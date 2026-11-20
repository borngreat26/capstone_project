"""Resumable, training-only P+T experiment; default: a stratified 150-gene pilot.

Use --stage fit to freeze cohort/models before --stage evaluate accesses test expression.
--all-genes is opt-in and is never launched by the pilot. Existing discovery is read only.
"""
import argparse
import json
import time
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.prediction import (THRESHOLDS, WORK, RESULTS, PRIMARY_RUN, digest, split_ids, load_expression,
                            chromosome_data, gene_data, training_associations, fit_gene,
                            scores_from_frame, metrics)
from src.utils import REPO_ROOT


def balanced_pick(frame, n, seed, sampler='chromosome-interleaved'):
    rng = np.random.default_rng(seed)
    groups = list(frame.groupby(['Gene_Chr', 'Length_Bin', 'Density_Bin', 'Variability_Bin'], sort=True))
    buckets = [rng.permutation(g.index).tolist() for _, g in groups]
    if sampler == 'chromosome-interleaved':
        # 'sorted' (pilot150) visited strata in string-sorted chromosome order, so a partial
        # pass drew only from chromosomes 1, 10, 11, ... Instead, shuffle strata within each
        # chromosome and take one stratum per chromosome in turn.
        chromosomes = [key[0] for key, _ in groups]
        order = rng.permutation(len(buckets))
        rank, seen = {}, {}
        for i in order:
            rank[i] = seen.get(chromosomes[i], 0); seen[chromosomes[i]] = rank[i] + 1
        buckets = [buckets[i] for i in sorted(order, key=lambda i: (rank[i], int(chromosomes[i])))]
    chosen = []
    while len(chosen) < min(n, len(frame)):
        for bucket in buckets:
            if bucket and len(chosen) < n:
                chosen.append(bucket.pop())
    return frame.loc[chosen].copy()


def prepare_genes(expression, groups):
    genes = pd.read_csv(REPO_ROOT/'data/quarter1_results/part1/protein_coding_genes.tsv',
                        sep='\t', dtype={'Gene_Chr': str})
    analyzed = pd.read_csv(REPO_ROOT/'results/part1/table_B_gene_summary.tsv',
                           sep='\t', usecols=['Gene_ID']).Gene_ID
    genes = genes.loc[genes.Gene_ID.isin(analyzed)].copy()
    assert len(genes) == 14661
    genes['Gene_Length'] = genes.Gene_End - genes.Gene_Start + 1
    # Supervised discovery p-values and full-sample expression variances are never read here.
    training = expression.loc[genes.Gene_ID, groups['train']].to_numpy(float)
    genes['Training_Expression_SD'] = training.std(axis=1, ddof=1)
    genes['Training_Expression_Mean'] = training.mean(axis=1)
    for chrom, group in genes.groupby('Gene_Chr'):
        positions = np.sort(pd.read_csv(REPO_ROOT/f'data/LDREF/1000G.EUR.{chrom}.bim',
                                       sep=r'\s+', header=None, usecols=[3])[3])
        genes.loc[group.index, 'Panel_Cis_SNPs'] = (
            np.searchsorted(positions, group.Gene_End+500000, side='right')
            - np.searchsorted(positions, np.maximum(1, group.Gene_Start-500000), side='left'))
    genes['SNPs_Per_Mb'] = genes.Panel_Cis_SNPs / ((genes.Gene_Length+1000000)/1e6)
    for source, destination in [('Gene_Length', 'Length_Bin'), ('SNPs_Per_Mb', 'Density_Bin'),
                                ('Training_Expression_SD', 'Variability_Bin')]:
        genes[destination] = pd.qcut(genes[source].rank(method='first'), 3, labels=False)
    genes = genes.sort_values('Gene_ID').reset_index(drop=True)
    return genes


def fit(args):
    out, work = RESULTS/args.run, WORK/args.run
    out.mkdir(parents=True, exist_ok=True); work.mkdir(parents=True, exist_ok=True)
    groups = split_ids()
    config = {'Schema': 1, 'Genes': args.genes, 'Candidates': args.candidates, 'Seed': args.seed,
              'All_Genes': args.all_genes, 'Thresholds': list(THRESHOLDS), 'LD_R2': .1, 'LD_KB': 250,
              'Evidence_Strata': {'Strong': 'training min p < 1e-8',
                                  'Moderate': '1e-8 <= training min p < 1e-3',
                                  'Weak': 'training min p >= 1e-3'},
              'Selection': 'Equal training-evidence quotas; balanced chromosome/length/density/training-variability bins; strongest training gene retained as case study',
              'Validation_Tie_Rule': 'Largest p-value cutoff on a predictive R2 tie',
              'Input_Expression_SHA256': digest(REPO_ROOT/'data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt'),
              'Split_SHA256': {k: digest(REPO_ROOT/f'data/week4_results/{k}.keep') for k in groups},
              'Discovery_Table_A_SHA256': digest(REPO_ROOT/'data/quarter1_results/part1/table_A_genome_wide_cis_eqtl.tsv.gz')}
    if args.sampler != 'sorted':  # pilot150 predates this key; omitting it keeps its frozen plan valid
        config['Sampler'] = args.sampler
    plan = out/'plan.json'
    if plan.exists() and json.loads(plan.read_text()) != config:
        raise ValueError('Run configuration/input changed. Use a new --run name; frozen runs are immutable.')
    if not plan.exists():
        plan.write_text(json.dumps(config, indent=2)+'\n')
    start = time.perf_counter()
    expression = load_expression()
    if not (out/'candidate_genes.tsv').exists():
        genes = prepare_genes(expression, groups)
        pool = genes if args.all_genes else balanced_pick(genes, args.candidates, args.seed, args.sampler)
        # User-proposed case-study candidate, declared before any validation or test outcomes.
        if not args.all_genes and 'RPS26' not in set(pool.Gene_Symbol):
            pool = pd.concat([pool.iloc[:-1], genes.loc[genes.Gene_Symbol.eq('RPS26')]])
        pool.sort_values(['Gene_Chr', 'Gene_Start']).to_csv(out/'candidate_genes.tsv', sep='\t', index=False)
    pool = pd.read_csv(out/'candidate_genes.tsv', sep='\t', dtype={'Gene_Chr': str})
    feature_rows = []
    for chrom, subset in pool.groupby('Gene_Chr', sort=True):
        variants, genotypes, columns = chromosome_data(chrom)
        for gene in subset.itertuples(index=False):
            association_path = work/'training_associations'/f'{gene.Gene_ID}.tsv.gz'
            association_path.parent.mkdir(exist_ok=True)
            if association_path.exists():
                associations = pd.read_csv(association_path, sep='\t')
            else:
                frame = gene_data(gene, variants, genotypes, columns)
                y = expression.loc[gene.Gene_ID, groups['train']].astype(float)
                associations = training_associations(frame, groups['train'], y)
                partial = association_path.with_name(association_path.name+'.partial')
                associations.to_csv(partial, sep='\t', index=False, compression='gzip')
                partial.replace(association_path)
            pmin = float(associations.P.min()) if len(associations) else np.nan
            evidence = 'Strong' if pmin < 1e-8 else 'Moderate' if pmin < 1e-3 else 'Weak'
            feature_rows.append({**gene._asdict(), 'Training_Min_P': pmin,
                                 'Training_Eligible_SNPs': len(associations),
                                 'Training_Nominal_Hits': int(associations.P.lt(.05).sum()),
                                 'Training_Gene_Bonferroni_Hits': int(associations.P.lt(.05/max(len(associations),1)).sum()),
                                 'Evidence': evidence,
                                 'Exclusion': '' if len(associations) and gene.Training_Expression_SD > 0 else 'No training-variable SNPs or expression'})
        print('Training candidate regressions complete: chromosome', chrom, flush=True)
    features = pd.DataFrame(feature_rows)
    features.to_csv(out/'candidate_training_features.tsv', sep='\t', index=False)
    if not (out/'cohort.tsv').exists():
        valid = features.loc[features.Exclusion.eq('')].copy()
        if args.all_genes:
            cohort = valid
        else:
            pieces = [balanced_pick(valid.loc[valid.Evidence.eq(tier)], args.genes//3, args.seed+i, args.sampler)
                      for i,tier in enumerate(['Strong','Moderate','Weak'])]
            cohort = pd.concat(pieces)
            if len(cohort) < args.genes:
                remaining = valid.loc[~valid.Gene_ID.isin(cohort.Gene_ID)]
                cohort = pd.concat([cohort, balanced_pick(remaining, args.genes-len(cohort), args.seed+4, args.sampler)])
            focus = valid.sort_values(['Training_Min_P','Gene_ID']).iloc[0]
            if focus.Gene_ID not in set(cohort.Gene_ID):
                same = cohort.index[cohort.Evidence.eq(focus.Evidence)]
                cohort = cohort.drop(same[-1] if len(same) else cohort.index[-1])
                cohort = pd.concat([cohort, focus.to_frame().T])
        cohort = cohort.sort_values(['Gene_Chr','Gene_Start']).reset_index(drop=True)
        cohort['Case_Study'] = cohort.Gene_ID.eq(cohort.sort_values(['Training_Min_P','Gene_ID']).iloc[0].Gene_ID)
        cohort.to_csv(out/'cohort.tsv', sep='\t', index=False)
    cohort = pd.read_csv(out/'cohort.tsv', sep='\t', dtype={'Gene_Chr': str})
    records = []; newly_fitted = 0
    genotype_provenance = {}
    for chrom, subset in cohort.groupby('Gene_Chr', sort=True):
        variants, genotypes, columns = chromosome_data(chrom)
        genotype_provenance[str(chrom)] = json.loads((WORK/'genotypes'/f'chr{chrom}.json').read_text())
        for gene in subset.itertuples(index=False):
            destination = out/'models'/gene.Gene_ID
            if (destination/'frozen_model.json').exists():
                model = json.loads((destination/'frozen_model.json').read_text())
                if digest(destination/'weights.tsv') != model['Weights_SHA256']:
                    raise ValueError('Frozen model weights changed')
            else:
                frame = gene_data(gene, variants, genotypes, columns)
                associations = pd.read_csv(work/'training_associations'/f'{gene.Gene_ID}.tsv.gz', sep='\t')
                model = fit_gene(gene, frame, associations,
                                 expression.loc[gene.Gene_ID, groups['train']].to_numpy(float),
                                 expression.loc[gene.Gene_ID, groups['validation']].to_numpy(float),
                                 groups, destination)
                newly_fitted += 1
            records.append(model)
            if args.max_models and newly_fitted >= args.max_models and len(records) < len(cohort):
                (out/'batch_checkpoint.json').write_text(json.dumps({
                    'New_Models_This_Batch': newly_fitted,
                    'Completed_Model_Files': len(list((out/'models').glob('*/frozen_model.json'))),
                    'Cohort_Genes': len(cohort), 'Ready_For_Test_Evaluation': False}, indent=2)+'\n')
                print('Batch checkpoint saved; resume the same fit command to continue.', flush=True)
                return False
        print('Training/validation model fits complete: chromosome', chrom, flush=True)
    pd.DataFrame(records).to_csv(out/'validation_selected_models.tsv', sep='\t', index=False)
    provenance_path = out/'genotype_provenance.json'
    if provenance_path.exists() and json.loads(provenance_path.read_text()) != genotype_provenance:
        raise ValueError('Frozen genotype sources have changed; use a new experiment run')
    provenance_path.write_text(json.dumps(genotype_provenance, indent=2)+'\n')
    frozen = {'Plan_SHA256': digest(plan), 'Cohort_SHA256': digest(out/'cohort.tsv'),
              'Models': {r['Gene_ID']: digest(out/'models'/r['Gene_ID']/'frozen_model.json') for r in records}}
    freeze_path = out/'frozen_experiment.json'
    if freeze_path.exists() and json.loads(freeze_path.read_text()) != frozen:
        raise ValueError('Frozen experiment cannot be changed after selection')
    freeze_path.write_text(json.dumps(frozen, indent=2)+'\n')
    timing_path = out/'fit_runtime.json'
    if not timing_path.exists():
        timing_path.write_text(json.dumps({'Fit_Seconds': time.perf_counter()-start,
                                          'Candidate_Genes': len(pool), 'Fitted_Genes': len(cohort),
                                          'Genotype_Cache_Bytes': sum(p.stat().st_size for p in (WORK/'genotypes').glob('*'))}, indent=2)+'\n')
    print('Models frozen before test evaluation:', len(cohort), flush=True)
    return True


def evaluate(args):
    out, work = RESULTS/args.run, WORK/args.run
    freeze_path = out/'frozen_experiment.json'
    frozen = json.loads(freeze_path.read_text()); signature = digest(freeze_path)
    assert digest(out/'plan.json') == frozen['Plan_SHA256']
    assert digest(out/'cohort.tsv') == frozen['Cohort_SHA256']
    config = json.loads((out/'plan.json').read_text())
    provenance = json.loads((out/'genotype_provenance.json').read_text())
    assert digest(REPO_ROOT/'data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt') == config['Input_Expression_SHA256']
    if (out/'test_evaluation.json').exists():
        assert json.loads((out/'test_evaluation.json').read_text())['Frozen_Experiment_SHA256'] == signature
        for gene_id, expected_hash in frozen['Models'].items():
            folder=out/'models'/gene_id
            assert digest(folder/'frozen_model.json')==expected_hash
            assert digest(folder/'weights.tsv')==json.loads((folder/'frozen_model.json').read_text())['Weights_SHA256']
        print('Using verified saved test evaluation; no repeated test evaluation.', flush=True)
        return
    groups = split_ids(); expression = load_expression()
    cohort = pd.read_csv(out/'cohort.tsv', sep='\t', dtype={'Gene_Chr': str})
    results = []; start = time.perf_counter()
    for chrom, subset in cohort.groupby('Gene_Chr', sort=True):
        variants, genotypes, columns = chromosome_data(chrom)
        if json.loads((WORK/'genotypes'/f'chr{chrom}.json').read_text()) != provenance[str(chrom)]:
            raise ValueError('Evaluation genotype source differs from the fitted model source')
        for gene in subset.itertuples(index=False):
            destination = out/'models'/gene.Gene_ID
            assert digest(destination/'frozen_model.json') == frozen['Models'][gene.Gene_ID]
            model = json.loads((destination/'frozen_model.json').read_text())
            assert digest(destination/'weights.tsv') == model['Weights_SHA256']
            prediction_path = work/'test_predictions'/f'{gene.Gene_ID}.tsv'
            metric_path = work/'test_predictions'/f'{gene.Gene_ID}.json'
            prediction_path.parent.mkdir(exist_ok=True)
            if metric_path.exists():
                record = json.loads(metric_path.read_text())
                assert record['Frozen_Experiment_SHA256'] == signature
                results.append(record)
                continue
            observed = expression.loc[gene.Gene_ID, groups['test']].to_numpy(float)
            baseline = np.repeat(model['Training_Mean'], len(observed))
            baseline_r2 = metrics(observed, baseline)['Predictive_R2']
            weights = pd.read_csv(destination/'weights.tsv', sep='\t')
            if len(weights):
                frame = gene_data(gene, variants, genotypes, columns)
                score = scores_from_frame(frame, weights, groups['test'])
                predicted = model['Intercept'] + model['Slope']*score
                measured = metrics(observed, predicted)
            else:
                score = np.zeros(len(observed)); predicted = baseline
                measured = dict(Predictive_R2=np.nan, Pearson_r=np.nan, Pearson_r_squared=np.nan)
            pd.DataFrame({'IID': groups['test'], 'Observed': observed, 'Genetic_Score': score,
                          'Predicted': predicted, 'Training_Mean_Baseline': baseline}).to_csv(prediction_path, sep='\t', index=False)
            record = {**gene._asdict(), **model, 'Test_N': len(observed),
                      'Test_Predictive_R2': measured['Predictive_R2'],
                      'Test_Pearson_r': measured['Pearson_r'], 'Test_Pearson_r_squared': measured['Pearson_r_squared'],
                      'Baseline_Test_R2': baseline_r2,
                      'Test_Beats_Baseline': bool(measured['Predictive_R2'] > baseline_r2) if len(weights) else False,
                      'Frozen_Experiment_SHA256': signature}
            partial=metric_path.with_suffix('.json.partial')
            partial.write_text(json.dumps(record, indent=2)+'\n')
            partial.replace(metric_path)
            results.append(record)
    pd.DataFrame(results).to_csv(out/'gene_prediction_results.tsv', sep='\t', index=False)
    (out/'test_evaluation.json').write_text(json.dumps({'Frozen_Experiment_SHA256': signature,
        'Genes_Evaluated': len(results), 'Test_N': len(groups['test']),
        'Test_Evaluation_Seconds': time.perf_counter()-start,
        'Protocol': 'Frozen models evaluated once per gene; reruns reuse saved predictions/metrics. Test results are not used to change selection.'}, indent=2)+'\n')
    print('Held-out evaluation complete:', len(results), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', default=PRIMARY_RUN)
    ap.add_argument('--sampler', choices=['chromosome-interleaved', 'sorted'], default='chromosome-interleaved',
                    help="'sorted' reproduces the superseded pilot150 draw")
    ap.add_argument('--genes', type=int, default=150)
    ap.add_argument('--candidates', type=int, default=660)
    ap.add_argument('--seed', type=int, default=20261120)
    ap.add_argument('--stage', choices=['fit','evaluate','all'], default='fit')
    ap.add_argument('--all-genes', action='store_true', help='Explicit opt-in after reviewing pilot resource estimates')
    ap.add_argument('--max-models', type=int, default=0, help='Fit at most this many new models per invocation; 0 means no batch limit')
    args = ap.parse_args()
    if args.all_genes and args.run in ('pilot150', PRIMARY_RUN):
        ap.error('Use a separate --run name for a genome-wide experiment')
    if Path(args.run).name != args.run or args.run in ('.','..') or args.max_models < 0:
        ap.error('Use a simple run name and a nonnegative model batch limit')
    complete = fit(args) if args.stage in ('fit','all') else True
    if args.stage in ('evaluate','all') and complete: evaluate(args)

if __name__ == '__main__': main()
