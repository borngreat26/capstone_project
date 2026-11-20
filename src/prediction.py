"""Training-only P+T expression prediction, using the existing eQTL and PRS modules."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, linregress

from src.eqtl import fit_cis
from src.prs import clump
from src.utils import REPO_ROOT

THRESHOLDS = (1., .05, .01, .001, .00001)
WORK = REPO_ROOT / 'data/prediction_work'
RESULTS = REPO_ROOT / 'results/expression_prediction'
# pilot150 is kept frozen but superseded: its sampler drew moderate/weak genes from few chromosomes.
PRIMARY_RUN = 'pilot150_v2'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def metrics(observed, predicted):
    observed, predicted = np.asarray(observed, float), np.asarray(predicted, float)
    if observed.shape != predicted.shape or not (np.isfinite(observed).all() and np.isfinite(predicted).all()):
        raise ValueError('Prediction/observation arrays must be finite and aligned')
    sst = np.square(observed - observed.mean()).sum()
    r = float(pearsonr(observed, predicted).statistic) if np.ptp(observed) > 0 and np.ptp(predicted) > 0 else np.nan
    return {'Predictive_R2': 1 - np.square(observed-predicted).sum()/sst if sst > 0 else np.nan,
            'Pearson_r': r, 'Pearson_r_squared': r*r}


def best_threshold(validation):
    valid = validation.loc[validation.SNPs.gt(0) & validation.Predictive_R2.notna()
                           & validation.Status.eq('Genetic score')]
    if valid.empty:
        return None
    return valid.sort_values(['Predictive_R2', 'Threshold'], ascending=[False, False]).iloc[0]


def split_ids():
    groups = {}
    for name in ['train', 'validation', 'test']:
        groups[name] = pd.read_csv(REPO_ROOT/f'data/week4_results/{name}.keep', sep=r'\s+',
                                   header=None, dtype=str)[1].tolist()
    if [len(groups[k]) for k in groups] != [275, 34, 35]:
        raise ValueError('Unexpected original Week 4 split')
    if len(set(sum(groups.values(), []))) != 344:
        raise ValueError('Participant groups overlap')
    return groups


def load_expression():
    return pd.read_csv(REPO_ROOT/'data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt',
                       sep='\t').set_index('TargetID')


def chromosome_data(chrom):
    """Reuse the discovery genotype export; cache dosages, never its full-sample betas."""
    source = REPO_ROOT/f'data/quarter1_results/part1/chromosomes/chr{chrom}.traw'
    folder = WORK/'genotypes'; folder.mkdir(parents=True, exist_ok=True)
    prefix = folder/f'chr{chrom}'
    stamp = {'Size': source.stat().st_size, 'Mtime_ns': source.stat().st_mtime_ns}
    marker = prefix.with_suffix('.json')
    if not marker.exists() or json.loads(marker.read_text())['Source'] != stamp:
        frame = pd.read_csv(source, sep=r'\s+')
        columns = frame.columns[6:].tolist()
        g = frame[columns].to_numpy(float)
        if not np.isin(g[np.isfinite(g)], [0, 1, 2]).all():
            raise ValueError('Unexpected hard-call dosage')
        g = np.where(np.isfinite(g), g, -1).astype(np.int8)
        np.save(prefix.with_suffix('.npy'), g)
        frame.iloc[:, :6].to_csv(prefix.with_suffix('.tsv'), sep='\t', index=False)
        marker.write_text(json.dumps({'Source': stamp, 'Source_SHA256': digest(source),
                                      'Columns': columns}, indent=2)+'\n')
    meta = json.loads(marker.read_text())
    variants = pd.read_csv(prefix.with_suffix('.tsv'), sep='\t')
    g = np.load(prefix.with_suffix('.npy'), mmap_mode='r')
    return variants, g, meta['Columns']


def gene_data(gene, variants, genotypes, columns):
    region = variants.POS.between(max(1, int(gene.Gene_Start)-500000), int(gene.Gene_End)+500000)
    frame = variants.loc[region].reset_index(drop=True).copy()
    g = np.asarray(genotypes[region], dtype=float)
    g[g < 0] = np.nan
    return pd.concat([frame, pd.DataFrame(g, columns=columns)], axis=1)


def training_associations(frame, ids, expression):
    # Explicitly resolve each participant's exported genotype column.
    mapping = {iid: iid+'_'+iid for iid in ids}
    result = fit_cis(frame, ids, expression, mapping).rename(columns={'COUNTED': 'A1'})
    # Week 4 uses complete calls. Incomplete training SNPs are not model candidates.
    return result.loc[result.N.eq(len(ids))].reset_index(drop=True)


def scores_from_frame(frame, weights, ids):
    if weights.empty:
        return np.zeros(len(ids))
    ordered = frame.set_index('SNP').loc[weights.SNP]
    if ordered.index.duplicated().any() or not np.array_equal(ordered.COUNTED, weights.A1):
        raise ValueError('Expression-score counted alleles do not match genotype dosages')
    g = ordered[[iid+'_'+iid for iid in ids]].to_numpy(float)
    if not np.isfinite(g).all():
        raise ValueError('Missing calls in scoring SNPs; no held-out imputation is performed')
    return weights.BETA.to_numpy() @ g


def fit_gene(gene, frame, associations, train_y, validation_y, groups, destination,
             bfile=None, keep=None):
    destination.mkdir(parents=True, exist_ok=True)
    bfile = bfile or REPO_ROOT/f'data/LDREF/1000G.EUR.{int(gene.Gene_Chr)}'
    keep = keep or REPO_ROOT/'data/week4_results/train.keep'
    indexes = clump(associations, bfile, destination/'training_clump', keep=keep)
    independent = associations.loc[associations.SNP.isin(indexes)].copy()
    rows = []
    baseline = metrics(validation_y, np.repeat(train_y.mean(), len(validation_y)))['Predictive_R2']
    for cutoff in THRESHOLDS:
        weights = independent.loc[independent.P.le(cutoff)]
        train_score = scores_from_frame(frame, weights, groups['train'])
        if len(weights) and np.std(train_score) > 0:
            calibration = linregress(train_score, train_y)
            intercept, slope, status = calibration.intercept, calibration.slope, 'Genetic score'
        else:
            intercept, slope = train_y.mean(), 0.
            status = 'Training-mean baseline (0 SNPs)' if weights.empty else 'Constant training genetic score'
        prediction = intercept + slope*scores_from_frame(frame, weights, groups['validation'])
        rows.append({'Threshold': cutoff, 'SNPs': len(weights), **metrics(validation_y, prediction),
                     'Intercept': intercept, 'Slope': slope, 'Status': status})
    validation = pd.DataFrame(rows)
    validation.to_csv(destination/'validation.tsv', sep='\t', index=False)
    selected = best_threshold(validation)
    model = {'Gene_ID': gene.Gene_ID, 'Gene': gene.Gene_Symbol,
             'Training_N': len(train_y), 'Validation_N': len(validation_y),
             'Eligible_SNPs': len(associations), 'Clumped_SNPs': len(independent),
             'Training_Mean': float(train_y.mean()), 'Baseline_Validation_R2': baseline,
             'Status': 'No valid nonempty genetic model' if selected is None else 'Genetic model',
             'SNPs': 0 if selected is None else int(selected.SNPs)}
    weights = independent.iloc[:0]
    if selected is not None:
        weights = independent.loc[independent.P.le(selected.Threshold)]
        model.update({'Threshold': float(selected.Threshold), 'Intercept': float(selected.Intercept),
                      'Slope': float(selected.Slope), 'Validation_Predictive_R2': float(selected.Predictive_R2),
                      'Validation_Beats_Baseline': bool(selected.Predictive_R2 > baseline)})
    weights.to_csv(destination/'weights.tsv', sep='\t', index=False)
    model['Weights_SHA256'] = digest(destination/'weights.tsv')
    # This is written before any test expression is passed to evaluation.
    marker = destination/'frozen_model.json.partial'
    marker.write_text(json.dumps(model, indent=2)+'\n')
    marker.replace(destination/'frozen_model.json')
    return model
