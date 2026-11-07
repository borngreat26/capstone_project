"""Build a self-contained LocusZoom.js page for one cis-eQTL locus.

The page holds its data and libraries, so it opens offline. scripts/render_locuszoom.py turns it into a PNG.
"""
import json
import urllib.request
from pathlib import Path

import numpy as np

VENDOR_FILES = {'d3.min.js': 'https://cdn.jsdelivr.net/npm/d3@5.16.0/dist/d3.min.js',
                'locuszoom.app.min.js': 'https://cdn.jsdelivr.net/npm/locuszoom@0.14.0/dist/locuszoom.app.min.js',
                'locuszoom.css': 'https://cdn.jsdelivr.net/npm/locuszoom@0.14.0/dist/locuszoom.css'}

PAGE_TEMPLATE = """<!doctype html>
<html><head><meta charset="utf-8"><title>__SYMBOL__ LocusZoom</title>
<style>__CSS__
body { margin: 0; background: white; font-family: Helvetica, Arial, sans-serif; color: #222; }
#figure { width: 1200px; padding: 16px 20px; box-sizing: border-box; background: white; }
h1 { font-size: 22px; margin: 0 0 4px; } .sub { font-size: 14px; margin: 0 0 10px; } .cap { font-size: 12px; color: #555; margin: 8px 0 0; }
</style>
<script>__D3__</script>
<script>__LZ__</script>
</head><body><div id="figure">
<h1>__SYMBOL__ cis-eQTL regional association</h1><p class="sub">__HEADER__</p>
<div id="lz"></div><p class="cap">__CAPTION__</p></div>
<script>
const D = __DATA__;
const sources = new LocusZoom.DataSources()
    .add('assoc', ['StaticSource', {data: D.assoc, prefix_namespace: false}])
    .add('ld', ['StaticSource', {data: D.ld}])
    .add('gene', ['StaticSource', {data: D.genes, prefix_namespace: false}]);
const layout = LocusZoom.Layouts.get('plot', 'standard_association', {
    width: 1160, responsive_resize: false, max_region_scale: D.end - D.start + 1,
    state: {chr: D.chr, start: D.start, end: D.end, genome_build: 'GRCh37'},
});
const panel = layout.panels.find((p) => p.id === 'association');
layout.toolbar.widgets = layout.toolbar.widgets.filter((w) => w.button_html !== 'LD Population: ');
panel.height = 400;
panel.margin.right = 30;
panel.axes.x.label = 'Chromosome __CHROM__ position, GRCh37 (Mb)';
panel.axes.y1.label = '-log10(cis-eQTL p-value)';
delete panel.axes.y2;
panel.legend.hidden = false;
panel.legend.origin = {x: 1000, y: 40};
panel.data_layers = panel.data_layers.filter((l) => l.id !== 'recombrate');
Object.assign(panel.data_layers.find((l) => l.id === 'significance'), {offset: D.cutoff_log10,
    style: {stroke: '#444', 'stroke-width': '1.5px', 'stroke-dasharray': '6px'}});
const points = panel.data_layers.find((l) => l.id === 'associationpvalues');
points.coalesce.active = false;
points.legend.push({shape: 'diamond', color: '#9632b8', size: 80, label: 'Lead SNP', label_size: 12});
points.label = {text: '{{assoc:rsid}}', spacing: 8, lines: {style: {stroke: '#333', 'stroke-width': '1px'}},
    filters: [{field: 'lz_is_ld_refvar', operator: '=', value: true}], style: {'font-size': '13px', 'font-weight': 'bold'}};
points.tooltip.html = '<strong>{{assoc:rsid|htmlescape}}</strong><br>Position: {{assoc:position}}<br>' +
    'P value: <strong>{{assoc:log_pvalue|logtoscinotation|htmlescape}}</strong><br>' +
    'r² with lead SNP: <strong>{{ld:correlation|scinotation|htmlescape}}</strong><br>' +
    '{{#if lz_is_ld_refvar}}<strong>Lead SNP (LD reference)</strong>{{/if}}';
const genePanel = layout.panels.find((p) => p.id === 'genes');
genePanel.margin.right = 30;
genePanel.height = 190;
const genes = genePanel.data_layers.find((l) => l.id === 'genes');
genes.namespace = {gene: 'gene'};
genes.data_operations = [{type: 'fetch', from: ['gene']}];
genes.tooltip.html = '<strong>{{gene_name|htmlescape}}</strong><br>{{gene_id|htmlescape}}<br>{{gene_type|htmlescape}}';
const lzPlot = LocusZoom.populate('#lz', sources, layout);
</script></body></html>"""


def download_vendor_files(vendor_dir):
    """Fetch the pinned LocusZoom.js and D3 files once."""
    vendor_dir = Path(vendor_dir)
    vendor_dir.mkdir(exist_ok=True)
    for name, url in VENDOR_FILES.items():
        if not (vendor_dir / name).exists():
            urllib.request.urlretrieve(url, vendor_dir / name)


def merge_exons(rows):
    """Collapse all transcripts' exons, because the gene track draws one model per gene."""
    merged = []
    for start, end in sorted(zip(rows[3], rows[4])):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def gene_track(gtf, annotation, chrom, start, end):
    """Gene track records for every annotated gene overlapping the window. `gtf` is a headerless GENCODE table."""
    exon_rows = gtf.loc[gtf[2].eq('exon')].assign(Gene_ID=lambda t: t[8].str.extract(r'gene_id "([^"]+)"')[0])
    track = []
    for gene, info in annotation.loc[(annotation['Start'] <= end) & (annotation['End'] >= start)].iterrows():
        exons = [{'exon_id': f'{gene}_{k}', 'start': int(s), 'end': int(e), 'strand': info['Strand']}
                 for k, (s, e) in enumerate(merge_exons(exon_rows.loc[exon_rows['Gene_ID'].eq(gene)]))]
        track.append({'gene_id': gene, 'gene_name': info['Symbol'], 'gene_type': info['Gene_type'], 'chrom': chrom,
                      'start': int(info['Start']), 'end': int(info['End']), 'strand': info['Strand'],
                      'transcripts': [{'transcript_id': f'{gene}_merged', 'exons': exons}]})
    return track


def write_locuszoom_page(path, vendor_dir, ld_table, track, chrom, start, end, symbol, gene_id, lead_snp, lead_position,
                         lead_p, cutoff, people):
    """Write the page. `ld_table` needs SNP, Position, Minus_Log10_P, Lead_SNP, and LD_R2."""
    vendor_dir = Path(vendor_dir)
    download_vendor_files(vendor_dir)
    data = {
        'chr': chrom, 'start': start, 'end': end, 'symbol': symbol, 'lead': lead_snp,
        'cutoff_log10': float(-np.log10(cutoff)),
        'assoc': [{'assoc:variant': r.SNP, 'assoc:rsid': r.SNP, 'assoc:position': int(r.Position), 'assoc:log_pvalue': r.Minus_Log10_P,
                   'lz_is_ld_refvar': bool(r.Lead_SNP)} for r in ld_table.itertuples()],
        'ld': [{'variant2': r.SNP, 'position2': int(r.Position), 'correlation': r.LD_R2} for r in ld_table.itertuples()],
        'genes': track,
    }
    header = (f'{symbol} ({gene_id}) · chr{chrom}:{start:,}–{end:,} (GRCh37; gene body ±500 kb) · '
              f'lead SNP {lead_snp}, chr{chrom}:{lead_position:,}, p = {lead_p:.2e}')
    caption = (f'{len(ld_table)} tested SNPs from additive cis-eQTL regressions in {people} people. '
               f'Color: unphased genotype r² with the lead SNP from PLINK2 in the same people. Dashed line: chr{chrom} Bonferroni p = {cutoff:.2e}. '
               'Genes: GENCODE v12, exons merged across transcripts. No recombination track.')
    page = PAGE_TEMPLATE.replace('__CHROM__', chrom)
    for key, value in {'__CSS__': (vendor_dir / 'locuszoom.css').read_text(), '__D3__': (vendor_dir / 'd3.min.js').read_text(),
                       '__LZ__': (vendor_dir / 'locuszoom.app.min.js').read_text(), '__SYMBOL__': symbol, '__HEADER__': header,
                       '__CAPTION__': caption, '__DATA__': json.dumps(data, allow_nan=False).replace('</', '<\\/')}.items():
        page = page.replace(key, value)
    Path(path).write_text(page)
