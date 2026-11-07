"""Render the LocusZoom.js page to a PNG and record exactly what was drawn.

Usage: python scripts/render_locuszoom.py page.html output.png drawn_points.json
Needs Playwright with Chromium: pip install playwright && playwright install chromium-headless-shell
"""
import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

html_path, png_path, points_path = sys.argv[1:4]

READ_DRAWN_POINTS = """() => {
    const panel = lzPlot.panels.association;
    const points = [];
    d3.selectAll('#lz path.lz-data_layer-scatter').each(function (d) {
        const [x, y] = this.getAttribute('transform').match(/-?[\\d.]+/g).map(Number);
        points.push({SNP: d['assoc:rsid'], Position: d['assoc:position'], Minus_Log10_P: d['assoc:log_pvalue'],
                     LD_R2: d['ld:correlation'], Lead_SNP: d.lz_is_ld_refvar === true, Fill: this.style.fill || this.getAttribute('fill'),
                     Shape_path: this.getAttribute('d'), Pixel_position: panel.x_scale.invert(x),
                     Pixel_minus_log10_p: panel.y1_scale.invert(y)});
    });
    return {points, state: lzPlot.state, x_extent: panel.x_extent, y1_extent: panel.y1_extent,
            genes: lzPlot.panels.genes.data_layers.genes.data.map((g) => g.gene_name)};
}"""

problems, outside_requests = [], []
with sync_playwright() as playwright:
    browser = playwright.chromium.launch()
    page = browser.new_page(viewport={'width': 1240, 'height': 1000}, device_scale_factor=2)
    page.on('pageerror', lambda error: problems.append(str(error)))
    page.on('console', lambda m: problems.append(m.text) if m.type == 'error' else None)
    page.on('request', lambda r: outside_requests.append(r.url) if r.url.startswith('http') else None)
    page.goto(Path(html_path).resolve().as_uri())
    page.wait_for_selector('#lz path.lz-data_layer-scatter', state='attached')
    page.wait_for_selector('#lz rect.lz-data_layer-genes', state='attached')
    page.wait_for_timeout(1500)
    drawn = page.evaluate(READ_DRAWN_POINTS)
    page.add_style_tag(content='[class*="lz-toolbar-button"] { visibility: hidden; }')  # hide buttons in the static copy only
    page.locator('#figure').screenshot(path=png_path)
    browser.close()

drawn.update(problems=problems, outside_requests=outside_requests)
Path(points_path).write_text(json.dumps(drawn))
print(f'Rendered {len(drawn["points"])} points; browser errors: {len(problems)}; outside requests: {len(outside_requests)}')
