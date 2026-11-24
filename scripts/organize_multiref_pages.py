"""Assemble the five English report pages from the existing results and plot assets."""
import argparse
from pathlib import Path
import re
import shutil

TEMPLATES = Path(__file__).parent/'multiref_web'
PAGES = {'index':'Overview', 'results':'Results', 'selection':'Step 1 & coincidences',
         'aggregation':'REF samples & offsets', 'diagnostics':'Diagnostics'}


def write_pages(source, output):
    output.mkdir(parents=True, exist_ok=True)
    # The full generator writes the single-page results first. The standalone
    # command also accepts a previously organized report as input.
    original = (source/'index.html').read_text()
    if 'id="comparison"' not in original:
        original = (source/'results.html').read_text()
    sections = re.findall(r'<section\b[^>]*>.*?</section>', original, re.S)
    def section(heading):
        return next(s for s in sections if f'<h2>{heading}</h2>' in s)
    results = ''.join(section(h) for h in (
        'Explore the four methods', 'Centered union compared with REF7',
        'Results at T1 = 10', 'Baseline and interpretation', 'Download the evidence'))
    results += '<p><a href="commands.html">Commands and files</a> · <a href="methodology.html">Estimator equations</a></p>'
    titles = {'index':'Timing resolution and REF selection', 'results':'Timing closure results',
              'selection':'Step 1 and the coincidence selection', 'aggregation':'Combining the REF samples',
              'diagnostics':'Timing diagnostics'}
    for page, title in titles.items():
        nav = ''.join(f'<a href="{slug}.html"'+(' aria-current="page"' if slug==page else '')+
                      f'>{label.replace("&","&amp;")}</a>' for slug,label in PAGES.items())
        scripts = ''
        if page == 'results':
            body = results
            scripts = '<script defer src="data.js"></script><script defer src="app.js?v=centered"></script>'
        else:
            body = (TEMPLATES/(('overview' if page=='index' else page)+'.html')).read_text()
            if page == 'diagnostics':
                scripts = '<script defer src="diagnostics.js"></script>'
            if page == 'index':
                scripts = """<script>
const anchor=location.hash.split('?')[0];
if (['#comparison','#comparisons','#results'].includes(anchor))
  location.replace('results.html'+(anchor==='#comparisons'?'#comparison':anchor));
if (anchor==='#diagnostics') location.replace('diagnostics.html'+location.hash);
</script>"""
        (output/(page+'.html')).write_text(f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>{title} · CMS BTL</title>
<link rel="stylesheet" href="style.css">{scripts}</head><body><header><div class="shell mast"><div class="identity">
<strong>CMS BTL</strong><span>September 2025 · REF selection study</span></div></div></header><main class="shell">
<div class="intro"><h1>{title}</h1></div><nav aria-label="Analysis pages">{nav}</nav>{body}</main>
<footer class="shell"><a href="../index.html">Historical single-REF report</a></footer></body></html>''')
    for page in ('methodology','commands'):
        content=(source/(page+'.html')).read_text().replace('href="index.html">Results','href="results.html">Results')
        (output/(page+'.html')).write_text(content)
    shutil.copyfile(TEMPLATES/'diagnostics.js',output/'diagnostics.js')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    write_pages(args.source,args.output)
