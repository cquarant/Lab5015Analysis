#!/usr/bin/env python3
"""Build an English, static HTML report and scientific plots from completed multi-REF stages."""
import argparse
import csv
import html
import json
import math
from pathlib import Path
import shutil
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
REPO=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(REPO/'main'),str(REPO/'macros')]
from moduleCharacterization_timingClosure import dump_json
from multiref_statistics import METHODS,REF_BARS,covariance_diagnostics
from compare_timingClosure import plot_methods,root_or_none
from organize_multiref_pages import write_pages
from package_multiref_diagnostics import package_all

RUNS=(4194,4257,4281,4291)
MODULES={4194:'DM9003 · SM0800',4257:'DM9001 · SM0691',4281:'FE4587 · SM4436',4291:'DM9004 · SM1365'}

def write_csv(path,rows,empty_fields=()):
    fields=list(dict.fromkeys(k for r in rows for k in r)) or list(empty_fields)
    if not fields: return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

def gap(metric):
    if not metric or min(metric['avg_mixed'],metric['direct_halfdiff'])<0: return None
    return math.sqrt(metric['avg_mixed'])-math.sqrt(metric['direct_halfdiff'])

def plot_ref_comparison(output, run, rows):
    samples = {
        'centered': ('Union, REF offsets removed', '#006b65'),
        'ref07': ('REF7, same selection', '#2876ad'),
        'published': ('Published REF7', '#969696'),
    }
    methods = {
        'direct_halfdiff': 'Direct',
        'avg_mixed': 'Mixed',
        'matched_halfdiff': 'Corrected difference',
        'catr': 'CatR',
    }
    bars = np.arange(16)
    for threshold in (3,5,7,10,15,20):
        for estimator in ('moments','gaussian'):
            selected = {
                (r['sample'], r['bar']): r for r in rows
                if r['run'] == run and r['threshold'] == threshold
                and r['stage'] == 'full' and r['estimator'] == estimator
                and r['sample'] in samples
            }

            def series(sample, field):
                return np.array([selected.get((sample, bar), {}).get(field) for bar in bars], dtype=float)

            title = f'Run {run} · T1 {threshold} · {estimator}'
            fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, constrained_layout=True)
            fig.suptitle(f'{title} · centered union vs REF7')
            for ax, (method, label) in zip(axes.flat, methods.items()):
                for sample in ('centered', 'ref07'):
                    name, color = samples[sample]
                    ax.plot(bars, series(sample, method+'_ps'), 'o-', color=color, label=name)
                ax.set(title=label, ylabel='Resolution [ps]', xticks=bars)
                ax.grid(alpha=.25)
            for ax in axes[-1]:
                ax.set_xlabel('DUT bar')
            axes[0,0].legend(fontsize=8)
            for ext in ('png','pdf'):
                fig.savefig(output/f'methods_centered_vs_ref7_{run}_{threshold}_{estimator}.{ext}')
            plt.close(fig)

            fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
            gaps = {}
            for sample, (name, color) in samples.items():
                gaps[sample] = series(sample, 'mixed_minus_direct_ps')
                axes[0].plot(bars, gaps[sample], 'o-', color=color, label=name)
            axes[1].plot(bars, gaps['centered']-gaps['ref07'], 'o-', color=samples['centered'][1])
            axes[0].set(ylabel='Mixed − direct resolution [ps]', title=title)
            axes[0].legend(fontsize=9)
            axes[1].set(xlabel='DUT bar', ylabel='Centered − REF7 gap [ps]', xticks=bars)
            for ax in axes:
                ax.axhline(0, color='black', lw=.7)
                ax.grid(alpha=.25)
            for ext in ('png','pdf'):
                fig.savefig(output/f'delta_{run}_{threshold}_{estimator}.{ext}')
            plt.close(fig)


def findings_html(summaries, statistics, offsets, historical, common_events, rows, covariances):
    central = [r for r in summaries if r['selection']=='central_DUT_4_11']
    sharing = [r for r in statistics if r['ref']=='union']
    pairs = sum(r['accepted'] for r in sharing)
    unique = sum(r['unique_events'] for r in sharing)
    shared = sum(r['shared_events'] for r in sharing)
    text = (f'<p>The union contains {pairs:,} coincidences representing {unique:,} distinct '
            f'DUT-bar observations, of which {shared:,} enter more than one REF sample.</p>')
    selected_offsets=[r['average_ref_mean_span_ps'] for r in offsets if r['stage']=='full']
    if selected_offsets:
        text+=(f'<p>The median span of the REF-dependent means of the average residual is '
               f'{np.median(selected_offsets):.0f} ps across DUT bars and thresholds. '
               'The uncentered union preserves these offsets. A Gaussian fitted to one peak of a '
               'mixture can respond differently for each residual, so a pooled Gaussian gap '
               'can be large or yield a negative extracted variance. Compare with the centered '
               'results before giving that gap a timing-resolution interpretation.</p>')
    high_thresholds=[r for r in central if r['threshold']>=5 and r['estimator']=='gaussian']
    centered_gaps=[abs(r['centered_median_gap_ps']) for r in high_thresholds
                   if r['centered_median_gap_ps'] is not None]
    if centered_gaps:
        union_valid=sum(r['union_valid_groups'] for r in high_thresholds)
        centered_valid=sum(r['centered_valid_groups'] for r in high_thresholds)
        total=8*len(high_thresholds)
        text+=(f'<p>Across the central DUT bars at T1 ≥ 5, the Gaussian mixed/direct comparison '
               f'has {union_valid}/{total} valid gaps in the uncentered union and '
               f'{centered_valid}/{total} after centering. The median absolute centered '
               f'run-by-threshold gap is {np.median(centered_gaps):.2f} ps. '
               'See the T1 = 3 rows separately for the low-threshold behavior.</p>')
    cells=[]
    for run in RUNS:
        for estimator in ('gaussian','moments'):
            gaps={(r['sample'],r['bar']):r['mixed_minus_direct_ps'] for r in rows
                  if r['run']==run and r['threshold']==10 and r['stage']=='full'
                  and r['estimator']==estimator and r['sample'] in ('centered','ref07')
                  and r['mixed_minus_direct_ps'] is not None}
            selected=[gaps['centered',bar]-gaps['ref07',bar] for bar in range(4,12)
                      if ('centered',bar) in gaps and ('ref07',bar) in gaps]
            r=next(r for r in central if r['run']==run and r['threshold']==10 and r['estimator']==estimator)
            values=[r['ref07_median_gap_ps'],r['centered_median_gap_ps'],float(np.median(selected)) if selected else None]
            formatted=''.join(f'<td>{v:+.2f}</td>' if v is not None else '<td>—</td>' for v in values)
            cells.append(f'<tr><td>{run}</td><td>{estimator}</td>{formatted}<td>{len(selected)}</td></tr>')
    text+=('''<p>At T1 = 10 the table shows the median mixed-minus-direct gap over DUT bars 4–11.
        The change column compares the centered union with REF7 within each DUT bar and then takes the median.</p>
        <div class="scroll"><table><thead><tr><th>Run</th><th>Estimator</th><th>REF7 gap [ps]</th>
        <th>Centered union gap [ps]</th><th>Change vs REF7 [ps]</th><th>DUT bars</th></tr></thead><tbody>'''
        + ''.join(cells) + '</tbody></table></div>')
    html_rows=[]
    for run in RUNS:
        selected=[r for r in common_events if r['run']==run and r['stage']=='full'
                  and r['threshold']==10 and 4<=r['bar']<=11 and r['mixed_direct_gap_change_ps'] is not None]
        if selected:
            median=np.median([r['mixed_direct_gap_change_ps'] for r in selected])
            html_rows.append(f'<tr><td>{run}</td><td>{len(selected)}</td><td>{median:+.2f}</td></tr>')
    if html_rows:
        text+=('''<h3>Comparison on the same events</h3><p>For each comparison below, the selected DUT
              observations are present in both REF samples. The change uses exact moments with
              independently trained corrections. This comparison omits REF pairs with fewer than
              100 shared observations.</p><div class="scroll"><table><thead><tr><th>Run</th>
              <th>REF–DUT groups</th><th>Median gap change vs REF7 [ps]</th></tr></thead><tbody>'''
              +''.join(html_rows)+'</tbody></table></div>')
    low_threshold=[]
    for run in RUNS:
        values=[]
        for estimator in ('gaussian','moments'):
            for stage in ('tw','full'):
                gaps=[r['mixed_minus_direct_ps'] for r in rows if r['run']==run and r['sample']=='centered'
                      and r['threshold']==3 and r['stage']==stage and r['estimator']==estimator
                      and 4<=r['bar']<=11 and r['mixed_minus_direct_ps'] is not None]
                values.append(f'<td>{np.median(gaps):+.2f} ({len(gaps)})</td>' if gaps else '<td>—</td>')
        low_threshold.append(f'<tr><td>{run}</td>'+''.join(values)+'</tr>')
    text+=('''<h3>Lowest threshold and correction stage</h3><p>At T1 = 3 the centered-union median
        mixed-minus-direct gap is shown before and after phase correction, over DUT bars 4–11.
        Compare the change in the Gaussian core with the moment gap, which already includes
        the tails before phase correction.</p><div class="scroll"><table><thead>
        <tr><th>Run</th><th>Gaussian, time walk [ps]</th><th>Gaussian, +phase [ps]</th>
        <th>Moments, time walk [ps]</th><th>Moments, +phase [ps]</th></tr></thead><tbody>'''
        +''.join(low_threshold)+'</tbody></table></div>')
    covariance_rows=[]
    for run in RUNS:
        selected=[r for r in covariances if r['run']==run and r['sample']=='centered'
                  and r['stage']=='full' and r['threshold']==10 and 4<=r['bar']<=11]
        fields=('residual_covariance_ps2','from_average_minus_matched_ps2','mixed_minus_direct_variance_ps2')
        values=''.join(f'<td>{np.median([r[key] for r in selected]):+,.1f}</td>' for key in fields)
        covariance_rows.append(f'<tr><td>{run}</td>{values}</tr>')
    text+=('''<h3>Covariance check</h3><p>Centered union at T1 = 10, medians over DUT bars 4–11.
        The direct covariance of the corrected DUT–REF residuals agrees with the variance identity.
        The mixed-minus-direct variance is a separate difference between correction methods.</p>
        <div class="scroll"><table><thead><tr><th>Run</th><th>Cov(X,Y) [ps²]</th>
        <th>Var(A) − Var(D)/4 [ps²]</th><th>Mixed − direct variance [ps²]</th></tr></thead><tbody>'''
        +''.join(covariance_rows)+'</tbody></table></div>')
    if historical:
        note=np.median([r['note_digitized_ratio'] for r in historical])
        text+=(f'<h3>Historical note comparison</h3><p>The {len(historical)} digitized DN-26-014 '
               f'CatR-to-direct ratios have median {note:.3f}. The table uses the same DUT bars '
               'where the corresponding fit is valid in our run 4194. The note includes a different '
               'run grouping and reconstruction, so its ratios are a scale for comparison.</p>')
        examples=[]
        for sample in ('published','ref07','union','centered'):
            formatted=[]
            for selection in ('fixed_t1','optimized'):
                values=[r[f'{sample}_catr_{selection}_ratio'] for r in historical
                        if r[f'{sample}_catr_{selection}_ratio'] is not None]
                formatted.append(f'<td>{np.median(values):.3f} ({len(values)})</td>' if values else '<td>—</td>')
            examples.append(f'<tr><td>{sample}</td>'+''.join(formatted)+'</tr>')
        text+=('''<div class="scroll"><table><thead><tr><th>Sample</th><th>CatR/direct at T1=10 (N)</th>
              <th>CatR/direct at independently minimized T1 (N)</th></tr></thead><tbody>'''
              +''.join(examples)+'</tbody></table></div>')
        centered=[r['centered_catr_fixed_t1_ratio'] for r in historical
                  if r['centered_catr_fixed_t1_ratio'] is not None]
        if centered:
            text+=(f'<p>After centering, the run 4194 ratio is {np.median(centered):.3f} '
                   f'at T1 = 10, compared with {note:.3f} in the historical note. The test '
                   'therefore gives no evidence that REF membership alone produces the improvement '
                   'reported there. Pooling adds a separate effect through REF-dependent offsets. '
                   'A quantitative reproduction of the note would also require its original '
                   'event selection, corrections, fitting and threshold choices.</p>')
    return text

def historical_comparison(rows, path):
    if not path.is_file():
        return []
    with path.open() as source:
        digitized={int(r['bar']):float(r['note_digitized_ratio']) for r in csv.DictReader(source)
                   if r['note_digitized_ratio']}
    output=[]
    for bar, ratio in digitized.items():
        row=dict(bar=bar,note_digitized_ratio=ratio)
        for sample in ('published','ref07','union','centered'):
            selected=[r for r in rows if r['run']==4194 and r['bar']==bar and r['sample']==sample
                      and r['stage']=='full' and r['estimator']=='gaussian']
            fixed=next((r for r in selected if r['threshold']==10),{})
            for target in ('catr','avg_mixed'):
                target_width=fixed.get(target+'_ps')
                direct=fixed.get('direct_halfdiff_ps')
                row[f'{sample}_{target}_fixed_t1_ratio']=target_width/direct if target_width is not None and direct else None
                target_choices=[r[target+'_ps'] for r in selected if r[target+'_ps'] is not None]
                direct_choices=[r['direct_halfdiff_ps'] for r in selected if r['direct_halfdiff_ps']]
                row[f'{sample}_{target}_optimized_ratio']=min(target_choices)/min(direct_choices) if target_choices and direct_choices else None
        output.append(row)
    return output

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work',type=Path,required=True); p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--style',type=Path,required=True)
    p.add_argument('--diagnostics-production',type=Path,help='Production directory containing archived step2 PNGs')
    a=p.parse_args(); out=a.output; out.mkdir(parents=True,exist_ok=True)
    (out/'plots').mkdir(exist_ok=True); (out/'data').mkdir(exist_ok=True)
    if a.style.resolve() != (out/'style.css').resolve():
        shutil.copy2(a.style,out/'style.css')
    rows=[]; covrows=[]; exclusions=[]; statistics=[]; common_events=[]; paired=[]; offsets=[]
    def add(run,sample,g,stage,metrics):
        for est in ('moments','gaussian'):
            m=metrics.get(est); boot=metrics.get('bootstrap',{}).get(est,{})
            row=dict(run=run,module=MODULES[run],sample=sample,bar=g['bar'],threshold=g['threshold'],vov=g['vov'],
                     stage=stage,estimator=est,accepted=g.get('accepted',g.get('counts',{}).get('accepted')),
                     valid=m is not None,bootstrap_requested=boot.get('requested'),bootstrap_successful=boot.get('successful'))
            for k in METHODS:
                row[k+'_variance_ps2']=m.get(k) if m else None
                row[k+'_ps']=root_or_none(m.get(k)) if m else None
                row[k+'_error_ps']=boot.get('resolution_errors_ps',{}).get(k)
            row['mixed_minus_direct_ps']=gap(m)
            row['closure_ps2']=m.get('closure_ps2') if m else None
            rows.append(row)
    for run in RUNS:
        for ref in REF_BARS:
            folder=a.work/'production_final'/str(run)/f'ref{ref:02d}'
            if not (folder/'complete.json').exists() and not (folder/'failure.json').exists(): raise RuntimeError(f'Pending production: {folder}')
            report=folder/'closure/summary.json'; sample=f'ref{ref:02d}'
            if not report.exists():
                exclusions.append(dict(run=run,ref=ref,bar='all',threshold='all',reason='No completed closure; see stage logs'))
                continue
            d=json.loads(report.read_text())
            statistics.append(dict(run=run,ref=ref,groups=len(d['groups']),skipped=len(d['skipped']),
                                   accepted=sum(g['counts']['accepted'] for g in d['groups']),
                                   phase_flags=sum(g['counts'].get('invalid_phase_bin',0) for g in d['groups'])))
            for g in d['skipped']: exclusions.append(dict(run=run,ref=ref,bar=g['bar'],threshold=g['threshold'],reason=g['reason'],**g['counts']))
            for g in d['groups']:
                flat=a.work/'flat_final'/str(run)/sample/(g['label']+'.npz')
                with np.load(flat) as data:
                    for stage in ('tw','full'):
                        add(run,sample,g,stage,g['stages'][stage])
                        covrows.append(dict(run=run,sample=sample,bar=g['bar'],threshold=g['threshold'],stage=stage,
                                            **covariance_diagnostics(data[stage],data['ref_diff'])))
            dest=out/'plots'/str(run)/sample; dest.mkdir(parents=True,exist_ok=True)
            for f in (folder/'comparison').glob('timing_*'):
                if f.suffix in ('.png','.pdf'): shutil.copy2(f,dest/f.name)
            shutil.copy2(report,out/'data'/f'{run}_{sample}_summary.json')
        published=json.loads((a.baseline/str(run)/'closure/summary.json').read_text())
        for g in published['groups']:
            for stage in ('tw','full'): add(run,'published',g,stage,g['stages'][stage])
        dest=out/'plots'/str(run)/'published'; dest.mkdir(parents=True,exist_ok=True)
        for f in (a.baseline/str(run)/'comparison').glob('timing_*'):
            if f.suffix in ('.png','.pdf'): shutil.copy2(f,dest/f.name)
        groups=[]
        for bar in range(16):
            for th in (3,5,7,10,15,20):
                label=f'bar{bar:02d}_Vov3.00_th{th:02d}_energyBin01'
                path=a.work/'aggregate'/str(run)/label/'summary.json'
                if not path.exists(): raise RuntimeError(f'Pending aggregation: {path}')
                g=json.loads(path.read_text()); groups.append(g)
                for stage in ('tw','full'):
                    data=g['stages'][stage]
                    averages={int(ref):v[2] for ref,v in data['mean_residual_by_ref_ps'].items()}
                    offsets.append(dict(run=run,bar=bar,threshold=th,stage=stage,
                                        average_ref_mean_span_ps=max(averages.values())-min(averages.values()),
                                        **{f'average_ref{ref:02d}_mean_ps':v for ref,v in averages.items()},
                                        **{f'{q}_between_ref_variance_ps2':v for q,v in data['between_ref_variance_ps2'].items()}))
                    for ref,diagnostic in data['common_events_vs_ref7'].items():
                        common_events.append(dict(run=run,bar=bar,threshold=th,stage=stage,ref=ref,**diagnostic))
                    for estimator,values in data['difference_to_ref7'].items():
                        uncertainty=data['union']['bootstrap']['paired_difference_to_ref7'][estimator]
                        paired.append(dict(run=run,bar=bar,threshold=th,stage=stage,estimator=estimator,**values,
                                           successful=uncertainty['successful'],
                                           **{key+'_error':v for key,v in uncertainty['errors'].items()}))
                    for sample in ('union','centered','single_ref_per_event'):
                        m=g['stages'][stage][sample]
                        n=g['unique_events'] if sample=='single_ref_per_event' else g['pairs']
                        add(run,sample,dict(g,accepted=n),stage,m)
                        covrows.append(dict(run=run,sample=sample,bar=bar,threshold=th,stage=stage,**m['covariance']))
                statistics.append(dict(run=run,ref='union',bar=bar,threshold=th,accepted=g['pairs'],
                                       unique_events=g['unique_events'],shared_events=g['shared_events'],
                                       duplicated_pairs=g['duplicated_pairs'],max_references_per_event=g['max_references_per_event']))
        for sample in ('union','centered','single_ref_per_event'):
            dest=out/'plots'/str(run)/sample; dest.mkdir(parents=True,exist_ok=True)
            plot_rows=[]; plot_groups=[]
            for g in groups:
                row={k:g[k] for k in ('label','bar','threshold','vov')}
                for est in ('moments','gaussian'):
                    m=g['stages']['full'][sample][est] or {}
                    for key in METHODS: row[f'{est}_{key}_ps']=root_or_none(m.get(key))
                plot_rows.append(row); plot_groups.append(dict(g,stages={'full':g['stages']['full'][sample]}))
            plot_methods(dest,plot_rows,dict(groups=plot_groups),title_prefix=f'Run {run} · {sample} · ',previous=False)
        dump_json(out/'data'/f'{run}_aggregate.json',dict(groups=groups))
        plot_ref_comparison(out/'plots', run, rows)
    comparison_keys=('run','bar','threshold','stage','estimator')
    baseline_rows={tuple(row[k] for k in comparison_keys):row for row in rows if row['sample']=='ref07'}
    for row in rows:
        baseline_row=baseline_rows.get(tuple(row[k] for k in comparison_keys),{})
        for method in METHODS:
            width=row[method+'_ps']
            reference=baseline_row.get(method+'_ps')
            row[method+'_delta_ref7_ps']=width-reference if width is not None and reference is not None else None
    write_csv(out/'data/all_results.csv',rows); write_csv(out/'data/covariance.csv',covrows)
    write_csv(out/'data/exclusions.csv',exclusions,('run','ref','bar','threshold','reason'))
    write_csv(out/'data/statistics.csv',statistics)
    write_csv(out/'data/common_events_vs_ref7.csv',common_events)
    write_csv(out/'data/paired_union_vs_ref7.csv',paired)
    write_csv(out/'data/ref_offsets.csv',offsets)
    historical=historical_comparison(rows,a.style.parent/'note_comparison_4194.csv')
    write_csv(out/'data/historical_comparison_4194.csv',historical)
    summaries=[]
    for run in RUNS:
        for th in (3,5,7,10,15,20):
            for est in ('moments','gaussian'):
                for selection,bars in [('all_DUT_bars',range(16)),('central_DUT_4_11',range(4,12))]:
                    r=dict(run=run,threshold=th,estimator=est,selection=selection)
                    for sample in ('published','ref07','union','centered','single_ref_per_event'):
                        selected=[x for x in rows if x['run']==run and x['sample']==sample and x['threshold']==th and x['estimator']==est and x['stage']=='full' and x['bar'] in bars and x['mixed_minus_direct_ps'] is not None]
                        r[sample+'_valid_groups']=len(selected)
                        r[sample+'_median_gap_ps']=float(np.median([x['mixed_minus_direct_ps'] for x in selected])) if selected else None
                    summaries.append(r)
    write_csv(out/'data/summary_by_run_threshold.csv',summaries)
    dump_json(out/'data/results.json',rows)
    shutil.copy2(a.work/'published_validation/validation.json',out/'data/baseline_validation.json')
    ui_fields=('run','sample','bar','threshold','stage','estimator','accepted')+tuple(
        method+suffix for method in METHODS for suffix in ('_ps','_error_ps','_delta_ref7_ps'))
    visible_rows=[{key:row[key] for key in ui_fields} for row in rows if row['stage']=='full']
    (out/'data.js').write_text('const RESULTS='+json.dumps(visible_rows,allow_nan=False,separators=(',',':'))+';\n')
    table_rows=[]
    for row in summaries:
        if row['selection']!='central_DUT_4_11':
            continue
        fields=[str(row[k]) for k in ('run','threshold','estimator')]
        for sample in ('centered','ref07','published'):
            value=row[sample+'_median_gap_ps']
            fields.append(f'{value:+.2f} ({row[sample+"_valid_groups"]})' if value is not None else '— (0)')
        table_rows.append('<tr>'+''.join(f'<td>{html.escape(field)}</td>' for field in fields)+'</tr>')
    table=''.join(table_rows)
    base='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CMS BTL · Multiple REF bars</title><link rel="stylesheet" href="style.css"><script defer src="data.js"></script><script defer src="app.js?v=centered"></script><body><header><div class="shell mast"><div class="identity"><strong>CMS BTL</strong><span>September 2025 · REF selection study</span></div></div></header><main class="shell"><div class="intro"><h1>Does the REF bar change the closure?</h1><p>Four runs, sixteen DUT bars, six thresholds and REF bars 4, 5, 7, 8, 9, 10, 11 and 12. TOFHIR-only energy intercalibration.</p></div><nav><a href="#results">Results</a><a href="#comparison">REF7 comparison</a><a href="methodology.html">Methods & limitations</a><a href="commands.html">Commands & files</a><a href="../index.html">Published single-REF analysis</a></nav>
<section><h2>Baseline and interpretation</h2><p>The published REF7 replay is reproduced exactly for all 384 groups. The current step1 fixes an absent-channel index check; it selects a slightly different sample. The REF-choice test therefore compares all references to REF7 with that same corrected selection, while retaining the published REF7 as a separate control.</p><p>CatR and the difference of consistently corrected channels agree identically for exact moments. This agreement is not evidence of zero intrinsic L–R covariance. The comparisons below pool the selected REF samples after aligning their residual means for each run, DUT bar and threshold.</p><p class="small">All uncertainties hold fitted corrections fixed. Bootstrap uncertainties for the uncentered union resample unique DUT events together with all their REF coincidences. The centered comparison displays point estimates.</p></section>
<section id="results"><h2>Explore the four methods</h2><div class="filters"><label>Run<select id="run">RUN_OPTIONS</select></label><label>Sample<select id="sample"><option value="centered" selected>Union · offsets removed</option><option value="union">Union · original offsets</option><option value="ref07">REF7 · same selection</option><option value="published">Published REF7</option><option value="single_ref_per_event">One REF per event · highest energy</option>REF_OPTIONS</select></label><label>Threshold T1<select id="threshold"><option>3</option><option>5</option><option>7</option><option selected>10</option><option>15</option><option>20</option></select></label><label>Estimator<select id="estimator"><option value="gaussian">Gaussian core</option><option value="moments">Exact moments</option></select></label></div><figure><img id="methods" alt="Timing methods versus DUT bar"><figcaption>Time-walk and phase corrected · invalid radicands and failed fits are left blank.<a id="pdf">PDF</a></figcaption></figure><div class="scroll"><table><thead><tr><th>DUT bar</th><th>Entries</th><th>Direct [ps]</th><th>Mixed [ps]</th><th>Corrected difference [ps]</th><th>CatR [ps]</th></tr></thead><tbody id="values"></tbody></table></div></section>
<section id="comparison"><h2>Centered union compared with REF7</h2><figure><img id="refComparison" alt="All four timing resolutions for the centered union and REF7"><figcaption>Four methods using the selected run, threshold and estimator.<a id="refComparisonPdf">PDF</a></figcaption></figure><figure><img id="delta" alt="Mixed minus direct after REF offset removal compared with REF7"><figcaption>The lower panel subtracts the current REF7 gap from the centered union gap for each DUT bar. It shows point estimates.<a id="deltaPdf">PDF</a></figcaption></figure><p class="small">The summary below shows medians across DUT bars 4–11. Full DUT-bar, threshold and correction-stage dependence is retained in the downloadable tables. A median across bars is a descriptive statistic, not a combined resolution measurement.</p><div class="scroll"><table><thead><tr><th>Run</th><th>T1</th><th>Estimator</th><th>Centered union</th><th>Current REF7</th><th>Published REF7</th></tr></thead><tbody>SUMMARY_TABLE</tbody></table></div></section><section><h2>Download the evidence</h2><div class="files"><a href="data/all_results.csv">All methods and errors (CSV)</a><a href="data/covariance.csv">Covariance diagnostics (CSV)</a><a href="data/statistics.csv">Event counts and overlap (CSV)</a><a href="data/exclusions.csv">Excluded groups (CSV)</a><a href="data/summary_by_run_threshold.csv">Threshold summary (CSV)</a><a href="data/baseline_validation.json">Exact baseline validation</a></div></section></main><footer class="shell">Output archive: /eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_selected_TOFHIR · Existing results preserved.</footer></body></html>'''
    base=base.replace('RUN_OPTIONS',''.join(f'<option value="{r}">{r} · {MODULES[r]}</option>' for r in RUNS)).replace('REF_OPTIONS',''.join(f'<option value="ref{r:02d}">REF bar {r}</option>' for r in REF_BARS if r!=7)).replace('SUMMARY_TABLE',table)
    base=base.replace('<a href="data/covariance.csv">', '<a href="data/common_events_vs_ref7.csv">Shared-event comparison</a><a href="data/paired_union_vs_ref7.csv">Paired uncentered union vs REF7</a><a href="data/ref_offsets.csv">REF offsets and between-bar variance</a><a href="data/covariance.csv">')
    base=base.replace('<section id="results">','<section><h2>Results at T1 = 10</h2>'+
                      findings_html(summaries,statistics,offsets,historical,common_events,rows,covrows)
                      +'</section><section id="results">')
    base=base.replace('The summary below shows medians across DUT bars 4–11.',
        'The summary below shows medians across DUT bars 4–11, with the number of valid gaps in parentheses.')
    if historical:
        base=base.replace('<a href="data/baseline_validation.json">','<a href="data/historical_comparison_4194.csv">Historical note comparison (4194)</a><a href="data/baseline_validation.json">')
    base=base.replace('<tbody id="values"></tbody></table></div>',
        '<tbody id="values"></tbody></table></div><p class="small">The value in parentheses is the change in resolution from REF7 with the same selection.</p>')
    (out/'index.html').write_text(base)
    page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><link rel="stylesheet" href="style.css"><title>{title} · CMS BTL</title><header><div class="shell mast"><div class="identity"><strong>CMS BTL</strong> September 2025</div></div></header><main class="shell"><div class="intro"><h1>{title}</h1></div><nav><a href="index.html">Results</a><a href="methodology.html">Methods</a><a href="commands.html">Commands &amp; files</a></nav><section>{body}</section></main></html>'''
    (out/'methodology.html').write_text(page.format(title='Methods and limitations',body="""
<h2>Selection and corrections</h2>
<p>Step 1 requires both ends of the selected REF bar and at most three active REF bars. For each DUT threshold it fits the calibrated REF mean energy and accepts energies between 0.70 times the fitted Landau location and 950 ADC units, with the original 20 ADC fallback. The DUT active-bar veto and channel thresholds follow the run configuration. TOFHIR energy factors use the DUT bar, side, overvoltage and discriminator setting. DM0348 uses the chosen REF bar and side at 3 V and T1 = 10.</p>
<p>REF characterization applies end-specific energy walk corrections trained on coincidences with DUT bar 7. Step 2 retrains DUT time walk and phase profiles for each REF selection. Replay applies the saved functions to the same samples for all estimators. The final selection requires all four raw timing residuals within 5000 ps, four ToTs in [0,25] ns, and mean DUT energy inside the saved range. A group needs valid corrections and at least 100 accepted coincidences. Negative extracted variances and invalid Gaussian tuples appear as blanks in width columns.</p>
<h2>Estimators and covariance</h2>
<p>Let X = L − REF and Y = R − REF after end-specific corrections, A = (X + Y)/2, D = X − Y, and Q be L − R after direct energy-ratio and mean-phase corrections. The four signed variances are Var(Q)/4, Var(A) − [Var(X) + Var(Y) − Var(Q)]/2, Var(D)/4, and [Var(X) + Var(Y) − 2 Var(A)]/2. Square roots give widths in ps. Gaussian-core widths use the existing 400-bin histograms and peak ± one RMS fit window.</p>
<p>For exact moments CatR equals Var(D)/4. The mixed variance is 2 Var(Q)/4 − Var(D)/4, so its difference from the direct variance measures a change in difference corrections. Cov(X,Y) = Var(A) − Var(D)/4 includes Var(REF) − Cov(L,REF) − Cov(R,REF) in addition to DUT covariance. Subtracting Var((REF<sub>L</sub> − REF<sub>R</sub>)/2) gives an exploratory proxy under independent REF ends and negligible DUT–REF covariance. Those assumptions are unverified. In pooled samples this REF half-difference variance also contains offsets between REF bars, including in the centered sensitivity check, which centers only the five reconstructed DUT observables. Gaussian fits use distinct windows, so the moment identities need not hold for their widths.</p>
<h2>Aggregation and statistical errors</h2>
<p>The union contributes one row per accepted DUT–REF coincidence with the corrections trained for that REF bar. The DUT event key uses both raw integer timestamps within a run, bar and threshold. Matches across REF samples must also share DUT energies and clock phases. Bootstrap resampling draws DUT event keys with replacement, keeping all of each event's coincidences together and pairing the union with REF7. It holds fitted corrections fixed. The shared-event table compares separately trained corrections on identical selected events for REF pairs with at least 100 common observations.</p>
<p>For the centered comparison, each of the five corrected DUT observables is shifted from its REF sample mean to the pooled mean before computing the resolutions. This is done independently for each run, DUT bar, threshold and correction stage. The centered and one-REF-per-event results provide point estimates. The one-REF-per-event sample uses the highest calibrated REF mean energy.</p>
"""))
    (out/'commands.html').write_text(page.format(title='Commands and files',body="""
<p>See PROGRESS.md in the EOS archive for execution paths, executable build and checkpoint history.</p>
<pre>python3 scripts/validate_multiref_baseline.py --help
python3 scripts/audit_ref_occupancy.py --help
python3 scripts/analyze_multiref.py --help
python3 scripts/extract_multiref_events.py --help
python3 main/aggregate_multiref.py --help
python3 scripts/report_multiref.py --help</pre>
<p>ROOT intermediates, provenance, logs, calibration configurations and analysis scripts are in <code>/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_selected_TOFHIR</code>.</p>
"""))
    (out/'app.js').write_text('''const ids=['run','sample','threshold','estimator'];
function render(){const [run,sample,threshold,estimator]=ids.map(x=>document.getElementById(x).value); const stem=`plots/${run}/${sample}/timing_${estimator}_Vov3.00_th${threshold.padStart(2,'0')}`;
document.getElementById('methods').src=stem+'.png';document.getElementById('pdf').href=stem+'.pdf';const d=`plots/delta_${run}_${threshold}_${estimator}`;document.getElementById('delta').src=d+'.png?v=centered';document.getElementById('deltaPdf').href=d+'.pdf?v=centered';
const comparison=`plots/methods_centered_vs_ref7_${run}_${threshold}_${estimator}`;
document.getElementById('refComparison').src=comparison+'.png';
document.getElementById('refComparisonPdf').href=comparison+'.pdf';
const rows=RESULTS.filter(r=>r.run==run&&r.sample===sample&&r.threshold==threshold&&r.estimator===estimator&&r.stage==='full').sort((a,b)=>a.bar-b.bar);
const f=v=>v===null||v===undefined?'—':Number(v).toFixed(2);
const delta=v=>v==null?'':` (${v>=0?'+':''}${f(v)})`;
document.getElementById('values').innerHTML=rows.map(r=>'<tr><td>'+r.bar+'</td><td>'+r.accepted+'</td>'+['direct_halfdiff','avg_mixed','matched_halfdiff','catr'].map(k=>'<td>'+f(r[k+'_ps'])+(r[k+'_error_ps']!=null?' ± '+f(r[k+'_error_ps']):'')+delta(r[k+'_delta_ref7_ps'])+'</td>').join('')+'</tr>').join('');}
ids.forEach(x=>document.getElementById(x).addEventListener('change',render));render();''')
    diagnostics=a.diagnostics_production or a.baseline.parent/'timing_closure_multiref_selected_TOFHIR/production_final'
    package_all(diagnostics, out/'diagnostic-images', out/'data')
    write_pages(out, out)
    print('Report written:',out)

if __name__=='__main__': main()
