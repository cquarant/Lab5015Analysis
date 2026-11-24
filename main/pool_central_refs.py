#!/usr/bin/env python3
"""Fit the union of separately calibrated DUT–REF coincidences for each DUT setting."""
import argparse
from contextlib import ExitStack
import csv
import json
from pathlib import Path
import sys
import tempfile

import matplotlib
matplotlib.use('Agg')
import numpy as np

from moduleCharacterization_timingClosure import QUANTITIES, dump_json, open_root, require, summarize_stage

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'macros'))
from compare_timingClosure import plot_methods, root_or_none


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("production", type=Path)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--fit-bootstrap", type=int, default=100)
    args = parser.parse_args()
    if min(args.bootstrap, args.fit_bootstrap) < 0:
        parser.error("Bootstrap counts must be nonnegative")
    import ROOT
    ROOT.gROOT.SetBatch(True)
    production = args.production.resolve()
    refs = [4, 5, 7, 8, 9, 10, 11, 12]
    stages = ('tw', 'full')
    labels = {}
    rows = []
    report = dict(reference_bars=refs, excluded_reference_bar=6, groups=[], source_uuids={},
                  bootstrap=args.bootstrap, fit_bootstrap=args.fit_bootstrap, root_version=ROOT.gROOT.GetVersion(),
                  pooling="One entry per accepted DUT–REF pair, with its own calibrated residuals. "
                          "No further centering. Bootstrap unit is one coincidence.")
    with tempfile.TemporaryDirectory(dir=production) as tmp, ExitStack() as stack:
        out = Path(tmp) / 'pooled'
        out.mkdir()
        for ref in refs:
            folder = production / f'ref{ref:02d}' / 'closure'
            single = json.loads((folder/'summary.json').read_text())
            source = open_root(ROOT, folder/'timing_closure.root')
            stack.callback(source.Close)
            report['source_uuids'][ref] = source.GetUUID().AsString()
            for g in single['groups'] + single['skipped']:
                if g['counts'].get('accepted', 0):
                    labels.setdefault(g['label'], {})[ref] = (g, source)
        destination = ROOT.TFile(str(out/'timing_closure.root'), 'CREATE')
        if destination.IsZombie():
            raise RuntimeError("Cannot create pooled ROOT file")
        with destination:
            for label, members in sorted(labels.items()):
                parts = {stage: [] for stage in stages}
                counts = {r: 0 for r in refs}
                means = {}
                for ref, (g, source) in members.items():
                    tree = require(source, 'events_'+label, 'TTree')
                    columns = [s+'_'+q for s in stages for q in QUANTITIES]
                    data = ROOT.RDataFrame(tree).AsNumpy(columns)
                    counts[ref] = int(tree.GetEntries())
                    if counts[ref] != g['counts']['accepted']:
                        raise RuntimeError(f'Incomplete REF {ref} sample for {label}')
                    means[ref] = {}
                    for stage in stages:
                        values = np.column_stack([data[stage+'_'+q] for q in QUANTITIES])
                        parts[stage].append(values)
                        means[ref][stage] = dict(zip(QUANTITIES, values.mean(axis=0)))
                group = {k: g[k] for k in ('label', 'bar', 'threshold', 'vov')}
                group.update(counts_by_ref=counts, mean_residual_ps_by_ref=means, accepted=sum(counts.values()), stages={})
                if group['accepted'] < 2:
                    report.setdefault('skipped', []).append(group)
                    continue
                row = {k: group[k] for k in ('label','bar','threshold','vov','accepted')}
                row.update({f'ref{r:02d}_count': counts[r] for r in refs})
                for stage, arrays in parts.items():
                    values = np.concatenate(arrays)
                    if not np.isfinite(values).all() or len(values) != group['accepted']:
                        raise RuntimeError(f'Invalid pooled sample for {stage} {label}')
                    destination.cd()
                    result = summarize_stage(ROOT, values, stage, label, args.bootstrap, args.fit_bootstrap,
                                             20261116 + group['bar'] + 100*group['threshold'])
                    group['stages'][stage] = result
                    for estimator in ('gaussian', 'moments'):
                        metrics = result[estimator] or {}
                        for key in ('direct_halfdiff','matched_halfdiff','catr','avg_mixed','avg_matched'):
                            prefix = '' if stage == 'full' else 'tw_'
                            row[f'{prefix}{estimator}_{key}_ps'] = root_or_none(metrics.get(key))
                rows.append(row)
                report['groups'].append(group)
                print(label, group['accepted'], flush=True)
        if not rows:
            raise RuntimeError("No pooled DUT groups")
        dump_json(out/'summary.json', report)
        with (out/'comparison.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        plot_methods(out, rows, report, title_prefix='REF 4–12 excluding 6 · ', previous=False)
        (out/'complete').touch()
        out.rename(production/'pooled')


if __name__ == '__main__':
    main()
