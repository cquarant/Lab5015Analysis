#!/usr/bin/env python3
"""Replay independently rerun REF7 corrections and compare with published baseline."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--production', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--published-inputs', action='store_true', help='Validate replay on the exact published inputs')
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    def check(run):
        dst = a.output/str(run)
        new = a.baseline/str(run) if a.published_inputs else a.production/str(run)/'ref07'
        if not (dst/'summary.json').exists():
            with (a.output/f'{run}.log').open('w') as log:
                subprocess.run([sys.executable, REPO/'main/moduleCharacterization_timingClosure.py',
                    '--reference', new/'refCharacterization.root', '--corrections', new/'step2.root',
                    '--output-dir', dst, '--bootstrap', '0', '--fit-bootstrap', '0'],
                    stdout=log, stderr=subprocess.STDOUT, check=True)
        old = json.loads((a.baseline/str(run)/'closure/summary.json').read_text())
        replay = json.loads((dst/'summary.json').read_text())
        previous = {g['label']:g for g in old['groups']}
        result = dict(run=run, groups=len(replay['groups']), failures=[], max_variance_delta_ps2=0.,
                      max_gaussian_width_delta_ps=0., accepted=0)
        if set(previous) != {g['label'] for g in replay['groups']}:
            result['failures'].append('group membership differs')
        for g in replay['groups']:
            b = previous[g['label']]
            result['accepted'] += g['counts']['accepted']
            for k in ('counts','source_entries_sha256'):
                if g[k] != b[k]:
                    result['failures'].append(g['label']+': '+k)
            for stage in ('tw','full'):
                for k, v in g['stages'][stage]['moments'].items():
                    delta = abs(v-b['stages'][stage]['moments'][k])
                    result['max_variance_delta_ps2'] = max(result['max_variance_delta_ps2'], delta)
                for k, v in g['stages'][stage]['fits'].items():
                    ov = b['stages'][stage]['fits'][k]['sigma']
                    if v['sigma'] is None or ov is None:
                        if v['sigma'] != ov:
                            result['failures'].append(g['label']+': fit validity '+k)
                    else:
                        result['max_gaussian_width_delta_ps'] = max(result['max_gaussian_width_delta_ps'], abs(v['sigma']-ov))
        result['passed'] = (not result['failures'] and result['groups']==96 and
                            result['max_variance_delta_ps2'] < 1e-6 and result['max_gaussian_width_delta_ps'] < 1e-4)
        print(json.dumps(result), flush=True)
        return result
    with ThreadPoolExecutor(a.jobs) as pool:
        results = list(pool.map(check, (4194,4257,4281,4291)))
    (a.output/'validation.json').write_text(json.dumps(results, indent=2)+'\n')
    if not all(r['passed'] for r in results):
        raise SystemExit('Baseline validation failed; inspect validation.json before production')


if __name__ == '__main__':
    main()
