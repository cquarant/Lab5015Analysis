#!/usr/bin/env python3
"""Plot a per-REF replay even when some original saved Gaussian fits are absent."""
import argparse
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'macros'))
from compare_timingClosure import root_or_none,plot_methods

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('summary',type=Path); p.add_argument('output',type=Path); a=p.parse_args()
    report=json.loads(a.summary.read_text()); rows=[]
    for g in report['groups']:
        row={k:g[k] for k in ('label','bar','threshold','vov')}; row['accepted']=g['counts']['accepted']
        for est in ('moments','gaussian'):
            m=g['stages']['full'][est] or {}
            for key in ('direct_halfdiff','avg_mixed','matched_halfdiff','catr','avg_matched'):
                row[f'{est}_{key}_ps']=root_or_none(m.get(key))
        rows.append(row)
    a.output.mkdir(parents=True,exist_ok=False)
    if rows:
        with (a.output/'comparison.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        plot_methods(a.output,rows,report,previous=False)

if __name__=='__main__': main()
