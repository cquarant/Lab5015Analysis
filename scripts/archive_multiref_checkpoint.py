#!/usr/bin/env python3
"""Archive only completed stages (never a live ROOT file), without deleting old outputs."""
import argparse
from pathlib import Path
import shutil

def copy(source,destination):
    if not source.is_file(): return
    if destination.exists() and destination.stat().st_size==source.stat().st_size and destination.stat().st_mtime_ns==source.stat().st_mtime_ns: return
    destination.parent.mkdir(parents=True,exist_ok=True)
    partial=destination.with_name(destination.name+'.partial')
    shutil.copy2(source,partial); partial.replace(destination)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path); p.add_argument('destination',type=Path); a=p.parse_args()
    for stage in ('production','production_final'):
        for run in (4194,4257,4281,4291):
            for ref in range(16):
                src=a.work/stage/str(run)/f'ref{ref:02d}'; dst=a.destination/stage/str(run)/f'ref{ref:02d}'
                if not src.is_dir(): continue
                for name in ('withref.cfg','import.json','import_rejected.json','failure.json'):
                    copy(src/name,dst/name)
                for exe,product in [('moduleCharacterization_step1','step1.root'),('refCharacterization','refCharacterization.root'),('moduleCharacterization_withRef_step2','step2.root')]:
                    if not (src/(exe+'.done.json')).exists(): continue
                    for name in (product,exe+'.log',exe+'.done.json'): copy(src/name,dst/name)
                if (src/'closure/summary.json').exists():
                    for name in ('summary.json','timing_closure.root'): copy(src/'closure'/name,dst/'closure'/name)
                    copy(src/'closure.log',dst/'closure.log')
                if (src/'complete.json').exists():
                    for f in (src/'comparison').iterdir(): copy(f,dst/'comparison'/f.name)
                    copy(src/'complete.json',dst/'complete.json')
        for name in ('production.json','status.json'): copy(a.work/stage/name,a.destination/stage/name)
    for section in ('published_validation','baseline_validation','occupancy_verified'):
        src=a.work/section
        if not src.is_dir(): continue
        for f in src.iterdir():
            if f.is_file(): copy(f,a.destination/section/f.name)
            elif f.is_dir() and (f/'summary.json').exists(): copy(f/'summary.json',a.destination/section/f.name/'summary.json')
    for run in (4194,4257,4281,4291):
        for group in (a.work/'aggregate'/str(run)).glob('bar*_energyBin01'):
            if not (group/'summary.json').is_file(): continue
            for name in ('histograms.root','summary.json'):
                copy(group/name,a.destination/'aggregate'/str(run)/group.name/name)
    for name in ('resume.log','aggregation.log'):
        copy(a.work/name,a.destination/name)
    for log in (a.work/'logs').glob('*.log'):
        copy(log,a.destination/'logs'/log.name)
    repo=Path(__file__).resolve().parents[1]
    for name in ('PROGRESS.md','TIMING_CLOSURE.md'):
        copy(repo/name,a.destination/name)
    scripts=['scripts/analyze_multiref.py','scripts/analyze_central_refs.py','scripts/validate_multiref_baseline.py',
             'scripts/audit_ref_occupancy.py','scripts/extract_multiref_events.py','scripts/archive_multiref_checkpoint.py',
             'main/moduleCharacterization_timingClosure.py','main/moduleCharacterization_withRef_step2.cpp',
             'main/multiref_statistics.py','main/aggregate_multiref.py','main/plot_multiref_single.py',
             'macros/compare_timingClosure.py','tests/test_multiref_statistics.py','tests/test_timing_closure.py',
             'scripts/report_multiref.py','scripts/verify_multiref_selected.py',
             'scripts/restore_multiref_inputs.py','scripts/resume_multiref_selected.sh']
    scripts.extend(['scripts/organize_multiref_pages.py','scripts/package_multiref_diagnostics.py'])
    scripts.extend(str(path.relative_to(repo)) for path in (repo/'scripts/multiref_web').iterdir() if path.is_file())
    for name in scripts: copy(repo/name,a.destination/'source'/name)
    print('Checkpoint archived to',a.destination)

if __name__=='__main__': main()
