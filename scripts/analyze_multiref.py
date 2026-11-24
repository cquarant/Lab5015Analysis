#!/usr/bin/env python3
"""Resumable, separate per-REF production using the existing three executables and replay.

Requires a successful published-baseline validation. No REF is excluded implicitly.
Runs at most --jobs executable processes; keeps intermediates and a per-stage receipt.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import fcntl
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from analyze_central_refs import option, run_logged

REPO=Path(__file__).resolve().parents[1]
RUNS=(4194,4257,4281,4291)
EXES=('moduleCharacterization_step1','refCharacterization','moduleCharacterization_withRef_step2')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def atomic_json(path,value):
    tmp=path.with_suffix(path.suffix+'.partial')
    tmp.write_text(json.dumps(value,indent=2)+'\n'); tmp.replace(path)

def settings(cfg):
    ignored={('Input','step1FileName'),('Output','outFileNameStep1'),('Output','outFileNameStep2'),
             ('Output','refFileName'),('Output','plotDir')}
    for section,key in ignored: cfg=option(cfg,section,key,'OUTPUT')
    for section,key in [('Input','inputDir')]: cfg=option(cfg,section,key,option(cfg,section,key).rstrip('/'))
    return cfg.strip()

def verified_relocated_input(checkpoint, previous, current, spills, run):
    # The previous central REF4/5 production used a local cache of the same
    # raw spills. Check its complete input list and processed count.
    raw_path=option(previous,'Input','inputDir')
    if not raw_path.startswith('/tmp/btl-central-ref.') or not raw_path.endswith('/reco'):
        return None
    if settings(option(previous,'Input','inputDir',option(current,'Input','inputDir')))!=settings(current):
        return None
    log=checkpoint/'moduleCharacterization_step1.log'
    if not log.is_file():
        return None
    content=log.read_text()
    added=Counter(re.findall(r'Adding file: (\S+)',content))
    expected=[row for row in spills if int(row['run'])==run and
              float(row['vov_first_entry'])==float(option(current,'Input','vov'))]
    entries=sum(int(row['entries']) for row in expected)
    tally=re.search(r'external bar loop: reading entry 0 / (\d+)',content)
    if not expected or not tally or int(tally[1])!=entries:
        return None
    if added!=Counter(Path(row['reco_path']).name for row in expected):
        return None
    return dict(raw_cache=raw_path,spill_count=len(expected),raw_entries=entries)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--validation',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bin-dir',type=Path,required=True)
    p.add_argument('--reuse',type=Path,action='append',default=[],
                   help='Try checkpoint directories in the given order')
    p.add_argument('--reuse-work',type=Path)
    p.add_argument('--runs',type=int,nargs='+',default=list(RUNS))
    p.add_argument('--ref-bars',type=int,nargs='+',default=[4,5,7,8,9,10,11,12])
    p.add_argument('--scratch',type=Path,help='Disposable ROOT and plotting caches')
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--bootstrap',type=int,default=200)
    p.add_argument('--fit-bootstrap',type=int,default=100)
    a=p.parse_args()
    if not 1<=a.jobs<=4 or any(b not in range(16) for b in a.ref_bars): p.error('Invalid jobs or reference bars')
    verified=json.loads(a.validation.read_text())
    if not all(any(v['run']==r and v['passed'] for v in verified) for r in a.runs): p.error('Baseline not verified')
    out=a.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    lock=(out/'.lock').open('a'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    for key,name in [('TMPDIR','tmp'),('MPLCONFIGDIR','mpl'),('CCACHE_DIR','ccache'),('CCACHE_TEMPDIR','ccache-tmp')]:
        local=(a.scratch.resolve() if a.scratch else out)/name
        local.mkdir(parents=True,exist_ok=True); os.environ[key]=str(local)
    code=['main/moduleCharacterization_timingClosure.py','main/moduleCharacterization_step1.cpp',
          'main/refCharacterization.cpp','main/moduleCharacterization_withRef_step2.cpp','src/referenceCharacterization_helper.cc']
    manifest=dict(code={f:sha(REPO/f) for f in code},bootstrap=a.bootstrap,fit_bootstrap=a.fit_bootstrap,
                  binaries={f:sha(a.bin_dir/(f+'.exe')) for f in EXES},library=sha(REPO/'lib/libTBStudies.so'),
                  baseline=str(a.baseline.resolve()))
    previous=out/'production.json'
    if previous.exists() and json.loads(previous.read_text())!=manifest: p.error('Production fingerprint changed; use a new output')
    atomic_json(previous,manifest)
    boot=['--bootstrap',str(a.bootstrap),'--fit-bootstrap',str(a.fit_bootstrap)]
    spill_manifest=a.baseline/'input_spills.csv'
    spills=list(csv.DictReader(spill_manifest.open())) if spill_manifest.is_file() else []
    def analyze(item):
        run,ref=item; folder=out/str(run)/f'ref{ref:02d}'; folder.mkdir(parents=True,exist_ok=True)
        cfg=(a.baseline/str(run)/'withref.cfg').read_text()
        mapping=list(map(int,option(cfg,'Channels','channelMapping').split())); asic=int(option(cfg,'Channels','REFasic'))
        for key,ch in zip(('chL','chR'),mapping[2*ref:2*ref+2]): cfg=option(cfg,'Coincidence',key,32*asic+ch)
        for section,key,name in [('Input','step1FileName','step1.root'),('Output','outFileNameStep1','step1.root'),
                                ('Output','refFileName','refCharacterization.root'),('Output','outFileNameStep2','step2.root'),
                                ('Output','plotDir','diagnostics')]: cfg=option(cfg,section,key,folder/name)
        config=folder/'withref.cfg'
        if config.exists() and config.read_text()!=cfg: raise RuntimeError(f'Configuration changed: {config}')
        config.write_text(cfg)
        if (folder/'complete.json').exists():
            receipt=json.loads((folder/'complete.json').read_text())
            if all((folder/k).is_file() and (folder/k).stat().st_size==v for k,v in receipt['sizes'].items()): return receipt
            raise RuntimeError(f'Incomplete completed output: {folder}')
        try:
            src=None
            relocation=None
            for candidate in a.reuse:
                checkpoint=candidate/str(run)/f'ref{ref:02d}'
                if not (checkpoint/'withref.cfg').is_file():
                    continue
                previous_config=(checkpoint/'withref.cfg').read_text()
                if settings(previous_config)!=settings(cfg):
                    relocation=verified_relocated_input(checkpoint,previous_config,cfg,spills,run)
                    if relocation is None:
                        atomic_json(folder/'import_rejected.json',dict(source=str(checkpoint),reason='Configuration/input path differs'))
                        continue
                if any((checkpoint/name).is_file() for name in ('step1.root','refCharacterization.root','step2.root')):
                    src=checkpoint
                    break
            if src and (src/'withref.cfg').exists() and not (folder/'import.json').exists():
                imports={}
                if relocation: imports['relocated_input_validation']=relocation
                for name in ('refCharacterization.root','step2.root','moduleCharacterization_step1.log',
                             'refCharacterization.log','moduleCharacterization_withRef_step2.log'):
                    if (src/name).exists(): shutil.copy2(src/name,folder/name); imports[name]=dict(source=str(src/name),sha256=sha(folder/name))
                legacy_work=(a.reuse_work/str(run)/f'ref{ref:02d}'
                             if a.reuse_work and src.parent.parent==a.reuse[-1] else None)
                if legacy_work or (src/'step1.root').exists():
                    step=src/'step1.root'
                    if legacy_work and not step.exists():
                        step=legacy_work/'step1.root'
                    if step.exists(): shutil.copy2(step,folder/'step1.root'); imports['step1.root']=dict(source=str(step),sha256=sha(step))
                # Reuse only genuinely completed closure, including its original bootstraps.
                summary=src/'closure/summary.json'
                if summary.exists():
                    old=json.loads(summary.read_text())
                    if old['arguments']['bootstrap']==a.bootstrap and old['arguments']['fit_bootstrap']==a.fit_bootstrap:
                        for name in ('closure','comparison'):
                            if (src/name).is_dir() and not (folder/name).exists(): shutil.copytree(src/name,folder/name)
                atomic_json(folder/'import.json',imports)
                for exe,product in zip(EXES,('step1.root','refCharacterization.root','step2.root')):
                    old_done=(src/(exe+'.done.json')).exists()
                    if exe!=EXES[0] and (src/'complete').exists():
                        old_done=True
                    if legacy_work:
                        old_done |= (legacy_work/(exe+'.done')).exists()
                    if old_done and (folder/product).is_file(): atomic_json(folder/(exe+'.done.json'),dict(imported=True))
            for exe in EXES:
                marker=folder/(exe+'.done.json')
                if marker.exists(): continue
                if exe==EXES[2] and (folder/(exe+'.log')).exists():
                    run_logged([a.bin_dir.resolve()/(EXES[1]+'.exe'),config],folder,folder/(EXES[1]+'.log'))
                start=time.time(); command=[a.bin_dir.resolve()/(exe+'.exe'),config]
                if exe==EXES[2]: command.append('--skip-uncalibratable')
                run_logged(command,folder,folder/(exe+'.log'))
                atomic_json(marker,dict(seconds=time.time()-start))
            with tempfile.TemporaryDirectory(dir=folder,prefix='replay-') as temp:
                work=Path(temp)
                if not (folder/'closure/summary.json').exists():
                    run_logged([sys.executable,REPO/'main/moduleCharacterization_timingClosure.py',
                                '--reference',folder/'refCharacterization.root','--corrections',folder/'step2.root',
                                '--output-dir',work/'closure','--skip-unusable',*boot],folder,folder/'closure.log')
                    (work/'closure').rename(folder/'closure')
                if not (folder/'comparison/comparison.csv').exists():
                    run_logged([sys.executable,REPO/'main/plot_multiref_single.py',folder/'closure/summary.json',
                                work/'comparison'],folder,folder/'comparison.log')
                    (work/'comparison').rename(folder/'comparison')
            summary=json.loads((folder/'closure/summary.json').read_text())
            files=['step1.root','refCharacterization.root','step2.root','closure/summary.json','closure/timing_closure.root','comparison/comparison.csv']
            result=dict(run=run,ref=ref,groups=len(summary['groups']),skipped=len(summary['skipped']),
                        accepted=sum(g['counts']['accepted'] for g in summary['groups']),
                        sizes={f:(folder/f).stat().st_size for f in files})
            atomic_json(folder/'complete.json',result)
            (folder/'failure.json').unlink(missing_ok=True)
            print(json.dumps(result),flush=True); return result
        except Exception as error:
            result=dict(run=run,ref=ref,error=str(error)); atomic_json(folder/'failure.json',result)
            print(json.dumps(result),flush=True); return result
    with ThreadPoolExecutor(a.jobs) as pool:
        results=list(pool.map(analyze,[(r,b) for r in a.runs for b in a.ref_bars]))
    atomic_json(out/'status.json',results)
    if any('error' in r for r in results): raise SystemExit('Some references failed; inspect status.json and per-REF logs')

if __name__=='__main__': main()
