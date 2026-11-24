#!/usr/bin/env python3
"""Checkpoint flat coincidence arrays with raw DUT identity and REF timing diagnostics."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'main'))
from moduleCharacterization_timingClosure import QUANTITIES,event_dictionary,open_root,require
from multiref_statistics import REF_BARS,cluster_ids

def extract(folder,output,ROOT):
    summary=json.loads((folder/'closure/summary.json').read_text())
    output.mkdir(parents=True,exist_ok=True)
    source=open_root(ROOT,folder/'refCharacterization.root'); replay=open_root(ROOT,folder/'closure/timing_closure.root')
    manifest=dict(reference_uuid=source.GetUUID().AsString(),closure_uuid=replay.GetUUID().AsString(),groups=[])
    corrections=open_root(ROOT,folder/'step2.root')
    actual=dict(reference=manifest['reference_uuid'],corrections=corrections.GetUUID().AsString())
    corrections.Close()
    if actual!=summary['input_uuids']:
        raise RuntimeError(f'Closure inputs differ from replay provenance in {folder}: {actual}')
    for g in summary['groups']+summary['skipped']:
        if g['counts'].get('accepted',0)<2 or 'source_entries_sha256' not in g: continue
        dst=output/(g['label']+'.npz')
        if dst.exists():
            with np.load(dst) as cached:
                if str(cached['reference_uuid'])!=manifest['reference_uuid'] or str(cached['closure_uuid'])!=manifest['closure_uuid']:
                    raise RuntimeError('Cached event inputs changed')
            continue
        events=require(replay,'events_'+g['label'],'TTree')
        columns=['source_entry','invalid_phase_bin']+[s+'_'+q for s in ('tw','full') for q in QUANTITIES]
        data=ROOT.RDataFrame(events).AsNumpy(columns)
        entries=data['source_entry']
        if len(entries)!=g['counts'].get('accepted',0): raise RuntimeError('Replay count mismatch')
        if hashlib.sha256(np.asarray(entries,dtype='<i8').tobytes()).hexdigest()!=g['source_entries_sha256']: raise RuntimeError('Replay entry hash mismatch')
        if len(entries)<2: continue
        tree=require(source,g['tree'],'TTree'); df=ROOT.RDataFrame(tree)
        definitions={'dutL':'event.timeL','dutR':'event.timeR','energyL':'event.energyL','energyR':'event.energyR',
                     'phaseL':'event.t1fineL','phaseR':'event.t1fineR',
                     'ref_diff':'double(event.timeL_ref_cor-event.timeR_ref_cor)',
                     'ref_energy':'0.5*(double(event.energyL_ref)+double(event.energyR_ref))'}
        for name,expr in definitions.items(): df=df.Define(name,expr)
        raw=df.AsNumpy(list(definitions)); raw={k:v[entries] for k,v in raw.items()}
        keys=np.column_stack([raw.pop('dutL'),raw.pop('dutR')]); _,counts=cluster_ids(keys)
        if np.any(counts!=1): raise RuntimeError(f'Ambiguous duplicated raw timestamp pair in {g["label"]}')
        values={s:np.column_stack([data[s+'_'+q] for q in QUANTITIES]) for s in ('tw','full')}
        for s,v in values.items():
            if np.max(np.abs(v[:,2]-(v[:,0]+v[:,1])/2))>1e-8 or np.max(np.abs(v[:,3]-(v[:,0]-v[:,1])))>1e-8:
                raise RuntimeError('Inconsistent corrected channel identities')
        temporary=dst.with_suffix('.partial.npz')
        np.savez_compressed(temporary,keys=keys,**raw,**values,invalid_phase=data['invalid_phase_bin'],
                            reference_uuid=manifest['reference_uuid'],closure_uuid=manifest['closure_uuid'])
        temporary.replace(dst)
    source.Close(); replay.Close()
    (output/'complete.json').write_text(json.dumps(manifest,indent=2)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('production',type=Path); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--published',action='store_true',help='Read run/closure rather than run/refXX/closure')
    p.add_argument('--runs',type=int,nargs='+',default=[4194,4257,4281,4291])
    a=p.parse_args(); import ROOT; ROOT.gROOT.SetBatch(True)
    with tempfile.TemporaryDirectory(dir='/tmp',prefix='multiref-dict-') as cache:
        event_dictionary(ROOT,Path(__file__).resolve().parents[1],cache)
        for run in a.runs:
            folders=[a.production/str(run)] if a.published else [a.production/str(run)/f'ref{b:02d}' for b in REF_BARS]
            for folder in folders:
                if not (folder/'closure/summary.json').exists(): continue
                output=a.output/str(run)/('published' if a.published else folder.name)
                if (output/'complete.json').exists(): continue
                extract(folder,output,ROOT); print(run,folder.name,'extracted',flush=True)

if __name__=='__main__': main()
