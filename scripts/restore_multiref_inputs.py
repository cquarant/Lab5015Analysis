#!/usr/bin/env python3
"""Restore the exact archived inputs for a completed replay, preserving alternatives."""
import argparse
import json
from pathlib import Path
import shutil


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('archive',type=Path)
    parser.add_argument('run',type=int)
    parser.add_argument('refs',type=int,nargs='+')
    args=parser.parse_args()
    import ROOT
    def uuid(path):
        root=ROOT.TFile.Open(str(path),'READ')
        if not root or root.IsZombie():
            raise RuntimeError(f'Invalid ROOT file: {path}')
        result=root.GetUUID().AsString()
        root.Close()
        return result
    quarantine=args.work/'before_replay_input_recovery'/str(args.run)
    quarantine.mkdir(parents=True,exist_ok=True)
    changed=False
    for ref in args.refs:
        name=f'ref{ref:02d}'
        target=args.work/'production_final'/str(args.run)/name
        source=args.archive/str(args.run)/name
        summary=json.loads((target/'closure/summary.json').read_text())
        receipt=json.loads((target/'complete.json').read_text())
        orig=quarantine/name
        updated={}
        for key,filename in [('reference','refCharacterization.root'),('corrections','step2.root')]:
            want=summary['input_uuids'][key]
            if uuid(target/filename)==want:
                continue
            if uuid(source/filename)!=want:
                raise RuntimeError(f'Archive UUID is incompatible: {source/filename}')
            orig.mkdir(exist_ok=True)
            if not (orig/filename).exists():
                shutil.copy2(target/filename,orig/filename)
            temp=target/(filename+'.restore')
            shutil.copy2(source/filename,temp)
            temp.replace(target/filename)
            receipt['sizes'][filename]=(target/filename).stat().st_size
            updated[filename]=dict(source=str(source/filename),root_uuid=want)
        if updated:
            changed=True
            receipt['restored_replay_inputs']=updated
            temp=target/'complete.restore.json'
            temp.write_text(json.dumps(receipt,indent=2)+'\n')
            temp.replace(target/'complete.json')
            flat=args.work/'flat_final'/str(args.run)/name
            if flat.is_dir():
                flat.rename(orig/'flat')
            print(args.run,ref,'restored',','.join(updated),flush=True)
    if changed:
        aggregate=args.work/'aggregate'/str(args.run)
        if aggregate.is_dir():
            aggregate.rename(quarantine/'aggregate')

if __name__=='__main__':
    main()
