#!/usr/bin/env python3
"""Pool calibrated coincidences; preserve REF membership and bootstrap raw DUT events."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import numpy as np
from moduleCharacterization_timingClosure import QUANTITIES,dump_json,summarize_stage,gauss_fit,gaussian_metrics,moments
from multiref_statistics import (METHODS,REF_BARS,cluster_ids,choose_single_reference,center_by_ref,
                                covariance_diagnostics,weighted_moments,error_summary)

def differences(a,b):
    out={k+'_variance_ps2':a[k]-b[k] for k in METHODS}
    for k in METHODS:
        out[k+'_resolution_ps']=float(np.sqrt(a[k])-np.sqrt(b[k])) if min(a[k],b[k])>=0 else float('nan')
    for prefix,m in [('sample',a),('ref7',b)]:
        out[prefix+'_mixed_direct_gap_ps']=float(np.sqrt(m['avg_mixed'])-np.sqrt(m['direct_halfdiff'])) if min(m['avg_mixed'],m['direct_halfdiff'])>=0 else float('nan')
    out['mixed_direct_gap_change_ps']=out['sample_mixed_direct_gap_ps']-out['ref7_mixed_direct_gap_ps']
    return out

def errors(rows):
    if len(rows)<2: return {}
    return {k:float(np.std([x[k] for x in rows],ddof=1)) for k in rows[0]
            if np.isfinite([x[k] for x in rows]).all()}

def cluster_bootstrap(ROOT,values,ids,refs,n_moments,n_fits,seed):
    rng=np.random.default_rng(seed); nclusters=int(ids.max())+1
    subsets={int(r):np.flatnonzero(refs==r) for r in np.unique(refs)}
    pooled_moments=[]; pooled_fits=[]; delta_moments=[]; delta_fits=[]
    per_ref={r:[] for r in subsets}; base_idx=subsets.get(7)
    if base_idx is None: raise RuntimeError('No safe REF7 baseline in group')
    for i in range(max(n_moments,n_fits)):
        cluster_weight=np.bincount(rng.integers(nclusters,size=nclusters),minlength=nclusters)
        weight=cluster_weight[ids]
        if i<n_moments:
            pooled=weighted_moments(values,weight); base=weighted_moments(values[base_idx],weight[base_idx])
            if pooled is not None and base is not None:
                pooled_moments.append(pooled); delta_moments.append(differences(pooled,base))
                for ref,idx in subsets.items():
                    metric=weighted_moments(values[idx],weight[idx])
                    if metric is not None: per_ref[ref].append(differences(metric,base))
        if i<n_fits:
            sample=values[np.repeat(np.arange(len(values)),weight)]
            fits={q:gauss_fit(ROOT,sample[:,j],'boot_union_'+q) for j,q in enumerate(QUANTITIES)}
            metric=gaussian_metrics(fits)
            if metric is not None: pooled_fits.append(metric)
            b=values[base_idx][np.repeat(np.arange(len(base_idx)),weight[base_idx])]
            bmetric=gaussian_metrics({q:gauss_fit(ROOT,b[:,j],'boot_ref7_'+q) for j,q in enumerate(QUANTITIES)})
            if metric is not None and bmetric is not None: delta_fits.append(differences(metric,bmetric))
    result={k:dict(error_summary(rows),requested=n) for k,rows,n in
            [('moments',pooled_moments,n_moments),('gaussian',pooled_fits,n_fits)]}
    result['paired_difference_to_ref7']={
        'moments':dict(successful=len(delta_moments),errors=errors(delta_moments)),
        'gaussian':dict(successful=len(delta_fits),errors=errors(delta_fits))}
    result['per_ref_paired_moment_differences']={r:dict(successful=len(rows),errors=errors(rows)) for r,rows in per_ref.items()}
    return result

def load_members(flat,run,label,reports):
    members={}; excluded={}
    for ref in REF_BARS:
        path=flat/str(run)/f'ref{ref:02d}'/(label+'.npz')
        if ref not in reports: excluded[ref]='No completed closure; inspect production failure/log'; continue
        groups,skipped=reports[ref]
        group=groups.get(label)
        if group is None:
            skip=skipped.get(label,{})
            excluded[ref]=skip.get('reason','No selected data group'); continue
        if not path.exists(): raise RuntimeError(f'Missing extracted events: {path}')
        with np.load(path) as data: arrays={k:data[k] for k in data.files}
        if len(arrays['keys'])!=group['counts']['accepted']: raise RuntimeError('Incomplete extracted events')
        members[ref]=(group,arrays)
    return members,excluded

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('production',type=Path); p.add_argument('--flat',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True); p.add_argument('--runs',type=int,nargs='+',default=[4194,4257,4281,4291])
    p.add_argument('--bootstrap',type=int,default=200); p.add_argument('--fit-bootstrap',type=int,default=100)
    p.add_argument('--thresholds',type=int,nargs='+',default=[3,5,7,10,15,20]); p.add_argument('--bars',type=int,nargs='+',default=list(range(16)))
    a=p.parse_args(); import ROOT; ROOT.gROOT.SetBatch(True)
    a.output.mkdir(parents=True,exist_ok=True)
    for run in a.runs:
        run_lock=(a.output/f'.run-{run}.lock').open('a')
        fcntl.flock(run_lock,fcntl.LOCK_EX)
        reports={}
        for ref in REF_BARS:
            folder=a.production/str(run)/f'ref{ref:02d}'
            if not (folder/'complete.json').exists() and not (folder/'failure.json').exists():
                raise RuntimeError(f'Production still pending: {folder}')
            if (folder/'closure/summary.json').is_file():
                data=json.loads((folder/'closure/summary.json').read_text())
                reports[ref]=tuple({g['label']:g for g in data[k]} for k in ('groups','skipped'))
        for bar in a.bars:
            for th in a.thresholds:
                label=f'bar{bar:02d}_Vov3.00_th{th:02d}_energyBin01'
                folder=a.output/str(run)/label; folder.mkdir(parents=True,exist_ok=True)
                members,excluded=load_members(a.flat,run,label,reports)
                identity={str(r):dict(reference_uuid=str(d['reference_uuid']),closure_uuid=str(d['closure_uuid'])) for r,(g,d) in members.items()}
                fingerprint=dict(inputs=identity,bootstrap=a.bootstrap,fit_bootstrap=a.fit_bootstrap,
                                 implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
                done=folder/'summary.json'
                if done.exists():
                    if json.loads(done.read_text())['fingerprint']!=fingerprint: raise RuntimeError('Aggregation inputs changed; use a new output')
                    continue
                if 7 not in members: raise RuntimeError(f'No usable REF7 group: {run} {label}')
                arrays=[d for g,d in members.values()]
                keys=np.concatenate([d['keys'] for d in arrays]); ids,counts=cluster_ids(keys)
                refs=np.concatenate([np.full(len(d['keys']),r,dtype=np.int16) for r,(g,d) in members.items()])
                refdiff=np.concatenate([d['ref_diff'] for d in arrays]); energy=np.concatenate([d['ref_energy'] for d in arrays])
                # Cross-reference matches must agree in DUT energy and phase, not just timestamps.
                for feature in ('energyL','energyR','phaseL','phaseR'):
                    v=np.concatenate([d[feature] for d in arrays]); lo=np.full(len(counts),np.inf); hi=-lo
                    np.minimum.at(lo,ids,v); np.maximum.at(hi,ids,v)
                    if np.any(hi!=lo): raise RuntimeError(f'Ambiguous raw event key ({feature}) in {run} {label}')
                single=choose_single_reference(ids,refs,energy)
                report=dict(run=run,label=label,bar=bar,threshold=th,vov=3.,fingerprint=fingerprint,
                            reference_bars=list(members),excluded=excluded,pairs=len(keys),unique_events=len(counts),
                            duplicated_pairs=int(len(keys)-len(counts)),shared_events=int(np.sum(counts>1)),
                            max_references_per_event=int(counts.max()),counts_by_ref={r:len(d['keys']) for r,(g,d) in members.items()},
                            stages={},uncertainty_scope='DUT-event cluster bootstrap with saved calibration functions held fixed; no training uncertainty')
                ref7_idx=np.flatnonzero(refs==7)
                baseline_position=np.full(len(counts),-1,dtype=np.int64)
                baseline_position[ids[ref7_idx]]=ref7_idx
                paired_indices={}
                for ref in members:
                    if ref==7: continue
                    indices=np.flatnonzero((refs==ref)&(baseline_position[ids]>=0))
                    if len(indices)>=100:
                        paired_indices[ref]=(indices,baseline_position[ids[indices]])
                root=ROOT.TFile(str(folder/'histograms.root'),'RECREATE')
                for stage in ('tw','full'):
                    values=np.concatenate([d[stage] for d in arrays]); centered,means=center_by_ref(values,refs)
                    result={}
                    for variant,v,rd in [('union',values,refdiff),('centered',centered,refdiff),
                                         ('single_ref_per_event',values[single],refdiff[single])]:
                        root.cd(); s=summarize_stage(ROOT,v,variant+'_'+stage,label,0,0,20261124)
                        s['covariance']=covariance_diagnostics(v,rd)
                        result[variant]=s
                    result['union']['bootstrap']=cluster_bootstrap(ROOT,values,ids,refs,a.bootstrap,a.fit_bootstrap,20261124+run+bar+100*th)
                    result['mean_residual_by_ref_ps']=means
                    result['common_events_vs_ref7']={
                        ref:dict(events=len(i),**differences(moments(values[i]),moments(values[j])))
                        for ref,(i,j) in paired_indices.items()}
                    result['between_ref_variance_ps2']=dict(zip(QUANTITIES,np.var(values,axis=0,ddof=1)-np.var(centered,axis=0,ddof=1)))
                    result['difference_to_ref7']={est:differences(result['union'][est],members[7][0]['stages'][stage][est])
                        for est in ('moments','gaussian') if result['union'][est] is not None and members[7][0]['stages'][stage][est] is not None}
                    report['stages'][stage]=result
                root.Close(); dump_json(folder/'summary.partial.json',report); (folder/'summary.partial.json').replace(done)
                print(run,label,'pairs',len(keys),'unique',len(counts),'shared',report['shared_events'],flush=True)
        run_lock.close()

if __name__=='__main__': main()
