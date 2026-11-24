#!/usr/bin/env python3
"""Event-cluster statistics for separately calibrated DUT–REF coincidences.

The primary estimand is the variance of the union of accepted pairs. A DUT event
may contribute multiple pairs, but bootstrap resampling is by raw DUT timestamp
pair. Centered and one-reference-per-event samples are explicit sensitivities.
"""
import numpy as np
from moduleCharacterization_timingClosure import QUANTITIES, derived_variances, moments

METHODS=('direct_halfdiff','avg_mixed','matched_halfdiff','catr')
REF_BARS=(4,5,7,8,9,10,11,12)

def cluster_ids(keys):
    keys=np.ascontiguousarray(keys,dtype=np.int64)
    if keys.ndim!=2 or keys.shape[1]!=2: raise ValueError('Expected two raw integer DUT timestamps')
    _, inverse, counts=np.unique(keys,axis=0,return_inverse=True,return_counts=True)
    return inverse,counts

def covariance_diagnostics(values, ref_diff):
    c=np.cov(values,rowvar=False,ddof=1)
    varref=np.var(ref_diff,ddof=1)/4
    m=moments(values)
    return dict(residual_covariance_ps2=float(c[0,1]),
                from_average_minus_matched_ps2=float(c[2,2]-c[3,3]/4),
                ref_halfdiff_variance_ps2=float(varref),
                intrinsic_covariance_proxy_ps2=float(c[0,1]-varref),
                mixed_minus_direct_variance_ps2=float(m['avg_mixed']-m['direct_halfdiff']),
                mixed_minus_matched_variance_ps2=float(m['avg_mixed']-m['matched_halfdiff']))

def center_by_ref(values, refs):
    centered=np.array(values,copy=True)
    target=np.mean(values,axis=0)
    means={}
    for ref in np.unique(refs):
        mask=refs==ref
        means[int(ref)]=values[mask].mean(axis=0)
        centered[mask]-=means[int(ref)]-target
    return centered,means

def choose_single_reference(ids, refs, ref_energy):
    # Highest calibrated REF mean energy; deterministic lower bar number breaks ties.
    order=np.lexsort((refs,-ref_energy,ids))
    return order[np.r_[True,ids[order][1:]!=ids[order][:-1]]]

def weighted_moments(values,weights):
    n=weights.sum()
    if n<2: return None
    # Center before squaring to avoid cancellation from ~1500 ps residual offsets.
    v=values-values.mean(axis=0)
    mean=np.sum(weights[:,None]*v,axis=0)/n
    variance=(np.sum(weights[:,None]*v*v,axis=0)-n*mean*mean)/(n-1)
    return derived_variances(dict(zip(QUANTITIES,variance)))

def error_summary(rows):
    out=dict(successful=len(rows),variance_errors_ps2={},resolution_errors_ps={})
    if len(rows)<2: return out
    for k in rows[0]:
        v=np.array([x[k] for x in rows],dtype=float)
        if np.isfinite(v).all():
            out['variance_errors_ps2'][k]=float(np.std(v,ddof=1))
            if k in METHODS and np.all(v>=0): out['resolution_errors_ps'][k]=float(np.std(np.sqrt(v),ddof=1))
    return out
