#!/usr/bin/env python3
"""Replay saved withRef corrections and compare timing estimators on identical events.

No calibration is fitted here: this deliberately isolates reconstruction/estimator
differences from calibration changes. See ../TIMING_CLOSURE.md.
"""
import argparse
from array import array
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np

QUANTITIES = ("left", "right", "average", "matched_diff", "direct_diff")
PREFIXES = {
    "left": "h1_deltaT_tL_tAveRefCor_eCor_phaseCor_",
    "right": "h1_deltaT_tR_tAveRefCor_eCor_phaseCor_",
    "average": "h1_deltaT_tAve_tAveRefCor_eCor_phaseCor_",
    "direct_diff": "h1_deltaT_tL_tR_eRatioCor_phaseMeanCor_",
}


def derived_variances(v):
    """Keep signed radicands: never turn an unphysical extraction into zero."""
    x, y, a, d, q = (v[k] for k in QUANTITIES)
    ref_matched = (x + y - d) / 2
    ref_mixed = (x + y - q) / 2
    return {
        "matched_halfdiff": d / 4,
        "direct_halfdiff": q / 4,
        "catr": (x + y - 2*a) / 2,
        "ref_matched": ref_matched,
        "ref_mixed": ref_mixed,
        "avg_matched": a - ref_matched,
        "avg_mixed": a - ref_mixed,
        "closure_ps2": (x + y - 2*a) / 2 - d / 4,
        "correction_gain_ps2": (q - d) / 4,
    }


def moments(values):
    v = dict(zip(QUANTITIES, np.var(values, axis=0, ddof=1)))
    return derived_variances(v)


def finite_json(value):
    if isinstance(value, dict):
        return {str(k): finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite_json(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def dump_json(path, data):
    Path(path).write_text(json.dumps(finite_json(data), indent=2, allow_nan=False) + "\n")


def open_root(ROOT, path):
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT input: {path}")
    return f


def require(f, name, kind):
    obj = f.Get(name)
    if not obj or not obj.InheritsFrom(kind):
        raise RuntimeError(f"Missing {kind} {name} in {f.GetName()}")
    return obj


def event_dictionary(ROOT, repository, cache):
    """Build only the event dictionary, using the repository's actual class schema."""
    header = (Path(repository) / "interface/AnalysisUtils.h").resolve()
    if not header.is_file():
        raise RuntimeError(f"Missing event schema: {header}")
    cache = Path(cache).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    old = Path.cwd()
    try:
        os.chdir(cache)
        if not ROOT.gInterpreter.GenerateDictionary("ModuleEventWithRefClass", str(header)):
            raise RuntimeError("Could not compile the ROOT event dictionary; check compiler/ROOT setup")
    finally:
        os.chdir(old)
    return hashlib.sha256(header.read_bytes()).hexdigest()


def gauss_fit(ROOT, values, name, save=False):
    """Match FitAndSaveHisto's 400-bin, peak +/- RMS Gaussian-core prescription."""
    is_diff = "diff" in name
    lo, hi = (-1000., 1000.) if is_diff else (500., 2500.)
    h = ROOT.TH1D(name, name + ";time residual [ps];Events", 400, lo, hi)
    h.SetDirectory(0)
    vals = np.asarray(values, dtype=np.float64)
    vals = np.ascontiguousarray(vals[np.isfinite(vals)])
    h.FillN(len(vals), vals, np.ones(len(vals), dtype=np.float64))
    result = {"sigma": None, "sigma_error": None, "status": -1,
              "entries": len(vals), "underflow": h.GetBinContent(0),
              "overflow": h.GetBinContent(401)}
    if len(vals) >= 10 and h.GetRMS() > 0:
        mean, rms = h.GetMean(), h.GetRMS()
        # ROOT otherwise lets the underflow/overflow bins seed GetMaximumBin.
        h.GetXaxis().SetRangeUser(max(lo, mean - 5*rms), min(hi, mean + 5*rms))
        peak, rms = h.GetBinCenter(h.GetMaximumBin()), h.GetRMS()
        xmin, xmax = max(lo, peak-rms), min(hi, peak+rms)
        if xmin < xmax:
            f = ROOT.TF1("fit_" + name, "gaus", xmin, xmax)
            f.SetParameters(h.GetMaximum(), peak, 20.)
            fit = h.Fit(f, "QRSN")
            status = int(fit)
            result.update(status=status, fit_min=xmin, fit_max=xmax,
                          covariance_status=int(fit.CovMatrixStatus()))
            if status == 0 and math.isfinite(f.GetParameter(2)) and f.GetParameter(2) > 0:
                result.update(sigma=f.GetParameter(2), sigma_error=f.GetParError(2))
            if save:
                f.Write()
    if save:
        h.Write()
    return result


def gaussian_metrics(fits):
    if any(fits[k]["sigma"] is None for k in QUANTITIES):
        return None
    return derived_variances({k: fits[k]["sigma"]**2 for k in QUANTITIES})


def bootstrap(ROOT, values, n_moments, n_fits, seed):
    """Paired event bootstrap, with calibration constants held fixed."""
    rng = np.random.default_rng(seed)
    results = {"moments": [], "gaussian": []}
    for i in range(max(n_moments, n_fits)):
        sample = values[rng.integers(0, len(values), len(values))]
        if i < n_moments:
            results["moments"].append(moments(sample))
        if i < n_fits:
            fits = {k: gauss_fit(ROOT, sample[:, j], "boot_" + k)
                    for j, k in enumerate(QUANTITIES)}
            metric = gaussian_metrics(fits)
            if metric is not None:
                results["gaussian"].append(metric)
    out = {}
    for estimator, rows in results.items():
        summary = {"requested": n_moments if estimator == "moments" else n_fits,
                   "successful": len(rows), "variance_errors_ps2": {}, "resolution_errors_ps": {}}
        if len(rows) >= 2:
            for key in rows[0]:
                v = np.array([r[key] for r in rows])
                summary["variance_errors_ps2"][key] = np.std(v, ddof=1)
                # Do not silently drop negative radicands in uncertainty estimation.
                if not key.endswith("ps2") and np.all(v >= 0):
                    summary["resolution_errors_ps"][key] = np.std(np.sqrt(v), ddof=1)
        out[estimator] = summary
    return out


def summarize_stage(ROOT, values, stage, label, n_moments, n_fits, seed):
    fits = {q: gauss_fit(ROOT, values[:, j], f"{stage}_{q}_{label}", save=True)
            for j, q in enumerate(QUANTITIES)}
    exact = moments(values)
    tolerance = 1e-10*max(1., exact["matched_halfdiff"])
    if abs(exact["closure_ps2"]) > tolerance:
        raise RuntimeError(f"Variance identity failed in {label}: {exact['closure_ps2']}")
    return {"moments": exact, "gaussian": gaussian_metrics(fits), "fits": fits,
        "bootstrap": bootstrap(ROOT, values, n_moments, n_fits, seed),
        "residual_covariance_ps2": np.cov(values, rowvar=False, ddof=1).tolist()}


def read_ranges(f):
    tree = require(f, "ranges", "TTree")
    out = {}
    for row in tree:
        # ROOT exposes a char[] as either a Python string or a LowLevelView,
        # depending on version. Read it through TLeafC in both cases.
        lr = str(tree.GetLeaf("LR").GetValueString())
        key = (int(row.bar), int(float(row.Vov)*100), int(row.vth))
        bounds = (float(row.rmin), float(row.rmax))
        if not all(math.isfinite(x) for x in bounds) or bounds[0] >= bounds[1]:
            raise RuntimeError(f"Invalid energy range for {key} {lr}: {bounds}")
        out.setdefault(key, {})[lr] = bounds
    return out


def load_corrections(f, label):
    spec = {
        "twL": ("f_p1_deltaT_tL_tAveRefCor_vs_eL_", "TF1"),
        "twR": ("f_p1_deltaT_tR_tAveRefCor_vs_eR_", "TF1"),
        "twD": ("f_p1_deltaT_tL_tR_vs_eRatio_", "TF1"),
        "phaseL": ("p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL_", "TProfile"),
        "phaseR": ("p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR_", "TProfile"),
        "phaseD": ("p1_deltaT_tL_tR_eRatioCor_vs_phaseMean_", "TProfile"),
        "meanL": ("h1_phaseL_", "TH1"), "meanR": ("h1_phaseR_", "TH1"),
        "meanD": ("h1_phaseAve_", "TH1"),
    }
    c = {k: require(f, prefix + label, kind) for k, (prefix, kind) in spec.items()}
    for side, stem, suffix in [("L", "h1_eL_", "landau"), ("R", "h1_eR_", "landau"),
                               ("D", "h1_eRatio_", "gaus")]:
        fit = require(f, "f_" + stem + label + "_" + suffix, "TF1")
        c["anchor" + side] = c["tw" + side].Eval(fit.GetParameter(1))
    return c


def corrected_row(e, c):
    # Subtract integer timestamps before conversion: preserve ps precision at large coarse times.
    tl, tr = int(e.timeL), int(e.timeR)
    fl, fr = int(e.timeL_ref_cor), int(e.timeR_ref_cor)
    x0, y0 = .5*((tl-fl)+(tl-fr)), .5*((tr-fl)+(tr-fr))
    el, er = float(e.energyL), float(e.energyR)
    tw = [c["twL"].Eval(el)-c["anchorL"],
          c["twR"].Eval(er)-c["anchorR"],
          c["twD"].Eval(el/er)-c["anchorD"]]
    phase = [float(e.t1fineL), float(e.t1fineR), .5*(int(e.t1fineL)+int(e.t1fineR))]
    pc, invalid = [], False
    for side, value in zip(("L", "R", "D"), phase):
        p = c["phase" + side]
        b, a = p.FindBin(value), p.FindBin(c["mean" + side].GetMean())
        invalid |= any(i < 1 or i > p.GetNbinsX() or p.GetBinEntries(i) == 0 for i in (a, b))
        # Exact legacy behavior retained, including empty-bin zero contents; flag separately.
        pc.append(p.GetBinContent(b)-p.GetBinContent(a))
    def pack(x, y, d):
        return [x, y, .5*(x+y), x-y, d]
    twrow = pack(x0-tw[0], y0-tw[1], (tl-tr)-tw[2])
    fullrow = pack(x0-tw[0]-pc[0], y0-tw[1]-pc[1], (tl-tr)-tw[2]-pc[2])
    return twrow, fullrow, invalid


def passes(e, ranges):
    key = (int(e.barID), int(float(e.Vov)*100), int(e.vth))
    if key not in ranges or not {"L", "R", "L-R"}.issubset(ranges[key]):
        raise RuntimeError(f"Missing selection ranges for {key}")
    tl, tr, fl, fr = map(int, (e.timeL, e.timeR, e.timeL_ref, e.timeR_ref))
    residuals = [.5*((fl-tl)+(fl-tr)), .5*((fr-tl)+(fr-tr)),
                 .5*((tl-fl)+(tl-fr)), .5*((tr-fl)+(tr-fr))]
    if any(abs(t) > 5000 for t in residuals):
        return "timing_cut"
    if any(not math.isfinite(float(v)) or v < 0 or v > 25
           for v in (e.totL, e.totR, e.totL_ref, e.totR_ref)):
        return "tot_cut"
    if not all(math.isfinite(float(v)) for v in (e.energyL, e.energyR)) or e.energyR == 0:
        return "invalid_energy"
    # Float32 addition matches the C++ energy-bin calculation (FindBin excludes upper edge).
    energy = .5*float(np.float32(e.energyL)+np.float32(e.energyR))
    lo, hi = ranges[key]["L-R"]
    if not lo <= energy < hi:
        return "energy_cut"
    return None


def write_events(ROOT, label, entries, tw, full, flags):
    tree = ROOT.TTree("events_" + label, "Same-event corrected residuals; times in ps")
    buffers = {}
    for name in [s + "_" + q for s in ("tw", "full") for q in QUANTITIES]:
        buffers[name] = array("d", [0.])
        tree.Branch(name, buffers[name], name + "/D")
    entrybuf, flagbuf = array("q", [0]), array("i", [0])
    tree.Branch("source_entry", entrybuf, "source_entry/L")
    tree.Branch("invalid_phase_bin", flagbuf, "invalid_phase_bin/I")
    for i, entry in enumerate(entries):
        entrybuf[0], flagbuf[0] = entry, int(flags[i])
        for stage, rows in (("tw", tw), ("full", full)):
            for j, q in enumerate(QUANTITIES):
                buffers[stage + "_" + q][0] = rows[i, j]
        tree.Fill()
    tree.Write()
    tree.SetDirectory(0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, help="Completed refCharacterization ROOT file (data_* + ranges)")
    parser.add_argument("--corrections", required=True, help="Matching withRef step2 ROOT output containing saved fits/profiles")
    parser.add_argument("--repository", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--output-dir", required=True, help="New directory; existing directories are refused")
    parser.add_argument("--bars", type=int, nargs="+")
    parser.add_argument("--vov", type=float)
    parser.add_argument("--threshold", type=int)
    parser.add_argument("--bootstrap", type=int, default=200, help="Paired moment bootstrap replicates")
    parser.add_argument("--fit-bootstrap", type=int, default=0, help="Optional paired Gaussian-refit replicates")
    parser.add_argument("--seed", type=int, default=20261116)
    parser.add_argument("--min-events", type=int, default=100)
    parser.add_argument("--no-event-tree", action="store_true")
    parser.add_argument("--skip-unusable", action="store_true", help="Record groups lacking fitted corrections instead of aborting the whole REF sample")
    args = parser.parse_args()
    if args.bootstrap < 0 or args.fit_bootstrap < 0 or args.min_events < 2:
        parser.error("Bootstrap counts must be nonnegative and min-events >= 2")
    out = Path(args.output_dir).resolve()
    if out.exists():
        parser.error(f"Output directory already exists: {out}")
    import ROOT
    ROOT.gROOT.SetBatch(True)
    with tempfile.TemporaryDirectory(prefix="timing-closure-dictionary-") as cache:
        schema_hash = event_dictionary(ROOT, args.repository, cache)
        reference = open_root(ROOT, args.reference)
        corrections = open_root(ROOT, args.corrections)
        ranges = read_ranges(reference)
        names = sorted({k.GetName() for k in reference.GetListOfKeys()
                        if k.GetName().startswith("data_") and k.GetClassName() == "TTree"})
        if not names:
            raise RuntimeError("No data_* trees in reference input")
        out.mkdir(parents=True)
        result_file = ROOT.TFile(str(out / "timing_closure.root"), "CREATE")
        if result_file.IsZombie():
            raise RuntimeError("Cannot create output ROOT file")
        report = {"schema_version": 1, "arguments": vars(args), "root_version": ROOT.gROOT.GetVersion(),
                  "event_schema_sha256": schema_hash,
                  "input_uuids": {"reference": reference.GetUUID().AsString(),
                                  "corrections": corrections.GetUUID().AsString()},
                  "uncertainty_scope": "Paired event bootstrap conditional on fixed saved corrections; no calibration uncertainty",
                  "groups": [], "skipped": []}
        try:
            for name in names:
                tree = reference.Get(name)
                if not tree.GetEntries():
                    continue
                tree.GetEntry(0)
                e = tree.event
                bar, vov, th = int(e.barID), float(e.Vov), int(e.vth)
                if ((args.bars is not None and bar not in args.bars) or
                    (args.vov is not None and abs(vov-args.vov) > 1e-4) or
                    (args.threshold is not None and th != args.threshold)):
                    continue
                label = f"bar{bar:02d}_Vov{vov:.2f}_th{th:02d}_energyBin01"
                try:
                    c = load_corrections(corrections, label)
                except RuntimeError as error:
                    if not args.skip_unusable:
                        raise
                    report['skipped'].append(dict(label=label, tree=name, bar=bar, vov=vov, threshold=th,
                        counts=dict(input=int(tree.GetEntries()), accepted=0), stages={},
                        reason='unusable corrections: '+str(error)))
                    continue
                entries, tw, full, flags = [], [], [], []
                counts = Counter()
                for i in range(tree.GetEntries()):
                    tree.GetEntry(i)
                    e = tree.event
                    counts["input"] += 1
                    if (int(e.barID), float(e.Vov), int(e.vth)) != (bar, vov, th):
                        raise RuntimeError(f"Mixed bar/Vov/threshold in {name}; refusing ambiguous grouping")
                    reason = passes(e, ranges)
                    if reason:
                        counts[reason] += 1
                        continue
                    a, b, invalid = corrected_row(e, c)
                    if not np.isfinite(a+b).all():
                        counts["nonfinite_correction"] += 1
                        continue
                    entries.append(i); tw.append(a); full.append(b); flags.append(invalid)
                    counts["accepted"] += 1
                    counts["invalid_phase_bin"] += int(invalid)
                tw, full = np.asarray(tw), np.asarray(full)
                group = {"label": label, "tree": name, "bar": bar, "vov": vov, "threshold": th,
                         "counts": dict(counts), "source_entries_sha256": hashlib.sha256(np.array(entries, dtype="<i8").tobytes()).hexdigest(),
                         "stages": {}}
                result_file.cd()
                if not args.no_event_tree:
                    write_events(ROOT, label, entries, tw, full, flags)
                if len(entries) < args.min_events:
                    group["reason"] = "too few events"
                    report["skipped"].append(group)
                    continue
                for stage, values in (("tw", tw), ("full", full)):
                    group["stages"][stage] = summarize_stage(
                        ROOT, values, stage, label, args.bootstrap, args.fit_bootstrap, args.seed + bar + 100*th)
                report["groups"].append(group)
                print(label, "accepted", len(entries), "matched/direct RMS half-diff [ps]",
                      math.sqrt(group["stages"]["full"]["moments"]["matched_halfdiff"]),
                      math.sqrt(group["stages"]["full"]["moments"]["direct_halfdiff"]), flush=True)
            if not report["groups"]:
                raise RuntimeError("No groups passed the requested filters/minimum event count")
            dump_json(out / "summary.json", report)
            ROOT.TNamed("provenance", json.dumps(finite_json({k:v for k,v in report.items() if k != "groups"}))).Write()
        finally:
            result_file.Close()
            reference.Close()
            corrections.Close()
    print("Wrote", out)


if __name__ == "__main__":
    main()
