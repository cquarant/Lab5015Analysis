#!/usr/bin/env python3
"""Compare closure output with saved previous withRef step-2 results."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "main"))
from moduleCharacterization_timingClosure import PREFIXES, derived_variances, dump_json, open_root


def previous_fit(f, name):
    h = f.Get(name)
    if not h or not h.InheritsFrom("TH1"):
        raise RuntimeError(f"Missing previous histogram: {name}")
    fit = f.Get("f_" + name + "_gaus")
    if not fit:
        fit = h.GetFunction("f_" + name + "_gaus")
    if not fit:
        raise RuntimeError(f"Missing saved Gaussian fit for {name}")
    sigma = float(fit.GetParameter(2))
    if not math.isfinite(sigma) or sigma <= 0:
        raise RuntimeError(f"Invalid previous Gaussian width for {name}: {sigma}")
    return {"sigma": sigma, "sigma_error": fit.GetParError(2),
            "entries": h.GetEntries(), "histogram": name,
            "underflow": h.GetBinContent(0), "overflow": h.GetBinContent(h.GetNbinsX()+1)}


def root_or_none(v):
    return math.sqrt(v) if v is not None and math.isfinite(v) and v >= 0 else None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--new", required=True, help="summary.json from timingClosure")
    p.add_argument("--previous", required=True, help="Previous withRef step-2 ROOT file")
    p.add_argument("--output-dir", required=True, help="New directory; existing directories refused")
    args = p.parse_args()
    out = Path(args.output_dir)
    if out.exists():
        p.error(f"Output directory already exists: {out}")
    import ROOT
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ROOT.gROOT.SetBatch(True)
    report = json.loads(Path(args.new).read_text())
    f = open_root(ROOT, args.previous)
    comparisons, rows = [], []
    for group in report["groups"]:
        label = group["label"]
        previous = {q: previous_fit(f, prefix + label) for q, prefix in PREFIXES.items()}
        v = {q: x["sigma"]**2 for q, x in previous.items()}
        # matched_diff is unavailable in old output; only quantities independent of it are used.
        v["matched_diff"] = float("nan")
        old = derived_variances(v)
        row = {"label": label, "bar": group["bar"], "vov": group["vov"], "threshold": group["threshold"],
               "accepted_new": group["counts"]["accepted"],
               "previous_direct_ps": root_or_none(old["direct_halfdiff"]),
               "previous_avg_mixed_ps": root_or_none(old["avg_mixed"]),
               "previous_catr_ps": root_or_none(old["catr"])}
        replay = group["stages"]["full"]
        for estimator in ("moments", "gaussian"):
            metrics = replay[estimator]
            for key in ("matched_halfdiff", "direct_halfdiff", "catr", "avg_mixed", "avg_matched"):
                row[f"{estimator}_{key}_ps"] = root_or_none(metrics[key]) if metrics else None
        for q in PREFIXES:
            new_sigma = replay["fits"][q]["sigma"]
            row[f"replay_minus_previous_{q}_width_ps"] = new_sigma-previous[q]["sigma"] if new_sigma else None
            row[f"previous_{q}_entries"] = previous[q]["entries"]
        rows.append(row)
        comparisons.append({"label": label, "previous": previous, "values": row})
    f.Close()
    if not rows:
        raise RuntimeError("No matching groups to compare")
    out.mkdir(parents=True)
    with (out / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    dump_json(out / "comparison.json", {"arguments": vars(args), "groups": comparisons,
        "notes": ["Previous outputs lack event IDs: equal populations cannot be established from histogram counts alone.",
                  "No independent-error pulls: old/new results can share events and calibrations.",
                  "Previous fit statuses were not persisted by the original helper.",
                  "Gaussian derived errors use paired refits only when requested in timingClosure; fixed calibration."]})
    for vov, threshold in sorted({(r["vov"], r["threshold"]) for r in rows}):
        selected = sorted([r for r in rows if r["vov"] == vov and r["threshold"] == threshold], key=lambda r:r["bar"])
        groups = {g["label"]: g for g in report["groups"]}
        for estimator in ("moments", "gaussian"):
            fig, ax = plt.subplots(figsize=(10, 6))
            series = [(f"{estimator}_direct_halfdiff_ps", "Direct ratio/phase correction", "direct_halfdiff"),
                      (f"{estimator}_matched_halfdiff_ps", "Difference of corrected channels", "matched_halfdiff"),
                      (f"{estimator}_catr_ps", "CatR formula, corrected channels", "catr"),
                      (f"{estimator}_avg_mixed_ps", "Current macro: mixed corrections", "avg_mixed")]
            if estimator == "gaussian":
                series.extend([("previous_direct_ps", "Previous direct difference", None),
                               ("previous_avg_mixed_ps", "Previous mixed extraction", None)])
            for column, title, metric in series:
                valid = [r for r in selected if r[column] is not None]
                if not valid:
                    continue
                x, y = [r["bar"] for r in valid], [r[column] for r in valid]
                errors = [groups[r["label"]]["stages"]["full"]["bootstrap"][estimator]
                          ["resolution_errors_ps"].get(metric) for r in valid] if metric else []
                if errors and all(e is not None for e in errors):
                    ax.errorbar(x, y, yerr=errors, marker="o", capsize=2, label=title)
                else:
                    ax.plot(x, y, marker="o", label=title)
            ax.set(xlabel="Bar", ylabel="Extracted width [ps]",
                   title=f"Vov {vov:.2f} V, threshold {threshold}: {estimator}")
            ax.grid(alpha=.3); ax.legend(fontsize=9)
            fig.tight_layout()
            for ext in ("png", "pdf"):
                fig.savefig(out / f"timing_{estimator}_Vov{vov:.2f}_th{threshold:02d}.{ext}")
            plt.close(fig)
    print("Wrote", out)


if __name__ == "__main__":
    main()
