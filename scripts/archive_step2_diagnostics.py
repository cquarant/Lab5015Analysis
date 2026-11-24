#!/usr/bin/env python3
"""Archive the diagnostic PNGs for one validated step2.root file."""
import argparse
import json
import re
import tarfile
import tempfile
from pathlib import Path

import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--previous", type=Path, default=Path(
        "/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_central_refs_TOFHIR"))
    args = parser.parse_args()
    root_path = args.sample / "step2.root"
    ROOT.gROOT.SetBatch(True)
    style = Path(__file__).resolve().parents[1] / "src/SetTDRStyle.cc"
    ROOT.gInterpreter.AddIncludePath(str(style.parent.parent))
    if ROOT.gROOT.LoadMacro(str(style)) != 0:
        raise RuntimeError(f"Cannot load {style}")
    ROOT.setTDRStyle()
    ROOT.gStyle.SetOptStat(1111)
    ROOT.gErrorIgnoreLevel = ROOT.kError

    source = ROOT.TFile.Open(str(root_path))
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open {root_path}")
    uuid = source.GetUUID().AsString()
    summary = json.loads((args.sample / "closure/summary.json").read_text())
    if uuid != summary["input_uuids"]["corrections"]:
        raise RuntimeError(f"Unexpected step2 UUID for {root_path}")
    keys = {k.GetName() for k in source.GetListOfKeys()
            if k.GetClassName() in {"TH1F", "TH2F", "TProfile"}}
    if not keys:
        raise RuntimeError(f"No diagnostic objects in {root_path}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(".partial")
    preserved = 0
    # Detect restored ROOT files whose neighboring diagnostics may be stale.
    existing = args.sample / "diagnostics"
    root_mtime = root_path.stat().st_mtime
    reuse = existing.is_dir() and abs(root_mtime - existing.stat().st_mtime) <= 5
    previous = args.previous / args.sample.parent.name / args.sample.name
    previous_tar = previous / "diagnostics.tar.gz"
    original_archive = None
    if not reuse and previous_tar.is_file():
        previous_root = ROOT.TFile.Open(str(previous / "step2.root"), "READ")
        if previous_root and not previous_root.IsZombie():
            if previous_root.GetUUID().AsString() == uuid:
                original_archive = tarfile.open(previous_tar, "r:gz")
            previous_root.Close()
    canvas = ROOT.TCanvas("step2_diagnostics", "", 600, 500)
    with tempfile.TemporaryDirectory() as scratch, tarfile.open(partial, "w:gz", compresslevel=1) as archive:
        originals = {}
        if original_archive:
            with original_archive:
                for member in original_archive:
                    name = Path(member.name).name
                    if member.isfile() and name.startswith("c_") and name[2:-4] in keys:
                        path = Path(scratch) / name
                        path.write_bytes(original_archive.extractfile(member).read())
                        originals[name] = path
        for name in sorted(keys):
            match = re.match(r"^[^_]+_(.+?)_bar", name)
            quantity = match[1].split("vs_", 1)[0].rstrip("_") if match else "UNKNOWN"
            relative = Path("diagnostics") / quantity / f"c_{name}.png"
            saved = args.sample / relative
            if reuse and saved.is_file() and saved.stat().st_mtime <= root_mtime:
                png = saved
                preserved += 1
            elif relative.name in originals:
                png = originals[relative.name]
                preserved += 1
            else:
                obj = source.Get(name)
                canvas.cd()
                obj.Draw("COLZ" if obj.InheritsFrom("TH2") else "")
                if obj.InheritsFrom("TProfile") and not obj.GetListOfFunctions().GetSize():
                    fit = source.Get(f"f_{name}")
                    if fit:
                        fit.Draw("same")
                png = Path(scratch) / relative.name
                canvas.Print(str(png))
                canvas.Clear()
            with png.open("rb") as image:
                if image.read(8) != b"\x89PNG\r\n\x1a\n":
                    raise RuntimeError(f"Invalid PNG {png}")
            archive.add(png, arcname=str(relative))
            if png != saved:
                png.unlink()
    receipt = dict(step2_uuid=uuid, plots=len(keys), original_pngs=preserved,
                   regenerated_pngs=len(keys) - preserved)
    with tarfile.open(partial, "r:gz") as archive:
        names = [member.name for member in archive]
    if len(names) != len(keys) or len(set(names)) != len(keys):
        raise RuntimeError(f"Invalid archive inventory {partial}")
    partial.replace(args.output)
    receipt_path = args.output.with_name("step2_diagnostics.json")
    temporary_receipt = receipt_path.with_suffix(".json.partial")
    temporary_receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    temporary_receipt.replace(receipt_path)
    print(args.sample, receipt, flush=True)


if __name__ == "__main__":
    main()
