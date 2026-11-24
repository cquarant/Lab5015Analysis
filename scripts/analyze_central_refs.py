#!/usr/bin/env python3
"""Run the September TB analysis separately for REF 4..12, then pool the coincidences."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import tarfile
import os

REPO = Path(__file__).resolve().parents[1]
REF_BARS = tuple(range(4, 13))
POOLED_REFS = tuple(b for b in REF_BARS if b != 6)


def option(config, section, key, value=None):
    block = re.search(rf"(?ms)^<{section}>\s*$.*?^</{section}>", config)
    if block is None:
        raise ValueError(f"Missing section {section}")
    pattern = rf"(?m)^([ \t]*{key}[ \t]+)([^#\r\n]+)"
    matches = list(re.finditer(pattern, block.group()))
    if len(matches) != 1:
        raise ValueError(f"Expected one {section}.{key}")
    match = matches[0]
    if value is None:
        return match[2].strip()
    start, end = match.span(2)
    return config[:block.start()+start] + str(value) + config[block.start()+end:]


def run_logged(command, folder, log):
    with Path(log).open("w") as stream:
        subprocess.run(command, cwd=folder, stdout=stream, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="Single-run withRef TB configuration")
    parser.add_argument("output", type=Path, help="Separate output directory for this run")
    parser.add_argument("--work-dir", type=Path, help="Temporary files and unpacked diagnostic plots for this run")
    parser.add_argument("--bin-dir", type=Path, default=REPO / "bin")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--fit-bootstrap", type=int, default=100)
    args = parser.parse_args()
    if args.jobs < 1 or min(args.bootstrap, args.fit_bootstrap) < 0:
        parser.error("Invalid job or bootstrap count")
    cfg = args.config.read_text()
    if (option(cfg, "Coincidence", "status"), option(cfg, "Input", "saveRefInfo"),
        option(cfg, "Input", "sourceName")) != ("yes", "1", "TB"):
        parser.error("Expected TB with REF coincidence and saved REF channels")
    if not option(cfg, "Input", "runs").isdigit():
        parser.error("Use one run per output directory")
    mapping = list(map(int, option(cfg, "Channels", "channelMapping").split()))
    asic = int(option(cfg, "Channels", "REFasic"))
    if len(mapping) != 32 or set(mapping) != set(range(32)) or not 0 <= asic < 8:
        parser.error("Invalid reference channel mapping")
    for section, key in [('Input', 'inputDir'), ('Input', 'discCalibration'), ('Input', 'interCalibration'),
                         ('Input', 'interCalibration_ref'), ('Cuts', 'minEnergiesFileName')]:
        path = option(cfg, section, key)
        if path != "0":
            cfg = option(cfg, section, key, Path(path).resolve())
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    scratch = args.work_dir.resolve() if args.work_dir else out
    scratch.mkdir(parents=True, exist_ok=True)
    if args.work_dir:
        for key, subdir in [('TMPDIR', 'tmp'), ('MPLCONFIGDIR', 'mpl')]:
            directory = scratch / subdir
            directory.mkdir(exist_ok=True)
            os.environ[key] = str(directory)
    # Hold the directory lock through the run, including while workers write ROOT files.
    lock = (scratch / ".production.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        parser.error("This production is already running")
    manifest = json.dumps(dict(config=cfg, bootstrap=args.bootstrap, fit_bootstrap=args.fit_bootstrap,
                               ref_bars=REF_BARS, pooled_refs=POOLED_REFS), indent=2)
    provenance = out / 'production.json'
    if provenance.exists() and provenance.read_text() != manifest:
        parser.error("Output already contains a different production")
    provenance.write_text(manifest)
    boot = ["--bootstrap", str(args.bootstrap), "--fit-bootstrap", str(args.fit_bootstrap)]
    def analyze(ref):
        folder = out / f"ref{ref:02d}"
        folder.mkdir(exist_ok=True)
        workdir = scratch / f"ref{ref:02d}"
        workdir.mkdir(exist_ok=True)
        marker = folder / 'complete'
        if marker.exists():
            return
        local = cfg
        for key, channel in zip(("chL", "chR"), mapping[2*ref:2*ref+2]):
            local = option(local, "Coincidence", key, 32*asic + channel)
        for section, key, filename in [('Input', 'step1FileName', 'step1.root'),
                ('Output', 'outFileNameStep1', 'step1.root'), ('Output', 'outFileNameStep2', 'step2.root'),
                ('Output', 'refFileName', 'refCharacterization.root'), ('Output', 'plotDir', 'diagnostics')]:
            base = workdir if filename in ("step1.root", "diagnostics") else folder
            local = option(local, section, key, base / filename)
        config = folder / 'withref.cfg'
        config.write_text(local)
        for exe in ('moduleCharacterization_step1', 'refCharacterization', 'moduleCharacterization_withRef_step2'):
            done = workdir / (exe + '.done')
            if not done.exists():
                # refCharacterization recreates step2.root. Restore its clean output before
                # retrying a failed Step 2, which otherwise appends to the partial file.
                if exe == 'moduleCharacterization_withRef_step2' and (folder / (exe + '.log')).exists():
                    run_logged([args.bin_dir.resolve() / 'refCharacterization.exe', config],
                               workdir, folder/'refCharacterization.log')
                run_logged([args.bin_dir.resolve() / (exe + '.exe'), config], workdir, folder/(exe + '.log'))
                done.touch()
        if workdir != folder and not (folder / "diagnostics.tar.gz").exists():
            archive = folder / "diagnostics.tar.gz.partial"
            with tarfile.open(archive, "w:gz", compresslevel=1) as stream:
                stream.add(workdir / "diagnostics", arcname="diagnostics")
            archive.replace(folder / "diagnostics.tar.gz")
        # Only the final result payload is staged on the destination filesystem for atomic rename.
        with tempfile.TemporaryDirectory(dir=folder) as tmp:
            work = Path(tmp)
            if not (folder / 'closure').exists():
                run_logged([sys.executable, REPO/'main/moduleCharacterization_timingClosure.py',
                    '--reference', folder/'refCharacterization.root', '--corrections', folder/'step2.root',
                    '--output-dir', work/'closure', *boot], workdir, folder/'closure.log')
                (work/'closure').rename(folder/'closure')
            if not (folder / 'comparison').exists():
                run_logged([sys.executable, REPO/'macros/compare_timingClosure.py', '--new', folder/'closure/summary.json',
                    '--previous', folder/'step2.root', '--output-dir', work/'comparison'], workdir, folder/'comparison.log')
                (work/'comparison').rename(folder/'comparison')
        marker.touch()
        print(f"REF {ref} complete", flush=True)

    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(analyze, REF_BARS))
    if not (out/'pooled'/'complete').exists():
        run_logged([sys.executable, REPO/'main/pool_central_refs.py', out, *boot], scratch, out/'pooled.log')


if __name__ == "__main__":
    main()
