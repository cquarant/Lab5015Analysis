#!/usr/bin/env bash
# Run from the analysis repository on lxplus. Outputs and receipts are persistent.
set -euo pipefail
base=/eos/home-c/cquarant/TB_Sept2025_Studies4DN
prior="$base/timing_closure_multiref_TOFHIR"
central=/afs/cern.ch/work/c/cquarant/btl_central_refs_20261124
out=${1:-/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124}
scratch=$(mktemp -d /tmp/btl-multiref-selected.XXXXXXXX)
trap 'rm -rf -- "$scratch"' EXIT
mkdir -p "$out/bin" "$out/logs" "$out/published_validation" "$scratch/ccache" "$scratch/ccache-tmp" "$scratch/mpl"
export TMPDIR="$scratch" CCACHE_DIR="$scratch/ccache" CCACHE_TEMPDIR="$scratch/ccache-tmp"
export MPLCONFIGDIR="$scratch/mpl" OPENBLAS_NUM_THREADS=1
exec 9>"$out/.resume.lock"
flock -n 9
cp "$prior/published_validation/validation.json" "$out/published_validation/validation.json"

if [[ ! -x "$out/bin/moduleCharacterization_withRef_step2" ]]; then
    rc=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.36.00/x86_64-almalinux9.5-gcc115-opt/bin/root-config
    g++ -O2 -I. $("$rc" --cflags) main/moduleCharacterization_withRef_step2.cpp \
        src/referenceCharacterization_helper.cc -L./lib -lTBStudies -L./DynamicTTree/lib -lDynamicTTree \
        -L./CfgManager/lib -lCfgManager $("$rc" --glibs) -lMinuit -lTreePlayer -lSpectrum \
        -o "$out/bin/moduleCharacterization_withRef_step2.partial" > "$out/logs/build.log" 2>&1
    mv "$out/bin/moduleCharacterization_withRef_step2.partial" "$out/bin/moduleCharacterization_withRef_step2"
fi
cp "$central/bin/moduleCharacterization_step1.exe" "$central/bin/refCharacterization.exe" "$out/bin/"
sed "s|$central/bin/moduleCharacterization_withRef_step2|$out/bin/moduleCharacterization_withRef_step2|" \
    "$central/bin/moduleCharacterization_withRef_step2.exe" > "$out/bin/moduleCharacterization_withRef_step2.exe"
chmod +x "$out/bin/moduleCharacterization_withRef_step2.exe"

extra=()
if [[ -n ${MULTIREF_RECOVER_LOCAL:-} && -d $MULTIREF_RECOVER_LOCAL ]]; then
    extra=(--reuse "$MULTIREF_RECOVER_LOCAL")
fi
for run in 4194 4257 4281 4291; do
    python3 scripts/analyze_multiref.py --baseline "$base/timing_closure_withref_TOFHIR" \
        --validation "$out/published_validation/validation.json" --output "$out/production_final" \
        --bin-dir "$out/bin" --scratch "$scratch" --runs "$run" --jobs 3 --ref-bars 7 4 5 8 9 10 11 12 \
        "${extra[@]}" --reuse "$prior/production_final" --reuse "$prior/production" \
        --reuse "$base/timing_closure_central_refs_TOFHIR" --reuse-work "$central/work" \
        >> "$out/logs/production-$run.log" 2>&1
    python3 scripts/extract_multiref_events.py "$out/production_final" --output "$out/flat_final" --runs "$run" \
        >> "$out/logs/extraction-$run.log" 2>&1
done
for run in 4194 4257 4281 4291; do
    python3 main/aggregate_multiref.py "$out/production_final" --flat "$out/flat_final" \
        --output "$out/aggregate" --runs "$run" >> "$out/logs/aggregation-$run.log" 2>&1
done
python3 scripts/report_multiref.py --work "$out" --baseline "$base/timing_closure_withref_TOFHIR" \
    --output "$out/report" --style "$base/../www/cms_btl_timing_closure/assets/style.css" >> "$out/logs/report.log" 2>&1
python3 scripts/verify_multiref_selected.py "$out" --safe-baseline "$prior/baseline_validation" \
    >> "$out/logs/verification.log" 2>&1
