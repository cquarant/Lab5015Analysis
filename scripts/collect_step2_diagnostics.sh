#!/usr/bin/env bash
set -euo pipefail
production=${1:?Provide the validated production_final directory}
output=${2:?Provide a staging directory}
jobs=${3:-4}
renderer=$(dirname "$(realpath "$0")")/archive_step2_diagnostics.py
mkdir -p "$output"
for run in 4194 4257 4281 4291; do
    for ref in 04 05 07 08 09 10 11 12; do
        printf '%s\n' "$run/ref$ref"
    done
done | xargs -P "$jobs" -I '{}' bash -c '
    set -euo pipefail
    sample=$1
    target=$3/$sample
    if [[ -s "$target/step2_diagnostics.tar.gz" && -s "$target/step2_diagnostics.json" ]]; then
        exit 0
    fi
    mkdir -p "$target"
    python3 "$4" "$2/$sample" "$target/step2_diagnostics.tar.gz" >"$target/render.log" 2>&1
' _ '{}' "$production" "$output" "$renderer"
