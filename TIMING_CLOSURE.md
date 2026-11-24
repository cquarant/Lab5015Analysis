Timing estimator closure for the September 2025 BTL analysis
==========================================================

These additive scripts leave all existing analysis scripts and outputs unchanged.
They diagnose the estimator/correction comparison discussed in DN-26-014.

The implementation deliberately **replays saved corrections** rather than fitting
new calibrations. It consumes the completed `refCharacterization` output (including
corrected REF times and the `ranges` tree) and the corresponding
`moduleCharacterization_withRef_step2` output (TF1 fits and phase profiles).
It is for the single-bar withRef analysis, not the separate two-hit energy-sharing
reconstruction. Do not pass the raw step-1 file as `--reference`.

Run from the analysis repository with Python 3, NumPy, PyROOT, and a C++ compiler
matching ROOT. The comparison additionally needs matplotlib. The existing
`source scripts/setup.sh` selects the site's ROOT environment. No Makefile change
or full rebuild is required: a temporary dictionary is compiled from the actual
`interface/AnalysisUtils.h` so the stored event class can be read.

```bash
python3 main/moduleCharacterization_timingClosure.py \
  --reference /path/to/refCharacterization_RUN.root \
  --corrections /path/to/moduleCharacterization_step2_RUN.root \
  --output-dir /path/to/new_timing_closure_RUN \
  --vov 3.0 --threshold 10 --bootstrap 200 --fit-bootstrap 100

python3 macros/compare_timingClosure.py \
  --new /path/to/new_timing_closure_RUN/summary.json \
  --previous /path/to/moduleCharacterization_step2_RUN.root \
  --output-dir /path/to/new_timing_comparison_RUN
```

For a quick replay use `--bootstrap 0 --fit-bootstrap 0`. Omit Vov/threshold
filters to process all settings, or select bars with `--bars 2 3 7`.
Each output directory must be new. Inputs are read-only. Missing correction
objects/ranges fail explicitly, rather than substituting zero corrections.
`--repository /path/to/Simona_TB_CERN_Sep25` supports running a staged copy elsewhere.

Analysis details
----------------

The event selection reproduces `PassSelection`: four raw time residuals within
5000 ps, four ToTs within [0,25], and mean DUT energy inside the stored L-R range.
The upper energy edge is excluded, matching `FindBin`. Nonfinite input energies
and zero denominators are rejected explicitly. The energy ranges are the two
bounds persisted by the current reference code (energyBin01).

Channel time walk, ratio time walk, their energy anchors, and the three separate
phase-profile corrections are replayed exactly from the saved objects. Both
time-walk-only and time-walk-plus-phase stages are recorded. The original code's
behavior for empty/overflow phase bins is retained for comparability and flagged
per event and in the cutflow. Such events warrant checking, not silently dropping
them from just one estimator. Long-long timestamp differences are formed before
conversion to floating point to preserve precision.

For the same accepted events define:

- X: corrected L minus corrected REF average;
- Y: corrected R minus corrected REF average;
- A = (X+Y)/2;
- D = X-Y: the previously missing difference of the externally corrected times;
- Q: the difference corrected separately using energy ratio and mean phase.

The output contains the following **signed variances** (ps²):

| Name | Definition |
|---|---|
| matched_halfdiff | Var(D)/4 |
| direct_halfdiff | Var(Q)/4 |
| catr | [Var(X)+Var(Y)−2 Var(A)]/2 |
| ref_matched | [Var(X)+Var(Y)−Var(D)]/2 |
| ref_mixed | [Var(X)+Var(Y)−Var(Q)]/2 |
| avg_matched | Var(A)−ref_matched |
| avg_mixed | Var(A)−ref_mixed |
| closure_ps2 | catr−matched_halfdiff; must vanish for event variances |
| correction_gain_ps2 | direct_halfdiff−matched_halfdiff |

`avg_mixed` reproduces the **latest summary macro**, which uses the separately
corrected difference in its REF triangulation. This differs from the CatR formula
written in the note. For exact variances, if H=Var(D)/4 and J=Var(Q)/4,
`avg_mixed = 2J−H`. Therefore changes in the separately corrected difference also
change the extracted reference and average resolutions. In particular, a smaller
J than H pushes this mixed extraction below even sqrt(J). This is not an
independent measurement of an improvement in the true DUT mean.

No extraction here establishes that intrinsic DUT/REF correlations are absent.
The reported covariance matrix is that of the **residuals**, which includes the
shared reference, not a measurement of intrinsic DUT-channel covariance.

Uncertainty and fitting
-----------------------

Exact unbinned sample variances use ddof=1 and all accepted events. Their closure
identity is checked numerically. Gaussian fits separately follow the existing
helper's 400-bin distributions, ±5 RMS axis range, and peak ±1 RMS fit interval,
with the original 20 ps seed. Left/right/average use [500,2500] ps and differences
[-1000,1000] ps. Underflow and overflow counts, fit status and covariance status
are saved. Fixed histogram windows can truncate tails differently; the moments
remain untruncated. Failed fits stay invalid; negative extracted variances are
retained, never silently clamped to zero.

Paired event bootstrap resamples all five quantities together, preserving their
correlations. By default 200 moment replicates run; Gaussian refits require
`--fit-bootstrap N` and can be expensive. Only successful complete fit tuples are
used, and their count is reported. A square-root uncertainty is not quoted when
bootstrap radicands cross zero. Derived Gaussian errors are absent unless the
paired-refit bootstrap was requested. Previous fit uncertainties are not assumed
independent of the new measurements; the comparison does not report naive pulls.

These are conditional uncertainties with **saved calibrations fixed**. They do
not include training uncertainty or constitute independent train/test validation.
Independent calibration files can be supplied, but their compatibility and
provenance must be established separately. The script does not claim to refit
REF calibration, optimize thresholds, or estimate the full systematic budget.

Outputs
-------

- `timing_closure.root`: all five residual histograms/fits for both stages,
  flat event trees with the input-tree entry ID and phase-validity flag, provenance.
  Use `--no-event-tree` to reduce storage.
- `summary.json`: settings, ROOT version, input ROOT UUIDs, schema hash, selected
  entry hashes, cutflows, moments, fitted widths, signed derived variances,
  residual covariance matrices, paired bootstrap uncertainties.
- Comparison directory: CSV/JSON with saved previous fits, event counts and
  replay-minus-previous widths, and PNG/PDF overlays grouped by Vov/threshold.

The previous step-2 files have no selected event IDs, so matching histogram counts
alone cannot prove equal samples. Validate replay against the same configuration
before interpreting differences against another run. Legacy outputs opened in
UPDATE mode may contain stale cycles; ROOT's latest cycle is used consistently.

The tests generate deterministic synthetic ROOT inputs with nonzero time-walk
and phase corrections, large coarse timestamps, selection failures, and a direct
correction that intentionally leaves an amplitude-dependent residual. They also
test an empty reference input and refusal to overwrite output. Synthetic inputs
must never be interpreted as test-beam results.
