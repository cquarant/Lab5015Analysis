# Multi-REF timing closure investigation

## Objective and constraints
Investigate REF choice for runs 4194, 4257, 4281, 4291, all 16 DUT bars and six thresholds; preserve the four estimators, TOFHIR-only energy calibration and existing selections. Keep old outputs. Validate REF7 before full production. Report in English. No delegated agents requested.

## Checkpoint 1 — workflow inventory (2026-11-24)
- Baseline: `/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_withref_TOFHIR/{run}` (ROOT, configurations, closure summaries and plots).
- Discovered pre-existing unfinished production: `/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_central_refs_TOFHIR/{run}`. REF4 and REF5 complete for all four runs; REF6 replay failed; REF7/8 partial; REF9–12 absent. Do not mistake this for a completed production.
- Its executables and scratch: `/afs/cern.ch/work/c/cquarant/btl_central_refs_20261124/{bin,work}`. Original run logs show replay failure on REF6. Repository binaries are older and must not be used accidentally.
- Existing uncommitted scripts include `scripts/analyze_central_refs.py`, `main/pool_central_refs.py`, and `tests/test_central_refs.py`. Runner hardcodes REF4–12; pooling hardcodes exclusion of REF6 without documenting a reason. Both need audit before reuse.
- Workflow: `main/moduleCharacterization_step1.cpp`: Coincidence.chL/chR, REF mapping/calibration, REF Landau selection and DUT selection. `main/refCharacterization.cpp`: REF correction trained on DUT bar 7 (`useAmpWalkPerBar=0`), then applied per DUT group. `main/moduleCharacterization_withRef_step2.cpp`: DUT corrections. `main/moduleCharacterization_timingClosure.py`: saved-correction replay and exact moments/Gaussian-core fits, paired bootstrap. `macros/compare_timingClosure.py`: comparisons.
- Mathematical issue: with X=L−REF, Y=R−REF, A=(X+Y)/2, D=X−Y, Q=directly corrected L−R, moment CatR=Var(D)/4 and matched triangularization also equals Var(D)/4 identically. Mixed extraction=Var(Q)/2−Var(D)/4. These are not independent measurements of intrinsic Cov(L,R). Residual covariance includes REF variance/cross-covariances.
- Current Python: PyROOT 6.40.04. Scratch `/tmp` has ~75 GB free.

## Commands executed
Bounded reads/listings of repository scripts, baseline and existing production; `git status --short`; `python3 -c 'import ROOT; print(ROOT.__version__)'`; `ldd bin/moduleCharacterization_step1.exe`; inspected run logs and completed closure summaries.

## Pending / exact next step
1. Inspect REF6 failure and existing binary wrapper; validate existing REF7 files against baseline, including counts, correction objects, moments and replay.
2. Audit calibration coverage and all REF bars' occupancy/quality; define usable selection from evidence, not hardcoded central membership.
3. Extend runner to configurable bars, checked resumability and explicit failures; retain/refit each REF's corrections. Save reusable output under a new multi-REF investigation directory.
4. Aggregate with documented event-overlap handling, centered/un-centered sensitivity, uncertainty and covariance diagnostics. Produce plots/tables/HTML and final physics assessment.

## Checkpoint 2 — baseline verified; important selection difference
- Exact published-input replay PASSED for all 384 groups, 3,750,136 accepted observations. Identical moments, Gaussian widths, counts and source-entry hashes. Validation: `/tmp/btl-multiref-20261124/published_validation/validation.json`.
- Independent existing central-production REF7 is NOT identical: accepted counts 848605 / 1040888 / 1176158 / 776485 for runs 4194 / 4257 / 4281 / 4291. Difference tracked in `/tmp/btl-multiref-20261124/baseline_validation/validation.json` (expected failed equality).
- Archived source hash of step1 matches `git show HEAD:main/moduleCharacterization_step1.cpp`. It lacks the `channelIdx < 0` guard when counting active REF bars, allowing out-of-bounds access for absent channels. Current source fixes it. Keep this fix; compare all REF choices to the NEW safe REF7, and separately quantify changes versus published REF7. Do not silently call them identical baselines.
- Original REF6 replay failures were AFS temporary-dictionary cleanup errors, NOT evidence of an unusable crystal. No exclusion of REF6 is justified by them.
- All 32 REF channel factors at 3 V / threshold 10 are finite positive. Existing source applies them by multiplication.
- New runner `scripts/analyze_multiref.py`: configurable all REF bars, stage receipts, production fingerprint, import of prior matching corrections, individual failure records; defaults to 200 moment / 100 fit bootstraps and two processes.
- Raw occupancy audit `scripts/audit_ref_occupancy.py`: one pass across explicit input_spills.csv, calibrated spectra and pair counts for all 16 REF bars, validates raw entry totals. Initial sandbox attempt had ROOT network errors; `/tmp/btl-multiref-20261124/occupancy` is INVALID, do not use. Rerun with network permission at `occupancy_verified` (tool session recorded in conversation).
- User was asked asynchronously how the independent analysis combines references; no answer yet. Proceed with explicit union of separately calibrated coincidences and overlap diagnostics, do not claim this reproduces an unspecified simultaneous-time estimator.

## Commands / resume
```
python3 scripts/validate_multiref_baseline.py --baseline /eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_withref_TOFHIR --production /eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_central_refs_TOFHIR --output /tmp/btl-multiref-20261124/published_validation --published-inputs --jobs 2
python3 scripts/audit_ref_occupancy.py --baseline /eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_withref_TOFHIR --output /tmp/btl-multiref-20261124/occupancy_verified
```
ROOT compilation and remote data may require escalation (ccache system path / network). Runner sets caches in its separate output. Do not use existing old `bin/*.exe`; use rebuilt wrappers in `/afs/cern.ch/work/c/cquarant/btl_central_refs_20261124/bin`.

## Modified/added files this task
`PROGRESS.md`, `scripts/validate_multiref_baseline.py`, `scripts/audit_ref_occupancy.py`, `scripts/analyze_multiref.py`. Existing dirty files preserved.

## Exact next step
Launch full configurable production after the successful published validation, importing completed/partial central-production ROOT inputs when cfg settings match. All 16 REF bars must be audited; low-population/failing groups must be listed, not silently omitted. Then implement event-identity-aware aggregation and covariance diagnostics, compare to both safe and published REF7, build English HTML report and publish separate outputs on EOS. Save new stage state after each milestone.

## Checkpoint 3 — all-bar production running, aggregation implemented
- User clarified the independent analysis **unites samples from different REF bars**. Primary comparison will do exactly this (separate corrections, uncentered union); no simultaneous averaged-REF-time hypothesis is needed.
- Raw audit completed in `occupancy_verified/{run}.{csv,root}`: all 16 REF bars have data in all four runs. Exact manifest raw entry counts checked. No blanket central-only selection justified.
- Running original production: `/tmp/btl-multiref-20261124/production`, log `production.log`, tool session **64266**. Two jobs. Command used all runs, ref order `7 4 5 6 8 0 1 2 3 9 10 11 12 13 14 15`, prior central production + work imported where matching. Running process holds its original Python code in memory. DO NOT resume this directory with modified runner: its source fingerprint is deliberately old. Copy v1 code snapshots are in `source_v1/`.
- Original imports REF4/5 rejected because their cfg inputDir points to a now-missing local raw cache, not central EOS. Recompute these instead of assuming raw equivalence.
- Some peripheral groups fail time-walk calibration (<3 usable profile bins). Legacy step2 aborts the entire REF sample. Added **opt-in** `--skip-uncalibratable` third CLI argument in `main/moduleCharacterization_withRef_step2.cpp`: marks unsupported group as `TNamed uncalibratable_<label>`, omits its events from later steps; unchanged successful fits and default failure behavior.
- Added opt-in `--skip-unusable` in replay to record missing-correction groups in `skipped`; `main/plot_multiref_single.py` plots valid replay results even if legacy Gaussian fit was absent. Existing comparison macro could fail for this reason (REF0 run4194); its closure is valid and reusable.
- Rebuilt new step2 in `/tmp/btl-multiref-20261124/bin/moduleCharacterization_withRef_step2` with **explicit absolute ROOT 6.36 root-config** and compiling `src/referenceCharacterization_helper.cc` directly. `source scripts/setup.sh` DOES NOT reliably override `/usr/bin/root-config` here (it stays ROOT6.40); previous failed builds left no valid binaries. Successful build uses `/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.36.00/x86_64-almalinux9.5-gcc115-opt/bin/root-config`. New bin wrappers reuse existing step1/ref binaries and use the new local step2.
- Recovery launched separately (one job, alongside original two) into `/tmp/btl-multiref-20261124/production_final`, log `recovery.log`, for run4194 refs `0 1 4 5`, reusing `production`. Tool exec cell56 yields session when waited. This is intended final production directory; after old production completes, import its completed stages for ALL other run/ref combinations into `production_final`. Runner now supports marker-aware import, failed-group recovery and skip flags. Its manifest includes current code; avoid semantic changes to fingerprinted code during the final run.
- Added `main/multiref_statistics.py`, `scripts/extract_multiref_events.py`, `main/aggregate_multiref.py`, `tests/test_multiref_statistics.py`. Flat arrays retain raw integer DUT L/R timestamps, energies, phases, corrected REF difference, both correction stages. Across references, matched event keys must have identical DUT energy+phase; within-reference duplicates fail. Primary union counts coincidences but bootstraps unique DUT event keys (200 moment / 100 Gaussian); paired REF7 differences use same cluster weights. Centered-per-REF and highest-energy-single-REF-per-event alternatives are sensitivity checks (central moments/fits; no fit errors claimed). Covariance diagnostics distinguish residual covariance, algebraic identity, mixed-correction difference and REF-halfdiff-based proxy (not identified intrinsic covariance).
- Extraction of currently completed groups ran in `flat/`, log `extract.log`, session31140. Repeat after more refs finish; completed per-REF receipts skip expensive work.
- Unit tests for clustering, weighted moments, centering and covariance PASSED (2 tests). Original algebra tests PASSED (2). Expanded synthetic ROOT integration tests running in session37689, log `integration.log` (also checks explicit missing-calibration skip vs strict default).
- Aggregation/report generation not yet run. No new EOS outputs published yet.

## Current exact next steps
1. Check sessions 64266 (old production), cell56/recovery, 37689 integration; inspect logs without flooding full ROOT carriage-return output.
2. Validate new step2 skip mode on REF1 recovery and confirm valid-group semantics against safe REF7 (at least selected group rerun), then fill `production_final` using old stages. Keep total executable workers <=4, preferably 2–3.
3. Run extraction on final production (and published baseline separately if needed), then `main/aggregate_multiref.py production_final --flat <flat_final> --output <aggregate>` once each run's 16 refs completed/explicitly failed. Exclusions are per group, not whole bars where recoverable.
4. Build plots, CSVs and English HTML report in existing style, compare safe REF7 and published REF7, document source-level selection bug and covariance identifiability. Archive new outputs separately on EOS and add a separate report page under existing site without replacing old results.

Example final-stage commands (ROOT cache variables should point into /tmp):
```
python3 scripts/analyze_multiref.py --baseline /eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_withref_TOFHIR --validation /tmp/btl-multiref-20261124/published_validation/validation.json --output /tmp/btl-multiref-20261124/production_final --bin-dir /tmp/btl-multiref-20261124/bin --reuse /tmp/btl-multiref-20261124/production --jobs 2
python3 scripts/extract_multiref_events.py /tmp/btl-multiref-20261124/production_final --output /tmp/btl-multiref-20261124/flat_final
python3 main/aggregate_multiref.py /tmp/btl-multiref-20261124/production_final --flat /tmp/btl-multiref-20261124/flat_final --output /tmp/btl-multiref-20261124/aggregate
```

## Checkpoint 4 — recovered from node-local scratch loss

- At continuation, the previous `/tmp/btl-multiref-20261124` and job handles were gone. Reused archived stages and the published validation under `/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_TOFHIR`. Its recovered run 4194 REF1 receipt verifies the skip mode retained 95 groups and omitted one, with 48,067 accepted observations.
- Rebuilt the optional-skip step2 locally with absolute CVMFS ROOT 6.36 `root-config`, explicitly compiling `src/referenceCharacterization_helper.cc`; build passed. Restored wrappers from the central production except for step2, which targets the rebuilt local executable.
- Runner accepts ordered, repeated `--reuse` options. Current production session **40622** on lxplus9104 resumed into `/tmp/btl-multiref-20261124/production_final`; see `resume.log`. It checks inputs in order: EOS `production_final`, EOS `production`, EOS `timing_closure_central_refs_TOFHIR` (with AFS work receipts). REF4/5 with unmatched raw paths are recomputed. By the first check, 13 of 16 run-4194 refs completed and the remaining three progressed through step2 and replay.
- Periodic EOS checkpoint session **7351** saves stage files and receipts once per minute until production status is written. Aggregation session **37440** waits for each run's 16 outcomes, extracts flat event records, and builds aggregate results in order 4194, 4257, 4281, 4291. See `aggregation.log`. Existing checkpoints are skipped. Re-run the archive after aggregation to capture results created after the production status appears.
- Aggregator now also reports moments on the exact intersection of each REF with REF7 when the shared sample has at least 100 events. This holds the selected sample fixed when comparing the separately fitted corrections. Aggregation also includes a central REF 4–12 sensitivity while retaining REF6. The report exports covariance, paired union-versus-REF7 differences and shared-event differences, and includes Methods and Commands pages.
- Statistical unit tests and Python syntax checks passed. Numeric multi-REF conclusions remain pending until the aggregate files exist.

### Resume commands and verification

Set `p=/tmp/btl-multiref-20261124`, `base=/eos/home-c/cquarant/TB_Sept2025_Studies4DN`, and `e=$base/timing_closure_multiref_TOFHIR`. First verify sessions 40622, 7351 and 37440 if available. On a different node, inspect EOS receipts, rebuild the local executable as above, and resume using:

```bash
python3 scripts/analyze_multiref.py --baseline "$base/timing_closure_withref_TOFHIR" --validation "$e/published_validation/validation.json" --output "$p/production_final" --bin-dir "$p/bin" --reuse "$e/production_final" --reuse "$e/production" --reuse "$base/timing_closure_central_refs_TOFHIR" --reuse-work /afs/cern.ch/work/c/cquarant/btl_central_refs_20261124/work --ref-bars 7 4 5 6 8 9 10 11 12 0 1 2 3 13 14 15 --jobs 3
python3 scripts/archive_multiref_checkpoint.py "$p" "$e"
```

For each completed run, run extraction and aggregation with `--runs RUN` as in the preceding checkpoint. Preserve completed EOS aggregate group files when changing nodes. After all 384 group summaries exist, generate the report, derive and add the numeric conclusions, verify tables/plots/static links, and copy the report to the new subdirectory of the existing CERN site and the analysis archive. The report generator expects `published_validation/validation.json` below `--work`; restore that small file from EOS if missing.

## Checkpoint 5 — current scope and persistent recovery

- User selected REF **4, 5, 7, 8, 9, 10, 11, 12**, excluding REF6 and peripherals for known detector performance. This supersedes previous REF lists, raw occupancy criteria, and the redundant central-only sensitivity. All 16 DUT bars, four runs and six thresholds remain. Expected coverage is **32 run–REF samples**, **3,072 individual** and **384 aggregate** DUT-bar/threshold groups before quality exclusions.
- Stopped prior temporary production, aggregation and archive sessions (all terminal). Verified compatible ROOT checkpoints on EOS and some still-visible files from the most recent temporary run. Retained genuine completed products.
- Persistent work: `/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124`. AFS quota had ~21 GiB available at start. Disposable ccache, ROOT dictionaries and Matplotlib caches use a fresh `/tmp/btl-multiref-selected.*`. Planned new EOS archive: `/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_selected_TOFHIR`.
- `scripts/resume_multiref_selected.sh` restores exact published-input validation from EOS, rebuilds the missing step2 executable against ROOT 6.36, runs the selected samples, extracts flat events, aggregates and builds the report. Repeat the same command on failure or a new node. AFS stage receipts and import provenance are persistent. The optional recovery variable is useful only when the preceding scratch directory still exists.
- The same-selection REF7 is independently stored per run in the selected production. Compare against the existing safe REF7 counts and the verified published-input values separately; preserve both in report.
- Synthetic omission smoke run, using actual 4194 REF6/7/8 arrays before narrowing the selection, successfully exercised aggregation, paired bootstrap, and event-key matching (24,060 pairs, 22,949 unique events for DUT7/T1=10). This is a software check only and is excluded from final physics results.

```bash
cd /afs/cern.ch/work/c/cquarant/MTDTB/Simona_TB_CERN_Sep25
bash scripts/resume_multiref_selected.sh /afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124
```

First resumed launch provided `MULTIREF_RECOVER_LOCAL=/tmp/btl-multiref-20261124/production_final` to reuse stage outputs still present. If a process is running, inspect `logs/production-RUN.log`, `logs/extraction-RUN.log`, `logs/aggregation-RUN.log`, and `logs/report.log` in persistent work before attempting restart. The workflow lock prevents a duplicate invocation.

**Next:** verify 32 completion receipts and the new REF7 counts. Finish aggregation, document measured union-minus-REF7 gaps and shared-event results, validate the report, publish separately and recheck artifact coverage before completing the goal.

### Persistent run evidence

- Execution session **90439** uses the standalone command above. Updated scripts compile and the shell syntax check passed. The rebuilt executable and build log are persistent.
- Selected run 4194 REF7 reproduces all 96 groups of the archived safe baseline, with identical counts, source-entry SHA256 and moment results at both correction stages (848,605 accepted). Archived published validation passed for all runs with counts 830,471 / 1,015,166 / 1,148,337 / 756,162.
- Validated the legacy REF4/5 relocated input logs against `input_spills.csv`: exact 3 V filename multisets and total raw entry counts for each run (45/7,381,902; 66/8,634,900; 68/10,829,215; 43/6,535,482). The runner will accept that specific legacy cache path when these checks and all other physics cfg options agree, recording the cross-check in its import receipt. The raw cache itself is unavailable for a byte-level checksum. Run 4194 REF4 independently reproduced the archived accepted count 616,532.
- Added a final verifier for exactly 32 selected sample receipts, exact REF7 fingerprints and fits against the distinct archived safe replay, 384 aggregate summaries, covariance algebra, CSV uniqueness and complete plot/link coverage. Run standalone with:

```bash
python3 scripts/verify_multiref_selected.py /afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124 \
  --safe-baseline /eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_TOFHIR/baseline_validation
```

### Checkpoint 6 — run 4194 ready for aggregate

- Eight selected REF receipts for run 4194 are complete with 96 groups each. REF4 and REF5 were verified by recomputation (616,532 and 836,834 accepted respectively). Event extraction finished for all eight.
- Run 4257 is in production, with REF4/5 using the audited relocated raw inventory. Aggregation worker **67080** runs concurrently with production, using a per-run lock also observed by the standalone workflow. First checked aggregate group (4194 DUT0, T1=7) has 25,365 pairs, 22,861 unique DUT observations and 2,125 observations shared across refs. Work is in the persistent AFS directory.
- To aggregate a specific already extracted run without waiting for the complete standalone sequence:

```bash
p=/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124
OPENBLAS_NUM_THREADS=1 python3 main/aggregate_multiref.py "$p/production_final" \
  --flat "$p/flat_final" --output "$p/aggregate" --runs 4194
```

- The English report and CSV now show safe-REF7 differences for each estimator separately, in addition to the mixed-minus-direct gap. No selection conclusion is drawn from the first incomplete aggregate.

### Checkpoint 7 — first offset diagnostic

- Safe REF7 reproduced exactly for both 4194 and 4257 in all 96 groups, including both correction stages, histogram fit results and event hashes. The final verifier repeats these checks for all four runs.
- First complete 4194 central groups at T1=10 show pronounced per-REF average-residual offsets. For DUT bars 4, 5 and 6, the span is about 0.87–0.89 ns and the between-REF contribution to the average-residual variance is 41,770–45,264 ps². Uncentered union gives Gaussian mixed-minus-direct gaps of invalid, +10.37 ps and +91.82 ps respectively. Centering gives −2.06, +0.05 and −0.16 ps. These are individual groups during production, not the full physics conclusion.
- The report now exports `ref_offsets.csv` and explains how one-Gaussian fits on offset mixtures can respond differently across observables. It also exports the historical 4194 note comparison at fixed and independently optimized thresholds using the existing digitized ratios.
- At last receipt check 4194 had 8/8 per-REF products, 8/8 flat extracts and 48/96 aggregates. 4257 had 4/8 per-REF products with unresolved refs progressing in the live standalone session 90439. No restart indicated. Exact resume commands remain in checkpoint 5.

### Checkpoint 8 — run 4194 aggregate complete

- All 96 selected-REF aggregate groups are present for run 4194. Across DUT bars 4–11, at T1=10 the uncentered Gaussian mixed-minus-direct gap is valid in only 3/8 groups (median +91.82 ps among those), versus 8/8 centered (median −0.07 ps). The exact-moment medians are +10.53 and −0.78 ps respectively. One REF per event without centering still shows a large pooled Gaussian gap and low validity. At T1=3 even the centered Gaussian median is +3.86 ps; for T1=5–20 it ranges between −1.37 and −0.07 ps. Comparison of shared events and external note remains to be summarized.
- 4257 selected production remains live in session 90439; REF9/10/11 replayed and REF12 is the last active raw production. Aggregation 67080 will start 4257 after its flat extraction. Do not restart based on elapsed time.
- Historical comparison for run 4194, T1=10: digitized CatR/direct median 0.9083 (13 bars); centered selected union 1.0057 (13); uncentered union 2.7959 (8 valid). Shared-event change to the moment mixed-minus-direct gap relative to safe REF7: median −0.65 ps across 36 eligible DUT–REF groups at T1=10, and +5.02 ps across 39 at T1=3. The report now renders both summaries and gives the validity denominators.

### Checkpoint 9 — second run extracted

- REF production and flat extraction for 4257 completed 8/8, all per-REF receipts with 96 groups. Its 96-group aggregation is running in session 67080. The initial 4257 T1=10 groups also show ~0.84 ns average-residual offset spans across REF bars.
- Session 90439 is live and is processing 4281; REF4/5 can reuse audited products and REF7 has reached replay. Continue to monitor the same sessions; the persistent working directory and resume command from checkpoint 5 remain authoritative.

### Checkpoint 10 — verify archived replay input UUIDs

- Detected a provenance bug when a legacy REF4/5 sample had a completed replay and original ROOT corrections, while its temporary step1 file/markers had been lost. The first resume could copy the old replay and then silently regenerate both of its inputs. ROOT UUID scan found only 4257 REF5 among the completed samples affected. 4281 REF4/5 were still actively writing regenerated inputs and may require identical recovery after their receipts appear.
- Stopped aggregate worker 67080. Recovery session **51372** restores 4257 REF5's two exact EOS files, matching both recorded replay input UUIDs. It preserves replaced ROOT files, that flat sample and all 4257 partial aggregates under `before_replay_input_recovery/4257`, rebuilds only the affected flat extraction and reaggregates 4257. Run 4194 passed the completed-sample UUID scan and is retained.
- Fixed the importer to recognize the old final `complete` receipt for its existing reference and correction files, regenerating a missing step1 independently. Extraction now rejects mismatched provenance before writing new event arrays. The final verifier checks both replay input UUIDs and the flat-cache UUIDs for every sample.
- The in-memory production of **4281** predates this importer fix. Let it finish the current run; its new extraction is expected to stop at a UUID mismatch. Verify first; then restore affected REF4/5, and resume the standalone command. The repair is a targeted data recovery, not a physics reprocessing:

```bash
p=/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124
central=/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_central_refs_TOFHIR
python3 scripts/restore_multiref_inputs.py "$p" "$central" 4281 4 5
bash scripts/resume_multiref_selected.sh "$p"
```

- Confirm recovery session 51372 is terminal before manually changing its 4257 files. The aggregate per-run lock prevents concurrent recalculation. Intermediate archived configurations and the four estimators have not been altered by the repair.

### Checkpoint 11 — 4257 aggregate revalidated and 4281 extraction

- 4257 REF5 was restored and freshly extracted; all 96 aggregate groups were rebuilt and are present. Central DUT bars, T1=10: Gaussian mixed-minus-direct gap has 3/8 valid uncentered, median +1.69 ps among those, and 8/8 valid centered, median −0.73 ps. Centered Gaussian T1=3 median is +3.08 ps. The moment gap changes by median −1.29 ps across 44 same-event REF–DUT groups at T1=10.
- 4281 REF4/5 were restored in session 99788 by matching both exact UUIDs before extraction. All eight sample receipts are present; UUID-checked extraction is progressing in session 90439. 4281 safe REF7 reproduced all 96 labels, selection counts, hashes, moments, Gaussian results and fit metadata exactly.
- Aggregate monitor **59983** is waiting for 4281 extraction and will then also handle 4291. The persistent standalone workflow will process 4291 and render the report. The exact commands in checkpoint 5 remain sufficient to resume if the session becomes terminal.

### Checkpoint 12 — all safe REF7 baselines verified

- Safe REF7 exact comparison now passed for all four runs (96 groups each, both stages). Accepted counts: 848,605 / 1,040,888 / 1,176,158 / 776,485. This is distinct from the already verified published replay with 3,750,136 events in total.
- All selected samples and flat receipts are complete for 4194, 4257 and 4281; 4194 and 4257 aggregates are complete (192 groups). 4281 aggregation and 4291 production remain active in sessions 59983 and 90439. At checkpoint 4291 REF4/5/7 were complete and no whole sample exclusions were recorded.
- Python and Bash syntax checks, two statistical unit tests, and `git diff --check` pass. Final report validation and publication await all 384 aggregate groups.

### Checkpoint 13 — all selected per-REF samples validated

- All 32 run–REF receipts are complete. Each has exactly 96 DUT/threshold groups (3,072 total; no selected sample or group exclusions). Independently verified exact 3 V spill filename lists and entry totals against the baseline manifest and both current ROOT input UUIDs against the saved replay for every sample.
- 4281 central DUT at T1=10: centered Gaussian mixed-minus-direct gap −0.18 ps across eight bars; at T1=3, +13.31 ps. The T1=3 centered Gaussian gap before phase correction was +3.40 ps, while the exact moment gap was +10.39 ps before phase and +9.27 ps after. Report includes the corresponding correction-stage table for every run.
- 4291 extraction and aggregation remain in the active workflow. Sessions 90439 and 59983 are authoritative when available. Remaining milestones are aggregate/report verification and publication.

### Checkpoint 14 — all flat caches ready, last aggregate running

- All 32 flat extraction receipts are complete. 4194, 4257 and 4281 each have 96 finalized aggregates. 4291 aggregation is running in session **59983**.
- Session **90439** completed all its production and extraction before exiting 2 at the next shell loop. The shell read a stale offset after the on-disk driver was updated during execution. The current driver passes `bash -n`; none of its completed stages need to be rerun.
- EOS archive copy started in session **74790**. A separate report and validation command is queued behind the 4291 aggregate lock. For standalone report generation after all 384 group receipts exist:

```bash
p=/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124
base=/eos/home-c/cquarant/TB_Sept2025_Studies4DN
export MPLCONFIGDIR="$(mktemp -d /tmp/btl-report-mpl.XXXXXXXX)"
python3 scripts/report_multiref.py --work "$p" --baseline "$base/timing_closure_withref_TOFHIR" \
  --output "$p/report" --style "$base/../www/cms_btl_timing_closure/assets/style.css"
python3 scripts/verify_multiref_selected.py "$p" --safe-baseline "$base/timing_closure_multiref_TOFHIR/baseline_validation"
```

### Checkpoint 15 — all 384 groups complete

- Session 59983 completed the selected aggregate. The union contains 18,496,110 pairs, 16,783,116 distinct DUT-bar observations and 1,459,113 observations shared between REF samples. Group-moment CatR identity has maximum absolute residual 1.13e-8 ps² across both stages and all variants.
- Median cross-REF span of average residuals is 847 ps across 384 groups. At T1 ≥ 5, on DUT bars 4–11, 88/160 uncentered and 160/160 centered Gaussian mixed/direct gaps are valid. At T1=10, centered Gaussian medians for 4194/4257/4281/4291 are −0.07 / −0.73 / −0.18 / −1.64 ps. The shared-event T1=10 moment gap changes relative to safe REF7 are −0.65 / −1.29 / −0.70 / −0.70 ps.
- The centered 4194 historical CatR/direct ratio median of 1.006 differs from 0.908 in the digitized DN-26-014 points. The observed uncentered discrepancy is consistent with Gaussian fitting of offset mixtures. The study gives no evidence that REF choice alone reproduces the reported improvement. The residual covariance identity is algebraic, and intrinsic DUT covariance is not identified without independent REF variance and correlation information.
- First EOS checkpoint copy completed. Report construction and full verification are running in session **64192**. The final archive refresh and separate CERN site deployment remain pending.

### Checkpoint 16 — report validated, publication started

- Session **64192** exited 0. `report/data/validation.json` verifies 32 run–REF samples, both ROOT replay input UUIDs per sample, both REF7 baselines, 384 aggregate groups, 18,432 unique result rows, and the complete PNG/PDF asset inventory. JavaScript syntax and actual rendering logic were checked in Node for all 576 run/sample/threshold/estimator selections. Each selection found 16 DUT rows and resolved all linked images and PDFs.
- Final report in `/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124/report`. The separate CERN website and EOS archive refresh have been started. The original single-REF site page is unchanged.

```bash
p=/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124
e=/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_selected_TOFHIR
site=/eos/home-c/cquarant/www/cms_btl_timing_closure/multiref_selected
rsync -rc --delay-updates -- "$p/report/" "$site/"
rsync -rc --delay-updates -- "$p/report/" "$e/report/"
python3 scripts/archive_multiref_checkpoint.py "$p" "$e"
rsync -rcni -- "$p/report/" "$site/"
```

### Checkpoint 17 — complete and published

- Publication session **64388** exited 0. Final checksum-only site comparison reported no differences. The analysis EOS archive contains eight completed samples, 96 aggregate summaries and 96 aggregate ROOT files per run. Its report HTML, validation JSON and all-results CSV match the CERN site copies byte for byte.
- Report URL: https://cquarant.webtest.cern.ch/cms_btl_timing_closure/multiref_selected/index.html . The deployed report directory inherited the existing Web Services EOS ACL. The HTTP request returned 302 to CERN SSO. An authenticated browser session was unavailable; the deployed files and interactive selectors were verified independently.
- Persistent ROOT outputs, reports and logs: `/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_selected_TOFHIR`. Flat event checkpoints and runnable binary are also retained in AFS work. The standalone driver and final verifier commands above remain available to reproduce the work.
- No remaining analysis steps. The findings are in the English report, with the raw union, centered and one-reference-per-event variants, both REF7 baselines and all four estimators.

### Checkpoint 18 — native step2 diagnostics recovery

- The multi-REF archive had retained `step2.root` for every sample. The report plots are separate from the original step2 PNGs. Collecting the stock diagnostic histograms and profiles for all 32 validated samples, reusing compatible original PNGs or redrawing stored ROOT objects and fit functions. Each output archive has a source UUID receipt. The four timing-resolution estimators are independent of this operation.
- Standalone recovery from the repository root:

```bash
p=/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124
e=/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_selected_TOFHIR
stage=/tmp/cquarant/step2-diagnostics-batch
bash scripts/collect_step2_diagnostics.sh "$p/production_final" "$stage" 6
rsync -rc --delay-updates --include='*/' --include='step2_diagnostics.tar.gz' --include='step2_diagnostics.json' --exclude='*' "$stage/" "$e/production_final/"
rsync -rcni --include='*/' --include='step2_diagnostics.tar.gz' --include='step2_diagnostics.json' --exclude='*' "$stage/" "$e/production_final/"
```

### Checkpoint 19 — step2 diagnostics published and verified

- All 32 `production_final/<run>/ref<NN>/step2_diagnostics.tar.gz` archives and their JSON receipts are on EOS. Each archive contains 3,936 unique PNGs with the standard quantity directories. The total is 125,952 PNGs. The 31,488 PNGs for run 4194 were redrawn using the saved histograms, profiles and fit functions. The other 94,464 are original PNGs.
- The checksum-only rsync check showed no differences. Each of the 32 EOS receipts matches the step2 UUID in its EOS closure summary. Representative restored and redrawn 1D, 2D and profile PNGs passed image validation. No step2 or closure recomputation was necessary.

### Checkpoint 20 — centered comparison report

- The four-method comparison overlays the centered union and REF7 with the corrected selection. The gap figure compares the centered mixed-minus-direct gap with both REF7 baselines, and the lower panel subtracts the current REF7 gap bar by bar. The centered results use point estimates. The index opens with the centered sample selected.
- Regenerate the report using the checkpoint 14 generator command with an output directory under `/tmp`. Use the original website style path so the digitized note table remains available. The analysis inputs and all-resolutions CSV are independent of the plot formatting.
- Publish these outputs from the generated report to the analysis archive report and the CERN website report directory. The selected files are `index.html`, `methodology.html`, `app.js`, `plots/delta_*.{png,pdf}`, and `plots/methods_centered_vs_ref7_*.{png,pdf}`.
- Published and verified both EOS copies. All 576 selector combinations resolved to 16 DUT rows and the expected PNG/PDF assets. All 192 comparison files and the three updated report documents matched the staged output by checksum. The analysis CSV and browser data matched the previously published copies byte for byte.

```bash
p=/afs/cern.ch/work/c/cquarant/btl_multiref_selected_20261124
base=/eos/home-c/cquarant/TB_Sept2025_Studies4DN
site=/eos/home-c/cquarant/www/cms_btl_timing_closure/multiref_selected
stage=$(mktemp -d /tmp/btl-centered-report.XXXXXXXX)
export MPLCONFIGDIR="$stage/mpl"
python3 scripts/report_multiref.py --work "$p" --baseline "$base/timing_closure_withref_TOFHIR" \
  --output "$stage/report" --style "$site/../assets/style.css"
filters=(--include='/index.html' --include='/methodology.html' --include='/app.js' \
         --include='/plots/' --include='/plots/delta_*.png' --include='/plots/delta_*.pdf' \
         --include='/plots/methods_centered_vs_ref7_*.png' --include='/plots/methods_centered_vs_ref7_*.pdf' --exclude='*')
for dest in "$base/timing_closure_multiref_selected_TOFHIR/report" "$site"; do
    rsync -rc --delay-updates "${filters[@]}" "$stage/report/" "$dest/"
    rsync -rcni "${filters[@]}" "$stage/report/" "$dest/"
done
```

### Checkpoint 21 — five-page multi-REF report and diagnostic viewer

- The overview maps all four plot labels to the timing estimators. Results, Step 1, REF aggregation and Diagnostics have dedicated pages. The aggregation page uses the existing before/after plots. The diagnostic viewer reads original PNG bytes for each run/REF/DUT/threshold/quantity. It displays the corresponding replay widths after time walk and after phase correction for Gaussian fits and exact moments.
- Verified 125,952 PNG entries in 512 per-bar bundles, including the exact historical example quantity across all 32 samples. Each sample contains all 41 quantities at all six thresholds and 16 DUT bars.
- Reproduce this presentation from an existing report and the step2 archives:

```bash
archive=/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_multiref_selected_TOFHIR
site=/eos/home-c/cquarant/www/cms_btl_timing_closure/multiref_selected
stage=$(mktemp -d /tmp/btl-report-pages.XXXXXXXX)
python3 scripts/package_multiref_diagnostics.py "$archive/production_final" "$stage/diagnostic-images"
python3 scripts/organize_multiref_pages.py "$archive/report" "$stage"
for dest in "$site" "$archive/report"; do
    rsync -rc --delay-updates "$stage/diagnostic-images/" "$dest/diagnostic-images/"
    rsync -rc --delay-updates "$stage/"*.html "$stage/diagnostics.js" "$dest/"
    rsync -rcni "$stage/" "$dest/"
done
```
- Published the five pages and original diagnostic image bundles on the CERN website and in the EOS report archive. All 576 Results combinations resolved to 16 DUT rows and their PNG/PDF links. Exercised the deployed diagnostics viewer for the requested 4194/REF7/DUT7/T1=10 final Gaussian example and 4281/REF12/DUT15/T1=3 time-walk moments. Both returned valid PNG bytes and matched the recorded stage widths.
