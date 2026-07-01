---
phase: 04-npe-training-advi-benchmark-ablation
plan: 01
subsystem: testing
tags: [benchmarktools, neuralestimators, random123, jld2, reproducibility, spike-env, pre-registration]

# Dependency graph
requires:
  - phase: 03-training-data-pipeline
    provides: "Loader.load_holdout / holdout.jld2 (reserved ≥20-stack set, D-10), keyed holdout_rng seeding, encode_d01 summary, generate_cache fixture builder"
  - phase: 01-environment-smoke-gate
    provides: "pinned spike env (NeuralEstimators v0.2.1) + UUID-keyed resolve-risk gate pattern in runtests.jl"
provides:
  - "BenchmarkTools v1.8.0 as a direct spike-env dep with the v0.2.1 pin provably intact (Turing still absent)"
  - "spike/npe/resimulate.jl — shared keyed holdout re-simulation path (resimulate_holdout / resimulate_holdout_pair) for run_advi.jl (04-04) and benchmark.jl (04-05)"
  - "spike/test/test_npe.jl — six pre-registered Phase-4 constants + SC1..SC5 skipped placeholders + a passing Wave-0 holdout reproducibility gate"
  - "Open Question 1 resolved: reserved holdout is byte-reproducible from (θ, global_index, master_seed) — no raw-image persistence needed"
affects: [04-02, 04-03, 04-04, 04-05, 04-06, 04-07, phase-5-sbc-bf-ood]

# Tech tracking
tech-stack:
  added: [BenchmarkTools v1.8.0]
  patterns: ["Phase-4 UUID-keyed resolve-risk gate extension", "shared keyed re-simulation helper", "pre-registered constants committed before any reported run", "bit-exact reproducibility round-trip gate"]

key-files:
  created:
    - spike/npe/resimulate.jl
    - spike/test/test_npe.jl
  modified:
    - spike/Project.toml
    - spike/Manifest.toml
    - spike/test/runtests.jl

key-decisions:
  - "resimulate_holdout recovers the holdout key from the STORED negative global_index (global_index[j] = -key) rather than assuming an identity mapping, so re-sim tracks whatever ordering the cache wrote"
  - "Wave-0 repro gate asserts bit-exact `==` (never `≈`): a mismatch is the signal to persist raw holdout images, not to loosen the tolerance"
  - "Turing/AdvancedVI/ForwardDiff asserted ABSENT from the spike env in runtests.jl gate (f) — they live only in the isolated spike/baseline/ env (D-01), never co-resolving the NeuralEstimators pin down"

patterns-established:
  - "resolve-risk gate (f): after each new spike dep, re-assert NeuralEstimators UUID 38f6df31-… == v0.2.1 with the new dep present"
  - "single shared re-sim path so the ADVI baseline and the NPE benchmark can never time divergent raw stacks"

requirements-completed: []

# Metrics
duration: 15min
completed: 2026-07-01
---

# Phase 4 Plan 01: Phase-4 Spike Foundations Summary

**BenchmarkTools added to the pinned spike env (v0.2.1 held), a shared keyed holdout re-simulation helper, and a test_npe.jl scaffold whose Wave-0 gate proves the reserved holdout is byte-reproducible (Open Question 1 resolved) with six pre-registered constants locked before any run.**

## Performance

- **Duration:** ~15 min
- **Started:** 2026-07-01T19:56Z (approx)
- **Completed:** 2026-07-01T20:11Z
- **Tasks:** 2
- **Files modified:** 5 (2 created, 3 modified)

## Accomplishments
- Added `BenchmarkTools v1.8.0` as a direct dep of the isolated spike env and re-froze `spike/Manifest.toml`; the `Pkg.add` did **not** move the `NeuralEstimators v0.2.1` pin (RESEARCH Pitfall 2), and `Turing` remains absent.
- Extended the `runtests.jl` resolve-risk gate with block **(f)**: asserts BenchmarkTools present, Turing absent, and re-asserts the v0.2.1 pin (UUID-keyed) — then wired `include(test_npe.jl)` so a single `runtests.jl` enumerates all Phase-4 SCs.
- Built `spike/npe/resimulate.jl` — `resimulate_holdout` / `resimulate_holdout_pair` replay the disjoint `holdout_rng(master_seed, key)` stream (D-10) to rebuild byte-identical raw `MultiChannelImage` stacks; this is the single shared re-sim path both `run_advi.jl` (04-04) and `benchmark.jl` (04-05) will consume.
- Created `spike/test/test_npe.jl` with the six pre-registered constants (`NPE_RMSE_TOLERANCE=1.2`, `SPEEDUP_GATE=100.0`, `ABL_REL_MARGIN=0.05`, `ABL_FOLD_CONSISTENCY=4`, `BENCH_THREADS=1`, `NPE_MASTER_SEED=0xC0FFEE`), SC1..SC5 named `@test_skip` placeholders, and a **passing** Wave-0 holdout reproducibility gate.
- **Open Question 1 RESOLVED:** `encode_d01(patch_summary(resimulate_holdout(...).mci_sample)) == load_holdout(...).summary_min[:, j]` holds bit-exact for the probed entries — no raw-image persistence fallback needed.

## Task Commits

Each task was committed atomically:

1. **Task 1: Add BenchmarkTools, extend resolve-risk gate, wire test_npe.jl** — `4001de5` (chore)
2. **Task 2: Keyed holdout re-sim helper + test_npe.jl SC scaffold + repro gate** — `67953e5` (feat)

## Files Created/Modified
- `spike/npe/resimulate.jl` (created) — `resimulate_holdout(dir, j; master_seed)` and `resimulate_holdout_pair(...)`; guarded includes of contract/prior/forward/seeding/encode/loader; uses `holdout_rng` (not `sample_rng`) for the reserved set.
- `spike/test/test_npe.jl` (created) — six pre-registered constants, SC1..SC5 skipped placeholders, the bit-exact Wave-0 holdout reproducibility gate + a bad-index `@test_throws ArgumentError` guard.
- `spike/Project.toml` (modified) — `BenchmarkTools` added to `[deps]`.
- `spike/Manifest.toml` (modified) — re-frozen with `BenchmarkTools v1.8.0` (+ `Profile`); NeuralEstimators unchanged at v0.2.1.
- `spike/test/runtests.jl` (modified) — resolve-risk gate (f) + `include(test_npe.jl)`.

## Decisions Made
- **Key recovery from stored index:** `resimulate_holdout` reads the holdout key back out of the stored negative `global_index` (`global_index[j] = -key`) instead of assuming `key == j`, making the helper robust to any holdout ordering the cache actually wrote.
- **Bit-exact, not approximate:** the repro gate uses `==` deliberately; per the plan, a failure must surface (signal to persist raw images), not be masked by a tolerance.
- **Boundary lock:** gate (f) also asserts `Turing` is *absent* from the spike project deps, encoding the D-01 isolation (Turing lives only in `spike/baseline/`) as an executable invariant.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- The worktree spawned at `0d9b87c` (v2.0-baseline), where `spike/` does not yet exist; the mandatory `<worktree_branch_check>` reset the worktree to the intended base `e9f376b`, at which the tracked `spike/` tree is present. Execution proceeded normally after the reset.

## Verification Evidence
- Task 1: `julia --project=spike -e '…'` → `OK` (BenchmarkTools present, Turing absent, NeuralEstimators == v0.2.1).
- Task 2: `julia --project=spike spike/test/test_npe.jl` → 11 pass, 5 broken (skipped), repro gate bit-exact.
- Full suite: `julia --project=spike spike/test/runtests.jl` green end-to-end — CPU smoke 11/11 (incl. gate (f)), SIM 1/2/3/4 all green, Phase-3 pipeline 69/69, Phase-4 scaffold 11 pass + 5 skipped.
- Manifest diff confirms `BenchmarkTools v1.8.0` added and `NeuralEstimators` still v0.2.1 (CUDA remains a weakdep-only extension).

## Known Stubs
The five SC1..SC5 `@test_skip true` testsets are intentional Wave-0 placeholders (not defects): the plan's explicit deliverable is to enumerate every Phase-4 success criterion as a skipped gate that later waves replace. Each is tagged in-file with the requirement (NPE-01/02/03, ABL-01/02) and the wave that fills it. No data-flow stub or hardcoded UI value exists.

## Next Phase Readiness
- BenchmarkTools is available for the >100× wall-clock claim (NPE-03); the pin is provably intact.
- `resimulate_holdout` is ready for the ADVI baseline (04-04) and the benchmark (04-05) to re-simulate identical raw stacks.
- The six pre-registered constants are committed before any reported run, discharging the "tune until calibrated" data-snooping blocker for Phase 4's own thresholds.
- SC1..SC5 are enumerated in the single `runtests.jl` gate, ready to be filled by subsequent Phase-4 waves.

---
*Phase: 04-npe-training-advi-benchmark-ablation*
*Completed: 2026-07-01*
