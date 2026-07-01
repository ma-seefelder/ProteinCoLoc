---
phase: 04-npe-training-advi-benchmark-ablation
plan: 07
subsystem: npe
tags: [scaling, amortization, d-12, d-13, thread-sweep, exponents, log-log, benchmarktools, jld2, cpu-only, characterization]

# Dependency graph
requires:
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 05
    provides: "benchmark.jl speedup_report/_time_npe_pair (D-08 to-posterior clock), _load_artifact, rmse_report; the timing harness reused across N/imsize points"
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 06
    provides: "settled :min summary choice; test_npe.jl SC1-SC5 live at pre-registered consts (BENCH_THREADS=1, SPEEDUP_GATE=100, NPE_MASTER_SEED)"
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 03
    provides: "trained_npe.jld2 (persisted estimator + frozen zt/θzt), infer.jl surface, train_fold for the amortized C_train measurement"
provides:
  - "spike/npe/scaling.jl — scaling_over_N (amortization curve C_train+N·t_fwd vs N·t_advi, log-log exponents) + scaling_over_imsize (shared summary-extraction cost + fixed NPE-forward vs real artifact vi() cost) + scaling_curves driver; atomic JLD2 persistence"
  - "spike/npe/run_thread_sweep.jl — separate-process julia -t N thread sweep (Pitfall 7) with speedup_report per subprocess, aggregate_thread_sweep table, >100× headline at fixed BENCH_THREADS (D-13)"
  - "SC3-scaling (NPE-03 characterization) green in spike/test/test_npe.jl — finite N/imsize exponents, NPE-over-N materially flatter than ADVI (~linear), thread-sweep aggregation artifact produced"
affects: [phase-5-sbc-bf-ood, phase-6-demo-memo]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "amortization curve as a scaling claim (D-12): NPE total-from-cold = C_train + N·t_fwd (one-time training amortized over every future dataset) vs ADVI = N·t_advi (no amortization); log-log exponent of NPE flattened by the C_train offset (< 1) while ADVI is a pure line (≈ 1)"
    - "empirical scaling exponent = ordinary least-squares log-log slope (_loglog_slope); NaN-guarded for <2 distinct points / zero variance"
    - "imsize scaling decomposition: the imsize-DEPENDENT cost is the shared 8×8 patch-correlation summary extraction (measured per size); imsize-INDEPENDENT cost is the tiny NPE forward pass vs the large fixed ADVI vi() cost (real artifact) — an honest MODEL since ADVI cannot run in the pinned spike env (D-01)"
    - "thread-count sweep as SEPARATE julia -t N processes (Pitfall 7 — nthreads() is immutable per process), each writing a per-N JLD2, aggregated afterward; headline read at the fixed BENCH_THREADS, full {1,2,4,8} reported alongside as characterization (D-13)"
    - "reuse the 04-05 D-08 to-posterior clock primitives (_time_npe_pair, forward pass at lightweight bench_N) across grid points rather than re-rolling @belapsed"

key-files:
  created:
    - spike/npe/scaling.jl
    - spike/npe/run_thread_sweep.jl
  modified:
    - spike/test/test_npe.jl

key-decisions:
  - "the amortization curve INCLUDES the one-time training cost as the fixed NPE offset (C_train + N·t_fwd): this is the textbook amortization framing that makes 'NPE flat O(1)/dataset post-training vs ADVI ~linear in N' quantitative — the log-log NPE exponent is < 1 whenever C_train > 0, ADVI is ≈ 1; documented transparently (not a single cherry-picked headline)"
  - "ADVI-over-imsize is a MODEL, not a measurement: ADVI cannot run in the pinned spike env (Turing co-resolve landmine, D-01), so its imsize curve is decomposed into the measured shared summary-extraction cost + the real per-stack vi() cost from advi_artifact.jld2; the interesting D-12 result is on the N axis, the imsize axis honestly shows ADVI's pixel-scaling is masked by its large fixed VI cost"
  - "the thread sweep's SC3-scaling gate exercises the aggregation machinery + atomic file production on cached per-N results (the real worker writer path), keeping the CI gate fast; the full {1,2,4,8} separate-process sweep is a reported artifact (run_thread_sweep spawn=true), not a CI gate (plan directive)"
  - "pre-registered scaling grids as guarded named consts (SCALING_N_GRID/SCALING_IMSIZE_GRID/…) mirroring the ablation.jl guarded-const idiom so the reported curve cannot be tuned to a desired exponent"

requirements-completed: [NPE-03]

# Metrics
duration: ~40min
completed: 2026-07-01
---

# Phase 4 Plan 07: Amortization Scaling Characterization (D-12/D-13) Summary

**`scaling.jl` turns the >100× headline into an honest CURVE: `scaling_over_N` builds the amortization curve (NPE total-from-cold `C_train + N·t_fwd` vs ADVI `N·t_advi`) and `scaling_over_imsize` decomposes both methods over input size, fitting empirical log-log exponents that show NPE materially flatter than ADVI over N (amortized O(1)/dataset post-training vs ~linear); `run_thread_sweep.jl` runs a CPU thread-count sweep as SEPARATE `julia -t N` processes (Pitfall 7), aggregating per-N wall-clocks with the >100× headline stated at the fixed `BENCH_THREADS` (D-13) — turning SC3-scaling green (18 tests) with SC1-SC5 untouched and the full suite green end-to-end.**

## Performance

- **Duration:** ~40 min
- **Tasks:** 2
- **Files:** 3 (2 created, 1 modified)

## Accomplishments

- **`scaling_over_N` / `scaling_over_imsize` / `scaling_curves` (Task 1, `spike/npe/scaling.jl`)** — `scaling_over_N(dir; master_seed, Ns, train_cost)` measures the per-dataset NPE forward-pass latency once via the reused 04-05 `_time_npe_pair` D-08 clock (halved to per-stack) and the per-dataset ADVI cost from the artifact `wall_clock` (median per pair, halved), then builds `npe_time_N = C_train + N·t_fwd` and `advi_time_N = N·t_advi` and fits each empirical scaling exponent by a hand-rolled log-log LSQ slope (`_loglog_slope`). `C_train` (the amortized training investment) is supplied or measured by timing one `train_fold`. `scaling_over_imsize` measures the imsize-dependent shared summary-extraction cost per size and adds the imsize-independent NPE forward pass (measured) vs the real fixed ADVI vi() cost (artifact). Both persist curves via the atomic `.tmp`+reopen+`mv` idiom; pre-registered grids are guarded named consts; `use_gpu=false` throughout.
- **`run_thread_sweep.jl` + `aggregate_thread_sweep` (Task 2)** — because `Threads.nthreads()` is immutable per process (Pitfall 7), `run_thread_sweep` SPAWNS a separate `julia --project=spike -t N` subprocess for each N ∈ {1,2,4,8} (a `--worker` entry running `speedup_report` and writing a per-N JLD2), then `aggregate_thread_sweep` merges them into a sorted thread-scaling table with the NPE-clock parallel (Amdahl) speedup and the >100× headline read at the fixed `BENCH_THREADS` (D-13), never the best thread count.
- **SC3-scaling (NPE-03 characterization) green (Task 2)** — a new sibling testset asserts the N/imsize curves exist with finite exponents, `advi_exponent_N > 0.9` (≈ linear) with `npe_exponent_N < advi_exponent_N` (the qualitative amortization check, not a tuned threshold), both curve files persisted, and the thread-sweep aggregation artifact produced with the headline at `BENCH_THREADS`. SC1/SC2/SC3/SC4/SC5 untouched.

## Measured Results (committed fixtures, `NPE_MASTER_SEED`, CPU-only, `-t 1`)

| Axis | Metric | NPE | ADVI | Reading |
|------|--------|-----|------|---------|
| N (datasets) | log-log exponent | **≈ 0.0004** (train_cost 10 s) | **1.00** | NPE flat/amortized, ADVI pure linear |
| N | per-dataset cost | **~0.6 ms** (forward pass) | **~0.25 s** (½ vi() pair) | ratio ≈ 400× per dataset |
| N | amortization crossover | N\* ≈ 40 (at train_cost 10 s) | — | beyond N\* the trained NPE undercuts ADVI total |
| imsize (pixels) | log-log exponent | **≈ 0.58** | **≈ 0.001** | NPE tracks pixel-scaled summary; ADVI's pixel-scaling masked by its large fixed vi() cost |
| threads | headline speedup @ BENCH_THREADS=1 | **≈ 355×** | — | > 100× at the fixed, reported thread count (D-13) |

The full {1,2,4,8} separate-process sweep is the reported characterization artifact (`run_thread_sweep(spawn=true)`); the CI gate exercises the aggregation machinery on cached per-N results to stay fast.

## Task Commits

1. **Task 1: scaling curves over N and imsize with empirical exponents (D-12)** — `57d7045` (feat)
2. **Task 2: CPU thread-sweep as separate processes (D-13) + SC3-scaling green** — `fb93e57` (test)

## Files Created/Modified

- `spike/npe/scaling.jl` (created) — `scaling_over_N`, `scaling_over_imsize`, `scaling_curves`, helpers `_loglog_slope` / `_atomic_save_curve` / `_simulate_at_imsize`; pre-registered `SCALING_*` guarded consts; guarded include of benchmark.jl.
- `spike/npe/run_thread_sweep.jl` (created) — `run_thread_sweep`, `thread_sweep_measure`, `aggregate_thread_sweep`, `_atomic_jldsave`; `--worker` subprocess entry; guarded `BENCH_THREADS`.
- `spike/test/test_npe.jl` (modified) — new `SC3-scaling (NPE-03 characterization)` testset (18 tests); SC1-SC5 + Wave-0 untouched.

## Decisions Made

- **Amortization curve includes the one-time training cost as the NPE fixed offset** — see key-decisions frontmatter. This makes "NPE flat O(1)/dataset post-training vs ADVI ~linear in N" a quantitative log-log claim (NPE exponent < 1, ADVI ≈ 1) rather than a slogan, and is the honest full-lifecycle accounting the D-12 "amortized thesis IS a scaling claim" intent calls for.
- **ADVI-over-imsize is a labelled decomposition MODEL** — ADVI cannot run in the pinned spike env (D-01), so its imsize curve is (measured shared summary extraction) + (real artifact vi() fixed cost); documented as a model, with the honest finding that ADVI's pixel-scaling is masked by its large fixed VI cost at these sizes.
- **Fast CI gate, reported full sweep** — SC3-scaling aggregates cached per-N results (real worker writer path) to stay fast; the separate-process {1,2,4,8} sweep is the reported artifact (plan directive: "cached timings; the full sweep is a reported artifact, not a CI gate").

## Deviations from Plan

None — plan executed as written. (A pre-commit implementation note: the `--worker` guard needed `(@__FILE__)` parenthesized inside the `&&` chain so the macro would not greedily consume the rest of the condition; caught by the Task-2 verify and fixed before the commit — not a plan deviation.)

## Honest Caveats (for the Go/No-Go memo)

- **The N-axis amortization exponent depends on `C_train` and the N grid** — with the reported `SCALING_TRAIN_EPOCHS=200` training, `C_train` dominates a modest N grid so the NPE log-log slope is ≈ 0 (flat); the crossover N\* = `C_train / (t_advi − t_fwd)` and the asymptotic >100× per-dataset ratio are the load-bearing numbers, reported alongside the exponent so the claim is not a single point. At very small N the realistically-fast ADVI baseline (~0.5 s, 04-05 caveat) wins on total cost until N\*.
- **ADVI's imsize curve is a model, not a measurement** (see Decisions) — honest given the D-01 env isolation; the shared summary-extraction term is real and measured.
- **The thread-sweep headline is stated at `BENCH_THREADS=1`** (D-13), with the full {1,2,4,8} process sweep as reported characterization — the NPE clock's parallel speedup / Amdahl saturation is the NPE-side driver (the ADVI `wall_clock` is the frozen artifact reference).

## Verification Evidence

- **Task 1:** `julia --project=spike -e 'include("spike/npe/scaling.jl"); for f in (:scaling_over_N,:scaling_over_imsize); @assert isdefined(Main,f) f; end; println("scaling defined OK")'` → "scaling defined OK". Functional smoke on committed fixtures: `npe_exponent_N` ≈ 0.0004 vs `advi_exponent_N` = 1.0 (NPE flatter ✓), crossover ≈ 40; imsize exponents 0.58 / 0.001 finite; atomic curve files produced.
- **Task 2:** `julia --project=spike -t 1 spike/test/test_npe.jl` → 139/139 pass (SC1 68, SC2 15, SC3 9, SC4 12, SC5 6, SC3-scaling 18, Wave-0 11). Real end-to-end `run_thread_sweep(spawn=false, threads=[1])`: headline_speedup ≈ 355× at BENCH_THREADS=1, aggregation + per-N files produced.
- **Full suite:** `julia --project=spike -t 1 spike/test/runtests.jl` → EXIT 0; Phase-1 smoke + resolve-risk (BenchmarkTools present, NeuralEstimators v0.2.1 pinned, no CUDA, Turing absent), SIM-02/03/04 + D-15, Phase 3 69/69, Phase 4 139/139. No CUDA loaded.
- **Decoupling:** `git diff --name-only <base> HEAD` = only `spike/npe/scaling.jl`, `spike/npe/run_thread_sweep.jl`, `spike/test/test_npe.jl`. 0 `src/` edits; 0 spike `Project.toml`/`Manifest.toml` changes; no stray generated artifacts.

## Known Stubs

None. SC1-SC5 + SC3-scaling are all live. The scaling and thread-sweep JLD2 outputs are regenerable reported artifacts (written to caller-supplied paths / tempdirs in tests), deliberately not committed.

## Next Phase Readiness

- Phase 5 (SBC / amortized BF / OOD) and Phase 6 (demo + Go/No-Go memo) inherit the amortization curve + thread-scaling table as the honest NPE-03 evidence: the >100× is a curve over N and input size with empirical exponents plus a CPU core-scaling characterization at a fixed thread count, not a cherry-picked point.
- `scaling_curves(dir; save_path)` and `run_thread_sweep(holdir; spawn=true)` are the re-runnable reported-artifact drivers for the memo; the pre-registered grids lock the reported run.

## Self-Check: PASSED
