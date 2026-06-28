---
phase: 03-training-data-pipeline
plan: 02
subsystem: data-pipeline
tags: [julia, random123, statsbase, threads, reproducibility, encoding, spike, decoupling]

# Dependency graph
requires:
  - phase: 03-01
    provides: JLD2 + Random123 in the spike env, resolve-risk gate, the five-SC MISSING test scaffold (SC-1 + D-11/D-12 placeholders this wave fills)
  - phase: 02-simulator
    provides: UNCHANGED sample_prior / simulate_pair / build_mci / patch_summary chain reached read-only via contract.jl
provides:
  - encode_d01 (LOCKED 128-dim D-01 vector = 64 imputed correlations + 64 binary mask)
  - encode_aug (D-02 augmented superset, AUG_DIM=142, + N_AUG_MOMENTS=14 scalar moments)
  - Philox4x keyed per-sample RNG (sample_rng / holdout_rng / fold_rng) + HOLDOUT_SALT/FOLD_SALT
  - cost-aware sample_imsize over IMSIZE_SET/IMSIZE_WEIGHTS (E[cost]≈4.68×, ≥1024² capped 10%)
  - generate_sample / generate_samples (RAW column-major θ/summary_min/summary_aug, thread==serial byte-identical)
  - SC-1 + D-11/D-12 testsets turned green in test_data_pipeline.jl
affects: [03-03 cache, 03-04 loader, 03-05 holdout, phase-4 summary-net input dim=128]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "per-sample Philox4x keyed RNG (master_seed, global_index) → order/thread-independent generation"
    - "guarded isdefined includes so a generator module loads both standalone and inside the runtests harness"
    - "column-major (d×K) RAW training matrices written per-column under Threads.@threads (no shared mutable state)"

key-files:
  created:
    - spike/data/encode.jl
    - spike/data/seeding.jl
    - spike/data/generate.jl
  modified:
    - spike/test/test_data_pipeline.jl
    - spike/NOTES.md

key-decisions:
  - "N_AUG_MOMENTS fixed at 14 (AUG_DIM=142): Manders M1/M2, whole-image Pearson, patch-grid median/IQR/mean/std/skewness/excess-kurtosis, fraction-missing, per-channel intensity median+IQR — Claude's discretion under D-02 (only 'cache D-01 + a superset' is locked)."
  - "IMSIZE_WEIGHTS = (0.55,0.35,0.05,0.03,0.02): ≥1024² fraction capped at 0.10, expected per-sample cost ≈4.68× the 256² baseline — the size-robustness (D-03) vs CPU-budget (D-06) / memory (D-12) lever."
  - "Generator includes are GUARDED with isdefined so contract.jl/forward.jl/prior.jl are not double-included (and structs not redefined) when test_simulator.jl already loaded them in runtests.jl — robust across the standalone and harness entry points (deviation Rule 3)."
  - "Byte-identity asserted with isequal (NaN-safe) rather than ==, so a degenerate (NaN) augmented moment could never spuriously fail the parallel==serial gate."

requirements-completed: [DATA-01]

# Metrics
duration: 12min
completed: 2026-06-28
---

# Phase 3 Plan 02: DATA-01 Generation Core Summary

**Three spike-local modules turn the calibrated prior + UNCHANGED Phase-2 simulator + frozen 8×8 summary contract into RAW, column-major training vectors in memory: `encode.jl` (locked 128-dim D-01 + AUG_DIM=142 D-02 superset), `seeding.jl` (Random123 Philox4x per-sample keyed RNG + disjoint salts + cost-aware imsize sampler), and `generate.jl` (order- and thread-independent generator whose parallel output is byte-identical to serial, proven across -t 1 and -t 4).**

## Performance

- **Duration:** ~12 min
- **Started:** 2026-06-28T09:38:35Z
- **Completed:** 2026-06-28T09:51:01Z
- **Tasks:** 3
- **Files modified:** 5 (3 created, 2 modified)

## Accomplishments

- **`spike/data/encode.jl`** — `encode_d01(M)` produces the LOCKED 128-dim vector (rows 1:64 = `vec(coalesce.(M,0.0))` imputed correlations in column-major order; rows 65:128 = `vec(Float64.(.!ismissing.(M)))` binary present/absent mask). A fully-missing 8×8 is KEPT (vals 0, mask 0), mirroring `induced_mu`'s NaN-not-throw degeneracy choice (D-13). `encode_aug(mci, M)` returns the `AUG_DIM=142` superset = `vcat(encode_d01(M), moments)` with `N_AUG_MOMENTS=14` empty-guarded scalar moments. `N_AUG_MOMENTS`/`AUG_DIM` are named constants so the layout is a contract.
- **`spike/data/seeding.jl`** — `sample_rng(master_seed, idx) = Philox4x(UInt64, (UInt64(master_seed), UInt64(idx)))` (an `AbstractRNG`, drops into the simulator unchanged); `holdout_rng`/`fold_rng` use `HOLDOUT_SALT`/`FOLD_SALT` XOR-salted disjoint key namespaces (D-10). `sample_imsize` draws from `IMSIZE_SET = {256²,512²,1024²,1376×1028,2048²}` under cost-aware `IMSIZE_WEIGHTS` (≥1024² capped 10%, expected cost ≈4.68×).
- **`spike/data/generate.jl`** — `generate_sample(master_seed, idx)` runs the full keyed chain (`sample_rng → sample_prior → sample_imsize → simulate_pair → build_mci → patch_summary → encode_d01/encode_aug`); `generate_samples(N; …)` fills pre-allocated column-major `theta`(7×N)/`summary_min`(128×N)/`summary_aug`(142×N) + `global_index`/`imsize` vectors, dispatching `Threads.@threads` vs serial. `simulate_pair`'s `ArgumentError`s propagate (ASVS V5). Includes are `isdefined`-guarded so the module is safe both standalone and in the harness.
- **Tests** — the SC-1 and D-11/D-12 placeholders in `test_data_pipeline.jl` are replaced with real gates: SC-1 asserts `(128,64)` shape, mask ∈ {0,1}, `mask=0 ⟺ value==0`, `all(isfinite)`, θ 7×64, AUG_DIM, and a propagated `ArgumentError`; D-11/D-12 asserts `parallel==serial` byte-identity and shuffled-order re-sort identity. Proven byte-identical under `-t 1` AND `-t 4`.
- **NOTES.md §5** — documents the 128-dim layout, `AUG_DIM`, the 14-moment list, the D-13 negative-tail caveat, the imsize table + expected cost, and the `-t auto`/`JULIA_NUM_THREADS` launch requirement (D-12).

## Task Commits

1. **Task 1: encode.jl — 128-dim D-01 vector + D-02 augmented superset** — `6928c1f` (feat)
2. **Task 2: seeding.jl — Random123 keyed RNG, disjoint salts, imsize sampler** — `d91934b` (feat)
3. **Task 3: generate.jl in-memory core + SC-1/D-11/D-12 gates** — `0b09564` (feat)

## Verification

- `julia --project=spike spike/test/test_data_pipeline.jl` → SC-1 (8 pass) + D-11/D-12 (9 pass), 17 pass / 3 broken (the remaining Wave-3 SC-2/3/4 placeholders).
- `julia --project=spike -t 4 spike/test/test_data_pipeline.jl` → identical 17 pass (thread-independence proven).
- `julia --project=spike spike/test/runtests.jl` → EXIT=0 (full suite incl. resolve-risk gate; guarded includes cause no double-include errors).
- `git diff --quiet f581d95 -- src Project.toml Manifest.toml` → CLEAN (src + root manifests byte-identical to baseline).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Guarded generator includes to prevent double-include in the harness**
- **Found during:** Task 3
- **Issue:** `runtests.jl` includes `test_simulator.jl` (which loads `contract.jl`/`forward.jl`) and then `test_data_pipeline.jl` (which includes `generate.jl`, which the plan specifies should include `contract.jl`/`prior.jl`/`forward.jl`). Unconditional re-includes would re-evaluate the `src/` struct definitions (e.g. `MultiChannelImage`) in the same `Main` module.
- **Fix:** wrapped the three Phase-2-chain includes in `isdefined(@__MODULE__, :symbol) || include(...)` guards; `seeding.jl`/`encode.jl` (functions/consts only, safe to redefine) stay unguarded. The module now loads cleanly both standalone and inside the harness.
- **Files modified:** spike/data/generate.jl
- **Commit:** 0b09564

**2. [Rule 1 - Robustness] isequal (not ==) for the parallel==serial byte-identity asserts**
- **Found during:** Task 3
- **Issue:** `==` returns `false` for `NaN`; a degenerate (fully-missing) augmented moment is documented NaN-safe, so a `==` byte-identity assert could spuriously fail even when both deterministic runs produced the identical NaN.
- **Fix:** used `isequal(...)` (treats `NaN===NaN`) for the θ/summary equality assertions; integer index/imsize vectors keep `==`.
- **Files modified:** spike/test/test_data_pipeline.jl
- **Commit:** 0b09564

## Known Stubs

None. SC-2 (cache round-trip), SC-3 (leak-free standardization), and SC-4 (k-fold + holdout) remain `@test_skip true` placeholders — these are Wave-3/Wave-4 scope by design (03-01 scaffold), not stubs introduced by this plan.

## Self-Check: PASSED

- FOUND: spike/data/encode.jl
- FOUND: spike/data/seeding.jl
- FOUND: spike/data/generate.jl
- FOUND: .planning/phases/03-training-data-pipeline/03-02-SUMMARY.md
- FOUND commit: 6928c1f (Task 1)
- FOUND commit: d91934b (Task 2)
- FOUND commit: 0b09564 (Task 3)

---
*Phase: 03-training-data-pipeline*
*Completed: 2026-06-28*
