---
phase: 04-npe-training-advi-benchmark-ablation
fixed_at: 2026-07-02T00:00:00Z
review_path: .planning/phases/04-npe-training-advi-benchmark-ablation/04-REVIEW.md
iteration: 1
findings_in_scope: 5
fixed: 4
skipped: 1
status: partial
---

# Phase 4: Code Review Fix Report

**Fixed at:** 2026-07-02
**Source review:** .planning/phases/04-npe-training-advi-benchmark-ablation/04-REVIEW.md
**Iteration:** 1

**Summary:**
- Findings in scope: 5 (WR-01..WR-05; the 7 IN-* Info items are out of scope for `critical_warning`)
- Fixed: 4
- Skipped: 1

All fixes were applied inside the decoupled `spike/` tree only (CLAUDE.md hard
constraint respected — no `src/`, root manifest, or manuscript-pipeline edits). The
full spike test gate `julia --project=spike spike/test/runtests.jl` exits 0 after the
fixes; the Phase-4 testset moved from 139/139 to 142/142 (WR-04 replaced two
tautological assertions with five measured-input assertions).

## Fixed Issues

### WR-01: Benchmark/scaling hardcode the `:min` summary while threading `m.variant`

**Files modified:** `spike/npe/benchmark.jl`, `spike/npe/scaling.jl`
**Commit:** 2064988
**Applied fix:** Added two variant-selection helpers in `benchmark.jl` —
`_holdout_summary(ho, variant)` (returns `ho.summary_aug` for `:aug`, `ho.summary_min`
for `:min`) and `_encode_variant(mci, variant)` (returns
`encode_aug(mci, patch_summary(mci))` for `:aug`, `encode_d01(patch_summary(mci))` for
`:min`, computing `patch_summary` once). Routed every data-reading site through them:
`rmse_report`'s stacked `Z` and its Δρ loop now use `_holdout_summary`; `_time_npe_pair`
and the GPU timing branch now use `_encode_variant`; `scaling_over_imsize`'s
summary-extraction timing and reference-forward-pass `Zref` now use `_encode_variant`.
`:min` behaviour is byte-identical, so the committed config and its tests are unchanged,
while an `:aug` model (a valid `choose_summary` ablation outcome) no longer hits a
`DimensionMismatch`/`142→width` feed error.

### WR-02: ADVI warm-up may not cover the real-data parameter shape

**Files modified:** `spike/baseline/run_advi.jl`
**Commit:** 1a3d201
**Applied fix:** Replaced the `smoke_mci_pair()` warm-up (fixed 128×128) in
`run_all_pairs` with a warm-up built from `resimulate_holdout(dir, 1)` /
`resimulate_holdout(dir, 2)` — the same re-simulation path the timed loop consumes — so
the warm-up `coloc_model` carries the real holdout's per-patch parameter dimension
(`num_control` = surviving-patch count from `_prepare_data`). This keeps the
ForwardDiff/AdvancedVI specialization the first timed pair reuses matched to the real
data shape, preventing a re-compile from leaking into `wall_clock[1]`. The `npairs >= 1`
assertion earlier in the function guarantees holdout entries 1 and 2 exist;
`resimulate_holdout` is already the loop's re-sim primitive, so no new dependency.

### WR-03: Asymmetric to-posterior clock — summary extraction timed for NPE but not ADVI

**Files modified:** `spike/npe/benchmark.jl`, `spike/npe/run_thread_sweep.jl`
**Commit:** a403be3
**Applied fix:** Documentation-accuracy fix (option (b) from the review). The NPE clock
deliberately includes shared summary extraction inside `@belapsed` while the ADVI
`wall_clock` excludes it (extraction runs in `build_coloc_model` before `advi_pair`'s
`t0`); the asymmetry only shrinks the reported speedup, so the >100× headline is a
conservative lower bound rather than a fraud risk. Rewrote the `_time_npe_pair` and
`speedup_report` docstrings and the `run_thread_sweep.jl` header comment to state the
clock is "NPE-conservative" (extraction charged to NPE, excluded for ADVI) instead of
claiming a "symmetric" to-posterior clock. No behavioural change; the measurement
semantics are intentionally left as-is (the header contract is that the NPE clock starts
from the raw stack, D-08), only the misleading "symmetric" wording is corrected.

### WR-04: Scaling "empirical exponents" are fit on analytic curves (tautological)

**Files modified:** `spike/npe/scaling.jl`, `spike/test/test_npe.jl`
**Commit:** 02a6b27
**Applied fix:** In `scaling_over_N`, `npe_time_N`/`advi_time_N` are closed-form
affine/linear functions of `N`, so `npe_exponent_N < advi_exponent_N` holds by
construction for any `C_train > 0` and the SC3 test `@test rN.npe_exponent_N <
rN.advi_exponent_N` (plus `@test rN.advi_exponent_N > 0.9`) could not fail. Labelled the
exponents in the `scaling_over_N` docstring and inline comments as slopes of an ANALYTIC
amortization MODEL (not measured empirical exponents), and replaced the two tautological
SC3 assertions with checks on the MEASURED inputs that actually make amortization real:
`train_cost > 0`, `t_fwd_per_dataset > 0`, `t_advi_per_dataset > 0`,
`t_fwd_per_dataset < t_advi_per_dataset` (the measured amortization premise), and
`isfinite(crossover_N) && crossover_N > 0`. These now fail if the measured per-dataset
NPE forward pass is not actually cheaper than the per-dataset ADVI cost. Test suite
confirms the new assertions pass (SC3-scaling 21/21).

## Skipped Issues

### WR-05: NPE reproducibility depends on global-RNG seeding

**File:** `spike/npe/benchmark.jl:134`, `spike/test/test_npe.jl:139-147`
**Reason:** skipped — the prescribed fix is not implementable without patching a pinned
dependency, which is out of scope (and outside `spike/`). The review's fix is to thread
an explicit `rng::AbstractRNG` (Random123 counter-based) into
`posterior_for`/`sampleposterior`. I verified the pinned NeuralEstimators **v0.2.1**
(the version the spike `Manifest.toml` resolves,
`packages/NeuralEstimators/gFxuZ`): `sampleposterior(estimator::PosteriorEstimator, Z;
N, device, use_gpu, kwargs...)` delegates to `sampleposterior(flow::NormalisingFlow, tz,
N; device)`, whose base draw is a bare `randn(Float32, d, N*K)` against the GLOBAL
default RNG — there is **no `rng` parameter anywhere in the call chain**. Threading a
Random123 stream into `posterior_for` would therefore be silently ignored by the flow;
genuine per-call determinism would require editing NeuralEstimators' `NormalisingFlow.jl`
in the package depot, which is (a) not a `spike/` edit and violates the decoupling
constraint's spirit, and (b) would desync the committed reproducible sampling path and
the SC1 `dρ ≈ mean(ρs .- ρc)` equality test that currently relies on identical global
seeding. The existing `Random.seed!` approach is the only reproducibility lever the
pinned library exposes. Recommend deferring to productionization, where a backend/version
that accepts an `rng` (or a vendored flow) can be adopted deliberately. (Note: the
related IN-01, `% 0xFFFFFFFF` seed narrowing on the same line, is an Info finding and out
of scope for this `critical_warning` pass.)

---

_Fixed: 2026-07-02_
_Fixer: Claude (gsd-code-fixer)_
_Iteration: 1_
