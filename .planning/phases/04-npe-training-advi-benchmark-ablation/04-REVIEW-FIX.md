---
phase: 04-npe-training-advi-benchmark-ablation
fixed_at: 2026-07-02T00:00:00Z
review_path: .planning/phases/04-npe-training-advi-benchmark-ablation/04-REVIEW.md
iteration: 2
findings_in_scope: 12
fixed: 11
skipped: 1
status: partial
---

# Phase 4: Code Review Fix Report

**Fixed at:** 2026-07-02
**Source review:** .planning/phases/04-npe-training-advi-benchmark-ablation/04-REVIEW.md
**Iteration:** 2 (`--all` scope; supersedes the iteration-1 `critical_warning` report)

**Summary (all 12 findings):**
- Findings in scope: 12 (5 warnings + 7 info)
- Fixed: 11 (4 warnings in iteration 1 + 7 info in iteration 2)
- Skipped: 1 (WR-05)

This iteration-2 report gives the FULL picture across both passes. The 5 warnings were
handled in the iteration-1 `critical_warning` pass (WR-01..WR-04 fixed, WR-05 skipped) and
are recorded below unchanged. Iteration 2 (`--all`) fixed the 7 remaining Info findings.

All edits stayed inside the decoupled `spike/` tree (CLAUDE.md hard constraint respected --
no `src/`, root `Project.toml`/`Manifest.toml`, or manuscript-pipeline edits). The phase
test gate `julia --project=spike spike/test/runtests.jl` exits 0 after every commit; the
Phase-4 testset moved from 142/142 (post-warning pass) to 150/150 (the IN-04 fixture-oracle
assertions and the IN-05 `rmse_ratio_ok` assertions add net coverage; no tautologies added).

## Fixed Issues

### WR-01: Benchmark/scaling hardcode the `:min` summary while threading `m.variant`

**Files modified:** `spike/npe/benchmark.jl`, `spike/npe/scaling.jl`
**Commit:** 2064988 (iteration 1)
**Applied fix:** Added `_holdout_summary(ho, variant)` and `_encode_variant(mci, variant)`
helpers and routed every data-reading site (rmse_report stacked `Z` + Δρ loop, the NPE
timing clocks, scaling_over_imsize) through them, so the summary matrix/encoding matches the
trained model's variant. `:min` behaviour is byte-identical; an `:aug` model no longer hits a
`DimensionMismatch`/`142→width` feed error.

### WR-02: ADVI warm-up may not cover the real-data parameter shape

**Files modified:** `spike/baseline/run_advi.jl`
**Commit:** 1a3d201 (iteration 1)
**Applied fix:** Replaced the fixed-size `smoke_mci_pair()` warm-up with a warm-up built from
`resimulate_holdout(dir, 1)`/`(dir, 2)` -- the same re-sim path the timed loop consumes -- so
the warm-up `coloc_model` carries the real per-patch parameter dimension and no ForwardDiff
re-compile leaks into `wall_clock[1]`.

### WR-03: Asymmetric to-posterior clock (summary extraction timed for NPE, not ADVI)

**Files modified:** `spike/npe/benchmark.jl`, `spike/npe/run_thread_sweep.jl`
**Commit:** a403be3 (iteration 1)
**Applied fix:** Documentation-accuracy fix (review option (b)). The NPE clock deliberately
includes shared summary extraction while ADVI excludes it; the asymmetry only shrinks the
reported speedup (a conservative lower bound). Rewrote the `_time_npe_pair`/`speedup_report`
docstrings and the `run_thread_sweep.jl` header to state the clock is "NPE-conservative"
rather than "symmetric". No behavioural change.

### WR-04: Scaling "empirical exponents" are fit on analytic curves (tautological)

**Files modified:** `spike/npe/scaling.jl`, `spike/test/test_npe.jl`
**Commit:** 02a6b27 (iteration 1)
**Applied fix:** Labelled the exponents as slopes of an ANALYTIC amortization MODEL (not
measured empirical exponents) and replaced the two tautological SC3 assertions with checks on
the MEASURED inputs (`train_cost > 0`, `t_fwd_per_dataset > 0`, `t_advi_per_dataset > 0`,
`t_fwd_per_dataset < t_advi_per_dataset`, finite positive `crossover_N`) that actually make
amortization real.

### IN-01: Odd `% 0xFFFFFFFF` narrowing for the seed

**Files modified:** `spike/npe/benchmark.jl`
**Commit:** 47a4d67
**Applied fix:** `rmse_report` seeded the global RNG with `UInt32(master_seed % 0xFFFFFFFF)`,
which reduces modulo 4294967295 (not 2^32) and maps `0xFFFFFFFF` (and any multiple of it) to
0. Changed to a bit mask -- `UInt32(master_seed & 0xFFFFFFFF)` -- so "take the low 32 bits" is
exact and surprise-free. The mask guarantees the value fits `UInt32` for any signed/unsigned
`master_seed`.

### IN-02: Duplicated RMSE/scale helpers across benchmark and ablation

**Files modified:** `spike/npe/_rmse_utils.jl` (new), `spike/npe/benchmark.jl`, `spike/npe/ablation.jl`
**Commit:** 7fcccd2
**Applied fix:** `benchmark.jl` (`_rmse_vector`/`_theta_scale`) and `ablation.jl`
(`_abl_rmse_vector`/`_abl_theta_scale`) held byte-identical copies -- a drift hazard. Extracted
both helpers into a new tiny shared `spike/npe/_rmse_utils.jl` that both files pull in via the
existing guarded include pattern (`isdefined(...) || include(...)`), and switched ablation's
two call sites to the shared names. Behaviour unchanged. (New file documented here per the
"do not create files unless the fix requires it" rule -- the review explicitly prescribed a
shared file.)

### IN-03: `_rmse_vector` hardcodes the 7-parameter axis

**Files modified:** `spike/npe/_rmse_utils.jl`
**Commit:** cfaa15c
**Applied fix:** The shared `_rmse_vector` built `[lut["θ$i"] for i in 1:7]`, a magic 7
decoupled from `NPE_D` (architecture.jl). Changed to iterate `1:NPE_D` so the per-parameter
RMSE axis tracks the true θ dimension. `NPE_D` is resolved at call time and the architecture
chain is always loaded before any rmse call, so include order is unaffected.

### IN-04: Tautological tests re-implement the unit under test

**Files modified:** `spike/test/test_npe.jl`
**Commit:** 332dbf1
**Applied fix:** SC5b re-derived `choose_summary`'s exact rule and SC4d re-derived `ablate`'s
`fold_wins_aug` count against the real fixture result -- both passed by construction with no
independent oracle. Replaced with hand-built inputs and fixed known outcomes: SC5b drives
`choose_summary` with four synthetic result tuples (clear win, sub-margin, inconsistent, and
the exact margin/consistency boundary) expecting `:aug`/`:min`; SC4d verifies the aug-win
counting predicate on three fixtures with known win counts (3/0/5), then cross-checks that
`ablate`'s returned `fold_wins_aug` agrees on its own `rho_per_fold`. Verified by running the
full suite -- the new assertions genuinely pass and test the decision/counting logic (not new
tautologies).

### IN-05: `rmse_ratio_ok` computed per thread count but never aggregated

**Files modified:** `spike/npe/run_thread_sweep.jl`, `spike/test/test_npe.jl`
**Commit:** a60497c
**Applied fix:** `thread_sweep_measure` wrote the per-thread-count D-09 RMSE-validity flag into
each per-N file, but `aggregate_thread_sweep` dropped it. Read `rmse_ratio_ok` from each row
and included it in both the returned NamedTuple and the persisted JLD2 (plus docstring return
list), and asserted it in the SC3 thread-sweep test (returned value + saved artifact).

### IN-06: Untyped `Any[...]` layer vector in the architecture builder

**Files modified:** `spike/npe/architecture.jl`
**Commit:** f538f1a
**Applied fix:** `build_estimator` accumulated layers in `Any[Dense(...)]`. Changed to a
concretely-typed `Dense[Dense(...)]` vector, matching the codebase's construct-time
concrete-typing convention. Harmless before (Chain splat recovers types); now deliberate.

### IN-07: Non-UTC, inconsistent artifact timestamp

**Files modified:** `spike/npe/train_npe.jl`, `spike/baseline/run_advi.jl`
**Commit:** 1af27c6
**Applied fix:** The NPE artifact used `Base.Libc.strftime(..., time())` (LOCAL, unlabeled),
while the ADVI artifact used `string(Dates.now())` (also local, unlabeled) -- two provenance
timestamps in different unlabeled bases. Switched BOTH to a UTC source with an explicit `Z`
suffix (`string(Dates.now(Dates.UTC)) * "Z"`) so cross-env ordering/repro audits are
unambiguous; added `using Dates` (a stdlib) to train_npe.jl. Existing committed artifacts keep
their old-format timestamp until regenerated; no test asserts the timestamp format.

## Skipped Issues

### WR-05: NPE reproducibility depends on global-RNG seeding

**File:** `spike/npe/benchmark.jl:134`, `spike/test/test_npe.jl:139-147`
**Reason:** skipped -- the prescribed fix (thread an explicit `rng::AbstractRNG` into
`posterior_for`/`sampleposterior`) is not implementable without patching a pinned dependency
outside `spike/`. The pinned NeuralEstimators **v0.2.1** `sampleposterior` has no `rng`
parameter anywhere in its call chain (the flow's base draw is a bare `randn` against the
global default RNG), so a threaded stream would be silently ignored; genuine per-call
determinism would require editing NeuralEstimators' `NormalisingFlow.jl` in the package depot,
violating the decoupling constraint, and would desync the SC1 equality test that relies on
identical global seeding. The existing `Random.seed!` approach is the only reproducibility
lever the pinned library exposes. Recommend deferring to productionization, where a
backend/version accepting an `rng` (or a vendored flow) can be adopted deliberately. (Note:
the related IN-01 `% 0xFFFFFFFF` narrowing on the same line WAS fixed this iteration.)

---

_Fixed: 2026-07-02_
_Fixer: Claude (gsd-code-fixer)_
_Iteration: 2_
