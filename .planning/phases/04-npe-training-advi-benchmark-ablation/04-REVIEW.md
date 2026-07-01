---
phase: 04-npe-training-advi-benchmark-ablation
reviewed: 2026-07-01T00:00:00Z
depth: standard
files_reviewed: 13
files_reviewed_list:
  - spike/baseline/model.jl
  - spike/baseline/run_advi.jl
  - spike/data/loader.jl
  - spike/npe/ablation.jl
  - spike/npe/architecture.jl
  - spike/npe/benchmark.jl
  - spike/npe/infer.jl
  - spike/npe/resimulate.jl
  - spike/npe/run_thread_sweep.jl
  - spike/npe/scaling.jl
  - spike/npe/train_npe.jl
  - spike/test/runtests.jl
  - spike/test/test_npe.jl
findings:
  critical: 0
  warning: 5
  info: 7
  total: 12
status: issues_found
---

# Phase 4: Code Review Report

**Reviewed:** 2026-07-01
**Depth:** standard
**Files Reviewed:** 13
**Status:** issues_found

## Summary

Reviewed the Phase-4 NPE training / ADVI-baseline / benchmark / ablation code at standard
depth, tracing the cross-file chains (`resimulate_holdout → advi_pair → build_coloc_model →
coloc_model`, and the NPE `standardize_summary → posterior_for → reconstruct → ρ` path).

**Decoupling constraint: respected.** No file edits `src/` or the root manifests; `src/`
is reached only through read-only `include()`s in `model.jl`, and the spawned-worker
subprocess uses Julia `Cmd` backticks (argv, no shell) with only fixed paths — no command
injection, no root-package coupling. The JLD2 atomic-write idiom (`.tmp` → reopen-assert →
`mv(...; force=true)`) is applied consistently.

No BLOCKER-tier defects were found in the committed happy path (`:min` summary, `:min`
trained model), and the test suite exercises it. However, several correctness-adjacent and
benchmark-honesty issues exist. The most important is a **latent variant mismatch**: the
benchmark and scaling code carry `m.variant` through their signatures but hardcode the
128-dim `:min` summary everywhere they actually read data — so if the phase's own ablation
(`choose_summary`) selects `:aug`, the entire NPE-vs-ADVI comparison breaks. The scaling
"amortization exponents" are also fit on analytic closed-form curves rather than measured
timings, making the headline amortization claim partly tautological.

No structural findings block was provided; the `## Structural Findings (fallow)` section is
therefore omitted.

## Narrative Findings (AI reviewer)

## Warnings

### WR-01: Benchmark/scaling hardcode the `:min` summary while threading `m.variant` — an `:aug` model (a valid ablation outcome) breaks the comparison

**File:** `spike/npe/benchmark.jl:160`, `spike/npe/benchmark.jl:187-188`, `spike/npe/benchmark.jl:246-249`, `spike/npe/scaling.jl:283`
**Issue:** `choose_summary` (`ablation.jl`) can legitimately return `:aug`, and the trained
model carries `variant`/`d_in` (142 for `:aug`). But every place the benchmark and scaling
code actually reads a summary, it uses the 128-dim `:min` matrix:
- `rmse_report` reads `ho.summary_min[:, stack_cols]` and passes `m.variant` to
  `standardize_summary` (`benchmark.jl:160`, `:187-188`).
- `_time_npe_pair` / GPU branch build the summary via `encode_d01(patch_summary(...))`,
  which yields the 128-dim `:min` vector (`benchmark.jl:246-249`, `:324-325`).
- `scaling_over_imsize` does the same at `scaling.jl:283`.

For an `:aug` model, `standardize_summary(S_min, zt_aug, :aug)` calls
`_row_partition(:aug, 128)` → `cont = vcat(1:64, 129:128) = 1:64`, then
`transform(zt, S[1:64, :])` where `zt` was fit on the 78 continuous `:aug` rows → a
`DimensionMismatch`. Even if it did not throw, feeding a 128-dim vector to a `Dense(142→width)`
summary net is a hard error. The design pretends to be variant-generic but is only correct
for `:min`. This is dormant in the committed config (`:min`) but is a real defect in a
configuration the phase itself is built to produce.
**Fix:** Select the summary matrix by variant instead of hardcoding `summary_min`, and make
the NPE clock extract the variant's summary. For example:

```julia
# in rmse_report / speedup_report / scaling
summ = m.variant === :aug ? ho.summary_aug : ho.summary_min
Z = standardize_summary(summ[:, stack_cols], m.zt, m.variant)

# and the on-the-fly path needs the aug encoder, e.g.:
raw = m.variant === :aug ? encode_aug(patch_summary(mci)) : encode_d01(patch_summary(mci))
```

Or, if `:aug` is explicitly out of scope for the benchmark, assert `m.variant === :min` at
the top of `rmse_report`/`speedup_report`/`scaling_over_*` so the failure is loud and
immediate rather than a downstream dimension crash.

### WR-02: ADVI warm-up may not cover the real-data parameter shape, so the first pair's timed `wall_clock` can include JIT compilation

**File:** `spike/baseline/run_advi.jl:207-229`
**Issue:** The warm-up fits `vi()` on `smoke_mci_pair()` (fixed `imsize=(128,128)`) to move
JIT compilation outside the timed loop. But `coloc_model`'s parameter dimension is
`filldist(..., num_control)` where `num_control = size(control, 1)` = the number of
per-patch correlation values `_prepare_data` returns, which depends on how many patches
survive the zero/background exclusion. The re-simulated holdout stacks use
`sample_imsize(rng)` (variable size) and can yield a different surviving-patch count than the
smoke pair. If it differs, ForwardDiff re-specializes (chunk size / array length) on the
first real pair, so `wall_clock[1]` still absorbs compile time — inflating `t_advi` and the
`t_advi / t_npe` speedup for that pair. The median headline is somewhat protected, but the
per-pair `speedups` and `t_advi` reported downstream are contaminated.
**Fix:** Warm up on a stack drawn from the SAME `resimulate_holdout` path (e.g.
`resimulate_holdout(dir, 1)` / `(dir, 2)`) rather than the synthetic `smoke_mci_pair`, so the
warm-up model has the exact parameter dimension the timed pairs use. Alternatively, run one
untimed `vi()` on pair 1's model before the timed loop.

### WR-03: Asymmetric to-posterior clock — summary extraction is timed for NPE but excluded for ADVI

**File:** `spike/npe/benchmark.jl:245-250` vs `spike/baseline/run_advi.jl:131-137`
**Issue:** `_time_npe_pair` puts `encode_d01(patch_summary($mci))` INSIDE the `@belapsed`
block, so the NPE clock is charged for summary extraction. `advi_pair` builds the model
(`build_coloc_model`, which runs `_prepare_data` → patch/correlation) BEFORE `t0`, so the
ADVI `wall_clock` EXCLUDES the identical summary-extraction cost. The header of `benchmark.jl`
and `scaling.jl` repeatedly assert a "symmetric to-posterior clock", but the shared
summary-extraction term is timed on only one side. The direction is conservative for the NPE
(it makes the NPE look slower, so the >100× is if anything understated), so this is not a
fraud risk — but the "symmetric" framing is inaccurate and the per-pair numbers are not
apples-to-apples.
**Fix:** Either (a) move summary extraction OUTSIDE the NPE timed block to match ADVI
(measuring forward pass only on both sides), or (b) keep it inside and update the header
comments to state the clock is deliberately conservative on the NPE side (summary extraction
included for NPE, excluded for ADVI) rather than "symmetric".

### WR-04: Scaling "empirical exponents" are fit on analytic curves, so the amortization result is tautological

**File:** `spike/npe/scaling.jl:208-212`, and its test `spike/test/test_npe.jl:336-337`
**Issue:** `npe_time_N = C_train .+ Ngrid .* t_fwd_per_dataset` and
`advi_time_N = Ngrid .* t_advi_per_dataset` are deterministic affine functions of `N`;
`_loglog_slope` is then applied to these closed-form values. `npe_exponent_N < advi_exponent_N`
is therefore true *by construction* for any `C_train > 0` — nothing about the network's
actual per-dataset scaling is measured. The test `@test rN.npe_exponent_N < rN.advi_exponent_N`
(with `train_cost = 5.0` supplied) cannot fail regardless of code correctness. The header
honestly labels this a "MODEL", but presenting a fit slope of an analytic line as an
"empirical scaling exponent" overstates what was measured.
**Fix:** Either measure `npe_time_N(N)` by actually running N forward passes end-to-end and
fitting the slope of the measured wall-clock, or rename `*_exponent_*` /
`npe_time_N`/`advi_time_N` to make explicit they are analytic amortization *models* and drop
the tautological `<` assertion in favor of a check on the measured inputs (`t_fwd_per_dataset`,
`t_advi_per_dataset`, `C_train`) themselves.

### WR-05: NPE reproducibility depends on global-RNG seeding and call order, contrary to the Random123 reproducibility mandate

**File:** `spike/npe/benchmark.jl:134`, `spike/test/test_npe.jl:139-147`
**Issue:** `rmse_report` establishes reproducible posterior draws with
`Random.seed!(UInt32(master_seed % 0xFFFFFFFF))`, i.e. by seeding the GLOBAL default RNG that
`sampleposterior` consumes. Determinism then depends on the exact order/quantity of RNG
consumption across `assess`, the interval loop, and the Δρ loop. The SC1 test's
`dρ ≈ mean(ρs .- ρc)` check (`test_npe.jl:139-147`) relies on `delta_rho` and the manual
recomputation drawing from the global stream in identical order. CLAUDE.md makes Random123
the reproducibility tool precisely because global-RNG seeding is fragile: any future
interleaved RNG consumer (a logging call, a reordered loop, a threaded sampler) silently
changes the results without changing the seed.
**Fix:** Thread an explicit `rng::AbstractRNG` (e.g. a Random123 counter-based stream keyed by
`master_seed`) into `posterior_for`/`sampleposterior` and pass it through
`rho_hat`/`rho_draws`/`delta_rho`, rather than mutating the process-global RNG. This also
removes the hidden coupling in the SC1 equality test.

## Info

### IN-01: Odd `% 0xFFFFFFFF` narrowing for the seed

**File:** `spike/npe/benchmark.jl:134`
**Issue:** `UInt32(master_seed % 0xFFFFFFFF)` reduces modulo 4294967295 (not 2^32). It maps
`0xFFFFFFFF → 0` and any multiple of 4294967295 to 0, which is a needless surprise; the
intent is clearly "take the low 32 bits".
**Fix:** Use a mask: `UInt32(master_seed & 0xFFFFFFFF)` (or `mod(master_seed, UInt32)`).

### IN-02: Duplicated RMSE/scale helpers across benchmark and ablation

**File:** `spike/npe/benchmark.jl:78-91` and `spike/npe/ablation.jl:86-100`
**Issue:** `_rmse_vector`/`_abl_rmse_vector` and `_theta_scale`/`_abl_theta_scale` are
byte-for-byte equivalent. The duplication is acknowledged in comments (to avoid pulling the
heavier benchmark chain into ablation), but it is a drift hazard: a fix to one will not track
the other.
**Fix:** Extract the two helpers into a tiny shared file (e.g. `spike/npe/_rmse_utils.jl`)
that both `include` via the existing guarded pattern, so there is a single definition.

### IN-03: `_rmse_vector`/`_abl_rmse_vector` hardcode the 7-parameter axis

**File:** `spike/npe/benchmark.jl:82`, `spike/npe/ablation.jl:90`
**Issue:** Both build `[lut["θ$i"] for i in 1:7]`, a magic `7` decoupled from `NPE_D`
(`architecture.jl:52`). If `NPE_D` ever changes, these silently mis-extract (or `KeyError`)
rather than tracking the true θ dimension.
**Fix:** Iterate `1:NPE_D` (or `1:length(values)`) instead of the literal `1:7`.

### IN-04: Tautological tests re-implement the unit under test

**File:** `spike/test/test_npe.jl:290-293` (SC5b), `spike/test/test_npe.jl:269-270` (SC4d)
**Issue:** SC5b re-derives the exact `choose_summary` boolean logic and asserts equality with
`choose_summary`'s output; SC4d re-derives `fold_wins_aug` the same way `ablate` computes it.
These pass by construction and give no independent oracle — a sign error copied into both the
function and the test would not be caught.
**Fix:** Assert against fixed expected outcomes for hand-constructed fixture inputs (e.g. a
result where aug clearly wins vs. clearly ties), not against a re-implementation of the
function's own formula.

### IN-05: `rmse_ratio_ok` is computed per thread count but never aggregated

**File:** `spike/npe/run_thread_sweep.jl:100-108`, `spike/npe/run_thread_sweep.jl:123-156`
**Issue:** `thread_sweep_measure` records `rmse_ratio_ok` into each per-N file, but
`aggregate_thread_sweep` reads `median_speedup`/`t_*_median` only and drops `rmse_ratio_ok`.
A thread count whose paired RMSE violated the D-09 tolerance would not surface anywhere in the
aggregated table or its persisted artifact.
**Fix:** Read and include `rmse_ratio_ok` in the aggregation NamedTuple and the saved file, so
the RMSE-validity of each swept thread count is visible alongside its speedup.

### IN-06: Untyped `Any[...]` layer vector in the architecture builder

**File:** `spike/npe/architecture.jl:80-85`
**Issue:** `layers = Any[Dense(...)]` is an abstractly-typed vector. It is harmless here
because the layers are splatted into `Chain(layers...)` (which recovers concrete types), but
an untyped container in the one arch-defining function reads as accidental rather than
deliberate. CONVENTION — the codebase otherwise favors concrete typing at construction.
**Fix:** Build a concretely-typed vector, e.g. `layers = Dense[Dense(d_in, width, gelu)]`, or
construct the `Chain` via a comprehension/`Tuple`.

### IN-07: Non-UTC, inconsistent artifact timestamp

**File:** `spike/npe/train_npe.jl:168`
**Issue:** `Dates_now()` uses `Base.Libc.strftime("%Y-%m-%dT%H:%M:%S", time())`, which formats
in LOCAL time with no timezone marker, whereas the ADVI artifact meta uses `string(Dates.now())`
(`run_advi.jl:267`). Provenance timestamps across the two cross-env artifacts are thus in
different, unlabeled time bases, complicating any later ordering/repro audit. CONVENTION.
**Fix:** Use a single UTC source for both (e.g. `string(Dates.now(Dates.UTC))`), accepting the
`Dates` dependency the header tried to avoid, or append an explicit `Z`/offset.

---

_Reviewed: 2026-07-01_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
