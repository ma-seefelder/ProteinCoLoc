---
phase: 02-forward-simulator-summary-contract
reviewed: 2026-06-26T00:00:00Z
depth: standard
files_reviewed: 6
files_reviewed_list:
  - spike/contract.jl
  - spike/simulator/forward.jl
  - spike/simulator/calibration.jl
  - spike/simulator/prior.jl
  - spike/02_simulator_demo.jl
  - spike/test/test_simulator.jl
findings:
  critical: 2
  warning: 7
  info: 3
  total: 12
status: issues_found
---

# Phase 02: Code Review Report

**Reviewed:** 2026-06-26
**Depth:** standard
**Files Reviewed:** 6
**Status:** issues_found

## Summary

Reviewed the complete Phase-2 deliverable: the 2-channel physics forward-simulator
(`simulate_pair`), its shared-latent generator, the offline induced-mu calibration
(`calibrate` / `emit_ghat`), the transform-sampling prior (`sample_prior`), the
CairoMakie plausibility demo, and the test suite.

The simulator physics, the PAVA isotonic-regression fit, the piecewise-linear ghat
inversion, the BG_FLOOR / _exclude_zero contract traps, and the rng-threading pattern
are all structurally sound. Two blockers require fixes before Phase 3 can depend on
this surface: a wrong-return-type extension of `Base.summary` that silently poisons
Julia's display machinery, and a crash path in `induced_mu` that is reachable via the
published API (`imsize = (8,8)` passes validation but produces an all-missing summary).
Seven warnings cover latent crashes, fragile one-sample tests, a shared mutable RNG
in the demo, and a reversed axis-convention in the warp call.

---

## Critical Issues

### CR-01: `Base.summary` extended with wrong return type — corrupts Julia display

**File:** `spike/contract.jl:83`

**Issue:** The function is defined as `function Base.summary(mci::MultiChannelImage)` and
returns `Matrix{Union{Float64,Missing}}`. `Base.summary(x)` is part of Julia's display
infrastructure: `show(io, x)` → `Base.summary(io, x)` → `print(io, Base.summary(x))`.
The expected return type is `String`. Returning a Matrix means:

1. Any REPL display of a `MultiChannelImage` (including in error context) executes the
   full `patch()`/`correlation()` pipeline and then prints the resulting numeric matrix
   as the "description" of the object — wrong semantics, expensive, and opaque.
2. Test.jl formats `@test` failure messages by calling `show()` on operands; a failing
   assertion involving a `MultiChannelImage` will embed a 8×8 correlation matrix in the
   error message body.
3. Any future code that calls `s = Base.summary(mci)` and treats `s` as a `String` will
   receive a `MethodError` on the first string operation.

The docstring's justification ("so the plan's `summary(mci)` call site works without
shadowing Base's generic") is specious: a module-local function named `patch_summary` or
`compute_summary` achieves the same call-site syntax without piggybacking on Base.

**Fix:**
```julia
# contract.jl — rename; do NOT extend Base.summary
"""
    patch_summary(mci::MultiChannelImage) -> Matrix{Union{Float64,Missing}}

The SIM-03 summary statistic: the fixed 8×8 per-patch Pearson correlation matrix ...
"""
function patch_summary(mci::MultiChannelImage)
    x = mci.data[1]; y = mci.data[2]
    xp, yp = patch.([x, y], 8)
    return correlation(xp, yp; method = :pearson)
end

# Update induced_mu and all three call sites in test_simulator.jl:
induced_mu(mci::MultiChannelImage) = Statistics.mean(skipmissing(patch_summary(mci)))
```

---

### CR-02: `induced_mu` crashes on all-missing summary — reachable via published API

**File:** `spike/contract.jl:98` (root cause: `spike/simulator/forward.jl:109`)

**Issue:** `induced_mu(mci) = Statistics.mean(skipmissing(summary(mci)))` throws
`ArgumentError: mean of empty collection` when every patch in the 8×8 grid is `missing`.
This is a reachable crash through the published API:

- `simulate_pair` validates `imsize[1] ≥ 8 && imsize[2] ≥ 8` (line 109) and accepts
  `imsize = (8, 8)`.
- At `(8, 8)` with an 8×8 patch grid, each patch is 1×1 = 1 pixel.
- `src/colocalization.jl:235` sets `ρ[i,j] = missing` when `length(a) ≤ 15`.
- All 64 patches contain 1 pixel, all fail the floor: `summary(mci)` returns a
  fully-missing 8×8 matrix.
- `Statistics.mean(skipmissing(fully_missing_matrix))` crashes.

The `@test_throws ArgumentError simulate_pair(... imsize=(4,4))` test (test_simulator.jl:115)
does NOT cover `imsize=(8,8)`, leaving this path untested and open.

**Fix:** Two complementary changes:

```julia
# forward.jl — tighten the minimum dimension so each patch has > 15 pixels:
# With an 8×8 grid, minimum patch size for ≥ 16 pixels is 4 px/axis → imsize ≥ 32.
# Use 32 as the practical floor (each patch then has 4×4 = 16 px, barely over the floor).
# Better: use 64 so each patch is 8×8 = 64 px with comfortable headroom.
(imsize[1] ≥ 64 && imsize[2] ≥ 64) ||
    throw(ArgumentError("imsize dims must be ≥ 64 for the 8×8 patch grid " *
                        "to produce non-missing patches (≥15 px each), got $imsize"))

# contract.jl — make induced_mu robust regardless:
function induced_mu(mci::MultiChannelImage)
    vals = collect(skipmissing(patch_summary(mci)))
    isempty(vals) && return NaN
    return Statistics.mean(vals)
end
```

---

## Warnings

### WR-01: `calibrate()` calls `mean`/`std` on potentially empty `μs`

**File:** `spike/simulator/calibration.jl:179-180`

**Issue:** `_induced_mu_samples` filters out non-finite results (`isfinite(m) && push!(μs, m)`).
If every one of the `n_per=100` simulations at a given ρ grid point returns a non-finite
`induced_mu` (possible under degenerate parameterisations during development), then `μs`
is empty and `mean(Float64[])` / `std(Float64[])` throw `ArgumentError`. The calibration
itself shows the correct guard pattern in `_w1`: `n < 5 && return NaN`.

**Fix:**
```julia
μs = _induced_mu_samples(rng, ρ, imsize, n_per)
if isempty(μs)
    @warn "all induced_mu non-finite at ρ=$ρ; skipping grid point"
    mu_mean[k] = NaN; mu_sd[k] = NaN; continue
end
mu_mean[k] = mean(μs); mu_sd[k] = std(μs)
```

---

### WR-02: D-15 monotone test is fragile — one simulation per ρ point

**File:** `spike/test/test_simulator.jl:151-161`

**Issue:** The D-15 test sweeps only three ρ values (`[-0.7, 0.0, 0.7]`) with a single
simulation per point at a fixed seed (Xoshiro(2026)). `@test μ_ind[1] < 0.0` and
`@test issorted(μ_ind)` both rely on a single stochastic draw landing on the "correct"
side of zero/ordering. A change in upstream `randn` dispatch, image size, or any
`include`-order side effect that shifts the rng state will silently invert the ordering,
converting a real regression into a false pass and vice versa. The demo sweep and
SIM-04(a) both correctly average N=3 replicates per point.

**Fix:**
```julia
# Average 3 replicates per ρ to suppress single-draw noise:
μ_ind = map(ρ_grid) do ρ
    vals = [induced_mu(build_mci(simulate_pair(Random.Xoshiro(2026 + i), _θ(ρ);
                                               imsize = (256, 256)))) for i in 0:2]
    mean(filter(isfinite, vals))
end
```

---

### WR-03: Shared mutable `DEMO_RNG` across all demo sweeps

**File:** `spike/02_simulator_demo.jl:56, 78, 83, 88, 93`

**Issue:** `const DEMO_RNG = Random.Xoshiro(2026)` is mutated sequentially by the
ρ-sweep (line 78), spillover-sweep (line 83), shift-sweep (line 88), and induced-μ
sample (line 93). Adding or removing any point from an earlier sweep changes all
subsequent results. Each sweep should be independent for the demo to be robust to
modification. Contrast with the test file, which correctly uses per-comparison fresh
seeds.

**Fix:** Give each sweep its own seed:
```julia
# Each panel's RNG is independent and reproducible:
ρ_rng     = Random.Xoshiro(2026)
spill_rng = Random.Xoshiro(2027)
shift_rng = Random.Xoshiro(2028)
mu_rng    = Random.Xoshiro(2029)

mean_corr = [mean_patch_corr(ρ_rng, _θ_med(ρ)) for ρ in ρ_grid]
corr_spill = [mean_patch_corr(spill_rng, merge(_θ_med(ρ_FIX), (spillover = s,))) ...]
# etc.
```

---

### WR-04: SIM-02 test N=120 — W1 estimator noise nearly equals tolerance

**File:** `spike/test/test_simulator.jl:196-211`

**Issue:** The SIM-02 gate draws N=120 samples at `imsize=(512,512)` and tests
`w1 < SIM02_W1_TOL = 0.10`. The Monte-Carlo standard error of the empirical W1
estimator for a heavy-tailed Cauchy-derived target scales as O(1/√n) ≈ 0.09 for
n=120. A true null W1 of 0.04 can plausibly appear as ≥0.10 by sampling variation,
giving a false failure. The calibration itself uses `metric_N=300` (SE ≈ 0.058); the
test should match or the tolerance should be adjusted for its smaller N.

**Fix:** Either increase N to 300 in the test (accepting the extra wall-clock cost) or
widen the tolerance proportionally:
```julia
N   = 300                  # match calibration metric_N
# OR, if 120 is the compute budget:
const SIM02_W1_TOL_TEST = 0.15   # = calibration tol + 2 SE buffer at N=120
@test w1 < SIM02_W1_TOL_TEST
```

---

### WR-05: `_w1_test` has no guard for empty input

**File:** `spike/test/test_simulator.jl:170-173`

**Issue:** The test-local `_w1_test` omits the `n < 5 && return NaN` guard present in
calibration.jl's `_w1`. If `induced` is empty (all 120 `induced_mu` calls return
non-finite), `n = 0`, `qs = Float64[]`, `quantile(a, []) = []`,
`mean(Float64[]) = ArgumentError`. While effectively unreachable with BG_FLOOR > 0,
the inconsistency between the two copies is a latent bug.

**Fix:**
```julia
_w1_test(a::Vector{Float64}, b::Vector{Float64}) = begin
    n = min(length(a), length(b))
    n < 5 && return NaN                     # match calibration.jl guard
    qs = ((1:n) .- 0.5) ./ n
    mean(abs.(quantile(a, qs) .- quantile(b, qs)))
end
```

---

### WR-06: `Translation(shift_dx, shift_dy)` axis convention is reversed

**File:** `spike/simulator/forward.jl:150`

**Issue:** `CoordinateTransformations.Translation(a, b)` shifts the **first array axis**
by `a` (rows = vertical = y in image coordinates) and the **second** by `b` (columns =
horizontal = x). The call `Translation(θ.shift_dx, θ.shift_dy)` therefore applies
`shift_dx` as the vertical (row) shift and `shift_dy` as the horizontal (column) shift —
the opposite of the conventional dx=horizontal, dy=vertical naming. Because both priors
are the identical `Uniform(-1, 1)` distribution the statistics are unaffected, but the
named fields in the θ NamedTuple, the docstring description, and the Turing prior they
mirror will all describe the wrong physical axis.

**Fix:**
```julia
# Swap to match convention: dx = column (horizontal), dy = row (vertical)
shifted = warp(ch2, Translation(θ.shift_dy, θ.shift_dx), axes(ch2);
               method = BSpline(Linear()), fillvalue = BG_FLOOR)
```

---

### WR-07: Nuisance parameter validation only checks `isfinite` — range violations pass silently

**File:** `spike/simulator/forward.jl:111-113`

**Issue:** The `ArgumentError` guard (lines 111–113) rejects non-finite θ fields but
accepts physically invalid nuisances: `spillover = -0.5` is finite and passes, then
`ch1 += -0.5 * ch2` subtracts signal from ch1 (negative bleed-through), corrupting
stage 4 silently. `label_efficiency = 1.5` passes, is clamped to 1.0 at line 131 with
no warning, so the caller believes they have a 150% label efficiency when the simulation
runs at 100%. `autofluorescence = -0.5` reduces background below BG_FLOOR potentially.
The documented contract states "untrusted parameter vector (T-02-IV)" — that contract
should extend to range bounds.

**Fix:** Add range checks in the validation block:
```julia
(0.0 ≤ θ.spillover ≤ 1.0) ||
    throw(ArgumentError("spillover must be in [0,1], got $(θ.spillover)"))
(0.0 ≤ θ.autofluorescence) ||
    throw(ArgumentError("autofluorescence must be ≥ 0, got $(θ.autofluorescence)"))
(0.0 ≤ θ.label_efficiency ≤ 1.0) ||
    throw(ArgumentError("label_efficiency must be in [0,1], got $(θ.label_efficiency)"))
(0.0 ≤ θ.noise) ||
    throw(ArgumentError("noise must be ≥ 0, got $(θ.noise)"))
```

---

## Info

### IN-01: Duplicate "(b)" section label in demo comments

**File:** `spike/02_simulator_demo.jl:82, 85`

**Issue:** Both the spillover-effect sweep (line 82) and the shift-effect sweep (line 85)
are labeled `# --- (b) ...`. The axis title at line 131 uses `"(b) |shift| ↑ → corr ↓"`,
which conflicts with both having the same letter. One should be `(c)` or `(b-ii)` to
match standard figure-panel convention.

**Fix:** Change line 85 to `# --- (c) sub-pixel |shift|-effect sweep`.

---

### IN-02: `ghat.jl` generated-artifact prerequisite not documented in test header

**File:** `spike/test/test_simulator.jl:167`

**Issue:** The `include(... "simulator/prior.jl")` at line 167 transitively includes
`spike/simulator/ghat.jl`. If calibration has not been run yet, the failure is a
bare `SystemError: opening file ... ghat.jl: No such file or directory` with no
guidance. A first-time contributor who runs the test suite before calibration will
spend time diagnosing this.

**Fix:** Add a one-line comment:
```julia
# Prior requires the frozen ghat.jl artifact — run `julia --project=spike
# spike/simulator/calibration.jl` once before running these tests.
include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
```

---

### IN-03: `eps()` guard in `_smooth_field` denominator is too small

**File:** `spike/simulator/forward.jl:86`

**Issue:** `(f .- mean(f)) ./ (std(f) + eps())` uses `eps(Float64) ≈ 2.2e-16` as the
zero-denominator floor. If `std(f)` were anomalously close to zero, dividing by eps()
produces values ≈ 4.5×10^15, which propagates to Poisson rates of 50×4.5e15 — far
outside the range of practical integer arithmetic. Unreachable in production (Gaussian
noise over a 256×256 grid will always have std ≈ 1), but the guard is fragile if the
function is ever called with a smaller or synthetic field.

**Fix:**
```julia
return (f .- mean(f)) ./ max(std(f), 1e-10)
```

---

_Reviewed: 2026-06-26_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
