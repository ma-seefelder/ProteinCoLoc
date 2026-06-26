---
phase: 02-forward-simulator-summary-contract
plan: 02
subsystem: simulation
tags: [julia, forward-simulator, sbi, microscopy, distributions, imagefiltering, warp, shared-latent, softplus]

# Dependency graph
requires:
  - phase: 02-01
    provides: "spike/contract.jl summary contract (build_mci/summary/induced_mu) over the UNCHANGED src/ patch()/correlation()"
  - phase: 01-environment-smoke-gate
    provides: "frozen spike env + baseline f581d95 + Phase-1 CUDA/NeuralEstimators-pin guard"
provides:
  - "spike/simulator/forward.jl: simulate_pair(rng, θ; imsize) — the 7-stage D-06 forward physics pipeline (SIM-01)"
  - "Shared-latent correlated smooth Gaussian-field stage-1 generator (D-15): standardized smooth latent + softplus mix with sign(ρ) flip; monotone ρ_true→induced-correlation incl. anti-correlation"
  - "SIM-01/SIM-03/D-15 testsets in spike/test/test_simulator.jl proving the simulator against the frozen summary contract"
affects: [02-03, prior-calibration, training-data-generation, npe-nre-training, sbc-calibration]

# Tech tracking
tech-stack:
  added: [ImageTransformations, CoordinateTransformations, Interpolations]
  patterns:
    - "Shared-latent standardized-smooth-field correlated generator (softplus(√|ρ|·L + √(1−|ρ|)·ε) with sign(ρ) channel-2 flip)"
    - "Explicit rng::AbstractRNG threaded through every stochastic stage (D-14)"
    - "Fixed-Gaussian-PSF nuisance via module constant, NOT a θ component (D-05)"
    - "θ/imsize validated at simulator entry (finite, ρ_true∈[-1,1], dims≥8)"

key-files:
  created: [spike/simulator/forward.jl]
  modified: [spike/test/test_simulator.jl, spike/Project.toml, spike/Manifest.toml]

key-decisions:
  - "σ_psf=1.3, STRUCT_σ=6.0 (≫σ_psf, the D-15 wash-out guard), GAIN=50, READ_NOISE=0.01, BG_FLOOR=0.02"
  - "Standardize each smoothed latent field to unit variance — required to keep corr=ρ through softplus and prevent induced-μ collapse under PSF+noise"
  - "Promote ImageTransformations/CoordinateTransformations/Interpolations to spike direct deps (RESEARCH A2) since Images does not re-export warp/Translation/BSpline; versions unchanged, NeuralEstimators stays v0.2.1"

patterns-established:
  - "Pattern: shared-latent standardized-smooth-field bivariate intensity generator (D-15), replacing the rejected puncta/spot-mixture"
  - "Pattern: 7-stage D-06 forward pipeline returning Vector{Matrix{Float64}} into the frozen Wave-1 contract"

requirements-completed: [SIM-01, SIM-03]

# Metrics
duration: 11min
completed: 2026-06-26
---

# Phase 2 Plan 02: Forward Physics Simulator (SIM-01) Summary

**`simulate_pair(rng, θ; imsize)` runs the 7-stage D-06 microscopy forward pipeline on a spike-validated shared-latent standardized-smooth-field generator (D-15), driving a monotone, sign-correct ρ_true→induced-patch-correlation (−0.43…+0.60 across the sweep) that flows unchanged through the Wave-1 summary contract (SIM-03).**

## Performance

- **Duration:** ~11 min
- **Started:** 2026-06-26T19:07:43Z
- **Completed:** 2026-06-26T19:18:54Z
- **Tasks:** 2
- **Files modified:** 4 (1 created, 3 modified)

## Accomplishments
- `spike/simulator/forward.jl`: `simulate_pair(rng::AbstractRNG, θ; imsize=(256,256))` executing all 7 D-06 stages in order — shared-latent correlated fields (D-15) → Bernoulli thinning → fixed Gaussian PSF → directional 2×2 spillover → autofluorescence+BG_FLOOR offset → sub-pixel `warp` shift on ch2 → Poisson(shot)+Gaussian(read) noise — returning a 2-element `Vector{Matrix{Float64}}`, all finite and ≥ 0 with small-positive background.
- Stage 1 is the spike-validated shared-latent **standardized** smooth Gaussian-field generator: `c1=softplus(√|ρ|·L+√(1−|ρ|)·ε₁)`, `c2=softplus(sign(ρ)·√|ρ|·L+√(1−|ρ|)·ε₂)`; the `sign(ρ)` flip makes anti-correlation reachable. The rejected puncta/spot-mixture is NOT used (D-15).
- Proven monotone ρ_true→induced-correlation through the REAL summary: induced-μ = −0.43 / −0.11 / +0.11 / +0.31 / +0.60 across ρ_true ∈ {−0.7,−0.3,0,0.3,0.7}, with the negative tail genuinely anti-correlated.
- SIM-01 + SIM-03-on-output + D-15-monotonicity testsets added; full runner `julia --project=spike spike/test/runtests.jl` exits 0 (5+14+10+10+5 = 44 passing tests) with the Phase-1 CUDA/NeuralEstimators-pin guard intact (NeuralEstimators v0.2.1, no CUDA).

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement simulate_pair — the 7-stage forward pipeline** - `caa8017` (feat)
2. **Task 2: SIM-01 + SIM-03 assertions on real simulator output** - `e11f8b2` (test)

**Plan metadata:** (final docs commit below)

## Files Created/Modified
- `spike/simulator/forward.jl` - The SIM-01 forward physics simulator: `simulate_pair`, the shared-latent generator (D-15), 7 D-06 stages, fixed-PSF nuisance, θ validation, explicit rng threading.
- `spike/test/test_simulator.jl` - Extended with SIM-01 (shape/finiteness/validation/determinism), SIM-03-on-output (valid MultiChannelImage, ≤2 missing), and D-15 (monotone + anti-correlation) testsets.
- `spike/Project.toml` - Promoted ImageTransformations/CoordinateTransformations/Interpolations to direct deps for `warp`/`Translation`/`BSpline`.
- `spike/Manifest.toml` - Re-resolved (project hash only; package versions unchanged, NeuralEstimators stays v0.2.1).

## Decisions Made
- **Standardize smoothed latent fields to unit variance.** Gaussian smoothing at STRUCT_σ=6 collapses the raw field variance to ≈1/(4π·36) (sd≈0.05), which buried the correlation signal under Poisson noise and flattened induced-μ to ≈0 with the negative tail unreachable. Standardizing each field restores full amplitude and makes the linear mix carry correlation exactly `sign(ρ)·a²=ρ` into softplus. This is the concrete mechanism behind the D-15 "STRUCT_σ ≫ σ_psf" guard.
- **Fixed nuisance constants** (Claude's discretion, D-05/D-07 ranges): σ_psf=1.3 px, STRUCT_σ=6.0 px, GAIN=50, READ_NOISE=0.01, BG_FLOOR=0.02.
- **warp** with `method=BSpline(Linear())`, `fillvalue=BG_FLOOR` (small-positive, not 0.0) — satisfies both the no-NaN-border (T-02-NAN) and background-not-zero (T-02-Z0) traps simultaneously.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Promoted three image packages to spike direct deps**
- **Found during:** Task 1 (warp sub-pixel shift stage)
- **Issue:** The plan's `warp(ch2, Translation(...); method=BSpline(Linear()))` failed — `Images` does not re-export `Translation`/`warp`/`BSpline` (`UndefVarError(:Translation)`), so they were not callable as written.
- **Fix:** Promoted `ImageTransformations`, `CoordinateTransformations`, and `Interpolations` (already pinned transitive deps of `Images` in the Manifest) to `spike/Project.toml` direct deps and added explicit `using` statements (RESEARCH A2's sanctioned fallback). `Pkg.resolve()` reported no package add/remove and no version change; the resolve-risk gate confirmed NeuralEstimators stays v0.2.1 with no CUDA.
- **Files modified:** spike/Project.toml, spike/Manifest.toml, spike/simulator/forward.jl
- **Verification:** Full `runtests.jl` exits 0; Phase-1 CUDA/NeuralEstimators-pin testset still passes.
- **Committed in:** caa8017 (Task 1 commit)

**2. [Rule 1 - Bug] Standardized the smoothed latent fields to prevent induced-μ collapse**
- **Found during:** Task 1 (stage-1 generator verification)
- **Issue:** A first pass without field standardization produced compressed, sign-wrong induced-μ (0.0…0.027 across the full ρ_true sweep; the negative-ρ point was 0.0, not anti-correlated) — exactly the D-15 wash-out failure mode, because smoothing collapses field variance below the noise floor.
- **Fix:** `_smooth_field` now standardizes each smoothed field to zero-mean/unit-variance before the softplus mix, restoring amplitude and preserving `corr=ρ`.
- **Files modified:** spike/simulator/forward.jl
- **Verification:** Induced-μ sweep now spans −0.43…+0.60, monotone, anti-correlated at negative ρ_true (D-15 testset green).
- **Committed in:** caa8017 (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (1 blocking, 1 bug)
**Impact on plan:** Both essential to deliver the plan's specified behavior (warp shift + the D-15 monotone/sign-correct generator). No scope creep; isolation preserved (only `spike/` touched, root clean vs f581d95).

## Issues Encountered
- None beyond the two deviations above. The NPE smoke in the shared runner takes ~3.5 min (Flux training) but passes; total suite green.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- The forward simulator is ready as the generative substrate for prior calibration (SIM-02, induced-μ → Turing μ-prior fit) and the SIM-04 plausibility gate (the remaining Phase-2 work in subsequent plans).
- `sample_prior` / `calibration.jl` and the CairoMakie plausibility figures are NOT yet built (out of scope for this plan).
- Decoupling invariant intact: `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` is clean.

---
*Phase: 02-forward-simulator-summary-contract*
*Completed: 2026-06-26*

## Self-Check: PASSED
- FOUND: spike/simulator/forward.jl
- FOUND: spike/test/test_simulator.jl
- FOUND: .planning/phases/02-forward-simulator-summary-contract/02-02-SUMMARY.md
- FOUND commit: caa8017 (Task 1, feat)
- FOUND commit: e11f8b2 (Task 2, test)
