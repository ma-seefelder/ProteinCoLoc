---
phase: 02-forward-simulator-summary-contract
plan: 03
subsystem: simulation
tags: [julia, sbi, prior-calibration, induced-mu, isotonic-regression, pchip, wasserstein, ks, distributions, statsbase, turing-prior]

# Dependency graph
requires:
  - phase: 02-02
    provides: "spike/simulator/forward.jl simulate_pair (shared-latent D-15 generator) — the generative substrate the calibration sweeps ρ_true through"
  - phase: 02-01
    provides: "spike/contract.jl summary contract (build_mci/summary/induced_mu) + read-only load_tiff over the UNCHANGED src/ patch()/correlation()"
  - phase: 01-environment-smoke-gate
    provides: "frozen spike env + baseline f581d95 + Phase-1 CUDA/NeuralEstimators-pin guard"
provides:
  - "spike/simulator/calibration.jl: offline induced-μ sweep + isotonic monotone ĝ fit + W1/KS evidence + size-invariance + faithful real-anchor (SIM-02 calibration engine)"
  - "spike/simulator/ghat.jl: frozen ĝ : μ ↦ ρ_true (knots/coeffs + ghat(μ) evaluator + GHAT_MU_MIN/MAX realized range), the reproducibility artifact"
  - "spike/simulator/prior.jl: MU_PRIOR + sample_prior(rng) drawing ρ_true=ghat(μ*) so induced μ matches the Turing μ-prior (SIM-02 prior)"
  - "spike/NOTES.md §3: the SIM-02 calibration evidence (ĝ method, span-past-±0.9 rationale, match metric, size-invariance, real anchor + negative-tail verdict, σ/τ/ν check)"
  - "SIM-02 prior-consistency testset (ĝ monotonicity + induced-μ W1 gate) in spike/test/test_simulator.jl"
affects: [02-04, plausibility-gate, training-data-generation, npe-nre-training, advi-baseline, sbc-calibration]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Offline induced-μ calibration: sweep ρ_true through the REAL summary → isotonic (PAVA) monotone fit → clamped piecewise-linear inverse ĝ, frozen as a plain-Julia literal artifact"
    - "Transform sampling for prior consistency: draw μ*~Turing-prior, set ρ_true=ĝ(μ*) so induced μ matches by construction"
    - "Pre-declared KS/Wasserstein tolerance fixed BEFORE the run (T-02-CAL data-snooping guard); consistency scoped to the physically-realized μ range (D-16)"

key-files:
  created: [spike/simulator/calibration.jl, spike/simulator/ghat.jl, spike/simulator/prior.jl]
  modified: [spike/test/test_simulator.jl, spike/NOTES.md]

key-decisions:
  - "Isotonic regression (PAVA) + clamped piecewise-linear inverse for ĝ (guaranteed-monotone, no extra dep, sanctioned over a bespoke optimizer / NormalizingFlows.jl)"
  - "Sweep ρ_true∈[-0.99,0.99] (Open Q1 RESOLVED, past ±0.9) — generator induces only μ≈±0.76 at knob ±0.9; realized induced-μ range [-0.680,+0.847]"
  - "Pre-declared SIM02_W1_TOL=0.10 / SIM02_KS_TOL=0.20; achieved W1=0.052, KS=0.101 over the realized range"
  - "Negative μ-prior tail documented PRIOR-ONLY: faithful load_tiff real anchor gives positive μ=0.329, negative μ=0.248 (both >0) → real anti-correlation not physically realized (D-16)"
  - "Nuisance ranges (D-03): spillover U(0,0.2), autofluor U(0,0.1), label_eff U(0.6,1), shift U(-1,1), noise U(0,1) — induced per-patch-corr pooled SD 0.317 inside the Turing scale prior"

patterns-established:
  - "Pattern: frozen offline-calibration artifact (calibration.jl emits ghat.jl) included read-only by prior.jl — no JLD2 (deferred to Phase 3)"
  - "Pattern: induced-μ consistency measured through the UNCHANGED summary and scoped to the physically-realized range with a pre-declared distance tolerance"

requirements-completed: [SIM-02]

# Metrics
duration: 32min
completed: 2026-06-26
---

# Phase 2 Plan 03: Induced-μ Calibration + Consistent Prior (SIM-02) Summary

**A frozen monotone induced-μ inverse ĝ (isotonic-PAVA fit, swept past ±0.9 on the shared-latent generator through the UNCHANGED summary) makes `sample_prior(rng)`'s induced per-patch-correlation mean match the Turing μ-prior `Truncated(Cauchy(0,0.3),-1,1)` to Wasserstein-1 = 0.052 over the physically-realized range [-0.68,+0.85], with the negative tail honestly flagged prior-only against a faithful `load_tiff` real anchor (D-16).**

## Performance

- **Duration:** ~32 min (incl. the ~9 min one-time offline calibration run)
- **Started:** 2026-06-26T21:18:00Z
- **Completed:** 2026-06-26T21:50:00Z
- **Tasks:** 2 (Task 2 = TDD RED→GREEN→docs)
- **Files modified:** 5 (3 created, 2 modified)

## Accomplishments
- `spike/simulator/calibration.jl`: an offline `calibrate(rng)` that sweeps `ρ_true ∈ [-0.99, 0.99]` (25 pts, n_per=100, 512²) through `simulate_pair` → `build_mci` → `induced_mu` (the UNCHANGED summary) with the 6 nuisances drawn from their priors, fits a **monotone isotonic (PAVA)** map and inverts it to a clamped piecewise-linear `ĝ : μ ↦ ρ_true`, and emits the frozen artifact + all evidence.
- `spike/simulator/ghat.jl` (auto-generated, committed reproducibility artifact): `GHAT_MU_KNOTS`/`GHAT_RHO_KNOTS`, `GHAT_MU_MIN/MAX = [-0.680, 0.847]`, and a `ghat(μ)` evaluator. `corspearman(ρ_grid, E[μ]) = 1.0`.
- `spike/simulator/prior.jl`: `MU_PRIOR` (verbatim Turing μ-prior) + `sample_prior(rng)` drawing `ρ_true = ghat(μ*)` and the 6 nuisances; deterministic under a fixed rng (D-14); pre-declared `SIM02_W1_TOL = 0.10`.
- **SIM-02 pass:** induced-μ Wasserstein-1 = **0.052 < 0.10** (KS = 0.101 < 0.20), evaluated over the physically-realized μ range (D-16). **Size-invariance:** max |Δ E[μ]| vs 512² across {256²,512²,1024²} = **0.032**.
- **Faithful real anchor (D-16):** via the package's own `load_tiff` (NOT a hand-rolled RGB→luminance reduction), real positive μ = **0.329**, negative μ = **0.248** — both positive → negative μ-prior tail documented **PRIOR-ONLY**; consistency scoped to the realized range.
- **σ/τ/ν consistency (D-03):** induced per-patch-corr pooled SD = 0.317 (inside `Truncated(Cauchy(0.1,0.3),1e-4,1)`); excess kurtosis ≈ -0.05 (near-Gaussian — honestly reported, consistent with a large-ν Exponential).
- `spike/NOTES.md §3` completed with the full SIM-02 evidence; `julia --project=spike spike/test/runtests.jl` exits 0 (SIM-02 testset +10, all prior testsets, Phase-1 CUDA/NeuralEstimators-pin guard intact).

## Task Commits

Each task was committed atomically:

1. **Task 1: Offline induced-μ sweep + monotone ĝ + real anchor + evidence** - `ef51011` (feat) — also ran `calibration.jl` to freeze `ghat.jl`
2. **Task 2: sample_prior + NOTES §3 + SIM-02 gate (TDD)**
   - RED: `edf3df1` (test) — failing SIM-02 prior-consistency testset (prior.jl absent)
   - GREEN: `0959b49` (feat) — `prior.jl` with `sample_prior` via frozen ĝ → testset green
   - docs: `3e45819` (docs) — completed NOTES §3 calibration evidence

**Plan metadata:** final docs commit (this SUMMARY + STATE + ROADMAP + REQUIREMENTS).

## Files Created/Modified
- `spike/simulator/calibration.jl` - SIM-02 calibration engine: sweep, PAVA isotonic fit + inverse, W1/KS metric, size-invariance, faithful real anchor, σ/τ/ν diagnostics, `emit_ghat` writer, script entry point.
- `spike/simulator/ghat.jl` - Frozen ĝ knots/coeffs + `ghat(μ)` + realized-range consts; auto-generated by calibration.jl (committed as the reproducibility artifact).
- `spike/simulator/prior.jl` - `MU_PRIOR`, `sample_prior(rng)`, nuisance prior consts (mirror of calibration.jl), `SIM02_W1_TOL`.
- `spike/test/test_simulator.jl` - Added the "SIM-02 prior consistency" testset (ĝ monotonicity + clamped tails, 7-field θ + determinism, induced-μ W1 < tol over the realized range).
- `spike/NOTES.md` - §3 rewritten from the Phase-1 seed into the full SIM-02 documentation deliverable.

## Decisions Made
- **ĝ fit = isotonic regression (PAVA) + clamped piecewise-linear inverse.** Guaranteed-monotone, no new dependency, sanctioned by 02-RESEARCH §Don't Hand-Roll over a bespoke optimizer or the excluded NormalizingFlows.jl. PAVA is the textbook isotonic algorithm, not bespoke code.
- **Sweep span [-0.99, 0.99] (past ±0.9, Open Q1 RESOLVED).** The shared-latent generator induces only μ≈±0.76 at knob ±0.9 post-degrade (D-15), so the knob extends past ±0.9 to reach the μ-prior tails. Realized induced-μ range is asymmetric [-0.680, +0.847] (softplus + BG_FLOOR push the positive tail further).
- **Pre-declared tolerance BEFORE the run (T-02-CAL).** `SIM02_W1_TOL=0.10` / `SIM02_KS_TOL=0.20`, justified by the transform-sampling argument (induced μ is an unbiased, conditional-noise-widened copy of the target). Not tuned to pass; achieved 0.052 / 0.101.
- **Negative tail prior-only (D-16).** The faithful `load_tiff` anchor shows real fluorescence never goes anti-correlated (positive control 0.33, negative control 0.25), so the D-01 consistency claim is scoped to the physically-realized μ range and the negative μ extreme is documented prior-only.

## Deviations from Plan

None - plan executed exactly as written. Both tasks, the offline calibration run, the frozen artifact, the TDD RED→GREEN cycle, the NOTES §3 evidence, and all acceptance criteria were delivered as specified; no auto-fixes (Rules 1-4) were required.

## Issues Encountered
- The one-time offline `calibration.jl` run takes ~9 min (2500+ `simulate_pair` calls at 512² for the sweep + metric + size-invariance). Expected for an offline freeze step; `ghat.jl` is committed so the gate (`runtests.jl`) never re-runs it.
- The induced per-patch-correlation distribution is near-Gaussian (excess kurtosis ≈ -0.05), not heavy-tailed. Reported honestly in NOTES §3: this does not contradict the finite-ν Exponential prior (which permits near-Gaussian at large ν); ν is "induced and checked", not set — and the SIM-02 gate is on μ, which passes.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- `sample_prior(rng)` + frozen `ĝ` are ready as the consistent generative prior for Plan 02-04 (the SIM-04 plausibility gate + CairoMakie figures) and downstream Phase-3 training-data generation.
- The induced-μ↔Turing-μ-prior consistency is proven and re-asserted in the test gate; σ/τ/ν consistency is documented.
- Decoupling invariant intact: `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` is clean (root byte-identical to baseline).

---
*Phase: 02-forward-simulator-summary-contract*
*Completed: 2026-06-26*

## Self-Check: PASSED
- FOUND: spike/simulator/calibration.jl
- FOUND: spike/simulator/ghat.jl
- FOUND: spike/simulator/prior.jl
- FOUND: spike/NOTES.md
- FOUND: spike/test/test_simulator.jl
- FOUND: .planning/phases/02-forward-simulator-summary-contract/02-03-SUMMARY.md
- FOUND commit: ef51011 (Task 1, feat — calibration + frozen ghat.jl)
- FOUND commit: edf3df1 (Task 2 RED, test)
- FOUND commit: 0959b49 (Task 2 GREEN, feat — prior.jl)
- FOUND commit: 3e45819 (Task 2 docs, NOTES §3)
- ROOT CLEAN vs f581d95 (decoupling invariant preserved)
