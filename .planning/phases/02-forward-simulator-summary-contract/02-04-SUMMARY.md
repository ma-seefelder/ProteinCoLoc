---
phase: 02-forward-simulator-summary-contract
plan: 04
subsystem: simulation
tags: [julia, sbi, plausibility-gate, cairomakie, spearman-monotonicity, paired-ttest, hypothesistests, statsbase, headless-figures, sim-04]

# Dependency graph
requires:
  - phase: 02-03
    provides: "spike/simulator/prior.jl sample_prior(rng) + frozen ĝ (induced-μ-consistent generative prior driving the plausibility sweeps)"
  - phase: 02-02
    provides: "spike/simulator/forward.jl simulate_pair (shared-latent D-15 generator) — the physics whose plausibility this gate proves"
  - phase: 02-01
    provides: "spike/contract.jl summary contract (build_mci/summary/induced_mu) over the UNCHANGED src/ patch()/correlation()"
  - phase: 01-environment-smoke-gate
    provides: "frozen spike env + baseline f581d95 + Phase-1 CUDA/NeuralEstimators-pin guard"
provides:
  - "spike/02_simulator_demo.jl: seeded headless demo regenerating the plausibility figures over the full sample_prior→simulate_pair→build_mci→summary chain"
  - "spike/figures/plausibility.{png,pdf}: rendered SIM-04 evidence (image panels, ρ_true monotonicity sweep, induced-μ vs μ-prior overlay, spillover/shift effect curves)"
  - "SIM-04 quantitative plausibility gate in spike/test/test_simulator.jl (monotonicity Spearman ≥ 0.95 + paired spillover/shift effects + valid MultiChannelImage)"
affects: [training-data-generation, npe-nre-training, go-no-go-memo, phase-3]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Headless CairoMakie evidence figures (CairoMakie.activate!(), Figure/heatmap/scatterlines/hist, save() to png+pdf) — spike-local, no OpenGL/display (D-13)"
    - "Quantitative plausibility gate: deterministic Spearman monotonicity + PAIRED perturbation effects (shared-seed pairing, OneSampleTTest p<0.05 beyond MC noise) instead of a visual check (D-12)"
    - "Script self-assert idiom (@assert isfile(png/pdf)) so the demo's automated verify is PowerShell-safe (no POSIX test -f)"

key-files:
  created: [spike/02_simulator_demo.jl, spike/figures/plausibility.png, spike/figures/plausibility.pdf]
  modified: [spike/test/test_simulator.jl, .gitignore]

key-decisions:
  - "Spearman monotonicity threshold fixed at 0.95 (RESEARCH default / A4 discretion); the deterministic 15-pt averaged sweep and the demo both reach Spearman 1.0, so 0.95 has full headroom and cannot flake"
  - "Perturbation effects asserted as PAIRED differences over shared Xoshiro seeds (only the spillover/shift param differs between the two runs) with OneSampleTTest p<0.05 — not single-seed flukes (T-02-GATE)"
  - "Effect-size floors δ_spill=0.02, δ_shift=0.01 (well below the observed +0.18 / +0.095 demo effects) to stay robustly beyond Monte-Carlo noise"
  - ".gitignore: scoped negation (!spike/figures/plausibility.png) to track the PNG deliverable past the repo-wide *.png ignore (D-13 evidence artifact); .gitignore is not src/ or a manifest, so decoupling is preserved"

patterns-established:
  - "Pattern: SIM-04 quantitative gate is deterministic under fixed Xoshiro seeds (D-14) so the plausibility proof is reproducible, not flaky"
  - "Pattern: 6 nuisances pinned at prior MEDIANS (spillover 0.1, autofluor 0.05, label_eff 0.8, shift 0, noise 0.5) for the monotonicity sweep to isolate ρ_true"

requirements-completed: [SIM-04]

# Metrics
duration: 22min
completed: 2026-06-26
---

# Phase 2 Plan 04: SIM-04 Plausibility Gate + CairoMakie Evidence Figures Summary

**A deterministic quantitative plausibility gate proves the shared-latent forward simulator behaves like real microscopy — ρ_true↑ → mean-patch-correlation↑ at Spearman 1.0 (≥0.95 bar) over a 15-point sweep, spillover RAISES (+0.18) and sub-pixel shift LOWERS (−0.095) correlation as paired effects significant beyond Monte-Carlo noise, all output stays a valid 8×8 MultiChannelImage — and a seeded headless `spike/02_simulator_demo.jl` regenerates the CairoMakie evidence figures under `spike/figures/`, closing Phase 2 with decoupling byte-identical to baseline f581d95.**

## Performance

- **Duration:** ~22 min
- **Started:** 2026-06-26T22:00:00Z
- **Completed:** 2026-06-26T22:22:00Z
- **Tasks:** 2
- **Files modified:** 5 (3 created, 2 modified)

## Accomplishments
- `spike/02_simulator_demo.jl`: a seeded (`Xoshiro(2026)`, D-14), CPU-only, **headless** demo that runs the full `sample_prior → simulate_pair → build_mci → summary` chain and renders a 2×4 `Figure`: high-ρ and low-ρ channel-pair `heatmap`s, the ρ_true→mean-patch-corr `scatterlines` sweep, the induced-μ `hist` overlaid on the Turing μ-prior density (`Truncated(Cauchy(0,0.3), GHAT_MU_MIN, GHAT_MU_MAX)`), and the spillover/|shift| effect curves. Saves `plausibility.{png,pdf}` and **self-asserts `@assert isfile(...)`** on both (PowerShell-safe).
- `spike/figures/plausibility.png` (raster) + `plausibility.pdf` (vector): the rendered SIM-04 evidence — verified visually (monotone sweep, induced-μ tracking the prior density, spillover-up/shift-down trends).
- **SIM-04 quantitative gate** added to `spike/test/test_simulator.jl` (21 assertions, all green):
  - **(a) monotonicity:** `corspearman(ρ_grid, mean_patch_corr) ≥ 0.95` over a deterministic 15-point sweep (N=3 MC replicates/point, nuisances at prior medians) — achieves **1.0**.
  - **(b) spillover:** paired `mean(Δ) > 0.02` and `OneSampleTTest p < 0.05` over 8 shared seeds (spillover 0.0 vs 0.2) — measurably RAISES correlation.
  - **(b) shift:** paired `mean(Δ) > 0.01` and `OneSampleTTest p < 0.05` over 8 shared seeds (no-shift vs |shift|≈1.27 px) — measurably LOWERS correlation.
  - **(c) validity:** both perturbation regimes still build a valid `MultiChannelImage` with ≤64 finite per-patch ρ and ≤2 missing (SIM-03 traps held).
- `julia --project=spike spike/test/runtests.jl` exits **0**: SIM-04 (+21), SIM-01/02/03, the D-15 monotonicity testset, and the Phase-1 CUDA-absence / NeuralEstimators-v0.2.1-pin resolve-risk guard all pass after CairoMakie/HypothesisTests entered the active env.
- Decoupling invariant preserved: `git diff --quiet f581d95 -- src/ Project.toml Manifest.toml` is clean (root byte-identical to baseline).

## Task Commits

Each task was committed atomically:

1. **Task 1: Seeded plausibility demo + CairoMakie figures** - `6eac6ad` (feat) — `spike/02_simulator_demo.jl` + `spike/figures/plausibility.{png,pdf}` + `.gitignore` negation
2. **Task 2: SIM-04 quantitative gate (monotonicity + perturbation + valid image)** - `55e330d` (test) — `@testset "SIM-04 plausibility"` in `spike/test/test_simulator.jl`

**Plan metadata:** final docs commit (this SUMMARY + STATE + ROADMAP + REQUIREMENTS).

## Files Created/Modified
- `spike/02_simulator_demo.jl` - Seeded headless demo: full sample→simulate→summary chain, 2×4 CairoMakie figure, save to png+pdf, `@assert isfile` self-check, console echo of Spearman/spillover/shift effect sizes.
- `spike/figures/plausibility.png` - Rendered SIM-04 evidence (raster).
- `spike/figures/plausibility.pdf` - Rendered SIM-04 evidence (vector).
- `spike/test/test_simulator.jl` - Added `using HypothesisTests` + the "SIM-04 plausibility" testset (3 sub-gates, 21 assertions); SIM-01/02/03 testsets and Phase-1 guards untouched.
- `.gitignore` - Scoped `!spike/figures/plausibility.png` negation so the evidence figure tracks past the repo-wide `*.png` ignore.

## Decisions Made
- **Spearman threshold = 0.95.** RESEARCH §Plausibility Gate default; A4/D-12 discretion. The deterministic averaged sweep and the demo both hit 1.0, so 0.95 has full headroom and — being seed-fixed — cannot flake. Recorded here per T-02-GATE.
- **Paired perturbation tests over shared seeds.** For each seed a single `Xoshiro(s)` drives both runs, so the ONLY difference is the spillover/shift parameter — a true paired difference. `OneSampleTTest` on the 8 differences gives `p<0.05` significance beyond MC noise (T-02-GATE), not a single-seed fluke. Direction asserted explicitly (spillover up, shift down).
- **Effect-size floors δ_spill=0.02, δ_shift=0.01.** Conservative — well under the observed +0.18 / +0.095 effects — so the gate fails only if the physics genuinely stops responding, not on noise.
- **Qualified `CairoMakie.Axis` and `DataAspect`.** `Axis` is exported by both Images and Makie (ambiguity error); qualifying the Makie symbols resolves it — the only adjustment to the planned CairoMakie recipe.

## Deviations from Plan

None - plan executed exactly as written. Both tasks, the headless demo, the rendered figures, the three SIM-04 sub-gates, and all acceptance criteria were delivered as specified. Two minor, non-architectural adjustments were made inline and are NOT scope deviations: (1) qualifying `CairoMakie.Axis` to resolve the Images/Makie name ambiguity (a Julia name-resolution fix, not a design change), and (2) a scoped `.gitignore` negation so the required PNG deliverable could be tracked past the repo-wide `*.png` ignore. Neither touches src/ or the manifests.

## Issues Encountered
- **`Axis` ambiguity:** `using Images` (via contract.jl) and `using CairoMakie` both export `Axis`, raising `UndefVarError: Axis not defined (ambiguity)`. Resolved by qualifying as `CairoMakie.Axis` at the four axis call sites; `DataAspect`/`heatmap`/`scatterlines`/`hist`/`save` are Makie-unique and needed no qualification.
- **`*.png` gitignore:** the repo's root `.gitignore` ignores all `*.png`. Since `plausibility.png` is a required `must_haves` artifact (D-13 evidence), added a scoped `!spike/figures/plausibility.png` negation rather than force-adding, so future regenerations track cleanly.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- **Phase 2 is complete.** SIM-01 (forward simulator), SIM-02 (induced-μ-consistent prior), SIM-03 (summary contract), and SIM-04 (quantitative plausibility gate) are all green and re-runnable via a single `julia --project=spike spike/test/runtests.jl` (exit 0).
- The plausibility evidence (`spike/figures/plausibility.{png,pdf}`) is regenerable and tracked — ready for the Go/No-Go memo and the manuscript.
- The simulator + consistent prior + frozen summary contract are the validated generative substrate for Phase 3 training-data generation and NPE/NRE training.
- Decoupling invariant intact across the whole phase: `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` is clean.

---
*Phase: 02-forward-simulator-summary-contract*
*Completed: 2026-06-26*

## Self-Check: PASSED
- FOUND: spike/02_simulator_demo.jl
- FOUND: spike/figures/plausibility.png
- FOUND: spike/figures/plausibility.pdf
- FOUND: spike/test/test_simulator.jl
- FOUND: .planning/phases/02-forward-simulator-summary-contract/02-04-SUMMARY.md
- FOUND commit: 6eac6ad (Task 1, feat — demo + figures + .gitignore negation)
- FOUND commit: 55e330d (Task 2, test — SIM-04 plausibility gate)
- ROOT CLEAN vs f581d95 (decoupling invariant preserved)
