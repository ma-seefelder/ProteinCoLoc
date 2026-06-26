# Phase 2: Forward Simulator + Summary Contract - Discussion Log

> **Audit trail only.** Not consumed by research/planning/execution agents — decisions live in CONTEXT.md.

**Date:** 2026-06-26
**Phase:** 02-forward-simulator-summary-contract
**Mode:** discuss
**Areas discussed:** Prior consistency mapping, Physics fidelity, Image geometry, Validation + plot backend

## Scouting findings that framed the discussion
- The Turing `@model` (`src/bayes.jl:251–316`) models the **per-patch correlation summary** as a Student-t sample (`control/sample[idx] ~ TDist(ν)·σ + μ`), not images — so SIM-02 "consistency" is a physics-θ → summary-stat-μ mapping problem.
- Real paper data in-repo: `test/test_images/{positive,negative}/*_c{1,2,3}.tif`, native **1376×1028, 16-bit, 3-channel**.
- `MultiChannelImage.pixel_size` stores the **pixel-dimension tuple** (`src/LoadImages.jl:444`), not µm.

## Decisions

### Prior consistency mapping (SIM-02)
- Options: Direct prior + induced-μ check (recommended) | Full induced-μ calibration | Direct adoption only
- **Selected:** **Full induced-μ calibration** → D-01, D-02, D-03. User chose the most rigorous path: fit/transform the ρ_true prior so the simulator's induced μ distribution matches the Turing μ-prior, documented with evidence in spike/NOTES.md.

### Physics fidelity (SIM-01)
- Options: 7-param θ, PSF fixed, Poisson+Gaussian (recommended) | Same θ, Gaussian-only noise | Add PSF width to θ
- **Selected:** **7-param θ, PSF fixed Gaussian nuisance, Poisson+Gaussian noise** → D-04..D-07. Directional 2×2 spillover, Bernoulli-thinning label efficiency, roadmap pipeline order.

### Image geometry
- Options: Match real size (recommended) | Small fixed 256² | Configurable default 512²
- **Selected (with extension):** **Match the real size 1376×1028 AND sweep other sensible microscopy sizes** → D-08, D-09, D-10, D-11. Configurable dims, 8×8 patch grid fixed, pixel_size = pixel-dim tuple. The "sweep all sensible sizes" instruction makes the common-microscopy-size set a research item.

### Validation + plot backend (SIM-04)
- Options: Quantitative gate + CairoMakie (recommended) | Quantitative gate + Plots.jl | Visual-only + CairoMakie
- **Selected:** **Quantitative gate + CairoMakie** → D-12, D-13. Monotonicity (Spearman) + measurable spillover/shift perturbation; CairoMakie figures, spike-local only (does not touch src/plot.jl).

## Claude's Discretion / research items
- PSF width, noise constants, Spearman threshold, calibration-fit method, the exact sensible-size sweep set, internal spike layout.
- RNG threaded explicitly (D-14): stdlib Random in Phase 2, Random123 deferred to Phase 3.

## Deferred
- Random123 seeding → Phase 3.
- src/plot.jl main-package backend migration → separate post-spike decision.
- Training-data scale / NPE / ADVI / BF → Phases 3, 4, 6.

## Connection to prior conversation
- The CairoMakie choice (D-13) partially addresses the user's earlier GLMakie/Makie→Plots compile-time question — but only for the **decoupled spike**, not the frozen main package. The src/ backend question stays open and deferred.
