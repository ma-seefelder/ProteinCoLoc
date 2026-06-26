# ProteinCoLoc v2.0 — AmortizedColoc

## What This Is

ProteinCoLoc is a Julia package for Bayesian colocalization analysis of multi-channel
fluorescence microscopy images (basis for the author's *Scientific Reports* publication).
**v2.0 ("AmortizedColoc")** extends it with **amortized simulation-based inference (SBI)**:
a 2-channel 2D physics forward-simulator trains a neural posterior estimator (NPE) and a
neural ratio estimator (NRE) on the package's *existing* patch-correlation summary statistic,
delivering calibrated colocalization posteriors and Bayes factors in **milliseconds per
dataset** instead of per-dataset ADVI minutes — with a Simulation-Based-Calibration (SBC)
coverage proof and an out-of-distribution (OOD) / misspecification flag.

The work begins as a **strictly decoupled spike** (`spike/`, own `Project.toml`, no edits to
`src/`) that proves the core principle and ends in a Go/No-Go memo. If the spike succeeds,
v2.0 **productionizes** the amortized inference into the main `src/` package.

## Core Value

A trained NPE/NRE produces **calibrated, amortized colocalization inference** (posterior +
Bayes factor) in a single forward pass, **>100× faster than the existing per-dataset ADVI**,
with a **demonstrated SBC/coverage calibration proof** and an honest OOD flag — the first
colocalization tool to combine amortized inference, calibration proof, and misspecification
detection.

## Requirements

### Validated

<!-- Inferred from the existing published codebase — these already work and are relied upon. -->

- ✓ Patch-based colocalization with hierarchical Student-t Turing model + ADVI (`src/colocalization.jl`) — existing
- ✓ Patch-correlation summary statistic: `correlation(x,y;method)`, `patch()`, `_prepare_data()` (`src/colocalization.jl`) — existing
- ✓ KDE-based Bayes factor on Δρ via `compute_BayesFactor()` (`src/bayes.jl`) — existing
- ✓ Image model + masking: `MultiChannelImage`, `_apply_mask!`, `otsu_thresholds` (`src/LoadImages.jl`) — existing
- ✓ `AbstractMultiChannelImage` interface + accessors (added 2026-02) — existing

### Active

<!-- v2.0 scope. The spike (AP1–AP7) is the gate; productionization is conditional on Go. -->

**Spike — amortized SBI proof of principle (decoupled in `spike/`):**

- [ ] Isolated `spike/` environment with own `Project.toml` + green NeuralEstimators smoke test (AP1)
- [ ] Physics forward-simulator `simulate_pair(θ) → MultiChannelImage` (ρ_true, spillover, autofluorescence, label efficiency, sub-pixel shift, PSF, noise) consistent with Turing priors (AP2)
- [ ] Training-data generator: θ~π → simulate → existing summary statistic → cached dataset (AP3)
- [ ] NPE trained, inferring **both ρ_true and Δρ**; benchmarked vs. real `colocalization()` ADVI (RMSE, interval width, **wall-clock speedup**) (AP4)
- [ ] **Summary-statistic ablation**: compare minimal patch-corr vector vs. + Manders/moments for efficacy (AP4-adjacent)
- [ ] **SBC / coverage** proof: rank histograms, KS/χ² uniformity, coverage curve, ECE/MCE traffic-light; **tune until calibrated** (AP5)
- [ ] **Amortized Bayes factor** via NRE/Evidence-Network, validated against `compute_BayesFactor()` (AP6)
- [ ] **OOD / misspecification flag** — triggers on misspecified images, quiet in-distribution (AP6)
- [ ] Reproducible `spike/demo.jl` (fixed seed) + Go/No-Go memo with metrics + full-build-out decision (AP7)

**Productionization — conditional on a Go decision:**

- [ ] Integrate amortized NPE/NRE inference into `src/` as a shipped feature of the main package
- [ ] Make `num_patches` (patch grid) **user-definable** in the productionized API (8×8 was the spike default)

**Backend evaluation (open consideration):**

- [ ] Evaluate migrating the Bayesian backend **Turing → RxInfer.jl** for performance — feasibility against the hierarchical Student-t model (message-passing/conjugacy constraints), benchmarked vs. current ADVI and vs. amortized SBI; decide whether it is complementary to, or made redundant by, amortized inference

### Out of Scope

<!-- Deliberate boundaries for the spike; some revisited only at full build-out. -->

- Full 4-level hierarchy (Pixel→Segment→Celltype→Neighborhood) — spike proves 2-channel 2D first
- 3D / Z-stacks — deferred to full build-out
- Higher-order chromatic aberration — beyond spike simulator realism
- GUI — not needed for the inference claim
- Multi→2-channel generalization — spike fixes 2 channels
- Edits to `src/bayes.jl`, `src/colocalization.jl`, or the main pipeline during the spike — strict decoupling
- Wet-lab data, AlphaFold3/Boltz/Chai dependencies — fully solo, simulation-only

## Context

- **Brownfield**: extends a published Julia package. The spike *reuses* existing functions
  (`correlation`, `patch`, `_prepare_data`, `compute_BayesFactor`, `MultiChannelImage`,
  `otsu_thresholds`) as the NPE input and ADVI baseline — no rewriting.
- **Provenance preserved**: `v1.0.2-scirep` tag (commit `65aad18`) freezes the Scientific
  Reports publication code state; `v2.0-baseline` tag (HEAD `0d9b87c`) freezes the starting
  point. Both are local until pushed.
- **Stack**: Julia-native — NeuralEstimators.jl (`PosteriorEstimator` + `RatioEstimator`),
  Flux.jl backend, Images.jl/ImageFiltering.jl for PSF/noise, HypothesisTests.jl for SBC
  uniformity. BayesFlow (Python via PythonCall/CondaPkg) is a documented fallback only.
- **Reference assets** from a sibling project (`BayesInteractomics`): `_bin_calibration` +
  `CalibrationResult` (ECE/MCE + traffic-light) and the parametric simulation engine
  (sweep/coverage, JLD2 cache) as templates.
- **Timeline**: ~4 weeks to spike Go/No-Go; lowest priority, pausable (must not conflict with
  manuscript submissions).

## Constraints

- **Decoupling**: Spike lives entirely in `ProteinCoLoc/spike/` with its own `Project.toml`/`Manifest.toml`; the main package, `src/`, and both manuscript pipelines must remain **provably untouched** during the spike.
- **Compute**: CPU-only is the portable, reproducible baseline (2D is small enough). A **GPU is available** — CUDA.jl is an *optional accelerator* for training, not a requirement; the spike must run without it.
- **Reproducibility**: Everything reproducible from `spike/demo.jl` with a fixed (Random123) seed.
- **Platform**: Windows-tauglich — Julia-native stack chosen partly to avoid Windows CUDA/Flux friction; GPU use must degrade gracefully to CPU.
- **Summary dimension**: Fixed (8×8 patch grid) during the spike; user-definable only after productionization.
- **Priors**: Simulator prior π(θ) must stay consistent with the existing Turing `@model` prior ranges (μ/ν/σ/τ) for ADVI comparability.

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Tag `v1.0.2-scirep` (65aad18) + `v2.0-baseline` (HEAD) before any v2.0 work | Freeze the Scientific Reports publication state and the dev baseline separately | ✓ Good (tags created, push pending) |
| v2.0 = spike **+ productionize** (conditional on Go) | User wants a shipped feature, not just research, if the spike validates | — Pending (gated at M5) |
| Infer **both** ρ_true and Δρ | Δρ for Bayes-factor comparison, ρ_true for SBC | — Pending |
| Build **both** NPE + NRE/Evidence-Network | Amortized BF is a core publishable claim | — Pending |
| **Compare** summary statistics (minimal patch-corr vs + Manders/moments) | Spike should measure summary-stat efficacy, not assume it | — Pending |
| `num_patches` fixed 8×8 in spike, user-definable in production | Fixed summary dim needed for the NPE; flexibility belongs in the shipped API | — Pending |
| Training budget: start 50k, scale to 200k as M2 accuracy demands | Faster first loop; GPU available if needed | — Pending |
| SBC miscalibration → **tune until calibrated** | Prioritize achieving nominal coverage before reporting | — Pending |
| NeuralEstimators.jl default, BayesFlow only fallback | Julia-native, Windows-tauglich, single package covers NPE+NRE | — Pending |
| Evaluate Turing→RxInfer.jl backend migration | Potential baseline-inference speedup; may be complementary to or redundant with amortized SBI | — Pending (research/feasibility) |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd:transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd:complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

## Current State

- **Phase 1 (Environment + Smoke Gate) — Complete (2026-06-26).** Isolated `spike/` env stands up the pre-1.0 NeuralEstimators v0.2.1 + Flux v0.16.10 stack on Julia 1.12.6; the CPU-only NPE smoke is green AND correct (recovered posterior mean within tolerance). Validated: ENV-01, ENV-02, ENV-03, ENV-04. Root baseline frozen at commit `f581d95`; `src/` and root manifests provably untouched. Parent coupling uses the D-01 `include()` fallback (Pkg.develop silently downgraded NeuralEstimators 0.2.1→0.1.4 against the parent's heavy tree — recorded for Phase 4).
- **Next:** Phase 2 — Forward Simulator + Summary Contract.

---
*Last updated: 2026-06-26 after Phase 1 completion*
