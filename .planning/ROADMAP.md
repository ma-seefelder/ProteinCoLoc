# Roadmap: ProteinCoLoc v2.0 — AmortizedColoc

## Overview

v2.0 proves and (conditionally) ships amortized simulation-based inference for colocalization. The work is a strictly decoupled spike (`spike/`, own `Project.toml`, no `src/` edits) built bottom-up as a 7-layer horizontal pipeline: a NeuralEstimators+Flux smoke gate de-risks the pre-1.0 API before any investment; a physics forward-simulator emits real `MultiChannelImage`s that the *unmodified* existing summary functions ingest; a version-guarded training cache feeds an NPE benchmarked at >100x ADVI speedup; a shared θ*~π→simulate→infer harness fans out into the SBC/BF/OOD validation trifecta; a reproducible demo closes with a Go/No-Go memo. Only on a Go decision does the final phase touch `src/`, integrating amortized inference behind a coexisting backend contract. Every phase is an independently runnable checkpoint, resumable across pauses.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Environment + Smoke Gate** - Isolated spike env, pinned Manifest, green NeuralEstimators/Flux toy NPE — the hard gate before all downstream investment (completed 2026-06-26)
- [ ] **Phase 2: Forward Simulator + Summary Contract** - `simulate_pair(θ) → MultiChannelImage` consumed unchanged by the existing summary functions, prior-consistent with the Turing model
- [ ] **Phase 3: Training-Data Pipeline** - θ~π → simulate → summary → standardized fixed-dim vectors, version-guarded JLD2 cache with leak-free split discipline
- [ ] **Phase 4: NPE Training + ADVI Benchmark + Ablation** - NPE for ρ_true and Δρ at >100x ADVI speedup, summary-statistic ablation gated on per-parameter RMSE
- [ ] **Phase 5: Validation Bundle (SBC + Amortized BF + OOD)** - The publishable trifecta off one simulate→infer harness: calibration proof, amortized log-BF, honest misspecification flag
- [ ] **Phase 6: Reproducible Demo + Go/No-Go Memo** - Seeded end-to-end `demo.jl`, decoupling proof, 2-3 page Go/No-Go decision memo
- [ ] **Phase 7: Productionization (conditional on Go)** - Integrate amortized inference into `src/` behind a coexisting backend contract with user-definable `num_patches`

## Phase Details

### Phase 1: Environment + Smoke Gate
**Goal**: A reproducible, isolated spike environment exists and the pre-1.0 NeuralEstimators + Flux stack is proven to work on this machine — the highest-risk unknown is retired before any other investment
**Depends on**: Nothing (first phase)
**Requirements**: ENV-01, ENV-02, ENV-03, ENV-04
**Success Criteria** (what must be TRUE):
  1. `spike/` activates its own `Project.toml` with the main package available read-only via `Pkg.develop`; main `Project.toml`/`Manifest.toml` and `src/` are untouched
  2. `spike/00_smoke.jl` trains a `PosteriorEstimator` (NormalisingFlow) on a 1-parameter Gaussian and runs `sampleposterior` green, CPU-only with no forced CUDA import
  3. `spike/Manifest.toml` is pinned and committed with exact NeuralEstimators/Flux versions as a reproducibility artifact; the green smoke test is the hard gate for all later phases
  4. A stack-decision note records NeuralEstimators.jl as default and BayesFlow (PythonCall) as fallback only
**Plans**: 4 plans
- [x] 01-01-PLAN.md - Reconcile the dirty-root precondition; record an agreed frozen baseline ref (ENV-01)
- [x] 01-02-PLAN.md - Minimal isolated spike env + green CPU-only NeuralEstimators/Flux NPE smoke gate (ENV-02)
- [x] 01-03-PLAN.md - Pin + commit Manifest, pin Julia version, write stack-decision note (ENV-03, ENV-04)
- [x] 01-04-PLAN.md - Pkg.develop coupling (include fallback) + decoupling proof vs baseline (ENV-01)

### Phase 2: Forward Simulator + Summary Contract
**Goal**: A physics forward-simulator emits real `MultiChannelImage` pairs from θ that the *unmodified* existing summary functions ingest, with a prior provably consistent with the Turing model so every downstream comparison stays valid
**Depends on**: Phase 1
**Requirements**: SIM-01, SIM-02, SIM-03, SIM-04
**Success Criteria** (what must be TRUE):
  1. `simulate_pair(θ) → MultiChannelImage` generates a 2-channel pair via correlated densities → Bernoulli thinning → PSF → 2×2 spillover → autofluorescence → sub-pixel shift → noise from θ = (ρ_true, spillover, autofluorescence, label efficiency, shift dx/dy, noise)
  2. `sample_prior()` and π(θ) mirror the existing Turing `@model` ranges (μ/ν/σ/τ), documented in `spike/NOTES.md`, so the ADVI benchmark and BF validation share one generative prior
  3. Output verifies as a valid `MultiChannelImage` so the existing `correlation()`/`patch()`/`_prepare_data()` apply unchanged (no `src/` edits)
  4. Plausibility plots confirm expected behavior: ρ_true ↑ → patch correlation ↑; spillover and sub-pixel shift visibly affect the pair
**Plans**: TBD

### Phase 3: Training-Data Pipeline
**Goal**: A reproducible, resumable generator turns the prior + simulator + reused summary functions into cached standardized training vectors, with leak-free split discipline baked into the loader structurally
**Depends on**: Phase 2
**Requirements**: DATA-01, DATA-02, DATA-03
**Success Criteria** (what must be TRUE):
  1. The generator maps θ~π → `simulate_pair` → existing summary functions → a fixed-dimension (8×8 = 64) standardized summary vector
  2. 50k–200k pairs write and reload from a version/hash-guarded JLD2 cache (BayesInteractomics pattern) under a Random123 seed; a changed simulator/prior auto-invalidates stale data
  3. Standardization statistics are fit on training folds only and applied to held-out folds — no standardization or parameter leakage across splits — enforced in the loader, not by per-script convention
  4. Leak-free k-fold cross-validation (e.g. 5-fold) is wired so reported NPE metrics are cross-validated, with a separate reserved set of ≥20 real-pipeline stacks held out solely for the ADVI benchmark
**Plans**: TBD

### Phase 4: NPE Training + ADVI Benchmark + Ablation
**Goal**: A trained NPE infers both ρ_true and Δρ in a single forward pass, demonstrably >100x faster than per-dataset ADVI at comparable accuracy, with summary-statistic sufficiency measured rather than assumed
**Depends on**: Phase 3
**Requirements**: NPE-01, NPE-02, NPE-03, ABL-01, ABL-02
**Success Criteria** (what must be TRUE):
  1. A `PosteriorEstimator` (summary net → NormalisingFlow) is trained to infer both ρ_true and Δρ and returns posteriors in a single forward pass for ≥20 held-out stacks
  2. NPE is benchmarked against real `colocalization()` ADVI on 20–30 stacks for RMSE and interval width, with metrics reported under leak-free k-fold cross-validation (per DATA-03)
  3. Measured NPE wall-clock is **>100x faster** than per-dataset ADVI **at comparable RMSE** (milliseconds vs minutes) — speedup is reported paired with accuracy, never alone
  4. The ablation scores **per-parameter RMSE** for minimal patch-correlation vs augmented (+ Manders/median/IQR moments) as a gating sufficiency diagnostic; a SBC-pass-but-high-RMSE outcome is treated as an insufficiency signal
  5. The chosen summary is justified by the ablation result, with its interaction with OOD detectability explicitly noted (couples to Phase 5)
**Plans**: TBD

### Phase 5: Validation Bundle (SBC + Amortized BF + OOD)
**Goal**: The publishable trifecta — calibration proof, amortized Bayes factor, and honest misspecification flag — built as sibling plans off one shared θ*~π→simulate→infer harness over the trained nets; the highest-scrutiny phase for scientific honesty
**Depends on**: Phase 4
**Requirements**: SBC-01, SBC-02, SBC-03, SBC-04, BF-01, BF-02, OOD-01, OOD-02
**Success Criteria** (what must be TRUE):
  1. SBC produces per-parameter rank histograms over M (≈2000) draws with KS/χ² uniformity, a coverage curve, and ECE/MCE via the `_bin_calibration` pattern with a traffic-light verdict
  2. M and the pass/fail threshold are **pre-registered before running**, and the reported SBC number comes from a **fresh, never-tuned-against, independently-seeded held-out run** (guards "tune until calibrated" against data-snooping); calibration is reported "under the simulator" and explicitly paired with the OOD result
  3. A `RatioEstimator`/Evidence-Network (l-POP loss) produces an amortized log-BF in one forward pass and reproduces `compute_BayesFactor()` in the well-specified regime, on identical Δρ and prior, without quadgk/KDE/shuffle
  4. The OOD flag fires on misspecified inputs and stays quiet in-distribution, validated as a controlled ROC experiment over a misspecification grid with **positive AND summary-orthogonal negative controls** plus a posterior-predictive channel; the structural blind spot (discrepancies orthogonal to the fixed summary) is measured and named, not hidden
  5. All preprocessing (standardization, OOD covariance/flow) is frozen from the training split only; misspecified test images are fully external
**Plans**: TBD

### Phase 6: Reproducible Demo + Go/No-Go Memo
**Goal**: The spike closes with a falsifiable, reproducible verdict — a single seeded script chains every layer and a memo states the metrics and a concrete full-build-out decision
**Depends on**: Phase 5
**Requirements**: DEMO-01, DEMO-02, DEMO-03
**Success Criteria** (what must be TRUE):
  1. `spike/demo.jl` chains the full pipeline reproducibly from a fixed Random123 seed, CPU-only with the pinned Manifest, and tabulates the success criteria
  2. The main repo and both manuscript pipelines are demonstrably untouched — `git status` on `src/`, `bayes.jl`, `colocalization.jl` is clean (decoupling proof)
  3. A 2–3-page Go/No-Go memo reports the metrics, frames SBC as "calibrated under the simulator" paired with the OOD result, and states a concrete full-build-out decision (hierarchy/3D/multi-channel)
**Plans**: TBD

### Phase 7: Productionization (conditional on Go)
**Goal**: On a Go decision from Phase 6, validated amortized inference is promoted into `src/` as a shipped feature coexisting with the existing Turing path — the only phase that edits `src/`
**Depends on**: Phase 6 (a Go decision in the memo); does not start on a No-Go
**Requirements**: PROD-01, PROD-02
**Success Criteria** (what must be TRUE):
  1. `colocalization_amortized()` exists in `src/`, coexisting with the existing Turing/ADVI path behind a shared `_prepare_data` input + `CoLocResult` output contract, so `compute_BayesFactor()` keeps working unchanged
  2. `num_patches` (patch grid) is **user-definable** in the productionized API via an estimator registry keyed by grid, rather than retraining the network ad hoc (8×8 was the spike default)
**Plans**: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3 → 4 → 5 → 6 → 7

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Environment + Smoke Gate | 4/4 | Complete   | 2026-06-26 |
| 2. Forward Simulator + Summary Contract | 0/TBD | Not started | - |
| 3. Training-Data Pipeline | 0/TBD | Not started | - |
| 4. NPE Training + ADVI Benchmark + Ablation | 0/TBD | Not started | - |
| 5. Validation Bundle (SBC + BF + OOD) | 0/TBD | Not started | - |
| 6. Reproducible Demo + Go/No-Go Memo | 0/TBD | Not started | - |
| 7. Productionization (conditional on Go) | 0/TBD | Not started | - |
