# Requirements: ProteinCoLoc v2.0 — AmortizedColoc

**Defined:** 2026-06-26
**Core Value:** A trained NPE/NRE produces calibrated, amortized colocalization inference (posterior + Bayes factor) in a single forward pass, >100× faster than per-dataset ADVI, with a demonstrated SBC/coverage calibration proof and an honest OOD flag.

## v1 Requirements

The v2.0 milestone: a decoupled spike (ENV–DEMO) that gates a conditional productionization (PROD). Each maps to a roadmap phase.

### Environment (decoupled spike scaffold)

- [ ] **ENV-01**: A `spike/` directory exists with its own `Project.toml`, isolated from the main package (`Pkg.activate`), with the main package available read-only via `Pkg.develop` — main `Project.toml`/`Manifest.toml` and `src/` untouched
- [x] **ENV-02**: A <30-line smoke test (`spike/00_smoke.jl`) trains a `PosteriorEstimator` (NormalisingFlow) on a 1-parameter Gaussian and runs `sampleposterior`, proving the NeuralEstimators + Flux backend works on this machine
- [x] **ENV-03**: The resolved `spike/Manifest.toml` is pinned and committed (exact NeuralEstimators/Flux versions) as a reproducibility artifact; smoke test is a hard gate before further investment
- [x] **ENV-04**: A stack-decision note records NeuralEstimators.jl as default and BayesFlow (PythonCall) as fallback only

### Simulator (forward physics model)

- [ ] **SIM-01**: `simulate_pair(θ) → MultiChannelImage` generates a 2-channel image pair from θ = (ρ_true, spillover, autofluorescence, label efficiency, sub-pixel shift dx/dy, noise level) via bivariate correlated densities → Bernoulli thinning → PSF convolution → 2×2 spillover mixing → autofluorescence → sub-pixel shift → Poisson/Gaussian noise
- [ ] **SIM-02**: The prior π(θ) and `sample_prior()` are consistent with the existing Turing `@model` prior ranges (μ/ν/σ/τ), documented in `spike/NOTES.md`
- [ ] **SIM-03**: Output verifies as a valid `MultiChannelImage` so the existing `correlation()`/`patch()` functions apply unchanged
- [ ] **SIM-04**: Plausibility plots confirm expected behavior (ρ_true ↑ → patch correlation ↑; spillover and shift visibly affect the pair)

### Training Data

- [ ] **DATA-01**: A generator maps θ~π → `simulate_pair` → existing summary functions → fixed-dimension (8×8 = 64) standardized summary vector
- [ ] **DATA-02**: 50k–200k pairs are produced with a version-guarded JLD2 cache (BayesInteractomics pattern) and a Random123-seeded loader; budget starts at 50k and scales as accuracy demands
- [ ] **DATA-03**: NPE accuracy/calibration is evaluated by **leak-free k-fold cross-validation** (e.g. 5-fold): standardization statistics are fit on training folds only and applied to the held-out fold, with no standardization or parameter leakage across folds; reported metrics are cross-validated. A separate reserved set of ≥20 real-pipeline stacks is held out solely for the ADVI speed/accuracy benchmark (NPE-02/03)

### Neural Posterior Estimation

- [ ] **NPE-01**: A `PosteriorEstimator` (summary net → NormalisingFlow) is trained to infer **both ρ_true and Δρ**
- [ ] **NPE-02**: For ≥20 held-out stacks, the NPE returns posteriors in a single forward pass, benchmarked against real `colocalization()` ADVI on 20–30 stacks for RMSE and interval width
- [ ] **NPE-03**: Measured NPE wall-clock is **>100× faster** than per-dataset ADVI at comparable RMSE (millisecond vs minutes)

### Summary-Statistic Ablation

- [ ] **ABL-01**: NPE accuracy is compared between the minimal patch-correlation summary and an augmented summary (+ Manders / median / IQR moments), scoring **per-parameter RMSE** as a gating sufficiency diagnostic (not decorative)
- [ ] **ABL-02**: The chosen summary is justified by the ablation result and its interaction with OOD detectability is noted

### Simulation-Based Calibration

- [ ] **SBC-01**: M (≈2000) draws θ*~π → simulate → L posterior draws → rank of θ* produce a rank histogram per parameter
- [ ] **SBC-02**: Uniformity is tested (KS/χ²) and a coverage curve (nominal vs empirical) computed; ECE/MCE via the `_bin_calibration` pattern with a traffic-light verdict
- [ ] **SBC-03**: M and the pass/fail threshold are **pre-registered before running**; the reported SBC is a fresh, never-tuned-against held-out run (guards "tune until calibrated" against data-snooping)
- [ ] **SBC-04**: Calibration is reported "under the simulator" and explicitly paired with the OOD result (SBC alone is not a real-data guarantee)

### Amortized Bayes Factor

- [ ] **BF-01**: A `RatioEstimator` / Evidence-Network (l-POP-style loss) produces an amortized log-BF (coloc vs null) in one forward pass, with Δρ and prior alignment identical to the KDE baseline
- [ ] **BF-02**: The amortized log-BF reproduces the existing `compute_BayesFactor()` in the well-specified regime, without quadgk/KDE/shuffle

### OOD / Misspecification Flag

- [ ] **OOD-01**: An OOD score (summary-density Mahalanobis / flow-likelihood + posterior-predictive mismatch) flags misspecified inputs and stays quiet in-distribution
- [ ] **OOD-02**: Validation is a controlled experiment (ROC over a misspecification grid) with **positive AND summary-orthogonal negative controls**; the structural blind spot (discrepancies orthogonal to the fixed summary) is measured and named, not hidden

### Demo & Decision

- [ ] **DEMO-01**: `spike/demo.jl` chains the full pipeline reproducibly from a fixed seed and tabulates success criteria
- [ ] **DEMO-02**: The main repo and both manuscript pipelines are demonstrably untouched (decoupling proof)
- [ ] **DEMO-03**: A 2–3-page Go/No-Go memo reports metrics and a concrete full-build-out decision (hierarchy/3D/multi-channel)

### Productionization (conditional on a Go decision)

- [ ] **PROD-01**: On Go, amortized inference is integrated into `src/` as `colocalization_amortized()` coexisting with the existing Turing path behind a shared input/output contract, so `compute_BayesFactor()` keeps working unchanged
- [ ] **PROD-02**: `num_patches` (patch grid) is **user-definable** in the productionized API (8×8 was the spike default), via an estimator registry rather than retraining the network ad hoc

## v2 Requirements

Deferred / optional; tracked but not gating the v2.0 milestone.

### Backend Evaluation

- **BACK-01**: Time-boxed feasibility evaluation of a Turing → RxInfer.jl backend migration (fixed-ν Gaussian scale-mixture formulation), benchmarked as an independent baseline cross-check; deferred to post-Go, never a spike dependency
- **BACK-02**: DeepSet permutation-invariant summary network as an upgrade if the MLP summary proves insufficient

## Out of Scope

| Feature | Reason |
|---------|--------|
| Full 4-level hierarchy (Pixel→Segment→Celltype→Neighborhood) | Spike proves 2-channel 2D first |
| 3D / Z-stacks | Deferred to full build-out (GPU territory) |
| Higher-order chromatic aberration | Beyond spike simulator realism |
| GUI | Not needed for the inference claim |
| Multi→2-channel generalization | Spike fixes 2 channels |
| Edits to `src/bayes.jl`, `src/colocalization.jl`, main pipeline (during spike) | Strict decoupling constraint |
| Wet-lab data, AlphaFold3/Boltz/Chai | Fully solo, simulation-only |
| Sequential SBI (SNPE/SNRE) | Breaks amortization; contradicts ms-inference claim |
| TARP / joint-coverage diagnostics | Marginal SBC clears the spike bar |

## Traceability

Mapped during roadmap creation (2026-06-26). Every v1 requirement maps to exactly one phase.

| Requirement | Phase | Status |
|-------------|-------|--------|
| ENV-01 | Phase 1 — Environment + Smoke Gate | Pending |
| ENV-02 | Phase 1 — Environment + Smoke Gate | Complete |
| ENV-03 | Phase 1 — Environment + Smoke Gate | Complete |
| ENV-04 | Phase 1 — Environment + Smoke Gate | Complete |
| SIM-01 | Phase 2 — Forward Simulator + Summary Contract | Pending |
| SIM-02 | Phase 2 — Forward Simulator + Summary Contract | Pending |
| SIM-03 | Phase 2 — Forward Simulator + Summary Contract | Pending |
| SIM-04 | Phase 2 — Forward Simulator + Summary Contract | Pending |
| DATA-01 | Phase 3 — Training-Data Pipeline | Pending |
| DATA-02 | Phase 3 — Training-Data Pipeline | Pending |
| DATA-03 | Phase 3 — Training-Data Pipeline | Pending |
| NPE-01 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Pending |
| NPE-02 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Pending |
| NPE-03 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Pending |
| ABL-01 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Pending |
| ABL-02 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Pending |
| SBC-01 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| SBC-02 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| SBC-03 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| SBC-04 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| BF-01 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| BF-02 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| OOD-01 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| OOD-02 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| DEMO-01 | Phase 6 — Reproducible Demo + Go/No-Go Memo | Pending |
| DEMO-02 | Phase 6 — Reproducible Demo + Go/No-Go Memo | Pending |
| DEMO-03 | Phase 6 — Reproducible Demo + Go/No-Go Memo | Pending |
| PROD-01 | Phase 7 — Productionization (conditional on Go) | Pending |
| PROD-02 | Phase 7 — Productionization (conditional on Go) | Pending |

**Coverage:**
- v1 requirements: 27 total
- Mapped to phases: 27 ✓
- Unmapped: 0
- v2 requirements (BACK-01, BACK-02): deferred, not gating v2.0 — intentionally unmapped

---
*Requirements defined: 2026-06-26*
*Last updated: 2026-06-26 after roadmap creation (traceability populated, 27/27 v1 mapped)*
