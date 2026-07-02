# Requirements: ProteinCoLoc v2.0 — AmortizedColoc

**Defined:** 2026-06-26
**Core Value:** A trained NPE/NRE produces calibrated, amortized colocalization inference (posterior + Bayes factor) in a single forward pass, >100× faster than per-dataset ADVI, with a demonstrated SBC/coverage calibration proof and an honest OOD flag.

## v1 Requirements

The v2.0 milestone: a decoupled spike (ENV–DEMO) that gates a conditional productionization (PROD). Each maps to a roadmap phase.

### Environment (decoupled spike scaffold)

- [x] **ENV-01**: A `spike/` directory exists with its own `Project.toml`, isolated from the main package (`Pkg.activate`), with the main package available read-only via `Pkg.develop` — main `Project.toml`/`Manifest.toml` and `src/` untouched
- [x] **ENV-02**: A <30-line smoke test (`spike/00_smoke.jl`) trains a `PosteriorEstimator` (NormalisingFlow) on a 1-parameter Gaussian and runs `sampleposterior`, proving the NeuralEstimators + Flux backend works on this machine
- [x] **ENV-03**: The resolved `spike/Manifest.toml` is pinned and committed (exact NeuralEstimators/Flux versions) as a reproducibility artifact; smoke test is a hard gate before further investment
- [x] **ENV-04**: A stack-decision note records NeuralEstimators.jl as default and BayesFlow (PythonCall) as fallback only

### Simulator (forward physics model)

- [x] **SIM-01**: `simulate_pair(θ) → MultiChannelImage` generates a 2-channel image pair from θ = (ρ_true, spillover, autofluorescence, label efficiency, sub-pixel shift dx/dy, noise level) via bivariate correlated densities → Bernoulli thinning → PSF convolution → 2×2 spillover mixing → autofluorescence → sub-pixel shift → Poisson/Gaussian noise
- [x] **SIM-02**: The prior π(θ) and `sample_prior()` are consistent with the existing Turing `@model` prior ranges (μ/ν/σ/τ), documented in `spike/NOTES.md`
- [x] **SIM-03**: Output verifies as a valid `MultiChannelImage` so the existing `correlation()`/`patch()` functions apply unchanged
- [x] **SIM-04**: Plausibility plots confirm expected behavior (ρ_true ↑ → patch correlation ↑; spillover and shift visibly affect the pair)

### Training Data

- [x] **DATA-01**: A generator maps θ~π → `simulate_pair` → existing summary functions → fixed-dimension (8×8 = 64) standardized summary vector
- [x] **DATA-02**: 50k–200k pairs are produced with a version-guarded JLD2 cache (BayesInteractomics pattern) and a Random123-seeded loader; budget starts at 50k and scales as accuracy demands
- [x] **DATA-03**: NPE accuracy/calibration is evaluated by **leak-free k-fold cross-validation** (e.g. 5-fold): standardization statistics are fit on training folds only and applied to the held-out fold, with no standardization or parameter leakage across folds; reported metrics are cross-validated. A separate reserved set of ≥20 real-pipeline stacks is held out solely for the ADVI speed/accuracy benchmark (NPE-02/03)

### Neural Posterior Estimation

- [x] **NPE-01**: A `PosteriorEstimator` (summary net → NormalisingFlow) is trained to infer **both ρ_true and Δρ**
- [x] **NPE-02**: For ≥20 held-out stacks, the NPE returns posteriors in a single forward pass, benchmarked against real `colocalization()` ADVI on 20–30 stacks for RMSE and interval width
- [x] **NPE-03**: Measured NPE wall-clock is **>100× faster** than per-dataset ADVI at comparable RMSE (millisecond vs minutes) — verified with an honestly-disclosed caveat: realized ADVI baseline is ~0.5s/pair (seconds, not minutes); >100× (median 325×) is scored on the forward-pass clock, ~16× on the full posterior-sample workload. Restate in Phase 6 Go/No-Go memo.

### Summary-Statistic Ablation

- [x] **ABL-01**: NPE accuracy is compared between the minimal patch-correlation summary and an augmented summary (+ Manders / median / IQR moments), scoring **per-parameter RMSE** as a gating sufficiency diagnostic (not decorative)
- [x] **ABL-02**: The chosen summary is justified by the ablation result and its interaction with OOD detectability is noted

### Simulation-Based Calibration

- [x] **SBC-01**: M (≈2000) draws θ*~π → simulate → L posterior draws → rank of θ* produce a rank histogram per parameter
- [x] **SBC-02**: Uniformity is tested (KS/χ²) and a coverage curve (nominal vs empirical) computed; ECE/MCE via the `_bin_calibration` pattern with a traffic-light verdict
- [x] **SBC-03**: M and the pass/fail threshold are **pre-registered before running**; the reported SBC is a fresh, never-tuned-against held-out run (guards "tune until calibrated" against data-snooping)
- [x] **SBC-04**: Calibration is reported "under the simulator" and explicitly paired with the OOD result (SBC alone is not a real-data guarantee)

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

### Cross-Method Comparator Harness (Phase 9 — feature expansion)

Spike-local harness that positions v2.0 as "knows when the classics are wrong." All new code under `spike/comparator/`; no `src/` or root-project edits; `spike/data/encode.jl` left byte-identical (hash-guarded).

- [ ] **CMP-01**: Classical battery (Costes-p, Manders M1/M2, Pearson, Spearman) runs on one shared `MultiChannelImage` input set and emits a per-method comparison table (SC1; D-01,D-04,D-06,D-09)
- [ ] **CMP-02**: Costes randomization significance p-value implemented new, spike-local, seeded (block-scramble null; p = (#{r ≥ r_obs}+1)/(n+1)); `manders(mci)` proven equal to the M1/M2 `encode_aug` computes for the same mci (SC1; D-04,D-05)
- [ ] **CMP-03**: Shared inputs are simulator-generated over a seeded θ grid with ground-truth regime labels; harness signature accepts any `Vector{MultiChannelImage}` (SC1; D-02,D-03)
- [ ] **CMP-04**: Divergence / traffic-light "knows when classics are wrong" column with **pre-declared** thresholds (SC1 positioning; D-10,D-14)
- [ ] **CMP-05**: Table persisted as tidy `DataFrame` + CSV + JLD2, content-addressed / seeded like the Phase-3 cache, over the comparator's OWN hash inputs (SC1; D-09,D-12)
- [ ] **CMP-06**: Tapqir bridge reproduces a published Tapqir tutorial example within tolerance **OR** skips-with-flag when the isolated Python env is absent (SC2; D-07,D-08)
- [ ] **CMP-07**: Reuse the BayesInteractomics comparator/audit pattern (`_bin_calibration`/`CalibrationResult`/traffic-light; `_ks_test_uniform`/report style) for table assembly + audit summary (SC3; D-11)
- [ ] **CMP-08**: Single seeded entry point `run_comparator.jl` (Random123 Philox); bit-reproducible across runs and thread counts (SC3; D-12)
- [ ] **CMP-09**: Spike test gate asserts determinism, finiteness, Manders-equality, traffic-light-by-oracle, and Tapqir anchor-or-skip (SC1,SC2; D-13)

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
| Tapqir as a live comparator column on our images | Wrong data regime (CoSMoS single-molecule spots vs diffuse 2-channel fields); anchor-only (D-07) |

## Traceability

Mapped during roadmap creation (2026-06-26). Every v1 requirement maps to exactly one phase. Feature-expansion Phase 9 (CMP-*) added 2026-07-02 during phase planning.

| Requirement | Phase | Status |
|-------------|-------|--------|
| ENV-01 | Phase 1 — Environment + Smoke Gate | Complete |
| ENV-02 | Phase 1 — Environment + Smoke Gate | Complete |
| ENV-03 | Phase 1 — Environment + Smoke Gate | Complete |
| ENV-04 | Phase 1 — Environment + Smoke Gate | Complete |
| SIM-01 | Phase 2 — Forward Simulator + Summary Contract | Complete |
| SIM-02 | Phase 2 — Forward Simulator + Summary Contract | Complete |
| SIM-03 | Phase 2 — Forward Simulator + Summary Contract | Complete |
| SIM-04 | Phase 2 — Forward Simulator + Summary Contract | Complete |
| DATA-01 | Phase 3 — Training-Data Pipeline | Complete |
| DATA-02 | Phase 3 — Training-Data Pipeline | Complete |
| DATA-03 | Phase 3 — Training-Data Pipeline | Complete |
| NPE-01 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Complete |
| NPE-02 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Complete |
| NPE-03 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Complete |
| ABL-01 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Complete |
| ABL-02 | Phase 4 — NPE Training + ADVI Benchmark + Ablation | Complete |
| SBC-01 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Complete |
| SBC-02 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Complete |
| SBC-03 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Complete |
| SBC-04 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Complete |
| BF-01 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| BF-02 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| OOD-01 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| OOD-02 | Phase 5 — Validation Bundle (SBC + BF + OOD) | Pending |
| DEMO-01 | Phase 6 — Reproducible Demo + Go/No-Go Memo | Pending |
| DEMO-02 | Phase 6 — Reproducible Demo + Go/No-Go Memo | Pending |
| DEMO-03 | Phase 6 — Reproducible Demo + Go/No-Go Memo | Pending |
| PROD-01 | Phase 7 — Productionization (conditional on Go) | Pending |
| PROD-02 | Phase 7 — Productionization (conditional on Go) | Pending |
| CMP-01 | Phase 9 — Cross-Method Comparator Harness | Pending |
| CMP-02 | Phase 9 — Cross-Method Comparator Harness | Pending |
| CMP-03 | Phase 9 — Cross-Method Comparator Harness | Pending |
| CMP-04 | Phase 9 — Cross-Method Comparator Harness | Pending |
| CMP-05 | Phase 9 — Cross-Method Comparator Harness | Pending |
| CMP-06 | Phase 9 — Cross-Method Comparator Harness | Pending |
| CMP-07 | Phase 9 — Cross-Method Comparator Harness | Pending |
| CMP-08 | Phase 9 — Cross-Method Comparator Harness | Pending |
| CMP-09 | Phase 9 — Cross-Method Comparator Harness | Pending |

**Coverage:**
- v1 requirements: 27 total
- Mapped to phases: 27 ✓
- Unmapped: 0
- Feature-expansion requirements: CMP-01..CMP-09 (Phase 9) mapped ✓; other Phase 8/10–16 IDs still TBD
- v2 requirements (BACK-01, BACK-02): deferred, not gating v2.0 — intentionally unmapped

---
*Requirements defined: 2026-06-26*
*Last updated: 2026-07-02 — added Phase 9 CMP-01..CMP-09 (cross-method comparator harness) during phase planning*
