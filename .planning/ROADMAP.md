# Roadmap: ProteinCoLoc v2.0 — AmortizedColoc

## Overview

v2.0 proves and (conditionally) ships amortized simulation-based inference for colocalization. The work is a strictly decoupled spike (`spike/`, own `Project.toml`, no `src/` edits) built bottom-up as a 7-layer horizontal pipeline: a NeuralEstimators+Flux smoke gate de-risks the pre-1.0 API before any investment; a physics forward-simulator emits real `MultiChannelImage`s that the *unmodified* existing summary functions ingest; a version-guarded training cache feeds an NPE benchmarked at >100x ADVI speedup; a shared θ*~π→simulate→infer harness fans out into the SBC/BF/OOD validation trifecta; a reproducible demo closes with a Go/No-Go memo. Only on a Go decision does the final phase touch `src/`, integrating amortized inference behind a coexisting backend contract. Every phase is an independently runnable checkpoint, resumable across pauses.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Environment + Smoke Gate** - Isolated spike env, pinned Manifest, green NeuralEstimators/Flux toy NPE — the hard gate before all downstream investment (completed 2026-06-26)
- [x] **Phase 2: Forward Simulator + Summary Contract** - `simulate_pair(θ) → MultiChannelImage` consumed unchanged by the existing summary functions, prior-consistent with the Turing model (completed 2026-06-26)
- [x] **Phase 3: Training-Data Pipeline** - θ~π → simulate → summary → standardized fixed-dim vectors, version-guarded JLD2 cache with leak-free split discipline (completed 2026-06-30)
- [x] **Phase 4: NPE Training + ADVI Benchmark + Ablation** - NPE for ρ_true and Δρ at >100x ADVI speedup, summary-statistic ablation gated on per-parameter RMSE (completed 2026-07-01)
- [x] **Phase 5: Validation Bundle (SBC + Amortized BF + OOD)** - The publishable trifecta off one simulate→infer harness: calibration proof, amortized log-BF, honest misspecification flag
 (completed 2026-07-02)
- [x] **Phase 6: Reproducible Demo + Go/No-Go Memo** - Seeded end-to-end `demo.jl`, decoupling proof, 2-3 page Go/No-Go decision memo
 (completed 2026-07-03)
- [ ] **Phase 7: Productionization (conditional on Go)** - Integrate amortized inference into `src/` behind a coexisting backend contract with user-definable `num_patches`

**v2.0 feature expansion (Phases 8–16) — all downstream of a Phase-6 Go:**
- [x] **Phase 8: External Physical Ground-Truth Corpus** - Non-circular validation anchors (physical 100%-coloc + segregated constructs); CBS ingested as labelled-simulated; versioned data contract (completed 2026-07-21)
- [x] **Phase 9: Cross-Method Comparator Harness** - Costes/Manders/Pearson/Spearman + Tapqir bridge on shared inputs; "knows when the classics are wrong" — parallelizable now
- [x] **Phase 10: Manuscript Skeleton + Related-Work Positioning** - Compiling Typst skeleton, explicit Tapqir/Costes/Manders delta, figure specs — parallelizable now
 (completed 2026-07-02)
- [ ] **Phase 11: Registration + Chromatic Uncertainty as Latent** - Promote dx/dy (+ chromatic warp) to inferred θ; posterior widens honestly under registration uncertainty
- [ ] **Phase 12: Spatial Colocalization Map (GP/CAR)** - Lattice prior over the correlation grid → amortized per-region Δρ map + uncertainty (descope-to-v2.1 candidate)
- [ ] **Phase 13: Three-Hypothesis Amortized Bayes Factor** - Evidence network extended to coloc/random/exclusion; replaces KDE+quadgk BF
- [ ] **Phase 14: Decision + Abstention Layer** - {coloc/not/ABSTAIN} at controlled Bayesian FDR; abstains on OOD/disagreement/ambiguity
- [ ] **Phase 15: Calibration Operating Envelope + CI Gate** - Adversarial nuisance sweep → domain-of-applicability; OOD-before-break; SBC regression gate in CI
- [ ] **Phase 16: External Validation + Manuscript Assembly** - Blind eval vs physical corpus + comparators; figures assembled; one-command repro + Zenodo

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
**Plans**: 4 plans (4 waves — serial: env → simulator → calibration → validation)
- [x] 02-01-PLAN.md - Extend spike env (StatsBase/Images/ImageFiltering/CairoMakie), re-freeze Manifest, establish the include() contract boundary, prove SIM-03 on a synthetic image (SIM-03)
- [x] 02-02-PLAN.md - `simulate_pair` 7-stage forward physics pipeline; SIM-03 on real simulator output (SIM-01, SIM-03)
- [x] 02-03-PLAN.md - Induced-μ calibration (sweep → monotone ĝ), `sample_prior` consistent with the Turing μ-prior, NOTES §3 evidence (SIM-02)
- [x] 02-04-PLAN.md - SIM-04 quantitative plausibility gate (monotonicity + paired perturbation) + CairoMakie figures (SIM-04)

### Phase 3: Training-Data Pipeline
**Goal**: A reproducible, resumable generator turns the prior + simulator + reused summary functions into cached standardized training vectors, with leak-free split discipline baked into the loader structurally
**Depends on**: Phase 2
**Requirements**: DATA-01, DATA-02, DATA-03
**Success Criteria** (what must be TRUE):
  1. The generator maps θ~π → `simulate_pair` → existing summary functions → a fixed-dimension (8×8 = 64) standardized summary vector
  2. 50k–200k pairs write and reload from a version/hash-guarded JLD2 cache (BayesInteractomics pattern) under a Random123 seed; a changed simulator/prior auto-invalidates stale data
  3. Standardization statistics are fit on training folds only and applied to held-out folds — no standardization or parameter leakage across splits — enforced in the loader, not by per-script convention
  4. Leak-free k-fold cross-validation (e.g. 5-fold) is wired so reported NPE metrics are cross-validated, with a separate reserved set of ≥20 real-pipeline stacks held out solely for the ADVI benchmark
**Plans**: 5 plans (5 waves — serial: env scaffold → generation core → cache/hash → leak-free loader → real-run gate)
- [x] 03-01-PLAN.md — Wave-0 env: add JLD2/Random123, re-freeze Manifest, extend resolve-risk gate, scaffold test_data_pipeline.jl (DATA-02)
- [x] 03-02-PLAN.md — DATA-01 generation core: 128-dim D-01 encode + D-02 augmented variant, Random123 keyed seeding, order/thread-independent generator (DATA-01)
- [x] 03-03-PLAN.md — DATA-02 cache: sharded atomic JLD2, content-hash version guard, resume-by-skip, reserved holdout, N-parameterized scale-up (DATA-02)
- [x] 03-04-PLAN.md — DATA-03 loader: leak-free k=5 fold standardization (fit-on-train-only, no global path), mask bypass, reserved-holdout exclusion (DATA-03)
- [x] 03-05-PLAN.md — Phase gate: real ≥50k generation run, cross-process thread-repro, full-suite green (DATA-02)

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
**Plans**: 7 plans (5 waves — W1 env scaffolds ∥ → W2 NPE train ∥ ADVI artifact → W3 benchmark → W4 ablation → W5 scaling/thread)
- [x] 04-01-PLAN.md — Wave-0 spike env: add BenchmarkTools, re-freeze Manifest + resolve-risk gate, test_npe.jl SC1..SC5 scaffold + pre-registered consts + keyed holdout-repro gate (NPE-01)
- [x] 04-02-PLAN.md — Isolated spike/baseline/ env (Turing/AdvancedVI + simulator deps) + read-only colocalization @model lift (NPE-02)
- [x] 04-03-PLAN.md — NPE architecture + CPU-only fixed-data training + inference surface (ρ̂/Δρ MC-diff/interval/ghat); SC1 (NPE-01)
- [x] 04-04-PLAN.md — ADVI port to modern AdvancedVI on raw re-simulated paired holdout → integrity-checked advi_artifact.jld2 (NPE-02, NPE-03)
- [x] 04-05-PLAN.md — ρ-space RMSE/interval benchmark + BenchmarkTools >100× speedup gate; SC2/SC3 (NPE-02, NPE-03)
- [x] 04-06-PLAN.md — :min vs :aug k=5 CV ablation + pre-registered decision rule + OOD-noted justification; SC4/SC5 (ABL-01, ABL-02)
- [x] 04-07-PLAN.md — Scaling curves over N/imsize + CPU thread-count sweep characterization; SC3-scaling (NPE-03)

### Phase 5: Validation Bundle (SBC + Amortized BF + OOD)
**Goal**: The publishable trifecta — calibration proof, amortized Bayes factor, and honest misspecification flag — built as sibling plans off one shared θ*~π→simulate→infer harness over the trained nets; the highest-scrutiny phase for scientific honesty
**Depends on**: Phase 4
**Requirements**: SBC-01, SBC-02, SBC-03, SBC-04, BF-01, BF-02, OOD-01, OOD-02
**Success Criteria** (what must be TRUE):
  1. SBC produces per-parameter rank histograms over M (≈2000) draws with KS/χ² uniformity, a coverage curve, and ECE/MCE via the `_bin_calibration` pattern with a traffic-light verdict
  2. M and the pass/fail threshold are **pre-registered before running**, and the reported SBC number comes from a **fresh, never-tuned-against, independently-seeded held-out run** (guards "tune until calibrated" against data-snooping); calibration is reported "under the simulator" and explicitly paired with the OOD result
  3. A `RatioEstimator`/Evidence-Network (NRE-as-model-comparison over a binary model index m∈{0,1}; the l-POP loss does NOT exist in NeuralEstimators v0.2.1 and is honestly re-labeled as this standard-NRE reduction, l-POP a documented fallback only) produces an amortized log-BF in one forward pass and reproduces `compute_BayesFactor()` in the well-specified regime, on identical Δρ and prior, without quadgk/KDE/shuffle
  4. The OOD flag fires on misspecified inputs and stays quiet in-distribution, validated as a controlled ROC experiment over a misspecification grid with **positive AND summary-orthogonal negative controls** plus a posterior-predictive channel; the structural blind spot (discrepancies orthogonal to the fixed summary) is measured and named, not hidden
  5. All preprocessing (standardization, OOD covariance/flow) is frozen from the training split only; misspecified test images are fully external
**Plans**: 4 plans (3 waves — W1 shared harness + pre-registration + SBC; W2 amortized BF ∥ OOD; W3 reported trifecta runs)
- [x] 05-01-PLAN.md — Shared θ*~π→simulate→infer harness + pre-registered consts + SBC machinery (rank histograms 7θ+Δρ, KS/χ², coverage, ECE/MCE traffic-light) (SBC-01..04)
- [x] 05-02-PLAN.md — Amortized Bayes factor: RatioEstimator model-index net + one-pass log-BF reproducing compute_BayesFactor over a held-out Δρ sweep (BF-01, BF-02)
- [x] 05-03-PLAN.md — OOD/misspecification flag: Mahalanobis + PP channels, ROC over the 4-family positive grid + summary-orthogonal negative control (named blind spot) (OOD-01, OOD-02)
- [x] 05-04-PLAN.md — Reported trifecta runs at locked pre-registered consts on the fresh reserved stream: run_sbc.jl (M=2000/L=999) + run_bf.jl (full Δρ sweep) + run_ood.jl (full 4-family grid), each gating on its thresholds (SBC-01..04, BF-01/02, OOD-01/02)

### Phase 6: Reproducible Demo + Go/No-Go Memo
**Goal**: The spike closes with a falsifiable, reproducible verdict — a single seeded script chains every layer and a memo states the metrics and a concrete full-build-out decision
**Depends on**: Phase 5
**Requirements**: DEMO-01, DEMO-02, DEMO-03
**Success Criteria** (what must be TRUE):
  1. `spike/demo.jl` chains the full pipeline reproducibly from a fixed Random123 seed, CPU-only with the pinned Manifest, and tabulates the success criteria
  2. The main repo and both manuscript pipelines are demonstrably untouched — `git status` on `src/`, `bayes.jl`, `colocalization.jl` is clean (decoupling proof)
  3. A 2–3-page Go/No-Go memo reports the metrics, frames SBC as "calibrated under the simulator" paired with the OOD result, and states a concrete full-build-out decision (hierarchy/3D/multi-channel)
**Plans**: 2 plans (2 waves — W1 demo.jl two-tier runner + decoupling proof ; W2 Go/No-Go memo + memo-gate)
- [x] 06-01-PLAN.md — Two-tier seeded demo.jl: fixture-scale NPE+BF+OOD chain proof (per-layer twin-run reproducibility) + git-status decoupling assertion + --full reported-gate dispatch (DEMO-01, DEMO-02)
- [x] 06-02-PLAN.md — 2-3 page Go/No-Go memo (Clean Go, both number sets, falsification, ship-gate, Phase 8-16 DAG) + demo.jl memo-content gate (DEMO-03)

### Phase 7: Productionization (conditional on Go)
**Goal**: On a Go decision from Phase 6, validated amortized inference is promoted into `src/` as a shipped feature coexisting with the existing Turing path — the only phase that edits `src/`
**Depends on**: Phase 6 (a Go decision in the memo); does not start on a No-Go
**Requirements**: PROD-01, PROD-02
**Success Criteria** (what must be TRUE):
  1. `colocalization_amortized()` exists in `src/`, coexisting with the existing Turing/ADVI path behind a shared `_prepare_data` input + `CoLocResult` output contract, so `compute_BayesFactor()` keeps working unchanged
  2. `num_patches` (patch grid) is **user-definable** in the productionized API via an estimator registry keyed by grid, rather than retraining the network ad hoc (8×8 was the spike default)
**Plans**: 11 plans (9 waves — W0 dep-surgery+resolve-gate+types ; W1 grid-param+registry ; W2 read-surfaces ; W3 training+persistence ; W4 gate-harness+GPU-smoke ; W5 8×8 ; W6 4×4 ∥ 16×16 ∥ sub-tile-map ; W7 32×32-conditional ; W8 public-API+release)
- [x] 07-00-PLAN.md — Dep surgery (Turing→ext, CUDA weakdep) + co-resolution HARD gate + version 2.0.0 + D-02 type hierarchy (PROD-01)
- [x] 07-01-PLAN.md — Grid-parametrize summary/encoder/loader/datagen + grid-keyed registry skeleton + train_and_register (PROD-02)
- [x] 07-02-PLAN.md — Amortized read surfaces: NPE infer + NRE bf (+non-clamped baseline) + OOD (+PP re-enable) (PROD-01)
- [x] 07-03-PLAN.md — NPE+NRE training GPU-plumbed + CPU-resident Flux.state persistence (PROD-01, PROD-02)
- [x] 07-04-PLAN.md — Per-grid CPU ship-gate harness + fresh-pre-reg gate_consts template + GPU-train smoke (PROD-02)
- [x] 07-05-PLAN.md — 8×8 pipeline: train → fresh CPU SBC/BF/OOD gate → register (proves gate machinery) (PROD-02)
- [x] 07-06-PLAN.md — 4×4 pipeline (coarse-robust): train → fresh CPU gate (PROD-02)
- [x] 07-07-PLAN.md — 16×16 pipeline (≥512² data, min-image-size caveat): train → fresh CPU gate (PROD-02) — honest FAIL recorded; grid 16 NOT registry-eligible by default
- [x] 07-08-PLAN.md — 8×8 sub-tile windowed coarse local colocalization map (inherits 8×8 gate; Phase-12 point clean) (PROD-02) — `local_coloc_map`/`LocalColocMap` + in-process `_ensure_grid_registered` bridge; full suite green
- [ ] 07-09-PLAN.md — 32×32 CONDITIONAL pipeline (≥1024² data) + user ship-with-caveat vs cap decision; 64×64 DROPPED (PROD-02)
- [ ] 07-10-PLAN.md — Public colocalization_amortized API + content-hashed Artifacts + register only gate-PASSED grids + integration test (PROD-01, PROD-02)

## Progress

**Execution Order:**
Spike (Phases 1–7) executes in numeric order: 1 → 2 → 3 → 4 → 5 → 6 → 7.
v2.0 feature expansion (8–16) is a DAG, not a chain, and all of it is downstream of a Phase-6 **Go**:
- Wave A {8, 9, 10} — no code dependency; parallelizable **now**, alongside Phases 5–7
- Wave B: after Phase 7 (the AP1 API spine), run 11 ∥ 13 as parallel workstreams; 12 follows 11 (shared estimator/training code — serialized to avoid a merge collision)
- Wave C {14, 15} — converge the features
- Wave D {16} — external validation + manuscript assembly

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Environment + Smoke Gate | 4/4 | Complete   | 2026-06-26 |
| 2. Forward Simulator + Summary Contract | 4/4 | Complete   | 2026-06-26 |
| 3. Training-Data Pipeline | 5/5 | Complete   | 2026-06-30 |
| 4. NPE Training + ADVI Benchmark + Ablation | 7/7 | Complete   | 2026-07-01 |
| 5. Validation Bundle (SBC + BF + OOD) | 4/4 | Complete   | 2026-07-02 |
| 6. Reproducible Demo + Go/No-Go Memo | 2/2 | Complete   | 2026-07-03 |
| 7. Productionization (conditional on Go) | 7/11 | In Progress|  |
| 8. External Physical Ground-Truth Corpus | 4/5 | In Progress|  |
| 9. Cross-Method Comparator Harness | 6/6 | Complete   | 2026-07-02 |
| 10. Manuscript Skeleton + Related-Work Positioning | 5/5 | Complete   | 2026-07-02 |
| 11. Registration + Chromatic Uncertainty as Latent | 6/11 | In Progress|  |
| 12. Spatial Colocalization Map (GP/CAR) | 0/TBD | Not started | - |
| 13. Three-Hypothesis Amortized Bayes Factor | 2/16 | In Progress|  |
| 14. Decision + Abstention Layer | 0/TBD | Not started | - |
| 15. Calibration Operating Envelope + CI Gate | 0/TBD | Not started | - |
| 16. External Validation + Manuscript Assembly | 0/TBD | Not started | - |

### Phase 8: External Physical Ground-Truth Corpus
**Goal**: Assemble a non-circular validation corpus whose truth does not come from the model's own simulator, so v2.0's calibration can be checked against external reality rather than mere self-consistency
**Depends on**: Nothing (parallelizable now, alongside Phases 5–7)
**Requirements**: SC1, SC2, SC3 (local — ROADMAP success criteria; Requirements = TBD in REQUIREMENTS.md)
**Success Criteria** (what must be TRUE):
  1. ≥1 physical 100%-coloc anchor (single-protein-two-channel / tandem fluorophore) and ≥1 segregated anchor (nuclear-vs-membrane) archived with provenance + content hashes
  2. The Colocalization Benchmark Source is ingested and explicitly flagged as simulated/secondary, never conflated with the physical anchors
  3. A versioned validation data contract + manifest (schema, provenance, split policy) is checked in
**Plans**: 5 plans (5 waves — W1 config+test-scaffold+fixture ; W2 hash+fetch ; W3 manifest+guards+conversion ; W4 CBS ingestion+manifest.csv ; W5 physical-anchor human-verify pinning)
- [x] 08-01-PLAN.md — Pre-declared config consts + scoped .gitignore + offline runtests.jl gate + synthetic two-channel TIFF fixture (SC1, SC3)
- [x] 08-02-PLAN.md — SHA-256 content-hash utility + fetch_verified (D-05 skip-vs-abort asymmetry, atomic commit, timeout) (SC1, SC3)
- [x] 08-03-PLAN.md — Manifest schema + D-06 tier guard + D-09 sealed-holdout guard + D-03 anchor predicate + content-addressed writer + read-only MultiChannelImage conversion (SC1, SC2, SC3)
- [x] 08-04-PLAN.md — Full CBS ingestion (simulated-secondary, seeded dev/eval split, zip-slip-safe) + committed manifest.csv + offline-skip fetch smoke (SC2, SC3)
- [x] 08-05-PLAN.md — Physical anchor pinning: positive tandem-FP + negative segregated construct, human-verified, archived sealed-holdout with provenance + SHA-256 (SC1) [autonomous:false]

**Completion note (2026-07-21)** — Phase 8 COMPLETE; offline gate `julia --project=. corpus/test/runtests.jl` 217/217, `src/` untouched, `git ls-files corpus/data` empty. Three recorded deviations carried into Phase 16:
  - **D-01 substitution (human-accepted):** no open-licensed, non-environment-quenched tandem-FP dataset exists in any public archive (every deposited tandem-FP set is a quenched autophagy reporter, disqualified by Pitfall 2). POSITIVE anchor is instead TetraSpeck 100 nm multicolor beads (RegiSTORM sample data, Zenodo `10.5281/zenodo.5509861`, CC-BY-4.0) — the same physical particle emits in both channels, so coloc is by construction AND state-independent.
  - **D-02 preference unmet:** NEGATIVE anchor is "Light My Cells" (BioImage Archive `S-BIAD1047`, CC-BY-4.0, nucleus vs mitochondria). Cross-study/cross-archive (`ANCHORS_MATCHED=false`) — imaging-condition confounds between the anchors are NOT controlled; Phase 16 must report this limitation.
  - **Hashes pending:** the bootstrap fetch was deliberately not run (the positive anchor is a ~6.3 GB archive — an explicit human decision). Both anchors carry the explicit `PENDING-FETCH` sentinel; no digest was fabricated. Run `bootstrap_anchor_hashes()` when online and authorized and replace the sentinel BEFORE Phase 16 opens the sealed holdout.

### Phase 9: Cross-Method Comparator Harness
**Goal**: A reproducible harness that runs the classical estimators and a Tapqir bridge on shared inputs, so v2.0 can be positioned as "knows when the classics are wrong," not merely "agrees with them"
**Depends on**: Nothing (parallelizable now)
**Requirements**: CMP-01, CMP-02, CMP-03, CMP-04, CMP-05, CMP-06, CMP-07, CMP-08, CMP-09
**Success Criteria** (what must be TRUE):
  1. Costes-p, Manders M1/M2, Pearson, Spearman all run on one shared input and emit a per-method comparison table
  2. A Tapqir bridge reproduces a published Tapqir example as a sanity anchor
  3. The harness reuses the BayesInteractomics comparator/audit pattern and is seeded/reproducible
**Plans**: 6 plans (4 waves: foundation -> {estimators, inputs, Tapqir} -> table/audit -> entry point + gate)
- [x] 09-01-PLAN.md - Wave-0 foundation: promote DataFrames/CSV, pre-declare D-14 consts (config.jl), scaffold test_comparator.jl + resolve-risk gate (CMP-04, CMP-09)
- [x] 09-02-PLAN.md - Classical estimator battery + seeded Costes block-scramble p-value (CMP-01, CMP-02)
- [x] 09-03-PLAN.md - Seeded regime-labelled shared-input builder, input-source-agnostic (CMP-03)
- [x] 09-04-PLAN.md - Isolated Tapqir sub-env + graceful skip-with-flag bridge + anchor capture (CMP-06)
- [x] 09-05-PLAN.md - Comparison table + divergence/traffic-light + content-addressed CSV/JLD2 + audit (CMP-01, CMP-04, CMP-05, CMP-07)
- [x] 09-06-PLAN.md - Seeded run_comparator entry point + filled D-13 test gate (CMP-08, CMP-09)

### Phase 10: Manuscript Skeleton and Related-Work Positioning
**Goal**: A compiling manuscript skeleton with related-work positioning drafted early — especially the explicit delta versus Tapqir/Costes/Manders — so experiments are shaped by the claims they must support
**Depends on**: Nothing (parallelizable now)
**Requirements**: TBD (claims derive from existing REQ outcomes; plans tag CONTEXT decisions D-01..D-12 + success criteria SC1/SC2/SC3)
**Success Criteria** (what must be TRUE):
  1. A Typst skeleton compiles with section scaffolding and a claim table
  2. Related work drafts the Tapqir differentiation (amortization + SBC + registration-UQ + spatial map)
  3. Figure specifications enumerate the panels each later phase must deliver
**Plans**: 5 plans (4 waves)
- [x] 10-01-PLAN.md — Compiling venue-neutral Typst skeleton: main.typ + section stubs + seeded refs.bib + claim_table machinery + compile gate (D-01..D-05, D-11, D-12)
- [x] 10-02-PLAN.md — Claim spine content: honestly-caveated claim table rows with per-row axis/phase/figure/status (D-06, D-07)
- [x] 10-03-PLAN.md — Related-work matrix + four-axis Tapqir differentiation, ASSUMED cells verified vs eLife 73860 (D-08, D-09)
- [x] 10-04-PLAN.md — Figure specifications enumerating per-phase panels linked to claim rows (D-10)
- [x] 10-05-PLAN.md — Integration regression: content-asserting compile gate + final DoD (D-12)

### Phase 11: Registration and Chromatic Uncertainty as Latent
**Goal**: Promote sub-pixel registration (and an optional chromatic warp) from a fixed simulator nuisance to an inferred latent, so the coloc posterior widens honestly under registration uncertainty instead of reporting false confidence
**Depends on**: Phase 7
**Requirements**: TBD
**Success Criteria** (what must be TRUE):
  1. dx/dy (+ optional 1-param chromatic warp) is added to θ (extending `spike/simulator/forward.jl` stage 6 + `prior.jl`) and the NPE is retrained on the extended prior
  2. Posterior width increases monotonically with injected registration uncertainty on a controlled sweep
  3. Deliberately mis-registered test images are handled without silent overconfidence
**Plans**: 11 plans (9 waves — W1 pre-registration+golden ; W2 the D-11/D-12 single commit ; W3 regression+provenance ∥ research-net scaffold ; W4 pre-flight probe+Tier-2 ; W5 probe-verdict checkpoint ; W6 λ-hierarchical datagen+training ; W7 SC2 ladder ; W8 SC3 breakdown+attenuation ∥ real-image ; W9 report+docs)
- [x] 11-01-PLAN.md — Tier-1 pre-registration consts, TOST/Wilson/Holm-direction stats, pre-edit golden fixture + PHASE11_BASE_SHA (D-01, D-04, D-07, D-08, D-10, D-14)
- [x] 11-02-PLAN.md — THE D-11/D-12 SINGLE COMMIT: chromatic ε as an 8th θ column, single composed affine stage 6 in both simulators, θ-arity ripple through src/ and the root suite, named limit #8 with the pinned pre-ε sha (D-02, D-09, D-10, D-11, D-12, D-15)
- [x] 11-03-PLAN.md — Exact-equality stage-6 regression vs the pre-edit golden, D-12 sha verification, the five §F20c decoupling commands, shipped-bundle load + colocalization_amortized regression (D-01, D-10, D-11, D-12, D-16)
- [x] 11-04-PLAN.md — Research-net scaffold: BoundedThetaTransform port + λ encoder/augmenter, d_in=129/D=8 smoke train, the λ-ablation tripwire (D-01, D-02, D-03, D-04)
- [x] 11-05-PLAN.md — D-06 pre-flight probe (paired, F5 mixture + 256² arm) and the append-only Tier-2 probe-derived constants (D-04, D-06, D-07, D-08, D-09)
- [x] 11-06-PLAN.md — Probe verdict: mechanical abort-criterion evaluation, blocking branch decision, iteration ledger (D-04, D-06)
- [ ] 11-07-PLAN.md — λ-hierarchical datagen (50k, F5 mixture) + research-net training + the λ-ablation tripwire firing green (D-01, D-02, D-03, D-09)
- [ ] 11-08-PLAN.md — SC2 ladder (GATED): λ-aware harness, permutation-Spearman monotonicity, per-rung TOST equivalence with the inverted Holm direction, per-rung shrinkage/vacuity (D-05, D-07, D-08, D-13)
- [ ] 11-09-PLAN.md — SC3 beyond-prior breakdown curve at λ=3.0 and λ=1.0 + the D-02 ρ_true attenuation measurement, both reported-not-gated (D-02, D-13, D-14, D-15)
- [ ] 11-10-PLAN.md — D-17 qualitative real-image transfer check on test/test_images/ with null-warp and integer-offset controls, sealed corpus holdout untouched (D-13, D-15, D-17)
- [ ] 11-11-PLAN.md — Figures, the Phase-11 results report, and the D-16 interpretation section in docs/amortized.md (D-01, D-05, D-12, D-13, D-14, D-16, D-17)

### Phase 12: Spatial Colocalization Map (GP/CAR)
**Goal**: Replace exchangeable patch pooling with a spatial lattice prior over the correlation grid, producing an amortized per-region Δρ map with calibrated per-region uncertainty — the feature that makes v2.0 "spatial" and differentiates it from Tapqir
**Depends on**: Phases 7, 11 (file overlap: both retrain/modify the shared PosteriorEstimator + simulator/training path — serialized to avoid a merge collision; spatial map trains on the registration-aware θ)
**Requirements**: TBD
**Success Criteria** (what must be TRUE):
  1. A lattice prior (CAR vs. AbstractGPs, chosen by mini-spike) is placed over the `correlation()` grid with a CNN/DeepSet summary
  2. `coloc_map(...)` returns a Δρ map + uncertainty map, amortized in a forward pass
  3. The spatial (CAR) model beats independent pooling in coverage on ≥1 real image
  *(Highest effort-risk phase; the natural descope-to-v2.1 candidate if amortization stalls.)*
**Plans**: TBD
- [ ] TBD (run /gsd:plan-phase 12 to break down)

### Phase 13: Three-Hypothesis Amortized Bayes Factor
**Goal**: Extend the amortized evidence network from two- to three-way model comparison — colocalized / random / mutually-exclusive — so segregation becomes a first-class testable hypothesis, replacing the fragile KDE+quadgk Bayes factor
**Depends on**: Phases 7, 11 (D-02: the three-way evidence net trains on Phase 11's registration-aware frozen `zt` and inherits its uncertainty conditioning input — execution blocks on Phase 11's research net existing; the original entry named Phase 7 only)
**Requirements**: TBD
**Success Criteria** (what must be TRUE):
**AMENDED** — SC1, SC2 and SC3 are superseded/scoped by `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md` (D-12/D-15, frozen before any Phase-13 result). The original text below is retained and must be cited alongside any amended result.
  1. A 3-way `RatioEstimator`/evidence network emits a log-BF simplex over {coloc, random, exclusion} in one forward pass (see amendment section 4)
  2. It reproduces `compute_BayesFactor()` (`src/bayes.jl:109`) in the overlapping 2-way regime without quadgk/KDE (see amendment sections 1-3)
  3. The exclusion hypothesis is validated on segregated ground-truth inputs (see amendment section 5 — three arms; the simulator ground truth gates, the alpha series and the `test/test_images/` check do not)
**Plans**: 16 plans (9 waves — W1 pre-registration ∥ SC2 amendment ; W2 labels ∥ α-series ∥ τ-probe code ; W3 two-head net ∥ real-image ingestion ; W4 D-07 closed-form verification ; W5 result type + suite wiring ; W6 Phase-11 preconditions ∥ τ-probe run [BLOCKED] ; W7 datagen + training [BLOCKED] ; W8 reported gate ∥ α-series run ∥ real-image run [BLOCKED] ; W9 phase report [BLOCKED])
- [x] 13-01-PLAN.md — Tier-1 pre-registration consts (seeds, τ probe spec, D-07 design + bars, gate floors, ECE band, α ladder) + literal-assertion test (D-01, D-04, D-05, D-06, D-07, D-09, D-10, D-12, D-13, D-14, D-15, D-16)
- [x] 13-02-PLAN.md — 13-SC2-AMENDMENT.md frozen before any result (two outcome-independent defects) + ROADMAP annotation (D-12, D-02, D-08, D-09, D-13, D-15)
- [ ] 13-03-PLAN.md — Three-way label surface: two-factor cut on ρ_sample level × control contrast, head_targets, measure_head_log_odds (D-05, D-07, D-08, D-11)
- [ ] 13-04-PLAN.md — α-graded disjoint-reassignment transform + its four invariants (bitwise α=0, no NEW zero (count(iszero) preserved), intensity conservation, frozen mask, substrate-agnostic (one code path for simulated and real)) (D-15, D-16)
- [ ] 13-05-PLAN.md — Fit-free UNPAIRED τ resolution probe code + abort criterion, needs no trained net (D-06, D-04)
- [ ] 13-06-PLAN.md — ThreeWayEvidenceNet: shared trunk verbatim + two BCE heads + masked joint loss + per-head read surface, with the A5 train smoke (D-08, D-09, D-10, D-11, D-03)
- [ ] 13-07-PLAN.md — D-07 correction VERIFIED against a closed-form Gaussian toy, with both negative controls (uncorrected logit, within-class reshape) (D-07, D-11)
- [ ] 13-08-PLAN.md — ThreeHypothesisColocResult in spike/ (log_bf_vs_random), empty-bin MCE trap demonstrated, resolve-risk clause (i) + suite wiring (D-08, D-13, D-14, D-04)
- [ ] 13-09-PLAN.md — Phase-11 precondition binding: loud block, no grid-8 fallback, derived conditioning length (D-02, D-03) [BLOCKED: Phase 11]
- [ ] 13-10-PLAN.md — REPORTED τ probe run + Tier-2 append of the measured τ (two-commit discipline) (D-06, D-04) [BLOCKED: Phase 11 simulator]
- [ ] 13-11-PLAN.md — Class-frequency-stratified conditioned datagen on the frozen Phase-11 zt + the single training run with measured per-head log-odds (D-02, D-03, D-07, D-10, D-11) [BLOCKED: Phase 11]
- [ ] 13-12-PLAN.md — REPORTED amended gate: per-class AUC + confusion matrix (descriptive) + per-head ECE with vacuous-pass guard; λ response and binary-NRE continuity reported not gated (D-12, D-13, D-14, D-03, D-05) [BLOCKED: Phase 11]
- [ ] 13-13-PLAN.md — REPORTED α-ladder run through the trained net: log-BF curves, m̄(α), crossing point α*, sealed holdout untouched (D-15, D-16) [BLOCKED: Phase 11]
- [ ] 13-14-PLAN.md — Phase-13 report: amended criteria beside originals, the D-05 trap in prose, four named limits incl. no labelled real-data validation, the real-image OOD + naming-correction honesty items, iteration ledger, deferrals, all 8 open questions closed (D-05, D-07, D-08, D-12, D-13, D-15, D-16, D-01, D-02, D-04) [BLOCKED: Phase 11]
- [ ] 13-15-PLAN.md — D-15 (AMENDED) real-image ingestion from `test/test_images/`: P13_REPO_ROOT/real_tif/load_real, ghat-anchor regression, derived grid truncation, the real α-ladder, and the five D-15 testsets (ingestion, real ladder, anti-snooping, read-only, not-a-gate) + A14 falsifier — Phase-11-INDEPENDENT (D-15, D-16, D-01, D-04)
- [ ] 13-16-PLAN.md — REPORTED real-image qualitative check: λ sweep with every log-BF printed beside its OOD verdict, Phase-13-vs-shipped OOD comparison, real α-ladder through the net, naming correction + target-substitution record, seal intact (D-15, D-16, D-03, D-08, D-01, D-04) [BLOCKED: Phase 11]

### Phase 14: Decision and Abstention Layer
**Goal**: Turn calibrated posteriors + the 3-way BF into an actionable batch decision {coloc / not / ABSTAIN} at a controlled Bayesian FDR, abstaining exactly when the tool should be silent
**Depends on**: Phases 11, 12, 13
**Requirements**: TBD
**Success Criteria** (what must be TRUE):
  1. `decide_coloc(...)` emits calibrated calls at a user-set Bayesian FDR across a batch, using conformal sets (ConformalPrediction.jl) + decision-risk
  2. Abstention triggers on OOD ∨ cross-method disagreement ∨ ambiguous conformal set
  3. A monotone risk-coverage curve shows abstention concentrates on hard/OOD cases
**Plans**: TBD
- [ ] TBD (run /gsd:plan-phase 14 to break down)

### Phase 15: Calibration Operating Envelope and CI Gate
**Goal**: Map the tool's domain of applicability by adversarially sweeping nuisances until coverage breaks, prove the OOD flag fires before it does, and lock calibration into CI as a regression gate
**Depends on**: Phases 11, 12
**Requirements**: TBD
**Success Criteria** (what must be TRUE):
  1. An adversarial sweep over spillover/PSF/autofluorescence/registration yields a domain-of-applicability map
  2. The OOD flag demonstrably fires before empirical coverage breaks ("OOD-before-break")
  3. An SBC/coverage regression test fails CI on calibration drift after a code change
**Plans**: TBD
- [ ] TBD (run /gsd:plan-phase 15 to break down)

### Phase 16: External Validation and Manuscript Assembly
**Goal**: Close v2.0 with a blind external evaluation against the physical corpus and comparator harness, and assemble the reproducible manuscript package
**Depends on**: Phases 8, 9, 10, 14, 15
**Requirements**: TBD
**Success Criteria** (what must be TRUE):
  1. v2.0 is blind-evaluated against the physical anchors (Phase 8) and comparator harness (Phase 9); results reported honestly, including the simulator-validated mid-range caveat
  2. All manuscript figures are assembled into the Phase-10 skeleton
  3. `run_v2.jl` reproduces the end-to-end result from a fixed seed; a Zenodo/DOI release is prepared
**Plans**: TBD
- [ ] TBD (run /gsd:plan-phase 16 to break down)
