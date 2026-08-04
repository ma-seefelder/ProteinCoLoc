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

- [x] **Phase 11: Registration + Chromatic Uncertainty as Latent** - Promote dx/dy (+ chromatic warp) to inferred θ; posterior widens honestly under registration uncertainty — **CLOSED NEGATIVE-BUT-USEFUL 2026-07-27**: registration at ≤3 px is not inferable from the 8×8 patch summary at any coloc level, AND does not need to be (Δρ RMSE flat in λ, ratio 1.0003). SC1g gate MIS-SPECIFIED, not failed. Plans 11-08…11-11 SUPERSEDED. See `11-CLOSURE.md`
- [x] **Phase 12: Spatial Colocalization Map (GP/CAR)** - Lattice prior over the correlation grid → amortized per-region Δρ map + uncertainty (descope-to-v2.1 candidate) — **CLOSED NEGATIVE 2026-08-03: THE GOAL WAS NOT ACHIEVED.** No spatial prior (CAR or GP) ever beat the neutralized ablation: `NONE-BEATS-ABLATION` fired in three independent mini-spike runs, so the shipped deliverable IS the ablation — the *same* network with spatial borrowing switched off — wearing the `SpatialColocResult` type, and **no trained spatial arm exists at any scale**. SPAT-06 was measured and FAILED its pre-registered band: pooled leave-region-out coverage **0.9707** against **[0.87, 0.93]**. SPAT-07 and SPAT-08 are DEFERRED TO v2.1 by recorded user rulings dated 2026-07-31 — not skipped, not forgotten. See `12-VERIFICATION.md`
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
| 12. Spatial Colocalization Map (GP/CAR) | 11/20 | In Progress|  |
| 13. Three-Hypothesis Amortized Bayes Factor | 16/17 | In Progress| 13-16: real-image arm REPORTED — **all its figures SUPERSEDED**, measured on the DAPI-counterstain pair (RANDOM on both pairs, both OOD-flagged 2.49x/2.42x, agreement). 13-17 corrected the operative pair to c2/c3 green/red per `13-D15-AMENDMENT.md`, **dropped** the redundancy arm (no second pair exists) and **re-ran the arm once**: on the amended pair the verdict is **COLOC / RANDOM** (was RANDOM / RANDOM) and both fixtures are OOD-flagged **4.198x / 5.297x — WORSE than the superseded 2.49x / 2.42x**. **The CONCLUSION is unchanged**: qualitative, n = 2, unlabelled, OOD-bound. **The gating arm is unaffected and no threshold moved**; D-04 allowance unspent |
| 14. Decision + Abstention Layer | 9/14 | In Progress|  |
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
**OUTCOME 2026-07-27 — the goal is ANSWERED NEGATIVELY, with evidence, and the answer is more useful than the one the phase expected.** SC1 was *delivered* (θ extended to 8 columns at `ca02b0e`, research net trained), but SC2/SC3 rest on a premise the phase disproved: registration at ≤3 px carries no recoverable information in the 8×8 patch-correlation summary. Read each criterion below against that.

  1. dx/dy (+ optional 1-param chromatic warp) is added to θ (extending `spike/simulator/forward.jl` stage 6 + `prior.jl`) and the NPE is retrained on the extended prior — **DONE** (research net only; the shipped bundle was never retrained)
  2. Posterior width increases monotonically with injected registration uncertainty on a controlled sweep — **MOOT.** The width the data demands is *flat* in λ (Δρ RMSE 0.13334 → 0.13338, ratio 1.0003). A monotone increase would have been *dishonest* width, not honest width
  3. Deliberately mis-registered test images are handled without silent overconfidence — **SATISFIED, by a different route than planned.** There is no false confidence to correct: the posterior already reports approximately the right width (1.021) for a quantity whose correct width ratio is ~1.00. Registration must instead be calibrated externally (fiducial beads), which is standard microscopy practice

**Plans**: 11 plans (9 waves — W1 pre-registration+golden ; W2 the D-11/D-12 single commit ; W3 regression+provenance ∥ research-net scaffold ; W4 pre-flight probe+Tier-2 ; W5 probe-verdict checkpoint ; W6 λ-hierarchical datagen+training ; W7 SC2 ladder ; W8 SC3 breakdown+attenuation ∥ real-image ; W9 report+docs)

- [x] 11-01-PLAN.md — Tier-1 pre-registration consts, TOST/Wilson/Holm-direction stats, pre-edit golden fixture + PHASE11_BASE_SHA (D-01, D-04, D-07, D-08, D-10, D-14)
- [x] 11-02-PLAN.md — THE D-11/D-12 SINGLE COMMIT: chromatic ε as an 8th θ column, single composed affine stage 6 in both simulators, θ-arity ripple through src/ and the root suite, named limit #8 with the pinned pre-ε sha (D-02, D-09, D-10, D-11, D-12, D-15)
- [x] 11-03-PLAN.md — Exact-equality stage-6 regression vs the pre-edit golden, D-12 sha verification, the five §F20c decoupling commands, shipped-bundle load + colocalization_amortized regression (D-01, D-10, D-11, D-12, D-16)
- [x] 11-04-PLAN.md — Research-net scaffold: BoundedThetaTransform port + λ encoder/augmenter, d_in=129/D=8 smoke train, the λ-ablation tripwire (D-01, D-02, D-03, D-04)
- [x] 11-05-PLAN.md — D-06 pre-flight probe (paired, F5 mixture + 256² arm) and the append-only Tier-2 probe-derived constants (D-04, D-06, D-07, D-08, D-09)
- [x] 11-06-PLAN.md — Probe verdict: mechanical abort-criterion evaluation, blocking branch decision, iteration ledger (D-04, D-06)
- [~] 11-07-PLAN.md — λ-hierarchical datagen (50k, F5 mixture) + research-net training + the λ-ablation tripwire firing green (D-01, D-02, D-03, D-09) — **RAN, then BLOCKED on the SC1g tripwire.** `11-07-SUMMARY.md` stays `status: blocked` as the honest record of where execution stopped; the tripwire that stopped it was subsequently shown to be mis-specified (`11-BAR-DERIVATION.md`)
- [S] 11-08-PLAN.md — SC2 ladder (GATED) — **SUPERSEDED, NEVER RUN.** Entry gate mis-specified; SC2 already answered analytically (≤1.134) and empirically (RMSE ratio 1.0003)
- [S] 11-09-PLAN.md — SC3 breakdown + D-02 attenuation — **SUPERSEDED IN SUBSTANCE, NEVER RUN.** The attenuation question *is* the flat-RMSE finding
- [S] 11-10-PLAN.md — D-17 real-image transfer check — **NOT APPLICABLE, NEVER RUN.** Nothing is inferred, so nothing transfers
- [S] 11-11-PLAN.md — Figures, report, docs — **SUPERSEDED BY `11-DIAGNOSIS.md`, NEVER RUN.** No `11-REPORT.md` exists or is owed; `docs/amortized.md` deliberately unedited (D-01 holds)

**Closure**: `11-CLOSURE.md` (2026-07-27) — evidence in `11-DIAGNOSIS.md`, corrected gate spec in
`11-BAR-DERIVATION.md`, structural findings in `deferred-items.md`. Legend: `[S]` = superseded,
`[~]` = ran but blocked. **No retraining, no re-seed, no reship**; `spike/validation/p11_consts.jl`
and the shipped `amended_v2/grid_8` bundle are byte-unchanged and the Phase-7 GO is untouched.
**Phase 15 must start from `11-BAR-DERIVATION.md`, not from `p11_consts.jl`** — the constant there is
correct as history and wrong as a specification, and a re-derived bar must be an equivalence-style
test, not a "must exceed" threshold, because the true value of the gated quantity is ~1.00.

### Phase 12: Spatial Colocalization Map (GP/CAR)

**Goal**: Replace exchangeable patch pooling with a spatial lattice prior over the correlation grid, producing an amortized per-region Δρ map with calibrated per-region uncertainty — the feature that makes v2.0 "spatial" and differentiates it from Tapqir
**Depends on**: Phase 7 — the substantive dependency. ~~Phases 7, 11 (file overlap: both retrain/modify the shared PosteriorEstimator + simulator/training path — serialized to avoid a merge collision; spatial map trains on the registration-aware θ)~~ **VOID** per the 2026-07-27 premise audit (V-1, V-2, V-3): no registration-aware trained model exists (plans 11-08…11-11 were superseded by the diagnosis) and Δρ RMSE was measured flat in λ (ratio 1.0003), so there is nothing to train on and nothing gained; Phase 11 is closed, so there is no concurrent writer to serialize against. What DOES stand (K-1): Phase 11's composed affine stage 6 landed at `ca02b0e` and is this phase's baseline, which D-06's field application stacks on.
**Requirements**: SPAT-01, SPAT-02, SPAT-03, SPAT-04, SPAT-05, SPAT-06, SPAT-07, SPAT-08, SPAT-09
**Success Criteria** (what must be TRUE):
**AMENDED** — SC1, SC2 and SC3 are superseded/scoped by `.planning/phases/12-spatial-colocalization-map/12-SC3-AMENDMENT.md` (D-02/D-09/D-10/D-11, frozen before any Phase-12 result). The original text below is retained and must be cited alongside any amended result.
**OUTCOME 2026-08-03 — the goal is ANSWERED NEGATIVELY. `12-VERIFICATION.md`'s headline is "NO — the phase did not achieve its stated goal."** The per-region Δρ map mechanism was built and works, but the two other things the goal sentence promises are false on the record. Read each criterion below against that. This is not read as an execution failure: every outcome came from a pre-registered, executed measurement.

  1. A lattice prior (CAR vs. AbstractGPs, chosen by mini-spike) is placed over the `correlation()` grid with a CNN/DeepSet summary — **MECHANISM BUILT; NO ARM WAS EVER SELECTED.** Both prior arms are implemented as dense linear algebra rather than `AbstractGPs.jl` (per the pre-registered `12-SC3-AMENDMENT.md` §2) — `car_sigma` at `spike/simulator/p12_lattice.jl:125`, `gp_sigma` at `:153` — and the lean CNN summary net exists (`build_p12_summary_net` at `spike/npe/p12_architecture.jl:253`). But the "chosen by mini-spike" clause resolved to **neither arm**: `select_prior` returned `NONE-BEATS-ABLATION` in all three runs (unseeded 18-epoch, seeded 18-epoch control, seeded 100-epoch treatment)
  2. `coloc_map(...)` returns a Δρ map + uncertainty map, amortized in a forward pass — **DELIVERED.** `p12_coloc_map` (`spike/validation/p12_coverage.jl:990`) returns the per-region Δρ map together with per-region uncertainty, each of the sample and control stacks one amortized forward pass, verified structurally. **Named caveat that travels with it (`12-VERIFICATION.md` §SPAT-05):** the model behind it is the non-spatial D-13 ablation, not a trained spatial arm, and it exists at spike scale (10,000 pairs), not the 50,000-pair version
  3. The spatial (CAR) model beats independent pooling in coverage on ≥1 real image — **NOT MET, AND UNREACHABLE AS WRITTEN.** Two independent reasons: (i) there is **no spatial model to compare against** — see criterion 1; (ii) the real-image arm (**12-19**) **never ran**, hard-blocked on `p12_train_full_report.jld2`, confirmed absent on disk, which only the never-executed **12-17** produces. The reduced substitute that DID run (12-16, on simulated data, ablation against a prior-only floor) failed its own pre-registered band — `12-16-SUMMARY.md`'s own words: **"SPAT-06 IS NOT MET: pooled leave-region-out coverage 0.9707 lies outside [0.87, 0.93]."** Not a boundary call: the interval on the independent unit is [0.9675, 0.9738], clear of the band by ≈ 0.037
  *(Highest effort-risk phase; the natural descope-to-v2.1 candidate if amortization stalls.)*
**Plans**: 20 plans exist and are committed. **16 executed (12-01 … 12-16); 12-17 … 12-20 are FORECLOSED and NEVER RAN.**

- [x] 12-01-PLAN.md — Freeze the Tier-1 Phase-12 pre-registration before a single Phase-12 number exists, and wire the Phase-12 test surface into `spike/test/runtests.jl` at the one position where it can actually run (SPAT-09)
- [x] 12-02-PLAN.md — Freeze the amended success criteria before any Phase-12 number can influence their wording, and mint the nine SPAT requirement IDs the ROADMAP had left unassigned — produced `12-SC3-AMENDMENT.md` (SPAT-09)
- [x] 12-03-PLAN.md — The lattice machinery: CAR and GP covariance arms, the D-05 per-cell rescale, the lag-1 reparametrization of the correlation length, the orthonormal DCT-II basis — `car_sigma` at `spike/simulator/p12_lattice.jl:125`, `gp_sigma` at `:153` (SPAT-01, SPAT-04)
- [x] 12-04-PLAN.md — Give the network a spatial *view* of the unchanged 128-row summary, topology pinned literally so nobody builds the measured-12.5 h variant — `build_p12_summary_net` at `spike/npe/p12_architecture.jl:253` (SPAT-03, SPAT-04)
- [x] 12-05-PLAN.md — Make the phase's three hard constraints executable rather than aspirational: `src/` untouched, the spike environment byte-frozen, the Phase-16 sealed holdout unburned (SPAT-09)
- [x] 12-06-PLAN.md — Answer cheaply, before anything expensive happens, whether `chromatic_eps` is identifiable from the patch-correlation summary. **Gates nothing.** MEASURED: `eps_ratio = 0.96307` — barely identifiable, 3.7 % better than the empirical prior-mean baseline — and `ridge_residual_shrinkage = 0.96870 > 0.90` ⇒ `vacuous = true` (SPAT-08 precursor, in scope, not deferred)
- [x] 12-07-PLAN.md — Turn 12-03's lattice arithmetic into a prior over ρ **fields** that leaves every region's marginal exactly where the Turing μ-prior puts it (D-05 copula), with a cheap route to pixel resolution (SPAT-01)
- [x] 12-08-PLAN.md — Capture the pre-edit simulator golden, then make stage 1 accept a spatially-varying ρ without moving a single byte on the path that does not use one — golden matched **12/12** exactly and a constant-valued field matched with `maximum(abs, Δ) = 0.0` against a `1e-12` tolerance (SPAT-02)
- [x] 12-09-PLAN.md — The field-aware training pool: draw a correlation length, draw a field conditional on it, simulate, encode the unchanged 128-row summary, and store the **drawn lattice** beside it; the single place the 72-row θ vector is assembled (SPAT-01, SPAT-02)
- [x] 12-10-PLAN.md — Prove at the pre-registered sample size, against a tolerance frozen before the field simulator existed, that D-05's copula leaves **every one of the 64 regions** with the `MU_PRIOR` marginal — max per-region W1 = **0.00675** at M = 20,000 across all 11 arm×rung configurations, against `P12_SIM02_W1_TOL_PERREGION = 0.10` (≈15× headroom) (SPAT-01)
- [x] 12-11-PLAN.md — The D-12 Stage-1 question the whole phase rests on — *do the other 63 regions tell you anything about the one you held out?* — answered by a closed-form leave-region-out ridge probe before a single network is trained. **COMPLETE IN FACT**: all artifacts on disk (`run_p12_stage1_ridge.jl`, `p12_stage1_report.jld2`, the blocked-adjudication record) and `12-STAGE1-VERDICT.md:3` reads `VERDICT: PROCEED` at commit `ee78cfe`. **ITS SUMMARY FRONTMATTER STILL READS `status: blocked`. THAT IS STALE, AND IT IS DELIBERATELY LEFT BYTE-UNCHANGED** per this phase's append-never-overwrite discipline — corrections are appended, records are never rewritten. **This is precisely why `init.execute-phase` counts 5 incomplete plans rather than 4.** **Do NOT "fix" `12-11-SUMMARY.md`**; the staleness is resolved here, in the roadmap, on purpose (SPAT-02)
- [x] 12-12-PLAN.md — Give the phase a real result type: a per-region Δρ map **with** a per-region uncertainty map, as a new subtype of the package's actual result hierarchy, without touching `src/` — `SpatialColocResultSpike` in `spike/p12/result.jl` (SPAT-05)
- [x] 12-13-PLAN.md — Turn D-08's hedge (*"expect this column may prove unidentifiable; report it as such if so"*) into a pre-registered measurement of the correlation length, using a probe strong enough that a null means something. **Gates nothing** (SPAT-04)
- [x] 12-14-PLAN.md — Build the one training surface every later Phase-12 plan uses, making a leaked standardizer and an out-of-distribution held-out encoding impossible rather than merely discouraged (arm as a keyword, everything else identical by construction). **Its Task 3 — the only coded producer of the D-13 descope bundle — is keyed on `p12_stage1_verdict() === :descope` and therefore RAN AS A NO-OP under the `:proceed` verdict.** That is the wiring defect that forced the standalone `12-D13-AUTHORISATION.md` (SPAT-03, SPAT-04, SPAT-09)
- [x] 12-15-PLAN.md — Decide by measurement rather than preference which lattice prior the phase uses, running both arms through a full-rank `K = P12_K_DEV = 63` head so the comparison measures the prior and not the head. **OUTCOME: `select_prior` returned `NONE-BEATS-ABLATION` in all three runs** — the original unseeded 18-epoch run, a seeded 18-epoch control, and a seeded 100-epoch treatment (`12-MINISPIKE-VERDICT.md`, `12-15-RERUN-VERDICT.md`). The qualifier *"at mini-spike scale"* stays on it, and a comparison between two failed arms is not evidence about priors: **the CAR-vs-GP question is UNRESOLVED, not answered** (SPAT-01, SPAT-04)
- [x] 12-16-PLAN.md — Build the D-09 leave-region-out predictive machinery, calibrate its one modelling constant, and measure the comparison on simulation, where ground truth still exists to check the machinery measures what it claims. **OUTCOME: pooled coverage 0.9707 outside the pre-registered [0.87, 0.93] — "SPAT-06 IS NOT MET"**, reported and not tuned; the ablation beats `prior_only_floor` by **+0.19195 nats/region** against `P12_STAGE2_LOGSCORE_MIN = 0.02`; `stage2_gate = :not_applicable_descope` (there is no spatial arm to compare) and `spat07_scope = spat08_scope = :deferred_to_v2_1` recorded in the artifact, not silently omitted (SPAT-05, SPAT-06)
- [S] 12-17-PLAN.md — Produce the two trained networks the rest of the phase scores: the spatial model on the prior 12-15 selected, and its D-10 matched ablation. **FORECLOSED, NEVER RAN.** It keys its 50,000-pair pool on `P12_CHOSEN_PRIOR` (`12-17-PLAN.md:118`, `p12_full_arms()` = `(P12_CHOSEN_PRIOR, :none)`), which is *asserted absent* at `spike/validation/p12_consts.jl:798`. Dispatching it throws. Superseded by the 2026-07-31 descope (SPAT-03, SPAT-04)
- [S] 12-18-PLAN.md — Report per-region calibration in the one space where this phase's own prior construction does not manufacture artifacts, with a test matched to what each row actually is. **DEFERRED TO v2.1** by user ruling dated 2026-07-31, `12-D13-AUTHORISATION.md` §9. **NEVER RAN**: `spike/test/test_p12_sbc.jl` remains a self-declaring `P12_PENDING_SCAFFOLD`, and no `p12_sbc.jl` runner exists anywhere in `spike/` (SPAT-07)
- [S] 12-19-PLAN.md — Score the leave-region-out predictive comparison on genuine microscopy data — the only part of SC3 a simulation cannot supply — and state the claim at exactly the strength the data support. **FORECLOSED, NEVER RAN.** `12-19-PLAN.md:202` resolves and hash-verifies **both bundles** from `p12_train_full_report.jld2`, which only 12-17 produces and which is confirmed absent on disk. Dead behind 12-17 (SPAT-06)
- [S] 12-20-PLAN.md — Test the explanation that would make this phase's headline result meaningless (the S-4 chromatic-radial confound), then record the D-12 Stage-2 verdict. **DEFERRED TO v2.1** by user ruling dated 2026-07-31, `12-D13-AUTHORISATION.md` §7. **NEVER RAN**, and the consequence is worth stating: **the phase therefore makes NO radial-confound claim at all**, and the S-4 warning carried forward from Stage 1 is **left standing and unaddressed** — nothing in the D-13 release should be read as having ruled a radial confound in or out (SPAT-08, SPAT-09)

**Closure**: `12-VERIFICATION.md` (2026-08-03, `status: gaps_found`; score **8/9 must-haves** — 6 verified,
2 accepted as deferred-to-v2.1 by recorded user ruling, **1 FAILED**) and `12-D13-AUTHORISATION.md` (the
standalone, explicitly-dated user authorisation that created the D-13 ablation as a NEW record without
touching any gate file). Legend: `[x]` = complete, `[S]` = superseded/foreclosed, `[~]` = ran but blocked
(not used in this phase).

**THE NAMED LIMIT, in the exact form `12-D13-AUTHORISATION.md` §10.5 fixes it for the manuscript:**
> **"We cannot exclude that longer training would have brought the spatial arms into the coverage band."**

**That sentence is a NAMED LIMIT of v2.0, not a caveat.** It is weaker than what §8 claimed, and it is
recorded in final form so nobody has to re-derive it later from a chain of appended corrections.

**SPAT-06 IS NOT TO BE RE-TUNED** by adjusting epochs, N, thresholds or the noise model (§10.1, §10.5).
Three specific moves were each offered and each REFUSED: **widening the admissibility loop** at
`spike/validation/run_p12_minispike.jl:249` (`for a in (:car, :gp)`) to include `:none` after seeing
`:none` fail it; **re-deriving the [0.87, 0.93] band** after seeing the arms miss it; and **relaxing the
early-stopping criterion** in order to train longer. **The loop stays wrong and documented as wrong.**
The user accepted a weaker negative rather than change a modelling constant after seeing the numbers.

**§10.1's honest statement — stronger than "the spatial priors were miscalibrated", and the sentence
that should survive:** *NO ARM IS CALIBRATED AT 100 EPOCHS, AND THE FALLBACK WAS RETAINED WITHOUT BEING
TESTED.* At 100 epochs all three arms lie outside [0.87, 0.93], `car` is closer to nominal than the
selected `none` (0.0492 against 0.0509), and both spatial arms beat the selected arm on RMSE — the arm
that won **won by default** and was never held to the criterion that disqualified the others.

**The `[x]` above means CLOSED, NOT SUCCEEDED — and the distinction is the whole point of this
block.** Phase 12 is finished in the only sense that matters operationally: nothing further will run,
all four remaining plans are foreclosed, and no measurement will be repeated. It is *not* finished in
the sense the goal sentence promises. **Five of nine requirements were never delivered** — SPAT-03,
SPAT-04 and SPAT-06 undelivered at the scope the goal means; SPAT-07 and SPAT-08 deferred to v2.1 —
and they will never be marked delivered. The 12-11 frontmatter staleness will never be repaired.

**A tally that counts this phase toward "13/16 complete" is counting a closed negative as a
success.** Any milestone summary, progress report or manuscript claim that leans on the phase count
must carry the qualifier, because the count alone cannot. **Phase 12 delivered 4 of 9 requirements
and answered its own goal question NO.**

*(Recorded 2026-08-03, superseding this block's original prediction that the phase would read
`partial` to the SDK permanently. That prediction was wrong: the SDK derives phase status from the
roadmap checkbox, not from plan/summary counts, so ticking `[x]` moved it straight to `complete` and
unblocked Phase 15. The checkbox is correct — the phase IS closed — and the prose is corrected here
rather than the checkbox reverted, since an unticked box misread the phase as unstarted, which is the
drift this amendment existed to end.)*

**Re-running `/gsd:execute-phase 12` is NOT indicated.** All four remaining plans are foreclosed by
rulings already on the record, not pending. `12-VERIFICATION.md` explicitly does not recommend executing
12-17, 12-18, 12-19 or 12-20, and does not recommend re-running, retuning or extending any existing
measurement.

### Phase 13: Three-Hypothesis Amortized Bayes Factor

**Goal**: Extend the amortized evidence network from two- to three-way model comparison — colocalized / random / mutually-exclusive — so segregation becomes a first-class testable hypothesis, replacing the fragile KDE+quadgk Bayes factor
**Depends on**: Phases 7, 11 (D-02: the three-way evidence net trains on Phase 11's registration-aware frozen `zt` and inherits its uncertainty conditioning input — execution blocks on Phase 11's research net existing; the original entry named Phase 7 only)
**Requirements**: TBD
**Success Criteria** (what must be TRUE):
**AMENDED** — SC1, SC2 and SC3 are superseded/scoped by `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md` (D-12/D-15, frozen before any Phase-13 result). The original text below is retained and must be cited alongside any amended result.

  1. A 3-way `RatioEstimator`/evidence network emits a log-BF simplex over {coloc, random, exclusion} in one forward pass (see amendment section 4)
  2. It reproduces `compute_BayesFactor()` (`src/bayes.jl:109`) in the overlapping 2-way regime without quadgk/KDE (see amendment sections 1-3)
  3. The exclusion hypothesis is validated on segregated ground-truth inputs (see amendment section 5 — three arms; the simulator ground truth gates, the alpha series and the `test/test_images/` check do not)

**Plans**: 17 plans (10 waves — W1 pre-registration ∥ SC2 amendment ; W2 labels ∥ α-series ∥ τ-probe code ; W3 two-head net ∥ real-image ingestion ; W4 D-07 closed-form verification ; W5 result type + suite wiring ; W6 Phase-11 preconditions ∥ τ-probe run [BLOCKED] ; W7 datagen + training [BLOCKED] ; W8 reported gate ∥ α-series run ∥ real-image run [BLOCKED] ; W9 phase report [BLOCKED] ; W10 D-15 channel-pair amendment + include-guard fix + real-arm re-run)

- [x] 13-01-PLAN.md — Tier-1 pre-registration consts (seeds, τ probe spec, D-07 design + bars, gate floors, ECE band, α ladder) + literal-assertion test (D-01, D-04, D-05, D-06, D-07, D-09, D-10, D-12, D-13, D-14, D-15, D-16)
- [x] 13-02-PLAN.md — 13-SC2-AMENDMENT.md frozen before any result (two outcome-independent defects) + ROADMAP annotation (D-12, D-02, D-08, D-09, D-13, D-15)
- [x] 13-03-PLAN.md — Three-way label surface: two-factor cut on ρ_sample level × control contrast, head_targets, measure_head_log_odds (D-05, D-07, D-08, D-11)
- [x] 13-04-PLAN.md — α-graded disjoint-reassignment transform + its four invariants (bitwise α=0, no NEW zero (count(iszero) preserved), intensity conservation, frozen mask, substrate-agnostic (one code path for simulated and real)) (D-15, D-16)
- [x] 13-05-PLAN.md — Fit-free UNPAIRED τ resolution probe code + abort criterion, needs no trained net (D-06, D-04)
- [x] 13-06-PLAN.md — ThreeWayEvidenceNet: shared trunk verbatim + two BCE heads + masked joint loss + per-head read surface, with the A5 train smoke (D-08, D-09, D-10, D-11, D-03)
- [x] 13-07-PLAN.md — D-07 correction VERIFIED against a closed-form Gaussian toy, with both negative controls (uncorrected logit, within-class reshape) (D-07, D-11)
- [x] 13-08-PLAN.md — ThreeHypothesisColocResult in spike/ (log_bf_vs_random), empty-bin MCE trap demonstrated, resolve-risk clause (i) + suite wiring (D-08, D-13, D-14, D-04)
- [x] 13-09-PLAN.md — Phase-11 precondition binding: loud block, no grid-8 fallback, derived conditioning length (D-02, D-03) — bound to `spike/npe/p11_research_npe.jld2`, n_cond = 1 derived, 74/74
- [x] 13-10-PLAN.md — REPORTED τ probe run + Tier-2 append of the measured τ (two-commit discipline) (D-06, D-04) — MEASURED τ = 0.15 at A(τ) = 0.916356 vs bar 0.9, reference λ = 3.0
- [x] 13-11-PLAN.md — Class-frequency-stratified conditioned datagen on the frozen Phase-11 zt + the single training run with measured per-head log-odds (D-02, D-03, D-07, D-10, D-11) — realized classes EXACTLY equal thirds; zt inherited; per-head corrections +0.0029420 / +0.0023543; run overfits early (best val at epoch 4), recorded not repaired
- [x] 13-12-PLAN.md — REPORTED amended gate: per-class AUC + confusion matrix (descriptive) + per-head ECE with vacuous-pass guard; λ response and binary-NRE continuity reported not gated (D-12, D-13, D-14, D-03, D-05) — **PASS, 6/6**: AUC 0.990169 / 0.988262 vs floor 0.9; ECE 0.0119808 / 0.012886 green vs 0.05, 0 empty bins, neither head vacuous; D-04 trigger did NOT fire (deep-tail |ρ|>0.9 exclusion AUC 0.997257), allowance UNSPENT; λ response flat but monotone (span 0.017 / 0.068 nats, Spearman 1.0); continuity corr 0.622 reported not gated
- [x] 13-13-PLAN.md — REPORTED α-ladder run through the trained net: log-BF curves, m̄(α), crossing point α*, sealed holdout untouched (D-15, D-16) — **α\* = 0.25** on the `:simulated` arm (64 items at ρ_true = 0, counter 4). The crossing lands exactly on the pre-registered boundary: `ghat(m̄)` gives ρ = −0.1443 at α = 0.125 (INSIDE the measured τ = 0.15 dead zone) and −0.2499 at α = 0.25 (first rung to clear it), so the one-rung lag behind the summary's own sign flip IS the τ dead zone, not an unexplained gap. All 8 invariants hold on all 64 REPORTED images; 0 absent patches at every rung; mask rows byte-identical across the ladder. Per-image exclusion-positive rate 0% → 27% → 86% → 98% → 100%. No gate (`P13_ALPHA_GATED = false`, 0 assertions); `consts.jl` byte-unchanged, allowance UNSPENT; seal untouched, physical anchor deferred to Phase 16
- [x] 13-14-PLAN.md — Phase-13 report: amended criteria beside originals, the D-05 trap in prose, four named limits incl. no labelled real-data validation, the real-image OOD + naming-correction honesty items, iteration ledger, deferrals, all 8 open questions closed (D-05, D-07, D-08, D-12, D-13, D-15, D-16, D-01, D-02, D-04) — `13-REPORT.md` written with the `## Scope of evidence` block and the three-arm standing table in the HEADER, before any result; every number re-read from the five persisted artifacts at report time; the epoch-4 checkpoint, the missed F5 `maxabs` bar and the DECLINED third amendment all recorded; four named limits (A–D) written in the liftable `docs/amortized.md` shape, with the anticipated prior-atom fifth entry evaluated and NOT added because the gate showed the opposite; §14 records five source discrepancies, incl. that the deep-tail band is the second-best of four rather than the best
- [x] 13-15-PLAN.md — D-15 (AMENDED) real-image ingestion from `test/test_images/`: P13_REPO_ROOT/real_tif/load_real, ghat-anchor regression, derived grid truncation, the real α-ladder, and the five D-15 testsets (ingestion, real ladder, anti-snooping, read-only, not-a-gate) + A14 falsifier — Phase-11-INDEPENDENT (D-15, D-16, D-01, D-04)
- [x] 13-16-PLAN.md — REPORTED real-image qualitative check: λ sweep with every log-BF printed beside its OOD verdict, Phase-13-vs-shipped OOD comparison, real α-ladder through the net, naming correction + target-substitution record, seal intact (D-15, D-16, D-03, D-08, D-01, D-04) — REPORTED, NOT GATED (`P13_REAL_IS_GATED = false`, zero `@test` lines, no threshold on any real-image quantity). Descriptive verdict RANDOM on both unmodified pairs at every λ rung in both read directions, coherent with the D-05 label rule (contrast ±0.1156, inside the measured τ = 0.15 dead zone). **Both fixtures flagged out-of-distribution by both detectors** — Phase-13 net 417.30 vs its own ID threshold 167.54 (2.49×), shipped reference 433.69 vs 179.14 (2.42×), reproducing the frozen constants exactly; **A12 outcome: agreement**, so the registration-aware basis did not move real microscopy back inside the training distribution. `alpha_star_real` = 0.875 (positive) / 0.625 (negative), never averaged with 13-13's `:simulated` arm. n = 2 specimens, no colocalization ground-truth label, so behaviour and never correctness. Seal intact, `consts.jl` byte-unchanged, allowance UNSPENT **[SUPERSEDED 2026-07-29 by `13-D15-AMENDMENT.md`: every figure in this bullet — `417.30 vs 167.54 (2.49×)`, `433.69 vs 179.14 (2.42×)`, `alpha_star_real 0.875 / 0.625`, contrast ±0.1156, "RANDOM on both unmodified pairs" — was measured on the pre-registered pairs `(1,2)` and `(1,3)`, BOTH of which contain `c1`, the DAPI/Hoechst nuclear counterstain (`test/runtests.jl:105`), so neither measured colocalization. The operative pair is the amended `(2,3)` green/red and the redundancy arm is DROPPED — once `c1` is excluded and the fixtures carry three channels, no second pair exists. The amended figures are measured by plan 13-17 and will be printed beside these. **No gating threshold moved and the phase's only gating arm cannot see this change** — `run_three_way_gate.jl` contains zero references to either pair constant. The record above is retained, not rewritten.]** **[AMENDED FIGURES, MEASURED BY 13-17 ON THE OPERATIVE `(2,3)` PAIR, printed beside their superseded `(1,2)` counterparts per `13-D15-AMENDMENT.md` §9: OOD **703.30 vs 167.54 (4.198×)** on `(2,3)` against **417.30 vs 167.54 (2.49×)** on `(1,2)`; shipped **948.97 vs 179.14 (5.297×)** on `(2,3)` against **433.69 vs 179.14 (2.42×)** on `(1,2)`; headline `log BF(C:R)` **+5.70311 / −8.65723** on `(2,3)` against **−0.4658 / −3.9432** on `(1,2)`, so the descriptive verdict is **COLOC (positive as sample) / RANDOM (negative as sample)** on `(2,3)` against **RANDOM / RANDOM** on `(1,2)`; `alpha_star_real` **`nothing` / `nothing`** on `(2,3)` against **0.875 / 0.625** on `(1,2)`; D-05 contrast **0.09819** on `(2,3)` against **0.11562** on `(1,2)`, both INSIDE τ = 0.15 — so on the corrected pair the label rule and the net's argmax now DISAGREE. `A12 = agreement` on both pairs. The redundancy arm no longer exists.]**
- [x] 13-17-PLAN.md — D-15 channel-pair amendment: correct the real arm from the DAPI counterstain to the c2/c3 green/red target pair, drop the redundancy arm, re-run 13-16 on the corrected substrate, and answer in the report whether the conclusion changes (D-15, D-16, D-01, D-04). **Buys no evidence**: unmasked, the corrected fixtures separate by 0.0788 against 0.0811 before — marginally *worse* — and the masked read of the same corrected pair, which separates them by 0.8576 and which `ghat.jl` itself calls the sharper separator, is FORBIDDEN because `patch_summary` applies no Otsu mask and the net was trained unmasked. It deletes an arm rather than adding one. **Mechanism RULED 2026-07-29: M1, amend in place** — `13-D15-AMENDMENT.md` §7 — so the plan carries **TWO separately disclosed changes to `spike/p13/consts.jl` in ONE edit** (§0.4). **CHANGE A** is the channel pair above. **CHANGE B** is the include-guard sentinel: `consts.jl:100` guards the whole Tier-1 body on `:P13_DEV_SEED`, which `spike/validation/p12_consts.jl:109` legitimately MIRRORS to assert seed disjointness, so a full-suite run skips the body and `runtests.jl:220` dies with 114 `UndefVarError`s (DEF-12-03). Re-pointed to `P13_DECLARED_DEVIATIONS` at **both ends plus six sibling callers**; `p12_consts.jl` is NOT touched and its mirror is correct. The two ride together because the byte-lock breaks either way: all three runners assert `h.consts_sha == p13_consts_sha()` against the 13-11 net, so the sha has to be re-derived regardless — the guard fix therefore costs **zero additional pre-registration**. The guard is re-derived as a **named, dated pre/post-amendment literal pair asserted on BOTH sides** (strictly stronger than the single equality it replaces; no widened `||` that would accept future drift), and the guard fix lands **before** the re-run so the suite can verify it. **No threshold moves; no seed literal changes; the gate result is untouched.** The **28 other poisoned include guards** the same sweep found are DOCUMENTED with remedies in `.planning/CONVENTIONS.md` C-01 and deliberately NOT fixed. **OUTCOME (2026-07-29): BOTH CHANGES LANDED IN ONE EDIT AND ARE DISCLOSED SEPARATELY** (`13-REPORT.md` §13, §14.6 CHANGE A, §14.7 CHANGE B). CHANGE B verified by targeted reproduction — **22 pass / 1 fail / 114 `UndefVarError` before, 137 pass / 0 error after**, under the exact condition the suite creates; the §5B.5 `UInt32`/`UInt64` redefinition knock-on was exercised and is a **non-event**; no seed edited; `p12_consts.jl` byte-unchanged; **DEF-12-03 closed by reference**. The `consts_sha` guard is re-derived as the named `P13_CONSTS_SHA` pre/post pair asserted on both sides at all three sites, with **no disjunction and no guard deleted**; `gate_report.jld2` / `alpha_report.jld2` byte-unchanged and the two gate-arm runners changed in **no hunk** other than that one. **The re-run answers the conclusion question: the numbers moved, one verdict moved, THE CONCLUSION IS UNCHANGED** — and the binding OOD limit got **worse** (4.198× against 2.491×), so the correction is reported as a correction and never as a rescue. **TWO CLAIMS ARE RETRACTED**: the α-ladder **no longer crosses zero** on the operative pair (`alpha_star_real` `nothing` / `nothing`), so *"negative induced μ is constructible from real microscopy pixels"* does **not** hold there; and the D-05 coherence check that previously *passed* now **disagrees**. The three `test_p13_real.jl` assertions encoding those claims are **left FAILING rather than rewritten**, which moves the suite abort to `runtests.jl:223` and masks seven Phase-13 files — all seven were therefore run individually and pass (correction retains its 2 known deliberate misses). **Whether to accept the new misses as named limits, amend the D-16 real-substrate expectation, or re-order the harness, is ESCALATED for a pre-registration ruling and NOT resolved in execution.** `P13_ITERATION_ALLOWANCE` remains **1 of 1 UNSPENT**

**Report:** `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md`
**Verdict:** The amended simulator gate (the one gating arm) **cleared all six pre-registered criteria on one run** — per-class one-vs-random AUC 0.990169 (coloc) / 0.988262 (exclusion) against the frozen floor 0.90, per-head ECE 0.0119808 / 0.012886 green against the frozen band 0.05, 0 empty bins, neither head vacuous. The D-04 iteration trigger did not fire and `P13_ITERATION_ALLOWANCE` remains 1 of 1, UNSPENT. The α-graded series and the real-image qualitative check are reported and are not gates; the report's `## Scope of evidence` header carries the standing of all three arms and the four named limits. **[AMENDED 2026-07-29, THE REAL-IMAGE CLAUSE ONLY: the real arm as reported by 13-16 read the DAPI counterstain, is SUPERSEDED by `13-D15-AMENDMENT.md`, and was re-measured on the c2/c3 pair by plan 13-17 with the redundancy arm dropped. **That re-run is DONE: the verdict moved to COLOC / RANDOM, the OOD limit worsened to 4.198× / 5.297×, and THE CONCLUSION IS UNCHANGED — qualitative, n = 2, unlabelled, OOD-bound. Two D-16 real-substrate claims are RETRACTED (§14.8).** **The gate verdict in the sentences above is untouched and cannot be reached by that amendment** — it is a simulator-ground-truth result, no threshold moved, and `P13_ITERATION_ALLOWANCE` is neither spent nor spendable on a real-image observation.]**

### Phase 14: Decision and Abstention Layer

**Goal**: Turn calibrated posteriors + the 3-way BF into an actionable batch decision {coloc / not / ABSTAIN} at a controlled Bayesian FDR, abstaining exactly when the tool should be silent
**Depends on**: Phases 11, 12, 13
**Requirements**: TBD
**Success Criteria** (what must be TRUE):

  1. `decide_coloc(...)` emits calibrated calls at a user-set Bayesian FDR across a batch, using conformal sets (ConformalPrediction.jl) + decision-risk
  2. Abstention triggers on OOD ∨ cross-method disagreement ∨ ambiguous conformal set
  3. A monotone risk-coverage curve shows abstention concentrates on hard/OOD cases

**Plans**: 14 plans in 8 waves (0-7)
> **SC1 and SC2 are AMENDED for this phase** — SC1 by D-02 (`ConformalPrediction.jl` is NOT added;
> conformal sets met in substance, hand-rolled, because the library breaks the running byte-frozen
> environment guard) and SC2 by D-05 (ASYMMETRIC fusion, not the literal `∨`: cross-method
> disagreement ALONE DECIDES and is recorded as a named field, because the literal OR would silence
> the tool exactly where Phase 9 says it should speak). Any citation of a Phase-14 result must cite
> the original wording above alongside the amendment. See `14-CONTEXT.md` D-02 and D-05.
> **Execute on the MAIN working tree — no git worktrees** (the OOD null fit needs the ~50 MB
> gitignored `spike/data/cache/p13/` pool).

Plans:
**Wave 0** *(the pre-registration freeze — must be committed before any result-producing commit)*

- [x] 14-01-PLAN.md — Wave 0: freeze the Tier-1 pre-registration (`spike/p14/consts.jl`) and assert it; opens with a blocking decision checkpoint on the five undERIVED bars

**Wave 1** *(blocked on Wave 0 completion)*

- [x] 14-02-PLAN.md — Wave 1: τ loaded with four-way provenance asserted, sha-pinned, divergence asserted (D-07)
- [x] 14-03-PLAN.md — Wave 1: `test_p14_decoupling.jl` — env/src/corpus guards plus the SC2-c, SC2-d, withdrawn-figure and forbidden-seed source greps

**Wave 2** *(blocked on Wave 1 completion)*

- [x] 14-04-PLAN.md — Wave 2: three-class posterior (anti-permutation fixture) + the running-MEAN Bayesian-FDR prefix rule
- [x] 14-05-PLAN.md — Wave 2: hand-rolled split conformal (LAC, order statistic) + the D-05 asymmetric fusion with D-06's three-valued OOD

**Wave 3** *(blocked on Wave 2 completion)*

- [x] 14-06-PLAN.md — Wave 3: `P14Result` / `P14BatchDecision` — the honesty commitments as machine-readable fields

**Wave 4** *(blocked on Wave 3 completion)*

- [x] 14-07-PLAN.md — Wave 4: `decide_coloc` — abstain-then-sort, `src/`-shaped signature, built in `spike/` and NOT shipped (D-01)
- [x] 14-08-PLAN.md — Wave 4: unstratified draw pools + the OOD density null on the Phase-13 basis

**Wave 5** *(blocked on Wave 4 completion)*

- [x] 14-09-PLAN.md — Wave 5: `run_p14_conformal.jl` (SC1-d) + the shared evaluation pool

**Wave 6** *(blocked on Wave 5 completion)*

- [ ] 14-10-PLAN.md — Wave 6: `run_p14_fdr_check.jl` (SC1-b gated, SC1-c reported-not-gated)
- [ ] 14-11-PLAN.md — Wave 6: `run_p14_riskcoverage.jl` (SC3-a/b/c)
- [ ] 14-12-PLAN.md — Wave 6: `run_p14_ood_arm.jl` (SC3-d)
- [ ] 14-13-PLAN.md — Wave 6: `run_p14_real_images.jl` (SC1-f, D-03a six-TIFF illustration + corpus record)

**Wave 7** *(blocked on Wave 6 completion)*

- [ ] 14-14-PLAN.md — Wave 7: `14-REPORT.md` + blocking human ratification of the honesty items

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
