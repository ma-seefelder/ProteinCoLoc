# Project Research Summary

**Project:** ProteinCoLoc v2.0 — AmortizedColoc
**Domain:** Amortized, calibrated simulation-based inference (NPE + NRE) for fluorescence-microscopy colocalization — a strictly decoupled Julia research spike (`spike/`) followed by conditional productionization into `src/`
**Researched:** 2026-06-26
**Confidence:** MEDIUM-HIGH

## Executive Summary

This is a **scientific-computing inference pipeline**, not a CRUD/web product. The v2.0 thesis is to replace per-dataset ADVI-minutes with **amortized millisecond inference**: a 2D 2-channel physics forward-simulator trains a neural posterior estimator (NPE, for rho_true and Delta-rho) and a neural ratio estimator (NRE, for an amortized Bayes factor) on the package's *existing* patch-correlation summary statistic, then certifies the result with Simulation-Based Calibration (SBC) and an out-of-distribution (OOD) flag. Experts in this space (post-2021 SBI literature) treat a trained network as untrustworthy until proven calibrated — so the *bundle* (calibrated amortized posterior + amortized BF + SBC/coverage + OOD), validated against the package's published KDE Bayes factor, is the publishable unit. The novelty is **"first to apply" SBI to colocalization, not "first to invent"** — every ingredient exists in cosmology/MRI/exoplanet SBI, but searches surface no prior colocalization application and none combining all four ingredients in this domain.

The research converges on a single, decisive stack decision: **NeuralEstimators.jl v0.2.1 is the one Julia-native engine that delivers both NPE (`PosteriorEstimator` + built-in `NormalisingFlow`) and the amortized Bayes factor (`RatioEstimator`/NRE)** from one API, on a **Flux v0.16 backend with CUDA-as-extension** (CPU-first, GPU optional and graceful). This de-risks the whole spike — *if the API holds*. NeuralEstimators is pre-1.0 and the Julia ML ecosystem is mid-migration Flux to Lux/Reactant, so the **AP1 smoke test is the hard gate** and a **pinned `spike/Manifest.toml` is a Definition-of-Done item**. RxInfer.jl (the open backend-migration question) is a clear **DEFER**: the non-conjugate hierarchical Student-t makes it costly to port, and it optimizes the very per-dataset axis amortized SBI is built to retire — keep it only as an optional, post-Go independent baseline cross-check, never a spike dependency.

The architecture is a **7-layer horizontal pipeline** (smoke to simulator to reused summary fns to training data to NPE/NRE to SBC/BF/OOD to demo), built strictly bottom-up with each layer an independently runnable script emitting a checkpoint. The load-bearing architectural fact is that **the existing summary statistic is the contract between simulator and network**: the simulator's only job is to emit a real `MultiChannelImage` that the *unmodified* `_prepare_data`/`patch`/`correlation` functions ingest. This keeps `src/` provably untouched, makes NPE directly comparable to ADVI (same input), and localizes `num_patches` flexibility to the post-Go productionization side (behind coexisting backends sharing one input/output contract). The **top risk is structural**: the OOD flag's detection power is bounded by the fixed summary — any misspecification orthogonal to patch-correlation (plausibly "double spillover"/"foreign PSF") is *provably undetectable*, so OOD validation must be reframed as a controlled ROC experiment with summary-orthogonal negative controls plus a posterior-predictive channel. The second risk is honesty: the chosen "tune until calibrated" policy is a data-snooping hazard unless M and the SBC threshold are pre-registered and the reported number comes from a fresh, never-tuned-against held-out SBC run.

## Key Findings

### Recommended Stack

NeuralEstimators.jl is the correct primary engine and the single decision that de-risks the spike — one package, one API, covering both amortized posterior and amortized ratio/Bayes-factor estimation on the existing summary vector. Flux is the default backend (most-documented, CPU-clean since CUDA became an extension in v0.14); Lux+Reactant is an *upgrade path*, not a competing choice. Do **not** pull in NormalizingFlows.jl or InvertibleNetworks.jl for the spike (NeuralEstimators owns the flow; note the British-spelling `NormalisingFlow` is distinct from the standalone `NormalizingFlows.jl`). See STACK.md for full detail.

**Core technologies:**
- **NeuralEstimators.jl v0.2.1** — NPE (`PosteriorEstimator`/`NormalisingFlow`) + NRE (`RatioEstimator`) + `train`/`sampleposterior`/`assess` — the only Julia-native package delivering both amortized posterior and amortized BF from one API
- **Flux.jl v0.16.x** — NN backend for the summary net + flow conditioner — most-documented backend, CPU-clean, GPU degrades gracefully (CUDA is an opt-in extension)
- **Distributions.jl + JLD2.jl + Random123.jl** — prior matched to the Turing `@model`, version/hash-guarded training-pair cache, counter-based reproducible seeding — the reproducibility-and-comparability spine
- **HypothesisTests.jl + StatsBase/CairoMakie** — KS/chi-square SBC uniformity, rank histograms, coverage/reliability plots — the calibration proof
- **Reuse existing deps** (Images, ImageFiltering for PSF/noise; KernelDensity, QuadGK for the KDE-BF baseline) — keeps the simulator and BF baseline consistent with `src/`
- **CUDA.jl (optional), Lux.jl (reserve), BenchmarkTools** — guarded GPU accelerator, throughput upgrade path, the >100x speedup measurement

**RxInfer.jl verdict:** MAYBE-feasible (needs a fixed-nu Gaussian-scale-mixture reformulation), recommend **DEFER** — largely redundant with the amortized goal; optional independent baseline only, post-Go.

### Expected Features

The audience is two-fold: the author making a Go/No-Go call, and SBI+microscopy reviewers judging credibility and novelty. See FEATURES.md.

**Must have (table stakes — reviewers reject without these):**
- Amortized NPE posterior for rho_true and Delta-rho, benchmarked vs ADVI (RMSE + interval width + >100x wall-clock)
- SBC rank histograms + KS/chi-square uniformity + coverage curve + ECE/MCE traffic-light (the calibration credibility anchor)
- Amortized log-BF (NRE/Evidence Network, l-POP loss) validated against the existing `compute_BayesFactor()`
- OOD/misspecification flag (fires off-distribution, quiet in-distribution)
- Prior consistency with the Turing `@model` ranges; held-out evaluation set; reproducible seeded `demo.jl`

**Should have (differentiators — the publishable combination):**
- First amortized SBI for colocalization, validated against a published Bayesian baseline ("first to apply")
- The trifecta bundled (calibrated posterior + amortized BF + OOD) as one shipped story
- Summary-statistic ablation as a *finding* (patch-corr vs + Manders/moments vs DeepSet), scored on per-parameter RMSE
- Honest pre-registered negative-result calibration reporting

**Defer (post-Go / v2+):**
- Productionize NPE/NRE into `src/`; user-definable `num_patches`; Turing to RxInfer evaluation
- 4-level hierarchy, 3D/Z-stacks, multi-channel, TARP joint-coverage, real wet-lab validation
- Anti-features explicitly excluded from the spike: GUI, GPU-as-requirement, sequential SBI (breaks amortization), BayesFlow on the core path, any edits to `src/`

### Architecture Approach

A horizontal-layers research pipeline built strictly bottom-up, each layer a standalone script emitting a checkpoint (plot, `.jld2`, or printed verdict), chained deterministically by `demo.jl`. The decoupling constraint (no `src/` edits) is an architectural invariant: `Pkg.develop` the parent package into the spike env and call qualified, read-only `ProteinCoLoc.correlation(...)`. Productionization (post-Go) folds validated estimators into `src/inference/` behind `colocalization_amortized()`, coexisting with the existing ADVI `colocalization()` and sharing the `_prepare_data` input + `CoLocResult` output contract. See ARCHITECTURE.md.

**Major components (spike, build order):**
1. **L0 smoke + isolated env** — toy NPE trains+samples on Flux; the hard gate before any investment
2. **L1 simulator** — `sample_prior()` + `simulate_pair(theta)::MultiChannelImage`, prior mirroring Turing priors
3. **L2 reused summary fns (unmodified)** — `patch`/`correlation`/`_prepare_data`: the simulator-network contract
4. **L3 training-data pipeline** — sample to sim to summary to standardize to version-guarded JLD2 cache (train 50-200k + held-out 5k)
5. **L4 neural inference** — MLP-first summary net (DeepSet as documented upgrade) to NPE + NRE; ADVI benchmark
6. **L5 validation siblings** — `sbc.jl`, `evidence_net.jl`, `ood.jl` (independent, share trained nets + sampler)
7. **L6 demo + Go/No-Go memo** — deterministic end-to-end from a fixed seed

### Critical Pitfalls

Top risks from PITFALLS.md, in priority order:

1. **OOD flag's structural blind spot (TOP RISK)** — a fixed patch-correlation summary can only detect misspecifications that *move* that summary; "double spillover"/"foreign PSF" may be summary-orthogonal and provably undetectable. *Avoid:* reframe OOD validation as a ROC experiment with positive AND summary-orthogonal negative controls, add a posterior-predictive-check channel, report blind spots honestly rather than shipping a binary "it fired."
2. **"Tune until calibrated" becomes data-snooping** — repeatedly tuning against the same SBC draws turns the KS p-value into a selection statistic that certifies nothing. *Avoid:* pre-register M and the threshold before looking at ranks; tune on one seed but report the number from a **fresh, never-tuned-against, independently-seeded SBC run**; log every tuning iteration.
3. **Summary-statistic insufficiency** — patch-corr was built for ADVI, not as a sufficient statistic under the full nuisance parameterization; SBC can *still pass* on a posterior over a lossy summary. *Avoid:* make the AP4 ablation diagnostic, not decorative — score per-parameter RMSE + interval width vs ADVI; SBC-pass + high-RMSE is the insufficiency smoking gun.
4. **NeuralEstimators/Flux API churn (pre-1.0, Flux to Lux/Reactant migration)** — a `Pkg.update` can silently break the spike mid-flight. *Avoid:* AP1 smoke gate + committed pinned `spike/Manifest.toml`; keep NormalizingFlows.jl / BayesFlow fallback notes warm but off the critical path.
5. **SBC-vs-real conflation + NRE-BF miscalibration/prior sensitivity** — SBC certifies calibration *under the simulator* only; quoting it as a real-data guarantee is the central honesty trap. The amortized BF is intrinsically prior-dependent and BCE-trained ratios are overconfident. *Avoid:* always pair SBC with OOD and state "under the simulator"; use l-POP loss (not vanilla BCE), pin identical Delta-rho + prior vs the KDE baseline, consider BNRE if ratios look overconfident.

Additional named traps: standardization/OOD-covariance leakage across train/held-out/SBC splits (fit on train split only), Windows GPU/XLA breakage (CPU-only default, guarded fall-through), and regenerating training data every run (version/hash-keyed cache).

## Implications for Roadmap

The spike's AP1-AP7 / M0-M5 structure already encodes a defensible phase plan; research strongly endorses it and sharpens the gates. Suggested phase structure:

### Phase 1: Environment + Smoke Gate (AP1 / M0)
**Rationale:** The single biggest unknown is whether the pre-1.0 NeuralEstimators flow API converges on Flux-on-Windows. Research is unanimous this is a *hard gate before any other investment* — building the simulator first risks weeks on a broken foundation.
**Delivers:** Isolated `spike/Project.toml` + **committed pinned `spike/Manifest.toml`**; green 20-line toy NPE (`PosteriorEstimator` + `NormalisingFlow`, `train`+`sampleposterior`); CPU-only, no forced CUDA import.
**Addresses:** Isolated env + green smoke test (table stakes), reproducibility.
**Avoids:** Pitfall 4 (API churn) and Pitfall 5 (Windows GPU breakage). DoD = pinned Manifest + smoke green on *this* machine.

### Phase 2: Forward Simulator + Summary Contract (AP2 / M1)
**Rationale:** Everything downstream consumes the simulator; prior drift here silently invalidates the ADVI benchmark and the BF validation. Build only after the smoke is green.
**Delivers:** `simulate_pair(theta) -> MultiChannelImage` (spillover, autofluorescence, label efficiency, PSF, sub-pixel shift, noise) that the *unmodified* `_prepare_data`/`patch`/`correlation` ingest; `sample_prior()` referencing Turing prior constants copied into `spike/NOTES.md` on Day 1.
**Uses:** Distributions.jl, ImageFiltering.jl; `Pkg.develop`-ed parent package (read-only).
**Implements:** L1 simulator + L2 reused-summary contract.
**Avoids:** Anti-Pattern 3 (prior drift), Anti-Pattern 1 (editing `src/`).

### Phase 3: Training-Data Pipeline (AP3)
**Rationale:** Caching and split discipline must exist before iterating on the network; this phase sets the leakage-prevention discipline structurally.
**Delivers:** sample to sim to summary to standardize to version/hash-guarded JLD2 cache; train 50-200k + held-out 5k; deterministic under Random123 seed.
**Uses:** JLD2, Random123, StatsBase; BayesInteractomics cache-guard template.
**Avoids:** Pitfall 8 (standardization/OOD-covariance leakage — fit scalers on train split only) and the regenerate-every-run performance trap.

### Phase 4: NPE Training + ADVI Benchmark + Ablation (AP4 / M2)
**Rationale:** The headline capability and the layer most likely to need a second iteration (MLP to DeepSet). The ablation lives here because summary sufficiency also drives OOD detection power.
**Delivers:** NPE for rho_true + Delta-rho; RMSE + interval width + >100x wall-clock vs real `colocalization()` ADVI; summary-statistic ablation (patch-corr vs + Manders/moments vs DeepSet) scored on **per-parameter RMSE**.
**Uses:** NeuralEstimators `PosteriorEstimator`, Flux, BenchmarkTools.
**Avoids:** Pitfall 3 (insufficiency) — make the ablation a gating diagnostic; confirm flow convergence before blaming the summary.

### Phase 5: Validation Bundle — SBC + Amortized BF + OOD (AP5-AP6 / M3-M4)
**Rationale:** The three validation tracks are architectural siblings — all consume the trained nets + sampler, none depend on each other — and together form the publishable bundle. Group them, but treat each gate seriously.
**Delivers:** SBC rank histograms + KS/chi-square + coverage + ECE/MCE traffic-light (pre-registered M + threshold; fresh held-out confirmatory run); amortized log-BF (l-POP loss) validated vs `compute_BayesFactor()` on identical Delta-rho + prior; OOD flag with ROC over a misspecification grid + summary-orthogonal negative controls + PPC channel.
**Uses:** HypothesisTests, KernelDensity/QuadGK (KDE-BF baseline), `RatioEstimator`.
**Avoids:** Pitfalls 1, 2, 6, 7 — the cluster of honesty/credibility risks. This is the highest-scrutiny phase.

### Phase 6: Reproducible Demo + Go/No-Go Memo (AP7 / M5)
**Rationale:** Closes the spike with a falsifiable, reproducible verdict.
**Delivers:** Seeded `demo.jl` chaining all layers CPU-only; metrics table; Go/No-Go memo framing SBC as "calibrated under the simulator," paired with the OOD result; full-build-out decision. Verify `git status` on `src/` is clean (decoupling proof).
**Avoids:** Pitfall 6 (SBC-vs-real conflation) in the memo language.

### Phase 7 (conditional on Go): Productionization
**Rationale:** Depends only on a Go memo; the spike is promoted (copied with provenance), never imported live.
**Delivers:** `colocalization_amortized()` in `src/inference/`, coexisting with ADVI behind the shared `_prepare_data`/`CoLocResult` contract; estimator registry keyed by `num_patches` (where the grid flexibility lives); consider a package extension/weakdep so the classic path loads without the NeuralEstimators tree.

### Phase Ordering Rationale

- **Strict bottom-up dependency chain:** smoke gates everything (de-risks the pre-1.0 API), simulator+prior-consistency gates the benchmark and BF validation, training cache gates the network, the network gates all three validation tracks. This is dictated by the dependency graph in FEATURES.md and the build-order table in ARCHITECTURE.md.
- **The validation trifecta is grouped** because the three scripts are mutually independent but share the trained nets and the simulate-infer harness — build the harness once, fan out.
- **The ablation sits with NPE (Phase 4), not validation,** because summary sufficiency couples to OOD detection power; deciding the summary late would force re-running SBC/BF/OOD.
- **Productionization is firewalled** behind the Go gate so the dependency-heavy Flux tree never leaks into the published manifest pre-decision.

### Research Flags

Phases likely needing deeper research during planning (`/gsd:plan-phase --research-phase`):
- **Phase 1 (AP1):** Exact NeuralEstimators v0.2.1 constructor/flow signatures must be verified live against `dev` docs (pre-1.0, churning); confirm the CPU-only Flux path is actually exercised by the package's own tests, not just the Reactant examples.
- **Phase 5 (AP5-AP6):** The OOD ROC experimental design (which axes are expected blind spots, how to construct summary-orthogonal perturbations) and the NRE loss choice (l-POP vs BNRE) are the highest-uncertainty design decisions and the highest reputational stakes.

Phases with standard/well-documented patterns (can skip research-phase):
- **Phase 3 (AP3):** JLD2 caching with version/hash guard is a directly reusable BayesInteractomics pattern.
- **Phase 6 (AP7):** Deterministic seeded orchestration is mechanical once the layers exist.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | MEDIUM-HIGH | Core SBI stack (NeuralEstimators/Flux versions) verified against official docs + GitHub releases (HIGH); the pre-1.0 API surface and the CPU-only Flux path are intentionally gated behind the AP1 smoke test (MEDIUM). RxInfer Student-t feasibility verified by mechanism, no published colocalization-shaped example (MEDIUM). |
| Features | HIGH | SBI table-stakes grounded in multiple recent peer-reviewed sources + method docs. The colocalization-novelty claim rests on *absence* of prior art in searches (negative evidence, MEDIUM). |
| Architecture | MEDIUM-HIGH | Component structure, data flow, and decoupling anchored in existing code + two concrete BayesInteractomics reference engines (HIGH); exact NeuralEstimators API surface deliberately re-verified at the smoke gate (MEDIUM). |
| Pitfalls | HIGH | SBI/SBC/misspecification methodology from arXiv primary sources (Schmitt, Talts, Delaunoy, Jeffrey & Wandelt). NeuralEstimators API specifics and Windows/Flux specifics are MEDIUM (pre-1.0, no project-specific repro). |

**Overall confidence:** MEDIUM-HIGH

### Gaps to Address

- **NeuralEstimators v0.2.1 exact API + CPU-Flux viability** — resolve in Phase 1 via the smoke test; do not write Phase 2+ until green and the Manifest is pinned. This is the single largest uncertainty and is structurally gated.
- **Summary-statistic sufficiency** — unknown until measured; resolve in Phase 4 via the ablation scored on per-parameter RMSE. Plan for a possible MLP to DeepSet + moments iteration.
- **OOD detection blind-spot structure** — which misspecification axes are summary-orthogonal is a research question, not a known; resolve in Phase 5 by constructing explicit negative controls and reporting ROC + named blind spots rather than a binary flag.
- **SBC pass criterion under "tune until calibrated"** — pre-register M and threshold during Phase 5 planning (a binding decision, not a runtime choice) to neutralize the data-snooping hazard.
- **RxInfer complementary-vs-redundant** — left as a post-Go, time-boxed experiment using the fixed-nu Gaussian-scale-mixture formulation; do not let it consume spike budget.

## Sources

### Primary (HIGH confidence)
- NeuralEstimators.jl official docs + GitHub releases (v0.2.1, 2026-04-07) — PosteriorEstimator/RatioEstimator/NormalisingFlow/train/sampleposterior/assess; Flux/Lux/SimpleChains/Reactant backends — https://msainsburydale.github.io/NeuralEstimators.jl/dev/
- Flux.jl docs — v0.16.x, CUDA-as-extension since v0.14 — http://fluxml.ai/Flux.jl/dev/guide/gpu/
- Talts et al. (2018), Validating Bayesian Inference Algorithms with SBC — rank histograms / uniformity — https://arxiv.org/pdf/1804.06788
- Schmitt et al. (2024), Detecting Model Misspecification in Amortized Bayesian Inference — summary-space MMD/Mahalanobis, detection-vs-inference trade-off, undetectable summary-sufficient misspecifications — https://arxiv.org/html/2406.03154v2
- Jeffrey & Wandelt (2024), Evidence Networks (l-POP-Exp loss) — fast amortized neural model comparison — https://arxiv.org/abs/2305.11241
- Delaunoy et al. (2022), Balanced Neural Ratio Estimation (BNRE) — NRE overconfidence, conservative loss — https://arxiv.org/abs/2208.13624
- Hermans et al. (2021), A Trust Crisis in SBI? — basis for the calibration-mandatory framing — https://arxiv.org/abs/2110.06581
- Existing code + BayesInteractomics reference engines — src/colocalization.jl, src/bayes.jl, src/LoadImages.jl; _bin_calibration/CalibrationResult, parametric simulation + JLD2 cache guard — read directly
- Project planning — .planning/PROJECT.md, .planning/10_plan_amortizedcoloc.md

### Secondary (MEDIUM confidence)
- RxInfer.jl docs (v1.17.0) + non-conjugate manual — ExponentialFamilyProjection runtime/accuracy warnings; basis for the DEFER verdict — https://docs.rxinfer.com/stable/
- Lux.jl GPU management; InvertibleNetworks.jl (JOSS) — reserve upgrade paths for GPU/image-sized flows
- Zaheer et al., Deep Sets — permutation-invariant summary networks — https://www.inference.vc/deepsets-modeling-permutation-invariance/
- Sailynoja et al. (2022), ECDF-over-histogram for SBC; State of the SciML Ecosystem 2025 (Flux to Lux migration signal)
- Colocalization domain references (Pearson/Manders/Costes, Coloc2, JACoP) — https://imagej.net/imaging/colocalization-analysis

### Tertiary (LOW confidence)
- Absence-of-prior-art for SBI in fluorescence colocalization — negative evidence underpinning the "first to apply" novelty framing; not a positive citation, validate during write-up
- Gaussian-scale-mixture Student-t reformulation as the recommended RxInfer route — standard result, applied judgment, not project-verified

---
*Research completed: 2026-06-26*
*Ready for roadmap: yes*
