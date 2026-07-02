---
gsd_state_version: 1.0
milestone: v2.0
milestone_name: milestone
status: planning
stopped_at: Phase 10 context gathered
last_updated: "2026-07-02T08:14:09.843Z"
last_activity: 2026-07-01
progress:
  total_phases: 16
  completed_phases: 4
  total_plans: 20
  completed_plans: 20
  percent: 25
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-06-26)

**Core value:** A trained NPE/NRE produces calibrated, amortized colocalization inference (posterior + Bayes factor) in a single forward pass, >100x faster than per-dataset ADVI, with a demonstrated SBC/coverage proof and an honest OOD flag.
**Current focus:** Phase 5 — validation bundle (sbc + amortized bf + ood)

## Current Position

Phase: 5
Plan: Not started
Status: Ready to plan
Last activity: 2026-07-01

Progress: [█████████░] 92%

## Performance Metrics

**Velocity:**

- Total plans completed: 20
- Average duration: - min
- Total execution time: 0.0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 4 | - | - |
| 02 | 4 | - | - |
| 03 | 5 | - | - |
| 04 | 7 | - | - |

**Recent Trend:**

- Last 5 plans: -
- Trend: -

*Updated after each plan completion*
| Phase 01 P02 | 22min | 2 tasks | 4 files |
| Phase 01 P03 | 9 | 2 tasks | 2 files |
| Phase 01 P04 | 13min | 2 tasks | 1 files |
| Phase 02 P01 | 15min | 2 tasks | 5 files |
| Phase 02 P02 | 11min | 2 tasks | 4 files |
| Phase 02 P03 | 32min | 2 tasks | 5 files |
| Phase 02 P04 | 22min | 2 tasks | 5 files |
| Phase 03 P01 | 8min | 2 tasks | 4 files |
| Phase 03 P02 | 12min | 3 tasks | 5 files |
| Phase 03 P03 | 26min | 3 tasks | 6 files |
| Phase 03 PP04 | 18min | 2 tasks tasks | 2 files files |

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 7-layer horizontal build order; ENV smoke gate (Phase 1) is the hard gate before all downstream investment (pre-1.0 NeuralEstimators API is the top unknown)
- [Roadmap]: SBC/BF/OOD grouped into one Phase 5 — sibling plans sharing one θ*~π→simulate→infer harness over the trained nets
- [Roadmap]: Phase 7 (productionization) is conditional on a Go decision in the Phase 6 memo and is the ONLY phase that edits `src/`; spike phases (1–6) keep `src/` provably untouched
- [Phase 01-02]: Spike env is minimal+isolated (NeuralEstimators 0.2.1, Flux 0.16.10, Distributions 0.25.128, no CUDA); CPU-only NPE smoke recovers theta within tol=0.3 under a re-runnable stdlib Test gate. CUDA-absence asserted via Pkg.dependencies (installed set) not a Manifest regex, since Flux/NNlib/Zygote/NeuralEstimators declare CUDA as inert weakdep extensions.
- [Phase ?]: [Phase 01-03]: Spike env frozen as reproducibility artifact — Julia 1.12.6 pinned via juliaup directory override (enforcing) + spike/.julia-version (doc-only); Manifest pins NeuralEstimators 0.2.1 + Flux 0.16.10, no top-level CUDA. NeuralEstimators default / BayesFlow fallback recorded in spike/NOTES.md (ENV-03, ENV-04).
- [Phase ?]: [Phase 01-04]: include() fallback chosen over Pkg.develop for ENV-01 — the parent's GLMakie/Makie 0.21 + Turing 0.44 + GraphNeuralNetworks 1.1.0 tree caps NeuralEstimators <0.2.1, silently downgrading the pinned 0.2.1 to 0.1.4 and breaking the v0.2.1-API smoke; later phases reach src/ read-only via include('../src/<file>.jl'). Package boundary deferred to Phase 4 (D-02). Decoupling proven byte-identical to baseline f581d95.
- [Phase ?]: SIM-03 summary contract proven on a synthetic strictly-positive image via UNCHANGED src/ patch()/correlation() at fixed 8x8 (D-10) before any physics exists
- [Phase ?]: CairoMakie/Images/ImageFiltering/HypothesisTests/StatsBase co-resolve cleanly with NeuralEstimators 0.2.1 (no downgrade, no CUDA); resolve-risk gate automated in runtests.jl
- [Phase ?]: Phase 2 simulator: shared-latent standardized-smooth-field generator (D-15) with sign(rho) flip drives monotone, sign-correct induced patch-correlation (-0.43 to +0.60)
- [Phase ?]: SIM-02 (02-03): frozen monotone ĝ (isotonic PAVA + clamped piecewise-linear inverse) maps μ↦ρ_true; induced μ matches Turing μ-prior W1=0.052<0.10 over realized range [-0.68,0.847]; negative μ tail prior-only (real anchor pos=0.33/neg=0.25)
- [Phase ?]: SIM-04 Spearman monotonicity threshold fixed at 0.95 (achieved 1.0); perturbation effects asserted as paired shared-seed differences with OneSampleTTest p<0.05 (T-02-GATE)
- [Phase ?]: CairoMakie plausibility figures are spike-local headless (D-13); .gitignore scoped negation tracks spike/figures/plausibility.png as the evidence deliverable
- [Phase ?]: [Phase 03-01]: JLD2 (cache backend, D-04) + Random123 (counter-based seeding, D-11) added to isolated spike env; JLD2 already transitive (promoted to direct), only Random123 v1.7.1 + RandomNumbers v1.6.0 newly installed; NeuralEstimators stays pinned v0.2.1 (no co-resolve downgrade); resolve-risk gate extended to cover both new deps; five-SC MISSING scaffold (test_data_pipeline.jl) wired into single runtests.jl gate
- [Phase ?]: [Phase 03-02]: DATA-01 generation core — encode_d01 (128-dim: 64 imputed corr + 64 binary mask, fully-missing kept per D-13), encode_aug (AUG_DIM=142, 14 moments, D-02), Philox4x per-sample keyed RNG keyed by (master_seed,idx) with disjoint HOLDOUT/FOLD salts (D-10/D-11), cost-aware imsize sampler (>=1024^2 capped 10%, E[cost]~4.68x, D-03); generate_samples parallel==serial byte-identical across -t 1 and -t 4 (D-12)
- [Phase ?]: [Phase 03-03]: DATA-02 sharded JLD2 cache — SHA-256 content hash over data-defining source bytes + canonical(config) names the cache dir and lives in meta.jld2 (D-05); atomic .tmp+integrity-check+mv shard writes with shard_done resume-by-skip (D-04/D-06); reserved >=20 ADVI holdout from disjoint XOR-salted holdout_rng into separate holdout.jld2 with negative global_index, structurally disjoint from main pool (D-10)
- [Phase ?]: [Phase 03-04]: DATA-03 leak-free loader — module Loader is the SOLE standardization path; ZScoreTransform fit on TRAIN columns only (D-07), no standardize_all symbol exists so leakage is impossible by construction (D-08), deterministic k=5 folds from fold_rng(master_seed XOR FOLD_SALT) (D-09), mask rows 65:128 bypass, holdout excluded structurally (D-10)

### Roadmap Evolution

- v2.0 feature expansion appended to the current milestone (Phases 8–16), derived from `next_project_analysis/11_plan_proteincoloc_v2.md` + scope decision log. Consolidated as ONE v2.0 paper (not split); multiplex/copula deferred to v2.1; expert-concordance dropped.
- Phase 8 added: External Physical Ground-Truth Corpus (Depends on: nothing — parallelizable now)
- Phase 9 added: Cross-Method Comparator Harness (Depends on: nothing — parallelizable now)
- Phase 10 added: Manuscript Skeleton + Related-Work Positioning (Depends on: nothing — parallelizable now)
- Phase 11 added: Registration + Chromatic Uncertainty as Latent (Depends on: Phase 7)
- Phase 12 added: Spatial Colocalization Map (GP/CAR) (Depends on: Phases 7, 11 — serialized behind 11 due to shared PosteriorEstimator/training code; descope-to-v2.1 candidate)
- Phase 13 added: Three-Hypothesis Amortized Bayes Factor (Depends on: Phase 7)
- Phase 14 added: Decision + Abstention Layer (Depends on: Phases 11, 12, 13)
- Phase 15 added: Calibration Operating Envelope + CI Gate (Depends on: Phases 11, 12)
- Phase 16 added: External Validation + Manuscript Assembly (Depends on: Phases 8, 9, 10, 14, 15)
- Dependency analysis (analyze-dependencies): sole correction vs. draft was serializing Phase 12 behind Phase 11 (file overlap on the estimator); Waves — A {8,9,10} now, B {11 ∥ 13, then 12}, C {14,15}, D {16}.

### Pending Todos

[From .planning/todos/pending/ — ideas captured during sessions]

None yet.

### Blockers/Concerns

[Issues that affect future work]

- [Phase 1]: NeuralEstimators v0.2.x is pre-1.0 and the Julia ML ecosystem is mid Flux→Lux/Reactant migration — the CPU-only Flux path must be confirmed exercised at the smoke gate; pinned Manifest is the mitigation
- [Phase 5]: OOD detection power is structurally bounded by the fixed patch-correlation summary — summary-orthogonal misspecifications are provably undetectable and must be named, not hidden
- [Phase 5]: "Tune until calibrated" is a data-snooping hazard — M and SBC threshold must be pre-registered before Phase 5 planning runs

## Deferred Items

Items acknowledged and carried forward from previous milestone close:

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| Backend (v2) | BACK-01 Turing→RxInfer.jl feasibility evaluation | Deferred (post-Go, never a spike dependency) | 2026-06-26 |
| Summary net (v2) | BACK-02 DeepSet permutation-invariant summary upgrade | Deferred (upgrade if MLP summary proves insufficient) | 2026-06-26 |

## Session Continuity

Last session: 2026-07-02T08:14:09.836Z
Stopped at: Phase 10 context gathered
Resume file: .planning/phases/10-manuscript-skeleton-and-related-work-positioning/10-CONTEXT.md
