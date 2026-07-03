---
gsd_state_version: 1.0
milestone: v2.0
milestone_name: milestone
status: planning
stopped_at: Phase 7 context gathered
last_updated: "2026-07-03T12:10:15.371Z"
last_activity: 2026-07-03
progress:
  total_phases: 16
  completed_phases: 8
  total_plans: 37
  completed_plans: 37
  percent: 50
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-06-26)

**Core value:** A trained NPE/NRE produces calibrated, amortized colocalization inference (posterior + Bayes factor) in a single forward pass, >100x faster than per-dataset ADVI, with a demonstrated SBC/coverage proof and an honest OOD flag.
**Current focus:** Phase 08 — external physical ground truth corpus

## Current Position

Phase: 08
Plan: Not started
Status: Ready to plan
Last activity: 2026-07-03

Progress: [██████████] 100%

## Performance Metrics

**Velocity:**

- Total plans completed: 22
- Average duration: - min
- Total execution time: 0.0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 4 | - | - |
| 02 | 4 | - | - |
| 03 | 5 | - | - |
| 04 | 7 | - | - |
| 06 | 2 | - | - |

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
| Phase 05 P01 | 22min | 3 tasks | 8 files |
| Phase 05 P05-02 | 19min | 2 tasks | 7 files |
| Phase 05 P05-03 | 70min | 2 tasks | 2 files |
| Phase 05 P05-04 | 60min | 3 tasks | 4 files |
| Phase 06 P01 | 40 | 3 tasks | 1 files |
| Phase 06 P02 | 20min | 2 tasks | 2 files |

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
- [Phase 05-01]: Phase-5 anti-snooping contract locked in committed spike/validation/consts.jl BEFORE any reported run — M=2000/L=999/bins=50, all SBC/BF/OOD thresholds, VAL_MASTER_SEED=0x5BC0FFEE disjoint from NPE_MASTER_SEED=0xC0FFEE (via VAL_SALT); fixture gates use a separate VAL_FIX_SEED so they never consume the reported stream (D-02/SBC-03)
- [Phase 05-01]: One shared harness.jl (draw_simulate_infer + paired Δρ path) composes the frozen Phase-4 read surface over the loaded-once net; sbc.jl delivers 7θ+Δρ ranks, KS/χ² uniformity, coverage curve, and ported _bin_calibration ECE/MCE traffic-light; figures.jl owns ROC/BF stubs so Wave-2 plans call, never edit it (SBC-01..04)
- [Phase 05-02]: amortized log-BF = logratio(m=1)-logratio(m=0)-measured_log_prior_odds via NRE model-comparison (v0.2.1 has no l-POP; honest relabel)
- [Phase 05-02]: KDE compute_BayesFactor math ported verbatim to isolated spike/baseline/ env; shared spike/ env untouched (NeuralEstimators 0.2.1)
- [Phase 05-02]: fixture reproduction corr=0.961 (D-08a PASS) but max|dlogBF|=3.18 (D-08b NOT met at fixture scale) - reported honestly, consts not tuned; reported gate is 05-04
- [Phase 05]: [Phase 05-03]: two-channel OOD flag (train-only Mahalanobis density + posterior-predictive per-feature discrepancy, OR-fused); ID-quantile operating point committed on a HELD-OUT ID pool (out-of-sample honest ~5% FPR); hand-rolled ROC/AUC, no new package
- [Phase 05]: [Phase 05-03]: fixed-summary OOD power structurally bounded — TWO measured named blind spots: D-04 correlation-preserving transforms (affine exact, block-permute residual) AND detector-noise mismatch (density+PP AUC~0); optics/PSF robustly detected; full four-family ROC is 05-04 reported job
- [Phase ?]: [Phase 06-01]: spike/demo.jl two-tier seeded runner — fast tier re-exercises NPE+BF/NRE+OOD at fixture scale on VAL_FIX_SEED with per-layer twin-run bit-reproducibility asserts (global RNG seeded alongside Philox since sampleposterior threads no rng); step-0 git-status assertion certifies src/ untouched; --full spawns reported gates as independent subprocesses; locked files byte-unchanged
- [Phase ?]: [Phase 06-02]: Go/No-Go memo renders a Clean Go — literal pre-registered gates did NOT pass as written, Go justified by reading each residual failure to a non-method cause (chi2-over-power at M=2000, residual data-scale gap, clamped-KDE tail artifact); credibility backed by the D-05 independent fresh-seed re-pre-registered confirmation ship-gate inside Phase 7
- [Phase ?]: [Phase 06-02]: BF max|d logBF|=10.9535 quoted after runtime-confirming against demo.jl printed bf_report max_abs_err (A1 resolved); memo reports BOTH number sets — Set 1 pre-registered all-FAIL (0.8861/0.730/2.81e-40) verbatim from 05-04-SUMMARY.md, Set 2 post-hoc (0.0164/0.9358/AUC 1.0) with data-snooping exposure + 3 mitigations (consts.jl vs e9c91d3, DEV seed 0xDE7C0DE, single confirmatory VAL run)
- [Phase ?]: [Phase 06-02]: demo.jl machine-gates DEMO-03 — assert isfile(MEMO_PATH) + occursin content asserts (Clean Go, falsification, 0.8861, 0.9358, ship-gate, Phase 11, Phase 13); SC3/DEMO-03 table rows now report real PASS

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
- [Phase 5-04]: ALL THREE reported pre-registered gates FAIL at locked consts on fresh VAL_MASTER_SEED (honest negatives, consts.jl untouched): SBC every-parameter KS/chi2 p~0 (overconfident/miscalibrated NPE); BF corr=0.886<0.95 AND max|dlogBF|=3.60>0.5; OOD density pooled AUC=0.730<0.80 (noise family AUC=0 blind spot). Feeds Phase-6 Go/No-Go as No-Go-or-iterate.
- [Phase 5-iter1 (2026-07-02)]: BOUNDED calibration iteration retrained the NPE on the SAME 50k cache with a higher-capacity flow (dstar 32→64, summary 2×128→3×256, coupling 6→10) + stabler recipe (LR 5e-4→2.5e-4, batch 64→128, 300 epochs/patience 40). Selected on a DISJOINT dev seed (0xDE7C0DE); consts.jl BYTE-UNCHANGED. Calibration improved dramatically — DEV SBC ρ_true KS D 0.158→0.038, Δρ 0.155→0.017; the frozen retrain overwrites trained_npe.jld2. CONFIRMATORY (VAL_MASTER_SEED, run once): SBC still FAIL overall but ALL 8 params now GREEN on ECE (ρ_true ECE=0.016) and 3/8 (spillover/autofluorescence/noise) pass KS — the strict all-8 KS+50bin-χ² conjunction rejects (label_efficiency/shift genuinely non-uniform, χ² hyper-sensitive at M=2000). BF FAIL corr=0.918 (↑ from 0.886) but max|dlogBF|=14.4 (↑ worse — sharper NPE Δρ tails inflate the KDE-baseline gap). OOD FAIL unchanged density pooled AUC=0.730 (Mahalanobis channel is NPE-independent; noise-family blind spot structural). Verdict: PARTIAL. Next levers: scale data toward 200k (SBC), retrain NRE (BF).
- [Phase 5-iter1]: OOD PP channel now VIABLE — finite-guard added in ood.jl (_theta_tuple maps non-finite posterior θ̂ to in-range fallbacks) so strong-misspec re-simulation no longer crashes; PP θ̂ finite on all 4 families. Still excluded from the reported OR-fusion (with_pp=false unchanged in run_ood.jl); available to re-enable in Phase-6/7.
- [Phase 5-iter2 (2026-07-02)]: TASK 1 — retrained the NRE on the SAME 50k (θ,summary) cache the NPE used (assemble_ratio_data_cache; 16× the prior 3k fresh-sim pairs, leak-free: cache gen-seed 0x134d8f3 disjoint from VAL/NPE seeds), higher capacity (num_summaries 32→64, conditioner 64→3×256), stabler recipe (LR 2.5e-4, batch 128, 300ep/pat40), + research-A7 difference encoding (pair_encode = concat + 64-dim continuous-row correlation CONTRAST, input dim 256→320; shared with bf.jl). Develop DEV 0xDE7C0DE (corr 0.862→0.903); CONFIRMATORY VAL_MASTER_SEED once (run_bf.jl unchanged): corr 0.918→0.936 (↑), max|Δ logBF| 14.4→10.95 (↓). BF gate still FAILS. VERDICT on max|Δ|: it is a KDE-BASELINE TAIL/CLAMP ARTIFACT, not an NRE deficiency — every large-|Δ| sweep point (|Δρ|≳0.4) is where the sharpened NPE Δρ posterior is fully one-sided, forcing the KDE baseline's P(Δρ>0) to its _clampp=1e-8 floor → logBF=±18.41, while the amortized NRE (correctly) stays bounded ~±8; mid-range (|Δρ|<0.4) agrees within ~1-2. D-08b (max|Δ|≤0.5) is structurally unclearable against a clamped-KDE baseline at the sweep tails — a pre-registration nuance, consts.jl untouched. BF fast test 17/17.
- [Phase 5-iter2]: TASK 2 — added an AUXILIARY image-noise OOD channel (ood.jl Channel 3: noise_features/fit_noise_null/noise_score) OR-fused with the density channel, closing the detector-noise blind spot. 10 scale/rotation/permutation-INVARIANT features (HF energy ratio, robust HF scale ratio, outlier fraction, HF excess kurtosis, lag-1 autocorr; ∇² via finite difference, NO FFT/new dep) designed from the SMOOTH training distribution (SC5, not the held-out test set). Null + robust-z fit on a fresh TRAIN-ONLY ID pool; reported detector = per-sample max of density/noise robust-z (continuous OR, D-05). maha_auc/combined_auc now report the FUSED detector; density maha_thr unchanged so run_ood's density-only neg-control flag is byte-identical. Develop DEV (noise-family AUC 0.0→1.0); CONFIRMATORY VAL once (run_ood.jl unchanged): OOD gate PASSES — all 4 families best-AUC=1.0 (noise 1.0 at every level), combined pooled AUC=1.0, ID fire-rate 0.05, neg-controls KS-invariant + density-quiet. OOD fast test 27/27. Updated blind-spot framing: affine+rotate stay quiet on ALL channels (true blind spots); block-permute stays density-quiet (frozen gate passes) but the noise channel correctly FIRES on its tile-seam HF artifacts (~0.93) — so block is no longer summary-orthogonal to the fused detector. Remaining blind spot is NARROWER: orthogonal to BOTH the 8×8 correlation summary AND the image-noise features (e.g. a pure positive-affine rescale, still exactly invariant). consts.jl BYTE-UNCHANGED (verified vs e9c91d3).

## Deferred Items

Items acknowledged and carried forward from previous milestone close:

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| Backend (v2) | BACK-01 Turing→RxInfer.jl feasibility evaluation | Deferred (post-Go, never a spike dependency) | 2026-06-26 |
| Summary net (v2) | BACK-02 DeepSet permutation-invariant summary upgrade | Deferred (upgrade if MLP summary proves insufficient) | 2026-06-26 |

## Session Continuity

Last session: 2026-07-03T12:10:15.361Z
Stopped at: Phase 7 context gathered
Resume file: .planning/phases/07-productionization-conditional-on-go/07-CONTEXT.md
