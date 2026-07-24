---
gsd_state_version: 1.0
milestone: v2.0
milestone_name: milestone
status: in-progress
stopped_at: GO decision taken (2026-07-24, Option A — GO with named limits, 07-GO-NO-GO-UPDATE.md). v2.0 ships/publishes as a calibrated coloc tool: TARGETS ρ_true/Δρ SBC-calibrated (randomized ranks), BF simulation-based AUC 0.994 (spike 014, KDE baseline replaced), OOD pass; nuisance marginal drift + twice-amended gate + atom handling documented as named limits. No further training/gate iteration. Option B (BF §6 gate integration + summary redesign) deferred.
last_updated: "2026-07-24T18:00:00.000Z"
last_activity: 2026-07-24 -- spikes 006-014 + post-hoc re-analysis + updated Go/No-Go; user chose GO with named limits (Option A)
progress:
  total_phases: 16
  completed_phases: 9
  total_plans: 53
  completed_plans: 51
  percent: 56
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-06-26)

**Core value:** A trained NPE/NRE produces calibrated, amortized colocalization inference (posterior + Bayes factor) in a single forward pass, >100x faster than per-dataset ADVI, with a demonstrated SBC/coverage proof and an honest OOD flag.
**Current focus (2026-07-24):** Phase 07 productionization is PAUSED at 9/11 plans (07-09 32×32 training never authorized). A calibration investigation (spikes 006-013) diagnosed the amended grid-8 ship-gate FAIL: coloc TARGETS ρ_true/Δρ are calibrated; ρ_true SBC-fail = prior-atom artifact (test-fixable via randomized ranks); nuisance fails = 0.06-SD marginal drift at the M=2000 over-power edge (test-fixable via a nuisance equivalence rule, spec draft `07-NUISANCE-SBC-SPEC-DRAFT.md`); BF attrition = KDE-baseline defect (Memo §5 incomplete). Neither model lever (capacity/truncation) earns the §6.4 retrain. **Open decision:** final bundled amendment on a retrained model vs Go/No-Go re-eval publishing with named limits. Findings in `Skill("spike-findings-proteincoloc")`.

## Current Position

**ACTIVE phase: 07 (productionization)** — 9 of 11 plans complete, wave 6 of 8 DONE; **07-08 is
COMPLETE**, next plan is **07-09** (wave 7 — `autonomous: false`, long 32×32 training run, MUST NOT
start without explicit human authorization).

**07-08 (2026-07-21):** the windowed sub-tile local colocalization map shipped
(`src/amortized/local_map.jl`, commits 25037e6 / 90f612e). `local_coloc_map(img, control,
channels; grid=8, tiles=(r,c))` cuts both images into r×c sub-tiles with the package's own
`patch()` tiler and scores each tile against its control tile through the FROZEN, already-gated
bundle — no training, no new gate, so the map inherits the 8×8 gate verdict verbatim (residual
FAILs included). Because 07-10 has not wired Artifacts yet, `_ensure_grid_registered(grid;
artifact_dir)` bridges in-process from the local `artifacts/grid_G/*.jld2` via
`_bundle_from_artifacts` + `register!` (touches neither `_lazy_load_from_artifact!` nor
Artifacts.toml). Degenerate tiles (below the grid size, zero survivors of the ≥15 floor, or a
non-finite read) take `LOCAL_MAP_SENTINEL = 0.0` AND a forced OOD flag — never a crash (T-7-04).
Phase-12 `SpatialColocResult`/`delta_rho_map`/`uncertainty_map` remain sketch-only in
`src/results.jl` (machine-asserted absent from the new module's executable code). `Pkg.test()`
green: windowed-map 30/30 plus a conditional real-8×8 smoke 5/5 (zero sentinels at 2×2 on 256²).

**Two carry-forward notes for 07-10:** (1) `_ensure_grid_registered` is a NO-OP when the grid is
already registered — an explicit registration always wins over a disk artifact, so superseding one
needs explicit invalidation, not a second `_ensure_*` call. (2) `_train_grid_pipeline` persists
ood_nulls with the `:density` channel ONLY (no `:noise`, no `:model`, **no `:thr`**), so
`ood_verdict` currently computes a score but has no threshold and returns `flag = false`; per-tile
OOD flagging is structural (degenerate tiles) only until the fitted `thr`/`zref` are persisted.
The 16×16 ship-gate RAN to completion (~16.5 min, `gate_report_16.jld2`, `status = :ran`) against the
byte-locked `gate_consts_16.jl` (7e2318b) on the fresh disjoint `PROD_SEED[16] = 0xb906f369f6cacf91`.
Recorded verdict is an **honest FAIL** (`gate-16x16.md`, commits 7f093a1 / e8f6a9d): SBC
`ks_pass=false` / `ece_pass=true` with **only 4/8 KS pass and the headline ρ_true itself rejected at
p = 7.0e-4**; BF corr 0.9150 < 0.95 AND max|Δ logBF| 12.72 > 0.5 (n=15 finite of 25); OOD
**NOT RUN / INCONCLUSIVE** (`auc = nothing`, `passed = nothing` — no `pos_sim` injected; ID operating
point 378.44 only, explicitly NOT a pass). No pre-registered constant was tuned.
**07-10 consequence:** grid 16 is recommended **NOT eligible** for default registry population — the
ρ_true rejection has no documented non-method cause, unlike the 8×8/4×4 residuals.

**Phase 08 (external-physical-ground-truth-corpus) — COMPLETE (2026-07-21).** 5 of 5 plans done.
Wave 5 (08-05) is unblocked and executed: both physical anchors human-verified, pinned, and sealed.
Offline gate `julia --project=. corpus/test/runtests.jl` → 217/217, ZERO network; `src/` provably
untouched; `git ls-files corpus/data` empty. NOTE: this Current Position was previously clobbered by a
parallel Phase-8 run — Phase 7, not Phase 8, is the active phase.

Last activity: 2026-07-21 -- Phase 07 plan 07-08 executed (windowed sub-tile local map, suite green); Phase 8 remains complete

Progress: [█████░░░░░] 56% of phases (9/16); 51/53 plans

## Resolved (2026-07-03): 07-00 CO-RESOLUTION GATE — GREEN

The Phase-7 Wave-0 co-resolution HARD GATE (07-00 Task 3) initially failed on a SECONDARY
GLMakie/Makie conflict (NOT the Finding-1 NeuralEstimators downgrade, which the Turing→weakdep
surgery had already cleared). Root cause: `GLMakie = "0.10.5"` forced Makie 0.21.18, while
NeuralEstimators 0.2.1 requires Makie 0.24.x.

**User decision (Option 2):** bump GLMakie compat (`0.10.5 → 0.13`, Makie 0.24.x) and keep
GLMakie as a normal CORE dependency (NOT an extension). Applied; the root `Pkg.resolve()` is now
GREEN and `Pkg.test()` passes fully.

Resolved versions (regenerated root Manifest.toml): NeuralEstimators **0.2.1** (not downgraded),
Flux **0.16.10**, Makie **0.24.12**, GLMakie **0.13.12**. Turing/CUDA remain weakdeps; the
`ext/ProteinCoLocTuringExt.jl` extension is unchanged. No `src/plot.jl` API changes were needed
at precompile/load (Makie 0.24 runtime plotting correctness deferred with the plotting tests).
spike/Project.toml + spike/Manifest.toml provably UNTOUCHED throughout.

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
| Phase 07 P07-02 | 45min | 3 tasks | 8 files |
| Phase 07 P07-03 | 50min | 4 tasks | 9 files |
| Phase 07 P07-04 | 55min | 3 tasks | 7 files |
| Phase 7 P6 | 40min | 3 tasks | 4 files |
| Phase 07 P07-08 | 55min | 2 tasks | 4 files |

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
- [Phase 07-01]: Grid coupling centralized ONCE in src/amortized/summary.jl — summary_dim(G)=2G^2 / cont_rows(G)=G^2 / ratio_input_dim(G)=5G^2 are the single source of truth (8x8 reproduces 128/64/320); patch_summary(mci,G) + encode_d01 + _summary_row_partition(:min,2G^2) grid-general for G in {4,8,16,32}. Every downstream amortized module derives dims from these helpers (no per-grid re-hardcoding)
- [Phase 07-01]: Estimator registry (PROD-02) keyed by patch grid — _SHIPPED_GRIDS=(4,8,16,32), 64 DROPPED (D-04); estimator_for validates grid and throws a train_and_register-pointing ArgumentError (T-7-05, no silent default-grid fallback); has_cuda_device() weakdep-safe (D-06 graceful CPU fallback); _lazy_load_from_artifact!/_train_grid_pipeline are honest hook-stubs for later plans
- [Phase 07-01]: Grid-parametric datagen (src/amortized/datagen.jl) — summary buffers Matrix(summary_dim(G),N), generating_config.summary_min_dim=summary_dim(G) so grids auto-separate into distinct content-hash cache dirs; imsize_set exposed for per-grid image-size bias (>=15-survivor floor, T-7-03); :min-only (spike :aug superset dropped); content hash uses Base.hash to avoid re-resolving the fragile Wave-0 Manifest; simulator chain referenced-but-promoted-later; spike/ byte-untouched
- [Phase 07-02]: Amortized READ surfaces promoted to src/amortized/{infer,bf,ood}.jl (PROD-01), grid-general + CPU-default (use_gpu=false everywhere, D-06). infer.jl: standardize_summary (frozen zt, mask bypass) + posterior_for/rho_draws/delta_rho (StatsBase.reconstruct BEFORE ρ read, Pitfall 5). bf.jl: amortized_log_bf (one NRE pass, measured log_prior_odds subtracted) + pair_encode (nc=G² derived ⇒ 5G²) + kde_log_bf_unclamped (Memo §5/T-7-06: compute_BayesFactor KDE math WITHOUT the 1e-8 _clampp floor so max|Δ logBF| is artifact-free). ood.jl: density (fit_ood_nulls continuous-rows+ridge) + noise (10 invariant features) + re-enabled posterior-predictive channel (ood_verdict defaults with_pp=true, Memo §5/T-7-04, kept crash-free by the iter1 _finite_or/_theta_tuple finite-guard) OR-fused by ood_verdict → OODVerdict. OOD ship-gate experiment (misspec families/ood_roc_over_grid, need simulator+ImageFiltering) deferred to 07-03+.
- [Phase 07-02]: Declared LinearAlgebra as a direct stdlib dep (cholesky/Symmetric/I for the OOD Mahalanobis) — a stdlib already in the Manifest with no version to resolve, so the co-resolution gate stays 4/4 green (NeuralEstimators 0.2.1 / Flux 0.16.10 pins intact); Pkg.test fully green (infer 11/11, bf 18/18, ood 31/31); spike/ byte-untouched.
- [Phase 07-04]: Per-grid CPU-reproducible ship-gate machinery (D-05) delivered to test/gate/{harness,sbc,gate_consts_template,run_gate}.jl + test/gpu_smoke.jl (PROD-02) — no grid trained yet; this is the gate the per-grid plans invoke. harness.jl draw_simulate_infer(m,rng;G,...) + paired-Δρ path is grid-parametrized (patch_summary(mci,G)), CPU-only (use_gpu=false on EVERY NeuralEstimators call), and loads a per-grid net via load_estimator (not the fixed spike trained_npe.jld2); frozen zt/θzt applied never re-fit (Pitfall 5). sbc.jl: M×8 rank table (7 θ + dedicated paired-draw Δρ column), KS+χ² uniformity via HypothesisTests (never hand-rolled), coverage curve, ported ECE/MCE CalibrationResult traffic-light, aggregated sbc_gate verdict. gate_consts_template.jl is the FRESH-per-grid pre-registration (SBC/BF/OOD consts) carrying a disjoint PROD_SEED[G] via a Philox stream salted (PROD_SALT) off the FORBIDDEN VAL_MASTER_SEED=0x5BC0FFEE and NPE_MASTER_SEED=0xC0FFEE — unit-asserted PROD_SEED[G]∉{those} for G∈{4,8,16,32} + 4 distinct seeds (anti-snooping Pitfall 3/T-7-08); BF gate uses the non-clamped kde_log_bf_unclamped baseline (Memo §5/T-7-06), OOD gate with_pp=true. run_gate.jl per-grid CLI --grid G [--sbc --bf --ood] loads the grid net + selected artifacts, runs the gates (use_gpu=false), writes an atomic .tmp→integrity→mv gate report; invocable before any grid is trained (missing NPE → :not_trained, no crash). bf_gate is self-contained (amortized_log_bf vs non-clamped KDE, no Turing); ood_gate computes the pre-registered id_threshold now + control-separability AUC when the per-grid plan injects pos_sim. gpu_smoke.jl: train_npe(...;use_gpu=true) runs on GPU when present and DEGRADES CLEANLY to CPU when CUDA absent (no error, CLAUDE.md graceful-fallback), asserts CPU-resident Flux.state persistence reloads/infers CPU-side; CUDA-present branch guarded by has_cuda_device() (==false here → fallback path exercised). Forward simulator NOT yet promoted (referenced-only in datagen/ood bodies) → harness takes an injectable `sim` seam (default_simulator resolves the promoted chain, errors clearly until then), mirroring 07-03's injectable datagen; the fixture SBC smoke injects a lightweight fake simulator. Pkg.test green (gate harness+SBC 29/29, GPU smoke 8/8, co-resolution 4/4 intact, no new external deps); spike/ byte-untouched.
- [Phase 07-07]: 16x16 (SHIP-WITH-CAVEAT fine grid) trained via _train_grid_pipeline(16) (80k pairs, 68k train, NPE d_in=512/dstar=64/10 coupling, NRE input_dim=1280) and CPU-gated on fresh disjoint PROD_SEED[16]=0xb906f369f6cacf91 with SBC_IMSIZE raised to 512^2. HONEST FAIL, constants byte-locked: SBC ks_pass=false/ece_pass=true with only 4/8 KS pass and the HEADLINE rho_true itself rejected (KS p=7.0e-4, ECE 0.0240) plus spillover 0.0027 -- both were comfortably uniform at 8x8 (0.567) and 4x4 (0.427), so the documented 'M=2000 over-sensitivity on summary-uninformative nuisance parameters' reading does NOT cover this; the paired Delta-rho column survives (KS 0.161, chi2 0.897, ECE 0.0060). BF corr 0.9150 (LOWEST of the three grids: 8x8 0.9472, 4x4 0.9410) AND max|d logBF| 12.72 (documented KDE-tail artifact), n=15 finite of 25. OOD NOT RUN / INCONCLUSIVE -- auc=nothing, passed=nothing (no pos_sim positive control injected); id_threshold=378.44 only, explicitly NOT a pass. CONSEQUENCE for 07-10: grid 16 recommended NOT eligible for default registry population (leave out of _SHIPPED_GRIDS or ship only behind explicit sign-off), and any entry must carry the >=512^2 minimum-image-size caveat (256 px/patch at 256^2 is marginal after background exclusion). PROVENANCE GAP: _train_grid_pipeline does not persist imsize_set, so the realized training image distribution is not recoverable from artifacts/grid_16/* alone.
- [Phase 07-08]: Local localisation ships as WINDOWED sub-tile inference over the frozen 8x8 bundle (src/amortized/local_map.jl), NOT a new fine-grid model — local_coloc_map returns a lightweight LocalColocMap (grid, tiles, delta_rho::Matrix, ood_flag::Matrix, meta), deliberately NOT an AbstractColocResult (no draws / no BF / no per-region uncertainty to back that interface). grid is a defaulted KEYWORD (8) so the map inherits whichever grid's gate the caller picks and CI can drive a tiny grid-4 bundle. Sub-tiles come from the package's OWN patch(img,nx,ny) (same trimming as the summary path) and tile MCIs carry the PARENT Otsu thresholds (tile-local Otsu would make tiles incomparable). Degenerate tiles take LOCAL_MAP_SENTINEL=0.0 AND a forced ood_flag=true so 'unscorable' can never read as 'measured no difference' (T-7-04). Wave-6 bridge _ensure_grid_registered(grid; artifact_dir) reuses _bundle_from_artifacts + register! rather than restating the three loads; it is a NO-OP on an already-registered grid.
- [Phase 07-03]: Amortized TRAINING layer promoted to src/amortized/{architecture,train_npe,train_ratio,persist,pipeline}.jl (PROD-01/02). build_estimator input-width-agnostic (q a NormalisingFlow INSTANCE positional; NPE_* Phase-5 consts); train_npe/train_ratio default use_gpu=has_cuda_device() with the spike use_gpu&&throw guards REMOVED (D-06) and LR/decay kept Float64 (AdamW-CosAnneal gotcha); ratio conditioner width from ratio_input_dim(G)=5G² (not 320), no custom loss (v0.2.1 hard-codes logit-BCE). SHIPPED PERSISTENCE migrated to CPU-resident Flux.state(cpu(est))+arch metadata → build_estimator/build_ratio_estimator + Flux.loadmodel! on load (Pitfall 4/T-7-07: device-independent, narrower deserialization surface T-7-01); OOD nulls persist directly; all through the atomic .tmp→reopen-integrity-@assert→mv(force=true) wrapper + schema_version + _estimator_ok/_ratio_ok/_ood_nulls_ok SKIP-IF-DONE predicates; seeded save→load→CPU inference bitwise-equal. _train_grid_pipeline(grid;...) factors datagen→zt→train_npe→train_ratio→fit_ood_nulls→persist→EstimatorBundle (SKIP-IF-DONE loads valid artifacts; injectable datagen seam makes it testable without the not-yet-promoted simulator; default_imsize_for 16→≥512²/32→≥1024²); train_and_register(grid) now runs it (D-04 on-ramp is real code, stub removed). Docstring documents -t auto datagen + Philox-per-index thread-count-independent byte-identical repro. GPU path plumbed but NOT exercised (has_cuda_device()==false, no CUDA loaded); CPU path fully green. Pkg.test green under -t auto (train_npe 11/11, train_ratio 12/12, persist 17/17, pipeline 19/19; co-resolution 4/4, no new external deps); spike/ byte-untouched.
- [Phase 07-06]: 4x4 bundle produced via the PUBLIC train_and_register(4) path (PROD-02/D-04 user-definable-grid happy path proven end-to-end; estimator_for(4) confirmed in-process), imsize_set constrained to ((256,256),) for the CPU budget; gate on fresh disjoint PROD_SEED[4]=0x8c0ad97b99bd6031 recorded honest FAIL (SBC ECE-green all 8, 7/8 KS pass vs 8x8 5/8, only shift_dx 0.0326; BF corr 0.941 near-miss + max|dlogBF| 12.40 KDE-tail artifact; OOD ID-op 31.79 only), gate_consts_4.jl byte-locked
- [Phase 07-06]: fixed a Rule-1 numerical crash in the shipped non-clamped KDE BF baseline (bf.jl kde_log_bf_unclamped) — QuadGK adaptive integral overshoots the tail probability a few ulp outside [0,1], crashing log(); clamped to the valid [0,1] domain so a saturated tail yields the honest ±Inf (dropped) NOT the forbidden 1e-8 finite floor; in-range values byte-identical
- [Phase 08-05]: BOTH physical ground-truth anchors human-verified and pinned (SC1). POSITIVE = TetraSpeck 100 nm multicolor fiducial beads (RegiSTORM v1.0.0 sample data, Zenodo DOI 10.5281/zenodo.5509861, CC-BY-4.0, Karlsson et al. 2023 BMC Bioinformatics 10.1186/s12859-023-05320-1): one physical bead emits in EVERY colour channel ⇒ coloc BY CONSTRUCTION and state-independent (multi-channel STORM .tif frame stacks + ThunderSTORM .csv; a mean/max projection is owed in Phase 16). NEGATIVE = "Light My Cells" (France-BioImaging / ISBI 2024, BioImage Archive S-BIAD1047, CC-BY-4.0, OME-TIFF): nucleus (DNA) vs mitochondria are disjoint BY BIOLOGY from the acquisition design (fields must be FILTERED in Phase 16 to those carrying BOTH channels). Neither label is a computed coloc score (Pitfall 1) and neither is an environment-quenched tandem (Pitfall 2).
- [Phase 08-05]: D-01 DEVIATION, human-accepted and recorded (never papered over) — the literal decision specified a CELLULAR tandem-fluorophore construct, but no open-licensed, NON-environment-quenched tandem-FP dataset could be verified in any public archive; the deposited tandem-FP sets are quenched mRFP/mCherry-EGFP-LC3 autophagy reporters, which the acceptance bar disqualifies. Substituted a multicolor-bead dataset: the SAME physical particle in both channels is a strictly STRONGER physical positive (no pH/quenching failure mode at all), but it IS a deviation from the literal wording.
- [Phase 08-05]: D-02 DEVIATION — the matched-same-study preference is NOT satisfiable; no qualifying single-study deposit carrying both a same-particle positive and a segregated negative could be verified. The anchors are cross-study / cross-archive (`ANCHORS_MATCHED = false`). KNOWN LIMITATION: imaging-condition confounds (microscope, objective, exposure, detector, prep) between positive and negative are NOT controlled; Phase 16 must report this.
- [Phase 08-05]: Hash discipline (T-08-16) — `sha256` uses a deliberately NON-hex `"PENDING-FETCH"` sentinel rather than an empty string (an empty hash is indistinguishable from an unfilled CBS row and could be silently accepted downstream), with `is_real_sha256`/`is_pending_hash` making "real 64-hex digest OR sentinel, never in between" a testable invariant. The bootstrap fetch was DELIBERATELY not executed — the positive anchor is a ~6.3 GB archive, so the bulk download stays an explicit human decision. `bootstrap_anchor_hashes()` is implemented and must be run online+authorized BEFORE Phase 16. No digest was fabricated.
- [Phase 08-05]: `finalized_manifest()` (corpus/anchor_rows.jl) SUPERSEDES the 08-04 placeholder `committed_manifest()`; `anchor_rows()` is parameterized on sha256/bytes so the pending path and the post-bootstrap path are the SAME code path and both are covered by the offline gate. All anchors remain `physical-primary` + `sealed_holdout`, reachable ONLY via `open_sealed_holdout(; reason)` (D-09) — asserted.

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
- [Phase 7 — F5 COVARIATE SHIFT (HIGHEST SEVERITY), 2026-07-21]: **The SBC calibration proof holds at 256² and is NOT demonstrated to transfer to the real 1376×1028 paper regime that Phase-2 decision D-08 anchors the design on.** Training runs 07-05/07-06 constrained `imsize_set` to `((256,256),)` as a documented compute-budget deviation, and the gates evaluate at `SBC_IMSIZE` = (256,256) for grids 4/8 and (512,512) for grid 16 — train and gate are self-consistent, so the recorded verdicts are internally valid, but they are valid *at 256²*. MEASURED on the grid-8 net: evaluating on the mixed imsize set instead of 256² moves `label_efficiency` mean(u) from 0.546 → **0.807** (z ≈ 16.8 at M=250), bias −0.130. At FIXED θ, changing image size shifts the mean per-patch correlation by roughly the ENTIRE `label_efficiency` prior range (G=8: 0.513 @2048² vs 0.470 @256² at le=0.95, Δ≈0.043, while le 0.65→0.95 moves it only 0.045) — the confound is as large as the signal. **Consequence: the chain "SBC passes ⇒ posteriors on real images are calibrated" does not currently hold, and the Phase-8 physical anchors and the Phase-16 blind evaluation both depend on exactly that chain.** Closing it requires training/gating at the real image dimensions, or demonstrating and gating summary invariance to image size — both need retraining and are OUT of the current scope. Full evidence: `.planning/phases/07-productionization-conditional-on-go/07-CALIBRATION-FINDINGS.md` (F5). Related F6: `_train_grid_pipeline` never persisted `imsize_set`, so grid 16's actual training image-size distribution is UNRECORDED; a zt forensic mildly favours 512²-only (0.4249 vs 0.4288 for a 512²-only pool, 0.405 for the mixture) but is explicitly NOT CONCLUSIVE at n=100 and the value is recorded as `unknown`, never guessed. Provenance is now persisted going forward (`training_imsize_provenance`); existing bundles honestly report `recorded=false`.
- [Phase 7 — F3 VACUOUS SBC PASS, 2026-07-21]: **grid 4's `label_efficiency` SBC pass is VACUOUS and must NOT be read as calibration evidence anywhere** (memo, manuscript, or any per-grid comparison table). MEASURED: post_sd/prior_sd = **0.994**, shrinkage centre = **0.8034** = the prior mean (0.80), tercile mean(u) = 0.18/0.49/0.82 — the textbook "posterior = prior" signature. Grid 4 passes SBC on this parameter (KS p = 0.567) **because it learns nothing about it**: a posterior that reproduces the prior yields uniform ranks BY CONSTRUCTION. Rank uniformity is NECESSARY BUT NOT SUFFICIENT for calibration, and a vacuous pass and an informative pass must never be tabulated as the same result. Related F1: where the SBC *does* fail on `label_efficiency` (grid 8 z=+7.4 p=1.84e-9; grid 16 z=+8.6 p=7.15e-13) the failure is a LOCATION BIAS, not overconfidence — tail mass 0.199/0.231/0.215 vs expected 0.20, histograms flat not U-shaped — so the standard overconfidence levers (widening, tempering, more capacity) are the WRONG remedy; the NPE systematically UNDERESTIMATES the parameter. F2 mechanism (INFERRED): the unbounded `NormalisingFlow` cannot represent the `Uniform(0.6,1.0)` box and the forward-KL objective constrains nothing about the implied marginal, so information gained, centre displacement and rank bias all grow monotonically 4→8→16 — **the bias grows with learning, not with ignorance**. MITIGATION SHIPPED (reporting only, no threshold or verdict changed): `sbc_gate` now reports per-parameter `shrinkage = post_sd/prior_sd` + a `vacuous` flag (cutoff `SBC_VACUOUS_SHRINKAGE = 0.95` in `test/gate/sbc.jl`, deliberately NOT in the frozen `gate_consts_<G>.jl`).
- [Phase 8-05 — RESOLVED 2026-07-21]: The human-verify blocker is CLEARED. The human confirmed both anchors and 08-05 Task 3 executed (commit 8eb9147): `corpus/anchor_rows.jl` + `corpus/test/test_anchors.jl` created, `corpus/manifest.csv` finalized (30 CBS + 2 sealed physical anchors), offline gate 217/217, `git ls-files corpus/data` empty, `src/` untouched. **Phase 8 is COMPLETE.** THREE residual items carried to Phase 16, all recorded honestly in `ANCHOR_PROVENANCE_NOTES`: (1) **sha256 = "PENDING-FETCH" on BOTH anchors** — the bootstrap fetch was DELIBERATELY not run because the positive anchor is a ~6.3 GB archive and a bulk download must stay an explicit human decision; no digest was fabricated (T-08-16). `bootstrap_anchor_hashes()` is implemented and must be run when online + authorized, and the sentinel replaced BEFORE Phase 16 opens the sealed holdout. (2) **D-01 substitution (human-accepted)** — no open-licensed, non-environment-quenched tandem-fluorophore dataset exists in any public archive (every deposited tandem-FP set is a quenched mRFP/mCherry-EGFP-LC3 autophagy reporter, disqualified by Pitfall 2). POSITIVE = TetraSpeck 100 nm multicolor beads (RegiSTORM sample data, Zenodo 10.5281/zenodo.5509861, CC-BY-4.0): the same physical particle emits in both channels, so coloc is by construction AND state-independent. (3) **D-02 preference unmet** — NEGATIVE = "Light My Cells" (BioImage Archive S-BIAD1047, CC-BY-4.0, nucleus vs mitochondria). Cross-study / cross-archive (`ANCHORS_MATCHED=false`); imaging-condition confounds between the two anchors are NOT controlled and Phase 16 MUST report this limitation. Phase-16 conversion work still owed: mean/max frame projection for the positive STORM stacks, and field-filtering the negative to fields carrying BOTH the nucleus and mitochondria channels.

## Deferred Items

Items acknowledged and carried forward from previous milestone close:

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| Backend (v2) | BACK-01 Turing→RxInfer.jl feasibility evaluation | Deferred (post-Go, never a spike dependency) | 2026-06-26 |
| Summary net (v2) | BACK-02 DeepSet permutation-invariant summary upgrade | Deferred (upgrade if MLP summary proves insufficient) | 2026-06-26 |

## Session Continuity

Last session: 2026-07-21T08:00:00.000Z
Stopped at: Phase 7 plan 07-07 COMPLETE (16x16 ship-gate ran; honest FAIL + unscored OOD recorded in gate-16x16.md)
Resume file: .planning/phases/07-productionization-conditional-on-go/07-08-PLAN.md
Resume action: continue the ACTIVE phase — Phase 7, plan 07-08, wave 6 of 8. Phase 8 remains COMPLETE.
NOTE: plan 07-09 (wave 7, the 32x32 grid) is `autonomous: false` and carries a long training run —
it MUST NOT be started without explicit human authorization. Separately, before Phase 16: run `bootstrap_anchor_hashes()` (corpus/anchor_rows.jl)
when online and authorized to replace the `PENDING-FETCH` sentinel on both physical anchors with real SHA-256
digests (the positive anchor is a ~6.3 GB download — an explicit human decision).
