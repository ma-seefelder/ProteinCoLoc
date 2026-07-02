# Phase 5: Validation Bundle (SBC + Amortized BF + OOD) - Research

**Researched:** 2026-07-02
**Domain:** Simulation-Based Calibration (SBC), amortized Bayes-factor / neural ratio estimation, OOD/misspecification detection — over a pinned Julia NeuralEstimators.jl v0.2.1 spike
**Confidence:** HIGH (API verified against installed v0.2.1 source; existing spike code read directly)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01 (SBC scope):** SBC covers all 7 per-stack θ (rank histograms from `sampleposterior` marginals) AND a **dedicated paired-draw Δρ SBC** — draw paired θ\* (sample θ\*_s, control θ\*_c), simulate both, rank the true Δρ\* = ρ\*_s − ρ\*_c among **differenced posterior draws** (MC difference of two independent single-stack passes, per Phase-4 D-03).
- **D-02 (anti-snooping):** M and pass/fail thresholds **pre-registered**; reported from a **fresh, independently-seeded, never-tuned-against held-out run** (Phase-4 pre-registered-consts + keyed holdout-repro-gate pattern). M ≈ 2000 (exact M + L = research). SBC reported **"calibrated under the simulator"** and explicitly paired with the OOD result (SBC-04).
- **D-03 (misspec positive controls):** ALL FOUR families in the ROC grid: (1) texture-model mismatch (puncta/granular/fibrillar), (2) noise-model mismatch (non-Gaussian/Poisson-shot/salt-pepper/structured), (3) optics/PSF mismatch (asymmetric/out-of-focus/aberrated/chromatic-warp), (4) background/illumination (gradients/vignetting/autofluorescence bleed/debris). Magnitudes + grid resolution = research.
- **D-04 (summary-orthogonal negative control):** correlation-preserving transforms — per-channel **affine intensity rescale (a·x+b)**, global **rotation/flip**, distribution-preserving spatial rearrangement. **Verify empirically** (KS test on 8×8 summary) that the summary is statistically unchanged so the OOD flag stays quiet by construction. KS-invariance ε = research.
- **D-05 (OOD fusion):** Report **both channels separately AND an OR-combined flag.** Summary-density channel (Mahalanobis OR normalizing-flow log-lik on the summary) + posterior-predictive mismatch channel; each its own ROC; combined flag fires if EITHER exceeds threshold. Density-vs-flow choice + PP discrepancy statistic = research.
- **D-06 (OOD threshold discipline):** Threshold-free **ROC/AUC headline** + **pre-registered ID-quantile operating point** (e.g. ~5% ID FPR) with density/flow + PP nulls **fit on TRAIN split only** and quantile committed **before** seeing misspecified data + **Youden-J as post-hoc reference only** (never gating). ROC on **fully external misspecified test images** (SC5).
- **D-07 (BF null):** **Read `compute_BayesFactor()` first, then mirror its null + prior EXACTLY.** (Body read — see §"BF Null, Confirmed".)
- **D-08 (BF agreement):** Pre-register **BOTH** (a) high rank/Pearson correlation between amortized and KDE log-BF across a held-out Δρ sweep AND (b) bounded absolute error in log-BF units within a pre-set tolerance over the decision-relevant Δρ range. Computed **without quadgk/KDE/shuffle** on the amortized side.

### Claude's Discretion / Research
Exact M/L; KS/χ² p-value thresholds; ECE/MCE traffic-light cutoffs; Mahalanobis-vs-flow choice + flow architecture; PP discrepancy statistic; perturbation magnitudes/grid resolution per family; correlation-preserving transform parameterization + KS ε; exact ID FPR quantile; BF correlation threshold + log-BF error tolerance + Δρ sweep range; shared-harness module layout; figure set. **All numeric thresholds pre-registered before the reported run.**

### Deferred Ideas (OUT OF SCOPE)
Demo + Go/No-Go memo (Phase 6); productionization / `src/` edits / user-definable `num_patches` (Phase 7); external physical corpus (Phase 8); cross-method comparator (Phase 9); three-hypothesis BF (Phase 13); DeepSet summary upgrade + adversarial nuisance/CI gate (Phase 15). This phase is "under the simulator," two-way coloc-vs-null only.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| SBC-01 | M(≈2000) θ\*~π → simulate → L posterior draws → rank of θ\* → per-parameter rank histogram | §SBC Mechanics: `sample_prior`→`simulate_pair`→`patch_summary`→`encode_d01`→`standardize_summary(zt)`→`posterior_for`(=`sampleposterior`,N=L)→`reconstruct(θzt)`→per-row rank. Δρ paired path (D-01) via `rho_draws` difference. |
| SBC-02 | Uniformity (KS/χ²) + coverage curve + ECE/MCE via `_bin_calibration` + traffic-light | §SBC Uniformity/Coverage: `ExactOneSampleKSTest`/`ChisqTest` (HypothesisTests, pinned); coverage curve hand-rolled or via `assess`+`coverage`; port `_bin_calibration`/`CalibrationResult` from BayesInteractomics (template only). |
| SBC-03 | M + pass/fail threshold pre-registered; fresh never-tuned held-out run | §Pre-Registration Discipline: mirror Phase-4 `test_npe.jl` const block; fresh Random123 disjoint counter stream (D-02). |
| SBC-04 | Report "under the simulator", paired with OOD | Reporting convention only (doc/figure caption + memo hand-off to Phase 6). |
| BF-01 | RatioEstimator/Evidence-Network amortized log-BF (coloc vs null) in one forward pass, Δρ + prior aligned to KDE baseline | §Amortized Bayes Factor: **l-POP loss NOT in v0.2.1** — use NRE-as-model-comparison via `RatioEstimator` + `logratio`, model-index parameter defined by sign(Δρ). Prior alignment confirmed. |
| BF-02 | Reproduces `compute_BayesFactor()` in well-specified regime, no quadgk/KDE/shuffle | §BF Reproduction: amortized log-BF = `logratio(Z;grid=1) − logratio(Z;grid=0)`; baseline BF read (odds-ratio of Δρ>0). Correlation + bounded-error gate (D-08). |
| OOD-01 | OOD score (summary-density + PP mismatch) fires on misspecified, quiet in-distribution | §OOD Channels: Mahalanobis on continuous summary rows (fit on train) + PP-mismatch; OR-combined (D-05). |
| OOD-02 | Controlled ROC over misspec grid with positive AND summary-orthogonal negative controls; blind spot measured/named | §Misspecification Grid: 4 positive families (D-03) + correlation-preserving negative control (D-04) with KS-invariance verification. |
</phase_requirements>

## Summary

This phase builds three sibling deliverables on **one shared `θ*~π → simulate → infer` harness** over the already-trained Phase-4 net (`spike/npe/trained_npe.jld2`), reusing the existing amortized read surface in `spike/npe/infer.jl` (`posterior_for`, `rho_draws`, `delta_rho`). All three consume the **frozen** train-split standardization (`m.zt` for the summary, `m.θzt` for θ) persisted inside `trained_npe.jld2` — nothing is re-fit on evaluation data (SC5/D-06).

**The single highest-impact finding:** the pinned **NeuralEstimators.jl v0.2.1 does NOT implement an l-POP / Evidence-Networks loss.** Its `RatioEstimator` is a standard Hermans-2020 neural ratio estimator trained with `logitbinarycrossentropy` (the loss is hard-coded: `_loss(::RatioEstimator, loss=nothing) = logitbinarycrossentropy` — a custom loss passed to `train` is *ignored* for this type). The BF requirement's parenthetical "l-POP-style loss" therefore cannot be met literally without hand-rolling a new estimator. The good news: the amortized coloc-vs-null Bayes factor is achievable **exactly** via the standard `RatioEstimator` used as a **neural model-comparison / Evidence Network** — encode a binary model index m∈{0,1} (coloc = Δρ>0, null = Δρ≤0) as the estimator's 1-D "parameter", and read `log-BF = logratio(Z; grid=[1]) − logratio(Z; grid=[0])` in one forward pass, with no quadgk/KDE/shuffle. This provably equals the baseline's posterior-odds/prior-odds quantity (proof in §Amortized Bayes Factor).

SBC is hand-rolled ranks (CLAUDE.md convention) tested with `HypothesisTests` (both `ExactOneSampleKSTest` and `ChisqTest` are available in the pinned env). The OOD channels use only pinned/stdlib tools (Mahalanobis via `LinearAlgebra`; ROC/AUC hand-rolled). **No new packages are needed** — which is also the safest path, because adding any dependency risks the Phase-1 co-resolve landmine that silently downgrades NeuralEstimators below the pinned 0.2.1.

**Primary recommendation:** Lay the harness in `spike/validation/` (sibling to `spike/npe/`) as a shared `harness.jl` (draw→simulate→infer, frozen-stats loader) consumed by `sbc.jl`, `bf.jl` (with a `train_ratio.jl` sibling to `train_npe.jl`), and `ood.jl`; land re-runnable gates + pre-registered consts in `spike/test/test_validation.jl` wired into the existing `spike/test/runtests.jl`; use only the pinned spike env.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| θ\*~π draw | Simulator (`spike/simulator/prior.jl`) | Harness | `sample_prior(rng)` is the sole prior; harness threads Random123 counter |
| Simulate image pair | Simulator (`spike/simulator/forward.jl`) | Harness | `simulate_pair(rng,θ;imsize)` unchanged; misspec grid perturbs *around* it |
| Summary extraction | Contract (`spike/contract.jl` → `src/` read-only) | Encode | `patch_summary` (8×8 pearson, UNCHANGED src) → `encode_d01` 128-dim |
| Standardization (freeze) | Loader stats via trained model (`m.zt`,`m.θzt`) | infer.jl | Frozen train-split transforms; never re-fit on eval (SC5) |
| Posterior sampling | NPE (`spike/npe/infer.jl` → NeuralEstimators) | — | `sampleposterior(PosteriorEstimator)` in one pass; reuse `posterior_for`/`rho_draws` |
| SBC ranks + uniformity | Validation (`spike/validation/sbc.jl`) | HypothesisTests | Hand-rolled ranks (CLAUDE.md), KS/χ² tests |
| Amortized BF | Validation (`spike/validation/bf.jl` + `train_ratio.jl`) | NeuralEstimators `RatioEstimator` | NRE model-comparison; `logratio` at m∈{0,1} |
| BF baseline (reproduction target) | `src/bayes.jl` (READ-ONLY via `include`) | baseline env | `compute_BayesFactor` KDE odds-ratio — never edited |
| OOD scoring | Validation (`spike/validation/ood.jl`) | LinearAlgebra/StatsBase | Mahalanobis + PP mismatch, ROC/AUC |
| Figures | Validation (`spike/validation/figures.jl`) | CairoMakie | rank hist, coverage, ECE, ROC, traffic-light |

## Standard Stack

**No new packages.** Every tool below is already in `spike/Project.toml` (pinned) or is a Julia stdlib. This is a deliberate constraint: the Phase-1 co-resolve landmine means any *added* dependency can silently downgrade NeuralEstimators from the pinned 0.2.1 to 0.1.4 and break the v0.2.1 API surface.

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| NeuralEstimators.jl | 0.2.1 (pinned) | `sampleposterior` (SBC draws), `RatioEstimator`/`logratio`/`train` (BF), `assess`/`coverage` | The trained net + the only ratio API; verified against installed source `[VERIFIED: installed pkg gFxuZ]` |
| Flux.jl | 0.16.10 (pinned) | NN backend for the ratio summary/inference nets (`Chain`,`Dense`,`gelu`,`Optimisers.AdamW`) | Same backend the NPE trained on |
| HypothesisTests.jl | pinned | `ExactOneSampleKSTest`, `ChisqTest` for SBC rank uniformity (SBC-02) | CLAUDE.md-mandated; already co-resolves cleanly (STATE.md) |
| StatsBase.jl | pinned | `ZScoreTransform`/`reconstruct` (frozen stats), rank/histogram helpers, `wsample` (used internally by ratio sampler) | Already the standardization lingua franca |
| Distributions.jl | pinned | `Uniform`, `Chisq`/`Cauchy` for prior draws + χ² reference; `Normal` for KS reference | Prior π(θ) source |
| Random123.jl | pinned | Counter-based disjoint-stream seeding for the fresh SBC run (D-02) | Phase-3 D-11 reproducibility discipline |
| JLD2.jl | pinned | Persist trained ratio net + pre-registered consts + frozen OOD nulls, atomic `.tmp`+`mv` | The cache.jl / save_npe idiom |
| CairoMakie.jl | pinned | Rank histograms, coverage curves, reliability/ECE, ROC, traffic-light (headless spike figures) | Phase-2 figure pattern |

### Supporting (stdlib — no dependency cost)
| Library | Purpose | When to Use |
|---------|---------|-------------|
| LinearAlgebra | Mahalanobis distance (`cholesky`/`\`, `Symmetric`), covariance inverse | OOD summary-density channel |
| Statistics | `mean`, `median`, `quantile`, `cor` | ranks, coverage, BF correlation gate |
| Test | Re-runnable gates | `spike/test/test_validation.jl` |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `RatioEstimator` (logit-BCE) for BF | Hand-rolled Evidence Network with l-POP loss (Jeffrey & Wandelt 2024) | l-POP is numerically nicer for extreme BFs but is **not in v0.2.1**; needs a bespoke `PointEstimator`+custom-loss `train` path (the `RatioEstimator` loss is hard-coded). Adds real scope for marginal benefit in the well-specified mid-range regime BF-02 targets. `[ASSUMED]` — confirm with user whether literal l-POP is required or the standard NRE model-comparison suffices. |
| Mahalanobis summary-density | NeuralEstimators `NormalisingFlow` fit on train summaries (log-lik score) | Flow catches non-Gaussian ID structure Mahalanobis misses, but needs a second trained density + tuning. D-05 says "pick by separability" — start Mahalanobis, escalate only if ROC AUC is weak. |
| Hand-rolled ROC/AUC | An ROC package (e.g. `ROCAnalysis`, `MLJ`) | Adding a package risks the co-resolve landmine; ROC/AUC from sorted scores is ~15 lines. Hand-roll. |
| Reuse Phase-4 NPE Δρ for BF | — | Would satisfy "no KDE" but violates BF-01's explicit `RatioEstimator` requirement. Train a ratio net. |

**Installation:** none. `Pkg.activate("spike")` uses the existing pinned Manifest. If ANY new package is ever proposed, it must first pass the `runtests.jl` resolve-risk gate (assert NeuralEstimators stays 0.2.1).

**Version verification:** `[VERIFIED: installed pkg]` NeuralEstimators 0.2.1 at `~/.julia/packages/NeuralEstimators/gFxuZ` (spike Manifest git-tree-sha1 `c780174…`); a stale 0.1.4 also present at `8G5tE` — the spike Manifest pins 0.2.1. Flux 0.16.10, HypothesisTests/CairoMakie/StatsBase/Distributions/Random123/JLD2 all present in `spike/Project.toml`.

## Package Legitimacy Audit

> This phase installs **no external packages**. All tools are pinned in the committed `spike/Manifest.toml` (already legitimacy-gated in Phases 1–4) or Julia stdlib.

| Package | Registry | Disposition |
|---------|----------|-------------|
| (none new) | — | Phase uses only the pinned spike env + stdlib |

**Packages removed due to slopcheck [SLOP] verdict:** none (none proposed).
**Packages flagged as suspicious [SUS]:** none.
**Guard:** the existing `runtests.jl` resolve-risk gate asserts NeuralEstimators==0.2.1 and no top-level CUDA — the planner should keep any BF/OOD code within these deps and NOT introduce a ROC/flow package.

## BF Null, Confirmed (D-07 — body of `compute_BayesFactor` read)

`src/bayes.jl:109` `compute_BayesFactor(posterior, prior; ρ_threshold=0.0)`:

```julia
Δρ_post  = posterior.posterior.μ_sample .- posterior.posterior.μ_control   # Δρ = μ_s − μ_c
p_post   = 1 - ∫_{-∞}^{ρ_threshold} kde(Δρ_post)          # = P(Δρ >  0 | data)
p_prior  = 1 - ∫_{-∞}^{ρ_threshold} kde(Δρ_prior)         # = P(Δρ >  0 | prior)
BF = (p_post/(1-p_post)) / (p_prior/(1-p_prior))          # posterior-odds / prior-odds
```

**Confirmed facts the amortized side must mirror exactly:**
- **Quantity:** Δρ = `μ_sample − μ_control` (the *induced per-condition mean* difference, a one-sided contrast). This is the difference the KDE is built on.
- **Hypotheses:** H1 (coloc) = {Δρ > 0}; H0 (null) = {Δρ ≤ 0}. It is a **one-sided interval hypothesis at threshold 0**, NOT a point null. `ρ_threshold=0.0` is the decision boundary.
- **BF form:** posterior-odds ÷ prior-odds for the event {Δρ>0}. `[VERIFIED: src/bayes.jl:109-136]`
- **Prior:** `μ_sample, μ_control` are each drawn i.i.d. from the Turing prior `Truncated(Cauchy(0,0.3),-1,1)` (`prior_chain` from `sample(m, Prior())`), which is symmetric about 0 ⇒ **prior odds of {Δρ>0} ≈ 1, log prior-odds ≈ 0** (verify numerically; the amortized side must subtract the *same* prior-odds it is trained under, not assume 0).

**Space alignment (ρ_true vs induced-μ):** the spike simulator sets `ρ_true = ghat(μ*)` with `ghat` **monotone increasing** (SIM-02), so `sign(ρ_s − ρ_c) == sign(μ_s − μ_c)`. The Phase-4 Δρ (D-03) is the MC difference of two single-stack ρ_true passes. Because the threshold is exactly 0 and ghat is monotone through 0, **{Δρ_ρ > 0} ≡ {Δμ > 0}** — the amortized model-index defined by sign(Δρ_ρ) reproduces the baseline's {Δμ>0} event. Document this equivalence in the BF plan; if `ghat(0) ≠ 0` numerically, define the coloc/null split on the induced-μ contrast to stay byte-faithful to the baseline (research: check `ghat(0)`).

## Amortized Bayes Factor — the l-POP gap and the correct v0.2.1 path

### What v0.2.1 actually provides
`RatioEstimator` estimates the **likelihood-to-evidence ratio** r(Z,θ)=p(Z|θ)/p(Z), trained on the log-scale via `logit(c*(Z,θ))` with a **hard-coded** `logitbinarycrossentropy` classifier loss (`RatioEstimator.jl:124`). Construction, training, and evaluation (verified verbatim, `[VERIFIED: RatioEstimator.jl gFxuZ]`):

```julia
# Source: RatioEstimator.jl:47,50,67 (installed v0.2.1)
summary_network = Chain(Dense(m, 64, gelu), Dense(64, 64, gelu), Dense(64, num_summaries))
estimator = RatioEstimator(summary_network, d; num_summaries = num_summaries)   # d = num_parameters
estimator = train(estimator, sampler, simulator; K = 1000)          # or fixed-data form (below)
logratio(estimator, z; grid = grid)          # K×G matrix of log r(z,θ) over grid columns
sampleposterior(estimator, z; grid = grid, logprior = θ -> 0f0)     # grid-based posterior
```

**No `pop`, `l-pop`, `evidence network`, or `bayesfactor` symbol exists anywhere in the package** (grepped). The `train` fixed-data form `train(est, θ_train, θ_val, Z_train, Z_val; ...)` works for `RatioEstimator` (the type's `_inputoutput` builds the contrastive dependent/independent pairs internally by shuffling θ). `[VERIFIED: train.jl:6, RatioEstimator.jl:94-122]`

### Recommended construction: NRE-as-model-comparison (Evidence Network)
Encode the **model index as the 1-D parameter**. Let m∈{0,1} with m=1 ⟺ Δρ>0 (coloc), m=0 ⟺ Δρ≤0 (null). Train a `RatioEstimator` with `num_parameters = 1` whose data `Z` is the **paired** summary (concatenate sample+control standardized summaries, 256-dim, OR a symmetric contrast encoding — research the encoding). Then:

```
log-BF(Z) = logratio(est, Z; grid = [1.0]) − logratio(est, Z; grid = [0.0])
```

**Why this equals the baseline BF (proof):** define models by restricting the prior to each region. Then p(Z|m=1)=∫p(Z|θ)π(θ|Δρ>0)dθ and p(Z|m=0)=∫p(Z|θ)π(θ|Δρ≤0)dθ. By Bayes, posterior-odds/prior-odds = [P(Δρ>0|Z)/P(Δρ≤0|Z)]/[P(Δρ>0)/P(Δρ≤0)] = p(Z|Δρ>0)/p(Z|Δρ≤0) = BF. The NRE ratio at m=1 vs m=0 estimates exactly p(Z|m=1)/p(Z|m=0). Hence `logratio(m=1)−logratio(m=0)` is the amortized log-BF, **in one forward pass, with no quadgk/KDE/shuffle** (the shuffle is a *training-time* NRE detail, not part of the amortized evaluation BF-02 forbids). `[CITED: Hermans et al. 2020 (linked in RatioEstimator docstring); ASSUMED for the model-index reduction — confirm framing with user]`

### Fallback (only if literal l-POP is mandated)
Hand-roll an Evidence Network: a plain binary classifier (`Chain(...)`→1 logit) trained with a custom l-POP-exponential loss; the logit at Z is the log posterior-odds; subtract log prior-odds. This bypasses `RatioEstimator` entirely (its loss is not swappable). Larger scope; reserve for a user "l-POP is required" decision. `[ASSUMED]`

### Training data for the ratio net
- Sibling to `train_npe.jl` → `spike/validation/train_ratio.jl` (CONTEXT.md: "the l-POP RatioEstimator likely lands as a sibling here").
- Draw paired (θ_s,θ_c)~π×π, simulate both via `simulate_pair`, encode both to 128-dim `:min`, **standardize with the frozen `m.zt`** (do NOT re-fit — SC5), concatenate, label by sign(Δρ). Use the **fixed-data** `train` form so the reserved-holdout and frozen-stats discipline hold (mirror `train_npe.jl` rationale). CPU-only (`use_gpu=false` on every call — Pitfall 1: `use_gpu` defaults to TRUE in NeuralEstimators).
- Reproduction gate (D-08): on a held-out Δρ sweep, compute amortized log-BF and the KDE baseline BF (from the isolated baseline env's ADVI posteriors OR from re-running the KDE odds-ratio on NPE Δρ draws), then assert (a) Pearson/Spearman correlation ≥ pre-registered threshold AND (b) `max|logBF_amortized − logBF_kde| ≤` pre-registered tolerance over the decision-relevant Δρ range.

## SBC Mechanics (SBC-01/02) — verified against `infer.jl` + NeuralEstimators

### Per-parameter marginal SBC (all 7 θ)
```
for i in 1:M:
    θ* = sample_prior(rng_sbc)                              # draw from π (prior.jl)
    mci = build_mci(simulate_pair(rng_sbc, θ*; imsize))     # simulate (forward.jl + contract.jl)
    S   = encode_d01(patch_summary(mci))                    # 128-dim raw summary
    Z   = standardize_summary(S, m.zt, :min)               # FROZEN train zt (infer.jl helper)
    draws_std = posterior_for(m.estimator, Z; N = L)        # 7×L, sampleposterior, one pass
    draws     = StatsBase.reconstruct(m.θzt, draws_std)     # 7×L un-standardized (physical θ)
    for p in 1:7: rank[p][i] = count(draws[p,:] .< θ*[p])   # rank ∈ 0:L
```
- `sampleposterior(::PosteriorEstimator, Z; N)` returns a d×N matrix for one dataset. `[VERIFIED: PosteriorEstimator.jl:130, inference.jl:56]`
- **Frozen stats:** `m.zt`/`m.θzt` come from `load_npe(trained_npe.jld2)` — the deployed model's train-split transforms (SC5). Never re-fit.
- **M ≈ 2000, L:** standard SBC needs L such that ranks span 0:L with a chosen bin count B where L+1 is a multiple of B (Talts et al. bin without remainder). Recommend **L = 999 (→ 1000 ranks 0..999), B = 20 or 50 bins**; M ≈ 2000 gives ~40–100 counts/bin (χ² valid). Pre-register M, L, B before the run (D-02/SBC-03). `[CITED: Talts et al. 2018 arXiv:1804.06788]`
- **Cost:** M×2 simulations for the paired path (~4000 `simulate_pair` calls) + M cheap forward passes. At imsize ~256² this is minutes on CPU; keep the *test-gate* M small (fixture) and the *reported* M ≈ 2000 in a script artifact (Phase-4 pattern: fast gate loads persisted results).

### Paired-draw Δρ SBC (D-01)
```
θ*_s, θ*_c = sample_prior(rng), sample_prior(rng)          # independent paired draw
Z_s, Z_c   = standardize_summary(encode_d01(patch_summary(build_mci(simulate_pair(...θ*_s)))), m.zt), ...
ρs = rho_draws(m.estimator, Z_s, m.θzt; N=L)               # L ρ_true draws (infer.jl)
ρc = rho_draws(m.estimator, Z_c, m.θzt; N=L)
Δρ_draws = ρs .- ρc                                        # MC difference (Phase-4 D-03)
Δρ*      = θ*_s.ρ_true − θ*_c.ρ_true
rank_Δρ  = count(Δρ_draws .< Δρ*)
```
`rho_draws` already exists in `infer.jl` — reuse it. The differenced draws are the correct posterior of Δρ under between-stack independence (the sample-vs-control design).

### Uniformity + coverage (SBC-02)
- **KS:** `ExactOneSampleKSTest((rank .+ 0.5) ./ (L+1), Uniform(0,1))` then `pvalue(test)`. `[CITED: HypothesisTests docs]`
- **χ²:** bin ranks into B equal bins → `ChisqTest(counts, fill(1/B, B))` (or `ChisqTest(counts)` against uniform expectation) → `pvalue`. Talts et al. recommend χ² on the rank histogram as the primary SBC test. `[CITED: Talts 2018; HypothesisTests docs]` `[ASSUMED: exact ChisqTest arg form — verify signature at implementation]`
- **Coverage curve (nominal vs empirical):** for a grid of central credible levels α, empirical coverage = fraction of the M draws whose θ\* falls inside the α-central interval of its L posterior draws; plot vs nominal α (should be diagonal). Can also use NeuralEstimators `assess(...)` + `coverage(assessment)` which computes MC coverage from `lower/truth/upper` columns (Hermans 2022 Def 2.1). `[VERIFIED: assess.jl:515-551]`
- **ECE/MCE traffic-light:** port `_bin_calibration`→`CalibrationResult` (ECE=Σ weight·|pred−obs|, MCE=max gap) from `~/Documents/GitHub/BayesInteractomics/src/diagnostics/calibration.jl` (template — copy the ~50-line function into `spike/validation/`, do NOT add BayesInteractomics as a dep). Map the SBC PIT/coverage into the (predicted-prob, empirical-positive) binning. Traffic-light cutoffs (green/yellow/red on ECE) pre-registered. `[VERIFIED: calibration.jl:22-72 read]`

## OOD / Misspecification (OOD-01/02)

### Summary-density channel (D-05)
- **Mahalanobis (recommended start):** fit mean μ_S and covariance Σ_S on the **standardized TRAIN summaries only** (`m.zt`-transformed). **Critical:** the 128-dim `:min` vector has 64 continuous rows (1:64) + 64 **binary mask** rows (65:128). The mask rows are near-constant (mostly all-1 for in-distribution non-degenerate images) → Σ singular. **Fit Mahalanobis on the continuous rows (1:64) only** (or the loader's `cont` partition), matching the loader's row split (`Loader._row_partition`). Score = (z−μ)'Σ⁻¹(z−μ) via `cholesky(Symmetric(Σ))`. `[VERIFIED: loader.jl:126-136 row partition; encode.jl mask layout]`
- **Normalizing-flow log-lik (escalation):** a NeuralEstimators `NormalisingFlow` trained on train summaries; score = −log q(Z). Use only if Mahalanobis ROC AUC is weak (D-05 "pick by separability"). Adds a trained artifact.
- Fit on TRAIN only, freeze, persist (JLD2). Quantile threshold pre-registered (D-06).

### Posterior-predictive mismatch channel (D-05)
- Infer θ̂ from Z_obs (`posterior_for` → posterior mean or draws), re-simulate P summaries from θ̂ via `simulate_pair`, compute a **discrepancy** between the observed summary and the posterior-predictive summary distribution. Candidate statistics (research): Mahalanobis of S_obs vs the PP-summary cloud; or a summary-space L2/energy distance; or per-feature z-score max. Rationale (D-05): catches inputs whose summary is marginally in-distribution but *inconsistent with the inferred θ* — reaches some summary-orthogonal cases the density channel misses.
- Template: `_ks_test_uniform` / `model_diagnostics` in `~/Documents/GitHub/BayesInteractomics/src/diagnostics/predictive_checks.jl` (copy pattern, not dep). `[VERIFIED: predictive_checks.jl:722-732 read]`

### Fusion + ROC (D-05/D-06)
- Report each channel's own ROC + AUC; combined flag fires if EITHER exceeds its pre-registered ID-quantile threshold (OR).
- **ROC/AUC:** hand-roll from sorted scores (ID negatives, misspec positives): sweep thresholds, compute TPR/FPR, AUC by trapezoid. ~15 lines, no new dep.
- **Operating point:** pre-registered ID quantile (e.g. 95th pct of ID scores → ~5% FPR), committed before seeing misspec data. **Youden-J** (max TPR−FPR) shown on the ROC as a labeled *post-hoc reference only*, never the gate.

### Misspecification grid (D-03 positive controls)
The simulator is a **shared-latent smooth-Gaussian-field** generator (`forward.jl`); OOD = things it cannot produce. Perturbations, per family:
1. **Texture-model (puncta/granular/fibrillar):** a NEW image generator (spot/blob mixture, fibre process) — NOT `simulate_pair`. Phase-2 explicitly rejected puncta as un-generatable, making it the canonical real-world OOD (real IF is punctate). Highest-value positive control.
2. **Noise-model:** corrupt `simulate_pair` output (or override stage 7) with non-Gaussian/salt-pepper/structured detector noise beyond the Poisson+Gaussian model.
3. **Optics/PSF:** apply an asymmetric/aberrated/out-of-focus PSF (anisotropic Gaussian, coma) or chromatic warp beyond the fixed `σ_psf=1.3` nuisance.
4. **Background/illumination:** add multiplicative illumination gradients, vignetting, autofluorescence bleed beyond `Uniform(0,0.1)`, bright debris blobs.
Sweep each family over a magnitude grid (research: 3–5 levels/family) so the ROC shows a dose-response. All applied to **fully external** images (SC5) — never the train/holdout set.

### Summary-orthogonal negative control (D-04) — the named blind spot
Apply, to in-distribution images, transforms that **provably/empirically preserve per-patch Pearson correlation**:
- **Per-channel affine `a·x+b`** (a>0): Pearson is exactly scale/shift invariant → 8×8 summary unchanged. Strongest theoretical guarantee.
- **Global rotation/flip:** permutes which patch is where but preserves the *set* of 8×8 correlation values → the summary *distribution* unchanged (KS on the 64 values passes) even though the vector is permuted. (Note: rotation by non-90° needs interpolation → verify empirically.)
- **Distribution-preserving spatial rearrangement:** block-permute patches; correlation values preserved as a set.
**Verify empirically (D-04):** KS test (`ExactOneSampleKSTest` or two-sample KS) on the 8×8 summary distribution before vs after transform; assert KS statistic < pre-registered ε. The OOD flag must stay **quiet** on these by construction — this *measures and names* the structural blind spot (OOD-02 honesty), paired with the SBC "under the simulator" caveat (SBC-04).

## Architecture Patterns

### System Architecture Diagram
```
                         trained_npe.jld2  (m.estimator, m.zt, m.θzt  — FROZEN train stats)
                                    │
          ┌─────────────────────────┼──────────────────────────────────────┐
          ▼                         ▼                                        ▼
   sample_prior(rng) ──► simulate_pair ──► build_mci ──► patch_summary ──► encode_d01 ──► standardize_summary(·, m.zt)
   (Random123 disjoint      (forward.jl)    (contract.jl, UNCHANGED src correlation 8×8)         │  = Z (128-dim)
    counter stream, D-02)                                                                        │
          │                                                                                      ▼
          │  θ*  (ground truth)                                              posterior_for(m.estimator, Z; N=L)  [one pass]
          │                                                                                      │  7×L std draws
          │                                                                    reconstruct(m.θzt, ·)  → physical θ
          ▼                                                                                      ▼
   ┌──────────────┐     ┌─────────────────────────┐     ┌──────────────────────────────────────────┐
   │ SBC (sbc.jl) │     │ Amortized BF (bf.jl +    │     │ OOD (ood.jl)                              │
   │ rank θ* among│     │ train_ratio.jl)          │     │ density: Mahalanobis(Z; μ_S,Σ_S train)    │
   │ L draws /7θ  │     │ RatioEstimator model-idx │     │ PP: discrepancy(S_obs, sim(θ̂))            │
   │ + Δρ paired  │     │ logBF = logratio(m=1)     │     │ ROC over misspec grid (D-03) + neg (D-04) │
   │ KS/χ²/cover  │     │        − logratio(m=0)    │     │ OR-fuse; ID-quantile op point             │
   │ ECE traffic  │     │ vs compute_BayesFactor    │     │ Youden-J post-hoc only                     │
   └──────┬───────┘     └───────────┬──────────────┘     └───────────────────┬──────────────────────┘
          └────────────────────────►│  figures.jl (CairoMakie) + pre-registered consts (test_validation.jl)
                                     ▼
                    spike/test/runtests.jl  (re-runnable gates, resolve-risk guard)
```

### Recommended Project Structure
```
spike/validation/          # NEW — sibling to spike/npe/
├── harness.jl             # shared draw→simulate→infer; load_frozen_model(); disjoint Random123 stream
├── sbc.jl                 # per-param + Δρ ranks, KS/χ², coverage curve, ECE/MCE (ported _bin_calibration)
├── train_ratio.jl         # RatioEstimator model-index training (sibling to npe/train_npe.jl)
├── bf.jl                  # amortized log-BF (logratio m=1 − m=0); baseline reproduction check
├── ood.jl                 # Mahalanobis + PP channels, misspec generators (D-03), neg control (D-04), ROC/AUC
├── figures.jl             # CairoMakie: rank hist, coverage, reliability/ECE, ROC, traffic-light
└── consts.jl              # (optional) shared pre-registered constants, or keep them in test_validation.jl
spike/test/
└── test_validation.jl     # SC1..SC5 gates + PRE-REGISTERED CONSTS block (mirror test_npe.jl), + wired into runtests.jl
```

### Pattern 1: Frozen-stats harness (reuse infer.jl surface)
**What:** every deliverable loads ONE model + its frozen transforms and never re-fits. **When:** all of SBC/BF/OOD.
```julia
# Source: composes existing spike/npe/infer.jl helpers (VERIFIED present)
m  = load_npe(joinpath(@__DIR__, "..", "npe", "trained_npe.jld2"))   # estimator, zt, θzt, d_in
Z  = standardize_summary(encode_d01(patch_summary(mci)), m.zt, :min) # frozen train zt
d  = posterior_for(m.estimator, Z; N = L)                            # one forward pass, use_gpu=false
θ  = StatsBase.reconstruct(m.θzt, d)                                 # physical θ (7×L)
```

### Pattern 2: Pre-registered consts + fresh disjoint seed (Phase-4 discipline)
**What:** all thresholds committed in `test_validation.jl` BEFORE the reported run; reported numbers come from a fresh Random123 counter range disjoint from training/holdout. **When:** SBC-03, D-02/D-06/D-08.
```julia
# Mirror spike/test/test_npe.jl:70-86 const block
const SBC_M          = 2000        # pre-registered draws (SBC-01/03)
const SBC_L          = 999         # posterior draws per SBC draw
const SBC_BINS       = 50          # rank-histogram bins (L+1 multiple of BINS)
const SBC_KS_ALPHA   = 0.05        # uniformity pass threshold
const OOD_ID_QUANTILE= 0.95        # pre-registered operating point (~5% ID FPR, D-06)
const BF_CORR_MIN    = 0.95        # D-08a rank/Pearson correlation gate
const BF_LOGBF_TOL   = 0.5         # D-08b bounded |Δ log-BF| over decision range
const VAL_MASTER_SEED= 0x5BC0FFEE  # DISJOINT from NPE_MASTER_SEED (D-02 fresh stream)
```
(Values illustrative — the planner/discuss step sets and locks the real numbers.)

### Anti-Patterns to Avoid
- **Re-fitting standardization or OOD nulls on eval data** — breaks SC5/D-06. Always use `m.zt` and train-only fits.
- **Fitting Mahalanobis on all 128 rows** — the 64 binary mask rows make Σ singular. Use continuous rows only.
- **Assuming `use_gpu=false`** — it defaults to TRUE in NeuralEstimators; pass `use_gpu=false` on EVERY `train`/`sampleposterior`/`logratio` call (CPU-only gate asserts no CUDA in `Base.loaded_modules`).
- **Passing a custom loss to `train(RatioEstimator, …)`** — silently ignored (loss hard-coded). If l-POP is truly required, hand-roll the estimator.
- **Adding a ROC/flow package** — risks the co-resolve landmine (NeuralEstimators→0.1.4). Hand-roll ROC; use the built-in `NormalisingFlow` only.
- **Editing `src/`** — reach `compute_BayesFactor`/`correlation` only via read-only `include()` (hard constraint until Phase 7).
- **Tuning M/thresholds to pass** — the data-snooping hazard the pre-registration exists to prevent (STATE.md Blockers).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Posterior sampling | Custom flow sampler | `sampleposterior` / `posterior_for` | Already wraps the flow; one pass |
| θ un-standardization | Manual (x·scale+mean) loops | `StatsBase.reconstruct(m.θzt, ·)` | Inverse of the frozen transform, matches training |
| KS / χ² uniformity | Manual CDF/statistic | `ExactOneSampleKSTest`, `ChisqTest` (HypothesisTests) | Pinned, exact p-values (CLAUDE.md mandate) |
| Interval / coverage | Manual quantiles | `interval`, `assess`+`coverage` | NeuralEstimators-native (infer.jl already uses `interval`) |
| ECE/MCE binning | Fresh binning code | Port `_bin_calibration`/`CalibrationResult` (BayesInteractomics) | Battle-tested pattern; traffic-light already defined |
| Ratio / BF classifier | Hand-rolled NRE from scratch | `RatioEstimator` + `logratio` (model-index) | The pinned, documented ratio API |
| Amortized log-BF integration | quadgk/KDE on amortized side | `logratio(m=1) − logratio(m=0)` | One forward pass; BF-02 forbids quadgk/KDE |

**Key insight:** almost the entire read side already exists in `spike/npe/infer.jl`; Phase 5 is *composition* of `posterior_for`/`rho_draws`/`delta_rho` + `HypothesisTests` + one new `RatioEstimator` training, not new inference machinery.

## Runtime State Inventory

> Phase 5 writes NEW artifacts; it does not rename/migrate. Included for completeness (no migration needed).

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | New JLD2 artifacts only: trained ratio net, frozen OOD nulls (μ_S/Σ_S), pre-registered consts, SBC/ROC result caches under `spike/validation/`+`spike/npe/` | write new; none migrated |
| Live service config | None — CPU-only local spike, no external services | None |
| OS-registered state | None | None |
| Secrets/env vars | None | None |
| Build artifacts | Consumes `spike/npe/trained_npe.jld2` (Phase-4) read-only; no rebuild of prior artifacts | None |

**Nothing found in categories 2–5:** verified — this is an additive, offline, CPU-only research phase.

## Common Pitfalls

### Pitfall 1: `use_gpu` defaults to TRUE
**What goes wrong:** a bare `sampleposterior`/`train`/`logratio` call tries CUDA. **Why:** NeuralEstimators default. **Avoid:** pass `use_gpu=false` everywhere (infer.jl already does). **Warning sign:** CUDA appears in `Base.loaded_modules` — the existing CPU-only assertions catch it.

### Pitfall 2: Singular covariance in Mahalanobis
**What goes wrong:** `cholesky` fails / NaN scores. **Why:** 64 binary mask rows are near-constant. **Avoid:** fit on continuous rows (1:64) only. **Warning sign:** non-PD covariance error at fit time.

### Pitfall 3: l-POP assumed present
**What goes wrong:** plan specifies l-POP loss that doesn't exist → dead end. **Why:** CLAUDE.md/CONTEXT.md say "l-POP" but v0.2.1 lacks it. **Avoid:** use NRE model-comparison; flag literal-l-POP as a user decision. **Warning sign:** searching for a `pop`/`evidence` symbol returns nothing (confirmed).

### Pitfall 4: SBC rank binning remainder
**What goes wrong:** χ² on rank histogram biased if (L+1) not divisible by bin count. **Why:** Talts et al. binning requirement. **Avoid:** pick L,B with (L+1) mod B == 0 (e.g. L=999, B=50 → 1000/50=20). **Warning sign:** non-uniform bin widths / spurious χ² rejection.

### Pitfall 5: prior-odds assumed 0
**What goes wrong:** subtracting log prior-odds=0 when it isn't → biased BF. **Why:** although μ prior is symmetric (odds≈1, log≈0), the *simulator* ρ-prior through ghat may not be perfectly symmetric. **Avoid:** compute prior-odds empirically from the training label balance / prior draws and subtract the measured value. **Warning sign:** amortized BF offset from KDE by a constant across the sweep.

### Pitfall 6: negative-correlation tail sparsity
**What goes wrong:** SBC ranks / OOD nulls unreliable at ρ<0. **Why:** Phase-2 D-16 / Phase-3 D-13 documented negative-ρ tail is prior-only-reachable, sparse in real anchors. **Avoid:** report SBC per-parameter and note the negative-Δρ regime as a documented caveat (do NOT stratify — distorts π). **Warning sign:** rank histogram spikes at extremes for ρ_true.

### Pitfall 7: co-resolve landmine on new deps
**What goes wrong:** adding ROC/flow package downgrades NeuralEstimators to 0.1.4. **Why:** Phase-1 documented co-resolve cap. **Avoid:** no new deps; hand-roll ROC; extend the resolve-risk gate. **Warning sign:** `NeuralEstimators v0.1.4` in a re-resolved Manifest.

## Code Examples

### Amortized log-BF in one forward pass (BF-02)
```julia
# Source: composes RatioEstimator.jl:160 logratio (VERIFIED v0.2.1)
# est_bf: RatioEstimator trained on paired summaries, model-index parameter m∈{0,1}
grid = reshape(Float32[0.0 1.0], 1, 2)          # two "parameter" configs: null, coloc
lr   = logratio(est_bf, Z_pair; grid = grid)     # 1×2 (one dataset): [logr(m=0) logr(m=1)]
log_prior_odds = log(p_coloc_prior / (1 - p_coloc_prior))   # measured from training labels (Pitfall 5)
log_BF = (lr[1,2] - lr[1,1]) - log_prior_odds    # amortized log Bayes factor — no quadgk/KDE/shuffle
```

### SBC rank uniformity (SBC-02)
```julia
# Source: HypothesisTests (pinned)
using HypothesisTests, Distributions
u   = (ranks .+ 0.5) ./ (SBC_L + 1)              # ranks ∈ 0:L → (0,1)
ks  = pvalue(ExactOneSampleKSTest(u, Uniform(0,1)))
cnts = [count(b -> (b-1) <= r*SBC_BINS/(SBC_L+1) < b, ...) for b in 1:SBC_BINS]  # rank histogram
chi = pvalue(ChisqTest(cnts))                    # vs uniform expectation (verify arg form at impl)
pass = ks > SBC_KS_ALPHA && chi > SBC_KS_ALPHA   # pre-registered thresholds
```

### Mahalanobis OOD score (OOD-01)
```julia
# Source: LinearAlgebra + StatsBase (stdlib/pinned). Fit on TRAIN continuous rows only.
using LinearAlgebra
cont, _ = Loader._row_partition(:min, 128)       # rows 1:64
μS  = vec(mean(Ztrain[cont, :]; dims=2))
ΣS  = cov(Ztrain[cont, :]; dims=2)
C   = cholesky(Symmetric(ΣS + 1e-6I))            # ridge for stability
maha(z) = (r = z[cont] .- μS; sum(abs2, C.L \ r))   # score; threshold at pre-registered ID quantile
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| KDE + quadgk Δρ Bayes factor (`src/bayes.jl`) | Amortized NRE model-comparison log-BF in one forward pass | this phase | The reproduction target (BF-02); baseline stays read-only |
| Per-dataset ADVI posterior | Amortized `sampleposterior` (Phase 4) | Phase 4 | SBC/BF/OOD all ride the amortized net |
| l-POP Evidence Networks (Jeffrey & Wandelt 2024) | Not in NeuralEstimators v0.2.1 | — | Standard Hermans-2020 NRE used instead; l-POP is an escape hatch only |

**Deprecated/outdated:** none in the pinned env; note NeuralEstimators 0.1.4 also on disk — ignore, spike pins 0.2.1.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | NRE model-index reduction (logratio m=1 − m=0) reproduces the baseline odds-ratio BF | Amortized BF | If user insists on literal l-POP, need a bespoke Evidence Network (larger scope) |
| A2 | Coloc/null split on sign(Δρ_ρ) ≡ baseline {Δμ>0} because ghat monotone through 0 | BF Null | If `ghat(0)≠0`, define split on induced-μ contrast instead (research check) |
| A3 | `ChisqTest(counts)` arg form for rank histogram | SBC Uniformity | Wrong arg form → compile error; verify HypothesisTests signature at impl |
| A4 | Mahalanobis on continuous rows gives usable separability | OOD density | May need `NormalisingFlow` log-lik escalation (D-05 allows) |
| A5 | Prior-odds ≈ but not exactly 1; measure empirically | Pitfall 5 | Constant BF offset if assumed 0 |
| A6 | Illustrative pre-registered const values (M=2000,L=999,etc.) | Patterns | Real values set/locked in discuss/plan step before reported run |
| A7 | Paired-summary encoding for the ratio net = concat(S_s,S_c) 256-dim | Amortized BF | A symmetric/contrast encoding may train better; research the encoding |

**All l-POP, ChisqTest-arg, ghat(0), and pre-registered-value items need confirmation before locking.**

## Open Questions

1. **Literal l-POP vs NRE model-comparison for BF-01.**
   - Known: v0.2.1 has no l-POP; `RatioEstimator` = logit-BCE NRE; model-comparison reduction is provably the same BF.
   - Unclear: whether the requirement/user mandates the literal l-POP loss.
   - Recommendation: use NRE model-comparison; surface as an explicit discuss-phase decision; keep hand-rolled l-POP as a documented fallback.

2. **Which BF baseline to correlate against (D-08).**
   - Known: `compute_BayesFactor` needs `CoLocResult` posteriors; the isolated Phase-4 baseline env produced ADVI posteriors.
   - Unclear: whether to reuse the ADVI artifact's Δρ posteriors or re-run the KDE odds-ratio on NPE Δρ draws for the sweep.
   - Recommendation: run the KDE odds-ratio (`compute_BayesFactor` logic) on a controlled Δρ sweep of re-simulated pairs so both sides see identical inputs (apples-to-apples, D-07).

3. **Paired-summary encoding for the ratio net.**
   - Known: BF is on a (sample,control) pair; the summary is per-stack.
   - Recommendation: start with 256-dim concat; consider a symmetric/difference encoding; research at plan time.

4. **PP-mismatch discrepancy statistic (D-05).**
   - Recommendation: Mahalanobis of S_obs vs PP-summary cloud as the default; compare separability against summary-density on the grid.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| NeuralEstimators.jl | SBC draws, ratio BF | ✓ | 0.2.1 (pinned) | — |
| Flux.jl | ratio net backend | ✓ | 0.16.10 | — |
| HypothesisTests.jl | KS/χ² SBC | ✓ | pinned | — |
| CairoMakie.jl | figures | ✓ | pinned | — |
| StatsBase / Distributions / Random123 / JLD2 | stats, priors, seeding, cache | ✓ | pinned | — |
| LinearAlgebra / Statistics / Test | Mahalanobis, ROC, gates | ✓ | stdlib | — |
| `spike/npe/trained_npe.jld2` | the net the harness infers with | ✓ (Phase-4 output) | — | re-run `train_npe.jl` at `NPE_MASTER_SEED` |
| `src/bayes.jl` `compute_BayesFactor` | BF reproduction target | ✓ (read-only include) | — | — |
| BayesInteractomics diagnostics | template for `_bin_calibration`/KS | ✓ (sibling repo, copy pattern) | — | write from scratch (~50 lines) |

**Missing dependencies with no fallback:** none. **With fallback:** trained_npe.jld2 (re-trainable); BayesInteractomics templates (re-implementable).

## Validation Architecture

> `workflow.nyquist_validation = true` (config.json). Included.

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (`@testset`/`@test`) |
| Config file | none — `spike/test/runtests.jl` includes per-phase test files |
| Quick run command | `julia --project=spike spike/test/runtests.jl` |
| Full suite command | same (single aggregated gate; add `test_validation.jl` include) |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SBC-01 | per-param + Δρ ranks produced (7×L, rank∈0:L) | unit | `julia --project=spike spike/test/runtests.jl` (SC1) | ❌ Wave 0 (`test_validation.jl`) |
| SBC-02 | KS/χ² p-values + coverage + ECE/MCE traffic-light finite | unit | same (SC2) | ❌ Wave 0 |
| SBC-03 | consts pre-registered; fresh disjoint seed | unit | same (assert const block + `VAL_MASTER_SEED≠NPE_MASTER_SEED`) | ❌ Wave 0 |
| SBC-04 | "under simulator" caption + OOD pairing present | doc | same (assert doc string) | ❌ Wave 0 |
| BF-01 | ratio net trained; log-BF = logratio(1)−logratio(0) finite, one pass | unit | same (SC3) | ❌ Wave 0 |
| BF-02 | corr ≥ BF_CORR_MIN AND max|Δ log-BF| ≤ BF_LOGBF_TOL vs `compute_BayesFactor`; no quadgk/KDE on amortized side | integration | same | ❌ Wave 0 |
| OOD-01 | Mahalanobis+PP scores; OR flag fires on misspec, quiet on ID | unit | same (SC4) | ❌ Wave 0 |
| OOD-02 | ROC/AUC over grid; neg-control KS<ε and flag quiet; blind spot named | integration | same (SC5) | ❌ Wave 0 |
| SC5/D-06 | all nulls fit on TRAIN only; misspec external | unit | assert fit inputs = train stats | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** quick gate on a SMALL fixture (tiny M/L, few grid points) — mirror `test_npe.jl` (persisted model loaded, not retrained; small `ABL_FIX_*`).
- **Per wave merge:** full `runtests.jl`.
- **Phase gate:** full suite green + the REPORTED SBC/BF/OOD numbers produced by a separate script at pre-registered M/L (Phase-4 pattern: fast gate loads persisted results; reported run is the artifact).

### Wave 0 Gaps
- [ ] `spike/test/test_validation.jl` — SC1..SC5 skipped scaffold + PRE-REGISTERED CONSTS block (mirror `test_npe.jl:70-96`), wired into `runtests.jl`.
- [ ] `spike/validation/harness.jl` — shared draw→simulate→infer + `load_frozen_model` + disjoint Random123 stream.
- [ ] Extend the `runtests.jl` resolve-risk gate to still assert NeuralEstimators==0.2.1 after any validation-code load.
- [ ] Fixture cache (small) reused across SBC/BF/OOD gates.

## Security Domain

> `security_enforcement` absent in config ⇒ treated as enabled. This is a CPU-only, offline, single-user scientific research spike with no network, auth, or persistence of untrusted data — most ASVS categories are N/A.

### Applicable ASVS Categories
| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | no users/sessions |
| V3 Session Management | no | — |
| V4 Access Control | no | local files only |
| V5 Input Validation | yes | θ/imsize already validated at `simulate_pair` entry (ArgumentError on out-of-range/non-finite); validation code must likewise reject bad grid/holdout indices (see `resimulate_holdout` throwing on bad index) |
| V6 Cryptography | no | no secrets; Random123 is for reproducibility, not security |

### Known Threat Patterns for this stack
| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Corrupt/half-written JLD2 artifact | Tampering | atomic `.tmp`+integrity-check+`mv` (existing `save_npe` idiom) |
| Silent stale-model / stale-stats reuse | Tampering | integrity-check estimator key on load; frozen-stats provenance from `trained_npe.jld2` meta |
| Data-snooping (tune-until-calibrated) | (scientific integrity) | pre-registered consts committed before the reported run; fresh disjoint seed (D-02) |
| Out-of-range perturbation crashing sim | DoS/robustness | reuse `simulate_pair` entry validation; guard misspec generators the same way |

## Sources

### Primary (HIGH confidence)
- `[VERIFIED]` Installed NeuralEstimators.jl v0.2.1 source (`~/.julia/packages/NeuralEstimators/gFxuZ/src/`): `Estimators/RatioEstimator.jl` (construction/loss/logratio/sampleposterior), `inference.jl` (sampleposterior/interval/posteriormean), `assess.jl` (assess/coverage/rmse), `train.jl` (train forms), `losses.jl` (no l-POP), `NeuralEstimators.jl` (exports).
- `[VERIFIED]` Spike code: `spike/npe/infer.jl`, `train_npe.jl`, `architecture.jl`; `spike/data/loader.jl`, `encode.jl`; `spike/simulator/prior.jl`, `forward.jl`; `spike/test/test_npe.jl` (pre-registration pattern); `spike/Project.toml`.
- `[VERIFIED]` `src/bayes.jl:109-136` `compute_BayesFactor` body (BF null confirmed).
- `[VERIFIED]` `~/Documents/GitHub/BayesInteractomics/src/diagnostics/calibration.jl` (`_bin_calibration`/`CalibrationResult`), `predictive_checks.jl` (`_ks_test_uniform`), `types.jl` (CalibrationResult struct).
- `[CITED]` NeuralEstimators.jl dev docs nav + methodology (msainsburydale.github.io/NeuralEstimators.jl/dev) — RatioEstimator/PosteriorEstimator, coverage (Hermans 2022).

### Secondary (MEDIUM confidence)
- `[CITED]` Talts et al. 2018, "Validating Bayesian Inference Algorithms with SBC" (arXiv:1804.06788) — rank-histogram + χ²/KS uniformity, binning.
- `[CITED]` Hermans et al. 2020 (NRE, linked in RatioEstimator docstring); Jeffrey & Wandelt 2024 (Evidence Networks / l-POP — for the fallback only).

### Tertiary (LOW confidence)
- `[ASSUMED]` Exact `ChisqTest` argument form; optimal paired-summary encoding for the ratio net; `ghat(0)` value — all flagged for verification at implementation.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — every tool verified in the pinned env; no new deps.
- SBC mechanics: HIGH — composes existing verified `infer.jl` surface + pinned HypothesisTests.
- Amortized BF: MEDIUM-HIGH — API verified; the l-POP gap and model-index reduction are the key judgment calls (flagged as assumptions).
- OOD: MEDIUM — approach is standard; density-vs-flow and PP statistic left to separability (per D-05).
- Pitfalls: HIGH — drawn from read source + Phase-1..4 documented landmines.

**Research date:** 2026-07-02
**Valid until:** 2026-08-01 (stable — pinned Manifest; API frozen at 0.2.1)
</content>
</invoke>
