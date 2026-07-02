# Phase 5: Validation Bundle (SBC + Amortized BF + OOD) - Context

**Gathered:** 2026-07-02
**Status:** Ready for planning

<domain>
## Phase Boundary

The publishable trifecta — a **calibration proof (SBC)**, an **amortized Bayes factor
(`RatioEstimator`/Evidence-Network)**, and an **honest OOD/misspecification flag** — built as
**sibling plans off ONE shared θ\*~π→simulate→infer harness** over the Phase-4 trained nets
(`spike/npe/trained_npe.jld2`). This is the **highest-scrutiny phase for scientific honesty**;
discussion clarified the honesty-critical implementation choices within this fixed boundary.

In scope: the shared SBC harness (rank histograms, KS/χ² uniformity, coverage curve, ECE/MCE
traffic-light via the `_bin_calibration` pattern) over all 7 θ **and** the Δρ contrast
(SBC-01..04); the amortized log-BF via a `RatioEstimator` (l-POP loss) reproducing
`compute_BayesFactor()` in the well-specified regime (BF-01/02); the OOD flag as a controlled ROC
over a misspecification grid with positive **and** summary-orthogonal negative controls plus a
posterior-predictive channel (OOD-01/02). **Out of scope:** the demo + Go/No-Go memo (Phase 6);
any `src/` edits or productionization (Phase 7); the external physical corpus (Phase 8) and
cross-method comparator (Phase 9). **No `src/` edits** — the harness reaches `src/`
(`compute_BayesFactor`, `correlation`) only via the read-only `include()` coupling.

Requirements: SBC-01, SBC-02, SBC-03, SBC-04, BF-01, BF-02, OOD-01, OOD-02.
</domain>

<decisions>
## Implementation Decisions

### SBC coverage scope (SBC-01, SBC-02)
- **D-01:** **SBC covers all 7 per-stack θ AND a dedicated Δρ-contrast SBC.** Rank histograms for
  the 7 sampled parameters (incl. ρ_true) are computed directly from `sampleposterior` marginals.
  In addition, a **dedicated paired-draw Δρ SBC**: draw a *paired* θ\* (sample θ\*_s, control θ\*_c),
  simulate both, and rank the true Δρ\* = ρ\*_s − ρ\*_c among **differenced posterior draws** (the
  MC difference of two independent single-stack passes, per Phase-4 D-03). Rationale: Δρ is the
  quantity the Bayes factor and the actual colocalization verdict rest on — calibrating it
  *directly* (not inheriting it by argument from ρ_true) closes the "you calibrated ρ but reported
  Δρ" gap. Costs an extra paired harness path.

### SBC anti-snooping discipline (SBC-03, SBC-04) — carried forward, not re-decided
- **D-02:** **M and pass/fail thresholds pre-registered; reported from a fresh, independently-seeded,
  never-tuned-against held-out run** — following the Phase-4 pre-registered-consts + keyed
  holdout-repro-gate pattern. M ≈ 2000 (exact value + L posterior draws/params = research). SBC is
  reported **"calibrated under the simulator"** and **explicitly paired with the OOD result** in all
  reporting (SBC-04) — SBC alone is not a real-data guarantee.

### Misspecification grid — positive controls (OOD-02)
- **D-03:** **All four positive-control families in the ROC grid (should-fire cases):**
  (1) **Texture-model mismatch** — puncta / granular / fibrillar spatial textures the shared-latent
  smooth-field simulator never produces (Phase-2 explicitly rejected puncta as un-generatable,
  which makes it the natural real-world OOD case; real IF images are punctate);
  (2) **Noise-model mismatch** — non-Gaussian / Poisson-shot / salt-pepper / structured detector
  noise vs the simulator's noise model;
  (3) **Optics/PSF mismatch** — asymmetric / out-of-focus / aberrated PSF, or chromatic warp beyond
  the simulator's fixed nuisance range;
  (4) **Background/illumination** — uneven illumination gradients, vignetting, autofluorescence
  bleed beyond prior range, bright artifacts/debris.
  Exact perturbation magnitudes + grid resolution per family = research/discretion.

### Misspecification grid — summary-orthogonal negative control (OOD-02, the named blind spot)
- **D-04:** **Negative control = correlation-preserving transforms.** Construct the
  summary-orthogonal perturbation from transforms provably/empirically invariant to per-patch
  Pearson correlation — per-channel **affine intensity rescale (a·x+b)**, global **rotation/flip**,
  and distribution-preserving spatial rearrangement. **Verify empirically** (KS test on the 8×8
  summary distribution) that the summary is statistically unchanged, so the OOD flag stays quiet
  **by construction** exactly where the fixed summary is blind. This *measures and names* the
  structural blind spot (discrepancies orthogonal to the fixed summary) rather than hiding it —
  the OOD-02 honesty requirement.

### OOD score fusion (OOD-01)
- **D-05:** **Report both channels separately AND an OR-combined flag.** The **summary-density
  channel** (Mahalanobis distance or normalizing-flow log-likelihood on the 8×8 summary vector) and
  the **posterior-predictive mismatch channel** (does inferred θ regenerate images whose summaries
  match the observed) each get their own ROC; the **combined flag fires if EITHER exceeds its
  threshold.** Rationale: diagnostic transparency (you see which channel caught what) + the OR
  captures complementary failure modes — summary-density catches summary-visible OOD, PP reaches
  some summary-orthogonal cases the density channel is blind to. Density-vs-flow choice = research
  (pick by separability); PP discrepancy statistic = research.

### OOD threshold discipline (OOD-01, SC5)
- **D-06:** **Threshold-free ROC/AUC headline + pre-registered ID-quantile operating point +
  Youden-J as post-hoc reference only.** Report the full **ROC curve + AUC** as the headline
  (score separability across all thresholds). The **reported operating point** is a
  **pre-registered in-distribution quantile** (e.g. ~5% ID false-positive rate) with the
  density/flow + PP nulls **fit on the TRAIN split only** and the quantile committed **before**
  seeing any misspecified data. The **Youden-J optimum** is shown on the ROC as a **descriptive
  post-hoc reference only — never the gating threshold** (clearly labeled), so best-case
  separability is visible without weakening the anti-snooping stance. ROC evaluated on **fully
  external misspecified test images** (SC5). Exact ID quantile = pre-registered before the run.

### Amortized Bayes factor — null definition (BF-01)
- **D-07:** **Read `compute_BayesFactor()` first, then mirror its null + prior EXACTLY.** Before
  defining the `RatioEstimator`'s coloc-vs-null hypotheses, read what the existing baseline actually
  tests on Δρ (`src/bayes.jl:109`, `compute_BayesFactor(posterior, prior; ρ_threshold=0.0)` — the
  `ρ_threshold=0.0` default hints at a Δρ≤0 / point-null-style test; **confirm by reading the
  function body**). Define the amortized null **identically** — same Δρ, same prior alignment
  (BF-01 mandates this) — whatever the baseline's null turns out to be. Guarantees BF-02 is an
  apples-to-apples *reproduction*, not two different tests that merely correlate.

### Amortized Bayes factor — agreement criterion (BF-02)
- **D-08:** **Pre-register BOTH correlation AND bounded log-BF error.** "Reproduces
  `compute_BayesFactor()` in the well-specified regime" is pre-registered as: (a) high rank/Pearson
  **correlation** between amortized and KDE log-BF across a held-out Δρ sweep, AND (b) **bounded
  absolute error in log-BF units** within a pre-set tolerance over the decision-relevant Δρ range.
  Analogous to Phase-4's pre-registered RMSE tolerance (D-09); catches both wrong ordering and
  right-ordering/wrong-magnitude. Exact correlation threshold + log-BF error tolerance + Δρ range =
  pre-registered before the reported run. Computed **without quadgk/KDE/shuffle** on the amortized
  side (BF-02).

### Claude's Discretion / Research
- Exact M (≈2000) and L posterior draws per SBC draw; exact KS/χ² p-value thresholds and ECE/MCE
  traffic-light cutoffs (all pre-registered, Phase-4 discipline).
- Summary-density channel: Mahalanobis vs normalizing-flow log-likelihood (choose by separability);
  flow architecture if used. Posterior-predictive mismatch statistic definition.
- Exact perturbation magnitudes / grid resolution per misspecification family (D-03); exact
  correlation-preserving transform parameterization + the KS-invariance ε (D-04).
- Exact ID false-positive quantile (D-06); exact BF correlation threshold, log-BF error tolerance,
  and Δρ sweep range (D-08) — all pre-registered before the reported run.
- Shared-harness module layout within `spike/` (e.g. `spike/validation/` alongside `spike/npe/`);
  figure set (CairoMakie: rank histograms, coverage curves, reliability/ECE, ROC curves,
  traffic-light).
</decisions>

<specifics>
## Specific Ideas

- **One shared θ\*~π→simulate→infer harness** feeds all three deliverables as sibling plans (roadmap
  goal) — SBC, BF, and OOD reuse the same draw→simulate→infer chain over the Phase-4 nets.
- **The blind spot is a feature of the honesty story, not a bug to hide** — D-04 constructs the
  summary-orthogonal negative control specifically so the paper can *name and measure* what the
  fixed 8×8 summary cannot see, paired with the SBC "under the simulator" caveat (SBC-04).
- **Everything frozen from the training split only** (SC5): OOD covariance/flow/PP nulls and all
  standardization are fit on TRAIN; misspecified test images are fully external.
- The `ρ_threshold=0.0` default on `compute_BayesFactor` is the anchor for the BF null (D-07) —
  read the body to confirm the exact hypothesis before mirroring.
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 5: Validation Bundle (SBC + Amortized BF + OOD)" — goal + 5 success
  criteria (SBC rank histograms/KS-χ²/coverage/ECE-MCE traffic-light; pre-registered M+threshold on
  a fresh held-out run, reported "under the simulator" paired with OOD; amortized log-BF reproducing
  `compute_BayesFactor()` without quadgk/KDE/shuffle; OOD ROC with pos + summary-orthogonal neg
  controls + posterior-predictive channel, blind spot measured; all preprocessing frozen from train).
- `.planning/REQUIREMENTS.md` — SBC-01, SBC-02, SBC-03, SBC-04, BF-01, BF-02, OOD-01, OOD-02.

### The nets + harness inputs (Phase-4 / Phase-3 outputs)
- `spike/npe/trained_npe.jld2` — the trained `PosteriorEstimator` (MLP → `NormalisingFlow` over
  7-dim θ) the harness infers with; `spike/npe/train_npe.jl` — its training entry (l-POP
  `RatioEstimator` for BF likely lands as a sibling here).
- `spike/data/loader.jl` (`Loader`) — the sole leak-free standardization path; **its train-split
  stats are what D-06's OOD nulls and all preprocessing must be frozen from** (SC5).
- `spike/data/encode.jl` — `encode_d01` (128-dim :min, the chosen summary per Phase-4 D-07); the
  8×8-derived vector the summary-density OOD channel scores.
- `spike/simulator/` — `sample_prior`, `simulate_pair`, `ghat`: the θ\*~π→simulate half of the
  shared harness (incl. the paired-draw Δρ path D-01 and the misspecification perturbations D-03/D-04).
- `spike/contract.jl` — `build_mci`, `patch_summary` (8×8 `correlation(...;:pearson)`, UNCHANGED
  `src/`), `induced_mu` — the read-only summary coupling.

### The BF baseline being reproduced (READ-ONLY — never edit in spike)
- `src/bayes.jl` — `compute_BayesFactor(posterior, prior; ρ_threshold=0.0)` (line 109) — the KDE
  Bayes factor the amortized log-BF must reproduce (BF-02); **read the body to fix the exact Δρ null
  + prior the `RatioEstimator` mirrors (D-07)**. `CoLocResult` (line 81). Reached via read-only
  `include()`.
- `src/colocalization.jl` — `correlation(x,y; method)`, `patch()`, `_exclude_zero` — the fixed 8×8
  summary the OOD blind-spot analysis (D-04) is defined relative to.

### Calibration reference implementation (sibling project — reuse the pattern)
- `~/Documents/GitHub/BayesInteractomics/src/diagnostics/calibration.jl` — `_bin_calibration` /
  `CalibrationResult` (ECE/MCE + traffic-light verdict) for SBC-02; `.../types.jl` — the result
  types; `.../predictive_checks.jl` — `_ks_test_uniform` / `model_diagnostics` templates for the
  SBC KS-uniformity + posterior-predictive channel (D-05). (Sibling repo, template only — not a
  spike dependency.)

### Prior-phase decisions this phase builds on
- `.planning/phases/04-npe-training-advi-benchmark-ablation/04-CONTEXT.md` — D-03 (Δρ = MC
  difference of two single-stack passes, the basis for D-01's paired Δρ SBC), D-07 (`:min` summary
  chosen, OOD-detectability the coupling), D-09 (pre-registered-tolerance pattern reused in D-08),
  pre-registered-consts + keyed holdout-repro-gate discipline (reused in D-02).
- `.planning/phases/03-training-data-pipeline/03-CONTEXT.md` — D-07/D-08 (loader fit-on-train-only,
  the freeze source for SC5/D-06), D-10 (reserved holdout discipline), D-11 (Random123 counter-based
  seeding for the fresh independently-seeded SBC run, D-02).

### Constraints & stack
- `CLAUDE.md` §Constraints (spike decoupling — no `src/` edits until Phase 7; CPU-only baseline;
  reproducible-from-fixed-seed; fixed 8×8 summary) + §Key API Notes (`RatioEstimator` = amortized
  log-BF in one forward pass; **l-POP loss**; `sampleposterior`/`assess`; SBC ranks hand-rolled then
  tested with **HypothesisTests** KS/χ²; reuse `_bin_calibration`/`CalibrationResult`) + §"What NOT
  to Use".
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `spike/npe/trained_npe.jld2` + `spike/npe/train_npe.jl` — the trained NPE the harness infers with;
  the training entry is the natural sibling home for the l-POP `RatioEstimator` (BF).
- `spike/simulator/` (`sample_prior`, `simulate_pair`, `ghat`) — drives every θ\*~π→simulate path:
  standard SBC draws, the paired Δρ draws (D-01), and the misspecification perturbations (D-03/D-04).
- `spike/data/loader.jl` train-split stats — the single freeze source for OOD nulls + preprocessing
  (SC5 / D-06).
- `spike/data/encode.jl` `encode_d01` — the 8×8-derived :min summary the summary-density OOD channel
  scores and the blind-spot control (D-04) is defined against.
- `src/bayes.jl` `compute_BayesFactor` / `CoLocResult` — the BF reproduction target (read-only).
- BayesInteractomics `src/diagnostics/calibration.jl` + `predictive_checks.jl` — ECE/MCE
  traffic-light + KS-uniformity + posterior-predictive templates (D-05, SBC-02).

### Established Patterns
- Spike reaches `src/` only via read-only `include()`/copy (Phase-1 D-01); root byte-identical to
  baseline `f581d95`. Same discipline for `compute_BayesFactor` here.
- Re-runnable quantitative gates under the spike `Test` harness (`spike/test/`, e.g.
  `test_npe.jl`) — SBC/BF/OOD land as re-runnable gates there with **pre-registered consts committed
  before the reported run** (Phase-4 pattern, D-02/D-06/D-08).
- Random123 counter-based reproducibility (Phase-3 D-11) — the fresh independently-seeded SBC run
  occupies its own disjoint counter range.

### Integration Points
- Trained Phase-4 NPE + `:min` summary → the shared θ\*~π→simulate→infer harness → SBC, BF, OOD as
  three sibling plans.
- Loader train-split stats → frozen OOD nulls / preprocessing (SC5).
- `compute_BayesFactor` (read-only) → BF-02 reproduction check.
- Phase-5 outputs (SBC verdict + amortized BF + OOD flag) → Phase-6 demo + Go/No-Go memo (SBC
  "under the simulator" framed with the OOD result and the Phase-4 speedup nuance).
</code_context>

<deferred>
## Deferred Ideas

- **Demo + Go/No-Go memo** — Phase 6 (this phase produces the metrics the memo reports).
- **Productionization / `src/` edits / user-definable `num_patches`** — Phase 7 (conditional on Go).
- **External physical ground-truth corpus** — Phase 8; **cross-method comparator (Costes-p etc.)** —
  Phase 9; both parallelizable and out of scope here. OOD/SBC here is "under the simulator"; external
  reality checks live in Phases 8/16.
- **Three-hypothesis BF (coloc/random/mutually-exclusive)** — Phase 13; this phase is the two-way
  coloc-vs-null BF only.
- **DeepSet summary upgrade / adversarial nuisance sweep + CI gate** — Phase 15; here the OOD grid is
  the honest first pass, not the exhaustive operating-envelope sweep.

### Reviewed Todos (not folded)
None — no pending todos matched this phase.
</deferred>

---
*Phase: 05-validation-bundle-sbc-amortized-bf-ood*
*Context gathered: 2026-07-02*
