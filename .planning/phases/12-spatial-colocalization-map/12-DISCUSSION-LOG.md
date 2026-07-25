# Phase 12: Spatial Colocalization Map (GP/CAR) - Discussion Log

> **Audit trail only.** Not consumed by researcher, planner, or executor agents — decisions live in
> `12-CONTEXT.md`. This log preserves the options presented and the reasoning offered.

**Date:** 2026-07-25
**Mode:** discuss (default, interactive)
**Areas discussed:** Spatial model + amortization · Spatially-varying simulator · Real-image coverage
claim (SC3) · Descope trigger and v2.1 fallback (all four selected by the user)

---

## Pre-discussion scouting

The defining constraint was surfaced before any question, from `local_map.jl:27-32` citing
`07-RESEARCH` Finding 2: *"a finer patch grid lengthens the summary vector but the NPE output stays a
single GLOBAL ρ, and windowing does not add a spatial prior, borrow strength between neighbours, or
propagate per-region uncertainty."* Phase 12 therefore needs both a different output head and a
spatially-varying simulator — it is not reachable by refining the grid.

Also noted up front: `SpatialColocResult` is already sketched at `src/results.jl:181-191`, so the
result *type* is decided (Phase 7 D-02) and Phase 12 fills it.

---

## Area 1: Spatial model + amortization

### Q1. How does the network emit a lattice-valued posterior?

Framing given: in SBI the lattice prior lives in the **simulator**, not the head — so "CAR vs GP" is
Area 2's question and Area 1's is how the network *emits* a lattice.

- **Hierarchical: global ρ + deviation field** *(recommended)* — global term stays the shipped
  quantity; calibration tractable (global + low-dim deviation summary, not 64 marginals); shrinkage IS
  the borrow-strength mechanism.
- Joint flow over all R×C dims — fully expressive, but 64 marginals against an SBC apparatus that
  strained at 8 parameters with 5 vacuous.
- Per-region conditional over region index — cheap and scalable, but no joint posterior, which is the
  independent-pooling baseline SC3 requires beating.

**Selected: Hierarchical.** → D-01

### Q2. How is the summary treated?

> **Correction issued mid-discussion.** An earlier framing in this same session claimed Phase 12 must
> spend the fixed-summary constraint. It must not. The summary is `2·G² = 128` dims at G=8 — 64
> continuous + 64 mask rows, **one entry per patch** — so per-region information is already present.
> Flattening discards the spatial *inductive bias*, not the *information*. SC1's CNN is an architecture
> choice, not an information necessity.

Also noted: SC1 offers "CNN/**DeepSet**", but a DeepSet is permutation-invariant over patches —
precisely the exchangeable pooling the phase goal says it replaces. Only the CNN half is coherent.

- **Reshape to G×G×2 and use a CNN** *(recommended)* — content byte-identical, so `patch_summary`, the
  OOD Mahalanobis reference and `DATAGEN_HASH_SRC_FILES` stay untouched, while gaining the locality bias.
- Keep flattened, structure in head only — identical input contract, but no inductive bias at all.
- New CNN summary over the raw grid — most faithful to SC1's text, but changes summary content and
  invalidates the shipped read surface.

**Selected: Reshape to G×G×2 and use a CNN.** → D-02

### Q3. What should the SC1 mini-spike decide?

Coupling raised: the low-rank hierarchical head bounds how rough a field it can represent, so a
full-rank CAR field would be systematically under-fit by it.

- Reduced: fix GP, spike rank + correlation length *(recommended)* — avoids reporting a CAR-vs-GP
  verdict that actually measures the head.
- **Full CAR-vs-GP comparison as SC1 specifies** — satisfies SC1 verbatim; CAR is lattice-native and
  might genuinely win.
- Skip the spike entirely — discards the one decision ROADMAP flagged as needing evidence.

**Selected: Full CAR-vs-GP comparison** (against the recommendation). Rationale acknowledged: given
Phase 7's amendment history, keeping a pre-registered criterion intact has real value. The confound was
then handled by the next question rather than accepted. → D-03

### Q4. How is the comparison de-confounded from the low-rank head?

- **Spike with a high-rank head, pick rank after** *(recommended)* — makes the comparison measure the
  prior; sequences the two decisions cleanly.
- Match a head to each prior — fairest head-to-head, but confounds in the opposite direction.
- Accept and name the confound — a CAR loss would be a head artifact reported as a prior verdict, the
  same misreading class that cost the Phase-7 SBC arm its credibility.

**Selected: Spike with a high-rank head, pick rank after.** → D-04

---

## Area 2: Spatially-varying simulator

Binding constraint stated: `prior.jl:71-82` draws `μ* ~ Truncated(Cauchy(0,0.3),-1,1)` — verbatim the
Turing μ-prior — and sets `ρ_true = ghat(μ*)`, so the induced per-patch correlation mean matches by
construction (SIM-02). CLAUDE.md makes π(θ)↔Turing consistency a hard constraint.

### Q1. How is the ρ field constructed?

- **Marginal-preserving copula** *(recommended)* — correlated Gaussian field → Φ → `MU_PRIOR` quantile
  → `ghat` elementwise, so every region's marginal μ is **exactly** `MU_PRIOR` and the lattice prior
  contributes only dependence.
- Field directly on ρ — simplest and stated in the units of interest, but breaks SIM-02 and ADVI
  comparability outright.
- Global μ* + zero-mean deviation field — mirrors the head most literally, but marginals become a
  convolution so SIM-02 holds only approximately.

**Selected: Marginal-preserving copula.** → D-05

### Q2. How is the field applied in the forward model?

- **Smooth interpolation to pixel resolution** *(recommended)* — no grid-aligned discontinuities for the
  network to exploit.
- Piecewise-constant per patch cell — unambiguous ground truth, but hard edges aligned exactly with the
  summary grid, so the net could learn block-edge cues absent from real images: passing SC2 in
  simulation and failing SC3 on real data.
- Continuous GP per pixel — most principled but substantially more expensive for essentially the same
  field.

**Selected: Smooth interpolation.** → D-06

### Q3. What is the per-region ground truth?

`prior.jl:30-32` warns *"Do NOT assert ρ_true == μ"*, so three candidate truths exist and are not the
same quantity.

- **The drawn G×G lattice values** *(recommended)* — truth is literally the sampled θ, which is what SBC
  ranks; interpolation is a rendering detail.
- Cell mean of the interpolated field — intuitive, but inserts a deterministic transform between what
  was sampled and what is scored, the structure that produced the ρ_true atom artifacts.
- Induced per-cell patch correlation — faithful to SIM-02 philosophy, but a measured noisy quantity
  rather than a parameter.

**Selected: The drawn G×G lattice values.** → D-07

### Q4. How is the spatial correlation length treated?

Distinction drawn from Phase 11 D-03 / Phase 13 D-03: registration uncertainty is something a **user
knows** (their instrument); the spatial correlation length of biology is not, so the conditioning answer
does not carry over.

- **Infer it as part of θ** *(recommended)* — otherwise the calibration claim is self-confirming, since
  matched simulations drawn at the assumed length would show good coverage regardless.
- Fix at pre-registered values — simplest, but bakes an unfalsifiable assumption into the phase's
  headline claim.
- Condition on it as an input — consistent with the sibling phases, but asks for a number the caller has
  no basis to supply.

**Selected: Infer it as part of θ.** → D-08

---

## Area 3: Real-image coverage claim (SC3)

Two problems with SC3 as written were laid out: (a) there is no ground-truth ρ field on a real image, so
"coverage" has no referent; (b) the independent-pooling baseline cannot do coverage at all —
`LocalColocMap` (`local_map.jl:56-79`) states in its own docstring that it "carries no posterior draws,
no Bayes factor and no per-region uncertainty".

### Q1. How is SC3 made measurable?

- **Leave-region-out predictive coverage** *(recommended)* — needs no ground truth; makes "beats
  independent pooling" well-defined since independent pooling cannot predict a held-out region beyond
  the prior; measures the borrowing mechanism directly.
- Semi-synthetic imposed field on real images — real texture with exact known truth, but the field is a
  construction rather than biology.
- Both — strongest evidence, two harnesses and two baselines.

**Selected: Leave-region-out predictive coverage.** SC3 amended, pre-declared. → D-09

### Q2. What is the "independent pooling" baseline?

- **Ablation: spatial prior neutralized** *(recommended)* — same net, same training, same summary, so a
  win is attributable to the spatial prior alone.
- Per-tile independent NPE reads + posterior SD — a product-level claim a reader understands, but
  confounded by a different net reading smaller windows.
- Both — answers the scientific and product questions separately.

**Selected: Ablation.** → D-10

### Q3. Which real images?

Corpus audit presented: only two `physical-primary` rows exist (`coloc`, and `segregated`/mitochondria);
the 30 CBS Red-Green rows are `simulated-secondary` — computer-generated, so they cannot support a "real
image" claim. Noted that coverage is scored per region, so at G=8 each image contributes ~64 regions and
the effective sample is regions, not images.

- **Both physical-primary anchors** *(recommended)* — ~128 regions across two different spatial regimes;
  a prior that helps on coloc but not segregation would be an important finding.
- Coloc anchor only — meets SC3 literally, but silent on segregated specimens where a smoothness prior
  is most likely to over-smooth.
- Physical anchors plus CBS as supplementary — more images, but CBS is simulated and needs careful
  labelling.

**Selected: Both physical-primary anchors.** → D-11

---

## Area 4: Descope trigger and v2.1 fallback

Framing: ROADMAP names Phase 12 the natural descope candidate but specifies no trigger — so in practice
the decision would be made late, under sunk-cost pressure. The CAR-vs-GP mini-spike (D-03) runs before
the expensive training and already computes coverage against the ablation.

### Q1. What triggers the descope?

- **Two-stage: spike gate, then calibration gate** *(recommended)* — fails fast at the cheap checkpoint;
  the spike already produces the needed numbers.
- Single post-training calibration gate — decided on the real model, but after the largest cost is sunk.
- Budget-based ceiling — objective, but measures effort rather than whether the approach works.

**Selected: Two-stage.** → D-12

### Q2. What ships on descope?

Noted that a fallback product already exists in the design: the ablation model (D-10) is built regardless
as the gate baseline, and produces per-region Δρ *with* uncertainty — strictly more than `LocalColocMap`
offers today.

- **Ship the ablation model** *(recommended)* — converts the descope from a write-off into a reduced
  release at no extra cost.
- Ship nothing new, document the negative result — clean and legitimate, but leaves users where Phase 7
  left them.
- `local_coloc_map` + per-tile uncertainty — most conservative, but extra work since the ablation rather
  than per-tile reads was chosen as the SC3 baseline.

**Selected: Ship the ablation model.** → D-13

---

## Corrections issued during this discussion

1. **The fixed-summary constraint does not have to be spent** (Area 1 Q2). The 128-dim summary already
   carries one entry per patch, so spatial information is present; only the inductive bias is lost by
   flattening. This changed the Q2 options materially and is recorded in D-02.

## ROADMAP clauses corrected

- **SC1's "CNN/DeepSet"** — the DeepSet half is incoherent with the phase goal (permutation-invariance
  over patches is the exchangeable pooling being replaced).
- **SC3's "coverage on ≥1 real image"** — undefined as written; amended to leave-region-out predictive
  coverage (D-09) with the amendment pre-declared.

## Deferred ideas raised during discussion

New CNN summary over the raw grid · DeepSet summary · joint 64-dim flow · per-region conditional head ·
per-tile independent baseline as a product comparison · semi-synthetic imposed-field parameter coverage ·
CBS benchmark as supplementary · head-matched-to-prior spike variant · productionizing the spatial map.

All carried into `12-CONTEXT.md` `<deferred>`.
