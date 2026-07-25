# Phase 13: Three-Hypothesis Amortized Bayes Factor — Research

**Researched:** 2026-07-25
**Domain:** Amortized neural model comparison (NRE / BCE-head evidence networks) in Julia;
NeuralEstimators v0.2.1 + Flux 0.16.10; simulation-based validation; semi-synthetic image ground truth
**Confidence:** HIGH on API surface, repo facts, prior masses, and the D-07 derivation.
MEDIUM on the τ probe design (design proposed, not yet measured) and the α-series behaviour.
LOW on anything downstream of Phase 11 (its research net does not exist yet).

---

<user_constraints>
## User Constraints (from 13-CONTEXT.md)

### Locked Decisions (verbatim from `13-CONTEXT.md` `<decisions>`)

**Retraining posture and lane**

- **D-01: Research lane.** Train the 3-way evidence net in `spike/` on a **fresh DEV seed** asserted
  disjoint from `PROD_SEED_V2`, `VAL_MASTER_SEED` (`0x5BC0FFEE`), `NPE_MASTER_SEED`
  (`0xC0FFEE`), `RATIO_PAIR_SEED` (`0x00000000004A7107`), `VAL_FIX_SEED` (`0xF1F7ED`),
  `CORPUS_MASTER_SEED` (`0x0000000000c05eed`) and all prior dev seeds. The shipped artifact, the
  shipped `bayes_factor` accessor, and the Phase-7 GO stay untouched.
  *Noted during discussion:* a reship would be materially cheaper here than in Phase 11/12 — Phase 13
  touches neither the NPE, the summary, nor any `DATAGEN_HASH_SRC_FILES` entry — but research lane
  was chosen to keep the GO intact.

- **D-02: Train on Phase-11's registration-aware `zt`,** not the shipped grid-8 `zt`.
  `train_ratio.jl:144` trains on pairs assembled from the FROZEN `zt` of an already-trained NPE, so
  the basis is inherited rather than created.
  **⚠ CONSEQUENCE FOR PLANNING: this adds a Phase 11 dependency that ROADMAP does not declare**
  (ROADMAP lists Phase 7 only). Phase 13 *planning* may proceed now; **execution blocks on Phase 11's
  research net existing.** It also weakens the SC2 overlap check (different input surface) — handled
  by D-12.

- **D-03: Condition the evidence net on the registration-uncertainty level** (Phase 11 D-03's extra
  input). The log-BF pair weakens as registration uncertainty grows. Rationale: Phase 11 D-07 commits
  to posterior width growing with uncertainty; a Bayes factor that stayed confident while the
  posterior widened would be an internal contradiction inside one codebase.

- **D-04: Pre-register thresholds, one documented iteration.** Carried forward from Phase 11 D-04 and
  **not re-asked**: lock every Phase-13 threshold in a Phase-13 consts file with the fresh DEV seed
  **before any run**, and declare the single-iteration allowance in advance.

**Hypothesis boundaries**

- **D-05: Two-factor cut — `ρ_sample` sign × control contrast.** Exclusion is claimed only when the
  sample is genuinely anti-correlated (`ρ_sample < -τ`) **and** stands apart from the control; coloc
  when `ρ_sample > +τ` and above control; random otherwise.
  **Why not the Δρ contrast:** `infer.jl:117-121` defines `Δρ = mean(ρ_sample_draws .- ρ_control_draws)`
  — a sample-vs-control contrast. A sample at `ρ = 0.8` under a control at `ρ = 0.9` has `Δρ < 0`.
  The binary net calls that "null" (harmless); a naive three-way extension on the same quantity would
  publish it as **"mutually exclusive"** while the sample is strongly colocalized. "Less colocalized
  than the control" is not segregation.

- **D-06: τ comes from measured summary resolution.** Probe the smallest ρ difference the fixed
  patch-correlation summary can reliably distinguish and set τ to it, so "random" means
  *indistinguishable from zero at this summary's resolution*. Same measure-don't-guess posture as
  Phase 11 D-06. **τ must be frozen in the Phase-13 consts file before any labelled data is
  generated.**

- **D-07: Stratify training coverage; restore scale with a 3-way prior-odds correction.** Build the
  evidence net's training set with stratified coverage across `ρ ∈ [-0.99, -τ]` rather than sampling
  straight from π(θ), then restore the correct evidence scale at read time via a **3-way
  generalization of `measure_log_prior_odds`** — the mechanism `bf.jl:74`
  (`amortized_log_bf(est_bf, Z_pair, log_prior_odds)`) already uses for the binary case.
  **π(θ) is untouched**, so the CLAUDE.md constraint that π(θ) match the Turing μ/ν/σ/τ ranges holds.
  **The 3-way correction must be DERIVED AND VERIFIED, not assumed to generalize.**
  Rejected: μ-prior truncation (spike 013) — VALIDATED *with cost*, it breaks that same
  prior-consistency constraint and degrades high-|ρ| coverage.

- **D-08: Output = log-BFs against the random reference.** Emit `log BF(coloc : random)` and
  `log BF(exclusion : random)`; the random entry is identically 0. Same semantics as the shipped
  binary `bayes_factor`, and **no model prior over the three hypotheses has to be invented**
  (turning evidence into a chosen hypothesis is Phase 14's job).
  **API note:** the sketched field name `log_bf_simplex` (`src/results.jl:172`) becomes a misnomer —
  it is not a simplex. Rename it when implementing `ThreeHypothesisColocResult`.

**Three-way head architecture**

- **D-09: Shared trunk, two BCE heads.** One Flux conditioner feeding two binary-cross-entropy heads
  (coloc-vs-random, exclusion-vs-random). Satisfies SC1's "one forward pass", preserves the logit/BCE
  semantics that make the D-07 prior-odds correction valid, and the shared conditioner makes the two
  evidences structurally consistent. **Cost accepted:** it lives outside `RatioEstimator`, inherits
  no gate lineage from the shipped binary net, and must be validated on its own terms (D-12/D-13).

- **D-10: Reuse the exact existing trunk.** `RATIO_SUMMARY_WIDTH = 256`, `RATIO_NUM_SUMMARIES = 64`,
  topology `Chain(Dense(in,256,gelu), Dense(256,256,gelu), Dense(256,256,gelu), Dense(256,64))`
  verbatim. Only the head and the input width change (`ratio_input_dim(G) = 5·G² = 320` at G=8,
  **plus 1** for the D-03 uncertainty input). Any difference in the 2-way overlap is then attributable
  to the head alone. Backed by spike 012 (**PARTIAL**): capacity *redistributes* error rather than
  removing it, so capacity is not a reliable lever in this codebase.

**Validation target (SC2 amended)**

- **D-12: Gate on simulation-based 3-way discrimination; report binary-NRE continuity.**
  **Gate:** per-class AUC plus a full confusion matrix over {exclusion, random, coloc} at known
  ground-truth labels — the method spike 014 validated and the GO adopted.
  **Reported, not gated:** correspondence with the **shipped binary NRE** in the overlapping 2-way
  regime at a fixed uncertainty level, as a continuity check.
  **ROADMAP SC2 is amended** — it must NOT be validated against `compute_BayesFactor()`, which named
  limit #3 established is invalid past |logBF| ≈ 6.9, precisely the confident-call tail where
  three-way verdicts matter most. **Document the amendment up front, before results.**

- **D-13: Also gate evidence calibration.** Reliability of the implied class probabilities
  (ECE/MCE via the existing `_bin_calibration` / `CalibrationResult` traffic-light), populating the
  `CalibrationMeta` field already sketched on `ThreeHypothesisColocResult`.
  **Why this does not contradict named limit #3:** that retreat was forced by a broken *reference*
  (the KDE baseline), not by a finding that NRE magnitudes are wrong — nothing was ever measured.
  Classifier calibration is checked directly against ground-truth labels with **no external
  reference**, so the obstacle does not apply. Phase 14's abstention layer needs validated magnitudes.

**Exclusion ground truth (SC3)**

- **D-15: Simulator gate + semi-synthetic graded series + the physical anchor as a named check.**
  **Corpus audit (`corpus/manifest.csv`):** exactly **one** real segregated anchor
  (`physical-primary | negative | segregated | 0.0 | mitochondria`), one positive coloc anchor, and
  30 `simulated-secondary` CBS Red-Green rows whose degrees are colocalization *percentages* — so the
  CBS `0.0` end is **random, not exclusion**. **The corpus has no graded segregation series at all.**
  Therefore: gate on simulator ground truth (`ρ_true < -τ`, exact labels, well-powered); add a
  semi-synthetic graded random→exclusion series as supporting evidence; report the single physical
  anchor as a named qualitative check.
  **Conceptual gap this addresses:** the simulator's "exclusion" is *negative intensity correlation
  across patches*; a biologist's "mutually exclusive" is *disjoint spatial localization*. These
  coincide often but not necessarily, and the semi-synthetic series is the only evidence that probes
  it.

- **D-16: Semi-synthetic series = mask-based disjoint reassignment.** Threshold ch1 into an object
  mask and progressively redistribute ch2 signal into its complement, graded by a mixing parameter
  `α ∈ [0, 1]` from random (α=0, unchanged) to fully disjoint (α=1). Segregation is constructed
  **spatially**, so the series tests the simulator's assumption rather than restating it.
  **Rejected: displacement of ch2 objects** — a global displacement is exactly the mis-registration
  nuisance Phase 11 marginalizes over, so the net could read the series as shift; a null result would
  be uninterpretable and a positive result unattributable.
  **Rejected: intensity anti-correlation** — reproduces the simulator's own mechanism, so it tests
  the net against its own training assumption (circular) and leaves the localization gap untouched.
  The segmentation/threshold rule is a pre-registered choice (D-04 applies).

### Claude's Discretion

- **D-11: Per-head class restriction and joint trunk training.** Each head trains **only on its own
  pair of classes** (coloc head sees coloc + random; exclusion head sees exclusion + random). One-vs-rest
  training would make the logit a coloc-vs-*mixture* ratio rather than `log BF(coloc : random)`,
  silently breaking D-08. The **trunk trains jointly** (both head losses summed) so it sees all three
  classes — that is what makes the two evidences structurally consistent per D-09.
  *Decided by Claude as following necessarily from D-08/D-09; user may override.*

- **D-14: No over-power workaround needed for the calibration gate.** `_bin_calibration` /
  `CalibrationResult` score ECE/MCE against a **traffic-light band**, not a p-value, so the D-13 gate
  cannot be over-powered the way Phase 7's M=2000 SBC point-null was. No TOST/equivalence machinery
  is required here (unlike Phase 11 D-08). *Decided by Claude; user may override.*

- **Sequencing.** Phase 13 planning proceeds now; execution blocks on Phase 11's research net
  existing (see D-02). *Decided by Claude as GSD-internal sequencing.*

### Deferred Ideas (OUT OF SCOPE — do not plan, do not research further)

- Reshipping the 3-way evidence net.
- Extending the corpus with real graded segregation anchors.
- μ-prior truncation to remove the atoms at source (spike 013).
- Full pairwise log-BF matrix including a direct coloc-vs-exclusion contrast.
- Posterior model probabilities over the three hypotheses.
- Training the 3-way net on the shipped grid-8 basis as a second, input-comparable net.
</user_constraints>

---

<phase_requirements>
## Phase Requirements

**No requirement IDs are mapped to this phase.** `.planning/ROADMAP.md:258` records
`**Requirements**: TBD` for Phase 13, and `grep` over `.planning/REQUIREMENTS.md` returns no
Phase-13 requirement block. [VERIFIED: repo grep 2026-07-25]

The phase is therefore governed by its three ROADMAP Success Criteria
(`.planning/ROADMAP.md:260-262`), **one of which D-12 amends**:

| SC | ROADMAP text (verbatim) | Phase-13 status |
|----|---|---|
| SC1 | "A 3-way `RatioEstimator`/evidence network emits a log-BF simplex over {coloc, random, exclusion} in one forward pass" | **Partially superseded.** D-09 rules out `RatioEstimator`; D-08 rules out "simplex". The surviving, bindable content is: *one network, one forward pass, three-way evidence.* |
| SC2 | "It reproduces `compute_BayesFactor()` (`src/bayes.jl:109`) in the overlapping 2-way regime without quadgk/KDE" | **AMENDED by D-12.** Replaced by simulation-based 3-way discrimination (per-class AUC + confusion matrix). Binary-NRE continuity reported, not gated. |
| SC3 | "The exclusion hypothesis is validated on segregated ground-truth inputs" | **Retained, scoped by D-15.** Simulator ground truth gates; α-series is load-bearing supporting evidence; physical anchor is a named qualitative check (see blocker in §J3). |

**Planner action:** because SC1 and SC2 as written are no longer literally achievable/desirable, the
plan must state the amended criteria explicitly and cite the amendment document (§H) — otherwise
`/gsd:verify-work` will score the phase against text the phase deliberately does not satisfy.
</phase_requirements>

---

## Project Constraints (from CLAUDE.md)

| Directive | Source | Phase-13 consequence |
|---|---|---|
| Spike lives entirely in `ProteinCoLoc/spike/` with its own `Project.toml`/`Manifest.toml`; `src/` and both manuscript pipelines remain **provably untouched** during the spike | CLAUDE.md §Constraints | D-01's lane. **Tension with `src/results.jl` — resolved in §G3.** |
| CPU-only is the portable reproducible baseline; GPU is an optional accelerator that must degrade gracefully | CLAUDE.md §Constraints | Train/eval CPU-default; any `use_gpu` must default `false` on the read path. |
| Everything reproducible from a fixed (Random123) seed | CLAUDE.md §Constraints | Philox-per-index; **plus** `Random.seed!` for anything that touches the global RNG (Flux `DataLoader(shuffle=true)`, `sampleposterior`). |
| Windows-tauglich; Julia-native stack | CLAUDE.md §Constraints | No new deps; JLD2 not HDF5. |
| Summary dimension fixed (8×8) during the spike | CLAUDE.md §Constraints | G = 8; `ratio_input_dim(8) = 320`. |
| Priors: simulator π(θ) must stay consistent with the Turing μ/ν/σ/τ ranges | CLAUDE.md §Constraints | **π(θ) must not be edited** — this is why D-07 needs a *correction*, not a prior change (§F). |
| NeuralEstimators API is young (v0.x) — verify signatures against docs/source, not memory | CLAUDE.md §Development Tools | Done: §A is verified against the installed 0.2.1 source. |
| Do not add `NormalizingFlows.jl` / `InvertibleNetworks.jl` / RxInfer / BayesFlow | CLAUDE.md §What NOT to Use | Phase 13 adds **zero** packages (§Standard Stack). |
| GSD workflow enforcement: no direct repo edits outside a GSD workflow | CLAUDE.md §GSD Workflow Enforcement | All Phase-13 edits happen under `/gsd:execute-phase`. |

---

## Summary

Phase 13 is, structurally, a **relabelling plus a head swap** — the CONTEXT's framing is correct and
I verified its premises. The existing evidence net already sees the entire `ρ ∈ [−0.99, +0.99]` range
(`GHAT_RHO_KNOTS`, `ghat.jl:44`) and already lumps exclusion into a single "null" class via
`RATIO_SPLIT_THRESHOLD = 0.0` (`train_ratio.jl:55`). No new physics, no new simulator, no new package.
What is genuinely new is (a) a **three-class label definition** that is *not* the Δρ contrast the
shipped net uses, (b) a **two-head Flux architecture outside `RatioEstimator`**, and (c) a
**per-head evidence-scale correction** whose derivation is materially different from the binary case.

Three findings reshape the plan beyond what CONTEXT already records.

**First — the architecture fork has a cheaper, better-supported route than a hand-rolled loop.**
`RatioEstimator`'s hard-coded loss is real and verified (`RatioEstimator.jl:124`), *but* the generic
fallback one file over is `_loss(estimator, loss) = loss` (`train.jl:654`) and the generic
input/output hook is `_inputoutput(estimator, Z, θ) = (Z, θ)` (`train.jl:657`). A **new subtype of
`NeuralEstimator`** therefore inherits NeuralEstimators' entire `train` loop — CosAnneal LR schedule,
early stopping, best-checkpoint-on-validation-risk, AdamW-Float64 handling — *and* accepts a custom
loss. That is the single highest-leverage finding in this document: it makes D-10's "any difference is
attributable to the head alone" defensible, because the optimizer, schedule and stopping rule stay
byte-for-byte the ones `train_ratio` used.

**Second — the D-07 correction does not generalize the way the binary formula suggests, and the
binary formula itself contains a (harmless) redundancy.** Under NeuralEstimators' shuffled-θ NRE
construction the ratio difference *already is* the log Bayes factor, so `amortized_log_bf`'s
`− log_prior_odds` is an extra ≈ 0.01-nat offset, invisible only because the binary labels are
balanced by construction (measured `P(Δρ>0) = 0.4993` [VERIFIED: 2 M-draw Monte Carlo, 2026-07-25]).
The D-09 heads are **plain BCE classifiers on class labels**, not shuffled-θ NREs — for them the
correction is genuinely required and equals `log(q_pos/q_neg)` **restricted to the head's own two
classes**. Crucially, that scalar is exact **if and only if** stratification changes class
*frequencies* and not the *within-class* θ-distribution. D-07's literal wording ("stratified coverage
across `ρ ∈ [−0.99, −τ]`") is within-class stratification, which a scalar cannot repair. §F derives
this and gives three admissible designs plus a runnable closed-form verification.

**Third — the three classes are not imbalanced, and the exclusion class is the *modal* one.**
`ghat` maps `ρ_true = 0` to `μ = 0.086485`, so `P(ρ_true ≤ 0) = 0.6097` under π
[VERIFIED: exact truncated-Cauchy CDF, 2026-07-25]. Under the D-05 two-factor cut at τ = 0.10 the
class masses are **E 34.8% / R 41.1% / C 24.1%**, with per-head π-implied log-odds of −0.53 (coloc)
and −0.17 (exclusion) — an order of magnitude larger than the binary net's −0.0102. So the
correction is decision-relevant even before any stratification, and "stratify because exclusion is
rare" is a false premise: stratification is about **coverage inside the exclusion range**, which is
dominated by the 4.85% point mass at ρ = −0.99 (vs 1.91% at +0.99, ratio 2.54× — CONTEXT's ≈2.5×
confirmed).

**Primary recommendation:** implement the two-head trunk as a `ThreeWayEvidenceNet <: NeuralEstimator`
with a masked two-term logit-BCE loss, trained through NeuralEstimators' own `train`; take the
scalar per-head correction route (D-07-i, class-frequency stratification with π-shaped within-class
draws) as the primary design and the importance-weighted route (D-07-ii) as the documented fallback;
verify the correction on a closed-form Gaussian toy before a single image is simulated; and measure
τ with a fit-free, unpaired AUC probe on the mean-patch-correlation readout, which is executable
against the Phase-11 *simulator* alone and therefore does not wait on the Phase-11 *net*.

---

## Architectural Responsibility Map

Phase 13 is a single-process research computation with no client/server tiers. The meaningful
"tiers" here are the project's own frozen-surface layers, so the map is stated against those.

| Capability | Primary tier | Secondary tier | Rationale |
|---|---|---|---|
| Image → 128-dim summary | Frozen `src/` summary path (`patch_summary` → `encode_d01`) | — | D-10/CLAUDE.md fix the summary; Phase 13 must not touch it. |
| Summary standardization (`zt`) | **Phase-11 research NPE (frozen)** | — | D-02: the ratio net *inherits* a frozen `zt`, never fits one (`train_ratio.jl:144` pattern). |
| Uncertainty level λ | Phase-11 conditioning input, appended **after** standardization | — | Phase-11 D-03 / 11-RESEARCH D13; λ is not a θ row and is not z-scored. |
| Class labelling (E/R/C) | **Phase-13 spike code**, from simulator ground-truth θ | — | Labels come from `ρ_true` and the control's `ρ_true`, never from a network output. |
| Three-way evidence | **Phase-13 two-head trunk** (`spike/`) | — | D-09; outside `RatioEstimator`. |
| Evidence-scale correction | **Phase-13 read surface** (per-head scalar) | — | D-07; must be applied at read time, exactly as `amortized_log_bf` does today. |
| Discrimination + calibration gate | **Phase-13 validation** (`spike/`), reusing `_bin_calibration` | spike-014 method | D-12/D-13. |
| Result type | **`spike/`** (see §G3 for the `src/` tension) | `src/results.jl` sketch (comment-only) | D-01 keeps `src/` untouched. |
| Physical anchor read | `corpus/` sealed holdout | — | **Blocked** (§J3). |

---

## Standard Stack

**Phase 13 adds zero packages.** Every capability it needs is already resolved in
`spike/Manifest.toml`. This is a deliberate constraint, not an accident: the spike environment's
`NeuralEstimators 0.2.1` pin is fragile under co-resolution (the Phase-1 GLMakie regression, spike
`runtests.jl` gates (d)–(h)), and the resolve-risk gate re-asserts the pin on every run.

### Core (versions read from the resolved manifests, not from memory)

| Library | Version | Purpose | Evidence |
|---|---|---|---|
| NeuralEstimators | **0.2.1** | `train` loop (CosAnneal + early stopping + best-checkpoint) reused by the custom estimator subtype; `logratio` for the shipped-binary continuity check | `spike/Manifest.toml` (`git-tree-sha1 = c780174a2d28714231c0976d5ca4f7c17ac69aa8`); `~/.julia/packages/NeuralEstimators/gFxuZ/Project.toml:3` [VERIFIED] |
| Flux | **0.16.10** | `Chain`/`Dense`/`gelu` trunk, `logitbinarycrossentropy`, `Optimisers.AdamW`, `@layer` | `spike/Manifest.toml`; `Flux/hrg9M/Project.toml:3` [VERIFIED] |
| Distributions | current (resolved) | π(θ) sampling; closed-form Gaussian toy for the §F verification | `spike/Project.toml` [VERIFIED] |
| Random123 | 1.7.1 (per STATE 03-01) | Philox counter-based per-index seeding | `spike/Project.toml` [VERIFIED] |
| JLD2 | current (resolved) | net + label-table + report persistence | `spike/Project.toml` [VERIFIED] |
| StatsBase | current (resolved) | `tiedrank`/`corspearman` for monotonicity, z-scoring helpers | `spike/Project.toml` [VERIFIED] |
| HypothesisTests | current (resolved) | only if a uniformity test is wanted; **not** required by D-13 (traffic-light band, not p-value) | `spike/Project.toml` [VERIFIED] |
| Images | **0.26.2** | `Images.otsu_threshold` for the D-16 ch1 object mask (same call `src/LoadImages.jl:445` and `spike/contract.jl:68` already use) | `spike/Manifest.toml` [VERIFIED] |
| ImageFiltering | **0.7.12** | available for optional mask smoothing/PSF in the α-series | `spike/Manifest.toml` [VERIFIED] |
| CairoMakie | current (resolved) | confusion-matrix / ROC / α-ladder figures | `spike/Project.toml` [VERIFIED] |

### Explicitly NOT needed / NOT to be added

| Package | Why not |
|---|---|
| `KernelDensity`, `QuadGK` | The KDE Bayes-factor baseline. **Neither is in `spike/Project.toml`** [VERIFIED] — the retired baseline is not even reachable from the spike env. D-12 declines to gate on it; do not add them. |
| `ROCAnalysis`, `MLJ`, `Metrics`-style packages | AUC/ROC is hand-rolled in this project (`src/amortized/ood.jl` `roc_auc`); `spike/test/runtests.jl` gate (h) **asserts `ROCAnalysis` and `MLJ` are absent**. Adding one fails the spike suite. |
| `NormalizingFlows.jl`, `InvertibleNetworks.jl` | Asserted absent by the same gate; CLAUDE.md "What NOT to Use". |
| `Turing` | Asserted absent from the spike env; caps NeuralEstimators below the pin. |
| `ImageSegmentation` / `ImageMorphology` as *new direct deps* | Already transitively available via `Images` 0.26.2; add nothing to `spike/Project.toml`. D-16's mask needs only `Images.otsu_threshold`. |

**Installation:** none. If any task proposes `Pkg.add` inside `spike/`, that is a planning error —
`spike/Manifest.toml` is a committed reproducibility artifact (STATE, Phase 01-03).

---

## Package Legitimacy Audit

**Not applicable — Phase 13 installs no external packages.**

Every library above is already a resolved node in the committed `spike/Manifest.toml` and has been
exercised by ten completed phases. No registry lookup, no `slopcheck` run, and no
`checkpoint:human-verify` install gate is required.

| Package | Registry | Disposition |
|---|---|---|
| (none) | — | No new installs in this phase |

**Packages removed due to slopcheck [SLOP] verdict:** none (no candidates).
**Packages flagged as suspicious [SUS]:** none (no candidates).

**Planner note:** the correct integrity control for this phase is not a package audit but the
existing **resolve-risk gate** in `spike/test/runtests.jl` (assertions (d)–(h)), which re-asserts
`NeuralEstimators == v"0.2.1"` by UUID after every dependency-touching change. Any Phase-13 task that
edits `spike/Project.toml` must extend that gate with a Phase-13 clause, mirroring clauses (e)–(h).

---

## A. Verified API Surface (D-09's premise, checked against installed source)

Everything in this section was read out of
`C:\Users\Manuel\.julia\packages\NeuralEstimators\gFxuZ\` — the exact tree pinned by
`spike/Manifest.toml`. Nothing here is from memory.

### A1. The hard-coded-loss landmine is REAL — with exact evidence

```julia
# src/Estimators/RatioEstimator.jl:124
_loss(estimator::RatioEstimator, loss = nothing) = logitbinarycrossentropy
```
```julia
# src/train.jl:654 — the GENERIC fallback, one file over
# For generic estimators, use the user-specified loss function
_loss(estimator, loss) = loss
```
```julia
# src/train.jl:223 (and :352, :538 — one per train method)
loss = _loss(estimator, loss)
```

[VERIFIED: installed NeuralEstimators 0.2.1 source, 2026-07-25]

`train_ratio.jl:38-39`'s comment is therefore accurate and the CONTEXT's finding #2 is confirmed:
passing a loss to `train` on a `RatioEstimator` is **silently discarded**, because the dispatch on
the concrete type wins over the generic fallback. There is no warning and no error.

### A2. `RatioEstimator` is structurally binary-only — and *why* matters

```julia
# src/Estimators/RatioEstimator.jl:92
RatioEstimator(summary_network, num_parameters::Integer; num_summaries::Integer, kwargs...) = ...
# :86 — the inference network always ends in a SINGLE output
inference_network = MLP(num_summaries + num_summaries_θ, 1; ..., output_activation = identity, ...)
```

The binariness is not in `num_parameters` (which is the *θ dimension*, here abused as a 1-dim model
index) — it is in the **fixed 1-unit output head** plus the logit-BCE loss. Even
`RatioEstimator(net, 3; ...)` would give a 1-output network scoring a 3-vector θ, not a 3-class
softmax. D-09's "architectural fork, not a tuning knob" is correct for a second, stronger reason
than the loss alone. [VERIFIED]

### A3. The NRE training construction — and the derivation consequence

```julia
# src/Estimators/RatioEstimator.jl:94-122 (abridged)
θ̃ = getobs(θ, shuffle(1:K))          # independent (marginal) pairs
θ  = hcat(θ, θ̃);  Z = hcat(Z, Z̃)
output = [ones(Float32,1,K)  zeros(Float32,1,K)]   # 1 = joint, 0 = product-of-marginals
```

This is the textbook Hermans-2020 NRE: the Bayes-optimal logit is
`log[ p(Z,θ) / (p(Z)p(θ)) ] = log[ p(Z|θ)/p(Z) ]`. Hence

```
logratio(Z; m=1) − logratio(Z; m=0) = log p(Z|m=1) − log p(Z|m=0) = log BF(1:0)   exactly.
```

**Consequence (report this to the planner):** `amortized_log_bf` (`bf.jl:74-79`) subtracts
`log_prior_odds` *on top of* that difference. Under the shipped net that term is
`−0.0102` (recorded in spike 014's README) because the labels are balanced — measured
`P(Δρ > 0) = 0.4993` over 2×10⁶ prior draws [VERIFIED: Monte Carlo, 2026-07-25] — so the shipped
numbers are unaffected at any reportable precision. But the term is **not** part of the NRE identity,
and D-07 explicitly forbids assuming the binary formula generalizes. §F derives the correct
per-head term from scratch rather than by analogy. *Do not "fix" `bf.jl` in this phase* — D-01
freezes the shipped accessor, and the offset is below the reporting floor.

### A4. The generic `train` hooks make a custom estimator a first-class citizen

```julia
# src/train.jl:657
_inputoutput(estimator, Z, θ) = (Z, θ)
# src/train.jl:660
data = _inputoutput(estimator, Z, _stripnames(_extractθ(θ)))
# src/Parameters.jl:80 — plain arrays pass through untouched
_stripnames(x::AbstractArray) = x
# ext/NeuralEstimatorsFluxExt.jl:39-41
_construct_train_state(estimator::NeuralEstimator, optimiser) =
    FluxTrainState(estimator, optimiser, Optimisers.setup(optimiser, estimator))
# ext/NeuralEstimatorsFluxExt.jl:67 and :80 — the ONLY contract the loss must satisfy
ls = loss(trainstate.model(input), output)
```

[VERIFIED: installed source]

So for `MyNet <: NeuralEstimator` holding Flux layers:
- `train(est, θ_train, θ_val, Z_train, Z_val; ...)` dispatches on `::NeuralEstimator` (`train.jl:101`) — works;
- the loss you pass is honoured (generic `_loss`);
- inputs/outputs pass through unchanged (generic `_inputoutput`);
- the loss is called as `loss(model(Z_batch), label_batch)` — a two-argument function of
  (network output, target array). **Both may be multi-row matrices**, which is exactly what a
  two-head masked loss needs.

**Requirement:** `Optimisers.setup(optimiser, estimator)` must be able to traverse the struct, so it
needs `Flux.@layer MyNet` (or `Functors.@functor`). `Flux.@layer` exists in 0.16.x
(`Flux/*/src/layers/macro.jl:55`) [VERIFIED].

**Side effect to know about:** `train` defaults `savepath = tempdir()` and unconditionally writes
`best_optimizer.bson` / `final_optimizer.bson` / `loss_per_epoch.csv` there
(`train.jl:168`, `ext/NeuralEstimatorsFluxExt.jl:45-59`). `train_ratio.jl` already lives with this.
Pass `savepath = nothing` if the plan wants a side-effect-free run, or an explicit dir to capture the
loss history as an artifact.

---

## B. The Two-Head Architecture (D-09, D-10, D-11)

### B1. Three implementation routes, with a recommendation

| Route | What it is | Pros | Cons | Verdict |
|---|---|---|---|---|
| **R1 — `<: NeuralEstimator` subtype + custom masked loss + NeuralEstimators `train`** | New struct wrapping the D-10 trunk + two `Dense(64,1)` heads; `Flux.@layer`; pass a 2-term masked logit-BCE to `train` | Inherits CosAnneal LR schedule, early stopping, best-val-risk checkpoint, AdamW-Float64 handling, batching — i.e. **the exact recipe `train_ratio` used**, which is what makes D-10's attribution argument hold | Depends on three *unexported internal* hooks (`_loss`, `_inputoutput`, `_construct_train_state`) of a pre-1.0 package | **RECOMMENDED** |
| R2 — hand-rolled `Flux` loop | Explicit `Optimisers.setup` + epoch loop + manual early stopping | No dependence on internals; total control | Must re-implement CosAnneal + patience + best-checkpoint **identically** or D-10's "difference attributable to the head alone" is no longer true; ~120 lines of gate-relevant code with no test lineage | **FALLBACK** (if R1 breaks on a future NeuralEstimators bump) |
| R3 — patch/fork `RatioEstimator` | Vendor a modified copy | — | Explicitly excluded by D-09; also duplicates the shuffled-θ `_inputoutput`, which is the *wrong* construction for class-labelled heads (§F1) | **REJECTED** |

**R1 is not exotic** — it is the documented extension surface: `_loss`, `_inputoutput` and `_risk` are
generic functions with per-type methods (`PosteriorEstimator.jl:85`, `QuantileEstimator.jl:350`,
`RatioEstimator.jl:124`), and the generic fallbacks exist precisely so a user type works. The risk is
version drift, mitigated by the pinned `spike/Manifest.toml` and by a Wave-0 smoke test that asserts
R1 trains for two epochs on a 200-sample toy before any real training starts.

### B2. The concrete architecture (D-10 numbers verified against source)

Every D-10 number checks out:

| D-10 claim | Source | Verified |
|---|---|---|
| `RATIO_NUM_SUMMARIES = 64` | `src/amortized/train_ratio.jl:53` | ✓ |
| `RATIO_SUMMARY_WIDTH = 256` | `src/amortized/train_ratio.jl:54` | ✓ |
| `Chain(Dense(in,256,gelu), Dense(256,256,gelu), Dense(256,256,gelu), Dense(256,64))` | `train_ratio.jl:77-78` | ✓ (exact) |
| `RatioEstimator(net, 1; num_summaries)` | `train_ratio.jl:79` | ✓ |
| `ratio_input_dim(G) = 5·G²` | `src/amortized/summary.jl:127` | ✓ |
| `ratio_input_dim(8) == 320` | `summary.jl:125` docstring + `5·64 = 320` | ✓ |
| `+1` for the D-03 uncertainty input ⇒ **321** | Phase-11 D-03; 11-RESEARCH D13 | ✓ (see §D2 for the *placement* trap) |
| `RATIO_SPLIT_THRESHOLD = 0.0` | `train_ratio.jl:55` | ✓ |
| `RATIO_PAIR_SEED = 0x00000000004A7107` | `train_ratio.jl:60` | ✓ |

**No mismatch found.** The trunk can be lifted verbatim.

### B3. D-11's masked joint loss — the exact construction

Encode the training targets as a **4-row** matrix so that one `loss(ŷ, y)` call satisfies both
D-11 clauses (per-head class restriction + joint trunk):

```
row 1 : y_C  ∈ {0,1}   label for the coloc head   (1 = coloc, 0 = random)
row 2 : w_C  ∈ {0,1}   participation weight       (1 iff the sample is coloc OR random)
row 3 : y_E  ∈ {0,1}   label for the exclusion head (1 = exclusion, 0 = random)
row 4 : w_E  ∈ {0,1}   participation weight       (1 iff the sample is exclusion OR random)
```

Class → (w_C, w_E): coloc → (1,0); exclusion → (0,1); random → (1,1).
**Random appears in BOTH heads' training sets** — that is what makes both logits ratios *against the
same reference*, which is D-08's whole premise. The trunk sees every sample through at least one
active head, so it trains jointly on all three classes (D-09/D-11).

Loss (weighted mean over active (sample, head) cells, so head imbalance does not silently reweight):

```julia
function masked_two_head_bce(ŷ, y)          # ŷ :: 2×n logits, y :: 4×n targets
    lC = Flux.logitbinarycrossentropy.(view(ŷ,1,:), view(y,1,:); agg = identity)
    lE = Flux.logitbinarycrossentropy.(view(ŷ,2,:), view(y,3,:); agg = identity)
    wC = view(y,2,:); wE = view(y,4,:)
    return (sum(wC .* lC) + sum(wE .* lE)) / (sum(wC) + sum(wE) + eps(Float32))
end
```

**Why not a 3-class softmax:** a softmax head's logit differences are
`log[q(c₁|Z)/q(c₂|Z)]`, which is a *posterior*-odds difference over three classes and requires a
model prior over all three to convert to a BF — precisely what D-08 refuses to invent. Two
independent binary heads give two independent 2-model comparisons, each with its own well-defined
correction. This is not a stylistic choice; it is what keeps §F's derivation valid.

**Why not one-vs-rest (the D-11 rationale, restated concretely):** with `w_C ≡ 1` for all samples, the
coloc head's negative class becomes the *mixture* `q_R·R + q_E·E`, so the corrected logit would be
`log[ p(Z|C) / (q̃_R p(Z|R) + q̃_E p(Z|E)) ]` — a mixture Bayes factor whose value depends on the
training mixture weights and which is **not** `log BF(coloc : random)`. It would silently violate D-08
while looking identical in the code.

### B4. Head shape and the "one forward pass" (SC1) claim

```julia
struct ThreeWayEvidenceNet <: NeuralEstimator
    trunk :: Any        # Chain(Dense(in,256,gelu), Dense(256,256,gelu), Dense(256,256,gelu), Dense(256,64))
    head  :: Any        # Dense(64, 2)   — two logits, one per hypothesis-vs-random comparison
end
Flux.@layer ThreeWayEvidenceNet
(m::ThreeWayEvidenceNet)(Z) = m.head(m.trunk(Z))     # 2×n logits, ONE pass
```

Using a single `Dense(64,2)` rather than two `Dense(64,1)` layers is mathematically identical
(the rows do not interact) and keeps the "one forward pass" claim literal and obvious. Either is
fine; the plan should pick one and state it.

---

## C. What the Shipped Path Actually Does (the surfaces Phase 13 mirrors)

| Surface | Location | Signature / value | Phase-13 use |
|---|---|---|---|
| `build_ratio_estimator` | `train_ratio.jl:72-80` | `(input_dim; num_summaries=64, width=256) -> RatioEstimator` | Trunk topology copied verbatim; the `RatioEstimator(...)` wrap is dropped |
| `measure_log_prior_odds` | `train_ratio.jl:91-96` | `(model_index) -> log(p/(1-p))`, `p = mean(model_index)`; **errors on degenerate balance** | Generalized per-head (§F3) |
| `assemble_ratio_pairs` | `train_ratio.jl:112-134` | pairs two distinct pool columns; label `Float32((ρ[i1]-ρ[i2]) > 0.0)` | Label rule **replaced** by the D-05 three-way cut; the pairing mechanics are reused |
| `train_ratio` | `train_ratio.jl:155-188` | `n=48_000, epochs=300, batchsize=128, lr=2.5e-4, weight_decay=1e-4, val_frac=0.15, stopping_epochs=40`, `AdamW(lr,(0.9,0.999),wd)` **Float64** | Recipe copied verbatim (D-10 attribution) |
| `amortized_log_bf` | `bf.jl:74-79` | `grid = [0 1]`, returns `lr[1,2]-lr[1,1]-log_prior_odds`; `use_gpu=false` | The read pattern to mirror; the shipped call itself is **frozen** |
| `pair_encode` | `bf.jl:91-100` | `vcat(Zs, Zc, Zs[1:G²]-Zc[1:G²])`; **asserts even length**, derives `nc = len÷2` | Reused unchanged — see §D2 for the λ trap it accidentally protects against |
| `rho_draws` / `delta_rho` | `infer.jl:101-121` | `Δρ = mean(ρ_s .- ρ_c)`; `reconstruct(θzt, ·)` applied **before** reading row 1 | D-05's counter-example; **not** the Phase-13 label |
| `api.jl` composition | `api.jl:140-141` | `amortized_log_bf(b.ratio.estimator, pair_encode(Zs,Zc), b.ratio.log_prior_odds)` | The eventual productionization shape. **Not touched in Phase 13.** |
| `_bin_calibration` / `CalibrationResult` | `spike/validation/sbc.jl:65-130` | see §I | D-13's gate |
| `CalibrationMeta` | `src/results.jl:117-126` | 8 fields incl. `grid::Int`, `gate::NamedTuple` | D-13 target shape |

---

## D. The Phase-11 Dependency (D-02, D-03) — What Must Exist Before Execution

### D1. The output contract Phase 11 must deliver

Phase 13 execution needs **all six** of these from Phase 11. Anything missing is a hard block.

| # | Artifact | Why Phase 13 needs it | Where Phase 11 produces it |
|---|---|---|---|
| 1 | **Trained Phase-11 research NPE**, persisted (JLD2, CPU-resident `Flux.state` + arch metadata) | D-02: the ratio net trains on pairs assembled from an already-trained NPE's basis | 11-RESEARCH §D12 (`spike/npe/train_npe.jl` extension) |
| 2 | **Frozen `zt`** (the 128-row summary `ZScoreTransform`, mask rows bypassed) | D-02: standardization must be *inherited*, never re-fit (`train_ratio.jl:144-146`; Pitfall 5) | carried in the persisted NPE handle, same as `npe_8.jld2` |
| 3 | **Frozen `θzt`** (the 8-row θ transform) | needed only to read `ρ_true` back out of posteriors for the D-12 continuity check and for label sanity checks | ditto |
| 4 | **`encode_lambda(λ)` (or its exact equivalent)** and the λ **range** `[LAMBDA_MIN, LAMBDA_MAX]` | D-03: Phase 13 must append λ in *exactly* the encoding the NPE was trained with, or the two nets disagree about what row 129 means | 11-RESEARCH §D13 |
| 5 | **The frozen Phase-11 `imsize_set`/weights** and the widened `SHIFT_PRIOR` + `ε` prior | Phase-13 label generation must simulate under the **same joint** the net was trained on (F5 covariate shift; spike-014 §"Image sizes" made this an explicit verified precondition) | Phase-11 consts file + `training_imsize_provenance` |
| 6 | **The θ arity and field order** (8 fields incl. `chromatic_eps`) | Phase 13 reads `ρ_true` positionally out of θ; an arity change silently shifts every row | `spike/simulator/prior.jl` post-Phase-11 |

### D2. How to express the precondition — and the λ placement trap

**The precondition check the planner should write** (a Wave-0 task that runs in seconds and fails
loudly, mirroring how `run_gate.jl` returns `:not_trained` rather than crashing):

```julia
# spike/p13/preconditions.jl — RUN FIRST; every later task depends on this being green.
const P11_NET = joinpath(@__DIR__, "..", "npe", "trained_npe_p11.jld2")   # exact name from Phase 11
isfile(P11_NET) || error("""
  Phase 13 is BLOCKED on Phase 11 (13-CONTEXT D-02).
  Expected the Phase-11 registration-aware research NPE at:
      $P11_NET
  Phase 13 must NOT fall back to the shipped grid-8 basis — that is a different input
  surface and would invalidate D-02, D-03 and the D-12 continuity framing.""")
h = load_npe(P11_NET)
@assert hasproperty(h, :zt) && hasproperty(h, :θzt)
@assert h.lambda_range isa Tuple{<:Real,<:Real}           # D-03 conditioning input exists
@assert size(h.θzt.mean, 1) == 8                          # ε joined θ (Phase-11 D-09)
@assert h.training_imsize_provenance.recorded == true      # F5: the joint is known
```

**The λ placement trap (high severity, easy to get wrong):**

Phase 11's network input is **129 rows** = `vcat(standardize_summary(128-row summary, zt), encode_lambda(λ))`
— λ appended *after* standardization (11-RESEARCH Pitfall 2). For Phase 13's pair encoding:

```
CORRECT:   Z_pair = vcat( pair_encode(Zs_128, Zc_128), λ )          # 320 + 1 = 321
WRONG:     Z_pair = pair_encode(Zs_129, Zc_129)                     # throws (see below)
ALSO WRONG: Z_pair = vcat(pair_encode(Zs_128,Zc_128), λ, λ)         # λ counted twice
```

The wrong form **fails loudly**, which is a gift: `pair_encode` asserts `iseven(length(Zs))`
(`bf.jl:95-96`) and 129 is odd, so it throws
`ArgumentError: pair_encode: summary length 129 must be 2·G² (even)`. The plan should still state
the correct form explicitly, and a unit test should assert `size(Z_pair,1) == 5*G^2 + 1`.

**Rationale for one λ, not two:** sample and control are two channels of the *same* acquisition under
the *same* stated registration uncertainty. Two λ values would imply the user holds different beliefs
about the two stacks, which the API does not express. State this as a pre-registered modelling choice.

### D3. What is genuinely unknowable until Phase 11 lands

| Unknown | Impact on Phase 13 | Mitigation available now |
|---|---|---|
| Exact `encode_lambda` form (raw λ? `log λ`? z-scored?) | Determines the 321st row's scale; a badly-scaled input can be ignored by the net | Plan reads it from Phase 11 at execution; **do not hard-code**. Unit-assert the value range. |
| Persisted-handle field names (`trained_npe_p11.jld2`? `h.zt`?) | The precondition script's literals | Write the precondition as a **single** task whose only job is to bind those names; every later task imports from it |
| Whether λ is a scalar or a small vector (11-RESEARCH flags a possible second ε-conditioning input in its Open Questions) | Input width 321 vs 322 | Derive width as `ratio_input_dim(G) + n_cond` where `n_cond = length(encode_lambda(λ))`. Never write 321 as a literal. |
| Realized Phase-11 `zt` statistics | The τ probe's *summary-level* numbers are `zt`-independent (raw correlations), so τ is **not** blocked; the *label generation* is | See §E4 — the τ probe runs on the simulator alone |
| Phase-11 net quality (11-RESEARCH H24 lists a "SC2 fails flat" branch) | If Phase 11 lands with a net whose λ input is ignored, D-03 conditioning in Phase 13 is vacuous too | Phase 13 must **measure and report** the λ-response of its own log-BF pair (it is D-03's stated claim). A flat response is an honest finding, not a bug to hide. |

---

## E. Hypothesis Boundaries (D-05, D-06) and the τ Probe

### E1. The D-05 cut, operationalized

D-05 fixes the *level* condition (`ρ_sample < −τ`, `ρ_sample > +τ`) but leaves "stands apart from the
control" to implementation. Two candidates:

| Variant | Definition | Property |
|---|---|---|
| **(a) τ-contrast** | C ⟺ `ρ_s > τ` **and** `ρ_s − ρ_c > τ`; E ⟺ `ρ_s < −τ` **and** `ρ_s − ρ_c < −τ`; R otherwise | One τ governs both factors; classes are symmetric under sign flip; "random" absorbs everything ambiguous |
| (b) sign-contrast | C ⟺ `ρ_s > τ` **and** `ρ_s > ρ_c`; E ⟺ `ρ_s < −τ` **and** `ρ_s < ρ_c`; R otherwise | Contrast factor has no dead zone, so a 0.001 contrast counts as "stands apart" |

**Recommend (a).** It is the only variant in which *both* factors are stated at the summary's measured
resolution — which is D-06's entire point. Under (b) the contrast factor is a strict-inequality
coin-flip near zero and adds label noise exactly in the region where the heads must be calibrated.
[ASSUMED — this is a design recommendation, not a measured result; it is a pre-registration choice
the user may override.]

**Measured class masses under π** (2×10⁶ i.i.d. draw pairs through `ghat(Truncated(Cauchy(0,0.3),-1,1))`,
`Random.seed!(20260725)`) [VERIFIED: computed 2026-07-25]:

| τ | variant (a) E / R / C | variant (b) E / R / C |
|---|---|---|
| 0.05 | 0.3850 / 0.3413 / 0.2738 | 0.4040 / 0.3143 / 0.2817 |
| 0.10 | 0.3482 / 0.4108 / 0.2410 | 0.3828 / 0.3634 / 0.2538 |
| 0.15 | 0.3119 / 0.4713 / 0.2168 | 0.3583 / 0.4078 / 0.2338 |

Implied π-level per-head log-odds under (a): at τ = 0.10, `log(q_C/q_R) = −0.5333`,
`log(q_E/q_R) = −0.1654`. Compare the shipped binary net's `−0.0102`. **The correction is 15–50×
larger than in the binary case and cannot be waved away as ≈0.**

### E2. The prior-atom asymmetry (CONTEXT finding #3) — confirmed exactly

Exact truncated-Cauchy CDF evaluation at the `ghat` clamp points [VERIFIED: 2026-07-25]:

| Quantity | Value |
|---|---|
| `P(μ* ≤ GHAT_MU_MIN = −0.67976)` → atom at `ρ = −0.99` | **0.048527** |
| `P(μ* ≥ GHAT_MU_MAX = 0.847149)` → atom at `ρ = +0.99` | **0.019108** |
| ratio (negative : positive) | **2.5396×** |
| `μ` such that `ghat(μ) = 0` | **0.086485** |
| `P(ρ_true ≤ 0)` | **0.60969** |

Two consequences the plan must absorb:

1. **`ρ_true = 0` is `μ = 0.0865`, not `μ = 0`.** The summary reads a *positive baseline* correlation at
   zero true colocalization (shared structure + background). Anyone reasoning about "the ρ = 0 point"
   in μ-space will be off by ~0.09 μ-units, which is larger than several of the `ghat` knot spacings.
2. **~9% of all exclusion-class draws sit exactly on the ρ = −0.99 rail** (0.0485 / 0.348 at τ=0.10).
   The exclusion head's training set therefore contains a large discrete spike at one θ value. That is
   both a within-class coverage problem (D-07's real motivation) and a hazard for any
   importance-weighting scheme (§F4 — the density is undefined at an atom).

### E3. What "reliably distinguish" should mean — proposed statistic

D-06 asks for *the smallest ρ difference the fixed patch-correlation summary can reliably distinguish*.
Recommended operationalization, chosen to be **fit-free** (hence honestly pre-registerable) and cheap:

> **Statistic.** `m̄(s) = mean of the PRESENT continuous summary rows` — i.e. the mean per-patch Pearson
> correlation over patches that survived the ≥15 valid-pixel floor. Rows `1:G²` only; the mask rows
> `G²+1:2G²` are excluded (they are the rows `_summary_row_partition` itself excludes from
> standardization, `summary.jl:89-98`).
>
> **Resolution measure.** `A(δ) = AUC( m̄ | ρ = 0  vs  m̄ | ρ = −δ )`, estimated from `R` independent
> draws per arm with **nuisances drawn from their full priors and `imsize` from the Phase-11 F5
> mixture**.
>
> **τ = the smallest δ on a pre-registered grid with `A(δ) ≥ AUC_τ`.**

**Why `m̄` and not something learned.** `m̄` *is* the quantity the frozen `ghat` calibration was built
on — `ghat.jl:34` records the sweep as `E[μ | ρ_true]` where μ is exactly the induced mean per-patch
correlation. So `m̄` is the canonical scalar readout of the summary with respect to ρ, it requires no
fitting (no snooping surface), and it is directly interpretable. A learned readout (LDA on 64 rows, or
the NPE's own `ρ̂`) would be *more* powerful and would therefore report a *smaller* τ — but it would
make τ a property of a trained model rather than of the summary, contradicting D-06's wording. Report
the NPE-based number as a secondary column if wanted; **gate on `m̄`**.

**Why UNPAIRED.** Phase 11's D-06 probe is deliberately paired (same latent field, only geometry
changes) because it measures *displacement*. τ measures *discriminability for a single unseen image*,
where the nuisances and the latent field are unknown. A paired design would report the resolution of a
counterfactual the deployed tool never has, and would understate τ by a large factor. **This is the
single most important difference between the Phase-11 probe and the Phase-13 probe and must be stated
in the consts file.**

**Recommended pre-registered probe specification (lock all of this in the consts file BEFORE running):**

| Knob | Recommended value | Rationale |
|---|---|---|
| `P13_TAU_DELTA_GRID` | `(0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20)` | brackets the plausible answer given the atoms and knot spacing; the τ=0.05–0.15 band is where class masses move fastest (§E1) |
| `P13_TAU_AUC` | `0.90` | the "reliably distinguish" bar. Pre-registered as a *number*, so "reliably" is not argued after the fact. Report the whole `A(δ)` curve regardless. |
| `P13_TAU_R` | `≥ 400` per arm | AUC standard error ≈ `0.03` at A≈0.9, R=400 → the 0.90 bar is resolved to ±0.03 |
| direction | measure **both** `A(0 vs −δ)` and `A(0 vs +δ)`; take `τ = max` | the atoms are asymmetric (§E2), so resolution may be too. A single symmetric τ is what D-05 assumes; taking the max is the conservative choice |
| `imsize` | the Phase-11 F5 mixture, **not** 256² | F5 covariate shift is a named limit; a 256²-only τ would not describe the deployment regime |
| λ | measure `τ(λ)` at every rung of the Phase-11 ladder; **freeze the scalar τ at a pre-registered reference λ** | resolution degrades with registration uncertainty; a single τ needs a stated λ. Recommend the ladder's *widest* rung (most conservative τ) or the shipped-equivalent rung — but the choice must be in the consts file before the run |
| statistic reporting | mean **and** median AUC over bootstrap resamples | one degenerate all-`missing` grid can move `m̄` by O(1) (the `encode_d01` missing→0 path) |
| abort criterion | if `A(δ) < AUC_τ` for **every** δ in the grid, τ is not measurable — **stop and report**, do not extend the grid | this is the one threshold that must be immune to its own result, or "below resolution" becomes unfalsifiable (Phase-11 §C10(c) precedent) |

### E4. **The τ probe does NOT block on the Phase-11 net** — sequencing win

The probe needs: `sample_prior`, `simulate_pair`, `patch_summary`, `encode_d01` — all simulator-side.
It needs **no** trained network, **no** `zt`, **no** λ conditioning input (λ enters only as the
*simulator's* shift half-width, which is a prior parameter, not a network input).

So the probe blocks on **Phase 11's simulator surgery** (D-09 `ε`, D-10 composed warp, D-02 widened
`SHIFT_PRIOR`) — which lands early in Phase 11 — **not** on Phase 11's training run. The planner should
sequence the τ probe as an early wave that can run the moment Phase 11's simulator commits merge, and
keep the *labelled data generation* behind the full Phase-11 net precondition.

**Caveat that must be stated:** if the probe is run against the *pre*-Phase-11 simulator
(`SHIFT_PRIOR = Uniform(-1,1)`, no `ε`), the measured τ will be **optimistic** — the deployment
distribution is wider. Running it early against the old simulator and freezing that τ would be a
pre-registration error. Either wait for the simulator merge, or freeze τ from a probe run against the
post-merge simulator. State which, in the consts file.

---

## F. D-07 Derived: The Three-Way Evidence-Scale Correction

This is the section D-07 demands be *derived and verified, not assumed*. It is derived from scratch.

### F1. The two constructions are different, and only one needs a correction

**Construction 1 — NeuralEstimators NRE (what the shipped binary net does).** Classifier separates
joint `(Z, m)` from product-of-marginals `(Z, m̃)` with `m̃` a shuffle of `m`, equal counts
(`RatioEstimator.jl:94-122`). Bayes-optimal logit:

```
s_NRE(Z, m) = log[ p(Z,m) / (p(Z) p(m)) ] = log[ p(Z|m) / p(Z) ]
⇒ s_NRE(Z,1) − s_NRE(Z,0) = log p(Z|1) − log p(Z|0) = log BF(1:0)      — no correction needed.
```

**Construction 2 — plain BCE classifier on class labels (what the D-09 heads do).** Head `h`
separates class A (positive) from class B (negative) with *training* class frequencies `q_A`, `q_B`
(`q_A + q_B = 1` within the head's restricted set, per D-11). Bayes-optimal logit:

```
s_h(Z) = log[ q(A|Z) / q(B|Z) ]
       = log[ q_A · p_q(Z|A) ] − log[ q_B · p_q(Z|B) ]
       = log[ p_q(Z|A) / p_q(Z|B) ] + log(q_A / q_B)
```

**⇒ the D-07 correction, derived:**

```
log BF_q(A : B) = s_h(Z) − log(q_A / q_B)
```

with `q_A`, `q_B` **measured on the head's own restricted training subset** — not on the full
three-class set, and not on π.

### F2. When the scalar is EXACT, and when it is not

`p_q(Z|A) = ∫ p(Z|θ) q(θ|A) dθ`, whereas the reportable quantity is
`p_π(Z|A) = ∫ p(Z|θ) π(θ|A) dθ`. Therefore:

> **Theorem (informal but exact).** The scalar correction `− log(q_A/q_B)` recovers the π-scale Bayes
> factor **if and only if** `q(θ|A) = π(θ|A)` and `q(θ|B) = π(θ|B)` — i.e. stratification changed only
> the class *frequencies*, never the within-class θ shape. No scalar can repair a within-class
> distributional change, because the change enters *inside* the integral and depends on `Z`.

D-07's literal instruction — "stratified coverage across `ρ ∈ [−0.99, −τ]`" — is a **within-class**
reshaping. So the naive reading of D-07 (stratify the ρ range, then apply a scalar) is **provably
insufficient**. This is exactly the failure D-07 anticipated when it said "must be DERIVED AND
VERIFIED, not assumed to generalize."

### F3. Three admissible designs (pick one, pre-register it)

| # | Design | Correction needed | π-consistency | Cost / risk |
|---|---|---|---|---|
| **D-07-i** | **Class-frequency stratification only.** Draw θ ~ π; oversample/undersample *whole classes* to target frequencies (accept/reject or replicate at class granularity). Within-class shape stays π\|class. | scalar `log(q_A/q_B)` — **exact** | ✓ exact | Cheapest, no new machinery. Does not increase coverage *inside* `[−0.99,−τ]` beyond what π gives — but §E1 shows π already gives 35% exclusion mass, so the coverage need is about the −0.99 atom, not the range. **RECOMMENDED as primary.** |
| **D-07-ii** | **Within-range ρ-stratification + importance-weighted BCE.** Draw ρ from a stratified proposal `g(ρ)` over the exclusion range; attach weight `w = π(ρ)/g(ρ)`; train with weighted BCE. | scalar, applied to the *weighted* frequencies | ✓ exact in the limit | Needs (1) weighted BCE — trivial, the masked loss already has weights; (2) the π-density of ρ (closed form, §F5); (3) **atom handling** — `π(ρ)` has a point mass at −0.99 where a density ratio is undefined. Must treat `{ρ = −0.99}` as its own stratum with a discrete weight. Weight variance inflates gradient noise. **DOCUMENTED FALLBACK.** |
| D-07-iii | **Redefine the hypotheses** as "ρ uniform over the exclusion range" etc. | scalar | ✗ — the reported BF is w.r.t. a different prior than π | Legitimate *if declared*, but it changes what "exclusion" means and breaks comparability with the shipped binary BF and with the Turing model. **Only with an explicit named limit.** |

**Recommendation:** plan D-07-i as the primary path and D-07-ii as a pre-declared single-iteration
fallback (which fits D-04's one-iteration allowance). If the D-12 confusion matrix shows the exclusion
class is under-resolved *specifically at high |ρ|*, that is the trigger to spend the iteration on
D-07-ii. Declare the trigger in advance.

### F4. The 3-way generalization of `measure_log_prior_odds`, concretely

```julia
"""
    measure_head_log_odds(class_labels; positive, negative) -> Float64

3-way generalization of `measure_log_prior_odds` (src/amortized/train_ratio.jl:91-96).

The binary version measures ONE scalar over ALL training labels. With two heads under D-11's
per-head class restriction, that is the WRONG quantity: each head's Bayes-optimal logit carries
the log-odds of ITS OWN restricted pair, so there are TWO scalars and each is measured only over
the samples that head actually trained on.

    n_pos = count(==(positive), class_labels)
    n_neg = count(==(negative), class_labels)
    return log(n_pos / n_neg)

DERIVATION (13-RESEARCH §F1): for a BCE head separating A from B with training frequencies
q_A, q_B, the Bayes-optimal logit is log[p_q(Z|A)/p_q(Z|B)] + log(q_A/q_B), so the log Bayes
factor is the logit MINUS log(q_A/q_B).

NOT the same as the binary NRE case: NeuralEstimators' RatioEstimator uses a shuffled-θ
construction whose ratio difference is ALREADY the log BF (RatioEstimator.jl:94-122), so its
`- log_prior_odds` term is a (numerically negligible, ~0.01 nat) redundancy. Do not copy it.
"""
function measure_head_log_odds(class_labels::AbstractVector; positive, negative)
    n_pos = count(==(positive), class_labels)
    n_neg = count(==(negative), class_labels)
    (n_pos > 0 && n_neg > 0) ||
        error("measure_head_log_odds: degenerate head balance n_pos=$n_pos n_neg=$n_neg")
    return log(n_pos / n_neg)
end
```

Read surface (mirrors `amortized_log_bf`'s shape and CPU default):

```julia
function three_way_log_bf(net, Z_pair, head_log_odds; use_gpu::Bool = false)
    Zc = Z_pair isa AbstractVector ? reshape(Z_pair, :, 1) : Z_pair
    s  = net(Zc)                                   # 2×1 logits, ONE forward pass
    return (coloc     = Float64(s[1,1]) - head_log_odds.coloc,
            exclusion = Float64(s[2,1]) - head_log_odds.exclusion,
            random    = 0.0)                       # D-08: the reference is identically 0
end
```

### F5. The runnable numerical verification (this is what "VERIFIED" must mean)

**Design a closed-form toy where the true three-way log-BFs are analytic**, then check the corrected
heads recover them. No quadrature, no new dependency — conjugate Gaussians only.

```
θ ~ N(0, s²)                       # the "ρ" analogue
Z | θ ~ N(θ, σ²)                   # the "summary" analogue (1-D or k-D i.i.d.)
classes:  E = {θ < −τ},  R = {|θ| ≤ τ},  C = {θ > τ}
```

The class-conditional evidence is closed form. With the marginal `Z ~ N(0, s²+σ²)` and the
posterior `θ|Z ~ N(m, v)`, `m = Z s²/(s²+σ²)`, `v = s²σ²/(s²+σ²)`:

```
p(Z | class) = N(Z; 0, s²+σ²) · [ Φ((b−m)/√v) − Φ((a−m)/√v) ] / [ Φ(b/s) − Φ(a/s) ]
log BF(C:R) = log p(Z|C) − log p(Z|R)          # exact, for every Z
```

**Verification protocol (a Wave-0/Wave-1 task, minutes on CPU):**

1. Sample a training set with **deliberately skewed** class frequencies (e.g. 60/25/15) so the
   correction is large and a missing correction is impossible to miss.
2. Train the D-09 two-head trunk (same `ThreeWayEvidenceNet`, smaller widths are fine) with the
   masked loss.
3. Measure `head_log_odds` with `measure_head_log_odds`.
4. On a held-out `Z` grid, compare `s_h(Z) − log(q_A/q_B)` against the analytic `log BF`.
5. **Pre-registered pass bars** (lock in the consts file): `cor ≥ 0.99` and
   `max |Δ log BF| ≤ 0.25` over the central 90% of the `Z` grid, plus a **negative control**:
   the *uncorrected* logit must FAIL the same bars by roughly `|log(q_A/q_B)|`. Without the negative
   control the test cannot distinguish "correction works" from "correction is unnecessary here."
6. Repeat step 1–5 with a **within-class reshaped** proposal (D-07-ii's scenario) and confirm the
   scalar-only correction **fails** — this is the empirical demonstration of §F2's theorem and it is
   what makes the D-07-i choice evidence-backed rather than asserted.

**The π-density of ρ (needed only for D-07-ii), closed form.** `ρ = ghat(μ)` with `ghat` piecewise
linear and strictly increasing on `[GHAT_MU_MIN, GHAT_MU_MAX]`, so on segment `i`

```
π_ρ(ρ) = f_μ(ghat⁻¹(ρ)) · (Δμ_i / Δρ_i),      f_μ = pdf of Truncated(Cauchy(0,0.3), −1, 1)
```
plus the two atoms `P(ρ = −0.99) = 0.048527`, `P(ρ = +0.99) = 0.019108` [VERIFIED, §E2].
`ghat⁻¹` is a 4-line inverse of the same knot table (`searchsortedlast` on `GHAT_RHO_KNOTS`).

---

## G. Output Semantics, the Rename, and the `src/` Boundary (D-08)

### G1. The rename

`src/results.jl:172` sketches `log_bf_simplex :: NTuple{3, Float64}  # {null, coloc, anti-coloc}`
with accessor `bayes_factor_simplex` (`:178`). D-08 is right that this is a misnomer: the three
values are `(log BF(C:R), 0, log BF(E:R))` — two free numbers plus a structural zero, not a point on
a 2-simplex. A reader who sees "simplex" will assume the entries sum to one, or that they are
probabilities.

**Recommended replacement:**

```julia
struct ThreeHypothesisColocResult <: AbstractColocResult
    grid             :: Int
    posterior        :: Matrix{Float64}
    log_bf_vs_random :: NamedTuple{(:coloc, :random, :exclusion), Tuple{Float64,Float64,Float64}}
    ood              :: OODVerdict
    calibration      :: CalibrationMeta
    meta             :: NamedTuple
end
log_bf_vs_random(r::ThreeHypothesisColocResult) = r.log_bf_vs_random
```

Why a `NamedTuple` and not an `NTuple{3}`: it makes the reference class self-documenting at every call
site and removes the positional-order hazard the sketch already has (the sketch's comment orders the
tuple `{null, coloc, anti-coloc}` — i.e. the reference class is *first* — which no reader would guess
from a name containing "simplex"). Cost: a slightly heavier type signature. Worth it.

**Phase 7 D-02 compliance** (`07-CONTEXT.md`, restated at `src/results.jl:19-33`): new variants must
be **new `AbstractColocResult` subtypes**, never bolted-on fields. `ThreeHypothesisColocResult`
satisfies this. It must also implement the four shared accessors or it will hit
`_iface_error` (`results.jl:48`):

| Accessor | Recommended implementation | Note |
|---|---|---|
| `posterior_draws` | `r.posterior` | trivial |
| `is_ood` | `r.ood.flag` | trivial |
| `bayes_factor` | `r.log_bf_vs_random.coloc` | **Must be documented**: the shared single-number accessor returns the *coloc-vs-random* entry, so a caller written against the binary result gets the semantically closest number rather than an error. State this in the docstring or a reader will assume it is a three-way summary. |
| `delta_rho` | `r.posterior`-derived Δρ draws **if** a control posterior is carried; otherwise `_iface_error` is the honest answer | Do not fabricate a Δρ. If the Phase-13 result does not carry the control draws, letting the interface error is correct behaviour. |

### G2. `bayes_factor_simplex` should not survive either

The sketch's accessor name inherits the same misnomer. Recommended: `log_bf_vs_random`. If backwards
compatibility is ever wanted, a deprecated alias is cheap — but nothing consumes the sketch today
(it is a comment), so a clean name costs nothing now and everything later.

### G3. **The `src/` boundary tension — reported, not guessed**

**The facts:**

1. CLAUDE.md §Constraints: *"the main package, `src/`, and both manuscript pipelines must remain
   **provably untouched** during the spike."*
2. D-01: *"All Phase-13 work lives in `spike/`. The shipped artifact, the shipped `bayes_factor`
   accessor and the Phase-7 GO stay untouched."* — note this enumerates the **artifact**, the
   **accessor** and the **GO**, not the whole directory.
3. Phase 7 was declared "the ONLY phase that edits `src/`" (STATE, Roadmap decisions) — but that
   posture has already been relaxed: **Phase 11 D-11 explicitly mirrors `ε` into
   `src/amortized/simulator.jl`**, with D-12 recording the provenance cost. So `src/` edits are not
   categorically forbidden post-Phase-7.
4. `src/results.jl` is **not** in `DATAGEN_HASH_SRC_FILES` (`src/amortized/datagen.jl:208-213`
   lists `colocalization.jl`, `LoadImages.jl`, `amortized/summary.jl`, `amortized/simulator.jl`)
   [VERIFIED]. So editing it has **zero** cache-provenance cost.
5. The Phase-13 block at `results.jl:161-192` is **pure comments** — editing it cannot change runtime
   behaviour, cannot affect the artifact (pinned by `git-tree-sha1` in `Artifacts.toml`), and cannot
   affect the Phase-7 GO.

**The tension is therefore between (1)'s literal wording and (2)'s enumerated scope.** They are not
reconcilable by inference; the planner should surface it rather than pick silently.

**Recommended resolution — two separable decisions:**

| Decision | Recommendation | Reason |
|---|---|---|
| Where does the **executable** `ThreeHypothesisColocResult` live? | **`spike/`** — e.g. `spike/p13/result.jl`, defined as `<: AbstractColocResult` after a read-only `include(joinpath(@__DIR__,"..","..","src","results.jl"))` (the same read-only reach `corpus/load.jl:44-45` and `spike/contract.jl` already use for `src/LoadImages.jl` / `src/colocalization.jl`) | Honours D-01 and CLAUDE.md fully, and still structurally satisfies Phase 7 D-02 (it *is* a new subtype of the real supertype). Zero `src/` bytes changed. |
| Does the **comment-only** rename land in `src/results.jl:172,178` now? | **Recommend NO for Phase 13; carry it as a one-line item for the productionization phase** — and record the misnomer in the Phase-13 report so the sketch cannot be quoted as authoritative in the meantime | It is a `src/` edit with no scientific benefit inside this phase, and D-01's plainest reading forbids it. If the user *wants* it now, it should be a separate, explicitly-labelled docs-only commit (`docs(results): rename Phase-13 sketch field ...`) that touches nothing else — cheap, but it is a user decision, not a planner inference. |

**Planner action:** put this in the plan as an explicit `checkpoint:human-verify` or a stated
assumption, not as a silent choice.

---

## H. Documenting the SC2 Amendment (D-12) — the Established Precedent

D-12 requires the amendment be documented **before results**. This project has a precise,
twice-used precedent, and Phase 13 should follow it verbatim.

### H1. Where amendments live in this project

| Artifact | What it is | Phase-13 analogue |
|---|---|---|
| `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md` | The canonical form: a standalone, dated, **`Status: FROZEN SPECIFICATION. Nothing has been run against it.`** document, with a §0 "Provenance disclosure — read this first", numbered defects (A1…), and per-item outcome-independence arguments | **`.planning/phases/13-.../13-SC2-AMENDMENT.md`** |
| `test/gate/gate_consts_8_v2.jl` | The *executable* companion: a **new, separately named** pre-registration file that supplements (never replaces) the original, with recomputed forbidden seeds and `@assert` self-checks | **`spike/p13/consts.jl`** (D-04) |
| `docs/amortized.md` §Named limits | Where amendments graduate into permanent published honesty (limits #1 and #3 are literally about the twice-amended gate and the retired KDE reference) | A Phase-13 entry is owed **only if** Phase 13 productionizes; while it stays research-lane, the amendment doc + phase report suffice |
| `.planning/ROADMAP.md` | Phase 13's SC list still says "reproduces `compute_BayesFactor()`" (`:261`) | **Should be annotated** with a pointer to the amendment, or `/gsd:verify-work` will score against superseded text |

### H2. The legitimacy test to reproduce (quoted from the precedent)

`07-GATE-AMENDMENT.md:22` states the test applied to every change:

> *"Would this change have been made, with this justification, by someone who had seen the gate's code
> but none of its outputs?"*

**Phase 13's answer is unusually strong and the amendment doc should say so plainly:** the
justification is not a Phase-13 result at all. It is `docs/amortized.md` named limit #3 — already
published, already part of the GO — establishing that the KDE reference is invalid past
`|logBF| ≈ log(L) ≈ 6.9` and that attrition is 100% baseline-side. Spike 014 quantified it on the
identical pairs: **34.75% non-finite (139/400), 32.5% saturated, dropped-pair median |Δρ| 0.736 vs
kept 0.283** [VERIFIED: `.planning/spikes/014-bf-sim-validation/README.md`]. Phase 13 is not amending
a gate because it failed; it is declining to adopt a reference already retired by a prior phase.
That is a materially cleaner position than Phase 7's, and the document should not blur the two.

### H3. What the amendment document must contain

1. **§0 Provenance disclosure** — state explicitly that **no Phase-13 result exists yet**, with the
   date and the commit sha, so the ordering is checkable.
2. **The defect in ROADMAP SC2**, argued only from the reference's validity conditions
   (L = 999 draws ⇒ smallest resolvable tail probability 1/999 ⇒ `|logBF| > 6.9` is kernel-tail
   extrapolation, not measurement — `references/bayes-factor-gate.md` §008).
3. **A second, independent defect D-02 introduces:** the Phase-13 net trains on the *Phase-11*
   registration-aware `zt` (129-row input surface), so it does not share an input space with
   `compute_BayesFactor()`'s Turing/ADVI path *or* with the shipped binary NRE. Even a valid KDE
   reference would not be a like-for-like comparison. This defect is architectural and
   outcome-independent.
4. **The replacement criterion, fully specified with numbers** (per-class AUC floor, confusion-matrix
   requirements, the calibration ECE band) — see §Validation Architecture.
5. **What is reported but not gated**: binary-NRE continuity at a fixed λ.
6. **An explicit statement that the original SC2 text is not erased** and that any report citing the
   amended criterion must cite the amendment alongside it (the `07-GATE-AMENDMENT.md:33-35` rule).

---

## I. The Calibration Gate (D-13) — Exact Signatures and Thresholds

### I1. `_bin_calibration` / `CalibrationResult` — two copies, both verified

| Copy | Location | Notes |
|---|---|---|
| Spike | `spike/validation/sbc.jl:65-72` (struct), `:82-130` (function) | **Use this one** — Phase 13 is spike-lane |
| Productionized | `test/gate/sbc.jl:44` (struct), `:60` (function) | Byte-equivalent pattern; `src/results.jl:117-126` carries the `CalibrationMeta` variant |

```julia
struct CalibrationResult          # spike/validation/sbc.jl:65-72
    bin_midpoints  :: Vector{Float64}
    predicted_rate :: Vector{Float64}
    observed_rate  :: Vector{Float64}
    bin_counts     :: Vector{Int}
    ece            :: Float64
    mce            :: Float64
end

_bin_calibration(posterior_probs :: AbstractVector{<:Real},
                 empirical_positive :: AbstractVector{Bool};
                 n_bins :: Int = 10) -> CalibrationResult   # :82-84
```

Semantics (read from `:85-130`): `n_bins` **equal-width** bins over `[0,1]`; last bin is
right-closed; empty bins contribute `predicted_rate = midpoint`, `observed_rate = 0.0`;
`ECE = Σ (count_i/total)·|pred_i − obs_i|`; `MCE = max_i |pred_i − obs_i|`.

**⚠ Empty-bin trap (high severity for D-13).** An empty bin is scored as `|midpoint − 0.0|` in the
**MCE** (it carries zero ECE weight but full MCE weight). With well-separated three-way classes many
of the 10 bins will be empty, so **MCE can be dominated by empty bins and is not a usable gate
statistic here.** Gate on **ECE**; report MCE with the empty-bin count beside it. Alternatively
pre-register `n_bins` low enough (e.g. 5) that every bin is populated — but do **not** hand-modify
`_bin_calibration`, which is shared, gate-lineage code.

### I2. Traffic-light thresholds

```julia
# spike/validation/consts.jl:49-50
const SBC_ECE_GREEN  = 0.05     # ECE ≤ GREEN  → :green
const SBC_ECE_YELLOW = 0.10     # ECE ≤ YELLOW → :yellow, else :red
# spike/validation/sbc.jl:139-143
sbc_traffic_light(ece) = ece <= SBC_ECE_GREEN ? :green :
                         ece <= SBC_ECE_YELLOW ? :yellow : :red
```

[VERIFIED: both files, 2026-07-25]. The same 0.05 / 0.10 values are re-declared in every
`test/gate/gate_consts_*.jl`. **Phase 13 should re-declare them in its own consts file rather than
importing them** — that is the established pattern (each gate carries its own frozen copy), and it
keeps the Phase-13 pre-registration self-contained.

### I3. What "implied class probabilities" means concretely

Each head's logit gives a two-class probability against random:
`p̂_C = σ(s_C(Z))`, `p̂_E = σ(s_E(Z))` — **the uncorrected logits**, since calibration is a
property of the classifier under its *training* frequencies. Two separate reliability curves:

| Curve | `posterior_probs` | `empirical_positive` | Evaluated on |
|---|---|---|---|
| Coloc head | `σ(s_C)` | `label == :coloc` | held-out samples with `label ∈ {:coloc, :random}` |
| Exclusion head | `σ(s_E)` | `label == :exclusion` | held-out samples with `label ∈ {:exclusion, :random}` |

Restricting each curve to its own head's class pair is not optional — it is the same D-11 argument:
feeding a coloc-head probability an exclusion-class sample asks the head about a class it never saw.

Spike 014's precedent for the analogous binary quantity: `reliability ECE of p̂(coloc) = 0.0188`
[VERIFIED: spike 014 README], comfortably green. That is the number Phase 13's coloc head should be
compared against in the continuity discussion.

**D-14 is correct:** ECE/MCE against a fixed band is not a hypothesis test, so there is no
over-power failure mode. No TOST machinery.

**One caveat D-14 does not cover, worth stating:** ECE is *under*-powered at small n and biased
downward when bins are sparse. Pre-register a minimum evaluation count (recommend `n ≥ 1000` per
head) so a green ECE is not a small-sample artifact.

---

## J. Exclusion Ground Truth (D-15, D-16)

### J1. The corpus audit — independently confirmed

`corpus/manifest.csv` has 36 lines (3 header comments + 1 column header + 32 data rows)
[VERIFIED: 2026-07-25]:

- **30** `simulated-secondary` CBS Red-Green rows (`cbs-RG-000` …), `truth_label = simulated-degree`,
  `coloc_ground_truth` a **colocalization percentage** — so the `0.0` end is *random*, not exclusion.
  D-15's reading is correct.
- **1** positive physical anchor: `pos-tetraspeck-01 | physical-primary | positive | coloc | 1.0`.
- **1** negative physical anchor: `neg-lightmycells-01 | physical-primary | negative | segregated | 0.0`,
  channels `2:nucleus-dna,mitochondria`.

**There is no graded segregation series.** D-15's premise holds and D-16's α-series is therefore
load-bearing evidence, exactly as CONTEXT says.

### J2. **BLOCKER — the physical anchor is not readable in Phase 13**

Three independent obstacles, all verified:

1. **`sha256 = "PENDING-FETCH"` on BOTH anchors** — the bootstrap fetch was deliberately not run
   (STATE Phase 08-05, T-08-16). `bootstrap_anchor_hashes()` exists but must be run online **and**
   authorized.
2. **`corpus/data/` is empty** [VERIFIED: `ls corpus/data` returns nothing] and `git ls-files corpus/data`
   is empty by design — the bytes are simply not on disk.
3. **Both anchors are `split = sealed_holdout`**, reachable **only** through
   `open_sealed_holdout(df; reason)` (`corpus/manifest.jl:149-157`), which is the **D-09 anti-snooping
   seal**. Opening it in Phase 13 spends a control that exists for Phase 16's blind evaluation.

**Recommendation:** Phase 13 should **NOT** open the seal. Plan the physical-anchor check as a
**declared deferral to Phase 16** with a one-paragraph note in the Phase-13 report, or as an explicit
`checkpoint:human-verify` if the user wants to authorize both the fetch and the seal break. Either
way the D-15 gate is unaffected: it rests on simulator ground truth, with the α-series as supporting
evidence. Both are fully available.

**Note on the download size:** the ~6.3 GB figure in STATE is the **positive** (TetraSpeck) anchor.
The segregated negative (`S-BIAD1047`, OME-TIFF from the EBI BioStudies FTP) is the one Phase 13
would want and is likely far smaller — but its size is recorded as `bytes = 0` (unfetched), so this is
**unverified**. If the user authorizes, fetching *only* the negative anchor is a much smaller ask than
STATE's headline number implies. [ASSUMED: relative sizes not measured.]

### J3. D-16's α-graded series — concrete design

**The trap that dominates this design.** `correlation()` calls `_exclude_zero`
(`src/colocalization.jl:154-169, 187-203`), which **drops any pixel where *either* channel is 0.0,
NaN, or missing**, and a patch with ≤ 15 survivors becomes `missing`
(`summary.jl:44`). Therefore:

> **Zeroing ch2 inside the ch1 mask does NOT create anti-correlation — it DELETES those pixels from
> the correlation entirely.** The measured patch correlation would then be computed over the
> complement only (≈ unchanged, i.e. *random*), and object-dense patches could fall below the
> 15-survivor floor and go `missing`, firing the summary's mask rows and plausibly the OOD flag.
> The α-ladder would look flat for a reason that has nothing to do with segregation.

**Recommended redistribution rule (intensity-preserving, zero-free, bitwise-exact at α = 0):**

```julia
# b :: strictly positive background floor. Recommended: b = quantile(vec(y), 0.05), asserted > 0.
#      (The spike simulator's own analogue is BG_FLOOR = 0.02, spike/simulator/forward.jl:66.)
# M :: the ch1 object mask, computed ONCE from the UNMODIFIED image and held FIXED across α.
function alpha_segregate(x, y, M, α; b)
    @assert b > 0 "background floor must be strictly positive (else _exclude_zero deletes pixels)"
    removed_px = α .* max.(y .- b, 0.0) .* M          # α = 0 ⇒ exactly 0.0, elementwise
    y_in       = y .- removed_px                       # α = 0 ⇒ y .- 0.0 == y BITWISE
    S_out      = sum(y[.!M]); @assert S_out > 0
    boost      = 1.0 + sum(removed_px) / S_out         # α = 0 ⇒ exactly 1.0
    y_out      = ifelse.(M, y_in, y .* boost)          # α = 0 ⇒ y .* 1.0 == y BITWISE
    return y_out
end
```

**Properties (each a testable assertion the plan should require):**

| Property | Assertion |
|---|---|
| `α = 0` reproduces the input **bitwise** | `alpha_segregate(x,y,M,0.0;b) == y` (exact `==`, not `≈`) — the Phase-11 D-10 "`ε=0` regression check" pattern |
| No pixel becomes zero | `all(alpha_segregate(...) .> 0)` given `b > 0` and `y[.!M] .> 0` |
| Total intensity preserved | `sum(y') ≈ sum(y)` to `rtol = 1e-10` — so the ladder is not confounded with a brightness change |
| Mask is α-invariant | `M` computed once from the unmodified `x` via `Images.otsu_threshold(x)`; **never** recomputed per α |
| Monotone | `m̄(α)` (mean patch correlation) should decrease monotonically; report the curve |

**Mask construction — reuse the frozen path.** `src/LoadImages.jl:235-241` already defines
`_calculate_mask(img) = [ch .> thr for (ch,thr) in zip(image_data(img), otsu_thresholds(img))]`, and
`Images.otsu_threshold` is the same call `LoadImages.jl:445` and `spike/contract.jl:68` use. Take
`M = _calculate_mask(mci)[1]` (the ch1 mask). This makes the segmentation rule a **frozen, already-shipped**
choice rather than a new one — which is the cheapest possible way to satisfy D-16's
"the segmentation/threshold rule is a pre-registered choice (D-04 applies)".

**Otsu is NOT re-consumed by the summary** — `patch_summary` → `patch` → `correlation` read
`mci.data` directly and never touch `otsu_threshold` [VERIFIED: grep over `src/colocalization.jl` and
`src/amortized/*.jl` shows Otsu used only in `api.jl:58`, `local_map.jl:208`, `simulator.jl:268`, i.e.
only at MCI *construction*]. So rebuilding the modified MCI (which recomputes ch2's Otsu) does not
perturb the summary. One trap fewer — but the plan should still hold the *mask* fixed across α, per
the `local_map.jl:208` precedent ("tile-local Otsu would make tiles incomparable").

**What to run the α-series on.** Two substrates, both worth doing:

| Substrate | Why | Availability |
|---|---|---|
| **Simulated images at ρ ≈ 0** (drawn through the Phase-11 simulator) | full control, arbitrary n, known nuisances, guaranteed positive background (`BG_FLOOR = 0.02`), and `α = 0` is *by construction* a random pair | ✓ available |
| **CBS `cbs-RG-000`** (the 0.0-degree Red-Green row) | real optics/noise; the manifest's own "random" end; `split = eval`, **not** sealed | ✓ available in principle, but `sha256` unfilled for CBS rows too — needs the fetch path. Verify before planning it as a hard deliverable. |

**The scientifically interesting readout** (CONTEXT §Specific Ideas): not merely "does exclusion get
detected", but **at what α the intensity-correlation notion starts to track disjoint localization**.
Report `log BF(E:R)(α)` and `log BF(C:R)(α)` as curves, plus `m̄(α)`, plus a **crossing point**
`α*` = the smallest α at which `log BF(E:R) > 0`. That single number is the phase's most quotable
result about the correlation-vs-localization gap.

---

## K. Pre-Registration: Seeds and the Phase-13 Consts File (D-01, D-04)

### K1. The complete forbidden-seed inventory (enumerated, as requested)

Assembled from `grep -rhoE "0x[0-9a-fA-F_]{4,}"` across `*.jl` plus the derived gate seeds. This is
the superset Phase 13 must assert against.

| Seed / salt | Value | Origin |
|---|---|---|
| `DEFAULT_MASTER_SEED` | `0x0000000000000001` | productionization datagen |
| `NPE_MASTER_SEED` | `0x0000000000c0ffee` | `spike/npe/train_npe.jl:65` |
| `VAL_MASTER_SEED` | `0x000000005bc0ffee` | `spike/validation/consts.jl:76` |
| `VAL_FIX_SEED` | `0x0000000000f1f7ed` | `spike/validation/consts.jl:79` (literal is `0xF1F7ED`) |
| `RATIO_PAIR_SEED` | `0x00000000004a7107` | `src/amortized/train_ratio.jl:60` |
| `CORPUS_MASTER_SEED` | `0x0000000000c05eed` | `corpus/manifest.csv` header |
| F2 dev seeds | `0x000000000de7c0de`, `0x00000000de7c0de2`, `0x000000000de7c0d3` | 07-CALIBRATION-FINDINGS |
| spike 006/008 dev | `0x00000000000b7a770/71/72/73` | `sources/006…`, `sources/008…` |
| spike 010 dev | `0x0000000000ca11b0/b1/b2/b3` | `sources/010…` |
| spike 013 train | `0x0000000000dec0de` | `sources/013-mu-prior-truncation/train_trunc.jl:52` |
| **spike 014 dev** | **`0x0000000000bf5014`** | `.planning/spikes/014-bf-sim-validation/run_sim_bf.jl:38` |
| misc burned | `0x0000000000a1ce01/02`, `0x00000000c0ffee5b`, `0x000000c0ffee5bc0`, `0x0000000000a70115`, `0x0000000000a70213`, `0x000000000134d8f3` | assorted spikes / gate scripts |
| **Phase-11 proposed** | **`0x000000000b11de71`** (`P11_DEV_SEED`), salt `0xa24baed4663ee121` | `11-RESEARCH.md:1651-1652` — **must be forbidden even though Phase 11 has not run yet** |
| 11-RESEARCH probe burns | `0xBEEF`, `0xCAFE`, `0xFEED`, `0xDEAD`, `Random.seed!(7)` | `11-RESEARCH.md:1672-1673` |
| Salts | `0x9e3779b97f4a7c15` (HOLDOUT), `0xd1b54a32d192ed03` (FOLD / CBS_SPLIT), `0xbf58476d1ce4e5b9` (VAL), `0x94d049bb133111eb` (PROD), `0xc4ceb9fe1a85ec53` (AMEND), `0xa5a5a5a5a5a5a5a5`, `0xa5a5a5a5deadbeef` (COSTES) | various |
| **Derived, must be recomputed** | `PROD_SEED[G]` and `PROD_SEED_V2[G]` for `G ∈ (4,8,16,32)` | `test/gate/gate_consts_8_v2.jl:275-285, 315-331` |

`PROD_SEED_V2` is the one D-01 names explicitly. It is **derived, not literal** — the only correct way
to forbid it is to recompute it, exactly as spike 014 did by loading the frozen consts into an
isolated module:

```julia
module _GC
    include(joinpath(@__DIR__, "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end   # READ-ONLY: no gate runs, the file is not modified (spike 014 precedent, run_sim_bf.jl:42-45)
```

### K2. Proposed Phase-13 seeds (verified absent from the whole repo)

```julia
const P13_DEV_SEED = 0x0000_0000_0B13_DE71   # "0B13" = phase 13, "DE71" = DEV-1 (mirrors P11)
const P13_SALT     = 0x2545_F491_4F6C_DD1D   # xorshift64* multiplier; ≠ every salt in the repo
```

[VERIFIED: `grep -ril` over `*.jl`, `*.md`, `*.toml` returns **zero** hits for `0b13de71`,
`b13de71`, `2545f4914f6cdd1d` in either case, 2026-07-25]

### K3. The Phase-13 consts file — required contents

Mirror `spike/validation/consts.jl` structurally (anti-snooping preamble, decoupling note,
`if !isdefined(@__MODULE__, :SENTINEL) … end` re-inclusion guard, grouped sections with inline
rationale) **and** adopt the three stronger `gate_consts_8_v2.jl` habits: recompute forbidden seeds
rather than trusting comments, carry a `_p13_forbidden()` function, and end with executable
`@assert` self-checks.

**Everything below must be frozen BEFORE the first run.** The τ probe consumes a stream, so its seed
and its own specification are Tier-1 (a seed chosen after seeing probe output is a snooped seed).

| Tier | Constant | Notes |
|---|---|---|
| **1 — before ANYTHING runs** | `P13_DEV_SEED`, `P13_SALT`, `p13_rng(counter)` | the probe itself consumes a stream |
| 1 | `P13_TAU_DELTA_GRID`, `P13_TAU_AUC`, `P13_TAU_R`, τ statistic (= `m̄`), unpaired design, reference λ, both-directions rule | else the τ measurement is selectable after the fact |
| 1 | **τ probe abort criterion** (`A(δ) < AUC_τ` for every δ ⇒ stop and report) | must be immune to its own result |
| 1 | `P13_CUT_VARIANT` (`:tau_contrast` or `:sign_contrast`, §E1) | changes every label |
| 1 | `P13_STRATIFICATION` (`:class_frequency` = D-07-i, or `:importance_weighted` = D-07-ii) | changes the correction's validity |
| 1 | §F5 verification pass bars (`cor ≥ 0.99`, `max|Δ| ≤ 0.25`) **and** the negative-control requirement | else the verification is unfalsifiable |
| 1 | D-12 gate numbers: per-class AUC floor, confusion-matrix requirements, `M` (label count) | the reported gate |
| 1 | D-13 numbers: `P13_ECE_GREEN = 0.05`, `P13_ECE_YELLOW = 0.10`, `n_bins`, min-n per head, **ECE-not-MCE gate rule** (§I1) | |
| 1 | D-16 α ladder, background-floor rule (`b = quantile(y, 0.05)`), mask rule (`_calculate_mask(mci)[1]`), the three invariants of §J3 | the "pre-registered segmentation choice" D-16 requires |
| 1 | `P13_ITERATION_ALLOWANCE = 1` **plus the pre-declared trigger** for spending it (recommend: "exclusion-class per-class AUC below floor *specifically at high \|ρ\|*" ⇒ switch to D-07-ii) | D-04's novelty; naming the trigger in advance is what stops it becoming "amended twice" |
| 2 — legitimately probe-derived | **`P13_TAU` itself** | this is the one number D-06 says must be *measured*. It is frozen into the consts file **after the probe and before any labelled data is generated** (D-06's exact wording), in a clearly-marked, separately-committed block |
| 2 | Training recipe values **copied, not chosen**: `epochs=300, batchsize=128, lr=2.5e-4, wd=1e-4, val_frac=0.15, stopping_epochs=40, n=48_000` | read from `train_ratio.jl:155-163`; copying them verbatim is what D-10's attribution argument rests on |

**The two-commit discipline for τ** is important and should be explicit in the plan: commit the consts
file with `P13_TAU = nothing` and the probe spec frozen → run the probe → commit `P13_TAU = <measured>`
with the probe log attached → *only then* generate labelled data. Any other order breaks D-06.

---

## Don't Hand-Roll

| Problem | Don't build | Use instead | Why |
|---|---|---|---|
| Training loop with early stopping + LR schedule | A bespoke Flux epoch loop | NeuralEstimators `train` via a `<: NeuralEstimator` subtype (§A4/§B1) | Reproducing CosAnneal + patience + best-val-checkpoint *identically* is the only way D-10's "difference attributable to the head alone" survives; a re-implementation is an uncontrolled variable |
| Binary cross-entropy on logits | `log(1+exp(-x))` by hand | `Flux.logitbinarycrossentropy` | numerically stable branchless form; hand-rolled versions overflow at \|logit\| > 30, which is exactly the confident-call tail this phase is about |
| ROC / AUC | A new package or a fresh implementation | `roc_auc` in `src/amortized/ood.jl` (hand-rolled, already gate-lineage) | `spike/test/runtests.jl` gate (h) **asserts `ROCAnalysis`/`MLJ` are absent**; adding one breaks the suite |
| Reliability diagram / ECE / MCE | A new binning routine | `_bin_calibration` (`spike/validation/sbc.jl:82`) | D-13 names it; it carries the gate lineage; a second implementation would need its own validation |
| Object mask from an image | A new thresholding rule | `_calculate_mask` / `Images.otsu_threshold` (`src/LoadImages.jl:235-241, 445`) | makes D-16's segmentation rule a *frozen, already-shipped* choice, satisfying D-04 at zero cost |
| Pair encoding | A new concat scheme | `pair_encode` (`bf.jl:91-100`) | grid-general, shared by trainer and reader, and its `iseven` assert catches the λ-placement bug (§D2) |
| Reproducible per-index RNG | `Random.seed!(i)` in a loop | `Philox4x(UInt64,(seed, counter))` via `spike/data/seeding.jl` | thread-count-independent byte-identical generation (STATE 03-02); `seed!`-in-loop is not |
| Truncated-Cauchy CDF / ρ-prior density | A numeric quadrature | closed forms in §E2 / §F5 | `QuadGK` is **not** in the spike env; the closed forms are exact and atom-aware |
| Tail probabilities of a posterior | KDE + `quadgk` | nothing — **do not** | Named limit #3; the reference is retired and its packages are not in the spike env |

**Key insight:** almost every "new" component this phase seems to need already exists in the repo with
gate lineage attached. The genuinely new code is small: the two-head struct, the masked loss, the
three-way label rule, `measure_head_log_odds`, the τ probe, and the α-series. Everything else is reuse.

---

## Common Pitfalls

### Pitfall 1: Copying `− log_prior_odds` from the binary formula
**What goes wrong:** the head's correction is applied with the wrong sign, the wrong denominator
(all three classes instead of the head's pair), or twice.
**Why it happens:** `amortized_log_bf` looks like the obvious template, and its `log_prior_odds` is
≈ 0 so a wrong copy is invisible in the binary case.
**How to avoid:** derive from §F1; use `measure_head_log_odds(labels; positive, negative)` with the
head's own pair; run the §F5 verification with a *deliberately skewed* toy and a negative control.
**Warning signs:** corrected log-BFs off by a constant ≈ 0.5 nat; the coloc and exclusion heads
disagreeing by exactly `log(q_E/q_C)` on random-class inputs.

### Pitfall 2: Zeroing ch2 inside the mask in the α-series
**What goes wrong:** `_exclude_zero` deletes those pixels, so the ladder is flat and patches may go
`missing`; the phase concludes "the net cannot see spatial segregation" when it never saw any.
**Why it happens:** "redistribute signal out of the mask" reads as "set to zero".
**How to avoid:** the strictly-positive background floor `b` in §J3; assert `all(y' .> 0)`.
**Warning signs:** rising `missing`-patch counts with α; mask rows of the summary changing with α;
the OOD flag firing on the α-ladder.

### Pitfall 3: Building the three-way label on Δρ
**What goes wrong:** a strongly colocalized sample under a *more* colocalized control is published as
"mutually exclusive."
**Why it happens:** `delta_rho` (`infer.jl:117-121`) is the natural handle and the shipped net uses
exactly that quantity.
**How to avoid:** D-05's two-factor cut on `ρ_sample` **level** × control contrast; assert in a unit
test that `(ρ_s, ρ_c) = (0.8, 0.9)` labels `:random` (or `:coloc` under variant (a) if the level
condition dominates) and **never** `:exclusion`.
**Warning signs:** exclusion-class members with `ρ_sample > 0`.

### Pitfall 4: Appending λ in the wrong place (or twice)
**What goes wrong:** input width 322 or 323; λ triple-counted through `pair_encode`'s contrast block;
or λ z-scored by a `zt` that never saw it.
**How to avoid:** §D2's `vcat(pair_encode(Zs_128, Zc_128), λ)`; assert
`size(Z_pair,1) == ratio_input_dim(G) + n_cond`.
**Warning signs:** a λ-response that is exactly flat (the net learned to ignore a nonsense row), or a
`pair_encode` `ArgumentError` about an odd summary length (this one is the *good* failure).

### Pitfall 5: Re-fitting `zt` instead of inheriting it
**What goes wrong:** the ratio net's input space silently differs from the NPE's; every downstream
comparison is meaningless.
**Why it happens:** the Phase-13 training pool is new, so "fit the standardizer on it" feels natural.
**How to avoid:** `train_ratio.jl:144-146`'s discipline — the ratio net *always* trains on an
already-frozen `zt` from a trained NPE. Assert the loaded `zt` object identity is the Phase-11 one.
**Warning signs:** standardized summaries with mean ≈ 0 / sd ≈ 1 on the Phase-13 pool (a frozen
transform applied to a *different* pool should NOT give exactly that).

### Pitfall 6: Global-RNG leakage breaking reproducibility
**What goes wrong:** two runs on the same Philox seed disagree.
**Why it happens:** `Flux.DataLoader(shuffle=true)` — used by NeuralEstimators' `_DataLoader`
(`train.jl:668-676`) — draws from the **global** RNG, and `sampleposterior` threads no `rng`
(documented at length in `test/gate/harness.jl:45-72`). This has already been a real regression in
this repo.
**How to avoid:** `Random.seed!(derived_from(P13_DEV_SEED))` immediately before every `train` call and
every posterior-draw loop, in addition to the Philox streams.
**Warning signs:** a twin-run bit-reproducibility assert failing only when threads > 1 or only on the
training arm.

### Pitfall 7: Gating on MCE
**What goes wrong:** empty reliability bins are scored `|midpoint − 0|`, so MCE ≈ 0.95 on a
well-separated three-way problem and the gate fails for a reason unrelated to calibration.
**How to avoid:** §I1 — gate on ECE, report MCE with the empty-bin count.
**Warning signs:** MCE ≫ ECE with several `bin_counts .== 0`.

### Pitfall 8: Measuring τ with a paired design
**What goes wrong:** τ comes out several-fold too small; "random" becomes a razor-thin band; the
random class is starved and the confusion matrix looks great for the wrong reason.
**How to avoid:** §E3's unpaired specification, nuisances and imsize sampled from their full priors.
**Warning signs:** measured τ below the `ghat` knot spacing near zero (≈ 0.08 in ρ) — that would be
claiming sub-knot resolution from a piecewise-linear calibration.

### Pitfall 9: Treating a `vacuous` clean result as evidence
**What goes wrong:** a head that learns nothing produces well-calibrated-looking probabilities by
reproducing the base rate — the exact F3 "vacuous SBC pass" failure, one level up.
**How to avoid:** report each head's AUC **beside** its ECE. A green ECE with AUC ≈ 0.5 is a vacuous
pass and must be labelled as such (`07-CALIBRATION-FINDINGS.md` F3 precedent).
**Warning signs:** ECE green, AUC ≈ 0.5, predicted probabilities clustered at the class base rate.

### Pitfall 10: Running the τ probe against the pre-Phase-11 simulator
**What goes wrong:** τ is measured under `SHIFT_PRIOR = Uniform(-1,1)` with no `ε`, i.e. an easier
world than the net will train in; the frozen τ is optimistic and the random band is too narrow.
**How to avoid:** §E4 — gate the probe on the Phase-11 *simulator merge*, and record which simulator
commit the probe ran against in the consts file.

---

## Code Examples

### 1. The custom estimator, end to end (D-09 / D-10 / D-11)

```julia
# spike/p13/net.jl
import NeuralEstimators: NeuralEstimator, train
import Flux
import Flux: Chain, Dense, gelu, logitbinarycrossentropy

# Trunk constants COPIED VERBATIM from src/amortized/train_ratio.jl:53-54,77-78 (D-10).
const P13_SUMMARY_WIDTH = 256
const P13_NUM_SUMMARIES = 64

struct ThreeWayEvidenceNet <: NeuralEstimator
    trunk::Any
    head::Any
end
Flux.@layer ThreeWayEvidenceNet
(m::ThreeWayEvidenceNet)(Z) = m.head(m.trunk(Z))      # 2×n logits — ONE forward pass (SC1)

function build_three_way_net(input_dim::Integer;
                             num_summaries::Integer = P13_NUM_SUMMARIES,
                             width::Integer = P13_SUMMARY_WIDTH)
    W = width
    trunk = Chain(Dense(input_dim, W, gelu), Dense(W, W, gelu),
                  Dense(W, W, gelu), Dense(W, num_summaries))   # verbatim train_ratio.jl:77-78
    return ThreeWayEvidenceNet(trunk, Dense(num_summaries, 2))
end

# D-11: per-head class restriction via participation weights; joint trunk via the summed loss.
function masked_two_head_bce(ŷ, y)          # ŷ :: 2×n, y :: 4×n = [y_C; w_C; y_E; w_E]
    lC = logitbinarycrossentropy.(view(ŷ, 1, :), view(y, 1, :); agg = identity)
    lE = logitbinarycrossentropy.(view(ŷ, 2, :), view(y, 3, :); agg = identity)
    wC = view(y, 2, :); wE = view(y, 4, :)
    return (sum(wC .* lC) + sum(wE .* lE)) / (sum(wC) + sum(wE) + eps(Float32))
end

# Recipe COPIED VERBATIM from train_ratio.jl:155-183. AdamW args stay Float64 (CosAnneal gotcha).
est = train(build_three_way_net(input_dim),
            y_tr, y_va, Z_tr, Z_va;
            loss = masked_two_head_bce,          # honoured: generic _loss (train.jl:654)
            epochs = 300, batchsize = 128, use_gpu = false,
            optimiser = Flux.Optimisers.AdamW(2.5e-4, (0.9, 0.999), 1e-4),
            stopping_epochs = 40, verbose = true)
```
*Sources: `src/amortized/train_ratio.jl:53-54,77-78,155-183`; NeuralEstimators
`src/train.jl:101,654,657`; `ext/NeuralEstimatorsFluxExt.jl:39-41,67,80`. All [VERIFIED].*

### 2. The three-way label rule (D-05, variant (a))

```julia
# spike/p13/labels.jl — P13_TAU is READ from the frozen consts file, never inlined.
@enum ThreeWayClass EXCLUSION RANDOM COLOC

function three_way_label(ρ_sample::Real, ρ_control::Real; τ::Real = P13_TAU)
    d = ρ_sample - ρ_control
    (ρ_sample >  τ && d >  τ) && return COLOC
    (ρ_sample < -τ && d < -τ) && return EXCLUSION
    return RANDOM
end

# D-05's counter-example, as an executable regression assertion:
@assert three_way_label(0.8, 0.9) != EXCLUSION   # "less colocalized than control" is NOT segregation

# D-11 participation weights.
function head_targets(c::ThreeWayClass)
    c === COLOC     && return Float32[1, 1, 0, 0]   # y_C=1, w_C=1, exclusion head inactive
    c === EXCLUSION && return Float32[0, 0, 1, 1]   # exclusion head only
    return Float32[0, 1, 0, 1]                       # RANDOM: negative for BOTH heads
end
```

### 3. The τ probe (D-06), fit-free and unpaired

```julia
# spike/p13/tau_probe.jl — needs NO trained network (see §E4).
mbar(Z128) = mean(@view(Z128[1:64])[@view(Z128[65:128]) .> 0.5])   # present continuous rows only

function summary_draws(ρ::Real, R::Integer, rng_base)
    map(1:R) do r
        rng = Philox4x(UInt64, (rng_base, UInt64(r)))
        θ   = merge(sample_prior(rng), (; ρ_true = ρ))       # nuisances ~ π, ρ pinned
        y   = simulate_pair(rng, θ; imsize = sample_imsize(rng))   # F5 mixture, NOT 256²
        mbar(encode_d01(patch_summary(build_mci(y))))
    end
end

function tau_auc(δ; R = P13_TAU_R, rng_base = p13_rng_key(:tau))
    a = summary_draws( 0.0, R, rng_base)
    b = summary_draws(  -δ, R, rng_base + 0x1000)   # DISJOINT stream — unpaired by design
    return roc_auc(vcat(a, b), vcat(falses(R), trues(R)))
end

τ_measured = findfirst(δ -> max(tau_auc(δ), tau_auc_pos(δ)) ≥ P13_TAU_AUC, P13_TAU_DELTA_GRID)
```

### 4. Closed-form verification of the D-07 correction (§F5)

```julia
import Distributions: Normal, cdf, logpdf
# Z|θ ~ N(θ,σ²), θ ~ N(0,s²); classes E={θ<-τ}, R={|θ|≤τ}, C={θ>τ}
function log_evidence(Z, a, b; s, σ)                  # log p(Z | θ ∈ [a,b])
    v = s^2*σ^2/(s^2+σ^2); m = Z*s^2/(s^2+σ^2); sd = sqrt(v)
    num = cdf(Normal(m,sd), b) - cdf(Normal(m,sd), a)
    den = cdf(Normal(0,s),  b) - cdf(Normal(0,s),  a)
    return logpdf(Normal(0, sqrt(s^2+σ^2)), Z) + log(num) - log(den)
end
analytic_log_bf_CR(Z; τ, s, σ) = log_evidence(Z, τ, Inf; s, σ) - log_evidence(Z, -τ, τ; s, σ)
# PASS bars are pre-registered; the negative control (uncorrected logit) must FAIL them.
```

---

## State of the Art

| Old approach (in this repo) | Current approach | When changed | Impact on Phase 13 |
|---|---|---|---|
| Bayes factor validated by magnitude agreement with a KDE + `quadgk` reference (`compute_BayesFactor`, `src/bayes.jl:109`) | **Simulation-based discrimination**: AUC + monotonicity + decision-calibration on labelled simulated pairs, no per-pair reference | Spikes 006–008 (diagnosis) → spike 014 (VALIDATED, AUC 0.994) → adopted by the 2026-07-24 GO as named limit #3 | D-12's gate *is* this method, extended from 2 classes to 3 |
| Raw SBC ranks | **Randomized ranks** for atom-bearing parameters | Spike 010 addendum; named limit #2 | Not directly used (Phase 13 has no SBC arm) but the *reason* — prior atoms — lands hardest on the exclusion end (§E2) |
| Clamped KDE tail (`_clampp = 1e-8`) | Non-clamped `kde_log_bf_unclamped` (`bf.jl:139-149`) | T-7-06 / Memo §5 | Historical only; **not used in Phase 13**, and its packages are absent from the spike env |
| Capacity as the calibration lever | Capacity **redistributes** error, does not remove it | Spike 012 (PARTIAL) | D-10's justification for reusing the trunk verbatim |
| `Pkg.develop` coupling of spike ↔ src | Read-only `include(joinpath(@__DIR__,"..","src",...))` | Phase 01-04 (the NeuralEstimators downgrade landmine) | The mechanism §G3 recommends for reaching `AbstractColocResult` from `spike/` |

**Deprecated / do not reintroduce:** the clamped KDE baseline; `compute_BayesFactor` as a validation
reference; μ-prior truncation as an atom fix (spike 013 — VALIDATED but costs π/Turing consistency,
and is explicitly deferred by CONTEXT).

---

## Assumptions Log

| # | Claim | Section | Risk if wrong |
|---|---|---|---|
| A1 | Variant (a) (τ on both the level and the contrast) is the right operationalization of "stands apart from the control" | §E1 | Class masses and every label shift; the D-12 confusion matrix is not comparable to a run under variant (b). **Mitigation: pre-register the choice; it is cheap to compute both mass tables before committing.** |
| A2 | `AUC_τ = 0.90` is the right "reliably distinguish" bar | §E3 | τ too small (random band razor-thin, random class starved) or too large (real signal swallowed by "random"). **Mitigation: report the whole `A(δ)` curve so a reader can re-read τ at any bar.** |
| A3 | `m̄` (mean present continuous rows) is the right τ statistic | §E3 | A more powerful readout would give a smaller τ; the reported τ is then a *summary-scalar* resolution, not the summary's full resolution. **Mitigation: report the NPE-based τ as a secondary column, labelled as model-dependent.** |
| A4 | D-07-i (class-frequency stratification) suffices for coverage | §F3 | If the exclusion head is under-resolved at high \|ρ\|, the phase spends its one iteration. **Mitigation: pre-declare that exact trigger.** |
| A5 | R1 (`<: NeuralEstimator` subtype) works with NeuralEstimators 0.2.1's internal hooks | §B1 | Falls back to R2 (hand-rolled loop), weakening D-10's attribution argument. **Mitigation: a two-epoch Wave-0 smoke test on a 200-sample toy, before any real training.** |
| A6 | The negative physical anchor (`S-BIAD1047`) is much smaller than the 6.3 GB positive anchor | §J2 | The physical check stays deferred (which is already the recommendation). Low impact. |
| A7 | Phase 11 will deliver a λ-conditioned net whose λ input is actually informative | §D3 | Phase 13's D-03 conditioning is vacuous too; the log-BF λ-response is flat. **Mitigation: measure and report the response honestly; a flat response is a finding.** |
| A8 | `P13_DEV_SEED = 0x0B13DE71` / `P13_SALT = 0x2545F4914F6CDD1D` remain unused when Phase 13 executes | §K2 | Seed collision → snooping exposure. **Mitigation: the executable `_p13_forbidden()` assert re-checks at runtime, including recomputed `PROD_SEED_V2`.** |
| A9 | Two λ values (one per stack) are not wanted; a single λ describes the acquisition | §D2 | Input width off by one and a modelling mismatch with Phase 11. **Mitigation: state it as a pre-registered modelling choice.** |

---

## Open Questions

1. **Does the comment-only `log_bf_simplex` rename land in `src/results.jl` during Phase 13?**
   - What we know: it is comment-only, zero provenance cost (`results.jl` ∉ `DATAGEN_HASH_SRC_FILES`),
     zero runtime impact; but it *is* a `src/` edit and CLAUDE.md says `src/` stays provably untouched
     during the spike.
   - What's unclear: whether D-01's enumerated scope (artifact / accessor / GO) was meant to be
     exhaustive.
   - Recommendation: **do not edit `src/` in Phase 13**; define the result type in `spike/`; record the
     misnomer in the Phase-13 report and carry the rename as a one-line item for productionization.
     Surface as a `checkpoint:human-verify` if the user wants it now.

2. **Which λ does the frozen scalar τ correspond to?**
   - What we know: resolution degrades with registration uncertainty, so τ is really `τ(λ)`.
   - What's unclear: whether D-05/D-06 intend a single scalar (their wording suggests yes).
   - Recommendation: measure `τ(λ)` across the ladder, **freeze the scalar at a pre-registered
     reference λ** (recommend the widest rung — most conservative), report the curve.

3. **Does Phase 13 need the CBS `cbs-RG-000` bytes for the α-series?**
   - What we know: CBS rows are `split = eval` (not sealed), but their `sha256` column is empty in the
     committed manifest and `corpus/data/` is empty.
   - What's unclear: whether the CBS fetch path is runnable without the anchor bootstrap.
   - Recommendation: plan the α-series on **simulated** substrate as the deliverable; treat the CBS
     arm as an optional extension gated on a cheap fetch check.

4. **Is `ε` (chromatic) also a conditioning input in Phase 11's net?**
   - What we know: 11-RESEARCH D13 says `ε` is deliberately *not* λ-scaled and flags a possible second
     conditioning input in its Open Questions.
   - Recommendation: derive Phase-13's input width as `ratio_input_dim(G) + length(encode_lambda(...))`;
     **never write 321 as a literal.**

5. **Should the confusion matrix be reported at `argmax` of the two corrected log-BFs, or at
   `logBF > 0` per head?**
   - What we know: with two heads and a structural zero for random, "declare the argmax of
     `(logBF_C, 0, logBF_E)`" is the natural rule and needs no extra threshold. Spike 014 used the
     analogous `logBF > 0` rule for the binary case.
   - What's unclear: the argmax rule can never declare `random` unless both log-BFs are negative,
     which is a *decision* rule — and decisions are Phase 14's job.
   - Recommendation: gate on **threshold-free** statistics (per-class AUC, one-vs-random) and report
     the argmax confusion matrix as a **descriptive** companion, explicitly labelled as not a decision
     rule. This keeps the Phase-13/Phase-14 boundary clean.

---

## Environment Availability

| Dependency | Required by | Available | Version | Fallback |
|---|---|---|---|---|
| Julia | everything | ✓ | 1.12.6 (pinned, `spike/Manifest.toml:3`) | — |
| NeuralEstimators | `train` loop, `logratio` | ✓ | 0.2.1 (pinned) | R2 hand-rolled loop |
| Flux | trunk, loss, optimiser | ✓ | 0.16.10 (pinned) | — |
| Images / ImageFiltering | D-16 Otsu mask | ✓ | 0.26.2 / 0.7.12 | — |
| Distributions / StatsBase / JLD2 / Random123 / HypothesisTests / CairoMakie | priors, stats, persistence, seeding, figures | ✓ | resolved in `spike/Manifest.toml` | — |
| CUDA | optional training acceleration | ✗ | — | CPU-only is the baseline (CLAUDE.md); `spike/test/runtests.jl` **asserts CUDA absent** — do not add it |
| KernelDensity / QuadGK | (only the retired KDE baseline) | ✗ in `spike/` | — | **Not needed** — D-12 declines to gate on it |
| **Phase-11 research NPE + `zt` + `encode_lambda`** | **D-02 / D-03 — every labelled-data task** | ✗ | — | **NONE. This is the hard block.** Do not substitute the shipped grid-8 basis. |
| Phase-11 simulator surgery (`ε`, composed warp, widened `SHIFT_PRIOR`) | τ probe (§E4) | ✗ | — | Probe may run against the old simulator **only** as a throwaway sanity run whose τ is NOT frozen |
| `corpus/data/` bytes (physical anchors) | D-15's named qualitative check | ✗ (`sha256 = "PENDING-FETCH"`, dir empty, `sealed_holdout`) | — | Defer to Phase 16 (recommended) or an authorized fetch + seal break |

**Missing dependencies with no fallback (block execution):**
- Phase-11 research NPE bundle (net + `zt` + `θzt` + λ encoding + imsize provenance).

**Missing dependencies with fallback:**
- Phase-11 simulator surgery → τ probe deferred one wave, not blocked outright.
- Physical anchors → declared deferral to Phase 16; gate is unaffected.

---

## Validation Architecture

### Test Framework

| Property | Value |
|---|---|
| Framework | Julia stdlib `Test` (`@testset` / `@test`) — no external test package |
| Config file | none (Julia convention); the spike gate entry point is `spike/test/runtests.jl` |
| Quick run command | `julia --project=spike spike/test/runtests.jl` |
| Full suite command | `julia --project=spike spike/test/runtests.jl` (spike lane) **and** `julia --project=. -e 'using Pkg; Pkg.test()'` (root — must stay green because `src/` is untouched) |
| Wiring convention | new files are `include`d at the bottom of `spike/test/runtests.jl` so **one command stays the gate** (precedent: `test_simulator.jl`, `test_data_pipeline.jl`, `test_npe.jl`, `test_sbc.jl`, `test_bf.jl`, `test_ood.jl`) |
| Reported-run convention | long/reported runs are **scripts**, not tests (`spike/validation/run_*.jl` pattern); the test suite carries fast fixtures on `VAL_FIX_SEED`-style streams that never consume the reported stream |

### Phase Criteria → Test Map

No requirement IDs exist (§Phase Requirements), so the map is keyed to the amended success criteria
and the decisions that carry testable content.

| Criterion | Behaviour | Test type | Automated command | File exists? |
|---|---|---|---|---|
| SC1 (amended) | `ThreeWayEvidenceNet` returns a 2×n logit matrix in one call; input width `== ratio_input_dim(G) + n_cond` | unit | `julia --project=spike spike/test/runtests.jl` (`@testset "P13 net shape"`) | ❌ Wave 0 |
| D-09/A5 | NeuralEstimators `train` accepts the custom estimator + custom loss; 2 epochs on a 200-sample toy reduce the validation risk | smoke | same, `@testset "P13 train smoke"` | ❌ Wave 0 |
| D-11 | `head_targets` gives `(1,1,0,0)/(0,0,1,1)/(0,1,0,1)`; random is negative for **both** heads; the masked loss ignores inactive cells (perturbing an inactive logit leaves the loss bit-identical) | unit | same, `@testset "P13 masked loss"` | ❌ Wave 0 |
| D-05 | `three_way_label(0.8, 0.9) != EXCLUSION`; sign symmetry; the τ dead-zone maps to `RANDOM` | unit | same, `@testset "P13 labels"` | ❌ Wave 0 |
| D-06 | τ probe is deterministic under a fixed key; unpaired streams are disjoint; `mbar` ignores mask-absent patches | unit (fixture-scale R) | same, `@testset "P13 tau probe"` | ❌ Wave 0 |
| **D-07** | corrected logits recover the **analytic** log BF on the closed-form Gaussian toy (`cor ≥ 0.99`, `max\|Δ\| ≤ 0.25`) **and** the uncorrected logit FAILS the same bars | unit (fast, ~seconds) | same, `@testset "P13 prior-odds correction"` | ❌ Wave 0 |
| D-07 (§F2) | scalar-only correction **fails** under a within-class reshaped proposal | unit | same, `@testset "P13 stratification theorem"` | ❌ Wave 0 |
| D-02/D-03 | precondition script errors with the Phase-11 pointer when the net is absent; asserts 8-row θ, λ range, imsize provenance | unit | same, `@testset "P13 preconditions"` | ❌ Wave 0 |
| D-03 | `Z_pair` width `== ratio_input_dim(G) + n_cond`; a naive 129-row `pair_encode` throws | unit | same | ❌ Wave 0 |
| D-13 | `_bin_calibration` on synthetic perfectly-calibrated input is green; on degenerate input is red; **empty-bin MCE trap** is asserted (documenting why ECE is the gate) | unit | same, `@testset "P13 calibration"` | ⚠ partial — `spike/test/test_sbc.jl:77-86` already covers the base function |
| **D-16** | `alpha_segregate(x,y,M,0.0;b) == y` **bitwise**; `all(y' .> 0)`; `sum(y') ≈ sum(y)`; mask α-invariant; `m̄(α)` monotone non-increasing on a fixture | unit | same, `@testset "P13 alpha series"` | ❌ Wave 0 |
| D-01/D-04 | `P13_DEV_SEED ∉ _p13_forbidden()` including **recomputed** `PROD_SEED` / `PROD_SEED_V2`; `P13_SALT` distinct from all repo salts; consts file re-include is a no-op | unit | same, `@testset "P13 seeds"` | ❌ Wave 0 |
| D-01 | `src/` provably untouched — `git diff --quiet HEAD -- src/` at the start of the reported run (precedent: `spike/demo.jl` step-0 assertion) | integration | same | ❌ Wave 0 |
| Decoupling | `spike/Project.toml`/`Manifest.toml` unchanged; `NeuralEstimators == v"0.2.1"` by UUID; `ROCAnalysis`/`MLJ`/`Turing`/`CUDA` absent | integration | same, existing gate clauses (d)–(h) **plus a Phase-13 clause (i)** | ⚠ extend existing |
| **D-12 (reported gate)** | per-class AUC over {E,R,C} + full confusion matrix at pre-registered thresholds, on the frozen consts | reported script (NOT a unit test) | `julia --project=spike -t auto spike/p13/run_three_way_gate.jl` | ❌ Wave N |
| D-12 (reported) | binary-NRE continuity at fixed λ — **reported, not gated** | reported script | same runner, separate section | ❌ Wave N |
| D-13 (reported gate) | per-head reliability ECE ≤ `P13_ECE_GREEN`, with per-head AUC reported beside it (vacuous-pass guard) | reported script | same runner | ❌ Wave N |
| D-15/D-16 (reported) | α-ladder log-BF curves + crossing point `α*` | reported script | `julia --project=spike spike/p13/run_alpha_series.jl` | ❌ Wave N |

### Sampling Rate

- **Per task commit:** `julia --project=spike spike/test/runtests.jl` — must stay green and fast
  (fixture-scale everywhere; no training, no reported streams consumed).
- **Per wave merge:** the same spike suite **plus** `Pkg.test()` on the root project (proves `src/`
  is untouched and the root Manifest still resolves).
- **Phase gate:** the reported D-12/D-13 runs execute **once**, on the frozen `spike/p13/consts.jl`,
  as standalone scripts; their reports are artifacts, not test assertions. Full spike suite green
  before `/gsd:verify-work`.

### Wave 0 Gaps

- [ ] `spike/p13/consts.jl` — the D-04 pre-registration (Tier-1 constants + `_p13_forbidden()` +
      executable `@assert` self-checks), committed **before** anything runs
- [ ] `spike/p13/preconditions.jl` — the Phase-11 hard block (§D2)
- [ ] `spike/p13/net.jl` — `ThreeWayEvidenceNet`, `build_three_way_net`, `masked_two_head_bce`
- [ ] `spike/p13/labels.jl` — `three_way_label`, `head_targets`, `measure_head_log_odds`
- [ ] `spike/p13/tau_probe.jl` — the D-06 probe
- [ ] `spike/p13/alpha_series.jl` — `alpha_segregate` + invariants
- [ ] `spike/p13/result.jl` — `ThreeHypothesisColocResult <: AbstractColocResult` (§G3)
- [ ] `spike/test/test_p13.jl` — every unit testset above, `include`d from `spike/test/runtests.jl`
- [ ] Extend `spike/test/runtests.jl` resolve-risk gate with a Phase-13 clause (i): no new deps,
      `NeuralEstimators == v"0.2.1"` re-asserted
- [ ] `.planning/phases/13-.../13-SC2-AMENDMENT.md` — **before any result exists** (§H)

*No framework install is needed — Julia stdlib `Test` is already the harness.*

---

## Security Domain

`security_enforcement` is not set in `.planning/config.json`; absent ⇒ enabled. Phase 13 is a local,
offline scientific computation with no network service, no user accounts, and no untrusted input in
the usual sense, so most ASVS categories are genuinely N/A. The two that are **not** N/A are real.

### Applicable ASVS categories

| ASVS category | Applies | Standard control |
|---|---|---|
| V2 Authentication | no | No auth surface; local scripts only |
| V3 Session management | no | No sessions |
| V4 Access control | **partly** | The `sealed_holdout` seal (`corpus/manifest.jl:156-157`, `open_sealed_holdout(df; reason)` with a non-empty audit reason) is a scientific-integrity access control. Phase 13 must **not** bypass it (§J2). |
| V5 Input validation | **yes** | Dimension/shape asserts at every boundary (`pair_encode`'s `iseven` guard, `ratio_input_dim` derivation, the precondition script's field asserts). Follow `api.jl:111-120`'s explicit-`ArgumentError` style. |
| V6 Cryptography | **yes (integrity only)** | SHA-256 for any fetched corpus bytes — `corpus/` already enforces "real 64-hex digest OR the `PENDING-FETCH` sentinel, never in between" (`is_real_sha256` / `is_pending_hash`, T-08-16). **Never** use Base `hash` for a download check. |

### Known threat patterns for this stack

| Pattern | STRIDE | Standard mitigation |
|---|---|---|
| **Untrusted deserialization** — JLD2/BSON load reconstructs arbitrary types | Tampering / Elevation | Phase 13 loads only artifacts it produced, from local paths. Follow T-7-01's precedent: verify before load, and keep the deserialization surface narrow (`Flux.state` + arch metadata → rebuild, rather than deserializing the estimator object). Do **not** load a `.jld2` from an unverified source. |
| Unverified download | Tampering | `fetch_verified(url, dest, sha256)` (`corpus/fetch.jl:80`) is the only sanctioned fetch path; a `PENDING-FETCH` sentinel must never be silently accepted as a digest. |
| Path traversal in artifact names | Tampering | Build all paths with `joinpath(@__DIR__, ...)`; never interpolate user strings into a path. |
| Silent overwrite of a frozen pre-registration | Repudiation (scientific integrity) | The atomic `.tmp` → integrity-check → `mv(force=true)` write wrapper already used across the project; plus git: consts files are **byte-locked** and their sha is quoted in the report (`gate_consts_16.jl` @ 7e2318b precedent). |
| Accidental consumption of a reserved seed stream | Repudiation (data snooping) | The executable `_p13_forbidden()` assert (§K1/K3) — the project's established control. |

---

## Sources

### Primary (HIGH confidence — read directly this session)

- `C:\Users\Manuel\.julia\packages\NeuralEstimators\gFxuZ\src\Estimators\RatioEstimator.jl` —
  `:78-92` constructors, `:94-122` `_inputoutput`, **`:124` the hard-coded loss**, `:160-175` `logratio`
- `…\NeuralEstimators\gFxuZ\src\train.jl` — `:101,107,113` train dispatch, `:164-176` kwargs,
  `:223` loss resolution, **`:654` generic `_loss`**, **`:657` generic `_inputoutput`**, `:660-676` dataloader
- `…\NeuralEstimators\gFxuZ\ext\NeuralEstimatorsFluxExt.jl` — `:39-41` train-state construction,
  `:45-59` checkpointing side effects, `:62-89` `_risk` / `_train_step` (**the `loss(model(input), output)` contract**)
- `…\NeuralEstimators\gFxuZ\Project.toml:3` — `version = "0.2.1"`
- `…\Flux\hrg9M\Project.toml:3` — `version = "0.16.10"`; `…\Flux\*\src\layers\macro.jl:55` — `@layer`
- `src/amortized/train_ratio.jl` (whole file), `src/amortized/bf.jl` (whole file),
  `src/amortized/infer.jl:70-134`, `src/amortized/summary.jl` (whole file),
  `src/amortized/api.jl:110-148`, `src/amortized/datagen.jl:199-224`
- `src/results.jl:19-60, 85-192`; `src/LoadImages.jl:78-130, 215-265, 430-455`;
  `src/colocalization.jl:140-204`
- `spike/simulator/ghat.jl` (whole file), `spike/simulator/prior.jl` (whole file),
  `spike/contract.jl:50-91`, `spike/validation/consts.jl` (whole file),
  `spike/validation/sbc.jl:55-160`, `spike/test/runtests.jl` (whole file),
  `spike/Project.toml`, `spike/Manifest.toml` (Flux / NeuralEstimators / Images / ImageFiltering blocks, `julia_version`)
- `corpus/manifest.csv`, `corpus/load.jl:40-95`, `corpus/manifest.jl:149-157`, `corpus/fetch.jl:50-80`
- `.planning/spikes/014-bf-sim-validation/README.md` + `run_sim_bf.jl:1-90` — the simulation-based
  BF method, its seed discipline, and the KDE-attrition numbers
- `.claude/skills/spike-findings-proteincoloc/references/bayes-factor-gate.md`,
  `references/sbc-calibration.md`
- `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md:1-60` — the amendment
  precedent and its legitimacy test
- `docs/amortized.md` §Named limits (limits #1–#7)
- `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-CONTEXT.md` and
  `11-RESEARCH.md` §C9, §D13, §G22-G23 — the upstream contract, probe pattern, and seed inventory
- `.planning/ROADMAP.md:255-264`; `.planning/REQUIREMENTS.md` (no Phase-13 IDs); `.planning/STATE.md`

### Computed this session (HIGH confidence — reproducible)

- Prior class masses, atom asymmetry, `ghat⁻¹(0)`: exact truncated-Cauchy CDF + `ghat` knot inversion
  (`scratchpad/p13probe.jl`) and a 2×10⁶-draw Monte Carlo under `Random.seed!(20260725)`
  (`scratchpad/p13joint.jl`). Both re-runnable with `julia --startup-file=no`.
- Seed-collision check for `0x0B13DE71` / `0x2545F4914F6CDD1D`: repo-wide case-insensitive `grep`,
  zero hits.

### Secondary (MEDIUM confidence)

- The D-07 derivation (§F1–F3) is standard density-ratio / class-prior-correction theory applied to
  this codebase's exact constructions. The algebra is elementary and stated in full; the **empirical**
  verification (§F5) has not been run and is a planned Phase-13 task.
- The observation that `amortized_log_bf`'s `− log_prior_odds` is redundant follows from the verified
  `_inputoutput` construction. It is a derivation, not a measurement — flagged rather than acted on.

### Tertiary (LOW confidence — flagged)

- Relative size of the `S-BIAD1047` negative anchor (§J2 A6) — not measured.
- Any statement about Phase 11's realized artifact names, `encode_lambda` form, or net quality —
  Phase 11 has not executed.

---

## Metadata

**Confidence breakdown:**

| Area | Level | Reason |
|---|---|---|
| NeuralEstimators / Flux API surface | **HIGH** | Read from the exact installed, manifest-pinned source; line numbers quoted |
| D-10 trunk numbers, D-09 premise | **HIGH** | Every constant checked against `train_ratio.jl`; no mismatch found |
| Prior masses, atom asymmetry, `ghat⁻¹(0)` | **HIGH** | Exact CDF + 2×10⁶ Monte Carlo, both re-runnable |
| D-07 correction derivation | **HIGH (algebra) / MEDIUM (empirics)** | Derivation is complete and elementary; the numerical verification is a planned task, not a result |
| Repo facts (corpus, seals, seeds, gates, precedents) | **HIGH** | Direct file reads |
| τ probe design | **MEDIUM** | Sound and cheap, but τ itself is unmeasured; the AUC bar is a recommendation (A2) |
| α-series design | **MEDIUM** | The `_exclude_zero` trap and the bitwise-α=0 property are verified from source; the ladder's *behaviour* is unmeasured |
| Anything downstream of Phase 11 | **LOW** | Phase 11 has not executed |

**Research date:** 2026-07-25
**Valid until:** ~2026-08-24 for the API/stack facts (pinned Manifest ⇒ stable).
**Invalidated earlier by:** Phase 11 landing (re-read its artifact contract), any `spike/Manifest.toml`
change, or a user override of D-05 variant (a) / the D-07 stratification design.
