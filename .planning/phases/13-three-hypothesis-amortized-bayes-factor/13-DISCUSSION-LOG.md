# Phase 13: Three-Hypothesis Amortized Bayes Factor - Discussion Log

> **Audit trail only.** Not consumed by researcher, planner, or executor agents — decisions live in
> `13-CONTEXT.md`. This log preserves the options presented and the reasoning offered.

**Date:** 2026-07-25
**Mode:** discuss `--analyze --all` (trade-off table before each question; all gray areas auto-selected)
**Areas:** Retraining posture and lane · Hypothesis boundaries · Three-way head architecture ·
Validation target (SC2 correction) · Exclusion ground truth (SC3)

`[--all] Selected all gray areas.`

---

## Pre-discussion scouting findings

Four findings were surfaced before any question, each of which changed how a question had to be posed:

1. **Exclusion is already in the training distribution** — `GHAT_RHO_KNOTS` spans `[-0.99, +0.99]`;
   `RATIO_SPLIT_THRESHOLD = 0.0` lumps random and exclusion into one "null" class. Three-way is
   largely a relabel plus a head swap.
2. **NeuralEstimators blocks a native 3-class head** — `RatioEstimator(net, 1; …)` is binary, and per
   `train_ratio.jl:38-39` a passed loss is *silently ignored* in v0.2.1.
3. **The prior atom lands on the exclusion end** — asymmetric `ghat` clamping puts ≈4.9% of prior mass
   at `ρ = -0.99` versus ≈1.9% at `+0.99` under `Truncated(Cauchy(0,0.3),-1,1)`.
4. **ROADMAP SC2 cites a retired reference** — `compute_BayesFactor()` was established invalid past
   `|logBF| ≈ 6.9` by named limit #3.

---

## Area 1: Retraining posture and lane

### Q1. What retraining posture should Phase 13 take?

Trade-off presented: research lane vs reship-the-ratio-net vs extend-in-place. Noted asymmetry — a
reship is materially cheaper here than in Phase 11/12 because Phase 13 touches neither the NPE, the
summary, nor any `DATAGEN_HASH_SRC_FILES` entry.

- **Research lane** *(recommended)* — fresh DEV seed disjoint from the pre-registered set; shipped
  artifact, `bayes_factor` accessor and the GO untouched.
- Reship the ratio net — new pre-registered gate for the evidence net alone.
- Research lane, reship-ready — build so a later reship is a pre-registration step, not a rebuild.

**Selected: Research lane.** → D-01

### Q2. Which frozen summary basis should the 3-way net train on?

Context: `train_ratio.jl:144` trains on pairs assembled from the FROZEN `zt` of an already-trained NPE.

- Shipped grid-8 frozen `zt` *(recommended)* — honors ROADMAP deps (Phase 7 only), keeps 13 parallel
  to 11 and 12, exact same-input 2-way comparison.
- **Phase-11 registration-aware `zt`** — exclusion robust to realistic misalignment; serializes 13
  behind 11 and makes the net non-input-comparable to the shipped binary BF.
- Shipped basis plus a registration sensitivity note.

**Selected: Phase-11 registration-aware `zt`** (against the recommendation). → D-02

> **Consequence recorded, not prompted on:** this adds a Phase 11 dependency ROADMAP does not declare.
> Phase 13 planning proceeds now; execution blocks on Phase 11's research net existing. It also
> weakens the SC2 overlap check, which fed directly into Area 4.

### Q3. Should the evidence net take Phase 11's registration-uncertainty level as an input?

- **Condition on uncertainty** *(recommended)* — the log-BF pair weakens as uncertainty grows,
  matching Phase 11 D-07; avoids a widening posterior beside an unmoved Bayes factor.
- Fixed uncertainty level — simplest, closest to like-for-like SC2, but bakes in an unstated assumption.
- Marginalize over uncertainty — honest average, but hides the dependence entirely.

**Selected: Condition on uncertainty.** → D-03

### Q4 (not asked — carried forward)

Pre-registration posture was **carried forward from Phase 11 D-04** rather than re-asked, as offering
back the obvious precedent would have been a rubber-stamp: pre-register thresholds in a Phase-13 consts
file with the fresh DEV seed before any run, one documented iteration allowed. User was offered the
chance to make Phase 13 stricter. → D-04

---

## Area 2: Hypothesis boundaries

### Q1. What quantity are the three hypotheses cut on?

Correctness trap surfaced from `infer.jl:117-121`
(`Δρ = mean(ρ_sample_draws .- ρ_control_draws)`): a sample at `ρ = 0.8` under a control at `ρ = 0.9`
has `Δρ < 0`. The binary net calls that "null" (harmless); a naive three-way extension on the same
quantity would publish it as "mutually exclusive".

- **Two-factor: `ρ_sample` sign × control contrast** *(recommended)* — segregation means segregation,
  control retained.
- Δρ contrast with a band — best SC2 continuity, but "exclusion" would mean "less colocalized than
  control".
- `ρ_sample` only with a band — matches biology, abandons the control comparison.

**Selected: Two-factor.** → D-05

### Q2. How is the null-band width τ set?

- **Measured resolution** *(recommended)* — τ = smallest ρ difference the summary can resolve; same
  posture as Phase 11 D-06.
- Fixed conventional value — arbitrary, no principled defence.
- Biologically-meaningful effect size — not established, varies by system.

**Selected: Measured resolution.** → D-06

### Q3. How is the exclusion class's prior-atom concentration handled?

Noted that the fix already exists: `measure_log_prior_odds` + `amortized_log_bf(…, log_prior_odds)`
means training class balance need not match π(θ) — only cover the classes, with prior odds restored
at read time.

- **Stratify + prior-odds correction** *(recommended)* — π(θ) untouched, so the Turing
  prior-consistency constraint holds; exclusion evidence becomes graded.
- Accept the atom, name it as a limit — the net plausibly learns "exclusion = exactly −0.99".
- Truncate the μ prior (spike 013) — VALIDATED with cost: breaks π(θ)/Turing consistency, degrades
  high-|ρ| coverage.

**Selected: Stratify + prior-odds correction.** → D-07

### Q4. What is the three-hypothesis output?

Context: `src/results.jl:172` sketches `log_bf_simplex :: NTuple{3,Float64}` without saying relative
to what — and Bayes factors vs posterior model probabilities are not interchangeable.

- **Log-BFs vs the random reference** *(recommended)* — same semantics as the shipped `bayes_factor`,
  no invented model prior, leaves hypothesis selection to Phase 14.
- Normalized model probabilities — requires model priors that do not exist; pre-empts Phase 14.
- Full pairwise log-BF matrix — only two of three independent; invites cherry-picking.

**Selected: Log-BFs vs the random reference.** → D-08 (with the `log_bf_simplex` rename)

---

## Area 3: Three-way head architecture

### Q1. How is three-way evidence built past a binary-only API?

- **Shared trunk, two BCE heads** *(recommended)* — one forward pass (SC1), preserves the logit/BCE
  semantics D-07 depends on, structurally consistent evidences. Cost: outside `RatioEstimator`, no
  inherited gate lineage.
- Two independent `RatioEstimator`s — maximum reuse of validated machinery; two passes, and the two
  evidences can be mutually inconsistent.
- Single 3-class softmax — reintroduces the invented model prior D-08 rejected.

**Selected: Shared trunk, two BCE heads.** → D-09

### Q2. What trunk capacity?

- **Reuse the exact 256/64 trunk** *(recommended)* — difference attributable to the head alone;
  backed by spike 012 (PARTIAL), which found capacity redistributes error rather than removing it.
- Re-tune capacity — confounds head change with capacity change.
- Larger trunk up front — contradicts spike 012 directly.

**Selected: Reuse the exact 256/64 trunk.** → D-10

### Claude's-discretion calls (reported, not asked)

Both follow necessarily from D-08/D-09, so they were settled rather than prompted: **each head trains
only on its own class pair** (one-vs-rest would make the logit a coloc-vs-mixture ratio, silently
breaking D-08), and **the trunk trains jointly** so it sees all three classes — the mechanism that
makes the two evidences structurally consistent. User offered an override; none requested. → D-11

---

## Area 4: Validation target (SC2 correction)

Two independent problems with SC2 were laid out: (a) it cites `compute_BayesFactor()`, which named
limit #3 established is invalid past `|logBF| ≈ 6.9`; (b) the Area-1 choices (Phase-11 basis +
uncertainty conditioning) mean the input surface no longer matches the shipped binary net's, so even a
valid overlap comparison would not be like-for-like.

### Q1. What does SC2 validate against?

- **3-way discrimination as the gate + binary-NRE continuity reported** *(recommended)* — mirrors the
  resolution the project already reached for the two-way BF; puts the continuity claim against the
  shipped net rather than the retired baseline.
- Binary-NRE reproduction as the gate — strongest backward compatibility, but gating on an approximate
  correspondence invites a Phase-7-style argument about what counts as agreement.
- Keep SC2 as written — validates against a reference proven invalid in exactly the tail that matters.

**Selected: 3-way discrimination + binary-NRE continuity.** → D-12 (SC2 amended, documented up front)

### Q2. Gate on evidence calibration too?

Tension stated plainly: named limit #3 records a deliberate retreat from magnitude claims, but Phase 14
will threshold on these magnitudes — an abstention threshold *is* a magnitude threshold. Argued that the
retreat was forced by a broken *reference*, not by a measurement that magnitudes are wrong, and that
classifier calibration needs no external reference at all.

- **Discrimination + calibration** *(recommended)* — ECE/MCE via existing `_bin_calibration` /
  `CalibrationResult` traffic-light into `CalibrationMeta`; gives Phase 14 usable magnitudes.
- Discrimination only — narrowest defensible claim, but pushes the problem downstream.
- Discrimination + calibration + exclusion-recall floor — a recall floor is gameable by widening τ.

**Selected: Discrimination + calibration.** → D-13

### Claude's-discretion call (reported, not asked)

`_bin_calibration` scores against a **traffic-light band**, not a p-value, so the calibration gate
cannot be over-powered the way Phase 7's M=2000 SBC point-null was. No TOST/equivalence machinery is
needed here, unlike Phase 11 D-08. → D-14

---

## Area 5: Exclusion ground truth (SC3)

Corpus audit presented (`corpus/manifest.csv`): exactly **one** real segregated anchor
(`physical-primary | negative | segregated | 0.0 | mitochondria`), one positive coloc anchor, and 30
`simulated-secondary` CBS Red-Green rows whose degrees are colocalization *percentages* — so the CBS
`0.0` end is **random, not exclusion**. The corpus therefore spans coloc→random and never
random→exclusion. Conceptual gap also named: the simulator's "exclusion" is negative intensity
correlation; a biologist's is disjoint spatial localization.

### Q1. What are the "segregated ground-truth inputs"?

- **Simulator gate + semi-synthetic graded series + physical anchor as a named check** *(recommended)*
  — the only option that fills the missing random→exclusion axis without reopening a completed phase.
- Simulator gate + physical anchor only — nothing semi-synthetic to defend, but n=1 and the
  correlation-vs-localization gap goes untested.
- Extend the corpus — the honest fix, but Phase 8 scope and Phase 8 is COMPLETE.

**Selected: Simulator gate + semi-synthetic series + physical anchor.** → D-15

### Q2. How is the semi-synthetic series constructed?

- **Mask-based disjoint reassignment** *(recommended)* — segregation constructed *spatially* via a
  mixing parameter α from random to fully disjoint; tests the simulator's assumption rather than
  restating it.
- Displacement of ch2 objects — confounds with the mis-registration nuisance Phase 11 marginalizes
  over; a null result would be uninterpretable.
- Intensity anti-correlation — reproduces the simulator's own mechanism; circular.

**Selected: Mask-based disjoint reassignment.** → D-16

---

## Deferred ideas raised during discussion

- Reshipping the 3-way evidence net (noted as materially cheaper here than in Phase 11/12)
- Extending the corpus with real graded segregation anchors (Phase 8 scope, COMPLETE)
- μ-prior truncation to remove atoms at source (spike 013, validated-with-cost)
- Full pairwise log-BF matrix including a direct coloc-vs-exclusion contrast
- Posterior model probabilities over the three hypotheses (Phase 14 territory)
- A second 3-way net on the shipped grid-8 basis for an exact 2-way overlap comparison

All carried into `13-CONTEXT.md` `<deferred>`.
