# 13 — PRE-REGISTRATION AMENDMENT to ROADMAP Phase-13 Success Criteria

**Date:** 2026-07-27
**Status:** **FROZEN SPECIFICATION. No Phase-13 result exists.**
**Implements:** `spike/p13/consts.jl` (Tier-1, frozen)
**Supplements — does NOT replace:** `.planning/ROADMAP.md` Phase-13 SC1, SC2 and SC3
**Depends on:** `docs/amortized.md` named limit 3; `.planning/spikes/014-bf-sim-validation/README.md`;
`13-CONTEXT.md` decisions **D-02**, **D-08**, **D-09**, **D-12**, **D-13**,
**D-15** (AMENDED 2026-07-25, commit `5b4da6d`), **D-16**

---

## 0. Provenance disclosure — read this first

**This amendment inverts Phase 7's disclosure, and the difference is the whole point.**

`07-GATE-AMENDMENT.md` was written **after** the original grid-8 ship-gate results were seen; its §0
says so, and every change in it had to be argued from test-design properties precisely because the
author already knew the verdicts. **Phase 13's position is not that.** At the time of writing, no
Phase-13 number of any kind exists: `spike/p13/` contains exactly one file — the frozen Tier-1
pre-registration `consts.jl` — there is no trained three-way net, no labelled evaluation set, no
gate report artifact, and no real-image or α-ladder run. Phase 13 is **not amending a criterion
because it failed**; it is declining to adopt a reference that a *prior* phase already retired and
already published as a named limit.

Recorded so the ordering is checkable afterwards:

| Field | Value |
|---|---|
| Date written | **2026-07-27** |
| Commit sha at time of writing | **`df16e76941cefd28147efb4d4c91278597836cbd`** |
| Phase-13 results in existence at that sha | **none** — no `spike/p13/*report*.jld2`, no trained net, no artifact of any kind |
| Phase-13 executable state at that sha | `spike/p13/consts.jl` (Tier-1 pre-registration) and `spike/test/test_p13_consts.jl` only |

The test of legitimacy applied throughout is the precedent's, quoted verbatim:

> *Would this change have been made, with this justification, by someone who had seen the gate's
> code but none of its outputs?*

**Phase 13's answer is stronger than "yes".** The justification is *not a Phase-13 result at all*.
It is `docs/amortized.md` named limit 3 — already published, already carried by the 2026-07-24 GO —
establishing that the KDE Bayes-factor reference is invalid past `|logBF| ≈ log(L) ≈ 6.9` and that
its attrition is 100% baseline-side. Both defects below are readable from source code with every
artifacts directory deleted.

**Do not blur the two positions.** A reader who files this document alongside
`07-GATE-AMENDMENT.md` as "the second time they moved a goalpost" has misread it. Phase 7 moved a
bar it had already missed; Phase 13 is refusing to be scored against an instrument its own prior
work proved cannot read the range in question. Those are different acts and this project's honesty
record depends on keeping them distinct.

### 0.1 A second disclosure, specific to §5

The SC3 scope written in §5 was itself **amended** (commit `5b4da6d`, 2026-07-25) after a
**correction to an earlier audit**. The first reading of D-15 named the `corpus/` physical
mitochondria anchor as the real-data check. That reading had counted the corpus anchors correctly
but had **not read their `split` and `sha256` columns**: both `physical-primary` rows are
`sha256 = "PENDING-FETCH"`, `bytes = 0`, `split = sealed_holdout` — unfetched, and reserved for the
**Phase-16** blind evaluation. The scope was therefore rewritten to substitute the six committed
`test/test_images/` TIFFs, with reduced standing.

Recording that correction here, rather than quietly adopting the new scope, is the same discipline
the rest of this document applies to SC2. The amendment is not improved by hiding that it needed
amending.

---

## 1. Defect A — the reference is invalid in exactly the regime that matters

### The defect

ROADMAP SC2 requires that the three-way evidence network *"reproduces `compute_BayesFactor()`
(`src/bayes.jl:109`) in the overlapping 2-way regime"*. That reference reads posterior tail mass
through a kernel density estimate over **L = 999** posterior draws. The smallest tail probability
that L draws can resolve is `1/L = 1/999`, so the largest ratio of two such probabilities that is
**measured** rather than **extrapolated** is bounded by

```
|log BF| ≤ log(999) ≈ 6.907
```

Past that bound the reference is reporting the shape of a Gaussian kernel's tail, not a property of
the posterior: there are no draws out there to weigh. Beyond it the estimator does one of two
things, both recorded — it returns a **non-finite** value (zero denominator mass), or it **saturates**
against the clamp rail.

Quantified on the identical pairs by spike 014
[`.planning/spikes/014-bf-sim-validation/README.md`, VERIFIED]:

| Quantity | Value |
|---|---|
| KDE baseline non-finite (attrited by the `isfinite` filter) | **34.75% (139/400)** |
| KDE baseline saturated at the ±18.42 rail (clamped variant) | **32.5% (130/400)** |
| median `|Δρ|` of **dropped** pairs | **0.736** |
| median `|Δρ|` of **kept** pairs | **0.283** |
| Attrition incurred by the amortized NRE on the same pairs | **0** |

Attrition is therefore **100% baseline-side**, and it is **not random**: the dropped pairs have
median `|Δρ|` 0.736 against 0.283 for the kept ones. The reference fails preferentially on the
high-signal pairs. Any "agreement" statistic computed on the surviving 261/400 is computed on a
survivorship-biased subset in which the reference happened to stay finite — which is close to the
definition of a circular check. See also `references/bayes-factor-gate.md` §008 (skill
`spike-findings-proteincoloc`, spikes 006–008), which is where this ceiling was first diagnosed.

`docs/amortized.md` named limit 3 states the settled position, quoted exactly:

> **The Bayes factor is validated by simulation-based DISCRIMINATION, not magnitude agreement with
> a per-pair reference.** The KDE Bayes-factor baseline was shown **invalid past |logBF| ≈ log(L) ≈
> 6.9** (attrition is 100% baseline-side; the NRE never goes non-finite). The BF is instead
> validated SBC-style (spike 014): **AUC(coloc vs null) = 0.994**, 0 attrition, monotone in evidence
> (Spearman 0.93), decision-calibration ECE 0.019, FPR 0.036 at logBF > 0.

### Outcome-independent justification

The ceiling is a property of **L**, not of any network. It is arithmetic on the reference's own
construction: count the draws, take the reciprocal, take the log. It is computable from
`src/bayes.jl` with every artifact on disk deleted and with no Phase-13 net in existence. The
finding is furthermore **already published** — named limit 3 in `docs/amortized.md` — and was
**already carried by the 2026-07-24 GO decision**, which adopted simulation-based discrimination as
the BF validation method for the shipped binary net. Phase 13 is inheriting a prior ruling, not
minting a new one.

*Would a blind reviewer have made this change?* **Yes, and more strongly than in Phase 7.** A
reviewer who had read `src/bayes.jl` and `docs/amortized.md` but had never seen a Phase-13 output
would object to SC2 on sight, because SC2 asks a phase to reproduce a reference the same repository
documents as invalid.

### Why this is not a convenience

The invalid band is not a peripheral corner. `|log BF| > 6.9` is the **confident-call tail** — the
region where a three-way verdict actually asserts something. A coloc-versus-exclusion decision at
`|log BF| ≈ 1` is a shrug; at `|log BF| ≈ 12` it is a claim. SC2 as written would validate the new
network exactly where nothing is at stake and fall silent exactly where the phase's contribution
lives.

### The new rule

SC2's reproduction requirement is **replaced** by the criterion in §3. Correspondence with the
**shipped binary NRE** in the overlapping 2-way regime is **reported, not gated**
(`P13_CONTINUITY_GATED = false`).

---

## 2. Defect B — the input surfaces are not shared (architectural, needs no results)

### The defect

D-02 trains the Phase-13 evidence net on **Phase 11's registration-aware frozen `zt`**, not on the
shipped grid-8 `zt`. Two things change at once:

1. the standardized summary **basis** is Phase 11's, fitted on registration-perturbed training data,
   not the shipped net's; and
2. an extra **registration-uncertainty conditioning input** (`λ`) is appended per D-03, so the input
   width is `ratio_input_dim(G) + n_cond`, derived and never written as a literal
   (`P13_INPUT_WIDTH_RULE = :ratio_input_dim_plus_ncond`,
   `P13_LAMBDA_PLACEMENT = :append_after_pair_encode`).

That input surface is shared with **neither** `compute_BayesFactor()`'s Turing/ADVI path (which
consumes the raw channel intensities through a probabilistic program, and has no `zt` and no `λ` at
all) **nor** the shipped binary NRE (a different frozen basis, no conditioning input). Therefore even
a **perfectly valid** KDE reference would not give a like-for-like comparison: a disagreement would
be a statement about the basis change and the conditioning input, not about either network's evidence.

### The subtler consequence — the two logits are different objects

The shipped binary NRE's logit is a **shuffled-θ likelihood-to-evidence ratio**: NeuralEstimators'
`RatioEstimator` separates the joint `(Z, m)` from the product of marginals `(Z, m̃)` with `m̃` a
shuffle of `m` at equal counts, so the Bayes-optimal logit is `log[p(Z|m)/p(Z)]` and the difference
of the two model indices *is* `log BF(1:0)` with no correction term (13-RESEARCH §A3).

The D-09 heads are **plain BCE classifiers on class labels**. Their Bayes-optimal logit carries the
head's own training class frequencies:

```
s_h(Z) = log[ p_q(Z|A) / p_q(Z|B) ] + log(q_A / q_B)
⇒ log BF_q(A : B) = s_h(Z) − log(q_A / q_B)
```

with `q_A`, `q_B` measured on that head's **restricted** training subset (13-RESEARCH §F1). So the
D-09 heads **need a correction the shipped net does not**, and the two quantities are **not
numerically interchangeable even at identical inputs**. That, and not any expected disagreement, is
why the binary-NRE continuity check is reported and not gated.

### Outcome-independent justification

This is a statement about which arrays enter which function, and about the algebra of two training
objectives. It requires reading `src/amortized/train_ratio.jl`, `src/amortized/summary.jl`,
`src/bayes.jl` and Phase 11's output contract — and nothing else. It is decidable before a single
weight is initialized, and it would be equally true if the Phase-13 net turned out to agree with the
KDE reference to four decimal places.

*Would a blind reviewer have made this change?* **Yes.** "These two estimators do not share an input
space, so magnitude agreement between them is not a validity check" is a design-level objection, not
a result-level one.

### Explicit non-change to shipped code

This amendment does **not** propose changing `src/amortized/bf.jl`. The shipped `amortized_log_bf`
accessor subtracts `− log_prior_odds` on top of the NRE difference, which is redundant under the
shuffled-θ identity; at the shipped net's measured balance (`P(Δρ > 0) = 0.4993`) that term is
`−0.0102`, which sits **below the reporting floor**. D-01 freezes the shipped accessor and the
Phase-7 GO. The redundancy is documented here and left alone.

---

## 3. The replacement criterion, fully specified

**Binding source of every number below: `spike/p13/consts.jl` (Tier-1, frozen, committed before this
document).** Constant names are given so the criterion cannot drift by paraphrase.

### 3.1 Gated — three-way discrimination at simulator ground truth

- **Per-class one-vs-random AUC** over {exclusion, random, coloc}, at simulator ground-truth labels.
- Floors: **`P13_AUC_FLOOR_COLOC = 0.90`** and **`P13_AUC_FLOOR_EXCLUSION = 0.90`**. These are
  floors, not targets — the shipped binary analogue reached AUC 0.994 in spike 014.
- Evaluated on **`P13_GATE_M = 4000`** labelled pairs, with at least
  **`P13_MIN_EVAL_PER_HEAD = 1000`** samples in each head's own restricted evaluation set (each head
  is scored only on its own class pair, per D-11).

### 3.2 Gated — per-head evidence calibration (D-13)

- **`ECE ≤ P13_ECE_GREEN = 0.05`** per head, computed with `P13_ECE_NBINS = 10` through the existing
  `_bin_calibration` / `CalibrationResult` traffic light. `P13_GATE_STATISTIC = :ece`.
- **MCE is reported, not gated**, and is reported **beside its empty-bin count**. Reason:
  `_bin_calibration` scores an empty bin as `predicted_rate = midpoint`, `observed_rate = 0.0`,
  which contributes **zero** ECE weight (ECE is count-weighted) but **full** MCE weight (MCE is an
  unweighted max over bins). On a well-separated three-way problem many of the ten bins will be
  empty, so MCE would be dominated by empty bins. `_bin_calibration` is shared, gate-lineage code
  and is **not** to be hand-modified to work around this.
- **Vacuous-pass guard:** each head's **AUC is reported beside its ECE**, and a green ECE on a head
  with AUC below **`P13_VACUOUS_AUC_FLOOR = 0.60`** must be **labelled a vacuous pass** in the
  report and never quoted as calibration evidence. A head that says "0.5, always" is perfectly
  calibrated and perfectly useless.
- No equivalence/TOST machinery is required (`P13_TOST_REQUIRED = false`): the traffic light is a
  **band**, not a p-value, so it cannot be over-powered the way Phase 7's M = 2000 SBC point-null was.

### 3.3 Reported, not gated

- The full **3×3 confusion matrix** at `argmax(logBF_C, 0, logBF_E)`
  (`P13_CONFUSION_RULE = :argmax_descriptive`), explicitly labelled **descriptive, not a decision
  rule**. Decisions and abstention are **Phase 14's** scope; a descriptive argmax quoted as a
  classifier would pre-empt that phase's design.
- **Binary-NRE continuity** in the overlapping 2-way regime at the reference λ
  (`P13_CONTINUITY_GATED = false`) — see §2 for why a disagreement would not be a failure.
- The **λ response** of the evidence pair (evidence must not stay confident as registration
  uncertainty grows).
- The **D-16 α-graded semi-synthetic series** and its crossing point α\* (`P13_ALPHA_GATED = false`).
- **Every real-image quantity, without exception.** See §5.

### 3.4 Where the labels come from

Labels are produced by the **D-05 two-factor cut**: the `ρ_sample` **level** (sign against τ)
**times** the contrast against the control. They are **never** taken from `Δρ`.

The counter-example that makes this non-negotiable: `infer.jl:117-121` defines
`Δρ = mean(ρ_sample_draws .- ρ_control_draws)`. A sample at `ρ = 0.8` under a control at `ρ = 0.9`
has `Δρ < 0`. The shipped binary net calls that "null", which is harmless. A naive three-way
extension cut on the same quantity would publish it as **"mutually exclusive"** while the sample is
strongly colocalized. **"Less colocalized than the control" is not segregation.**

---

## 4. What SC1's superseded wording becomes

ROADMAP SC1 reads: *"A 3-way `RatioEstimator`/evidence network emits a log-BF simplex over
{coloc, random, exclusion} in one forward pass."* Two of its terms are not achievable and are not
desirable:

- **`RatioEstimator` is ruled out by D-09.** `RatioEstimator(net, 1; num_summaries)` takes
  `num_parameters = 1` — a *binary* model index — and its training loss is **hard-coded**
  `logitbinarycrossentropy`; **a passed loss is silently ignored in v0.2.1**. A three-way head is not
  reachable by passing a loss. This is an architectural fork, not a tuning knob.
- **"log-BF simplex" is ruled out by D-08.** The output is `log BF(coloc : random)` and
  `log BF(exclusion : random)` with the random entry **identically 0** — two free numbers plus a
  structural zero. That is not a point on a simplex, and a field named "simplex" would be quoted as
  one.

**The surviving, bindable criterion:** *one network, one forward pass, three-way evidence emitted as
`log BF(coloc : random)` and `log BF(exclusion : random)`, with the random entry identically 0.*

One shared Flux trunk (reused verbatim per D-10) feeding two BCE heads satisfies "one forward pass"
and preserves the logit/BCE semantics that make the §3.4 / §2 correction valid; a softmax head would
have broken them.

---

## 5. SC3 scope (AMENDED — three arms, one of which gates)

This section **replaces the pre-amendment single-sentence deferral entirely.** SC3 reads *"The
exclusion hypothesis is validated on segregated ground-truth inputs."* It is **retained** and
**scoped** into three arms of deliberately unequal standing:

| Arm | Substrate | Standing |
|---|---|---|
| **Simulator ground truth** | `ρ_true < −τ`, exact labels, well-powered, `P13_GATE_M = 4000` | **GATES** — the D-12/D-13 floors of §3 |
| **D-16 α-graded series** | mask-based disjoint reassignment on simulated **and** real pairs | Supporting evidence — **not a gate** (`P13_ALPHA_GATED = false`) |
| **`test/test_images/` qualitative check** | six committed 1028 × 1376 microscopy TIFFs | Required deliverable — **explicitly not a gate** |

The ordering is stated **here, before any curve exists**, because the arm that will attract the most
reader attention (real microscopy) is the weakest one.

### (a) Why the α series is load-bearing

The corpus holds **one** coloc anchor (`pos-tetraspeck-01`), **one** segregated anchor
(`neg-lightmycells-01`), and **30** `simulated-secondary` CBS Red-Green rows whose degrees are
colocalization **percentages** — so the CBS `0.0` end is **random, not exclusion**. **There is no
graded segregation series anywhere in the corpus.** The D-16 mask-based construction is therefore the
**only** graded segregation evidence that will exist for this phase, which is why it is load-bearing
supporting evidence rather than a nice-to-have.

It is also the only arm that probes the **conceptual gap**: the simulator's "exclusion" is *negative
intensity correlation across patches*, while a biologist's "mutually exclusive" is *disjoint spatial
localization*. D-16 constructs segregation **spatially**, so the series tests the simulator's
assumption instead of restating it.

### (b) Why the corpus physical anchor is not used

Both `physical-primary` rows in `corpus/manifest.csv` are `sha256 = "PENDING-FETCH"`, `bytes = 0`,
and `split = sealed_holdout` — **unfetched**, not on disk, and reserved for the **Phase-16 blind
evaluation**. Consuming the segregated anchor here would irreversibly burn Phase 16 on the very
hypothesis Phase 16 exists to evaluate.

**The seal-break escape hatch an earlier draft floated is WITHDRAWN.** Concretely, for the whole of
Phase 13:

- there is **no `checkpoint:human-verify` bypass** and no flag that opens the seal;
- **no fetch** is performed — `bootstrap_anchor_hashes()` is not run;
- **`open_sealed_holdout` is not called** in any Phase-13 script, test, or flag path — every
  appearance of that name in Phase-13 material, including this document, is inside prose stating it is
  not called;
- no Phase-13 **executable** artifact references `corpus/` at all. `corpus/manifest.csv` may be cited
  in **prose** (this document, plans, the report) for provenance only.

This mirrors the ruling Phase 11's planner made for its D-17 and the amended Phase 12 D-11.

### (c) What the real-image arm is, and what it is not

**The six `test/test_images/` TIFFs carry no colocalization ground-truth label.** There is no
recorded ρ, no recorded Δρ, no segregation degree, and no provenance metadata asserting one.

The folder names `positive/` and `negative/` refer to the original package's **biological** test
conditions — they are **not** colocalization labels — and the measured summary confirms this
directly: the "negative" pair has **positive** mean patch correlation **`m̄ = +0.2481`**
(`ρ_true ≈ +0.215`). It is **not** an anti-correlated pair and **is not an exclusion example**.
Nothing in this phase may treat `negative/` as one.

**Therefore the real-data check is QUALITATIVE ONLY.** It can show that the three-way exclusion
verdicts behave *sensibly* on real microscopy — that ingestion works, that the evidence responds
monotonically to a constructed α-ladder, that the sign of the verdict flips where the summary's sign
flips, that nothing degenerates. **It cannot show the verdicts are correct**, because correctness
requires a label the data does not have.

**Real segregation validation remains deferred to Phase 16's blind evaluation** — which is precisely
what the sealed holdout exists for, and precisely why Phase 13 does not open it.

### (d) Two findings that make the honesty statement stronger than "no ground truth"

Both must be **reported**.

**(i) The shipped OOD detector flags both fixtures.** Mahalanobis summary-density score **433.69**
against the frozen in-distribution threshold **179.14** — **2.4× over**, `is_ood = true` for both
read directions. Real microscopy at this size and these statistics sits outside the simulator's
training distribution. That is **not a bug and not a reason to skip the check**; it is the
misspecification channel doing its job, and it is arguably the most honest single number the arm will
produce. *(Prior measurement of the **shipped** 8×8 bundle, 13-RESEARCH §J4.6 — not a Phase-13
result. The Phase-13 net is a different net and may score differently; measuring and reporting that
comparison is itself part of the arm.)*

**Standing requirement:** **every real-image number is reported beside its OOD verdict and its λ** —
in every printed table, every JLD2 field, and every figure caption.

**(ii) The frozen simulator calibration already records that real anti-correlation was never
observed.** `spike/simulator/ghat.jl:40` states:

> `negative tail : physically reachable in real fluorescence = false ⇒ negative μ-prior tail
> documented PRIOR-ONLY; consistency scoped to realized range`

This is a pre-existing, frozen Phase-02 finding, and it lands directly on the hypothesis Phase 13
promotes to first-class: **the exclusion hypothesis has no observed real-data instance anywhere in
this project.** Every unmodified real image the project has ever measured sits at **positive** induced
μ — even the "negative" control, at `+0.2481`.

This is not a reason to weaken the phase; it is the reason the α series matters. The D-16
construction produces `m̄ = −0.168` / `−0.252` (`ρ_true ≈ −0.35 / −0.44`) **from real microscopy
pixels** — the first negative induced μ this project has obtained from real data. It does **not**
refute the Phase-02 verdict; it **sharpens** it: **negative μ is constructible from real images but
has not been observed to occur naturally in the images this project holds.** *(Prior measurement of
the frozen simulator calibration and the committed fixtures, 13-RESEARCH §J5 — not a Phase-13
result.)*

### (e) A real image has no known λ

A real image's true registration error is **unmeasured and unmeasurable from the file alone**, so
there is no correct value to plug in. The pre-registered choice
(`P13_REAL_LAMBDA_READS_RULE = :full_phase11_ladder`,
`P13_REAL_LAMBDA_HEADLINE_RULE = :widest_rung`,
`P13_REAL_LAMBDA_HEADLINE_EXPECTED = 3.0`) is to **report the whole sweep** across the Phase-11
ladder and **headline the widest, most conservative rung** — because a conservative λ makes the
evidence *weaker*, so an exclusion verdict that survives it is the credible one. **A real image is
never read at an invented single λ**, and λ is never estimated from the images (that would be a new
estimator with its own validation burden, inside a phase that is not about registration).

### (f) The scope limit of the α ladder

The measured real ladder spans `ρ_true ∈ [−0.44, +0.33]` — a **densely-covered interior region** of
`GHAT_RHO_KNOTS`. It therefore probes the **sign transition and the near-exclusion regime**, and it
**does not probe the deep-exclusion tail** near the ±0.99 clamp atoms, where the prior-atom asymmetry
(≈4.9% of prior mass at the negative clamp vs ≈1.9% at the positive one, about **2.5 : 1 against the
exclusion end**) lands hardest. A reader who assumes otherwise will over-read the result.

Note also the substrate asymmetry, which must be handled in reporting and never averaged away: on
simulated substrate the ladder runs **random → exclusion**; on real substrate it runs **moderately
colocalized → near-exclusion**, because the unmodified real pairs measure `m̄ = +0.3292 / +0.2481`.
Both are informative; they are **not the same experiment**.

### (g) The declarative close

**No pass/fail threshold is defined for any real-image quantity.**

Stated in exactly those words so it is greppable rather than merely intended. If no threshold exists,
the arm cannot accidentally become a gate, and no Phase-13 verdict can turn on it. The two documented
exemptions are **ingestion regressions, not bars on a result**:
`P13_REAL_ANCHOR_MBAR = (positive = 0.3292, negative = 0.2481)` with
`P13_REAL_ANCHOR_TOL = 1e-3` checks that the *loader reproduces the frozen lineage anchors*, and the
`P13_ALPHA_*` invariant bounds check that the *transform did what it says*. Neither scores an
evidence number.

---

## 6. The original is not erased

**ROADMAP's SC1, SC2 and SC3 remain byte-present and citable.** In particular SC2 keeps its literal
wording — *"reproduces `compute_BayesFactor()` (`src/bayes.jl:109`) in the overlapping 2-way regime
without quadgk/KDE"* — in `.planning/ROADMAP.md`. This amendment adds a **pointer**, not a rewrite.

The amended criteria are a **second, separately named** pre-registration (this document plus
`spike/p13/consts.jl`), scored on a **fresh, disjoint DEV seed** (`P13_DEV_SEED`, asserted absent
from the whole forbidden-seed inventory).

**Any report that cites an amended criterion must cite the original alongside it.** Reporting only
the amended criterion would be exactly the thing this document exists to prevent. The phase report
(plan 13-14) is required to print the original SC1/SC2/SC3 text beside the amended criteria, so a
reader can see what was changed and judge the change for themselves.

---

## 7. Iteration allowance

**`P13_ITERATION_ALLOWANCE = 1`**, authorised **here**, before any result exists.

**The one pre-declared trigger.** The D-12 evaluation shows the **exclusion-class per-class AUC below
`P13_AUC_FLOOR_EXCLUSION` specifically at high `|ρ|`** (the deep-exclusion tail near the −0.99 atom)
while the mid-range is resolved. The allowance is then spent on switching to
`P13_STRATIFICATION_FALLBACK = :importance_weighted` (D-07-ii: within-range ρ-stratification with
importance-weighted BCE, with the `ρ = −0.99` atom as its own stratum) and **re-running once**.

**It is never spent on relaxing a threshold after training.** A **second** amendment is **not
authorised by this document**, and cannot be authorised by amending this document. Phase 7 was
amended twice after seeing results; the credibility cost of that is the reason this clause exists,
and named limit 1 in `docs/amortized.md` is the receipt.

**The allowance is scoped to the GATE.** It can **never** be spent in response to a real-image
observation, because no real-image quantity has a bar to miss (§5(g)). An α-ladder curve or a real
TIFF verdict cannot trigger a re-run under this document.

---

## Appendix — provenance of every number quoted above

**No measured Phase-13 number appears anywhere in this file.** Every number quoted is a prior
measurement, labelled here so the FROZEN status is unambiguous:

| Number | Source | Kind |
|---|---|---|
| `|logBF| ≤ log(999) ≈ 6.9` | `src/bayes.jl` L = 999, arithmetic | property of the reference's code |
| 34.75% (139/400), 32.5% (130/400), 0.736 / 0.283, AUC 0.994 | `.planning/spikes/014-bf-sim-validation/README.md` | prior spike, **shipped binary net** |
| `−0.0102`, `P(Δρ > 0) = 0.4993` | 13-RESEARCH §A3 Monte Carlo | prior measurement, **shipped net** |
| `m̄ = +0.3292 / +0.2481`, ladder `ρ_true ∈ [−0.44, +0.33]` | 13-RESEARCH §J5, executed 2026-07-25 against the committed fixtures | prior measurement, **frozen simulator calibration** |
| OOD density 433.69 vs ID threshold 179.14 | 13-RESEARCH §J4.6 | prior measurement, **shipped 8×8 bundle** |
| ≈4.9% / ≈1.9% prior clamp mass (≈2.5 : 1) | 13-RESEARCH §E2, `spike/simulator/ghat.jl` + `prior.jl` | property of the frozen prior |
| every `P13_*` constant | `spike/p13/consts.jl` (Tier-1, frozen) | pre-registered choice, not a measurement |

`P13_TAU` is deliberately **absent** from the Tier-1 pre-registration: D-06 requires it be
**measured**, and it is appended as Tier-2 after the probe run. Its absence at this sha is a further
check that this document precedes the phase's results.
