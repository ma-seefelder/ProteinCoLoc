# 12 — PRE-REGISTRATION AMENDMENT to ROADMAP Phase-12 Success Criteria

**Date:** 2026-07-29
**Status:** **FROZEN SPECIFICATION. No Phase-12 result exists.**
**Implements:** `spike/validation/p12_consts.jl` (Tier-1, frozen)
**Supplements — does NOT replace:** `.planning/ROADMAP.md` Phase-12 SC1, SC2 and SC3
**Depends on:** `12-CONTEXT.md` decisions **D-02**, **D-05**, **D-06**, **D-07**, **D-08**, **D-09**,
**D-10**, **D-11** (AMENDED 2026-07-25), **D-12**, **D-13**; orchestrator resolutions **R-4**,
**R-7**, **R-8**; the 2026-07-27 premise audit **V-1 … V-3**, **K-1**, **S-4**;
`12-CONSTANTS-FOR-CONFIRMATION.md` (CLOSED 2026-07-29)

---

## 0. Provenance disclosure — read this first

**This document is written before Phase 12 has produced a single number.** Three of the ROADMAP's
Phase-12 success criteria are defective as *written*, and all three defects are readable from source
code and from decisions already recorded — none of them from a result. Correcting them after a result
existed would be indistinguishable from tuning a criterion to an outcome, which is the failure this
project has already paid for once (`07-GATE-AMENDMENT.md`, amended twice **after** the grid-8 ship-gate
verdicts were seen, and carried as named limit 1 in `docs/amortized.md`).

Recorded so the ordering is checkable afterwards:

| Field | Value |
|---|---|
| Date written | **2026-07-29** |
| Commit sha at time of writing | **`af8e1ca905438ebfbd1f5fb96535c95947b6a533`** |
| Phase-12 results in existence at that sha | **none** — verified, see below |
| Phase-12 executable state at that sha | `spike/validation/p12_consts.jl` (Tier-1 pre-registration) plus 11 `spike/test/test_p12_*.jl` scaffolding files — and **nothing else**: no `spike/simulator/p12_*`, no `spike/npe/p12_*`, no `spike/data/p12_*`, no `spike/p12/`, no `run_p12_*` runner |

**The "none" row was verified, not asserted.** At the sha above, all four of the following are empty:

| Checked | Result at `af8e1ca9` |
|---|---|
| `spike/validation/p12_*report*.jld2` | absent — no such file |
| `spike/data/cache/p12/**` | absent — `spike/data/cache/` holds `fixture`, `p11`, `p13` and one content-addressed pool, and no `p12` directory |
| `spike/npe/p12_*.jld2` | absent — no such file |
| any trained Phase-12 net under any name | absent — `spike/npe/` holds only `trained_npe.jld2` (shipped lane) and `p11_research_npe.jld2` (Phase 11) |

No Phase-12 network exists, no Phase-12 pool exists, and no Phase-12 artifact of any kind exists.
**No Phase-12 result exists.**

The test of legitimacy applied throughout is the precedent's (`13-SC2-AMENDMENT.md` §0), quoted
verbatim:

> *Would this change have been made, with this justification, by someone who had seen the gate's
> code but none of its outputs?*

For every section below the answer is **yes**, and the reason is the same in each case: the
justification is a property of the *criterion's own wording* set against code and decisions that
already exist, never a Phase-12 measurement.

**One disclosure the reader is owed.** §5 records a Δρ semantics decision taken on 2026-07-28, and
§5 also **withdraws** — rather than deletes — an earlier draft clause of that decision which
contradicted the clause immediately preceding it. The retraction is written into §5 for the same
reason this whole document exists: an amendment is not improved by hiding that it needed amending.

---

## 1. SC1 correction — the DeepSet is dropped; the summary net is a CNN

ROADMAP SC1 reads *"A lattice prior (CAR vs. AbstractGPs, chosen by mini-spike) is placed over the
`correlation()` grid with a **CNN/DeepSet** summary."*

**The DeepSet half is incoherent with the phase's own goal and is dropped. Only the CNN is
retained.** A DeepSet is *permutation-invariant over patches* — which **is** the exchangeable pooling
this phase's goal sentence says it replaces. A criterion that permits the very structure the phase
exists to remove could be satisfied by a model that fails the phase. This is a contradiction between
SC1 and the ROADMAP goal line, not a preference between architectures.

The amended clause: *a lattice prior is placed over the `correlation()` grid and the summary is read
by a **CNN** over the `G×G×2` reshape.*

**D-02 is recorded here so the scope of the change is unmistakable: summary CONTENT is unchanged.**
The summary is `2·G² = 128` dims at G = 8 (64 continuous rows + 64 mask rows, one entry per patch), so
per-region information is **already present**; flattening discards the spatial *inductive bias*, not
the spatial *information*. Only the network's **view** of the same 128 numbers changes. Therefore
`patch_summary`, the OOD summary-density Mahalanobis reference, `local_coloc_map`, the shipped
content-hashed artifact and `DATAGEN_HASH_SRC_FILES` all stay valid and untouched. What does **not**
transfer for free is comparability against the shipped MLP conditioner, because the input contract
differs; any such comparison must say so.

**R-7 is binding on the CNN and is not a detail.** The topology is **lean** and is written literally
into the plans: a 1×1 channel bottleneck (`2→16→32`, `1×1→8`, flatten 512, `Dense(512,128)`), measured
at ~1.8 h at D = 72 / 50k, against ~12.5 h for a fat `2→32→64→64` / `Dense(4096,256)` variant doing the
same job. A criterion that says only "a CNN" silently selects the 12.5 h version.

*Outcome-independent?* Yes. "Permutation-invariance over patches is the thing being replaced" is a
statement about the goal sentence and the architecture's definition. It is decidable with every
artifact directory deleted.

---

## 2. SC1 correction — the two prior arms, named concretely (`AbstractGPs.jl` is NOT used)

SC1 names the comparison as *"CAR vs. AbstractGPs"*. **`AbstractGPs.jl` is not used, and no new
package is added.** Both arms are implemented as **dense 64×64 covariance/precision algebra using
stdlib `LinearAlgebra`** at G = 8 — a `G² × G²` matrix is 4096 entries, which is small enough that a GP
library buys nothing.

**The reason is a hard constraint, not a taste call.** Adding the package forces a `Pkg.resolve()` on
the frozen `spike/Manifest.toml`, which could dislodge the **NeuralEstimators 0.2.1** pin the entire
milestone rests on. The Phase-7 Wave-0 co-resolution gate already failed once on exactly this class of
transitive conflict. Stdlib `using` needs no `Project.toml` entry (Julia's default `LOAD_PATH` carries
`@stdlib`), and `spike/test/runtests.jl` asserts the **exact** 16-name dependency set, so an addition
fails the suite by name.

The amended clause: *a lattice prior — **CAR or a dense GP kernel**, chosen by the pre-registered
mini-spike — is placed over the `correlation()` grid; neither arm introduces a package dependency.*

**R-8 is recorded with it, because it is what makes the comparison mean anything.** BOTH arms are
parametrized by **induced lag-1 correlation r₁**, not by the CAR α or the GP length-scale ℓ. Measured,
lag-1 correlation runs 0.136 → 0.947 across α ∈ [0.5, 0.999], so a uniform prior on α would put
roughly 90 % of its mass on "no spatial structure" and would make D-08's correlation length
unidentifiable **by parametrization rather than by physics** — the S-1 prior-echo trap this milestone
already met once in Phase 11. Parametrizing both arms by r₁ makes the mini-spike a question about the
kernels' **shape** rather than about two incomparable coordinate systems. `P12_R1_PRIOR =
Uniform(0.05, 0.95)` is Tier-1 and is bound once, so no plan can construct its own.

R-8 also carries a correction D-05 as written does not mention and which this document therefore
states: CAR marginal variances are **not** uniform (measured sd 0.658 at an edge against 0.992 in the
interior at α = 0.95), so each cell must be divided by `sqrt(diag(Σ))` before Φ. Without that rescale
the per-region SIM-02 claim is simply false at the lattice boundary.

*Outcome-independent?* Yes. It is a dependency-resolution constraint plus a measured property of the
kernels' parametrization, both established before any network exists.

---

## 3. SC2 scoping — research lane, a spike-local `coloc_map`

ROADMAP SC2 reads *"`coloc_map(...)` returns a Δρ map + uncertainty map, amortized in a forward
pass."* Read literally that is a **public-API deliverable in `src/`**, inside a phase this project has
put in the research lane. It is **scoped**, not weakened: SC2 is satisfied by a **spike-local**
`coloc_map` returning a spike-local `SpatialColocResult`-shaped struct, with the comment-only sketch at
`src/results.jl:180-191` left **byte-unchanged** as the sketch it is.

Three reasons, all pre-existing:

1. **The 2026-07-24 GO locks "no further training, no gate iteration"** on the shipped lane. A
   Phase-12 net is a new net; shipping it through the public API is a gate action the GO explicitly
   closed.
2. **A Phase-12 net is not drop-in comparable to the shipped read surface.** The shipped path is a
   128-row summary vector → MLP → 7-marginal flow. Phase 12 is a 4-D `(G,G,2,K)` array → CNN →
   72-marginal flow. Those are different input and output contracts; a `coloc_map` in `src/` would
   have to reconcile them, which is a productionization phase's job and is listed as deferred in
   `12-CONTEXT.md`.
3. **The shipped bundle is content-hashed and Release-hosted.** Reshipping is a *release* action, not
   a phase action.

This follows the **Phase-13 13-08 precedent** exactly: the executable result subtype lives in
`spike/`, reaching its supertype `AbstractColocResult` by a **read-only** include of `src/results.jl`,
so Phase 7 D-02 (variants are new subtypes, never bolted-on fields) is satisfied while `src/` stays
provably untouched — the CLAUDE.md hard constraint.

The amended clause: *a **spike-local** `coloc_map(...)` returns a per-region Δρ map and a per-region
uncertainty map from **one** amortized forward pass, as a new `AbstractColocResult` subtype defined in
`spike/`.* "Amortized in a forward pass" is retained unchanged — it is the phase's actual claim.

*Outcome-independent?* Yes. It is a lane decision plus two contract mismatches, all readable today.

---

## 4. SC3 amendment — leave-region-out predictive coverage against a matched retrained ablation

ROADMAP SC3 reads *"The spatial (CAR) model beats independent pooling in coverage on ≥1 real image."*
**As written it is not evaluable**, for two independent reasons:

- there is **no ground-truth ρ field on a real image**, so "coverage" in the parameter sense has
  nothing to be coverage *of*; and
- the named baseline has **no uncertainty to compare**. `src/amortized/local_map.jl:55-79` states in
  its own docstring that `LocalColocMap` *"carries no posterior draws, no Bayes factor and no
  per-region uncertainty"*. A coverage comparison against an object with no intervals is undefined.

SC3 is **amended** to: *leave-region-out **PREDICTIVE** coverage (D-09) of the spatial model against a
**MATCHED RETRAINED ABLATION** (D-10, R-4), scored on real microscopy TIFFs (D-11 AMENDED), on both
calibration and a proper scoring rule.* The following eight clauses are stated **in advance** and are
part of the criterion.

### (a) Predictive coverage is not parameter coverage

The claim is *"the model predicts held-out regions honestly"*, and **no ground truth is required or
claimed anywhere in it**. The construction masks region *r* in an image's summary, predicts region *r*
from its neighbours and the shared terms, and checks the predictive interval against **that image's own
observed entry** at *r*. Nothing outside the image is consulted. A reader must not restate this as
"the recovered ρ field is correct" — that is a different claim, it is not made here, and on real data
it is not makeable.

### (b) The baseline is a retrained ablation, not a read-time pin

The "independent pooling" arm is **the same Phase-12 network, same training procedure, same summary,
RETRAINED on an r₁ → 0 prior** (`P12_ABLATION_R1 = 0.0`, Tier-1, asserted to lie outside the prior
bracket by design). Per R-4: a read-time pin of the correlation length would be an
**out-of-distribution read**, not a matched ablation, and it would confound the very gate D-10 exists to
make attributable. It is explicitly **not** `LocalColocMap`, for the reason in (a)'s second bullet.
A consequence worth stating: this makes D-13's descope deliverable a genuinely trained model rather
than a mis-specified read of another one.

### (c) The ablation is NOT a prior-only baseline

With r₁ → 0 the regions become **conditionally independent**, but the other 63 regions still inform the
shared nuisances — especially `chromatic_eps`, whose effect is global — and the global term `c₀`.
The correct description, and the words that must be used in every report:
**"no spatial borrowing, full nuisance and global borrowing"**.

Because a reader will otherwise conflate the two contributions, a **third, genuinely prior-only floor
arm** is reported alongside, so the gap between "prior only" and "no spatial borrowing" and "spatial
borrowing" can be read off directly. Only the second of those three is the gate baseline.

### (d) The criterion is TWO-SIDED, and a coverage-only win is not a pass

Both of the following are required:

1. **Calibration** — empirical coverage ≈ `P12_COVERAGE_NOMINAL = 0.90`, required of **BOTH** arms. A
   miscalibrated arm is **DISQUALIFIED**, never "better". An arm with too-wide intervals can beat a
   calibrated one on nominal coverage while being strictly less useful.
2. **Sharpness given calibration** — a **proper scoring rule**, with the spatial arm required to
   improve the mean log predictive density over the matched ablation by at least
   **`P12_STAGE2_LOGSCORE_MIN` = 0.02 nats per held-out region**.

**Pre-declared here, before any curve exists: a `coverage-only win` with no scoring-rule improvement
is NOT a pass.** Stated in those words so it is greppable rather than merely intended.

### (e) The Fisher-z observation-noise model is a declared modelling assumption

The predictive check compares an interval against an **observed patch correlation**, which is a noisy
measurement, so it needs an observation-noise model. The pre-registered choice is Fisher-z with a
single calibrated `n_eff` (`P12_FISHERZ_NEFF`, a Tier-2 append with provenance). **D-09 therefore
validates the JOINT of (posterior + noise model), not the posterior alone.** A failure is attributable
to either component, and a pass is a statement about the pair. This is declared now because it cannot
be honestly declared afterwards.

### (f) The data are the committed TIFFs, and the independent unit is the IMAGE

The real-image data are `test/test_images/{positive,negative}/*_c{1,2,3}.tif` — **six FILES, two
DATASETS** (1028 × 1376, genuine microscopy, already exercised by the root suite).

They are **NOT** the `corpus/` physical anchors. Both `physical-primary` rows in `corpus/manifest.csv`
carry `sha256 = PENDING-FETCH`, `bytes = 0` and `split = sealed_holdout`; they are unfetched and
reserved for the **Phase-16 blind evaluation**, and scoring them here would **irreversibly burn Phase
16** on the question Phase 16 exists to answer. Phase 12 consumes nothing from `corpus/`; the manifest
is cited in prose for provenance only.

**The honest n, stated in advance: 2 images × ~64 regions, and the INDEPENDENT UNIT IS THE IMAGE.**
The effective independent n is therefore **2**, and the 128 region-draws are **pseudo-replicates**
(see §5). `12-01-PLAN.md:198-204` pre-registers exactly this unit, in these words: *"the 64 regions of
one image share every nuisance and the global term, so counting region-draws as independent would
overstate the power of the gate by up to 8× in sd terms."* Any interval on a real-arm coverage number
comes from the **between-image spread**, never from a binomial over 128.

**Do NOT write "the effective sample is regions."** That phrasing contradicts the pre-registered unit
above, and it is the same unit-mismatch this phase has already had to correct twice.

### (g) The S-4 chromatic-radial confound is named in advance

`CHROMATIC_PRIOR = Uniform(-0.02, 0.02)` applies a **radial** magnification difference to channel 2 in
stage 6 of the simulator. The effective per-patch correlation therefore varies **radially in every
training draw, with no biological field present**. A spatial map can score well by learning that
radial chromatic signature — structurally the same failure class as D-06's grid-alignment artifact.

**This is named here, before any result, as the leading alternative explanation of any SC3 win**, with
its three guards:

1. a **radial-energy report** of the recovered deviation field;
2. an **ε = 0 evaluation ablation** (performance must survive `chromatic_eps` held at zero);
3. a **radially-orthogonalized field arm** (the advantage measured with the radial component projected
   out).

All four Phase-12 guards are **reporting-only**; Guard 4 in particular is **purely descriptive for
v2.0** — it reports the retained fraction of the headline log-score advantage with **no pass/fail
reading**, because no repo artifact records how much of a per-region log-score advantage is radial
(which is precisely why the guard exists), and because the orthogonal projection removes *real* radial
signal along with the artifact, so a strict fraction would penalise a truthful model with genuinely
radial biology. The judgement is carried instead by a **written attribution argument** in the Stage-2
verdict: crossing a reporting bar without that argument makes the verdict ITERATE or DESCOPE.

### (h) The real-image arm is REPORTED, not gated

**The real-image arm is `reported, not gated`, and the Stage-2 gate rests on the SIMULATED arm.**

The derivation, recorded so it cannot be re-argued afterwards:

- **`P12_STAGE2_N_MIN = 271` is *sized* with the Wald normal-approximation half-width** for 271
  *independent* trials: `z·√(p(1−p)/N) = 1.644854·√(0.09)/√N = 0.49346/√N ≤ 0.03 ⇒ √N ≥ 16.449 ⇒
  N ≥ 270.6 ⇒ N_min = 271`.
- **The interval actually reported and tested is the Wilson score interval**
  (`spike/validation/p11_stats.jl:49-64`), whose own docstring records that Wald is deliberately **not**
  used near p = 0.90 with n in the low hundreds, because its actual coverage there is poor and it can
  extend past 1.
- **Both labels are stated because an earlier draft said "Wilson-derived".** That names the *reporting*
  construction for the *sizing* arithmetic. Wilson is a different construction and yields a different
  N, so the phrase is wrong in a way a referee who knows the difference will catch — and this document
  is frozen and never edited, so the label is fixed here, before it freezes. Equally, the phase is not
  "Wald-based"; that would contradict `p11_stats.jl`. The number 271 is correct for what the formula
  computes.
- **`12-01-PLAN.md:198-204` fixes the independent unit as the DATASET**, not the region-draw, "because
  the 64 regions of one image share every nuisance and the global term."
- **D-11 supplies two images.** Applying `P12_STAGE2_COVERAGE_TOST_DELTA = 0.03` to an arm with
  effective independent n = 2 would be a criterion decided by noise — a coin flip presented as a
  pre-registration.

Therefore: **the simulated arm (12-16) is where N ≥ 271 is actually achieved, and it is what gates.**
The real arm establishes that the method runs on genuine microscopy and predicts held-out regions
sensibly — which is what D-11 meant by *"this costs Phase 12 very little"*, since D-09 had already made
the criterion ground-truth-free.

**A reader must not read the real-image numbers as a gate outcome.** They are evidence that the method
transfers, at n = 2, and nothing stronger.

---

## 5. Δρ semantics — all three maps ship, and Δρ is the primary deliverable

**Decision taken 2026-07-28; recorded here as the frozen form.** In this milestone Δρ means *sample
minus control*: `src/amortized/infer.jl:107-121` defines
`delta_rho(est, Zsample, Zcontrol, θzt)` as the Monte-Carlo difference of **two independent
single-stack passes** (D-03), `src/amortized/local_map.jl:48` documents `LOCAL_MAP_SENTINEL` as *"no
evidence of a sample-vs-control difference"*, and `local_map.jl:30` names `SpatialColocResult` /
`delta_rho_map` as that map's calibrated successor. **Phase 12 ships the sample map, the control map
AND the difference.** The reasoning: a reader benefits from seeing where a difference *comes from*
rather than being handed a contrast they cannot interrogate.

**(1) `region_delta_rho` remains the PRIMARY named deliverable.** Per the `SpatialColocResult` sketch
at `src/results.jl:180-191`, computed as the per-region Monte-Carlo difference mirroring
`src/amortized/infer.jl:117-121`. `region_rho_sample` and `region_rho_control` are *additional exposed
artifacts*, never a replacement. The package keeps exactly one meaning of Δρ.

**(2) The read cost is `zero forward passes` extra, and datagen is unchanged.** A Δρ read already
performs both single-stack passes internally — `infer.jl:117-121` is
`mean(ρ_sample_draws .- ρ_control_draws)` over two `rho_draws` calls — so exposing
`region_rho_sample` and `region_rho_control` adds **zero forward passes** at read time; they are
values the computation already held and threw away. And the **training pool stays SINGLE-STACK**: the
net maps one summary to θ, and the pairing happens at *scoring* time via
`spike/validation/harness.jl:124-136` (`draw_simulate_infer_paired`), so datagen does **not** double.
This is stated explicitly because a paired pool would have pushed the measured 65–80 min / 50k to
~130–160 min against an append-only `P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150`.

**(3) The per-region Δρ\* ground truth comes from `harness.jl:124-136` at SCORING time, not from the
pool.** 12-16's parameter-coverage cross-check and 12-18's SBC obtain it by drawing two independent
prior samples and differencing their fields at read time — exactly as `draw_simulate_infer_paired`
already does. **12-09 does NOT gain a paired control field draw.**

> **WITHDRAWN, not deleted.** An earlier draft of this clause said 12-09 *did* gain such a draw. It is
> withdrawn here on the record, because it contradicted clause (2) immediately above, and freezing both
> into this document would have created a permanent, citable self-contradiction in a file that is never
> edited. 12-09's own acceptance criteria forbid any `_control` pool field or `rho_field_delta` column.

Note the scope of this clause, so it is not read against §4(a): Δρ\* is a **simulator** quantity used
only in the simulated parameter-coverage cross-check and in SBC. The SC3 predictive-coverage criterion
of §4 is **ground-truth-free and single-stack on both arms**, and consumes no Δρ\* at all.

**(4) 12-12's range guard is REPLACED, not tightened.** The old `[-2, 2]` bound admitted ρ ∈ [-1, 1]
silently — the one check that could have caught a ρ map mislabelled as a Δρ map is what hid it.
Tightening cannot repair that: with a **true** Δρ the range genuinely **is** [-2, 2], so a range check
simply cannot distinguish the two quantities. The check that can is a **STRUCTURAL IDENTITY**: assert
the result carries all three maps and that

    region_delta_rho ≈ region_rho_sample − region_rho_control

**elementwise**, to a stated tolerance. That is unfalsifiable-by-accident in a way no range ever is, and
it becomes available *precisely because* all three maps are now shipped. The `[-1, 1]` and `[-2, 2]`
bounds are retained as **sanity** checks explicitly labelled as such, and are never the guard.

**(5) Sentinel semantics are defined PER MAP.** `local_map.jl:48`'s wording is about the *difference*
and is meaningless for a single-stack map, so:

| Map | Sentinel means |
|---|---|
| `region_delta_rho` | "no evidence of a sample-vs-control difference" — the inherited meaning |
| `region_rho_sample` | "this region was not scorable in THIS stack" — a statement about measurability, not about correlation being zero |
| `region_rho_control` | as above, for the control stack |

The invariant is preserved in **all three**: a sentinel is **never** silently indistinguishable from a
measured zero, and every sentinel region carries a forced OOD flag. A region unscorable in *either*
stack is unscorable in the difference, so the difference's sentinel set is the **union** of the two
single-stack sets. A result built with no control stack at all is a legal state of the type, announces
itself through `meta.delta_rho_available = false`, and **must never be summarized as a Δρ of zero**.

### What Δρ is on the real arm: it is NOT computed there

The real-image data is, and stays, the six committed TIFFs = **two specimens** (six FILES, two
DATASETS). Two facts settle what Δρ can mean on them, and the answer is a **named limit**, not an
omission.

**Fact 1 — the SC3 criterion does not need a control.** D-09's leave-region-out construction masks
region *r* in ONE image's summary and scores that image's own observed entry at *r*. Every step is
single-stack. Predictive coverage is therefore computable on each specimen independently: **n = 2
images, ~64 regions each.** In the same breath: the **effective independent n is 2, not 128**. The
region-draws are **pseudo-replicates**, because `12-01-PLAN.md:198-204` pre-registers the dataset as
the independent unit — *"the 64 regions of one image share every nuisance and the global term"*, and
counting region-draws as independent would overstate power *"by up to 8× in sd terms"*. Any interval on
the real-arm coverage comes from the between-image spread, never from a binomial over 128.

**Fact 2 — Δρ needs two EXCHANGEABLE draws, and the two specimens are not exchangeable.** The argument
is exchangeability, and it is **not** an argument about shared nuisances: the simulated pair does not
share nuisances either — `harness.jl:124-136` draws two independent priors, exactly as the shipped Δρ
does. What the simulated pair has is **exchangeability**: both draws come from the same π(θ), so Δρ\*
is a **within-population** difference, and that is the quantity the net is calibrated against. The two
real anchors are `role = positive` / `truth = coloc` and `role = negative` / `truth = segregated` —
different specimen types, deliberately built to differ. Their contrast is **between-population**.
Scoring a between-population contrast against a gate calibrated on a within-population one would be two
arms measuring different quantities under one gate: the same defect class this section fixes.

**The counted corpus fact, and it was COUNTED, not inherited.** The Phase-8 external corpus supplies
**zero** images today. Counted this session from `corpus/manifest.csv` and `corpus/data/`:

| Counted | Value |
|---|---|
| `corpus/data/` | **empty** — no fetched content |
| Manifest data rows | **32** |
| Rows with `bytes = 0` | **32 of 32** |
| `physical-primary` rows | **2** — both `sha256 = PENDING-FETCH`, both `split = sealed_holdout` (Phase-16, untouchable here) |
| `simulated-secondary` rows | **30** — empty `sha256`, `split ∈ {eval, dev}` (16 eval, 14 dev) |
| Rows actually fetched | **0** |

The two `physical-primary` rows are `sealed_holdout` and reserved for Phase 16; the 30 CBS rows are
`simulated-secondary` (computer-generated) and cannot support a "real image" claim at any n. The
manifest's `role` column is `positive` / `negative` / `benchmark` — those are *experimental* positive
and negative controls, i.e. **different specimens with opposite ground truth**, not sample/control
stacks of the same specimen. Fetching would therefore not change this section's conclusion.

**So, as a named limit:** the per-region Δρ map is validated **on simulation**, where the sample and
control are two **exchangeable** draws from one prior so that Δρ\* is a within-population difference;
on real data the phase demonstrates per-region predictive calibration on **two specimens**, and
**no exchangeable sample/control pair is available** to validate Δρ itself. 12-19 records
`real_delta_rho_computed = false`. A specimen contrast between the two anchors may appear as an
illustrative figure, labelled a **specimen difference** and explicitly **not** the scored Δρ. This
milestone already carries several named limits and they have cost it nothing; a quiet pairing of the two
anchors would cost it a great deal.

**Reporting the two single-stack maps is therefore NOT the "different quantities" defect.** All three
maps are first-class deliverables everywhere, the SC3 criterion is single-stack on **both** arms and is
unaffected, and the one quantity that genuinely cannot transfer is named as not transferring.

---

## The legitimacy test, and the freeze

Applied to every section above, quoted verbatim from the precedent:

> *Would this change have been made, with this justification, by someone who had seen the gate's
> code but none of its outputs?*

**Yes, section by section.** §1 is the goal sentence set against the definition of a DeepSet. §2 is a
dependency-resolution constraint plus a measured property of the kernels' parametrization. §3 is a lane
decision plus two input/output contract mismatches. §4 is the observation that a real image has no
ground-truth ρ field and that the named baseline has no uncertainty — both readable from
`local_map.jl:55-79`. §5 is the definition of Δρ in `infer.jl:117-121` set against a counted corpus
audit. **Not one of them is a Phase-12 measurement, because at the sha in §0 there are none.**

**ROADMAP's SC1, SC2 and SC3 remain byte-present and citable.** This document adds a **pointer**, not a
rewrite. Any report that cites an amended criterion **must cite the original alongside it**; reporting
only the amended criterion would be exactly what this document exists to prevent.

**THIS DOCUMENT IS FROZEN. Later Phase-12 plans CITE it; none EDIT it.** A diff against this file after
its first commit is, by construction, a pre-registration breach. If a later finding contradicts
something written here, the correct action is a **separately named, separately dated** document that
records the contradiction — never an edit to this one.
