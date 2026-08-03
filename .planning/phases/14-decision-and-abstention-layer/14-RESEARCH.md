# Phase 14: Decision and Abstention Layer — Research

**Researched:** 2026-08-03
**Domain:** Bayesian multiple-testing (posterior-expected FDR), split conformal prediction, selective
prediction / risk-coverage, wired onto a frozen three-way evidence net in Julia, CPU-only, in a
byte-frozen package environment.
**Confidence:** HIGH on what is on disk (every wiring claim below was read out of a file or loaded
out of a `.jld2` in this session). HIGH on the FDR / conformal / selective-prediction method. MEDIUM
on the two judgement calls flagged as such (§C.3 monotonicity bar shape, §E.4 prevalence prior).

---

<user_constraints>

## User Constraints (from 14-CONTEXT.md)

### Locked Decisions

- **D-01: Bridge — build in `spike/`, against a `src/`-shaped API.** `decide_coloc(...)` is built in
  the spike research lane consuming Phase 13's three-way net (`spike/p13/three_way_net.jld2`,
  `spike/p13/net.jl`), but given the signature it would have in `src/` — taking a bundle + images and
  returning an `AbstractColocResult` subtype. This keeps SC1's three-way claim honest, keeps the
  Phase-7 GO untouched (matching Phase 13's D-01 research-lane posture), and makes later promotion
  mechanical rather than a rewrite.
  **Consequence for planning:** the decision layer is NOT shipped in v2.0. Say so plainly in any
  report; do not let "src-shaped" be read as "in src".

- **D-02: Hybrid — Bayesian FDR primary, hand-rolled split conformal as the hedge.**
  Bayesian FDR on the calibrated posterior is the primary rule in-distribution (this is the natural
  fit: a calibrated `P(H₀ | data)` gives direct expected-FDR control by sorting and thresholding, and
  SBC already backs the calibration). A **small hand-rolled split-conformal set** is the
  distribution-free hedge under misspecification. **`ConformalPrediction.jl` is NOT added.**
  **SC1 is AMENDED:** "conformal sets (ConformalPrediction.jl)" → conformal sets met *in substance*,
  hand-rolled; the library name is dropped. Rationale: the named library is the one option with a
  concrete verifiable cost (breaks the frozen-env guard) and the least payoff — split conformal is a
  sorted-score quantile, and the library wraps MLJ models rather than a NeuralEstimators posterior.

- **D-03: Conformal calibrates on held-out SIMULATOR draws; the corpus is a REPORTED CHECK only.**
  Fit conformal scores on held-out simulator draws (tight α, costs no corpus data), then **report
  observed coverage on the 30 open corpus rows beside it** as a named check.
  **This is deliberately not the strongest available claim.** Corpus-only calibration would be
  genuinely non-circular but caps α ≥ 1/31 ≈ **0.032** (n=30) and would **consume Phase 16's blind
  evaluation data**. Reporting both numbers shows the sim-vs-real gap instead of hiding it behind
  either source alone.
  **Named limit to carry into the manuscript:** the headline conformal guarantee is simulator-derived
  and therefore inherits the simulator's misspecification; the corpus check is what bounds that.

- **D-04: FDR is controlled PER IMAGE PAIR. Tiles are displayed but explicitly NOT FDR-controlled.**
  Per-pair is the only unit where the calibration provably holds — the posterior, the three-way BF,
  and the SBC proof are all per-pair. The per-tile map from `local_coloc_map` is still surfaced, but
  **labelled in the result as uncontrolled display, not a controlled call.**
  **Rationale, stated so it is not re-litigated:** claiming per-region FDR after Phase 12 returned NO
  would be exactly the "counting a closed negative as a success" failure corrected in ROADMAP at
  `e29c846`. `LocalColocMap` carries `delta_rho` and `ood_flag` but **no uncertainty field** — there
  is nothing to control against.

- **D-05: ASYMMETRIC fusion. SC2 is AMENDED.**
  - `OOD` fires → **ABSTAIN**
  - ambiguous conformal set → **ABSTAIN**
  - cross-method disagreement **alone** → **DECIDE**, and record the disagreement as a named field on
    the result. Disagreeing with Costes/Manders is the product thesis, not a failure mode.
  - cross-method disagreement **∧** OOD → **ABSTAIN** (independent reason to distrust ourselves)

  **SC2's literal `∨` is replaced by this rule.** Rationale: the literal OR silences the tool exactly
  where Phase 9 says it should speak, and three OR'd triggers compound — a tool that is usually
  silent is not useful. This preserves SC2's protective intent without gutting Phase 9.
  *Considered and not chosen:* a four-state output adding `DECIDED-CONTRA-CLASSICAL` — see Deferred.

- **D-06: A `false` OOD flag that means NOT CHECKED must never read as "in distribution."**
  The OOD input is **three-valued**: fired / clear / not-checked. `_recorded_ood_threshold`
  (`src/amortized/local_map.jl:93`) returns `nothing` when the gate report is absent, unreadable, or
  records a non-finite threshold, and `local_map.jl:70-72` states outright that *"a `false` flag on a
  non-sentinel tile means NOT CHECKED, not in distribution."*
  **Not-checked ABSTAINS by default**, with an explicit opt-out keyword the caller must set
  deliberately. The reason for silence is recorded in the result so it is auditable.
  **Rationale:** the failure this prevents is silent and scientific; the cost it imposes is loud and
  operational.

- **D-07: τ is INHERITED from Phase 13's artifact by loading it, with provenance asserted.**
  Load τ from the Phase-13 report rather than copying the literal into a Phase-14 constants file;
  **fail loudly** if the artifact is missing or its hash does not match. Copying the value is how two
  numbers silently diverge; re-deriving guarantees a second, different cut on the same hypothesis
  space.
  - **τ** (three-way cut) — inherited, measured by Phase 13, never tuned here.
  - **α_FDR** — a **user parameter**, per SC1's "user-set Bayesian FDR".
  - **α_conformal** — the miscoverage level for the hedge, **pre-registered separately**. Distinct
    from α_FDR and must not be conflated in code or prose.

### Claude's Discretion

- Exact conformal score function (posterior-probability-based vs BF-based) — pick during planning and
  justify against the calibration evidence.
- The concrete Bayesian-FDR estimator form (direct posterior-probability sorting vs a decision-risk
  formulation) — both satisfy D-02; pick the one that composes better with the three-way output.
- Result type shape and field names, provided it subtypes `AbstractColocResult` per the established
  hierarchy.
- Whether the risk-coverage curve (SC3) is computed on simulator draws, corpus rows, or both — but
  D-03's circularity caveat applies to whichever is chosen and must be stated.

### Deferred Ideas (OUT OF SCOPE)

- **Per-region / per-tile FDR control** — blocked on Phase 12's NO (no calibrated per-region
  uncertainty). Revisit in v2.1 alongside SPAT-07, which was itself deferred there.
- **Re-partitioning the corpus for a non-circular conformal guarantee** — genuinely better science
  (removes D-03's circularity caveat) but shrinks Phase 16's blind evaluation set and needs a
  pre-registered split ratio. A pre-registration matter, not a Phase-14 decision.
- **Promoting the three-way net into `src/`** — would make the decision layer shippable, but reopens
  the Phase-7 GO and pulls a DEV-seed research artifact onto the production path. Deliberately not
  done here (D-01).
- **Four-state output adding `DECIDED-CONTRA-CLASSICAL`** — considered during discussion as an
  alternative to D-05's asymmetric rule. Rejected for v2.0 as a non-standard output state needing
  extra defence, but it is the more transparent design if reviewers push back on the asymmetry.
- **Adding `ConformalPrediction.jl`** — would require amending or exempting the frozen-env guard and
  accepting the MLJ dependency tree. Revisit only if the hand-rolled conformal proves insufficient.

</user_constraints>

---

<phase_requirements>

## Phase Requirements

ROADMAP.md records `**Requirements**: TBD` for Phase 14 (`.planning/ROADMAP.md:377`) and
`.planning/REQUIREMENTS.md` contains **no** `DEC-*` / decision-layer requirement family
[VERIFIED: grep of `.planning/REQUIREMENTS.md` for `DEC-|Decision|abstain|FDR` returns only the
`### Demo & Decision` heading at line 58, whose three items DEMO-01..03 are all `[x]` and belong to
Phase 6].

**Consequence for the planner:** there are no requirement IDs to map. Coverage must be expressed
against the three amended success criteria (SC1-amended, SC2-amended, SC3) and against D-01..D-07,
exactly as Phase 13 did (`13-VERIFICATION.md` "Requirements Coverage": *"Coverage is instead
expressed as the D-01 through D-16 decision log"*). Do **not** invent requirement IDs — that is a
`.planning/REQUIREMENTS.md` edit and a user decision.

</phase_requirements>

---

## Summary

This phase has an unusually favourable evidence position and one unusually unfavourable one, and both
need to reach the planner before task decomposition.

**Favourable.** Everything Phase 14 consumes is on disk, loadable, and self-describing. The trained
net (`spike/p13/three_way_net.jld2`, 953 KB), its measured per-head correction, τ = 0.15, the
`consts.jl` sha both sides of the D-15 amendment, the Phase-11 basis it rides, and its 50 MB training
pool are all present and were opened in this session. More usefully still, the persisted gate report
`spike/p13/three_way_gate_report.jld2` already carries **4,000 labelled evaluation items with their
`logbf_coloc`, `logbf_exclusion`, `prob_coloc`, `prob_exclusion`, `class`, `rho_sample`, `rho_control`
and `lambda` columns** — a free, zero-recompute development substrate for every algorithm in this
phase. And the missing conceptual piece D-08 deliberately declined to supply — a prior over the three
hypotheses — was *measured* anyway and persisted in the net's own metadata as
`pi_class_masses = (exclusion = 0.31272, random = 0.47073, coloc = 0.21655)` over 200,000 prior draws.
That single NamedTuple is what turns two log Bayes factors into a proper three-class posterior, and it
is measured rather than chosen, which is exactly this project's standard.

**Unfavourable, and it must not be discovered during execution.** D-03's "report observed coverage on
the 30 open corpus rows beside it" **cannot be executed today.** `corpus/data/` is empty (0 entries),
every one of the 30 open rows carries `bytes = 0` and an empty `sha256` field, and Phase 12 already
recorded this as `corpus_images_available = 0`. Separately, the 30 open rows are single-image CBS
benchmark frames, while this net's input is a *paired* (sample, control) encoding — so even fetched,
they would not directly form the pairs `decide_coloc` consumes. The only real microscopy on disk is
`test/test_images/{positive,negative}/*_c{1,2,3}.tif` — **two specimens, no colocalization label,
both flagged OOD at ~4.2× their own threshold** by Phase 13. The honest Phase-14 posture is therefore
to *report the check as unexecutable with its reason*, run the n=2 fixture pair as a behaviour note,
and carry the gap as a named limit — not to quietly drop it and not to substitute a number that looks
like a coverage estimate.

**Primary recommendation:** Compute a three-class posterior by softmax over
`(log π_C + logBF_C, log π_R, log π_E + logBF_E)`; define the null as the composite
`{random ∪ exclusion}` so `v = P(H₀|Z) = 1 − P(coloc|Z)`; apply the direct posterior-probability
Bayesian-FDR rule (sort `v` ascending, accept the largest prefix whose **running mean** stays ≤ α_FDR)
**after** the abstention filter, and claim FDR over the DECIDED subset only; hedge with a hand-rolled
split-conformal set using the **LAC score `s = 1 − p̂(y|Z)`** at the finite-sample quantile
`sort(s)[⌈(n+1)(1−α_conf)⌉]`, abstaining on both `|set| ≥ 2` and `|set| = 0`; and validate SC3 with a
**normalized selective-skill statistic against oracle and random reference curves**, not with a raw
monotonicity assertion the empirical curve cannot satisfy.

---

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Three-class posterior from two log BFs | Spike decision layer (new) | — | D-08 explicitly declined to invent the prior; the prior *masses* were measured and persisted, so the composition belongs here and nowhere else. |
| Bayesian FDR over a batch | Spike decision layer (new) | — | Batch-level; no per-image object owns it. Pure arithmetic over a vector of `v`. |
| Split-conformal calibration + set construction | Spike decision layer (new) | — | Needs a held-out simulator draw stream and a frozen quantile; nothing existing does this. |
| Class-conditional evidence (log BF) | Frozen Phase-13 net (`spike/p13/net.jl`) | — | READ-ONLY. Phase 14 never retrains, never re-tunes, never re-derives τ. |
| Three-way class boundary | Frozen `three_way_label` (`spike/p13/labels.jl:121`) | — | Loaded, not reimplemented; it carries the D-05 two-factor cut and its counter-example. |
| OOD verdict | Shipped `src/amortized/ood.jl` (`ood_verdict`, `maha_score`, `id_threshold`) | Phase-13-lane density null (`p13_real_id_reference`) | The *shipped* function is consumed unchanged; the *nulls* must be fitted on the Phase-13 basis, because the shipped bundle's `zt` is a different standardizer (see Pitfall 4). |
| Not-checked detection | Spike decision layer (new), pattern from `src/amortized/local_map.jl:93` | — | D-06. `_recorded_ood_threshold` returning `nothing` is the pattern, but it reads a *shipped grid-8 gate report*; the spike lane needs its own three-valued resolver. |
| Classical comparator statistics | Frozen `spike/comparator/classical.jl` | — | READ-ONLY. `manders`, `costes_p`, `patch_correlation` consumed as-is. |
| μ → ρ_true rescaling for the cross-method comparison | Frozen `ghat` (`spike/simulator/ghat.jl:110`) | — | The only in-repo bridge between the classical per-patch-correlation scale and the ρ_true scale τ is defined on. Prevents a unit mismatch (Pitfall 5). |
| Result type | New `<: AbstractColocResult` subtype in `spike/p14/` | — | Phase-7 D-02 pattern; precedent is `ThreeHypothesisColocResult` (`spike/p13/result.jl:142`). |
| Per-tile display (uncontrolled) | Shipped `LocalColocMap` (`src/amortized/local_map.jl:73`) | — | D-04. Carried as a labelled field, never entering the FDR arithmetic. |

---

## Standard Stack

**No package may be added.** `spike/test/test_p12_decoupling.jl:176-179` asserts
`spike/Project.toml` and `spike/Manifest.toml` are byte-unchanged vs `HEAD` *and* untracked-clean,
and `:189-192` asserts the dependency name set is *exactly* a named sixteen. A `Pkg.add` fails all
three.

### Available (verified by reading `spike/Project.toml`, 16 entries)

| Package | Used for in this phase | Why sufficient |
|---------|------------------------|----------------|
| `Statistics` (stdlib) | `mean`, `quantile`, `median` | Conformal needs a sorted order statistic (see Pitfall 6), FDR needs a running mean. |
| `StatsBase` | `corspearman`, `sample`, `denserank` if needed | Already the house rank/robust-stats package. |
| `JLD2` | load the net / gate report / τ probe; persist the Phase-14 report | The house artifact format; `save_three_way`/`load_three_way` (`spike/p13/net.jl:459,496`) is the atomic-write idiom to copy. |
| `SHA` (stdlib) | D-07 provenance assertions | `p13_consts_sha()` (`spike/p13/net.jl:524`) is the existing pattern. |
| `Random123` | `Philox4x` counter-based per-index streams | Mandatory: every draw in this project is index-keyed (CONVENTIONS C-04). |
| `Flux` / `NeuralEstimators` v0.2.1 | forward pass only, via `load_three_way` | Never `train`. |
| `Distributions` | prior draws inside the simulator chain | Transitively required by `sample_prior`. |
| `CairoMakie` | risk-coverage curve, reliability plot, abstention-stratification figure | Already used for `spike/figures/p13_*.png`. |
| `HypothesisTests` | *optional* — a binomial CI on realized conformal coverage | Not required; `Distributions.Binomial` quantiles suffice. |
| `Images` / `ImageFiltering` / `ImageTransformations` / `Interpolations` / `CoordinateTransformations` | transitively required by the simulator and the read-only `src/` contract | — |
| `CSV` / `DataFrames` | reading `corpus/manifest.csv` metadata (metadata-only, allowlisted) | See Pitfall 8 on the corpus token ban. |
| `BenchmarkTools` | not needed here | — |

### Explicitly NOT to be added

| Package | Why it is tempting | Why it must not be added |
|---------|--------------------|--------------------------|
| `ConformalPrediction.jl` | SC1 names it | D-02 forbids it; it breaks three running assertions; and it wraps MLJ models, not a Flux/NeuralEstimators posterior. |
| `MLJ`, `ROCAnalysis` | AUC / ROC | `spike/test/runtests.jl` **asserts both are ABSENT** [VERIFIED: cited at `spike/p13/tau_probe.jl:98-103`]. The in-repo tie-aware `roc_auc` (`src/amortized/ood.jl:277`, twin at `spike/validation/ood.jl:319`) is the one to reuse. |
| `MultipleTesting.jl`, `Distances.jl`, any FDR package | the FDR rule | The rule is `sort` + `cumsum` + `findlast`. Three lines. |
| `CUDA` | — | `spike/test/test_p12_decoupling.jl:352` asserts no loaded module name contains `CUDA`. |

**Installation:** none. Any `Pkg` operation on the spike environment is a phase failure.

---

## Package Legitimacy Audit

**Not applicable — this phase installs zero external packages.** The environment is byte-frozen by a
running test (`spike/test/test_p12_decoupling.jl:176-179`) and the dependency name set is asserted to
be exactly sixteen (`:189-192`). No registry lookup, slopcheck run, or postinstall inspection is
required because no package is fetched. The planner should carry an executable guard re-asserting the
freeze from a Phase-14 test file (the house pattern: each phase re-states its own claim in its own
file — `test_p12_decoupling.jl:186-188` says so explicitly).

---

## A. Bayesian FDR on the three-way output

### A.1 What the net actually gives you

Read out of the artifacts in this session:

| Object | Where | Value / shape |
|--------|-------|---------------|
| Corrected log BFs | `three_way_log_bf(net, Z_pair, head_log_odds)` — `spike/p13/net.jl:349-358` | `(coloc = logit_C − hlo.coloc, exclusion = logit_E − hlo.exclusion, random = 0.0)` |
| `head_log_odds` | persisted in `three_way_net.jld2`, returned by `load_three_way` | `(coloc = 0.0029420438990930215, exclusion = 0.0023543271616353607)` [VERIFIED: JLD2 load] |
| Uncorrected head probabilities | `three_way_probs(net, Z_pair)` — `spike/p13/net.jl:372-376` | `(coloc = σ(logit_C), exclusion = σ(logit_E))` |
| π-level class masses | `three_way_net.jld2` → `meta.pi_class_masses` | `(exclusion = 0.31272, random = 0.47073, coloc = 0.21655)` [VERIFIED: JLD2 load] |
| π-level head log-odds | `meta.pi_head_log_odds` | `(coloc = −0.7764632135784509, exclusion = −0.40897645975438723)` |

**The semantics, stated precisely because everything downstream depends on it.** Each head is a
logit-BCE classifier trained *only on its own class pair* (D-11 participation weights,
`spike/p13/labels.jl:153`, loss at `spike/p13/net.jl:239-245`). Its Bayes-optimal logit is
`log[p(Z|A)/p(Z|B)] + log(q_A/q_B)` where `q` are the **training-subset** frequencies; subtracting
`head_log_odds` removes the frequency term and leaves a pure likelihood ratio
`[CITED: 13-REPORT.md §5, and the derivation is reproduced at spike/p13/net.jl:339-344]`.

Crucially, the training pool was stratified **at class granularity only** — `stratify_by_class`
(`spike/p13/datagen.jl:583`) accepts or rejects *whole items* on the basis of class alone and never
touches θ, which is stated at `:540-560` as the exactness condition: `q(θ|class) == π(θ|class)`. So
`log BF(C:R)` and `log BF(E:R)` are **π-scale likelihood ratios**, and the training frequencies have
already been divided out. **Do not subtract or re-apply any frequency term a second time.**

### A.2 Getting a scalar null posterior probability — the recommendation

Because both numbers are likelihood ratios against the same `random` reference (the whole point of
D-08's structural zero, enforced at `spike/p13/result.jl:159-161`), a prior over the three hypotheses
composes exactly:

```julia
# π over {exclusion, random, coloc}: MEASURED, read off the net's own metadata.
# NEVER retyped -- h = load_three_way(path); pi = h.meta.pi_class_masses
function three_class_posterior(logbf, pi_masses)
    a = (log(pi_masses.coloc)     + logbf.coloc,      # H_C
         log(pi_masses.random)    + 0.0,              # H_R  (structural zero)
         log(pi_masses.exclusion) + logbf.exclusion)  # H_E
    m = maximum(a)                                    # log-sum-exp, never exp() first
    w = exp.(a .- m)
    s = sum(w)
    return (coloc = w[1]/s, random = w[2]/s, exclusion = w[3]/s)
end
```

This is the *only* step in the phase that introduces information Phase 13 refused to introduce, and
the prior it introduces is measured (200,000 prior draws, recorded in `meta.pi_class_masses`), not
chosen. That distinction should be stated in the report in exactly those words.

**Is the null `{random}` or `{random ∪ exclusion}`? — RECOMMENDATION: `{random ∪ exclusion}`.**

The decision the phase must emit is `{coloc / not / ABSTAIN}` (ROADMAP SC1, 14-CONTEXT domain
paragraph). "Not coloc" is therefore the complement of coloc, and an item that is genuinely
*segregated* is unambiguously not a colocalization discovery — calling it one would be the worst
possible error, and it is the exact error the confusion matrix shows the net never makes (structural
zeros: 0 exclusion items called coloc, 0 coloc items called exclusion, `13-REPORT.md` §6b). So:

```
v_i = P(H₀ | Z_i) = 1 − P(coloc | Z_i)
```

Three reasons this composite pooling is right rather than merely convenient:

1. **Posterior probability is additive over a partition**, so `P({R ∪ E}|Z) = P(R|Z) + P(E|Z)` needs
   no approximation. The composite-null objection that bites *frequentist* FDR (which null value do
   you compute the p-value at?) does not arise for a posterior probability.
2. **The alternative — null = `{random}` alone — would leave exclusion items uncontrolled**, i.e. an
   exclusion item wrongly called coloc would count as neither a true nor a false discovery. That is a
   hole in the guarantee, not a tighter guarantee.
3. **It keeps the FDR unit and the decision unit identical** (D-04: per image pair, one call). A
   three-way FDR (one per hypothesis) would need three α's and a multiple-comparison story across
   hypotheses, which nothing in this phase's scope supports.

**Report the exclusion split beside the FDR number, always.** `v_i` is a sum of two terms; a batch
whose false discoveries are mostly `P(E|Z)` mass is scientifically different from one whose are
mostly `P(R|Z)` mass. Carry `(p_random, p_exclusion)` on the result so the composite is decomposable
by a reader.

### A.3 The estimator — exact form, and what it controls

**The direct posterior-probability procedure** [CITED: Newton, Noueiry, Sarkar & Ahlquist 2004,
*Biostatistics* 5(2):155-176; Müller, Parmigiani, Robert & Rousseau 2004 *JASA*; Müller, Parmigiani &
Rice 2007, "FDR and Bayesian Multiple Comparisons Rules", https://web.ma.utexas.edu/users/pmueller/pap/MPR07.pdf]:

Let `v_i = P(H₀ᵢ | data)` for `i = 1..n`. Sort ascending: `v₍₁₎ ≤ v₍₂₎ ≤ … ≤ v₍ₙ₎`. Define

```
FDR̂(k) = (1/k) · Σ_{j=1}^{k} v₍ⱼ₎        # RUNNING MEAN, not running max
k*      = max { k : FDR̂(k) ≤ α_FDR }     # empty ⇒ k* = 0, reject nothing
R       = the k* smallest-v items
```

Three details the planner must encode as tests:

- **Running MEAN, not running max.** The controlled quantity is the *posterior expected false
  discovery proportion*, `E[FDP | data] = (1/|R|)·Σ_{i∈R} P(H₀ᵢ|data)`, which is a mean by
  construction. A running max would be a different (and much more conservative) rule that controls
  the per-comparison posterior error probability instead. Write the test as: on a hand-built
  `v = [0.01, 0.02, 0.30]` with α = 0.10, the mean rule accepts `k* = 2` (means 0.010, 0.015, 0.110)
  and the max rule would accept `k* = 2` too — so **use a discriminating fixture**:
  `v = [0.01, 0.02, 0.15]`, α = 0.10 → mean rule `k* = 3` (means 0.010, 0.015, 0.060), max rule
  `k* = 2`. That fixture separates them.
- **The running mean of an ascending sequence is non-decreasing in `k`**, so `k*` is a prefix and the
  rule is exactly equivalent to a threshold `t*` on `v` (accept iff `v_i ≤ v₍ₖ*₎`). Assert the prefix
  property executably; it is what makes the rule a *rule* and not a subset search.
- **What is controlled.** `E[FDP | data] ≤ α_FDR`, *conditional on the model being right*. This is
  the Bayesian / posterior-expected FDR, **not** a frequentist long-run FDR guarantee. Do not write
  "controls the FDR at α" without the qualifier. Storey's pFDR is the frequentist quantity that
  coincides with it under a two-groups model; that connection is worth a sentence in the report and
  no more.

### A.4 Direct posterior sorting vs the decision-risk formulation — RECOMMENDATION

**Recommendation: implement the direct posterior-probability rule as the primary, and *report* the
implied decision-risk threshold beside it as a derived quantity. Do not make the loss-function form
the primary.**

The Müller et al. decision-theoretic result is that under a bivariate loss over (false discoveries,
false negatives) — e.g. `L_c = c·FD + FN` — the optimal Bayes rule is a **threshold on the marginal
posterior probability of the null**, specifically reject iff `v_i < c/(1+c)`
[CITED: Müller, Parmigiani & Rice 2007]. The two formulations are therefore structurally the *same
rule*; they differ only in what the user supplies.

Why the direct form wins here, in order of weight:

1. **SC1's own wording is "a user-set Bayesian FDR"** (`.planning/ROADMAP.md:379`). That is α, not a
   cost ratio. Making the user supply `c` and then reporting the induced FDR would satisfy the
   criterion only by translation.
2. **The loss form needs a 3×3 loss matrix this project has no basis to elicit.** With three
   hypotheses, "false discovery" is not one event: calling an *exclusion* item coloc and calling a
   *random* item coloc are different scientific errors with different costs, and a biologist's
   relative cost for them is not in this repo. Inventing one would be exactly the chosen-not-measured
   number the project's discipline forbids (`14-CONTEXT.md` specifics: *"measured rather than
   chosen"*). The direct form needs no such matrix — it needs only `P(coloc|Z)`.
3. **Composability with abstention.** The FDR rule operates on a vector of `v` and is indifferent to
   how that vector was filtered (see §E.2). A loss-minimisation formulation would have to price
   abstention into the loss as a third action, which is a genuinely larger design (Chow's reject
   option) and would make α_FDR, α_conformal and the abstention cost three interacting knobs instead
   of two independent ones. D-07 explicitly demands α_FDR and α_conformal stay unconflated; a shared
   loss would conflate them.

**What to report anyway, because it costs nothing:** the realized threshold `t* = v₍ₖ*₎` and its
implied cost ratio `c* = t*/(1 − t*)`. One line, and it makes the decision-theoretic reading
available to a referee without committing the implementation to it.

---

## B. Hand-rolled split conformal, zero new dependencies

### B.1 The construction, exactly

Split (inductive) conformal for classification, with a calibration set of `n` exchangeable labelled
points `(Z_i, Y_i)` drawn independently of the trained net
[CITED: Vovk/Papadopoulos inductive conformal; standard modern statement in Angelopoulos & Bates,
"A Gentle Introduction to Conformal Prediction"]:

1. Choose a nonconformity score `s(Z, y)` — higher = worse fit.
2. Compute calibration scores `s_i = s(Z_i, Y_i)`, `i = 1..n`.
3. Take the **order statistic**
   ```
   q̂ = sort(s)[ ceil(Int, (n + 1) * (1 - α_conf)) ]
   ```
   i.e. the `⌈(n+1)(1−α)⌉ / n` empirical quantile.
4. For a new `Z`, emit `C(Z) = { y ∈ {coloc, random, exclusion} : s(Z, y) ≤ q̂ }`.

**Guarantee.** Under exchangeability of the calibration points with the test point,
`P(Y_test ∈ C(Z_test)) ≥ 1 − α_conf`, marginally over the draw of both the calibration set and the
test point; and if the scores are almost surely distinct, also `≤ 1 − α_conf + 1/(n+1)`.

**Feasibility bound.** Step 3 requires `⌈(n+1)(1−α)⌉ ≤ n`, i.e. `α ≥ 1/(n+1)`. At `n = 30` this is
`α ≥ 1/31 ≈ 0.03226` — **which reproduces D-03's recorded 0.032 exactly**, confirming the CONTEXT's
arithmetic. [VERIFIED: derivation and D-03's number agree.]

**Dependencies:** `sort` (Base) and integer indexing. Nothing else. Not even `Statistics.quantile` —
and it must *not* be used (Pitfall 6).

### B.2 Score function — RECOMMENDATION: LAC, `s(Z, y) = 1 − p̂(y | Z)`

Three candidates were evaluated.

| # | Score | Form | Verdict |
|---|-------|------|---------|
| (i) | APS / adaptive | cumulative sorted posterior mass down to and including `y` | REJECT for this phase |
| (ii) | BF-margin | e.g. `−(logBF_y − max_{y'≠y} logBF_{y'})` | **REJECT — repo evidence is decisive** |
| (iii) | LAC / THR | `1 − p̂(y|Z)` on the three-class posterior of §A.2 | **RECOMMEND** |

**Why (iii).** The decisive evidence is in the gate report, and it is unusually strong for this
project:

- Both heads are **ECE-green on the very probabilities the score would use**: `ece_coloc = 0.011981`,
  `ece_exclusion = 0.012886` against `P13_ECE_GREEN = 0.05`, with **0 of 10 empty bins on both heads**
  and `vacuous_coloc = vacuous_exclusion = false` at head AUCs 0.9902 / 0.9883
  [VERIFIED: JLD2 load of `three_way_gate_report.jld2` in this session; matches `13-REPORT.md` §6c].
  LAC's known weakness is poor *conditional* coverage when `p̂` is miscalibrated
  [CITED: Sadinle, Lei & Wasserman, "Least Ambiguous Set-Valued Classifiers with Bounded Error
  Levels", JASA 2019]. That weakness is minimised precisely here.
- LAC's known strength is that it produces the **smallest average set size** among procedures with
  the marginal coverage guarantee [CITED: Sadinle et al. 2019]. Smaller sets ⇒ fewer `|set| ≥ 2`
  ambiguities ⇒ fewer abstentions. D-05's stated rationale is that *"a tool that is usually silent is
  not useful"*; the score choice should serve that, and LAC is the one that does.
- LAC is computed on **the same posterior the FDR rule uses**, so the two layers are coherent by
  construction. A conformal set built on a different object than the FDR rule's `v` can contradict it
  (set = {coloc} while `v` says reject), which would be an incoherence a referee would find.

**Why NOT (ii), the BF-margin score — and this is the load-bearing rejection.** `13-REPORT.md` §11
item 5 hands Phase 14 a named finding *by name*:

> **The F5 confident-tail precision shortfall → Phase 14.** … corrected log Bayes factors are
> unbiased in the bulk (mean deviation < 0.07 nat) and lose sub-nat precision in the confident tail,
> where the *sign* of the verdict is not in doubt but the *magnitude* is least determined. Phase 14's
> abstention layer consumes these magnitudes; **thresholds placed in the confident tail would rest on
> the least-determined part of the estimate.**

A margin score is a monotone function of the log-BF magnitudes, and its `(1−α)` quantile sits, by
construction, out in the confident tail. That is exactly the placement Phase 13 warned against. The
posterior probability `p̂` compresses the same tail into `[0,1]`, so its quantile sits where the
probability is *calibrated* (ECE-green, 0 empty bins), not where the nat-scale magnitude is least
determined. **This is the single strongest score-selection argument available and it should be quoted
in the plan.**

**Why NOT (i), APS.** APS buys better conditional coverage at the cost of systematically larger sets.
With `p̂` already ECE-green, the conditional-coverage premium is small; the set-size cost is not. It
is the right choice if the calibration were poor, which it is not. Record it as the documented
fallback if the realized `|set| ≥ 2` rate on the calibration set is implausibly low (which would
indicate over-confident `p̂`), and pre-register that fallback *before* looking.

### B.3 "Ambiguous conformal set" — operational definition for 3 classes

LAC applies **one global threshold `q̂` to `1 − p̂(y)`** across all classes. Hence
`C(Z) = { y : p̂(y|Z) ≥ 1 − q̂ }`, and three regimes exist:

| `|C(Z)|` | Meaning | D-05 action |
|---|---|---|
| **1** | exactly one class is plausible at level α_conf | **DECIDE** — singleton `{coloc}` → coloc; singleton `{random}` or `{exclusion}` → not coloc |
| **≥ 2** | two or more classes are plausible; the evidence does not separate them at α_conf | **ABSTAIN** (`trigger = :ambiguous_set`) |
| **0** | *no* class reaches the plausibility floor — `max_y p̂(y) < 1 − q̂` | **ABSTAIN** (`trigger = :empty_set`) |

**The empty set is not a corner case and it must not be collapsed into "ambiguous".** An empty
conformal set means the observed point is *more nonconforming than `(1−α)` of everything the
calibration set contained* — that is a distribution-free misspecification signal, arriving through a
completely different channel from the Mahalanobis/PP/noise OOD detector. It is arguably the most
valuable single output of the conformal hedge, because it is the one signal in this phase that does
not depend on the simulator being right about *densities*, only about *exchangeability*. Give it its
own trigger symbol, count it separately in the report, and cross-tabulate it against `ood_state` —
if the two agree, that is corroboration worth reporting; if they disagree, that is the more
interesting finding.

Note the arithmetic: with 3 classes and `p̂` summing to 1, `|C| = 0` requires `1 − q̂ > 1/3`, i.e.
`q̂ < 2/3`. Since `q̂` is the `(1−α)` quantile of `1 − p̂(true class)`, an empty set can only occur
when the calibration scores are concentrated near 0 (a confident, well-separated classifier) — which
is exactly this net. So the empty case is *reachable here*, unlike in many textbook settings. Confirm
this against the calibration draws and record the realized `q̂`.

### B.4 Sample-size arithmetic

Two distinct quantities; conflating them is a classic error.

**(a) Feasibility.** `α_conf ≥ 1/(n+1)`.

**(b) Tightness of the realized coverage.** Conditional on the calibration set, the realized coverage
of split conformal is exactly `Beta(n + 1 − ℓ, ℓ)`-distributed with `ℓ = ⌊(n+1)·α⌋`
[CITED: Vovk 2012, "Conditional validity of inductive conformal predictors"]. Hence
`sd ≈ sqrt(α(1−α)/(n+2))`:

| n | α = 0.10 → sd of realized coverage | α ≥ 1/(n+1) |
|---|---|---|
| 30 | **0.053** | 0.0323 |
| 200 | 0.021 | 0.0050 |
| 1,000 | 0.0095 | 0.0010 |
| **2,000** | **0.0067** | 0.0005 |
| 4,000 | 0.0047 | 0.00025 |

**Recommendation: `n_cal = 2000` fresh π-distributed simulator draws** at a new Phase-14 counter.
At α_conf = 0.10 that pins realized coverage to 0.90 ± 0.007 (1 sd), which is tight enough that the
reported number is a property of the method and not of the calibration draw. Cost is ~half the gate
run (§ Environment Availability: the 4,000-item gate set cost ~11 min CPU), i.e. ≈ 5-6 min.

**Contrast with the n=30 corpus floor.** Even ignoring that the images do not exist (§ Pitfall 1),
n=30 gives realized coverage 0.90 ± 0.053 — a ±10 % two-sigma band. That is not a guarantee, it is a
rumour. D-03's choice is correct and the sd number above is the concrete way to say so in the report.

**A separate `n_eval` set is needed for the *reported* coverage**, and its precision is binomial:
`sd = sqrt(α(1−α)/m)`. At `m = 2000`, ±0.0067. Recommend `n_eval = 2000` at a *third* counter.

### B.5 Feasibility confirmed

Split conformal as specified needs: `sort`, `ceil`, integer indexing, and a forward pass through an
already-loaded net. All Base/stdlib. **Zero new dependencies. Confirmed against the 16 entries of
`spike/Project.toml` read in this session.**

---

## C. Risk-coverage curves (SC3)

### C.1 The standard construction

Selective prediction [CITED: El-Yaniv & Wiener, JMLR 11 (2010); Geifman & El-Yaniv, NeurIPS 2017,
"Selective Classification for Deep Neural Networks"; Geifman, Uziel & El-Yaniv, ICLR 2019,
"Bias-Reduced Uncertainty Estimation for Deep Neural Classifiers" (AURC / E-AURC)]:

Given a confidence score `κ(Z)` and a loss `ℓ(ŷ, y)`, for threshold `θ`:

```
coverage(θ)       = P(κ(Z) ≥ θ)
selective_risk(θ) = E[ ℓ · 1{κ ≥ θ} ] / coverage(θ)
```

Sweep `θ` over the observed `κ` values, plot `selective_risk` against `coverage`.
`AURC = ∫ selective_risk d(coverage)`; `E-AURC = AURC − AURC_oracle`.

### C.2 What makes it monotone, and what breaks it

**Population statement (derivable, not cited):** the population risk-coverage curve is non-decreasing
in coverage **iff** `κ` is a monotone function of the conditional risk `E[ℓ | Z]` (up to ties).
Reasoning: increasing coverage adds the next-lowest-confidence points; the running average rises iff
every added point has expected loss ≥ the current average, which holds exactly when `κ` orders points
by expected loss. [ASSUMED — this is a one-line derivation, not a quoted theorem.]

**What breaks monotonicity in practice, all four live here:**

1. **Finite-sample noise at low coverage.** The denominator is tiny; a single hard item moves the
   selective risk by `1/k`. The empirical curve is a step function that jitters badly below ~10-20 %
   coverage. This alone will falsify any literal "the curve is monotone" assertion.
2. **`κ` not monotone in the true error probability** — i.e. miscalibration. Here `p̂` is ECE-green,
   so this is the *least* of the four risks, which is worth saying.
3. **Ties in `κ`.** With a sigmoid output ties are rare but the sweep must define a tie rule; use
   descending unique thresholds, the same convention `roc_auc` already uses (`src/amortized/ood.jl:281`).
4. **Multiple incommensurate abstention triggers.** *This is the specific hazard of D-05.* An OOD
   flag is a **binary** trigger; a conformal set size is an **integer**; only `p̂` is continuous. A
   risk-coverage *curve* requires a single scalar to sweep. **A curve cannot be drawn over the D-05
   fusion itself.** The planner must not conflate the two objects.

### C.3 RECOMMENDATION: what SC3 should actually measure

SC3 asks for two things at once — *"a monotone risk-coverage curve"* AND *"shows abstention
concentrates on hard/OOD cases"*. Separate them, and pre-register three falsifiable numbers before
any result exists. This project has now been shown four times that a bar can measure something other
than what it names (`.planning/STATE.md`, the `12-STAGE1-VERDICT.md` §7 families), so each bar below
is stated with what it would take to falsify it.

**Object 1 — the curve.** Sweep the scalar confidence `κ(Z) = max_y p̂(y|Z)`, which is *the same
quantity LAC thresholds*, so the curve and the conformal layer describe one ordering. Loss `ℓ = 1`
if the singleton call disagrees with the true three-way class, else 0. Plot alongside two reference
curves computed on the same data:

- **Oracle:** sort by actual correctness (all correct items first). Lower envelope.
- **Random:** shuffle. A straight line at the base error rate.

**Statistic S1 — selective skill (pre-register a floor).**
```
skill = (AURC_random − AURC) / (AURC_random − AURC_oracle)   ∈ [0, 1]
```
Normalized, so it is immune to the base error rate — which matters because the base error rate here
depends on the class mix and would otherwise let a favourable draw carry the bar. **Suggested
pre-registered floor: `skill ≥ 0.60`.** Falsified if the confidence ordering is no better than 60 %
of the way from random to oracle. *This bar is a judgement call and must be frozen before the run;
it is the one number in this section with no derivation behind it, and the report must say so.*
[ASSUMED]

**Statistic S2 — monotone trend, stated honestly.** Do **not** assert pointwise monotonicity. Assert:
```
spearman( coverage, selective_risk ) ≥ 0.95   over coverage ∈ [0.20, 1.00]
```
with the **0.20 floor pre-registered before the run** and its reason recorded (§C.2 item 1: below it
the denominator is too small for the ratio to mean anything). Report the raw curve over the *full*
range including the noisy part — showing it, not hiding it — and apply the statistic only on the
declared range. Anything else (monotonizing the curve, or picking the floor after seeing where it
becomes monotone) is the failure mode this project has paid for repeatedly.

**Statistic S3 — "concentrates on hard cases", and it must not be circular.**
```
hard_i := 1{ argmax_y p̂(y|Z_i) ≠ true_class_i }      # "would have been wrong"
AUC_hard := roc_auc( abstention_score | hard = 0 , abstention_score | hard = 1 )
```
using the in-repo tie-aware `roc_auc` (`src/amortized/ood.jl:277`), with
`abstention_score := 1 − max_y p̂(y|Z)`. `hard` uses the simulator's true label, which exists, so
there is no circularity. **Suggested pre-registered floor: `AUC_hard ≥ 0.80`.** [ASSUMED]

**Statistic S4 — "concentrates on OOD cases", measured on a separate arm.** Run the batch through
the four in-repo misspecification families — `misspec_texture`, `misspec_noise`, `misspec_optics`,
`misspec_background`, collected as `OOD_FAMILIES` at `spike/validation/ood.jl:499`, with
`ood_roc_over_grid` at `:612` — and report:
```
abstain_rate(misspecified at the strongest level)  −  abstain_rate(in-distribution)
```
**Suggested pre-registered floor: ≥ 0.50 (absolute).** Falsified if abstention is no more frequent
under a strong, deliberately injected misspecification than in-distribution. This is a *different*
claim from S3 and both are needed, because SC3's "hard/OOD" is a slash, not a synonym.

### C.4 Simulator draws, corpus rows, or both — RECOMMENDATION

**Recommendation: simulator draws only for every curve and every statistic; the real substrate gets a
reported n=2 behaviour note and nothing that looks like a rate.**

Grounds, in descending order of hardness:

1. **The corpus has no bytes.** `corpus/data/` is empty; all 32 rows record `bytes = 0`; the 30 open
   rows have an *empty* `sha256` field and the 2 sealed rows carry the `PENDING-FETCH` sentinel
   [VERIFIED: `ls corpus/data` → 0 entries, and `corpus/manifest.csv` read in full]. Phase 12 recorded
   `corpus_images_available = 0` for exactly this reason.
2. **The unit does not match even if fetched.** `decide_coloc` consumes a *pair* — `p13_pair_encode(Zs, Zc)`
   (`spike/p13/net.jl:398`) requires two summaries, and the D-05 cut is a two-factor rule over
   `(ρ_sample, ρ_control)` (`spike/p13/labels.jl:121`). The 30 CBS rows are single frames at a graded
   `coloc_ground_truth ∈ {0.0 … 0.9}` with no control acquisition. Pairing them would require a
   pairing convention that does not exist and would be invented here.
3. **Labels.** The curve needs a three-way class label; CBS supplies a scalar degree, not a
   `{coloc, random, exclusion}` class, and mapping degree → class needs a threshold — a *second* τ on
   a *different* scale, which D-07 forbids in spirit and §Pitfall 5 forbids in units.
4. **The 2 real fixtures are OOD.** `13-REPORT.md` limit D: both fixtures score density 703.30 against
   their own threshold 167.54 (4.198×). Under D-05 the correct output for both is ABSTAIN. That is a
   *result worth reporting* — "on the only real microscopy this project holds, the decision layer
   abstains, and that is the designed behaviour" — but it is n=2 and bounds nothing.

**The D-03 circularity caveat, stated for the record and required in any Phase-14 report:** the
risk-coverage curve, the conformal quantile and the FDR check are all computed on draws from the same
simulator the net was trained on. Train-joint equals eval-joint by construction (`13-REPORT.md` §6d).
**Every number in SC1/SC3 is therefore a well-specified-regime number and inherits the simulator's
misspecification in full.** It is a statement about the decision layer's internal coherence, not about
its performance on real microscopy. The bounding evidence D-03 intended — the corpus check — is
unavailable, so the bound is *absent*, not merely loose. Say that in those words.

---

## D. Wiring into what actually exists

### D.1 Loading the net — concrete call sequence

```julia
# spike/p14/... ; run as: julia --project=spike spike/p14/<script>.jl
isdefined(@__MODULE__, :ThreeWayEvidenceNet) || include(joinpath(@__DIR__, "..", "p13", "net.jl"))
isdefined(@__MODULE__, :three_way_label)     || include(joinpath(@__DIR__, "..", "p13", "labels.jl"))
isdefined(@__MODULE__, :ThreeHypothesisColocResult) ||
    include(joinpath(@__DIR__, "..", "p13", "result.jl"))
isdefined(@__MODULE__, :load_p13_basis)      || include(joinpath(@__DIR__, "..", "p13", "preconditions.jl"))

h = load_three_way(joinpath(@__DIR__, "..", "p13", "three_way_net.jld2"))
# h :: NamedTuple with fields:
#   net, head_log_odds, input_dim, num_summaries, width, tau, consts_sha,
#   cut_variant, schema_version, meta
```
[VERIFIED: `spike/p13/net.jl:496-514` (function body) and a live `JLD2.load` of the artifact in this
session.]

**Measured values in the artifact** [VERIFIED: JLD2 load, 2026-08-03]:

| Key | Value |
|---|---|
| `tau` | `0.15` |
| `cut_variant` | `:tau_contrast` |
| `head_log_odds` | `(coloc = 0.0029420438990930215, exclusion = 0.0023543271616353607)` |
| `input_dim` / `num_summaries` / `width` | `321` / `64` / `256` |
| `consts_sha` | `100e97a37bb470bb7cb4fbdbd8e33f219d83adcf61fe398a90cf3ad3391f3fa6` |
| `meta.grid` / `meta.n_cond` | `8` / `1` |
| `meta.pi_class_masses` | `(exclusion = 0.31272, random = 0.47073, coloc = 0.21655)` |
| `meta.pool_dir` | `…/spike/data/cache/p13/7c65a1a9…a29a1` (50 MB, present, **gitignored**) |
| `meta.phase11_sha` | `659cda09112c4f20d82e7fd1bc6dc22d57f3853ee245f3c6f781407eb68eb571` |

`input_dim = 321 = 5·8² + 1` = `three_way_input_dim(8, 1)` (`spike/p13/net.jl:132`). Never write it
as a literal — derive it, as that function's banner comment insists.

**Read surface, per pair:**
```julia
Z    = encode_conditioned_pair(Zs, Zc, [lambda])   # spike/p13/net.jl:429
lbf  = three_way_log_bf(h.net, Z, h.head_log_odds) # :349  -> (coloc, exclusion, random=0.0)
prob = three_way_probs(h.net, Z)                   # :372  -> UNCORRECTED head sigmoids
```
Note the **key-order trap the codebase already documents**: `three_way_log_bf` returns
`(coloc, exclusion, random)` while `ThreeWayLogBF` stores `(coloc, random, exclusion)`. The
normalizer `_p13_as_three_way_logbf` (`spike/p13/result.jl:179`) accepts any order and throws on a
key mismatch. **Always go through it; never build the NamedTuple positionally.**

### D.2 The `ThreeHypothesisColocResult` inner constructor, and what it constrains

```julia
struct ThreeHypothesisColocResult <: AbstractColocResult   # spike/p13/result.jl:142-164
    grid             :: Int
    posterior        :: Matrix{Float64}
    log_bf_vs_random :: ThreeWayLogBF        # NamedTuple{(:coloc,:random,:exclusion), NTuple{3,Float64}}
    ood              :: OODVerdict
    calibration      :: CalibrationMeta
    meta             :: NamedTuple
end
# inner constructor throws ArgumentError naming D-08 unless lbf.random === 0.0  (:159-161)
```

**Constraint on a downstream Phase-14 type.** Two options:

- **(a) Compose, do not inherit.** Carry a `ThreeHypothesisColocResult` as a *field* of the Phase-14
  result. Then the `random === 0.0` invariant is enforced once, at construction, and the Phase-14
  type does not have to re-state it.
- **(b) Re-state the invariant.** If the Phase-14 type stores the triple directly, it must call
  `_p13_as_three_way_logbf` and re-assert `random === 0.0`, or it becomes a second place where a
  differently-referenced triple can enter the system silently.

**RECOMMENDATION: (a).** It reuses the guard rather than duplicating it, matches the "new subtype,
never bolted-on fields" pattern of Phase-7 D-02 (`src/results.jl:161-166`), and makes the Phase-14
result strictly *additive* over Phase 13's — a reader can see exactly what the decision layer added.

### D.3 τ — exactly where it lives, and what can be hashed (D-07)

**τ = 0.15 appears in FOUR places on disk, and all four agree** [VERIFIED: all four loaded/read in
this session]:

| # | Location | Key / field | Value |
|---|----------|-------------|-------|
| 1 | `spike/p13/three_way_net.jld2` | `"tau"` | `0.15` |
| 2 | `spike/p13/tau_probe_report.jld2` | `"tau"` | `0.15` (with `tau_status = "measured"`, `tau_measured = true`) |
| 3 | `spike/p13/three_way_gate_report.jld2` | `"P13_TAU"` | `0.15` |
| 4 | `spike/p13/consts.jl:854` | `const P13_TAU` (Tier-2 append block, `:853-872`) | `0.15` |

**RECOMMENDATION for D-07's "loaded, not hardcoded":** treat **(1)** as the source of truth — τ
travels with the net that was trained under it — and cross-assert against (2), (3) and (4). The
constant at (4) is in scope anyway (`labels.jl` and `net.jl` both `include` `consts.jl`), so
including it as a fourth cross-check is free and strictly stronger than reading it as primary.

**Provenance material available to hash / assert** [all VERIFIED by loading in this session]:

| Handle | Value | Assertion to write |
|---|---|---|
| `h.consts_sha` (net, train-time) | `100e97a3…3fa6` | `== P13_CONSTS_SHA.pre_amendment` (`spike/p13/net.jl:557`) |
| `p13_consts_sha()` (file, now) | `d6a6e63b…e247` | `== P13_CONSTS_SHA.post_amendment` (`:558`) |
| `tau_probe_report["repo_sha"]` | `cecb69c5…9f80` | `== P13_TAU_PROBE_SHA` (`consts.jl:856`) |
| `tau_probe_report["simulator_sha"]` | `78dc37f5…5f63` | `== P13_TAU_SIMULATOR_SHA` (`consts.jl:857`) |
| `tau_probe_report["measurement_sha256"]` | `06eb1bc8…c1e6` | self-consistency of the persisted curve |
| `tau_probe_report["tau_bar"]` | `0.9` | `== P13_TAU_AUC` |
| `P13_TAU_MEASURED_AUC` | `0.91635625` | `≥ P13_TAU_AUC` and `== curve[idx].auc` |
| `h.meta.phase11_sha` | `659cda09…b571` | net rode the expected Phase-11 basis |
| `gate_report["net_consts_sha"]` | `100e97a3…3fa6` | `== h.consts_sha` |

**⚠ A TRAP THE PLANNER MUST BE TOLD ABOUT.** `tau_probe_report.jld2` records
`consts_sha256 = 640ec90edb80ef62eacbbf64c7b40cb9ca87c59d28a3a7438ef9d382b0aca4cf` — a **THIRD**
value, matching *neither* literal in `P13_CONSTS_SHA`. This is **correct and expected**: the probe
ran on 2026-07-27 *before* the Tier-2 τ block was appended to `consts.jl` (which is precisely the
auditable ordering `consts.jl:842-846` cites as the evidence τ was measured before it was locked). A
naive Phase-14 assertion `probe["consts_sha256"] ∈ values(P13_CONSTS_SHA)` **would fail**, and an
executor "fixing" it by widening the check would destroy the integrity argument. Write the assertion
as: `probe["consts_sha256"] ∉ values(P13_CONSTS_SHA)` — i.e. **assert the divergence**, with the
reason in the message. That is the same shape as the existing pinned-both-sides assertion at
`run_three_way_gate.jl:369-370` and for the same reason.

Also record honestly: `tau_probe_report["repo_dirty_at_run"] == true`. The probe ran on a dirty tree.
Surface it in the Phase-14 report; do not assert it false.

### D.4 The OOD channel — signatures and the three-valued input

**Shipped surface** (`src/amortized/ood.jl`, reached read-only):

```julia
fit_ood_nulls(Ztrain::AbstractMatrix; variant = :min)         # :77  -> (cont, μS, C)
maha_score(nulls, z::AbstractVector) -> Float64                # :103
cov_cols(Zc::AbstractMatrix) -> Matrix                          # :88
pp_mismatch_score(m, Z_obs; reps, rng, imsize, N) -> Float64    # :166
noise_features(pair) -> Vector{Float64}                         # :211 (10 features)
fit_noise_null(F::AbstractMatrix) -> (μ, σ, C)                  # :243
noise_score(nn, pair) -> Float64                                # :259
roc_auc(neg, pos) -> (fpr, tpr, auc)                            # :277  (tie-aware Mann-Whitney U)
id_threshold(id_scores; q = OOD_ID_QUANTILE) -> Float64         # :302  (OOD_ID_QUANTILE = 0.95, :61)
youden_j(fpr, tpr) -> (j, fpr, tpr, idx)                        # :310  NEVER the gate
ood_flag(nulls, m, Z; maha_thr, pp_thr, rng, reps, imsize, N) -> Bool   # :329
ood_verdict(ood_nulls, Zs, Zc = nothing; pair = nothing, with_pp = true,
            rng = nothing, reps = OOD_PP_REPS, imsize = (256,256)) -> OODVerdict  # :363
```

`OODVerdict` is `struct OODVerdict; score::Float64; flag::Bool; per_channel::NamedTuple; end`
(`src/results.jl:93-97`).

**The critical mechanic for D-06, read directly out of the source (`src/amortized/ood.jl:394-396`):**

```julia
thr  = haskey(ood_nulls, :thr) ? ood_nulls.thr : nothing
flag = thr === nothing ? false : fused > thr
return OODVerdict(fused, flag, per)
```

**So `flag == false` is genuinely two-valued at the source.** When `:thr` is absent from `ood_nulls`,
`ood_verdict` returns `flag = false` **for every input**, unconditionally. This is precisely the
hazard `src/amortized/local_map.jl:110-119` names in prose (*"Without one, `ood_verdict`
short-circuits to `flag = false` and the per-tile flag is VACUOUS… which silently understates
risk"*), and the `nothing`-returning resolver at `:93-108` is the shipped detection pattern.

**RECOMMENDED three-valued resolver for the spike lane:**

```julia
# :fired | :clear | :not_checked   -- NEVER a Bool
function p14_ood_state(ood_nulls, verdict::OODVerdict)
    haskey(ood_nulls, :thr) || return :not_checked
    t = ood_nulls.thr
    (t isa Real && isfinite(t)) || return :not_checked   # mirrors _recorded_ood_threshold's guard
    return verdict.flag ? :fired : :clear
end
```

Mirroring `_recorded_ood_threshold`'s exact acceptance predicate (`t isa Real && isfinite(t)`,
`local_map.jl:107`) matters: a non-finite threshold is *also* not-checked there, and a Phase-14
resolver that accepted `Inf` would silently disagree with the shipped semantics.

**Where does `:thr` come from in the spike lane?** Not from the shipped bundle — see Pitfall 4. The
Phase-13 precedent is `p13_real_id_reference(pool_dir)` (`spike/p13/run_p13_realimage.jl`), which
fits the density null on the standardized summaries of *the net's own training pool* (both `Zs` and
`Zc` halves) and takes `id_threshold(scores; q = 0.95)` — the copied shipped quantile. `pool_dir` is
read off `h.meta.pool_dir`, never guessed. **Reuse that construction verbatim.** It gives a
single-channel (density) null on the correct basis, which yields `:fired`/`:clear`; and the honest
statement is that the noise and PP channels are **not** wired in this lane, which should be recorded
rather than papered over.

**Accessor contract implication.** `AbstractColocResult` requires `is_ood(r) -> Bool`
(`src/results.jl:68-72`). A three-valued state cannot round-trip through it. **RECOMMENDATION:**
```julia
is_ood(r::P14Result)    = r.ood_state === :fired      # interface-conformant, DOCUMENTED as lossy
ood_state(r::P14Result) = r.ood_state                 # NEW accessor, the one the decision layer uses
```
and document on `is_ood` that a `false` return collapses `:clear` and `:not_checked`, that the
decision layer never routes through it, and that a consumer who cares must call `ood_state`. This is
the same "documented rather than silent" treatment `bayes_factor(::ThreeHypothesisColocResult)`
received (`spike/p13/result.jl:205-219`), and the phase should follow that precedent explicitly.

### D.5 Cross-method disagreement — exact signatures and a proposed definition

**Available classical estimators** (`spike/comparator/classical.jl`, frozen, read-only):

```julia
manders(mci)                       -> (M1, M2)   # :66   NaN-safe, Otsu thresholds from the mci
pearson_whole(mci)                 -> Float64    # :85
spearman_whole(mci)                -> Float64    # :95
patch_correlation(mci; method = :spearman) -> Float64  # :108  mean per-patch ρ over the FROZEN 8×8 grid
costes_rng(master_seed, idx)       -> Philox4x   # :135
costes_p(master_seed, idx, mci; n = COSTES_N_SCRAMBLE, block = COSTES_BLOCK_PX) -> Float64  # :181
```
Constants (`spike/comparator/config.jl`): `COSTES_N_SCRAMBLE = 200` (:37), `COSTES_BLOCK_PX = 7`
(:38), `COSTES_SALT = 0xA5A5A5A5DEADBEEF` (:39), `COSTES_ALPHA = 0.05` (:43),
`DIVERGENCE_WARN = 0.30` / `DIVERGENCE_FAIL = 0.60` (:50-51), `MASTER_SEED = 0x00C0FFEE` (:56).

**PROPOSED OPERATIONAL DEFINITION — two channels, both on one scale, no new constant.**

**Channel 1 — classical three-way call, via `ghat`.**
```julia
rho_hat_s = ghat(patch_correlation(mci_s))        # spike/simulator/ghat.jl:110
rho_hat_c = ghat(patch_correlation(mci_c))
classical_call = three_way_label(rho_hat_s, rho_hat_c; tau = tau_loaded)   # labels.jl:121
disagree_1 = (classical_call !== amortized_call)
```
Why this is the right construction, and each link was verified:

- `patch_correlation` is the **mean per-patch correlation over the frozen 8×8 grid**
  (`classical.jl:108-114`), reducing `skipmissing` entries with `_safe(mean, ...)`.
- `mbar_from_summary` — the statistic τ was *measured on* — is `sum(vals .* mask)/sum(mask)`
  (`spike/p13/tau_probe.jl:135-147`), i.e. the mean over *present* patches. **The same quantity.**
- `ghat(μ) -> ρ_true` is the frozen, monotone, clamped piecewise-linear map from *induced mean
  per-patch correlation* to ρ_true (`spike/simulator/ghat.jl:44-46, 103-112`; knots
  `GHAT_MU_KNOTS`/`GHAT_RHO_KNOTS` at `:93,:95`; `corspearman(ρ_grid, E[μ]) = 1.0` at `:53`).
  **It already runs in the direction needed** — μ in, ρ out — so no inversion is written here.
- `three_way_label` then applies **the same τ, on the same scale**, giving one cut and one number.
  No second threshold is invented, which is what D-07's spirit requires.

**Channel 2 — Costes significance.**
```julia
p_s = costes_p(P14_DEV_SEED, idx, mci_s)          # NOT comparator MASTER_SEED -- see Pitfall 7
costes_call = p_s < COSTES_ALPHA ? :coloc : :not_coloc
disagree_2 = (costes_call === :coloc) != (amortized_call === COLOC)
```
Genuinely independent of channel 1: it is a block-permutation test, not a correlation threshold.

**Record as a NamedTuple field, not a Bool** (D-05 says "record the disagreement as a named field"):
```julia
cross_method = (classical_call, costes_p_sample, costes_call,
                rho_hat_sample = rho_hat_s, rho_hat_control = rho_hat_c,
                disagree_label = disagree_1, disagree_costes = disagree_2,
                disagree_any = disagree_1 || disagree_2,
                basis = :ghat_rho_true_scale, tau = tau_loaded)
```
`basis` is not decoration — it records *which scale the comparison was made on*, which is the thing
Pitfall 5 says will otherwise be forgotten.

**Alternative considered and NOT recommended:** thresholding `|M1 − M2|` or a Manders cut. Manders
coefficients are fraction-of-signal-above-threshold quantities on `[0,1]` with no map to ρ_true, so
any cut on them is a new, chosen constant — exactly what this phase must avoid. Report `manders(mci)`
as a descriptive column; do not derive a call from it.

### D.6 The `AbstractColocResult` accessor contract

`src/results.jl:35-79`. Supertype `AbstractColocResult` (`:43`); `_iface_error(r, f)` (`:47-49`)
names the accessor, the type, and the four required accessors. Every subtype must implement:

| Accessor | Contract | Phase-14 recommendation |
|---|---|---|
| `delta_rho(r)` | Δρ (sample − control) draws | **Fall through to `_iface_error` unless the run genuinely carried them.** This is the precedent at `spike/p13/result.jl:240-243`, whose docstring gives the reason: fabricating Δρ from a three-way path publishes a number no measurement backs, and the D-05 counter-example (ρ_s = 0.8 under ρ_c = 0.9) is why. Do the same. |
| `bayes_factor(r)` | "the colocalization Bayes factor evidence" (`:59-65`) | Return `log BF(coloc : random)`, documented as lossy, exactly as `spike/p13/result.jl:205-220` does. |
| `is_ood(r)` | Bool | `r.ood_state === :fired`, documented as lossy (§D.4). |
| `posterior_draws(r)` | parameter × draw matrix | Pass through from the carried `ThreeHypothesisColocResult`. |

New accessors layered on top (never bolted-on fields — Phase-7 D-02):
`decision(r) -> Symbol` in `{:coloc, :not_coloc, :abstain}`; `abstain_reason(r) -> Symbol`;
`ood_state(r)`; `conformal_set(r)`; `null_posterior(r)`; `cross_method(r)`;
`local_map(r)` (labelled uncontrolled, D-04).

**Batch object.** `decide_coloc` is inherently batch-level (FDR is a batch quantity). Recommend a
`P14BatchDecision` carrying the vector of per-pair results plus the batch fields
`(alpha_fdr, k_star, threshold_t_star, implied_cost_ratio, n_decided, n_abstained,
abstain_reason_counts, posterior_expected_fdp, class_prior_used, fdr_scope = :decided_subset_only)`.
The last field is not cosmetic — it is the machine-readable form of the honest claim in §E.2.

### D.7 The decoupling guards, quoted verbatim

The planner must not break these, and 14-CONTEXT asks for analogous guards for D-03 and D-06.

**Environment freeze — `spike/test/test_p12_decoupling.jl:170-183`:**
```julia
@testset "the spike environment is byte-frozen" begin
    if _HAS_GIT
        @test success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                          dir = P12_REPO_ROOT))
        @test isempty(readchomp(Cmd(`git status --porcelain -- spike/Project.toml spike/Manifest.toml`;
                                    dir = P12_REPO_ROOT)))
    else
        @info "git unavailable — skipping the executable spike-environment freeze assertions"
    end
end
```

**Exactly-sixteen dependency set — `:185-193`:**
```julia
@test Set(keys(Pkg.project().dependencies)) == Set([
    "BenchmarkTools", "CSV", "CairoMakie", "CoordinateTransformations", "DataFrames",
    "Distributions", "Flux", "HypothesisTests", "ImageFiltering", "ImageTransformations",
    "Images", "Interpolations", "JLD2", "NeuralEstimators", "Random123", "StatsBase"])
@test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
```

**Sealed holdout — `:210-297`.** The property is *"no byte of a sealed corpus IMAGE is ever read"*,
**not** *"the string `corpus` appears nowhere"* (`:211-213`). The structure is four checks:
(a) no image load resolves onto the corpus, in ANY phase source, allowlisted or not (`:236-244`);
(b) non-allowlisted sources may not mention the corpus at all, comment-stripped first, so prose is
free and a code reference is not (`:246-256`); (c) allowlisted metadata readers get the narrower
property — metadata yes, the sealed image directory never, *"neither spelled as a path nor reached
through the constant `corpus/config.jl` binds to it"* (`:258-270`); (d) `corpus/` is byte-unchanged
and untracked-clean (`:291-296`), plus a **positive** assertion that the permitted real-image path
(`test/test_images`) is actually exercised, so the three negatives cannot be satisfied by a phase
that reads no images at all (`:272-290`).

**Recommended Phase-14 analogues** (14-CONTEXT explicitly asks for these):

1. **D-03 guard** — a `test_p14_decoupling.jl` re-stating (a)-(d) over the Phase-14 source set, with
   its own `_P14_CORPUS_METADATA_ALLOWED` allowlist, plus the positive `test/test_images` assertion.
   *Add one Phase-14-specific check:* assert that no Phase-14 source contains a literal `0.032` or
   `1/31` outside a comment — i.e. that the corpus-derived α floor never became an executable
   default (which would be the corpus silently calibrating the hedge).
2. **D-06 guard** — assert executably that `p14_ood_state` returns `:not_checked` when `:thr` is
   absent AND when `:thr` is non-finite, and that `decide_coloc` with `ood_state = :not_checked` and
   the opt-out keyword left at its default returns `:abstain` with `abstain_reason = :ood_not_checked`.
   Also grep the executable source for a `!verdict.flag` or `verdict.flag == false` pattern used as a
   proxy for "in distribution" — that is the exact bug D-06 exists to prevent, and it should fail the
   suite if reintroduced (the same grep-the-source technique `spike/test/test_p13_tau.jl` uses to keep
   a grid-extension path from ever appearing).
3. **D-07 guard** — the four-way τ agreement plus the sha table of §D.3, including the *assert-the-
   divergence* form for the probe's third `consts_sha256`.
4. **Seed guard** — a `test_p14_consts.jl` mirroring `test_p13_consts.jl`'s disjointness block
   (`spike/p13/consts.jl:730-741`): P14 seeds ∉ forbidden set, P14 counters pairwise distinct, and —
   the Phase-14-specific one — **the full cross-product of Philox key words is pairwise distinct**
   between the P14 family and every P13 key word (`P13_DEV_SEED ⊻ P13_SALT ⊻ c` for
   `c ∈ {1,2,3,4,5,99}` and `P13_DEV_SEED ⊻ P13_SALT`). See Pitfall 2.

---

## E. Pitfalls specific to this phase

### E.1 Where the Bayesian FDR guarantee silently fails

| Failure | Mechanism | Detection / mitigation |
|---|---|---|
| **Prevalence mismatch** — *the biggest one* | `v_i` depends on `π`. `pi_class_masses` is the **simulator prior's** class mix (E 0.313 / R 0.471 / C 0.217), not any real batch's prevalence. A user batch that is 90 % coloc makes every `v_i` too large and the rule too conservative; a batch that is 1 % coloc makes them too small and the FDR claim false. | Make the class prior an **explicit keyword** defaulting to the measured `h.meta.pi_class_masses`. Report an FDR **sensitivity curve** over a pre-registered grid of `π_C` (e.g. 0.05, 0.10, 0.217, 0.40, 0.70) rather than one number. Do **not** default to empirical-Bayes estimation of `π` from the batch — it is a second estimated quantity with its own failure modes and would make the guarantee circular. Offer it as an opt-in and label it as such. |
| **Miscalibration** | The guarantee is an identity given correct posteriors; nothing rescues it if `p̂` is wrong. | Carry the Phase-13 `CalibrationMeta` (ECE 0.0120 / 0.0129, green, 0 empty bins, AUC 0.99, `vacuous = false`) on the result, so every FDR number ships with the calibration evidence beside it. That is `p13_calibration_meta`'s whole design (`spike/p13/result.jl:337-352`). |
| **Composite null pooling done on the wrong object** | Summing `P(R\|Z) + P(E\|Z)` is exact; summing *Bayes factors* or *logits* is not. | Assert `p_coloc + p_random + p_exclusion ≈ 1` to `1e-12` on every item. Cheap, catches every arithmetic slip. |
| **The prior-atom asymmetry** | The θ prior carries clamp atoms — ~4.85 % of prior mass at ρ = −0.99 vs ~1.91 % at +0.99, a ~2.5:1 asymmetry against the exclusion end (`13-REPORT.md` §9, the evaluated non-finding; `spike/p13/tau_probe.jl:335-340`). The measured `pi_class_masses` inherit it. | Phase 13 evaluated the condition for a named limit and it *did not hold for discrimination* (deep-tail exclusion AUC 0.997257). But it enters the **prior term** of §A.2 directly, which Phase 13 never used. Record it as a stated property of the prior, and it is one more reason to report the π-sensitivity curve. |
| **Independence assumed across the batch** | `E[Σ 1{H₀ᵢ}] = Σ v_i` needs only linearity of expectation, **not** independence — so this one does *not* bite. | Say so explicitly; it is the standard reviewer question and the answer is clean. |
| **Log-sum-exp omitted** | `exp(logbf)` on a `logbf` range of roughly ±22 nats (`three_way_gate_report.jld2` continuity block shows `max_abs_diff = 22.02`) overflows/underflows to 0 or Inf. | Always subtract the max before `exp`. Test with a synthetic `logbf = 800.0`. |

### E.2 Does routing to ABSTAIN before the FDR sort change what is controlled?

**Short answer: yes, it changes the *scope* of the claim, and no, it does not invalidate it — but
only because the quantity is Bayesian.** State this in the report; it is the subtlest thing in the
phase.

**The argument.** The controlled quantity is
`E[FDP | data] = E[ (Σ_{i∈R} 1{H₀ᵢ}) / |R| | data ]`. The rejection set `R` is a deterministic
function of the observed data (the abstention filter and the FDR sort are both data-measurable). So
`R` is σ(data)-measurable and pulls out of the conditional expectation:
`E[FDP | data] = (1/|R|)·Σ_{i∈R} P(H₀ᵢ | data) = (1/|R|)·Σ_{i∈R} v_i`.
**This holds for *any* data-dependent selection of `R`, including one that filtered on OOD status,
conformal ambiguity, or anything else.** Unlike frequentist FDR — where a data-dependent
pre-selection is a genuine selective-inference problem requiring correction — the Bayesian posterior
quantity is immune, because you already conditioned on the data. [ASSUMED — this is the standard
conditioning argument, stated here as a derivation rather than a quoted theorem; a referee should be
able to check it in two lines.]

**Therefore: ABSTAIN FIRST, then FDR-sort the survivors.** Ordering matters and this order is
strictly better:

- Abstaining first *improves the premise*. The guarantee needs the posteriors to be right; abstention
  routes away precisely the items where they are least likely to be (OOD-flagged, not-checked,
  conformally ambiguous or empty). The retained set is the set on which the model is most defensible.
- The reverse order — FDR-sort the whole batch, then abstain from some accepted items — leaves the
  residual accepted set with a running mean that may **exceed** α, because removing items changes
  both numerator and denominator in an uncontrolled way. That ordering is wrong and the plan should
  say so explicitly so it is not re-derived.

**The honest claim, in the exact words the report should use:**

> The posterior expected false discovery proportion is controlled at α_FDR **over the DECIDED subset
> only** — the batch minus abstentions. It is **not** controlled over the whole batch: an abstained
> item is neither a discovery nor a non-discovery and appears in neither the numerator nor the
> denominator. The guarantee is conditional on the model (calibrated posteriors and a correct class
> prior) and is a posterior expectation, not a frequentist long-run rate.

**And the mandatory companion number.** *Always report the decided fraction (coverage) beside the FDR
level.* A rule that abstains on 95 % of a batch hits any α trivially. `(α_FDR, n_decided/n_total)` is
one quantity, not two, and quoting either alone is misleading. This should be enforced structurally:
make the report writer take both or neither, the same way `p13_calibration_meta` makes `auc` a
required keyword so a calibration verdict can never ship without its discrimination number
(`spike/p13/result.jl:313-315`).

### E.3 Live landmines from this repository's history

| # | Landmine | Why it fires here | Guard |
|---|---|---|---|
| 1 | **The corpus has no images** | D-03's reported check is unexecutable; `corpus/data/` is empty, `bytes = 0` on all 32 rows, `corpus_images_available = 0`. | Report the check as unexecutable with its reason; never substitute a number. Plan a task whose deliverable is *the recorded absence*, not a coverage figure. |
| 2 | **Phase 12's permutation artifact** (user memory `phase12-permutation-artifact`) | An orthogonal permutation preserves norms and spectra, so a wrong-basis bug is O(1) yet invisible to structural verification — it silently inverted a phase's headline finding. Phase 14 has an exactly analogous exposure: **the class order.** `three_way_probs` returns `(coloc, exclusion)`; `three_way_log_bf` returns `(coloc, exclusion, random)`; `ThreeWayLogBF` stores `(coloc, random, exclusion)`; the gate's `confusion_class_order` is `["exclusion","random","coloc"]`; `pi_class_masses` is keyed `(exclusion, random, coloc)`. **Five orderings.** A positional index anywhere silently permutes the posterior, and every aggregate statistic (ECE, AUC, FDR) is invariant to a consistent relabelling — so it would pass structural verification. | Everything by NAME. Never `[1]`/`[2]`/`[3]` on any of these. Add a directed test: build a fixture where the three classes have *deliberately unequal* known posteriors (e.g. 0.7/0.2/0.1) and assert each named field lands on its intended value. An equal-mass fixture would not catch a permutation. |
| 3 | **The `SPEEDUP_GATE` masks later includes** (`13-REPORT.md` §11 item 6, `deferred-items.md` D-13-A) | The full `spike/test/runtests.jl` exits 1 at the Phase-4 gate, so the whole Phase-13 include block never runs; Phase-13 tests were run per-file. The same will happen to Phase 14. | Do **not** plan any verification step that depends on the aggregate suite being green. Every Phase-14 test file must be independently runnable: `julia --project=spike spike/test/test_p14_*.jl`. Record per-file pass counts in the report, as Phase 13 did. |
| 4 | **The KDE-baseline attrition** (skill `references/bayes-factor-gate.md`; user memory `v2-go-decision-named-limits`) | The KDE Bayes factor is invalid past \|log BF\| ≈ log(999) ≈ 6.9, and Phase 13 *declined* it as a validation target. | Do not reach for `compute_BayesFactor()` / KernelDensity / QuadGK as any kind of reference here. They are not direct deps and `spike/test/runtests.jl:108-150` asserts their absence. There is no baseline comparison in this phase's scope. |
| 5 | **"The bar was wrong, not the model" — four times** (`.planning/STATE.md`, `12-STAGE1-VERDICT.md` §7) | Every pre-registered number in this phase (α_conformal, the S1/S2/S3/S4 floors, the coverage floor of §C.3) is a candidate for the fifth. | Each bar must ship with (a) a derivation or an explicit "this is a judgement call" label, (b) the unit it is stated in, and (c) what would falsify it. The §C.3 floors are flagged `[ASSUMED]` above for exactly this reason. |
| 6 | **The gitignored 50 MB pool** (user memory `worktree-cleanup-destroys-gitignored-bulk-data`) | `h.meta.pool_dir` → `spike/data/cache/p13/7c65a1a9…` is 50 MB, gitignored, and does **not exist inside a fresh worktree**. `p13_real_id_reference` needs it to fit the OOD density null. A worktree-based execution would find it missing, or a cleanup would destroy it. | **Execute Phase 14 on the main working tree, no worktrees** — the same ruling Phases 12 and 13 took, for the same reason. State it in the plan so an orchestrator does not default to worktrees. |
| 7 | **The `12-15`/`12-17` pool-overlap hazard (CONVENTIONS C-04)** | Index-keyed generation means two pools at the same key share samples `1:min(n,m)` byte-identically. | Phase 13's arithmetic already defends against it by folding the counter into the **key word**, not the counter word (`spike/p13/datagen.jl:161-181` explains why). Phase 14 must copy that construction and assert the key words distinct — see Pitfall 2 below. |

### E.4 Ranked pitfalls, with detection

**Pitfall 1 — the calibration set must be π-distributed, NOT class-balanced.**
The Phase-13 gate set is stratified to equal thirds (`stratify_by_class`, realized E 1334 / R 1333 /
C 1333) while π at τ = 0.15 is E 0.3119 / R 0.4713 / C 0.2168
(`spike/p13/datagen.jl:559-563`, and `meta.pi_class_masses`). **Split conformal's guarantee requires
the calibration points to be exchangeable with the test point.** A class-balanced calibration set and
a π-distributed deployment stream are not exchangeable, and the realized coverage will miss `1−α` by
an amount driven by the class-conditional score distributions. This is silent: the code runs, the
quantile computes, the number looks fine.
*Fix:* draw the calibration set **unstratified** — `p13_simulate_indices(1:n, basis; rng_for = p14_cal_rng)`
(`spike/p13/datagen.jl:502`), not `stratify_by_class`. Same for the risk-coverage and FDR evaluation
sets, because FDR depends on prevalence.
*Detection:* assert the realized class masses of the calibration set are within Monte-Carlo error of
`pi_class_masses` (a χ² or a simple ±3 sd band on each), and **fail** if they are equal thirds.

**Pitfall 2 — Philox key-word collision between the P14 and P13 streams.**
The house construction is `Philox4x((SEED ⊻ SALT ⊻ COUNTER, idx))` (`spike/p13/datagen.jl:180-181`),
guarded by explicit inequality asserts at `:197-199`. If Phase 14 uses a new `P14_DEV_SEED` with the
same `P13_SALT`, then `P14_DEV_SEED ⊻ P13_SALT ⊻ c` can collide with `P13_DEV_SEED ⊻ P13_SALT ⊻ c'`
whenever `P14_DEV_SEED ⊻ P13_DEV_SEED == c ⊻ c'` — and the counters are small integers, so this is
not a negligible probability if the seeds are chosen by the usual hex-mnemonic pattern
(`P13_DEV_SEED = 0x0B13DE71`, `P13_FIX_SEED = 0x0B13F1F7`).
*Fix:* assert the **full cross-product** of key words pairwise distinct, over
`c ∈ {P13_TAU_COUNTER, P13_DATAGEN_COUNTER, P13_GATE_COUNTER, P13_ALPHA_COUNTER,
P13_CONTINUITY_COUNTER, P13_FIXTURE_COUNTER}` and the plain `p13_rng` key. Also assert P14 seeds ∉
`_p13_forbidden()` (`spike/p13/consts.jl:730-741`).

**Pitfall 3 — reusing the gate's 4,000 items as the conformal calibration set.**
Tempting (free, on disk, structurally disjoint from training since `P13_GATE_COUNTER ≠
P13_DATAGEN_COUNTER`, asserted at `run_three_way_gate.jl:358`). But it is (a) class-balanced
(Pitfall 1) and (b) **the very set the Phase-13 ECE gate was reported on** — using it to calibrate a
Phase-14 hedge means the hedge is calibrated on the set whose calibration was the headline claim.
*Recommendation:* use it freely as a **development substrate** (getting the code right, sanity plots,
fixture construction) and as a **reported cross-check**; draw fresh, unstratified sets at new P14
counters for every reported number. Cost ≈ 5-6 min per 2,000 items.

**Pitfall 4 — the shipped OOD null is on a DIFFERENT basis than the Phase-13 net.**
`artifacts/amended_v2/grid_8/ood_nulls_8.jld2` + `gate_report_8.jld2` carry the shipped grid-8
density null and its recorded `report.ood.id_threshold`. But the shipped bundle's `zt` standardizer
is the Phase-7 one, while the Phase-13 net rides **Phase 11's** `zt`, inherited and never re-fit
(`h.meta.zt_provenance = "INHERITED from the Phase-11 research NPE, never re-fit (D-02)"`,
`h.meta.phase11_sha`). Scoring a Phase-11-standardized summary against a Phase-7-fitted Mahalanobis
null produces a number with no interpretation. Phase 13 handled this by keeping the shipped bundle as
a **read-only comparison reference only** and building its own ID reference on its own basis
(`run_p13_realimage.jl` docstrings quoted in §D.4).
*Fix:* fit the Phase-14 density null on the Phase-13 net's own pool (`h.meta.pool_dir`) at
`q = OOD_ID_QUANTILE = 0.95`. If the shipped bundle is also read, label it a comparison reference and
never fuse the two thresholds.

**Pitfall 5 — unit mismatch between the classical scale and τ.**
τ = 0.15 is defined on **ρ_true**, measured by a probe that stops at `encode_d01(patch_summary(...))`
and never applies the standardization (`spike/p13/tau_probe.jl:68-73`). `patch_correlation` returns
the **induced mean per-patch correlation** μ. Applying τ directly to μ is a criterion applied in a
unit it was not derived for — the exact failure `P12_STAGE1_RATIO_CEILING` carries an explicit
recorded warning about, and one of the four "the bar was wrong" families.
*Fix:* `ghat(μ)` first (§D.5), and record `basis = :ghat_rho_true_scale` on the result so the choice
is auditable rather than implicit.
*Detection:* assert `ghat` is applied on every path that compares a classical statistic to τ — a
one-line grep test over the Phase-14 sources, in the style of `test_p13_tau.jl`'s grid-extension grep.

**Pitfall 6 — `Statistics.quantile` is the WRONG function for the conformal quantile.**
`quantile(s, 1-α)` uses linear interpolation between order statistics (type-7 by default). The
conformal guarantee is stated for the **order statistic** `sort(s)[⌈(n+1)(1−α)⌉]`. Interpolating
gives a slightly smaller `q̂` and breaks the finite-sample coverage bound. It fails *silently* — the
realized coverage is off by O(1/n), invisible at n = 2000 without a directed test.
*Fix:* write the order statistic explicitly. *Detection:* a small-n unit test — at n = 9, α = 0.1,
`⌈10·0.9⌉ = 9`, so `q̂ = max(s)`; `quantile(s, 0.9)` gives something strictly smaller. Assert the
exact index.

**Pitfall 7 — `costes_p` and the forbidden seed.**
`spike/comparator/config.jl:56` sets `MASTER_SEED = 0x00000000_00C0FFEE`, which **is**
`NPE_MASTER_SEED` in Phase 13's forbidden-seed inventory (`spike/p13/consts.jl:117`). `costes_p`
takes `master_seed` as its first positional argument (`classical.jl:181`), so the fix is trivial —
pass `P14_DEV_SEED`. But an executor copying the Phase-9 call site verbatim would consume a forbidden
seed's XOR-derived stream, and the seed-disjointness assertion in `test_p14_consts.jl` would not
catch it (it checks constants, not call sites).
*Detection:* a grep test asserting no Phase-14 source calls `costes_p(MASTER_SEED, …)`.

**Pitfall 8 — the corpus token ban is a code ban, not a prose ban.**
`test_p12_decoupling.jl:246-256` strips comments before scanning, so *prose about the holdout is
free and a CODE reference is not*, and `:214-218` explicitly requires the constants file to name
`CORPUS_MASTER_SEED` and the decoupling test itself to contain the tokens. An executor who tries to
satisfy the guard by removing every mention of the corpus will break the guard's own preconditions.
Copy the allowlist structure, do not invent a stricter one.

**Pitfall 9 — vacuous conformal sets from an over-confident `p̂`.**
If `q̂` comes out very small (a highly separable calibration set), `|C(Z)| = 1` almost always and the
hedge never fires — it looks like a pass but measures nothing. This is the same shape as the
`vacuous_pass` guard Phase 13 built for ECE (`spike/p13/result.jl:284-298`): *"a head that learns
nothing produces well-calibrated-looking probabilities"*.
*Fix:* report `q̂`, the realized set-size distribution, and the realized `|set| ≥ 2` / `|set| = 0`
rates on the calibration set **beside** the coverage number, and pre-register a "vacuity label"
(not a failure) if the ambiguous rate is below some floor. Follow `vacuous_pass`'s design: a label,
not a verdict.

**Pitfall 10 — the net is an epoch-4 checkpoint of an early-overfitting run.**
`13-REPORT.md` limit C and §2: best-validation risk 0.13625 at epoch 4, rising to 0.38632, with a
~10× train/validation gap, and the recipe was deliberately not tuned. Every Phase-14 number inherits
this. It does not invalidate anything (the gate was scored on a fresh disjoint set and passed 6/6),
but it must travel with the report and must not be discovered by a reviewer.

---

## Runtime State Inventory

Phase 14 is additive (new files under `spike/p14/`), not a rename or refactor. Recorded for
completeness because two categories are non-empty and matter for execution.

| Category | Items Found | Action Required |
|---|---|---|
| Stored data | `spike/data/cache/p13/7c65a1a9…` — the 50 MB, **gitignored** Phase-13 training pool that `h.meta.pool_dir` points at and that the OOD density null must be fitted on. Also `spike/data/cache/p13/57014f9e…` (456 KB, a fixture-scale pool). | **Read-only.** Never regenerate, never write into `cache/p13/`. Phase-14 caches (if any) go under a new `cache/p14/` root, mirroring the `p11`/`p12` separation asserted at `test_p12_decoupling.jl:321-349`. |
| Live service config | None — no external service. Verified: no daemon, no scheduler, no remote config in this repo. | None. |
| OS-registered state | None. Verified: no Task Scheduler / systemd / pm2 registration anywhere in the repo. | None. |
| Secrets / env vars | None. Verified: no `.env`, no SOPS, no CI secret referenced by the spike. | None. |
| Build artifacts | `artifacts/amended_v2/grid_8/` (present locally, gitignored, Release-hosted): `npe_8.jld2`, `ood_nulls_8.jld2`, `ratio_8.jld2`, `gate_report_8.jld2`. Phase 13 opened them **read-only** as a comparison reference. | **Read-only.** `test_p12_decoupling.jl:299-319` guards the shipped bundle; the executable half is skipped when `artifacts/` is untracked, so the guard is *weaker than it looks* — Phase 14 must not rely on it and must simply never write there. |

---

## Environment Availability

| Dependency | Required by | Available | Version / evidence | Fallback |
|---|---|---|---|---|
| Julia + `spike/Project.toml` env | everything | ✓ | 16 deps, NeuralEstimators pinned `v0.2.1` | — |
| `spike/p13/three_way_net.jld2` | the whole phase | ✓ | 953,027 bytes; loaded in this session | **None. Hard precondition — fail loudly.** |
| `spike/p13/tau_probe_report.jld2` | D-07 provenance | ✓ | 15,343 bytes; `tau = 0.15`, `tau_status = "measured"` | None. |
| `spike/p13/three_way_gate_report.jld2` | dev substrate + cross-check | ✓ | 425,996 bytes; 4,000 items with all columns | Not blocking (can redraw). |
| `spike/npe/p11_research_npe.jld2` | the frozen `zt` basis via `load_p13_basis()` | ✓ | 4,933,362 bytes, sha `659cda09…` matches `h.meta.phase11_sha` | **None. Hard precondition.** |
| `spike/data/cache/p13/7c65a1a9…` (50 MB pool) | OOD density null fit | ✓ | present, gitignored | Regenerable via `p13_generate_pool` but that is a large spend; **do not plan on it**. Execute on the main tree (Pitfall 6 above). |
| `artifacts/amended_v2/grid_8/*` | optional shipped comparison reference | ✓ | 4 files present locally | Degrade to "frozen constants quoted alone" — Phase 13's own fallback (`p13_real_shipped_reference`). |
| `corpus/data/*` — the 30 open rows' images | **D-03's reported check** | **✗** | `corpus/data/` = **0 entries**; all rows `bytes = 0`; `sha256` empty on the 30 open rows, `PENDING-FETCH` on the 2 sealed | **NO FALLBACK. Report the check as unexecutable, with the reason and the evidence.** A fetch is a human decision (licence + bandwidth) and out of Phase-14 scope. |
| `test/test_images/{positive,negative}/*_c{1,2,3}.tif` | n=2 real behaviour note | ✓ | 6 TIFFs, 2 specimens × 3 channels, **no coloc labels**, both OOD at ~4.2× | This is the only permitted real-image path (asserted positively at `test_p12_decoupling.jl:284-290`). |
| `git` | the decoupling guards | ✓ | guards degrade to `@info` skips without it (`test_p12_decoupling.jl:165,180,294`) | Guards skip, do not fail. |
| GPU / CUDA | — | not needed | `P13_USE_GPU = false`; `train_three_way` hard-throws on `use_gpu = true` | CPU-only is the baseline. |

**Missing dependency with no fallback:** the corpus images. This is the one item that changes what
Phase 14 can claim, and it must be surfaced in planning, not in execution.

**Compute budget reference:** the 4,000-item gate evaluation cost ~11 min CPU end-to-end
(`13-REPORT.md` §6d), with ~1.58× rejection overhead under stratification (which an *unstratified*
draw avoids entirely). So each 2,000-item unstratified Phase-14 set should cost ≲ 4 min. Three sets
(calibration, evaluation, OOD-arm) ≲ 15 min total. Training time is zero — nothing is trained.

---

## Code Examples

### The three-class posterior and the FDR rule (the two core primitives)

```julia
"""
    p14_class_posterior(logbf, prior) -> NamedTuple

Softmax over (log prior + log BF against the `random` reference). `logbf` must carry the keys
:coloc, :random, :exclusion with `random === 0.0` (D-08). `prior` is the MEASURED
`h.meta.pi_class_masses`, never a chosen value.
"""
function p14_class_posterior(logbf, prior)
    logbf.random === 0.0 || throw(ArgumentError(
        "p14_class_posterior: logbf.random must be exactly 0.0 (D-08 structural zero); got $(logbf.random)"))
    a = (log(prior.coloc) + logbf.coloc, log(prior.random), log(prior.exclusion) + logbf.exclusion)
    m = maximum(a)                       # log-sum-exp: NEVER exp() a raw ±22-nat logbf
    w = exp.(a .- m); s = sum(w)
    p = (coloc = w[1]/s, random = w[2]/s, exclusion = w[3]/s)
    @assert abs(p.coloc + p.random + p.exclusion - 1.0) < 1e-12
    return p
end

"""
    p14_bayes_fdr(v, alpha) -> NamedTuple

The direct posterior-probability Bayesian-FDR rule (Newton et al. 2004; Mueller et al. 2004/2007).
`v[i] = P(H0 | Z_i)`. Accept the LARGEST PREFIX of the ascending order whose RUNNING MEAN stays at
or below `alpha`. Returns the accepted ORIGINAL indices, k*, the realized posterior expected FDP,
the induced threshold t* and its implied decision-risk cost ratio.

RUNNING MEAN, NOT RUNNING MAX -- the controlled quantity E[FDP | data] is a mean by construction.
The running mean of an ascending sequence is non-decreasing, so the accepted set is genuinely a
prefix and the rule reduces to a threshold on v.
"""
function p14_bayes_fdr(v::AbstractVector{<:Real}, alpha::Real)
    0.0 < alpha < 1.0 || throw(ArgumentError("p14_bayes_fdr: alpha must be in (0,1), got $alpha"))
    ord  = sortperm(v)                      # ascending
    sv   = v[ord]
    run  = cumsum(sv) ./ (1:length(sv))     # running MEAN
    @assert issorted(run)                   # the prefix property, asserted not assumed
    k    = findlast(<=(alpha), run)
    k === nothing && return (accepted = Int[], k_star = 0, fdp = 0.0,
                             t_star = NaN, cost_ratio = NaN, n = length(v))
    t = sv[k]
    return (accepted = ord[1:k], k_star = k, fdp = run[k],
            t_star = t, cost_ratio = t / (1 - t), n = length(v))
end
```

### Hand-rolled split conformal (LAC), zero dependencies

```julia
const P14_CLASS_KEYS = (:coloc, :random, :exclusion)   # by NAME, never positional (Pitfall/landmine 2)

"LAC / THR nonconformity score: s(Z, y) = 1 - p_hat(y | Z)."
p14_lac_score(p::NamedTuple, y::Symbol) = 1.0 - getproperty(p, y)

"""
    p14_conformal_quantile(scores, alpha) -> Float64

The FINITE-SAMPLE conformal quantile: the ceil((n+1)(1-alpha))-th ORDER STATISTIC.

DO NOT USE `Statistics.quantile`. Its default type-7 linear interpolation between order statistics
returns a strictly smaller value and silently breaks the >= 1-alpha coverage bound by O(1/n).
Requires alpha >= 1/(n+1); at n = 30 that is alpha >= 1/31 ~ 0.0323 (the D-03 corpus floor).
"""
function p14_conformal_quantile(scores::AbstractVector{<:Real}, alpha::Real)
    n = length(scores)
    k = ceil(Int, (n + 1) * (1 - alpha))
    k <= n || throw(ArgumentError(
        "p14_conformal_quantile: alpha = $alpha is infeasible at n = $n; " *
        "split conformal requires alpha >= 1/(n+1) = $(1/(n+1))"))
    return sort(collect(Float64, scores))[k]
end

"""
    p14_conformal_set(p, qhat) -> (set::Tuple{Vararg{Symbol}}, status::Symbol)

`status` is :singleton (DECIDE), :ambiguous (|set| >= 2 -> ABSTAIN) or :empty
(|set| == 0 -> ABSTAIN). The EMPTY case is kept separate on purpose: it is a distribution-free
misspecification signal arriving through a different channel from the Mahalanobis OOD detector,
and collapsing it into :ambiguous would discard that.
"""
function p14_conformal_set(p::NamedTuple, qhat::Real)
    set = Tuple(y for y in P14_CLASS_KEYS if p14_lac_score(p, y) <= qhat)
    st  = length(set) == 1 ? :singleton : (isempty(set) ? :empty : :ambiguous)
    return (set = set, status = st)
end
```

### The D-05 asymmetric fusion, with D-06's three-valued OOD

```julia
"""
    p14_fuse(; ood_state, conformal_status, disagree, allow_unchecked_ood = false) -> (Symbol, Symbol)

D-05 ASYMMETRIC fusion (NOT the literal OR) plus D-06's not-checked default.
Returns (action, reason) with action in (:decide, :abstain).

Order is load-bearing and each branch cites its decision:
  * OOD fired                          -> ABSTAIN                       (D-05)
  * OOD not-checked, no opt-out        -> ABSTAIN                       (D-06)
  * conformal set ambiguous or empty   -> ABSTAIN                       (D-05)
  * disagreement AND OOD not clear     -> ABSTAIN                       (D-05)
  * disagreement ALONE                 -> DECIDE, recorded as a field   (D-05: the Phase-9 thesis)
"""
function p14_fuse(; ood_state::Symbol, conformal_status::Symbol, disagree::Bool,
                  allow_unchecked_ood::Bool = false)
    ood_state in (:fired, :clear, :not_checked) || throw(ArgumentError(
        "p14_fuse: ood_state must be :fired, :clear or :not_checked (D-06); got $ood_state"))
    ood_state === :fired && return (:abstain, :ood_fired)
    if ood_state === :not_checked && !allow_unchecked_ood
        return (:abstain, :ood_not_checked)      # D-06: a false flag that means NOT CHECKED
    end
    conformal_status === :ambiguous && return (:abstain, :conformal_ambiguous)
    conformal_status === :empty     && return (:abstain, :conformal_empty)
    if disagree && ood_state !== :clear
        return (:abstain, :disagreement_and_ood) # only reachable under the explicit opt-out
    end
    return (:decide, disagree ? :decided_contra_classical : :decided)
end
```

---

## State of the Art

| Old approach | Current approach | When changed | Impact here |
|---|---|---|---|
| Frequentist BH FDR on p-values | Direct posterior-probability Bayesian FDR on `P(H₀\|data)` | Newton et al. 2004 / Müller et al. 2004 | The right family for a Bayesian posterior; and it is the reason pre-filtering by abstention does not require a selective-inference correction (§E.2). |
| Full (transductive) conformal | Split / inductive conformal | Papadopoulos et al. 2002; Lei et al. 2018 | Makes the hedge a single sorted-array lookup, which is why zero new dependencies is achievable. |
| Softmax-threshold reject option (Chow) | Selective prediction with risk-coverage / AURC | El-Yaniv & Wiener 2010; Geifman & El-Yaniv 2017; Geifman et al. 2019 | Gives SC3 its standard object and its standard normalization (E-AURC / selective skill). |
| APS as the default conformal classification score | LAC/THR when the classifier is well calibrated; APS when conditional coverage matters | Sadinle et al. 2019 vs Romano et al. 2020 | The ECE-green evidence here selects LAC (§B.2). |

**Deprecated / not applicable in this repo:**
- `compute_BayesFactor()` / KDE+quadgk as any kind of reference — declined by Phase 13, invalid past
  \|log BF\| ≈ 6.9, and both packages are absent from the direct deps by assertion.
- `youden_j` as an operating point — present in `src/amortized/ood.jl:310` but explicitly *"A LABELED
  REFERENCE only — NEVER the gate"*. If any Phase-14 threshold is tempted toward it, that is a
  post-hoc bar and forbidden.

---

## Validation Architecture

### Test framework

| Property | Value |
|---|---|
| Framework | Julia stdlib `Test` (`@testset` / `@test`) |
| Config file | none — plain scripts under `spike/test/` |
| Quick run command | `julia --project=spike spike/test/test_p14_<name>.jl` |
| Full suite command | `julia --project=spike spike/test/runtests.jl` — **known to exit 1 at the Phase-4 `SPEEDUP_GATE`, masking every later include block** (`deferred-items.md` D-13-A, `13-REPORT.md` §11 item 6). Per-file runs are the only reliable signal. |
| Fixture stream | `p13_fix_rng(P13_FIXTURE_COUNTER)` pattern — Phase 14 needs its own `P14_FIX_SEED` so a quick gate never pre-observes a reported stream |

### SC → validation map

| SC | Behaviour | Measured on | Pre-registered number | Falsified if | Test type | Command |
|---|---|---|---|---|---|---|
| **SC1-a** (FDR rule correctness) | `p14_bayes_fdr` implements the running-**mean** prefix rule and reduces to a threshold | hand-built fixtures, no simulation | `v = [0.01, 0.02, 0.15]`, α = 0.10 ⇒ `k* = 3`, `fdp = 0.06`; `issorted(run)` holds on 10⁴ random ascending vectors | any discriminating fixture returns the running-max answer, or `run` is not sorted | unit | `julia --project=spike spike/test/test_p14_fdr.jl` |
| **SC1-b** (realized FDP) | On a π-distributed simulator batch with known labels, the realized false-discovery proportion among DECIDED-and-accepted items tracks α_FDR | fresh unstratified draws, `n_eval = 2000`, at `P14_EVAL_COUNTER` | over α_FDR ∈ {0.01, 0.05, 0.10, 0.20}: realized FDP within the binomial 95 % band of the *predicted* `fdp` returned by the rule | realized FDP exceeds the upper band at any α on the grid | integration | `julia --project=spike spike/p14/run_p14_fdr_check.jl` |
| **SC1-c** (prior sensitivity, named limit) | The FDR estimate's dependence on the assumed class prior is quantified, not hidden | same batch, π_C swept over the pre-registered grid {0.05, 0.10, 0.217, 0.40, 0.70} | reported curve; **no bar** — this is a REPORTED, NOT GATED quantity | n/a (reporting obligation; failing to report it is the failure) | integration | same runner |
| **SC1-d** (conformal coverage) | Split-conformal marginal coverage ≥ 1 − α_conf | calibrate on `n_cal = 2000` at `P14_CAL_COUNTER`; measure on `n_eval = 2000` at `P14_EVAL_COUNTER`; both **unstratified** | realized coverage ≥ 1 − α_conf − 3·sqrt(α(1−α)/n_eval); with α_conf = 0.10, n = 2000 ⇒ ≥ 0.880 | realized coverage falls below the band | integration | `julia --project=spike spike/p14/run_p14_conformal.jl` |
| **SC1-e** (conformal quantile arithmetic) | `p14_conformal_quantile` returns the **order statistic**, not an interpolated quantile | fixture | at n = 9, α = 0.10 ⇒ `q̂ == maximum(s)`, and `q̂ != Statistics.quantile(s, 0.9)` for a non-degenerate `s`; at n = 30 the α = 0.03 call **throws** | either equality flips | unit | `test_p14_conformal.jl` |
| **SC1-f** (D-03 corpus check) | The corpus-side coverage check is attempted and its outcome recorded | `corpus/manifest.csv` metadata only | `corpus_images_available == 0` and the check is recorded as **UNEXECUTABLE with reason**; no coverage number is emitted from corpus data | a coverage number is emitted from corpus data, or the check is silently omitted | integration | `test_p14_decoupling.jl` + the report |
| **SC2-a** (fusion truth table) | `p14_fuse` implements D-05 asymmetry exactly | fixture — all 3 × 3 × 2 × 2 = 36 input combinations | full truth table asserted cell by cell; in particular `(:clear, :singleton, true, false) ⇒ (:decide, :decided_contra_classical)` | any cell differs | unit | `test_p14_fuse.jl` |
| **SC2-b** (D-06 not-checked) | A missing or non-finite `:thr` yields `:not_checked` and abstains by default | fixture `ood_nulls` with `:thr` absent / `Inf` / `NaN` / a finite value | all three non-finite/absent cases ⇒ `:not_checked` ⇒ `:abstain` with `reason = :ood_not_checked`; opt-out flips to `:decide` | any case reads as `:clear` | unit | `test_p14_fuse.jl` |
| **SC2-c** (D-06 anti-regression) | No Phase-14 source uses `verdict.flag == false` as a proxy for in-distribution | grep of comment-stripped Phase-14 sources | zero hits for `!.*\.flag` / `\.flag == false` outside `p14_ood_state` | any hit | unit (source grep) | `test_p14_decoupling.jl` |
| **SC2-d** (cross-method scale) | Every classical-vs-amortized comparison passes through `ghat` | grep of comment-stripped Phase-14 sources | every `patch_correlation(` call site compared to τ is wrapped in `ghat(` | any bare comparison | unit (source grep) | `test_p14_decoupling.jl` |
| **SC3-a** (selective skill) | Confidence ordering beats random | `n_eval = 2000` unstratified | `skill = (AURC_random − AURC)/(AURC_random − AURC_oracle) ≥ 0.60` | skill below floor | integration | `run_p14_riskcoverage.jl` |
| **SC3-b** (monotone trend) | Selective risk rises with coverage on the declared range | same | `spearman(coverage, selective_risk) ≥ 0.95` over coverage ∈ **[0.20, 1.00]**, floor frozen before the run | ρ below 0.95, **or** the floor is changed after seeing the curve | integration | same runner |
| **SC3-c** (concentrates on hard) | Abstention score discriminates would-be errors | same | `AUC_hard ≥ 0.80` via the in-repo `roc_auc` | AUC below floor | integration | same runner |
| **SC3-d** (concentrates on OOD) | Abstention rate rises under injected misspecification | `OOD_FAMILIES` at the strongest level vs in-distribution, matched n | `abstain_rate(misspec) − abstain_rate(ID) ≥ 0.50` | difference below floor | integration | `run_p14_ood_arm.jl` |
| **D-07** (τ provenance) | τ is loaded, four-way agreed, and sha-pinned | the four artifacts | `net.tau == probe.tau == gate.P13_TAU == P13_TAU == 0.15`; `net.consts_sha == P13_CONSTS_SHA.pre_amendment`; `p13_consts_sha() == .post_amendment`; **`probe.consts_sha256 ∉ values(P13_CONSTS_SHA)`** (assert the divergence, §D.3) | any equality fails, or the divergence assertion is widened to a disjunction | unit | `test_p14_provenance.jl` |
| **D-01/env** | `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `corpus/`, `artifacts/` untouched; 16 deps exactly; no CUDA | git + `Pkg` | all assertions from `test_p12_decoupling.jl:152-208, 291-296, 351-353`, re-stated in the Phase-14 file | any diff | unit | `test_p14_decoupling.jl` |
| **Seed discipline** | P14 streams provably disjoint from every P13 and forbidden stream | constants | P14 seeds ∉ `_p13_forbidden()`; P14 counters pairwise distinct; **full cross-product of Philox key words pairwise distinct** vs the P13 family (Pitfall 2) | any collision | unit | `test_p14_consts.jl` |
| **Class-order anti-permutation** | Every named field lands on its intended class | fixture with deliberately UNEQUAL posteriors (0.7 / 0.2 / 0.1) | each of `p.coloc`, `p.random`, `p.exclusion` equals its intended value | any permutation | unit | `test_p14_posterior.jl` |

### Sampling rate

- **Per task commit:** the relevant `test_p14_*.jl` unit file(s) — each < 10 s, no simulation.
- **Per wave merge:** all `test_p14_*.jl` unit files, run per-file (never via `runtests.jl`).
- **Phase gate:** all unit files green **plus** the three integration runners
  (`run_p14_conformal.jl`, `run_p14_riskcoverage.jl`, `run_p14_ood_arm.jl`) producing persisted
  `.jld2` artifacts with recorded seeds, counters and sha provenance, before `/bm:verify-work`.

### Wave 0 gaps

- [ ] `spike/p14/consts.jl` — the Tier-1 pre-registration: `P14_DEV_SEED`, `P14_FIX_SEED`,
      `P14_CAL_COUNTER`, `P14_EVAL_COUNTER`, `P14_OOD_COUNTER`, `P14_FIXTURE_COUNTER`,
      `P14_ALPHA_CONFORMAL`, `P14_N_CAL`, `P14_N_EVAL`, `P14_SKILL_FLOOR`,
      `P14_SPEARMAN_FLOOR`, `P14_COVERAGE_FLOOR`, `P14_AUC_HARD_FLOOR`, `P14_OOD_MARGIN_FLOOR`,
      `P14_ITERATION_ALLOWANCE`. **Frozen and committed in a separate commit before any result
      exists**, with the Tier-1/Tier-2 guard-block structure of `spike/p13/consts.jl` and the
      append-never-overwrite rule.
- [ ] `spike/test/test_p14_consts.jl` — literal-value + disjointness assertions (mirrors
      `test_p13_consts.jl`, 137 tests).
- [ ] `spike/test/test_p14_decoupling.jl` — the environment/corpus/src guards (mirrors
      `test_p12_decoupling.jl`, with the Phase-14 allowlist and the two new source-grep checks).
- [ ] `spike/test/test_p14_provenance.jl` — the D-07 four-way τ + sha table.
- [ ] `spike/test/fixtures/` — the unequal-posterior class-order fixture and the FDR discriminating
      fixture.
- [ ] No framework install needed (`Test` is stdlib).

---

## Assumptions Log

| # | Claim | Section | Risk if wrong |
|---|---|---|---|
| A1 | The population risk-coverage curve is non-decreasing in coverage **iff** the confidence score is monotone in conditional risk | §C.2 | Low. It is a two-line derivation; if wrong, the S2 statistic's *rationale* changes but the statistic itself (a Spearman bar on a declared range) still measures the right thing. |
| A2 | Bayesian posterior-expected FDP is immune to data-dependent pre-selection because `R` is σ(data)-measurable | §E.2 | **High if wrong** — it is the whole justification for abstain-then-sort and for the "controlled over the decided subset" claim. The argument is standard conditioning, but it should be stated as a derivation in the report so a referee can check it rather than take it on trust. Worth an explicit sentence in the plan asking the user to sanity-check. |
| A3 | `skill ≥ 0.60`, `spearman ≥ 0.95`, `AUC_hard ≥ 0.80`, `ood margin ≥ 0.50`, coverage floor `0.20` | §C.3, Validation | **Judgement calls with no derivation.** This project has been burned four times by an underived bar. Each must be frozen before any result and labelled in the report as a judgement call, or replaced with a user-ruled number. **Recommend the planner surface these four numbers to the user before freezing.** |
| A4 | `n_cal = n_eval = 2000` is the right size | §B.4 | Low. The Beta-distribution sd arithmetic is exact; 2000 is a cost/precision choice. Larger is strictly better and cheap. |
| A5 | The simulator's per-item cost is ≈ the gate's (~11 min / 4000 items, minus the 1.58× stratification overhead) | Environment Availability | Low. If wrong, only the compute estimate moves; nothing structural. |
| A6 | Both fixtures in `test/test_images/` will produce ABSTAIN under D-05 | §C.4 | Low, and it is a *prediction to be tested*, not an input. If they do not abstain despite 4.198× OOD scores, that is a finding about the OOD wiring worth chasing. |
| A7 | The `pi_class_masses` in `h.meta` were computed under τ = 0.15 with the D-05 two-factor cut and the post-Phase-11 prior | §A.2 | Medium. It matches `datagen.jl:559-563`'s stated table for τ = 0.15 (E 0.3119 / R 0.4713 / C 0.2168) to 4 decimals, which is strong corroboration — but the planner should add an executable re-derivation from `sample_prior` + `three_way_label` on a fixture stream as a cheap confirmation, since the FDR guarantee rests on it. |

---

## Open Questions (RESOLVED — resolutions appended 2026-08-03 during `/bm:plan-phase 14`)

> **Resolution index.** All five questions below were carried into planning and resolved there. The
> original question text is left unedited (corrections are appended, never overwritten); each
> resolution names the plan that carries it.
>
> | Q | Resolution | Carried by |
> |---|---|---|
> | Q1 — are the SC3 bars acceptable, or should the user rule them? | **Escalated, not decided.** `14-01` opens with a **blocking** decision checkpoint (`autonomous: false`) offering: accept as proposed / user sets them / freeze as REPORTED-NOT-GATED. They may be ruled **before** the freeze commit and never after. | `14-01-PLAN.md` |
> | Q2 — corpus: report the absence, or fetch? | **Neither, as originally framed.** Superseded by **D-03a**: the corpus holds no unsealed physical truth (30/32 rows are `simulated-secondary`; the 2 `physical-primary` rows are Phase 16's sealed holdout) and its bytes are unfetched *by design*. The real-data check uses the **six committed microscopy TIFFs**, reported as an illustration, not a coverage claim. **No fetch is planned.** The `α ≥ 1/31 ≈ 0.032` figure is withdrawn. | `14-13-PLAN.md`, `14-CONTEXT.md` D-03a |
> | Q3 — default the class prior, or force the caller to supply one? | **Default to the measured `h.meta.pi_class_masses`** (measured, not chosen), exposed as a keyword, with the prior-sensitivity curve a **required** output rather than an option. | `14-04`, `14-06`, `14-07` |
> | Q4 — whose `CalibrationMeta` does the Phase-14 result carry? | **Carry Phase 13's forward unchanged**, and record Phase-14 conformal coverage as its own separate field, so no second differently-sourced ECE can be conflated with the gated one. | `14-07-PLAN.md` |
> | Q5 — is a `P14_ITERATION_ALLOWANCE` needed, and what is its trigger? | **Yes — `P14_ITERATION_ALLOWANCE = 1`** with exactly one pre-declared trigger, frozen in `spike/p14/consts.jl` and explicitly never spendable on relaxing a bar after a result. | `14-01-PLAN.md` |

1. **Are the four SC3 bars (A3) acceptable to the user, or should they be ruled?**
   - What we know: this project treats every bar as a pre-registration matter and has recorded four
     separate cases of a bar measuring something other than what it named.
   - What's unclear: whether the user wants to set these numbers personally.
   - **Recommendation:** surface all four (plus the 0.20 coverage floor) to the user during
     `/bm:plan-phase` review, freeze them in `spike/p14/consts.jl` in a commit that precedes every
     result commit, and record the freeze ordering the way `consts.jl:842-846` records τ's.

2. **What should `decide_coloc` do about the 30 corpus rows — report the absence, or attempt a fetch?**
   - What we know: `corpus/data/` is empty, `bytes = 0`, `corpus_images_available = 0`; the rows are
     single frames, not pairs; a fetch is a licence + bandwidth decision explicitly flagged as a human
     one in `13-REPORT.md` limit A.
   - **Recommendation:** report the absence with evidence, run the n=2 fixture behaviour note, and
     carry "the D-03 corpus bound is absent, not loose" as a named limit. Do **not** plan a fetch task.

3. **Should the class prior default to `pi_class_masses`, or should the caller be forced to supply one?**
   - What we know: the FDR number is prevalence-sensitive (§E.1), and the measured masses describe the
     *simulator*, not any real batch.
   - **Recommendation:** default to `h.meta.pi_class_masses` (so the default is measured, not chosen),
     make it a keyword, and make the prior-sensitivity curve a **required** output rather than an
     option — the same structural trick `p13_calibration_meta` uses for `auc`.

4. **Does the Phase-14 result carry a `CalibrationMeta`, and if so, whose?**
   - What we know: `ThreeHypothesisColocResult` requires one; `p13_calibration_meta` builds it from a
     spike `CalibrationResult` plus a required `auc` (`spike/p13/result.jl:337`).
   - What's unclear: whether Phase 14 should re-run `_bin_calibration` on its own draws or carry
     Phase 13's gate calibration forward.
   - **Recommendation:** carry Phase 13's forward *unchanged* (it is the calibration the decision layer
     relies on) and add a **separate** Phase-14 conformal-coverage record as its own field. Re-running
     `_bin_calibration` on Phase-14 draws would produce a second, differently-sourced ECE that a reader
     would conflate with the gated one.

5. **Is a `P14_ITERATION_ALLOWANCE` needed, and what would its single pre-declared trigger be?**
   - What we know: Phase 13 declared one, defined its trigger mechanically, and did not spend it
     (`three_way_gate_report.jld2` → `iteration_trigger_fired = false`).
   - **Recommendation:** declare `P14_ITERATION_ALLOWANCE = 1` with exactly one pre-declared trigger —
     e.g. *"realized conformal coverage falls below the SC1-d band, in which case the allowance is spent
     on switching the score from LAC to APS and re-running ONCE"* — and state in the constants file that
     it is never spent on relaxing a bar after a result. This is the house pattern
     (`consts.jl` iteration_trigger_text) and it pre-empts the exact argument that would otherwise arise.

---

## Sources

### Primary (HIGH confidence — read on disk in this session)

- `spike/Project.toml` — the exact 16 dependencies; NeuralEstimators pinned v0.2.1.
- `spike/p13/net.jl` — `ThreeWayEvidenceNet` (:157), `three_way_input_dim` (:132),
  `masked_two_head_bce` (:239), `three_way_log_bf` (:349), `three_way_probs` (:372),
  `save_three_way`/`load_three_way` (:459/:496), `p13_consts_sha` (:524), `P13_CONSTS_SHA` (:556).
- `spike/p13/result.jl` — `ThreeHypothesisColocResult` (:142), inner constructor D-08 guard (:159),
  `_p13_as_three_way_logbf` (:179), the four accessors (:194-243), `log_bf_vs_random` (:257),
  `vacuous_pass` (:297), `p13_calibration_meta` (:337).
- `spike/p13/labels.jl` — `three_way_label` (:121), `head_targets` (:153), `class_masses` (:301).
- `spike/p13/tau_probe.jl` — the D-06 probe, `mbar_from_summary` (:135), `tau_from_curve` (:548),
  the no-grid-extension rule (:526-537).
- `spike/p13/consts.jl` — seeds (:117-122, :197-198), counters (:213-218), `p13_tau()` (:272),
  `P13_TAU_*` Tier-1 (:289-331), `P13_GATE_M` (:419), disjointness asserts (:730-741),
  the Tier-2 τ block (:805-872).
- `spike/p13/datagen.jl` — `p13_datagen_rng` (:180) + key-word asserts (:197-199),
  `p13_simulate_indices` (:502), `stratify_by_class` (:583) and the class-mass table (:559-563).
- `spike/p13/run_three_way_gate.jl` — `p13_gate_rng` (:155), the stream-hygiene and sha assertions
  (:358-370), `P13_GATE_CLASS_ORDER` (:131).
- `spike/p13/preconditions.jl` — `P13_PHASE11_NET` (:213), `p13_require_phase11` (:275),
  `load_p13_basis` (:411).
- `spike/p13/run_p13_realimage.jl` — `p13_real_shipped_reference`, `p13_real_id_reference`, the
  copied `OOD_ID_QUANTILE = 0.95` and ridge.
- `src/results.jl` — `AbstractColocResult` (:43), `_iface_error` (:47), the four accessors (:51-79),
  `OODVerdict` (:93), `CalibrationMeta` (:117), the Phase-13/12 sketch block (:161-192).
- `src/amortized/ood.jl` — every signature listed in §D.4, incl. the `thr === nothing ⇒ flag = false`
  short-circuit (:394-396).
- `src/amortized/local_map.jl` — `LOCAL_MAP_SENTINEL` (:53), `LocalColocMap` + the NOT-CHECKED prose
  (:55-79), `_recorded_ood_threshold` (:93-108), the vacuous-flag warning (:110-119).
- `spike/comparator/classical.jl` — `manders` (:66), `pearson_whole` (:85), `spearman_whole` (:95),
  `patch_correlation` (:108), `costes_rng` (:135), `costes_p` (:181).
- `spike/comparator/config.jl` — `COSTES_*` (:37-43), `DIVERGENCE_*` (:50-51), `MASTER_SEED` (:56).
- `spike/simulator/ghat.jl` — `ghat` (:110), knots (:93,:95), monotonicity record (:44-53).
- `spike/validation/sbc.jl` — `CalibrationResult` (:65-72), `_bin_calibration` (:82).
- `spike/validation/ood.jl` — `OOD_FAMILIES` (:499), `ood_roc_over_grid` (:612), `roc_auc` (:319).
- `spike/test/test_p12_decoupling.jl` — the guards at :152-208, :210-297, :299-319, :321-349, :351-353.
- `corpus/manifest.csv` (all 37 lines), `corpus/anchor_rows.jl` (:85 `PENDING_FETCH`, :201/:217
  `split = "sealed_holdout"`, :268-276 the 30+2 = 32 contract), and `ls corpus/data` → **0 entries**.
- **Loaded artifacts (live `JLD2.load` in this session):** `spike/p13/three_way_net.jld2`,
  `spike/p13/tau_probe_report.jld2`, `spike/p13/three_way_gate_report.jld2` — every numeric value
  quoted in §A.1, §D.1 and §D.3 comes from these loads, not from a document.
- `.planning/phases/13-*/13-REPORT.md` (§2, §5, §6, §6b-d, §9 limits A-E, §11), `13-VERIFICATION.md`,
  `.planning/phases/14-*/14-CONTEXT.md`, `.planning/ROADMAP.md:374-383`, `.planning/STATE.md`,
  `.planning/REQUIREMENTS.md`, `.claude/skills/spike-findings-proteincoloc/SKILL.md`, `CLAUDE.md`.

### Secondary (MEDIUM confidence — method literature, verified by web search this session)

- Newton, Noueiry, Sarkar & Ahlquist (2004), *Biostatistics* — the direct posterior-probability FDR
  approach. https://web.ma.utexas.edu/users/pmueller/pap/MPR07.pdf (Müller, Parmigiani & Rice 2007,
  "FDR and Bayesian Multiple Comparisons Rules") states the decision-theoretic equivalence: the
  optimal rule under (FD, FN) losses is a threshold on the marginal posterior probability of the null.
- https://cran.r-project.org/package=bayefdr — an independent implementation confirming the
  sort-then-running-average form is the standard operationalization.
- Sadinle, Lei & Wasserman (2019), *JASA*, "Least Ambiguous Set-Valued Classifiers with Bounded Error
  Levels" — the LAC score `1 − p̂(y|x)` and its minimal-average-set-size optimality.
- https://www.emergentmind.com/topics/split-conformal-prediction-cp and
  https://par.nsf.gov/servlets/purl/10292056 (Romano, Sesia & Candès, "Classification with Valid and
  Adaptive Coverage") — the `⌈(n+1)(1−α)⌉/n` quantile, the ≥ 1−α marginal guarantee, and APS as the
  adaptive alternative.
- Vovk (2012), "Conditional validity of inductive conformal predictors" — the Beta distribution of
  realized coverage, used for the §B.4 sample-size table.
- El-Yaniv & Wiener (JMLR 2010); Geifman & El-Yaniv (NeurIPS 2017); Geifman, Uziel & El-Yaniv
  (ICLR 2019) — risk-coverage, AURC, E-AURC.

### Tertiary (LOW confidence — flagged for validation)

- The four SC3 bars and the 0.20 coverage floor (§C.3). No literature sets these; they are judgement
  calls and are labelled as such in the Assumptions Log (A3). **Do not present them as derived.**

---

## Metadata

**Confidence breakdown:**

- **Wiring / on-disk facts:** HIGH — every path, signature, line number, constant and artifact value
  was read or loaded in this session, not recalled. The four τ locations, the `head_log_odds`, the
  `pi_class_masses`, the sha triple and the empty `corpus/data/` were all verified directly.
- **Bayesian FDR method:** HIGH — the running-mean form and the "controls posterior expected FDP"
  reading were cross-verified against the Müller/Parmigiani/Rice PDF and an independent CRAN
  implementation.
- **Split conformal:** HIGH — the quantile formula, the coverage bound, the `α ≥ 1/(n+1)` feasibility
  condition (which independently reproduces D-03's 0.032) and the LAC/APS tradeoff are all confirmed.
- **Risk-coverage / SC3:** MEDIUM — the standard construction and AURC/E-AURC are well established;
  the *specific* bars are judgement calls (A3) and the monotonicity "iff" is a derivation, not a
  citation.
- **Pitfalls:** HIGH for those grounded in a file or an artifact (1-8, 10); MEDIUM for the
  prevalence-sensitivity claim (E.1), which is a modelling argument rather than a measurement.
- **The abstain-then-sort argument (E.2 / A2):** MEDIUM-HIGH — the conditioning argument is standard
  and short, but it carries the phase's central honesty claim and should be independently checked.

**Research date:** 2026-08-03
**Valid until:** the on-disk facts are valid until any Phase-13 artifact or `spike/Project.toml`
changes (both are guarded by running tests, so change is loud). The method sections are stable
literature and do not expire. **Re-verify `corpus/data/` before writing any Phase-14 report** — it is
the one input whose state could change by a human action outside this phase.
