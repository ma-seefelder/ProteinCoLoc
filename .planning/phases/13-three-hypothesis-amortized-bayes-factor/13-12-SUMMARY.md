---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 12
status: complete
subsystem: amortized-inference
tags: [reported-gate, discrimination, calibration, pre-registration, three-way-evidence]
requires:
  - "spike/p13/three_way_net.jld2 (13-11) — the trained two-head evidence net + its MEASURED per-head log-odds"
  - "spike/p13/consts.jl (13-01/13-10) — the byte-locked pre-registration incl. MEASURED P13_TAU = 0.15"
  - "spike/p13/datagen.jl (13-11) — stratify_by_class / assemble_conditioned_pairs / p13_encode_pair"
  - "spike/p13/result.jl (13-08) — p13_calibration_meta / p13_traffic_light / vacuous_pass"
  - "spike/validation/sbc.jl — _bin_calibration (shared, gate-lineage, NOT modified)"
  - "spike/validation/ood.jl — the hand-rolled tie-aware roc_auc"
  - "spike/npe/p11_research_npe.jld2 (Phase 11) — the bound research basis"
provides:
  - "spike/p13/run_three_way_gate.jl — the REPORTED D-12/D-13 gate runner (persist-before-gate)"
  - "spike/p13/three_way_gate_report.jld2 — the reported artifact: every statistic, every threshold, the seed, the counter and the consts fingerprint"
  - "spike/figures/p13_confusion.png — the confusion/ROC figure (script artifact, gitignored under the house *.png rule)"
affects:
  - "13-13, 13-14, 13-16 — the phase report cites this verdict; the per-item scores are persisted so no plan needs to re-run the gate"
tech-stack:
  added: []
  patterns:
    - "Persist-before-gate: the artifact is written, and the figure rendered, before any assertion can throw"
    - "Reported evaluation stream at a DISJOINT counter: the datagen key construction at P13_GATE_COUNTER"
    - "Exact round trip between two frozen z-score bases via StatsBase.reconstruct, for a reported-not-gated cross-net comparison"
key-files:
  created:
    - spike/p13/run_three_way_gate.jl
    - spike/p13/three_way_gate_report.jld2
    - spike/figures/p13_confusion.png
  modified: []
decisions:
  - "The figure is rendered BEFORE the gate assertions, not after as the plan's step 14 ordered, because a failing @testset throws and an after-the-gate figure would exist only for a passing run"
  - "Continuity reads the spike-lane frozen binary NRE, because src/amortized/bf.jl imports KernelDensity/QuadGK (absent from the lean spike env) and consts.jl section I2 reserves the shipped grid-8 bundle for the OOD comparison only"
metrics:
  duration: ~50 min wall clock (of which 10.86 min was the reported evaluation-set generation)
  completed: 2026-07-29
---

# Phase 13 Plan 12: Reported Three-Way Gate Summary

**VERDICT: PASS.** The amended D-12/D-13 gate ran **once**, on a fresh 4,000-pair evaluation set
drawn at a counter disjoint from the training pool, and cleared **all six** pre-registered criteria:
per-class one-vs-random AUC **0.990169** (coloc) and **0.988262** (exclusion) against the frozen
floor **0.9**; per-head ECE **0.0119808** and **0.012886** against the frozen green band **0.05**;
and **neither head is a vacuous pass**. No threshold was touched, no seed was rerolled, the
evaluation set was not redrawn, and `P13_ITERATION_ALLOWANCE` remains **1 of 1, UNSPENT** — its one
pre-declared trigger was evaluated mechanically and did **not** fire.

This also answers, for these two criteria specifically, the open question plan 13-11 carried
forward: the persisted **epoch-4 best checkpoint of an early-overfitting run is sufficient** for
D-12 discrimination and D-13 calibration. It does not answer it for 13-13 or 13-16.

## Citing this result — the amended criterion never travels alone

`13-SC2-AMENDMENT.md` section 6 requires any report quoting an amended criterion to quote the
original beside it. The original ROADMAP Phase-13 success criteria, verbatim and still byte-present
at `.planning/ROADMAP.md:271-273`:

> 1. A 3-way `RatioEstimator`/evidence network emits a log-BF simplex over {coloc, random,
>    exclusion} in one forward pass
> 2. It reproduces `compute_BayesFactor()` (`src/bayes.jl:109`) in the overlapping 2-way regime
>    without quadgk/KDE
> 3. The exclusion hypothesis is validated on segregated ground-truth inputs

What this plan gated is the **amended** criterion: per-class one-vs-random AUC plus per-head ECE at
simulator ground truth. SC2's reproduction requirement was replaced — before any Phase-13 result
existed — because the KDE reference is invalid past `|logBF| ~ log(999) ~ 6.9` with 100%
baseline-side attrition, and because D-02 gives the two nets different input surfaces. SC1's
"`RatioEstimator`" and "log-BF simplex" wording was superseded by D-09/D-08 for architectural
reasons. **SC3 is not closed by this plan**: only its first arm (simulator ground truth) gates, and
the α-ladder and real-image arms are 13-13's and 13-16's, explicitly not gates.

## 1. Discrimination — the gate (D-12)

Per-class one-vs-random AUC at **simulator ground-truth labels**, each head scored only on its own
class pair (D-11), on the corrected log Bayes factors. Threshold-free by construction: an AUC is a
ranking statistic, so no decision rule and no operating point enters the gated number.

| Head | Statistic | Measured | Frozen floor | n (pos vs neg) | Verdict |
|------|-----------|----------|--------------|----------------|---------|
| coloc | AUC(log BF(C:R) \| coloc vs random) | **0.990169** | `P13_AUC_FLOOR_COLOC` = 0.9 | 1333 vs 1333 | **PASS** |
| exclusion | AUC(log BF(E:R) \| exclusion vs random) | **0.988262** | `P13_AUC_FLOOR_EXCLUSION` = 0.9 | 1334 vs 1333 | **PASS** |

Realized class counts E 1334 / R 1333 / C 1333 (frequencies 0.3335 / 0.33325 / 0.33325 against the
equal-thirds target), from 6,322 draws for 4,000 accepted items — acceptance overhead **1.58×**,
consistent with the class masses under π at τ = 0.15. Each head's restricted evaluation set holds
**2,666** (coloc) and **2,667** (exclusion) items, both far above the pre-registered
`P13_MIN_EVAL_PER_HEAD = 1000`.

The exclusion AUC resolved over |ρ_sample| bands. **These are reporting bands and carry no bar of
their own**; the frozen floor is applied unchanged inside each. They exist so the D-04 trigger,
which is worded about the deep-exclusion tail, can be evaluated mechanically rather than argued.

| \|ρ_sample\| band | n | exclusion AUC |
|---|---|---|
| (0.15, 0.50] | 637 | 0.979025 |
| (0.50, 0.75] | 325 | 0.995709 |
| (0.75, 0.90] | 133 | 0.998144 |
| (0.90, 1.00] | 239 | **0.997257** |

Discrimination **improves** toward the deep-exclusion tail; the weakest band is the one nearest the
τ boundary, which is what the design predicts and the opposite of the failure mode the iteration
trigger was reserved for.

## 2. Confusion matrix — DESCRIPTIVE, NOT A DECISION RULE

Reported at `argmax(log BF(C:R), 0, log BF(E:R))`, with **no assertion attached to any cell**. The
argmax can only declare `random` when *both* log Bayes factors are negative, which is a **decision**
rule; decisions and abstention are **Phase 14's** scope, and a descriptive argmax quoted as a
classifier would pre-empt that phase's design.

| true \ predicted | exclusion | random | coloc |
|---|---|---|---|
| **exclusion** | 1265 | 69 | 0 |
| **random** | 101 | 1148 | 84 |
| **coloc** | 0 | 66 | 1267 |

The two structural zeros are the informative cells: **not one** exclusion item was called coloc and
**not one** coloc item was called exclusion. Every error is an adjacent-class confusion with the
shared `random` reference.

## 3. Calibration — ECE gates, MCE does not (D-13)

Computed per head on that head's own restricted class pair, on the **UNCORRECTED** head
probabilities against the binary label, through the shared gate-lineage `_bin_calibration` at
`P13_ECE_NBINS = 10`. `_bin_calibration` was **not** modified.

| Head | ECE (GATED) | Band | Green cutoff | MCE (reported, never gated) | Empty bins | Head AUC | Vacuous? |
|------|-------------|------|--------------|-----------------------------|------------|----------|----------|
| coloc | **0.0119808** | `:green` | `P13_ECE_GREEN` = 0.05 | 0.0816278 | **0 of 10** | 0.990169 | **no** |
| exclusion | **0.012886** | `:green` | `P13_ECE_GREEN` = 0.05 | 0.116682 | **0 of 10** | 0.988262 | **no** |

Both passes are **substantive, not vacuous**: each head's AUC sits far above
`P13_VACUOUS_AUC_FLOOR = 0.60`, so neither is the "0.5, always" head that is perfectly calibrated
and perfectly useless. The AUC beside each ECE is computed on the very quantity the ECE was computed
on (the uncorrected probability) and coincides exactly with the gated AUC above, because
subtracting a constant is rank-preserving.

**The empty-bin trap did not bite this run** — 0 of 10 bins are empty on both heads
(coloc counts `[1111, 72, 40, 44, 48, 29, 49, 45, 100, 1128]`, exclusion
`[1063, 84, 54, 48, 51, 44, 48, 65, 110, 1100]`). MCE is still reported and still not gated: the
reason it is not the gate statistic is structural (an empty bin carries zero ECE weight and full MCE
weight), and a rule that happens not to bite on one run is not thereby a good gate.

**Context, not a threshold:** the shipped *binary* net's decision-calibration ECE in spike 014 was
**0.0188**. It is quoted beside the coloc head's 0.0119808 as a scale reference only — a different
net on a different input surface — and no Phase-13 verdict turns on it.

## 4. Reported observations — no assertion attached to anything below

### 4a. λ response (D-03)

Identical summaries, at otherwise matched inputs, re-encoded at each rung of Phase 11's own frozen
`SC2_RUNGS` ladder — only the appended conditioning value moves. The ladder was **read** from
`spike/validation/p11_consts.jl`, not chosen here.

| λ | mean \|logBF_C\| | mean \|logBF_E\| | mean logBF_C \| coloc | mean logBF_E \| exclusion | mean logBF_C \| random | mean logBF_E \| random |
|---|---|---|---|---|---|---|
| 0.25 | 7.36726 | 8.26806 | 5.93705 | 6.28207 | −6.19158 | −5.97236 |
| 0.50 | 7.36882 | 8.27439 | 5.95553 | 6.27987 | −6.18161 | −5.98289 |
| 1.00 | 7.37196 | 8.28694 | 5.99219 | 6.27519 | −6.16166 | −6.00391 |
| 1.50 | 7.37506 | 8.29935 | 6.02845 | 6.27014 | −6.14166 | −6.02489 |
| 2.00 | 7.37811 | 8.31161 | 6.06430 | 6.26472 | −6.12162 | −6.04582 |
| 2.50 | 7.38114 | 8.32371 | 6.09974 | 6.25892 | −6.10154 | −6.06671 |
| 3.00 | 7.38413 | 8.33566 | 6.13477 | 6.25274 | −6.08141 | −6.08755 |

Span of `mean|logBF|` across the whole ladder: **0.0169 nats** (coloc head) and **0.0676 nats**
(exclusion head). Spearman(rung, `mean|logBF|`) = **1.0** on both heads.

### 4b. Binary-NRE continuity (reported, `P13_CONTINUITY_GATED = false`)

Read at the reference λ = 3.0 over the **2,666** coloc/random evaluation pairs — the overlapping
two-way regime.

| Quantity | Value |
|---|---|
| corr(three-way log BF(C:R), binary log BF) | 0.62216 |
| max \|difference\| | 22.02 |
| mean difference | −4.2066 |
| n | 2666 |
| gated | **false** |

## 5. The D-04 iteration trigger, evaluated MECHANICALLY

`P13_ITERATION_TRIGGER` fires on exactly one condition: the exclusion-class AUC below
`P13_AUC_FLOOR_EXCLUSION` **specifically at high |ρ|** (the deep-exclusion tail near the −0.99 atom)
**while the mid-range is resolved**. Both halves must hold.

| Quantity | Measured | Floor | Half of the trigger satisfied? |
|---|---|---|---|
| exclusion AUC, deep tail \|ρ\| > 0.9 (n = 239) | **0.997257** | 0.9 | **no** — it is *above* the floor |
| exclusion AUC, mid range \|ρ\| ≤ 0.9 (n = 1095) | 0.986299 | 0.9 | yes (resolved) |

**THE TRIGGER DID NOT FIRE** (`iteration_trigger_fired = false`, persisted in the artifact). The
number the decision rests on is the deep-tail exclusion AUC **0.997257**, which is not below its
floor — it is the *best* of the four |ρ| bands.

**`P13_ITERATION_ALLOWANCE` is therefore UNSPENT (1 of 1) and NO RETRAINING IS AUTHORISED**,
independently of the gate verdict. `P13_STRATIFICATION_FALLBACK = :importance_weighted` is not
switched on, `consts.jl` is byte-unchanged, and nothing in the runner was adjusted in response to
anything it printed.

## 6. The λ response, read honestly

The response is **directionally correct and numerically negligible, and both halves of that sentence
matter.** Evidence magnitude moves **monotonically** across the ladder — Spearman 1.0 on both heads,
with every one of the seven rungs in order — so the conditioning row is demonstrably *wired* and the
net is *reading* it. But the total movement across the entire trained range, from λ = 0.25 to
λ = 3.0, is **0.017 nats** on the coloc head and **0.068 nats** on the exclusion head, against
evidence magnitudes of 7.4 and 8.3 nats. That is a **0.2% and 0.8% change** in the evidence over a
twelve-fold change in stated registration uncertainty.

**This is a finding about the conditioning input, not a defect, and it was pre-authorised as one.**
Phase 11 closed on precisely this result one level up: registration error at ≤3 px is not inferable
from the fixed 8×8 patch-correlation summary, and the posterior width the data demands is flat in λ
(RMSE ratio 1.0003). A three-way evidence net trained on that same summary cannot manufacture a
sensitivity the summary does not carry. The honest reading is that **the 8×8 summary does not encode
enough registration information for the evidence to depend materially on λ** — not that the λ input
is broken, and not that the evidence is wrongly confident.

Two directional details are worth recording rather than averaging away. The coloc head grows
*slightly more* confident as λ widens (5.937 → 6.135 on coloc items) while the exclusion head grows
*slightly less* (6.282 → 6.253); and the random-class means move in mirror image. Neither movement
is large enough to support an interpretation, and neither should be quoted as one. What can be said
is that the evidence pair does **not** stay flat by being disconnected — it is connected, and the
connection is small.

## 7. The continuity numbers, and why they are not a criterion

`corr = 0.622`, `max|Δ| = 22.02`, `mean Δ = −4.21`. **These are not a criterion and no Phase-13
verdict turns on them** (`P13_CONTINUITY_GATED = false`, frozen before any result existed). The two
architectural reasons, both decidable from source before a single weight was initialized:

1. **The two nets do not share an input surface (D-02).** This net reads Phase 11's
   registration-aware frozen `zt` *plus* an appended λ conditioning row; the binary net reads a
   different frozen `zt` and has no conditioning input at all. A disagreement is therefore a
   statement about the basis change and the conditioning input, not about either net's evidence.
   (The summaries were round-tripped exactly between the two frozen bases — `StatsBase.reconstruct`
   then re-standardize, an affine map on the continuous rows with the mask rows bypassed by both —
   so the two nets genuinely saw the *same acquisitions*. That is what makes the comparison
   meaningful as an observation; it is not what would make it a criterion.)
2. **The two logits are different objects.** The binary net's logit is a shuffled-θ
   likelihood-to-evidence ratio, so its model-index difference *is* the log Bayes factor with no
   correction term. The D-09 heads are plain BCE classifiers on class labels, whose Bayes-optimal
   logit carries `log(q_A/q_B)` and therefore **needs** a correction the binary net does not. They
   are not numerically interchangeable even at identical inputs.

The mean difference of −4.21 nats is consistent with exactly that: a systematic offset between two
differently-constructed quantities, not a disagreement about which pairs are colocalized.

## What this gate does NOT establish

Stated plainly so the pass is not over-read:

- **Simulator ground truth only.** Every label here comes from the D-05 two-factor cut on drawn θ.
  No real microscopy was scored, and no real image in this project has a colocalization ground-truth
  label. Labelled real segregation validation remains deferred to Phase 16's sealed holdout, which
  this phase does not open.
- **The well-specified regime.** The evaluation joint equals the training joint by construction
  (F5). This is a calibration-and-discrimination result under the simulator, not a misspecification
  result.
- **Not a decision rule.** The confusion matrix is descriptive. Turning this evidence into calls,
  with abstention, is Phase 14's job.
- **The exclusion hypothesis still has no observed real-data instance** anywhere in this project
  (`spike/simulator/ghat.jl:40`): every unmodified real image measured here sits at positive induced
  μ, including the "negative" control at +0.2481. That is why 13-13's α-ladder matters.

## Verification performed

| Check | Result |
|-------|--------|
| `julia --project=spike -t auto spike/p13/run_three_way_gate.jl` | **exit 0** — `Reported three-way pre-registered gate (D-12/D-13): 6 pass / 6 total` |
| Task-1 `<verify>` (artifact has `auc_coloc`, `confusion`, `P13_AUC_FLOOR_COLOC`) | prints `gate artifact ok` |
| Task-2 `<verify>` (nine outcome keys) | all present |
| All fifteen documented keys present | **15/15** |
| Thresholds recorded IN the artifact equal the frozen consts | AUC floors 0.9/0.9, ECE green 0.05, nbins 10, M 4000, min/head 1000, vacuous floor 0.6 — all match |
| `spike/p13/three_way_gate_report.jld2` exists, no `.tmp` sibling | 425,996 bytes, confirmed |
| `spike/figures/p13_confusion.png` rendered | 157,477 bytes |
| `git diff --quiet HEAD -- spike/p13/consts.jl` | **exit 0 — byte-unchanged** |
| `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | **exit 0 — byte-unchanged** (and asserted again at run time by step 0) |
| `grep -v '^\s*#' … \| grep -c '0\.90'` | **0** (thresholds referenced, never inlined) |
| `grep -v '^\s*#' … \| grep -c 'kde(\|quadgk('` | **0** |
| `grep -c 'DESCRIPTIVE, NOT A DECISION RULE'` | **3** |
| `@test` lines in the runner | **6**, none mentioning `mce` or `continuity` |
| `consts_git_blob_sha` in the artifact | `5a4ea2223ce11bcae2f974c643ec1e523458e8ac` — identical to the one 13-11 recorded |
| `net_consts_sha` == live `p13_consts_sha()` | `100e97a3…` both sides (asserted at run time) |
| Post-commit deletion check | no deletions |

**Test files run.** Per-file runs are the available gate signal for Phase 13: the full spike suite
aborts at `spike/test/runtests.jl:169` on the Phase-4 `SPEEDUP_GATE` (68.35 vs the pre-registered
bar 100.0), which masks the entire Phase-13 include block. That is pre-existing, already logged as
`deferred-items.md` D-13-A, and is **not** this plan's to fix — it is a pre-registered threshold in
a Phase-4 file. This plan's own deliverable is a **runner, not a test**, and its verification is the
run itself plus the artifact checks above; no `test_p13_*` file was modified, so none needed
re-running.

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 3 — Blocking] The figure is rendered BEFORE the gate, not after (plan step 14)**
- **Found during:** Task 1, writing the runner.
- **Issue:** the plan orders the figure at step 14, *after* the `@testset`. A failing `@testset`
  throws, so nothing after it runs — an after-the-gate figure would therefore exist only for a
  **passing** run. That contradicts the plan's own binding rule that the report must survive a
  failing gate, and the analog it names (`run_bf.jl:118-120`) emits its figure *before* the gate for
  exactly this reason.
- **Fix:** the figure is rendered immediately after the artifact is persisted and before the
  headline and the assertions, with the ordering and its reason written into the source comment.
  It carries no assertion either way.
- **Files:** `spike/p13/run_three_way_gate.jl`. **Commit:** `e5f9569`.

**2. [Rule 3 — Blocking] The confusion figure is composed in the runner, not routed through a
`figures.jl` helper**
- **Found during:** Task 1.
- **Issue:** the plan asks for the figure "through the existing `spike/validation/figures.jl`
  helpers into `spike/figures/p13_confusion.png`". Two things make that impossible as written:
  `figures.jl` has **no** confusion-matrix helper, and every helper it does have routes its output
  through `_val_fig_path`, which hard-joins to `spike/validation/figures/` — not the path the plan
  requires.
- **Fix:** `p13_gate_figure` composes the confusion heatmap and both ROC curves in one PNG using the
  CairoMakie surface, headless backend activation and caption discipline that `figures.jl`
  establishes and that the runner inherits by including it. The reason is documented in the
  function's docstring.
- **Files:** `spike/p13/run_three_way_gate.jl`. **Commit:** `e5f9569`.

**3. [Rule 3 — Blocking] Continuity reads the spike-lane binary NRE, not `src/amortized/bf.jl` +
`artifacts/grid_8/`**
- **Found during:** Task 1.
- **Issue:** the plan names `src/amortized/bf.jl:62-79`'s `amortized_log_bf` "used unchanged and
  read-only". That file `import`s **KernelDensity** and **QuadGK**, neither of which is a dependency
  of `spike/Project.toml`, so including it from the spike environment cannot work — which is the
  same reason `spike/validation/bf.jl` exists and says so in its own banner. Separately,
  `consts.jl` section I2 sanctions exactly **one** read of the shipped `artifacts/grid_8/` bundle,
  as a frozen **OOD comparison reference**, and this is not that read.
- **Fix:** the continuity section reads the frozen spike-lane binary NRE
  (`spike/validation/trained_ratio.jld2`, the artifact `run_bf.jl` itself reports against) through
  `spike/validation/bf.jl`'s `amortized_log_bf` — the same shuffled-θ `RatioEstimator` construction,
  carrying its own frozen `zt` and its own measured `log_prior_odds`. The choice, and which net was
  read, are printed in the run and recorded in the artifact's `continuity.note`. If the artifact is
  ever absent (it is gitignored and regenerable) the runner reports "NOT COMPUTED" with the reason
  rather than failing — continuity is reported-not-gated, so no verdict can depend on it.
- **Files:** `spike/p13/run_three_way_gate.jl`. **Commit:** `e5f9569`.

### Disclosures (not deviations, but the ordering must be checkable)

**4. A load smoke was run on the gate stream BEFORE the reported run, and the runner was not edited
afterwards.** To prove the runner executes end to end before spending ~11 minutes of CPU, it was run
once at `m = 90, min_per_head = 30`. That run consumed the reported gate stream and printed numbers
(AUC 0.99 / 0.984, ECE 0.0738 / 0.0789), so they were seen before the reported run. Three things
bound what that means:

- The runner refuses to write to the reported artifact path in that mode (it writes
  `*_smoke.jld2`), prints `SMOKE MODE -- NOT A REPORTED RUN`, and **does not run the gate
  assertions at all**. Both smoke outputs were deleted before the reported run.
- **Not one byte of the runner was changed after the smoke.** The file that produced the reported
  numbers is byte-identical to the file that produced the smoke numbers. No threshold, no statistic
  definition, no band edge and no seed moved in response to anything the smoke printed — including
  its ECE, which read *yellow*.
- Every threshold was already byte-locked in `consts.jl` and committed long before either run.

**5. The smoke's yellow ECE is itself evidence that `P13_MIN_EVAL_PER_HEAD` is load-bearing.** At
60 items per head the ECE read 0.0738 / 0.0789 (`:yellow`); at 2,666 / 2,667 it reads 0.0120 /
0.0129 (`:green`). Recorded with a caveat: the pre-registration warns that ECE is biased *downward*
at small n, and the direction observed here is the opposite — consistent with per-bin noise
dominating when ten bins share sixty samples. Both are reasons the minimum exists, and the small-n
number is not a measurement of anything.

**6. The figure is gitignored under the house `*.png` rule** (`.gitignore:401`), like every other
reproducible spike figure. It is regenerable from the committed runner and the committed artifact.
No untracked file is left behind. `spike/figures/plausibility.png` is negated back in because it is
a tracked *deliverable*; this one is a script artifact, which the plan itself says.

## Assumption Drift (advisory)

**A. "A flat λ response" was the pre-authorised expectation; what was measured is flat *and*
perfectly monotone.**
- **Planned:** the plan, 13-RESEARCH §D3 and the 13-09 ruling all pre-authorise a **flat** λ response
  as an honest finding, framed as the outcome to expect given Phase 11's negative result.
- **Actual:** the response is negligible in magnitude (0.017 / 0.068 nats over the full ladder) but
  **Spearman 1.0 on both heads** — all seven rungs in order, on both, in opposite directions.
- **Why:** "flat" and "disconnected" are not the same thing. The conditioning row is wired, appended
  once, and read; the summary simply does not carry enough registration information for it to move
  the evidence materially.
- **Materiality:** it changes the sentence the phase report should write. "The λ input has no
  effect" would be **false**; "the λ input has a monotone but negligible effect, because the 8×8
  summary does not carry the information" is what was measured. It also means a future summary
  redesign has something to move, rather than a dead input.

**B. The gate was framed as the place where an epoch-4 checkpoint might fail; it passed comfortably.**
- **Planned:** 13-11 recorded early overfitting as an honest finding and explicitly deferred "is a
  best-at-epoch-4 checkpoint sufficient?" to the downstream gates.
- **Actual:** both AUCs land at ~0.99 and both ECEs at ~0.012 — well inside the floors, not marginal.
- **Why:** the persisted artifact is the *best-validation* checkpoint, so the overfitting consumed
  budget rather than weights; and the D-11 masked loss gives the trunk all three classes, which is a
  strong signal for a small MLP on a 321-dim input.
- **Materiality:** the "overfits early" finding stays on the record and stays honest, but it must
  **not** be carried forward as a caveat on the D-12/D-13 numbers. It remains open for 13-13 and
  13-16.

## Known Stubs

None. Every number reported here was computed from the persisted artifact; no placeholder, mock or
hardcoded value flows into any table above. The one "not computed" branch in the runner (continuity
when the binary net artifact is absent) did not fire — 2,666 pairs were scored.

## Requirements satisfied

D-12 (per-class AUC gate + descriptive confusion matrix + reported-not-gated continuity), D-13
(per-head ECE gated, MCE reported beside its empty-bin count, every ECE beside its AUC, vacuous-pass
guard asserted), D-14 (band not p-value; no TOST machinery), D-08 (evidence emitted as two log Bayes
factors against a structural-zero random reference), D-05 (labels from the two-factor cut, never
from Δρ), D-03 (λ response measured and reported honestly), D-04 (persist-before-gate; thresholds
written into the artifact; iteration trigger evaluated mechanically; allowance unspent), D-01
(reserved disjoint stream, run-time decoupling proof, CPU-only, no package installed).

## Commits

| Commit | Task | Description |
|--------|------|-------------|
| `e5f9569` | Task 1 | run the amended three-way gate once — it passes on all six locked criteria |
| _(this)_ | Task 2 | record the outcome and the mechanical D-04 trigger evaluation |

Task 2 produced **no new artifact and no code change**, as the plan requires: it read
`spike/p13/three_way_gate_report.jld2` and wrote this SUMMARY. The artifact is byte-identical to
what Task 1 committed.

A parallel agent is executing Phase 12 on this branch. All commits here used explicit paths only; no
race occurred on the index, and nothing was rebased or amended.

## Self-Check: PASSED

- `spike/p13/run_three_way_gate.jl` — FOUND (859 lines)
- `spike/p13/three_way_gate_report.jld2` — FOUND (425,996 bytes, no `.tmp` sibling)
- `spike/figures/p13_confusion.png` — FOUND (157,477 bytes)
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-12-SUMMARY.md` — FOUND
- commit `e5f9569` — FOUND in history
- `git diff --quiet HEAD -- spike/p13/consts.jl` — exit 0
- `git diff --quiet HEAD -- spike/p13/three_way_gate_report.jld2` — exit 0 (read, not rewritten)
- `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` — exit 0
