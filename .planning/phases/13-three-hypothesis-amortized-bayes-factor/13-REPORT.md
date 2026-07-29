# Phase 13 — Three-Hypothesis Amortized Bayes Factor: Results Report

**Date:** 2026-07-29
**Status:** COMPLETE — the amended simulator gate cleared all six pre-registered criteria on one
run (`6 pass / 6 total`, exit 0), and both real microscopy fixtures are simultaneously flagged
out-of-distribution by both amortized detectors at ~2.5× their own in-distribution thresholds.

**Lane:** research (`spike/`) — **not shipped** (D-01). The shipped `artifacts/amended_v2/grid_8/`
bundle, the shipped `bayes_factor` accessor (`src/amortized/bf.jl`, `src/amortized/api.jl`) and the
2026-07-24 Phase-7 GO decision are **untouched by this phase**. Nothing in this report changes a
published claim about the shipped 8×8 estimator. `src/` is byte-unchanged across all sixteen plans
(§13).

**Every number in this report was read from a persisted artifact** — `spike/p13/tau_probe_report.jld2`,
`spike/p13/three_way_net.jld2`, `spike/p13/three_way_gate_report.jld2`, `spike/p13/alpha_report.jld2`,
`spike/p13/realimage_report.jld2` — or from `spike/p13/consts.jl`, and re-read at report time rather
than retyped from a console log or copied from a SUMMARY. Where a SUMMARY and an artifact disagree,
the artifact wins and the discrepancy is recorded in §14.

---

## Scope of evidence

> **Scope of evidence.** The Phase-13 gate is simulator ground truth (`rho_true < -tau`, exact labels,
> well-powered, `P13_GATE_M = 4000`). Two further arms are reported and **neither is a gate**: the
> D-16 semi-synthetic alpha-graded random-to-exclusion series (`P13_ALPHA_GATED = false`), and a
> qualitative check on six committed real microscopy TIFFs (`P13_REAL_IS_GATED = false`,
> `P13_REAL_QUALITATIVE_ONLY = true`). The TIFFs carry no colocalization label — the folder names
> `positive`/`negative` are the original package's biological test conditions, and the "negative" pair
> in fact has mean patch correlation **+0.2481** (`alpha0_mbar` = +0.24805233, `realimage_report.jld2`).
> The real-data arm can show the three-way verdicts behave sensibly on real microscopy; **it cannot
> show they are correct.** It covers two specimens, and both of them are flagged out-of-distribution
> by both the Phase-13 and the shipped detector, so every real-image number in this report is printed
> beside its OOD verdict and its lambda. A real image has no known registration uncertainty, so the
> whole pre-registered lambda ladder is reported and the widest, most conservative rung
> (`lambda = 3.0`) is the headline; lambda was never estimated from the images. The sealed-holdout
> rows of the provenance manifest were **not** opened.

| Arm | Substrate | Standing |
|---|---|---|
| Simulator ground truth | `rho_true < -tau`, exact labels, well-powered | **GATES** |
| D-16 alpha-graded series | mask-based disjoint reassignment, simulated **and** real | Supporting evidence — **not a gate** |
| `test/test_images/` qualitative check | six committed 1028 × 1376 microscopy TIFFs | Required deliverable — **explicitly not a gate** |

*No pass/fail threshold is defined for any real-image quantity.*

> **AMENDED AFTER THIS REPORT WAS FIRST WRITTEN — `13-D15-AMENDMENT.md` (2026-07-29), executed by
> plan 13-17.** The real-image arm's channel pair was **mis-designated**: `c1` is the DAPI/Hoechst
> nuclear counterstain, so both originally pre-registered pairs measured counterstain-versus-protein
> overlap rather than colocalization. The operative pair is now **`(2, 3)` = green/red**, the anchors
> are the unmasked **0.4603 / 0.3815**, and the redundancy arm is **dropped** — with `c1` excluded
> there is no second pair. The arm was **re-run once** on the corrected substrate.
> **The conclusion is unchanged**; §8c answers that question in words and prints every amended number
> beside its superseded counterpart. **No gating threshold moved**, no seed changed, no net was
> retrained, `P13_ITERATION_ALLOWANCE` is still **1 of 1 UNSPENT**, and **the phase's only gating arm
> (13-12, simulator ground truth) does not read the channel pair at all**, so §6a–§6d are untouched
> by this amendment. A **second, separately justified** change rode in the same edit — the
> `consts.jl` include-guard sentinel, which changed **no value** — and it is disclosed on its own in
> §13 and §14.7, not folded into the first. The figures quoted in the tables of §8b and §8c below are
> the **superseded** run and are retained as the record of what the counterstain pair produced.

This block is stated here, before every result, because the arm that will attract the most reader
attention (real microscopy) is the weakest one. Sections 8b, 8c and named limit D refer back to it;
they do not restate it as if it were a local caveat.

---

## 1. What this phase claims

The original ROADMAP Phase-13 success criteria are **byte-present and citable** at
`.planning/ROADMAP.md`; they were supplemented, never rewritten, by
`.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md`. Each amended
criterion is printed here beside its original, because reporting only the amended criterion is the
exact failure the amendment document exists to prevent (13-SC2-AMENDMENT §6).

| # | Original ROADMAP wording (verbatim) | Amended criterion | Pointer | Verdict |
|---|---|---|---|---|
| SC1 | *"A 3-way `RatioEstimator`/evidence network emits a log-BF simplex over {coloc, random, exclusion} in one forward pass"* | One network, one forward pass, three-way evidence emitted as `log BF(coloc : random)` and `log BF(exclusion : random)`, with the random entry identically 0. `RatioEstimator` is ruled out by D-09; "simplex" is ruled out by D-08. | `13-SC2-AMENDMENT.md` §4 | **MET** — `ThreeWayEvidenceNet` emits a `2 × n` logit matrix in one pass (§2); the structural zero is enforced by an inner constructor with `===`. |
| SC2 | *"It reproduces `compute_BayesFactor()` (`src/bayes.jl:109`) in the overlapping 2-way regime without quadgk/KDE"* | Reproduction **replaced** by per-class one-vs-random AUC (floors 0.90) plus per-head ECE (green band 0.05) at simulator ground truth. Correspondence with the shipped binary NRE is **reported, not gated** (`P13_CONTINUITY_GATED = false`). | `13-SC2-AMENDMENT.md` §§1–3 | **MET on the amended criterion** — 6/6, §6. The original reproduction requirement was **not attempted** and is not claimed. |
| SC3 | *"The exclusion hypothesis is validated on segregated ground-truth inputs"* | Retained and **scoped into three arms of deliberately unequal standing** (the table above). Only the simulator arm gates. | `13-SC2-AMENDMENT.md` §5 | **MET for the gating arm only** (§6). The alpha series (§8) and the real-image check (§8c) are supporting/qualitative and close no correctness claim. |

**The amendment was frozen before any Phase-13 result existed, and the ordering is checkable.** Its
§0 records the repository sha at time of writing, `df16e76941cefd28147efb4d4c91278597836cbd`
(2026-07-27), at which `spike/p13/` contained exactly one file — the Tier-1 pre-registration
`consts.jl` — with no trained net, no labelled evaluation set, no gate report and no alpha or
real-image run of any kind. `P13_TAU` is deliberately **absent** from the Tier-1 block, which is a
further check that the document precedes the phase's measurements.

**This is not Phase 7's position and must not be filed alongside it.** `07-GATE-AMENDMENT.md` was
written after its gate's outputs were seen. Phase 13 declined to be scored against a reference its
own prior work had already retired and already published as `docs/amortized.md` named limit 3 — a
ruling inherited, not minted, and readable from `src/bayes.jl` with every artifact deleted.

---

## 2. The architecture

One shared Flux trunk, copied **verbatim** from the shipped binary ratio estimator, feeding a single
two-output linear head, so both evidences come out of one forward pass and both retain logit/BCE
semantics (a softmax head would have broken the D-07 correction).

| Item | Value | Source |
|---|---|---|
| Trunk | `Chain(Dense(d_in, 256, gelu), Dense(256, 256, gelu), Dense(256, 256, gelu), Dense(256, 64))` | D-10, verbatim from `src/amortized/train_ratio.jl` (`RATIO_SUMMARY_WIDTH = 256`, `RATIO_NUM_SUMMARIES = 64`) |
| Head | one `Dense(64, 2)` → a `2 × n` logit matrix | D-09 |
| Input width | `three_way_input_dim(8, 1) = 5·G² + n_cond = 321` — **derived, never written as a literal** (`P13_INPUT_WIDTH_RULE = :ratio_input_dim_plus_ncond`; source grep for the literal returns 0) | `net_meta.input_dim = 321`, `net_meta.n_cond = 1` |
| `n_cond` | **1**, derived as `length(encode_lambda(3.0))` from Phase 11's own encoder, closing `d_in = 129 = 128 + 1` on the Phase-11 side | 13-09 |
| Lambda placement | appended **once**, after the pair encoding (`P13_LAMBDA_PLACEMENT = :append_after_pair_encode`); the `iseven` tripwire rejects the wrong placement | 13-06 |
| Loss | masked two-term logit-BCE; each head trains only on its own class pair (D-11), the trunk on all three | 13-06 |
| Training recipe | 300 epochs max, batch 128, lr 2.5e-4, weight decay 1e-4, val_frac 0.15, patience 40 — **copied verbatim** from the shipped binary net (D-10) so any difference is attributable to the head alone | `net_meta` |

**Why `RatioEstimator` could not express this.** `RatioEstimator(net, 1; num_summaries)` takes
`num_parameters = 1` — a *binary* model index — and its training loss is hard-coded
`logitbinarycrossentropy`; **a passed loss is silently ignored in v0.2.1**. A three-way head is not
reachable by passing a loss. This is an architectural fork, not a tuning knob, and it is why D-09
accepted the cost that this net **inherits no gate lineage from the shipped binary net** (limit C).

**Route R1 was used; the pre-declared route-R2 fallback was not needed.** Plan 13-06 proved, with a
real two-epoch training pass, that NeuralEstimators 0.2.1's own `train` accepts the custom
`ThreeWayEvidenceNet <: NeuralEstimator` subtype and **honours** the passed masked loss (the
two-different-losses check is the proof, since two freshly-initialized nets would differ anyway).
That matters beyond convenience: it means the optimizer, the Float64 CosAnneal schedule, the patience
rule and the best-validation checkpoint are byte-for-byte the shipped binary net's, which is exactly
what D-10's attribution argument requires.

**The shipped artifact is an epoch-4 checkpoint of an early-overfitting run, and that belongs here
rather than only in a SUMMARY.** The single training run reached its best validation risk
**0.13625102 at epoch 4** of a 300-epoch budget, then rose monotonically to **0.38632473** while
training risk fell to **0.036551084** — a ~10× train/validation gap — with early stopping firing at
epoch 45. NeuralEstimators' `train` returns the **best-validation** state, so the overfitting consumed
budget rather than weights, but the fact stands on the record. **The recipe was not tuned**, because
re-tuning even one value would have converted D-10's attribution argument into a coincidence and
would have spent the iteration allowance that `P13_ITERATION_TRIGGER` reserves for something else
(§10). Whether an epoch-4 checkpoint is *sufficient* was left to the downstream reads; §6, §8a and
§8c each answer it for their own arm, and none of them found it to be the binding constraint.

**The training pool.** 48,000 labelled items, realized class frequencies **exactly** equal thirds
(0.333333 / 0.333333 / 0.333333 against the `1//3` target), by whole-class accept/reject to an
integer quota of 16,000 each; acceptance overhead 1.545× (48,000 accepted from 74,149 draws). The
standardizer is Phase 11's, **inherited and provably not re-fit**: verdict `inherited`, with
`max|per-row mean| = 0.057353463` and `max|per-row sd − 1| = 0.035997152` over 96,000 standardized
columns, against a negative control that confirms a genuine re-fit lands at machine zero (< 1e-8,
i.e. more than 1e3× smaller than the observed residual). Training took **2.37 minutes**; pool
generation took ~120 minutes, so the spend is ~98% data generation.

---

## 3. The hypothesis boundary, and the trap

Labels come from the **D-05 two-factor cut**: the `rho_sample` **level** (its sign against `tau`)
**times** the contrast against the control. `COLOC` iff `rho_s > tau` **and** `rho_s − rho_c > tau`;
`EXCLUSION` iff `rho_s < −tau` **and** `rho_s − rho_c < −tau`; `RANDOM` otherwise. Labels are
**never** taken from `Δρ`.

### 3a. The trap, in plain prose

`src/amortized/infer.jl:117-121` defines `Δρ = mean(ρ_sample_draws .- ρ_control_draws)` — a
sample-versus-control contrast, and the quantity the shipped binary net splits on. A sample at
`ρ = 0.8` under a control at `ρ = 0.9` has `Δρ < 0`. The binary net calls that "null", which is
harmless. **A three-way test cut on the same quantity would publish it as "mutually exclusive" while
the sample is strongly colocalized.** `less colocalized than the control` is not segregation.

A reader may well assume the naive extension — it is the obvious generalization of a shipped,
working, published rule, and nothing in the binary net's own code warns against it. **Naming that
trap is itself a contribution of this phase**, independent of any number in this report, and it is
the reason D-05 is a two-factor rule rather than a sign test on one number.

The counter-example is not prose only: `three_way_label(0.8, 0.9; tau = TAU_FIX) != EXCLUSION` and
`=== RANDOM` are executable regression assertions at `spike/test/test_p13_labels.jl:111-112`, with
the mirrored form `three_way_label(-0.8, -0.9) === RANDOM` at line 148. The trap cannot be
reintroduced without turning the suite red.

---

## 4. tau, measured

`tau` was **measured, not chosen**. The whole `A(δ)` curve was computed and persisted **before**
`tau_from_curve` was ever evaluated, so no verdict could steer the measurement. Design: **unpaired**
(`P13_TAU_DESIGN = :unpaired`, disjoint reference and contrast arm keys, so the statistic is a genuine
two-sample separation rather than a paired difference), `P13_TAU_R = 400` draws per arm at every grid
point, statistic `:mbar_mask_weighted`, scored in **both** directions with the max rule
(`P13_TAU_BOTH_DIRECTIONS = true`), against the frozen bar `P13_TAU_AUC = 0.90`.

| δ | A_neg | A_pos | **A (max)** | bootstrap median | degenerate draws | clears 0.90 |
|---|---|---|---|---|---|---|
| 0.020 | 0.529025 | 0.592806 | 0.592806 | 0.592806 | 0 | no |
| 0.030 | 0.566481 | 0.625869 | 0.625869 | 0.625184 | 0 | no |
| 0.050 | 0.639375 | 0.689925 | 0.689925 | 0.688956 | 0 | no |
| 0.075 | 0.723788 | 0.760544 | 0.760544 | 0.761053 | 0 | no |
| 0.100 | 0.800338 | 0.819356 | 0.819356 | 0.821203 | 0 | no |
| **0.150** | **0.916356** | 0.905794 | **0.916356** | 0.917831 | 0 | **YES** |
| 0.200 | 0.978556 | 0.953613 | 0.978556 | 0.978541 | 0 | yes |

**`P13_TAU = 0.15`** — grid index 6 of 7, `tau_status = :measured`. The probe did **not** return
`:below_resolution`; the abort criterion `P13_TAU_ABORT_EXTEND_GRID = false` was never approached,
because δ = 0.20 also clears and the curve is monotone in δ in both directions.

**Reference lambda.** `P13_TAU_REFERENCE_LAMBDA = 3.0`, the realization of the pre-registered rule
`:widest_rung` against Phase-11's `LAMBDA_MAX`. tau is really tau(λ), measured at the **widest**
registration uncertainty, i.e. the hardest rung; it must never be quoted as an unconditional
resolution, which would overstate it at narrower λ.

**The both-directions rule earned its place.** At δ = 0.15 the two directions differ (0.916 negative
vs 0.906 positive); both cleared, so the max rule did not change the verdict, but the positive
direction alone would have cleared by only 0.006 and a single-direction probe would have reported a
near-miss as a comfortable pass.

**Sub-knot check:** `ghat_knot_spacing_near_zero = 0.0825`, and 0.15 > 0.0825, so
`tau_sub_knot = false` — the dead zone is representable in the simulator's own ρ mapping.

**Provenance, and the two-commit ordering as evidence.** `consts.jl` sha256 at run time
`640ec90edb80ef62eacbbf64c7b40cb9ca87c59d28a3a7438ef9d382b0aca4cf`; runner sha256
`d05057be5d3fc2bdfd491b2a932a5ad5a19fe81f0fd784974f9901a2d11586bd`; measurement sha256
`06eb1bc8fce8e9bf479cf20f314f100f1d0e1bcdbf39cc1fe01ea8d43921c1e6`; simulator-directory sha
`78dc37f517ad1b8e1e71afa963691f7ddefe5f63` with the runtime guard confirming the post-Phase-11
surgery (θ arity 8 with `chromatic_eps`, shift-prior half-width 3.0). The git history carries the
audit trail: `581602a` persists the artifact and **contains no tau constant**, `bfba6ac` appends the
Tier-2 constant, `45f5a44` flips the downstream assertions. The artifact provably existed at a commit
containing no tau value.

**Disclosed rather than hidden:** `repo_dirty_at_run = true`. The probe ran with an unclean tree
(concurrent Phase-12 planning work in `.planning/`, disjoint from `spike/`). The recorded `repo_sha`
therefore does not fully describe the tree; the two things tau depends on are pinned independently
and more tightly, by the `consts.jl` sha256 and the simulator sha plus the runtime guard.

**The honest trade-off this number carries.** A 0.15 dead zone is wide relative to the coloc scale:
it means "random" includes every pair whose sample-versus-control separation is under 0.15, and both
real fixtures land inside it (§8c). That is a property of the fixed 8×8 patch-correlation summary,
measured rather than assumed, and the full curve is persisted so any reader can re-read tau at a
different bar.

---

## 5. The D-07 correction, derived and verified

**The derivation, in five lines.** A BCE head's Bayes-optimal logit over its own restricted class
pair is `s_h(Z) = log[p_q(Z|A)/p_q(Z|B)] + log(q_A/q_B)`, where `q_A, q_B` are that head's **training
subset** class frequencies. Therefore `log BF_q(A:B) = s_h(Z) − log(q_A/q_B)` — the logit minus the
head's own log-odds. **The shipped binary net does not need this term**: its `RatioEstimator` loss
separates the joint `(Z, m)` from the product of marginals with a shuffle at equal counts, so its
Bayes-optimal logit already *is* `log[p(Z|m)/p(Z)]` and the model-index difference is the log Bayes
factor with no correction. **Copying the binary formula would therefore be wrong**, and wrong
invisibly. Finally, the scalar is exact **only if** stratification changed class *frequencies* and
not the within-class *shape*.

**Measured corrections** (`three_way_net.jld2`, training subset only, D-11 restricted):
`head_log_odds = (coloc = +0.0029420438990930215, exclusion = +0.0023543271616353607)`, against
π-level references `(coloc = −0.7764632135784509, exclusion = −0.40897645975438723)` computed from
200,000 prior class draws (π class masses 0.31272 / 0.47073 / 0.21655). The training-subset pair is
near zero **because stratification removed the prior imbalance** — but it is measured, not assumed:
the contiguous split does not divide a stratified pool into exactly equal thirds. The π-level terms
are 40–75× larger than the shipped binary net's own `−0.0102`, which is the concrete reason copying
by analogy was forbidden.

**The verification, against a closed-form conjugate-Gaussian toy** (plan 13-07; bars
`P13_F5_CORR_MIN = 0.99`, `P13_F5_MAXABS_TOL = 0.25` over the central `P13_F5_CENTRAL_FRAC = 0.90`
of each head's own evaluation band, `head_log_odds = (coloc = −0.49833, exclusion = +0.91050)`):

| Arm | Head | `corr` (bar ≥ 0.99) | `maxabs` (bar ≤ 0.25) | `meandev` | `maxabs` after removing the best constant | Reading |
|---|---|---|---|---|---|---|
| CORRECTED | coloc | 0.99751 ✅ | 0.5498 ❌ | −0.0693 | 0.4805 | corr clears, max misses |
| CORRECTED | exclusion | 0.99880 ✅ | 0.3783 ❌ | +0.0668 | 0.4452 | corr clears, max misses |
| UNCORRECTED (negative control 1) | coloc | 0.99751 | 1.0482 (must miss) ✅ | **−0.5676** | — | control fired |
| UNCORRECTED (negative control 1) | exclusion | 0.99880 | 1.1726 (must miss) ✅ | **+0.9773** | — | control fired |
| RESHAPED (negative control 2) | exclusion | 0.99831 | 1.0289 (must miss) ✅ | −0.7873 | **0.6745** | control fired, irreparably |
| RESHAPED (negative control 2) | coloc *(its classes untouched)* | 0.99873 | 0.3607 | −0.0535 | 0.3071 | unaffected, as predicted |

**The negative controls are what make this a verification rather than an assertion, and both fired
with the specific predicted magnitude.** Control 1: the uncorrected signed mean deviation lands on
each head's *own* measured log-odds to within 0.07 nat, and the identity
`uncorrected − corrected == log(q_pos/q_neg)` holds to 1e-9 — which rules out wrong sign, wrong
denominator and double application, each of which would also miss the bars but would *not* land on
that number. Control 1 also produced a finding worth carrying: **correlation cannot detect a missing
correction at all**, because a missing correction is a constant shift and correlation is
shift-invariant (uncorrected `corr` 0.99751, bit-comparable to corrected). A verification gating on
`P13_F5_CORR_MIN` alone would have passed a completely uncorrected net. Control 2: under a
within-class reshape at *identical* class frequencies the exclusion head's deviation grows to 1.0289
and remains **0.6745 after removing the best constant in hindsight** — the theorem as a number, not
an anecdote: no scalar of the D-07 form could have repaired it.

**The `maxabs` bar was MISSED on both heads, and that is recorded, not repaired.** The miss is **not
under-training** — measured, not argued: the analytic Bayes-optimal logits scored under the same
masked loss on the same held-out data give `net = 0.408983` vs `Bayes-optimal = 0.407191`, an excess
of **0.001792 nat (0.4%)**. The mechanism is confident-tail flatness of logit-BCE: both worst
deviations sit at the outer edge of the retained band where the true log BF is large (coloc
`Z = 1.765`, true +4.564, predicted +4.014; exclusion `Z = −2.188`, true +6.953, predicted +6.574),
and where a head is already right with probability near one a half-nat logit error costs almost
nothing in loss. About 90% of each band already meets the 0.25 bar; the bar is a maximum and the top
decile carries it over. The miss is robust to the choices the executor was free to make — the
deliberately data-starved coloc head reads 0.5498 at n = 6,000 and 0.509 at n = 48,000, an 8× data
increase for no improvement — so more data is not the lever. `spike/test/test_p13_correction.jl`
carries the two failing assertions **as-is**, uncommented and unsoftened, and neither the toy, the
bars, the widths nor the epochs were adjusted after the miss was seen. Amending the F5 statistic from
a maximum to a high quantile would have been a **third** amendment in this project's history; it was
**not** done (§"What was DECLINED"). The consequence is carried to Phase 14 in §11.

---

## 6. The gate

**VERDICT: PASS — all six pre-registered criteria cleared, on ONE run.**
`julia --project=spike -t auto spike/p13/run_three_way_gate.jl` → exit 0,
`Reported three-way pre-registered gate (D-12/D-13): 6 pass / 6 total`. Fresh 4,000-pair evaluation
set drawn at `P13_GATE_COUNTER = 3`, disjoint from the training-pool counter; realized classes
E 1334 / R 1333 / C 1333 (frequencies 0.3335 / 0.33325 / 0.33325) from 6,322 draws, acceptance
overhead 1.58×. No threshold was touched, no seed rerolled, and the evaluation set was not redrawn.

### 6a. Discrimination against the frozen floors (D-12) — GATED

| Head | Statistic | Measured | Frozen floor | n (pos vs neg) | Verdict |
|---|---|---|---|---|---|
| coloc | AUC(`log BF(C:R)` \| coloc vs random) | **0.9901687725007021** | `P13_AUC_FLOOR_COLOC` = 0.90 | 1333 vs 1333 | **PASS** |
| exclusion | AUC(`log BF(E:R)` \| exclusion vs random) | **0.9882624329245729** | `P13_AUC_FLOOR_EXCLUSION` = 0.90 | 1334 vs 1333 | **PASS** |

Each head's restricted evaluation set holds **2,666** (coloc) and **2,667** (exclusion) items, both
far above the pre-registered `P13_MIN_EVAL_PER_HEAD = 1000`. The gated statistic is threshold-free by
construction: an AUC is a ranking statistic, so no decision rule and no operating point enters it.

Resolved over `|rho_sample|` bands. **These are reporting bands and carry no bar of their own**; the
frozen floor is applied unchanged inside each. They exist so the D-04 trigger, which is worded about
the deep-exclusion tail, can be evaluated mechanically rather than argued.

| `|rho_sample|` band | n | exclusion AUC |
|---|---|---|
| (0.15, 0.50] | 637 | 0.9790253685870447 |
| (0.50, 0.75] | 325 | 0.9957089272318079 |
| (0.75, 0.90] | 133 | **0.9981442729103328** |
| (0.90, 1.00] | 239 | 0.9972566363348159 |

The weakest band is the one nearest the tau boundary, which is what the design predicts and the
opposite of the failure mode the iteration trigger was reserved for. Discrimination rises
monotonically across the first three bands and the deepest band sits a hair below the third
(0.997257 vs 0.998144) — see §14, which corrects a SUMMARY sentence on this point.

### 6b. Confusion matrix — **DESCRIPTIVE, NOT A DECISION RULE**

Reported at `argmax(log BF(C:R), 0, log BF(E:R))` (`P13_CONFUSION_RULE = :argmax_descriptive`), with
**no assertion attached to any cell**. The reason it cannot be a decision rule: the argmax can only
declare `random` when *both* log Bayes factors are negative, which is itself a decision — and
decisions and abstention are **Phase 14's** scope. A descriptive argmax quoted as a classifier would
pre-empt that phase's design.

| true \ predicted | exclusion | random | coloc |
|---|---|---|---|
| **exclusion** | 1265 | 69 | 0 |
| **random** | 101 | 1148 | 84 |
| **coloc** | 0 | 66 | 1267 |

The two structural zeros are the informative cells: **not one** exclusion item was called coloc and
**not one** coloc item was called exclusion. Every error is an adjacent-class confusion with the
shared `random` reference.

### 6c. Calibration (D-13) — ECE GATED, MCE reported beside its empty-bin count

Computed per head on that head's own restricted class pair, on the **uncorrected** head probabilities
against the binary label, through the shared gate-lineage `_bin_calibration` at
`P13_ECE_NBINS = 10`. `spike/validation/sbc.jl` was **not** modified.

| Head | ECE (GATED) | Band | Green cutoff | MCE (reported, never gated) | Empty bins | Head AUC | Vacuous? |
|---|---|---|---|---|---|---|---|
| coloc | **0.011980787564747793** | `:green` | `P13_ECE_GREEN` = 0.05 | 0.08162779293277045 | **0 of 10** | 0.9901687725007021 | **no** |
| exclusion | **0.012885988131763223** | `:green` | `P13_ECE_GREEN` = 0.05 | 0.1166815446181731 | **0 of 10** | 0.9882624329245729 | **no** |

Neither pass is vacuous: each head's AUC sits far above `P13_VACUOUS_AUC_FLOOR = 0.60`, so neither is
the "0.5, always" head that is perfectly calibrated and perfectly useless. The AUC beside each ECE is
computed on the very quantity the ECE was computed on, and coincides exactly with the gated AUC
because subtracting a constant is rank-preserving.

**The empty-bin trap did not bite this run** — 0 of 10 bins empty on both heads. MCE is still
reported and still not gated, because the reason it is not the gate statistic is *structural* (an
empty bin carries zero ECE weight and full MCE weight, demonstrated by an executable testset in plan
13-08), and a rule that happens not to bite on one run is not thereby a good gate.
`P13_TOST_REQUIRED = false`: a traffic light is a band, not a p-value, so it cannot be over-powered
the way Phase 7's M = 2000 SBC point-null was (D-14).

### 6d. What this gate does NOT establish

- **Simulator ground truth only.** Every label comes from the D-05 cut on drawn θ. No real microscopy
  was scored here.
- **The well-specified regime.** The evaluation joint equals the training joint by construction. This
  is a discrimination-and-calibration result under the simulator, not a misspecification result.
- **Not a decision rule.** §6b is descriptive; turning evidence into calls, with abstention, is
  Phase 14's job.
- **A disclosed load smoke preceded the reported run.** The runner was executed once at
  `m = 90, min_per_head = 30` to prove it runs end to end before spending ~11 minutes of CPU; it
  wrote to a `*_smoke.jld2` path, printed `SMOKE MODE -- NOT A REPORTED RUN`, ran **no** gate
  assertion, and both smoke outputs were deleted. **Not one byte of the runner changed afterwards**,
  and every threshold was byte-locked in `consts.jl` long before either run. The smoke's ECE read
  *yellow* (0.0738 / 0.0789 at 60 items per head) and nothing moved in response — which is also
  direct evidence that `P13_MIN_EVAL_PER_HEAD` is load-bearing.

---

## 7. Reported, not gated

### 7a. The lambda response (D-03)

Identical summaries at otherwise matched inputs, re-encoded at each rung of Phase 11's own frozen
`SC2_RUNGS` ladder (read from `spike/validation/p11_consts.jl`, not chosen here); only the appended
conditioning value moves.

| λ | mean \|logBF_C\| | mean \|logBF_E\| | mean logBF_C \| coloc | mean logBF_E \| exclusion | mean logBF_C \| random | mean logBF_E \| random |
|---|---|---|---|---|---|---|
| 0.25 | 7.367255 | 8.268062 | 5.937045 | 6.282073 | −6.191578 | −5.972361 |
| 0.50 | 7.368818 | 8.274391 | 5.955528 | 6.279871 | −6.181615 | −5.982888 |
| 1.00 | 7.371958 | 8.286943 | 5.992191 | 6.275192 | −6.161657 | −6.003908 |
| 1.50 | 7.375064 | 8.299347 | 6.028449 | 6.270142 | −6.141658 | −6.024885 |
| 2.00 | 7.378114 | 8.311609 | 6.064300 | 6.264719 | −6.121618 | −6.045818 |
| 2.50 | 7.381137 | 8.323709 | 6.099742 | 6.258921 | −6.101536 | −6.066705 |
| 3.00 | 7.384127 | 8.335663 | 6.134773 | 6.252744 | −6.081411 | −6.087547 |

`span_coloc = 0.016871622583837897` nats, `span_exclusion = 0.06760069735650376` nats,
`spearman_coloc = spearman_exclusion = 1.0`.

**Reported honestly: the response is directionally correct and numerically negligible, and both
halves of that sentence matter.** Evidence magnitude moves **monotonically** across the ladder —
Spearman 1.0 on both heads, all seven rungs in order — so the conditioning row is demonstrably wired
and read. But the total movement across a twelve-fold change in stated registration uncertainty is
**0.017 nats (0.2%)** on the coloc head and **0.068 nats (0.8%)** on the exclusion head, against
evidence magnitudes of 7.4 and 8.3 nats. **"The λ input has no effect" would be FALSE**; what was
measured is a monotone but negligible effect. A flat response was **pre-authorised** as an expected
outcome before the run, because Phase 11 closed on precisely this result one level up: registration
error at ≤ 3 px is not inferable from the fixed 8×8 patch-correlation summary, and the posterior
width the data demands is flat in λ (RMSE ratio 1.0003). A three-way evidence net trained on that same
summary cannot manufacture a sensitivity the summary does not carry. The honest reading is that the
**8×8 summary does not encode enough registration information for the evidence to depend materially
on λ** — not that the input is broken, and not that the evidence is wrongly confident. Two directional
details, recorded rather than averaged away and neither large enough to support an interpretation: the
coloc head grows slightly *more* confident as λ widens while the exclusion head grows slightly *less*,
and the random-class means move in mirror image.

### 7b. Binary-NRE continuity (`P13_CONTINUITY_GATED = false`)

Read at the reference λ = 3.0 over the 2,666 coloc/random evaluation pairs — the overlapping two-way
regime.

| Quantity | Value |
|---|---|
| corr(three-way `log BF(C:R)`, binary log BF) | 0.6221587375775006 |
| max \|difference\| | 22.02011533415178 |
| mean difference | −4.206613427934039 |
| n | 2666 |
| gated | **false** |

**These are not a criterion and no Phase-13 verdict turns on them**, for two architectural reasons
that were both decidable from source before a single weight was initialized. (1) **The two nets do
not share an input surface (D-02):** this net reads Phase 11's registration-aware frozen `zt` plus an
appended λ row; the binary net reads a different frozen `zt` with no conditioning input at all, so a
disagreement is a statement about the basis change and the conditioning input, not about either net's
evidence. (The summaries were round-tripped exactly between the two frozen bases via
`StatsBase.reconstruct` then re-standardize — an affine map on the continuous rows, mask rows bypassed
by both — so the two nets genuinely saw the same acquisitions. That is what makes the comparison
meaningful as an *observation*; it is not what would make it a criterion.) (2) **The two logits are
different objects** — see §5: the binary logit is a shuffled-θ likelihood-to-evidence ratio needing no
correction, the D-09 heads are BCE classifiers that *do*. They are not numerically interchangeable
even at identical inputs, and the mean offset of −4.21 nats is consistent with exactly that. The
comparison read the spike-lane frozen binary NRE (`spike/validation/trained_ratio.jld2`), because
`src/amortized/bf.jl` imports KernelDensity and QuadGK, neither of which is in the lean spike
environment, and because `consts.jl` §I2 reserves the shipped bundle for the OOD comparison only.

---

## 8. The alpha-graded series (two substrates, reported separately and NEVER averaged)

`P13_ALPHA_SUBSTRATE = (:simulated, :real)` is two-valued and pre-registered. What differs between
the arms is **what `alpha = 0` MEANS**: a genuinely random pair on the simulated substrate, a
moderately colocalized pair on the real one. A merged curve would be a statement about neither.
Figures: `spike/figures/p13_alpha_ladder.png` (simulated), `spike/figures/p13_realimage.png` (real),
`spike/figures/p13_confusion.png` (gate). All three are gitignored under the house `*.png` rule and
regenerable byte-identically from the committed runners and artifacts.

### 8a — simulated substrate (`alpha_report.jld2`)

64 images drawn at `rho_true` pinned **exactly** 0.0, one control acquisition each held fixed at
`rho = 0` across the whole ladder while only the sample is graded, `lambda = 3.0`, `G = 8`, reserved
stream at `P13_ALPHA_COUNTER = 4`.

| alpha | mean log BF(E:R) | median | mean log BF(C:R) | median | m-bar | `ghat(m-bar)` = induced rho | vs `tau` = 0.15 | absent patches |
|---|---|---|---|---|---|---|---|---|
| 0.000 | **−4.073865** | −4.094041 | −3.238762 | −3.491486 | +0.067767 | −0.02645 | inside the dead zone | 0 / 4096 |
| 0.125 | −1.137295 | −1.095068 | −5.693727 | −5.934632 | −0.018156 | −0.14430 | **inside, by 0.006** | 0 / 4096 |
| 0.250 | **+1.505085** | +1.438862 | −7.309971 | −7.499742 | −0.104272 | −0.24989 | **CLEARS tau** | 0 / 4096 |
| 0.375 | +3.678397 | +3.575812 | −8.067960 | −8.080470 | −0.186036 | −0.36684 | clears | 0 / 4096 |
| 0.500 | +5.446412 | +5.384720 | −8.403060 | −8.387613 | −0.259751 | −0.45129 | clears | 0 / 4096 |
| 0.625 | +6.842829 | +6.802206 | −8.633947 | −8.618346 | −0.323327 | −0.53094 | clears | 0 / 4096 |
| 0.750 | +7.912789 | +7.851003 | −8.852493 | −8.819391 | −0.376358 | −0.60285 | clears | 0 / 4096 |
| 0.875 | +8.724826 | +8.695825 | −9.061059 | −9.019144 | −0.419661 | −0.66511 | clears | 0 / 4096 |
| 1.000 | +9.344724 | +9.318827 | −9.250164 | −9.188022 | −0.454642 | −0.70471 | clears | 0 / 4096 |

**`alpha_star = 0.25`** — the smallest rung of the frozen `P13_ALPHA_LADDER` at which the mean
`log BF(exclusion : random)` first exceeds 0. Reported as a **RUNG**, never interpolated: the ladder's
own resolution is the resolution of the statistic, and the true crossing lies somewhere in
`(0.125, 0.25]` with this design unable to say where inside it. Per-image, the exclusion-positive
fraction runs 0/64 → 17/64 → 55/64 → 63/64 → 64/64 and the spread *narrows* monotonically from 1.88
to 0.82 nats; the coloc head moves the other way and never disagrees (negative on 64/64 at every rung
above `alpha = 0`).

**The interpretation — the correlation-versus-localization reading.** The simulator's "exclusion"
mechanism *is* negative intensity correlation across patches; a biologist's "mutually exclusive" is
*disjoint spatial localization*. This ladder is built **spatially** (ch2's above-background mass is
moved out of the ch1 object mask and returned, intensity-conserving, to the complement) and is scored
by a net trained only on the anti-correlation notion, so it asks at what alpha the second notion
starts to track the first. **The answer is that the two notions are separated by exactly the
pre-registered hypothesis boundary and by nothing else**: `alpha = 0.125` induces `rho = −0.1443`,
inside the measured tau dead zone by 0.006, and `alpha = 0.25` induces `rho = −0.2499`, the first
rung to clear it. The one-rung lag behind the summary's own sign flip **is** the tau dead zone, not an
unexplained gap. So: **roughly a quarter of ch2's above-background mass must be moved out of the ch1
objects before disjoint spatial localization registers as exclusion under the intensity-correlation
notion** — and below that, the segregation is genuinely present in the pixels and genuinely invisible
to the *summary*, which is a more useful place to have located the limitation than in the net.

**THE CAVEAT, STATED PLAINLY: alpha is a CONSTRUCTION parameter, not a physical quantity.** No
microscope has an alpha, no specimen has an alpha, and nothing in this project converts an alpha into
a biological degree of segregation. `alpha_star` therefore **CALIBRATES** the correlation-versus-
localization gap — how much of *this particular constructed transform* is needed before the evidence
responds — it does **not MEASURE** the gap as a property of real biology. A reader who quotes "0.25"
as "25% segregation is detectable" has over-read it. Two further scope limits: the ladder spans
induced rho `[−0.705, −0.026]`, so it does **not** probe the deep-exclusion tail near the −0.99 prior
atom; and one transform under one Otsu mask rule and one redistribution rule was used, so a different
segmentation or redistribution rule would give a different `alpha_star` and nothing here bounds by how
much.

**The ladder is not flat or shaped for a mechanical reason, and the evidence is in a table rather than
an intention.** `verify_alpha_invariants` ran on the **reported** substrate — all 64 items, not a
fixture — **before any curve was read**, with an `INVARIANT VIOLATION` early-return path that did not
fire.

| Invariant | Images failing | Measured |
|---|---|---|
| `bitwise_alpha0` (exact `==`, never `isapprox`) | 0 / 64 | — |
| `no_new_zeros` (`count(iszero)` preserved) | 0 / 64 | source zeros 0 … 10,659 per image |
| `intensity_conserved` (rtol 1e-10) | 0 / 64 | worst residual **9.63e-16** |
| `mask_invariant` (ch1 untouched, recomputed mask equals ladder mask at every rung) | 0 / 64 | — |
| `mask_fraction_ok` inside (0.01, 0.40) | 0 / 64 | realized 0.2840 … 0.3427, mean 0.3070 |
| `max_value_ok` | 0 / 64 | **VACUOUSLY TRUE — see below** |
| `mbar_monotone` | 0 / 64 | — |
| `missing_nonincreasing` | 0 / 64 | — |

| Warning sign | Measured | Reading |
|---|---|---|
| absent patches at `alpha = 0` / `alpha = 1` | **0 / 4,096** both | — |
| did the absent-patch count rise across the ladder? | **NO** | pixels are MOVED, not deleted |
| did the summary's MASK ROWS change with alpha? | **NO — byte-identical at every rung, all 64 images** | the zero set is untouched |

Had masked pixels been zeroed rather than pushed to a strictly-positive background floor, the measured
correlation would have been computed over the untouched complement only — approximately unchanged,
i.e. random — and the ladder would have read flat for a reason with nothing to do with segregation.
That failure did not occur.

**DECLARED LIMIT — `max_value_ok` is vacuously true on this arm and must never be quoted as
evidence.** The simulated arm declares an **unbounded** dynamic range (`declared_max_value = Inf`)
rather than the frozen `P13_ALPHA_MAX_VALUE_BOUND = 1.0`, which is a Gray-TIFF property of the *real*
substrate; the simulator's softplus intensities measured **4.079 … 17.964** here before any alpha, so
the frozen bound would have thrown on the first image for a reason that is a property of the substrate
rather than of the transform. Nothing was clamped; the realized maximum is recorded per image. And
**at `alpha = 1` the mask region is not EMPTY, it is AT the background floor** — the endpoint is "ch2
reduced to background wherever ch1 has objects", not "ch2 absent". Do not write "fully disjoint"
without that qualifier.

### 8b — real substrate (`realimage_report.jld2`)

Both conditions, at the headline `lambda = 3.0`, on the frozen `P13_REAL_ALPHA_GRID` (= the simulated
arm's ladder, so the rungs are comparable rung for rung), with the control half held fixed at the
*other* unmodified real fixture. All eight invariants hold on both conditions (worst intensity
residual 2.86e-16 / 1.72e-16; mask fractions 0.134933 / 0.229189; realized maxima 0.773286 / 0.586862,
both below 1.0 and recorded rather than clamped). **Every row of this arm is OOD-flagged; see §8c.**

| alpha | positive: log BF(E:R) | log BF(C:R) | m-bar | negative: log BF(E:R) | log BF(C:R) | m-bar |
|---|---|---|---|---|---|---|
| 0.000 | −10.27271 | −0.46577 | +0.32916 | −9.05337 | −3.94316 | +0.24805 |
| 0.125 | −9.56408 | −1.14050 | +0.27312 | −7.55256 | −5.17191 | +0.18357 |
| 0.250 | −8.65155 | −1.88244 | +0.21231 | −5.55352 | −6.38048 | +0.11407 |
| 0.375 | −7.49962 | −2.61873 | +0.14958 | −3.19712 | −7.33669 | +0.04554 |
| 0.500 | −6.02779 | −3.32545 | +0.08690 | −0.79561 | −8.01521 | −0.01866 |
| 0.625 | −4.16637 | −4.06013 | +0.02465 | **+1.54869** | −8.53514 | −0.07918 |
| 0.750 | −1.90658 | −4.93801 | −0.03851 | +3.92881 | −9.02240 | −0.13912 |
| 0.875 | **+0.70625** | −5.94974 | −0.10423 | +6.38335 | −9.57032 | −0.20059 |
| 1.000 | +3.40055 | −6.83244 | −0.16787 | +8.24502 | −10.10820 | −0.25198 |

**`alpha_star_real = (positive = 0.875, negative = 0.625)`** — each a RUNG, never interpolated.

**On real substrate `alpha = 0` is a moderately COLOCALIZED pair, not a random one** (measured
`m-bar` **+0.32916** / **+0.24805**, i.e. `+0.3292` / `+0.2481` against the frozen anchors to
3.7e-5 / 4.8e-5). So this ladder runs **moderately colocalized → near-exclusion**, while §8a runs
**random → exclusion**. The difference between `alpha_star = 0.25` there and 0.875 / 0.625 here is a
difference of **starting point**, not of net quality, and the two arms are never averaged.

**Scope limit.** The real ladder spans `ghat(m-bar)` from **+0.33042 down to −0.44184**, i.e.
`rho_true` in roughly `[−0.44, +0.33]` — a densely covered interior region of `GHAT_RHO_KNOTS`. It
probes the sign transition and the near-exclusion regime and **does NOT probe the deep-exclusion
tail** near the ±0.99 clamp atoms, where the prior-atom asymmetry (≈4.9% of prior mass at the negative
clamp against ≈1.9% at the positive one, about 2.5 : 1 against the exclusion end) lands hardest. That
tail was tested, and cleared, only by the D-12 gate on simulator ground truth (exclusion AUC 0.997257
at `|rho| > 0.9`).

**Grid truncation, in one sentence:** `patch()` truncates each frame to the largest exact multiple of
the grid, so at G = 8 the last **4 of 1028 rows (0.39%)** and 0 of 1376 columns are not used; each
patch holds 22,016 px and zero patches go missing. No crop, resize or pad step exists anywhere in the
arm.

---

## 8c. The real-image qualitative check (D-15 AMENDED TWICE) — reported, never gated

**Read the `## Scope of evidence` block in this report's header before reading this section.** It is
not restated here, deliberately, so its framing cannot be mistaken for a local caveat of this section.

### DOES THE CONCLUSION CHANGE? — the corrected channel pair, answered in words, first

**The channel pair was wrong, it has been corrected under an authorised amendment, the arm was
re-run once — and the conclusion is unchanged.** The numbers moved, one verdict moved, and the
overall reading did not: this arm is still **qualitative, n = 2 specimens, no ground-truth label,
OOD-bound, weak evidence — coherence and never correctness.**

That the limitation **survived its own correction is the point of reporting it.** It tells a reader
the weakness is **structural** — a property of two unlabelled specimens sitting outside the
simulator's training distribution — rather than an artefact of having read the counterstain
channel. A limitation that was never tested is a weaker statement than one that was.

**This correction is not a rescue, and three measured facts say so rather than a promise.**

1. **The binding OOD limit got WORSE, not better.** 703.2995 against a 167.5446 threshold —
   **4.198×** over, where the superseded pair read 2.491×. The shipped detector agrees and also
   worsens: 5.297× against 2.421×. The real fixtures now sit *further* outside the training
   distribution, not closer to it.
2. **The amendment adopted the narrower separation.** Unmasked, the corrected fixtures separate by
   **0.0788**, against **0.0811** for the superseded pair — and against the **0.8576** the masked
   read of the *same corrected pair* would have given, which is permitted by every prior source and
   was refused anyway (§14.6).
3. **The phase now reports one arm where it reported two.** The redundancy arm is **dropped, not
   deferred**: with `c1` excluded and three channels present there is no second pair.

**What DID change, stated as a change and not restated from the old run.** At the headline
`lambda = 3.0`, reading the **positive** fixture as sample, `log BF(C:R)` moved from **−0.4658**
(superseded pair `(1, 2)`) to **+5.70311** (amended pair `(2, 3)`), so the descriptive argmax moved
from **RANDOM to COLOC** in that direction. Reading the **negative** fixture as sample it stays
**RANDOM** (−8.65723 amended, against −3.9432 superseded). The headline OOD figure moved from
417.2974 / 167.5446 / 2.491× to **703.2995 / 167.5446 / 4.198×**. That figure was already known to
be strongly pair-dependent — the dropped redundancy pair `(1, 3)` had landed at 607.5728 (3.63×), a
46% move — so it is reported here **as a change**, never as a restatement.

**And two things this report previously claimed do NOT survive the correction** (§14.8). The
α-ladder **no longer crosses zero** on either fixture, so `alpha_star_real` is `nothing / nothing`
where it was 0.875 / 0.625; and the D-05 contrast of **0.09819** sits inside `tau = 0.15`, so the
phase's own label rule says RANDOM while the net says COLOC — the coherence check that previously
*passed* now *disagrees*. **Neither is repaired, retuned or argued away.**

### Every quantity that moved, each labelled with the pair it was measured on

| Quantity | SUPERSEDED — pair `(1, 2)` | **AMENDED — pair `(2, 3)`** |
|---|---|---|
| `m-bar` positive / negative | +0.32916 / +0.24805 | **+0.46027 / +0.38147** |
| headline `log BF(C:R)`, positive as sample | −0.4658 | **+5.70311** |
| headline `log BF(E:R)`, positive as sample | −10.2727 | **−12.9744** |
| headline `log BF(C:R)`, negative as sample | −3.9432 | **−8.65723** |
| headline `log BF(E:R)`, negative as sample | −9.0534 | **−12.8627** |
| descriptive argmax, positive as sample | RANDOM | **COLOC** |
| descriptive argmax, negative as sample | RANDOM | **RANDOM** |
| Phase-13 OOD density / threshold / ratio | 417.2974 / 167.5446 / 2.491× | **703.2995 / 167.5446 / 4.198×** |
| per-acquisition densities (positive / negative) | 417.2974 / 321.1505 | **703.2995 / 398.7601** |
| shipped OOD density / threshold / ratio | 433.6884 / 179.1368 / 2.421× | **948.9745 / 179.1368 / 5.297×** |
| A12 outcome | `agreement` | **`agreement`** |
| `alpha_star_real` positive / negative | 0.875 / 0.625 | **`nothing` / `nothing`** — never crosses zero |
| ladder `m-bar` span, positive | +0.32916 → **−0.16787** | **+0.46027 → +0.13806** (stays positive) |
| ladder `m-bar` span, negative | +0.24805 → **−0.25198** | **+0.38147 → +0.20414** (stays positive) |
| realized mask fraction positive / negative | 0.134933 / 0.229189 (**on c1**) | **0.051150 / 0.032858** (**on c2**) |
| eight α invariants | all hold, worst residual 2.86e-16 | **all hold, worst residual 1.80e-16** |
| ρ through the frozen `ghat`, positive / negative | +0.33042 / +0.21480 | **+0.48253 / +0.38434** |
| D-05 contrast vs `tau = 0.15` | 0.11562 — inside | **0.09819 — inside** |
| redundancy arm | reported at 607.5728 (3.63×) | **DROPPED — no second pair exists** |

**The mask fraction stayed inside its pre-registered band** (`P13_ALPHA_MASK_FRACTION_BOUNDS =
(0.01, 0.40)`) on the c2 mask, so the runner did **not** take its `INVARIANT VIOLATION` branch and
the ladder curves are readable. The band was **not** widened and was never at risk of being.

### The sweep, with every log-BF pair on the same row as its OOD verdict

The tables below are the **superseded** run, retained as the record of what the counterstain pair
produced. The amended readings are in the delta table above and in
`spike/p13/realimage_report.jld2`, which additionally persists `superseded_*` fields so the
artifact alone reconstructs the change.

Primary arm, **superseded** `P13_REAL_CHANNEL_PAIR = (1, 2)`. The read rule is
`P13_REAL_LAMBDA_READS_RULE = :full_phase11_ladder`, realized as Phase 11's
`SC2_RUNGS = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)`; the headline rule is
`P13_REAL_LAMBDA_HEADLINE_RULE = :widest_rung`, realized at **3.0** and asserted equal to both
Phase 11's `LAMBDA_MAX` and the frozen `P13_REAL_LAMBDA_HEADLINE_EXPECTED`. **A real image has no
known lambda, so the whole sweep is reported and the widest, most conservative rung is the headline**
— a conservative lambda makes the evidence weaker, so a verdict that survives it is the credible one.
Lambda was **never estimated from the images**; that would be a new estimator with its own validation
burden inside a phase that is not about registration. The motivating domain judgement — that widefield
chromatic aberration plus filter-cube/stage repeatability on these particular microscopes is ≥ 1 px —
is labelled **`[ASSUMED]`** in the artifact and gates nothing.

| direction | lambda | log BF(C:R) | log BF(E:R) | OOD density | OOD threshold | is_ood |
|---|---|---|---|---|---|---|
| positive as sample | 0.25 | −0.725239 | −10.261014 | 417.2974 | 167.5446 | **true** |
| positive as sample | 0.50 | −0.701851 | −10.262200 | 417.2974 | 167.5446 | **true** |
| positive as sample | 1.00 | −0.654947 | −10.264497 | 417.2974 | 167.5446 | **true** |
| positive as sample | 1.50 | −0.607885 | −10.266683 | 417.2974 | 167.5446 | **true** |
| positive as sample | 2.00 | −0.560664 | −10.268779 | 417.2974 | 167.5446 | **true** |
| positive as sample | 2.50 | −0.513290 | −10.270785 | 417.2974 | 167.5446 | **true** |
| **positive as sample** | **3.00 (HEADLINE)** | **−0.465771** | **−10.272712** | **417.2974** | **167.5446** | **true** |
| negative as sample | 0.25 | −4.118168 | −8.994916 | 417.2974 | 167.5446 | **true** |
| negative as sample | 0.50 | −4.102732 | −9.000272 | 417.2974 | 167.5446 | **true** |
| negative as sample | 1.00 | −4.071577 | −9.010961 | 417.2974 | 167.5446 | **true** |
| negative as sample | 1.50 | −4.040039 | −9.021613 | 417.2974 | 167.5446 | **true** |
| negative as sample | 2.00 | −4.008124 | −9.032229 | 417.2974 | 167.5446 | **true** |
| negative as sample | 2.50 | −3.975831 | −9.042814 | 417.2974 | 167.5446 | **true** |
| **negative as sample** | **3.00 (HEADLINE)** | **−3.943164** | **−9.053374** | **417.2974** | **167.5446** | **true** |

At the headline rung both log Bayes factors are negative in both read directions, so the descriptive
argmax lands on the RANDOM reference for both specimens. The density is constant across the sweep
because the lambda row is appended *after* the pair encoding and therefore never enters the summary
the density is computed on; it is printed on every row regardless, because the standing D-15
requirement is that no real-image quantity is ever quotable without its verdict beside it. The
per-acquisition scores behind the pair max are 417.2974 (positive) and 321.1505 (negative).

**Why RANDOM is the coherent answer here, by the phase's own label rule.** Mapping the measured
`m-bar` through the frozen `spike/simulator/ghat.jl` map (a read of frozen code, not a fit), the
positive fixture sits at `rho = +0.33042` and the negative at `+0.21480`, so the D-05 contrast is
**±0.11562 — inside the measured `tau = 0.15` dead zone**. Both specimens are individually correlated
(`rho_s > tau` on both) but to *similar degrees*, and `less colocalized than the control` is not
segregation — precisely the counter-example of §3a. **This is a consistency check that passed, not a
correctness check**: what was checked is that the net agrees with a rule derived from the same summary
it reads.

**SUPERSEDED BY THE AMENDED READING, AND THE COHERENCE ARGUMENT DOES NOT SURVIVE.** On the amended
pair `(2, 3)` the contrast is **0.09819** — still inside `tau = 0.15`, so the D-05 label rule still
assigns RANDOM — but the net's descriptive argmax is **COLOC** on the positive direction
(`log BF(C:R) = +5.70311`). **What was a consistency check that passed is, on the corrected pair, a
consistency check that disagrees.** It is recorded as a finding (§14.8) and is not repaired: the
net is not retrained, `tau` is not moved, and the label rule is not re-derived.

A redundancy arm through the second pre-registered channel pair `(1, 3)` re-read the **same two
specimens** and reached the same qualitative outcome (both log Bayes factors negative at every rung,
verdict RANDOM) at a *higher* density, 607.5728 — 3.63× over threshold. **THAT ARM NO LONGER
EXISTS.** `13-D15-AMENDMENT.md` CHANGE A **dropped** it — not deferred, not repointed — because
`(1, 3)` was DAPI-versus-red and, with `c1` excluded and the fixtures carrying exactly three
channels, **there is no second pair to re-read**. It measured −0.0925 masked on the *positive*
fixture, which is noise rather than corroboration, and a redundancy arm that corroborates nothing is
worse than none because it gets quoted as a second observation. **This phase now reports one arm
where it reported two.**

### The OOD finding — this arm's binding limit

**Both fixtures are flagged out-of-distribution by BOTH detectors, in both read directions, on both
pre-registered channel pairs.**

| detector | density | threshold | ratio | is_ood |
|---|---|---|---|---|
| **Phase-13** (this net, Phase-11 registration-aware basis, fit on its own 96,000-acquisition training pool) | 417.29737445899303 | 167.54463293780395 | **2.4907×** | true |
| **shipped** (frozen comparison reference, `artifacts/amended_v2/grid_8/`) | 433.68840347024536 | 179.13677813 | **2.4210×** | true |
| frozen pre-registration constants | `P13_REAL_OOD_SHIPPED_DENSITY` = **433.69** | `P13_REAL_OOD_SHIPPED_THRESHOLD` = **179.14** | 2.42× | true |

The shipped read **reproduced the frozen constants exactly** (`reproduced_frozen = true`), which is
what makes the comparison trustworthy rather than merely adjacent. The Phase-13 channel is the shipped
recipe applied to this net's own inputs — continuous-rows-only Mahalanobis, same 1e-6 ridge, operating
point at the copied `OOD_ID_QUANTILE = 0.95` — so no number was chosen. Realized in-distribution
quantiles: q50 = 40.17, q90 = 149.57, q95 = **167.54**, q99 = 199.47, max = 290.69. **The real
fixtures score above the maximum of 96,000 simulated acquisitions.** (Disclosed rather than buried:
that pool is class-frequency stratified, so its rho marginal is not the prior's; each detector is
referenced to its own net's training distribution, which is the comparison this arm is making.)

**This is not a bug and it is not a reason to skip the check.** It is the net's own misspecification
channel reporting, correctly, that real microscopy at this frame size and these statistics sits
outside the simulator's training distribution. Suppressing, softening or "fixing" it would be the
single most damaging thing this arm could do, which is why every number in every table above carries
its density, its threshold and its boolean verdict on the same row, in the printed output and in the
persisted artifact alike. **It is arguably the most honest single number the real-data arm produced.**

**The A12 outcome observed was `agreement`, and what it buys is narrow.** The registration-aware
Phase-11 basis did **NOT** move real microscopy back inside the training distribution: 2.49× over its
own threshold against the shipped net's 2.42× over its own. The two magnitudes are not directly
comparable (different standardizer, different training pool, different threshold), but the *verdict*
is the same, and that agreement **corroborates** named limit D rather than weakening it. Had the two
diverged, that divergence would itself have been the finding — a statement about what the
registration-aware basis buys on real data — and the named limit's OOD sentence would have needed
rewording. It did not.

### The naming correction (Open Question 8, answered YES)

**`positive`/`negative` are the original package's BIOLOGICAL test conditions, not colocalization
labels.** The "negative" pair measures mean patch correlation **+0.2481** (`rho_true` ≈ +0.215) — a
positively correlated pair, not an anti-correlated one. **A reader who assumes `negative` means "not
colocalized" will misread every real-image figure in this phase**: they would read the RANDOM verdict
on the "negative" specimen as a miss, when the specimen is in fact moderately colocalized and RANDOM
is what the pre-registered label rule assigns it. `P13_REAL_NAMING_CORRECTION` is printed verbatim in
the runner's banner before any code runs and is persisted verbatim into the artifact.

**Stated explicitly: this phase therefore holds NO real exclusion example at all.**

**A second naming trap — NAMED first, then CORRECTED under an authorised amendment.** Both originally
pre-registered channel pairs — `(1, 2)` and `(1, 3)` — include channel 1, which `test/runtests.jl:105`
records as the DAPI/Hoechst **nuclear counterstain**, not a target protein. Plans 13-15 and 13-16 both
detected this and both **declined to fix it**, naming the limit instead — **which was the right call
under the pre-registration as it then stood**, and is precisely why the correction required a dated,
argued amendment rather than a patch.

`13-D15-AMENDMENT.md` (2026-07-29) then made that correction: the operative pair is **`(2, 3)` =
green/red, the two target proteins**, the anchors are the unmasked **0.4603 / 0.3815**, and the
redundancy arm is dropped. **This is not an un-pre-registered third pair.** It is the pair
`spike/simulator/ghat.jl:59` — this project's own frozen calibration source — had **already** adopted
on 2026-07-25, 39 h before `consts.jl` was written; Phase 13 had drifted from its own upstream and the
amendment closes that drift. The full timeline, including the 41 minutes between the 13-16 run and the
user's disposition, is in §14.6 and is not softened.

The frozen c1/c2 anchor regression **survives, unchanged in meaning**, as a **read-chain identity
check** — it proves this ingestion is bit-for-bit the load path the frozen calibration used — and it
is **not** evidence that the c1/c2 pair measures colocalization. `spike/test/test_p13_real.jl` now
pins **both sides** of that check to `channels = (1, 2)` and to the literals 0.3292 / 0.2481, because
otherwise both sides would have moved together under the amendment, the assertion would still have
passed, and the identity proof would have been silently destroyed while its own comment went on
claiming to test it.

**The naming correction is SHARPENED, not retired, by the amendment.** On the operative pair the
"negative" biological control measures **+0.3815** — *more* positively correlated than the +0.2481
previously quoted, not less. **Both fixtures moved further into positive correlation.** This limit
belongs in the manuscript beside the folder-name correction.

### The qualitative-only statement

The six committed TIFFs carry **no colocalization ground-truth label** — no recorded rho, no recorded
Δρ, no segregation degree, no provenance metadata asserting one. So this arm can show that the
three-way verdicts behave *sensibly* on real microscopy (ingestion works, the evidence responds
monotonically to a constructed ladder, the verdict tracks the pre-registered hypothesis boundary,
nothing degenerates) and **it cannot show they are correct**, because correctness requires a label the
data does not have. **No pass/fail threshold is defined for any real-image quantity**, so this arm
cannot accidentally become a gate and no Phase-13 verdict turns on it.

**The n, stated plainly: two specimens.** Six TIFF files, because each specimen carries three
channels; the two pre-registered channel pairs give four reads, but they are four reads of the *same
two specimens*, so the effective independent n is **2**, not 4 and not 6. Two specimens cannot support
a rate, a coverage claim, a confidence interval or a comparison of conditions. **Labelled real
segregation validation remains deferred to Phase 16's blind evaluation**, which is what the sealed
holdout exists for and precisely why Phase 13 did not open it.

### The Phase-02 lineage note

`spike/simulator/ghat.jl:40` and `spike/NOTES.md:238-247` already record, as a frozen Phase-02
finding, that *"negative tail : physically reachable in real fluorescence = false ⇒ negative μ-prior
tail documented PRIOR-ONLY; consistency scoped to realized range"*. Real anti-correlation was never
observed, and every unmodified real image measured in this phase confirms it (`m-bar` +0.32916 and
+0.24805, both positive). **The exclusion hypothesis, which Phase 13 promotes to first class,
therefore has no observed real-data instance anywhere in this project.**

On the **superseded** pair the D-16 construction appeared to **sharpen** that verdict rather than
refute it: at `alpha = 1` it produced `m-bar` = −0.16787 / −0.25198 (`rho` ≈ −0.347 / −0.442) **from
real microscopy pixels**, which was reported as the first negative induced μ this project had obtained
from real data. The sentence was intended for the manuscript, verbatim:

> **Negative induced μ is CONSTRUCTIBLE from real microscopy pixels via the D-16 mask-based
> reassignment, but it has NOT been observed to occur naturally in the images this project holds.**

> **⚠ THAT SENTENCE IS RETRACTED. It was measured on the counterstain pair and DOES NOT HOLD on the
> operative one.** On the amended pair `(2, 3)` the α-ladder **never reaches negative `m-bar` at
> all**: +0.46027 → **+0.13806** (positive fixture) and +0.38147 → **+0.20414** (negative fixture),
> with `alpha_star_real` = `nothing` on both. Strict monotone decrease still holds; only the reach
> into negative correlation is gone.
>
> **Why, and it is not mysterious.** The α transform's mask is the Otsu mask of the pair's **first**
> channel. That was `c1`, the **nuclear counterstain**, and is now `c2`, the **green target**.
> Moving green intensity out of *nuclei* — two largely disjoint compartments — drives the
> correlation negative; moving red out of *green objects*, when green and red co-occur and share a
> bright background, does not. The measured c2 mask fractions are also much smaller (0.0511 /
> 0.0329 against 0.1349 / 0.2292 on c1), so far less mass is redistributed.
>
> **The corrected replacement claim, stated no more strongly than the data supports:**
>
> > **On the two committed specimens, the D-16 mask-based reassignment does NOT construct negative
> > induced μ from the green/red target pair, and negative induced μ has still never been observed
> > to occur naturally in the images this project holds. Whether it is constructible from real
> > pixels at all therefore remains OPEN on n = 2 unlabelled specimens — the earlier affirmative
> > was an artefact of segregating a target channel against a nuclear counterstain.**
>
> The `spike/test/test_p13_real.jl` assertion encoding the old claim is **left failing** rather than
> rewritten (§14.8). **This is a reduction in what the phase claims.**

### Read-only and seal evidence

`sha256(git ls-files -s test/test_images)` recorded by the runner **before and after** the run:
`eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276` both times — and matching the
digest plan 13-15 recorded independently, so the committed fixtures are byte-identical across both
plans. `git status --porcelain test/test_images/` is empty. **`open_sealed_holdout` was never called**
— not in a script, not in a test, not behind a flag — no Phase-13 executable references the sealed
provenance directory (the reported runner contains the string zero times, comments included, and
asserts that about its own source text at load), no fetch was performed, and the earlier
`checkpoint:human-verify` escape-hatch proposal was **withdrawn** by the amended D-15 rather than left
dormant. This arm consumed **no RNG stream at all**.

---

## What was DECLINED

Declines are first-class results of this phase and are listed here so they are not mistaken for
omissions.

1. **The KDE + `quadgk` Bayes-factor reference (`compute_BayesFactor()`, `src/bayes.jl:109`) was
   declined as a validation target**, and nothing in Phase 13 was scored against it. It reads
   posterior tail mass through a kernel density estimate over L = 999 draws, so the largest ratio it
   can *measure* rather than extrapolate is `|log BF| ≤ log(999) ≈ 6.907` — arithmetic on the
   reference's own construction. Beyond that it goes non-finite or saturates: on identical pairs,
   34.75% (139/400) non-finite and 32.5% (130/400) saturated, with dropped pairs at median `|Δρ|`
   0.736 against 0.283 for kept ones, and **zero** attrition on the amortized side. Attrition is 100%
   baseline-side and *not* random, so agreement computed on the survivors is survivorship-biased. The
   invalid band is exactly the confident-call tail where a three-way verdict asserts something. No
   `kde(` or `quadgk(` call appears in any Phase-13 executable (source grep: 0).
2. **Continuity with the shipped binary NRE was declined as a criterion** (§7b),
   `P13_CONTINUITY_GATED = false`, frozen before any result — for two architectural reasons, not
   because a disagreement was expected.
3. **MCE was declined as the gate statistic** (§6c) for a structural reason, and
   `spike/validation/sbc.jl` was **not** hand-modified to work around the empty-bin behaviour.
4. **The confusion matrix was declined as a decision rule** (§6b). Decisions and abstention are
   Phase 14's scope.
5. **No confusion matrix was computed on real data, and no real-image quantity was given a
   threshold.** Both refusals are structural: the runner contains zero `@test` lines.
6. **The single authorised iteration was declined** — not unavailable, *declined*, because its one
   pre-declared trigger did not fire (§10).
7. **A second amendment was declined.** The F5 `maxabs` shortfall (§5) could have been answered by
   amending the statistic from a maximum to a high quantile, which q90 already passes on both heads.
   That would have been a **third** amendment in this project's history, and
   `13-SC2-AMENDMENT.md` §7 states plainly that a second amendment is not authorised by that document
   and cannot be authorised by amending it. The shortfall is reported as measured instead.
8. **The `src/results.jl` rename was declined in-phase** (§11 item 1), and the sketch is not quoted as
   authoritative in the meantime.
9. **The corpus physical anchor and the CBS rows were declined**: the physical anchor is sealed for
   Phase 16, and CBS was dropped from scope entirely (its `0.0` end is *random*, not exclusion) rather
   than deferred, because the amended D-15 supplies a real substrate needing no fetch and no hash gate.
10. **A third channel pair was declined** after the fixtures had been measured (§8c), and that
    decline **still stands** — but it must be read precisely, because `13-D15-AMENDMENT.md` did
    subsequently change the pair. **The amendment does NOT add an un-pre-registered third pair.** It
    **corrects** the *first* pair to `(2, 3)` — the pair this project's own frozen calibration source
    `spike/simulator/ghat.jl:59` had **already** adopted on 2026-07-25, before `spike/p13/consts.jl`
    existed — and it **removes** the second pair rather than adding one, leaving the phase with
    **one** channel pair where it previously had two. What was declined during execution was
    *choosing a new pair on the executor's own initiative after seeing the numbers*; what happened
    instead was a dated, argued, separately-disclosed amendment that closes a drift from an upstream
    correction. §14.6 states the full timeline, including the 41 minutes between the 13-16 run and
    the user's disposition, without softening it.

---

## 9. Named limits

Written in the `docs/amortized.md` entry shape so a future productionization can lift them verbatim.
**`docs/amortized.md` was NOT edited** — a published named-limits entry is owed only if Phase 13
productionizes; while it stays research-lane, the amendment plus this report suffice.

**Limit A — the physical segregation anchor is n = 0 in practice.**
> The provenance corpus contains exactly **ONE** real segregated anchor — `neg-lightmycells-01`,
> `physical-primary`, condition `negative`, truth label `segregated`, channels
> `2:nucleus-dna,mitochondria` — and Phase 13 did not read it. It is **triply unavailable**, and each
> obstacle is independently sufficient: (1) its recorded `sha256` is the unfilled sentinel
> `"PENDING-FETCH"`, because the bootstrap fetch was deliberately never run, so there is no digest to
> verify any bytes against; (2) the corpus data directory is **empty** and the row records
> `bytes = 0`, which states exactly that nothing was fetched; and (3) the row's split is
> `sealed_holdout`, reachable only through the anti-snooping accessor. **Phase 13 deliberately did not
> open that seal.** The seal is not an inconvenience to be routed around — it is a control that exists
> for **Phase 16's** blind evaluation, and consuming the one segregated anchor here would irreversibly
> burn Phase 16 on the very hypothesis Phase 16 exists to evaluate. No seal-break escape hatch exists
> anywhere in Phase 13. The n = 1 physical segregation check is therefore **DEFERRED TO PHASE 16**
> (`P13_PHYSICAL_ANCHOR_DEFERRED_TO = "Phase 16"`, frozen in the pre-registration). No download-size
> claim attaches to this deferral: `bytes = 0` records only that nothing was fetched. **The exclusion
> claim therefore rests on simulator ground truth plus a semi-synthetic construction.**

**Limit B — "exclusion" is defined as negative intensity correlation, not disjoint localization.**
> The simulator's exclusion mechanism is a sign flip on channel 2's shared latent component, i.e.
> **negative intensity correlation across patches**; a biologist's "mutually exclusive" is **disjoint
> spatial localization**. These coincide often but not necessarily, and the D-16 alpha series is the
> only evidence in this project that probes the gap. On the simulated substrate it quantifies the gap
> as `alpha_star = 0.25`: roughly a quarter of channel 2's above-background mass must be moved out of
> the channel-1 object mask before disjoint localization registers as exclusion under the
> intensity-correlation notion — and the crossing lands exactly at the first rung whose induced rho
> clears the measured `tau = 0.15`, so the separation is the pre-registered hypothesis boundary and
> nothing else. **`alpha` is a CONSTRUCTION parameter, not a physical quantity**: no microscope and no
> specimen has an alpha, so `alpha_star` **calibrates** the gap for one particular constructed
> transform and does **not measure** it as a property of real biology. One mask rule (frozen Otsu), one
> redistribution rule; a different rule would give a different `alpha_star` and nothing here bounds by
> how much. The ladders span induced rho `[−0.705, −0.026]` (simulated) and roughly `[−0.44, +0.33]`
> (real), so neither probes the deep-exclusion tail near the ±0.99 clamp atoms.

**Limit C — this net is research-lane and inherits no gate lineage from the shipped binary net.**
> The three-way evidence network lives in `spike/`, trains on **Phase 11's** registration-aware frozen
> `zt` with an appended registration-uncertainty conditioning input (D-02/D-03), and sits **outside**
> `RatioEstimator` because that type takes a binary model index and hard-codes its loss (D-09). It
> therefore inherits **no** gate lineage from the shipped binary net and was validated entirely on its
> own terms, against thresholds frozen before it existed. Its logit is a different object from the
> shipped net's — a BCE class logit needing a per-head `log(q_A/q_B)` correction, against a shuffled-θ
> likelihood-to-evidence ratio needing none — so the two are **not numerically interchangeable even at
> identical inputs**, and the reported continuity numbers (corr 0.622, mean offset −4.21 nats) are a
> statement about that, not a disagreement about which pairs are colocalized. **Nothing here changes
> the shipped artifact, the shipped `bayes_factor` accessor, or the Phase-7 GO.** The shipped artifact
> was opened read-only, once, purely as an OOD comparison reference and never as an evidence basis.
> The persisted net is additionally the **epoch-4 best-validation checkpoint of an early-overfitting
> run** whose recipe was deliberately not tuned.

**Limit D — the exclusion hypothesis has no labelled real-data validation.**
> The three-way evidence network is gated on **simulator ground truth**. The real-image arm runs on
> six committed microscopy TIFFs (`test/test_images/`) that carry **no colocalization ground-truth
> label**; it is a qualitative behaviour check, not a correctness check, and it is **explicitly not a
> gate** — no pass/fail threshold is defined for any quantity in it. It covers **two specimens**, so it
> can support no rate and no coverage claim. Both specimens are additionally flagged
> **out-of-distribution** by both amortized detectors, in both read directions and on every channel
> pair ever read: on the **amended** pair `(2, 3)` the Phase-13 net scores density **703.30 against
> its own in-distribution threshold 167.54 (4.198×)** and the shipped 8×8 bundle **948.97 against
> 179.14 (5.297×)** — **worse than the superseded pair's 417.30 / 167.54 (2.49×) and 433.69 / 179.14
> (2.42×)**, and still agreement rather than divergence, so the registration-aware basis did **not**
> move real microscopy back inside the training distribution. The folder names `positive`/`negative`
> are **biological** conditions, not colocalization labels, and on the operative pair the "negative"
> pair measures `m-bar` **+0.3815** (**+0.2481** on the superseded pair) — *more* positively
> correlated, not less — so this phase holds **no real exclusion example at all**.
> **The channel pair was MIS-DESIGNATED in the pre-registration and has been CORRECTED, not merely
> named.** Both originally frozen pairs contained `c1`, the nuclear counterstain;
> `13-D15-AMENDMENT.md` (2026-07-29) amends the operative pair to `(2, 3)` and **drops the redundancy
> arm because no second pair exists**, and the arm was re-run once. **The corrected arm separates the
> two fixtures marginally *less* well (0.0788 against 0.0811), and its OOD limit is worse** — the
> correction bought no evidence. On the corrected pair the D-16 α-ladder **no longer reaches negative
> induced μ at all**, so the claim that negative μ is constructible from real pixels is **retracted**
> (§8c). The provenance manifest holds exactly **one** physically-segregated anchor; it
> is `split = sealed_holdout`, `sha256 = "PENDING-FETCH"`, `bytes = 0`, and reserved for the
> **Phase-16** blind evaluation, so it was deliberately **not** consumed here. **Labelled real
> segregation validation is deferred to Phase 16.**

**A fifth entry was considered and NOT added, and the reason is a measurement.** The pre-registration
anticipated a prior-atom limit — about 4.9% of prior mass at the `rho = −0.99` clamp against about
1.9% at the positive one, roughly 2.5 : 1 against the exclusion end — to be written as a named limit
**if the gate showed exclusion-end degradation**. The gate showed the opposite: exclusion AUC
**0.997257** in the deep tail `|rho| > 0.9` (n = 239), against 0.979025 in the band nearest tau. The
condition for the fifth entry was evaluated and did not hold, so it is recorded here as an evaluated
non-finding rather than silently dropped. The atom asymmetry remains a documented property of the
frozen prior and is the reason both alpha ladders' scope limits are stated in §8.

---

## 10. Iteration ledger (D-04)

**`P13_ITERATION_ALLOWANCE = 1`, authorised in `13-SC2-AMENDMENT.md` §7 before any result existed. It
was NOT spent. It remains 1 of 1, UNSPENT.**

**The one pre-declared trigger, verbatim from `three_way_gate_report.jld2`
(`iteration_trigger_text`):**

> THE ONE PRE-DECLARED CONDITION for spending P13_ITERATION_ALLOWANCE: the D-12 evaluation shows the
> exclusion-class per-class AUC below P13_AUC_FLOOR_EXCLUSION SPECIFICALLY AT HIGH |rho| (i.e. the
> deep-exclusion tail near the -0.99 atom), while the mid-range is resolved. The allowance is then
> spent on switching to P13_STRATIFICATION_FALLBACK (D-07-ii) and re-running ONCE. It is NEVER spent
> on relaxing a threshold after training. A second iteration is not authorised by this file.

**Evaluated mechanically, on the number the trigger is worded about:**

| Quantity | Measured | Floor | Half of the trigger satisfied? |
|---|---|---|---|
| exclusion AUC, deep tail `|rho| > 0.9` (n = 239) | **0.9972566363348159** | 0.90 | **no** — it is far *above* the floor |
| exclusion AUC, mid range `|rho| ≤ 0.9` (n = 1095) | 0.986299314554666 | 0.90 | yes (resolved) |

**`iteration_trigger_fired = false`**, persisted in the artifact. Both halves must hold and the first
does not. `P13_STRATIFICATION_FALLBACK = :importance_weighted` was **not** switched on, `consts.jl` is
byte-unchanged, and nothing in any runner was adjusted in response to anything it printed.

**Two clarifications that keep this ledger from being over- or under-read.**

- **The allowance is scoped to the GATE and could NOT have been spent in response to a real-image
  observation**, because no real-image quantity has a bar to miss (`13-SC2-AMENDMENT.md` §5(g)). An
  alpha-ladder curve or a real TIFF verdict cannot trigger a re-run under that document.
- **The discipline is about PRE-REGISTRATION INTEGRITY, not about compute budget.** Training took
  **2.37 minutes**; pool generation took ~120 minutes, so the marginal cost of a second training run
  on the existing pool is *minutes*. **It is cheap; that is not the reason it is forbidden.** Nobody
  should argue for a retrain on cost grounds.

---

## 11. Deferred to a later phase

1. **Rename `log_bf_simplex` → `log_bf_vs_random` and `bayes_factor_simplex` → `log_bf_vs_random` in
   the `src/results.jl` sketch (`src/results.jl:172` and `:178`).** *Owner: a productionization phase
   (user decision).* The sketch documents the field as `NTuple{3,Float64}` ordered
   `{null, coloc, anti-coloc}` — reference class FIRST — so a caller indexing `[1]` for "the
   colocalization evidence" would silently read the structural zero, and a field named "simplex" that
   is not a simplex will be quoted as one. The change is **comment-only, zero provenance cost**
   (`results.jl` is not in `DATAGEN_HASH_SRC_FILES`) and zero runtime impact. It was deliberately not
   done here because D-01's plainest reading forbids a `src/` edit during the spike and because
   whether D-01's enumerated scope (artifact / accessor / GO) was meant to be exhaustive is a user
   decision, not a planner inference. The executable rename already exists in the spike lane:
   `spike/p13/result.jl` defines `ThreeHypothesisColocResult <: AbstractColocResult` with a
   **name-keyed** evidence triple whose inner constructor enforces `random === 0.0`. **Until the rename
   lands, the `src/results.jl` sketch must not be quoted as authoritative.**
   `P13_RESULTS_RENAME_DEFERRED = true` is the frozen record.
2. **The n = 1 physical segregated anchor → Phase 16's blind evaluation.** *Owner: Phase 16.* See
   limit A. The `PENDING-FETCH` sentinels must be replaced with real digests before Phase 16 opens the
   holdout; that fetch is an explicit human decision and is not a Phase-13 one.
3. **Reshipping the three-way evidence net.** *Owner: a future productionization phase.* Materially
   cheaper here than in Phase 11/12 — Phase 13 touches neither the NPE, the summary, nor any
   `DATAGEN_HASH_SRC_FILES` entry — but deferred to keep the Phase-7 GO intact (D-01). Any reship must
   carry limit C's caveat that the persisted artifact is an epoch-4 checkpoint of an untuned,
   early-overfitting run, and should re-open the recipe question rather than inherit it.
4. **Extending the corpus with real graded segregation anchors.** *Owner: a roadmap change, not a
   Phase-13 decision.* This is the honest fix for limit A's n = 0, but Phase 8 is COMPLETE and
   reopening it is a roadmap-level call.
5. **The F5 confident-tail precision shortfall → Phase 14.** *Owner: Phase 14.* §5 measured that
   corrected log Bayes factors are unbiased in the bulk (mean deviation < 0.07 nat) and lose sub-nat
   precision in the confident tail, where the *sign* of the verdict is not in doubt but the
   *magnitude* is least determined. Phase 14's abstention layer consumes these magnitudes; thresholds
   placed in the confident tail would rest on the least-determined part of the estimate.
6. **The Phase-4 `SPEEDUP_GATE` red, and its drift.** *Owner: user decision, logged as
   `deferred-items.md` D-13-A.* The full spike suite exits 1 at the Phase-4 gate, which masks the
   entire Phase-13 include block, so Phase-13 test files were run directly. Two options already stand
   (confirm the flow-marginal hypothesis and record it as a named limit scoped to the research net, or
   explicitly retire the legacy gate as not meaningful against an 8-column prior), and a third question
   is now open: the measured median speedup has drifted **92.50 → 83.97 → 68.35** across the project's
   history (with a 50.402 reading recorded by plan 13-08 on the same machine), and whether that is
   machine-load noise or a real trend needs repeated measurement to answer. Re-deriving or relaxing a
   pre-registered threshold is not an executor call.

---

## 12. Open questions closed

All eight of `13-RESEARCH.md`'s open questions are closed.

1. **Does the comment-only `log_bf_simplex` rename land in `src/results.jl` during Phase 13?** —
   **NO, DEFERRED.** `src/` stays byte-unchanged (D-01); the executable rename lives in
   `spike/p13/result.jl`; the one-line productionization item is recorded in §11.1 with its reasoning.
2. **Which λ does the frozen scalar tau correspond to?** — **The widest rung.**
   `P13_TAU_REFERENCE_LAMBDA = 3.0`, the realization of the pre-registered `:widest_rung` rule against
   Phase 11's `LAMBDA_MAX`, stored as its own constant so the rule cannot silently re-point. tau must
   never be quoted without it.
3. **Does Phase 13 need the CBS `cbs-RG-000` bytes?** — **NO**, resolved in the research. CBS was
   dropped from scope (its `0.0` end is *random*, not exclusion), no fetch was triggered, and nothing
   under the provenance corpus directory was checked, fetched or referenced by any Phase-13 executable.
4. **Is the chromatic ε a second conditioning input?** — **NO.** `n_cond = 1`, derived at run time as
   `length(encode_lambda(3.0))`; ε is an inferred θ field (row 8 of the frozen box, bounded ±0.02), not
   a conditioning row. `d_in = 129 = 128 + 1` closes the arithmetic on the Phase-11 side and
   `three_way_input_dim(8, 1) = 321` on this one.
5. **argmax versus per-head threshold for the confusion matrix?** — **Gated threshold-free; argmax
   reported descriptively.** Per-class AUC is the gated statistic; the 3×3 matrix is a labelled
   descriptive companion (§6b); the decision rule is left to Phase 14.
6. **At which λ is the real-image check read, and which rung is the headline?** — **The full
   pre-registered ladder is reported and the widest rung is the headline.**
   `:full_phase11_ladder` realized as `(0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)`; `:widest_rung` realized
   at **3.0**, asserted at run time against both Phase-11's `LAMBDA_MAX` and the frozen
   `P13_REAL_LAMBDA_HEADLINE_EXPECTED`. Headline values: positive-as-sample `log BF(C:R) = −0.465771`,
   `log BF(E:R) = −10.272712`; negative-as-sample `−3.943164` / `−9.053374`. Both rules were locked in
   `spike/p13/consts.jl` before any real read, and **λ was deliberately NOT estimated from the
   images**.
7. **Should the Phase-13 net's OOD verdict be compared to the shipped net's?** — **YES**, and it was.
   `P13_REAL_OOD_COMPARISON = true`: Phase-13 **417.30 / 167.54 (2.49×)** against shipped **433.69 /
   179.14 (2.42×)**, the shipped read reproducing the frozen constants exactly. **A12 outcome:
   `agreement`.** The shipped bundle was opened read-only as a comparison reference only and was never
   an evidence basis; every log Bayes factor in this report came from the Phase-13 net on Phase 11's
   basis, with no fallback path.
8. **Does the `positive`/`negative` naming need an explicit correction?** — **YES.** Stated once, with
   the measured **+0.2481**, in §8c, and flagged here for the manuscript's data description alongside
   the counterstain-channel caveat.

---

## 13. Decoupling evidence

> **THIS SECTION EXISTS TO PROVE THAT NOTHING WAS TAMPERED WITH — AND THIS PHASE THEN AMENDED THE
> FILE IT WAS PROVING.** `13-D15-AMENDMENT.md` (2026-07-29), executed by plan 13-17, edited
> `spike/p13/consts.jl`. Three claims this section used to make therefore became **verifiably
> false**: that `consts.jl` had only two commits, that its blob sha was `5a4ea222…`, and that its
> sha256 `100e97a3…` "matches `consts_sha` and `net_consts_sha` in every artifact". Rather than
> quietly overwrite them, the section is **corrected here with the pre-amendment values kept beside
> the post-amendment ones**, and the deliberate sha divergence is disclosed as the designed
> property it is. A decoupling section that silently re-baselined itself after an amendment would
> be worth nothing.

Commands **re-run at 13-17 time** on the main working tree. Nothing below is copied from the
previous version of this table; every value was re-derived from live output.

| Command | Pre-amendment (as this section read before 13-17) | Observed now (post-amendment) |
|---|---|---|
| `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | exit 0 | **exit 0 — byte-unchanged** |
| `git diff --quiet HEAD -- test/ docs/ corpus/` | exit 0 | **exit 0 — byte-unchanged** |
| `git status --porcelain test/test_images/` | empty | **empty** |
| `git log --oneline -- spike/p13/consts.jl` | two commits: `c42cc8e` (13-01, Tier-1), `bfba6ac` (13-10, Tier-2 tau append) | **THREE commits.** A third now exists: **`9fe9b02`** (13-17) — *"amend consts.jl — TWO separately justified changes in one edit"*. **It postdates 13-16, and that is disclosed rather than hidden:** `13-D15-AMENDMENT.md` §0 opens on exactly this ordering. See §14.6 and §14.7. |
| `git rev-parse HEAD:spike/p13/consts.jl` | `70fe66df33832af2e648ed3f81cc3ad356ce4cd8` | **`d82caea64e19e8509f3abc1c09140cc7c2e2d497`** |
| `git hash-object --no-filters spike/p13/consts.jl` | `5a4ea2223ce11bcae2f974c643ec1e523458e8ac` | **`ef117e6ebc42095b777320ae0ceb680c13d7c6b9`** |
| `sha256(read("spike/p13/consts.jl"))` | `100e97a37bb470bb7cb4fbdbd8e33f219d83adcf61fe398a90cf3ad3391f3fa6` | **`d6a6e63bf19bb13ad5c61eeaaf9fc18fab27c5aa7ba879bad5c493a5c71ae247`** |
| CHANGE B — the include-guard sentinel | `spike/p13/consts.jl:100` guarded the whole Tier-1 body on `:P13_DEV_SEED` | **guards on `:P13_DECLARED_DEVIATIONS`**, at the body wrapper and all eight callers. `spike/validation/p12_consts.jl` is **byte-unchanged** and its `P13_DEV_SEED` mirror — which is **correct**, since asserting seed disjointness requires naming the seeds — is intact. `P13_DEV_SEED` itself is byte-unchanged. |
| `spike/p13/three_way_gate_report.jld2` | 425,996 bytes | **byte-unchanged — the gate was not re-run** |
| `spike/p13/alpha_report.jld2` | 53,391 bytes, 65 keys | **byte-unchanged — not re-run** |
| `spike/p13/realimage_report.jld2` | 81,649 bytes, 57 keys | **re-run by 13-17 on the amended pair; 65 keys** |
| `spike/p13/tau_probe_report.jld2` | 15,343 bytes | unchanged |
| `spike/p13/three_way_net.jld2` | 953,027 bytes | **unchanged — no net was retrained** |

### `consts_sha` and `net_consts_sha` now DIFFER, by design — the amendment's most load-bearing disclosure

The retired claim was that `100e97a3…` "matches `consts_sha` and `net_consts_sha` in every
artifact". After the amendment that is **false, deliberately**, and the four values are:

| Artifact | `net_consts_sha` (recorded at TRAIN time) | `consts_sha` (recomputed at RUN time) |
|---|---|---|
| `realimage_report.jld2` (re-run by 13-17) | `100e97a3…` — **pre-amendment** | `d6a6e63b…` — **post-amendment** |
| `gate_report.jld2` / `alpha_report.jld2` (not re-run) | `100e97a3…` | `100e97a3…` — both pre-amendment, because neither was re-run |

**Why this is legitimate:** training reads **no** `P13_REAL_*` constant (grep-verified,
`13-D15-AMENDMENT.md` §3), so the amended values are ones the trained net never consumed. No net
was retrained and no weight changed.

**The guard was re-derived to mean MORE, not less** (`13-D15-AMENDMENT.md` §7.3). The old
assertion `h.consts_sha == p13_consts_sha()` only required the two sides to **agree with each
other** — a retrain plus an undisclosed edit would have satisfied it silently. It is replaced, at
all three runner sites, by two assertions against a named, dated literal pair `P13_CONSTS_SHA`
(`spike/p13/net.jl`): the artifact side must equal `pre_amendment` and the file side must equal
`post_amendment`. A **third** sha256 on `consts.jl` — including a future undisclosed edit — now
fails, which the original form would have accepted. The widened
`… || h.consts_sha == PRE_AMENDMENT_SHA` disjunction was **rejected in a code comment**, because
once it exists the artifact side alone satisfies the first branch and any future drift passes
silently. **No guard was deleted or blanket-disabled.**

### The full-suite signature, BEFORE and AFTER CHANGE B

| | BEFORE (pre-13-17) | AFTER (post-13-17) |
|---|---|---|
| where `runtests.jl` aborts | **`runtests.jl:174`** → `test_p12_suite.jl:86` → `test_p12_decoupling.jl:150` | **`runtests.jl:223`** → `test_p13_real.jl` |
| `test_p13_consts.jl` | **never reached at all** | **137 pass / 137 total, 0 errors** |
| `test_p13_labels.jl` / `test_p13_alpha.jl` | never reached | **164/164** and **129/129** — executing and reporting |
| suite exit code | 1 | 1 |

**Two honest qualifications, because the raw comparison overstates what CHANGE B did.**

1. **The BEFORE run did not abort where the amendment predicted.** It aborted in the **Phase-12**
   block, before Phase 13 was reached, on `test_p12_decoupling.jl`'s *"the Phase-11 pool is
   read-only to Phase 12"* (25 pass / 1 fail). A **concurrently executing Phase-12 agent** fixed
   that between the two runs (`8880c29`), which is what let the AFTER run reach Phase 13 at all.
   **That movement is Phase 12's, not CHANGE B's**, and attributing it here would be false.
2. **The clean evidence for CHANGE B is a targeted reproduction, not the suite.** Loading
   `spike/validation/p12_consts.jl` (which defines the mirrored `P13_DEV_SEED`) and then
   `spike/test/test_p13_consts.jl`, in one process — exactly the suite's condition —
   produced **22 pass / 1 fail / 114 error, every error an `UndefVarError`** before the fix and
   **137 pass / 0 error** after it. `P13_DECLARED_DEVIATIONS` was confirmed **not** defined by the
   Phase-12 mirror, which is what makes it a sound sentinel.

**CHANGE B un-masks test files; it silences nothing.** The suite still exits 1, and it now exits 1
*later*, on real assertions rather than on 114 undefined constants.

**Run-time assertions, not just after-the-fact greps.** Every reported runner asserts the `src/` and
manifest decoupling at step 0, before it does any work, and asserts the `P13_CONSTS_SHA` pair on
both sides, so a net trained under a different pre-registration cannot be read **and** `consts.jl`
at any unauthorised sha fails.

**Resolve-risk clause (j)** in `spike/test/runtests.jl` asserts `ROCAnalysis`, `MLJ`,
`NormalizingFlows`, `InvertibleNetworks`, `Turing`, `ImageSegmentation`, `ImageMorphology`, `QuadGK`
and `KernelDensity` are absent as direct dependencies, and re-asserts NeuralEstimators v0.2.1 by UUID
`38f6df31-6b4a-4144-b2af-7ace2da57606`. All ten assertions ran and passed in a real suite run.
**No package was installed anywhere in this phase**, on any plan.

**Test status, stated exactly as observed rather than as hoped.**

- **This plan (13-14) wrote two Markdown files and ran no Julia test suite of its own.** The
  Phase-13 test files were not modified by it, so none needed re-running.
- The **root** package suite `julia --project=. -e 'using Pkg; Pkg.test()'` was re-run at report time
  and observed **green**: `Testing ProteinCoLoc tests passed`, exit 0 (its last testset,
  `amended gate machinery (07-GATE-AMENDMENT required code changes)`, reads 254/254). That is the
  executable proof that the read-only include of `src/results.jl` changed nothing in the shipped
  package. Plan 13-11 observed the same result independently.
- **UPDATED BY 13-17.** The **full spike suite** `julia --project=spike spike/test/runtests.jl` still
  **exits 1**, and that is the correct state: the Phase-4 `SPEEDUP_GATE` is `@test_skip`-**paused**
  (`a494247`, `DEFERRED-NPE-03-WALLCLOCK`, bar byte-unchanged) and `test_p13_correction.jl` carries two
  deliberate, committed, pre-registered MEASURED MISSES. What **changed** is that the Phase-13 block now
  **executes instead of being masked**: `test_p13_consts.jl` reports **137 pass / 0 error** where it
  previously reported 22 pass / 1 fail / **114 `UndefVarError`**, and `test_p13_labels.jl` (164/164) and
  `test_p13_alpha.jl` (129/129) run and report. That is the signal — **not the exit code**.
- **The suite now aborts at `runtests.jl:223`, on `test_p13_real.jl`, and 13-17 caused that.** Three
  assertions there are **MEASURED MISSES of the corrected channel pair**, left failing deliberately
  rather than rewritten (§8b, §14.8). Because a thrown testset aborts `include`, this masks the **seven**
  Phase-13 files after it. **All seven were therefore run individually and are recorded here:**
  `test_p13_tau.jl` 99/99 exit 0; `test_p13_net.jl` 93/93 exit 0; `test_p13_result.jl` 43/43 exit 0;
  `test_p13_calibration.jl` 95/95 exit 0; `test_p13_preconditions.jl` 74/74 exit 0;
  `test_p13_datagen.jl` 706/706 exit 0; `test_p13_correction.jl` **exit 1** with 88 pass and its **2
  deliberate** F5 `maxabs` failures. **Whether to accept the three new misses as named limits, to
  amend the D-16 real-substrate expectation, or to re-order the harness so the misses stop masking
  siblings, is a pre-registration decision that 13-17 escalated rather than resolved.**
- `test_p13_real.jl` reports **273 pass / 3 fail** (it read 34 pass / 0 broken before 13-17; the file
  gained the operative-pair regression and the pinned c1/c2 lineage check).

---

## 14. Discrepancies found while writing this report

Recorded because "the artifact wins" is only meaningful if the disagreements are printed.

1. **A SUMMARY sentence about the deep-tail band is wrong, and the artifact corrects it.**
   `13-12-SUMMARY.md` §5 states the deep-tail exclusion AUC 0.997257 "is the *best* of the four |ρ|
   bands". Read from `three_way_gate_report.jld2`, `band_auc = [0.9790253685870447, 0.9957089272318079,
   0.9981442729103328, 0.9972566363348159]`: the **highest** band is (0.75, 0.90] at **0.998144**, and
   the deep tail is the **second** highest. The dispatch briefing for this plan repeated the SUMMARY's
   wording. **The substantive claims are unaffected** — the deep tail is far above its 0.90 floor, the
   weakest band is the one nearest tau, and `iteration_trigger_fired = false` either way — but
   "the best of the four" is not what was measured, and §6a states the ranking as measured.
2. **`spike/p13/consts.jl` has two stable fingerprints and they are not in conflict.** The tracked git
   blob is `70fe66df…` (git applies its CRLF→LF filter before hashing); the value every runner recorded
   as `consts_git_blob_sha` is `5a4ea222…`, which is `git hash-object --no-filters`, i.e. the hash of
   the raw working-tree bytes. Both were verified at report time and both describe the **same file**,
   whose sha256 is `100e97a3…`. A reader comparing the two without knowing which is which would think
   the pre-registration had changed. It has not.
3. **The `runtests.jl` line number in `deferred-items.md` has drifted.** D-13-A records the abort at
   `spike/test/runtests.jl:169`; on the current tree `test_npe.jl` is included at **line 189**, because
   a concurrently-executing Phase-12 agent wired its aggregator in first (`261a0cc`, `afae172`). The
   abort itself — Phase-4 `SPEEDUP_GATE` at `test_npe.jl:230` — is unchanged, and so is the masking
   consequence for the Phase-13 block (now at lines 220–239).
4. **An advisory `consts.jl` COMMENT does not reproduce, and nothing is scored against it.** The
   `P13_REAL_ALPHA_GRID` comment block states the measured summary-level sign crossing sits at
   `alpha ≈ 0.49` (positive) and `≈ 0.46` (negative). Measured on the committed fixtures through the
   frozen ingestion, `m-bar` crosses zero at `alpha ≈ 0.674` on the positive fixture and `≈ 0.464` on
   the negative one — the negative figure reproduces, the positive one does not. It is a **comment, not
   a constant**; the executable pre-registration is unaffected, and the claim the comment was making
   ("strictly interior to this grid rather than pinned at an endpoint") remains true on both fixtures.
   `consts.jl` was **not** edited: the pre-registration is append-only and a comment is not a bar.
5. **Two Phase-13 SUMMARY files carry no `status:` frontmatter field** (`13-06-SUMMARY.md`,
   `13-07-SUMMARY.md`). Both plans completed; 13-07's own verdict line records PARTIAL (one bar
   cleared, one missed, both controls fired). Recorded here because the plan index counts a plan
   complete on a file's mere existence, and a missing `status:` has cost this project before.

---

**Discrepancies 6–8 were added by plan 13-17 (2026-07-29), which amended the pre-registration.
They are numbered separately, and CHANGE A and CHANGE B are given SEPARATE entries on purpose:
"we opened the file for X and also did Y" belongs on the page, not inferable from a diff.**

6. **CHANGE A — THE REAL-IMAGE CHANNEL PAIR NAMED THE WRONG PHYSICAL OBJECT, AND THE CORRECTION
   POSTDATES THE RUN IT AFFECTS.** `spike/p13/consts.jl` froze
   `P13_REAL_CHANNEL_PAIR = (1, 2)` and `P13_REAL_REDUNDANCY_PAIR = (1, 3)`, and
   `test/runtests.jl:105` records the fixture channels as `["blue", "green", "red"]` — so `c1` is
   the DAPI/Hoechst **nuclear counterstain**, not a target protein. Both pre-registered pairs
   contained it, so the arm as executed measured counterstain-versus-protein overlap, which is not
   colocalization. `13-D15-AMENDMENT.md` (2026-07-29) amends the pair to `(2, 3)`, **drops** the
   redundancy arm with no replacement, and amends the anchors to the unmasked `0.4603 / 0.3815`.

   **The timeline, stated without softening and without overstatement.** The factual correction was
   on the record **before Phase 13 executed**: `spike/simulator/ghat.jl` was corrected to the c2/c3
   pair at **2026-07-25 23:28** (`78dc37f`), and a STATE.md blocker naming Phase 13 specifically —
   *"MUST be corrected before Phase 13 executes"* — landed at **2026-07-25 23:56:54** (`ea4a7d3`),
   which also fixed `11-10-PLAN.md`. `consts.jl` was nevertheless locked at **2026-07-27 14:43:12**
   (`c42cc8e`) still quoting the superseded figures — **38 h 46 min 18 s after the blocker**.
   13-16 then ran at **2026-07-29 17:07** (`2112bed`) on the wrong pair, and the user's specific
   disposition (pair `(2,3)`, redundancy dropped, unmasked values, handle as an amendment) arrived
   at **17:48**, **41 minutes after that run finished** (`f411ee7`). **Neither half of that is
   obscured here.** What was settled before any Phase-13 result existed is *which channels are the
   target proteins*; what was ruled 41 minutes after the run is *how to dispose of an
   already-known defect*.

   **No gating threshold moved**, the phase's only gating arm (13-12, simulator ground truth) does
   not read `P13_REAL_CHANNEL_PAIR` at all, and `P13_ITERATION_ALLOWANCE` remains **1 of 1
   UNSPENT**. The amendment also **took the worse number**: unmasked separation 0.0788 against the
   superseded pair's 0.0811, and against the **0.8576** the masked read of the same corrected pair
   would have given — a figure `ghat.jl:69-71` itself calls the sharper separator and which no
   prior source forbids. It is refused on a code property (`patch_summary` applies no Otsu mask, so
   the net was trained on unmasked summaries), and the masked figures appear nowhere in this arm.

7. **CHANGE B — A LATENT LOADING DEFECT WAS FIXED IN THE SAME EDIT, AND IT IS A SEPARATE ACT WITH A
   SEPARATE JUSTIFICATION.** `spike/p13/consts.jl:100` guarded its **entire Tier-1 body** on
   `:P13_DEV_SEED`, and `spike/validation/p12_consts.jl:109` **legitimately mirrors** that name in
   order to assert seed disjointness. In a full-suite run the wrapper saw the mirror, skipped its
   own body, and ~104 Tier-1 constants were never defined — `test_p13_consts.jl` reporting 22 pass
   / 1 fail / **114 `UndefVarError`** and aborting eight sibling files. **Per-file runs passed,
   which is exactly why it stayed invisible.** The body wrapper and all eight callers now guard on
   `:P13_DECLARED_DEVIATIONS` — a name `consts.jl` alone declares, that is not a seed and not a
   prior bound. **`spike/validation/p12_consts.jl` was NOT touched: its mirror is correct, and the
   defect is on the guard side, every time.**

   **CHANGE B changes no value at all** — only which symbol an `if` tests — and it makes the suite
   report **more** failures, not fewer. **Why it rode in this amendment rather than a later one:**
   `consts.jl`'s sha256 changes for CHANGE A regardless, and the train-time `consts_sha` guard has
   to be re-derived in all three runners regardless, so CHANGE B cost **zero additional
   pre-registration integrity**; holding it back would have meant opening the byte-locked file a
   second time under a second amendment. This closes **DEF-12-03** (recorded in
   `.planning/phases/12-spatial-colocalization-map/deferred-items.md`), closed **by reference**
   here because a live Phase-12 executor owns that file. It is **1 of 29** poisoned guards a
   mechanical sweep found; **the other 28 are documented with a remedy each in
   `.planning/CONVENTIONS.md` C-01 and were deliberately NOT fixed.**

8. **THE CORRECTED PAIR RETRACTS TWO CLAIMS THIS REPORT PREVIOUSLY MADE, AND THE ASSERTIONS THAT
   ENCODE THEM ARE LEFT FAILING RATHER THAN REWRITTEN.** On the amended pair the D-16 α-ladder
   **never crosses zero** (`alpha_star_real` is `nothing` on both conditions, where it was
   0.875 / 0.625), so the sentence *"negative induced μ is CONSTRUCTIBLE from real microscopy
   pixels via the D-16 mask-based reassignment"* — offered in §8c for the manuscript — **does not
   hold on the operative pair**; it was measured on the counterstain pair. And the D-05 contrast is
   **0.09819**, inside `tau = 0.15`, so the phase's own label rule assigns RANDOM while the net's
   descriptive argmax says **COLOC** on the positive direction — so the *"consistency check that
   passed"* is, on the corrected pair, a consistency check that **disagrees**.
   `spike/test/test_p13_real.jl`'s three affected assertions were **NOT** rewritten to the
   newly-measured values: rewriting a substrate expectation to match what was measured *after* the
   measurement is the act the amendment exists not to be. They are reported, not adjusted, and the
   disposition is escalated.

---

*Phase: 13-three-hypothesis-amortized-bayes-factor*
*Report written: 2026-07-29 (plan 13-14). Sources: `spike/p13/{tau_probe_report,three_way_net,three_way_gate_report,alpha_report,realimage_report}.jld2`, `spike/p13/consts.jl`, `13-SC2-AMENDMENT.md`, and the fifteen preceding plan SUMMARYs.*
