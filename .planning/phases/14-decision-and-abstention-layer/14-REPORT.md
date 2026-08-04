# Phase 14 — Decision and Abstention Layer: Results Report

**Written:** 2026-08-04
**Repository HEAD at assembly:** `71e748782b6f19c9fda6a148c1f7a1219a45e5f0`
**Pre-registration freeze commit:** `a6c825867dbc786f7c3927df6760295ee0930c77`
(*"pre-register(14-01): freeze the Phase-14 Tier-1 pre-registration"*, Mon Aug 3 20:39:29 2026 +0200)

---

## Scope of evidence, and how every number below was obtained

**Every number in this report was read out of a `.jld2` artifact with `JLD2.load` and pasted.** No
value was retyped from a plan summary or from memory. Where a summary and an artifact could have
disagreed, the artifact is authoritative and the artifact key is named beside the number.

| Artifact | sha256 (on-disk at report time) | bytes | Sections it backs |
|---|---|---|---|
| `spike/p14/p14_conformal_report.jld2` | `0b1338a55fafc0a10bdac6c6bf694a613ec8eac7b0540b118ac22d93fe3ee7c9` | 31,539 | §5, §3 (SC1-d, SC1-e) |
| `spike/p14/p14_fdr_report.jld2` | `104cc9b430897f4b0d2d4f601214c9408acc5e28963b0528d73874d68e041f8c` | 48,609 | §4, §3 (SC1-a, SC1-b, SC1-c) |
| `spike/p14/p14_riskcoverage_report.jld2` | `4a7f293e514af89deb88dca72837787a2d81c9d86cb2a407b1927de6c3baab4e` | 416,102 | §6, §3 (SC3-a/b/c) |
| `spike/p14/p14_ood_arm_report.jld2` | `5d1549535b39186a557aea77222ba1c2e1f8091e98a8fd73aea60ca84320aaaf` | 738,407 | §6, §3 (SC3-d, SC2-b) |
| `spike/p14/p14_real_images_report.jld2` | `95da68c27e99a25ae35f807f9202419b1e16d8a22d4062a2d9821b86be4c3cd8` | 78,410 | §7, §3 (SC1-f) |
| `spike/p14/p14_eval_pool.jld2` | `4840aafca16a0618f3842f3bee0f91cc69a8201cf92bb866979919d20f806464` | 687,719 | the shared substrate for §4, §5, §6 |

All six are **gitignored, regenerable** outputs (`.gitignore: spike/p14/*.jld2`, added by plan 14-09).
The sha256 above is of the file as it stands on disk; the runners are the committed deliverables and
the counter-based Philox streams make every one of them reproducible bit-for-bit.

**The pre-registration these numbers were scored against.** `spike/p14/consts.jl` is
**byte-unchanged against HEAD** — `git diff --exit-code HEAD -- spike/p14/consts.jl` exits 0, verified
after every reported run and again while writing this report. Its content digest, recorded inside all
five artifacts, is `consts_sha = 4e8ee4be102ccbadd526d01302d1587279884e26bddfc9157336e03086d86ca5`.
Not one constant was edited, before or after seeing any result.

---

## 1. What this phase built, and what it is not

`decide_coloc` takes a batch of image pairs and a user-set Bayesian FDR level and emits, per pair,
one of **{`:coloc` / `:not_coloc` / `:abstain`}** — with every abstention carrying the named trigger
that produced it. Underneath it sits Phase 13's three-hypothesis evidence net, a three-class
posterior, a running-**mean** Bayesian-FDR prefix rule, a hand-rolled split-conformal hedge, and the
D-05 asymmetric fusion of the OOD channel, the conformal set and the classical cross-method
comparison.

**D-01, stated plainly and early.** The layer is built in `spike/`, against a **src-shaped**
signature: it takes a bundle plus images and returns an `AbstractColocResult` subtype, so a later
promotion into `src/` is mechanical rather than a rewrite. **"src-shaped" is not "in src". The
decision layer is NOT shipped in v2.0.** Nothing in this phase is reachable from the public API, and
no user of ProteinCoLoc v2.0 can call it.

That is enforced, not asserted. `spike/test/test_p14_decoupling.jl` (51 assertions, all passing)
checks, on every run and inside every runner at step 0:

| Guarded tree | Assertion | Status at report time |
|---|---|---|
| `src/` | `git diff --quiet HEAD -- src` **and** `git status --porcelain -- src` empty (the untracked case a diff cannot see) | clean |
| `spike/Project.toml`, `spike/Manifest.toml` | byte-unchanged; dependency name set is exactly the frozen sixteen; `ConformalPrediction` / `MLJ` / `ROCAnalysis` / `MultipleTesting` each asserted absent **by name** | clean |
| `corpus/` | byte-unchanged; sealed holdout never opened or named | clean |
| `artifacts/` | byte-unchanged | clean |
| `test/` | byte-unchanged; the six TIFFs' tree digest taken before **and** after the real-image run and asserted equal | clean, digest `eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276` both sides |
| CUDA | no CUDA module loaded in any Phase-14 process | clean |

`git status --porcelain -- src spike/Project.toml spike/Manifest.toml corpus test artifacts` is
**empty** at the time this report was written.

---

## 2. The two ROADMAP amendments, quoted beside the originals

Any citation of a Phase-14 result must carry **both** amendments **and** the original wording. Every
artifact carries the amendment string in a required key so the pairing travels with the number.

### SC1

> **AS ORIGINALLY WRITTEN (`.planning/ROADMAP.md`, Phase 14, criterion 1):**
> *"`decide_coloc(...)` emits calibrated calls at a user-set Bayesian FDR across a batch, using
> conformal sets (ConformalPrediction.jl) + decision-risk"*

**SC1 is AMENDED by D-02.** Conformal sets are met **in substance**, hand-rolled; the library name is
dropped. `conformal_library` is persisted as `"none -- hand-rolled; SC1 AMENDED by D-02"`
(`p14_conformal_report.jld2`).

**The reason, in full.** (i) `ConformalPrediction.jl` appears nowhere in this repository except in the
ROADMAP sentence itself — it is aspirational, not inherited. (ii) Adding it would modify
`spike/Project.toml` and `spike/Manifest.toml`, which a **running** test asserts byte-unchanged; the
frozen-environment guard would have to be amended or exempted to accommodate a library. (iii) The
library wraps **MLJ** models, not a NeuralEstimators posterior, so the integration is not a drop-in.
(iv) Split conformal is a sorted-array order statistic — `sort`, then index `ceil((n+1)(1-α))`. The
cost was concrete and the payoff was a name.

### SC2

> **AS ORIGINALLY WRITTEN (`.planning/ROADMAP.md`, Phase 14, criterion 2):**
> *"Abstention triggers on OOD ∨ cross-method disagreement ∨ ambiguous conformal set"*

**SC2 is AMENDED by D-05** to an asymmetric rule:

- OOD fires → **ABSTAIN**
- ambiguous conformal set → **ABSTAIN**
- cross-method disagreement **alone** → **DECIDE**, and record the disagreement as a named field
- cross-method disagreement **∧** OOD → **ABSTAIN**

**The reason.** Phase 9 exists so v2.0 can be positioned as *"knows when the classics are wrong, not
merely agrees with them."* The literal `∨` makes disagreement with Costes/Manders a trigger for
**silence** — so the tool would go quiet in exactly the cases that demonstrate its headline
differentiator. Separately, three OR'd triggers compound: a tool that is usually silent is not
useful. The asymmetry preserves SC2's protective intent without gutting Phase 9.

**Considered and not chosen, recorded here because it is the more transparent design if reviewers
push back on the asymmetry:** a four-state output adding **`DECIDED-CONTRA-CLASSICAL`** as an
explicit state, rather than folding the disagreement into a field on a `:decide`. It was rejected for
v2.0 as a non-standard output state needing extra defence, not because it is worse.

### SC3

SC3 (*"A monotone risk-coverage curve shows abstention concentrates on hard/OOD cases"*) is **not
amended**. One qualification of measurement, not of criterion: SC3-b is tested as a **rank
correlation on a declared coverage range**, never as a pointwise monotonicity assertion, because the
empirical curve is a step function whose denominator is the number of selected items —
`pointwise_monotonicity_asserted = false` is a persisted key.

---

## 3. The SC table — every row, every verdict

Rows are `14-VALIDATION.md`'s SC → validation map in its own order. **No row is without a verdict.**

| ID | Behaviour | Pre-registered number | MEASURED | Verdict | Artifact key | Reproduce |
|---|---|---|---|---|---|---|
| **SC1-a** | `p14_bayes_fdr` is the running-**mean** prefix rule | `v = [0.01, 0.02, 0.15]`, α = 0.10 ⇒ `k* = 3`, `fdp = 0.06`; `issorted(run)` on 10⁴ ascending vectors | fixture reproduces `k* = 3`, `fdp = 0.06`; substituting a running **max** turns the gate red in 11 places (watched, not argued) | **PASS** | `test_p14_fdr.jl` (55/55) | `julia --project=spike spike/test/test_p14_fdr.jl` |
| **SC1-b** | Realized FDP among decided-and-accepted items tracks α_FDR | realized FDP inside the binomial 95 % band of the rule's predicted `fdp` at every α ∈ {0.01, 0.05, 0.10, 0.20} | within band at **all four** levels; `sc1b_met = true`; decided fraction **0.888 (1776 of 2000)** at every level | **PASS** | `p14_fdr_report.jld2` → `alpha_table`, `sc1b_met` | `julia --project=spike -t auto spike/p14/run_p14_fdr_check.jl` |
| **SC1-c** | Prior sensitivity quantified, not hidden | **no bar — REPORTED, NOT GATED** | 20 rows persisted; `pi_sensitivity_gated = false`. At π_coloc = 0.70, α = 0.20 the rule predicts **0.1997520196632641** and realizes **0.3466666666666667** | **REPORTED, NOT GATED** | `p14_fdr_report.jld2` → `pi_sensitivity` | same runner |
| **SC1-d** | Split-conformal marginal coverage ≥ 1 − α_conf − 3·√(α(1−α)/n) | **≥ 0.8798753882025019** | **0.9115** (1823 of 2000); gap **+0.03162461179749809**; `q̂ = 0.4498001628978108` | **PASS** | `p14_conformal_report.jld2` → `coverage`, `band_lower`, `qhat` | `julia --project=spike -t auto spike/p14/run_p14_conformal.jl` |
| **SC1-e** | `p14_conformal_quantile` is the **order statistic**, not an interpolated quantile | n = 9, α = 0.10 ⇒ `q̂ == maximum(s)` and `q̂ != Statistics.quantile(s, 0.9)`; at n = 30 the α = 0.03 call **throws** | both equalities hold; the throw fires; substituting `Statistics.quantile` turns the gate red in 6 places (watched) | **PASS** | `test_p14_conformal.jl` (83/83) | `julia --project=spike spike/test/test_p14_conformal.jl` |
| **SC1-f** | Real-data check on the six committed TIFFs; corpus recorded unavailable-**by design** | six recorded decisions + triggers; `corpus_images_available == 0` as unfetched by design with the tier breakdown | 6 rows / 2 decision units, **both ABSTAIN** (`:ood_fired`); `corpus_images_available = 0`, `corpus_status = :unfetched_by_design`, tiers 30 `simulated-secondary` / 2 `physical-primary`, splits dev 14 / eval 16 / sealed_holdout 2 | **PASS** (structurally — SC1-f carries no threshold) | `p14_real_images_report.jld2` → `images`, `decision_units`, `corpus_status`, `tier_breakdown` | `julia --project=spike -t auto spike/p14/run_p14_real_images.jl` |
| **SC2-a** | `p14_fuse` implements the D-05 asymmetry exactly | all 3 × 3 × 2 × 2 = 36 cells; in particular `(:clear, :singleton, true, false) ⇒ (:decide, :decided_contra_classical)` | full truth table cell by cell; restoring SC2's literal `∨` in branch 5 turns the gate red in 6 places, two of them the cell that *is* D-05 | **PASS** | `test_p14_fuse.jl` (222/222) | `julia --project=spike spike/test/test_p14_fuse.jl` |
| **SC2-b** | D-06: missing / non-finite OOD threshold ⇒ `:not_checked` ⇒ abstain by default | absent / `Inf` / `NaN` all ⇒ `:not_checked` ⇒ `:abstain` with `reason = :ood_not_checked`; opt-out flips to `:decide` | all three cases hold; opt-out behaves; and in the reported runs `any_not_checked = false` — no item ever defaulted | **PASS** | `test_p14_fuse.jl`; `p14_ood_arm_report.jld2` → `any_not_checked` | `julia --project=spike spike/test/test_p14_fuse.jl` |
| **SC2-c** | No source uses `verdict.flag == false` as an in-distribution proxy | zero hits for `!.*\.flag` / `\.flag == false` outside `p14_ood_state` | zero hits, lane-wide, on comment-stripped Phase-14 sources | **PASS** | `test_p14_decoupling.jl` (51/51) | `julia --project=spike spike/test/test_p14_decoupling.jl` |
| **SC2-d** | Every classical-vs-amortized comparison passes through `ghat` | every `patch_correlation(` compared to τ is wrapped in `ghat(` | zero bare comparisons; the real-image rows record `basis = :ghat_rho_true_scale` | **PASS** | `test_p14_decoupling.jl`; `p14_real_images_report.jld2` → `images[*].cross_method.basis` | same |
| **SC3-a** | Confidence ordering beats random (selective skill) | `skill ≥ P14_SKILL_FLOOR = 0.60` | **0.9017113386798696** | **PASS** | `p14_riskcoverage_report.jld2` → `skill`, `sc3a_met` | `julia --project=spike -t auto spike/p14/run_p14_riskcoverage.jl` |
| **SC3-b** | Selective risk rises with coverage | `corspearman ≥ 0.95` on coverage ∈ [0.20, 1.00], floor frozen before the run | **0.9701942302111639** on 1601 of 2000 curve points | **PASS** | same → `spearman`, `spearman_n_points`, `sc3b_met` | same |
| **SC3-c** | Abstention concentrates on hard cases | `AUC_hard ≥ 0.80` via the in-repo `roc_auc` | **0.9117431863910738** (154 hard items of 2000) | **PASS** | same → `auc_hard`, `n_hard`, `sc3c_met` | same |
| **SC3-d** | Abstention concentrates on OOD cases | `min over families of (abstain_rate(family) − abstain_rate(ID)) ≥ P14_OOD_MARGIN_FLOOR = 0.50` | **−0.129**, weakest family `:noise` | **FAIL — NOT MET** | `p14_ood_arm_report.jld2` → `margin`, `margin_per_family`, `sc3d_met = false` | `julia --project=spike -t auto spike/p14/run_p14_ood_arm.jl` |
| **D-07** | τ loaded, four-way agreed, sha-pinned | `net.tau == probe.tau == gate.P13_TAU == P13_TAU == 0.15`; net sha == pre-amendment; `p13_consts_sha()` == post-amendment; **probe sha ∉ both** | `tau = 0.15`, `four_way_agreed = true`; `net_consts_sha = 100e97a3…`, `file_consts_sha = d6a6e63b…`, `probe_consts_sha = 640ec90e…` — the probe **diverges**, asserted rather than papered over | **PASS** | every artifact → `provenance` (`tau`, `four_way_agreed`, the seven shas) | `julia --project=spike spike/test/test_p14_provenance.jl` (77/77) |
| **D-01/env** | `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `corpus/`, `artifacts/` untouched; 16 deps exactly; no CUDA | all assertions restated for Phase 14 | all clean; dependency set exactly the frozen sixteen; no CUDA module loaded | **PASS** | `test_p14_decoupling.jl` (51/51) | `julia --project=spike spike/test/test_p14_decoupling.jl` |
| **Seeds** | P14 streams provably disjoint from every P13 and forbidden stream | P14 seeds ∉ `_p13_forbidden()`; counters pairwise distinct; **full Philox key-word cross-product** pairwise distinct vs the P13 family | 175 cross-product assertions pass; `master_seed = 0x000000000b14de71`, `salt = 0x5851f42d4c957f2d`; counters cal 1 / eval 2 / ood 3 / real 4 | **PASS** | `test_p14_consts.jl` (273/273) | `julia --project=spike spike/test/test_p14_consts.jl` |
| **Class order** | Every named field lands on its intended class | fixture with **unequal** posteriors 0.7 / 0.2 / 0.1; `p.coloc` / `p.random` / `p.exclusion` each equal their intended value | all three hold; a deliberately permuted class assignment fails the gate in 9 places (watched, not argued) | **PASS** | `test_p14_posterior.jl` (48/48) | `julia --project=spike spike/test/test_p14_posterior.jl` |
| **Prior re-derivation** | `pi_class_masses` is re-derivable, not taken on trust | re-derived masses match to the fixture's MC error | re-derived; the reported prior in use is `class_prior_used = (exclusion = 0.31272, random = 0.47073, coloc = 0.21655)` and is carried in every artifact | **PASS** | `test_p14_posterior.jl`; every artifact → `class_prior_used` | same |

**Tally: 17 PASS, 1 REPORTED-NOT-GATED, 1 FAIL.** The FAIL is SC3-d and it is stated at the same
volume as the passes, in §6 and in §8. Its disposition was **ruled by the user on 2026-08-04**
(`14-SC3D-RULING.md`, commit `f4934a9`): the row is recorded FAILED and **scoped to the density-only
wiring**, not recorded as a limit of the shipped OOD flag. The ruling's six required statements are
implemented in §6, *The SC3-d disposition*. §6 also carries a **verdict-sensitivity strip** — a
read-only table of what each SC3 verdict would have been at other bar values, added on the same
day's ruling about how the judgement-call bars are presented. **No bar was changed.**

---

## 4. FDR, stated at the scope it actually holds

*Source: `p14_fdr_report.jld2`. Substrate: the shared `p14_eval_pool.jld2`, n = 2000 at
`P14_EVAL_COUNTER = 2`, unstratified.*

### The α table — with its decided fraction on every row

| α_FDR | k\* | n_decided | **decided_fraction** | n_accepted | n_false_disc. | predicted FDP | realized FDP | 95 % band | within band |
|---|---|---|---|---|---|---|---|---|---|
| 0.01 | 311 | 1776 | **0.888** | 311 | 1 | 0.009917363921582958 | 0.003215434083601286 | [−0.0010957396597117783, 0.020930467502877692] | yes |
| 0.05 | 389 | 1776 | **0.888** | 389 | 20 | 0.049434169389304713 | 0.05141388174807198 | [0.027892170912349745, 0.07097616786625968] | yes |
| 0.10 | 426 | 1776 | **0.888** | 426 | 42 | 0.09898741069175962 | 0.09859154929577464 | [0.07062738085259579, 0.12734744053092345] | yes |
| 0.20 | 491 | 1776 | **0.888** | 491 | 106 | 0.1985881503174122 | 0.2158859470468432 | [0.16330073324987912, 0.23387556738494528] | yes |

The band is **derived**, not chosen: `predicted ± 1.96·√(predicted(1−predicted)/k*)` with
`P14_SC1B_BAND_Z = 1.96`. The partition is `n_decided = 1776`, `n_abstained = 224`
(`abstain_reason_counts = (ood_fired = 183, conformal_empty = 41)`). Sanity, persisted:
`max_abs_dev_sum_to_one = 2.220446049250313e-16` and `max_abs_dev_null_pooling =
2.3554969823605387e-16`, both against a 1e-12 tolerance.

### The honest claim, verbatim (`fdr_claim`)

> The posterior expected false discovery proportion is controlled at α_FDR over the **DECIDED SUBSET
> ONLY** — the batch minus abstentions. It is **not** controlled over the whole batch: an abstained
> item is neither a discovery nor a non-discovery and appears in neither the numerator nor the
> denominator. The guarantee is conditional on the model — calibrated posteriors and a correct class
> prior — and is a **posterior expectation, not a frequentist long-run rate**. A rule that abstains on
> most of a batch meets any level trivially, so the level and the decided fraction are ONE quantity
> and are never quoted apart.

`fdr_scope = :decided_subset_only` is a persisted symbol, and `P14FDRRow` is one declared concrete
`NamedTuple` type whose key set is asserted — a row missing `decided_fraction` **cannot be
constructed**, so the pairing is structural rather than a convention a reviewer must check.

### Assumption A2, written out as a derivation a referee can check

This is the **load-bearing assumption of the phase**. It is the entire justification for
abstain-then-sort and for the decided-subset scope, and `14-RESEARCH.md`'s Assumptions Log marks it
**High risk if wrong**. Two checkable lines:

1. The controlled quantity is
   `E[FDP | data] = (1/|R|) · Σ_{i ∈ R} P(H0_i | data)`, where `R` is the rejection set.
2. `R` is a **deterministic function of the observed data** — both the abstention filter and the
   prefix sort are computed from the data alone — hence `R` is **σ(data)-measurable** and pulls
   straight out of the conditional expectation. The identity therefore holds for **any**
   data-dependent selection of `R`, including one that filtered on OOD status or conformal ambiguity.

Unlike the frequentist case, where a data-dependent pre-selection is a genuine selective-inference
problem requiring correction, the Bayesian posterior quantity is immune **because the conditioning
has already happened**. (Persisted verbatim as `assumption_a2_statement`.)

**No independence across the batch is assumed.** `E[Σ 1{H0_i}] = Σ v_i` follows from linearity of
expectation alone. That is the standard reviewer question, and the answer travels inside the artifact
as `independence_note`.

### The ordering rule, and why the reverse order is wrong

**ABSTAIN FIRST, THEN SORT** (`ordering_note`). The reverse — sort the whole batch, then abstain from
some of the accepted items — is wrong: removing items from a set whose running mean was computed over
all of them changes both numerator and denominator in an uncontrolled way, so the residual accepted
set's running mean can **exceed** the level. Abstaining first also improves the premise, because it
routes away precisely the items whose posteriors are least likely to be right. The ordering is bound
by a test that a reordering turns red (`test_p14_decide.jl`, 146/146).

### SC1-c — prior sensitivity, REPORTED AND NOT GATED

`pi_sensitivity_gated = false` is persisted so this table can never later be read as a criterion that
was met. Every rung is a full re-run of the whole rule through the same `_p14_rule_at_prior` code path
`decide_coloc` uses; the reconstruction proof at the pool's own prior reproduced every one of the 2000
items with `max_abs_dev_posterior = 0.0` and `max_abs_dev_null = 0.0`.

| π_coloc | α | k\* | n_decided (frac) | predicted | realized |
|---|---|---|---|---|---|
| 0.05 | 0.01 | 240 | 1777 (0.8885) | 0.009871791417894052 | 0.0 |
| 0.05 | 0.05 | 317 | 1777 (0.8885) | 0.04988652159765778 | 0.0031545741324921135 |
| 0.05 | 0.10 | 357 | 1777 (0.8885) | 0.09967752941532872 | 0.025210084033613446 |
| 0.05 | 0.20 | 414 | 1777 (0.8885) | 0.19881838140767902 | 0.0966183574879227 |
| 0.10 | 0.01 | 271 | 1776 (0.888) | 0.009887673523070992 | 0.0 |
| 0.10 | 0.05 | 352 | 1776 (0.888) | 0.04958812640612884 | 0.022727272727272728 |
| 0.10 | 0.10 | 388 | 1776 (0.888) | 0.09868442291574546 | 0.059278350515463915 |
| 0.10 | 0.20 | 449 | 1776 (0.888) | 0.19990599183368657 | 0.1492204899777283 |
| **0.217** | 0.01 | 311 | 1777 (0.8885) | 0.009892142202187133 | 0.003215434083601286 |
| **0.217** | 0.05 | 389 | 1777 (0.8885) | 0.04933271013421822 | 0.05141388174807198 |
| **0.217** | 0.10 | 427 | 1777 (0.8885) | 0.09966332630692931 | 0.10070257611241218 |
| **0.217** | 0.20 | 492 | 1777 (0.8885) | 0.19892887488431868 | 0.21747967479674796 |
| 0.40 | 0.01 | 352 | 1772 (0.886) | 0.009917787526747676 | 0.022727272727272728 |
| 0.40 | 0.05 | 425 | 1772 (0.886) | 0.04969586188965847 | 0.1011764705882353 |
| 0.40 | 0.10 | 462 | 1772 (0.886) | 0.09901277414147118 | 0.17316017316017315 |
| 0.40 | 0.20 | 533 | 1772 (0.886) | 0.19874272068695037 | 0.28142589118198874 |
| 0.70 | 0.01 | 401 | 1770 (0.885) | 0.009858782852563825 | 0.06982543640897755 |
| 0.70 | 0.05 | 477 | 1770 (0.885) | 0.049625657130969995 | 0.18238993710691823 |
| 0.70 | 0.10 | 521 | 1770 (0.885) | 0.09944241408694038 | 0.2495201535508637 |
| 0.70 | 0.20 | 600 | 1770 (0.885) | 0.1997520196632641 | 0.3466666666666667 |

**The direction matters, and it fails in the dangerous direction.** *Below* the measured prevalence
the rule is merely conservative (π = 0.05, α = 0.10: predicts ~0.0997, realizes ~0.0252 — it discovers
less than it could). *Above* it, the rule is **anti-conservative and the FDR claim becomes false**
(π = 0.70, α = 0.10: predicts ~0.0994, realizes ~0.2495, two and a half times the level). **Nothing
throws.** The number simply stops meaning what its name says. The rung nearest the measured prior
(0.217 against 0.21655) reproduces the headline to the third decimal, and the decided fraction barely
moves across the whole sweep (0.885 to 0.8885) — so none of this is an artifact of the rule abstaining
its way to a friendly number.

**The prevalence a real batch actually has is not measured anywhere in this project.** Empirical-Bayes
estimation of π from the batch is deliberately **not** offered as a default: estimating the prior from
the same batch the guarantee is quoted over would make the guarantee circular.

### The prior-atom asymmetry

The θ prior carries clamp atoms and they are asymmetric: roughly 4.85 % of mass at ρ = −0.99 against
roughly 1.91 % at ρ = +0.99, about 2.5 to 1 **against** the exclusion end. The measured
`pi_class_masses` inherit those atoms, and **in Phase 14 they enter the PRIOR TERM of the class
posterior directly — a term Phase 13 never used**, because Phase 13 published log Bayes factors and
declined to introduce a prior at all. An asymmetry that was inert upstream is load-bearing here.
(`prior_atom_note`.)

### The unit

FDR is controlled **per image pair** (D-04). Per-region and per-tile FDR is out of scope — see §8.

---

## 5. The conformal hedge

*Source: `p14_conformal_report.jld2`.*

| Quantity | Value | Key |
|---|---|---|
| `q̂` (LAC, α_conf = 0.10) | **0.4498001628978108** | `qhat` |
| Realized coverage | **0.9115** (1823 covered of 2000) | `coverage`, `n_covered` |
| Pre-registered band (lower) | **0.8798753882025019** | `band_lower` |
| Gap | **+0.03162461179749809** | `coverage_gap` |
| Calibration set / evaluation set | n = 2000 @ counter 1 / n = 2000 @ counter 2, **both unstratified** | `n_cal`, `n_eval` |
| Calibration scores (min / median / max) | 9.873126405324228e-8 / 0.01705689278977429 / 0.9956241391144685 | `scores_summary` |
| `iteration_trigger_fired` | **false** — `P14_ITERATION_ALLOWANCE = 1` remains **unspent** | `iteration_trigger_fired` |
| Elapsed | 624.4879999160767 s | `elapsed_s` |

The band is derived from `P14_ALPHA_CONFORMAL = 0.1`, `P14_N_EVAL = 2000` and
`P14_SC1D_BAND_SD_MULTIPLIER = 3`; recomputed by hand it matches the persisted value to < 1e-12.
Coverage exceeding nominal (0.9115 > 0.90) is expected, not anomalous: split conformal guarantees
coverage ≥ 1 − α and the `ceil((n+1)(1−α))` order statistic makes it mildly conservative.

### Set-size distribution — and the qualifier that belongs beside the coverage number

| Status | Count | Rate |
|---|---|---|
| singleton | 1954 | 0.977 |
| **ambiguous** | **0** | **0.000** |
| empty | 46 | 0.023 |

**The hedge never actually hedges.** `ambiguous_rate` is **exactly 0.0** — at this `q̂` the LAC score
admits no 2-class set anywhere in 2000 items. The hedge's only departure from a bare point prediction
is the empty set (2.3 %). `vacuous_hedge = false` only because the label tests
`(ambiguous_rate + empty_rate) ≥ P14_AMBIGUOUS_RATE_FLOOR = 0.01`, and the empty sets alone carry it
past that floor. **On the ambiguous channel the hedge is silent, and coverage here is delivered by
classifier accuracy rather than by set-valued caution.** This is the singleton-dominant regime
14-RESEARCH §B.3 named in advance.

### Class masses — both draws unstratified (`p14_assert_unstratified` passed on each)

| Draw | exclusion | random | coloc |
|---|---|---|---|
| calibration (n = 2000) | 0.3145 | 0.4650 | 0.2205 |
| evaluation (n = 2000) | 0.3315 | 0.4545 | 0.2140 |
| expected prior π (`class_prior_used`) | 0.31272 | 0.47073 | 0.21655 |

Neither draw is equal thirds: the `random` and `coloc` masses sit 11–12 binomial SE from 1/3. Both
tracked the measured prior.

### Conformal status × OOD state (`cross_tab_status_by_ood`)

| | fired | clear | not_checked |
|---|---|---|---|
| singleton | 178 | 1776 | 0 |
| ambiguous | 0 | 0 | 0 |
| empty | 5 | 41 | 0 |

`ood_available = true` (threshold 167.5446329378043, `:density` wired), so no item resolved to
`:not_checked` and the D-06 all-abstain outcome did not arise. The shared pool's decisions are
`:not_coloc` 1381 / `:coloc` 395 / `:abstain` 224, with abstain reasons `:ood_fired` 183 and
`:conformal_empty` 41 — the 5 items that were both empty and OOD-fired are attributed to OOD under the
ABSTAIN-FIRST ordering.

### The guarantee is simulator-derived — what that means, exactly

`guarantee_basis = :simulator_derived_held_out_draws`, persisted in every artifact.

The calibration draws and the evaluation draws come from the **same forward simulator the evidence net
was trained on**, at reserved counters disjoint from the training stream. **The training joint
therefore equals the evaluation joint by construction**, so 0.9115 — and every other SC1 and SC3
number in this report — is a **well-specified-regime number that inherits the simulator's
misspecification in full**. These numbers say the arithmetic, the ordering and the calibration are
consistent with each other. They say **nothing** about how a real microscopy batch degrades any of it.
A skill of 0.90 on simulator draws is not a claim about a bench.

**The bound D-03 intended is ABSENT, not merely loose.** Those are the words, and the difference
matters: *loose* invites a reader to discount a weak bound, whereas nothing in this phase bounds the
layer's behaviour on real microscopy at all. The intended substrate turned out to hold no unsealed
physical ground truth (§7). The bound is **absent, not merely loose**.

---

## 6. Risk-coverage and abstention concentration

*Sources: `p14_riskcoverage_report.jld2` (SC3-a/b/c) and `p14_ood_arm_report.jld2` (SC3-d).*

### The three met statistics

| Statistic | Measured | Frozen bar | Verdict |
|---|---|---|---|
| SC3-a selective skill | **0.9017113386798696** | `P14_SKILL_FLOOR = 0.6` | MET |
| SC3-b `corspearman(coverage, selective_risk)` on [0.20, 1.00] | **0.9701942302111639** | `P14_SPEARMAN_FLOOR = 0.95` | MET |
| SC3-c `AUC_hard` | **0.9117431863910738** | `P14_AUC_HARD_FLOOR = 0.8` | MET |

Supporting quantities, all persisted: `aurc = 0.010420724333625658`,
`aurc_oracle = 0.0030436576034615016`, `aurc_random = 0.07809877211088907`,
`e_aurc = 0.007377066730164156`, `skill_denominator = 0.07505511450742756`,
`skill_degenerate = false`, `base_error_rate = 0.077` (154 wrong argmax calls of 2000),
`n_hard = 154`, `spearman_n_points = 1601`, `oracle_is_lower_envelope = true`,
`max_abs_dev_confidence = 0.0` (the swept scalar is pinned, by name, to the column 14-09 persisted).

**The RAW curve is shown over its full range, including the part the statistic excludes.**
`curve_is_raw_full_range = true`; the 2000 curve points span coverage **[0.0005, 1.0]** — the left
edge is a one-item denominator where the selective risk is exactly 0 or 1. SC3-b is computed on the
1601 points with coverage ≥ 0.20 and the excluded noisy region is persisted, plotted and shaded in
`spike/figures/p14_riskcoverage.png` (rendered **before** the gate was evaluated) so a reader can judge
the exclusion rather than take it on trust.

### THE FIVE BARS ARE JUDGEMENT CALLS WITH NO DERIVATION

`P14_JUDGEMENT_CALL_BARS = [:P14_SKILL_FLOOR, :P14_SPEARMAN_FLOOR, :P14_COVERAGE_FLOOR,
:P14_AUC_HARD_FLOOR, :P14_OOD_MARGIN_FLOOR]`. `bars_are_judgement_calls = true` and
`bar_is_judgement_call = true` are persisted; the banner in every runner prints the label
**mechanically from that frozen tuple**, so a bar cannot quietly stop being labelled one.

| Bar | Value | Status |
|---|---|---|
| `P14_SKILL_FLOOR` | 0.6 | **JUDGEMENT CALL** — no derivation |
| `P14_SPEARMAN_FLOOR` | 0.95 | **JUDGEMENT CALL** — no derivation |
| `P14_COVERAGE_FLOOR` | 0.2 | **JUDGEMENT CALL** — has a stated reason, no derivation of the value |
| `P14_AUC_HARD_FLOOR` | 0.8 | **JUDGEMENT CALL** — no derivation |
| `P14_OOD_MARGIN_FLOOR` | 0.5 | **JUDGEMENT CALL** — no derivation |

**The freeze evidence.** All five were frozen in commit
`a6c825867dbc786f7c3927df6760295ee0930c77` — a commit that contains **no Phase-14 result of any
kind**, checkable with `git show`. `coverage_floor_frozen_before_run = true` and
`coverage_floor_freeze_commit` record it inside the artifact. Plan 14-01 opened with a **blocking
decision checkpoint** on exactly these five bars and the user ruled them *accept-as-proposed* on
2026-08-03, before the freeze commit and before any Phase-14 code ran. **They were never relaxed.**
`forbidden_action` names the two moves that were not available: monotonizing the curve (running-max or
isotonic fit) and re-choosing the coverage floor after seeing where it becomes monotone.

**Only two of the five carry any stated reasoning at all.** The selective-skill statistic is
normalized between the random and oracle curves computed on the same items, which makes the bar immune
to the base error rate — but *the normalization is derived and the floor is not*. The coverage floor
has a recorded reason (at coverage k/n one hard item moves the selective risk by 1/k, so a trend
statistic at the left edge measures sampling noise) — but *that reason does not derive 0.20*.

**The licensing standard, recorded from `.planning/STATE.md`.** This project has now recorded **four**
separate cases of a bar measuring something other than what it named — SC1g's component-vs-total, the
n = 2-against-271 real arm, the Wald-labelled-Wilson sizing (DEF-12-04), and the Phase-12 Stage-1
control ceiling. The honest framing is that pre-registration is what *surfaces* these, and all four
were caught before corrupting a result — but the count is now high enough that *"we pre-registered
it"* no longer settles an argument on its own. **The standard any amendment is held to is evidence,
independent of the machinery under audit, that the bar measures something other than what it names.**
No such evidence exists for any of the five bars here, so none of them moves, and **the SC3-d
shortfall below is reported, not re-tuned.**

### Verdict sensitivity — what each verdict would have been at a different bar (NO BAR WAS CHANGED)

**Read this framing before the tables.** **No bar was changed, and none can be.** The frozen values
in `spike/p14/consts.jl` — byte-unchanged against HEAD, re-asserted at the moment this strip was
computed — are the **only** values that gate anything in Phase 14, and every verdict reported
anywhere in this document is the verdict at those frozen values. The rungs below are
**hypothetical**. They exist so that a reader can see whether a conclusion **hinges on an underived
number or is robust to it**, which is the actual defence against a judgement-call bar. **This is
presentation, not relaxation.** It is the user's ruling of 2026-08-04 on how the five bars are
framed, and it changes nothing measured: nothing was re-simulated, nothing re-trained, and the four
gated statistics were **loaded** from the artifacts above and compared against a sweep.

*Source: `p14_bar_sensitivity.jld2`, written by `spike/p14/run_p14_bar_sensitivity.jl` — a runner
that gates nothing (`gates_nothing = true`, `no_bar_changed = true`,
`recomputed_anything_gated = false`) and reads the frozen bars from `consts.jl` by name. Its verdict
at each frozen rung is **pinned by assertion** to the `sc3a_met` / `sc3b_met` / `sc3c_met` /
`sc3d_met` flags the gated runners persisted, so the strip cannot be scoring under a different rule
from the one Phase 14 was scored under.*

**What a one-sided bar's sensitivity actually is, said plainly rather than dressed up.** Every SC3
gate has the form `statistic ≥ bar`, so the bar value at which the verdict flips **is** the measured
statistic, by construction. The informative content of the strip is therefore the **distance from
the frozen bar to that flip point** — the room the verdict had. It is not new evidence and must not
be quoted as any.

`*` marks the frozen bar. `MET` / `NOT` is the verdict that hypothetical bar would have produced.

| Gate | Measured | 0.30 | 0.40 | 0.50 | 0.60 | 0.70 | 0.80 | 0.85 | 0.90 | 0.91 | 0.95 | Frozen → verdict | Flips at | Distance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **SC3-a** selective skill | 0.9017 | MET | MET | MET | **MET\*** | MET | MET | MET | MET | NOT | NOT | `P14_SKILL_FLOOR = 0.60` → **MET** | 0.9017 | **+0.3017** |

| Gate | Measured | 0.80 | 0.85 | 0.90 | 0.93 | 0.95 | 0.96 | 0.97 | 0.98 | 0.99 | 1.00 | Frozen → verdict | Flips at | Distance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **SC3-b** Spearman | 0.9702 | MET | MET | MET | MET | **MET\*** | MET | MET | NOT | NOT | NOT | `P14_SPEARMAN_FLOOR = 0.95` → **MET** | 0.9702 | **+0.0202** |

| Gate | Measured | 0.60 | 0.70 | 0.75 | 0.80 | 0.85 | 0.90 | 0.91 | 0.92 | 0.95 | 0.99 | Frozen → verdict | Flips at | Distance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **SC3-c** `AUC_hard` | 0.9117 | MET | MET | MET | **MET\*** | MET | MET | MET | NOT | NOT | NOT | `P14_AUC_HARD_FLOOR = 0.80` → **MET** | 0.9117 | **+0.1117** |

| Gate | Measured | 0.00 | 0.05 | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 | 0.70 | 0.80 | Frozen → verdict | Flips at | Distance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **SC3-d** abstain margin | −0.1290 | NOT | NOT | NOT | NOT | NOT | NOT | **NOT\*** | NOT | NOT | NOT | `P14_OOD_MARGIN_FLOOR = 0.50` → **NOT MET** | −0.1290 | **−0.6290** |

**What the strip says, row by row.**

- **SC3-a and SC3-c are far from any boundary.** Skill would have passed a bar anywhere up to 0.90 —
  half the admissible range above the frozen 0.60 — and `AUC_hard` up to 0.91. Neither verdict hinges
  on the choice of its underived number; a referee substituting a different-but-reasonable floor gets
  the same answer.
- **SC3-b is the one PASS with modest room, and this is the honest place to say so.** The Spearman
  clears the frozen 0.95 by **0.0202**. A floor of 0.98 — not an unreasonable number for a
  monotone-trend criterion — would have turned SC3-b red. The verdict is real, but it is not robust
  to a materially stricter choice, and any citation of SC3-b should carry that.
- **SC3-d does not fail marginally, and no bar choice rescues it.** The measured margin is
  **negative**: the `noise` arm abstains *less* than the in-distribution arm. It fails at **every**
  rung of the sweep including **0.00** — i.e. it would fail a bar that asked only for "not worse than
  in-distribution". Nothing about the frozen 0.50 produced this result, and lowering it would not
  have avoided it. (Why it fails, and the scope of that failure, is the next subsection.)

**`P14_COVERAGE_FLOOR` is different in kind, and is treated differently.** It is not a pass/fail bar:
it is the coverage at which the raw risk-coverage curve is **truncated** before SC3-b's Spearman is
taken. Its sensitivity is therefore a genuine **re-derivation** of the statistic at other
truncations, computed from the raw full-range curve that `p14_riskcoverage_report.jld2` already
persists (`curve_coverage`, `curve_selective_risk`). It **was** computable from the persisted
artifact; nothing here is fabricated. The re-derivation is **pinned by assertion** at the frozen
truncation to the persisted `spearman` (agreement to < 1e-12) and to `spearman_n_points = 1601`
before any other rung is reported, so every other row is produced by the same computation SC3-b was
scored with.

| Truncation | Curve points in range | Re-derived Spearman | vs the frozen `P14_SPEARMAN_FLOOR = 0.95` |
|---|---|---|---|
| 0.02 | 1961 | 0.9368 | **NOT MET** |
| 0.05 | 1901 | 0.9427 | **NOT MET** |
| 0.10 | 1801 | 0.9523 | MET |
| 0.15 | 1701 | 0.9616 | MET |
| **0.20 \*** (frozen) | **1601** | **0.9702** | **MET** |
| 0.25 | 1501 | 0.9780 | MET |
| 0.30 | 1401 | 0.9846 | MET |
| 0.40 | 1201 | 0.9927 | MET |
| 0.50 | 1001 | 0.9919 | MET |
| 0.60 | 801 | 0.9977 | MET |
| 0.70 | 601 | 0.9989 | MET |

**And this row is the one a referee should look at hardest, so it is stated without hedging.** The
Spearman is **monotonically increasing in the truncation** across the whole sweep bar one inversion
(0.40 → 0.50), which is exactly what the coverage floor's stated rationale predicts: the lower the
truncation, the more of the one-item-denominator region enters the correlation, and that region is
sampling noise. **SC3-b would have been NOT MET at a floor of 0.05 or 0.02**, and it is MET at every
floor from 0.10 upward. The frozen 0.20 is therefore **inside the passing region but not at its
edge** — one rung of margin below it (0.15, 0.10) and every rung above it also pass. A reader is
entitled to ask whether 0.20 was chosen to land there; the answer is on the record and checkable:
`P14_COVERAGE_FLOOR` was frozen in `a6c825867dbc786f7c3927df6760295ee0930c77`, a commit that contains
no Phase-14 result of any kind, so it could not have been. That is the whole value of the
pre-registration, and this table is what makes the value legible rather than asserted.

### SC3-d — NOT MET

*Source: `p14_ood_arm_report.jld2`. Five matched arms of n = 1000 at `P14_OOD_COUNTER = 3`; one
in-distribution, one per in-repo `OOD_FAMILIES` member at `level = :strongest`
(`family_level = 4`, `family_level_is_strongest_existing_rung = true`). `arms_matched_n = true`.
`qhat_recalibrated_here` — the 14-09 `q̂` was **loaded**, and the density threshold
167.5446329378043 is **pinned by an executable assertion** to 14-09's
(`ood_threshold_pinned_to_14_09 = true`).*

| Arm | abstain rate | n_abstain | ood_fired | ood_not_checked | conformal_empty | conformal_ambiguous | disagreement_and_ood |
|---|---|---|---|---|---|---|---|
| in_distribution | 0.129 | 129 | 112 | 0 | 17 | 0 | 0 |
| texture | 0.999 | 999 | 999 | 0 | 0 | 0 | 0 |
| **noise** | **0.000** | **0** | **0** | 0 | 0 | 0 | 0 |
| optics | 0.387 | 387 | 376 | 0 | 11 | 0 | 0 |
| background | 0.989 | 989 | 989 | 0 | 0 | 0 | 0 |

| Family | margin = rate(family) − rate(ID) | vs floor 0.5 |
|---|---|---|
| texture | **+0.870** | MET |
| **noise** | **−0.129** | **NOT MET** |
| optics | **+0.258** | NOT MET |
| background | **+0.860** | MET |

**Gated headline: `margin = −0.129` (minimum over families, weakest family `:noise`) against
`P14_OOD_MARGIN_FLOOR = 0.5`. `sc3d_met = false`. SC3-d IS NOT MET.**

The minimum, not the mean, is gated **by design** (`gated_reading = :minimum_over_families`): the mean
of the four is +0.465, which would also have missed the floor, but quoting it would have hidden that
one family produces *zero* detections. A report quoting a mean would be quoting a different statistic
under this statistic's name.

**Was the margin driven by `ood_fired` or by `ood_not_checked`?** By **`ood_fired`** — and for `noise`
by neither, because nothing fired at all. `any_not_checked = false` is persisted: **no item in any arm
resolved `:not_checked`**, so every OOD-triggered abstention here is a *fired detector*, never a
default abstention. The `margin_driver` field reports `dominant = :ood_fired` for texture, optics and
background. For `noise` the reported `dominant` symbol is `:ood_not_checked`, and that is an **argmax
artefact, not a finding**: all five per-reason deltas for that family are ≤ 0
(`ood_fired = −0.112`, `ood_not_checked = 0.0`, `conformal_empty = −0.017`, rest 0.0), so `argmax`
selects the only non-negative one. The substantive fact is that the `noise` arm abstained on **nothing
at all**: 0 fired, 0 empty sets, 1000 singletons, 1000 `:clear`.

**The empty-set × `ood_state` cross-tab** (rows `singleton / ambiguous / empty`, columns
`fired / clear / not_checked`):

| Arm | empty ∧ fired | empty but CLEAR | fired but not empty | n_empty | n_fired |
|---|---|---|---|---|---|
| in_distribution | 1 | 17 | 111 | 18 | 112 |
| texture | 19 | 0 | 980 | 19 | 999 |
| noise | 0 | 0 | 0 | 0 | 0 |
| optics | 14 | 11 | 362 | 25 | 376 |
| background | 23 | 0 | 966 | 23 | 989 |

`:ambiguous` is empty in every arm — 14-09's LAC finding carried forward. The two signals share no
assumption (the density null needs the simulator to be right about densities; an empty conformal set
needs only exchangeability), and **both readings are reported without editorialising one into the
other**: where the density channel fires, *every* empty set sits on an item it also flagged (texture
19 of 19, background 23 of 23) — genuine corroboration. In-distribution, **17 of 18 empty sets land on
items the density null called CLEAR** — the exchangeability channel is seeing something the density
channel is not, on data the density channel considers ordinary.

**What SC3-d's failure actually says.** The claim is *"abstention concentrates on data the model was
not trained on."* Measured against four deliberately injected misspecifications, that is **true for
three and false for the fourth**, and the fourth is the one named `noise`. That is not a coincidence:
`channels_wired = (:density,)` and `channels_not_wired = (:pp, :noise)`. **The one family that slips
past is the family the unwired channel exists for.** The single-wired-channel limit was declared in
`spike/p14/pools.jl` and printed in the runner banner **before any arm was built**; SC3-d has now
converted it from a stated limit into a measured one with a number attached.

**Two things this does not license, and neither was done.** (1) Wiring the noise channel now to clear
the bar — that is authoring the experiment after seeing its result. `P14_ITERATION_ALLOWANCE = 1` has
exactly one pre-declared trigger, split-conformal coverage falling *below* the SC1-d band, and 14-09
measured 0.9115 *above* a band of 0.8798753882025019, so the allowance is unspent and unavailable
here (`iteration_allowance_applies_here = false`). (2) Relaxing `P14_OOD_MARGIN_FLOOR`, dropping the
`noise` family, or lowering the grid rung. `spike/p14/consts.jl` is byte-unchanged against HEAD, and
**nothing was tuned to try to clear this gate.**

**No misspecification was authored for this phase.** The families are the in-repo `OOD_FAMILIES`
(texture, noise, optics, background) from `spike/validation/ood.jl` at the strongest existing rung;
an acceptance-criterion grep proves no family is defined in the runner. Arms are matched item by item
— identical θ, λ, image size and stream position, differing only in which generator produced the
pixels.

### The SC3-d disposition — ruled by the user on 2026-08-04

*Binding record: `.planning/phases/14-decision-and-abstention-layer/14-SC3D-RULING.md`, committed
`f4934a9`. The ruling answers the `resume-signal`'s question — "if any SC row FAILED, say how you
want it recorded" — and it is reproduced in force here.*

**Phase 14 closes NEGATIVE-BUT-USEFUL on this row. SC3-d is recorded as a named limit SCOPED TO THE
DENSITY-ONLY WIRING — not as a limit of the shipped OOD flag.** The six points below are the ruling's
own instruction for how this row must be stated, and each is implemented rather than summarized.

**1. The failure, unsoftened.** **SC3-d FAILED.** Margin **−0.129** against the frozen floor
**0.50**; `sc3d_met = false`; weakest family `:noise` with **0 of 1000** items flagged at the
strongest existing rung against an in-distribution abstain rate of 0.129. The abstain-rate table
above stands unchanged, the per-family margins stand unchanged, and no bar was relaxed, re-derived,
re-scoped or swapped for a neighbouring statistic. The verdict-sensitivity strip above shows the
failure is not a boundary effect: the margin is negative, so **no** non-negative bar would have
passed it.

**2. The scope, stated immediately beside the failure.** The measurement was taken on a
**density-only** detector: `ood_channels_wired = (:density,)`, `ood_channels_not_wired = (:pp,
:noise)` (`spike/p14/decide.jl:583-584`, `:979-980`; `spike/p14/pools.jl:487,506,563-564`;
`decide.jl:286,293` switches the posterior-predictive channel off with `with_pp = false`). **One of
the shipped flag's three channels.** The shipped **fused** detector attains **AUC 1.0 on the same
`noise` family** — read from `artifacts/amended_v2/grid_8/gate_report_8.jld2`, key **`report[:ood]`**,
loaded programmatically by `spike/p14/run_p14_bar_sensitivity.jl` and persisted to
`p14_bar_sensitivity.jld2` (`gate8_density_auc`, `gate8_noise_auc`, `gate8_fused_auc`). Four
misspecification levels per family; `levels = 4`, `n = 200`, `id_threshold = 179.1367781263274`,
`fused_threshold = 7.128118439466407`:

| Family | `density_auc` (**what Phase 14 wired**) | `noise_auc` | `fused_auc` (**shipped**) |
|---|---|---|---|
| `texture` | 1.0 / 1.0 / 1.0 / 0.99475 | 1.0 / 1.0 / 1.0 / 1.0 | 1.0 / 1.0 / 1.0 / 1.0 |
| **`noise`** | **0.0 / 0.0 / 0.0 / 0.0** | 1.0 / 1.0 / 1.0 / 1.0 | **1.0 / 1.0 / 1.0 / 1.0** |
| `optics` | 0.58 / 0.6755 / 0.656875 / 0.66975 | 1.0 / 1.0 / 1.0 / 1.0 | 1.0 / 1.0 / 1.0 / 1.0 |
| `background` | 0.55925 / 0.95875 / 0.0 / 0.960125 | 1.0 / 1.0 / 1.0 / 1.0 | 1.0 / 1.0 / **0.0** / 1.0 |

The density channel is **exactly blind** to the `noise` family — AUC 0.0 at every level, i.e. worse
than chance. That is the *by-design* blind spot the noise channel exists to close:
`test/gate/misspec.jl:32-33` states the reported detector is the OR-fused density∨noise pair because
*"the image-noise channel closes the detector-noise blind spot of a correlation-only summary."* **The
scope is DENSITY-ONLY, not NOISE-FAMILY-ONLY**: `optics` shows the same effect more mildly — density
0.58–0.66975 against a fused 1.0 — and its 0.387 abstain rate above is consistent with that. The
limit is a property of the **wiring**, not of one family.

**3. The honest one-line form, quoted from the ruling.** **"the Phase-14 decision layer's abstention
is blind to detector-noise misspecification because it wires one of the shipped flag's three
channels — not because the shipped flag is blind to it."**

**4. What is NOT claimed, stated explicitly.** **The fused detector's SC3-d margin is UNMEASURED.**
Phase 14 never ran the fused detector through `decide_coloc`, so it is **unknown** whether SC3-d would
pass with all channels wired, and **this report does not assert that it would**
(`fused_sc3d_margin_measured = false` is persisted in `p14_bar_sensitivity.jld2`). What *is* known is
that the channel SC3-d needed was present in the shipped artifact and absent from the decision layer.
**Phase 15 tests the fused detector independently** — its OOD arm runs `gate_ood_roc`, i.e. the fused
density∨noise pair, and sweeps a `noise` axis (Phase 15 **D-02 / D-06 / D-07**). Phase 15's
three-valued verdict for that axis is the number that settles it, not anything in this report.

**5. `P14_ITERATION_ALLOWANCE` is UNSPENT, and why.** The allowance is 1 and its **single**
pre-declared trigger is *conformal coverage falling below the SC1-d band*. That trigger **did not
fire**: SC1-d measured **0.9115** against a band of **0.8798753882025019** — above, not below
(`iteration_allowance_applies_here = false`, `iteration_trigger_fired = false`). The allowance is
therefore unavailable here and remains unspent. It could not have been spent on SC3-d in any case:
`P14_ITERATION_TRIGGER` says in terms that it is *never* spent on relaxing a bar after a result and
specifically never on any member of `P14_JUDGEMENT_CALL_BARS`.

**6. The two rejected alternatives, recorded with their reasons, so the refusal to re-wire is visible
as a choice rather than an omission.**

| Rejected alternative | Reason on the record |
|---|---|
| **Route SC3-d to Phase 15 as a pre-registered rebuild with a fresh bar** | Phase 15's OOD arm runs `gate_ood_roc` — the **fused** density∨noise detector — so a fresh bar there would measure a *differently equipped* detector and answer a different question. It would also collide with Phase 15's **D-14**, which lists "add a detector channel for the failing mechanism" under Deferred and forbids it as an in-phase response. |
| **Re-wire `(:pp, :noise)` into the Phase-14 decision layer and re-measure SC3-d** | It would very likely clear the bar — on the AUC table above, near-certainly — **and that is precisely why it is refused.** Improving the instrument after seeing the result, in order to pass a frozen bar, is the pattern this project has already had to correct four times (`.planning/STATE.md`, the four-amendments observation). |
| *(also rejected)* **Record SC3-d as a general limit of the OOD flag** | Factually too pessimistic and falsifiable by opening the shipped gate report: `fam_best_auc[:noise] = 1.0`, `fused_auc[:noise] = 1.0` at all four levels. The limit belongs to the wiring, not to the flag. |

---

## 7. The real-data check (D-03a)

*Source: `p14_real_images_report.jld2`. Substrate: `test/test_images`, read-only.*

Six committed microscopy TIFFs, two conditions, opened through Phase 13's frozen `load_real` /
`real_tif`. Operative pair `P13_REAL_CHANNEL_PAIR = [2, 3]` (green/red, the two target proteins),
passed explicitly at every call site; c1 is the DAPI/Hoechst nuclear counterstain and enters no
comparison, though it is still opened and digested because the unit of record is the file. Two
decision units: `positive_as_sample` (control = negative) and `negative_as_sample`.

| file | ch | in pair | decision | trigger | OOD state | score / threshold | conformal | P(coloc) |
|---|---|---|---|---|---|---|---|---|
| `test/test_images/positive/positive_c1.tif` | 1 | false | **abstain** | `:ood_fired` | `:fired` | 703.2995398773389 / 167.5446329378043 | singleton `[:coloc]` | 0.9928014238079867 |
| `test/test_images/positive/positive_c2.tif` | 2 | true | **abstain** | `:ood_fired` | `:fired` | 703.2995398773389 / 167.5446329378043 | singleton `[:coloc]` | 0.9928014238079867 |
| `test/test_images/positive/positive_c3.tif` | 3 | true | **abstain** | `:ood_fired` | `:fired` | 703.2995398773389 / 167.5446329378043 | singleton `[:coloc]` | 0.9928014238079867 |
| `test/test_images/negative/negative_c1.tif` | 1 | false | **abstain** | `:ood_fired` | `:fired` | 703.2995398773389 / 167.5446329378043 | singleton `[:random]` | 7.997639204997099e-5 |
| `test/test_images/negative/negative_c2.tif` | 2 | true | **abstain** | `:ood_fired` | `:fired` | 703.2995398773389 / 167.5446329378043 | singleton `[:random]` | 7.997639204997099e-5 |
| `test/test_images/negative/negative_c3.tif` | 3 | true | **abstain** | `:ood_fired` | `:fired` | 703.2995398773389 / 167.5446329378043 | singleton `[:random]` | 7.997639204997099e-5 |

**The prediction held.** `predicted_behaviour` was written into the runner **before** the run and
persisted verbatim beside the outcome: *"ABSTAIN on both specimens (13-REPORT limit D: density 703.30
vs threshold 167.54, 4.198x)"* → `observed_decisions = [:abstain, :abstain]` → `prediction_held =
true`. Abstaining here is the **designed** behaviour under D-05/D-06, not a failure.

**The two misspecification signals disagree on this substrate.** Both conformal sets are singletons —
the exchangeability channel saw nothing unusual — while the density null fires at 4.198× its
operating point. Recorded, not editorialised. Cross-method channel on both units: `classical_call =
RANDOM`, `costes_p_sample = 0.004975124378109453` (`costes_call = :coloc`),
`basis = :ghat_rho_true_scale`. The fused OOD trigger fires first, so the classical disagreement never
reaches the decision — the D-05 asymmetry doing its job.

### The corpus record — computed from the manifest at run time, not transcribed

| field | value |
|---|---|
| `corpus_images_available` | **0** (total bytes 0 across all 32 rows) |
| `corpus_status` | **`:unfetched_by_design`** |
| `tier_breakdown` | `simulated-secondary` = 30, `physical-primary` = 2 |
| `split_breakdown` | `dev` = 14, `eval` = 16, `sealed_holdout` = 2 |
| `physical_truth_available` | **false** (0 unsealed physical-primary rows with bytes) |
| `d03_bound_status` | **`:absent_not_loose`** |
| `sealed_holdout_read` / `corpus_images_opened` | false / 0 |
| `manifest sha256` | `0dfc7e50e5680befd5436a77a01ecce170e1579c0a832841cb1a24c76b615701` |

The corpus is **unfetched by design** — never "gone" and never "empty". The ignore rules were
**located by name at run time** (`.gitignore` lines 459, 461, 463, 465, 466, with `corpus/data/*` at
461) and their text persisted, so the "by design" reading rests on the file rather than on a claim.
The line range 14-CONTEXT's D-03a quotes had already drifted; the runner does not depend on it.

**The only `physical-primary` rows are Phase 16's sealed holdout, and they were not read** — not in a
script, not in a test, not behind a flag. There is no second blind set, so consuming one would
irreversibly burn the blind evaluation on the very hypothesis it exists to evaluate. The remaining 30
rows are `simulated-secondary` (the CBS benchmark — *another simulator*), so **a completed fetch would
still not have produced the physical check D-03 wanted.**

### The D-03a substitution and its rationale

The check moved from the manifest's physical-primary anchors to the six committed TIFFs, following
the precedent Phase 13 set at `spike/p13/real_images.jl:71,145` — adopted *"precisely so that Phase
16's blind corpus stays blind"* (`13-D15-AMENDMENT.md`). `substitution_record` carries the from, the
to, the precedent, the amendment, the rationale, `sealed_holdout_read = false` and
`metadata_only = true`.

### The withdrawn figure

D-03 recorded a miscoverage floor derived from the manifest — a ratio of one to thirty-one, roughly
three percent. **It is WITHDRAWN, not merely superseded.** It was computed from the wrong n (the
manifest is partitioned dev 14 / eval 16 / sealed_holdout 2, not "30 open") **and** on rows that are
not physical truth. It appears as **no numeral** anywhere in the Phase-14 runners, artifacts or this
report; `spike/test/test_p14_decoupling.jl` bans the literal from Phase-14 executable code, which is
why it is written here in prose.

### What the six TIFFs are

Six TIFFs in two conditions, **with no colocalization ground-truth labels**, is an **illustration, not
a coverage claim.** It bounds nothing tightly. `is_coverage_claim = false` and `no_rate_computed =
true` are required keys of the atomic save — an artifact that failed to say so is never written at all
— and no rate, fraction or percentage is computed over these six images anywhere in the runner.

**A naming correction that must travel with any quote of this table:** `positive/` and `negative/` are
the original package's *biological test conditions*, **not colocalization labels**. The `negative`
pair is a positively correlated pair, not an anti-correlated one; nothing here may treat it as an
exclusion example.

**Read-only discipline:** `readonly_digest_before == readonly_digest_after ==
eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276`, asserted equal **after** the
artifact was persisted.

---

## 8. Named limits, carried into the manuscript

**Enumerated from code — `p14_named_limits()` in `spike/p14/result.jl`, eight items, printed by
running it rather than transcribed.** Every artifact embeds the same tuple (`named_limits`,
`named_limits_count = 8`), so a limit cannot be silently dropped by a report writer.

1. The decision layer is built in the spike research lane and is **NOT shipped in v2.0** (D-01): it
   has a src-shaped signature and subtypes the shipped `AbstractColocResult`, but src-shaped is not
   in-src, and no byte of `src/` is edited by this phase.
2. The conformal guarantee is **SIMULATOR-DERIVED** (D-02/D-03) and inherits the simulator's
   misspecification in full; the intended real-data bound is **ABSENT, not merely loose**, because the
   substrate it was to be computed on holds no unsealed physical ground truth and its bytes are
   unfetched (D-03a).
3. The real-data check is six committed microscopy TIFFs in two conditions (D-03a) — an
   **ILLUSTRATION**, not a coverage claim. It bounds nothing tightly and must never be quoted as
   though it did.
4. FDR is controlled over the **DECIDED SUBSET ONLY** — the batch minus abstentions — and the decided
   fraction must be quoted beside α, because a rule that abstains on most of a batch hits any α
   trivially. The guarantee is a posterior expectation conditional on the model, not a frequentist
   long-run rate.
5. The per-tile map is **UNCONTROLLED DISPLAY**, never an FDR-controlled call (D-04): Phase 12
   returned NO on calibrated per-region uncertainty and the shipped per-tile map carries no
   uncertainty field, so there is nothing per-tile to control against.
6. **SC1 is AMENDED** by D-02 and **SC2 is AMENDED** by D-05. Any report citing a Phase-14 result must
   cite the original ROADMAP wording alongside the amendment.
7. The four SC3 bars and the 0.20 coverage floor are **JUDGEMENT CALLS** with no derivation
   (`P14_JUDGEMENT_CALL_BARS`). They were ratified by the user before any Phase-14 result existed, and
   a missed bar is REPORTED, not re-tuned.
8. The underlying net is an **epoch-4 checkpoint of an early-overfitting run** (13-REPORT limit C),
   and every Phase-14 number inherits that limit unchanged.

**Four further limits, added by this report:**

**A. The net is an epoch-4 checkpoint of an untuned, early-overfitting run.** 13-REPORT limit C:
best validation risk **0.13625 at epoch 4**, rising to **0.38632**, roughly a **10× train/validation
gap**, with the recipe deliberately not tuned. **Every number in this report inherits it** — the
posteriors the FDR rule sorts, the confidence scalar the risk-coverage curve sweeps, and the LAC
scores the conformal quantile is taken over are all outputs of that checkpoint.

**B. Only the summary-density OOD channel is wired in this lane.** `channels_wired = (:density,)`,
`channels_not_wired = (:pp, :noise)`. The shipped flag is an OR-fusion over three channels; here a
`:clear` means **one** channel says clear, not three. This cuts both ways and is recorded as such: it
is not an excuse for the SC3-d shortfall, because the layer as measured is the layer as it exists —
but a reader comparing this margin to a three-channel detector's is not comparing like with like.
**Per the user's SC3-d ruling of 2026-08-04 this is the limit SC3-d's failure is scoped to** — see
§6, *The SC3-d disposition*, for the loaded `density_auc` / `noise_auc` / `fused_auc` table from
`artifacts/amended_v2/grid_8/gate_report_8.jld2 → report[:ood]` that measures it, and for the
explicit statement that the **fused** detector's SC3-d margin is **unmeasured**.

**C. The τ probe ran on a dirty tree.** `provenance.probe_repo_dirty_at_run = true`, inherited from
Phase 13 and **surfaced rather than asserted away**. τ = 0.15 itself is four-way agreed
(`four_way_agreed = true`) across the net artifact, the constants file, the probe and the gate, and
the probe's own `consts_sha` (`640ec90e…`) is asserted to **diverge** from both pinned literals rather
than being widened into a disjunction that would hide the divergence.

**D. Per-region FDR is out of scope, and not merely unimplemented.** Phase 12 returned NO on
calibrated per-region uncertainty, and `LocalColocMap` carries `grid`, `tiles`, `delta_rho`,
`ood_flag`, `meta` — **no uncertainty field at all**. There is nothing to control against, so per-tile
output is **uncontrolled display**. A grep asserts the runners produce no such number anywhere.
Claiming per-region FDR would repeat the counting-a-closed-negative-as-a-success failure this project
has already corrected once.

---

## 9. Forward-pointing items

1. **The conformal hedge's ambiguous channel is silent, and nothing in this phase tests whether it
   would speak under misspecification.** *Owner: Phase 15.* `ambiguous_rate = 0.0` across 2000
   in-distribution items and across all five OOD arms; the only non-singleton output the hedge ever
   emitted was the empty set. Evidence: `p14_conformal_report.jld2` → `set_size_counts`,
   `p14_ood_arm_report.jld2` → `set_size_counts_per_arm`. The pre-declared remedy
   (`P14_ITERATION_TRIGGER`: switch LAC → APS) is **unspent** and was never triggered, because it fires
   only on a coverage *shortfall* and coverage was above band. A future phase choosing to exercise APS
   must do so as a *new* pre-registration, not by spending this allowance retroactively.
2. **The empty conformal set is the one signal in this phase that does not depend on the simulator
   being right about densities.** *Owner: Phase 15/16.* Evidence: the in-distribution cross-tab, where
   17 of 18 empty sets land on items the density null called CLEAR. It is the natural candidate for a
   misspecification channel that would survive the D-03 absence, and it is currently un-gated.
3. **SC3-d's shortfall names a concrete build item: wire the `(:pp, :noise)` channels.** *Owner: Phase
   15, as a pre-registered change with its own bar.* Evidence: `channels_not_wired`, and the `noise`
   arm's 0-of-1000 flag rate. Doing it inside Phase 14 would have been authoring the experiment after
   seeing its result; doing it in Phase 15 against a freshly frozen bar is legitimate.
4. **The class prior is the simulator's mix, and no real batch's prevalence is measured anywhere in
   this project.** *Owner: Phase 16.* Evidence: the SC1-c table, where π = 0.70 turns a nominal 0.10
   into a realized 0.2495. A Phase-16 external evaluation is the first opportunity to observe a real
   prevalence; until then every quoted FDR level must travel with both its decided fraction and this
   table.
5. **The D-03 real-data bound remains ABSENT.** *Owner: Phase 16, or a roadmap change.* Evidence:
   `physical_truth_available = false`, `d03_bound_status = :absent_not_loose`. Either Phase 16's
   sealed holdout supplies it, or the corpus is extended with real graded anchors — which is a roadmap
   decision, not an executor's.
6. **The Phase-4 `SPEEDUP_GATE` red still masks every later include block.** *Owner: user decision,
   `deferred-items.md` D-13-A.* Carried forward unchanged from 13-REPORT §11 item 6. See §10.
7. **A CRLF/LF blob-sha discrepancy in the recorded `consts_git_blob_sha`.** *Owner: a future
   provenance plan; cosmetic, recorded so nobody mistakes it for tampering.* The artifacts record
   `consts_git_blob_sha = 26ac1e7b2b5e85badd3e2dc9f99238fff21bc034`, which is
   `git hash-object --no-filters` over the working-tree bytes. The blob actually committed is
   `c1dbd185435f975ff237dae6e9d9eca48cd57134`, after git's CRLF→LF normalization on Windows. The two
   describe the same content; a reader running `git rev-parse HEAD:spike/p14/consts.jl` and comparing
   to the artifact would otherwise see a mismatch and suspect an edit. The sha256 of the raw file,
   `4e8ee4be102ccbadd526d01302d1587279884e26bddfc9157336e03086d86ca5`, is the unambiguous digest and
   agrees everywhere.
8. **A PROTECTIVE prediction for Phase 15's `noise` axis, recorded BEFORE that phase runs.** *Owner:
   Phase 15 (D-02 / D-06 / D-07).* Phase 15 sweeps a `noise` axis on the **fused** detector and
   assigns it PROTECTIVE / LATE / SILENT-BUT-SAFE. On the AUC evidence in §6 —
   `gate_report_8.jld2 → report[:ood]`, `noise_auc[:noise] = 1.0` and `fused_auc[:noise] = 1.0` at all
   four levels — **PROTECTIVE is expected.** This is written down *now*, from the SC3-d ruling of
   2026-08-04, so that a `LATE` outcome there would be a genuine surprise with recorded evidence
   against it, rather than a result read after the fact.
9. **`fused_auc[:background]` drops to 0.0 at level 3 while `noise_auc[:background]` is 1.0 there —
   flagged, not diagnosed.** *Owner: Phase 15, before the `background` and the new pure-offset
   autofluorescence axes (D-02a) are read.* Evidence: the §6 scope table, and
   `fire_rate[:background] = [1.0, 1.0, 0.0, 1.0]` in the same artifact. The z-score OR-fusion is
   **not monotone in its channels**, so a fused AUC below its best channel is possible — but a clean
   **zero** exactly at one rung, with 1.0 either side of it, is worth a look before that family is
   reported. **Not diagnosed here; recorded so it is not discovered inside a reported run.**

---

## 10. Reproduction

### The standing warning about the aggregate suite

**`julia --project=spike spike/test/runtests.jl` exits 1 at the Phase-4 `SPEEDUP_GATE`, and that
failure masks every later include block — including all ten Phase-14 test files.** The aggregate suite
**must not be used as a Phase-14 gate**. Per-file runs are the only reliable signal
(`deferred-items.md` D-13-A; `13-REPORT.md` §11 item 6). Every Phase-14 verification in this report
was performed per-file.

### Per-file test runs — measured, all ten exit 0

| File | Top-level testset | Pass / Total | Exit |
|---|---|---|---|
| `spike/test/test_p14_consts.jl` | P14 Tier-1 pre-registration (D-07) | 273 / 273 | 0 |
| `spike/test/test_p14_decoupling.jl` | P14 decoupling: `src/`, deps, the sealed holdout and the phase's own greps | 51 / 51 | 0 |
| `spike/test/test_p14_provenance.jl` | P14 tau provenance (D-07) | **78 / 78** (was 77 before the sensitivity runner landed — see below) | 0 |
| `spike/test/test_p14_posterior.jl` | P14 three-class posterior (D-04, class order, prior re-derivation) | 48 / 48 | 0 |
| `spike/test/test_p14_fdr.jl` | P14 Bayesian-FDR prefix rule (SC1-a) | 55 / 55 | 0 |
| `spike/test/test_p14_conformal.jl` | P14 hand-rolled split conformal (SC1-e, D-02) | 83 / 83 | 0 |
| `spike/test/test_p14_fuse.jl` | P14 D-05 asymmetric fusion and D-06 three-valued OOD (SC2-a, SC2-b) | 222 / 222 | 0 |
| `spike/test/test_p14_result.jl` | P14 result and batch types (D-01, D-04, D-05, D-06) | 121 / 121 | 0 |
| `spike/test/test_p14_decide.jl` | P14 decide: the composed pipeline, abstain-first, at fixture scale | 146 / 146 | 0 |
| `spike/test/test_p14_pools.jl` | P14 pools: unstratified draws, the Phase-13 basis null, and D-06 end to end | 102 / 102 | 0 |
| **Total** | | **1179 / 1179** | **all 0** |

```bash
for f in consts decoupling provenance posterior fdr conformal fuse result decide pools; do
  julia --project=spike spike/test/test_p14_$f.jl || echo "FAILED: $f"
done
```

**Why the provenance count moved from 77 to 78, recorded rather than quietly restated.** The first
assembly of this report measured 77. `test_p14_provenance.jl:73` builds `P14_LANE_FILES` by
**globbing** `spike/p14/*.jl` rather than by an enumerated list, precisely so that a file added in a
later wave is scanned automatically instead of needing the list edited — the comment there says so.
Adding `spike/p14/run_p14_bar_sensitivity.jl` (the verdict-sensitivity strip, §6) therefore adds
**exactly one** assertion: the lane-wide guard that no Phase-14 source writes τ as a hardcoded
literal, now applied to the new file, which passes. **No test changed, no test was added by hand, and
the delta is +1 by construction.**

### The four reported runners

```bash
julia --project=spike -t auto spike/p14/run_p14_conformal.jl     # SC1-d, SC1-e; writes p14_conformal_report.jld2 AND the shared p14_eval_pool.jld2 (624.5 s)
julia --project=spike -t auto spike/p14/run_p14_fdr_check.jl     # SC1-b, SC1-c; consumes the shared pool (8.3 s)
julia --project=spike -t auto spike/p14/run_p14_riskcoverage.jl  # SC3-a/b/c; consumes the shared pool (33.3 s)
julia --project=spike -t auto spike/p14/run_p14_ood_arm.jl       # SC3-d; needs the ~52 MB gitignored spike/data/cache/p13 pool (3947.5 s)
julia --project=spike -t auto spike/p14/run_p14_real_images.jl   # SC1-f; read-only on the six committed TIFFs (22.7 s)
```

**And one runner that is NOT a reported runner and gates nothing** — the verdict-sensitivity strip of
§6, added on the user's 2026-08-04 ruling about how the judgement-call bars are presented:

```bash
julia --project=spike spike/p14/run_p14_bar_sensitivity.jl       # §6 strip; READ-ONLY, gates nothing; writes p14_bar_sensitivity.jld2 (6.2 s)
```

It re-simulates nothing and re-measures nothing: it loads the already-measured `skill`, `spearman`,
`auc_hard` and SC3-d `margin` and compares them against a sweep of hypothetical bar values, re-derives
SC3-b's Spearman at alternative truncations of the persisted raw curve (pinned to the persisted value
at the frozen truncation), and reads `report[:ood]` out of the shipped
`artifacts/amended_v2/grid_8/gate_report_8.jld2` **read-only**. It asserts, before computing anything,
that `spike/p14/consts.jl` is byte-unchanged against HEAD.

**Order matters:** `run_p14_conformal.jl` produces `p14_eval_pool.jld2`, the shared n = 2000
evaluation set that the FDR and risk-coverage runners **consume rather than redraw**, so the three
score the *same items at the same counter*. A missing pool fails loudly naming 14-09.

**Execute on the main working tree — no git worktrees.** `run_p14_ood_arm.jl` fits its density null on
the ~52 MB gitignored `spike/data/cache/p13/` pool, which does not exist inside a fresh worktree and
was destroyed once already this milestone.

**Every runner is persist-before-assert.** The artifact (and, for risk-coverage, the figure) is
written **before** any gate is evaluated, so a failing gate still leaves a complete self-describing
report on disk. This was confirmed empirically, not designed and hoped for: `run_p14_ood_arm.jl`
exited non-zero on the SC3-d `AssertionError` and the 738,407-byte artifact — including
`sc3d_met = false` — was on disk regardless. **The evidence of an honest failure survived the
failure.**

### Environment

Julia 1.12.6, 32 threads, CPU-only (`julia_version` and `nthreads` are persisted in every artifact).
No package was installed at any point in this phase; the frozen sixteen dependencies are asserted by
name on every test run.

---

## 11. One-paragraph summary for a reader in a hurry

Phase 14 built `decide_coloc`: a batch decision layer emitting {coloc / not / ABSTAIN} at a user-set
Bayesian FDR, in the spike research lane, against a **src-shaped** signature. It is **NOT shipped in
v2.0**. Of nineteen pre-registered validation rows, **seventeen PASS, one is REPORTED-NOT-GATED, and
one FAILS.** The FDR rule tracks its own predicted false discovery proportion inside a binomial 95 %
band at all four pre-registered levels, over a **DECIDED SUBSET ONLY** of 1776 of 2000 items
(decided fraction 0.888). Split-conformal coverage is 0.9115 against a derived band of
0.8798753882025019 — but the hedge emits **zero** two-class sets, so that coverage is delivered by
classifier accuracy rather than by set-valued caution. Abstention orders and concentrates well
(skill 0.9017, Spearman 0.9702, AUC_hard 0.9117, all against underived judgement-call bars). **SC3-d
is NOT MET at −0.129 against a floor of 0.50**, because the one wired OOD channel is blind to the one
injected misspecification family the unwired channel exists for — a limit declared in code before any
arm was built, now measured. On the only real microscopy this project may touch, the layer abstains,
as predicted in writing beforehand — an **illustration, not a coverage claim**. And every gated number
above is a **well-specified-regime** number: train-joint equals eval-joint by construction on a shared
simulator pool, so the real-data bound D-03 intended is **absent, not merely loose**.
