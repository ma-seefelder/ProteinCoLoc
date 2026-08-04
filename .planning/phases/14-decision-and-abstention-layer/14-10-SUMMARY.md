---
phase: 14-decision-and-abstention-layer
plan: 10
subsystem: decision-and-abstention-layer
status: complete
tags: [bayesian-fdr, sc1-b, sc1-c, prior-sensitivity, decided-subset, pre-registered-gate]
requires:
  - "14-09 (spike/p14/p14_eval_pool.jld2 — the SHARED n_eval=2000 evaluation set at P14_EVAL_COUNTER)"
  - "spike/p14/fdr.jl (p14_bayes_fdr, p14_fdr_over_decided with its required n_total)"
  - "spike/p14/decide.jl (_p14_rule_at_prior, _p14_prior_at, _p14_no_cross_method, _p14_headline_line)"
  - "spike/p14/consts.jl (P14_ALPHA_FDR_GRID, P14_PI_COLOC_GRID, P14_N_EVAL, P14_EVAL_COUNTER — frozen before any Phase-14 result existed)"
provides:
  - "SC1-b: a measured realized-vs-predicted false discovery proportion at four pre-registered levels, each with its decided fraction"
  - "SC1-c: the 20-row prior-sensitivity table, REPORTED and explicitly not gated"
  - "spike/p14/p14_fdr_report.jld2 — the alpha table, the pi table, the honesty block and the provenance"
affects:
  - "any report quoting a Phase-14 FDR number: the scope, the prior sensitivity and the per-pair unit now travel inside the artifact"
tech-stack:
  added: []
  patterns: ["persist-before-assert", "atomic .jld2 save with reopen integrity check", "run-time decoupling proof", "locked-threshold banner", "one declared row type so decided_fraction cannot be omitted", "honesty block as REQUIRED artifact keys"]
key-files:
  created:
    - "spike/p14/run_p14_fdr_check.jl (961 lines)"
    - "spike/p14/p14_fdr_report.jld2 (gitignored, reproducible)"
  modified: []
decisions: [D-04, D-05, D-07]
validates: ["SC1-b", "SC1-c", "SC1-a"]
metrics:
  duration: "~25 min wall (the reported run itself: 10.9 s)"
  completed: 2026-08-04
---

# Phase 14 Plan 10: Bayesian-FDR Control Check (SC1-b) and Prior Sensitivity (SC1-c) Summary

**SC1-b PASSES at all four pre-registered levels.** On the shared 2000-item evaluation set the
realized false discovery proportion among decided-and-accepted items falls inside the binomial 95 %
band of the rule's own predicted posterior-expected FDP at every level of `P14_ALPHA_FDR_GRID`. The
decided fraction is **0.888 (1776 / 2000)** at every level and travels with every number below.

**SC1-c is REPORTED and is not gated.** The 20-row prior sweep shows the guarantee is genuinely
prior-sensitive: at a batch that is 70 % coloc the realized FDP is **0.347** where the rule predicts
0.200. That is the documented failure mode, measured rather than argued.

## SC1-b: the four-row table (the plan's required reproduction, with the decided fraction column)

| α_FDR | k\* | n_decided | **decided_fraction** | predicted FDP | realized FDP | band (95 %) | within band |
|---|---|---|---|---|---|---|---|
| 0.01 | 311 | 1776 | **0.888** | 0.009917 | 0.003215 (1/311) | [−0.001096, 0.020930] | yes |
| 0.05 | 389 | 1776 | **0.888** | 0.049434 | 0.051414 (20/389) | [0.027892, 0.070976] | yes |
| 0.10 | 426 | 1776 | **0.888** | 0.098987 | 0.098592 (42/426) | [0.070627, 0.127347] | yes |
| 0.20 | 491 | 1776 | **0.888** | 0.198588 | 0.215886 (106/491) | [0.163301, 0.233876] | yes |

**Verdict, plainly: SC1-b is MET.** No row exceeds its upper band, which is the pre-registered
falsification condition. The margin is not uniform and that is worth saying: at α = 0.20 the
realized 0.2159 sits at about 82 % of the way from the predicted value to the upper band, so the
level at which the rule is most permissive is also the level where realized and predicted are
furthest apart. The band is derived — `predicted ± 1.96·√(predicted(1−predicted)/k*)` — not chosen.

Acceptance thresholds and the decision-theoretic reading, reported and not gated:
`t* = 0.0739 / 0.4135 / 0.7410 / 0.9189` and `cost_ratio = 0.0798 / 0.7051 / 2.861 / 11.333`.

### The composite-null decomposition of each accepted set (reported beside the number)

| α_FDR | mean p_random | mean p_exclusion | share random | share exclusion |
|---|---|---|---|---|
| 0.01 | 0.009917 | 0.000000 | 1.0000 | 0.0000 |
| 0.05 | 0.049430 | 0.000005 | 0.9999 | 0.0001 |
| 0.10 | 0.098970 | 0.000017 | 0.9998 | 0.0002 |
| 0.20 | 0.198509 | 0.000079 | 0.9996 | 0.0004 |

**Essentially all of the accepted null mass is `random`, not `exclusion`.** The rule's residual risk
is confusing noise with signal, not confusing segregation with colocalization. That is the more
benign of the two, and it is only visible because `v` was reported split rather than as a total.

### The section-E.1 sanity assertion

`p_coloc + p_random + p_exclusion = 1` to a maximum absolute deviation of **2.22e-16** across all
2000 items, and the stored composite null equals `p_random + p_exclusion` to **2.36e-16**. Both are
persisted (`max_abs_dev_sum_to_one`, `max_abs_dev_null_pooling`) with the 1e-12 tolerance they were
scored against.

### The partition

`decided = 1776`, `abstained = 224` (`:ood_fired` 183, `:conformal_empty` 41). This reconciles
exactly with 14-09's recorded pool. Abstention happened FIRST; only the 1776 survivors were sorted.

## SC1-c: the prior-sensitivity sweep — REPORTED, NOT GATED (no bar, and none is implied)

`pi_sensitivity_gated = false` is persisted in the artifact so the table cannot later be read as a
criterion that was met. All 20 rows (5 rungs of `P14_PI_COLOC_GRID` × 4 levels of
`P14_ALPHA_FDR_GRID`):

| π_coloc | α | k\* | n_decided (frac) | predicted | realized |
|---|---|---|---|---|---|
| 0.05 | 0.01 | 240 | 1777 (0.8885) | 0.009872 | 0.000000 |
| 0.05 | 0.05 | 317 | 1777 (0.8885) | 0.049887 | 0.003155 |
| 0.05 | 0.10 | 357 | 1777 (0.8885) | 0.099678 | 0.025210 |
| 0.05 | 0.20 | 414 | 1777 (0.8885) | 0.198818 | 0.096618 |
| 0.10 | 0.01 | 271 | 1776 (0.8880) | 0.009888 | 0.000000 |
| 0.10 | 0.05 | 352 | 1776 (0.8880) | 0.049588 | 0.022727 |
| 0.10 | 0.10 | 388 | 1776 (0.8880) | 0.098684 | 0.059278 |
| 0.10 | 0.20 | 449 | 1776 (0.8880) | 0.199906 | 0.149220 |
| **0.217** | 0.01 | 311 | 1777 (0.8885) | 0.009892 | 0.003215 |
| **0.217** | 0.05 | 389 | 1777 (0.8885) | 0.049333 | 0.051414 |
| **0.217** | 0.10 | 427 | 1777 (0.8885) | 0.099663 | 0.100703 |
| **0.217** | 0.20 | 492 | 1777 (0.8885) | 0.198929 | 0.217480 |
| 0.40 | 0.01 | 352 | 1772 (0.8860) | 0.009918 | 0.022727 |
| 0.40 | 0.05 | 425 | 1772 (0.8860) | 0.049696 | 0.101176 |
| 0.40 | 0.10 | 462 | 1772 (0.8860) | 0.099013 | 0.173160 |
| 0.40 | 0.20 | 533 | 1772 (0.8860) | 0.198743 | 0.281426 |
| 0.70 | 0.01 | 401 | 1770 (0.8850) | 0.009859 | 0.069825 |
| 0.70 | 0.05 | 477 | 1770 (0.8850) | 0.049626 | 0.182390 |
| 0.70 | 0.10 | 521 | 1770 (0.8850) | 0.099442 | 0.249520 |
| 0.70 | 0.20 | 600 | 1770 (0.8850) | 0.199752 | 0.346667 |

**Read this table before quoting any FDR number.** Three things it says:

1. **The guarantee is only as good as the prior's match to the batch, and the direction matters.**
   Below the measured prevalence the rule is merely CONSERVATIVE — at π = 0.05 it predicts 0.0999
   and realizes 0.0252, so it discovers less than it could. Above it, the rule is
   **ANTI-CONSERVATIVE and the FDR claim becomes false**: at π = 0.70 and a nominal 0.10 the
   realized proportion is 0.2495, two and a half times the level. Nothing in the run throws; the
   number simply stops meaning what its name says.
2. **The rung nearest the measured prior reproduces the headline.** π = 0.217 against the measured
   0.21655 gives 0.003215 / 0.051414 / 0.100703 / 0.217480 versus the headline 0.003215 / 0.051414
   / 0.098592 / 0.215886. The two agree to the third decimal, which is the corroboration the sweep
   exists to provide and not a coincidence — see the reconstruction proof below.
3. **The decided fraction barely moves** (0.885 to 0.8885 across the whole sweep), so none of the
   above is an artifact of the rule abstaining its way to a friendly number at some rungs.

The prevalence a real batch actually has is **not measured anywhere in this project**, so this table
is the honest statement of what the FDR claim depends on. Empirical-Bayes estimation of π from the
batch is deliberately not offered as a default: estimating the prior from the same batch the
guarantee is quoted over would make the guarantee circular.

## The honesty block travels inside the artifact

Eight string keys are REQUIRED keys of the atomic save, so an artifact missing one is never written
at all: `fdr_claim`, `assumption_a2_statement`, `ordering_note`, `independence_note`, `prior_note`,
`prior_atom_note`, `unit_note`, `amendment`. Their substance:

- **`fdr_claim`** — control is over the DECIDED SUBSET ONLY; an abstained item is in neither the
  numerator nor the denominator; the guarantee is conditional on the model and is a posterior
  expectation, not a frequentist long-run rate.
- **`assumption_a2_statement`** — the two-line derivation, stated as a derivation a referee can
  check: `E[FDP | data] = (1/|R|)·Σ_{i∈R} P(H0_i | data)`, and `R` is a deterministic function of the
  observed data, hence σ(data)-**measurable**, hence pulls out of the conditional expectation — so a
  data-dependent pre-selection needs no selective-inference correction, unlike the frequentist case.
- **`ordering_note`** — abstain first, then sort. The reverse order can leave the residual accepted
  set's running mean ABOVE α, because removing items moves both numerator and denominator.
- **`independence_note`** — `E[Σ 1{H0_i}] = Σ v_i` needs only linearity of expectation, not
  independence across the batch. The standard reviewer question, with a clean answer.
- **`prior_note`** — the simulator's class mix is not any real batch's prevalence; the sensitivity
  table above is the answer; empirical Bayes is deliberately not the default.
- **`prior_atom_note`** — the θ prior carries clamp atoms (~4.85 % of mass at ρ = −0.99 against
  ~1.91 % at +0.99, ~2.5:1 against the exclusion end) which the measured `pi_class_masses` inherit,
  and which in Phase 14 enter the PRIOR TERM of the posterior directly — a term Phase 13 never used.
- **`unit_note`** — FDR is controlled **PER IMAGE PAIR** (D-04). Per-region and per-tile FDR is out
  of scope: the tile map carries no uncertainty field and Phase 12 returned NO on calibrated
  per-tile uncertainty, so there is nothing to control against. The runner produces no such number
  anywhere, and a grep asserts it.
- **`amendment`** — SC1 amended by D-02, SC2 by D-05; any report citing this must cite the original
  ROADMAP wording alongside it.

## Tasks

**Task 1 — realized-vs-predicted FDP over the frozen α grid.** Commit `f1a764a`. Created
`spike/p14/run_p14_fdr_check.jl` on the 14-09 skeleton: banner, step-0 run-time decoupling proof,
locked-threshold banner, atomic save with reopen integrity check, smoke redirect, guarded entry
point. The shared pool is CONSUMED, never redrawn: absence fails loudly naming 14-09, and the
recorded counter, size and class masses are asserted against the frozen constants (the masses are
re-derived from the pool's own stored labels through `p14_assert_unstratified`, not read off the
stored summary).

**Task 2 — the π sweep and the honest-claim block.** Commit `cc1cf19`. Every rung is a full re-run
of the whole rule through `_p14_rule_at_prior`, the same code path `decide_coloc` uses, so the sweep
cannot drift from the headline as a second simpler loop would.

## Deviations from Plan

### Auto-fixed / structural choices, none requiring a decision

**1. [Rule 2 - Missing critical check] The π sweep is preceded by an executable reconstruction proof**

- **Found during:** Task 2
- **Issue:** The SC1-b table is computed from the pool's RECORDED per-item fields, while the π sweep
  must recompute posteriors and re-run the rule from the stored log Bayes factors. Those are two
  routes to one set of numbers. A sweep that quietly disagreed with the table at the very rung where
  they should coincide would still produce a perfectly well-formed curve, with nothing saying which
  route was wrong. That is precisely the Phase-12 class of failure (an O(1) inconsistency invisible
  to structural verification).
- **Fix:** `_p14_fdr_assert_reconstruction` re-runs `_p14_rule_at_prior` at the pool's OWN prior and
  asserts, on every one of the 2000 items, that the re-run reproduces the three class posteriors and
  the composite null to 1e-12 and the conformal status and abstain/decide partition exactly. It
  **passed** with `max|Δp| = 0.0` and `max|Δv| = 0.0` and `n_decided = 1776`, and the deviations are
  persisted as `pi_sensitivity_reconstruction`.
- **Files modified:** `spike/p14/run_p14_fdr_check.jl`
- **Commit:** `cc1cf19`

**2. [Rule 3 - Structural] A single declared row type for the α table**

- **Found during:** Task 1
- **Issue:** T-14-06 requires that no FDR number can exist without its decided fraction. Building
  rows as ad-hoc NamedTuples makes that a convention a reviewer must check.
- **Fix:** `P14FDRRow` is one declared concrete NamedTuple type; `_p14_fdr_row` asserts the key set
  and converts. A row missing `decided_fraction`, or carrying an extra key, cannot be constructed.
  The band fields are `Union{Missing,Float64}` so a level with no discoveries records `missing`
  rather than a `NaN` that would read as "computed and undefined".
- **Files modified:** `spike/p14/run_p14_fdr_check.jl`
- **Commit:** `f1a764a`

**3. [Rule 3 - Blocking] `_p14_save_report` was not importable**

- **Found during:** Task 1
- **Issue:** The plan says to reuse 14-09's `_p14_save_report`. It is defined inside
  `run_p14_conformal.jl`, which also defines a zero-argument `main`; including one runner from
  another would put two `main` methods in one namespace where the last include silently wins.
- **Fix:** A runner-local `_p14_fdr_save_report` carrying the identical `_p13_gate_save_report`
  idiom, with the reason for the duplication documented in its docstring. This matches the Phase-13
  precedent, where the save helper is carried per runner. The eval-pool PATH is likewise named
  locally and then pinned by asserting the pool's own recorded counter and size.
- **Files modified:** `spike/p14/run_p14_fdr_check.jl`
- **Commit:** `f1a764a`

### Assumption Drift (advisory)

**The plan assumed the LAC hedge would be the interesting part of the abstention story; the measured
result is that OOD carries it.** Planned framing: abstention is a fused signal. Actual: of 224
abstentions, 183 are `:ood_fired` and 41 are `:conformal_empty`, and 14-09 already recorded that the
hedge emits zero 2-class sets at this q̂. Why it matters: the "decided subset" whose FDR is
controlled is defined mostly by the density OOD channel, which is the one channel wired in this
lane. Non-blocking, and it does not touch any SC1-b or SC1-c number.

## Verification Performed

All commands were run and their real output observed.

- `julia --project=spike -t auto spike/p14/run_p14_fdr_check.jl` — **exit 0**, 10.9 s (Task 1 run:
  exit 0, 9.8 s).
- `spike/p14/p14_fdr_report.jld2` exists with mtime **09:28:09Z**, newer than the Task-1 run start
  **09:27:46Z**; rewritten again by the Task-2 run started 09:31:19Z.
- `alpha_table` has exactly 4 rows and every row carries `decided_fraction` — "rows ok".
- `fdr_scope` is `:decided_subset_only` — verified by assertion.
- `pi_sensitivity` has exactly **20** rows; `pi_sensitivity_gated == false`; all ten SC1-c keys
  present — "SC1-c keys ok".
- `assumption_a2_statement` is a non-empty `String` containing `measurable` — "a2 measurable ok".
- `unit_note` contains `PER IMAGE PAIR` — "unit_note ok".
- The printed headline carries `FDR SCOPE: DECIDED SUBSET ONLY` together with the α, the decided
  fraction as `1776/2000 (88.8%)` and the measured realized FDP, on ONE line per level, emitted
  through the same `_p14_headline_line` a batch object would print.
- `grep -v '^\s*#' spike/p14/run_p14_fdr_check.jl | grep -cE '\b0\.01\b|\b0\.05\b|\b0\.20\b'` — **0**.
- `grep -v '^\s*#' spike/p14/run_p14_fdr_check.jl | grep -ci 'per_region\|per_tile\|local_coloc_map'`
  — **0**.
- `git status --porcelain -- src spike/Project.toml spike/Manifest.toml corpus artifacts` — **empty**.
- `git diff HEAD -- spike/p14/consts.jl` — **empty**; no constant was modified, before or after
  seeing any number.
- All 10 `spike/test/test_p14_*.jl` files run **individually** — 10/10 PASS (never via the aggregate
  suite, which exits 1 early at the Phase-4 speedup gate). This includes
  `test_p14_decoupling.jl`, whose lane-wide source scan now covers the new runner.

## Known Stubs

None.

## Threat Flags

None. This runner opens no network path, reads no image, writes only under `spike/p14/`, adds no
dependency, and consumes no random stream at all.

## Self-Check: PASSED

- `spike/p14/run_p14_fdr_check.jl` (961 lines, ≥ 200 required; contains `P14_ALPHA_FDR_GRID` 11×,
  `p14_eval_pool` 2×, `p14_fdr_over_decided` 2×) — FOUND
- `spike/p14/p14_fdr_report.jld2` — FOUND
- `spike/p14/p14_eval_pool.jld2` (consumed, unmodified) — FOUND
- commit `f1a764a` (Task 1) — FOUND
- commit `cc1cf19` (Task 2) — FOUND
