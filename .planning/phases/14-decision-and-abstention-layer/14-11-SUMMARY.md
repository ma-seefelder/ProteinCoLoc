---
phase: 14-decision-and-abstention-layer
plan: 11
subsystem: decision-and-abstention-layer
status: complete
tags: [risk-coverage, selective-prediction, sc3-a, sc3-b, sc3-c, judgement-call-bars, pre-registered-gate]
requires:
  - "14-09 (spike/p14/p14_eval_pool.jld2 — the SHARED n_eval=2000 evaluation set at P14_EVAL_COUNTER)"
  - "spike/p14/posterior.jl (p14_confidence — the SAME scalar the LAC hedge thresholds)"
  - "spike/p14/decide.jl (p14_load_bundle — the D-07 inherited tau, prior and provenance)"
  - "spike/validation/ood.jl (roc_auc — the in-repo tie-aware Mann-Whitney AUC, reused not reimplemented)"
  - "spike/validation/figures.jl (the headless CairoMakie backend)"
  - "StatsBase.corspearman (never hand-rolled, PATTERNS S6)"
  - "spike/p14/consts.jl (P14_SKILL_FLOOR, P14_SPEARMAN_FLOOR, P14_COVERAGE_FLOOR, P14_AUC_HARD_FLOOR — frozen in a6c8258 before this runner existed)"
provides:
  - "SC3-a: selective skill normalized between an oracle and a random reference curve computed on the SAME data"
  - "SC3-b: a Spearman trend statistic on the pre-registered coverage range, with the RAW full-range curve persisted beside it"
  - "SC3-c: AUC_hard — abstention score against would-have-been-wrong items, scored on the simulator's true label"
  - "spike/p14/p14_riskcoverage_report.jld2 — the three curves, the three statistics, the four-bar honesty block and the provenance"
  - "spike/figures/p14_riskcoverage.png — the raw curve and both references, rendered BEFORE the gate"
affects:
  - "any report quoting a Phase-14 SC3 number: the bars-are-judgement-calls label, the coverage-floor freeze commit and the well-specified-regime caveat now travel inside the artifact"
tech-stack:
  added: []
  patterns: ["persist-before-assert", "figure-before-gate", "atomic .jld2 save with reopen integrity check", "run-time decoupling proof", "locked-threshold banner printed mechanically from P14_JUDGEMENT_CALL_BARS", "honesty block as REQUIRED artifact keys", "degenerate denominator guarded before the division"]
key-files:
  created:
    - "spike/p14/run_p14_riskcoverage.jl (1008 lines)"
    - "spike/p14/p14_riskcoverage_report.jld2 (gitignored, reproducible)"
    - "spike/figures/p14_riskcoverage.png (gitignored, reproducible)"
  modified: []
decisions: [D-03, D-05, D-07]
validates: ["SC3-a", "SC3-b", "SC3-c"]
metrics:
  duration: "~35 min wall (the reported run itself: 40.0 s)"
  completed: 2026-08-04
---

# Phase 14 Plan 11: Selective Risk-Coverage (SC3-a / SC3-b / SC3-c) Summary

**All three pre-registered gates are MET on the shared 2000-item evaluation set.** Selective skill
**0.9017** against a floor of 0.60; Spearman trend **0.9702** against 0.95; AUC_hard **0.9117**
against 0.80. Every one of those four floors is a **judgement call with no derivation** — the
artifact says so mechanically, in a required key, next to the numbers.

## The measured numbers

| Statistic | Measured | Frozen bar | Verdict |
|---|---|---|---|
| **SC3-a** selective skill | **0.9017113386798696** | `P14_SKILL_FLOOR = 0.60` | **MET** |
| **SC3-b** `corspearman(coverage, selective_risk)` on coverage ∈ [0.20, 1.00] | **0.9701942302111639** | `P14_SPEARMAN_FLOOR = 0.95` | **MET** |
| **SC3-c** `AUC_hard` | **0.9117431863910738** | `P14_AUC_HARD_FLOOR = 0.80` | **MET** |

Supporting quantities, all persisted:

| Quantity | Value |
|---|---|
| `n_eval` | 2000 (the SHARED 14-09 pool, consumed, never redrawn) |
| base error rate | 0.077 — 154 of 2000 argmax calls wrong |
| `AURC` (measured) | 0.010420724333625658 |
| `AURC_oracle` | 0.0030436576034615016 |
| `AURC_random` | 0.07809877211088907 |
| `E_AURC` = AURC − AURC_oracle | 0.007377066730164156 |
| curve points / coverage range | 2000 points over **[0.0005, 1.0]** — the RAW full range |
| Spearman subset | 1601 of 2000 points, coverage ∈ [0.20, 1.00] |
| hard items | 154 of 2000 |
| oracle below measured | `true` (the envelope behaves as an envelope) |

**One sentence, as the plan's output section requires: all four SC3 bars — 0.60, 0.95, 0.20 and
0.80 — are judgement calls with no derivation, frozen in `spike/p14/consts.jl` at commit
`a6c825867dbc786f7c3927df6760295ee0930c77` before this runner existed, and they may not be relaxed
after a result.**

## What was actually verified, and how

Every claim below was produced by running the command, not by reading the code.

- `julia --project=spike -t auto spike/p14/run_p14_riskcoverage.jl` → **exit 0**, all three
  assertions passed after the artifact and the figure were already on disk.
- The plan's Task-2 automated verify (the twelve required artifact keys +
  `bars_are_judgement_calls`) → **exit 0**, printing
  `skill=0.9017113386798696 spearman=0.9701942302111639 auc_hard=0.9117431863910738`.
- A second artifact check asserted from disk: `d["bars"] == (skill_floor = 0.60,
  spearman_floor = 0.95, coverage_floor = 0.20, auc_hard_floor = 0.80)`; `minimum(curve_coverage)
  = 0.0005 < 0.20` (the raw curve spans below the floor); `spearman_n_points == count(c -> c >=
  0.20, curve_coverage)` (the Spearman really was computed on the declared subset only);
  `bars_note` contains both `judgement` and `four`; `pointwise_monotonicity_asserted == false`;
  the freeze commit matches. → **exit 0**.
- Grep acceptance criteria on `spike/p14/run_p14_riskcoverage.jl`:
  `function roc_auc|function corspearman` → **0**; `roc_auc(` → **1**; `corspearman` → **4**;
  comment-stripped `\b0\.60\b|\b0\.95\b|\b0\.80\b|\b0\.20\b` → **0** (no bar inlined);
  `JUDGEMENT CALL` → **10** occurrences, including one beside each of the four bars in the
  locked-threshold banner.
- `spike/figures/p14_riskcoverage.png` and `spike/p14/p14_riskcoverage_report.jld2` both exist with
  mtimes newer than the run start, and both are written **before** the first gate is evaluated —
  confirmed by the runner's own step ordering ([6/7] figure, [7/7] persist, then headline, then
  assertions) and by inspecting the rendered figure.
- All ten `spike/test/test_p14_*.jl` files run individually → **exit 0** each.
- `git status --porcelain -- src spike/Project.toml spike/Manifest.toml corpus artifacts` → empty.
- `git diff --quiet HEAD -- spike/p14/consts.jl` → **byte-unchanged**.

## How each statistic was constructed (and what was deliberately not done)

**The sweep is over one scalar, and it is the hedge's own.** `κ_i = p14_confidence(p_i) = max_y
p̂(y|Z_i)` — recomputed BY NAME from the three stored class posteriors and then pinned to the
column 14-09 persisted (max absolute deviation **0.0**). That is the same quantity the LAC hedge
thresholds, so the curve and the conformal layer describe **one** ordering rather than two. A
risk-coverage curve **cannot** be drawn over the D-05 fusion — an OOD flag is binary and a
conformal set size is an integer, so the fusion admits no total order to sweep — and the artifact
carries that sentence as `score_note` (T-14-37).

**The curve is raw and full-range.** Swept over the descending unique values of κ with `>=` at each
threshold, matching the tie convention `roc_auc` already uses in this repository. The first point
sits at coverage 0.0005, where the denominator is one item and the selective risk is exactly 0 or
1. That noise is **shown, not hidden**: it is persisted, it is plotted, and the figure shades the
[0.20, 1.00] region so a reader can see precisely what the SC3-b statistic excludes.

**Both reference curves are computed on the same data.** Oracle = every correct item first (the
lower envelope); random = a **recorded** permutation from `shuffle(p14_rng(P14_EVAL_COUNTER),
1:n)`, whose seed, salt, counter and recipe are all persisted, so the "random" reference is
regenerable bitwise rather than being an unseeded draw. Normalizing between both is what makes
SC3-a immune to the base error rate — a bar on a raw AURC would mostly be measuring how easy the
batch was. The three integrals are only comparable over a common coverage span, so the runner
**asserts** that the measured curve and the references start at the same coverage (within `5/n`)
rather than assuming it.

**SC3-b is a rank correlation, never a pointwise claim.** `pointwise_monotonicity_asserted = false`
is a persisted key. The empirical curve is a step function whose denominator is the number of
selected items, so at coverage k/n one hard item moves the risk by 1/k; asserting pointwise
monotonicity would be asserting something the object cannot satisfy (T-14-36).

**SC3-c is not circular.** `hard_i` is *the argmax call would have been wrong* against the
**simulator's** true label, which exists for every item — not against the abstention, the conformal
set, or anything else the decision layer emits. It is literally the same vector as the curve's loss
`ℓ`, computed once, because two definitions could only create a way for them to disagree. The
abstention score is `1 − κ`, matching `roc_auc`'s higher-is-more-positive convention, and the AUC
comes from the in-repo tie-aware `roc_auc` behind a one-line `p14_rc_auc` wrapper copied in shape
from `p13_gate_auc`.

**The degenerate denominator is guarded before the division.** `skill` normalizes by
`AURC_random − AURC_oracle`; if the two references coincided (every item correct, or every item
wrong, so that ordering the batch cannot change anything) the honest answer is `NaN` with a note,
never a number obtained by dividing by ~0. It did not trigger here (denominator 0.0751), and the
guard plus its note are persisted regardless.

## The anti-snooping record

The most snoopable number in this plan is the **coverage floor**: a floor chosen *after* seeing
where the curve becomes monotone is invisible in the resulting statistic and total in the claim.
Three things are on disk against that:

- `coverage_floor_frozen_before_run = true` and `coverage_floor_freeze_commit =
  "a6c825867dbc786f7c3927df6760295ee0930c77"` — a commit that contains no Phase-14 result of any
  kind, checkable by any reader with `git show`.
- `coverage_floor_rationale` — §C.2 item 1: below roughly this coverage the selective-risk
  denominator is too small for the ratio to mean anything. **That reason is real; it does not
  derive the specific value**, and the artifact says exactly that.
- `forbidden_action` — names the two moves that are not available here: monotonizing the curve
  (running-max or isotonic fit) and re-choosing the floor after seeing the curve.

`P14_ITERATION_ALLOWANCE = 1` could not have been spent on any of these three gates in any case:
its single pre-declared trigger is split-conformal coverage falling **below** the SC1-d band, and
14-09 measured 0.9115 against a band of 0.8799, so the trigger did not fire. Both facts are
persisted (`iteration_trigger_fired = false`, `iteration_allowance_applies_here = false`).

## The honesty block, written INTO the artifact

Seven string keys, each a **required** key of the atomic save — an artifact missing one is never
written at all, so the caveats cannot be dropped by a report writer (T-14-32 / T-14-33):

| Key | What it carries |
|---|---|
| `bars_note` | the four floors have no derivation; the four recorded cases of an underived bar measuring something other than what it named; the licensing standard for any amendment |
| `coverage_floor_rationale` | why a floor exists, and that the reason does not derive the value |
| `forbidden_action` | monotonizing the curve, and post-hoc floor selection, are breaches |
| `assumption_a1_note` | a two-line derivation that the population curve is non-decreasing **iff** the score is monotone in conditional risk — a property of the score, not of the estimator, and it can fail on a real head |
| `hard_note` | why SC3-c is not circular |
| `circularity_note` | **every SC1/SC3 number is a well-specified-regime number**; train-joint equals eval-joint by construction; the bounding evidence D-03 intended is **absent, not merely loose** |
| `amendment` | `P14_AMENDMENT_NOTICE` — SC1 amended by D-02, SC2 by D-05 |

## What these numbers do NOT say

The curve, the conformal quantile and the FDR check are all computed on draws from the **same
simulator the evidence net was trained on**. Training joint equals evaluation joint by
construction. So SC3-a/b/c say that the confidence ordering, the arithmetic and the calibration are
consistent with each other in the well-specified regime. They say **nothing** about how a real
microscopy batch degrades this, and they inherit the simulator's misspecification in full. The
held-out real-data arm that would bound that degradation is absent, not loose. A skill of 0.90 on
simulator draws is not a claim about a bench.

## Deviations from Plan

**None — the plan executed as written.** Two implementation choices were made inside the plan's
latitude and are recorded here because they affect how a number should be read:

1. **Reference-curve grid.** The plan mandates the descending-unique-κ sweep for the *measured*
   curve. The oracle and random references are **orderings**, not scores — the oracle's natural
   score takes two values, and a two-point curve would under-resolve the very envelope it exists to
   draw — so they are evaluated at every coverage `k/n`. Because the three AURCs must be integrals
   over a common interval for their difference to be a skill, the runner **asserts** span agreement
   (measured 0.0005 vs reference 0.0005) rather than assuming it, and persists the spans and the
   tolerance.
2. **`AURC` is a raw integral, not a span-normalized mean risk.** Any common normalization cancels
   exactly in the skill ratio; the span is persisted so a reader can convert without re-running
   anything.

### Assumption Drift (advisory)

None material. The plan's expectations (a step-function curve, a noisy low-coverage tail, and a
Spearman well clear of the floor on the declared range) matched what the run produced.

## Known Stubs

None. Every quantity the artifact names was computed and persisted.

## Threat Flags

None. This plan adds no network endpoint, auth path, file-access pattern or schema at a trust
boundary; it reads one existing `.jld2` and writes one `.jld2` plus one `.png`, both gitignored
regenerable artifacts.

## Self-Check

- `spike/p14/run_p14_riskcoverage.jl` — FOUND (1008 lines, committed)
- `spike/p14/p14_riskcoverage_report.jld2` — FOUND (gitignored bulk artifact, regenerable)
- `spike/figures/p14_riskcoverage.png` — FOUND (gitignored bulk artifact, regenerable)
- commit `4854eed` — FOUND
- commit `a22713d` — FOUND

## Self-Check: PASSED
