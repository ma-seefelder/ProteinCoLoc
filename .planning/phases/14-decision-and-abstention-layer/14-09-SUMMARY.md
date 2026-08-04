---
phase: 14-decision-and-abstention-layer
plan: 09
subsystem: decision-and-abstention-layer
status: complete
tags: [conformal, coverage, sc1-d, calibration, shared-eval-pool, pre-registered-gate]
requires:
  - "14-07 (p14_decide_one, the composed per-pair path)"
  - "14-08 (p14_draw_pool, p14_assert_unstratified)"
  - "spike/p14/conformal.jl (p14_lac_score, p14_conformal_calibrate, p14_hedge_diagnostics)"
  - "spike/p14/consts.jl (every threshold, frozen before any Phase-14 result existed)"
provides:
  - "SC1-d: a measured realized split-conformal marginal coverage with a recorded verdict"
  - "spike/p14/p14_conformal_report.jld2 — q-hat, coverage, band, set-size distribution, provenance"
  - "spike/p14/p14_eval_pool.jld2 — the SHARED n_eval=2000 evaluation set consumed by the wave-6 runners"
affects:
  - "the SC1-b (realized FDP), SC3-a/b/c (risk-coverage) and SC3-d (OOD in-distribution) runners, which now score the SAME items at the SAME counter"
tech-stack:
  added: []
  patterns: ["persist-before-assert", "atomic .jld2 save with reopen integrity check", "run-time decoupling proof", "locked-threshold banner", "smoke-mode path redirect"]
key-files:
  created:
    - "spike/p14/p14_conformal_report.jld2 (gitignored, reproducible)"
    - "spike/p14/p14_eval_pool.jld2 (gitignored, reproducible)"
  modified:
    - ".gitignore"
decisions: [D-02, D-03, D-03a, D-07]
validates: ["SC1-d", "SC1-e"]
metrics:
  duration: "~35 min wall (the reported run itself: 627.8 s)"
  completed: 2026-08-04
---

# Phase 14 Plan 09: Split-Conformal Coverage (SC1-d) Summary

**SC1-d PASSES.** Realized split-conformal marginal coverage is **0.9115** (1823 of 2000) against
the pre-registered lower band **0.8799**, a gap of **+0.0316**, with `q̂ = 0.4498` and a
non-vacuous hedge. `iteration_trigger_fired = false`, so the single `P14_ITERATION_ALLOWANCE`
remains unspent.

## The Reported Numbers

| Quantity | Value |
|---|---|
| `q̂` (LAC, at α = 0.10) | 0.4498001628978108 |
| Realized coverage | **0.9115** (1823 / 2000) |
| Pre-registered band (lower) | **0.8798753882025019** (0.880 to 3 dp) |
| Gap (coverage − band) | +0.03162461179749809 |
| Nominal 1 − α | 0.90 |
| `iteration_trigger_fired` | **false** |
| `vacuous_hedge` | false |
| Calibration scores (min / median / max) | 9.873e-8 / 0.01706 / 0.99562 |
| Elapsed | 627.8 s |

The band is **derived, not chosen**: `1 − α − 3·√(α(1−α)/n_eval)` recomputed by hand from the
stored `P14_ALPHA_CONFORMAL` and `P14_N_EVAL` matches the persisted `band_lower` to < 1e-12.

### Set-size distribution (reported beside the coverage number, Pitfall 9)

| Status | Count | Rate |
|---|---|---|
| singleton | 1954 | 0.977 |
| ambiguous | 0 | 0.000 |
| empty | 46 | 0.023 |

### Class masses — Pitfall 1 held (both draws unstratified)

| Draw | exclusion | random | coloc |
|---|---|---|---|
| calibration (n = 2000) | 0.3145 | 0.4650 | 0.2205 |
| evaluation (n = 2000) | 0.3315 | 0.4545 | 0.2140 |
| expected prior π | 0.31272 | 0.47073 | 0.21655 |

Neither draw is equal thirds: the `random` and `coloc` masses sit **11-12 binomial SE** from 1/3
(calibration 12.49 and 10.70 SE; evaluation 11.49 and 11.32 SE). Both draws tracked the measured
prior, and `p14_assert_unstratified` passed on each.

### Conformal status × OOD state

| | fired | clear | not_checked |
|---|---|---|---|
| singleton | 178 | 1776 | 0 |
| ambiguous | 0 | 0 | 0 |
| empty | 5 | 41 | 0 |

The OOD reference was **available** (`thr = 167.54`, `:density` wired; `:pp` and `:noise` not
wired), so no item resolved to `:not_checked` and the D-06 all-abstain outcome did not arise.

### The shared evaluation pool

`spike/p14/p14_eval_pool.jld2` holds all 2000 items. Decisions: `:not_coloc` 1381, `:coloc` 395,
`:abstain` 224. Abstention reasons reconcile exactly with the cross-tab: `:ood_fired` 183 and
`:conformal_empty` 41 — the 5 items that were *both* empty and OOD-fired are attributed to OOD
under the ABSTAIN-FIRST ordering (46 empty total − 5 OOD-attributed = 41).

## Honest Reading of the Result

Three things qualify the pass and belong in any report that quotes it:

1. **The hedge never actually hedges.** `ambiguous_rate` is **exactly 0.0** — at this `q̂` the LAC
   score never admits a 2-class set. The hedge's only departure from a bare point prediction is
   the empty set (2.3%). The run therefore sits squarely in 14-RESEARCH §B.3's *singleton-dominant*
   regime. `vacuous_hedge` is `false` only because the label is
   `(ambiguous_rate + empty_rate) < P14_AMBIGUOUS_RATE_FLOOR` and the empty sets carry it past the
   0.01 floor; on the ambiguous channel alone the hedge is silent. Coverage here is being delivered
   by classifier accuracy, not by set-valued caution.

2. **Coverage exceeds nominal (0.9115 > 0.90).** This is expected, not anomalous: split conformal's
   finite-sample guarantee is coverage ≥ 1 − α, and the `ceil((n+1)(1−α))` order statistic makes it
   mildly conservative.

3. **The guarantee is simulator-derived.** Persisted as
   `guarantee_basis = :simulator_derived_held_out_draws`. Calibration and evaluation draws come from
   the same forward simulator the evidence net was trained on, at reserved counters disjoint from
   the training stream — so train-joint equals eval-joint *by construction* and 0.9115 is a
   **well-specified-regime** number inheriting the simulator's misspecification in full. The
   real-data bound D-03 intended is **ABSENT, not merely loose** (D-03a). The full `circularity_note`
   is persisted inside the artifact so a reader who never opens a planning document still meets it.

`conformal_library` is persisted as `"none -- hand-rolled; SC1 AMENDED by D-02"`.

## Tasks

**Task 1 — the runner skeleton: ALREADY COMPLETE, not re-executed and not re-committed.** Delivered
by the pre-existing commit `b288a36`, which contained both the skeleton and a fully implemented
`main`. Its acceptance criteria were re-verified on disk rather than assumed: the inlined-threshold
scan returns 0, `_p14_save_report`/`jldsave` appears 6 times, the `git diff --quiet HEAD -- src`
guard appears twice, `abspath(PROGRAM_FILE)` appears exactly once, and the include-is-a-no-op
command exits 0.

**Task 2 — the reported run: executed.** No code fix was needed; the `main` committed in `b288a36`
already satisfied the Task 2 spec in full (required artifact keys, persist-before-assert ordering,
the shared eval pool). The only commit this plan added is the `.gitignore` deviation below.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] The reported `.jld2` outputs were not actually gitignored**

- **Found during:** Task 2, before the run
- **Issue:** `git check-ignore` returned exit 1 for `spike/p14/p14_conformal_report.jld2` and
  `spike/p14/p14_eval_pool.jld2`. They were being treated as gitignored bulk data but no rule
  covered them, so the run would have left ~700 KB of untracked generated output in the tree.
- **Fix:** Added `spike/p14/*.jld2` to `.gitignore`, following the established Phase-5 (05-04)
  reported-artifact pattern at `.gitignore:434-440` — the runner is the deliverable, the report is a
  reproducible output whose numbers are recorded verbatim in this SUMMARY. The glob also covers the
  `_smoke` redirect targets.
- **Files modified:** `.gitignore`
- **Commit:** `c2da70b`

### Not a deviation, recorded for the next executor

A smoke run at `n = 120` tripped `p14_assert_unstratified` with "THE DRAW IS CLASS-BALANCED". This
is a **small-n power artifact of the guard, not a defect**: three binomial SE around 1/3 is 0.129 at
n = 120, so even a genuinely π-distributed draw lands inside the "balanced" envelope. At the reported
n = 2000 the SE is 0.0105 and `random = 0.471` sits ~13 SE from 1/3, so the guard regains full power
and passed. **Nothing was changed** — smoke below roughly n = 400 simply cannot exercise this guard.

## Verification Performed

All commands were run and their real output observed.

- `julia --project=spike -t auto spike/p14/run_p14_conformal.jl` — **exit 0**, 627.8 s.
- Both artifacts exist with mtime 11:10, newer than the 10:59 run start.
- All 10 required report keys present (`qhat`, `coverage`, `band_lower`, `cal_masses`,
  `eval_masses`, `vacuous_hedge`, `guarantee_basis`, `circularity_note`, `amendment`,
  `iteration_trigger_fired`) — "keys ok".
- All 9 required pool keys present, `n = 2000`.
- `band_lower` recomputed by hand matches to < 1e-12 and equals 0.880 to three decimals.
- `cal_masses` / `eval_masses` are not approximately equal thirds (11-12 SE away).
- `git status --porcelain -- src spike/Project.toml spike/Manifest.toml corpus artifacts` — **empty**.
- `git diff HEAD -- spike/p14/consts.jl` — **empty**; no constant was modified.
- All 10 `spike/test/test_p14_*.jl` files run **individually** — 10/10 PASS (never via `runtests.jl`,
  which exits 1 early at the Phase-4 speedup gate).

## Known Stubs

None.

## Self-Check: PASSED

- `spike/p14/run_p14_conformal.jl` — FOUND
- `spike/p14/p14_conformal_report.jld2` — FOUND
- `spike/p14/p14_eval_pool.jld2` — FOUND
- `.planning/phases/14-decision-and-abstention-layer/14-09-SUMMARY.md` — FOUND
- commit `b288a36` (Task 1, pre-existing) — FOUND
- commit `c2da70b` (.gitignore) — FOUND
