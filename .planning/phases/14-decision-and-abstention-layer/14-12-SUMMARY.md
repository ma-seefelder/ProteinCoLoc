---
phase: 14-decision-and-abstention-layer
plan: 12
subsystem: decision-and-abstention-layer
status: complete
tags: [sc3-d, ood-arm, misspecification, abstention, judgement-call-bars, pre-registered-gate, gate-not-met, named-limit]
requires:
  - "14-09 (spike/p14/p14_conformal_report.jld2 — the q-hat = 0.4498001628978108 the hedge was calibrated at, LOADED not re-calibrated)"
  - "spike/p14/pools.jl (p14_draw_pool, p14_misspec_pool, p14_ood_reference, p14_assert_unstratified)"
  - "spike/p14/decide.jl (p14_load_bundle, p14_decide_one, p14_ood_input — the D-07 inherited tau, prior and provenance)"
  - "spike/validation/ood.jl (OOD_FAMILIES and OOD_GRID_LEVELS — the four in-repo positive controls, CONSUMED not re-authored)"
  - "spike/p14/consts.jl (P14_OOD_MARGIN_FLOOR, P14_N_OOD, P14_OOD_COUNTER — frozen before this runner existed)"
  - "spike/data/cache/p13/<pool> (~51 MB gitignored; the Phase-13 net's OWN pool the density null is fitted on)"
provides:
  - "SC3-d: the per-family abstain-rate margin between a strongly misspecified arm and a matched in-distribution arm, with a conservative minimum-over-families headline"
  - "spike/p14/p14_ood_arm_report.jld2 — five matched arms, per-arm trigger breakdowns, the per-family margin table, the empty-set x OOD-state cross-tab and the honesty block"
  - "the recorded finding that the wired density channel is BLIND to misspec_noise at the strongest rung (0 of 1000 items flagged)"
affects:
  - "the Phase-14 verdict: SC3-d is the first Phase-14 gate to be NOT MET; SC1-b, SC1-d and SC3-a/b/c were all met"
  - "any report quoting a Phase-14 OOD claim: the single-wired-channel limit is now a measured limit rather than a stated one"
tech-stack:
  added: []
  patterns: ["persist-before-assert", "atomic .jld2 save with reopen integrity check", "run-time decoupling proof", "locked-threshold banner printed mechanically from P14_JUDGEMENT_CALL_BARS", "honesty block as REQUIRED artifact keys", "matched-arm design with item-by-item equality assertion", "minimum-over-families rather than mean", "ood_fired and ood_not_checked never summed"]
key-files:
  created:
    - "spike/p14/run_p14_ood_arm.jl (930 lines)"
    - "spike/p14/p14_ood_arm_report.jld2 (gitignored, reproducible; 738 KB)"
  modified: []
decisions: [D-05, D-06, D-07]
validates: ["SC3-d", "SC2-b"]
metrics:
  duration: "~2 h 40 min wall (the reported run itself: 3947.5 s = 65.8 min; one earlier attempt was killed at ~55 min before it could persist)"
  completed: 2026-08-04
---

# Phase 14 Plan 12: SC3-d — Does Abstention Concentrate on OOD Data? Summary

**SC3-d is NOT MET: the gated minimum-over-families abstain-rate margin is −0.129 against the frozen
P14_OOD_MARGIN_FLOOR = 0.50, and the shortfall is carried entirely by one family — the wired density
channel flags 0 of 1000 `misspec_noise` items at the strongest grid rung, while flagging 999/1000
`texture` and 989/1000 `background`.**

## What was measured

Five matched arms of n = 1000 each on the reserved `P14_OOD_COUNTER` stream: one in-distribution and
one per in-repo `OOD_FAMILIES` member at `level = :strongest` (rung 4 of 4 of the existing grid). No
misspecification was authored here. Every item of every arm went through the same composed
`p14_decide_one` path, at the same LOADED 14-09 q-hat (0.4498001628978108), against the same density
null (threshold 167.5446329378043), with `allow_unchecked_ood = false` (D-06).

### Per-arm abstain rates and trigger breakdown

| Arm | abstain rate | n_abstain | ood_fired | ood_not_checked | conformal_empty | conformal_ambiguous | disagreement_and_ood |
|---|---|---|---|---|---|---|---|
| in_distribution | 0.129 | 129 | 112 | 0 | 17 | 0 | 0 |
| texture | 0.999 | 999 | 999 | 0 | 0 | 0 | 0 |
| noise | 0.000 | 0 | 0 | 0 | 0 | 0 | 0 |
| optics | 0.387 | 387 | 376 | 0 | 11 | 0 | 0 |
| background | 0.989 | 989 | 989 | 0 | 0 | 0 | 0 |

(The in-distribution `ood_fired` count is 112 by subtraction from the persisted breakdown; the
persisted table is authoritative. Realized in-distribution class masses: exclusion 0.322 / random
0.474 / coloc 0.204, which passed `p14_assert_unstratified` against the measured prior.)

### Per-family margins

| Family | margin = rate(family) − rate(ID) | vs floor 0.50 | dominant channel | driven by NOT-CHECKED |
|---|---|---|---|---|
| texture | **+0.870** | MET | ood_fired | false |
| noise | **−0.129** | NOT MET | (no channel moved; see below) | false in substance |
| optics | **+0.258** | NOT MET | ood_fired | false |
| background | **+0.860** | MET | ood_fired | false |

**Gated headline: `margin = min over families = −0.129` (weakest family `noise`) < 0.50 → SC3-d NOT
MET.** The minimum, not the mean, is gated by design (T-14-41): the mean of these four is +0.465,
which would still have missed the floor, but quoting it would have hidden that one family produces
*zero* detections.

## Which channel drove the margin — the D-06 question, answered

**No item in any arm resolved `:not_checked`.** `any_not_checked = false` is persisted. Every
OOD-triggered abstention above is a *fired* detector, not a default abstention, so nothing here is
the failure mode T-14-39 exists to catch. The `margin_driver` field reports `ood_fired` as dominant
for texture, optics and background.

For `noise` the reported `dominant` symbol is `:ood_not_checked` — this is an **argmax artefact and
not a finding**, and the honest reading is stated here rather than left to the field: all five deltas
for that family are ≤ 0 (`ood_fired = −0.112`, `ood_not_checked = 0.0`, `conformal_empty = −0.017`,
the rest 0.0), so `argmax` selects the only non-negative one. The substantive fact is that the
`noise` arm abstained on *nothing at all*: 0 fired, 0 empty sets, 1000 singletons, 1000 `:clear`.

## The empty-set × OOD-state cross-tab — two independent signals

Rows `(:singleton, :ambiguous, :empty)` × columns `(:fired, :clear, :not_checked)`. `:ambiguous` is
empty in every arm, which is 14-09's LAC finding carried forward (this q̂ emits no 2-class sets).

| Arm | empty AND fired (agreement) | empty but CLEAR (disagreement) | fired but NOT empty (disagreement) | n_empty | n_fired |
|---|---|---|---|---|---|
| in_distribution | 1 | 17 | 111 | 18 | 112 |
| texture | 19 | 0 | 980 | 19 | 999 |
| noise | 0 | 0 | 0 | 0 | 0 |
| optics | 14 | 11 | 362 | 25 | 376 |
| background | 23 | 0 | 966 | 23 | 989 |

**Both readings, neither editorialised into the other.** The two channels share no assumption — the
density null needs the simulator to be right about densities, the empty conformal set needs only
exchangeability.

- **Agreement (corroboration):** in the three families the density channel detects, *every* empty
  set is on an item the null also flagged (texture 19/19, background 23/23, optics 14/25). Two
  channels that share no assumption pointing at the same items is real corroboration.
- **Disagreement, and it is the more interesting direction:** in-distribution, **17 of 18 empty sets
  land on items the density null called CLEAR.** The exchangeability channel is seeing something the
  density channel is not, on data the density channel considers ordinary. Optics shows the same at
  11/25.
- **Disagreement, the other way:** `fired but not empty` is large everywhere the detector fires — the
  density null flags hundreds of items whose conformal sets are perfectly ordinary singletons. The
  two signals are genuinely not redundant in either direction.
- **The `noise` arm produces neither signal.** Not one empty set and not one flag: both wired
  channels are silent on it.

## The finding, stated plainly

The claim SC3-d makes is *"abstention concentrates on data the model was not trained on."* Measured
against four deliberately injected misspecifications, that is **true for three of them and false for
the fourth**, and the fourth is the one named `noise`.

That is not a coincidence, and the artifact records why rather than inferring it after the fact:
**this lane wires only the summary-density channel** (`channels_wired = (:density,)`,
`channels_not_wired = (:pp, :noise)`). The unwired channels include a *noise* channel. So the one
family that slips past is the family the missing channel is for. The single-channel limit was
declared in `pools.jl` before this runner existed and printed in the banner before any arm was
built; SC3-d has now converted it from a stated limit into a **measured** one with a number attached.

Two things this does **not** license, both named in the artifact's `bar_note` and in the assertion
message:

1. **Wiring the noise channel now to clear the bar.** That would be authoring the experiment after
   seeing its result. `P14_ITERATION_ALLOWANCE = 1` has exactly one pre-declared trigger — the
   split-conformal coverage band, which 14-09 measured *above* the band — so the allowance is unspent
   and is not available here.
2. **Relaxing P14_OOD_MARGIN_FLOOR, dropping the `noise` family, or lowering the grid rung.**
   `spike/p14/consts.jl` is byte-unchanged against HEAD and was verified as such after the run.

A separate honest note on the bar itself: **P14_OOD_MARGIN_FLOOR = 0.50 has no derivation.** It is
the fifth member of `P14_JUDGEMENT_CALL_BARS` and is printed and persisted from that tuple
mechanically. Even the best-performing families would be scored against a number nobody derived.
That is a reason to report the per-family table as the primary result and the pass/fail as
secondary — which is what the artifact does — not a reason to move the number.

## Tasks and commits

| Task | Name | Commit | Files |
|---|---|---|---|
| 1 | The matched in-distribution vs misspecified arms (SC3-d) | `df01d6a` | `spike/p14/run_p14_ood_arm.jl` |
| 2 | The empty-set cross-tab, the honesty block, persist then assert | `53f18c6` | `spike/p14/run_p14_ood_arm.jl` |

## Verification performed (observed, not assumed)

- **Task 2 `<verify>` command run verbatim, exit 0:** printed `margin=-0.129 floor=0.5`, and every
  one of the thirteen required keys is present with `bar_is_judgement_call == true`.
- `length(d["arms"]) == 5`; all five arms carry `n == P14_N_OOD == 1000`; every
  `trigger_breakdown`'s key set is a subset of the closed `p14_abstain_reasons()`.
- `length(d["margin_per_family"]) == 4`.
- `d["channels_not_wired"] == (:pp, :noise)`.
- Every arm's `empty_set_crosstab` sums to exactly 1000.
- `grep -c 'function misspec\|OOD_FAMILIES = ' spike/p14/run_p14_ood_arm.jl` → **0**.
- `grep -v '^\s*#' spike/p14/run_p14_ood_arm.jl | grep -cE '\b0\.50\b|\b1000\b'` → **0**.
- **All ten `spike/test/test_p14_*.jl` files run individually, all PASS** (conformal, consts, decide,
  decoupling, fdr, fuse, pools, posterior, provenance, result).
- `git diff --quiet HEAD -- spike/p14/consts.jl` → clean. `git status --porcelain -- src
  spike/Project.toml spike/Manifest.toml corpus artifacts spike/data/cache/p13` → **empty**.
- **PERSIST-BEFORE-ASSERT confirmed empirically:** the run exited non-zero on the SC3-d
  `AssertionError`, and `spike/p14/p14_ood_arm_report.jld2` (738 KB) exists on disk with the complete
  report, including `sc3d_met = false`. The evidence of the failure survived the failure.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 — Blocking] `OOD_FAMILIES` was not in scope under the plan's prescribed guard**

- **Found during:** Task 1, first smoke run (`UndefVarError: OOD_FAMILIES not defined in Main`).
- **Issue:** `spike/p14/pools.jl` guards its `spike/validation/ood.jl` include on `fit_ood_nulls`,
  but `src/amortized/ood.jl` — pulled in transitively by `decide.jl`, which every P14 runner loads
  first — *also* defines `fit_ood_nulls`. The guard therefore short-circuits and the validation file
  is never reached, so `OOD_FAMILIES` (which exists **only** in the validation file) is absent. My
  first attempt guarded on `OOD_GRID_LEVELS`, which has the same problem via
  `spike/validation/consts.jl`.
- **Fix:** guard on `OOD_FAMILIES` itself — the symbol this runner actually needs — and place the
  include **first**, before `decide.jl`. Ordering is load-bearing, not tidy: `spike/validation/ood.jl`
  carries spike-lane twins of `fit_ood_nulls`, `maha_score` and `roc_auc`, and in Julia the last
  definition wins. Plans 14-09/14-10/14-11 all fitted their density nulls with the **shipped**
  definitions (their guards named symbols `src/` had already defined). Including the validation file
  *after* `decide.jl` would silently swap the implementation the null is fitted with, moving the
  threshold and making this runner's abstain rates incomparable with 14-09's. Loading it first leaves
  the shipped definitions on top, exactly as in 14-09.
- **Files modified:** `spike/p14/run_p14_ood_arm.jl`. **Commit:** `df01d6a`.

**2. [Rule 2 — Missing critical functionality] The operating point is now PINNED to 14-09's**

- **Found during:** Task 1, while reasoning about deviation 1.
- **Issue:** the include-ordering argument above was an *argument*. A wrong ordering would move the
  fitted Mahalanobis threshold silently, because a slightly different threshold is still a perfectly
  plausible threshold, and the ID-arm abstain rate quoted beside 14-09's would no longer be
  comparable.
- **Fix:** added `"pool_provenance"` to the required 14-09 keys and an executable assertion that the
  threshold fitted here equals the one 14-09 recorded. It passed: **167.5446329378043 both sides.**
  The ordering argument is now proved at run time rather than argued.
- **Files modified:** `spike/p14/run_p14_ood_arm.jl`. **Commit:** `df01d6a`.

### Process deviations

**3. Task 1 was verified by a SMOKE run, not by the reported run.**

Task 1's `<verify>` names the full reported command. Running it would have cost ~66 minutes and
produced an artifact Task 2 immediately overwrites (Task 1's file cannot persist the Task-2 keys).
Task 1 was therefore verified with the runner's own designed non-reported path (`main(n_ood = 6)` →
`_smoke.jld2`, gate not asserted), which exercised every code path including the arm-matching
assertions, and the **reported** run was executed once, in Task 2, against the complete file. The
smoke artifact's schema was also validated against the full acceptance list before the reported run
was launched. Stated here so no reader has to infer it.

**4. The first reported run was killed by the harness at ~55 minutes, before it could persist.**

No artifact was produced and nothing was left half-written (the atomic `.tmp` + `mv` save means a
kill mid-write leaves a discardable `.tmp`, not a torn artifact). The run was relaunched detached and
completed in 3947.5 s. This is a harness/compute note, not a result note: the counter-based Philox
stream makes the run reproducible, so the relaunch is a re-execution, not a re-draw.

## Assumption Drift (advisory)

**1. "Only the density channel is wired" moved from a stated limit to the cause of the gate result.**

- **Found during:** Task 2, reading the reported numbers.
- **Planned:** the single-wired-channel limit was to be *recorded* in the artifact as a named limit
  (plan action item 1 and the `channels_not_wired` honesty key) — framed as scope disclosure.
- **Actual:** it is the *mechanism* of the SC3-d shortfall. The gated minimum is set by the `noise`
  family, and the unwired channels are `(:pp, :noise)`. The limit did not merely accompany the
  result; it produced it.
- **Why it matters:** a reader of `margin = −0.129` alone would conclude the abstention layer does
  not track misspecification. The correct reading is that the layer tracks three of four injected
  misspecifications strongly (margins 0.870 / 0.860 / 0.258) and is blind to the fourth *through the
  channel this lane did not wire.* Advisory only; it changes no threshold and gates nothing.

## Known Stubs

None. No placeholder values, no unwired data paths, no TODOs.

## Threat Flags

None. This plan adds no network endpoint, auth path, file-write surface or schema at a trust
boundary. The runner is read-only with respect to `src/`, the frozen spike environment, `corpus/`,
`artifacts/` and the gitignored `spike/data/cache/p13` pool, and all four were asserted clean after
the run.

## Self-Check: PASSED

- `spike/p14/run_p14_ood_arm.jl` — FOUND (930 lines).
- `spike/p14/p14_ood_arm_report.jld2` — FOUND (738,407 bytes), all thirteen required keys present.
- Commit `df01d6a` — FOUND. Commit `53f18c6` — FOUND.
- `spike/p14/consts.jl` — byte-unchanged against HEAD.
