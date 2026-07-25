---
phase: 11-registration-and-chromatic-uncertainty-as-latent
plan: 06
subsystem: validation
tags: [pre-registration, decision-record, d-06, d-04, abort-criterion, jld2, re-derivation, no-compute]

# Dependency graph
requires:
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 01
    provides: "the Tier-1 abort criterion (P11_PROBE_S_FLOOR = 0.9, P11_PROBE_SPAN_FLOOR = 0.02) and P11_ITERATION_ALLOWANCE = 1, frozen at d336699 before the probe existed"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 05
    provides: "spike/validation/p11_probe_report.jld2 and the Tier-2 measurements P11_PROBE_S_MEASURED / P11_PROBE_SPAN_MEASURED, appended at 2780a06"
provides:
  - ".planning/phases/11-.../11-PROBE-VERDICT.md — the mechanically evaluated D-06 abort criterion, the three pre-registered branches with their downstream consequences, and the recorded decision"
  - "VERDICT: ABOVE RESOLUTION — the abort criterion did NOT fire on either leg"
  - "DECISION: proceed — plans 11-07 through 11-11 authorised as planned"
  - "P11_ITERATION_ALLOWANCE = 1 carried forward UNSPENT (iterations spent: 0 of 1)"
affects: [11-07-research-trainer, 11-08-ladder, 11-09-breakdown, 11-10-real-image, 11-11-report]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "A pre-registered decision is re-derived from the raw artifact before it is recorded, so the verdict rests on recomputation rather than on transcription from a prior SUMMARY"
    - "An autonomously resolved blocking checkpoint records WHY it was resolvable without a human (the rule was locked and mechanical, and the alternative branches' preconditions were false) rather than merely that it was"

key-files:
  created:
    - .planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-PROBE-VERDICT.md
  modified: []

key-decisions:
  - "DECISION: proceed — neither abort leg fired (S_probe 1.0 vs floor 0.9; span 0.046885 vs floor 0.02), so branches 2 and 3, whose entry condition is 'the criterion fired', were not available to be chosen"
  - "The iteration allowance was deliberately left UNSPENT: spending it on a passing criterion is irreversible, leaving it is not"
  - "The blocking checkpoint was resolved autonomously, but made conditional on an independent re-derivation from the raw artifact tables — a disagreement would have been recorded as BLOCKED instead"

requirements-completed: [D-04, D-06]

# Metrics
duration: ~20min
completed: 2026-07-25
---

# Phase 11 Plan 06: The D-06 Probe Verdict Summary

**The pre-registered abort criterion was evaluated from the frozen constants and independently
re-derived from the raw probe artifact, and it did NOT fire on either leg — `S_probe = 1.0` against
a floor of 0.9 and a ladder span of 0.046885 Δρ_eq against a floor of 0.02 — so the recorded
decision is PROCEED, plans 11-07 through 11-11 are authorised as planned, and
`P11_ITERATION_ALLOWANCE = 1` is carried forward UNSPENT with no threshold altered.**

## Verdict: **ABOVE RESOLUTION** — abort criterion NOT FIRED

| Quantity | Measured (F5 mixture arm) | Pre-registered floor | Leg fires? |
|---|---|---|---|
| `S_probe` | **1.0** | needs >= 0.9 (`P11_PROBE_S_FLOOR`) | no |
| ladder span (Δρ_eq) | **0.046884564903982** | needs >= 0.02 (`P11_PROBE_SPAN_FLOOR`) | no |

## DECISION: proceed

**Downstream authorisation:** Plans 11-07 through 11-11 are authorised as planned.
**Iterations spent: 0 of 1.**

Two pre-registered conditions attach and are recorded in the verdict file:

- `P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150` is a **BLOCKER** threshold, not a downgrade trigger.
  Exceeding it must be recorded as a blocker; it is explicitly not licence to downgrade the arm to
  256².
- The SC1g λ-ablation tripwire now bars at `P11_LAMBDA_ABLATION_FACTOR = 2.502` (up from the `1.15`
  placeholder). A failure there in plan 11-07 is a **result**, not licence to relax the constant.

## Performance

- **Duration:** ~20 min wall
- **Tasks:** 3 (Task 2 was the blocking decision checkpoint, resolved autonomously — see below)
- **Files created:** 1; **modified:** 0; **deleted:** 0
- **Compute spent:** none — no simulation, no training, no seed consumed

## Task Commits

1. **Task 1: evaluate the criterion and draft the verdict** — `4f7daa8` (docs)
2. **Task 2: the blocking decision checkpoint** — no commit; resolved autonomously (below)
3. **Task 3: record the decision, the ledger and the authorisation** — `6faafd4` (docs)

## The independent re-derivation (the substantive work of this plan)

The plan asks for the criterion to be *computed* from the frozen file rather than transcribed. I
went one step further as instructed and recomputed **both statistics from the raw
per-`(θ, replicate, rung)` tables** inside `spike/validation/p11_probe_report.jld2` — `shift_norm_f5`
(dims `(5, 32, 12)`) and `drho_norm_f5` (dims `(5, 32, 4)`) — using per-rung means, an OLS fit, the
Δρ_eq inversion, a midrank Spearman and a span, all written for this verdict rather than by calling
the probe script's own helpers. No RNG, no simulation, no stream consumed.

| Quantity | Re-derived from raw tables | Tier-2 frozen | Artifact headline field |
|---|---|---|---|
| `S_probe` | 1.0 | 1.0 | 1.0 |
| ladder span | 0.046884564903982295 | 0.046884564903982365 | 0.046884564903982365 |
| Δρ_eq slope (F5) | 6.4585054837724964 | 6.458505483772493 | — |
| Δρ_eq intercept (F5) | 0.018126542853519334 | 0.018126542853519556 | — |
| λ_max/λ_min ratio | 5.003943268109935 | 5.003943268109952 | — |

Residuals are at the `1e-16` level — float summation-order artifacts of a differently written
reduction, five orders of magnitude below the `1e-12` agreement tolerance and incapable of moving
either leg. The script printed `REDERIVATION_AGREES` and
`abort_criterion_fired (re-derived) = false`, matching `abort_criterion_fired (frozen) = false`.

**This mattered procedurally, not just as diligence:** the autonomous resolution below was made
*conditional* on this agreement. Had it disagreed, the verdict would have been written as BLOCKED.

## The blocking checkpoint, and how it was resolved

Task 2 is a `checkpoint:decision` with `gate="blocking"`. It was resolved **autonomously by the
orchestrator under the pre-registered rule, with no user present**, because the phase ran in a
background non-interactive session.

The justification recorded in the verdict file is narrow and specific: the D-06 criterion is
pre-registered and mechanical, so the verdict is a deterministic evaluation of an already-locked
rule rather than a fresh judgement. The floors were frozen in Tier 1 at `d336699` *before the probe
existed*; the measurements were frozen in Tier 2 at `2780a06`. Both legs clear their floor, so
**branches 2 and 3 — whose entry condition is "the criterion fired" — were not available to be
chosen.** Branch 1 is the only branch whose stated precondition holds.

The one genuinely discretionary sliver — whether to spend the allowance *despite* a passing
criterion, e.g. to buy a wider ladder — was **not** exercised on the user's behalf. It was left
unspent, which is the reversible option: the allowance survives for plan 11-07 or 11-08, whereas
spending it here could not have been undone.

## Verification — observed results

Every command below was run and its real output observed.

| Check | Command | Observed |
|---|---|---|
| Task 1 verify (criterion from frozen consts) | `julia --project=spike -e 'include(...p11_consts.jl); fired = ...; println("abort_criterion_fired=", fired)'` | `abort_criterion_fired=false` |
| Plan verification 1 (same, bare boolean form) | `julia --project=spike -e 'include(...); println((...) \|\| (...))'` | `false` |
| Independent re-derivation from the raw artifact | ad-hoc script over `p11_probe_report.jld2` | `S_probe=1.0`, `span=0.046884564903982295`, `RESULT: REDERIVATION_AGREES`, re-derived verdict `false` |
| Pre-registration intact (no threshold moved) | `julia --project=spike spike/test/test_p11_consts.jl` | **84 Tier-1 + 32 Tier-2 = 116 pass / 0 fail**, incl. the testset "the TIER-1 block is still exactly what plan 11-01 locked" (11/11) |
| Task 1 acceptance: allowance named | `grep -c 'P11_ITERATION_ALLOWANCE' 11-PROBE-VERDICT.md` | `2` (bar: >= 1) |
| Task 1 acceptance: three branches present | `grep -c 'PROCEED\|SPEND THE ITERATION\|STOP SC2'` | `6` (bar: >= 3) |
| Task 1 acceptance: verdict leads the doc | first `^## ` heading | `## Verdict: **ABOVE RESOLUTION**` (line 13) |
| Task 1 acceptance: scope annotations | `grep -n 'interpolation onset\|registration ladder only'` | both present (lines 83, 182) |
| Task 3 verify: placeholder gone | `grep -c 'DECISION: pending'` | `0` |
| Task 3 acceptance: exactly one decision line | `grep -cE 'DECISION: (proceed\|spend-iteration\|stop-sc2)'` | `1` (line 256, `DECISION: proceed`) |
| Task 3 acceptance: ledger as integers | `grep -c 'of 1'` | `3` (bar: >= 1) |
| No spike/src byte moved | `git diff --stat 9e99b69..HEAD -- spike/ src/` | empty |
| No deletions in either commit | `git diff --diff-filter=D --name-only HEAD~1 HEAD` per commit | empty each time |

**NOT verified here, and deliberately not claimed.** This plan trains nothing, evaluates no net,
and reads no posterior.

- `S_probe = 1.0` is a property of the **forward model's summary displacement**, not of any
  posterior. Whether posterior width tracks λ is unknown and is measured for the first time in
  plans 11-07/11-08. Passing the probe buys the right to try; it does not predict SC2.
- No calibration, coverage, monotonicity or identifiability claim is made or supported.
- The shift and `chromatic_eps` θ columns remain **vacuous by design** (D-05).
- `spike/test/runtests.jl` was **not** run to green — the pre-existing `SPEEDUP_GATE` blocker
  stands and was **not** lowered (see Deferred Issues).

## STATE.md — proposed text for the orchestrator (NOT written by this executor)

Per the parallel-execution contract I did not touch `.planning/STATE.md`. The plan's Task 3 asks
for a dated Current Position line and, on a non-`proceed` branch, a Blockers entry. The branch is
`proceed`, so **no Blockers entry is required**. Proposed line for `## Current Position`:

> - **2026-07-25 — Phase 11 probe verdict (plan 11-06): ABOVE RESOLUTION, DECISION `proceed`.**
>   The D-06 abort criterion did not fire: `S_probe = 1.0` against `P11_PROBE_S_FLOOR = 0.9`, and a
>   Δρ_eq ladder span of 0.046885 against `P11_PROBE_SPAN_FLOOR = 0.02` (F5 mixture arm; both
>   independently re-derived from `p11_probe_report.jld2`). Plans 11-07 through 11-11 are
>   authorised as planned. `P11_ITERATION_ALLOWANCE = 1` is UNSPENT (0 of 1); no threshold was
>   altered. The ~1.5 h datagen-plus-training spend is authorised, and
>   `P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150` is a BLOCKER threshold, not a downgrade trigger.
>   Checkpoint resolved autonomously under the pre-registered rule (no user present), conditional
>   on the re-derivation agreeing. See `11-PROBE-VERDICT.md`.

## Decisions Made

1. **PROCEED, and the reason is that the alternatives were unavailable rather than unattractive.**
   Branches 2 and 3 both open with "the criterion fired". It did not, on either leg, so branch 1 is
   the only branch whose precondition holds. Recording it this way keeps the abort branch
   falsifiable: it did not fire *this time*, which is a different statement from "we chose not to
   fire it".
2. **The iteration allowance is left UNSPENT, deliberately.** Spending an allowance on a passing
   criterion is irreversible; leaving it is not. It stays available for 11-07/11-08.
3. **The autonomous resolution was made conditional, not assumed.** The verdict file records that a
   disagreement between the re-derivation and the frozen values would have produced a BLOCKED
   verdict, so the autonomy is bounded by a check rather than by assertion.
4. **The Task-1 draft was committed with `DECISION: pending` before the decision commit.** Two
   commits rather than one, so the git history shows the criterion was evaluated *before* a branch
   was selected — the same audit-trail logic the two-tier consts file uses.

## Deviations from Plan

### Adjustments required by the execution contract

**1. `.planning/STATE.md` not written (orchestrator-owned)**
- **Found during:** Task 3
- **Issue:** The plan lists `.planning/STATE.md` in `files_modified` and Task 3 directs an edit to
  it, but the parallel-execution contract reserves STATE.md to the orchestrator.
- **Resolution:** The file was left untouched and the exact proposed text is reproduced above under
  "STATE.md — proposed text for the orchestrator". Task 3's STATE.md acceptance criteria
  (`grep -c 'probe verdict' .planning/STATE.md` >= 1, and the Blockers-entry criterion) are
  therefore **not evaluated by this executor**; the Blockers criterion is moot in any case because
  it is conditioned on a non-`proceed` branch.
- **Files modified:** none

**2. The blocking checkpoint was resolved without a user**
- **Found during:** Task 2
- **Issue:** Task 2 is `checkpoint:decision`, `gate="blocking"`, and normal protocol is to STOP.
  The session is background and non-interactive, so no prompt could be answered.
- **Resolution:** Resolved under the pre-registered rule per the orchestrator's explicit
  instruction, and — as required — the resolution was verified rather than trusted by re-deriving
  both statistics from the frozen artifact. The autonomy, its justification, the unspent allowance
  and the unchanged thresholds are all recorded in `11-PROBE-VERDICT.md` under
  "How this checkpoint was resolved", so a later reader cannot mistake it for a human decision.
- **Files modified:** `11-PROBE-VERDICT.md`

### Auto-fixed Issues

**None.** No Rule 1/2/3 situation arose: this plan writes documents and evaluates frozen constants.
No Rule 4 (architectural) situation arose. No package was installed. No dependency, no seed, no
threshold and no byte of `spike/` or `src/` was touched.

## Assumption Drift (advisory)

**1. "The verdict document's numbers can be carried over from the 11-05 SUMMARY."**
- **Found during:** Task 1
- **Planned:** the plan directs reproducing the sweep tables "from the 11-05 SUMMARY", which frames
  the SUMMARY as the source of record for the numbers.
- **Actual:** the two *decision-bearing* numbers were instead recomputed from the raw artifact
  tables, and the SUMMARY was used only for the descriptive sweep tables. The re-derived slope and
  intercept differ from the frozen Tier-2 constants in the 16th significant figure.
- **Why it matters:** it changes what the verdict rests on. A transcribed verdict inherits any
  error in the transcription chain; a re-derived one does not. It also establishes that the
  artifact's raw tables — not just its headline fields — are sufficient to reconstruct the
  criterion, which is what makes the pre-registration auditable by a third party.

## Issues Encountered

- The worktree spawned at a commit **behind** the declared base (`0d9b87c`, an ancestor of
  `9e99b69`). The startup protocol's base assertion caught it and `git reset --hard 9e99b69`
  corrected it before any work; the tree was clean, so nothing was lost. Worth noting only because
  the plan's inputs (`p11_consts.jl` Tier 2, `p11_probe_report.jld2`) land in `2780a06`, which is
  reachable only from the corrected base — executing from the stale base would have read a
  Tier-2-less consts file and failed on `P11_PROBE_S_MEASURED` being undefined.

## Deferred Issues

- **`spike/test/runtests.jl` still exits 1** on the pre-existing `SPEEDUP_GATE` (NPE-03) blocker
  recorded in `deferred-items.md` (92.5× / 84.0× against a `> 100×` bar). Not this plan's to fix;
  the gate was **not** lowered. This plan added no file to the aggregate runner and ran no code
  outside `test_p11_consts.jl` and a read-only artifact re-derivation.

## Known Stubs

**None.** The verdict document is complete: the criterion is evaluated, the verdict leads, all
three branches are laid out with concrete downstream consequences, the ledger states spent and
remaining as integers, and the decision is recorded with its reason and authorisation.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no schema at a trust boundary, and no code.
The registered threats were handled as planned:

| Threat ID | Handling | Evidence |
|---|---|---|
| T-11-24 (unrecorded or ambiguous branch decision) | `11-PROBE-VERDICT.md` carries the branch id, the one-sentence reason, the date and the explicit downstream authorisation; the STATE.md mirror is drafted above for the orchestrator, and no Blockers entry is due on a `proceed` branch | `grep -cE 'DECISION: (proceed\|spend-iteration\|stop-sc2)'` = `1` |
| T-11-25 (relaxing a floor or re-seeding for a nicer verdict) | Both named as forbidden under every branch in the document; `p11_consts.jl` byte-unchanged and its literal-value gate re-run green | `git diff --stat 9e99b69..HEAD -- spike/` empty; `test_p11_consts.jl` 116/116 incl. the Tier-1-unchanged testset |
| T-11-26 (spending a second, unauthorised iteration) | The ledger states 0 spent / 1 remaining against `P11_ITERATION_ALLOWANCE = 1`, and records that a second cannot be authorised by amending the file | `grep -c 'of 1'` = `3` |
| T-11-SC (package-manager installs) | None attempted; no code and no packages | no `Project.toml`/`Manifest.toml` in the diff |

## User Setup Required

None.

## Next Phase Readiness

- **Ready — plan 11-07 (research trainer).** Training is authorised. Two pre-registered conditions
  travel with the authorisation: the 150-minute datagen ceiling is a BLOCKER (not licence to
  downgrade to 256²), and the λ-ablation tripwire bars at 2.502, which is materially harder than
  the retired 1.15 placeholder.
- **Ready — plan 11-08 (the SC2 ladder).** `SC2_SPEARMAN_FLOOR = 0.5` and the Δρ_eq coefficients are
  frozen; the ladder's Δρ_eq axis is the **F5 mixture** fit (`6.458505 · Δρ + 0.018127`), not the
  256² one.
- **For the report (plan 11-11).** Three report-facing items are restated in the verdict file: the
  interpolation-onset reading of the 0 → 0.25 px step, the negative-zero-row artifact of the affine
  inversion, and the fact that the two misalignment axes are **not** comparable in magnitude on the
  realistic arm (chromatic 0.1895 vs shift 0.0586 at their prior edges, reversing the 256²
  relation). Also that SC2 is a **registration** ladder only.
- **For the orchestrator.** Apply the proposed STATE.md Current Position line above. No Blockers
  entry is due.
- **Left red for someone else:** the pre-existing `SPEEDUP_GATE` blocker in `spike/test/runtests.jl`.

## Self-Check: PASSED

- Created file exists on disk:
  `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-PROBE-VERDICT.md`.
- Both task commits exist in `git log`: `4f7daa8`, `6faafd4`.
- `git diff --stat 9e99b69..HEAD -- spike/ src/` is empty — no pre-registration byte moved, no
  simulator or package code touched.
- `.planning/STATE.md` and `.planning/ROADMAP.md` are untouched, as the parallel-execution contract
  requires.

---
*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Completed: 2026-07-25*
