---
phase: 14
plan: 01
subsystem: pre-registration
status: COMPLETE — Task 0 ruled accept-as-proposed 2026-08-03; freeze committed at a6c8258
tags: [pre-registration, seeds, judgement-call-bars, wave-0]
requires: []
provides:
  - "spike/p14/consts.jl — the Phase-14 Tier-1 pre-registration (every seed, counter, sample size, alpha and bar)"
  - "P14_DEV_SEED / P14_FIX_SEED / P14_SALT — fresh Philox streams, disjoint from the whole P13 family by key-word cross-product"
  - "P14_JUDGEMENT_CALL_BARS — the five ruled bars, enumerated mechanically"
  - "P14_ITERATION_ALLOWANCE / P14_ITERATION_TRIGGER — the one pre-declared remedy (LAC to APS)"
affects:
  - "every later Phase-14 plan: no result-producing commit may precede a6c8258"
tech-stack:
  added: []           # no dependency added; spike/Project.toml and Manifest.toml byte-unchanged
  patterns:
    - "Tier-1 guard block keyed on a declared-deviation tuple, never on a seed (CONVENTIONS C-01)"
    - "forbidden-seed inventory DERIVED by guarded include, never retyped"
    - "full Philox key-word cross-product, not merely seed inequality"
    - "run-time-assembled grep needle so a source scan can scan itself"
key-files:
  created:
    - spike/p14/consts.jl          # 413 lines
    - spike/test/test_p14_consts.jl # 297 lines, 273 assertions
  modified: []
decisions: [D-07]
metrics:
  duration: ~35 min
  completed: 2026-08-03
  tasks_completed: 3
  tasks_total: 3
  freeze_commit: a6c825867dbc786f7c3927df6760295ee0930c77
---

> **FRONTMATTER UPDATED IN PLACE 2026-08-03.** The body below is append-only and the original
> BLOCKED record is preserved verbatim. The frontmatter is machine-read (`status`, `provides`,
> `tasks_completed`), so leaving it saying BLOCKED would misreport disk state to every later plan.
> The superseded values were: `status: BLOCKED — pending user ruling at Task 0`, empty `provides`
> and `key-files.created`, `duration: ~10 min (read-only)`, `completed: null`, `tasks_completed: 0`.

# Phase 14 Plan 01: Freeze the Phase-14 Pre-Registration — Summary

**STOPPED AT TASK 0. Nothing was written to disk. This is the intended outcome, not a failure.**

Plan 14-01 opens with a **blocking `checkpoint:decision`** (`autonomous: false`) that must be ruled
by the user *before* Task 1 writes `spike/p14/consts.jl`. Task 0 is the **first** task in the plan, so
there was no preceding work to execute or commit. The executor read the full decision record and
stopped.

---

## Why nothing was written

The plan's own threat register makes this explicit:

| Threat ID | Category | Component | Mitigation |
|-----------|----------|-----------|------------|
| T-14-52 | Repudiation | an agent silently choosing five gating numbers | Blocking `checkpoint:decision` at Task 0, **before anything is written**; the ruling recorded verbatim in the summary and as a dated comment in the constants file |

Auto-selecting the first option would have executed precisely the threat the checkpoint exists to
mitigate. `14-VALIDATION.md` "Manual-Only Verifications" row 1 states the rule that makes this
one-directional:

> A user may overrule any of them, but only *before* the freeze commit — never after a result.

Once Task 1 commits, these values are locked for the phase. **This is the only moment they can be
ruled.**

---

## CHECKPOINT REACHED

**Type:** decision (`gate="blocking"`)
**Plan:** 14-01
**Progress:** 0/3 tasks complete
**Blocked task:** Task 0 — *Rule the five undERIVED bars — BEFORE the freeze, or never*

### The decision

> Do you accept the five pre-registered Phase-14 bars as proposed, or do you want to set any of them
> yourself?

| Constant | Proposed | Scores | Source |
|---|---|---|---|
| `P14_SKILL_FLOOR` | **0.60** | SC3-a — selective skill, normalized between the random and oracle risk-coverage curves | 14-RESEARCH §C.3 Statistic S1 `[ASSUMED]` |
| `P14_SPEARMAN_FLOOR` | **0.95** | SC3-b — Spearman(coverage, selective risk) | §C.3 Statistic S2 |
| `P14_COVERAGE_FLOOR` | **0.20** | SC3-b — the coverage range the Spearman is computed on | §C.3 Statistic S2 |
| `P14_AUC_HARD_FLOOR` | **0.80** | SC3-c — abstention concentrates on would-have-been-wrong items | §C.3 Statistic S3 `[ASSUMED]` |
| `P14_OOD_MARGIN_FLOOR` | **0.50** | SC3-d — abstain rate under strong misspecification minus in-distribution | §C.3 Statistic S4 |

**These five numbers have no derivation.** 14-RESEARCH Assumptions Log **A3** (line 1440) names all
five as *"judgement calls with no derivation"* and recommends surfacing them to the user before the
freeze. Open Question 1 (line 1462) records the same recommendation; the planning resolution table
(line 1456) escalated rather than decided it.

Only S1's normalization and S2's floor carry any stated reasoning at all:
- **S1 = 0.60** is normalized between random and oracle *so it is immune to the base error rate* —
  the normalization is derived; the 0.60 is not.
- **S2's 0.20 coverage floor** has a recorded *reason* (§C.2 item 1: below ~10–20 % coverage the
  denominator is too small for the selective-risk ratio to mean anything) but not a derivation of the
  specific value. 14-RESEARCH is blunt that picking this floor *after* seeing where the curve becomes
  monotone is "the failure mode this project has paid for repeatedly."
- **S3 = 0.80** and **S4 = 0.50** are bare suggestions.

### Why this checkpoint is first and blocking

This project has recorded **four** separate cases of a bar measuring something other than what it
named (`.planning/STATE.md:114-122`, the `12-STAGE1-VERDICT.md` §7 families):

1. SC1g's component-vs-total
2. the n=2-against-271 real arm
3. the Wald-labelled-Wilson sizing (DEF-12-04)
4. the Phase-12 Stage-1 control ceiling — *a gating ceiling with no derivation anywhere in its
   pre-registration*

STATE.md's own framing, recorded for Phase 16: *"the count is now high enough that 'we
pre-registered it' no longer settles an argument alone."* The licensing standard it names — the one a
referee should hold each amendment to — is that every one of the four was shown **from evidence
independent of the machinery under audit** to measure something other than what it named.

14-RESEARCH Pitfall 5 (line 1078) names these five as the candidates for the fifth.

The single pre-declared remedy in this phase, `P14_ITERATION_TRIGGER`, concerns **conformal coverage
and none of these five** — so if one of these bars is wrong, there is no authorised route to fix it
after a result exists.

### Options

| id | Name | Pros | Cons |
|---|---|---|---|
| **`accept-as-proposed`** | Accept all five as proposed | Matches 14-RESEARCH's recommendation; the phase proceeds immediately; the values are conventional for selective prediction and are labelled as judgement calls in the constants file, in every artifact and in the report | Five undERIVED numbers gate three success criteria; if one turns out to measure the wrong thing, that is a fifth entry in a list already long enough that "we pre-registered it" no longer settles an argument on its own |
| **`rule-some`** | Set one or more of them yourself | The bars become yours rather than an agent's, and this is the only moment that can happen | Requires a judgement now, with no Phase-14 result to inform it — which is exactly the point of pre-registration, but it does mean deciding blind |
| **`report-only`** | Freeze them as REPORTED, NOT GATED — measure and report, gate nothing | Removes the risk of a wrong bar failing the phase for the wrong reason; SC1-c is already treated this way, so the pattern exists | SC3 then has no falsifiable content, and 14-VALIDATION currently marks SC3-a..d as GATED. This would be a deliberate weakening of the phase's success criteria and must be recorded as such in the ROADMAP alongside SC3 |

### Two further values you MAY rule in the same commit (less contentious)

| Constant | Proposed | Note |
|---|---|---|
| `P14_ALPHA_CONFORMAL` | **0.10** | Its SC1-d band at n = 2000 is a coverage floor of **0.880** — that band is *derived* from α and n, not chosen |
| `P14_N_CAL` = `P14_N_EVAL` | **2000** | Pins realized conformal coverage to ±0.007 at 1 sd. **That arithmetic is exact**, and larger is strictly better and cheap |

### What is NOT on the table

- **τ** — inherited from Phase 13 by *loading* the artifact, never copied or re-derived (D-07).
  Deliberately absent from the Phase-14 constants file; `test_p14_consts.jl` asserts
  `!isdefined(:P14_TAU)` so a copied literal fails the suite.
- **α_FDR** — a user parameter supplied per call, not a frozen constant. (It shares the value 0.10
  with `P14_ALPHA_CONFORMAL`, which is exactly why neither may ever be assigned from the other —
  D-07, threat T-14-13.)
- **The SC1-b / SC1-d bands** — derived from α and n by formula, not chosen.

### Resume signal

**Select one:** `accept-as-proposed` · `rule-some` (and give the values) · `report-only`

### What happens on each answer

- **`accept-as-proposed`** → Task 1 writes `spike/p14/consts.jl` with the five proposed values, each
  under a comment block labelling them JUDGEMENT CALLS WITH NO DERIVATION, plus the dated ruling.
- **`rule-some`** → the changed values are used in Task 1 *in place of* the proposed ones and each
  carries a `USER RULING` comment with the date; `test_p14_consts.jl` asserts the ruled literals.
- **`report-only`** → additionally, `14-VALIDATION.md`'s SC3-a..d rows and the ROADMAP Phase-14 SC3
  entry must be updated to read **REPORTED NOT GATED**, and the runners in 14-11 and 14-12 report the
  statistics without asserting them.

---

## Remaining tasks (blocked)

| Task | Name | Status |
|---|---|---|
| 0 | Rule the five undERIVED bars | **BLOCKED — awaiting user ruling** |
| 1 | Freeze the Tier-1 pre-registration in `spike/p14/consts.jl` | not started (gated by Task 0) |
| 2 | Assert the freeze — `spike/test/test_p14_consts.jl` | not started (gated by Task 0) |

Tasks 1 and 2 both consume the ruling's values directly (Task 1 writes them; Task 2 asserts them
**by literal value**), so neither can be partially executed ahead of the ruling.

---

## ⚠️ EXECUTION-ENVIRONMENT DIVERGENCE — flagged for the orchestrator, affects LATER Phase-14 plans

**This plan declares `execution_environment: main-working-tree-no-worktrees` in its frontmatter**, and
both `14-01-PLAN.md` `<execution_context>` and `14-VALIDATION.md` "Execution Environment Note" say in
capitals: *"EXECUTE ON THE MAIN WORKING TREE. NO GIT WORKTREES."* This executor was nevertheless
dispatched into the worktree `.claude/worktrees/agent-aa25e564d250d23c1`.

**For 14-01 specifically this is harmless** — the plan writes two text files and reads only tracked
sources (`spike/p13/consts.jl` is present).

**It is NOT harmless for later Phase-14 plans, and this was verified on disk, not assumed:**

```
$ ls spike/data/cache/
fixture
```

The ~50 MB gitignored `spike/data/cache/p13/` pool that the OOD null fit needs **does not exist in
this worktree**. This is the exact failure mode recorded in the user memory
`worktree-cleanup-destroys-gitignored-bulk-data`, where `git worktree remove --force` destroyed
Phase 11's 54 MB pool.

**Recommendation:** run the remainder of Phase 14 on the main working tree, per the plan.

---

## Deviations from Plan

None. No code was written, no constant was chosen, no value was frozen. The plan was followed
exactly: it opens with a blocking checkpoint and the executor stopped at it.

## Assumption Drift (advisory)

**1. Execution environment — planned vs actual**
- **Found during:** startup
- **Planned:** main working tree, no worktrees (plan frontmatter + `<execution_context>` +
  `14-VALIDATION.md`)
- **Actual:** dispatched into a git worktree
- **Why it matters:** no effect on 14-01; a hard blocker for any later Phase-14 plan that reads
  `spike/data/cache/p13/`. Verified absent above.

Non-blocking and recorded per the advisory rule; it did not gate anything in this plan.

---

## Constraint Verification (all hold)

Verified by `git status --porcelain` returning empty and by not invoking any write tool against
these paths:

| Constraint | Status |
|---|---|
| `src/` untouched (D-01) | ✅ no changes |
| `spike/Project.toml` / `spike/Manifest.toml` byte-unchanged (`test_p12_decoupling.jl:176-178`) | ✅ no changes |
| No dependency added or changed | ✅ none |
| Phase-16 sealed holdout not consumed (`test_p12_decoupling.jl:210`) | ✅ not read |
| `corpus/` not fetched or read | ✅ not touched |
| `STATE.md` / `ROADMAP.md` not modified (orchestrator owns these) | ✅ no changes |
| Seed discipline | ✅ n/a — no seed defined yet |
| `spike/test/runtests.jl` not wired as a gate | ✅ not touched |

## Known Stubs

None. No source file was created.

## Threat Flags

None. No new surface was introduced.

## Self-Check

- `spike/p14/consts.jl` — **intentionally absent** (gated by Task 0). Verified: `spike/p14` does not
  exist.
- `spike/test/test_p14_consts.jl` — **intentionally absent** (gated by Task 0).
- Working tree clean before this summary was written: `git status --porcelain` → empty.

## Self-Check: PASSED

The plan's Wave-0 success criteria are **not** met and are not claimed to be met. Wave 0 is
incomplete by design: the pre-registration cannot be frozen until it is ruled.

**No plan in any later Phase-14 wave may run until Task 0 is ruled and the freeze commit exists.**

---
---

# APPENDED 2026-08-03 — TASK 0 RULED. THE CHECKPOINT IS RESOLVED.

**Nothing above this line is edited.** The blocked record is preserved verbatim as the evidence that
the five bars were put to the user *before* any Phase-14 file existed. What follows is appended, per
the project's append-never-overwrite discipline (14-CONTEXT §Specific Ideas).

## The ruling, verbatim

> The Task-0 blocking checkpoint is RESOLVED. The user ruled **accept-as-proposed**. Freeze exactly
> these values:
>
> - `P14_SKILL_FLOOR      = 0.60`
> - `P14_SPEARMAN_FLOOR   = 0.95`
> - `P14_COVERAGE_FLOOR   = 0.20`
> - `P14_AUC_HARD_FLOOR   = 0.80`
> - `P14_OOD_MARGIN_FLOOR = 0.50`
>
> Also ruled in the same freeze: `P14_ALPHA_CONFORMAL = 0.10`, and `P14_N_CAL = P14_N_EVAL = 2000`.
>
> NOT ruled here and must follow the plan as written: τ (inherited by loading with provenance
> asserted, D-07), α_FDR (per-call user parameter), and the SC1-b / SC1-d bands (derived by formula).
>
> Record in `consts.jl`, at the point of definition, that these five bars are `[ASSUMED]` —
> pre-registered by explicit user ruling on 2026-08-03 without an analytic derivation, frozen BEFORE
> any Phase-14 result existed. State plainly in the same comment that if a bar is missed, the result
> is REPORTED and NOT re-tuned — matching the Phase-11 and Phase-12 precedent. Do not soften that.

**Option selected: `accept-as-proposed`.** No proposed value was changed, so no value carries a
"USER RULING replaced the proposal" note — but the *act of ruling* is itself recorded as a dated
`USER RULING` in `spike/p14/consts.jl`, because "the user ratified these five" is the fact that
distinguishes them from five numbers an agent picked.

**Consequences carried into Task 1:**

- `report-only` was NOT selected, so `14-VALIDATION.md`'s SC3-a..d rows stay **GATED** and neither
  that file nor the ROADMAP is touched by this plan.
- The ruling adds one clause the proposal did not contain: **a missed bar is REPORTED, never
  re-tuned.** That is stronger than the plan's own wording and is written into the constants file as
  binding text, not as commentary.
- `P14_ITERATION_ALLOWANCE` still concerns conformal coverage only. It cannot be spent on any of
  these five, and the file says so.

---

# APPENDED 2026-08-03 — PLAN 14-01 COMPLETE. WAVE 0 IS DONE.

## The freeze commit — the sha every later plan cites as ordering evidence

```
a6c825867dbc786f7c3927df6760295ee0930c77
pre-register(14-01): freeze the Phase-14 Tier-1 pre-registration
```

**It contains exactly one file** — `spike/p14/consts.jl`, 413 insertions. No `.jld2`, no
`spike/p14/run_p14_*.jl`, no measured number. Verified with `git show --stat HEAD`. That is the
acceptance criterion this task exists for: **every result-producing commit in Phase 14 must be a
descendant of `a6c8258`,** and the ordering is a matter of git history rather than trust.

## Commits, in order

| # | Task | Commit | Files |
|---|---|---|---|
| 0 | Rule the five undERIVED bars | `de9a6db` | `14-01-SUMMARY.md` (the ruling, recorded **before** anything was written to `spike/`) |
| 1 | Freeze the Tier-1 pre-registration | **`a6c8258`** | `spike/p14/consts.jl` (413 lines) |
| 2 | Assert the freeze | `0ff4669` | `spike/test/test_p14_consts.jl` (297 lines) |

Task 0 was committed separately and first, on purpose: the ruling has to exist on disk before the
values it ruled, or the record cannot distinguish "the user ratified these" from "an agent wrote
these and the user was told afterwards."

## Verification — actually run, with observed output

| Check | Command | Observed |
|---|---|---|
| Tier-1 includes cleanly, every `@assert` passes | `julia --project=spike -e 'include("spike/p14/consts.jl"); …'` | printed `tier-1 ok`, exit 0 |
| The freeze is executably asserted | `julia --project=spike spike/test/test_p14_consts.jl` | **273 pass, 0 fail, 0 error**, 0.2 s test time / 2.6 s wall including startup |
| τ is absent from non-comment source | `grep -v '^\s*#' spike/p14/consts.jl \| grep -c 'P14_TAU\|const P13_TAU'` | `0` |
| Guard sentinel is not a seed | `grep -c 'isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS)' …` | `2` |
| All 20 required constant names present as literal text | per-name `grep -q` | 20/20 ok |
| AGPL header byte-identical to `test_p12_decoupling.jl:1-19` | `diff <(sed -n '1,19p' …)` | identical |
| CPU-only testset present | `grep -c 'CUDA' spike/test/test_p14_consts.jl` | `1` |
| `src`, `spike/Project.toml`, `spike/Manifest.toml`, `corpus` byte-unchanged | `git diff --exit-code HEAD -- …` | exit 0 |

**Falsification check, run and reverted.** Flipping one digit of `P14_DEV_SEED`
(`…0B14_DE71` → `…0B14_DE72`) in place made `test_p14_consts.jl` exit **1**. The flip was chosen so
it does *not* collide with a forbidden seed — so the file's own self-checks still passed and it was
the **test's literal assertion** that caught it, which is what the criterion is actually about. The
scratch change was reverted with `git checkout --` and `git diff --exit-code` confirmed clean before
anything was staged.

Test breakdown as reported by Julia:

```
P14 Tier-1 pre-registration (D-07)                                               |  273    273  0.2s
  the locked literals                                                            |   22
  seeds are fresh and disjoint                                                   |   31
  streams are counter-separated                                                  |    9
  the Philox key-word cross-product is pairwise distinct vs the whole P13 family |  175
  tau is NOT Tier-1 here (D-07)                                                  |    5
  the two alphas are distinct bindings (D-07)                                    |    7
  the judgement-call bars are declared as such                                   |    8
  the iteration allowance is pre-declared and non-trivial                        |    6
  declared deviations and decoupling are on the record                           |    9
  P14 ran CPU-only                                                               |    1
```

## What is frozen

**Streams.** `P14_DEV_SEED = 0x…0B14_DE71`, `P14_FIX_SEED = 0x…0B14_F1F7`,
`P14_SALT = 0x5851_F42D_4C95_7F2D` (the MMIX LCG multiplier). The forbidden inventory is **derived**,
never retyped: `_p14_forbidden_list()` collects `_p13_forbidden()` — which itself *recomputes* both
frozen ship-gate seed families rather than trusting a comment — and adds Phase 13's own two seeds,
which Phase 13 could not forbid to itself.

**The disjointness proof goes past seed identity.** Two seeds that differ are not enough, because
what actually keys Philox in this project's runners is `seed ⊻ salt ⊻ counter`. The file asserts the
**full key-word cross-product**: 12 Phase-14 key words (dev × 5 counters, fix × 5 counters, plus both
bare) against 14 Phase-13 key words — 168 pairs, plus pairwise distinctness within Phase 14. The test
re-runs all of it and adds three draw-a-number checks against the P13 gate, datagen and fixture
streams.

**Counters.** `CAL = 1`, `EVAL = 2`, `OOD = 3`, `REAL = 4`, `FIXTURE = 99`. `P14_EVAL_COUNTER` is
recorded in the file as the **shared** reported set that SC1-b, SC1-d and SC3-a/b/c all read — shared
deliberately, and written down so the multiplicity is visible now rather than discovered later.

**Sizes and alphas.** `N_CAL = N_EVAL = 2000`, `N_OOD = 1000` per arm. `P14_ALPHA_CONFORMAL = 0.10`;
`P14_ALPHA_FDR_GRID = (0.01, 0.05, 0.10, 0.20)`. The test asserts the sd arithmetic
(`sqrt(0.10·0.90/2000) ≈ 0.0067`) rather than restating it in a comment.

**The five ruled bars** — `0.60 / 0.95 / 0.20 / 0.80 / 0.50` — sit under a comment block that says
they are JUDGEMENT CALLS WITH NO DERIVATION, names the four recorded cases of an underived bar
measuring the wrong thing, records the user's 2026-08-03 ruling, and states the binding clause: **if
a bar is missed the result is REPORTED and NOT re-tuned.** `test_p14_consts.jl` asserts the strings
`JUDGEMENT CALL`, `USER RULING`, `2026-08-03` and `NOT RE-TUNED` against the **raw** (un-stripped)
source, so the words cannot be deleted while the numbers stay.

**τ is not here.** `spike/p14/consts.jl` has no Tier-2 block and must never grow one: D-07 makes τ
*inherited by loading*, not appended. Its absence is enforced three ways — an `@assert` in the file,
`!isdefined` in the test, and a comment-stripped source grep.

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 3 — Blocking] The τ-absence assertion could not name τ literally**

- **Found during:** Task 1, immediately after the first draft of `consts.jl`.
- **Issue:** The plan asks (Task 1 acceptance) that
  `grep -v '^\s*#' spike/p14/consts.jl | grep -c 'P14_TAU\|const P13_TAU'` return **0**, and
  `test_p14_consts.jl` testset 5 asserts `!occursin("P14_TAU", P14_CONSTS_CODE)` over the
  comment-stripped source. Mirroring Phase 13's Tier-1 self-check literally —
  `@assert !isdefined(@__MODULE__, :P14_TAU) …` — writes that exact string on a **code** line, so the
  guard would have failed the very file it protects. The two requirements are contradictory in their
  literal form.
- **Fix:** assemble the needle at run time — `@assert !isdefined(@__MODULE__, Symbol("P14_", "TAU"))`
  — which is this repository's own idiom for exactly this situation:
  `spike/test/test_p12_decoupling.jl` builds three needles by concatenation "so the contiguous
  literal never appears in this source at all" (14-PATTERNS §3.11). The reason is written into the
  file beside the assertion rather than left for a reader to reconstruct.
- **Files modified:** `spike/p14/consts.jl`
- **Commit:** `a6c8258`

### Structural additions the plan implies but does not name

None of these is a choice about a value; each exists only so an assertion the plan *does* demand can
be expressed. Recorded so the file's surface is fully declared:

| Name | Why it exists |
|---|---|
| `P14_REPORTED_COUNTERS`, `P14_ALL_COUNTERS` | the fixture-counter and pairwise-distinctness assertions need the sets by name |
| `P14_P13_COUNTERS` | the Phase-13 counter family, held in the `P14_` namespace so this file never shadows a P13 name |
| `_p14_key_words()`, `_p14_p13_key_words()` | the cross-product the plan requires, expressed once and reused by the test instead of duplicated |
| `P14_SRC_UNTOUCHED` | mirrors `P13_SRC_UNTOUCHED`; makes D-01 machine-checked rather than asserted in prose |
| `@assert P14_ALPHA_CONFORMAL in P14_ALPHA_FDR_GRID` | a tripwire on the D-07 warning itself: if the two ever stop sharing 0.10, the comment explaining why they must never be assigned from one another has gone stale and must be re-read |

### Not done, deliberately

- **`14-VALIDATION.md` and `ROADMAP.md` untouched.** Their SC3-a..d rows stay **GATED**. That edit is
  conditional on the `report-only` option, which the user did not select.
- **`STATE.md` / `ROADMAP.md` not written** — the orchestrator owns those. (Both show as modified in
  the working tree; those edits are the orchestrator's, made concurrently with this run, and were
  neither made nor staged by this executor.)

## Assumption Drift (advisory)

**1. Execution environment — the divergence flagged by the previous run is RESOLVED**

- **Found during:** startup
- **Planned:** main working tree, no worktrees
- **Actual:** main working tree, no worktrees — this run was dispatched correctly
- **Why it matters:** the earlier BLOCKED record above flagged a worktree dispatch and warned it
  would break later Phase-14 plans that read the gitignored `spike/data/cache/p13/` pool. That
  advisory stands for the *later* plans; it no longer applies to this run. `git worktree` was never
  invoked.

Nothing else drifted materially. Every value written matches the ruling and the plan.

## Constraint Verification (all hold, all checked on disk)

| Constraint | Check | Status |
|---|---|---|
| `src/` untouched (D-01) | `git diff --exit-code HEAD -- src` | exit 0 |
| `spike/Project.toml` + `Manifest.toml` byte-unchanged (`test_p12_decoupling.jl:176-178`) | `git diff --exit-code HEAD -- …` | exit 0 |
| No dependency added or changed | only `Test` (stdlib) + `Random123` (already used by p13) | ok |
| Phase-16 sealed holdout not consumed (`test_p12_decoupling.jl:210`) | `corpus/` never read; `git diff --exit-code HEAD -- corpus` | exit 0 |
| DEV seed disjoint from every forbidden stream | executable, in the file **and** in the test | 31 assertions pass |
| Corrections APPENDED, never overwritten | this summary appends; `consts.jl` is a first write | ok |
| `runtests.jl` not wired as a gate | mentioned once, in a comment saying it is *not* a gate | ok |
| `STATE.md` / `ROADMAP.md` not modified by this executor | not staged, not written | ok |
| Files modified outside the allowed three | none | ok |

## Known Stubs

None. Neither file contains a placeholder, a TODO, or a value awaiting a later run. τ's absence is
not a stub — it is the D-07 requirement, asserted three ways.

## Threat Flags

None. No network surface, no file write at runtime, no new trust boundary. The two threats this plan
was written to mitigate are both discharged:

- **T-14-52** (an agent silently choosing five gating numbers) — the blocking checkpoint was honoured
  by the previous run, the ruling was committed at `de9a6db` **before** `spike/p14/` existed, and the
  ruling is reproduced as a dated `USER RULING` comment at the point of definition.
- **T-14-02** (post-hoc tampering with a bar) — every bar is asserted by literal value, and the
  falsification check above confirms the assertions bite.

## Self-Check

- `spike/p14/consts.jl` — **FOUND** (413 lines)
- `spike/test/test_p14_consts.jl` — **FOUND** (297 lines)
- commit `de9a6db` — **FOUND**
- commit `a6c8258` — **FOUND**, single file, no artifact
- commit `0ff4669` — **FOUND**
- `git status --porcelain` limited to the three permitted paths — **clean** (the only other modified
  files are `.planning/STATE.md` and `.planning/config.json`, both written by the orchestrator, not
  by this executor)

## Self-Check: PASSED

**Wave 0 is complete.** The pre-registration exists, was ruled by the user before it was written, is
committed at `a6c8258`, and is executably asserted. No result-producing code exists yet, by
construction. The five judgement-call bars are on disk, labelled as judgement calls, before any
number they will score.

Later Phase-14 plans may now run. Every one of them must commit as a descendant of `a6c8258`.

