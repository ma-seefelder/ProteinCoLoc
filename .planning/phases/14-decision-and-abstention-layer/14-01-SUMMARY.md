---
phase: 14
plan: 01
subsystem: pre-registration
status: BLOCKED — pending user ruling at Task 0 (blocking checkpoint:decision)
tags: [pre-registration, seeds, judgement-call-bars, wave-0]
requires: []
provides: []          # nothing yet — Task 0 gates Tasks 1 and 2
affects: []
tech-stack:
  added: []           # no dependency added; spike/Project.toml and Manifest.toml byte-unchanged
  patterns: []
key-files:
  created: []         # spike/p14/consts.jl NOT written — deliberately
  modified: []
decisions: [D-07]
metrics:
  duration: ~10 min (read-only)
  completed: null
  tasks_completed: 0
  tasks_total: 3
---

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
