---
phase: 14-decision-and-abstention-layer
plan: 14
subsystem: decision-and-abstention-layer
status: paused
pause_reason: "blocking checkpoint:human-verify (Task 2) — the two honesty items are the user's to ratify, not an executor's"
tags: [report, honesty-items, sc-table, amendments-cited, judgement-call-bars, gate-not-met, blocking-checkpoint]
requires:
  - "14-09 (p14_conformal_report.jld2, p14_eval_pool.jld2)"
  - "14-10 (p14_fdr_report.jld2)"
  - "14-11 (p14_riskcoverage_report.jld2)"
  - "14-12 (p14_ood_arm_report.jld2)"
  - "14-13 (p14_real_images_report.jld2)"
  - "spike/p14/result.jl (p14_named_limits — the limit list enumerated from code, not from memory)"
provides:
  - ".planning/phases/14-decision-and-abstention-layer/14-REPORT.md — every SC row with its measured value, verdict and artifact key; both amendments cited beside the originals; the manuscript-bound honesty items"
affects:
  - "any manuscript claim citing a Phase-14 number: the scope, the amendment pairing, the simulator provenance and the judgement-call labels now travel in one document"
tech-stack:
  added: []
  patterns: ["every number loaded via JLD2.load and pasted, never retyped", "artifact sha256 recorded beside every section", "named limits enumerated from code", "amendment quoted beside the original wording", "the failed gate stated at the same volume as the passes"]
key-files:
  created:
    - ".planning/phases/14-decision-and-abstention-layer/14-REPORT.md (781 lines)"
  modified: []
decisions: [D-01, D-02, D-03, D-03a, D-04, D-05, D-06, D-07]
validates: ["SC1-a", "SC1-b", "SC1-c", "SC1-d", "SC1-e", "SC1-f", "SC2-a", "SC2-b", "SC2-c", "SC2-d", "SC3-a", "SC3-b", "SC3-c", "SC3-d", "D-07", "D-01/env", "Seeds", "Class order", "Prior re-derivation"]
metrics:
  duration: "~55 min wall"
  completed: null
---

# Phase 14 Plan 14: The Phase-14 Report — Summary (PAUSED AT CHECKPOINT)

**Task 1 is complete and committed (`d57e39c`). Task 2 is a BLOCKING `checkpoint:human-verify` and
this plan stops there.** The orchestrator is running unattended, so there is no user available to
approve inline. Nothing was self-approved and no user response was invented.

`14-REPORT.md` exists at 781 lines. Every number in it was loaded out of a `.jld2` artifact with
`JLD2.load` and pasted; no value was retyped from a plan summary. Every section names the artifact it
was read from and every artifact's sha256 is recorded at the top.

## The headline the report carries

**Seventeen SC rows PASS, one is REPORTED-NOT-GATED, one FAILS.**

| Row | Measured | Verdict |
|---|---|---|
| SC1-b realized-vs-predicted FDP | within the binomial 95 % band at all four α; decided fraction **0.888** (1776 of 2000) | PASS |
| SC1-c prior sensitivity | at π_coloc = 0.70, α = 0.20 the rule predicts 0.1997520196632641 and realizes 0.3466666666666667 | REPORTED, NOT GATED |
| SC1-d conformal coverage | **0.9115** vs band **0.8798753882025019**; `q̂ = 0.4498001628978108` | PASS |
| SC1-f real-image check | both specimens ABSTAIN (`:ood_fired`), density 703.2995398773389 vs 167.5446329378043 | PASS (structural) |
| SC3-a / b / c | skill **0.9017113386798696** / Spearman **0.9701942302111639** / AUC_hard **0.9117431863910738** | PASS |
| **SC3-d abstain-rate margin** | **−0.129** vs the frozen floor **0.5**; `sc3d_met = false` | **FAIL — NOT MET** |

**SC3-d is stated in the report as plainly as the passes**, in §3 (the table), §6 (its own subsection
with the per-family breakdown, the driver analysis and the cross-tab) and §8 (limit B). No constant in
`spike/p14/consts.jl` was relaxed, reinterpreted or touched — it is **byte-unchanged against HEAD**,
verified with `git diff --exit-code` after the report was written.

## The report's ten sections, and what each is for

| § | Content | Backing artifact |
|---|---|---|
| 1 | What was built; D-01 stated early — src-shaped is not in-src, the layer is NOT shipped in v2.0; the six guarded-tree assertions | git + `test_p14_decoupling.jl` |
| 2 | SC1 and SC2 quoted **as originally written** beside their amendments, with the reason for each and the considered-and-rejected `DECIDED-CONTRA-CLASSICAL` four-state design | ROADMAP + `amendment` key |
| 3 | The 19-row SC table: ID, behaviour, pre-registered number, measured value, verdict, artifact key, reproduce command. **No row without a verdict.** | all five |
| 4 | FDR at the scope it holds: the four-row α table **with `decided_fraction` on every row**, the verbatim claim, Assumption A2 as a two-line derivation, the ordering rule, the 20-row π table, the prior-atom asymmetry | `p14_fdr_report.jld2` |
| 5 | The hedge: `q̂`, coverage vs band, the set-size distribution with the zero-ambiguous qualifier, both class-mass tables, the cross-tab, and the simulator-derived / absent-not-loose statement | `p14_conformal_report.jld2` |
| 6 | Risk-coverage, the raw full-range curve, all five bars labelled JUDGEMENT CALLS with the freeze sha, the licensing standard from STATE.md, and SC3-d's failure with driver analysis and cross-tab | `p14_riskcoverage_report.jld2`, `p14_ood_arm_report.jld2` |
| 7 | The six TIFFs as an illustration, the corpus as unfetched by design, the D-03a substitution, the withdrawn figure in prose | `p14_real_images_report.jld2` |
| 8 | `p14_named_limits()` enumerated **from code** (8 items) plus four further limits: the epoch-4 checkpoint, the single wired OOD channel, the dirty-tree τ probe, and per-region FDR out of scope | `spike/p14/result.jl` |
| 9 | Seven forward-pointing items for Phases 15/16, each naming its evidence | — |
| 10 | Per-file reproduction with measured pass counts and the standing `runtests.jl` warning | measured |

## Verification performed (observed, not assumed)

- **The plan's automated honesty-token check ran verbatim and exited 0** — all eleven required tokens
  present (`DECIDED SUBSET ONLY`, `src-shaped`, `NOT shipped in v2.0`, `JUDGEMENT CALL`,
  `illustration, not a coverage claim`, `absent, not merely loose`, `SC1 is AMENDED`,
  `SC2 is AMENDED`, `Assumption A2`, `unfetched by design`, `epoch-4`) and **both withdrawn-figure
  numerals absent**. Printed `report honesty tokens ok`.
- `wc -l` → **781 lines** (≥ 220 required).
- `grep -c "p14_named_limits"` → 1; `grep -cE "p14_.*_report\.jld2"` → 22. Both `key_links` patterns
  satisfied.
- `git diff --exit-code HEAD -- src spike/Project.toml spike/Manifest.toml corpus test` → **exit 0**.
- `git diff --exit-code HEAD -- spike/p14/consts.jl` → **exit 0, byte-unchanged**.
- **All ten `test_p14_*.jl` files run individually — all exit 0, 1178 assertions total**
  (consts 273, decoupling 51, provenance 77, posterior 48, fdr 55, conformal 83, fuse 222, result 121,
  decide 146, pools 102). `spike/test/runtests.jl` was **not** used.
- Every α row in §4 carries a populated `decided_fraction` column (hand-scanned: 0.888 × 4).
- Every SC row in §3 has a non-empty verdict cell; no row reads `TBD`.

## Deviations from Plan

### Auto-fixed / recorded

**1. [Rule 1 — Bug found while assembling] The `1/311` fraction would have tripped the report's own
honesty grep.**

- **Found during:** Task 1, writing §4.
- **Issue:** The α = 0.01 row's realized FDP is one false discovery among 311 accepted. Written the
  way 14-09-SUMMARY writes it — `0.003215 (1/311)` — the literal string contains `1/31`, which the
  plan's own `<verify>` asserts **absent** because it is the withdrawn D-03 figure. The check would
  have failed on a number that has nothing to do with the withdrawn figure.
- **Fix:** the α table carries `n_accepted` and `n_false_discoveries` as separate columns, so the
  count is fully recoverable without ever writing the fraction.
- **Files modified:** `14-REPORT.md`. **Commit:** `d57e39c`.

**2. [Rule 2 — Missing disclosure] A CRLF/LF blob-sha discrepancy that would read as tampering.**

- **Found during:** Task 1, cross-checking the recorded provenance against git.
- **Issue:** every artifact records `consts_git_blob_sha = 26ac1e7b2b5e85badd3e2dc9f99238fff21bc034`,
  but `git rev-parse HEAD:spike/p14/consts.jl` returns `c1dbd185435f975ff237dae6e9d9eca48cd57134`. A
  reader checking the provenance would see a mismatch on the pre-registration file and reasonably
  suspect an edit.
- **Diagnosis (measured, not argued):** `git hash-object --no-filters` on the working file returns the
  artifact's value; `git hash-object` (filtered) returns the committed value. The runner's helper
  hashes raw working-tree bytes, and git normalizes CRLF→LF on commit. The content is identical; the
  sha256 `4e8ee4be…` agrees on both sides.
- **Fix:** recorded as forward-pointing item 7 in §9 rather than silently omitted.
- **Files modified:** `14-REPORT.md`. **Commit:** `d57e39c`.

### Assumption Drift (advisory)

**1. The plan frames §7 around "whether the ABSTAIN prediction held"; the more interesting measured
fact is that the two misspecification channels DISAGREE on that substrate.**

- **Planned:** report the six decisions, the abstention reasons, and whether the prediction held.
- **Actual:** the prediction held exactly, *and* both conformal sets are singletons while the density
  null fires at 4.198× its operating point. The exchangeability channel saw nothing unusual on the
  only real microscopy this project holds.
- **Why it matters:** a reader taking "both ABSTAIN" as "the layer detected that this is real data"
  would be over-reading — one wired channel detected it and the other did not. Recorded in §7 rather
  than editorialised.
- Advisory only. Nothing gated, no threshold touched.

## Known Stubs

None. The report is a document; every section is populated from a loaded artifact.

## Threat Flags

None. This plan adds no network endpoint, auth path, file-write surface or schema at a trust boundary.
It reads six `.jld2` files and writes one Markdown file under `.planning/`.

## Self-Check: PASSED

- `.planning/phases/14-decision-and-abstention-layer/14-REPORT.md` — FOUND (781 lines)
- `spike/p14/p14_conformal_report.jld2` — FOUND
- `spike/p14/p14_fdr_report.jld2` — FOUND
- `spike/p14/p14_riskcoverage_report.jld2` — FOUND
- `spike/p14/p14_ood_arm_report.jld2` — FOUND
- `spike/p14/p14_real_images_report.jld2` — FOUND
- `spike/p14/p14_eval_pool.jld2` — FOUND
- commit `d57e39c` (Task 1) — FOUND

---

## CHECKPOINT REACHED — Task 2 is BLOCKING and was NOT self-approved

**Type:** `checkpoint:human-verify`, `gate="blocking"`
**Progress:** 1 of 2 tasks complete.

The plan's acceptance criteria require the user's answers to two questions to be recorded **verbatim**
in this file. They are **not recorded** because they have not been given. This section is the record
of what is outstanding.

### Question 1 — does §4's Assumption A2 derivation convince you as written?

`14-RESEARCH.md`'s Assumptions Log flags A2 **High risk if wrong**. It is the entire justification for
abstain-then-sort and for the "controlled over the decided subset" claim. §4 states it as two
checkable lines: the rejection set is a deterministic function of the observed data, therefore
σ(data)-measurable, therefore it pulls out of the conditional expectation, so the identity holds for
*any* data-dependent selection — unlike the frequentist case, where pre-selection is a genuine
selective-inference problem. **Your call: does that read as a derivation a referee can check, or does
it need a referee-facing rewrite before it goes near the manuscript?**

### Question 2 — are you content with how the five judgement-call bars are presented?

`skill ≥ 0.6`, `spearman ≥ 0.95`, `AUC_hard ≥ 0.8`, `ood margin ≥ 0.5`, and the `0.2` coverage floor.
All five were frozen in `a6c825867dbc786f7c3927df6760295ee0930c77` (a commit containing no Phase-14
result), you ruled them accept-as-proposed on 2026-08-03, and **none was relaxed** — including after
SC3-d missed. §6 labels all five JUDGEMENT CALLS WITH NO DERIVATION and records the licensing
standard. **They cannot be changed now — a bar changed after a result is a pre-registration breach —
but how they are *presented* in the manuscript is yours to rule.** Given this project has recorded
four cases of an underived bar measuring the wrong thing, are you content with the current framing?

### Question 3 — how do you want the FAILED gate recorded?

SC3-d is NOT MET at **−0.129** against **0.5**, carried entirely by the `noise` family (0 of 1000
flagged) — the family the unwired `(:pp, :noise)` channel exists for. Nothing was tuned to try to
clear it. `P14_ITERATION_ALLOWANCE`'s single pre-declared trigger (conformal coverage *below* the
SC1-d band) did not fire and the allowance is unspent, so no in-phase remedy is available. **Say how
you want this recorded when you close the phase. No bar may be relaxed.**

### Mechanical checks the user was asked to run — already run by the executor, all green

| Check | Result |
|---|---|
| ten `test_p14_*.jl` individually | all exit 0, 1178 assertions |
| `spike/test/runtests.jl` used? | **no** |
| `git status --porcelain -- src spike/Project.toml spike/Manifest.toml corpus test` | empty |
| report contains every honesty token | yes |
| report contains the two withdrawn numerals | **no** |
| `spike/p14/consts.jl` byte-unchanged | yes |

### Resume signal

Type **"approved"** to close the phase, or describe what must change. If changes are requested they
are applied and this checkpoint is re-run; the phase does not close on an unresolved objection.
