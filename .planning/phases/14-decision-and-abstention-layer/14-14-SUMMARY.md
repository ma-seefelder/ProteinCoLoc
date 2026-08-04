---
phase: 14-decision-and-abstention-layer
plan: 14
subsystem: decision-and-abstention-layer
status: complete
pause_reason: ""
tags: [report, honesty-items, sc-table, amendments-cited, judgement-call-bars, gate-not-met, checkpoint-resolved, verdict-sensitivity-strip, sc3d-ruling-applied]
requires:
  - "14-09 (p14_conformal_report.jld2, p14_eval_pool.jld2)"
  - "14-10 (p14_fdr_report.jld2)"
  - "14-11 (p14_riskcoverage_report.jld2)"
  - "14-12 (p14_ood_arm_report.jld2)"
  - "14-13 (p14_real_images_report.jld2)"
  - "spike/p14/result.jl (p14_named_limits — the limit list enumerated from code, not from memory)"
  - "14-SC3D-RULING.md (f4934a9 — the user's separately-ruled disposition of the FAILED SC3-d row)"
  - "artifacts/amended_v2/grid_8/gate_report_8.jld2 (READ-ONLY: report[:ood], the shipped fused detector's per-family AUCs)"
provides:
  - ".planning/phases/14-decision-and-abstention-layer/14-REPORT.md — every SC row with its measured value, verdict and artifact key; both amendments cited beside the originals; the verdict-sensitivity strip; the six-point SC3-d disposition; the manuscript-bound honesty items"
  - "spike/p14/run_p14_bar_sensitivity.jl + p14_bar_sensitivity.jld2 — the READ-ONLY verdict-sensitivity strip; gates nothing, changes no bar"
affects:
  - "any manuscript claim citing a Phase-14 number: the scope, the amendment pairing, the simulator provenance and the judgement-call labels now travel in one document"
  - "Phase 15's noise axis: a PROTECTIVE prediction is now on the record BEFORE that phase runs, and fused_auc[:background] = 0.0 at level 3 is flagged for inspection first"
tech-stack:
  added: []
  patterns: ["every number loaded via JLD2.load and pasted, never retyped", "artifact sha256 recorded beside every section", "named limits enumerated from code", "amendment quoted beside the original wording", "the failed gate stated at the same volume as the passes", "a sensitivity strip that gates nothing, with its verdict at the frozen bar PINNED by assertion to the gated runner's own persisted flag"]
key-files:
  created:
    - ".planning/phases/14-decision-and-abstention-layer/14-REPORT.md (781 lines at Task 1; 992 after the checkpoint rulings were applied)"
    - "spike/p14/run_p14_bar_sensitivity.jl (the READ-ONLY verdict-sensitivity strip; gates nothing)"
  modified:
    - ".planning/phases/14-decision-and-abstention-layer/14-REPORT.md (§3 tally note, §6 strip + SC3-d disposition, §8 limit B, §9 items 8–9, §10 runner list and pass counts)"
decisions: [D-01, D-02, D-03, D-03a, D-04, D-05, D-06, D-07]
validates: ["SC1-a", "SC1-b", "SC1-c", "SC1-d", "SC1-e", "SC1-f", "SC2-a", "SC2-b", "SC2-c", "SC2-d", "SC3-a", "SC3-b", "SC3-c", "SC3-d", "D-07", "D-01/env", "Seeds", "Class order", "Prior re-derivation"]
metrics:
  duration: "~55 min wall (Task 1) + ~35 min wall (checkpoint resolution)"
  completed: "2026-08-04T15:47:29Z"
---

# Phase 14 Plan 14: The Phase-14 Report — Summary (COMPLETE)

**Task 1 was completed and committed (`d57e39c`). Task 2 was a BLOCKING `checkpoint:human-verify`
and this plan stopped there — nothing was self-approved and no user response was invented.** The
user has now answered, and this summary records those answers verbatim and what changed as a result.
Task 2 is **RESOLVED**; the plan is **complete**.

`14-REPORT.md` was 781 lines at Task 1 and is **992 lines** after the two rulings were applied. Every
number in it was loaded out of a `.jld2` artifact with `JLD2.load` and pasted; no value was retyped
from a plan summary, from memory, or from a ruling document. Every section names the artifact it was
read from and every artifact's sha256 is recorded at the top.

---

## Checkpoint resolved — the user's answers, verbatim

**Question 1 — the Assumption A2 derivation in §4: does it convince as written, or does it need a
referee-facing rewrite before it goes near the manuscript?**

> **USER'S ANSWER:** "Accept as written. The two steps are the whole argument and both are checkable;
> the independence question — the one a referee actually asks — is answered in place. It goes into
> the phase record unchanged, and any manuscript-level polish happens in Phase 16 when the prose is
> written for a referee rather than for a phase report."

**What changed as a result: nothing.** §4's A2 derivation is **unchanged, byte for byte**, and no
rewrite was attempted. The manuscript-level polish is deferred to Phase 16 by the user's own
instruction, which is where the prose is written for a referee.

**Question 2 — the framing of the five judgement-call bars in §6: is the current presentation enough,
given four recorded cases of an underived bar measuring the wrong thing?**

> **USER'S ANSWER:** "Add a verdict-sensitivity strip. Keep every bar frozen and every verdict as-is,
> but add a read-only table showing what each SC3 verdict would have been across a range of nearby
> bar values. This changes nothing measured — it lets a reader see whether a conclusion hinges on the
> arbitrary number or is robust to it, which is the actual defence against an underived bar. It is
> presentation, not relaxation."

**What changed as a result: Work Item A below.** A new §6 subsection, *"Verdict sensitivity — what
each verdict would have been at a different bar (NO BAR WAS CHANGED)"*, computed by a new read-only
runner. **Every bar is still frozen and every verdict is still the verdict at the frozen bar.**

**Question 3 — how the FAILED SC3-d row is recorded — was ruled SEPARATELY by the user on
2026-08-04** and is on disk at `.planning/phases/14-decision-and-abstention-layer/14-SC3D-RULING.md`
(commit `f4934a9`). Its ruling: **Phase 14 closes NEGATIVE-BUT-USEFUL on that row; SC3-d is recorded
as a named limit SCOPED TO THE DENSITY-ONLY WIRING, not as a limit of the shipped OOD flag.** The
ruling's six-point instruction under *"HOW `14-REPORT.md` AND `14-14-SUMMARY.md` MUST STATE IT"* is
implemented as **Work Item B** below. **All three checkpoint items are now closed.**

---

## Work Item A — the verdict-sensitivity strip (Question 2's ruling), EXECUTED

**Commit `6ad9fab`.** New runner `spike/p14/run_p14_bar_sensitivity.jl`; new §6 subsection in
`14-REPORT.md`; artifact `spike/p14/p14_bar_sensitivity.jld2` (gitignored, regenerable).

**Nothing measured changed.** The runner **gates nothing** (`gates_nothing = true`,
`no_bar_changed = true`, `recomputed_anything_gated = false`), re-simulates nothing, re-trains
nothing and re-measures nothing. It loads the already-measured `skill = 0.9017113386798696`,
`spearman = 0.9701942302111639`, `auc_hard = 0.9117431863910738` and SC3-d `margin = −0.129` and
compares each against a sweep of **hypothetical** bar values. The frozen bars are read from
`consts.jl` **by name** and never retyped. `spike/p14/consts.jl` is **byte-unchanged** — the runner
asserts that itself, at run time, before computing anything.

**Three assertions make the strip trustworthy rather than merely plausible:**
- each sweep grid is asserted to **contain** its frozen bar, so the gating value always appears
  inside its own strip;
- the strip's verdict at each frozen rung is **pinned** to the `sc3a_met` / `sc3b_met` / `sc3c_met` /
  `sc3d_met` flags the gated runners persisted — if the strip scored under a different rule it fails
  loudly;
- the coverage-floor re-derivation is **pinned** at the frozen truncation to the persisted
  `spearman` (agreement < 1e-12) and `spearman_n_points = 1601`.

**What the strip actually shows** (this is the point of it, so it is stated and not buried):

| Gate | Measured | Frozen bar | Flips at | Distance | Reading |
|---|---|---|---|---|---|
| SC3-a skill | 0.9017 | 0.60 | 0.9017 | **+0.3017** | far from any boundary |
| SC3-b Spearman | 0.9702 | 0.95 | 0.9702 | **+0.0202** | **the one PASS with modest room** — a floor of 0.98 would have turned it red |
| SC3-c `AUC_hard` | 0.9117 | 0.80 | 0.9117 | **+0.1117** | far from any boundary |
| SC3-d margin | −0.1290 | 0.50 | −0.1290 | **−0.6290** | fails at **every** rung including **0.00** — no bar choice rescues it |

`P14_COVERAGE_FLOOR` differs in kind — it is a **truncation**, not a pass/fail bar — so its
sensitivity is a genuine re-derivation of SC3-b's Spearman at other truncations from the persisted
raw full-range curve. **It was computable from the persisted artifact** (`curve_coverage`,
`curve_selective_risk`); nothing was fabricated. Re-derived: 0.9368 at truncation 0.02, 0.9427 at
0.05, 0.9523 at 0.10, 0.9616 at 0.15, **0.9702 at the frozen 0.20**, rising to 0.9989 at 0.70.
**SC3-b would have been NOT MET at a truncation of 0.05 or 0.02** and is MET from 0.10 upward — the
frozen 0.20 sits inside the passing region with one rung of margin below it. That the floor could not
have been chosen to land there is checkable: it was frozen in
`a6c825867dbc786f7c3927df6760295ee0930c77`, a commit containing no Phase-14 result of any kind.

**Honest statement of what a one-sided strip can and cannot be**, recorded in the artifact and in the
report: every SC3 gate has the form `statistic ≥ bar`, so the flip point **is** the measured
statistic by construction. The strip's informative content is the **distance** from the frozen bar to
it. It is not new evidence and is not presented as any.

---

## Work Item B — the SC3-d ruling applied to `14-REPORT.md`, EXECUTED

**Commit `6ad9fab`.** New §6 subsection *"The SC3-d disposition — ruled by the user on 2026-08-04"*,
plus §3's tally note, §8 limit B and §9 items 8–9. All six numbered points of the ruling implemented:

1. **FAILED, unsoftened.** Margin **−0.129** vs the frozen **0.50**, `sc3d_met = false`, `noise`
   **0 of 1000**, the abstain-rate table intact. The sensitivity strip adds that the margin is
   *negative*, so no non-negative bar would have passed it.
2. **The scope, immediately beside it.** `ood_channels_wired = (:density,)`,
   `ood_channels_not_wired = (:pp, :noise)`, cited to `decide.jl:583-584,979-980`,
   `pools.jl:487,506,563-564` and `decide.jl:286,293`. The `density_auc` / `noise_auc` / `fused_auc`
   table is **LOADED** from `artifacts/amended_v2/grid_8/gate_report_8.jld2`, key **`report[:ood]`**,
   by the Work-Item-A runner and persisted — **not retyped from the ruling document.** Loading it
   produced *more precise* values than the ruling's rounded table (`texture` density level 4 is
   **0.99475**, not 0.995; `background` is **0.55925 / 0.95875 / 0.0 / 0.960125**; `optics` is
   **0.58 / 0.6755 / 0.656875 / 0.66975**), which is exactly why the instruction said to load it.
   `density_auc[:noise] = 0.0` at **all four** levels — the density channel is worse than chance on
   that family — against `fused_auc[:noise] = 1.0` at all four.
3. **The honest one-line form, quoted from the ruling:** *"the Phase-14 decision layer's abstention
   is blind to detector-noise misspecification because it wires one of the shipped flag's three
   channels — not because the shipped flag is blind to it."*
4. **The fused detector's SC3-d margin is UNMEASURED**, stated explicitly, and the report **does not
   assert it would pass** (`fused_sc3d_margin_measured = false` persisted). Phase 15 tests the fused
   detector independently (D-02 / D-06 / D-07) and its three-valued verdict is what settles it.
5. **`P14_ITERATION_ALLOWANCE` is UNSPENT**, with the reason: its single pre-declared trigger is
   conformal coverage *below* the SC1-d band, and SC1-d measured **0.9115** against a band of
   **0.8798753882025019** — above, not below. `iteration_allowance_applies_here = false`.
6. **The rejected alternatives, with reasons**, as a table so the refusal to re-wire is visible as a
   choice: routing to Phase 15 with a fresh bar (measures a differently-equipped detector; collides
   with Phase 15's D-14), re-wiring `(:pp, :noise)` and re-measuring (**would likely clear the bar,
   and that is precisely why it is refused**), and recording SC3-d as a general limit of the OOD flag
   (falsifiable by opening the shipped gate report).

**Carried into §9 as forward-pointing items:** item **8**, the **PROTECTIVE** prediction for Phase
15's `noise` axis, recorded **before** that phase runs so a `LATE` outcome would be a genuine
surprise; item **9**, the flagged `fused_auc[:background] = 0.0 at level 3` against
`noise_auc[:background] = 1.0` there — **flagged, not diagnosed**, so it is not discovered inside a
reported run.

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
| 6 | Risk-coverage, the raw full-range curve, all five bars labelled JUDGEMENT CALLS with the freeze sha, the licensing standard from STATE.md, **the verdict-sensitivity strip (added by Work Item A)**, SC3-d's failure with driver analysis and cross-tab, and **the six-point SC3-d disposition (Work Item B)** | `p14_riskcoverage_report.jld2`, `p14_ood_arm_report.jld2`, `p14_bar_sensitivity.jld2`, `gate_report_8.jld2` |
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

### Re-verification after the checkpoint rulings were applied (observed, not assumed)

Every check above was re-run against the amended report. All green:

- **The honesty-token check ran verbatim again and exited 0** — all eleven tokens still present, and
  **both withdrawn numerals still absent.** The new sensitivity strip emits no value whose printed
  form contains `0.032` or `1/31`; the check was **not** weakened to accommodate it and nothing was
  reformatted to dodge it — it simply does not collide. (The nearest neighbour anywhere in the
  document remains SC1-d's pre-existing coverage gap `+0.03162461179749809`, which does not match.)
- `git diff --exit-code HEAD -- spike/p14/consts.jl` → **exit 0, byte-unchanged.** The new runner
  additionally asserts this itself, at run time, before it computes anything.
- `git status --porcelain -- src spike/Project.toml spike/Manifest.toml corpus test` → **empty.**
- **All ten `test_p14_*.jl` files run individually — all exit 0, 1179 assertions total**
  (consts 273, decoupling 51, **provenance 78**, posterior 48, fdr 55, conformal 83, fuse 222,
  result 121, decide 146, pools 102). `spike/test/runtests.jl` was **not** used.
- **The +1 in the provenance count is explained rather than restated.** It was 77; it is 78.
  `test_p14_provenance.jl:73` builds `P14_LANE_FILES` by **globbing** `spike/p14/*.jl` precisely so a
  file added later is scanned automatically. Adding `run_p14_bar_sensitivity.jl` therefore adds
  **exactly one** assertion — the lane-wide guard that no Phase-14 source writes τ as a hardcoded
  literal, now applied to the new file, which passes. No test was edited or added by hand. The delta
  and its cause are recorded in §10 of the report.
- The strip runner's own internal pins all held: the verdict at each frozen bar matched the gated
  runners' persisted `sc3a/b/c/d_met`, and the Spearman re-derived at the frozen truncation matched
  the persisted `0.9701942302111639` to < 1e-12 on the same 1601 curve points.
- `artifacts/amended_v2/grid_8/gate_report_8.jld2` **exists** and was opened read-only; the six
  `.jld2` reads and the one `artifacts/` read produced every number in Work Items A and B.
  `git status --porcelain -- artifacts` is untouched by this plan.

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

**3. [Recorded, not auto-fixed] This plan ran on the MAIN WORKING TREE, not in an isolated
worktree.**

- **Found during:** plan start (both the Task-1 run and the checkpoint-resolution run).
- **Deviation:** GSD's default execution isolation is a git worktree. This plan ran directly on
  `gsd/v2.0-milestone` in the main checkout instead, as its own frontmatter declares
  (`execution_environment: main-working-tree-no-worktrees`).
- **Why:** every `.jld2` this plan reads is **gitignored** (`.gitignore:448`, `spike/p14/*.jld2`) and
  exists only in this checkout. A fresh worktree would contain none of them, and the plan's central
  discipline — *every number loaded from an artifact, never retyped* — would have been impossible to
  honour there. The same applies to `spike/data/cache/p13/` (~52 MB), and this milestone has already
  lost a 54 MB pool once to a worktree cleanup.
- **Containment:** the guarded tree (`src`, `spike/Project.toml`, `spike/Manifest.toml`, `corpus`,
  `test`) was verified clean before and after every commit, `spike/p14/consts.jl` is byte-unchanged,
  and each runner asserts the D-01 lane guard at run time rather than only in review.
- **Files modified:** none by the deviation itself. **Commits:** `d57e39c`, `6ad9fab`, and this
  summary's commit.

**4. [Recorded] The provenance test's assertion count moved 77 → 78 as a direct, expected
consequence of adding a file to `spike/p14/`.**

- **Found during:** the post-ruling re-verification.
- **Issue:** the report's §10 recorded 77 for `test_p14_provenance.jl`; after Work Item A it measures
  78. A silently restated count would look like a number that drifted.
- **Diagnosis (measured, not argued):** `test_p14_provenance.jl:73` globs `spike/p14/*.jl` rather
  than enumerating them, by design, so the lane-wide "no hardcoded τ literal" guard arms itself on
  any new file. One new file ⇒ exactly one new assertion, which passes.
- **Fix:** §10 now states **78 / 78** and **1179 / 1179** with the cause spelled out beneath the
  table. **Commit:** `6ad9fab`.

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
It reads seven `.jld2` files — the six Phase-14 artifacts plus the shipped
`artifacts/amended_v2/grid_8/gate_report_8.jld2`, the latter strictly **read-only** — and writes one
Markdown file under `.planning/`, one runner under `spike/p14/`, and that runner's own gitignored
`spike/p14/p14_bar_sensitivity.jld2`. Nothing under `artifacts/`, `src/` or `corpus/` is written.

## Self-Check: PASSED

- `.planning/phases/14-decision-and-abstention-layer/14-REPORT.md` — FOUND (992 lines)
- `.planning/phases/14-decision-and-abstention-layer/14-SC3D-RULING.md` — FOUND
- `spike/p14/run_p14_bar_sensitivity.jl` — FOUND
- `spike/p14/p14_bar_sensitivity.jld2` — FOUND (gitignored, regenerable)
- `spike/p14/p14_conformal_report.jld2` — FOUND
- `spike/p14/p14_fdr_report.jld2` — FOUND
- `spike/p14/p14_riskcoverage_report.jld2` — FOUND
- `spike/p14/p14_ood_arm_report.jld2` — FOUND
- `spike/p14/p14_real_images_report.jld2` — FOUND
- `spike/p14/p14_eval_pool.jld2` — FOUND
- `artifacts/amended_v2/grid_8/gate_report_8.jld2` — FOUND (read-only)
- commit `d57e39c` (Task 1) — FOUND
- commit `f4934a9` (the SC3-d ruling) — FOUND
- commit `6ad9fab` (Work Items A and B) — FOUND

---

## CHECKPOINT RESOLVED — Task 2 was BLOCKING, was NOT self-approved, and has now been answered

**Type:** `checkpoint:human-verify`, `gate="blocking"`
**Progress:** **2 of 2 tasks complete.**

The plan's acceptance criteria require the user's answers to be recorded **verbatim** in this file.
They are recorded above under *"Checkpoint resolved — the user's answers, verbatim"*. **The history
of what was asked is kept below, unedited in substance, with each question's outcome appended** — a
resolved checkpoint should still show what the executor refused to decide for itself.

### Question 1 — does §4's Assumption A2 derivation convince you as written?

`14-RESEARCH.md`'s Assumptions Log flags A2 **High risk if wrong**. It is the entire justification for
abstain-then-sort and for the "controlled over the decided subset" claim. §4 states it as two
checkable lines: the rejection set is a deterministic function of the observed data, therefore
σ(data)-measurable, therefore it pulls out of the conditional expectation, so the identity holds for
*any* data-dependent selection — unlike the frequentist case, where pre-selection is a genuine
selective-inference problem. **Your call: does that read as a derivation a referee can check, or does
it need a referee-facing rewrite before it goes near the manuscript?**

> **ANSWERED — accept as written** (verbatim above). **Outcome: §4 is unchanged, byte for byte.** No
> rewrite was attempted; manuscript-level polish is deferred to Phase 16 by the user's instruction.

### Question 2 — are you content with how the five judgement-call bars are presented?

`skill ≥ 0.6`, `spearman ≥ 0.95`, `AUC_hard ≥ 0.8`, `ood margin ≥ 0.5`, and the `0.2` coverage floor.
All five were frozen in `a6c825867dbc786f7c3927df6760295ee0930c77` (a commit containing no Phase-14
result), you ruled them accept-as-proposed on 2026-08-03, and **none was relaxed** — including after
SC3-d missed. §6 labels all five JUDGEMENT CALLS WITH NO DERIVATION and records the licensing
standard. **They cannot be changed now — a bar changed after a result is a pre-registration breach —
but how they are *presented* in the manuscript is yours to rule.** Given this project has recorded
four cases of an underived bar measuring the wrong thing, are you content with the current framing?

> **ANSWERED — add a verdict-sensitivity strip** (verbatim above). **Outcome: Work Item A.** §6 gains
> a read-only strip computed by a new runner that gates nothing. **Every bar stayed frozen, every
> verdict stayed as-is, and `consts.jl` is byte-unchanged.** Presentation, not relaxation.

### Question 3 — how do you want the FAILED gate recorded?

SC3-d is NOT MET at **−0.129** against **0.5**, carried entirely by the `noise` family (0 of 1000
flagged) — the family the unwired `(:pp, :noise)` channel exists for. Nothing was tuned to try to
clear it. `P14_ITERATION_ALLOWANCE`'s single pre-declared trigger (conformal coverage *below* the
SC1-d band) did not fire and the allowance is unspent, so no in-phase remedy is available. **Say how
you want this recorded when you close the phase. No bar may be relaxed.**

> **ANSWERED SEPARATELY on 2026-08-04**, in
> `.planning/phases/14-decision-and-abstention-layer/14-SC3D-RULING.md` (commit `f4934a9`), which is
> the binding text and is not paraphrased here. **Outcome: Work Item B.** SC3-d is recorded
> **FAILED** and **scoped to the density-only wiring**, not to the shipped OOD flag; the fused
> detector's SC3-d margin is stated as **unmeasured**; `P14_ITERATION_ALLOWANCE` stays **unspent**;
> and the rejected alternatives are on the record with their reasons. **No bar was relaxed and
> nothing was re-wired.**

### Mechanical checks the user was asked to run — run by the executor, all green

Recorded twice: as measured at Task 1, and as re-measured after the rulings were applied.

| Check | At Task 1 | After the rulings |
|---|---|---|
| ten `test_p14_*.jl` individually | all exit 0, 1178 assertions | **all exit 0, 1179 assertions** (+1, cause recorded) |
| `spike/test/runtests.jl` used? | **no** | **no** |
| `git status --porcelain -- src spike/Project.toml spike/Manifest.toml corpus test` | empty | **empty** |
| report contains every honesty token | yes | **yes** |
| report contains the two withdrawn numerals | **no** | **no** |
| `spike/p14/consts.jl` byte-unchanged | yes | **yes** |

### Resume signal — CONSUMED

The user answered all three questions rather than typing "approved" bare: Question 1
*accept-as-written*, Question 2 *add a verdict-sensitivity strip*, Question 3 ruled separately in
`14-SC3D-RULING.md`. The two changes requested were applied (Work Items A and B, commit `6ad9fab`)
and the mechanical checks were re-run against the amended report. **No objection remains
unresolved, so the plan closes.**
