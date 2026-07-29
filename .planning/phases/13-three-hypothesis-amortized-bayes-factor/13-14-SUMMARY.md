---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 14
status: complete
subsystem: reporting / phase close-out
tags: [phase-report, amended-criteria, named-limits, iteration-ledger, scope-of-evidence, honesty]
requires:
  - "13-SC2-AMENDMENT.md — the amended criteria, frozen 2026-07-27 at sha df16e769, before any Phase-13 result existed"
  - "spike/p13/tau_probe_report.jld2 (13-10) — the whole A(delta) curve, the measured tau, all three fingerprints"
  - "spike/p13/three_way_net.jld2 (13-11) — head_log_odds, the recipe, the Phase-11 provenance, the training trajectory"
  - "spike/p13/three_way_gate_report.jld2 (13-12) — every gated and reported gate statistic, the thresholds, the mechanical trigger evaluation"
  - "spike/p13/alpha_report.jld2 (13-13) — the :simulated ladder, alpha_star, the eight invariants, the warning signs"
  - "spike/p13/realimage_report.jld2 (13-16) — the lambda sweep with paired OOD verdicts, both detectors, alpha_star_real, the read-only digest pair"
  - "spike/p13/consts.jl — the byte-locked pre-registration, READ ONLY"
  - "the fifteen preceding Phase-13 plan SUMMARYs"
provides:
  - ".planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md — the Phase-13 results report (1,086 lines): the Scope-of-evidence header, the amended criteria beside the originals, the D-05 trap, the D-07 derivation with its negative controls, the gate, the two alpha arms kept separate, the real-image arm with every number beside its OOD verdict, what was DECLINED, four named limits, the iteration ledger, the deferrals, all eight open questions closed, the decoupling evidence, and five recorded source discrepancies"
  - ".planning/ROADMAP.md — the Phase-13 checklist completed across all sixteen plan lines, with a Report pointer and a simulator-gate-only Verdict line"
affects:
  - "Phase 14 — inherits the F5 confident-tail precision caveat and the descriptive-not-decision-rule boundary"
  - "Phase 16 — inherits the sealed-holdout deferral, the PENDING-FETCH prerequisite and named limit A"
  - "a future productionization phase — inherits the log_bf_simplex rename and the four named limits in liftable shape"
tech-stack:
  added: []
  patterns:
    - "Report-time re-read: every quoted number was loaded from its .jld2 at write time rather than copied from a SUMMARY, which is what surfaced the deep-tail band discrepancy"
    - "Positional honesty control: the scope-of-evidence block is asserted to PRECEDE section 1, so a qualitative arm cannot be met without its framing"
key-files:
  created:
    - .planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md
  modified:
    - .planning/ROADMAP.md
decisions:
  - "One section beyond the plan's enumerated list was added — `## What was DECLINED` — placed between 8c and the named limits, because the orchestrator's honesty bar requires declines as first-class content and the plan's ordering is otherwise preserved exactly."
  - "The conditional FIFTH named limit (the prior-atom asymmetry) was NOT added, because its stated condition — that the gate show exclusion-end degradation — was evaluated against the artifact and did not hold (deep-tail exclusion AUC 0.997257 against 0.979025 nearest tau). It is recorded as an evaluated non-finding rather than silently dropped."
  - "docs/amortized.md was NOT edited; the four limits are written in its entry shape so a future productionization can lift them verbatim."
metrics:
  duration: ~55 min
  completed: 2026-07-29
  tasks: 2
  commits: 2
---

# Phase 13 Plan 14: The Phase-13 Report Summary

**PHASE VERDICT, IN ONE LINE: the three-way amortized Bayes factor works on simulated data — the
amended gate cleared all six pre-registered criteria on one run, on thresholds frozen before the net
existed, without spending the iteration allowance — and it is not shown to work on real microscopy,
because both real fixtures are flagged out-of-distribution at ~2.5x their own in-distribution
thresholds and the registration-aware basis did not fix that.**

Both halves of that sentence are in the report's opening, not one in the headline and one in a
limitations appendix.

## The `## Scope of evidence` block, reproduced verbatim from the report header

> **Scope of evidence.** The Phase-13 gate is simulator ground truth (`rho_true < -tau`, exact labels,
> well-powered, `P13_GATE_M = 4000`). Two further arms are reported and **neither is a gate**: the
> D-16 semi-synthetic alpha-graded random-to-exclusion series (`P13_ALPHA_GATED = false`), and a
> qualitative check on six committed real microscopy TIFFs (`P13_REAL_IS_GATED = false`,
> `P13_REAL_QUALITATIVE_ONLY = true`). The TIFFs carry no colocalization label — the folder names
> `positive`/`negative` are the original package's biological test conditions, and the "negative" pair
> in fact has mean patch correlation **+0.2481** (`alpha0_mbar` = +0.24805233, `realimage_report.jld2`).
> The real-data arm can show the three-way verdicts behave sensibly on real microscopy; **it cannot
> show they are correct.** It covers two specimens, and both of them are flagged out-of-distribution
> by both the Phase-13 and the shipped detector, so every real-image number in this report is printed
> beside its OOD verdict and its lambda. A real image has no known registration uncertainty, so the
> whole pre-registered lambda ladder is reported and the widest, most conservative rung
> (`lambda = 3.0`) is the headline; lambda was never estimated from the images. The sealed-holdout
> rows of the provenance manifest were **not** opened.

| Arm | Substrate | Standing |
|---|---|---|
| Simulator ground truth | `rho_true < -tau`, exact labels, well-powered | **GATES** |
| D-16 alpha-graded series | mask-based disjoint reassignment, simulated **and** real | Supporting evidence — **not a gate** |
| `test/test_images/` qualitative check | six committed 1028 × 1376 microscopy TIFFs | Required deliverable — **explicitly not a gate** |

*No pass/fail threshold is defined for any real-image quantity.*

The block sits in the report **HEADER**, before section 1 and therefore before the gate verdict, the
tau curve and every table. That placement is machine-asserted (`HEADER_OK`), because a reader who
meets the real-image section without having met this framing will read a behaviour check as a
validation.

## The four named-limit titles

- **Limit A — the physical segregation anchor is n = 0 in practice.** One anchor exists, it is triply
  unavailable (`sha256 = "PENDING-FETCH"`, `bytes = 0`, `split = sealed_holdout`), and Phase 13
  deliberately left the anti-snooping seal intact for Phase 16's blind evaluation.
- **Limit B — "exclusion" is defined as negative intensity correlation, not disjoint localization.**
  `alpha_star = 0.25` calibrates the gap for one constructed transform; alpha is a construction
  parameter, not a physical quantity.
- **Limit C — this net is research-lane and inherits no gate lineage from the shipped binary net.**
  It trains on the Phase-11 basis, lives outside `RatioEstimator`, was validated on its own terms, and
  its persisted artifact is an epoch-4 checkpoint of an untuned early-overfitting run. Nothing here
  changes the shipped artifact or the Phase-7 GO.
- **Limit D — the exclusion hypothesis has no labelled real-data validation.** Two specimens, no
  colocalization label, both out-of-distribution on both detectors (417.30 / 167.54 = 2.49x; 433.69 /
  179.14 = 2.42x), both channel pairs include the nuclear counterstain, no real exclusion example
  exists anywhere in this phase, and the labelled check is deferred to Phase 16.

A conditional **fifth** entry (the ~2.5 : 1 prior-atom asymmetry against the exclusion end) was
**evaluated and NOT added**: its condition was "if the gate showed exclusion-end degradation", and the
gate showed the opposite. Recorded in the report as an evaluated non-finding rather than dropped.

## Iteration-ledger outcome

**`P13_ITERATION_ALLOWANCE = 1`. NOT SPENT. It remains 1 of 1, UNSPENT.**

The one pre-declared trigger was evaluated **mechanically**, on the number it is worded about, and
`iteration_trigger_fired = false` is persisted in `three_way_gate_report.jld2`:

| Quantity | Measured (artifact) | Floor | Trigger half satisfied? |
|---|---|---|---|
| exclusion AUC, deep tail `\|rho\| > 0.9` (n = 239) | **0.9972566363348159** | 0.90 | **no** — far above the floor |
| exclusion AUC, mid range `\|rho\| <= 0.9` (n = 1095) | 0.986299314554666 | 0.90 | yes (resolved) |

Both halves must hold; the first does not. `P13_STRATIFICATION_FALLBACK = :importance_weighted` was
not switched on, nothing was retrained, and `spike/p13/consts.jl` is byte-unchanged. Two clarifications
the report states explicitly: the allowance is **scoped to the gate** and could not have been spent on
a real-image observation (no real-image quantity has a bar to miss); and the discipline is about
**pre-registration integrity, not compute budget** — training took 2.37 minutes, so a retrain is cheap,
and that is precisely not the reason it is forbidden.

## Decoupling evidence — commands run, outputs observed

| Command | Observed |
|---|---|
| `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | **exit 0 — byte-unchanged** (checked before Task 1's commit and again after Task 2's) |
| `git diff --quiet HEAD -- test/ docs/ corpus/` | **exit 0 — byte-unchanged**; `docs/amortized.md` was not edited |
| `git status --porcelain test/test_images/` | empty |
| `git log --oneline -- spike/p13/consts.jl` | exactly two commits: `c42cc8e` (13-01 Tier-1) and `bfba6ac` (13-10 Tier-2 tau append). **No commit that produced a Phase-13 result touched it.** |
| `git rev-parse HEAD:spike/p13/consts.jl` | `70fe66df33832af2e648ed3f81cc3ad356ce4cd8` (tracked git blob, CRLF→LF filtered) |
| `git hash-object --no-filters spike/p13/consts.jl` | `5a4ea2223ce11bcae2f974c643ec1e523458e8ac` — the value every runner recorded as `consts_git_blob_sha` |
| `sha256(read("spike/p13/consts.jl"))` | `100e97a37bb470bb7cb4fbdbd8e33f219d83adcf61fe398a90cf3ad3391f3fa6` — matches `consts_sha` / `net_consts_sha` in every artifact |
| `julia --project=. -e 'using Pkg; Pkg.test()'` | **exit 0 — `Testing ProteinCoLoc tests passed`** (last testset `amended gate machinery (07-GATE-AMENDMENT required code changes)` 254/254) |
| `git diff --diff-filter=D --name-only HEAD~1 HEAD` (both commits) | empty — no deletions |
| `git diff -U0 -- .planning/ROADMAP.md \| grep '^[+-]### Phase'` | none — no other phase block was touched |

**Not satisfiable, and not this plan's to fix:** the plan's acceptance criterion
`git diff --quiet HEAD -- src/ test/ docs/ spike/ corpus/` exits **1**, because a **concurrently
executing Phase-12 agent** has `spike/test/test_p12_architecture.jl` modified in the shared working
tree. This plan touched no file under `spike/`; both of its commits used explicit paths
(`git commit -- <path>`), so nothing of that agent's work was staged or committed here. The hard
constraint that actually binds Phase 13 — `src/`, `spike/Project.toml`, `spike/Manifest.toml` — holds
at exit 0.

## Verification performed (every command was RUN; the observed result is what is reported)

| Check | Observed |
|---|---|
| Task-1 `<verify>` (six `occursin` assertions on the report) | prints **`REPORT_OK`** |
| Task-1 positional check (scope-of-evidence precedes section 1) | prints **`HEADER_OK`** |
| `grep -c 'GATES'` / `grep -c 'not a gate'` | **1** / **2** (bars: ≥1, ≥2) |
| `grep -c '^## '` | **17** (bar: ≥14) |
| `grep -c 'Limit A\|Limit B\|Limit C\|Limit D'` | **4** (bar: ≥4) |
| `grep -c '433.69\|433.7'` / `'0.2481'` / `'biological'` | **3** / **5** / **3** (bars: ≥1 each) |
| `grep -c 'qualitative'` / `'Phase 16'` | **7** / **8** (bars: ≥1, ≥2) |
| Pitfall 2d — `grep -ci 'validated on real\|validated against test_images'` | **0** |
| No gate language on a real-image line — `grep -cE '(real-image\|test_images).*\b(PASS\|FAIL)\b'` | **0** |
| `grep -c 'Open questions closed'` / `'less colocalized than the control'` / `'P13_ITERATION_ALLOWANCE'` | **1** / **2** / **2** |
| Report length | **1,086 lines** (bar: ≥220) |
| Artifact key-links present — `three_way_gate_report` / `realimage_report` | **5** / **5** occurrences |
| Task-2 `<verify>` (Phase-13 block, sixteen plan lines, Report + Verdict + amendment + TBD present) | prints **`ROADMAP_OK`** |
| `awk '/^### Phase 13:/,/^### Phase 14:/' \| grep -c '13-[0-9][0-9]-PLAN.md'` | **16** |
| same pipeline for `13-15-PLAN.md\|13-16-PLAN.md` | **2** |
| unticked lines with no explicit reason | **0** |
| `Verdict:` line containing a real-image number | **0** (caught and fixed once — see deviation 2) |
| Phase-13-scoped `grep -c "reproduces \`compute_BayesFactor()\`"` | **1** (the original SC2 wording is still present) |

**Spot-checks of quoted numbers against their artifacts** (all read at report time, not copied):
`auc_coloc = 0.9901687725007021` and `auc_exclusion = 0.9882624329245729`
(`three_way_gate_report.jld2`); `alpha_star = 0.25` (`alpha_report.jld2`); headline real
`log BF(C:R) = -0.465771` / `log BF(E:R) = -10.272712` at `lambda = 3.0`, `ood_density = 417.29737`,
`ood_threshold = 167.54463`, `is_ood = true` (`realimage_report.jld2`).

**What this plan did NOT run.** It wrote two Markdown files and modified no executable, so no Phase-13
test file needed re-running and none was. The **full spike suite was not re-run**: it exits 1 at the
Phase-4 `SPEEDUP_GATE` (`test_npe.jl:230`, measured 68.35 against the pre-registered bar 100.0),
masking the entire Phase-13 include block. That is pre-existing, logged as `deferred-items.md` D-13-A,
and re-deriving or relaxing a pre-registered threshold is a user decision. The plan's
`<verification>` item 7 ("`spike/test/runtests.jl` still exits 0") is therefore **unsatisfiable for
reasons that predate this plan**, and saying so is more useful than a run that would tell us what
three prior plans already observed.

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 2 — Missing critical] One section beyond the plan's enumerated list: `## What was
DECLINED`**
- **Found during:** Task 1, drafting.
- **Issue:** the plan enumerates sections 1–13 and the header, but the dispatch brief's honesty bar
  requires the phase's *declines* to be first-class content, not scattered asides. Ten of them exist
  (the KDE reference, continuity as a criterion, MCE as the gate statistic, the confusion matrix as a
  decision rule, any real-data confusion matrix, the unspent iteration, a **third** amendment, the
  `src/` rename, the corpus anchor and CBS, and a third channel pair) and each was a deliberate act
  rather than an omission. Left scattered, a reader would read them as gaps.
- **Fix:** one extra section, placed between 8c and the named limits, so the plan's own ordering and
  numbering are preserved exactly.
- **Files:** `13-REPORT.md`. **Commit:** `0c1a0bd`.

**2. [Rule 1 — Bug] The first `Verdict:` line tripped the plan's own real-image-number assertion**
- **Found during:** Task 2 verification.
- **Issue:** the draft verdict line ended "the α-graded series and the `test/test_images/` check are
  reported and are not gates", which matches the acceptance criterion
  `grep 'Verdict:' | grep -c '0.2481\|433.69\|test_images'` → must be `0`. The criterion is right and
  the line was wrong: a roadmap verdict line should not point at the real-image fixture directory at
  all.
- **Fix:** rephrased to "the real-image qualitative check". The criterion now outputs `0`.
- **Files:** `.planning/ROADMAP.md`. **Commit:** `48fb808`.

**3. [Rule 2 — Missing critical] The conditional fifth named limit was evaluated, not assumed**
- **Found during:** Task 1, §9.
- **Issue:** the plan says to add the prior-atom note as a fifth entry "if the gate showed
  exclusion-end degradation". Writing it unconditionally would have asserted a degradation that was
  not measured; omitting it silently would have hidden that the question was asked.
- **Fix:** the condition was evaluated against the artifact (deep-tail exclusion AUC **0.997257** vs
  **0.979025** in the band nearest tau — the opposite of degradation), the entry was not added, and
  the evaluation is recorded in §9 as an evaluated non-finding.
- **Files:** `13-REPORT.md`. **Commit:** `0c1a0bd`.

### Recorded, not deviations

**4. `.planning/STATE.md` was written ADDITIVELY by hand, and no mutating state handler was run.**
`## Current Position` is a hand-authored block carrying three concurrent phases (11 closed, 12
planned, 13 executing). `bm-sdk query state.advance-plan` / `state.update-progress` rewrite that block
from disk-derived counts and would have flattened the Phase-11 and Phase-12 lines while a Phase-12
agent is live on this branch. Only the Phase-13 lines and the `last_activity` scalar were touched;
Phase 11's and Phase 12's lines are byte-unchanged. `roadmap.update-plan-progress` was likewise not
run, because the plan scopes Task 2 to the Phase-13 block only.

**5. No `requirements.mark-complete` was run.** The plan's `requirements:` frontmatter lists D-IDs
(context decisions), not REQ-IDs, and ROADMAP Phase 13 records `**Requirements**: TBD`. There is
nothing in `REQUIREMENTS.md` to tick.

**6. A parallel agent is executing Phase 12 on this branch.** Both commits used explicit paths; the
ROADMAP was re-read immediately before editing (its last touching commit was still `2ce97a8`, so no
race occurred), and nothing was rebased or amended.

## Assumption Drift (advisory)

**A. A number the briefing and a SUMMARY both assert is not what the artifact says.**
- **Planned:** `13-12-SUMMARY.md` §5 and the dispatch brief both state the deep-tail exclusion AUC
  0.997257 is "the **best** of the four |ρ| bands", and use that phrasing to argue the trigger did not
  fire.
- **Actual:** `three_way_gate_report.jld2` records
  `band_auc = [0.9790253685870447, 0.9957089272318079, 0.9981442729103328, 0.9972566363348159]`. The
  highest band is **(0.75, 0.90] at 0.998144**; the deep tail is the **second** highest.
- **Why:** discrimination rises monotonically across the first three bands and dips by 0.0009 in the
  fourth — a difference far below anything this design can resolve, but "the best of the four" is
  still not what was measured.
- **Materiality:** nothing substantive moves — the deep tail is 0.097 above its floor, the weakest
  band is the one nearest tau, and `iteration_trigger_fired = false` either way. But this is exactly
  the class of drift the "read every number from the artifact" rule exists to catch, and the report
  states the ranking as measured and records the correction in §14.

**B. The briefing's framing of the two `consts.jl` fingerprints is inverted.**
- **Planned:** the brief describes `70fe66df…` as "the git blob" and `5a4ea222…` as "the runners' own
  **normalized** fingerprint".
- **Actual:** it is the other way round. `70fe66df…` is the git blob **because git normalizes** CRLF→LF
  before hashing; `5a4ea222…` is `git hash-object --no-filters`, i.e. the hash of the **raw,
  un-normalized** working-tree bytes, and that is what the runners record.
- **Why:** on Windows with `core.autocrlf` active the two always differ, and neither is wrong.
- **Materiality:** both were verified at report time to describe the same file (sha256
  `100e97a3…`), so the integrity claim is unaffected — but a reader comparing them without knowing
  which is which would conclude the pre-registration had changed. §14.2 states it.

**C. `deferred-items.md`'s `runtests.jl:169` line reference has drifted to 189.**
- **Planned:** D-13-A (and three SUMMARYs) cite the abort at `spike/test/runtests.jl:169`.
- **Actual:** `test_npe.jl` is now included at **line 189**, because the concurrent Phase-12 agent
  wired its aggregator in first (`261a0cc`, `afae172`). The Phase-13 include block has moved to lines
  220–239.
- **Why:** a shared file under concurrent edit.
- **Materiality:** the abort (`test_npe.jl:230`, `SPEEDUP_GATE`) and its masking consequence are
  unchanged; only the citation is stale. Recorded in §14.3 so a future reader looking at line 169
  does not conclude the record is wrong.

## Known Stubs

None. Every table in the report is populated from a persisted artifact or from the frozen
`consts.jl`; no placeholder, mock or hardcoded value flows into any number. No new executable was
written, so no data source is left unwired.

## Threat Flags

None. This plan wrote two Markdown files. It adds no network endpoint, no auth path, no file-access
pattern and no schema at a trust boundary. **No package was installed.** The register's own threats
were handled as designed: T-13-41 (a number drifting from its artifact) fired for real and is
recorded above as drift A; T-13-10 (reporting only the amended criterion) is mitigated by §1's
side-by-side table and the Phase-13-scoped grep asserting the original SC2 wording still present;
T-13-42, T-13-54, T-13-55, T-13-57 and T-13-26 are mitigated by the machine-checked criteria in the
verification table; T-13-05 stays `accept (deferred to Phase 16)` — `git diff --quiet HEAD -- corpus/`
exits 0 and nothing in this plan reached the sealed corpus.

## Requirements satisfied

**D-12** (the amended criteria stated with the original ROADMAP wording verbatim beside them and the
amendment cited), **D-15** (the Scope-of-evidence block in the HEADER before any result, with the
three-arm standing table and the greppable no-threshold sentence; the real-image arm reported with
its OOD verdicts, its naming correction and its qualitative-only statement; the four named limits
written), **D-05** (the correctness trap named in prose, with its executable regression assertion
cited), **D-07** (the derivation, why the binary formula does not generalize, and why no scalar can
repair a within-class reshape — with both negative controls and the missed `maxabs` bar reported as
measured), **D-08** (the `log_bf_simplex` misnomer recorded as a deferred one-line productionization
item with `src/` proven byte-unchanged), **D-13/D-03** (each head's AUC beside its ECE and its
empty-bin count; the lambda response reported honestly, flat and monotone), **D-04** (the iteration
ledger, explicit and mechanical), **D-01** (research lane; the shipped artifact, the shipped
`bayes_factor` accessor and the Phase-7 GO all recorded untouched), **D-02** (the different input
surface stated as the reason continuity is not a criterion), **D-16** (both alpha arms reported
separately and never averaged), **D-15/D-16** (the ROADMAP tick-off covers all sixteen plan lines
including 13-15 and 13-16).

## Commits

| Commit | Task | Description |
|---|---|---|
| `0c1a0bd` | Task 1 | the Phase-13 report — the gate passed 6/6 and both real fixtures are OOD-flagged |
| `48fb808` | Task 2 | tick the Phase-13 checklist across all sixteen plans and record the gate verdict |
| _(this)_ | close-out | this SUMMARY plus the additive STATE.md write |

## Self-Check: PASSED

- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md` — FOUND (1,086 lines)
- `.planning/ROADMAP.md` — FOUND, Phase-13 block edited only (`git diff -U0` shows no other
  `### Phase` heading in its hunks), 4 insertions / 2 deletions
- commit `0c1a0bd` — FOUND in history
- commit `48fb808` — FOUND in history
- `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` — exit 0
- `git diff --quiet HEAD -- test/ docs/ corpus/` — exit 0
- `git diff --quiet HEAD -- spike/p13/consts.jl` — exit 0, byte-unchanged
- post-commit deletion scan over both commits — no deletions
