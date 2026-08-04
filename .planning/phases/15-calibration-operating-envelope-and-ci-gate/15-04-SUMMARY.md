---
phase: 15-calibration-operating-envelope-and-ci-gate
plan: 04
subsystem: gate-decision-rules
tags: [pre-registration, ece, break-criterion, ood, sc2, source-guards, fixtures]
requires:
  - test/gate/p15_consts.jl (Tier-1 bars, from 15-01)
  - test/gate/sbc.jl (sbc_calibration / sbc_coverage / SBC_PARAM_LABELS — the production statistics)
  - test/gate/misspec.jl (gate_ood_roc's fire_rate and MEASURED id_fire_rate — consumed, not called)
  - src/amortized/ood.jl (the shipped detector — read-only, never assigned)
provides:
  - test/gate/p15_envelope.jl (the break criterion and the SC2 verdict as pure functions)
  - test/gate/test_p15_break.jl (129 assertions; both source guards shown to fire)
  - p15_ece_null as the bit-reproducible regenerator for P15_ECE_MARGIN
affects:
  - 15-07 (calls p15_axis_break / p15_report_only per rung; measures the rung-0 anchor)
  - 15-08 (calls p15_ood_rung / p15_axis_verdicts to produce the SC2 verdict)
  - 15-09 (reports the domain map; consumes the degenerate field and the mechanism labels)
tech-stack:
  added: []
  patterns:
    - pure decision reductions separated from the measurement that feeds them
    - source guards whose banned token is ASSEMBLED at load time so the guard cannot flag itself
    - structural bans checked on comment-stripped source (CONVENTIONS C-05), trailing comments included
    - branch fixtures placed at half-margin, never on the frontier (CONVENTIONS C-06)
key-files:
  created:
    - test/gate/p15_envelope.jl
    - test/gate/test_p15_break.jl
  modified: []
decisions:
  - "p15_ece_null defaults to Philox4x, not the plan's Xoshiro: only the Philox stream that p15_consts.jl RECORDS reproduces the pinned bars, and it does so exactly"
  - "The mce guard strips TRAILING comments as well as whole-line ones, because p15_consts.jl declares the prohibition in a trailing comment on a frozen const line"
  - "The detector guard is a zero-hit ban on decision code but VALUE EQUALITY on p15_consts.jl, which byte-forks the shipped operating point on purpose; equality is the stronger check"
metrics:
  duration_minutes: 34
  tasks: 3
  files_created: 2
  files_modified: 0
  tests_added: 129
  completed: 2026-08-04
---

# Phase 15 Plan 04: The Break Criterion and the SC2 Verdict Summary

Turned D-04a's anchor arithmetic and D-06's three-valued SC2 verdict into pure, self-asserting
functions with 129 fixtures, before any rung exists to argue with them — and found that the ECE
margin's regenerator reproduces the frozen Tier-1 bars **bit-for-bit**, not merely within
Monte-Carlo error.

## What Was Built

### `test/gate/p15_envelope.jl` (447 lines)

Pure reductions only. The file performs **no inference, no forward simulation and no net load** —
`grep -c 'posterior_for\|sampleposterior\|simulate_pair'` returns **0** on the raw source.

| Function | What it decides |
|---|---|
| `p15_ece(ranks_col, L)` | The closed-form ECE identity. Documented as the **null-model twin only**: the production path is `sbc_calibration(...).ece`, and a second implementation on the reported path would make the envelope incomparable to the shipped gate report. |
| `p15_ece_null(M, L; B, q, seed, rng)` | The Monte-Carlo null and therefore the **regenerator for `P15_ECE_MARGIN`** — no net, no simulation, no outcome. |
| `p15_read_shipped_anchor(path)` | Returns `nothing` rather than throwing when the git-ignored report is absent. Docstring states in full that this is **reporting context and explicitly NOT the break anchor**. |
| `p15_break_threshold_for(ece_rung0)` | Thin delegate to the Tier-1 `p15_break_threshold`, so no caller inlines `rung0 + margin`. |
| `p15_break_rung(eces, thr)` | First crossing, else `nothing` — D-13's UNBOUNDED reading, which is not a pass. |
| `p15_axis_break(rows, thr)` | Earliest crossing across the gating columns; asserts `P15_BREAK_COLUMNS == (1, 8)` at call time and rejects rung 0 as a ladder rung. |
| `p15_report_only(row)` | The non-gating record. The maximum-calibration-error field is absent **even as a pass-through**. |
| `p15_ood_rung(rates, id_fire_rate)` | Crossing against the **MEASURED** baseline plus the frozen binomial margin. |
| `p15_axis_verdict(L_break, L_ood)` | `:protective` / `:late` / `:silent_but_safe`, with `:late` named as the failure SC2 exists to catch. |
| `p15_sc2_pass(verdicts)` | Passes iff no axis is `:late`; docstring records D-14's two forbidden responses. |
| `p15_axis_verdicts(break, ood)` | Per-axis records carrying mechanism label and prior boundary rung, plus D-13's `degenerate` field. **Errors** rather than guessing when an axis is missing from either input — a missing key would otherwise read as "no break" and fake an unbounded envelope. |

### `test/gate/test_p15_break.jl` (469 lines, 129 assertions, 9 testsets)

Standalone-runnable and picked up automatically by `runtests.jl`'s Phase-15 block (15-01 wired the
fixed five-file list; **`runtests.jl` was not edited by this plan**). Loads `p15_envelope.jl` into
the isolated module `P15BreakEnv`, because `p15_consts.jl` defines the same const names the
template already loaded into the test module and a direct include would be silently skipped.

**The two source guards assemble their banned token at load time** (`"m" * "ce"`,
`"OOD_ID_" * "QUANTILE"`, `"id_" * "threshold"`) so the guard file does not flag itself — including
in its failure messages and in its own violating fixture. `grep -c 'mce' test/gate/test_p15_break.jl`
returns **0** on the raw source.

## Verification (observed, not assumed)

| Check | Command | Observed |
|---|---|---|
| Task 1 verify (import order corrected — Deviation 1) | the plan's one-liner | printed `ok` |
| Task 2 verify | the plan's one-liner, verbatim | printed `ok` |
| Task 3 standalone | `julia --project=. --threads=auto test/gate/test_p15_break.jl` | **129 pass, 0 fail, 0 error**, `EXIT=0` |
| Task 3 wall-clock | `time` around the same command | **18.8 s** total (testset body 6.6 s) — bar is 60 s |
| 15-01 testset still green | `julia --project=. --threads=auto test/gate/test_p15_consts.jl` | `EXIT=0` |
| Full suite | `julia --project=. --threads=auto -e 'using Pkg; Pkg.test()'` | `Testing ProteinCoLoc tests passed`; the new testset ran inside it at 129/129 |
| Envelope reads the forbidden field | `grep -v '^\s*#' … \| grep -c 'mce'` | `0` |
| Envelope performs inference | `grep -c 'posterior_for\|sampleposterior\|simulate_pair'` | `0` |
| Envelope retunes the detector | `grep -c 'OOD_ID_QUANTILE\s*='` | `0` |
| Read-only files untouched | `git status --porcelain src/amortized/ood.jl test/gate/sbc.jl test/gate/misspec.jl` | empty |
| Accidental deletions | `git diff --diff-filter=D --name-only HEAD~3 HEAD` | empty |

`test_p15_families.jl` could not be run: it is sibling plan 15-03's file and does not exist in this
worktree. The suite reported it as absent-and-skipped by name, which is the wired behaviour.

## The margin's reproducibility, on the record (required by the plan's `<output>`)

`p15_ece_null(1000, 999)` at the **recorded** design (`B = P15_ECE_NULL_B = 20 000`,
`seed = P15_ECE_NULL_SEED = 20260804`, `Philox4x(UInt64, (seed, 0))`, `L = SBC_L = 999`):

| Quantity | Regenerated | Tier-1 pin | Agreement |
|---|---|---|---|
| `E[ECE₀]` | `0.010328326315789473` | `P15_ECE_NULL_MEAN = 0.01033` | exact at the pin's 5 dp |
| `q95(ECE₀)` | `0.019263157894736843` | `P15_ECE_NULL_Q95 = 0.01926` | exact at the pin's 5 dp |
| `q95 − mean` | `0.008934831578947370` | `P15_ECE_MARGIN = 0.00893` | exact at the pin's 5 dp |
| closed form `0.336/√M` | `0.010625252938165755` | — | pin is 2.8 % below, as its comment states |

**This is a bit-reproduction, not a statistical agreement**, and it settles a question the plan left
open: the pinned bars were measured on the Philox stream the pre-registration records, so the margin
can be re-derived exactly by anyone with the artifacts tree deleted. Regeneration cost: **0.47 s**.

The cheaper `B = 2000` figures the testset actually uses, for reference:
mean `0.010370763`, q95 `0.019371053`, margin `0.009000289` — deviations from the pins of
`4.1e-5`, `1.1e-4` and `7.0e-5` against tolerances of `0.0015`, `0.0030` and `0.0015`
(≈ 6 SE, derived from sd(ECE₀) ≈ 0.0045 at M = 1000, not chosen).

Closed-form check at three M, all inside the 5 % bar: M = 250 → 3.0 %, M = 1000 → 2.4 %,
M = 2000 → 2.0 %.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Task 1's verify one-liner used `Random` before importing it**
- **Found during:** Task 1, first verify run
- **Issue:** The command reads `… r = rand(Random.Xoshiro(7), 0:999, 500); import Random; …`.
  `julia -e` evaluates top-level statements in order, so the observed failure was
  `ERROR: UndefVarError: `Random` not defined in `Main``. The envelope file itself had already
  loaded cleanly — the defect is in the verification command, not in the code under test.
- **Fix:** Moved `import Random` to the front of the same one-liner. **Nothing else was changed**,
  and every assertion in it ran and passed.
- **Commit:** n/a (verification command, not a repository file)

**2. [Rule 1 - Bug] `p15_ece_null` defaults to `Philox4x`, not the plan's `Random.Xoshiro`**
- **Found during:** Task 1
- **Issue:** The plan directs drawing "from a seeded `Random.Xoshiro(seed)`" while calling the
  function "the regenerator for the Tier-1 margin: the bar can be re-derived from the test design
  alone". Those two statements are incompatible. `p15_consts.jl:628` records the design as
  `Philox4x(UInt64, (seed, 0))`, and 15-01's Deviation 3 records that it *measured* the pins under
  exactly that design. A Xoshiro-based function is a second **estimate** of the same quantity, not a
  **regeneration** of the pinned one — the distinction this phase's licensing standard turns on.
- **Fix:** Default `rng = Philox4x(UInt64, (UInt64(seed), UInt64(0)))`, with the `rng` keyword kept
  so any stream can be supplied. The docstring states why. **Evidence the fix is right:** the
  regenerator now reproduces all three pinned bars exactly (table above); under Xoshiro it could
  only have agreed to Monte-Carlo error.
- **Files modified:** `test/gate/p15_envelope.jl`
- **Commit:** `b677297`

**3. [Rule 1 - Bug] The forbidden-field guard must strip TRAILING comments, not only whole-line ones**
- **Found during:** Task 3, pre-scanning the guard targets before writing the testset
- **Issue:** The plan specifies "DROP every line matching `^\s*#`, then assert no surviving line
  matches `\bmce\b`". `test/gate/p15_consts.jl:698` is
  `const P15_MCE_FORBIDDEN = true      # reading `.mce` anywhere in Phase 15 is a defect` — a
  **trailing** comment on a code line. A whole-line-only filter leaves it, so the guard would flag
  **the very file that declares the prohibition**, and that file is frozen, append-only and one of
  the eight SHA-256-pinned files. The guard would have been unsatisfiable without a
  pre-registration breach.
- **Fix:** `_p15brk_code_of(line)` returns everything up to the first `#` that is not inside a
  double-quoted string. This is the correct reading of CONVENTIONS C-05 in both directions — a ban a
  comment can *satisfy* is not a ban, and a ban a comment can *violate* is not a ban either. The
  firing fixture in testset 9 covers **both** the whole-line and the trailing case explicitly, so
  the filter is demonstrated rather than trusted. Its per-line string tracking is documented as a
  stated limitation.
- **Files modified:** `test/gate/test_p15_break.jl`
- **Commit:** `28f3db6`

**4. [Rule 3 - Blocking] The detector-retuning guard is split: zero-hit ban on decision code, VALUE EQUALITY on the pre-registration**
- **Found during:** Task 3, same pre-scan
- **Issue:** `test/gate/p15_consts.jl:268` carries `const OOD_ID_QUANTILE = 0.95` as part of its
  deliberate **byte-fork** of the amended-v2 rule set — 15-01's testset 4 asserts shared constants
  per-name equal to `gate_consts_8_v2.jl`. A zero-hit ban over that file is unsatisfiable without
  editing a frozen, pinned pre-registration, which is exactly the breach this phase's discipline
  exists to prevent. `misspec.jl:390`'s `id_threshold = maha_thr` shows the same regex also matches
  NamedTuple **field names**, which are not assignments to a detector parameter at all.
- **Fix:** Two-part guard. **(a)** Zero-hit ban over the decision code — `p15_misspec.jl`,
  `p15_envelope.jl`, `ci_golden.jl` and every existing `test_p15_*.jl`. **(b)** For
  `p15_consts.jl`, assert the recorded value **equals** the shipped `ProteinCoLoc.OOD_ID_QUANTILE`,
  that it is 0.95, that its only hits are `const` declarations, and that it assigns no threshold at
  all. **This is strictly stronger than the ban it replaces on that file:** absence of the line
  would not catch a *changed* number, and equality does.
- **Files modified:** `test/gate/test_p15_break.jl`
- **Commit:** `28f3db6`

**5. [Rule 3 - Blocking] Task 2's action text and its own acceptance criterion contradict each other**
- **Found during:** Task 2
- **Issue:** The action requires `p15_ood_rung`'s docstring to state "`OOD_ID_QUANTILE = 0.95`
  (`src/amortized/ood.jl:61`)", while the same task's acceptance criterion requires
  `grep -c 'OOD_ID_QUANTILE\s*=' test/gate/p15_envelope.jl` to return **0**. The docstring the
  action asks for is precisely the string the criterion forbids.
- **Fix:** Kept the substance and dropped the literal form: the docstring reads "The detector's
  in-distribution operating quantile is 0.95 (`src/amortized/ood.jl:61`), so about 5 % of
  IN-DISTRIBUTION images fire BY CONSTRUCTION". D-07's full reasoning — why "any image fires" and
  why a chosen fixed rate were both rejected — is present verbatim in substance. The executable
  criterion returns 0.
- **Files modified:** `test/gate/p15_envelope.jl`
- **Commit:** `26cb03f`

### Additions not named in the plan

Three, all additive and all asserted; none changes a value or a rule the plan specified.
`P15_ENVELOPE_REPO_ROOT` (so the default report path resolves from the file rather than from the
process working directory), the `rng` keyword on `p15_ece_null` (Deviation 2), and the
missing-axis error in `p15_axis_verdicts` (a `get(..., nothing)` fallback would have let an
incomplete sweep be read as an unbounded envelope, which is one of D-13's two pre-declared
conclusions — too consequential to reach by accident).

## Assumption Drift (advisory)

**A. The regenerator's RNG family is a free choice.** *Planned:* the plan specifies
`Random.Xoshiro(seed)` for `p15_ece_null` and sets purely statistical tolerances, implying the bar
is only re-derivable to Monte-Carlo error. *Actual:* the pinned bars are **bit-reproducible** under
the Philox stream `p15_consts.jl` records, and are so reproduced here. *Why it matters:* the ECE
margin is a stronger artifact than the plan assumed — 15-09 can state that the margin regenerates
exactly from the recorded design in 0.47 s, rather than that it agrees within error.

**B. "Comment filtering" means dropping whole-line comments.** *Planned:* `^\s*#` line dropping,
stated as sufficient for CONVENTIONS C-03/C-05 compliance. *Actual:* this repository declares its
prohibitions in **trailing** comments on `const` lines, so the whole-line reading makes a guard
flag the file that states the rule. *Why it matters:* 15-06 (the golden re-bless discipline) and
15-09 are likely to write source guards of the same shape; they should reuse `_p15brk_code_of`'s
treatment rather than re-derive the whole-line filter and hit the same wall.

**C. The shipped gate report is readable during Phase-15 work.** *Planned:* testset 3 is written as
though `artifacts/amended_v2/grid_8/gate_report_8.jld2` may be present. *Actual:* `artifacts/` is
git-ignored, so it is absent in **every worktree**, not merely in a fresh clone — the value
assertions in testset 3 have therefore never executed anywhere in this plan's run. This is a
tolerated, `@info`-announced skip, but 15-07 must not assume the anchor reader has been exercised
against a real file.

## Known Stubs

None.

## Threat Flags

None. This plan adds no network endpoint, no auth path and no schema at a trust boundary. Its only
file access is reading already-committed repository source inside a test plus one `mktempdir`
fixture. `T-15-17`, `T-15-18` and `T-15-19` are mitigated as planned: the threshold form, margin and
column scope stay in the frozen Tier-1 block with `P15_BREAK_COLUMNS == (1, 8)` asserted at call
time; both source guards exist and have been shown to fire; and `p15_report_only` does not pass the
artifact statistic through at all.

## Notes for Later Plans

- **`p15_ece` must never appear on the reported path.** The sweep computes ECE with
  `sbc_calibration(...).ece`. `p15_ece` exists to drive the null and to be cross-checked; using it
  in the sweep would make the envelope incomparable to the shipped gate report.
- **`p15_axis_break` rejects rung 0.** Rung 0 is the shared anchor, measured once per net; feeding
  it in as a ladder rung raises an `AssertionError` by design.
- **`p15_axis_verdicts` requires every one of `P15_AXES` in both inputs** and errors otherwise, so a
  partially-completed sweep cannot be silently reported as `:unbounded_envelope`.
- **Any new Phase-15 source file is scanned automatically** by both guards if it is named
  `test_p15_*.jl`, or is one of `p15_misspec.jl` / `p15_envelope.jl` / `ci_golden.jl`. Assembling a
  banned token from fragments is the pattern to copy if another guard is added.

## Self-Check: PASSED

- `test/gate/p15_envelope.jl` — FOUND (447 lines)
- `test/gate/test_p15_break.jl` — FOUND (469 lines; 100 source lines carry a `@test` macro, and the
  runner reported 129 assertions because several sit inside loops)
- `.planning/phases/15-calibration-operating-envelope-and-ci-gate/15-04-SUMMARY.md` — FOUND
- commits `b677297`, `26cb03f`, `28f3db6` — all present in `git log`
- `git status --porcelain` on `src/amortized/ood.jl`, `test/gate/sbc.jl`, `test/gate/misspec.jl` — empty
- `STATE.md` and `ROADMAP.md` — NOT modified by this plan (the orchestrator owns those writes)
