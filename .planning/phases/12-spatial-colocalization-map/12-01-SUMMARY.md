---
phase: 12-spatial-colocalization-map
plan: 01
status: complete
subsystem: pre-registration + test wiring
tags: [tier-1, seeds, thresholds, include-ordering, r-5]
requires: []
provides:
  - "spike/validation/p12_consts.jl — Tier-1 pre-registration (seeds, salts, counters, ladders, gate thresholds, two repo roots, descope routing, ITERATION_ALLOWANCE)"
  - "spike/test/test_p12_consts.jl — literal re-assertion + executable seed disjointness + wiring manifest"
  - "spike/test/test_p12_suite.jl — the single Phase-12 include point, with printed per-include markers"
  - "nine green pending scaffolds, one per later Phase-12 plan"
affects:
  - "spike/test/runtests.jl (one added include line + comments; nothing removed or reordered)"
tech-stack:
  added: []
  patterns:
    - "two-tier append-only pre-registration (p11_consts.jl idiom), extended to THREE Tier-2 sentinels"
    - "isolated read-only module (GateV2P12) to READ the frozen gate constants, never retype them"
    - "printed P12-SUITE-RAN markers as the only detector of a skipped include"
key-files:
  created:
    - spike/validation/p12_consts.jl
    - spike/test/test_p12_consts.jl
    - spike/test/test_p12_suite.jl
    - spike/test/test_p12_lattice.jl
    - spike/test/test_p12_prior.jl
    - spike/test/test_p12_architecture.jl
    - spike/test/test_p12_datagen.jl
    - spike/test/test_p12_train.jl
    - spike/test/test_p12_result.jl
    - spike/test/test_p12_sbc.jl
    - spike/test/test_p12_coverage.jl
    - spike/test/test_p12_decoupling.jl
  modified:
    - spike/test/runtests.jl
decisions:
  - "The forbidden-seed LIST carries P13_BURNED_DEV_SEEDS once and PROVES that the three Phase-11 burned tuples are subsumed by it, rather than concatenating all four — which would have double-counted twelve keys and made the plan's own no-duplicate assertion unsatisfiable."
  - "The aggregator is wired FIRST among the phase-test includes, not merely before the include that throws today. 'Before the throwing file' is contingent on knowing which file throws; first position is structural. Orchestrator ruling, 2026-07-29."
  - "Testset 9 asserts the aggregator precedes EVERY sibling include, with sibling names derived from the source text rather than hardcoded, so the assertion cannot go stale the way the original premise did."
  - "The Phase-4 SPEEDUP_GATE shortfall is NOT touched: re-deriving that bar is a pre-registration decision the user owns (standing project ruling)."
metrics:
  duration: ~85 min
  completed: 2026-07-29
  tasks: 3
  commits: 5
  files: 13
---

# Phase 12 Plan 01: Tier-1 Pre-registration and Phase-12 Test Wiring — Summary

Every Phase-12 threshold, seed, salt and counter is now frozen as a literal before a single
Phase-12 number exists, with the gated/reported split machine-checkable, the DEV stream
executably proven disjoint from all 43 reserved, burned and recomputed streams in the repo, and
all ten Phase-12 testsets provably *executing* under `spike/test/runtests.jl` rather than merely
wired into it.

## The wiring blocker, and how it was resolved

Executed as three tasks, then a fourth commit resolving a blocker this plan's own verify exposed.
Kept in full below because the failure mode is instructive and the fix changed the plan's
mechanism.

### What disagreed with what

The plan (and `12-PATTERNS.md` §0.3) both state that `test_p13_correction.jl` is the include
whose thrown testset "aborts the remaining includes", and derive from that the instruction to
insert the Phase-12 aggregator *immediately before it* so that Phase-12 testsets are guaranteed
to run.

**That premise is false as of 2026-07-29.** A strictly earlier include already throws:

```
runtests.jl:169  include(joinpath(@__DIR__, "test_npe.jl"))
  SC3 (NPE-03): Test Failed at spike/test/test_npe.jl:230
    Expression: sr.median_speedup > SPEEDUP_GATE
     Evaluated: 84.43856733198933 > 100.0
  ERROR: LoadError: Some tests did not pass: 149 passed, 1 failed, 0 errored, 0 broken.
  in expression starting at spike/test/runtests.jl:169
```

Because that throw aborts every later include, **nothing from `runtests.jl:174` onward runs
today** — not the Phase-5 validation bundle, not the Phase-9 comparator, not any of the eleven
Phase-13 includes, and not the Phase-12 aggregator this plan wires in at line 219.

Measured twice, on two independent runs (84.44 and 89.77 against a gate of 100.0), the second
with no competing processes — so it is not a load artefact. The `>100×` speedup gate genuinely
misses on this machine.

### Evidence (the failing state, before the fix)

`Task 3 <verify>` run verbatim against the tree as the plan specified it:

```
SUITE_EXIT=1
MISSING Phase-12 include: test_p12_consts.jl
MISSING Phase-12 include: test_p12_lattice.jl
MISSING Phase-12 include: test_p12_prior.jl
MISSING Phase-12 include: test_p12_architecture.jl
MISSING Phase-12 include: test_p12_datagen.jl
MISSING Phase-12 include: test_p12_train.jl
MISSING Phase-12 include: test_p12_result.jl
MISSING Phase-12 include: test_p12_sbc.jl
MISSING Phase-12 include: test_p12_coverage.jl
MISSING Phase-12 include: test_p12_decoupling.jl
VERIFY_RESULT: MISSING=1 (0 == pass)
```

All ten markers are absent. Run standalone (`include("spike/test/test_p12_suite.jl")`), the same
aggregator prints all ten markers and exits 0 with zero failures — so the aggregator is correct;
only its *position in the suite* fails to deliver execution.

The plan's `must_haves.truths` entry **"A Phase-12 testset that is added later cannot silently
never run"** was therefore NOT delivered by the specified insertion point. It was delivered
against `test_p13_correction.jl` and against nothing else. *(It is delivered now — see the
resolution below.)*

### Why the executor escalated instead of auto-fixing

The mechanical fix looked like a one-line move: put `include(joinpath(@__DIR__, "test_p12_suite.jl"))`
before `include(joinpath(@__DIR__, "test_npe.jl"))` at `runtests.jl:169` instead. That would have
satisfied every stated acceptance criterion (one added include line, nothing removed or reordered,
lower line index than `test_p13_correction.jl`, testset 9 still passes) *and* the verify.

It was **not** applied by the executor, for two reasons — and the escalation was upheld, since the
ruling both hardened the fix beyond that one-line version and supplied a standing project ruling
the executor did not have:

1. It hops Phase-12's testsets in front of a pre-existing, unaddressed abort that also kills
   Phases 5, 9 and 13. Doing that for Phase 12 only is a cross-phase process decision of the same
   class `runtests.jl:191-199` explicitly refuses to make in passing ("Resolving that shortfall is
   a phase-level scientific decision ... and is NOT resolved by this wiring").
2. A Phase-13 executor is live on this branch and its entire test surface is currently dead in
   the suite. That is a finding the orchestrator needs before anyone reshuffles include order.

### The three options, for the orchestrator

| Option | Action | Cost |
|---|---|---|
| A | Move the one aggregator line before `runtests.jl:169` | one line; Phase-12 tests run today; Phases 5/9/13 stay dead |
| B | Fix or re-scope the Phase-4 `SPEEDUP_GATE` (currently 84.44 vs 100.0) | the whole suite runs again; this is a v2.0 headline-claim question, not a wiring one |
| C | Accept that Phase-12 testsets are run via `test_p12_suite.jl` directly, not via `runtests.jl` | zero code; the "cannot silently never run" guarantee is downgraded to a convention |

### RESOLUTION — orchestrator ruling, 2026-07-29: **Option A, hardened** (commit `afae172`)

Two facts decided it, neither of which the executor had:

1. **The placement was deliberate-but-stale, not deliberate-and-load-bearing.** `.planning/STATE.md`
   already recorded this exact finding, measured during 13-11 on the same day: *"the suite never
   REACHES line 213. It aborts EARLIER, at runtests.jl:169 -> test_npe.jl:230, on the Phase-4
   SPEEDUP_GATE (median speedup 68.35 against the pre-registered bar 100.0)."* The project had
   corrected its model of the suite **in writing, before 12-01 was dispatched**; the plan's
   insertion point simply predated that correction. Fixing it therefore serves the plan's stated
   intent rather than overriding a live design choice.
2. **Option B is out of scope by standing ruling.** STATE.md: *"Re-deriving that bar is a
   pre-registration decision and is left for the user."* `SPEEDUP_GATE` and `test_npe.jl` were
   not touched.

Option A was taken **hardened**, because plain Option A repeats the defect at one remove:

**(a) The aggregator is now FIRST among the phase-test includes**, ahead of `test_simulator.jl`,
not merely ahead of whichever file happens to throw. "Before the throwing include" is contingent
on knowing which file throws *today* — that is precisely the construction that just failed. First
position is *structural*: no other phase's failure can mask Phase 12, whichever file throws. The
comment block moved with it and now names the real abort (Phase-4 `SPEEDUP_GATE`, pre-existing,
user-owned) instead of `test_p13_correction.jl`. The same stale claim was also corrected in
`test_p12_suite.jl`'s own header, which had asserted the old placement in prose.

**(b) Testset 9 now encodes the guarantee instead of the stale premise.** The old assertion
(aggregator offset < `test_p13_correction.jl` offset) stayed *true* after the move but was
insufficient for the same reason it was insufficient before — it was true throughout the period
the aggregator never ran. It is replaced by an assertion that the aggregator precedes **every**
other phase-test include, with the sibling filenames **derived from the source text** by regex
rather than hardcoded, since a literal list would go stale exactly as the original premise did.
The `test_p13_correction.jl` check is retained as an explicitly-labelled corollary because it is
a stated acceptance criterion of this plan.

**Proof the strengthening bites:** reverting the aggregator to its original placement in a scratch
copy — a placement that *passed* the old assertion — now fails the new one six times over, naming
each include that jumped ahead: `test_simulator.jl`, `test_data_pipeline.jl`, `test_npe.jl`,
`test_sbc.jl`, `test_bf.jl`, `test_ood.jl` …

**Task 3 `<verify>`, re-run against the committed tree:**

```
SUITE_EXIT=1
VERIFY_RESULT: MISSING=0 (0 == pass)
  6:P12-SUITE-RAN: test_p12_consts.jl        35:P12-SUITE-RAN: test_p12_train.jl
 19:P12-SUITE-RAN: test_p12_lattice.jl       39:P12-SUITE-RAN: test_p12_result.jl
 23:P12-SUITE-RAN: test_p12_prior.jl         43:P12-SUITE-RAN: test_p12_sbc.jl
 27:P12-SUITE-RAN: test_p12_architecture.jl  47:P12-SUITE-RAN: test_p12_coverage.jl
 31:P12-SUITE-RAN: test_p12_datagen.jl       51:P12-SUITE-RAN: test_p12_decoupling.jl
```

**Both facts recorded, because only one of them is ours.** All ten markers are present and every
Phase-12 testset reports green (143 + 9 = 152 passes, 0 failures). **The suite still exits
non-zero**, at `test_npe.jl` (`runtests.jl:189` after the insertion shifted line numbers), for a
Phase-4 reason Phase 12 does not own and did not touch.

**Net effect on `runtests.jl` versus its pre-Phase-12 state: 20 insertions, 0 deletions** — one
`include(` line plus its comment block. No Phase-13 include was touched or reordered.

## What Was Built

**Task 1 — `spike/validation/p12_consts.jl`** (commit `ee86ab8`, 801 lines)

Tier 1 in the `p11_consts.jl` shape: the two-tier header naming **three** reserved Tier-2
sentinels with the plan that opens each; an isolated `module GateV2P12` that READS the frozen
imsize mixture and both derived ship-gate seed families rather than retyping them; the forbidden
seed set; the salt inventory; the fresh DEV/FIXTURE streams and their nine reserved reported
counters; lattice geometry; the bound `P12_R1_PRIOR`; SIM-02; the D-12 Stage-1 and Stage-2 gating
thresholds with their arithmetic recorded verbatim; the MCAR masking set; three reporting-only
bars; the budget; the iteration allowance; two repo roots; the shared `_p12_interval` helper;
executable descope routing; five declared deviations; the gated-vs-reported tuples; and the
Stage-1 golden-regression triple.

**Task 2 — `spike/test/test_p12_consts.jl`** (commit `2d620b6`, 242 lines, 10 testsets, 122 tests)

**Task 3 — the aggregator, nine scaffolds and the wiring line** (commit `261a0cc`, 11 files)

## Verification Results — what was actually run and observed

| Check | Command | Observed |
|---|---|---|
| Task 1 verify | the plan's `julia -e 'include(...)'` one-liner | **PASS** — printed `72 1 tier1-ok` |
| `spike/Project.toml` / `Manifest.toml` unchanged | `git diff --quiet HEAD --` | **PASS** |
| `src/` and `corpus/` unchanged | `git diff --quiet HEAD -- src corpus` | **PASS** |
| non-comment `P12_` line count > 60 | `grep -v '^\s*#' \| grep -c` | **PASS** — 138 |
| all 17 required literals present | `grep -F` each | **PASS** — 17/17 |
| `P12_R1_PRIOR` is a bound `Uniform` | `params(...) == (0.05, 0.95)` | **PASS** |
| double include, no redefinition | included twice in one session | **PASS** |
| Tier-2 absent (3 assertions) | in-file + testset 8 | **PASS** |
| `P12_PRIMARY_CHECKOUT != P12_REPO_ROOT` inside a worktree | loaded from a throwaway `git worktree add --detach` tree | **PASS** — worktree root vs `.../ProteinCoLoc`; equal in the primary checkout |
| git-absent fallback warns and names the degradation | `withenv("PATH" => "")` | **PASS** — warning emitted, enclosing root returned |
| `p12_require_proceed` throws with no verdict file | temp dir | **PASS** (also exercised PROCEED, DESCOPE, and malformed-file throw) |
| Task 2 verify | `include("spike/test/test_p12_consts.jl")` | **PASS** — 143/143 after the ruling (122 before), 0 failed, 0 errored, exit 0 |
| mutation `0.95 → 0.94` breaks testset 5 | scratch copy | **PASS** — 2 failures in "recorded derivations reproduce" |
| moving the aggregator after the throwing arm breaks testset 9 | scratch copy | **PASS** — `4887 < 4833` failed |
| **strengthened** testset 9 rejects the *original* placement | scratch copy | **PASS** — 6 failures naming `test_simulator.jl`, `test_data_pipeline.jl`, `test_npe.jl`, `test_sbc.jl`, `test_bf.jl`, `test_ood.jl` |
| each scaffold: one non-comment `P12_PENDING_SCAFFOLD`, no `@test false`/`@test_broken` | `grep` | **PASS** — 9/9 |
| `runtests.jl` net diff vs pre-Phase-12 | `git diff 2d620b6 --` | **PASS** — 20 insertions (1 include + comments), 0 deletions |
| aggregator standalone | `include("spike/test/test_p12_suite.jl")` | **PASS** — all ten markers, exit 0 |
| **Task 3 verify (markers in the full suite log)** | the plan's loop over all ten | **FAIL before the ruling** (all ten missing) → **PASS after** (`MISSING=0`, all ten present, every Phase-12 testset green) |
| suite exit code | `julia --project=spike spike/test/runtests.jl` | **non-zero, expected and not ours** — aborts at `test_npe.jl` on the Phase-4 `SPEEDUP_GATE`, unchanged from the pre-plan baseline |

**Not run:** `julia --project=. -e 'using Pkg; Pkg.test()'`. No pre-plan baseline for it was
captured, so "unchanged from its pre-plan state" is not something this run can assert. `src/` is
byte-unchanged and this plan adds nothing to the main package, so the main-package suite cannot
have been affected by it. Recorded rather than implied.

## Deviations from Plan

### 1. [Rule 3 — Blocking] The forbidden-seed list cannot contain all four burned tuples *and* satisfy the no-duplicate assertion

- **Found during:** Task 1.
- **Issue:** The plan asks for `F2_DEV_SEEDS`, `SPIKE_DEV_SEEDS`, `P11_RESEARCH_BURNED` **and**
  the 25-entry `P13_BURNED_DEV_SEEDS`, and separately for
  `@assert length(unique(_p12_forbidden_list())) == length(_p12_forbidden_list())`. All twelve
  entries of the first three tuples are *also* entries of `P13_BURNED_DEV_SEEDS` — Phase 13
  re-collected the same historical burns. A naive concatenation therefore duplicates twelve keys
  and the assertion fires at include time, so the file would not load at all.
- **Fix:** All four constants are declared exactly as specified, read from their defining sites.
  `_p12_forbidden_list()` lists the *superset* once, and three new `issubset` assertions PROVE
  that nothing is lost by doing so. The forbidden SET is identical under either construction;
  only the double-count is removed, and the duplicate check stays live and meaningful.
  The plan text specifies which constants to declare, not the body of `_p12_forbidden_list()`, so
  this is a reading that satisfies both instructions rather than a change to either.
- **Commit:** `ee86ab8`.

### 2. Extra in-file self-checks and extra literal assertions beyond the plan's enumerated list

- **Found during:** Tasks 1 and 2.
- **What:** `p12_consts.jl` carries roughly thirty foot-of-block `@assert`s where the plan's item
  (d) lists five; `test_p12_consts.jl` asserts the *inputs* to each recorded derivation as
  literals inside testset 5 and the three reporting-only bars inside testset 6.
- **Why:** `p11_consts.jl` — the file R-5 says to mirror — carries about twenty such self-checks,
  so an enumerated subset reads as a floor, not a ceiling; and the test file's own purpose header
  claims "every locked value is asserted AS A LITERAL", which would over-claim if only the four
  seeds were. Asserting the derivation inputs also stops testset 5 passing by both sides drifting.
- **Consequence worth knowing:** mutating `P12_STAGE1_RATIO_CEILING` now fails at *include* time
  (in-file assert) as well as in testset 5. Demonstrated both ways in a scratch copy.
- **No testsets were added or removed**; the plan's list of ten is exactly what shipped.

### 3. Task 2's verify was run after Task 3's files existed

- **Why:** testset 9 reads `test_p12_suite.jl` and `runtests.jl`, both of which Task 3 creates —
  so Task 2 cannot verify green before Task 3 exists. Files were written for both tasks, both
  verifies were run, then the two commits were made in plan order. Inherent to the plan's task
  ordering, not a change to it.

## Assumption Drift (advisory)

**The suite's throwing include is not the one the phase modelled.** Planned: `runtests.jl`'s only
abort-causing include is `test_p13_correction.jl`, so "before it" means "guaranteed to run".
Actual: `test_npe.jl` at `runtests.jl:169` throws first and has done so for some time, so
"before `test_p13_correction.jl`" and "guaranteed to run" are not the same position any more.
Recorded here as well as in the blocker because it changes how a reader should interpret every
"the Phase-12 testsets run in the same harness" claim in this phase's documents.

**A second, smaller drift.** `12-PATTERNS.md` §0.3 pins the insertion point at "between line 207
and line 208". Those line numbers had already drifted to 218/219 by execution time. The plan
anticipated this and told me to locate by content, which is what was done.

## Known Stubs

Nine deliberate pending scaffolds — `test_p12_{lattice,prior,architecture,datagen,train,result,
sbc,coverage,decoupling}.jl`. Each carries the literal `P12_PENDING_SCAFFOLD` in a comment and in
an `@info`, asserts only the CPU-only invariant, and contains no failing assertion **by design**:
a red scaffold would throw and abort the remaining Phase-12 includes, which is precisely the trap
this wiring exists to close. Owners: lattice → 12-03, architecture → 12-04, decoupling → 12-05,
prior → 12-07, datagen → 12-09, result → 12-12, train → 12-14, coverage → 12-16, sbc → 12-18.

## Threat Flags

| Flag | File | Description |
|------|------|-------------|
| threat_flag: availability | `spike/validation/p12_consts.jl` | `p12_primary_checkout()` falls back **only** when `Sys.which("git") === nothing`, as the plan specifies. Observed consequence: loading the file from a directory tree that is not a git checkout (an export, a tarball, a sandbox) makes git exit 128 and the include **throw**, taking every Phase-12 runner with it. Encountered for real while building the scratch mutation harness. The `test_p13_result.jl:199-201` precedent guards on "is this a checkout at all" for exactly this reason. Tier 1 is append-only, so this is **reported, not repaired**. |

## Carried-forward advisories — accepted as reported, deliberately UNREPAIRED

Both are consequences of Tier 1 being append-only: they are now documented facts rather than
open work items. Orchestrator ruling 2026-07-29 accepts both as reported and carries them upward.

1. **`p12_primary_checkout()` throws outside a git checkout.** The threat-flag row above, in full.
   The plan specifies fallback *only* on `Sys.which("git") === nothing`, so a tree where git
   exists but the directory is not a checkout produces exit 128 and a thrown include rather than a
   degraded-but-working root. Hit for real, not hypothesised.

2. **`Z_TWO_SIDED_90` is an unprefixed name that `spike/validation/p11_consts.jl:209` already
   binds to `1.644853627`; Phase 12 binds it to `1.644854`** — more and fewer digits of the same
   quantile Φ⁻¹(0.95). It is safe today for two independent reasons, both checked rather than
   assumed: the two consts files are never loaded into one module by the test suite (nothing
   `runtests.jl` reaches includes `p11_consts.jl` into `Main` — `p13/preconditions.jl` reads it
   through the isolated `module _P11C` precisely to avoid the sentinel collision), and under
   Julia 1.12 a same-name `const` rebinding is legal rather than an error. Were both ever loaded
   into one module, the later binding would silently change the default `z` of
   `p11_stats.jl`'s `wilson_ci`. Recorded in `p12_consts.jl` §7 as a fact; **both files are
   frozen pre-registrations and neither may be edited**, so no repair is available or attempted.

## For the Next Plan

- Every Tier-1 name later plans consume is defined and verified. `P12_CACHE_ROOT`/`P12_POOL_SCHEMA`
  (12-09), `P12_THETA_ROWS`/`P12_ZSCORE_ARM` (12-04), `P12_SPATIAL_SENTINEL` (12-12) and the
  `P12_*_LOAD_ONLY` flags are owned elsewhere, as the plan specifies. Nothing is left undefined.
- `P12_GOLDEN_IMSIZE`/`P12_GOLDEN_KEYS` are deliberately **absent** from Tier 1; 12-08 owns them
  locally, following `capture_p11_golden.jl:51-52`.
- Tier 1 is now **append-only**. The three reserved sentinels are `:P12_CHOSEN_PRIOR` and
  `:P12_N_LOW` (12-15) and `:P12_FISHERZ_NEFF` (12-16). Each append removes its own negative
  assertion in `p12_consts.jl` **and** its mirror in testset 8 of `test_p12_consts.jl`, in the
  same commit, and only then.
- **Add a Phase-12 testset by adding a line to `test_p12_suite.jl` and to `P12_TEST_FILES` in
  `test_p12_consts.jl` — never to `runtests.jl`.** Testset 9 fails if a file is missing from
  either. The aggregator is first among the phase-test includes and testset 9 asserts it stays
  there, so a new Phase-12 test cannot be masked by another phase's failure.
- **The suite still exits non-zero** at `test_npe.jl` (Phase-4 `SPEEDUP_GATE`, 84.44 vs 100.0).
  Phase-12 testsets all run and report *before* that point. Do not read a red suite as a Phase-12
  failure; read the `P12-SUITE-RAN` markers and the Phase-12 testset summaries.

## Self-Check: PASSED

All 13 created/modified files verified present on disk; all five commits verified present in
`git log`: `ee86ab8` (Tier 1), `2d620b6` (literal-assertion test), `261a0cc` (aggregator,
scaffolds, wiring), `afae172` (orchestrator ruling: aggregator first + structural assertion),
`f32a998`/this update (summary). `spike/data/cache/p11` intact at 54 MB; `src/`, `corpus/`,
`spike/Project.toml` and `spike/Manifest.toml` byte-unchanged; `artifacts/amended_v2/grid_8`
untouched; `.planning/STATE.md` and `.planning/ROADMAP.md` untouched; `test_npe.jl` and
`SPEEDUP_GATE` untouched.
