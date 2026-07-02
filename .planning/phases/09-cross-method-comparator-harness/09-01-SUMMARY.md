---
phase: 09-cross-method-comparator-harness
plan: 01
subsystem: testing
tags: [julia, dataframes, csv, comparator, pre-registration, resolve-risk-gate, costes, traffic-light]

# Dependency graph
requires:
  - phase: 03-training-data-pipeline
    provides: Random123/Philox keyed-salt seeding convention (HOLDOUT_SALT/FOLD_SALT) the new COSTES_SALT stays disjoint from
  - phase: 04-npe-training-advi-benchmark-ablation
    provides: the test_npe.jl pre-registered-const (D-14) idiom and the runtests.jl resolve-risk gate shape (d/e/f) copied for gate (g)
provides:
  - DataFrames + CSV promoted to direct spike deps with zero version drift (NeuralEstimators still v0.2.1)
  - spike/comparator/config.jl — nine pre-declared D-14 comparator thresholds committed before any table exists
  - spike/test/test_comparator.jl — CMP-01..08 skipped-scaffold testset wired into the single spike gate
  - runtests.jl resolve-risk gate (g) asserting DataFrames/CSV present, PythonCall/CondaPkg absent, NeuralEstimators pinned
affects: [09-cross-method-comparator-harness later waves, tapqir-bridge, comparator-table-artifact]

# Tech tracking
tech-stack:
  added: [DataFrames v1.8.2 (promoted to direct), CSV v0.10.16 (promoted to direct)]
  patterns:
    - "Pre-registration (D-14): comparator thresholds committed as const before any table run"
    - "Resolve-risk gate keyed by NeuralEstimators UUID re-asserted after every dep promotion"
    - "Skipped-scaffold testset enumerating every SC in the single spike gate"

key-files:
  created:
    - spike/comparator/config.jl
    - spike/test/test_comparator.jl
  modified:
    - spike/Project.toml
    - spike/Manifest.toml
    - spike/test/runtests.jl
    - .gitignore

key-decisions:
  - "DataFrames/CSV promotion is inert: both already resolved transitively at 1.8.2/0.10.16; Pkg.add added no packages, only flipped them to [deps]"
  - "COSTES_SALT = 0xA5A5A5A5DEADBEEF chosen distinct from HOLDOUT_SALT/FOLD_SALT so the Costes null stream never collides with data/CV streams"
  - "TAPQIR_PUBLISHED_VALUE = NaN documents 'not yet captured' rather than pinning a guessed number; plan 09-04 pins the real value"
  - "Gate (g) forbids PythonCall/CondaPkg in the MAIN env; the Tapqir Python stack is isolated to spike/comparator/tapqir_env only (T-09-03 boundary)"

patterns-established:
  - "Inertness check: snapshot Manifest to .manifest_baseline.toml, assert post-promotion versions equal the snapshot"
  - "Every new spike dep is guarded by a resolve-risk gate that re-asserts the NeuralEstimators v0.2.1 pin"

requirements-completed: [CMP-04, CMP-09]

# Metrics
duration: 12min
completed: 2026-07-02
---

# Phase 9 Plan 01: Comparator Wave-0 Foundation Summary

**DataFrames/CSV promoted to direct spike deps (inert, no NeuralEstimators downgrade), nine D-14 comparator thresholds pre-declared in config.jl, and a green CMP-01..08 skipped-scaffold testset wired into the single spike gate behind a Python-forbidding resolve-risk gate.**

## Performance

- **Duration:** ~12 min
- **Completed:** 2026-07-02
- **Tasks:** 3
- **Files modified:** 6 (2 created, 4 modified)

## Accomplishments
- Promoted DataFrames v1.8.2 + CSV v0.10.16 to `spike/Project.toml [deps]` with a proven-inert re-resolve (no version drift, NeuralEstimators stays pinned v0.2.1); root Project.toml/Manifest byte-unchanged.
- Committed `spike/comparator/config.jl` with all nine pre-declared D-14 thresholds (Costes randomization params, divergence traffic-light bands, harness master seed, Tapqir anchor placeholders) BEFORE any comparison table exists — the anti-data-snooping guarantee.
- Stood up `spike/test/test_comparator.jl` (CMP-01..08 named skipped child testsets, one per requirement) and wired it into `runtests.jl` after `test_npe.jl`; the whole suite runs green (exit 0) with the comparator scaffold skipped.
- Added resolve-risk gate (g) to `runtests.jl`: asserts DataFrames + CSV present, PythonCall + CondaPkg absent from the main env, and re-asserts NeuralEstimators == v0.2.1.

## Task Commits

Each task was committed atomically:

1. **Task 1: Promote DataFrames + CSV to direct deps and re-freeze the Manifest** - `7b97369` (chore)
2. **Task 2: Pre-declare the comparator constants in config.jl (D-14)** - `46c6a13` (feat)
3. **Task 3: Scaffold test_comparator.jl and wire the resolve-risk gate (g)** - `19f21d4` (test)

## Files Created/Modified
- `spike/Project.toml` - added DataFrames + CSV to `[deps]`
- `spike/Manifest.toml` - re-frozen (only project_hash changed; no version deltas)
- `spike/comparator/config.jl` - nine pre-declared D-14 comparator constants with rationale comments
- `spike/test/test_comparator.jl` - CMP-01..08 skipped-scaffold testset including config.jl consts
- `spike/test/runtests.jl` - include of test_comparator.jl + resolve-risk gate (g)
- `.gitignore` - ignore future `spike/comparator/tapqir_env/.CondaPkg/` and the `spike/.manifest_baseline.toml` inertness scratch

## Decisions Made
- The promotion is inert by design: DataFrames/CSV were already resolved transitively, so `Pkg.add` added zero packages and only marked them direct. The inertness check (baseline snapshot + version equality) proves no GLMakie-class downgrade occurred.
- `COSTES_SALT` set to a value provably distinct from `HOLDOUT_SALT`/`FOLD_SALT` (asserted in the config verify), keeping the Costes scramble stream disjoint from the data and CV RNG streams (D-10).
- `TAPQIR_PUBLISHED_VALUE` left as `NaN` with a comment pointing to plan 09-04 for pinning — honest "not yet captured" instead of a fabricated anchor.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None. Git reported the usual LF→CRLF warnings on the Windows checkout (cosmetic, no content impact).

## Threat Model Compliance
- **T-09-01 (Manifest resolve tampering):** mitigated — inertness check + gate (g) re-assert NeuralEstimators v0.2.1 after promotion.
- **T-09-02 (root project tampering):** mitigated — all `Pkg.add` scoped to `--project=spike`; `git diff` confirms root Project.toml/Manifest byte-unchanged.
- **T-09-03 (PythonCall/CondaPkg leaking into main env):** mitigated — gate (g) asserts both absent from the main spike env.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- The comparator scaffold, pre-registered thresholds, and tabular deps are in place; later Phase-9 waves can implement the classical estimators (09-02+), the table/artifact, the run_comparator entry point, and the Tapqir bridge (09-04, which pins TAPQIR_PUBLISHED_VALUE) by filling the corresponding CMP-* skipped testsets.
- No blockers.

## Self-Check: PASSED

---
*Phase: 09-cross-method-comparator-harness*
*Completed: 2026-07-02*
