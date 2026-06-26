---
phase: 01-environment-smoke-gate
plan: 01
subsystem: infra
tags: [git, baseline, reproducibility, decoupling]

requires: []
provides:
  - "Frozen root baseline ref f581d95 for the untouched-root proof (ENV-01 / DEMO-02)"
  - "01-BASELINE.md recording the reconciliation decision and the verify command"
affects: [04-coupling-decoupling-proof, demo, env-01]

tech-stack:
  added: []
  patterns:
    - "Untouched-root proof asserts byte-identity against a recorded baseline commit ref"

key-files:
  created:
    - .planning/phases/01-environment-smoke-gate/01-BASELINE.md
  modified:
    - Project.toml
    - Manifest.toml

key-decisions:
  - "User (package author) chose commit-as-baseline: the GLMakie 0.10.5 [compat] pin + Manifest re-resolve were intentional and are frozen as the v2.0 spike baseline"
  - "Baseline ref = commit f581d95; later phases assert root-untouched against it"

patterns-established:
  - "Pattern: protected files = Project.toml, Manifest.toml, src/; .planning/ docs are not root drift"

requirements-completed: [ENV-01]

duration: 4min
completed: 2026-06-26
---

# Phase 1 / Plan 01: Root Baseline Reconciliation Summary

**Froze the intentionally-dirty root manifests (GLMakie 0.10.5 [compat] pin + full re-resolve) as the v2.0 spike baseline commit `f581d95`, with a recorded verify command for all later untouched-root proofs**

## Performance

- **Duration:** ~4 min
- **Tasks:** 2 (1 checkpoint:decision, 1 auto)
- **Files modified:** 3

## Accomplishments
- Surfaced the live `git diff` of the two protected root files to the user and captured an explicit reconciliation decision (`commit-as-baseline`).
- Committed `Project.toml` + `Manifest.toml` as `f581d95` — `chore(baseline): freeze root manifests as v2.0 spike baseline`.
- Recorded the baseline ref, decision, and the reusable `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` verify command in `01-BASELINE.md`.

## Task Commits

1. **Task 1: Decide the frozen baseline (checkpoint:decision)** — resolved interactively; user selected `commit-as-baseline`.
2. **Task 2: Execute reconciliation + record baseline** — root commit `f581d95`; `01-BASELINE.md` committed separately with this summary.

## Files Created/Modified
- `Project.toml` — committed the `[compat] GLMakie = "0.10.5"` pin (now part of baseline `f581d95`)
- `Manifest.toml` — committed the full re-resolve (now part of baseline `f581d95`)
- `.planning/phases/01-environment-smoke-gate/01-BASELINE.md` — baseline record + verify command

## Decisions Made
- **commit-as-baseline** (user choice): the dirty root was intentional; freezing it avoids discarding the deliberate GLMakie pin. Baseline ref is the resulting commit `f581d95`.

## Deviations from Plan
None — plan executed exactly as written for the chosen option. (Note: the plan's `<facts>` cited HEAD `9055a13` from planning time; HEAD had since advanced via planning commits, but the commit-as-baseline path and its `git diff --quiet HEAD` verify are HEAD-relative and unaffected.)

## Issues Encountered
- Git emitted LF→CRLF warnings on the Windows checkout when staging the manifests — cosmetic only, no content impact.

## User Setup Required
None.

## Next Phase Readiness
- Baseline `f581d95` is frozen; `src/` and root manifests are clean against it.
- Plan 01-02 (minimal `spike/` env + CPU-only NPE smoke) is unblocked.

---
*Phase: 01-environment-smoke-gate*
*Completed: 2026-06-26*
