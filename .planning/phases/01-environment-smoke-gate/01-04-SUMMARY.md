---
phase: 01-environment-smoke-gate
plan: 04
subsystem: infra
tags: [julia, pkg, decoupling, neuralestimators, flux, sbi]

# Dependency graph
requires:
  - phase: 01-01
    provides: agreed root baseline ref f581d95 + untouched-root verify command
  - phase: 01-02
    provides: green CPU-only NeuralEstimators v0.2.1 smoke (PosteriorEstimator/NormalisingFlow)
  - phase: 01-03
    provides: committed known-good minimal spike Manifest (safe revert point) + NOTES.md
provides:
  - "Parent-package coupling resolved via documented include() fallback (D-01)"
  - "Verbatim evidence that Pkg.develop(path='..') downgrades NeuralEstimators 0.2.1 => 0.1.4 and breaks the smoke"
  - "Decoupling proof: root Project.toml/Manifest.toml/src byte-identical to baseline f581d95"
  - "ENV-01 satisfied; the DEMO-02 decoupling-proof pattern recorded for Phase 6"
affects: [phase-02-simulator, phase-04-advi-baseline, phase-06-go-no-go-memo]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "include() fallback over Pkg.develop when the parent dep tree conflicts with the pinned SBI stack"
    - "Decoupling proof via git diff --quiet against a fixed baseline ref (DEMO-02 pattern)"

key-files:
  created: []
  modified:
    - spike/NOTES.md

key-decisions:
  - "include() fallback chosen over Pkg.develop: the parent's GLMakie/Makie 0.21 + Turing 0.44 + GraphNeuralNetworks 1.1.0 tree caps NeuralEstimators at <0.2.1, downgrading the pinned 0.2.1 to 0.1.4 and breaking the v0.2.1-API smoke (D-01)"
  - "Package boundary deferred to Phase 4 (D-02): the ADVI baseline will use a dedicated parent-aware env or a targeted include of bayes.jl, not by dragging the parent into the lean SBI env"

patterns-established:
  - "Pattern: attempt clean package boundary first, revert on smoke-break, document conflict + fall back to include() — root never edited"
  - "Pattern: assert publication-integrity via git diff --quiet <baseline-ref> -- Project.toml Manifest.toml src/"

requirements-completed: [ENV-01]

# Metrics
duration: 13min
completed: 2026-06-26
---

# Phase 1 Plan 04: Parent-Package Coupling + Decoupling Proof Summary

**ENV-01 satisfied via the D-01 include() fallback after Pkg.develop downgraded NeuralEstimators 0.2.1 to 0.1.4 and broke the smoke; root Project.toml/Manifest.toml/src proven byte-identical to baseline f581d95.**

## Performance

- **Duration:** 13 min
- **Started:** 2026-06-26T15:11:43Z
- **Completed:** 2026-06-26T15:25:18Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments
- Attempted the clean package boundary (`Pkg.develop(path="..")` + `Pkg.resolve()`) and captured the verbatim conflict: it silently downgraded **NeuralEstimators v0.2.1 ⇒ v0.1.4** to satisfy the parent's heavy tree (Turing 0.44.5, GLMakie 0.10.18/Makie 0.21, Images 0.26, GraphNeuralNetworks 1.1.0, PackageCompiler 2.4.0).
- Confirmed the downgrade **broke the green smoke** (`MethodError: no method matching NormalisingFlow(::Int64; num_summaries::Int64)` — v0.1.4 ships the old ApproximateDistributions API), then reverted `spike/Project.toml` + `spike/Manifest.toml` to the Plan 01-03 known-good state and re-verified the smoke **GREEN** (4/4 pass, `mu_hat ≈ 0.799` vs `theta_true = 0.7`).
- Adopted and documented the **D-01 `include()` fallback** in `spike/NOTES.md`: later phases reach `src/colocalization.jl` (correlation/patch/_prepare_data) etc. via `include("../src/<file>.jl")` read-only; Phase 2 adds StatsBase+Statistics (lightweight) for that, Phase 4 decides the heavier Turing-model path.
- Proved the **decoupling guarantee**: `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` exits 0 — the frozen publication state is byte-identical to the agreed Plan-01 baseline; `src/` has zero modifications.

## Task Commits

Each task was committed atomically:

1. **Task 1: Attempt Pkg.develop coupling; fall back to include() on conflict** - `c2717d9` (docs)
2. **Task 2: Prove the decoupling guarantee (root + src/ untouched vs baseline)** - `b4ae100` (docs)

## Files Created/Modified
- `spike/NOTES.md` - Added §4 (Parent-Package Coupling Outcome — verbatim develop→downgrade→smoke-break evidence + include() fallback decision and Phase-2/4 wiring) and §5 (Decoupling Proof against baseline f581d95).

## Decisions Made
- **include() fallback over Pkg.develop (D-01):** The parent package transitively constrains NeuralEstimators to `<v0.2.1`. Forcing the develop would abandon the pinned, proven v0.2.1 stack — unacceptable — so the package boundary is deferred and `src/` is reached read-only via `include()`.
- **Package boundary deferred to Phase 4 (D-02):** Phase 4's ADVI baseline (which needs the real Turing `colocalization()` model) will use a dedicated parent-aware environment or a targeted include of `bayes.jl` with a minimal Turing-only dep set, rather than dragging GLMakie/Turing into the lean SBI env.

## Deviations from Plan

None - plan executed exactly as written. The `include()` fallback is the pre-authorized D-01 branch, not a deviation; the develop attempt legitimately conflicted (downgrade + smoke break), which is the documented trigger for the fallback.

## Issues Encountered
- `Pkg.develop(path="..")` did **not** raise an exception (so a naive "did it throw?" check would have falsely reported success) — it silently downgraded NeuralEstimators and only the smoke re-run surfaced the break. Resolution: treat "smoke must stay green" as the real success criterion (per the plan), capture the downgrade + `MethodError` as the verbatim conflict evidence, revert the two spike files via `git checkout --`, and adopt the include() fallback.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- ENV-01 complete: parent reachable read-only via the documented include() fallback; smoke green; decoupling proven. Phase 1 (environment + smoke gate) is fully satisfied across all four plans.
- **Phase 2 note:** to use the summary-statistic contract, add `StatsBase` + `Statistics` to the spike env (lightweight, no GLMakie/Turing) and `include("../src/colocalization.jl")`.
- **Phase 4 note:** the Turing ADVI baseline must NOT pull the parent into the lean SBI env (it caps NeuralEstimators <0.2.1); use a separate parent-aware env or a targeted `bayes.jl` include with a minimal Turing dep set.
- **Phase 6 note:** repeat the decoupling proof `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` for DEMO-02.

## Self-Check: PASSED

- FOUND: `.planning/phases/01-environment-smoke-gate/01-04-SUMMARY.md`
- FOUND: `spike/NOTES.md`
- FOUND commit: `c2717d9` (Task 1 — coupling outcome)
- FOUND commit: `b4ae100` (Task 2 — decoupling proof)
- FOUND commit: `5898b9d` (SUMMARY)

---
*Phase: 01-environment-smoke-gate*
*Completed: 2026-06-26*
