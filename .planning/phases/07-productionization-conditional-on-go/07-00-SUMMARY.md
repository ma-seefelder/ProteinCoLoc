---
phase: 07-productionization-conditional-on-go
plan: 00
subsystem: packaging + result-type hierarchy
status: BLOCKED (co-resolution hard gate failed — user decision required)
tags: [breaking-release, package-extension, dependency-surgery, type-hierarchy, co-resolution-gate]
requirements: [PROD-01]
dependency_graph:
  requires: []
  provides:
    - AbstractColocResult hierarchy + shared accessor interface (delta_rho/bayes_factor/is_ood/posterior_draws)
    - AmortizedColocResult shipped subtype
    - ext/ProteinCoLocTuringExt.jl (internal weakdep-gated ADVI reference path)
    - root Project.toml v2.0.0 with hard amortized deps + Turing/CUDA weakdeps
  affects:
    - all downstream Phase-7 plans (dispatch on AbstractColocResult; BLOCKED on the co-resolution gate)
tech_stack:
  added: [NeuralEstimators (dep, 0.2.1), Flux (dep, 0.16.10), JLD2 (dep), HypothesisTests (dep), CUDA (weakdep), Turing (moved dep -> weakdep)]
  patterns: [Julia package extension (weakdep-gated), generic-function stubs + abstract-supertype dispatch, D-02 result-type hierarchy]
key_files:
  created:
    - src/results.jl
    - ext/ProteinCoLocTuringExt.jl
    - .planning/phases/07-productionization-conditional-on-go/07-00-SUMMARY.md
  modified:
    - Project.toml
    - src/ProteinCoLoc.jl
    - src/plot.jl
    - src/utils.jl
    - src/bayes.jl
    - test/runtests.jl
decisions:
  - D-02 hierarchy: AbstractColocResult supertype + accessor interface + AmortizedColocResult; AdviColocResult internal-only (ext), not exported
  - D-03 surgery: Turing/KDE/QuadGK-requiring ADVI path moved verbatim into ext/ProteinCoLocTuringExt.jl; core precompiles without Turing
  - CO-RESOLUTION GATE FAILED on a GLMakie/Makie conflict (NOT the Finding-1 NeuralEstimators downgrade, which is cleared)
metrics:
  tasks_completed: 2
  tasks_total: 3
  files_created: 3
  files_modified: 6
  completed_date: 2026-07-03
---

# Phase 7 Plan 00: Breaking-Release Foundation + Co-Resolution Gate Summary

**One-liner:** Established the D-02 `AbstractColocResult` type hierarchy and moved the entire
Turing/ADVI reference path into a weakdep package extension (clearing the Finding-1
NeuralEstimators-0.1.4 downgrade), but the root `Pkg.resolve()` HARD GATE **failed on a
secondary GLMakie(0.10.5)→Makie-0.21 vs NeuralEstimators-0.2.1→Makie-0.24 conflict** that
requires a user decision before any downstream Phase-7 inference plan can run.

## Status: BLOCKED

Tasks 1 and 2 are complete and committed. Task 3 (the co-resolution HARD GATE) FAILED and, per
the gate protocol, execution STOPPED without forcing or working around the conflict. **The plan
is NOT marked complete.** STATE.md carries the full `## Blocker` record; ROADMAP/plan-counter
were deliberately NOT advanced.

## What Was Built

### Task 1 — D-02 result-type hierarchy (`src/results.jl`) — COMPLETE (commit 2b98ee7)
- `abstract type AbstractColocResult` + the shared accessor interface (`delta_rho`,
  `bayes_factor`, `is_ood`, `posterior_draws`) defined on the supertype with an `_iface_error`
  fallback so unimplemented subtypes fail loudly.
- Shipped concrete subtype `AmortizedColocResult <: AbstractColocResult`
  (`grid`, `posterior::Matrix{Float64}` 7×N row-1=ρ_true, `delta_rho_draws`,
  `log_bayes_factor`, `ood::OODVerdict`, `calibration::CalibrationMeta`, `meta::NamedTuple`)
  with all four accessors.
- Supporting value types `OODVerdict` and `CalibrationMeta` (the spike's `CalibrationResult`
  fields + grid/gate provenance).
- Phase-12 `SpatialColocResult`/`delta_rho_map`/`uncertainty_map` and Phase-13
  `ThreeHypothesisColocResult`/`bayes_factor_simplex` present ONLY as commented extension-point
  sketches (not implemented).
- Exported the D-02 surface; `AdviColocResult` deliberately NOT exported (D-01).

### Task 2 — Dependency surgery → v2.0.0 (commit c832119) — COMPLETE
- `Project.toml`: added hard deps NeuralEstimators (compat `0.2.1`), Flux (compat `0.16.10`),
  JLD2, HypothesisTests; moved **Turing → `[weakdeps]`**; added **CUDA → `[weakdeps]`**
  (optional GPU training, D-06); added `[extensions] ProteinCoLocTuringExt = "Turing"`; bumped
  `version` 1.0.1 → **2.0.0**.
- `ext/ProteinCoLocTuringExt.jl` (new): internal, weakdep-gated ADVI reference — the demoted
  `AdviColocResult <: ProteinCoLoc.AbstractColocResult` (fields verbatim from `CoLocResult`),
  the Turing `@model` + ADVI `colocalization()`, the KDE/quadgk `compute_BayesFactor()`,
  `_prepare_data()`, `convert_posterior_samples()`, the D-02 accessors for the ADVI result,
  and the three ADVI-result plotters. Not exported.
- `src/ProteinCoLoc.jl`: dropped `using Turing`/`Variational`, `import KernelDensity: kde`,
  `import QuadGK: quadgk`, and the `include("bayes.jl")`; declared generic-function stubs for
  the ADVI-path entry points; removed `colocalization`/`compute_BayesFactor`/`CoLocResult`
  from the export list.
- `src/bayes.jl`: gutted to a MOVED-marker (content relocated to the extension).

### Task 3 — Co-resolution HARD GATE — FAILED (see Blocker)
- Ran an actual root `Pkg.resolve()` (root project activated, spike never touched).
- Added the resolve-assertion gate + D-02 unit tests to `test/runtests.jl`; retired the
  now-invalid ADVI-path public-API tests (that path is internal reference, validated by the
  later per-grid D-05 ship-gate).

## CO-RESOLUTION GATE RESULT: **FAILED**

```
Unsatisfiable requirements detected for package Makie [ee78f7c6]:
 Makie [ee78f7c6] log:
 ├─possible versions are: 0.9.0 - 0.24.12 or uninstalled
 ├─restricted to versions 0.21.18 by an explicit requirement, leaving only versions: 0.21.18
 └─restricted by compatibility requirements with NeuralEstimators [38f6df31] to versions: 0.24.0 - 0.24.12 or uninstalled — no versions left
   └─NeuralEstimators [38f6df31] log:
     ├─possible versions are: 0.1.0 - 0.2.1 or uninstalled
     └─restricted to versions 0.2.1 - 0.2 by ProteinCoLoc [12345678], leaving only versions: 0.2.1
       └─ProteinCoLoc [12345678] is fixed to version 2.0.0
```

**Interpretation:** The Turing→extension surgery succeeded at its stated purpose —
NeuralEstimators is pinned to **0.2.1** and is NOT downgraded to 0.1.4 (Finding 1 cleared). The
remaining blocker is a **different** package: `GLMakie = "0.10.5"` (root compat) forces
**Makie 0.21.18**, while **NeuralEstimators 0.2.1 requires Makie 0.24.x** — incompatible.

**Candidate fixes (07-00 Task 3 checkpoint — user decision required):**
1. `glmakie-to-ext` — move GLMakie into a weakdep extension too (removes Makie from the core
   resolve). Cost: core plotting (`plot`/`local_correlation_plot`/`plot_mask`) becomes
   weakdep-gated; the Plots test set moves behind the extension.
2. Bump GLMakie — drop the `0.10.5` compat pin to a GLMakie that ships Makie 0.24.x. Cost:
   possible GLMakie API touch-ups in `src/plot.jl`.
3. Escalate a harder blocker.

## Deviations from Plan

### [Rule 3 - Blocking] Turing-path code in `src/plot.jl` and `src/utils.jl` moved/adapted
- **Found during:** Task 2.
- **Issue:** The plan's Task-2 file list named only Project.toml / bayes.jl / ProteinCoLoc.jl /
  ext, but `CoLocResult` was ALSO dispatched on in `src/plot.jl` (three ADVI-result plotters
  using kde/quadgk) and `src/utils.jl` (`generate_txt`). Because the key_links contract requires
  `AdviColocResult` to be defined IN the extension, those symbols become undefined in core, so
  the core module would fail to precompile without Turing.
- **Fix:** Moved the three ADVI-result plotters into `ext/ProteinCoLocTuringExt.jl`; declared
  exported generic-function stubs for them in core; switched `generate_txt` to dispatch on the
  core-defined `AbstractColocResult` supertype (runtime duck-typing on the ADVI result's fields).
- **Files modified:** src/plot.jl, src/utils.jl, ext/ProteinCoLocTuringExt.jl, src/ProteinCoLoc.jl.
- **Commit:** c832119.

### [Design] `bayes_factor` accessor for the ADVI result is a 2-arg delegation
- The ADVI Bayes factor is intrinsically a posterior-vs-prior comparison; a single
  `AdviColocResult` cannot produce it. The interface's 1-arg `bayes_factor(::AdviColocResult)`
  raises an informative error pointing at `bayes_factor(posterior, prior)` /
  `compute_BayesFactor`, and a 2-arg method delegates to `compute_BayesFactor` (as the plan's
  "delegating to compute_BayesFactor" intent requires).

## Authentication Gates
None.

## Known Stubs
- The core generic-function stubs `colocalization`/`compute_BayesFactor`/`plot_posterior`/
  `bayesplot`/`bayes_rangeplot` intentionally have zero methods until Turing is loaded (they are
  the extension interface, D-01). This is by design, not an unfinished stub.

## Self-Check

Created files:
- FOUND: src/results.jl
- FOUND: ext/ProteinCoLocTuringExt.jl
- FOUND: .planning/phases/07-productionization-conditional-on-go/07-00-SUMMARY.md

Commits:
- FOUND: 2b98ee7 (Task 1 — D-02 hierarchy)
- FOUND: c832119 (Task 2 — dependency surgery, v2.0.0)

Verifications:
- Task 2 TOML structure verify: PASS (version 2.0.0; NeuralEstimators/Flux/JLD2/HypothesisTests
  in [deps]; Turing/CUDA in [weakdeps]; ProteinCoLocTuringExt extension; compat pins present).
- Task 3 co-resolution gate: FAILED (GLMakie/Makie conflict — documented; NOT worked around).
- spike/Project.toml + spike/Manifest.toml: provably UNTOUCHED. Root Manifest.toml: unmodified
  (resolve threw before writing).

## Self-Check: PASSED (with Task 3 BLOCKED as designed)
