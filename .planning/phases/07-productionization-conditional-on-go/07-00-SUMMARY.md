---
phase: 07-productionization-conditional-on-go
plan: 00
subsystem: packaging + result-type hierarchy
status: COMPLETE (co-resolution hard gate GREEN)
tags: [breaking-release, package-extension, dependency-surgery, type-hierarchy, co-resolution-gate]
requirements: [PROD-01]
dependency_graph:
  requires: []
  provides:
    - AbstractColocResult hierarchy + shared accessor interface (delta_rho/bayes_factor/is_ood/posterior_draws)
    - AmortizedColocResult shipped subtype
    - ext/ProteinCoLocTuringExt.jl (internal weakdep-gated ADVI reference path)
    - root Project.toml v2.0.0 with hard amortized deps + Turing/CUDA weakdeps
    - GREEN root Pkg.resolve (NeuralEstimators 0.2.1, Flux 0.16.10, Makie 0.24.12, GLMakie 0.13.12)
  affects:
    - all downstream Phase-7 plans (dispatch on AbstractColocResult; unblocked by the green gate)
tech_stack:
  added: [NeuralEstimators (dep, 0.2.1), Flux (dep, 0.16.10), JLD2 (dep), HypothesisTests (dep), CUDA (weakdep), Turing (moved dep -> weakdep), GLMakie (bumped 0.10.5 -> 0.13)]
  patterns: [Julia package extension (weakdep-gated), generic-function stubs + abstract-supertype dispatch, D-02 result-type hierarchy]
key_files:
  created:
    - src/results.jl
    - ext/ProteinCoLocTuringExt.jl
    - .planning/phases/07-productionization-conditional-on-go/07-00-SUMMARY.md
  modified:
    - Project.toml
    - Manifest.toml
    - src/ProteinCoLoc.jl
    - src/plot.jl
    - src/utils.jl
    - src/bayes.jl
    - test/runtests.jl
decisions:
  - D-02 hierarchy: AbstractColocResult supertype + accessor interface + AmortizedColocResult; AdviColocResult internal-only (ext), not exported
  - D-03 surgery: Turing/KDE/QuadGK-requiring ADVI path moved verbatim into ext/ProteinCoLocTuringExt.jl; core precompiles without Turing
  - Co-resolution blocker resolved via user Option 2 - bump GLMakie 0.10.5 -> 0.13 (Makie 0.24.x), GLMakie kept in CORE (not an extension)
  - Gate GREEN - NeuralEstimators pinned 0.2.1 (Finding 1 cleared); Makie 0.24.12 via GLMakie 0.13.12; no src/plot.jl API changes needed at precompile/load
metrics:
  tasks_completed: 3
  tasks_total: 3
  files_created: 3
  files_modified: 7
  completed_date: 2026-07-03
---

# Phase 7 Plan 00: Breaking-Release Foundation + Co-Resolution Gate Summary

**One-liner:** Established the D-02 `AbstractColocResult` type hierarchy, moved the entire
Turing/ADVI reference path into a weakdep package extension, and cleared the Finding-1
co-resolution HARD GATE — bumping GLMakie 0.10.5→0.13 (Makie 0.24.x) resolved the secondary
GLMakie/Makie conflict while keeping NeuralEstimators pinned at 0.2.1, so the root package now
resolves, precompiles, loads, and passes `Pkg.test()` fully green.

## Status: COMPLETE — co-resolution gate GREEN

All three tasks complete. The gate initially failed on a secondary GLMakie/Makie conflict
(escalated as a checkpoint:decision). The user chose **Option 2** — bump GLMakie compat, keep
plotting in the core — and the gate then went green.

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

### Task 3 — Co-resolution HARD GATE — GREEN (commits 44095c8, dbf9d3c)
- Initial resolve FAILED on a secondary GLMakie/Makie conflict (escalated as a
  checkpoint:decision — recorded verbatim in git commit 63420f4).
- **User decision (Option 2):** raised `[compat] GLMakie = "0.10.5" → "0.13"` (Makie 0.24.x),
  keeping GLMakie as a normal CORE dependency (NOT an extension). Regenerated the root
  `Manifest.toml` (fresh full resolve; the stale v1.0.1-era Manifest could not upgrade in place).
- Added the resolve-assertion gate + D-02 unit tests to `test/runtests.jl`; retired the
  now-invalid ADVI-path public-API tests (that path is internal reference, validated later by
  the per-grid D-05 ship-gate under `using Turing`).

## CO-RESOLUTION GATE RESULT: **GREEN**

Resolved versions (from the regenerated `Manifest.toml`, confirmed via `Pkg.dependencies()`):

| Package | Version | Note |
|---------|---------|------|
| NeuralEstimators | **0.2.1** | pinned; NOT downgraded to 0.1.4 — Finding 1 cleared |
| Flux | **0.16.10** | pinned |
| Makie | **0.24.12** | now satisfies NeuralEstimators 0.2.1 |
| GLMakie | **0.13.12** | bumped from 0.10.x; kept as a normal CORE dep |
| Turing / CUDA | (weakdeps) | resolvable weakdeps; extension unchanged |

- `using ProteinCoLoc` **precompiles and loads cleanly** under Makie 0.24 (76 deps precompiled;
  ProteinCoLoc + its extension config build without error).
- **No `src/plot.jl` API changes were needed** at precompile/load level: all GLMakie calls are
  qualified (`GLMakie.*`) and inside function bodies, so the Makie 0.21→0.24 jump does not break
  precompilation. Runtime plotting-attribute correctness under Makie 0.24 is unverified here and
  is deferred with the plotting tests (see Deferred).
- `Pkg.test()` is **fully green**: co-resolution gate 4/4, D-02 hierarchy 10/10, LoadImages
  32/32, patch 3/3, Colocalization 29/29.

## Deviations from Plan

### [Rule 3 - Blocking] Turing-path code in `src/plot.jl` and `src/utils.jl` moved/adapted
- **Found during:** Task 2.
- **Issue:** `CoLocResult` was ALSO dispatched on in `src/plot.jl` (three ADVI-result plotters
  using kde/quadgk) and `src/utils.jl` (`generate_txt`), beyond the plan's named Task-2 files.
  Because the key_links contract requires `AdviColocResult` to live IN the extension, those
  symbols are undefined in core, so the core module would fail to precompile without Turing.
- **Fix:** Moved the three ADVI-result plotters into `ext/ProteinCoLocTuringExt.jl`; declared
  exported generic-function stubs for them in core; switched `generate_txt` to dispatch on the
  core-defined `AbstractColocResult` supertype.
- **Files modified:** src/plot.jl, src/utils.jl, ext/ProteinCoLocTuringExt.jl, src/ProteinCoLoc.jl.
- **Commit:** c832119.

### [Rule 1 - Bug] `cor` regression from removing `using Turing`
- **Found during:** Task 3 verification (`Pkg.test`).
- **Issue:** `colocalization.jl:222` uses a bare `cor` that was reaching scope via
  `using Turing`'s re-export of `Statistics.cor`. Removing `using Turing` broke it
  (`UndefVarError: cor`). It is the ONLY bare Statistics re-export used in core.
- **Fix:** `import Statistics: cor` in `src/ProteinCoLoc.jl`.
- **Commit:** dbf9d3c.

### [Rule 3 - Blocking] Test-harness modernization for `Pkg.test()`
- **Found during:** Task 3 verification.
- **Issue:** The legacy `using .ProteinCoLoc` (leading dot) and package-root-relative fixture
  paths never ran under standard `Pkg.test()`; and `import Pkg` (used by the co-resolution gate)
  needs `Pkg` declared as a test dependency.
- **Fix:** `using ProteinCoLoc`; `cd(dirname(@__DIR__))` to anchor fixture paths; added
  `[extras]`/`[targets]` declaring `Pkg` + `Test`. Also corrected two PRE-EXISTING stale
  assertions to current reality (the `MultiChannelImage` 4-type-parameter signature, and
  `_exclude_zero`'s non-mutating return-value contract) — unrelated to Phase 7 but blocking a
  green suite.
- **Commit:** dbf9d3c.

### [Design] `bayes_factor` accessor for the ADVI result is a 2-arg delegation
- The ADVI Bayes factor is intrinsically a posterior-vs-prior comparison; a single
  `AdviColocResult` cannot produce it. `bayes_factor(::AdviColocResult)` raises an informative
  error pointing at `bayes_factor(posterior, prior)`, and a 2-arg method delegates to
  `compute_BayesFactor`.

## Authentication Gates
None.

## Known Stubs
- The core generic-function stubs `colocalization`/`compute_BayesFactor`/`plot_posterior`/
  `bayesplot`/`bayes_rangeplot` intentionally have zero methods until Turing is loaded (they are
  the extension interface, D-01). This is by design.

## Deferred Items
- **Runtime plotting under Makie 0.24:** the ADVI-result plotters (in the ext) and the core
  `plot`/`local_correlation_plot`/`plot_mask` are precompile/load-clean but their Makie
  0.21→0.24 runtime-attribute correctness is unverified (the ADVI-path plotting tests were
  retired to the D-05 ship-gate; the GLMakie plotting tests are not exercised headless here).
  Logged for a later plan / the ship-gate.

## Self-Check

Created files:
- FOUND: src/results.jl
- FOUND: ext/ProteinCoLocTuringExt.jl
- FOUND: .planning/phases/07-productionization-conditional-on-go/07-00-SUMMARY.md

Commits:
- FOUND: 2b98ee7 (Task 1 — D-02 hierarchy)
- FOUND: c832119 (Task 2 — dependency surgery, v2.0.0)
- FOUND: 63420f4 (Task 3 — gate assertion + initial blocker record)
- FOUND: 44095c8 (Task 3 — gate GREEN, GLMakie bump + regenerated Manifest)
- FOUND: dbf9d3c (Task 3 — green Pkg.test: cor fix + harness modernization)

Verifications:
- Root `Pkg.resolve()`: GREEN. NeuralEstimators 0.2.1, Flux 0.16.10, Makie 0.24.12, GLMakie 0.13.12.
- `using ProteinCoLoc`: loads; AbstractColocResult abstract; AmortizedColocResult <: it;
  colocalization/AdviColocResult NOT exported.
- `Pkg.test()`: fully green (co-resolution 4/4, D-02 10/10, LoadImages 32/32, patch 3/3,
  Colocalization 29/29).
- spike/Project.toml + spike/Manifest.toml: provably UNTOUCHED.

## Self-Check: PASSED
