---
phase: 04-npe-training-advi-benchmark-ablation
plan: 02
subsystem: infra
tags: [julia, turing, advancedvi, dynamicppl, adam, isolated-env, decoupling, read-only-include]

# Dependency graph
requires:
  - phase: 01-environment-smoke-gate
    provides: read-only include() coupling discipline (spike/contract.jl), decoupling baseline f581d95
  - phase: 02-forward-simulator-summary-contract
    provides: build_mci pattern, MultiChannelImage summary contract, ghat(μ)→ρ
provides:
  - "Isolated, pinned spike/baseline/ Julia env (Turing 0.45.0, AdvancedVI 0.6.2, DynamicPPL 0.41.8) — the only safe home for the ADVI ground-truth baseline (D-01)"
  - "spike/baseline/model.jl: read-only lift of the hierarchical colocalization() @model as a top-level coloc_model, plus build_coloc_model(...) constructor"
  - "const DynamicPPL = Turing.DynamicPPL alias pattern that lets src/bayes.jl include cleanly in a minimal env without a dep add or src/ edit"
affects: [04-04-advi-baseline-port, 04-05-benchmark, 05-sbc-bf-ood]

# Tech tracking
tech-stack:
  added: [Turing 0.45.0, AdvancedVI 0.6.2, DataFrames, JLD2, ForwardDiff, StatsBase, Images, ImageFiltering, Random123 (all in the ISOLATED spike/baseline env only)]
  patterns: [isolated-per-baseline Julia env with pinned Manifest, read-only include() one level deeper (../../src), verbatim @model lift, transitive-module aliasing, guarded-optional cross-plan include]

key-files:
  created:
    - spike/baseline/Project.toml
    - spike/baseline/Manifest.toml
    - spike/baseline/model.jl
  modified: []

key-decisions:
  - "Include src/bayes.jl read-only AND lift the nested @model verbatim as top-level coloc_model — the nested @model is unreachable (local to colocalization()'s body), so a verbatim copy is the only way to instantiate it directly (spike 003 MODEL-1 confirmed it captures no enclosing var)"
  - "Alias const DynamicPPL = Turing.DynamicPPL instead of adding DynamicPPL as a direct dep — bayes.jl references the bare module name; aliasing keeps the Task-1 dep set intact and touches neither src/ nor the env"
  - "smoke_mci_pair synthetic fixture as the self-contained construction smoke — the physics simulator (spike/simulator/forward.jl) cannot load here (its ImageTransformations/CoordinateTransformations/Interpolations deps are deliberately excluded, D-11 minimalism); re-simulation via 04-01 resimulate.jl is a guarded-optional hook"

patterns-established:
  - "Isolated baseline env (spike/baseline/) with its own pinned Project/Manifest — mirrors the spike/ isolation but carries the Turing stack the spike env must never see (Phase-1 co-resolve landmine)"
  - "Read-only include(joinpath(@__DIR__, '..', '..', 'src', ...)) coupling from one directory deeper than contract.jl"
  - "Guarded cross-plan include (isfile + try/catch) so a parallel-wave dependency (04-01 resimulate.jl) lights up when merged but never blocks load"

requirements-completed: [NPE-02]

# Metrics
duration: ~13min
completed: 2026-07-01
---

# Phase 4 Plan 02: Isolated ADVI Baseline Env + Read-only @model Lift Summary

**Pinned, fully isolated spike/baseline/ env (Turing 0.45.0 / AdvancedVI 0.6.2 / DynamicPPL 0.41.8) with a read-only lift of the hierarchical colocalization() @model that instantiates a DynamicPPL.Model from a MultiChannelImage pair — de-risking the ADVI port (04-04) with Turing provably absent from the spike env and src/ untouched.**

## Performance

- **Duration:** ~13 min
- **Started:** 2026-07-01T18:02Z (approx, first baseline dir creation)
- **Completed:** 2026-07-01T18:10:30Z
- **Tasks:** 2
- **Files modified:** 3 created

## Accomplishments
- Resolved and pinned a fresh isolated `spike/baseline/` env with all 9 required direct deps; KernelAbstractions deliberately kept out of `[deps]` (D-11).
- Confirmed the Phase-1 landmine is respected: `Turing` is absent from the spike env (`julia --project=spike` assertion green) and `src/` stays byte-clean.
- Lifted the hierarchical `@model` read-only into `spike/baseline/model.jl` and exposed `build_coloc_model(...)`, which constructs a `DynamicPPL.Model` from a `MultiChannelImage` sample/control pair via the UNCHANGED `_prepare_data`.
- Recorded the resolved baseline versions for 04-04: **Turing 0.45.0, AdvancedVI 0.6.2, DynamicPPL 0.41.8** (RESEARCH warned the `vi` return type differs across Turing 0.43→0.45 and expected AdvancedVI ~0.7 — the resolved 0.6.2 is the exact version the port must target).

## Task Commits

Each task was committed atomically:

1. **Task 1: Resolve and pin the isolated baseline environment** - `ede3f69` (chore)
2. **Task 2: Lift the colocalization @model read-only into the baseline env** - `216e9ad` (feat)

## Files Created/Modified
- `spike/baseline/Project.toml` - isolated baseline `[deps]`: Turing, AdvancedVI, DataFrames, JLD2, ForwardDiff, StatsBase, Images, ImageFiltering, Random123
- `spike/baseline/Manifest.toml` - pinned reproducibility artifact (Turing 0.45.0 stack)
- `spike/baseline/model.jl` - read-only include() of src/LoadImages.jl + src/colocalization.jl + src/bayes.jl; verbatim top-level `coloc_model`; `build_coloc_model(...)`; `smoke_mci_pair` fixture; guarded `spike/npe/resimulate.jl` hook

## Decisions Made
- **Include bayes.jl AND lift the @model verbatim:** the `@model` is nested inside `colocalization()` and unreachable as a top-level name, so including bayes.jl (for `_prepare_data`/`CoLocResult`/`compute_BayesFactor` and to satisfy the key_links `include.*src.*bayes` coupling) is combined with a verbatim top-level copy of the model block for direct instantiation.
- **`const DynamicPPL = Turing.DynamicPPL` alias:** bayes.jl uses the bare `DynamicPPL` name (its `convert_posterior_samples` signature). `using Turing` loads but does not bind that name, and `using DynamicPPL` fails (transitive-only). Aliasing from Turing is the minimal, decoupling-faithful fix — no dep added, no src/ edit.
- **Synthetic `smoke_mci_pair` for the construction smoke:** the calibrated physics simulator can't load in this env by design (its transform deps are excluded), so a dependency-light correlated fixture proves `build_coloc_model` end-to-end; faithful re-simulation is deferred to 04-01's `resimulate_holdout` via a guarded hook.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] `DynamicPPL` name unbound when including src/bayes.jl**
- **Found during:** Task 2 (lifting the @model)
- **Issue:** `src/bayes.jl:24` references the bare module `DynamicPPL` in a method signature (evaluated at include time). `using Turing` loads DynamicPPL but does not bind the name in `Main`; `using DynamicPPL` errors because it is a transitive (not direct) dep of the baseline env.
- **Fix:** Added `const DynamicPPL = Turing.DynamicPPL` before the `src/` includes.
- **Files modified:** spike/baseline/model.jl
- **Verification:** `include("spike/baseline/model.jl")` loads; `build_coloc_model` returns a `DynamicPPL.Model`.
- **Committed in:** 216e9ad (Task 2 commit)

**2. [Rule 3 - Blocking] Physics simulator cannot load in the minimal baseline env**
- **Found during:** Task 2 (planned "include the simulator + resimulate.jl for a re-simulated smoke pair")
- **Issue:** `spike/simulator/forward.jl` `using`s ImageTransformations / CoordinateTransformations / Interpolations, which are intentionally NOT in the Task-1 dep set (D-11 minimalism); a direct include breaks model.jl load. 04-01's `spike/npe/resimulate.jl` is also absent (same parallel wave, `depends_on: []`).
- **Fix:** Guarded the `resimulate.jl` include (`isfile` + `try/catch`, exposed via `RESIMULATE_AVAILABLE`) and added a self-contained `smoke_mci_pair` fixture as the construction smoke input. No new deps added (kept the pinned env minimal).
- **Files modified:** spike/baseline/model.jl
- **Verification:** model.jl loads with `RESIMULATE_AVAILABLE[] == false`; `build_coloc_model(smoke_mci_pair()..., [1,2], 8)` constructs a model.
- **Committed in:** 216e9ad (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (both Rule 3 - blocking).
**Impact on plan:** Both fixes were necessary to complete Task 2 in a genuinely isolated, minimal env without editing `src/` or expanding the pinned dep set. No scope creep — the acceptance criteria (model loads, constructs from an MCI pair, no inference call, src/ clean) are all met. The re-simulation smoke path is preserved as a guarded hook that will activate when 04-01 merges.

## Issues Encountered
- **AdvancedVI resolved to 0.6.2, not the RESEARCH-anticipated ~0.7.** Not a blocker for this plan (no inference run here) but flagged for 04-04: verify the `vi(model, family, max_iter; adtype=...)` signature and `VIResult`/`.q` return shape against AdvancedVI 0.6.2 specifically, not 0.7 docs.

## User Setup Required
None - no external service configuration required (offline, CPU-only research spike).

## Next Phase Readiness
- `spike/baseline/` env is pinned and reproducible; `build_coloc_model(...)` is ready for 04-04 to drive the modern `vi()` ADVI port (`vi(m, q_meanfield_gaussian, ITER; adtype=AutoForwardDiff())` → `VIResult.q`).
- **04-04 handoff notes:** (1) target AdvancedVI **0.6.2** API; (2) `convert_posterior_samples`/`colocalization` in the included bayes.jl still use the REMOVED `DynamicPPL.syms` / `vi(m, ADVI(n,iter))` API — 04-04 must port the inference call in the baseline (not in src/), per D-01; (3) re-simulation of raw holdout images depends on 04-01's `resimulate_holdout`, wired here as a guarded hook.

## Threat Surface
No new threats. Both threat-model entries satisfied: T-04-SC2 (all deps flagship JuliaStats/TuringLang, pinned Manifest committed) and T-04-DECOUPLE (no src/ writes — `git status src/` clean; coupling via read-only include() only).

## Self-Check: PASSED

- Files verified present: spike/baseline/Project.toml, spike/baseline/Manifest.toml, spike/baseline/model.jl, 04-02-SUMMARY.md
- Commits verified: ede3f69 (Task 1), 216e9ad (Task 2), 3d9b883 (SUMMARY)
- src/ and spike/ env (Project.toml/Manifest.toml) byte-clean — decoupling and isolation intact

---
*Phase: 04-npe-training-advi-benchmark-ablation*
*Completed: 2026-07-01*
