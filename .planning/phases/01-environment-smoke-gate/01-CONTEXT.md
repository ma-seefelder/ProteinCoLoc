# Phase 1: Environment + Smoke Gate - Context

**Gathered:** 2026-06-26
**Status:** Ready for planning

<domain>
## Phase Boundary

Stand up a reproducible, isolated `spike/` environment with its own `Project.toml`/`Manifest.toml`, make the main package reachable read-only, and prove the pre-1.0 NeuralEstimators + Flux stack works on this machine via a green smoke test — retiring the highest-risk unknown before any further investment. The main `Project.toml`/`Manifest.toml` and `src/` stay provably untouched. Simulator, training data, and inference belong to later phases.

Requirements: ENV-01, ENV-02, ENV-03, ENV-04.
</domain>

<decisions>
## Implementation Decisions

### Parent-package coupling (ENV-01)
- **D-01:** The spike reaches the existing `src/` functions via `Pkg.develop(path="..")` **first** (clean package boundary). If co-resolving the parent dep tree (GLMakie 0.10.5, Turing, PackageCompiler) with Flux 0.16 conflicts or is unreasonably heavy, **document the conflict and fall back to `include()`** of the specific source files for the spike. Either path must leave the main `Project.toml`/`Manifest.toml` and `src/` untouched.
- **D-02:** Rationale for preferring `Pkg.develop`: Phase 4's ADVI baseline needs the *real* Turing `colocalization()` model through the package, so keeping the package boundary now avoids rework later; the `include()` fallback exists only to keep Phase 1 unblocked if the resolve fails.

### Smoke test (ENV-02)
- **D-03:** The smoke gate is **green + correctness**: it must `train` a `PosteriorEstimator` (NormalisingFlow) on a 1-parameter Gaussian, run `sampleposterior`, **and assert the recovered posterior mean is within a tolerance of the known value** — not merely "runs without error." Catches a silently-wrong API/backend.
- **D-04:** CPU-only, with **no forced CUDA import** (ROADMAP SC2). The core stays close to the ~30-line target; the correctness assertion adds a few lines.

### Reproducibility artifact (ENV-03)
- **D-05:** Pin **both** the package versions (`spike/Manifest.toml`, exact NeuralEstimators/Flux + transitive versions captured after the first `Pkg.resolve()`) **and the Julia version** (record exact `julia --version`; add a `.julia-version` / juliaup-channel note). Commit as the reproducibility artifact. The main repo declares no Julia compat bound and the SBI stack is version-sensitive, so the Julia version is part of "reproducible from a fixed seed."

### Gate durability (ENV-03 / cross-phase)
- **D-06:** Wrap the smoke in a **re-runnable harness** `spike/test/runtests.jl` (stdlib `Test`), so "green smoke is the hard gate for all later phases" can be re-checked on demand before each subsequent phase. `spike/00_smoke.jl` remains the readable entry script. **No external CI** in scope.

### Stack-decision note (ENV-04, locked to default)
- **D-07:** Record NeuralEstimators.jl as default and BayesFlow (PythonCall) as fallback-only in a spike note. Default location `spike/NOTES.md` (also seeds the place where Phase 2 will copy the Turing prior ranges); a dedicated `spike/STACK-DECISION.md` is acceptable at planner discretion.

### Claude's Discretion
- Exact tolerance value for the correctness assertion.
- File/dir layout within `spike/` (beyond the named `00_smoke.jl`, `test/runtests.jl`, `NOTES.md`, `Project.toml`/`Manifest.toml`).
- Exact juliaup/`.julia-version` mechanism for recording the Julia version.
- `NOTES.md` vs a dedicated stack-decision file for D-07.
</decisions>

<specifics>
## Specific Ideas

- The user is the package author and weights the **decoupling / publication-integrity guarantee** highly — the publication state is frozen at tags `v1.0.2-scirep` (65aad18) and `v2.0-baseline` (HEAD, 0d9b87c), both pushed to origin. Phase 1 must visibly preserve that: nothing under `src/` or the root manifests changes.
- Reproducibility is treated as a first-class Definition-of-Done concern, hence pinning the Julia version, not just packages.
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Stack & API
- `.planning/research/STACK.md` — NeuralEstimators v0.2.1 + Flux v0.16 versions, the exact `Pkg.add` install list, AP1 smoke API notes (`PosteriorEstimator`/`NormalisingFlow`/`train`/`sampleposterior`), CUDA-as-extension behavior, and the "What NOT to Use" list (NormalizingFlows.jl/InvertibleNetworks.jl excluded). §"Key API Notes (AP1 smoke target)" and §"Version Compatibility" are directly relevant.
- `.planning/research/ARCHITECTURE.md` — L0 smoke layer, the 7-layer horizontal build order, and the decoupling-boundary recommendation (`Pkg.develop` parent into spike env vs raw `include`).
- `.planning/research/PITFALLS.md` — P4 (NeuralEstimators/Flux pre-1.0 churn → pinned-Manifest + smoke-gate mitigation) and P5 (Windows/GPU breakage → CPU-only).

### Requirements & scope
- `.planning/REQUIREMENTS.md` — ENV-01..04 (the four Phase-1 requirements).
- `.planning/ROADMAP.md` §"Phase 1: Environment + Smoke Gate" — the four success criteria.
- `.planning/PROJECT.md` §Constraints — decoupling (hard), CPU-first / GPU-optional, reproducibility-from-fixed-seed.
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets (exposed via the D-01 coupling, used in Phase 2+)
- `src/colocalization.jl` — `correlation(x,y;method)`, `patch()`, `_prepare_data()` (the summary-statistic contract).
- `src/LoadImages.jl` — `MultiChannelImage`, `AbstractMultiChannelImage`, `_apply_mask!`, `otsu_thresholds`.
- `src/bayes.jl` — `compute_BayesFactor()` (KDE baseline, Phase 6).
- These are NOT used by the Phase-1 smoke itself (a synthetic 1-param Gaussian), but the coupling mechanism chosen here must expose them for later phases.

### Established Patterns / constraints on the resolve
- Main `Project.toml` deps that bear on D-01's resolve weight: `Turing`, `GLMakie` (pinned `0.10.5`), `PackageCompiler`, `Distributions`, `Images`, `KernelDensity`, `QuadGK`, `Random123`, `StatsBase`. GLMakie + Turing are the heavy items the spike does not need for the smoke.
- Sibling project `BayesInteractomics` provides JLD2-cache and `_bin_calibration`/`CalibrationResult` patterns for later phases (not Phase 1).

### Integration Points
- `spike/` is a sibling directory at the repo root with its own activated environment; `Pkg.develop(path="..")` targets the repo root (the ProteinCoLoc package).
</code_context>

<deferred>
## Deferred Ideas

- **CUDA / GPU acceleration** — only revisited if Phase 4 training is slow on CPU (PROJECT constraint); kept out of the Phase-1 env.
- **External CI** for the smoke gate — out of scope; the re-runnable local harness (D-06) is the chosen guard.
- **RxInfer baseline** — v2/deferred (BACK-01), not a spike dependency.

### Reviewed Todos (not folded)
None — no pending todos matched this phase.
</deferred>

---
*Phase: 01-environment-smoke-gate*
*Context gathered: 2026-06-26*
