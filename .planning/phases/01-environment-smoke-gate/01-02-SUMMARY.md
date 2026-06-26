---
phase: 01-environment-smoke-gate
plan: 02
subsystem: infra
tags: [julia, neuralestimators, flux, sbi, npe, normalising-flow, cpu-only, testing, spike]

# Dependency graph
requires:
  - phase: 01-01
    provides: frozen root baseline (commit f581d95) + recorded baseline ref for the untouched-root proof
provides:
  - isolated spike/ environment (own Project.toml + resolved Manifest.toml) with exactly NeuralEstimators 0.2.1, Flux 0.16.10, Distributions 0.25.128 (no installed CUDA)
  - spike/00_smoke.jl — CPU-only NPE smoke that trains a PosteriorEstimator(NormalisingFlow) on a 1-param Gaussian and asserts recovered posterior mean within tol of a known theta
  - spike/test/runtests.jl — re-runnable stdlib Test gate (correctness + CUDA-absence) usable as the hard gate before every later phase
  - verified-against-installed v0.2.1 API call sequence (resolves the three pre-1.0 drift hotspots)
affects: [01-03, 01-04, 02-simulator, 04-advi-baseline, 05-sbc-bf-ood]

# Tech tracking
tech-stack:
  added: [NeuralEstimators v0.2.1, Flux v0.16.10, Distributions v0.25.128]
  patterns: [isolated Pkg.activate("spike") env, verify-against-installed-docstrings for pre-1.0 deps, green+correctness smoke gate, CPU-only via use_gpu=false + installed-set CUDA guard]

key-files:
  created: [spike/Project.toml, spike/Manifest.toml, spike/00_smoke.jl, spike/test/runtests.jl]
  modified: []

key-decisions:
  - "q is passed POSITIONALLY as a constructed NormalisingFlow instance: PosteriorEstimator(network, q). The convenience constructor's q= keyword expects a TYPE, not an instance — the planned q=<instance> keyword form was wrong."
  - "CUDA-absence is asserted via Pkg.dependencies()/project deps/loaded modules, NOT a Manifest text regex — Flux/NNlib/Zygote/NeuralEstimators declare CUDA as inert weakdep extensions that false-positive a \\bCUDA\\b regex."
  - "m=200 replicates per data set (informative posterior) + K=8000 / epochs=60 gives a 3x margin (|delta|=0.099 < tol=0.3); NPE mu_hat matches the analytic Bayesian posterior mean to within 0.004."

patterns-established:
  - "Verify-against-installed-docstrings: read ?PosteriorEstimator/?NormalisingFlow/?sampleposterior/?train on the resolved env before writing code against a pre-1.0 dependency."
  - "Green+correctness smoke gate: assert a recovered quantity within tolerance, not merely runs-without-error."
  - "CPU-only contract enforced on two fronts: no installed CUDA package + train(...; use_gpu=false), guarded by a re-runnable Test assertion."

requirements-completed: [ENV-02]

# Metrics
duration: 22min
completed: 2026-06-26
---

# Phase 1 Plan 02: Environment + Smoke Gate Summary

**Isolated CPU-only spike env (NeuralEstimators 0.2.1 + Flux 0.16.10) whose green stdlib-Test gate trains a NormalisingFlow PosteriorEstimator on a 1-param Gaussian and recovers theta within a 3x-margin tolerance, with CUDA proven absent.**

## Performance

- **Duration:** ~22 min
- **Completed:** 2026-06-26
- **Tasks:** 2
- **Files modified:** 4 created

## Accomplishments
- Stood up the minimal, provably isolated `spike/` environment (own `Project.toml` + resolved `Manifest.toml`) with exactly NeuralEstimators 0.2.1, Flux 0.16.10, Distributions 0.25.128 — and confirmed no installed CUDA package node (only inert weakdep extension references remain).
- Wrote `spike/00_smoke.jl` against the INSTALLED v0.2.1 docstrings, resolving all three pre-1.0 drift hotspots; it trains a `PosteriorEstimator(NormalisingFlow)` CPU-only and asserts the recovered posterior mean (mu_hat=0.79911894) is within tol=0.3 of theta_true=0.7 (|delta|=0.099 — a comfortable 3x margin; mu_hat matches the analytic posterior mean to 0.004).
- Wrapped the smoke in `spike/test/runtests.jl`, a re-runnable stdlib Test gate that asserts both correctness and CUDA-absence; `julia --project=spike spike/test/runtests.jl` exits 0 deterministically across consecutive runs.
- Retired the phase's two top risks: P4 (NeuralEstimators/Flux pre-1.0 API churn) and P5 (Windows/GPU breakage), with the root baseline (f581d95) provably untouched.

## Task Commits

Each task was committed atomically:

1. **Task 1: Create + resolve minimal isolated spike env** - `56f85c3` (chore)
2. **Task 2: CPU-only NPE smoke + re-runnable Test gate** - `6124c62` (feat)

**Plan metadata:** (this commit) (docs: complete plan)

## Files Created/Modified
- `spike/Project.toml` - Minimal `[deps]`: NeuralEstimators, Flux, Distributions (no heavy root deps, no CUDA).
- `spike/Manifest.toml` - Resolved reproducibility artifact (NeuralEstimators 0.2.1, Flux 0.16.10, Distributions 0.25.128 + transitive set); no installed CUDA node.
- `spike/00_smoke.jl` - Readable CPU-only NPE smoke with AGPL header, `Random.seed!(2026)`, and the posterior-mean correctness assertion (D-03/D-04).
- `spike/test/runtests.jl` - Stdlib Test hard gate (D-06): includes the smoke, asserts the tolerance, and asserts CUDA is neither a dependency nor loaded.

## Decisions Made
- **q passed positionally as an instance.** Installed v0.2.1 exposes `PosteriorEstimator(summary_network, q::ApproximateDistribution)` (instance, positional) and a separate convenience form `PosteriorEstimator(summary_network, d; num_summaries, q=<Type>)` where `q` is a TYPE. Used the positional-instance form with `q = NormalisingFlow(d; num_summaries=...)`.
- **m=200 replicates per data set.** A single observation (the plan skeleton's m=1) makes the posterior mean track one noisy draw, not theta. Using 200 sorted replicates makes the data informative so both the analytic and NPE posterior means sit near theta_true with a 3x tolerance margin.
- **Tolerance tol=0.3** (Claude's discretion per D-03) — generous enough to separate a correct backend from a broken one, with the actual run at |delta|=0.099.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Plan's CUDA-absence check (`\bCUDA\b` Manifest regex) is a false-positive**
- **Found during:** Task 1 (env verify) — the plan's verify command `@assert !occursin(r"\bCUDA\b", manifest)` failed on a correct, CUDA-free environment.
- **Issue:** Flux/NNlib/Zygote/NeuralEstimators/MLDataDevices/Atomix declare CUDA as a *weak* dependency (optional package extension — the CUDA-as-extension design since Flux v0.14). Those `[deps.X.weakdeps]`/`[deps.X.extensions]` entries are inert text in `Manifest.toml` and match the regex even though CUDA is never installed or resolved.
- **Fix:** Verified the env-level CPU-only guarantee with the correct check: no installed package node named CUDA (`Pkg.dependencies()`), no `[[deps.CUDA]]` stanza, CUDA not in the project's direct deps. The Test gate's D-04 assertions use `Pkg.dependencies()` (installed set), `Pkg.project().dependencies` (direct deps), and `Base.loaded_modules` (runtime) instead of a Manifest text regex.
- **Files modified:** spike/test/runtests.jl (CUDA-absence asserted on the installed/loaded set, not Manifest text).
- **Verification:** Confirmed `NONE INSTALLED`, no `[[deps.CUDA]]` stanza, `haskey(...,"CUDA")==false`; gate's 3 CUDA assertions pass.
- **Committed in:** 6124c62 (Task 2 commit).

**2. [Rule 1 - Bug] Plan interface note `PosteriorEstimator(network, d; ..., q=<instance>)` mis-specified**
- **Found during:** Task 2 (docstring verification — exactly the drift this step exists to catch, P4).
- **Issue:** In installed v0.2.1, the `q=` keyword of the convenience constructor expects a TYPE (default NormalisingFlow); passing a constructed instance there is the wrong overload.
- **Fix:** Used the 2-arg positional form `PosteriorEstimator(summary_network, q)` with `q` a constructed `NormalisingFlow` instance, per the canonical installed docstring example.
- **Files modified:** spike/00_smoke.jl.
- **Verification:** Smoke trains and samples successfully; recovered mean within tolerance.
- **Committed in:** 6124c62 (Task 2 commit).

---

**Total deviations:** 2 auto-fixed (both Rule 1 — correctness of plan-supplied checks/signatures).
**Impact on plan:** Both were the pre-1.0 API-drift and weakdep-vs-installed nuances this smoke gate exists to surface. No scope creep; the deliverables match the plan's `files_modified` exactly.

## Issues Encountered
- First-call precompilation of NeuralEstimators/Flux/Zygote took ~60s and the full first gate run ~3m02s (training included); subsequent runs are faster. No failures — patience per the plan note paid off.

## Known Stubs
None — the smoke is a complete, correctness-asserting gate; nothing is stubbed.

## User Setup Required
None - no external service configuration required. (Optional, per D-05/Plan 01-03: a directory-scoped Julia pin via `juliaup override set 1.12.6` and a recorded `.julia-version`/NOTES.md — that pinning + stack-decision note is Plan 01-03 scope.)

## Next Phase Readiness
- The hard gate is established: `julia --project=spike spike/test/runtests.jl` is green and re-runnable, the durable precondition for all later phases.
- `spike/Manifest.toml` exists and is committed; Plan 01-03 will freeze/pin it further and add the Julia-version record + NOTES.md (ENV-03/ENV-04).
- Plan 01-04 will attempt `Pkg.develop(path="..")` coupling (with the `include()` fallback) — independent of this gate.
- Root `Project.toml`/`Manifest.toml`/`src/` remain provably untouched vs baseline f581d95.

## Self-Check: PASSED
- Files: spike/Project.toml, spike/Manifest.toml, spike/00_smoke.jl, spike/test/runtests.jl — all FOUND.
- Commits: 56f85c3, 6124c62 — both FOUND.
- Root untouched vs f581d95: CLEAN.

---
*Phase: 01-environment-smoke-gate*
*Completed: 2026-06-26*
