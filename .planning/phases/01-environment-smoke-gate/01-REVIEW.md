---
phase: 01-environment-smoke-gate
reviewed: 2026-06-26T00:00:00Z
depth: standard
files_reviewed: 3
files_reviewed_list:
  - spike/00_smoke.jl
  - spike/test/runtests.jl
  - spike/Project.toml
findings:
  critical: 0
  warning: 2
  info: 3
  total: 5
status: issues_found
---

# Phase 01: Code Review Report

**Reviewed:** 2026-06-26
**Depth:** standard
**Files Reviewed:** 3
**Status:** issues_found (2 warnings, 3 info — no blockers)

## Summary

Reviewed the Phase 1 environment + smoke gate: `spike/00_smoke.jl` (NPE train + correctness
assert), `spike/test/runtests.jl` (Test harness wrapping the smoke), and `spike/Project.toml`
(spike dependency declarations).

The spike environment is structurally sound. `spike/Manifest.toml` is committed and pins Julia
1.12.6, NeuralEstimators 0.2.1, and Flux 0.16.10. `git diff f581d95 -- Project.toml Manifest.toml`
returns clean — no root-level contamination. The CUDA check in `runtests.jl` is correctly
designed: it queries `Pkg.dependencies()` (which excludes weakdeps) and `Base.loaded_modules`
rather than string-matching `Manifest.toml` text, where CUDA appears only as inert weakdep
metadata. `posteriormean` is a confirmed NeuralEstimators.jl v0.2.1 export (equivalent to
`mean(draws; dims=2)`), verified against `spike/NOTES.md` §2 and `01-RESEARCH.md`.

Two warnings require attention before this gate is used as a reference pattern for later phases.

---

## Warnings

### WR-01: `@assert` inside `include()` converts test failures to unhandled errors

**File:** `spike/test/runtests.jl:46` / `spike/00_smoke.jl:78`

**Issue:** `runtests.jl` line 46 calls `include("00_smoke.jl")` inside the outer `@testset`.
That included file executes `@assert abs(mu_hat - theta_true) < tol ...` at line 78 before
control returns to the testset body. When the assertion fails, it throws `AssertionError` —
not a `Test.Fail` — inside the `@testset` block. Julia's `@testset` catches unhandled exceptions
and records them as `Error During Test`, not `Test Failed`. The two nested testsets (`ENV-02
correctness` and `D-04 CPU-only`) never execute because the exception exits the enclosing body.

Concrete consequence: a broken NPE (API regression, bad training, version mismatch) causes
`runtests.jl` to exit with `Error 1 Total 1` rather than `Fail N Total N`. The CPU-only checks
never run, so a broken build also hides whether the CPU-only constraint held. The NOTES.md §4
already shows this exact failure mode verbatim (`Got exception outside of a @test`).

The `@assert` is the correct mechanism for standalone execution
(`julia --project=spike spike/00_smoke.jl`) where no Test framework is active. The conflict
arises only when `include()`'d into a testset.

**Fix:** Replace the `@assert` in `00_smoke.jl` with a conditional that works in both contexts,
or remove it and rely on the test harness for correctness:

Option A — remove `@assert` from `00_smoke.jl` entirely (simplest; standalone callers see
`println` output rather than a hard abort):
```julia
# Remove line 78 of 00_smoke.jl; println on line 80-81 already reports the result.
```

Option B — guard the assert so it is only active in standalone (non-Test) mode:
```julia
# spike/00_smoke.jl  line 78 replacement
if !isdefined(Main, :Test) || !isdefined(Main.Test, :get_testset)
    @assert abs(mu_hat - theta_true) < tol \
        "smoke FAILED: recovered $mu_hat vs true $theta_true (tol $tol)"
end
```

Option C — in `runtests.jl`, catch `AssertionError` and convert it to a test failure:
```julia
# runtests.jl — replace include() with:
try
    include(joinpath(@__DIR__, "..", "00_smoke.jl"))
catch e
    e isa AssertionError ? (@test false) : rethrow(e)
end
```

---

### WR-02: `Random.seed!(2026)` reproducibility is not thread-count-robust

**File:** `spike/00_smoke.jl:41`

**Issue:** `Random.seed!(2026)` seeds the current task's RNG (Julia 1.7+: task-local RNG). If
`train()` in NeuralEstimators internally spawns tasks via `Threads.@spawn` (e.g., for data
loading or minibatch generation), each spawned task starts with a deterministically derived RNG
state — but the ORDER in which spawned tasks consume that state is scheduling-dependent and
changes with `Threads.nthreads()`. The training result is therefore only byte-reproducible in
the single-threaded case (`--threads=1`, which is the interactive default but not the
`--threads=auto` CI default on multi-core machines).

The 0.3 tolerance absorbs natural variation, so a thread-count change will not flip the gate
from green to red under normal convergence. However, CLAUDE.md mandates "everything reproducible
from `demo.jl` with a fixed seed," which is not strictly satisfied under `--threads=auto`.

**Fix:** Add an explicit assertion or at minimum a comment that documents the single-thread
assumption. The assertion prevents silent reproducibility breakage if Julia's default ever
changes:

```julia
# spike/00_smoke.jl  — after line 41
@assert Threads.nthreads() == 1 """
    smoke seed (Random.seed!(2026)) is only task-local.
    Run with `julia --threads=1 --project=spike spike/00_smoke.jl`
    (or `julia --project=spike ...`) to guarantee reproducibility.
    Detected $(Threads.nthreads()) threads — aborting.
    """
```

If multi-threaded training is desirable for performance, replace `Random.seed!` with
`Random123` (already planned for a later phase) or use `Flux.seed!` / `rng = Xoshiro(2026)` and
thread it explicitly through `prior_sampler` and `gaussian_simulator`.

---

## Info

### IN-01: Correctness tolerance (0.3) is generous but not trivially passing

**File:** `spike/00_smoke.jl:77`

**Issue:** With `m = 200` i.i.d. observations and `sigma = 1.0`, the Bayesian posterior mean for
`theta ~ Normal(0,1)` concentrates sharply: posterior SD ≈ `1/sqrt(m + 1) ≈ 0.07`. A well-
trained NPE recovers `mu_hat ≈ theta_true ± 0.07`, well inside `tol = 0.3`. The question for
adversarial review is whether a pathological but technically-passing result exists.

Analysis: a degenerate estimator that returns the prior mean (0.0) gives
`|0.0 - 0.7| = 0.7 > 0.3` — the gate FAILS. An estimator that outputs random prior samples
(`N(0,1)`, mean ≈ 0 over 1000 draws) also fails. The only way to pass accidentally is to bias
the output toward the data likelihood, which implies the NPE has learned a non-trivial mapping.
The tolerance is not trivially passing.

**Note:** The NOTES.md §4 records the actual green run: `mu_hat ≈ 0.799`, `|delta| = 0.099`.
This is well within 0.3 and consistent with a trained NPE, confirming the gate is discriminating.

---

### IN-02: Implicit `sort()` preprocessing contract is not documented

**File:** `spike/00_smoke.jl:57`, `spike/00_smoke.jl:71`

**Issue:** Both the training simulator (line 57: `sort(randn(m))`) and the test observation
(line 71: `sort(randn(m))`) sort their observation vectors before passing them to the network.
This is internally consistent and deliberately replaces the permutation-invariant DeepSet
(planned for Phase 2) with a simpler `Chain(Dense, Dense, Dense)` that expects a fixed-order
input. The gate is not broken.

However, this creates an undocumented preprocessing contract: any data fed to the trained
estimator in later phases must also be sorted. If Phase 2 passes unsorted patch vectors to the
same summary network architecture, the estimator will silently receive out-of-distribution
inputs. The comment on line 60-61 does not mention the sorted-input requirement.

**Fix:** Add a one-line comment:
```julia
# simulator returns m SORTED replicates; summary network expects sorted input (non-DeepSet).
```
And add the analogous comment next to line 71.

---

### IN-03: `verbose = false` in `train()` suppresses training-divergence diagnostics

**File:** `spike/00_smoke.jl:67`

**Issue:** `verbose = false` is correct for non-cluttered test output. However, if training
diverges (NaN loss, gradient explosion), the failure is invisible until the correctness assertion
fires at the end of the script. `NaN < 0.3` evaluates to `false` in IEEE 754, so a NaN `mu_hat`
still triggers the assertion correctly — there is no false-pass risk here.

The practical concern is debugging: a failed smoke gives no indication of whether the failure
came from diverged training vs. a bad inference call. Running with `verbose = true` fixes this
but requires editing the source.

**Fix (optional):** Document `verbose = false` as deliberate and note the debug form:
```julia
est = train(est, prior_sampler, gaussian_simulator;
            K = 8000, epochs = 60, use_gpu = false,
            verbose = false)   # set to true to see per-epoch loss for debugging
```

---

_Reviewed: 2026-06-26_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
