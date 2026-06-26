# Phase 1: Environment + Smoke Gate - Pattern Map

**Mapped:** 2026-06-26
**Files analyzed:** 6 new files (all under `spike/`)
**Analogs found:** 2 strong in-repo / 6 total (the SBI smoke + the stack docs have NO in-repo analog — external API only)

> **Context for the planner:** Phase 1 is greenfield inside a NEW isolated `spike/` directory. Only TWO of the six deliverables have a genuine in-repo analog (`spike/Project.toml` ← root `Project.toml`; `spike/test/runtests.jl` ← root `test/runtests.jl`). The NeuralEstimators NPE smoke, the pinned Manifest, the stack-decision note, and the Julia-version record have **no local precedent** — for those, mirror the canonical external API sequence captured verbatim in `01-RESEARCH.md` §"Code Examples" rather than any invented local analog. This is called out honestly per file below.

---

## File Classification

| New File | Role | Data Flow | Closest Analog | Match Quality |
|----------|------|-----------|----------------|---------------|
| `spike/Project.toml` | config (Pkg manifest) | declarative | `Project.toml` (root) | role-match (must be MINIMAL, not a copy) |
| `spike/Manifest.toml` | config (generated artifact) | batch / generated | `Manifest.toml` (root) | structural only (DO NOT hand-write — `Pkg.resolve()` emits it) |
| `spike/00_smoke.jl` | script (smoke harness) | transform / request-response | `src/*.jl` (license header only) | **NO behavioral analog** — external NeuralEstimators API |
| `spike/test/runtests.jl` | test | request-response | `test/runtests.jl` (root) | role-match (structure only; subject is synthetic, not images) |
| `spike/NOTES.md` | doc (stack-decision + Julia-version record) | declarative | — | **NO analog** — content from CLAUDE.md / STACK.md |
| `spike/.julia-version` | config (version record) | declarative | — | **NO analog** — single literal line `1.12.6` |

---

## Pattern Assignments

### `spike/Project.toml` (config, declarative)

**Analog:** `C:\Users\Manuel\Documents\GitHub\ProteinCoLoc\Project.toml` (lines 1-22)

**Structure to mirror** (the `name`/`uuid`/`version` + `[deps]` + `[compat]` skeleton):
```toml
name = "ProteinCoLoc"
uuid = "12345678-1234-1234-1234-123456789abc"
version = "1.0.1"

[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
...
[compat]
GLMakie = "0.10.5"
```

**What to COPY:** the `[deps]` UUID-per-line table format and the trailing `[compat]` block convention. The exact UUIDs for `Distributions` (`31c24e10-a181-5473-b8eb-7969acd0382f`), `Random` (`9a3f8284-a2c9-5f02-9a11-845980a1fd5c`), and `Statistics` (`10745b16-79ce-11e8-11f9-7d13ad32a3b2`) are already present in the root and can be reused verbatim.

**What to CHANGE / NOT copy (critical per RESEARCH §Standard Stack, line 74):**
- Phase-1 `spike/Project.toml` must be **MINIMAL**: only `NeuralEstimators`, `Flux`, `Distributions` (+ `Test` is stdlib, may be omitted from `[deps]`). Do NOT copy the root's heavy deps (`Turing`, `GLMakie`, `PackageCompiler`, `Images`, `CSV`, `DataFrames`, `KernelDensity`, `QuadGK`) — those widen the resolve and the Windows-fragility surface the smoke exists to isolate.
- A `spike/` env normally does NOT need its own `name`/`uuid`/`version` header (it is an environment, not a published package) — a bare `[deps]`/`[compat]` file is sufficient. If `Pkg.test()` package-context is later desired (RESEARCH line 336), a header can be added at planner discretion.
- Use `Pkg.add([...])` to populate this (RESEARCH lines 90-97), then let `Pkg.resolve()` write the `[compat]` bounds — do not hand-pin Flux/NeuralEstimators versions in `Project.toml`; pinning lives in the Manifest (ENV-03).

---

### `spike/Manifest.toml` (config, GENERATED artifact)

**Analog:** root `Manifest.toml` (structural reference only — 114 KB, machine-generated).

**Pattern:** This file is **NEVER hand-written**. It is emitted by `Pkg.resolve()` / `Pkg.add()` inside the activated `spike/` env (RESEARCH Pattern 3, lines 191-201). The executor's job is to (a) run the resolve, (b) confirm the smoke is green against it, then (c) **commit it verbatim** as the reproducibility artifact (ENV-03 / D-05).

**Verification assertion the planner must schedule** (Pitfall 2, RESEARCH lines 229-233): after resolve, confirm the Manifest contains **no `CUDA` / `CUDA_jll` entry**. This is the env-level half of the CPU-only guarantee (the `use_gpu=false` flag is the other half).

**No code excerpt applies** — there is nothing to copy; the analog is purely "a resolved Manifest looks like the root one and is committed."

---

### `spike/00_smoke.jl` (script, transform) — NO IN-REPO BEHAVIORAL ANALOG

**Analog (behavior):** NONE. There is no existing NeuralEstimators / Flux / normalising-flow code anywhere in this repo (verified: `src/` is image-loading + Turing ADVI + Makie plotting; no NN backend). Inventing a local analog would be dishonest.

**Authoritative source instead:** `01-RESEARCH.md` §"Code Examples" lines 273-304 (the verified-against-v0.2.1 smoke skeleton) and lines 252-271 (the verbatim API signatures). **Mirror that sequence**, not memory or stale STACK.md notes.

**The one in-repo convention that DOES apply — the license header** (from `src/colocalization.jl` lines 1-19, identical block atop every `src/*.jl`):
```julia
#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
...
You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=#
```
Open the smoke script with this same AGPL `#= ... =#` header for repo consistency (note: `test/runtests.jl` uses a shorter MIT-style 5-line `#####` banner — the executor should pick ONE; the `#=` AGPL block is the dominant `src/` convention).

**Core API pattern to mirror** (from RESEARCH lines 274-304 — verify each call against installed `?` docstrings per Pattern 1, lines 169-172):
```julia
using NeuralEstimators, Flux, Distributions, Random
Random.seed!(2026)                       # determinism — Pitfall 5

d = 1; σ = 1.0f0; m = 1
sample(K)  = Float32.(randn(1, K))                                   # prior θ~N(0,1), d×K
simulate(θ) = [Float32.(θ[1,k] .+ σ.*randn(1,m)) for k in 1:size(θ,2)]

num_summaries = 16
network = Chain(Dense(1 => 32, relu), Dense(32 => num_summaries))
q   = NormalisingFlow(d, num_summaries)
est = PosteriorEstimator(network, d; num_summaries = num_summaries, q = q)
est = train(est, sample, simulate; epochs = 200, use_gpu = false)    # CPU-only, D-04

θ_true = 0.7f0
Z_obs  = Float32.(θ_true .+ σ .* randn(1, m))
θ̂      = sampleposterior(est, Z_obs; N = 1000)                        # d×N matrix
μ̂      = posteriormean(θ̂)[1]
tol = 0.3
@assert abs(μ̂ - θ_true) < tol  "smoke FAILED: recovered $μ̂ vs true $θ_true"
```

**Three drift hotspots the executor MUST re-verify against the installed v0.2.1 docstrings** (RESEARCH Pitfall 1 + Assumptions A3, lines 223-227 / 390):
1. `PosteriorEstimator(network, d; num_summaries=d, q=...)` — NOT the stale `PosteriorEstimator(q, network)`.
2. `sampleposterior(est, Z; N=1000)` — `N` is now a **keyword**, not positional.
3. Whether `q` accepts a constructed `NormalisingFlow` instance vs a type (A3, the single most likely construction snag).

---

### `spike/test/runtests.jl` (test, request-response)

**Analog:** `C:\Users\Manuel\Documents\GitHub\ProteinCoLoc\test\runtests.jl` (lines 1-17 — the harness preamble + first `@testset`).

**Harness preamble pattern to mirror** (root lines 6-17):
```julia
import Images
import Statistics: cor
using .ProteinCoLoc
using Random123
using Test
Random123.seed!(1234)

@testset "LoadImages" verbose = true begin
    ...
end
```

**What to COPY (structure):**
- The `using Test` + top-level `@testset "name" verbose = true begin ... end` wrapper (root line 17).
- The deterministic seed-at-top discipline (root line 12). **BUT:** Phase 1 uses stdlib `Random.seed!`, NOT `Random123` — RESEARCH Pitfall 5 (line 248) explicitly defers Random123 to later phases. Replace `using Random123 / Random123.seed!(1234)` with `using Random / Random.seed!(2026)`.
- `@test <condition>` assertion style (root lines 22, 27, 30, etc.).

**What to CHANGE / NOT copy:**
- The root test loads real TIFF images and exercises `MultiChannelImage` (lines 17-62). The spike testset's subject is the **synthetic 1-param Gaussian smoke** — none of the image machinery applies. Do NOT pull in `Images` or `using .ProteinCoLoc` (the smoke needs nothing from the parent package; coupling is verified separately per Pattern 2).
- Wrap the smoke via `include(joinpath(@__DIR__, "..", "00_smoke.jl"))` then assert (RESEARCH lines 305-313), OR refactor the smoke body into a callable function — planner's choice.

**Two assertions this harness must contain** (RESEARCH Test Map, lines 338-344):
```julia
@testset "NeuralEstimators CPU smoke" begin
    include(joinpath(@__DIR__, "..", "00_smoke.jl"))
    @test abs(μ̂ - θ_true) < tol                       # ENV-02 correctness (D-03)
    @test !haskey(Pkg.project().dependencies, "CUDA")  # D-04 CPU-only guard (adjust check)
end
```
**Run command** (RESEARCH line 335): `julia --project=spike spike/test/runtests.jl`.

---

### `spike/NOTES.md` (doc) — NO IN-REPO ANALOG

**Analog:** NONE — there is no existing `NOTES.md` / decision-doc pattern in the repo.

**Content source:** CLAUDE.md §"What NOT to Use" (NeuralEstimators default; BayesFlow/PythonCall fallback-only; NormalizingFlows.jl & InvertibleNetworks.jl excluded) + `01-RESEARCH.md` D-07 (line 17) and the fallback ladder (lines 82-87). Must record:
1. **Stack decision (ENV-04 / D-07):** NeuralEstimators.jl = default; BayesFlow (PythonCall) = fallback only, invoked only if the Julia flow won't converge after tuning.
2. **Julia-version record (D-05):** literal `1.12.6` + the `juliaup override set 1.12.6` mechanism note (RESEARCH lines 192-201). Note explicitly that `.julia-version` is documentation-only (NOT auto-read by juliaup) and the *enforcing* pin is the directory override.
3. **Seed for Phase 2 (D-07, CONTEXT line 32):** this file is where Phase 2 will copy the Turing prior ranges. For convenience, the executor MAY pre-seed the exact ranges from `src/bayes.jl` lines 274-291 (e.g. `μ ~ Truncated(Cauchy(0,0.3),-1,1)`, `ν ~ Exponential()`, `σ ~ Truncated(Cauchy(0.1,0.3),0.0001,1)`, `τ ~ Truncated(Cauchy(0.1,0.3),0.0001,1)`) — but that is Phase-2 scope, optional now.

---

### `spike/.julia-version` (config) — NO IN-REPO ANALOG

**Analog:** NONE. Single-line literal file:
```
1.12.6
```
Documentation record only (RESEARCH line 201); the real pin is `juliaup override set 1.12.6` run inside `spike/`. Planner discretion on the exact mechanism (CONTEXT line 37).

---

## Shared Patterns

### License / file header
**Source:** `src/colocalization.jl` lines 1-19 (and identical atop every `src/*.jl`: `ProteinCoLoc.jl`, `LoadImages.jl`, `bayes.jl`).
**Apply to:** `spike/00_smoke.jl` (and optionally `spike/test/runtests.jl`).
**Pattern:** AGPL `#= ... =#` block, author "Dr. rer. nat. Manuel Seefelder". (Note the inconsistency: `test/runtests.jl` lines 1-5 use a shorter `#####`-banner MIT header instead — choose one convention per file; `#=` AGPL is dominant.)

### Deterministic seed at top of executable file
**Source:** `test/runtests.jl` line 12 (`Random123.seed!(1234)`) and the package's stated Random123 reproducibility intent.
**Apply to:** `spike/00_smoke.jl` and `spike/test/runtests.jl`.
**Phase-1 deviation:** use stdlib `Random.seed!(2026)`, NOT `Random123` (Random123 is a later-phase concern — RESEARCH Pitfall 5, line 248).

### `@testset "..." verbose = true begin ... end` wrapper
**Source:** `test/runtests.jl` line 17.
**Apply to:** `spike/test/runtests.jl`.

### Isolated-env discipline (hard constraint)
**Source:** CLAUDE.md §Constraints + RESEARCH Pattern 2 (lines 174-189).
**Apply to:** all `Pkg` operations — always `Pkg.activate("spike")` first; NEVER `Pkg.add`/`resolve` against the root project.

---

## No Analog Found

Files with no close in-repo match — planner must use the external API sequence in `01-RESEARCH.md` (and CLAUDE.md stack notes), NOT a fabricated local pattern:

| File | Role | Data Flow | Reason / Authoritative Source |
|------|------|-----------|-------------------------------|
| `spike/00_smoke.jl` | script | transform | No NeuralEstimators/Flux/NN code exists in repo. Use RESEARCH §Code Examples lines 273-304 + verify `?` docstrings. |
| `spike/Manifest.toml` | config | generated | Machine-emitted by `Pkg.resolve()`; not authored. Commit verbatim (ENV-03). |
| `spike/NOTES.md` | doc | declarative | No decision-doc precedent. Content from CLAUDE.md "What NOT to Use" + RESEARCH D-07. |
| `spike/.julia-version` | config | declarative | No precedent. One literal line `1.12.6`. |

---

## Planning Preconditions (surfaced by RESEARCH — not pattern-derived, but blocking)

These are NOT pattern assignments but the planner must schedule them as they gate ENV-01's "untouched root" proof:

1. **Dirty-root reconciliation (Pitfall 4, RESEARCH lines 241-244 / Q1 lines 397-400):** `git status` CONFIRMED at mapping time shows `M Project.toml` and `M Manifest.toml` — the two files the decoupling guarantee protects are **already modified vs HEAD/`v2.0-baseline`**. The "untouched" proof needs an agreed clean baseline (commit/restore OR snapshot-hash) decided with the user BEFORE any `spike/` work. Tag `v2.0-baseline` exists as a candidate frozen reference (confirmed: `git tag` lists `v1.0.2-scirep`, `v2.0-baseline`).
2. **`Pkg.develop(path="..")` co-resolve (Pitfall 3 / D-01):** discrete gated task AFTER the smoke is green; `include()` fallback pre-authorized. Smoke green-light does not depend on it.

---

## Metadata

**Analog search scope:** repo root (`Project.toml`, `Manifest.toml`), `test/` (`runtests.jl`), `src/` (all 11 `.jl` modules — headers + module structure + Turing `@model` in `bayes.jl`).
**Files scanned:** root `Project.toml`, root `Manifest.toml` (size only), `test/runtests.jl`, `src/ProteinCoLoc.jl`, `src/colocalization.jl`, `src/bayes.jl` (grep), plus dir listings of `src/` and `test/`.
**Sibling `BayesInteractomics`:** referenced in CLAUDE.md for later-phase JLD2 / `_bin_calibration` patterns; NOT reachable in this working tree and NOT needed for Phase 1 — no analog drawn from it.
**Pattern extraction date:** 2026-06-26
