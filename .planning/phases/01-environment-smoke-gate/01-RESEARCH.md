# Phase 1: Environment + Smoke Gate - Research

**Researched:** 2026-06-26
**Domain:** Reproducible isolated Julia environment scaffolding + a correctness-asserting NeuralEstimators.jl/Flux.jl NPE smoke test (Windows, CPU-only, pre-1.0 SBI stack)
**Confidence:** HIGH (NeuralEstimators v0.2.1 API verified against source docstrings on `main`; Julia/juliaup mechanics verified on this machine; decoupling mechanics verified against the actual root `Project.toml`)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** The spike reaches `src/` via `Pkg.develop(path="..")` **first** (clean package boundary). If co-resolving the parent dep tree (GLMakie 0.10.5, Turing, PackageCompiler) with Flux 0.16 conflicts or is unreasonably heavy, **document the conflict and fall back to `include()`** of the specific source files. Either path must leave the main `Project.toml`/`Manifest.toml` and `src/` untouched.
- **D-02:** Rationale for preferring `Pkg.develop`: Phase 4's ADVI baseline needs the *real* Turing `colocalization()` model through the package, so keeping the package boundary now avoids rework later; the `include()` fallback exists only to keep Phase 1 unblocked if the resolve fails.
- **D-03:** The smoke gate is **green + correctness**: it must `train` a `PosteriorEstimator` (NormalisingFlow) on a 1-parameter Gaussian, run `sampleposterior`, **and assert the recovered posterior mean is within a tolerance of the known value** — not merely "runs without error."
- **D-04:** CPU-only, with **no forced CUDA import** (ROADMAP SC2). Core stays close to the ~30-line target; the correctness assertion adds a few lines.
- **D-05:** Pin **both** the package versions (`spike/Manifest.toml`, exact versions captured after the first `Pkg.resolve()`) **and the Julia version** (record exact `julia --version`; add a `.julia-version` / juliaup-channel note). Commit as the reproducibility artifact.
- **D-06:** Wrap the smoke in a **re-runnable harness** `spike/test/runtests.jl` (stdlib `Test`). `spike/00_smoke.jl` remains the readable entry script. **No external CI** in scope.
- **D-07:** Record NeuralEstimators.jl as default and BayesFlow (PythonCall) as fallback-only. Default location `spike/NOTES.md`; a dedicated `spike/STACK-DECISION.md` is acceptable at planner discretion.

### Claude's Discretion
- Exact tolerance value for the correctness assertion.
- File/dir layout within `spike/` (beyond the named `00_smoke.jl`, `test/runtests.jl`, `NOTES.md`, `Project.toml`/`Manifest.toml`).
- Exact juliaup/`.julia-version` mechanism for recording the Julia version.
- `NOTES.md` vs a dedicated stack-decision file for D-07.

### Deferred Ideas (OUT OF SCOPE)
- **CUDA / GPU acceleration** — only revisited if Phase 4 training is slow on CPU; kept out of the Phase-1 env.
- **External CI** for the smoke gate — out of scope; the re-runnable local harness (D-06) is the chosen guard.
- **RxInfer baseline** — v2/deferred (BACK-01), not a spike dependency.
- Simulator, training data, real inference, SBC, OOD, BF — all later phases. Phase 1 uses a synthetic 1-param Gaussian only.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| ENV-01 | `spike/` with own `Project.toml`, isolated (`Pkg.activate`), main package read-only via `Pkg.develop`; root `Project.toml`/`Manifest.toml` and `src/` untouched | Decoupling mechanics + the **dirty-root blocker** (root manifests already modified vs HEAD — see Open Questions Q1); `Pkg.develop(path="..")` vs `include()` fallback verified against actual root deps |
| ENV-02 | <30-line smoke (`spike/00_smoke.jl`) trains a `PosteriorEstimator` (NormalisingFlow) on a 1-param Gaussian, runs `sampleposterior` | Verified v0.2.1 API call sequence (Code Examples §); `use_gpu=false` requirement; posterior-mean tolerance assertion (d×N layout) |
| ENV-03 | Resolved `spike/Manifest.toml` pinned + committed (exact NeuralEstimators/Flux); smoke is the hard gate | Pin-after-resolve workflow; Julia version pin via `juliaup override set 1.12.6` + recorded `.julia-version`; re-runnable `runtests.jl` |
| ENV-04 | Stack-decision note: NeuralEstimators.jl default, BayesFlow (PythonCall) fallback only | Content sourced from CLAUDE.md "What NOT to Use" + STACK.md fallback ladder |
</phase_requirements>

## Summary

Phase 1 is environment scaffolding plus a single high-value de-risking test. There is **no novel algorithm** here; the entire value is proving that the pre-1.0 NeuralEstimators v0.2.1 + Flux v0.16 stack trains a normalising-flow posterior estimator and recovers a known parameter **on this exact Windows machine, CPU-only**, then freezing that environment (package Manifest + Julia version) so every later phase inherits a proven, reproducible base. The top project risk (P4: pre-1.0 API churn; P5: Windows/GPU breakage) is retired here or not at all.

The single most important technical finding is **API signature drift** in NeuralEstimators between v0.2.0 and v0.2.1. Older examples (the CRAN/R vignette, some blog snippets, and the project's own STACK.md notes) show `PosteriorEstimator(q, network)` and a positional `sampleposterior(estimator, Z, N)`. The **current v0.2.1 source** uses `PosteriorEstimator(network, d; num_summaries = d, q = ...)` and `sampleposterior(estimator, Z; N = 1000)` (N is now a keyword). The smoke MUST be written against the installed v0.2.1 docstrings (`?PosteriorEstimator`, `?NormalisingFlow`, `?sampleposterior` in the REPL on the resolved env), not from memory or stale examples. This is exactly the churn the smoke gate exists to catch.

The second load-bearing finding: `train(...)` defaults to **`use_gpu::Bool = true`**. To satisfy D-04 ("no forced CUDA import"), the smoke must pass `use_gpu = false` (or `device = cpu_device()`) explicitly. With Flux as the backend and no CUDA in the Manifest, `use_gpu=false` guarantees the CPU path and prevents any CUDA artifact from being pulled or imported.

**Primary recommendation:** Build the smoke FIRST against a minimal `spike/Project.toml` (NeuralEstimators, Flux, Distributions, Test only) and get it green CPU-only; THEN, as a *separate, independently-verifiable step*, attempt `Pkg.develop(path="..")` and record whether the parent co-resolves with Flux 0.16 — falling back to `include()` per D-01 if it does not. The smoke's success does not depend on the coupling, so the two highest-risk unknowns (does the SBI stack work? does the parent co-resolve?) are isolated rather than entangled.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Isolated dependency resolution | Build / Pkg environment (`spike/Project.toml`+`Manifest.toml`) | — | Decoupling invariant lives in the env boundary, not in code |
| Parent-package read-only access | Pkg environment (`Pkg.develop`) | Source include (`include("../src/...")` fallback) | D-01: package boundary preferred; raw include is the unblock path |
| NPE training + sampling | Library (NeuralEstimators/Flux, CPU) | — | All numerics owned by the SBI stack; spike only orchestrates |
| Correctness assertion | Test harness (stdlib `Test`) | Readable script (`00_smoke.jl`) | D-03/D-06: gate is a re-runnable test; script is the human-readable twin |
| Reproducibility pinning | Build artifacts (`Manifest.toml` + Julia version record) | VCS (git commit) | D-05: env freeze is a committed artifact, not runtime behaviour |
| CPU-only guarantee | Library config (`use_gpu=false`) + env (no CUDA dep) | — | D-04: enforced both by not installing CUDA and by the train flag |

## Standard Stack

### Core (Phase 1 only — the minimal smoke set)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| NeuralEstimators.jl | **v0.2.1** (released 2026-04-07) `[VERIFIED: GitHub releases API]` | `PosteriorEstimator`, `NormalisingFlow`, `train`, `sampleposterior`, `posteriormean` | The chosen SBI engine (CLAUDE.md locked stack); v0.2.1 added Lux backend-agnosticism, Flux still first-class |
| Flux.jl | **v0.16.x** `[CITED: CLAUDE.md stack table]` | NN backend for the flow's coupling layers | Default NeuralEstimators backend; CUDA is a separate extension since v0.14, so CPU installs stay clean |
| Distributions.jl | current `[CITED: CLAUDE.md]` | Define the toy Gaussian prior + likelihood for the smoke | Lingua franca; trivial for a 1-param Gaussian |
| Test (stdlib) | bundled with Julia 1.12.6 | `@test`/`@testset` wrapper for `runtests.jl` (D-06) | Zero-dependency, ships with Julia |

> Keep the **Phase-1** `spike/Project.toml` minimal (the four above). The full spike stack (JLD2, Random123, Images, ImageFiltering, StatsBase, HypothesisTests, KernelDensity, QuadGK, BenchmarkTools, CairoMakie) from STACK.md belongs to **later phases** and should NOT be added in Phase 1 — every extra dep widens the resolve and the Windows-fragility surface the smoke is meant to isolate. Add them per-phase as needed.

### Supporting (decision only — do not install in Phase 1)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| ProteinCoLoc (the parent, via `Pkg.develop`) | local path `..` | Read-only access to `correlation`/`patch`/`_prepare_data` for Phase 2+ | D-01: develop into spike env; not used BY the Phase-1 smoke |
| CUDA.jl | — | Optional GPU accelerator | **Never in Phase 1.** Deferred until Phase 4 training proves too slow |

### Alternatives Considered (locked by CLAUDE.md / CONTEXT — no exploration needed)
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| NeuralEstimators built-in `NormalisingFlow` | NormalizingFlows.jl (TuringLang) | Separate flow stack; explicitly on the "What NOT to Use" list for the spike |
| NeuralEstimators (Julia) | BayesFlow via PythonCall/CondaPkg | **Fallback only** (D-07) — invoke only if the Julia flow won't converge after tuning |
| Flux backend | Lux.jl + Reactant | Reserve for GPU/throughput; not a Phase-1 concern |

**Installation (Phase 1 minimal env):**
```julia
# Run from the repo root; activates spike/ — NEVER touches root Project.toml/Manifest.toml
using Pkg
Pkg.activate("spike")
Pkg.add(["NeuralEstimators", "Flux", "Distributions"])   # Test is stdlib
# (Pkg.develop(path = "..") attempted SEPARATELY — see Architecture Patterns §Pattern 2)
Pkg.resolve()                                            # then commit spike/Manifest.toml
```

**Version verification performed this session:**
- `NeuralEstimators` v0.2.1, published 2026-04-07 — confirmed via GitHub releases API. `[VERIFIED: github releases]`
- Julia `1.12.6` (juliaup `release` channel, default) — confirmed on this machine via `julia --version` / `juliaup status`. `[VERIFIED: local shell]`
- Flux v0.16.x compat with NeuralEstimators v0.2.1 — per CLAUDE.md compat table; **exact transitive versions to be captured from the first `spike/Pkg.resolve()`** and committed (ENV-03). `[CITED: CLAUDE.md]`

## Package Legitimacy Audit

> Julia packages resolve through the **Julia General registry**, not npm/PyPI/crates. `slopcheck` (a PyPI/npm/crates tool) does not apply to the Julia ecosystem and was not run. Legitimacy was instead verified via the General registry path (official maintainer GitHub repo + release history + a maintained companion CRAN package).

| Package | Registry | Age / Signal | Source Repo | Verdict | Disposition |
|---------|----------|--------------|-------------|---------|-------------|
| NeuralEstimators.jl | Julia General | 7 releases, 605+ commits, v0.2.1 (2026-04-07); companion CRAN R pkg | github.com/msainsburydale/NeuralEstimators.jl | Legitimate, but **pre-1.0 / API-churning** | Approved — pin exact version in Manifest |
| Flux.jl | Julia General | FluxML org; ubiquitous, mature | github.com/FluxML/Flux.jl | Legitimate | Approved |
| Distributions.jl | Julia General | JuliaStats org; ubiquitous | github.com/JuliaStats/Distributions.jl | Legitimate | Approved |

**Packages removed (slop):** none.
**Packages flagged suspicious:** none (NeuralEstimators is legitimate; its risk is *version churn*, not provenance — mitigated by the pinned Manifest, ENV-03).

## Architecture Patterns

### System Architecture Diagram (Phase 1 data + control flow)

```
                       spike/Project.toml  (activate; isolated env)
                                │
                 Pkg.add(NeuralEstimators, Flux, Distributions)
                                │
                          Pkg.resolve()  ──────────►  spike/Manifest.toml  (PIN + commit, ENV-03)
                                │                              ▲
                                ▼                              │ exact versions
   ┌──────────────────────  00_smoke.jl  (readable script, ENV-02) ─────────────┐
   │  d = 1                                                                       │
   │  prior θ ~ N(0,1)         ──► sampler(K)  -> 1×K parameter matrix            │
   │  likelihood Z|θ ~ N(θ,σ²) ──► simulator(θ) -> data conditioned on θ          │
   │  net  = summary/inner Flux Chain (output dim = num_summaries)                │
   │  q    = NormalisingFlow(d=1, num_summaries)                                  │
   │  est  = PosteriorEstimator(net, d=1; num_summaries, q)                       │
   │  est  = train(est, sampler, simulator; epochs=…, use_gpu=FALSE)  ◄── CPU-only│
   │  θ̂s   = sampleposterior(est, Z_obs; N=…)   -> d×N matrix                      │
   │  μ̂    = posteriormean(θ̂s)  == mean(θ̂s; dims=2)                                │
   │  @test abs(μ̂ - θ_true) < tol           ◄── correctness assertion (D-03)      │
   └──────────────────────────────────────────────────────────────────────────────┘
                                │  same logic, wrapped in @testset
                                ▼
                spike/test/runtests.jl  (re-runnable hard gate, D-06)
                                │
                   julia --project=spike spike/test/runtests.jl  → GREEN

   Decoupling boundary (verified separately, ENV-01):
       Pkg.develop(path="..")  ──► parent ProteinCoLoc available read-only
            └─ if resolve conflicts with Flux 0.16 → fall back to include("../src/…")
       root Project.toml / Manifest.toml / src/  ── PROVABLY UNTOUCHED  (git status clean)
```

### Recommended Project Structure (Phase 1 deliverables only)
```
ProteinCoLoc/
├── src/                 # UNTOUCHED (read-only)
├── Project.toml         # root — UNTOUCHED (see Open Questions Q1: currently dirty!)
├── Manifest.toml        # root — UNTOUCHED (see Open Questions Q1: currently dirty!)
└── spike/
    ├── Project.toml      # minimal: NeuralEstimators, Flux, Distributions (+ parent via develop)
    ├── Manifest.toml      # PINNED + committed reproducibility artifact (ENV-03)
    ├── .julia-version     # records "1.12.6" (documentation pin; see D-05 mechanism below)
    ├── 00_smoke.jl        # readable smoke script (ENV-02)
    ├── NOTES.md           # stack-decision note (ENV-04 / D-07) + Julia-version record
    └── test/
        └── runtests.jl    # stdlib Test wrapper of the smoke (hard gate, D-06)
```

### Pattern 1: Verify-against-installed-docstrings, not memory (pre-1.0 churn discipline)
**What:** Before writing the smoke, the implementer opens the REPL on the resolved spike env and reads `?PosteriorEstimator`, `?NormalisingFlow`, `?train`, `?sampleposterior`. The smoke is written from those docstrings.
**When to use:** Always, for any pre-1.0 dependency. NeuralEstimators is v0.2.x; signatures demonstrably changed between v0.2.0 and v0.2.1.
**Why:** Catches the exact failure mode (P4) where a copy-pasted example targets a different version and silently mis-constructs the estimator.

### Pattern 2: Decouple the smoke from the coupling (build order)
**What:** Stand up the minimal env and get the smoke green *before* attempting `Pkg.develop(path="..")`. The smoke uses a synthetic Gaussian and needs **nothing** from `src/`.
**When to use:** Phase 1 specifically — two independent risks (SBI stack works / parent co-resolves) should not be entangled.
**Example sequence:**
```julia
# Step A — minimal env, prove the stack (no parent yet)
Pkg.activate("spike"); Pkg.add(["NeuralEstimators","Flux","Distributions"]); Pkg.resolve()
# → write & run 00_smoke.jl → GREEN  → commit spike/Manifest.toml

# Step B — attempt the clean package boundary (D-01)
Pkg.develop(path = "..")     # brings ProteinCoLoc deps (Turing, GLMakie 0.10.5, PackageCompiler, Images…)
Pkg.resolve()                # ← the risk point: does GLMakie 0.10.5 co-resolve with Flux 0.16?
#   success → re-run smoke (must stay green), commit updated Manifest
#   conflict/too-heavy → DOCUMENT it, `Pkg.rm` the develop, fall back to include("../src/colocalization.jl")
```
**Why:** `Pkg.develop` pulls the parent's full `[deps]` as resolve constraints. The root pins **GLMakie = "0.10.5"** (an old, tightly-bound Makie-ecosystem version) alongside Turing and PackageCompiler. Co-resolving that with Flux 0.16 may conflict on shared transitive deps. The fallback (`include`) sidesteps the resolve entirely and keeps Phase 1 unblocked, at the cost of re-doing the boundary in Phase 4 (acceptable per D-02).

### Pattern 3: Reproducibility freeze as a committed artifact (ENV-03 / D-05)
**What:** After the smoke is green, commit `spike/Manifest.toml` (exact NeuralEstimators/Flux + every transitive version) AND record the Julia version. Pin the Julia version per-directory with juliaup.
**Mechanism (verified on this machine):**
```bash
juliaup status                     # → release = 1.12.6 (default)
juliaup override set 1.12.6        # run inside spike/ : directory-scoped Julia pin (juliaup feature, confirmed)
# also write the literal version into a committed file for humans + non-juliaup users:
#   spike/.julia-version  ->  "1.12.6"
#   and a line in spike/NOTES.md
```
**Note:** `.julia-version` is **not** auto-read by Julia/juliaup (that is an rbenv/pyenv convention, not a Julia one). Use it purely as a committed documentation record; the *enforcing* mechanism is `juliaup override set` (directory override) — both together satisfy D-05.

### Anti-Patterns to Avoid
- **Writing the simulator/data pipeline before the smoke is green** (P4/Anti-Pattern 4): the smoke is the gate; nothing downstream is built until it passes.
- **Adding the full spike dependency list in Phase 1**: widens the resolve and the Windows-fragility surface unnecessarily. Add per-phase.
- **Trusting `use_gpu`'s default**: it defaults to `true`; on a CUDA-less machine it falls back, but explicitly passing `use_gpu=false` is the contract that satisfies D-04 and prevents accidental CUDA artifact pulls.
- **Editing or "lightly refactoring" `src/` or the root manifests** to make coupling easier (Anti-Pattern 1): violates the decoupling invariant and the publication-integrity guarantee.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Normalising-flow posterior | A custom flow / coupling layers | `NormalisingFlow(d, num_summaries)` | Built into NeuralEstimators; the whole point of the chosen stack |
| Training loop / early stopping / LR schedule | A manual Flux training loop | `train(estimator, sampler, simulator; epochs, stopping_epochs, ...)` | Handles batching, validation risk, CosAnneal schedule, device placement |
| Posterior summary | Manual mean over samples | `posteriormean(θ̂)` (`= mean(θ̂; dims=2)`) | Correct d×N layout handling; matches package conventions |
| Julia version pinning | A custom shell wrapper | `juliaup override set` | Native juliaup directory-override feature |
| Test harness | A bespoke pass/fail script | stdlib `Test` (`@testset`/`@test`) | Zero-dep, re-runnable, the D-06 requirement |

**Key insight:** Phase 1 writes almost no logic of its own — it *wires* a proven library and *freezes* the result. The only first-party code is ~30 lines of glue (prior, likelihood, network shape) and one assertion.

## Common Pitfalls

### Pitfall 1: API signature drift (NeuralEstimators v0.2.0 → v0.2.1)
**What goes wrong:** Copy a `PosteriorEstimator(q, network)` / positional `sampleposterior(est, Z, N)` example (CRAN R vignette, STACK.md notes) and it errors or mis-binds on v0.2.1.
**Why:** v0.2.1 uses `PosteriorEstimator(network, d; num_summaries=d, q=...)` and `sampleposterior(est, Z; N=1000)` (N is now keyword). `[VERIFIED: NeuralEstimators main source — inference.jl, dev quick-start]`
**How to avoid:** Read `?PosteriorEstimator`/`?NormalisingFlow`/`?sampleposterior` on the resolved env and write from those. Pin the Manifest so the docstrings can't shift under you.
**Warning signs:** `MethodError` on the constructor; `sampleposterior` returning an unexpected shape; examples that only run after un-pinning.

### Pitfall 2: GPU default pulls/forces CUDA (violates D-04)
**What goes wrong:** `train` defaults to `use_gpu=true`; a stray CUDA artifact or import sneaks onto the "CPU-only" path.
**Why:** `train(...; use_gpu::Bool = true)` `[VERIFIED: train.jl source]`.
**How to avoid:** Pass `use_gpu=false` (or `device = cpu_device()`); never `Pkg.add("CUDA")` in the Phase-1 env; assert in the harness that no CUDA package is loaded (e.g. check `CUDA` is not a key in the active project, or that the Manifest has no CUDA entry).
**Warning signs:** Manifest contains `CUDA`/`CUDA_jll`; runtime asks for a CUDA driver during a CPU run.

### Pitfall 3: `Pkg.develop` co-resolve conflict (GLMakie 0.10.5 + Flux 0.16)
**What goes wrong:** `Pkg.develop(path="..")` + `Pkg.resolve()` fails or drags a huge, fragile tree into the spike env.
**Why:** The root pins GLMakie 0.10.5 (old Makie ecosystem) with Turing + PackageCompiler; shared transitive deps may be incompatible with Flux 0.16.
**How to avoid:** Build the smoke first WITHOUT the parent (Pattern 2). Attempt develop separately; on conflict, document and fall back to `include()` (D-01). The Phase-1 smoke needs nothing from `src/`, so a failed develop does not block the gate.
**Warning signs:** `Unsatisfiable requirements` from the resolver mentioning Makie/Colors/StaticArrays; resolve pulling hundreds of packages.

### Pitfall 4: The "untouched root" claim against an already-dirty tree
**What goes wrong:** ENV-01 requires the root `Project.toml`/`Manifest.toml` untouched, but they are **already modified vs HEAD** in the current working tree (see Open Questions Q1). A naive `git status`-clean proof will fail for reasons unrelated to the spike.
**How to avoid:** Establish a clean baseline reference BEFORE creating `spike/` — either commit/restore the root files, or snapshot their current content hash and assert byte-identity against that snapshot after the spike work. Decide this with the user up front.
**Warning signs:** `git status` shows `M Project.toml` / `M Manifest.toml` before any spike work is done (it does, right now).

### Pitfall 5: Non-deterministic smoke flakes the gate
**What goes wrong:** A tolerance assertion on a stochastically-trained flow occasionally fails, undermining "hard gate."
**How to avoid:** Seed the RNG (Julia `Random.seed!` is sufficient for Phase 1; Random123 is a later-phase concern), train enough epochs for the 1-param Gaussian to converge comfortably, and choose a **generous tolerance** (the assertion's job is to catch a silently-wrong backend, not to certify accuracy). A 1-param Gaussian posterior mean within, e.g., a few ×0.1 of truth is ample; tune so a correct run passes with margin and a broken backend fails clearly.

## Code Examples

### Verified v0.2.1 API surface (from NeuralEstimators `main` source)
```julia
# Source: github.com/msainsburydale/NeuralEstimators.jl  src/ApproximateDistributions/NormalisingFlow.jl
#   NormalisingFlow(d::Integer, num_summaries::Integer; num_coupling_layers = 6, use_act_norm = true, backend = nothing)
#     d            = dimension of the parameter vector (1 for the smoke)
#     num_summaries = dimension of the summary statistics = output width of the neural network
#
# Source: dev docs quick-start (v0.2.1)
#   PosteriorEstimator(network, d; num_summaries = d, q = GaussianMixture)   # q may be a NormalisingFlow instance
#
# Source: src/inference.jl (verbatim docstring)
#   sampleposterior(estimator::PosteriorEstimator, Z; N::Integer = 1000, kwargs...)
#     → single data set: returns a  d × N  matrix (parameters in ROWS, draws in COLUMNS)
#   posteriormean(θ::AbstractMatrix) = mean(θ; dims = 2)        # → d × 1
#
# Source: src/train.jl (verbatim signature + kwargs)
#   train(estimator, sampler, simulator; epochs = 100, stopping_epochs = 5,
#         batchsize = 32, optimiser = Adam(5e-4), use_gpu = true, device = nothing, ...)
#   ⚠ use_gpu defaults to TRUE — pass use_gpu = false for the CPU-only smoke (D-04)
```

### Smoke skeleton (illustrative — MUST be re-verified against installed `?` docstrings)
```julia
# spike/00_smoke.jl  — 1-parameter Gaussian NPE smoke (ENV-02), CPU-only (D-04)
using NeuralEstimators, Flux, Distributions, Random
Random.seed!(2026)

d = 1                                   # infer the mean θ
σ = 1.0f0                               # known noise
m = 1                                   # iid replicates per data set (1 for the smoke)

# prior sampler: returns a d×K matrix of θ draws  (θ ~ N(0,1))
sample(K) = Float32.(randn(1, K))
# simulator: given a d×K parameter matrix, return K data sets Z|θ ~ N(θ, σ²)
simulate(θ) = [Float32.(θ[1, k] .+ σ .* randn(1, m)) for k in 1:size(θ, 2)]

# summary/inner network → output width = num_summaries
num_summaries = 16
network = Chain(Dense(1 => 32, relu), Dense(32 => num_summaries))
q   = NormalisingFlow(d, num_summaries)
est = PosteriorEstimator(network, d; num_summaries = num_summaries, q = q)

est = train(est, sample, simulate; epochs = 200, use_gpu = false)   # CPU-only

# inference on data generated from a KNOWN θ
θ_true = 0.7f0
Z_obs  = Float32.(θ_true .+ σ .* randn(1, m))
θ̂      = sampleposterior(est, Z_obs; N = 1000)   # d×N matrix
μ̂      = posteriormean(θ̂)[1]                      # mean over columns

tol = 0.3                                          # generous; catches a broken backend (Claude's discretion, D-03)
@assert abs(μ̂ - θ_true) < tol  "smoke FAILED: recovered $μ̂ vs true $θ_true"
```
```julia
# spike/test/runtests.jl  — re-runnable hard gate (D-06)
using Test
@testset "NeuralEstimators CPU smoke" begin
    include(joinpath(@__DIR__, "..", "00_smoke.jl"))   # or refactor the smoke into a function and call it
    @test abs(μ̂ - θ_true) < tol
    @test !haskey(Pkg.project().dependencies, "CUDA")  # CPU-only guard (D-04); adjust to chosen check
end
```
> The exact `sample`/`simulate` return shapes, the `m` keyword on `train`, and whether `q` is passed positionally or as the `q=` keyword **must be confirmed against the installed v0.2.1 docstrings** — these are the precise spots where pre-1.0 drift bites. Treat the skeleton as structure, not gospel.

## Runtime State Inventory

> Phase 1 creates a new isolated environment; it is greenfield within `spike/`. The only "existing runtime state" that matters is whatever must stay UNTOUCHED.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — no databases or datastores involved in env scaffolding | None |
| Live service config | None — no external services | None |
| OS-registered state | **juliaup directory override** for `spike/` will be registered in juliaup's config (`juliaup override set 1.12.6`) — this is intended, not incidental | Record the override in `spike/NOTES.md` so it is reproducible on another machine |
| Secrets/env vars | None | None |
| Build artifacts | **Root `Project.toml` + `Manifest.toml` are currently modified vs HEAD** (working-tree dirty — full Manifest re-resolve + a `[compat] GLMakie = "0.10.5"` addition to Project.toml). These are the artifacts the "untouched" proof depends on. | **Reconcile the dirty root state BEFORE spike work** (commit/restore or snapshot-hash baseline) — see Open Questions Q1. `spike/artifacts/` (later phases) is already covered by `.gitignore`'s `artifacts/` rule. |

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (ships with Julia 1.12.6) |
| Config file | none — `spike/test/runtests.jl` (created Phase 1, D-06) |
| Quick run command | `julia --project=spike spike/test/runtests.jl` |
| Full suite command | `julia --project=spike -e 'using Pkg; Pkg.test()'` (once `spike/Project.toml` declares the package context) — or just the runtests invocation above for Phase 1 |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| ENV-02 | NPE trains + samples + recovers θ within tol | smoke/unit | `julia --project=spike spike/test/runtests.jl` | ❌ Wave 0 — create `00_smoke.jl` + `test/runtests.jl` |
| ENV-04 / D-04 | No CUDA on the CPU path | unit (guard `@test`) | same harness | ❌ Wave 0 — add CUDA-absence assertion |
| ENV-01 | Root `Project.toml`/`Manifest.toml`/`src/` untouched | integration (git check) | `git status --porcelain src/ Project.toml Manifest.toml` clean vs agreed baseline | ❌ Wave 0 — needs baseline reconciliation (Q1) |
| ENV-03 | `spike/Manifest.toml` pinned + committed; Julia version recorded | manual/artifact | inspect committed `spike/Manifest.toml` + `spike/.julia-version` | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** `julia --project=spike spike/test/runtests.jl` (the smoke is fast — seconds to low minutes on CPU for a 1-param Gaussian).
- **Per wave merge:** same (single gate).
- **Phase gate:** green `runtests.jl` + clean `git status` on root files vs the agreed baseline + committed `spike/Manifest.toml`.

### Wave 0 Gaps
- [ ] `spike/00_smoke.jl` — covers ENV-02 (the readable script).
- [ ] `spike/test/runtests.jl` — covers ENV-02 + D-04 CUDA-absence guard (the hard gate).
- [ ] `spike/Project.toml` + resolved/committed `spike/Manifest.toml` — covers ENV-03.
- [ ] `spike/.julia-version` + `juliaup override set 1.12.6` + NOTES.md record — covers ENV-03 Julia pin.
- [ ] `spike/NOTES.md` (or `STACK-DECISION.md`) — covers ENV-04.
- [ ] Root-baseline reconciliation step (Q1) — precondition for ENV-01's untouched proof.

## Security Domain

> This phase is offline, simulation-only environment scaffolding: no authentication, no network service, no user/untrusted input, no PII. The ASVS application-security categories (V2 auth, V3 session, V4 access control, V5 input validation, V6 crypto) **do not apply**. `[ASSUMED — based on phase scope]`

The one genuinely relevant integrity axis is **software supply chain / reproducibility**, and it is already addressed by the phase's own requirements:

| Concern | STRIDE | Mitigation (already in scope) |
|---------|--------|-------------------------------|
| Dependency substitution / silent version drift | Tampering | Pinned, committed `spike/Manifest.toml` (ENV-03); legitimacy audit above; Julia General registry provenance |
| Irreproducible environment | Repudiation | Julia version pin (`juliaup override`) + recorded `.julia-version` (D-05) |
| Accidental CUDA artifact on the "clean" path | Tampering | `use_gpu=false` + no CUDA dep + harness guard (D-04) |

No separate security tasks are needed beyond the existing ENV-03/D-04 requirements.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `PosteriorEstimator(q, network)` (v0.2.0 / CRAN vignette) | `PosteriorEstimator(network, d; num_summaries=d, q=...)` (v0.2.1) | v0.2.1 (2026-04-07) | Stale examples error; write from installed docstrings |
| `sampleposterior(est, Z, N)` positional N (STACK.md note) | `sampleposterior(est, Z; N=1000)` keyword N | by v0.2.1 | Positional N call fails |
| Flux-only NeuralEstimators | Backend-agnostic (Flux **or** Lux), `use_gpu`/`device` kwargs, Reactant device option | v0.2.1 (PR #76) | Flux still works and is the chosen backend; just confirm CPU device handling |

**Deprecated/outdated for this phase:**
- The STACK.md "Key API Notes" `train(estimator, sampler, simulator; ...)` shape is broadly correct, but its `sampleposterior(estimator, Z, N)` is stale (N is now keyword) — minor but real.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Flux v0.16.x is the version that co-resolves with NeuralEstimators v0.2.1 (exact transitive set unknown until first `Pkg.resolve()`) | Standard Stack | Low — resolve will pick a compatible set; just capture & pin it |
| A2 | `juliaup override set 1.12.6` is the correct directory-pin mechanism (verified `juliaup override` exists with set/unset/status; exact behavior assumed standard) | Pattern 3 / D-05 | Low — fallback is the committed `.julia-version` + NOTES record |
| A3 | The `q=` argument of `PosteriorEstimator` accepts a constructed `NormalisingFlow` instance (dev quick-start shows `q = GaussianMixture` as a type; instance-passing inferred from older `PosteriorEstimator(q, network)` form) | Code Examples | Medium — **verify against installed `?PosteriorEstimator`**; this is the single most likely smoke-construction snag |
| A4 | Security/ASVS categories don't apply (offline sim-only scaffolding) | Security Domain | Low — phase has no network/auth/input surface |
| A5 | A generous tolerance (~0.3 on a unit-scale 1-param Gaussian) reliably separates a correct backend from a broken one | Pitfall 5 / Code Examples | Low — tune empirically; widen if a correct run is marginal |
| A6 | The root `Project.toml`/`Manifest.toml` dirty state is incidental (not a deliberate uncommitted spike pre-step) and should be reconciled | Open Questions Q1 | Medium — if the user intended those changes, the baseline definition changes; must confirm with user |

## Open Questions (RESOLVED)

1. **The root `Project.toml`/`Manifest.toml` are already modified vs HEAD (`v2.0-baseline`).** `Project.toml` gained a `[compat] GLMakie = "0.10.5"` block and `Manifest.toml` is a full re-resolve (1162 insert / 1073 delete). ENV-01 / DEMO-02's "untouched" proof needs a clean reference.
   - What we know: working tree is dirty on exactly the two files the decoupling guarantee protects; tags `v1.0.2-scirep` and `v2.0-baseline` exist as frozen references.
   - What's unclear: whether these edits are intentional (and should be committed as the new baseline) or accidental (and should be `git restore`d).
   - **Recommendation:** Before any `spike/` work, the plan must include a task to reconcile this with the user — commit or restore the root files — and record the agreed baseline (tag or commit hash). The untouched-proof then asserts byte-identity against that reference. Flag this as a **planning precondition**, not a mid-phase surprise.
   - **RESOLVED (Plan 01-01):** Reconciliation is a discrete blocking decision checkpoint (commit-as-baseline | restore-to-head | adopt-dirty-snapshot); the chosen baseline ref (commit short-hash, HEAD, or per-file `git hash-object` hashes) is recorded in 01-BASELINE.md, and every later untouched-root proof asserts state-matches-baseline against it.

2. **Does `Pkg.develop(path="..")` co-resolve with Flux 0.16?** Unknown until run. Recommendation: make it a discrete, gated task with the `include()` fallback pre-authorized (D-01); the smoke green-light does not depend on its outcome.
   - **RESOLVED (Plan 01-04 T1):** `Pkg.develop(path="..")` co-resolve is a discrete gated task with the `include()` fallback pre-authorized per D-01; the smoke gate (Plan 01-02) does not depend on its outcome.

3. **Exact `m`/replicate handling and `q` passing in v0.2.1.** Resolve by reading installed docstrings during implementation (Pattern 1). Low risk to planning; the planner should allocate a "verify API against REPL docstrings" step inside the ENV-02 task.
   - **RESOLVED (Plan 01-02 T2):** before writing the smoke, a docstring-verification step opens a REPL on the installed env and reads `?PosteriorEstimator` / `?NormalisingFlow` / `?sampleposterior` to confirm the exact v0.2.1 signatures (PosteriorEstimator positional order, N-as-keyword, q-as-instance, and m/replicate handling).

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Julia | Everything | ✓ | 1.12.6 (juliaup `release`, default) | — |
| juliaup | Julia version pin (D-05) | ✓ | present; `override` subcommand confirmed | committed `.julia-version` record only |
| git | Decoupling proof + commits | ✓ | repo active (tags present) | — |
| NeuralEstimators.jl | ENV-02 | ✗ (not yet installed) | target v0.2.1 | BayesFlow/PythonCall fallback (D-07) — never in Phase 1 |
| Flux.jl | ENV-02 | ✗ (not yet installed) | target v0.16.x | Lux backend (deferred) |
| CUDA / GPU | — | not required | — | CPU-only is the baseline (D-04) |
| Internet (Pkg registry) | First `Pkg.add`/`resolve` | assumed ✓ | — | none — needed once to populate the spike Manifest |

**Missing dependencies with no fallback:** none that block planning (NeuralEstimators/Flux are installed *by* the phase; that is the work, not a blocker).
**Missing dependencies with fallback:** SBI stack convergence failure → BayesFlow (D-07), explicitly out of Phase-1 scope unless the smoke can't be made green after tuning.

## Sources

### Primary (HIGH confidence)
- github.com/msainsburydale/NeuralEstimators.jl — `src/inference.jl` (`sampleposterior` d×N layout, `posteriormean = mean(;dims=2)`), `src/train.jl` (train signatures + `use_gpu=true` default, epochs/stopping_epochs/batchsize/optimiser kwargs), `src/ApproximateDistributions/NormalisingFlow.jl` (`NormalisingFlow(d, num_summaries; num_coupling_layers=6, ...)`). Read verbatim this session.
- GitHub releases API — NeuralEstimators v0.2.1 published 2026-04-07 (PR #76: Lux backend-agnosticism).
- Local shell — `julia --version` = 1.12.6; `juliaup status`/`juliaup override --help`; root `Project.toml` deps + GLMakie 0.10.5 pin; `git tag`/`git diff` confirming the dirty root manifests.
- msainsburydale.github.io/NeuralEstimators.jl/dev/ — quick-start `PosteriorEstimator(network, d; num_summaries=d, q=GaussianMixture)`, `train(estimator, sampler, simulator)`, `sampleposterior(estimator, Z)`.

### Secondary (MEDIUM confidence)
- CLAUDE.md project Technology Stack table + "What NOT to Use" + Version Compatibility (Flux v0.16.x ↔ NeuralEstimators v0.2.1; CUDA-as-extension since Flux v0.14).
- `.planning/research/STACK.md`, `ARCHITECTURE.md`, `PITFALLS.md` (P4 churn, P5 Windows/GPU) — project research inputs.
- CRAN NeuralEstimators vignette / RDocumentation v0.2.0 — source of the *older* `PosteriorEstimator(q, network)` / `NormalisingFlow(d,d)` form (used here to identify the drift, not as current truth).

### Tertiary (LOW confidence)
- General Julia-ecosystem knowledge that `.julia-version` is not auto-read by juliaup (treated as documentation-only; the enforcing pin is `juliaup override`).

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — versions and the minimal set verified against source + local machine.
- Architecture / API call sequence: HIGH on shapes (read from source), MEDIUM on the exact `q`-passing form and `m` handling (flagged A3 — verify against installed docstrings).
- Pitfalls: HIGH — drift and `use_gpu` default verified in source; dirty-root blocker verified via git.
- Decoupling co-resolve outcome: MEDIUM — cannot be known until `Pkg.resolve()` runs; deliberately gated.

**Research date:** 2026-06-26
**Valid until:** ~2026-07-10 (pre-1.0 NeuralEstimators; re-verify API if the spike Manifest is ever re-resolved to a newer version)
