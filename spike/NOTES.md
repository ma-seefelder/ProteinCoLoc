# Spike Notes — ProteinCoLoc v2.0 (AmortizedColoc)

Reproducibility + stack-decision record for the decoupled `spike/` environment.
Created by Phase 1 (environment-smoke-gate), Plan 01-03. Satisfies ENV-03 (Julia-version
record half) and ENV-04 (stack decision).

---

## 1. Stack Decision (ENV-04 / D-07)

**Default SBI engine: `NeuralEstimators.jl` (v0.2.1).**
The amortized inference stack for this spike — and for v2.0 productionization if the spike
succeeds — is NeuralEstimators.jl on the **Flux.jl (v0.16.10)** backend, CPU-only. It provides
both the amortized posterior estimator (`PosteriorEstimator` + built-in `NormalisingFlow`) and
the amortized ratio/Bayes-factor estimator (`RatioEstimator`) from one API. This is the proven,
green-smoke path (Plan 01-02).

**Fallback only: `BayesFlow` (Python, via `PythonCall`/`CondaPkg`).**
BayesFlow is the documented **fallback-only** SBI engine. It is invoked **solely** if the
Julia normalising-flow fails to converge or calibrate after genuine tuning effort. It is NOT
on the core path: it carries two-language complexity, Conda/Windows friction, and a dual-language
reproducibility/seed burden. Do not reach for it unless the Julia flow is demonstrably stuck.

### Explicitly excluded (CLAUDE.md "What NOT to Use")

These packages are **deliberately excluded** for the spike — adding them is a regression:

| Excluded | Reason | Use instead |
|----------|--------|-------------|
| `NormalizingFlows.jl` (TuringLang) | Separate flow stack; duplicates NeuralEstimators' built-in `NormalisingFlow`; GPU-compat caveats. Zero spike benefit. | NeuralEstimators' built-in `NormalisingFlow` |
| `InvertibleNetworks.jl` | Solves image-sized (>1024²) flow scalability we do not have on an 8×8 / 16-summary vector; pure overhead now. | NeuralEstimators flow (revisit only for 3D/hierarchy build-out) |
| `RxInfer.jl` *as a Turing replacement* | Non-conjugate Student-t needs reformulation; optimizes the per-dataset axis that amortized SBI is built to retire. | Keep Turing ADVI as the accuracy baseline; defer RxInfer (BACK-01) |
| `BayesFlow` *on the core path* | Two-language complexity, Conda/Windows friction, reproducibility burden. | Julia-native NeuralEstimators; BayesFlow strictly last-resort fallback (above) |
| `CUDA.jl` *as a requirement* | Windows Flux/CUDA friction is the stated risk; 2D summary-vector training is small. CPU-only is the portable baseline. | CPU-only baseline; CUDA.jl only as an optional Flux extension that degrades gracefully (deferred to Phase 4 if training is slow) |

**Naming gotcha:** NeuralEstimators' internal flow type is `NormalisingFlow` (British spelling),
distinct from the standalone `NormalizingFlows.jl` package above. Do not conflate them.

---

## 2. Julia-Version Record (ENV-03 / D-05)

**Pinned Julia version: `1.12.6`.**

- **Enforcing mechanism (the real pin):** a juliaup **directory override** scoped to `spike/`:
  ```
  juliaup add 1.12.6
  juliaup override set --path <repo>/spike 1.12.6
  juliaup override status      # → <repo>/spike  →  1.12.6
  ```
  This is the mechanism that actually selects Julia 1.12.6 when working inside `spike/`.
  Note: `1.12.6` was added as an explicit juliaup channel (the machine default is the floating
  `release` channel, which currently resolves to 1.12.6 but may move). The override pins the
  exact version, not the floating channel.

- **Documentation record (NOT enforcing):** `spike/.julia-version` contains the single literal
  line `1.12.6`. This file is **documentation-only** — Julia/juliaup do **not** auto-read it
  (that is an rbenv/pyenv convention, not a Julia one). It exists so humans and non-juliaup
  users can see the intended version. The enforcing pin is the juliaup directory override above.

- **To reproduce on another machine:** install juliaup, run `juliaup add 1.12.6`, then
  `juliaup override set --path <clone>/spike 1.12.6`. Then instantiate the pinned environment:
  `julia --project=spike -e 'using Pkg; Pkg.instantiate()'`.

### Package pin (ENV-03 / D-05)

The exact package versions are frozen in the committed `spike/Manifest.toml` (the reproducibility
artifact). Key pins captured after the first green resolve:

| Package | Version |
|---------|---------|
| NeuralEstimators | 0.2.1 |
| Flux | 0.16.10 |
| Distributions | 0.25.128 |
| (Manifest `julia_version`) | 1.12.6 |

**CPU-only guarantee:** the Manifest contains **no top-level `[[deps.CUDA]]` installed package**.
CUDA appears only as inert *weakdep / extension* declarations under Flux / NNlib / Zygote /
NeuralEstimators (e.g. `FluxCUDAExt = "CUDA"`); these are not installed and load nothing on the
CPU path. The CPU-only contract is enforced both by this absence and by `use_gpu = false` at the
`train` call (D-04).

---

## 4. Parent-Package Coupling Outcome (ENV-01 / D-01) — Plan 01-04

**Decision: `include()` fallback chosen.** `Pkg.develop(path="..")` was attempted first
(D-01 develop-first), co-resolved a parent dep tree that is **incompatible with the pinned
NeuralEstimators v0.2.1**, and **broke the green smoke**. Per D-01 the develop was reverted and
the `include()` fallback is adopted. The root `Project.toml`/`Manifest.toml`/`src/` were never
edited; only `spike/Project.toml` + `spike/Manifest.toml` were touched and then reverted to the
committed Plan 01-03 known-good state.

### Verbatim conflict evidence

`Pkg.develop(path="..")` + `Pkg.resolve()` did **not** raise an exception, but it silently
**downgraded** NeuralEstimators to satisfy the parent's heavy tree (Turing 0.44.5, GLMakie
0.10.18 → Makie 0.21, Images 0.26, GraphNeuralNetworks 1.1.0, PackageCompiler 2.4.0):

```
⌃ [38f6df31] ↓ NeuralEstimators v0.2.1 ⇒ v0.1.4
  [12345678] + ProteinCoLoc v1.0.1 `...\ProteinCoLoc`
```

`Pkg.status --outdated -m` confirms 0.2.1 is unreachable alongside the parent tree:

```
⌃ [38f6df31] NeuralEstimators v0.1.4 (<v0.2.1)
⌅ [09ab397b] StructArrays v0.6.21 (<v0.7.3): GeometryBasics, ShaderAbstractions
```

Re-running the smoke against the developed env then **errored** — v0.1.4 ships the OLD
ApproximateDistributions API, incompatible with the v0.2.1-written smoke:

```
NeuralEstimators CPU smoke: Error During Test
  Got exception outside of a @test
  LoadError: MethodError: no method matching NormalisingFlow(::Int64; num_summaries::Int64)
  Closest candidates are:
    NormalisingFlow(::D, ::Vector{<:NeuralEstimators.CouplingLayer})   # v0.1.4 form
  at spike/00_smoke.jl:62
Test Summary: NeuralEstimators CPU smoke | Error 1 Total 1   → exit 1
```

**Root cause:** the parent package transitively constrains NeuralEstimators to `<v0.2.1`
(the GLMakie 0.10.x / Makie 0.21 / GraphNeuralNetworks 1.1.0 ecosystem the parent pins does
not co-resolve with NeuralEstimators 0.2.1 / Flux 0.16 / StructArrays 0.7). This is exactly
Pitfall 3 from 01-RESEARCH.md. Forcing the develop would mean abandoning the pinned, proven
v0.2.1 stack — unacceptable — so the package boundary is deferred (D-02: re-evaluate at Phase 4,
which may use a separate parent-aware env for the ADVI baseline rather than dragging the parent
into the lean SBI env).

### Resolution: revert + `include()` fallback

1. `git checkout -- spike/Project.toml spike/Manifest.toml` → restored the Plan 01-03
   minimal env (NeuralEstimators 0.2.1, Flux 0.16.10, Distributions 0.25.128, no parent).
2. `Pkg.instantiate()` + re-ran the smoke → **GREEN** (4/4 pass,
   recovered `mu_hat ≈ 0.799` vs `theta_true = 0.7`, |delta| < tol = 0.3).
3. **`include()` fallback (the chosen coupling path):** later phases reach the reusable `src/`
   assets by `include`-ing the specific source files into the spike env rather than via the
   package boundary, e.g.:
   ```julia
   # Phase 2+ (NOT Phase 1): bring in the summary-statistic contract read-only.
   include(joinpath(@__DIR__, "..", "src", "colocalization.jl"))  # correlation, patch, _prepare_data
   ```
   - `src/colocalization.jl` (`correlation`, `patch`, `_prepare_data`) — the summary-statistic
     contract Phase 2 needs. It references `corspearman`/`corkendall` (StatsBase) and
     `mean`/`median`/`quantile` (Statistics), so Phase 2 must add **StatsBase + Statistics** to
     the spike env (lightweight — does NOT pull GLMakie/Turing) before `include`-ing it.
   - `src/LoadImages.jl` (`MultiChannelImage`) and `src/bayes.jl` (Turing `colocalization()`
     model, `compute_BayesFactor`) carry heavier deps (Images, Turing, GLMakie via the module).
     Phase 4's ADVI baseline will decide between a dedicated parent-aware environment and a
     targeted `include` of just `bayes.jl` with a minimal Turing-only dep set — out of Phase-1
     scope.
   - The `src/` files retain their AGPL module-level imports; `include` is **read-only** — no
     source file is edited. This keeps the publication baseline frozen (Task 2 proof below).

**Net for ENV-01:** the parent is reachable read-only via the documented `include()` fallback;
the smoke stays green; the root is byte-identical to the baseline (see §5). Either coupling
path was a SUCCESS per D-01 — the fallback is pre-authorized, not a failure.

## 3. Phase-2 seed — Turing prior ranges (OPTIONAL, Phase-2 scope)

Pre-seeded here as a convenience for Phase 2, which will copy these into the simulator prior π(θ)
so it stays consistent with the existing Turing `@model` ranges (μ/ν/σ/τ) for ADVI comparability.
Source: `src/bayes.jl` (~lines 274–291). **Phase-2 scope — not used by the Phase-1 smoke.**

```julia
mu    ~ Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)
nu    ~ Exponential()
sigma ~ Truncated(Cauchy(0.1, 0.3), 1e-4, 1.0)
tau   ~ Truncated(Cauchy(0.1, 0.3), 1e-4, 1.0)
```
