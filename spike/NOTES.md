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
