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

## 5. Decoupling Proof (ENV-01 / DEMO-02 pattern) — Plan 01-04 Task 2

**Result: CLEAN.** The frozen publication state under the root manifests and `src/` is
**byte-identical** to the agreed Plan-01 baseline. `Pkg.develop` targeted only the spike env
and never wrote to the root; the develop attempt + revert left the protected paths **untouched**.

- **Baseline ref:** `f581d95` (`f581d95dbc3596c5d9a064d6325ed70ca43a79c0`,
  "chore(baseline): freeze root manifests as v2.0 spike baseline"), recorded in 01-BASELINE.md.
- **Protected paths:** `Project.toml`, `Manifest.toml`, `src/`.
- **Proof command (verbatim from 01-BASELINE.md, run from repo root):**
  ```sh
  git diff --quiet f581d95 -- Project.toml Manifest.toml src/ && echo "root untouched" \
    || { echo "ROOT DRIFT vs baseline f581d95"; exit 1; }
  ```
  → exited **0** (clean): `root untouched`. Additionally `git diff --quiet HEAD -- src/` → `0`
  (zero `src/` modifications).

This is the decoupling-proof pattern DEMO-02 repeats at Phase 6. Any future ENV-01 / DEMO-02
check MUST assert byte-identity against `f581d95` using the command above.

## 3. SIM-02 Induced-μ Calibration — the prior-consistency proof (Plan 02-03)

The simulator prior π(θ) is made **provably consistent** with the Turing `@model` (D-01/D-02).
This section documents the fitted map ĝ, the calibration evidence, the faithful real anchor
(D-16), and the σ/τ/ν consistency check (D-03). Source target: `src/bayes.jl` (~lines 274–291).

### Turing target ranges (the consistency target)

```julia
mu    ~ Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)   # ← the SIM-02 calibration target (induced μ)
nu    ~ Exponential()                            # ← tail heaviness, induced & checked (D-03)
sigma ~ Truncated(Cauchy(0.1, 0.3), 1e-4, 1.0)   # ← per-patch-corr spread, induced & checked
tau   ~ Truncated(Cauchy(0.1, 0.3), 1e-4, 1.0)
```

### The fitted map ĝ : μ ↦ ρ_true

- **Method:** sweep `ρ_true` → measure E[μ | ρ_true] through the UNCHANGED `patch()`/`correlation()`
  summary → **isotonic regression (PAVA)** to enforce monotonicity → invert to a **clamped
  piecewise-linear** ĝ (monotone by construction; "Don't Hand-Roll" §02-RESEARCH). Not a bespoke
  optimizer; not the standalone NormalizingFlows.jl.
- **Where the knots live:** frozen in `spike/simulator/ghat.jl` (`GHAT_MU_KNOTS`, `GHAT_RHO_KNOTS`,
  + `GHAT_MU_MIN/MAX`), auto-generated by `spike/simulator/calibration.jl`. `prior.jl` `include`s it.
- **Reproduce:** `julia --project=spike spike/simulator/calibration.jl` (offline, ≈9 min, one-time).

### Sweep + span-past-±0.9 rationale (D-15)

- **Generator:** the spike-validated **shared-latent correlated smooth Gaussian-field** generator
  (`forward.jl`, D-15) — the calibration sweeps `θ.ρ_true` into induced μ through exactly this generator.
- **Span:** `ρ_true ∈ [−0.99, 0.99]`, 25 grid points, `n_per = 100` pairs/point, calibrated at **512²**.
  The knob span **extends past ±0.9** because the generator induces only **μ ≈ ±0.76 at knob ±0.9
  post-degrade** (D-15) — the extension to ±0.99 lets ĝ reach the Turing μ-prior tails.
- **Realized induced-μ range:** **[−0.680, +0.847]** (asymmetric: the softplus + background floor
  make the positive tail reach further than the negative). `corspearman(ρ_grid, E[μ]) = 1.0` (monotone).

### Induced-vs-target μ match (the SIM-02 pass metric)

- **Pre-declared tolerance (T-02-CAL, fixed BEFORE the run, NOT tuned to pass):**
  `SIM02_W1_TOL = 0.10` (Wasserstein-1, primary) / `SIM02_KS_TOL = 0.20` (KS, secondary).
- **Result (n = 300 at 512², draw μ*~MU_PRIOR → ρ_true=ĝ(μ*) → simulate → measure induced μ):**
  **Wasserstein-1 = 0.052 < 0.10 ✓**, KS = 0.101 < 0.20 ✓.
- **Scope (D-16):** evaluated over the **physically-realized μ range [−0.680, +0.847]** — i.e. against
  `Truncated(Cauchy(0,0.3), GHAT_MU_MIN, GHAT_MU_MAX)`. The μ-prior tails beyond this range are
  **prior-only** (see real anchor); ĝ clamps ρ_true to its swept endpoints there.
- **Re-asserted in the gate:** `spike/test/test_simulator.jl` "SIM-02 prior consistency" (W1 < tol at
  512², N=120) — green under `julia --project=spike spike/test/runtests.jl`.

### Size-invariance (T-02-SI, D-08/D-09)

Re-measured E[μ | ρ] at `ρ ∈ {−0.7, 0, 0.7}` across `{256², 512², 1024²}`. **Max |Δ E[μ]| vs 512²
= 0.032** — the 512² ĝ transfers across sizes (Pearson is sample-size-robust above the 15-px floor),
so it applies to the 1376×1028 anchor and 2048² regimes without re-fitting.

### Real anchor (D-16) — faithful extraction + negative-tail reachability

- **Extraction:** the package's own loader path `load_tiff` (`Float64.(Images.Gray.(Images.load(p)))`,
  `src/LoadImages.jl:214-216`), reached read-only through the `contract.jl` `include()` coupling — **NOT
  a hand-rolled RGB→luminance reduction**. Applied to `test/test_images/{positive,negative}/*_c{1,2}.tif`.
- **Measured real induced μ:** **positive = 0.329**, **negative = 0.248** (both run through the UNCHANGED
  summary). Both are positive, in a narrow moderate band.
- **Negative-tail verdict:** the **negative μ-prior tail is NOT physically reachable** in this real
  fluorescence data (real anti-correlation never observed; even the "negative"/non-colocalized control
  sits at +0.25). It is therefore documented **PRIOR-ONLY**: the simulator *can* generate negative
  induced μ (sign-flip, D-15), but the D-01 consistency claim is **scoped to the physically-realized μ
  range** where real data lives. This honours the D-16 caveat that consistency need only hold over the
  realized range, with the negative extreme flagged as prior-only.

### Chosen nuisance ranges (Claude's discretion, D-03)

`calibration.jl` (source of truth) and `prior.jl` use IDENTICAL ranges; the calibration measured
E[μ | ρ] averaging over exactly these:

| Nuisance | Prior | Rationale |
|----------|-------|-----------|
| `spillover` | `Uniform(0.0, 0.2)` | directional bleed-through, modest |
| `autofluorescence` | `Uniform(0.0, 0.1)` | additive background offset |
| `label_efficiency` | `Uniform(0.6, 1.0)` | Bernoulli keep-probability |
| `shift_dx`, `shift_dy` | `Uniform(-1.0, 1.0)` | sub-pixel registration error |
| `noise` | `Uniform(0.0, 1.0)` | Poisson + Gaussian noise scale |

### σ / τ / ν consistency (D-03 — a CHECK, not a fit; the SIM-02 gate is on μ)

- **σ / τ (per-patch-correlation spread):** induced **pooled SD ≈ 0.317** of the per-patch ρ values
  under the prior — **inside** the Turing scale-prior support `Truncated(Cauchy(0.1,0.3),1e-4,1)`
  (i.e. `1e-4 ≤ 0.317 ≤ 1`). Consistent ✓.
- **ν (tail heaviness, `Exponential()`):** induced per-patch-correlation **excess kurtosis ≈ −0.05**
  — i.e. the well-specified induced distribution is **approximately Gaussian** (very slightly
  light-tailed), not heavy-tailed. This does **not** contradict a finite-ν Exponential prior (which
  *permits* near-Gaussian as ν grows large); heavier tails would only arise from outlier patches
  (few-effective-pixel / spillover-driven), which are modest in the clean regime. ν is **induced and
  checked, honestly near-Gaussian here**, not set.

## 5. DATA-01 Encoding Layout (D-01 / D-02) — Plan 03-02

The generation core (`spike/data/encode.jl`) emits two fixed-dimension layouts from the FROZEN 8×8
`patch_summary(mci)` matrix. The cache stores both RAW (standardization is the loader's job, D-07).

### The 128-dim D-01 vector (`encode_d01`, LOCKED)

| Rows | Content | Encoding |
|------|---------|----------|
| 1:64 | per-patch Pearson correlations, **column-major `vec()` order** of the 8×8 grid | `vec(coalesce.(M, 0.0))` — `missing` imputed to `0.0` |
| 65:128 | the parallel binary present/absent **mask**, same column-major ordering | `vec(Float64.(.!ismissing.(M)))` — `1.0` present, `0.0` missing |

- **Phase-4 input dim is unambiguously 128.** The mask channel (65:128) makes a missing patch
  self-describing: `mask=0 ⟺ value row==0`. The loader (Wave 4) z-scores ONLY rows 1:64 (and the
  continuous D-02 moment rows) and passes the binary mask through UNCHANGED (z-scoring a 0/1 mask
  re-couples folds through the mask mean).
- **Degeneracy is KEPT, never dropped (D-13).** A fully-missing 8×8 (every patch below the ≥15-px
  floor) → vals all 0, mask all 0. This mirrors `induced_mu`'s NaN-not-throw choice (contract.jl:104);
  the sample is cached (the mask encodes the degeneracy). Dropping it would distort the π(θ)-faithful
  training distribution. It is rare here because `BG_FLOOR` (forward.jl:66) keeps every pixel > 0.

### The D-02 augmented superset (`encode_aug`, `AUG_DIM = 142`)

`encode_aug(mci, M) = vcat(encode_d01(M), moments)`; rows 1:128 are exactly `encode_d01(M)`, rows
129:142 are `N_AUG_MOMENTS = 14` scalar moments (named constant = the layout contract). In append
order: **Manders M1, M2** (thresholds = `mci.otsu_threshold`), **whole-image Pearson**, patch-grid
**median, IQR, mean (=induced_μ), std, skewness, excess-kurtosis**, **fraction-missing**, and
per-channel intensity **median + IQR** (ch1, ch2). Each `skipmissing` reduction is empty-guarded
(NaN-safe) the way `induced_mu` is. The augmented set is cached so Phase-4 ABL-01 can slice extra
features with ZERO re-simulation; the exact moment list is Claude's discretion (D-02), only the
commitment to cache both variants is locked.

### D-13 negative-correlation-tail caveat (documented, NOT engineered away)

θ is sampled i.i.d. from π(θ) with **no stratification**. Per the SIM-02 real anchor (§3 above), the
**negative induced-μ tail is prior-only** — real fluorescence data never showed anti-correlation
(positive control μ≈0.33, negative control μ≈0.25). The simulator *can* generate negative induced μ
(sign-flip, D-15), so the training pool will be **sparse in the strongly-negative tail**. This is
documented as a caveat, NOT fixed by stratifying θ (stratification would distort π(θ) and break SBC,
D-13). Phase-5 SBC/coverage is therefore evaluated honestly over the realized range.
