# Phase 5: Validation Bundle (SBC + Amortized BF + OOD) - Pattern Map

**Mapped:** 2026-07-02
**Files analyzed:** 8 planned (7 new under `spike/validation/`, 1 new gate under `spike/test/`)
**Analogs found:** 8 / 8 (7 exact/role-match in-repo, 1 partial + sibling-repo template)

All new code lives under `spike/` (hard decoupling constraint, CLAUDE.md). Every planned
file inherits the same in-repo house style: GNU-AGPL license header block, a top-of-file
`# spike/<path> --- <phase-tag> (<req IDs / decisions>)` doc comment, guarded idempotent
`include(joinpath(@__DIR__, ...))` chains, `use_gpu = false` on every NeuralEstimators call,
and reaching `src/` only via read-only `include()`.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `spike/validation/harness.jl` | service (shared draw→sim→infer) | transform / request-response | `spike/npe/infer.jl` | exact |
| `spike/validation/sbc.jl` | service (calibration diagnostics) | batch / transform | `spike/npe/ablation.jl` + `spike/npe/infer.jl` | role-match |
| `spike/validation/train_ratio.jl` | service (training entry) | batch | `spike/npe/train_npe.jl` | exact |
| `spike/validation/bf.jl` | service (amortized BF + reproduction) | transform / request-response | `spike/npe/infer.jl` + `src/bayes.jl` (read-only target) | role-match |
| `spike/validation/ood.jl` | service (OOD scoring + misspec grid) | batch / transform | `spike/npe/ablation.jl` + `spike/data/loader.jl` + `spike/simulator/forward.jl` | role-match |
| `spike/validation/figures.jl` | utility (CairoMakie figures) | transform | `spike/02_simulator_demo.jl` (Phase-2 figure pattern) | partial |
| `spike/validation/consts.jl` *(optional)* | config (pre-registered constants) | — | `spike/test/test_npe.jl:70-96` const block | role-match |
| `spike/test/test_validation.jl` | test (SC1..SC5 gates) | — | `spike/test/test_npe.jl` + `spike/test/runtests.jl` | exact |

## Pattern Assignments

### `spike/validation/harness.jl` (service, transform)

**Analog:** `spike/npe/infer.jl` (composes the entire amortized read surface; Phase 5 is
*composition* of this surface, not new inference machinery — RESEARCH "Key insight").

**License + file-header + decoupling doc-comment pattern** (`infer.jl:1-56`): copy the
AGPL header block verbatim, then the `# spike/validation/harness.jl --- Phase-5 shared
harness (...)` banner naming the reqs and the CENTRAL frozen-stats / CPU-only / decoupling
invariants, then the guarded includes.

**Guarded, order-dependent includes** (`infer.jl:55-56`):
```julia
# ORDER MATTERS: trainer first (load_npe + loader + architecture), then simulator halves.
isdefined(@__MODULE__, :load_npe)     || include(joinpath(@__DIR__, "..", "npe", "train_npe.jl"))
isdefined(@__MODULE__, :posterior_for)|| include(joinpath(@__DIR__, "..", "npe", "infer.jl"))
isdefined(@__MODULE__, :sample_prior) || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair)|| include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)    || include(joinpath(@__DIR__, "..", "contract.jl"))
```

**Frozen-stats load + standardize + one-pass infer chain** — the exact composition the
harness wraps (RESEARCH Pattern 1; `infer.jl:76-99` `standardize_summary` + `posterior_for`):
```julia
m  = load_npe(joinpath(@__DIR__, "..", "npe", "trained_npe.jld2"))   # estimator, zt, θzt, d_in
Z  = standardize_summary(encode_d01(patch_summary(build_mci(sim))), m.zt, :min)  # FROZEN train zt
d  = posterior_for(m.estimator, Z; N = L)                            # 7×N, one pass, use_gpu=false
θ  = StatsBase.reconstruct(m.θzt, d)                                 # physical θ (7×N)
```
Reuse `rho_draws`/`delta_rho` (`infer.jl:120-138`) for the paired-Δρ path (D-01) directly —
do NOT re-derive; `delta_rho` already IS the MC difference of two single-stack passes (D-03).

**Draw half** — `sample_prior(rng)` is the SOLE prior (`prior.jl:71-82`); the harness only
threads a Random123 disjoint counter stream around it (D-02). `simulate_pair(rng, θ; imsize)`
is the forward half (`forward.jl:110`). Both are called UNCHANGED; the misspec grid perturbs
*around* them (see ood.jl).

---

### `spike/validation/train_ratio.jl` (service, batch)

**Analog:** `spike/npe/train_npe.jl` — CONTEXT.md/RESEARCH: "the l-POP `RatioEstimator`
likely lands as a sibling here." Mirror its every discipline; only the estimator type
(`RatioEstimator`, model-index parameter) and the paired-summary data assembly differ.

**Fixed-data `train` form + CPU-only guard** (`train_npe.jl:102-128`) — use the same
fixed-data form so the reserved-holdout + frozen-stats discipline holds; note the AdamW
Float64 args gotcha:
```julia
use_gpu && throw(ArgumentError("train_ratio: use_gpu=true is out of scope (D-10 CPU-only gate)"))
est = train(est, θtr_std, θva_std, Ztr, Zva;         # fixed-data form; θ here = model index m∈{0,1}
            epochs = epochs, batchsize = batchsize, use_gpu = false,
            optimiser = Flux.Optimisers.AdamW(5e-4, (0.9, 0.999), 1e-4),
            stopping_epochs = 10, verbose = verbose)
```
NOTE (RESEARCH Anti-Pattern): a custom loss passed to `train(::RatioEstimator, …)` is
silently ignored (loss hard-coded `logitbinarycrossentropy`). Do not attempt to pass l-POP.

**Frozen-stats reuse (SC5):** standardize BOTH paired summaries with the FROZEN `m.zt`
(`infer.jl:76` `standardize_summary`), never re-fit — the RatioEstimator must see summaries
in the same input space as the deployed NPE. Label by `sign(Δρ)` (BF-01 null, see bf.jl).

**Atomic JLD2 persistence** (`train_npe.jl:147-186` `save_npe`/`load_npe`): copy the
`.tmp` + reopen-integrity-check + `mv(...; force=true)` idiom for the trained ratio net and
the measured prior-odds constant (Pitfall 5):
```julia
tmp = path * ".tmp"
jldsave(tmp; schema_version = RATIO_MODEL_SCHEMA, estimator = est, zt = zt,
        log_prior_odds = measured_lpo, meta = (...))
JLD2.jldopen(tmp, "r") do f
    @assert haskey(f, "estimator") "save: integrity check failed"
end
mv(tmp, path; force = true)   # filesystem-atomic commit
```

---

### `spike/validation/bf.jl` (service, transform / request-response)

**Analog:** `spike/npe/infer.jl` (amortized single-pass read; `delta_rho`) for the amortized
side; **`src/bayes.jl:109` `compute_BayesFactor` (READ-ONLY reproduction target, never edit).**

**Amortized log-BF in one forward pass** (RESEARCH Code Examples; composes
`RatioEstimator`/`logratio`, no quadgk/KDE/shuffle — BF-02):
```julia
grid = reshape(Float32[0.0 1.0], 1, 2)          # null (m=0), coloc (m=1)
lr   = logratio(est_bf, Z_pair; grid = grid, use_gpu = false)   # 1×2: [logr(m=0) logr(m=1)]
log_BF = (lr[1,2] - lr[1,1]) - log_prior_odds    # subtract MEASURED prior-odds (Pitfall 5)
```

**Null mirrored EXACTLY from the baseline** (D-07; confirmed body, RESEARCH §BF Null):
Δρ = `μ_sample − μ_control`; H1 = {Δρ>0}, H0 = {Δρ≤0} (one-sided interval at threshold 0,
NOT a point null); BF = posterior-odds ÷ prior-odds for {Δρ>0}. Because `ghat` is monotone
through 0 (`ghat.jl:59`, SIM-02), the model index `sign(Δρ_ρ)` reproduces the baseline
`{Δμ>0}` event — but verify `ghat(0)` numerically (A2) and define the split on the induced-μ
contrast if `ghat(0) ≠ 0`.

**Reproduction gate (D-08)** — mirror the `ablation.jl`/`benchmark.jl` "assess-then-compare-
to-pre-registered-tolerance" shape: compute amortized log-BF and the KDE baseline BF over a
held-out Δρ sweep, assert `cor ≥ BF_CORR_MIN` AND `max|Δ logBF| ≤ BF_LOGBF_TOL`. The
tolerance-gate pattern is exactly Phase-4's D-09 RMSE gate (`test_npe.jl:73`).

---

### `spike/validation/sbc.jl` (service, batch / transform)

**Analog:** `spike/npe/ablation.jl` (the batch-loop + pre-registered-decision-rule shape) and
`spike/npe/infer.jl` (`rho_draws` for the Δρ paired path). **Template (copy pattern, NOT a
dep):** BayesInteractomics `calibration.jl` `_bin_calibration`/`CalibrationResult` and
`predictive_checks.jl` `_ks_test_uniform`.

**Batch loop + guarded-const structure** (`ablation.jl:60-72, 143-158`): copy the
`for f in 1:K, …` accumulation-into-preallocated-array shape and the guarded pre-registered
const declarations:
```julia
if !isdefined(@__MODULE__, :SBC_M)
    const SBC_M = 2000        # pre-registered draws
end
rank_table = Matrix{Int}(undef, SBC_M, 7)       # per-param ranks; +1 col for Δρ paired path
for i in 1:SBC_M
    θ★  = sample_prior(rng_sbc)                  # disjoint Random123 stream (D-02)
    Z   = standardize_summary(encode_d01(patch_summary(build_mci(simulate_pair(rng_sbc, θ★)))), m.zt, :min)
    d   = StatsBase.reconstruct(m.θzt, posterior_for(m.estimator, Z; N = SBC_L))  # 7×L physical
    for p in 1:7; rank_table[i, p] = count(<(θ★[p]), @view d[p, :]); end          # rank ∈ 0:L
end
```
Δρ paired path (D-01): reuse `rho_draws(m.estimator, Z_s, m.θzt; N=L)` and `Z_c`, then
`Δρ_draws = ρs .- ρc`, `rank = count(<(Δρ★), Δρ_draws)` (RESEARCH §Paired-draw Δρ SBC).

**KS / χ² uniformity via HypothesisTests** (RESEARCH Code Examples; Don't-Hand-Roll):
```julia
u  = (rank_table[:, p] .+ 0.5) ./ (SBC_L + 1)
ks = pvalue(ExactOneSampleKSTest(u, Uniform(0, 1)))
chi = pvalue(ChisqTest(cnts))                 # verify ChisqTest arg form at impl (A3)
pass = ks > SBC_KS_ALPHA && chi > SBC_KS_ALPHA # pre-registered thresholds
```

**ECE/MCE + traffic-light** — port `_bin_calibration` (BayesInteractomics
`calibration.jl:22-72`) verbatim into `spike/validation/` (~50 lines; do NOT add the sibling
repo as a dep) with the `CalibrationResult` struct (`types.jl:227-234`):
```julia
struct CalibrationResult
    bin_midpoints::Vector{Float64}; predicted_rate::Vector{Float64}
    observed_rate::Vector{Float64}; bin_counts::Vector{Int}; ece::Float64; mce::Float64
end
# ECE = Σ (bin_count/total)·|pred−obs|;  MCE = max gap  (calibration.jl:60-69)
```
Map the SBC PIT/coverage into the (predicted-prob, empirical-positive) binning; traffic-light
cutoffs on ECE pre-registered.

---

### `spike/validation/ood.jl` (service, batch / transform)

**Analog:** `spike/npe/ablation.jl` (batch scoring + pre-registered decision rule),
`spike/data/loader.jl` (`_row_partition` + fit-on-train-only), `spike/simulator/forward.jl`
+ `prior.jl` (the misspec generators are perturbed siblings of these).

**Fit nulls on TRAIN continuous rows only (SC5/D-06, Pitfall 2)** — reuse the loader's row
split (`loader.jl:126-136`, exposed as `infer.jl:65` `_summary_row_partition`); fit
Mahalanobis on rows 1:64 only (the 64 binary mask rows 65:128 make Σ singular):
```julia
cont, _ = Loader._row_partition(:min, 128)       # rows 1:64
μS = vec(mean(Ztrain[cont, :]; dims = 2));  ΣS = cov(Ztrain[cont, :]; dims = 2)
C  = cholesky(Symmetric(ΣS + 1e-6I))             # ridge for stability
maha(z) = (r = z[cont] .- μS; sum(abs2, C.L \ r))
```

**PP-mismatch channel** — infer θ̂ via `posterior_for`, re-simulate via `simulate_pair`,
compute discrepancy (Mahalanobis of S_obs vs the PP-summary cloud is the default, D-05);
template: `_ks_test_uniform` (`predictive_checks.jl:729-736`, copy pattern).

**Misspec generators (D-03)** — build as spike-local functions that perturb `simulate_pair`
output or replace stage-7 noise, guarded with the same `simulate_pair` entry validation
(ArgumentError on out-of-range, mirroring `train_fold`'s `use_gpu` guard, `train_npe.jl:109`).
**Negative control (D-04):** correlation-preserving transforms (per-channel affine a·x+b,
rotation/flip, block-permute); verify KS-invariance on the 8×8 `patch_summary` distribution
with `ExactOneSampleKSTest` and assert statistic < pre-registered ε.

**ROC/AUC** — hand-roll from sorted scores (~15 lines, no new dep — Pitfall 7 co-resolve
landmine); operating point = pre-registered ID quantile; Youden-J shown as post-hoc label
only (D-06).

---

### `spike/validation/figures.jl` (utility, transform) — partial analog

**Analog:** `spike/02_simulator_demo.jl` (Phase-2 CairoMakie figure pattern) — no exact
in-`validation/` analog. Planner should follow the RESEARCH figure set (rank histograms,
coverage curves, reliability/ECE, ROC, traffic-light) and the Phase-2 headless-save idiom.
Keep figure generation OUT of the fast test gate (Phase-4 pattern: figures are script
artifacts, not gate assertions).

---

### `spike/test/test_validation.jl` (test) + `runtests.jl` wiring

**Analog:** `spike/test/test_npe.jl` (structure) + `spike/test/runtests.jl` (wiring +
resolve-risk gate).

**Pre-registered const block committed BEFORE the reported run** (`test_npe.jl:70-96`,
SBC-03/D-02/D-06/D-08) — the anti-snooping contract; note the DISJOINT master seed:
```julia
const SBC_M           = 2000        # pre-registered draws (SBC-01/03)
const SBC_L           = 999         # posterior draws/SBC draw (L+1 divisible by SBC_BINS, Pitfall 4)
const SBC_BINS        = 50
const SBC_KS_ALPHA    = 0.05
const OOD_ID_QUANTILE = 0.95        # ~5% ID FPR operating point (D-06)
const BF_CORR_MIN     = 0.95        # D-08a
const BF_LOGBF_TOL    = 0.5         # D-08b
const VAL_MASTER_SEED = 0x5BC0FFEE  # DISJOINT from NPE_MASTER_SEED 0xC0FFEE (D-02 fresh stream)
```

**Fast-gate-loads-persisted-artifact + tiny fixture** (`test_npe.jl:80-134`, SC1): the quick
gate loads `trained_npe.jld2` and the trained ratio net (does NOT retrain), runs SBC/BF/OOD
on a SMALL fixture (tiny M/L, few grid points); the REPORTED numbers come from a separate
script at pre-registered M/L. Copy the `generate_cache(mktempdir(); …)` fixture idiom and the
guarded includes of the units under test (`test_npe.jl:54-68`).

**Wire into `runtests.jl`** (`runtests.jl:92-103`): add
`include(joinpath(@__DIR__, "test_validation.jl"))`. **Extend the resolve-risk gate**
(`runtests.jl:52-86`) to re-assert `NeuralEstimators == v"0.2.1"` (UUID
`38f6df31-6b4a-4144-b2af-7ace2da57606`) AFTER validation code loads and that no ROC/flow
package or Turing entered the env (Pitfall 7).

## Shared Patterns

### Frozen train-split standardization (applies to harness, sbc, train_ratio, bf, ood — SC5/D-06)
**Source:** `spike/npe/infer.jl:76-84` (`standardize_summary`) + `train_npe.jl:180-186`
(`load_npe` → `m.zt`, `m.θzt`).
```julia
m = load_npe(".../trained_npe.jld2")                       # frozen estimator + zt + θzt
Z = standardize_summary(rawS, m.zt, :min)                  # continuous rows z-scored, mask 65:128 passthrough
θ = StatsBase.reconstruct(m.θzt, draws)                    # un-standardize before any ρ read (Pitfall 5)
```
NEVER re-fit standardization or OOD nulls on eval data (Anti-Pattern). All OOD nulls fit on
TRAIN continuous rows only (`loader.jl:181-182` fit-on-train + `_row_partition`).

### Pre-registered consts + fresh disjoint seed (applies to sbc, bf, ood, test_validation — D-02)
**Source:** `spike/test/test_npe.jl:70-96` (locked const block, committed before the reported
run) + `train_npe.jl:64-66` (guarded `const` so redefinition to the same value under
`runtests.jl` is a silent no-op). Reported numbers come from a Random123 counter range
disjoint from training/holdout (`VAL_MASTER_SEED ≠ NPE_MASTER_SEED`).

### Atomic JLD2 persistence (applies to train_ratio, ood nulls)
**Source:** `spike/npe/train_npe.jl:147-186` (`save_npe`/`load_npe`): `jldsave` to `.tmp` →
reopen integrity-check the required key → `mv(...; force=true)`. A crash before `mv` leaves
only a discardable `.tmp`, never a half-written artifact.

### CPU-only discipline (applies to ALL files — D-10/Pitfall 1)
**Source:** every NeuralEstimators call in `infer.jl`/`train_npe.jl`/`ablation.jl` passes
`use_gpu = false` (it defaults to TRUE). `train_fold`/`fold_rmse` also hard-`throw` on
`use_gpu=true` (`train_npe.jl:109`, `ablation.jl:101`) — mirror this guard. CUDA is never
imported; the `runtests.jl` gate asserts it stays out of `Base.loaded_modules`.

### Read-only `src/` coupling (applies to bf.jl — hard constraint until Phase 7)
**Source:** `spike/contract.jl:41-47` — `using StatsBase/Statistics` BEFORE
`include("../src/...")` (correlation() builds its method Dict at call time). Reach
`compute_BayesFactor`/`correlation` ONLY via read-only `include()`; never edit `src/`.

### Guarded, order-dependent includes (applies to ALL files)
**Source:** `spike/npe/infer.jl:55-56`, `ablation.jl:62`, `train_npe.jl:59-60`:
`isdefined(@__MODULE__, :sym) || include(joinpath(@__DIR__, ...))` — keeps each file loadable
standalone AND idempotent inside `runtests.jl`.

## No Analog Found

| File | Role | Data Flow | Reason |
|------|------|-----------|--------|
| `spike/validation/figures.jl` | utility (viz) | transform | No CairoMakie figure module exists under `spike/npe/` or `spike/validation/`; only the Phase-2 demo script (`02_simulator_demo.jl`) is a loose precedent. Planner should follow the RESEARCH figure set + Phase-2 headless-save idiom; keep figures out of the fast gate. |

Note: `spike/validation/consts.jl` is optional — the RESEARCH structure allows keeping the
pre-registered constants inline in `test_validation.jl` (the `test_npe.jl` precedent). Only
factor out `consts.jl` if the same constants are needed by both the scripts and the gate.

## Conventions

Convention derivation skipped (`gsd-tools verify conventions --derive --scope spike` returned
`{ "skipped": true, "reason": "no-readable-files" }` — the deterministic module could not
enumerate the scope in this environment). The axis table below is therefore stated from the
directly-read `spike/` analogs rather than the tool's majority-vote output; treat it as
descriptive, not tool-derived.

| Axis | Dominant | Share | Entropy | Status |
|------|----------|-------|---------|--------|
| File-name casing | `snake_case.jl` (`train_npe.jl`, `test_npe.jl`, `loader.jl`) | ~100% (spike) | low | named contract |
| Identifier casing | `snake_case` funcs (`sample_prior`, `load_fold`, `posterior_for`); `SCREAMING_SNAKE` consts (`NPE_MASTER_SEED`, `SBC_M`); Unicode math (`θzt`, `ρ_true`, `Δρ`, `μS`) for math objects | ~90% | low-med | named contract |
| Export / module style | flat top-level functions per file (no module) EXCEPT `loader.jl` which uses `module Loader` + `export` + trailing `using .Loader` to control the public surface | mixed | med | contested hotspot |
| Include / dep style | guarded `isdefined(@__MODULE__, :sym) || include(joinpath(@__DIR__, ...))`; `using <Pkg>` with a trailing purpose comment; `src/` only via read-only `include()` | ~100% | low | named contract |

**Contested hotspots (author's choice).** The **module-vs-flat** split is the one genuinely
contested axis and it is *intentional*: `spike/data/loader.jl` wraps itself in `module Loader`
with an explicit `export` list precisely to make `names(Loader)` a controlled public surface
(the D-08 leak-free guarantee — no `standardize_all` symbol can exist), then re-exports via
`using .Loader`. Every other spike file uses flat top-level functions pulled in by guarded
`include`. Each half is internally consistent for its purpose (controlled-surface module vs.
composed-include utility); the split is contested only when viewed repo-wide. New Phase-5
files should match the LOCAL style of what they compose: `harness.jl`/`sbc.jl`/`bf.jl`/
`ood.jl`/`train_ratio.jl` follow the flat `spike/npe/*.jl` include-composition style (they
are siblings of `infer.jl`/`train_npe.jl`/`ablation.jl`), and should NOT introduce a new
`module` wrapper unless they need loader-style surface control. Identifier casing follows the
math/code split already in `infer.jl` (Unicode for math objects θ/ρ/Δ/μ, snake_case for
plumbing, SCREAMING_SNAKE for pre-registered constants).

## Metadata

**Analog search scope:** `spike/npe/`, `spike/data/`, `spike/simulator/`, `spike/test/`,
`spike/contract.jl`, `src/bayes.jl` (read-only target), and the sibling
`~/Documents/GitHub/BayesInteractomics/src/diagnostics/` (template only).
**Files scanned (read in full or targeted):** `spike/npe/infer.jl`, `train_npe.jl`,
`ablation.jl`; `spike/data/loader.jl`, `encode.jl`; `spike/simulator/prior.jl`, `forward.jl`
(signatures), `ghat.jl` (signatures); `spike/contract.jl`; `spike/test/test_npe.jl`,
`runtests.jl`; BayesInteractomics `calibration.jl`, `predictive_checks.jl`, `types.jl`.
**Pattern extraction date:** 2026-07-02
