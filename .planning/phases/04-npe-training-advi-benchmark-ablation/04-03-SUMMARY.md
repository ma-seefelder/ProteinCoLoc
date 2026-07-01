---
phase: 04-npe-training-advi-benchmark-ablation
plan: 03
subsystem: npe
tags: [neuralestimators, posteriorestimator, normalisingflow, npe, flux, cpu-only, zscore, jld2, delta-rho, ghat]

# Dependency graph
requires:
  - phase: 03-training-data-pipeline
    provides: "Loader.load_fold (leak-free Ztr/θtr/zt), load_holdout (reserved ≥20-stack set), _row_partition (mask-bypass split), encode_d01/patch_summary"
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 01
    provides: "test_npe.jl scaffold + six pre-registered constants (NPE_MASTER_SEED), resimulate helper, resolve-risk gate"
  - phase: 02-forward-simulator-summary-contract
    provides: "spike/simulator/ghat.jl frozen μ→ρ inverse (SIM-02)"
provides:
  - "spike/npe/architecture.jl — build_estimator(d_in, D=7; dstar,depth,width) → PosteriorEstimator (MLP summary net → NormalisingFlow); locked NPE_DSTAR/DEPTH/WIDTH so ablation/benchmark share the arch"
  - "spike/npe/train_npe.jl — train_fold(...) CPU-only fixed-data train on leak-free loader tensors + leak-free θ ZScoreTransform + save_npe/load_npe atomic JLD2"
  - "spike/npe/infer.jl — posterior_for / rho_hat / delta_rho (D-03) / interval_width / rho_from_mu (ghat μ→ρ axis) / standardize_summary"
  - "spike/npe/trained_npe.jld2 — persisted CPU-only estimator + frozen θ/zt transforms for the 04-05 benchmark and 04-06 ablation"
  - "SC1 (NPE-01) green in spike/test/test_npe.jl"
affects: [04-05, 04-06, phase-5-sbc-bf-ood]

# Tech tracking
tech-stack:
  added: []
  patterns: ["fixed-data NeuralEstimators train on the leak-free loader cache", "leak-free θ ZScoreTransform (fit-on-θtr, reconstruct draws) mirroring the loader Z path", "single-source build_estimator arch shared across trainer/benchmark/ablation", "ρ-space benchmark alignment via the frozen ghat (rho_from_mu)", "fast SC gate loads the persisted model rather than retraining"]

key-files:
  created:
    - spike/npe/architecture.jl
    - spike/npe/train_npe.jl
    - spike/npe/infer.jl
    - spike/npe/trained_npe.jld2
  modified:
    - spike/test/test_npe.jl

key-decisions:
  - "AdamW optimiser args are Float64 (5e-4,(0.9,0.999),1e-4): a Float32-parameterised AdamW state trips Optimisers.adjust! against NeuralEstimators' Float64 CosAnneal lr_schedule (blocking-issue fix)"
  - "A single-dataset summary MUST be a d_in×1 matrix, not a bare Vector (NeuralEstimators reads a Vector as many 1-element datasets); posterior_for coerces vectors"
  - "θ is standardized leak-free (fit on θtr only) and posterior draws are un-standardized via StatsBase.reconstruct before any ρ read (Pitfall 5); the θ-transform is frozen alongside the model"
  - "SC1 loads the persisted trained_npe.jld2 (produced by train_fold @ NPE_MASTER_SEED) and evaluates on the reserved holdout — a fast, deterministic, non-flaky gate ('load from cache where practical')"
  - "Architecture consts are NPE_-prefixed to avoid clobbering 00_smoke.jl's `const d`/`num_summaries` when both load under runtests.jl"

requirements-completed: [NPE-01]

# Metrics
duration: ~40min
completed: 2026-07-01
---

# Phase 4 Plan 03: NPE Training + Inference Surface Summary

**A CPU-only `PosteriorEstimator` (MLP summary net → NeuralEstimators `NormalisingFlow` over the 7-dim θ, ρ_true row 1) trains fixed-data on the Phase-3 leak-free loader tensors with a frozen leak-free θ-standardization, is persisted to `trained_npe.jld2`, and drives an inference surface (single-pass posteriors, ρ̂, the Δρ MC-difference (D-03), interval widths, and the frozen `ghat` μ→ρ benchmark axis) — turning SC1 (NPE-01) green with ρ̂ recovering holdout ρ_true at corr ≈ 0.98.**

## Performance

- **Duration:** ~40 min
- **Tasks:** 3
- **Files modified:** 5 (4 created, 1 modified)

## Accomplishments
- **`architecture.jl`** — `build_estimator(d_in, D=7; dstar=32, depth=2, width=128)` returns a `PosteriorEstimator(Chain(Dense(d_in,128,gelu), Dense(128,128,gelu), Dense(128,32)), NormalisingFlow(7; num_summaries=32))`. `q` is constructed as an INSTANCE passed POSITIONALLY (the v0.2.1 gotcha), and the topology is exposed as locked `NPE_*` consts so the ablation (04-06) is a fair arch-controlled comparison (D-05).
- **`train_npe.jl`** — `train_fold(dir, fold; master_seed, variant=:min, use_gpu=false, epochs=200, batchsize=64, savepath=nothing)` loads the leak-free fold tensors, fits a leak-free θ `ZScoreTransform` on θtr only (Pitfall 5), standardizes θtr/θva, builds the estimator, and calls the FIXED-DATA `train(est, θtr_std, θva_std, Ztr, Zva; use_gpu=false, optimiser=AdamW(5e-4,(0.9,0.999),1e-4), stopping_epochs=10)`. `save_npe`/`load_npe` persist/reload {estimator, θ-transform, zt, meta} via the atomic `.tmp`→reopen-assert→`mv(force=true)` cache.jl idiom. `use_gpu=false` on every call (D-10).
- **`infer.jl`** — `posterior_for` (7×N single pass; vector→d_in×1 coercion), `rho_hat`/`rho_draws` (un-standardize with `reconstruct` then read row 1), `delta_rho` (mean of two independent single-stack passes, D-03), `interval_width` (via NeuralEstimators `interval`), `rho_from_mu = mean(ghat.(μ_draws))` (the frozen ρ-space benchmark axis, Pattern 2/A1), and `standardize_summary` (apply the frozen loader `zt` to a raw holdout summary with mask-row bypass).
- **`trained_npe.jld2`** — a genuine CPU-only estimator trained on a reproducible N=256 fold-1 cache at `NPE_MASTER_SEED`, with both frozen transforms, reloading cleanly for downstream 04-05/04-06.
- **SC1 (NPE-01) green** — 7×N single-pass posteriors for all 20 reserved holdout stacks, every ρ̂ ∈ [-1,1], a finite Δρ equal to the seeded MC-diff of two passes, CPU-only. Empirical signal: ρ̂ vs known holdout ρ_true correlates at **0.983**.

## Task Commits

1. **Task 1: NPE architecture + CPU-only fixed-data training** — `5460529` (feat)
2. **Task 2: Inference surface (rho_hat/delta_rho/ghat map)** — `8139ab3` (feat)
3. **Task 3: Turn SC1 (NPE-01) green** — `4022e57` (test)

## Files Created/Modified
- `spike/npe/architecture.jl` (created) — `build_estimator` + `NPE_D/DSTAR/DEPTH/WIDTH` consts.
- `spike/npe/train_npe.jl` (created) — `train_fold`, `fit_theta_transform`, `save_npe`, `load_npe`.
- `spike/npe/infer.jl` (created) — the five inference functions + `rho_draws` + `standardize_summary`.
- `spike/npe/trained_npe.jld2` (created) — persisted trained estimator + frozen transforms (~2.2 MB).
- `spike/test/test_npe.jl` (modified) — SC1 `@test_skip` replaced with the live NPE-01 gate; guarded `include` of infer.jl added; SC2..SC5 + Wave-0 gate untouched.

## Decisions Made
- **Float64 AdamW args:** a Float32-parameterised AdamW state raised an `Optimisers.adjust!` type error against NeuralEstimators' Float64 `CosAnneal` lr_schedule; matching the RESEARCH example's Float64 literals fixes it.
- **Single-dataset = d_in×1 matrix:** a bare Vector is interpreted by NeuralEstimators as many 1-element datasets (a DataLoader over scalars); `posterior_for` coerces a Vector to a column matrix so one stack yields one 7×N posterior.
- **SC1 loads the persisted model:** rather than retrain (cache-gen is imsize-draw-variable and training adds ~45 s), SC1 loads `trained_npe.jld2` and reads the reserved holdout from `NPE_REPRO_DIR`. The holdout stream depends only on `(master_seed, n_holdout)`, so it is byte-identical to the model's held-out set and was excluded from training by construction (D-10) — fast, deterministic, non-flaky, and it exercises the persisted artifact downstream consumes.
- **`NPE_`-prefixed arch consts:** avoids colliding with `00_smoke.jl`'s `const d`/`num_summaries` when both files load in `runtests.jl`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] AdamW/lr-schedule dtype mismatch**
- **Found during:** Task 1 (first training run)
- **Issue:** `Flux.Optimisers.AdamW(5f-4,(0.9f0,0.999f0),1f-4)` (Float32) tripped `Optimisers.adjust!` with an `eta::Float64` from NeuralEstimators' `CosAnneal` lr_schedule → training aborted.
- **Fix:** Use Float64 AdamW args `(5e-4,(0.9,0.999),1e-4)` (matches RESEARCH Pattern 1 example).
- **Files modified:** spike/npe/train_npe.jl
- **Commit:** 5460529

**2. [Rule 3 - Blocking] Bare-vector single-dataset input**
- **Found during:** Task 2 (end-to-end functional check)
- **Issue:** `sampleposterior(est, Zvec)` on a `Vector{Float32}` was read as a DataLoader of 1-element datasets, erroring in the summary network's `Dense` (dimension mismatch).
- **Fix:** `posterior_for` coerces an `AbstractVector` to a `d_in×1` matrix (one dataset) before sampling; documented the gotcha.
- **Files modified:** spike/npe/infer.jl
- **Commit:** 8139ab3

Both are Rule-3 blocking-issue fixes made inline (no package installs, no architectural change).

## Verification Evidence
- Task 1 verify: `include(train_npe.jl)` loads, CUDA absent from `Base.loaded_modules`, `build_estimator(128)` returns a `PosteriorEstimator`.
- Artifact: cache gen 7 s + train 45 s (N=256, CPU-only); `trained_npe.jld2` written atomically and reloads to a usable estimator (d_in=128, variant=:min).
- Task 2 verify: all five functions defined; `rho_from_mu` in range; full functional pass on the 20-stack holdout (ρ̂∈[-1,1], Δρ finite, widths≥0, corr(ρ̂,true)=0.983), CPU-only.
- Task 3: `julia --project=spike spike/test/test_npe.jl` → SC1 68 pass, SC2..SC5 skipped, Wave-0 gate 11 pass (13 s).
- Full suite: `julia --project=spike spike/test/runtests.jl` → EXIT 0; CPU smoke 11/11 (D-04 CPU-only 10/10), SIM-02/04 + D-15 green, Phase-3 69/69, Phase-4 79 pass + 4 skipped. No CUDA loaded.

## Known Stubs
SC2 (NPE-02), SC3 (NPE-03), SC4 (ABL-01), SC5 (ABL-02) remain intentional `@test_skip` placeholders filled by later Phase-4 waves (04-05 benchmark, 04-06 ablation) — the enumerated-gate pattern from 04-01, not defects. No data-flow or hardcoded-value stub exists in the delivered NPE code.

## Next Phase Readiness
- `trained_npe.jld2` + `load_npe` give 04-05/04-06 a persisted CPU-only estimator with both frozen transforms.
- `rho_from_mu` (frozen `ghat`) is the ready cross-method ρ-space axis the 04-05 benchmark maps ADVI μ-posteriors onto.
- `standardize_summary` reproduces the loader's leak-free standardization for any raw holdout stack.
- `train_fold(...; variant=:min|:aug)` is arch-identical across variants, ready for the 04-06 k=5 CV ablation.

## Self-Check: PASSED
