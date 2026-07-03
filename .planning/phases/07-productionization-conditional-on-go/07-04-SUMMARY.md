---
phase: 07-productionization-conditional-on-go
plan: 04
subsystem: per-grid CPU-reproducible ship-gate machinery (SBC/BF/OOD harness + fresh pre-registration + GPU-train smoke)
status: COMPLETE
tags: [ship-gate, d-05, sbc, bayes-factor, ood, fresh-pre-registration, prod-seed-disjoint, cpu-reproducible, use-gpu-false, gpu-train-smoke, graceful-cpu-fallback, injectable-simulator]
requirements: [PROD-02]
dependency_graph:
  requires:
    - 07-00 (co-resolution GREEN: NeuralEstimators 0.2.1 / Flux 0.16.10; HypothesisTests/Distributions/JLD2 direct deps)
    - 07-01 (grid-parametric summary.jl: patch_summary(mci,G)/encode_d01/_summary_row_partition; registry.jl has_cuda_device)
    - 07-02 (read surfaces: standardize_summary/posterior_for/rho_draws (infer.jl); amortized_log_bf/pair_encode/kde_log_bf_unclamped (bf.jl); maha_score/id_threshold/roc_auc/ood_verdict (ood.jl))
    - 07-03 (train_npe/fit_summary_transform; save_estimator/load_estimator/load_ratio/load_ood_nulls CPU-resident; default_artifacts_root)
  provides:
    - "test/gate/harness.jl: draw_simulate_infer(m,rng;G,...) + draw_simulate_infer_paired (grid-parametrized, use_gpu=false, per-grid net via load_estimator) + default_simulator injection seam"
    - "test/gate/sbc.jl: sbc_ranks (M×8, 7 θ + dedicated Δρ column), sbc_uniformity (KS+χ² via HypothesisTests), sbc_coverage, sbc_calibration (ECE/MCE CalibrationResult), sbc_traffic_light, aggregated sbc_gate verdict"
    - "test/gate/gate_consts_template.jl: fresh-per-grid pre-registration (SBC/BF/OOD consts) + disjoint PROD_SEED[G] (Philox salted off VAL_MASTER_SEED/NPE_MASTER_SEED) + prod_rng(G)"
    - "test/gate/run_gate.jl: per-grid CLI --grid G [--sbc --bf --ood], bf_gate (amortized vs non-clamped KDE baseline), ood_gate (ID operating point + injected-control AUC), atomic write_gate_report; invocable pre-training (:not_trained)"
    - "test/gpu_smoke.jl: train_npe(...;use_gpu=true) GPU-or-graceful-CPU-fallback + CPU-resident persistence assertion (has_cuda_device()-guarded CUDA branch)"
  affects:
    - "the per-grid D-05 ship-gate plans (07-05+) invoke run_gate.jl per grid, copy gate_consts_template.jl → gate_consts_<G>.jl with fresh values, and inject the promoted simulator + OOD control families"
    - "the placeholder EstimatorBundle.calibration is filled by a passed sbc_gate/bf_gate/ood_gate run per grid"
tech_stack:
  added: []
  patterns: [fresh disjoint PROD_SEED[G] via Philox XOR-salt (anti-snooping T-7-08), use_gpu=false CPU-resident gate inference (T-7-07 reproducibility), injectable simulator seam (default_simulator, testable before forward-model promotion), non-clamped KDE BF baseline (Memo §5/T-7-06), atomic .tmp→reopen-integrity→mv gate-report writer, standalone-loadable guarded-include test-file convention, has_cuda_device()-guarded GPU-train graceful CPU fallback (D-06)]
key_files:
  created:
    - test/gate/harness.jl
    - test/gate/sbc.jl
    - test/gate/gate_consts_template.jl
    - test/gate/run_gate.jl
    - test/gpu_smoke.jl
    - .planning/phases/07-productionization-conditional-on-go/07-04-SUMMARY.md
  modified:
    - test/runtests.jl
decisions:
  - "Gate inference is CPU-only everywhere: every harness/sbc/bf/ood call passes use_gpu=false against a CPU-resident frozen net loaded via load_estimator (T-7-07) — the pre-registered numbers are deterministic regardless of training hardware"
  - "PROD_SEED[G] is a FULL-WIDTH Philox draw keyed (PROD_MASTER ⊻ PROD_SALT, G), asserted ∉ {VAL_MASTER_SEED=0x5BC0FFEE, NPE_MASTER_SEED=0xC0FFEE} for G∈{4,8,16,32} (Pitfall 3/T-7-08); PROD_SALT distinct from the spike VAL_SALT"
  - "The forward simulator (sample_prior/simulate_pair/build_mci) is NOT yet promoted into src (only referenced in datagen/ood bodies), so the harness takes an INJECTABLE `sim` seam (default_simulator resolves the promoted chain, errors clearly until then) — mirroring 07-03's injectable datagen; the fixture SBC smoke injects a lightweight fake simulator"
  - "gate_consts_template.jl IS a valid, includable pre-registration file (the canonical values); each per-grid plan copies it to gate_consts_<G>.jl and re-expresses M/L/imsize + PROD_SEED FRESH before that grid's gate runs"
  - "run_gate.jl is invocable before any grid is trained: a missing NPE artifact returns status=:not_trained (no crash) — the full reported run happens in the per-grid plans. The BF gate is self-contained (amortized_log_bf vs non-clamped kde_log_bf_unclamped, no Turing); the OOD gate computes the ID operating point now and the control-separability AUC only when the per-grid plan injects pos_sim/neg_sim"
  - "GPU-train smoke asserts the graceful use_gpu=true→CPU fallback path (this machine has has_cuda_device()==false) and CPU-resident persistence; the CUDA-present assertion is guarded by has_cuda_device() so CPU-only runs pass"
  - "runtests.jl wiring for BOTH smokes landed in the Task-3 commit (after every referenced file exists) so each intermediate per-task commit keeps the suite green"
metrics:
  tasks_completed: 3
  tasks_total: 3
  files_created: 6
  files_modified: 1
  completed_date: 2026-07-03
---

# Phase 7 Plan 04: Per-grid CPU-reproducible Ship-Gate Machinery (D-05) Summary

**One-liner:** Built the per-grid ship-gate the shipped grids will each pass — a
grid-parametrized CPU (`use_gpu=false`) SBC/BF/OOD harness against a CPU-resident frozen net, a
fresh-per-grid `gate_consts` pre-registration template carrying a disjoint `PROD_SEED[G]` (salted
off the spike's `VAL_MASTER_SEED`/`NPE_MASTER_SEED`, unit-asserted distinct), a `run_gate.jl`
per-grid CLI, and a GPU-train smoke proving the D-06 `use_gpu=true` → graceful CPU fallback with
CPU-resident persistence — all green under `Pkg.test()`, no grid trained yet.

## What was built

### Task 1 — grid-parametrized CPU gate harness + SBC (`test/gate/{harness,sbc}.jl`)
Promoted `spike/validation/{harness,sbc}.jl` with the three D-05 deltas:
- **`draw_simulate_infer(m, rng; G, imsize, N, sim)`** + **`draw_simulate_infer_paired`** — the
  θ*~π → simulate → CPU-infer chain, now **grid-parametrized** (`patch_summary(mci, G)`), every
  inference call **`use_gpu = false`** against a **per-grid** net loaded via
  `ProteinCoLoc.load_estimator` (not the fixed spike `trained_npe.jld2`). Frozen-stats discipline:
  `m.zt`/`m.θzt` applied, never re-fit; ρ read after `StatsBase.reconstruct` (Pitfall 5). The paired
  path also returns `Zs`/`Zc` for the BF gate.
- **`sbc_ranks`** — the M×8 rank table (7 θ + a dedicated paired-draw Δρ column), rank =
  `count(<(θ*), draws)`. **`sbc_uniformity`** — KS (`ExactOneSampleKSTest`) + χ² (`ChisqTest`)
  p-values via HypothesisTests (never hand-rolled). **`sbc_coverage`**, **`sbc_calibration`** →
  the ported `CalibrationResult` ECE/MCE, **`sbc_traffic_light`**, and an aggregated **`sbc_gate`**
  verdict over all 8 columns scored against `SBC_KS_ALPHA`/`SBC_ECE_GREEN`.

### Task 2 — fresh pre-registration template + `run_gate.jl` CLI (`test/gate/{gate_consts_template,run_gate}.jl`)
- **`gate_consts_template.jl`** — the FRESH-per-grid pre-registration: `SBC_M/L/BINS`,
  `SBC_KS_ALPHA`, `SBC_ECE_GREEN/YELLOW`, `SBC_IMSIZE`, `BF_CORR_MIN`/`BF_LOGBF_TOL`/`BF_SWEEP_*`,
  `OOD_ID_QUANTILE`/`OOD_AUC_MIN`/`OOD_GRID_LEVELS`/`OOD_PP_REPS`, plus **`PROD_SEED[G]`** derived
  by a Philox stream salted (`PROD_SALT`) **disjoint** from the forbidden `VAL_MASTER_SEED`
  (`0x5BC0FFEE`) and `NPE_MASTER_SEED` (`0xC0FFEE`), and **`prod_rng(G)`**. Header documents:
  values are re-expressed fresh per grid, the BF gate uses the non-clamped baseline, the OOD gate
  runs `with_pp = true`.
- **`run_gate.jl`** — the CLI `--grid G [--sbc] [--bf] [--ood] [--artifacts-root PATH]`, loading
  the grid consts (else the template), the CPU-resident net, and the ratio/OOD artifacts for the
  selected gates. **`bf_gate`** scores amortized `amortized_log_bf` vs the **non-clamped**
  `kde_log_bf_unclamped` baseline (self-contained, no Turing). **`ood_gate`** computes the density
  ID-score distribution + pre-registered `id_threshold`, and the control-separability ROC AUC when
  the per-grid plan injects `pos_sim`. Reports written through the atomic `.tmp`→reopen-integrity→
  `mv(force=true)` **`write_gate_report`**. Invocable before any grid is trained → `:not_trained`.

### Task 3 — GPU-train smoke, graceful CPU fallback (`test/gpu_smoke.jl`)
`train_npe(...; use_gpu = true)` runs on GPU when present and **degrades cleanly to CPU** when
CUDA is absent (no error — CLAUDE.md graceful-fallback constraint), asserts the frozen net persists
**CPU-resident** (`Flux.state` key, not the whole object) and reloads/infers CPU-side; the
CUDA-present branch is guarded by `has_cuda_device()` so CPU-only CI passes on the fallback path.

## Verification

- `julia --project -t 1 -e 'using Pkg; Pkg.test()'` — **green**: co-resolution gate 4/4 intact
  (NeuralEstimators 0.2.1 / Flux 0.16.10), new **`ship-gate harness + SBC (gate/, 07-04)` 29/29**,
  **`GPU-train smoke (graceful CPU fallback, D-06)` 8/8**; all pre-existing testsets still pass.
- `julia --project -e 'include("test/gpu_smoke.jl")'` — 8/8 standalone.
- The fixture SBC smoke drives the harness with an injected fake simulator (forward model not yet
  promoted), asserts `PROD_SEED[G] ∉ {VAL_MASTER_SEED, NPE_MASTER_SEED}` for G∈{4,8,16,32} + 4
  distinct seeds, and that `run_gate(4)` returns `:not_trained` on an untrained grid.

## Deviations from Plan

None — plan executed as written. Three scope-consistent judgment calls (all within the plan's
"full run happens in the per-grid plans" framing):
1. **Injectable simulator seam (`default_simulator`/`sim=`).** The forward model
   (`sample_prior`/`simulate_pair`/`build_mci`) is not yet promoted into `src` (only referenced in
   `datagen.jl`/`ood.jl` bodies), so the harness takes an injectable `sim` (mirroring 07-03's
   injectable `datagen`). This is what makes the fixture SBC smoke runnable today; the real gate
   resolves the promoted chain automatically once a per-grid plan lands it.
2. **`run_gate.jl` `:not_trained` honest exit.** No grid is trained in this plan, so a missing NPE
   artifact returns `:not_trained` rather than erroring — satisfying "`--grid 8 --sbc` invocable".
3. **runtests wiring committed in Task 3.** The single `runtests.jl` edit wires BOTH smokes; it
   landed in the Task-3 commit (after `gpu_smoke.jl` and all gate files exist) so every
   per-task commit keeps `Pkg.test()` green.

## Idempotency / reproducibility

- **CPU-reproducible:** every gate inference is `use_gpu=false` against a CPU-resident net;
  `prod_rng(G)` (Philox, keyed) is bit-identical per grid, so a re-run reproduces the numbers.
- **Fresh disjoint seeds:** `PROD_SEED[G]` never reuses the training/validation streams (asserted).
- **Atomic writes:** the gate-report writer reuses the `.tmp`→reopen-integrity→`mv(force=true)`
  wrapper; a crash leaves only a discardable `.tmp`.
- **Thread-independence:** gates run `-t 1`; data-gen (per-grid plans) runs `-t auto` and is
  byte-identical (Philox-per-index).

## spike/ decoupling

`spike/Project.toml` and `spike/Manifest.toml` **provably untouched** (`git status --short spike/`
empty). Only `spike/validation/{harness,sbc,consts}.jl` were READ to port logic.

## Self-Check: PASSED
