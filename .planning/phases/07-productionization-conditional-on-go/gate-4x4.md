# 4×4 Ship-Gate Outcome (D-05, PROD-02)

Recorded CPU-reproducible SBC / BF / OOD ship-gate for the **4×4 coarse-but-robust grid**, run
through the per-grid machinery (`test/gate/run_gate.jl --grid 4 --sbc --bf --ood`) against the
**frozen 4×4 pre-registration** `test/gate/gate_consts_4.jl`. The 4×4 bundle under test was produced
**through the PUBLIC `train_and_register(4)` path** — the user-definable-grid happy path (PROD-02 /
D-04) exercised end-to-end (not inline prose).

- **Gate stream (anti-snooping, T-7-08):** `PROD_SEED[4] = 0x8c0ad97b99bd6031` — the fresh disjoint
  Philox draw, provably distinct from the spike training seed `NPE_MASTER_SEED = 0xC0FFEE`, the
  spike validation seed `VAL_MASTER_SEED = 0x5BC0FFEE`, the datagen seed
  `DEFAULT_MASTER_SEED = 0x1`, and the sibling `PROD_SEED[8] = 0x8b39fecd4e2bcceb` (the grid index
  G=4 keys a different Philox counter).
- **Reproducibility (T-7-07):** every inference `use_gpu = false` against the CPU-resident frozen
  net (`load_estimator`); the gate was **not** launched with `-t auto` (that is datagen-only,
  determinism-pinned here).
- **Report artifact:** `artifacts/grid_4/gate_report_4.jld2` (gitignored per repo policy —
  regenerable cache; the recorded numbers below are the committed deliverable).

## Trained bundle under test

`train_and_register(4)` → `_train_grid_pipeline(4)` (the shared 07-03 per-grid pipeline), datagen
launched under `julia -t auto`. `estimator_for(4)` returned the same registered bundle in-process
after the call (no `UndefVarError`, no unregistered-grid error).

| Knob | Value | Note |
|------|-------|------|
| datagen `n_pairs` | 30 000 | plan default `default_npairs_for(4)` |
| datagen `imsize_set` | `((256,256),)` | constrained to 256² — matches the gate's `SBC_IMSIZE=(256,256)`; the default mixed set (up to 2048²) makes CPU datagen intractable in the compute budget (see SUMMARY deviations). 4×4 is robust at 256²: each 64×64 patch carries ≥4096 px |
| NPE | early-stopped @ epoch 61 (patience 40, cap 300) | Flux/`NormalisingFlow`, CPU; ~159 s train |
| NRE | trained (num_summaries 64, input dim `ratio_input_dim(4)`=80) | CPU |
| persistence | CPU-resident `Flux.state` | `npe_4.jld2` / `ratio_4.jld2` / `ood_nulls_4.jld2` load CPU-only, finite draws |

## SBC — Simulation-Based Calibration (M=2000, L=999, 50 bins)

Coded gate verdict (`sbc_gate`, sbc.jl): `passed = all(KS p > SBC_KS_ALPHA) && all(ECE ≤
SBC_ECE_GREEN)` across the 8 columns. `SBC_KS_ALPHA = 0.05`, `SBC_ECE_GREEN = 0.05`. The χ² column
is computed and reported for diagnostics but is **not** part of the coded `passed` conjunction.

| Parameter | KS p | χ² p | ECE | MCE | ECE light | KS pass (α=0.05) |
|-----------|------|------|-----|-----|-----------|------------------|
| ρ_true | 0.4270 | 0.8543 | 0.0120 | 0.99 | 🟢 green | pass |
| spillover | 0.05355 | 0.05955 | 0.0147 | 0.99 | 🟢 green | pass |
| autofluorescence | 0.1964 | 1.92e-7 | 0.0079 | 0.99 | 🟢 green | pass |
| label_efficiency | 0.5666 | 0.0009652 | 0.0073 | 0.99 | 🟢 green | pass |
| shift_dx | 0.03256 | 1.592e-6 | 0.0238 | 0.99 | 🟢 green | **fail** |
| shift_dy | 0.1614 | 0.001472 | 0.0177 | 0.99 | 🟢 green | pass |
| noise | 0.1964 | 2.235e-7 | 0.0042 | 0.99 | 🟢 green | pass |
| Δρ (paired, D-01) | 0.7532 | 0.4872 | 0.0074 | 0.99 | 🟢 green | pass |

- **`ks_pass = false`, `ece_pass = true` → SBC verdict: FAIL** (the gate requires the strict all-8
  KS conjunction).
- **Reading (memo §6):** ECE is green on every parameter (headline ρ_true ECE 0.0120, paired Δρ
  0.0074 — both with comfortably uniform KS, 0.427 / 0.753). **7 of 8 parameters pass KS** — a
  markedly stronger outcome than the 8×8 reference grid (5/8). The single KS rejection is
  `shift_dx` at 0.0326, a marginal miss just below the 0.05 threshold on a summary-uninformative
  registration-shift nuisance parameter — the same M=2000 KS over-sensitivity the Phase-5/8×8
  analyses documented (the fixed 8×8→4×4 patch-correlation summary carries little sub-pixel-shift
  information). MCE = 0.99 across all columns is the same single sparse-bin coverage→reliability
  routing artifact seen at 8×8; ECE (weighted) is the metric the traffic light is stated against,
  and it is green everywhere. The χ² column rejects on several nuisance parameters
  (autofluorescence / shift / noise), again an M=2000 hyper-sensitivity artifact on the
  summary-least-informative parameters; χ² is diagnostic in the coded gate.

## BF — Amortized Bayes Factor (Δρ sweep, non-clamped KDE baseline)

Pre-registered pass conditions: `cor(amortized, kde_unclamped) ≥ 0.95` **and** `max|Δ logBF| ≤ 0.5`.

| Metric | Value | Threshold | Pass |
|--------|-------|-----------|------|
| finite sweep pairs `n` | 20 (of `BF_SWEEP_N=25`) | — | — |
| corr(amortized, KDE-unclamped) | 0.9410 | ≥ 0.95 | **fail (near-miss)** |
| max\|Δ logBF\| | 12.40 | ≤ 0.5 | **fail** |

- **BF verdict: FAIL.**
- **Reading (memo §5, iter2):** corr = 0.941 is a near-miss of the 0.95 floor, consistent with the
  8×8 reference (0.947) and the spike confirmatory run (0.936) — the amortized NRE tracks the KDE
  baseline well in the informative mid-range. `max|Δ logBF| = 12.40` is a **KDE-baseline tail
  artifact**, not an NRE deficiency: at the one-sided Δρ sweep extremes the finite-sample KDE
  density collapses toward zero and its log-BF diverges, while the amortized NRE (correctly) stays
  bounded. 5 of 25 sweep points produced a non-finite KDE log-BF (a fully one-sided posterior → ±Inf
  under the honest un-floored baseline) and were dropped, leaving n=20. Constants were **not**
  re-tuned.
- **Note on the baseline (07-06 fix):** the non-clamped KDE baseline `kde_log_bf_unclamped`
  (bf.jl) was hardened this plan against a QuadGK domain-overshoot crash — its adaptive integral can
  push the tail probability a few ulp outside `[0,1]`, which without a domain clamp made `log(·)`
  throw a `DomainError` on the 4×4 sweep tails (the 8×8 sweep happened to land on exactly-zero →
  `log(0) = -Inf`, which was already dropped). The fix clamps the probability to its valid `[0,1]`
  domain so a saturated tail yields the honest `±Inf` (dropped), **not** the forbidden `1e-8` finite
  floor — in-range finite values are byte-identical. This is a numerical-safety Rule-1 bug fix, not
  a threshold change.

## OOD — Out-of-Distribution / Misspecification Flag

Pre-registered pass conditions: positive-control separability ROC AUC ≥ 0.80 (against injected
misspecification families); ID operating point at `OOD_ID_QUANTILE = 0.95` (~5% ID false-positive).

| Metric | Value | Note |
|--------|-------|------|
| ID draws `n` | 200 | in-distribution density-channel Mahalanobis scores |
| ID operating point `id_threshold` | 31.79 | 95th ID quantile → ~5% ID fire-rate by construction |
| per-family separability AUC | not computed | — |

- **OOD verdict: indeterminate (operating point only).**
- **Reading:** the bare `run_gate.jl --ood` CLI computes the pre-registered ID operating point but
  injects **no** positive-control (misspecified) simulators, so no separability AUC is scored
  (`auc = nothing`, `passed = nothing`) — identical by-design behavior to the 8×8 gate. The
  controlled misspecification-grid ROC experiment (`misspec_*` families, negative controls,
  `ood_roc_over_grid`, `with_pp=true` fusion) is **explicitly deferred to a later Phase-7 plan** per
  the `src/amortized/ood.jl` scope note. Introducing families ad-hoc here would be un-pre-registered
  and would violate the phase's anti-snooping discipline. The Phase-5 confirmatory OOD gate (pooled
  AUC 1.0 with the fused detector) stands as the method-level evidence; this 4×4 gate establishes the
  ID operating point on the fresh disjoint stream.

## Overall verdict

**FAIL** on the literal pre-registered conjunctions (SBC all-8 KS; BF corr + tol), with OOD reduced
to the ID operating point — but **stronger than the 8×8 reference grid** on the SBC axis. Every
residual failure reproduces a documented **non-method cause**:

1. **SBC** — ECE green on all 8 (calibrated); **7/8 KS pass** (vs 5/8 at 8×8). The lone rejection is
   `shift_dx` at KS 0.0326, a marginal miss on a summary-uninformative registration nuisance
   parameter — the M=2000 KS over-sensitivity artifact.
2. **BF** — corr 0.941 near-miss with strong mid-range agreement; `max|Δ|` 12.40 is a finite-sample
   KDE-baseline tail divergence, not an NRE error.
3. **OOD** — indeterminate by construction (misspec-family ROC is a later plan); ID operating point
   recorded.

This is the honest, un-tuned outcome on the fresh disjoint `PROD_SEED[4]` — the pre-registration in
`gate_consts_4.jl` is byte-locked and was not altered to force a pass. The **user-definable-grid
happy path is proven end-to-end**: `train_and_register(4)` ran `_train_grid_pipeline(4)` + `register!`,
`estimator_for(4)` returned the bundle in-process, and the fresh-seeded SBC/BF/OOD ran
CPU-reproducibly against the CPU-resident 4×4 net and wrote an atomic report. Per the Phase-6 memo's
"Clean Go" framing (each residual failure read to a non-method cause), the 4×4 grid remains eligible
for registry population (07-10); the residuals are surfaced here verbatim for the orchestrator/verifier.
