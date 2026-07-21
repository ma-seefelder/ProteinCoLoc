# 16×16 Ship-Gate Outcome (D-05, PROD-02)

Recorded CPU-reproducible SBC / BF / OOD ship-gate for the **16×16 fine grid — SHIP-WITH-CAVEAT**,
run through the per-grid machinery (`test/gate/run_gate.jl --grid 16 --sbc --bf --ood`) against the
**frozen 16×16 pre-registration** `test/gate/gate_consts_16.jl` (committed 7e2318b, *before* the run).

- **Gate stream (anti-snooping, T-7-08):** `PROD_SEED[16] = 0xb906f369f6cacf91`
  (= 13 332 611 383 314 534 289) — the fresh disjoint Philox draw, provably distinct from the spike
  training seed `NPE_MASTER_SEED = 0xC0FFEE`, the spike validation seed `VAL_MASTER_SEED = 0x5BC0FFEE`,
  the datagen seed `DEFAULT_MASTER_SEED = 0x1`, and the sibling `PROD_SEED[8] = 0x8b39fecd4e2bcceb` /
  `PROD_SEED[4] = 0x8c0ad97b99bd6031` (the grid index G=16 keys a different Philox counter).
- **Reproducibility (T-7-07):** every inference `use_gpu = false` against the CPU-resident frozen net
  (`load_estimator`); the gate was **not** launched with `-t auto` (that is datagen-only —
  determinism-pinned here). Wall clock ≈ 16.5 min for the full `--sbc --bf --ood` run.
- **Report artifact:** `artifacts/grid_16/gate_report_16.jld2`, `status = :ran` (gitignored per repo
  policy — regenerable cache; the recorded numbers below are the committed deliverable).

## Minimum-image-size caveat (the registry entry MUST carry this)

**16×16 is informative only on images ≥ 512².** At the coarse grids' `SBC_IMSIZE = (256,256)` each
16×16 patch spans only 256 px, which is marginal after background exclusion — the ≥15-survivor floor
is not reliably cleared. The 16×16 pre-registration therefore raises `SBC_IMSIZE` to **(512, 512)**
(each 32×32 patch carries 1024 px), and this is the larger reported simulate_pair image size relative
to the 8×8 and 4×4 gates *by design*. Any registry/user-facing entry for grid 16 must state the
≥512² minimum-image-size requirement; applying the 16×16 estimator to 256² data is out of its
validated domain.

## Trained bundle under test

`_train_grid_pipeline(16)` (the shared 07-03 per-grid pipeline), datagen launched under `julia -t auto`.

| Knob | Value | Source |
|------|-------|--------|
| datagen `n_pairs` | 80 000 | `npe_16.jld2` meta (`= default_npairs_for(16)`) |
| train/val split | 68 000 train / 12 000 val | `ood_nulls_16.jld2` meta `n_train` (`val_frac = 0.15`) |
| datagen `imsize_set` | **not recorded in the persisted metadata** | `_train_grid_pipeline` does not persist `imsize_set`; the pipeline default `default_imsize_for(16)` is the ≥512²-biased set `((512,512),(1024,1024),(1376,1028),(2048,2048))` (Pitfall 2). This is a *provenance gap* — the realized training image distribution cannot be re-derived from the artifacts alone. |
| NPE arch | `d_in = 512`, `D = 7`, `dstar = 64`, depth 3 × width 256, 10 coupling layers (flow 2×128) | `npe_16.jld2` arch |
| NPE device | `use_gpu = false` (CPU) | `npe_16.jld2` meta |
| NRE | `input_dim = 1280` (`= ratio_input_dim(16) = 5·16²`), `num_summaries = 64`, summary width 256, `ratio_n = 80 000` | `ratio_16.jld2` |
| persistence | CPU-resident `Flux.state` | `npe_16.jld2` / `ratio_16.jld2` / `ood_nulls_16.jld2` load CPU-only, finite draws |

## SBC — Simulation-Based Calibration (M = 2000, L = 999, 50 bins, `SBC_IMSIZE = (512,512)`)

Coded gate verdict (`sbc_gate`, sbc.jl): `passed = all(KS p > SBC_KS_ALPHA) && all(ECE ≤ SBC_ECE_GREEN)`
across the 8 columns. `SBC_KS_ALPHA = 0.05`, `SBC_ECE_GREEN = 0.05`. The χ² column is computed and
reported for diagnostics but is **not** part of the coded `passed` conjunction.

| Parameter | KS p | χ² p | ECE | MCE | ECE light | KS pass (α = 0.05) |
|-----------|------|------|-----|-----|-----------|--------------------|
| ρ_true | 0.000703 | 0.01483 | 0.02403 | 0.99 | 🟢 green | **fail** |
| spillover | 0.002749 | 0.001024 | 0.01068 | 0.99 | 🟢 green | **fail** |
| autofluorescence | 0.3363 | 5.904e-5 | 0.008921 | 0.99 | 🟢 green | pass |
| label_efficiency | 7.149e-13 | 1.208e-9 | 0.01134 | 0.99 | 🟢 green | **fail** |
| shift_dx | 0.01671 | 0.001438 | 0.01458 | 0.99 | 🟢 green | **fail** |
| shift_dy | 0.1964 | 0.002462 | 0.008553 | 0.99 | 🟢 green | pass |
| noise | 0.2835 | 1.033e-5 | 0.008184 | 0.99 | 🟢 green | pass |
| Δρ (paired, D-01) | 0.1614 | 0.8971 | 0.005974 | 0.99 | 🟢 green | pass |

- **`ks_pass = false`, `ece_pass = true` → SBC verdict: FAIL** (the gate requires the strict all-8 KS
  conjunction). `sbc.passed = false`.
- Recorded caption from the report: *"calibrated under the simulator; pair with the OOD result"*.

## BF — Amortized Bayes Factor (Δρ sweep, non-clamped KDE baseline)

Pre-registered pass conditions: `cor(amortized, kde_unclamped) ≥ BF_CORR_MIN = 0.95` **and**
`max|Δ logBF| ≤ BF_LOGBF_TOL = 0.5`.

| Metric | Value | Threshold | Pass |
|--------|-------|-----------|------|
| finite sweep pairs `n` | 15 (of `BF_SWEEP_N = 25`) | — | — |
| corr(amortized, KDE-unclamped) | 0.9150 | ≥ 0.95 | **fail** |
| max\|Δ logBF\| | 12.72 | ≤ 0.5 | **fail** |

- **BF verdict: FAIL** (`corr_pass = false`, `tol_pass = false`, `passed = false`).
- Constants were **not** re-tuned. 10 of 25 sweep points produced a non-finite KDE log-BF under the
  honest un-floored baseline and were dropped, leaving n = 15.

## OOD — Out-of-Distribution / Misspecification Flag

Pre-registered pass conditions: positive-control separability ROC AUC ≥ `OOD_AUC_MIN = 0.80`;
ID operating point at `OOD_ID_QUANTILE = 0.95` (~5 % ID false-positive rate).

| Metric | Value | Note |
|--------|-------|------|
| ID draws `n` | 200 | in-distribution density-channel Mahalanobis scores |
| ID operating point `id_threshold` | 378.44 | 95th ID quantile → ~5 % ID fire-rate by construction |
| separability AUC | `nothing` | **not computed** |
| OOD `passed` | `nothing` | **no verdict** |

- **OOD verdict: NOT RUN / INCONCLUSIVE. This is NOT a pass.** The gate produced `auc = nothing` and
  `passed = nothing` — there is no evidence here that the 16×16 OOD flag separates misspecified from
  in-distribution data. Only the ID operating point was established.
- **Cause:** the bare `run_gate.jl --ood` CLI computes the pre-registered ID operating point but the
  invocation injected **no** `pos_sim` positive control, so `ood_gate` had nothing to score an AUC
  against. `run_gate.jl` supports the `pos_sim` seam (07-04); it was simply not supplied.
- **Recorded as a gap:** the 07-07 plan's Task 3 asked for a `--ood` verdict, and by not injecting a
  positive control the run cannot deliver one. The same absence occurred at 8×8 and 4×4, where it was
  attributed to the deliberate deferral of the controlled misspecification-grid ROC experiment
  (`misspec_*` families, negative controls, `ood_roc_over_grid`, `with_pp = true` fusion) to a later
  Phase-7 plan per the `src/amortized/ood.jl` scope note. That deferral is a legitimate reason for the
  absence, but it does not convert the absence into a pass: **the 16×16 OOD gate is unscored** and must
  be scored before the 16×16 registry entry can claim a validated misspecification flag. The Phase-5
  confirmatory OOD result (pooled AUC 1.0 with the fused detector) is method-level evidence on the
  spike 8×8 net only — it is not a 16×16 result.

## Overall verdict — FAIL (literal pre-registration)

| Gate | Coded verdict | Detail |
|------|---------------|--------|
| SBC | **FAIL** | `ks_pass = false` (4/8 KS rejections), `ece_pass = true` |
| BF | **FAIL** | `corr_pass = false` (0.9150 < 0.95), `tol_pass = false` (12.72 > 0.5) |
| OOD | **NO VERDICT** | `auc = nothing`, `passed = nothing` — unscored, not a pass |

The 16×16 grid **does not clear its frozen pre-registration**. `test/gate/gate_consts_16.jl` is
byte-locked and was not altered, relaxed, or reinterpreted to force a pass, and no threshold was
re-tuned after seeing the numbers. On the SBC axis this is the **weakest of the three gated grids**
(4/8 KS pass, vs 7/8 at 4×4 and 5/8 at 8×8), including a rejection on the headline `ρ_true`
(KS p = 7.0e-4) — which the 4×4 and 8×8 grids both passed. On the BF axis corr = 0.9150 is the
**lowest** of the three (4×4 0.9410, 8×8 0.9472) and is no longer a marginal near-miss of the 0.95
floor.

## Interpretation (separate from the verdict above)

The following is analysis, **not** a re-scoring of the gate. The verdict is FAIL as recorded.

1. **ECE is green on all 8 parameters** (max 0.0240 on ρ_true, Δρ at 0.0060) — on the weighted
   reliability metric the 16×16 net is well-calibrated, and the paired Δρ column is uniform on both
   KS (0.161) and χ² (0.897). The strict all-8 KS conjunction is what rejects. `MCE = 0.99` across all
   columns is the same single sparse-bin coverage→reliability routing artifact seen at 4×4 and 8×8;
   ECE (weighted) is the metric the traffic light is stated against.
2. **The KS rejections concentrate where expected but now include ρ_true.** `label_efficiency`
   (7.1e-13) and `shift_dx` (0.0167) reproduce the documented Phase-5 / 8×8 pattern — M = 2000 KS
   hyper-sensitivity on the parameters the fixed patch-correlation summary carries least information
   about. `ρ_true` at 7.0e-4 and `spillover` at 0.0027 are **new** and are not explained by that
   pattern: they are consistent with the 16×16 summary (512-dim, 256 patches at 32×32 px on 512²
   images) being noisier per patch than the coarse grids, i.e. a genuine grid-resolution effect
   rather than a test-power artifact. This should be read as a real weakness of 16×16, not dismissed.
3. **`max|Δ logBF| = 12.72` carries the known clamped-KDE tail signature** documented in Phase 5
   (iter2) and Go/No-Go memo §5: at the one-sided Δρ sweep extremes the finite-sample KDE density
   collapses and its log-BF diverges while the amortized NRE stays bounded; 10/25 sweep points went
   non-finite and were dropped. That reading applies to `max|Δ|`. It does **not** explain
   corr = 0.9150, which is computed on the 15 *finite* pairs and is the lowest of the three grids —
   the mid-range agreement is genuinely weaker at 16×16.
4. **OOD is simply missing**, see above. No interpretation can substitute for the unscored AUC.

## Recommendation for 07-10 (registry population)

Unlike 4×4 and 8×8 — where every residual failure read to a documented non-method cause and the
grids stayed eligible under the memo's "Clean Go" framing — the 16×16 residuals are **not fully
attributable to pre-registration/baseline design**: the ρ_true KS rejection and the corr = 0.9150
mid-range BF gap are grid-resolution effects. Combined with the unscored OOD gate and the
unrecorded training image distribution, the honest position is:

- **16×16 should be registered, if at all, only as SHIP-WITH-CAVEAT and explicitly flagged as
  weaker than 4×4/8×8**, carrying (a) the ≥512² minimum-image-size requirement, (b) the ρ_true
  SBC-KS rejection, (c) the un-validated OOD flag.
- The decision is deferred to 07-10 / the phase verifier; this document states the evidence, not the
  disposition.
