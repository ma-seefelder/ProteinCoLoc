# 16×16 Ship-Gate Outcome (D-05, PROD-02) — SHIP-WITH-CAVEAT grid

Recorded CPU-reproducible SBC / BF / OOD ship-gate for the **16×16 fine grid**, run through the
per-grid machinery (`test/gate/run_gate.jl --grid 16 --sbc --bf --ood`) against the **frozen 16×16
pre-registration** `test/gate/gate_consts_16.jl` (committed at `7e2318b`, byte-unchanged through the
run).

- **Gate stream (anti-snooping, T-7-08):** `PROD_SEED[16] = 0xb906f369f6cacf91` — the fresh disjoint
  Philox draw, provably distinct from the spike training seed `NPE_MASTER_SEED = 0xC0FFEE`, the
  spike validation seed `VAL_MASTER_SEED = 0x5BC0FFEE`, the datagen seed
  `DEFAULT_MASTER_SEED = 0x1`, and the sibling `PROD_SEED[8] = 0x8b39fecd4e2bcceb` /
  `PROD_SEED[4] = 0x8c0ad97b99bd6031` (the grid index G=16 keys a different Philox counter).
- **Reproducibility (T-7-07):** every inference `use_gpu = false` against the CPU-resident frozen
  net (`load_estimator`); the gate was **not** launched with `-t auto` (that is datagen-only,
  determinism-pinned here).
- **Report artifact:** `artifacts/grid_16/gate_report_16.jld2` (gitignored per repo policy —
  regenerable cache; the recorded numbers below are the committed deliverable).

## Trained bundle under test

The 16×16 bundle was produced by `_train_grid_pipeline(16)` (the shared 07-03 per-grid pipeline) in a
**prior session** (artifacts timestamped 2026-07-20 18:22); this plan resumed at the gate run. The
knobs below are read **directly from the persisted artifact metadata** — the values that are not
persisted are marked as such rather than reconstructed.

| Knob | Value | Source |
|------|-------|--------|
| datagen `n_pairs` | 80 000 | `npe_16.jld2` `meta.n_pairs` |
| training device | CPU (`use_gpu = false`) | `npe_16.jld2` `meta.use_gpu` |
| NPE architecture | `d_in = 512`, `D = 7`, `dstar = 64`, depth 3 × width 256, 10 coupling layers, flow 2×128 | `npe_16.jld2` `arch` |
| NRE | `ratio_n = 80 000`, `num_summaries = 64` | `ratio_16.jld2` `meta` + gate runtime log |
| OOD nulls | fitted on `n_train = 68 000` TRAIN-ONLY ID rows | `ood_nulls_16.jld2` `meta` |
| persistence | CPU-resident `Flux.state` | all three load CPU-only, finite draws |
| datagen `imsize_set` | **not recorded** in artifact metadata | `_train_grid_pipeline` persists only `grid`/`n_pairs`/`use_gpu` (pipeline.jl:200–202); the plan specified `default_imsize_for(16)` (≥512²-biased) |
| NPE/NRE early-stop epochs, wall time | **not recorded** | no training log retained |

`summary_dim(16) = 512`, `cont_rows(16) = 256`, `ratio_input_dim(16) = 1280` — the `d_in = 512` in the
persisted `arch` confirms the bundle is genuinely the 16×16 grid, not a mis-keyed coarse-grid net.

## SBC — Simulation-Based Calibration (M=2000, L=999, 50 bins, `SBC_IMSIZE = (512,512)`)

Coded gate verdict (`sbc_gate`, sbc.jl): `passed = all(KS p > SBC_KS_ALPHA) && all(ECE ≤
SBC_ECE_GREEN)` across the 8 columns. Pre-registered `SBC_KS_ALPHA = 0.05`, `SBC_ECE_GREEN = 0.05`.
The χ² column is computed and reported for diagnostics but is **not** part of the coded `passed`
conjunction.

Note the raised gate image size: `SBC_IMSIZE = (512,512)` (vs `(256,256)` for the 8×8 and 4×4 gates)
is itself part of the frozen 16×16 pre-registration — each 32×32 patch carries 1024 px, clearing the
≥15-survivor floor.

| Parameter | KS p | threshold | KS pass | χ² p | ECE | threshold | MCE | ECE light |
|-----------|------|-----------|---------|------|-----|-----------|-----|-----------|
| ρ_true | 0.000703 | > 0.05 | **fail** | 0.01483 | 0.02403 | ≤ 0.05 | 0.99 | 🟢 green |
| spillover | 0.002749 | > 0.05 | **fail** | 0.001024 | 0.01068 | ≤ 0.05 | 0.99 | 🟢 green |
| autofluorescence | 0.3363 | > 0.05 | pass | 5.904e-5 | 0.008921 | ≤ 0.05 | 0.99 | 🟢 green |
| label_efficiency | 7.149e-13 | > 0.05 | **fail** | 1.208e-9 | 0.01134 | ≤ 0.05 | 0.99 | 🟢 green |
| shift_dx | 0.01671 | > 0.05 | **fail** | 0.001438 | 0.01458 | ≤ 0.05 | 0.99 | 🟢 green |
| shift_dy | 0.1964 | > 0.05 | pass | 0.002462 | 0.008553 | ≤ 0.05 | 0.99 | 🟢 green |
| noise | 0.2835 | > 0.05 | pass | 1.033e-5 | 0.008184 | ≤ 0.05 | 0.99 | 🟢 green |
| Δρ (paired, D-01) | 0.1614 | > 0.05 | pass | 0.8971 | 0.005974 | ≤ 0.05 | 0.99 | 🟢 green |

- **`ks_pass = false`, `ece_pass = true` → SBC verdict: FAIL.**
- **KS pass count: 4 of 8** — weaker than 4×4 (7/8) and weaker than the 8×8 reference (5/8).
- **The material finding, stated plainly: `ρ_true` itself FAILS KS at p = 7.0e-4.** This is *not* the
  documented "M=2000 over-sensitivity on a summary-uninformative nuisance parameter" reading that
  covered the 8×8 and 4×4 residuals. `ρ_true` is the **headline** parameter and it is the one the
  colocalization claim rests on. At 8×8 its KS p was 0.567 and at 4×4 0.427 — comfortably uniform;
  at 16×16 it is rejected. Its ECE also roughly doubles (0.0240 vs 0.0120 at 4×4 and 0.0062 at 8×8),
  the largest ECE in the whole table, though still inside the green band. `spillover` likewise flips
  from clearly-uniform at both coarse grids (0.641 / 0.054) to rejected here (0.0027).
- **What survives:** the paired **Δρ** column — the quantity the shipped `delta_rho` accessor
  actually reports — remains uniform (KS 0.161, χ² 0.897) with the lowest ECE in the table
  (0.00597). ECE is green on all 8 columns. `MCE = 0.99` across every column is the same single
  sparse-bin coverage→reliability routing artifact recorded at 8×8 and 4×4; ECE (count-weighted) is
  the metric the traffic light is stated against.
- **Honest reading:** raising the training/gate image size to 512² was the mitigation for the
  fine-grid pixel budget, and it did not rescue the marginal-parameter calibration. The most
  economical explanation consistent with the 07-RESEARCH feasibility analysis is that 512 summary
  dimensions estimated from 256 patches of ~1024 px each is a materially noisier regression problem
  than 128 dims from 64 patches — the flow is honest on the *derived* Δρ contrast but drifts on the
  marginal ρ_true rank. That is a hypothesis, **not** a demonstrated non-method cause, and it is not
  offered as an excuse for the failure.

## BF — Amortized Bayes Factor (Δρ sweep, non-clamped KDE baseline)

Pre-registered pass conditions: `cor(amortized, kde_unclamped) ≥ BF_CORR_MIN = 0.95` **and**
`max|Δ logBF| ≤ BF_LOGBF_TOL = 0.5`, over `BF_SWEEP_N = 25` points on `[BF_SWEEP_LO, BF_SWEEP_HI] =
[-0.6, 0.8]`.

| Metric | Value | Threshold | Pass |
|--------|-------|-----------|------|
| finite sweep pairs `n` | 15 (of `BF_SWEEP_N = 25`) | — | — |
| corr(amortized, KDE-unclamped) | 0.9150 | ≥ 0.95 | **fail** |
| max\|Δ logBF\| | 12.72 | ≤ 0.5 | **fail** |

- **BF verdict: FAIL.**
- **Reading:** corr = 0.9150 is the **weakest** of the three grids (8×8: 0.9472, 4×4: 0.9410, spike
  confirmatory: 0.936) — the near-miss framing that applied at the coarse grids is thinner here.
  `max|Δ logBF| = 12.72` reproduces the documented KDE-baseline tail divergence (at the one-sided Δρ
  sweep extremes the finite-sample KDE density collapses and its log-BF diverges while the amortized
  NRE stays bounded); 10 of 25 sweep points produced a non-finite KDE log-BF under the honest
  un-floored baseline and were dropped, leaving n=15. Constants were **not** re-tuned.

## OOD — Out-of-Distribution / Misspecification Flag

Pre-registered pass conditions: positive-control separability ROC AUC ≥ `OOD_AUC_MIN = 0.80`
(against injected misspecification families); ID operating point at `OOD_ID_QUANTILE = 0.95`
(~5% ID false-positive).

| Metric | Value | Threshold | Note |
|--------|-------|-----------|------|
| ID draws `n` | 200 | — | in-distribution density-channel Mahalanobis scores |
| ID operating point `id_threshold` | 378.44 | — | 95th ID quantile → ~5% ID fire-rate by construction |
| per-family separability AUC | not computed (`nothing`) | ≥ 0.80 | no positive-control simulator injected |

- **OOD verdict: NOT RUN / INCONCLUSIVE (operating point only) — this is explicitly NOT a pass.**
  `passed = nothing` is an *absent* verdict, not a green one: there is no evidence in this gate that
  the 16×16 misspecification flag separates OOD from ID inputs.
- **Reading:** the bare `run_gate.jl --ood` CLI computes the pre-registered ID operating point but
  injects **no** positive-control (misspecified) simulators, so no separability AUC is scored
  (`auc = nothing`, `passed = nothing`) — identical by-design behavior to the 8×8 and 4×4 gates. The
  controlled misspecification-grid ROC experiment (`misspec_*` families, negative controls,
  `ood_roc_over_grid`, `with_pp = true` fusion) is **explicitly deferred to a later Phase-7 plan** per
  the `src/amortized/ood.jl` scope note. Introducing families ad-hoc here would be un-pre-registered
  and would violate the phase's anti-snooping discipline. The `id_threshold` of 378.44 (vs 103.86 at
  8×8, 31.79 at 4×4) scales as expected with the Mahalanobis dimension (`cont_rows(16) = 256`).

## Minimum-image-size caveat (MUST be carried by the 07-10 registry entry)

Per 07-RESEARCH §FEASIBILITY VERDICT, 16×16 is **SHIP-WITH-CAVEAT** and the caveat is a hard
physical constraint, independent of this gate's verdict:

> **16×16 is informative only on images ≥512².** `patch(img, 16)` tiles into 16×16 patches of
> `(W÷16)×(H÷16)` px, and `correlation(...)` sets a patch to `missing` when ≤15 pixels survive
> `_exclude_zero` (background exclusion). Pixels per patch: **256 at 256² — marginal after
> background exclusion**; 1024 at 512² (good); 4096 at 1024²; 5504 at the 1376×1028 real anchor.
> Applying a 16×16 estimator to 256² inputs therefore risks a summary vector dominated by `missing`
> patches, which the estimator will silently impute — producing a confident-looking posterior driven
> by the imputation, not the data.

This gate itself ran at the raised `SBC_IMSIZE = (512,512)` precisely because of this constraint, so
the numbers above are the **best case** for 16×16, not a typical-image case. Any registry entry for
grid 16 must surface a minimum-image-size warning (≥512², recommended ≥1024²).

## Overall verdict

**FAIL** — on every scored component, against the byte-locked `gate_consts_16.jl`:

| Component | Verdict | Headline |
|-----------|---------|----------|
| SBC (M=2000, L=999) | **FAIL** (`ks_pass = false`, `ece_pass = true`) | ECE green on all 8, but **KS 4/8** and **ρ_true itself rejected at p = 7.0e-4** |
| BF (n = 15 finite / 25) | **FAIL** | corr 0.9150 < 0.95; max\|Δ logBF\| 12.72 > 0.5 |
| OOD (n = 200) | indeterminate | ID operating point 378.44; per-family AUC deferred to the later misspec-ROC plan |

### Registry eligibility (07-10): **NOT eligible**

16×16 is **NOT** eligible for registry population in 07-10 on this evidence. The reason is specific,
and it is *not* merely "the literal conjunction failed" — 8×8 and 4×4 also failed literally and were
carried forward under the Phase-6 memo's "Clean Go" framing, which requires **every residual failure
to read to a documented non-method cause**. That framing does not extend to this gate:

1. **The headline parameter is miscalibrated.** `ρ_true` KS p = 7.0e-4 (and `spillover` 0.0027).
   Both were comfortably uniform at 8×8 and 4×4. The "M=2000 KS over-sensitivity on
   summary-uninformative nuisance parameters" explanation covered `label_efficiency` / `shift_*` /
   `noise`; it does **not** cover ρ_true, the parameter the colocalization claim is stated in. There
   is no documented non-method cause for this rejection.
2. **The BF axis degrades rather than holds.** corr 0.9150 is below both sibling grids and below the
   spike confirmatory run — the "near-miss of 0.95" reading is materially weaker here.
3. The `max|Δ logBF|` tail divergence and the indeterminate OOD **are** documented non-method causes
   and are not counted against 16×16.

Recommendation for 07-10: leave grid 16 **out** of `_SHIPPED_GRIDS`, or ship it only behind an
explicit user sign-off that also carries the minimum-image-size caveat above. Do not populate it by
default on the strength of the 8×8/4×4 precedent — the failure mode is different in kind.

This is the honest, un-tuned outcome on the fresh disjoint `PROD_SEED[16] = 0xb906f369f6cacf91`.
`test/gate/gate_consts_16.jl` is byte-identical to commit `7e2318b` and was not altered before,
during, or after the run; no threshold, seed, M, L, or bin count was touched, and no iteration toward
a pass was attempted. Per this project's Phase-5 precedent, a recorded negative is a valid and
publishable outcome.
