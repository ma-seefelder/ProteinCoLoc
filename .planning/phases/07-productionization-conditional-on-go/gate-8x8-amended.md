# Grid 8×8 — Amended Ship-Gate, Confirmatory Run

**Run:** 2026-07-21 17:11–17:42 (30.9 min), CPU-only, single binding invocation
**Net:** `artifacts/amended_v2/grid_8/` (bounded logit-θ space, trained on the realistic imsize mixture, 50k pairs)
**Pre-registration:** `test/gate/gate_consts_8_v2.jl`, frozen at commit `c69d381` — **unchanged between freezing and this run** (verified by `git diff`)
**Seed:** `PROD_SEED_V2` → 1135605683775656488 (fresh, asserted disjoint from all prior seeds)
**Protocol:** amendment §6 — ONE run, no retry, no re-seed, no threshold change. Honoured.

---

## Overall verdict: **FAIL**

| Arm | v2 (binding) | v1 (same data, comparison) |
|---|---|---|
| **SBC** | **FAIL** (`ks_pass=false`, `chi2_pass=false`, `ece_pass=true`) | **FAIL** (`ks_pass=false`, `chi2_pass=false`) |
| **BF** | **FAIL — `corr_verdict = :invalid`** | FAIL (`corr_pass=true`, `tol_pass=false`) |
| **OOD** | **PASS** (pooled AUC 1.0, ID fire rate 0.05) | — |

The gate FAILED. This is reported as the primary result; interpretation is separated below.

---

## 1. SBC — FAIL under **both** rule sets

This is the single most important line in this report: **the multiplicity correction did not rescue SBC.** The failure is not an artifact of the uncorrected all-8 conjunction (defect A1). It survives Holm–Bonferroni at FWER 0.05.

Holm-adjusted KS p-values — four rejections survive:

| parameter | KS p | Holm-adj p | χ² p | ECE | shrinkage | vacuous |
|---|---|---|---|---|---|---|
| **ρ_true** | **1.12e-7** | **8.96e-7** | 3.53e-26 | 0.0371 | **0.130** | no |
| **spillover** | 6.82e-6 | 4.77e-5 | 1.72e-4 | 0.0192 | 1.010 | **YES** |
| autofluorescence | 0.0761 | 0.228 | 0.765 | 0.0032 | 0.997 | **YES** |
| **label_efficiency** | **1.92e-4** | **1.15e-3** | 1.51e-3 | 0.0330 | 0.923 | no |
| shift_dx | 0.883 | 0.883 | 0.909 | 0.0028 | 0.976 | **YES** |
| **shift_dy** | 1.67e-3 | 8.37e-3 | 0.0433 | 0.0153 | 0.987 | **YES** |
| noise | 0.0192 | 0.0767 | 0.157 | 0.0172 | 1.025 | **YES** |
| Δρ | 0.196 | 0.393 | 0.0528 | 0.0043 | 0.136 | no |

### 1a. Five of eight parameters are VACUOUS

`vacuous_params = [spillover, autofluorescence, shift_dx, shift_dy, noise]` — shrinkage (post_sd/prior_sd) between 0.92 and 1.03, i.e. **the posterior is essentially the prior**. Per finding F3, a uniform rank distribution from such a parameter is *not* calibration evidence.

Only three parameters carry real information: **ρ_true (0.130)**, **Δρ (0.136)**, and marginally **label_efficiency (0.923)**.

Consequence for reading this table: `shift_dx`'s clean KS p = 0.883 is meaningless — it is a vacuous column. Conversely, `spillover`'s rejection at a shrinkage of 1.01 is a rejection on a column that learned nothing, which is itself odd and worth a look.

### 1b. ρ_true calibration BROKE under the realistic image-size regime

The decisive comparison, with the θ-space held constant:

| regime | ρ_true KS p |
|---|---|
| 256²-only, unbounded θ (original v1 gate, grid 8) | 0.567 |
| 256²-only, bounded θ (DEV diagnostic) | z = −1.25 (clean) |
| **realistic mixture, bounded θ (this run)** | **1.12e-7** |

Because the DEV diagnostic isolated the θ-space effect at 256² and found ρ_true unharmed, the degradation here is attributable to the **image-size regime**, not the bounded-θ fix. Finding **F5 is confirmed on the reported gate**: the calibration demonstrated at 256² does not transfer to realistic microscopy dimensions.

The 256² result was flattering. That is exactly what F5 warned about, and it is now measured rather than inferred.

### 1c. The bounded-θ fix did help, but not enough

`label_efficiency` KS p improved from **1.84e-9** (v1 gate, 256², unbounded) to **1.92e-4** (here) — roughly five orders of magnitude — while its shrinkage stayed informative (0.923, non-vacuous). The mechanism fix (F2) works. It is not sufficient to pass.

### 1d. Realised image-size mixture

`realised_imsize_counts` = 512²:1631, 1024²:1017, 1376×1028:961, 2048²:391 (n=4000 draws) against targets 1600/1000/1000/400. The mixture sampled as specified.

---

## 2. BF — FAIL, and the correlation verdict is `:invalid`, not merely failing

| quantity | value |
|---|---|
| n (finite pairs) | **58** |
| n attempted | 100 |
| `n_min` required | 69 |
| r̂ | 0.9516 |
| one-sided 95% lower bound | 0.9255 (needs ≥ 0.96394) |
| two-sided 95% CI | [0.9192, 0.9712] |
| q50 / q95 / max \|Δ logBF\| | 1.489 / **12.406** / 13.544 |

**42% attrition.** The gate attempted 100 pairs and obtained 58 finite ones, below the amendment's `n_min = 69` floor, so `corr_verdict = :invalid` — the run did not collect enough valid data to render a correlation judgement at the pre-registered precision. This is a **process failure needing investigation**, not a clean statistical fail, and it must not be reported as one. The amendment set `n_min` precisely so attrition could not silently return the gate to the underpowered regime; that guard fired as designed.

**Rules attribution (the free comparison):** under v1 rules the same data gives `corr_pass = true` (r̂ = 0.9516 > 0.95); under v2 it gives `corr_pass = false`. The amendment's stricter lower-confidence-bound rule is doing exactly what it was written to do — it refuses a point estimate that its own uncertainty does not support. Had we kept v1 rules, this run would have recorded a passing correlation on 58 pairs.

`q95|Δ logBF| = 12.4` against a tolerance of 0.5 — the tail disagreement with the KDE baseline persists and is not explained by clamping alone (this baseline is the non-clamped one).

---

## 3. OOD — PASS

| quantity | value |
|---|---|
| pooled AUC | **1.0** |
| per-family best AUC | texture 1.0, noise 1.0, optics 1.0, background 1.0 |
| ID fire rate | 0.05 (matches the pre-registered 1 − `OOD_ID_QUANTILE`) |
| ID threshold / fused threshold | 179.14 / 7.128 |

Honest detail behind the pooled 1.0:
- **Density channel alone scores AUC 0.0 on the `noise` family at every level** — worse than chance. The image-noise channel (1.0 at every level) is what carries it. The documented detector-noise blind spot of the correlation summary reproduces here.
- **`background` level 3 fused AUC = 0.0** — non-monotonic in level, as previously observed on grids 4/8/16. Reported as-is; the gate reads the best level.
- Density on `optics` is weak throughout (0.58–0.67); again the noise channel carries the family.

**Negative controls did NOT all pass:** `neg_pass` = block **false** (fire rate 0.10), affine true (0.0), rotate **false** (0.133). These are not folded into `passed` (they never were), but a 10–13% fire rate on transforms that should be summary-invariant is above the 5% ID rate and deserves scrutiny.

---

## 4. Limitation the amendment conceded, restated

Amendment §5: this run changes **both** the decision rules and the model relative to the original v1 gate, so an outcome cannot be attributed to either alone. Two mitigations were built in and both paid off:
- the **v1/v2 side-by-side on identical data** isolates the rules effect (decisive on the BF arm: v1 pass → v2 fail);
- the **DEV-seed θ-diagnostic at 256²** isolates the model effect (decisive on ρ_true: unharmed by bounded θ, therefore the regression here is the imsize regime).

Because SBC fails under *both* rule sets, no attribution ambiguity remains for the headline result: **the SBC failure is a property of the model/regime, not of the amended rules.**

---

## 5. What this run establishes

1. **The amended gate FAILED.** One run, as pre-registered. No retry, no re-seed, no threshold adjustment.
2. **F5 is confirmed, not merely suspected.** ρ_true calibration collapses (KS 0.567 → 1.12e-7) when the image-size regime moves from 256² to realistic microscopy dimensions. Any calibration claim for this tool must state the image-size regime it was demonstrated at.
3. **The bounded-θ fix (F2) works but is insufficient** — label_efficiency improved ~5 orders of magnitude while staying informative.
4. **5/8 parameters are vacuous.** The correlation-only summary constrains ρ_true, Δρ and (partially) label_efficiency; it learns essentially nothing about spillover, autofluorescence, shift_dx, shift_dy and noise. Any future SBC "pass" on those columns must be read as vacuous.
5. **The BF arm could not be evaluated at the pre-registered precision** (58 < 69 finite pairs). The attrition cause is unknown and is the most concrete open technical question.
6. **OOD passes** on the fused detector, with the density channel's noise-family blind spot and two failing negative controls documented.

## 6. Open items (none actioned in this run)

- Diagnose the 42% BF pair attrition (non-finite pairs) — blocks a valid BF verdict.
- Investigate why `spillover` rejects at shrinkage 1.01 (a rejection on a column that learned nothing).
- Negative controls `block` (0.10) and `rotate` (0.133) fire above the ID rate.
- ρ_true miscalibration under the realistic regime is now the primary scientific blocker for shipping.

**The original v1 pre-registration and its FAIL record remain intact and must be co-cited with this report** (`gate-8x8.md`, `gate-4x4.md`, `gate-16x16.md`, `gate_consts_{4,8,16}.jl`).
