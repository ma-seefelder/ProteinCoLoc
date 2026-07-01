<!--
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder — GNU AGPL v3 or later.
-->

# Summary-Statistic Choice — Ablation Justification (ABL-02, D-07)

**Deliverable:** the written justification consumed by Phase 5 (SBC / amortized BF / OOD).
**Decision:** **`:min`** (the minimal 128-dim patch-correlation summary).
**Rule:** pre-registered, rho_true-led, advisory-with-rule (D-07 / RESEARCH A3).
**Produced by:** `spike/npe/ablation.jl` — `ablate(...)` + `choose_summary(...)`.

---

## What was ablated (ABL-01, D-06)

The **same** NPE architecture (`build_estimator`: MLP summary net -> NeuralEstimators
`NormalisingFlow` over the 7-dim theta, rho_true row 1) was trained on the **two
already-cached summary variants** under the loader's leak-free **k=5 CV**, with **zero
re-simulation** — only the `load_fold(...; variant=...)` toggle differs:

| Variant | Dim | Content |
|---------|-----|---------|
| `:min`  | 128 | D-01 encoding: 64 imputed patch correlations + 64 binary present/absent mask rows |
| `:aug`  | 142 | `:min` **superset** + 14 continuous moments (Manders M1/M2, whole-image Pearson, patch-grid median/IQR/mean/std/skew/excess-kurtosis, fraction-missing, per-channel intensity median+IQR) |

Architecture, `master_seed`, fold membership, the leak-free theta-standardization recipe
and every training hyper-parameter were held **identical** across arms, so the comparison
isolates the **summary choice**, not the network (D-06). Per-parameter RMSE (all 7 theta,
rho_true headline) was scored on each held-out fold via NeuralEstimators `assess`/`rmse`,
scaled to physical units by `theta_zt.scale` (Pitfall 5). CPU-only (`use_gpu=false`, D-10).

**Reported run:** N=200 main-pool stacks, K=5, epochs=200, batch 64, N_draws=2000,
`master_seed = NPE_MASTER_SEED (0xC0FFEE)`, CPU-only.

---

## Results

### Headline — rho_true RMSE per fold (physical units)

| Fold | `:min` | `:aug` | aug beats min? |
|------|--------|--------|----------------|
| 1 | 0.1240 | 0.0889 | yes |
| 2 | 0.1059 | 0.1148 | no |
| 3 | 0.1148 | 0.0894 | yes |
| 4 | 0.1444 | 0.0830 | yes |
| 5 | 0.1201 | 0.1207 | no |
| **mean** | **0.1218** | **0.0993** | **3 / 5 folds** |

### All-7 theta mean per-parameter RMSE

| Parameter | `:min` | `:aug` |
|-----------|--------|--------|
| rho_true (headline) | 0.1218 | 0.0993 |
| spillover | 0.0598 | 0.0572 |
| autofluorescence | 0.0294 | 0.0290 |
| label efficiency | 0.1207 | 0.0765 |
| sub-pixel shift | 0.6078 | 0.6143 |
| noise | 0.5914 | 0.6003 |
| theta7 | 0.3004 | 0.3050 |

---

## The pre-registered decision rule (D-07)

`choose_summary` returns `:aug` **iff both**:

1. **Relative margin** — `rho_rmse_aug <= (1 - ABL_REL_MARGIN) * rho_rmse_min`
   with `ABL_REL_MARGIN = 0.05`; threshold `= (1 - 0.05) * 0.1218 = 0.1158`.
   -> aug `0.0993 <= 0.1158` => **TRUE** (aug improves rho_true by ~18% on the fold mean).
2. **Fold consistency** — `fold_wins_aug >= ABL_FOLD_CONSISTENCY` with
   `ABL_FOLD_CONSISTENCY = 4`.
   -> `fold_wins_aug = 3 < 4` => **FALSE** (aug wins only 3 of 5 folds).

Both conditions are required. Condition 2 fails, so the rule returns **`:min`**.

**Why this is the right call, not a technicality.** The augmented moments deliver a real
mean-RMSE gain, but that gain is **not consistent across folds** (aug is *worse* on folds
2 and 5). A summary that helps on some data partitions and hurts on others is a
fold-to-fold-**unstable** signal — the classic symptom of a few extra features letting the
flow fit partition-specific idiosyncrasy rather than the shared physics. The two-part rule
is designed exactly to reject "wins on average, unstable per fold": the margin alone would
have chosen `:aug`; the consistency guard vetoes it. Defaulting to `:min` under this
condition is the pre-registered, parsimony-preferring outcome (D-07), decided **before**
the run and not tuned to it.

**Insufficiency vs. calibration (D-07).** A downstream SBC-pass on `:min` must **not** be
read as "the summary is sufficient." Calibration (nominal coverage) and sufficiency (low
RMSE / high information) are distinct: an under-informative summary can be honestly
calibrated yet wide. RMSE is the sufficiency axis here, and `:min`'s rho_true RMSE (0.122)
is the number Phase 5 inherits as the accuracy floor — a SBC-pass-but-high-RMSE outcome is
an **insufficiency signal**, not a pass.

---

## OOD-detectability interaction — the explicit Phase-5 coupling (ABL-02)

The summary choice is not accuracy-only; it couples to Phase-5 OOD/misspecification power,
and that coupling is **why the tie/instability default is `:min`, not `:aug`**:

- **Parsimony => tighter, more orthogonal discrepancy channel.** The minimal summary is a
  pure patch-correlation vector (+ mask). Its discrepancy from the training distribution is
  concentrated in a single, well-understood channel, which an OOD flag can monitor with
  fewer degrees of freedom. The 14 augmented moments **add breadth** to the summary, and
  each extra dimension is another axis along which an in-distribution-looking-but-wrong
  input could "hide" by matching moments while the underlying physics is misspecified —
  diluting the per-dimension OOD signal.

- **Structural OOD blind spot — named, not hidden.** Both variants share a hard limit: the
  spike's summary is a **fixed 8x8 patch-correlation statistic** (CLAUDE.md constraint).
  Any misspecification that is **orthogonal to the patch-correlation (and moment) subspace**
  is **structurally undetectable** by an OOD flag built on this summary — the flag can only
  see discrepancies that project onto what the summary measures. `:aug` widens that
  measured subspace slightly (more moments => a few more detectable directions) but does
  **not** remove the blind spot; it trades a marginal, fold-unstable accuracy gain for a
  broader, harder-to-calibrate discrepancy space. Given the accuracy gain fails the
  consistency gate, the parsimonious `:min` is the better base for the Phase-5 OOD story.

**Consequence for Phase 5.** Phase 5 builds SBC, the amortized Bayes factor, and the OOD
flag on the **`:min`** summary and the `trained_npe.jld2` estimator. The OOD flag's power
is bounded to discrepancies expressible in the patch-correlation subspace; misspecifications
orthogonal to it will not trigger it. This bound is inherited, documented, and must be
stated honestly in the Go/No-Go memo (Phase 6).

---

## Reproduce

```julia
using Pkg; Pkg.activate("spike")
include("spike/data/generate.jl")
include("spike/npe/ablation.jl")
dir = generate_cache(mktempdir(); N=200, master_seed=NPE_MASTER_SEED, n_holdout=20)
r   = ablate(dir; master_seed=NPE_MASTER_SEED, K=5, N=2000, epochs=200, batchsize=64)
choose_summary(r)   # => :min
```

The SC4/SC5 gates in `spike/test/test_npe.jl` exercise `ablate`/`choose_summary` on a
small fixture cache (fast, deterministic); this document reports the full-scale run.
