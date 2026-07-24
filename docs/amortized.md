# Amortized colocalization inference (v2.0 "AmortizedColoc")

v2.0 is a **breaking release** (D-01). The primary public API is amortized simulation-based
inference: a trained neural posterior estimator (NPE) and neural ratio estimator (NRE) produce a
**calibrated colocalization posterior, a Δρ contrast, an amortized Bayes factor, and an
out-of-distribution (OOD) flag in milliseconds per dataset** — a single forward pass, not the
per-dataset ADVI minutes of v1.0. The Turing/ADVI path remains as an **internal, weakdep-gated
reference implementation** (`ext/ProteinCoLocTuringExt.jl`) and is deliberately **not exported**;
`colocalization` / `compute_BayesFactor` are internal only.

> **Read the "Named limits" section below before using this tool for a scientific claim.** The v2.0
> GO is honest but **conditional** ("GO with named limits", `07-GO-NO-GO-UPDATE.md`, Option A). The
> limits are non-negotiable and must accompany any published result.

## Public entry point

```julia
using ProteinCoLoc

r = colocalization_amortized(img, control, channels;
                             num_patches = 8, N = 2000, use_gpu = false)
```

- `img::MultiChannelImage`, `control::MultiChannelImage` — the sample and control images. Δρ is the
  Monte-Carlo difference `ρ(img) − ρ(control)` (D-03): two independent single-stack posterior passes.
- `channels::AbstractVector{<:Integer}` — the **two** channels to compare, e.g. `[1, 2]`.
- `num_patches` — the registry key (see below). **v2.0 ships only `8`.**
- `N` — posterior draws per stack (default 2000).
- `use_gpu` — CPU-reproducible shipped default `false` (D-06); a plumbed opt-in only.

Returns an `AmortizedColocResult <: AbstractColocResult`. Read it through the shared **accessor
interface** (dispatch on the interface, never on concrete fields — D-02):

| Accessor | Returns |
|---|---|
| `posterior_draws(r)` | `7×N` physical-θ posterior draws (row 1 = `ρ_true`) |
| `delta_rho(r)` | length-`N` Δρ (sample − control) draw cloud |
| `bayes_factor(r)` | scalar **log** Bayes factor (coloc vs null), from the NRE |
| `is_ood(r)` | `Bool` — OOD / misspecification flag |

`r.calibration::CalibrationMeta` carries the grid's ship-gate provenance; `r.meta` carries the run
metadata (`N`, `channels`, `num_patches`, `use_gpu`).

An unregistered / unshipped `num_patches` raises an `ArgumentError` that names the shipped grids and
points at `train_and_register` — **never a silent fallback to a default grid** (T-7-05).

## The grid registry (`num_patches` as key, PROD-02)

Each patch grid `G` has its own calibrated `EstimatorBundle` (NPE + NRE + OOD nulls + frozen
standardizers), keyed by `num_patches` in the registry:

- `estimator_for(G)` — return the bundle for grid `G`: a registered bundle, else a **lazy-loaded**
  shipped-grid artifact, else a clear error.
- `train_and_register(G)` — the documented on-ramp for **any other grid**: runs the full per-grid
  pipeline (`_train_grid_pipeline`) and registers the result. This is real, user-runnable code — a
  grid outside the shipped family is a training run away, not a hard-coded limit.

### Shipped grids: **`{8}` only**

`_SHIPPED_GRIDS = (8,)`. The v2.0 GO rests **entirely on the 8×8 reference grid**. The other grids
are excluded, each for a recorded reason:

| Grid | Status | Reason |
|---|---|---|
| **8×8** | **SHIPPED** | The GO rests on it (targets calibrated, BF sim-validated, OOD passes). |
| 4×4 | excluded | Never post-hoc re-analysed with the corrected SBC/BF statistics the GO relies on (`gate-4x4.md`). |
| 16×16 | excluded | Gate **FAILED** with no non-method cause for the ρ_true rejection (`gate-16x16.md`). |
| 32×32 | excluded | **CAPPED** — never trained; a gate would run the apparatus known defective (`gate-32x32.md`). |
| 64×64 | dropped | D-04 (finer grids lengthen the input but the NPE output stays a single global ρ). |

### Content-hashed, lazy artifact loading

The shipped 8×8 bundle is **not committed to the repository**. It is resolved from a
content-hashed `Artifacts.toml` entry (`git-tree-sha1`, `lazy = true`), hosted as a GitHub-Release
asset. `_lazy_load_from_artifact!(8)`:

1. uses the artifact store if the content-addressed artifact is already installed (verified at
   install time);
2. else verifies the in-repo dev bundle's **git-tree-sha1 against `Artifacts.toml`** and loads it
   (dev/CI, no network);
3. else lazily downloads the Release asset and verifies it.

Every path is **tree-sha1 verified before load** (T-7-01 — never load a content-mismatched model).
The shipped bundle is the retrained **bounded-θ / realistic-image-size** 8×8 net the GO rests on
(`artifacts/amended_v2/grid_8/`).

## CPU-reproducible by default (D-06)

Every shipped inference call defaults `use_gpu = false`. Training may use a GPU
(`train_and_register` defaults `use_gpu = has_cuda_device()`), but the persisted nets are
CPU-resident (`Flux.state`) and the pre-registered ship-gate + shipped inference run CPU-only, so
what is pre-registered and what a user reproduces is **device-independent**. GPU degrades gracefully
to CPU when CUDA is absent.

## Windowed sub-tile local map (coarse)

`local_coloc_map(img, control, channels; grid = 8, tiles = (r, c))` runs the frozen, already-gated
bundle over image sub-tiles and returns a `LocalColocMap` (a per-tile Δρ matrix + per-tile OOD
flag). It is a **coarse windowed readout, not a calibrated per-region posterior**: each tile is an
independent global-ρ read on a smaller window, with no spatial prior, no borrowing of strength
between neighbours, and no per-region uncertainty. Degenerate tiles (below the grid size, below the
≥15-survivor floor, or a non-finite read) take a finite sentinel (`0.0`) **and** a forced OOD flag,
so "unscorable" can never read as "measured no difference" (T-7-04). The calibrated per-region
map (GP/CAR lattice, `SpatialColocResult`) is Phase 12 and is **not** implemented here.

## Named limits (non-negotiable honesty — carry these into the manuscript)

The v2.0 GO is **GO with named limits** (`07-GO-NO-GO-UPDATE.md §5`). These are properties of the
shipped 8×8 estimator and its calibration evidence; they must be stated with any published result.

1. **The gate was amended twice; the final grid-8 verdict is post-hoc, not a clean
   pre-registration.** Pre-registration credibility on grid 8 is spent, and the memo says so. The
   literal amended ship-gate (one shot, `PROD_SEED_V2`) **FAILED** (`gate-8x8-amended.md`); the GO
   rests on a corrected **post-hoc re-analysis** (`posthoc_reanalysis.jl`).

2. **ρ_true SBC calibration depends on randomized-rank atom handling.** The prior places mass on
   **non-realizable μ** (a documented simulator/prior inconsistency), which puts atoms in the ρ_true
   rank distribution and makes the raw KS test undefined (raw KS ≈ 5.9e-4 → **randomized-rank KS ≈
   0.72**, the textbook fix). Δρ is calibrated (KS ≈ 0.20). Holm over the targets: PASS.

3. **The Bayes factor is validated by simulation-based DISCRIMINATION, not magnitude agreement with
   a per-pair reference.** The KDE Bayes-factor baseline was shown **invalid past |logBF| ≈ log(L) ≈
   6.9** (attrition is 100% baseline-side; the NRE never goes non-finite). The BF is instead
   validated SBC-style (spike 014): **AUC(coloc vs null) = 0.994**, 0 attrition, monotone in evidence
   (Spearman 0.93), decision-calibration ECE 0.019, FPR 0.036 at logBF > 0. This is a concept-proof
   (single seed); full gate integration is deferred (Option B, `07-GO-NO-GO-UPDATE.md §4`).

4. **Nuisance-parameter marginals carry a ~0.08-SD residual drift** not certifiable as negligible at
   M=2000. Under a pre-registered δ=0.10-SD equivalence (TOST) test, 4/6 nuisances pass; 2/6
   (autofluorescence, label_efficiency) fail — not because the drift is large (point drift ~0.08 SD
   < 0.10) but because the 90% CI edges just past 0.10. **No model lever removes it** (capacity
   redistributes it, spike 012; μ-truncation fixes only the atoms at a real cost to |ρ|≥0.95
   inference and ADVI comparability, spike 013). This is a **marginal-reproduction limit on
   parameters the fixed patch-correlation summary cannot constrain** (5/8 columns are vacuous —
   posterior ≈ prior), **not** a failure of the coloc inference. It is a named calibration limit.

5. **The one-shot gate carries ±1.5-z per-parameter seed variance.** Single-seed verdicts are
   fragile; the recorded numbers are one draw of a noisy statistic.

6. **The image-size regime matters (F5).** The 8×8 calibration is demonstrated on the realistic
   image-size mixture the shipped net was trained on (`imsize_source = :generate_samples`, provenance
   recorded). A calibration claim does **not** automatically transfer to a different image-size
   regime; the shipped net records its training image-size distribution for exactly this reason.

7. **The shipped OOD flag is the summary-density (Mahalanobis) channel** at the recorded ID
   operating point (~5% ID false-positive rate). The **fused** detector that reached gate AUC 1.0
   also used a noise channel and a posterior-predictive channel; those are **gate-time constructs**
   (they need the raw image pair / re-simulation) and are not carried in the frozen bundle. The
   density channel alone has a documented blind spot on detector-noise misspecification. Known OOD
   blind spot: correlation-preserving affine transforms.

**Findings and scripts:** `Skill("spike-findings-proteincoloc")`, `07-GO-NO-GO-UPDATE.md`,
`.planning/phases/07-productionization-conditional-on-go/` (gate reports, `posthoc_reanalysis.jl`),
`.planning/spikes/006-014`.
