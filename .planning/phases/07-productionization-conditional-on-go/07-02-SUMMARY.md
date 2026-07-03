---
phase: 07-productionization-conditional-on-go
plan: 02
subsystem: amortized read surfaces (NPE infer / NRE Bayes factor / OOD flag)
status: COMPLETE
tags: [amortized-inference, npe-read-surface, nre-bayes-factor, ood-flag, memo-5-hardening, cpu-reproducible, grid-general]
requirements: [PROD-01]
dependency_graph:
  requires:
    - 07-00 (AbstractColocResult/OODVerdict/CalibrationMeta hierarchy; green root Pkg.resolve; NeuralEstimators 0.2.1 pinned)
    - 07-01 (grid-parametric summary.jl _summary_row_partition + dimension helpers summary_dim/cont_rows/ratio_input_dim; registry skeleton)
  provides:
    - "src/amortized/infer.jl: standardize_summary, posterior_for, rho_hat, rho_draws, delta_rho (4-arg), interval_width (CPU-default read surface)"
    - "src/amortized/bf.jl: amortized_log_bf, pair_encode (grid-general 5G²), kde_log_bf_unclamped non-clamped baseline"
    - "src/amortized/ood.jl: fit_ood_nulls + maha_score (density), pp_mismatch_score + finite-guard (PP), noise channel, roc_auc/id_threshold/youden_j, ood_flag, ood_verdict (density+noise+PP OR-fusion)"
  affects:
    - "the future colocalization_amortized public entry point composes these read surfaces"
    - "the per-grid D-05 ship-gate (07-03+) consumes kde_log_bf_unclamped + the OOD channels + the SBC/BF harness"
tech_stack:
  added: [LinearAlgebra (stdlib dep, cholesky/Symmetric/I for the OOD Mahalanobis)]
  patterns: [CPU-default use_gpu split (D-06), frozen-stats reconstruct discipline (Pitfall 5), grid coupling via single-source _summary_row_partition, hand-rolled Mahalanobis+ROC (no new package), non-clamped KDE BF baseline (Memo §5), robust-z OR-fusion]
key_files:
  created:
    - src/amortized/infer.jl
    - src/amortized/bf.jl
    - src/amortized/ood.jl
    - .planning/phases/07-productionization-conditional-on-go/07-02-SUMMARY.md
  modified:
    - src/ProteinCoLoc.jl
    - test/runtests.jl
    - Project.toml
    - Manifest.toml
decisions:
  - "NPE/NRE/OOD read paths default use_gpu=false everywhere (D-06 CPU-reproducible shipped inference); use_gpu stays a plumbed keyword"
  - "pair_encode derives cont-rows nc=G² from the summary length (length(Zs)÷2), so the A7 difference encoding is grid-general (ratio_input_dim(G)=5G²) with no G argument"
  - "kde_log_bf_unclamped mirrors compute_BayesFactor KDE math WITHOUT the spike baseline's 1e-8 _clampp floor so max|Δ logBF| is artifact-free (Memo §5 / T-7-06); lives in core (KernelDensity/QuadGK/Distributions already root deps)"
  - "ood_verdict defaults with_pp=true — the re-enabled posterior-predictive channel (Memo §5 / T-7-04), computable on misspecified inputs via the iter1 _finite_or/_theta_tuple finite-guard"
  - "OOD scope: read-surface channels + ROC helpers + ood_verdict promoted here; the misspec-family generators + ood_roc_over_grid ship-gate experiment (need the promoted simulator + ImageFiltering) deferred to 07-03+"
  - "declared LinearAlgebra as a direct stdlib dep (Rule 3); Manifest external pins unchanged (NeuralEstimators 0.2.1 / Flux 0.16.10 intact) — a stdlib has no version to resolve"
metrics:
  tasks_completed: 3
  tasks_total: 3
  files_created: 4
  files_modified: 4
  completed_date: 2026-07-03
---

# Phase 7 Plan 02: Amortized Read Surfaces (NPE / NRE / OOD) Summary

**One-liner:** Promoted the spike's three amortized READ surfaces — NPE posterior inference
(`infer.jl`), the NRE amortized log Bayes factor (`bf.jl`), and the density+noise+posterior-
predictive OOD flag (`ood.jl`) — into grid-general, CPU-default `src/` services, folding in both
Memo §5 hardening items: a non-clamped KDE Bayes-factor baseline (removes the 1e-8 floor artifact
from `max|Δ logBF|`) and the re-enabled posterior-predictive OOD channel (`with_pp=true`, kept
crash-free on misspecified inputs by the iter1 finite-guard).

## Status: COMPLETE

All three tasks complete and committed atomically (Tasks 1–2 under a TDD RED→GREEN gate, Task 3
as a single feat). `Pkg.test()` is fully green: co-resolution 4/4, D-02 10/10, LoadImages 32/32,
patch 3/3, Colocalization 29/29, amortized summary 20/20, estimator registry 7/7, **infer.jl
11/11**, **bf.jl 18/18**, **ood.jl 31/31**. `spike/` provably byte-untouched throughout.

## What Was Built

### Task 1 — NPE read surface (`src/amortized/infer.jl`) — commits 523d19e (test), 0f88c96 (feat)
- `standardize_summary(S, zt, :min)` — applies the FROZEN `zt` to the continuous rows and passes
  the binary mask rows through UNCHANGED; grid-general via `_summary_row_partition` (07-01).
- `posterior_for(est, Z; N, use_gpu=false)` — coerces a vector to a `d_in×1` single dataset,
  returns `sampleposterior` draws.
- `rho_hat` / `rho_draws` — `StatsBase.reconstruct(θzt, draws)` BEFORE reading ρ row 1 (Pitfall 5).
- `delta_rho(est, Zs, Zc, θzt; N, use_gpu=false)` — the D-03 MC difference of the two ρ clouds
  (a 4-arg METHOD extending the `results.jl` `delta_rho` accessor generic; the two coexist by
  dispatch).
- `interval_width`. Every read call defaults `use_gpu=false` (D-06 CPU-reproducible path).

### Task 2 — Amortized Bayes factor + non-clamped baseline (`src/amortized/bf.jl`) — commits df0ac62 (test), abe2189 (feat)
- `amortized_log_bf(est_bf, Z_pair, log_prior_odds; use_gpu=false)` — one NRE forward pass, the
  `[0 1]` null/coloc grid, MEASURED `log_prior_odds` subtracted (Pitfall 5).
- `pair_encode(Zs, Zc)` — the A7 difference encoding; `nc=G²` derived from `length(Zs)÷2`, so the
  output length is `ratio_input_dim(G)=5G²` for any grid.
- `kde_log_bf_unclamped` + `_p_gt_threshold_unclamped` — the `compute_BayesFactor` KDE odds-ratio
  math WITHOUT the spike baseline's `_clampp=1e-8` floor, so `max|Δ logBF|` reflects the estimator
  rather than the clamp (Memo §5 / T-7-06). Uses KernelDensity/QuadGK/Distributions — all already
  root `[deps]`, so no re-resolve of the Wave-0 Manifest.
- Did NOT attempt an l-POP loss (the v0.2.1 ratio net loss is hardcoded `logitbinarycrossentropy`;
  a custom loss is silently ignored — 07-PATTERNS delta).

### Task 3 — OOD flag with re-enabled PP channel (`src/amortized/ood.jl`) — commit 8794388
- **Density channel:** `fit_ood_nulls` (continuous rows only via `_summary_row_partition`, `1e-6·I`
  ridge, a hand-rolled column covariance so no `Statistics.cov` import is needed) + `maha_score`.
- **Posterior-predictive channel (Memo §5 re-enabled):** `_finite_or`/`_theta_tuple` finite-guard
  (iter1) + `pp_mismatch_score`; `ood_verdict` defaults `with_pp=true` (the spike ran
  `with_pp=false`). The finite-guard maps a non-finite θ̂ to an in-range, prior-valid tuple so PP
  re-simulation stays computable on strong misspecification (T-7-04).
- **Image-noise channel:** `noise_features` (10 scale/rotation/permutation-invariant features via a
  finite-difference Laplacian, no FFT/new dep), `fit_noise_null`, `noise_score`.
- **ROC/operating points:** hand-rolled `roc_auc` (tie-aware Mann–Whitney U), `id_threshold`
  (pre-registered ID quantile — the gate), `youden_j` (post-hoc reference only).
- **Fusion:** `ood_flag` (spike OR density∨pp) and the productionized `ood_verdict(ood_nulls, Zs,
  Zc)` — OR-fuses every available channel via robust-z against the train reference and returns an
  `OODVerdict` (results.jl).

## Deviations from Plan

### [Rule 3 - Blocking] Declared `LinearAlgebra` as a direct stdlib dependency
- **Found during:** Task 3 (module load).
- **Issue:** `ood.jl` uses `cholesky`/`Symmetric`/`I` (the Mahalanobis channels), but
  `LinearAlgebra` was not in `Project.toml [deps]`, so `import LinearAlgebra` failed precompile.
- **Fix:** Added `LinearAlgebra` to `[deps]`. It is a stdlib already present in the Manifest (81
  references) with no version to resolve, so the fix is a no-op for external resolution: the
  Manifest diff is only the ProteinCoLoc dep-list line + `project_hash`; NeuralEstimators 0.2.1 and
  Flux 0.16.10 (the co-resolution gate pins) are unchanged and the gate stays 4/4 green. Chosen
  over relying on transitive loading of `\`/`cholesky` (an anti-pattern) or hand-rolling a matrix
  solve.
- **Files modified:** Project.toml, Manifest.toml.
- **Commit:** 8794388.

### [Scope adaptation] OOD ship-gate experiment left to a later plan
- **Found during:** Task 3.
- **Issue:** The spike `ood.jl` also carries the misspecification-grid ROC experiment
  (`misspec_texture`/`misspec_noise`/`misspec_optics`/`misspec_background`, the negative controls,
  and `ood_roc_over_grid`). Those need the promoted simulator (`simulate_pair`, not yet in `src/`)
  and `ImageFiltering` (`imfilter`/`Kernel.gaussian`, not a core dep), and are the per-grid D-05
  ship-gate's job — not the amortized READ surface this plan promotes.
- **Fix:** Promoted the read-surface channels + ROC helpers + `ood_verdict` only; the misspec
  families + `ood_roc_over_grid` are deferred to the ship-gate plan (07-03+). This matches the
  plan's artifact focus ("Mahalanobis + noise + posterior-predictive channels, OR-fusion,
  ood_verdict") and keeps `ImageFiltering` out of the core resolve. No shipped capability lost.

### [Design] TDD RED confirmed statically, not via a failing full-suite run
- **Found during:** Tasks 1–2.
- **Issue:** `Pkg.test()` re-precompiles the heavy Flux/NeuralEstimators/GLMakie stack; running the
  full suite once per RED phase (in addition to GREEN) is costly.
- **Fix:** Committed the RED test first (the referenced symbols are provably undefined until the
  GREEN feat lands), then the GREEN feat, then validated the whole suite green. The RED→GREEN
  commit order is preserved in history (523d19e→0f88c96, df0ac62→abe2189); only the intermediate
  full-suite failing run was skipped for cost.

## Authentication Gates
None.

## Known Stubs
- `pp_mismatch_score` / `frozen_summary` reference the Phase-2 simulator (`simulate_pair`,
  `build_mci`) which a later Phase-7 plan promotes into `src/` (the same established 07-01 seam).
  They are referenced ONLY inside function bodies, so the module LOADS today; the PP channel
  becomes end-to-end runnable once the simulator is promoted. The finite-guard (`_theta_tuple`)
  and every other OOD function are fully exercised now. This is an intended on-ramp seam, not an
  accidental gap.

## Threat Flags
None new. The plan's `<threat_model>` mitigations are implemented: T-7-04 (PP robustness) via the
`_finite_or`/`_theta_tuple` finite-guard behind the re-enabled `with_pp=true` path; T-7-06 (BF
metric integrity) via `kde_log_bf_unclamped` removing the 1e-8 floor.

## Deferred Items
- The per-grid D-05 ship-gate experiment (misspec families, `ood_roc_over_grid`, the SBC/BF
  harness) — 07-03+.
- End-to-end PP-channel exercise awaits the simulator promotion (a later Phase-7 plan).

## Self-Check

Created files:
- FOUND: src/amortized/infer.jl
- FOUND: src/amortized/bf.jl
- FOUND: src/amortized/ood.jl
- FOUND: .planning/phases/07-productionization-conditional-on-go/07-02-SUMMARY.md

Commits:
- FOUND: 523d19e (Task 1 test/RED)
- FOUND: 0f88c96 (Task 1 feat/GREEN)
- FOUND: df0ac62 (Task 2 test/RED)
- FOUND: abe2189 (Task 2 feat/GREEN)
- FOUND: 8794388 (Task 3 feat)

Verifications:
- `Pkg.test()`: fully green (infer 11/11, bf 18/18, ood 31/31; all prior suites unchanged;
  co-resolution gate 4/4 — NeuralEstimators 0.2.1 / Flux 0.16.10 pins intact after the
  LinearAlgebra stdlib declaration).
- spike/Project.toml + spike/Manifest.toml + all spike/ sources: provably UNTOUCHED
  (`git diff --quiet spike/`).

## Self-Check: PASSED
