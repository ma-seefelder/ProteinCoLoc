---
phase: 07-productionization-conditional-on-go
plan: 07
subsystem: testing
tags: [neural-estimators, flux, sbc, bayes-factor, ood, ship-gate, npe, nre, amortized-inference, reproducibility, fine-grid, min-image-size]

# Dependency graph
requires:
  - phase: 07-03
    provides: _train_grid_pipeline + default_imsize_for(16) (≥512²-biased) + default_npairs_for(16)
  - phase: 07-04
    provides: per-grid CPU ship-gate machinery (run_gate.jl, harness.jl, sbc.jl, gate_consts_template.jl)
  - phase: 07-05
    provides: proven 8×8 pipeline shape + reference gate-8x8.md verdict to interpret against
provides:
  - "16×16 NPE+NRE+OOD-nulls bundle produced by _train_grid_pipeline(16), CPU-resident (artifacts/grid_16/*.jld2, on-disk, gitignored)"
  - "Fresh 16×16 ship-gate pre-registration test/gate/gate_consts_16.jl (PROD_SEED[16]=0xb906f369f6cacf91, SBC_IMSIZE=(512,512))"
  - "Recorded 16×16 CPU SBC/BF/OOD ship-gate outcome (gate-16x16.md) — honest FAIL + unscored OOD"
  - "Documented ≥512² minimum-image-size caveat the registry entry must carry"
affects: [07-10, registry-population]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Fine-grid gate raises SBC_IMSIZE per grid so patches clear the ≥15-survivor floor (16×16 → 512²)"
    - "Ship-gate reads the grid's OWN locked gate_consts_G.jl + fresh disjoint PROD_SEED[G] (anti-snooping)"

key-files:
  created:
    - test/gate/gate_consts_16.jl
    - .planning/phases/07-productionization-conditional-on-go/gate-16x16.md
    - artifacts/grid_16/npe_16.jld2 (on-disk, gitignored)
    - artifacts/grid_16/ratio_16.jld2 (on-disk, gitignored)
    - artifacts/grid_16/ood_nulls_16.jld2 (on-disk, gitignored)
    - artifacts/grid_16/gate_report_16.jld2 (on-disk, gitignored)
  modified: []

key-decisions:
  - "Recorded the honest FAIL on the byte-locked gate_consts_16.jl — SBC ks_pass=false (4/8 KS rejections INCLUDING the headline ρ_true), BF corr 0.9150 < 0.95 AND max|Δ logBF| 12.72 > 0.5. No constant was re-tuned, relaxed, or reinterpreted."
  - "Recorded the OOD gate as NOT RUN / INCONCLUSIVE (auc=nothing, passed=nothing) — no pos_sim positive control was injected, so no separability AUC exists. Explicitly NOT presented as a pass; logged as a gap that must be closed before the 16×16 registry entry can claim a validated misspecification flag."
  - "Flagged 16×16 as materially weaker than 4×4/8×8 rather than folding it into the memo's 'every residual reads to a non-method cause' framing: the ρ_true KS rejection (7.0e-4) and the lowest-of-three BF corr (0.9150, on the FINITE pairs) are grid-resolution effects, not M=2000 test-power or clamped-KDE-tail artifacts."
  - "Did NOT force-commit artifacts/grid_16/*.jld2 — .gitignore:63 ignores artifacts/; the recorded markdown numbers are the deliverable (07-05/07-06 precedent)."

patterns-established:
  - "When a gate's residuals are NOT attributable to pre-registration/baseline design, say so explicitly rather than reusing the prior grids' exculpatory reading"

requirements-completed: [PROD-02]

# Metrics
duration: 16.5min (gate run) + training (earlier session)
completed: 2026-07-21
---

# Phase 7 Plan 07: 16×16 Per-Grid Ship-Gate Summary

**16×16 (SHIP-WITH-CAVEAT fine grid; summary_dim 512, 256 patches) trained via
`_train_grid_pipeline(16)` on 80 000 pairs and gated CPU-reproducibly on the fresh disjoint
`PROD_SEED[16] = 0xb906f369f6cacf91` against the byte-locked `gate_consts_16.jl` (SBC_IMSIZE raised
to 512²). Verdict: honest FAIL — SBC `ks_pass=false` (4/8 KS rejections, including the headline
ρ_true at p = 7.0e-4), BF corr 0.9150 and max|Δ logBF| 12.72 both failing, and the OOD gate
UNSCORED (`auc = nothing`, `passed = nothing`) because no positive control was injected. This is the
weakest of the three gated grids.**

## Performance

- **Duration:** ~16.5 min gate run (training completed in an earlier session)
- **Completed:** 2026-07-21
- **Tasks:** 3 (Task 1 training, Task 2 pre-registration, Task 3 gate run + record)
- **Files committed:** `test/gate/gate_consts_16.jl` (7e2318b, before the run), `gate-16x16.md`

## Accomplishments
- Produced the 16×16 bundle via `_train_grid_pipeline(16)` (n_pairs 80 000 = `default_npairs_for(16)`,
  68 000 train / 12 000 val, `use_gpu=false`, NPE `d_in=512` / `dstar=64` / 10 coupling layers, NRE
  `input_dim=1280` = `ratio_input_dim(16)`), persisted CPU-resident; all three artifacts load CPU-only.
- Committed the fresh 16×16 pre-registration `gate_consts_16.jl` **before** the reported run
  (anti-snooping, T-7-08), with `SBC_IMSIZE` raised to (512,512) so each 32×32 patch carries 1024 px.
- Ran the determinism-pinned (single-thread, `use_gpu=false`) SBC/BF/OOD ship-gate keyed by
  `PROD_SEED[16]`; `gate_report_16.jld2` written with `status = :ran`.
- Recorded the verdict verbatim in `gate-16x16.md`, including the ≥512² minimum-image-size caveat
  the registry entry must carry, and the explicit OOD NOT-RUN record.

## Task Commits

1. **Task 2: Fresh 16×16 pre-registration (gate_consts_16.jl)** — `7e2318b` (feat) — committed before the gate run.
2. **Task 1: Produce the 16×16 bundle** — no commit (outputs are `artifacts/grid_16/*.jld2`, gitignored; on-disk + load-verified).
3. **Task 3: Run + record the 16×16 CPU SBC/BF/OOD gate + caveat** — this plan's commit (gate-16x16.md; `gate_report_16.jld2` gitignored).

## Gate Verdict (honest, un-tuned)

| Gate | Verdict | Headline numbers |
|------|---------|------------------|
| SBC (M=2000, L=999) | **FAIL** (`ks_pass=false`, `ece_pass=true`) | ECE green on all 8 (ρ_true 0.0240, Δρ 0.0060); **KS 4/8 pass** — rejections on ρ_true (7.0e-4), spillover (2.7e-3), label_efficiency (7.1e-13), shift_dx (0.0167) |
| BF (n=15 finite/25) | **FAIL** (`corr_pass=false`, `tol_pass=false`) | corr 0.9150 (< 0.95; lowest of the three grids), max\|Δ logBF\| 12.72 (> 0.5) |
| OOD (n=200) | **NOT RUN / INCONCLUSIVE** | `id_threshold = 378.44` only; `auc = nothing`, `passed = nothing` — no `pos_sim` injected, no separability evidence. **Not a pass.** |

Interpretation is recorded in `gate-16x16.md` and kept strictly separate from the verdict: ECE is
green on all 8 and `max|Δ logBF|` carries the documented clamped-KDE tail signature (Phase 5 iter2 /
memo §5), but the ρ_true KS rejection and the corr = 0.9150 mid-range gap are **not** covered by
those readings and are treated as genuine 16×16 grid-resolution weaknesses.

## Deviations from Plan

**1. [Gap] OOD gate produced no verdict — no positive control injected**
- **Found during:** Task 3 (gate run)
- **Issue:** The plan's Task 3 asked for an OOD verdict, but `run_gate.jl --grid 16 --ood` was invoked
  without a `pos_sim` positive-control simulator, so `ood_gate` computed only the pre-registered ID
  operating point (`id_threshold = 378.44`) and returned `auc = nothing`, `passed = nothing`.
- **Disposition:** Recorded as NOT RUN / INCONCLUSIVE, **not** as a pass. The same absence occurred at
  8×8 and 4×4 and is partly explained by the deliberate deferral of the controlled
  misspecification-grid ROC experiment (`misspec_*` families, `ood_roc_over_grid`) to a later Phase-7
  plan per the `src/amortized/ood.jl` scope note — but the deferral does not convert absence into
  evidence. Injecting families ad-hoc here would be un-pre-registered and would violate the phase's
  anti-snooping discipline, so no ad-hoc control was added.
- **Carried to:** the later misspec-ROC plan / 07-10; the 16×16 registry entry cannot claim a
  validated OOD flag until scored.

**2. [Provenance gap] Training `imsize_set` not recoverable from the artifacts**
- **Found during:** Task 3 (writing the report)
- **Issue:** `_train_grid_pipeline` persists `grid`, `n_pairs`, `use_gpu`, `n_train`, `ratio_n` in the
  artifact meta but **not** `imsize_set`, and the datagen pool is in-memory (no cache dir left behind).
  The realized training image distribution for the 16×16 bundle therefore cannot be re-derived from
  `artifacts/grid_16/*` alone. The pipeline default is the ≥512²-biased `default_imsize_for(16)`.
- **Disposition:** Recorded honestly as a provenance gap in `gate-16x16.md` rather than asserting a
  value that was not verified. A future plan should add `imsize_set` to the persisted meta.

**3. [Repo convention] Did not commit the trained/report .jld2 artifacts**
- `.gitignore:63` ignores `artifacts/`; the recorded markdown numbers are the committed deliverable
  (07-05 / 07-06 precedent).

---

**Total deviations:** 3 (1 unscored gate recorded as a gap, 1 provenance gap, 1 repo convention).
No pre-registration constant was altered.

## Next Phase Readiness
- The 16×16 bundle exists on-disk (`artifacts/grid_16/`, gitignored); 07-10 can regenerate via
  `train_and_register(16)` (SKIP-IF-DONE loads if present).
- **07-10 caution:** unlike 4×4 and 8×8, the 16×16 residuals are not fully attributable to
  pre-registration/baseline design. If registered at all it must be SHIP-WITH-CAVEAT and carry
  (a) the ≥512² minimum-image-size requirement, (b) the ρ_true SBC-KS rejection, (c) the
  un-validated OOD flag. The disposition decision belongs to 07-10 / the phase verifier.

## Known Stubs
None. The trained bundle is real (finite CPU draws), the gate report is real (`status = :ran`), and
the recorded numbers are the honest gate output. The OOD `auc = nothing` is an **unscored gate**, not
a stub and not a pass.

## Self-Check: PASSED

- On-disk artifacts: `npe_16.jld2`, `ratio_16.jld2`, `ood_nulls_16.jld2`, `gate_report_16.jld2` — all FOUND (gitignored).
- Committed deliverables: `test/gate/gate_consts_16.jl` TRACKED (7e2318b); `.planning/.../gate-16x16.md` TRACKED.
- `spike/` byte-untouched; `src/` untouched by this plan.

---
*Phase: 07-productionization-conditional-on-go*
*Completed: 2026-07-21*
