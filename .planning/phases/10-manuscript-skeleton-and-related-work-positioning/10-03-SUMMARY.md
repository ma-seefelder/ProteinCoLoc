---
phase: 10-manuscript-skeleton-and-related-work-positioning
plan: 03
subsystem: manuscript
tags: [related-work, positioning, tapqir, typst, bibliography]
requires:
  - "manuscript skeleton (10-01): main.typ, refs.bib, related_work.typ stub, tiered-bib wiring"
provides:
  - "Populated related-work section: seven-axis capability matrix + four-axis Tapqir differentiation"
  - "Verified Tapqir comparison cells (A1 SBC / A2 BF / A3 registration) against Ordabayev et al. 2022"
affects:
  - "manuscript H-tier reference list (first real citations now render)"
tech-stack:
  added: []
  patterns:
    - "Wide table figure spanning both columns via place(top, scope: parent, float: true) + kind: table"
    - "Per-cell verification provenance recorded as a // source comment block before locking prose"
key-files:
  created: []
  modified:
    - manuscript/sections/related_work.typ
decisions:
  - "All three ASSUMED Tapqir cells (A1/A2/A3) confirmed by the paper; none required correction, framed conservatively"
  - "Reused existing refs.bib keys (ordabayev2022, costes2004, manders1993, talts2018); no refs.bib edits needed"
metrics:
  duration: ~10 min
  completed: 2026-07-02
  tasks: 2
  files: 1
---

# Phase 10 Plan 03: Related-Work Positioning Summary

Drafted `manuscript/sections/related_work.typ` as a seven-axis capability comparison matrix over six colocalization methods plus a mandatory four-axis Tapqir differentiation, with the three flagged-ASSUMED Tapqir cells verified against Ordabayev et al. 2022 (eLife 2022;11:e73860, DOI 10.7554/eLife.73860) before the prose was locked; the manuscript still compiles green (SC2/D-12).

## What Was Built

- **Task 1 — Tapqir cell verification (commit f3d9898).** Fetched the eLife full text and verified the three ⚠ cells, recording a per-cell verdict + source as a `//` provenance comment block at the top of `related_work.typ`:
  - **A1 (formal SBC):** CONFIRMED — Tapqir validates via posterior-predictive checking (Gelman 2013) and simulation parameter-recovery; the article contains no simulation-based-calibration rank-uniformity or coverage test ("simulation-based" = 0 hits, "coverage" = 0 hits in the full text). Cell kept **Partial** (Bayesian per-spot probabilities, no formal SBC).
  - **A2 (explicit Bayes factor):** CONFIRMED — Tapqir emits per-spot `p(specific)` = `p(m=1)` posterior probabilities, not a dataset-level model-comparison Bayes factor ("bayes factor" = 0 hits). Cell kept **Partial**.
  - **A3 (registration-UQ):** CONFIRMED — cross-channel mapping and microscope-drift correction are external preprocessing steps (Friedman & Gelles 2015) done before inference; registration accuracy bounds localization precision but is not propagated as an inferred latent into the coloc call. Cell kept **No** (external channel mapping + drift correction). (Tapqir does model within-AOI spot xy — a distinct quantity.)
- **Task 2 — matrix + prose (commit 89582fa).** Rendered the comparison matrix as a parent-spanning `table` figure `<tab:related-work>` (rows: ProteinCoLoc v2.0, Tapqir, Costes randomization, Manders M1/M2, Pearson/Spearman, ProteinCoLoc v1 ADVI; columns: amortized inference, calibrated UQ/SBC, Bayes factor/model comparison, registration UQ as latent, spatial per-region map, OOD/misspec flag, speed). Wrote the four-axis Tapqir differentiation prose covering amortization, SBC/calibration, registration-UQ, and spatial map, plus a classical-methods lead-in.

## Honesty Posture (as required)

- v2.0 calibration stated as holding **"under the simulator"** (appears 5×), paired with the SBC framing "v2.0 adds a formal SBC coverage proof" rather than "Tapqir is uncalibrated."
- The **summary-orthogonal OOD blind spot** is explicitly named: misspecifications leaving the fixed patch-correlation summary unchanged are provably invisible to the flag.
- Speedup caveated (median ~325× forward-pass / ~16× full posterior workload vs ~0.5 s/pair ADVI at matched RMSE) — no bare ">100×".
- Tapqir vs v2.0 framed as **overlapping-but-distinct data regimes** (sparse single-molecule CoSMoS spots vs dense/diffuse patch-correlation microscopy), not blanket superiority.
- v2.0 reserved-slot cells carry phase annotations (Ph 5/11/12/13) so they read as pending, not proven.

## Verification

- `bash manuscript/build.sh` exits **0** (manuscript compiles; hard DoD met).
- Four axes present (`amortiz`, `SBC/calibrat`, `registration`, `spatial` all grep-match).
- Citations `@ordabayev2022`, `@costes2004`, `@manders1993`, `@talts2018` all resolve (no dangling citation → build would otherwise fail).
- Task-1 acceptance: DOI reference present (2 hits); `formal SBC|Bayes factor|registration` (12 hits, ≥3).
- No content-overflow or error output from the compiler.

## Deviations from Plan

None — plan executed as written. All three ASSUMED cells were confirmable against the primary source (no conservative-fallback framing needed beyond the honesty posture already required).

## Known Non-Blocking Warning

`typst compile` now prints `warning: layout did not converge within 5 attempts` (build still exits 0). Root cause verified: the template's tiered-bibliography (`manuscript/lib/tiered-bib.typ`) uses `state()` + `query(selector(ref))` (line 136) to assign each citation to its tier; this query-over-references iterates the moment the first real `@cite` enters the H tier — which this plan is the first to do. It is a template-level behavior, not caused by the wide float table (no overflow warning is emitted), and it will appear for any plan that adds citations. Not a compile failure; no degradation applied. Flagged for potential template hardening in a later phase if desired.

## Bibliography Note

No `refs.bib` edits were required: all four cited keys were already seeded in Wave 1 (`ordabayev2022`, `costes2004`, `manders1993`, `talts2018`). No new `@citekey` needed to be introduced.

## Self-Check

- File exists: `manuscript/sections/related_work.typ` — FOUND
- Commit f3d9898 (Task 1) — recorded
- Commit 89582fa (Task 2) — recorded

## Self-Check: PASSED
