---
phase: 10-manuscript-skeleton-and-related-work-positioning
plan: 02
subsystem: manuscript
tags: [typst, claim-table, honesty-posture, manuscript-skeleton]
requires:
  - "manuscript/claims.typ placeholder + claim_table() from Plan 10-01"
  - "manuscript/main.typ importing claims, claim_table (Plan 10-01)"
provides:
  - "Full seeded claim spine: 9-row array-of-dicts of headline v2.0 claims (C1-C9)"
  - "claim_table() render fn extended with a Notes/honest-caveat column"
  - "Honesty slots by construction: speedup caveat, calibration-under-simulator, OOD blind spot"
affects:
  - "manuscript/main.typ <tab:claims> figure (renders the seeded rows)"
tech-stack:
  added: []
  patterns:
    - "RESEARCH Pattern 1: data-driven claim table (array-of-dicts + render fn, one source)"
    - "Pitfall 4 avoidance: claims.typ stays #let-only so #import emits nothing"
key-files:
  created: []
  modified:
    - "manuscript/claims.typ"
decisions:
  - "Added a notes field + Notes column (rather than extra rows) to carry verbatim honest caveats, so honesty is visible in the rendered PDF and diffable in source"
  - "Used fr column widths (auto,2fr,auto,auto,auto,auto,3fr) so long claim/notes text wraps instead of overflowing the A4 page"
metrics:
  duration: ~6 min
  completed: 2026-07-02
  tasks: 1
  files: 1
---

# Phase 10 Plan 02: Claim Spine Seeding Summary

Seeded `manuscript/claims.typ` with the full 9-row honestly-caveated claim spine (C1–C9), each row mapping a v2.0 headline claim to its differentiator axis, owning/supporting phase, figure spec ID, status, and a verbatim honest-caveat note — replacing the single Plan 10-01 placeholder while keeping the file `#let`-only and the manuscript build green.

## What Was Built

Task 1 replaced the placeholder `C0` row in `#let claims` with the complete headline claim set derived from CONTEXT D-06 and the RESEARCH §"Phase Requirements" mapping table, and extended `claim_table()` to render a seventh **Notes / honest caveat** column:

| ID | Claim | Axis | Phase | Fig | Status |
|----|-------|------|-------|-----|--------|
| C1 | Amortized ms-scale inference (single forward pass) | Amortization | 4,7 | F1 | supported (measured) |
| C2 | >100x speedup at comparable RMSE (caveated) | Speed | 4 | F1 | supported (caveat) |
| C3 | Formal SBC rank-uniformity + coverage proof | Calibration/SBC | 5 | F2 | reserved slot |
| C4 | Honest OOD / misspecification flag | OOD | 5 | F3 | reserved slot; blind-spot named |
| C5 | Amortized Bayes factor (2-way → 3-way) | Model comparison | 5,13 | F7 | reserved slot |
| C6 | Registration uncertainty as inferred latent | Registration-UQ | 11 | F5 | reserved slot |
| C7 | Amortized spatial per-region map + per-region UQ | Spatial map | 12 | F6 | reserved slot |
| C8 | Cross-method: knows when the classics are wrong | Comparator | 9 | F4 | reserved slot |
| C9 | Non-circular external validation | External validation | 8 | F4 | reserved slot |

**Honesty posture baked in by construction:**
- **C2 notes** carry the verbatim Phase-4 caveat: median ~325× forward-pass (min pair ~125×), ~16× on the full posterior-sample (N=2000) workload, realized ADVI baseline ~0.5 s/pair, paired with RMSE parity (NPE 0.1271 vs ADVI 0.1270) and wider NPE 90% intervals (0.884 vs 0.605). No bare ">100× vs minutes".
- **C3 notes** state calibration is "under the simulator" and always paired with the OOD result.
- **C4 notes** name the summary-orthogonal OOD blind spot ("provably undetectable … named, not hidden").

## Acceptance Criteria

| Criterion | Result |
|-----------|--------|
| `#let claims` array ≥ 9 rows | 9 rows (C1–C9) |
| Speedup row text contains "325", "16", "0.5" | present (`325`, `16x`, `0.5 s/pair`) |
| ≥1 status "supported" and ≥1 "reserved slot" | 4 "supported", 8 "reserved slot" |
| Row/notes reference "under the simulator" + OOD blind spot | "under the simulator" ×4; "blind spot" ×2 |
| Every row has non-empty `fig` and `phase` (D-07) | all 9 rows populated |
| `grep -c "reserved slot"` ≥ 5 | 8 |
| `bash manuscript/build.sh` exits 0 | exit 0, `manuscript/main.pdf` rebuilt |

## Deviations from Plan

None — plan executed exactly as written. The plan's `<action>` explicitly permitted "a notes field" for the honesty caveats; a Notes column was added to `claim_table()` so those caveats render in the PDF and remain diffable in source. This is within the stated schema latitude, not a deviation.

## Threat Surface

No new surface. Static Typst source only; no runtime, network, input, or data flow. The T-10-02 disposition (mitigate overclaiming) is satisfied by construction: the status column plus verbatim Phase-4 caveats prevent asserting pending results as done.

## Self-Check: PASSED

- FOUND: `manuscript/claims.typ` (modified, 9-row spine, `claim_table` present)
- FOUND commit `8b24252` (feat(10-02): seed full honestly-caveated claim spine)
- `bash manuscript/build.sh` → exit 0
- No files modified outside `manuscript/claims.typ`; STATE.md / ROADMAP.md untouched
