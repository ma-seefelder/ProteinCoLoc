# Phase 7: Productionization - Discussion Log

> **Audit trail only.** Not consumed by downstream agents. Decisions live in 07-CONTEXT.md.

**Date:** 2026-07-03
**Phase:** 07-productionization-conditional-on-go
**Mode:** discuss (interactive)
**Areas discussed:** API framing (reframe) · Type hierarchy · Dependency/packaging · Result types · Estimator registry · Ship-gate

## Reframe (user-initiated, via "Other" on area selection)

User: "You do not have to care for compatibility with the old code as this will be a new
breaking release; so changes to the API are completely acceptable."

Clarified in plain text → user chose **Option 3**: amortized-only public API; Turing/ADVI kept
internal (ship-gate + future validation); plus directive: **API must be extensible and
maintainable — clear types, understandable type hierarchy.**

This relaxes ROADMAP SC1 / PROD-01 backward-compat clause. Captured as D-01/D-02.

## Questions & answers

| Area | Options offered | User's answer |
|------|-----------------|---------------|
| Turing disposition | replace / coexist-redesigned / **amortized-only, Turing internal** | **Option 3** (D-01) |
| Packaging | **hard deps** / package extension / split subpackage | **Hard deps in Project.toml** (D-03) |
| Result types | **abstract supertype + concrete** / single rich struct / parametric | **Abstract supertype + concrete + accessor interface** (D-02) |
| Estimator registry | 8×8 + train API / 8×8 only / **small grid family** | **5-grid family (4×4,8×8,16×16,32×32,64×64)** + train_and_register (D-04) |
| Ship-gate | gate 8×8 (others provisional) / **gate all 5** / gate 8×8 block others | **Gate all 5 grids before ship** (D-05) |

## Consequences flagged to user

- D-04 + D-05 make Phase 7 large: 5 grids × (data-gen → NPE+NRE train → fresh per-grid
  ship-gate). Likely needs splitting into per-grid plans. Captured as a scope flag.
- 64×64 (4096 patches) is a feasibility risk — tiny noisy patches; calibration may not
  transfer from 8×8. Captured as a researcher risk flag.

## Deferred ideas captured

- Phase-7 splitting into per-grid plans (planner sizing).
- 64×64 feasibility assessment before committing to gate it (researcher).
- Re-enabling OOD PP channel in production.
- RxInfer cross-check (post-Go, paper nice-to-have, BACK-01).
