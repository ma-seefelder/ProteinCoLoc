# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-06-26)

**Core value:** A trained NPE/NRE produces calibrated, amortized colocalization inference (posterior + Bayes factor) in a single forward pass, >100x faster than per-dataset ADVI, with a demonstrated SBC/coverage proof and an honest OOD flag.
**Current focus:** Phase 1 — Environment + Smoke Gate

## Current Position

Phase: 1 of 7 (Environment + Smoke Gate)
Plan: 0 of TBD in current phase
Status: Ready to plan
Last activity: 2026-06-26 — Roadmap created (7 phases, 27/27 v1 requirements mapped)

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity:**
- Total plans completed: 0
- Average duration: - min
- Total execution time: 0.0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| - | - | - | - |

**Recent Trend:**
- Last 5 plans: -
- Trend: -

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 7-layer horizontal build order; ENV smoke gate (Phase 1) is the hard gate before all downstream investment (pre-1.0 NeuralEstimators API is the top unknown)
- [Roadmap]: SBC/BF/OOD grouped into one Phase 5 — sibling plans sharing one θ*~π→simulate→infer harness over the trained nets
- [Roadmap]: Phase 7 (productionization) is conditional on a Go decision in the Phase 6 memo and is the ONLY phase that edits `src/`; spike phases (1–6) keep `src/` provably untouched

### Pending Todos

[From .planning/todos/pending/ — ideas captured during sessions]

None yet.

### Blockers/Concerns

[Issues that affect future work]

- [Phase 1]: NeuralEstimators v0.2.x is pre-1.0 and the Julia ML ecosystem is mid Flux→Lux/Reactant migration — the CPU-only Flux path must be confirmed exercised at the smoke gate; pinned Manifest is the mitigation
- [Phase 5]: OOD detection power is structurally bounded by the fixed patch-correlation summary — summary-orthogonal misspecifications are provably undetectable and must be named, not hidden
- [Phase 5]: "Tune until calibrated" is a data-snooping hazard — M and SBC threshold must be pre-registered before Phase 5 planning runs

## Deferred Items

Items acknowledged and carried forward from previous milestone close:

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| Backend (v2) | BACK-01 Turing→RxInfer.jl feasibility evaluation | Deferred (post-Go, never a spike dependency) | 2026-06-26 |
| Summary net (v2) | BACK-02 DeepSet permutation-invariant summary upgrade | Deferred (upgrade if MLP summary proves insufficient) | 2026-06-26 |

## Session Continuity

Last session: 2026-06-26
Stopped at: ROADMAP.md and STATE.md written; REQUIREMENTS.md traceability updated
Resume file: None
