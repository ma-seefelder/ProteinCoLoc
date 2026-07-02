---
id: SEED-002
status: dormant
planted: 2026-07-02
planted_during: v2.0 — AmortizedColoc (Phase 5)
trigger_when: at the start of the v2.1 milestone (once v2.0 — AmortizedColoc is completed/archived)
scope: medium
---

# SEED-002: Copula / rank-based nonlinear dependence layer for colocalization

## Why This Matters

v2.0's colocalization signal rests on Pearson (and optionally Spearman/Kendall) correlation via `correlation()` in `src/colocalization.jl`. Pearson captures only linear monotone comovement and is sensitive to the zero-exclusion threshold (`length(a) <= 15 → missing`). A **copula-based dependence layer** would capture nonlinear, saturating, or tail-dependent association and report the *form* of dependence, not just its magnitude — and there is existing copula machinery in BayesInteractomics to reuse.

## When to Surface

**Trigger:** at the start of the v2.1 milestone, once v2.0 — AmortizedColoc is completed/archived. Surface during `/gsd:new-milestone`.

Deferred from v2.0 (not dropped): it is a genuine enhancement but **not load-bearing** for the v2.0 claim, and switching off Pearson would break comparability to the published v1 tool. Keep v2.0 on the Pearson/Spearman summary.

## Scope Estimate

**Medium** — a phase or two. Add a copula dependence estimator, integrate with the amortized/summary path, validate that it recovers nonlinear dependence the Pearson summary misses. Lower risk than SEED-001.

## Breadcrumbs

- `next_project_analysis/11_plan_proteincoloc_v2.md` — Scope-Guard (copula deferred to v2.1)
- `next_project_analysis/decisions/proteincoloc_v2_scope_reasoning_log.md` — functionality POSTERIOR: copula = v2.1 DEFER (not load-bearing; breaks Pearson-comparability to v1)
- `src/colocalization.jl` — `correlation()` (Pearson/Spearman/Kendall; the `<=15` zero-exclusion threshold) is what a copula layer would supplement
- BayesInteractomics — existing copula/EM machinery to reuse
- Related: SEED-001 (multiplex partial correlation — the other v2.1 modeling deferral)

## Notes

Consider bundling with SEED-001 into a single v2.1 "richer dependence structure" paper (partial correlation + copula), or keeping copula as a lighter standalone enhancement.
