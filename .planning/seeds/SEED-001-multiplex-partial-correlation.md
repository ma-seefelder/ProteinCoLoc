---
id: SEED-001
status: dormant
planted: 2026-07-02
planted_during: v2.0 — AmortizedColoc (Phase 5)
trigger_when: at the start of the v2.1 milestone (once v2.0 — AmortizedColoc is completed/archived)
scope: large
---

# SEED-001: N-channel / multiplex joint colocalization via a Bayesian Gaussian graphical model — direct vs. indirect via partial correlation

## Why This Matters

v2.0 is deliberately restricted to two channels (`_prepare_data` in `src/colocalization.jl` is hardwired to a channel pair). Multiplex imaging (CODEX/IMC/multiplexed IF) needs *joint* colocalization across ≥3 channels, and — more importantly — the ability to separate **direct** co-occurrence from **indirect** association (A–B appears colocalized only because both track C). A Bayesian Gaussian graphical model with shrinkage over the N-channel correlation matrix yields **partial correlations**, which is exactly that direct-vs-indirect distinction with calibrated uncertainty.

This is the imaging analog of the spoke-vs-matrix / direct-contact problem from the interactome work — a methodological through-line across the portfolio, and the highest-ceiling extension because of the multiplex-imaging audience.

## When to Surface

**Trigger:** at the start of the v2.1 milestone, i.e. once the v2.0 — AmortizedColoc milestone is completed/archived. Surface during `/gsd:new-milestone` when the next milestone's scope is defined.

Deferred from v2.0 (not dropped) to keep the v2.0 tool paper focused and avoid overloading a solo release.

## Scope Estimate

**Large** — a full milestone. New generative structure (N-channel simulator), a graphical-model prior with shrinkage, retrained/extended estimators, and its own validation. Needs explicit positioning against the incumbent.

## Positioning / Risk

- **Incumbent to beat:** Jaqaman lab conditional colocalization analysis in 3-color microscopy (JCB 2022, e202106129). The space is neither empty nor trivially ours.
- **Open angle:** calibrated UQ + direct-vs-indirect via partial correlation (Gaussian graphical model), which the incumbent does not provide.

## Breadcrumbs

- `next_project_analysis/11_plan_proteincoloc_v2.md` — Scope-Guard (multiplex explicitly deferred to v2.1)
- `next_project_analysis/decisions/proteincoloc_v2_scope_reasoning_log.md` — functionality POSTERIOR: multiplex = v2.1 DEFER; SPIKE 3 (incumbent found)
- `src/colocalization.jl` — `_prepare_data` (2-channel hardwire to generalize), `correlation()` (grid the graphical model extends)
- Related: SEED-002 (copula / nonlinear dependence — the other v2.1 modeling deferral)

## Notes

The other half of the v2.0→v2.1 deferral. Consider whether partial correlation (linear, Gaussian graphical model) and copula dependence (SEED-002) should be one combined v2.1 "dependence-structure" paper or two.
