# Phase 9: Cross-Method Comparator Harness - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-07-02
**Phase:** 09-cross-method-comparator-harness
**Mode:** --auto (recommended defaults auto-selected; no interactive questions)
**Areas discussed:** Shared input source & format, Classical estimators (reuse vs new), Tapqir bridge scope, Comparison table schema & positioning, Reproducibility & harness structure

---

## Shared input source & format

| Option | Description | Selected |
|--------|-------------|----------|
| Simulator-generated MultiChannelImages over a seeded θ grid | Reuse `simulate_pair`; ground-truth regime per row; identical container for classics + NPE | ✓ |
| External corpus / CBS as primary input | Would take a Phase-8 dependency; Phase 9 is parallelizable-now | |
| Bespoke synthetic arrays | Bypasses the shared `MultiChannelImage` contract, breaks input parity | |

**Auto choice:** Simulator-generated `MultiChannelImage` set, input-source-agnostic API so the Phase-8 corpus can be fed later without redesign (D-01/D-02/D-03).
**Notes:** "Shared inputs" (SC1) means all estimators + NPE see byte-identical `MultiChannelImage`s; simulator θ supplies the ground-truth regime label needed to show *when* a classic is wrong.

---

## Classical estimators — reuse vs new

| Option | Description | Selected |
|--------|-------------|----------|
| Reuse existing Manders/Pearson/Spearman, implement only Costes-p | Factor `encode_aug` Manders + `src` `correlation`; add missing Costes randomization | ✓ |
| Reimplement all four fresh in the harness | Risks diverging from ablation summary math; wasteful | |
| Pull estimators from an external package | New dependency; Windows/decoupling friction | |

**Auto choice:** Reuse `spike/data/encode.jl:encode_aug` (Manders M1/M2, whole-image Pearson) and `src/colocalization.jl:correlation` (Pearson/Spearman); implement Costes-p new in a spike-local module (D-04/D-05/D-06).
**Notes:** Costes randomization p-value is missing everywhere and must be built. All new code stays in `spike/comparator/` (no `src/` edits).

---

## Tapqir bridge scope

| Option | Description | Selected |
|--------|-------------|----------|
| Minimal published-example sanity anchor, optional/skippable | Reproduce a published Tapqir example; PythonCall/CondaPkg; graceful skip if env absent | ✓ |
| Full Tapqir comparator column on our simulator images | Invalid — Tapqir models CoSMoS single-molecule spots, not diffuse 2-channel coloc | |
| Drop Tapqir entirely | Violates SC2 | |

**Auto choice:** Minimal PythonCall/CondaPkg bridge reproducing a published Tapqir example as a sanity anchor, degrading gracefully to skip-with-flag (D-07/D-08).
**Notes:** SC2 only asks for a published-example anchor; Tapqir is an anchor, not a comparator column, because of the data-regime mismatch. Not a hard gate on the classical table.

---

## Comparison table schema & "knows when classics are wrong"

| Option | Description | Selected |
|--------|-------------|----------|
| Tidy per-input DataFrame + divergence flag + traffic-light | Rows=inputs, cols=methods + NPE + OOD; derived divergence column; BayesInteractomics audit style | ✓ |
| Bare per-method scalar table | Meets SC1 letter but misses the "when classics are wrong" positioning | |
| Free-form report only | Not reproducible/queryable | |

**Auto choice:** Tidy DataFrame with ground-truth regime + all methods + NPE ρ̂/OOD + a divergence indicator, persisted CSV+JLD2, presented with the traffic-light pattern (D-09/D-10/D-11).
**Notes:** Reuses BayesInteractomics `_bin_calibration`/`CalibrationResult`/traffic-light and `model_diagnostics` reporting per SC3.

---

## Reproducibility & harness structure

| Option | Description | Selected |
|--------|-------------|----------|
| Single seeded entry point + spike test + content-addressed artifact | `run_comparator.jl`, Random123/Philox, determinism + finiteness + Tapqir-anchor test | ✓ |
| Ad-hoc script, no test | Not reproducible; fails SC3 | |

**Auto choice:** `spike/comparator/run_comparator.jl` seeded with Philox; artifact content-addressed like Phase-3 cache; spike test asserts determinism, finite classical outputs, and Tapqir-anchor-or-clean-skip (D-12/D-13). Pre-declared divergence/traffic-light thresholds to avoid data-snooping (D-14).
**Notes:** Carries the roadmap's Phase-5 anti-data-snooping blocker into this phase's threshold discipline.

---

## Claude's Discretion

- Costes randomization scheme (pixel-scramble vs block/Van Steensel shift) and scramble count.
- Table column ordering, optional Arrow output, figure rendering.
- Module/file layout under `spike/comparator/`.

## Deferred Ideas

- Classics + NPE on the Phase-8 physical corpus / CBS for the final head-to-head → Phase 16.
- Costes regression-based auto-threshold beyond Otsu.
- Full Tapqir integration as a live comparator column (rejected as invalid for our data regime).
- Manuscript prose framing of the delta → Phase 10.
