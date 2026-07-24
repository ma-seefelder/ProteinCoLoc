---
name: spike-findings-proteincoloc
description: >-
  Persisted findings from 13 ProteinCoLoc spikes, packaged as build blueprints. Use when working on
  (a) the v2.0 AmortizedColoc Bayes-factor ship-gate / KDE-vs-NRE baseline, (b) SBC calibration of
  the NPE (atoms, nuisance drift, over-power, capacity/truncation levers), or (c) the analysis-only
  src/ refactoring roadmap (hot-path perf, type hierarchy, bayes @model modularity, headless viz,
  provenance). Read the relevant references/ blueprint before implementing any of these.
---

# Spike Findings — ProteinCoLoc

Thirteen spikes across two decoupled efforts, packaged so a future build conversation can pick up the
verified conclusions without re-deriving them. **Read the matching `references/` blueprint before
building; open the `sources/NNN-*/` READMEs + scripts for full evidence.**

Two independent tracks:
- **Spikes 001–005** — a *separate, earlier* session: an **analysis-only** refactoring roadmap for
  the v1.0 `src/` core. Proposals only; nothing implemented; gated by the "don't touch `src/` during
  the spike" rule.
- **Spikes 006–013** — v2.0 AmortizedColoc **ship-gate diagnostics**: why the amended grid-8
  confirmatory gate's BF and SBC arms failed, and whether any lever fixes them.

## Requirements (project-level non-negotiables)

From `.planning/spikes/MANIFEST.md` and the v2.0 CLAUDE.md, plus the gate constraints these spikes
established:

- **Decoupling.** The v2.0 spike lives entirely in `spike/` with its own `Project.toml`/`Manifest`.
  `src/`, the main package, and both manuscript pipelines stay **provably untouched** during the
  spike. The `src/` refactors (001–005) are a **separate effort** — sequence them before resuming
  the spike or at the productionization boundary; never interleave.
- **Simulator prior ↔ Turing `@model`.** π(θ) must stay consistent with the existing Turing μ/ν/σ/τ
  prior ranges for ADVI comparability. (Truncating the μ-prior — spike 013 — *breaks* this; it is a
  known, separately-decided cost, not a free fix.)
- **Reproducibility & seed discipline.** Everything reproducible from a fixed (Random123) seed. Any
  spike touching the ship-gate runs on a **DEV seed asserted disjoint** from `PROD_SEED_V2`,
  `VAL_MASTER_SEED` (0x5BC0FFEE), `NPE_MASTER_SEED` (0xC0FFEE) and prior dev seeds. Never consume a
  pre-registered seed in a spike; never write a gate report from one.
- **Gate changes are pre-registration matters, not spike decisions.** Every BF/SBC remedy below is a
  *proposal*; **none were applied.** They bundle into one final amendment on a retrained model, per
  `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md` §6 (and §6.4 blocks
  iterating the current amended_v2 gate).
- **BF-gate reality (spikes 006–008).** BF-arm attrition is 100% baseline-side (KDE `p_post` hits
  {0,1}→±Inf); the NRE never goes non-finite. The KDE Bayes factor is **not a valid reference past
  |log-BF| ≈ log(L) ≈ 6.9** (L=999). Do not gate magnitude agreement on a tail the reference cannot
  resolve.
- **SBC reality (spikes 009–013).** The coloc **targets ρ_true/Δρ are calibrated.** ρ_true's SBC
  failure is a **prior-atom artifact** (fix test-side with randomized ranks). Residual nuisance
  failures are a ~0.06-SD marginal location drift on **non-identified** nuisances at the M=2000
  over-power edge — no realistic model lever (capacity, truncation) reliably removes it; the correct
  response is a nuisance-appropriate equivalence test.
- **Compute / platform.** CPU-only is the reproducible baseline; GPU (CUDA.jl) is an optional
  accelerator that must degrade gracefully. Windows-tauglich. Fixed 8×8 summary grid during the
  spike.

## Feature Areas

| Area | Blueprint | Spikes | What it covers |
|------|-----------|--------|----------------|
| **Bayes-factor ship-gate** | `references/bayes-factor-gate.md` | 006, 007, 008 | Attrition is baseline-side & boundary-triggered; survivorship bias in r̂; the clamp masks a real NRE-vs-KDE tail disagreement; log-space tail baseline recipe; the KDE reference is invalid past \|log-BF\|≈log(L). **008 PARTIAL** (remedy failed — that is the finding). |
| **SBC calibration of the NPE** | `references/sbc-calibration.md` | 009, 010, 011, 012, 013 | Shrinkage≠location; ρ_true "overconfidence" **retracted** → prior-atom artifact; randomized-rank recipe; targets calibrated; nuisance drift = M=2000 over-power; capacity redistributes (**012 PARTIAL**); μ-truncation costs high-\|ρ\| + ADVI. Recommended: randomized ranks + equivalence test + strict point-null for targets only. |
| **`src/` refactoring roadmap** | `references/src-refactoring-analysis.md` | 001, 002, 003, 004, 005 | Analysis-only (nothing implemented). Δρ bug (4/5 dimensions); hot-path function barriers (3.4–6.0×, 135–1394× mem); real type hierarchy; lift the `@model`; CairoMakie headless viz; seed+manifest provenance. Detail in `.planning/spikes/PROPOSAL.md`. |

## Processed Spikes (001–013)

| # | Name | Type | Verdict | Area |
|---|------|------|---------|------|
| 001 | hotpath-perf-memory | analysis | VALIDATED | src refactor |
| 002 | type-hierarchy-extensibility | analysis | VALIDATED | src refactor |
| 003 | bayes-model-modularity | analysis | VALIDATED | src refactor |
| 004 | visualization-reuse | analysis | VALIDATED | src refactor |
| 005 | provenance-reproducibility | analysis | VALIDATED | src refactor |
| 006 | bf-attrition-mechanism | standard | VALIDATED | BF gate |
| 007 | bf-survivorship-bias | standard | VALIDATED | BF gate |
| 008 | bf-baseline-remedies | comparison | **PARTIAL** | BF gate |
| 009 | shrinkage-vs-location | standard | VALIDATED | SBC |
| 010 | calibration-vs-imsize | standard | VALIDATED (parent README overconfidence claim **retracted** by ADDENDUM-atoms.md) | SBC |
| 011 | sbc-power-nuisance | standard | VALIDATED | SBC |
| 012 | capacity-lever | standard | **PARTIAL** | SBC |
| 013 | mu-prior-truncation | standard | VALIDATED | SBC |

## How to use this skill

1. Identify the area (BF gate / SBC / src refactor) and open the matching `references/*.md` — each is
   a build blueprint (Requirements / How to Build It / What to Avoid / Constraints / Origin).
2. For a specific claim, open `sources/NNN-*/README.md` (and the `.jl` scripts) for full evidence and
   run commands. **Note `sources/010-calibration-vs-imsize/ADDENDUM-atoms.md` retracts its parent
   README's conclusion — read the addendum as the settled finding.**
3. Honor the constraints above: PARTIAL verdicts (008, 012) are *not* solved build steps; every gate
   remedy is unapplied and a pre-registration matter; every result is single-DEV-seed unless noted.

## Related

- `.planning/spikes/MANIFEST.md`, `PROPOSAL.md`, `CONVENTIONS.md`, `WRAP-UP-SUMMARY.md`
- `.planning/phases/07-productionization-conditional-on-go/` — `07-GATE-AMENDMENT.md`,
  `07-NUISANCE-SBC-SPEC-DRAFT.md`
- User auto-memory: `phase7-amended-gate-diagnosis`, `src-refactoring-proposal-spike`,
  `phase5-gates-failed-npe-overconfident` (partially superseded).
