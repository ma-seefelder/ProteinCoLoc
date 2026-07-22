# Spike Manifest

## Idea

Identify refactoring targets across the **ProteinCoLoc v1.0 `src/` core** with three goals:
**(1)** improve performance (runtime + memory), **(2)** improve extensibility & code reuse
(e.g. a proper structure hierarchy), **(3)** ease maintenance and improve data visualization
and provenance tracking. This is an **analysis-only** spike: every target is proposed with
concrete before/after sketches and evidence, but **nothing is implemented**.

Method: each of the five dimensions below was analyzed by an independent agent that read the
in-scope source and produced structured findings; **every finding was then adversarially
verified against the actual code** by a second agent (46 agents total). Performance claims for
the hot path are backed by a standalone micro-benchmark (`001-hotpath-perf-memory/`).

## Scope

- **In scope (src/ core):** `ProteinCoLoc.jl`, `LoadImages.jl`, `colocalization.jl`, `bayes.jl`,
  `utils.jl`, `main.jl`, `plot.jl`, `script.jl`.
- **Out of scope:** `gui.jl`, `gui_css.jl`, `compile.jl`, and the `spike/` AmortizedColoc (v2.0) code.

## Requirements (constraints the proposals must honor)

- **Analysis only** — propose, do not implement. (User directive.)
- **Behavior-preserving** — every proposal preserves numerical behavior unless it is explicitly
  flagged as a correctness-bug fix.
- **Correct file map** (the analysis brief initially mislabeled this; verifiers corrected it):
  - `colocalization.jl` (≤ line 240) contains **only** `patch` / `unpatch` / `_exclude_zero` / `correlation`.
  - `bayes.jl` contains `convert_posterior_samples`, `CoLocResult` (struct **and** construction),
    `compute_BayesFactor`, `_prepare_data`, the `colocalization()` function, and the inline `@model`.
- **Decoupling caveat** — CLAUDE.md forbids touching `src/` during the active v2.0 spike. These
  `src/` refactors are a **separate effort**: sequence them either *before* resuming the spike or
  during the v2.0 **productionization** phase; do not interleave with the decoupled spike work.
- **Stack reality** — the resolved/installed Turing stack (0.43–0.45) differs from the legacy
  ADVI/`vi`/`DynamicPPL.syms` API the posterior path is written against; see Spike 003 / INFER-1.

## Spikes

| #   | Name | Type | Validates | Verdict | Tags |
|-----|------|------|-----------|---------|------|
| 001 | hotpath-perf-memory | analysis | Hot path `patch→correlation→_prepare_data` + KDE/BF: pinpoint allocation/type hotspots & propose redesign | ✅ VALIDATED | performance, memory, type-stability |
| 002 | type-hierarchy-extensibility | analysis | Image/result type system: propose a real interface + structure hierarchy | ✅ VALIDATED | extensibility, types, reuse |
| 003 | bayes-model-modularity | analysis | Inline `@model` + inference flow: propose modular, runnable, testable inference | ✅ VALIDATED | bayes, modularity, correctness |
| 004 | visualization-reuse | analysis | `plot.jl`+`utils.jl` plotting: propose reusable, themed, headless-safe layer | ✅ VALIDATED | visualization, dry, maintenance |
| 005 | provenance-reproducibility | analysis | Outputs + result struct: propose seed/version/param capture & provenance | ✅ VALIDATED | provenance, reproducibility, output |
| 006 | bf-attrition-mechanism | standard | BF-gate attrition: which side produces the non-finite log-BF, and is it random? | ✅ VALIDATED | bayes-factor, ship-gate, kde, numerics |
| 007 | bf-survivorship-bias | standard | Does the survivors-only filter bias the gate's correlation estimate? | ✅ VALIDATED | bayes-factor, survivorship-bias, statistics |
| 008 | bf-baseline-remedies | comparison | clamped vs unclamped vs log-space KDE baseline: is there one that is finite AND unsaturated? | ⚠ PARTIAL | bayes-factor, kde, baseline, remedy |
| 009 | shrinkage-vs-location | standard | Is the spillover rejection at shrinkage 1.01 a bug, or explainable? | ✅ VALIDATED | sbc, calibration, diagnostics |
| 010 | calibration-vs-imsize | standard | Why does ρ_true calibration collapse under the image-size mixture? | ✅ VALIDATED | sbc, overconfidence, imsize |
| 011 | sbc-power-nuisance | standard | Is the M=2000 SBC KS test over-powered for non-identified nuisances? | ✅ VALIDATED | sbc, power-analysis, over-power |

**Verdict legend:** VALIDATED = the refactoring opportunity is real, evidence-backed, and worth doing.

## Correctness bugs found in passing (separate from refactoring)

| Bug | Where | Confirmed by | Severity |
|-----|-------|--------------|----------|
| Δρ mixes posterior μ_sample with **prior** μ_control → wrong effect size in `result.txt` | `utils.jl:368` | 4 of 5 dimensions ✓✓✓✓ | **High** |
| Posterior path uses removed Turing `ADVI(n,iter)`/`vi(m,alg)`/`DynamicPPL.syms` API | `bayes.jl:323-329`, `bayes.jl:44` | bayes ✓ | **High** (may not run) |
| `convert_posterior_samples` slices magic `1:10`, mislabels 2 local params → posterior CSV has 10 mislabeled cols vs prior's 8 | `bayes.jl:42-50` | bayes ✓ | Medium |
| `plot()` `minmax_norm!` mutates shared image data → inference silently depends on plot flags | `plot.jl:149` | perf ✓ | Low |
| Fisher-z docstring claims a transform that is never applied (dead `#tanh` line) | `bayes.jl:40,52-53` | bayes ✓ (reclassified: stale-doc/trap, not a live bug) | Low |

See `PROPOSAL.md` for the consolidated, prioritized roadmap and `00N-*/README.md` for full detail.
