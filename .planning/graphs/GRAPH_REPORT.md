# Graph Report - ProteinCoLoc  (2026-06-28)

## Corpus Check
- 25 files · ~111,548 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 155 nodes · 180 edges · 20 communities (17 shown, 3 thin omitted)
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `e06414e4`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 14|Community 14]]

## God Nodes (most connected - your core abstractions)
1. `ProteinCoLoc` - 13 edges
2. `Pkg.add("CUDA")           # loaded as a Flux extension, separate from CPU path` - 13 edges
3. `calibrate()` - 9 edges
4. `3. SIM-02 Induced-μ Calibration — the prior-consistency proof (Plan 02-03)` - 9 edges
5. `Random` - 8 edges
6. `Distributions` - 6 edges
7. `Statistics` - 6 edges
8. `ProteinCoLoc` - 6 edges
9. `Spike Notes — ProteinCoLoc v2.0 (AmortizedColoc)` - 6 edges
10. `StatsBase` - 4 edges

## Surprising Connections (you probably didn't know these)
- None detected - all connections are within the same source files.

## Communities (20 total, 3 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.11
Nodes (15): Base, CairoMakie, CSV, DataFrames, GLMakie, Images, KernelDensity, QuadGK (+7 more)

### Community 1 - "Community 1"
Cohesion: 0.12
Nodes (7): Distributions, Flux, HypothesisTests, NeuralEstimators, Pkg, Random, Test

### Community 2 - "Community 2"
Cohesion: 0.12
Nodes (15): 1. Stack Decision (ENV-04 / D-07), 2. Julia-Version Record (ENV-03 / D-05), 4. Parent-Package Coupling Outcome (ENV-01 / D-01) — Plan 01-04, 5. Decoupling Proof (ENV-01 / DEMO-02 pattern) — Plan 01-04 Task 2, code:block1 (juliaup add 1.12.6), code:block2 (⌃ [38f6df31] ↓ NeuralEstimators v0.2.1 ⇒ v0.1.4), code:block3 (⌃ [38f6df31] NeuralEstimators v0.1.4 (<v0.2.1)), code:block4 (NeuralEstimators CPU smoke: Error During Test) (+7 more)

### Community 3 - "Community 3"
Cohesion: 0.27
Nodes (9): AbstractMultiChannelImage, _apply_mask!(), _calculate_mask(), channel_names(), image_data(), num_channels(), otsu_thresholds(), shuffle_blocks!() (+1 more)

### Community 4 - "Community 4"
Cohesion: 0.15
Nodes (13): Alternatives Considered, Architecture, Conventions, Developer Profile, GSD Workflow Enforcement, Key API Notes (AP1 smoke target), Pkg.add("CUDA")           # loaded as a Flux extension, separate from CPU path, Project Skills (+5 more)

### Community 5 - "Community 5"
Cohesion: 0.36
Nodes (11): calibrate(), emit_ghat(), _fmt(), _induced_mu_samples(), _isotonic_increasing(), _ks_stat(), _patchcorr_diagnostics(), _real_anchor() (+3 more)

### Community 6 - "Community 6"
Cohesion: 0.17
Nodes (11): Constraints, Core Technologies, Development Tools, Executive Verdict (read first), In ProteinCoLoc/spike/ — own environment, do NOT touch the root project, Installation, Optional GPU accelerator (only if training is slow; spike must run without it), Project (+3 more)

### Community 7 - "Community 7"
Cohesion: 0.2
Nodes (10): 3. SIM-02 Induced-μ Calibration — the prior-consistency proof (Plan 02-03), Chosen nuisance ranges (Claude's discretion, D-03), code:julia (mu    ~ Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)   # ← the SIM), σ / τ / ν consistency (D-03 — a CHECK, not a fit; the SIM-02 gate is on μ), Induced-vs-target μ match (the SIM-02 pass metric), Real anchor (D-16) — faithful extraction + negative-tail reachability, Size-invariance (T-02-SI, D-08/D-09), Sweep + span-past-±0.9 rationale (D-15) (+2 more)

### Community 8 - "Community 8"
Cohesion: 0.29
Nodes (6): CoordinateTransformations, ImageFiltering, ImageTransformations, Interpolations, simulate_pair(), _smooth_field()

### Community 9 - "Community 9"
Cohesion: 0.29
Nodes (6): Citation, Contact, Description, Features, License, ProteinCoLoc

## Knowledge Gaps
- **62 isolated node(s):** `NeuralEstimators`, `Flux`, `CairoMakie`, `ImageFiltering`, `ImageTransformations` (+57 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Random` connect `Community 1` to `Community 0`, `Community 8`, `Community 5`?**
  _High betweenness centrality (0.058) - this node is a cross-community bridge._
- **Why does `ProteinCoLoc` connect `Community 0` to `Community 1`?**
  _High betweenness centrality (0.044) - this node is a cross-community bridge._
- **Why does `Statistics` connect `Community 0` to `Community 8`, `Community 5`?**
  _High betweenness centrality (0.030) - this node is a cross-community bridge._
- **What connects `NeuralEstimators`, `Flux`, `CairoMakie` to the rest of the system?**
  _62 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.11 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.12 - nodes in this community are weakly interconnected._
- **Should `Community 2` be split into smaller, more focused modules?**
  _Cohesion score 0.12 - nodes in this community are weakly interconnected._