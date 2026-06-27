# Phase 3: Training-Data Pipeline - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-27
**Phase:** 03-training-data-pipeline
**Mode:** discuss (`--all` — all gray areas auto-selected, discussed interactively)
**Areas discussed:** Summary schema & ablation provisioning, Cache format / resumability / invalidation, Leak-free split discipline (loader), Seeding / parallelism / θ-sampling

---

## Summary schema & ablation provisioning (DATA-01)

### Missing-patch encoding
| Option | Description | Selected |
|--------|-------------|----------|
| Impute + 64-dim mask channel | Fill missing→0, append 64-dim present/absent mask (128-dim total); net distinguishes true-zero from missing | ✓ |
| Impute only, keep 64-dim | Fill missing→0; honors literal 64-dim contract; net can't tell missing from zero | |
| Require complete (reject missing) | Discard pairs with any missing patch; biases toward dense images | |

**User's choice:** Impute + 64-dim mask channel → D-01 (summary becomes 128-dim).

### Ablation provisioning
| Option | Description | Selected |
|--------|-------------|----------|
| Cache both summary variants now | Store minimal + augmented moments in one pass; Phase-4 ablation needs no re-simulation | ✓ |
| Cache minimal only | Smaller cache; Phase 4 re-runs generator for augmented arm | |

**User's choice:** Cache both summary variants now → D-02.

### Image size
| Option | Description | Selected |
|--------|-------------|----------|
| Fixed 256², spot-check real size | Fast, CPU-feasible; lean on 8×8 scale-invariance | |
| Fixed 1376×1028 (real anchor) | Maximal realism; ~27× pixel cost, CPU-prohibitive at scale | |
| Sample across sizes | Draw imsize from a range (256²–2048²); size-robust net; highest cost | ✓ |

**User's choice:** Sample across sizes → D-03 (imsize is a per-sample sampled nuisance; budget/parallelism must absorb the cost).

---

## Cache format, resumability & invalidation (DATA-02)

### Cache layout
| Option | Description | Selected |
|--------|-------------|----------|
| Sharded chunks (e.g. 10k/shard) | Resume = skip completed shards; parallel writers; bounded memory | ✓ |
| Single file + periodic checkpoint | Simpler load; coarser resume; single writer | |
| One JLD2 per sample | Finest resume; 200k tiny files, Windows-hostile | |

**User's choice:** Sharded chunks → D-04.

### Auto-invalidation guard
| Option | Description | Selected |
|--------|-------------|----------|
| Hash source + config | Content-hash simulator/prior/ghat source + config; any code or config change invalidates | ✓ |
| Hash config only | Misses source edits | |
| Manual version integer | Relies on remembering to bump | |

**User's choice:** Hash source + config → D-05 (satisfies SC-2).

### Sample budget
| Option | Description | Selected |
|--------|-------------|----------|
| Default 50k, N parameterized, sharded scale-up | Start at 50k; add shards to extend without regenerating | ✓ |
| Generate 200k upfront | Heavy cost before knowing 50k suffices | |
| Adaptive to Phase-4 accuracy | Most sample-efficient; couples Phase 3 to Phase-4 loop | |

**User's choice:** Default 50k, N parameterized, sharded scale-up → D-06.

---

## Leak-free split discipline — the loader (DATA-03)

### Standardization location
| Option | Description | Selected |
|--------|-------------|----------|
| Store raw; loader standardizes per fold | Cache standardization-free; z-stats fit on train folds only | ✓ |
| Store globally-standardized summaries | Bakes global stats into every fold — the prohibited leakage | |
| Store raw + global-stats sidecar | Invites accidental global standardization | |

**User's choice:** Store raw; loader standardizes per fold → D-07.

### Structural leak-guard mechanism
| Option | Description | Selected |
|--------|-------------|----------|
| Loader owns the split; no global path | Single API, stats fit-on-train/applied-to-val; no global-standardize function to misuse | ✓ |
| Convention + leakage regression test | Global path exists, guarded by a test; footgun remains | |

**User's choice:** Loader owns the split; no global path → D-08 (+ D-09 k=5-fold, deterministic from seed).

### Reserved ADVI-benchmark holdout
| Option | Description | Selected |
|--------|-------------|----------|
| Separate seed stream, i.i.d. from prior | Disjoint Random123 range, excluded from all folds, sampled from π(θ) | ✓ |
| Separate seed stream, stratified ρ grid | Even RMSE coverage; not prior-representative | |
| Carved from same pool post-hoc | Risks fold overlap; not structurally isolated | |

**User's choice:** Separate seed stream, i.i.d. from prior → D-10 (≥20 stacks).

---

## Seeding, parallelism & θ-sampling (reproducibility)

### Per-sample seeding
| Option | Description | Selected |
|--------|-------------|----------|
| Counter-based per-sample stream | RNG = f(master_seed, sample_index); reproducible regardless of order/thread count | ✓ |
| Sequential Xoshiro jumps | Order-dependent; parallel scheduling changes results | |
| Per-shard seed only | Samples not independently addressable | |

**User's choice:** Counter-based per-sample stream → D-11 (Random123, adopting Phase-2 D-14 deferral).

### Parallelism
| Option | Description | Selected |
|--------|-------------|----------|
| Multithreaded, serial fallback | Threads.@threads; results identical regardless of thread count | ✓ |
| Distributed multiprocess | Scales beyond one node; heavier setup | |
| Serial only | Simplest; slow at 50k–200k | |

**User's choice:** Multithreaded, serial fallback → D-12 (CPU-only portable baseline; Distributed reserved).

### θ / ρ_true sampling
| Option | Description | Selected |
|--------|-------------|----------|
| i.i.d. from prior π(θ) | Training distribution stays π(θ) (SBC depends on it); negative-tail sparsity documented as caveat | ✓ |
| Stratified over ρ_true | Even tail coverage; distorts π(θ), needs importance weighting | |
| i.i.d. prior + weighted tail oversample | Coverage + calibration; weight-handling complexity | |

**User's choice:** i.i.d. from prior π(θ) → D-13 (D-16 negative-tail caveat documented, not engineered away).

---

## Claude's Discretion

- Exact shard size (D-04); imsize set / sampling distribution over sizes (D-03 / D-09).
- Exact augmented-moment composition for the second summary variant (D-02).
- Internal generator/loader module layout within `spike/` and the JLD2 key schema per shard.
- Exact hash function and source-file set feeding the D-05 content hash.
- Mask dtype and 128-dim tensor layout handed to NeuralEstimators (D-01).

## Deferred Ideas

- NPE/NRE training, ADVI benchmark, summary ablation — Phase 4 (this phase provisions D-02 + D-10).
- SBC / amortized BF / OOD — Phase 5.
- Distributed.jl multiprocess generation — reserved (D-12).
- Fixing negative-correlation-tail sparsity (stratification / importance weighting) — explicitly not done (D-13).
