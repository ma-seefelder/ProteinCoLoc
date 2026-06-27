# Phase 3: Training-Data Pipeline - Context

**Gathered:** 2026-06-27
**Status:** Ready for planning

<domain>
## Phase Boundary

A **reproducible, resumable generator** that turns the calibrated prior π(θ) + the Phase-2
`simulate_pair` simulator + the **unchanged** existing summary functions into a **cached,
standardized** set of training vectors, with **leak-free split discipline baked into the
loader structurally** (not by per-script convention).

In scope: the θ~π → `simulate_pair` → `patch_summary`/summary → fixed-dimension summary-vector
generator (DATA-01); the version/hash-guarded, resumable, Random123-seeded JLD2 cache at
50k–200k scale (DATA-02); and the leak-free k-fold loader + reserved ADVI benchmark holdout
(DATA-03). **Out of scope:** NPE/NRE training, the ADVI benchmark itself, and the
summary-statistic ablation (Phase 4); SBC / amortized BF / OOD (Phase 5). **No `src/` edits** —
the generator reaches `src/` only via the read-only `include()` coupling in `spike/contract.jl`.

Requirements: DATA-01, DATA-02, DATA-03.
</domain>

<decisions>
## Implementation Decisions

### Summary schema — what gets cached per sample (DATA-01)
- **D-01:** **Impute + 64-dim mask channel.** The 8×8 `patch_summary` matrix (which can contain
  `missing` from `_exclude_zero`d background patches) is encoded as a **128-dim** summary
  vector: 64 correlation values with missing→0 imputation, plus a **64-dim binary present/absent
  mask** so the net can distinguish a genuine zero correlation from a missing patch. This
  deliberately extends the literal "8×8 = 64" contract to 64 values + 64 mask = 128; the
  **8×8 patch grid itself stays fixed** (Phase-2 D-10 / CLAUDE.md). Document the mask convention
  explicitly so Phase-4's summary-net input dimension is unambiguous.
- **D-02:** **Cache BOTH summary variants now.** In one generation pass, store per sample the
  minimal patch-correlation summary (D-01) **and** the augmented-moment summary
  (+Manders / median / IQR moments) that Phase-4's ablation (ABL-01) needs. Simulation is the
  expensive step; extra summary floats are cheap — so Phase 4 runs its ablation with **zero
  re-simulation**. (Exact augmented-moment set is a research/planning detail; the *commitment to
  cache both* is locked.)
- **D-03:** **Sample image size across a range.** `imsize` is a **per-sample sampled nuisance**
  drawn from a range of sensible microscopy sizes (e.g. 256²–2048², the D-09 size sweep), not
  fixed at the 256² `simulate_pair` default nor pinned to the 1376×1028 real anchor. This trains
  a **size-robust** net and exercises the 8×8 grid's near-scale-invariance directly. Cost
  implication is real (large sizes dominate sim time) — the budget/parallelism decisions below
  must absorb it. (Exact size set / sampling distribution over sizes = Claude's discretion /
  researcher, per D-09.)

### Cache format, resumability & invalidation (DATA-02)
- **D-04:** **Sharded JLD2 chunks** (e.g. ~10k samples/shard, BayesInteractomics pattern).
  **Resume = skip completed shards, restart only the in-flight one.** Enables parallel writers
  and bounded memory. (Exact shard size = Claude's discretion / tuneable.)
- **D-05:** **Version/hash guard hashes source + config.** A content hash over the
  simulator/prior/`ghat` **source** PLUS the generating **config** (imsize-range, summary spec
  incl. the D-01 mask, θ-dim, master seed) is stored with the cache; **any code OR config change
  auto-invalidates** stale data. Directly satisfies Success Criterion 2 ("a changed
  simulator/prior auto-invalidates stale data"). A config-only or hand-bumped-integer guard was
  rejected as too weak (would silently reuse stale vectors after a source edit).
- **D-06:** **Default budget 50k, N parameterized, sharded scale-up.** Generate **50k by
  default**; `N` is a single config knob; **adding more shards extends the dataset without
  regenerating** existing ones (composes with D-04). Scale toward 200k when Phase-4 accuracy
  demands. Not a one-shot 200k upfront; not auto-coupled to the Phase-4 training loop.

### Leak-free split discipline — the loader (DATA-03)
- **D-07:** **Store RAW summaries + θ; the loader standardizes per fold.** The cache is
  **standardization-free**: it holds raw summary vectors (both variants, D-02) and raw θ. The
  loader fits z-scoring mean/std **on the training folds only** and applies them to the held-out
  fold. This is the leak-safe foundation DATA-03 requires. Storing globally-standardized
  summaries (or a global-stats sidecar) is **rejected** — it bakes global statistics into every
  fold, the exact leakage the requirement prohibits.
- **D-08:** **The loader owns the split; there is NO public global-standardize path.** A single
  loader API returns per-fold standardized tensors with stats fit-on-train / applied-to-val.
  Because no global-standardize function is exposed to misuse, **cross-fold leakage is impossible
  by construction** — satisfying SC-3's "enforced in the loader, not by per-script convention."
  (A convention-plus-regression-test approach was rejected: it leaves the footgun in place.)
- **D-09:** **k-fold cross-validation, k = 5** (per ROADMAP SC-4 "e.g. 5-fold"). Fold assignment
  is **deterministic from the master seed** so reported NPE metrics are reproducibly
  cross-validated. (k is configurable; 5 is the default.)
- **D-10:** **Reserved ADVI-benchmark holdout: separate seed stream, i.i.d. from prior.** A
  reserved set of **≥20 stacks** for the Phase-4 ADVI speed/accuracy benchmark (NPE-02/03) is
  drawn from its **own disjoint Random123 counter range**, **excluded from every CV fold by
  construction**, and sampled **i.i.d. from π(θ)** (honest typical-case benchmark). Not carved
  post-hoc from the main pool (which would risk fold overlap), and i.i.d.-from-prior rather than
  a stratified grid (keeps it prior-representative).

### Seeding, parallelism & θ-sampling (reproducibility)
- **D-11:** **Random123 counter-based per-sample seeding** (adopting Phase-2 D-14's deferral).
  Each sample's RNG = `f(master_seed, global_sample_index)` — independent and reproducible
  **regardless of execution order or thread count**; the reserved set (D-10) occupies a disjoint
  counter range. This also addresses the Phase-1 thread-reproducibility note (WR-02). The
  simulator already accepts an explicit `rng::AbstractRNG` (Phase-2 D-14), so this slots in
  without simulator changes.
- **D-12:** **Multithreaded generation with serial fallback.** `Threads.@threads` over
  samples/shards; counter-based seeding (D-11) makes results **identical regardless of thread
  count**, degrading to serial on a single thread. Portable CPU-only baseline (CLAUDE.md);
  Distributed.jl multiprocess is reserved for if single-node threading proves too slow.
- **D-13:** **θ sampled i.i.d. from the calibrated prior π(θ).** The training distribution stays
  **exactly π(θ)**, which SBC calibration (Phase 5) depends on. The **negative-correlation tail
  sparsity flagged in Phase-2 D-16 is documented as a known caveat, NOT engineered away** by
  stratification (which would distort π(θ) and require importance weighting to stay calibrated).

### Claude's Discretion
- Exact shard size (D-04) and the imsize set / sampling distribution over sizes (D-03 / D-09).
- The exact augmented-moment set for the second summary variant (D-02) — Manders/median/IQR
  composition is indicative; the locked commitment is "cache both variants."
- Internal module/file layout of the generator + loader within `spike/` (e.g.
  `spike/data/`), and the precise JLD2 key schema inside each shard.
- The exact hash function and which source files feed the D-05 content hash.
- Mask dtype / how the 128-dim vector is laid out in the tensor handed to NeuralEstimators.
</decisions>

<specifics>
## Specific Ideas

- The summary contract is already proven end-to-end in `spike/contract.jl`:
  `patch_summary(mci)` → 8×8 `correlation(...; method=:pearson)` matrix (UNCHANGED `src/`),
  `induced_mu(mci)` → `mean(skipmissing(...))`. The D-01 128-dim encoding wraps this matrix; it
  does not modify it.
- The simulator API the generator drives: `sample_prior(rng)` → θ, `simulate_pair(rng, θ; imsize)`
  → `MultiChannelImage` pair, with `ghat(μ)` (in generated `spike/simulator/ghat.jl`) the
  calibrated μ↔ρ_true inverse. θ is the Phase-2 D-04 7-vector
  `(ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise)`.
- BayesInteractomics is the cited reference implementation for the JLD2 sharded-cache pattern
  (Phase-1 code_context) and for `_bin_calibration`/`CalibrationResult` (Phase-5, not here).
- Carry-forward (hard constraint): spike isolation — nothing under `src/` or the root manifests
  may change (decoupling baseline `f581d95`). The generator reaches `src/` only through the
  read-only `include()` coupling.
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 3: Training-Data Pipeline" — goal + 4 success criteria
  (generator mapping, version/hash-guarded resumable JLD2 cache, fit-on-train-only
  standardization enforced in the loader, 5-fold CV + ≥20-stack reserved ADVI holdout).
- `.planning/REQUIREMENTS.md` — DATA-01, DATA-02, DATA-03 (and the downstream NPE-02/03,
  ABL-01 these decisions provision for).

### The simulator + summary the generator drives (Phase 2 outputs)
- `spike/contract.jl` — the frozen summary contract: `build_mci`, `patch_summary`
  (8×8 `correlation(...; :pearson)`, UNCHANGED `src/`), `induced_mu`. The D-01 128-dim encoding
  wraps `patch_summary`'s matrix.
- `spike/simulator/forward.jl` — `simulate_pair(rng, θ; imsize)` (the per-sample physics);
  `spike/simulator/prior.jl` — `sample_prior(rng)` / π(θ);
  `spike/simulator/calibration.jl` + generated `spike/simulator/ghat.jl` — `ghat(μ)` the
  calibrated μ↔ρ_true inverse.
- `.planning/phases/02-forward-simulator-summary-contract/02-CONTEXT.md` — Phase-2 decisions,
  esp. D-04 (7-param θ), D-09 (configurable image size / size sweep), D-10 (fixed 8×8 grid),
  D-14 (explicit `rng`, Random123 deferred to here), D-16 (negative-tail reachability caveat).

### Summary & data contract in src/ (READ-ONLY, never edit)
- `src/colocalization.jl` — `patch()`, `correlation(x,y; method)` (line 221), `_exclude_zero`
  (source of the missing patches D-01 handles).
- `src/LoadImages.jl` — `MultiChannelImage` struct/constructor + `pixel_size = size(data[1])`
  convention (line 444); what a valid simulated image must satisfy.

### Constraints & stack
- `CLAUDE.md` §Constraints + §Technology Stack — spike decoupling (hard), CPU-only baseline /
  GPU-optional, reproducibility-from-fixed-seed, **JLD2** for the cache, **Random123** for
  counter-based seeding, **StatsBase** for standardization/standardize, fixed 8×8 summary during
  the spike, and the "What NOT to Use" list.
- `.planning/phases/01-environment-smoke-gate/01-CONTEXT.md` — spike isolation, decoupling
  baseline `f581d95`, the `include()` coupling, BayesInteractomics JLD2-cache reference.
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets (used UNCHANGED — no src/ edits)
- `spike/contract.jl`: `build_mci`, `patch_summary`, `induced_mu` — the summary the generator
  caches (D-01 wraps these).
- `spike/simulator/`: `sample_prior`, `simulate_pair`, `ghat` — the per-sample generation chain.
- `src/colocalization.jl` `correlation`/`patch`/`_exclude_zero`, `src/LoadImages.jl`
  `MultiChannelImage` — reached read-only via `contract.jl`'s `include()`.

### Established Patterns
- Spike reaches `src/` via `include()` (Phase-1 D-01 fallback); spike-local only, byte-identical
  baseline `f581d95`.
- Explicit-`rng`-threaded reproducibility (Phase-2 D-14) → now upgraded to Random123
  counter-based (D-11).
- JLD2 sharded cache + version guard modeled on BayesInteractomics (Phase-1 code_context).
- Validation = re-runnable quantitative gate under the spike `Test` harness (`spike/test/`).

### Integration Points
- Generator output (raw summary vectors + θ, both variants) → sharded JLD2 → the leak-free
  k-fold loader (D-07/D-08) → consumed by Phase-4 NPE `train`. The loader's per-fold standardized
  tensors are the direct input to NeuralEstimators.
- The reserved ADVI holdout (D-10) is consumed by Phase-4 NPE-02/03 (NPE-vs-ADVI benchmark) and
  must never enter a CV fold.
</code_context>

<deferred>
## Deferred Ideas

- **NPE/NRE training, the ADVI benchmark, the summary-statistic ablation** — Phase 4 (this phase
  only *provisions* the cached augmented summary, D-02, and the reserved holdout, D-10).
- **SBC / amortized Bayes factor / OOD flag** — Phase 5 (the `_bin_calibration`/`CalibrationResult`
  reuse lands there, not here).
- **Distributed.jl multiprocess generation** — reserved (D-12) unless single-node threading is
  too slow.
- **Fixing the negative-correlation tail sparsity** (stratification / importance weighting) —
  explicitly NOT done (D-13); documented as a caveat, revisit only if Phase-4/5 calibration
  demands it.

### Reviewed Todos (not folded)
None — no pending todos matched this phase.
</deferred>

---
*Phase: 03-training-data-pipeline*
*Context gathered: 2026-06-27*
