# Phase 4: NPE Training + ADVI Benchmark + Ablation - Context

**Gathered:** 2026-06-30
**Status:** Ready for planning

<domain>
## Phase Boundary

Train a `PosteriorEstimator` (MLP summary net → `NormalisingFlow`) on the Phase-3 cached,
leak-free-standardized vectors to infer **ρ_true and Δρ** in a single forward pass; benchmark it
against the existing `colocalization()` ADVI for **RMSE / interval width** and a **>100× speedup
at comparable accuracy** (NPE-01/02/03); and run the **minimal (D-01) vs augmented (D-02)**
summary-statistic **ablation** scoring per-parameter RMSE as a gating sufficiency diagnostic
(ABL-01/02).

In scope: NPE architecture + training on the Phase-3 loader output; the dedicated-env ADVI
ground-truth baseline and its JLD2 hand-off; the speedup/accuracy benchmark (incl. complexity
scaling, CPU core-scaling, and an optional GPU add-on); the summary ablation and its
summary-choice justification. **Out of scope:** SBC, the amortized Bayes factor (`RatioEstimator`),
and the OOD flag (Phase 5 — built off one shared simulate→infer harness over these trained nets);
the demo + Go/No-Go memo (Phase 6); any `src/` edits or productionization (Phase 7). **No `src/`
edits** — the spike and the baseline reach `src/` only by read-only copy/`include()`.

Requirements: NPE-01, NPE-02, NPE-03, ABL-01, ABL-02.
</domain>

<decisions>
## Implementation Decisions

### ADVI ground-truth baseline (NPE-02/03)
- **D-01:** **Dedicated isolated baseline env, ported API — NOT the spike env, NOT `src/`.** The
  ADVI baseline runs in its **own isolated Julia environment** (e.g. `spike/baseline/` with its
  own `Project.toml`) carrying a **modern Turing/AdvancedVI**; the hierarchical `@model` is
  **copied/`include()`d read-only** and the **removed `vi(m, ADVI(int,int))` / `syms` calls are
  ported to the current API there**. It is NOT added to the spike env (Turing co-resolve would cap
  NeuralEstimators below the pinned 0.2.1 — the Phase-1 landmine) and NOT fixed in `src/` (spike
  decoupling holds until Phase 7).
- **D-02:** **ADVI runs offline → JLD2 artifact the spike consumes.** The baseline runs ADVI on
  the reserved holdout pairs and serializes `{per-parameter posterior summaries, interval widths,
  wall-clock}` to a **JLD2 artifact**; the spike's NPE benchmark reads that artifact for the RMSE
  and >100× comparison. This decouples the broken-Turing baseline from the pinned-NeuralEstimators
  spike. **The parent project's own env is also broken** (Manifest pins uninstalled Turing; `src`
  uses the removed API) — the baseline must NOT depend on repairing it.

### NPE targets — ρ_true and Δρ (NPE-01)
- **D-03:** **One NPE over per-stack θ; Δρ by differencing two single-stack passes.** A single
  `PosteriorEstimator` infers the per-stack θ (including ρ_true). The **Δρ posterior = Monte-Carlo
  difference of two independent single-stack posterior passes** (sample − control). Reuses the
  Phase-3 single-stack cache unchanged, stays fully amortized (2 forward passes), and assumes
  between-stack independence (matches the sample-vs-control experimental design). Avoids the
  paired-summary re-generation a joint estimator would require.
- **D-04:** **Δρ benchmark = paired reserved-holdout stacks, symmetric on both methods.** Sample/
  control **pairs are formed from the reserved holdout**; `Δρ_true = ρ_true(sample) −
  ρ_true(control)`. **Both** the NPE (difference of two passes) and ADVI (joint on the pair) run on
  **identical pairs** → apples-to-apples Δρ RMSE. The per-stack ρ_true benchmark stays per-stack.

### Summary network
- **D-05:** **MLP summary net (DeepSet deferred).** An **MLP conditioner** over the fixed
  standardized loader output (128-dim D-01 / 142-dim D-02) → NeuralEstimators `NormalisingFlow`.
  DeepSet (BACK-02) is **deferred** until MLP sufficiency is *measured*, not assumed. The D-01 mask
  rows (65:128) are part of the input vector and are bypassed by loader standardization (Phase-3
  D-08). Depth/width/activation = research/discretion.

### Summary ablation (ABL-01/02)
- **D-06:** **Ablate both cached variants under leak-free k=5 CV, zero re-simulation.** Train the
  NPE on the **minimal (D-01)** and **augmented (D-02)** summaries — both already cached (Phase-3
  D-02) — and compare **per-parameter RMSE** under the leak-free 5-fold CV the loader enforces
  (DATA-03).
- **D-07:** **ρ_true/Δρ-led, advisory-with-rule gate.** The headline decision is on **ρ_true (and
  Δρ) RMSE**, with **all 7 θ reported**. Choose augmented **only if it materially beats minimal on
  ρ_true RMSE by a pre-set margin**; otherwise keep minimal for **parsimony + OOD detectability**
  (couples to Phase 5). A **SBC-pass-but-high-RMSE outcome is treated as an insufficiency signal**,
  not a pass. Research sets the exact margin/tie-break.

### Speedup benchmark protocol (NPE-03)
- **D-08:** **Raw-stack → posterior on identical inputs; summary time counted; training excluded.**
  Time both methods from the **same raw `MultiChannelImage` stack to posterior**. NPE clock =
  **summary extraction + forward pass** (training is **excluded** — the amortized claim); ADVI
  clock = the full `vi()` run. CPU-only, `BenchmarkTools`, over the ≥20 reserved holdout stacks.
  **Speedup is always reported paired with RMSE, never alone.**
- **D-09:** **"Comparable RMSE" = pre-registered tolerance.** NPE per-parameter RMSE must fall
  within a **pre-set, pre-registered tolerance** of ADVI's (e.g. ≤1.2×; research fixes the exact
  value before running).
- **D-10:** **CPU gates; GPU is an optional bonus.** **All pass/fail criteria (>100×, RMSE) are
  decided on CPU.** GPU timings are reported **only if a CUDA device is present**, skipping
  gracefully otherwise (Flux CUDA extension, **no forced import**). Honors the CLAUDE.md hard
  CPU-only-baseline / GPU-optional / graceful-degrade constraint; keeps the pinned spike env clean.
- **D-11:** **KernelAbstractions.jl only if non-perturbing, summary/simulator kernels only.** KA is
  evaluated for **cpu/gpu auto-dispatch of the summary/simulator kernels only** (the network's
  dispatch is already handled by Flux), and adopted **only if it does not perturb the pinned spike
  env / NeuralEstimators 0.2.1**. Otherwise CPU/GPU dispatch stays Flux-native.
- **D-12:** **Complexity scaling — the >100× as a curve, both methods.** Characterize wall-clock vs
  **(a) number of datasets N** (the amortization story: NPE ~flat O(1)/dataset post-training, ADVI
  ~linear in N) and **(b) input size** (imsize / patch grid), for **both NPE inference and ADVI**;
  report empirical scaling exponents. The speedup is presented as a curve, not a single point.
- **D-13:** **CPU core-scaling — reported characterization, not a gate.** Sweep thread counts (e.g.
  1/2/4/8) reporting NPE (summary+forward) and ADVI wall-clock vs `Threads.nthreads()` on the
  reserved holdout; identify parallel speedup / Amdahl saturation. This is **characterization, not
  pass/fail** — the >100× headline is stated at a **fixed, reported thread count**.

### Claude's Discretion
- MLP depth/width/activation; flow type/coupling/transform stack; optimizer/LR/epochs/early-stopping/batch size.
- Exact comparable-RMSE tolerance (D-09), ablation margin/tie-break (D-07), thread-count sweep (D-13), and the N / imsize ranges for the scaling study (D-12) — all pre-registered before the reported run.
- Baseline env layout/location, the ADVI `iter`/`num_samples` defaults used for the baseline, and the JLD2 artifact schema (D-01/D-02).
- The exact KA evaluation depth (D-11) and how the optional GPU path is wired (D-10).
- Whether other θ parameters join ρ_true/Δρ in the headline table (all 7 are reported regardless).
</decisions>

<specifics>
## Specific Ideas

- **The amortized thesis IS a scaling claim** — present >100× as a *curve* over datasets N and
  input size (D-12) plus a CPU core-scaling sweep (D-13), not a single number.
- **KernelAbstractions.jl** was suggested for automatic cpu/gpu dispatch; scoped to the
  summary/simulator kernels (the net is already Flux-dispatched) and gated on not disturbing the
  pinned env (D-11).
- The ADVI baseline lives in a **self-contained ported env precisely because neither the spike env
  (NeuralEstimators 0.2.1 pin) nor the parent env (broken Turing) can run it directly** (D-01/D-02).
- Δρ is treated as a *paired contrast* derived by differencing independent single-stack posteriors
  (D-03), with the benchmark pairs drawn symmetrically for both methods (D-04).
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 4: NPE Training + ADVI Benchmark + Ablation" — goal + 5 success
  criteria (single-pass NPE for ρ_true & Δρ; leak-free CV RMSE/interval benchmark vs ADVI; >100×
  paired with comparable RMSE; per-parameter ablation as a gating sufficiency diagnostic; chosen
  summary justified with OOD interaction noted).
- `.planning/REQUIREMENTS.md` — NPE-01, NPE-02, NPE-03, ABL-01, ABL-02; BACK-02 (DeepSet deferral).

### The ADVI baseline being benchmarked (READ-ONLY — the comparison target, never edit in spike)
- `src/bayes.jl` — `colocalization()` (line 234, the hierarchical Turing `@model`);
  `vi(m, ADVI(num_latent, iter))` (line 325 — **REMOVED pre-0.35 Turing/AdvancedVI API, must be
  ported in the isolated baseline env**); `CoLocResult` (line 81); `_prepare_data` (line 164);
  `compute_BayesFactor` (line 109). Copy/`include()` read-only into the baseline env.
- `src/colocalization.jl` — `patch()` / `correlation(x,y; method)` / `_exclude_zero` — the summary
  the NPE input derives from (fixed 8×8 grid).
- **Known issue (port before the baseline runs):** the ADVI path uses removed Turing/AdvancedVI API
  AND the Manifest pins uninstalled versions. Analysis + fix sketch:
  `.planning/spikes/003-bayes-model-modularity/README.md` (INFER-1) and `.planning/spikes/PROPOSAL.md`
  (B2). Verify current line numbers before editing the baseline copy.

### Phase-3 inputs the NPE consumes
- `spike/data/loader.jl` (module `Loader`) — `load_fold(dir, fold; K=5, master_seed, variant=:min|:aug)`,
  `load_holdout(dir)`, `load_main_pool(dir)`, `all_folds` — the **sole** leak-free standardization
  path; per-fold standardized tensors are the direct NPE `train` input.
- `spike/data/encode.jl` — `encode_d01` (128-dim: 64 corr + 64 mask) / `encode_aug` (AUG_DIM=142)
  summary encodings (the two ablation variants).
- `spike/data/cache.jl` + the cache dir — `holdout.jld2` (reserved ≥20-stack ADVI benchmark set,
  disjoint from every fold), `shard_*.jld2` (CV pool), `meta.jld2` (content-hash guard).
- `.planning/phases/03-training-data-pipeline/03-CONTEXT.md` — D-01 (128-dim mask encoding),
  D-02 (both variants cached → zero-re-sim ablation), D-07/D-08 (loader leak-free, no global path),
  D-10 (reserved holdout disjoint), D-09 (k=5 deterministic folds).

### Constraints & stack
- `CLAUDE.md` §Constraints — spike **decoupling** (no `src/` edits until Phase 7), **CPU-only
  baseline / GPU-optional / graceful degrade**, reproducible-from-fixed-seed, **fixed 8×8 summary**
  — plus §Technology Stack / §Key API Notes (`PosteriorEstimator` + `NormalisingFlow` +
  `train`/`sampleposterior`/`assess`; Flux backend; **BenchmarkTools** for the >100× claim;
  "What NOT to Use").
- `.planning/phases/01-environment-smoke-gate/01-CONTEXT.md` — spike isolation, decoupling baseline
  `f581d95`, the read-only `include()` coupling, and the **NeuralEstimators<0.2.1 co-resolve
  landmine** that forbids adding Turing to the spike env (the reason D-01 isolates the baseline).
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `spike/data/loader.jl` (`Loader`): per-fold standardized `(theta=7×N, summary=128|142×N)` tensors
  — feed straight into NeuralEstimators `train`. `load_holdout` returns the reserved benchmark set.
- `spike/data/encode.jl`: the two summary variants (`:min` 128-dim, `:aug` 142-dim) the ablation
  (D-06) compares — already cached, no re-simulation.
- `spike/simulator/` (`sample_prior`, `simulate_pair`, `ghat`): drives the paired-holdout Δρ
  benchmark inputs (D-04) and any scaling-study input-size sweep (D-12).
- `src/bayes.jl` `colocalization()` / `CoLocResult` / `_prepare_data`: the ADVI baseline target,
  copied read-only into the isolated baseline env (D-01).

### Established Patterns
- Spike reaches `src/` only via read-only `include()`/copy; byte-identical baseline `f581d95`
  (Phase-1 D-01). The baseline env follows the same read-only-coupling discipline (D-01).
- Re-runnable quantitative gates under the spike `Test` harness (`spike/test/`) — the benchmark and
  ablation should land as re-runnable gates there.
- Random123 counter-based reproducibility (Phase-3 D-11) — the benchmark/scaling runs stay
  seed-reproducible.

### Integration Points
- Phase-3 loader/holdout → NPE `train` (CV folds) + NPE/ADVI benchmark (reserved holdout).
- Isolated baseline env → JLD2 ADVI artifact → spike benchmark (the only cross-env hand-off; D-02).
- The trained NPE(s) + chosen summary (D-07) feed **Phase 5** (SBC / amortized BF / OOD over one
  shared simulate→infer harness); the OOD-detectability note on the summary choice is the explicit
  coupling.
</code_context>

<deferred>
## Deferred Ideas

- **SBC, amortized Bayes factor (`RatioEstimator`), OOD flag** — Phase 5. In particular **Δρ
  correctness vs `compute_BayesFactor()` is reproduced in BF-02 (Phase 5)**, not here.
- **DeepSet permutation-invariant summary (BACK-02)** — upgrade only if the MLP (D-05) proves
  insufficient by the ablation.
- **RxInfer.jl independent cross-check baseline (BACK-01)** — post-Go, never a spike dependency.
- **Productionization / user-definable `num_patches`** — Phase 7 (conditional on a Go).
- **Making GPU a required deliverable / putting CUDA in the reproducibility path** — out of scope;
  CPU is the gating baseline, GPU is an optional bonus (D-10).

### Reviewed Todos (not folded)
None — no pending todos matched this phase.
</deferred>

---
*Phase: 04-npe-training-advi-benchmark-ablation*
*Context gathered: 2026-06-30*
