# Phase 7: Productionization (conditional on Go) - Context

**Gathered:** 2026-07-03
**Status:** Ready for planning

<domain>
## Phase Boundary

Promote the spike's validated amortized inference into `src/` as the package's **primary,
public colocalization API**, on the strength of the Phase-6 **Clean Go**. This is the **only
phase that edits `src/`**.

**Scope reframe locked with the user (overrides ROADMAP SC1 / PROD-01 backward-compat
clause):** v2.0 is a **new breaking release**. There is **no obligation to preserve the old
public API or keep `compute_BayesFactor` working unchanged.** The API may be redesigned
freely. The single hard design directive is that the new API be **extensible and
maintainable — clear types, an understandable type hierarchy** (see D-02).

**Turing/ADVI disposition — Option 3 (user-chosen):** amortized inference is the **only
supported public entry point**; the Turing/ADVI model and `compute_BayesFactor` are **kept in
the repo as internal reference code** (used by the D-05 ship-gate and future validation), not
as supported public API.

In scope: the `src/` amortized inference entry point (`colocalization_amortized()` or its
redesigned equivalent), the result type hierarchy, the estimator registry (5 grids), the
promotion/adaptation of spike inference code into `src/`, adding the amortized deps to the
root package, and running the all-grid ship-gate. **Out of scope:** the feature-DAG work
(Phases 8, 11–16); any new scientific capability beyond promoting what the spike validated;
retraining the *summary statistic* itself (the patch-correlation summary is reused).

Requirements: PROD-01, PROD-02.
</domain>

<decisions>
## Implementation Decisions

### API framing & Turing disposition **[USER-DECIDED — reframe]**
- **D-01:** **Breaking release, amortized-only public API (Option 3).** The public surface is
  the amortized path. The Turing/ADVI `colocalization()` + `compute_BayesFactor()` are
  demoted to **internal reference code** — retained in the repo, not exported as supported
  API — serving the D-05 ship-gate (as the non-amortized accuracy reference) and future
  validation/cross-checks. No backward-compatibility obligation to the pre-v2.0 API.

### Type hierarchy (the governing design directive) **[USER-DECIDED]**
- **D-02:** **Design an extensible, maintainable type hierarchy with clear types** — this is a
  first-class deliverable, not incidental. Concretely:
  - An **abstract `AbstractColocResult`** supertype with concrete subtypes
    **`AmortizedColocResult`** (posterior draws, NRE/amortized Bayes factor, OOD flag,
    calibration/diagnostic metadata) and a kept-internal **`AdviColocResult`** (the Turing
    path's result).
  - A **shared accessor interface** over the supertype — e.g. `delta_rho(r)`,
    `bayes_factor(r)`, `is_ood(r)`, `posterior_draws(r)` — so downstream/plotting code is
    **type-agnostic** and dispatches on the interface, not concrete fields.
  - The hierarchy must **leave clean extension points** for the roadmap's future result
    variants: the three-hypothesis BF (Phase 13) and the spatial per-region Δρ map
    (Phase 12) should slot in as **new `AbstractColocResult` subtypes / interface
    extensions**, not as bolted-on fields. The planner/researcher should sketch these
    extension points even though they are out of scope to implement here.

### Dependency / packaging **[USER-DECIDED]**
- **D-03:** **Add NeuralEstimators, Flux, and JLD2 (+ any transitive spike inference deps) as
  hard direct dependencies in the root `Project.toml`.** Amortized is the primary path, so it
  must work on a clean install with no optional-load machinery. Accepts the heavier base
  install (Flux load time) as the correct tradeoff for an amortized-first package. The spike's
  isolated `spike/Project.toml`/`Manifest.toml` remain untouched; this is a *new* root-package
  dependency set for the productionized code.

### Estimator registry (PROD-02) **[USER-DECIDED]**
- **D-04:** **Ship a five-grid pre-trained family — 4×4, 8×8, 16×16, 32×32, 64×64 —** in a
  registry keyed by patch grid, PLUS a documented **`train_and_register(grid)`** entry point
  so users can add further grids without ad-hoc retraining. The registry maps a requested
  `num_patches`/grid to its trained estimator; an unregistered grid gives a clear, actionable
  error pointing at `train_and_register`.
  - **RATIONALE (user):** the finer grids exist to enable **"local localisation" analysis** —
    a finer patch grid yields **per-region / spatially-resolved** colocalization, not merely a
    higher-dimensional summary. This is the on-ramp to the Phase-12 spatial Δρ map. The
    feasibility question for the fine grids is therefore **"does local colocalization
    resolution hold up at this grid?"**, not just "is the summary trainable" — the researcher
    must evaluate the fine grids against the *local-resolution* goal.
  - **Consequence for the planner:** the summary-vector dimension is **coupled to the grid**
    (8×8 → 64 continuous + 64 mask rows; other grids scale accordingly), so **each grid needs
    its own training-data generation, its own NPE (and NRE) architecture + training run**, not
    a reuse of the 8×8 net. This is a substantial, repeated pipeline — see the scope flag in
    Deferred Ideas.
  - Decide (research/discretion) **where the shipped `.jld2` estimators live** — Julia
    `Artifacts.toml` lazy-download vs. bundled-in-repo — given five grids up to 64×64 (4096
    patches) may be large.

### Ship-gate operationalization (D-05 carried from Phase 6) **[USER-DECIDED]**
- **D-05:** **Gate ALL five grids before ship.** Every bundled estimator must pass its **own
  independent, fresh-seed, re-pre-registered SBC/BF/OOD confirmation run** before the release
  ships — no provisional/unvalidated estimators in the public registry. The gate re-expresses
  the pre-registration **fresh** (new locked consts + new disjoint seed per the Phase-6 memo's
  ship-gate definition), it does **not** silently reuse the spike's `VAL_MASTER_SEED` (which
  the net has now been iterated against). Results are recorded per grid; a grid that fails its
  gate does not ship (and is not merged into `src/` as public).

### GPU acceleration for training **[USER-RAISED; recommended default, confirmable]**
- **D-06:** **GPU MAY accelerate the 5-grid training.** The CPU-only rule was a *spike*
  constraint (portability/reproducibility de-risking) and does **not** bind Phase 7. Per
  CLAUDE.md, GPU is an optional accelerator for training with graceful CPU fallback.
  **Recommended split (told to the planner; user confirmation pending):**
  - **Train on GPU** (NeuralEstimators `use_gpu=true` / CUDA.jl as a Flux extension) — most
    valuable for the heavy 32×32 / 64×64 grids; add CUDA.jl with graceful CPU fallback.
  - **Keep the frozen nets, the per-grid ship-gate SBC/BF/OOD confirmation runs, and the
    shipped default inference path CPU-reproducible** (graceful GPU→CPU), so the
    pre-registered ship-gate numbers stay deterministic/reproducible regardless of training
    hardware. GPU changes training-time, NOT the local-localisation *statistical* calibration
    question (that stands on its own merits, D-04 rationale).
  - If the user later prefers GPU-everywhere or GPU-training-only-ephemeral, the planner
    adjusts — this decision affects hardware/repro plumbing, not the type hierarchy or registry.

### Claude's Discretion
- Exact naming of the public entry point and accessor functions (within the D-02 hierarchy).
- Estimator storage mechanism (Artifacts vs bundled) and registry file format.
- How the internal Turing path is structured/namespaced to stay usable by the ship-gate
  without being exported.
- Whether Phase 7 is split into per-grid plans (strongly implied by D-04 + D-05 — see flag).
</decisions>

<specifics>
## Specific Ideas

- The user explicitly values **long-term maintainability** — the type hierarchy (D-02) is the
  headline concern, weighted above minimizing effort. Prefer explicit clarity over clever
  abstraction (the user picked an abstract-supertype hierarchy over a parametric
  `ColocResult{Backend}` precisely for readability).
- "New breaking release" means the researcher/planner should feel free to **name and shape the
  API from scratch** rather than contorting to the old `colocalization()` signature.
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 7: Productionization" — SC1–2, PROD-01/02. **NOTE:** SC1's
  "coexisting … so `compute_BayesFactor` keeps working unchanged" clause is **RELAXED by D-01**
  (breaking release; Turing path internal-only). PROD-02 (`num_patches` user-definable via
  registry) stands and is realized by D-04.
- `.planning/REQUIREMENTS.md` — PROD-01, PROD-02 definitions.

### The Go decision & the ship-gate this phase must honor
- `.planning/phases/06-reproducible-demo-go-no-go-memo/06-GO-NO-GO-MEMO.md` — the Clean Go
  verdict authorizing this phase and the **D-05 independent-confirmation ship-gate** definition.
- `.planning/phases/06-reproducible-demo-go-no-go-memo/06-CONTEXT.md` §D-05 — ship-gate intent.

### Existing src/ surface (the integration target — read to redesign, not to preserve)
- `src/ProteinCoLoc.jl` — module + `export` list to be **redesigned** for the breaking release.
- `src/bayes.jl` — `CoLocResult` struct (→ becomes internal `AdviColocResult`),
  `compute_BayesFactor` (internal reference), `_prepare_data` (shared summary input),
  `colocalization()` (internal Turing/ADVI path).
- `src/colocalization.jl`, `src/LoadImages.jl` — `MultiChannelImageStack`, patch/correlation
  summary the amortized path consumes.

### Spike inference code being promoted/adapted into src/
- `spike/npe/train_npe.jl`, `spike/npe/infer.jl`, `spike/npe/architecture.jl` — NPE training +
  amortized read surface + summary/DeepSet architecture (per-grid).
- `spike/validation/train_ratio.jl`, `spike/validation/bf.jl` — NRE / amortized Bayes factor.
- `spike/validation/ood.jl` — OOD flag (density + noise channels).
- `spike/validation/harness.jl`, `spike/validation/sbc.jl` — the θ*~π→simulate→infer +
  SBC harness the per-grid ship-gate re-uses.
- `spike/data/generate.jl`, `spike/data/cache.jl`, `spike/data/loader.jl` — per-grid
  training-data generation + leak-free loading.
- `spike/validation/consts.jl` — the pre-registration to **re-express fresh** per grid (D-05).
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **Promotable spike inference stack** (proven in Phases 1–6): NPE train/infer, NRE/BF, OOD,
  and the SBC harness — Phase 7 adapts these from `spike/` into `src/`, it does not rebuild them.
- **Existing summary contract in src/:** `_prepare_data(img, channels, num_patches; cor_method)`
  + the patch/correlation summary already produce the exact input the amortized nets were
  trained on — the amortized path reuses this input surface.
- **`CoLocResult` → `AdviColocResult`:** the existing struct (`img, control, channels,
  num_patches, posterior::DataFrame, advi_result`) becomes the internal Turing result subtype
  under the new `AbstractColocResult` hierarchy.

### Established Patterns
- **Grid ↔ summary-dimension coupling:** the patch grid determines the summary vector length
  (8×8 → 128-row summary with a 64/64 continuous/mask split). Each registry grid therefore
  needs its own architecture + training + calibration — the core reason D-04+D-05 make this a
  large phase.
- **Atomic JLD2 persistence + pre-registered-consts + disjoint-seed discipline** (from the
  spike) carry directly into the estimator registry storage and the per-grid ship-gate.

### Integration Points
- Amortized entry point consumes `_prepare_data` output; returns an `AmortizedColocResult`.
- The internal Turing path stays wired to the ship-gate as the non-amortized reference.
- New root-`Project.toml` deps (NeuralEstimators/Flux/JLD2) must co-resolve with the existing
  Turing/GLMakie/Images stack — a resolve-risk the researcher should check early.
</code_context>

<deferred>
## Deferred Ideas

- **⚠ Scope flag (planner): Phase 7 is large.** D-04 (5 grids) × D-05 (per-grid fresh ship-gate)
  = five full pipelines: per-grid training-data generation → NPE + NRE training → fresh
  re-pre-registered SBC/BF/OOD confirmation. This likely warrants **splitting into per-grid
  plans** (or a shared-infrastructure plan + five grid plans). The planner should size this
  realistically and consider whether some grids are sequenced rather than all-at-once.
- **⚠ Feasibility risk (researcher): the fine grids (esp. 64×64 = 4096 patches) for LOCAL
  localisation.** The fine grids are meant to deliver **per-region colocalization** (D-04
  rationale). At high resolution each patch is tiny and its per-patch correlation estimate is
  noisy, so the question is whether **calibrated local Δρ resolution** survives — not just
  whether the summary is trainable. The researcher must assess, per fine grid, whether local
  colocalization is genuinely calibratable before the planner commits to gating it, and surface
  any grid that can't meet the local-resolution bar to the user rather than shipping a weak
  estimator. (Note the design overlap with Phase 12's spatial map — flag if the local-grid
  approach here should inform or be informed by that phase.)
- **Re-enabling the OOD posterior-predictive channel** in the productionized OOD flag (the
  θ̂ finite-guard made it viable in iter1 but the reported OR-fusion still runs `with_pp=false`)
  — a hardening item to fold in here or defer.
- **RxInfer independent cross-check** of the ADVI reference — still a post-Go paper nice-to-have
  (BACK-01), not a Phase-7 dependency.

### Reviewed Todos (not folded)
None — `todo.match-phase 7` not run (no todo backlog matched in prior phases).
</deferred>

---

*Phase: 07-productionization-conditional-on-go*
*Context gathered: 2026-07-03 (interactive; breaking release, amortized-only public API, 5-grid all-gated)*
