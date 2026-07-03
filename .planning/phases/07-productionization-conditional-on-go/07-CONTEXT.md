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

### Estimator registry (PROD-02) **[USER-DECIDED — REVISED after feasibility verdict]**
- **D-04 (revised):** **Ship a CAPPED grid family, plus a windowed local-map on-ramp.** The
  original 5-grid ambition was reframed after the research feasibility verdict (07-RESEARCH.md
  §FEASIBILITY VERDICT), which found that (a) a finer grid lengthens the *input* but the NPE
  output stays a **single GLOBAL ρ_true** — finer grids do **not** produce a per-region map
  (that is Phase 12's job), and (b) tiny fine-grid patches fall below `src/`'s `≥15-survivor`
  floor on normal-size images → calibrated-but-uninformative. The user's chosen scope:
  - **Ship 4×4, 8×8, 16×16** as calibrated **global** estimators in a registry keyed by patch
    grid, plus a documented **`train_and_register(grid)`** entry point (unregistered grid →
    clear error pointing at it).
  - **32×32** ships only **with a loud large-image caveat** (informative only on big, dense
    images), conditional on its ship-gate.
  - **64×64 is DROPPED** from the shipped family — not shipped as a "local map." The fine-grid
    ambition (true per-region Δρ maps) is **deferred to Phase 12** (GP/CAR spatial lattice).
  - **Local-localisation on-ramp shipped NOW:** run the validated **8×8 estimator on image
    sub-tiles (windowed inference)** to produce a genuine (coarse) local colocalization map —
    **no fine-grid training required**. This is the honest local-localisation feature for
    Phase 7; it dovetails with, rather than pre-empts, Phase 12.
  - The summary-vector dimension is **coupled to the grid**, so each shipped grid still needs
    its own data-gen + NPE + NRE + training run. Storage mechanism (Artifacts lazy-download vs
    bundled; recommend `Flux.state`+`loadmodel!` over whole-object `jldsave` for version
    robustness) is research/discretion.

### Ship-gate operationalization (D-05 carried from Phase 6) **[USER-DECIDED — scope follows D-04]**
- **D-05 (revised):** **Gate every SHIPPED grid before ship.** Each estimator that ships
  (4×4, 8×8, 16×16, and 32×32 if included) must pass its **own independent, fresh-seed,
  re-pre-registered SBC/BF/OOD confirmation run** (new locked consts + new disjoint
  `PROD_SEED[grid]` per the Phase-6 memo's ship-gate definition — **not** the spike's
  `VAL_MASTER_SEED`). A grid that fails its gate does not ship. 32×32 is explicitly
  **conditional** on its gate outcome + the large-image caveat; 64×64 is not gated because it
  is not shipped. The windowed 8×8-sub-tile local map inherits the 8×8 estimator's gate (no
  separate net to gate).

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

- **RESOLVED — feasibility verdict (07-RESEARCH.md §FEASIBILITY VERDICT).** The fine grids do
  NOT deliver local localisation with the current architecture (finer input → still a single
  global ρ; per-region maps are Phase 12). 64×64 additionally collapses to uninformative on
  normal-size images. **User's decision (2026-07-03): cap the family (ship 4×4/8×8/16×16, 32×32
  caveated, DROP 64×64) + ship the 8×8-sub-tile windowed coarse local map as the Phase-7
  local-localisation feature; defer true per-region maps to Phase 12.** Captured in D-04/D-05.
- **⚠ Scope note (planner): still a multi-pipeline phase.** Each shipped grid (4×4/8×8/16×16,
  +32×32 conditional) needs its own data-gen → NPE + NRE training → fresh re-pre-registered
  ship-gate. Warrants **shared-infra plan + per-grid plans** (researcher proposed: Wave-0
  co-resolution+types+registry gated on a clean `Pkg.resolve`, then per-grid ascending-risk
  plans, then API assembly registering only what passed). Size realistically.
- **⚠ Co-resolution HARD GATE (Wave-0):** adding NeuralEstimators/Flux to the root `Project.toml`
  alongside Turing/GLMakie may not `Pkg.resolve` cleanly. Must be de-risked FIRST with an actual
  resolve; if it fails, escalate to the user (candidate fixes: Turing or GLMakie → extension).
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
