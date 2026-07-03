# Phase 7: Productionization (conditional on Go) — Research

**Researched:** 2026-07-03
**Domain:** Julia package API design (type hierarchy + registry), amortized SBI promotion (NeuralEstimators/Flux), dependency co-resolution, fine-grid patch-correlation feasibility
**Confidence:** HIGH on codebase facts and the co-resolution risk (verified from Phase-1 record + Project.toml); MEDIUM on fine-grid feasibility (reasoned extrapolation from the ≥15-survivor floor + spike's validated 8×8); MEDIUM on artifact-size/compute estimates.

## Summary

Phase 7 promotes the spike's proven amortized inference (NPE posterior + NRE Bayes factor + OOD flag + SBC harness) into `src/` as the new **primary, breaking-release public API**, behind an extensible `AbstractColocResult` type hierarchy (D-02), served by a **grid-keyed estimator registry** carrying a 5-grid family (D-04), with **every grid gated by its own fresh, re-pre-registered SBC/BF/OOD confirmation run** (D-05). The spike code is adapted, not rebuilt — the NPE/NRE/OOD/harness surfaces already exist and are proven by six completed phases of working code (the single strongest form of API verification for the young NeuralEstimators v0.2.1).

Two findings dominate planning and must reach the user before any grid is committed to.

**Finding 1 — the dependency co-resolution risk is real and already documented.** Phase 1 recorded (STATE line 89) that the root package's `Turing + GLMakie(old) + GraphNeuralNetworks` tree **caps NeuralEstimators below 0.2.1 and silently downgrades it to 0.1.4, breaking the v0.2.1 API.** D-03 requires adding NeuralEstimators/Flux/JLD2 as *hard root deps* — which will re-trigger exactly this conflict. The clean, idiomatic fix follows directly from D-01: since the Turing/ADVI path is **internal-reference-only** (ship-gate + future validation), it does **not** need to be a core *runtime* dependency. Move it to a **Julia package extension (weakdep)** or the test environment, so the core resolve is just `NeuralEstimators + Flux + Images + …` — the combination the spike *proved* co-resolves cleanly. This should be the first plan (Wave 0), gated before any `src/` inference code is written.

**Finding 2 — the fine grids do NOT deliver "local localisation" as currently conceived, and calibration ≠ resolution.** The current NPE maps the *whole* G×G correlation summary → **7 global θ** (one global ρ_true). A finer grid enriches the *input* but the *output stays global* — it produces no per-region Δρ map. Genuine per-region output is precisely Phase 12's job (GP/CAR lattice). Separately, at fine grids each patch is tiny: on the spike's dominant 256² training images, a 64×64 grid gives **4×4 = 16-pixel patches**, which after the src/ `≥15-survivor` background floor (`colocalization.jl:235`) collapse to nearly all-`missing`. SBC calibration can still *pass* at fine grids (the flow learns honest, wide posteriors), but the estimator becomes **calibrated-but-uninformative** — the exact opposite of the "sharper local resolution" the D-04 rationale assumes. This must be surfaced loudly; the recommendation is to cap the *genuinely useful* shipped family and reframe fine grids honestly.

**Primary recommendation:** Wave-0 plan does dependency co-resolution (Turing→extension) + type hierarchy + registry + summary/encoder grid-parametrization, gated on a clean root `Pkg.resolve`. Then one plan per grid (data-gen → NPE → NRE → fresh re-pre-registered SBC/BF/OOD gate), sequenced by ascending risk (8×8 → 4×4 → 16×16 → 32×32 → 64×64). A final plan assembles the public API and registers only the grids that passed. Ship 4×4/8×8/16×16 as calibrated global estimators; ship 32×32 with a large-image caveat; do **not** ship 64×64 as a "local map" (global-only or defer to Phase 12).

## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01: Breaking release, amortized-only public API (Option 3).** Public surface is the amortized path. Turing/ADVI `colocalization()` + `compute_BayesFactor()` demoted to **internal reference code** — retained, not exported — serving the D-05 ship-gate and future validation. No backward-compat obligation to the pre-v2.0 API.
- **D-02: Design an extensible, maintainable type hierarchy with clear types** (first-class deliverable). Abstract `AbstractColocResult` supertype; concrete `AmortizedColocResult`; kept-internal `AdviColocResult`. A shared accessor interface (`delta_rho`/`bayes_factor`/`is_ood`/`posterior_draws`) so downstream code is type-agnostic. Leave clean extension points for Phase 13 (three-hypothesis BF) and Phase 12 (spatial map) as **new subtypes / interface extensions**, not bolted-on fields. Prefer explicit clarity over clever abstraction (user chose abstract-supertype over parametric `ColocResult{Backend}`).
- **D-03: Add NeuralEstimators, Flux, JLD2 (+ transitive spike inference deps) as hard direct deps in the root `Project.toml`.** Amortized is the primary path; must work on a clean install with no optional-load machinery. Heavier base install accepted. Spike's isolated `spike/Project.toml`/`Manifest.toml` remain untouched.
- **D-04: Ship a five-grid pre-trained family — 4×4, 8×8, 16×16, 32×32, 64×64 —** in a registry keyed by patch grid, plus a documented `train_and_register(grid)`. Registry maps requested `num_patches`/grid to its trained estimator; an unregistered grid gives a clear, actionable error pointing at `train_and_register`. Finer grids exist to enable **"local localisation"** (per-region/spatially-resolved coloc), the on-ramp to Phase 12. Summary-vector dimension is **coupled to the grid**, so each grid needs its own data-gen + NPE + NRE + training run. Storage mechanism (Artifacts lazy-download vs bundled) is research/discretion.
- **D-05: Gate ALL five grids before ship.** Every bundled estimator must pass its own independent, **fresh-seed, re-pre-registered** SBC/BF/OOD confirmation run — new locked consts + new disjoint seed, **NOT** the spike's `VAL_MASTER_SEED`. Results recorded per grid; a grid that fails its gate does not ship and is not merged as public. No provisional estimators in the public registry.

### Claude's Discretion
- Exact naming of the public entry point and accessor functions (within the D-02 hierarchy).
- Estimator storage mechanism (Artifacts vs bundled) and registry file format.
- How the internal Turing path is structured/namespaced to stay usable by the ship-gate without being exported.
- Whether Phase 7 is split into per-grid plans (strongly implied by D-04 + D-05).

### Deferred Ideas (OUT OF SCOPE)
- Feature-DAG work (Phases 8, 11–16); any new scientific capability beyond promoting what the spike validated; retraining the *summary statistic itself* (the patch-correlation summary is reused).
- **⚠ Scope flag:** Phase 7 is large — 5 grids × per-grid fresh ship-gate = five full pipelines; likely split into shared-infra + per-grid plans.
- **⚠ Feasibility risk:** fine grids (esp. 64×64) for LOCAL localisation — does calibrated local Δρ resolution survive? (This RESEARCH's headline verdict below.)
- Re-enabling the OOD posterior-predictive channel in the productionized OOD flag (iter1 finite-θ̂ guard made it viable; reported OR-fusion still runs `with_pp=false`) — hardening item.
- RxInfer independent cross-check of the ADVI reference (BACK-01) — post-Go paper nice-to-have, not a Phase-7 dependency.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| PROD-01 | Amortized inference integrated into `src/` as the primary public API. *(SC1's "coexisting so `compute_BayesFactor` keeps working unchanged" clause is RELAXED by D-01: breaking release, Turing path internal-only.)* | Promotable spike surfaces (`spike/npe/*`, `spike/validation/*`), the `AbstractColocResult` hierarchy sketch (§Architecture Patterns), the Turing→extension co-resolution fix (§Dependency Co-Resolution). Public entry consumes `_prepare_data`/`patch_summary` output, returns `AmortizedColocResult`. |
| PROD-02 | `num_patches` (patch grid) **user-definable** via an estimator registry, not ad-hoc retraining. | Grid-keyed registry + `train_and_register(grid)` design (§Registry Design), grid↔summary-dim coupling analysis (§Grid Coupling), the FEASIBILITY VERDICT gating which grids can actually ship. |
</phase_requirements>

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Public `colocalization_amortized(...)` entry | `src/` package API (library) | — | The new primary surface; consumes the existing summary contract, returns `AmortizedColocResult`. |
| Result type hierarchy + accessor interface | `src/` types module | Phase 12/13 extension subtypes | D-02 governing deliverable; dispatch target for all downstream (plotting, Phase 14 decision layer). |
| Estimator registry (grid → bundle) | `src/` registry module | Artifacts.toml (storage) | PROD-02; keyed by `num_patches`; lazy-loads shipped `.jld2`. |
| NPE posterior inference | `src/` (promoted from `spike/npe`) | NeuralEstimators/Flux (runtime dep) | Core amortized read; per-grid trained flow. |
| NRE / amortized Bayes factor | `src/` (promoted from `spike/validation/bf.jl`) | NeuralEstimators | Model-comparison ratio; per-grid trained. |
| OOD / misspecification flag | `src/` (promoted from `spike/validation/ood.jl`) | LinearAlgebra/Statistics (hand-rolled, no dep) | Paired guard for SBC-under-simulator honesty. |
| Turing/ADVI reference (`colocalization`, `compute_BayesFactor`) | **Package extension (weakdep) OR test env** | Turing (moved off core runtime) | D-01 internal-only; drives the ship-gate; kept off the core resolve to fix Finding 1. |
| Per-grid ship-gate (SBC/BF/OOD, fresh pre-reg) | `test/` or `src/validation` (promoted harness) | HypothesisTests | D-05; the phase's validation strategy. |
| Shipped estimator storage | Artifacts.toml (lazy, content-hashed) | GitHub Release assets (hosting) | D-04 discretion; large binaries out of git history. |

---

## FEASIBILITY VERDICT — Fine Grids / Local Localisation

> **This is the headline deliverable. Read before planning any grid.**

### The two things the D-04 rationale conflates

**(A) Higher-resolution INPUT ≠ per-region OUTPUT.** `[VERIFIED: spike/npe/architecture.jl + infer.jl]` The current `PosteriorEstimator` maps the **entire** G×G correlation summary vector through an MLP summary net → a `NormalisingFlow` over **7 global θ** (`NPE_D = 7`; ρ_true is θ row 1). `rho_hat` reads row 1 of the *global* posterior; `delta_rho` is the MC difference of two *global* single-stack posteriors. A finer grid makes the input vector longer (2·G² rows) but the **output is still a single global ρ_true** — there is **no per-region Δρ map** anywhere in the promoted stack. "Finer grid → per-region colocalization" is therefore a **misconception about the current architecture**: the finer grid enriches the global estimate, it does not spatialize the output. Genuine per-region output (a Δρ map + per-region uncertainty) is exactly **Phase 12's** deliverable (GP/CAR lattice prior, `coloc_map(...)` returning a map — ROADMAP §Phase 12). **This must be surfaced to the user: shipping 32×32/64×64 estimators does not, by itself, produce a spatial map.**

**(B) Calibration ≠ resolution.** SBC checks that posterior *width is honest*, not that it is *narrow*. At fine grids the per-patch correlation is estimated from very few pixels, so it is high-variance; the flow can still learn an honest (wide) posterior and **pass SBC** — while being **uninformative** about local structure. Passing the D-05 gate at a fine grid does **not** prove the fine grid is scientifically useful for local resolution.

### Pixels-per-patch budget (the hard physical constraint)

`[VERIFIED: src/colocalization.jl:37-60, 221-239]` `patch(img, G)` tiles the image into G×G patches of `(W÷G)×(H÷G)` pixels; `correlation(...)` computes per-patch Pearson **after `_exclude_zero`** (drops zero/NaN/background pixels) and sets the patch to **`missing` when ≤15 pixels survive**. So the usable signal per patch is `(W÷G)²` pixels *minus* background exclusion.

Pixels per patch `(W÷G)²` (before background exclusion):

| Grid | 256² (55% of spike data) | 512² (35%) | 1024² (5%) | 1376×1028 (real anchor) | 2048² (2%) |
|------|------|------|------|------|------|
| 4×4 | 4096 | 16384 | 65536 | 88408 | 262144 |
| 8×8 | 1024 | 4096 | 16384 | 22016 | 65536 |
| 16×16 | 256 | 1024 | 4096 | 5504 | 16384 |
| 32×32 | 64 | 256 | 1024 | 1376 | 4096 |
| 64×64 | **16** | 64 | 256 | 336 | 1024 |

`[VERIFIED: spike/data/seeding.jl:52-62]` The spike's `IMSIZE_SET` is weighted **0.55 → 256², 0.35 → 512²**, i.e. 90% of training stacks are ≤512². The `sample_imsize` docstring states the ≥64 floor exists specifically "so the **8×8** patch grid clears the ≥15-px floor." Nothing in the current training distribution supports grids finer than ~16×16.

Approximate per-patch correlation noise (Fisher-z SE ≈ `1/√(n_eff−3)`, with `n_eff` the *surviving* pixel count, realistically ~30–60% of the raw count after background exclusion):

- n≈1024 → SE(ρ)≈0.03 (tight, informative)
- n≈256 → SE(ρ)≈0.06 (usable)
- n≈64 → SE(ρ)≈0.13 (noisy)
- n≈16 → SE(ρ)≈0.29 **or `missing`** (uninformative; correlation swings ±0.3 by noise alone)

### Per-grid verdict

Assumes the planner **regenerates training data per grid at an image-size distribution adequate for that grid** (the existing 256²-heavy cache cannot train grids finer than 16×16 — see Grid Coupling). "Ship" = a calibrated **global** estimator; the spatial-map claim is addressed separately in the Phase-12 overlap note.

| Grid | Summary dim (2·G²) | Verdict | Reason |
|------|------|---------|--------|
| **4×4** (16 patches) | 32 | **SHIP** | Robust at any image size (≥4096 px/patch even at 256²). Very coarse "local" resolution but calibratable and cheap. |
| **8×8** (64 patches) | 128 | **SHIP** | The exact grid the spike validated end-to-end (Phases 4–6). Proven. The reference grid. |
| **16×16** (256 patches) | 512 | **SHIP-WITH-CAVEAT** | Needs data biased to ≥512² (px/patch 256 at 256² is marginal after background exclusion). Calibratable if data regenerated; document the minimum-image-size caveat in the registry. |
| **32×32** (1024 patches) | 2048 | **SHIP-WITH-CAVEAT / global-only** | Needs ≥1024² images (px/patch 64 at 512² is noisy; 64 at 256² collapses to `missing`). Calibration likely holds with wide posteriors; local resolution degraded. Ship only if its fresh D-05 gate passes AND registry warns "requires large images." |
| **64×64** (4096 patches) | 8192 | **DO-NOT-SHIP as a local map** (global-only or DEFER) | On the dominant 256² stack, 16 px/patch ⇒ **≤15-survivor floor ⇒ nearly all `missing`** ⇒ mask rows ~0 ⇒ posterior collapses to prior. Requires ≥2048² for genuine signal (1024 px/patch). 8192-dim summary + massive missingness + highest compute. Even trained on huge images it yields a *global* ρ, not a map. **Recommend: do not ship as "local localisation"; either ship as an experimental global estimator with a loud caveat, or defer the fine-grid ambition to Phase 12.** |

### Design overlap with Phase 12 (Spatial Colocalization Map, GP/CAR) — flag

The D-04 fine grids are **not the right primitive** for per-region resolution and **partially pre-empt Phase 12**. Phase 12 (ROADMAP §Phase 12) explicitly replaces exchangeable patch pooling with a **spatial lattice (CAR/GP) prior over the correlation grid** and a `coloc_map(...)` that returns a **Δρ map + uncertainty map** — the actual spatialized output. The Phase-7 fine grids only lengthen the *input* to a still-global estimator. Recommendation:

- **Ship 7's fine grids (if at all) as global estimators**, framed honestly as "higher-dimensional summary for large images," not as a spatial map.
- **Let the fine grids INFORM Phase 12** (they exercise the data-gen/training machinery at high input dimension), but **defer the per-region-map claim to Phase 12**, where the output architecture actually produces one.
- A cheap honest "on-ramp" that needs **no** 32×32/64×64 estimator: run the shipped **8×8** estimator on image **sub-tiles** (windowed inference) to get a coarse map. This gives a genuine (coarse) local readout without the fine-grid feasibility cliff. Consider recommending this over shipping 64×64.

### Compute cost reality (per-grid, so the planner can size the phase)

`[ASSUMED — extrapolated from the spike's 8×8 run at 50k pairs]` Each grid is a full pipeline: **data-gen → NPE train → NRE train → fresh re-pre-registered SBC/BF/OOD gate**. **Phase 7 relaxes the spike's CPU-only training rule — a GPU is available and SHOULD accelerate NPE+NRE training** (see §GPU Acceleration & Reproducibility Split). The two cost columns below separate the **image-simulation** axis (data-gen; CPU-bound, GPU does not help — it is Julia image ops, not tensor math) from the **neural training** axis (GPU-accelerable):

| Grid | Data-gen (CPU, sim-bound) | Summary dim | NPE+NRE training (CPU baseline) | NPE+NRE training (GPU-accelerated) |
|------|------|------|------|------|
| 4×4 | 50k @ 256²-heavy (~1×) | 32 | low | trivial |
| 8×8 | 50k (spike cache reusable) | 128 | baseline (already trained) | — |
| 16×16 | 50k–100k @ ≥512² (~4× sim) | 512 | ~2–3× | ~1× (GPU absorbs the bigger net) |
| 32×32 | 100k @ ≥1024² (~16× sim) | 2048 | ~4–6× | ~1.5–2× (GPU absorbs the net; **sim now dominates**) |
| 64×64 | 100k–200k @ ≥2048² (~64× sim) | 8192 | ~10–20× | ~2–4× (GPU absorbs the net; **data-gen is the bottleneck**) |

**Key GPU insight:** GPU eases the fine-grid **training-time** axis (the bigger 2048/8192-dim nets), but the **data-generation** axis (simulating 100k–200k image pairs at ≥1024²/≥2048²) is **CPU/image-ops-bound and GPU does not help it** — it becomes the fine-grid bottleneck. So GPU materially de-risks 32×32/64×64 *training feasibility* but not the *data-gen wall clock*. The **statistical** feasibility (does calibrated local Δρ resolution hold up) is entirely independent of GPU and is **unchanged** by this — the DO-NOT-SHIP-as-local-map / global-θ-output findings above stand on their own merits.

Plus the **D-05 gate per grid**: SBC at M=2000 = 2000 fresh simulate+infer passes at that grid's image size, the BF Δρ sweep, and the OOD ROC grid. The gate's *simulation* cost scales with image size (CPU-bound); its *inference* passes run on the shipped **CPU-reproducible** path (§GPU section) so the pre-registered numbers stay deterministic. **Sequencing recommendation:** do 8×8 first (re-validate the proven grid on a fresh seed to establish the D-05 machinery), then ascending compute/risk. If 32×32/64×64 remain infeasible even with GPU training — because of data-gen cost **or** the statistical feasibility cliff — that is a **user decision point** (cap the family), not a silent descope.

---

## GPU Acceleration & Reproducibility Split

> **User update (2026-07-03):** the CPU-only rule was a **spike constraint only** and does **not** bind Phase 7. A GPU is available and MAY accelerate training. CLAUDE.md stance holds: **GPU is an optional accelerator with graceful CPU fallback.** This section revises the training plan accordingly; it does **not** soften any feasibility verdict (the local-localisation calibration question is statistical, GPU-independent).

### The mechanism (NeuralEstimators v0.2.1) — verified

`[VERIFIED: spike/npe/train_npe.jl:29-33 + spike/validation/train_ratio.jl:49-51, execution-proven across Phases 4–6]` NeuralEstimators' `train`, `sampleposterior`, and `logratio` all take a **`use_gpu` keyword that DEFAULTS TO `true`** — the spike disables it (`use_gpu = false`) on *every* call precisely because the spike was CPU-only. Phase 7 simply **stops forcing `use_gpu = false` in the training path**. `[CITED: CLAUDE.md Technology Stack + fluxml.ai GPU guide]` Since Flux ≥0.14, **CUDA loads as a package extension** (`CUDA.jl` added separately; Flux picks it up automatically), so:

- **Training:** call `train(est, …; use_gpu = true)` for NPE and NRE. With `CUDA.jl` present and a device visible, tensors move to GPU; with no CUDA device, NeuralEstimators/Flux **degrade to CPU** without code change (graceful fallback — the CLAUDE.md constraint).
- **Dependency:** add `CUDA.jl` as an **optional** dependency — a **weakdep/extension of the productionized package** (mirrors the Turing→extension move), NOT a hard runtime dep. This keeps a clean-CPU install lean (a core resolve without CUDA) and honors "spike/degrade must run without GPU." A `train_and_register(grid; use_gpu = has_cuda_device())` default auto-selects.
- **Verify at Wave 0:** confirm the installed NeuralEstimators 0.2.1 honors `use_gpu=true` on a tiny GPU smoke train (mirror the ENV-02 CPU smoke), and that the CPU path is byte-reproducible with `use_gpu=false`.

### What GPU changes — and what it does NOT

| Axis | GPU effect | Verdict impact |
|------|-----------|----------------|
| NPE/NRE **training time** (heavy 2048/8192-dim nets at 32×32/64×64) | **Eased** — GPU absorbs the bigger nets (see revised compute table) | Fine-grid **training feasibility** improves; sequencing risk drops. |
| **Data generation** (100k–200k image sims at ≥1024²/≥2048²) | **No help** — CPU/image-ops bound (Julia `imfilter`/patch/correlation), not tensor math | Remains the fine-grid wall-clock bottleneck. |
| **Local-localisation calibration** (does per-region Δρ resolution hold up) | **None** — a statistical property of the summary + prior, independent of hardware | **Feasibility verdict UNCHANGED.** 64×64 stays DO-NOT-SHIP-as-a-local-map; global-θ-output finding stands. |

**Do not let GPU availability soften Finding 2.** GPU makes it *cheaper to train* a 64×64 net; it does nothing to make that net produce a per-region map (it still outputs global θ) or to make 16-pixel patches informative. Those are architecture/statistics facts.

### Reproducibility split (plan as the default, pending final user confirmation)

**Train on GPU; keep the shipped/gated/inference paths CPU-reproducible.** Concretely:

1. **Shipped frozen nets** — trained however is fastest (GPU), but **persisted as CPU-resident arrays** (move parameters to CPU before `Flux.state` + atomic save). The artifact carries no device state, so it loads identically on any machine.
2. **Per-grid ship-gate (SBC/BF/OOD, D-05)** — runs the **CPU inference path** (`use_gpu = false` on every `sampleposterior`/`logratio`), keyed by the fresh `PROD_SEED[grid]` (Random123). This is what makes the **pre-registered numbers deterministic and reproducible** regardless of training hardware. The gate is the pre-registration; the gate is CPU.
3. **Shipped default inference** (`colocalization_amortized`) — **CPU by default** (`use_gpu = false`), so end users reproduce results without a GPU. Optionally expose `use_gpu` for batch throughput, but the *documented, reproducible* default is CPU.

**Determinism caveat (state it plainly):** GPU *training* is typically **non-deterministic** (non-associative float reductions, cuDNN algorithm selection) — this is **acceptable and by design here**, because *what is pre-registered and reproduced is the **frozen net + CPU inference/gate**, not the training run*. The training seed influences which net you get; once frozen, the net's CPU inference and the CPU gate are deterministic. If two GPU training runs yield slightly different frozen nets, each must still **independently pass its own fresh CPU gate** — the D-05 contract is on the gate outcome, not on training-run bit-reproducibility. (If the user later wants bit-reproducible *training* too, that forces CPU training or a deterministic-GPU config and a large slowdown — flag as a trade, do not assume it.)

### Compute-table caveat

The revised GPU column in the feasibility §Compute cost table reflects **training** acceleration only. The data-gen column is unchanged (CPU/sim-bound). Re-estimate both against the real 8×8 GPU-train timing captured in Wave 0 / the first grid.

---

## Standard Stack

All packages are **already proven** in the spike Manifest and by six completed phases (the strongest verification for young NeuralEstimators v0.2.1). **Pin the productionized deps to the exact versions the spike validated** — the shipped estimators are tied to them; do not bump.

### Core (add to root `Project.toml` as hard deps — D-03)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| NeuralEstimators | **0.2.1** | `PosteriorEstimator`/`RatioEstimator`/`NormalisingFlow`/`train`/`sampleposterior`/`logratio`/`assess`/`interval` | Only Julia-native package delivering both amortized posterior AND amortized ratio from one API. `[VERIFIED: spike Manifest + 6 phases of working code]` `[CITED: msainsburydale.github.io/NeuralEstimators.jl/dev]` |
| Flux | **0.16.10** | NN backend for the MLP summary net + flow conditioner | Windows-tauglich; since v0.14 CUDA is a separately-loaded extension so CPU installs stay lean, GPU degrades gracefully. `[VERIFIED: spike Manifest]` |
| JLD2 | current (spike-pinned) | Persist trained estimators + frozen transforms + calibration reports | Pure-Julia atomic persistence; already the spike idiom (`save_npe`/`save_ratio`). `[VERIFIED: spike Manifest]` |
| Random123 | 1.7.1 | Counter-based reproducible seeding of data-gen + gates | Stateless/splittable; the spike's reproducibility primitive. `[VERIFIED: spike Manifest]` |
| HypothesisTests | current | KS/χ² rank-uniformity for the per-grid SBC gate | Standard; the spike's SBC uniformity test. `[VERIFIED: spike Manifest]` |
| StatsBase | current | `ZScoreTransform` fit/transform/reconstruct (frozen standardization) | Already root + spike dep. `[VERIFIED: both Project.toml]` |

### Supporting (already root deps or spike-proven)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Images / ImageFiltering | current | `MultiChannelImage` assembly + PSF ops for the simulator used by data-gen + gates | Already root + spike. |
| Distributions | current | Prior π(θ) sampling consistent with the Turing `@model` ranges | Already root + spike. |
| LinearAlgebra / Statistics | stdlib | Hand-rolled Mahalanobis OOD + ROC/AUC (no new dep) | OOD flag; no ROC package enters the resolve. |
| CUDA | current (optional) | **GPU acceleration of NPE/NRE `train(...; use_gpu=true)`** | Optional accelerator (§GPU section); add as a **weakdep/extension**, NOT a hard core dep, so clean-CPU installs stay lean and degrade gracefully. `[CITED: fluxml.ai GPU guide — CUDA is a Flux extension since v0.14]` |

### Turing — move OFF the core runtime (the co-resolution fix)
| Library | Version | Disposition | Rationale |
|---------|---------|-------------|-----------|
| Turing | current | **weakdep (package extension) OR test-only dep** — NOT a core runtime dep | D-01 makes the ADVI path internal-reference-only (ship-gate + validation). Keeping Turing off the core resolve is what lets NeuralEstimators 0.2.1 resolve (Finding 1). See §Dependency Co-Resolution. |
| GLMakie | current (root pins 0.10.5) | Consider moving plotting to a weakdep extension | Old GLMakie/Makie is part of the cap that downgrades NeuralEstimators; shrinking the plot stack from the core resolve de-risks. |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Whole-object `jldsave(estimator)` (current spike idiom) | `Flux.state(est)` + architecture rebuild + `Flux.loadmodel!` | For a *shipped* package whose estimators must survive Flux/Julia patch bumps, the Flux-recommended `state`+`loadmodel!` path is far more robust than whole-object serialization (which breaks on internal type-layout changes). `[CITED: fluxml.ai/Flux.jl/stable/guide/saving]` **Recommend migrating the registry's persistence to `Flux.state` while keeping the atomic `.tmp`+integrity+`mv` idiom.** |
| Artifacts lazy-download | bundle `.jld2` in-repo | Bundling bloats clone + git history permanently (5 grids up to 8192-dim flows). Lazy artifacts are the ecosystem standard for large binaries. `[CITED: pkgdocs.julialang.org/v1/artifacts]` |
| Turing as weakdep | Turing as hard core dep | Hard dep re-triggers the NeuralEstimators<0.2.1 cap (Finding 1). |

**Installation (root, after Turing→extension):**
```julia
# In the ROOT project — new hard deps (D-03)
Pkg.add(["NeuralEstimators", "Flux", "JLD2", "Random123", "HypothesisTests"])
# Pin to the spike-validated versions before committing Manifest:
#   NeuralEstimators 0.2.1, Flux 0.16.10
```

**Version verification:** The spike Manifest already pins NeuralEstimators 0.2.1 / Flux 0.16.10 as a reproducibility artifact (ENV-03), and these versions are proven by execution across Phases 1–6. `[VERIFIED: spike/Project.toml + STATE Phase 01-03]`. Do **not** upgrade during Phase 7 — the trained estimators are version-coupled.

## Package Legitimacy Audit

> Julia ecosystem (Julia General registry), not npm/PyPI — slopcheck (npm/PyPI-oriented) is **not applicable**. No new or hallucinated packages are introduced: every package is a pre-existing spike dependency, registered in the Julia General registry, and exercised by six completed phases.

| Package | Registry | Provenance | Disposition |
|---------|----------|------------|-------------|
| NeuralEstimators 0.2.1 | Julia General | In spike Manifest; companion CRAN R package; active repo | Approved (execution-verified) |
| Flux 0.16.10 | Julia General | In spike Manifest; FluxML org | Approved |
| JLD2 | Julia General | In spike + root; JuliaIO | Approved |
| Random123, HypothesisTests, StatsBase, Images, Distributions | Julia General | In spike + root | Approved |

**Packages removed due to slopcheck [SLOP] verdict:** none (tool N/A for Julia).
**Packages flagged [SUS]:** none.
**Supply-chain note (real risk, see Security Domain):** shipped `.jld2` estimators are **deserialized Julia objects** — loading an untrusted `.jld2` can reconstruct arbitrary types. Ship via **content-hashed Artifacts (tree-sha1 verified)** so the registry only loads integrity-checked model files.

## Architecture Patterns

### System Architecture Diagram

```
                    ┌─────────────────────────────────────────────┐
 user image stacks  │  colocalization_amortized(img, control,     │
 (MultiChannelImage │      channels; num_patches = 8)   [PUBLIC]   │
  Stack)  ──────────▶                                              │
                    └───────────────┬─────────────────────────────┘
                                    │ validate grid
                                    ▼
                    ┌───────────────────────────────┐   unregistered grid
                    │  estimator_for(num_patches)    │──────────────────────▶ ArgumentError
                    │  (grid-keyed registry)         │   "call train_and_register(grid)"
                    └───────────────┬────────────────┘
                                    │ EstimatorBundle{npe, ratio, ood_nulls, zt, θzt, calib}
                                    │  (lazy-loaded from Artifacts on first use)
                                    ▼
   _prepare_data / patch_summary(img, G)  →  encode_d01_G  →  standardize_summary(·, zt)
                                    │ (UNCHANGED src summary, grid-parametrized)
              ┌─────────────────────┼──────────────────────┬───────────────────┐
              ▼                     ▼                      ▼                   ▼
        NPE posterior         NRE log-BF            OOD flag            (Δρ = MC diff of
     sampleposterior(npe)   logratio(ratio)     Mahalanobis+noise       two NPE passes)
        → 7×N draws        → amortized log-BF     OR-fusion → is_ood
              │                     │                      │                   │
              └─────────────────────┴──────────┬───────────┴───────────────────┘
                                               ▼
                              ┌────────────────────────────────────┐
                              │  AmortizedColocResult               │
                              │   <: AbstractColocResult            │
                              │  accessor interface:                │
                              │  delta_rho / bayes_factor /         │
                              │  is_ood / posterior_draws           │
                              └────────────────────────────────────┘

  Internal reference path (NOT exported; weakdep extension / test env):
      colocalization() ADVI  →  AdviColocResult  →  used by the per-grid ship-gate
      compute_BayesFactor() (KDE)  →  BF reference for the gate
```

### Component Responsibilities
| Component | New file (suggested) | Promoted from | Responsibility |
|-----------|----------------------|---------------|----------------|
| Result types + accessors | `src/results.jl` | new (D-02) | `AbstractColocResult`, `AmortizedColocResult`, internal `AdviColocResult`, accessor interface |
| Registry | `src/registry.jl` | new (D-04) | grid → `EstimatorBundle`; `estimator_for`, `register!`, `train_and_register`, Artifacts loading |
| NPE inference | `src/amortized/infer.jl` | `spike/npe/infer.jl` | `posterior_for`, `rho_draws`, `delta_rho`, `standardize_summary` |
| NPE training | `src/amortized/train_npe.jl` | `spike/npe/train_npe.jl` | `build_estimator`, `train_fold`, `save_npe`/`load_npe` (→ `Flux.state`) |
| NRE / BF | `src/amortized/bf.jl` + `train_ratio.jl` | `spike/validation/{bf,train_ratio}.jl` | `amortized_log_bf`, ratio training |
| OOD | `src/amortized/ood.jl` | `spike/validation/ood.jl` | Mahalanobis + noise channels, OR-fusion, PP re-enable (hardening) |
| Grid-parametrized summary | edit `patch_summary` + `encode_d01` | `spike/contract.jl` + `spike/data/encode.jl` | make grid `G` a parameter (§Grid Coupling) |
| Ship-gate harness | `test/gate/` or `src/validation/` | `spike/validation/{harness,sbc}.jl` + `consts.jl` | per-grid fresh-seed SBC/BF/OOD, fresh pre-registration |
| Turing reference | `ext/ProteinCoLocTuringExt.jl` | existing `src/bayes.jl` body | ADVI `colocalization` + `compute_BayesFactor`, weakdep-gated |

### Pattern 1: Abstract result type + accessor interface (D-02)
```julia
# src/results.jl — the governing D-02 deliverable
abstract type AbstractColocResult end

# --- accessor interface (defined on the ABSTRACT type; downstream dispatches on these) ---
"""Posterior Δρ (sample − control). Concrete subtypes must implement."""
delta_rho(r::AbstractColocResult)       = _iface_error(r, :delta_rho)
bayes_factor(r::AbstractColocResult)    = _iface_error(r, :bayes_factor)
is_ood(r::AbstractColocResult)          = _iface_error(r, :is_ood)
posterior_draws(r::AbstractColocResult) = _iface_error(r, :posterior_draws)
_iface_error(r, f) = error("$(typeof(r)) does not implement the AbstractColocResult accessor `$f`")

# --- primary shipped subtype ---
struct AmortizedColocResult <: AbstractColocResult
    grid            :: Int                       # num_patches (registry key)
    posterior       :: Matrix{Float64}           # 7×N physical-θ draws (row 1 = ρ_true)
    delta_rho_draws :: Vector{Float64}           # MC difference cloud (D-03)
    log_bayes_factor:: Float64                   # amortized NRE log-BF
    ood             :: OODVerdict                 # score, flag, per-channel breakdown
    calibration     :: CalibrationMeta           # SBC/ECE provenance for THIS grid (D-05)
    meta            :: NamedTuple
end
delta_rho(r::AmortizedColocResult)       = mean(r.delta_rho_draws)
bayes_factor(r::AmortizedColocResult)    = exp(r.log_bayes_factor)
is_ood(r::AmortizedColocResult)          = r.ood.flag
posterior_draws(r::AmortizedColocResult) = r.posterior

# --- internal reference subtype (NOT exported; lives with the Turing extension) ---
struct AdviColocResult <: AbstractColocResult    # was src/bayes.jl `CoLocResult`
    img; control; channels; num_patches
    posterior :: DataFrame
    advi_result
end
delta_rho(r::AdviColocResult) = r.posterior.μ_sample .- r.posterior.μ_control
# bayes_factor(r::AdviColocResult) delegates to compute_BayesFactor (KDE reference)
```
**Extension points (sketch only — out of scope to implement, D-02):**
```julia
# Phase 13 — three-hypothesis BF: new subtype + interface extension, NOT a new field
struct ThreeHypothesisColocResult <: AbstractColocResult ... end
bayes_factor_simplex(r::ThreeHypothesisColocResult)      # {coloc, random, exclusion}
bayes_factor(r::ThreeHypothesisColocResult) = _reduce_to_twoway(r)   # honors the 2-way accessor

# Phase 12 — spatial map: new subtype + map accessors
struct SpatialColocResult <: AbstractColocResult ... end
delta_rho_map(r::SpatialColocResult)      # per-region Δρ (the actual spatial output)
uncertainty_map(r::SpatialColocResult)
delta_rho(r::SpatialColocResult) = mean(delta_rho_map(r))   # scalar summary honors the interface
```

### Pattern 2: Grid-keyed registry + `train_and_register` (PROD-02, D-04)
```julia
# src/registry.jl
struct EstimatorBundle
    grid       :: Int
    npe                                   # PosteriorEstimator (+ θzt)
    ratio                                 # RatioEstimator (+ log_prior_odds)
    ood_nulls                             # frozen Mahalanobis + noise nulls
    zt                                    # frozen summary ZScoreTransform
    θzt                                   # frozen θ ZScoreTransform
    calibration :: CalibrationMeta        # the grid's PASSED D-05 gate report
end

const _REGISTRY = Dict{Int, EstimatorBundle}()
const _SHIPPED_GRIDS = (4, 8, 16, 32, 64)   # only those that PASSED their gate get populated

function estimator_for(grid::Integer)
    haskey(_REGISTRY, grid) && return _REGISTRY[grid]
    grid in _SHIPPED_GRIDS && (return _lazy_load_from_artifact!(grid))   # Artifacts on first use
    throw(ArgumentError(
        "No estimator registered for a $(grid)×$(grid) patch grid. " *
        "Shipped grids: $(_SHIPPED_GRIDS). " *
        "Train and register a new grid with `train_and_register($grid)`."))
end

register!(b::EstimatorBundle) = (_REGISTRY[b.grid] = b)
function train_and_register(grid::Integer; kwargs...)
    b = _train_grid_pipeline(grid; kwargs...)   # data-gen → NPE → NRE → OOD nulls
    register!(b); return b
end
```

### Pattern 3: Turing/ADVI as a weakdep extension (D-01 internal-only, Finding-1 fix)
```toml
# root Project.toml
[weakdeps]
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"

[extensions]
ProteinCoLocTuringExt = "Turing"
```
```julia
# ext/ProteinCoLocTuringExt.jl — loaded ONLY when Turing is present (ship-gate/test env)
module ProteinCoLocTuringExt
using ProteinCoLoc, Turing, ...
# defines the ADVI colocalization() + compute_BayesFactor() reference path,
# returning an AdviColocResult. Not exported; the per-grid gate calls it as the
# non-amortized accuracy/BF reference.
end
```
This keeps `Turing` out of the **core** resolve (fixing Finding 1) while satisfying D-01 (the reference code is retained and usable by the gate) and D-03 (NeuralEstimators/Flux/JLD2 are the hard runtime deps).

### Anti-Patterns to Avoid
- **Do NOT `export` the Turing path or `AdviColocResult`.** D-01 makes it internal reference only. Exporting it re-creates a "supported" second public API.
- **Do NOT reuse the 8×8 net for other grids.** The summary dim is `2·G²` — each grid needs its own architecture + training (D-04). `[VERIFIED: encode.jl hardcodes 128; architecture builds from `size(Ztr,1)`]`.
- **Do NOT reuse `VAL_MASTER_SEED` or the spike's `consts.jl` for any grid's gate.** D-05 mandates fresh consts + fresh disjoint seed per grid (the net was iterated against `VAL_MASTER_SEED`).
- **Do NOT ship whole-object-serialized estimators for cross-version robustness.** Migrate to `Flux.state`+`loadmodel!`.
- **Do NOT present 32×32/64×64 as a "spatial map."** It is a higher-D global estimate (Finding 2).

## Grid Coupling (the mechanical work each grid needs)

`[VERIFIED: spike/data/encode.jl, spike/contract.jl, spike/data/generate.jl]`

| Coupled site | Current (8×8 hardcoded) | Generalize to grid `G` |
|--------------|-------------------------|------------------------|
| `patch_summary` (contract.jl:88) | `patch.([x,y], 8)` | `patch.([x,y], G)` — `patch`/`correlation` already accept any `num_patches` |
| `encode_d01` (encode.jl:71) | 128 = 64 vals + 64 mask | `2·G²` (rows 1:G² vals, G²+1:2G² mask) |
| `generate.jl` buffers | `summary_min = Matrix(128, N)` | `Matrix(2G², N)` |
| loader `_row_partition` | continuous 1:64, mask 65:128 | continuous 1:G², mask G²+1:2G² |
| NRE `pair_encode`/`RATIO_INPUT_DIM` | 320 = 128+128+64 | `2·(2G²)+G²` = `5G²` |
| NPE `build_estimator` | reads `size(Ztr,1)`; input-agnostic already | no change (data-driven) ✓ |
| cache `generating_config` | `summary_min_dim = 128` | `2G²` (invalidates cross-grid caches — good) |

**Implication:** the Wave-0 plan must parametrize the summary + encoder by `G` **once**, then each grid plan instantiates it. The content-hash cache key already includes `summary_min_dim`, so distinct grids get distinct cache dirs automatically (no collision).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Posterior sampling / density flow | Custom normalizing flow | NeuralEstimators `NormalisingFlow` (built-in) | Proven in spike; British-spelling type, distinct from `NormalizingFlows.jl` (do not add the latter). |
| Amortized Bayes factor | quadgk/KDE/shuffle | NeuralEstimators `logratio` on a `RatioEstimator` | One forward pass; the whole point of BF-01/02. |
| Per-parameter RMSE / intervals | manual error loops | NeuralEstimators `assess`/`rmse`/`interval` | Spike convention; correct scaling. |
| Estimator persistence across versions | whole-object `jldsave` | `Flux.state` + `loadmodel!` + atomic `.tmp` idiom | Forward-compatible loading of shipped artifacts. `[CITED: Flux docs]` |
| Large model-file distribution | commit `.jld2` to git | lazy content-hashed `Artifacts.toml` | Ecosystem standard; integrity + lean clones. `[CITED: Pkg artifacts docs]` |
| Reproducible seeding | `Random.seed!` | Random123 `Philox4x` keyed streams | Spike's disjoint-stream discipline (D-02/D-05). |
| SBC rank uniformity | custom stats | HypothesisTests KS/χ² | Spike's `consts.jl` gate. |
| Dependency conflict "just pin it" | ad-hoc `[compat]` bounds | early `Pkg.resolve` spike + Turing→extension | Finding 1 is structural, not a compat-bound tweak. |

**Key insight:** Almost nothing needs building — the spike already built and validated the inference machinery. Phase 7's real work is **API design (types + registry), dependency surgery (Turing→extension), grid generalization, and 5 fresh gates.** The risk is in resolve + fine-grid feasibility, not in new algorithms.

## Runtime State Inventory

> Rename/refactor-adjacent (promotion from `spike/` to `src/` + a breaking public API). Runtime/state considerations:

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | Trained spike artifacts `spike/npe/trained_npe.jld2`, `spike/validation/trained_ratio.jld2` are the **8×8** estimators. They are frozen at the spike env's Flux 0.16.10 layout. | For the shipped 8×8 grid, **re-persist via `Flux.state`** into the registry's artifact format; do not just copy the whole-object `.jld2`. Other grids are trained fresh. |
| Live service config | None — no external services, DBs, or daemons. | None — verified: the project is a Julia library, no runtime services. |
| OS-registered state | None. | None — verified: no scheduled tasks / daemons. |
| Secrets/env vars | None. `NPE_MASTER_SEED`/`VAL_MASTER_SEED` are config constants, not secrets (per file comments). | None. New per-grid gate seeds are new committed consts (D-05), not secrets. |
| Build artifacts / installed packages | Root `Manifest.toml` will change when NeuralEstimators/Flux/JLD2 are added; `version = 1.0.1` → **bump to 2.0.0** (breaking release, D-01). Spike Manifest stays untouched (D-03). | Update root `Project.toml` version to 2.0.0; regenerate root `Manifest.toml` via the resolve spike; leave `spike/*` byte-identical. |

**Decoupling note:** D-03 requires the spike env stay untouched. Phase 7 edits **only** the root package. The promotion is a **copy-and-adapt** (spike files remain as the frozen reproducibility record), not a move.

## Common Pitfalls

### Pitfall 1: Adding NeuralEstimators to root silently downgrades it to 0.1.4
**What goes wrong:** `Pkg.add("NeuralEstimators")` into the current root (Turing + GLMakie 0.10.5) resolves to **0.1.4**, whose API differs from the 0.2.1 the spike/estimators require — everything "works" but against the wrong API.
**Why it happens:** `[VERIFIED: STATE Phase 01-04]` The parent's Turing/GLMakie/GraphNeuralNetworks tree caps NeuralEstimators <0.2.1.
**How to avoid:** Wave-0 resolve spike **before** any src code: move Turing to a weakdep extension (and consider GLMakie→extension), then assert the resolved Manifest pins NeuralEstimators **0.2.1** and Flux **0.16.10**. Automate the assertion (mirror the spike's resolve-risk gate).
**Warning signs:** `Pkg.status` shows NeuralEstimators 0.1.4; `logratio`/`NormalisingFlow` kwargs error.

### Pitfall 2: Training a grid on the 256²-heavy cache
**What goes wrong:** A 32×32 or 64×64 net trained on the spike's `IMSIZE_SET` (90% ≤512²) sees mostly-`missing` summaries → learns the prior → passes SBC trivially but is useless.
**Why it happens:** `[VERIFIED: seeding.jl:52-62 + colocalization.jl:235]` ≥15-survivor floor vs 16-px patches.
**How to avoid:** Each fine grid regenerates data with an image-size distribution adequate for its patch size (see Feasibility table). The cache content-hash key already forces a fresh cache per summary-dim.
**Warning signs:** High fraction-missing in `encode_aug`'s `frac_missing`; posterior width ≈ prior width; ECE green but interval width ≈ prior.

### Pitfall 3: Reusing the spike's pre-registration for a grid's gate
**What goes wrong:** Scoring a grid against `VAL_MASTER_SEED` / the spike `consts.jl` is data-snooping — the net was iterated against that stream (memo §4).
**Why it happens:** Copy-paste of `spike/validation/consts.jl`.
**How to avoid:** Per grid, commit a **new** consts file (M/L/bins/thresholds) and a **new disjoint** seed **before** the gate runs (D-05). Recommend a `PROD_SEED[grid]` family provably disjoint from `VAL_MASTER_SEED`, `NPE_MASTER_SEED`, and each other.
**Warning signs:** Gate imports `spike/validation/consts.jl`; any grid's seed equals `0x5BC0FFEE`.

### Pitfall 4: Whole-object estimator serialization in the shipped registry
**What goes wrong:** A shipped `.jld2` holding a whole `PosteriorEstimator` object fails to load after a Flux patch bump (internal struct layout changed).
**How to avoid:** Persist `Flux.state(est)` + enough metadata to rebuild the architecture, load with `loadmodel!`. Keep the atomic `.tmp`+integrity+`mv` idiom.

### Pitfall 5: Exporting the Turing path "for convenience"
**What goes wrong:** Users depend on `colocalization()`/`compute_BayesFactor()`, re-creating the maintenance burden D-01 explicitly retired.
**How to avoid:** Keep them unexported, in the weakdep extension, documented as "internal reference for the ship-gate."

## Code Examples

### Grid-parametrized summary (the one shared edit)
```julia
# promoted patch_summary, now grid-parametrized (was hardcoded `8`)
function patch_summary(mci::MultiChannelImage, G::Integer)
    x, y = mci.data[1], mci.data[2]
    xp, yp = patch.([x, y], G)                     # src/colocalization.jl UNCHANGED, G-parametric
    return correlation(xp, yp; method = :pearson)  # G×G Union{Float64,Missing}
end

function encode_d01(M::AbstractMatrix)             # dim = 2·G²
    vals = vec(coalesce.(M, 0.0))                  # G² imputed correlations (column-major)
    mask = vec(Float64.(.!ismissing.(M)))          # G² present/absent mask
    return vcat(vals, mask)
end
```

### Public entry point (shape — name is discretion)
```julia
"""
    colocalization_amortized(img, control, channels; num_patches=8, N=2000) -> AmortizedColocResult

Amortized colocalization inference in a forward pass. Looks up the trained estimator
family for `num_patches` in the registry (errors with a `train_and_register` hint if the
grid is unregistered), computes the frozen summary, and returns posterior draws, the
amortized log Bayes factor, and the OOD flag as an `AmortizedColocResult`.
"""
function colocalization_amortized(img, control, channels; num_patches::Integer = 8, N = 2000)
    b   = estimator_for(num_patches)                       # registry (PROD-02)
    Zs  = standardize_summary(encode_d01(patch_summary(img,     num_patches)), b.zt)
    Zc  = standardize_summary(encode_d01(patch_summary(control, num_patches)), b.zt)
    post = posterior_for(b.npe, Zs; N = N)                 # 7×N standardized → reconstruct
    Δρ   = rho_draws(b.npe, Zs, b.θzt; N=N) .- rho_draws(b.npe, Zc, b.θzt; N=N)
    logbf = amortized_log_bf(b.ratio, pair_encode(Zs, Zc), b.ratio_log_prior_odds)
    ood   = ood_verdict(b.ood_nulls, Zs, Zc)
    return AmortizedColocResult(num_patches, reconstruct(b.θzt, post), Δρ, logbf, ood,
                                b.calibration, (; N))
end
```

## State of the Art

| Old Approach | Current Approach | When | Impact |
|--------------|------------------|------|--------|
| Per-dataset ADVI (`colocalization()`, minutes/pair) | Amortized NPE forward pass (ms) | v2.0 | The core thesis; ADVI demoted to internal reference (D-01). |
| KDE+quadgk Bayes factor | NRE `logratio` amortized log-BF | v2.0 | `compute_BayesFactor` internal-only; Phase 13 extends to 3-way. |
| Whole-object Flux serialization | `Flux.state`+`loadmodel!` | Flux ≥0.14 guidance | Robust cross-version artifact loading. `[CITED: Flux docs]` |
| Fixed 8×8 grid | Grid-keyed registry + `train_and_register` | Phase 7 (PROD-02) | User-definable `num_patches`; feasibility-capped. |

**Deprecated/outdated:**
- Root `version = 1.0.1` → bump to **2.0.0** (breaking, D-01).
- Root `GLMakie = "0.10.5"` compat pin — part of the resolve cap; revisit (extension or bump).

## Validation Architecture

> `workflow.nyquist_validation = true` → this section is authoritative. The **per-grid ship-gate IS the validation strategy** (D-05).

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (spike precedent: `runtests.jl`, `test_*.jl` fixtures vs `run_*.jl` reported) |
| Config file | root `test/runtests.jl` (new); per-grid gate under `test/gate/` |
| Quick run command | `julia --project -e 'using Pkg; Pkg.test()'` (fixture-scale, fast) |
| Full suite command | per-grid reported gate: `julia --project test/gate/run_gate.jl --grid G` (M=2000, full sweeps) |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| PROD-01 | `colocalization_amortized` returns `AmortizedColocResult`; accessors dispatch | unit | `Pkg.test()` (result-type + accessor tests) | ❌ Wave 0 |
| PROD-01 | Root resolves with NeuralEstimators 0.2.1 (no downgrade) | integration | resolve-risk assertion in `runtests.jl` | ❌ Wave 0 |
| PROD-02 | Registry returns bundle for shipped grid; errors w/ hint for unregistered | unit | `Pkg.test()` (registry tests) | ❌ Wave 0 |
| PROD-02 (per grid) | Fresh re-pre-registered **SBC** rank-uniformity (KS/χ²/ECE) | reported gate | `run_gate.jl --grid G --sbc` | ❌ per-grid |
| PROD-02 (per grid) | Fresh **BF** reproduction vs (non-clamped) KDE reference | reported gate | `run_gate.jl --grid G --bf` | ❌ per-grid |
| PROD-02 (per grid) | Fresh **OOD** ROC (incl. re-enabled PP channel) | reported gate | `run_gate.jl --grid G --ood` | ❌ per-grid |

### Sampling Rate
- **Per task commit:** `Pkg.test()` fixture-scale (types, registry, resolve, tiny-M SBC smoke).
- **Per wave merge:** the merged grid's full `run_gate.jl --grid G` (reported scale).
- **Phase gate:** every **shipped** grid green on its fresh-seed gate; unregistered/failed grids excluded from `_SHIPPED_GRIDS`.

### CPU-reproducibility of the ship-gate (mandatory)
`[per §GPU Acceleration & Reproducibility Split]` Training may use the GPU, but the **D-05 ship-gate and the shipped default inference are CPU-reproducible**: every gate `sampleposterior`/`logratio` call passes `use_gpu = false`, keyed by the fresh `PROD_SEED[grid]` (Random123), against a **CPU-resident frozen net** (parameters moved to CPU before persistence). The pre-registered numbers are therefore deterministic and reproducible regardless of training hardware. GPU-training non-determinism is acceptable because the pre-registered/reproduced object is the *frozen net + CPU gate*, not the training run; a net that trains slightly differently must still independently pass its own fresh CPU gate.

### Fresh pre-registration per grid (D-05)
- One committed `gate_consts_G.jl` per grid **before** its gate runs: `SBC_M`, `SBC_L`, `SBC_BINS`, KS/χ²/ECE thresholds, BF corr/tol + sweep, OOD quantile/AUC — mirroring `spike/validation/consts.jl` but with **new values re-expressed fresh** and a **new disjoint `PROD_SEED[G]`** (CPU-keyed).
- Falsification/pass criteria inherit the memo §6 conditions (ECE green, BF mid-range agreement, OOD family detectable by a summary-orthogonal channel).
- **Ride-along hardening (memo §5):** re-enable the OOD posterior-predictive channel (iter1 finite-θ̂ guard) in the reported OR-fusion; add a **non-clamped BF baseline** so `max|Δ logBF|` is testable without the KDE `1e-8` floor artifact.

### Wave 0 Gaps
- [ ] `test/runtests.jl` — root test entry (resolve-risk + type/registry units)
- [ ] `test/gate/run_gate.jl` + `gate_consts_G.jl` template — per-grid fresh **CPU** gate
- [ ] Promote `spike/validation/{harness,sbc,bf,ood}.jl` → `src/validation` or `test/gate` (grid-parametrized, `use_gpu=false` on gate/inference)
- [ ] GPU-train smoke (mirror ENV-02 CPU smoke): assert `train(...; use_gpu=true)` works with CUDA present AND degrades to CPU when absent; assert frozen net persists CPU-resident
- [ ] Framework: no new install — stdlib `Test` + HypothesisTests (already a dep); `CUDA.jl` as optional weakdep for training only

## Security Domain

> `security_enforcement` absent = enabled. This is a scientific Julia library with no auth/web/user-account surface, so classic ASVS categories are largely N/A; the live risks are **supply-chain + artifact integrity + numerical robustness**.

### Applicable ASVS Categories
| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | No accounts/sessions. |
| V3 Session Management | no | N/A. |
| V4 Access Control | no | Library, no multi-tenant boundary. |
| V5 Input Validation | **yes** | Validate `num_patches` against the registry (clear error, no silent fallback); validate `channels`/image shape as the existing `_prepare_data` already does. |
| V6 Cryptography | no (but integrity yes) | No crypto; use **content-hash (tree-sha1) Artifacts** for model-file integrity. |

### Known Threat Patterns for this stack
| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Untrusted `.jld2`/Flux deserialization reconstructs arbitrary Julia types | Tampering / Elevation | Ship estimators **only** via content-hashed lazy Artifacts (tree-sha1 verified before load); never load user-supplied `.jld2` in the public path without hash verification. |
| Dependency confusion / silent version downgrade (NeuralEstimators 0.1.4) | Tampering | Pin exact versions in the root Manifest; automated resolve-risk assertion (Pitfall 1). |
| Heavy transitive surface (Flux/Zygote/ChainRules) enlarges attack/breakage surface | — | Accepted tradeoff (D-03); mitigated by pinned Manifest + reproducible resolve. |
| Non-finite θ̂ / NaN correlations at fine grids crash inference | DoS (robustness) | Reuse the spike's finite-guards (`induced_mu` NaN-not-throw; OOD `_theta_tuple` finite fallback); validate finiteness at the public boundary. |

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Fine-grid per-patch correlation noise (Fisher-z SE) and the ≤15-survivor collapse make 32×32/64×64 calibrated-but-uninformative on typical images | FEASIBILITY VERDICT | If wrong (e.g. real images are large + dense), 32×32/64×64 could be more useful than stated; verify by a data-gen + tiny-SBC probe per fine grid before committing its full pipeline. |
| A2 | Fine grids require regenerating training data at larger image sizes; the 256²-heavy spike cache cannot train >16×16 | Feasibility / Pitfall 2 | If wrong, fine grids could reuse existing cache (cheaper). Low risk — grounded in the ≥15-floor code. |
| A3 | Moving Turing to a weakdep extension resolves the NeuralEstimators<0.2.1 cap | Dependency Co-Resolution / Finding 1 | If Flux/NeuralEstimators still conflict with another root dep (Images/GLMakie transitive), the resolve spike surfaces it early; worst case GLMakie also → extension, or a harder blocker escalated to user. **Confirm via the Wave-0 resolve spike before any src code.** |
| A4 | Per-grid compute/artifact-size estimates (10–20× for 64×64; 100s of MB total) | Compute cost / Storage | Rough; the Wave-0/8×8 timing calibrates the real numbers. Planner should re-estimate after the first grid. |
| A5 | Julia General registry + execution-proof suffices for package legitimacy (slopcheck N/A for Julia) | Package Legitimacy Audit | Low — all packages are pre-existing spike deps. |
| A6 | Migrating to `Flux.state`+`loadmodel!` is needed for robust shipped-artifact loading | Standard Stack / Pitfall 4 | If the package pins an exact Manifest forever, whole-object jldsave *might* suffice; but a shipped library should not assume that. |
| A7 | NeuralEstimators 0.2.1 `train(...; use_gpu=true)` + CUDA-as-Flux-extension gives GPU training with graceful CPU fallback; frozen-net CPU inference/gate stays deterministic | GPU Acceleration & Reproducibility Split | Mechanism is execution-verified from spike code (`use_gpu` defaults true); the *GPU* path itself is unexercised in the spike (CPU-only) — confirm via a Wave-0 GPU-train smoke. Reproducibility split assumes CPU inference/gate is what's pre-registered (design choice, pending final user confirmation). |

> **Note:** A3 (co-resolution) is the single assumption most likely to reshape the plan and must be de-risked **first** via an actual `Pkg.resolve` spike. A1 (fine-grid feasibility) is the one most likely to change what ships.

## Open Questions (RESOLVED)

1. **Does the Turing→extension move fully clear the root resolve?**
   - **Resolution:** decided by the Wave-0 hard resolve gate (07-00 Task 3) — `Pkg.resolve` must pin NeuralEstimators 0.2.1 or execution halts and escalates (GLMakie→extension is the standby fix).
   - Known: Turing + old GLMakie caps NeuralEstimators (Phase-1 verified). Spike proved NeuralEstimators 0.2.1 co-resolves with Images/ImageFiltering/HypothesisTests/StatsBase/CairoMakie.
   - Unclear: whether GLMakie (0.10.5) or another root transitive still conflicts once Turing is out.
   - Recommendation: **Wave-0 resolve spike is a hard gate.** If it fails, escalate to the user (candidate: GLMakie→extension too, or a Makie bump) before writing inference code.

2. **How many grids actually ship?**
   - **Resolution:** revised D-04 capped family — ship 4/8/16 (07-05/06/07), 32×32 conditional on its gate + user sign-off (07-09), 64×64 DROPPED; only gate-PASSED grids populate `_SHIPPED_GRIDS` in 07-10.
   - Known: 8×8 proven; 4×4 robust; 16×16 workable with ≥512² data.
   - Unclear: whether 32×32/64×64 pass a fresh gate at a CPU-tractable cost and whether they are scientifically worth shipping (Finding 2).
   - Recommendation: gate them; if infeasible or uninformative, **cap the shipped family with the user's sign-off** rather than shipping weak estimators. Reframe fine-grid "local localisation" as Phase-12 work.

3. **Public entry-point + accessor naming (discretion).**
   - **Resolution:** `colocalization_amortized` (07-10) with the D-02 accessors `delta_rho`/`bayes_factor`/`is_ood`/`posterior_draws` on `AbstractColocResult` (07-00).
   - Recommendation: `colocalization_amortized` (matches ROADMAP SC1 wording) with accessors `delta_rho`/`bayes_factor`/`is_ood`/`posterior_draws` exactly as D-02 names them.

4. **Should the shipped 8×8 estimator be re-trained or re-persisted?**
   - **Resolution:** re-persisted via `Flux.state` and re-gated on a fresh `PROD_SEED[8]` in 07-05 (the spike whole-object `.jld2` is not copied; the spike numbers were iterated against VAL_MASTER_SEED).
   - The spike's `trained_npe.jld2`/`trained_ratio.jld2` are whole-object at Flux 0.16.10. Re-persist via `Flux.state`; the D-05 gate must run **fresh-seed** regardless (the spike's numbers were iterated against `VAL_MASTER_SEED`).

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Julia | everything | ✓ (spike pins) | 1.12.6 | — |
| NeuralEstimators | amortized core | ✓ (spike Manifest) | 0.2.1 | none — hard requirement |
| Flux | NN backend | ✓ | 0.16.10 | Lux (reserve, not needed) |
| JLD2/Random123/HypothesisTests/StatsBase/Images/Distributions | pipeline + gate | ✓ | spike-pinned | — |
| Turing | internal ADVI reference (ship-gate) | ✓ (root dep today) | current | move to weakdep/test env |
| CUDA | **GPU-accelerated NPE/NRE training** (Phase 7, `train(...; use_gpu=true)`) | optional (GPU available per user) | — | graceful CPU fallback — training, gate, and shipped inference all run without a GPU |
| slopcheck | package legitimacy | N/A (Julia, not npm/PyPI) | — | Julia General registry + execution proof |

**Missing dependencies with no fallback:** none — all inference deps are spike-proven.
**Blocking concern (not a missing tool):** the **root co-resolution** of these deps with Turing/GLMakie (Finding 1) — an environment-configuration risk, resolved by the Turing→extension move + Wave-0 resolve gate.

## Recommended Plan Decomposition (Claude's discretion, strongly implied by D-04+D-05)

| Plan | Wave | Scope | Gate |
|------|------|-------|------|
| **7-00 Shared infrastructure** | 0 | Root dep add + **Turing→weakdep extension** + **CUDA→optional weakdep (GPU training)** + `Pkg.resolve` risk gate; version→2.0.0; `AbstractColocResult` hierarchy + accessors (D-02); registry skeleton + `train_and_register(grid; use_gpu=has_cuda_device())` + unregistered-grid error (D-04); **grid-parametrize** `patch_summary`/`encode_d01`/loader/NRE dims; per-grid **CPU** gate harness + `gate_consts` template (D-05); Flux.state persistence (CPU-resident); OOD-PP re-enable + non-clamped BF baseline (memo §5). | Root resolves w/ NeuralEstimators 0.2.1; GPU-train smoke passes AND degrades to CPU; types+registry unit tests green; 8×8 round-trips through the registry. |
| **7-01 … 7-05 Per-grid pipelines** | 1..n | One grid each: data-gen (adequate imsize, CPU/sim-bound) → **GPU-accelerated** NPE train → NRE train → OOD nulls → **fresh re-pre-registered CPU SBC/BF/OOD gate** → register if pass. Sequence **8×8 first** (establish gate machinery on the proven grid), then 4×4, 16×16, 32×32, 64×64 (ascending risk/compute). GPU eases training for 32×32/64×64; data-gen stays the fine-grid bottleneck. | Grid's fresh CPU gate green (or excluded from `_SHIPPED_GRIDS` with recorded reason). |
| **7-06 Public API + release** | final | `colocalization_amortized` public entry; finalize accessor interface + docs; populate registry from Artifacts; write `Artifacts.toml` (lazy, content-hashed, GitHub-Release-hosted); integration test; register only passed grids. | PROD-01/02 acceptance; only gated grids shipped. |

**Sequencing note:** 32×32/64×64 are **conditional** — their inclusion depends on their gate outcome and a user decision on the feasibility caveats. Treat "cap the family" as an explicit, expected branch, not a failure.

## Sources

### Primary (HIGH confidence)
- Codebase (execution-verified across Phases 1–6): `spike/npe/{architecture,infer,train_npe}.jl`, `spike/validation/{bf,train_ratio,ood,harness,consts}.jl`, `spike/data/{generate,encode,loader,seeding}.jl`, `spike/contract.jl`, `src/{ProteinCoLoc,bayes,colocalization}.jl`, root + spike `Project.toml`.
- `.planning/STATE.md` (Phase 01-04 co-resolution record — the Finding-1 anchor; Phase 5 iter1/iter2; Phase 6 Clean Go).
- `.planning/phases/06-*/06-GO-NO-GO-MEMO.md` + `06-CONTEXT.md` (ship-gate definition, D-05).
- `.planning/phases/07-*/07-CONTEXT.md` (D-01..D-05).
- `.planning/ROADMAP.md` §Phase 7/12/13.
- `[CITED]` fluxml.ai/Flux.jl/stable/guide/saving — `Flux.state`+`loadmodel!` recommended persistence.
- `[CITED]` pkgdocs.julialang.org/v1/artifacts + docs.julialang.org Artifacts — lazy content-hashed artifacts for large binaries.

### Secondary (MEDIUM confidence)
- CLAUDE.md Technology Stack table (NeuralEstimators 0.2.1 / Flux 0.16.x sourcing) — cross-checked against the spike Manifest.
- WebSearch (Flux saving discourse; Julia artifacts best-practice) — verified against official docs above.

### Tertiary (LOW confidence)
- Per-grid compute/time and total artifact-size estimates (A4) — extrapolated from the 8×8 spike; to be calibrated by the first grid.

## Metadata

**Confidence breakdown:**
- Standard stack / API: HIGH — proven by six completed phases + spike Manifest; versions pinned.
- Dependency co-resolution (Finding 1) + fix: HIGH on the risk (Phase-1 verified), MEDIUM on the exact fix clearing fully (needs the Wave-0 resolve spike).
- Type hierarchy / registry design: HIGH — grounded in existing structs + D-02/D-04.
- Fine-grid feasibility (Finding 2): MEDIUM — reasoned from verified code (≥15 floor, global-output architecture, IMSIZE distribution) but not empirically run at fine grids.
- Compute/artifact estimates: LOW-MEDIUM — extrapolated.

**Research date:** 2026-07-03
**Valid until:** ~2026-08-03 (stable Julia stack; re-check only if NeuralEstimators/Flux are bumped — which Phase 7 should NOT do).
