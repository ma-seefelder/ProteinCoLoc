# Phase 12: Spatial Colocalization Map (GP/CAR) — Research

**Researched:** 2026-07-27
**Domain:** Lattice (CAR/GP) priors as *simulator* priors; amortized hierarchical NPE with a CNN
summary over a G×G correlation grid; leave-region-out predictive coverage
**Confidence:** HIGH on the API and the simulator arithmetic (measured this session against the
pinned v0.2.1 source and the live spike environment); MEDIUM on compute estimates (single-shot
timings on a machine concurrently running Phase 13); MEDIUM on the D-09 gate design (novel
construction, no prior art in this repo).

---

<user_constraints>
## User Constraints (from 12-CONTEXT.md)

### Locked Decisions

Copied verbatim from `12-CONTEXT.md` `<decisions>`. All thirteen stand after the 2026-07-27
premise audit.

- **D-01: Hierarchical head — global ρ + low-rank per-region deviation field.** The global term stays
  exactly the quantity the shipped estimator reports, so continuity is preserved; calibration stays
  tractable (SBC on the global term plus a low-dimensional deviation summary, **not** 64 independent
  marginals); and shrinkage of deviations toward zero **is** the borrow-strength mechanism the CAR/GP
  prior is meant to supply.
  Rejected: a joint flow over all 64 lattice dims (the SBC apparatus already strained at 8 parameters
  with 5 vacuous); a per-region conditional over region index (yields no joint posterior — that *is*
  the independent-pooling baseline SC3 requires the model to beat).

- **D-02: Reshape the summary to G×G×2 and read it with a CNN. Summary CONTENT is unchanged.**
  The summary is `2·G² = 128` dims at G=8 (64 continuous + 64 mask rows, **one entry per patch**), so
  per-region information is **already present** — flattening discards the spatial *inductive bias*,
  not the spatial *information*. SC1's CNN is an architecture choice, not an information necessity.
  Consequence: `patch_summary`, the OOD density Mahalanobis reference, `local_coloc_map`, the shipped
  artifact and `DATAGEN_HASH_SRC_FILES` all stay valid. Only the network's *view* of the same numbers
  changes. Comparability claims against the shipped MLP conditioner still need care, since the input
  contract differs.

- **D-03: Run the full CAR-vs-GP mini-spike, SC1 verbatim.** Chosen over a reduced spike specifically
  to keep a pre-registered criterion intact given Phase 7's amendment history.

- **D-04: The mini-spike runs with a deliberately HIGH-RANK deviation head; production head rank is
  chosen afterwards.** A low-rank head bounds how *rough* a field it can represent, so a full-rank CAR
  field would be systematically under-fit by it — a CAR loss would then be a **head artifact reported
  as a prior verdict**. Cost accepted: the spike's head is not the production head.

- **D-05: Marginal-preserving copula construction.** Draw a correlated standard-Gaussian field from
  the CAR/GP lattice prior → push each cell through Φ → then the `MU_PRIOR` quantile function → apply
  `ghat` **elementwise** to obtain the ρ-field. Under the copula, **every region's marginal μ is
  exactly `MU_PRIOR`**, so SIM-02 holds per region and the lattice prior contributes **only spatial
  dependence**. Rejected: a field directly on ρ; global μ* plus a zero-mean deviation field.

- **D-06: Apply the field by smooth interpolation to pixel resolution** (bilinear/bicubic upsample of
  the G×G field), not piecewise-constant per cell. Rationale is an artifact risk: piecewise-constant
  cells create hard discontinuities aligned **exactly** with the patch grid the summary reads.

- **D-07: Per-region ground truth = the drawn G×G lattice values.** Interpolation (D-06) is a
  rendering detail with no bearing on scoring. Making the scored quantity anything else inserts a
  deterministic transform between what was sampled and what is scored. Rejected: the cell mean of the
  interpolated field; the induced per-cell patch correlation.

- **D-08: The spatial correlation length is INFERRED as part of θ,** not fixed and not conditioned on.
  If smoothness is fixed wrongly, the per-region uncertainty is wrong in a way that
  matched-simulation calibration **cannot** reveal. Not conditioned on because the spatial correlation
  length of biological colocalization is not something a caller can supply. Expect this column may
  prove unidentifiable; report it as such if so.

- **D-09: SC3 becomes leave-region-out PREDICTIVE coverage.** Hold out lattice regions, predict them
  from their neighbours, and check the predictive interval against the held-out region's observed
  summary. **No ground truth is required anywhere.** SC3's wording is amended — pre-declare the
  amendment before results.

- **D-10: The "independent pooling" baseline is a MATCHED ABLATION** — the same Phase-12 network, same
  training, same summary, with the spatial prior neutralized (correlation length driven to zero).
  Rejected as the *gate*: per-tile independent NPE reads. `LocalColocMap` cannot serve as a coverage
  baseline as-is.

- **D-11 (AMENDED 2026-07-25): Score on the committed real microscopy TIFFs in `test/test_images/`,
  NOT on the corpus physical anchors.** Both anchor rows are `split = sealed_holdout`,
  `sha256 = PENDING-FETCH`, reserved for the Phase-16 blind evaluation. Use
  `test/test_images/{positive,negative}/*_c{1,2,3}.tif` (1028×1376). Coverage is scored **per region**,
  so at G=8 each image contributes ~64 regions. **Do NOT touch `corpus/` sealed_holdout rows.**

- **D-12: Two-stage pre-registered trigger.** Stage 1 (cheap, before the training spend): at the
  CAR-vs-GP mini-spike, require the better prior to show a stated coverage improvement over the
  neutralized-prior ablation on simulated data. Stage 2 (after training): the per-region calibration
  gate must pass within one documented iteration.

- **D-13: On descope, ship the ABLATION MODEL** — a per-region Δρ map **with** per-region uncertainty,
  deferring only the spatial borrowing to v2.1. Caveat: if the Stage-1 trigger fires, the ablation
  exists only at spike scale.

### Carried forward (not re-asked)

- **Research-lane posture.** Train in `spike/` on a fresh DEV seed asserted disjoint from
  `PROD_SEED_V2`, `VAL_MASTER_SEED` (`0x5BC0FFEE`), `NPE_MASTER_SEED` (`0xC0FFEE`),
  `RATIO_PAIR_SEED` (`0x00000000004A7107`), `VAL_FIX_SEED` (`0xF1F7ED`), `CORPUS_MASTER_SEED`
  (`0x0000000000c05eed`) and prior dev seeds; shipped artifact and GO untouched.
  **Flagged in CONTEXT as an assumption** — see §"Research-lane posture: the evidence" below, where
  this research finds the assumption **supported**.
- **Pre-registration discipline** (Phase 11 D-04): lock thresholds in a Phase-12 consts file with the
  fresh DEV seed before any run; one documented iteration allowed.
- **Phase 7 D-02:** new result variants slot in as new `AbstractColocResult` subtypes, never as
  bolted-on fields.

### Claude's Discretion

**Nothing was delegated to Claude's discretion in this discussion. Every gray area was decided by
the user.** Every recommendation below that is not a direct consequence of a locked decision is
therefore flagged as an **open question for the user**, not as a settled choice.

### Deferred Ideas (OUT OF SCOPE)

- A new CNN summary over the raw `correlation()` grid, unconstrained by the current 128-dim encoding.
- A DeepSet summary (named in SC1, incoherent with the phase goal).
- Joint flow over all 64 lattice dimensions.
- Per-region conditional head over region index.
- Per-tile independent NPE reads as a product-level baseline (reportable, never the gate).
- Semi-synthetic imposed-field series on real images for *parameter* coverage.
- CBS Red-Green benchmark as a graded supplementary series.
- Matching a head to each prior in the mini-spike (full-rank for CAR, low-rank for GP).
- Productionizing the spatial map.

### VOID (from the 2026-07-27 premise audit) — do not plan around

- **V-1** "The spatial map trains on the registration-aware θ." No such model exists, and Δρ RMSE is
  flat in λ (ratio 1.0003).
- **V-2** "Sequencing matters: Phase 11 lands first." Phase 11 closed 2026-07-27.
- **V-3** "File overlap … serialized to avoid a merge collision." No concurrent Phase-11 writer.

</user_constraints>

---

## Phase Requirements

`REQUIREMENTS.md` maps **no requirement IDs to Phase 12** (`Requirements: TBD` in ROADMAP;
`phase_req_ids` is null). The traceability table lists IDs only through CMP-09 (Phase 9), and states
"other Phase 8/10–16 IDs still TBD". [VERIFIED: `.planning/REQUIREMENTS.md:156`]

**Planner action required:** either (a) mint `SPAT-01…SPAT-0n` IDs in `REQUIREMENTS.md` covering the
three amended success criteria, or (b) explicitly record "no new REQ-IDs" as a ruling — which is the
precedent Phase 11 set ("no new REQ-IDs" is listed among its eight orchestrator rulings,
`.planning/STATE.md`). Do not leave it implicit.

---

## Summary

Phase 12's premise survives the Phase-11 negative result intact, but this research changes three
things a planner would otherwise get wrong.

**First, the stack is smaller than SC1 implies.** At G=8 the lattice has 64 cells, so both candidate
priors are 64×64 dense linear algebra: CAR is `Q = D − αW` inverted and Cholesky-factored; GP is a
Matérn/SE correlation matrix Cholesky-factored. Both are ~15 lines of `LinearAlgebra` and cost
**40.4 µs per draw** (measured). `AbstractGPs.jl` / `KernelFunctions.jl` / `GaussianRandomFields.jl`
are all unnecessary, and adding any of them forces a `Pkg.resolve()` on `spike/Manifest.toml`, which
is the one thing that could dislodge the `NeuralEstimators 0.2.1` / `Flux 0.16.10` pin the whole
milestone rests on. **Recommendation: zero new dependencies. `LinearAlgebra` and `SparseArrays` are
stdlibs, loadable via `@stdlib` on `LOAD_PATH` without touching `spike/Project.toml`** — the existing
`run_p11_recovery.jl` already does exactly this with `LinearAlgebra`/`Statistics`/`Dates`.

**Second, the copula (D-05) works, and the CAR α parameter does not.** Measured: after the
mandatory per-cell rescale by `sqrt(diag(Σ))`, every region's induced μ matches `MU_PRIOR` with
Wasserstein-1 ≤ **0.00732** against the pre-declared `SIM02_W1_TOL = 0.10` — SIM-02 holds per region
with two orders of magnitude of headroom. Per-region `ghat` atom mass is **6.79 %**, essentially
identical to the current global 6.78 %, so D-05 multiplies the *number* of atoms by 64 without
inflating their *rate*. But CAR's α is a terrible correlation-length knob at G=8: lag-1 correlation
is 0.136 at α=0.5, 0.438 at α=0.95, 0.692 at α=0.99, 0.947 at α=0.999. Ninety percent of α's range
produces a nearly independent field. A uniform prior on α would put almost all mass on "no spatial
structure", making D-08's inferred correlation length unidentifiable **by parametrization, not by
physics** — the exact prior-echo trap Phase 11's S-1 warns about. The GP lengthscale ℓ behaves far
better (lag-1 correlation 0.135 → 0.992 over ℓ ∈ [0.5, 8]).

**Third, and most consequentially, there is a large pre-existing radial confound and a modest
per-region signal.** Measured on the real image size (1376×1028): the per-patch correlation carries
a per-draw noise sd of **0.0662** in correlation units, while a *constant-ρ* image already shows an
across-region spread of **0.0847** — a ratio of only **1.28**. And holding everything fixed except
`chromatic_eps` 0 → 0.02 (the prior edge) moves the summary by **‖Δ‖₂ = 0.72–2.59**, against an
independent-draw noise floor of 0.9765, with **cor(|Δ|, radius) = +0.48 to +0.83** and a
linear-in-radius model explaining **18–69 %** of the perturbation. Chromatic aberration therefore
injects an above-noise, predominantly *radial*, per-region correlation gradient into **every training
draw** — larger than the maximal 3 px registration shift by a factor of ~3. A spatial map can score
well on D-09's leave-region-out coverage by learning that radial signature and never touching biology.

**Primary recommendation:** Build the lattice prior with stdlib dense linear algebra (no new deps),
parametrize the prior by its *induced lag-1 correlation* rather than by α, express the deviation field
in the 2-D **DCT-II** basis (which is the lattice Laplacian's eigenbasis, is an exact orthonormal
bijection so D-07 is satisfied, and is smoothness-ordered so D-04's "high rank" means simply K = 63),
score SBC in *Gaussian* field space (atom-free) and read ρ through the frozen monotone map at read
time, and gate Stage 1 (D-12) on a **leave-region-out ridge** that needs no trained network at all —
predict region r's ρ from the summary **with row r excluded**, against a predict-the-per-region-prior
baseline, with the row-r-included version as the positive control. That single cheap experiment is the
linear analogue of the whole phase and can kill it for ~20 minutes of CPU.

---

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Lattice prior draw (CAR/GP) → correlated Gaussian field | Simulator prior (`spike/simulator/`) | — | It is π(θ); it belongs beside `sample_prior`, not in the net |
| Copula transform Φ → `MU_PRIOR` quantile → `ghat` elementwise | Simulator prior | — | D-05 is a prior construction; `ghat` is frozen and lives here |
| Field → pixel-resolution ρ map (bilinear upsample) | Forward simulator stage 1 | — | It changes the *generative physics*, not the parameter |
| Composed affine warp (registration ∘ chromatic) | Forward simulator stage 6 | — | Already on HEAD (`ca02b0e8a696269a2b8710b11ca3073edc5ae6c7`); D-06 stacks on it |
| G×G×2 reshape of the encoded summary | Spike-local adapter beside `_summary_row_partition` semantics | `src/amortized/summary.jl` (READ-ONLY) | D-02: summary *content* unchanged; only the net's view changes. Grid coupling stays stated once |
| CNN summary net | Spike net (`spike/npe/p12_architecture.jl`) | — | Deviation from the shipped MLP; must not propagate to `src/` |
| Hierarchical θ (global + deviation coefficients + nuisances + ℓ) | `NormalisingFlow` d-vector | — | Verified: no custom `ApproximateDistribution` needed |
| Per-region ρ / Δρ read-out | Read surface (spike `infer`) | — | Deterministic monotone post-processing of flow draws |
| Leave-region-out predictive coverage | Validation harness (`spike/validation/`) | — | Extends `harness.jl`, per 12-CONTEXT `<code_context>` |
| Result container (`SpatialColocResult`-shaped) | Spike-local struct | `src/results.jl` sketch (READ-ONLY) | Mirrors Phase 13's 13-08 pattern; `src/` stays untouched |

---

## Research-lane posture: the evidence

12-CONTEXT flags research-lane as an **assumption, not a locked answer**, and asks for a verdict.
**The evidence supports it, strongly.** Four independent reasons:

1. **The Phase-7 GO explicitly locks it.** "No further training/gate iteration; Option B … deferred"
   (`.planning/STATE.md`, Current focus 2026-07-24). [VERIFIED: STATE.md]
2. **A Phase-12 net is not drop-in comparable to the shipped read surface.** The shipped bundle reads
   a 128-**row vector** through an MLP and a 7-marginal flow (`NPE_D = 7`,
   `src/amortized/architecture.jl:169`). D-02's CNN reads a **4-D `(G,G,2,K)` array** and D-01's head
   is a 24–72-marginal flow. This is the identical situation Phase 11 recorded as its honesty item #4
   ("Scope every number to the research net"). [VERIFIED: source read]
3. **The shipped artifact is content-hashed and released.** `Artifacts.toml [grid_8]` carries
   `git-tree-sha1 = 90e6b63a…` verified before load (T-7-01). Reshipping means re-cutting a tarball
   and a GitHub Release — a release action, not a phase action. [VERIFIED: STATE.md 07-10 entry]
4. **12-CONTEXT itself defers productionization** ("Productionizing the spatial map — explicitly a
   future phase").

**One caveat the planner must resolve.** ROADMAP SC2 says "`coloc_map(...)` returns a Δρ map +
uncertainty map", which reads as a public-API deliverable. In the research lane that becomes a
**spike-local** `coloc_map` returning a spike-local `SpatialColocResult`-shaped struct, with
`src/results.jl:180-192` left as the untouched sketch. **This is exactly the pattern Phase 13 already
uses** — plan `13-08-PLAN.md` is titled "ThreeHypothesisColocResult **in `spike/`**". Recommend
following it and stating the SC2 scoping in the phase report as a named limit.

---

## Standard Stack

### Core — already present, zero new dependencies

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `LinearAlgebra` (stdlib) | Julia 1.12.6 | Dense `cholesky`/`inv`/`eigen` on the 64×64 lattice precision and covariance | At n = 64 dense is exact, instant and dependency-free. `SparseArrays` buys nothing at this size. [VERIFIED: measured, 40.4 µs/draw] |
| `SparseArrays` (stdlib) | Julia 1.12.6 | Convenient construction of the 4-neighbour adjacency `W` (`spzeros` + fill) | Stdlib; loadable via `@stdlib` without editing `spike/Project.toml`. [VERIFIED: ran under `--project=spike`] |
| `Distributions` | in `spike/Project.toml` | `Truncated(Cauchy(0,0.3),-1,1)` = `MU_PRIOR`; `Normal()` CDF for the copula; `quantile` for the inverse-CDF step | Already the SIM-02 prior surface. [VERIFIED: `spike/simulator/prior.jl:40`] |
| `Interpolations` | in `spike/Project.toml` | Reference bilinear/bicubic upsample (D-06) — but see the matmul note below | Already a dep (used by `forward.jl` stage 6). [VERIFIED: Project.toml] |
| `Flux` | 0.16.10 (pinned) | `Conv`, `SamePad`, `flatten`, `Dense` — the D-02 CNN summary net | Verified to work as a NeuralEstimators summary network on a 4-D array. [VERIFIED: live smoke this session] |
| `NeuralEstimators` | 0.2.1 (pinned) | `PosteriorEstimator`, `NormalisingFlow`, `train`, `sampleposterior`, `assess` | Verified against the installed source at `~/.julia/packages/NeuralEstimators/gFxuZ`. [VERIFIED: source read + live run] |
| `HypothesisTests` | in `spike/Project.toml` | KS / χ² on SBC ranks | Existing SBC apparatus. |
| `Random123` | in `spike/Project.toml` | Counter-based Philox DEV stream | Existing seed discipline. |
| `JLD2`, `StatsBase`, `Statistics`, `Printf`, `Dates` | present / stdlib | Artifacts, z-scoring, reporting | Existing. |

**Installation:**
```bash
# NONE. spike/Project.toml and spike/Manifest.toml stay byte-unchanged.
```

### Alternatives Considered — and rejected

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Hand-written 64×64 dense GP kernel | `KernelFunctions.jl` / `AbstractGPs.jl` | Would give `Matern32Kernel`/`SEKernel` and `kernelmatrix` for free. But it requires adding a dep and running `Pkg.resolve()` on the frozen `spike/Manifest.toml`. The Phase-7 Wave-0 co-resolution gate already failed **once** on a transitive Makie conflict that forced NeuralEstimators toward 0.1.4 (`.planning/STATE.md`, "Resolved 2026-07-03"). The kernel is 1 line; the resolve risk is a milestone-level hazard. **Reject.** |
| Hand-written CAR precision | `GaussianRandomFields.jl` | Generates GRFs on grids with several covariance families. Same resolve risk; also oriented at *large* grids where circulant embedding pays off. At n = 64 it is pure overhead. **Reject.** [ASSUMED: package capabilities from training data, not verified this session] |
| Dense `cholesky(Σ)` | Sparse `cholesky(Q)` + `L' \ z` sampling | The textbook GMRF route and correct for large lattices. At n = 64 the dense route is simpler and equally fast, and — decisively — the **per-cell marginal rescale D-05 requires needs `diag(Q⁻¹)` anyway**, which means forming the dense inverse regardless. **Reject sparse.** |
| DCT-II deviation basis | Eigenvectors of the realized `Q` at a reference α | Exact for CAR but α-dependent, so the head's basis would change with the inferred correlation length. DCT-II is fixed and α-free. **Prefer DCT.** |
| DCT-II deviation basis | PCA of simulated fields | Data-dependent; introduces a fitting step and a leakage surface the Phase-3 loader discipline exists to prevent. **Reject.** |

**Version verification.** `NeuralEstimators` v0.2.1 was verified by reading
`~/.julia/packages/NeuralEstimators/gFxuZ/Project.toml` (`version = "0.2.1"`) and by cross-checking
`spike/Manifest.toml` (`version = "0.2.1"`, `git-tree-sha1 = c780174a2d28714231c0976d5ca4f7c17ac69aa8`).
No registry query was needed and none was made. [VERIFIED: local depot + manifest]

---

## Package Legitimacy Audit

**No external packages are installed by this phase.** The recommendation is stdlib-only
(`LinearAlgebra`, `SparseArrays`), plus packages already pinned in the frozen `spike/Manifest.toml`.

| Package | Registry | Age | Downloads | Source Repo | slopcheck | Disposition |
|---------|----------|-----|-----------|-------------|-----------|-------------|
| `LinearAlgebra` | Julia stdlib | ships with Julia 1.12.6 | n/a | JuliaLang/julia | n/a | Approved — stdlib, no install |
| `SparseArrays` | Julia stdlib | ships with Julia 1.12.6 | n/a | JuliaLang/julia | n/a | Approved — stdlib, no install |
| `NeuralEstimators` | Julia General (already pinned) | n/a | n/a | msainsburydale/NeuralEstimators.jl | n/a | Already installed; tree-sha1 pinned in Manifest |
| `Flux`, `Distributions`, `Interpolations`, `JLD2`, `HypothesisTests`, `Random123`, `StatsBase` | Julia General (already pinned) | n/a | n/a | various JuliaLang/FluxML orgs | n/a | Already installed; no resolve |

**Packages removed due to slopcheck [SLOP] verdict:** none — no package installation is recommended.
**Packages flagged as suspicious [SUS]:** none.

**Honest caveat on tooling.** `slopcheck` was **not run**. It targets npm/PyPI namespaces; the Julia
General registry is outside its coverage, and this phase installs nothing from any registry. The
*rejected* alternatives (`AbstractGPs.jl`, `KernelFunctions.jl`, `GaussianRandomFields.jl`) are named
from training knowledge and were **not** verified against the Julia General registry this session —
they are tagged `[ASSUMED]` above and, since they are rejected, no verification is owed. **If the
planner overrules the zero-dependency recommendation, each proposed package must be gated behind a
`checkpoint:human-verify` task AND behind a `Pkg.resolve()` dry-run that proves NeuralEstimators
stays at 0.2.1** — the co-resolution gate is a Phase-7 precedent, not a formality.

---

## Architecture Patterns

### System Architecture Diagram

```
                        ┌───────────────────────────────────────────┐
  DEV seed (Philox) ───► │  p12_consts.jl  (Tier-1 pre-registration) │
  counters 1..n          │  seeds, ladders, thresholds, ITER=1       │
                        └───────────────┬───────────────────────────┘
                                        │
      ┌─────────────────────────────────▼──────────────────────────────────┐
      │  PRIOR π(θ)   — spike/simulator/p12_prior.jl                       │
      │                                                                    │
      │  ℓ or α  ──► lattice kernel  ──► Σ (64×64)  ──► chol L             │
      │                                    │                               │
      │  z ~ N(0,I₆₄) ──► L·z ──► ÷ sqrt(diag Σ)  ──► standardized field   │
      │                                    │                               │
      │                            DCT-II rotation (orthonormal, exact)    │
      │                                    │                               │
      │                        c = (c₀ global, c₁..c₆₃ deviations)  ◄── θ rows
      │                                    │                               │
      │             Φ(z)  ──►  quantile(MU_PRIOR, ·)  ──►  μ-field         │
      │                                    │                               │
      │                        ghat(·) elementwise  ──►  ρ-field (G×G)     │
      │                                    │            ▲ per-region ATOMS │
      │  7 nuisances (spillover…chromatic_eps) ─────────┴───────► θ rows   │
      └────────────────────┬───────────────────────────────────────────────┘
                           │  ρ-field (G×G)   +  θ nuisances
                           ▼
      ┌────────────────────────────────────────────────────────────────────┐
      │  FORWARD MODEL — spike/simulator/forward.jl (7 stages)             │
      │  stage 1: bilinear upsample ρ-field → pixel res (D-06)             │
      │           a(x,y)=√|ρ(x,y)|, b=√(1−|ρ|), sign(ρ(x,y)) per pixel     │
      │  stages 2..5 unchanged                                             │
      │  stage 6: composed affine (registration ∘ chromatic)  ◄── UNCHANGED│
      │  stage 7 unchanged                                                 │
      └────────────────────┬───────────────────────────────────────────────┘
                           │  [ch1, ch2]
                           ▼
      ┌────────────────────────────────────────────────────────────────────┐
      │  SUMMARY — UNCHANGED src/ path: patch() → correlation() →          │
      │  patch_summary → encode_d01 → 128 rows (64 cont + 64 mask)         │
      └────────────────────┬───────────────────────────────────────────────┘
                           │
        ┌──────────────────┴───────────────────┐
        │  RANDOM REGION MASKING (train-time)  │  ← makes "region unobserved"
        │  zero value + zero mask for k cells  │    in-distribution (see Pitfall 4)
        └──────────────────┬───────────────────┘
                           │
                  reshape to (G, G, 2, K)   ← D-02, column-major, exact
                           │
        ┌──────────────────▼──────────────────┐        ┌────────────────────┐
        │ CNN summary net (Flux Chain)         │──tz──► │ NormalisingFlow(d) │
        │ Conv 3×3 pad=Same ×2 → Conv 1×1 → …  │ dstar  │ d = 1+K+7+1        │
        └─────────────────────────────────────┘        └─────────┬──────────┘
                                                                  │ draws
        ┌─────────────────────────────────────────────────────────▼────────┐
        │ READ SURFACE: inverse-DCT → z-field → Φ → MU_PRIOR quantile →    │
        │ ghat → per-region ρ; paired sample/control pass → per-region Δρ  │
        │ + per-region sd  →  SpatialColocResult-shaped struct (spike)     │
        └─────────────────────────────────────┬────────────────────────────┘
                                              │
   ┌──────────────────────────────────────────┼──────────────────────────────┐
   │                                          │                              │
   ▼                                          ▼                              ▼
 SBC (global c₀ + deviation-energy      LEAVE-REGION-OUT PREDICTIVE     RADIAL-ENERGY
 summary; randomized ranks on any        COVERAGE on 6 real TIFFs        GUARD (S-4)
 ρ-space quantity)                       vs MATCHED ABLATION (ℓ→0)       + offset-grid
                                                                          guard (D-06)
```

### Recommended file structure (spike-local, mirrors the Phase-11/13 layout)

```
spike/
├── validation/
│   ├── p12_consts.jl              # Tier-1 pre-registration: DEV seed, counters, ladders,
│   │                              #   D-12 Stage-1/Stage-2 thresholds, ITERATION_ALLOWANCE=1
│   ├── p12_stage1_ridge.jl        # the cheap leave-region-out ridge gate (NO trained net)
│   ├── p12_sbc.jl                 # global + deviation-summary SBC, randomized ranks
│   ├── p12_coverage.jl            # leave-region-out predictive coverage (sim + real)
│   ├── p12_guards.jl              # radial-energy + offset-grid + eps=0 ablation checks
│   └── run_p12_*.jl               # one reported runner per artifact
├── simulator/
│   ├── p12_lattice.jl             # CAR/GP kernels, dense chol, DCT-II basis, marginal rescale
│   └── p12_prior.jl               # copula: field → Φ → MU_PRIOR quantile → ghat elementwise
├── npe/
│   ├── p12_architecture.jl        # CNN summary net + wide NormalisingFlow + θ transform
│   └── train_p12_npe.jl
├── data/
│   └── p12_generate.jl            # field-aware datagen; stores the drawn G×G lattice (D-07)
└── test/
    ├── test_p12_consts.jl
    ├── test_p12_lattice.jl
    ├── test_p12_prior.jl
    ├── test_p12_architecture.jl
    └── test_p12_coverage.jl
```

### Pattern 1: CAR / GP lattice prior as a *sampling* prior (n = 64, dense)

```julia
# Source: derived from the standard GMRF construction (Rue & Held); VERIFIED by running the
# arithmetic below in the spike environment on 2026-07-27.
using LinearAlgebra, SparseArrays

const G = 8
idx(i, j) = (j - 1) * G + i          # MUST match vec() column-major ordering of the G×G matrix

function lattice_adjacency(G)
    W = spzeros(Float64, G*G, G*G)
    for j in 1:G, i in 1:G, (di, dj) in ((1,0), (-1,0), (0,1), (0,-1))
        i2, j2 = i + di, j + dj
        (1 <= i2 <= G && 1 <= j2 <= G) || continue
        W[idx(i,j), idx(i2,j2)] = 1.0
    end
    W
end

# --- CAR arm: Q = D − αW, α ∈ (0, 1) strictly, SPD guaranteed for α < 1 ------------
function car_sigma(G, α)
    W = lattice_adjacency(G); D = Diagonal(vec(sum(W; dims = 2)))
    Σ = inv(Symmetric(Matrix(D - α * W)))
    return Symmetric(Σ)
end

# --- GP arm: squared-exponential (or Matérn) correlation over lattice coordinates ---
function gp_sigma(G, ℓ; jitter = 1e-8)
    P = [(float(i), float(j)) for j in 1:G for i in 1:G]   # column-major, matches idx()
    K = [exp(-((p[1]-q[1])^2 + (p[2]-q[2])^2) / (2ℓ^2)) for p in P, q in P]
    return Symmetric(K + jitter * I)
end

# --- THE STEP THAT MAKES D-05 CORRECT --------------------------------------------
# Σ does NOT have a unit diagonal (CAR edge cells have smaller marginal variance).
# Without this rescale the copula's per-region marginal is WRONG and SIM-02 fails.
function field_sampler(Σ)
    sd = sqrt.(diag(Σ))
    L  = cholesky(Σ + 1e-12I).L
    return rng -> (L * randn(rng, size(Σ, 1))) ./ sd      # marginally N(0,1) per cell
end
```

**Measured behaviour of this exact code** (G = 8, 4-neighbour lattice, this session):

| CAR α | eigmin(Q) | marginal sd range | r(lag 1) | r(lag 2) | r(lag 4) |
|---|---|---|---|---|---|
| 0.50 | 1.4915 | [0.518, 0.741] | 0.136 | 0.019 | 0.001 |
| 0.90 | 0.3449 | [0.606, 0.909] | 0.351 | 0.143 | 0.039 |
| 0.95 | 0.1738 | [0.658, 0.992] | 0.438 | 0.231 | 0.101 |
| 0.99 | 0.0350 | [0.903, 1.235] | 0.692 | 0.556 | 0.419 |
| 0.999 | 0.0035 | [2.200, 2.370] | 0.947 | 0.922 | 0.883 |

| GP ℓ | posdef | r(lag 1) | r(lag 2) | r(lag 4) | cond(K) |
|---|---|---|---|---|---|
| 0.5 | true | 0.135 | 0.000 | 0.000 | 2.8e0 |
| 1.0 | true | 0.607 | 0.135 | 0.000 | 1.5e3 |
| 2.0 | true | 0.882 | 0.607 | 0.135 | 1.3e9 |
| 4.0 | true | 0.969 | 0.882 | 0.607 | 3.8e9 |
| 8.0 | true | 0.992 | 0.969 | 0.882 | 5.5e9 |

[VERIFIED: measured 2026-07-27, `julia --project=spike`]

**Two consequences the planner must act on:**

1. **Do not put a uniform prior on α.** Ninety percent of α ∈ (0,1) gives lag-1 correlation below
   0.44; the interesting range is α ∈ [0.99, 1). **Parametrize both arms by the same physically
   meaningful quantity — the induced lag-1 correlation r₁ ∈ [0.05, 0.95] — and solve for α or ℓ
   numerically inside the sampler.** This makes the CAR-vs-GP comparison a fair test of *shape*
   (D-03's actual question) instead of a test of two incomparable parametrizations, and it makes
   D-08's inferred correlation length a quantity with a flat, interpretable prior.
2. **The GP arm needs jitter and should probably be Matérn, not SE.** cond(K) reaches 5.5e9 at ℓ = 8
   with an SE kernel. A Matérn ν = 3/2 kernel is far better conditioned at the same effective range.
   Pre-register the jitter value; do not tune it after a failed Cholesky.

**`_smooth_field` is NOT reusable for this.** `spike/simulator/forward.jl:94-99` generates a
*pixel-resolution* Gaussian field by low-pass-filtering white noise at `STRUCT_σ = 6.0` and
standardizing — it is the *texture* generator for stage 1, not a lattice prior, it has no tunable
lattice correlation matrix, and its correlation length is a frozen constant that stage 1's physics
depends on. **Do not repurpose it and do not change `STRUCT_σ`** — `STRUCT_σ ≫ σ_psf` is the D-15
condition the whole ρ→induced-μ monotonicity rests on (`forward.jl:66-70`).

### Pattern 2: The copula (D-05), and how to verify it

```julia
# Source: D-05, implemented and VERIFIED 2026-07-27.
using Distributions
include(joinpath(@__DIR__, "ghat.jl"))                       # ghat, frozen SIM-02
const MU_PRIOR = Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)      # prior.jl:40, verbatim Turing μ-prior

function rho_field(rng, sampler)
    z  = sampler(rng)                       # 64-vector, marginally N(0,1)  (after the rescale!)
    u  = cdf.(Normal(), z)                  # → Uniform(0,1) per region
    μ  = quantile.(Ref(MU_PRIOR), u)        # → MU_PRIOR per region, EXACTLY
    return reshape(ghat.(μ), G, G), reshape(μ, G, G), reshape(z, G, G)
end
```

**Measured verification (M = 20 000 draws, CAR α = 0.95):**

| Check | Measured | Bar |
|---|---|---|
| per-cell empirical sd after rescale | [0.9852, 1.0115] | 1.0 (target) |
| per-region W1(induced μ*, `MU_PRIOR`) — **max over 64 regions** | **0.00732** | `SIM02_W1_TOL = 0.10` |
| per-region W1 — mean over 64 regions | 0.00344 | — |
| per-region `ghat` atom mass (ρ ∈ {−0.99, +0.99}) | **0.0679** | — |
| current GLOBAL 1-D prior atom mass (reference) | 0.0678 | — |

[VERIFIED: measured 2026-07-27]

**Reading:** the copula preserves the marginal essentially exactly, so the SIM-02 claim holds
**per region** with 13× headroom on the pre-declared tolerance. And per-region atom *rate* is
unchanged at 6.79 % — D-05 does not make the atom problem worse per region; it makes it **64× more
frequent per draw**, which is exactly why S-3's randomized-rank budget must be in the plan from the
start rather than discovered.

**The SIM-02 verification pattern to reuse:** `spike/simulator/prior.jl:42-46` pre-declares
`SIM02_W1_TOL = 0.10` *before* the calibration run. Mirror it exactly: declare a per-region W1
tolerance in `p12_consts.jl` **before** the field simulator runs, and assert `maximum(w1_per_region)
≤ tol` — not the mean. (The numbers above say this will pass comfortably; declaring the bar in
advance is nonetheless the discipline, and the measured 0.00732 must **not** be used to set the
tolerance — that would be reverse-engineering a bar from a result.)

### Pattern 3: D-06 field application — use two matmuls, not a scalar interpolant

```julia
# Separable bilinear upsample of a G×G field to (m, n) pixels: F_full = A * F * B'
function bilinear_op(G, m)
    A  = zeros(m, G)
    for (r, x) in enumerate(range(1, G; length = m))
        i = clamp(floor(Int, x), 1, G-1); t = x - i
        A[r, i] = 1 - t; A[r, i+1] = t
    end
    A
end
# stage 1 then uses per-pixel a/b/sign instead of scalars:
#   ρpx = A * ρfield * B'
#   a   = sqrt.(abs.(ρpx));  b = sqrt.(1 .- abs.(ρpx))     # a² + b² = 1 holds PER PIXEL
#   ch1 = _softplus.(a .* L .+ b .* ε1)
#   ch2 = _softplus.(sign.(ρpx) .* a .* L .+ b .* ε2)
```

**Measured cost** (this session):

| Output size | scalar `Interpolations` comprehension | separable matmul |
|---|---|---|
| 512×512 | — | 0.0003 s |
| 1024×1024 | — | 0.0018 s |
| 1376×1028 | **0.3622 s** | **0.0252 s** |
| 2048×2048 | — | 0.0273 s |

[VERIFIED: measured 2026-07-27]

**This is a 14× difference and it is load-bearing.** Naively calling a scalar interpolant per pixel
costs 0.36 s per sample at the real image size; over 50 000 pairs on the F5 mixture that is
**~5 hours added to datagen**, versus ~8 minutes for the matmul form. A plan that says "bilinear
upsample" without specifying the implementation will silently produce the 5-hour version.

**Offset-grid guard (D-06 / `<specifics>`).** The same construction makes the guard trivial: build
`A`/`B` from `range(1 - 0.5, G + 0.5; length = m)` instead of `range(1, G; length = m)` and the field
lattice is offset by half a cell relative to the summary's patch boundaries. Pre-register: the trained
model's per-region RMSE and predictive coverage must not degrade by more than a declared amount on the
offset-grid evaluation arm. If it does, the model is reading grid-aligned block edges.

### Pattern 4: The hierarchical head — verified API, no custom distribution needed

**This was the single highest-uncertainty question and it is now settled.** Verified two ways:
reading the pinned v0.2.1 source, and running it.

```julia
# Source: ~/.julia/packages/NeuralEstimators/gFxuZ/src/Estimators/PosteriorEstimator.jl:1-8
#         and .../ApproximateDistributions/NormalisingFlow.jl:1-32   [VERIFIED: local depot, v0.2.1]
#
#   PosteriorEstimator(summary_network, q::ApproximateDistribution)
#   PosteriorEstimator(summary_network, num_parameters::Integer; num_summaries::Integer, q = nothing, kwargs...)
#   NormalisingFlow(d::Integer, dstar::Integer; num_coupling_layers = 6, use_act_norm = true,
#                   backend = nothing, kwargs...)
#   # ApproximateDistributions.jl:29 also accepts num_summaries as a KEYWORD:
#   #   (::Type{T})(num_parameters::Integer; num_summaries::Integer, kwargs...)
#
# The summary network is called directly:  _summarystatistics(estimator, Z) = estimator.summary_network(Z)
#   (Estimators.jl:76)  ⇒ ANY Flux model accepting Z works, including a CNN on a 4-D array.

using Flux, NeuralEstimators
G, DSTAR = 8, 128
D = 1 + K_dev + 7 + 1        # global c₀ + K deviation coefficients + 7 nuisances + correlation length

cnn = Chain(
    Conv((3,3), 2  => 16, gelu; pad = SamePad()),   # 8×8×2  → 8×8×16
    Conv((3,3), 16 => 32, gelu; pad = SamePad()),   # 8×8×16 → 8×8×32
    Conv((1,1), 32 =>  8, gelu),                    # 8×8×32 → 8×8×8   (channel bottleneck)
    Flux.flatten,                                   # → 512
    Dense(512, DSTAR),
)
q   = NormalisingFlow(D; num_summaries = DSTAR, num_coupling_layers = 10, depth = 2, width = 128)
est = PosteriorEstimator(cnn, q)                    # q passed POSITIONALLY as an INSTANCE

# Training data: Z as a 4-D array, K data sets along the LAST dimension.
# train.jl:668 wraps MLUtils.DataLoader(f32(data); batchsize, shuffle) which batches the last dim.
est = train(est, θ_train, θ_val, Z_train, Z_val;      # Z_* :: Array{Float32,4} == (G,G,2,K)
            epochs = ..., batchsize = 64, use_gpu = false)

draws = sampleposterior(est, Z[:, :, :, 1:1]; N = 2000)   # → Matrix{Float32}, D × N
```

**Live smoke result (this session):**
- `cnn(randn(Float32, 8, 8, 2, 512))` → `(128, 512)`. ✅
- `NormalisingFlow(72; num_summaries = 128, num_coupling_layers = 10, depth = 2, width = 128)` builds;
  `numdistributionalparams = 1440`. ✅
- `train(...)` on `(8,8,2,512)` completes. ✅
- `sampleposterior(est, Z[:,:,:,1:1]; N = 200)` → `Matrix{Float32}`. ✅

**Answer to "does the flow just get a wider θ vector, or is a custom distribution needed?"** — the
flow just gets a wider θ vector. `NormalisingFlow` is dimension-generic (`CouplingLayer` splits `d`
into ⌊d/2⌋ / ⌈d/2⌉ and works for any `d ≥ 2`). **No custom `ApproximateDistribution` is required.**
[VERIFIED: source + live run]

**Note on the docs' spatially-varying example.** `msainsburydale.github.io/…/examples/data_gridded_nonstationary`
does exactly this phase's shape (image-in, field-out) — but with `PointEstimator` and a `Unet`, i.e.
a *point* estimator with no posterior. `Unet` is **not exported by v0.2.1**
(`grep -rn "Unet" src/` → no matches; exports list at `NeuralEstimators.jl:49-64`). So there is no
built-in image-to-image *distributional* head available. **D-01's finite-dimensional low-rank vector
head is not just a calibration convenience — it is the only route the pinned library supports.**
[VERIFIED: source grep + CITED: msainsburydale.github.io/NeuralEstimators.jl/dev/examples/data_gridded_nonstationary]

### Pattern 5: The deviation basis — 2-D DCT-II, and why it resolves D-01 vs D-07

D-01 wants a **low-rank** deviation field; D-07 wants the scored quantity to be **the sampled
parameter, with no deterministic transform in between**. These look contradictory. They are not,
provided the basis change is an **exact orthonormal bijection**.

- The 4-neighbour lattice graph Laplacian on a G×G grid with free (Neumann) boundaries is a Kronecker
  sum of two path-graph Laplacians, whose eigenvectors are exactly the **DCT-II basis**
  `φ_p(i) ∝ cos(π p (i − ½) / G)`. So the 2-D DCT-II basis is the natural smoothness-ordered
  eigenbasis of the lattice. [ASSUMED: standard spectral-graph result from training knowledge; the
  *numerical* claims in this document do not depend on it — the DCT basis is a valid orthonormal
  basis regardless of whether it exactly diagonalizes `Q`]
- **Caveat, stated rather than glossed:** for `Q = D − αW` with α < 1 the degree matrix `D` is *not*
  a multiple of `I` (boundary cells have degree 2 or 3), so DCT-II diagonalizes `Q` **exactly only in
  the α → 1 (intrinsic CAR) limit**. For α < 1 it is an excellent but approximate diagonalization.
  This costs nothing: the basis only has to be orthonormal and smoothness-ordered, which it is
  exactly.
- **Because the DCT is orthonormal and exactly invertible, ranking a DCT coefficient IS ranking a
  sampled parameter.** D-07's prohibition is specifically about *lossy* transforms — the `ghat` clamp
  (which manufactures atoms) and "the cell mean of the interpolated field" (a projection that discards
  information). An orthogonal rotation discards nothing. **Sample the field, express it in the DCT
  basis, and let the DCT coefficients BE the θ rows.** Then D-01 and D-07 are simultaneously satisfied.
- **"High-rank" (D-04) then has a concrete meaning: K = G² − 1 = 63**, i.e. every non-constant DCT
  mode. The head is then *full rank on the deviation space* and cannot under-fit a rough CAR field by
  construction — which is precisely what D-04 asks for. Reconstruction is `F = idct(c)`, exact.
- **The production head (chosen after the mini-spike) is then a truncation to the first K modes** in
  the smoothness ordering, and the truncation error is directly measurable on simulated fields — so
  "what rank is enough?" becomes an arithmetic question, not a guess.

### Anti-Patterns to Avoid

- **Sampling the CAR field without rescaling by `sqrt(diag(Σ))`.** Edge cells have marginal sd 0.658
  vs interior 0.992 at α = 0.95. Skipping the rescale silently breaks D-05's marginal-preservation
  claim, breaks SIM-02 per region, and puts a *position-dependent* prior on ρ — which is exactly the
  radial/edge artifact class the phase is trying to guard against.
- **A uniform prior on the CAR α.** See Pattern 1.
- **A `GlobalMeanPool` or `DeepSet` anywhere in the summary net.** Permutation-invariant pooling over
  patches *is* the exchangeable pooling this phase exists to replace (12-CONTEXT records this as an
  SC1 correction and a deferred item). A flatten-then-Dense head after the convs is the point.
- **Flattening 8·8·64 = 4096 into a wide Dense.** Measured 7× slower with no representational
  argument at G = 8 (see Compute Budget). Use a 1×1 channel bottleneck first.
- **Transposing the G×G reshape.** `encode_d01` builds `vec(coalesce.(M, 0.0))` column-major from the
  G×G `correlation()` matrix (`src/amortized/summary.jl:72-76`), so `reshape(rows[1:G^2], G, G)`
  recovers `M` exactly with **no transpose**. Any radius/adjacency computation must use the same
  column-major indexing (`idx(i,j) = (j-1)*G + i`). An off-by-transpose here would silently rotate the
  lattice relative to the image and would be nearly undetectable in a symmetric test.
- **Repurposing `_smooth_field` or changing `STRUCT_σ`.** See Pattern 1.
- **Appending Phase-12 testsets after `test_p13_correction.jl` in `spike/test/runtests.jl`.** That
  testset contains committed measured misses and **throws**, aborting every subsequent `include`
  (documented in the runtests.jl comment block). Phase-12 includes must go **before** it or they will
  silently never run.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Rank-uniformity testing | A custom KS / χ² | `HypothesisTests.ExactOneSampleKSTest` / `ChisqTest` | Already the project's rule; already wired in `spike/validation/sbc.jl` |
| Atom handling in SBC | Truncating the μ-prior | **Randomized ranks** (recipe in the skill's `references/sbc-calibration.md`) | Spike 013 measured that truncation doubles high-\|ρ\| bias and breaks ADVI comparability; randomized ranks cost nothing |
| Per-region posterior read | A new inference path | Extend `spike/validation/harness.jl` | 12-CONTEXT `<code_context>`: "the per-region coverage runs should extend it rather than build a parallel harness" |
| Degenerate-region handling | Silent zeros | The `LOCAL_MAP_SENTINEL` + forced-OOD-flag discipline (`src/amortized/local_map.jl:42-53`) | T-7-04: an unscorable region must never read as "measured no difference" |
| Pre-registration file | Ad-hoc constants in scripts | The two-tier `p11_consts.jl` pattern (`spike/validation/p11_consts.jl`) | Append-never-edit makes the git history itself the audit trail |
| Seed disjointness | A comment asserting it | `_p11_forbidden()`-style **executable** `@assert`, with derived gate seeds RECOMPUTED | `p11_consts.jl:104-135`: "a comment naming a derived seed is not evidence" |
| Cholesky / eigen on 64×64 | Custom decompositions | `LinearAlgebra.cholesky` / `eigen` | Dense, exact, instant at n = 64 |
| GP kernel matrix | `AbstractGPs.jl` | A 1-line comprehension | The dependency's resolve risk exceeds its value at n = 64 (see Alternatives) |
| Bilinear upsample | A scalar interpolant loop | Two matmuls (`A * F * B'`) | 14× measured, and it makes the offset-grid guard a one-line change |

**Key insight:** in this phase the usual "don't hand-roll" instinct **inverts** for the spatial
prior. The library-shaped solutions (`AbstractGPs`, `GaussianRandomFields`, sparse GMRF machinery) all
exist to make *large* spatial problems tractable. At G = 8 the problem is a 64×64 dense matrix, and the
binding constraint is not implementation effort — it is the frozen `spike/Manifest.toml` and the
NeuralEstimators 0.2.1 pin that the entire v2.0 calibration story is scoped to. Writing 15 lines of
`LinearAlgebra` is strictly cheaper than a resolve.

---

## Runtime State Inventory

Phase 12 is not a rename/refactor phase, but it does write persistent artifacts and consume seed
streams, so the equivalent inventory is worth stating.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | Phase-11 pool at `spike/data/cache/<hash>/` (50 000 pairs) is **reusable read-only** for the ε-identifiability check but **cannot** supply per-region truth (it has a single global ρ per sample). A new Phase-12 pool is required. | New content-hashed cache dir; do not overwrite the P11 pool |
| Live service config | None — no external service. | None |
| OS-registered state | None. | None |
| Secrets / env vars | None. | None |
| Build artifacts | `artifacts/amended_v2/grid_8/` (shipped, content-hashed via `Artifacts.toml`, tree-sha1 `90e6b63a…`) must stay **byte-unchanged**. Phase-12 artifacts go under a new spike-local root. | Assert byte-unchanged at phase close, as Phase 11 did |
| Reserved seed streams | `NPE_MASTER_SEED 0xC0FFEE`, `VAL_MASTER_SEED 0x5BC0FFEE`, `VAL_FIX_SEED 0xF1F7ED`, `DEFAULT_MASTER_SEED 0x1`, `CORPUS_MASTER_SEED 0xc05eed`, `P11_DEV_SEED 0x0B11DE71`, `P11_RESEARCH_BURNED` (4 keys), `F2_DEV_SEEDS` (2), `SPIKE_DEV_SEEDS` (6), derived `PROD_SEED`/`PROD_SEED_V2` (8), plus **Phase 13's** DEV seed and burned keys in `spike/p13/`. | `p12_consts.jl` must recompute and forbid ALL of them, **including Phase 13's**, which `p11_consts.jl` predates and therefore does not list |
| Concurrent writers | **Phase 13 is EXECUTING** (`.planning/STATE.md`: "Phase: 13 … EXECUTING, Plan 1 of 16"). Phase 11's closure records that concurrent agents on one branch race on the **git index**, not only on files. | Use pathspec-scoped commits (`git commit -- <paths>`); do not `git commit` without a pathspec |

**Explicitly checked and empty:** no database, no OS task registration, no secret, no external service
carries a Phase-12-relevant string.

---

## Common Pitfalls

### Pitfall 1: The chromatic radial confound can fake a pass (12-CONTEXT S-4) — MEASURED, and it is large

**What goes wrong:** `CHROMATIC_PRIOR = Uniform(-0.02, 0.02)` (`spike/simulator/prior.jl:57`) applies
a **radial** magnification difference to channel 2 in stage 6. The effective per-patch correlation
therefore varies **radially in every training draw** with no biological field present. A spatial map
learns the radial signature, predicts held-out regions from neighbours partly *because the whole image
has a deterministic radial gradient*, and passes D-09.

**Measured magnitude.** Two independent measurements agree that this is a first-order effect.

*From Phase 11's own pre-registered probe artifact* (`spike/validation/p11_probe_report.jld2`, F5
mixture, paired design, metric `paired_l2_rows_1_64`) — compared against the independent-key noise
floor **0.9765** that `11-DIAGNOSIS.md §1` measured on the same metric:

| chromatic ε | paired ‖Δs‖₂ | ÷ noise floor | Δρ-equivalent |
|---|---|---|---|
| 0.0025 | 0.0940 | 0.096 | 0.0117 |
| 0.005 | 0.2172 | 0.222 | 0.0308 |
| 0.010 | 0.5513 | 0.565 | 0.0826 |
| **0.020 (prior edge)** | **1.2417** | **1.272** | **0.1895** |
| 0.030 (beyond prior) | 1.7440 | 1.786 | 0.2672 |

Compare the **entire** registration ladder: 3.0 px (`LAMBDA_MAX`) gives ‖Δs‖₂ = 0.3966, ratio 0.406.
**Chromatic ε at its prior edge moves the summary ~3× more than the maximal registration shift, and
it is the only nuisance in the model that clears the per-draw noise floor.**

*Direct measurement of the spatial structure* (this session, 6 prior draws, 1376×1028, paired same-key
ε = 0 vs ε = 0.02):

| base | ρ_true | ‖Δ‖₂ (64 rows) | cor(\|Δ\|, radius) | radial-linear R² |
|---|---|---|---|---|
| 1 | −0.471 | 1.4867 | +0.718 | 0.520 |
| 2 | −0.021 | 0.7166 | +0.480 | 0.182 |
| 3 | +0.664 | 2.5920 | +0.818 | 0.674 |
| 4 | +0.381 | 1.7895 | +0.792 | 0.630 |
| 5 | +0.161 | 1.2366 | +0.644 | 0.416 |
| 6 | −0.820 | 2.0295 | +0.832 | 0.693 |

[VERIFIED: measured 2026-07-27]

**Reading:** the perturbation is unambiguously radial (mean cor +0.71; a single linear-in-radius term
explains a mean 52 % of the 64-dimensional perturbation), and it grows with |ρ| exactly as physics
predicts (magnification mismatch degrades an existing correlation; at ρ ≈ 0 there is nothing to
degrade — base 2, ρ = −0.021, is the smallest effect).

**How to avoid — three pre-registered guards, all cheap:**

1. **Radial-energy report (mandatory).** Regress the recovered per-region deviation field on a radial
   basis (constant + radius, or the first two radial DCT-like modes) and report the fraction of field
   energy that is radial, on both simulated and real evaluations. Pre-register a ceiling. A map whose
   deviation field is >X% radial is reporting chromatic aberration.
2. **ε = 0 evaluation ablation.** Re-run the D-09 leave-region-out coverage on an evaluation pool with
   `chromatic_eps` held at exactly 0. If the coverage/score advantage over the D-10 matched ablation
   *vanishes*, the spatial win was the chromatic signature.
3. **Radially-orthogonalized field arm.** Draw the CAR/GP field, project out its radial component, and
   check the map still recovers it. This is the positive control that separates "reads radial nuisance
   structure" from "reads the biological field".

Additionally, **report `chromatic_eps`'s own shrinkage/vacuity** through the existing F3 diagnostic
(`shrinkage = post_sd / prior_sd` + boolean `vacuous`). Because ε clears the noise floor, it is
plausibly the *first identified nuisance* in this project's history — and if it is identified, the net
can explain the radial gradient away and the confound is mitigated. If it is vacuous, it cannot, and
the confound is live. Either answer is publishable; not measuring it is not.

### Pitfall 2: The per-region signal is much weaker than the global one

**What goes wrong:** Phase 11's ridge recovered *global* ρ from the 128-row summary at ratio 0.157
(84 % error reduction), and it is tempting to read that as licence for per-region recovery.
12-CONTEXT S-2 correctly says it "does **not** prove per-region recovery". Measured, the per-region
picture is much tighter.

**Measured** (6 draws, 1376×1028, independent-key paired):

| Quantity | Measured |
|---|---|
| per-patch correlation, per-draw noise sd (raw correlation units) | **0.0662** |
| within-image across-region SD of the 64 patch correlations, under a **constant-ρ** simulator | **0.0847** |
| ratio (apparent across-region spread ÷ pure noise) | **1.28** |

[VERIFIED: measured 2026-07-27]

**Reading:** under the *current* (spatially constant ρ) simulator, most of the apparent
region-to-region variation in an image is estimator noise. Per-region recovery is therefore possible
only when the injected field's per-region contrast is comparable to or larger than ~0.066 in
correlation units. Under the D-05 copula every region's marginal is `MU_PRIOR` regardless of the
correlation length, so the **within-image contrast is controlled entirely by the correlation length**:
short length → large contrast, easily recoverable; long length → the field is nearly constant and
there is nothing per-region to recover (and, correctly, nothing to lose).

**Consequence for the plan:** (a) the D-12 Stage-1 gate must be run *across a correlation-length
ladder*, not at one value, and the report must state the length at which per-region recovery
crosses the noise floor — that is a genuine scientific deliverable; (b) the noise sd is
image-size dependent (patches at 512² carry ~1/7 the pixels of 1376×1028), so the ladder must be run
on the F5 mixture the net trains on, not on a single convenient size.

### Pitfall 3: "Alive and used" ≠ "learned" for the correlation length (12-CONTEXT S-1)

**What goes wrong:** the exact Phase-11 failure, repeated. Phase 11 read λ-marginal widening of
4.84–9.05× as the conditioning being "alive and used" when **12.0 = `LAMBDA_MAX/LAMBDA_MIN` is the
prior-echo ceiling** reachable by a network that ignores the image entirely.

**Why it happens here:** if D-08's correlation length ℓ is inferred and its posterior is wide,
that is trivially consistent with the flow reproducing its prior. Reporting "ℓ's posterior responds"
proves nothing.

**How to avoid — pre-register the Phase-11 ridge, per S-1:** ridge-regress ℓ on the 64 continuous
summary rows (there is **no** conditioning row, since D-08 infers rather than conditions), against the
baseline of predicting the prior mean of ℓ, with a `ρ_true` positive control on the same predictors
and the same split. **Add one engineered feature that makes the test fair:** the empirical **lag-1
spatial autocorrelation (Moran's I) of the 64 continuous rows**, which is the natural near-sufficient
statistic for ℓ. A ridge on raw rows alone is a weak probe for a *second-order* quantity; a ridge on
raw rows + Moran's I is the honest one. Report both.

### Pitfall 4: The masked-region encoding may be out of distribution

**What goes wrong:** D-09 holds out lattice regions. The natural encoding is exactly what
`encode_d01` already produces for an unusable patch — value 0 in the continuous row, 0 in the mask row
(`src/amortized/summary.jl:60-76`). So the mask row **does** give a principled "this region
unobserved" encoding, and `local_map.jl`'s sentinel discipline is the precedent for treating it as
"unscorable, not zero".

**But:** in training, that encoding only arises when a patch falls below the ≥15-surviving-pixel floor
of `correlation`/`_exclude_zero`. On the F5 image-size mixture with the softplus background floor
(`BG_FLOOR = 0.02`, which exists precisely so `_exclude_zero` does *not* thin patches), masked patches
are **rare**. A net that has essentially never seen a masked region will behave arbitrarily when D-09
masks one, and the coverage number will measure out-of-distribution behaviour rather than borrowing
strength.

**How to avoid:** add **random region masking (MCAR) as a training-time augmentation** — for each
training sample, mask k ~ some small distribution of randomly chosen regions (set value 0 and mask 0)
before the reshape. Pre-register the masking distribution. This is the single change that makes D-09
measure what it claims to measure, and it is cheap. Report the realized masked-region rate in
training so a reader can check that D-09's held-out configuration is in-distribution.

**Warning sign:** if leave-region-out predictive coverage is wildly off nominal for *both* the spatial
model and the D-10 ablation, suspect this before suspecting the prior.

### Pitfall 5: The D-09 gate may not discriminate

**What goes wrong:** the predictive interval for a held-out region's *observed* correlation entry is
the convolution of (posterior uncertainty about the region's ρ) with (per-region observation noise,
measured at **0.0662**). If observation noise dominates, both the spatial model and the ℓ→0 ablation
will show approximately nominal coverage, and "beats independent pooling in coverage" will be a
coin flip.

**How to avoid:** gate on **two** things, as most predictive-checking practice does —
*calibration* (coverage ≈ nominal, required of **both** arms; a miscalibrated arm is disqualified, not
"better") and *sharpness given calibration* (a **proper scoring rule** — mean log predictive density
or CRPS of the held-out observation). The spatial prior's claim is that it makes held-out predictions
*sharper at the same calibration*, which is exactly what a proper scoring rule measures and what raw
coverage cannot. Pre-register both, and pre-register that a coverage-only win with no scoring-rule
improvement is **not** a pass.

### Pitfall 6: D-10's baseline is not a "prior-only" baseline

**What goes wrong:** D-10 says "independent pooling has no mechanism to predict a held-out region
beyond the prior". Taken literally this is **not true**, and stating it in the report would be a
falsifiable overclaim.

**Why:** with the correlation length driven to zero the *regions* become conditionally independent —
but the other 63 observed regions still inform (a) the shared **nuisances** (noise level, label
efficiency, spillover, and especially `chromatic_eps`, whose radial gradient is global), and (b) any
global term in the head. So the ablation's predictive for the held-out region is the prior
**conditional on the shared nuisances and the global term**, which is strictly narrower than the
marginal prior.

**How to avoid:** this is actually the *right* baseline — it isolates the spatial-borrowing mechanism
exactly, which is what D-10 wants. Just describe it correctly: "no spatial borrowing, full nuisance
and global borrowing". Do not write "beyond the prior". A third, genuinely prior-only arm (predict the
held-out region from its marginal prior alone) is worth reporting as a floor, so the reader can see how
much of the ablation's performance is nuisance/global borrowing versus spatial borrowing.

### Pitfall 7: Per-region atoms, unrandomized, cost ~15 coverage points

**What goes wrong:** applying `ghat` elementwise (D-05) puts a prior atom at ρ = ±0.99 on **every
region**, at a measured rate of **6.79 %** per region. SBC assumes a continuous prior; θ* on an atom
is forced to a rank extreme.

**Measured cost of ignoring it:** Phase 11's coverage run (no atom randomization) reported empirical
coverage **0.723–0.767 against nominal 0.90**, z-sd ≈ 1.4 — intervals ~40 % too narrow — and
`11-DIAGNOSIS` attributes this to the known prior-atom artifact rather than a calibration defect.

**How to avoid:** budget randomized ranks from the start, using the recipe in the skill's
`references/sbc-calibration.md`:
```julia
u = (R[:, p] .+ 0.5) ./ (L + 1)
x = theta[:, p]; lo, hi = minimum(x), maximum(x)
atom = (x .== lo) .| (x .== hi)
ur = copy(u); ur[atom] = rand(count(atom))
pvalue(ExactOneSampleKSTest(ur, Uniform(0, 1)))
```
**Better still: avoid the atoms entirely for the deviation field.** If the θ rows are the *Gaussian*
field's DCT coefficients rather than the ρ-field, there are **no atoms at all** in those rows — the
clamp lives only in the read-time `ghat` map, which is a monotone elementwise transform and therefore
preserves credible intervals exactly. Randomized ranks are then needed only for whatever ρ-space
quantity is additionally reported. This is a strictly better parameterization and costs nothing.

### Pitfall 8: Over-powered SBC on non-identified rows (spike 011)

**What goes wrong:** M = 2000 point-null SBC resolves a 0.059-SD marginal drift, and a neural flow's
irreducible drift on a **non-identified** parameter is ~0.04–0.07 SD — no capacity lever removes it
(spike 012 PARTIAL). With D-01's head the parameter count *rises* (to 24–72), so the number of
non-identified rows rises with it, and a naive per-row point-null gate will fail loudly and
uninformatively.

**How to avoid:** apply the project's settled policy from `references/sbc-calibration.md` —
strict point-null SBC **only** for the targets (the global term and Δρ), and a **nuisance-appropriate
equivalence test** (TOST, per the `07-NUISANCE-SBC-SPEC-DRAFT.md` and Phase-11 D-08 pattern) for the
nuisances and for the high-index deviation coefficients, which are non-identified almost by
construction.

---

## Code Examples

### Reshape the encoded summary to (G, G, 2, K) — exact, no transpose

```julia
# Source: derived from src/amortized/summary.jl:60-76 (encode_d01) and :89-98
# (_summary_row_partition). VERIFIED by construction: encode_d01 uses column-major vec().
"""
    reshape_summary(Zraw::AbstractMatrix, G::Integer) -> Array{Float32,4}

Zraw is (2G²) × K, exactly what encode_d01 produces (rows 1:G² imputed correlations,
rows G²+1:2G² the binary present-mask). Returns (G, G, 2, K):
channel 1 = the correlation lattice, channel 2 = the present-mask lattice.
`reshape(view, G, G)` recovers the original G×G `patch_summary` matrix EXACTLY — no transpose.
"""
function reshape_summary(Zraw::AbstractMatrix, G::Integer)
    n = G^2
    size(Zraw, 1) == 2n || throw(DimensionMismatch("expected $(2n) rows, got $(size(Zraw,1))"))
    K = size(Zraw, 2)
    out = Array{Float32,4}(undef, G, G, 2, K)
    @views out[:, :, 1, :] .= reshape(Zraw[1:n,       :], G, G, K)
    @views out[:, :, 2, :] .= reshape(Zraw[n+1:2n,    :], G, G, K)
    return out
end
```

**Test this against a round-trip, not by eye:** for a random `M::Matrix{Union{Float64,Missing}}`,
assert `reshape_summary(encode_d01(M), G)[:, :, 1, 1] == Float32.(coalesce.(M, 0.0))`. A transposed
implementation passes visual inspection and fails this.

**Standardization note (open question, see OQ-3):** the frozen `zt` z-scores each of the 64 continuous
rows *independently* — a per-pixel affine map, which fights a convolution's weight sharing. A single
shared scalar mean/sd across the 64 continuous rows would preserve translation equivariance. Both are
defensible; **whichever is used must be recorded as a declared deviation**, and the mask rows must
never be z-scored (`_summary_row_partition`'s reason: "z-scoring a 0/1 mask would re-couple folds
through the mask mean").

### The D-12 Stage-1 gate: leave-region-out ridge, no trained network

```julia
# Source: the design of spike/validation/run_p11_recovery.jl, adapted to the per-region question.
# The Phase-11 template's three load-bearing features are preserved EXACTLY:
#   (1) the predictor set EXCLUDES the row that would answer the question trivially,
#   (2) the baseline is "predict the conditional prior mean", computed EMPIRICALLY,
#   (3) a positive control proves the harness is live.
#
# Predictors : the 64 continuous rows (+64 mask rows), with ROW r ZEROED AND MASKED.
# Target     : the drawn lattice value at region r (D-07 ground truth).
# Baseline   : predict the per-region prior mean, from the realized draws.
# +control   : the SAME ridge WITH row r included (must be far below 1.0, else the harness is dead).
#
# Run across a correlation-length ladder. The decisive number is the ratio at each rung:
#   ratio(r-excluded) < 1  ⇒ neighbours carry information ⇒ spatial borrowing exists
#   ratio(r-excluded) ≈ 1  ⇒ NO borrowing available ⇒ DESCOPE before the training spend (D-12 Stage 1)
```

This costs one modest field-simulator pool (a few thousand samples at 512²) plus closed-form ridge
solves, i.e. tens of minutes. It is the **linear analogue of the entire phase** and it is falsifiable
before any network exists — which is precisely what D-12 Stage 1 asks for and what Phase 11 proved is
worth building first.

### The leave-region-out predictive construction (D-09)

```
For each real image (6 committed TIFFs) and each region r in 1:64:
  1.  Z_full = encode_d01(patch_summary(mci, 8))          # the observed 128 rows
  2.  Z_r    = mask_region(Z_full, r)                     # value 0, mask 0 at r
  3.  draws  = sampleposterior(est, reshape_summary(Z_r, 8); N)   # amortized, one pass
  4.  ρ_r^(s) = read_region(draws, r)                     # posterior for the HELD-OUT region
  5.  predictive for the OBSERVED entry:
         y_r^(s) = ρ_r^(s) + e^(s),   e^(s) ~ observation-noise model
  6.  score: (a) is the observed Z_full[r] inside the nominal predictive interval?  [calibration]
             (b) log predictive density / CRPS of Z_full[r] under {y_r^(s)}          [sharpness]
  7.  repeat for the D-10 matched ablation (ℓ→0) and (optionally) a marginal-prior floor arm
```

**The observation-noise model is a declared modelling assumption and must be pre-registered and
stated as a limit.** Two defensible forms:

- *Empirical:* calibrate the per-region noise sd from simulation as a function of image size and
  |ρ| — this session measured 0.0662 pooled at 1376×1028, with a clear |ρ| dependence
  (0.0799 at ρ ≈ 0.16, 0.0393 at ρ ≈ −0.82), consistent with the classical `(1 − r²)/√n_eff` form.
- *Fisher-z:* transform to `z = atanh(r)`, where `Var(z) ≈ 1/(n_eff − 3)` with a single `n_eff`
  calibrated once from simulation per image size and frozen in `p12_consts.jl`.

Prefer Fisher-z: it needs one calibrated constant instead of a surface, it is the standard treatment,
and it makes the |ρ| dependence automatic. **State plainly in the report that D-09 validates the
*joint* of (posterior + this noise model), not the posterior alone** — that is the honest scope, and
it is the same class of caveat Phase 11 attached to its hybrid quadrature column.

---

## State of the Art

| Old approach (in this repo) | Current approach (Phase 12) | Why it changed |
|---|---|---|
| Global scalar ρ, exchangeable patch pooling | Lattice prior + per-region deviation field | ROADMAP Phase 12 goal; `local_map.jl:27-32` states a finer grid cannot get there |
| `LocalColocMap`: per-tile Δρ point estimates, no uncertainty | `SpatialColocResult`-shaped: per-region Δρ **+ per-region sd** | `local_map.jl:55-79` docstring: "no posterior draws, no Bayes factor and no per-region uncertainty" |
| MLP over a 128-row vector | CNN over a (G, G, 2) lattice | D-02: same information, added spatial inductive bias |
| θ = 7 rows (shipped) / 8 rows (spike lane) | θ = 1 + K + 7 + 1 rows | D-01 + D-08 |
| SC3 "beats independent pooling in coverage on ≥1 real image" | SC3 amended: leave-region-out **predictive** coverage vs a **matched ablation** | D-09/D-10: no ground-truth ρ field exists on a real image, and the old baseline had no uncertainty |
| SC1 "CNN/DeepSet summary" | CNN only | 12-CONTEXT: a DeepSet is permutation-invariant over patches — the very pooling being replaced |

**Deprecated / not applicable here:**
- **The finer-grid family** is ruled out with evidence, not opinion (`11-CLOSURE.md` #6): 8×8→32×32
  multiplies signal by 8.9× and noise floor by 10.3×. Do not propose G > 8 as a route to per-region
  resolution.
- **`Unet` / image-to-image estimators** are in the NeuralEstimators *docs* but not in the pinned
  v0.2.1 exports, and the documented example uses `PointEstimator` (no posterior). Not a route.

---

## Assumptions Log

| # | Claim | Section | Risk if wrong |
|---|---|---|---|
| A1 | The 2-D DCT-II basis is the eigenbasis of the free-boundary lattice graph Laplacian (and hence the intrinsic-CAR precision) | Pattern 5 | Low. The basis only needs to be orthonormal and smoothness-ordered, which is exact regardless. Only the *interpretive* claim "these are the CAR modes" would need softening. |
| A2 | `AbstractGPs.jl` / `KernelFunctions.jl` / `GaussianRandomFields.jl` provide the described capabilities | Alternatives | None — all three are **rejected**; no plan depends on them. If overruled, each needs registry verification + a `Pkg.resolve()` dry-run. |
| A3 | Phase-11 training took ~18 epochs at 50k (back-derived from 16.08 min ÷ measured 0.557 s/epoch/512-sample MLP+D=8 cost) | Compute Budget | Medium. All epoch-scaled estimates move proportionally. Mitigate by measuring one real epoch before committing to a full run. |
| A4 | The Phase-11 pool's stored fields (`:summary_min`, `:theta`, `:lambda`) suffice for a read-only ε-identifiability ridge | Stage-1 gate | Low. Verified the field names exist in `run_p11_recovery.jl`; not executed. If wrong, the ε check needs its own small pool (minutes). |
| A5 | `chromatic_eps` is drawn independently of the spatial field and remains a per-image global | Pitfall 1 | Low — read directly from `prior.jl:103`. Would only change if Phase 12 alters the ε prior, which nothing asks for. |
| A6 | The observation-noise model needed by D-09 can be reduced to a single frozen `n_eff` per image size (Fisher-z) | D-09 construction | Medium. If the residual is strongly non-Gaussian on the z-scale, coverage will be biased for *both* arms equally, so the *comparison* survives but the absolute calibration claim does not. Check the z-residual normality on simulation before freezing. |
| A7 | Timings measured this session are representative | Compute Budget | Medium — the machine was concurrently running Phase 13. Treat all wall-clock figures as **upper-ish bounds with ±50 % uncertainty**, and re-measure one epoch before sizing the full run. |

---

## Open Questions

1. **What exactly is "the global term" under the copula? (D-01 vs D-05)**
   - *What we know:* D-01 says "the global term stays exactly the quantity the shipped estimator
     reports" (a scalar ρ_true). D-05's copula gives **every region its own** `MU_PRIOR` marginal —
     there is no scalar the field is a deviation *from*. D-05 explicitly **rejected** "global μ* plus a
     zero-mean deviation field" because it breaks the marginal.
   - *What's unclear:* whether the global θ row should be (a) the DCT DC coefficient `c₀` — a genuinely
     *sampled* parameter under the orthonormal rotation, D-07-clean, but not the shipped quantity; or
     (b) a read-time-derived scalar ρ (e.g. `ghat` of the pooled μ) reported for continuity but not a
     θ row.
   - *Recommendation:* **use `c₀` as the θ row and report the shipped-comparable scalar ρ as a
     derived read-time quantity**, stating clearly that it is derived. This satisfies D-07 exactly and
     satisfies D-01's *intent* (continuity of what is reported) rather than its literal wording.
     **This needs the user's confirmation — 12-CONTEXT delegated nothing to Claude's discretion.**

2. **Does the head parameterize the field in Gaussian space or ρ space?**
   - *What we know:* Gaussian space is atom-free (measured 6.79 % atoms per region in ρ space);
     `ghat ∘ F⁻¹ ∘ Φ` is monotone elementwise so credible intervals transform exactly; the deliverable
     (Δρ map) is in ρ units either way.
   - *What's unclear:* whether reporting SBC in Gaussian space is acceptable to the manuscript's
     framing, given the phase's headline is a Δρ map.
   - *Recommendation:* parameterize and score in Gaussian space; report *both* the Gaussian-space SBC
     (clean) and the randomized-rank ρ-space SBC (comparable to the milestone's existing tables).
     Costs one extra table.

3. **Per-row z-scoring vs a single shared scalar for the 64 continuous rows.**
   - *What we know:* the existing frozen `zt` is per-row; a CNN's weight sharing assumes translation
     equivariance, which per-row affine breaks.
   - *What's unclear:* whether it matters empirically at G = 8, where a radial chromatic gradient makes
     the per-position statistics genuinely different anyway.
   - *Recommendation:* keep per-row (matching the existing loader) as the default, add shared-scalar as
     a declared alternative arm only if the CNN underperforms the MLP control. Record which was used.

4. **How does the correlation length enter the D-10 ablation, mechanically?**
   - *What we know:* D-10 says "correlation length driven to zero". At ℓ → 0 (or α → 0) the lattice
     covariance → I, the regions become i.i.d. under the prior, and the copula still gives each region
     `MU_PRIOR` exactly.
   - *What's unclear:* whether the ablation retrains on an ℓ = 0 prior (a genuinely different training
     joint — the clean "matched ablation") or reuses the spatial net with ℓ pinned at 0 at read time
     (cheaper, but a mis-specified read).
   - *Recommendation:* **retrain**. D-10's whole point is "same network, same training, same summary,
     spatial prior neutralized" — a read-time pin is not that, and it would confound the gate with an
     out-of-distribution read. It also makes D-13's descope deliverable a real trained model.

5. **Does the D-12 Stage-1 threshold get set before or after the ridge?**
   - *What we know:* D-12 requires "a stated coverage improvement" pre-registered; Phase 11's Tier-1 /
     Tier-2 structure is the established way to have a measurement inform a threshold without
     data-snooping (judgment factor Tier-1, measured multiplicand Tier-2, floor written as a product).
   - *Recommendation:* mirror it exactly. Tier 1 in `p12_consts.jl` before anything runs; Tier 2
     appended (never edited) with one-line provenance per constant naming the artifact it came from.

6. **Is `chromatic_eps` identifiable?** Not previously measured, cheap to answer (ridge on the existing
   Phase-11 pool, ε excluded from nothing since there is no ε conditioning row). Because the ε effect
   clears the per-draw noise floor (ratio 1.272 at the prior edge), the answer is plausibly *yes* —
   which would be the project's first identified nuisance and would materially soften Pitfall 1. Worth
   a Wave-0 task.

---

## Environment Availability

| Dependency | Required by | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Julia | everything | ✓ | 1.12.6 (juliaup directory override, enforcing) | — |
| `spike/` environment | everything | ✓ | resolves; `NeuralEstimators 0.2.1`, `Flux 0.16.10` | — |
| `NeuralEstimators` v0.2.1 source | API verification | ✓ | `~/.julia/packages/NeuralEstimators/gFxuZ`, tree-sha1 `c780174a…` | — |
| `LinearAlgebra`, `SparseArrays`, `Statistics`, `Random`, `Printf`, `Dates` | lattice prior, ridge | ✓ | stdlib, via `@stdlib` on `LOAD_PATH` | — |
| CUDA / GPU | optional accelerator | ✗ (asserted absent by `spike/test/runtests.jl` D-04 check) | — | CPU-only; this is the *required* baseline, not a fallback |
| `test/test_images/{positive,negative}/*_c{1,2,3}.tif` | D-11 real-image arm | ✓ | committed, 1028×1376, exercised by `test/runtests.jl` | — |
| `corpus/` physical anchors | — | ✗ (`sha256 = PENDING-FETCH`, `bytes = 0`, `split = sealed_holdout`) | — | **Not a fallback — forbidden by D-11 AMENDED** |
| Phase-11 50k pool (`spike/data/cache/…`) | optional ε-identifiability ridge | likely ✓ (restored, verified by exact image-size counts in `11-DIAGNOSIS`) | — | Regenerate (~56 min) or draw a smaller fresh pool (~10 min) |
| `slopcheck` | package legitimacy | ✗ / not applicable | — | No package install is recommended; see the audit section |
| `ctx7` CLI / Context7 MCP | live API docs | ✗ (`command -v ctx7` → not found; no `mcp__context7__*` tools available) | — | **Used instead:** the pinned v0.2.1 source in the local depot (strictly higher confidence than docs) + `WebFetch` on the official docs site |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** Context7/ctx7 (fallback used and it is a *better* source);
`corpus/` anchors (deliberately excluded, not needed).

**Concurrency note (not a dependency, but it blocks):** Phase 13 is executing on the same branch.
`11-CLOSURE.md` records that concurrent agents race on the **git index**, and a Phase-12 planning
commit already swept two Phase-11 files into itself (`b82a02d`). **Every Phase-12 commit must be
pathspec-scoped.**

---

## Compute Budget

Anchored on Phase 11's measured baseline (`11-DIAGNOSIS`, "Rough compute estimate"): datagen
**56.37 min / 50k pairs**, training **16.08 min**, both CPU-only on 32 threads.

**What Phase 12 adds — measured this session:**

| Added cost | Measured | Over 50k pairs |
|---|---|---|
| Lattice draw + copula + `ghat` (n = 64, Cholesky precomputed) | **40.4 µs/draw** | **2.0 s** — negligible |
| Field → pixel bilinear upsample, **separable matmul** | 0.0003 s (512²) … 0.0273 s (2048²) | **~8 min** on the F5 mixture (realized weights 0.401/0.252/0.246/0.101 from `11-DIAGNOSIS`) |
| Field → pixel upsample, **naive scalar interpolant** | 0.3622 s at 1376×1028 | **~5 h** — do not do this |
| Per-pixel `a`/`b`/`sign` in stage 1 (≈4 extra elementwise passes over the image) | not measured | est. **+5–15 min** [ASSUMED] |
| **Datagen total** | | **≈ 65–80 min** (vs `P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150`) |

**Training cost — measured, and this is where the surprise is** (3 epochs at K = 512, extrapolated
to 50k × 18 epochs; see A3/A7 for the caveats):

| Summary net | Summary params | θ dim D | s / 3 epochs @ K=512 | ≈ min/epoch @ 50k | ≈ h @ 18 epochs |
|---|---|---|---|---|---|
| **fat CNN** (2→32→64→64, flatten 4096, Dense 4096→256) | 1 137 760 | 72 | 76.5 | 41.5 | **12.5** |
| **lean CNN** (2→16→32, 1×1 →8, flatten 512, Dense 512→128) | 70 872 | 72 | 11.1 | 6.0 | **1.8** |
| lean CNN | 70 872 | 24 | 2.6 | 1.4 | 0.4 |
| lean CNN | 70 872 | 9 | 2.4 | 1.3 | 0.4 |
| MLP control (shipped shape) | 197 504 | 24 | 8.9 | 4.8 | 1.5 |
| MLP control (shipped shape) | 197 504 | 8 | 2.2 | 1.2 | 0.4 |

[VERIFIED: measured 2026-07-27; MEDIUM confidence on absolute values — single-shot timings on a
machine concurrently running Phase 13. The **relative** figures are the robust part.]

**Two decision-relevant readings:**
1. **The fat CNN is out.** A `Dense(4096, 256)` after flattening 8·8·64 is 92 % of the summary net's
   parameters and costs **7×**. A 1×1 channel bottleneck before the flatten recovers essentially all
   of it. Specify the lean topology in the plan; do not leave "a CNN" underspecified.
2. **The θ width is the other lever:** D = 8 → 72 costs ~4–5×. That is D-04's high-rank mini-spike
   head, and it is affordable **only at reduced pool size**.

**Recommended sizing:**

| Activity | Pool | Config | Est. datagen | Est. training | Est. total |
|---|---|---|---|---|---|
| D-12 **Stage-1 ridge** (no net) | ~5k @ 512² | — | ~6 min | 0 | **~10 min** |
| D-03 **mini-spike**, per arm (CAR / GP / neutralized) | 10k @ 512² | lean CNN, D = 72 (K = 63) | ~13 min | ~22 min | **~35 min** |
| D-03 mini-spike, all 3 arms | | | | | **≈ 1.8 h** |
| **Full run** (chosen prior) | 50k, F5 mixture | lean CNN, D = 1+K+7+1 with K ≈ 15–24 | ~70 min | ~30 min | **~1.7 h** |
| SBC + coverage + guards | M = 2000, F5 | | ~30 min sim | ~15 min | **~45 min** |
| **Per full iteration** | | | | | **≈ 2.5–3.5 h** |

This gives D-12's two-stage trigger real numbers: **Stage 1 costs ~10 minutes and can kill the phase
before the ~1.8 h mini-spike; Stage 2 sits after a further ~2.5 h.** Declare a wall-clock ceiling in
`p12_consts.jl` the way `P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150` does, and **measure one real epoch
before committing to the full run** (A3/A7).

---

## Validation Architecture

`workflow.nyquist_validation` is `true` in `.planning/config.json`. [VERIFIED: config read]

### Test Framework

| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (`@testset`), invoked directly (no Pkg test target for the spike) |
| Config file | none — `spike/test/runtests.jl` is the single entry point |
| Quick run command | `julia --project=spike spike/test/runtests.jl` |
| Full suite command | `julia --project=spike -t auto spike/test/runtests.jl` (root package: `julia --project=. -e 'using Pkg; Pkg.test()'` — must stay green and must show `src/` untouched) |
| Reported-run commands | `julia --project=spike -t auto spike/validation/run_p12_*.jl` (each writes one `.jld2` artifact) |

**CRITICAL ordering constraint.** `spike/test/runtests.jl` ends with
`include("test_p13_correction.jl")`, whose outer `@testset` **throws** on two committed measured
misses; the file's own comment records that "a thrown testset aborts the remaining includes, so any
sibling placed after it would silently never run." **All Phase-12 `include`s must be inserted before
that line.** [VERIFIED: `spike/test/runtests.jl` tail]

### Phase Requirements → Test Map

No REQ-IDs exist for Phase 12 (see §Phase Requirements). The map below is keyed on the amended
success criteria and the decisions.

| Behaviour | Source | Test type | Automated command | File exists? |
|---|---|---|---|---|
| Tier-1 pre-registration constants are literal, seeds provably disjoint (incl. Phase-13 seeds) | carried-forward, D-12 | unit | `julia --project=spike spike/test/runtests.jl` (`test_p12_consts.jl`) | ❌ Wave 0 |
| CAR/GP kernels SPD; per-cell marginal sd = 1 after rescale | D-03, D-05 | unit | `test_p12_lattice.jl` | ❌ Wave 0 |
| DCT-II basis is orthonormal; `idct(dct(F)) == F` to float tolerance | D-04, D-07 | unit | `test_p12_lattice.jl` | ❌ Wave 0 |
| Prior parametrized by induced lag-1 correlation; r₁ monotone in the parameter | Pattern 1 | unit | `test_p12_lattice.jl` | ❌ Wave 0 |
| **SIM-02 per region**: max over 64 regions of W1(induced μ, `MU_PRIOR`) ≤ pre-declared tol | D-05 | integration | `run_p12_sim02.jl` → `test_p12_prior.jl` asserts the artifact | ❌ Wave 0 |
| Per-region `ghat` atom mass measured and recorded (not gated) | D-05, S-3 | integration | same artifact | ❌ Wave 0 |
| `reshape_summary` round-trips `encode_d01` exactly (no transpose) | D-02 | unit | `test_p12_architecture.jl` | ❌ Wave 0 |
| Mask rows are never z-scored; mask channel is channel 2 | D-02 | unit | `test_p12_architecture.jl` | ❌ Wave 0 |
| CNN summary net + wide `NormalisingFlow` builds, trains 1 epoch, `sampleposterior` returns D×N | D-01, D-02 | smoke | `test_p12_architecture.jl` | ❌ Wave 0 |
| Field at constant ρ reproduces the current simulator **bit-for-bit** (regression) | D-06 | unit | `test_p12_prior.jl` — mirrors the `P11_STAGE6_EXACT` golden-fixture pattern | ❌ Wave 0 |
| Separable upsample equals a reference interpolant to tolerance | D-06 | unit | `test_p12_prior.jl` | ❌ Wave 0 |
| **D-12 Stage-1**: leave-region-out ridge beats the per-region prior baseline, with positive control | D-12 | reported | `run_p12_stage1_ridge.jl` | ❌ Wave 1 |
| Correlation-length identifiability ridge (raw rows + Moran's I) vs prior-mean baseline | D-08, S-1 | reported | `run_p12_ell_ridge.jl` | ❌ Wave 1 |
| `chromatic_eps` identifiability ridge on the Phase-11 pool | Pitfall 1 | reported | `run_p12_eps_ridge.jl` | ❌ Wave 1 |
| SBC: global term + deviation-energy summary; randomized ranks; equivalence test on nuisances | D-01, D-07, Pitfalls 7–8 | reported | `run_p12_sbc.jl` | ❌ Wave 3 |
| **Leave-region-out predictive coverage** on simulation, spatial vs matched ablation, coverage **and** proper scoring rule | D-09, D-10, Pitfall 5 | reported | `run_p12_coverage.jl` | ❌ Wave 3 |
| Same on the 6 real TIFFs; per-region n reported; sealed holdout untouched | D-11 | reported | `run_p12_coverage.jl --real` | ❌ Wave 3 |
| **Radial-energy guard**: fraction of deviation-field energy that is radial | S-4 | reported | `run_p12_guards.jl` | ❌ Wave 3 |
| **ε = 0 ablation**: spatial advantage survives with chromatic held at 0 | S-4 | reported | `run_p12_guards.jl` | ❌ Wave 3 |
| **Offset-grid guard**: no degradation when the field lattice is offset half a cell | D-06 | reported | `run_p12_guards.jl` | ❌ Wave 3 |
| Decoupling: `src/`, `spike/Project.toml`, `spike/Manifest.toml` byte-unchanged; `corpus/` untouched | CLAUDE.md, D-11 | integration | `test_p12_decoupling.jl` (git-status/sha assertions, Phase-11 §F20c pattern) | ❌ Wave 0 |
| Root suite still green and `src/` provably untouched | CLAUDE.md | integration | `julia --project=. -e 'using Pkg; Pkg.test()'` | ✅ exists |

### Sampling Rate

- **Per task commit:** `julia --project=spike spike/test/runtests.jl` (unit + smoke tiers only;
  reported runners are never invoked from the suite).
- **Per wave merge:** `julia --project=spike -t auto spike/test/runtests.jl` **plus**
  `julia --project=. -e 'using Pkg; Pkg.test()'` (the decoupling proof).
- **Phase gate:** full spike suite green + every reported `.jld2` artifact present and re-derivable +
  the D-12 Stage-2 verdict recorded, before `/gsd:verify-work`.

### Wave 0 Gaps

- [ ] `spike/validation/p12_consts.jl` — Tier-1 pre-registration (fresh DEV seed, reserved counters,
      Stage-1/Stage-2 thresholds, `P12_ITERATION_ALLOWANCE = 1`, per-region W1 tolerance, jitter,
      radial-energy ceiling, offset-grid tolerance, wall-clock ceiling)
- [ ] `spike/test/test_p12_consts.jl` — literal re-assertion of every Tier-1 constant + executable
      seed-disjointness including **Phase 13's** seeds
- [ ] `spike/simulator/p12_lattice.jl` + `spike/test/test_p12_lattice.jl`
- [ ] `spike/simulator/p12_prior.jl` + `spike/test/test_p12_prior.jl` (incl. the constant-ρ
      bit-for-bit golden fixture, captured **before** the first simulator edit)
- [ ] `spike/npe/p12_architecture.jl` + `spike/test/test_p12_architecture.jl`
- [ ] `spike/test/test_p12_decoupling.jl`
- [ ] Wire all Phase-12 includes into `spike/test/runtests.jl` **before** `test_p13_correction.jl`

*(No framework install is needed — Julia stdlib `Test` is already the harness.)*

---

## Security Domain

`security_enforcement` is not set in `.planning/config.json` (absent ⇒ enabled), so the section is
included. This is an offline, local, scientific-computing phase with no network, no user input surface,
no authentication and no persistence of untrusted data.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard control |
|---------------|---------|------------------|
| V2 Authentication | no | No auth surface; no service |
| V3 Session Management | no | No sessions |
| V4 Access Control | **yes, in a project-specific sense** | The `corpus/` sealed holdout is access-controlled by convention (`open_sealed_holdout(; reason)`); D-11 forbids consumption. Enforce with an executable assertion, not a comment |
| V5 Input Validation | **yes** | The simulator already validates θ at entry (`forward.jl:122-158`, `ArgumentError` on non-finite / out-of-range / `chromatic_eps ≤ −1`). Extend the same discipline to the field: assert the ρ-field is finite and in [−1, 1] elementwise, and that the correlation-length parameter is inside its declared box |
| V6 Cryptography | no (but see below) | No secrets. SHA-256 is used for **content addressing** (cache dirs, artifact integrity), never for security. Do not hand-roll — reuse the Phase-3 `hashguard.jl` pattern |
| V12 File Handling | **yes** | Artifact reads/writes use the established atomic `.tmp` → reopen-integrity-check → `mv(force=true)` wrapper. Deserialization surface: JLD2 loads of *this project's own* artifacts only (T-7-01's trusted-input framing) |

### Known Threat Patterns for this stack

| Pattern | STRIDE | Standard mitigation |
|---|---|---|
| Untrusted JLD2 deserialization | Tampering / EoP | Load only in-repo, self-produced artifacts; verify integrity on reopen (existing wrapper). Never load a downloaded `.jld2` without the content-hash check `_lazy_load_from_artifact!` already performs |
| Silent corruption of a frozen pre-registration | Tampering (scientific integrity) | Append-never-edit constants + literal re-assertion in a test + git history as the audit trail (`p11_consts.jl` pattern) |
| Consuming a reserved seed stream | Tampering (scientific integrity) | Executable `_forbidden()` assertion with **derived** gate seeds recomputed, not trusted from comments |
| Burning the sealed holdout | Information disclosure (of the blind test set) | D-11: an executable assertion that no `corpus/` path is opened by any Phase-12 code path |
| Concurrent agents racing the git index | Tampering | Pathspec-scoped commits (`git commit -- <paths>`); recorded as a structural item in `deferred-items.md` |
| Non-finite propagation crashing a reported run | DoS (of the run) | The `LOCAL_MAP_SENTINEL` + forced-OOD discipline: a degenerate region must be finite-sentinelled and flagged, never crash and never masquerade as a measurement (T-7-04) |

---

## Sources

### Primary (HIGH confidence)

- **`~/.julia/packages/NeuralEstimators/gFxuZ/`** — the *installed, pinned* v0.2.1 source
  (`Project.toml: version = "0.2.1"`; `spike/Manifest.toml` tree-sha1 `c780174a…`). Read:
  `src/Estimators/PosteriorEstimator.jl` (constructors, docstring example),
  `src/ApproximateDistributions/NormalisingFlow.jl` (constructor + `CouplingLayer`),
  `src/ApproximateDistributions/ApproximateDistributions.jl:28-30` (keyword `num_summaries` form),
  `src/Estimators/Estimators.jl:76-89` (`_summarystatistics` calls the network directly),
  `src/train.jl:6-8, 164-230, 661-671` (`train` signatures, `_dataloader`, `MLUtils.DataLoader`),
  `src/inference.jl:56-72` (`sampleposterior`), `src/NeuralEstimators.jl:40-98` (full export list;
  **no `Unet`**).
- **Live execution in the spike environment** (`julia --project=spike`, Julia 1.12.6, 2026-07-27):
  CAR/GP spectra and marginal-sd tables; copula marginal preservation and per-region atom mass;
  per-draw timings; separable-vs-scalar upsample timings; per-region noise floor and the chromatic
  radial decomposition; the CNN + wide-flow API smoke and the six-config training benchmark.
- **Repository sources read directly:** `spike/simulator/{prior,forward,ghat}.jl`,
  `src/amortized/{summary,architecture,local_map}.jl`, `src/results.jl:150-192`,
  `spike/validation/{consts,p11_consts}.jl`, `spike/validation/run_p11_recovery.jl`,
  `spike/npe/{architecture,p11_architecture}.jl`, `spike/test/runtests.jl`, `spike/Project.toml`.
- **`spike/validation/p11_probe_report.jld2`** — Phase 11's frozen, pre-registered probe artifact
  (`generated = 2026-07-25T21:26:36.787Z`); the `eps_*` and `shift_*` ladders quoted in Pitfall 1.
- **`.planning/phases/11-…/11-DIAGNOSIS.md`, `11-CLOSURE.md`** — the noise floor 0.9765, the ridge
  harness design, the flat-λ RMSE result, the coverage 0.72–0.77, the 56.37/16.08 min baseline.
- **Skill `spike-findings-proteincoloc` → `references/sbc-calibration.md`** — randomized-rank recipe,
  M = 2000 over-power (0.059 SD at 50 % power), spikes 012/013 verdicts.
- **`git`** — `ca02b0e8a696269a2b8710b11ca3073edc5ae6c7` confirmed an ancestor of `HEAD`
  (`cecb69c58770a5b45bec8151c657d9024b999f80`) via `git merge-base --is-ancestor`.

### Secondary (MEDIUM confidence)

- **https://msainsburydale.github.io/NeuralEstimators.jl/dev/** — navigation/API page index (fetched);
  **https://msainsburydale.github.io/NeuralEstimators.jl/dev/examples/data_gridded_nonstationary** —
  the spatially-varying-parameters example (I2I `Unet` + `PointEstimator`, 4-D `grid×grid×1×K` data).
  MEDIUM because the site tracks `dev`, not the pinned v0.2.1 — corroborated against the local source,
  which is why the `Unet` absence is stated as a v0.2.1 fact and the docs example as a `dev` pattern.
- **Compute extrapolations** — measured 3-epoch timings scaled to 50k × 18 epochs; the 18-epoch figure
  is back-derived from Phase 11's 16.08 min (A3), and the machine was concurrently busy (A7).

### Tertiary (LOW confidence — flagged, nothing depends on them)

- **`AbstractGPs.jl` / `KernelFunctions.jl` / `GaussianRandomFields.jl` capability claims** — training
  knowledge only, not verified this session. All three are **rejected**, so no plan depends on them.
- **The DCT-II ↔ lattice-Laplacian eigenbasis identity** — standard spectral-graph result recalled from
  training, not re-derived. Only the interpretive framing depends on it (A1).

### Not consulted, and why

- **Context7 / `ctx7` CLI** — `command -v ctx7` returned not-found and no `mcp__context7__*` tools are
  available in this agent's toolset. The documented fallback (official docs via `WebFetch`) was used,
  **and** superseded by reading the pinned source in the local depot, which is strictly more
  authoritative for a version-pinned question. CLAUDE.md's instruction to "verify signatures against
  `dev` docs, not memory" is satisfied more strongly than it asks.
- **Knowledge graph** (`.planning/graphs/graph.json`) — queried; **stale by 706 hours / 453 commits**
  (`built_at_commit e06414e` vs `HEAD cecb69c`) and returned zero nodes for "spatial lattice
  correlation grid". Treated as unavailable; all structural knowledge in this document comes from
  direct file reads.

---

## Metadata

**Confidence breakdown:**

| Area | Level | Reason |
|------|-------|--------|
| NeuralEstimators v0.2.1 API (head, CNN summary, 4-D data, training, sampling) | **HIGH** | Read the pinned source AND ran it end-to-end this session |
| Lattice prior construction, CAR/GP spectra, marginal rescale necessity | **HIGH** | Measured directly; the numbers are reproducible from the tables |
| D-05 copula marginal preservation and per-region atom mass | **HIGH** | Measured at M = 20 000; W1 max 0.00732 vs tol 0.10 |
| Chromatic radial confound (S-4) | **HIGH** | Two independent measurements agree (Phase-11 frozen probe artifact + fresh paired simulation) |
| Per-region noise floor | **MEDIUM-HIGH** | Measured, but n = 6 base draws at one image size; the |ρ| dependence is visible but not characterized |
| Compute estimates | **MEDIUM** | Single-shot timings on a concurrently-busy machine; relative ratios robust, absolutes ±50 % |
| D-09 predictive-coverage construction | **MEDIUM** | No prior art in this repo; the observation-noise model is a declared assumption (A6) and the discrimination risk (Pitfall 5) is a genuine design hazard |
| DCT-II basis choice | **MEDIUM-HIGH** | Orthonormality and invertibility are exact and testable; the "these are the CAR modes" interpretation is approximate for α < 1 and is stated as such |
| Descope-trigger sizing (D-12) | **MEDIUM** | Follows from the compute estimates, inherits their uncertainty |

**Research date:** 2026-07-27
**Valid until:** ~2026-08-26 for the API and stack findings (the environment is version-pinned and
frozen, so these do not decay); **re-verify the compute numbers before sizing any full run** — they
were measured under concurrent load.
