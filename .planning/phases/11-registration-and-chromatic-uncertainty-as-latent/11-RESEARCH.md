# Phase 11: Registration and Chromatic Uncertainty as Latent — Research

**Researched:** 2026-07-25
**Domain:** Simulator surgery (affine warp composition), amortized SBI conditioning, equivalence-test coverage validation, provenance/content-hash discipline
**Confidence:** HIGH for the repo-mechanics findings (all measured or code-cited); MEDIUM for the training-budget and threshold-derivation proposals (arithmetic is sound, but the numbers depend on a probe that has not yet run under pre-registration).

---

<user_constraints>
## User Constraints (from 11-CONTEXT.md)

**These are LOCKED. This research answers HOW, never WHETHER.**

### Locked Decisions (verbatim from `11-CONTEXT.md` `<decisions>`)

- **D-01: Research lane only.** Train a Phase-11 research net in `spike/` on a **fresh DEV seed**
  asserted disjoint from `PROD_SEED_V2`, `VAL_MASTER_SEED` (`0x5BC0FFEE`), `NPE_MASTER_SEED`
  (`0xC0FFEE`), `VAL_FIX_SEED` (`0xF1F7ED`), `CORPUS_MASTER_SEED` (`0x0000000000c05eed`) and all
  prior dev seeds. The shipped `amended_v2/grid_8` artifact and the Phase-7 GO (Option A,
  2026-07-24) stay untouched — no reship, no gate reopened. Science advances without spending
  pre-registration credibility.

- **D-02: Widen `SHIFT_PRIOR` to `Uniform(-3, 3)` px** for the research net (from `Uniform(-1, 1)`).
  Covers realistic multi-channel error — chromatic aberration plus stage/filter-cube drift — and
  gives the SC2 ladder real headroom (0.25 → 3 px). **A wider nuisance prior adds marginal spread
  that may attenuate `ρ_true` inference; this must be MEASURED against the frozen net, not assumed.**

- **D-03: Amortize over the uncertainty level.** Append the registration-uncertainty level (the
  shift prior's half-width) to the network input so one training run yields a continuous, controlled
  sweep — feed the same image at level 0.25 → 3 px and read posterior width directly. **Consequence:
  the input dimension changes, so this net is NOT drop-in comparable to the shipped read surface.**

- **D-04: Pre-register, allow one iteration.** Lock the SC2 monotonicity and SC3 no-overconfidence
  thresholds in a Phase-11 consts file (pattern: `spike/validation/consts.jl`) with the fresh DEV
  seed **BEFORE any run**, and state up front that this is a research-lane result permitting **one
  documented iteration**.

- **D-05: Marginalize only — no identification claim.** The patch-correlation summary stays
  unchanged. `gate-8x8-amended.md:42` lists
  `vacuous_params = [spillover, autofluorescence, shift_dx, shift_dy, noise]` at shrinkage 0.92–1.03
  — the posterior IS the prior. Widening the shift prior propagates registration uncertainty into
  `ρ_true`/`Δρ` as honest marginal spread; the shift columns (and `ε`) remain vacuous and are
  **reported as such**. Per `gate-8x8-amended.md:46`, a clean KS p on a vacuous column is
  meaningless and must not be quoted as calibration evidence.
  **SC2 is a claim about the COLOC posterior widening, not about recovering dx/dy/ε.**

- **D-06: Simulator-only pre-flight probe BEFORE training.** Probe = sweep shift `0 → 3` px and `ε`
  through `forward.jl` at fixed θ and measure how far the 128-dim summary actually moves, relative
  to the movement induced by a known `Δρ`. Minutes, no training. **The probe supplies the SC2
  threshold; the threshold is not guessed.**

- **D-07: SC2 pass criterion is two-part.** (a) Posterior SD (or 90% HDI width) of `Δρ` is
  **monotone non-decreasing** across the uncertainty ladder, AND (b) empirical coverage of the 90%
  credible interval stays within nominal at **every** rung.

- **D-08: Per-rung coverage tested by equivalence, not by a point null.** TOST-style test for
  equivalence to nominal with a **pre-registered tolerance band** (e.g. ±3 percentage points on the
  90% CI), **Holm-corrected across rungs**. **The tolerance band must be justified from the D-06
  probe, not chosen after seeing results.**

- **D-09: Add a 1-parameter radial chromatic scale `ε`,** prior `Uniform(-0.02, 0.02)`, symmetric.
  `ε` joins θ with its own entry in `prior.jl`. Note: `ε` is expected to be another **vacuous**
  column (D-05 applies to it too).

- **D-10: Single composed affine warp.** Build one `AffineMap = Translation(dy, dx) ∘ LinearMap(
  scale about image centre)` and pass it to a **single** `warp` call, replacing the current
  `Translation`-only call at `forward.jl:168-174`. One interpolation pass keeps the `ε = 0,
  shift = 0` case byte-comparable with today's simulator. The `WR-06` axis-order comment must be
  **re-derived** for the composed map.

- **D-11: Mirror the `ε` extension into `src/amortized/simulator.jl` now** (not spike-only).

- **D-12: Record the provenance consequence of D-11 as a named limit** in `docs/amortized.md`,
  pinning the pre-`ε` commit sha. Provenance is repaired **by reference rather than by bytes**.

- **D-13: SC3 is a coverage claim, not a flagging claim.** No OOD flag requirement. **Structural
  tension:** the wider the prior, the less mis-registration is detectable as OOD — SC2's mechanism
  directly erodes any density-channel detector. **Accepted consequence:** badly mis-registered data
  returns a wide-but-unflagged answer.

- **D-14: Report the breakdown curve.** Test coverage well past the training prior (0 → 8 px,
  `|ε| → 0.05`) and report the point where 90% coverage first drops below nominal. **SC3 passes if
  coverage holds throughout the training prior; the breakdown point beyond it is REPORTED, not
  gated.**

- **D-16: Deliverable = report + a `docs/amortized.md` interpretation note.** **No API change, no
  artifact change.**

- **D-17: Validation is simulation + a qualitative real-image check.** Deliberately shift one
  channel by a known amount and confirm the posterior widens as the simulated sweep predicts.
  **This cannot be a coverage test.** **Read-only access to those images; the manuscript pipelines
  must be provably unmodified, and the report must say so.**

### Claude's Discretion

- **D-15: Misalignment for SC3 is simulator-injected, not post-hoc.** Force `shift`/`ε` to the test
  value and run the single composed affine warp (D-10); do **not** post-hoc warp already-simulated
  images.

### Deferred Ideas (OUT OF SCOPE — do not plan)

- Summary extension for registration identifiability (lag / cross-correlation-peak features)
- A dedicated registration OOD channel
- Staged marginalize-then-identify comparison
- Off-axis optical centre for the chromatic warp
- Dedicated bead / registration-target acquisition
- API-level width advisory
- Productionizing the Phase-11 research net
</user_constraints>

---

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| — | **No requirement IDs are mapped.** `.planning/ROADMAP.md:235` records `**Requirements**: TBD` for Phase 11, and a grep of `.planning/REQUIREMENTS.md` finds no `REG-*` / `CHROM-*` / Phase-11 identifiers. | The phase is governed by the three ROADMAP Success Criteria (`ROADMAP.md:237-239`) as refined by the 17 locked decisions. The planner should either (a) plan against SC1/SC2/SC3 directly, or (b) mint Phase-11 requirement IDs into `REQUIREMENTS.md` as a Wave-0 task. This research assumes (a). |

**ROADMAP Success Criteria (`.planning/ROADMAP.md:237-239`), as reframed by `11-CONTEXT.md`:**

- **SC1** — dx/dy (+ 1-param chromatic warp) in θ, NPE retrained. **dx/dy portion ALREADY SATISFIED**
  (`spike/simulator/prior.jl:56,78-79`; `spike/simulator/forward.jl:172-174`;
  `src/amortized/simulator.jl:94,116-117,236-238`; read back at `src/amortized/ood.jl:129-130`).
  Phase-11 delta = the `ε` term + the research-lane retrain.
- **SC2** — posterior width increases monotonically with injected registration uncertainty. **Never run.**
- **SC3** — deliberately mis-registered images handled without silent overconfidence. **Unaddressed.**
</phase_requirements>

---

## Summary

Phase 11's real technical content is far smaller than the ROADMAP text implies and far riskier in
exactly one place. The `shift_dx`/`shift_dy` latents already exist end-to-end, so SC1 reduces to
adding one θ column (`ε`) plus a research-lane retrain. The genuinely hard parts are (i) the
composed-affine surgery in stage 6 of two mirrored simulators without perturbing existing behaviour,
(ii) the θ-arity ripple from D-11 into `src/`, which touches a content-hashed file, a hard-coded
`theta_dim = 7`, a hard-coded `NPE_D = 7`, and roughly a dozen literal `7`s in `test/runtests.jl`,
and (iii) designing a per-rung equivalence test whose thresholds are genuinely probe-derived rather
than reverse-engineered.

**The highest-risk unknown named in the brief — "does the 8th θ column break the shipped bundle
load?" — is answered NO, definitively, with code evidence.** `src/amortized/persist.jl:97` rebuilds
the estimator with `build_estimator(a.d_in, a.D; …)` where `a` is the **architecture metadata read
off disk**, not the module constant; and the θ un-standardizer `θzt` is a persisted
`BoundedThetaTransform` carrying its own 7-element `lo`/`hi` vectors (`architecture.jl:100-104`),
so `theta_prior_bounds()` is never re-consulted at load time. The shipped `Artifacts.toml`
`git-tree-sha1` pin is likewise independent of the datagen digest. What D-11 *does* break is the
`src/` **training** path (`datagen.jl:170,274,393` hard-code `7`; `train_npe.jl:136` passes no `D`,
so `build_estimator` silently defaults to `NPE_D = 7`) and roughly a dozen assertions in
`test/runtests.jl`. Those are loud, fixable failures, not silent corruption — but they are real
work the CONTEXT does not enumerate.

I ran an **indicative** (not pre-registered) simulator-only probe and a **working prototype of the
D-10 composed warp**. Two results de-risk the phase substantially. First, the composed
`Translation(dy,dx) ∘ recenter(LinearMap(s·I), c)` at `ε = 0` produces output that is **bit-for-bit
identical** to today's `simulate_pair` — verified through the whole seven-stage pipeline, not just
the `warp` call, for four (dy,dx) pairs including `(0,0)`. D-10's regression check is achievable at
exact `==`, not merely `isapprox`. Second, the summary is **not** below resolution: at 256², a 3 px
misalignment moves the 64 continuous summary rows by as much as a `Δρ` of ≈ 0.19, and `ε = 0.02`
moves it by as much as `Δρ ≈ 0.10` (≈ 0.34 at the 1376×1028 anchor). The D-06 "SC2 fails flat"
risk looks low. But the probe also exposed a design trap: the between-seed Monte-Carlo floor
(‖Δs‖ ≈ 2.65) is **larger than the entire shift effect at 3 px**, so an unpaired probe would be
swamped — the probe *must* be paired on the RNG key.

**Primary recommendation:** Sequence the phase as `pre-registration consts → capture pre-ε sha +
golden fixture → composed-warp surgery (spike + src + docs limit + test updates, one commit) →
pre-flight probe → freeze probe-derived thresholds → train research net → SC2 ladder → SC3 breakdown
→ real-image check → report + docs note`. Train in `spike/` with a bespoke trainer (not
`_train_grid_pipeline`, which cannot express `D = 8` or a 129-row input), and append the
uncertainty level λ as a 129th input row **after** the frozen 128-row standardization.

---

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Chromatic `ε` forward physics | Simulator (`spike/simulator/forward.jl` stage 6) | Mirror in `src/amortized/simulator.jl` | D-10/D-11; stage 6 already owns geometric misalignment |
| `ε` prior + support box | Prior (`spike/simulator/prior.jl`, `src/amortized/simulator.jl:140-148`) | — | `theta_prior_bounds()` is the declared single source of truth (`simulator.jl:127`) |
| Uncertainty-level λ conditioning | Spike training/inference only | — | D-03 explicitly must not reach the shipped read surface |
| Coverage / monotonicity statistics | Spike validation (`spike/validation/`) | Reuse `test/gate/sbc.jl` reductions | D-07/D-08 say EXTEND, not fork |
| Provenance / named limit | `docs/amortized.md` §Named limits | git (pinned sha) | D-12; repaired by reference |
| Shipped read surface | `src/amortized/{api,registry,persist}.jl` | — | **Must be byte-unchanged in behaviour** (D-16) |
| Real-image check | `test/test_images/` via `src/LoadImages.jl` (read-only) | — | D-17; see §F for why NOT `corpus/` |

---

## Standard Stack

Everything Phase 11 needs is already resolved in both manifests. **No new dependency is required.**

### Core (already present, versions verified from the resolved manifests)

| Library | Version | Purpose | Evidence |
|---------|---------|---------|----------|
| `CoordinateTransformations` | **0.6.4** | `Translation`, `LinearMap`, `recenter`, `∘` | `Manifest.toml` and `spike/Manifest.toml` both pin 0.6.4 [VERIFIED: local manifest read] |
| `ImageTransformations` | **0.10.3** | `warp`, `center` | both manifests [VERIFIED] |
| `Interpolations` | **0.16.3** | `BSpline(Linear())` | both manifests [VERIFIED] |
| `ImageFiltering` | **0.7.12** | PSF `imfilter`/`Kernel.gaussian` (stages 1,3 — untouched) | both manifests [VERIFIED] |
| `StaticArrays` | **1.9.18** | transitive; `SVector`/`SMatrix` used inside `warp` | both manifests [VERIFIED] |
| `NeuralEstimators` | **0.2.1** | `PosteriorEstimator`, `NormalisingFlow`, `train`, `sampleposterior` | `~/.julia/packages/NeuralEstimators/gFxuZ/Project.toml:3` [VERIFIED] |
| `Flux` | 0.16.10 | NN backend | root `Project.toml:38` compat [VERIFIED] |
| `Random123` | — | `Philox4x` keyed RNG | `spike/Project.toml:16` [VERIFIED] |
| `HypothesisTests` | — | KS/χ² (existing SBC arm); **not** needed for TOST if hand-rolled | `spike/Project.toml:9` [VERIFIED] |

### Explicitly NOT needed

| Considered | Verdict | Why |
|------------|---------|-----|
| `MultipleTesting.jl` | **Do not add** | Holm is already hand-rolled twice in-repo: `test/gate/sbc.jl:352-362` (`sbc_holm_adjusted`) and `test/gate/gate_consts_8_v2.jl:361-371` (`holm_adjusted`), with a unit test asserting they agree (`sbc.jl:349-350`). Reuse. Neither is in `spike/Project.toml`, so a spike-local copy or a `include` of the gate file is required. |
| `LinearAlgebra` (`I`) in `spike/` | **Avoid** | `spike/Project.toml` does not list it, so `using LinearAlgebra` would require editing spike deps. Use a literal `[s 0.0; 0.0 s]` `Matrix{Float64}` instead — `ImageTransformations.try_static` (`autorange.jl:99-100`) converts a plain-matrix `AffineMap` to `SMatrix`/`SVector` automatically. [VERIFIED: source read + measured equality] |
| `Rotations.jl` | Not needed | A radial scale is a diagonal `LinearMap`; no rotation group required. |
| `Isotonic` / monotone-regression package | Not needed | The recommended SC2 statistic (§E17) is a Spearman + permutation test, computable with `StatsBase.corspearman` (already a dep). |

---

## Package Legitimacy Audit

**Not applicable — Phase 11 installs no external packages.** Every library named above is already
resolved and pinned in the committed `Manifest.toml` / `spike/Manifest.toml`. The only dependency
change even contemplated (`LinearAlgebra` into `spike/Project.toml`) is explicitly avoided by the
recommended design.

If the planner nonetheless adds a dependency, the Package Legitimacy Gate must run first; the Julia
ecosystem equivalent of the registry check is `Pkg.Registry` resolution plus a General-registry
`Package.toml` `repo` field check.

---

## A. Simulator Surgery (D-09, D-10, D-11, D-15)

### A1. What stage 6 does today — and are the two copies identical?

**`spike/simulator/forward.jl:168-174` (verbatim):**

```julia
    # --- (6) sub-pixel registration shift on channel 2 only (fillvalue > 0) ------
    # WR-06: Translation(a, b) shifts the FIRST array axis (rows = vertical = dy) by a
    # and the SECOND (columns = horizontal = dx) by b. Pass (dy, dx) so the named
    # fields map to their conventional physical axes (dx horizontal, dy vertical).
    shifted = warp(ch2, Translation(θ.shift_dy, θ.shift_dx), axes(ch2);
                   method = BSpline(Linear()), fillvalue = BG_FLOOR)
    ch2 = Matrix{Float64}(collect(shifted))
```

**`src/amortized/simulator.jl:234-238` (verbatim):**

```julia
    # --- (6) sub-pixel registration shift on channel 2 only (fillvalue > 0) ------
    # WR-06: Translation(dy, dx) — first array axis (rows = vertical = dy), second (cols = dx).
    shifted = warp(ch2, Translation(θ.shift_dy, θ.shift_dx), axes(ch2);
                   method = BSpline(Linear()), fillvalue = BG_FLOOR)
    ch2 = Matrix{Float64}(collect(shifted))
```

**Behaviourally identical: YES. Byte-identical: NO** — the executable lines are character-for-character
the same (`warp(ch2, Translation(θ.shift_dy, θ.shift_dx), axes(ch2); method = BSpline(Linear()),
fillvalue = BG_FLOOR)` and the `collect` line), but the WR-06 comment above them is abridged in
`src/`. The imports differ in style only: `forward.jl:49-51` uses `using ImageTransformations` /
`using CoordinateTransformations` / `using Interpolations`; `simulator.jl:46-48` uses
`import ImageTransformations: warp` / `import CoordinateTransformations: Translation` /
`import Interpolations: BSpline, Linear`. `BG_FLOOR = 0.02` and every other stage constant match
(`forward.jl:57-66` ↔ `simulator.jl:152-156`).

**Correction the plan must absorb before re-deriving WR-06.** `warp` is **backward-mode**. From the
installed source, `ImageTransformations/ET6hs/src/warp.jl:4`:

> `Transform the coordinates of img, returning a new imgw satisfying imgw[I] = img[tform(I)].`

and `warp.jl:167-176` implements exactly `out[I] = _getindex(img, tform(SVector(I.I)))`. So `tform`
maps **destination index → source index**, and `Translation(dy, dx)` moves image *content* by
`−(dy, dx)`. The existing WR-06 comment describes the *argument-slot* convention correctly (first
slot = rows = dy, second = cols = dx) but its "shifts … by a" phrasing implies forward mode. Because
`SHIFT_PRIOR` is symmetric about 0 this is observationally irrelevant for the marginal law of the
simulated data — but the re-derived comment must state the backward-mode fact explicitly, because
**it determines the sign of the `ε` scale factor** (see A2).

### A2. The exact composed-affine expression

**Recommended, verified-working form** (uses only symbols already imported or importable from the
pinned deps; introduces no new dependency):

```julia
    # --- (6) ONE composed affine warp: radial chromatic scale then registration shift ---------
    # WR-06 (RE-DERIVED for the composed map, D-10). `warp` is BACKWARD-mode:
    #   out[I] = ch2[A(I)]   (ImageTransformations warp.jl:4, :167-176)
    # so `A` maps a DESTINATION index to the SOURCE index it samples, and the image CONTENT
    # undergoes A⁻¹. Two consequences:
    #   • Translation slot order is unchanged: slot 1 = first array axis (rows = vertical = dy),
    #     slot 2 = second axis (cols = horizontal = dx). Pass (dy, dx) exactly as before.
    #   • To MAGNIFY channel 2's content by (1 + ε) about the image centre, the backward map must
    #     SHRINK coordinates by that factor: s = 1 / (1 + ε).
    # `recenter(t, c) == Translation(c) ∘ t ∘ Translation(-c)` (CoordinateTransformations
    # core.jl:103-105), so the scale is applied about `c`, the geometric centre of `axes(ch2)`.
    # The whole thing collapses to ONE `AffineMap` (compose methods, affine.jl:139-165), so `warp`
    # performs exactly ONE interpolation pass — the D-10 correctness requirement.
    s = 1.0 / (1.0 + θ.chromatic_eps)
    c = map(ax -> (first(ax) + last(ax)) / 2, axes(ch2))     # (128.5, 128.5) at 256²
    A = Translation(θ.shift_dy, θ.shift_dx) ∘ recenter(LinearMap([s 0.0; 0.0 s]), c)
    shifted = warp(ch2, A, axes(ch2); method = BSpline(Linear()), fillvalue = BG_FLOOR)
    ch2 = Matrix{Float64}(collect(shifted))
```

**Symbol provenance (all pinned, both environments):**

| Symbol | Package | Version | Source line |
|--------|---------|---------|-------------|
| `Translation` | CoordinateTransformations | 0.6.4 | `src/affine.jl:17-22` |
| `LinearMap` | CoordinateTransformations | 0.6.4 | `src/affine.jl:50-52` |
| `recenter` | CoordinateTransformations | 0.6.4 | `src/core.jl:103-106` (exported, `CoordinateTransformations.jl:7`) |
| `∘` (compose) | Base, methods added by CoordinateTransformations | 0.6.4 | `src/core.jl:56` (`const compose = ∘`), affine methods `src/affine.jl:139-165` |
| `warp` | ImageTransformations | 0.10.3 | `src/warp.jl:162-183` |
| `BSpline`, `Linear` | Interpolations | 0.16.3 | (unchanged from today's call) |
| `try_static` | ImageTransformations (internal) | 0.10.3 | `src/autorange.jl:97-105` — converts plain-`Matrix` `AffineMap` → `SMatrix`/`SVector` |

**Import deltas required.** `spike/simulator/forward.jl` already does `using CoordinateTransformations`
(line 50), which exports `LinearMap` and `recenter` — **no import change needed in the spike**.
`src/amortized/simulator.jl:47` does `import CoordinateTransformations: Translation` and must become
`import CoordinateTransformations: Translation, LinearMap, recenter`.

**Why `c` is a plain `Tuple`, not `ImageTransformations.center`.** `recenter(trans, origin::Tuple)`
is a defined method (`core.jl:106`) that converts to `SVector` internally, so no `StaticArrays`
import is needed. `map(ax -> (first(ax)+last(ax))/2, axes(ch2))` reproduces
`ImageTransformations.center` exactly (`ImageTransformations.jl:80-81`:
`center(img) = SVector{N}(map(_center, axes(img)))`, `_center(ind) = (first(ind)+last(ind))/2`).
Avoiding the `center` import also avoids any name collision inside the `ProteinCoLoc` module.

**Measured composition result.** At `ε = 0.02`, `axes = (1:256, 1:256)`, `shift = 0`:

```
A1 = AffineMap([0.9803921568627451 0.0; 0.0 0.9803921568627451], [2.5196078431372655, 2.5196078431372655])
corner (1,1) maps to source [3.5000000000000107, 3.5000000000000107]  displacement = 3.5355 px
```

which is exactly `|c − x| · ε/(1+ε) = 127.5·√2 · 0.02/1.02 = 3.5355`. This matches the CONTEXT's
D-09 magnitude table (`~3.6 px` at ε = 0.02 on a 256² field) to within the `1/(1+ε)` factor.
[VERIFIED: executed `julia` against `spike/Project.toml`, 2026-07-25]

### A3. Proving the `ε = 0, shift = 0` regression check — **exact equality is achievable**

**Result: `==` holds, bit for bit, through the entire `simulate_pair` pipeline.** I built a
prototype of the composed stage 6 and compared its 64 continuous summary rows against the current
`simulate_pair` for 12 independent Philox keys:

```
eps=0,shift=0 composed == legacy simulate_pair : true
```

and compared the raw warped matrices directly for four shift pairs:

```
dy=0.0  dx=0.0  : exact_equal=true  maxabsdiff=0.0
dy=0.37 dx=-0.82: exact_equal=true  maxabsdiff=0.0
dy=2.5  dx=1.25 : exact_equal=true  maxabsdiff=0.0
dy=-3.0 dx=3.0  : exact_equal=true  maxabsdiff=0.0
```

[VERIFIED: executed 2026-07-25, `CoordinateTransformations 0.6.4` / `ImageTransformations 0.10.3` /
Julia 1.12.6]

**Why it is exact, not lucky.** Tracing the compose methods:

1. `recenter(LinearMap(M), c)` = `Translation(c) ∘ LinearMap(M) ∘ Translation(-c)`.
   Inner: `compose(::LinearMap, ::Translation) = AffineMap(M, M*(-c))` (`affine.jl:143-145`). At
   `M = [1 0; 0 1]`, `M*(-c) == -c` exactly (`1.0*(-c₁) + 0.0*(-c₂)`).
   Outer: `compose(::Translation, ::AffineMap) = AffineMap(t2.linear, t1.translation + t2.translation)`
   (`affine.jl:163-165`) ⇒ `AffineMap(M, c + (-c)) = AffineMap(M, [0.0, 0.0])` — exact, since
   `x - x == 0.0` for any finite `x`.
2. `Translation(dy,dx) ∘ AffineMap(M, 0)` ⇒ `AffineMap(M, [dy, dx])`.
3. `try_static` promotes to `AffineMap(SMatrix{2,2}(I), SVector(dy,dx))`.
4. Evaluation `(::AffineMap)(x)` = `M*x + v` (`affine.jl:110-114`). With `M` the exact identity and
   `x::SVector{2,Int}`, `1.0*x₁ + 0.0*x₂ == Float64(x₁)` exactly (and identically under FMA), so the
   result equals `Translation(dy,dx)(x) = x + SVector(dy,dx)` bit for bit.

Identical source coordinates ⇒ identical B-spline interpolation ⇒ identical output. The rest of the
pipeline is untouched and RNG-driven from the same key.

**The concrete test the plan must specify — two parts:**

- **T1 (cheap, no fixture, same-commit):** in `spike/test/`, assert
  `warp(X, A(ε=0,dy,dx), axes(X); method=BSpline(Linear()), fillvalue=BG_FLOOR) ==
   warp(X, Translation(dy,dx), axes(X); method=BSpline(Linear()), fillvalue=BG_FLOOR)`
  for `X` a fixed `Random.seed!`-generated 256² matrix and `(dy,dx) ∈ {(0,0), (0.37,−0.82),
  (2.5,1.25), (−3,3)}`. Uses `==`, not `isapprox`. Requires no golden bytes.
- **T2 (pipeline-level, needs a PRE-EDIT fixture):** persist, **before touching stage 6**, a JLD2
  golden of `simulate_pair(Philox4x(UInt64,(P11_FIXTURE_SEED, k)), θ_fixed; imsize=(256,256))` for
  `k = 1:4` and an explicitly constructed 7-field `θ_fixed` (**not** `sample_prior`, whose RNG
  consumption changes the instant `ε` is added — see A4). After the edit, re-run with
  `θ_fixed ∪ (chromatic_eps = 0.0,)` and assert `==` on both channel matrices.

**Tolerance fallback.** Exact equality is the expected and measured outcome; a failure is a genuine
signal (e.g. a StaticArrays/Interpolations version bump changing the reduction order). If the
planner wants a documented escape hatch, the tightest defensible one is
`maximum(abs.(a .- b)) ≤ 1e-13` (a few ULP on values of order 0.02–3.0), and any use of it must be
recorded in the report as a deviation.

### A4. Where θ arity changes ripple — complete enumeration

| # | Site | Line(s) | Today | Required change | Blast radius |
|---|------|---------|-------|-----------------|--------------|
| 1 | `spike/simulator/prior.jl` `SHIFT_PRIOR` | :56 | `Uniform(-1.0, 1.0)` | **D-02**: `Uniform(-3.0, 3.0)` — but see §D13, λ-conditioned sampling replaces the constant draw in the research trainer | spike |
| 2 | `spike/simulator/prior.jl` new `CHROMATIC_PRIOR` | after :57 | absent | `const CHROMATIC_PRIOR = Uniform(-0.02, 0.02)` (D-09) | spike |
| 3 | `spike/simulator/prior.jl` `sample_prior` NamedTuple | :73-81 | 7 fields | 8th field `chromatic_eps = rand(rng, CHROMATIC_PRIOR)` — **changes the RNG stream for every downstream draw** | spike |
| 4 | `spike/simulator/forward.jl` validation `all(isfinite, …)` | :119-120 | 7 fields | add `θ.chromatic_eps` (or read defensively, see A5) | spike |
| 5 | `spike/simulator/forward.jl` stage 6 | :168-174 | `Translation` only | composed `AffineMap` (A2) | spike |
| 6 | `spike/simulator/forward.jl` header docstrings | :38, :44-45, :94 | "7-field NamedTuple", stage-6 description | update to 8-field + composed warp | doc only |
| 7 | `src/amortized/simulator.jl` `SHIFT_PRIOR` | :94 | `Uniform(-1.0, 1.0)` | **Leave at ±1** unless the planner wants the src mirror to change the *training distribution* — D-11 mirrors **`ε`**, not D-02's widening. See §A5 note. | **decision point** |
| 8 | `src/amortized/simulator.jl` `CHROMATIC_PRIOR` | after :95 | absent | add (D-11) | src, hashed |
| 9 | `src/amortized/simulator.jl` `sample_prior` | :109-120 | 7 fields | 8th field | src, hashed |
| 10 | `src/amortized/simulator.jl` `theta_prior_bounds()` | **:140-148** | `NTuple{7,…}` | append `(minimum(CHROMATIC_PRIOR), maximum(CHROMATIC_PRIOR))` ⇒ `NTuple{8,…}`; docstring :123, :125-126 | src |
| 11 | `src/amortized/simulator.jl` validation | :190-191 | 7 fields | add `chromatic_eps` (defensively) | src |
| 12 | `src/amortized/simulator.jl` stage 6 | :234-238 | `Translation` only | composed `AffineMap` | src, hashed |
| 13 | `src/amortized/simulator.jl` imports | :47 | `import CoordinateTransformations: Translation` | `…: Translation, LinearMap, recenter` | src, hashed |
| 14 | **`src/amortized/datagen.jl` θ buffer** | **:170** | `Matrix{Float64}(undef, 7, N)` | `8` — otherwise `theta[:, j] = s.theta` throws `DimensionMismatch` (`:178`) because `generate_sample` returns `collect(values(θ))` (`:143`) | **src — breaks datagen** |
| 15 | **`src/amortized/datagen.jl` `theta_dim`** | **:274** | `theta_dim = 7,` inside `generating_config` | `8` — this is a **hashed config field**, so it re-digests the cache dir independently of the source-byte change | src |
| 16 | **`src/amortized/datagen.jl` holdout buffer** | **:393** | `Matrix{Float64}(undef, 7, H)` | `8` | src |
| 17 | `src/amortized/datagen.jl` docstrings | :132, :157 | "7-vector", "theta::Matrix 7×N" | update | doc only |
| 18 | **`src/amortized/architecture.jl` `NPE_D`** | **:169** | `const NPE_D = 7` | See A5 — **do NOT change** unless the src training path is also being retargeted; `build_estimator`'s `D` default flows from here | src |
| 19 | `src/amortized/architecture.jl` header comment | :167-168 | "ρ_true stays θ ROW 1 and D stays 7" | update if :169 changes | doc only |
| 20 | `src/amortized/train_npe.jl` `build_estimator` call | **:136** | `build_estimator(d_in; …)` — **no `D` argument** | must become `build_estimator(d_in, length(theta_prior_bounds()); …)` if src-side 8-θ training is ever wanted; `train_npe` exposes **no `D` keyword** (:114-124) | src |
| 21 | `src/amortized/train_npe.jl` `arch` record | :147 | `D = NPE_D` | follows :136 | src |
| 22 | **`src/amortized/ood.jl` `_theta_tuple`** | **:124-132** | builds a 7-field θ from `v[1..7]` | must emit `chromatic_eps` (defensively: `length(v) ≥ 8 ? _finite_or(v[8], 0.0) : 0.0`) — this tuple is fed straight to `simulate_pair` at `:169` | **src — the shipped-compat hinge** |
| 23 | `src/amortized/ood.jl` `_theta_tuple` docstring | :121-123 | lists 7 ranges | update | doc only |
| 24 | **`test/runtests.jl`** literal 7s | **:495, :507, :513, :517, :535, :550, :556, :558, :562, :624, :700, :734, :796, :811, :813, :1185, :1399** | e.g. `@test length(b) == 7` (:535), `randn(7, n)` (:507) | every θ-shaped literal must become 8 (or `length(theta_prior_bounds())`) | **test suite goes red without this** |
| 25 | `test/gate/sbc.jl` `SBC_PARAM_LABELS` | :123-124 | 8 labels (7 θ + Δρ) | becomes 9 if the src-side gate is ever re-run; **not needed for Phase 11** (no gate re-run, D-01) | deferred |
| 26 | `spike/validation/sbc.jl` `SBC_PARAM_LABELS` | :245-246 | 8 labels | Phase-11 research reporting needs 9 (8 θ + Δρ) — spike-local | spike |

**Hard-coded `7` audit result:** `grep -n "theta_dim\|== 7\|randn(7" src/ test/` yields exactly the
sites above. There is no other literal-7 θ coupling in `src/`.

### A5. The mirroring asymmetry — does the 8th θ column break the shipped bundle? **NO.**

This is the brief's designated highest-risk unknown. Definitive answer with code evidence.

**The shipped read path, traced end to end:**

1. `colocalization_amortized(...)` → `estimator_for(8)` (`src/registry.jl:170-177`)
2. → `_lazy_load_from_artifact!(8)` (`registry.jl:126-161`) — resolves `Artifacts.toml`'s
   `grid_8` entry via `Artifacts.artifact_hash(name, toml)` (`:133`), verifies
   `Pkg.GitTools.tree_hash(dir)` against the pinned `git-tree-sha1` (`:101-108`)
3. → `_bundle_from_artifacts(8, npe_path, ratio_path, ood_path)` (`registry.jl:155`,
   `src/amortized/pipeline.jl:104`)
4. → `load_estimator(path)` (`src/amortized/persist.jl:92-103`):

```julia
    a   = d["arch"]
    est = build_estimator(a.d_in, a.D; dstar = a.dstar, depth = a.depth, width = a.width,
                          num_coupling_layers = a.num_coupling_layers,
                          flow_depth = a.flow_depth, flow_width = a.flow_width)
    Flux.loadmodel!(est, d["model_state"])
    return (estimator = est, θzt = d["theta_transform"], zt = d["zt"], …)
```

**`a.D` is read from the JLD2 file, not from `NPE_D`.** It was written at training time by
`train_npe.jl:147` (`arch = (d_in = …, D = NPE_D, …)`) and frozen into the artifact. So even if
`NPE_D` were changed to 8, the shipped bundle would still rebuild a **7-marginal** `NormalisingFlow`
and `Flux.loadmodel!` would succeed. `build_estimator`'s only guard is `dstar >= D`
(`architecture.jl:203`) — satisfied at `dstar = 64, D = 7`.

**`θzt` is likewise frozen data, not a recomputation.** It is a persisted
`BoundedThetaTransform(lo, hi, zt)` (`architecture.jl:100-104`) whose `lo`/`hi` are 7-element
`Vector{Float64}`s baked in at `fit_theta_transform` time (`train_npe.jl:66-77`).
`theta_prior_bounds()` appears **only** as the *default value* of `fit_theta_transform`'s `bounds`
keyword (`train_npe.jl:66`) — a training-time call. `StatsBase.reconstruct(t::BoundedThetaTransform, Y)`
(`architecture.jl:157-159`) reads `t.lo`/`t.hi`, never the module function. **Confirmed: changing
`theta_prior_bounds()` to return 8 entries cannot affect loading or reading the shipped bundle.**

**No version/compat guard is needed for the artifact.** But **one compat guard IS needed in code**,
and it is the single genuine shipped-path hazard:

> **`src/amortized/ood.jl:124-132` `_theta_tuple(v)` constructs a 7-field NamedTuple from a
> 7-row posterior-mean vector and hands it to `simulate_pair` (`ood.jl:169`).** If the mirrored
> `simulate_pair` unconditionally reads `θ.chromatic_eps`, this path throws.

**Mitigation (mandatory, both files):** read `ε` defensively in `simulate_pair`, so a legacy 7-field
θ is a valid input meaning "no chromatic aberration":

```julia
    ε = hasproperty(θ, :chromatic_eps) ? θ.chromatic_eps : 0.0
```

and extend `_theta_tuple` to emit the field when the vector is long enough:

```julia
    chromatic_eps = length(v) >= 8 ? clamp(_finite_or(v[8], 0.0), -1.0, 1.0) : 0.0,
```

With that guard, `simulate_pair(rng, θ_7field)` is **byte-identical to today** (A3), so the PP-channel
re-simulation and every other legacy call site keep working unchanged.

**How exposed is that path in practice?** Low, but non-zero. `registry.jl:157-159` records that
"The persisted `ood_nulls` carries only the `:density` channel", and `ood_verdict` only invokes
`pp_mismatch_score` when `haskey(ood_nulls, :model) && ood_nulls.model !== nothing`
(`ood.jl:371-374`), which the shipped bundle does not satisfy. So the shipped `is_ood(r)` path never
reaches `simulate_pair`. But `test/gate/misspec.jl`, `src/amortized/local_map.jl` and any future gate
run do — the guard is cheap insurance and should be non-negotiable.

**What must NOT cross the spike/src boundary (D-03):**

| Crosses into `src/` (D-11) | Stays spike-only (D-03) |
|---|---|
| `CHROMATIC_PRIOR` + the `chromatic_eps` θ field | The uncertainty level λ as a network input |
| The composed affine stage 6 | The 129-row input surface |
| The 8-entry `theta_prior_bounds()` | The λ prior / hierarchical sampler |
| The defensive `ε` read + `_theta_tuple` guard | `D = 8` flow construction / `NPE_D` |
| `theta_dim = 8` in `generating_config` | Any change to `build_estimator`'s `D` default |

**Recommendation on `SHIFT_PRIOR` in `src/` (row 7 of A4):** D-11 says "mirror the **`ε`** extension",
and D-02 scopes the ±3 px widening to "**the research net**". Mirroring the widening into `src/`
would change the *shipped simulator's* nuisance prior, which is a strictly larger provenance event
than D-12 describes and would make `theta_prior_bounds()[5:6]` disagree with the shipped bundle's
frozen `θzt.lo/hi` — a silent semantic drift on the reconstruct path if anything ever refits.
**Leave `src/amortized/simulator.jl:94` at `Uniform(-1.0, 1.0)`.** Flagged in Open Questions.

---

## B. Provenance (D-12)

### B6. What the datagen content hash covers, and what editing `simulator.jl` costs

**`src/amortized/datagen.jl:208-213` — the hashed source list:**

```julia
const DATAGEN_HASH_SRC_FILES = String[
    joinpath(@__DIR__, "..", "colocalization.jl"),   # patch() / correlation() (defines the summary)
    joinpath(@__DIR__, "..", "LoadImages.jl"),        # MultiChannelImage container
    joinpath(@__DIR__, "summary.jl"),                 # patch_summary(mci,G) / encode_d01
    joinpath(@__DIR__, "simulator.jl"),               # sample_prior / simulate_pair / build_mci (07-05)
]
```

**`datagen.jl:234-242` — how the digest is formed:**

```julia
function cache_hash(config::NamedTuple; src_files = DATAGEN_HASH_SRC_FILES)
    h = zero(UInt)
    for f in src_files                                     # FIXED order = deterministic
        h = hash(read(f), h)
        h = hash(0x00, h)                                  # inter-file separator
    end
    h = hash(_canonical(config), h)
    return string(h; base = 16, pad = 16)
end
```

**Coverage:** the *raw bytes* of all four files (so **comments and whitespace count** — the D-12
named-limit comment inside `simulator.jl` would itself flip the digest), plus a sorted, `repr`-based
canonical serialization of the `generating_config` NamedTuple (`_canonical`, `:221-224`). That config
carries `N`, `shard_size`, `master_seed`, `k`, `grid`, **`theta_dim`**, `imsize_set`,
`imsize_weights`, `summary_min_dim`, `schema_version` (`:268-279`).

**So D-11 flips the digest twice, independently:** once via the source bytes of `simulator.jl`, and
once via `theta_dim: 7 → 8` in the hashed config (A4 row 15).

**Cache directory naming:** `open_or_invalidate(cache_root, config)` (`:368-385`) resolves
`joinpath(cache_root, cache_hash(config))` — the digest **is** the directory name. Contents per dir:
`shard_0001.jld2 …` at `SHARD_SIZE = 10_000` samples/shard (`:54`, `:288`), plus `holdout.jld2` and
`meta.jld2`.

**Observable consequence of the edit:**

- A run with the new sources resolves to a **new, empty** directory and regenerates from scratch.
  It does **not** error — `open_or_invalidate` only errors when an *existing* dir's stored
  `meta.jld2` hash disagrees with the recomputation (`:379-381`), which cannot happen for a
  fresh dir.
- The **old** directories become orphaned: never read, never garbage-collected.

**Where they live and how big:** `cache_root` is a caller argument with no default in `datagen.jl`;
the pipeline supplies it. `artifacts/` is the repo-local root, and it is **gitignored** — a
`git ls-files artifacts` returns nothing while `artifacts/amended_v2/grid_8/` exists on disk with
`{npe,ratio,ood_nulls,gate_report}_8.jld2`. The shipped 8×8 pool was 50 000 pairs
(`scripts/train_grid8_amended_v2.jl:11`); at `summary_min` 128×Float64 + θ 7×Float64 + index +
imsize per sample that is on the order of **~55 MB/50k pool** plus JLD2 overhead — i.e. tens of MB,
not GB.

**Should they be deleted?** **No, and the plan should say so explicitly.** They are the only on-disk
record of the training pool the shipped net came from. Deleting them would destroy exactly the
provenance D-12 is trying to preserve by reference. Recommended plan wording: *"orphaned caches are
retained, not pruned; they are the byte-level fallback if the git-pinned reconstruction is ever
disputed."* If disk pressure ever forces a prune, `subhashes()` (`:250-252`) in the stale
`meta.jld2` records which source files produced them, so the pruning is auditable.

### B7. The `Artifacts.toml` pin is independent of the datagen digest — code path

**Evidence, in resolution order:**

1. `Artifacts.toml` (whole file, 8 lines) contains only:
   ```toml
   [grid_8]
   git-tree-sha1 = "90e6b63a8a234d067b407fefd7914f2ae4845448"
   lazy = true
       [[grid_8.download]]
       sha256 = "17162904c69b30e2a99a6dd934434553d08ebaadbfd55dea664b36efcf92b6d5"
       url = "https://github.com/ma-seefelder/ProteinCoLoc/releases/download/v2.0.0/grid_8.tar.gz"
   ```
   No datagen digest appears anywhere in it.
2. `src/registry.jl:133` — `expected = Artifacts.artifact_hash(name, toml)` reads that literal
   `git-tree-sha1` straight from the TOML. No call to `cache_hash`, `subhashes`, or
   `generating_config`.
3. `src/registry.jl:101-108` — `_verify_tree_sha1(dir, expected, grid)` computes
   `Base.SHA1(Pkg.GitTools.tree_hash(dir))` over the **bundle directory**
   (`artifacts/amended_v2/grid_8/`, `registry.jl:96-97`), i.e. over the four `.jld2` files. `src/`
   source bytes are not in that tree.
4. `src/registry.jl:137-147` — all three resolution branches (installed store / in-repo dev
   fallback / `ensure_artifact_installed` download) key on `expected` only.
5. `datagen.jl` is not imported by `registry.jl` and `cache_hash` appears nowhere in the
   `colocalization_amortized` call chain (`api.jl` → `registry.jl` → `pipeline.jl:104` →
   `persist.jl:92`).

**Conclusion (HIGH confidence):** after the `simulator.jl` edit, `colocalization_amortized(img,
control, channels; num_patches = 8)` and the artifact download/verification are **unaffected**. The
CONTEXT's D-12 correction is confirmed by code.

**The one caveat worth a sentence in the report:** `_train_grid_pipeline`'s `skip_if_done` branch
(`pipeline.jl:163-166`) checks only *artifact file integrity*, not the datagen digest — so a
post-edit `train_and_register(8)` against the existing `artifacts/amended_v2/` root would **load the
old bundle and not retrain**. That is the desired behaviour under D-01 (no reship), but it means the
digest change is invisible at that call site. Worth naming so nobody later mistakes silence for
"nothing changed".

### B8. The shape of an existing named limit, and how to get the pinned sha

**Verbatim, `docs/amortized.md:132-135` (limit #4 — the direct antecedent of D-05):**

> 4. **Nuisance-parameter marginals carry a ~0.08-SD residual drift** not certifiable as negligible at
>    M=2000. Under a pre-registered δ=0.10-SD equivalence (TOST) test, 4/6 nuisances pass; 2/6
>    (autofluorescence, label_efficiency) fail — not because the drift is large (point drift ~0.08 SD
>    < 0.10) but because the 90% CI edges just past 0.10. **No model lever removes it** (capacity
>    redistributes it, spike 012; μ-truncation fixes only the atoms at a real cost to |ρ|≥0.95
>    inference and ADVI comparability, spike 013). This is a **marginal-reproduction limit on
>    parameters the fixed patch-correlation summary cannot constrain** (5/8 columns are vacuous —
>    posterior ≈ prior), **not** a failure of the coloc inference. It is a named calibration limit.

**The house voice, distilled:** numbered item → **bold one-sentence claim** → the measured numbers
inline → what was ruled out and where the evidence lives (`spike NNN`, a filename) → a closing
sentence that bounds the scope of the damage ("…**not** a failure of the coloc inference"). All
limits sit under `## Named limits (non-negotiable honesty — carry these into the manuscript)`
(`docs/amortized.md:106`), introduced by `docs/amortized.md:108-109`. The list currently runs
**1–7**, not 1–4 (the CONTEXT's "#1–#4" undercounts); **D-12's limit is #8** and D-16's
interpretation note is a separate `##` section, not a numbered limit.

**Draft skeleton for limit #8 (planner to fill the sha):**

> 8. **The shipped `grid_8` bundle's training distribution is pinned by git reference, not by the
>    current tree.** `src/amortized/simulator.jl` gained a chromatic-aberration parameter `ε` in
>    Phase 11 (research lane, no reship). The shipped 8×8 net was trained by
>    `src/amortized/simulator.jl` **as of commit `<SHA>`**, before that parameter existed; the
>    current file is a strict superset whose `ε = 0` behaviour is byte-identical (regression-tested).
>    Because `simulator.jl` is a `DATAGEN_HASH_SRC_FILES` member (`src/amortized/datagen.jl:208-213`),
>    the edit re-digests the training-data cache directory and orphans the pool the shipped net was
>    trained from. **The artifact itself is unaffected** — it is pinned by `git-tree-sha1` in
>    `Artifacts.toml`, not by the datagen digest — so download, verification and
>    `colocalization_amortized(...)` are unchanged. The training distribution remains exactly
>    recoverable via `git show <SHA>:src/amortized/simulator.jl`. This is a provenance-by-reference
>    repair, **not** a change to the shipped estimator.

**Command to capture the sha — and the ordering constraint.** The sha must be `HEAD` **immediately
before** the ε commit, so it must be captured as a *task output* before any edit:

```bash
git rev-parse HEAD          # e.g. 3e92430ca1f8ca59dc08d6912abd68684bd2d575
git rev-parse --short HEAD  # 3e92430
```

**Verification that the captured sha is the right one** (run *after* the ε commit lands):

```bash
# The pinned sha must be the FIRST PARENT of the ε commit, and its simulator.jl must lack `chromatic_eps`.
git rev-parse "<EPS_COMMIT>^"                                  # must equal <SHA>
git show "<SHA>:src/amortized/simulator.jl" | grep -c chromatic_eps   # must be 0
```

**D-12 sequencing constraint (from the CONTEXT, restated for the planner):** the `simulator.jl` edit
and the `docs/amortized.md` limit must land **in the same commit**. That creates a chicken-and-egg
with the sha: the sha is `HEAD` *before* the commit, which is knowable at authoring time. So the
task order inside the single commit is: `git rev-parse HEAD` → write the limit text with that sha →
edit `simulator.jl` → `git add` both → commit. The plan must state this explicitly or the executor
will reach for `HEAD` after committing and pin the wrong sha.

---

## C. Pre-flight Probe (D-06)

> **What follows includes an INDICATIVE probe I ran on 2026-07-25.** It is **not** the
> pre-registered D-06 probe: single fixed θ, `R = 12` paired replicates, 256² only, no θ averaging,
> no imsize mixture, and it consumed an ad-hoc key (`0xBEEF`/`0xCAFE`/`0xFEED`) that must therefore
> be added to the Phase-11 forbidden-seed set. Its purpose is to de-risk the phase and to show the
> planner that the metric and the arithmetic work.

### C9. Probe design

**The 128-dim summary producer.** Two equivalent chains, one per environment:

| Environment | Chain | Sites |
|---|---|---|
| spike | `build_mci(simulate_pair(rng, θ; imsize)) → patch_summary(mci) → encode_d01(M)` | `spike/contract.jl:60-70`, `spike/contract.jl:86-91` (fixed 8×8), `spike/data/encode.jl:72-76` |
| src | `build_mci(...) → patch_summary(mci, G) → encode_d01(M)` | `src/amortized/simulator.jl:260-270`, `src/amortized/summary.jl` (`patch_summary`, `encode_d01:72-76`) |

Rows `1:64` are the imputed per-patch correlations, rows `65:128` the binary present-mask
(`encode.jl:60-66`). **The probe must operate on rows 1:64 only** — the mask rows are near-constant
and are precisely the rows `_summary_row_partition` excludes from standardization
(`src/amortized/summary.jl:89-98`); including them would dilute the statistic with a
near-zero-variance block.

**Holding θ fixed while sweeping.** Do **not** call `sample_prior`: adding `ε` changes its RNG
consumption (A4 row 3), so a `sample_prior`-based probe is not comparable across the edit. Construct
θ explicitly:

```julia
mk(ρ, dx, dy, ε) = (ρ_true = ρ, spillover = 0.1, autofluorescence = 0.05,
                    label_efficiency = 0.8, shift_dx = dx, shift_dy = dy,
                    noise = 0.5, chromatic_eps = ε)
```

**The PAIRED design is load-bearing, not a nicety.** Each replicate `r` reuses the *same*
`Philox4x(UInt64, (P11_PROBE_SEED, r))`, so stages 1–5 (the latent fields, the thinning draws, the
PSF) are bit-identical across the sweep and only the stage-6 geometry differs. Measured
justification:

```
independent-seed ||Δs||₂ (same θ, different key) = 2.653  (sd 0.132)
|shift| = 3.0 px, PAIRED                          = 1.281  (sd 0.162)
```

The between-seed Monte-Carlo floor is **2.1× larger than the entire 3 px effect**. An unpaired probe
would report noise. [VERIFIED: measured 2026-07-25]

**The sensitivity metric — recommendation and justification.**

> **Metric: `Δρ_eq(x) = the Δρ that displaces the summary as far as perturbation x does`,**
> obtained by inverting a calibration curve `‖Δs‖₂(Δρ)` measured in the same paired design, where
> `Δs = s(θ_perturbed) − s(θ_base)` restricted to the 64 continuous rows.

**Why this and not the alternatives:**

- **Why not raw `‖Δs‖₂`?** It has no units. D-06's own wording demands "relative to the movement
  induced by a known `Δρ`", so the Δρ-referenced form is what the decision literally asks for, and
  it is directly interpretable by a reviewer ("3 px of misalignment looks like a Δρ of 0.19").
- **Why not Mahalanobis under the training-summary covariance?** It is the natural OOD-flavoured
  choice and it is *available* (`fit_ood_nulls`/`maha_score`, `src/amortized/ood.jl:77-106`), but it
  requires a fitted `Σ` from a training pool — which does not exist until after datagen, violating
  D-06's "simulator-only, before training, minutes". It also answers a different question ("is this
  summary unusual?") than SC2's ("how much information about Δρ does misalignment destroy?").
  **Recommend reporting it as a secondary column if a pool happens to be available, never as the
  gate input.**
- **Why Euclidean on the 64 raw correlation rows and not the z-scored ones?** All 64 rows are the
  same physical quantity (a Pearson correlation in [−1,1]) on the same scale, so z-scoring adds a
  fitted nuisance for no gain. Stated as a pre-registration choice.
- **Robustness:** report the **median** over replicates alongside the mean, because a single
  degenerate patch (`missing` → 0 via `encode_d01`) can move `‖Δs‖₂` by O(1).

**Full pre-registered probe specification (what the plan should lock):**

| Knob | Value | Rationale |
|---|---|---|
| θ base points | ≥ 5 draws from `sample_prior` on `P11_PROBE_SEED`, then held fixed | avoids over-fitting the metric to one θ; the indicative run used 1 and is therefore weaker |
| Replicates R | 32 per (θ, rung) | `sd/√R` ≈ 0.03 on ‖Δs‖ at the measured sd ≈ 0.16 |
| `imsize` | the F5 mixture `((512,512),(1024,1024),(1376,1028),(2048,2048))` w = `(.40,.25,.25,.10)`, **plus** a 256² arm for comparability with the indicative run | F5: the probe must speak about the regime the net will train in (`gate_consts_8_v2.jl:152-153`) |
| shift rungs | `0, 0.25, 0.5, 1, 1.5, 2, 3` px along the diagonal (`dx = dy = |s|/√2`) | matches the D-07 ladder |
| shift extension rungs (D-14) | `4, 5, 6, 8` px | SC3 breakdown range |
| ε rungs | `0, 0.0025, 0.005, 0.01, 0.02` | matches D-09's magnitude table |
| ε extension rungs (D-14) | `0.03, 0.05` | SC3 breakdown range |
| Δρ calibration rungs | `0.02, 0.05, 0.10, 0.20` at fixed base ρ | the reference scale |
| statistic | mean and median `‖Δs₁:₆₄‖₂`, plus `Δρ_eq` | see above |

### C9b. INDICATIVE probe results (256², ρ = 0.5, R = 12, paired)

**Shift sweep** (`dx = dy = |s|/√2`):

| \|shift\| (px) | mean ‖Δs‖₂ | sd | **Δρ_eq** |
|---|---|---|---|
| 0.00 | 0.000 | 0.000 | 0.000 |
| 0.25 | 0.380 | 0.075 | 0.052 |
| 0.50 | 0.391 | 0.069 | 0.054 |
| 1.00 | 0.452 | 0.057 | 0.063 |
| 1.50 | 0.806 | 0.131 | 0.118 |
| 2.00 | 0.872 | 0.118 | 0.128 |
| 3.00 | 1.281 | 0.162 | **0.191** |
| 5.00 | 1.910 | 0.194 | 0.288 |
| 8.00 | 2.938 | 0.227 | 0.447 |

**Δρ calibration curve** (shift = 0, same keys):

| Δρ | mean ‖Δs‖₂ | sd |
|---|---|---|
| 0.02 | 0.201 | 0.020 |
| 0.05 | 0.365 | 0.029 |
| 0.10 | 0.686 | 0.047 |
| 0.20 | 1.337 | 0.063 |

Ordinary-least-squares over the last three points: `‖Δs‖₂ ≈ 6.48·Δρ + 0.041`, hence
`Δρ_eq = (‖Δs‖₂ − 0.041)/6.48`. (The curve is near-linear over the range used; the plan should
re-fit on the pre-registered probe and report R².)

**ε sweep** (shift = 0, using the D-10 prototype):

| ε | corner displacement (px, 256²) | mean ‖Δs‖₂ | **Δρ_eq** |
|---|---|---|---|
| 0.0025 | 0.45 | 0.136 | 0.015 |
| 0.005 | 0.90 | 0.180 | 0.021 |
| 0.010 | 1.79 | 0.324 | 0.044 |
| 0.020 | 3.54 | 0.657 | **0.095** |
| 0.030 | 5.25 | 1.026 | 0.152 |
| 0.050 | 8.59 | 1.785 | 0.269 |

**ε at the real-data anchor 1376×1028** (R = 4): ε = 0.005 → 0.318; ε = 0.02 → **2.206**
(≈ 3.4× the 256² response, as expected from the larger half-diagonal).

**Verdict from the indicative run:** the summary is **comfortably above resolution** on both axes.
SC2 has ≈ 3.7× headroom in `Δρ_eq` across the intended λ ladder (0.052 → 0.191), and `ε` at its
prior edge is worth ≈ 0.10 `Δρ_eq` at 256² and considerably more at realistic sizes. **The D-06
"fails flat / below resolution" branch looks unlikely** — but the pre-registered probe still governs,
and §C10 specifies what to do if it disagrees.

**One design trap the indicative run exposed.** `‖Δs‖₂` jumps from 0.000 at shift = 0 to 0.380 at
shift = 0.25 and then is nearly flat to 1.0 px. That is not misalignment sensitivity — it is the
**onset of interpolation smoothing**: at exactly integer offsets the backward map samples on-grid
and `BSpline(Linear())` is an identity lookup, whereas any sub-pixel offset engages the
interpolation kernel. This is precisely the smoothing-decorrelates effect D-10 exists to avoid
duplicating, showing up as a step at the origin. **Consequences for the plan:**
(i) the reported probe curve must annotate the 0 → 0.25 step as interpolation onset, not signal;
(ii) the SC2 ladder must **not** include a literal `λ = 0` rung — start at `λ = 0.25`; under
D-03's λ-conditioning each rung draws `shift ~ Uniform(−λ, λ)`, which is continuous in λ and
smooths the step out, but a λ = 0 rung would be a qualitatively different (no-interpolation)
regime.

### C10. From probe output to the SC2 threshold and the D-08 band

**(a) The SC2 monotonicity threshold.** The probe supplies the *effect scale* that justifies the
pre-registered floor; the *noise scale* comes from the pre-registered `N` per rung and is handled by
a permutation null (see §E17). Concretely:

1. Compute the probe's own rank correlation across the ladder,
   `S_probe = corspearman(λ_rungs, mean_Δρ_eq_per_rung)`. On the indicative data
   `Δρ_eq` is strictly increasing across all seven rungs ⇒ `S_probe = 1.0`.
2. Pre-register **`SC2_SPEARMAN_FLOOR = 0.5 · S_probe`**, i.e. `0.5` on the indicative numbers.
   The 0.5 is an **attenuation allowance**: the posterior-width response to λ is a *composition* of
   (summary displacement → λ) with (posterior width → summary), and the second map is a learned,
   lossy stage. Halving is a round, conventional allowance fixed before the run — the same move
   `07-NUISANCE-SBC-SPEC-DRAFT.md:91-94` makes for δ = 0.10 SD.
3. Pre-register the **abort branch**: if the probe returns `S_probe < 0.9`, or if
   `Δρ_eq(λ_max) − Δρ_eq(λ_min) < 0.02` (i.e. the whole ladder is worth less than the smallest
   calibrated Δρ step), declare **"SC2 below resolution"** and take the §C10(c) branch.

**(b) The D-08 TOST tolerance band — the arithmetic.** The band is on *coverage*, but the probe
measures *summary displacement*, so the link must be made explicitly. Two steps:

*Step 1 — convert a coverage band into a width-error band.* For a symmetric posterior, if the
reported 90% interval is `±1.644854·c·σ` while the true scale is `σ`, empirical coverage is
`2Φ(1.644854·c) − 1`. Solving both directions:

| reported-width factor `c` | coverage | deviation from 0.90 |
|---|---|---|
| 0.92 | 0.8698 | **−3.0 pp** |
| 1.00 | 0.9000 | 0 |
| 1.10 | 0.9260 | **+2.6 pp** |

So **a ±3 pp coverage band ≈ a −8 % / +10 % error in the reported interval width.**

*Step 2 — check that band against the effect the phase is demonstrating.* The probe says the SC2
ladder should widen the posterior by a factor on the order of
`Δρ_eq(λ_max)/Δρ_eq(λ_min) = 0.191/0.052 ≈ 3.7×` (≈ +270 %). The equivalence band is therefore
**≈ 1/30 of the effect being claimed** — small enough that "coverage is provably close to nominal"
cannot be confused with "coverage moved with the ladder". That is an outcome-independent
justification anchored in the probe, exactly what D-08 requires.

*Step 3 — pre-register.* `SC2_COVERAGE_NOMINAL = 0.90`, `SC2_TOST_DELTA = 0.03`,
`SC2_TOST_ALPHA = 0.05` (per side), and record the Step-1/Step-2 derivation verbatim in the consts
file the way `gate_consts_8_v2.jl:189-198` records the Fisher-z derivation.

**(c) What "SC2 fails flat / below resolution" looks like numerically, and the required branch.**

*Numerically:* `Δρ_eq` differences between adjacent rungs at or below the paired replicate noise —
concretely `max_r Δρ_eq(λ_r) − min_r Δρ_eq(λ_r) < 0.02`, or `S_probe < 0.9`, or per-rung
`sd(‖Δs‖₂)/√R` exceeding the whole rung-to-rung increment. Note the indicative data already shows a
*local* flat stretch (0.25 → 0.5 → 1.0 px: 0.052, 0.054, 0.063) — so the criterion must be on the
**ladder span**, not on every adjacent pair, or the design fails on the interpolation-onset plateau.

*Required plan branch (the CONTEXT calls this a real possible outcome, so a plan without it is
incomplete):*

1. **Do not train.** The probe is explicitly the cheap gate before the expensive step.
2. **Re-parameterise the ladder, not the criterion.** Raise `λ_max` toward the D-14 extension range
   (up to 8 px) and/or move the ladder to the F5 realistic image sizes, where the same pixel shift
   is a *smaller* fraction of a patch but the ε response is much larger (measured: ε = 0.02 gives
   0.657 at 256² vs 2.206 at 1376×1028). Re-run the probe once. This is the "one documented
   iteration" D-04 pre-authorises, and spending it here — before training — is the cheapest possible
   place to spend it.
3. **If the re-run also fails:** report SC2 as **NOT DEMONSTRABLE at the fixed 8×8 summary
   resolution**, with the probe curve as the evidence, and route the remaining budget to SC3 (which
   is a coverage claim and does not require a resolvable width gradient) plus the D-16 docs note.
   Explicitly link the finding to the deferred "summary extension for registration identifiability"
   idea, since a below-resolution probe is direct evidence *for* that deferred work. **Do not**
   widen the ladder a second time or relax `SC2_SPEARMAN_FLOOR` — that is the Phase-7
   amend-twice pattern D-04 exists to avoid.

### C11. Runtime

**Measured on this machine** (Julia 1.12.6, `spike/Project.toml`, single-threaded; `nproc` = 32):

| `imsize` | `simulate_pair` |
|---|---|
| 256² | **44.9 ms/call** |
| 512² | **164.3 ms/call** |

[VERIFIED: measured 2026-07-25, 10 and 3 calls respectively after a warm-up call]

Cross-check against the repo's own recorded figures: `gate_consts_8_v2.jl:140-148` records
`512×512 → 0.0274 s/pair` and `1376×1028 → 0.0974 s/pair` "measured s/pair (32 thr)". My
single-thread 512² of 164 ms against their 32-thread 27.4 ms implies ≈ 6× effective parallel
speedup (memory-bandwidth-bound, not 32×) — consistent, and the right factor to use for budgeting.

**Probe cost estimate (pre-registered spec from §C9):**

- Rungs: 7 shift + 4 shift-extension + 5 ε + 2 ε-extension + 4 Δρ + 1 base = 23
- Evaluations: 23 rungs × 5 θ base points × 32 replicates = **3 680 simulate+summary calls**
- At 256² single-threaded: 3 680 × ~50 ms ≈ **3 min**
- At the F5 mixture (E[cost] ≈ 0.0835 s/pair at 32 threads, `gate_consts_8_v2.jl:148`):
  3 680 × 0.0835 s ≈ **5 min** threaded
- Running **both** arms: **< 10 minutes wall clock.** D-06's "minutes, no training" is accurate.

Add `patch_summary` cost: it is `patch()` + `correlation()` on the two channels
(`spike/contract.jl:86-91`); on the indicative run the full simulate+summary loop at 256² came in
under ~60 ms/call end to end, so the summary is a small fraction. Budget 1.3× the simulate figures.

---

## D. Training the Research Net (D-02, D-03)

### D12. `spike/npe/train_npe.jl` end to end, and what must change

**Current control flow (`spike/npe/train_npe.jl`):**

| Line(s) | What happens |
|---|---|
| :59-60 | guarded includes: `architecture.jl` (→ `build_estimator`, `NPE_*` consts), `../data/loader.jl` (→ `load_fold`) |
| :64-66 | `const NPE_MASTER_SEED = 0xC0FFEE` |
| :80 | `fit_theta_transform(θtr) = fit(ZScoreTransform, θtr; dims = 2)` — **plain z-score, not the bounded/logit transform**; the F2 `BoundedThetaTransform` lives only in `src/` (`src/amortized/architecture.jl:100`, `src/amortized/train_npe.jl:66`) |
| :102-116 | `train_fold(dir, fold; …)`; **`use_gpu && throw(...)` hard CPU-only guard at :114** |
| :116 | `fold_data = load_fold(dir, fold; master_seed, variant)` — reads a pre-built JLD2 cache; `Ztr` is **already standardized by the loader's train-only `zt`** |
| :119-121 | leak-free θ z-score fit on `θtr` only, applied to both |
| :123-126 | `d_in = size(fold_data.Ztr, 1)`; `est = build_estimator(d_in; dstar, depth, width, num_coupling_layers, flow_depth, flow_width)` — **no `D` argument ⇒ `D = NPE_D = 7`** (`spike/npe/architecture.jl:60,87`) |
| :138-141 | `train(est, θtr_std, θva_std, fold_data.Ztr, fold_data.Zva; epochs, batchsize, use_gpu=false, optimiser = AdamW(lr, (0.9,0.999), wd), stopping_epochs, verbose)` |
| :143-148 | result tuple + optional `save_npe` |
| :160-180 | `save_npe` — atomic `.tmp` → reopen-assert → `mv`; persists the **whole `estimator` object** (unlike src's `Flux.state` form) |

**The three required deltas:**

**(1) Widened `SHIFT_PRIOR` (D-02).** Not a `train_npe.jl` change at all — it is a datagen-time
change (`spike/simulator/prior.jl:56`). But under D-03 the constant prior is *replaced* by the
λ-conditional draw (see D13), so `SHIFT_PRIOR`'s role changes from "the prior" to "the widest rung"
(`SHIFT_LAMBDA_MAX = 3.0`). Recommendation: keep `SHIFT_PRIOR = Uniform(-3.0, 3.0)` in `prior.jl`
so `theta_prior_bounds()`-equivalent support boxes stay honest and so a non-λ code path still works,
and have the Phase-11 sampler override the two shift fields.

**(2) The 8th θ column.** `build_estimator(d_in; …)` must become
`build_estimator(d_in, 8; …)` — the positional `D` already exists (`spike/npe/architecture.jl:87`),
so **no signature change is needed**, only an argument. Guard `dstar >= D` is satisfied
(64 ≥ 8, `architecture.jl:95`). `fit_theta_transform` (a plain `ZScoreTransform`) is arity-agnostic.
Downstream, `spike/npe/infer.jl` docstrings say "7×N" at :92, :102, :109 — comments only, the code
is generic. `spike/validation/sbc.jl:159-171` hardcodes `Matrix{Int}(undef, M, 8)` and `for p in 1:7`
— must become 9 and `1:8`, and `SBC_PARAM_LABELS` (:245-246) gains `"ε"`.

**(3) The λ conditioning input (D-03) — recommended mechanism.**

> **Append λ as a 129th input ROW, after the frozen 128-row standardization.**

```julia
# spike-local; NEVER crosses into src/ (D-03)
const LAMBDA_MIN = 0.25
const LAMBDA_MAX = 3.0
# Fixed, invertible, data-independent map to ~[0,1] — NOT a fitted transform, so there is no
# leak and no frozen statistic to persist beyond these two constants.
encode_lambda(λ) = (λ - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN)

# Z128 is the OUTPUT of the existing frozen path: standardize_summary(encode_d01(M), zt, :min)
augment_input(Z128, λ) = vcat(Z128, fill(Float32(encode_lambda(λ)), 1, size(Z128, 2)))
```

**Why this and not a `DeepSet` / a second network branch:**

- `build_estimator` is **already documented as input-width-agnostic** — `src/amortized/architecture.jl:193`
  and `spike/npe/architecture.jl:76-79` both say it reads `d_in` and needs no per-grid change. Setting
  `d_in = 129` requires **zero** architecture code change.
- The summary network is a plain `Chain(Dense(d_in, width, gelu), …)`
  (`spike/npe/architecture.jl:101-106`), so a 129th input feature is a first-class conditioning
  variable that the MLP can mix into every learned summary, which then conditions **every coupling
  layer** of the flow. From the installed NeuralEstimators docs
  (`Estimators/PosteriorEstimator.jl` docstring): *"`NormalisingFlow` uses `t` as a conditioning
  input at each coupling layer."* [VERIFIED: installed v0.2.1 source, `Project.toml:3`]
- The fixed-data `train(estimator, θ_train, θ_val, Z_train, Z_val; …)` method
  (`NeuralEstimators/gFxuZ/src/train.jl:101`) is generic in `T`, so a `129×n` matrix is accepted with
  no API gymnastics. **No v0.2.1 feature is required that does not already exist in this repo's
  working code.**
- A `DeepSet` is the wrong tool: it exists for permutation-invariant *replicates*, and λ is a
  single scalar covariate of one dataset, not a replicate.

**The one trap.** `standardize_summary` delegates to `_summary_row_partition`
(`spike/npe/infer.jl:65` → `Loader._row_partition`, `spike/data/loader.jl:125-135`), which
**errors** on a non-128 row count for `:min` (`loader.jl:126`: `nrows == 128 || error(...)`). The
src twin errors on an odd count (`src/amortized/summary.jl:91-92`, `iseven`). Both are *good* — they
make "append λ then standardize" fail loudly. **The order must therefore be: standardize the 128
rows first (frozen path, untouched), then `vcat` the λ row.** State this in the plan as an
invariant; it also guarantees `zt` stays exactly the shipped-shape 64-row transform.

### D13. Hierarchical sampling of λ — the crux of D-03

**The generative story the net must be trained against:**

```
for each training sample i = 1 … N:
    λᵢ            ~ Uniform(LAMBDA_MIN, LAMBDA_MAX)          # the UNCERTAINTY LEVEL
    shift_dxᵢ     ~ Uniform(-λᵢ, +λᵢ)                        # shift drawn CONDITIONAL on λᵢ
    shift_dyᵢ     ~ Uniform(-λᵢ, +λᵢ)
    εᵢ            ~ Uniform(-0.02, +0.02)                    # D-09: FIXED prior, NOT λ-scaled
    ρ_trueᵢ       = ghat(μ*ᵢ),  μ*ᵢ ~ MU_PRIOR               # unchanged (prior.jl:72,74)
    (other 4 nuisances) ~ their existing priors              # unchanged (prior.jl:75-77,80)
    θᵢ            = (ρ_trueᵢ, spillover, autofl, label_eff, shift_dxᵢ, shift_dyᵢ, noiseᵢ, εᵢ)
    imsizeᵢ       ~ the F5 categorical (sample_imsize)
    yᵢ            = simulate_pair(rngᵢ, θᵢ; imsize = imsizeᵢ)
    Zᵢ            = vcat( standardize_summary(encode_d01(patch_summary(build_mci(yᵢ))), zt, :min),
                          encode_lambda(λᵢ) )                # 129 rows
    train on (θᵢ, Zᵢ)        # 8 θ rows; λ is NOT a θ row
```

**Why this exact structure and not the alternatives — three properties that must all hold:**

1. **λ is a CONDITIONING INPUT, never a θ row.** The flow has `D = 8` marginals
   (`ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise, ε`). λ is
   not inferred, so it must not be in θ. If λ were a θ row, the net would be asked to *estimate*
   the user's own stated uncertainty from the image — which is both meaningless and would make the
   read-at-a-chosen-λ sweep impossible.

2. **λ is drawn FIRST, and the shift is drawn CONDITIONAL on it.** This is what makes the sweep
   valid. The NPE's objective is a forward-KL fit to `p(θ | Z)` under whatever joint it is trained
   on. Here the joint is `π(λ)·π(θ|λ)·p(y|θ)`, so the learned `q(θ | Z₁₂₈, λ)` targets
   `p(θ | y, λ)` — the posterior a user would hold *given they believe their registration
   uncertainty is λ*. That is exactly the quantity SC2 sweeps. If instead λ were drawn *after* the
   shift (or independently of it), λ would carry no information about the shift and the net would
   correctly learn to ignore the 129th row — **SC2 would fail flat for a reason that has nothing to
   do with resolution.** This is the single most important implementation detail in the phase.

3. **`Uniform(-λ, λ)`, not a Gaussian of scale λ.** The half-width parameterisation is what D-03
   names verbatim ("the shift prior's half-width"), it keeps `λ = 3` exactly equal to D-02's widened
   `SHIFT_PRIOR = Uniform(-3, 3)` (so the widest rung *is* the stated prior), and it keeps the
   marginal support bounded, which matters because `spike`'s `fit_theta_transform` is a plain
   z-score with no support box.

**The marginal-vs-conditional consequence the report must state.** Integrating λ out, the *marginal*
prior on `shift_dx` is **not** `Uniform(-3,3)`: it is
`p(s) = ∫ 1/(2λ)·1[|s|≤λ] dπ(λ)`, a peaked, heavier-at-zero density. Therefore:

- the `theta_prior_bounds()`-style support box `[-3, 3]` remains correct (support is unchanged),
- but any SBC/rank statistic on the `shift_*` columns must be computed **conditional on λ** (i.e.
  ranks pooled within a rung), never against a `Uniform(-3,3)` reference. Given D-05 declares those
  columns vacuous and non-evidential anyway, the practical rule is: **report `shift_*`/`ε` ranks
  per-rung, labelled non-evidence** (see §E19).

**ε is deliberately NOT λ-scaled.** D-09 states the prior as literally `Uniform(-0.02, 0.02)`. Tying
ε to λ would make its marginal a mixture, contradicting D-09's literal text. Consequence, stated
plainly: **the SC2 ladder is a registration ladder only**; the chromatic term contributes a constant
marginal spread at every rung, and its own dose-response is measured by SC3/D-14's fixed-ε injection
(which needs no conditioning input). D-03's wording — "the registration-uncertainty level (the shift
prior's half-width)", singular — supports this reading. Flagged in Open Questions in case the
planner wants a second conditioning input.

**Reading at inference (the SC2 sweep):** for a fixed image `y`, compute `Z₁₂₈` **once**, then
evaluate `sampleposterior(est, vcat(Z₁₂₈, encode_lambda(λ_r)))` for each rung `λ_r`. The image is
never re-simulated across rungs, so the sweep is *exactly* the controlled comparison D-03 promises,
with zero simulation noise between rungs.

**Seeding caveat inherited from the harness.** `posterior_for` → `sampleposterior` **threads no
`rng`** and draws the flow's base samples from the **global** RNG — documented at length in
`test/gate/harness.jl:45-72` and fixed there by `seed_gate_global!(G)` (`:92-96`). The Phase-11
evaluation loops must do the same (`Random.seed!` derived deterministically from `P11_DEV_SEED`)
or the reported ladder will not be reproducible. This was a real regression once already; the plan
must wire it into every arm.

### D14. Training budget

**Recorded in-repo timings:**

| Source | Figure |
|---|---|
| `07-05-SUMMARY.md:90` | 50 000 pairs at **256² only**: datagen ≈ 6 min; whole pipeline ≈ **21 min** |
| `07-05-SUMMARY.md:100,120` | 50 000 pairs at the **full default mixture (up to 2048²)**: killed at **~45 min still in datagen** |
| `07-05-SUMMARY.md:120` | steady-state NPE training **≈ 1 s/epoch at ~10k train, ≈ 3.5 s/epoch at ~42.5k train** (an early 14.7 s/epoch probe was compile-dominated) |
| `07-06-SUMMARY.md:71` | grid-16 run: train ≈ **2.7 min** compute + gate ≈ 6 min |
| `gate_consts_8_v2.jl:148` | F5 mixture `E[cost] = 0.08346 s/pair` ⇒ **50 000 pairs ≈ 4 173 s ≈ 1.16 h** (32 threads) |
| measured here | `simulate_pair` 256² = 44.9 ms, 512² = 164.3 ms single-thread; `nproc` = 32 |

**Recommended budget:**

| Item | Setting | Cost |
|---|---|---|
| Pool size | `n_pairs = 50 000` | matches every prior run (`train_grid8_amended_v2.jl:11`), so the comparison to the shipped net is capacity-controlled |
| `imsize` mixture | the F5 set/weights, read from `test/gate/gate_consts_8_v2.jl:152-153` (**never retyped** — the same discipline as `scripts/train_grid8_amended_v2.jl:34-35`) | **≈ 1.2 h datagen** at `-t auto` on 32 threads |
| NPE epochs | 300, `stopping_epochs = 40`, `batchsize = 128`, `lr = 2.5e-4`, `weight_decay = 1e-4` — the Phase-5-calibrated recipe (`spike/npe/train_npe.jl:102-113,132-137`) | **≈ 15–20 min** at 3.5 s/epoch × ~300 epochs, less with early stopping |
| Architecture | `dstar = 64, depth = 3, width = 256, coupling = 10, flow_depth = 2, flow_width = 128` (`spike/npe/architecture.jl:60-66`), `D = 8`, `d_in = 129` | unchanged capacity ⇒ comparable to the shipped net |
| **Total** | | **≈ 1.5 h wall clock, CPU, single run** |

**GPU.** CLAUDE.md makes CPU the reproducible baseline and CUDA an optional accelerator. The spike
trainer has a hard `use_gpu && throw` guard at `spike/npe/train_npe.jl:114`; the src trainer removed
it (`src/amortized/train_npe.jl:107-108`, defaults to `has_cuda_device()`). Since training is only
~20 min of the 1.5 h and **datagen dominates** (and datagen is CPU-bound regardless), **GPU buys
almost nothing here**. Recommend keeping the CPU-only guard and stating so in the report. If the
planner wants GPU anyway, note the AdamW-Float64 gotcha (`spike/npe/train_npe.jl:129-131`).

**Do NOT route through `_train_grid_pipeline`.** It cannot express this net:
`src/amortized/train_npe.jl:136` calls `build_estimator(d_in; …)` with **no `D`** and `train_npe`
exposes **no `D` keyword** (:114-124), so `D` is pinned to `NPE_D = 7`; and its `datagen` seam
returns a `7×N` θ that `generate_samples` allocates at `datagen.jl:170`. The `datagen`-injection
pattern that spike 013 used (`.planning/spikes/013-mu-prior-truncation/trunc_datagen.jl`) works for
a *prior* change but **not** for an arity change. **Write a spike-local trainer**, which is what
D-01 wants anyway. It may still reuse `spike/data/generate.jl` + `spike/data/cache.jl` +
`spike/data/loader.jl` with the θ-buffer width parameterised.

---

## E. Evaluation: SC2 Ladder + SC3 Breakdown (D-07, D-08, D-13, D-14)

### E15. What `spike/validation/harness.jl` gives you, and the minimal extension

**Reusable as-is:**

| Function | Site | Role |
|---|---|---|
| `load_frozen_model(path)` | `harness.jl:78-80` | one-time frozen load via `load_npe`; **rename/param for the Phase-11 artifact path** |
| `val_rng(master_seed)` | `harness.jl:93-94` | `Philox4x(UInt64, (seed ⊻ VAL_SALT, 0))` — the salted-disjoint-stream idiom to copy |
| `draw_simulate_infer(m, rng; imsize, N)` | `harness.jl:106-113` | the θ*~π → simulate → standardize → `posterior_for` → `reconstruct` chain |
| `draw_simulate_infer_paired(m, rng; imsize, N)` | `harness.jl:124-136` | the **paired Δρ** path — returns `(θs, θc, ρs, ρc)`; consumers form `Δρ_draws = ρs .- ρc` and `Δρ* = θs.ρ_true - θc.ρ_true` |
| `rho_draws(est, Z, θzt; N, use_gpu)` | `spike/npe/infer.jl:120-124` | N un-standardized ρ draws — the Δρ building block |
| `interval_width(draws; probs)` | `spike/npe/infer.jl:148-151` | per-parameter CI width via NeuralEstimators `interval` — **use this, do not hand-roll quantiles** |
| `sbc_shrinkage(post_sd, prior_draws)` | `test/gate/sbc.jl:243-250` | the F3 vacuity diagnostic (see E19) |
| `sbc_holm_adjusted(p)` | `test/gate/sbc.jl:352-362` | Holm, valid under arbitrary dependence |
| `sbc_ranks_and_spread` | `test/gate/sbc.jl:169-210` | rank table **plus** `post_sd` and `prior_draws` — the shape the shrinkage diagnostic needs |
| `seed_gate_global!(G)` | `test/gate/harness.jl:92-96` | the global-RNG pin `sampleposterior` requires |

**The minimal extension (D-07/D-08 say EXTEND, not fork) — four small changes:**

1. **`draw_simulate_infer_paired` gains a `λ` argument** and, at the two `rho_draws` call sites
   (`harness.jl:133-134`), passes `augment_input(Zs, λ)` / `augment_input(Zc, λ)`. Everything else
   is untouched. Give λ a default of `nothing` meaning "no augmentation", so the 128-row path is
   preserved bit-for-bit for any legacy caller.
2. **A `forced` keyword on the θ draw** so SC3 can inject a fixed `shift`/`ε` (D-15) while leaving
   the other five fields prior-drawn:
   `θ = merge(sample_prior(rng), forced)` where `forced = (shift_dx = …, shift_dy = …, chromatic_eps = …)`.
   `merge` on NamedTuples preserves field order for existing keys, so `collect(values(θ))` stays
   row-aligned with `theta_prior_bounds()`. **The plan must assert that** — a one-line unit test.
3. **A new `p11_ladder(m, rng; rungs, N_per_rung, …)`** that loops rungs × draws and returns, per
   rung, `(widths, covered, post_sd, Δρ_star)`. ~40 lines. It composes the two functions above; it
   introduces no inference machinery, exactly as `harness.jl:24-26` describes the original design.
4. **Reuse, don't re-derive, the reductions:** `sbc_coverage` (`test/gate/sbc.jl:285-290`) and
   `interval_width` already compute the two quantities SC2 needs.

**Important reuse subtlety.** `spike/validation/harness.jl:109` hardcodes `patch_summary(mci)` (the
fixed 8×8 spike form, `spike/contract.jl:91`), whereas `test/gate/harness.jl:197` uses the
grid-parametric `patch_summary(mci, G)`. Phase 11 is 8×8 only (CLAUDE.md fixed-summary constraint),
so either is fine; prefer the **spike** harness as the base since D-01 puts everything in `spike/`
and it avoids depending on `test/gate/gate_consts_*.jl` load order.

### E16. The TOST / equivalence test, concretely in Julia

**The statistic.** Per rung `r`, over `N` independent paired draws, the indicator
`cᵢ = 1[ Δρ*ᵢ ∈ CI₉₀(Δρ_drawsᵢ) ]` where `CI₉₀` is the central 90% interval from `interval(draws;
probs = [0.05, 0.95])`. Empirical coverage `p̂_r = mean(c)`.

**Distributional assumption.** The `cᵢ` are i.i.d. Bernoulli(`p_r`) across draws — justified because
each draw uses an independent prior θ pair, an independent forward simulation, and (given the global
seed pin) an independent posterior sample. So `N·p̂_r ~ Binomial(N, p_r)`. **Use a Wilson score
interval, not Wald** — at `p = 0.9` and `N` in the low hundreds Wald's coverage is poor near the
boundary, and Wilson is the standard remedy.

**The two one-sided tests.** With nominal `p₀ = 0.90` and band `δ = 0.03`:

- `H₀₁ : p_r ≤ p₀ − δ = 0.87` vs `H₁₁ : p_r > 0.87`
- `H₀₂ : p_r ≥ p₀ + δ = 0.93` vs `H₁₂ : p_r < 0.93`

Equivalence is declared iff **both** are rejected at `α = 0.05`. The TOST p-value is
`p_TOST = max(p₁, p₂)`. Operationally identical (and easier to report) is the **CI form**: the
`(1−2α) = 90%` two-sided Wilson interval for `p_r` must lie **entirely inside `[0.87, 0.93]`**.

```julia
# Wilson score interval — no new dependency (Φ⁻¹ constants pre-registered, as gate_consts_8_v2.jl:351 does)
const Z_TWO_SIDED_90 = 1.644853627      # Φ⁻¹(0.95); the 90% two-sided limb = TOST at α = 0.05/side
function wilson_ci(k::Integer, n::Integer; z::Real = Z_TWO_SIDED_90)
    p̂ = k / n
    d  = 1 + z^2 / n
    c  = (p̂ + z^2 / (2n)) / d
    h  = (z / d) * sqrt(p̂ * (1 - p̂) / n + z^2 / (4n^2))
    return (lo = c - h, hi = c + h)
end

# One-sided normal-approximation TOST p-values (report BOTH these and the CI; they agree by construction)
_Φ(x) = 0.5 * (1 + erf(x / sqrt(2)))            # or use Distributions.cdf(Normal(), x)
function tost_pvalue(k::Integer, n::Integer; p0 = 0.90, δ = 0.03)
    p̂  = k / n
    se = sqrt(max(p̂ * (1 - p̂), 1e-12) / n)
    p1 = 1 - _Φ((p̂ - (p0 - δ)) / se)            # H₀₁: p ≤ p0 − δ
    p2 =     _Φ((p̂ - (p0 + δ)) / se)            # H₀₂: p ≥ p0 + δ
    return max(p1, p2)
end
```

**How the ±3 pp band enters:** as `δ`, justified by the §C10(b) arithmetic (±3 pp ↔ −8 %/+10 %
reported-width error, ≈ 1/30 of the ~270 % effect the ladder demonstrates), pre-registered in the
consts file **before** the probe's coverage numbers exist.

**Holm across rungs — and the direction trap.**

> **The Holm usage here is the MIRROR IMAGE of the existing SBC usage. Get the inequality right.**

`test/gate/sbc.jl:374` defines `sbc_holm_pass(p; fwer) = all(>(fwer), sbc_holm(p))` — "pass iff Holm
rejects **nothing**", correct for a *uniformity* null where rejection is bad. For equivalence,
rejection is **good**: the phase must reject every rung's non-equivalence null. So:

```julia
# Reuse the frozen adjuster; INVERT the comparison.
p11_tost_pass(p_tost::AbstractVector; fwer::Real = 0.05) =
    all(<=(fwer), sbc_holm_adjusted(p_tost))       # note <=, not >
```

A plan that copies `sbc_holm_pass` verbatim would compute the exact opposite verdict. This deserves
its own unit test with a hand-computed vector.

**An honesty note the report should carry (does not contradict D-08, strengthens it).** "All rungs
equivalent" is an **intersection–union** claim: the union-intersection principle says that requiring
every component TOST to reject at level α already has family-wise size ≤ α, so **no multiplicity
correction is required** and Holm makes the test *strictly harder*. D-08's choice is therefore
conservative in the right direction (the same direction `gate_consts_8_v2.jl:213-218` argues for the
BF lower-confidence-bound rule). Applying it is defensible and should be described as deliberately
conservative rather than as necessary.

**N per rung.** Requiring the 90% Wilson interval at `p̂ = 0.90` to fit inside a ±0.03 band:

```
half-width ≈ z·√(p(1−p)/N) = 1.644854·√(0.09/N) = 0.49346/√N ≤ 0.03
⇒ √N ≥ 16.449  ⇒  N ≥ 270.6  ⇒  N_min = 271
```

That is the *knife-edge* case (`p̂` exactly nominal). For slack when `p̂` lands a little off nominal,
**pre-register `N_PER_RUNG = 500`** (half-width 0.0221, so equivalence stays decidable for
`p̂ ∈ [0.878, 0.922]`). Record `N_min = 271` as the derivation, mirroring how
`gate_consts_8_v2.jl:193-198` records `n_min = 69` beside `BF_GATE_N = 100`.

**Reconciling with the M = 2000 over-power lesson.** The Phase-7 finding (skill
`references/sbc-calibration.md`, spike 011) is that an **M = 2000 point-null** KS test resolves a
0.059-SD drift and therefore rejects on drift too small to matter. **TOST inverts the incentive:
larger N makes equivalence *easier* to establish, because the CI shrinks toward the point estimate.
There is no over-power pathology for an equivalence test — only a cost ceiling.** So the Phase-7
lesson does not cap `N` here; it is the *reason* D-08 chose TOST in the first place. The report
should say this in one sentence, because a reader who remembers "M = 2000 was too big" will
otherwise object to `N = 500 × R rungs`.

**Cost.** `N = 500` × 8 rungs × 2 sims per paired draw = 8 000 `simulate_pair` calls. At the F5
mixture (0.0835 s/pair, 32 threads) ≈ **11 min**; posterior draws add a forward pass each and are
milliseconds. The SC3 breakdown ladder (§E18) adds ~10 more rungs ≈ **14 min**. **Evaluation is
cheap; training is the cost.**

### E17. The SC2 monotonicity test

**Statistic (pre-registered):** `W_r` = the mean over the rung's `N` draws of the **90% HDI width of
the Δρ draw cloud**, computed with `interval_width(hcat(Δρ_draws); probs = [0.05, 0.95])`
(`spike/npe/infer.jl:148-151`). Report posterior SD as a secondary column; the HDI width is
preferred because it is the quantity a user reads and because it is what D-07 names first.

**Test (pick ONE — this is the recommendation):**

> **One-sided permutation test on the Spearman correlation between λ and the per-draw width,
> plus a pre-registered floor on the point estimate.**
>
> 1. Pool all `R × N` pairs `(λ_r, w_{r,i})`. Compute `S = corspearman(λ, w)`
>    (`StatsBase.corspearman`, already a dependency).
> 2. Permutation null: shuffle the rung labels `B = 10 000` times, recompute `S`, and take
>    `p_perm = (1 + #{S_b ≥ S}) / (B + 1)`.
> 3. **PASS iff `p_perm ≤ 0.05` AND `S ≥ SC2_SPEARMAN_FLOOR`** (the probe-derived floor from
>    §C10(a), `0.5 · S_probe`).
> 4. Report, but do not gate on: the rung means `W_r`, the maximum downward rung-to-rung step
>    `max_r (W_r − W_{r+1})₊`, and an isotonic (pool-adjacent-violators) fit for the figure.

**Why this over the alternatives:**

- **vs. pairwise "every step ≥ −τ":** requires a τ with no outcome-independent source. The
  indicative probe already shows a *local* plateau (0.25 → 0.5 → 1.0 px gives 0.052, 0.054, 0.063
  `Δρ_eq`), driven by the interpolation-onset artifact — a strict pairwise rule would fail on a
  known simulator artifact rather than on the science.
- **vs. isotonic-fit goodness-of-fit:** the null distribution of an isotonic residual statistic is
  awkward to pre-register and needs a bespoke calibration run. The permutation null is exact under
  "λ has no effect", needs no distributional assumption, and is computable from the ladder data
  itself.
- **vs. Spearman with an asymptotic p-value:** the `R × N` widths are **not** independent within a
  rung in the way an asymptotic Spearman assumes (they share a rung). The permutation over *rung
  labels* is the honest null for exactly that structure.

**How the probe parameterises it:** the probe supplies (i) `S_probe`, hence the floor
`SC2_SPEARMAN_FLOOR = 0.5 · S_probe`; (ii) the ladder itself (which λ rungs are far enough apart in
`Δρ_eq` to be worth spending draws on); and (iii) the §C10(c) abort criterion. All three are locked
into the consts file *after the probe and before training* — see §G23.

**Pre-registered rung ladder (SC2):** `λ ∈ {0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0}` — 7 rungs, all
inside the training prior `[LAMBDA_MIN, LAMBDA_MAX] = [0.25, 3.0]`, no `λ = 0` rung (§C9b trap).

### E18. The SC3 breakdown curve (D-14)

**Two structurally different regimes — the plan must not blur them:**

| | **In-prior (GATED, SC3 pass/fail)** | **Beyond-prior (REPORTED only)** |
|---|---|---|
| Mechanism | λ-rungs; `shift ~ Uniform(-λ, λ)` | **fixed-magnitude injection** (D-15), λ held at `LAMBDA_MAX` |
| Range | λ ∈ {0.25 … 3.0}; ε ~ its prior | shift ∈ {4, 5, 6, 8} px; ε ∈ {0.03, 0.05} |
| Claim type | conditional-on-λ coverage — a *valid* SBC-style claim, because λ is a training conditioning variable | frequentist conditional coverage at a fixed out-of-prior θ — **not** implied by SBC |
| Verdict | TOST + Holm (§E16) | curve reported with CIs; **no gate** |

The in-prior arm is *the same runs as SC2(b)*: D-07(b) and SC3-in-prior are literally the same
coverage statistic on the same ladder. **Do not run them twice.** The plan should state that SC3's
in-prior pass is satisfied by the SC2(b) TOST result, and that the incremental SC3 work is only the
beyond-prior extension.

**The exact, pre-registerable definition of "the breakdown point":**

> **`BREAKDOWN_SHIFT` := the smallest rung magnitude `x` in the pre-registered beyond-prior ladder
> at which the upper limb of the 90% two-sided Wilson interval for empirical 90%-CI coverage falls
> strictly below the nominal 0.90 — i.e. the first `x` with `wilson_ci(k_x, N).hi < 0.90`.
> Reported as "coverage holds to < `x` px; degrades at and beyond `x`". If no rung on the ladder
> satisfies the condition, report `BREAKDOWN_SHIFT > 8 px` (right-censored). `BREAKDOWN_EPS` is
> defined identically on the ε ladder.**

Why the *upper* limb and not the point estimate: it makes "coverage has genuinely dropped" an
evidentiary claim rather than a coin flip on Monte-Carlo noise, matching the burden-of-proof logic
`gate_consts_8_v2.jl:211-218` already uses for the BF arm. Why *first* and not *last*: monotone
degradation is expected but not guaranteed; "first" is unambiguous under non-monotonicity, and the
full curve is reported anyway.

**Beyond-prior ladder (pre-registered):**
`shift ∈ {3.5, 4, 5, 6, 8}` px (diagonal, `dx = dy = x/√2`) with ε ~ prior;
`|ε| ∈ {0.025, 0.03, 0.04, 0.05}` with shift ~ `Uniform(-3, 3)`;
`N_PER_RUNG = 500` each; λ pinned at `LAMBDA_MAX = 3.0`.

**The λ-at-read-time choice must be pre-registered and stated.** Pinning λ = 3.0 answers "the user
already declares maximal registration uncertainty and *still* gets mis-registered data" — the most
favourable honest setting, so a breakdown there is a genuinely damning number. Recommend
additionally reporting a λ = 1.0 secondary curve ("the user under-declares"), which will break
earlier and is the more realistic user story. Flagged in Open Questions since D-14 does not settle it.

### E19. Vacuous-column reporting (D-05, F3)

**Where `shrinkage` is computed today:** `sbc_shrinkage(post_sd, prior_draws)` in
**`test/gate/sbc.jl:243-250`**:

```julia
function sbc_shrinkage(post_sd::AbstractMatrix, prior_draws::AbstractMatrix)
    n = size(post_sd, 2)
    psd_mean = [Statistics.mean(view(post_sd, :, p)) for p in 1:n]
    prior_sd = [Statistics.std(view(prior_draws, :, p)) for p in 1:n]
    shrink   = [prior_sd[p] > 0 ? psd_mean[p] / prior_sd[p] : NaN for p in 1:n]
    vacuous  = [isfinite(s) && s >= SBC_VACUOUS_SHRINKAGE for s in shrink]
    return (shrinkage = shrink, post_sd_mean = psd_mean, prior_sd = prior_sd, vacuous = vacuous)
end
```

with `const SBC_VACUOUS_SHRINKAGE = 0.95` at `:229` — deliberately defined **in the gate file, not
in the frozen consts** (`:226-228`), precisely so it is a reporting label and not a pass condition.
Its inputs come from `sbc_ranks_and_spread` (`:169-210`), which returns `post_sd` and `prior_draws`
as `M×8` matrices; the results are threaded into each `per_param` entry at `:455-456` and summarised
as `vacuous_params` at `:491`. `sbc_gate`'s docstring is explicit (`:431-435`): *"They are REPORTING
ONLY and feed no pass/fail decision."*

**There is no `src/` copy** — `grep shrinkage src/` returns only an unrelated comment at
`registry.jl:157`. So Phase 11 must either `include` `test/gate/sbc.jl` or (cleaner, given D-01's
spike containment) copy the ~8-line function into the Phase-11 spike validation module with an
attribution comment.

**How `ε` and the widened `shift_*` flow through it:**

- The θ column count goes 7 → 8, so the rank/spread matrices go `M×8` → `M×9` (8 θ + Δρ) and
  `SBC_PARAM_LABELS` gains `"ε"`. `sbc_shrinkage` is already generic in `n = size(post_sd, 2)`
  (`:244`) — **no change needed**.
- **Critical:** `prior_sd` must be the SD of the θ* values **actually drawn** — which, under D-13's
  hierarchical sampler, means the shift columns' `prior_sd` **within a rung** is
  `sd(Uniform(-λ_r, λ_r)) = λ_r/√3`, not the pooled marginal. Computing shrinkage against the pooled
  marginal would inflate `prior_sd` and make the shift columns *look* informative. The docstring's
  own defence — *"an empirical prior sample — no analytic form is assumed"* (`:236`) — holds only if
  `prior_draws` is restricted to the same rung. **Compute shrinkage per-rung.**
- `ε`'s prior is λ-independent, so its `prior_sd = 0.04/√12 = 0.01155` pooled or per-rung alike.

**The exact wording pattern for labelling a uniform rank on a vacuous column as non-evidence.**
The house voice already exists verbatim in three places; reuse it rather than inventing a phrasing:

- `gate-8x8-amended.md:46`: *"`shift_dx`'s clean KS p = 0.883 is meaningless — it is a vacuous
  column."*
- `gate-8x8-amended.md:125`: *"Any future SBC 'pass' on those columns must be read as vacuous."*
- `07-CALIBRATION-FINDINGS.md:86-88`: *"**This must NOT be read as calibration evidence anywhere** —
  not in the Go/No-Go memo, not in the manuscript, not in any per-grid comparison table. A vacuous
  pass and an informative pass are not the same result and must never be tabulated as if they
  were."*
- And the width/location caveat from the skill (`references/sbc-calibration.md:54-56`):
  *"`vacuous = shrinkage ≥ 0.95` means *uninformative about spread*, NOT *posterior = prior*. Both
  readings must be stated wherever the flag is used."*

**Recommended Phase-11 template (one sentence per vacuous column, in the report's SBC table
caption):**

> *`shift_dx`, `shift_dy` and `ε` carry per-rung shrinkage `s = <value>` (≥ 0.95), i.e. the posterior
> retains essentially all of the prior spread on those columns. Their rank statistics are therefore
> **vacuous by construction and are not calibration evidence** — they are tabulated for completeness
> only. Per spike 009, a vacuous column can still legitimately reject on a location shift, so a
> rejection here is also uninformative about the coloc claim. The Phase-11 claim rests on `ρ_true`
> and `Δρ`, whose shrinkage is `<value>`/`<value>`.*

---

## F. Real-Image Check (D-17)

### F20. The read-only path — and a scope correction the planner must absorb

**D-17 names `corpus/` as the target. `corpus/` cannot serve as the target.** Three independent
blockers, all verified:

1. **The bytes are not there.** `corpus/data/` is empty (`ls -la corpus/data` → only `.` and `..`),
   and `.gitignore` excludes `corpus/data/*` (`corpus/config.jl:109-112`: *"This tree is gitignored…
   large/licensed bytes never enter git"*).
2. **The physical anchors have no integrity hashes and are unfetchable as-is.** Both
   `physical-primary` rows carry `sha256 = PENDING-FETCH` (`corpus/anchor_rows.jl:85`,
   confirmed by `grep -c PENDING-FETCH corpus/manifest.csv` → 2), and `fetch_verified(url, dest,
   expected_sha256)` (`corpus/fetch.jl:80`) requires a real expected hash — a mismatch is a **hard
   error** (`corpus/config.jl:97-100`). Bootstrapping them (`bootstrap_anchor_hashes`,
   `anchor_rows.jl:247`) means a network fetch of a Zenodo zip and an EBI BioStudies tree.
3. **Decisive: both physical anchors are `split = sealed_holdout`.** `corpus/config.jl:55-58`:
   *"`sealed_holdout` is the blind physical-anchor split opened ONLY at Phase-16 evaluation (D-09);
   the sealed split is unreachable from every default manifest accessor by construction."* Opening
   it requires `open_sealed_holdout(df; reason)` (`corpus/manifest.jl:156`) and would **burn the
   Phase-16 blind evaluation** — a far larger credibility cost than anything Phase 11 buys.
   The 30 non-sealed rows are all `simulated-secondary` CBS benchmark images (and also have empty
   `sha256` / `bytes = 0`), i.e. **not real microscopy**, which defeats the purpose of D-17.

**What the repo actually has, and what D-17 should target.** `test/test_images/` holds **real,
committed, two-condition, three-channel microscopy TIFFs** used by the v1.0 manuscript test suite:

| Path | Measured | Used by |
|---|---|---|
| `test/test_images/positive/positive_c{1,2,3}.tif` | **1028 × 1376**, `RGB{N0f16}`, 8 488 000 B each | `test/runtests.jl:86`, `:166` |
| `test/test_images/negative/negative_c{1,2,3}.tif` | 1028 × 1376, `RGB{N0f16}` | (parallel negative condition) |

[VERIFIED: `Images.load` executed 2026-07-25; sizes and eltypes as shown]

**1028 × 1376 is exactly the D-08 real-data anchor size** carried in the F5 training mixture
(`gate_consts_8_v2.jl:152`), so these images sit squarely inside the regime the research net will be
trained on. They are tracked in git (`git ls-files` lists all six), require **no network**, are
**not** in any sealed split, and are literally "the existing real images already used by the
manuscript pipelines" — which is D-17's own phrasing. **Recommend retargeting D-17 to
`test/test_images/`.** This is a *target substitution within D-17's stated intent*, not a scope
change; it is flagged in Open Questions for the planner to confirm.

**The read chain (read-only, no `src/` edits):**

```julia
# src/LoadImages.jl:432-446 — the frozen convenience constructor
img = MultiChannelImage("positive",
        ["test/test_images/positive/positive_c1.tif",
         "test/test_images/positive/positive_c2.tif"],
        ["ch1", "ch2"])
# load_tiff (src/LoadImages.jl:214-217) = Float64.(Images.Gray.(Images.load(path)))
#   ⇒ img.data :: Vector{Matrix{Float64}}, each 1028×1376, values in [0,1]
#   ⇒ pixel_size = size(data[1]) (:444), otsu_threshold = Images.otsu_threshold.(data) (:445)
```

Then the *unchanged* summary path: `patch_summary(img)` → `encode_d01` → `standardize_summary(·,
m.zt, :min)` → `augment_input(·, λ)` → `sampleposterior`. `corpus/load.jl:75-79` (`summary_grid`)
is the analogous corpus-side helper and documents the same read-only-include discipline
(`corpus/load.jl:29-33`) — worth citing as the pattern even though the corpus is not the target.

### F20b. Injecting a KNOWN shift on a real image without a second interpolation pass

**The honest position: on real data a second interpolation pass is unavoidable, and the report must
say so.** D-15's "simulator-injected, not post-hoc" is achievable for *simulated* data because the
simulator owns stage 6. A real image arrives already-formed; any shift we impose is by definition
post-hoc and costs one resampling pass.

**Recommended handling — make the artifact measurable rather than pretend it away:**

1. **Use the same single composed `AffineMap`** on the real channel 2 that stage 6 uses, so the
   resampling kernel, `fillvalue` and axis convention are identical to the simulator:
   ```julia
   s = 1.0 / (1.0 + ε);  c = map(ax -> (first(ax)+last(ax))/2, axes(ch2))
   A = Translation(dy, dx) ∘ recenter(LinearMap([s 0.0; 0.0 s]), c)
   ch2w = Matrix{Float64}(collect(warp(ch2, A, axes(ch2);
                          method = BSpline(Linear()), fillvalue = BG_FLOOR)))
   ```
   **Exactly one pass**, not two — the composition still buys what D-10 bought.
2. **Include an ε = 0, shift = 0 "null warp" arm.** By §A3 this is a bit-exact identity, so it
   isolates *zero* smoothing and gives the reference posterior width. Any width increase at the
   smallest nonzero rung that also appears here would be a bug, not physics.
3. **Include a `shift = (1.0, 0.0)` integer-offset arm.** At exactly integer offsets the backward
   map samples on-grid and `BSpline(Linear())` is an identity lookup — so this arm has
   **misalignment without interpolation smoothing**. Comparing it against a `0.5 px` arm separates
   "the posterior widened because of misalignment" from "the posterior widened because resampling
   decorrelates". This is the cleanest available control and it costs one extra arm.
4. **State the residual limitation verbatim in the report:** *"On real images the injected
   misalignment necessarily involves one resampling pass that the native acquisition did not have.
   The integer-offset control arm bounds that contribution; the simulated sweep, where stage 6 owns
   the only pass, remains the primary evidence."*

**BG_FLOOR on real data.** `fillvalue = BG_FLOOR = 0.02` is a simulator constant chosen so
`_exclude_zero` does not thin patches (`forward.jl:63-66`). On a real image the natural background
differs. Recommend `fillvalue = median(ch2)` **or** cropping the analysis to the interior so no
fill-region pixel enters a patch; whichever is chosen must be pre-registered, because the fill
region sits at the image edge and the 8×8 grid's border patches would otherwise see an artificial
constant.

### F20c. Proving the manuscript pipelines are unmodified

**The repo's own definition of the claim** (`.planning/ROADMAP.md:126`):

> *"The main repo and both manuscript pipelines are demonstrably untouched — `git status` on `src/`,
> `bayes.jl`, `colocalization.jl` is clean (decoupling proof)"*

**Phase 11 cannot make that claim literally**, because D-11 deliberately edits
`src/amortized/simulator.jl`. The honest, checkable claim is narrower and must be worded precisely:

> **The v1.0 manuscript pipeline files and the real image bytes are byte-unchanged; only
> `src/amortized/` (the v2.0 subtree) is edited, and within it only `simulator.jl`, plus a
> defensive read in `ood.jl` and the θ-arity constants in `datagen.jl`.**

**Exact verification commands (run at report time; every one is deterministic and reviewable):**

```bash
# 1. The v1.0 manuscript pipeline sources are untouched since the Phase-11 base commit.
git diff --stat <PHASE11_BASE_SHA>..HEAD -- src/bayes.jl src/colocalization.jl \
                                             src/LoadImages.jl src/plot.jl
#    EXPECT: empty output.

# 2. Nothing outside src/amortized/ changed inside src/.
git diff --name-only <PHASE11_BASE_SHA>..HEAD -- src/ | grep -v '^src/amortized/'
#    EXPECT: empty output.

# 3. The real image bytes are untouched (belt and braces — they are read-only inputs).
git diff --stat <PHASE11_BASE_SHA>..HEAD -- test/test_images/
#    EXPECT: empty output.
git ls-files -s test/test_images/ | sha256sum      # record the digest in the report

# 4. The shipped artifact pin is unchanged.
git diff --stat <PHASE11_BASE_SHA>..HEAD -- Artifacts.toml
#    EXPECT: empty output.

# 5. The v1 manuscript tests still pass against the edited tree.
julia --project -e 'using Pkg; Pkg.test()'
```

Command 5 is the substantive one: `test/runtests.jl:86` and `:166` exercise the real-image
`MultiChannelImage` path, so a green suite is direct evidence the v1 read path is intact. (Recall
from §A4 row 24 that the suite **will** go red on θ-arity literals until those are updated — that is
expected work, not evidence of pipeline breakage, and the report must distinguish the two.)

### F21. Which images, how many, and the expected signature

**Recommendation: both conditions, channels `[1, 2]`, 2 images total.**

| Image | Role | Why |
|---|---|---|
| `test/test_images/positive/` ch1 vs ch2 | high-coloc arm | the width response should be largest where there is real correlation to destroy |
| `test/test_images/negative/` ch1 vs ch2 | low-coloc arm | near-zero `ρ`; misalignment has little correlation left to destroy, so the width response should be **markedly weaker** |

Two images is the right number: D-17 is explicitly **qualitative** ("this cannot be a coverage
test"), and the scientific content is the *shape of the width-vs-shift curve*, not a population
estimate. A third arm — the same positive image with channels `[1, 3]` — is cheap and gives a
second, independent high-coloc curve if the planner wants redundancy.

**Design: `Δρ` needs a sample and a control.** The amortized read surface computes
`Δρ = ρ(img) − ρ(control)` from two independent single-stack passes (`docs/amortized.md:24-25`;
`spike/npe/infer.jl:134-138`). For the width sweep use **positive as sample, negative as control**,
apply the injected warp to the **sample's channel 2 only**, and read the Δρ cloud width. That
mirrors the simulated ladder exactly (which also warps ch2 of one stack).

**Pre-registered rungs (mirror the simulated ladder so the curves are comparable):**
`shift ∈ {0 (null-warp), 0.5, 1.0 (integer control), 1.5, 2.0, 3.0}` px along the diagonal, plus
`ε ∈ {0, 0.005, 0.01, 0.02}` at shift = 0. Read at `λ = LAMBDA_MAX = 3.0` and again at
`λ = 1.0`. Total 2 images × 10 warp arms × 2 λ = 40 posterior passes — **seconds**, since the
expensive part (`patch_summary` on a 1028×1376 image) runs once per warp arm.

**Expected width-response signature (state this as a pre-registered prediction, so the check can
fail):**

1. **Monotone non-decreasing** Δρ-cloud width with injected shift magnitude on the positive arm.
2. **Steeper on the positive arm than the negative arm** — misalignment destroys correlation, and
   there is more of it to destroy at high `ρ`.
3. **Comparable in `Δρ_eq` terms to the simulated ladder.** From the probe, at the 1376×1028 anchor
   `ε = 0.02` moves the summary by `‖Δs‖ ≈ 2.21`; a matching magnitude on real data is the strong
   transfer evidence D-17 is after. A real-image response that is *much smaller* would say the
   simulator overstates how badly misalignment degrades patch correlation on real texture — a
   finding worth reporting either way.
4. **The null-warp arm (ε = 0, shift = 0) reproduces the un-warped posterior exactly** (§A3), and
   the integer-offset arm shows misalignment-without-smoothing.

**Non-signature (what would falsify transfer):** width flat in shift on both arms, or the negative
arm widening as much as the positive arm (which would indicate the width response is driven by
resampling smoothing rather than by decorrelation).

---

## G. Pre-registration File (D-04)

### G22. The `spike/validation/consts.jl` pattern, and the proposed Phase-11 seed

**Structure to mirror (`spike/validation/consts.jl`, 80 lines):**

| Element | Line(s) | What to copy |
|---|---|---|
| Anti-snooping preamble | :21-40 | The explicit contract: *"Every M/L/bin-count and every pass/fail threshold … is LOCKED in this ONE file and committed BEFORE any reported run. … Changing any value below after a reported run would be 'tune until calibrated' data-snooping"* |
| Seed-discipline paragraph | :30-34 | Names the disjointness requirement and *why* (never reuse the training stream) |
| Decoupling note | :36 | *"spike-local constants; touches no src/"* |
| Re-inclusion guard | :42, :80 | `if !isdefined(@__MODULE__, :SBC_M) … end` — one block keyed on a sentinel const, so a guarded re-include is a silent no-op |
| Grouped consts with inline rationale | :43-70 | One `# ---` section per arm, every constant carrying its justification as a trailing comment |
| Seed block with `≠` assertions in comments | :72-79 | `const VAL_MASTER_SEED = 0x5BC0FFEE   # ≠ 0xC0FFEE (NPE_MASTER_SEED): fresh disjoint stream` |

**The stronger pattern to also copy — `test/gate/gate_consts_8_v2.jl`:** it goes further and is the
better model for D-04 because it (a) *recomputes* the forbidden seeds rather than trusting comments
(`:275-285` recomputes every v1 `PROD_SEED[G]` purely to forbid them), (b) carries a
**`_forbidden_seeds()`** function and a **redraw loop** (`:298-322`), and (c) ends with
**executable `@assert` self-checks** on the frozen constants (`:434-453`). Phase 11 should adopt all
three.

**Complete inventory of `0x` seed/salt literals in the repo** (from
`grep -rhoE "0x[0-9a-fA-F_]{4,}" --include=*.jl .`, normalised and de-duplicated):

```
0x0000000000000001  (DEFAULT_MASTER_SEED)      0x0000000000c05eed  (CORPUS_MASTER_SEED)
0x00000000004a7107                             0x0000000000c0ffee  (NPE_MASTER_SEED)
0x0000000000a1ce01 / 02                        0x0000000000ca11b0 / b1 / b2 / b3
0x0000000000b7a770 / 71 / 72 / 73              0x0000000000dec0de
0x0000000000bf5014                             0x0000000009e3779b  (PROD_MASTER)
0x000000000de7c0d3 / 0de7c0de  (DEV_SEEDS)     0x000000005bc0ffee  (VAL_MASTER_SEED)
0x00000000c0de2026  (AMEND_MASTER)             0x00000000c0ffee5b
0x00000000de7c0de2  (DEV_SEEDS)                0x000000c0ffee5bc0
0x00a70115 / 0x00a70213                        0x000f1f7ed        (VAL_FIX_SEED, as 0xf1f7ed)
0x0134d8f3                                     0xa5a5a5a5a5a5a5a5 / 0xa5a5a5a5deadbeef
SALTS: 0x9e3779b97f4a7c15 (HOLDOUT) · 0xd1b54a32d192ed03 (FOLD / CBS_SPLIT)
       0xbf58476d1ce4e5b9 (VAL_SALT) · 0x94d049bb133111eb (PROD_SALT)
       0xc4ceb9fe1a85ec53 (AMEND_SALT)
```

Plus the **derived** (non-literal) seeds that must be recomputed and forbidden: `PROD_SEED[G]` for
`G ∈ (4,8,16,32)` (`gate_consts_8_v2.jl:275-285`) and `PROD_SEED_V2[G]` (`:315-331`).
Plus, **because this research session burned them**, the ad-hoc probe keys `0xBEEF`, `0xCAFE`,
`0xFEED`, `0xDEAD` and `Random.seed!(7)`.

**Proposed values (verified absent from the inventory above):**

```julia
const P11_DEV_SEED  = 0x0000_0000_0B11_DE71   # "0B11" = phase 11, "DE71" = DEV-1
const P11_SALT      = 0xA24B_AED4_663E_E121   # mx3 mixing constant; ≠ every salt in the repo
```

**The disjointness assertion code (executable, not a comment):**

```julia
if !isdefined(@__MODULE__, :P11_DEV_SEED)
    import Random123: Philox4x

    # ---- The FORBIDDEN set: every seed a Phase-11 run must not consume -----------------
    const NPE_MASTER_SEED     = 0x0000_0000_00C0_FFEE   # spike TRAINING stream
    const VAL_MASTER_SEED     = 0x0000_0000_5BC0_FFEE   # spike VALIDATION stream
    const VAL_FIX_SEED        = 0x0000_0000_000F_1F7E   # NOTE: literal is 0xF1F7ED — see below
    const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001   # productionization DATAGEN
    const CORPUS_MASTER_SEED  = 0x0000_0000_00C0_5EED   # corpus (D-17 target data)
    const F2_DEV_SEEDS        = (0x0000_0000_0DE7_C0DE, 0x0000_0000_DE7C_0DE2)
    const SPIKE_DEV_SEEDS     = (0x0000_0000_00CA_11B0, 0x0000_0000_00CA_11B1,
                                 0x0000_0000_00CA_11B2, 0x0000_0000_00CA_11B3,
                                 0x0000_0000_00DE_C0DE, 0x0000_0000_0DE7_C0D3)
    # Keys consumed by the 11-RESEARCH indicative probe (2026-07-25) — burned, not reusable.
    const P11_RESEARCH_BURNED = (0x0000_0000_0000_BEEF, 0x0000_0000_0000_CAFE,
                                 0x0000_0000_0000_FEED, 0x0000_0000_0000_DEAD)

    # v1 and v2 gate seeds are DERIVED, not literal — recompute them so they can be forbidden.
    const PROD_SALT   = 0x94D0_49BB_1331_11EB
    const PROD_MASTER = 0x0000_0000_09E3_779B
    const AMEND_SALT   = 0xC4CE_B9FE_1A85_EC53
    const AMEND_MASTER = 0x0000_0000_C0DE_2026
    _derive(master, salt, G) = rand(Philox4x(UInt64, (UInt64(master) ⊻ salt, UInt64(G))), UInt64)
    const PROD_SEED    = Dict(G => _derive(PROD_MASTER,  PROD_SALT,  G) for G in (4, 8, 16, 32))
    const PROD_SEED_V2 = Dict(G => _derive(AMEND_MASTER, AMEND_SALT, G) for G in (4, 8, 16, 32))

    _p11_forbidden() = (UInt64(0), NPE_MASTER_SEED, VAL_MASTER_SEED, 0x0000_0000_00F1_F7ED,
                        DEFAULT_MASTER_SEED, CORPUS_MASTER_SEED,
                        F2_DEV_SEEDS..., SPIKE_DEV_SEEDS..., P11_RESEARCH_BURNED...,
                        values(PROD_SEED)..., values(PROD_SEED_V2)...)

    # ---- The FRESH Phase-11 DEV stream ------------------------------------------------
    const P11_SALT     = 0xA24B_AED4_663E_E121   # mx3 constant; ≠ PROD/AMEND/VAL/HOLDOUT/FOLD salts
    const P11_DEV_SEED = 0x0000_0000_0B11_DE71   # phase 11, DEV-1

    p11_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P11_DEV_SEED) ⊻ P11_SALT, UInt64(counter)))

    # ---- EXECUTABLE self-checks (cheap; no inference, no simulation) -------------------
    @assert !(UInt64(P11_DEV_SEED) in _p11_forbidden())
    @assert P11_SALT ∉ (0x9E37_79B9_7F4A_7C15,   # HOLDOUT_SALT
                        0xD1B5_4A32_D192_ED03,   # FOLD_SALT / CBS_SPLIT_SALT
                        0xBF58_476D_1CE4_E5B9,   # VAL_SALT
                        PROD_SALT, AMEND_SALT)
    @assert length(unique(_p11_forbidden())) == length(_p11_forbidden()) - 0  # no accidental dupes
end
```

**Note on `VAL_FIX_SEED`.** The literal in `spike/validation/consts.jl:79` is `0xF1F7ED`
(6 hex digits = `0x0000_0000_00F1_F7ED`). The stub above shows a deliberately wrong-looking
placeholder to force the planner to read the source rather than trust this document — **use
`0x0000_0000_00F1_F7ED`**.

**Also lock in the same file** (this is D-04's substance, not just the seed): `LAMBDA_MIN`,
`LAMBDA_MAX`, the SC2 rung ladder, the SC3 beyond-prior ladder, `N_PER_RUNG`,
`SC2_COVERAGE_NOMINAL`, `SC2_TOST_DELTA`, `SC2_TOST_ALPHA`, `P11_HOLM_FWER`, `PERM_B`,
`SBC_L`, `P11_IMSIZE_SET`/`WEIGHTS` (read from `gate_consts_8_v2.jl`, never retyped), and — the
D-04 novelty — an explicit `const P11_ITERATION_ALLOWANCE = 1` with a comment stating that the
allowance is declared in advance and that a second iteration is not authorised by this file.

### G23. What must be locked BEFORE any run vs. what is legitimately probe-derived

This distinction drives the plan's task ordering, so it is drawn explicitly.

**Tier 1 — locked before ANYTHING runs (before even the probe):**

| Constant | Why it cannot wait |
|---|---|
| `P11_DEV_SEED`, `P11_SALT`, `p11_rng` | The probe itself consumes a stream; if the seed were chosen after seeing probe output it would be a snooped seed |
| `P11_PROBE_SEED` (may be `p11_rng(1)`) | ditto |
| The probe's own specification (§C9: θ base count, R, rungs, `imsize` arms, metric = paired `‖Δs₁:₆₄‖₂` and `Δρ_eq`) | Otherwise the metric could be selected to make the threshold convenient |
| **The §C10(c) abort criterion** (`S_probe < 0.9` or ladder span `< 0.02` `Δρ_eq`) | This is the one threshold that must be immune to the probe's own result, or the "below resolution" branch is unfalsifiable |
| `SC2_COVERAGE_NOMINAL = 0.90`, `SC2_TOST_ALPHA = 0.05`, `P11_HOLM_FWER = 0.05`, `PERM_B = 10_000` | Conventional, derivable with no data |
| `LAMBDA_MIN = 0.25`, `LAMBDA_MAX = 3.0` | Fixed by D-02/D-07's stated ladder |
| SC2 rung ladder, SC3 beyond-prior ladder | Fixed by D-07/D-14's stated ranges |
| `N_PER_RUNG = 500` and its `N_min = 271` derivation | Pure arithmetic on the band (§E16) |
| `P11_ITERATION_ALLOWANCE = 1` | D-04's whole point is declaring it in advance |
| `P11_IMSIZE_SET`/`WEIGHTS` (by reference to `gate_consts_8_v2.jl:152-153`) | F5 validity condition, no data needed |
| The `ε = 0` regression tolerance (`==`, fallback `1e-13`) | §A3 arithmetic |

**Tier 2 — legitimately probe-derived, locked at the END of the probe step, BEFORE training:**

| Constant | Source |
|---|---|
| `SC2_SPEARMAN_FLOOR = 0.5 · S_probe` | §C10(a) — the 0.5 attenuation factor is Tier 1; only `S_probe` is measured |
| `SC2_TOST_DELTA` (the ±3 pp band) | §C10(b) — the *arithmetic* is Tier 1; the probe supplies the effect scale that justifies calling 3 pp negligible |
| The `Δρ_eq` calibration coefficients (slope/intercept, R²) | §C9b — reporting only, but frozen so the ladder's `Δρ_eq` axis is fixed |
| Whether the ladder runs at 256², the F5 mixture, or both | §C10(c) branch 2 |

**Tier 3 — NOT pre-registerable, reported as measured:**

everything downstream of training — the ladder widths, coverages, TOST p-values, the breakdown
point, the real-image curves, the `ρ_true` attenuation measurement (D-02).

**The ordering rule this implies, stated for the planner:**

> The Tier-1 consts file is **committed** (its own commit, before the simulator edit). The probe
> runs against it. A **second, additive commit** appends the Tier-2 block to the *same* file, with
> the probe's numbers and a one-line provenance comment per constant. Training does not start until
> that second commit exists. No constant is ever *edited* — only appended — so the git history is
> itself the pre-registration audit trail.

---

## H. Risks and Landmines

### H24. Failure modes, ordered by probability × impact

**R1 — `test/runtests.jl` goes red on θ-arity literals. (probability HIGH ≈ 0.95 · impact MEDIUM)**
`grep` finds ~17 literal-`7` θ couplings, including `@test length(b) == 7` (`:535`),
`@test length(t.zt.mean) == 7` (`:550`), `randn(7, n)` (`:507, :624, :796, :1185`),
`randn(7, 500)` (`:556`), `size(draws,1) == 7` (`:495, :517, :734, :813`) and two synthetic datagen
stubs returning `randn(7, N)` (`:700, :1399`). **Early warning:** the very first `Pkg.test()` after
the D-11 commit. **Mitigation:** treat the test update as a first-class task *inside the D-11
commit*, and prefer `length(ProteinCoLoc.theta_prior_bounds())` over a new literal `8` so this never
recurs. **This is certain, cheap, and must be planned — not discovered.**

**R2 — `src/amortized/datagen.jl` throws `DimensionMismatch` on the first src-side datagen.
(probability HIGH ≈ 0.9 if src datagen is ever run · impact MEDIUM)**
`generate_samples` allocates `Matrix{Float64}(undef, 7, N)` (`:170`) and assigns
`theta[:, j] = s.theta` (`:178`) where `s.theta = collect(values(θ))` (`:143`) is now length 8.
Same at `:393` for the holdout. **Early warning:** any `_train_grid_pipeline` / `generate_cache`
call. **Mitigation:** update `:170`, `:274` (`theta_dim`), `:393` in the D-11 commit. Note the
failure is **loud**, which is the good case.

**R3 — the PP-channel `_theta_tuple` → `simulate_pair` break. (probability MEDIUM ≈ 0.4 ·
impact HIGH)** `ood.jl:124-132` builds a 7-field θ and `ood.jl:169` feeds it to `simulate_pair`.
The *shipped* path never reaches it (`registry.jl:157-159`: only `:density` is persisted, and
`ood_verdict` gates PP on `haskey(ood_nulls, :model)`, `ood.jl:371`), but `test/gate/misspec.jl`,
`src/amortized/local_map.jl` and any gate re-run do. **Early warning:** a `MethodError`/`type has no
field chromatic_eps` in any OOD or local-map test. **Mitigation (mandatory):** the defensive read
`ε = hasproperty(θ, :chromatic_eps) ? θ.chromatic_eps : 0.0` in **both** `simulate_pair`s, plus the
`length(v) >= 8` guard in `_theta_tuple`. With those, legacy 7-field θ is bit-identical (§A3).

**R4 — D-03's λ is drawn independently of the shift, so SC2 fails flat for the wrong reason.
(probability MEDIUM ≈ 0.3 · impact CRITICAL)** If the sampler draws `shift ~ Uniform(-3,3)` and λ
separately, λ carries no information and the net correctly learns to ignore the 129th input. The
result is a null SC2 that looks exactly like "below resolution" but is an implementation bug.
**Early warning (cheap, mandatory):** a **λ-ablation unit test** — feed one fixed `Z₁₂₈` at
`λ = 0.25` and `λ = 3.0` and assert the Δρ posterior SDs differ by more than Monte-Carlo noise.
If they do not, the conditioning is dead. Run this **immediately after training, before the ladder**.
Also worth an assertion inside the sampler that `abs(shift_dx) <= λ` for every generated sample.

**R5 — the probe returns "below resolution". (probability LOW ≈ 0.15 · impact HIGH)** Indicative
data says the summary moves plenty (3 px ≈ `Δρ_eq` 0.19; ε = 0.02 ≈ 0.095 at 256² and ≈ 0.34 at the
anchor), so the risk is lower than the CONTEXT feared — but the indicative run used one θ and one
image size. **Early warning:** the probe itself, ~10 min, before any training. **Mitigation:** the
pre-registered §C10(c) branch — do not train, re-parameterise the ladder once (D-04's single
allowance), and if it still fails, report SC2 as not demonstrable at the fixed 8×8 resolution and
route the finding to the deferred summary-extension idea.

**R6 — `ρ_true` attenuation from the wider nuisance prior. (probability MEDIUM ≈ 0.5 that it is
measurable · impact MEDIUM)** D-02 says explicitly this "must be MEASURED against the frozen net,
not assumed". **How to measure it, concretely — a paired, like-for-like comparison:**

1. Generate a **single** evaluation set of `K = 500` (θ*, image) pairs from the **shipped** prior
   (`SHIFT_PRIOR = Uniform(-1,1)`, `ε = 0`) on a fresh `p11_rng(...)` counter — one set, used by
   both nets, so the data is identical.
2. Read it through the **frozen shipped 8×8 net** (`estimator_for(8)`, 128-row input) → `ρ̂`, and
   through the **Phase-11 research net** at `λ = 1.0` (the level matching the shipped prior) →
   `ρ̂'`.
3. Report the three quantities that separate the possible stories:
   - `RMSE(ρ̂, ρ*)` vs `RMSE(ρ̂', ρ*)` — accuracy loss;
   - `mean(post_sd)` for each — **this is the attenuation**, i.e. how much wider the research net's
     `ρ_true` posterior is at matched conditions;
   - 90% coverage for each — whether the extra width is honest or slack.
4. Ratio-report: `attenuation = mean(post_sd') / mean(post_sd)`. A value near 1.0 at `λ = 1.0` is
   the good outcome (the conditioning did its job and the wider prior costs nothing when the user
   declares low uncertainty). A value well above 1.0 at `λ = 1.0` means the widening leaked into
   the low-λ regime and must be reported as a cost.

   **Caveat that must accompany the number:** the two nets differ in *two* ways (prior width and
   input surface), so the comparison is not a clean single-factor experiment — the same "the
   retrain is a confound" honesty `07-GATE-AMENDMENT` applies to itself (`gate_consts_8_v2.jl:28-45`).
   State it.

**R7 — Holm applied in the wrong direction on the TOST family. (probability MEDIUM ≈ 0.35 ·
impact HIGH)** `sbc_holm_pass` is `all(>(fwer), …)` (`test/gate/sbc.jl:374`) — correct for
uniformity, **exactly backwards for equivalence**. A copy-paste inverts the verdict silently.
**Early warning:** a unit test on a hand-computed p-vector (e.g. `[0.001, 0.002, 0.04]` must PASS
equivalence at FWER 0.05 under `all(<=, holm_adjusted(·))` and FAIL under `all(>, ·)`).
**Mitigation:** name the Phase-11 function `p11_tost_pass` (never reuse the name) and unit-test it.

**R8 — the `ε = 0` regression check silently degraded to `isapprox`. (probability LOW ≈ 0.2 ·
impact MEDIUM)** D-10's regression value comes from `==`. An executor hitting one failing case may
"fix" it with a tolerance. **Early warning:** grep the diff for `isapprox`/`atol` in the new test.
**Mitigation:** state in the plan that exact equality is *measured to hold* (§A3), that a failure is
a signal not a nuisance, and that any tolerance must be recorded as a deviation in the report.

**R9 — the pinned pre-ε sha captured after the commit. (probability MEDIUM ≈ 0.4 ·
impact MEDIUM)** D-12 requires the sha of `HEAD` *before* the ε commit, but the limit text ships
*in* that commit. **Early warning:** the post-hoc check `git rev-parse "<EPS_COMMIT>^"` disagreeing
with the pinned value. **Mitigation:** make "capture pre-ε sha" an explicit, ordered task producing
a recorded value, and add the two verification commands from §B8 to the plan's verification steps.

**R10 — NeuralEstimators conditioning-input API does not support D-03. (probability VERY LOW
≈ 0.05 · impact HIGH)** Assessed and effectively retired: `build_estimator` is documented
input-width-agnostic in both copies (`spike/npe/architecture.jl:76-79`,
`src/amortized/architecture.jl:193`); `train(estimator, θ_train, θ_val, Z_train, Z_val; …)` is
generic in the data type (`NeuralEstimators/gFxuZ/src/train.jl:101`); and the installed v0.2.1
`PosteriorEstimator` docstring states the summary vector conditions the flow "at each coupling
layer". **Early warning:** a dimension error on the very first `train` call — minutes, not hours.
**Mitigation:** a `d_in = 129` smoke train on 200 synthetic samples as a Wave-0 task, before the
1.2 h datagen.

**R11 — training non-convergence with the widened prior. (probability LOW ≈ 0.2 ·
impact MEDIUM)** Spike 012 found the current architecture is near its useful capacity ceiling
(val risk drifted 0.55 → 0.59; skill `references/sbc-calibration.md:104-107`), so the instinct to
"add capacity" is contraindicated. **Early warning:** validation risk not improving by epoch ~50,
or early stopping firing before ~80 epochs. **Mitigation:** keep the Phase-5-calibrated recipe
(lr 2.5e-4, batch 128, 300 epochs, patience 40 — `spike/npe/train_npe.jl:132-137`) unchanged; if it
stalls, **do not raise capacity** (012 PARTIAL: capacity *redistributes* drift, it does not reduce
it) — instead reduce `LAMBDA_MAX` and report the reduced ladder.

**R12 — Windows/CPU runtime blowout. (probability LOW ≈ 0.15 · impact MEDIUM)** The recorded
precedent is real: a 50k datagen on the full default mixture was **killed at ~45 min still in
datagen** (`07-05-SUMMARY.md:100,120`). **Early warning:** datagen not reaching shard 2 of 5 within
20 min. **Mitigation:** the sharded cache already supports **resume-by-skip**
(`datagen.jl:445-446`, `shard_done`), so a killed run resumes; launch with `julia -t auto` (32
threads available); and pre-register a fallback to a cheaper mixture with the F5 limit recorded via
`training_imsize_provenance` (`persist.jl:143-158`). Do **not** silently fall back to 256² — that is
the exact 07-05 artifact `gate_consts_8_v2.jl:130-134` was written to prevent.

**R13 — D-17's stated target (`corpus/`) is not executable. (probability HIGH ≈ 0.9 ·
impact LOW, once known)** Corpus bytes absent, anchor hashes `PENDING-FETCH`, and both physical
anchors are `sealed_holdout` (§F20). **Early warning:** immediate — `ls corpus/data`.
**Mitigation:** retarget to `test/test_images/` (real, 1028×1376, committed, non-sealed, already
used by the manuscript test suite). Flagged in Open Questions for confirmation, since it is a
target substitution.

**R14 — the interpolation-onset step is misread as SC2 signal. (probability MEDIUM ≈ 0.3 ·
impact MEDIUM)** `‖Δs‖₂` jumps 0.000 → 0.380 between shift 0 and 0.25 px purely because integer
offsets bypass the B-spline kernel (§C9b). A ladder including `λ = 0` would show a large, spurious
first step. **Early warning:** the probe curve's shape. **Mitigation:** ladder starts at
`λ = 0.25`; the integer-offset control arm (§F20b item 3) is carried into the real-image check;
the report annotates the step explicitly.

### H25. The SC2/SC3 structural tension (D-13) — stated in this document's own words

**The tension.** SC2 and SC3 are not independent criteria; SC2's mechanism is what makes SC3's
detection impossible.

The shipped OOD flag is a **summary-density Mahalanobis** test: it fires when an incoming summary
sits far from the cloud of summaries the net was *trained* on (`src/amortized/ood.jl:77-106`,
operating point `OOD_ID_QUANTILE = 0.95` at `:61`; the shipped bundle carries **only** that channel,
`src/registry.jl:157-159`, `docs/amortized.md:145-150`). Detection is therefore, by construction,
detection of *unfamiliarity*.

SC2 works by making mis-registration **familiar**. D-02 widens the shift prior from ±1 to ±3 px and
D-09 adds `ε ∈ [-0.02, 0.02]`, precisely so the net has *seen* mis-registered data and can widen its
`Δρ` posterior honestly in response. But a summary the net has seen is, definitionally, a summary
the density channel scores as in-distribution. **Every pixel of misalignment that SC2 teaches the
net to handle is a pixel of misalignment the density channel stops flagging.** The two criteria pull
on the same rope in opposite directions, and the harder SC2 succeeds, the more thoroughly SC3's
would-be detector is disabled.

D-13 resolves this by **changing what SC3 asserts**: not "mis-registered data is flagged" but
"mis-registered data is *covered*" — the `Δρ` credible interval still contains the truth. That is a
coherent and defensible resolution, because a covered interval **is** the operational absence of
silent overconfidence: the user is not told something false, they are told something appropriately
uncertain.

**The accepted cost, which the report must state plainly rather than bury:** badly mis-registered
data returns a **wide but unflagged** answer. The interval is honest, but the user is not told *why*
it is wide, and cannot distinguish "these two proteins genuinely have an ambiguous relationship"
from "your channels are misaligned by 3 px". The failure is honest but **silent** — and D-14's
breakdown curve is what converts that silence into a quotable operating limit ("coverage holds to
X px; beyond that, neither the interval nor the flag protects you").

**Does this create an unresolvable planning contradiction? No.** Three reasons:

1. The two criteria are measured by the **same statistic on the same runs** — per-rung coverage.
   SC2(b) (D-07b) and SC3-in-prior (D-13/D-14) are literally the same TOST result read twice (§E18).
   There is no experiment that must simultaneously satisfy two incompatible tests.
2. The only *gated* SC3 claim is **in-prior** coverage, which the SC2 ladder already produces. The
   beyond-prior region — the only place where the erosion of OOD detection actually bites — is
   explicitly **reported, not gated** (D-14).
3. The deferred "dedicated registration OOD channel" idea is recorded as the thing that *would*
   break the tension, and is deferred for an independent reason (it needs fitted nulls in a
   re-shipped bundle, contradicting D-01's no-reship posture).

The residual is a **documentation obligation, not a planning conflict**: the tension must appear in
the Phase-11 report and in the D-16 `docs/amortized.md` interpretation note, in plain language,
because a reader who knows the tool has an OOD flag will otherwise assume misalignment is covered by
it. Recommended one-liner for the docs note: *"A wide `Δρ` posterior under registration uncertainty
is the honest answer, not a flagged error — the OOD flag does not fire on misalignment the model was
trained to expect, and by design cannot."*

---

## Validation Architecture

> `workflow.nyquist_validation` is `true` in `.planning/config.json`, so this section is required.

### Test Framework

| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (`@testset`/`@test`) — root suite `test/runtests.jl` (1400+ lines), spike suite `spike/test/` |
| Config file | `Project.toml:42-46` (`[extras] Test`, `[targets] test = ["Test"]`); spike has no separate test target, files are `include`d |
| Quick run command (src) | `julia --project -e 'using Pkg; Pkg.test()'` |
| Quick run command (spike) | `julia --project=spike spike/test/<file>.jl` |
| Full suite command | `julia --project -t auto -e 'using Pkg; Pkg.test()'` |
| Existing gate machinery reused | `test/gate/sbc.jl`, `test/gate/harness.jl`, `spike/validation/{harness,sbc}.jl` |

### Success Criteria → Test Map

| SC | Behaviour | Test type | Automated command | Exists? |
|----|-----------|-----------|-------------------|---------|
| **SC1a** (dx/dy already latent) | `shift_dx`/`shift_dy` are prior-drawn, applied in stage 6, read off the NPE | regression (already covered) | `julia --project -e 'using Pkg; Pkg.test()'` (`runtests.jl` θ round-trip at :541-563) | ✅ |
| **SC1b** (`ε` in θ) | `sample_prior` returns 8 fields; `theta_prior_bounds()` returns 8; `ε` reaches stage 6 | unit | `julia --project -e 'using Pkg; Pkg.test()'` → new `@testset "chromatic ε (D-09)"` | ❌ **Wave 0** |
| **SC1c** (composed warp is one pass) | the stage-6 `tform` is a single `AffineMap`, not a `ComposedTransformation`, and `warp` is called once | unit | same suite; `@test A isa CoordinateTransformations.AffineMap` | ❌ **Wave 0** |
| **SC1d** (ε = 0 regression) | `ε = 0, shift = s` composed warp `==` legacy `Translation(s)` warp, and full `simulate_pair` `==` a pre-edit golden | **regression, exact `==`** | `julia --project=spike spike/test/test_stage6_regression.jl` | ❌ **Wave 0** (golden must be captured **pre-edit**) |
| **SC1e** (backward compat) | `simulate_pair(rng, θ_7field)` still runs and is bit-identical | unit | same | ❌ **Wave 0** |
| **SC1f** (net trains) | `d_in = 129`, `D = 8` estimator trains on 200 synthetic samples without error | smoke | `julia --project=spike spike/test/test_p11_smoke.jl` | ❌ **Wave 0** |
| **SC1g** (λ conditioning is live) | same `Z₁₂₈` at `λ = 0.25` vs `λ = 3.0` gives materially different Δρ posterior SD | **behavioural, blocks the ladder** | `julia --project=spike spike/test/test_lambda_ablation.jl` | ❌ **Wave 0** |
| **SC2a** (monotone widening) | Spearman(λ, Δρ-width) permutation test passes and `S ≥ SC2_SPEARMAN_FLOOR` | experiment (reported) | `julia --project=spike -t auto spike/validation/run_p11_ladder.jl` | ❌ **Wave 2** |
| **SC2b** (coverage honest) | per-rung TOST equivalence to 0.90 within ±3 pp, Holm across 7 rungs, **`all(<=(fwer), …)`** | experiment (reported) | same run | ❌ **Wave 2** |
| **SC3a** (in-prior coverage) | identical statistic to SC2b — **not a separate run** | experiment | same run | ❌ **Wave 2** |
| **SC3b** (breakdown curve) | first beyond-prior rung with `wilson_ci(...).hi < 0.90` | experiment (reported, not gated) | `julia --project=spike -t auto spike/validation/run_p11_breakdown.jl` | ❌ **Wave 2** |
| **D-02** (ρ_true attenuation) | paired frozen-vs-research RMSE / post_sd / coverage on one shared eval set | experiment (reported) | `julia --project=spike spike/validation/run_p11_attenuation.jl` | ❌ **Wave 2** |
| **D-17** (real-image transfer) | Δρ width monotone in injected shift on the positive arm; steeper than the negative arm | qualitative (reported) | `julia --project=spike spike/validation/run_p11_realimage.jl` | ❌ **Wave 3** |
| **D-12** (provenance) | pinned sha is the ε commit's first parent and its `simulator.jl` lacks `chromatic_eps` | shell assertion | `git rev-parse "<EPS>^"` + `git show "<SHA>":src/amortized/simulator.jl \| grep -c chromatic_eps` | ❌ **Wave 1** |
| **D-16/D-01** (no shipped change) | `Artifacts.toml` unchanged; `src/bayes.jl`/`colocalization.jl`/`LoadImages.jl`/`plot.jl` unchanged; suite green | shell + suite | the five commands in §F20c | ❌ **Wave 3** |
| **TOST direction** | `p11_tost_pass` uses `<=`, verified on a hand-computed vector | unit | root or spike suite | ❌ **Wave 0** |

### Sampling Rate

- **Per task commit:** `julia --project -e 'using Pkg; Pkg.test()'` for any task touching `src/`;
  `julia --project=spike spike/test/<touched>.jl` for spike-only tasks. Target < 30 s for the spike
  unit files (the root suite is minutes — run it on `src/`-touching commits only).
- **Per wave merge:** full root suite **plus** all `spike/test/` files.
- **Phase gate (before `/gsd:verify-work`):** full root suite green **and** the five §F20c
  provenance/decoupling commands producing empty output **and** every Wave-2/3 experiment script
  having written its report artifact.
- **Experiments are NOT tests.** The ladder/breakdown/attenuation/real-image runs are reported
  measurements with pre-registered pass criteria, executed once each on the frozen consts. They must
  not be wired into CI (a stochastic multi-minute experiment in a test suite is a flaky gate), but
  their *harness functions* must have fast fixture-scale unit tests (`SBC_FIX_*` pattern,
  `spike/validation/consts.jl:53-56`, `gate_consts_8_v2.jl:169-171`).

### Wave 0 Gaps

- [ ] `spike/validation/p11_consts.jl` — Tier-1 pre-registration (seed, ladders, N, band, α, Holm
      FWER, iteration allowance) + executable disjointness `@assert`s (§G22)
- [ ] `spike/test/test_p11_consts.jl` — asserts `P11_DEV_SEED ∉ _p11_forbidden()`, salt distinctness,
      ladder monotonicity, `N_PER_RUNG ≥ 271`, `(SBC_L + 1) % bins == 0`
- [ ] **Pre-edit golden fixture** `spike/test/fixtures/p11_stage6_golden.jld2` — captured **before**
      any stage-6 change (SC1d); a task whose ordering is load-bearing
- [ ] `spike/test/test_stage6_regression.jl` — SC1c/SC1d/SC1e, exact `==`
- [ ] Root-suite `@testset "chromatic ε (D-09)"` in `test/runtests.jl` — SC1b
- [ ] **Update ~17 literal-`7` θ assertions in `test/runtests.jl`** to
      `length(ProteinCoLoc.theta_prior_bounds())` (§A4 row 24) — R1
- [ ] `spike/test/test_p11_smoke.jl` — `d_in = 129, D = 8` train on 200 synthetic samples (SC1f);
      **must run before the 1.2 h datagen** (R10)
- [ ] `spike/test/test_lambda_ablation.jl` — SC1g, the R4 tripwire
- [ ] `spike/test/test_p11_tost.jl` — Wilson CI, TOST p-value, and the Holm **direction** (R7)
- [ ] `spike/test/test_p11_forced_theta.jl` — `merge(sample_prior(rng), forced)` preserves field
      order so `collect(values(θ))` stays row-aligned (§E15 item 2)
- [ ] Holm adjuster available in `spike/` — copy `sbc_holm_adjusted` (`test/gate/sbc.jl:352-362`)
      with attribution, or `include` the gate file; **no new dependency**

---

## Recommended Task Sequence

**Legend:** ⛔ = hard ordering constraint (violating it invalidates the result) · ∥ = parallelizable

| # | Step | Ordering | Parallel? |
|---|------|----------|-----------|
| **W0-1** | Write + commit `spike/validation/p11_consts.jl` (Tier 1 only) and `test_p11_consts.jl` | ⛔ **before every run**, D-04 | — |
| **W0-2** | Capture `PHASE11_BASE_SHA` and the **pre-ε sha** (`git rev-parse HEAD`) as recorded values | ⛔ **before any `src/` or simulator edit**, D-12 (R9) | ∥ with W0-1 |
| **W0-3** | Capture the **pre-edit golden fixture** for `simulate_pair` at 4 keys × fixed 7-field θ | ⛔ **before the stage-6 edit**, SC1d | after W0-2 |
| **W0-4** | Write the TOST/Wilson/Holm-direction helpers + `test_p11_tost.jl` | none | ∥ with W0-1..3 |
| **W1-1** | **THE D-11/D-12 COMMIT (single commit).** Compose the affine warp in `spike/simulator/forward.jl` and `src/amortized/simulator.jl`; add `CHROMATIC_PRIOR` + the 8th θ field to both `prior.jl` and `simulator.jl`; extend `theta_prior_bounds()`; defensive `ε` read in both `simulate_pair`s; `_theta_tuple` guard in `ood.jl`; `datagen.jl` 7→8 at `:170,:274,:393`; update ~17 literal-7s in `test/runtests.jl`; re-derive the WR-06 comment; **and** append named limit #8 to `docs/amortized.md` with the W0-2 sha | ⛔ **D-12 same commit as D-11**; ⛔ after W0-3 | — |
| **W1-2** | `spike/test/test_stage6_regression.jl` (exact `==` vs the W0-3 golden and vs legacy `Translation`) + full root `Pkg.test()` | ⛔ after W1-1 | — |
| **W1-3** | Verify the pinned sha: `git rev-parse "<EPS>^"` and `git show "<SHA>":…\|grep -c chromatic_eps` | ⛔ after W1-1 | ∥ with W1-2 |
| **W1-4** | `spike/test/test_p11_smoke.jl` — `d_in=129, D=8` trains on 200 synthetic samples (retires R10) | after W1-1 | ∥ with W1-2/3 |
| **W2-1** | **RUN THE PRE-FLIGHT PROBE** (D-06), ~10 min, simulator only | ⛔ **strictly before any training**, D-06 | — |
| **W2-2** | Append the **Tier-2** probe-derived block to `p11_consts.jl` (`SC2_SPEARMAN_FLOOR`, band justification, `Δρ_eq` coefficients, chosen imsize arm) — **append, never edit** | ⛔ after W2-1, ⛔ before W3-1 | — |
| **W2-3** | **Abort check:** if `S_probe < 0.9` or ladder span `< 0.02 Δρ_eq`, take the §C10(c) branch (one re-parameterised probe; then stop or proceed) | ⛔ gate on W2-1 | — |
| **W3-1** | λ-hierarchical datagen (50k, F5 mixture read from `gate_consts_8_v2.jl`) | ⛔ after W2-2 | — |
| **W3-2** | Train the research net (`d_in = 129, D = 8`), persist with `imsize` provenance | ⛔ after W3-1 | — |
| **W3-3** | **`test_lambda_ablation.jl` — the R4 tripwire** | ⛔ after W3-2, ⛔ **before W4-\*** | — |
| **W4-1** | SC2 ladder run (7 rungs × N=500) → widths, coverage, TOST, Holm, permutation Spearman. **Also yields SC3-in-prior** | ⛔ after W3-3 | — |
| **W4-2** | SC3 beyond-prior breakdown run (shift + ε ladders, λ pinned at 3.0 and 1.0) | after W3-3 | ∥ with W4-1 |
| **W4-3** | D-02 `ρ_true` attenuation measurement (frozen vs research on one shared eval set) | after W3-3 | ∥ with W4-1/2 |
| **W4-4** | Per-rung shrinkage / vacuity table for `shift_*` and `ε` (§E19) | after W4-1 | ∥ with W4-2/3 |
| **W5-1** | D-17 real-image check on `test/test_images/` (null-warp + integer-offset controls) | after W3-3 | ∥ with W4-\* |
| **W5-2** | Figures (ladder, breakdown curve, real-image width response, probe curves) | after W4-1/2, W5-1 | ∥ |
| **W6-1** | Write the Phase-11 report under `.planning/phases/11-…/` (probe, ladder, breakdown, attenuation, real-image, the D-13 tension in plain language, the vacuous-column caption) | after W4-\*, W5-\* | — |
| **W6-2** | D-16 interpretation section in `docs/amortized.md` (separate `##` section; scrupulous about which net each number came from) | after W6-1 | ∥ with W6-3 |
| **W6-3** | Run and record the five §F20c decoupling/provenance verification commands | after W1-1 | ∥ |

**The two ordering constraints that must not be traded away:** ⛔ **W2-1 (probe) strictly before
W3-1/W3-2 (training)** — D-06's entire value is being the cheap gate before the expensive step; and
⛔ **W1-1 is one commit containing both the `ε` edit and the `docs/amortized.md` limit** — D-12.
A third, easily missed: ⛔ **W0-3 (golden fixture) before W1-1**, or the D-10 regression check
becomes unverifiable.

---

## Common Pitfalls

### Pitfall 1: Drawing λ independently of the shift
**What goes wrong:** SC2 shows no width response; the phase concludes "below resolution" and
possibly burns its single D-04 iteration on a re-parameterised ladder that cannot help.
**Why it happens:** the sampler reads naturally as "draw θ, draw λ, train" — the conditional
dependence is easy to drop.
**How to avoid:** `λ` first, then `shift ~ Uniform(-λ, λ)`; assert `abs(shift) <= λ` per sample.
**Warning signs:** the λ-ablation test (W3-3) shows equal posterior SD at λ = 0.25 and λ = 3.0.

### Pitfall 2: Standardizing after appending λ
**What goes wrong:** either a loud error (`Loader._row_partition` requires exactly 128 rows,
`spike/data/loader.jl:126`; `src`'s twin requires an even count, `summary.jl:91-92`) or, if someone
"fixes" the partition, a silently z-scored λ that couples the conditioning variable to the training
pool's λ distribution.
**How to avoid:** standardize the 128 rows with the frozen `zt`, **then** `vcat` the λ row encoded
by a fixed, data-independent affine map.

### Pitfall 3: Copying `sbc_holm_pass` for the equivalence family
**What goes wrong:** the verdict inverts — the phase reports PASS exactly when coverage is *not*
equivalent to nominal.
**How to avoid:** §E16's `all(<=(fwer), sbc_holm_adjusted(p_tost))`, a distinct function name, and a
hand-computed unit test.

### Pitfall 4: Reading the interpolation-onset step as misalignment sensitivity
**What goes wrong:** the probe's 0 → 0.25 px jump (0.000 → 0.380 measured) inflates the apparent
effect and a λ = 0 rung makes the ladder look strongly monotone for the wrong reason.
**How to avoid:** ladder starts at λ = 0.25; annotate the step; carry the integer-offset control
arm into the real-image check.

### Pitfall 5: Computing shrinkage against the pooled marginal shift prior
**What goes wrong:** under hierarchical λ the pooled `sd(shift_dx)` exceeds `λ_r/√3` at low rungs,
so `shrinkage = post_sd/prior_sd` comes out well below 0.95 and the shift columns *look* informative
— contradicting D-05 and the F3 discipline.
**How to avoid:** compute `prior_sd` **within each rung**.

### Pitfall 6: Pinning the sha after committing
**What goes wrong:** the D-12 limit points at the ε commit itself, so `git show <SHA>:…/simulator.jl`
already contains `chromatic_eps` and the provenance-by-reference claim is false.
**How to avoid:** W0-2 captures it as a recorded value; W1-3 verifies with
`git rev-parse "<EPS>^"`.

### Anti-Patterns to Avoid

- **Two sequential `warp` calls.** The whole point of D-10; a second interpolation pass widens the
  posterior for a reason unrelated to misalignment and contaminates SC2 directly.
- **Routing the research net through `_train_grid_pipeline`.** It cannot express `D = 8`
  (`src/amortized/train_npe.jl:136` passes no `D`; no `D` keyword exists at :114-124) or a 129-row
  input, and its `datagen` seam allocates a `7×N` θ (`datagen.jl:170`).
- **Mirroring D-02's widened `SHIFT_PRIOR` into `src/`.** D-11 mirrors `ε`, not the widening; a
  widened src prior would desynchronise `theta_prior_bounds()[5:6]` from the shipped bundle's frozen
  `θzt.lo/hi`.
- **Deleting orphaned datagen caches.** They are the byte-level record of the shipped net's training
  pool — exactly what D-12 is preserving by reference.
- **Opening the corpus `sealed_holdout` split.** It would burn the Phase-16 blind evaluation for a
  qualitative check that `test/test_images/` serves better.
- **Adding capacity if training stalls.** Spike 012 is PARTIAL: capacity redistributes the drift and
  overfits near its ceiling.
- **Quoting a clean KS/rank result on `shift_dx`, `shift_dy` or `ε` as calibration evidence.**
  `gate-8x8-amended.md:46,125`; `07-CALIBRATION-FINDINGS.md:86-88`.

---

## Don't Hand-Roll

| Problem | Don't build | Use instead | Why |
|---|---|---|---|
| Affine composition / recentring | manual 2×2 matrix + offset arithmetic | `Translation ∘ recenter(LinearMap(M), c)` (CoordinateTransformations 0.6.4) | `compose` methods collapse to one `AffineMap` (`affine.jl:139-165`), which is what guarantees the single `warp` pass **and** the bit-exact ε = 0 identity (§A3) |
| Image centre | a literal `128.5` or `size(x)/2` | `map(ax -> (first(ax)+last(ax))/2, axes(x))` (== `ImageTransformations.center`, `ImageTransformations.jl:80-81`) | `size/2` is wrong for `OneTo` axes (gives 128.0, not 128.5) and silently off-centres the scale |
| Credible-interval width | manual `quantile` calls | `interval_width(draws; probs)` → NeuralEstimators `interval` (`spike/npe/infer.jl:148-151`) | already the house rule ("Never hand-roll quantiles", `infer.jl:146`) |
| Holm adjustment | a fresh implementation | `sbc_holm_adjusted` (`test/gate/sbc.jl:352-362`) / `holm_adjusted` (`gate_consts_8_v2.jl:361-371`) | already unit-tested to agree; valid under arbitrary dependence |
| Shrinkage / vacuity flag | a new diagnostic | `sbc_shrinkage` (`test/gate/sbc.jl:243-250`) + `SBC_VACUOUS_SHRINKAGE` | F3-mandated, already generic in the column count |
| Reproducible per-sample RNG | `Random.seed!` in a loop | `Philox4x(UInt64, (seed ⊻ SALT, idx))` (`spike/data/seeding.jl:75-96`) | thread- and order-independent; the project's reproducibility contract |
| Atomic artifact write | ad-hoc `jldsave` | the `.tmp` → reopen-assert → `mv(...; force=true)` idiom (`spike/npe/train_npe.jl:160-180`, `src/amortized/persist.jl:66-82`) | a crash leaves a discardable `.tmp`, never a torn artifact |
| Image-size mixture | retyping the F5 tuple | read `SBC_IMSIZE_SET`/`WEIGHTS` from `test/gate/gate_consts_8_v2.jl:152-153` | the binding invariant is train-joint == eval-joint; `scripts/train_grid8_amended_v2.jl:34-35` establishes the never-retype rule |
| ROC/AUC, if the breakdown curve wants separability | a new implementation | `roc_auc` / `id_threshold` / `youden_j` (`src/amortized/ood.jl:268-305`) | hand-rolled, dependency-free, tie-aware Mann–Whitney |

**Key insight:** every primitive Phase 11 needs already exists in this repo, usually with a comment
explaining why the obvious hand-rolled version is wrong. The dangerous novelty is not any single
function — it is the *composition* (λ-conditional sampling, inverted Holm direction, per-rung
shrinkage) where reusing a house function with the wrong semantics is easy and silent.

---

## Environment Availability

| Dependency | Required by | Available | Version | Fallback |
|---|---|---|---|---|
| Julia | everything | ✓ | 1.12.6 | — |
| CPU threads | datagen (`-t auto`) | ✓ | **32** (`nproc`) | serial (≈ 6× slower) |
| `CoordinateTransformations` | D-10 composed warp | ✓ | 0.6.4 (both manifests) | — |
| `ImageTransformations` | `warp`, `center` | ✓ | 0.10.3 | — |
| `Interpolations` | `BSpline(Linear())` | ✓ | 0.16.3 | — |
| `NeuralEstimators` | NPE | ✓ | 0.2.1 | — |
| `Flux` | backend | ✓ | 0.16.10 (compat pin) | — |
| `Random123` | Philox | ✓ | in both `Project.toml`s | — |
| `StatsBase` (`corspearman`) | SC2 statistic | ✓ | `spike/Project.toml:17` | — |
| Real microscopy images (D-17) | real-image check | ✓ | `test/test_images/{positive,negative}/*_c{1,2,3}.tif`, 1028×1376 `RGB{N0f16}`, committed | — |
| Shipped 8×8 bundle (D-02 attenuation baseline) | frozen-vs-research comparison | ✓ | `artifacts/amended_v2/grid_8/{npe,ratio,ood_nulls,gate_report}_8.jld2` present on disk (gitignored); tree-sha1 `90e6b63a…` | lazy download from the GitHub Release |
| GPU / CUDA | optional accelerator | not probed | — | CPU is the baseline (CLAUDE.md); datagen dominates, so GPU buys little |
| Network | **not required** | — | — | — |
| `corpus/` image bytes | D-17 **as literally written** | ✗ | `corpus/data/` empty; anchor `sha256 = PENDING-FETCH`; anchors are `sealed_holdout` | **`test/test_images/`** (§F20) |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** `corpus/` bytes → `test/test_images/` (see Open Question 1).

---

## Security Domain

> `security_enforcement` is not set in `.planning/config.json`; treated as enabled.

### Applicable ASVS categories

| Category | Applies | Standard control |
|---|---|---|
| V2 Authentication | no | no auth surface; local scientific computation |
| V3 Session management | no | no sessions |
| V4 Access control | **partially** | the corpus `sealed_holdout` split is an *epistemic* access control (`corpus/config.jl:55-58`, `manifest.jl:156`). Phase 11 must not open it — §F20 blocker 3 |
| V5 Input validation | **yes** | `simulate_pair` validates every θ field at entry (`forward.jl:111-134`, `simulator.jl:186-202`, tagged `T-02-IV: untrusted parameter vector`). **The new `chromatic_eps` must join those guards** — at minimum `isfinite` and `-1 < ε` (so `1 + ε > 0` and `s = 1/(1+ε)` is finite and positive; a non-positive `1+ε` would mirror or collapse the image) |
| V6 Cryptography | **yes, indirectly** | `Artifacts.toml` integrity is a SHA-256 download hash + a `git-tree-sha1` verified before load (`registry.jl:101-108`, T-7-01). Phase 11 changes neither and must not weaken either. `corpus/config.jl:97-100` mandates cryptographic SHA-256 (never Base `hash`) for download integrity — relevant only if the planner overrides Open Question 1 and fetches the corpus |
| V12 File handling | **yes** | D-17 reads TIFFs. The path is the frozen `load_tiff` (`src/LoadImages.jl:214-217`) on **committed, in-repo** bytes — no user-supplied path, no upload, no archive extraction |
| V14 Configuration | **yes** | the pre-registration consts file is a security-of-the-scientific-claim artifact: `@assert`-enforced seed disjointness (§G22) is the control that prevents a snooped stream |

### Known threat patterns for this stack

| Pattern | STRIDE | Mitigation in place / required |
|---|---|---|
| Loading a tampered or content-mismatched model artifact | Tampering | `_verify_tree_sha1` before every load (`registry.jl:101-108`); Phase 11 **must not** add a bypass |
| Deserialization of an arbitrary object graph | Tampering / RCE | `src` persists `Flux.state` + arch metadata, **not** whole objects (`persist.jl:20-36`, T-7-01). **Note:** the *spike* trainer still persists the whole estimator (`spike/npe/train_npe.jl:163-174`); acceptable for a research-lane artifact that is never shipped, but the report should say the Phase-11 net is **not** a distributable artifact |
| Reading a stale/orphaned training cache as if current | Tampering | `open_or_invalidate` errors on a stored-vs-recomputed digest mismatch and names the diverging file via `subhashes` (`datagen.jl:368-385`) |
| Non-finite / adversarial θ reaching the simulator | DoS / Tampering | entry validation (`forward.jl:111-134`); **extend to `ε`** |
| Seed reuse invalidating a pre-registered claim | Repudiation | executable `@assert !(P11_DEV_SEED in _p11_forbidden())` (§G22) |
| Unverified network fetch of corpus bytes | Tampering | avoided entirely by targeting `test/test_images/`; if overridden, `fetch_verified` hard-errors on hash mismatch (`corpus/fetch.jl:80`, `config.jl:97-100`) |

**No new attack surface is introduced by Phase 11 as scoped** (no API change, no artifact change, no
network dependency).

---

## State of the Art

| Old approach (in this repo) | Current approach | When changed | Impact on Phase 11 |
|---|---|---|---|
| Strict point-null SBC at M = 2000 on all 8 columns | Equivalence (TOST) for non-identified nuisances; strict point null for targets only | spike 011 → `07-NUISANCE-SBC-SPEC-DRAFT.md` | D-08 follows it; §E16 implements it, and inverts the Holm direction accordingly |
| Uncorrected all-8 α = 0.05 conjunction (FWER 33.7 %) | Holm–Bonferroni at FWER 0.05 | `07-GATE-AMENDMENT.md` A1 / `gate_consts_8_v2.jl:94-100` | reuse `sbc_holm_adjusted`; **direction inverts** for equivalence |
| Raw KS on `ρ_true` (undefined under prior atoms) | **Randomized ranks** for atom-carrying columns | spike 010 ADDENDUM | `ρ_true` is still the only atom-carrying column; the Phase-11 SBC-style reporting must apply it |
| Scalar `SBC_IMSIZE = (256,256)` | F5 realistic mixture, train-joint == eval-joint, provenance persisted | `gate_consts_8_v2.jl:126-166`, `persist.jl:143-158` | Phase 11 must read the mixture from the frozen file and persist `imsize` provenance |
| Plain `ZScoreTransform` on raw θ | `BoundedThetaTransform` (logit-of-prior-box ∘ z-score) | F2 remedy, `src/amortized/architecture.jl:44-104` | **spike still uses the plain z-score** (`spike/npe/train_npe.jl:80`). Open Question 4 |
| KDE + `quadgk` Bayes factor as the reference | simulation-based discrimination (AUC 0.994) | spike 014 | not touched by Phase 11 |
| Whole-object model persistence | `Flux.state` + arch metadata | Pitfall 4 / T-7-07, `persist.jl` | the spike trainer still uses the old form; fine for a non-shipped research net, state it |

**Deprecated / superseded readings to avoid repeating:**
- "The NPE is overconfident on `ρ_true`" — **retracted** (spike 010 ADDENDUM); it is a prior-atom artifact.
- "Raising flow capacity fixes the nuisance drift" — **PARTIAL/false** (spike 012): it redistributes.
- "Truncate the μ-prior to fix `ρ_true` SBC" — costs high-|ρ| accuracy and ADVI comparability (spike 013).

---

## Assumptions Log

| # | Claim | Section | Risk if wrong |
|---|---|---|---|
| A1 | Appending λ as a 129th input row is the cleanest D-03 mechanism | §D12 | LOW — verified that `build_estimator` is width-agnostic and `train` is generic; an alternative (separate conditioning branch) exists if it underperforms |
| A2 | `0.5 · S_probe` is a defensible attenuation allowance for `SC2_SPEARMAN_FLOOR` | §C10(a) | MEDIUM — the factor is judgment, not derivation. Must be pre-registered before the run; the report must name it as judgment (the same honesty `07-NUISANCE-SBC-SPEC-DRAFT.md:95-101` applies to δ = 0.10) |
| A3 | `N_PER_RUNG = 500` is enough for a decidable ±3 pp band | §E16 | LOW — `N_min = 271` is exact arithmetic; 500 gives slack for `p̂ ∈ [0.878, 0.922]` |
| A4 | Training ≈ 1.5 h wall clock at 50k on the F5 mixture | §D14 | MEDIUM — extrapolated from `gate_consts_8_v2.jl:148` (0.0835 s/pair at 32 thr) and `07-05-SUMMARY.md`; the 07-05 run at the *default* mixture was killed at 45 min, so the mixture choice matters |
| A5 | `ε` is NOT λ-scaled (fixed `Uniform(-0.02,0.02)`) | §D13 | MEDIUM — follows D-09's literal text and D-03's singular "the shift prior's half-width", but it means SC2's ladder is registration-only. Open Question 3 |
| A6 | `test/test_images/` are the "real images used by the manuscript pipelines" | §F20 | MEDIUM — they are real, 1028×1376, committed, and exercised by `runtests.jl:86,166`, but D-17 names `corpus/`. Open Question 1 |
| A7 | `src/amortized/simulator.jl:94` `SHIFT_PRIOR` stays at ±1 | §A5 | MEDIUM — D-11 says mirror `ε`; D-02 scopes the widening to the research net. Open Question 2 |
| A8 | The indicative probe's `Δρ_eq` scale transfers to the pre-registered probe | §C9b | MEDIUM — one θ, one imsize, R = 12. The real probe re-measures everything; the indicative numbers are directional only |
| A9 | Exact `==` holds for the ε = 0 regression under future dependency bumps | §A3 | LOW — measured true today; a StaticArrays/Interpolations bump could change reduction order. Documented fallback `1e-13` |
| A10 | GPU buys little because datagen dominates | §D14 | LOW — datagen is ~1.2 h of a ~1.5 h budget and is CPU-bound |
| A11 | λ = 3.0 is the right read-time level for the beyond-prior breakdown curve | §E18 | MEDIUM — D-14 does not specify it. Open Question 5 |

---

## Open Questions for the Planner

**1. D-17's stated target (`corpus/`) is not executable; `test/test_images/` is. Confirm the
substitution.**
Evidence: `corpus/data/` is empty; both `physical-primary` rows carry `sha256 = PENDING-FETCH`
(`corpus/anchor_rows.jl:85`); and — decisively — both are `split = sealed_holdout`, which
`corpus/config.jl:55-58` reserves for Phase-16 blind evaluation and which
`corpus/manifest.jl:156` gates behind `open_sealed_holdout(df; reason)`. The 30 non-sealed rows are
`simulated-secondary` CBS images with empty hashes. Meanwhile `test/test_images/{positive,negative}/`
holds six **real** 1028×1376 `RGB{N0f16}` microscopy TIFFs, committed to git, already used by
`test/runtests.jl:86,166` — i.e. literally D-17's phrase "the existing real images already used by
the manuscript pipelines". **I could not resolve from code whether the user meant `corpus/`
specifically or "the real images we already have"; the discussion log (`11-DISCUSSION-LOG.md:300-302`)
says the corpus was "identified during the answer", which suggests an in-conversation guess rather
than a verified target.** Recommend `test/test_images/`; escalate if the planner disagrees.

**2. Does D-11's mirror include D-02's widened `SHIFT_PRIOR`, or only `ε`?**
D-11's text is "Mirror the `ε` extension into `src/amortized/simulator.jl`"; D-02 scopes ±3 px to
"the research net". I recommend leaving `src/amortized/simulator.jl:94` at `Uniform(-1.0, 1.0)`,
because widening it would desynchronise `theta_prior_bounds()[5:6]` from the shipped bundle's frozen
`θzt.lo`/`θzt.hi` (a semantic drift that is currently harmless only because `θzt` is persisted) and
would be a larger provenance event than D-12's named limit describes. **Not resolvable from the
decision text alone.**

**3. Should there be a SECOND conditioning input for the chromatic level?**
D-03 says "the registration-uncertainty level (**the shift prior's half-width**)" — singular, and
explicitly about the shift. So I designed one conditioning input, with `ε` as a fixed marginalized
nuisance (§D13). The consequence is that **SC2's monotone-widening claim is about registration only**;
`ε`'s dose-response appears only in the SC3/D-14 fixed-injection curve. D-09's rationale ("the ±0.02
range matches the ±3 px translation headroom so both misalignment axes are comparable") could be read
as wanting them on one joint ladder, which would need either a second input or an λ-scaled `ε` prior
(contradicting D-09's literal `Uniform(-0.02, 0.02)`). **Flagging rather than choosing.**

**4. Should the Phase-11 spike trainer adopt the `BoundedThetaTransform` (F2 remedy)?**
`spike/npe/train_npe.jl:80` still uses a plain `ZScoreTransform`; the F2 bounded/logit transform
exists only in `src/` (`src/amortized/architecture.jl:88-104`). The measured F2 defect was ~2 % of
posterior mass leaking outside `label_efficiency`'s box with a sign-matched bias
(`architecture.jl:50-52`). A research net trained without it inherits that defect, which would
contaminate the `ρ_true` attenuation comparison (R6) against a shipped net that *has* it. Porting it
is ~50 lines (the struct + two `StatsBase` methods) and would make the comparison cleaner — but it is
a model change not named in any decision, and D-13's precedent is that unnamed model changes are how
credibility gets spent. **Recommend porting it and stating the port explicitly**, but this is the
planner's call.

**5. At what λ should the beyond-prior breakdown curve be read?**
D-14 says "the point where 90% coverage first drops below nominal" but does not fix the read-time λ,
and the answer differs materially: at λ = 3.0 (user declares maximal uncertainty) coverage will hold
further out than at λ = 1.0 (user under-declares). I recommend λ = 3.0 as the primary curve — a
breakdown there is the more damning, more quotable number — with λ = 1.0 reported as a secondary,
more realistic user story. **Not settled by the decisions.**

**6. Are Phase-11 requirement IDs wanted in `REQUIREMENTS.md`?**
`ROADMAP.md:235` says `**Requirements**: TBD` and no Phase-11 IDs exist. The plan can proceed against
SC1/SC2/SC3 directly, but every other phase in this project has requirement IDs and the verification
tooling may expect them. **Cheap either way; flagging so it is a decision rather than an omission.**

**7. Does the SC2 ladder run at 256², the F5 mixture, or both?**
F5's binding invariant is train-joint == eval-joint, so the ladder must match whatever the research
net trains on. Training on the F5 mixture is the scientifically right answer (transferable to the
real-data regime, and the ε response is ~3.4× larger at the anchor) and costs ~1.2 h datagen; a 256²
run is ~6 min but re-creates exactly the compute-budget artifact `gate_consts_8_v2.jl:130-134` was
written to eliminate. **Recommend the F5 mixture**, with the choice recorded in the Tier-2 consts
block after the probe (§G23). Flagged because it is the single largest budget decision in the phase.

**8. What is the exact meaning of "provably unmodified manuscript pipelines" under D-11?**
`ROADMAP.md:126` defines the decoupling proof as `git status` clean on `src/`, `bayes.jl`,
`colocalization.jl`. D-11 deliberately edits `src/amortized/simulator.jl`, so the literal claim is no
longer available. §F20c proposes the narrower, checkable claim (v1.0 pipeline files byte-unchanged;
only `src/amortized/` touched) with five concrete commands. **The planner should confirm this
narrower wording is acceptable, because the report has to state it and a reviewer will compare it
against the ROADMAP text.**

---

## Sources

### Primary (HIGH confidence — read directly in this session)

- **Repository sources** (all cited with file:line throughout): `spike/simulator/{prior,forward}.jl`,
  `src/amortized/{simulator,datagen,architecture,train_npe,persist,ood,summary,pipeline}.jl`,
  `src/{registry,LoadImages}.jl`, `Artifacts.toml`, `spike/{validation,npe,data,contract}.jl` files,
  `test/{runtests.jl,gate/*.jl}`, `corpus/{config,load,manifest,anchor_rows,fetch}.jl`,
  `corpus/manifest.csv`, `docs/amortized.md`, `scripts/train_grid8_amended_v2.jl`,
  `Project.toml`, `Manifest.toml`, `spike/Project.toml`, `spike/Manifest.toml`.
- **Installed package sources**: `CoordinateTransformations 0.6.4` (`src/{core,affine}.jl`),
  `ImageTransformations 0.10.3` (`src/{warp,autorange,ImageTransformations}.jl`),
  `NeuralEstimators 0.2.1` (`Project.toml`, `src/train.jl`, `src/Estimators/PosteriorEstimator.jl`).
- **Executed measurements (this session, 2026-07-25, Julia 1.12.6, `spike/Project.toml`):**
  `simulate_pair` timings (44.9 ms @ 256², 164.3 ms @ 512²); composed-affine ↔ legacy `Translation`
  bit-equality at ε = 0 (4 shift pairs, plus 12-key full-pipeline summary equality); the paired shift
  / ε / Δρ sensitivity sweeps and the independent-seed noise floor; `nproc = 32`; TIFF dimensions and
  eltypes; the repo-wide `0x` seed-literal inventory.
- **Planning artifacts:** `11-CONTEXT.md`, `11-DISCUSSION-LOG.md`, `.planning/{ROADMAP,STATE,REQUIREMENTS,config.json}`,
  `07-{CALIBRATION-FINDINGS,GATE-AMENDMENT,NUISANCE-SBC-SPEC-DRAFT,CONTEXT}.md`,
  `gate-8x8-amended.md`, `test/gate/gate_consts_8_v2.jl`.
- **Project skill:** `.claude/skills/spike-findings-proteincoloc/SKILL.md` and
  `references/sbc-calibration.md` (randomized ranks, M = 2000 over-power, capacity/truncation levers).

### Secondary (MEDIUM confidence)

- `07-05-SUMMARY.md` / `07-06-SUMMARY.md` training timings — recorded by prior runs on (presumably)
  this machine; used as budget anchors and cross-checked against my own measurements.
- The `E[cost] = 0.08346 s/pair` F5 figure (`gate_consts_8_v2.jl:148`) — a documented derivation, not
  re-measured here.

### Tertiary (LOW confidence — flagged, not relied upon)

- The `0.5` attenuation factor for `SC2_SPEARMAN_FLOOR` and the λ = 3.0 breakdown read-level are
  **my proposals**, not repo facts or established practice. Both are in the Assumptions Log and the
  Open Questions.

**No external web sources were consulted.** Everything above is either read from this repository,
read from the installed package sources, or measured by execution — which is the appropriate
evidence standard for a phase whose entire content is repo-internal mechanics.

---

## Metadata

**Confidence breakdown:**

| Area | Level | Reason |
|---|---|---|
| Simulator surgery (§A) | **HIGH** | Composition semantics read from installed source; ε = 0 bit-equality *measured* through the full pipeline; the arity ripple enumerated by exhaustive grep |
| Shipped-bundle compatibility (§A5) | **HIGH** | Traced the whole load path; `arch.D` and `θzt.lo/hi` are demonstrably read from disk, not recomputed |
| Provenance (§B) | **HIGH** | `cache_hash`/`Artifacts.artifact_hash` code paths read in full; independence established by absence of any call edge |
| Probe design + feasibility (§C) | **MEDIUM-HIGH** | Metric and arithmetic are sound and the indicative run is real, but it used one θ, one image size, R = 12 |
| Training mechanics (§D) | **MEDIUM-HIGH** | The 129-row/`D=8` route is verified against the installed API and the repo's own width-agnostic contract; the hierarchical sampler is a design proposal (correct, but untested here) |
| Evaluation statistics (§E) | **MEDIUM-HIGH** | TOST/Wilson/permutation arithmetic is exact; the Holm-direction trap is verified against source; `SC2_SPEARMAN_FLOOR`'s attenuation factor is judgment |
| Real-image check (§F) | **MEDIUM** | The corpus blockers are verified facts; the retarget is a recommendation the user has not confirmed |
| Pre-registration (§G) | **HIGH** | Seed inventory is exhaustive; the pattern is copied from two existing frozen files |
| Risks (§H) | **MEDIUM-HIGH** | Probabilities are calibrated judgment; every mechanism is code-cited |

**Research date:** 2026-07-25
**Valid until:** ~2026-08-24 (30 days). The repo-mechanics findings are stable as long as the pinned
manifests are unchanged; re-verify the ε = 0 bit-equality (§A3) if `CoordinateTransformations`,
`ImageTransformations`, `Interpolations` or `StaticArrays` is bumped.
