# Phase 11: Registration and Chromatic Uncertainty as Latent - Discussion Log

**Date:** 2026-07-25
**Mode:** discuss (default, interactive)
**Areas discussed:** Retraining posture, Identifiability strategy, Chromatic warp scope,
Mis-registration handling (all four selected by the user)

---

## Pre-discussion scouting finding

Before any question was asked, scouting established that **ROADMAP SC1 is already nominally
satisfied**: `shift_dx`/`shift_dy` are sampled from `SHIFT_PRIOR = Uniform(-1, 1)`
(`spike/simulator/prior.jl`), applied in `spike/simulator/forward.jl` stage 6, mirrored in
`src/amortized/simulator.jl`, and emitted by the NPE (`src/amortized/ood.jl:129-130`). This was
surfaced to the user rather than planning redundant work, and reframed the phase around three
genuine gaps: no chromatic term, the SC2 sweep never run, SC3 unaddressed.

Two structural tensions were identified as the basis for the gray areas:

1. **Retraining vs. the GO decision** — SC1 demands retraining; the Phase-7 GO (Option A,
   2026-07-24) locked "no further training, no gate iteration", and the shipped bundle is a frozen
   content-hashed Artifact.
2. **Identifiability** — `docs/amortized.md` named limit #4 establishes that the fixed
   patch-correlation summary structurally *cannot* constrain the shift nuisances.

---

## Area 1: Retraining posture

### Q1. What is the retraining posture, given SC1 demands retraining but the Phase-7 GO locked "no further training"?

Options presented:
- **Retrain in research lane** — train a Phase-11 research net in `spike/` on a fresh DEV seed
  asserted disjoint from `PROD_SEED_V2` / `VAL_MASTER_SEED` / `NPE_MASTER_SEED`; prove SC1-SC3
  there; shipped artifact and GO untouched.
- Reopen the ship gate — retrain and re-pre-register a new shipped bundle.
- No retraining — characterise the frozen net's behaviour under registration error only.

**Selected: Retrain in research lane.** → D-01

### Q2. What shift range should the Phase-11 research net be trained on? (current `SHIFT_PRIOR` = `Uniform(-1, 1)` px)

Options presented:
- **Widen to ±3 px** — covers chromatic aberration plus stage/filter-cube drift; gives the SC2 sweep
  real headroom (0.25 → 3). Cost: a wider nuisance prior adds marginal spread that may attenuate
  `ρ_true` inference — must be measured against the frozen net, not assumed.
- Keep ±1 px — minimal delta, comparable to the frozen net.
- Widen to ±5 px — maximal headroom, far outside realistic registration error.

**Selected: Widen to ±3 px.** → D-02

### Q3. How should the SC2 monotone-widening sweep be constructed?

Options presented:
- **Condition on uncertainty level** — append the registration-uncertainty level (shift prior
  half-width) to the network input and amortize over it; one training run gives a continuous
  controlled sweep. Cost: input dimension changes, so the net is not drop-in comparable to the
  shipped read surface.
- Train a ladder of nets — exactly comparable input surface, N training runs, discrete sweep.
- Post-hoc injection only — cheapest, but confounds prior width with data mismatch.

**Selected: Condition on uncertainty level.** → D-03

### Q4. How should Phase 11's SC2/SC3 thresholds be handled, given Phase 7's gate was amended twice?

Options presented:
- **Pre-register, allow one iteration** — lock thresholds in a Phase-11 consts file with a fresh DEV
  seed *before* any run, and state up front that this is a research-lane result permitting one
  documented iteration. Names the allowance in advance instead of amending after the fact.
- Strict pre-registration, no iteration.
- Exploratory, no pre-registration.

**Selected: Pre-register, allow one iteration.** → D-04

---

## Area 2: Identifiability strategy

**Framing established before questioning:** `gate-8x8-amended.md:42` names
`vacuous_params = [spillover, autofluorescence, shift_dx, shift_dy, noise]` at shrinkage 0.92-1.03 —
the posterior IS the prior. Line 46: `shift_dx`'s clean KS p = 0.883 is *meaningless*. Line 125: any
future SBC pass on those columns must be read as vacuous.

### Q1. Should Phase 11 make registration identifiable?

Options presented:
- **Marginalize-only, no ID claim** — keep the summary unchanged; widening the shift prior propagates
  registration uncertainty into `ρ_true`/`Δρ` as honest marginal spread; shift columns stay vacuous
  and are reported as such. Cost: Phase 11 can never say "we estimate the registration error", only
  "we are honest about not knowing it".
- Extend summary to identify shift (lag / cross-correlation features) — breaks summary reuse,
  invalidates the shipped 128-dim read surface, OOD reference, `local_coloc_map` and artifact.
- Both, staged — marginalize-only first, then lag features; ~double effort.

**Selected: Marginalize-only, no ID claim.** → D-05

### Q2. Should a pre-flight sensitivity probe run before the training run?

Risk raised: at grid 8 on 256² images each patch is ~32 px, so ±3 px is ~9% displacement; after PSF
blur the summary may barely move and SC2 would fail *flat* — below resolution, not wrong design.

Options presented:
- **Yes — simulator-only probe first** — sweep shift 0→3 px through `forward.jl` at fixed θ and
  measure summary displacement vs that induced by a known `Δρ`; minutes not hours; gives SC2 a
  defensible threshold instead of a guess.
- No — train first, measure after.
- Probe the frozen net instead — costs no training but confounds widening with misspecification
  (the frozen net was trained on ±1 px).

**Selected: Yes — simulator-only probe first.** → D-06

### Q3. What is the SC2 pass criterion?

Options presented:
- **Width increases AND coverage holds** — monotone non-decreasing `Δρ` width across the ladder AND
  90% CI coverage within nominal at every rung. Part (b) is what makes widening honest rather than
  merely larger.
- Width increases only — cannot distinguish honest widening from variance inflation.
- Width tracks the oracle — strongest claim, but the oracle and its tolerance both need defending.

**Selected: Width increases AND coverage holds.** → D-07

### Q4. How should the per-rung coverage check be tested, given the M=2000 over-power failure?

Context raised: Phase-7 spikes found "a ~0.06-SD marginal location drift on non-identified nuisances
at the M=2000 over-power edge" — a literal p-gate at large M rejects on drift too small to matter,
and a per-rung check multiplies that exposure.

Options presented:
- **Equivalence test per rung + Holm** — TOST-style equivalence to nominal with a pre-registered
  tolerance band, Holm-corrected across rungs; the exact remedy the Phase-7 spikes recommended.
  Passing means "coverage is provably close to nominal", not "we failed to reject".
- Coverage curve with CI band — still a failure-to-reject; weakens as M shrinks.
- Strict point-null, small M — reads as tuning M until the test passes.

**Selected: Equivalence test per rung + Holm.** → D-08

---

## Area 3: Chromatic warp scope

**Framing established:** `forward.jl` stage 6 is translation-only; there is no chromatic term
anywhere in the simulator. Lateral chromatic aberration is predominantly a radial magnification
difference — zero at image centre, maximal at corners — so a global `dx`/`dy` **cannot** represent
it. It is a distinct failure mode, not a reparameterisation, and it bites hardest in exactly the
regime ProteinCoLoc targets: patch-wise correlation over a field of view.

### Q1. Add the chromatic warp, or stay translation-only?

Options presented:
- **Add 1-param radial scale** — channel 2 scaled by (1+ε) about the image centre; the only
  spatially-varying misalignment in the model, and the natural bridge to Phase 12's per-region map.
  Cost: one more vacuous θ column, one more simulator stage, two-axis probe.
- Translation-only, defer chromatic — smallest scope; ROADMAP SC1's parenthetical goes unmet.
- Add ε, probe-gated — decides on evidence but leaves Phase 11 scope unfixed until the probe lands.

**Selected: Add 1-param radial scale.** → D-09

### Q2. What prior should ε get?

Magnitude table presented (256² field, half-diagonal ~181 px, corner displacement ≈ ε × 181):
0.005 → ~0.9 px (well-corrected apochromat); 0.01 → ~1.8 px (typical achromat); 0.02 → ~3.6 px
(poorly corrected). Sign raised as a real question — which channel magnifies more depends on which
is the longer emission wavelength, and users image in either order.

Options presented:
- **`Uniform(-0.02, 0.02)`, symmetric** — matches the ±3 px translation headroom so both axes are
  comparable; symmetric because channel ordering is the user's choice, not the model's.
- `Uniform(-0.01, 0.01)`, symmetric — attenuates `ρ_true` less; half the SC2 chromatic headroom.
- `Uniform(0, 0.02)`, one-sided — most informative per unit spread, but a silent correctness trap
  when users swap channels.

**Selected: `Uniform(-0.02, 0.02)`, symmetric.** → D-09

### Q3. How should the chromatic scale compose with the existing translation?

Correctness issue raised: every resampling pass applies interpolation smoothing, and **smoothing
decorrelates** — a second `warp` would widen the posterior for a reason unrelated to misalignment,
silently contaminating exactly what SC2 measures.

Options presented:
- **Single composed affine warp** — one `AffineMap = Translation(dy, dx) ∘ LinearMap(scale about
  centre)` in a single `warp` call; one interpolation pass, so the `ε = 0, shift = 0` case is
  byte-comparable with today's simulator and gives a regression check. Cost: the `WR-06` axis-order
  comment must be re-derived.
- Two sequential warps — minimal diff, `WR-06` stays valid verbatim; two interpolation passes
  inflate posterior width for the wrong reason.
- Composed, with centre as a named constant.

**Selected: Single composed affine warp.** → D-10

### Q4. Does the ε extension mirror into `src/amortized/simulator.jl`, or stay spike-only?

Constraint surfaced: `src/amortized/datagen.jl:208-213` lists `simulator.jl` in
`DATAGEN_HASH_SRC_FILES` — its source *bytes* feed the content-hash digest.

Options presented:
- spike/ only — src untouched, digest unchanged, provenance chain intact; two simulators drift apart.
- **Mirror into src/ now** — the two simulators never diverge; changes the datagen digest, orphans
  caches, severs reproducibility of the shipped net's training distribution from the current tree.
- spike/ only, plus a divergence test asserting the digest is unchanged.

**Selected: Mirror into src/ now.** → D-11

> **Correction issued during discussion.** The initial framing overstated the cost. The shipped
> artifact is pinned by `git-tree-sha1` in `Artifacts.toml` (`90e6b63a…`), **not** by the datagen
> digest. Editing `simulator.jl` does not break the artifact download or
> `colocalization_amortized(...)`. It orphans training-data caches and severs the claim that current
> `src/` reproduces the shipped net's training distribution — a narrower, documentation-shaped cost.

### Q5. How should the severed provenance claim be recorded?

Options presented:
- **Named limit + pinned pre-ε commit** — record in `docs/amortized.md` that the shipped `grid_8`
  bundle was trained by `src/amortized/simulator.jl` as of commit `<sha>` (pre-ε), so the training
  distribution stays exactly recoverable via git. Provenance repaired by reference rather than by
  bytes; consistent with how named limits #1-#4 are handled.
- Version the simulator in-file (`SIMULATOR_VERSION`) — strongest reproducibility; dead code in a
  hashed file forever and two paths to maintain.
- Regenerate and re-pin — cleanest end state, but reopens the ship lane and contradicts D-01.

**Selected: Named limit + pinned pre-ε commit.** → D-12

---

## Area 4: Mis-registration handling

**Framing established:** `src/registry.jl:158` — "The persisted `ood_nulls` carries only the
`:density` channel", so the shipped flag is Mahalanobis on the continuous summary rows even though
`ood.jl` implements density + posterior-predictive + noise.

**Tension named before questioning:** the density channel fires when a summary looks unlike the
*training* summaries. Widening the shift prior (D-02) and adding ε (D-09) deliberately pull
mis-registered summaries *into* the training distribution. The wider the prior, the less
mis-registration is detectable as OOD — SC2's mechanism directly erodes SC3's detector. They are not
independent criteria.

### Q1. What does SC3 "handled without silent overconfidence" mean operationally?

Options presented:
- Widen in-domain, flag beyond-prior — partition at the training prior; needs a pre-registered
  boundary magnitude and a check that the density channel actually separates beyond-prior shift.
- **Widening alone satisfies SC3** — purely a no-overconfidence claim: the `Δρ` credible interval
  must still cover the truth. No flag requirement. A covered interval *is* the absence of silent
  overconfidence. Cost: badly mis-registered data returns a wide-but-unflagged answer — honest but
  silent.
- Add a registration OOD channel — breaks the tension outright, but needs its own fitted nulls and
  is not in the shipped artifact.

**Selected: Widening alone satisfies SC3.** → D-13

### Q2. Over what misalignment range must SC3 coverage hold?

Raised: "any injected misalignment" cannot be literal — a 50 px shift breaks coverage on any model,
so where SC3 stops testing is what the manuscript can actually claim.

Options presented:
- **Report the breakdown curve** — test to 0 → 8 px and |ε| → 0.05, report where 90% coverage first
  drops below nominal; gate in-prior, report beyond. Yields a quotable limit and pre-empts the
  reviewer question. Costs evaluation rungs only, no extra training.
- In-prior only — cleanest pass/fail, says nothing about the regime users most need to know about.
- In-prior + 2× extrapolation — arbitrary stopping point.

**Selected: Report the breakdown curve.** → D-14

### Claude's-discretion call (reported, not asked)

**Misalignment for SC3 is simulator-injected, not post-hoc.** Force `shift`/`ε` to the test value and
run the single composed affine warp (D-10). Post-hoc warping of already-simulated images would add
the second interpolation pass D-10 exists to eliminate, and its smoothing would decorrelate
independently of misalignment, contaminating the coverage curve. User offered an override; none
requested. → D-15

### Q3. What does Phase 11 actually deliver, given the research-lane / no-reship posture?

Options presented:
- **Report + docs interpretation note** — results report under `.planning/` plus a short
  interpretation section in `docs/amortized.md` quoting the breakdown point; no artifact or API
  change. Cost: the wording must be scrupulous about which model each number came from.
- Report only — cleanest separation; the finding stays where no user will read it.
- Report + docs + API width warning — most useful to users, but touches the frozen public API and
  contradicts D-01.

**Selected: Report + docs interpretation note.** → D-16

### Q4. Is there real mis-registered microscopy data to validate against?

Options presented:
- Simulation-only — exact ground truth, fully controlled; wholly conditional on the forward model
  being right about how misalignment degrades patch correlation.
- **Simulation + existing manuscript images** — deliberately shift one channel of the real images by
  a known amount and confirm the posterior widens as the simulated sweep predicts. No ground-truth
  ρ, so not a coverage test, but a matching width response is strong evidence the mechanism
  transfers. Needs read-only access and a documented note that the pipelines were not modified.
- Simulation + a dedicated acquisition — strongest evidence; instrument-access-blocked.

**Selected: Simulation + existing manuscript images.** → D-17

Target data identified during the answer: `corpus/` (pre-registered manifest, D-08,
`CORPUS_MASTER_SEED = 0x0000000000c05eed`), fetched read-only via `corpus/fetch.jl` /
`corpus/load.jl`.

---

## Deferred ideas raised during discussion

- Summary extension for registration identifiability (lag / cross-correlation features)
- A dedicated registration OOD channel
- Staged marginalize-then-identify comparison
- Off-axis optical centre for the chromatic warp
- Dedicated bead / registration-target acquisition
- API-level posterior-width advisory
- Productionizing the Phase-11 research net

All carried into `11-CONTEXT.md` `<deferred>`.
