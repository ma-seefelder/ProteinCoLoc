# Phase 15: Calibration Operating Envelope and CI Gate — Research

**Researched:** 2026-08-04
**Domain:** Julia (1.12.6) — SBC/ECE measurement harness extension, misspecification generators, GitHub Actions CI for a Julia package
**Confidence:** HIGH for everything measured on disk in this session; MEDIUM for the compute extrapolation to 32 threads; MEDIUM for CI (verified action versions, unverified against an actual run because no `.github/` exists yet)

---

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01: `src`-side, in `test/gate/` — NOT the spike research lane.**
  This deliberately breaks with Phases 13 and 14, which built in `spike/` (their D-01s). Three
  reasons, all structural: SC3's gate must run against the **shipped package** to be a regression
  gate at all; the envelope is a claim about the **shipped net** (D-10); and the three machines this
  phase needs — `test/gate/misspec.jl` (families), `test/gate/sbc.jl` (ranks/ECE/coverage),
  `test/gate/harness.jl` (seeded draw→simulate→infer) — **already live in `test/gate/`**.
  **Consequence for planning:** new sweep generators extend `misspec.jl`'s existing signature; the
  runner extends `run_gate.jl`'s CLI rather than forking a parallel one.

- **D-02: Six axes, UNION of SC1's named four and the shipped four, each TYPED BY MECHANISM.**
  Axes: `spillover` (new generator), `registration` (new generator), `optics` (= SC1's PSF, exists),
  `background` (= SC1's autofluorescence, exists), `texture` (exists), `noise` (exists).
  Every axis carries a mechanism label: **in-prior extrapolation** (spillover, registration,
  background/autofluorescence) vs **out-of-model misspecification** (optics/PSF, texture, noise).
  **The typing is the finding, not bookkeeping.** An in-prior axis that breaks coverage is a
  *training* failure — the net saw those θ and still miscalibrates. An out-of-model axis that breaks
  is a *scope* limit — expected, honest, and precisely what the OOD flag exists for. Reporting one
  number for both would produce a domain map nobody can act on.
  Rationale for the union over SC1's literal four: dropping `noise` would discard the channel the
  Phase-5 OOD pass rests on (`test/gate/misspec.jl:32-33` — the image-noise channel closes the
  correlation-only summary's blind spot). Dropping spillover/registration would discard the two axes
  a microscopist actually controls. **SC1 is EXTENDED, not amended** — all four named axes are
  swept; two more are added.

- **D-03: Each in-prior axis sweeps BOTH within and beyond prior support, in one ladder, with the
  prior boundary as a MARKED RUNG.**
  Directly answers "does calibration hold where we trained, and where does it stop." A break *inside*
  prior support is a materially different (and more serious) result than a break outside it, and a
  ladder that starts at the prior edge cannot distinguish them.
  **Rationale for not assuming in-prior coverage:** the Phase-7 amended gate FAILED on nuisance
  marginals, so "already established" is not free. Out-of-model axes (optics/texture/noise) have no
  prior support to mark and sweep from their fixed/absent baseline outward.

- **D-04: ECE is the BREAK statistic. Coverage band and KS/χ² p-value are computed and reported at
  every rung but are explicitly NOT gating.**
  The break threshold is **anchored, not chosen**: set at the ECE the shipped `amended_v2/grid_8`
  net attains at its own in-prior operating point, so "break" means *measurably worse than what we
  shipped*.
  **Rationale, recorded so it is not re-litigated.** A significance-based break measures `SBC_M`, not
  the tool — at M = 2000 the χ² rejects on drift far below what matters, which is the documented
  "over-powered-χ² artifact" behind part of the Phase-7 gate failure. And a pre-registered coverage
  band is structurally the criterion Phase 12 concluded was *"a criterion whose pass/fail tracks
  something other than what it was meant to measure"* (STATE.md, 2026-07-31 ruling). ECE is an effect
  size and does not grow teeth with M. Reporting all three keeps the other two auditable rather than
  hidden.

- **D-05: The envelope is defined by Δρ and ρ_true ONLY. The six nuisance marginals are measured and
  reported per rung but do not gate.**
  These two are the SBC-calibrated targets (randomized-rank atom handling) and the quantities a user
  acts on. The nuisance marginal drift (~0.08 SD) is an **already-disclosed named v2.0 limit**
  (`docs/amortized.md`); gating on it would re-fail the tool for a reason already on the record and
  would collapse the envelope to near-zero on every axis, producing no map at all.

- **D-06: SC2 is declared PER AXIS and is THREE-VALUED. All three states are defined and
  pre-registered BEFORE any rung is run.**
  - **PROTECTIVE** — `L_ood ≤ L_break`: the flag fires at or before the rung where coverage breaks.
  - **LATE** — `L_ood > L_break`: coverage breaks while the flag is still quiet. **This is the
    failure SC2 exists to catch.**
  - **SILENT-BUT-SAFE** — coverage never breaks within the swept ladder, so no warning was owed.
  **SC2 passes iff no axis is LATE.**
  **Rationale:** "the flag never fired" has two opposite meanings — working correctly where nothing
  was wrong, versus failing to warn where something was. A two-valued pass/fail collapses them, and
  collapsing them is what would force a fifth after-the-fact amendment. Pre-declaring
  `SILENT-BUT-SAFE` before results exist is the cheapest available protection against that.
  Note the asymmetry that motivates it: spillover, registration and background are **in-prior**, so
  the OOD detector was fit on ID data *containing* them and may correctly never fire on them.

- **D-07: "The flag fires" at a rung means the per-rung FIRE-RATE exceeds the MEASURED
  in-distribution baseline by a pre-registered margin.**
  Not "any image fires": `OOD_ID_QUANTILE = 0.95` (`src/amortized/ood.jl:61`) means ~5% of
  *in-distribution* images fire by construction, so an any-image rule fires at rung 1 on every axis
  and proves nothing. Not a chosen fixed rate (e.g. 50%) either — an undserived bar is exactly the
  Stage-1-ceiling pattern that had to be adjudicated. The ID baseline is measured on this net, on the
  Phase-15 stream, and the margin is frozen in `p15_consts.jl` (D-12).

- **D-08: GitHub Actions, TWO-TIER.** New `.github/workflows/`:
  - **Fast tier** — runs on every push and PR (ubuntu-latest, minutes). The deterministic
    calibration golden check (D-09). This is what makes SC3 literally true.
  - **Slow tier** — release-triggered and manually dispatchable (`workflow_dispatch`). Full-scale
    SBC/coverage at reported scale. This is what makes SC3 *mean* something.
  The slow tier needs an explicit trigger discipline recorded in the workflow file, or it silently
  stops running and the fast tier quietly becomes the whole gate.

- **D-09: The fast gate is a GOLDEN-VALUE regression with tight tolerance.**
  Fixed seed + small M (the `SBC_FIX_M = 8` fixture is the right shape); assert the calibration
  statistic matches a committed golden to tight tolerance.
  **This works here specifically because of counter-based Philox seeding** — pools regenerate
  byte-identically (the property that let Phase 11's lost 54 MB pool be rebuilt at zero iteration
  cost). Near-zero flake, seconds not hours.
  **Re-blessing the golden MUST be an explicit, committed, reasoned act** — same append-never-
  overwrite discipline as the constants files. A golden that can be silently regenerated is a rubber
  stamp, and a flaky gate gets disabled. The re-bless procedure is a required deliverable, not
  documentation garnish.

- **D-10: `artifacts/amended_v2/grid_8` carries the REPORTED envelope; `grid_16` is a single
  SECONDARY contrast.**
  The envelope is a claim about the shipped tool, and v2.0 ships the 8×8 reference grid only, so the
  shipped net must carry the reported number. One contrast at *higher* summary dimension is the
  cheapest evidence the envelope is not a grid-8 artifact — and higher is the informative direction:
  if the envelope widens with dimension, that is a concrete v2.1 lever rather than a curiosity.
  `grid_4` is not swept.

- **D-11: A FRESH `PROD_SEED_P15`, derived by the same salted-Philox redraw rule, with disjointness
  ASSERTED. All draws go through `test/gate/harness.jl` on that stream.**
  Asserted disjoint from: all four `PROD_SEED[G]` (G ∈ 4,8,16,32), **`PROD_SEED_V2[8]`** (spent —
  see the domain-section constraint), both `DEV_SEEDS`, `NPE_MASTER_SEED`, `VAL_MASTER_SEED`,
  `DEFAULT_MASTER_SEED`, and `0`. Copy the redraw-loop construction at
  `test/gate/gate_consts_8_v2.jl:278,301-305`.
  **This is structural immunity to C-04, not an index-range check.** A different Philox *key* gives a
  different stream; a different index range on the same key does not. Rejected alternatives, both
  named as "correct constructions" in STATE.md but wrong here: carving from the production pool via
  `pool_indices` inherits the val-block model-selection contamination *while passing every index
  check*; the `VAL_MASTER_SEED` validation stream drove early stopping for the shipped net and is
  contaminated for exactly this use.
  Also assert `seed_gate_global!` is wired into every Phase-15 arm — `test/gate/harness.jl:50-59`
  warns that `sampleposterior` draws from the **global** RNG and not the passed `rng`, and that
  numbers produced without it are not retroactively reproducible.

- **D-12: A full `p15_consts.jl`, two-tier, mirroring `p12_consts.jl` / `spike/p13/consts.jl`.**
  **Tier 1 (append-only bars, frozen before any rung runs):** the ECE break threshold and its anchor
  derivation (D-04), the OOD fire-rate margin (D-07), the ladder endpoints for all six axes (D-03),
  the CI golden tolerance (D-09), and `P15_ITERATION_ALLOWANCE` (D-13).
  **Tier 2 (sentinel-keyed appends):** every measurement, and every later correction — appended,
  never overwritten.
  **Why this weight is warranted here.** STATE.md records this milestone has concluded *"the bar was
  wrong, not the model"* **four times**, and that *"the count is now high enough that 'we
  pre-registered it' no longer settles an argument alone."* The stated standard each amendment was
  held to — shown **from evidence independent of the machinery under audit** to measure something
  other than what it named — applies to every Phase-15 bar and should be recorded next to it.
  Not chosen: extending `gate_consts_8_v2.jl`, which is bound to one spent run and would entangle
  this phase with a closed amendment.

- **D-13: BOTH degenerate outcomes are pre-declared; ladder endpoints are FROZEN; exactly ONE
  documented iteration allowance.**
  Before any rung runs, write down what each degenerate result *means*:
  - **Empty envelope** (breaks at rung 1 on every axis) — the tool has no usable operating range at
    the swept resolution; a reportable negative.
  - **Unbounded envelope** (never breaks within the ladder) — the envelope exceeds the swept range;
    the honest report is the bound, not a wider search.
  `P15_ITERATION_ALLOWANCE = 1`, **spendable on ladder extension ONLY in the "never breaks" case**,
  and recorded as spent when used. Any other extension is an amendment and must be argued as one.
  **Rationale for the freeze:** "extend until it fails" is searching for the failure you want. The
  Phase-12 ruling refused the symmetric move in as many words — *"that would be tuning a model until
  it passes its own calibration gate — a gate that has then stopped measuring anything."*

- **D-14: A LATE axis is REPORTED AS A NAMED LIMIT and the phase CLOSES on it. No retuning.**
  "The OOD flag does not protect against X" is a real, publishable finding and follows the Phase
  11/12 negative-but-useful precedent directly.
  **Explicitly forbidden in response to a LATE axis:** lowering `OOD_ID_QUANTILE` until the flag
  fires early enough (tuning a detector until it passes its own gate), and building a new detector
  channel for the failing mechanism (new capability, its own phase — see Deferred).

### Claude's Discretion

- Exact functional form of the `spillover` and `registration` misspecification generators, provided
  they match the existing `misspec.jl` signature `(rng, θ; imsize, level, G)` and are monotone in
  `level`.
- Ladder resolution — number of rungs per axis and their spacing (linear vs geometric) — subject to
  the endpoints being frozen per D-13.
- Which calibration statistic instantiates the D-09 golden (rank-histogram digest vs ECE scalar vs
  both), and whether the net-artifact hash is asserted alongside it.
- Report/artifact schema and file naming, following the established `*_report.jld2` + atomic
  `.tmp` → reopen-integrity → `mv(...; force=true)` pattern in `run_gate.jl`.
- Whether the `grid_16` contrast (D-10) sweeps all six axes or a justified subset, provided the
  reduction is stated and the reported `grid_8` envelope is complete.

### Deferred Ideas (OUT OF SCOPE)

- **A new OOD detector channel targeting a LATE axis** — scientifically legitimate (the detector is
  multi-channel by design) but it is new capability, reopens the shipped OOD surface, and would let
  this phase tune its way to a pass. Its own phase, v2.1. Explicitly forbidden as a Phase-15
  response by D-14.
- **Retuning `OOD_ID_QUANTILE`** — would make a LATE axis pass. Refused for the same reason the
  Phase-12 ruling refused to hunt the epoch count that lands coverage in the band.
- **Sweeping `grid_4`, and a full three-grid envelope trend** — D-10 takes one contrast, not a trend.
  Revisit if the `grid_16` contrast shows the envelope is strongly grid-dependent.
- **Per-region / per-tile envelope mapping** — blocked on Phase 12's NO, exactly as Phase 14's D-04
  is. v2.1, alongside SPAT-07.
- **Windows CI runner for the fast tier** — the Windows-tauglich constraint is currently asserted in
  prose rather than executably. Adding `windows-latest` to the fast-tier matrix would test it for
  real, at ~2× CI minutes. Worth doing before v2.0 release if Actions minutes allow.
- **Making the slow tier a scheduled (nightly/weekly) job** rather than release-triggered — better
  drift detection, but consumes minutes continuously and needs a failure-routing policy.
</user_constraints>

---

<phase_requirements>
## Phase Requirements

`.planning/ROADMAP.md` Phase 15 records `**Requirements**: TBD` — **no requirement IDs are mapped to
this phase** [VERIFIED: read `.planning/ROADMAP.md`, Phase-15 detail block]. Per the orchestrator's
instruction, CONTEXT.md's D-01…D-14 are the requirement set. Traceability map:

| ID | Description | Research Support (section below) |
|----|-------------|----------------------------------|
| D-01 | src-side lane, in `test/gate/` | §Architecture Patterns — Integration Surface; every machine's exact signature is documented |
| D-02 | Six axes, typed by mechanism | §Generator Contract; §Prior Support vs Level Units; **§Pitfall 2 (`background` typing conflict)**; **§Pitfall 1 (do not mutate `OOD_FAMILIES`)** |
| D-03 | Ladder spans prior boundary, boundary is a marked rung | §Prior Support vs Level Units gives the exact numeric boundary per axis; §Pattern 2 gives the ladder construction that makes the boundary an exact rung |
| D-04 | ECE is the break statistic, anchored | §The ECE Statistic, Exactly; **§The ECE Null Distribution (measured)**; anchors read from disk |
| D-05 | Envelope defined by Δρ and ρ_true only | §The ECE Statistic — per-column; column indices 1 and 8 |
| D-06 | Three-valued SC2 per axis | §OOD Arm — `gate_ood_roc` returns exactly the needed fields |
| D-07 | Fire-rate vs measured ID baseline | §OOD Arm — `id_fire_rate` and `fire_rate[fam][lvl]` are already returned |
| D-08 | Two-tier GitHub Actions | §CI for This Julia Package (verified action versions, Julia pin, GLMakie/OpenGL, artifact fetch) |
| D-09 | Golden-value fast gate | §The D-09 Golden — feasibility, fragility, and what to assert |
| D-10 | grid_8 reported, grid_16 contrast | **§Pitfall 3 — the grid_16 net is a v1 bundle with `recorded = false` provenance** |
| D-11 | Fresh `PROD_SEED_P15`, asserted disjoint | §Seed Derivation — exact construction + exact assertion set |
| D-12 | Two-tier `p15_consts.jl` | §Pattern 4 — the two-tier template, verbatim structure from `p12_consts.jl` |
| D-13 | Frozen endpoints, both degenerate outcomes pre-declared | §Compute Budget (what a frozen endpoint costs); §Open Question 2 |
| D-14 | LATE axis is a named limit, no retuning | No research needed — nothing in this document proposes a detector change |
</phase_requirements>

---

## Summary

This phase needs **no new dependency, no new library, and no new statistical method**. Every machine
it requires already exists in `test/gate/` and every number it needs to anchor against is recoverable
from disk. The research therefore concentrates on six things the planner cannot guess: the exact
generator contract, the exact numeric prior boundaries in each generator's `level` units, the exact
mechanics of the ECE statistic (including its **noise floor**, which turns out to be decisive), the
compute budget (measured, then extrapolated), the CI surface (which does not exist and has three
concrete obstacles), and the seed construction.

**Three findings change what can be planned.**

1. **The ECE break threshold cannot be anchored at the shipped Δρ value as written.** The shipped
   `amended_v2/grid_8` net's in-prior ECEs are **ρ_true = 0.03713** and **Δρ = 0.00429**
   [VERIFIED: read from `artifacts/amended_v2/grid_8/gate_report_8.jld2` this session]. But the
   *expected* ECE of a **perfectly calibrated** net at M = 2000 is **0.00751** (measured by Monte
   Carlo, 2000 replicates of uniform ranks, this session). Δρ's anchor sits **below its own null
   mean** — it is a lucky draw from the noise distribution, not a signal. A break threshold anchored
   literally at 0.00429 declares a break on ~90 % of *perfectly calibrated* rungs. The break
   threshold must be `ECE(rung 0, same M, same stream) + margin`, where the margin is derived from
   the null distribution of the statistic at the chosen M — a quantity computable in closed form
   before any rung runs, from M and the level grid alone, with no net and no simulation.

2. **The D-10 `grid_16` contrast is not like-for-like and the shipped invariant refuses it.**
   `artifacts/grid_16/npe_16.jld2` reports `training_imsize_provenance = (…, recorded = false)`
   [VERIFIED: `ProteinCoLoc.training_imsize_provenance` called on the bundle this session]. It is a
   **v1** net — trained at 256², unbounded-θ, before the F5 mixture retrain — so `run_gate` returns
   `:provenance_mismatch` and runs no arm, and even bypassing that guard confounds *grid dimension*
   with *image-size regime* and *θ-space parameterisation*. There is no `artifacts/amended_v2/grid_16`.
   This needs a user ruling before planning (see Open Questions).

3. **`MCE` is structurally broken in this pipeline and must not be used.** Every one of the eight
   columns records `mce = 0.99` [VERIFIED, all 8 columns, this session]. `_bin_calibration` walks all
   50 reliability bins including the 31 empty ones, where `predicted_rate = midpoint` and
   `observed_rate = 0.0`, so the maximum gap is pinned at the highest bin midpoint forever. ECE is
   unaffected (empty bins carry weight 0).

**Primary recommendation:** Build the sweep as a thin new arm (`--envelope`) over the *existing*
`sbc_gate` / `gate_ood_roc` with a Phase-15-local `P15_FAMILIES` NamedTuple and a Phase-15-local
seed; run it at **M = 1000 per rung, 5 rungs per axis, rung 0 shared** (≈ 7.3 h wall-clock for
grid_8 at 32 threads); anchor the break threshold on rung 0 at the same M on the same stream, with a
margin derived from the Monte-Carlo null of the ECE statistic; and build CI as two workflows that pin
Julia to **1.12.6**, instantiate the committed `Manifest.toml`, install the OpenGL/xvfb system
libraries GLMakie needs, and load the net through `ProteinCoLoc.estimator_for(8)` (which resolves the
shipped bundle from `Artifacts.toml`, hash-verified, no network needed on the dev path).

---

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Prior draw + forward simulation | `src/amortized/simulator.jl` | — | Shipped, frozen; Phase 15 must not touch it (no retraining/retuning, §6.4) |
| Misspecified forward models | `test/gate/misspec.jl` | — | Already the home of the four families; D-01 puts the two new ones here |
| Frozen summary read chain | `test/gate/harness.jl` `gate_summary` / `draw_simulate_infer` | `src/amortized/summary.jl` | The harness owns the *frozen-stats* discipline; `src` owns the transform maths |
| Posterior draws | `src/amortized/infer.jl` `posterior_for` | NeuralEstimators `sampleposterior` | Takes no `rng`; the **global** RNG must be pinned by the arm (harness.jl:50-59) |
| Rank → ECE / coverage / KS-χ² | `test/gate/sbc.jl` | HypothesisTests.jl | Never hand-rolled; ECE is `_bin_calibration` (ported reliability diagram) |
| OOD scoring + fire rate | `src/amortized/ood.jl` (read surface) | `test/gate/misspec.jl` `gate_ood_roc` (experiment) | The split is stated in `ood.jl`'s header; Phase 15 consumes both, changes neither |
| Pre-registration constants | `test/gate/p15_consts.jl` (**new**) | — | D-12; selected via the existing `run_gate.jl --consts` override |
| Sweep orchestration + report | `test/gate/run_gate.jl` (**new arm**) | — | D-01; atomic `.tmp` → reopen → `mv` write pattern already there |
| CI orchestration | `.github/workflows/` (**new**) | `test/gate/` golden script | The workflow owns triggers/matrix; the science owns the assertion |

---

## Standard Stack

### Core — everything is already a dependency; **add nothing**

| Library | Version (resolved) | Purpose | Why Standard |
|---------|-------------------|---------|--------------|
| `Random123` | in `Manifest.toml` | `Philox4x` counter-based seeding for `PROD_SEED_P15` | Project-wide pattern; C-04's guarantee half rests on it |
| `HypothesisTests` | in `Manifest.toml` | KS + χ² rank uniformity (reported, not gating, per D-04) | `sbc_uniformity` already calls it; never hand-roll |
| `Statistics` / `StatsBase` | stdlib / manifest | ECE arithmetic, `reconstruct` for the frozen θ transform | Already the whole gate's substrate |
| `JLD2` | in `Manifest.toml` | Report persistence, same schema as `gate_report_*.jld2` | Established atomic write pattern in `run_gate.jl` |
| `ImageFiltering` | in `Manifest.toml` | `imfilter` / `Kernel.gaussian` in the misspec generators | Already used by `misspec.jl` and `simulator.jl` |

**Installation:** none. `Project.toml` and `Manifest.toml` must stay byte-unchanged.

**Why zero new dependencies matters here beyond tidiness** [VERIFIED: read `Project.toml`,
`Manifest.toml`, `test/runtests.jl:32-49`]: `Project.toml` declares
`[compat] NeuralEstimators = "0.2.1"`, `Flux = "0.16.10"` — Julia caret semantics, so **0.2.2 and
0.16.11 are both admissible**. But `test/runtests.jl:47-48` asserts the *exact strings*
`"0.2.1"` / `"0.16.10"`. Any change that forces a re-resolve can therefore fail the co-resolution
hard gate on an upstream patch release that nobody in this repo did anything to invite. The committed
`Manifest.toml` (`julia_version = "1.12.6"`, `NeuralEstimators 0.2.1`, `Flux 0.16.10`) is what holds
that together, and every CI job must instantiate it rather than resolve afresh.

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Extending `run_gate.jl` with `--envelope` | A standalone `test/gate/run_envelope.jl` | D-01 says extend, not fork. Extending also inherits `--consts`, `--artifacts-root`, the atomic writer and the `:not_trained` honesty path for free. **Follow D-01.** |
| Adding two families to `OOD_FAMILIES` | A Phase-15-local `P15_FAMILIES = merge(OOD_FAMILIES, (; spillover, registration))` | **Strongly prefer the local merge** — see Pitfall 1. `gate_ood_roc` and `ood_gate` both take `families` as a keyword, so this is zero-touch. |
| `julia-actions/julia-runtest` for the fast tier | `julia --project=. <golden script>` after explicit `Pkg.instantiate()` | Prefer the explicit script: `Pkg.test()` sandboxes and re-resolves; `Pkg.instantiate()` on the committed manifest does not. Also avoids running the 1550-line full suite on every push. |
| Bit-exact golden equality | Tolerance-based golden with a pinned Julia patch version | Prefer tolerance + pin — see Pitfall 5 (Julia's global RNG stream is not a cross-version contract). |

---

## Package Legitimacy Audit

**Not applicable — this phase installs no external packages.** Every library it uses is already a
declared dependency in the committed `Project.toml`/`Manifest.toml`, resolved and in use by shipped
code. No `npm`/`pip`/`cargo`/`Pkg.add` step appears anywhere in the recommended plan, and the hard
constraint "prefer zero new dependencies" is satisfiable in full.

The only registry interaction in the plan is CI's `Pkg.instantiate()` against the **committed
manifest**, which installs exactly the pinned `git-tree-sha1`s and resolves nothing.

---

## Architecture Patterns

### System Architecture Diagram

```
                    p15_consts.jl  (Tier 1 bars, frozen; Tier 2 appends)
                            |  --consts
                            v
   PROD_SEED_P15 ---> prod_seed/prod_rng aliases ---> seed_gate_global!(G)
        |                        |                            |
        |                        v                            v
        |          Philox stream per (axis, rung)     GLOBAL RNG (flow base draws)
        |                        |                            |
        v                        v                            |
  sample_prior(π) ------> [ MISSPEC GENERATOR ] --------------+
   (unchanged)             gen(rng, θ; imsize, level, G)
                                 |
                                 |  2 x Matrix{Float64}
                                 v
                          build_mci -> patch_summary(·, G)
                                 |
                                 v
                     encode_d01 -> standardize_summary(·, m.zt, :min)     [FROZEN]
                                 |
                     +-----------+------------------------+
                     |                                    |
                     v                                    v
        posterior_for (CPU) -> reconstruct(m.θzt)    maha_score(density,·)
                     |                                noise_score(nnull,pair)
                     v                                    |
             M x 8 RANK TABLE                             v
                     |                            zfuse -> > fused_thr ?
        +------------+-----------+                        |
        v            v           v                        v
   sbc_calibration  sbc_uniformity  sbc_coverage    fire_rate[axis][rung]
     -> ECE          -> KS/chi2 p    -> curve        vs id_fire_rate
        |  (GATING, cols 1 & 8)  (reported)               |
        v                                                 v
   L_break(axis) = first rung with ECE > thr      L_ood(axis) = first rung with
        |                                          fire_rate > id + margin
        +---------------------+---------------------------+
                              v
                 SC2 verdict per axis: PROTECTIVE / LATE / SILENT-BUT-SAFE
                              |
                              v
              envelope_report_8.jld2   (atomic .tmp -> reopen -> mv)
                              |
                              v
                   [ CI FAST TIER ]  golden(ECE, rank digest, net hash)
                   [ CI SLOW TIER ]  full-scale re-run on release/dispatch
```

### Recommended File Layout

```
test/gate/
├── p15_consts.jl          # NEW — two-tier pre-registration (D-12), selected by --consts
├── p15_misspec.jl         # NEW — misspec_spillover / misspec_registration + P15_FAMILIES
├── misspec.jl             # UNCHANGED (OOD_FAMILIES not mutated — see Pitfall 1)
├── sbc.jl                 # UNCHANGED (consumed as-is)
├── harness.jl             # UNCHANGED (consumed as-is)
├── run_gate.jl            # +1 arm: --envelope, + envelope_gate(...)
└── ci_golden.jl           # NEW — the D-09 fast-tier script (also the golden re-bless entry point)
.github/workflows/
├── ci.yml                 # NEW — fast tier (push + pull_request)
└── calibration-slow.yml   # NEW — slow tier (release + workflow_dispatch)
artifacts/amended_v2/grid_8/
└── envelope_report_8.jld2 # NEW output (git-ignored, like every other artifact)
```

### Pattern 1: The misspecification-generator contract (VERIFIED, exact)

Every family in `OOD_FAMILIES` is a function with this signature
[VERIFIED: `test/gate/misspec.jl:104-105, 125-126, 150-151, 166-167`]:

```julia
gen(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS, G::Integer = 8)
    -> Vector{Matrix{Float64}}     # exactly 2 elements, each `imsize`, all finite and >= 0
```

- `rng` is positional and **untyped** — a `Philox4x` is passed in practice.
- `θ` is positional and **untyped** — the 8-field NamedTuple from `sample_prior`.
- All three of `imsize`, `level`, `G` are **keywords** and all three are always supplied by the
  consumer [VERIFIED: `misspec.jl:371` — `img = gen(rng, θ; imsize = isz, level = lvl, G = G)`].
  The declared defaults are therefore decorative at the call site but must be present for
  signature-compatibility.
- Every family begins with `_guard_misspec_imsize(imsize, G)` (`misspec.jl:76-83`) — throws unless
  both dims ≥ 64 **and** `imsize[d] ÷ G ≥ 4`.
- The return is `[ch1, ch2]`, `Matrix{Float64}`, clamped with `max.(·, 0.0)`.

**What `level` scales, per shipped family** [VERIFIED: read line by line]:

| Family | `level` drives | Value at `level = 1 … 4` |
|--------|----------------|--------------------------|
| `texture` | spot count `30·level` (shared), `20·level` (private ×2); smoothing `σ = max(0.8, 3.0 − 0.5·level)` | 30/60/90/120 spots; σ = 2.5/2.0/1.5/1.0 (**saturates at 0.8 for level ≥ 4.4**) |
| `noise` | spike fraction `0.02·level`; amplitude `2.0·level` | 2/4/6/8 % pixels; amp 2/4/6/8 × max |
| `optics` | `σy = σ_psf + 2.0·level`, `σx = σ_psf + 0.2·level` | σy = 3.3/5.3/7.3/9.3 px against a fixed `σ_psf = 1.3` |
| `background` | gradient `0.5·level`/`0.3·level`, vignette `0.2·level`, **bleed `0.1 + 0.3·level`** | bleed = 0.4/0.7/1.0/1.3 |

**How `G` is threaded:** only into `_guard_misspec_imsize` inside the generators; the *summary* read
uses `G` separately in `gate_summary(m, img, G)` (`misspec.jl:64-66`) and in
`negctrl_block_permute(...; blocks = G)`. A generator that ignores `G` beyond the guard is contract-
compliant (three of the four do exactly that).

**How `OOD_FAMILIES` is consumed** [VERIFIED: `misspec.jl:316-396`]: `gate_ood_roc` takes
`families = OOD_FAMILIES` as a **keyword**, then uses only `keys(families)` and `families[fam]`. It
loops `for lvl in 1:levels` with `levels::Integer = OOD_GRID_LEVELS`, also a keyword. `ood_gate`
(`run_gate.jl:290-295`) likewise forwards `families` and `levels` as keywords. **A new family
therefore slots in with zero consumer change — by passing a different NamedTuple, not by editing the
frozen one.**

### Pattern 2: The two new generators (recommended concrete form)

The design goal is that the *prior boundary is an exact rung* (D-03) and that rung 0 is a
byte-identity with the unmisspecified simulator (an executable correctness guard, the project's
established style — cf. `spike/test/test_p12_decoupling.jl`).

```julia
# test/gate/p15_misspec.jl
"""
    misspec_spillover(rng, θ; imsize = SBC_IMSIZE, level = 0, G = 8)

Axis: SPILLOVER, mechanism = IN-PRIOR EXTRAPOLATION.
`level = 0` is the pass-through in-prior rung (byte-identical to `simulate_pair`); levels
1..length(P15_SPILLOVER_LADDER) OVERRIDE θ.spillover with the frozen ladder magnitude.
`P15_SPILLOVER_LADDER[2] == maximum(SPILLOVER_PRIOR) == 0.2` — the MARKED PRIOR BOUNDARY rung.
"""
function misspec_spillover(rng, θ; imsize = SBC_IMSIZE, level::Integer = 0, G::Integer = 8)
    _guard_misspec_imsize(imsize, G)
    level == 0 && return ProteinCoLoc.simulate_pair(rng, θ; imsize = imsize)
    s = P15_SPILLOVER_LADDER[level]
    return ProteinCoLoc.simulate_pair(rng, merge(θ, (; spillover = s)); imsize = imsize)
end

function misspec_registration(rng, θ; imsize = SBC_IMSIZE, level::Integer = 0, G::Integer = 8)
    _guard_misspec_imsize(imsize, G)
    level == 0 && return ProteinCoLoc.simulate_pair(rng, θ; imsize = imsize)
    r = P15_SHIFT_LADDER[level]                    # px, |(dx, dy)|
    d = r / sqrt(2.0)                              # fixed 45° direction: |shift| == r exactly
    return ProteinCoLoc.simulate_pair(rng, merge(θ, (; shift_dx = d, shift_dy = d));
                                      imsize = imsize)
end

const P15_FAMILIES = merge(OOD_FAMILIES, (spillover    = misspec_spillover,
                                          registration = misspec_registration))
```

Four properties this buys, each of them load-bearing:

1. **Exact prior-boundary rung.** `P15_SPILLOVER_LADDER[2] = 0.2` is *literally*
   `maximum(SPILLOVER_PRIOR)`; assert that in `p15_consts.jl` rather than writing the number twice.
2. **Level-independent RNG consumption.** `simulate_pair` draws the same number of values regardless
   of `spillover`, `shift_dx`, `shift_dy` [VERIFIED: stages 4 and 6 of `simulate_pair` consume no
   `rng`]. So with the same stream every rung on these two axes sees **the same θ\*, the same image
   sizes and the same latent fields** — a matched-pairs ladder in which the only varying quantity is
   the misspecification magnitude. Rung-to-rung ECE differences are far less noisy than independent
   runs would be. (This does **not** hold for `texture`/`noise`, whose generators draw extra values.)
3. **Rung-0 identity is assertable.** `misspec_spillover(rng, θ; level = 0)` must equal
   `simulate_pair(rng′, θ)` bit-for-bit for identically-seeded `rng`/`rng′`. That is a cheap unit
   test that catches a whole class of wiring errors before any rung burns compute.
4. **Monotone in `level`** by construction (the ladders are increasing tuples), satisfying the
   discretion constraint.

**Guards the plan must respect:** `simulate_pair` validates `0.0 ≤ spillover ≤ 1.0` (`simulator.jl:232-233`)
— so `P15_SPILLOVER_LADDER` must not exceed 1.0. There is **no** range guard on `shift_dx`/`shift_dy`
(only finiteness), so the shift ladder is unbounded from the simulator's side; large shifts fill the
vacated border with `BG_FLOOR` (`simulator.jl:290`), which is itself part of the misspecification.

### Pattern 3: The sweep arm — reuse, do not rebuild

`misspec_simulator` (`misspec.jl:194-199`) already wraps a generator into the exact
`(sample_prior, simulate_pair, build_mci)` contract that `sbc_gate` consumes as `sim`. So **an SBC
run at a misspecification rung is one existing call**:

```julia
sbc_gate(m; G = 8, M = P15_SWEEP_M, L = SBC_L,
         imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET, imsize_weights = GATE_IMSIZE_WEIGHTS,
         sim = misspec_simulator(P15_FAMILIES[axis]; level = rung, G = 8),
         rng = p15_rng(8, axis, rung))
```

and the OOD fire-rate arm at all rungs of all axes is one existing call:

```julia
gate_ood_roc(m, ood_nulls.density; G = 8, rng = p15_rng_ood(8),
             families = P15_FAMILIES, levels = P15_RUNGS,
             n_id = …, n_pos = …, imsize_set = …, imsize_weights = …)
```

The new code in `run_gate.jl` is therefore an orchestration loop, a break/OOD-crossing reduction, and
a report writer — not a new harness. This is the single largest reason D-01 chose the src-side lane
and it holds up on inspection.

### Pattern 4: The two-tier constants file (structure to copy verbatim)

From `spike/validation/p12_consts.jl:21-95` [VERIFIED, read]:

```julia
# Guarded as ONE Tier-1 block keyed on :P15_… so a re-include is a silent no-op.
if !isdefined(@__MODULE__, :P15_SWEEP_M)
    import Random123: Philox4x
    # --- 1. FORBIDDEN seeds --------------------------------------------------
    # --- 2. FRESH stream (distinct master ⊻ salt, redraw loop) ---------------
    # --- 3. TIER 1 BARS (append-only; each carries its derivation inline) ----
    # --- 4. GATING vs REPORTING-ONLY constant-name tuples, intersection asserted empty
    # --- 5. Self-checks (cheap; no inference, no simulation) -----------------
end
# --- TIER 2 blocks, each keyed on its OWN sentinel const, appended never edited ---
if !isdefined(@__MODULE__, :P15_ECE_ANCHOR_MEASURED)  # opened by plan 15-xx
    …
end
```

Two details that are easy to get wrong and are stated explicitly in the source:
- **Tier 2 is a SEQUENCE of blocks, each with its OWN sentinel** — a second append reusing an earlier
  sentinel is silently skipped once that block exists (`p12_consts.jl:33-37`).
- **`P15_GATING_CONSTANTS` / `P15_REPORTING_ONLY_CONSTANTS` as two disjoint name tuples with an
  in-file `@assert isempty(intersect(...))`** makes D-04/D-05's "reported is not gated" split
  machine-checkable rather than prose (`p12_consts.jl:52-57`).

### Anti-Patterns to Avoid

- **Mutating `OOD_FAMILIES`** — see Pitfall 1.
- **Reusing `PROD_SEED_V2[8]`** — spent. `run_gate.jl:47-54` records the one-run binding, and the
  report it produced is on disk with `seed = 1135605683775656488`.
- **Selecting `gate_consts_8_v2.jl` for a Phase-15 run** — it is bound to a closed amendment; D-12
  mirrors it into a new file for exactly this reason. Use `--consts p15_consts.jl`.
- **Treating `mce` as a statistic** — it is 0.99 for every column (Pitfall 4).
- **Gating on KS or χ² p-values** — D-04 forbids it, and the M=2000 over-power artifact is on the
  record; but *do* compute and report them (they are free — one more reduction of one rank table).

---

## Prior Support vs `level` Units — the D-02/D-03 typing, made numeric

[All VERIFIED from `src/amortized/simulator.jl:85-98, 134-146, 166-175, 179-183`]

| θ field | Prior | Support | Fixed? |
|---------|-------|---------|--------|
| `ρ_true` | `ghat(μ*)`, `μ* ~ Truncated(Cauchy(0, 0.3), -1, 1)` | `[-0.99, +0.99]` (extrema of `GHAT_RHO_KNOTS`) | no |
| `spillover` | `Uniform(0.0, 0.2)` | **`[0, 0.2]`** | no |
| `autofluorescence` | `Uniform(0.0, 0.1)` | **`[0, 0.1]`** | no |
| `label_efficiency` | `Uniform(0.6, 1.0)` | `[0.6, 1.0]` | no |
| `shift_dx`, `shift_dy` | `Uniform(-1.0, 1.0)` | **`[-1, 1]` px each** | no |
| `noise` | `Uniform(0.0, 1.0)` | `[0, 1]` | no |
| `chromatic_eps` | `Uniform(-0.02, 0.02)` | `[-0.02, 0.02]` | no |
| **`σ_psf`** | — | **none — `const σ_psf = 1.3`** (`simulator.jl:179`) | **YES** |

**Per-axis prior boundary, in the units each generator's `level` operates in:**

| Axis | Mechanism (D-02) | Quantity the ladder moves | Prior boundary, in that quantity | Boundary expressible as a rung? |
|------|------------------|---------------------------|----------------------------------|---------------------------------|
| `spillover` | in-prior extrapolation | `θ.spillover` (override) | **0.2** | **YES** — make it `P15_SPILLOVER_LADDER[2]` |
| `registration` | in-prior extrapolation | `|(shift_dx, shift_dy)|` px (override) | **1.0 px** | **YES** — make it `P15_SHIFT_LADDER[2]` |
| `background` | in-prior extrapolation *(per D-02)* | bleed offset `0.1 + 0.3·level` | **0.1** (`maximum(AUTOFLUORESCENCE_PRIOR)`) | **NO — see Pitfall 2.** `level = 1` already gives 0.4 = 4× the boundary |
| `optics` | out-of-model | `σy = 1.3 + 2.0·level` | **none** — `σ_psf` is fixed, not a prior | N/A; sweep from the fixed baseline outward |
| `texture` | out-of-model | spot count / sharpness | **none** — a generator the simulator cannot express at any θ | N/A |
| `noise` | out-of-model | spike fraction + heavy tail | **none** — beyond the Poisson+Gaussian model at any `θ.noise` | N/A |

Note for the report's honesty section: `misspec_optics` also perturbs `σx = 1.3 + 0.2·level`, so the
axis moves *both* anisotropy and overall blur; it is "aberrated/out-of-focus", not pure anisotropy.

---

## The ECE Statistic, Exactly

**Which function yields ECE:** `sbc_calibration(ranks; levels, L, n_bins) -> CalibrationResult`
(`sbc.jl:302-315`), whose `.ece` field `sbc_gate` records per column (`sbc.jl:451-456`).

**Exact signatures and return shapes** [VERIFIED, read `test/gate/sbc.jl` in full]:

| Function | Signature | Returns |
|----------|-----------|---------|
| `sbc_ranks` | `(m; G, M, L, imsize, imsize_set, imsize_weights, sim, rng)` | `Matrix{Int}` `M×8`, each entry ∈ `0:L` |
| `sbc_ranks_and_spread` | same | `(ranks::M×8 Int, post_sd::M×8 Float64, prior_draws::M×8 Float64, imsizes::Vector{Tuple{Int,Int}} (len M), imsizes_delta::same)` |
| `sbc_uniformity` | `(ranks::AbstractVector{<:Integer}; L, bins)` | `(ks_p, chi2_p)` — **one column at a time** |
| `sbc_coverage` | `(ranks::Vector{<:Integer}; levels = 0.05:0.05:0.95, L)` | `(nominal::Vector, empirical::Vector)`, length 19 |
| `sbc_calibration` | `(ranks::Vector{<:Integer}; levels, L, n_bins)` | `CalibrationResult(bin_midpoints, predicted_rate, observed_rate, bin_counts, ece, mce)` — 50-element vectors |
| `sbc_traffic_light` | `(ece::Real)` | `:green` / `:yellow` / `:red` against `SBC_ECE_GREEN`/`_YELLOW` |
| `sbc_gate` | `(m; G, M, L, bins, chi2_bins, imsize…, sim, rng)` | large NamedTuple; `per_param` is an 8-element `Vector{NamedTuple}` with `(label, ks_p, chi2_p, ece, mce, light, shrinkage, post_sd, prior_sd, vacuous)` |

**Column order** (`SBC_PARAM_LABELS`, `sbc.jl:123-124`):
`1 = ρ_true`, 2 = spillover, 3 = autofluorescence, 4 = label_efficiency, 5 = shift_dx, 6 = shift_dy,
7 = noise, **`8 = Δρ`**. **D-05's envelope is therefore columns 1 and 8.**

### What the ECE number actually is — derived and then verified numerically

`sbc_calibration` forms, for each of the 19 nominal levels α ∈ {0.05, 0.10, …, 0.95} and each of the
M ranks, the pair `(predicted = α, positive = [|u − 0.5| ≤ α/2])` with `u = (rank+0.5)/(L+1)`. It
then hands 19·M pairs to `_bin_calibration` with 50 equal-width bins over [0,1]. Bin width is 0.02
and the α values are spaced 0.05, so **each α lands in its own bin and 31 of the 50 bins are empty**.
Empty bins get `bin_counts = 0`, hence `weight = 0` in the ECE sum. Each occupied bin holds exactly M
samples, `predicted_rate = α` exactly and `observed_rate = ĉ(α)` = empirical central-interval
coverage. Therefore:

> **ECE = (1/19) · Σ_{α ∈ {0.05,…,0.95}} |α − ĉ(α)|** — a plain mean absolute coverage error over
> the 19 nominal levels. It does not depend on `SBC_BINS` at all (as long as no two α share a bin).

**Verified numerically against the shipped report this session** — recomputing the identity by hand
from the stored rank table reproduces the stored value to the last float digit:

| column | stored `ece` | hand-computed `mean|α − ĉ(α)|` |
|--------|--------------|-------------------------------|
| 1 (ρ_true) | `0.03713157894736839` | `0.037131578947368425` |
| 8 (Δρ) | `0.0042894736842104715` | `0.004289473684210523` |

### The in-prior anchors, read from the shipped report

`artifacts/amended_v2/grid_8/gate_report_8.jld2` — **the ECE anchor IS recoverable, nothing needs
re-measuring to obtain it** [VERIFIED, opened this session].

- JLD2 top-level keys: `["schema_version", "report"]`.
- `report` fields: `(:status, :grid, :seed, :gate_consts_version, :imsize_set, :imsize_weights,
  :realised_imsize_counts, :training_imsize_provenance, :sbc, :bf, :ood)`.
- `report.sbc` fields: `(:grid, :M, :L, :ranks, :per_param, :post_sd, :prior_draws, :bins,
  :chi2_bins, :imsize_set, :imsize_weights, :realised_imsize_counts, :imsizes, :imsizes_delta,
  :vacuous_cutoff, :vacuous_params, :rules_version, :v1_verdict, :v2_verdict, :ks_pass, :chi2_pass,
  :ece_pass, :passed, :caption)`.
- `report.sbc.ranks` is the **full 2000×8 Int rank table** — so any rank-derived statistic can be
  recomputed post-hoc without re-running anything.
- `report.seed = 1135605683775656488` (= `PROD_SEED_V2[8]`, the spent one).
- `report.training_imsize_provenance = (imsize_set = ((512,512),(1024,1024),(1376,1028),(2048,2048)),
  imsize_weights = (0.4,0.25,0.25,0.1), imsize_source = :generate_samples, recorded = true)`.

Per-column values at M = 2000, L = 999:

| column | ece | mce | ks_p | shrinkage | light |
|--------|-----|-----|------|-----------|-------|
| **ρ_true** | **0.03713** | 0.99 | 1.12e-7 | 0.1302 | green |
| spillover | 0.01921 | 0.99 | 6.82e-6 | 1.0101 | green |
| autofluorescence | 0.00321 | 0.99 | 0.0761 | 0.9974 | green |
| label_efficiency | 0.03300 | 0.99 | 1.92e-4 | 0.9226 | green |
| shift_dx | 0.00284 | 0.99 | 0.8834 | 0.9759 | green |
| shift_dy | 0.01532 | 0.99 | 0.0017 | 0.9869 | green |
| noise | 0.01724 | 0.99 | 0.0192 | 1.0250 | green |
| **Δρ** | **0.00429** | 0.99 | 0.1964 | 0.1363 | green |

### The ECE Null Distribution — **measured**, and it is decisive

Because ECE is a mean of 19 absolute deviations of binomial proportions with n = M each, it has a
**positive expectation even for a perfectly calibrated estimator**. This is computable with no net,
no simulation and no data: draw M uniform ranks, compute the statistic, repeat. Measured this session
(2000 replicates per M, uniform ranks on `0:999`, `L = 999`, `n_bins = 50` — i.e. `_bin_calibration`
reproduced exactly):

| M | E[ECE₀] | sd | q95 | q99 |
|---|---------|-----|-----|-----|
| 100 | 0.03248 | 0.01457 | 0.06105 | 0.07684 |
| 250 | 0.02059 | 0.00896 | 0.03789 | 0.04884 |
| 500 | 0.01431 | 0.00651 | 0.02737 | 0.03484 |
| **1000** | **0.01038** | 0.00449 | **0.01884** | 0.02463 |
| 2000 | 0.00751 | 0.00338 | 0.01424 | 0.01821 |

Closed form that reproduces the table to ~3 %: `E[ECE₀] ≈ 0.336 / √M`
(derivation: `E|ĉ−α| ≈ √(2/(πM))·√(α(1−α))`, averaged over the 19 α giving mean `√(α(1−α)) = 0.4085`;
`0.4085 · √(2/π) = 0.326`). Put this derivation next to the bar in `p15_consts.jl` per D-12.

**Three consequences the planner must design around.**

1. **The Δρ anchor (0.00429) is below its own null mean (0.00751).** Anchoring the break threshold at
   the literal shipped Δρ value would make a *perfectly calibrated* rung break with probability
   ≈ 0.9. That is the D-04 failure mode ("a criterion whose pass/fail tracks something other than
   what it was meant to measure") reappearing in a new place, and it is catchable before any rung
   runs — which is exactly the licensing standard STATE.md sets for a bar.
2. **The fix is a rung-0 anchor at the sweep's own M on the P15 stream, plus a null-derived margin.**
   Recommended, and derivable before any result:
   `P15_ECE_BREAK_THRESHOLD(col) = ECE_rung0(col) + P15_ECE_MARGIN`, with
   `P15_ECE_MARGIN = q95(ECE₀ at P15_SWEEP_M) − E[ECE₀ at P15_SWEEP_M]`
   (at M = 1000: `0.01884 − 0.01038 = 0.00846`). This is an *effect size* margin — it does not grow
   teeth with M, and it is calibrated so that a genuinely-still-calibrated rung crosses it ~5 % of
   the time. Both terms come from arithmetic on the test design, never from an outcome.
   Rung 0 must be measured once per net at the sweep's M and appended as a Tier-2 constant.
3. **Reducing M inflates every rung's ECE, including rung 0** — but because the anchor is measured at
   the same M, the inflation is *common-mode and cancels*. What M actually costs is **resolution**:
   the smallest true ECE degradation the ladder can detect is ≈ the margin. At M = 1000 that is
   ≈ 0.0085; at M = 500, ≈ 0.0131; at M = 250, ≈ 0.0173; at M = 2000, ≈ 0.0067. Since the shipped
   ρ_true ECE is 0.0371, a resolution of 0.0085 detects a ~23 % relative degradation — adequate. A
   resolution of 0.0173 (M = 250) detects only a ~47 % degradation, which would blur the ladder.

---

## Compute Budget

**This is the planning blocker the phase description flags, and it resolves favourably.**

### Measured on this machine, this session (Julia 1.12.6, **1 thread**, CPU)

| Image size | `simulate_pair` (s) | `patch_summary` G=8 (s) | `patch_summary` G=16 (s) |
|------------|--------------------|-------------------------|--------------------------|
| 512×512 | 0.332 | (0.014 after JIT) | 0.014 |
| 1024×1024 | 2.448 | 0.086 | 0.119 |
| 1376×1028 | 2.066 | 0.177 | 1.622 † |
| 2048×2048 | 5.928 | 0.396 | 1.875 † |

`posterior_for(estimator, Z; N = 999, use_gpu = false)` = **0.028 s**, independent of image size.
`load_estimator` (cold, per process) = 21 s. `using ProteinCoLoc` (cold) = 91 s including
precompilation of the package itself.

† Single un-replicated samples with visible GC/JIT contamination (the 1024² sim time is also high
relative to a linear-in-pixels fit). **The G=16 summary numbers in particular should be re-measured
with `BenchmarkTools` before the grid_16 budget is committed.** A linear-in-pixels fit to the G=8
column gives ≈ 0.094 µs/px; to the simulate column, ≈ 1.41 µs/px at 1 thread.

### Cross-check against the one reported gate run — the model reconciles

The amended confirmatory run (`gate-8x8-amended.md`, `confirmatory-run.log`) took **30.9 min**
wall-clock (17:11:02 → 17:41:59) for SBC(M=2000, L=999) + BF(n=100) + OOD, and
`07-GATE-AMENDMENT.md:140-149` records a **measured** simulate cost table "at 32 threads" giving
`E[cost] = 0.08346 s/pair` over the same mixture.

Per SBC draw the harness does **3 simulate_pair, 3 patch_summary, 3 posterior calls** — one
`draw_simulate_infer` (1 pair) plus one `draw_simulate_infer_paired` (2 pairs, one shared imsize)
[VERIFIED: `sbc.jl:187-207`]. With `E[patch_summary G=8] ≈ 0.107 s` (mixture-weighted from the table
above; not thread-parallel) and `posterior_for = 0.028 s`:

```
t_draw(32 thr) ≈ 3·0.0835 + 3·0.107 + 3·0.028 = 0.656 s
M = 2000  ->  1311 s = 21.9 min  (SBC arm)
   + BF (n=100 -> 200 pairs + KDE)      ≈ 1.5 min
   + OOD (~1130 pairs, no posterior)    ≈ 4   min
   + process + net load                 ≈ 2   min
   ----------------------------------------------
                                        ≈ 29.4 min  vs 30.9 min observed  ✓
```

The model reproduces the only wall-clock datum this project has, from independent measurements. Use
it. **The 1-thread figure for the same arm is ≈ 2.9 h** — so `JULIA_NUM_THREADS` is not an
optimisation here, it is a 13× budget factor and must be pinned in the runner and in the plan.

### Per-rung and full-grid cost

Misspecification generators cost more than `simulate_pair` (extra `imfilter` at large σ for `optics`,
two full-image `rand` arrays for `noise`); `spillover`/`registration` cost exactly the same, and
`texture` replaces `simulate_pair` entirely. Take an axis-averaged overhead factor **f = 1.3**, so
`t_rung(M) ≈ M · 0.85 s` at 32 threads:

| M | per rung | rung 0 (f = 1.0) | grid_8 total: 1 rung-0 + 6 axes × 5 rungs |
|---|----------|------------------|--------------------------------------------|
| 2000 | 28.4 min | 21.9 min | **14.6 h** |
| **1000** | **14.2 min** | 10.9 min | **7.3 h** |
| 500 | 7.1 min | 5.5 min | 3.6 h |
| 250 | 3.5 min | 2.7 min | 1.8 h |

The OOD fire-rate arm is essentially free: `gate_ood_roc` runs **no posterior draws** (only
`maha_score` / `noise_score`), so 6 families × 5 levels × `n_pos = 40` + `n_fit = n_id = 200` ≈ 1600
pairs ≈ **5.6 min per net** for the whole D-06/D-07 experiment.

### Recommendation — concrete

> **`P15_SWEEP_M = 1000`, `P15_RUNGS = 5` per axis, one shared rung-0 anchor per net,
> `SBC_L = 999` unchanged, run with `JULIA_NUM_THREADS ≥ 16`.**
> grid_8 reported envelope ≈ **7.3 h**; add ≈ 6 min for the OOD arm.

**Justification, stated so it survives review.** M = 2000 is *affordable* (14.6 h) but buys only a
1.27× improvement in break resolution (0.0067 vs 0.0085) for a 2× compute cost, and the reported
scale of `SBC_M = 2000` is not a scientific requirement of *this* experiment — it was frozen in
`gate_consts_8_v2.jl:68` explicitly because "changing M would change power in a direction the author
knew he needed", a concern about a *significance* gate. D-04 has deliberately moved the break
statistic off significance and onto an effect size, which removes exactly that objection: ECE's
sensitivity to M is a *precision* property, disclosed in the table above, not a power-hacking lever.
M = 1000 gives a break resolution of 0.0085 ECE against a shipped ρ_true ECE of 0.0371 — it resolves
a ~23 % relative degradation, which is the granularity a "domain of applicability" map is read at.
Going below M = 500 is not recommended: the resolution (0.013+) starts to approach the shipped
anchor's own magnitude and the ladder stops discriminating.

**If the wall-clock must be halved:** cut *rungs*, not M (4 rungs × 6 axes = 5.8 h at M = 1000).
Cutting M degrades every rung's resolution; cutting rungs degrades only the ladder's granularity,
which D-13 already treats as a reportable bound ("the honest report is the bound, not a wider
search").

---

## The OOD Arm (D-06 / D-07)

`gate_ood_roc` already returns **exactly** the two quantities D-07 names, with no new code
[VERIFIED: `misspec.jl:387-395`]:

- `fire_rate :: Dict{Symbol, Vector{Float64}}` — `fire_rate[family][level]` = fraction of the
  `n_pos` positives at that rung whose OR-fused robust-z exceeds `fused_thr`.
- `id_fire_rate :: Float64` = `mean(id_fused .> fused_thr)` — the **measured** in-distribution
  baseline on this net, on this stream, which is what D-07 requires (it is ≈ `1 − OOD_ID_QUANTILE` =
  0.05 by construction, but measured, not assumed).

So `L_ood(axis) = first rung r with fire_rate[axis][r] > id_fire_rate + P15_OOD_FIRE_MARGIN`.

Sizing note for the margin's derivation: `fire_rate` at a rung is a binomial proportion with
`n = n_pos`. At the shipped `n_pos = 40` its standard error at p ≈ 0.05 is 0.034 — so a margin below
~0.10 is inside one SE of the baseline and is not a detection. Either raise `n_pos` (cheap: the OOD
arm has no posterior draws — `n_pos = 200` for 6 families × 5 levels is ≈ 6000 extra pairs ≈ 20 min)
or derive the margin from the binomial null at the chosen `n_pos`. **Recommend `n_pos = 200` and a
margin derived as the one-sided binomial 95th percentile of `Binomial(n_pos, id_fire_rate)/n_pos`
minus `id_fire_rate`** — again pure design arithmetic, no outcome.

The negative controls (`gate_negctrls`, `verify_summary_invariance`) are **not** needed for Phase 15
and re-running them would re-open a measurement recorded as failing (`neg_pass`: block 0.10, rotate
0.133) in a phase that has no mandate to touch it. Leave them out of the envelope arm.

---

## CI for This Julia Package (D-08 / D-09)

**Starting state: `.github/` does not exist** [VERIFIED: `ls .github` → absent]. The whole directory
is new.

### Verified current action versions (checked this session, not from memory)

| Action | Current major tag | Notes |
|--------|-------------------|-------|
| `julia-actions/setup-julia` | **`@v3`** | inputs `version` (default `'1'`), `arch`, `show-versioninfo`, `include-all-prereleases`; named versions `lts`, `pre`, `nightly`, `min`, `min-minor`, `min-patch` [CITED: github.com/julia-actions/setup-julia] |
| `julia-actions/cache` | **`@v3`** | inputs `cache-name`, `cache-packages`, `cache-artifacts`, `cache-registries`, `cache-compiled`, `cache-scratchspaces`, `delete-old-caches`, `include-matrix`, `depot`, `save-always`, `token` — all default `true` for the five depot subdirs [CITED: github.com/julia-actions/cache] |
| `julia-actions/julia-buildpkg` | `@v1` | inputs `localregistry`, `git_cli`; sets `JULIA_PKG_SERVER_REGISTRY_PREFERENCE=eager` [CITED: github.com/julia-actions/julia-buildpkg] |
| `julia-actions/julia-runtest` | `@v1` | inputs `check_bounds` (`yes`), `coverage` (`true`), `depwarn` (`yes`), `force_latest_compatible_version` (`auto`), `inline` (`yes`), `prefix` (`''`), `project` (`@.`), `annotate` (`false`) [CITED: raw action.yml, julia-actions/julia-runtest master] |

`v3` for both `setup-julia` and `cache` is newer than the `v2` that most blog posts and training data
still show — **do not write `@v2`.**

### Three concrete obstacles, all verified on disk

**(1) The co-resolution hard gate vs. dependency resolution.** `test/runtests.jl:32-49` asserts
`string(NeuralEstimators version) == "0.2.1"` and `string(Flux version) == "0.16.10"` — exact
strings. `Project.toml` `[compat]` permits `0.2.x ≥ 0.2.1` and `0.16.x ≥ 0.16.10`. `Manifest.toml`
**is committed** (git-tracked at the repo root) and pins `julia_version = "1.12.6"`,
`NeuralEstimators 0.2.1`, `Flux 0.16.10` [all VERIFIED].

Implications for the workflow, in order of importance:
- **Pin Julia exactly: `version: '1.12.6'`.** A different minor makes the committed manifest
  unusable (`Pkg.instantiate` will re-resolve or error) and the gate becomes a lottery. Do **not**
  use `'1'` or `'lts'`.
- **Prefer explicit `Pkg.instantiate()` over `Pkg.test()` for the fast tier.** `Pkg.test()` builds a
  sandbox environment; with the legacy `[extras]`/`[targets]` layout used here (there is **no**
  `test/Project.toml`) it re-resolves in that sandbox with tiered preservation — usually preserving
  the manifest versions, but not by contract. `Pkg.instantiate()` on the committed manifest installs
  the pinned tree-hashes and resolves nothing.
- Keep `force_latest_compatible_version` at its default (`auto` → false on ordinary pushes) if
  `julia-runtest` is used at all; **never set it to `true`** — it would deliberately break the gate.
- Optional hardening for the planner to weigh: tighten `[compat]` to `NeuralEstimators = "=0.2.1"`,
  `Flux = "=0.16.10"` so the declared contract matches the asserted one. This adds no dependency and
  removes a whole class of CI flake, but it edits `Project.toml`, which the phase constraints ask to
  leave alone. **Flag for the user; do not do it silently.**

**(2) GLMakie is a hard, unconditional dependency and CI is headless.** `src/ProteinCoLoc.jl:28`
does `import GLMakie` at module load, so **every** `using ProteinCoLoc` — including the fast tier's
golden script — loads GLMakie 0.13.12 (an OpenGL backend, with a `PrecompileTools` workload)
[VERIFIED: read `src/ProteinCoLoc.jl`, `Manifest.toml`]. On `ubuntu-latest` the job needs, before any
Julia step:

```yaml
- run: |
    sudo apt-get update
    sudo apt-get install -y xorg-dev mesa-utils xvfb libgl1 freeglut3-dev \
                            libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libxext-dev
```
and the Julia invocation should be prefixed with `xvfb-run -a` (`julia-runtest` exposes this as its
`prefix` input; a plain `run:` step just writes it inline). This mirrors Makie's own headless CI
[CITED: MakieOrg/Makie.jl issue #1953; GLMakie README "GLMakie's CI has no GPU"].
**This is the highest-risk unknown in the CI work and should be the first thing a plan proves**, with
a throwaway `using ProteinCoLoc; println("ok")` job before anything else is built on top of it.

**(3) The net is not in git — but it is fetchable, and the resolution path is already written.**
`artifacts/` is entirely git-ignored (`git ls-files artifacts/` is empty). `Artifacts.toml` declares
one lazy artifact `grid_8`, `git-tree-sha1 = 90e6b63a…`, downloaded from
`https://github.com/ma-seefelder/ProteinCoLoc/releases/download/v2.0.0/grid_8.tar.gz` (5.75 MB).
**Verified this session that the release tarball contains the amended_v2 net:** the tarball's
`npe_8.jld2` and `artifacts/amended_v2/grid_8/npe_8.jld2` share sha256
`198bb0783d35023fc4be6c04740589f6acb7d72fa047f3867fb91358b9042fed`, while the v1
`artifacts/grid_8/npe_8.jld2` hashes differently. `ProteinCoLoc.estimator_for(8)` →
`_lazy_load_from_artifact!(8)` resolves in three integrity-checked steps: artifact store →
in-repo `artifacts/amended_v2/grid_$(G)` dev dir with tree-sha1 verification → lazy download
(`registry.jl:96-147`). **So CI's golden script should call `ProteinCoLoc.estimator_for(8)`, not
`load_gate_model(path)`** — it gets the shipped net with hash verification and works both locally
(dev dir) and on a runner (download), and `julia-actions/cache@v3` caches the downloaded artifact by
default (`cache-artifacts: true`).

Note that `run_gate`'s default `artifacts_root` is `<pkgroot>/artifacts` — i.e. the **v1** `grid_8`
directory. Every Phase-15 invocation against the shipped net must pass
`--artifacts-root artifacts/amended_v2`, exactly as the amended confirmatory run did.

### Recommended two-tier shape (D-08)

```yaml
# .github/workflows/ci.yml            — FAST TIER
on:
  push: { branches: [main] }
  pull_request:
concurrency: { group: ci-${{ github.ref }}, cancel-in-progress: true }
jobs:
  calibration-golden:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: sudo apt-get update && sudo apt-get install -y xorg-dev mesa-utils xvfb libgl1 …
      - uses: julia-actions/setup-julia@v3
        with: { version: '1.12.6' }
      - uses: julia-actions/cache@v3
      - run: julia --project=. -e 'using Pkg; Pkg.instantiate()'
      - run: xvfb-run -a julia --project=. test/gate/ci_golden.jl        # D-09
```

```yaml
# .github/workflows/calibration-slow.yml — SLOW TIER
on:
  release: { types: [published] }
  workflow_dispatch:
    inputs:
      sweep_m: { description: 'SBC draws per rung', default: '1000' }
jobs:
  full-sbc:
    runs-on: ubuntu-latest
    timeout-minutes: 350        # GitHub-hosted job cap is 6 h — see Pitfall 6
    env: { JULIA_NUM_THREADS: '4' }   # ubuntu-latest is 4-core
    …
```

**Trigger discipline, recorded in the workflow file per D-08:** put a header comment in
`calibration-slow.yml` stating that it is release-triggered plus `workflow_dispatch`, that no
schedule is set on purpose (Deferred), and that its last successful run's report path is the evidence
the slow tier is alive. A `workflow_dispatch`-only workflow that nobody dispatches is silently dead,
which is the exact failure D-08 names.

### The D-09 Golden — feasible, with two named fragilities

**Feasible.** `SBC_FIX_M = 8`, `SBC_FIX_L = 15`, `SBC_FIX_BINS = 4` (`gate_consts_8_v2.jl:169-171`;
note `(15+1) % 4 == 0`, so the rank-bin evenness assertion holds). 8 draws × 3 pairs = 24
`simulate_pair` calls; at a small fixed `imsize` (pass a concrete tuple, which takes the scalar branch
of `gate_imsize` and consumes **zero** rng — `harness.jl:158-171`) this is seconds. Determinism comes
from both halves being pinned: the Philox `rng` for prior/simulation and `seed_gate_global!` for the
flow's base draws.

**Fragility 1 — the global RNG is not a cross-version contract.** `seed_gate_global!` calls
`Random.seed!` on Julia's default task-local RNG. Julia has changed the default RNG and the
`seed!`/task-local semantics across minor versions. A bit-exact golden is therefore a *Julia-patch*
assertion, not a *calibration* assertion. **Recommendation:** pin `version: '1.12.6'` in the workflow
(already required by the manifest), assert with a tolerance `P15_GOLDEN_TOL` rather than `==`, and
record in the re-bless procedure that a Julia version bump is a legitimate re-bless reason **and must
be named as such in the commit message**.

**Fragility 2 — a golden at M = 8 measures determinism, not calibration.** At M = 8 the ECE null mean
is ≈ 0.12; the number carries no calibration information. That is fine and should be stated in the
script's header: *the fast gate detects change, not miscalibration; the slow tier detects
miscalibration.* Anything else oversells SC3.

**What to assert (the discretion item), recommended:** all three, cheapest first —
1. `sha256` of `artifacts/…/npe_8.jld2` (or the resolved `Artifacts.artifact_hash("grid_8")`) — a net
   swap must fail loudly and instantly rather than as a numeric drift;
2. a **rank-table digest** (`hash` of the 8×8 Int matrix) — catches any change in the draw order,
   the summary chain or the flow forward pass, with zero tolerance ambiguity;
3. the **ECE scalar per column** against the committed golden within `P15_GOLDEN_TOL` — the
   human-readable signal, and the one that degrades gracefully across a Julia bump.

**Re-bless procedure (required deliverable, D-09):** a `--rebless` flag on `test/gate/ci_golden.jl`
that (a) refuses to run unless `P15_REBLESS_REASON` is passed as a non-empty string, (b) **appends**
a new sentinel-keyed Tier-2 block to `p15_consts.jl` rather than editing the previous golden, and
(c) prints the old and new values side by side for the commit message. Then guard it the way this
project guards constraints — with a running test that asserts the golden constants file has no
*modified* lines relative to `origin/main` for the golden block, mirroring
`spike/test/test_p12_decoupling.jl`'s byte-unchanged assertions and CONVENTIONS C-03's
"gate on changed **declaration** lines only" `git diff -U0 | grep -E '^[+-]\s*const '` recipe.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SBC ranks at a misspecified forward model | A parallel sweep harness | `sbc_gate(m; …, sim = misspec_simulator(gen; level, G))` | `misspec_simulator` already produces the exact `(sample_prior, simulate_pair, build_mci)` contract `sbc_gate` consumes |
| Per-rung OOD fire rate | A new scoring loop | `gate_ood_roc(…; families = P15_FAMILIES, levels = R)` | Returns `fire_rate[fam][lvl]` **and** the measured `id_fire_rate` — literally D-07's two inputs |
| ECE / coverage / KS / χ² | Any new statistic | `sbc_calibration` / `sbc_coverage` / `sbc_uniformity` | Already the shipped definitions; a second implementation makes the envelope incomparable to the gate report |
| Multiplicity correction | A new correction | `holm_adjusted` (frozen in the consts) / `sbc_holm` | Already frozen and unit-tested; but note D-04 makes p-values reported-only here |
| Atomic report writing | `JLD2.jldsave` directly | `write_gate_report(path, report)` | `.tmp` → reopen-integrity → `mv(...; force=true)`; a crash leaves a discardable `.tmp` |
| Loading the shipped net | `load_gate_model("artifacts/…")` | `ProteinCoLoc.estimator_for(8)` | Tree-sha1 verified, three-path resolution, works on a CI runner with no repo artifacts |
| Image-size draws | A second sampler | `gate_imsize` → `ProteinCoLoc.sample_imsize` | `harness.jl:126-130` — "the gate joint and the training joint cannot drift apart in implementation even if the constants agree" |
| Seed derivation | A literal constant | The salted-Philox redraw loop | D-11; a different **key** gives a different stream, a different index range does not (C-04) |

**Key insight:** this phase's entire novelty is the *sweep and its pre-registration*, not any of its
machinery. Every line of new statistical code is a line the plan-checker should question.

---

## Common Pitfalls

### Pitfall 1: Appending to `OOD_FAMILIES` silently rewrites a PASSED, frozen gate arm
**What goes wrong:** `OOD_FAMILIES` (`misspec.jl:181-184`) is consumed by `gate_ood_roc`'s
`pooled_pos` (strongest-level scores over **all** families) and by `ood_gate`'s
`passed = auc_ok && all(values(fam_pass))`. Adding `spillover` and `registration` changes
`combined_auc` and adds two entries to the pass conjunction — so any future run of the *shipped* OOD
arm produces a different, non-comparable number, and the recorded `OOD — PASS (pooled AUC 1.0)`
result stops being reproducible from the code.
**Why it happens:** CONTEXT.md's code-context section says "two new families … appended to
`OOD_FAMILIES`", which reads as a literal instruction.
**How to avoid:** define `P15_FAMILIES = merge(OOD_FAMILIES, (; spillover, registration))` in the
Phase-15 file and pass it through the existing `families =` keyword. `misspec.jl` stays byte-
unchanged; D-02's intent ("two more axes are swept") is fully satisfied; nothing frozen moves.
**Warning signs:** a diff that touches `test/gate/misspec.jl:181-184`.

### Pitfall 2: `background`'s ladder starts *outside* the prior, so D-03's marked boundary rung is not attainable with the shipped generator
**What goes wrong:** D-02 types `background` as in-prior extrapolation (SC1's autofluorescence) and
D-03 requires the ladder to span the prior boundary with the boundary as a marked rung. But
`misspec_background` sets `bleed = 0.1 + 0.3·level`, so **`level = 1` already gives 0.4 — four times
`maximum(AUTOFLUORESCENCE_PRIOR) = 0.1`** — and the generator additionally applies a multiplicative
illumination gradient and radial vignette that the simulator has **no parameter for at any θ**. On
mechanism, the shipped `background` family is a *compound* axis, closer to out-of-model than to
in-prior extrapolation.
**Why it happens:** the family was designed as an OOD positive control (make it clearly out of
distribution), not as a prior-boundary ladder.
**How to avoid — two admissible routes, both must be stated in the report, not chosen silently:**
 (a) **Honest limit (recommended, no new code):** sweep `background` at its native levels 1..R, type
 it **out-of-model / compound** in the domain map with the reason, and record that its prior boundary
 lies *below rung 1* so D-03's marked rung for this axis is rung 0 (the in-prior anchor) rather than
 an interior rung. Note the offset value at each rung so the compound nature is auditable.
 (b) **Pure-offset ladder:** add a `misspec_autofluorescence` generator to the Phase-15 file with
 `merge(θ, (; autofluorescence = P15_AF_LADDER[level]))` and `P15_AF_LADDER[2] == 0.1` — clean
 typing, exact boundary rung, level-independent rng consumption (same trick as spillover). Cost: it
 is a **seventh** axis or a replacement for `background`, and D-02 froze six. **This needs a user
 ruling, not a planner decision.**
**Warning signs:** a domain map that reports "background breaks at rung 1" and reads it as an
in-prior training failure.

### Pitfall 3: The `grid_16` contrast is blocked by a shipped invariant and is triply confounded
**What goes wrong:** `run_gate` calls `assert_imsize_provenance` and, under
`SBC_REQUIRE_IMSIZE_PROVENANCE = true`, returns `(status = :provenance_mismatch, …)` **running no
arm** when the net's persisted provenance is absent. Verified this session:
`artifacts/grid_16/npe_16.jld2` → `(imsize_set = :unknown, imsize_weights = :unknown,
imsize_source = :unknown, recorded = false)`. There is no `artifacts/amended_v2/grid_16`.
**Why it happens:** grids 4/8/16 were trained on 2026-07-20 under the v1 regime (256², unbounded θ);
only grid 8 was retrained on 2026-07-21 for the amendment.
**How to avoid:** do not plan a grid_16 sweep until this is ruled. Bypassing the guard
(`require_provenance = false`) is *possible* but the resulting contrast varies **grid dimension ×
image-size regime × θ-space parameterisation** simultaneously — it cannot support D-10's stated
inference ("evidence the envelope is not a grid-8 artifact"). See Open Question 1.
**Warning signs:** a plan task that invokes the sweep at `--grid 16` without naming the confound.

### Pitfall 4: `MCE` is pinned at 0.99 by empty reliability bins and means nothing
**What goes wrong:** `_bin_calibration` (`sbc.jl:99-104`) computes `mce` over **all** `n_bins`
including empty ones, where `predicted_rate = midpoint` and `observed_rate = 0.0`. With 19 α values
in 50 bins, 31 bins are empty and the top one has midpoint 0.99. Verified: **all eight columns of the
shipped report record `mce = 0.99`.**
**Why it happens:** the function was ported from a reliability-diagram context where predicted
probabilities populate the whole [0,1] range.
**How to avoid:** never use `mce` as a Phase-15 statistic; if a worst-level statistic is wanted,
compute `max|α − ĉ(α)|` over the 19 levels directly (report-only, since D-04 fixes ECE as the break).
**Warning signs:** a table where MCE is identical across every column and every rung.

### Pitfall 5: Forgetting the second half of reproducibility
**What goes wrong:** `posterior_for` → `NeuralEstimators.sampleposterior` accepts **no `rng`**; the
flow's base draws come from the **global** RNG (`harness.jl:49-60`). The passed Philox stream drives
only `sample_prior`, `gate_imsize` and `simulate_pair`. A Phase-15 arm that constructs its own loop
without calling `seed_gate_global!` produces numbers that are not reproducible — silently, and only
discoverable by trying to re-run.
**Why it happens:** the promotion of the harness into `test/gate/` already dropped this once
(recorded as a fixed regression in the same comment block), and the grid-4/8/16 v1 rank tables are
permanently non-reproducible because of it.
**How to avoid:** call `seed_gate_global!(G)` at the top of the new `--envelope` arm, and — since
`gate_global_seed` derives from `prod_seed(G)` — define `prod_seed`/`prod_rng` aliases in
`p15_consts.jl` exactly as `gate_consts_8_v2.jl:342-343` does, so the global half rides the P15
stream automatically. Assert in a unit test that `gate_global_seed(8)` under `p15_consts.jl` differs
from its value under `gate_consts_8_v2.jl`.
**Warning signs:** an arm that constructs `Philox4x` inline instead of going through `prod_rng`.

### Pitfall 6: The slow tier can exceed a GitHub-hosted job's wall-clock cap
**What goes wrong:** a full-scale sweep is 7.3 h at 32 threads; `ubuntu-latest` GitHub-hosted runners
are **4-core** and jobs are capped at **6 hours**. Naively porting the sweep into the slow tier
guarantees a timeout, on a job that has already burned hours.
**How to avoid:** the slow tier should run the **reported-scale SBC/coverage on the shipped net**
(the ~22–30 min arm the confirmatory run demonstrates), **not** the 6-axis envelope sweep. The
envelope sweep is a one-off local/self-hosted run whose artifact is committed as a report, exactly
like `gate-8x8-amended.md`. Say this explicitly in the workflow header so nobody later "fixes" it by
moving the sweep into CI.
**Warning signs:** a `calibration-slow.yml` with a matrix over axes.

### Pitfall 7: Cross-rung ECEs are not independent replicates
**What goes wrong:** `seed_gate_global!(G)` depends only on `G`, so every rung's flow base-draw
stream is identical; and for `spillover`/`registration` the *entire* rng consumption is
level-independent, so those ladders are matched-pairs. Treating rung ECEs as independent samples —
e.g. applying a multiplicity correction across rungs, or quoting an independent CI per rung — is
wrong in both directions.
**How to avoid:** state the pairing in the report; use it as the *reason* the ladder is readable at
M = 1000 (differences are far less noisy than levels), and do **not** correct for multiplicity across
rungs. If independent rungs are wanted instead, pass a per-rung `rng = p15_rng(G, axis, rung)`
(distinct Philox counters) — but the matched design is the better one for a monotone ladder.

---

## Runtime State Inventory

Phase 15 is additive (new files + one new arm + a new workflow directory); it is not a rename,
refactor or migration. The inventory is included anyway because the phase writes artifacts and
consumes seeds.

| Category | Items Found | Action Required |
|----------|-------------|-----------------|
| Stored data | `artifacts/amended_v2/grid_8/*.jld2` (read-only input); the phase writes a **new** `envelope_report_8.jld2` alongside. `artifacts/` is git-ignored — verified `git ls-files artifacts/` is empty. | none — write a new file, never overwrite `gate_report_8.jld2` |
| Live service config | **None** — verified: no `.github/`, no `.git/hooks` beyond samples, no `Makefile`, no external service configuration in the repo. | none (the phase *creates* the first `.github/`) |
| OS-registered state | **None** — verified: no scheduled tasks, daemons or service registrations referenced anywhere in the repo. | none |
| Secrets / env vars | **None required.** The fast tier needs only the default `GITHUB_TOKEN` (used by `julia-actions/cache` for cache deletion). `JULIA_NUM_THREADS` is a runner env var, not a secret. | none |
| Build artifacts | Julia precompile cache (`~/.julia/compiled`) — cached by `julia-actions/cache@v3` (`cache-compiled: true`); the `grid_8` Julia artifact (5.75 MB) — cached by `cache-artifacts: true`. | none; both regenerate |
| **Seeds (project-specific)** | `PROD_SEED_V2[8] = 1135605683775656488` is **SPENT** (recorded in `gate_report_8.jld2`). `PROD_SEED[4,8,16,32]` spent. `DEV_SEEDS`, `VAL_MASTER_SEED`, `NPE_MASTER_SEED`, `DEFAULT_MASTER_SEED`, `VAL_FIX_SEED`, `RATIO_PAIR_SEED`, `CORPUS_MASTER_SEED`, `P11/P12/P13 DEV+FIX` all burned. | derive `PROD_SEED_P15` fresh and assert disjoint from **all** of them (below) |

---

## Seed Derivation (D-11) — exact construction and exact assertion set

### The construction to copy

Verbatim shape from `gate_consts_8_v2.jl:266-343`, with fresh mixing constants:

```julia
# --- forbidden set (strictly larger than v2's) ---
const NPE_MASTER_SEED     = 0x0000_0000_00C0_FFEE
const VAL_MASTER_SEED     = 0x0000_0000_5BC0_FFEE
const VAL_FIX_SEED        = 0x0000_0000_00F1_F7ED
const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001
const RATIO_PAIR_SEED     = 0x0000_0000_004A_7107
const CORPUS_MASTER_SEED  = 0x0000_0000_00C0_5EED
const DEV_SEEDS           = (0x0000_0000_0DE7_C0DE, 0x0000_0000_DE7C_0DE2)
const P11_DEV_SEED        = 0x0000_0000_0B11_DE71
const P12_DEV_SEED        = 0x0000_0000_0B12_DE71
const P12_FIX_SEED        = 0x0000_0000_0B12_F1F7
const P13_DEV_SEED        = 0x0000_0000_0B13_DE71
const P13_FIX_SEED        = 0x0000_0000_0B13_F1F7

# v1 and v2 gate seeds are DERIVED, so RECOMPUTE them here (never trust a comment):
const PROD_SALT   = 0x94D0_49BB_1331_11EB;  const PROD_MASTER  = 0x0000_0000_09E3_779B
const AMEND_SALT  = 0xC4CE_B9FE_1A85_EC53;  const AMEND_MASTER = 0x0000_0000_C0DE_2026
_derive(master, salt, G) = (rng = Philox4x(UInt64, (UInt64(master) ⊻ salt, UInt64(G)));
                            s = rand(rng, UInt64);
                            while s in (VAL_MASTER_SEED, NPE_MASTER_SEED, UInt64(0))
                                s = rand(rng, UInt64) end; s)
const PROD_SEED    = Dict(G => _derive(PROD_MASTER,  PROD_SALT,  G) for G in (4,8,16,32))
const PROD_SEED_V2 = Dict(G => _derive(AMEND_MASTER, AMEND_SALT, G) for G in (4,8,16,32))

_p15_forbidden() = (UInt64(0), NPE_MASTER_SEED, VAL_MASTER_SEED, VAL_FIX_SEED,
                    DEFAULT_MASTER_SEED, RATIO_PAIR_SEED, CORPUS_MASTER_SEED, DEV_SEEDS...,
                    P11_DEV_SEED, P12_DEV_SEED, P12_FIX_SEED, P13_DEV_SEED, P13_FIX_SEED,
                    values(PROD_SEED)..., values(PROD_SEED_V2)...)

# --- the FRESH stream: distinct (master, salt) pair + redraw loop ---
const P15_SALT   = <fresh 64-bit, NOT in P15_REPO_SALTS below>
const P15_MASTER = <fresh 64-bit>
function _derive_prod_seed_p15(G::Integer)
    rng = Philox4x(UInt64, (UInt64(P15_MASTER) ⊻ P15_SALT, UInt64(G)))
    s = rand(rng, UInt64)
    while s in _p15_forbidden(); s = rand(rng, UInt64); end
    return s
end
const PROD_SEED_P15 = Dict{Int,UInt64}(G => _derive_prod_seed_p15(G) for G in (4, 8, 16, 32))
prod_seed_p15(G) = get(() -> _derive_prod_seed_p15(G), PROD_SEED_P15, Int(G))
prod_rng_p15(G)  = Philox4x(UInt64, (prod_seed_p15(G), UInt64(0)))
# The aliases that make the EXISTING gate machinery ride the P15 stream (v2 pattern, :342-343):
prod_seed(G::Integer) = prod_seed_p15(G)
prod_rng(G::Integer)  = prod_rng_p15(G)
```

**`P15_SALT` must not be any of these repo salts** [VERIFIED: `p12_consts.jl:234-245`]:
`0x9E3779B97F4A7C15` (HOLDOUT), `0xD1B54A32D192ED03` (FOLD/CBS_SPLIT), `0xBF58476D1CE4E5B9` (VAL),
`0x94D049BB133111EB` (PROD), `0xC4CEB9FE1A85EC53` (AMEND), `0xA24BAED4663EE121` (P11),
`0xA5A5A5A5A5A5A5A5` and `0xA5A5A5A5DEADBEEF` (COSTES), `0x2545F4914F6CDD1D` (P13 == P11_DATAGEN),
`0xFF51AFD7ED558CCD` (P12), `0xC2B2AE3D27D4EB4F` (P12_DATAGEN).

### The exact assertion set

Both in-file `@assert`s (the `gate_consts_8_v2.jl:450-453` pattern) **and** an independent testset in
`test/runtests.jl` loading the file into an **isolated module** — mandatory, because `p15_consts.jl`
defines the same const names (`SBC_M`, `prod_seed`, …) that `run_gate.jl` already loaded into the
test module via the template (`runtests.jl:878`). Copy the `module GateConstsV2 … end` wrapper at
`runtests.jl:1147-1149`.

```julia
module P15Consts
    include(joinpath(@__DIR__, "gate", "p15_consts.jl"))
end
@testset "Phase-15 pre-registration (p15_consts.jl)" begin
    P = P15Consts; fresh = P.PROD_SEED_P15[8]
    @test fresh != 0
    @test fresh != P.NPE_MASTER_SEED
    @test fresh != P.VAL_MASTER_SEED
    @test fresh != P.DEFAULT_MASTER_SEED
    @test !(fresh in P.DEV_SEEDS)
    @test fresh != P.PROD_SEED_V2[8]                      # THE SPENT SEED — D-11 names it
    for G in (4, 8, 16, 32)
        @test fresh != P.PROD_SEED[G]
        @test fresh != P.PROD_SEED_V2[G]
        @test P.PROD_SEED_P15[G] != P.PROD_SEED[G]
        @test P.PROD_SEED_P15[G] != P.PROD_SEED_V2[G]
    end
    @test isempty(intersect(Set(values(P.PROD_SEED_P15)), Set(values(P.PROD_SEED))))
    @test isempty(intersect(Set(values(P.PROD_SEED_P15)), Set(values(P.PROD_SEED_V2))))
    @test length(unique(values(P.PROD_SEED_P15))) == 4
    @test !(fresh in P._p15_forbidden())
    @test !(P.P15_SALT in P.P15_REPO_SALTS)               # distinct KEY, not just a distinct index
    # the stream itself
    @test P.prod_seed(8) == fresh
    @test rand(P.prod_rng(8), UInt64, 4) == rand(P.prod_rng(8), UInt64, 4)   # reproducible
    @test rand(P.prod_rng(8), UInt64, 4) != rand(P.prod_rng(4), UInt64, 4)   # per-grid disjoint
    # the GLOBAL half (Pitfall 5) rides the P15 stream, not the spent one
    @test P.prod_seed(8) != GateConstsV2.PROD_SEED_V2[8]
    # C-04 immunity is a KEY difference, shown by construction:
    @test (UInt64(P.P15_MASTER) ⊻ P.P15_SALT) != (UInt64(GateConstsV2.AMEND_MASTER) ⊻ GateConstsV2.AMEND_SALT)
end
```

**Show the assertion fires.** CONVENTIONS C-04 rule 3: *"an `isempty(intersect(...))` never shown
non-empty is an assurance, not a test."* Add one negative fixture — e.g. a locally-constructed
`_derive_prod_seed_p15` variant seeded with `AMEND_MASTER`/`AMEND_SALT`, asserted to land inside
`_p15_forbidden()` — so the guard is demonstrated to catch a real collision.

---

## Code Examples

### Reading the ECE anchor from the shipped report (no re-measurement needed)

```julia
# Source: verified against artifacts/amended_v2/grid_8/gate_report_8.jld2 this session
import JLD2
rep = JLD2.jldopen("artifacts/amended_v2/grid_8/gate_report_8.jld2", "r") do f
    f["report"]                                   # top-level keys: schema_version, report
end
anchor = Dict(p.label => p.ece for p in rep.sbc.per_param)
# anchor["ρ_true"] == 0.03713157894736839 ; anchor["Δρ"] == 0.0042894736842104715
rep.sbc.ranks           # the FULL 2000x8 Int rank table — any rank statistic is recomputable
rep.sbc.M, rep.sbc.L    # 2000, 999   <- the M the anchor was measured at
```

### The ECE identity (use this, not a re-derivation)

```julia
# Source: derived from test/gate/sbc.jl:302-315 + :60-108, verified numerically this session
ece(ranks, L; levels = 0.05:0.05:0.95) = Statistics.mean(
    abs(α - Statistics.mean(abs.((ranks .+ 0.5) ./ (L + 1) .- 0.5) .<= α / 2)) for α in levels)
```

### The null-ECE margin (computable before any rung runs)

```julia
# Source: this session's Monte-Carlo probe; reproduces E[ECE_0] ~= 0.336/sqrt(M)
function ece_null_quantiles(M, L; B = 20_000, q = 0.95, rng = Random.Xoshiro(20260804))
    v = [ece(rand(rng, 0:L, M), L) for _ in 1:B]
    (mean = Statistics.mean(v), q = Statistics.quantile(v, q))
end
# M=1000, L=999 -> (mean = 0.01038, q95 = 0.01884)  =>  P15_ECE_MARGIN = 0.00846
```

### The sweep loop (orchestration only — every call already exists)

```julia
# Source: composed from test/gate/sbc.jl:437, misspec.jl:194-199, misspec.jl:316
rows = NamedTuple[]
for (axis, gen) in pairs(P15_FAMILIES), rung in 0:P15_RUNGS
    rung == 0 && axis != first(keys(P15_FAMILIES)) && continue      # rung 0 measured ONCE
    g = sbc_gate(m; G = 8, M = P15_SWEEP_M, L = SBC_L,
                 imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
                 imsize_weights = GATE_IMSIZE_WEIGHTS,
                 sim = rung == 0 ? default_simulator() :
                       misspec_simulator(gen; level = rung, G = 8),
                 rng = prod_rng(8))
    push!(rows, (axis = axis, rung = rung,
                 mechanism   = P15_AXIS_MECHANISM[axis],       # :in_prior | :out_of_model
                 prior_edge  = (rung == P15_PRIOR_BOUNDARY_RUNG[axis]),
                 ece_rho     = g.per_param[1].ece,             # D-05 GATING
                 ece_drho    = g.per_param[8].ece,             # D-05 GATING
                 ks_p        = [x.ks_p   for x in g.per_param], # reported, NOT gating
                 chi2_p      = [x.chi2_p for x in g.per_param], # reported, NOT gating
                 coverage    = sbc_coverage(view(g.ranks, :, 1); L = SBC_L),
                 shrinkage   = [x.shrinkage for x in g.per_param],
                 realised    = g.realised_imsize_counts))
end
```

---

## State of the Art

| Old approach | Current approach | When changed | Impact on this phase |
|--------------|------------------|--------------|----------------------|
| `julia-actions/setup-julia@v1`/`@v2` | **`@v3`** | current major | Write `@v3`; `@v2` still works but is not current |
| `julia-actions/cache@v1`/`@v2` | **`@v3`** | current major | Write `@v3` |
| `actions/cache` hand-rolled for `~/.julia` | `julia-actions/cache@v3` (five depot subdirs, incl. compiled + artifacts) | — | Caches the 5.75 MB `grid_8` artifact and the precompile cache for free |
| `setup-julia` numeric versions only | named `lts`/`pre`/`min`/`min-minor`/`min-patch` (since v2.2.0) | 2024 | **Do not use them here** — the committed manifest pins 1.12.6 |
| Legacy `[extras]`/`[targets]` test deps | `test/Project.toml`, and workspaces on Julia 1.12+ | Julia 1.12 | This repo still uses the legacy form; that is *why* `Pkg.instantiate()` is preferred over `Pkg.test()` for the deterministic tier |

**Deprecated / not applicable:**
- `Pkg.test()` as the fast tier's entry point — it sandboxes and re-resolves; the co-resolution gate
  asserts exact versions that `[compat]` does not pin.
- `mce` from `CalibrationResult` — structurally degenerate here (Pitfall 4).

---

## Assumptions Log

| # | Claim | Section | Risk if wrong |
|---|-------|---------|---------------|
| A1 | The 32-thread cost extrapolation (`t_draw ≈ 0.656 s`) holds on the machine that will run the sweep | Compute Budget | Wall-clock off by the thread-count ratio (up to 13× if run single-threaded). **Mitigation: the plan should time one rung at M = 50 and reconcile before committing the grid.** Reconciles to the one observed run at 30.9 min, so the model itself is sound. |
| A2 | Misspecification generators average ≈ 1.3× `simulate_pair` cost | Compute Budget | Budget error up to ~1.5× on the `optics` axis (large-σ `imfilter`). Not measured; measurable in minutes. |
| A3 | The `patch_summary` G=16 timings (1.6–1.9 s) are GC/JIT-contaminated rather than real | Compute Budget | If real, the grid_16 contrast costs ~4× the grid_8 sweep. **Re-measure with BenchmarkTools before budgeting grid_16.** |
| A4 | `Pkg.instantiate()` on the committed manifest installs the pinned versions without re-resolving | CI | If wrong, the co-resolution gate could fail on an upstream patch release. Directly testable in the first CI run. |
| A5 | GLMakie 0.13.12 precompiles and loads on headless `ubuntu-latest` given xvfb + the listed apt packages | CI, Pitfall list | If wrong, the **entire** CI tier is blocked (GLMakie is imported unconditionally at package load). Highest-risk CI unknown; prove it first with a trivial `using ProteinCoLoc` job. |
| A6 | `n_pos = 40` is too small for a fire-rate margin and should be raised to ~200 | OOD Arm | Under-powered D-07 detection (SE 0.034 at p = 0.05). The cost of raising it is ~20 min, so the downside of being wrong is trivial. |
| A7 | `julia-runtest`'s `test_args` input exists on the current `v1` tag | CI | The README documents it; the fetched `action.yml` did not list it. Immaterial under the recommended design, which does not use `julia-runtest` for the fast tier. |
| A8 | The GitHub-hosted job wall-clock cap is 6 h and `ubuntu-latest` is 4-core | Pitfall 6 | If the cap is different, the slow-tier scoping argument shifts but the conclusion (do not run the 7.3 h sweep in CI) does not. |

---

## Open Questions

1. **The D-10 `grid_16` contrast — the only planning blocker that needs a user ruling.**
   - *What we know:* `artifacts/grid_16/npe_16.jld2` exists but reports `recorded = false` training
     provenance; `run_gate` refuses (`:provenance_mismatch`) under the amended pre-registration; there
     is no `amended_v2/grid_16`; and retraining is forbidden (07-GATE-AMENDMENT §6.4, and D-10 itself
     says the envelope is a claim about the *shipped* net).
   - *What's unclear:* whether a confounded contrast is worth running at all.
   - *Recommendation:* present three options rather than picking one.
     (i) **Drop the grid_16 contrast** and report the grid_8 envelope alone, with the reason recorded
     as a named limit — consistent with Phase 11/12's negative-but-useful precedent and costs nothing.
     (ii) **Run it with `require_provenance = false`** and report it explicitly as a *three-way
     confounded* contrast (grid × imsize regime × θ-space), which is honest but supports a much
     weaker inference than D-10 states.
     (iii) **Defer to v2.1** alongside the already-deferred grid_4 trend.
     My reading of the evidence favours (i) or (iii); (ii) risks a number that a referee will
     correctly discount, and D-10's stated purpose ("evidence the envelope is not a grid-8 artifact")
     is not served by a contrast that varies three things at once.

2. **`background`'s prior boundary (Pitfall 2)** — route (a) *honest limit, no new code* or route (b)
   *pure-offset `misspec_autofluorescence` generator*. Route (b) gives D-03 exactly what it asks for
   on this axis at the cost of a seventh axis (or replacing `background`), which is a D-02 change.
   Needs a ruling before `p15_consts.jl` freezes the ladder endpoints.

3. **`P15_ECE_MARGIN`'s derivation form.** Recommended: `q95(ECE₀) − E[ECE₀]` at `P15_SWEEP_M`
   (= 0.00846 at M = 1000), which gives a ~5 % per-rung false-break rate under exact calibration. An
   equally defensible alternative is `q99 − E` (= 0.01425), giving ~1 %. Both are pure design
   arithmetic; the choice is a false-break/miss tradeoff and should be made *and stated with its
   operating characteristic* before any rung runs, not after.

4. **Whether to tighten `[compat]` to `"=0.2.1"` / `"=0.16.10"`.** It would make the declared contract
   match the asserted one and remove a CI flake class at zero dependency cost, but it edits
   `Project.toml`, which the phase constraints ask to leave alone. Flagged, not recommended either way.

---

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|-------------|-----------|---------|----------|
| Julia | everything | ✓ | **1.12.6** (matches `Manifest.toml` `julia_version`) | — |
| `ProteinCoLoc` package env | sweep + golden | ✓ | resolves; `using` costs ~91 s cold, ~20 s warm | — |
| `artifacts/amended_v2/grid_8/` | the reported envelope | ✓ | npe/ratio/ood_nulls/gate_report present; tree matches `Artifacts.toml` `grid_8` | lazy download from the v2.0.0 GitHub release |
| `artifacts/grid_16/` | D-10 contrast | ✓ present, **✗ usable** | v1 bundle, `recorded = false` provenance | see Open Question 1 |
| `NeuralEstimators` / `Flux` | inference | ✓ | 0.2.1 / 0.16.10 (manifest) | — |
| `Random123`, `HypothesisTests`, `JLD2`, `ImageFiltering` | sweep | ✓ | in manifest | — |
| GPU / CUDA | — | not required | — | CPU-only is the contract (`use_gpu = false` throughout the harness) |
| `.github/` + Actions runner | SC3 | ✗ **absent** | — | none — this phase creates it |
| OpenGL / xvfb on the runner | GLMakie at package load | unknown (no CI yet) | — | none; the apt package list above is the only route |
| Multi-core host for the sweep | 7.3 h budget | ✓ (this session measured `Threads.nthreads() == 1` by default) | — | single-threaded is ~13× slower — **set `JULIA_NUM_THREADS`** |

**Missing dependencies with no fallback:** `.github/` (created by this phase); a proven headless
GLMakie load on `ubuntu-latest` (must be established by the first CI task, before anything is built on it).

**Missing dependencies with fallback:** the shipped net on a runner — resolved via
`Artifacts.toml` lazy download, cached by `julia-actions/cache@v3`.

---

## Validation Architecture

`.planning/config.json` does not set `workflow.nyquist_validation` to `false`, so this section is
included.

### Test Framework

| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (`@testset`/`@test`), driven by `test/runtests.jl` (1550 lines) |
| Config file | none — legacy `[extras] Test` + `[targets] test = ["Test"]` in `Project.toml`; **no `test/Project.toml`** |
| Quick run command | `julia --project=. -e 'include("test/gate/<file>.jl")'` for a single file, or a dedicated script; there is **no** per-testset selector today |
| Full suite command | `julia --project=. -e 'using Pkg; Pkg.test()'` (heavy: loads Turing ext, Images, GLMakie, runs `test_integration.jl` and `test_local_map.jl`) |
| Isolation pattern for a consts file | `module P15Consts; include(".../p15_consts.jl"); end` — required, mirroring `runtests.jl:1147-1149` |

### Phase Requirements → Test Map

| Req | Behaviour | Test type | Automated command | File exists? |
|-----|-----------|-----------|-------------------|--------------|
| D-11 | `PROD_SEED_P15` disjoint from all burned streams; assertion demonstrably fires | unit | `julia --project=. -e 'include("test/gate/p15_consts_test.jl")'` | ❌ Wave 0 |
| D-11 | `seed_gate_global!` rides the P15 stream, not `PROD_SEED_V2` | unit | same file | ❌ Wave 0 |
| D-02 | New generators satisfy the `(rng, θ; imsize, level, G)` contract, return 2×`Matrix{Float64}` of the right size, finite, ≥ 0 | unit | `.../p15_misspec_test.jl` | ❌ Wave 0 |
| D-02 | **Rung-0 identity**: `misspec_spillover(rng, θ; level=0) === simulate_pair(rng′, θ)` bit-for-bit | unit | same file | ❌ Wave 0 |
| D-02 | Generators are monotone in `level` (effective magnitude strictly increasing) | unit | same file | ❌ Wave 0 |
| D-02 | `test/gate/misspec.jl` is **byte-unchanged** (Pitfall 1) | structural | `git diff --quiet HEAD~ -- test/gate/misspec.jl` in the plan's verification | ❌ Wave 0 |
| D-03 | `P15_SPILLOVER_LADDER[2] == maximum(SPILLOVER_PRIOR)` and `P15_SHIFT_LADDER[2] == maximum(SHIFT_PRIOR)` — the boundary is derived, never a duplicated literal | unit (in-file `@assert` + test) | `.../p15_consts_test.jl` | ❌ Wave 0 |
| D-04 | ECE identity: `sbc_calibration(...).ece == mean|α − ĉ(α)|` on a fixture rank vector | unit | `.../p15_ece_test.jl` | ❌ Wave 0 |
| D-04 | Null-ECE margin reproduces `≈ 0.336/√M` at M ∈ {250, 1000} | unit (fast MC, B = 2000) | same file | ❌ Wave 0 |
| D-04 | Anchor is read from `gate_report_8.jld2`, not hard-coded | integration | `.../p15_anchor_test.jl` (skips honestly if the artifact is absent, mirroring `:not_trained`) | ❌ Wave 0 |
| D-05 | Break decision reads **only** columns 1 and 8 | unit | `.../p15_envelope_test.jl` with a synthetic rank table | ❌ Wave 0 |
| D-06 | Three-valued SC2 classifier returns PROTECTIVE / LATE / SILENT-BUT-SAFE on three synthetic fixtures | unit | same file | ❌ Wave 0 |
| D-07 | Fire-rate rule compares against the **measured** `id_fire_rate`, not `1 − OOD_ID_QUANTILE` | unit | same file | ❌ Wave 0 |
| D-09 | Golden script is deterministic across two runs in one process and across two processes | smoke | `julia --project=. test/gate/ci_golden.jl` ×2, diff | ❌ Wave 0 |
| D-09 | `--rebless` refuses without a reason string and **appends** rather than edits | unit | `.../p15_golden_test.jl` | ❌ Wave 0 |
| D-12 | `P15_GATING_CONSTANTS ∩ P15_REPORTING_ONLY_CONSTANTS == ∅`; every Tier-1 bar names its derivation | unit | `.../p15_consts_test.jl` | ❌ Wave 0 |
| D-13 | `P15_ITERATION_ALLOWANCE == 1` and the ladder endpoints are present and frozen | unit | same file | ❌ Wave 0 |
| SC3 | The fast workflow actually fails on injected calibration drift | manual-only (needs a real PR to observe) | — | n/a — justified: a workflow's failure behaviour is only observable in Actions |

### Sampling Rate

- **Per task commit:** the Phase-15 unit files only — seconds, no net, no simulation:
  `julia --project=. -e 'for f in ("p15_consts_test","p15_misspec_test","p15_ece_test","p15_envelope_test"); include("test/gate/$f.jl"); end'`
- **Per wave merge:** the above plus the D-09 golden (`ci_golden.jl`, seconds) plus a **single**
  `M = 50` rung as a live smoke of the sweep arm (~40 s at 32 threads).
- **Phase gate:** full `Pkg.test()` green (which includes the co-resolution gate), the fast CI
  workflow green on a real PR, and the envelope sweep report written.

### Wave 0 Gaps

- [ ] `test/gate/p15_consts.jl` — the pre-registration itself (D-12), must land **before** any rung
- [ ] `test/gate/p15_consts_test.jl` — seed disjointness (incl. a firing negative fixture), ladder/prior-boundary derivation, tier split
- [ ] `test/gate/p15_misspec.jl` + `p15_misspec_test.jl` — the two generators, `P15_FAMILIES`, rung-0 identity
- [ ] `test/gate/p15_ece_test.jl` — ECE identity + null-margin derivation
- [ ] `test/gate/p15_envelope_test.jl` — break/OOD-crossing reduction on synthetic fixtures
- [ ] `test/gate/ci_golden.jl` + `p15_golden_test.jl` — the fast gate and its re-bless discipline
- [ ] Wiring: a `@testset` block in `test/runtests.jl` that loads `p15_consts.jl` into an isolated module
- [ ] No framework install needed — `Test` is already the target.

---

## Security Domain

`security_enforcement` is not set to `false` in `.planning/config.json`, so this section is included.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard control |
|---------------|---------|------------------|
| V2 Authentication | no | No user-facing auth surface; CI uses the default `GITHUB_TOKEN` only |
| V3 Session Management | no | No sessions |
| V4 Access Control | **yes (CI)** | `permissions:` block in each workflow, least-privilege (`contents: read`; `actions: write` only if `julia-actions/cache`'s `delete-old-caches` is kept) |
| V5 Input Validation | **yes** | `_guard_misspec_imsize` on every generator; `simulate_pair`'s θ-range `throw`s; `workflow_dispatch` inputs must be validated (parse `sweep_m` as an integer with bounds, never interpolate into a shell string) |
| V6 Cryptography | **yes (integrity, not secrecy)** | `Artifacts.toml` `git-tree-sha1` + `sha256` verification via `_verify_tree_sha1` / `ensure_artifact_installed` — already implemented, do not bypass |
| V14 Configuration | **yes** | Pin every action to a major tag at minimum (`@v3`), pin Julia to `1.12.6`, pin the runner image explicitly |

### Known Threat Patterns for a Julia GitHub Actions workflow

| Pattern | STRIDE | Standard mitigation |
|---------|--------|---------------------|
| Untrusted PR code executing with write permissions | Elevation of Privilege | Use `pull_request` (not `pull_request_target`) for the fast tier; set `permissions: contents: read` at workflow level |
| Script injection via `${{ github.event.* }}` into `run:` | Tampering | Never interpolate event data into shell; pass through `env:` and quote |
| Supply-chain drift in a floating action tag | Tampering | Pin to major tags now (`@v3`/`@v1`); consider SHA-pinning before v2.0 release |
| Model-artifact substitution (a swapped `npe_8.jld2`) | Tampering | `_verify_tree_sha1` already refuses a mismatched bundle; D-09's golden should additionally assert the artifact hash |
| Cache poisoning across branches | Tampering | `julia-actions/cache` scopes by workflow+job by default; do not widen `cache-name` across trust boundaries |
| Unbounded compute from `workflow_dispatch` | Denial of Service | `timeout-minutes` on the slow-tier job; validate `sweep_m` bounds |

Nothing in this phase handles personal data, credentials, or network input beyond the pinned,
hash-verified artifact download.

---

## Sources

### Primary (HIGH confidence — read or executed this session)

- `test/gate/misspec.jl` (396 lines, read in full) — generator contract, `OOD_FAMILIES`, `gate_ood_roc` return shape
- `test/gate/sbc.jl` (496 lines, read in full) — `sbc_ranks`/`_and_spread`/`uniformity`/`coverage`/`calibration`/`gate`, `_bin_calibration`, `SBC_PARAM_LABELS`
- `test/gate/harness.jl` (299 lines, read in full) — `draw_simulate_infer`, `gate_global_seed`/`seed_gate_global!`, the global-RNG caveat, `gate_imsize`, `assert_imsize_provenance`
- `test/gate/run_gate.jl` (541 lines, read in full) — CLI, `--consts` override, the `PROD_SEED_V2[8]` one-run binding, `bf_gate`, `ood_gate`, `write_gate_report`
- `test/gate/gate_consts_8_v2.jl` (454 lines, read in full) — seed derivation + redraw loop, `SBC_*` constants, F5 mixture and its cost table, `SBC_FIX_*`
- `src/amortized/simulator.jl` (read in full) — all prior objects and bounds, `σ_psf = 1.3`, `simulate_pair` stage order and validation
- `src/amortized/ood.jl` (relevant sections) — `OOD_ID_QUANTILE`, `id_threshold`, `roc_auc`, `ood_verdict`
- `src/registry.jl:88-160` — `_dev_bundle_dir` = `artifacts/amended_v2/grid_$(G)`, `_lazy_load_from_artifact!` three-path resolution, `_verify_tree_sha1`
- `spike/validation/p12_consts.jl:21-95, 99-260, 723-725` — the two-tier template, the forbidden-seed and repo-salt lists
- `test/runtests.jl:32-49, 878, 1137-1210` — the co-resolution hard gate, the template include, the isolated-module pattern
- `Project.toml`, `Manifest.toml`, `Artifacts.toml`, `.gitignore`, `git ls-files` — dependency and artifact provenance
- **Executed:** `ProteinCoLoc` timing probe (simulate/summary/posterior at four image sizes, `nthreads = 1`)
- **Executed:** `gate_report_8.jld2` structural probe (keys, per-column ECE/MCE, ECE identity check, `training_imsize_provenance` for three bundles)
- **Executed:** Monte-Carlo null-ECE distribution, 2000 replicates at M ∈ {100, 250, 500, 1000, 2000}
- **Executed:** `sha256` comparison proving `artifacts/grid_8.tar.gz` == the `amended_v2` net
- `.planning/phases/07-productionization-conditional-on-go/gate-8x8-amended.md` + `confirmatory-run.log` — the 30.9 min wall-clock datum and the recorded verdicts
- `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md` §6 — the no-iteration protocol and the measured per-pair cost table
- `.planning/CONVENTIONS.md` C-04, C-05 — index-keyed generation's two faces; raw-vs-stripped source assertions
- `.planning/STATE.md:55-81, 84-160` — the Phases 14/15/16 pool-overlap warning; the 2026-07-31 Phase-12 rulings and the four-amendments standard
- `.claude/skills/spike-findings-proteincoloc/SKILL.md` — seed discipline; §6.4 blocks iterating the amended gate

### Secondary (verified against official repositories)

- https://github.com/julia-actions/setup-julia — `@v3`, inputs, named versions
- https://github.com/julia-actions/cache — `@v3`, inputs, depot subdirs cached
- https://github.com/julia-actions/julia-buildpkg — `@v1`, `localregistry`/`git_cli`
- https://raw.githubusercontent.com/julia-actions/julia-runtest/master/action.yml — full input list with defaults
- https://pkgdocs.julialang.org/v1/creating-packages/ — legacy `[extras]`/`[targets]` vs `test/Project.toml` vs workspaces

### Tertiary (single-source, flagged)

- GLMakie headless CI apt package list — from search results referencing MakieOrg/Makie.jl issue #1953
  and the GLMakie README's "GLMakie's CI has no GPU" note. **Not verified by an actual CI run** — this
  is assumption A5 and should be proven by the first CI task.

---

## Metadata

**Confidence breakdown:**

- Machinery contracts (generator signature, `sbc.jl` shapes, `OOD_FAMILIES` consumption, seed
  construction, JLD2 report schema): **HIGH** — every claim read from the source file or executed
  against the artifact this session.
- ECE mechanics, anchors, and the null distribution: **HIGH** — the identity was verified numerically
  to the last float digit against the stored rank table, and the null distribution was measured, not
  assumed.
- Prior bounds and the in-prior/out-of-model typing: **HIGH** for the numbers; **MEDIUM** for the
  `background` typing conflict's resolution, which needs a user ruling (Open Question 2).
- Compute budget: **MEDIUM-HIGH** — the model reconciles the one observed 30.9 min run from
  independent measurements, but the per-rung extrapolation to misspecified generators (A2) and the
  G=16 summary cost (A3) are unmeasured.
- grid_16 blocker: **HIGH** — `recorded = false` read directly off the bundle.
- CI action versions and inputs: **HIGH** — fetched from the official repositories, not memory.
- CI viability end-to-end (GLMakie headless, `Pkg.instantiate` behaviour): **MEDIUM** — reasoned from
  verified repository facts but never executed on a runner, because no CI exists yet.

**Research date:** 2026-08-04
**Valid until:** 2026-09-03 for the in-repo findings (nothing here moves unless the repo moves);
**2026-08-18** for the GitHub Actions version claims (action major tags move faster than that window
suggests — re-verify `setup-julia`/`cache` tags at plan time if more than two weeks have passed).
