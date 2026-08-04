# Phase 15: Calibration Operating Envelope and CI Gate - Context

**Gathered:** 2026-08-04 (discuss --analyze --all)
**Status:** Ready for planning

<domain>
## Phase Boundary

Map the tool's domain of applicability by adversarially sweeping nuisances until coverage breaks,
prove the OOD flag fires before it does, and lock calibration into CI as a regression gate.

**Six scouting findings, verified on disk. Four change what can be built.**

1. **There is no CI.** No `.github/` directory, no active git hooks (`.git/hooks/` has samples only),
   no `Makefile`. SC3 says "fails CI on calibration drift" — the CI it would fail does not exist.
   This phase builds one (D-08).

2. **SC1's four named axes do not match the four that exist.** `test/gate/misspec.jl:181-183` ships
   `OOD_FAMILIES = (texture, noise, optics, background)`. Mapping onto SC1: PSF→`optics` ✓,
   autofluorescence→`background` ✓, **spillover → no generator**, **registration → no generator**.
   Two named axes must be built; two shipped families are unnamed in SC1. Resolved by D-02.

3. **Three of SC1's four axes are INSIDE the training prior, not misspecifications.**
   `src/amortized/simulator.jl:92,138-142` — `spillover ~ Uniform(0, 0.2)`, `autofluorescence`, and
   `shift_dx`/`shift_dy ~ Uniform(-1, 1)` are all *inferred nuisances the net was trained on*. Only
   `σ_psf = 1.3` (`simulator.jl:179`) is fixed. Sweeping an in-prior axis is **extrapolation past
   prior support**; sweeping PSF is **genuine out-of-model misspecification**. These are two
   different experiments and D-02/D-03 keep them typed apart.

4. **Full SBC does not fit in any CI.** `SBC_M = 2000` × `SBC_L = 999`
   (`test/gate/gate_consts_8_v2.jl:68`, `spike/validation/consts.jl:44-45`) over the
   `SBC_IMSIZE_SET` mixture of 512²/1024²/1376×1028/2048² images. `SBC_FIX_M = 8`
   (`gate_consts_8_v2.jl:169`) exists as a smoke fixture. That gap is the whole SC3 design problem,
   resolved by the two-tier D-08.

5. **The C-04 pool-overlap hazard fires in this phase by name.** `.planning/STATE.md` flags Phases
   14/15/16 explicitly: *"a later phase generates 'a fresh evaluation pool' at
   `arm = P12_CHOSEN_PRIOR` — under C-04 its first `min(n, 50000)` samples ARE the production
   training set. Any coverage, SBC or BF number so obtained is measured on training data and is not
   a calibration result."* Resolved structurally by D-11.

6. **Both dependency phases closed NEGATIVE.** Phase 11 never reshipped a registration-aware net
   (`11-CLOSURE.md`); Phase 12 delivered no calibrated per-region uncertainty (`12-VERIFICATION.md`,
   goal not achieved). The envelope is therefore mapped for the **shipped** `amended_v2/grid_8` net,
   not for a Phase-11/12 successor (D-10).

**A HARD CONSTRAINT ON SEEDS, verified on disk, that binds planning.** `test/gate/run_gate.jl:47-52`
records that `gate_consts_8_v2.jl` *"binds it to exactly ONE run on the fresh `PROD_SEED_V2[8]`
(protocol §6.3), so consuming it must be an explicit, auditable invocation … and never a side
effect."* **That single run is spent** — it produced `artifacts/amended_v2/grid_8/gate_report_8.jld2`.
Phase 15 must NOT draw on `PROD_SEED_V2[8]`. See D-11.

**Out of scope:** retraining or retuning any net; per-region/per-tile calibration (blocked on Phase
12's NO); promoting the Phase-13 three-way net or the Phase-14 decision layer into `src/`; the blind
external evaluation itself (Phase 16); adding a new OOD detector channel (see Deferred).

</domain>

<decisions>
## Implementation Decisions

### Lane and integration surface

- **D-01: `src`-side, in `test/gate/` — NOT the spike research lane.**
  This deliberately breaks with Phases 13 and 14, which built in `spike/` (their D-01s). Three
  reasons, all structural: SC3's gate must run against the **shipped package** to be a regression
  gate at all; the envelope is a claim about the **shipped net** (D-10); and the three machines this
  phase needs — `test/gate/misspec.jl` (families), `test/gate/sbc.jl` (ranks/ECE/coverage),
  `test/gate/harness.jl` (seeded draw→simulate→infer) — **already live in `test/gate/`**.
  **Consequence for planning:** new sweep generators extend `misspec.jl`'s existing signature; the
  runner extends `run_gate.jl`'s CLI rather than forking a parallel one.

### Sweep design

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

### Break criterion

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

### OOD-before-break (SC2)

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

### CI regression gate (SC3)

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

### Nets and draw provenance

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

### Pre-registration and failure discipline

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

### AMENDMENTS — appended 2026-08-04 after Phase-15 research (`15-RESEARCH.md`, `919be48`)

**Additive block. Nothing above is edited.** Four findings from research were ruled by the user
before planning. Each is recorded here as a correction to the decision it supersedes, with the
evidence that forced it.

- **D-10a (SUPERSEDES D-10's contrast half; D-10's reported-net half stands).** **The `grid_16`
  contrast is DROPPED. The reported envelope is `amended_v2/grid_8` ONLY.**
  Verified directly: `artifacts/grid_16/npe_16.jld2`'s `meta` is `(:grid, :n_pairs, :use_gpu)` —
  `training_imsize_provenance` is **absent entirely**, not merely unrecorded, so
  `SBC_REQUIRE_IMSIZE_PROVENANCE = true` rejects it and `run_gate` returns `:provenance_mismatch`
  and runs no arm. It is a v1-era net (256², unbounded θ); bypassing the check would confound
  **grid × imsize-regime × θ-space** in a single comparison, so the resulting number could not
  answer the grid-dependence question the contrast exists to answer.
  **Carry as a named limit, not a silent scope cut:** *grid-dependence of the operating envelope is
  UNTESTED; the contrast was attempted and blocked by missing provenance on the only available
  higher-dimension net.* Halves the compute budget.

- **D-02a (EXTENDS D-02; the six axes and the mechanism typing are unchanged).** **A SEVENTH axis is
  added: a pure-offset autofluorescence generator**, whose `level` maps directly onto the
  simulator's autofluorescence offset parameter so its prior boundary is markable per D-03.
  The existing `background` family is **reclassified as out-of-model ONLY** and keeps no boundary
  rung. Evidence: `background`'s `level = 1` bleed is already **0.4 against a prior max of 0.1**, and
  it composes a radial vignette and gradient the simulator has **no parameter for** — so it cannot
  express "at the prior edge" in any units, and D-03 was unsatisfiable on it as written.
  **D-03 now holds on every in-prior axis without exception**, which is why this route was chosen
  over amending D-03 to "where markable".

- **D-04a (SUPERSEDES D-04's ANCHOR, not its statistic).** **The ECE break threshold is anchored at
  RUNG 0 measured at the sweep's own M on the Phase-15 stream — NOT at the shipped gate report's
  ECE.** The margin above it is **`q95 − E[ECE₀]`** from the measured null.
  **This is arithmetic, not preference.** Read from `artifacts/amended_v2/grid_8/gate_report_8.jld2`:
  ρ_true ECE = **0.03713**, Δρ ECE = **0.00429**. The null distribution of ECE for a *perfectly
  calibrated* net at M = 2000 has **E[ECE₀] ≈ 0.00751** (closed form ≈ `0.336/√M`). **Δρ's shipped
  anchor sits BELOW its own null mean** — anchoring there would declare a break on ~90 % of
  perfectly-calibrated rungs. That is D-04's own failure mode reappearing inside the anchor. Anchoring
  at rung 0 at the sweep's M makes the noise floor **common-mode** and removes it.
  **D-04's substance is unchanged:** ECE is still the break statistic, band and p-value are still
  reported and still non-gating.
  **`MCE` MUST NOT BE USED anywhere in this phase** — it is pinned at exactly **0.99 on all eight
  columns** by empty reliability bins. It is an artifact, not a measurement.
  **q95 was chosen over q99 knowing the cost:** ~5 % false-break rate ⇒ ≈1.5 expected false breaks
  across ~30 axis×rung tests. **The map must therefore be read PER AXIS, not per rung** — state this
  wherever the map is reported.

- **D-08a (REFINES D-08; two-tier structure unchanged).** **CI installs via `Pkg.instantiate()` from
  the committed `Manifest.toml`, not by resolving.** `[compat]` is left **untouched**.
  Reason: `test/runtests.jl`'s co-resolution hard gate asserts `NeuralEstimators == "0.2.1"` and
  `Flux == "0.16.10"` as **exact strings**, while `[compat]` permits a wider range — so a CI resolve
  could satisfy `[compat]` and still fail the gate. Instantiating from the manifest makes the pinned
  versions what CI actually gets.
  Not chosen: hard-pinning `[compat]` to `"=0.2.1"`/`"=0.16.10"`, which would block downstream
  co-installation of an AGPL package others may depend on.

- **Compute budget adopted from research (not a user ruling; recorded so the number is traceable):**
  **M = 1000 per rung, 5 rungs per axis**, ≈ **7.3 h** for the full `grid_8` map. The cost model
  reconciles the one observed gate run to within 5 % (predicted 29.4 min vs 30.9 min observed).
  M = 2000 costs 2× for a 1.27× resolution gain. **`JULIA_NUM_THREADS` is a budget factor, not an
  optimisation** — single-threaded is ~13× slower and must be set in every CI job and run script.

- **Implementation constraint found by research and binding on planning:** **do NOT append to
  `OOD_FAMILIES`.** Mutating it would silently rewrite the shipped, *passed* OOD arm's `combined_auc`
  and its pass conjunction. Use a local `P15_FAMILIES = merge(OOD_FAMILIES, ...)` passed through
  `gate_ood_roc`'s existing `families =` keyword — zero-touch on the shipped gate.

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

</decisions>

<specifics>
## Specific Ideas

- **The mechanism typing in D-02 is a deliverable, not an annotation.** The domain map should be
  readable as two overlaid stories: where training is insufficient, and where the model runs out.
- **Every bar in `p15_consts.jl` should carry its derivation next to it.** The Stage-1 ceiling
  needed adjudication precisely because it *"carries no derivation anywhere in the pre-registration"*
  (STATE.md). A bar with no derivation is a liability this project has already paid for.
- **A negative result here is a good outcome and should be written that way.** Phases 11 and 12 both
  closed negative-but-useful; an empty or unbounded envelope, or a LATE axis, is publishable and
  should not be framed as a phase failure.
- Follow the project discipline throughout: **measured rather than chosen**, pre-registered before
  results exist, **corrections appended rather than overwritten**.

</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### The machinery this phase extends (all already in `test/gate/`)
- `test/gate/misspec.jl` — `OOD_FAMILIES = (texture, noise, optics, background)` at :181-183;
  generator signature `(rng, θ; imsize, level, G)`; `gate_ood_roc` at :298-393. D-02 adds two
  families here.
- `test/gate/sbc.jl` — `sbc_ranks`, `sbc_uniformity`, `sbc_coverage`, `sbc_calibration`,
  `sbc_traffic_light`. Supplies the D-04 statistics.
- `test/gate/harness.jl` — `draw_simulate_infer`; `gate_global_seed`/`seed_gate_global!` at :75-95;
  **:50-59 is the global-RNG caveat D-11 requires honouring**; `gate_imsize` at :146.
- `test/gate/run_gate.jl` — the per-grid gate CLI (`--grid`/`--sbc`/`--bf`/`--ood`/`--consts`);
  **:47-52 records the one-run binding of `PROD_SEED_V2[8]` that D-11 must not violate**.
- `test/gate/gate_consts_8_v2.jl` — the pre-registration to MIRROR (not extend): `SBC_M` :68,
  `SBC_IMSIZE_SET`/`_WEIGHTS` :152-153, `SBC_FIX_M` :169, forbidden-seed list :301-305, seed redraw
  loop :278. D-12 copies this structure into a new `p15_consts.jl`.

### The net and the simulator under test
- `artifacts/amended_v2/grid_8/gate_report_8.jld2` — the shipped net's gate report; source of the
  D-04 in-prior ECE anchor.
- `src/amortized/simulator.jl` — `SPILLOVER_PRIOR` :92, `SHIFT_PRIOR` :95, θ draw :138-142, prior
  bounds :171-172, **fixed `σ_psf = 1.3` :179**. The in-prior/out-of-model split D-02 rests on.
- `src/amortized/ood.jl` — `OOD_PP_REPS` :60, **`OOD_ID_QUANTILE = 0.95` :61** (the ~5% ID fire rate
  D-07 exists to handle), `id_threshold` :296-302, `roc_auc`, `ood_flag` OR-fusion.
- `test/runtests.jl` — the co-resolution hard gate (NeuralEstimators 0.2.1 / Flux 0.16.10) that any
  CI job must satisfy first.

### Constraints and hazards that bind this phase
- `.planning/CONVENTIONS.md` **C-04** (:246 ff.) — index-keyed generation has two faces; a "fresh
  pool" at the same arm is a superset of the training pool. **Grounds D-11.**
- `.planning/STATE.md` — the "⚠️ FOR PHASES 14/15/16" block naming this phase; the 2026-07-31
  Phase-12 ruling (coverage band tracking training duration; the refusal to tune to a gate); the
  four-amendments observation and its licensing standard.
- `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md` — §6.4 blocks
  iterating the current amended gate. Phase 15 measures; it does not retrain or retune.
- `docs/amortized.md` — the named v2.0 limits, including the ~0.08 SD nuisance marginal drift that
  D-05 declines to re-gate on.

### Why the dependency phases constrain rather than enable
- `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-CLOSURE.md` — closed
  negative; no registration-aware net exists to sweep.
- `.planning/phases/12-spatial-colocalization-map/12-VERIFICATION.md` — goal not achieved; no
  calibrated per-region uncertainty, so the envelope is per-pair only.
- `.planning/phases/14-decision-and-abstention-layer/14-CONTEXT.md` — D-04's per-pair unit and D-06's
  three-valued NOT-CHECKED semantics; the same asymmetry logic recurs in D-06 here.
- `.planning/ROADMAP.md` Phase 15 entry — SC1/SC2/SC3 as written. **SC1 is EXTENDED by D-02 (four
  named axes all swept, two added); SC2's evidential form is FIXED by D-06.** Both must be cited
  alongside any Phase-15 result.

### Project-level discipline
- `./.claude/skills/spike-findings-proteincoloc/SKILL.md` — seed discipline and the rule that gate
  changes are pre-registration matters rather than in-phase decisions.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`test/gate/` is ~2,650 lines of exactly the right machinery** — `misspec.jl` (396),
  `run_gate.jl` (541), `sbc.jl` (496), `harness.jl` (299), plus five `gate_consts_*` files. This
  phase extends it; it does not build a parallel harness. **This is the single biggest reason for
  D-01's src-side lane.**
- `gate_ood_roc` (`misspec.jl:298`) — already computes fused/density/noise AUC and **`fire_rate`**
  per family × level against a TRAIN-only fitted null. `fire_rate` is directly what D-07 needs.
- `OOD_GRID_LEVELS` + the family × level loop (`misspec.jl:352-386`) — the ladder structure already
  exists; D-03 adds prior-boundary marking and two new families to it.
- `run_gate.jl`'s `--consts FILE` override (:44-58) — the clean, auditable way to point the runner at
  `p15_consts.jl` without touching any committed pre-registration.
- `src/amortized/ood.jl` — `ood_flag`, `ood_verdict`, `OODVerdict`, `is_ood`, `id_threshold`,
  `roc_auc`. Consume directly; **do not reimplement and do not retune** (D-14).

### Established Patterns
- **Two-tier constants files** — Tier 1 append-only bars, Tier 2 sentinel-keyed measurement appends
  (`p12_consts.jl`, `spike/p13/consts.jl`). D-12 follows this exactly.
- **Salted-Philox seed derivation with a forbidden-seed redraw loop** — `gate_consts_8_v2.jl:278,
  301-305`. D-11 copies it.
- **Atomic report writes** — `.tmp` → reopen-integrity check → `mv(...; force=true)`, in
  `run_gate.jl`.
- **Executable constraint guards** — this project encodes hard constraints as running tests
  (`spike/test/test_p12_decoupling.jl`), not prose. Phase 15 should guard D-11's seed disjointness
  and D-09's re-bless discipline the same way.
- **Honest `:not_trained` returns** — `run_gate.jl` returns an honest status rather than crashing
  when an artifact is absent (:16-19). New Phase-15 arms should match.

### Integration Points
- `test/gate/misspec.jl` — two new families (`spillover`, `registration`) appended to
  `OOD_FAMILIES`, matching the existing generator signature.
- `test/gate/p15_consts.jl` (new) — the two-tier pre-registration, selected via `--consts`.
- `test/gate/run_gate.jl` — a new sweep arm (e.g. `--envelope`) alongside `--sbc`/`--bf`/`--ood`.
- `.github/workflows/` (new, whole directory) — the two-tier CI of D-08.
- A committed golden artifact + its documented re-bless procedure (D-09).

</code_context>

<deferred>
## Deferred Ideas

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

</deferred>

---

*Phase: 15-calibration-operating-envelope-and-ci-gate*
*Context gathered: 2026-08-04*
