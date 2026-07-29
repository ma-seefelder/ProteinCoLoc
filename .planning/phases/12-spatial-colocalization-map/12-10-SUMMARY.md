---
phase: 12-spatial-colocalization-map
plan: 10
status: complete
subsystem: reported validation (the per-region SIM-02 claim)
tags: [sim02, marginal-preservation, pre-registered-bar, reported-artifact, mutation-tested, atom-mass]
requires: ["12-01", "12-03", "12-07", "12-08"]
provides:
  - "spike/validation/run_p12_sim02.jl — the REPORTED per-region SIM-02 runner: M = 20 000 induced-mu W1 vs MU_PRIOR per region, both spatial arms across every P12_R1_LADDER rung plus the D-10 ablation, gated on the MAX over the 64 regions"
  - "spike/validation/p12_sim02_report.jld2 — TRACKED artifact: w1_per_region / w1_max / w1_mean, the full-support companion, corner regions, per-cell field sd and atom_mass_per_region, for all 11 arm x rung configurations"
  - "spike/test/test_p12_prior.jl testset 17 — holds the committed artifact to the same Tier-1 bar and asserts the ladder sweep, so the claim cannot silently decay"
affects: []
tech-stack:
  added: []
  patterns:
    - "the pre-registered bar is READ from Tier 1; the runner carries no numeric tolerance literal of its own and writes `tolerance_source` into the artifact"
    - "the reported artifact is TRACKED and re-asserted by the suite, so a 20 000-draw claim survives without re-running it"
    - "an analytic reference distribution instead of a second finite sample, so the reported statistic carries no avoidable Monte-Carlo noise"
    - "mutation-verified gate: the rescale was deleted in memory and the runner's own @assert observed to fire, with the committed artifact proved untouched"
key-files:
  created:
    - spike/validation/run_p12_sim02.jl
    - spike/validation/p12_sim02_report.jld2
  modified:
    - spike/test/test_p12_prior.jl
decisions:
  - "PASS: max per-region W1 = 0.006750 (scoped, gp@r1=0.25) and 0.007615 (full support, car@r1=0.5) against the Tier-1 bar 0.10, at M = 20 000, on all 11 arm x rung configurations. ~15x headroom."
  - "The tolerance P12_SIM02_W1_TOL_PERREGION = 0.10 was fixed in Tier 1 BEFORE the field simulator existed, is identical to prior.jl's SIM02_W1_TOL, and was NOT derived from this run. p12_consts.jl is byte-unchanged by this plan."
  - "ADDITION, not a threshold change: the same estimator is also reported over the FULL MU_PRIOR support and held to the SAME pre-registered bar. Measured on the rescale-deleted mutant at M = 20 000, the pre-registered SCOPED statistic PASSES a broken copula at three of the five CAR rungs (0.0775 / 0.0394 / 0.0905 at r1 = 0.50 / 0.75 / 0.95) while the full-support statistic fails at two of them (0.1138 / 0.3241)."
  - "The per-region atom mass is REPORTED and gates nothing: 0.067641 mean over configurations, 0.07355 max over regions. It reproduces the analytic MU_PRIOR mass outside [GHAT_MU_MIN, GHAT_MU_MAX] = 0.0676354 to 5.3e-6 — an independent confirmation of the marginal claim."
metrics:
  duration: ~80 min
  completed: 2026-07-29
  tasks: 2
  commits: 3
  files: 3
---

# Phase 12 Plan 10: The Per-Region SIM-02 Claim — Summary

**The claim is established, at the pre-registered sample size, against a bar frozen before the
field simulator existed.** Under D-05's copula every one of the 64 regions carries the `MU_PRIOR`
marginal: the max per-region Wasserstein-1 distance over **all 11 arm × rung configurations** is
**0.006750**, against `P12_SIM02_W1_TOL_PERREGION = 0.10` — roughly **15× headroom**. The per-region
`ghat` atom mass is on record at **0.067641** as a number that gates nothing, and the suite now
enforces both without anyone re-running a 20 000-draw job.

## What Was Built

**Task 1 — `spike/validation/run_p12_sim02.jl`** (commit `61c293a`, 351 lines, new file)

Prior-only: no image simulation, no `simulate_pair`, no network, no `src/` file, no package added,
no constant mutated. `using JLD2, Dates, Statistics` plus two guarded repo-local includes;
`Distributions` and `MU_PRIOR` arrive transitively through `p12_prior.jl` → `prior.jl`, so this file
*cannot* fork π(θ) by restating it. All randomness rides one stream on one counter,
`p12_rng(P12_SIM02_COUNTER)`, consumed sequentially — the run is deterministic and independent of
thread count.

The scoped reference is **copied, not re-derived**: `Truncated(Cauchy(0.0, 0.3), GHAT_MU_MIN,
GHAT_MU_MAX)`, verbatim from `test_p12_prior.jl:77`, which in turn implements `prior.jl:42-46`'s D-16
scoping. Because a restatement is a fork risk, the file asserts at load that its base distribution is
`MU_PRIOR.untruncated` — a Tier-1 or `prior.jl` change that moved the μ-prior would fail this runner
at load rather than leave a reported W1 measured against a stale target.

The W1 estimator is a deliberate second copy of `calibration.jl:137-143` (that file's `_w1` is not
factored out into anything includable; reaching it means including `calibration.jl` whole, which
re-binds its own `MU_PRIOR` and `SIM02_W1_TOL` and pulls `contract.jl` → `src/` and `forward.jl` into
a prior-only runner). The reference is **analytic**, not a second finite sample — 12-07's precedent —
so the reported number carries no avoidable Monte-Carlo noise.

**Task 2 — the run, the tracked artifact, and testset 17** (commit `6775b68`)

`spike/validation/p12_sim02_report.jld2` is **tracked** (`git check-ignore` exits 1), matching the
`p11_recovery_report.jld2` precedent. `spike/test/test_p12_prior.jl` gains **one** testset (17),
appended below 12-07's eleven and 12-08's five, none of which was touched.

## The measurement — re-derived by loading the `.jld2`, not copied from the console

`M = 20000`, `counter = 1`, `tolerance = 0.1`,
`tolerance_source = "Tier 1, p12_consts.jl, fixed before the field simulator existed"`,
`schema_version = 1`, `generated = 2026-07-29T18:26:39.704Z`, `elapsed_min = 0.149`.

| config | w1_max | w1_mean | w1full_max | w1full_mean | sd_min | sd_max | atom_mean | atom_max |
|---|---|---|---|---|---|---|---|---|
| car@r1=0.05 | 0.005858 | 0.002729 | 0.006992 | 0.003130 | 0.9856 | 1.0124 | 0.067847 | 0.07100 |
| car@r1=0.25 | 0.005559 | 0.002767 | 0.006675 | 0.003237 | 0.9904 | 1.0127 | 0.067373 | 0.07320 |
| car@r1=0.50 | 0.005050 | 0.002744 | **0.007615** | 0.003269 | 0.9897 | 1.0114 | 0.067659 | 0.07095 |
| car@r1=0.75 | 0.004271 | 0.002472 | 0.005030 | 0.002819 | 0.9906 | 1.0110 | 0.067852 | 0.07140 |
| car@r1=0.95 | 0.004051 | 0.002294 | 0.003431 | 0.002477 | 0.9985 | 1.0084 | 0.069157 | 0.07085 |
| gp@r1=0.05 | 0.006432 | 0.003027 | 0.006575 | 0.003499 | 0.9895 | 1.0176 | 0.067782 | 0.07355 |
| gp@r1=0.25 | **0.006750** | 0.002870 | 0.006520 | 0.003456 | 0.9845 | 1.0088 | 0.067389 | 0.07125 |
| gp@r1=0.50 | 0.005210 | 0.002651 | 0.005738 | 0.003036 | 0.9905 | 1.0105 | 0.067706 | 0.07260 |
| gp@r1=0.75 | 0.004415 | 0.002557 | 0.006988 | 0.003054 | 0.9908 | 1.0114 | 0.067819 | 0.07145 |
| gp@r1=0.95 | 0.004571 | 0.002557 | 0.004006 | 0.002802 | 0.9896 | 1.0044 | 0.066248 | 0.06945 |
| none@r1=0.0 | 0.005158 | 0.002560 | 0.006482 | 0.003249 | 0.9899 | 1.0079 | 0.067216 | 0.07005 |

**Verdict: PASS on every arm × rung, on both scopings, against the unchanged pre-registered bar.**

- **GATED** — max over all configurations of the per-region W1: **0.0067503** (scoped, `gp@r1=0.25`)
  and **0.0076153** (full support, `car@r1=0.50`), against `P12_SIM02_W1_TOL_PERREGION` = **0.10**.
- Mean of the per-configuration means: **0.0026570** (scoped). Reported; gates nothing.
- Per-cell field sd across all 11 configurations × 64 cells: **[0.98448, 1.01759]** — the
  `./ sqrt.(diag(Σ))` rescale, made visible.
- Corner regions `(1,1) / (1,8) / (8,1) / (8,8)`, at `car@r1=0.50`: 0.002645 / 0.002185 / 0.004281 /
  0.002107. Across all configurations the corner values run 0.001511 – 0.006280 — **no corner
  anomaly**, and the per-config worst region is scattered (55, 61, 52, 42, 39, 6, 42, 17, 10, 9, 3),
  i.e. Monte-Carlo noise rather than a spatial pattern.
- **Atom mass — REPORTED, GATES NOTHING**: mean over configurations **0.067641**, per-configuration
  range 0.066248 – 0.069157, per-region max **0.07355**.

### Does it match the research forecast? Yes, and one number matches the companion statistic better

`12-RESEARCH.md` Pattern 2 measured its row at **CAR α = 0.95**, which is induced r₁ ≈ 0.44
(`p12_lattice.jl:225-230`), so the comparable rung is **`car@r1=0.50`** (α = 0.947559, 12-07).

| Quantity | Research forecast | Measured here | Reading |
|---|---|---|---|
| per-cell sd after rescale | [0.9852, 1.0115] | [0.9897, 1.0114] at `car@r1=0.50` | **confirmed** |
| per-region W1, **max** | **0.00732** | 0.005050 scoped / 0.007615 full support | **confirmed** |
| per-region W1, mean | 0.00344 | 0.002744 scoped / 0.003269 full support | **confirmed** |
| per-region atom mass | **0.0679** | **0.067659** | **confirmed** |
| global 1-D reference atom mass | 0.0678 | analytic 0.0676354 (see below) | **confirmed** |

Two honest observations, stated plainly rather than framed away:

1. **The forecast sits between the two scopings.** At the comparable rung the scoped max (0.005050)
   is 31 % *below* the forecast and the full-support max (0.007615) is 4 % *above* it; the same
   holds for the means (0.002744 / 0.003269 against 0.00344). The likeliest explanation is that the
   research measured over the full support, or against a second finite sample rather than an
   analytic reference — both of which add noise the scoped analytic estimator here removes. This is
   a note about how the forecast was computed, not a claim about which number is right: **both**
   measured statistics clear the pre-registered bar by more than an order of magnitude.
2. **The atom mass is not merely close to the forecast — it is analytically predicted.** Under the
   copula the atoms are exactly the `MU_PRIOR` mass falling outside `[GHAT_MU_MIN, GHAT_MU_MAX]`,
   which is `cdf(MU_PRIOR, GHAT_MU_MIN) + 1 - cdf(MU_PRIOR, GHAT_MU_MAX)` = **0.06763542**. Measured
   mean: **0.06764070**. Agreement to **5.3e-6**. Because the atom mass is a *functional of the
   marginal*, this is an independent confirmation of the per-region marginal claim that does not go
   through the W1 estimator at all. It changes nothing and gates nothing; it is an input to 12-18's
   randomized-rank budget and to the phase report's named limits, as pre-registered.

### The tolerance was not reverse-engineered, in as many words

`P12_SIM02_W1_TOL_PERREGION = 0.10` was declared in the **Tier-1 pre-registration**
(`p12_consts.jl` §5), is **identical to `prior.jl:46`'s `SIM02_W1_TOL`**, and was fixed **before the
Phase-12 field simulator existed** — before `p12_lattice.jl`, before `p12_prior.jl`, before this
runner. **It was not derived from this run, and no threshold was moved, widened, added or
re-derived by this plan.** The runner carries **no numeric tolerance literal of its own** (verified:
after stripping comment lines the only float literals in the file are the copied
`Cauchy(0.0, 0.3)` reference parameters and the `0.5` quantile-grid offset of the W1 estimator), and
`spike/validation/p12_consts.jl` is **byte-unchanged** (`git diff --quiet HEAD` clean).

## Verification Results — what was actually run and observed

| Check | Command | Observed |
|---|---|---|
| Task 1 verify (verbatim from the plan) | the plan's `julia -e` one-liner | **PASS**, exit 0 — printed `sim02-load-ok` |
| the reported run | `julia --project=spike -t auto spike/validation/run_p12_sim02.jl` | **PASS**, exit 0 — 11 configurations, all `PASS`, 19.6 s wall (`elapsed_min` 0.149) |
| Task 2 verify (verbatim from the plan) | the plan's `JLD2.load` one-liner | **PASS**, exit 0 — printed `0.0067503000646160684 0.1 20000` |
| testset 17 green in isolation | `julia --project=spike -e 'using Test; include("spike/test/test_p12_prior.jl")'` | **PASS** — 17 testsets, **146 assertions**, 0 failed, 0 errored, 22.2 s |
| 12-07's and 12-08's assertions still green | same run | **PASS** — 5/13/4/9/4/4/5/3/5/10/1 (12-07, **63**) and 11/20/5/10/15 (12-08, **61**) exactly as their summaries record; testset 17 adds **22** |
| testset 17 green under the suite | `julia --project=spike spike/test/runtests.jl` | **PASS** — 22/22 after the `P12-SUITE-RAN: test_p12_prior.jl` marker |
| every earlier Phase-12 file still green under the suite | same run | **PASS** — `test_p12_consts.jl`, `test_p12_lattice.jl`, `test_p12_architecture.jl` … all zero failures |
| decoupling (12-05's live self-scanning test) | same run | **PASS 26/26** — `src/` byte-unchanged, spike env byte-frozen, the frozen sixteen dependencies, sealed holdout not consumed, Phase-11 pool read-only |
| artifact is tracked | `git check-ignore -q spike/validation/p12_sim02_report.jld2` | **exit 1** (not ignored) — matches the `p11_recovery_report.jld2` precedent |
| `p12_consts.jl` byte-unchanged | `git diff --quiet HEAD -- spike/validation/p12_consts.jl` | **PASS** — this plan appends no Tier-2 constant and edits no Tier-1 line |
| frozen surfaces byte-unchanged | `git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml src corpus spike/simulator` | **PASS** |
| no pool written | `git status --porcelain -- spike/data/cache` | **empty** — prior-only, as designed |
| `spike/data/cache/p11` intact | `du -sh` | **PASS** — 54 MB |
| `.planning/STATE.md`, `.planning/ROADMAP.md` untouched | `git diff --quiet HEAD --` | **PASS** |

**Not run:** `julia --project=. -e 'using Pkg; Pkg.test()'`. No pre-plan baseline exists (12-01,
12-03 and 12-07 recorded the same gap), so "unchanged" is not something this run can assert. `src/`
is byte-unchanged by `git diff --quiet` and by 12-05's live decoupling testset, and this plan adds
nothing to the main package and loads nothing from it. Recorded rather than implied.

### The mutation proofs — the gate was watched failing, and the artifact was proved untouched

Both proofs were applied **in memory** (redefining `field_sampler` in `Main` after the include —
12-07's idiom) or on a **byte-exact backup**; nothing on disk was left changed.

**(1) The runner's own `@assert`.** With the `./ sqrt.(diag(Σ))` rescale deleted, the REAL `main()`
was called at the REAL `M = 20 000`. It aborted with the runner's own message:

> `AssertionError: SIM-02 per-region FAILURE (scoped): arm = car, r1 = 0.05, region 55 reads
> W1 = 0.12203516427760797 against the Tier-1 bar P12_SIM02_W1_TOL_PERREGION = 0.1.`

The assertion loop sits **before** the `jldsave`, so the committed artifact's mtime was verified
**unchanged** by the mutant run. The mutant's full table, at M = 20 000:

| CAR rung | scoped w1_max | verdict on the pre-registered bar | full-support w1_max | verdict |
|---|---|---|---|---|
| 0.05 | **0.122035** | **FAIL** | **0.161973** | **FAIL** |
| 0.25 | **0.108917** | **FAIL** | **0.148447** | **FAIL** |
| 0.50 | 0.077527 | *passes a broken copula* | **0.113751** | **FAIL** |
| 0.75 | 0.039385 | *passes a broken copula* | 0.086825 | *passes* |
| 0.95 | 0.090524 | *passes a broken copula* | **0.324124** | **FAIL** |

Three readings follow, and all three are measurements rather than opinions:
- **Sweeping the ladder is what makes this runner's assertion able to fail.** A single-rung run at
  r₁ = 0.50, 0.75 or 0.95 would have certified a demonstrably broken copula. This is 12-07's
  recorded advice, confirmed at the reported M.
- **The full-support companion catches two of the three rungs the scoped statistic misses**, for
  the reason 12-07 recorded: the D-16 scoping removes exactly the tails into which a compressed
  marginal moves its mass.
- **Neither statistic fires at r₁ = 0.75** — a named limit, not a defect: the pre-registered bar was
  sized as a "distributions are close" criterion for a *scalar* claim with design headroom, never as
  a mutation detector. It is stated here so no later document reads "the SIM-02 gate is green" as
  "no marginal defect of any size could be present at any rung".
- The GP arm and the `:none` ablation are **unaffected by the mutation** (identical to six decimal
  places), because the Matérn-3/2 kernel and the identity both carry a unit diagonal, so
  `./ sqrt.(diag(Σ))` is a no-op there. The rescale is a **CAR-arm** property; a GP-only screen
  cannot see it.

**(2) Testset 17.** A tampered copy of the artifact (`w1_max[3] → 0.5`, `config_labels` truncated to
one entry) was written over the real path after a byte-exact backup. The testset went **red: 10
passed, 12 failed** — the breached max *and* the truncated sweep both caught. The original was then
restored from the backup and the SHA-256 verified identical
(`d41d5da6fbd7744c58850759132560903cf16bca39b6d306919088adfbff8819` before and after), and the
testset re-run green at 22/22. No `git checkout`, `stash`, `clean`, `reset` or `amend` was used.

## Deviations from Plan

### 1. [Rule 2 — missing critical verification] A full-support companion statistic, held to the SAME pre-registered bar

- **Found during:** Task 1, before the reported run, while probing the mutant on the fixture stream.
- **Issue:** the plan gates on the W1 **scoped** to `[GHAT_MU_MIN, GHAT_MU_MAX]`. 12-07's recorded
  Assumption Drift states that this scoping *damps the very failure mode the claim is about*,
  because a compressed marginal moves most of its mass into the tails the scoping removes. Measured
  (above), that is exactly what happens: the scoped statistic passes a rescale-deleted copula at
  three of five CAR rungs.
- **Fix:** the runner computes the **same estimator over the full `MU_PRIOR` support** and asserts it
  against the **same** `P12_SIM02_W1_TOL_PERREGION`. This is an **addition, not a threshold change**:
  no new constant exists, no bar was moved, relaxed or invented, and the pre-registered scoped
  statistic is still computed, still reported and still gated exactly as specified. `w1_max` in the
  artifact remains the scoped quantity the plan's verify command reads.
- **Files:** `spike/validation/run_p12_sim02.jl`, `spike/test/test_p12_prior.jl`.
- **Commits:** `61c293a`, `6775b68`.

### 2. [Rule 2] Testset 17 asserts the ladder sweep, which the plan's assertion list does not name

- **Found during:** Task 2, from the mutation table above.
- **Issue:** the plan's seven listed assertions would all pass on an artifact regenerated from a
  single rung — and such an artifact would certify a broken copula at three of five CAR rungs.
- **Fix:** testset 17 asserts that `config_labels` contains every `P12_R1_LADDER` rung for both
  spatial arms plus the ablation, and that there are exactly `2 × 5 + 1 = 11` of them. Also added:
  the full-support max under the same bar, `all(length == P12_G^2)` over the per-region vectors
  (rather than only `first(values(...))`), and that `tolerance_source` records the bar as
  pre-existing the field simulator. Nothing was removed: all seven of the plan's assertions are
  present verbatim in substance.

### 3. Atom counting is by exact equality against the frozen knots, not a `±0.99` literal

- **Why:** the plan says "fraction of `rho_field` entries at `±0.99`". `ghat` *returns*
  `GHAT_RHO_KNOTS[1]` / `[end]` unchanged when μ falls outside the realized range, so the atoms are
  those `Float64`s exactly. Reading the frozen knots keeps the count from drifting from the map and
  keeps a hand-written numeric literal out of a file whose acceptance criterion forbids one. The
  measured value (0.067641) matches the analytic prediction to 5.3e-6, so nothing is being missed.

### 4. One descriptive inaccuracy in the plan, reported rather than silently corrected

- **`read_first` describes `calibration.jl`'s W1 machinery as something to "reuse … if one exists".**
  `_w1` exists (`calibration.jl:137-143`) but is a file-scope function in a script that also binds
  its **own** `MU_PRIOR` and `SIM02_W1_TOL` and includes `contract.jl` → `src/` and `forward.jl`.
  It is guarded by `abspath(PROGRAM_FILE) == @__FILE__`, so including it would *not* re-run the
  calibration — but it would fork π(θ)'s bindings and pull the whole forward model into a
  prior-only runner. The plan's fallback ("write one local helper and say in the comment that it
  duplicates a scoped computation deliberately, with the scoping copied, not re-derived") is what
  shipped, with the reason stated accurately in the file rather than as "it would re-run".

## Assumption Drift (advisory)

**The plan assumed the pre-registered bar, swept across the ladder, is what makes this runner able
to fail. Measured, the ladder sweep is necessary but the scoped statistic alone is not sufficient.**

- **Planned:** *"This runner is one of the two that may fail the phase"* — i.e. `maximum(w1) ≤
  P12_SIM02_W1_TOL_PERREGION`, scoped, across the five rungs, is the protective assertion.
- **Actual:** across the five CAR rungs the scoped statistic on a rescale-deleted copula reads
  0.1220 / 0.1089 / 0.0775 / 0.0394 / 0.0905 — it fires at two rungs and passes at three. The full
  ladder sweep is therefore load-bearing (a single-rung runner would be inert at three of five
  rungs), and the full-support companion is what covers two of the three misses. At r₁ = 0.75
  nothing fires.
- **Why:** the bar is `SIM02_W1_TOL`'s value, sized as a "distributions are close" criterion for the
  *scalar* induced-μ claim with 13× design headroom, and rightly not reverse-engineered from any
  measurement. It was never sized as a mutation detector.
- **Why this is advisory:** nothing gated moved. The pre-registered bar is unchanged and is asserted
  exactly as specified; the reported claim passes it with ~15× headroom on both scopings. Recorded
  so no later document reads "SIM-02 is green" as "no marginal defect of any magnitude can exist".

## Deferred Issues

**`runtests.jl` exits 1 with 3 failures in `test_p13_real.jl` ("P13 real alpha ladder", lines 226 and
292) — Phase 13's, not this plan's.** The Phase-12 abort point has now moved a **fourth** time
(Phase-4 `SPEEDUP_GATE` → `test_sbc.jl` → `test_p13_consts.jl`, which 12-07 recorded, → now
`test_p13_real.jl`), exactly as the standing briefing warns. `test_p13_real.jl` includes only
`spike/p13/real_images.jl`; nothing in Phase 13 references `run_p12_sim02.jl`,
`p12_sim02_report.jld2` or `test_p12_prior.jl`, and the Phase-13 files are owned by a live
concurrent executor on this branch. Every Phase-12 testset runs and reports green before the abort
(`P12-SUITE-RAN` markers for consts, lattice, prior, architecture, datagen, train, result, sbc,
coverage, decoupling). Do not read a red suite as a Phase-12 failure.

DEF-12-03 (the `:P13_DEV_SEED` include-guard poisoning) is unchanged and remains the user's open
ruling; it is not touched here.

## Known Stubs

None. The runner, the artifact and the testset are complete; nothing in this plan is deferred to a
later one.

## Threat Flags

None. No new network, auth or file-access surface: the runner samples in memory, reads only the
frozen pre-registration, the frozen `ghat` and the committed copula, and writes exactly one
artifact inside `spike/validation/` by the shared atomic `.tmp` → reopen-and-check → `mv(force)`
discipline. The plan's five dispositions are mitigated and asserted: **T-12-36** (the runner carries
no tolerance of its own; `tolerance_source` is written into the artifact; `p12_consts.jl` asserted
byte-unchanged), **T-12-37** (the gate is `maximum` over 64 regions, the corners are stored by name,
and no `mean(` appears in either gating assertion), **T-12-08** (`gates nothing` in the caption,
asserted by testset 17; the atom mass changed no threshold), **T-12-38** (the artifact is tracked
and re-asserted by the suite — proved by tampering it and watching the suite go red), **T-12-07**
(three pathspec-scoped commits; no index.lock race occurred despite the concurrent Phase-13
executor).

## For the Next Plan

- **The per-region SIM-02 claim is settled evidence.** `12-11`, `12-15`, `12-16` and the phase report
  may cite `p12_sim02_report.jld2` directly; the numbers above are re-derivable by
  `JLD2.load` and are enforced by `test_p12_prior.jl` testset 17.
- **Quote the atom mass as 0.067641 (measured) / 0.0676354 (analytic).** 12-18's randomized-rank
  budget sizes against a per-region atom rate of ≈ 6.76 %, occurring on **all 64 regions per draw**.
  It is materially *at* the forecast, not above it, so no budget revision is implied.
- **Do not quote "0.658 edge / 0.992 interior"** (12-07's standing correction) and do not quote the
  research's per-region W1 forecast as if it were this run's measurement — the two differ by scoping
  (see the forecast table).
- **If a later plan re-generates this artifact, it must keep the 11-configuration sweep.** Testset 17
  asserts it, and the mutation table above is why.
- `spike/data/cache/` is untouched by this plan: the runner is prior-only and writes no pool.

## Self-Check: PASSED

- `spike/validation/run_p12_sim02.jl` — **FOUND** (351 lines).
- `spike/validation/p12_sim02_report.jld2` — **FOUND**, tracked, SHA-256
  `d41d5da6fbd7744c58850759132560903cf16bca39b6d306919088adfbff8819`.
- `spike/test/test_p12_prior.jl` — **FOUND** (684 lines), 17 testsets / 146 assertions green.
- Commits **FOUND** in `git log`: `61c293a` (the runner), `6775b68` (the artifact and testset 17).
- `spike/Project.toml`, `spike/Manifest.toml`, `src/`, `corpus/`, `spike/simulator/`,
  `spike/validation/p12_consts.jl` byte-unchanged; `.planning/STATE.md` and `.planning/ROADMAP.md`
  untouched; `spike/data/cache/p11` intact at 54 MB; `spike/data/cache` free of new files.
- All commits used explicit pathspecs; no `git stash`, `rebase`, `amend`, `reset`, `clean` or
  `checkout --` was run at any point.
