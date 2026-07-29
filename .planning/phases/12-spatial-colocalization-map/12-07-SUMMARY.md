---
phase: 12-spatial-colocalization-map
plan: 07
status: complete
subsystem: spatial prior (the D-05 copula)
tags: [copula, marginal-preservation, generative-order, separable-upsample, offset-grid, mutation-tested]
requires: ["12-01", "12-03"]
provides:
  - "spike/simulator/p12_prior.jl — rho_field copula (z -> Phi -> MU_PRIOR quantile -> ghat elementwise), r1-first field-aware prior draw, derived scalar theta view, separable bilinear upsample with a keyword offset-grid arm, ASVS-V5 field validation"
  - "spike/test/test_p12_prior.jl — eleven testsets, 63 assertions, replacing the 12-01 pending scaffold"
affects: []
tech-stack:
  added: []
  patterns:
    - "the pre-registered bar is READ, never re-derived; where it cannot falsify, a companion falsifier is added rather than the bar moved"
    - "mutation-verified invariants: every marginal/fidelity claim was broken in memory and observed to fail before being committed"
    - "an analytic reference distribution instead of a second finite sample, so a screening statistic carries no avoidable Monte-Carlo noise"
key-files:
  created:
    - spike/simulator/p12_prior.jl
  modified:
    - spike/test/test_p12_prior.jl
    - .planning/phases/12-spatial-colocalization-map/deferred-items.md
decisions:
  - "MEASURED CORRECTION: the CAR marginal-sd spread runs 0.9863 at a CORNER against 0.6534 in the INTERIOR — the reverse of the '0.658 edge / 0.992 interior' pairing recorded in 12-RESEARCH Pattern 1 and repeated in p12_lattice.jl:41-48. The plan's testset 2 was therefore specified on the four cells a missing rescale does NOT break, and would have been a test that cannot fail."
  - "The pre-registered bar P12_SIM02_W1_TOL_PERREGION = 0.10 does NOT falsify a deleted rescale at r1 = 0.50 (measured max W1 0.0789). It DOES at the bottom ladder rung (0.1264). The bar was not touched; the screen was extended to the rung where the pre-registered bar itself bites."
  - "Testsets 1-3 hoist the sampler and call `rho_field` directly rather than `sample_p12_prior` 2 000 times: the object under test is identical and the alternative re-solves the r1 bisection on every draw."
metrics:
  duration: ~75 min
  completed: 2026-07-29
  tasks: 2
  commits: 2
  files: 2
---

# Phase 12 Plan 07: The D-05 Copula, the r₁-First Draw and the Separable Upsample — Summary

A ρ field can now be drawn whose every one of the 64 regions carries exactly the `MU_PRIOR`
marginal, conditional on a correlation length drawn **first**, with the prior arm recorded on the
draw and a pixel-resolution rendering that costs two matmuls and reaches its own half-cell guard
arm through a single keyword — and the marginal, the order and the upsample are each proved by
**mutation**, not asserted.

## What Was Built

**Task 1 — `spike/simulator/p12_prior.jl`** (commit `6c77dfd`, 305 lines, new file)

Flat top-level functions, `#= =#` AGPL header, two guarded includes (`prior.jl` first — it binds
`MU_PRIOR`, `SIM02_W1_TOL` and all seven nuisance priors, which Phase 12 reuses **verbatim** and
never restates; then `p12_lattice.jl`, which transitively brings the Tier-1 pre-registration).

| Function | What it delivers |
|---|---|
| `validate_rho_field(F)` | ASVS V5: rejects non-square, non-finite and out-of-`[-1,1]` fields, naming the offending value; called from inside `rho_field` |
| `rho_field(rng, sampler; G)` | the D-05 copula, returning `(rho, mu, z)` — `z` is kept because it is atom-free and is what 12-09's θ rows and 12-18's clean SBC arm are built from |
| `sample_p12_prior(rng; G, arm, r1)` | r₁ **first**, field conditional on it, then the seven nuisances in `prior.jl`'s own relative order; carries `arm` and `r1` on the returned draw |
| `theta_scalar_view(draw)` | the legacy 8-field θ `simulate_pair` accepts, with `ρ_true` labelled **DERIVED** (R-1) |
| `p12_bilinear_op(G, m; offset)` | the `m×G` operator; `offset = true` **is** the D-06 half-cell guard arm — one keyword, not a second code path |
| `p12_upsample(F, imsize; offset)` | `A * F * B'` |
| `p12_upsample_ops(G, imsize; offset)` | the two operators, so a datagen loop can hoist them out of the inner loop |

`P12_R1_PRIOR` is **read** through the guarded includes and is deliberately not re-declared: a
local `const` with the same value would be benign in Julia and therefore *worse* — it silently
defeats the guard. Verified: after stripping comment lines the file contains no `P12_R1_PRIOR =`,
no `MU_PRIOR =` and no `*_PRIOR =` assignment.

**Task 2 — `spike/test/test_p12_prior.jl`** (commit `c78f753`, 344 lines, scaffold replaced)

Eleven testsets, **63 assertions**, all randomness on `p12_fix_rng(P12_FIXTURE_COUNTER)`.
`p12_rng` — the reported stream — is not touched anywhere in the file (`grep` count 0).
`P12_PENDING_SCAFFOLD` is gone (`grep -c` returns 0).

## Verification Results — what was actually run and observed

| Check | Command | Observed |
|---|---|---|
| Task 1 verify (verbatim from the plan) | the plan's `julia -e` one-liner | **PASS**, exit 0 — printed `0.116 prior-ok` |
| frozen/read-only paths unchanged | `git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml src corpus spike/simulator/{prior,ghat,forward}.jl` | **PASS** — in particular `forward.jl` is byte-unchanged by this plan |
| no re-declared prior constants | `grep -vE '^\s*#' \| grep -E 'P12_R1_PRIOR=\|MU_PRIOR=\|SIM02_W1_TOL=\|_PRIOR='` | **PASS** — zero matches |
| `ghat` applied elementwise / lattice reached | `grep -c 'ghat\.'` / `grep -c 'lattice_sigma'` | **PASS** — 3 and 4 |
| no per-pixel interpolant | `grep -c 'Interpolations'` | **PASS** — 0 |
| generative order, read from the source | line indices inside `sample_p12_prior`, comments stripped | **PASS** — r₁ at 3, field at 5, first nuisance at 13 |
| offset arm is a real arm | `p12_bilinear_op(8,64;offset=true) != (...;offset=false)` | **PASS**; every offset row still sums to exactly 1.0 |
| Task 2 verify | `julia --project=spike -e 'using Test; include("spike/test/test_p12_prior.jl")'` | **PASS** — 11 testsets, 63 passed, 0 failed, 0 errored, exit 0 |
| file runtime budget | same run, `time` | **9.85 s** cold (incl. compilation) against the plan's 60 s ceiling |
| eleven testsets green under `runtests.jl` | `julia --project=spike spike/test/runtests.jl` | **PASS** — 5 / 13 / 4 / 9 / 4 / 4 / 5 / 3 / 5 / 10 / 1 = **63, zero failures, zero errors**, after the `P12-SUITE-RAN: test_p12_prior.jl` marker |
| every earlier Phase-12 testset still green | same run | **PASS** — `test_p12_consts.jl`, `test_p12_lattice.jl` (137), `test_p12_architecture.jl`, `test_p12_decoupling.jl` (24/24) all report zero failures |
| decoupling (12-05's live self-scanning test) | same run | **PASS 24/24** — `src/` byte-unchanged, spike env byte-frozen, dependency name set exactly the frozen sixteen, sealed holdout not consumed, Phase-11 pool read-only |
| `spike/data/cache/p11` intact | `du -sh` | **PASS** — 54 MB |
| `.planning/STATE.md`, `.planning/ROADMAP.md`, `p12_consts.jl` unchanged | `git diff --quiet HEAD --` | **PASS** (this plan opens no Tier-2 sentinel) |

### Measured values the plan asked to record

| Quantity | Measured (M = 2 000, `arm = :car`, `r1 = 0.50`, fixture stream) | Bar |
|---|---|---|
| per-region W1 vs `MU_PRIOR`, **max** over 64 regions | **0.016056** (region 42 = `p12_idx(2,6)`) | `P12_SIM02_W1_TOL_PERREGION` = 0.10 |
| per-region W1, mean | 0.008277 | — |
| corner regions `(1,1) / (1,8) / (8,1) / (8,8)` | 0.01217 / 0.00773 / 0.00599 / 0.01355 | 0.10 |
| per-region W1 max at the **bottom** rung `r1 = 0.05` | 0.020338 | 0.10 |
| per-region `ghat` atom mass (ρ at ±0.99) | **0.065563** | reported, not gated (research 0.0679; global 1-D reference 0.0678) |
| per-cell sd of the Gaussian field the copula consumes | within `[0.93, 1.07]` on all 64 regions | — |
| empirical lattice lag-1, pinned `r1 = 0.90` vs `0.05` (200 draws each) | **0.89386** vs **0.04932** | — |
| ablation `arm = :none`, 1 000 draws: lag-1 / W1 max | −0.0017 / 0.0281 | \|lag-1\| < 0.05 ; W1 ≤ 0.10 |
| separable upsample vs the inline per-pixel reference at (129, 97) | **1.11e-16** | < 1e-12 |
| constant field 0.42 rendered to 512² | max deviation **0.0** (exact, this size) | < 1e-12 |

### The mutation proofs — every invariant was broken and observed to fail

Both mutations were applied **in memory** (redefining the function in `Main` after the include), so
no file on disk was edited and nothing was reverted. Verified `git diff --quiet` clean on
`p12_prior.jl` afterwards.

| Mutation | Testset 1 | Testset 2 | Testset 6 |
|---|---|---|---|
| `./ sqrt.(diag(Σ))` deleted from `field_sampler` | **PASSES** — max W1 0.0789 ≤ 0.10 | **FAILS 2 of 13** — the per-cell sd band (44 of 64 regions out of band, observed range [0.644, 1.002]) and the bottom-rung bar (**0.12639 > 0.10**) | not reached (the thrown testset aborts the file) |
| `A * F * B'` → `A * F * B` | — | — | **ERRORS** — `DimensionMismatch`, at (129, 97) *and* at a square (256, 256); the missing transpose can never silently return a wrong image |

The unrescaled **corner** W1s measure 0.0112 / 0.0080 / 0.0069 / 0.0115 — all comfortably inside
the bar. That is the executable proof that the plan's corner-only falsifier could not have failed.

**Not run:** `julia --project=. -e 'using Pkg; Pkg.test()'`. No pre-plan baseline exists (12-01 and
12-03 recorded the same gap), so "unchanged" is not something this run can assert. `src/` is
byte-unchanged by `git diff --quiet` and by 12-05's live decoupling testset, and this plan adds
nothing to the main package and loads nothing from it, so the main-package suite cannot have been
affected. Recorded rather than implied.

**The suite exit code is 1 and it is not ours** — see Deferred Issues below. Every Phase-12
testset runs and reports green before the abort.

## Deviations from Plan

### 1. [Rule 1 — Bug] Testset 2's premise is inverted, and as written the testset could not fail

- **Found during:** Task 2, probing the mutation before writing the assertions.
- **Issue:** The plan specifies testset 2 as *"the falsifier for a missing `sqrt.(diag(Σ))`
  rescale: without it these four **corner** regions are the ones that break"*, citing *"the
  measured edge/interior sd pair (0.658 / 0.992 at α = 0.95)"* — which is 12-RESEARCH.md Pattern 1
  and `p12_lattice.jl:41-48`, both of which read *"0.658 at an edge cell against 0.992 in the
  interior"*. **Measured this session on the committed `car_sigma`, the pairing is the other way
  round:**

  | Cell | α = 0.95 | at `r1 = 0.50` (α = 0.947559) |
  |---|---|---|
  | corner (1,1) | **0.9922** | **0.9863** |
  | edge-mid (1,4) | 0.7961 | 0.7910 |
  | interior (4,4) | **0.6575** | **0.6534** |

  `argmax(sqrt.(diag(Σ)))` is always a corner and `argmin` always an interior cell — as it must be,
  since `Σ = (D − αW)^{-1}` gives LOW-degree cells the LARGEST variance and corners have degree 2.
  The consequence is not cosmetic: the corners are the cells whose *unrescaled* marginal is already
  nearest 1, so they are exactly the regions that do **not** break when the rescale is deleted
  (measured corner W1 without the rescale: 0.0112 against 0.0122 with it). A testset asserting the
  falsifier only on the four corners is a test that cannot fail.
- **Compounding:** the pre-registered bar does not falsify the mutation at `r1 = 0.50` either
  (0.0789 ≤ 0.10), so testset 1 could not have covered for it.
- **Fix:** testset 2 keeps the four corner assertions (they hold, and the boundary claim is worth
  stating) and adds the three things that actually bite — (i) the structural sd assertions
  `argmax ∈ corners`, `argmin ∉ boundary`, `sd(1,1) > sd(1,4) > sd(4,4)`, which make the correction
  executable rather than a footnote; (ii) the per-cell sd band on the field the copula consumes,
  paired with the falsifier that `Σ` itself is *not* in that band; (iii) the **pre-registered bar at
  the bottom `P12_R1_LADDER` rung**, where the same deletion measures 0.1264 and breaches it. **No
  threshold was invented, re-derived or moved** — only the rung the pre-registered bar is also
  screened at. Both the correction and its measurement are written into the test file's header.
- **Files modified:** `spike/test/test_p12_prior.jl`.
- **Commit:** `c78f753`.

### 2. [Rule 3 — Blocking] `using Statistics` added to `p12_prior.jl`

- **Found during:** Task 1. `theta_scalar_view` needs `mean(draw.z_field)`. The plan lists only
  `using Distributions` and `using Random`. `Distributions` does re-export `Statistics.mean`, so
  the call would have resolved — but silently, through a re-export, which is exactly the kind of
  dependency nobody can see.
- **Fix:** `using Statistics`, declared explicitly with a header note that it is a stdlib reached
  through the default `@stdlib` LOAD_PATH exactly as `run_p11_recovery.jl:50-53` reaches it.
  **Nothing was added to `spike/Project.toml`** — verified byte-unchanged, and the frozen
  sixteen-name dependency assertion still passes.

### 3. Testsets 1-3 hoist the sampler; the ablation testset draws 1 000 instead of 200

- **Found during:** Task 2.
- **Hoisting:** `sample_p12_prior` re-solves the r₁ bisection on *every* call, which is right for a
  datagen loop of independent r₁ draws and pointless for a fixed-r₁ Monte Carlo. Testsets 1-3 build
  `field_sampler(lattice_sigma(:car, r1))` once and call `rho_field` — the object under test is
  identical. Testsets 4 and 5 go through the full `sample_p12_prior` entry point, so that path is
  still exercised.
- **1 000 vs 200 (testset 5):** the `:none` arm is the identity covariance, so it costs no bisection
  and 1 000 draws take 0.05 s. At M = 200 the screening W1 measures **0.0685** against a 0.10 bar —
  a two-thirds-of-the-bar reading that is Monte-Carlo noise, not signal, and a latent flake. At
  M = 1 000 it is **0.0281**. This is an *increase*, taken for the same reason
  `test_simulator.jl:213-216` raised N rather than loosening `SIM02_W1_TOL`: keep the honest
  pre-declared bar intact.

### 4. Assertions beyond the plan's enumerated list, in six testsets

- **What:** testset 1 adds the M-is-a-screen assertions (`P12_SIM02_M == 20_000`, `_P12P_M <
  P12_SIM02_M`); testset 3 adds a tighter informative band and asserts the atoms really are the
  `ghat` knot endpoints; testset 4 adds `d.arm === :car`, the `:gp` arm round-trip, the pinned-r₁
  round-trip and the determinism check `a == c`; testset 5 adds `P12_ABLATION_R1 < P12_R1_MIN` and
  the `:car`-at-ablation-r₁ identity; testset 6 adds the offset arm against the same reference and
  the hoisted-operator equality; testset 9 adds a positive case so the guard is not always-throw.
- **Why:** each covers a surface the plan's item list creates but does not assert. Same shape as
  12-03's deviation 2. No testset was added or removed — the plan's eleven are exactly what shipped.

### 5. Three descriptive inaccuracies in the plan, reported rather than silently corrected

- **Task 1 acceptance criterion** states the verify command *"prints `(8, 8)`, an `r1` inside
  `[P12_R1_MIN, P12_R1_MAX]`, `true`, and `(512, 512)`"*. The command as written prints only
  `round(d.r1;digits=3), " prior-ok"`; the other three facts are `@assert`s, not prints. The
  criterion's *printing* clause is therefore unsatisfiable as written, while the substance it
  describes is fully covered by the assertions. Observed output: `0.116 prior-ok`, exit 0.
- **`read_first` calls `prior.jl:37` a "guarded `include(joinpath(@__DIR__, "ghat.jl"))`."** It is
  at line 37 and it is **unguarded**. This matters slightly and in `p12_prior.jl`'s favour: the
  guard on *this* file's include of `prior.jl` is what keeps `ghat`'s knot vectors from being
  re-bound on a second load, and that is now stated in the include comment.
- **`read_first` cites `sample_prior` at `prior.jl:87-105`.** It is at `:93-105`; `:75-92` is its
  docstring. No consequence.

## Assumption Drift (advisory)

**Testset 1's pre-registered bar was assumed to be the thing that catches a broken copula in
seconds. It is not, at the rung the plan chose.**

- **Planned:** *"this testset exists so a broken copula fails in 120 s rather than in the reported
  run"* — i.e. the screening assertion `maximum(w1) ≤ P12_SIM02_W1_TOL_PERREGION` at
  `arm = :car, r1 = 0.5` is the protective assertion.
- **Actual, measured by mutation:** with the `./ sd` rescale deleted the max per-region W1 at
  `r1 = 0.50` is **0.0789** — badly wrong (the mean rises from 0.0083 to 0.0527, and the worst
  region is interior cell 27) yet still **inside** the 0.10 bar. Testset 1 passes on a broken
  copula. What catches it is testset 2, and only through the additions described in Deviation 1.
- **Why:** the bar is `SIM02_W1_TOL`'s value, fixed for a *scalar* induced-μ claim with a 13×
  headroom design margin, and it was (rightly) not reverse-engineered from any measured value. It
  was never sized as a mutation detector. The W1 scoping to `[GHAT_MU_MIN, GHAT_MU_MAX]` also damps
  the signal, because a compressed marginal moves the most mass in exactly the tails that scoping
  removes.
- **Why this is advisory:** nothing gated moved. The pre-registered bar is unchanged, is still
  asserted exactly as specified, and the reported claim (`run_p12_sim02.jl` at M = 20 000) is
  untouched. Recorded so no later document reads "testset 1 is green" as "the rescale is present".

## Deferred Issues

**`runtests.jl` exits 1 at `test_p13_consts.jl` with 114 `UndefVarError`s — pre-existing, not
this plan's, logged as DEF-12-03.**

`spike/p13/consts.jl:100` guards its entire Tier-1 block on `:P13_DEV_SEED`, and
`spike/validation/p12_consts.jl:109` legitimately binds that same name (R-5 obliges Phase 12 to
forbid the Phase-13 stream *by name*; the value is identical, so no seed changed and no
pre-registration was breached). The Phase-12 aggregator loads first, so
`test_p13_consts.jl:42`'s guarded include silently no-ops and every Phase-13 constant is missing.

Reproduced with **neither 12-07 file loaded** — including `test_p12_consts.jl` alone already leaves
`P13_DEV_SEED` defined and `P13_SALT` undefined. It surfaced only now because commit `a494247`
("PAUSE the NPE-03 wall-clock gate") converted the Phase-4 `SPEEDUP_GATE` to `@test_broken`, so the
suite reached the Phase-13 block for the first time (Phase 4 now reports 149 passed / 1 broken
instead of throwing).

Not fixed here: `p12_consts.jl` is frozen Tier-1 append-only, and the two Phase-13 files are owned
by a **live concurrent Phase-13 executor on this branch**. This is the fourth instance of one
structural class (13-09, DEF-12-01, DEF-12-02, DEF-12-03); the full mechanism, the reproduction and
the preferred one-line fix are written up in
`.planning/phases/12-spatial-colocalization-map/deferred-items.md`.

## Known Stubs

None. Both files are complete implementations; nothing in this plan is deferred to a later one.

## Threat Flags

None. No new network, auth or file-access surface: this plan adds in-memory sampling and linear
algebra plus a test file, writes no artifact, and reads only the frozen pre-registration and the
frozen `ghat`. The threat register's seven Phase-12 dispositions for this plan are all mitigated
and asserted — T-12-15 by testset 2 (in its corrected, falsifiable form), T-12-24 by the guarded
includes plus the no-re-declaration grep, T-12-25 by `arm` on the returned draw (testset 4),
T-12-17 by the separable operators and the absent per-pixel comprehension, T-12-26 by testset 4's
two behavioural order checks, T-12-27 by `validate_rho_field` and testset 9, and T-12-07 by
pathspec-scoped commits throughout.

## For the Next Plan

- **12-08** edits `forward.jl` stage 1. `forward.jl` is **byte-unchanged** by this plan, so the
  golden capture that must precede that edit is still capturable against pristine code.
  `p12_upsample` and `p12_upsample_ops` are the two matmuls stage 1 needs; hoist the operators.
  The exact-equality golden criterion belongs on the **legacy scalar-θ path**, which does not go
  through the upsample at all — a constant *field* agrees only to a few ULP, which is what
  `P12_STAGE1_CONSTFIELD_TOL` is for.
- **12-09** owns the 72-row θ assembly and is the only file that includes both this and
  `p12_architecture.jl`. Build the rows from `draw.z_field` (atom-free), not from `rho_field`, and
  carry `draw.r1` and `draw.arm` into the pool.
- **12-10 and 12-15 must not quote the "0.658 edge / 0.992 interior" pairing.** It is the wrong way
  round (Deviation 1). The corrected, measured numbers are in this summary and are asserted
  executably in `test_p12_prior.jl` testset 2. `12-RESEARCH.md` Pattern 1 and `p12_lattice.jl:41-48`
  still carry the inverted wording; `p12_lattice.jl` is not frozen, but editing it is out of this
  plan's scope and its *conclusion* (the rescale is mandatory) is unaffected and independently
  proved here.
- **12-10's reported SIM-02 runner** should use `P12_SIM02_M = 20 000` and the same scoping this
  file uses. Note that the screening bar alone does not detect a deleted rescale at `r1 = 0.50`
  (Assumption Drift): the reported runner should sweep the ladder, not a single rung.
- **A note for anyone adding a testset here:** 12-08 and 12-10 also append to
  `test_p12_prior.jl`. The shared 2 000-draw fixture is `const _P12P_F05` / `_P12P_W05`, computed
  once at file scope; reuse it rather than drawing again.
- The Phase-4 `SPEEDUP_GATE` no longer aborts `runtests.jl` (it is `@test_broken` since `a494247`),
  so the suite now runs to Phase 13 and dies there instead (DEF-12-03). Read the `P12-SUITE-RAN`
  markers and the Phase-12 testset summaries; do not read a red suite as a Phase-12 failure.

## Self-Check: PASSED

Both files verified present on disk (`spike/simulator/p12_prior.jl`, 305 lines;
`spike/test/test_p12_prior.jl`, 344 lines). Both task commits verified present in `git log`:
`6c77dfd` (the copula, the r₁-first draw, the separable upsample) and `c78f753` (the eleven
testsets). `spike/Project.toml`, `spike/Manifest.toml`, `src/`, `corpus/`, `spike/simulator/prior.jl`,
`spike/simulator/ghat.jl` and `spike/simulator/forward.jl` byte-unchanged;
`spike/validation/p12_consts.jl` byte-unchanged (this plan opens no Tier-2 sentinel);
`.planning/STATE.md` and `.planning/ROADMAP.md` untouched; `spike/data/cache/p11` intact at 54 MB;
`artifacts/amended_v2/grid_8` untouched. Both commits used explicit pathspecs; no `git stash`,
`rebase`, `amend`, `reset` or `clean` was run, and no index.lock race occurred despite the
concurrent Phase-13 executor on the same branch.
