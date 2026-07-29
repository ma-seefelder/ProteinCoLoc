---
phase: 12-spatial-colocalization-map
plan: 08
status: complete
subsystem: forward simulator (D-06 stage-1 spatial rho)
tags: [golden-fixture, capture-before-edit, byte-preservation, stage-1, bilinear-upsample, mutation-tested]
requires: ["12-01", "12-03", "12-07"]
provides:
  - "spike/test/capture_p12_golden.jl — one-shot PRE-EDIT capture: dirty-tree refusal on spike/simulator, capture-time sha embedded as p12_base_sha, explicit 8-field theta literal with NO rho_field, atomic .tmp -> integrity reopen -> mv"
  - "spike/test/fixtures/p12_stage1_golden.jld2 — TRACKED pre-edit stage-1 golden, six FIXTURE-stream keys at 128^2"
  - "spike/simulator/forward.jl — stage 1 accepts an OPTIONAL theta.rho_field through one defensive hasproperty binding; _p12_check_field entry guard; _p12_upsample_local (A*F*B'); scalar path byte-unchanged"
  - "spike/test/test_p12_prior.jl — testsets 12-16 appended (61 new assertions; 16 testsets / 124 assertions total)"
affects: ["12-09", "12-11", "12-20"]
tech-stack:
  added: []
  patterns:
    - "capture-BEFORE-you-edit: the fixture is committed in its own commit, at a sha whose tree carries the pre-edit file, and the sha is asserted to be an ancestor of HEAD"
    - "a criterion that measures 0.0 is not automatically a falsifier: prove the new branch ran by making it produce a DIFFERENT answer, not a slightly different one"
    - "where two plans' validators are asserted to agree, measure first — the one divergence is a finding to record, not a test to force"
key-files:
  created:
    - spike/test/capture_p12_golden.jl
    - spike/test/fixtures/p12_stage1_golden.jld2
  modified:
    - spike/simulator/forward.jl
    - spike/test/test_p12_prior.jl
decisions:
  - "MEASURED CORRECTION: a constant-VALUED rho field does NOT differ from the scalar path 'in the last few ULP'. It is EXACTLY equal — max|delta| = 0.0 across all six keys, because _p12_upsample_local(fill(0.7,8,8), sz) is exactly 0.7 in every pixel at 128^2/256^2/512^2. The pre-registered criterion (<= P12_STAGE1_CONSTFIELD_TOL) therefore holds; the plan's paired falsifier `0 < max|delta|` is measurably FALSE and its stated reasoning is inverted. Mutation-proved: with the field branch DISABLED the difference is also 0.0, so zero difference proves nothing either way."
  - "The replacement falsifier is an OFF-VALUE constant field: rho_field = fill(0.1) against theta.rho_true = 0.7 must override the scalar. Measured 3.04 against the golden, and 0.0 against a scalar theta at 0.1. With the field branch disabled it collapses to 0.0 and the assertion goes red."
  - "MEASURED DIVERGENCE: validate_rho_field (12-07) ACCEPTS a 1x1 field; simulate_pair's inlined _p12_check_field REJECTS it. Not a bug in either — the >=2-cells-per-axis rule lives in p12_bilinear_op's G>=2 guard on the p12_prior.jl side, so p12_upsample refuses it too. Testset 16 asserts the four-way agreement AND this one divergence, rather than asserting a falsehood."
  - "P12_GOLDEN_IMSIZE and P12_GOLDEN_KEYS are bound in the capture script, NOT in p12_consts.jl: they are fixture geometry, nothing gates on them, and Tier 1 is closed."
metrics:
  duration: ~95 min
  completed: 2026-07-29
  tasks: 3
  commits: 4
  files: 4
---

# Phase 12 Plan 08: The Pre-Edit Golden and the Optional ρ Field — Summary

The forward simulator now mixes ρ **per pixel** when a θ hands it a G×G field, and a θ that does
not is **byte-identical** to the simulator as it stood before the edit — proved against a fixture
captured and committed **one commit earlier**, at a sha whose tree still carries the pre-edit file.

## The three numbers the plan asked for, copied from command output

| Quantity | Value |
|---|---|
| `p12_base_sha` (embedded in the fixture, printed by the capture script) | **`a889dc4f715b3046399c4b69d7ff00649678a709`** |
| Exact arm — a θ with no `rho_field` vs the golden bytes, six keys, both channels | **12 / 12 `==`**, zero failures |
| Constant-field arm — `maximum(abs, Δ)` vs the golden bytes, six keys | **0.0** (bar `P12_STAGE1_CONSTFIELD_TOL` = 1e-12) |

## What Was Built

**Task 1 — `capture_p12_golden.jl` + its fixture** (commit `c8d128e`, 154 lines + 1.51 MB fixture)

Structural mirror of `capture_p11_golden.jl`: the timing-discipline header in the analog's voice,
`git status --porcelain -- spike/simulator` refusal (it **passed**, the tree was clean at capture
time), the 40-char sha assertion, the explicit `THETA_GOLDEN_P12` literal with in-file
`!hasproperty(..., :rho_field)` and `length(...) == 8` assertions, and the atomic
`.tmp → jldsave → reopen-and-`@assert haskey`` → `mv(...; force = true)` write.

The θ is the Phase-11 literal plus an **explicit** `chromatic_eps = 0.0`, so the fixture is about
stage 1 alone and the stage-6 term is pinned inert. Six keys at 128², riding
`p12_fix_rng(P12_FIXTURE_COUNTER + k)` — the FIXTURE stream on the FIXTURE counter.

**Task 2 — `forward.jl` stage 1** (commit `1ccf5cf`, +95 / −4)

| Piece | What it is |
|---|---|
| `rho_field_val = hasproperty(θ, :rho_field) ? θ.rho_field : nothing` | ONE defensive binding, immediately after `chromatic_eps_val`, copying that line's discipline exactly |
| `_p12_check_field(F)` | entry guard: matrix, square, both dims ≥ 2, finite, elementwise in `[-1,1]`; throws `ArgumentError` naming the offending value, **before any pixel work** |
| `_p12_upsample_local(F, imsize)` | `A * F * B'`, `offset = false` only; the half-cell guard arm stays in `p12_prior.jl` alone |
| the stage-1 branch | the plan's code verbatim; the three `_smooth_field` calls stay **before** it, so RNG consumption is identical on both paths |
| docstring | one added sentence naming the optional `rho_field` |

`forward.jl` gained **zero** `include`s and **zero** `using`s
(`git diff HEAD -- spike/simulator/forward.jl | grep -c '^+using\|^+include'` → `0`).

**Task 3 — testsets 12-16 appended to `test_p12_prior.jl`** (commits `50e5d0a`, `d200391`, +274)

12-07's eleven testsets and 63 assertions are untouched. The five new ones add **61** assertions:
fixture-precedes-edit (11), the exact arm plus its source-level anti-softening tripwire (20), the
constant-field criterion plus the override falsifier (5), the spatial check (10), the two-validator
cross-check (15).

## Verification Results — what was actually run and observed

| Check | Command | Observed |
|---|---|---|
| Task 1 verify (verbatim from the plan) | the plan's `julia -e` one-liner | **PASS**, exit 0 — printed `a889dc4f715b3046399c4b69d7ff00649678a709` |
| capture-time dirty-tree refusal | inside the script | **PASS** — `git status --porcelain -- spike/simulator` was empty |
| fixture is TRACKED, not gitignored | `git check-ignore -q …` | **PASS** — exit 1 (not ignored); committed in `c8d128e` |
| re-capture on an unchanged tree is byte-identical | captured twice, `a["channels"] == b["channels"]` | **PASS** — `true` (Philox is counter-based); only `generated` differs; **the first capture was restored and is what is committed** |
| Task 2 verify | `julia --project=spike -e 'using Test; include("spike/test/test_stage6_regression.jl")'` | **PASS** — **51 / 51**, exit 0, incl. the source-level `count("warp(", code) == 1` tripwire |
| stage 6 + the constants block byte-identical to the pre-edit blob | line-block compare against `git show a889dc4:…` | **PASS** — both `true`; every diff hunk sits at old line ≤ 171 (stage 1 and above) |
| zero new `using` / `include` | `git diff … \| grep -c` | **PASS** — `0` |
| Task 3 verify | `julia --project=spike -e 'using Test; include("spike/test/test_p12_prior.jl")'` | **PASS** — **16 testsets, 124 assertions, 0 failed, 0 errored**, exit 0, 21.3 s cold |
| all sixteen green under `runtests.jl` | `julia --project=spike spike/test/runtests.jl` | **PASS** — 5/13/4/9/4/4/5/3/5/10/1 + **11/20/5/10/15**, zero failures, after the `P12-SUITE-RAN: test_p12_prior.jl` marker |
| every earlier Phase-12 testset still green | same run | **PASS** — `test_p12_consts.jl`, `test_p12_lattice.jl`, `test_p12_architecture.jl` all zero failures |
| decoupling (12-05's live self-scanning test) | same run | **PASS 24 / 24** — `src/` byte-unchanged, spike env byte-frozen, the frozen sixteen deps, sealed holdout unconsumed, Phase-11 pool read-only |
| `P12_PENDING_SCAFFOLD` gone | `grep -c` | **PASS** — `0` |
| the reported stream is not pre-observed | `grep -c 'p12_rng('` in the test file | **PASS** — `0`; all randomness rides `p12_fix_rng` |
| frozen / read-only paths | `git diff --quiet HEAD -- src corpus spike/Project.toml spike/Manifest.toml spike/validation/p12_consts.jl` | **PASS** |
| `.planning/STATE.md`, `.planning/ROADMAP.md` | `git diff --quiet HEAD --` | **PASS** — untouched (orchestrator-owned) |
| `spike/data/cache/p11` intact | `du -sh` | **PASS** — 54 MB |

### Fixture provenance, both ids recorded as the plan's verification asks

```
git show a889dc4…:spike/simulator/forward.jl | git hash-object --stdin
    -> 28c3ce6f9d57b6e7efa9e8dce8d11844aa5756b0
git rev-parse a889dc4…:spike/simulator/forward.jl
    -> 28c3ce6f9d57b6e7efa9e8dce8d11844aa5756b0     (identical — the fixture's tree is pre-edit)
git rev-parse HEAD:spike/simulator/forward.jl
    -> baa797cde18920a3c9b22fcec765afb6638f567a     (post-edit)
git merge-base --is-ancestor a889dc4… HEAD          -> exit 0
```

The pre-edit blob was also identical in the **working tree** at capture time, so the fixture
records the committed bytes and not a local variant.

### Measured values

| Quantity | Measured | Bar / reference |
|---|---|---|
| exact arm, six keys × two channels | **12/12 `==`** | `P12_STAGE1_EXACT = true` |
| constant field at `ρ_true = 0.7`, `max|Δ|` vs golden, six keys | **0.0** | ≤ 1e-12 |
| `_p12_upsample_local(fill(0.7,8,8), sz) .- 0.7`, `sz ∈ {128²,256²,512²}` | **0.0**, every entry `==` | — |
| off-value field `fill(0.1)` vs golden (the branch-ran falsifier) | **3.0399999999999996** | > 0.5 |
| the same off-value field vs a **scalar** θ at 0.1 | **0.0** | ≤ 1e-12 |
| across-region patch-correlation spread, key 1: control → field | **0.15500 → 0.34057** (2.197×) | directional, not a gate |
| across-region spread, key 2: control → field | **0.21973 → 0.33278** (1.514×) | directional |
| the same on all six keys (probe, not asserted) | field larger on **6 / 6** | — |
| field-arm weakest region vs control's weakest, key 1 | **−0.597 vs +0.212** | the control never goes negative |
| worst `max|ρ_px|` over adversarial `|F| = 1` fields at six sizes incl. 1376×1028 | **1.0 exactly**, never > 1 | `sqrt(1 − |ρ_px|)` cannot `DomainError` |

The last row is a safety probe the plan did not ask for: `_p12_check_field` admits an entry of
exactly `±1`, so if the upsample could overshoot by one ULP, `sqrt(1.0 - abs(ρpx))` would throw
`DomainError` inside a 50 000-pair datagen loop. It cannot — the weights are a non-negative convex
combination and the measured worst case is exactly 1.0. **No clamp was added**, so the plan's stage-1
code block stands verbatim.

### The mutation proofs — every claim broken and observed to fail

Applied to **scratch copies** in the session temp directory. The tracked `forward.jl` was read only;
`git hash-object` afterwards equals `git rev-parse HEAD:spike/simulator/forward.jl` (`baa797c…`).

| Mutant | Testset 13 (exact arm) | Testset 14 criterion | Testset 14 falsifier |
|---|---|---|---|
| *(baseline, committed code)* | GREEN, 0/12 failing | GREEN, `Δ = 0.0` | GREEN, `3.04 > 0.5` |
| **A** — the three `_smooth_field` calls reordered | **RED, 12/12 failing** | RED, `Δ = 5.541` | GREEN |
| **B** — the field branch never runs (`rho_field_val = nothing`) | GREEN *(correct: the no-field path is unaffected)* | **GREEN, `Δ = 0.0`** | **RED, `0.0`** |

Mutant A is the plan's requested proof that the exact arm is sensitive to what it guards.
**Mutant B is the decisive one**: with the field branch disabled the constant-field difference is
*also* exactly 0.0. So a zero difference is consistent with the branch running **and** with it not
running, which is why the plan's `0 < max|Δ|` reasoning fails and why the falsifier had to move to
the off-value form — the only assertion in the pair that separates the two worlds.

## Deviations from Plan

### 1. [Rule 1 — Bug] Testset 14's falsifier `0 < maximum(abs, Δ)` is measurably false, and its stated reasoning is inverted

- **Found during:** Task 3, probing before writing the assertions.
- **Plan text:** *"The falsifiable form is that the constant-field path is within tolerance but NOT
  byte-identical to the scalar path: `@test 0 < maximum(abs, Δ) ≤ P12_STAGE1_CONSTFIELD_TOL`. The
  strict lower bound is the load-bearing half — if the field branch were skipped the result would be
  byte-identical … so a difference of precisely zero here proves the branch did NOT run."* The same
  premise appears in the plan's `<context>` and, in gentler form, in `p12_consts.jl:710-714` and
  `12-VALIDATION.md`.
- **Measured:** the constant-field difference is **exactly 0.0**, on all six keys, in both channels —
  and the branch **did** run. `_p12_upsample_local(fill(0.7, 8, 8), sz)` returns exactly `0.7` in
  every pixel at 128², 256² and 512² (`all(==(0.7), P)` is `true`), so `a`, `b` and `sgn` are
  bit-identical to the scalar branch's and all seven stages reproduce the golden bytes. 12-07's
  testset 8 measured the same thing independently at 0.42 / 512². As written the assertion would
  have **failed on correct code**.
- **And the reasoning is inverted, not merely the constant:** mutant B (field branch disabled) also
  measures 0.0. Zero difference therefore proves *nothing*; the plan's load-bearing half was
  load-bearing in the wrong direction.
- **Fix:** the pre-registered criterion `max|Δ| ≤ P12_STAGE1_CONSTFIELD_TOL` is asserted **unchanged
  and unmoved**, the measured value is recorded with `@info`, and the branch-ran proof is done where
  it actually bites — an **off-value** constant field (`fill(0.1)` against `ρ_true = 0.7`) must
  override the scalar: measured 3.04 against the golden, 0.0 against a scalar θ at 0.1, and 0.0
  under mutant B, where the assertion goes red. **No threshold was invented, re-derived or moved.**
  The whole correction is written into the test file's header.
- **Files modified:** `spike/test/test_p12_prior.jl`. **Commit:** `50e5d0a`.

### 2. [Rule 1 — Bug] Testset 16's shared bad-input list contains one input the two validators do not agree on

- **Found during:** Task 3, before writing the testset.
- **Plan text:** *"for a shared list of bad inputs (an entry > 1, an entry < −1, a NaN, a non-square
  matrix, **a 1×1 matrix**), assert **both** `validate_rho_field` and a `simulate_pair` call carrying
  that field throw."*
- **Measured:**

  | input | `validate_rho_field` | `simulate_pair(rho_field=…)` | `p12_upsample` |
  |---|---|---|---|
  | entry > 1 | `ArgumentError` | `ArgumentError` | no throw |
  | entry < −1 | `ArgumentError` | `ArgumentError` | no throw |
  | `NaN` | `ArgumentError` | `ArgumentError` | no throw |
  | non-square | `ArgumentError` | `ArgumentError` | no throw |
  | **1×1** | **no throw** | `ArgumentError` | `ArgumentError` |

  A 1×1 matrix is square, finite and in range, so 12-07's validator accepts it by design. As
  specified the assertion would have **failed**.
- **This is not a defect in either file.** The "≥ 2 cells per axis" rule is a constraint of the
  **renderer**, and on the `p12_prior.jl` side it lives in `p12_bilinear_op`'s `G ≥ 2` guard rather
  than in the validator; `forward.jl`'s inlined `_p12_check_field` merges validation and rendering
  into one entry guard because it has no separate upsample to fail in. The two sides together reject
  exactly what the merged guard rejects.
- **Fix:** testset 16 asserts the four-way agreement on the four inputs where the validators agree,
  and asserts the 1×1 divergence **explicitly and in all three directions** (validator accepts,
  `p12_upsample` throws, `simulate_pair` throws), so a future edit to either side is caught rather
  than papered over. It also adds `_p12_upsample_local == p12_upsample` at two sizes, which is the
  upsample half of the same bounded-duplication claim.
- **Files modified:** `spike/test/test_p12_prior.jl`. **Commit:** `50e5d0a`.

### 3. [Rule 3 — Blocking] `using JLD2` and a guarded `include` of `contract.jl` added to the test file

- Testsets 12-14 read the JLD2 fixture; testset 15 needs the real SIM-03 patch summary
  (`build_mci` / `patch_summary`), which lives in `spike/contract.jl`. Both are additions the plan
  implies but does not list. The include is **guarded** on `:build_mci` (already loaded inside the
  aggregated suite) and `contract.jl`'s two `src/` includes are read-only. Measured cost when not
  already loaded: **0.84 s**. **Nothing was added to `spike/Project.toml`** — verified byte-unchanged,
  and 12-05's frozen sixteen-name dependency assertion still passes 5/5.

### 4. Assertions beyond the plan's enumerated list

Testset 12 adds the fixture's stream/counter provenance (`fixture_counter`, `seed`, `salt`) and
`chromatic_eps == 0.0`; testset 13 adds `@test P12_STAGE1_EXACT == true` and the four
non-vacuity assertions on the source tripwire (`!isempty`, `occursin("_P12G_CHANS")`); testset 14
adds the off-value-vs-scalar equality and `_P12G_OFFVAL != ρ` so the falsifier cannot pass by
comparing `ρ_true` with itself; testset 15 adds the weakest-region assertion that separates
"spatially varying" from "merely noisier"; testset 16 adds the upsample equality and a CPU-only
check for the half of the file that loads `contract.jl`. Same shape as 12-07's deviation 4.

### 5. Two descriptive inaccuracies in the plan, reported rather than silently corrected

- **The `<verification>` block asks that `runtests.jl` show "the pre-existing
  `test_stage6_regression` testsets green".** `test_stage6_regression.jl` is **not included in
  `runtests.jl` at all** (the aggregated suite's include list runs `test_p12_suite.jl`,
  `test_simulator.jl`, `test_data_pipeline.jl`, `test_npe.jl`, `test_sbc.jl`, `test_bf.jl`,
  `test_ood.jl`, `test_comparator.jl`, then the eleven Phase-13 files). It was run **per file**
  instead — 51/51 green — which is the briefing's sanctioned gate signal. Nothing is wrong with the
  test; the verification line is not satisfiable as written.
- **`read_first` cites `test_stage6_regression.jl:88-100` for `_strip_comment_lines`.** The helper is
  at `:79-80`; `:88-100` is inside the `SC1c` testset, which is where the *tripwire that uses it*
  lives (`:95-96`). No consequence — both were read.

## Assumption Drift (advisory)

**The plan assumed the constant-field arm and the exact arm are two different strengths of the same
claim. Measured, at these operating points they are the same claim.**

- **Planned** (`<context>`, and `p12_consts.jl:710-714` behind it): the exact `==` arm holds only on
  the no-field path, while a constant *field* "cannot be bit-identical" and agrees only to a few ULP,
  which is why `P12_STAGE1_CONSTFIELD_TOL` exists as *a criterion and not an escape hatch*.
- **Actual:** a constant field is bit-identical too, at 128², 256² and 512², at ρ = 0.7 and at
  ρ = 0.1. The tolerance is never consumed.
- **Why:** the row weights `1-t` and `t` do round, but on a *constant* field the products
  `(1-t)·c` and `t·c` recombine to exactly `c` at every sampled position for these `(G, m)` pairs, so
  the rounding never materializes. It is a property of the geometry, not a guarantee — a different
  `(G, m)` or a different constant could consume some of the tolerance.
- **Why this is advisory:** nothing gated moved. `P12_STAGE1_CONSTFIELD_TOL` is Tier-1 frozen, is
  asserted exactly as pre-registered, and is the right bound to keep — it is precisely the headroom
  that makes the claim robust to the `(G, m)` pairs 12-09 and 12-11 will actually use. Recorded so no
  later document reads "the constant-field arm is a tolerance arm" as evidence that the upsample is
  lossy here.

## Deferred Issues

**`runtests.jl` still exits 1 at `test_p13_consts.jl` with 114 `UndefVarError`s — pre-existing, not
this plan's, already logged as DEF-12-03 by 12-07.** Reproduced unchanged this session (22 passed /
1 failed / 114 errored, aborting at `test_p13_consts.jl:52`). The mechanism is the include-guard
poisoning class: `spike/p13/consts.jl:100` guards its Tier-1 block on `:P13_DEV_SEED`, which
`p12_consts.jl:109` legitimately binds. Every Phase-12 testset runs and reports green **before** the
abort. Not fixed here: `p12_consts.jl` is frozen Tier-1 append-only and the Phase-13 files are owned
by a live concurrent executor on this branch.

## Known Stubs

None. All four files are complete; nothing in this plan is deferred to a later one.

## Threat Flags

None. No new network, auth or file-access surface. The plan's six dispositions are all mitigated and
asserted: **T-12-28** by `p12_base_sha` + the dirty-tree refusal + testset 12's ancestry and
no-`rho_field` assertions; **T-12-29** by the preserved scalar branch, the pre-branch
`_smooth_field` calls, testset 13's six-key `==` and its source-level anti-softening tripwire —
mutation-proved red under mutant A; **T-12-30** by `test_stage6_regression.jl`'s
`count("warp(", code) == 1` still green and no `warp` in the new code; **T-12-31** by zero new
`include`/`using` (grep-asserted) and testset 16's cross-check; **T-12-27** by `_p12_check_field`
throwing before any pixel work; **T-12-07** by four pathspec-scoped commits, no index.lock race
despite the concurrent Phase-13 executor.

## For the Next Plan

- **12-09 (datagen):** stage 1 now consumes `θ.rho_field` directly, so build the 72-row θ and hand
  the field through unchanged. **Hoist the operators** — `_p12_upsample_local` rebuilds `A` and `B`
  on every call by design (it is a private copy with no hoisting API); the pool loop should use
  `p12_upsample_ops` where it renders anything itself. Fields must be `G×G` with **both dims ≥ 2**;
  a 1×1 is accepted by `validate_rho_field` and rejected by the simulator (Deviation 2).
- **12-11 / 12-20:** the D-06 half-cell guard arm is reached by **pre-rendering an offset field into
  a θ** (`p12_upsample(F, imsize; offset = true)` → feed the result's lattice as `rho_field`), never
  by a second branch in `forward.jl`. There is exactly one code path in the simulator.
- **12-10, which also appends to `test_p12_prior.jl`:** the file is now 16 testsets / 124 assertions
  and ~21 s cold. Reuse `_P12G`, `_P12G_THETA`, `_P12G_CHANS`, `_P12G_KEYS` and `_p12g_spread` rather
  than reloading the fixture, and keep the `_p12g_*` / `_P12G_*` phase prefix.
- **Do not quote "a constant field differs in the last few ULP" as a measured fact.** It is 0.0 here
  (Assumption Drift). `P12_STAGE1_CONSTFIELD_TOL` remains the right bound and remains Tier-1 frozen.
- **`test_stage6_regression.jl` is not in `runtests.jl`.** Any plan that needs the Phase-11 stage-6
  guarantee must run that file explicitly; a green aggregated suite does not cover it.

## Self-Check: PASSED

All four files verified present on disk: `spike/test/capture_p12_golden.jl` (154 lines),
`spike/test/fixtures/p12_stage1_golden.jld2` (1.51 MB, tracked, `git check-ignore` exit 1),
`spike/simulator/forward.jl` (+95/−4, blob `baa797cde18920a3c9b22fcec765afb6638f567a`),
`spike/test/test_p12_prior.jl` (618 lines). All four commits verified present in `git log`:
`c8d128e` (capture + fixture), `1ccf5cf` (stage 1), `50e5d0a` (testsets 12-16), `d200391` (header
call-count correction). Every commit used explicit pathspecs; no `git stash`, `rebase`, `amend`,
`reset`, `clean` or `checkout --` was run, and the two mutations were applied to scratch copies
outside the repository. `src/`, `corpus/`, `spike/Project.toml`, `spike/Manifest.toml` and
`spike/validation/p12_consts.jl` byte-unchanged; `.planning/STATE.md` and `.planning/ROADMAP.md`
untouched; `spike/data/cache/p11` intact at 54 MB; `artifacts/amended_v2/grid_8` untouched.

**Not run:** `julia --project=. -e 'using Pkg; Pkg.test()'`. No pre-plan baseline exists (12-01,
12-03 and 12-07 recorded the same gap), so "unchanged" is not something this run can assert. The
substance is covered by stronger evidence: the root package has its own `src/amortized/simulator.jl`
and contains **no include of `spike/`** (grep: the only `spike/` mentions in `src/` are provenance
comments), `src/` is byte-unchanged by `git diff --quiet` and by 12-05's live 24/24 decoupling
testset, and this plan adds nothing to the main package and loads nothing from it. Recorded rather
than implied.
