---
phase: 11-registration-and-chromatic-uncertainty-as-latent
plan: 02
subsystem: simulator
tags: [julia, affine-warp, coordinate-transformations, theta-arity, provenance, content-hash, named-limits]

# Dependency graph
requires:
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 01
    provides: "spike/test/fixtures/p11_stage6_golden.jld2 (the PRE-EDIT stage-6 golden bytes) and PHASE11_BASE_SHA = 17ebd1e, without which the D-10 regression claim is unverifiable"
  - phase: 07-productionization-conditional-on-go
    provides: "the shipped amended_v2/grid_8 bundle, DATAGEN_HASH_SRC_FILES provenance machinery, and the docs/amortized.md named-limits list this extends"
provides:
  - "chromatic_eps as an 8th theta field (Uniform(-0.02, 0.02)) in BOTH mirrored simulators"
  - "stage 6 as ONE composed AffineMap (Translation o recenter(LinearMap(scale))) — a single interpolation pass"
  - "theta_prior_bounds() returning 8 derived rows; DATAGEN_THETA_DIM as the derived theta-arity single source of truth"
  - "_theta_tuple emitting chromatic_eps only when the vector is long enough (shipped 7-row compat hinge)"
  - "docs/amortized.md named limit #8 — provenance-by-reference, pinning 17ebd1e"
  - "EPS_COMMIT_SHA = ca02b0e8a696269a2b8710b11ca3073edc5ae6c7"
affects: [11-03-stage6-regression, 11-05-probe, 11-07-research-trainer, 11-ladder-and-breakdown, 11-11-interpretation-note]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Composed affine collapsed to a single AffineMap so `warp` performs exactly one interpolation pass"
    - "Defensive NamedTuple field read (`hasproperty`) as a structural backward-compatibility hinge, bound once and shared by guard and use site"
    - "Theta arity derived from `length(theta_prior_bounds())` at every buffer/assert site instead of a literal"

key-files:
  created:
    - .planning/phases/11-registration-and-chromatic-uncertainty-as-latent/deferred-items.md
  modified:
    - docs/amortized.md
    - spike/simulator/prior.jl
    - spike/simulator/forward.jl
    - spike/test/test_p11_forced_theta.jl
    - src/amortized/simulator.jl
    - src/amortized/datagen.jl
    - src/amortized/ood.jl
    - src/amortized/train_npe.jl
    - test/runtests.jl
    - test/gpu_smoke.jl
    - test/test_local_map.jl
    - test/test_integration.jl

key-decisions:
  - "Named limit #8 quotes BOTH the fixture-pinned PHASE11_BASE_SHA (17ebd1e) and the epsilon commit's immediate parent (a2ced26), with the byte-identity of simulator.jl across that span stated — satisfying the orchestrator's pin and the plan's first-parent criterion at once"
  - "train_npe now derives the flow's marginal count from `size(theta, 1)` instead of NPE_D; NPE_D and build_estimator's default are untouched, so the shipped bundle's frozen D is unaffected"
  - "test/test_integration.jl keeps its literal 7 — it pins the SHIPPED bundle's frozen marginal count, which the prior extension must NOT move"
  - "test/runtests.jl:495 asserts against ProteinCoLoc.NPE_D rather than the theta arity, because build_estimator called without D uses the shipped default"

requirements-completed: [D-02, D-09, D-10, D-11, D-12, D-15]

# Metrics
duration: ~85min
completed: 2026-07-25
---

# Phase 11 Plan 02: Chromatic ε as an 8th θ Column and a Single Composed Affine Stage 6 Summary

**`chromatic_eps` is now an 8th θ field in both mirrored simulators and stage 6 is ONE composed `AffineMap` — measured bit-identical to the pre-edit golden fixture at `chromatic_eps = 0` — landed in a single commit together with `docs/amortized.md` named limit #8, which pins the shipped `grid_8` bundle's training distribution to the pre-ε sha `17ebd1e`.**

## The two shas

```
PRE_EPS_SHA (pinned, fixture-embedded PHASE11_BASE_SHA)  17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348
PRE_EPS_SHA (immediate first parent of the ε commit)     a2ced26d93c7409a465aa50320e188868b3ab462
EPS_COMMIT_SHA                                            ca02b0e8a696269a2b8710b11ca3073edc5ae6c7
```

**Verified parent relationship (commands run, output observed):**

| Check | Command | Observed |
|---|---|---|
| First parent of the ε commit | `git rev-parse HEAD^` | `a2ced26d93c7409a465aa50320e188868b3ab462` |
| Parent is quoted in the limit | `grep -c a2ced26… docs/amortized.md` | `1` |
| Pinned sha is quoted in the limit | `grep -c 17ebd1e… docs/amortized.md` | `2` |
| Pinned sha predates the parameter | `git show 17ebd1e…:src/amortized/simulator.jl \| grep -c chromatic_eps` | `0` |
| Parent predates the parameter | `git show a2ced26…:src/amortized/simulator.jl \| grep -c chromatic_eps` | `0` |
| `simulator.jl` unchanged across the span | `git diff --stat 17ebd1e…a2ced26… -- src/amortized/simulator.jl` | empty |
| Fixture agrees with the pinned sha | `JLD2.load(...)["phase11_base_sha"]` | `17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348` |

Both quoted shas are genuinely pre-ε and `src/amortized/simulator.jl` is byte-identical at both, so the
`git show <SHA>:src/amortized/simulator.jl` reconstruction in limit #8 is correct against either.

## Accomplishments

- **One commit carries the provenance event and its repair.** `ca02b0e` contains the simulator
  edit, the θ-arity ripple and the `docs/amortized.md` limit together. No window exists in which
  `src/` claims a training distribution it no longer reproduces, and the cache was not re-digested
  twice.
- **Stage 6 is a single interpolation pass, proven bit-exact.** Replayed against the plan-11-01
  golden fixture at all four captured Philox keys: the composed affine at `chromatic_eps = 0`
  reproduces the pre-edit bytes with `==` (no tolerance), for BOTH an explicit 8-field θ and a
  legacy 7-field θ. The `WR-06` comment is re-derived for backward-mode `warp`, including why the
  scale is `1/(1 + chromatic_eps)` and why the composition collapses to one `AffineMap`.
- **The shipped-compat hinge holds on both sides.** `simulate_pair` reads `chromatic_eps` through
  one defensive `hasproperty` binding shared by the guard and stage 6, and `_theta_tuple` emits the
  8th field only at `length(v) >= 8`. A 7-row posterior mean — what the shipped bundle produces —
  therefore still round-trips into a valid θ.
- **θ arity is derived, not retyped.** `DATAGEN_THETA_DIM = length(theta_prior_bounds())` is the
  single source of truth for datagen's buffers and the hashed `theta_dim`; every test site uses
  `length(ProteinCoLoc.theta_prior_bounds())` rather than a fresh literal `8`.
- **The two deliberate non-mirrors are asserted, not just commented.** The src `SHIFT_PRIOR` stays
  `Uniform(-1.0, 1.0)` (ruling Q2) and `NPE_D` stays 7; both are pinned by assertions in the new
  testset and by the untouched shipped-bundle assertions in `test/test_integration.jl`.

## Verification — observed results

Every command below was run and its real output observed.

| Check | Command | Observed |
|---|---|---|
| Root suite (baseline, pre-edit) | `julia --project -t auto -e 'using Pkg; Pkg.test()'` | **passed**, exit 0 |
| Root suite (final) | `julia --project -t auto -e 'using Pkg; Pkg.test()'` | **passed**, exit 0 |
| New testset ran | `chromatic ε (D-09)` line in the suite output | **25 pass / 25 total**, 3.0s |
| D-10 golden regression | replay of the 4 fixture keys, `==` on both channels | legacy 7-field **true**, 8-field ε=0 **true** |
| Spike prior contract | the plan's `sample_prior`/`SHIFT_PRIOR`/`CHROMATIC_PRIOR` assert chain | prints `ok` |
| Spike legacy θ | `simulate_pair` on an explicit 7-field θ | prints `legacy 7-field ok` |
| src θ contract | the plan's `theta_prior_bounds`/`SHIFT_PRIOR`/`_theta_tuple` assert chain | prints `ok` |
| Spike forced-θ | `julia --project=spike spike/test/test_p11_forced_theta.jl` | **28 pass / 0 fail**, exit 0 |
| Spike pre-registration | `julia --project=spike spike/test/test_p11_consts.jl` | exit 0 |
| Spike equivalence stats | `julia --project=spike spike/test/test_p11_tost.jl` | exit 0 |
| One warp in spike | `grep -v '^#' spike/simulator/forward.jl \| grep -c 'warp('` | `1` |
| Old translation call gone | `grep -v '^ *#' … \| grep -c 'Translation(θ.shift_dy, θ.shift_dx), axes(ch2)'` | `0` |
| `recenter` present / no LinearAlgebra | `grep -c recenter` / `grep -c LinearAlgebra` | `4` / `0` |
| One warp in src | `grep -v '^ *#' src/amortized/simulator.jl \| grep -c 'warp('` | `1` |
| src shift prior untouched (Q2) | `grep -c 'SHIFT_PRIOR *= *Uniform(-1.0, 1.0)' src/amortized/simulator.jl` | `1` |
| datagen literals gone | `grep -c 'undef, 7,'` / `grep -c 'theta_dim *= *7'` | `0` / `0` |
| `NPE_D` untouched | `grep -n NPE_D src/amortized/architecture.jl` | `const NPE_D          = 7` (line 169) |
| `architecture.jl` not in diff | `git diff --name-only -- src/amortized/architecture.jl` | empty |
| Only `src/amortized/` touched (Q8) | `git diff --name-only 17ebd1e… -- src/ \| grep -v '^src/amortized/'` | empty |
| Manuscript pipelines byte-unchanged | `git diff --stat 17ebd1e… -- src/bayes.jl src/colocalization.jl src/LoadImages.jl src/plot.jl` | empty |
| Artifact + fixtures byte-unchanged | `git diff --stat 17ebd1e… -- Artifacts.toml test/test_images/` | empty |
| No dependency change | `git diff --stat -- Project.toml Manifest.toml spike/Project.toml spike/Manifest.toml` | empty |
| No deletions in the commit | `git diff --diff-filter=D --name-only HEAD~1 HEAD` | empty |
| No orphaned cache pruned | `git status --porcelain -- artifacts/` | empty (`artifacts/` is gitignored; nothing deleted) |

**Not verified here (out of this plan's scope):** nothing in this plan trains a net, runs the D-06
probe, runs an SC2 ladder or evaluates coverage, so **no calibration, width, monotonicity or
identifiability claim is made or checked**. In particular the new `chromatic_eps` column is
*expected* to be vacuous (D-05) and this plan produced no evidence either way. The shipped
`amended_v2/grid_8` bundle was neither retrained nor re-gated.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] `train_npe` was left internally inconsistent by the θ-arity change**
- **Found during:** Task 3 (first attempt to update the root suite)
- **Issue:** The plan states "do NOT change `src/amortized/train_npe.jl:136`". That instruction
  rests on a premise the change invalidates. `fit_theta_transform` defaults `bounds =
  theta_prior_bounds()` and hard-requires `size(θtr, 1) == length(bounds)`, so after the extension
  it demands **8** θ rows — while `build_estimator(d_in; …)` at `:136` still built a **7**-marginal
  flow from `NPE_D`. Measured directly: 7 rows → `DimensionMismatch: θ has 7 rows but the prior box
  has 8`; 8 rows → `AssertionError: d == flow.d`. **Every** input threw. The whole `src` training
  path (and therefore `_train_grid_pipeline`) was dead, and `datagen` now emits 8-row θ, so it was
  guaranteed to be hit.
- **Fix:** Derived the flow's marginal count from the θ actually passed in —
  `D_theta = size(θtr, 1)`, forwarded to `build_estimator(d_in, D_theta; …)` — and made `arch.D`
  record `D_theta` rather than `NPE_D`. The second half is not cosmetic: `persist.load_estimator`
  rebuilds from `arch.D`, so a stale `arch.D` would have persisted a net that cannot be reloaded.
  **`NPE_D` (still 7), `build_estimator`'s default `D`, and `architecture.jl` are untouched**, so
  the shipped bundle — whose `D` is read off disk — is unaffected. Both are asserted.
- **Files modified:** `src/amortized/train_npe.jl` (2 code lines + comments/docstring)
- **Verification:** the round-trip is now asserted in two places —
  `@test res.arch.D == Dθ` (`test/runtests.jl`) and the save→load→`posterior_for` assertion in
  `test/gpu_smoke.jl`, which caught the stale-`arch.D` case empirically before the fix.
- **Deviates from:** the plan's explicit "do NOT change `train_npe.jl:136`" and its acceptance
  criterion `git diff --name-only -- src/amortized/architecture.jl src/amortized/train_npe.jl` is
  empty. `architecture.jl` **is** empty in the diff; `train_npe.jl` is not. Flagged rather than
  silently satisfied.
- **Committed in:** `ca02b0e`

**2. [Rule 3 - Blocking] Three θ-arity sites outside the plan's enumeration**
- **Found during:** Task 3 (suite runs 1 and 2)
- **Issue:** `11-RESEARCH.md` §A4 row 24 enumerates the literal-7 sites in `test/runtests.jl` only;
  its "hard-coded 7 audit" grepped `src/` and `test/` but the plan's file list carried only
  `test/runtests.jl`. Three further files feed θ to `train_npe` and broke the suite:
  `test/gpu_smoke.jl:38` (`DimensionMismatch`, aborted the run), `test/gpu_smoke.jl:63` (the
  save→load round-trip, which correctly reported 8 vs 7), and `test/test_local_map.jl:30,177`.
- **Fix:** Replaced each with `length(ProteinCoLoc.theta_prior_bounds())`, matching the
  derive-never-literal rule used everywhere else.
- **Files modified:** `test/gpu_smoke.jl`, `test/test_local_map.jl`
- **Verification:** root suite exit 0 on the third run.
- **Committed in:** `ca02b0e`

**3. [Rule 3 - Blocking] `spike/test/test_p11_forced_theta.jl` pinned the PRE-edit θ arity**
- **Found during:** Task 3 (plan verification item 2)
- **Issue:** Plan verification asserts this file "still exit 0 after the arity change". It did not
  (exit 1, 4 failures) — by 11-01's own design. The file states in-line: *"TODAY the last field is
  `noise`. Asserting that here makes the append-at-end contract for `chromatic_eps` checkable
  BEFORE and AFTER the θ-arity edit"*. It is a declared tripwire that this edit was supposed to
  trip, not a frozen pre-registration (that is `test_p11_consts.jl`, which stayed green).
- **Fix:** Updated it to exactly the post-edit form it documents — `last(keys(θ)) ==
  :chromatic_eps`, rows 1..7 asserted unchanged — and **kept it a real tripwire** rather than a
  restatement of the prior's current shape by adding a synthetic in-place-overwrite case
  (`merge(θ, (chromatic_eps = 0.017,))` must not grow or reorder) and a synthetic new-key case
  (`:future_param` must append at the end).
- **Files modified:** `spike/test/test_p11_forced_theta.jl`
- **Verification:** `julia --project=spike spike/test/test_p11_forced_theta.jl` — **28 pass / 0
  fail**, exit 0 (was 17 pass / 4 fail).
- **Committed in:** `ca02b0e`

**4. [Rule 2 - Missing critical] The pinned sha: fixture anchor vs first parent**
- **Found during:** Task 1 step 1
- **Issue:** Two authorities disagreed. The plan requires `PRE_EPS_SHA = git rev-parse HEAD` at
  task start, with the acceptance criterion that it equals the ε commit's **first parent**
  (`a2ced26`). The orchestrator's non-negotiable requires the sha to be `17ebd1e`, the
  fixture-embedded `PHASE11_BASE_SHA`, verified against the fixture first.
- **Fix:** Verified the fixture reads `17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348`, verified
  `src/amortized/simulator.jl` is **byte-identical** between `17ebd1e` and `a2ced26`
  (`git diff --stat` empty), and quoted **both** in limit #8 with that byte-identity stated
  explicitly. The `git show <SHA>:…` reconstruction is therefore correct against either, and both
  acceptance readings are satisfied simultaneously.
- **Files modified:** `docs/amortized.md`
- **Verification:** the parent-relationship table above.
- **Committed in:** `ca02b0e`

**5. [Rule 1 - Bug] Two plan acceptance criteria are unsatisfiable as literally written**
- **Found during:** Task 2 and Task 3
- **Issue:** (a) `grep -c 'NPE_D = 7' src/amortized/architecture.jl` outputs `0`, not `1` — the
  file uses column-aligned spacing (`const NPE_D          = 7`), so the single-space pattern cannot
  match. (b) The criterion "`git show --stat HEAD --name-only` lists exactly these seven paths"
  cannot hold given deviations 1-3. (c) `test/runtests.jl:495` is listed as a site to convert to
  the θ arity, but it asserts the output of `build_estimator(d_in; …)` called **without** `D`,
  which uses `NPE_D = 7` by design — converting it to 8 would have made it fail.
- **Fix:** (a) Verified the substantive claim instead: `NPE_D` is line 169, unchanged, and
  `architecture.jl` is absent from the diff. (b) The commit carries twelve paths, listed above.
  (c) Left as 7 but rewritten to assert `ProteinCoLoc.NPE_D` with a comment naming why, so it can
  never be "corrected" to the θ arity by a later editor.
- **Files modified:** `test/runtests.jl`
- **Verification:** root suite exit 0; `git diff --name-only -- src/amortized/architecture.jl`
  empty.
- **Committed in:** `ca02b0e`

---

**Total deviations:** 5 auto-fixed (1 bug, 3 blocking, 1 missing-critical). **No Rule 4
(architectural) situation arose; no checkpoint was hit.**
**Impact on scope:** Deviation 1 is the only one that changes behaviour in `src/`, and it is
strictly a repair of a break this plan introduced; it does not touch the shipped read surface, the
public API, `NPE_D`, `Artifacts.toml` or the artifact. Deviations 2, 3 and 5 are test-side. All are
inside the D-16 boundary (no API change, no artifact change).

## Assumption Drift (advisory)

**1. "Only `test/runtests.jl` couples to θ arity in the test tree."**
- **Found during:** Task 3
- **Planned:** `11-RESEARCH.md` §A4 row 24's list of 17 line numbers in `test/runtests.jl`, cited
  as the complete blast radius, with "there is no other literal-7 θ coupling in `src/`".
- **Actual:** `test/gpu_smoke.jl` and `test/test_local_map.jl` also construct θ for `train_npe`,
  and `scripts/diag_bounded_theta_dev.jl` hardcodes 7 in four places.
- **Why it matters:** a reader auditing the ripple against §A4 alone would conclude it is closed.
  It is closed for the *suite*; `scripts/` is logged in `deferred-items.md`, not fixed.

**2. "`src/`'s training path was going to keep working."**
- **Found during:** Task 3
- **Planned:** the plan treats `train_npe` as out of scope and untouched, implying the src-side
  trainer is simply "not retargeted at 8-θ".
- **Actual:** after the extension the src trainer accepted **no** θ arity at all — not 7, not 8.
  "Do not retarget it" and "leave it working" were not compatible.
- **Why it matters:** any statement that Phase 11 left `src/` training untouched must be read as
  "left it *derived*". `build_estimator`'s default `D` — the thing the shipped bundle actually
  depends on — genuinely is untouched, and that is the claim worth making.

**3. "The 7 vs 8 distinction is one number."**
- **Found during:** Tasks 2-3
- **Planned:** implicitly, a single θ arity moving 7 → 8.
- **Actual:** three distinct quantities that merely coincided at 7 — the **prior's** θ arity (now
  8), the **shipped bundle's** frozen flow marginal count (still 7, read off disk), and
  `NPE_D`/`build_estimator`'s default (still 7). `test/test_integration.jl` and
  `test/runtests.jl:495` assert the latter two and must NOT move.
- **Why it matters:** the Phase-11 report should not describe "θ is now 8-dimensional" without
  saying that the shipped estimator's posterior is still 7-row — that separation is precisely what
  makes D-16's "no artifact change" true.

## Issues Encountered

- **`train_npe` threw on every θ arity** after the prior extension (both failure modes reproduced
  by direct execution before any fix). Resolved as deviation 1.
- **The suite needed three runs** to reach green: run 1 errored in `test/gpu_smoke.jl:34`
  (`DimensionMismatch`), run 2 failed at `test/gpu_smoke.jl:63` (save→load round-trip reporting
  8 vs 7 — which is the *correct* new behaviour and confirmed the `arch.D` half of deviation 1 was
  necessary), run 3 passed. A pre-edit baseline run was taken first, so the green→green transition
  is established rather than assumed.
- **`spike/test/test_p11_forced_theta.jl` was red by design** after the edit; the file's own
  comments specified the post-edit form. Resolved as deviation 3.

## Known Stubs

None. Every field, guard and buffer introduced is fully wired: `chromatic_eps` is drawn from its
prior, validated at entry, and reaches stage 6 in both simulators; `DATAGEN_THETA_DIM` is consumed
at all three buffer sites and by the hashed config; `_theta_tuple`'s 8th field is exercised in the
new testset in all four branches (short vector, in-range, out-of-range clamp, non-finite).

No placeholder value, empty container, mock data path or TODO was introduced.

**Intentionally deferred, not stubbed** (logged in `deferred-items.md`): `test/gate/sbc.jl`'s
`SBC_PARAM_LABELS` remains 8 labels (7 θ + Δρ) and `sbc_ranks_and_spread` still loops `for p in
1:7`. This is `11-RESEARCH.md` §A4 row 25's explicit deferral — the src ship gate is not re-run in
Phase 11 (D-01) and the frozen pre-registration must not be edited. The loop is arity-safe: it
reads the first 7 rows of an 8-row draw matrix, which is why the gate testsets stayed green.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no new file-access pattern and no schema at
a trust boundary. The registered threats were handled as planned:

| Threat ID | Handling | Evidence |
|---|---|---|
| T-11-07 (`chromatic_eps` reaching stage 6) | `isfinite` joined the existing wall and a new `-1 < chromatic_eps` bound added, in **both** `simulate_pair` copies | `@test_throws ArgumentError` over `(NaN, Inf, -1.0, -2.5)` in the new testset |
| T-11-08 (7-row vector into an 8-field `simulate_pair`) | `length(v) >= 8` guard in `_theta_tuple` + `hasproperty` read in both simulators | legacy-7-field θ asserted bit-identical to the 8-field ε=0 call |
| T-11-09 (training provenance severed) | Named limit #8 landed in the SAME commit, pinning two full 40-char shas | `git rev-parse HEAD^` = `a2ced26`; both shas `grep -c chromatic_eps` = `0` |
| T-11-10 (orphaned caches deleted) | Do-not-prune recorded in `open_or_invalidate`'s docstring and in limit #8; nothing deletes | no deletions in the commit; `artifacts/` untouched |
| T-11-11 (artifact integrity) | `Artifacts.toml` and `registry.jl` untouched | `git diff --stat 17ebd1e… -- Artifacts.toml` empty |
| T-11-12 (stale cache read as current) | Unchanged digest-mismatch error; the double flip (source bytes **and** `theta_dim`) auto-separates the new pool | `theta_dim` is now `DATAGEN_THETA_DIM`, a hashed field |
| T-11-SC (package installs) | None attempted; `LinearAlgebra` avoided via the plain `[s 0.0; 0.0 s]` matrix | `grep -c LinearAlgebra` = `0`; all four manifests byte-unchanged |

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- **Ready:** plan 11-03 can formalize the stage-6 regression as a committed test file. The
  regression itself is already **measured green** here (exact `==` against all four golden keys,
  both θ arities), so 11-03 is codifying a verified result, not discovering one.
- **Ready:** the D-06 probe can now sweep `chromatic_eps` through `forward.jl` at fixed θ.
- **Ready:** D-15 simulator-injected misalignment works — `merge(θ, (shift_dx = …, shift_dy = …,
  chromatic_eps = …))` overwrites in place without reordering, asserted in the spike forced-θ test.
- **Constraint for the research trainer (11-07):** it is spike-local. The λ-conditioned 129-row
  input surface and any `D = 8` flow construction must NOT cross into `src/`; `NPE_D` stays 7 and
  `build_estimator`'s default stays `NPE_D`.
- **Constraint for the report:** the vacuity of `chromatic_eps` and the widened `shift_*` columns
  is **unmeasured** by this plan and must not be asserted from it in either direction.

## Self-Check: PASSED

- `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/deferred-items.md` exists
  on disk.
- All twelve modified paths appear in `git show --name-only ca02b0e`.
- Commit `ca02b0e8a696269a2b8710b11ca3073edc5ae6c7` exists in `git log`, and its first parent is
  `a2ced26d93c7409a465aa50320e188868b3ab462` — the sha quoted in `docs/amortized.md`.
- `docs/amortized.md` limit `8.` sits at line 152, after limit `7.` (line 145) and before
  `**Findings and scripts:**` (line 172).

---
*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Completed: 2026-07-25*
