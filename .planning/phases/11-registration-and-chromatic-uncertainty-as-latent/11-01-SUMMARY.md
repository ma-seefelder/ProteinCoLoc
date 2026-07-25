---
phase: 11-registration-and-chromatic-uncertainty-as-latent
plan: 01
subsystem: testing
tags: [julia, random123, philox, pre-registration, tost, holm, wilson, jld2, sbc]

# Dependency graph
requires:
  - phase: 07-productionization-conditional-on-go
    provides: "the frozen amended pre-registration test/gate/gate_consts_8_v2.jl (F5 imsize mixture, holm_adjusted, PROD_SEED/PROD_SEED_V2) and the vacuous-column reporting discipline"
  - phase: 05-validation-bundle-sbc-amortized-bf-ood
    provides: "spike/validation/consts.jl (the anti-snooping consts file shape, VAL_MASTER_SEED/VAL_FIX_SEED) and the spike test/harness idioms"
provides:
  - "spike/validation/p11_consts.jl — the Tier-1 anti-snooping contract: DEV seed/salt, forbidden-seed families, lambda and SC2/SC3 ladders, TOST band, probe spec, result-immune abort criterion, F5 imsize read-by-reference, datagen ceiling, stage-6 regression contract, one-iteration allowance, four declared net deviations"
  - "spike/validation/p11_stats.jl — wilson_ci / tost_pvalue / p11_holm_adjusted / p11_tost_pass (inverted direction) / p11_perm_spearman / p11_shrinkage / breakdown_point"
  - "spike/test/fixtures/p11_stage6_golden.jld2 — the PRE-EDIT simulate_pair golden bytes plus PHASE11_BASE_SHA"
  - "PHASE11_BASE_SHA = 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348 (the decoupling/provenance diff anchor)"
affects: [11-02-simulator-chromatic-term, 11-03-stage6-regression, 11-05-probe-tier2-append, 11-ladder-and-breakdown, 11-report]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Two-tier append-only pre-registration (Tier 1 committed before the probe; Tier 2 appended, never edited)"
    - "Isolated `module GateV2` read of the frozen gate constants, declared at file top level because a module expression cannot sit inside a guard block"
    - "Copy-with-attribution of gate helpers that cannot be depended on from the spike environment"

key-files:
  created:
    - spike/validation/p11_consts.jl
    - spike/validation/p11_stats.jl
    - spike/test/test_p11_consts.jl
    - spike/test/test_p11_tost.jl
    - spike/test/test_p11_forced_theta.jl
    - spike/test/capture_p11_golden.jl
    - spike/test/fixtures/p11_stage6_golden.jld2
  modified: []

key-decisions:
  - "P11_DEV_SEED = 0x0000_0000_0B11_DE71 with P11_SALT = 0xA24B_AED4_663E_E121, executably proven disjoint from 22 forbidden seeds including both recomputed ship-gate families"
  - "The forbidden-seed recomputation is cross-checked against GateV2.PROD_SEED / GateV2.PROD_SEED_V2, proving the single-draw derivation reproduces the gate's redraw-loop rule"
  - "The Holm-agreement test asserts against GateV2.holm_adjusted rather than test/gate/sbc.jl, because sbc.jl cannot load standalone under --project=spike (measured)"
  - "module GateV2 is declared outside the guard block: Julia rejects a module expression that is not at top level"
  - "P11_VACUOUS_SHRINKAGE lives in p11_stats.jl, not in the pre-registration, because it is a reporting label and not a gate input"

patterns-established:
  - "Locked constants are asserted AS LITERALS in a companion test, so a post-hoc edit breaks the suite loudly"
  - "A golden fixture embeds its own capture-time commit sha and refuses to be written against a dirty tree for the files it precedes"
  - "Statistics whose comparison direction can be silently inverted get a named direction-trap comment plus a hand-computed fixture vector"

requirements-completed: [D-01, D-04, D-07, D-08, D-10, D-14]

# Metrics
duration: 17min
completed: 2026-07-25
---

# Phase 11 Plan 01: Pre-registration, Equivalence Statistics and the Pre-Edit Golden Summary

**The Phase-11 anti-snooping contract is committed and executably self-checked: a fresh DEV seed proven disjoint from 22 forbidden streams, the SC2/SC3 ladders and the ±3 pp TOST band with their derivations recorded, an equivalence Holm whose `<=` direction is pinned by a hand-computed vector, and the pre-edit `simulate_pair` golden bytes captured at `PHASE11_BASE_SHA` before a single simulator byte moved.**

## PHASE11_BASE_SHA

```
17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348
```

Short form: `17ebd1e`.

This is the commit `HEAD` pointed at when `spike/test/capture_p11_golden.jl` ran, verified at
capture time to have a clean working tree for `spike/simulator/` and `src/amortized/simulator.jl`.
It is stored inside the fixture as `phase11_base_sha` / `phase11_base_sha_short`, and it is the
anchor the decoupling and provenance diffs consume:

```bash
git diff --stat 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- \
    src/bayes.jl src/colocalization.jl src/LoadImages.jl src/plot.jl
git diff --name-only 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- src/ | grep -v '^src/amortized/'
git diff --stat 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- test/test_images/
git diff --stat 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- Artifacts.toml
```

## Performance

- **Duration:** ~17 min
- **Started:** 2026-07-25T18:44Z (worktree base `e645cde`)
- **Completed:** 2026-07-25T19:01Z
- **Tasks:** 3
- **Files created:** 7 (no file modified, no file deleted)

## Accomplishments

- **Tier-1 pre-registration is committed before anything runs.** `spike/validation/p11_consts.jl`
  (379 lines, exactly one guard block) locks the λ range, the seven-rung SC2 ladder, both
  beyond-prior ladders, `N_PER_RUNG = 500` with its `N_min = 271` derivation, the ±3 pp TOST band
  with the coverage-to-width algebra recorded verbatim, the D-06 probe specification, and the
  abort criterion that must outrank the probe's own result.
- **Seed disjointness is executable, not prose.** `_p11_forbidden()` spans 22 seeds — the four
  named reserved streams, the corpus stream, eight prior DEV/burned keys, and all eight
  ship-gate seeds, which are **recomputed** (never trusted from a comment) and then
  cross-checked against the frozen gate file's own `PROD_SEED` / `PROD_SEED_V2`. The
  cross-check also proves the single-draw derivation reproduces the gate's redraw-loop rule,
  so nothing is silently missing from the forbidden set.
- **The Holm direction trap is closed.** `p11_tost_pass` uses `all(<=(fwer), ...)`, the mirror
  image of `sbc_holm_pass`'s `all(>(fwer), ...)`. The test asserts, on the same hand-computed
  vector `[0.001, 0.002, 0.04] → [0.003, 0.004, 0.04]`, that the equivalence reading passes and
  the SBC reading fails — so the two directions are shown to genuinely disagree there, not merely
  to be differently spelled. The name `holm_pass` is never reused in the spike file.
- **The F5 mixture is read, never retyped.** `P11_IMSIZE_SET`/`P11_IMSIZE_WEIGHTS` bind
  `GateV2.SBC_IMSIZE_SET`/`WEIGHTS` through the isolated-module idiom; the consts file contains
  zero occurrences of the anchor size literal, and the test asserts identity (`===`) with the
  gate constants rather than equality with a retyped copy.
- **The pre-edit golden exists and is self-dating.** The fixture holds 4 × 2 `256²` channel
  matrices from an explicit 7-field θ literal, the DEV seed/salt, the fixture counter, a UTC
  `generated` stamp, and the capture-time sha. The capture script refuses to run on a dirty tree
  for the two simulator files it claims to precede.

## Task Commits

1. **Task 1: Tier-1 pre-registration consts + literal-value test** — `d336699` (feat)
2. **Task 2: equivalence statistics with the inverted Holm direction** — `17ebd1e` (feat)
3. **Task 3: pre-edit stage-6 golden fixture + base sha capture** — `9392897` (test)

## Files Created

- `spike/validation/p11_consts.jl` — the Tier-1 anti-snooping contract (379 lines, one guard block).
- `spike/validation/p11_stats.jl` — Wilson CI, TOST, Holm (copied, then inverted), rung-label
  permutation Spearman, shrinkage vacuity diagnostic, breakdown point.
- `spike/test/test_p11_consts.jl` — 84 assertions pinning every locked value as a literal.
- `spike/test/test_p11_tost.jl` — 36 passing assertions + 1 deliberate skip; every statistic is
  exercised in both directions.
- `spike/test/test_p11_forced_theta.jl` — 21 exact (no-tolerance) assertions on the `merge`
  field-order and θ row-alignment invariants.
- `spike/test/capture_p11_golden.jl` — the one-shot pre-edit capture script.
- `spike/test/fixtures/p11_stage6_golden.jld2` — the golden bytes + `phase11_base_sha`.

No file was modified and no file was deleted. `spike/Project.toml`, `spike/Manifest.toml` and all
of `src/` are byte-unchanged since the worktree base (`git diff --stat` empty). No package was
installed.

## Verification — observed results

Every command below was run and its real output observed.

| Check | Command | Observed |
|---|---|---|
| Tier-1 consts gate | `julia --project=spike spike/test/test_p11_consts.jl` | **84 pass / 0 fail**, exit 0 |
| Equivalence statistics | `julia --project=spike spike/test/test_p11_tost.jl` | **36 pass / 0 fail / 1 broken (deliberate `@test_skip`)**, exit 0 |
| θ row alignment | `julia --project=spike spike/test/test_p11_forced_theta.jl` | **21 pass / 0 fail**, exit 0 |
| Golden fixture readback | `JLD2.load(...)`, 4 channels, 40-char sha, no `:chromatic_eps` | prints `ok` |
| Fixture is tracked | `git ls-files --error-unmatch spike/test/fixtures/p11_stage6_golden.jld2` | exit 0 |
| No dependency change | `git diff --stat -- spike/Project.toml spike/Manifest.toml` | empty |
| No `src/` change | `git diff --stat e645cde..HEAD -- src/` | empty |
| One guard block | `grep -c 'if !isdefined' spike/validation/p11_consts.jl` | `1` |
| Mixture read, not retyped | `grep -c '(1376, 1028)' spike/validation/p11_consts.jl` | `0`; `module GateV2` present |
| Inverted direction present | `grep -c 'all(<=(' spike/validation/p11_stats.jl` | `2` |
| `holm_pass` never reused | `grep -c 'function holm_pass\|^holm_pass' spike/validation/p11_stats.jl` | `0` |
| Attribution on copied fns | `grep -c 'test/gate/sbc.jl' spike/validation/p11_stats.jl` | `5` |
| No tolerance in row-order test | `grep -c 'atol\|isapprox' spike/test/test_p11_forced_theta.jl` | `0` |
| Fixture θ is a literal | `grep -c 'sample_prior' spike/test/capture_p11_golden.jl` | `0` |
| Clean tree at capture | `git diff --stat -- spike/simulator/ src/amortized/simulator.jl` | empty |

**Not verified here (out of this plan's scope):** nothing in this plan trains, simulates a ladder,
or evaluates a net, so no calibration, coverage or width claim is made or checked. The
`SC2_SPEARMAN_FLOOR` and the `Δρ_eq` calibration coefficients are deliberately absent — the test
asserts `!isdefined(@__MODULE__, :SC2_SPEARMAN_FLOOR)` to prove Tier 2 has not been pre-empted.

## Decisions Made

1. **`module GateV2` is declared at file top level, outside the guard block.** Julia rejects a
   `module` expression that is not at top level (`syntax: "module" expression not at top level` —
   reproduced directly). Idempotency is instead carried by the guarded-include idiom every caller
   already uses. This is documented in-file so it does not read as an oversight.
2. **The forbidden-seed recomputation is cross-checked against the frozen file's own values.**
   `@assert PROD_SEED_V2 == GateV2.PROD_SEED_V2` (and the v1 twin) turns "we recomputed the derived
   gate seeds" from a claim into a proof, and incidentally proves the single-draw `_derive` matches
   the gate's redraw-loop rule.
3. **`P11_VACUOUS_SHRINKAGE = 0.95` lives in `p11_stats.jl`, not in the pre-registration.** It is a
   reporting label, and putting a reporting label in the frozen consts invites it being read as a
   pass condition — the same reasoning `test/gate/sbc.jl:226-229` gives.
4. **`p11_perm_spearman` takes `rng` as a required keyword** (no default), so no caller can
   accidentally ride an implicit global stream.
5. **`P11_PERM_ALPHA = 0.05` was added to Tier 1** beyond the plan's enumeration: SC2(a) is gated on
   `p_perm ≤ α` and that α is exactly the kind of threshold D-04 requires to be locked in advance.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] `module GateV2` moved outside the guard block**
- **Found during:** Task 1 (Tier-1 consts file)
- **Issue:** The plan specifies the isolated-module read *inside* the `if !isdefined(...)` Tier-1
  block. Julia rejects that: `syntax: "module" expression not at top level` (reproduced with a
  minimal `julia -e` case before writing the file).
- **Fix:** Declared `module GateV2` at file top level immediately above the guard block, with an
  in-file comment naming the language restriction and the guarded-include idiom that supplies
  idempotency instead. The acceptance criteria are unaffected — one `if !isdefined`, `module GateV2`
  present, anchor literal absent.
- **Files modified:** `spike/validation/p11_consts.jl`
- **Verification:** file loads and all self-asserts pass; `grep -c 'if !isdefined'` = 1.
- **Committed in:** `d336699`

**2. [Rule 2 - Missing critical] The Holm-agreement test asserts against a frozen adjuster that actually loads**
- **Found during:** Task 2 (equivalence statistics)
- **Issue:** The plan directs the agreement testset to load `test/gate/sbc.jl` into an isolated
  `module GateSBC`, falling back to `@test_skip` if it cannot load. It cannot: it transitively
  includes `test/gate/harness.jl:35`, which does `using ProteinCoLoc`, and that package does not
  resolve under `--project=spike` (measured — `ArgumentError: Package ProteinCoLoc not found in
  current path`). A bare skip would have left the copied Holm with **no** agreement evidence.
- **Fix:** Kept the documented `@test_skip` naming the unavailable route, and additionally asserted
  agreement against `GateV2.holm_adjusted` (`test/gate/gate_consts_8_v2.jl:361-371`) — the
  authoritative frozen implementation that `sbc_holm` itself defers to, already loaded for the
  imsize mixture, and loadable under `--project=spike`. Also asserted that `GateV2.holm_pass`
  returns the *opposite* verdict on the hand-computed vector, which strengthens the direction claim
  against the real frozen function rather than against a local restatement.
- **Files modified:** `spike/test/test_p11_tost.jl`, `spike/validation/p11_stats.jl` (attribution
  comment records the measured load failure)
- **Verification:** `julia --project=spike spike/test/test_p11_tost.jl` — 36 pass, 1 broken, exit 0.
- **Committed in:** `17ebd1e`

**3. [Rule 2 - Missing critical] `P11_PERM_ALPHA` added to Tier 1**
- **Found during:** Task 1 (Tier-1 consts file)
- **Issue:** The plan locks `PERM_B` and the Tier-2 Spearman floor but not the permutation test's
  own α, which is a pass/fail threshold for SC2(a) and therefore squarely inside D-04's scope. Left
  unlocked it could have been chosen after seeing `p_perm`.
- **Fix:** Added `const P11_PERM_ALPHA = 0.05` with its rationale comment; asserted in the test.
- **Files modified:** `spike/validation/p11_consts.jl`, `spike/test/test_p11_consts.jl`
- **Verification:** consts test passes; the value is asserted as a literal.
- **Committed in:** `d336699`

---

**Total deviations:** 3 auto-fixed (1 blocking, 2 missing-critical). **No Rule 4 (architectural)
situation arose; no checkpoint was hit.**
**Impact on plan:** None on scope. Deviation 1 is a language constraint with an equivalent
implementation. Deviations 2 and 3 both *strengthen* the pre-registration rather than relax it —
neither adds a dependency, changes a locked value, or touches `src/`.

## Assumption Drift (advisory)

**1. "The plan's Tier-1 enumeration is the complete set of thresholds SC2 gates on."**
- **Found during:** Task 1
- **Planned:** the enumerated Tier-1 list (band, α per side, FWER, `PERM_B`, ladders, N).
- **Actual:** SC2(a)'s permutation-test α was not in the list although it is a gating threshold.
- **Why it matters:** a reader auditing D-04 compliance against the plan's enumeration alone would
  have found a gap. Recorded as deviation 3 above; the consts file is now the complete set.

**2. "`test/gate/sbc.jl` is the frozen Holm source."**
- **Found during:** Task 2
- **Planned:** copy from and agree with `test/gate/sbc.jl:352-362`.
- **Actual:** that file's own docstring says it is the *fallback*; the frozen implementation is
  `holm_adjusted` in `gate_consts_8_v2.jl`, and only the latter is reachable from the spike
  environment. The copied bytes are identical in both, so the copy itself is unaffected — only the
  attribution and the agreement target changed.
- **Why it matters:** the report's provenance sentence for the copied helper should name
  `gate_consts_8_v2.jl` as the authority, with `sbc.jl` as the byte-identical fallback.

## Issues Encountered

- **`test/gate/sbc.jl` is not loadable from the spike environment.** Diagnosed by direct execution
  (`ArgumentError: Package ProteinCoLoc not found in current path`, raised at
  `test/gate/harness.jl:35`). Resolved as deviation 2 — no dependency was added to force the load,
  which would have violated the spike decoupling constraint.
- **`grep -c 'if !isdefined'` initially returned 2** because a prose comment quoted the guard form.
  The comment was reworded; no code changed.

## Known Stubs

None. Every function in `p11_stats.jl` is fully implemented and exercised in both directions;
no placeholder value, empty container or TODO was introduced. `SC2_SPEARMAN_FLOOR` and the
`Δρ_eq` coefficients are **absent by design**, not stubbed — they are Tier-2 constants that must
not exist until the probe has run, and their absence is asserted.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no file-access pattern beyond a single
`spike/test/fixtures/` write, and no schema at a trust boundary. The registered threats were
handled as planned:

| Threat ID | Handling |
|---|---|
| T-11-01 (snooped seed) | `@assert !(P11_DEV_SEED in _p11_forbidden())` over 22 seeds, re-asserted pairwise in the test |
| T-11-02 (post-hoc threshold edit) | every locked value asserted as a literal; Tier 2 is a future append, never an edit |
| T-11-03 (post-edit fixture) | `phase11_base_sha` embedded in the JLD2; capture refuses on a dirty simulator tree |
| T-11-04 (torn write) | `.tmp` → reopen-and-`@assert haskey` → `mv(...; force = true)` |
| T-11-05 (retyped mixture) | read via `GateV2`; anchor literal count in the consts file is 0; test asserts `===` identity |
| T-11-SC (package installs) | none attempted; `spike/Project.toml` and `spike/Manifest.toml` byte-unchanged |

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- **Ready:** the simulator edit may now proceed. `PHASE11_BASE_SHA` is recorded, the pre-edit golden
  exists and is tracked, and the θ row-alignment invariants the append-at-end layout depends on are
  asserted and currently green (`last(keys(θ)) == :noise`).
- **Ready:** the probe may consume `p11_rng(P11_PROBE_COUNTER)` against a locked specification and a
  locked abort criterion.
- **Sequencing constraint for the Tier-2 append:** it must be a **second guard block keyed on
  `SC2_SPEARMAN_FLOOR`**, appended to `p11_consts.jl`, never an edit to the Tier-1 block. The
  Tier-1 test asserts `SC2_SPEARMAN_FLOOR` is currently undefined; that assertion must be removed
  in the same commit that appends Tier 2, and only then.
- **Not started, by design:** no training, no ladder, no evaluation. `P11_ITERATION_ALLOWANCE = 1`
  is unspent.

## Self-Check: PASSED

- All 7 created artifact paths exist on disk (`ls -1` returned every one).
- All 3 task commits exist in `git log --oneline --all`: `d336699`, `17ebd1e`, `9392897`.
- The `PHASE11_BASE_SHA` recorded above matches `phase11_base_sha` inside the committed fixture.

---
*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Completed: 2026-07-25*
