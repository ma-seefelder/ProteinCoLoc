---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 01
subsystem: pre-registration (spike-local, Phase-13 Tier 1)
tags: [pre-registration, anti-snooping, seeds, thresholds, d-04, d-01, spike]
requires:
  - "test/gate/gate_consts_8_v2.jl (READ-ONLY: recomputed PROD_SEED / PROD_SEED_V2, F5 imsize mixture)"
  - "Random123 (already resolved in spike/Project.toml; no dependency added)"
provides:
  - "spike/p13/consts.jl — the Phase-13 Tier-1 pre-registration (89 locked constants + 2 RNG constructors + _p13_forbidden())"
  - "p13_rng(counter) / p13_fix_rng(counter) — the reported and fixture Philox streams"
  - "_p13_forbidden() / _p13_forbidden_list() — the 41-entry executable forbidden-seed proof"
  - "p13_tau() — the Tier-2 tripwire that fails loudly until plan 13-10 appends the measured tau"
  - "spike/test/test_p13_consts.jl — 132 literal-value assertions guarding the pre-registration"
affects:
  - "every later Phase-13 plan (13-02 .. 13-16) reads its thresholds, ladders and seeds from spike/p13/consts.jl"
  - "plan 13-10 must APPEND a second guard block (P13_TAU) and RELOCATE the Tier-2-absent assertion"
tech-stack:
  added: []
  patterns:
    - "ONE guard block keyed on a phase sentinel (spike/validation/consts.jl:42)"
    - "isolated read-only `module _GC` include of a frozen pre-registration, at top level because Julia rejects a module expression inside an `if` (spike/validation/p11_consts.jl:66-77)"
    - "recompute-don't-trust-comments for DERIVED seeds (test/gate/gate_consts_8_v2.jl:263-331)"
    - "salt-XOR-into-the-first-Philox-key-word (spike/data/seeding.jl:75-96)"
    - "assert-the-pre-registration-is-locked testset with literal values (spike/test/test_bf.jl:91-101)"
    - "comment-stripped source-grep gate (spike/test/test_bf.jl:67-68, 83-89)"
key-files:
  created:
    - spike/p13/consts.jl
    - spike/test/test_p13_consts.jl
  modified: []
decisions:
  - "P13_DEV_SEED = 0x0000_0000_0B13_DE71 and P13_FIX_SEED = 0x0000_0000_0B13_F1F7, executably proven disjoint from a 41-entry forbidden set that includes the RECOMPUTED PROD_SEED and PROD_SEED_V2 families and P11_DEV_SEED (forbidden before Phase 11 has even run)"
  - "P13_TAU is deliberately absent (Tier 2); p13_tau() errors with a pointer to the post-Phase-11 probe, so labelled data cannot be generated before tau is measured (D-06)"
  - "The isolated gate-consts module sits at file top level, not inside the guard block — Julia rejects a `module` expression that is not at top level; idempotency comes from the callers' guarded-include idiom (matches spike/validation/p11_consts.jl)"
  - "The not-a-gate source assertion uses a word-boundary-anchored name pattern and asserts SET EQUALITY with the two documented exemptions, because the unanchored pattern matches the 'MIN' substring inside the plan-mandated P13_REAL_NAMING_CORRECTION"
metrics:
  duration: ~35 min
  completed: 2026-07-27
  tasks: 2
  files: 2
  commits: 2
  tests: 132 pass / 0 fail / 0 error
---

# Phase 13 Plan 01: Tier-1 Pre-Registration Summary

Locked the Phase-13 pre-registration as 89 executable constants plus a 41-entry
recomputed forbidden-seed proof, with `P13_TAU` deliberately absent behind a
loudly-failing Tier-2 tripwire — so every threshold, ladder, seed and design choice a
reported Phase-13 number will be scored against is committed before any Phase-13 code
runs and before a single labelled datum exists.

## What Was Built

**`spike/p13/consts.jl` (720 lines, Tier 1 only, ONE guard block keyed on `:P13_DEV_SEED`)**

- **Section A — seeds, salts, streams (D-01).** Seven named reserved seeds read from their
  defining sites, `P11_DEV_SEED` forbidden *even though Phase 11 has not run yet*, 25 burned
  spike/probe keys, the 8-salt repository inventory, and both frozen ship-gate seed families
  **recomputed** through an isolated read-only module rather than trusted from a comment.
  `_p13_forbidden()` returns 41 distinct `UInt64`s; `_p13_forbidden_list()` keeps the
  declaration order so a duplicate stays visible (asserted). The fresh `p13_rng` /
  `p13_fix_rng` streams, six reserved counters (fixtures on 99, never a reported counter),
  and `P13_GLOBAL_RNG_DISCIPLINE` recording the `Random.seed!`-before-every-`train` rule that
  Flux's `DataLoader(shuffle = true)` makes mandatory.
- **Section B — F5 mixture.** `P13_IMSIZE_SET` / `P13_IMSIZE_WEIGHTS` **read** from the frozen
  gate pre-registration, never retyped (train-joint == eval-joint invariant).
- **Section C — the D-05 cut.** `:tau_contrast` (variant (a)) with the full class definition and
  the counter-example that makes the decision load-bearing: `(ρ_s, ρ_c) = (0.8, 0.9)` has
  `Δρ < 0` and a naive three-way `Δρ` extension would publish it as "mutually exclusive".
- **Section D — the D-06 probe spec.** δ grid, the 0.90 AUC bar as a *number*, R = 400 per arm,
  the mask-weighted `m̄` statistic, the **unpaired** design (the single most important
  difference from Phase 11's probe), both-directions max, bootstrap B, the reference-λ *rule*
  plus expected value, the grid-never-extended abort criterion, and the
  post-Phase-11-simulator requirement.
- **Section E — the D-07 design.** `:class_frequency` primary with the exactness theorem, the
  pre-declared `:importance_weighted` fallback, `P13_ITERATION_ALLOWANCE = 1` with its ONE
  named trigger, and the F5 verification bars including the mandatory negative control and the
  deliberately skewed toy frequencies.
- **Section F — the D-09/D-10 recipe.** All nine numbers copied verbatim from
  `src/amortized/train_ratio.jl`, the smoke-gate pair, the λ-placement rule with both wrong
  forms named, and `P13_INPUT_WIDTH_RULE` (never write 321 as a literal).
- **Section G/H — the D-12/D-13 gates.** AUC floors, evaluation count, argmax-is-descriptive,
  continuity-reported-not-gated; ECE band, bin count, **ECE-is-the-gate / MCE-is-not** with the
  empty-bin trap in full, the minimum per-head evaluation count, the vacuous-pass AUC floor,
  and `P13_TOST_REQUIRED = false`.
- **Section I — the D-16 alpha ladder.** Nine rungs, the strictly-positive background floor,
  the frozen ch1 Otsu mask computed once, the **two-valued** `P13_ALPHA_SUBSTRATE`, the
  **corrected** `:count_iszero_preserved` invariant, the mask-fraction guard, the
  assert-and-record max-value bound, and the `_exclude_zero` trap.
- **Section I2 — the D-15 real-image arm.** The three not-a-gate sentinels, the pre-registered
  conditions/channel pairs, the grid-truncation math and transposition footnote, the frozen
  `ghat` anchors, the shared alpha grid, the λ read/headline *rules* with the `[ASSUMED]` A11
  domain note, the OOD comparison with the frozen shipped reference numbers, and the two
  verbatim honesty strings (naming correction, target-substitution record).
- **Section J + 30 executable `@assert`s.** Decoupling, deferrals, the 4-tuple of declared
  deviations, and self-checks that run at include time — ending with
  `@assert !isdefined(@__MODULE__, :P13_TAU)`.

**`spike/test/test_p13_consts.jl` (287 lines, 13 nested testsets, 132 assertions)** re-asserts
every locked value **as a literal**, re-runs the seed-disjointness proof against the recomputed
gate families, proves counter separation and stream reproducibility, proves the Tier-2 slot is
empty (`@test_throws ErrorException p13_tau()`), and runs the comment-stripped source greps.

## Verification Performed

All five plan-level verification steps were **run** and their real output observed:

| # | Command | Observed |
|---|---------|----------|
| 1 | `julia --project=spike -e 'include("spike/p13/consts.jl"); println("tier1 ok")'` | printed `tier1 ok`, exit 0 |
| 2 | `julia --project=spike spike/test/test_p13_consts.jl` | `132 Pass / 132 Total`, 0 fail, 0 error, exit 0 |
| 3 | `git diff --quiet HEAD -- src/` | exit 0 (`src/` byte-unchanged) |
| 4 | `git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml` | exit 0 (manifests byte-unchanged) |
| 5 | `julia --project=spike -e '... @assert !isdefined(Main, :P13_TAU) ...'` | printed `tier2 absent`, exit 0 |

Task-level acceptance criteria, all **run**:

- guard blocks outside comments: `1` (expected 1)
- `gate_consts_8_v2` references: `4` (expected ≥ 1)
- `(1376, 1028)` outside comments: `0` — the F5 mixture is read, never retyped
- `all(... .> 0)` outside comments: `0` — the forbidden real-data invariant is never written
- `P13_REAL_*` bar-shaped names outside comments, word-boundary-anchored: exactly
  `P13_REAL_ANCHOR_TOL` and `P13_REAL_OOD_SHIPPED_THRESHOLD` (see Deviation 2)
- locked-value spot-check (`P13_ITERATION_ALLOWANCE`, `P13_CUT_VARIANT`, `P13_STRATIFICATION`,
  `P13_GATE_STATISTIC`): printed `ok`
- tau tripwire: printed `tau tripwire ok`
- substrate/zero-invariant/shared-grid: printed `substrate ok`
- real-arm sentinels: printed `real arm ok`
- `grep -c 'p13_tau' spike/test/test_p13_consts.jl` → `2` (≥ 1)
- `grep -c 'startswith(strip' spike/test/test_p13_consts.jl` → `1` (≥ 1)
- `julia --project=spike -e 'include("spike/test/test_p13_consts.jl")'` → exit 0 (includable
  from a harness, not only runnable standalone)

**Not verified (and deliberately out of this plan's scope):** the full spike suite
(`spike/test/runtests.jl`) and the root `Pkg.test()` were NOT run. `runtests.jl` is owned by
the concurrently-running Phase-11 executor and is not in this plan's `files_modified`, so
`test_p13_consts.jl` is not yet wired into the harness — it is proven includable and
standalone-runnable instead. The wave-level suite run belongs to the orchestrator.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 — Blocking] The isolated `module _GC` cannot live inside the guard block**

- **Found during:** Task 1
- **Issue:** The plan's action places `module _GC ... end` inside
  `if !isdefined(@__MODULE__, :P13_DEV_SEED)`. Julia rejects that outright:
  `ERROR: syntax: "module" expression not at top level` (reproduced before writing the file).
- **Fix:** `module _GC` sits at file top level, immediately above the guard block, with a
  comment recording *why* it is outside it and noting that idempotency comes from the callers'
  guarded-include idiom. This is byte-for-byte the pattern the house twin
  `spike/validation/p11_consts.jl:66-77` already uses for the same file and the same reason —
  so this is a convergence on the established analog, not an invention. The isolation,
  read-only posture, two-level path and the `gate_consts_8_v2` `key_links` pattern are all
  preserved, and `grep -c 'if !isdefined'` outside comments is still `1`.
- **Files modified:** `spike/p13/consts.jl`
- **Commit:** `c42cc8e`

**2. [Rule 1 — Bug] The plan's not-a-gate grep is over-broad against a name the same plan mandates**

- **Found during:** Task 1 acceptance checking
- **Issue:** The plan's literal criterion
  `grep -oE 'P13_REAL_[A-Z_]*(FLOOR|THRESHOLD|MIN|MAX|TOL)'` emits a **third** hit,
  `P13_REAL_NAMIN`, because `MIN` occurs as a substring inside `NAMING` — i.e. inside
  `P13_REAL_NAMING_CORRECTION`, a constant the same plan mandates by name and whose contents
  Task 2 is required to assert. The hit is a regex artifact, not a threshold.
- **Fix:** The constant keeps its mandated name (downstream plan 13-16 and 13-VALIDATION both
  reference the naming correction). The **machine-checked** gate — the source assertion in
  `spike/test/test_p13_consts.jl`, which is what 13-VALIDATION row T-13-54 actually scores —
  anchors the pattern on both ends (`\bP13_REAL_[A-Z_]*(?:FLOOR|THRESHOLD|MIN|MAX|TOL)\b`) and
  asserts **set equality** with the two documented exemptions. Under the anchored pattern the
  comment-stripped source yields exactly `P13_REAL_ANCHOR_TOL` and
  `P13_REAL_OOD_SHIPPED_THRESHOLD` and nothing else — verified.
- **Note for a verifier:** running the plan's *unanchored* one-liner will show three lines. The
  substantive requirement ("no `P13_REAL_*` constant is a pass/fail threshold outside the two
  documented exemptions") holds and is asserted; only the plan's shell mechanism was
  substring-unsafe.
- **Files modified:** `spike/test/test_p13_consts.jl`
- **Commit:** `2e13552`

**3. [Rule 1 — Bug] Case mismatch broke the file's own self-check**

- **Found during:** Task 1 verification (the include failed on its own `@assert`)
- **Issue:** `P13_REAL_NAMING_CORRECTION` was first written with `BIOLOGICAL` in caps for
  emphasis, while the plan-mandated self-check is the case-sensitive
  `occursin("biological", ...)`. The include aborted — the self-check did its job.
- **Fix:** lowercased to `biological test conditions`; the emphasis now rides on the
  surrounding `NOT colocalization labels`.
- **Files modified:** `spike/p13/consts.jl`
- **Commit:** `c42cc8e` (fixed before the commit)

### Conservative additions (Rule 2, strictly strengthening)

- `0x0000_0000_0000_0007` added to `P13_BURNED_DEV_SEEDS`: 13-RESEARCH §K1 lists
  `Random.seed!(7)` among the 11-RESEARCH probe burns but the plan's enumeration named only the
  four hex keys. Forbidding one more observed key can only strengthen the disjointness claim.
- `P13_REPO_SALTS` includes `P11_SALT = 0xA24B_AED4_663E_E121` read from
  `spike/validation/p11_consts.jl:143` (the plan listed the value; the provenance is now in the
  comment).
- The four `0x0000_0000_00B7_A77{0,1,2,3}` keys are the widened form of the repo's actual
  `0x00B7A771` / `0x00B7A772` literals (read from
  `.planning/spikes/006-bf-attrition-mechanism/measure_attrition.jl:25` and
  `.planning/spikes/008-bf-baseline-remedies/logspace_baseline.jl:75`). 13-RESEARCH §K1 renders
  them as a malformed 17-nibble literal; the defining sites were used instead, per the
  read-the-literal-out-of-its-defining-site rule.

## Assumption Drift (advisory)

**1. `spike/validation/p11_consts.jl` already exists and already solved two of this plan's problems**

- **Found during:** Task 1 (looking for a collision-safe way to declare `NPE_MASTER_SEED`)
- **Planned:** 13-PATTERNS.md and the plan present `spike/validation/consts.jl` +
  `test/gate/gate_consts_8_v2.jl` as a *two-part composite* analog, and flag the
  "measured-after-the-probe Tier-2 slot" as a **Phase-13-specific addition with no analog** that
  "the executor must invent".
- **Actual:** Phase 11 has already landed `spike/validation/p11_consts.jl`, which is a
  single-file, near-exact analog: same two-tier structure with two differently-keyed guard
  blocks, the same isolated read-only gate-consts module (with the module-not-at-top-level
  problem already solved and commented), the same recompute-the-derived-seeds discipline, and an
  actual **Tier-2 appended block**. The Tier-2 marker idiom did not need inventing; it needed
  copying.
- **Why it matters:** the resulting file is closer to house style than the plan anticipated, and
  a reader comparing `p13/consts.jl` with `validation/p11_consts.jl` should see deliberate
  parallel structure rather than two independent designs.

**2. The Phase-11 constants the real-image arm defers to are already committed, not pending**

- **Planned:** `P13_TAU_REFERENCE_LAMBDA_EXPECTED` and `P13_REAL_LAMBDA_HEADLINE_EXPECTED` are
  framed as expectations to be checked "once Phase 11 lands".
- **Actual:** `LAMBDA_MIN = 0.25` and `LAMBDA_MAX = 3.0` are already committed literals at
  `spike/validation/p11_consts.jl:166-172`, and `SC2_RUNGS` already ends at 3.0. The expected
  value 3.0 is therefore already consistent with the landed constant.
- **Why it matters:** it does not change what is locked — the *rule* stays Tier-1 and the
  realized number must still be asserted against Phase-11's constant at read time (plans 13-10
  and 13-16) — but the "Phase 11 has not landed" framing in the file comments is now only
  partially true (its consts have landed; its trained net has not, which is what D-02 actually
  blocks on).

## Known Stubs

None. This plan produces a pre-registration and its assertion suite; both are complete and
executable. `P13_TAU` is **not** a stub — its absence is the pre-registration guarantee, it is
asserted absent, and calling `p13_tau()` fails loudly with a pointer to plan 13-10.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no schema and no new file-access
pattern. It reads exactly one file outside `spike/` (`test/gate/gate_consts_8_v2.jl`) and reads
it read-only through an isolated module, which is the disposition the plan's own register
(T-13-03) prescribes. All paths are `joinpath(@__DIR__, ...)` (T-13-06). No package was
installed (T-13-SC): `spike/Project.toml` and `spike/Manifest.toml` are byte-unchanged.

## Locked Tier-1 Constants

`spike/p13/consts.jl` @ commit **`c42cc8e`** (720 lines, one guard block on `:P13_DEV_SEED`).

### A. Seeds, salts and streams (D-01)

| Constant | Committed value |
|---|---|
| `P13_DEV_SEED` | `0x000000000b13de71` |
| `P13_FIX_SEED` | `0x000000000b13f1f7` |
| `P13_SALT` | `0x2545f4914f6cdd1d` |
| `DEFAULT_MASTER_SEED` | `0x0000000000000001` (forbidden) |
| `NPE_MASTER_SEED` | `0x0000000000c0ffee` (forbidden) |
| `VAL_MASTER_SEED` | `0x000000005bc0ffee` (forbidden) |
| `VAL_FIX_SEED` | `0x0000000000f1f7ed` (forbidden) |
| `RATIO_PAIR_SEED` | `0x00000000004a7107` (forbidden) |
| `CORPUS_MASTER_SEED` | `0x0000000000c05eed` (forbidden) |
| `P11_DEV_SEED` | `0x000000000b11de71` (forbidden before Phase 11 runs) |
| `P13_BURNED_DEV_SEEDS` | 25 keys |
| `P13_REPO_SALTS` | 8 salts |
| `_p13_forbidden()` | 41 distinct `UInt64` (incl. 4 recomputed `PROD_SEED` + 4 `PROD_SEED_V2`) |
| `P13_TAU_COUNTER` / `P13_DATAGEN_COUNTER` / `P13_GATE_COUNTER` | `1` / `2` / `3` |
| `P13_ALPHA_COUNTER` / `P13_CONTINUITY_COUNTER` / `P13_FIXTURE_COUNTER` | `4` / `5` / `99` |
| `P13_GLOBAL_RNG_DISCIPLINE` | doc string (two-layer RNG rule) |

### B. F5 image-size mixture (READ, never retyped)

| Constant | Committed value |
|---|---|
| `P13_IMSIZE_SET` | `((512, 512), (1024, 1024), (1376, 1028), (2048, 2048))` |
| `P13_IMSIZE_WEIGHTS` | `(0.4, 0.25, 0.25, 0.1)` |

### C/D. The D-05 cut and the D-06 tau probe

| Constant | Committed value |
|---|---|
| `P13_CUT_VARIANT` | `:tau_contrast` |
| `P13_TAU` | **absent (Tier 2)** — `p13_tau()` errors |
| `P13_TAU_DELTA_GRID` | `(0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2)` |
| `P13_TAU_AUC` | `0.9` |
| `P13_TAU_R` | `400` |
| `P13_TAU_STATISTIC` | `:mbar_mask_weighted` |
| `P13_TAU_DESIGN` | `:unpaired` |
| `P13_TAU_BOTH_DIRECTIONS` | `true` |
| `P13_TAU_BOOTSTRAP_B` | `1000` |
| `P13_TAU_REFERENCE_LAMBDA_RULE` | `:widest_rung` |
| `P13_TAU_REFERENCE_LAMBDA_EXPECTED` | `3.0` |
| `P13_TAU_ABORT_EXTEND_GRID` | `false` |
| `P13_TAU_REQUIRE_POST_P11_SIMULATOR` | `true` |

### E. D-07 stratification and F5 verification bars

| Constant | Committed value |
|---|---|
| `P13_STRATIFICATION` | `:class_frequency` |
| `P13_TARGET_CLASS_FREQ` | `(exclusion = 1//3, random = 1//3, coloc = 1//3)` |
| `P13_STRATIFICATION_FALLBACK` | `:importance_weighted` |
| `P13_ITERATION_ALLOWANCE` | `1` |
| `P13_ITERATION_TRIGGER` | string; one named condition (exclusion AUC below floor at high \|ρ\|) |
| `P13_F5_CORR_MIN` | `0.99` |
| `P13_F5_MAXABS_TOL` | `0.25` |
| `P13_F5_CENTRAL_FRAC` | `0.9` |
| `P13_F5_SKEW_FREQ` | `(0.6, 0.25, 0.15)` |
| `P13_F5_NEGCTRL_REQUIRED` | `true` |

### F. D-09/D-10 architecture and copied recipe

| Constant | Committed value |
|---|---|
| `P13_SUMMARY_WIDTH` | `256` |
| `P13_NUM_SUMMARIES` | `64` |
| `P13_TRAIN_N` | `48000` |
| `P13_EPOCHS` | `300` |
| `P13_BATCHSIZE` | `128` |
| `P13_LR` | `0.00025` |
| `P13_WEIGHT_DECAY` | `0.0001` |
| `P13_VAL_FRAC` | `0.15` |
| `P13_STOPPING_EPOCHS` | `40` |
| `P13_USE_GPU` | `false` |
| `P13_SMOKE_N` / `P13_SMOKE_EPOCHS` | `200` / `2` |
| `P13_LAMBDA_PLACEMENT` | `:append_after_pair_encode` |
| `P13_INPUT_WIDTH_RULE` | `:ratio_input_dim_plus_ncond` |

### G/H. D-12 gate and D-13 calibration

| Constant | Committed value |
|---|---|
| `P13_GATE_M` | `4000` |
| `P13_AUC_FLOOR_COLOC` | `0.9` |
| `P13_AUC_FLOOR_EXCLUSION` | `0.9` |
| `P13_CONFUSION_RULE` | `:argmax_descriptive` |
| `P13_CONTINUITY_GATED` | `false` |
| `P13_ECE_GREEN` | `0.05` |
| `P13_ECE_YELLOW` | `0.1` |
| `P13_ECE_NBINS` | `10` |
| `P13_GATE_STATISTIC` | `:ece` |
| `P13_MIN_EVAL_PER_HEAD` | `1000` |
| `P13_VACUOUS_AUC_FLOOR` | `0.6` |
| `P13_TOST_REQUIRED` | `false` |

### I. D-16 alpha ladder (two substrates)

| Constant | Committed value |
|---|---|
| `P13_ALPHA_LADDER` | `(0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0)` |
| `P13_ALPHA_BG_QUANTILE` | `0.05` |
| `P13_ALPHA_MASK_RULE` | `:calculate_mask_ch1_once` |
| `P13_ALPHA_N_IMAGES` | `64` |
| `P13_ALPHA_INVARIANT_RTOL` | `1.0e-10` |
| `P13_ALPHA_GATED` | `false` |
| `P13_ALPHA_SUBSTRATE` | `(:simulated, :real)` |
| `P13_ALPHA_ZERO_INVARIANT` | `:count_iszero_preserved` |
| `P13_ALPHA_MASK_FRACTION_BOUNDS` | `(0.01, 0.4)` |
| `P13_ALPHA_MAX_VALUE_BOUND` | `1.0` |

### I2. D-15 real-image qualitative arm

| Constant | Committed value |
|---|---|
| `P13_REAL_IS_GATED` | `false` |
| `P13_REAL_QUALITATIVE_ONLY` | `true` |
| `P13_REAL_READ_ONLY` | `true` |
| `P13_REAL_CONDITIONS` | `("positive", "negative")` |
| `P13_REAL_SAMPLE` / `P13_REAL_CONTROL` | `"positive"` / `"negative"` |
| `P13_REAL_CHANNEL_PAIR` | `(1, 2)` |
| `P13_REAL_REDUNDANCY_PAIR` | `(1, 3)` |
| `P13_REAL_IMSIZE` | `(1028, 1376)` |
| `P13_REAL_GRID_TRUNCATION_ROWS` | `4` |
| `P13_REAL_ANCHOR_MBAR` | `(positive = 0.3292, negative = 0.2481)` |
| `P13_REAL_ANCHOR_TOL` | `0.001` (exemption 1: ingestion regression tolerance) |
| `P13_REAL_ALPHA_GRID` | `== P13_ALPHA_LADDER` |
| `P13_REAL_LAMBDA_READS_RULE` | `:full_phase11_ladder` |
| `P13_REAL_LAMBDA_HEADLINE_RULE` | `:widest_rung` |
| `P13_REAL_LAMBDA_HEADLINE_EXPECTED` | `3.0` |
| `P13_REAL_OOD_COMPARISON` | `true` |
| `P13_REAL_OOD_SHIPPED_DENSITY` | `433.69` |
| `P13_REAL_OOD_SHIPPED_THRESHOLD` | `179.14` (exemption 2: shipped-net reference value) |
| `P13_REAL_NAMING_CORRECTION` | string (contains `biological`, `0.2481`) |
| `P13_REAL_SUBSTITUTION_RECORD` | string (contains `Phase 16`) |

### J. Decoupling, deferrals, declared deviations

| Constant | Committed value |
|---|---|
| `P13_SRC_UNTOUCHED` | `true` |
| `P13_PHYSICAL_ANCHOR_DEFERRED_TO` | `"Phase 16"` |
| `P13_RESULTS_RENAME_DEFERRED` | `true` |
| `P13_DECLARED_DEVIATIONS` | 4 strings (custom estimator subtype; two-logit head; masked two-term BCE; per-head evidence-scale correction) |

## Handoff Notes for Later Phase-13 Plans

1. **Plan 13-10 (tau probe)** must APPEND a second guard block keyed on a sentinel other than
   `:P13_DEV_SEED`, and must **relocate** — never delete —
   `@assert !isdefined(@__MODULE__, :P13_TAU)` so it stays above the appended block. It must
   also flip the `"tau is Tier-2 and absent"` testset to assert the measured value, in a
   separately-committed edit, and assert the realized reference λ against Phase-11's
   `LAMBDA_MAX`.
2. **No plan may EDIT a Tier-1 line.** Constants are appended only; the git history of
   `spike/p13/consts.jl` is the audit trail, and `spike/test/test_p13_consts.jl` will fail
   loudly on any modification.
3. **`spike/test/test_p13_consts.jl` is not yet in `spike/test/runtests.jl`.** Wiring it in
   belongs to whichever plan owns the harness edit (the file was concurrently held by the
   Phase-11 executor during this plan). It is proven includable and standalone-runnable.
4. **Two names to avoid re-declaring** if a Phase-13 file is ever loaded into the same module as
   `spike/validation/p11_consts.jl`: both files declare `NPE_MASTER_SEED`, `VAL_MASTER_SEED`,
   `VAL_FIX_SEED`, `DEFAULT_MASTER_SEED`, `CORPUS_MASTER_SEED` and `P11_DEV_SEED` — with
   **identical values**, so a co-load is at worst a redefinition warning, never a value
   conflict. This was checked deliberately.

## Self-Check: PASSED

Created files (checked with `[ -f ... ]`):

- `FOUND: spike/p13/consts.jl`
- `FOUND: spike/test/test_p13_consts.jl`

Commits (checked with `git log --oneline --all | grep -q`):

- `FOUND: c42cc8e` — `feat(13-01): lock the Phase-13 Tier-1 pre-registration`
- `FOUND: 2e13552` — `test(13-01): assert every locked Tier-1 constant as a literal`

No file deletions in either commit (`git diff --diff-filter=D` empty for both). No untracked
files remain beyond this SUMMARY.
