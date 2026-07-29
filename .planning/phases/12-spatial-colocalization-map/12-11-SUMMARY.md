---
phase: 12-spatial-colocalization-map
plan: 11
status: blocked
subsystem: D-12 Stage-1 gate (leave-region-out ridge across the r1 ladder)
tags: [gate, pre-registration, ridge, leave-region-out, positive-control, information-limit, adjudication-blocked]
requires: ["12-01", "12-03", "12-04", "12-07", "12-09", "12-10"]
provides:
  - "spike/validation/run_p12_stage1_ridge.jl — the D-12 Stage-1 gate: 64 leave-region-out ridges per rung across the five pinned-r1 rungs, with an included-row positive control, an empirical per-region prior baseline, a reported rho-space arm, a reported radial-R2, and a reported own-row information limit"
  - "spike/validation/p12_stage1_report.jld2 — per-rung ratios (mean/min/max/per-region), control ratios, rho arm, radial R2, within-image contrasts, selected penalties, liveness diagnostics, and the two gate components"
  - ".planning/phases/12-spatial-colocalization-map/12-STAGE1-ADJUDICATION-BLOCKED.md — the measurement, the rule, the evidence, and the one question a human must rule on"
affects: ["12-14", "12-15", "12-16", "12-17", "12-18", "12-19", "12-20"]
tech-stack:
  added: []
  patterns:
    - "a positive control is only a liveness proof if its bar is reachable; audit it against an information limit computed WITHOUT the estimator it judges"
    - "a two-state banner is a defect when the governing rule has three branches — a failed control may not print DESCOPE"
    - "an unpenalised solve is unavailable by construction on a masked design: zeroing a region's rows makes the Gram exactly rank-deficient"
key-files:
  created:
    - spike/validation/run_p12_stage1_ridge.jl
    - spike/validation/p12_stage1_report.jld2
    - .planning/phases/12-spatial-colocalization-map/12-STAGE1-ADJUDICATION-BLOCKED.md
    - .planning/phases/12-spatial-colocalization-map/12-11-SUMMARY.md
  modified: []
decisions:
  - "NO VERDICT RECORDED, AND 12-STAGE1-VERDICT.md DELIBERATELY NOT CREATED. borrowing_ok = TRUE at 5 of 5 rungs; control_live = FALSE. The frozen interpretation rule forbids BOTH verdicts when the control fails, and its third branch ('harness not live') is contradicted by the artifact. `p12_stage1_verdict()` therefore returns `:absent` — whose own docstring defines it as 'the Stage-1 gate has not been adjudicated yet', which is exactly true — and `p12_require_proceed` refuses 12-17/18/19/20."
  - "THE CONTROL IS AT THE INFORMATION LIMIT OF THE DATA, not broken. The own-row ridge attains sqrt(1 - corr^2) to within 5e-4 at every rung, and that limit is 0.766/0.712/0.637 at r1 = 0.05/0.25/0.50 — above the 0.5 ceiling, so no estimator of any kind can clear it there. The Phase-11 benchmark also reproduces: the global-level control measures 0.19687 at r1 = 0.95 against Phase 11's 0.157."
  - "P12_STAGE1_CONTROL_CEILING CARRIES NO DERIVATION IN THE PRE-REGISTRATION, unlike P12_STAGE1_RATIO_CEILING which carries one verbatim plus its own mis-scaled-unit warning. Whether 0.5 is a live bar or a mis-scaled one is a PRE-REGISTRATION question and was not answered here."
  - "THE FINDING STANDS REGARDLESS: neighbours do inform a held-out region, monotonically in r1 (ratio_mean 0.926 -> 0.336), and the crossing at which neighbours beat the region's OWN measurement lies between r1 = 0.75 (parity, delta 0.0063) and r1 = 0.95."
  - "NO ITERATION SPENT. P12_ITERATION_ALLOWANCE stands at 1. No Tier-1 constant edited, no Tier-2 constant appended, no threshold moved."
metrics:
  duration: ~55 min
  completed: 2026-07-29
  tasks: 2
  commits: 3
  files: 3
---

# Phase 12 Plan 11: The D-12 Stage-1 Gate — Summary

**Status: BLOCKED on a pre-registration ruling.** The measurement is complete, reproducible and
committed. The gate's two halves disagreed: the half the phase is about **passed at every rung**,
and the half that exists to prove the instrument works **did not clear a ceiling that the data
makes unreachable at the three shortest correlation lengths**. The frozen interpretation rule
forbids recording either verdict in that case, and its third branch asserts something the artifact
contradicts. Resolving that is a user decision, so no verdict was written.

Full detail, with every number re-derived from the committed artifact:
`.planning/phases/12-spatial-colocalization-map/12-STAGE1-ADJUDICATION-BLOCKED.md`.

---

## 1. What was measured

Five pools of 5,000 samples, `arm = :car`, image size pinned to `(512, 512)`, `r1` **PINNED** per
rung and reaching the content hash (five distinct hash directories, asserted before generation;
every realized `r1` and every realized image size asserted after loading). Realized `n_test = 1250`
per rung, equal to `P12_STAGE1_N_TEST_PER_RUNG_DERIVED`, asserted at every rung.

| r₁ | `ratio_mean` | `ratio_min` | `ratio_max` | `control_ratio_mean` | `rho_arm_ratio_mean` | `radial_r2` | vs 0.95 |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.92557 | 0.84632 | 1.01020 | 0.73188 | 0.92439 | 0.34522 | **PASS** |
| 0.25 | 0.84660 | 0.75726 | 0.97470 | 0.67740 | 0.84641 | 0.37843 | **PASS** |
| 0.50 | 0.71870 | 0.61900 | 0.85285 | 0.59827 | 0.71985 | 0.25225 | **PASS** |
| 0.75 | 0.56184 | 0.48333 | 0.70133 | 0.49659 | 0.56276 | 0.41096 | **PASS** |
| 0.95 | 0.33569 | 0.29381 | 0.42932 | 0.32369 | 0.32311 | 0.74890 | **PASS** |

```
borrowing_ok = 5 of 5 rungs at or below 0.95   (1 required)   -> TRUE
control_live = max(control) 0.73188 <= 0.5                    -> FALSE
stage1_pass  = borrowing_ok && control_live                   -> FALSE
```

The gated target is the Gaussian `z_field[r]` (R-2, atom-free); the ρ arm is reported only and
tracks it to within 0.002 at four of five rungs, so nothing depends on which is quoted.

## 2. Why no verdict

The rule, fixed in the plan before the numbers existed: *"`control_live` false ⇒ the harness is not
live and the whole run is uninformative. Do not record a PROCEED or a DESCOPE. Record 'harness not
live', name the most likely causes in order (a broken standardizer, a target/predictor
misalignment, a pool with the wrong `r1`), and spend the single `P12_ITERATION_ALLOWANCE` on fixing
the harness."*

It forbids both verdicts. Its third branch asserts the harness is not live — and the artifact says
otherwise. All three named causes are ruled out executably (masked rows verified exactly zero raw
and standardized; column-major `p12_idx` alignment confirmed by an own-row correlation rising 0.635
→ 0.881 with r₁; every realized `r1` asserted equal to its pin in five distinct hash directories).
And the decisive evidence, which the rule does not ask for:

| r₁ | own-row ridge | its information limit `sqrt(1-corr²)` | Δ | full-128 control | global-level control |
|---|---|---|---|---|---|
| 0.05 | 0.76588 | 0.76564 | 0.00024 | 0.73188 | 0.77173 |
| 0.25 | 0.71236 | 0.71201 | 0.00035 | 0.67740 | 0.61400 |
| 0.50 | 0.63716 | 0.63667 | 0.00049 | 0.59827 | 0.38141 |
| 0.75 | 0.55590 | 0.55551 | 0.00039 | 0.49659 | 0.23234 |
| 0.95 | 0.47199 | 0.47183 | 0.00016 | 0.32369 | 0.19687 |

The ridge **attains** the analytic limit at every rung, and the limit itself is **above the 0.5
ceiling** at r₁ ≤ 0.50 — so at those rungs no estimator, correct or broken, can reach the bar. The
Phase-11 benchmark reproduces too: 0.19687 against Phase 11's 0.157 for global ρ_true, at the rung
where the comparison is meaningful and at the noisiest single size in the F5 set.

So all three writable contents the plan specifies — `PROCEED`, `DESCOPE`, `harness not live` — are
untrue. Per the dispatch's standing instruction, no verdict was picked and `status: blocked` is
recorded instead.

**The open question, in one sentence:** does a positive control that is provably at the information
limit of the data, but above `P12_STAGE1_CONTROL_CEILING = 0.5` at r₁ ≤ 0.50, mean the harness is
not live — or that the ceiling is unreachable in that unit at those rungs? Both readings are
defensible from the record (§6 of the blocked document lays out each and its consequence).
`P12_STAGE1_CONTROL_CEILING` carries **no derivation** in the pre-registration, unlike the ratio
ceiling which carries one verbatim together with its own mis-scaled-unit warning; and 12-RESEARCH's
qualitative form of the same bar — *"must be far below 1.0"* — is met at every rung (a 27–68 %
error reduction). That is exactly why it is not the executing agent's call.

## 3. The finding, which stands whatever the ruling

Neighbouring regions **do** carry information about a held-out region, monotonically in the
correlation length: `ratio_mean` 0.926 → 0.847 → 0.719 → 0.562 → 0.336, clearing the pre-registered
0.95 ceiling at all five rungs. The crossing at which the other 63 regions become **more**
informative than the region's own measurement lies between r₁ = 0.75 (parity, Δ = 0.0063) and
r₁ = 0.95 (0.336 vs 0.472) — the deliverable Pitfall 2 named in advance, and it has the shape
Pitfall 2 predicted.

Quoted against the measured own-row limit rather than against the research's 0.0662 per-draw noise
sd: that figure was measured at 1376×1028 while this ladder is pinned to 512², and Pitfall 2 itself
records the noise sd as image-size dependent, so the two are not comparable and were not compared.

**S-4 early warning, as the plan requires:** `radial_r2` = 0.345 / 0.378 / 0.252 / 0.411 /
**0.749**. Three quarters of the per-region recoverability spread is radial at r₁ = 0.95. Carry
into 12-20's guard interpretation **whatever the ruling** — together with the competing benign
reading, that at long correlation lengths the spread is small and dominated by lattice geometry
(12-07 measured CAR marginal sd largest at low-degree corner/edge cells, itself a radial pattern
with no chromatic content). This run does not separate the two; 12-20's ε = 0 and offset-grid arms
do.

## 4. Verification actually run

| Check | Command | Observed |
|---|---|---|
| Task-1 verify, verbatim | `julia --project=spike -e 'P12_STAGE1_LOAD_ONLY = true; include(...)'` | `stage1-load-ok` |
| Task-2 verify, verbatim | the plan's `JLD2.load` + mandatory-key assertion | `false true false [0.9256, …] [0.7319, …]` — all eight keys present, `ratio_mean` a vector |
| header strings | source read | contains `GATES`; does **not** contain `gates nothing` |
| single masking implementation | comment-stripped source | `mask_regions(` appears **once**; no second implementation |
| no threshold literal of its own | source read | references all three Tier-1 gate constants; the only numerics are the ridge penalty grid, the split fractions and the analog's `1e-12` zero-variance guard |
| no abort on the gate outcome | source read | 13 `@assert`, 0 `throw`; every assertion structural (pool completeness, dims, θ arity, r1 pin, imsize, n_test, distinct dirs, artifact integrity) |
| routing while unadjudicated | `p12_stage1_verdict(; dir)` | `:absent`; `p12_require_proceed("12-17")` **threw**, as required |
| blocked doc carries no verdict | `eachmatch(P12_STAGE1_VERDICT_PATTERN, src)` | **0** matching lines |
| runner not in the suite | `grep -rn run_p12_stage1 spike/test/` | no references — the suite is untouched by this plan |
| frozen surfaces | `git diff --quiet HEAD -- spike/validation/p12_consts.jl src corpus spike/Project.toml spike/Manifest.toml` | **CLEAN** |
| `.planning` files I must not touch | `git status --porcelain -- .planning/STATE.md .planning/ROADMAP.md` | empty — neither touched, neither committed |
| commit deletions | `git diff --diff-filter=D --name-only HEAD~1 HEAD` | empty |

### The Phase-11 pool, before and after

Byte-identical — sizes **and** mtimes, all six files (`diff` of the two listings is empty):

```
meta.jld2         11 764 bytes   mtime 1785165536
shard_0001.jld2   11 206 103     mtime 1785163029
shard_0002.jld2   11 206 103     mtime 1785163648
shard_0003.jld2   11 206 103     mtime 1785164275
shard_0004.jld2   11 206 103     mtime 1785164906
shard_0005.jld2   11 206 103     mtime 1785165536
du -sh spike/data/cache/p11 -> 54M      git status --porcelain -- spike/data/cache/p11 -> empty
```

### Compute, measured

Datagen 2.01 / 1.90 / 1.94 / 2.01 / 2.08 min per rung = **9.94 min** for all five (well under the
briefing's ~7.3 min/rung projection, because the rungs are pinned to 512² rather than the F5
mixture the projection was measured on). The 64-region ridge stage costs **~7 s per rung**. First
run end-to-end: 10.06 min.

**The five pools are gitignored and unbacked.** They regenerate byte-identically from the
counter-based Philox stream keyed per global index — but regenerating costs ~10 minutes of CPU.
Nothing here is available for reuse at no further compute cost.

## 5. Deviations from plan

### `[Rule 1 — Bug] An unpenalised solve that is singular by construction on the arm the gate is about`

- **Found during:** Task 1, in a deliberate 64-sample smoke run made **before** spending the
  datagen — which is what the smoke was for.
- **Issue:** my `_p12s1_select_and_fit` seeded its search with `best_w = A \ rhs`, an OLS fit at
  zero penalty. On the **masked** design two rows are exactly zero by construction, so the Gram is
  exactly rank-deficient (126 of 128) no matter how many fit samples there are. `SingularException`,
  on the gated arm, at every region.
- **Fix:** seed at `RIDGE_GRID[1]` instead, so every solve carries a penalty; the reason is recorded
  in the function's docstring so it cannot be reintroduced.
- **Verified:** smoke green; the Gram rearrangement then reproduces the analog's literal `ridge_fit`
  with `max|Δw| = 0.0` (bit-identical).
- **Commit:** `e88b97f`

### `[Rule 2 — Missing critical functionality] A liveness diagnostic, because the rule's premise was untested`

- **Found during:** Task 2, after the first run returned `control_live = false`.
- **Issue:** the frozen rule directs the phase's single iteration allowance at "fixing the harness"
  on the strength of an untested premise. Nothing in the artifact could distinguish a broken harness
  from an unreachable bar, so the run could not be honestly reported either way.
- **Fix:** three REPORTED keys — `own_row_bound` (`sqrt(1 - corr²)`, an information limit computed
  with no ridge, no split and no standardizer, so it cannot share the failure modes of what it
  audits), `ownrow_ratio`, and `global_control_ratio` (the Phase-11-comparable benchmark).
- **Discipline:** they enter no gate component; **every pre-existing artifact key is bit-identical**
  across the two runs (verified key-by-key: the only differences are `generated`, `elapsed_min`, and
  the six new keys). Added after the result was seen, which is recorded at the site, in the artifact
  and here. No threshold moved.
- **Commit:** `986d518`

### `[Rule 1 — Bug] A two-state banner in a three-branch rule`

- **Issue:** the plan specifies "print a PASS/DESCOPE banner" from `stage1_pass`. But `stage1_pass`
  is an AND, and the plan's own interpretation rule forbids reading a failed **control** as a
  descope. The first run therefore printed `DESCOPE` for a run from which a descope may not be read
  — and the banner is what a human sees.
- **Fix:** three states; `!control_live` prints `NO VERDICT — the positive control did not clear its
  ceiling`, followed by the liveness numbers and an explicit note that which reading applies is a
  pre-registration question. Artifact values unchanged.
- **Commit:** `986d518`

### `[Deliberate] 12-STAGE1-VERDICT.md not created; a differently-named document instead`

The plan's artifact spec requires `12-STAGE1-VERDICT.md` to `contain: "PROCEED"`. It cannot be
written truthfully, and writing it with no matching `VERDICT:` line would make `p12_stage1_verdict`
throw a *malformed-file* error at four later plans rather than report the true state. Leaving it
**absent** makes the same function return `:absent`, documented as *"the Stage-1 gate has not been
adjudicated yet"* — precisely correct. The measurement, the rule, the evidence and the open question
are recorded in `12-STAGE1-ADJUDICATION-BLOCKED.md` in the `11-PROBE-VERDICT.md` shape, carrying
**zero** lines matching the verdict pattern.

### `[Deliberate] Task 2(d)'s STATE.md note not written`

Conditional on a DESCOPE, which was not recorded. The dispatch also forbids touching `STATE.md`.
Neither `STATE.md` nor `ROADMAP.md` was modified or committed.

### Falsification proofs (mutate → red → revert byte-exactly)

Per §6 of the briefing. No committed file was mutated: each mechanism was falsified by driving it
with inputs that break its precondition, so nothing needed reverting.

| Probe | Assertion under test | Result |
|---|---|---|
| Resolve all five rungs at ONE pin inside one cache root (the directory collapse a dropped `r1_pin` produces) | the five-distinct-directories assertion | 1 of 5 distinct → **would fire**; with the pin, 5 of 5 |
| Mask region r, then standardize | the exclusion is real | masked rows exactly `0.0` raw **and** exactly `0.0` standardized, while the unmasked row is `0.804` |
| Purely radial input vector | `_p12s1_radial_r2` | `1.0`; a non-radial index vector gives `0.0` |
| Gram-cached solve vs the analog's literal `ridge_fit` | the rearrangement is not a different estimator | `max|Δw| = 0.0` |

**One weakness reported rather than smoothed over:** the `n_test` assertion compares a **rounded**
count, so it cannot distinguish n = 4999 / 5000 / 5001 (all round to 1250). It catches any n
differing by ≳ 10, which is the failure mode it exists for (a wrong-sized or foreign pool), but it
is not a byte-level guard and should not be quoted as one.

## 6. Assumption Drift (advisory)

- **The ladder image size.** *Planned / pre-registered:* `P12_STAGE1_IMSIZE = (512, 512)`, pinned,
  as Tier 1 fixes and as the plan's `key_links` require. *12-RESEARCH Pitfall 2 consequence (b)*
  says the opposite in its own words: *"the ladder must be run on the F5 mixture the net trains on,
  not on a single convenient size"*, because the per-region noise sd is image-size dependent.
  *Actual:* the pinned size, per Tier 1. *Why:* Tier 1 is frozen and append-only, and its own §6
  compute note shows the single size was chosen knowingly ("the research sized the Stage-1 ridge at
  ~6 min for 5 000 samples at 512-squared"). **Consequence a reader must carry:** 512² is the
  noisiest single size in the F5 set, so every ratio here is a **conservative** (pessimistic)
  estimate of what the same probe would show on the training joint — which cuts in favour of the
  borrowing result and against the control ceiling being reachable. Not resolved here; recorded.
- **Rung independence.** *Implicitly assumed by "five rungs":* five independent experiments.
  *Actual:* `p12_datagen_rng` is keyed by global index alone and pinning `r1` skips the r1 draw, so
  sample i of every rung shares its stream from the field draw onward — the ladder is
  common-random-number **paired**, not independent. *Why it matters:* the comparison is
  variance-reduced and no rung-to-rung difference can be a seed artifact, but the five rungs must
  not be counted as five independent replications. Recorded in the runner header.

## 7. Threat Flags

None. No network endpoint, no auth path, no file access outside the gitignored Phase-12 cache root
and the two committed artifacts, no schema change at a trust boundary. Dispositions applied:
T-12-39 (`mask_regions` zeroes both rows; the included-row design is reported as the control only),
T-12-21 (`control_live` is a mandatory component and its failure blocked the verdict — the threat
this register names is exactly what happened, and the mitigation held),
T-12-40 (no throw on the gate outcome, so the artifact was written and committed),
T-12-41 (`n_test` asserted against the derived count at every rung),
T-12-22 (standardizer fit on FIT only; penalty selected on VALIDATION; the test block touched once),
T-12-06 (`radial_r2` reported per rung and carried forward),
T-12-11 (own cache root; before/after Phase-11 shard listing identical),
T-12-42 (`STATE.md` not touched at all),
T-12-07 (three pathspec-scoped commits on a branch two other agents are committing to).

## 8. Known Stubs

None.

---

## Self-Check: PASSED

- `spike/validation/run_p12_stage1_ridge.jl` — FOUND
- `spike/validation/p12_stage1_report.jld2` — FOUND (tracked, committed)
- `.planning/phases/12-spatial-colocalization-map/12-STAGE1-ADJUDICATION-BLOCKED.md` — FOUND
- `12-STAGE1-VERDICT.md` — deliberately ABSENT; `p12_stage1_verdict()` returns `:absent`
- `e88b97f` `feat(12-11): the D-12 Stage-1 gate …` — FOUND
- `986d518` `feat(12-11): report the Stage-1 liveness limit …` — FOUND
- `2b96ecb` `docs(12-11): Stage-1 measured; the pre-registered rule yields neither verdict` — FOUND
- `spike/data/cache/p11` — 54 MB, six files, sizes and mtimes unchanged
- `src`, `corpus`, `spike/Project.toml`, `spike/Manifest.toml`, `spike/validation/p12_consts.jl` — byte-unchanged
- `.planning/STATE.md`, `.planning/ROADMAP.md` — untouched
