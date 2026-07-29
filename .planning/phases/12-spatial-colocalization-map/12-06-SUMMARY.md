---
phase: 12-spatial-colocalization-map
plan: 06
status: complete
subsystem: spike-validation
tags: [identifiability, chromatic-eps, ridge, reported-only, read-only-pool, S-4-confound]
requires:
  - "spike/data/cache/p11/bff8550…350a (the EXISTING Phase-11 50k pool, read-only)"
  - "spike/data/p11_generate.jl (p11_pool_dir, load_p11_pool)"
  - "spike/validation/p12_consts.jl (P12_EPS_RIDGE_COUNTER, P12_VACUOUS_SHRINKAGE_FLOOR)"
  - "spike/simulator/prior.jl (CHROMATIC_PRIOR, transitively via p11_generate.jl)"
provides:
  - "spike/validation/run_p12_eps_ridge.jl — read-only identifiability ridge on theta row 8"
  - "spike/validation/p12_eps_ridge_report.jld2 — the reported measurement (TRACKED)"
affects:
  - "12-20: S-4's radial chromatic confound is measured LIVE, so the guards carry the SC3 claim"
  - "the phase report's named limit #4 (vacuous nuisance columns): NOT revised — see below"
tech-stack:
  added: []
  patterns:
    - "read-only pool access: p11_pool_dir + load_p11_pool only; shard-presence assertion because p11_pool_dir mkpaths"
    - "empirically-computed prior baseline (fit-split mean/sd), never the analytic value"
    - "validation-selected ridge penalty, deterministic tail-block test split touched once"
    - "ridge_residual_shrinkage named apart from 12-18's shrinkage (F3, two estimators)"
key-files:
  created:
    - "spike/validation/run_p12_eps_ridge.jl"
    - "spike/validation/p12_eps_ridge_report.jld2"
  modified: []
decisions:
  - "Include order INVERTED from the plan's: p11 chain FIRST, p12_consts.jl SECOND (measured sentinel collision)"
  - "A shard-presence assertion was added because the plan's isdir check provably cannot fire"
  - "eps_ratio 0.9631 is BELOW 1.0 but the pre-registered F3 flag says vacuous=true; both reported, neither reconciled"
  - "The two predictor arms are provably the same regression on this pool: the 64 mask rows are identically 1.0"
metrics:
  tasks: 2
  commits: 2
  tests_added: 0
  duration_min: 34
  completed: 2026-07-29
requirements: [SPAT-08]
---

# Phase 12 Plan 06: Chromatic-ε Identifiability — Summary

Measured, on the existing Phase-11 50 000-sample pool and without touching a byte of it, whether
`chromatic_eps` (θ row 8) is linearly recoverable from the 128-row patch-correlation summary. It
is — barely. A 128-predictor ridge beats the empirical prior-mean baseline by **3.7 %**
(`eps_ratio = 0.96307`), with a live positive control (`control_rho_ratio = 0.15742`, effectively
identical to Phase 11's 0.157). Against the pre-registered F3 bar, **96.9 % of the prior spread
survives**, so `ridge_residual_shrinkage = 0.96870 > P12_VACUOUS_SHRINKAGE_FLOOR = 0.90` and the
flag reads **`vacuous = true`**. S-4's radial chromatic confound is therefore **LIVE**.

**This number is REPORTED and gated nothing.** No threshold was changed, no plan was gated, no
descope was authorised, and **nothing was appended to `spike/validation/p12_consts.jl`** — that
file is byte-unchanged at HEAD (asserted below). The phase's single iteration allowance was not
spent.

---

## The measurement, re-derived by loading the artifact

Every number below was read back out of `spike/validation/p12_eps_ridge_report.jld2` with
`JLD2.load`, not copied from the console.

| Key | Value |
|---|---|
| `eps_ratio["all128"]` | **0.9630683102676035** |
| `eps_ratio["continuous64"]` | **0.9630683102676035** |
| `control_rho_ratio` (all128) | **0.15742124368198485** |
| `control_rho_ratio_by_arm["continuous64"]` | 0.15742124368198485 |
| `ridge_residual_shrinkage` (all128) | **0.9687016173325083** |
| `ridge_residual_shrinkage_by_arm["continuous64"]` | 0.9687016173325083 |
| `vacuous` | **true** (both arms) |
| `vacuous_shrinkage_floor` | 0.9 |
| `n_total` / `n_fit` / `n_val` / `n_test` | 50000 / 30000 / 7500 / **12500** |
| selected penalty, ε (both arms) | **1.0e-6** |
| selected penalty, ρ_true control (both arms) | **10.0** |
| `ridge_grid` | [1.0e-6, 1.0e-4, 1.0e-2, 1.0e-1, 1.0, 10.0, 100.0, 1000.0] |
| `elapsed_min` | **0.09196666479110718** (5.5 s) |
| `zscore_arm` | `per_row` |
| `counter` | 4 (`P12_EPS_RIDGE_COUNTER`) |
| `schema_version` / `p11_pool_readonly` | 1 / `true` |
| `generated` | 2026-07-29T16:18:57.439Z |

Supporting quantities:

| Key | Value |
|---|---|
| `eps_prior_mean_fit` | -4.7904723649332935e-5 |
| `eps_prior_sd_fit` (the shrinkage denominator) | 0.011520577079042955 |
| `eps_prior_sd_analytic` (`std(CHROMATIC_PRIOR)`) | 0.011547005383792516 |
| `arms["all128"]["eps_ridge_rmse"]` | 0.011160001649072734 |
| `arms["all128"]["eps_prior_rmse"]` | 0.01158796476853418 |
| `arms["all128"]["control_rho_ridge_rmse"]` | 0.07411825414389875 |
| `arms["all128"]["control_rho_prior_rmse"]` | 0.47082752245071213 |

### Per-|ρ| breakdown (5 bins, `all128`; `continuous64` is bit-identical)

| \|ρ\| bin | n | ridge RMSE | prior RMSE | ratio |
|---|---|---|---|---|
| [0.0, 0.2) | 4346 | 0.010989038374988172 | 0.011648636081467732 | **0.9433755418345552** |
| [0.2, 0.4) | 3229 | 0.01102610018820681 | 0.011526865724464787 | **0.956556660914758** |
| [0.4, 0.6) | 1966 | 0.01120636880556201 | 0.011498656715017739 | **0.9745806908841806** |
| [0.6, 0.8) | 1322 | 0.011200777863048409 | 0.01155860819147627 | **0.9690420920494789** |
| [0.8, 1.0) | 1637 | 0.011767613131494587 | 0.011676850609968417 | **1.0077728596997453** |

---

## Verdict, under the interpretation rules frozen before the number was seen

**Rule 1 — is the harness live?** `control_rho_ratio = 0.15742` is comfortably below 1.0 and
matches Phase 11's measured 0.157 to three decimals (11-DIAGNOSIS.md §2, an 84 % error reduction).
The harness is **LIVE**. The ε result is therefore informative, and a verdict may be given.

**Rules 2 and 3 do not resolve to a single branch, and that is the honest finding.** The plan's
two ε branches were `eps_ratio < 1.0` ⇒ first identified nuisance, and `eps_ratio ≈ 1.0` ⇒
vacuous. The measurement landed **between them, and it landed there on both sides of the same
data**:

- `eps_ratio = 0.96307` is strictly below 1.0, and it is below 1.0 by more than noise. The
  sampling sd of an RMSE ratio at `n_test = 12500` is ≈ `1/sqrt(2·12500) = 0.0063`, so the
  0.0369 shortfall is ≈ 5.8 sd. There is a **real, reproducible linear signal about ε in the
  128-row summary**. That is a genuine first: no nuisance in this project has previously produced
  a ratio distinguishable from 1.0.
- The pre-registered F3 bar, declared in Tier 1 before anything ran, is not on the ratio but on
  **what fraction of the prior spread survives**: `ridge_residual_shrinkage = 0.96870`, against a
  floor of 0.90. **96.9 % of the ε prior spread survives**, so `vacuous = true`.

The two are not in conflict — they answer different questions. "Is there signal?" is yes. "Is
there enough signal to matter?" is, by the bar this phase committed to in advance, **no**. Under
the standard that governs the phase, **ε is VACUOUS**.

Consequently:

- **S-4's radial chromatic confound is LIVE.** The net cannot be expected to explain the radial
  gradient away from a channel that removes 3 % of ε's spread, so **the 12-20 guards carry the
  weight of the SC3 claim**, exactly as the plan's `eps_ratio ≈ 1.0` branch specifies.
- **Named limit #4 (vacuous nuisance columns) does NOT need revision.** ε does not become "the
  project's first identified nuisance" on this evidence; it becomes the first nuisance whose
  vacuity is measured rather than assumed, with a named non-zero-but-immaterial signal. That is a
  strengthening of limit #4's wording, not a retraction of it, and the phase report should say so
  in those terms.
- Nothing here licenses a re-run, a wider grid, or a lower floor. The floor was set in advance and
  is a member of `P12_REPORTING_ONLY_CONSTANTS`; it decides nothing and it was not touched.

### Two structural findings a reader must not misread

**1. The two predictor arms are bit-identical because the mask block is degenerate, not because
the mask rows are irrelevant to a general summary.** `eps_ratio`, `ridge_residual_shrinkage` and
`control_rho_ratio` agree to all 16 digits across `all128` and `continuous64`. Measured cause:
over all 50 000 pool columns the present-mask rows 65:128 are **identically 1.0** (min = max = 1.0,
per-row sd = 0.0 on all 64 rows; all 64 value rows have sd > 1e-12). The train-only standardizer
forces sd = 1 on zero-variance rows, so those 64 standardized rows are **exactly 0** and
contribute exactly nothing to the ridge. The arms are therefore *provably the same regression* on
this pool. The plan's branch "if only (a) works, the mask rows are carrying the signal" is
**unreachable here** — no patch is ever missing at the Phase-11 image sizes — and the second arm
is confirmatory rather than informative. It is still reported, because the degeneracy itself is
the fact.

**2. The per-|ρ| trend runs OPPOSITE to the physics the plan predicted.** The plan's stated reason
for binning was that "the ε perturbation grows with |ρ| by physics (magnification mismatch
degrades an existing correlation; at ρ ≈ 0 there is nothing to degrade)", which predicts the
*lowest* ratio at high |ρ|. The measurement is monotone the other way: 0.9434 at |ρ| ∈ [0, 0.2)
rising to **1.0078** at |ρ| ∈ [0.8, 1.0) — i.e. at high |ρ| the ridge is *worse than the prior
mean*, which is what an overfit direction with no real signal looks like out of sample. Recorded
as an anomaly, not explained: this plan measures, it does not model. Two candidate readings worth
one sentence each for the phase report — (i) the largest bin (n = 4346) has the most statistical
power and the high-|ρ| bins the least (n = 1322–1637, sampling sd ≈ 0.019, so 1.0078 is within
noise of 1.0), and (ii) at high |ρ| the summary rows are dominated by the correlation itself,
leaving less residual variance for a weak radial term. **No follow-up run was made to distinguish
them.**

---

## Phase-11 pool: provably unmutated

Resolved directory (from the artifact's `pool_dir`):
`spike/data/cache/p11/bff8550959d517d796692bae26ed61765bd923c06c010b82e7b187942368350a`

`stat -c '%n %s %y'` on all six files, **BEFORE** the run:

```
meta.jld2        11764     2026-07-27 17:18:56.656029000 +0200
shard_0001.jld2  11206103  2026-07-27 16:37:09.127512100 +0200
shard_0002.jld2  11206103  2026-07-27 16:47:28.396512800 +0200
shard_0003.jld2  11206103  2026-07-27 16:57:55.666234000 +0200
shard_0004.jld2  11206103  2026-07-27 17:08:26.636385300 +0200
shard_0005.jld2  11206103  2026-07-27 17:18:56.096202400 +0200
```

**AFTER** the run:

```
meta.jld2        11764     2026-07-27 17:18:56.656029000 +0200
shard_0001.jld2  11206103  2026-07-27 16:37:09.127512100 +0200
shard_0002.jld2  11206103  2026-07-27 16:47:28.396512800 +0200
shard_0003.jld2  11206103  2026-07-27 16:57:55.666234000 +0200
shard_0004.jld2  11206103  2026-07-27 17:08:26.636385300 +0200
shard_0005.jld2  11206103  2026-07-27 17:18:56.096202400 +0200
```

All six sizes and all six mtimes are **identical to the nanosecond**; the five shards are
11 206 103 bytes each as pre-registered. `du -sb spike/data/cache/p11` reads **56 042 279** bytes
before and after. `ls -1 spike/data/cache/p11/` lists exactly one hash directory before and after,
so no stray empty pool directory was created either. `git status --porcelain -- spike/data/cache`
is **empty**.

This plan generated **no pool and no cache**. It only read one.

## Decoupling assertions (all run, all passed)

| Check | Result |
|---|---|
| `git diff --quiet HEAD -- src spike/Project.toml spike/Manifest.toml corpus spike/validation/p12_consts.jl` | **byte-unchanged** |
| `git diff --quiet HEAD -- spike/simulator spike/npe` | **byte-unchanged** |
| `git diff --quiet HEAD -- .planning/STATE.md .planning/ROADMAP.md` | **untouched** |
| `git status --porcelain -- spike/data/cache` | **empty** |
| stripped-source forbidden tokens (`generate_p11_pool`, `write_p11_shard`, `Pkg.add`, `corpus`) | **zero occurrences** |
| `using` set | `JLD2, Dates, Statistics, LinearAlgebra` only (stdlib via `@stdlib` LOAD_PATH) |
| test suite | **unchanged** — no testset added; grep confirms only `run_p12_eps_ridge.jl` itself and a `p12_consts.jl` counter comment mention the runner, so it is never invoked from `runtests.jl` |

## Verification commands actually run

Task 1 verify (the plan's `<automated>` block, verbatim) → printed `eps-ridge-load-ok`.
Additionally run, to close the gap the acceptance criteria describe but the automated block omits:
the same stripping with `corpus` added to the forbidden list, plus a source read confirming
`best_v` is compared on `Xval` → printed `stripped-source-clean`.

Task 2 verify (the plan's `<automated>` block, verbatim) → exited 0 and printed
`Dict("continuous64" => 0.9630683102676035, "all128" => 0.9630683102676035) 0.15742124368198485
0.9687016173325083 true`.

---

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] The plan's include order does not load — sentinel collision between
`p12_consts.jl` and `p11_consts.jl`**

- **Found during:** Task 1
- **Issue:** The plan specifies "Guarded includes: `p12_consts.jl`, then `../data/p11_generate.jl`."
  That order throws. `p11_consts.jl:79` guards its entire Tier-1 block on
  `if !isdefined(@__MODULE__, :P11_DEV_SEED)`, and `p12_consts.jl:108` legitimately **binds**
  `P11_DEV_SEED` (it forbids the Phase-11 stream by name, as R-5 requires). Loading the Phase-12
  constants first therefore makes the Phase-11 Tier-1 block a silent no-op — `LAMBDA_MIN` is never
  defined — and `p11_consts.jl`'s Tier-2 block, which is guarded on a *different* sentinel and so
  still runs, dies with `UndefVarError: SC2_SPEARMAN_ATTENUATION not defined`. Measured, both
  orders, not inferred.
- **Fix:** Inverted the order — Phase-11 chain first (it pulls `p11_consts.jl`, the frozen prior
  and the cache layer), Phase-12 pre-registration second — and documented the reason in a header
  block so the next `run_p12_*` runner does not rediscover it. The re-binding of the shared seed
  names is benign: identical literals in both files, and `p12_consts.jl:424-428` already records
  the one deliberate exception (`Z_TWO_SIDED_90`) as legal under Julia 1.12. Verified on Julia
  1.12.6. **No constants file was edited.**
- **Files modified:** `spike/validation/run_p12_eps_ridge.jl`
- **Commit:** a878ce1

**2. [Rule 2 - Missing critical check] `@assert isdir(dir)` provably cannot fire, so an absent
pool would not have been caught**

- **Found during:** Task 1
- **Issue:** The plan instructs "assert `isdir(dir)`" as the refusal that protects against a
  missing pool. But `p11_pool_dir` → `open_or_invalidate` (`spike/data/cache.jl:162`) ends in
  `isdir(dir) || mkpath(dir)`: the resolver **creates** the directory when it is absent, so the
  assertion is satisfied by the very call it is meant to check. An absent pool presents as a
  freshly-created **empty** hash directory, not as a missing one, and the failure would only have
  surfaced later as `load_p11_pool`'s generic "no shard_*.jld2 found".
- **Fix:** Kept the `isdir` assertion as instructed (with an honest short message) and added the
  load-bearing check next to it: enumerate `shard_*.jld2` under the resolved directory and refuse,
  with the plan's script-path-not-function-name wording, if the list is empty. The comment records
  why the first assertion cannot fire, so the redundancy is not mistaken for belt-and-braces.
- **Files modified:** `spike/validation/run_p12_eps_ridge.jl`
- **Commit:** a878ce1

### Additions beyond the plan's literal key list (no behaviour changed)

The artifact carries a few keys the plan did not name, all reported-only: `eps_ratio_primary`,
`primary_arm`, `control_rho_ratio_by_arm`, `ridge_residual_shrinkage_by_arm`, `vacuous_by_arm`,
`vacuous_shrinkage_floor`, `eps_prior_mean_fit`, `eps_prior_sd_fit`, `eps_prior_sd_analytic`,
`rho_bin_edges`, `theta_row_eps`, `theta_row_rho`, `test_fraction`, and a per-arm `arms` dict.
They exist because the plan requires *both arms* reported while the verify reads *scalar*
`ridge_residual_shrinkage` / `vacuous` / `control_rho_ratio` keys; the scalars are the primary
(`all128`) arm and the `_by_arm` dicts carry both. `eps_ratio` is itself the two-arm Dict, which
the verify prints. **No key named `shrinkage` exists** (asserted by the plan's own verify).

## Assumption Drift (advisory)

**1. The plan predicted the ε signal would strengthen with |ρ|; it weakens.**
- **Found during:** Task 2
- **Planned:** "the ε perturbation grows with |ρ| by physics (magnification mismatch degrades an
  existing correlation; at ρ ≈ 0 there is nothing to degrade)" — the stated rationale for the
  5-bin breakdown.
- **Actual:** the ratio is *lowest* (best recovery) at |ρ| ∈ [0, 0.2) at 0.9434 and rises
  monotonically to 1.0078 at |ρ| ∈ [0.8, 1.0), i.e. worse than the prior mean in the top bin.
- **Why it matters:** a reader taking the plan's rationale at face value would expect the
  strongest ε evidence exactly where the measurement is weakest. Recorded, not explained; the two
  candidate readings are listed above and neither was tested. Non-blocking; nothing was re-run.

**2. The "second arm isolates the mask rows" assumption is void on this pool.**
- **Found during:** Task 2
- **Planned:** "Report both — if only (a) works, the mask rows are carrying the signal, which is
  itself informative."
- **Actual:** the mask rows are identically 1.0 across all 50 000 columns, standardize to exactly
  zero, and cannot carry anything. The arms are the same regression by construction, so the
  branch cannot fire on Phase-11 data.
- **Why it matters:** the bit-identical arm columns would otherwise read as a copy-paste bug.
  They are not; the degeneracy was measured directly (per-row sd = 0.0 on all 64 mask rows).

## Deferred Issues

None.

## Known Stubs

None.

## Threat Flags

None. This plan opened no network path, no auth path, no file-write path beyond one tracked
artifact under `spike/validation/`, and no schema at a trust boundary.

## Self-Check: PASSED

- `spike/validation/run_p12_eps_ridge.jl` — FOUND
- `spike/validation/p12_eps_ridge_report.jld2` — FOUND (tracked, committed)
- commit `a878ce1` — FOUND in `git log`
