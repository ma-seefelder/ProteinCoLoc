---
status: complete
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 10
subsystem: pre-registration (Tier-2 probe-derived constant)
tags: [tau, resolution-probe, d-06, d-04, d-01, two-commit-discipline, measured-not-chosen]

# Dependency graph
requires:
  - phase: 13-three-hypothesis-amortized-bayes-factor
    plan: 01
    provides: "spike/p13/consts.jl Tier-1 block — the frozen probe spec (P13_TAU_DELTA_GRID, P13_TAU_AUC, P13_TAU_R, P13_TAU_STATISTIC, P13_TAU_DESIGN, P13_TAU_BOTH_DIRECTIONS, P13_TAU_BOOTSTRAP_B, P13_TAU_REFERENCE_LAMBDA_RULE, P13_TAU_ABORT_EXTEND_GRID) committed with P13_TAU deliberately ABSENT"
  - phase: 13-three-hypothesis-amortized-bayes-factor
    plan: 05
    provides: "spike/p13/tau_probe.jl — tau_curve, tau_from_curve, simulator_provenance_guard"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 02
    provides: "the simulator surgery the probe MUST measure against: theta arity 8 with chromatic_eps, composed warp in stage 6, widened SHIFT_PRIOR half-width 3.0"
provides:
  - "spike/p13/run_tau_probe.jl — the reported runner: locked-spec banner, runtime post-Phase-11 hard stop, whole-curve computation, atomic persist-before-verdict"
  - "spike/p13/tau_probe_report.jld2 — the complete A(delta) curve in BOTH directions with bootstrap medians, the entire locked spec, the seed/counter, both provenance shas and the consts fingerprint"
  - "MEASURED: P13_TAU = 0.15, the first grid delta clearing the pre-registered bar, at A(tau) = 0.916356 against P13_TAU_AUC = 0.9"
  - "P13_TAU_REFERENCE_LAMBDA = 3.0 — tau is really tau-of-lambda, so the reference travels with the number"
  - "the Tier-2 guard block in spike/p13/consts.jl, appended in a SEPARATE commit after the artifact"
affects: [13-11-datagen, 13-12-gate, 13-13-alpha-series, 13-16-realimage, 13-14-report]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Two-commit discipline for a probe-derived constant: the specification commit must contain NO value, the artifact commit must contain NO constant, and the constant commit must add only the appended Tier-2 block — git history IS the pre-registration audit trail"
    - "A Tier-2 append never edits Tier 1: the retired Tier-1 self-check is CONVERTED INTO A COMMENT quoting the commit sha at which it held, rather than deleted, so the original intent survives in the file"
    - "A live provenance assertion placed OUTSIDE the Tier-1 guard (so it re-runs on every include) catches a constant set from anywhere other than its own appended block"

key-files:
  created:
    - spike/p13/run_tau_probe.jl
    - spike/p13/tau_probe_report.jld2
  modified:
    - spike/p13/consts.jl
    - spike/test/test_p13_consts.jl

key-decisions:
  - "tau = 0.15 was MEASURED, not chosen: the whole 7-point curve was computed and persisted BEFORE tau_from_curve was ever evaluated, so no verdict could steer the measurement"
  - "The grid was NOT extended and the bar was NOT relaxed — neither was needed (the bar cleared at grid point 6 of 7), and P13_TAU_ABORT_EXTEND_GRID = false stood unused"
  - "The max-over-directions rule mattered: at tau the negative direction gave 0.916356 and the positive 0.905794. Both cleared 0.9, so the asymmetry did not change the verdict, but the reported tau is the max-direction value as pre-registered"
  - "tau = 0.15 sits ABOVE the ghat knot spacing near zero (0.0825), so tau_sub_knot = false and the sub-knot caution did not fire"

patterns-established:
  - "Persist-before-verdict: the artifact is written and integrity-checked before the pass/fail statistic is computed, so a crash cannot lose an expensive measurement and a verdict cannot retroactively shape what was recorded"
  - "The reported artifact is self-describing: it carries its own locked spec key-by-key, so a reader can re-read tau at any bar without trusting the runner"

requirements-completed: [D-06, D-04, D-01]

# Metrics
duration: ~40min (28.6 min of it the probe curve itself)
completed: 2026-07-27
---

# Phase 13 Plan 10: Reported τ Resolution Probe and Tier-2 Freeze Summary

**The pre-registered D-06 resolution probe ran once against a runtime-verified post-Phase-11
simulator and MEASURED τ = 0.15 — the smallest ρ difference the frozen 8×8 patch-correlation
summary can separate at the pre-registered bar (A(τ) = 0.916356 against P13_TAU_AUC = 0.9) — and the
value was frozen in an appended Tier-2 block in a commit that lands strictly after the artifact
commit, which is the auditable evidence that τ was measured before it was locked.**

## Performance

- **Duration:** ~40 min (the A(δ) curve alone took 1713.49 s = 28.6 min: 7 grid points × 2 directions × 400 draws, each a full simulate-and-summarize at an F5-mixture image size)
- **Completed:** 2026-07-27T21:53:41+02:00
- **Tasks:** 3 (Task 1 checkpoint, Tasks 2-3 auto)
- **Files modified:** 4

## Probe provenance

Recorded here as Task 1 required, and independently re-read from
`spike/p13/tau_probe_report.jld2` at close-out.

| Item | Value |
|---|---|
| Simulator guard | `post_p11 = true` |
| Guard reason (verbatim) | `POST-Phase-11 simulator confirmed: theta arity 8 with chromatic_eps present, shift-prior half-width 3.0 == expected 3.0.` |
| θ arity measured | 8 (with `chromatic_eps` present) |
| Shift-prior half-width measured | 3.0 (== expected 3.0) |
| Repo HEAD sha at the probe run | `cecb69c58770a5b45bec8151c657d9024b999f80` |
| Simulator-directory commit sha | `78dc37f517ad1b8e1e71afa963691f7ddefe5f63` |
| `consts.jl` fingerprint (sha256) | `640ec90edb80ef62eacbbf64c7b40cb9ca87c59d28a3a7438ef9d382b0aca4cf` |
| Runner fingerprint (sha256) | `d05057be5d3fc2bdfd491b2a932a5ad5a19fe81f0fd784974f9901a2d11586bd` |
| Measurement fingerprint (sha256) | `06eb1bc8fce8e9bf479cf20f314f100f1d0e1bcdbf39cc1fe01ea8d43921c1e6` |
| Artifact generated | `2026-07-27T19:39:58.238Z` |
| Reported stream | `P13_DEV_SEED = 185851505` at `P13_TAU_COUNTER = 1` |
| Arm keys (disjoint, unpaired by construction) | reference `2685821657552978796`, contrast `13507013436369108857` |
| `repo_dirty_at_run` | `true` — recorded, not hidden (see Issues) |

## The measured A(δ) curve, both directions

The whole curve is a reported artifact: it was computed and persisted in full before
`tau_from_curve` was evaluated, so a reader can re-read τ at any bar. Every grid point drew the
pre-registered `P13_TAU_R = 400` per arm, with **zero degenerate draws** at every δ.

| δ | A_neg | A_pos | **A (max)** | A bootstrap median | clears 0.9 |
|---|---|---|---|---|---|
| 0.02 | 0.529025 | 0.592806 | 0.592806 | 0.592806 | no |
| 0.03 | 0.566481 | 0.625869 | 0.625869 | 0.625184 | no |
| 0.05 | 0.639375 | 0.689925 | 0.689925 | 0.688956 | no |
| 0.075 | 0.723788 | 0.760544 | 0.760544 | 0.761053 | no |
| 0.10 | 0.800338 | 0.819356 | 0.819356 | 0.821203 | no |
| **0.15** | **0.916356** | 0.905794 | **0.916356** | 0.917831 | **YES** |
| 0.20 | 0.978556 | 0.953613 | 0.978556 | 0.978541 | yes |

**τ = 0.15** — grid index 6 of 7, `tau_status = :measured`, `tau_measured = true`.

Two things worth reading off this curve:

1. **The bar cleared with room, and one point to spare.** δ = 0.20 also clears, so τ is not pinned
   at the edge of the grid; the abort criterion (`P13_TAU_ABORT_EXTEND_GRID = false`) was never
   approached. The curve is monotone in δ in both directions, which is the sanity property a
   resolution curve must have.
2. **The both-directions rule earned its place.** The prior atoms are asymmetric by roughly 2.5×,
   and the two directions genuinely differ (at δ = 0.15, 0.916 vs 0.906; at δ = 0.20, 0.979 vs
   0.954). Here both directions clear the bar so the max-rule did not change the verdict — but at
   δ = 0.15 the positive direction alone would have cleared by only 0.006, and a single-direction
   probe would have been reporting a near-miss as a comfortable pass.

**Reference λ:** `P13_TAU_REFERENCE_LAMBDA = 3.0`, the realization of
`P13_TAU_REFERENCE_LAMBDA_RULE = :widest_rung` against Phase-11's `LAMBDA_MAX = 3.0`. τ is really
τ(λ) — measured at the WIDEST registration uncertainty, i.e. the hardest rung — so the reference
travels with the number and downstream code must not quote τ without it.

**Sub-knot check:** `ghat_knot_spacing_near_zero = 0.0825`, and τ = 0.15 > 0.0825, so
`tau_sub_knot = false`. The printed caution did not fire: τ is not finer than the calibration
knots that define ĝ near zero, so the dead zone is representable in the simulator's own ρ mapping.

## Task Commits

The two-commit discipline is the whole point of this plan, and it is visible in history:

1. **Task 1: Confirm the Phase-11 simulator surgery has merged** — checkpoint, no commit. The guard
   output and both shas are recorded above; they were consumed as literals by Task 3.
2. **Task 2: Write and run the reported probe, persisting the whole curve** — `581602a`
   (`feat(13-10): run the reported D-06 resolution probe and persist the whole A(delta) curve`),
   2026-07-27 21:41:39 +0200. Adds `run_tau_probe.jl` (330 lines) and the 15343-byte artifact.
   **Contains no τ constant.**
3. **Task 3: Append the Tier-2 τ block and update the pre-registration assertions** — two commits,
   because the constant and its downstream assertions are separable:
   - `bfba6ac` (`feat(13-10): append the Tier-2 block with the MEASURED P13_TAU = 0.15`),
     2026-07-27 21:46:26 +0200 — the appended guard block.
   - `45f5a44` (`test(13-10): flip the three downstream tau-is-absent assertions to the measured
     tau`), 2026-07-27 21:53:41 +0200 — the test-side flip.

Ordering verified: `581602a` (artifact, no constant) → `bfba6ac` (constant) → `45f5a44` (tests).
The artifact provably existed at a commit containing no τ value, so no reader has to take the
measurement order on trust.

## Files Created/Modified

- `spike/p13/run_tau_probe.jl` — the reported runner: `="^78` banner echoing every locked knob by
  name and value, stream-disjointness asserts against `P13_FIX_SEED` and `_p13_forbidden()`, a
  runtime `error` from `simulator_provenance_guard()` quoting Pitfall 10, whole-curve computation,
  atomic `.tmp` → reopen-integrity-check → `mv(force = true)` persist **before** any verdict.
- `spike/p13/tau_probe_report.jld2` — the self-describing artifact: per-δ `auc_neg`/`auc_pos`/`auc`/
  `auc_median_boot`/`n_*`, the whole locked spec key-by-key, seed + counter, both shas, three
  sha256 fingerprints, and a UTC `generated` stamp.
- `spike/p13/consts.jl` — the appended Tier-2 guard block (`P13_TAU`, `P13_TAU_MEASURED_AUC`,
  `P13_TAU_PROBE_SHA`, `P13_TAU_SIMULATOR_SHA`, `P13_TAU_PROBE_ARTIFACT`,
  `P13_TAU_REFERENCE_LAMBDA`) plus its two executable self-asserts, and the live
  provenance assertion at line 745 that sits OUTSIDE the Tier-1 guard.
- `spike/test/test_p13_consts.jl` — `"tau is Tier-2 and absent"` replaced by
  `"tau is measured and locked (Tier 2)"` (8 assertions).

## Decisions Made

- **The retired Tier-1 self-check was converted, not deleted.** `@assert !isdefined(@__MODULE__,
  :P13_TAU)` is now a comment block (consts.jl:716-727) recording that the assertion HELD and
  naming the commit that proves it: `git show c42cc8e:spike/p13/consts.jl` carries the live line.
  The historical statement survives in the file where a reader will find it, which a deletion would
  not have achieved. No other Tier-1 line changed.
- **A live provenance tripwire replaced it** (consts.jl:740-746): `P13_TAU` defined *without*
  `P13_TAU_PROBE_ARTIFACT` now errors. It sits outside the Tier-1 guard so it re-runs on every
  include, catching a τ hand-edited in, set from a REPL, or injected by a caller — the failure mode
  the deleted assertion used to cover, in the only form that still works once Tier 2 legitimately
  supplies a value.
- **`P13_TAU_REFERENCE_LAMBDA = 3.0` is stored as its own constant** rather than left implicit in
  `:widest_rung`. The rule resolves against a Phase-11 constant; if Phase 11's `LAMBDA_MAX` ever
  moved, an unstored rule would silently re-point τ at a different measurement.

## Deviations from Plan

**1. Task 3 landed as two commits rather than one**

- **Found during:** Task 3
- **Issue:** The plan specified one commit for the Tier-2 append plus the assertion update. The
  consts change and the downstream test flip are independently reviewable, and three downstream
  assertions (not the one the plan anticipated) asserted τ-is-absent.
- **Fix:** Split into `bfba6ac` (consts) and `45f5a44` (tests).
- **Impact:** Strictly stronger than the plan. The binding requirement was that the constant land
  *after* the artifact; both commits do. No acceptance criterion is weakened — the criterion reads
  `git log --oneline -2 -- spike/p13/` for artifact-before-constant ordering, which holds.

---

**Total deviations:** 1 (commit granularity)
**Impact on plan:** None on substance. The pre-registration audit trail is intact and, if anything,
finer-grained.

## Issues Encountered

- **`repo_dirty_at_run = true`.** The probe ran with an unclean working tree (the concurrent
  Phase-12 planning work). This is **recorded in the artifact rather than suppressed**, which is
  the honest handling: the recorded `repo_sha` (`cecb69c`) does not fully describe the tree the
  probe ran in. It does not compromise τ — the two things τ depends on are pinned independently and
  more tightly than by the repo sha: the `consts.jl` sha256 fingerprint
  (`640ec90e…`) pins the spec byte-exactly, and the simulator sha (`78dc37f5…`) plus the runtime
  guard pin the measured world. The dirty paths were `.planning/` documentation, disjoint from
  `spike/`.
- **No `.tmp` residue.** The atomic write completed; `spike/p13/tau_probe_report.jld2.tmp` does not
  exist.

## Verification at close-out

Re-run in the primary checkout on 2026-07-28, at merge commit `a90036f`:

- `julia --project=spike spike/test/test_p13_consts.jl` → **137/137 pass, 0 fail**, including the
  8-assertion `"tau is measured and locked (Tier 2)"` testset.
- `p13_tau() === P13_TAU === 0.15`; `P13_TAU in P13_TAU_DELTA_GRID`;
  `P13_TAU_MEASURED_AUC = 0.916356 >= P13_TAU_AUC = 0.9`.
- The artifact's recorded grid equals the frozen `P13_TAU_DELTA_GRID` exactly — 7 points, no
  extension: `(0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2)`.
- `length(curve) == length(P13_TAU_DELTA_GRID) == 7`.
- `src/`, `spike/Project.toml` and `spike/Manifest.toml` byte-unchanged.

## Why this SUMMARY was written at close-out rather than by the executing agent

The executing agent committed all three of this plan's commits inside worktree
`agent-aca63027b2785bf16` and then hit a session limit before writing the SUMMARY. Those commits
were merged into `gsd/v2.0-milestone` as `a90036f`. The resuming orchestrator hit the
`safe_resume_gate` Case A condition (production commits present, SUMMARY absent) and chose the
`close out manually` recovery path rather than re-dispatching, because **re-running the probe would
have been a second reported measurement of a pre-registered quantity** and would have had to be
declared as spending the D-04 iteration allowance that is reserved for the stratification switch.
Every number above is re-read from the committed artifact and the committed consts file, not
reconstructed from the agent's narration. `status: complete` is set deliberately and explicitly.

## Next Phase Readiness

- **τ exists, so labelled data can now be generated.** `three_way_label` reads τ through
  `p13_tau()`, which no longer errors — this was the gate that made plans 13-11 onward
  unrunnable.
- **Downstream code must quote τ with its reference λ = 3.0.** τ was measured at the widest
  registration rung; presenting it as an unconditional resolution would overstate it at narrower λ.
- **The D-04 iteration allowance is UNSPENT.** One reported probe run, one counter consumed
  (`P13_TAU_COUNTER = 1` on `P13_DEV_SEED`).
- **A note for the report (13-14):** the dead zone τ = 0.15 is wide relative to the coloc scale.
  Whether that is a limitation of the phase or of the 8×8 summary is a reporting question, and the
  full A(δ) curve is persisted precisely so 13-14 can state the trade-off with numbers rather than
  adjectives.

---
*Phase: 13-three-hypothesis-amortized-bayes-factor*
*Completed: 2026-07-27 (SUMMARY written at close-out 2026-07-28)*
