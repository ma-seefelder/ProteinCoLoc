---
phase: quick-260725-vl8
plan: 01
subsystem: spike/simulator (SIM-02 induced-μ calibration diagnostic)
tags: [calibration, diagnostic, honesty-fix, d-16, real-anchor, otsu-mask]
requires:
  - src/LoadImages.jl::_calculate_mask / _apply_mask! (called read-only via contract.jl)
  - spike/contract.jl::build_mci / patch_summary / induced_mu / load_tiff (unchanged)
provides:
  - "_real_anchor(): masked AND unmasked induced μ + surviving-patch count per fixture"
  - "neg_reachable derived from the MASKED negative-fixture μ"
  - "corrected, two-sided real-anchor provenance in the frozen ghat.jl header"
affects:
  - spike/NOTES.md §3 (still repeats the withdrawn reading — out of scope, see Follow-ups)
  - spike/test/test_simulator.jl:225 (same — out of scope)
  - .planning/phases/02-forward-simulator-summary-contract/02-03-SUMMARY.md (same — out of scope)
tech-stack:
  added: []
  patterns:
    - "diagnostic reports BOTH estimators rather than silently choosing one"
    - "in-place mutation ordering made explicit (_apply_mask! mutates mci.data)"
key-files:
  created: []
  modified:
    - spike/simulator/calibration.jl
    - spike/simulator/ghat.jl
decisions:
  - "neg_reachable is computed from the MASKED negative-fixture μ (the estimator comparable to the v1.0 analysis pipeline), and is explicitly documented as NOT evidence of physical reachability"
  - "the unqualified 'negative tail physically reachable = false ⇒ PRIOR-ONLY' claim is WITHDRAWN; reachability is NOT ESTABLISHED on n=2 fixtures with two oppositely-biased estimators"
  - "ghat.jl's anchor block is hand-patched, not regenerated — re-running calibration.jl would re-fit and change the frozen knots"
metrics:
  duration: ~25 min
  completed: 2026-07-25
  tasks: 3
  files: 2
---

# Quick 260725-vl8: Fix the Real-Anchor Diagnostic in spike/simulator Summary

The SIM-02 real anchor now reports both the unmasked estimator (what `patch_summary`, and therefore
the v2.0 training and inference path, actually sees) and the masked estimator (what the v1.0 analysis
pipeline sees via `_apply_mask!` + `_exclude_zero`) with surviving-patch counts, and the
negative-tail reachability overclaim that rested on the unmasked pair alone is withdrawn rather than
inverted.

## What Was Measured

Run on the two real fixtures in `test/test_images/` via the package's own `load_tiff` (D-16 —
faithful load, not an RGB→luminance reduction):

| fixture  | unmasked μ | patches | masked μ  | patches |
|----------|-----------:|--------:|----------:|--------:|
| positive | 0.3292     | 64/64   | −0.0729   | 21/64   |
| negative | 0.2481     | 64/64   |  0.1149   | 20/64   |

Exact values observed in the Task-1 verify run:
`positive = (unmasked_mu = 0.329163298494035, unmasked_n = 64, masked_mu = -0.07293570883861245, masked_n = 21)`,
`negative = (unmasked_mu = 0.24805233208940947, unmasked_n = 64, masked_mu = 0.11485341998017469, masked_n = 20)`.

These reproduce the pre-dispatch reference measurement exactly.

**`neg_reachable` = `false`** — but now derived from `anchor.negative.masked_mu` (0.1149, which is
> 0), not from the unmasked 0.2481, and now carrying an explicit note that the predicate is not by
itself evidence of physical reachability.

Two facts drove the honesty fix:
1. The **positive control measures NEGATIVE when masked** (−0.0729) and the ordering between the
   controls **inverts** between the two estimators (unmasked: positive > negative; masked:
   positive < negative).
2. Only ~20 of 64 patches clear the ≥15-survivor floor under masking, so the masked number rests on
   roughly a third of the image.

## What Was Built

**Task 1 — `spike/simulator/calibration.jl`, `_real_anchor()` (commit `a891ca5`)**
Each fixture is loaded once and measured twice off the same `MultiChannelImage`: unmasked first,
then `_apply_mask!(mci, _calculate_mask(mci))` (mirroring `src/utils.jl:86` verbatim, using the Otsu
thresholds `build_mci` already stored — not a re-thresholding), then measured again. The ordering is
load-bearing because `_apply_mask!` mutates `mci.data` in place; that is stated in a comment. A
local `_anchor_measure(mci)` returns `(μ, n)` where `n` counts the 8×8 summary entries that are
non-`missing` **and** finite — i.e. exactly the patches entering the mean, out of 64 —
by calling `patch_summary`, never reimplementing it. Return shape is a NamedTuple per fixture:
`(unmasked_mu, unmasked_n, masked_mu, masked_n)`.

`_anchor_na()` supplies the same nested shape (`NaN` μ, `0` counts) for the `catch` fallback in
`calibrate()`, so a fixture load failure cannot throw downstream. `neg_reachable` now reads
`anchor.negative.masked_mu`.

**Task 2 — emit template + printed report (commit `231e8d7`)**
The single `real anchor (D-16)` line and the `negative tail` ternary in `emit_ghat`'s template are
replaced by a block carrying all four μ values, all four patch counts, both bias directions, and the
NOT-ESTABLISHED verdict; the `println` report at the script entry mirrors the same content compactly.
Every other evidence line (generator, sweep, induced-μ range, monotonicity, SIM-02 match,
size-invariance, σ/τ, ν) and the `GHAT_*` / `ghat()` emission are byte-identical.

**Task 3 — `spike/simulator/ghat.jl` header (commit `151ad79`)**
Comment-only hand patch. The corrected anchor block is copied verbatim from the Task-2 template, so
a future regeneration is a no-op on that block. A `HAND-PATCH PROVENANCE (2026-07-25)` note sits
directly under the `AUTO-GENERATED — DO NOT EDIT BY HAND` banner recording that the edit was
deliberate, why it could not be regenerated (re-running `calibration.jl` re-runs a stochastic sweep
and re-fits the map, changing the frozen knots), and that the knots and `ghat()` are byte-identical.

## The Corrected Claim (honesty)

The new text does **not** replace one overclaim with its opposite. It records that **both**
estimators are biased, in **opposite** directions:

- **UNMASKED is biased POSITIVE** — background pixels are dark in both channels and co-vary, so the
  shared background alone lifts the per-patch correlation.
- **MASKED is biased NEGATIVE** — each channel is thresholded independently and `_exclude_zero` then
  keeps only pixels bright in **both**, a selection on both variables that restricts the joint range
  and induces spurious negative correlation. This is the known weakness of thresholded Pearson and
  the reason Manders/Costes coefficients exist. Masking also leaves only ~20 of 64 patches above the
  ≥15-survivor floor.

Verdict recorded: **negative-tail reachability is NOT ESTABLISHED** either way — n = 2 fixtures and
two oppositely-biased estimators cannot settle it. The former unqualified
`"physically reachable = false ⇒ PRIOR-ONLY"` rested on the unmasked numbers alone and is
**WITHDRAWN**. The SIM-02 consistency claim stays scoped to `[GHAT_MU_MIN, GHAT_MU_MAX]`; that
scoping is unchanged and does not depend on this anchor.

## Verification (all run, output observed)

| Check | Command | Observed |
|-------|---------|----------|
| Task 1 | `julia --project=spike -e 'include("spike/simulator/calibration.jl"); a = _real_anchor(); …'` | printed the four μ / four counts above; all asserts held (64/64 unmasked, 21 and 20 masked, `positive.masked_mu < 0 < positive.unmasked_mu`); **`REAL-ANCHOR OK`** |
| Task 2 (template) | probe script rendering `emit_ghat(r; path=tempdir()/ghat_probe.jl)` with a stub `r` | rendered block contains `UNMASKED`, `MASKED`, `0.3292`, `0.2481`, `-0.0729`, `0.1149`, `64/64`, `21/64`, `20/64`, `NOT ESTABLISHED`, `WITHDRAWN`, `Manders`; the old sentence absent; **`EMIT-TEMPLATE OK`** |
| Task 2 (report) | probe evaluating the 5 new `println` source lines against a stub `r` | all five rendered; same 12 substrings present; **`REPORT-LINES OK`** |
| Task 3 (diff) | `git diff -U0 -- spike/simulator/ghat.jl \| grep -E '^[+-][^+-]' \| grep -vE '^[+-]\s*#'` | **empty** — every added/removed line is a comment |
| Task 3 (load) | `julia --project=spike -e 'include("spike/simulator/ghat.jl"); …'` | `GHAT_MU_MIN == -0.67976`, `GHAT_MU_MAX == 0.847149`, 25 + 25 knots, `ghat(0.0) = -0.12004360606355068`; **`GHAT-HEADER OK`** |
| Decoupling | `git status --porcelain src/ spike/contract.jl` | **empty** |
| Scope | `git diff --name-only HEAD~3 HEAD` | exactly `spike/simulator/calibration.jl`, `spike/simulator/ghat.jl` |

The full `calibrate()` sweep was deliberately **NOT** run (it would re-fit and change the frozen
knots); every verify `include`s the file, and the `PROGRAM_FILE` guard keeps the script entry inert.
The script-entry `println` block is therefore not exercised as a whole run — the five new lines were
verified individually by evaluating their actual source text against a stub `r`, which is recorded
above rather than implied.

## Deviations from Plan

**None affecting scope or behavior.** Two implementation details worth recording:

1. **[Rule 3 — blocking] Task 2's verify command was a sketch.** The plan's command referenced a
   nonexistent `GHAT_MU_KNOTS_STUB()` and explicitly authorized defining "whatever minimal stub `r`
   the current `emit_ghat` signature needs". A probe script was written to the scratchpad supplying
   a literal 3-element `mu_knots`/`rho_knots` stub and the real `_real_anchor()` output, writing to
   `tempdir()`. `spike/simulator/ghat.jl` was never written by the probe.
2. **Added verification beyond the plan.** The plan's Task-2 verify covered only the emit template,
   leaving the printed `println` report unexercised. A second probe evaluates those five source
   lines directly against a stub `r`, so the "both numbers in the printed report" must-have is
   backed by observed output rather than by inspection.

## Assumption Drift (advisory)

None material. The plan predicted masked patch counts of "≈21/20" and masked μ of −0.0729 / 0.1149;
the run reproduced 21/20 and those μ values exactly, so no re-numbering of the ghat.jl header was
needed.

## Out-of-Scope Follow-ups (identified by the planner, NOT edited)

Three artifacts still repeat the now-withdrawn PRIOR-ONLY reading and should be corrected in a
separate change:

1. **`spike/NOTES.md` §3 "Real anchor (D-16)"** — lines 225, 242–247 and 309 state "the negative
   μ-prior tail is NOT physically reachable", "documented PRIOR-ONLY", and "real fluorescence data
   never showed anti-correlation", all resting on the unmasked pair.
2. **`spike/test/test_simulator.jl:225`** — comment: "the negative μ-prior tail past `GHAT_MU_MIN`
   is prior-only (real anchor > 0)".
3. **`.planning/phases/02-forward-simulator-summary-contract/02-03-SUMMARY.md`** — lines 39, 55, 70
   and 97 record "positive μ=0.329, negative μ=0.248 (both >0) → real anti-correlation not
   physically realized" as a frozen decision and headline claim.

Note that (3) is a historical phase summary; correcting it is a records question (amend vs. append a
correction note), not a code change.

## Known Stubs

None. The two probe scripts live in the session scratchpad, not in the repo.

## Threat Flags

None. No new network endpoint, auth path, file-access pattern, or schema change. `src/` was not
touched (`git status --porcelain src/` empty); the shipped v2.0 summary/training/inference path was
not touched; the frozen fit is byte-identical.

## Commits

| Task | Commit | Message |
|------|--------|---------|
| 1 | `a891ca5` | `fix(vl8-01): measure the real anchor masked AND unmasked, derive neg_reachable from the masked mu` |
| 2 | `231e8d7` | `docs(vl8-02): carry both anchor estimators and the two-sided bias into the emit template and report` |
| 3 | `151ad79` | `docs(vl8-03): correct the frozen ghat.jl real-anchor provenance header (comment-only)` |

## Self-Check: PASSED

- Files exist: `spike/simulator/calibration.jl`, `spike/simulator/ghat.jl`, this SUMMARY.
- Commits exist: `a891ca5`, `231e8d7`, `151ad79`.
- `git status --porcelain src/ spike/contract.jl` empty; the only other dirty paths belong to the
  concurrent Phase-11 executor (`spike/data/generate.jl`, `spike/npe/*`, `spike/test/*`) and the
  pre-existing `.planning/HANDOFF.json` deletion — none were staged by this plan.
