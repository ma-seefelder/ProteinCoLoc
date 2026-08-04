---
phase: 14-decision-and-abstention-layer
plan: 13
subsystem: decision-and-abstention-layer
status: complete
tags: [sc1-f, real-images, d-03a, illustration-not-coverage, corpus-unfetched-by-design, sealed-holdout-untouched, read-only, abstain, prediction-held]
requires:
  - "14-09 (spike/p14/p14_conformal_report.jld2 — q-hat = 0.4498001628978108, LOADED not re-calibrated)"
  - "spike/p14/decide.jl (decide_coloc per-pair path, p14_load_bundle, p14_ood_input — the D-07 inherited tau, prior and provenance)"
  - "spike/p14/pools.jl (p14_ood_reference — the density null fitted on the Phase-13 net's own pool)"
  - "spike/p13/real_images.jl (load_real, real_tif, verify_real_readonly_digest — the FROZEN real-image read path, reused verbatim)"
  - "test/test_images/{positive,negative}/*_c{1,2,3}.tif (the six committed microscopy TIFFs, READ-ONLY)"
  - "corpus/manifest.csv (METADATA ONLY — tier, split and byte columns; no image row opened)"
provides:
  - "SC1-f on the D-03a substrate: six recorded per-file decisions with their abstention triggers, on the only real microscopy this project may touch"
  - "spike/p14/p14_real_images_report.jld2 — the six-row images table, the two decision units, the pre-stated prediction beside its outcome, and the computed corpus availability record"
  - "the recorded result that BOTH specimens ABSTAIN (:ood_fired) — the designed behaviour, and the prediction held"
  - "the computed corpus record: 30 simulated-secondary / 2 physical-primary, dev = 14 / eval = 16 / sealed_holdout = 2, 0 image bytes, :unfetched_by_design"
affects:
  - "the Phase-14 verdict: SC1-f is satisfied structurally (no coverage number, sealed holdout untouched) rather than by clearing a threshold — there is no threshold"
  - "spike/test/test_p14_decoupling.jl: the POSITIVE `test/test_images` assertion now fires (51 passing, up from 50), so the three sealed-holdout negatives are no longer vacuous"
  - "any manuscript claim about real-data behaviour: the D-03 bound is recorded as ABSENT, not merely loose"
tech-stack:
  added: []
  patterns: ["persist-before-print", "atomic .jld2 save with reopen integrity check", "run-time decoupling proof extended to test/ and corpus/", "read-only digest taken before and after and asserted equal", "prediction written down before the run and persisted beside its outcome", "breakdown counts DERIVED from the manifest, never hard-coded", "ignore-rule evidence located at run time rather than transcribed", "withdrawn figure expressed in prose so no numeral exists to leak"]
key-files:
  created:
    - "spike/p14/run_p14_real_images.jl (956 lines)"
    - "spike/p14/p14_real_images_report.jld2 (gitignored, reproducible)"
  modified: []
decisions: [D-03, D-03a, D-05, D-06]
validates: ["SC1-f", "SC2-a", "SC2-b"]
metrics:
  duration: "~50 min wall (each reported run of the runner: ~23-26 s)"
  completed: 2026-08-04
---

# Phase 14 Plan 13: SC1-f — The Decision Layer on the Only Real Microscopy Summary

**On both committed specimens the decision layer ABSTAINS with `abstain_reason = :ood_fired`, at a
summary density of 703.2995 against the fitted operating point 167.5446 (4.198×) — exactly the
behaviour predicted in writing before the run — and the provenance corpus is recorded as
`:unfetched_by_design` with a tier breakdown computed from its manifest (30 `simulated-secondary`,
2 `physical-primary`, both of them the untouched Phase-16 sealed holdout). This is an illustration,
not a coverage claim.**

## What ran

`spike/p14/run_p14_real_images.jl`, on the main working tree, on the house runner shape (banner,
step-0 run-time decoupling proof, locked-inheritance banner, atomic `_p14_real_save_report`,
smoke-mode redirect, persist-before-print, guarded `abspath(PROGRAM_FILE)` entry point).

- **Substrate:** the six committed microscopy TIFFs under `test/test_images/`, opened through
  Phase 13's frozen `load_real` / `real_tif` — no second real-image read path was authored.
- **Operative pair:** `P13_REAL_CHANNEL_PAIR = (2, 3)` (green/red, the two target proteins),
  passed EXPLICITLY at every call site. c1 is the DAPI/Hoechst nuclear counterstain and enters no
  comparison; it is still opened and digested, because the unit of record is the file.
- **Decision path:** `decide_coloc(bundle, img, control, [2,3]; ood_nulls, qhat, ...)` — the
  src-shaped per-pair entry point — at the 14-09 q-hat **loaded, never re-calibrated**, with
  `allow_unchecked_ood = false` passed explicitly, on the density null fitted on the Phase-13 net's
  own pool.
- **Two decision units**, matching Phase 13's real arm: `positive_as_sample` (control = negative)
  and `negative_as_sample` (control = positive).

## The six recorded decisions

The decision is a property of a two-channel acquisition compared against a control, so each file row
carries the decision of the unit in which **its condition was the sample**, plus `in_operative_pair`.

| file | ch | in pair | decision | trigger | OOD state | score / threshold | conformal | P(coloc) |
|---|---|---|---|---|---|---|---|---|
| `test/test_images/positive/positive_c1.tif` | c1 | false | **abstain** | `:ood_fired` | `:fired` | 703.2995 / 167.5446 | singleton `[:coloc]` | 0.99280 |
| `test/test_images/positive/positive_c2.tif` | c2 | true  | **abstain** | `:ood_fired` | `:fired` | 703.2995 / 167.5446 | singleton `[:coloc]` | 0.99280 |
| `test/test_images/positive/positive_c3.tif` | c3 | true  | **abstain** | `:ood_fired` | `:fired` | 703.2995 / 167.5446 | singleton `[:coloc]` | 0.99280 |
| `test/test_images/negative/negative_c1.tif` | c1 | false | **abstain** | `:ood_fired` | `:fired` | 703.2995 / 167.5446 | singleton `[:random]` | 0.00008 |
| `test/test_images/negative/negative_c2.tif` | c2 | true  | **abstain** | `:ood_fired` | `:fired` | 703.2995 / 167.5446 | singleton `[:random]` | 0.00008 |
| `test/test_images/negative/negative_c3.tif` | c3 | true  | **abstain** | `:ood_fired` | `:fired` | 703.2995 / 167.5446 | singleton `[:random]` | 0.00008 |

Cross-method channel (D-05), on both units: `classical_call = RANDOM`, `costes_p = 0.004975`
(`costes_call = :coloc`), `basis = :ghat_rho_true_scale`, `costes_seed = P14_DEV_SEED`. The fused
OOD trigger fires first, so the classical disagreement never reaches the decision — which is the
D-05 asymmetry doing its job.

Both units score the identical density 703.2995 because the shipped `ood_verdict` fuses over the
sample and the control summaries with `max`, and the two units are the same pair of acquisitions in
swapped roles.

## The prediction held

Written into the runner **before** the run, from 13-REPORT limit D, and persisted verbatim beside
the outcome:

> `predicted_behaviour = "ABSTAIN on both specimens (13-REPORT limit D: density 703.30 vs threshold
> 167.54, 4.198x)"` → `observed_decisions = [:abstain, :abstain]` → `prediction_held = true`

On the only real microscopy this project holds, the decision layer abstains — and that is the
**designed** behaviour under D-05/D-06, not a failure. The measured density reproduces 13-REPORT's
703.30 / 167.54 to four decimal places through an entirely independent composition path, which is
itself evidence that the Phase-13 and Phase-14 read chains are the same chain.

## The corpus availability record (metadata only)

Every number below was **computed from `corpus/manifest.csv`** at run time, not transcribed:

| field | value |
|---|---|
| `corpus_images_available` | 0 (total bytes 0 across all 32 rows) |
| `corpus_status` | `:unfetched_by_design` |
| `tier_breakdown` | `simulated-secondary` = 30, `physical-primary` = 2 |
| `split_breakdown` | `dev` = 14, `eval` = 16, `sealed_holdout` = 2 |
| `physical_truth_available` | `false` (0 unsealed physical-primary rows with bytes) |
| `d03_bound_status` | `:absent_not_loose` |
| `sealed_holdout_read` / `corpus_images_opened` | `false` / `0` |

The ignore-rule evidence is located at run time by name (`.gitignore:459, 461, 463, 465, 466`) and
the rule text persisted, so the "by design" reading rests on the file rather than on a claim. The
withdrawn D-03 miscoverage figure is recorded as **withdrawn** in prose and appears as **no numeral**
anywhere in the runner or the artifact.

## This is an illustration, not a coverage claim

Six TIFFs in two conditions, with no colocalization ground-truth labels, bounds nothing tightly.
`is_coverage_claim = false` and `illustration_note` are REQUIRED keys of the atomic save, so a report
that failed to say so is never written at all — and no rate, fraction or percentage is computed over
these six images anywhere in the runner.

## Structural guarantees (SC1-f has no threshold)

SC1-f's four failure conditions are all structural, and each is prevented by construction:

| failure condition | prevention | evidence |
|---|---|---|
| a coverage number is emitted | no rate is computed; `is_coverage_claim = false` is a required key | artifact |
| the sealed holdout is read | metadata only; the sealed image directory is never named, spelled or opened | `test_p14_decoupling.jl` checks (a)/(c) pass |
| the corpus reported as gone/empty | `corpus_status_note` contains `unfetched`, and neither `gone` nor `empty` | asserted, verified |
| six TIFFs presented as a coverage claim | banner, `illustration_note`, and the plan's own wording | artifact + printout |

Additionally: `git ls-files -s test/test_images` digest taken **before and after** the run and
asserted equal (`eaee22f9…`, unchanged); the six per-file sha256 digests recorded; `git status
--porcelain -- test` and `-- corpus` empty; `git diff HEAD -- corpus test src spike/Project.toml
spike/Manifest.toml` clean; `spike/p14/consts.jl` byte-unchanged.

## Verification (all observed, not inferred)

| check | result |
|---|---|
| `julia --project=spike -t auto spike/p14/run_p14_real_images.jl` | exit 0, artifact written, 6 rows |
| Task-2 artifact key/value assertion (11 keys + 4 values) | exit 0, printed the computed breakdowns |
| `julia --project=spike spike/test/test_p14_decoupling.jl` | **51 / 51 pass** (was 50 — the positive `test/test_images` assertion armed) |
| `git status --porcelain -- test` | empty |
| `git diff --exit-code HEAD -- corpus` / `git status --porcelain -- corpus` | exit 0 / empty |
| `grep -v '^\s*#' … \| grep -cE '0\.032\|1/31'` | `0` |
| `grep -v '^\s*#' … \| grep -ci 'corpus/data\|CORPUS_DATA_DIR\|open_sealed_holdout'` | `0` |
| `d["corpus_status_note"]` contains `unfetched`, not `gone`/`empty` | `true` / `false` / `false` |
| `git diff --exit-code HEAD -- spike/p14/consts.jl` | exit 0 — byte-unchanged |

## Deviations from Plan

### Auto-fixed / design resolutions

**1. [Rule 3 — resolved gap] The "six rows" vs "two specimens" reconciliation**
- **Found during:** Task 1.
- **Issue:** The plan requires a **6-row** `images` table each carrying a decision, and separately
  predicts the outcome for **both specimens**. A decision is a property of a two-channel acquisition
  compared against a control, so there are exactly two decision units and six files.
- **Resolution:** the artifact carries BOTH — a 6-row `images` table keyed on the file (condition,
  channel, repo-relative path, sha256, `in_operative_pair`) where each row carries the decision of
  the unit in which its condition was the sample, AND a 2-row `decision_units` table. A required
  `decision_unit_note` states the distinction in the artifact so the duplication cannot be misread
  as six independent measurements.
- **Commit:** `74266aa`.

**2. [Rule 2 — missing critical check] `corpus/` byte-equality asserted at run time**
- **Found during:** Task 2.
- **Issue:** the plan's acceptance criteria check the provenance tree in review; nothing checked it
  while the metadata read happened (the T-14-11 hole the lane guard exists to close for `src/`).
- **Fix:** step 0 now also asserts `git diff --quiet HEAD -- corpus` and an empty
  `git status --porcelain -- corpus`, plus the same pair for `test/` before any fixture is opened.
- **Commit:** `8b08743`.

**3. [Rule 2] The read-only digest asserted, not merely recorded**
- **Found during:** Task 1.
- **Issue:** the plan asks for sha256 digests of the six files; it does not ask for the tree-level
  proof to be *checked*.
- **Fix:** `verify_real_readonly_digest()` is taken before and after and asserted equal AFTER the
  artifact is persisted, so a violated read-only discipline fails loudly with the evidence intact.
- **Commit:** `74266aa`.

### Assumption Drift (advisory)

**A. The `.gitignore` line range the plan cites is stale.**
- **Planned:** record `.gitignore:451-453` as the evidence that the image bytes are unfetched by
  design (the range 14-CONTEXT D-03a also quotes).
- **Actual:** the rules now live at `.gitignore:459-466`, with the ignore rule itself at 461;
  lines 451-453 are Phase-13 training side effects.
- **Why it matters:** a transcribed line number goes stale the first time an unrelated rule is
  inserted above it, and this one already had. The runner therefore **locates the rules by name at
  run time** and persists their text and their line numbers, so the record moves with the file.
  Nothing was relaxed; the evidence got stronger.
- Advisory only — no gate, no decision, no threshold touched.

**B. Every classical channel and every conformal set is a singleton here, and the abstention is
carried entirely by the OOD channel.**
- **Planned:** the plan anticipates ABSTAIN and names the OOD score as the reason.
- **Actual:** confirmed, and additionally: the conformal hedge emits a **singleton** on both units
  (`[:coloc]` and `[:random]`), so the hedge itself saw nothing unusual. The two independent
  misspecification signals **disagree** on this substrate — the density null fires hard while the
  exchangeability-based conformal set does not. That is the more interesting reading (14-12's
  `crosstab_note` says so for the simulated arms) and it is recorded rather than editorialised.
- Advisory only.

## Named limits carried by this result

- The decision layer is built in the spike research lane and is **NOT shipped in v2.0** (D-01);
  `src/` is byte-unchanged and untracked-clean, asserted on every run.
- Only the **summary-density** OOD channel is wired; `(:pp, :noise)` are not. A `:fired` here is one
  channel firing, not three — and 14-12 measured that channel blind to `misspec_noise`.
- The conformal guarantee is **simulator-derived** and inherits the simulator's misspecification in
  full; the real-data bound D-03 intended is **absent, not merely loose**.
- The underlying net is an epoch-4 checkpoint of an early-overfitting run.
- SC1 and SC2 are AMENDED (D-02, D-05); the original wording must be cited alongside any number here.

## Known Stubs

None.

## Threat Flags

None. No new network endpoint, auth path, file-access pattern or schema at a trust boundary was
introduced; the one new read surface is a committed text file read as metadata, and the sealed
holdout was neither opened nor named.

## Self-Check: PASSED

- `spike/p14/run_p14_real_images.jl` — FOUND (956 lines).
- `spike/p14/p14_real_images_report.jld2` — FOUND (gitignored bulk artifact, reproducible).
- Commit `74266aa` (Task 1) — FOUND.
- Commit `8b08743` (Task 2) — FOUND.
