---
phase: 12-spatial-colocalization-map
plan: 02
status: complete
subsystem: planning-artifacts
tags: [pre-registration, success-criteria, requirements, roadmap, frozen-spec]
requires:
  - ".planning/ROADMAP.md Phase-12 entry (SC1/SC2/SC3 originals)"
  - ".planning/phases/12-spatial-colocalization-map/12-CONTEXT.md (D-02, D-05..D-13, R-4, R-7, R-8, V-1..V-3, K-1, S-4)"
  - ".planning/phases/12-spatial-colocalization-map/12-CONSTANTS-FOR-CONFIRMATION.md (CLOSED 2026-07-29)"
  - ".planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md (the §0 shape and the legitimacy test)"
  - "spike/validation/p12_consts.jl (Tier-1, for the constants quoted)"
provides:
  - "12-SC3-AMENDMENT.md — frozen pre-declaration of amended SC1/SC2/SC3 plus the §5 Δρ decision"
  - "SPAT-01 … SPAT-09 minted and mapped to Phase 12"
  - "ROADMAP Phase-12 AMENDED annotation + corrected Depends-on line"
affects:
  - "every later Phase-12 plan (they cite the amendment; none may edit it)"
  - "12-12 (structural-identity guard), 12-16 (Stage-2 gate arm), 12-18 (SBC), 12-19 (real arm), 12-20 (guards)"
tech-stack:
  added: []
  patterns: ["frozen pre-registration amendment (13-SC2-AMENDMENT precedent)", "counted-not-inherited corpus audit"]
key-files:
  created:
    - ".planning/phases/12-spatial-colocalization-map/12-SC3-AMENDMENT.md"
  modified:
    - ".planning/REQUIREMENTS.md"
    - ".planning/ROADMAP.md"
decisions:
  - "SC1's DeepSet half dropped: permutation-invariance over patches IS the exchangeable pooling the phase replaces"
  - "AbstractGPs.jl not used; both prior arms are dense 64×64 stdlib LinearAlgebra at G=8, parametrized by induced lag-1 r₁"
  - "SC2 scoped to a spike-local coloc_map; src/results.jl:180-191 stays a byte-unchanged sketch (13-08 precedent)"
  - "SC3 amended to leave-region-out predictive coverage vs a matched RETRAINED ablation, two-sided (calibration AND log-score)"
  - "The real-image arm is reported, not gated; the Stage-2 gate rests on 12-16's simulated arm at N ≥ 271"
  - "All three maps ship (sample, control, difference); region_delta_rho stays primary; the guard is a structural identity, not a range"
  - "Δρ is NOT computed on the real arm — the two specimens are not exchangeable; SC3 is single-stack and therefore unaffected"
requirements: [SPAT-06, SPAT-09]
metrics:
  duration_min: 34
  tasks_completed: 2
  files_changed: 3
  completed: 2026-07-29
---

# Phase 12 Plan 02: Freeze the SC Amendment and Mint the SPAT IDs — Summary

Froze `12-SC3-AMENDMENT.md` at sha `af8e1ca9` with a counted, verified §0 disclosure that **no
Phase-12 result exists**, correcting three defective ROADMAP success criteria and encoding the
2026-07-28 Δρ decision as §5; minted SPAT-01 … SPAT-09 so every Phase-12 plan's `requirements` field
resolves.

## What Was Built

### Task 1 — `12-SC3-AMENDMENT.md` (commit `ac9afdb`)

459 lines. Six `##` headings: §0 provenance plus **exactly five** numbered amendment sections
(§1–§5), plus one unnumbered close carrying the legitimacy test and the freeze statement.

- **§0** — provenance table with the write-time sha, plus a second table recording that all four
  emptiness checks were **run, not assumed** (see Verification Evidence below).
- **§1** — SC1: the DeepSet is dropped (permutation-invariance over patches *is* the pooling the
  phase replaces); CNN only; D-02 recorded so the scope is unmistakable — summary **content** is
  unchanged, so `patch_summary`, the OOD Mahalanobis reference, `local_coloc_map`, the shipped
  artifact and `DATAGEN_HASH_SRC_FILES` all stay valid. R-7's lean topology written in literally.
- **§2** — SC1: `AbstractGPs.jl` is not used; both arms are dense 64×64 stdlib `LinearAlgebra`,
  because adding the package forces a `Pkg.resolve()` that could dislodge the NeuralEstimators 0.2.1
  pin (the Phase-7 Wave-0 co-resolution gate already failed once on this class). R-8 recorded: both
  arms parametrized by induced lag-1 r₁, plus the `sqrt.(diag(Σ))` rescale D-05 omits.
- **§3** — SC2 scoped to a spike-local `coloc_map`, three reasons (the 2026-07-24 GO; the input/output
  contract mismatch — 128-row + MLP + 7-marginal flow versus `(G,G,2,K)` + 72-marginal flow; the
  content-hashed Release-hosted bundle), following the 13-08 precedent.
- **§4** — SC3 amended to leave-region-out predictive coverage against a matched **retrained**
  ablation, with clauses (a)–(h) all pre-declared: not parameter coverage; retrained not pinned;
  "no spatial borrowing, full nuisance and global borrowing" plus a third prior-only floor arm;
  two-sided with `P12_STAGE2_LOGSCORE_MIN` = 0.02 nats and an explicit "a coverage-only win is NOT a
  pass"; the Fisher-z noise model declared as a joint-validation assumption; six FILES / two DATASETS
  with the **image** as the independent unit; the S-4 chromatic-radial confound named in advance with
  its three guards; and (h) the real arm **reported, not gated**.
- **§5** — the Δρ decision: all three maps ship, `region_delta_rho` primary, **zero forward passes**
  extra at read time and a single-stack pool, the withdrawn-not-deleted 12-09 clause, the
  structural-identity guard replacing 12-12's range check, per-map sentinel semantics, the counted
  corpus audit, and the named limit that Δρ is not computed on the real arm because the two specimens
  are not **exchangeable**.

### Task 2 — SPAT IDs and the ROADMAP annotation (commit `092e12e`)

- `.planning/REQUIREMENTS.md`: a `### Spatial Colocalization Map (Phase 12 — feature expansion)`
  block with SPAT-01 … SPAT-09 in the existing checkbox format, nine matching rows in the
  ID-to-phase traceability table (status `Pending`), an updated Coverage line, and an **appended**
  `*Last updated*` footer (the 2026-07-02 line is retained, not overwritten).
- `.planning/ROADMAP.md` Phase-12 block only: `**Requirements**: TBD` → the nine IDs; an `**AMENDED**`
  line in the exact Phase-13 shape; and the `**Depends on**` line with the registration-aware-θ and
  serialization clauses struck as VOID (V-1, V-2, V-3), Phase 7 left as the substantive dependency,
  and `ca02b0e` recorded as this phase's baseline (K-1).

## Verification Evidence — what was actually run

Every claim below is an observed command result, not an inspection.

**§0 emptiness checks, at sha `af8e1ca905438ebfbd1f5fb96535c95947b6a533`:**

| Checked | Observed |
|---|---|
| `ls spike/validation/p12_*report*.jld2` | `No such file or directory` |
| `ls spike/data/cache/p12` | `No such file or directory` (cache holds `fixture`, `p11`, `p13`, one content-addressed pool) |
| `ls spike/npe/p12_*.jld2` | `No such file or directory` |
| any trained Phase-12 net under any name | `spike/npe/` holds only `trained_npe.jld2` and `p11_research_npe.jld2` |
| `git ls-files \| grep -i p12` | `spike/validation/p12_consts.jl` + 11 `spike/test/test_p12_*.jl` — no simulator, npe, data, `spike/p12/`, or `run_p12_*` file |

**Corpus audit — COUNTED this session, not inherited** (`awk -F, ... | sort | uniq -c` on
`corpus/manifest.csv`, plus `ls -la corpus/data`):

| Counted | Observed |
|---|---|
| `corpus/data/` | empty |
| data rows | 32 |
| `bytes = 0` | 32 of 32 |
| `physical-primary` | 2 — both `sha256 = PENDING-FETCH`, both `split = sealed_holdout` |
| `simulated-secondary` | 30 — empty `sha256`, 16 `eval` + 14 `dev` |
| fetched | 0 |

**Task 1 `<verify>` (run):** `test -f <file> && grep -q "No Phase-12 result exists" <file>` → **PASS**.

**Task 1 acceptance greps (run):** all eight required whole-file literals present
(`FROZEN SPECIFICATION` 1, `No Phase-12 result exists` 2, `sealed_holdout` 3, `PENDING-FETCH` 2,
`P12_STAGE2_LOGSCORE_MIN` 1, `no spatial borrowing, full nuisance and global borrowing` 1,
`coverage-only win` 2, `reported, not gated` 1); the 40-char sha present and equal to
`git rev-parse HEAD` at write time; exactly five numbered `##` amendment sections plus §0; all six
§5-scoped literals present (`region_rho_sample` 4, `region_rho_control` 4, `region_delta_rho` 3,
`zero forward passes` 2, `two specimens` 3, `no exchangeable sample/control pair` 1) with
`exchangeable` present (3) and `matched pair` absent (0 in §5, 0 in the whole file).

**Task 2 `<verify>` (run):** `[ "$(grep -o 'SPAT-0[1-9]' .planning/REQUIREMENTS.md | sort -u | wc -l)" -eq 9 ] && grep -q "12-SC3-AMENDMENT" .planning/ROADMAP.md` → **PASS**.
Each of SPAT-01…SPAT-09 verified present in **both** the requirement table and the traceability table
(1/1 each). `registration-aware θ` survives in ROADMAP only inside the struck/VOID annotation
(line 255). `git diff HEAD -- .planning/ROADMAP.md` touched 3 lines, all inside the Phase-12 block.

**Post-commit:** `git diff --diff-filter=D --name-only HEAD~2 HEAD` → empty (no deletions).
Files changed by the two commits: exactly the three in `files_modified`. `.planning/STATE.md`,
`src/` and `spike/` untouched. `git log --oneline -- 12-SC3-AMENDMENT.md` → **one** commit, as the
freeze requires.

**Constants quoted in the amendment, verified against `spike/validation/p12_consts.jl`:**
`P12_STAGE2_N_MIN = 271` (:406), `P12_STAGE2_COVERAGE_TOST_DELTA = 0.03` (:405),
`P12_STAGE2_LOGSCORE_MIN = 0.02` (:434), `P12_COVERAGE_NOMINAL = 0.90` (:404),
`P12_ABLATION_R1 = 0.0` (:313), `P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150` (:476). The Wald-sizing /
Wilson-reporting distinction in §4(h) matches `p12_consts.jl:409-422` and `p11_stats.jl:49-64` (whose
`wilson_ci` docstring does carry "WALD IS DELIBERATELY NOT USED") — the plan's citation was correct
and is reproduced unchanged.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 — Bug] Three stale line citations in the plan's context prose, corrected before freezing**

The plan is a write-once document's source, so a citation frozen with the wrong line number is
permanent. In all three cases the **cited content exists and was verified verbatim**; only the line
numbers had drifted, so the numbers were corrected rather than the substance changed. Recorded rather
than silently improved, per the wave-1 instruction.

| Plan cited | Actual | Verified |
|---|---|---|
| `12-01:172-178` (dataset as the independent unit) | `12-01-PLAN.md:198-204` | both quoted strings exact at 198-204: *"the 64 regions of one image share every nuisance and the global term"*, *"by up to 8× in sd terms"*. 172-178 is the `THE r₁ PRIOR ITSELF` block. |
| `infer.jl:107` / `infer.jl:108` | `src/amortized/infer.jl:107-121` (docstring + method) and `:117-121` (the formula) | `:108` is the signature line inside the docstring; `mean(ρs .- ρc)` is at `:120`. 13-SC2-AMENDMENT sets the precedent by citing `infer.jl:117-121`. |
| `src/results.jl:181` (`region_delta_rho`) | `src/results.jl:180-191` (the sketch block) | `:181` is the `struct SpatialColocResult` line; `region_delta_rho` is at `:183`. The plan itself uses `180-191` in §3. |

- **Found during:** Task 1, pre-write citation verification
- **Files modified:** `12-SC3-AMENDMENT.md` (the frozen text carries the corrected ranges)
- **Commit:** `ac9afdb`

**2. [Rule 1 — Bug] §0's "Phase-12 executable state" row understated what is on disk**

The first draft wrote *"`p12_consts.jl` and its literal-assertion test only"*, mirroring
13-SC2-AMENDMENT. `git ls-files | grep -i p12` shows **11** `spike/test/test_p12_*.jl` scaffolding
files, not one. A §0 row is the one part of this document whose whole value is being literally true,
so it was corrected before commit to name the 11 test files and to state positively what is
**absent** (no simulator, npe, data, `spike/p12/` or `run_p12_*` file). The substantive claim —
no result, no pool, no trained net — is unchanged and stronger for being specific.

- **Found during:** Task 1, pre-commit re-verification
- **Commit:** `ac9afdb`

**3. [Rule 3 — Blocking] A markdown line-wrap broke a required literal**

`**no\nexchangeable sample/control pair**` split the acceptance-grep literal
`no exchangeable sample/control pair` across a newline and a `**` marker; the grep returned 0. Fixed
by moving the bold boundary and rewrapping. Caught by running the acceptance criterion rather than
by reading — which is the point of it being greppable.

- **Found during:** Task 1, acceptance verification (pre-commit)
- **Commit:** `ac9afdb`

### Assumption Drift (advisory)

None material. The plan's §5 source text, the settled decisions in the dispatch brief, and
`12-CONSTANTS-FOR-CONFIRMATION.md` agreed on every point checked (Δρ semantics, the structural
identity, single-stack pool, `n_low` as a Tier-2 append, Guard 4 descriptive, the real arm reported
not gated). No contradiction between clauses was found that required blocking.

### One thing checked and deliberately NOT changed

§4(a) states SC3 requires no ground truth, while §5(3) describes a per-region Δρ\* ground truth from
`harness.jl:124-136`. These read as contradictory on a fast pass and are not: Δρ\* is a **simulator**
quantity consumed only by 12-16's parameter-coverage cross-check and 12-18's SBC, while the SC3
predictive-coverage criterion is ground-truth-free and single-stack on both arms. Rather than
soften either clause, an explicit scoping sentence was added to §5(3) naming the boundary, so a
future reader cannot file the pair as a self-contradiction.

## Concurrency Notes

Two other agents were live on `gsd/v2.0-milestone` during execution (a Phase-13 executor and the
12-01 executor finishing a `spike/test/` wiring fix). No `index.lock` contention or race was
encountered; both commits used explicit pathspecs (`git commit -- <paths>`), and `HEAD` was
`af8e1ca9` at both the write and the commit of Task 1, so the §0 sha is the true parent of the
amendment's own commit. `spike/test/test_p12_consts.jl` was modified in the working tree by the other
agent throughout and was never staged here.

## Known Stubs

None. This plan produces documents only; it runs no simulation and trains nothing.

## Self-Check: PASSED

- `FOUND: .planning/phases/12-spatial-colocalization-map/12-SC3-AMENDMENT.md`
- `FOUND: .planning/REQUIREMENTS.md` (SPAT-01…SPAT-09, 9 unique)
- `FOUND: .planning/ROADMAP.md` (`**AMENDED**` + `12-SC3-AMENDMENT.md` in the Phase-12 block)
- `FOUND: ac9afdb` — `docs(12-02): freeze 12-SC3-AMENDMENT.md before any Phase-12 result exists`
- `FOUND: 092e12e` — `docs(12-02): mint SPAT-01..SPAT-09 and annotate the Phase-12 ROADMAP block`
