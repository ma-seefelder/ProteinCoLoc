---
phase: 08-external-physical-ground-truth-corpus
plan: 03
subsystem: corpus
tags: [corpus, data-contract, manifest, tier-guard, sealed-holdout, content-addressing, offline-gate, decoupling]
requires:
  - corpus/config.jl consts (MANIFEST_COLUMNS, TIERS, SPLITS, ROLES, MANIFEST_SCHEMA_VERSION, CORPUS_MASTER_SEED, MIN_PATCH_VALID_PIXELS) from 08-01
  - corpus/hash.jl (manifest_hash, canonical_config) from 08-02
  - corpus/test/runtests.jl single offline gate from 08-01
  - corpus/test/fixtures/fixture_ch1.tif + fixture_ch2.tif committed fixture from 08-01
  - src/LoadImages.jl (MultiChannelImage ctor, load_tiff) — read-only reuse
  - src/colocalization.jl (patch, correlation) — read-only reuse
provides:
  - corpus/manifest.jl (build_manifest/validate_manifest schema; physical_anchors/cbs_benchmark typed accessors; corpus_default split guard; open_sealed_holdout; mixed_tier; validate_anchor_row D-03 predicate; write_manifest content-addressed atomic writer)
  - corpus/load.jl (load_anchor two-TIFF -> MultiChannelImage; summary_grid; has_signal Pitfall-6 guard)
  - corpus/test/test_manifest.jl + corpus/test/test_load.jl wired into the offline gate
affects:
  - corpus/test/runtests.jl
tech-stack:
  added: []
  patterns:
    - "tidy DataFrame(col=...) assembly (analog spike/comparator/table.jl:build_table)"
    - "content-addressed atomic .tmp -> verify -> mv(force=true) CSV+JLD2 writer (analog table.jl:write_table)"
    - "structural tier/split separation — no default both-tier or sealed accessor exists (analog datagen.jl no-standardize_all discipline)"
    - "guarded isdefined(...) || include(...) read-only reach into frozen src/ + sibling corpus files"
    - "pre-registration CSV comment header quoting schema version + master seed"
key-files:
  created:
    - corpus/manifest.jl
    - corpus/load.jl
    - corpus/test/test_manifest.jl
    - corpus/test/test_load.jl
  modified:
    - corpus/test/runtests.jl
decisions:
  - "write_manifest content-addresses on MANIFEST_SRC_FILES (config.jl+hash.jl+manifest.jl bytes) AND a SHA-256 digest of the manifest CSV content, folded into the canonical_config — so BOTH a contract-code edit and a manifest-data edit resolve to a new dir (D-08 detectability over code + data)"
  - "physical_anchors filters physical-primary AND split!=sealed_holdout, so on a valid manifest it returns ZERO rows by construction — the anchors are unreachable via any default accessor; open_sealed_holdout(;reason) is the sole path (D-09 anti-snooping)"
  - "CSV.write repositions an IOStream and clobbers a pre-written comment header; write_manifest serializes the table to a dedicated IOBuffer first, then writes header + table bytes to the file"
metrics:
  duration: ~15 min
  tasks: 2
  files: 5
  completed: 2026-07-20
---

# Phase 8 Plan 03: Corpus Data-Contract + Conversion Layer Summary

Builds the versioned validation data contract (SC3) and the structural guards that make
never-conflated tiers (SC2/D-06) and the sealed physical-anchor holdout (D-09) true BY
CONSTRUCTION — mirroring the Phase-3 "no `standardize_all` symbol exists" leakage-impossible
discipline. `corpus/manifest.jl` enforces the 14-column D-08 schema, exposes SEPARATE typed
tier accessors (no default both-tier or sealed accessor exists), enforces the D-03 anchor
acceptance predicate, and content-addresses the committed contract; `corpus/load.jl` converts a
two-single-channel-TIFF anchor into a `MultiChannelImage` through the UNCHANGED frozen `src/`
summary path. Fully offline against the committed fixture, zero `src/` edits.

## What Was Built

**Task 1 — Manifest schema + tier/split guards + D-03 predicate + content-addressed writer (RED 6bde9c0, GREEN 6d7cd14)**
- `corpus/manifest.jl` (guarded includes of config.jl + hash.jl; `using DataFrames/CSV/JLD2`):
  - `build_manifest(rows)` — tidy assembly of NamedTuple/DataFrameRow/Dict rows into a DataFrame
    with EXACTLY the 14 `MANIFEST_COLUMNS`, in declared order.
  - `validate_manifest(df)` — throws on a missing column, a tier/split/role outside the pre-declared
    `config.jl` vocabulary, or the D-09 invariant violation (a `physical-primary` row whose split
    != `sealed_holdout`); returns `df` on success.
  - D-06/SC2 structural separation: `physical_anchors(df)` (physical-primary AND not sealed ⇒ ZERO
    rows on a valid manifest), `cbs_benchmark(df)` (simulated-secondary only), `corpus_default(df)`
    (asserts it yields no sealed_holdout row), and `mixed_tier(df)` — the SOLE both-tier path.
  - D-09 split guard: `open_sealed_holdout(df; reason)` `@assert !isempty(reason)` and returns the
    sealed physical anchors — documented as Phase-16-only, the sole reachable path.
  - D-03 predicate: `validate_anchor_row(row)` returns false unless truth_label ∈ {coloc,segregated}
    AND license + accession are non-empty AND coloc_ground_truth is an ASSERTED 1.0/0.0 (a
    Pearson/Manders-derived score is rejected — Pitfall 1).
  - D-08 writer: `write_manifest(df; outdir)` validates, then content-addresses into
    `joinpath(outdir, manifest_hash(MANIFEST_SRC_FILES, canonical_config(config)))` (config folds
    schema version + master seed + a SHA-256 digest of the manifest content), writing `manifest.csv`
    (pre-registration comment header quoting `MANIFEST_SCHEMA_VERSION` + `CORPUS_MASTER_SEED`) and a
    `manifest.jld2` twin, each via `.tmp` → integrity-check → `mv(force=true)`.
- `corpus/test/test_manifest.jl` — asserts every behavior bullet: exact schema, validate passes;
  bad tier/split/role each throw; physical+dev throws (D-09); missing column throws; corpus_default
  yields zero sealed rows; physical_anchors empty; cbs_benchmark only simulated; open_sealed_holdout
  reason guard; mixed_tier both tiers; D-03 accept/reject cases; write_manifest header + JLD2 twin +
  deterministic dir.

**Task 2 — Two-channel TIFF → MultiChannelImage loader via read-only src reuse (a97be92)**
- `corpus/load.jl` (guarded READ-ONLY includes of `src/LoadImages.jl` + `src/colocalization.jl`):
  - `load_anchor(name, ch1_path, ch2_path; channels)` — calls the FROZEN convenience
    `MultiChannelImage` constructor (auto pixel_size + Otsu threshold).
  - `summary_grid(img; num_patches=8)` — `patch(·,8)` per channel + `correlation(·,·)` through the
    UNCHANGED summary path.
  - `has_signal(grid)` — Pitfall-6 guard: true iff any non-missing entry AND finite mean.
  - Documents that the RGB/composite single-TIFF (CBS) layout is handled separately in 08-04.
- `corpus/test/test_load.jl` — builds a MultiChannelImage from the committed fixture, asserts the
  8×8 grid is non-all-missing with a finite mean, has_signal true, and false on an all-missing grid.

Both testsets wired into `corpus/test/runtests.jl` — the single offline gate stays one command.

## Verification Results

- `julia --project=. corpus/test/runtests.jl` → **92/92 pass** in ~17s, zero network
  (fixture 12 + hashing 17 + fetch 20 + manifest 34 + load 9).
- `git diff --name-only <base> HEAD -- src/` → empty; `git status --porcelain src/` → empty
  (src/ provably untouched, hard CLAUDE.md constraint; load.jl reuses src/ read-only via include).
- No untracked/generated files left; `corpus/data/*` stays gitignored.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] CSV.write clobbered the pre-registration comment header**
- **Found during:** Task 1 (GREEN run — 3 header assertions failed).
- **Issue:** `CSV.write(io, df)` on the installed CSV.jl repositions an open `IOStream`, overwriting
  the comment header lines written before it (the file emerged with the raw table and header
  fragments), so the D-08 pre-registration header was absent from `manifest.csv`.
- **Fix:** Serialize the table to a dedicated `IOBuffer` first (`CSV.write(tbuf, df)` → `take!`),
  then write the comment header + table bytes to the file. Header now precedes the table intact.
- **Files modified:** corpus/manifest.jl
- **Commit:** 6d7cd14

Two implementation choices the plan left open, resolved within its stated latitude (not deviations):
(1) `write_manifest`'s content-address config folds a SHA-256 digest of the manifest CSV content in
addition to the source-file bytes, so D-08 detectability covers both contract-code and manifest-data
edits; (2) `build_manifest`/`validate_anchor_row` accept NamedTuple, DataFrameRow, or Dict rows via a
small `_rowget` helper, widening acceptance with no behavioral change.

## Note on scope

Per the plan objective, this plan does NOT write the committed `corpus/manifest.csv` — CBS rows in
08-04 write it first. `write_manifest` is exercised here only against temp dirs in the test.

## Known Stubs

None. The remaining commented `include(...)` lines in `runtests.jl` (test_cbs.jl @ 08-06,
test_anchors.jl @ 08-07) are intentional wiring stubs for downstream Phase-8 plans, carrying no
runtime behavior — the designed extension point for the single-gate pattern, not data stubs.

## Threat Flags

None. All new surface maps to the plan's `<threat_model>` mitigations (T-08-07 sealed-holdout guard,
T-08-08 tier separation, T-08-09 content-addressed manifest, T-08-10 anchor-row predicate).

## Self-Check: PASSED

- FOUND: corpus/manifest.jl
- FOUND: corpus/load.jl
- FOUND: corpus/test/test_manifest.jl
- FOUND: corpus/test/test_load.jl
- FOUND commit: 6bde9c0 (Task 1 RED)
- FOUND commit: 6d7cd14 (Task 1 GREEN)
- FOUND commit: a97be92 (Task 2)
