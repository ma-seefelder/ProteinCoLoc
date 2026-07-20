---
phase: 08-external-physical-ground-truth-corpus
plan: 01
subsystem: corpus
tags: [corpus, data-contract, offline-gate, fixture, anti-snooping, decoupling]
requires: []
provides:
  - corpus/config.jl pre-declared consts (CORPUS_MASTER_SEED, MANIFEST_SCHEMA_VERSION, TIERS, SPLITS, ROLES, MANIFEST_COLUMNS, CBS_SPLIT_SALT, DOWNLOAD_TIMEOUT_S, MIN_PATCH_VALID_PIXELS, CORPUS_DATA_DIR)
  - corpus/test/runtests.jl single offline test gate (downstream include stubs wired)
  - corpus/test/make_fixture.jl deterministic two-channel TIFF fixture generator
  - committed corpus/test/fixtures/ (fixture_ch1.tif, fixture_ch2.tif)
  - scoped .gitignore keeping corpus/data bytes out of git
affects:
  - .gitignore
tech-stack:
  added: []
  patterns:
    - "pre-declared SCREAMING_SNAKE const block (analog spike/comparator/config.jl)"
    - "Random123 Philox keyed determinism (analog src/amortized/datagen.jl)"
    - "single-runner include-wiring test gate (analog spike/test/runtests.jl)"
    - "scoped-negation .gitignore (analog .gitignore spike/data/cache)"
    - "read-only include() reach into frozen src/ summary path"
key-files:
  created:
    - corpus/config.jl
    - corpus/test/runtests.jl
    - corpus/test/make_fixture.jl
    - corpus/test/fixtures/.gitkeep
    - corpus/test/fixtures/fixture_ch1.tif
    - corpus/test/fixtures/fixture_ch2.tif
  modified:
    - .gitignore
decisions:
  - "CORPUS_MASTER_SEED = 0x0000000000C05EED — distinct from every other repo master seed (spike VAL 0x5BC0FFEE / NPE 0xC0FFEE / comparator 0x00C0FFEE) so corpus draw streams cannot collide"
  - "runtests.jl includes src/LoadImages.jl + src/colocalization.jl read-only rather than `using ProteinCoLoc`, avoiding the heavy display-bound GLMakie load while proving the unchanged summary path"
  - "fixture conversion test regenerates the committed fixture if a fresh checkout lacks it (deterministic), so the gate is green from a clean clone"
metrics:
  duration: ~7 min
  tasks: 2
  files: 7
  completed: 2026-07-20
---

# Phase 8 Plan 01: Corpus Wave-0 Foundation Summary

Establishes the decoupled `corpus/` module's Wave-0 gate: pre-declared data-contract consts, a single OFFLINE test runner, a deterministic synthetic two-channel TIFF fixture, and a scoped `.gitignore` that keeps downloaded image bytes out of git — proving the fixture converts through the UNCHANGED `src/` patch/correlation summary path to a finite 8x8 grid with zero network and zero `src/` edits.

## What Was Built

**Task 1 — Pre-declared config + scoped .gitignore (commit 607cc0a)**
- `corpus/config.jl`: a pure declarative `const` block mirroring the `spike/comparator/config.jl` pre-registration idiom. Declares `CORPUS_MASTER_SEED` (0x0000000000C05EED, distinct namespace), `MANIFEST_SCHEMA_VERSION = 1`, the `TIERS`/`SPLITS`/`ROLES` vocabulary, the 14-field `MANIFEST_COLUMNS` D-08 schema, `CBS_SPLIT_SALT`, `DOWNLOAD_TIMEOUT_S`, `MIN_PATCH_VALID_PIXELS` (documents the `_exclude_zero` >15-px floor), and `CORPUS_DATA_DIR`. Includes a `# NOTE` mandating cryptographic SHA-256 (never Base `hash`) for download integrity, and that mismatch = hard error vs unreachable = skip-with-flag.
- `.gitignore`: scoped-negation block ignoring `corpus/data/*` (downloaded bytes, D-04) while keeping the manifest, tests, and `corpus/test/fixtures/` tracked.

**Task 2 — Offline test gate + deterministic fixture (commit 0373d29)**
- `corpus/test/make_fixture.jl`: `make_fixture(dir; seed=CORPUS_MASTER_SEED)` deterministically writes two 64x64 `Float64` TIFFs via a Random123 Philox stream. Channels share a period-16 smooth positive field (varies within every 8x8 patch ⇒ non-zero variance ⇒ finite Pearson) plus per-channel jitter, all intensities strictly in (0.05, 0.95) so every patch keeps all 64 pixels above the >15 floor.
- `corpus/test/runtests.jl`: the single offline gate. Imports the Statistics/StatsBase/Images names the frozen `src/` files resolve at call time, includes `config.jl` + `src/LoadImages.jl` + `src/colocalization.jl` read-only, then runs two testsets — round-trip + byte-identical determinism, and fixture conversion (build `MultiChannelImage`, `patch(·,8)` + `correlation`, assert grid is not all-`missing` and mean is finite). Commented include stubs reserve where 08-02…08-07 wire their suites so ONE command stays the gate.
- `corpus/test/fixtures/`: committed `fixture_ch1.tif` / `fixture_ch2.tif` + `.gitkeep`.

## Verification Results

- `julia --project=. corpus/test/runtests.jl` → 12/12 pass in ~14s, zero network.
- `julia -e 'include("corpus/config.jl"); ...'` → distinct-seed, tier/split, and 14-column assertions all hold (`config OK`).
- `git diff --name-only <base> HEAD -- src/` → empty (src/ provably untouched, hard CLAUDE.md constraint).
- `git status --porcelain corpus/data` → empty; `git check-ignore corpus/data/foo.tif` confirms the bytes dir is ignored while `corpus/test/fixtures/*.tif` are NOT ignored.

## Deviations from Plan

None - plan executed exactly as written. (The plan's specified `CBS_SPLIT_SALT = 0xD1B54A32D192ED03` shares its hex with datagen's `FOLD_SALT`; this is not a collision in practice because `CORPUS_MASTER_SEED` is a distinct key namespace from the datagen master seed. Documented inline in config.jl; used verbatim as the plan specified.)

## Known Stubs

The commented `include(...)` lines in `runtests.jl` are intentional wiring stubs for downstream Phase-8 plans (08-02 hash, 08-03 manifest, 08-04 fetch, 08-05 load, 08-06 CBS, 08-07 anchors). They carry no runtime behavior yet and are the designed extension point for the single-gate pattern — not unresolved data stubs.

## Self-Check: PASSED

- FOUND: corpus/config.jl
- FOUND: corpus/test/runtests.jl
- FOUND: corpus/test/make_fixture.jl
- FOUND: corpus/test/fixtures/.gitkeep
- FOUND: corpus/test/fixtures/fixture_ch1.tif
- FOUND: corpus/test/fixtures/fixture_ch2.tif
- FOUND commit: 607cc0a (Task 1)
- FOUND commit: 0373d29 (Task 2)
