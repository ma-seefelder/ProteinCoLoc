---
phase: 08-external-physical-ground-truth-corpus
plan: 04
subsystem: corpus
tags: [corpus, cbs, simulated-secondary, tier-separation, seeded-split, zip-slip, rgb-split, manifest-contract, offline-gate, decoupling]
requires:
  - corpus/config.jl consts (CORPUS_MASTER_SEED, CBS_SPLIT_SALT, CORPUS_DATA_DIR, TIERS, MANIFEST_COLUMNS) from 08-01
  - corpus/manifest.jl (build_manifest, validate_manifest, _manifest_header) from 08-03
  - corpus/fetch.jl (fetch_verified skip-with-flag / bootstrap) from 08-02
  - corpus/test/runtests.jl single offline gate from 08-01
provides:
  - corpus/cbs.jl (cbs_rows full-corpus enumeration; assign_cbs_split seeded dev/eval; safe_extract zip-slip guard; _split_rgb/cbs_rgb_to_channels; committed_manifest/write_committed_manifest; fetch_cbs_smoke)
  - corpus/manifest.csv (committed versioned contract — 30 CBS simulated-secondary rows + 2 sealed physical-anchor placeholders)
  - corpus/test/test_cbs.jl wired into the offline gate
affects:
  - corpus/test/runtests.jl
requirements: [SC2, SC3]
tech-stack:
  added: []
  patterns:
    - "hard-coded tier=simulated-secondary/role=benchmark per row so validate_manifest enforces SC2 structurally (never conflated with anchors)"
    - "Philox4x keyed by (master_seed xor SALT, row_index) for an order-independent bit-reproducible dev/eval split (analog datagen.jl holdout/fold RNGs)"
    - "path-containment check (normpath+abspath, path-separator-aware prefix) rejecting zip-slip BEFORE any write"
    - "injectable entry-provider seam (entries= keyword) keeping safe_extract security logic testable offline with NO zip dependency"
    - "disk-I/O-free channel-selection core (_split_rgb) so RGB->two-channel logic is unit-testable without saving a TIFF"
    - "flat committed manifest.csv via serialize-to-buffer-then-header write (avoids CSV.write IOStream header clobber, analog write_manifest)"
    - "bootstrap fetch_verified(url, dest, nothing) smoke returning :skipped offline (analog tapqir_bridge :skipped contract)"
key-files:
  created:
    - corpus/cbs.jl
    - corpus/manifest.csv
    - corpus/test/test_cbs.jl
  modified:
    - corpus/test/runtests.jl
decisions:
  - "safe_extract exposes an injectable `entries=` provider seam instead of hard-wiring a zip reader: the T-08-SC constraint forbids new deps and Julia has no stdlib zip reader, so the offline gate injects crafted entries directly (testing the real path-safety code) while the default `_zip_entries` errors clearly, deferring a concrete online reader to a future/Phase-16 caller. The `:present` online branch (never hit offline) is guarded, so no reader is fabricated — honoring the network constraint (no invented artifacts)."
  - "cbs_rows emits `split=\"dev\"` as a valid-enum placeholder so validate_manifest(build_manifest(cbs_rows())) passes standalone; assign_cbs_split then overwrites it with the seeded dev/eval assignment. Keying the Philox draw by the 1-based row index (not iteration order) makes the split independent of order and thread count."
  - "source_url is the VERIFIED colocalization-benchmark.com/download/ prefix ONLY; the exact per-set zip filename is left as a documented fetch-time TODO (RESEARCH A6) rather than inventing a false-precise filename — sha256/bytes stay empty/zero until the first authorized fetch records them."
  - "The two physical-anchor rows are PENDING-provenance placeholders (physical-primary + sealed_holdout, asserted 1.0/0.0 biological truth) so the committed contract exercises the schema/D-09 guards with anchor rows present; 08-05 finalizes their accession/url/hash behind a human-verify."
metrics:
  duration: ~12 min
  tasks: 2
  files: 4
  completed: 2026-07-20
---

# Phase 8 Plan 04: CBS Ingestion + Committed Manifest Contract Summary

Ingests the Colocalization Benchmark Source (CBS) — the one fully-verified external dataset —
into the corpus strictly as the `simulated-secondary` tier (SC2/D-06), never conflated with the
physical anchors. `corpus/cbs.jl` enumerates the FULL corpus (3 channel pairs × 10 labelled
degrees = 30 rows, D-07), assigns each a seeded bit-reproducible dev/eval split keyed by
`CORPUS_MASTER_SEED ⊻ CBS_SPLIT_SALT` (D-09, never `sealed_holdout`), provides a zip-slip-safe
extractor (T-08-11) and an RGB/composite→two-channel splitter, and writes the committed versioned
`corpus/manifest.csv` (SC3) = the 30 CBS rows plus two sealed physical-anchor placeholders. The
live CBS fetch is a bootstrap smoke that returns `:skipped` offline (D-05). The whole layer is
buildable and green with ZERO network, and `src/` is provably untouched.

## What Was Built

**Task 1 — CBS enumeration + seeded split + zip-slip-safe extraction (ddf2e95)**
- `corpus/cbs.jl` (guarded includes of config.jl + manifest.jl + fetch.jl; `import Images` +
  `Random123.Philox4x` + `Random`):
  - `cbs_rows()` — the FULL corpus (D-07): for each pair in `("Red-Green","Red-Blue","Green-Blue")`
    × each degree in `0:10:90`, a row with `anchor_id` `cbs-RG-000`…`cbs-GB-090`,
    `tier="simulated-secondary"`, `role="benchmark"`, `truth_label="simulated-degree"`,
    `coloc_ground_truth = degree/100` (the ARCHIVE-LABELLED value, never computed),
    `accession="CBS/<pair>/<degree>"`, `source_url` = the verified `…/download/` PREFIX,
    CBS citation + `CC-BY-NC-SA-4.0`, `format="tiff"`, `channels` = the pair (`"R,G"` …),
    `sha256=""` / `bytes=0` (filled on first authorized fetch), placeholder `split="dev"`.
  - `assign_cbs_split(rows; seed=CORPUS_MASTER_SEED)` — `Philox4x` keyed by
    `(seed ⊻ CBS_SPLIT_SALT, row_index)` → `"dev"`/`"eval"`; bit-reproducible, order-independent,
    NEVER `sealed_holdout` (D-09). Returns new rows (inputs untouched).
  - `safe_extract(zip, dir; entries=…)` + `_safe_target_path` — validates every entry's normalized
    (`normpath`+`abspath`, path-separator-aware) path stays within `dir` and `error(...)`s on any
    escape (`../evil.tif`, absolute paths) BEFORE writing a byte (T-08-11). `entries` is an
    injectable `(name, getbytes)` provider seam so the guard is testable offline with no zip dep.
  - `_split_rgb(img, pair)` / `cbs_rgb_to_channels(path, pair)` — split an RGB/composite image into
    two `Matrix{Float64}` channels via `Images.red`/`green`/`blue`; the disk-I/O-free core is unit
    tested without saving a TIFF.
- `corpus/test/test_cbs.jl` (sections a–e) — 30-row enumeration + tier/role/label asserts;
  labelled ground truth; validate_manifest pass; deterministic dev/eval-only split; zip-slip
  rejection + safe-entry admission; per-pair channel split. Wired into `runtests.jl`.

**Task 2 — Committed manifest.csv + offline-skip fetch smoke (03c0484)**
- `corpus/cbs.jl` (added): `anchor_placeholder_rows()` (2 `physical-primary`/`sealed_holdout`
  PENDING placeholders, asserted 1.0/0.0), `committed_manifest()` (32-row validated contract),
  `write_committed_manifest()` (flat `corpus/manifest.csv` with the D-08 pre-registration comment
  header, serialize-to-buffer-then-write to dodge the CSV.write header clobber, atomic `.tmp→mv`),
  and `fetch_cbs_smoke(; nsets=1)` (bootstrap `fetch_verified(prefix, dest, nothing)` →
  `:skipped` offline, never throws, never requires network).
- `corpus/manifest.csv` — the committed contract: 30 CBS `simulated-secondary` rows (seeded
  dev/eval) + 2 sealed anchor placeholders, header stamping `MANIFEST_SCHEMA_VERSION` +
  `CORPUS_MASTER_SEED`.
- `corpus/test/test_cbs.jl` (sections f–i) — committed-contract shape (32 rows, 30/2 tier split,
  anchors sealed), header stamp, on-disk `manifest.csv` present + validates, offline-skip smoke
  (`status ∈ (:skipped,:present)`), and the D-04 staging guard (`git ls-files corpus/data` empty).

## Verification Results

- `julia --project=. corpus/test/runtests.jl` → **138/138 pass** in ~37s, ZERO network
  (fixture 12 + hashing 17 + fetch 20 + manifest 34 + load 9 + **CBS 46**).
- `git ls-files corpus/data` → empty (no `.tif`/`.zip` bytes staged, D-04/Pitfall 3).
- `git check-ignore corpus/manifest.csv` → not ignored (the committed contract is trackable).
- `git diff --name-only 842aea9 HEAD -- src/` → empty; `git status --porcelain src/` → empty
  (src/ provably untouched, hard CLAUDE.md constraint; cbs.jl reaches src/ only transitively via
  manifest.jl's read-only includes).
- No unexpected deletions across either commit; no untracked generated files left.

## Deviations from Plan

None — plan executed exactly as written.

Two implementation choices the plan left open, resolved within its stated latitude (not
deviations): (1) `safe_extract` takes an injectable `entries=` provider (defaulting to a clear
"no reader configured" error) instead of hard-wiring a zip reader — the T-08-SC no-new-deps
constraint plus the absence of a Julia stdlib zip reader make the seam the honest offline-testable
choice; the online `:present` extraction branch is guarded and never fabricated (network
constraint). (2) `cbs_rows` emits a valid `split="dev"` placeholder so it validates standalone,
with `assign_cbs_split` overwriting it — keeping the two responsibilities separable as the plan
splits them.

## Known Stubs

- `corpus/manifest.csv` physical-anchor rows carry `accession`/`source_url`/`citation`/`license` =
  `"PENDING"` and `sha256=""`/`bytes=0`. This is intentional and plan-directed: the two anchors are
  finalized in **08-05** behind a `checkpoint:human-verify` once their accession/license/hash are
  confirmed. They are already correctly tier/split/role-tagged (`physical-primary`,
  `sealed_holdout`) so the schema and D-09 guards are exercised now.
- CBS rows carry `sha256=""`/`bytes=0` and a prefix-only `source_url` by design — the trusted
  baseline hash + exact zip filename + byte size are recorded on the first authorized online fetch
  (RESEARCH A6). No false-precise provenance is invented.
- `safe_extract`'s default `_zip_entries` errors until a concrete online zip reader is wired; the
  offline gate injects entries directly and the online branch is guarded. No runtime behavior is
  stubbed on the offline path.

## Threat Flags

None. All new surface maps to the plan's `<threat_model>`: T-08-11 (zip-slip) mitigated by
`safe_extract`/`_safe_target_path`; T-08-12 (tier conflation) mitigated by hard-coded
simulated-secondary/benchmark tagging + `validate_manifest`; T-08-06 (offline CI) satisfied by the
`:skipped` smoke; T-08-01 (bytes into git) asserted empty; T-08-SC (supply chain) honored — no new
deps.

## Self-Check: PASSED

- FOUND: corpus/cbs.jl
- FOUND: corpus/test/test_cbs.jl
- FOUND: corpus/manifest.csv
- FOUND commit: ddf2e95 (Task 1)
- FOUND commit: 03c0484 (Task 2)
