---
phase: 08-external-physical-ground-truth-corpus
plan: 02
subsystem: corpus
tags: [corpus, integrity, sha256, fetch, atomic-write, skip-with-flag, offline-gate, decoupling]
requires:
  - corpus/config.jl consts (DOWNLOAD_TIMEOUT_S, MANIFEST_SCHEMA_VERSION) from 08-01
  - corpus/test/runtests.jl single offline gate from 08-01
  - corpus/test/fixtures/fixture_ch1.tif committed fixture from 08-01
provides:
  - corpus/hash.jl (file_sha256 streaming, manifest_hash fixed-order fold, canonical_config, subhashes)
  - corpus/fetch.jl (fetch_verified with two distinct D-05 paths + atomic commit + bounded timeout)
  - corpus/test/test_hash.jl + corpus/test/test_fetch.jl wired into the offline gate
affects:
  - corpus/test/runtests.jl
tech-stack:
  added: []
  patterns:
    - "SHA-256 content-hash fold over fixed-order NUL-separated bytes (analog spike/comparator/table.jl:comparator_hash)"
    - "canonical sort-by-key config serializer (analog table.jl:_cmp_canonical)"
    - "atomic .tmp -> verify -> mv(force=true) durable write (analog datagen.jl:write_shard)"
    - "skip-with-flag NamedTuple on unreachable source (analog tapqir_bridge.jl:tapqir_anchor)"
    - "hard-error-on-mismatch, distinct from the skip path (analog datagen.jl:open_or_invalidate)"
    - "overridable network seam via a Ref{Function} + invokelatest (world-age-safe test injection)"
    - "guarded isdefined(...) || include(...) sibling loading (analog table.jl:47-50)"
key-files:
  created:
    - corpus/hash.jl
    - corpus/fetch.jl
    - corpus/test/test_hash.jl
    - corpus/test/test_fetch.jl
  modified:
    - corpus/test/runtests.jl
decisions:
  - "Network I/O isolated in a Ref{Function} `_DOWNLOAD_HOOK` seam invoked through `Base.invokelatest`, so tests swap in a local-copy stub by ASSIGNING the ref (runtime mutation, not a method redefinition) — sidesteps Julia world-age when the stub closure is defined in the same top-level testset block that calls fetch_verified. This is the offline-testability mechanism (no live network)."
  - "file_sha256 streams the file in 1 MiB blocks through SHA256_CTX rather than read()-ing whole, so a >100 MB CBS archive is never fully resident in memory; the digest is identical to bytes2hex(SHA.sha256(read(path)))."
  - "subhashes re-expressed over cryptographic SHA-256 (file_sha256), NOT datagen's Base-hash variant — corpus download integrity (D-04/D-05) requires a cryptographic digest (08-PATTERNS L91-96 caveat)."
  - "manifest_hash documents an explicit FIXED-ORDER contract: the fold is order-sensitive by design so identical contract inputs content-address identically and any source/config edit is detectable."
metrics:
  duration: ~8 min
  tasks: 2
  files: 5
  completed: 2026-07-20
---

# Phase 8 Plan 02: Corpus Integrity + Fetch Spine Summary

Builds the D-04/D-05 integrity spine of the `corpus/` module: a stdlib-`SHA` content-hash
utility (`corpus/hash.jl`) and a `Downloads`-based `fetch_verified` (`corpus/fetch.jl`) whose
two D-05 code paths never conflate — an UNREACHABLE source degrades to a green skip-with-flag,
while a content-hash MISMATCH is a HARD `error(...)` abort. Every artifact is committed
atomically (`.tmp -> verify -> mv`), the download is bounded by `DOWNLOAD_TIMEOUT_S`, and the
whole layer is buildable and testable fully OFFLINE against the committed fixture.

## What Was Built

**Task 1 — SHA-256 content-hash utility (RED 95f644a, GREEN 494afe5)**
- `corpus/hash.jl` (guarded `isdefined || include("config.jl")`, `using SHA` stdlib):
  - `file_sha256(path)` — streams 1 MiB blocks through `SHA256_CTX` (large archives never fully
    in memory); returns 64-char lowercase hex equal to `bytes2hex(SHA.sha256(read(path)))`.
  - `manifest_hash(src_files, canonical_config)` — folds `SHA.update!` over fixed-order,
    NUL-separated file bytes then the canonical config string (the `comparator_hash` idiom);
    order-sensitive by an explicit documented contract.
  - `canonical_config(nt)` — sort-by-key `"k=repr(v)"` serializer (the `_cmp_canonical` idiom),
    stable across insertion order.
  - `subhashes(src_files)` — `Dict path=>file_sha256(f)` for which-source-changed diagnosis,
    over cryptographic SHA-256 (not Base hash).
- `corpus/test/test_hash.jl` — determinism, same-bytes-equal, flipped-byte-differs, canonical
  key-order stability, `manifest_hash` order/config sensitivity, and 64-char subhash assertions.

**Task 2 — fetch_verified with two distinct D-05 paths (RED 6627901, GREEN 957a32e)**
- `corpus/fetch.jl` (guarded includes of config.jl + hash.jl, `using Downloads` stdlib):
  - `_DOWNLOAD_HOOK::Ref{Function}` — the ONE network seam; default closure calls
    `Downloads.download(url, dest; timeout = DOWNLOAD_TIMEOUT_S)`. Invoked via
    `Base.invokelatest` so a test stub is world-age-safe.
  - `fetch_verified(url, dest, expected_sha256::Union{AbstractString,Nothing})`:
    Path A (skip) — download in the ONLY try/catch; unreachable ⇒ `(status=:skipped,
    value=missing, reason="unreachable: ...")`, no `.tmp` left, `InterruptException` rethrown.
    Bootstrap — `nothing`/`""` expected ⇒ record hash without comparing, commit atomically,
    `(status=:present, sha256=got, bootstrapped=true)`.
    Path B (abort, OUTSIDE any try/catch) — digest `!=` ⇒ `rm(.tmp)` + `error(...)`.
    Success ⇒ atomic `mv(.tmp -> dest)`, `(status=:present, sha256=got)`.
- `corpus/test/test_fetch.jl` — OFFLINE, network seam stubbed with a local fixture copy:
  `:present` atomic commit (no `.tmp`), uppercase-hash case-insensitivity, bootstrap, wrong-hash
  `@test_throws ErrorException` (`.tmp`+dest removed), and unreachable `:skipped` (never throws).

Both testsets wired into `corpus/test/runtests.jl` — the single offline gate remains one command.

## Verification Results

- `julia --project=. corpus/test/runtests.jl` → **49/49 pass** in ~7s, zero network
  (fixture 12 + hashing 17 + fetch 20).
- Structural asymmetry (grep-confirmed, plan acceptance): the download `_download(url, tmp)`
  is the ONLY statement inside the try/catch (L87-89); the mismatch `error(...)` sits OUTSIDE
  it (L108); `Downloads.download(...; timeout = DOWNLOAD_TIMEOUT_S)` present (L51).
- `git diff --name-only <base> HEAD -- src/` → empty (src/ provably untouched, hard CLAUDE.md
  constraint).
- No untracked/generated files; no unexpected deletions across either feat commit.

## Deviations from Plan

None — plan executed exactly as written.

Two implementation choices the plan left open, resolved within its stated latitude (not
deviations): (1) the plan offered "`file://` URL OR monkeypatch the download seam" for offline
tests — chose the `_DOWNLOAD_HOOK` seam (the plan's own fallback for Windows `file://`
unavailability), hardened with `Base.invokelatest` for world-age safety; (2) `manifest_hash`'s
`canonical_config` parameter typed `AbstractString` (accepts `SubString`/`String`) rather than
the plan's literal `String`, widening acceptance with no behavioral change.

## TDD Gate Compliance

Both tasks followed RED -> GREEN. Git log shows, per task, a `test(08-02)` commit (failing test,
gate red on missing implementation) immediately followed by a `feat(08-02)` commit (implementation,
gate green). No REFACTOR commit needed. Gate sequence intact:
`test 95f644a -> feat 494afe5` (hash), `test 6627901 -> feat 957a32e` (fetch).

## Self-Check: PASSED

- FOUND: corpus/hash.jl
- FOUND: corpus/fetch.jl
- FOUND: corpus/test/test_hash.jl
- FOUND: corpus/test/test_fetch.jl
- FOUND commit: 95f644a (Task 1 RED)
- FOUND commit: 494afe5 (Task 1 GREEN)
- FOUND commit: 6627901 (Task 2 RED)
- FOUND commit: 957a32e (Task 2 GREEN)
