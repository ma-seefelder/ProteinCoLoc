# Phase 8: External Physical Ground-Truth Corpus - Pattern Map

**Mapped:** 2026-07-20
**Files analyzed:** 9 new/modified corpus-module files
**Analogs found:** 9 / 9 (all exact or role-match; every analog VERIFIED against real source)

> Scope note: this is a **decoupled data-acquisition + data-contract** phase. The recommended
> home is a new top-level `corpus/` module (or `test/corpus/`), NOT `src/` model code. It reaches
> `src/LoadImages.jl` / `src/colocalization.jl` **read-only** via `include()`. Every mechanic
> already exists in Phases 3 and 9 and should be reused near-verbatim — the real work is curation,
> not new infrastructure. Concrete file names below are the planner's to finalize; classifications
> and analogs are load-bearing.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `corpus/config.jl` | config | declarative | `spike/comparator/config.jl` | exact |
| `corpus/hash.jl` (or folded into fetch) | utility | transform | `spike/comparator/table.jl` (`comparator_hash`) + `src/amortized/datagen.jl` (`cache_hash`/`subhashes`) | exact |
| `corpus/fetch.jl` | service | file-I/O + request-response (network) | `datagen.jl:write_shard` (atomic) + `tapqir_bridge.jl` (skip-with-flag, `_wait_timeout`) | role-match (composed) |
| `corpus/manifest.jl` (schema + typed accessors + tier/split guards) | model + guard | CRUD / transform + access-control | `table.jl:build_table`/`write_table` + `datagen.jl` structural-holdout discipline | role-match |
| `corpus/manifest.csv` (the committed contract artifact) | config/data | transform (tabular) | `table.jl:write_table` (DataFrame → CSV with header) | exact |
| `corpus/load.jl` (two-channel TIFF → `MultiChannelImage`) | loader/adapter | file-I/O + transform | `src/LoadImages.jl:load_tiff` (L214) + `MultiChannelImage` ctor (L126 / convenience L432) | exact (read-only reuse) |
| `corpus/test/runtests.jl` | test | — | `spike/test/runtests.jl` | exact |
| synthetic two-channel TIFF fixture | test fixture | file-I/O | `spike/data/cache/fixture/` committed-fixture idiom | role-match |
| `.gitignore` (scoped-negation for corpus data dir) | config | — | existing `.gitignore` L404 / L406-409 scoped-negation | exact |

## Pattern Assignments

### `corpus/config.jl` (config, declarative)

**Analog:** `spike/comparator/config.jl` (VERIFIED) — pure pre-declared `const` block, no logic.

Copy the pre-registration pattern: a `CORPUS_MASTER_SEED`, a schema-version tag, the tier/split
enums, and any byte-size/threshold consts declared HERE before any manifest exists (anti-snooping).

**Master-seed + pre-declared const pattern** (`config.jl:53-56`):
```julia
# --- Harness master seed (D-12) -------------------------------------------------
# Random123 (Philox) master seed for the whole comparator harness, so the emitted
# table is bit-reproducible across runs and thread counts, exactly like the Phase-3 cache.
const MASTER_SEED = 0x00000000_00C0FFEE  # comparator harness Philox master seed
```
For the corpus, declare e.g. `const CORPUS_MASTER_SEED = 0x...` plus the enums as `const`
tuples/sets so the tier/split guards validate against a single source (mirror the header
provenance discipline at `config.jl:21-31`).

**Disjoint-salt pattern for split streams** (`datagen.jl:56-60`) — reuse when seeding a CBS
dev/eval split so its stream never collides with any other draw stream:
```julia
const HOLDOUT_SALT = 0x9E3779B97F4A7C15   # reserved holdout stream (D-10)
const FOLD_SALT    = 0xD1B54A32D192ED03   # k-fold permutation stream
```

---

### `corpus/hash.jl` (utility, transform)

**Analog:** `spike/comparator/table.jl:comparator_hash` (VERIFIED) — SHA-256 stdlib content hash;
`src/amortized/datagen.jl:cache_hash`/`subhashes` (VERIFIED) — the diagnosable per-file variant.

**Manifest content-hash pattern** (`table.jl:225-233`) — fold `SHA.update!` over fixed-order,
NUL-separated source bytes then a canonical config string:
```julia
function comparator_hash(config::NamedTuple; src_files = CMP_SRC_FILES)
    ctx = SHA.SHA256_CTX()
    for f in src_files                                     # FIXED order = deterministic
        SHA.update!(ctx, read(f))
        SHA.update!(ctx, UInt8[0x00])                      # inter-file separator
    end
    SHA.update!(ctx, Vector{UInt8}(_cmp_canonical(config)))
    return bytes2hex(SHA.digest!(ctx))
end
```

**Canonical (sort-by-field) config serializer** (`table.jl:211-214`, mirrored in `datagen.jl:221-224`):
```julia
function _cmp_canonical(config::NamedTuple)::String
    ks = sort(collect(keys(config)))
    return join(("$(k)=$(repr(getproperty(config, k)))" for k in ks), ";")
end
```

**Per-file sub-hashes for which-source-changed diagnosis** (`datagen.jl:250-252`):
```julia
function subhashes(src_files = DATAGEN_HASH_SRC_FILES)
    return Dict{String,String}(f => string(hash(read(f)); base = 16, pad = 16) for f in src_files)
end
```

> NOTE (verified nuance to flag to the planner): `datagen.jl:cache_hash` uses **Base `hash`**
> (not SHA-256) deliberately, to avoid touching a fragile Wave-0 Manifest (see its comment at
> `datagen.jl:199-205`). For the corpus's **download-integrity** hash you MUST use real SHA-256
> (`SHA.sha256` / `SHA256_CTX`, as in `table.jl`), NOT Base `hash` — D-04/D-05 integrity requires
> a cryptographic digest. Use the `table.jl` SHA-256 idiom for file+manifest hashing; borrow only
> the `subhashes` *diagnosis structure* from `datagen.jl`, re-expressed over `SHA.sha256`.

A `file_sha256(path) = bytes2hex(SHA.sha256(read(path)))` one-liner (stream in chunks for large
files) is the download-verification primitive (RESEARCH §Fetch + Integrity Mechanics).

---

### `corpus/fetch.jl` (service, file-I/O + network request-response)

**Analog (composed from three VERIFIED sources):**
- Atomic `.tmp → integrity-check → mv(force=true)`: `datagen.jl:write_shard` + `table.jl:write_table`.
- Skip-with-flag NamedTuple on failure: `tapqir_bridge.jl:tapqir_anchor`.
- Subprocess/stall timeout: `tapqir_bridge.jl:_wait_timeout`.

**Atomic durable-write idiom** (`datagen.jl:323-332`) — the exact `.tmp → reopen-integrity → mv` shape:
```julia
function write_shard(dir, n::Integer, theta, summary_min, global_index, imsize)
    path = shard_path(dir, n)
    tmp  = path * ".tmp"
    JLD2.jldsave(tmp; theta, summary_min, global_index, imsize, schema_version = SCHEMA_VERSION)
    JLD2.jldopen(tmp, "r") do f                            # integrity check before commit
        @assert haskey(f, "theta") "integrity check failed: $tmp missing theta"
    end
    mv(tmp, path; force = true)                            # atomic commit
    return path
end
```
The CSV-writer variant (`table.jl:277-284`) shows the same idiom for a text artifact with a
non-empty integrity check (`@assert filesize(csv_tmp) > 0`).

**Skip-with-flag contract (D-05 UNREACHABLE → green)** (`tapqir_bridge.jl:180-205`) — return a
`(status=..., value/reason=...)` NamedTuple, NEVER throw, on an environment/network failure:
```julia
function tapqir_anchor(; tol = TAPQIR_TOL)
    avail = try
        _tapqir_env_available()
    catch
        false
    end
    avail || return (status = :skipped, value = missing,
                     reason = "Tapqir/Python env unavailable")
    try
        # ... work ...
    catch e
        return (status = :skipped, value = missing, reason = sprint(showerror, e))
    end
end
```
Apply as: wrap `Downloads.download` in `try/catch` → `(status=:skipped, reason="unreachable: ...")`.

**Stall-proof timeout** (`tapqir_bridge.jl:61-72`) — reuse if a download can hang (ASVS DoS
mitigation, RESEARCH §Security):
```julia
function _wait_timeout(proc::Base.Process, timeout_s::Real)
    t0 = time()
    while process_running(proc) && (time() - t0) < timeout_s
        sleep(0.25)
    end
    if process_running(proc)
        try; kill(proc); catch; end
        sleep(0.5)
        try; kill(proc, Base.SIGKILL); catch; end
    end
    return proc
end
```

> ASYMMETRY GUARD (D-05, Pitfall 4): the skip-with-flag path is for **unreachable** only. A
> **content-hash MISMATCH must `error(...)`/abort** — it is tamper/corruption, not recoverable.
> The precedent is `datagen.jl:open_or_invalidate` (below), which `error`s on a stale-cache hash
> mismatch rather than silently reusing. Keep two DISTINCT code paths: `catch` on the download →
> `:skipped`; explicit `!=` on the digest → `error`. Never wrap both in one `try/catch`.

**Hard-error-on-mismatch precedent** (`datagen.jl:368-385`):
```julia
function open_or_invalidate(cache_root, config)
    h         = cache_hash(config)
    dir       = joinpath(cache_root, h)
    meta_path = joinpath(dir, "meta.jld2")
    if isfile(meta_path)
        meta   = JLD2.load(meta_path)
        stored = get(meta, "hash", nothing)
        if stored != h
            old_sub = get(meta, "subhashes", Dict{String,String}())
            new_sub = subhashes()
            changed = [k for k in keys(new_sub) if get(old_sub, k, nothing) != new_sub[k]]
            error("stale cache at $dir: content-hash mismatch " *
                  "(stored=$(stored), recomputed=$(h)); changed sources: $(changed)")
        end
    end
    isdir(dir) || mkpath(dir)
    return dir
end
```

---

### `corpus/manifest.jl` (model + guard, CRUD/transform + access-control)

**Analog:** `table.jl:build_table` (VERIFIED, DataFrame assembly) + `datagen.jl` structural-holdout
discipline (VERIFIED) for the tier/split guards.

**Tidy DataFrame assembly** (`table.jl:175-189`) — one row per image/channel-file, named columns:
```julia
    return DataFrame(
        input_id     = input_id,
        regime       = collect(regime),
        rho_true     = collect(rho_true),
        Costes_p     = Costes_p,
        # ...
        divergence   = divc,
        light        = lightc,
    )
```
Map the RESEARCH §Manifest Schema columns (`anchor_id, tier, role, truth_label,
coloc_ground_truth, accession, source_url, citation, license, format, channels, sha256, split,
bytes`) into this shape, then persist via the `write_table` content-addressed atomic writer.

**Content-addressed atomic CSV+JLD2 writer with pre-registration header** (`table.jl:268-296`) —
copy verbatim for the committed manifest so any contract edit is detectable:
```julia
function write_table(df, config; outdir)
    dir = joinpath(outdir, comparator_hash(config))
    isdir(dir) || mkpath(dir)
    header = _threshold_header(config)
    csv_path = joinpath(dir, "table.csv"); csv_tmp = csv_path * ".tmp"
    open(csv_tmp, "w") do io
        for line in split(header, '\n'; keepempty = false)
            println(io, "# ", line)                    # CSV comment header (pre-registration)
        end
        CSV.write(io, df)
    end
    @assert filesize(csv_tmp) > 0 "csv integrity check failed: $csv_tmp is empty"
    mv(csv_tmp, csv_path; force = true)
    # ... jld2 twin with reopen-integrity, then mv ...
end
```

**Structural tier/split guards (D-06/D-09 — the anti-snooping mechanism):** the discipline to copy
is `datagen.jl`'s "there is NO `standardize_all` symbol, so leakage is impossible by construction"
(see the module header note and the negative-index holdout below). Concretely for the manifest:
- Expose SEPARATE typed accessors — `physical_anchors(m)` vs `cbs_benchmark(m)` — and make any
  tier-blending path a distinct, explicitly-named `mixed_tier(...)` the caller must opt into.
- `sealed_holdout` rows are unreachable from every default accessor; provide an explicit
  `open_sealed_holdout(; reason)` called ONLY by Phase-16 code.

**Structural-disjointness precedent** (`datagen.jl:_write_holdout`, L387-417) — the holdout is a
SEPARATE artifact drawn from an XOR-salted disjoint RNG namespace with **negative** global indices
so it shares no value with the main pool. This is the "structural, not by-value/convention" bar the
split guard must meet:
```julia
        rng = holdout_rng(master_seed, j)                  # disjoint key namespace (D-10)
        # ...
        @inbounds global_index[j]   = -j                   # negative ⇒ never a main-pool index
```
A unit test must assert the default corpus iterator never yields a `sealed_holdout` row.

---

### `corpus/load.jl` (loader/adapter, file-I/O + transform) — READ-ONLY reach into `src/`

**Analog:** `src/LoadImages.jl` (VERIFIED). Reach it via `include()` read-only; do NOT edit `src/`.

**TIFF → `Matrix{Float64}` channel** (`LoadImages.jl:214-217`):
```julia
function load_tiff(path::S) where {S<:AbstractString}
    # Fused operations to reduce allocations
    return Float64.(Images.Gray.(Images.load(path)))
end
```

**Primary constructor** (`LoadImages.jl:126-130`) — asserts channel/threshold length agreement:
```julia
function MultiChannelImage(data::Vector{Matrix{T}}, channels::Vector{S}, name::S, path::Vector{S}, pixel_size::Tuple{I, I}, otsu_threshold::Vector{F}) where {T <: Union{Missing, Float64}, S <: AbstractString, F <: AbstractFloat, I <: Int}
    @assert length(data) == length(channels) "Length of data and channels vectors must match"
    @assert length(otsu_threshold) == length(channels) "Length of otsu_threshold vector must match number of channels"
    new{T, S, F, I}(data, channels, name, path, pixel_size, otsu_threshold)
end
```

**Convenience constructor from paths** (`LoadImages.jl:432-455`) — the SIMPLEST anchor path: give it
`name` + a `Vector` of the two single-channel TIFF paths and it loads, fills `pixel_size =
size(data[1])`, and computes `otsu_threshold = Images.otsu_threshold.(data)`:
```julia
function MultiChannelImage(name::S, path::Vector{S}, channels::Vector{S} =[]) where {S <: AbstractString}
    data = load_tiff.(path)
    if isempty(channels)
        @warn "No channel names provided, using default names"
        channels = ["channel_$i" for i in 1:length(data)]
    end
    pixel_size = size(data[1])
    otsu_threshold = Images.otsu_threshold.(data)
    return MultiChannelImage(data, channels, name, path, pixel_size, otsu_threshold)
end
```
- **Two single-channel TIFFs** (typical anchor layout): pass both paths to this constructor directly.
- **One RGB/composite TIFF** (typical CBS layout): split planes (`Images.channelview` /
  `Images.red`/`green`/`blue` → `Float64`), pick the two per the CBS pair (RG/RB/GB), then call the
  primary constructor.

**Downstream (UNCHANGED) summary path** — once a `MultiChannelImage` exists, `patch(·, 8)` +
`correlation(·)` from `src/colocalization.jl` apply verbatim (VERIFIED signatures):
- `patch(img::AbstractMatrix, num_patches::Integer)` — `colocalization.jl:37`.
- `correlation(x::Array{T,4}, y::Array{T,4}; method=:pearson)` — `colocalization.jl:221`.

> Pitfall 6 (VERIFIED in source): `correlation` calls `_exclude_zero` (`colocalization.jl:187-203`)
> which drops `0.0`/`NaN`/`missing`, and any patch with `length ≤ 15` valid pixels → `missing`
> (`colocalization.jl:235`). Verify on the REAL anchor bytes (not just the fixture) that patches
> retain enough non-zero signal, else the anchor reduces to an all-`missing` 8×8 grid.

---

### `corpus/test/runtests.jl` (test)

**Analog:** `spike/test/runtests.jl` (VERIFIED) — stdlib `Test`, `@testset`/`@test`, sub-suites wired
in via `include`.

**Harness + include-wiring pattern** (`runtests.jl:40-46, 120`):
```julia
using Test
using Pkg

@testset "NeuralEstimators CPU smoke" verbose = true begin
    include(joinpath(@__DIR__, "..", "00_smoke.jl"))
    @testset "ENV-02 correctness (recovered posterior mean)" begin
        @test abs(mu_hat - theta_true) < tol
    end
end
# ...sub-suites wired into the same runner so ONE command is the gate:
include(joinpath(@__DIR__, "test_simulator.jl"))
```
Wire the RESEARCH §Validation Architecture testsets (manifest schema, tier guard, split guard,
integrity hard-fail, skip-with-flag, deterministic hashing, fixture conversion, no-bytes-staged)
into a single `corpus/test/runtests.jl`. The live CBS-fetch smoke must **skip cleanly offline**
(reuse the `tapqir_anchor` `:skipped` contract), so `runtests.jl` never requires network.

---

### synthetic two-channel TIFF fixture (test fixture)

**Analog:** the committed `spike/data/cache/fixture/` idiom (see `.gitignore` L406-409) — a tiny
tracked fixture that lets round-trip/integrity/conversion tests run with ZERO network. Generate the
two-channel TIFF in-test (or commit a tiny tracked file) with enough non-zero per-8×8-patch signal
to survive the `_exclude_zero` ≥15-valid-pixel floor.

---

### `.gitignore` (config) — scoped-negation for the corpus data dir

**Analog:** existing `.gitignore` scoped-negation lines (VERIFIED) — ignore the bulk data dir,
negate the committed fixture:
```gitignore
# Phase-3 DATA-02 bulk training cache (large, regenerable, content-hash-named) — gitignored.
spike/data/cache/*
# …except the tiny committed round-trip / on-disk-schema fixture (D-04 smoke).
!spike/data/cache/fixture/
```
Mirror for the corpus: ignore the downloaded-bytes dir (`corpus/data/*`), negate the manifest,
fetch script, tests, and the tiny fixture. A test asserts no `.tif`/`.zip` under the corpus data
dir is staged (D-04, Pitfall 3).

## Shared Patterns

### SHA-256 content hashing (D-04/D-05/D-08)
**Source:** `spike/comparator/table.jl:comparator_hash` (L225-233) + `_cmp_canonical` (L211-214).
**Apply to:** `corpus/hash.jl`, `corpus/fetch.jl` (download verification), `corpus/manifest.jl`
(manifest content hash). Use stdlib `SHA` (`SHA.sha256` / `SHA256_CTX`) — never hand-roll; never
substitute Base `hash` for the download-integrity digest.

### Atomic durable write (`.tmp → integrity → mv(force=true)`)
**Source:** `src/amortized/datagen.jl:write_shard` (L323-332) + `spike/comparator/table.jl:write_table` (L268-296).
**Apply to:** every persisted corpus artifact — the fetched file commit, the manifest CSV, any JLD2
cache. A crash must never leave a half-written artifact.

### Skip-with-flag graceful degradation (D-05 unreachable = green)
**Source:** `spike/comparator/tapqir_bridge.jl:tapqir_anchor` (L180-205), probe `_tapqir_env_available`
(L114-132), timeout `_wait_timeout` (L61-72).
**Apply to:** `corpus/fetch.jl` (unreachable source → `(status=:skipped, reason=...)`), the live
CBS-fetch test smoke (skips offline). Distinct from the hard-error mismatch path.

### Hard-error on content-hash mismatch (D-05 tamper/corruption = abort)
**Source:** `src/amortized/datagen.jl:open_or_invalidate` (L368-385).
**Apply to:** `corpus/fetch.jl` digest check, `corpus/manifest.jl` contract-hash check. `error(...)`,
never a swallowed skip.

### Seeded reproducibility (Random123 Philox + disjoint salts)
**Source:** `spike/comparator/config.jl:MASTER_SEED` (L56) + `src/amortized/datagen.jl` salts (L56-60)
and keyed RNGs `sample_rng`/`holdout_rng`/`fold_rng` (L81-100).
**Apply to:** any CBS subset selection / dev-eval split assignment — keyed by a corpus
`MASTER_SEED`, bit-reproducible.

### Structural (not conventional) separation for anti-snooping
**Source:** `src/amortized/datagen.jl` holdout discipline — disjoint XOR-salted RNG + negative
global indices (`_write_holdout`, L387-417) and the "no `standardize_all` symbol exists" note.
**Apply to:** the manifest tier guard (D-06, physical-primary vs simulated-secondary) and the split
guard (D-09, sealed_holdout unreachable via default accessors). Enforce by construction + a unit
test, not a comment.

### Read-only reach into frozen `src/`
**Source:** `spike/comparator/table.jl:47-50` guarded-`include` pattern; `datagen.jl` references
`src/` summary functions only inside function bodies.
**Apply to:** `corpus/load.jl` — `include("../src/LoadImages.jl")` / `colocalization.jl` read-only;
NO edits to `src/` (hard CLAUDE.md decoupling constraint).

## No Analog Found

None. Every file this phase needs maps to a VERIFIED in-repo analog. The only genuinely new
primitive is `Downloads.download` (stdlib) for network I/O — it has no repo precedent (not yet used
anywhere), but it is a one-call stdlib primitive composed INTO the verified `fetch_verified`
skeleton (RESEARCH §Fetch + Integrity Mechanics), so no pattern is missing.

## Analogs Named in RESEARCH but Requiring a Caveat (all located, all verified)

| RESEARCH claim | Verification result |
|----------------|---------------------|
| `table.jl:comparator_hash`/`write_table` SHA-256 + atomic | VERIFIED at L225-233 / L268-296 (SHA-256 via `SHA256_CTX`; `.tmp`→check→`mv`). |
| `datagen.jl:cache_hash`/`subhashes`/`write_shard`/`write_meta`/`open_or_invalidate` | VERIFIED L234-242 / L250-252 / L323-332 / L341-357 / L368-385. **Caveat:** `cache_hash` uses Base `hash`, NOT SHA-256 (by design, to avoid Manifest churn) — do NOT reuse it for download integrity; use the `table.jl` SHA-256 idiom there. |
| `tapqir_bridge.jl` skip-with-flag `(status=:skipped\|:passed\|:failed)` + `_wait_timeout` | VERIFIED L180-205 (NamedTuple) / L61-72 (timeout). Note actual return field is `status` with `value`/`reason`, and the passed/failed pair is `:passed`/`:failed`. |
| `config.jl:MASTER_SEED`, pre-declared consts | VERIFIED L56 (+ full declarative-const file). |
| `src/LoadImages.jl` `MultiChannelImage` ctor ~L126, `load_tiff` ~L214 | VERIFIED exactly at L126 and L214. **Bonus:** a convenience path-constructor at L432 is the simplest anchor entry (auto pixel_size + otsu). |
| `spike/test/runtests.jl` `@testset`/`@test` | VERIFIED L40-116 + include-wiring L120-143. |
| DataFrames.jl / CSV.jl already deps | VERIFIED (used in `table.jl`; `runtests.jl` gate (g) asserts both resolved). |

## Conventions

**Convention derivation skipped (no-readable-files):** the shared `gsd-tools.cjs verify conventions
--derive` module returned `{ "skipped": true, "reason": "no-readable-files" }` for `spike/comparator`,
`src`, and repo-wide — it does not parse Julia sources, so no automated axis table is available.

Hand-observed conventions from the verified Julia analogs (planner guidance, non-authoritative):

| Axis | Observed dominant | Notes |
|------|-------------------|-------|
| File-name casing | `snake_case.jl` in `spike/`, `src/amortized/` (e.g. `tapqir_bridge.jl`, `datagen.jl`); PascalCase for a few legacy `src/` files (`LoadImages.jl`) | New `corpus/` files should follow the directory-local `snake_case.jl` majority. |
| Identifier casing | `snake_case` functions (`generate_sample`, `write_shard`, `load_tiff`); `SCREAMING_SNAKE` module consts (`MASTER_SEED`, `HOLDOUT_SALT`, `DIVERGENCE_WARN`) | Match: functions snake_case, pre-declared consts SCREAMING_SNAKE. |
| Module/export style | Include-based composition (`include(...)`), NO `module`/`export` in `spike/`; guarded `isdefined(...) \|\| include(...)` so files load standalone or after a sibling | Reuse the guarded-include idiom (`table.jl:47-50`) for `corpus/` sibling files. |
| Doc/comment style | Triple-quoted docstrings above each function + a file-header `#= license =#` block + a leading `# path --- purpose` line | Match for every new corpus file. |

**Contested hotspots (author's choice):** file-name casing is the one repo-wide contested axis —
`src/` carries PascalCase legacy modules (`LoadImages.jl`) while `spike/` and `src/amortized/` are
uniformly `snake_case.jl`. Each directory is internally consistent; the split is contested only
repo-wide. New `corpus/` files must match the **directory-local** style (snake_case, matching the
`spike/comparator/` prototype this phase mirrors), not deviate to author's choice.

## Metadata

**Analog search scope:** `spike/comparator/`, `src/amortized/`, `src/`, `spike/test/`, `.gitignore`.
**Files scanned (read):** `table.jl`, `tapqir_bridge.jl`, `config.jl`, `datagen.jl`,
`runtests.jl`, `LoadImages.jl`, `colocalization.jl`, `.gitignore`.
**Pattern extraction date:** 2026-07-20
