#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=#

#############################################################################################
# src/amortized/datagen.jl --- grid-parametric training-data generation + atomic cache.
#
# Promoted (grid-generalized) from the proven spike data pipeline:
#   spike/data/seeding.jl   -> per-sample Philox4x keyed RNGs + cost-aware imsize sampler
#   spike/data/generate.jl  -> in-memory column-major generation core + sharded cache driver
#   spike/data/cache.jl     -> atomic (.tmp + reopen-integrity + mv) sharded JLD2 persistence
#
# THE GRID DELTA (D-04): every summary buffer is `Matrix(summary_dim(G), N)` (NOT literal 128),
# and the cache `generating_config` carries `summary_min_dim = summary_dim(G)`, so distinct grids
# hash to distinct content-addressed cache dirs — a 4×4 and a 16×16 cache never collide, and
# raising the grid auto-invalidates cross-grid caches. `imsize_set`/`imsize_weights` are exposed
# so per-grid plans can bias image size upward (fine grids need bigger images to clear the
# ≥15-survivor patch floor).
#
# :MIN-ONLY (scope, D-04): the shipped estimators consume the `:min` D-01 encoding only. The
# spike's `:aug` moment-superset (encode_aug / summary_aug / AUG_DIM) was an ABL-01 slicing
# convenience and is NOT part of the shipped path, so it is dropped here.
#
# SIMULATOR CHAIN: `sample_prior` / `simulate_pair` / `build_mci` are the Phase-2 forward model,
# promoted into src/ by a later Phase-7 plan. They are referenced ONLY inside function bodies
# here, so this skeleton LOADS today and becomes runnable once the simulator is promoted.
#
# REPRODUCIBILITY (D-11/D-12): each sample is a pure function of `(master_seed, grid, idx)`; the
# parallel path is byte-identical to serial for any thread count, and a shuffled generation order
# re-sorts to identical columns.
#############################################################################################

import JLD2
import Random123: Philox4x
import Random: AbstractRNG

# --- Cache constants ------------------------------------------------------------------------
const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001   # a config value, NOT a secret; overridable
const SHARD_SIZE          = 10_000                   # ~10-20 MB/shard; bounds resume + write unit
const SCHEMA_VERSION      = 2                         # on-disk schema tag (2 = grid-parametric, :min-only)

# --- Disjoint Random123 key-namespace salts (D-10) -----------------------------------------
# Fixed nonzero UInt64s XOR-ed into the first key word to carve provably disjoint draw streams
# (SplitMix64 / golden-ratio mixing constants; only nonzero-and-distinct matters).
const HOLDOUT_SALT = 0x9E3779B97F4A7C15   # reserved holdout stream (D-10)
const FOLD_SALT    = 0xD1B54A32D192ED03   # k-fold permutation stream

# --- imsize categorical (D-03) --------------------------------------------------------------
# Discrete size set: all dims ≥ 64 and 8-divisible-friendly; 1376×1028 is the real-data anchor.
# Cost ∝ W·H. `IMSIZE_SET`/`IMSIZE_WEIGHTS` are the DEFAULTS — per-grid plans override via the
# `imsize_set`/`imsize_weights` keyword so fine grids can bias larger (clear the ≥15-px floor).
const IMSIZE_SET = (
    (256, 256),     # 1×    -- the cheap workhorse
    (512, 512),     # 4×
    (1024, 1024),   # 16×
    (1376, 1028),   # ~21.6× -- real-data anchor
    (2048, 2048),   # 64×   -- the expensive tail
)
const IMSIZE_WEIGHTS = (0.55, 0.35, 0.05, 0.03, 0.02)   # Σ = 1.0; ≥1024² fraction = 0.10

"""
    sample_rng(master_seed, idx) -> Philox4x

The main-pool per-sample keyed RNG (D-11): the draw stream depends ONLY on `(master_seed, idx)`,
independent of thread or execution order. An `AbstractRNG` usable directly in the simulator.
"""
sample_rng(master_seed::Integer, idx::Integer) =
    Philox4x(UInt64, (UInt64(master_seed), UInt64(idx)))

"""
    holdout_rng(master_seed, idx) -> Philox4x

The reserved-holdout per-sample keyed RNG (D-10): `HOLDOUT_SALT` XOR-ed into the first key word,
so its stream is PROVABLY DISJOINT from `sample_rng(master_seed, idx)` for the same `idx`.
"""
holdout_rng(master_seed::Integer, idx::Integer) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ HOLDOUT_SALT, UInt64(idx)))

"""
    fold_rng(master_seed) -> Philox4x

The k-fold permutation RNG, keyed by `(master_seed ⊻ FOLD_SALT, 0)`, reproducible and
independent of the generation key stream.
"""
fold_rng(master_seed::Integer) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ FOLD_SALT, UInt64(0)))

"""
    sample_imsize(rng; imsize_set=IMSIZE_SET, imsize_weights=IMSIZE_WEIGHTS) -> Tuple{Int,Int}

Draw one image size from `imsize_set` under the cost-aware `imsize_weights` (D-03), reproducible
under a fixed `rng`. Hand-rolled inverse-CDF over the fixed categorical (no extra dependency).
"""
function sample_imsize(rng::AbstractRNG; imsize_set = IMSIZE_SET,
                       imsize_weights = IMSIZE_WEIGHTS)
    r = rand(rng)
    c = 0.0
    @inbounds for i in 1:length(imsize_set)
        c += imsize_weights[i]
        r <= c && return imsize_set[i]
    end
    return imsize_set[end]   # float-rounding guard: r just over Σw -> the last bin
end

# ============================ in-memory generation core =====================================

"""
    generate_sample(master_seed, idx, grid; imsize_set, imsize_weights) -> NamedTuple

Generate ONE training sample for a `grid×grid` summary, keyed by `(master_seed, grid, idx)`:

  rng = sample_rng(master_seed, idx)
  θ   = sample_prior(rng)                                    # Phase-2 prior (7-field NamedTuple)
  isz = sample_imsize(rng; imsize_set, imsize_weights)       # per-sample image size (D-03)
  mci = build_mci(simulate_pair(rng, θ; imsize = isz))       # Phase-2 forward model
  M   = patch_summary(mci, grid)                             # grid-parametric GxG summary

Returns `(theta, s_min, idx, imsize)` where `theta = collect(values(θ))` (7-vector, prior field
order) and `s_min = encode_d01(M)` (length `summary_dim(grid) = 2·grid²`).
"""
function generate_sample(master_seed::Integer, idx::Integer, grid::Integer;
                         imsize_set = IMSIZE_SET, imsize_weights = IMSIZE_WEIGHTS)
    rng = sample_rng(master_seed, idx)
    θ   = sample_prior(rng)
    isz = sample_imsize(rng; imsize_set = imsize_set, imsize_weights = imsize_weights)
    mci = build_mci(simulate_pair(rng, θ; imsize = isz))
    M   = patch_summary(mci, grid)                          # grid-parametric summary
    return (
        theta  = collect(values(θ)),          # 7-vector, prior field order
        s_min  = encode_d01(M),               # summary_dim(grid)-dim D-01 encoding
        idx    = Int(idx),
        imsize = isz,
    )
end

"""
    generate_samples(N; grid, master_seed, indices, imsize_set, imsize_weights, parallel) -> NamedTuple

Generate `N` samples for `grid` into pre-allocated, column-major matrices (one sample = one
column). The summary buffer is `Matrix(summary_dim(grid), N)` — the grid coupling threaded from
the single-source dimension helper, NOT a literal 128.

Returns `(theta::Matrix 7×N, summary_min::Matrix summary_dim(grid)×N, global_index::Vector{Int},
imsize::Vector)`. Every column is a pure function of `(master_seed, grid, indices[j])`, so the
`parallel` and serial paths are BYTE-IDENTICAL for any thread count (D-11/D-12).
"""
function generate_samples(N::Integer; grid::Integer,
                          master_seed::Integer = DEFAULT_MASTER_SEED,
                          indices = 1:N,
                          imsize_set = IMSIZE_SET, imsize_weights = IMSIZE_WEIGHTS,
                          parallel::Bool = Threads.nthreads() > 1)
    idxv = collect(indices)
    @assert length(idxv) == N "indices length $(length(idxv)) != N=$N"

    d            = summary_dim(grid)                        # grid coupling (single source)
    theta        = Matrix{Float64}(undef, 7, N)
    summary_min  = Matrix{Float64}(undef, d, N)            # summary_dim(grid)×N, NOT 128×N
    global_index = Vector{Int}(undef, N)
    imsize       = Vector{Tuple{Int,Int}}(undef, N)

    fill_col! = function (j::Int)
        s = generate_sample(master_seed, idxv[j], grid;
                            imsize_set = imsize_set, imsize_weights = imsize_weights)
        @inbounds theta[:, j]       = s.theta
        @inbounds summary_min[:, j] = s.s_min
        @inbounds global_index[j]   = s.idx
        @inbounds imsize[j]         = s.imsize
        return nothing
    end

    if parallel
        Threads.@threads for j in 1:N
            fill_col!(j)
        end
    else
        for j in 1:N
            fill_col!(j)
        end
    end

    return (theta = theta, summary_min = summary_min,
            global_index = global_index, imsize = imsize)
end

# ============================ content-hash version guard ====================================
# A hash over the data-DEFINING src source bytes PLUS a canonical serialization of the config.
# The digest names the cache directory, so ANY change to the summary source OR to a hashed config
# field (incl. `summary_min_dim = summary_dim(grid)`) yields a different digest and auto-separates
# stale/cross-grid data. Uses Base `hash` (no new stdlib dependency) so the fragile Wave-0
# co-resolution Manifest stays untouched.

# FIXED-order data-defining src files (@__DIR__ = src/amortized/). The simulator files join this
# list when a later plan promotes them; the list is a keyword so callers can extend it.
const DATAGEN_HASH_SRC_FILES = String[
    joinpath(@__DIR__, "..", "colocalization.jl"),   # patch() / correlation() (defines the summary)
    joinpath(@__DIR__, "..", "LoadImages.jl"),        # MultiChannelImage container
    joinpath(@__DIR__, "summary.jl"),                 # patch_summary(mci,G) / encode_d01
]

"""
    _canonical(config::NamedTuple) -> String

Deterministic, sort-by-field-name serialization of `config` so field ORDER can never perturb the
hash. Each value is rendered with `repr` (stable for the Int/UInt64/Tuple/Symbol fields carried).
"""
function _canonical(config::NamedTuple)::String
    ks = sort(collect(keys(config)))
    return join(("$(k)=$(repr(getproperty(config, k)))" for k in ks), ";")
end

"""
    cache_hash(config; src_files=DATAGEN_HASH_SRC_FILES) -> String

The content hash: fold Base `hash` over the raw bytes of each `src_files` entry (FIXED order,
NUL-separated) then over `_canonical(config)`, returning a padded hex digest. Identical
`(config, source bytes)` ⇒ identical digest; perturbing ANY source byte OR ANY hashed config
field (incl. the grid via `summary_min_dim`) flips it. This digest names the cache directory.
"""
function cache_hash(config::NamedTuple; src_files = DATAGEN_HASH_SRC_FILES)
    h = zero(UInt)
    for f in src_files                                     # FIXED order = deterministic
        h = hash(read(f), h)
        h = hash(0x00, h)                                  # inter-file separator
    end
    h = hash(_canonical(config), h)
    return string(h; base = 16, pad = 16)
end

"""
    subhashes(src_files=DATAGEN_HASH_SRC_FILES) -> Dict{String,String}

Per-file digests keyed by source path, so a `cache_hash` mismatch is DIAGNOSABLE — the caller can
diff stored vs recomputed sub-hashes to name exactly WHICH source file changed.
"""
function subhashes(src_files = DATAGEN_HASH_SRC_FILES)
    return Dict{String,String}(f => string(hash(read(f)); base = 16, pad = 16) for f in src_files)
end

# ============================ sharded atomic JLD2 cache =====================================

"""
    generating_config(N; grid, master_seed, k, shard_size, imsize_set, imsize_weights) -> NamedTuple

The canonical generating config hashed by `cache_hash`. Carries every field whose change must
invalidate the cache: N, shard_size, master_seed, k, grid, θ-dim, the imsize-spec, the
grid-coupled `summary_min_dim = summary_dim(grid)` (so distinct grids auto-separate into distinct
content-hashed dirs), and the on-disk schema version. Field ORDER is irrelevant (`_canonical`
sorts by name).
"""
function generating_config(N::Integer; grid::Integer, master_seed::Integer, k::Integer,
                           shard_size::Integer = SHARD_SIZE,
                           imsize_set = IMSIZE_SET, imsize_weights = IMSIZE_WEIGHTS)
    return (
        N               = Int(N),
        shard_size      = Int(shard_size),
        master_seed     = UInt64(master_seed),
        k               = Int(k),
        grid            = Int(grid),
        theta_dim       = 7,
        imsize_set      = imsize_set,
        imsize_weights  = imsize_weights,
        summary_min_dim = summary_dim(grid),      # grid coupling ⇒ per-grid cache separation
        schema_version  = SCHEMA_VERSION,
    )
end

"""
    shard_path(dir, n) -> String

Deterministic on-disk path of shard `n` (1-based): `shard_0001.jld2`, … (zero-padded so a lexical
sort equals a numeric sort).
"""
shard_path(dir, n::Integer) = joinpath(dir, "shard_" * lpad(n, 4, '0') * ".jld2")

"""
    _loads_ok(path) -> Bool

True iff `path` exists and reopens with the full column-major key set present — the resume
integrity predicate. Any open/read error (a torn `.tmp` from a crash) is reported not-done.
"""
function _loads_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "theta") && haskey(f, "summary_min") && haskey(f, "schema_version")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end

"""
    shard_done(path) -> Bool

Resume-by-skip predicate: a shard is "done" iff its FINAL file exists and passes the reopen
integrity check. Done shards are skipped on re-run.
"""
shard_done(path) = _loads_ok(path)

"""
    write_shard(dir, n, theta, summary_min, global_index, imsize) -> String

Atomically commit shard `n`: `jldsave` the column-major buffers to a `.tmp`, reopen to
integrity-check, then `mv(...; force=true)` (a filesystem-atomic rename). A crash before the `mv`
leaves only the discardable `.tmp`, never a half-written final shard. Returns the final path.
"""
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

"""
    write_meta(dir, config) -> String

Persist the cache manifest into `meta.jld2` (atomically): the full `cache_hash`, per-file
`subhashes` (for which-source-changed diagnosis), the raw `config`, N, SHARD_SIZE, and
SCHEMA_VERSION. The stored `hash` is the integrity boundary `open_or_invalidate` re-checks.
"""
function write_meta(dir, config)
    path = joinpath(dir, "meta.jld2")
    tmp  = path * ".tmp"
    JLD2.jldsave(tmp;
        hash           = cache_hash(config),
        subhashes      = subhashes(),
        config         = config,
        N              = config.N,
        shard_size     = SHARD_SIZE,
        schema_version = SCHEMA_VERSION,
    )
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "hash") "meta integrity check failed: $tmp missing hash"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    open_or_invalidate(cache_root, config) -> String

Resolve the content-hash-named cache directory `joinpath(cache_root, cache_hash(config))`,
creating it if absent. On an EXISTING dir with a stored `meta.jld2`, recompute the hash and
compare; on mismatch ERROR with the diverging source sub-hashes rather than reusing stale data. A
changed config/source resolves to a DIFFERENT dir (the hash names it), so versions — and grids —
auto-separate.
"""
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

# Reserved holdout (D-10): drawn from the PROVABLY DISJOINT, XOR-salted `holdout_rng` namespace
# into a SEPARATE holdout.jld2, never indexed into the main pool. `global_index` is NEGATIVE so
# it shares no value with the positive 1:N main pool — disjointness is structural, not by-value.
function _write_holdout(dir, master_seed::Integer, H::Integer, grid::Integer;
                        imsize_set = IMSIZE_SET, imsize_weights = IMSIZE_WEIGHTS)
    d            = summary_dim(grid)
    theta        = Matrix{Float64}(undef, 7, H)
    summary_min  = Matrix{Float64}(undef, d, H)
    global_index = Vector{Int}(undef, H)
    imsize       = Vector{Tuple{Int,Int}}(undef, H)
    for j in 1:H
        rng = holdout_rng(master_seed, j)                  # disjoint key namespace (D-10)
        θ   = sample_prior(rng)
        isz = sample_imsize(rng; imsize_set = imsize_set, imsize_weights = imsize_weights)
        mci = build_mci(simulate_pair(rng, θ; imsize = isz))
        M   = patch_summary(mci, grid)
        @inbounds theta[:, j]       = collect(values(θ))
        @inbounds summary_min[:, j] = encode_d01(M)
        @inbounds global_index[j]   = -j                   # negative ⇒ never a main-pool index
        @inbounds imsize[j]         = isz
    end
    path = joinpath(dir, "holdout.jld2")
    tmp  = path * ".tmp"
    JLD2.jldsave(tmp; theta, summary_min, global_index, imsize,
                 schema_version = SCHEMA_VERSION, namespace = "holdout")
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "theta") "holdout integrity check failed: $tmp"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    generate_cache(cache_root; grid, N, master_seed, k, n_holdout, shard_size,
                   imsize_set, imsize_weights) -> String

The sharded generation driver for a `grid×grid` summary. Resolves the content-hash-named cache
dir (`open_or_invalidate`), partitions `1:N` into `shard_size` chunks, and `Threads.@threads`
over SHARDS (each fills column-major buffers via keyed `generate_samples` and writes atomically),
SKIPPING shards already done (resume-by-skip). Then draws the reserved holdout from the disjoint
`holdout_rng` namespace and writes `meta.jld2`. Distinct grids resolve to distinct cache dirs
(the grid enters the content hash via `summary_min_dim`). Returns the cache directory.

Launch with `julia -t auto` for parallel generation; the result is BYTE-IDENTICAL for any thread
count (per-global-index keying, D-12).
"""
function generate_cache(cache_root; grid::Integer, N::Integer,
                        master_seed::Integer = DEFAULT_MASTER_SEED,
                        k::Integer = 5, n_holdout::Integer = 20,
                        shard_size::Integer = SHARD_SIZE,
                        imsize_set = IMSIZE_SET, imsize_weights = IMSIZE_WEIGHTS)
    config = generating_config(N; grid = grid, master_seed = master_seed, k = k,
                               shard_size = shard_size,
                               imsize_set = imsize_set, imsize_weights = imsize_weights)
    dir    = open_or_invalidate(cache_root, config)

    nshards = cld(N, shard_size)
    Threads.@threads for s in 1:nshards
        path = shard_path(dir, s)
        shard_done(path) && continue                       # resume-by-skip
        lo  = (s - 1) * shard_size + 1
        hi  = min(s * shard_size, N)
        out = generate_samples(hi - lo + 1; grid = grid, master_seed = master_seed,
                               indices = lo:hi, imsize_set = imsize_set,
                               imsize_weights = imsize_weights, parallel = false)
        write_shard(dir, s, out.theta, out.summary_min, out.global_index, out.imsize)
    end

    hpath = joinpath(dir, "holdout.jld2")
    _loads_ok(hpath) || _write_holdout(dir, master_seed, max(n_holdout, 20), grid;
                                       imsize_set = imsize_set, imsize_weights = imsize_weights)

    write_meta(dir, config)
    return dir
end
