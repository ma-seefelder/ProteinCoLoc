#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
E-Mail: manuel.seefelder@uni-ulm.de
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany

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

# spike/data/cache.jl --- DATA-02 sharded atomic JLD2 cache (D-04 / D-06).
#
# The persistence half of the pipeline: RAW, column-major (θ, summary_min,
# summary_aug) shards written ATOMICALLY (.tmp + reopen-integrity-check + mv) so
# a crash/resume can never leave a half-written shard in the dataset (T-03-07).
# Resume = skip shards whose final file exists and reloads (`shard_done`),
# regenerating only the in-flight one; raising N just appends new shards without
# touching existing ones (D-06). The cache directory is NAMED by the D-05 content
# hash, and `meta.jld2` stores that hash + per-file sub-hashes + the raw config;
# `open_or_invalidate` recomputes the hash and ERRORs on a stored-vs-recomputed
# mismatch rather than silently reusing stale data (T-03-06).
#
# COLUMN-MAJOR SCHEMA (D-04, matches NeuralEstimators d×K): one sample = one
# column. theta 7×n, summary_min 128×n (D-01), summary_aug AUG_DIM×n (D-02).
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local file I/O; the only src/
# reach is transitive through hashguard.jl's read-only `read` of frozen sources.

using JLD2   # jldsave / jldopen / load (pure-Julia, Windows-clean cache backend)

# hashguard.jl supplies cache_hash / subhashes / HASH_SRC_FILES. Guarded so this
# file is safe to load standalone and after generate.jl already pulled it in.
isdefined(@__MODULE__, :cache_hash) || include(joinpath(@__DIR__, "hashguard.jl"))

# --- Cache constants (D-04) -----------------------------------------------------
# ~10k samples/shard: a shard is (7+128+AUG_DIM) floats × 10k ≈ 10-20 MB in
# memory, trivially fits and writes in one `jldsave` (NO in-file append). Bounds
# the resume granularity (re-doing one shard on crash is cheap) and the parallel
# write unit (one shard per thread, no shared mutable state).
const SHARD_SIZE     = 10_000
const SCHEMA_VERSION = 1   # on-disk schema tag; bump iff the shard key set changes

"""
    shard_path(dir, n::Integer) -> String

The deterministic on-disk path of shard `n` (1-based) inside `dir`:
`shard_0001.jld2`, `shard_0002.jld2`, … (zero-padded so lexical sort == numeric).
"""
shard_path(dir, n::Integer) = joinpath(dir, "shard_" * lpad(n, 4, '0') * ".jld2")

"""
    _loads_ok(path) -> Bool

True iff `path` exists and reopens with the full column-major key set present --
the resume integrity predicate. Any open/read error (a torn `.tmp` left from a
crash) is caught and reported as not-done so the shard is regenerated.
"""
function _loads_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "theta") && haskey(f, "summary_min") &&
                haskey(f, "summary_aug") && haskey(f, "schema_version")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end

"""
    shard_done(path) -> Bool

Resume-by-skip predicate (D-04/D-06): a shard is "done" iff its FINAL file exists
and passes the reopen integrity check. Done shards are skipped on re-run.
"""
shard_done(path) = _loads_ok(path)

"""
    write_shard(dir, n, theta, summary_min, summary_aug, global_index, imsize) -> String

Atomically commit shard `n`: `jldsave` the column-major buffers to a `.tmp`,
reopen to integrity-check, then `mv(...; force=true)` (a filesystem-atomic
rename). A crash before the `mv` leaves only the discardable `.tmp`, never a
half-written final shard (T-03-07). Returns the final shard path.
"""
function write_shard(dir, n::Integer, theta, summary_min, summary_aug,
                     global_index, imsize)
    path = shard_path(dir, n)
    tmp  = path * ".tmp"
    jldsave(tmp; theta, summary_min, summary_aug, global_index, imsize,
            schema_version = SCHEMA_VERSION)
    JLD2.jldopen(tmp, "r") do f                            # integrity check before commit
        @assert haskey(f, "theta") "integrity check failed: $tmp missing theta"
    end
    mv(tmp, path; force = true)                            # atomic commit
    return path
end

"""
    write_meta(dir, config) -> String

Persist the cache manifest into `meta.jld2` (atomically): the full D-05
`cache_hash`, the per-file `subhashes` (for which-source-changed diagnosis), the
raw `config`, N, SHARD_SIZE, and SCHEMA_VERSION. The stored `hash` is the
integrity boundary `open_or_invalidate` re-checks on every reopen.
"""
function write_meta(dir, config)
    path = joinpath(dir, "meta.jld2")
    tmp  = path * ".tmp"
    jldsave(tmp;
        hash           = cache_hash(config),
        subhashes      = subhashes(),
        config         = config,
        N              = config.N,
        shard_size     = SHARD_SIZE,
        schema_version = SCHEMA_VERSION,
    )
    JLD2.jldopen(tmp, "r") do f                            # integrity check before commit
        @assert haskey(f, "hash") "meta integrity check failed: $tmp missing hash"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    open_or_invalidate(cache_root, config) -> String

Resolve the content-hash-named cache directory `joinpath(cache_root,
cache_hash(config))`, creating it if absent. On an EXISTING dir with a stored
`meta.jld2`, recompute the hash and compare to the stored one; on mismatch ERROR
with the diverging source sub-hashes rather than reusing stale data (D-05 /
T-03-06). A changed config or source instead resolves to a DIFFERENT dir (the
hash names it), so versions auto-separate; the in-dir check additionally catches
a corrupted/tampered manifest.
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
