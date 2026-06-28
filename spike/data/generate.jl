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

# spike/data/generate.jl --- DATA-01 in-memory generation core (D-12 / D-03 / D-13).
#
# Composes the UNCHANGED Phase-2 chain (sample_prior -> simulate_pair ->
# build_mci -> patch_summary) with the Wave-2 keyed seeding + encoding into RAW,
# column-major training matrices (one sample = one column, matching the
# NeuralEstimators d×K convention). Persistence (sharding/hashing) is Wave 3;
# standardization is Wave 4 -- here the cache is RAW and in-memory.
#
# REPRODUCIBILITY (D-11/D-12): each sample is keyed by its GLOBAL INDEX, so the
# parallel (`Threads.@threads`) path is byte-identical to the serial fallback for
# ANY thread count, and a shuffled generation order re-sorts to identical columns.
# No shared mutable RNG, no locks -- distinct columns are written independently.
#
# DISTRIBUTION (D-13): θ is sampled i.i.d. from π(θ), no stratification; the
# simulate_pair ArgumentErrors are PROPAGATED (ASVS V5), never swallowed.
#
# DECOUPLING (hard constraint, CLAUDE.md): reaches src/ ONLY transitively through
# contract.jl's read-only include(). All code lives under spike/. The includes
# below are GUARDED (isdefined) so this file is safe to load both standalone and
# inside runtests.jl after test_simulator.jl already loaded the Phase-2 chain.

# --- ORDER MATTERS: contract.jl FIRST (brings build_mci/patch_summary + the
#     read-only src/ coupling into scope), then the simulator, then the siblings.
isdefined(@__MODULE__, :build_mci)     || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :sample_prior)  || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair) || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
include(joinpath(@__DIR__, "seeding.jl"))
include(joinpath(@__DIR__, "encode.jl"))
# Wave-3 persistence layer (sharded cache + D-05 content-hash guard). Guarded so
# generate.jl stays loadable both standalone and after the harness pulled it in.
isdefined(@__MODULE__, :write_shard) || include(joinpath(@__DIR__, "cache.jl"))

# Default master seed (a config value, NOT a secret). Overridable per call.
const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001

"""
    generate_sample(master_seed::Integer, idx::Integer) -> NamedTuple

Generate ONE training sample, keyed entirely by `(master_seed, idx)`:

  rng = sample_rng(master_seed, idx)        # Philox4x keyed stream (D-11)
  θ   = sample_prior(rng)                    # UNCHANGED prior.jl (7-field NamedTuple)
  isz = sample_imsize(rng)                   # D-03 per-sample image size
  mci = build_mci(simulate_pair(rng, θ; imsize=isz))   # UNCHANGED simulator + contract
  M   = patch_summary(mci)                   # UNCHANGED 8×8 summary

Returns `(theta, s_min, s_aug, idx, imsize)` where `theta = collect(values(θ))`
is the 7-vector in prior field order (ρ_true, spillover, autofluorescence,
label_efficiency, shift_dx, shift_dy, noise), `s_min = encode_d01(M)` (128-dim),
`s_aug = encode_aug(mci, M)` (AUG_DIM). `simulate_pair`'s ArgumentErrors propagate.
"""
function generate_sample(master_seed::Integer, idx::Integer)
    rng = sample_rng(master_seed, idx)
    θ   = sample_prior(rng)
    isz = sample_imsize(rng)
    mci = build_mci(simulate_pair(rng, θ; imsize = isz))
    M   = patch_summary(mci)
    return (
        theta  = collect(values(θ)),   # 7-vector, prior field order
        s_min  = encode_d01(M),        # 128-dim D-01
        s_aug  = encode_aug(mci, M),   # AUG_DIM D-02 superset
        idx    = Int(idx),
        imsize = isz,
    )
end

"""
    generate_samples(N::Integer; master_seed=DEFAULT_MASTER_SEED,
                     indices=1:N, parallel=Threads.nthreads()>1) -> NamedTuple

Generate `N` samples into pre-allocated, column-major matrices (one sample = one
column). Column `j` is filled from `generate_sample(master_seed, indices[j])`.

Returns `(theta::Matrix 7×N, summary_min::Matrix 128×N,
summary_aug::Matrix AUG_DIM×N, global_index::Vector{Int} N, imsize::Vector N)`.

Because every column is a pure function of `(master_seed, indices[j])` and each is
written independently, the `parallel` (`Threads.@threads`) and serial paths are
BYTE-IDENTICAL for any thread count, and generating `indices` in any order then
re-sorting by `global_index` yields identical columns (D-11/D-12). θ is i.i.d.
from π(θ) (no stratification, D-13).
"""
function generate_samples(N::Integer; master_seed::Integer = DEFAULT_MASTER_SEED,
                          indices = 1:N, parallel::Bool = Threads.nthreads() > 1)
    idxv = collect(indices)
    @assert length(idxv) == N "indices length $(length(idxv)) != N=$N"

    theta        = Matrix{Float64}(undef, 7, N)
    summary_min  = Matrix{Float64}(undef, 128, N)
    summary_aug  = Matrix{Float64}(undef, AUG_DIM, N)
    global_index = Vector{Int}(undef, N)
    imsize       = Vector{Tuple{Int,Int}}(undef, N)

    # Fill a single column -- no shared mutable state across j (distinct columns).
    fill_col! = function (j::Int)
        s = generate_sample(master_seed, idxv[j])
        @inbounds theta[:, j]       = s.theta
        @inbounds summary_min[:, j] = s.s_min
        @inbounds summary_aug[:, j] = s.s_aug
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

    return (theta = theta, summary_min = summary_min, summary_aug = summary_aug,
            global_index = global_index, imsize = imsize)
end

# ============================ Wave-3 sharded cache driver ========================

"""
    generating_config(N; master_seed, k, shard_size=SHARD_SIZE) -> NamedTuple

The canonical generating config hashed by D-05 (`cache_hash`). Carries every
field whose change must invalidate the cache: N, shard_size, master_seed, k, the
θ-dim, the imsize-spec (`IMSIZE_SET` + `IMSIZE_WEIGHTS`), the summary-spec (D-01
128-dim + D-02 `AUG_DIM`/`N_AUG_MOMENTS`), and the on-disk schema version. Field
ORDER is irrelevant -- `canonical` sorts by name.
"""
function generating_config(N::Integer; master_seed::Integer, k::Integer,
                           shard_size::Integer = SHARD_SIZE)
    return (
        N               = Int(N),
        shard_size      = Int(shard_size),
        master_seed     = UInt64(master_seed),
        k               = Int(k),
        theta_dim       = 7,
        imsize_set      = IMSIZE_SET,
        imsize_weights  = IMSIZE_WEIGHTS,
        summary_min_dim = 128,
        summary_aug_dim = AUG_DIM,
        n_aug_moments   = N_AUG_MOMENTS,
        schema_version  = SCHEMA_VERSION,
    )
end

# Reserved ADVI holdout (D-10): drawn from the PROVABLY DISJOINT, XOR-salted
# `holdout_rng` namespace into a SEPARATE holdout.jld2, NEVER indexed into the
# main pool. Stored `global_index` is NEGATIVE (-1,-2,…) so it shares no value
# with the positive 1:N main pool -- disjointness is structural, not by-value.
function _write_holdout(dir, master_seed::Integer, H::Integer)
    theta        = Matrix{Float64}(undef, 7, H)
    summary_min  = Matrix{Float64}(undef, 128, H)
    summary_aug  = Matrix{Float64}(undef, AUG_DIM, H)
    global_index = Vector{Int}(undef, H)
    imsize       = Vector{Tuple{Int,Int}}(undef, H)
    for j in 1:H
        rng = holdout_rng(master_seed, j)          # disjoint key namespace (D-10)
        θ   = sample_prior(rng)                     # i.i.d. π(θ), UNCHANGED prior.jl
        isz = sample_imsize(rng)
        mci = build_mci(simulate_pair(rng, θ; imsize = isz))
        M   = patch_summary(mci)
        @inbounds theta[:, j]       = collect(values(θ))
        @inbounds summary_min[:, j] = encode_d01(M)
        @inbounds summary_aug[:, j] = encode_aug(mci, M)
        @inbounds global_index[j]   = -j            # negative ⇒ never a main-pool index
        @inbounds imsize[j]         = isz
    end
    path = joinpath(dir, "holdout.jld2")
    tmp  = path * ".tmp"
    jldsave(tmp; theta, summary_min, summary_aug, global_index, imsize,
            schema_version = SCHEMA_VERSION, namespace = "holdout")
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "theta") "holdout integrity check failed: $tmp"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    generate_cache(cache_root; N, master_seed=DEFAULT_MASTER_SEED, k=5,
                   n_holdout=20, shard_size=SHARD_SIZE) -> String

The DATA-02 sharded generation driver. Resolves the content-hash-named cache dir
(`open_or_invalidate`, D-05), partitions `1:N` into `shard_size` chunks, and
`Threads.@threads` over SHARDS (each shard fills column-major buffers via the
keyed `generate_samples` and writes atomically) -- SKIPPING shards already done
(resume-by-skip, D-04). Then draws the reserved ≥20-stack ADVI holdout from the
disjoint `holdout_rng` namespace into a separate `holdout.jld2` (D-10, also
resume-skipped) and writes `meta.jld2`. Raising `N` appends new shards without
regenerating existing ones (D-06). Returns the cache directory.

Launch with `julia -t auto` (or `JULIA_NUM_THREADS`) for parallel generation;
the result is BYTE-IDENTICAL for any thread count (per-global-index keying, D-12).
"""
function generate_cache(cache_root; N::Integer,
                        master_seed::Integer = DEFAULT_MASTER_SEED,
                        k::Integer = 5, n_holdout::Integer = 20,
                        shard_size::Integer = SHARD_SIZE)
    config = generating_config(N; master_seed = master_seed, k = k,
                               shard_size = shard_size)
    dir    = open_or_invalidate(cache_root, config)

    nshards = cld(N, shard_size)
    Threads.@threads for s in 1:nshards
        path = shard_path(dir, s)
        shard_done(path) && continue               # resume-by-skip (D-04/D-06)
        lo  = (s - 1) * shard_size + 1
        hi  = min(s * shard_size, N)
        out = generate_samples(hi - lo + 1; master_seed = master_seed,
                               indices = lo:hi, parallel = false)
        write_shard(dir, s, out.theta, out.summary_min, out.summary_aug,
                    out.global_index, out.imsize)
    end

    # Reserved disjoint-namespace holdout (D-10); skip if already materialized.
    hpath = joinpath(dir, "holdout.jld2")
    _loads_ok(hpath) || _write_holdout(dir, master_seed, max(n_holdout, 20))

    write_meta(dir, config)
    return dir
end
