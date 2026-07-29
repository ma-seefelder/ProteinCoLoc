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

# spike/data/p12_generate.jl --- Phase-12 FIELD-AWARE training-pool generation (D-05, D-07, D-08).
#
# THE GENERATIVE STORY, IN EXACTLY THIS ORDER. FOR EVERY TRAINING SAMPLE i:
#     r1_i        ~ Uniform(P12_R1_MIN, P12_R1_MAX)      <- THE CORRELATION LENGTH, DRAWN **FIRST**
#     z_field_i   ~ N(0, Sigma(arm, r1_i)) then rescaled <- THE FIELD, CONDITIONAL ON r1_i
#     rho_field_i  = ghat.(quantile.(MU_PRIOR, Phi(z)))  <- D-05 copula, marginal-preserving
#     7 nuisances ~ their existing priors                <- UNCHANGED from prior.jl
#     imsize_i    ~ THE F5 CATEGORICAL                   <- READ FROM THE FROZEN PRE-REGISTRATION
#     y_i          = the 7-stage forward model at imsize_i, with rho supplied PER PIXEL
#
# THIS ORDER IS THE ENTIRE POINT OF THE FILE. IT IS `p11_generate.jl:20-40` TRANSPOSED ONTO A
# FIELD. IF THE FIELD WERE DRAWN FIRST AND A CORRELATION LENGTH REPORTED FOR IT AFTERWARDS, r1
# WOULD CARRY NO INFORMATION ABOUT THE FIELD, THE NETWORK WOULD **CORRECTLY** LEARN TO IGNORE THE
# r1 ROW, AND D-08 WOULD FAIL FLAT FOR A REASON THAT HAS NOTHING TO DO WITH RESOLUTION -- A NULL
# THAT LOOKS EXACTLY LIKE "BELOW THE FIXED 8x8 SUMMARY'S RESOLUTION" WHILE ACTUALLY BEING AN
# IMPLEMENTATION BUG. THOSE TWO READINGS HAVE COMPLETELY DIFFERENT CONSEQUENCES. The order is
# enforced upstream inside `sample_p12_prior` (p12_prior.jl:184-207), which is the ONLY place a
# Phase-12 draw is made; this file never re-implements a draw, so it cannot get the order wrong
# by a second route.
#
# THE STORED GROUND TRUTH IS THE **DRAWN G x G LATTICE** (D-07). Not the cell mean of the
# interpolated pixel field, and not the induced per-cell patch correlation. D-07 scores the drawn
# lattice values because those are the SAMPLED parameters: a cell mean is a lossy projection, and
# an induced per-cell correlation is a MEASURED, NOISY quantity -- SBC ranks parameters, never
# measurements. Both rejected alternatives are asserted against in `test_p12_datagen.jl` testset 2.
#
# WHAT IS STORED, AND WHY BOTH FIELDS. `z_field` is the GAUSSIAN field: it is what the theta rows
# are built from (R-2), and in Gaussian space there are no `ghat` atoms, so it is the clean SBC
# target. `rho_field` is its monotone copula image: it is the deliverable's UNITS and 12-19 needs
# it, but it carries a measured 6.79 % clamp atom per region and is NOT a theta row. `mu_field` is
# kept so the copula step is checkable after the fact rather than only at generation time.
#
# THIS POOL IS SINGLE-STACK, AND THAT IS CORRECT -- Delta-rho IS A READ-TIME CONSTRUCTION, NOT A
# POOL PROPERTY. The trained net maps ONE summary to theta; Delta-rho is formed AFTERWARDS from two
# independent single-stack passes of that same net (`src/amortized/infer.jl:108`, D-03), exactly as
# `spike/validation/harness.jl:124-136` already does for the scalar lane: it draws two priors,
# encodes both, and lets the consumer form the difference with truth
# `Delta-rho* = theta_s.rho - theta_c.rho`. The pairing therefore happens at SCORING time, in 12-16
# and 12-18, which draw FRESH datasets through the harness rather than reading this pool. Two
# consequences, recorded so nobody re-adds a paired control draw here:
#   * DATAGEN DOES NOT DOUBLE. A paired pool would mean two forward simulations and two summary
#     passes per dataset, taking the research's measured 65-80 min / 50k to ~130-160 min against
#     `P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150` -- straddling a ceiling that Tier 1 makes
#     APPEND-ONLY and that this file declares a BLOCKER rather than a downgrade trigger.
#   * NOTHING IS LOST. The per-region Delta-rho* truth 12-18 ranks is `rho*_sample - rho*_control`
#     from the paired draw 12-18 makes ITSELF, so it exists without this pool storing anything
#     extra.
#
# THE THETA ROWS ARE THE **GAUSSIAN** FIELD'S DCT COEFFICIENTS (R-2), NOT THE rho FIELD'S. Row 1 is
# the DC coefficient c0 and rows 2..64 are the deviations in the frozen `p12_dct_order` smoothness
# order. Delta-rho is NOT a theta row either: theta is 72 rows, and Delta-rho is a DERIVED quantity
# carried as a separate appended rank column by its consumer.
#
# WHAT THIS FILE DELIBERATELY DOES **NOT** DO. It stores the RAW 128-row summary and the RAW drawn
# lattice. It does not z-score anything, it does not re-view the 128 rows as a (G, G, 2, K) array,
# and it does not mask regions. All three are the trainer's TRAIN-ONLY jobs: fitting the frozen
# summary transform here would leak the validation fold into the fit (leak-freedom), and the MCAR
# region-masking augmentation (Pitfall 4, `P12_MASK_K_SET`) is applied by 12-14 so that ONE pool can
# serve both the masked and the unmasked arm. The realized masking RATE is therefore recorded at
# TRAINING time, not here.
#
# AND THE TRAINER MASKS THE **RAW** ROWS, BEFORE ANY Z-SCORING. Recorded here because this header is
# where a wrong ordering would get enshrined two waves before 12-14 runs. A patch whose correlation
# cannot be scored reaches the net as raw `0.0` with mask `0`; after per-row z-scoring that becomes
# `(0 - mu_r)/sd_r`, which is NOT zero. Masking AFTER the z-score would instead place the row MEAN
# at a held-out region, so the net would read "perfectly average" where inference time means
# "absent", and D-09 would measure out-of-distribution behaviour instead of borrowing strength.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local generation; reaches src/ only READ-ONLY,
# transitively through `contract.jl`'s include and through the content-hash guard's `read` of
# frozen source bytes. Nothing here writes to, or imports from, the ProteinCoLoc module, and
# nothing here adds a package to spike/Project.toml.
#
# Run (self-check):
#     julia --project=spike -t auto -e 'include("spike/data/p12_generate.jl")'

# --- ORDER MATTERS. The Tier-1 pre-registration FIRST (it binds P12_DEV_SEED, P12_DATAGEN_SALT,
#     P12_G, P12_K_DEV, P12_D_MINISPIKE, the r1 bracket, the F5 mixture and the wall-clock ceiling,
#     and it imports Philox4x at UInt64 width), then the frozen summary contract, then the
#     Phase-12 prior (which transitively pulls prior.jl, ghat.jl and p12_lattice.jl), the forward
#     model, the row LAYOUT, the encoder, and the sharded-cache persistence layer.
#
#     A NOTE ON GUARD SENTINELS, BECAUSE THIS FILE IS THE ONE THAT LOADS THE MOST OF THEM.
#     `p12_consts.jl` itself BINDS `P11_DEV_SEED` and `P13_DEV_SEED` (its forbidden-seed
#     inventory), which are the sentinels `p11_consts.jl` and `spike/p13/consts.jl` guard their own
#     Tier-1 blocks on. A process that wants a FOREIGN pre-registration as well must therefore load
#     the foreign one FIRST. This file needs none of them, so it guards on Phase-12-local names
#     only and never on a name another phase owns.
isdefined(@__MODULE__, :P12_DEV_SEED)      || include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))
isdefined(@__MODULE__, :build_mci)         || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :sample_p12_prior)  || include(joinpath(@__DIR__, "..", "simulator", "p12_prior.jl"))
isdefined(@__MODULE__, :simulate_pair)     || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :p12_dct_vec)       || include(joinpath(@__DIR__, "..", "simulator", "p12_lattice.jl"))
isdefined(@__MODULE__, :p12_theta_index)   || include(joinpath(@__DIR__, "..", "npe", "p12_architecture.jl"))
isdefined(@__MODULE__, :encode_d01)        || include(joinpath(@__DIR__, "encode.jl"))
isdefined(@__MODULE__, :write_shard)       || include(joinpath(@__DIR__, "cache.jl"))

# =============================================================================================
# 1. The reserved Phase-12 datagen stream (R-5)
# =============================================================================================
#
# `P12_DATAGEN_SALT` is Tier-1 (p12_consts.jl:256) rather than declared here, which is the ONE
# structural difference from the Phase-11 analog: Phase 12's pre-registration was written after
# Phase 11 shipped, so the salt inventory and the salt itself live in the same frozen file and the
# disjointness is asserted there BEFORE any generation code exists. The assertions below are
# deliberately redundant with that file's own: they prove the property AT THE POINT OF USE, so a
# future edit that reached for a different salt here would fail on load.
#
# The stream is carved by its SALT, not by a reserved counter, because the counter word must carry
# the per-sample GLOBAL INDEX -- that is what makes the pool byte-identical for any thread count.
# `P12_DATAGEN_COUNTER` is therefore recorded in the pool's generating config as the audit-trail
# tag for this activity rather than consumed as a key word.
@assert P12_DATAGEN_SALT ∉ P12_REPO_SALTS "P12_DATAGEN_SALT reuses an existing repository salt"
@assert P12_DATAGEN_SALT != P12_SALT
@assert (UInt64(P12_DEV_SEED) ⊻ P12_DATAGEN_SALT) != (UInt64(P12_DEV_SEED) ⊻ P12_SALT)
@assert !(UInt64(P12_DEV_SEED) in _p12_forbidden())

"""
    p12_datagen_rng(idx::Integer) -> Philox4x

The per-sample keyed RNG of the Phase-12 training pool:
`Philox4x((P12_DEV_SEED xor P12_DATAGEN_SALT, idx))`. Counter-based and keyed by the GLOBAL INDEX,
so the draw stream is a pure function of the index -- the pool is BYTE-IDENTICAL regardless of
thread count, execution order, or how many times a killed run is resumed. Mirrors
`p11_datagen_rng` (`spike/data/p11_generate.jl:113-114`).
"""
p12_datagen_rng(idx::Integer) =
    Philox4x(UInt64, (UInt64(P12_DEV_SEED) ⊻ P12_DATAGEN_SALT, UInt64(idx)))

# On-disk schema tag for a Phase-12 pool shard. Bump iff the shard key set changes.
const P12_POOL_SCHEMA = 1

# The RAW encoded-summary height, DERIVED from the frozen grid rather than retyped as 128: rows
# 1:G^2 are the imputed patch correlations and rows G^2+1:2G^2 the binary present-mask.
const P12_SUMMARY_MIN_DIM = 2 * P12_G^2

@assert P12_SUMMARY_MIN_DIM == 128

# =============================================================================================
# 2. THE ONLY theta ASSEMBLER IN THE PHASE
# =============================================================================================

"""
    p12_theta_column(draw; G = P12_G, K = P12_K_DEV) -> Vector{Float64}

Turn one `sample_p12_prior` draw into ONE column of the theta matrix, in `P12_THETA_ROWS` order:

    [ c0 ; c1..cK ; spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise,
      chromatic_eps ; r1 ]

`c` is the 2-D DCT-II of the **GAUSSIAN** field, permuted into the frozen `p12_dct_order`
smoothness ordering, so `c[1]` is the DC coefficient and `c[2:1+K]` are the deviations. The rho
field is NOT a theta row: in rho space the elementwise `ghat` clamp puts a measured 6.79 % atom
mass on every region, while the Gaussian field is atom-free (R-2).

THIS FUNCTION IS THE ONLY PLACE IN PHASE 12 THAT BUILDS A theta COLUMN, and that is deliberate.
`p12_architecture.jl` (12-04) owns the row LAYOUT, `p12_prior.jl` (12-07) owns the DRAW, and this
file is the only one that includes both -- so the layout is honoured exactly once. A second
assembler anywhere would be a second chance to go off by one row, and an off-by-one here makes
every SBC table, every coverage number and every ablation comparison in this phase silently about
the wrong parameter while every other test in the suite still passes.

The guard below is therefore keyed off `p12_theta_index` (the layout) rather than off arithmetic:
a layout change cannot silently pass here, because the assertion reads the layout it is checking.
"""
function p12_theta_column(draw; G::Integer = P12_G, K::Integer = P12_K_DEV)
    G >= 2 || throw(ArgumentError("p12_theta_column: G must be ≥ 2, got $G"))
    (0 <= K <= G^2 - 1) || throw(ArgumentError(
        "p12_theta_column: K must satisfy 0 ≤ K ≤ G²-1 = $(G^2 - 1), got $K"))
    size(draw.z_field) == (G, G) || throw(DimensionMismatch(
        "p12_theta_column: expected a $G×$G z_field, got $(size(draw.z_field))"))

    c = p12_dct_vec(vec(draw.z_field); G = G)[p12_dct_order(G)]

    nuis = (draw.spillover, draw.autofluorescence, draw.label_efficiency,
            draw.shift_dx, draw.shift_dy, draw.noise, draw.chromatic_eps)

    out = Vector{Float64}(undef, 1 + K + 7 + 1)
    out[1] = c[1]                                   # the global DC term
    @inbounds for k in 1:K
        out[1 + k] = c[1 + k]                       # the deviations, in p12_dct_order
    end
    @inbounds for (t, v) in enumerate(nuis)
        out[1 + K + t] = v                          # the 7 nuisances, EXISTING relative order
    end
    out[1 + K + 7 + 1] = draw.r1                    # APPENDED LAST (D-08)

    # --- The structural guard, read off the LAYOUT and not off arithmetic ---------------------
    @assert length(out) == 1 + K + 7 + 1
    K == P12_K_DEV && @assert length(out) == P12_D_MINISPIKE "p12_theta_column: at the " *
        "pre-registered K = $(P12_K_DEV) the column must be P12_D_MINISPIKE = $(P12_D_MINISPIKE) " *
        "rows, got $(length(out))"
    @assert out[p12_theta_index(:c0; K = K)] == c[1]
    @assert out[p12_theta_index(:r1; K = K)] == draw.r1
    @assert out[p12_theta_index(:chromatic_eps; K = K)] == draw.chromatic_eps
    @assert out[p12_theta_index(:spillover; K = K)]        == draw.spillover
    @assert out[p12_theta_index(:autofluorescence; K = K)] == draw.autofluorescence
    @assert out[p12_theta_index(:label_efficiency; K = K)] == draw.label_efficiency
    @assert out[p12_theta_index(:shift_dx; K = K)]         == draw.shift_dx
    @assert out[p12_theta_index(:shift_dy; K = K)]         == draw.shift_dy
    @assert out[p12_theta_index(:noise; K = K)]            == draw.noise
    return out
end

# =============================================================================================
# 3. The image-size mixture (F5)
# =============================================================================================

"""
    sample_p12_imsize(rng::AbstractRNG) -> Tuple{Int,Int}

Draw one image size from the F5 mixture, `P12_IMSIZE_SET` under `P12_IMSIZE_WEIGHTS`.

BOTH ARE **READ** FROM `spike/validation/p12_consts.jl`, which itself reads them from the frozen
amended pre-registration through its isolated `GateV2P12` module. THE TUPLE IS NEVER RETYPED HERE
-- a retyped mixture is a mixture that can silently drift, and F5's binding invariant is
train-joint == eval-joint.

WHY A MIXTURE AND NOT ONE CONVENIENT SIZE. Pitfall 2 measured the per-region noise sd as
IMAGE-SIZE DEPENDENT: a patch at 512² carries roughly one seventh of the pixels of the same patch
at 1376×1028, so its correlation estimate is correspondingly noisier. A single-size pool would
therefore put the per-region signal-to-noise ratio somewhere the net will never see at read time,
and every per-region interval would be calibrated for the wrong noise level. Hand-rolled
inverse-CDF over the fixed categorical, mirroring `sample_p11_imsize`
(`spike/data/p11_generate.jl:202-210`); no extra dependency.
"""
function sample_p12_imsize(rng::AbstractRNG)
    r = rand(rng)
    c = 0.0
    @inbounds for i in eachindex(P12_IMSIZE_SET)
        c += P12_IMSIZE_WEIGHTS[i]
        r <= c && return P12_IMSIZE_SET[i]
    end
    return P12_IMSIZE_SET[end]   # float-rounding guard: r just over the weight sum -> last bin
end

# =============================================================================================
# 4. One sample, one block
# =============================================================================================

"""
    generate_p12_sample(idx; imsize_sampler, arm, imsize, r1) -> NamedTuple

Generate ONE pool row, keyed entirely by its global index `idx`:

    rng    = p12_datagen_rng(idx)
    draw   = a Phase-12 prior draw -- r1 FIRST, then the field conditional on it
    isz    = imsize_sampler(rng)                 # the F5 categorical, unless `imsize` is pinned
    Z      = the RAW 128-row encoding of the patch grid of the simulated pair

`r1` IS THREADED THROUGH TO THE DRAW, NEVER DROPPED. A pinned `r1` is what 12-11's D-12 Stage-1
ladder and the D-10 ablation are built on.

THE DRAW COMES BEFORE THE IMAGE SIZE, AND THE HEADER'S "DRAWN **FIRST**" IS MEANT LITERALLY. The
first numbers this key produces are r1 and the field; the image size is a nuisance of the
OBSERVATION, not of the parameter, so it is consumed afterwards -- exactly the ordering the
Phase-11 analog uses (`p11_generate.jl:226-238`, theta then imsize). Pinning `imsize` skips the
sampler entirely and therefore consumes less of `rng` than the mixture path does; that is why a
pinned-size caller passes a sampler closure instead whenever the two paths must stay comparable.

Returns `(summary_min, theta, z_field, rho_field, mu_field, r1, arm, imsize, global_index)`.
Storing `z_field` **and** `rho_field` is deliberate: the first is what SBC ranks (atom-free), the
second is what the Delta-rho map reports. `ArgumentError`s from the forward model propagate
(ASVS V5), never swallowed.
"""
function generate_p12_sample(idx::Integer; imsize_sampler = sample_p12_imsize,
                             arm::Symbol = :car, imsize = nothing, r1 = nothing)
    rng  = p12_datagen_rng(idx)
    draw = sample_p12_prior(rng; arm = arm, r1 = r1)
    isz  = imsize === nothing ? imsize_sampler(rng) : imsize
    # The legacy 8-field theta view PLUS the field: `simulate_pair` reads `rho_field` through one
    # defensive binding and mixes stage 1 PER PIXEL from its bilinear upsample (D-06).
    θs      = merge(theta_scalar_view(draw), (rho_field = draw.rho_field,))
    channels = simulate_pair(rng, θs; imsize = isz)
    mci     = build_mci(channels)
    return (
        summary_min  = encode_d01(patch_summary(mci)),   # RAW 128-row encoding, never z-scored
        theta        = p12_theta_column(draw),           # the ONE assembler
        z_field      = draw.z_field,                     # D-07 GROUND TRUTH (atom-free, SBC)
        rho_field    = draw.rho_field,                   # the deliverable's units (12-19)
        mu_field     = draw.mu_field,                    # the copula's intermediate, checkable
        r1           = draw.r1,
        arm          = draw.arm,
        imsize       = isz,
        global_index = Int(idx),
    )
end

"""
    generate_p12_block(indices; imsize_sampler, arm, imsize, r1, parallel) -> NamedTuple

Fill column-major buffers for `indices` (one sample = one column, the NeuralEstimators d-by-K
convention; the lattices are a `G x G x n` array rather than a vector of matrices, so downstream
slicing is cheap and type-stable). Because every column is a pure function of its global index and
each is written independently -- no shared mutable RNG, no locks -- the parallel and serial paths
are BYTE-IDENTICAL for any thread count.

PARALLEL INSIDE A SHARD, SERIAL OVER SHARDS. This is the analog's choice
(`p11_generate.jl:250-282`) and it is made for the analog's reason: at these pool sizes there are
only a handful of shards, so threading OVER shards would leave most of the available threads idle;
serial-over-shards also makes progress monotone, which is what the wall-clock guard needs.
"""
function generate_p12_block(indices; imsize_sampler = sample_p12_imsize, arm::Symbol = :car,
                            imsize = nothing, r1 = nothing,
                            parallel::Bool = Threads.nthreads() > 1)
    idxv  = collect(indices)
    n     = length(idxv)
    G     = P12_G
    theta = Matrix{Float64}(undef, P12_D_MINISPIKE, n)
    smin  = Matrix{Float64}(undef, P12_SUMMARY_MIN_DIM, n)
    zf    = Array{Float64,3}(undef, G, G, n)
    rf    = Array{Float64,3}(undef, G, G, n)
    mf    = Array{Float64,3}(undef, G, G, n)
    r1v   = Vector{Float64}(undef, n)
    armv  = Vector{Symbol}(undef, n)
    imsz  = Vector{Tuple{Int,Int}}(undef, n)
    gidx  = Vector{Int}(undef, n)

    fill_col! = function (j::Int)
        s = generate_p12_sample(idxv[j]; imsize_sampler = imsize_sampler, arm = arm,
                                imsize = imsize, r1 = r1)
        @inbounds theta[:, j]  = s.theta
        @inbounds smin[:, j]   = s.summary_min
        @inbounds zf[:, :, j]  = s.z_field
        @inbounds rf[:, :, j]  = s.rho_field
        @inbounds mf[:, :, j]  = s.mu_field
        @inbounds r1v[j]       = s.r1
        @inbounds armv[j]      = s.arm
        @inbounds imsz[j]      = s.imsize
        @inbounds gidx[j]      = s.global_index
        return nothing
    end

    if parallel
        Threads.@threads for j in 1:n
            fill_col!(j)
        end
    else
        for j in 1:n
            fill_col!(j)
        end
    end

    return (theta = theta, summary_min = smin, z_field = zf, rho_field = rf, mu_field = mf,
            r1 = r1v, arm = armv, imsize = imsz, global_index = gidx)
end

# =============================================================================================
# 5. Sharded, atomic, resume-by-skip persistence
# =============================================================================================
#
# The shard NAMING (`shard_path`), the shard SIZE, the content-hash directory resolution
# (`open_or_invalidate`) and the hash itself (`cache_hash` / `subhashes`) are reused UNCHANGED from
# `spike/data/cache.jl`. Only the integrity predicate, the writer and the manifest writer are
# Phase-12-local, because a Phase-12 shard carries the three lattices, `r1` and `arm` and carries
# NO `summary_aug` column -- `cache.jl`'s `_loads_ok` requires `summary_aug` and would therefore
# declare every finished Phase-12 shard NOT done, silently defeating resume-by-skip on a run that
# is budgeted in HOURS. The write idiom itself is the identical sequence:
# mkpath -> .tmp -> jldsave -> reopen-and-assert -> mv(...; force = true).

"""
    P12_CACHE_ROOT

The Phase-12 pool cache root. Sits under the gitignored bulk-cache tree (`.gitignore`).

**THIS IS A DIFFERENT DIRECTORY FROM `cache/p11`, AND THAT IS A HARD REQUIREMENT.** Phase 11's
50 000-pair pool is 54 MB of gitignored, unbacked data that costs roughly 56 minutes to rebuild,
has ALREADY been destroyed once by a worktree cleanup, and is still needed by 12-15 and 12-17. It
is READ-ONLY to Phase 12: nothing in this phase may write into, overwrite, move or regenerate it.
The disjointness is asserted by `test_p12_decoupling.jl` (12-05) and again by testset 3 of
`test_p12_datagen.jl`.
"""
const P12_CACHE_ROOT = joinpath(@__DIR__, "cache", "p12")

# THE DISJOINTNESS IS STRUCTURAL, NOT ASSERTED HERE, AND THAT IS 12-05's RULE RATHER THAN A
# WEAKENING OF IT. `test_p12_decoupling.jl:337` requires that the Phase-11 cache-root literal appear
# NOWHERE in this file's comment-stripped source -- so no code path here can construct that path at
# all, which is a stronger guarantee than a runtime check that the constant is not it. An earlier
# draft carried exactly such a check and FAILED that test, on its own defensive assertion. The
# positive half of the rule is asserted here; the negative half is asserted from the test files,
# where naming the other phase's root is what they are for.
@assert basename(P12_CACHE_ROOT) == "p12"

"""
    p12_shard_loads_ok(path) -> Bool

True iff `path` exists and reopens with the full Phase-12 key set present. Any open/read error (a
torn file left by a kill) is caught and reported as not-done, so the shard is regenerated rather
than silently trusted.
"""
function p12_shard_loads_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "theta") && haskey(f, "summary_min") && haskey(f, "z_field") &&
                haskey(f, "rho_field") && haskey(f, "mu_field") && haskey(f, "r1") &&
                haskey(f, "arm") && haskey(f, "imsize") && haskey(f, "global_index") &&
                haskey(f, "schema_version")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end

"Resume-by-skip predicate: a Phase-12 shard is done iff its FINAL file passes the reopen check."
p12_shard_done(path) = p12_shard_loads_ok(path)

"""
    write_p12_shard(dir, n, out) -> String

Atomically commit shard `n` from a `generate_p12_block` result: `jldsave` the column-major buffers
to a `.tmp`, reopen to integrity-check, then `mv(...; force = true)` (a filesystem-atomic rename).
A crash before the `mv` leaves only the discardable `.tmp`, never a half-written shard, so a killed
run always resumes from a consistent pool (T-12-17).
"""
function write_p12_shard(dir, n::Integer, out)
    path = shard_path(dir, n)
    tmp  = path * ".tmp"
    jldsave(tmp; theta = out.theta, summary_min = out.summary_min,
            z_field = out.z_field, rho_field = out.rho_field, mu_field = out.mu_field,
            r1 = out.r1, arm = out.arm, imsize = out.imsize,
            global_index = out.global_index, schema_version = P12_POOL_SCHEMA)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "z_field") "write_p12_shard: integrity check failed, $tmp missing z_field"
        @assert haskey(f, "r1")      "write_p12_shard: integrity check failed, $tmp missing r1"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end

"""
    p12_generating_config(n; shard_size, arm, r1, imsize_tag) -> NamedTuple

The canonical generating config the D-05 content hash names the pool directory by. Carries every
field whose change must invalidate the pool.

**`arm` AND `r1_pin` ARE IN HERE FOR TWO DIFFERENT REASONS AND BOTH ARE LOAD-BEARING.**

`arm` separates the CAR pool from the GP pool, so 12-15's controlled CAR-vs-GP comparison cannot
read one arm's pool for both (T-12-34).

`r1_pin` -- the pinned value, or `nothing` when r1 is DRAWN -- separates the five rungs of 12-11's
D-12 Stage-1 ladder. Without it, `n`, `arm` and the r1 PRIOR BOUNDS are identical across rungs, so
all five would hash to the SAME directory; resume-by-skip would then serve rung 1's data for rungs
2-5, 12-11's `n_test` assertion would pass at every rung because the COUNT is right, and the gate
that authorises the entire training spend would report a five-rung ladder that is one rung repeated
five times. "The correlation length at which recovery crosses the noise floor" -- which 12-11 calls
a scientific deliverable of this phase -- would be fabricated, silently, with every assertion green.

**`imsize_tag` IS THE ONE FIELD A CALLER MUST NOT FORGET.** The image-size SAMPLER is a closure and
cannot be hashed, so the tag is the only thing that distinguishes an F5-mixture pool from a pool
pinned to one size. A caller that pins `imsize_sampler` without passing a matching `imsize_tag`
gets a pool that hashes to the F5-mixture directory and is indistinguishable from it. The realized
counts are written into the manifest (F6) so the mistake is at least DIAGNOSABLE after the fact,
but the tag is what prevents it.
"""
p12_generating_config(n::Integer; shard_size::Integer = SHARD_SIZE, arm::Symbol = :car,
                      r1 = nothing, imsize_tag::Symbol = :f5_mixture) = (
    N               = Int(n),
    shard_size      = Int(shard_size),
    master_seed     = UInt64(P12_DEV_SEED),
    datagen_salt    = UInt64(P12_DATAGEN_SALT),
    datagen_counter = P12_DATAGEN_COUNTER,   # audit tag for this activity (R-5 ledger)
    theta_dim       = P12_D_MINISPIKE,
    G               = P12_G,
    K_dev           = P12_K_DEV,
    arm             = arm,
    r1_min          = P12_R1_MIN,            # the r1 PRIOR bounds ...
    r1_max          = P12_R1_MAX,
    r1_pin          = r1 === nothing ? nothing : Float64(r1),   # ... and the PIN, separately
    imsize_set      = P12_IMSIZE_SET,
    imsize_weights  = P12_IMSIZE_WEIGHTS,
    imsize_tag      = imsize_tag,
    summary_min_dim = P12_SUMMARY_MIN_DIM,
    schema_version  = P12_POOL_SCHEMA,
)

"""
    write_p12_meta(dir, config; imsize_counts, r1_quantiles) -> String

Persist the cache manifest into `meta.jld2` (atomically). Carries everything `cache.jl`'s
`write_meta` carries -- the D-05 `hash` that `open_or_invalidate` re-checks, the per-file
`subhashes` for which-source-changed diagnosis, the raw `config`, `N`, the shard size and the
schema tag -- PLUS the two REALIZED distributions.

THE REALIZED FIELDS ARE THE WHOLE REASON THIS WRITER IS PHASE-12-LOCAL (F6). The INTENDED prior is
not evidence of what a pool contains; only the realized values are. They cannot live in `config`,
because `config` is what the directory name is hashed from and a realized count is not known until
after the pool exists.
"""
function write_p12_meta(dir, config; imsize_counts, r1_quantiles)
    path = joinpath(dir, "meta.jld2")
    tmp  = path * ".tmp"
    jldsave(tmp;
        hash                    = cache_hash(config),
        subhashes               = subhashes(),
        config                  = config,
        N                       = config.N,
        shard_size              = config.shard_size,
        schema_version          = P12_POOL_SCHEMA,
        realized_imsize_counts  = imsize_counts,
        realized_r1_quantiles   = r1_quantiles,
    )
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "hash") "write_p12_meta: integrity check failed, $tmp missing hash"
        @assert haskey(f, "realized_r1_quantiles") "write_p12_meta: integrity check failed, " *
            "$tmp missing realized_r1_quantiles"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    p12_pool_dir(n; arm = :car, r1 = nothing, cache_root = P12_CACHE_ROOT,
                 shard_size = SHARD_SIZE, imsize_tag = :f5_mixture) -> String

Resolve (creating if absent) the content-hash-named pool directory for a pool of `n` samples,
WITHOUT generating anything. The trainer and every reported runner use this to find the pool the
datagen step wrote.

THE `arm` AND `r1` KEYWORDS ARE WRITTEN OUT IN FULL HERE ON PURPOSE, rather than deferred to the
analog with an ellipsis: `r1` is a HARD REQUIREMENT of 12-11, not an option, and a keyword a
sibling plan depends on must be VISIBLE in the signature that declares it. See
[`p12_generating_config`](@ref) for what each one prevents.
"""
p12_pool_dir(n::Integer; arm::Symbol = :car, r1 = nothing, cache_root = P12_CACHE_ROOT,
             shard_size::Integer = SHARD_SIZE, imsize_tag::Symbol = :f5_mixture) =
    open_or_invalidate(cache_root, p12_generating_config(n; shard_size = shard_size, arm = arm,
                                                         r1 = r1, imsize_tag = imsize_tag))

"""
    p12_pool_complete(dir, n; shard_size = SHARD_SIZE) -> Bool

True iff `dir` holds a COMPLETE pool of `n` samples: every one of the `cld(n, shard_size)` shards
passes the reopen integrity check, and the manifest exists.

**THIS, NOT `isdir`, IS THE PRESENCE TEST.** `p12_pool_dir` resolves the content-hash directory by
calling `open_or_invalidate`, which **CREATES the directory it resolves** (`cache.jl:162`). So
`isdir(p12_pool_dir(n; arm = arm))` is TRUE the first time it is ever asked, for a pool that does
not exist and has never been generated -- an assertion that cannot fail, and a "generate only if
absent" branch that never fires. A caller that wants "is the pool already there?" must ask this
function; a caller that wants "make sure it is there" should just call
[`generate_p12_pool`](@ref), whose resume-by-skip already costs nothing when the shards are
complete.
"""
function p12_pool_complete(dir, n::Integer; shard_size::Integer = SHARD_SIZE)
    isdir(dir) || return false
    isfile(joinpath(dir, "meta.jld2")) || return false
    return all(s -> p12_shard_done(shard_path(dir, s)), 1:cld(n, shard_size))
end

"""
    generate_p12_pool(n; shard_size = SHARD_SIZE, imsize_sampler = sample_p12_imsize,
                      arm = :car, r1 = nothing, cache_root = P12_CACHE_ROOT,
                      imsize_tag = :f5_mixture, parallel, wallclock_ceiling_min,
                      verbose = true) -> String

Generate the Phase-12 field-aware training pool of `n` samples into a sharded, atomic,
content-hash-named cache, and return the directory.

Shards are processed SERIALLY while the samples INSIDE a shard are filled in parallel -- see
[`generate_p12_block`](@ref) for why.

Already-finished shards are SKIPPED (`p12_shard_done`), so a killed run RESUMES rather than
restarting. Because every column is a pure function of its global index, a resumed pool is
BYTE-IDENTICAL to one generated in a single run.

THE WALL-CLOCK GUARD IS A **BLOCKER, NOT A DOWNGRADE TRIGGER.** After every shard the run projects
its own total against `P12_DATAGEN_WALLCLOCK_CEILING_MIN` and THROWS on projected exceedance. It
does not, and must not, shrink the image-size arm to fit the budget.

`r1` IS PASSED TO BOTH THE DRAW AND THE CONTENT HASH. When it is pinned, every realized `r1` in the
written pool is asserted equal to it before the manifest is written -- so a pin that is ACCEPTED but
not APPLIED fails at generation rather than at interpretation, hours later.
"""
function generate_p12_pool(n::Integer; shard_size::Integer = SHARD_SIZE,
                           imsize_sampler = sample_p12_imsize,
                           arm::Symbol = :car, r1 = nothing,
                           cache_root = P12_CACHE_ROOT,
                           imsize_tag::Symbol = :f5_mixture,
                           parallel::Bool = Threads.nthreads() > 1,
                           wallclock_ceiling_min::Real = P12_DATAGEN_WALLCLOCK_CEILING_MIN,
                           verbose::Bool = true)
    config  = p12_generating_config(n; shard_size = shard_size, arm = arm, r1 = r1,
                                    imsize_tag = imsize_tag)
    dir     = open_or_invalidate(cache_root, config)
    nshards = cld(n, shard_size)

    if verbose
        println("="^78)
        println("Phase-12 field-aware training pool (D-05, D-07, D-08)")
        println("  N = $n   shards = $nshards x $shard_size   threads = $(Threads.nthreads())")
        println("  arm = $arm   r1 " *
                (r1 === nothing ? "~ Uniform($P12_R1_MIN, $P12_R1_MAX) (DRAWN)" : "= $r1 (PINNED)"))
        println("  field -> copula -> rho; theta = $(P12_D_MINISPIKE) rows in P12_THETA_ROWS order")
        println("  imsize arm = $imsize_tag   set = $(P12_IMSIZE_SET)")
        println("  stream = P12_DEV_SEED $(repr(UInt64(P12_DEV_SEED))) xor " *
                "P12_DATAGEN_SALT $(repr(P12_DATAGEN_SALT)), keyed per global index")
        println("  wall-clock ceiling = $wallclock_ceiling_min min (BLOCKER, not a downgrade trigger)")
        println("  dir = $dir")
        println("="^78)
    end

    t0        = time()
    n_done    = 0                    # shards this RUN actually generated (resumed ones excluded)
    n_todo    = count(s -> !p12_shard_done(shard_path(dir, s)), 1:nshards)
    n_skipped = nshards - n_todo
    verbose && n_skipped > 0 &&
        println("[resume] $n_skipped of $nshards shards already complete — skipping them")

    for s in 1:nshards
        path = shard_path(dir, s)
        if p12_shard_done(path)
            continue                                   # resume-by-skip
        end
        ts  = time()
        lo  = (s - 1) * shard_size + 1
        hi  = min(s * shard_size, n)
        out = generate_p12_block(lo:hi; imsize_sampler = imsize_sampler, arm = arm, r1 = r1,
                                 parallel = parallel)
        write_p12_shard(dir, s, out)
        n_done += 1

        elapsed_min = (time() - t0) / 60
        shard_min   = (time() - ts) / 60
        verbose && println("[shard $s/$nshards] $(hi - lo + 1) samples in " *
                           "$(round(shard_min; digits = 2)) min " *
                           "(cumulative $(round(elapsed_min; digits = 2)) min)")

        # --- The pre-registered wall-clock BLOCKER. Projected from the work THIS RUN did, so a
        #     resumed run is judged on its own remaining work rather than on a stale total. The
        #     projection includes first-shard compile overhead, which makes it conservative -- it
        #     fires early rather than late, which is the safe direction for a blocker.
        projected_min = elapsed_min * (n_todo / n_done)
        if projected_min > wallclock_ceiling_min
            println("!"^78)
            println("PHASE-12 DATAGEN WALL-CLOCK BLOCKER")
            println("!"^78)
            throw(ErrorException(
                "p12 datagen projected total $(round(projected_min; digits = 1)) min exceeds the " *
                "pre-registered ceiling of $wallclock_ceiling_min min " *
                "($(n_done)/$(n_todo) shards in $(round(elapsed_min; digits = 2)) min). " *
                "RECORD THIS AS A BLOCKER, naming the measured rate and the shortfall. " *
                "P12_DATAGEN_WALLCLOCK_CEILING_MIN is APPEND-ONLY Tier-1 and may NOT be raised. " *
                "SILENTLY DOWNGRADING THE IMAGE-SIZE ARM TO A SINGLE CHEAP SIZE IS FORBIDDEN: " *
                "that is the exact 07-05 compute-budget artifact the F5 amendment exists to " *
                "correct, and it would break the F5 binding invariant that the joint an estimator " *
                "is TRAINED on must equal the joint it is EVALUATED on. The pool is sharded and " *
                "resumable — re-launching continues from the last completed shard."))
        end
    end

    # --- The realized record (F6), and the pin check -------------------------------------------
    # Read back from DISK rather than from the in-memory buffers: on a resumed run most shards were
    # never in memory at all, and it is the WRITTEN pool the assertion has to be about.
    files = _p12_shard_files(dir)
    r1_realized = reduce(vcat, [JLD2.load(f, "r1")     for f in files])
    imsz_realized = reduce(vcat, [JLD2.load(f, "imsize") for f in files])
    if r1 !== nothing
        pin = Float64(r1)
        bad = count(!=(pin), r1_realized)
        @assert bad == 0 "generate_p12_pool: r1 was PINNED to $pin but $bad of " *
            "$(length(r1_realized)) realized draws differ (range " *
            "$(extrema(r1_realized))). A pin that is accepted and not applied would give 12-11 a " *
            "five-rung ladder that is one rung repeated five times."
    end
    imsize_counts = realized_imsize_counts(imsz_realized)
    r1_quantiles  = realized_r1_quantiles(r1_realized)

    write_p12_meta(dir, config; imsize_counts = imsize_counts, r1_quantiles = r1_quantiles)
    total_min = (time() - t0) / 60
    verbose && println("pool complete in $(round(total_min; digits = 2)) min " *
                       "($n_done generated, $n_skipped resumed) -> $dir")
    return dir
end

"Shard files of `dir` in lexical (== numeric, zero-padded) order."
_p12_shard_files(dir) =
    sort(filter(f -> startswith(basename(f), "shard_") && endswith(f, ".jld2"),
                readdir(dir; join = true)))

"""
    load_p12_pool(dir) -> NamedTuple

Load the RAW (never z-scored) Phase-12 pool by concatenating every shard in `dir` in lexical
(== numeric, zero-padded) order. Returns

    (theta = D×N, summary_min = 128×N, z_field = G×G×N, rho_field = G×G×N, mu_field = G×G×N,
     r1 = N, arm = N, imsize = N, global_index = N)

SINGLE-STACK, per the file header: the nine fields above are ALL of them. There is no paired
control stack and no stored per-region difference field. A reader expecting either should read the
header's explanation of why Delta-rho is built at READ time, from two independent passes of the
trained net, rather than stored here.

The lattices come back as a 3-D array rather than a vector of matrices so downstream slicing is
cheap and type-stable.
"""
function load_p12_pool(dir)
    files = _p12_shard_files(dir)
    isempty(files) && error("load_p12_pool: no shard_*.jld2 found in $dir")
    theta        = reduce(hcat, [JLD2.load(f, "theta")        for f in files])
    summary_min  = reduce(hcat, [JLD2.load(f, "summary_min")  for f in files])
    z_field      = cat([JLD2.load(f, "z_field")   for f in files]...; dims = 3)
    rho_field    = cat([JLD2.load(f, "rho_field") for f in files]...; dims = 3)
    mu_field     = cat([JLD2.load(f, "mu_field")  for f in files]...; dims = 3)
    r1           = reduce(vcat, [JLD2.load(f, "r1")           for f in files])
    arm          = reduce(vcat, [JLD2.load(f, "arm")          for f in files])
    imsize       = reduce(vcat, [JLD2.load(f, "imsize")       for f in files])
    global_index = reduce(vcat, [JLD2.load(f, "global_index") for f in files])
    return (theta = theta, summary_min = summary_min, z_field = z_field, rho_field = rho_field,
            mu_field = mu_field, r1 = r1, arm = arm, imsize = imsize, global_index = global_index)
end

# =============================================================================================
# 6. The realized record (F6)
# =============================================================================================
#
# THE INTENDED PRIOR IS NOT EVIDENCE OF WHAT A POOL CONTAINS, ONLY THE REALIZED VALUES ARE. F6
# established an unrecorded training image-size distribution as a provenance DEFECT; the same
# argument applies verbatim to the r1 prior, which is why both are reported and both are written
# into the manifest.

# `realized_imsize_counts` is DEFINED UNDER A GUARD because `spike/data/p11_generate.jl:503-509`
# defines the same one-argument function, with identical semantics, and the two files are both
# included into `Main` during a full-suite run. The guard makes "whichever loaded last wins"
# unnecessary rather than merely harmless -- a phase-namesake silently overwriting another phase's
# method by an identical zero-keyword signature is a defect this phase has already found once.
if !isdefined(@__MODULE__, :realized_imsize_counts)
    # The REALIZED image-size distribution of a generated pool, as counts per size.
    function realized_imsize_counts(imsize)
        counts = Dict{Tuple{Int,Int},Int}()
        for z in imsize
            counts[z] = get(counts, z, 0) + 1
        end
        return counts
    end
end

"""
    realized_r1_quantiles(r1; probs = (0.0, 0.25, 0.5, 0.75, 1.0)) -> NamedTuple

The REALIZED correlation-length distribution of a generated pool, as the quantiles at `probs`.
Returns `(probs = collect(probs), values = ...)`, so the monotonicity of `values` is directly
checkable and the pair round-trips through JLD2 unambiguously.

On a PINNED pool every entry is the pin, so the vector is constant -- `issorted` still holds, and
that degenerate shape is itself the evidence that the pin reached the draw.
"""
function realized_r1_quantiles(r1; probs = (0.0, 0.25, 0.5, 0.75, 1.0))
    v = collect(float.(r1))
    isempty(v) && throw(ArgumentError("realized_r1_quantiles: empty r1 vector"))
    return (probs = collect(float.(probs)), values = [quantile(v, p) for p in probs])
end
