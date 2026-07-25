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

# spike/data/p11_generate.jl --- Phase-11 lambda-hierarchical training-pool generation (D-03, D-09).
#
# THE GENERATIVE STORY, IN EXACTLY THIS ORDER. FOR EVERY TRAINING SAMPLE i:
#     lambda_i        ~ Uniform(LAMBDA_MIN, LAMBDA_MAX)   <- THE UNCERTAINTY LEVEL, DRAWN **FIRST**
#     shift_dx_i      ~ Uniform(-lambda_i, +lambda_i)     <- DRAWN **CONDITIONAL ON** lambda_i
#     shift_dy_i      ~ Uniform(-lambda_i, +lambda_i)     <- DRAWN **CONDITIONAL ON** lambda_i
#     chromatic_eps_i ~ CHROMATIC_PRIOR                   <- D-09: ITS OWN FIXED PRIOR, NOT lambda-SCALED
#     rho_true_i      = ghat(mu*_i), mu*_i ~ MU_PRIOR     <- UNCHANGED
#     the other four nuisances ~ their existing priors    <- UNCHANGED
#     imsize_i        ~ THE F5 CATEGORICAL                <- READ FROM THE FROZEN PRE-REGISTRATION
#     y_i             = simulate_pair(rng_i, theta_i; imsize = imsize_i)
#
# THIS ORDER IS THE ENTIRE POINT OF THE FILE. IT EXISTS TO PREVENT **PITFALL 1 / RISK R4**: IF
# lambda WERE DRAWN INDEPENDENTLY OF THE SHIFT (OR AFTER IT), lambda WOULD CARRY NO INFORMATION
# ABOUT THE SHIFT, THE NETWORK WOULD **CORRECTLY** LEARN TO IGNORE THE 129th INPUT ROW, AND SC2
# WOULD FAIL FLAT FOR A REASON THAT HAS NOTHING TO DO WITH RESOLUTION -- A NULL THAT LOOKS EXACTLY
# LIKE "BELOW THE FIXED 8x8 SUMMARY'S RESOLUTION" WHILE ACTUALLY BEING AN IMPLEMENTATION BUG. THOSE
# TWO READINGS HAVE COMPLETELY DIFFERENT CONSEQUENCES. THE PER-SAMPLE INVARIANT ASSERTION IN
# `sample_p11_theta` AND THE BEHAVIOURAL TRIPWIRE `spike/test/test_lambda_ablation.jl` ARE THE TWO
# THINGS THAT TELL THEM APART; NEITHER MAY BE WEAKENED.
#
# lambda IS **NOT** A THETA ROW. It is a conditioning INPUT (D-03). theta stays 8 rows in
# `sample_prior` field order; lambda rides alongside in its own column-parallel vector.
#
# WHAT THIS FILE DELIBERATELY DOES **NOT** DO. It stores the RAW 128-row summary and the raw
# lambda separately. It does not z-score anything and it does not build the 129-row input.
# Fitting the frozen summary transform is the loader's train-only job (leak-freedom), and the
# lambda row is appended only AFTER that fit has been applied -- see the ordering contract in
# `spike/npe/p11_architecture.jl`. Doing either here would either leak validation statistics into
# the fit or z-score the user-declared conditioning variable against the training pool.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local generation; reaches src/ only READ-ONLY,
# transitively through `contract.jl`'s include and through the content-hash guard's `read` of
# frozen source bytes. Nothing here writes to, or imports from, the ProteinCoLoc module.
# Guarded includes keep the file loadable standalone AND inside a test harness.

# --- ORDER MATTERS: the Tier-1/Tier-2 pre-registration FIRST (it binds NPE_MASTER_SEED,
#     VAL_MASTER_SEED, VAL_FIX_SEED and DEFAULT_MASTER_SEED at UInt64 width; a later file that
#     rebinds one of them at a narrower integer width is a hard `const` error), then the frozen
#     summary contract, the prior, the forward model, the seeding primitives whose salt idiom
#     `p11_datagen_rng` mirrors, the encoder, and the sharded-cache persistence layer.
isdefined(@__MODULE__, :LAMBDA_MIN)    || include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :build_mci)     || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :sample_prior)  || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair) || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :HOLDOUT_SALT)  || include(joinpath(@__DIR__, "seeding.jl"))
isdefined(@__MODULE__, :encode_d01)    || include(joinpath(@__DIR__, "encode.jl"))
isdefined(@__MODULE__, :write_shard)   || include(joinpath(@__DIR__, "cache.jl"))

# =============================================================================================
# 1. The reserved Phase-11 datagen stream (D-01)
# =============================================================================================

"""
    P11_DATAGEN_SALT

The fresh XOR salt that carves the Phase-11 TRAINING-POOL draw stream. A distinct salt gives a
distinct Philox KEY WORD, hence a provably disjoint stream family.

DISTINCT FROM, and asserted below against, every other salt in the repository:
`HOLDOUT_SALT` (0x9E37…, the reserved ADVI holdout), `FOLD_SALT` (0xD1B5…, the loader's k-fold
permutation), `VAL_SALT` (0xBF58…, the Phase-5 validation stream -- written as a literal here
because `spike/validation/consts.jl` cannot be included alongside `p11_consts.jl`), `PROD_SALT`
and `AMEND_SALT` (the two frozen ship-gate families), and `P11_SALT` (this phase's own
`p11_rng` family, which the probe and every evaluation counter ride).

Value is the xorshift64* multiplier -- good avalanche, well away from zero; only
nonzero-and-distinct is load-bearing.
"""
const P11_DATAGEN_SALT = 0x2545_F491_4F6C_DD1D

# The Phase-11 datagen stream is carved by its SALT, not by a reserved counter, because the
# counter word must carry the per-sample GLOBAL INDEX (that is what makes the pool byte-identical
# for any thread count). `P11_DATAGEN_COUNTER` is therefore recorded in the pool's generating
# config as the audit-trail tag for this activity rather than consumed as a key word; the
# disjointness it exists to guarantee is delivered here by the salt and asserted below.
@assert P11_DATAGEN_SALT ∉ (HOLDOUT_SALT,                # reserved ADVI holdout stream
                            FOLD_SALT,                   # loader k-fold permutation stream
                            0xBF58_476D_1CE4_E5B9,       # VAL_SALT (spike/validation/consts.jl)
                            PROD_SALT, AMEND_SALT,       # the two frozen ship-gate families
                            P11_SALT)                    # this phase's own p11_rng family
@assert (UInt64(P11_DEV_SEED) ⊻ P11_DATAGEN_SALT) != (UInt64(P11_DEV_SEED) ⊻ P11_SALT)

"""
    p11_datagen_rng(idx::Integer) -> Philox4x

The per-sample keyed RNG of the Phase-11 training pool: `Philox4x((P11_DEV_SEED ⊻
P11_DATAGEN_SALT, idx))`. Counter-based and keyed by the GLOBAL INDEX, so the draw stream is a
pure function of the index -- the pool is BYTE-IDENTICAL regardless of thread count, execution
order, or how many times a killed run is resumed. Mirrors `sample_rng`/`holdout_rng`
(`spike/data/seeding.jl:75-87`).
"""
p11_datagen_rng(idx::Integer) =
    Philox4x(UInt64, (UInt64(P11_DEV_SEED) ⊻ P11_DATAGEN_SALT, UInt64(idx)))

# =============================================================================================
# 2. The hierarchical sampler -- lambda FIRST, shift CONDITIONAL on lambda
# =============================================================================================

# theta arity, DERIVED from the prior rather than retyped (the repo-wide derive-never-retype
# convention). A future theta extension cannot silently desynchronise the buffer width.
const P11_THETA_DIM = length(sample_prior(Random.Xoshiro(0)))

# On-disk schema tag for a Phase-11 pool shard. Bump iff the shard key set changes.
const P11_POOL_SCHEMA = 1

"""
    sample_p11_theta(rng::AbstractRNG) -> (theta::NamedTuple, lambda::Float64)

Draw one Phase-11 training sample's parameters under the lambda-hierarchical joint
`pi(lambda) * pi(theta | lambda)`.

Draw order (load-bearing, see the file header):

  1. `lambda ~ Uniform(LAMBDA_MIN, LAMBDA_MAX)` -- the registration-uncertainty LEVEL, FIRST;
  2. `shift_dx, shift_dy ~ Uniform(-lambda, +lambda)` -- CONDITIONAL on that draw;
  3. every remaining field from its own existing prior object, exactly as `sample_prior` does,
     including `chromatic_eps` from `CHROMATIC_PRIOR` (D-09 / ruling Q3: a FIXED marginalized
     nuisance, deliberately NOT lambda-scaled and carrying no conditioning input of its own).

The returned NamedTuple is assembled in the SAME field order `sample_prior` uses, so
`collect(values(theta))` stays row-aligned with `p11_theta_prior_bounds()` and with every
positional theta read downstream.

THE MARGINAL-VERSUS-CONDITIONAL CONSEQUENCE, WHICH THE REPORT MUST STATE. INTEGRATING lambda
OUT, THE **MARGINAL** PRIOR ON `shift_dx` IS **NOT** `Uniform(-3, 3)`: IT IS
`p(s) = INTEGRAL 1/(2*lambda) * 1[|s| <= lambda] dpi(lambda)`, A PEAKED, HEAVIER-AT-ZERO DENSITY.
THE SUPPORT BOX `[-3, 3]` IS UNCHANGED AND STAYS CORRECT, BUT ANY RANK, COVERAGE OR SHRINKAGE
STATISTIC ON THE SHIFT COLUMNS MUST BE COMPUTED **CONDITIONAL ON lambda** -- I.E. POOLED WITHIN A
RUNG -- NEVER AGAINST A `Uniform(-3, 3)` REFERENCE. (D-05 already declares those columns vacuous
and non-evidential, so the practical rule is: report them per-rung, labelled non-evidence.)
"""
function sample_p11_theta(rng::AbstractRNG)
    # (1) the UNCERTAINTY LEVEL, drawn FIRST.
    lambda = rand(rng, Uniform(LAMBDA_MIN, LAMBDA_MAX))
    # (2) the registration error, drawn CONDITIONAL on it. `Uniform(-lambda, lambda)` (the
    #     half-width parameterisation D-03 names verbatim), never a Gaussian of scale lambda:
    #     it keeps lambda = LAMBDA_MAX exactly equal to the stated SHIFT_PRIOR, so the widest
    #     rung IS the prior and the ladder never leaves it.
    shift_dx = rand(rng, Uniform(-lambda, lambda))
    shift_dy = rand(rng, Uniform(-lambda, lambda))
    # (3) everything else from its own unchanged prior object.
    μ_star           = rand(rng, MU_PRIOR)
    spillover        = rand(rng, SPILLOVER_PRIOR)
    autofluorescence = rand(rng, AUTOFLUORESCENCE_PRIOR)
    label_efficiency = rand(rng, LABEL_EFFICIENCY_PRIOR)
    noise            = rand(rng, NOISE_PRIOR)
    chromatic_eps    = rand(rng, CHROMATIC_PRIOR)   # D-09 / Q3: fixed prior, NOT lambda-scaled

    theta = (
        ρ_true           = ghat(μ_star),
        spillover        = spillover,
        autofluorescence = autofluorescence,
        label_efficiency = label_efficiency,
        shift_dx         = shift_dx,
        shift_dy         = shift_dy,
        noise            = noise,
        chromatic_eps    = chromatic_eps,
    )

    # THE CHEAP STRUCTURAL GUARD AGAINST PITFALL 1 / R4. If a future edit ever reorders the draws
    # so the shift stops being conditional on lambda, this fires on the FIRST sample rather than
    # surfacing hours later as an unexplained null SC2. It is asserted per sample on purpose.
    @assert abs(theta.shift_dx) <= lambda && abs(theta.shift_dy) <= lambda "sample_p11_theta: " *
        "|shift| exceeds lambda ($(theta.shift_dx), $(theta.shift_dy) vs $lambda) — the shift is " *
        "no longer drawn CONDITIONAL on lambda (Pitfall 1 / R4)"

    return (theta, lambda)
end

"""
    sample_p11_imsize(rng::AbstractRNG) -> Tuple{Int,Int}

Draw one image size from the F5 mixture, `P11_IMSIZE_SET` under `P11_IMSIZE_WEIGHTS`.

BOTH ARE **READ** FROM `spike/validation/p11_consts.jl`, WHICH ITSELF READS THEM FROM THE FROZEN
AMENDED PRE-REGISTRATION THROUGH THE ISOLATED `GateV2` MODULE. THE TUPLE IS NEVER RETYPED HERE --
a retyped mixture is a mixture that can silently drift, and F5's binding invariant is
train-joint == eval-joint. Hand-rolled inverse-CDF over the fixed categorical, mirroring
`sample_imsize` (`spike/data/seeding.jl:106-114`); no extra dependency.
"""
function sample_p11_imsize(rng::AbstractRNG)
    r = rand(rng)
    c = 0.0
    @inbounds for i in eachindex(P11_IMSIZE_SET)
        c += P11_IMSIZE_WEIGHTS[i]
        r <= c && return P11_IMSIZE_SET[i]
    end
    return P11_IMSIZE_SET[end]   # float-rounding guard: r just over the weight sum -> last bin
end

"""
    generate_p11_sample(idx::Integer; imsize_sampler = sample_p11_imsize) -> NamedTuple

Generate ONE pool row, keyed entirely by its global index `idx`:

    rng    = p11_datagen_rng(idx)
    theta, lambda = sample_p11_theta(rng)        # lambda FIRST, shift conditional
    isz    = imsize_sampler(rng)                 # the F5 categorical
    s_min  = encode_d01(patch_summary(build_mci(simulate_pair(rng, theta; imsize = isz))))

Returns `(theta::Vector, lambda::Float64, s_min::Vector, idx::Int, imsize)`. `theta` is
`collect(values(...))` in prior field order; `s_min` is the RAW (un-z-scored) 128-row encoding.
`simulate_pair`'s `ArgumentError`s propagate (ASVS V5), never swallowed.
"""
function generate_p11_sample(idx::Integer; imsize_sampler = sample_p11_imsize)
    rng            = p11_datagen_rng(idx)
    theta, lambda  = sample_p11_theta(rng)
    isz            = imsize_sampler(rng)
    mci            = build_mci(simulate_pair(rng, theta; imsize = isz))
    return (
        theta  = collect(values(theta)),   # P11_THETA_DIM-vector, prior field order
        lambda = lambda,                   # NOT a theta row -- a conditioning level (D-03)
        s_min  = encode_d01(patch_summary(mci)),   # RAW 128-row encoding, un-standardized
        idx    = Int(idx),
        imsize = isz,
    )
end

"""
    generate_p11_block(indices; imsize_sampler, parallel) -> NamedTuple

Fill column-major buffers for `indices` (one sample = one column, the NeuralEstimators d-by-K
convention). Because every column is a pure function of its global index and each is written
independently -- no shared mutable RNG, no locks -- the parallel and serial paths are
BYTE-IDENTICAL for any thread count.

The theta-buffer width is `P11_THETA_DIM` (derived from the prior), never a literal.
"""
function generate_p11_block(indices; imsize_sampler = sample_p11_imsize,
                            parallel::Bool = Threads.nthreads() > 1)
    idxv   = collect(indices)
    n      = length(idxv)
    theta  = Matrix{Float64}(undef, P11_THETA_DIM, n)
    lambda = Vector{Float64}(undef, n)
    s_min  = Matrix{Float64}(undef, 128, n)
    gidx   = Vector{Int}(undef, n)
    imsz   = Vector{Tuple{Int,Int}}(undef, n)

    fill_col! = function (j::Int)
        s = generate_p11_sample(idxv[j]; imsize_sampler = imsize_sampler)
        @inbounds theta[:, j] = s.theta
        @inbounds lambda[j]   = s.lambda
        @inbounds s_min[:, j] = s.s_min
        @inbounds gidx[j]     = s.idx
        @inbounds imsz[j]     = s.imsize
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

    return (theta = theta, lambda = lambda, summary_min = s_min,
            global_index = gidx, imsize = imsz)
end

# =============================================================================================
# 3. Sharded, atomic, resume-by-skip persistence
# =============================================================================================
#
# The shard NAMING (`shard_path`), the shard SIZE, the content-hash directory resolution
# (`open_or_invalidate`) and the manifest write (`write_meta`) are reused UNCHANGED from
# `spike/data/cache.jl`. Only the integrity predicate and the writer are Phase-11-local, because
# a Phase-11 shard carries a `lambda` column and carries NO `summary_aug` column -- `cache.jl`'s
# `_loads_ok` requires `summary_aug` and would therefore declare every finished Phase-11 shard
# NOT done, silently defeating resume-by-skip on a 1.2 h run. The write idiom itself is the
# identical S4 sequence: mkpath -> .tmp -> jldsave -> reopen-and-assert -> mv(...; force = true).

"The Phase-11 pool cache root. Sits under the gitignored bulk-cache tree (`.gitignore`)."
const P11_CACHE_ROOT = joinpath(@__DIR__, "cache", "p11")

"""
    p11_shard_loads_ok(path) -> Bool

True iff `path` exists and reopens with the full Phase-11 key set present. Any open/read error
(a torn file left by a kill) is caught and reported as not-done, so the shard is regenerated
rather than silently trusted.
"""
function p11_shard_loads_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "theta") && haskey(f, "lambda") && haskey(f, "summary_min") &&
                haskey(f, "global_index") && haskey(f, "imsize") && haskey(f, "schema_version")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end

"Resume-by-skip predicate: a Phase-11 shard is done iff its FINAL file passes the reopen check."
p11_shard_done(path) = p11_shard_loads_ok(path)

"""
    write_p11_shard(dir, n, out) -> String

Atomically commit shard `n` from a `generate_p11_block` result: `jldsave` the column-major
buffers to a `.tmp`, reopen to integrity-check, then `mv(...; force = true)` (a
filesystem-atomic rename). A crash before the `mv` leaves only the discardable `.tmp`, never a
half-written shard, so a killed run always resumes from a consistent pool (T-11-30).
"""
function write_p11_shard(dir, n::Integer, out)
    path = shard_path(dir, n)
    tmp  = path * ".tmp"
    jldsave(tmp; theta = out.theta, lambda = out.lambda, summary_min = out.summary_min,
            global_index = out.global_index, imsize = out.imsize,
            schema_version = P11_POOL_SCHEMA)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "lambda") "write_p11_shard: integrity check failed, $tmp missing lambda"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end

"""
    p11_generating_config(n; shard_size, imsize_tag) -> NamedTuple

The canonical generating config the D-05 content hash names the pool directory by. Carries every
field whose change must invalidate the pool: the size, the shard size, the DEV seed and datagen
salt, the reserved-counter audit tag, the derived theta arity, the lambda range (the
hierarchical joint itself), the image-size mixture, the summary dimension and the schema tag.
"""
p11_generating_config(n::Integer; shard_size::Integer = SHARD_SIZE,
                      imsize_tag::Symbol = :f5_mixture) = (
    N               = Int(n),
    shard_size      = Int(shard_size),
    master_seed     = UInt64(P11_DEV_SEED),
    datagen_salt    = UInt64(P11_DATAGEN_SALT),
    datagen_counter = P11_DATAGEN_COUNTER,   # audit tag for this activity (D-01 ledger)
    theta_dim       = P11_THETA_DIM,
    lambda_min      = LAMBDA_MIN,
    lambda_max      = LAMBDA_MAX,
    imsize_set      = P11_IMSIZE_SET,
    imsize_weights  = P11_IMSIZE_WEIGHTS,
    imsize_tag      = imsize_tag,
    summary_min_dim = 128,
    schema_version  = P11_POOL_SCHEMA,
)

"""
    p11_pool_dir(n; cache_root, shard_size, imsize_tag) -> String

Resolve (creating if absent) the content-hash-named pool directory for a pool of `n` samples,
without generating anything. The trainer uses this to find the pool the datagen step wrote.
"""
p11_pool_dir(n::Integer; cache_root = P11_CACHE_ROOT, shard_size::Integer = SHARD_SIZE,
             imsize_tag::Symbol = :f5_mixture) =
    open_or_invalidate(cache_root, p11_generating_config(n; shard_size = shard_size,
                                                         imsize_tag = imsize_tag))

"""
    generate_p11_pool(n; cache_root, shard_size, imsize_sampler, imsize_tag,
                      parallel, wallclock_ceiling_min, verbose) -> String

Generate the Phase-11 lambda-hierarchical training pool of `n` samples into a sharded, atomic,
content-hash-named cache, and return the directory.

Shards are processed SERIALLY while the samples INSIDE a shard are filled in parallel. This is
deliberate and differs from `generate_cache`, which threads over shards: at the pre-registered
`P11_N_PAIRS` there are only five shards, so threading over shards would leave 27 of the 32
available threads idle and turn a ~1.2 h run into a multi-hour one. Serial-over-shards also
makes progress monotone, which is what the wall-clock guard below needs.

Already-finished shards are SKIPPED (`p11_shard_done`), so a killed run RESUMES rather than
restarting -- the R12 mitigation.

THE WALL-CLOCK GUARD IS A BLOCKER, NOT A DOWNGRADE TRIGGER. After every shard the run projects
its own total against `P11_DATAGEN_WALLCLOCK_CEILING_MIN` and THROWS on projected exceedance.
It does not, and must not, shrink the image-size arm to fit the budget.
"""
function generate_p11_pool(n::Integer; cache_root = P11_CACHE_ROOT,
                           shard_size::Integer = SHARD_SIZE,
                           imsize_sampler = sample_p11_imsize,
                           imsize_tag::Symbol = :f5_mixture,
                           parallel::Bool = Threads.nthreads() > 1,
                           wallclock_ceiling_min::Real = P11_DATAGEN_WALLCLOCK_CEILING_MIN,
                           verbose::Bool = true)
    config  = p11_generating_config(n; shard_size = shard_size, imsize_tag = imsize_tag)
    dir     = open_or_invalidate(cache_root, config)
    nshards = cld(n, shard_size)

    if verbose
        println("="^78)
        println("Phase-11 lambda-hierarchical training pool (D-03, D-09)")
        println("  N = $n   shards = $nshards x $shard_size   threads = $(Threads.nthreads())")
        println("  lambda ~ Uniform($LAMBDA_MIN, $LAMBDA_MAX); shift ~ Uniform(-lambda, lambda)")
        println("  imsize arm = $imsize_tag   set = $(P11_IMSIZE_SET)")
        println("  stream = P11_DEV_SEED $(repr(UInt64(P11_DEV_SEED))) xor " *
                "P11_DATAGEN_SALT $(repr(P11_DATAGEN_SALT)), keyed per global index")
        println("  wall-clock ceiling = $wallclock_ceiling_min min (BLOCKER, not a downgrade trigger)")
        println("  dir = $dir")
        println("="^78)
    end

    t0        = time()
    n_done    = 0                    # shards this RUN actually generated (resumed ones excluded)
    n_todo    = count(s -> !p11_shard_done(shard_path(dir, s)), 1:nshards)
    n_skipped = nshards - n_todo
    verbose && n_skipped > 0 &&
        println("[resume] $n_skipped of $nshards shards already complete — skipping them")

    for s in 1:nshards
        path = shard_path(dir, s)
        if p11_shard_done(path)
            continue                                   # resume-by-skip (R12)
        end
        ts  = time()
        lo  = (s - 1) * shard_size + 1
        hi  = min(s * shard_size, n)
        out = generate_p11_block(lo:hi; imsize_sampler = imsize_sampler, parallel = parallel)
        write_p11_shard(dir, s, out)
        n_done += 1

        elapsed_min = (time() - t0) / 60
        shard_min   = (time() - ts) / 60
        verbose && println("[shard $s/$nshards] $(hi - lo + 1) samples in " *
                           "$(round(shard_min; digits = 2)) min " *
                           "(cumulative $(round(elapsed_min; digits = 2)) min)")

        # --- The pre-registered wall-clock BLOCKER. Projected from the work THIS RUN did, so a
        #     resumed run is judged on its own remaining work rather than on a stale total. The
        #     projection includes first-shard compile overhead, which makes it conservative --
        #     it fires early rather than late, which is the safe direction for a blocker.
        projected_min = elapsed_min * (n_todo / n_done)
        if projected_min > wallclock_ceiling_min
            println("!"^78)
            println("PHASE-11 DATAGEN WALL-CLOCK BLOCKER")
            println("!"^78)
            throw(ErrorException(
                "p11 datagen projected total $(round(projected_min; digits = 1)) min exceeds the " *
                "pre-registered ceiling of $wallclock_ceiling_min min " *
                "($(n_done)/$(n_todo) shards in $(round(elapsed_min; digits = 2)) min). " *
                "RECORD THIS AS A BLOCKER in .planning/STATE.md, naming the measured rate and the " *
                "shortfall. SILENTLY DOWNGRADING THE IMAGE-SIZE ARM TO A SINGLE CHEAP SIZE IS " *
                "FORBIDDEN: that is the exact 07-05 compute-budget artifact the F5 amendment " *
                "exists to correct, and it would leave the D-17 real-image check with no basis " *
                "for transfer. The pool is sharded and resumable — re-launching continues."))
        end
    end

    write_meta(dir, config)
    total_min = (time() - t0) / 60
    verbose && println("pool complete in $(round(total_min; digits = 2)) min " *
                       "($n_done generated, $n_skipped resumed) -> $dir")
    return dir
end

"""
    load_p11_pool(dir) -> NamedTuple

Load the RAW (un-standardized) Phase-11 pool by `hcat`-ing every shard in `dir` in lexical
(== numeric, zero-padded) order. Returns
`(theta = P11_THETA_DIM x N, lambda = N, summary_min = 128 x N, global_index = N, imsize = N)`.
"""
function load_p11_pool(dir)
    files = sort(filter(f -> startswith(basename(f), "shard_") && endswith(f, ".jld2"),
                        readdir(dir; join = true)))
    isempty(files) && error("load_p11_pool: no shard_*.jld2 found in $dir")
    theta        = reduce(hcat, [JLD2.load(f, "theta")        for f in files])
    summary_min  = reduce(hcat, [JLD2.load(f, "summary_min")  for f in files])
    lambda       = reduce(vcat, [JLD2.load(f, "lambda")       for f in files])
    global_index = reduce(vcat, [JLD2.load(f, "global_index") for f in files])
    imsize       = reduce(vcat, [JLD2.load(f, "imsize")       for f in files])
    return (theta = theta, lambda = lambda, summary_min = summary_min,
            global_index = global_index, imsize = imsize)
end

"""
    realized_imsize_counts(imsize) -> Dict{Tuple{Int,Int},Int}

The REALIZED image-size distribution of a generated pool, as counts per size. F6 established
that an unrecorded training image-size distribution is a provenance defect, and the intended
WEIGHTS are not evidence of what a pool actually contains -- only the counts are.
"""
function realized_imsize_counts(imsize)
    counts = Dict{Tuple{Int,Int},Int}()
    for z in imsize
        counts[z] = get(counts, z, 0) + 1
    end
    return counts
end
