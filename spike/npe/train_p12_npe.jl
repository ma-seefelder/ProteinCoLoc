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

# spike/npe/train_p12_npe.jl --- THE single Phase-12 training surface.
#
# SCOPE. Every Phase-12 net is trained through this one file, and the PRIOR ARM IS A KEYWORD. CAR,
# GP and the D-10 `arm = :none` ablation are therefore the same code path on different pools, which
# is what makes SC3's matched-ablation claim ATTRIBUTABLE: a coverage difference between two arms
# cannot come from a different split, a refitted transform or a different masking rate, because
# there is only one of each. If 12-15, 12-16, 12-17, 12-18, 12-19 and 12-20 each wired their own
# loop, the arms would stop being comparable and the claim would collapse quietly.
#
# THE MASK ORDERING CORRECTION -- A DELIBERATE DEPARTURE FROM A SIBLING PLAN'S PROSE.
# `12-09-PLAN.md:69` says the MCAR masking is applied "to the standardized input". THAT IS WRONG
# AND THIS FILE DOES NOT FOLLOW IT. Masking happens on the RAW 128 rows, BEFORE `zt`, because the
# entire point of Pitfall 4 is that a held-out region must look EXACTLY like a region
# `correlation()` could not score. A naturally-unusable patch reaches the net as raw value `0.0`
# with mask `0`, which after per-row z-scoring becomes `(0 - mu_r)/sd_r` -- NOT zero. Masking after
# standardization would instead place the ROW MEAN at the held-out region, so the net would read
# "this region is perfectly average" exactly where read time says "this region is absent" -- an
# out-of-distribution encoding, and a silent one. D-09's own construction in 12-RESEARCH.md
# ("The leave-region-out predictive construction") masks `Z_full = encode_d01(patch_summary(mci))`,
# i.e. RAW, and confirms the ordering. This is a deviation from a sibling plan's PROSE, not from
# any decision. The order is `augment -> standardize -> reshape`, and step 6 computes it through a
# named local `Zraw_masked` so the ordering is visible in one expression rather than implied.
#
# `augment_mask!` IS **NOT** DEFINED HERE, AND THAT IS A CORRECTION TO THIS PLAN'S OWN TEXT.
# `12-14-PLAN.md:154` specifies `augment_mask!(Zraw, rng; G, kset)` as one of this file's functions.
# THAT FUNCTION ALREADY EXISTS -- `spike/npe/p12_architecture.jl:190`, written by 12-04, with the
# keyword spelled `k_set` -- and this file `include`s that one into the same scope. Defining it
# again here would SILENTLY OVERWRITE the pre-registered 12-04 augmentation with an identically
# named method, and whichever definition won would depend on include order. That is exactly the
# defect class 12-12 already found once in this phase (a bare helper whose zero-positional
# signature overwrote its Phase-13 namesake) and which `test_p12_architecture.jl` would then be
# testing the WRONG function for. So: 12-04's `augment_mask!` is CALLED, never redefined.
#
# The plan wanted the realized per-column `k` back from it (F6: "the intended distribution is not
# evidence of what a run contained"). 12-04's returns `Zraw`, so this file instead MEASURES what the
# net saw, from the DATA, via `_p12train_masked_counts` before and after the call. That is strictly
# better evidence than the drawn `k` would be -- it counts what the net actually saw -- and it keeps
# the augmentation single-sourced.
#
# WHAT THE MEASURED NUMBERS ARE, NAMED FOR WHAT THEY COUNT AND NOT FOR WHAT THEY ARE NEAR.
# NONE OF THE THREE IS THE DRAWN `k`, and calling any of them "realized k" would repeat, in
# miniature, the ridge-versus-posterior `shrinkage` collision this phase already named apart: two
# quantities that coincide in ONE regime, labelled as if they coincide always.
#   `natural_mask_rate`       -- zeroed mask rows per region BEFORE augmentation: the rate the
#                                forward model produced on its own (a patch below the >=15-surviving-
#                                pixel floor).
#   `realized_mask_rate`      -- zeroed mask rows per region AFTER augmentation. THE TOTAL the net
#                                saw, not a delta. This is the key name 12-14 mandates.
#   `augmented_only_mask_rate` -- the difference, i.e. NEWLY zeroed rows. This UNDER-COUNTS the drawn
#                                `k` by any overlap: if `augment_mask!` selects a region the forward
#                                model had already masked, that draw zeroes nothing new.
# `realized_mask_rate * G^2` equals the mean drawn `k` ONLY when the pre-augmentation mask has no
# zeros at all. On the F5 mixture natural masking is rare, so the two are usually close -- which is
# exactly why the distinction has to be written down rather than left to be inferred from a number
# that "looks about right". `test_p12_train.jl` testset 4 therefore exercises the k-set property on
# an ALL-ONES-mask fixture, where no natural zeros exist and the recount IS the drawn k.
#
# LEAK-FREEDOM. Both transforms are fitted in `fit_p12_transforms` and NOWHERE else -- one
# `fit(ZScoreTransform` in the whole file, asserted by this plan's acceptance criteria -- on the
# TRAIN split only, and are FROZEN into the saved bundle. Nothing downstream ever refits.
#
# CPU-ONLY. `use_gpu = false` is passed EXPLICITLY on every NeuralEstimators call: v0.2.1 defaults
# it to TRUE, which the executor briefing records as having already bitten 12-04.
#
# GATE ORDER IS A HARD PRECONDITION. D-12 Stage 1 exists to make the descope call BEFORE the
# training spend, so `_assert_spend_allowed` refuses a real run unless the Stage-1 verdict is on
# disk and reads PROCEED. The 1-epoch smoke path is exempt -- it is a unit test, not a spend -- and
# says so in its own message. On DESCOPE this file still runs, but ONLY with `arm = :none`; that
# asymmetry is 12-11's routing and is enforced, not documented.
#
# `p12_stage1_verdict` IS NOT DEFINED HERE EITHER. It is a Tier-1 function (12-01) reached through
# the guarded `p12_consts.jl` include, for the same overwrite reason as `augment_mask!` above.
# `_assert_spend_allowed` IS this file's own, because the arm asymmetry is not Tier 1's to own.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. Reaches `src/` not at all. THE TRAINER READS
# POOLS AND NEVER MAKES THEM -- `generate_p12_pool` does not appear in this file, and the pool
# resolution asserts completeness and names the generation command rather than running it.
#
# This file is a LIBRARY: it defines functions and calls none.

using NeuralEstimators   # train, PosteriorEstimator
using Flux               # optimiser
using StatsBase          # ZScoreTransform, fit, transform
using JLD2               # atomic bundle persistence
using Dates              # UTC-labeled artifact timestamp
using Statistics         # mean, std
using Random123          # the reserved counter-based stream
using Random             # seed! -- the WEIGHT-INIT seeding fix (see _p12_train_seed below)

# --- ORDER MATTERS: the pre-registration FIRST (the `p11_generate.jl:53-56` reason -- consts bind
#     reserved seeds at UInt64 width and a later narrower rebinding is a hard `const` error), then
#     the model surface, then the pool reader. Guarded for idempotency. Only Phase-12-local
#     sentinels are guarded on, so the include-guard poisoning class cannot arise here.
isdefined(@__MODULE__, :P12_DEV_SEED)        || include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))
isdefined(@__MODULE__, :build_p12_estimator) || include(joinpath(@__DIR__, "p12_architecture.jl"))
isdefined(@__MODULE__, :load_p12_pool)       || include(joinpath(@__DIR__, "..", "data", "p12_generate.jl"))

"Persistence schema for a Phase-12 bundle (bump on any breaking layout change)."
const P12_NET_SCHEMA = 1

"UTC, explicitly Z-labeled (IN-07): an unlabeled local timestamp was an audit defect."
_p12_now() = string(Dates.now(Dates.UTC)) * "Z"

# =============================================================================================
# 1. The Stage-1 gate guard -- the arm asymmetry, which Tier 1 does not own
# =============================================================================================

"""
    _assert_spend_allowed(arm, epochs, n; verdict_dir = nothing)

Refuse a real training spend unless the D-12 Stage-1 verdict is on disk and permits this arm.

THE SMOKE PATH IS EXEMPT. `epochs <= 1 && n <= 512` is a unit test, not a spend, and the skip
message says so -- otherwise every test run would depend on the repository's live verdict state.

`verdict_dir` THREADS STRAIGHT THROUGH to 12-01's Tier-1 `p12_stage1_verdict(; dir = ...)`, and
`nothing` means "take that function's own default", so no production call site changes. It exists
because 12-14's testset 7 must exercise the absent / DESCOPE / PROCEED branches WITHOUT writing
into the real `.planning/` tree -- and because 12-11 (wave 6) always writes a verdict before this
plan (wave 7) runs, so a bare call can never reach the `:absent` branch and that branch would
otherwise be untestable. Tier 1 is append-only, so this keyword and its Tier-1 counterpart had to
be declared together, before 12-01 executed.

The three branches:
  `:absent`   -> throw, naming BOTH missing paths and the command that produces them.
  `:descope`  -> throw UNLESS `arm === :none`, quoting 12-11's routing sentence verbatim.
  `:proceed`  -> return.
"""
function _assert_spend_allowed(arm::Symbol, epochs::Integer, n::Integer; verdict_dir = nothing)
    if epochs <= 1 && n <= 512
        @info "train_p12_npe: SMOKE PATH (epochs = $epochs <= 1, n = $n <= 512) — the Stage-1 " *
              "gate guard is SKIPPED. This is a unit test, not a compute spend; a real run at " *
              "any larger size or epoch count is gated."
        return :smoke
    end

    v = verdict_dir === nothing ? p12_stage1_verdict() : p12_stage1_verdict(; dir = verdict_dir)

    if v === :absent
        report = joinpath(@__DIR__, "..", "validation", "p12_stage1_report.jld2")
        error("""
            train_p12_npe: REFUSING A REAL TRAINING SPEND — the D-12 Stage-1 gate has not been
            adjudicated. Both of these must exist first:
                $(normpath(report))
                $(P12_STAGE1_VERDICT_FILE)  (under the phase's planning directory)
            Produce the first with:
                julia --project=spike -t auto spike/validation/run_p12_stage1_ridge.jl
            The second is the ADJUDICATION of that artifact and is written by a human ruling, not
            by this or any other script. D-12 Stage 1 exists precisely to make the descope call
            BEFORE the training spend; running the spend first would make the gate decorative.""")
    elseif v === :descope
        arm === :none || error("""
            train_p12_npe: the Stage-1 verdict is DESCOPE, on which `arm = :$(arm)` may not be
            trained. 12-11's routing sentence, verbatim:
                "If the D-12 Stage-1 gate does not clear its pre-registered bars the phase
                 DESCOPES: the mini-spike is not trained, and 12-17, 12-18, 12-19 and 12-20 do
                 not run."
            `arm = :none` is the ONE legal arm on this route: 12-11 routes the descope through
            12-12, 12-14 and 12-16 to build D-13's ablation fallback, and this trainer is the only
            thing on that route that trains anything.""")
        @info "train_p12_npe: Stage-1 verdict is DESCOPE; arm = :none is permitted (the D-13 " *
              "ablation fallback)."
    end
    return v
end

# =============================================================================================
# 2. Masking bookkeeping, the transforms, and the theta-row truncation operator
# =============================================================================================

"""
    _p12train_masked_counts(Zraw; G = P12_G) -> Vector{Int}

Per-column count of regions whose MASK row reads 0, i.e. regions the net will see as ABSENT.

Measured from the DATA rather than from the drawn `k`, which is what lets `realized_mask_rate` be
evidence of what a run CONTAINED rather than of what its sampler intended (F6). Called before and
after 12-04's `augment_mask!` so the augmentation's own contribution is separable from the rate
the forward model produced naturally.
"""
function _p12train_masked_counts(Zraw::AbstractMatrix; G::Integer = P12_G)
    n = G^2
    size(Zraw, 1) == 2n || throw(DimensionMismatch(
        "_p12train_masked_counts: expected $(2n) rows (2·G² at G = $G), got $(size(Zraw, 1))"))
    return [count(iszero, view(Zraw, (n + 1):(2n), c)) for c in axes(Zraw, 2)]
end

"""
    fit_p12_transforms(Zraw_train, theta_train) -> (zt, theta_zt)

THE ONLY PLACE EITHER TRANSFORM IS FITTED, on the TRAIN split only.

Per `P12_ZSCORE_ARM === :per_row` the summary transform is fitted over the 64 CONTINUOUS rows
only; the 64 mask rows are never fitted and never applied to, mirroring
`src/amortized/summary.jl:85-98`'s recorded reason -- z-scoring a 0/1 mask would re-couple folds
through the mask mean.

R-3's shared-scalar alternative is reachable by changing `P12_ZSCORE_ARM`; it is NOT taken by
default, and whichever arm ran is persisted into the bundle as `zscore_arm` so a later comparison
is readable rather than inferred.

The zero-variance guard is `run_p11_recovery.jl:75-84`'s: a constant row gets scale 1.0 rather
than dividing by zero. Still leak-free -- the fit saw only train columns.
"""
function fit_p12_transforms(Zraw_train::AbstractMatrix, theta_train::AbstractMatrix;
                            G::Integer = P12_G)
    n = G^2
    size(Zraw_train, 1) == 2n || throw(DimensionMismatch(
        "fit_p12_transforms: expected $(2n) raw summary rows, got $(size(Zraw_train, 1))"))
    P12_ZSCORE_ARM === :per_row || error(
        "fit_p12_transforms: only the pre-registered :per_row arm is implemented; " *
        "P12_ZSCORE_ARM is :$(P12_ZSCORE_ARM). R-3's shared-scalar arm is a declared " *
        "alternative and adopting it is a recorded deviation, not a silent branch.")

    zt = StatsBase.fit(StatsBase.ZScoreTransform, Zraw_train[1:n, :]; dims = 2)
    @inbounds for i in eachindex(zt.scale)
        (isfinite(zt.scale[i]) && zt.scale[i] > 0) || (zt.scale[i] = one(eltype(zt.scale)))
    end
    theta_zt = StatsBase.fit(StatsBase.ZScoreTransform, theta_train; dims = 2)
    @inbounds for i in eachindex(theta_zt.scale)
        (isfinite(theta_zt.scale[i]) && theta_zt.scale[i] > 0) ||
            (theta_zt.scale[i] = one(eltype(theta_zt.scale)))
    end
    return zt, theta_zt
end

"""
    standardize_p12(Zraw, zt; G = P12_G) -> Matrix{Float64}

Apply the FROZEN `zt` to the 64 continuous rows and pass the 64 mask rows through UNTOUCHED.
NEVER calls `fit` -- that lives in `fit_p12_transforms` and only there.
"""
function standardize_p12(Zraw::AbstractMatrix, zt; G::Integer = P12_G)
    n = G^2
    size(Zraw, 1) == 2n || throw(DimensionMismatch(
        "standardize_p12: expected $(2n) rows (2·G² at G = $G), got $(size(Zraw, 1))"))
    out = Matrix{Float64}(undef, size(Zraw))
    out[1:n, :]       = StatsBase.transform(zt, Float64.(Zraw[1:n, :]))
    out[(n + 1):2n, :] = Float64.(Zraw[(n + 1):2n, :])   # mask rows: never re-z-scored
    return out
end

"""
    p12_truncate_theta(theta, K) -> Matrix

THE theta-row truncation operator, and the ONE place it exists.

The pool's theta is always `P12_D_MINISPIKE` = 72 rows, laid out by `p12_theta_column` as
`[c0; deviations 1:K_dev; the 7 nuisances; r1]`. A production head at `K < P12_K_DEV` is that same
layout with the DEVIATION BLOCK SHORTENED, so the operator is a contiguous row selection returning
`1 + K + 7 + 1` rows. `K == P12_K_DEV` returns the input unchanged, so the default path is provably
a no-op.

WHY THIS EXISTS, recorded because its absence survived three plan-checker passes. 12-15
pre-registers `P12_D_PROD = 1 + P12_K_PROD + 7 + 1` and 12-17 trains "the production head", but no
plan said how 72 pool rows become `D_PROD` rows. Without it 12-17 is UNBUILDABLE: the trainer would
fit `D = size(theta, 1) = 72` and write a 72-entry `theta_rows`, and `load_p12_npe` -- which
REFUSES, never repairs, a bundle whose `theta_rows` disagrees with its `D` -- would reject the
production bundle for every downstream consumer. A gap that fires only at `K_PROD < 63` is exactly
the kind that reaches execution.
"""
function p12_truncate_theta(theta::AbstractMatrix, K::Integer)
    K <= P12_K_DEV || throw(ArgumentError(
        "p12_truncate_theta: K = $K exceeds P12_K_DEV = $(P12_K_DEV). The production head is a " *
        "TRUNCATION of the mini-spike layout, never wider (12-15's P12_D_PROD <= P12_D_MINISPIKE " *
        "in row terms)."))
    K >= 1 || throw(ArgumentError("p12_truncate_theta: K must be ≥ 1, got $K"))
    size(theta, 1) == P12_D_MINISPIKE || throw(DimensionMismatch(
        "p12_truncate_theta: expected $(P12_D_MINISPIKE) pool theta rows, got $(size(theta, 1))"))
    K == P12_K_DEV && return theta                       # provable no-op on the default path
    out = vcat(theta[1:1, :], theta[2:(1 + K), :], theta[(end - 7):end, :])
    @assert size(out, 1) == 1 + K + 7 + 1 "p12_truncate_theta: produced $(size(out, 1)) rows, " *
        "expected $(1 + K + 7 + 1)"
    rows = P12_THETA_ROWS(K)
    @assert length(rows) == size(out, 1)
    @assert first(rows) === :c0 && last(rows) === :r1 "p12_truncate_theta: the truncated layout " *
        "must still start at c0 and end at r1"
    return out
end

# =============================================================================================
# 3. The training surface
# =============================================================================================

"""
    train_p12_npe(; arm, n, K_dev, pool_dir, epochs, batchsize, val_fraction, out_path,
                  verbose, verdict_dir, rng) -> NamedTuple

Train ONE Phase-12 net and persist it atomically. The arm is a keyword; everything else is
identical across arms by construction.

`K_dev` IS THE ONLY KNOB THAT SETS THE HEAD WIDTH. It defaults to `P12_K_DEV`, so the mini-spike is
unchanged; 12-17 passes `K_dev = P12_K_PROD`. There is deliberately NO separate `D` keyword -- `D`
is DERIVED from `K_dev` by the layout, and two independent width arguments could disagree. The name
is `K_dev` rather than bare `K` because `K` already means the number of datasets in a batch in
`reshape_summary`'s `(G, G, 2, K)` shape, and because `K_dev` is also the bundle key this value is
written to, so one name serves the argument, the bundle and the docs.

`verdict_dir` is forwarded verbatim to `_assert_spend_allowed`; see there.
"""
function train_p12_npe(; arm::Symbol = :car,
                       n::Integer = P12_MINISPIKE_N,
                       K_dev::Integer = P12_K_DEV,
                       pool_dir = nothing,
                       pool_indices = nothing,
                       epochs::Integer = 18,
                       batchsize::Integer = 64,
                       val_fraction::Real = 0.2,
                       out_path = nothing,
                       verbose::Bool = true,
                       verdict_dir = nothing,
                       rng = p12_rng(P12_MINISPIKE_COUNTER))
    t0 = time()

    # 1. THE GATE, FIRST. Nothing below runs before the spend is permitted.
    verdict = _assert_spend_allowed(arm, epochs, n; verdict_dir = verdict_dir)

    # 2. Resolve the pool. THE TRAINER READS POOLS; IT DOES NOT MAKE THEM.
    dir = pool_dir === nothing ? p12_pool_dir(n; arm = arm) : pool_dir
    isdir(dir) && p12_pool_complete(dir, n) || error("""
        train_p12_npe: no COMPLETE pool for arm = :$(arm), n = $n at
            $dir
        This trainer READS pools and refuses to generate one -- generating from inside a training
        call is how an arm silently acquires a pool the comparison assumed it shared. Produce it
        through the ONE datagen path, `spike/data/p12_generate.jl`, whose pool-building entry point
        takes `n` and the `arm` keyword; the plans that own that spend are 12-15 (per-arm
        mini-spike pools) and 12-17 (the full production pool).

        THE GENERATING SYMBOL IS NAMED BY FILE RATHER THAN SPELLED HERE, DELIBERATELY. 12-14's own
        acceptance criteria require that this file contain ZERO occurrences of that identifier once
        comments are stripped, while also asking this message to name the generation command -- two
        requirements that cannot both be met literally. Naming the file and the owning plans
        satisfies the reader without putting the symbol in the trainer's source, which is what
        makes "the trainer never generates" checkable by a text scan rather than by review.

        NOTE `isdir` ALONE IS NOT A PRESENCE TEST: `p12_pool_dir` resolves through
        `open_or_invalidate`, which CREATES the directory it resolves (12-09 section 5), so an
        absent pool appears as a freshly-created EMPTY hash directory. `p12_pool_complete` is the
        real check and is the one that fired here.""")

    if verbose
        println("="^78)
        println("Phase-12 NPE training — arm = :$arm   n = $n   K_dev = $K_dev")
        println("  pool     : $dir")
        println("  verdict  : $verdict")
        println("  threads  : $(Threads.nthreads())   CPU-ONLY (use_gpu = false explicitly)")
        println("="^78)
    end

    pool  = load_p12_pool(dir)
    Zpool = Float64.(pool.summary_min)          # 128 x N, RAW (never z-scored at rest)
    N     = size(Zpool, 2)
    @assert size(Zpool, 1) == 2 * P12_G^2 "train_p12_npe: expected $(2 * P12_G^2) summary rows, " *
        "got $(size(Zpool, 1))"

    # 3. Deterministic TAIL-BLOCK validation split, no shuffle, so two arms at the same `n` get
    #    the SAME index partition -- which is what makes the D-10 matched ablation matched.
    # THE AVAILABLE COLUMNS. `pool_indices === nothing` means ALL of them, which is the behaviour
    # every existing call site had before this keyword existed -- so the default path is provably
    # unchanged, and `avail == 1:N` is asserted below rather than assumed.
    #
    # WHY THIS KEYWORD EXISTS, AND WHY THE OBVIOUS ALTERNATIVE IS A LEAK. 12-15 needs a scoring
    # block disjoint from everything training saw. It cannot get one from a second pool:
    # `generate_p12_sample(idx)` keys its RNG on `p12_datagen_rng(idx)`, the GLOBAL INDEX ALONE,
    # and `generate_p12_pool(n)` always generates indices 1:n -- so two pools at the same
    # arm/imsize share samples 1:min(n,m) BYTE-IDENTICALLY. A "separate scoring pool" is a
    # SUPERSET of the training pool, never a held-out set, and there is no index-offset keyword to
    # escape that. Verified executably before this keyword was added. Carving the block out of the
    # ONE pool and handing training the remainder is therefore the only leak-free construction,
    # and it is also 12-15's own primary instruction.
    avail = pool_indices === nothing ? collect(1:N) : collect(pool_indices)
    if pool_indices === nothing
        @assert avail == collect(1:N)
    else
        allunique(avail) || error("train_p12_npe: pool_indices contains duplicates; a repeated " *
            "column would be trained on twice and would corrupt the val split's disjointness")
        (minimum(avail) >= 1 && maximum(avail) <= N) || error(
            "train_p12_npe: pool_indices range $(extrema(avail)) is outside the pool's 1:$N")
        length(avail) >= 2 || error("train_p12_npe: pool_indices must select at least 2 columns")
        sort!(avail)
    end
    Navail = length(avail)

    nval      = round(Int, val_fraction * Navail)
    val_idx   = avail[(Navail - nval + 1):Navail]      # deterministic TAIL block of the available
    train_idx = avail[1:(Navail - nval)]               # columns, no shuffle
    @assert isempty(intersect(train_idx, val_idx)) "train_p12_npe: train and val index sets overlap"
    verbose && println("[1/6] split: $(length(train_idx)) train / $(length(val_idx)) val " *
                       "of $Navail available column(s)" *
                       (pool_indices === nothing ? " (whole pool)" :
                        " (SUBSET of a $N-column pool; the complement is the caller's held-out set)") *
                       " — deterministic tail block, no shuffle")

    # 4. TRUNCATE theta FIRST, so `theta_zt` is fitted over exactly the rows the head predicts and
    #    never over rows that are then discarded.
    theta = p12_truncate_theta(Float64.(pool.theta), K_dev)
    D     = size(theta, 1)
    @assert D == 1 + K_dev + 7 + 1
    verbose && println("[2/6] theta truncated to K_dev = $K_dev -> D = $D " *
                       "(P12_K_DEV = $(P12_K_DEV), P12_D_MINISPIKE = $(P12_D_MINISPIKE))")

    # 5. Fit BOTH transforms on the TRAIN split ONLY. The one and only fit site.
    zt, theta_zt = fit_p12_transforms(view(Zpool, :, train_idx), view(theta, :, train_idx))
    verbose && println("[3/6] transforms fitted TRAIN-ONLY " *
                       "($(length(zt.mean)) continuous rows z-scored, $(P12_G^2) mask rows " *
                       "bypassed; theta over $D rows)")

    # 6. AUGMENT ON RAW -> STANDARDIZE -> RESHAPE. The ordering correction, computed through the
    #    named local `Zraw_masked` so the order is visible in ONE expression rather than implied.
    #    `augment_mask!` is 12-04's (p12_architecture.jl:190) and is CALLED, never redefined.
    rng_train = p12_rng(P12_MINISPIKE_COUNTER)
    rng_val   = p12_rng(P12_MINISPIKE_COUNTER + 1)   # a DISTINCT stream: train and val masks must
                                                      # be independent, or the validation loss is
                                                      # measured on the training masks.
    Zraw_train_masked = copy(Zpool[:, train_idx])
    Zraw_val_masked   = copy(Zpool[:, val_idx])
    natural_train = _p12train_masked_counts(Zraw_train_masked)
    natural_val   = _p12train_masked_counts(Zraw_val_masked)
    augment_mask!(Zraw_train_masked, rng_train)
    augment_mask!(Zraw_val_masked,   rng_val)
    realized_train = _p12train_masked_counts(Zraw_train_masked)
    realized_val   = _p12train_masked_counts(Zraw_val_masked)

    Ztrain = reshape_summary(standardize_p12(Zraw_train_masked, zt), P12_G)
    Zval   = reshape_summary(standardize_p12(Zraw_val_masked,   zt), P12_G)

    # MEASURED, not intended (F6). See the header for what each of these three IS: none is the
    # drawn `k`, and `realized_mask_rate` is the POST-augmentation TOTAL rather than a delta.
    nreg = P12_G^2
    realized_mask_rate       = mean(realized_train) / nreg
    natural_mask_rate        = mean(natural_train)  / nreg
    augmented_only_mask_rate = mean(realized_train .- natural_train) / nreg
    verbose && println("[4/6] MCAR augmentation on RAW rows (Pitfall 4): zeroed-row rate AFTER = " *
                       "$(round(realized_mask_rate; digits = 5)) of $nreg regions " *
                       "(natural rate BEFORE = $(round(natural_mask_rate; digits = 5)); " *
                       "newly zeroed = $(round(augmented_only_mask_rate; digits = 5)) — this " *
                       "under-counts the drawn k by any overlap with naturally masked regions)")

    theta_train_std = Float32.(StatsBase.transform(theta_zt, theta[:, train_idx]))
    theta_val_std   = Float32.(StatsBase.transform(theta_zt, theta[:, val_idx]))

    # --- AUDITABILITY FIX 1 of 2: SEED THE WEIGHT INITIALISATION -------------------------------
    # RECORDED DEFECT, found while harvesting 12-15: `build_p12_estimator` takes no RNG and
    # `NeuralEstimators.train` was called without one, so Flux drew its weight init from the
    # UNSEEDED global RNG. Only the mask streams and the deterministic tail-block split were
    # reproducible, which made two runs at identical settings non-comparable -- and left 12-15's
    # decisive 0.33 % arm margin unattributable.
    #
    # THIS IS AN AUDITABILITY FIX, NOT TUNING. It does not change WHAT the model can learn; it makes
    # what it did learn reproducible. Recorded as the fix to a defect, never as an adjustment made
    # in pursuit of a better number.
    #
    # COUNTER +2 OFF THE EXISTING RESERVED STREAM -- so NO new Tier-1 constant and NO
    # pre-registration edit. `P12_MINISPIKE_COUNTER + 0` and `+ 1` are the train/val mask streams
    # above; `+ 2` was free.
    train_seed = rand(p12_rng(P12_MINISPIKE_COUNTER + 2), UInt64)
    Random.seed!(train_seed)

    # 7. D comes from the TRUNCATED theta, never a literal and never from the raw pool.
    est = build_p12_estimator(; G = P12_G, D = D)
    verbose && println("[5/6] estimator built (G = $(P12_G), D = $D); training $epochs epoch(s) " *
                       "at train_seed = $train_seed ...")

    # --- AUDITABILITY FIX 2 of 2: PERSIST THE PER-EPOCH RISK TRACE -----------------------------
    # RECORDED DEFECT: the bundle stored `epochs` and `batchsize` but NO per-epoch risk, and
    # `train` returns only the estimator. 12-15's convergence evidence -- the train/val ratios that
    # showed the run ended on its EPOCH BUDGET rather than on convergence -- survived only because
    # that run's stdout happened to be redirected to a scratch log. A training run whose
    # convergence cannot be audited afterwards cannot support a claim about convergence.
    #
    # `train` ALREADY computes this and writes `loss_per_epoch.csv` under `savepath` (headerless,
    # (epochs+1) x 2: training risk, validation risk; row 1 is the initial validation risk twice).
    # Nothing ever read it back. Parsed by hand rather than through CSV.jl so this adds no dependency.
    lossdir = mktempdir()
    est = NeuralEstimators.train(est, theta_train_std, theta_val_std, Ztrain, Zval;
                                 epochs = epochs,
                                 batchsize = batchsize,
                                 use_gpu = false,          # EXPLICIT: v0.2.1 defaults this to TRUE
                                 savepath = lossdir,
                                 verbose = verbose)

    risk_trace = let f = joinpath(lossdir, "loss_per_epoch.csv")
        isfile(f) ? reduce(vcat, [permutedims(parse.(Float64, split(strip(l), ',')))
                                  for l in readlines(f) if !isempty(strip(l))]) :
                    Matrix{Float64}(undef, 0, 2)
    end
    @assert size(risk_trace, 2) == 2 || isempty(risk_trace) "train_p12_npe: the risk trace must " *
        "be (epochs+1) x 2 -- training risk, validation risk -- got $(size(risk_trace))"

    elapsed_min = (time() - t0) / 60
    bundle = (schema_version = P12_NET_SCHEMA,
              estimator = est,
              zt = zt,                      # these two field names deliberately match the shipped
              theta_zt = theta_zt,          # bundle so harness.jl consumers need no special case
              arm = arm,
              r1_prior = arm === :none ? P12_ABLATION_R1 : (P12_R1_MIN, P12_R1_MAX),
              zscore_arm = P12_ZSCORE_ARM,
              theta_rows = P12_THETA_ROWS(K_dev),   # the ARGUMENT IS MANDATORY: a bare
                                                    # P12_THETA_ROWS() writes 72 names into a
                                                    # bundle whose D may be P12_D_PROD, and
                                                    # load_p12_npe refuses exactly that bundle
              D = D,
              G = P12_G,
              K_dev = K_dev,
              n_pool = N,
              n_train = length(train_idx),
              n_val = length(val_idx),
              # THE REALIZED INDEX SETS, so a scoring runner can ASSERT its held-out block is
              # disjoint from everything training saw rather than re-derive the arithmetic and
              # trust it. Both are stored because the val block drives early stopping inside
              # `NeuralEstimators.train` and is therefore just as contaminated as the train block
              # for scoring purposes -- a scorer must exclude BOTH.
              train_indices = collect(train_idx),
              val_indices = collect(val_idx),
              pool_indices_given = pool_indices === nothing ? nothing : collect(avail),
              epochs = epochs,
              batchsize = batchsize,
              # THE TWO AUDITABILITY FIELDS. `train_seed` makes the run reproducible; `risk_trace`
              # makes its convergence checkable AFTER the fact instead of only in a live stdout.
              # Both are ADDED fields -- `schema_version` is deliberately NOT bumped, because
              # `load_p12_npe` asserts `schema_version == 1` and checks for the PRESENCE of its
              # required keys, so an added key is backward-compatible while a bumped version would
              # make every existing bundle unloadable.
              train_seed = train_seed,
              risk_trace = risk_trace,          # (epochs+1) x 2: [training risk  validation risk]
              mask_k_set = collect(P12_MASK_K_SET),
              realized_mask_rate = realized_mask_rate,
              natural_mask_rate = natural_mask_rate,
              augmented_only_mask_rate = augmented_only_mask_rate,
              mask_rate_naming_note =
                  "NONE of these three is the drawn k. `realized_mask_rate` is the fraction of " *
                  "regions whose MASK ROW READS ZERO AFTER augmentation (a TOTAL, not a delta); " *
                  "`natural_mask_rate` is the same count BEFORE augmentation, produced by the " *
                  "forward model's >=15-surviving-pixel floor; `augmented_only_mask_rate` is the " *
                  "difference, i.e. NEWLY zeroed rows, which UNDER-COUNTS the drawn k by any " *
                  "overlap with an already-masked region. realized_mask_rate * G^2 equals the mean " *
                  "drawn k ONLY when the pre-augmentation mask has no zeros. Do not compare any of " *
                  "them against P12_MASK_K_SET without that condition holding — that would be the " *
                  "same unit error as conflating a ridge residual ratio with a posterior width.",
              realized_mask_counts_train = realized_train,
              realized_mask_counts_val = realized_val,
              natural_mask_counts_train = natural_train,
              natural_mask_counts_val = natural_val,
              pool_dir = string(dir),
              stage1_verdict = verdict,
              elapsed_min = elapsed_min,
              generated = _p12_now(),
              julia_version = string(VERSION),
              caption = "Phase-12 NPE, arm = :$arm, head width K_dev = $K_dev (D = $D). MCAR " *
                        "region masking applied to the RAW rows BEFORE standardization so a " *
                        "held-out region is encoded byte-identically to a naturally unusable " *
                        "patch (Pitfall 4). Both transforms fitted TRAIN-ONLY and frozen here. " *
                        "CPU-only.")

    if out_path !== nothing
        save_p12_npe(out_path, bundle)
        verbose && println("[6/6] bundle -> $out_path  ($(round(elapsed_min; digits = 2)) min)")
    else
        verbose && println("[6/6] not persisted (out_path = nothing)  " *
                           "($(round(elapsed_min; digits = 2)) min)")
    end
    return bundle
end

"""
    save_p12_npe(path, bundle) -> String

Atomic persistence: `jldsave` to a `.tmp`, reopen and integrity-check, then `mv(...; force = true)`
(a filesystem-atomic rename). A crash before the `mv` leaves only a discardable `.tmp`, never a
half-written model.
"""
function save_p12_npe(path, bundle::NamedTuple)
    mkpath(dirname(abspath(path)))
    tmp = string(path) * ".tmp"
    jldsave(tmp; (k => getfield(bundle, k) for k in propertynames(bundle))...)
    JLD2.jldopen(tmp, "r") do f
        for k in ("estimator", "zt", "theta_zt", "theta_rows", "D", "arm", "realized_mask_rate")
            haskey(f, k) || error("save_p12_npe: integrity check failed, $tmp missing $k")
        end
    end
    mv(tmp, path; force = true)
    return string(path)
end

"""
    load_p12_npe(path) -> NamedTuple

The read side. Asserts `schema_version == 1` and that `theta_rows` has length `D`.

A BUNDLE WHOSE `theta_rows` DISAGREES WITH ITS `D` IS REFUSED, NOT REPAIRED. That refusal is what
makes `p12_truncate_theta` load-bearing: without the truncation operator a production head at
`K_PROD < 63` would be written with 72 row names against a smaller `D`, and every downstream
consumer would be handed a bundle whose row labels do not describe its own output.
"""
function load_p12_npe(path)
    isfile(path) || error("load_p12_npe: no bundle at $path")
    d = JLD2.load(path)
    haskey(d, "estimator") || error("load_p12_npe: $path missing estimator key")
    d["schema_version"] == P12_NET_SCHEMA || error(
        "load_p12_npe: $path carries schema_version $(d["schema_version"]), expected $(P12_NET_SCHEMA)")
    length(d["theta_rows"]) == d["D"] || error(
        "load_p12_npe: $path carries $(length(d["theta_rows"])) theta_rows against D = $(d["D"]). " *
        "A bundle whose row labels do not describe its own output is REFUSED, not repaired — it " *
        "is the signature of a production head written without p12_truncate_theta.")
    return NamedTuple(Symbol(k) => v for (k, v) in d)
end
