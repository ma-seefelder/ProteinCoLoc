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

# spike/p13/datagen.jl --- Phase-13 labelled pool on the Phase-11 basis (D-02, D-03, D-07, D-11).
#
# =============================================================================================
# (a) THE FROZEN-`zt` DISCIPLINE. BOTH summaries in EVERY pair are standardized with the FROZEN
#     PHASE-11 `zt`, INHERITED through `load_p13_basis()` (`preconditions.jl`). NO
#     Z-SCORE TRANSFORM IS FITTED ANYWHERE IN THIS FILE OR ANYWHERE IN THE PHASE-13 DATA PATH.
#     Re-fitting is catastrophic AND INVISIBLE: the net still receives numbers of the right shape
#     and roughly the right scale, it trains, it converges, and every claim about it silently
#     refers to a basis the net was never trained on.
#     THE WARNING SIGN (Pitfall 5) is standardized PHASE-13 summaries whose per-row moments come
#     out at mean ~= 0 AND sd ~= 1 ON THE PHASE-13 POOL -- a transform frozen on the PHASE-11
#     pool should not reproduce exactly that on a different pool. `check_frozen_zt` below MEASURES
#     that warning sign and RETURNS it as a number, so it is reported rather than hoped for.
#     READ THE CAVEAT IN `check_frozen_zt`'s DOCSTRING BEFORE INTERPRETING THE VERDICT: D-02
#     deliberately draws the Phase-13 pool from the SAME joint the Phase-11 net trained on
#     (train-joint == eval-joint, the F5 invariant), so the two pools are draws from the same
#     marginal and the SOFT moment heuristic loses much of its power here. The HARD half of the
#     check -- `zt` OBJECT IDENTITY against the memoized single load -- is the unambiguous one.
#
# (b) STRATIFICATION IS BY WHOLE CLASS, NEVER WITHIN CLASS (D-07-i,
#     `P13_STRATIFICATION = :class_frequency`). Whole labelled items are ACCEPTED or REJECTED to
#     reach the pre-registered `P13_TARGET_CLASS_FREQ`; the within-class theta distribution is
#     left EXACTLY as pi gives it. THE ONE-LINE REASON: the scalar per-head correction
#     `log(q_pos/q_neg)` is exact IFF stratification changed only class FREQUENCIES, because a
#     WITHIN-class distributional change enters INSIDE the evidence integral and depends on Z --
#     so no scalar can repair it (13-RESEARCH F2, demonstrated empirically by plan 13-07's
#     NEGATIVE CONTROL 2). Reshaping rho inside the exclusion range is therefore NOT a stronger
#     version of this design; it is a DIFFERENT design (`P13_STRATIFICATION_FALLBACK =
#     :importance_weighted`) that is spendable only on the pre-declared D-04 trigger.
#
# (c) LAMBDA IS APPENDED ONCE, AFTER `pair_encode` (D-03,
#     `P13_LAMBDA_PLACEMENT = :append_after_pair_encode`). Every conditioned column here is built
#     by exactly one call to `p13_encode_pair`, which is the single place the append happens.
#     Appending BEFORE the pair encoding pushes the conditioning value through the contrast block
#     (and trips `p13_pair_encode`'s `iseven` guard, which is the GOOD failure); appending TWICE
#     is SILENT, which is why the realized width is re-asserted against
#     `three_way_input_dim(G, p13_conditioning_length())` on the way out.
#     ONE lambda per pair: the sample stack and the control stack are two channels of the SAME
#     acquisition under the SAME stated registration uncertainty.
#
# (d) TWO-LAYER RNG DISCIPLINE (`P13_GLOBAL_RNG_DISCIPLINE`). LAYER 1 lives here: every pool
#     column is a pure function of its GLOBAL INDEX through a counter-based Philox stream
#     (`p13_datagen_rng(idx)`), so the pool is BYTE-IDENTICAL for any thread count, any execution
#     order and any number of resumed runs. `Random.seed!` IS NEVER CALLED IN THIS FILE -- not at
#     the top, and above all not inside a generation loop, which would correlate the draws and
#     destroy exactly the reproducibility the counter-based stream buys. LAYER 2 (a
#     `Random.seed!` immediately before `train`, because Flux's shuffling `DataLoader` draws from
#     the GLOBAL RNG) belongs to the TRAINER, `spike/p13/train_three_way.jl`, and only there.
#
# (e) DECOUPLING (hard constraint, CLAUDE.md): spike-local. Reaches `src/` only transitively and
#     READ-ONLY through the simulator/summary contract chain and through the content-hash guard's
#     `read` of frozen source bytes. It writes only into the gitignored bulk cache
#     `spike/data/cache/p13/`, imports no `ProteinCoLoc` symbol, adds no package, and never
#     imports CUDA. `src/`, `spike/Project.toml` and `spike/Manifest.toml` stay byte-unchanged.
#
# THE COST SHAPE, STATED SO A SLOW RUN IS NOT MISTAKEN FOR A HANG. The class label is a function
# of THETA ALONE (`three_way_label(rho_s, rho_c)`), so the whole-class accept/reject decision is
# made BEFORE any image is simulated: a rejected draw costs two `sample_prior` calls and nothing
# else. Only accepted items are simulated. That is not an optimisation of the design, it is the
# design read literally -- the item is rejected, and a rejected item has no images.
#
# Flat top-level functions (sibling style of `net.jl` / `preconditions.jl` -- no module wrapper).
# Guarded, order-dependent includes keep the file loadable standalone AND idempotent under
# `runtests.jl`.

using Statistics         # mean, std (the realized-frequency and moment reporting)
using Dates              # UTC-labelled artifact timestamp
using Random123          # Philox4x -- the counter-based per-index stream
using Random             # AbstractRNG, Xoshiro (the derived theta arity)
using Distributions      # Uniform (the lambda level and the conditional shift)
using JLD2               # jldsave / jldopen / load -- pure-Julia, Windows-clean persistence
using SHA                # stdlib: provenance sha of the frozen basis and the frozen consts

# ORDER MATTERS. `preconditions.jl` FIRST: it carries the WHOLE preamble this file needs -- the
# isolated `_P11C` read of Phase 11's Tier-1 pre-registration (LAMBDA_MIN/LAMBDA_MAX), the
# Phase-11 architecture (so JLD2 cannot hand back a reconstructed placeholder for the frozen
# theta transform), the Phase-13 pre-registration, the NPE loader, the standardization surface,
# the net/encoding surface, the simulator halves, the src/ summary contract and the encoder.
# Then the label surface, then the net (already pulled in by preconditions, guarded anyway so the
# dependency is stated), then the sharded content-hash cache layer.
isdefined(@__MODULE__, :p13_require_phase11)  || include(joinpath(@__DIR__, "preconditions.jl"))
isdefined(@__MODULE__, :three_way_label)      || include(joinpath(@__DIR__, "labels.jl"))
isdefined(@__MODULE__, :ThreeWayEvidenceNet)  || include(joinpath(@__DIR__, "net.jl"))
isdefined(@__MODULE__, :shard_path)           || include(joinpath(@__DIR__, "..", "data", "cache.jl"))

if !isdefined(@__MODULE__, :P13_POOL_SCHEMA)

    # =========================================================================================
    # LOCAL KNOBS. NOT PRE-REGISTERED BARS -- no reported number depends on any of them.
    # =========================================================================================
    # They live here rather than in `consts.jl` because `consts.jl` IS the pre-registration and
    # must not grow a knob no claim depends on. Everything a claim DOES depend on -- P13_TRAIN_N,
    # the recipe, tau, the target class frequencies, the imsize mixture -- is READ from there.

    "On-disk schema tag for a Phase-13 labelled pool shard. Bump iff the shard key set changes."
    const P13_POOL_SCHEMA = 1

    """
        P13_POOL_CACHE_ROOT

    Root of the Phase-13 labelled-pool cache. Sits under the gitignored bulk-cache tree
    (`.gitignore`: `spike/data/cache/*`), in its OWN `p13/` subtree.

    DELIBERATELY A DIFFERENT DIRECTORY FROM `spike/data/cache/p11/`. That one holds the 50k
    Phase-11 pool a Phase-12 plan depends on; it is READ-ONLY to this phase and is never pruned,
    regenerated or invalidated from here. The `p13/` prefix makes a content-hash collision with
    it structurally impossible rather than merely unlikely.
    """
    const P13_POOL_CACHE_ROOT = joinpath(@__DIR__, "..", "data", "cache", "p13")

    """
        P13_POOL_SHARD_SIZE

    Items per shard. Bounds the RESUME granularity: a killed run re-does at most one shard.
    Sized so one shard is a few minutes of wall clock at the F5 image-size mixture, which makes
    progress legible and a kill cheap. Smaller than `cache.jl`'s `SHARD_SIZE = 10_000` because a
    Phase-13 item is TWO simulations, not one.
    """
    const P13_POOL_SHARD_SIZE = 4_000

    """
        P13_STRAT_MAX_DRAW_FACTOR

    Safety bound on the whole-class accept/reject search, as a multiple of the requested pool
    size. The EXPECTED factor is ~1.54 at tau = 0.15 (see `stratify_by_class`); this bound exists
    only so a mis-specified target or a broken label rule fails with a diagnosis instead of
    spinning forever. It is NOT a tuning knob and no reported number depends on it.
    """
    const P13_STRAT_MAX_DRAW_FACTOR = 25

    # --- Default filenames inside a resolved pool directory ---------------------------------
    const P13_POOL_META_FILE = "meta.jld2"

end

# =============================================================================================
# 1. The reserved Phase-13 datagen stream (D-01, layer 1)
# =============================================================================================

"""
    p13_datagen_rng(idx::Integer) -> Philox4x

The per-item keyed RNG of the Phase-13 labelled pool:
`Philox4x((P13_DEV_SEED xor P13_SALT xor P13_DATAGEN_COUNTER, idx))`.

WHY THE COUNTER GOES IN THE **KEY** WORD AND THE INDEX IN THE COUNTER WORD. `p13_rng(counter)`
puts the activity counter in the COUNTER word, which is right for an activity that consumes ONE
sequential stream (the tau probe, the gate). A POOL cannot do that: the counter word must carry
the per-item GLOBAL INDEX, because that is exactly what makes a column a pure function of its
index and therefore the pool byte-identical for any thread count. So the reserved
`P13_DATAGEN_COUNTER` is folded into the KEY word instead, which is the same trick Phase 11 used
with a dedicated salt (`spike/data/p11_generate.jl:90-114`) -- here expressed WITHOUT adding a
constant to the frozen pre-registration, because every ingredient is already in it.

DISJOINTNESS IS ASSERTED, NOT ASSUMED (below): the resulting key word differs from every
`p13_rng(c)` key word (those all share `P13_DEV_SEED xor P13_SALT`) and from the whole fixture
family (which rides `P13_FIX_SEED`). A repeated key is a repeated stream, and a reported pool
sharing a stream with a reported evaluation is the D-01 failure this arithmetic prevents.
"""
p13_datagen_rng(idx::Integer) =
    Philox4x(UInt64, (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_DATAGEN_COUNTER),
                      UInt64(idx)))

"""
    p13_fix_datagen_rng(idx::Integer) -> Philox4x

The FIXTURE twin of `p13_datagen_rng`, riding `P13_FIX_SEED` at `P13_FIXTURE_COUNTER`. Fixtures
must never consume -- and so never pre-observe -- a reported stream, which is why the test file
draws from here and the reported pool never does.
"""
p13_fix_datagen_rng(idx::Integer) =
    Philox4x(UInt64, (UInt64(P13_FIX_SEED) ⊻ P13_SALT ⊻ UInt64(P13_FIXTURE_COUNTER),
                      UInt64(idx)))

# The three key words that must differ. Stated as executable assertions rather than as a comment,
# because "obviously distinct" is how two lanes end up sharing a Philox key.
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_DATAGEN_COUNTER)) !=
        (UInt64(P13_DEV_SEED) ⊻ P13_SALT) "the datagen key word collides with the p13_rng family"
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_DATAGEN_COUNTER)) !=
        (UInt64(P13_FIX_SEED) ⊻ P13_SALT ⊻ UInt64(P13_FIXTURE_COUNTER)) "the datagen key word collides with the fixture datagen family"

# =============================================================================================
# 2. The lambda-hierarchical sampler -- ONE lambda per PAIR, shift CONDITIONAL on it
# =============================================================================================

"""
    p13_sample_imsize(rng::AbstractRNG) -> Tuple{Int,Int}

Draw one image size from the F5 mixture, `P13_IMSIZE_SET` under `P13_IMSIZE_WEIGHTS`.

BOTH ARE **READ** FROM THE FROZEN PRE-REGISTRATION, never retyped -- `consts.jl` reads them in
turn from the frozen amended gate file through its isolated `_GC` module, and
`p13_require_phase11` asserts they equal Phase 11's own `P11_IMSIZE_SET`/`P11_IMSIZE_WEIGHTS`.
F5's binding invariant is train-joint == eval-joint, and a retyped mixture is a mixture that can
silently drift out of agreement with the net it scores.

Hand-rolled inverse-CDF over the fixed categorical, mirroring `sample_imsize`
(`spike/data/seeding.jl:106-114`) and `sample_p11_imsize`
(`spike/data/p11_generate.jl:202-210`); no extra dependency.
"""
function p13_sample_imsize(rng::AbstractRNG)
    r = rand(rng)
    c = 0.0
    @inbounds for i in eachindex(P13_IMSIZE_SET)
        c += P13_IMSIZE_WEIGHTS[i]
        r <= c && return P13_IMSIZE_SET[i]
    end
    return P13_IMSIZE_SET[end]   # float-rounding guard: r just over the weight sum -> last bin
end

"""
    p13_sample_theta_given_lambda(rng::AbstractRNG, lambda::Real) -> NamedTuple

Draw one theta from `pi(theta | lambda)`: the registration error CONDITIONAL on the supplied
uncertainty level, everything else from its own unchanged prior object.

MIRRORS `sample_p11_theta` (`spike/data/p11_generate.jl:153-189`) IN DRAW ORDER AND IN JOINT,
with ONE difference that is the entire point: `lambda` is an ARGUMENT here, not a draw. The
Phase-11 pool drew lambda per SAMPLE because a Phase-11 sample is one acquisition; a Phase-13
item is a sample stack AND its control stack, and D-03 fixes ONE lambda for the pair, so the
level is drawn once by the caller and passed to both halves.

    shift_dx, shift_dy ~ Uniform(-lambda, +lambda)   <- CONDITIONAL on the level
    rho_true            = ghat(mu*), mu* ~ MU_PRIOR  <- unchanged
    the other five nuisances ~ their existing priors <- unchanged, chromatic_eps LAST

`Uniform(-lambda, lambda)` is the HALF-WIDTH parameterisation, never a Gaussian of scale lambda,
so `lambda = LAMBDA_MAX` is exactly the stated `SHIFT_PRIOR` and the widest rung IS the prior.

WHY THE CONDITIONAL DRAW MATTERS EVEN THOUGH PHASE 13 DOES NOT INFER THE SHIFT. The Phase-11 net
reads the conditioning row as the level the shift was drawn under. Drawing the shift from the
MARGINAL `Uniform(-3, 3)` while telling the net `lambda = 0.25` would feed it a conditioning
value that is a lie about the same column's own contents -- a covariate shift against the very
input row D-03 exists to exercise, and nothing would throw.

The returned NamedTuple is assembled in `sample_prior`'s field order, so `collect(values(theta))`
stays row-aligned with `p11_theta_prior_bounds()` and with every positional theta read.
"""
function p13_sample_theta_given_lambda(rng::AbstractRNG, lambda::Real)
    lambda > 0 || throw(ArgumentError(
        "p13_sample_theta_given_lambda: lambda must be > 0 (got $lambda); a zero level would " *
        "collapse the registration prior to a point mass, which is not a rung of the ladder"))
    shift_dx = rand(rng, Uniform(-lambda, lambda))
    shift_dy = rand(rng, Uniform(-lambda, lambda))
    μ_star           = rand(rng, MU_PRIOR)
    spillover        = rand(rng, SPILLOVER_PRIOR)
    autofluorescence = rand(rng, AUTOFLUORESCENCE_PRIOR)
    label_efficiency = rand(rng, LABEL_EFFICIENCY_PRIOR)
    noise            = rand(rng, NOISE_PRIOR)
    chromatic_eps    = rand(rng, CHROMATIC_PRIOR)   # fixed prior, NOT lambda-scaled (P11 D-09)

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

    # The cheap structural guard Phase 11 wrote for the same reason (its Pitfall 1 / R4): if a
    # future edit ever stops drawing the shift CONDITIONAL on lambda, this fires on the FIRST
    # item instead of surfacing hours later as an unexplained flat lambda response.
    @assert abs(theta.shift_dx) <= lambda && abs(theta.shift_dy) <= lambda "p13_sample_theta_given_lambda: |shift| exceeds lambda ($(theta.shift_dx), $(theta.shift_dy) vs $lambda) -- the shift is no longer drawn CONDITIONAL on lambda"

    return theta
end

"""
    p13_draw_labelled_item(idx::Integer; tau = p13_tau(), imsize_sampler = p13_sample_imsize,
                           rng_for = p13_datagen_rng) -> NamedTuple

The CHEAP half of one training item: everything that does NOT require a simulation. Returns
`(rng, lambda, theta_s, theta_c, imsize, rho_s, rho_c, class, idx)`.

THIS EXISTS BECAUSE THE LABEL IS A FUNCTION OF THETA ALONE. `three_way_label` reads only
`rho_s` and `rho_c`, so the whole-class accept/reject decision of `stratify_by_class` can be made
here, before any image exists. A rejected item costs two prior draws. `p13_simulate_labelled`
calls this function and then CONTINUES FROM THE RETURNED `rng`, so the simulated item is
bit-identical to what a single-pass implementation would have produced -- the split is an
evaluation order, not a different generative story.

`rng` IS RETURNED ON PURPOSE. Rebuilding it and re-drawing would also be deterministic, but
returning it makes it impossible for the prescreen and the simulation to drift apart.
"""
function p13_draw_labelled_item(idx::Integer; tau::Real = p13_tau(),
                                imsize_sampler = p13_sample_imsize,
                                rng_for = p13_datagen_rng)
    rng     = rng_for(idx)
    lambda  = rand(rng, Uniform(LAMBDA_MIN, LAMBDA_MAX))   # the LEVEL, drawn FIRST, once per PAIR
    theta_s = p13_sample_theta_given_lambda(rng, lambda)
    theta_c = p13_sample_theta_given_lambda(rng, lambda)
    isz     = imsize_sampler(rng)
    rho_s   = theta_s.ρ_true
    rho_c   = theta_c.ρ_true
    return (rng = rng, lambda = lambda, theta_s = theta_s, theta_c = theta_c, imsize = isz,
            rho_s = rho_s, rho_c = rho_c,
            class = three_way_label(rho_s, rho_c; tau = tau), idx = Int(idx))
end

"""
    p13_simulate_labelled(idx::Integer, basis; tau = p13_tau(),
                          imsize_sampler = p13_sample_imsize,
                          rng_for = p13_datagen_rng) -> NamedTuple

ONE labelled training item, complete: `(Zs, Zc, lambda, rho_s, rho_c, class, imsize, idx)`.

    rng     = rng_for(idx)                       # counter-based, keyed by the GLOBAL INDEX
    lambda ~ Uniform(LAMBDA_MIN, LAMBDA_MAX)     # the LEVEL, drawn FIRST
    theta_s, theta_c ~ pi(theta | lambda)        # BOTH under the SAME lambda
    imsize ~ the F5 mixture                      # ONE size for the acquisition
    Zs, Zc  = standardize_summary(encode_d01(patch_summary(build_mci(simulate_pair(...)))),
                                  basis.zt, basis.variant)
    class   = three_way_label(rho_s, rho_c; tau = tau)

SAMPLE AND CONTROL SHARE ONE LAMBDA AND ONE IMAGE SIZE because they are two channels of ONE
acquisition: the user states a single registration uncertainty for the instrument and the field
size is a property of that instrument, not of the specimen. Two lambdas would encode the belief
that the user holds different beliefs about the alignment of the two stacks, which is not a claim
the read surface expresses and not a claim any acquisition supports (D-03).

THE LABEL COMES FROM SIMULATOR GROUND TRUTH ONLY -- NEVER FROM A NETWORK OUTPUT. It is
`three_way_label(theta_s.rho_true, theta_c.rho_true)`, read off the drawn theta. A label derived
from a posterior read, an estimated rho-hat, or the net's own logit would make the training
target a function of the model being trained; the class boundary would then move with the model
and the reported evidence would be circular. `tau` is the summary's MEASURED resolution (D-06)
and defaults through `p13_tau()`, which ERRORS if tau has not been frozen -- so no labelled datum
can exist before tau is measured.

BOTH summaries are put into the FROZEN Phase-11 input space through `standardize_summary` with
`basis.zt`. NOTHING IS FITTED HERE (banner (a)). `basis` must be the memoized handle from
`load_p13_basis()`; passing a freshly loaded one is a different object and fails
`assert_frozen_zt` on purpose.
"""
function p13_simulate_labelled(idx::Integer, basis; tau::Real = p13_tau(),
                               imsize_sampler = p13_sample_imsize,
                               rng_for = p13_datagen_rng)
    d   = p13_draw_labelled_item(idx; tau = tau, imsize_sampler = imsize_sampler,
                                 rng_for = rng_for)
    rng = d.rng
    Zs  = standardize_summary(
              encode_d01(patch_summary(build_mci(
                  simulate_pair(rng, d.theta_s; imsize = d.imsize)))), basis.zt, basis.variant)
    Zc  = standardize_summary(
              encode_d01(patch_summary(build_mci(
                  simulate_pair(rng, d.theta_c; imsize = d.imsize)))), basis.zt, basis.variant)
    return (Zs = Zs, Zc = Zc, lambda = d.lambda, rho_s = d.rho_s, rho_c = d.rho_c,
            class = d.class, imsize = d.imsize, idx = d.idx)
end

# =============================================================================================
# 3. Whole-class accept/reject stratification (D-07-i)
# =============================================================================================

"""
    p13_class_quota(n::Integer, target = P13_TARGET_CLASS_FREQ) -> NamedTuple

Integer per-class item counts summing EXACTLY to `n`, by largest remainder in the fixed order
`(exclusion, random, coloc)`.

The fixed tie-break order is what makes the quota a deterministic function of `n` alone. `target`
holds `Rational`s, so the products are exact and the remainder distribution is not a floating
point coin flip.
"""
function p13_class_quota(n::Integer, target = P13_TARGET_CLASS_FREQ)
    n >= 1 || throw(ArgumentError("p13_class_quota: n must be >= 1 (got $n)"))
    keys_in_order = (:exclusion, :random, :coloc)
    raw   = [n * getproperty(target, k) for k in keys_in_order]        # exact Rationals
    base  = [floor(Int, r) for r in raw]
    short = n - sum(base)
    # Largest remainder first; ties broken by the fixed key order above (`sortperm` is stable).
    order = sortperm([-(raw[i] - base[i]) for i in eachindex(base)])
    for j in 1:short
        base[order[((j - 1) % length(base)) + 1]] += 1
    end
    @assert sum(base) == n "p13_class_quota: quota $(base) does not sum to $n"
    return (exclusion = base[1], random = base[2], coloc = base[3])
end

"""
    p13_select_indices(n::Integer; target = P13_TARGET_CLASS_FREQ, tau = p13_tau(),
                       imsize_sampler = p13_sample_imsize, rng_for = p13_datagen_rng,
                       first_index = 1, max_draw_factor = P13_STRAT_MAX_DRAW_FACTOR,
                       shuffle = true)
        -> (accepted, classes, n_draws, quota)

THE whole-class accept/reject search, and the ONLY implementation of it. Walks
`first_index, first_index+1, ...` in order, draws the CHEAP theta-only item at each index, and
ACCEPTS the whole item iff its class has quota left. Returns the accepted global indices, their
classes in the same order, how many indices were examined, and the integer quota.

SERIAL BY NECESSITY, NOT BY OVERSIGHT. The decision at index `i` depends on which quotas were
already full when `i` was examined, so the search cannot be threaded without changing the answer.
It is also cheap -- the label is a function of theta alone, so no image is simulated here.

BOTH `stratify_by_class` AND `p13_generate_pool` CALL THIS. Two copies of an accept/reject rule is
two chances for the pool on disk to disagree with the pool a test builds, so there is one.

=================================================================================================
WHY THE ACCEPTED ORDER IS PERMUTED BEFORE IT IS RETURNED (`shuffle = true`) -- MEASURED, NOT
PRECAUTIONARY
=================================================================================================
ACCEPT/REJECT DESTROYS EXCHANGEABILITY ACROSS COLUMNS, and every downstream consumer of this list
assumes it. The classes fill their quotas at different times because pi does not produce them at
equal rates: at tau = 0.15 the masses are E 0.3119 / R 0.4713 / C 0.2168, so RANDOM fills first,
then EXCLUSION, and COLOC -- the rarest -- is still being accepted long after both others have
stopped. The TAIL of the accepted list is therefore nearly pure COLOC.

MEASURED at n = 48000, tau = 0.15, walking from index 1: the contiguous `P13_VAL_FRAC = 0.15`
validation split (the last 7200 of 48000) came out **E 1266 / R 0 / C 5934** -- with ZERO RANDOM
items. RANDOM is the SHARED NEGATIVE CLASS OF BOTH HEADS (D-11), so a validation set without it
gives both heads positives and no negatives: the validation risk would then be minimised by
predicting "positive" everywhere, and it is that risk which drives early stopping and the
best-checkpoint rule. The trained net would have been selected against a degenerate signal.

The analog's comment -- "contiguous, the stream is already i.i.d. across columns"
(`spike/validation/train_ratio.jl:296`) -- is TRUE THERE and FALSE HERE, for exactly the reason
above: that pool is unstratified. Rather than special-case the split, the ORDER is repaired at
the source, so a contiguous split, a shard boundary and a `head_log_odds` measured on a training
subset are all valid again for the same reason they are valid in the analog.

THE PERMUTATION CHANGES NO DISTRIBUTION. It reorders columns; it does not add, drop, reweight or
reshape a single item, so `q(theta | class)` is still exactly `pi(theta | class)` and the F2
exactness condition is untouched. It is drawn from the SAME counter-based family as the pool
(`rng_for(0)` -- index 0 is never a sample index, since the walk starts at `first_index >= 1`),
so it is a pure function of the seed and reproduces bitwise on any thread count and on any resume.
"""
function p13_select_indices(n::Integer; target = P13_TARGET_CLASS_FREQ, tau::Real = p13_tau(),
                            imsize_sampler = p13_sample_imsize, rng_for = p13_datagen_rng,
                            first_index::Integer = 1,
                            max_draw_factor::Integer = P13_STRAT_MAX_DRAW_FACTOR,
                            shuffle::Bool = true)
    n >= 1 || throw(ArgumentError("p13_select_indices: n must be >= 1 (got $n)"))
    quota  = p13_class_quota(n, target)
    want   = Dict(EXCLUSION => quota.exclusion, RANDOM => quota.random, COLOC => quota.coloc)
    filled = Dict(EXCLUSION => 0, RANDOM => 0, COLOC => 0)

    accepted = Int[];            sizehint!(accepted, n)
    classes  = ThreeWayClass[];  sizehint!(classes, n)
    idx       = Int(first_index) - 1
    n_draws   = 0
    max_draws = max_draw_factor * n
    while length(accepted) < n
        idx     += 1
        n_draws += 1
        n_draws <= max_draws || error("""
            p13_select_indices: examined $n_draws draws without filling the quota $(quota) for
            n = $n (filled so far: exclusion $(filled[EXCLUSION]), random $(filled[RANDOM]),
            coloc $(filled[COLOC])). The bound is max_draw_factor = $max_draw_factor times n,
            against an EXPECTED overhead of about 1.54x at tau = 0.15. Either `target` asks for a
            class pi almost never produces, or the label rule is degenerate at tau = $tau. This is
            a diagnosis, not a budget to raise.""")
        d = p13_draw_labelled_item(idx; tau = tau, imsize_sampler = imsize_sampler,
                                   rng_for = rng_for)
        filled[d.class] < want[d.class] || continue          # WHOLE-ITEM reject; theta untouched
        filled[d.class] += 1
        push!(accepted, d.idx)
        push!(classes, d.class)
    end

    # Repair exchangeability across columns -- see the docstring's measured E/R/C = 1266/0/5934
    # finding. `rng_for(0)` rides the SAME counter-based family as the pool, and index 0 is never
    # a sample index, so the permutation is reproducible and cannot collide with an item stream.
    if shuffle
        perm     = randperm(rng_for(0), length(accepted))
        accepted = accepted[perm]
        classes  = classes[perm]
    end

    return (accepted = accepted, classes = classes, n_draws = n_draws, quota = quota)
end

"""
    p13_simulate_indices(indices, basis; tau, imsize_sampler, rng_for, parallel) -> Vector

Simulate the accepted items for `indices`. Threaded, because each item is a pure function of its
global index -- no shared mutable RNG, no locks -- so the parallel and serial paths are
BYTE-IDENTICAL for any thread count.
"""
function p13_simulate_indices(indices, basis; tau::Real = p13_tau(),
                              imsize_sampler = p13_sample_imsize,
                              rng_for = p13_datagen_rng,
                              parallel::Bool = Threads.nthreads() > 1)
    items = Vector{Any}(undef, length(indices))
    fill_one! = function (j::Int)
        items[j] = p13_simulate_labelled(indices[j], basis; tau = tau,
                                         imsize_sampler = imsize_sampler, rng_for = rng_for)
        return nothing
    end
    if parallel
        Threads.@threads for j in eachindex(indices)
            fill_one!(j)
        end
    else
        for j in eachindex(indices)
            fill_one!(j)
        end
    end
    return items
end

"""
    p13_assert_classes_agree(items, classes, indices)

The prescreen and the simulation must agree on every label. If they ever disagree the pool is
silently mislabelled, which is undetectable downstream -- so it is an outright error here.
"""
function p13_assert_classes_agree(items, classes, indices)
    for j in eachindex(items)
        items[j].class === classes[j] || error(
            "p13: item $(indices[j]) simulated as $(items[j].class) but was accepted as " *
            "$(classes[j]) -- the whole-class prescreen and the simulation disagree")
    end
    return nothing
end

"""
    stratify_by_class(n::Integer, basis; target = P13_TARGET_CLASS_FREQ, tau = p13_tau(),
                      imsize_sampler = p13_sample_imsize, rng_for = p13_datagen_rng,
                      first_index = 1, parallel = Threads.nthreads() > 1,
                      max_draw_factor = P13_STRAT_MAX_DRAW_FACTOR, verbose = false)
        -> (items, realized, n_draws, accepted_indices)

Generate `n` labelled items whose CLASS FREQUENCIES match `target`, by ACCEPTING OR REJECTING
WHOLE ITEMS and never reshaping anything within a class.

THIS IS DESIGN D-07-i (`P13_STRATIFICATION = :class_frequency`), and the accept/reject is at
CLASS GRANULARITY PRECISELY SO `q(theta | class)` STAYS `pi(theta | class)`. That equality is the
EXACTNESS CONDITION for the scalar per-head correction `log(q_pos/q_neg)` (13-RESEARCH F2): the
scalar recovers the pi-scale Bayes factor if and only if stratification changed the class
FREQUENCIES only. A within-class change (for instance stratifying rho across the exclusion range,
which is D-07's own literal wording) enters INSIDE the evidence integral and depends on Z, so no
scalar can repair it -- plan 13-07's NEGATIVE CONTROL 2 demonstrated exactly that failure on a
closed-form toy rather than asserting it. So: an item is kept or dropped as a whole, on the sole
basis of its class, and its theta is never touched.

THE EXPECTED REJECTION OVERHEAD, SO A SLOW RUN IS NOT MISTAKEN FOR A HANG. Under pi the three
classes are NOT equally likely; the measured masses (13-RESEARCH E1, variant (a)) are

    tau = 0.05   E 0.3850 / R 0.3413 / C 0.2738
    tau = 0.10   E 0.3482 / R 0.4108 / C 0.2410
    tau = 0.15   E 0.3119 / R 0.4713 / C 0.2168

so with an equal-thirds target the BINDING class is COLOC and the expected number of DRAWS is
about `(n/3) / q_C`, i.e. **~1.54 x n at tau = 0.15** (~1.38x at 0.10, ~1.22x at 0.05). Roughly a
third of all draws are rejected, and the rejected ones are cheap: the label is a function of
theta alone, so a rejection costs two `sample_prior` calls and NO simulation. Only the accepted
`n` items are simulated. `max_draw_factor` bounds the search at
`max_draw_factor * n` draws and errors with a diagnosis, so a mis-specified target or a broken
label rule fails loudly instead of spinning.

DETERMINISM. The acceptance decision for index `i` depends on which quotas were already full when
`i` was examined, so the SEARCH IS DELIBERATELY SERIAL and walks `first_index, first_index+1, ...`
in order. The SIMULATION of the accepted indices is order-independent (each item is a pure
function of its index) and is therefore threaded. Same accepted index list, same `Z` matrix, any
thread count.

Returns the items in ACCEPTED-INDEX order together with `realized` (a `class_masses` NamedTuple),
`n_draws` (how many indices were examined) and `accepted_indices`.
"""
function stratify_by_class(n::Integer, basis; target = P13_TARGET_CLASS_FREQ,
                           tau::Real = p13_tau(),
                           imsize_sampler = p13_sample_imsize,
                           rng_for = p13_datagen_rng,
                           first_index::Integer = 1,
                           parallel::Bool = Threads.nthreads() > 1,
                           max_draw_factor::Integer = P13_STRAT_MAX_DRAW_FACTOR,
                           verbose::Bool = false)
    # --- PHASE A (serial, cheap): whole-item accept/reject on the theta-only label -----------
    sel = p13_select_indices(n; target = target, tau = tau, imsize_sampler = imsize_sampler,
                             rng_for = rng_for, first_index = first_index,
                             max_draw_factor = max_draw_factor)

    verbose && println("[stratify] $(length(sel.accepted)) items accepted from " *
                       "$(sel.n_draws) draws " *
                       "(overhead $(round(sel.n_draws / n; digits = 3))x); quota $(sel.quota)")

    # --- PHASE B (threaded): simulate ONLY the accepted items --------------------------------
    items = p13_simulate_indices(sel.accepted, basis; tau = tau,
                                 imsize_sampler = imsize_sampler, rng_for = rng_for,
                                 parallel = parallel)
    p13_assert_classes_agree(items, sel.classes, sel.accepted)

    return (items = items,
            realized = class_masses(sel.classes),
            n_draws = sel.n_draws,
            accepted_indices = sel.accepted)
end

# =============================================================================================
# 4. Conditioned pair assembly (D-03, D-11)
# =============================================================================================

"""
    assemble_conditioned_pairs(items) -> NamedTuple

Assemble the net's training input and target from a vector of labelled items:

    (Z = three_way_input_dim(G, n_cond) x n,
     targets = 4 x n,
     classes = n-vector of ThreeWayClass,
     lambdas = n-vector of Float64)

Each column of `Z` is `p13_encode_pair(item.Zs, item.Zc, item.lambda)` -- ONE call, so LAMBDA IS
APPENDED EXACTLY ONCE AND EXACTLY AFTER THE PAIR ENCODING (D-03,
`P13_LAMBDA_PLACEMENT = :append_after_pair_encode`). Each column of `targets` is
`head_targets(item.class)` via `target_matrix`, rows `[y_C; w_C; y_E; w_E]`.

TWO ASSERTIONS, ON ENTRY AND ON EXIT, BOTH ArgumentErrors:

  * ON ENTRY every `Zs` and `Zc` is `2*G^2` rows -- EQUAL length across the pair and across the
    pool, and EVEN. An odd length is the signature of lambda folded INTO the summary (the wrong
    form): `p13_pair_encode`'s own `iseven` guard would catch it, but catching it here names the
    cause instead of reporting an arithmetic surprise 48000 columns later.
  * ON EXIT `size(Z, 1) == three_way_input_dim(G, p13_conditioning_length())` and
    `size(targets, 1) == 4`. The exit width check is what catches DOUBLE-APPENDED conditioning,
    which is otherwise SILENT -- it produces a taller matrix that trains perfectly happily.

`G` is DERIVED from the summary length (`G = isqrt(length(Zs) / 2)`), never assumed to be 8, and
the realized width is never written as a literal.
"""
function assemble_conditioned_pairs(items)
    isempty(items) && throw(ArgumentError(
        "assemble_conditioned_pairs: items must be non-empty (a 4x0 target is never what the " *
        "caller meant)"))
    n   = length(items)
    len = length(items[1].Zs)
    for (j, it) in enumerate(items)
        length(it.Zs) == len && length(it.Zc) == len || throw(ArgumentError(
            "assemble_conditioned_pairs: item $j has summary lengths " *
            "($(length(it.Zs)), $(length(it.Zc))) but the pool is $len rows -- both halves of " *
            "every pair must be the same 2*G^2 rows"))
    end
    iseven(len) || throw(ArgumentError(
        "assemble_conditioned_pairs: summary length $len must be 2*G^2 (even). An ODD length is " *
        "the signature of the conditioning value folded INTO the summary instead of appended " *
        "after the pair encoding (D-03)."))
    G    = isqrt(len ÷ 2)
    2 * G^2 == len || throw(ArgumentError(
        "assemble_conditioned_pairs: summary length $len is not 2*G^2 for any integer G"))

    n_cond   = p13_conditioning_length()
    expected = three_way_input_dim(G, n_cond)
    Z        = Matrix{Float32}(undef, expected, n)
    lambdas  = Vector{Float64}(undef, n)
    classes  = Vector{ThreeWayClass}(undef, n)
    for j in 1:n
        it          = items[j]
        Z[:, j]     = p13_encode_pair(it.Zs, it.Zc, it.lambda)
        lambdas[j]  = Float64(it.lambda)
        classes[j]  = it.class
    end
    targets = target_matrix(classes)

    size(Z, 1) == expected || throw(ArgumentError(
        "assemble_conditioned_pairs: built width $(size(Z, 1)) but " *
        "three_way_input_dim($G, $n_cond) is $expected -- the conditioning vector is misplaced " *
        "or double-counted (D-03)"))
    size(targets, 1) == 4 || throw(ArgumentError(
        "assemble_conditioned_pairs: targets have $(size(targets, 1)) rows, not the 4 rows " *
        "[y_C; w_C; y_E; w_E] the masked two-head loss indexes positionally (D-11)"))

    return (Z = Z, targets = targets, classes = classes, lambdas = lambdas)
end

"""
    check_frozen_zt(items, basis) -> NamedTuple

Run `assert_frozen_zt` on the ASSEMBLED STANDARDIZED pool (both halves of every pair, so `2n`
columns) and return its verdict together with the measured moments, so the Pitfall-5 warning sign
is a REPORTED NUMBER rather than a hope.

THE HARD HALF IS UNAMBIGUOUS. `assert_frozen_zt` first requires `basis.zt === load_p13_basis().zt`
-- object identity against the memoized single load. A caller who fitted their own standardizer
and packed it into a look-alike handle ERRORS here. That is the check that actually enforces D-02.

THE SOFT HALF IS A HEURISTIC, AND ITS POWER IS LIMITED **HERE IN PARTICULAR**. The heuristic's
premise is that a transform frozen on the PHASE-11 pool, applied to a DIFFERENT pool, should not
reproduce mean ~= 0 and sd ~= 1. But D-02 requires the Phase-13 pool to be drawn from the SAME
joint the Phase-11 net trained on -- same priors, same lambda hierarchy, same F5 image-size
mixture (that IS the F5 train-joint == eval-joint invariant) -- so the per-column marginal is by
construction the same one `zt` was fitted against. At large `n` the standardized moments may
therefore land close to 0 and 1 LEGITIMATELY, and the verdict may read `:refit_suspected` while
nothing was re-fit. THAT READING IS REPORTED, NOT SUPPRESSED: it is a limitation of the
heuristic, and the honest handling is to say so next to the number. What would be evidence of an
actual re-fit is the HARD identity failing, or a z-score transform being CONSTRUCTED anywhere in
the Phase-13 data path -- both of which are checked, the second by source grep in
`spike/test/test_p13_datagen.jl`.

`:indeterminate` is returned below `P13_FROZEN_ZT_MIN_COLS` columns, where the per-row moments are
too noisy to judge, rather than a false clean bill.
"""
function check_frozen_zt(items, basis)
    isempty(items) && throw(ArgumentError("check_frozen_zt: items must be non-empty"))
    len = length(items[1].Zs)
    M   = Matrix{Float32}(undef, len, 2 * length(items))
    for (j, it) in enumerate(items)
        M[:, 2j - 1] = it.Zs
        M[:, 2j]     = it.Zc
    end
    return assert_frozen_zt(basis, M)
end

# =============================================================================================
# 5. Persistence -- content-hash-named, sharded, atomic, resumable
# =============================================================================================

"""
    p13_pool_config(n; shard_size, tau, target, imsize_tag) -> NamedTuple

The canonical generating config the content hash names the pool directory by. Carries every field
whose change MUST invalidate the pool: the size and shard size, the seed/salt/counter that define
the draw stream, the tau in force and the frozen cut variant (both of which change every label),
the target class frequencies and the stratification design, the lambda range (the hierarchical
joint itself), the F5 image-size mixture, the theta arity, the summary width, the conditioning
length and the schema tag.

`cache_hash` folds these over the RAW BYTES of the frozen simulator/summary sources as well
(`HASH_SRC_FILES`, which includes the two frozen `src/` files as a decoupling-breach tripwire), so
an edited simulator resolves to a DIFFERENT directory instead of silently reusing a stale pool.
"""
function p13_pool_config(n::Integer; shard_size::Integer = P13_POOL_SHARD_SIZE,
                         tau::Real = p13_tau(), target = P13_TARGET_CLASS_FREQ,
                         imsize_tag::Symbol = :f5_mixture)
    return (
        N                = Int(n),
        shard_size       = Int(shard_size),
        master_seed      = UInt64(P13_DEV_SEED),
        salt             = UInt64(P13_SALT),
        datagen_counter  = Int(P13_DATAGEN_COUNTER),
        tau              = Float64(tau),
        cut_variant      = P13_CUT_VARIANT,
        stratification   = P13_STRATIFICATION,
        target_exclusion = string(target.exclusion),
        target_random    = string(target.random),
        target_coloc     = string(target.coloc),
        lambda_min       = Float64(LAMBDA_MIN),
        lambda_max       = Float64(LAMBDA_MAX),
        imsize_set       = P13_IMSIZE_SET,
        imsize_weights   = P13_IMSIZE_WEIGHTS,
        imsize_tag       = imsize_tag,
        theta_dim        = length(sample_prior(Random.Xoshiro(0))),
        summary_dim      = 2 * length(load_p13_basis().zt.mean),
        n_cond           = p13_conditioning_length(),
        schema_version   = P13_POOL_SCHEMA,
    )
end

"""
    p13_pool_dir(n; cache_root, kwargs...) -> String

Resolve (creating if absent) the content-hash-named pool directory for a pool of `n` items,
WITHOUT generating anything. The trainer uses this to find the pool the datagen step wrote, and
`open_or_invalidate` re-checks the stored manifest hash on reopen rather than trusting the name.
"""
p13_pool_dir(n::Integer; cache_root = P13_POOL_CACHE_ROOT, kwargs...) =
    open_or_invalidate(cache_root, p13_pool_config(n; kwargs...))

"""
    p13_items_to_pool(items) -> NamedTuple

Pack a vector of labelled items into column-major buffers (one item = one column, the
NeuralEstimators d-by-K convention): `(Zs, Zc, lambda, rho_s, rho_c, class, global_index,
imsize)`. `class` is stored as `Int` rather than as the `@enum` value, so reading a shard back
never depends on `ThreeWayClass` being reconstructible by the deserializer.
"""
function p13_items_to_pool(items)
    isempty(items) && throw(ArgumentError("p13_items_to_pool: items must be non-empty"))
    n   = length(items)
    len = length(items[1].Zs)
    Zs  = Matrix{Float32}(undef, len, n)
    Zc  = Matrix{Float32}(undef, len, n)
    lam = Vector{Float64}(undef, n)
    rs  = Vector{Float64}(undef, n)
    rc  = Vector{Float64}(undef, n)
    cl  = Vector{Int}(undef, n)
    gi  = Vector{Int}(undef, n)
    isz = Vector{Tuple{Int,Int}}(undef, n)
    for (j, it) in enumerate(items)
        Zs[:, j] = it.Zs
        Zc[:, j] = it.Zc
        lam[j]   = Float64(it.lambda)
        rs[j]    = Float64(it.rho_s)
        rc[j]    = Float64(it.rho_c)
        cl[j]    = Int(it.class)
        gi[j]    = Int(it.idx)
        isz[j]   = it.imsize
    end
    return (Zs = Zs, Zc = Zc, lambda = lam, rho_s = rs, rho_c = rc, class = cl,
            global_index = gi, imsize = isz)
end

"""
    p13_pool_to_items(pool) -> Vector

Unpack column-major pool buffers back into the item NamedTuples the rest of this file consumes.
The `Zs`/`Zc` entries are column VIEWS, so no summary is copied.
"""
function p13_pool_to_items(pool)
    n = length(pool.lambda)
    return [(Zs = view(pool.Zs, :, j), Zc = view(pool.Zc, :, j),
             lambda = pool.lambda[j], rho_s = pool.rho_s[j], rho_c = pool.rho_c[j],
             class = ThreeWayClass(pool.class[j]), imsize = pool.imsize[j],
             idx = pool.global_index[j]) for j in 1:n]
end

"""
    save_pool(path, data; meta = (;)) -> String

Atomically persist one labelled-pool file: write to `path * ".tmp"`, REOPEN to integrity-check
the required keys, then `mv(...; force = true)`. This is the filesystem-atomic idiom of
`save_ratio` (`spike/validation/train_ratio.jl:331-347`) and `write_p11_shard`
(`spike/data/p11_generate.jl:330-341`): a crash mid-write leaves only a discardable `.tmp`, never
a torn file that a later run would happily read as complete.

`data` may be either a vector of items or the column-major buffers `p13_items_to_pool` produces.

PERSISTED PROVENANCE, so a pool can be tied to the exact conditions it was produced under:
`schema_version`, the realized class frequencies, the `tau` in force, the Phase-11 artifact path
and its sha256 when it is on disk, the seed / salt / counter, the git blob sha of
`spike/p13/consts.jl`, and a UTC `generated` stamp. Anything the caller adds through `meta` is
merged on top.
"""
function save_pool(path, data; meta = (;))
    pool = data isa NamedTuple ? data : p13_items_to_pool(data)
    d    = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp;
        schema_version  = P13_POOL_SCHEMA,
        Zs              = pool.Zs,
        Zc              = pool.Zc,
        lambda          = pool.lambda,
        rho_s           = pool.rho_s,
        rho_c           = pool.rho_c,
        class           = pool.class,
        global_index    = pool.global_index,
        imsize          = pool.imsize,
        realized_freq   = class_masses(ThreeWayClass.(pool.class)),
        tau             = Float64(p13_tau()),
        cut_variant     = P13_CUT_VARIANT,
        stratification  = P13_STRATIFICATION,
        phase11_net     = P13_PHASE11_NET,
        phase11_sha     = p13_phase11_sha(),
        master_seed     = UInt64(P13_DEV_SEED),
        salt            = UInt64(P13_SALT),
        datagen_counter = Int(P13_DATAGEN_COUNTER),
        consts_sha      = p13_consts_sha(),
        meta            = merge((generated = string(Dates.now(Dates.UTC)) * "Z",), meta))
    JLD2.jldopen(tmp, "r") do f
        for k in ("Zs", "Zc", "class", "lambda", "global_index", "realized_freq")
            @assert haskey(f, k) "save_pool: integrity check failed, $tmp missing $k"
        end
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end

"""
    load_pool(path) -> NamedTuple

Paired reader for `save_pool`. Errors when the file is absent, when a required key is missing, or
when the stored schema version is not the one this file writes -- a silently accepted stale layout
would train a head against labels whose semantics are not the ones documented here. Returns the
column-major buffers plus the persisted provenance block.
"""
function load_pool(path)
    isfile(path) || error("load_pool: no pool file at $path")
    d = JLD2.load(path)
    for k in ("Zs", "Zc", "lambda", "rho_s", "rho_c", "class", "global_index")
        haskey(d, k) || error("load_pool: $path missing $k key")
    end
    schema = get(d, "schema_version", nothing)
    schema == P13_POOL_SCHEMA ||
        error("load_pool: $path has schema_version $schema, expected $P13_POOL_SCHEMA")
    return (Zs = d["Zs"], Zc = d["Zc"], lambda = d["lambda"],
            rho_s = d["rho_s"], rho_c = d["rho_c"], class = d["class"],
            global_index = d["global_index"], imsize = get(d, "imsize", nothing),
            realized_freq = get(d, "realized_freq", nothing),
            tau = get(d, "tau", nothing), cut_variant = get(d, "cut_variant", nothing),
            stratification = get(d, "stratification", nothing),
            phase11_net = get(d, "phase11_net", nothing),
            phase11_sha = get(d, "phase11_sha", nothing),
            consts_sha = get(d, "consts_sha", nothing),
            schema_version = schema, meta = get(d, "meta", nothing))
end

"""
    p13_phase11_sha() -> Union{String,Nothing}

sha256 of the bound Phase-11 artifact, or `nothing` when it is absent. Recorded in every
persisted pool so a pool can be tied to the exact frozen basis it was standardized against. The
artifact is opened READ-ONLY; nothing here moves, rewrites or prunes it.
"""
p13_phase11_sha() =
    isfile(P13_PHASE11_NET) ? bytes2hex(SHA.sha256(read(P13_PHASE11_NET))) : nothing

"""
    p13_pool_shard_done(path) -> Bool

Resume-by-skip predicate: a Phase-13 pool shard is done iff its FINAL file exists and reopens
with the full key set. Any open/read error (a torn file left by a kill) is reported as not-done so
the shard is regenerated rather than silently trusted.

A Phase-13 shard carries `Zs`/`Zc`/`class` and carries NO `theta`/`summary_aug` column, so
`cache.jl`'s own `_loads_ok` would declare every finished shard NOT done and silently defeat
resume on a multi-hour run. Hence a local predicate, exactly as Phase 11 needed one.
"""
function p13_pool_shard_done(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "Zs") && haskey(f, "Zc") && haskey(f, "class") &&
                haskey(f, "lambda") && haskey(f, "global_index") && haskey(f, "schema_version")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end

"""
    p13_generate_pool(n = P13_TRAIN_N; basis = load_p13_basis(), cache_root, shard_size, tau,
                      target, imsize_sampler, parallel, verbose) -> String

Generate the Phase-13 labelled, class-stratified pool of `n` items into a sharded, atomic,
content-hash-named cache, and return the directory.

THE ORDER OF OPERATIONS, WHICH IS THE WHOLE FILE IN FIVE LINES:

  1. `p13_require_phase11()` through `load_p13_basis()` -- the hard precondition, FIRST;
  2. the whole-class accept/reject search over the theta-only label, SERIAL and cheap;
  3. the accepted indices simulated in shards, threaded inside a shard, SERIAL over shards so
     progress is monotone and legible;
  4. each shard committed atomically and SKIPPED on a re-run (`p13_pool_shard_done`), so a killed
     run RESUMES instead of restarting;
  5. the manifest written last, which is what makes the directory readable by
     `p13_load_pool_dir`.

Because the accepted index list is a pure function of `(n, tau, target, first_index)` and the
draw stream, step 2 is recomputed cheaply on every resume rather than persisted -- there is no
second source of truth to fall out of sync.

NO WALL-CLOCK BLOCKER IS IMPOSED. Phase 11 carried one because its budget was pre-registered;
this plan authorises ONE training run and one pool, and a throwing ceiling here would abort the
phase's only authorised spend on a projection. The per-shard rate and the projected total are
PRINTED instead, so a run that is going to be slow says so early and the decision stays with the
operator.
"""
function p13_generate_pool(n::Integer = P13_TRAIN_N; basis = load_p13_basis(),
                           cache_root = P13_POOL_CACHE_ROOT,
                           shard_size::Integer = P13_POOL_SHARD_SIZE,
                           tau::Real = p13_tau(), target = P13_TARGET_CLASS_FREQ,
                           imsize_sampler = p13_sample_imsize,
                           rng_for = p13_datagen_rng,
                           parallel::Bool = Threads.nthreads() > 1,
                           verbose::Bool = true)
    config  = p13_pool_config(n; shard_size = shard_size, tau = tau, target = target)
    dir     = open_or_invalidate(cache_root, config)
    quota   = p13_class_quota(n, target)
    nshards = cld(n, shard_size)

    if verbose
        println("="^78)
        println("Phase-13 labelled three-way training pool (D-02, D-03, D-07-i, D-11)")
        println("  N = $n   shards = $nshards x $shard_size   threads = $(Threads.nthreads())")
        println("  tau = $tau  (reference lambda = $(P13_TAU_REFERENCE_LAMBDA); tau is tau(lambda))")
        println("  cut = $P13_CUT_VARIANT   stratification = $P13_STRATIFICATION (WHOLE CLASS)")
        println("  target class freq = $target  ->  quota = $quota")
        println("  lambda ~ Uniform($LAMBDA_MIN, $LAMBDA_MAX); shift ~ Uniform(-lambda, lambda); " *
                "ONE lambda per PAIR")
        println("  imsize mixture = $(P13_IMSIZE_SET) w = $(P13_IMSIZE_WEIGHTS)")
        println("  frozen basis = $(basename(P13_PHASE11_NET)) (zt INHERITED, never re-fit)")
        println("  input width = three_way_input_dim(G, $(p13_conditioning_length())) = " *
                "$(three_way_input_dim(isqrt(length(basis.zt.mean)), p13_conditioning_length()))")
        println("  stream = P13_DEV_SEED $(repr(UInt64(P13_DEV_SEED))) xor P13_SALT xor " *
                "P13_DATAGEN_COUNTER $(P13_DATAGEN_COUNTER), keyed per global index")
        println("  dir = $dir")
        println("="^78)
    end

    # --- Step 2: the accepted index list (serial, cheap, recomputed on every resume) ---------
    # RECOMPUTED, NEVER PERSISTED. The list is a pure function of (n, tau, target, first_index)
    # and the draw stream, so re-deriving it on a resume costs seconds and removes a second
    # source of truth that could fall out of sync with the shards on disk.
    t_sel    = time()
    sel      = p13_select_indices(n; target = target, tau = tau,
                                  imsize_sampler = imsize_sampler, rng_for = rng_for)
    accepted = sel.accepted
    classes  = sel.classes
    n_draws  = sel.n_draws
    realized = class_masses(classes)
    verbose && println("[select] $n accepted from $n_draws draws " *
                       "(overhead $(round(n_draws / n; digits = 3))x) in " *
                       "$(round((time() - t_sel) / 60; digits = 2)) min")
    verbose && println("[select] realized class freq = $realized  vs target $target")

    # --- Steps 3-4: simulate in shards, atomically, resuming by skip -------------------------
    t0        = time()
    n_done    = 0
    n_todo    = count(s -> !p13_pool_shard_done(shard_path(dir, s)), 1:nshards)
    n_skipped = nshards - n_todo
    verbose && n_skipped > 0 &&
        println("[resume] $n_skipped of $nshards shards already complete -- skipping them")

    for s in 1:nshards
        path = shard_path(dir, s)
        p13_pool_shard_done(path) && continue               # resume-by-skip
        ts    = time()
        lo    = (s - 1) * shard_size + 1
        hi    = min(s * shard_size, n)
        idxs  = accepted[lo:hi]
        items = p13_simulate_indices(idxs, basis; tau = tau, imsize_sampler = imsize_sampler,
                                     rng_for = rng_for, parallel = parallel)
        p13_assert_classes_agree(items, classes[lo:hi], idxs)
        save_pool(path, items;
                  meta = (shard = s, nshards = nshards, lo = lo, hi = hi,
                          n_draws = n_draws, realized_freq = realized))
        n_done += 1

        elapsed_min = (time() - t0) / 60
        shard_min   = (time() - ts) / 60
        if verbose
            projected = elapsed_min * (n_todo / n_done)
            println("[shard $s/$nshards] $(hi - lo + 1) items in " *
                    "$(round(shard_min; digits = 2)) min (cumulative " *
                    "$(round(elapsed_min; digits = 2)) min, projected total " *
                    "$(round(projected; digits = 1)) min)")
        end
    end

    # --- Step 5: the manifest, last ---------------------------------------------------------
    write_meta(dir, config)
    verbose && println("pool complete in $(round((time() - t0) / 60; digits = 2)) min " *
                       "($n_done generated, $n_skipped resumed) -> $dir")
    return dir
end

"""
    p13_load_pool_dir(dir) -> NamedTuple

Load a generated pool by `hcat`-ing every shard in `dir` in lexical (== numeric, zero-padded)
order. Returns the same column-major shape `load_pool` returns for a single file, with the
provenance block taken from the FIRST shard (every shard is written under one config, and the
directory name is the content hash of that config).
"""
function p13_load_pool_dir(dir)
    isdir(dir) || error("p13_load_pool_dir: no pool directory at $dir")
    files = sort(filter(f -> startswith(basename(f), "shard_") && endswith(f, ".jld2"),
                        readdir(dir; join = true)))
    isempty(files) && error("p13_load_pool_dir: no shard_*.jld2 found in $dir")
    parts = [load_pool(f) for f in files]
    return (Zs = reduce(hcat, (p.Zs for p in parts)),
            Zc = reduce(hcat, (p.Zc for p in parts)),
            lambda = reduce(vcat, (p.lambda for p in parts)),
            rho_s = reduce(vcat, (p.rho_s for p in parts)),
            rho_c = reduce(vcat, (p.rho_c for p in parts)),
            class = reduce(vcat, (p.class for p in parts)),
            global_index = reduce(vcat, (p.global_index for p in parts)),
            imsize = reduce(vcat, (p.imsize for p in parts)),
            tau = parts[1].tau, cut_variant = parts[1].cut_variant,
            stratification = parts[1].stratification,
            phase11_net = parts[1].phase11_net, phase11_sha = parts[1].phase11_sha,
            consts_sha = parts[1].consts_sha,
            schema_version = parts[1].schema_version, nshards = length(files), dir = dir)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so `include`-ing this file (the trainer, the test suite, a later gate) can NEVER
# trigger generation. Run the pool deliberately:
#
#     julia --project=spike -t auto spike/p13/datagen.jl
if abspath(PROGRAM_FILE) == @__FILE__
    basis = load_p13_basis()
    dir   = p13_generate_pool(P13_TRAIN_N; basis = basis, verbose = true)
    pool  = p13_load_pool_dir(dir)
    items = p13_pool_to_items(pool)
    real  = class_masses(ThreeWayClass.(pool.class))
    zt    = check_frozen_zt(items, basis)
    println("="^78)
    println("REALIZED CLASS FREQUENCIES vs TARGET (D-07-i, whole-class accept/reject)")
    for k in (:exclusion, :random, :coloc)
        r = getproperty(real, k)
        t = Float64(getproperty(P13_TARGET_CLASS_FREQ, k))
        println("  $(rpad(k, 10)) realized $(round(r; digits = 6))   target " *
                "$(round(t; digits = 6))   delta $(round(r - t; digits = 6))")
    end
    println("FROZEN-zt CHECK (Pitfall 5): verdict = $(zt.verdict)  " *
            "max|mean| = $(zt.max_abs_mean)  max|sd-1| = $(zt.max_abs_sd_dev)  " *
            "over $(zt.n_cols) columns x $(zt.n_rows) continuous rows")
    println("  read check_frozen_zt's docstring before interpreting the SOFT verdict: the " *
            "Phase-13 pool is drawn from the SAME joint by D-02, which limits the heuristic")
    println("pool dir = $dir")
    println("="^78)
end
