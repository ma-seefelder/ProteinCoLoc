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

# spike/p13/alpha_series.jl --- D-16 semi-synthetic random-to-exclusion series (D-15, D-16).
#
# THE TRAP THAT DOMINATES THIS FILE, STATED IN FULL BEFORE ANY CODE.
# `correlation()` calls `_exclude_zero` (src/colocalization.jl:154-169, 187-203), which DROPS
# any pixel where EITHER channel is 0.0, NaN or missing, and a patch left with <= 15 survivors
# becomes `missing` (summary.jl:44). So the obvious reading of "redistribute ch2 signal out of
# the ch1 object mask" -- set ch2 to zero inside the mask -- does NOT create anti-correlation.
# It DELETES exactly those pixels from the correlation:
#
#   * the measured patch correlation is then computed over the COMPLEMENT ONLY, i.e. over the
#     region the transform did not touch, so it comes back approximately UNCHANGED (random);
#   * object-dense patches can fall below the 15-survivor floor and go `missing`, firing the
#     summary's own mask rows and plausibly the OOD flag;
#   * the ladder therefore reads FLAT, and it reads flat for a reason that has nothing
#     whatsoever to do with segregation.
#
# The same insight is already recorded in the opposite direction by `negctrl_affine`
# (spike/validation/ood.jl:513-528): a blanket `+b` there would turn read-noise zeros INTO
# signal and change which pixels `_exclude_zero` drops, perturbing ~3/64 patch correlations.
# Which pixels are zero is part of the summary's definition, not an implementation detail.
#
# HENCE THE TWO NON-NEGOTIABLE PROPERTIES OF THE TRANSFORM BELOW.
#   (1) A STRICTLY-POSITIVE BACKGROUND FLOOR `b`. Masked ch2 pixels are pushed DOWN TO `b`,
#       never to zero, so the zero SET is untouched and every pixel keeps participating in its
#       patch correlation. The in-repo precedent is `BG_FLOOR = 0.02`
#       (spike/simulator/forward.jl:66): "background stays small-positive, never hard 0.0".
#       HONEST CAVEAT ON THE ENDPOINT: at alpha = 1 the mask region is not EMPTY, it is AT `b`.
#       The endpoint is "ch2 reduced to background wherever ch1 has objects", NOT "ch2 absent".
#       Do not write "fully disjoint" without that qualifier.
#   (2) TOTAL ch2 INTENSITY IS CONSERVED. Whatever is removed from inside the mask is added
#       back over its complement as a single multiplicative boost, so an alpha-response can
#       never be a brightness change wearing a segregation costume.
#
# WHY SPATIALLY AND NOT BY ANTI-CORRELATION (D-16). The simulator's own "exclusion" mechanism
# IS negative intensity correlation across patches (spike/simulator/forward.jl stage 1, the
# sign(rho) flip on channel 2's shared latent component). Constructing the ladder that way
# would test the net against its own training assumption -- circular, and it would leave the
# actual question untouched. This ladder instead constructs DISJOINT SPATIAL LOCALIZATION: ch2
# is moved out of where ch1's objects are. That is what a biologist means by "mutually
# exclusive", and the gap between the two notions is the only thing this series exists to probe.
# Displacement of ch2 objects is likewise rejected: a global displacement is precisely the
# mis-registration nuisance Phase 11 marginalizes over, so a null result would be
# uninterpretable and a positive result unattributable.
#
# ONE CODE PATH, TWO SUBSTRATES (D-15). The interface boundary is `Vector{Matrix{Float64}}`,
# which `simulate_pair` (spike/simulator/forward.jl:103) and `load_tiff.(paths)`
# (src/LoadImages.jl:433) both already emit. `alpha_segregate` is typed on
# `AbstractMatrix{Float64}` and knows NOTHING about where its input came from -- there is no
# substrate branch anywhere in this file, because a provenance branch would make the two arms
# different experiments and their rungs incomparable. The substrate-dependent knobs (the
# background floor, the frame size, what alpha = 0 MEANS, and the declared dynamic range the
# value bound is stated in) are supplied by the CALLER as provenance, never decided here.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. `src/` is reached ONLY read-only,
# transitively through spike/contract.jl's include chain, for the already-shipped
# `_calculate_mask`/`otsu_thresholds` mask rule and the unchanged patch/correlation summary. No
# src/ byte, no manifest byte and no `corpus/` path is touched; the sealed holdout is never
# opened by this phase.

using Statistics     # quantile (the background floor), mean (the mask fraction)

# ORDER MATTERS. The pre-registration first (P13_ALPHA_LADDER and the frozen bounds), then
# spike/contract.jl -- which puts `using StatsBase; using Statistics; using Images` in scope
# BEFORE it read-only-includes src/LoadImages.jl and src/colocalization.jl, the load-order
# landmine documented at spike/contract.jl:34-39 -- then the 128-row encoding the m-bar
# statistic is read out of. Guarded for idempotency under a runner, the house include idiom.
isdefined(@__MODULE__, :P13_ALPHA_LADDER) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :build_mci)        || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)       || include(joinpath(@__DIR__, "..", "data", "encode.jl"))

# --- (1) The strictly-positive background floor ----------------------------------------------

"""
    alpha_background_floor(y::AbstractMatrix; q = P13_ALPHA_BG_QUANTILE) -> Float64

The strictly-positive background floor `b` masked ch2 pixels are pushed DOWN TO (never to
zero), computed as `quantile(vec(y), q)` with `q = P13_ALPHA_BG_QUANTILE` = 0.05.

The in-repo precedent is `const BG_FLOOR = 0.02` (spike/simulator/forward.jl:66), commented
there as "background stays small-positive, never hard 0.0" for exactly this reason: a hard 0.0
would let `_exclude_zero` (src/colocalization.jl:154-169) thin patches below the 15-survivor
floor. This is its data-driven analogue, so the same rule serves both substrates.

Throws an `ArgumentError` naming `_exclude_zero` when the result is not strictly positive. That
is a runtime guard rather than a comment because it is REACHABLE: an image with more than
`100q` per cent exact-zero pixels silently yields `b = 0` and re-opens the trap in full. Not
observed on the committed real fixtures (b = 0.00624 / 0.00591, both far below the Otsu
threshold), which is the point -- it must fail loudly on the first image where it is not true.
"""
function alpha_background_floor(y::AbstractMatrix; q::Real = P13_ALPHA_BG_QUANTILE)
    (0.0 < q < 1.0) ||
        throw(ArgumentError("alpha_background_floor: q must be in (0, 1), got $q"))
    b = quantile(vec(Float64.(y)), q)
    b > 0.0 || throw(ArgumentError(
        "alpha_background_floor: the background floor must be STRICTLY POSITIVE, got b = $b " *
        "at q = $q. A floor of 0.0 re-opens the _exclude_zero trap " *
        "(src/colocalization.jl:154-169): _exclude_zero DROPS any pixel where either channel " *
        "is zero, so pushing masked ch2 pixels to 0.0 would DELETE them from the correlation " *
        "instead of anti-correlating them, and the alpha ladder would read flat for a reason " *
        "unrelated to segregation. More than $(100q)% of this channel is exactly zero."))
    return b
end

# --- (2) The ch1 object mask: the frozen, already-shipped Otsu rule --------------------------

"""
    ch1_object_mask(pair) -> BitMatrix

The ch1 object mask, `_calculate_mask(build_mci(pair))[1]` -- the FROZEN, ALREADY-SHIPPED Otsu
rule of src/LoadImages.jl:235-241 (`ch .> thr` against `otsu_thresholds`), reached read-only
through spike/contract.jl. Recorded in the pre-registration as
`P13_ALPHA_MASK_RULE = :calculate_mask_ch1_once`. Reusing the shipped helper is what makes
D-16's "the segmentation/threshold rule is a pre-registered choice (D-04 applies)" true at zero
cost; do NOT invent a new thresholding rule here.

COMPUTED ONCE FROM THE UNMODIFIED ch1 AND HELD FIXED ACROSS THE WHOLE LADDER. The precedent is
src/local_map.jl:208, which rejects a tile-local Otsu because a per-tile threshold makes tiles
incomparable; a per-alpha threshold makes RUNGS incomparable in exactly the same way, and the
comparability of the rungs is the entire content of the series.

VERIFIED FINDING, recorded so nobody re-derives it: `patch_summary` never re-consumes Otsu.
`patch_summary` -> `patch` -> `correlation` read `mci.data` directly; Otsu appears only at MCI
CONSTRUCTION (src/api.jl:58, src/local_map.jl:208, src/amortized/simulator.jl:268). So
rebuilding the MCI per rung -- which does recompute ch2's threshold -- does NOT perturb the
summary. The mask is held fixed for COMPARABILITY, not for correctness of the summary.
"""
function ch1_object_mask(pair)
    return _calculate_mask(build_mci(pair))[1]
end

# --- (3) The transform ------------------------------------------------------------------------

"""
    alpha_segregate(x, y, M, alpha; b, max_value = P13_ALPHA_MAX_VALUE_BOUND) -> Matrix{Float64}

Redistribute ch2 (`y`) intensity OUT of the ch1 object mask `M` and into its complement, graded
by `alpha` in `[0, 1]`. `x` (ch1) is never modified; it is taken only to assert the frame
shapes agree. Typed on `AbstractMatrix{Float64}` so ONE code path serves the simulated and the
real substrate (D-15) -- there is no substrate branch.

    removed = alpha .* max.(y .- b, 0.0) .* M     # alpha = 0 => exactly 0.0, elementwise
    y_in    = y .- removed                        # alpha = 0 => y .- 0.0 == y BITWISE
    boost   = 1.0 + sum(removed) / sum(y[.!M])    # alpha = 0 => exactly 1.0
    y_out   = ifelse.(M, y_in, y .* boost)        # alpha = 0 => y .* 1.0 == y BITWISE

THE ELEMENTWISE FORM IS LOAD-BEARING AND MUST NOT BE RESTRUCTURED INTO A BRANCH ON A ZERO
ALPHA. The whole point of the bitwise-identity property is that it holds STRUCTURALLY, with no
special case: a special case would prove only that the special case works.

# The contract (the seven invariants, each an executable assertion in test_p13_alpha.jl)

| # | Invariant | Form |
|---|---|---|
| 1 | `alpha = 0` reproduces the input BITWISE | exact `==`, never `isapprox` |
| 2 | alpha introduces NO NEW zero | `count(iszero, y_out) == count(iszero, y)` |
| 3 | total ch2 intensity preserved | `isapprox(sum(y_out), sum(y); rtol = P13_ALPHA_INVARIANT_RTOL)` |
| 4 | the mask is alpha-invariant | `M` computed once from the UNMODIFIED ch1 |
| 5 | the mask fraction is in band | `P13_ALPHA_MASK_FRACTION_BOUNDS` |
| 6 | the realized maximum is in the declared range | `maximum(y_out) <= max_value`, RECORDED |
| 7 | `m-bar(alpha)` is monotone non-increasing | reported by `alpha_ladder_summaries` |

INVARIANT 2 IS NOT AN ALL-PIXELS-STRICTLY-POSITIVE CHECK, AND THAT NAIVE FORM IS DELIBERATELY
ABSENT FROM THIS FILE -- NOT EVEN QUOTED, so a source-grep gate cannot be defeated by prose
that names what it forbids. The all-positive property is FALSE on real data: the TIFFs already
carry exact-zero pixels BEFORE any alpha -- 2 in `positive_c2`, 3 in `negative_c2`, of
1 414 528 -- so `minimum(y_out) == 0.0` at every alpha INCLUDING alpha = 0, where `y_out` IS
`y`. It is also false on simulated substrate once the noise stage has run. Asserting it would
turn a CORRECT algorithm into a red suite. The pre-registered rule is
`P13_ALPHA_ZERO_INVARIANT = :count_iszero_preserved`, i.e. "alpha introduces no NEW zero", and
it is checked against the count on the UNMODIFIED channel. What actually drives `_exclude_zero`
is the zero SET, and pushing masked pixels down to `b > 0` leaves that set untouched.

# Guards (named errors, never a silent uninterpretable ladder)

`alpha` outside `[0, 1]`; a non-positive `b` (the `_exclude_zero` message); a mask fraction
outside `P13_ALPHA_MASK_FRACTION_BOUNDS`, because outside that band the complement cannot
absorb the redistributed mass and the boost blows up (MEASURED: a 50% mask gives boost 2.678,
while the real Otsu masks at 13.49% and 22.92% give 1.509 and 1.598); a mask covering the whole
field, leaving nowhere to redistribute to; and the post-condition `maximum(y_out) <= max_value`,
whose error message RECORDS the realized maximum. NOTHING IS EVER SILENTLY CLAMPED.

`max_value` defaults to the pre-registered `P13_ALPHA_MAX_VALUE_BOUND` (= 1.0), whose measured
basis is the real Gray TIFFs, whose values live in `[0, 1]` and whose post-boost maximum stayed
at or below 0.75. It is a KEYWORD because a declared dynamic range is a property of the
SUBSTRATE, not of this algorithm: the simulator's softplus intensities are unnormalized and
reach ~6 on a 128x128 draw, so the simulated arm declares its own range at the call site while
the realized maximum is still RECORDED by `alpha_ladder_summaries` and
`verify_alpha_invariants`. Deciding the range in here would be exactly the provenance branch
D-15 forbids. The summary is a per-patch CORRELATION and is invariant to a global affine
rescale, so the boost bites only through the mask/complement contrast WITHIN a patch -- which
is the intended mechanism, not a side effect.
"""
function alpha_segregate(x::AbstractMatrix{Float64}, y::AbstractMatrix{Float64},
                         M::AbstractMatrix{Bool}, alpha::Real;
                         b::Real, max_value::Real = P13_ALPHA_MAX_VALUE_BOUND)
    (0.0 <= alpha <= 1.0) ||
        throw(ArgumentError("alpha_segregate: alpha must be in [0, 1], got $alpha"))
    b > 0.0 || throw(ArgumentError(
        "alpha_segregate: the background floor b must be STRICTLY POSITIVE, got b = $b. " *
        "_exclude_zero (src/colocalization.jl:154-169) DROPS any pixel where either channel " *
        "is zero, so a floor of 0.0 DELETES the masked pixels from the correlation instead " *
        "of anti-correlating them and the ladder reads flat for the wrong reason."))
    (size(x) == size(y) == size(M)) || throw(ArgumentError(
        "alpha_segregate: ch1, ch2 and the mask must share a frame, got " *
        "$(size(x)) / $(size(y)) / $(size(M))"))

    # The PRE-REGISTERED mask-fraction guard. Outside this band the ladder is not merely
    # imprecise, it is uninterpretable: the complement is too small to absorb the redistributed
    # mass, so the boost blows up and the "segregation" signal is confounded with a global
    # rescale of a shrinking region.
    frac = mean(M)
    (first(P13_ALPHA_MASK_FRACTION_BOUNDS) <= frac <= last(P13_ALPHA_MASK_FRACTION_BOUNDS)) ||
        throw(ArgumentError(
            "alpha_segregate: mask fraction $frac is outside the pre-registered band " *
            "$(P13_ALPHA_MASK_FRACTION_BOUNDS). Outside it the mask complement cannot absorb " *
            "the redistributed mass and the boost blows up (measured: a 50% mask gives boost " *
            "2.678, while the real Otsu masks at 13.49% and 22.92% give 1.509 and 1.598), so " *
            "the resulting ladder is uninterpretable rather than merely noisy."))

    removed = alpha .* max.(y .- b, 0.0) .* M     # alpha = 0 -> exactly 0.0, elementwise
    y_in    = y .- removed                        # alpha = 0 -> y .- 0.0 == y BITWISE
    S_out   = sum(y[.!M])
    S_out > 0.0 || throw(ArgumentError(
        "alpha_segregate: the mask complement carries no intensity (sum = $S_out), so the ch1 " *
        "mask covers the whole field and there is nowhere to redistribute the removed ch2 " *
        "mass to. A segregation ladder needs a complement to segregate INTO."))
    boost = 1.0 + sum(removed) / S_out            # alpha = 0 -> exactly 1.0
    y_out = ifelse.(M, y_in, y .* boost)          # alpha = 0 -> y .* 1.0 == y BITWISE

    realized = maximum(y_out)
    realized <= max_value || throw(ArgumentError(
        "alpha_segregate: realized maximum $realized exceeds the declared range bound " *
        "$max_value at alpha = $alpha (boost = $boost, mask fraction = $frac). The value is " *
        "REPORTED, never clamped: silently clamping would change which pixels _exclude_zero " *
        "sees and would make the rung incomparable with the others. Either the substrate's " *
        "declared dynamic range is wrong for this input or the boost is too large."))
    return Matrix{Float64}(y_out)
end

# --- (4) The ladder ---------------------------------------------------------------------------

"""
    alpha_ladder(pair; ladder = P13_ALPHA_LADDER, mask = nothing, floor = nothing,
                 max_value = P13_ALPHA_MAX_VALUE_BOUND)
        -> (pairs, mask, floor, ladder)

Build the whole alpha ladder from one `pair = [ch1, ch2]`, returning one two-channel
`Vector{Matrix{Float64}}` per rung in `ladder` order, alongside the mask and floor actually
used. ch1 is passed through UNTOUCHED in every rung.

THE MASK AND THE FLOOR ARE LADDER-LEVEL CONSTANTS, computed ONCE from the UNMODIFIED pair when
not supplied. Recomputing either per rung is the comparability bug this signature exists to
prevent: a per-alpha Otsu threshold or a per-alpha quantile would make each rung a different
experiment, and a monotone curve across rungs would then mean nothing (the src/local_map.jl:208
"tile-local Otsu makes tiles incomparable" precedent, applied to alpha instead of tiles). Pass
`mask`/`floor` explicitly only to SHARE ladder-level constants across arms, never to vary them.

`max_value` is forwarded to `alpha_segregate` and is the substrate's DECLARED dynamic range;
see that docstring for why the range is a caller-supplied knob rather than a branch in here.
"""
function alpha_ladder(pair; ladder = P13_ALPHA_LADDER, mask = nothing, floor = nothing,
                      max_value::Real = P13_ALPHA_MAX_VALUE_BOUND)
    length(pair) == 2 || throw(ArgumentError(
        "alpha_ladder: expected a 2-channel pair, got $(length(pair)) channels"))
    x = Matrix{Float64}(pair[1])
    y = Matrix{Float64}(pair[2])
    M = mask === nothing ? ch1_object_mask([x, y]) : mask
    b = floor === nothing ? alpha_background_floor(y) : floor
    pairs = [Vector{Matrix{Float64}}([x, alpha_segregate(x, y, M, a; b = b,
                                                         max_value = max_value)])
             for a in ladder]
    return (pairs = pairs, mask = M, floor = b, ladder = ladder)
end

# --- (5) The per-rung summary readout ---------------------------------------------------------

"""
    alpha_ladder_summaries(pair; ladder = P13_ALPHA_LADDER,
                           max_value = P13_ALPHA_MAX_VALUE_BOUND)
        -> Vector of (alpha, mbar, n_missing, max_value, sum_ratio, n_zero)

The per-rung readout of the ladder, one NamedTuple per alpha.

`mbar` is the MASK-WEIGHTED mean of the PRESENT continuous summary rows: rows 1 to G^2 of
`encode_d01(patch_summary(build_mci(rung)))` weighted by the parallel present/absent rows
G^2+1 to 2G^2 (spike/data/encode.jl:56-75). The weighting is not cosmetic. `encode_d01` imputes
a `missing` patch as 0.0, so an UNWEIGHTED mean would silently pull `m-bar` toward 0 exactly in
the regime where patches start dropping out -- i.e. it would manufacture the appearance of a
decaying ladder out of missingness. Returns `NaN` when every patch is missing, mirroring
`induced_mu`'s NaN-not-throw degeneracy idiom (spike/contract.jl:104).

`n_missing` COUNTS THE ABSENT PATCHES AND MUST NOT RISE WITH ALPHA. A rising count is the
Pitfall-2 warning sign that pixels are being DELETED rather than moved, i.e. that the
strictly-positive floor has stopped protecting the zero set. It is reported per rung rather
than only aggregated so a partial breakdown is visible at the alpha where it starts.

`max_value`, `sum_ratio` and `n_zero` are the realized maximum, the
`sum(y_alpha) / sum(y_0)` intensity-conservation residual and the exact-zero count of the
transformed channel. They exist so the real arm (plan 13-16) and the simulated arm (plan 13-13)
can REPORT the measured maximum, the conservation residual and the zero count per rung without
re-deriving any of them -- the pre-registration requires the maximum recorded, not clamped.
"""
function alpha_ladder_summaries(pair; ladder = P13_ALPHA_LADDER,
                                max_value::Real = P13_ALPHA_MAX_VALUE_BOUND)
    L = alpha_ladder(pair; ladder = ladder, max_value = max_value)
    y0_sum = sum(L.pairs[1][2])
    out = NamedTuple[]
    for (a, rung) in zip(L.ladder, L.pairs)
        enc  = encode_d01(patch_summary(build_mci(rung)))
        n    = length(enc) ÷ 2                     # G^2; vals in 1:n, present/absent in n+1:2n
        vals = @view enc[1:n]
        msk  = @view enc[(n + 1):(2n)]
        wsum = sum(msk)
        mbar = wsum == 0.0 ? NaN : sum(vals .* msk) / wsum
        push!(out, (alpha       = a,
                    mbar        = mbar,
                    n_missing   = Int(n - wsum),
                    max_value   = maximum(rung[2]),
                    sum_ratio   = sum(rung[2]) / y0_sum,
                    n_zero      = count(iszero, rung[2])))
    end
    return out
end

# --- (6) The invariant verifier ---------------------------------------------------------------

"""
    verify_alpha_invariants(pair; ladder = P13_ALPHA_LADDER,
                            rtol = P13_ALPHA_INVARIANT_RTOL,
                            max_value = P13_ALPHA_MAX_VALUE_BOUND)
        -> NamedTuple of booleans and measured residuals

Report -- never enforce -- every D-16 invariant on one pair, in the shape of
`verify_summary_invariance` (spike/validation/ood.jl:567-581) but with assertions STRONGER than
statistical, because at alpha = 0 the output must be the input bit for bit.

Fields: `bitwise_alpha0`, `no_new_zeros`, `intensity_conserved`, `mask_invariant`,
`mask_fraction_ok`, `max_value_ok`, `mbar_monotone`, `missing_nonincreasing`,
`max_rel_intensity_err`, `source_zero_count`, `realized_mask_fraction`, `realized_max_value`.

  * THE FIELD IS `no_new_zeros`, AND NEVER AN ALL-POSITIVE FIELD -- whose name is not even
    quoted here, so the source-grep gate cannot be defeated by prose. It compares
    `count(iszero, y_out)`
    against the count on the UNMODIFIED ch2, which is what makes the check valid on BOTH
    substrates: the committed real TIFFs carry 2-3 source zeros before any alpha, so an
    all-positive field would report a correct transform as broken. `source_zero_count` is
    returned so a reader sees the reference the comparison was made against rather than having
    to trust it.
  * `bitwise_alpha0` uses exact `==`, never `isapprox`, mirroring the Phase-11 zero-epsilon
    regression pattern: an approximate check at alpha = 0 would pass a transform that quietly
    perturbs every pixel by a rounding step, and a perturbed zero set is a changed summary.
  * `mask_invariant` RECOMPUTES `ch1_object_mask` on each transformed pair and compares it to
    the ladder mask. ch1 is untouched by construction, so this is a real leak check, not a
    tautology -- it fails the moment a future edit writes to channel 1.
  * `mbar_monotone` is NON-STRICT non-increase at tolerance `<= previous + 1e-9`, a float-noise
    allowance and not a loosened bar; a genuinely rising rung exceeds it by orders of magnitude.
  * `mask_fraction_ok` and `max_value_ok` are REPORTED here rather than thrown, which is why the
    ladder is built with the guards' reporting path (an out-of-band mask short-circuits to a
    report with the realized fraction, and the range check is evaluated against `max_value`
    after the fact). Enforcement lives in `alpha_segregate`; this function's job is to say what
    was measured.
"""
function verify_alpha_invariants(pair; ladder = P13_ALPHA_LADDER,
                                 rtol::Real = P13_ALPHA_INVARIANT_RTOL,
                                 max_value::Real = P13_ALPHA_MAX_VALUE_BOUND)
    x = Matrix{Float64}(pair[1])
    y = Matrix{Float64}(pair[2])
    M = ch1_object_mask([x, y])
    b = alpha_background_floor(y)
    frac = mean(M)
    src_zeros = count(iszero, y)

    # An out-of-band mask fraction is REPORTED, not thrown: a verifier that throws cannot tell a
    # caller which invariant failed. The remaining fields are false because no interpretable
    # ladder exists to measure them on.
    in_band = first(P13_ALPHA_MASK_FRACTION_BOUNDS) <= frac <= last(P13_ALPHA_MASK_FRACTION_BOUNDS)
    if !in_band
        return (bitwise_alpha0 = false, no_new_zeros = false, intensity_conserved = false,
                mask_invariant = false, mask_fraction_ok = false, max_value_ok = false,
                mbar_monotone = false, missing_nonincreasing = false,
                max_rel_intensity_err = NaN, source_zero_count = src_zeros,
                realized_mask_fraction = frac, realized_max_value = NaN)
    end

    # Built with the range check DISABLED so the report can state the realized maximum and
    # compare it against `max_value` itself, instead of dying inside the transform.
    L = alpha_ladder([x, y]; ladder = ladder, mask = M, floor = b, max_value = Inf)
    S = alpha_ladder_summaries([x, y]; ladder = ladder, max_value = Inf)

    y0 = L.pairs[1][2]
    sum_y = sum(y)
    rel_errs = [abs(sum(p[2]) - sum_y) / abs(sum_y) for p in L.pairs]

    bitwise   = first(ladder) == 0.0 && y0 == y            # exact ==, never isapprox
    no_new    = all(count(iszero, p[2]) == src_zeros for p in L.pairs)
    conserved = all(isapprox(sum(p[2]), sum_y; rtol = rtol) for p in L.pairs)
    ch1_fixed = all(p[1] == x for p in L.pairs)
    mask_inv  = ch1_fixed && all(ch1_object_mask(p) == M for p in L.pairs)

    realized_max = maximum(s.max_value for s in S)
    max_ok = realized_max <= max_value

    mbars = [s.mbar for s in S]
    monotone = all(mbars[i] <= mbars[i - 1] + 1e-9 for i in 2:length(mbars))
    miss = [s.n_missing for s in S]
    miss_ok = all(miss[i] <= miss[1] for i in 2:length(miss))

    return (bitwise_alpha0 = bitwise, no_new_zeros = no_new, intensity_conserved = conserved,
            mask_invariant = mask_inv, mask_fraction_ok = in_band, max_value_ok = max_ok,
            mbar_monotone = monotone, missing_nonincreasing = miss_ok,
            max_rel_intensity_err = maximum(rel_errs), source_zero_count = src_zeros,
            realized_mask_fraction = frac, realized_max_value = realized_max)
end
