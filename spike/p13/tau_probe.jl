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

# spike/p13/tau_probe.jl --- D-06 summary-resolution probe (fit-free, UNPAIRED).
#
# (a) WHAT TAU IS, AND WHY IT IS MEASURED RATHER THAN CHOSEN.
# Phase 13 promotes "mutually exclusive" to a first-class hypothesis, and the whole
# three-way cut hangs off ONE number. D-05 says a sample is colocalized when
# rho_s > +tau and segregated when rho_s < -tau; everything between is RANDOM. So tau
# is not a convenience threshold, it is the DEFINITION of the random class -- and if it
# were chosen, "random" would mean whatever made the confusion matrix look best.
# D-06 fixes that by defining random operationally: INDISTINGUISHABLE FROM ZERO AT THIS
# SUMMARY'S RESOLUTION. This file measures that resolution.
#
#   A(delta) = AUC( mbar | rho_true = 0   vs   mbar | rho_true = -+delta )
#   tau      = the SMALLEST delta on the frozen P13_TAU_DELTA_GRID with A(delta) >= P13_TAU_AUC
#
# Read that as a claim about the SUMMARY, not about any network: it says the fixed 8x8
# patch-correlation statistic cannot tell rho = 0 from rho = tau/2 more reliably than a
# coin, so calling anything inside +-tau "colocalized" or "segregated" would be
# publishing a distinction the measurement instrument does not support.
#
# (b) WHY THE DESIGN IS UNPAIRED -- THE SINGLE MOST IMPORTANT DIFFERENCE FROM PHASE 11.
# Phase 11's D-06 probe is DELIBERATELY PAIRED: same latent field, only the geometry
# changes, because it measures a DISPLACEMENT and a displacement is only defined against
# the un-displaced counterfactual. This probe measures something else entirely --
# DISCRIMINABILITY FOR A SINGLE UNSEEN IMAGE, where the latent field, the spillover, the
# label efficiency, the registration error, the chromatic term, the noise level and the
# frame size are ALL unknown. Pairing here would hold every one of those fixed across the
# two arms and report the resolution of a counterfactual the deployed tool never has;
# it would UNDERSTATE tau by a large factor, the random band would come out razor-thin,
# the random class would be starved, and the confusion matrix would then look excellent
# for a reason that has nothing to do with the net (13-RESEARCH Pitfall 8). The two arms
# below therefore draw from DISJOINT PHILOX KEYS: their nuisances and their latent fields
# differ, which is exactly what makes the measured number describe one unseen image.
# `P13_TAU_DESIGN = :unpaired` records this in the pre-registration; it is asserted here.
#
# (c) WHY THE STATISTIC IS mbar AND NOT SOMETHING LEARNED.
# mbar is the MASK-WEIGHTED mean of the PRESENT continuous summary rows -- the canonical
# scalar readout of the summary with respect to rho, and literally the quantity the frozen
# ghat calibration was built on (spike/simulator/ghat.jl:34 records the sweep as
# E[mu | rho_true], where mu IS the induced mean per-patch correlation). It is FIT-FREE,
# which is what makes it honestly pre-registerable: there is no fitting surface to snoop.
# A LEARNED readout -- an LDA on the 64 rows, or the NPE's own rho-hat -- would be MORE
# powerful and would therefore report a SMALLER tau. That is precisely why it is refused:
# it would make tau a property of a trained MODEL rather than of the SUMMARY, which
# contradicts D-06's wording, and it would re-introduce the Phase-11 net dependency this
# probe exists to avoid. A model-dependent number may be reported as a SECONDARY column.
# `P13_TAU_STATISTIC = :mbar_mask_weighted` is the frozen record of the choice.
#
# (d) THE PROBE STOPS AT encode_d01(patch_summary(...)).
# It never applies the frozen summary standardization, so it needs NO `zt` and NO trained
# network -- the sequencing win of 13-RESEARCH E4. There is deliberately no `zt` argument
# anywhere in this file: tau is a property of the RAW per-patch correlations, which is why
# it can be measured off the simulator alone the moment the Phase-11 simulator merges,
# while the labelled data generation stays behind the full Phase-11 net precondition.
#
# (e) DECOUPLING (hard constraint, CLAUDE.md): spike-local. `src/` is reached ONLY
# read-only and transitively, through spike/contract.jl's include chain, for the UNCHANGED
# patch/correlation summary. No src/ byte, no manifest byte, no new dependency, and no
# `corpus/` path -- the sealed holdout is never opened by this phase.
#
# WHAT THIS FILE IS *NOT*. It is the probe MACHINERY only. The REPORTED run and the Tier-2
# freeze of P13_TAU are plan 13-10's `spike/p13/run_tau_probe.jl`, which must first confirm
# `simulator_provenance_guard().post_p11`. Including this file simulates nothing.

using Statistics             # median (the bootstrap column), mean is not used on the mask path
import Random123: Philox4x   # counter-based per-index keys: the two disjoint unpaired arms

# ORDER MATTERS, and every include is guarded for idempotency under a runner (house idiom).
# The pre-registration first (the frozen grid, bar, R, statistic, design and the F5 mixture),
# then the prior (sample_prior + the frozen ghat knots), the forward simulator, the read-only
# src/ summary contract, the 128-row encoding the statistic is read out of, and the seeding
# primitives whose salt-XOR-into-the-first-key-word idiom the two arms copy.
isdefined(@__MODULE__, :P13_DEV_SEED) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :sample_prior) || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair) || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)    || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)   || include(joinpath(@__DIR__, "..", "data", "encode.jl"))
isdefined(@__MODULE__, :sample_rng)   || include(joinpath(@__DIR__, "..", "data", "seeding.jl"))
# The include below exists SOLELY to REUSE the in-repo hand-rolled tie-aware AUC, never to
# re-implement it and never to add an AUC package. spike/test/runtests.jl asserts the two
# obvious AUC/ML packages (ROCAnalysis, MLJ) are ABSENT from the pinned environment, so a
# fresh dependency here would turn the suite red on purpose. The guard is keyed on the
# function itself so the heavier validation chain is pulled at most once.
isdefined(@__MODULE__, :roc_auc)      || include(joinpath(@__DIR__, "..", "validation", "ood.jl"))

# =====================================================================================
# 1. The statistic
# =====================================================================================

"""
    mbar_from_summary(Z128::AbstractVector) -> Float64

The MASK-WEIGHTED mean of the PRESENT continuous summary rows -- the D-06 statistic,
frozen as `P13_TAU_STATISTIC = :mbar_mask_weighted`.

With `G2 = length(Z128) ÷ 2`, rows `1:G2` are the per-patch correlations and rows
`G2+1:2G2` the parallel present/absent indicator (spike/data/encode.jl:56-75), so

    mbar = sum(vals .* mask) / sum(mask)

Returns `NaN` when every patch is absent, mirroring `induced_mu`'s NaN-not-throw
degeneracy idiom (spike/contract.jl:104) rather than reporting a silent zero.

# Why the weighting is load-bearing, not cosmetic

`encode_d01` imputes a `missing` patch as `0.0`. An UNWEIGHTED `mean(vals)` would
therefore let the missing-to-zero imputation pull the statistic toward 0, and ONE
degenerate all-absent grid can move the value by order 1. On a resolution probe that is
not noise, it is BIAS WITH A DIRECTION: it drags both arms toward the same point and so
DEFLATES the measured discriminability, which inflates tau, which widens the random band.
The mask rows are exactly the rows `_summary_row_partition` (src/amortized/summary.jl:89-98)
itself excludes from the standardization path -- so weighting by them is the house reading
of the encoding, not a new convention. `spike/p13/alpha_series.jl` reads its ladder the
same way, and the two must not diverge.
"""
function mbar_from_summary(Z128::AbstractVector)
    G2 = length(Z128) ÷ 2
    G2 > 0 || throw(ArgumentError(
        "mbar_from_summary: expected a 2*G^2-row encoding, got length $(length(Z128))"))
    iseven(length(Z128)) || throw(ArgumentError(
        "mbar_from_summary: the encoding must have an EVEN row count (G^2 value rows plus " *
        "G^2 present/absent rows), got length $(length(Z128))"))
    vals = view(Z128, 1:G2)
    mask = view(Z128, (G2 + 1):(2 * G2))
    w = sum(mask)
    w == 0.0 && return NaN     # every patch absent: report the degeneracy, never a silent 0
    return sum(vals .* mask) / w
end

# =====================================================================================
# 2. The image-size mixture (F5)
# =====================================================================================

"""
    sample_imsize(rng::AbstractRNG, set::Tuple, weights::Tuple) -> Tuple{Int,Int}

Draw one image size from an EXPLICIT categorical `(set, weights)` by hand-rolled
inverse CDF, consuming exactly ONE `rand(rng)` so the surrounding draw stream advances
identically whatever mixture is passed. No extra dependency.

# Why this is a THREE-ARGUMENT method and not a keyword one

`spike/data/seeding.jl:106` already owns `sample_imsize(rng::AbstractRNG)` over the
generation-time `IMSIZE_SET`, and that file is in scope here (the validation chain pulls
it). Keyword arguments do NOT participate in Julia dispatch, so a
`sample_imsize(rng; set, weights)` method would have the SAME positional signature and
would silently OVERWRITE the generation sampler -- or be overwritten by it, depending on
include order, in which case this probe would quietly measure under the 256-squared-inclusive
generation mixture instead of the frozen F5 one. Taking the mixture POSITIONALLY makes this
an ADDITIVE method: both samplers coexist and neither can clobber the other.

# Why the mixture and not 256 squared

`P13_IMSIZE_SET` / `P13_IMSIZE_WEIGHTS` are READ from the frozen gate file through
`spike/p13/consts.jl`, never retyped. F5 covariate shift is a NAMED LIMIT of this project:
rank, coverage and calibration claims hold only under the joint the estimator was trained
on. A 256-squared-only probe would report the resolution of a regime the deployed tool never
sees, so the measured tau would not describe the deployment distribution it is used to
partition.
"""
function sample_imsize(rng::AbstractRNG, set::Tuple, weights::Tuple)
    length(set) == length(weights) || throw(ArgumentError(
        "sample_imsize: set and weights must agree in length, got " *
        "$(length(set)) / $(length(weights))"))
    r = rand(rng)
    c = 0.0
    @inbounds for i in 1:length(set)
        c += weights[i]
        r <= c && return set[i]
    end
    return set[end]     # float-rounding guard: r just over the weight sum -> the last bin
end

# =====================================================================================
# 3. The two UNPAIRED arms
# =====================================================================================

"""
    tau_probe_arm_keys(; seed = P13_DEV_SEED) -> (reference, contrast)

The two Philox FIRST-KEY-WORDS the probe's arms ride, as `UInt64`s:

    reference = seed xor P13_SALT                  # the reported Phase-13 key (p13_rng's key)
    contrast  = seed xor P13_SALT xor HOLDOUT_SALT  # a PROVABLY DISJOINT second key

This is the `spike/data/seeding.jl:75-96` idiom verbatim -- XOR a fixed nonzero salt into
the first key word to carve a second stream that can never collide with the first for any
counter. `HOLDOUT_SALT` is reused rather than a fresh salt invented because carving a
disjoint arm is EXACTLY what that constant already means in this repository, and inventing
a new salt in a non-pre-registration file would put a stream-defining value outside the
frozen record.

**The disjointness of these two keys is what makes the design UNPAIRED** (`P13_TAU_DESIGN`).
Because the arms differ in their KEY, every nuisance and the whole latent field differ
between them, so `A(delta)` answers "can the summary tell a rho = 0 image from a
rho = -+delta image when nothing else is held fixed" -- the question the deployed tool
actually faces. A paired design would key both arms identically and report a
counterfactual resolution (13-RESEARCH Pitfall 8).

Pass `seed = P13_FIX_SEED` for fixture-scale work so a quick gate never consumes -- and so
never pre-observes -- the reported stream (D-01).
"""
function tau_probe_arm_keys(; seed = P13_DEV_SEED)
    reference = UInt64(seed) ⊻ P13_SALT
    contrast  = reference ⊻ HOLDOUT_SALT
    # A defensive equality check rather than a comment: if a future edit made the two arms
    # share a key the probe would silently become PAIRED and would understate tau, and
    # nothing about its output shape would change.
    reference == contrast && error(
        "tau_probe_arm_keys: the two arms collapsed onto ONE key, which would make the " *
        "probe PAIRED (P13_TAU_DESIGN = :unpaired) and would understate tau by a large " *
        "factor. Refusing to return a paired arm pair.")
    return (reference = reference, contrast = contrast)
end

"""
    summary_draws(rho, R, stream_key; counter_base = 0,
                  imsize_set = P13_IMSIZE_SET, imsize_weights = P13_IMSIZE_WEIGHTS)
        -> Vector{Float64}

`R` independent `mbar` values at a PINNED `rho_true = rho`, with every nuisance drawn from
its OWN full prior and the frame size drawn from the frozen F5 mixture. One arm of the
probe. Deterministic under `(stream_key, counter_base)`, and bitwise-reproducible across
processes and thread counts because the RNG is keyed PER INDEX:

    rng   = Philox4x(UInt64, (UInt64(stream_key), UInt64(counter_base + r)))
    theta = merge(sample_prior(rng), (; rho_true = rho))     # nuisances ~ pi, rho PINNED
    imsize -> sample_imsize(rng, imsize_set, imsize_weights)
    mbar_from_summary(encode_d01(patch_summary(build_mci(simulate_pair(rng, theta; imsize)))))

# Two details that are easy to get wrong and silent when wrong

**The merge must OVERWRITE, not append.** `sample_prior` names the colocalization knob with
the Unicode field the simulator uses, and `merge` on a NamedTuple replaces an existing key
IN PLACE while preserving field ORDER -- which matters because the house theta-to-vector map
is positional (`collect(values(theta))`; `src/amortized/ood.jl` `_theta_tuple` reads
`v[1]..v[7]`). Merging under a DIFFERENT spelling would APPEND a ninth field and leave the
knob at its prior draw, so both arms would be distributed identically, `A(delta)` would come
back at chance for every delta, and the probe would report `:below_resolution` for a reason
that has nothing to do with the summary's resolution. The arity and the pinned value are
therefore ASSERTED below rather than trusted.

**The RNG is threaded per index, never seeded inside the loop.** `Random.seed!(r)` in a loop
gives correlated low-entropy streams and is forbidden by `P13_GLOBAL_RNG_DISCIPLINE`; the
counter-based key is the whole reason the arms are reproducible and provably disjoint.

Entries can be `NaN` when a draw yields an all-absent patch grid. They are RETURNED as `NaN`
rather than dropped here, so the caller can COUNT the degeneracies instead of silently
shrinking an arm; `tau_auc_at` filters and reports the count.
"""
function summary_draws(rho::Real, R::Integer, stream_key::Integer;
                       counter_base::Integer = 0,
                       imsize_set::Tuple = P13_IMSIZE_SET,
                       imsize_weights::Tuple = P13_IMSIZE_WEIGHTS)
    R >= 1 || throw(ArgumentError("summary_draws: R must be >= 1, got $R"))
    (-1.0 <= rho <= 1.0) ||
        throw(ArgumentError("summary_draws: rho must be in [-1, 1], got $rho"))
    out = Vector{Float64}(undef, R)
    for r in 1:R
        rng = Philox4x(UInt64, (UInt64(stream_key), UInt64(counter_base + r)))
        theta = merge(sample_prior(rng), (; ρ_true = rho))
        # The overwrite-not-append assertions (see the docstring): a merge under the wrong
        # spelling would leave the knob at its prior draw and flatten the whole curve.
        length(theta) == 8 || error(
            "summary_draws: pinning rho changed the theta arity to $(length(theta)); the " *
            "merge APPENDED a field instead of OVERWRITING the colocalization knob, so " *
            "the arm would be distributed identically to the reference arm.")
        theta.ρ_true == rho || error(
            "summary_draws: rho was not pinned (got $(theta.ρ_true), wanted $rho).")
        imsize = sample_imsize(rng, imsize_set, imsize_weights)
        out[r] = mbar_from_summary(
            encode_d01(patch_summary(build_mci(simulate_pair(rng, theta; imsize = imsize)))))
    end
    return out
end

# =====================================================================================
# 4. A(delta) -- the resolution measure at one delta
# =====================================================================================

"""
    _auc_magnitude(reference, contrast) -> Float64

`max(A, 1 - A)` on the third element of the REUSED in-repo tie-aware `roc_auc`
(spike/validation/ood.jl:309-335, the exact Mann-Whitney U form). Returns `NaN` when either
arm is empty.

The magnitude, not the raw AUC, is the reported quantity. `mbar` DECREASES with rho on the
negative side, so the raw AUC of the `rho = -delta` arm against the `rho = 0` arm runs
toward 0, not toward 1: MEASURED at fixture scale, `A_neg = 0.0` at `delta = 0.5`, which is
PERFECT separation reading as "worse than chance". A bar applied to the raw value would
score a perfectly resolved delta as unresolved, and the direction of the readout is not part
of what D-06 asks. Discriminability is symmetric under a sign flip of the score; the bar is
applied to `max(A, 1 - A)`.
"""
function _auc_magnitude(reference::AbstractVector, contrast::AbstractVector)
    (isempty(reference) || isempty(contrast)) && return NaN
    a = roc_auc(reference, contrast)[3]     # REUSED, never re-implemented
    return isnan(a) ? NaN : max(a, 1.0 - a)
end

_finite_only(v::AbstractVector) = collect(Float64, Iterators.filter(isfinite, v))

"""
    tau_auc_at(delta; R = P13_TAU_R, arm_keys = tau_probe_arm_keys(),
               both_directions = P13_TAU_BOTH_DIRECTIONS,
               bootstrap = P13_TAU_BOOTSTRAP_B, boot_rng = p13_rng(P13_TAU_COUNTER),
               reference = nothing,
               imsize_set = P13_IMSIZE_SET, imsize_weights = P13_IMSIZE_WEIGHTS)
        -> NamedTuple

`A(delta)` in BOTH directions plus the reported magnitude:

    (neg, pos, auc = max(neg, pos), auc_median_boot, n_reference, n_neg, n_pos, n_degenerate)

`neg` separates `mbar | rho = 0` from `mbar | rho = -delta`, `pos` the same at `+delta`, and
`auc` is the MAX. Taking the max is the CONSERVATIVE reading and it is frozen as
`P13_TAU_BOTH_DIRECTIONS = true` for a measured reason: the prior atoms are ASYMMETRIC by a
factor of about 2.5 (0.048527 of the prior mass on the `rho = -0.99` rail against 0.019108 on
`+0.99`, exact truncated-Cauchy CDF at the ghat clamps), so the resolution may be asymmetric
too, while D-05 assumes ONE symmetric tau. The larger of the two is the honest single number.

**The two arms are UNPAIRED**: the reference arm rides `arm_keys.reference` and both contrast
arms ride `arm_keys.contrast`, two PROVABLY DISJOINT Philox keys, so the nuisances and the latent
field differ across a comparison. That disjointness IS the unpaired design
(`P13_TAU_DESIGN`); see `tau_probe_arm_keys`. Within `arm_keys.contrast` the two directions take
NON-OVERLAPPING counter blocks (`1:R` and `R+1:2R`) so `-delta` and `+delta` never share a
draw either.

`reference` accepts a PRE-COMPUTED reference arm so `tau_curve` can draw `rho = 0` ONCE and
read the whole grid against it. That is deliberate common-random-numbers across grid points:
it costs one arm instead of one per delta, and it makes the curve's shape a statement about
delta rather than about reference-arm noise. Within any single comparison the arms are still
unpaired, which is the property that matters.

Non-finite draws (an all-absent patch grid) are filtered and COUNTED in `n_degenerate`
rather than silently treated as scores -- `NaN` compares false against everything, so an
unfiltered `NaN` would quietly contribute a zero to the U statistic.

`bootstrap` resamples both arms `bootstrap` times and returns the MEDIAN magnitude, because
one degenerate grid can move `mbar` by order 1 through the missing-to-zero imputation path
and the mean is not robust to that. `boot_rng` makes the resampling reproducible.

Throws an `ArgumentError` on `delta <= 0`: `A(0)` is chance by construction and would report
"the summary cannot resolve zero from zero" as a finding.
"""
function tau_auc_at(delta::Real;
                    R::Integer = P13_TAU_R,
                    arm_keys = tau_probe_arm_keys(),
                    both_directions::Bool = P13_TAU_BOTH_DIRECTIONS,
                    bootstrap::Integer = P13_TAU_BOOTSTRAP_B,
                    boot_rng = p13_rng(P13_TAU_COUNTER),
                    reference = nothing,
                    imsize_set::Tuple = P13_IMSIZE_SET,
                    imsize_weights::Tuple = P13_IMSIZE_WEIGHTS)
    delta > 0.0 || throw(ArgumentError(
        "tau_auc_at: delta must be STRICTLY POSITIVE, got $delta. A(0) is chance by " *
        "construction, so a zero delta would report 'the summary cannot resolve zero from " *
        "zero' as though it were a measurement."))

    ref_raw = reference === nothing ?
        summary_draws(0.0, R, arm_keys.reference; counter_base = 0,
                      imsize_set = imsize_set, imsize_weights = imsize_weights) :
        collect(Float64, reference)
    neg_raw = summary_draws(-delta, R, arm_keys.contrast; counter_base = 0,
                           imsize_set = imsize_set, imsize_weights = imsize_weights)
    pos_raw = both_directions ?
        summary_draws(delta, R, arm_keys.contrast; counter_base = R,
                      imsize_set = imsize_set, imsize_weights = imsize_weights) :
        Float64[]

    ref = _finite_only(ref_raw)
    neg = _finite_only(neg_raw)
    pos = _finite_only(pos_raw)
    n_degenerate = (length(ref_raw) - length(ref)) + (length(neg_raw) - length(neg)) +
                   (length(pos_raw) - length(pos))

    a_neg = _auc_magnitude(ref, neg)
    a_pos = both_directions ? _auc_magnitude(ref, pos) : NaN
    auc = both_directions ? max(a_neg, a_pos) : a_neg

    return (neg = a_neg, pos = a_pos, auc = auc,
            auc_median_boot = _boot_auc_median(ref, neg, pos, bootstrap, boot_rng;
                                               both_directions = both_directions),
            n_reference = length(ref), n_neg = length(neg), n_pos = length(pos),
            n_degenerate = n_degenerate)
end

"""
    _boot_auc_median(ref, neg, pos, B, rng; both_directions) -> Float64

The median of the reported magnitude over `B` bootstrap resamples of the arms (each arm
resampled with replacement to its own length). Returns `NaN` for `B <= 0` or an empty arm.

The MEDIAN is the point of the exercise: `encode_d01` imputes an absent patch as `0.0`, so a
single degenerate grid can move `mbar` by order 1, and a mean AUC over resamples inherits
that leverage while a median does not. `P13_TAU_BOOTSTRAP_B` freezes `B`, and the mean is
still reported alongside by the runner (plan 13-10) so the two can be compared rather than
one being quietly preferred.
"""
function _boot_auc_median(ref::AbstractVector, neg::AbstractVector, pos::AbstractVector,
                          B::Integer, rng; both_directions::Bool = true)
    B <= 0 && return NaN
    (isempty(ref) || isempty(neg)) && return NaN
    nr = length(ref); nn = length(neg); np = length(pos)
    out = Vector{Float64}(undef, B)
    for b in 1:B
        r = [ref[rand(rng, 1:nr)] for _ in 1:nr]
        an = _auc_magnitude(r, [neg[rand(rng, 1:nn)] for _ in 1:nn])
        if both_directions && np > 0
            ap = _auc_magnitude(r, [pos[rand(rng, 1:np)] for _ in 1:np])
            out[b] = max(an, ap)
        else
            out[b] = an
        end
    end
    return median(out)
end

# =====================================================================================
# 5. The whole curve
# =====================================================================================

"""
    tau_curve(; grid = P13_TAU_DELTA_GRID, R = P13_TAU_R,
              bootstrap = P13_TAU_BOOTSTRAP_B, arm_keys = tau_probe_arm_keys(),
              both_directions = P13_TAU_BOTH_DIRECTIONS,
              boot_rng = p13_rng(P13_TAU_COUNTER),
              imsize_set = P13_IMSIZE_SET, imsize_weights = P13_IMSIZE_WEIGHTS)
        -> Vector{NamedTuple}

`A(delta)` over the WHOLE frozen grid, one NamedTuple per delta:

    (delta, auc_neg, auc_pos, auc, auc_median_boot, n_reference, n_neg, n_pos, n_degenerate)

**The whole curve is always computed and always reported, never short-circuited at the first
delta that clears the bar.** Two reasons, both about the reader rather than the runner:
a reader can re-read tau at ANY bar from the reported curve without re-running the probe (so
the pre-registered `P13_TAU_AUC = 0.90` stops being a black box), and a curve that is FLAT or
NON-MONOTONE is diagnostic in itself -- it says the statistic is not tracking rho at all,
which is a different and more serious finding than "tau is large".

The `rho = 0` reference arm is drawn ONCE and shared across the grid (common random numbers);
see `tau_auc_at` for why that is deliberate and why it does not make any single comparison
paired.

Defaults ARE the frozen knobs: `P13_TAU_DELTA_GRID`, `P13_TAU_R`, `P13_TAU_BOOTSTRAP_B` and
`P13_TAU_BOTH_DIRECTIONS`. The keywords exist so a FIXTURE-scale suite can run the same code
path off `P13_FIX_SEED` at a small `R`; they are not there to be retuned for the reported run.
"""
function tau_curve(; grid = P13_TAU_DELTA_GRID,
                   R::Integer = P13_TAU_R,
                   bootstrap::Integer = P13_TAU_BOOTSTRAP_B,
                   arm_keys = tau_probe_arm_keys(),
                   both_directions::Bool = P13_TAU_BOTH_DIRECTIONS,
                   boot_rng = p13_rng(P13_TAU_COUNTER),
                   imsize_set::Tuple = P13_IMSIZE_SET,
                   imsize_weights::Tuple = P13_IMSIZE_WEIGHTS)
    isempty(grid) && throw(ArgumentError(
        "tau_curve: the delta grid is empty. The grid is pre-registered " *
        "(P13_TAU_DELTA_GRID) and is never constructed at run time."))
    ref = summary_draws(0.0, R, arm_keys.reference; counter_base = 0,
                        imsize_set = imsize_set, imsize_weights = imsize_weights)
    out = NamedTuple[]
    for d in grid
        a = tau_auc_at(d; R = R, arm_keys = arm_keys, both_directions = both_directions,
                       bootstrap = bootstrap, boot_rng = boot_rng, reference = ref,
                       imsize_set = imsize_set, imsize_weights = imsize_weights)
        push!(out, (delta = d, auc_neg = a.neg, auc_pos = a.pos, auc = a.auc,
                    auc_median_boot = a.auc_median_boot,
                    n_reference = a.n_reference, n_neg = a.n_neg, n_pos = a.n_pos,
                    n_degenerate = a.n_degenerate))
    end
    return out
end

# =====================================================================================
# 6. Reading tau off the curve -- and the abort criterion
# =====================================================================================

"""
    ghat_knot_spacing_near_zero() -> Float64

The spacing of the frozen `GHAT_RHO_KNOTS` adjacent to zero (0.0825 on the shipped knots).
Derived from the knots, never retyped, so it cannot drift out of agreement with the
calibration it describes.
"""
function ghat_knot_spacing_near_zero()
    k = GHAT_RHO_KNOTS
    i = argmin(abs.(k))
    lo = i > 1 ? abs(k[i] - k[i - 1]) : Inf
    hi = i < length(k) ? abs(k[i + 1] - k[i]) : Inf
    return min(lo, hi)
end

"""
    tau_from_curve(curve; bar = P13_TAU_AUC, warn_sub_knot = true) -> NamedTuple

Read tau off a reported curve: `(tau, status, delta_index, bar, sub_knot)`.

    tau = the FIRST delta in `curve` order whose `auc` reaches `bar`   -> status = :measured
    tau = nothing, when no delta on the curve reaches it               -> status = :below_resolution

`curve` need only carry `delta` and `auc` per entry, so a synthetic curve can exercise the
stopping rule without spending a single simulation.

# THE GRID IS NEVER EXTENDED ON A MISS. THIS IS NOT NEGOTIABLE.

`P13_TAU_ABORT_EXTEND_GRID = false` is the frozen constant that says so, and it is read here
rather than described. A `:below_resolution` return is a REPORTABLE FINDING, not a failure to
be worked around: it says the fixed patch-correlation summary cannot reliably distinguish ANY
pre-registered rho difference, and therefore that Phase 13 cannot draw a defensible three-way
cut at this summary's resolution. Adding a larger delta AFTER SEEING that result would make
"below resolution" UNFALSIFIABLE -- the grid would always eventually contain a delta that
passes, so the abort criterion would never be able to fire and would stop carrying
information. THE ONE THRESHOLD THAT MUST BE IMMUNE TO ITS OWN RESULT IS THIS ONE. There is
deliberately no code path in this file that appends to, widens or regenerates the grid, and
`spike/test/test_p13_tau.jl` greps the executable source to keep it that way.

# The sub-knot caution

When the measured tau falls BELOW the `ghat` knot spacing near zero the result is PRINTED
with a caution (`sub_knot = true`) and NOT turned into an error. ghat is a piecewise-linear
interpolation of a measured sweep, so a tau finer than its knot spacing would be claiming a
resolution the calibration map itself does not carry, which is 13-RESEARCH Pitfall 8's stated
warning sign. It is a caution rather than a failure because the finding may be real and the
right response is to say so in the report, not to silently round tau up to a knot.
"""
function tau_from_curve(curve; bar::Real = P13_TAU_AUC, warn_sub_knot::Bool = true)
    isempty(curve) && throw(ArgumentError(
        "tau_from_curve: the curve is empty; nothing was measured."))
    spacing = ghat_knot_spacing_near_zero()
    for (i, row) in enumerate(curve)
        a = row.auc
        if isfinite(a) && a >= bar
            sub_knot = row.delta < spacing
            if sub_knot && warn_sub_knot
                println("NOTE (13-RESEARCH Pitfall 8): the measured tau = $(row.delta) is " *
                        "BELOW the ghat knot spacing near zero ($(spacing)). That claims a " *
                        "resolution finer than the piecewise-linear calibration map itself " *
                        "carries. Reported, not corrected -- state it in the report.")
            end
            return (tau = row.delta, status = :measured, delta_index = i, bar = bar,
                    sub_knot = sub_knot)
        end
    end
    # THE ABORT PATH. Read the frozen constant rather than describing it, so a future edit
    # that flipped it would change behaviour visibly instead of leaving a stale comment.
    P13_TAU_ABORT_EXTEND_GRID && error(
        "tau_from_curve: P13_TAU_ABORT_EXTEND_GRID is true, but no grid-extension path " *
        "exists in this file and none may be added (D-04 / T-13-17). Extending the grid " *
        "after seeing a miss would make :below_resolution unfalsifiable.")
    println("ABORT CRITERION FIRED: no delta on the pre-registered grid reached the bar " *
            "$(bar). tau is BELOW THIS SUMMARY'S RESOLUTION. The grid is NOT extended " *
            "(P13_TAU_ABORT_EXTEND_GRID = false); this is a reportable finding, not a " *
            "failure to be worked around.")
    return (tau = nothing, status = :below_resolution, delta_index = nothing, bar = bar,
            sub_knot = false)
end

# =====================================================================================
# 7. Simulator provenance
# =====================================================================================

"""
    simulator_provenance_guard() -> NamedTuple

`(post_p11, reason, has_chromatic_eps, theta_arity, shift_halfwidth, expected_halfwidth)`.

REPORTS -- never throws, and never stops anything. The CALLER that turns a `false` into a
hard stop is plan 13-10's reported runner, immediately before it spends the reserved stream;
keeping the check reportable is what lets a unit suite exercise it and lets the runner quote
the reason in its own error message.

# Why freezing a tau measured against the PRE-Phase-11 simulator is a pre-registration error

Phase 11 widens the registration prior and adds a chromatic term, so the POST-Phase-11
deployment distribution is strictly WIDER than the pre-Phase-11 one. A wider nuisance joint
can only make two rho values HARDER to tell apart, so a tau measured against the old
simulator is OPTIMISTIC: it would set the random band too narrow, starve the random class,
and hand the confusion matrix a flattering result that came from measuring an easier world
than the net trains in (13-RESEARCH Pitfall 10). `P13_TAU_REQUIRE_POST_P11_SIMULATOR = true`
is the frozen record that the probe blocks on the SIMULATOR merge -- not on Phase 11's
training run, which the probe does not need at all.

# The two markers

  * `sample_prior` returns EIGHT fields including the chromatic epsilon (Phase 11 D-09
    appends it LAST so existing positional theta indices stay valid).
  * the registration-shift prior half-width equals the widened value, checked against
    `P13_TAU_REFERENCE_LAMBDA_EXPECTED` -- the SAME frozen number the reference-lambda rule
    is stated in, so the two cannot drift apart. The pre-Phase-11 half-width was 1.0.

The probe theta is drawn from `p13_fix_rng(P13_FIXTURE_COUNTER)`, so asking the question
never consumes a reported stream.
"""
function simulator_provenance_guard()
    try
        theta = sample_prior(p13_fix_rng(P13_FIXTURE_COUNTER))
        has_eps = hasproperty(theta, :chromatic_eps)
        arity = length(theta)
        hw = (maximum(SHIFT_PRIOR) - minimum(SHIFT_PRIOR)) / 2
        expected = float(P13_TAU_REFERENCE_LAMBDA_EXPECTED)
        widened = isapprox(hw, expected; rtol = 1e-9)
        post = has_eps && arity == 8 && widened
        reason = post ?
            "POST-Phase-11 simulator confirmed: theta arity $arity with chromatic_eps " *
            "present, shift-prior half-width $hw == expected $expected." :
            "PRE-Phase-11 simulator markers (13-RESEARCH Pitfall 10): " *
            "chromatic_eps present = $has_eps, theta arity = $arity (want 8), " *
            "shift-prior half-width = $hw (want $expected). Freezing a tau measured here " *
            "would be OPTIMISTIC because the deployment distribution is wider."
        return (post_p11 = post, reason = reason, has_chromatic_eps = has_eps,
                theta_arity = arity, shift_halfwidth = hw, expected_halfwidth = expected)
    catch err
        # A guard that threw would be useless to its caller: it must be able to SAY that the
        # simulator surface is not the expected one, including when asking crashes.
        return (post_p11 = false,
                reason = "simulator_provenance_guard could not read the simulator surface " *
                         "at all: $(sprint(showerror, err)). Treat as PRE-Phase-11.",
                has_chromatic_eps = false, theta_arity = -1,
                shift_halfwidth = NaN,
                expected_halfwidth = float(P13_TAU_REFERENCE_LAMBDA_EXPECTED))
    end
end

# =====================================================================================
# 8. Script entry point -- GUARDED, and it runs NOTHING
# =====================================================================================
# GUARDED so `include`-ing this file (from the unit suite, from a runner, or from a REPL)
# never triggers a probe run. The reported probe is a separate, single-purpose script that
# also records the simulator provenance and persists the whole curve; running the probe as a
# side effect of loading its machinery would spend the reserved stream by accident.
if abspath(PROGRAM_FILE) == @__FILE__
    println("spike/p13/tau_probe.jl holds the D-06 probe machinery and runs nothing. " *
            "The REPORTED probe is spike/p13/run_tau_probe.jl (plan 13-10), which must " *
            "first confirm simulator_provenance_guard().post_p11.")
end
