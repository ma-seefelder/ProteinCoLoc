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

# spike/p13/labels.jl --- Phase-13 three-way class boundary and per-head targets (D-05, D-07, D-11).
#
# THIS FILE HOLDS THE TWO PLACES THIS PHASE IS MOST LIKELY TO BE SILENTLY WRONG. Both
# traps look identical to correct code at a glance, both are cheap to encode as an
# executable assertion, and both are expensive to discover after a training run.
#
# (a) WHY THE CUT IS ON THE rho_sample LEVEL TIMES THE CONTROL CONTRAST, AND NOT ON THE
#     SAMPLE-MINUS-CONTROL CONTRAST ALONE (D-05).
#     The natural handle is the quantity the shipped binary net already uses: delta_rho =
#     mean(rho_sample_draws .- rho_control_draws), src/amortized/infer.jl:117-121. It is a
#     sample-vs-control CONTRAST and nothing more. THE COUNTER-EXAMPLE THAT MAKES THIS
#     DECISION LOAD-BEARING: a sample at rho_s = 0.8 under a control at rho_c = 0.9 has
#     delta_rho = -0.1 < 0. The binary net calls that "null", which is harmless. A naive
#     three-way extension of the SAME quantity would publish it as "MUTUALLY EXCLUSIVE"
#     while the sample is strongly COLOCALIZED. "Less colocalized than the control" is not
#     segregation. The fix is the two-factor cut of variant (a) (P13_CUT_VARIANT =
#     :tau_contrast): the LEVEL factor (rho_s < -tau) has to agree with the CONTRAST factor
#     (rho_s - rho_c < -tau) before "exclusion" may be claimed. Both factors are stated at
#     ONE tau -- the summary's MEASURED resolution -- which is D-06's entire point, and the
#     rule is sign-symmetric, so the coloc and exclusion classes are mirror images.
#
# (b) WHY EACH HEAD TRAINS ONLY ON ITS OWN CLASS PAIR (D-11).
#     One-vs-rest is the default reflex and it is wrong here. If the coloc head also trained
#     on the exclusion class as a negative, its Bayes-optimal logit would be a ratio against
#     the MIXTURE q_R*R + q_E*E -- a mixture Bayes factor -- not log BF(coloc : random). The
#     code would look identical and D-08's published semantics would be silently broken. The
#     four-row target vector from `head_targets` encodes the restriction as PARTICIPATION
#     WEIGHTS: the random class is the negative for BOTH heads (which is what makes the two
#     logits ratios against the SAME reference, D-08), the coloc class switches the exclusion
#     head OFF, and the exclusion class switches the coloc head OFF. The TRUNK still trains
#     jointly on all three classes; only the head losses are restricted.
#
#     Nothing here invents a model prior over the three hypotheses (D-08). The random class
#     is a shared REFERENCE, not a prior mass, and `measure_head_log_odds` measures a
#     TRAINING-SET frequency ratio -- never a belief about which hypothesis is true.
#
# (c) TAU IS TIER-2 AND MAY NOT BE INLINED (D-01/D-06). Every function here takes tau as a
#     keyword defaulting to `p13_tau()`, which ERRORS until plan 13-10 appends the MEASURED
#     value to spike/p13/consts.jl. So labelled three-way data cannot be generated before
#     tau has been measured and frozen: the tripwire is the default argument, not a comment.
#     Tests pass tau explicitly, which is why the suite is runnable before the Tier-2 commit.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reaches `src/` only read-only. The
# analog is src/amortized/train_ratio.jl:82-134 and its spike-lane twin
# spike/validation/train_ratio.jl:231-262; neither is edited by this phase.

using Statistics    # mean (empirical class fractions)
using Random        # AbstractRNG (the threaded counter-based stream)

# ORDER MATTERS: the pre-registration first (p13_tau, P13_CUT_VARIANT, P13_DATAGEN_COUNTER,
# p13_rng), then the simulator prior (sample_prior). Guarded for idempotency under
# runtests.jl, the house guarded-include idiom.
isdefined(@__MODULE__, :P13_DEV_SEED) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :sample_prior) || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))

# The ordering is exclusion / random / coloc so that `Int` order matches INCREASING rho,
# which is what makes the D-12 confusion matrix read left-to-right (most anti-correlated in
# the first row/column, most colocalized in the last) instead of in an arbitrary order.
if !isdefined(@__MODULE__, :ThreeWayClass)
    @enum ThreeWayClass EXCLUSION RANDOM COLOC
end

# --- (a) The D-05 boundary ------------------------------------------------------------------
#
# The full rule, restated in the pre-registration's own words (spike/p13/consts.jl section C):
#     COLOC      iff  rho_s >  tau  AND  rho_s - rho_c >  tau
#     EXCLUSION  iff  rho_s < -tau  AND  rho_s - rho_c < -tau
#     RANDOM     otherwise
# The rejected variant (b) tests the contrast with a strict inequality and NO dead zone, so a
# 0.001 contrast would count as "stands apart" -- label noise placed exactly where the heads
# must be calibrated. The frozen choice is recorded as P13_CUT_VARIANT = :tau_contrast; this
# function implements that constant and nothing else.
"""
    three_way_label(rho_sample, rho_control; tau = p13_tau()) -> ThreeWayClass

The D-05 three-way class boundary, variant (a) ("tau-contrast"): the LEVEL of `rho_sample`
AND its CONTRAST against `rho_control` must BOTH clear `tau`, with the same `tau`, and with
the same sign.

    COLOC      iff  rho_sample >  tau  and  rho_sample - rho_control >  tau
    EXCLUSION  iff  rho_sample < -tau  and  rho_sample - rho_control < -tau
    RANDOM     otherwise (the dead zone absorbs everything ambiguous)

DELIBERATELY NOT USED: the sample-minus-control contrast Δρ of
`src/amortized/infer.jl:117-121`. On its own it is not a segregation statement.

    three_way_label(0.8, 0.9; tau = 0.10)  ===  RANDOM        # NEVER EXCLUSION

A sample at ρ = 0.8 under a control at ρ = 0.9 has Δρ = -0.1 < 0. Publishing that as
"mutually exclusive" while the sample is strongly colocalized is the failure this two-factor
cut exists to prevent: "less colocalized than the control" is not segregation.

The rule is sign-symmetric: negating BOTH arguments maps COLOC to EXCLUSION, EXCLUSION to
COLOC and RANDOM to RANDOM. `tau` is the summary's MEASURED resolution (D-06) and is Tier-2 --
the default `p13_tau()` ERRORS until plan 13-10 appends the measured value, so no labelled
datum can exist before tau is frozen. `P13_CUT_VARIANT` records the frozen variant.

Errors with an `ArgumentError` on non-positive `tau` (a zero or negative dead zone would
silently turn variant (a) into the rejected variant (b)).
"""
function three_way_label(rho_sample::Real, rho_control::Real; tau::Real = p13_tau())
    tau > 0 || throw(ArgumentError(
        "three_way_label: tau must be > 0 (got $tau); a non-positive dead zone degenerates " *
        "variant (a) into the rejected sign-contrast variant (b)"))
    d = rho_sample - rho_control
    (rho_sample >  tau && d >  tau) && return COLOC
    (rho_sample < -tau && d < -tau) && return EXCLUSION
    return RANDOM
end

# --- (b) The D-11 per-head participation weights -------------------------------------------
"""
    head_targets(c::ThreeWayClass) -> Vector{Float32}

The 4-row per-sample training target of D-11, with row legend `[y_C; w_C; y_E; w_E]`:
`y_C`/`y_E` are the BCE targets of the coloc and exclusion heads, `w_C`/`w_E` their
PARTICIPATION WEIGHTS (1 = this head trains on this sample, 0 = this head ignores it).

    COLOC     -> Float32[1, 1, 0, 0]    # coloc head positive; exclusion head INACTIVE
    EXCLUSION -> Float32[0, 0, 1, 1]    # exclusion head positive; coloc head INACTIVE
    RANDOM    -> Float32[0, 1, 0, 1]    # the shared NEGATIVE class of BOTH heads

WHY THE RESTRICTION IS THE POINT (D-11). Because RANDOM is the negative class for BOTH heads,
each head's Bayes-optimal logit is a ratio against the SAME reference, which is exactly what
makes the two published numbers `log BF(coloc : random)` and `log BF(exclusion : random)`
(D-08). Under a one-vs-rest scheme (`w_C == 1` everywhere) the coloc head's negative class
would be the MIXTURE `q_R*R + q_E*E` and its logit a MIXTURE Bayes factor -- a different
quantity, published under D-08's name, produced by code that looks identical.

The TRUNK still sees all three classes (both head losses are summed), which is what makes the
two evidences structurally consistent (D-09). Only the head losses are restricted.
"""
function head_targets(c::ThreeWayClass)
    c === COLOC     && return Float32[1, 1, 0, 0]
    c === EXCLUSION && return Float32[0, 0, 1, 1]
    return Float32[0, 1, 0, 1]
end

"""
    target_matrix(classes) -> Matrix{Float32}

Assemble the per-batch training target from a vector of `ThreeWayClass` labels.

CONTRACT: `size(target_matrix(classes)) == (4, length(classes))`, column `j` being
`head_targets(classes[j])` and the rows being `[y_C; w_C; y_E; w_E]` in that order. This is
the layout the masked two-term BCE loss of `spike/p13/net.jl` indexes positionally, so the
row order is part of the contract, not a convention.

Errors with an `ArgumentError` on an empty `classes` (a 4x0 target is never what the caller
meant, and `reduce(hcat, [])` would fail with an unhelpful message).
"""
function target_matrix(classes::AbstractVector{ThreeWayClass})
    isempty(classes) && throw(ArgumentError("target_matrix: classes must be non-empty"))
    return reduce(hcat, head_targets.(classes))::Matrix{Float32}
end

# --- (c) The D-07 per-head evidence-scale correction ---------------------------------------
#
# THE ONE THING NOT TO COPY (13-RESEARCH Pitfall 1). The binary read surface subtracts a
# single measured term -- `measure_log_prior_odds(model_index)`, src/amortized/train_ratio.jl
# :91-96, consumed at src/amortized/bf.jl:74 -- over ALL training labels. That term belongs to
# a DIFFERENT construction: NeuralEstimators' `RatioEstimator` separates the joint from the
# product of marginals using a SHUFFLED theta with equal counts (RatioEstimator.jl:94-122), so
# its logit DIFFERENCE is already the log Bayes factor and the subtracted `log_prior_odds` is a
# numerically negligible (~0.01 nat) redundancy. Because it is ~0 in the binary case, a wrong
# copy -- wrong sign, wrong denominator, or applied twice -- is INVISIBLE there. It is not
# invisible here: at tau = 0.10 the pi-level per-head terms are -0.5333 and -0.1654, i.e.
# 15-50x larger. The quantity below is derived from scratch in 13-RESEARCH F1/F4 and counted
# over the head's OWN pair.
"""
    measure_head_log_odds(class_labels; positive, negative) -> Float64

The MEASURED per-head evidence-scale correction of D-07: `log(n_pos / n_neg)`, counted
**only** over the head's OWN two classes.

    n_pos = count(==(positive), class_labels)
    n_neg = count(==(negative), class_labels)

DERIVATION (13-RESEARCH F1). For a plain BCE head separating class A from class B with
TRAINING frequencies q_A, q_B (q_A + q_B = 1 within that head's restricted set, per D-11),
the Bayes-optimal logit is

    s_h(Z) = log[ p_q(Z|A) / p_q(Z|B) ] + log(q_A / q_B)

so the log Bayes factor is the logit MINUS `log(q_A / q_B)`. That is the scalar this function
measures. It is MEASURED, never assumed 0: the three-way class masses under π are far from
balanced, so assuming 0 would shift every published log BF by up to ~0.53 nat.

NOT THE BINARY CASE, AND THE BINARY FORM MUST NOT BE COPIED (Pitfall 1). `RatioEstimator`'s
shuffled-θ construction already CONTAINS the ratio, so the binary read surface's
`− log prior-odds` term (`src/amortized/train_ratio.jl:91-96`, applied at
`src/amortized/bf.jl:74`) is a numerically negligible redundancy over a DIFFERENT denominator:
ONE scalar over ALL labels. Here there are TWO scalars, each over its head's restricted pair.
Copying the binary form gives the wrong denominator (all three classes instead of the head's
pair) and the error is ~0.5 nat, not ~0.01.

EXACTNESS CONDITION (13-RESEARCH F2). This scalar recovers the π-scale Bayes factor if and
only if stratification changed the class FREQUENCIES only, never the WITHIN-class θ shape --
because a within-class distributional change enters INSIDE the evidence integral and depends
on Z, so no scalar can repair it. `P13_STRATIFICATION = :class_frequency` is the design that
satisfies the condition; `P13_STRATIFICATION_FALLBACK = :importance_weighted` is the
pre-declared route if it must be relaxed.

Returns a finite scalar; errors on degenerate head balance (either class absent), in the exact
style of the binary analog it generalizes.
"""
function measure_head_log_odds(class_labels::AbstractVector{ThreeWayClass};
                               positive::ThreeWayClass, negative::ThreeWayClass)
    positive === negative && throw(ArgumentError(
        "measure_head_log_odds: positive and negative must differ (both $positive); a head " *
        "cannot be its own reference"))
    n_pos = count(==(positive), class_labels)
    n_neg = count(==(negative), class_labels)
    (n_pos > 0 && n_neg > 0) ||
        error("measure_head_log_odds: degenerate head balance n_pos=$n_pos n_neg=$n_neg " *
              "(head $positive vs $negative; need both > 0)")
    return log(n_pos / n_neg)
end

"""
    head_log_odds(class_labels) -> NamedTuple

Both D-07 corrections at once, in the EXACT shape the read surface of `spike/p13/net.jl`
consumes:

    (coloc     = measure_head_log_odds(class_labels; positive = COLOC,     negative = RANDOM),
     exclusion = measure_head_log_odds(class_labels; positive = EXCLUSION, negative = RANDOM))

RANDOM is the negative of both entries -- the shared reference of D-08 -- and each entry is
counted over its own head's pair only, so the coloc entry's denominator is the RANDOM count
and never the coloc complement.
"""
head_log_odds(class_labels::AbstractVector{ThreeWayClass}) = (
    coloc     = measure_head_log_odds(class_labels; positive = COLOC,     negative = RANDOM),
    exclusion = measure_head_log_odds(class_labels; positive = EXCLUSION, negative = RANDOM),
)

# --- Cheap prior-only labelling -------------------------------------------------------------
"""
    prior_class_draws(n; rng = p13_rng(P13_DATAGEN_COUNTER), tau = p13_tau()) -> Vector{ThreeWayClass}

`n` three-way class labels drawn from π WITHOUT simulating a single image -- the three-way
generalization of `prior_label_draws` (`spike/validation/train_ratio.jl:246-262`).

The D-05 label depends ONLY on θ (two independent `sample_prior(rng).ρ_true` draws, sample and
control), so it needs no forward simulation. That is how the per-head log-odds of
`head_log_odds` and the class-mass table are measured at large `n` for negligible cost.

`rng` is THREADED explicitly through the whole loop (one counter-based stream, consumed
sequentially). `Random.seed!` is NEVER called inside the loop: re-seeding per index would
correlate the draws and destroy the reproducibility the counter-based stream provides
(P13_GLOBAL_RNG_DISCIPLINE).

The `tau` default is the Tier-2 tripwire -- this function ERRORS before drawing anything until
plan 13-10 has appended the MEASURED tau, which is what makes "no labelled datum before tau
is frozen" (D-06) a property of the code rather than of a reviewer's diligence.
"""
function prior_class_draws(n::Integer;
                           rng::AbstractRNG = p13_rng(P13_DATAGEN_COUNTER),
                           tau::Real = p13_tau())
    n >= 1 || throw(ArgumentError("prior_class_draws: n must be >= 1 (got $n)"))
    classes = Vector{ThreeWayClass}(undef, n)
    for j in 1:n
        rho_s = sample_prior(rng).ρ_true
        rho_c = sample_prior(rng).ρ_true
        classes[j] = three_way_label(rho_s, rho_c; tau = tau)
    end
    return classes
end

"""
    class_masses(classes) -> NamedTuple

Empirical class fractions `(exclusion = ..., random = ..., coloc = ...)`, summing to 1.

Used by `spike/test/test_p13_labels.jl` to reproduce the frozen class-mass table of
13-RESEARCH E1 in-repo, and by the datagen plan to verify the realized frequencies against
the `P13_TARGET_CLASS_FREQ` stratification target. Errors with an `ArgumentError` on an empty
`classes` (a mass table over zero samples is meaningless, and the fractions would be `NaN`).
"""
function class_masses(classes::AbstractVector{ThreeWayClass})
    isempty(classes) && throw(ArgumentError("class_masses: classes must be non-empty"))
    return (exclusion = mean(c -> c === EXCLUSION, classes),
            random    = mean(c -> c === RANDOM,    classes),
            coloc     = mean(c -> c === COLOC,     classes))
end

# D-05's counter-example, as an executable regression assertion (13-RESEARCH Pitfall 3: the
# warning sign is an exclusion-class member with rho_sample > 0). Evaluated at an EXPLICIT tau
# so this file stays loadable before the Tier-2 tau commit.
@assert three_way_label(0.8, 0.9; tau = 0.10) != EXCLUSION  # "less colocalized than the control" is NOT segregation
