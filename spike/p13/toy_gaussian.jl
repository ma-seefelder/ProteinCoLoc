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

# spike/p13/toy_gaussian.jl --- closed-form verification of the D-07 per-head correction.
#
# WHY THIS FILE EXISTS. D-07 does not say "apply the binary prior-odds correction to two heads".
# It says the three-way correction must be DERIVED AND VERIFIED, and it says so because the
# binary formula does NOT generalize by analogy: the shipped binary net's subtracted term belongs
# to NeuralEstimators' shuffled-theta construction, where the logit DIFFERENCE already contains
# the ratio and the term is a ~0.01-nat redundancy. The D-09 heads are plain BCE classifiers,
# for which the correction is genuinely load-bearing (15-50x larger) AND is exact only under a
# condition -- 13-RESEARCH F2's theorem -- that a scalar cannot check on itself.
#
# Everything else in this phase compares one ESTIMATE against another ESTIMATE. The
# conjugate-Gaussian construction below is the only object in this repository whose three-way log
# Bayes factors are ANALYTIC, so it is the only place where the implemented correction can be
# checked against TRUTH. That is the whole reason to carry a toy.
#
# THE TWO NEGATIVE CONTROLS, AND WHAT EACH ONE RULES OUT.
#
#   CONTROL 1 -- the UNCORRECTED logit must FAIL the same bars.
#     Rules out a VACUOUS PASS. A verification in which the corrected numbers clear the bars but
#     the uncorrected numbers would have cleared them too proves only that the toy's class
#     frequencies were nearly balanced; it cannot distinguish "the correction works" from "the
#     correction was unnecessary here". `P13_F5_SKEW_FREQ` makes the frequencies deliberately
#     lopsided so the omission is impossible to miss, and the control asserts the failure has the
#     SPECIFIC magnitude of the head's own log-odds -- so a correction with the wrong sign, the
#     wrong denominator, or applied twice fails the control too.
#
#   CONTROL 2 -- the scalar must FAIL under a WITHIN-CLASS reshaped proposal.
#     Rules out an ASSERTED design choice. `P13_STRATIFICATION = :class_frequency` was chosen
#     over within-class stratification on the strength of 13-RESEARCH F2's theorem, which says
#     the scalar recovers the pi-scale Bayes factor IFF stratification changed class FREQUENCIES
#     and never the within-class shape. `toy_reshaped_dataset` is that theorem's counter-example
#     made runnable: identical class frequencies, one class's within-class shape replaced, and
#     the same scalar correction applied. Its failure is EXPECTED and is the empirical evidence
#     for the pre-registered design, not a bug.
#
# NO QUADRATURE AND NO KDE. The class-conditional evidence of a Gaussian likelihood under an
# interval-truncated Gaussian prior is closed form (two normal CDFs over two normal CDFs), so
# nothing here integrates numerically. `QuadGK` and `KernelDensity` are NOT in the spike
# environment and the KDE Bayes-factor baseline is retired (named limit 3); reaching for either
# would mean installing a package, which this phase does not do.
#
# THE TOY IS NOT A MODEL OF ANYTHING. theta is the rho analogue and Z the summary analogue, and
# the resemblance stops there: there is no control channel, no image, no nuisance. Its only job
# is to be a construction in which the correction's target value is KNOWN.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. Touches no src/ byte, adds no dependency,
# installs nothing, never imports CUDA. Guarded includes keep the file loadable standalone AND
# idempotent under a runner. This file has no script entry point: it is included by
# spike/test/test_p13_correction.jl and by nothing else.

import Distributions: Normal, cdf, logpdf   # the closed-form evidence: two CDFs over two CDFs
using Statistics                            # quantile (the central-band restriction), cor, mean
using Random                                # AbstractRNG threading + randperm (the class shuffle)

# ORDER MATTERS: the frozen pre-registration first (P13_F5_*, P13_FIX_SEED, P13_VAL_FRAC), then
# the label surface (ThreeWayClass, target_matrix, head_log_odds), then the net (the REAL
# architecture and the REAL masked loss this verification must exercise). Guarded for
# idempotency, the house guarded-include idiom.
isdefined(@__MODULE__, :P13_DEV_SEED)       || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :three_way_label)    || include(joinpath(@__DIR__, "labels.jl"))
isdefined(@__MODULE__, :ThreeWayEvidenceNet) || include(joinpath(@__DIR__, "net.jl"))

# --- 1. The closed-form class-conditional evidence ------------------------------------------
"""
    toy_log_evidence(Z, a, b; s, sigma) -> Float64

`log p(Z | theta in [a, b])` for the conjugate toy `theta ~ N(0, s^2)`, `Z | theta ~ N(theta, sigma^2)`.

DERIVATION, in three lines. Marginally `Z ~ N(0, s^2 + sigma^2)` and the posterior is
`theta | Z ~ N(m, v)` with `m = Z*s^2/(s^2+sigma^2)` and `v = s^2*sigma^2/(s^2+sigma^2)`.
Restricting theta to `[a, b]` multiplies the joint by the prior mass of that interval and the
marginal by the POSTERIOR mass of it, so

    p(Z | theta in [a,b]) = N(Z; 0, s^2+sigma^2) * P(theta in [a,b] | Z) / P(theta in [a,b])

Open-ended intervals are passed as `-Inf` / `Inf` and handled by `cdf` directly; no special case
and no truncation constant is needed.

THE IDENTITY THIS SATISFIES, and the sanity check the suite exercises: the numerators
`P(theta in [a,b] | Z)` of a PARTITION of the real line sum to one, so the class evidences mix
back to the marginal,

    sum_over_classes  P(class) * p(Z | class)  ==  N(Z; 0, s^2 + sigma^2)

exactly, for every `Z`. A closed form that merely looks plausible would not satisfy that.

Errors on a non-positive `s` or `sigma` and on an empty or inverted interval.
"""
function toy_log_evidence(Z::Real, a::Real, b::Real; s::Real, sigma::Real)
    s > 0 || throw(ArgumentError("toy_log_evidence: s must be > 0 (got $s)"))
    sigma > 0 || throw(ArgumentError("toy_log_evidence: sigma must be > 0 (got $sigma)"))
    b > a || throw(ArgumentError("toy_log_evidence: need b > a (got a=$a, b=$b)"))
    t  = s^2 + sigma^2
    m  = Z * s^2 / t
    sd = sqrt(s^2 * sigma^2 / t)
    post = Normal(m, sd)
    pri  = Normal(0.0, float(s))
    num = cdf(post, b) - cdf(post, a)     # P(theta in [a,b] | Z)
    den = cdf(pri, b)  - cdf(pri, a)      # P(theta in [a,b])
    return logpdf(Normal(0.0, sqrt(t)), Z) + log(num) - log(den)
end

"""
    toy_prior_class_masses(; tau, s) -> NamedTuple

The exact prior class masses `(exclusion, random, coloc)` of the toy's level cut under
`theta ~ N(0, s^2)`: `P(theta < -tau)`, `P(|theta| <= tau)`, `P(theta > tau)`. They sum to one.

These are the mixing weights of the identity documented on `toy_log_evidence`, and they are also
the reference the DELIBERATELY SKEWED `P13_F5_SKEW_FREQ` departs from -- the whole point of the
toy being that the training frequencies are NOT these.
"""
function toy_prior_class_masses(; tau::Real, s::Real)
    tau > 0 || throw(ArgumentError("toy_prior_class_masses: tau must be > 0 (got $tau)"))
    s > 0 || throw(ArgumentError("toy_prior_class_masses: s must be > 0 (got $s)"))
    pri = Normal(0.0, float(s))
    lo  = cdf(pri, -tau)
    hi  = cdf(pri, tau)
    return (exclusion = lo, random = hi - lo, coloc = 1.0 - hi)
end

# --- 2. The analytic three-way log Bayes factors ---------------------------------------------
"""
    toy_analytic_log_bf(Z; tau, s, sigma) -> NamedTuple

The EXACT three-way log Bayes factors against the random reference at summary `Z`:

    (coloc     = log p(Z | theta >  tau) - log p(Z | |theta| <= tau),
     exclusion = log p(Z | theta < -tau) - log p(Z | |theta| <= tau),
     random    = 0.0)

The shape is DELIBERATELY the one `three_way_log_bf` (`spike/p13/net.jl`) returns, field for
field and in the same order, so the comparison against the network is elementwise on named
fields. A truth vector in a different container would invite a reshape step, and a silent
transpose inside a reshape step is exactly the class of bug this verification exists to catch.

`random` is identically `0.0` for the same structural reason it is in the read surface (D-08):
both reported numbers are ratios AGAINST the random reference, so the reference against itself
is zero by construction, never by measurement.

SIGN-SYMMETRY IDENTITY. The toy's prior and likelihood are both symmetric about zero and the cut
is symmetric in `tau`, so for every `Z`

    toy_analytic_log_bf(Z).coloc  ==  toy_analytic_log_bf(-Z).exclusion

-- reflecting `Z` swaps the coloc and exclusion evidences and leaves the random one alone. It is
an equality between the two HEADS at reflected arguments, NOT a negation of either: at `tau =
0.3, s = 1, sigma = 0.5` and `Z = 2` the coloc entry is about `+5.83` while the exclusion entry
at `Z = -2` is the same `+5.83`, and neither is the other's negative (the exclusion entry at
`Z = +2` is about `-5.62`, a different number). The suite asserts the equality form.
"""
function toy_analytic_log_bf(Z::Real; tau::Real, s::Real, sigma::Real)
    tau > 0 || throw(ArgumentError("toy_analytic_log_bf: tau must be > 0 (got $tau)"))
    lR = toy_log_evidence(Z, -tau,  tau; s = s, sigma = sigma)
    lC = toy_log_evidence(Z,  tau,  Inf; s = s, sigma = sigma)
    lE = toy_log_evidence(Z, -Inf, -tau; s = s, sigma = sigma)
    return (coloc = lC - lR, exclusion = lE - lR, random = 0.0)
end

# --- 3. The toy's class boundary --------------------------------------------------------------
"""
    toy_class_of(theta; tau) -> ThreeWayClass

The toy's class of a latent `theta` under the plain level cut

    theta >  tau -> COLOC
    theta < -tau -> EXCLUSION
    otherwise    -> RANDOM

THIS IS THE LEVEL FACTOR OF D-05 ONLY, AND DELIBERATELY SO. The real cut
(`three_way_label`, `spike/p13/labels.jl`) is TWO-factor -- the level of `rho_sample` AND its
contrast against a control channel must both clear `tau` -- because "less colocalized than the
control" is not segregation. The toy has no control channel at all, so the contrast factor has
no meaning here and is not modelled. The two-factor cut is verified in its own suite
(`spike/test/test_p13_labels.jl`); this file verifies the EVIDENCE-SCALE CORRECTION and nothing
else, and it borrows only the class ENUM from the label surface.

Returns the same `ThreeWayClass` enum the real labeller returns, so `target_matrix` and
`head_log_odds` are consumed unmodified rather than re-implemented for the toy.
"""
function toy_class_of(theta::Real; tau::Real)
    tau > 0 || throw(ArgumentError("toy_class_of: tau must be > 0 (got $tau)"))
    theta >  tau && return COLOC
    theta < -tau && return EXCLUSION
    return RANDOM
end

# Split `n` into per-class counts at the requested frequencies, giving the remainder to the
# largest class so the realized frequencies stay as close to the request as integer counts allow.
function _toy_class_counts(n::Integer, freq)
    length(freq) == 3 || throw(ArgumentError(
        "_toy_class_counts: freq must be a 3-tuple (exclusion, random, coloc), got $freq"))
    all(f -> f > 0, freq) || throw(ArgumentError(
        "_toy_class_counts: every class frequency must be > 0 (got $freq); an absent class " *
        "would make measure_head_log_odds degenerate"))
    isapprox(sum(freq), 1.0; atol = 1e-9) || throw(ArgumentError(
        "_toy_class_counts: freq must sum to 1 (got $(sum(freq)))"))
    k = collect(floor.(Int, n .* collect(freq)))
    k[argmax(freq)] += n - sum(k)
    return k
end

# Rejection-draw a single theta from pi(theta | class): draw theta ~ N(0, s^2) and keep it only
# if it lands in the requested class. THE ACCEPT/REJECT IS THE POINT, not an implementation
# detail -- it is what makes the within-class shape EXACTLY pi(theta | class) while the class
# COUNTS are set by the caller. That is design D-07-i, and it is the condition under which the
# scalar correction is exact (13-RESEARCH F2).
function _toy_draw_theta(rng::AbstractRNG, c::ThreeWayClass, tau::Real, s::Real;
                         maxtries::Integer = 100_000)
    for _ in 1:maxtries
        th = s * randn(rng)
        toy_class_of(th; tau = tau) === c && return th
    end
    error("_toy_draw_theta: no draw landed in class $c after $maxtries tries " *
          "(tau=$tau, s=$s) -- the class has negligible prior mass at these parameters")
end

# --- 4. The primary dataset: class FREQUENCIES changed, within-class shape preserved ----------
"""
    toy_draw_dataset(n; tau, s, sigma, freq = P13_F5_SKEW_FREQ, rng) -> NamedTuple

`n` toy observations at the requested class frequencies, returned as

    (Z = 1-by-n Float32 matrix, classes = Vector{ThreeWayClass}, theta = Vector{Float64})

`freq` is ordered `(exclusion, random, coloc)` -- the enum's own increasing-rho order, the same
order `class_masses` reports -- and defaults to the frozen `P13_F5_SKEW_FREQ`.

THIS IS DESIGN D-07-i, THE PRIMARY. Whole classes are drawn to the target COUNTS by accept/reject
on `theta ~ N(0, s^2)`, so `pi(theta | class)` is preserved EXACTLY and only the class
frequencies move. That is precisely the condition under which the scalar per-head correction
recovers the pi-scale Bayes factor (13-RESEARCH F2), which is why this is the arm expected to
PASS the pre-registered bars.

WHY THE FREQUENCIES ARE DELIBERATELY LOPSIDED. Under `P13_F5_SKEW_FREQ` the per-head log-odds are
several times the `P13_F5_MAXABS_TOL` tolerance, so a MISSING correction cannot slip through as
noise -- it is the negative control's entire leverage. A near-balanced toy would make the
verification vacuous, which is the failure mode `P13_F5_NEGCTRL_REQUIRED` exists to forbid.

THE RETURNED SAMPLE IS SHUFFLED. Classes are generated in blocks, so an unshuffled return would
put a single class in any contiguous validation split -- a head would then be scored on a
validation set containing none of its own positives. `rng` is THREADED through the theta draws,
the observation noise and the permutation; `Random.seed!` is never called here.
"""
function toy_draw_dataset(n::Integer; tau::Real, s::Real, sigma::Real,
                          freq = P13_F5_SKEW_FREQ, rng::AbstractRNG)
    n >= 3 || throw(ArgumentError("toy_draw_dataset: n must be >= 3 (got $n)"))
    counts  = _toy_class_counts(n, freq)
    theta   = Vector{Float64}(undef, n)
    classes = Vector{ThreeWayClass}(undef, n)
    j = 0
    for (c, k) in zip((EXCLUSION, RANDOM, COLOC), counts)
        for _ in 1:k
            j += 1
            theta[j]   = _toy_draw_theta(rng, c, tau, s)
            classes[j] = c
        end
    end
    noise = sigma .* randn(rng, n)
    perm  = randperm(rng, n)
    theta, classes, noise = theta[perm], classes[perm], noise[perm]
    Z = reshape(Float32.(theta .+ noise), 1, n)
    return (Z = Z, classes = classes, theta = theta)
end

# --- 5. The negative-control dataset: within-class shape RESHAPED ----------------------------
"""
    toy_reshaped_dataset(n; tau, s, sigma, freq = P13_F5_SKEW_FREQ, reshape_halfwidth, rng)

THE SUBSTRATE OF NEGATIVE CONTROL 2. IDENTICAL CLASS FREQUENCIES TO `toy_draw_dataset`, AND A
DELIBERATELY DIFFERENT WITHIN-CLASS SHAPE: the exclusion class's theta is drawn UNIFORMLY over
`[-reshape_halfwidth*s, -tau]` instead of from `pi(theta | EXCLUSION)`, which is the
`(-Inf, -tau)`-truncated normal. The random and coloc classes are untouched, so the two arms
differ in EXACTLY ONE respect and the failure below is attributable to it alone.

THE SCALAR CORRECTION IS EXPECTED TO FAIL ON THIS DATASET, AND THAT FAILURE IS THE RESULT, NOT A
BUG. 13-RESEARCH F2 states the exactness condition as an IFF: the scalar `-log(q_A/q_B)` recovers
the pi-scale Bayes factor if and only if `q(theta|A) == pi(theta|A)` and `q(theta|B) ==
pi(theta|B)`. Here `q(theta|EXCLUSION) != pi(theta|EXCLUSION)`, and the discrepancy enters INSIDE
the evidence integral `p_q(Z|A) = integral p(Z|theta) q(theta|A) dtheta`, where it is a function
of `Z`. No CONSTANT can repair a function of `Z`. The exclusion head's corrected output is
therefore wrong by a Z-dependent amount, while the coloc head -- whose two classes were not
reshaped -- keeps clearing the bars.

THIS IS WHY `P13_STRATIFICATION = :class_frequency` IS THE PRE-REGISTERED PRIMARY DESIGN AND
D-07's literal within-range rho-stratification is NOT. The pre-declared fallback
`P13_STRATIFICATION_FALLBACK = :importance_weighted` exists precisely because a reshaped proposal
needs a per-sample weight `pi(theta)/g(theta)`, not a scalar. Without this dataset that ruling
would rest on the derivation alone; with it, the ruling is evidenced.

`reshape_halfwidth` is expressed in units of `s` and is a property of the CONTROL, not a tuned
knob: it only has to make the uniform proposal visibly unlike the truncated normal.
"""
function toy_reshaped_dataset(n::Integer; tau::Real, s::Real, sigma::Real,
                              freq = P13_F5_SKEW_FREQ, reshape_halfwidth::Real = 4.0,
                              rng::AbstractRNG)
    n >= 3 || throw(ArgumentError("toy_reshaped_dataset: n must be >= 3 (got $n)"))
    reshape_halfwidth * s > tau || throw(ArgumentError(
        "toy_reshaped_dataset: reshape_halfwidth*s must exceed tau, else the reshaped " *
        "exclusion range is empty (got $(reshape_halfwidth*s) vs tau=$tau)"))
    counts  = _toy_class_counts(n, freq)
    theta   = Vector{Float64}(undef, n)
    classes = Vector{ThreeWayClass}(undef, n)
    lo = -reshape_halfwidth * s
    j = 0
    for (c, k) in zip((EXCLUSION, RANDOM, COLOC), counts)
        for _ in 1:k
            j += 1
            # THE ONE DIFFERENCE FROM toy_draw_dataset, and it is confined to this branch.
            theta[j]   = c === EXCLUSION ? lo + (-tau - lo) * rand(rng) :
                                           _toy_draw_theta(rng, c, tau, s)
            classes[j] = c
        end
    end
    noise = sigma .* randn(rng, n)
    perm  = randperm(rng, n)
    theta, classes, noise = theta[perm], classes[perm], noise[perm]
    Z = reshape(Float32.(theta .+ noise), 1, n)
    return (Z = Z, classes = classes, theta = theta)
end

# --- 6. The Bayes-optimal reference logits ---------------------------------------------------
"""
    toy_ideal_logits(Zgrid, head_log_odds; tau, s, sigma) -> Matrix{Float32}

The `2 x n` logit matrix a PERFECT head would emit on `Zgrid`: the analytic log Bayes factor of
each head PLUS that head's measured log-odds, which is the inverse of the D-07 correction

    ideal_logit_h(Z) = log BF_h(Z) + log(q_pos / q_neg)

from 13-RESEARCH F1's `s_h(Z) = log[p_q(Z|A)/p_q(Z|B)] + log(q_A/q_B)`.

WHAT THIS IS FOR: SEPARATING "THE CORRECTION IS WRONG" FROM "THE NETWORK IS NOT CONVERGED".
Scored under `masked_two_head_bce` against the same targets, these logits are the minimum
attainable loss -- the Bayes risk of the classification problem. Comparing a trained net's
held-out loss against it answers, with a number rather than an opinion, whether a residual
deviation from the analytic log Bayes factor is a training shortfall or a property of the
objective. It matters because logit-BCE is nearly FLAT in the confident tail: where the head is
already right with probability near one, a half-nat logit error costs almost nothing in loss, so a
net can sit at the Bayes risk and still be visibly off in nats out there.

Callers must drop any column where a returned logit is not finite: the closed-form evidence
underflows to `-Inf` far enough into a tail, and one such column would poison an aggregate.
"""
function toy_ideal_logits(Zgrid, head_log_odds; tau::Real, s::Real, sigma::Real)
    Zv = collect(Float64.(vec(Zgrid)))
    out = Matrix{Float32}(undef, 2, length(Zv))
    for (j, z) in enumerate(Zv)
        t = toy_analytic_log_bf(z; tau = tau, s = s, sigma = sigma)
        out[1, j] = Float32(t.coloc     + head_log_odds.coloc)
        out[2, j] = Float32(t.exclusion + head_log_odds.exclusion)
    end
    return out
end

# --- 6b. Fitting the REAL net with the REAL masked loss --------------------------------------
"""
    toy_fit_two_head(dataset; input_dim, width, num_summaries, epochs, batchsize,
                     val_frac, seed, verbose) -> NamedTuple

Train the ACTUAL `ThreeWayEvidenceNet` with the ACTUAL `masked_two_head_bce` through the ACTUAL
`train_three_way` on the toy, and return

    (net, head_log_odds, n_train, n_val, classes_train)

THE OBJECT UNDER TEST MUST BE THE OBJECT THAT SHIPS. A stand-in classifier fitted here would
verify an algebraic identity about a different network than the one the phase gates, and the
correction's whole risk (13-RESEARCH Pitfall 1: wrong sign, wrong denominator, applied twice)
lives in the seam between the head and the read surface. `width` and `num_summaries` are smaller
than the pre-registered `P13_SUMMARY_WIDTH`/`P13_NUM_SUMMARIES` because the toy's input is
one-dimensional -- and they are PASSED EXPLICITLY by the caller for that reason, never hard-coded
in here, so nothing in this file can be mistaken for a change to the D-10 recipe.

THE CORRECTION IS MEASURED ON THE TRAINING SUBSET ONLY. `head_log_odds` counts `classes[1:n_train]`
and nothing else. Measuring it on the evaluation set instead would be a leak in the exact sense
that matters here: the correction is a property of the frequencies the head was TRAINED at, so an
evaluation-set count would silently fold an evaluation-set artifact into every reported Bayes
factor. That is threat T-13-23's shape and it is cheap to get wrong.

WHY A `seed` AND NOT AN `rng`. Both stochastic parts of a fit reach for the GLOBAL RNG and cannot
be threaded: Flux's `glorot_uniform` initializer draws from it when the layers are constructed,
and NeuralEstimators builds its `DataLoader` with `shuffle = true`, which draws from it every
epoch. So layer 2 of `P13_GLOBAL_RNG_DISCIPLINE` is the only mechanism available, and it is
applied twice -- `Random.seed!` immediately before the net is BUILT, and again inside
`train_three_way` immediately before `train`. `seed` defaults to the FIXTURE seed, so a unit test
never touches a reported Phase-13 stream.

The validation split is taken from the tail of the (already shuffled) dataset and is floored at
`batchsize`, because NeuralEstimators' loader is built with `partial = false` and DROPS an
incomplete final batch -- a validation set smaller than one batch would yield zero batches.
"""
function toy_fit_two_head(dataset; input_dim::Integer = 1, width::Integer = 32,
                          num_summaries::Integer = 8, epochs::Integer, batchsize::Integer,
                          val_frac::Real = P13_VAL_FRAC, seed = P13_FIX_SEED,
                          verbose::Bool = false)
    Z, classes = dataset.Z, dataset.classes
    size(Z, 1) == input_dim || throw(ArgumentError(
        "toy_fit_two_head: dataset.Z has $(size(Z, 1)) rows but input_dim is $input_dim"))
    n = size(Z, 2)
    n_va = max(batchsize, round(Int, val_frac * n))
    n_tr = n - n_va
    n_tr >= batchsize || throw(ArgumentError(
        "toy_fit_two_head: training split $n_tr is smaller than batchsize $batchsize"))

    classes_tr = classes[1:n_tr]
    classes_va = classes[n_tr+1:end]
    tg_tr = target_matrix(classes_tr)
    tg_va = target_matrix(classes_va)

    Random.seed!(UInt64(seed))          # layer 2, before the Flux initializer draws
    net = build_three_way_net(input_dim; num_summaries = num_summaries, width = width)

    trained = train_three_way(net, tg_tr, tg_va, Z[:, 1:n_tr], Z[:, n_tr+1:end];
                              epochs = epochs, batchsize = batchsize, seed = seed,
                              savepath = nothing, verbose = verbose)

    return (net = trained, head_log_odds = head_log_odds(classes_tr),
            n_train = n_tr, n_val = n_va, classes_train = classes_tr)
end

# --- 7. The per-head evaluation grid -----------------------------------------------------------
"""
    toy_head_grid(dataset, positive::ThreeWayClass) -> Vector{Float64}

The `Z` values of a drawn dataset that belong to ONE head's own class pair: `positive` together
with `RANDOM`, the shared negative of both heads.

EACH HEAD IS SCORED ON ITS OWN RESTRICTED PAIR, AND THAT IS THE PHASE'S PRE-REGISTERED
CONVENTION, not a convenience. It follows from D-11 -- the coloc head never sees an exclusion
sample, so `log BF(coloc : random)` is a statement about the coloc-versus-random comparison and
is only estimable where that comparison has support -- and the frozen consts file already commits
to it for the reported gate, where `P13_GATE_M` is sized so that "each head is scored only on its
own class pair, per D-11" puts enough samples in each RESTRICTED evaluation set.

The consequence for this verification is worth stating plainly: scoring the coloc head on the
pooled three-class sample would evaluate it deep in the exclusion region, where it has no
training support at all and is extrapolating. That would measure the extrapolation behaviour of
a network, not the correctness of a scalar correction, which is what this file is for.
"""
function toy_head_grid(dataset, positive::ThreeWayClass)
    positive === RANDOM && throw(ArgumentError(
        "toy_head_grid: RANDOM is the shared negative, not a head's positive class"))
    idx = findall(c -> c === positive || c === RANDOM, dataset.classes)
    isempty(idx) && throw(ArgumentError("toy_head_grid: no $positive or RANDOM samples present"))
    return collect(Float64.(vec(dataset.Z)[idx]))
end

# --- 8. The comparison, with the corrected/uncorrected switch as a PARAMETER ------------------
"""
    toy_compare(net, head_log_odds, Zgrid; tau, s, sigma, corrected = true,
                central_frac = P13_F5_CENTRAL_FRAC) -> NamedTuple

Compare the network's two log Bayes factors against `toy_analytic_log_bf` over `Zgrid`, returning

    (coloc     = (corr = ..., maxabs = ..., meandev = ..., maxabs_demeaned = ...),
     exclusion = (corr = ..., maxabs = ..., meandev = ..., maxabs_demeaned = ...),
     n = <points scored>, lo = ..., hi = ...)

`corr` is the correlation against the analytic truth, `maxabs` the largest absolute deviation,
and `meandev` the SIGNED mean deviation -- which is what lets the negative control assert the
specific offset rather than merely "something failed".

`maxabs_demeaned` is the largest absolute deviation AFTER the best constant has been removed,
and it exists to state 13-RESEARCH F2's theorem as a number. The theorem says a within-class
distributional change cannot be repaired by ANY scalar, and `maxabs_demeaned` is exactly the
residual that the best scalar in hindsight would leave: if it still exceeds the tolerance, no
correction of the D-07 form -- the measured one, a better-measured one, or an oracle one --
could have brought that head inside the bar. Quoting `maxabs` alone would leave open the reply
"then measure a better constant".

THE `corrected` SWITCH IS A PARAMETER OF THIS FUNCTION, NOT A SECOND COPY OF THE CODE. That is
deliberate and it is the only reason negative control 1 proves anything: the corrected and
uncorrected arms are then PROVABLY the same forward pass, the same grid, the same restriction and
the same statistics, differing in exactly one subtraction. Two hand-written near-duplicates could
drift, and a drifting duplicate makes the control's failure unattributable.

NOTE WHAT THE UNCORRECTED ARM CANNOT FAIL ON. A missing correction is a CONSTANT shift, and
correlation is shift-invariant, so the uncorrected `corr` is identical to the corrected one to
floating-point noise. The uncorrected arm therefore fails on `maxabs` and on `meandev`, and a
verification that gated on correlation alone would have missed the entire defect.

THE CENTRAL-BAND RESTRICTION IS BY QUANTILE OF `Zgrid`, so it presumes `Zgrid` is a SAMPLE from
the evaluation distribution rather than an evenly-spaced ruler: the pre-registered
`P13_F5_CENTRAL_FRAC` then means "the central 90 percent of the observations", trimming the tails
where neither head has training support and where any classifier is extrapolating. Passing an
evenly-spaced grid would silently redefine the band.
"""
# The four deviation statistics of one head, computed once so the corrected and uncorrected arms
# cannot diverge in how they are summarized.
function _toy_dev_stats(pred::AbstractVector, truth::AbstractVector)
    d = pred .- truth
    md = mean(d)
    return (corr = cor(pred, truth), maxabs = maximum(abs.(d)), meandev = md,
            maxabs_demeaned = maximum(abs.(d .- md)))
end

function toy_compare(net::ThreeWayEvidenceNet, head_log_odds, Zgrid;
                     tau::Real, s::Real, sigma::Real, corrected::Bool = true,
                     central_frac::Real = P13_F5_CENTRAL_FRAC)
    0 < central_frac <= 1 || throw(ArgumentError(
        "toy_compare: central_frac must be in (0, 1] (got $central_frac)"))
    Zv = collect(Float64.(vec(Zgrid)))
    length(Zv) >= 2 || throw(ArgumentError("toy_compare: Zgrid needs at least 2 points"))
    tail = (1 - central_frac) / 2
    lo, hi = quantile(Zv, tail), quantile(Zv, 1 - tail)
    Zk = filter(z -> lo <= z <= hi, Zv)

    arch = three_way_arch(net)
    arch.input_dim == 1 || throw(ArgumentError(
        "toy_compare: expects the 1-D toy net (input_dim 1), got $(arch.input_dim)"))
    S = net(reshape(Float32.(Zk), 1, :))

    shiftC = corrected ? head_log_odds.coloc     : 0.0
    shiftE = corrected ? head_log_odds.exclusion : 0.0
    predC = Float64.(S[1, :]) .- shiftC
    predE = Float64.(S[2, :]) .- shiftE

    truth = [toy_analytic_log_bf(z; tau = tau, s = s, sigma = sigma) for z in Zk]
    trueC = [t.coloc for t in truth]
    trueE = [t.exclusion for t in truth]

    return (coloc     = _toy_dev_stats(predC, trueC),
            exclusion = _toy_dev_stats(predE, trueE),
            n = length(Zk), lo = lo, hi = hi)
end
