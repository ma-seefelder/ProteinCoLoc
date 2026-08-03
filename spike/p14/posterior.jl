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

# spike/p14/posterior.jl --- the three-class posterior, composed BY NAME from the MEASURED prior.
#
# THIS FILE IS THE ONE PLACE IN PHASE 14 THAT INTRODUCES INFORMATION PHASE 13 DECLINED TO
# INTRODUCE. Phase 13 published two log Bayes factors against a shared `random` reference and
# deliberately stopped there (D-08): a Bayes factor is evidence, and turning evidence into a
# posterior requires a prior over the three hypotheses, which Phase 13 refused to invent. Phase 14
# needs that posterior -- the Bayesian-FDR rule is a sort over `P(H0 | Z)` and nothing else -- so
# it supplies the prior. THE PRIOR IS MEASURED, NEVER CHOSEN: it is `h.meta.pi_class_masses`, the
# realized class fractions of 200,000 prior draws under the inherited three-way cut, read off the
# net's own metadata. That distinction belongs in the report in exactly those words.
#
# THE PERMUTATION HAZARD, AND WHY EVERY READ BELOW IS BY NAME.
# FIVE class orderings coexist in this codebase: `three_way_probs` returns (coloc, exclusion);
# `three_way_log_bf` returns (coloc, exclusion, random); `ThreeWayLogBF` stores (coloc, random,
# exclusion); the Phase-13 gate's confusion order is ["exclusion","random","coloc"]; and
# `pi_class_masses` is keyed (exclusion, random, coloc). ECE, AUC and FDR are all invariant to a
# CONSISTENT relabelling, which is exactly what makes an inconsistent one so dangerous -- it
# produces numbers that look right. This milestone has already paid for that once: an orthogonal
# permutation preserves norms and spectra, so the Phase-12 wrong-basis bug was O(1) and yet
# invisible to structural verification, and it silently inverted a phase's headline finding.
# THERE IS NO POSITIONAL INDEX ANYWHERE IN THIS FILE. Not one `[1]`, `[2]` or `[3]`; the
# log-sum-exp maximum is a three-term `max` over named locals rather than a `maximum` over a
# collection, precisely so no reordering can change a number.
#
# LOG-SUM-EXP IS NOT AN OPTIMIZATION HERE. The raw log Bayes factors this net emits span roughly
# +-22 nats (the Phase-13 gate report's continuity block records max_abs_diff = 22.02), and the
# decision layer must stay finite on the tails beyond that. `exp` of a raw value overflows to Inf
# or underflows to 0, and Inf/Inf is NaN -- a silent NaN posterior would propagate into the FDR
# sort as an unordered element. The maximum is subtracted BEFORE any `exp`.
#
# DECOUPLING (hard constraint, CLAUDE.md and D-01): spike-local. Touches no src/ byte, adds no
# dependency, writes nothing, opens no image and consumes no reported stream. It reaches
# spike/p13/ READ-ONLY through guarded includes, and it obtains the three-way cut ONLY through
# `p14_load_tau()` (D-07) -- no threshold literal appears here under any name.
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three. Stated here
# because the same reminder is carried in three separate spike/p13/ files and the mistake is
# silent.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :p14_class_posterior) || include(joinpath(@__DIR__, "posterior.jl"))

# --- Guarded includes, in dependency order ------------------------------------------------------
# The Phase-14 Tier-1 pre-registration: P14_CLASS_KEYS, P14_FIX_SEED, P14_FIXTURE_COUNTER.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
# The Phase-13 evidence net: `load_three_way`, so a caller can reach `meta.pi_class_masses` from
# the same namespace it composes the posterior in.
isdefined(@__MODULE__, :ThreeWayEvidenceNet) || include(joinpath(@__DIR__, "..", "p13", "net.jl"))
# The Phase-13 result surface: P13_LOG_BF_KEYS, the frozen key set of the evidence triple.
isdefined(@__MODULE__, :ThreeHypothesisColocResult) ||
    include(joinpath(@__DIR__, "..", "p13", "result.jl"))
# The Phase-13 label surface: `prior_class_draws` and `class_masses`, REUSED rather than
# re-implemented -- a local re-implementation is how two versions of one statistic diverge.
isdefined(@__MODULE__, :three_way_label) || include(joinpath(@__DIR__, "..", "p13", "labels.jl"))
# The D-07 threshold loader. The ONLY sanctioned route to the inherited three-way cut.
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED) || include(joinpath(@__DIR__, "provenance.jl"))

# THE GUARD COVERS ONLY THE `const`s. Julia 1.12 DROPS every docstring written inside an
# `if ... end` block (the parser emits the `Core.@doc` call but the docsystem never registers it),
# so the documented functions sit at top level; method redefinition under a re-include is silent
# and harmless, while `const` redefinition is what would warn. Same split, same reason, as
# spike/p13/result.jl and spike/p14/provenance.jl.
if !isdefined(@__MODULE__, :P14_POSTERIOR_LOADED)
    # The canonical Float64 three-class posterior type, in the frozen P14_CLASS_KEYS order. The
    # Phase-14 analogue of `ThreeWayLogBF` (spike/p13/result.jl:115), and it exists for the same
    # reason: a NamedTuple removes the positional hazard outright.
    const P14ClassPosterior = NamedTuple{P14_CLASS_KEYS, Tuple{Float64, Float64, Float64}}

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else
    # can make this body skip.
    const P14_POSTERIOR_LOADED = true
end

# --- Validation ---------------------------------------------------------------------------------

"""
    _p14_check_class_keys(nt, expected, what) -> nothing

Assert that `nt` carries EXACTLY the keys `expected`, in any order, and throw an `ArgumentError`
naming the caller and the observed keys otherwise.

Missing or extra keys are rejected rather than silently dropped, which is the same contract
`_p13_as_three_way_logbf` enforces on the evidence triple (`spike/p13/result.jl:179-184`). Key
ORDER is deliberately not constrained: every read in this file is by name, so accepting any order
is safe, and asserting the key SET is what actually catches a caller who built the tuple from a
differently-ordered upstream surface.
"""
function _p14_check_class_keys(nt::NamedTuple, expected::NTuple{3, Symbol}, what::AbstractString)
    Set(keys(nt)) == Set(expected) || throw(ArgumentError(
        "$what: expected exactly the keys $expected (in any order); got $(keys(nt))."))
    return nothing
end

# --- The posterior ------------------------------------------------------------------------------

"""
    p14_class_posterior(logbf, prior) -> NamedTuple{P14_CLASS_KEYS}

The three-class posterior `P(class | Z)` over `{coloc, random, exclusion}`, composed from the two
Phase-13 log Bayes factors and a class prior.

`logbf` must carry exactly the keys `P13_LOG_BF_KEYS` -- `(:coloc, :random, :exclusion)`, in any
order -- where `coloc` is `log BF(coloc : random)`, `exclusion` is `log BF(exclusion : random)`,
and `random` is the reference class scored against itself and therefore EXACTLY `0.0` (D-08). The
zero is checked with `===` and not `==`, because `-0.0 == 0.0` is true while the two are different
bit patterns, and the claim is that the reference was scored against itself.

`prior` must carry exactly the keys `P14_CLASS_KEYS`, each strictly positive.

# The prior is MEASURED, never chosen

The value to pass is `h.meta.pi_class_masses` from `load_three_way` -- the realized class fractions
of 200,000 prior draws under the inherited three-way cut, `(exclusion = 0.31272, random = 0.47073,
coloc = 0.21655)`. **This is the only step in Phase 14 that introduces information Phase 13 declined
to introduce**, and the report must say so in those words. `p14_rederive_class_masses` re-derives
those masses from the simulator prior and the frozen label rule, so the number is corroborated
rather than taken on trust (Assumptions Log A7).

The prior is an explicit ARGUMENT, not a hidden default, because the guarantee downstream of it is
only as good as the prior's match to the batch: `pi_class_masses` is the SIMULATOR prior's class
mix, not any real batch's prevalence, and a batch that is 90 % coloc makes every null posterior too
large (the rule too conservative) while a batch that is 1 % coloc makes them too small (the FDR
claim false). That sensitivity is reported over the frozen `P14_PI_COLOC_GRID` rather than hidden
behind one number.

# The arithmetic

Both log Bayes factors are likelihood ratios against the SAME `random` reference (that is what
D-08's structural zero buys), so a prior over the three hypotheses composes exactly:

    a_k  = log(prior_k) + logbf_k          (logbf_random being the structural zero)
    p_k  = exp(a_k - max a) / sum(exp(a - max a))

The maximum is subtracted BEFORE any `exp`. This is not an optimization: raw log Bayes factors span
roughly +-22 nats on this net, so `exp` on a raw value under/overflows and `Inf/Inf` is `NaN`. The
returned masses are asserted to sum to 1 within 1e-12 -- cheap, and it catches every arithmetic
slip, including the one that matters most (pooling a composite null by summing Bayes factors or
logits rather than probabilities, which is not exact and is not what happens here).

Every class is read and written BY NAME. There is no positional index in this function, and the
three-term `max` is by name for the same reason.

# Returns

`(coloc = ..., random = ..., exclusion = ...)` in the frozen `P14_CLASS_KEYS` order, summing to 1.
"""
function p14_class_posterior(logbf::NamedTuple, prior::NamedTuple)
    # VALIDATE FIRST, then compute. A key-set failure reported from inside the arithmetic would
    # surface as an obscure `type NamedTuple has no field` rather than as the caller's mistake.
    _p14_check_class_keys(logbf, P13_LOG_BF_KEYS, "p14_class_posterior: `logbf`")
    _p14_check_class_keys(prior, P14_CLASS_KEYS, "p14_class_posterior: `prior`")

    # D-08's STRUCTURAL ZERO, with `===` semantics. A non-zero reference entry means the triple was
    # built against some other reference (or a third measurement was fabricated), and every number
    # downstream of it would silently be on a different scale.
    logbf.random === 0.0 || throw(ArgumentError(
        "p14_class_posterior (D-08): logbf.random must be exactly 0.0 -- the random class is the " *
        "REFERENCE, scored against itself. Got $(logbf.random)."))

    # A non-positive prior mass has no logarithm. Left unchecked it becomes a silent -Inf that
    # zeroes a whole class, which is a decision, not an arithmetic accident.
    (prior.coloc > 0 && prior.random > 0 && prior.exclusion > 0) || throw(ArgumentError(
        "p14_class_posterior: every prior class mass must be strictly positive; got " *
        "(coloc = $(prior.coloc), random = $(prior.random), exclusion = $(prior.exclusion))."))

    # THE UNNORMALIZED LOG WEIGHTS, ONE NAMED LOCAL PER CLASS. Named locals rather than a tuple
    # because a tuple would invite `a[1]`, and a positional read here is the whole hazard.
    a_coloc     = log(prior.coloc)     + logbf.coloc
    a_random    = log(prior.random)                      # + logbf.random, the D-08 structural zero
    a_exclusion = log(prior.exclusion) + logbf.exclusion

    # LOG-SUM-EXP. Subtract the maximum BEFORE `exp`: the raw log BFs span roughly +-22 nats on
    # this net, so `exp` on a raw value overflows to Inf (or underflows to 0) and Inf/Inf is NaN.
    # A three-term `max` by name, never `maximum` over a positional collection.
    m = max(a_coloc, a_random, a_exclusion)
    w_coloc     = exp(a_coloc     - m)
    w_random    = exp(a_random    - m)
    w_exclusion = exp(a_exclusion - m)
    s = w_coloc + w_random + w_exclusion

    # WRITTEN BY NAME, exactly as it was read. The literal keyword form is the frozen
    # P14_CLASS_KEYS layout, so this IS a `P14ClassPosterior` without a positional construction.
    p = (coloc = w_coloc / s, random = w_random / s, exclusion = w_exclusion / s)
    @assert abs(p.coloc + p.random + p.exclusion - 1.0) < 1e-12 "p14_class_posterior: the three-class posterior does not sum to 1 (got $(p.coloc + p.random + p.exclusion)); the composite null is only additive over a genuine partition"
    return p
end

"""
    p14_null_posterior(p) -> Float64

The scalar null posterior `v = P(H0 | Z) = 1 - p.coloc` that the Bayesian-FDR rule sorts on.

# The composite null is H0 = {random OR exclusion}, and that choice is deliberate

The decision this phase emits is `{coloc / not / ABSTAIN}`, so "not coloc" IS the complement of
coloc. Four reasons this pooling is right rather than merely convenient:

 1. **It needs no approximation.** Posterior probability is additive over a partition, so
    `P({random OR exclusion} | Z) = P(random | Z) + P(exclusion | Z)` is exact. The composite-null
    objection that bites *frequentist* FDR -- at which null value do you compute the p-value? --
    simply does not arise for a posterior probability.
 2. **A null of {random} alone would leave exclusion items UNCONTROLLED.** An exclusion item
    wrongly called coloc would then count as neither a true nor a false discovery. That is a hole
    in the guarantee, not a tighter guarantee.
 3. **It keeps the FDR unit identical to the decision unit** (D-04: per image pair, one call). A
    three-way FDR would need three levels and a multiple-comparison story across hypotheses, which
    nothing in this phase's scope supports.
 4. **Calling a segregated item a colocalization discovery is the worst available error**, and it
    is precisely the error the composite null prices in.

Use [`p14_null_split`](@ref) to carry the decomposition beside the total, because a batch whose
false discoveries are mostly exclusion mass is scientifically different from one whose are mostly
random mass, and the report must be able to say which.
"""
p14_null_posterior(p::NamedTuple) = 1.0 - p.coloc

"""
    p14_null_split(p) -> NamedTuple

The composite null of [`p14_null_posterior`](@ref), decomposed:
`(p_random = p.random, p_exclusion = p.exclusion, v = 1 - p.coloc)`.

`v` is a SUM OF TWO TERMS and the two are not scientifically interchangeable. A batch whose false
discoveries sit mostly on `p_exclusion` is one where the tool is confusing segregation with
colocalization; a batch whose sit mostly on `p_random` is one where it is confusing noise with
signal. Reporting only the total hides which. Every read is by name.
"""
p14_null_split(p::NamedTuple) =
    (p_random = p.random, p_exclusion = p.exclusion, v = p14_null_posterior(p))

"""
    p14_confidence(p) -> Float64

The maximum class posterior, `max(p.coloc, p.random, p.exclusion)`, taken as a THREE-TERM `max`
over named fields -- never `maximum` over a positional collection.

This is the confidence score the risk-coverage curve sweeps and the same quantity the conformal
hedge thresholds, so the two layers describe ONE ordering of the batch rather than two that happen
to agree on the fixtures. Recorded here rather than re-derived at each call site for exactly that
reason.
"""
p14_confidence(p::NamedTuple) = max(p.coloc, p.random, p.exclusion)

# --- The prior, re-derived rather than trusted ---------------------------------------------------

"""
    p14_rederive_class_masses(; n, rng = p14_fix_rng(P14_FIXTURE_COUNTER)) -> NamedTuple{P14_CLASS_KEYS}

Re-derive the class prior from the simulator prior and the frozen three-way label rule, and return
it re-keyed BY NAME into the `P14_CLASS_KEYS` order.

# Why this exists

The whole Bayesian-FDR guarantee rests on `h.meta.pi_class_masses` being the class mix of the prior
the net was trained under. Assumptions Log A7 records that claim as corroborated only by a table in
a docstring, and recommends an executable re-derivation as cheap confirmation. This is it. The
recommendation was cheap to follow because the three-way label depends ONLY on theta -- two
independent prior draws of rho -- so nothing is simulated, and 20,000 items cost a fraction of a
second.

# What it reuses, and what it refuses to re-implement

`prior_class_draws` (`spike/p13/labels.jl:278`) draws the pairs and applies `three_way_label`;
`class_masses` (`:301`) counts them. Neither is re-implemented here -- a local copy of a labelling
rule is how two versions of one cut silently diverge, and the lane guard greps for exactly that.
The three-way cut itself is obtained ONLY through `p14_load_tau()`, the D-07 route: no threshold
literal appears in this file under any name.

# The stream

`rng` defaults to the FIXTURE stream, never a reported one, so a re-derivation can never pre-observe
draws that a reported number will later ride. `n` is a REQUIRED keyword rather than a defaulted one:
the agreement band is `3 * sqrt(p * (1 - p) / n)`, so a caller who does not state `n` has not stated
what the comparison means, and no sample size for this activity is frozen in the pre-registration.

# Returns

`(coloc = ..., random = ..., exclusion = ...)`. Note that `class_masses` returns its fractions in
the (exclusion, random, coloc) order -- one of the five orderings in this codebase -- and the
re-keying below is by NAME, which is the point.
"""
function p14_rederive_class_masses(; n::Integer,
                                     rng::AbstractRNG = p14_fix_rng(P14_FIXTURE_COUNTER))
    n >= 1 || throw(ArgumentError("p14_rederive_class_masses: n must be >= 1 (got $n)"))
    tau = p14_load_tau().tau      # D-07: inherited by LOADING, never a literal and never re-derived
    classes = prior_class_draws(n; rng = rng, tau = tau)
    m = class_masses(classes)
    return (coloc = m.coloc, random = m.random, exclusion = m.exclusion)
end
