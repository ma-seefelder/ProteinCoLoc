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

# spike/p14/fdr.jl --- the direct Bayesian-FDR prefix rule, and the wrapper that refuses to
# report it naked.
#
# NO PACKAGE WAS ADDED, AND NONE WAS NEEDED. The rule is `sortperm` + `cumsum` + `findlast`.
# Bayesian FDR is a sort over posterior probabilities, not a frequentist p-value adjustment, so a
# multiple-testing package would be the wrong tool as well as a new dependency in an environment
# a running test asserts byte-frozen. The file shape -- docstring naming the returned tuple and
# the convention, an explicit degenerate branch BEFORE the arithmetic, pure Base and Statistics --
# is the house hand-rolled numeric-helper shape of `roc_auc` (spike/validation/ood.jl:309-335).
#
# THE ONE SUBSTITUTION THIS FILE IS MOST LIKELY TO SUFFER IS A RUNNING MAXIMUM IN PLACE OF THE
# RUNNING MEAN. The two rules agree on most vectors, which is precisely what makes the swap
# survivable in casual testing: it is the more conservative rule, so it never produces an obviously
# wrong answer -- it silently rejects fewer items and controls a DIFFERENT quantity (the
# per-comparison posterior error probability instead of the posterior expected false discovery
# proportion). `spike/test/test_p14_fdr.jl` carries a fixture chosen because it SEPARATES them.
#
# THE UNIT IS THE IMAGE PAIR (D-04). Per-region/per-tile FDR is out of scope and is not merely
# unimplemented: Phase 12 returned NO on calibrated per-region uncertainty, and `LocalColocMap`
# carries no uncertainty field at all, so there is nothing to control against. A per-tile map may
# be displayed beside a decision; it is never FDR-controlled, and claiming otherwise would be the
# "counting a closed negative as a success" failure this project has already corrected once.
#
# DECOUPLING (hard constraint, CLAUDE.md and D-01): spike-local. Touches no src/ byte, adds no
# dependency, writes nothing, reads no image and consumes no random stream at all.
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :p14_bayes_fdr) || include(joinpath(@__DIR__, "fdr.jl"))

# --- Guarded includes ---------------------------------------------------------------------------
# The Phase-14 Tier-1 pre-registration (P14_ALPHA_FDR_GRID is the frozen SWEEP; the level a call
# runs at is a USER PARAMETER and is never read from a constant here -- D-07).
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
# The composite null this rule sorts on: `v = p14_null_posterior(p) = 1 - P(coloc | Z)`, the
# pooled {random OR exclusion} mass.
isdefined(@__MODULE__, :p14_class_posterior) || include(joinpath(@__DIR__, "posterior.jl"))

# --- Validation ---------------------------------------------------------------------------------

"""
    _p14_check_null_posteriors(v, what) -> nothing

Assert that every entry of `v` is a finite probability, and throw an `ArgumentError` naming the
caller, the offending index and the offending value otherwise.

NOT DECORATIVE. `v` is a vector of posterior probabilities by construction
(`p14_null_posterior`), so a violation means a caller built it some other way. The concrete damage
is silent: a `NaN` sorts to the END of the ascending order and turns every running mean after it
into a `NaN`, so the sequence the prefix property is asserted on is no longer ordered and the
`issorted` assertion below would fire with a message about sorting rather than about the caller's
input.
"""
function _p14_check_null_posteriors(v::AbstractVector{<:Real}, what::AbstractString)
    for (i, x) in enumerate(v)
        (isfinite(x) && 0 <= x <= 1) || throw(ArgumentError(
            "$what: every entry of `v` must be a finite posterior probability in [0,1]; " *
            "entry $i is $x."))
    end
    return nothing
end

# --- The rule -----------------------------------------------------------------------------------

"""
    p14_bayes_fdr(v, alpha) -> NamedTuple

The direct posterior-probability Bayesian-FDR rule (Newton, Noueiry, Sarkar & Ahlquist 2004;
Müller, Parmigiani, Robert & Rousseau 2004).

`v[i] = P(H0 | Z_i)` is the null posterior of item `i` -- for this phase, the composite null
`{random OR exclusion}` produced by [`p14_null_posterior`](@ref). Sort ascending, take the running
MEAN, and accept the LARGEST PREFIX whose running mean stays at or below `alpha`:

    FDP_hat(k) = (1/k) * sum of the k smallest v
    k*         = max { k : FDP_hat(k) <= alpha }      (empty => k* = 0, reject nothing)

# The running MEAN, and why a running maximum is a different rule

The controlled quantity is `E[FDP | data] = (1/|R|) * sum over R of P(H0_i | data)`, the posterior
EXPECTED false discovery proportion -- a mean by construction. A running maximum would control the
PER-COMPARISON posterior error probability instead: a different, and much more conservative, rule.
The two agree on most vectors, so the substitution is survivable in casual testing; the fixture in
`spike/test/test_p14_fdr.jl` is chosen because it separates them.

# What is controlled, stated with its qualifier

`E[FDP | data] <= alpha`, **conditional on the model being right** -- calibrated posteriors and a
class prior that matches the batch. This is the BAYESIAN, posterior-expected FDR. It is **not** a
frequentist long-run FDR guarantee, and the phrase "controls the FDR at alpha" must never be
written without that qualifier. (Storey's pFDR is the frequentist quantity that coincides with it
under a two-groups model; that connection is worth one sentence in a report and no more.)

Note what is NOT needed: `E[sum of 1{H0_i}] = sum of v_i` follows from linearity of expectation
alone, so the guarantee does **not** assume independence across the batch. That is the standard
reviewer question and the answer is clean.

# The prefix property

The running mean of an ASCENDING sequence is non-decreasing, so `k*` really is a prefix and the
rule reduces exactly to a threshold `t*` on `v` (accept iff `v_i <= t*`). That is what makes it a
RULE rather than a subset search, and it is `@assert`ed below rather than assumed.

# The reported cost ratio

`cost_ratio = t*/(1 - t*)` is the decision-risk cost ratio the accepted threshold implies: under a
bivariate loss `c * (false discoveries) + (false negatives)`, the optimal Bayes rule is exactly a
threshold on the null posterior at `c/(1 + c)` (Müller, Parmigiani & Rice 2007). It is REPORTED so
the decision-theoretic reading is available to a referee, and the implementation does not commit to
it: SC1 asks for a user-set FDR level, not a cost ratio, and a three-hypothesis loss matrix is a
quantity this project has no basis to elicit and would therefore have to invent.

# Returns

`(accepted, k_star, fdp, t_star, cost_ratio, n)` where `accepted` holds ORIGINAL indices into `v`
(not sorted positions), `fdp` is the realized posterior expected false discovery proportion of the
accepted set, and `t_star` / `cost_ratio` are `NaN` when nothing is accepted.
"""
function p14_bayes_fdr(v::AbstractVector{<:Real}, alpha::Real)
    # VALIDATE FIRST (the house order): the function name, the requirement, then the observed
    # value.
    0.0 < alpha < 1.0 || throw(ArgumentError(
        "p14_bayes_fdr: alpha must lie strictly inside (0,1); got $alpha."))
    _p14_check_null_posteriors(v, "p14_bayes_fdr")

    # THE DEGENERATE BRANCH, BEFORE THE ARITHMETIC. An empty batch is a real answer -- it is what a
    # batch that abstained on everything hands to this rule -- and it must not be reached through
    # an empty `cumsum` and an empty `findlast`.
    n = length(v)
    n == 0 && return (accepted = Int[], k_star = 0, fdp = 0.0,
                      t_star = NaN, cost_ratio = NaN, n = 0)

    ord = sortperm(v)                       # ascending; ORIGINAL indices, kept to the end
    sv  = v[ord]
    run = cumsum(sv) ./ (1:n)               # the running MEAN -- see the docstring
    @assert issorted(run) "p14_bayes_fdr: the running mean of an ascending vector is not non-decreasing, so the accepted set would not be a prefix. This is asserted rather than assumed: without it, `findlast` would be selecting a k whose predecessors are not all admissible, and the rule would silently be a subset search."

    k = findlast(<=(alpha), run)
    k === nothing && return (accepted = Int[], k_star = 0, fdp = 0.0,
                             t_star = NaN, cost_ratio = NaN, n = n)

    t = sv[k]
    return (accepted = ord[1:k], k_star = k, fdp = run[k],
            t_star = t, cost_ratio = t / (1 - t), n = n)
end

"""
    p14_fdr_over_decided(v_decided, alpha; n_total) -> NamedTuple

[`p14_bayes_fdr`](@ref) over the DECIDED subset of a batch, merged with the scope fields that make
the resulting number a claim rather than a fragment.

`v_decided` holds the null posteriors of the decided items ONLY -- the batch minus abstentions.
`n_total` is the size of the batch they were drawn from and is a REQUIRED keyword.

# Why `n_total` is required rather than optional

`(alpha, n_decided/n_total)` is ONE quantity, not two: a rule that abstains on 95 % of a batch
hits any alpha trivially, so quoting either number alone is misleading. Making the keyword
required is the structural version of that sentence -- the same trick as `p13_calibration_meta`,
where `auc` is required precisely so a calibration verdict can never ship without its
discrimination number (`spike/p13/result.jl:313-315`).

# The honest claim, in the words the report must use

> The posterior expected false discovery proportion is controlled at the user-set level **over the
> DECIDED subset only** -- the batch minus abstentions. It is **not** controlled over the whole
> batch: an abstained item is neither a discovery nor a non-discovery and appears in neither the
> numerator nor the denominator. The guarantee is conditional on the model (calibrated posteriors
> and a correct class prior) and is a posterior expectation, not a frequentist long-run rate.

`fdr_scope = :decided_subset_only` carries that sentence in machine-readable form, so an artifact
cannot be read as claiming batch-wide control.

# ABSTAIN FIRST, THEN SORT -- and why the reverse order is wrong (Assumption A2)

The rejection set `R` is a deterministic function of the observed data, hence sigma(data)-measurable,
so it pulls straight out of the conditional expectation:

    E[FDP | data] = (1/|R|) * sum over R of P(H0_i | data)

**for ANY data-dependent selection of `R`** -- including one that filtered on OOD status, conformal
ambiguity, or anything else. Unlike frequentist FDR, where a data-dependent pre-selection is a
genuine selective-inference problem needing correction, the Bayesian posterior quantity is immune
because the conditioning has already happened.

The REVERSE order -- sort the whole batch, then abstain from some of the accepted items -- is
wrong, and is recorded here so it is not re-derived: removing items changes both the numerator and
the denominator of the running mean in an uncontrolled way, so the residual accepted set may
EXCEED the level. Abstaining first also improves the premise, since it routes away precisely the
items whose posteriors are least likely to be right.

This derivation is **Assumption A2** in `14-RESEARCH.md`. It is LOAD-BEARING -- it is the whole
justification for abstain-then-sort and for the decided-subset scope -- and it is a two-line
derivation rather than a quoted theorem, so the report must state it explicitly and let a referee
check it rather than assert it.

# Returns

The `p14_bayes_fdr` tuple merged with
`(n_total, n_decided, decided_fraction, fdr_scope = :decided_subset_only)`.
"""
function p14_fdr_over_decided(v_decided::AbstractVector{<:Real}, alpha::Real; n_total::Integer)
    n_total > 0 || throw(ArgumentError(
        "p14_fdr_over_decided: n_total must be > 0; got $n_total. A decided fraction over an " *
        "empty batch is not a number, and the fraction is not optional."))
    length(v_decided) <= n_total || throw(ArgumentError(
        "p14_fdr_over_decided: the decided subset ($(length(v_decided))) is larger than the " *
        "batch it came from (n_total = $n_total); that is a bookkeeping error, not a decided " *
        "fraction above 1."))

    inner = p14_bayes_fdr(v_decided, alpha)
    # SCOPE MERGED OVER THE RULE, never under it: a caller cannot overwrite the scope label that
    # says what the number does and does not cover.
    return merge(inner, (n_total          = Int(n_total),
                         n_decided        = length(v_decided),
                         decided_fraction = length(v_decided) / n_total,
                         fdr_scope        = :decided_subset_only))
end
