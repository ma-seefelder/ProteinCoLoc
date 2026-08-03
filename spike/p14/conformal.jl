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

# spike/p14/conformal.jl --- the distribution-free hedge: split conformal, HAND-ROLLED (D-02).
#
# SC1 IS AMENDED BY D-02, AND THE AMENDMENT IS CITED HERE BECAUSE A READER OF THIS FILE IS EXACTLY
# THE READER WHO NEEDS IT. The ROADMAP's SC1 names a library by name. That name is DROPPED and the
# criterion is met IN SUBSTANCE, hand-rolled, for two reasons that were both verified rather than
# assumed:
#   1. Adding it would force a `Pkg.resolve()` against an environment a RUNNING TEST asserts
#      byte-frozen (spike/test/test_p14_decoupling.jl, "the spike environment is byte-frozen"), and
#      this milestone has already paid for that once -- the Phase-7 Wave-0 co-resolution gate FAILED
#      on a transitive conflict that pushed NeuralEstimators back below the pin the whole spike is
#      built on.
#   2. It wraps MLJ models, not a NeuralEstimators posterior, so it would have to be adapted anyway
#      -- while the construction it would supply is `sort` + `ceil` + integer indexing.
# ANY REPORT QUOTING A NUMBER FROM THIS FILE MUST CITE THE ORIGINAL SC1 WORDING ALONGSIDE THE
# AMENDMENT. "Met in substance" is a claim a reader is entitled to check, not one to be told.
#
# THE ONE FUNCTION THIS FILE MUST NOT CALL IS `Statistics.quantile` (Pitfall 6). The finite-sample
# conformal guarantee is stated for an ORDER STATISTIC. The default type-7 rule INTERPOLATES
# LINEARLY between two neighbouring order statistics and so returns a strictly smaller value for
# any non-degenerate sample; a smaller qhat admits fewer classes, and the >= 1 - alpha coverage
# bound breaks by O(1/n). It breaks SILENTLY: at n = 2000 the shift is invisible beside the binomial
# noise of the evaluation set. `spike/test/test_p14_conformal.jl` carries the directed n = 9 test
# that separates the two, and the plan's acceptance grep keeps the interpolating call out of this
# source. Note that the name IS in scope here regardless -- spike/p13/result.jl does `using
# Statistics` and the spike lane is one flat namespace -- so the guard is the test and the grep,
# not the import. `Statistics` is imported QUALIFIED below for the median in the score summary,
# which is the honest minimum.
#
# THE SCORE IS LAC / THR: s(Z, y) = 1 - p_hat(y | Z). The rejected alternative is recorded because
# the rejection is load-bearing, not stylistic. A BF-MARGIN score is a monotone function of the log
# Bayes-factor magnitudes, so its (1 - alpha) quantile sits BY CONSTRUCTION out in the confident
# tail -- and 13-REPORT.md hands Phase 14 a named finding that the corrected log Bayes factors are
# unbiased in the bulk but lose sub-nat precision exactly there, "where the sign of the verdict is
# not in doubt but the magnitude is least determined". A threshold placed there would rest on the
# least-determined part of the estimate. The posterior probability compresses that same tail into
# [0,1], so its quantile sits where the probability is CALIBRATED (both Phase-13 heads are
# ECE-green with zero empty bins on the very probabilities this score reads). APS is the
# pre-registered fallback and the ONE authorised iteration is reserved for it
# (P14_ITERATION_TRIGGER); it is not chosen here because it buys conditional coverage this net does
# not need at a set-size cost D-05 explicitly does not want.
#
# DECOUPLING (hard constraint, CLAUDE.md and D-01): spike-local. Touches no src/ byte, ADDS NO
# DEPENDENCY, writes nothing, reads no image, opens no artifact and consumes no reported stream.
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :p14_conformal_quantile) || include(joinpath(@__DIR__, "conformal.jl"))

# --- Guarded includes, in dependency order ------------------------------------------------------
# The Phase-14 Tier-1 pre-registration: P14_CLASS_KEYS and P14_AMBIGUOUS_RATE_FLOOR. NOTE that
# P14_ALPHA_CONFORMAL is NOT read anywhere in this file: alpha is an ARGUMENT here, so a caller
# cannot silently inherit a level it did not state, and alpha_FDR can never be derived from it
# (D-07 -- the two share the value 0.10 at one grid point, which is exactly why).
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
# The posterior this score is computed ON -- the SAME object the Bayesian-FDR rule sorts on, so the
# hedge and the primary rule cannot contradict each other. Also supplies `_p14_check_class_keys`,
# reused rather than re-implemented.
isdefined(@__MODULE__, :p14_class_posterior) || include(joinpath(@__DIR__, "posterior.jl"))

import Statistics       # QUALIFIED: only `Statistics.median`, for the score summary. See the banner
                        # for why the interpolating quantile is nonetheless reachable in this flat
                        # namespace and why the guard is the directed test, not the import.

# THE GUARD COVERS ONLY THE `const`s. Julia 1.12 DROPS every docstring written inside an
# `if ... end` block, so the documented functions sit at top level; method redefinition under a
# re-include is silent and harmless, while `const` redefinition is what would warn. Same split,
# same reason, as spike/p13/result.jl and spike/p14/posterior.jl.
#
# THE STATUS VOCABULARY IS DELIBERATELY DUPLICATED IN spike/p14/fuse.jl, AND THE DUPLICATION IS
# SELF-POLICING. `fuse.jl` must stay loadable WITHOUT this file, because this file pulls the whole
# Phase-13 stack (Flux, NeuralEstimators, Images) through posterior.jl and the fusion rule is pure
# symbol arithmetic that a caller must be able to test in a second, not in twenty. So both files
# declare the same three symbols under the same name, each behind the same isdefined-guard, and
# whichever loads SECOND asserts agreement rather than silently losing. A textual divergence
# therefore fails loudly the moment both files are in one process -- which `decide.jl` guarantees.
if isdefined(@__MODULE__, :P14_CONFORMAL_STATUSES)
    @assert P14_CONFORMAL_STATUSES === (:singleton, :ambiguous, :empty) "spike/p14/conformal.jl: P14_CONFORMAL_STATUSES was already defined as $(P14_CONFORMAL_STATUSES), which disagrees with this file's (:singleton, :ambiguous, :empty). The conformal status vocabulary is duplicated in spike/p14/fuse.jl on purpose (see the note above); a divergence must fail loudly here rather than resolve to whichever file loaded first."
else
    const P14_CONFORMAL_STATUSES = (:singleton, :ambiguous, :empty)
end

if !isdefined(@__MODULE__, :P14_CONFORMAL_LOADED)
    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else can
    # make this body skip.
    const P14_CONFORMAL_LOADED = true
end

# --- Validation ---------------------------------------------------------------------------------

"""
    _p14_check_scores(scores, what) -> nothing

Assert that every calibration score is finite, and throw an `ArgumentError` naming the caller, the
offending index and the offending value otherwise.

NOT DECORATIVE, and the damage is silent in a specific way. `sort` places `NaN` at the END of the
ascending order, so a single `NaN` anywhere in a calibration set makes the top order statistic
`NaN`, and `NaN <= qhat` is `false` for every class — the conformal set comes back EMPTY for every
input, and an empty set is a signal this phase reports as a misspecification finding. A silent
arithmetic fault would then be indistinguishable from the most interesting result the hedge can
produce.
"""
function _p14_check_scores(scores::AbstractVector{<:Real}, what::AbstractString)
    for (i, x) in enumerate(scores)
        isfinite(x) || throw(ArgumentError(
            "$what: every calibration score must be finite; entry $i is $x."))
    end
    return nothing
end

# --- The score ------------------------------------------------------------------------------------

"""
    p14_lac_score(p, y) -> Float64

The LAC / THR nonconformity score `s(Z, y) = 1 - p_hat(y | Z)`. HIGHER means a WORSE fit, which is
the convention every conformal construction assumes and the same "higher ⇒ worse" convention
`roc_auc` states for the OOD detectors.

`p` is the three-class posterior from [`p14_class_posterior`](@ref) and `y` a member of
`P14_CLASS_KEYS`. The class is read BY NAME via `getproperty` — never by position — because five
class orderings coexist in this codebase and every headline metric is invariant to a *consistent*
relabelling, which is exactly what makes an inconsistent one produce numbers that look right.

# Why this score and not a Bayes-factor margin

A margin score is a monotone function of the log Bayes-factor magnitudes, so its `(1 - alpha)`
quantile sits by construction in the confident tail — and Phase 13 hands this phase a named finding
that the corrected log Bayes factors are unbiased in the bulk and lose sub-nat precision exactly
there. The posterior probability compresses that tail into `[0,1]`, so its quantile sits where the
probability is calibrated. Both Phase-13 heads are ECE-green with zero empty bins on precisely the
probabilities this score reads.
"""
function p14_lac_score(p::NamedTuple, y::Symbol)
    y in P14_CLASS_KEYS || throw(ArgumentError(
        "p14_lac_score: `y` must be one of the frozen class keys $(P14_CLASS_KEYS); got $y."))
    _p14_check_class_keys(p, P14_CLASS_KEYS, "p14_lac_score: `p`")
    return 1.0 - getproperty(p, y)
end

# --- The quantile ---------------------------------------------------------------------------------
#
# WHY THE FEASIBILITY BOUND IS COMPUTED FROM THE ARGUMENT `n` AND NEVER WRITTEN AS A LITERAL.
# D-03a WITHDREW the figure "alpha >= 1/31 ~ 0.032" on two counts, both verified directly against
# corpus/manifest.csv: it was computed from the wrong n (the corpus is partitioned dev = 14 /
# eval = 16 / sealed = 2, not "30 open"), AND it was computed on rows that are not physical ground
# truth (30 of 32 are tier: simulated-secondary -- the CBS benchmark, i.e. another simulator). A
# withdrawn figure is not a superseded one: it must never become an executable default, because a
# default is how a corpus-derived floor would silently calibrate the hedge that D-03a just removed
# it from. spike/test/test_p14_decoupling.jl bans the literal and exempts the general expression.
#
# THIS NOTE LIVES IN A COMMENT AND NOT IN THE DOCSTRING BELOW, DELIBERATELY. The same guard file
# strips whole-line comments before scanning and then bans the corpus tokens outright in every
# non-allowlisted Phase-14 source -- so prose about the holdout is free HERE and a docstring
# mentioning it is not. That is Pitfall 8 ("the token ban is a code ban, not a prose ban") landing
# on this file: the first draft of the docstring below named the token, and this guard caught it.

"""
    p14_conformal_quantile(scores, alpha) -> Float64

The FINITE-SAMPLE conformal quantile: the `ceil(Int, (n + 1) * (1 - alpha))`-th ORDER STATISTIC of
the calibration scores.

# DO NOT USE `Statistics.quantile`

Its default type-7 rule interpolates LINEARLY between two neighbouring order statistics and returns
a strictly smaller value for any non-degenerate sample. A smaller quantile admits fewer classes into
the set, so the finite-sample `>= 1 - alpha` coverage bound breaks by `O(1/n)`. **It breaks
silently**: at `n = 2000` the shift is entirely invisible beside the binomial noise of the
evaluation set, and no aggregate number would look wrong. The conformal guarantee is stated for the
ORDER STATISTIC and for nothing else. `spike/test/test_p14_conformal.jl` carries the directed
`n = 9` test that separates the two answers, and it asserts BOTH equalities because either one
flipping is the failure.

# It THROWS rather than clamps, and that is the deliberate difference from `ghat`

`ghat` (`spike/simulator/ghat.jl:110-121`) is the closest in-repo analog — a frozen, order-statistic-
style lookup — and it CLAMPS out-of-range input to the nearest swept endpoint, which is right for a
calibration curve. Here it would be wrong: `ceil((n + 1)(1 - alpha)) > n` means the requested level
is not achievable from `n` points AT ALL, and clamping to the maximum would return a number that
looks like a quantile while carrying no guarantee whatsoever. Split conformal requires
`alpha >= 1/(n + 1)`, and the message computes that bound FROM THE ARGUMENT `n` — never from a
hard-coded figure. D-03a WITHDREW one such hard-coded floor (see the note above this function);
`test_p14_decoupling.jl` bans the withdrawn literal outright, and the general expression in `n` is
the exemption it carves out.

# Dependencies

`sort`, `ceil`, integer indexing. That is the whole of it (D-02, and 14-RESEARCH §B.5).
"""
function p14_conformal_quantile(scores::AbstractVector{<:Real}, alpha::Real)
    # VALIDATE FIRST (the house order): the function name, the requirement, then the observed value.
    0.0 < alpha < 1.0 || throw(ArgumentError(
        "p14_conformal_quantile: alpha must lie strictly inside (0,1); got $alpha."))
    n = length(scores)
    n >= 1 || throw(ArgumentError(
        "p14_conformal_quantile: the calibration set is empty (n = 0); a conformal quantile over " *
        "no points is not a number, and returning one would be a guarantee with nothing behind it."))
    _p14_check_scores(scores, "p14_conformal_quantile")

    k = ceil(Int, (n + 1) * (1 - alpha))
    k <= n || throw(ArgumentError(
        "p14_conformal_quantile: alpha = $alpha is infeasible at n = $n. Split conformal takes " *
        "the ceil((n + 1) * (1 - alpha))-th order statistic, which here is index $k > n = $n; " *
        "feasibility requires alpha >= 1/(n + 1) = $(1 / (n + 1)) at this n. This THROWS rather " *
        "than clamping: a clamped value would look like a quantile and carry no coverage " *
        "guarantee at all."))

    # `sort`, never `sort!` -- a caller's calibration vector is not this function's to reorder.
    return sort(collect(Float64, scores))[k]
end

# --- The set --------------------------------------------------------------------------------------

"""
    p14_conformal_set(p, qhat) -> NamedTuple

The LAC prediction set `C(Z) = { y : s(Z, y) <= qhat }` and its size regime, returned as
`(set, status, qhat)`. LAC applies ONE global threshold to `1 - p_hat(y)` across all classes, so
equivalently `C(Z) = { y : p_hat(y | Z) >= 1 - qhat }`.

`status` is one of `P14_CONFORMAL_STATUSES` and the three regimes are kept SEPARATE:

| `length(set)` | `status`     | meaning                                                     | D-05 |
|---------------|--------------|-------------------------------------------------------------|------|
| 1             | `:singleton` | exactly one class is plausible at level `alpha`              | DECIDE  |
| >= 2          | `:ambiguous` | the evidence does not separate two or more classes           | ABSTAIN |
| 0             | `:empty`     | *no* class reaches the plausibility floor                    | ABSTAIN |

# Why the EMPTY case has its own symbol and is never folded into `:ambiguous`

Both abstain, so collapsing them would change no decision — and would destroy the most informative
single output of this hedge. An empty set means the observed point is **more nonconforming than
`(1 - alpha)` of everything the calibration set contained**. That is a distribution-free
misspecification signal arriving through a **completely different channel** from the
Mahalanobis / posterior-predictive / noise OOD detector: it is the one signal in this phase that does
not depend on the simulator being right about *densities*, only about *exchangeability*. Counted
separately it can be cross-tabulated against the OOD state — agreement is corroboration worth
reporting, and disagreement is the more interesting finding. Folded in, it is unrecoverable.

Note the arithmetic, because in many textbook settings this case is unreachable and here it is not:
with three classes and `p_hat` summing to 1, an empty set requires `1 - qhat > 1/3`, i.e.
`qhat < 2/3`. Since `qhat` is the `(1 - alpha)` quantile of `1 - p_hat(true class)`, that happens
when the calibration scores concentrate near 0 — a confident, well-separated classifier, which is
what this net is. So the empty regime is REACHABLE for this net and must be measured, not assumed
away.

`qhat` is echoed back on the result so a set can never be read without the threshold it was cut at.
"""
function p14_conformal_set(p::NamedTuple, qhat::Real)
    _p14_check_class_keys(p, P14_CLASS_KEYS, "p14_conformal_set: `p`")
    isfinite(qhat) || throw(ArgumentError(
        "p14_conformal_set: qhat must be finite; got $qhat. A non-finite threshold silently " *
        "admits every class or none, and both look like a legitimate regime."))

    # Built by iterating the FROZEN key tuple, so the emitted set is always in P14_CLASS_KEYS order
    # regardless of the order the caller's posterior happened to be declared in.
    set = Tuple(y for y in P14_CLASS_KEYS if p14_lac_score(p, y) <= qhat)
    status = length(set) == 1 ? :singleton : (isempty(set) ? :empty : :ambiguous)
    return (set = set, status = status, qhat = Float64(qhat))
end

"""
    _p14_set_status(s) -> Symbol

The status of one hedge outcome, accepting either a full [`p14_conformal_set`](@ref) tuple or a bare
status symbol, and validating the result against `P14_CONFORMAL_STATUSES`.

Two shapes, ONE vocabulary. The validation is the point: a mistyped symbol would otherwise be
counted into none of the three rates, so the rates would still sum to less than 1 and nothing would
say why.
"""
function _p14_set_status(s)
    st = s isa Symbol ? s : s.status
    st in P14_CONFORMAL_STATUSES || throw(ArgumentError(
        "_p14_set_status: a hedge outcome must carry one of $(P14_CONFORMAL_STATUSES); got $st."))
    return st
end

# --- The diagnostics ------------------------------------------------------------------------------

"""
    p14_hedge_diagnostics(sets, qhat) -> NamedTuple

The realized set-size distribution of a batch of conformal sets, returned as
`(qhat, n, singleton_rate, ambiguous_rate, empty_rate, vacuous_hedge)`.

# `vacuous_hedge` is a LABEL, not a verdict

If `qhat` comes out very small — a highly separable calibration set — then `length(C(Z)) == 1`
almost always and the hedge never fires. **That looks like a pass and measures nothing.** This is
the same shape as the `vacuous_pass` guard Phase 13 built for ECE
(`spike/p13/result.jl:284-298`: *"a head that learns nothing produces well-calibrated-looking
probabilities"*), and it gets the same treatment: `vacuous_hedge = true` when
`ambiguous_rate + empty_rate` falls below the pre-registered `P14_AMBIGUOUS_RATE_FLOOR`, reported
BESIDE the coverage number and travelling with every quoted figure.

It is emphatically **not** a failure condition, and that is deliberate. Abstaining rarely can also be
the CORRECT behaviour on well-specified data, so making it a gate would penalise the right answer.
What must not happen is a coverage number quoted without it.

# The empty batch is refused, not answered

Zero sets would give `0 + 0 < floor`, i.e. `vacuous_hedge = true` — the right word for entirely the
wrong reason: a hedge that was never exercised is not a hedge that never fires. The degenerate case
is therefore thrown before the arithmetic rather than returned through it.
"""
function p14_hedge_diagnostics(sets, qhat::Real)
    n = length(sets)
    n >= 1 || throw(ArgumentError(
        "p14_hedge_diagnostics: the outcome collection is empty. A vacuity LABEL over zero sets " *
        "would report `vacuous_hedge = true` for a hedge that was never exercised, which is the " *
        "right word for the wrong reason."))
    isfinite(qhat) || throw(ArgumentError(
        "p14_hedge_diagnostics: qhat must be finite; got $qhat."))

    st = [_p14_set_status(s) for s in sets]
    singleton_rate = count(==(:singleton), st) / n
    ambiguous_rate = count(==(:ambiguous), st) / n
    empty_rate     = count(==(:empty), st) / n
    return (qhat           = Float64(qhat),
            n              = n,
            singleton_rate = singleton_rate,
            ambiguous_rate = ambiguous_rate,
            empty_rate     = empty_rate,
            vacuous_hedge  = (ambiguous_rate + empty_rate) < P14_AMBIGUOUS_RATE_FLOOR)
end

# --- Calibration ----------------------------------------------------------------------------------

"""
    p14_conformal_calibrate(posteriors, true_classes, alpha) -> NamedTuple

Fit the split-conformal threshold from labelled calibration pairs and return
`(qhat, n_cal, alpha_conformal, scores_summary)`.

The whole construction is: score each pair with [`p14_lac_score`](@ref), sort, take the
`ceil((n + 1)(1 - alpha))`-th order statistic. **`ConformalPrediction.jl` is NOT a dependency** —
SC1's named mechanism is AMENDED by D-02 and met in SUBSTANCE, and this docstring is one of the
three places (source, artifact, report) the amendment is required to be cited. A report quoting any
number produced here must carry the original SC1 wording alongside the amendment; "met in substance"
is a claim a reader is entitled to check.

`alpha_conformal` is echoed back on the result so a threshold can never be read without the level it
was cut at. It is an ARGUMENT of this function, never read from a constant here and never derived
from `alpha_FDR` — the two happen to share the value 0.10 at one grid point of the pre-registered
sweep, which is precisely why neither may ever be assigned from the other (D-07).

`scores_summary` is `(min, median, max)`. It is reported because a `qhat` quoted alone cannot be
sanity-checked: a near-zero `qhat` beside near-zero scores is a separable classifier, while a
near-zero `qhat` beside a wide score spread is a bug.

# The exchangeability precondition this function CANNOT check

The guarantee requires the calibration points to be exchangeable with the test point. A
CLASS-BALANCED calibration set drawn against a prior-distributed deployment stream is NOT
exchangeable with it, and the realized coverage will then miss `1 - alpha` by an amount driven by the
class-conditional score distributions — silently, because the code runs and the number looks fine.
Nothing in this signature can detect that: it sees scores, not how they were drawn. The calibration
draws must be UNSTRATIFIED, and the runner that produces them is where that is asserted.
"""
function p14_conformal_calibrate(posteriors, true_classes, alpha::Real)
    length(posteriors) == length(true_classes) || throw(ArgumentError(
        "p14_conformal_calibrate: `posteriors` and `true_classes` must be the same length; got " *
        "$(length(posteriors)) and $(length(true_classes)). A zip over mismatched lengths would " *
        "silently truncate the calibration set and shift the quantile."))
    n = length(posteriors)
    n >= 1 || throw(ArgumentError(
        "p14_conformal_calibrate: the calibration set is empty (n = 0)."))

    s = [p14_lac_score(p, y) for (p, y) in zip(posteriors, true_classes)]
    qhat = p14_conformal_quantile(s, alpha)
    ss = sort(s)
    return (qhat            = qhat,
            n_cal           = n,
            alpha_conformal = Float64(alpha),
            scores_summary  = (min    = first(ss),
                               median = Statistics.median(ss),
                               max    = last(ss)))
end
