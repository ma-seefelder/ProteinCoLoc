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

# spike/p14/result.jl --- the types the decision layer emits: `P14Result` (per image pair) and
# `P14BatchDecision` (per batch).
#
# THIS FILE IS WHERE THIS PHASE'S HONESTY COMMITMENTS STOP BEING PROSE AND BECOME FIELDS. Every
# one of them was, at some point in the planning documents, a sentence someone could forget to
# copy into a report. Here each is a name a compiler and a test can see:
#
#     fdr_scope = :decided_subset_only   the claim is controlled over the DECIDED SUBSET ONLY
#     local_map_uncontrolled             the per-tile map is display, never a controlled call (D-04)
#     cross_method                       the disagreement is RECORDED, not suppressed (D-05)
#     ood_state                          the OOD input is THREE-valued, never a Bool (D-06)
#     abstain_reason                     every silence names its trigger (D-05/D-06)
#     decided_fraction                   the denominator alpha is meaningless without (RESEARCH E.2)
#
# D-01, SAID PLAINLY AND REPEATED IN EVERY REPORT THAT QUOTES THIS TYPE. `P14Result` has a
# `src/`-shaped signature and subtypes the SHIPPED `AbstractColocResult`, but "SRC-SHAPED" IS NOT
# "IN SRC": the decision layer is built in the spike research lane, it is NOT shipped in v2.0, and
# this file makes ZERO edits to src/. The three reaches into src/ below are read-only `include`s,
# exactly as spike/p13/result.jl:36-45 established for the same reason.
#
# COMPOSITION, NOT RESTATEMENT (14-RESEARCH section D.2, recommendation (a)). The Phase-13 evidence
# triple is carried as a FIELD of type `ThreeHypothesisColocResult` rather than re-stored here. Two
# consequences, both wanted: D-08's `random === 0.0` structural zero is enforced ONCE, at Phase
# 13's own inner constructor, so there is no second place a differently-referenced triple can enter
# the system; and the Phase-14 type is strictly ADDITIVE over Phase 13's, so a reader can see
# exactly what the decision layer added and nothing else.
#
# WHAT `p14_named_limits()` IS FOR. The manuscript's limit list is enumerable FROM CODE. A report
# that drops one is then a report that failed to iterate a tuple, which is visible, rather than a
# report that forgot a sentence, which is not.
#
# A NOTE ON THE WORDING OF THOSE LIMITS, RECORDED SO IT IS NOT READ AS EVASION. This file is
# scanned by spike/test/test_p14_decoupling.jl, whose corpus-token ban is applied to the
# COMMENT-STRIPPED source -- a docstring and a string literal are CODE to that scan, and 14-05
# already paid for that once. `spike/p14/result.jl` is deliberately NOT added to the guard's
# allowlist (adding an entry there is a decision, not a fix), so the named limits say "real-data"
# where the planning documents say the substrate's name. Nothing is softened: the substrate D-03a
# withdrew is the one whose 30 of 32 rows are `tier: simulated-secondary` (another simulator) and
# whose only two physical-primary rows are the sealed reservation for the Phase-16 blind
# evaluation, with zero bytes on disk. The full record lives in 14-CONTEXT.md D-03a and in
# P14_DECLARED_DEVIATIONS (spike/p14/consts.jl section H), which IS allowlisted and names it.
#
# WHAT THIS FILE DELIBERATELY DOES NOT DO. It never computes anything. It validates, it stores and
# it reports. No threshold literal appears here under any name, no random stream is consumed, and
# nothing is written to disk. The decision itself is `spike/p14/fuse.jl`'s and the arithmetic is
# `spike/p14/fdr.jl`'s and `posterior.jl`'s; a result type that recomputed either would be a second
# implementation of a rule that already has one.
#
# DECOUPLING (hard constraint, CLAUDE.md and D-01): spike-local. Reaches src/ READ-ONLY through
# guarded includes, touches no src/ byte, adds no dependency, writes nothing and opens no image.
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :P14Result) || include(joinpath(@__DIR__, "result.jl"))

# ORDER MATTERS (spike/contract.jl:34-39, the documented load order): StatsBase + Statistics must
# be in scope before src/colocalization.jl is reached transitively, and Images before
# src/LoadImages.jl (its convenience constructor calls Images.otsu_threshold). Both are reached
# through spike/p13/result.jl -> spike/validation/sbc.jl -> harness.jl -> contract.jl below, so the
# `using` lines come FIRST and the includes after.
using StatsBase      # corspearman/corkendall for the transitively-reached correlation() Dict
using Statistics     # mean, quantile
using Images         # otsu_threshold (transitive src/LoadImages.jl requirement)

# --- Guarded includes, in dependency order ------------------------------------------------------
# READ-ONLY reach into src/ for the supertype, `_iface_error`, `OODVerdict` and `CalibrationMeta`.
isdefined(@__MODULE__, :AbstractColocResult) ||
    include(joinpath(@__DIR__, "..", "..", "src", "results.jl"))
# The Phase-14 Tier-1 pre-registration: P14_CLASS_KEYS and P14_PI_COLOC_GRID. No BAR and no
# threshold is read from it here -- this file has no number in it to tune.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
# The Phase-13 result type, COMPOSED as a field (14-RESEARCH D.2 (a)). Pulls src/results.jl,
# spike/p13/consts.jl and spike/validation/sbc.jl transitively.
isdefined(@__MODULE__, :ThreeHypothesisColocResult) ||
    include(joinpath(@__DIR__, "..", "p13", "result.jl"))
# The fusion vocabulary: P14_OOD_STATES, P14_ABSTAIN_REASONS, P14_CONFORMAL_STATUSES and
# `p14_abstain_reasons()`. Pure symbol arithmetic; it loads in a fraction of a second.
isdefined(@__MODULE__, :p14_fuse) || include(joinpath(@__DIR__, "fuse.jl"))
# READ-ONLY reach into src/ for `LocalColocMap`. `src/registry.jl` comes first because
# local_map.jl annotates a method with `::EstimatorBundle` and a type annotation is EVALUATED at
# definition time -- without it the include raises UndefVarError rather than defining anything.
# `MultiChannelImage` (its other definition-time requirement) is already in scope by here, through
# spike/contract.jl on the Phase-13 chain above; that is why this pair is LAST.
isdefined(@__MODULE__, :EstimatorBundle) ||
    include(joinpath(@__DIR__, "..", "..", "src", "registry.jl"))
isdefined(@__MODULE__, :LocalColocMap) ||
    include(joinpath(@__DIR__, "..", "..", "src", "amortized", "local_map.jl"))

# THE GUARD BLOCK COVERS ONLY THE `const`s AND THE `struct`s, NOT THE FUNCTIONS, AND THAT SPLIT IS
# DELIBERATE. Julia 1.12 DROPS every docstring written inside an `if ... end` block -- the parser
# emits the `Core.@doc` call but the docsystem never registers it, so `@doc f` returns `nothing`
# for a function documented inside a guard (verified on 1.12.6). The two lossy accessors below
# (`bayes_factor` and `is_ood`) are REQUIRED to be documented and the accompanying test asserts the
# docstring is retrievable from the running system, so every function lives at TOP LEVEL where the
# docsystem can see it. Method redefinition under a re-include is silent and harmless in Julia;
# only `const` and `struct` redefinition would warn or throw, and those are what the guard protects.
if !isdefined(@__MODULE__, :P14Result)

# The CLOSED set of decisions. `{coloc / not / ABSTAIN}` is the phase boundary's own vocabulary and
# a fourth state is not reachable: the four-state design adding an explicit DECIDED-CONTRA-CLASSICAL
# action was considered and rejected for v2.0 (14-CONTEXT Deferred Ideas), and that information is
# carried instead as `cross_method.disagree_any` beside a two-state decision.
const P14_DECISIONS = (:coloc, :not_coloc, :abstain)

# WHERE THE CLASS PRIOR CAME FROM, as a closed vocabulary rather than free text. The DEFAULT is the
# MEASURED `h.meta.pi_class_masses` -- the realized class fractions of the prior draws the net was
# trained under -- so the default is measured, not chosen. Empirical-Bayes estimation of pi from
# the batch is deliberately NOT offered as a third member: it is a second estimated quantity with
# its own failure modes, and estimating the prior from the same batch the guarantee is quoted over
# would make that guarantee circular (14-RESEARCH E.1, "Prevalence mismatch").
const P14_CLASS_PRIOR_SOURCES = (:measured_pi_class_masses, :caller_supplied)

# The ONLY legal scope of a Phase-14 FDR claim, as a symbol rather than a sentence.
const P14_FDR_SCOPE = :decided_subset_only

# The token every headline carries. UPPERCASE and literal, so a grep over a report finds it and a
# reader skimming one cannot miss it.
const P14_HEADLINE_SCOPE_TOKEN = "FDR SCOPE: DECIDED SUBSET ONLY"

"""
    P14Result

The per-image-pair output of the Phase-14 decision layer: a `{coloc / not_coloc / abstain}` call,
the three-class posterior it was made from, the conformal hedge's verdict, the three-valued OOD
state, the recorded cross-method comparison, and -- explicitly labelled as uncontrolled -- the
per-tile display map.

A NEW subtype of the real `AbstractColocResult` (Phase-7 D-02), defined in the spike lane per D-01.
`src/results.jl` is byte-unchanged and no field is bolted onto any shipped struct. **"src-shaped"
is not "in src": the decision layer is NOT shipped in v2.0.**

# Fields
- `three_way::ThreeHypothesisColocResult`: the Phase-13 evidence result, COMPOSED rather than
  re-stated, so D-08's `random === 0.0` structural zero is enforced once at its own constructor.
- `class_posterior::NamedTuple`: `P(class | Z)` over `P14_CLASS_KEYS`, summing to 1.
- `null_posterior::Float64`: `v = 1 - P(coloc | Z)`, the COMPOSITE null `{random OR exclusion}`
  that the Bayesian-FDR rule sorts on.
- `null_split::NamedTuple`: `(p_random, p_exclusion, v)` -- the composite null decomposed, because
  a batch whose false discoveries sit on exclusion mass is scientifically different from one whose
  sit on random mass.
- `decision::Symbol`: one of `P14_DECISIONS`.
- `abstain_reason::Union{Symbol,Nothing}`: the trigger that fired, from `p14_abstain_reasons()`,
  and `nothing` for a decided item. Present if and only if `decision === :abstain`.
- `ood_state::Symbol`: one of `P14_OOD_STATES` -- `:fired`, `:clear` or `:not_checked` (D-06).
- `conformal::NamedTuple`: `(set, status, qhat)` exactly as `p14_conformal_set` returns it.
- `cross_method::NamedTuple`: the classical comparison, carried as a NAMED FIELD rather than
  suppressed (D-05). Must carry `basis` -- which SCALE the comparison was made on.
- `local_map_uncontrolled::Union{Nothing,LocalColocMap}`: the per-tile map, or `nothing`. **Not an
  FDR-controlled output** -- see [`local_map`](@ref) and D-04.
- `meta::NamedTuple`: run metadata. If -- and ONLY if -- the run actually carried control posterior
  draws, put them in `meta.delta_rho_draws`; see [`delta_rho`](@ref).

The inner constructor rejects, with an `ArgumentError` naming `P14Result`, the decision it enforces
and the observed value: an unknown decision or OOD state, an abstention without a reason, a decided
call carrying one, a class posterior on the wrong keys or not summing to 1, a null posterior that
disagrees with its own class posterior, a conformal tuple on an unknown status, and a cross-method
record that does not say which scale it was computed on.
"""
struct P14Result <: AbstractColocResult
    three_way              :: ThreeHypothesisColocResult
    class_posterior        :: NamedTuple
    null_posterior         :: Float64
    null_split             :: NamedTuple
    decision               :: Symbol
    abstain_reason         :: Union{Symbol, Nothing}
    ood_state              :: Symbol
    conformal              :: NamedTuple
    cross_method           :: NamedTuple
    # THE LABEL IS IN THE TYPE, NOT ONLY IN A DOCSTRING (D-04). A caller reading a field list --
    # which is what a caller reads -- cannot mistake this for a controlled call.
    local_map_uncontrolled :: Union{Nothing, LocalColocMap}
    meta                   :: NamedTuple

    function P14Result(three_way::ThreeHypothesisColocResult, class_posterior::NamedTuple,
                       null_posterior::Real, null_split::NamedTuple, decision::Symbol,
                       abstain_reason::Union{Symbol, Nothing}, ood_state::Symbol,
                       conformal::NamedTuple, cross_method::NamedTuple,
                       local_map_uncontrolled::Union{Nothing, LocalColocMap},
                       meta::NamedTuple)
        decision in P14_DECISIONS || throw(ArgumentError(
            "P14Result: `decision` must be one of $(P14_DECISIONS); got $decision. The phase " *
            "boundary emits three states and a fourth was considered and rejected for v2.0."))

        # D-06. NOT a Bool, and it cannot be coerced to one: `ood_verdict` returns a false flag for
        # EVERY input when no threshold was recorded, so a `false` there means NOT CHECKED rather
        # than in distribution, and storing a Bool would erase that distinction at the boundary.
        ood_state in P14_OOD_STATES || throw(ArgumentError(
            "P14Result (D-06): `ood_state` must be one of $(P14_OOD_STATES); got $ood_state. " *
            "The OOD input is THREE-valued -- a `false` flag on an unchecked channel means NOT " *
            "CHECKED, not in distribution."))

        # D-05/D-06: EVERY ABSTAIN CARRIES WHICH TRIGGER FIRED. A silent abstention is nearly as
        # bad as a wrong call, because a reader cannot tell a tool that DECLINED from a tool that
        # BROKE, nor an OOD abstention (fix the input) from a conformal-empty one (question the
        # simulator). The converse is enforced too: a decided call carrying a reason would let an
        # abstention trigger be recorded on an item that was never silenced.
        if decision === :abstain
            (abstain_reason isa Symbol && abstain_reason in p14_abstain_reasons()) ||
                throw(ArgumentError(
                    "P14Result (D-05/D-06): an :abstain decision must carry an `abstain_reason` " *
                    "from $(p14_abstain_reasons()); got $(repr(abstain_reason)). Every ABSTAIN " *
                    "carries which trigger fired -- a silent abstention is nearly as bad as a " *
                    "wrong call."))
        else
            abstain_reason === nothing || throw(ArgumentError(
                "P14Result (D-05): a $decision decision must carry `abstain_reason === nothing`; " *
                "got $(repr(abstain_reason)). A reason on a decided item records a silence that " *
                "never happened."))
        end

        # THE CLASS-ORDER CONTRACT. Five class orderings coexist in this codebase and every one of
        # ECE, AUC and FDR is invariant to a CONSISTENT relabelling -- which is exactly what makes
        # an inconsistent one so dangerous: it produces numbers that look right. The key SET is
        # asserted; the order is not, because every read here is BY NAME.
        Set(keys(class_posterior)) == Set(P14_CLASS_KEYS) || throw(ArgumentError(
            "P14Result: `class_posterior` must carry exactly the keys $(P14_CLASS_KEYS) (in any " *
            "order); got $(keys(class_posterior))."))
        s = class_posterior.coloc + class_posterior.random + class_posterior.exclusion
        abs(s - 1.0) < 1e-12 || throw(ArgumentError(
            "P14Result: `class_posterior` must sum to 1 within 1e-12; got $s. The composite null " *
            "is additive only over a genuine partition, so a posterior that does not sum to 1 " *
            "makes every null posterior downstream of it wrong."))

        # The composite null and the class posterior are ONE quantity stored twice; a disagreement
        # means they came from different items or different arithmetic.
        v = Float64(null_posterior)
        abs(v - (1.0 - class_posterior.coloc)) < 1e-12 || throw(ArgumentError(
            "P14Result: `null_posterior` must equal 1 - class_posterior.coloc within 1e-12; got " *
            "$v against $(1.0 - class_posterior.coloc). The composite null is {random OR " *
            "exclusion}, i.e. the exact complement of coloc."))

        # The hedge's own vocabulary, single-sourced with spike/p14/conformal.jl. `:empty` is a
        # status of its own and is never reported as ambiguous.
        (haskey(conformal, :set) && haskey(conformal, :status) && haskey(conformal, :qhat)) ||
            throw(ArgumentError(
                "P14Result: `conformal` must carry (:set, :status, :qhat), exactly as " *
                "`p14_conformal_set` returns them; got $(keys(conformal)). A set read without " *
                "the threshold it was cut at is not interpretable."))
        conformal.status in P14_CONFORMAL_STATUSES || throw(ArgumentError(
            "P14Result: `conformal.status` must be one of $(P14_CONFORMAL_STATUSES); got " *
            "$(conformal.status)."))

        # PITFALL 5, MADE STRUCTURAL. `basis` records WHICH SCALE the classical comparison was made
        # on, and that is precisely the thing that is otherwise forgotten -- this project has four
        # recorded cases of a criterion applied in a unit it was not derived for. A cross-method
        # record without its basis is a comparison whose scale nobody can reconstruct later.
        haskey(cross_method, :basis) || throw(ArgumentError(
            "P14Result (D-05, Pitfall 5): `cross_method` must carry the key `:basis` -- which " *
            "SCALE the classical comparison was made on; got $(keys(cross_method))."))

        return new(three_way, class_posterior, v, null_split, decision, abstain_reason,
                   ood_state, conformal, cross_method, local_map_uncontrolled, meta)
    end
end

end # if !isdefined(@__MODULE__, :P14Result) -- consts + structs only

# --- The four REQUIRED accessors (the AbstractColocResult interface) -----------------------------

"""
    posterior_draws(r::P14Result)

The posterior draws backing the result (parameter x draw layout), passed straight through from the
composed [`ThreeHypothesisColocResult`](@ref).
"""
posterior_draws(r::P14Result) = posterior_draws(r.three_way)

"""
    bayes_factor(r::P14Result) -> Float64

Return `log BF(coloc : random)` from the composed three-hypothesis result.

DOCUMENTED RATHER THAN SILENT, for the same reason as
`bayes_factor(::ThreeHypothesisColocResult)`: a single-number accessor on a three-way result is
inherently lossy and the reader deserves to know which number they got. The shared interface's
`bayes_factor` has exactly one shipped meaning -- "the colocalization Bayes factor evidence"
(`src/results.jl:59-65`) -- and the coloc-vs-random entry is the semantically closest quantity this
result carries, so a consumer written against the binary `AmortizedColocResult` keeps working and
gets the number it meant.

A caller who wants THE THREE-WAY EVIDENCE must use [`log_bf_vs_random`](@ref), which returns all
three entries by name. This accessor deliberately says nothing about the exclusion hypothesis, and
it says nothing at all about the DECISION: the call is [`decision`](@ref), the evidence is this.
"""
bayes_factor(r::P14Result) = bayes_factor(r.three_way)

"""
    is_ood(r::P14Result) -> Bool

`true` exactly when `ood_state(r) === :fired`.

**LOSSY ON PURPOSE, AND LOSSY IN A DIRECTION THAT MATTERS.** A `false` return COLLAPSES `:clear`
and `:not_checked`, which are different states (D-06): the first means a detector answered and said
the input looks in-distribution; the second means no detector ran at all, because `ood_verdict`
returns a false flag for every input when no operating point was recorded. Reading the second as
the first makes the tool speak confidently in exactly the regime where it has no evidence it
should.

**The decision layer NEVER routes through this accessor.** It reads [`ood_state`](@ref), and
`p14_fuse` takes the three-valued symbol. This method exists so that a consumer written against the
shipped `AbstractColocResult` interface keeps working -- and its docstring exists so that such a
consumer knows what it lost. A caller who cares about the difference must call `ood_state(r)`.

Same treatment, same reason, as `bayes_factor(::ThreeHypothesisColocResult)`.
"""
is_ood(r::P14Result) = r.ood_state === :fired

"""
    delta_rho(r::P14Result)

Return the carried Delta rho draws if -- and only if -- the run actually produced them, i.e. if
`meta.delta_rho_draws` is present. Otherwise fall through to `_iface_error`.

**DO NOT FABRICATE A DELTA RHO.** Delta rho is a sample-minus-control contrast that requires a
CONTROL posterior; the three-way decision path does not need one and in general does not carry one.
Synthesising a plausible-looking value (zeros, the coloc draws, a difference against a notional
zero control) would put a number in front of a reader that no measurement backs. Letting the
documented interface error fire is the honest answer, and it is exactly what `_iface_error`
(`src/results.jl:47-49`) exists for: it names the accessor, the type, and the four accessors the
interface requires.

D-05's own counter-example is why this matters: a sample at rho = 0.8 under a control at rho = 0.9
has Delta rho < 0, which a naive three-way reading would publish as "mutually exclusive" while the
sample is strongly colocalized. The three-way call and the sign of a contrast are not the same
question, and a fabricated Delta rho would make them look like one.
"""
function delta_rho(r::P14Result)
    haskey(r.meta, :delta_rho_draws) || _iface_error(r, :delta_rho)
    return r.meta.delta_rho_draws
end

"""
    log_bf_vs_random(r::P14Result) -> ThreeWayLogBF

The full evidence triple by name, passed through from the composed three-hypothesis result:
`(coloc = log BF(C:R), random = 0.0, exclusion = log BF(E:R))`.

Defined on `P14Result` so that [`bayes_factor`](@ref)'s docstring is ACTIONABLE rather than merely
true: it tells a caller who wants the three-way evidence to use this accessor, and a caller holding
a `P14Result` must therefore be able to.
"""
log_bf_vs_random(r::P14Result) = log_bf_vs_random(r.three_way)

# --- The Phase-14 accessors, LAYERED ON TOP and never bolted-on fields ---------------------------

"""
    three_way(r::P14Result) -> ThreeHypothesisColocResult

The composed Phase-13 evidence result. The Phase-14 type is strictly ADDITIVE over Phase 13's, and
this accessor is how a consumer reaches the part that is not new.
"""
three_way(r::P14Result) = r.three_way

"""
    decision(r::P14Result) -> Symbol

The call: one of `P14_DECISIONS` -- `:coloc`, `:not_coloc` or `:abstain`.

An `:abstain` is not a failure and not a missing value; it is an ANSWER, and it always carries its
trigger in [`abstain_reason`](@ref).
"""
decision(r::P14Result) = r.decision

"""
    abstain_reason(r::P14Result) -> Union{Symbol,Nothing}

Which trigger fired, from the closed set `p14_abstain_reasons()`, or `nothing` for a decided item.

The constructor enforces the biconditional -- a reason is present if and only if the decision is
`:abstain` -- so a caller may branch on either one and never on both. `:ood_not_checked` is a
DISTINCT reason from `:ood_fired` (D-06) and `:conformal_empty` a distinct one from
`:conformal_ambiguous` (D-05); collapsing either pair would destroy the only information that
separates "fix the input" from "question the simulator".
"""
abstain_reason(r::P14Result) = r.abstain_reason

"""
    ood_state(r::P14Result) -> Symbol

The THREE-valued OOD state: `:fired`, `:clear` or `:not_checked` (D-06).

**This is the accessor the decision layer uses, and the one a consumer who cares must use.**
[`is_ood`](@ref) cannot round-trip three states through a `Bool` and collapses the last two.
"""
ood_state(r::P14Result) = r.ood_state

"""
    conformal_set(r::P14Result) -> Tuple

The classes the split-conformal hedge admits at `alpha_conformal`, in the frozen `P14_CLASS_KEYS`
order. Read it beside [`conformal_status`](@ref) and `r.conformal.qhat`: a set can never be
interpreted without the threshold it was cut at.
"""
conformal_set(r::P14Result) = r.conformal.set

"""
    conformal_status(r::P14Result) -> Symbol

One of `P14_CONFORMAL_STATUSES` -- `:singleton`, `:ambiguous` or `:empty`.

`:empty` is kept as its own status and never folded into `:ambiguous`: both abstain, so folding
them would change no decision and would destroy the most informative single output of the hedge --
an empty set means the point is more nonconforming than `(1 - alpha)` of the whole calibration set,
a distribution-free misspecification signal on a completely different channel from the OOD detector.
"""
conformal_status(r::P14Result) = r.conformal.status

"""
    class_posterior(r::P14Result) -> NamedTuple

`P(class | Z)` over `P14_CLASS_KEYS`, summing to 1 (asserted at construction). Every read is BY
NAME; there is no positional index on this type anywhere.
"""
class_posterior(r::P14Result) = r.class_posterior

"""
    null_posterior(r::P14Result) -> Float64

`v = P(H0 | Z) = 1 - P(coloc | Z)`, the COMPOSITE null `{random OR exclusion}` the Bayesian-FDR
rule sorts on.

The pooling is exact rather than approximate -- posterior probability is additive over a partition
-- and it is the right null rather than merely the convenient one: a null of `{random}` alone would
leave exclusion items UNCONTROLLED, so an exclusion item wrongly called coloc would count as
neither a true nor a false discovery. That is a hole in the guarantee, not a tighter guarantee.
"""
null_posterior(r::P14Result) = r.null_posterior

"""
    null_split(r::P14Result) -> NamedTuple

The composite null decomposed: `(p_random, p_exclusion, v)`.

`v` is a SUM OF TWO TERMS and the two are not scientifically interchangeable. A batch whose false
discoveries sit mostly on `p_exclusion` is one where the tool is confusing segregation with
colocalization; a batch whose sit mostly on `p_random` is one where it is confusing noise with
signal. Reporting only the total hides which.
"""
null_split(r::P14Result) = r.null_split

"""
    cross_method(r::P14Result) -> NamedTuple

The classical comparison, RECORDED AS A NAMED FIELD rather than suppressed (D-05).

**Disagreeing with Costes/Manders is the product thesis, not a failure mode.** SC2 as originally
written makes cross-method disagreement a trigger for SILENCE, which would mute the tool exactly in
the cases that demonstrate its headline differentiator -- Phase 9 exists so that v2.0 can be
positioned as *"knows when the classics are wrong, not merely agrees with them"*. D-05 AMENDS SC2:
disagreement ALONE decides, and is recorded here and as the `:decided_contra_classical` reason.
Disagreement TOGETHER WITH a non-clear OOD state abstains, because then there is an independent
reason to distrust ourselves.

Always carries `basis` -- which SCALE the comparison was made on -- and the constructor rejects a
record without it (Pitfall 5).
"""
cross_method(r::P14Result) = r.cross_method

"""
    local_map(r::P14Result) -> Union{Nothing,LocalColocMap}

The per-tile map, carried for display.

**THE PER-TILE MAP IS UNCONTROLLED DISPLAY, NOT A CONTROLLED CALL (D-04).** FDR is controlled PER
IMAGE PAIR, because that is the only unit where the calibration provably holds: the posterior, the
three-way Bayes factor and the SBC proof are all per-pair quantities. `LocalColocMap` carries
`delta_rho` and `ood_flag` but **no uncertainty field** (`src/amortized/local_map.jl:73-79`), so
there is nothing per-tile to control AGAINST -- an FDR statement needs a per-item null posterior and
none exists at tile resolution.

Phase 12 returned NO on calibrated per-region uncertainty (`12-VERIFICATION.md`: goal NOT achieved;
the shipped deliverable is the neutralized ablation), so claiming per-region FDR after that would be
exactly the "counting a closed negative as a success" failure corrected in ROADMAP at `e29c846`.
The field is therefore NAMED `local_map_uncontrolled`, and
[`p14_local_map_is_controlled`](@ref) is the executable, greppable form of the same statement.

Per-region FDR is deferred to v2.1 alongside SPAT-07, and it is deferred because its foundation
does not exist -- not because it was out of budget.
"""
local_map(r::P14Result) = r.local_map_uncontrolled

"""
    p14_local_map_is_controlled(::P14Result) -> Bool

Always `false`. The EXECUTABLE form of D-04, so the claim is greppable and assertable rather than
only readable: the per-tile map is uncontrolled display and no Phase-14 FDR arithmetic may consume
it.

It takes the result rather than no argument so that it reads as a property OF A RESULT at every
call site, which is where the question is actually asked.
"""
p14_local_map_is_controlled(::P14Result) = false
