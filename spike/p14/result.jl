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

# --- The batch object ---------------------------------------------------------------------------

# The keys every rung of the prior-sensitivity curve must carry. Enumerated so a curve computed
# with a missing column fails at construction rather than at report time.
const P14_PI_SENSITIVITY_KEYS = (:pi_coloc, :k_star, :fdp, :n_accepted)

"""
    P14BatchDecision

The batch-level output of the decision layer. Bayesian FDR is a BATCH quantity -- the rule is a
sort over null posteriors and a prefix cut -- so the per-pair results alone cannot express it.

**THE FIELD THIS TYPE EXISTS FOR IS `fdr_scope`.** It is always `:decided_subset_only`, and that is
not cosmetic: it is the machine-readable form of the phase's honest claim.

> The posterior expected false discovery proportion is controlled at the user-set level **over the
> DECIDED subset only** -- the batch minus abstentions. It is **not** controlled over the whole
> batch: an abstained item is neither a discovery nor a non-discovery and appears in neither the
> numerator nor the denominator. The guarantee is conditional on the model (calibrated posteriors
> and a correct class prior) and is a posterior expectation, not a frequentist long-run rate.

**AND THE MANDATORY COMPANION NUMBER.** A rule that abstains on 95 % of a batch hits any alpha
trivially, so `(alpha_fdr, n_decided/n_total)` is ONE quantity and quoting either alone is
misleading. `n_total` is therefore a REQUIRED keyword of [`p14_batch_decision`](@ref) rather than
`length(results)`, `decided_fraction` is a stored field, and [`p14_headline`](@ref) emits both in a
single string so a copy-paste cannot separate them.

# Fields
- `results::Vector{P14Result}`: the WHOLE batch, abstentions included.
- `alpha_fdr::Float64`: the user-set FDR level this batch was cut at (SC1: a per-call USER
  parameter, never read from a constant, and never derived from the conformal miscoverage level).
- `k_star::Int`: the accepted prefix length.
- `threshold_t_star::Float64`: the implied threshold on the null posterior; `NaN` when nothing was
  accepted.
- `implied_cost_ratio::Float64`: `t*/(1 - t*)`, REPORTED so the decision-theoretic reading is
  available to a referee; the implementation does not commit to it. `NaN` when nothing was accepted.
- `posterior_expected_fdp::Float64`: the realized `E[FDP | data]` of the accepted set.
- `accepted::Vector{Int}`: **positions INTO THE DECIDED SUBSET**, exactly as `p14_bayes_fdr`
  returns them -- NOT positions into `results`. `results[accepted]` is therefore WRONG whenever
  anything abstained, and the constructor asserts `1 <= i <= n_decided` to catch the confusion.
  A caller that needs the mapping back to batch positions must carry its own decided-index vector
  in `meta`; this type does not invent one.
- `n_total::Int`, `n_decided::Int`, `n_abstained::Int`, `decided_fraction::Float64`: the
  denominator, its two parts and their ratio, all cross-checked against one another.
- `abstain_reason_counts::NamedTuple`: how many items each trigger silenced. Keys are a subset of
  `p14_abstain_reasons()` and the values sum to `n_abstained`.
- `class_prior_used::NamedTuple`: the class prior the null posteriors were composed from, on
  `P14_CLASS_KEYS`.
- `class_prior_source::Symbol`: one of `P14_CLASS_PRIOR_SOURCES`, so a silently switched prior is
  not possible.
- `pi_sensitivity::Vector{<:NamedTuple}`: the REQUIRED prior-sensitivity curve, one rung per
  `P14_PI_COLOC_GRID` value, each carrying `P14_PI_SENSITIVITY_KEYS`.
- `conformal::NamedTuple`: the batch's hedge diagnostics (`qhat` and the realized set-size rates,
  including the `vacuous_hedge` LABEL).
- `fdr_scope::Symbol`: always `:decided_subset_only`; there is no other legal value in this phase.
- `meta::NamedTuple`: run metadata and provenance.

Build it with [`p14_batch_decision`](@ref). The inner constructor validates the invariants that are
properties of the stored values, so they hold no matter which entry point was used.
"""
struct P14BatchDecision
    results                :: Vector{P14Result}
    alpha_fdr              :: Float64
    k_star                 :: Int
    threshold_t_star       :: Float64
    implied_cost_ratio     :: Float64
    posterior_expected_fdp :: Float64
    accepted               :: Vector{Int}
    n_total                :: Int
    n_decided              :: Int
    n_abstained            :: Int
    decided_fraction       :: Float64
    abstain_reason_counts  :: NamedTuple
    class_prior_used       :: NamedTuple
    class_prior_source     :: Symbol
    pi_sensitivity         :: Vector{<:NamedTuple}
    conformal              :: NamedTuple
    fdr_scope              :: Symbol      # ALWAYS :decided_subset_only -- validated below
    meta                   :: NamedTuple

    function P14BatchDecision(results::Vector{P14Result}, alpha_fdr::Real, k_star::Integer,
                              threshold_t_star::Real, implied_cost_ratio::Real,
                              posterior_expected_fdp::Real, accepted::Vector{Int},
                              n_total::Integer, n_decided::Integer, n_abstained::Integer,
                              decided_fraction::Real, abstain_reason_counts::NamedTuple,
                              class_prior_used::NamedTuple, class_prior_source::Symbol,
                              pi_sensitivity::Vector{<:NamedTuple}, conformal::NamedTuple,
                              fdr_scope::Symbol,   # only :decided_subset_only is legal
                              meta::NamedTuple)
        0.0 < alpha_fdr < 1.0 || throw(ArgumentError(
            "P14BatchDecision: `alpha_fdr` must lie strictly inside (0,1); got $alpha_fdr."))

        # THE SCOPE IS NOT A CALLER CHOICE. Compared against the literal so the value a reader
        # greps for is the value the type enforces; `P14_FDR_SCOPE` is a name for the same literal
        # and is pinned to it by the assertion below this struct.
        fdr_scope === :decided_subset_only || throw(ArgumentError(
            "P14BatchDecision: `fdr_scope` must be :decided_subset_only; got $fdr_scope. There " *
            "is no other legal value in this phase -- an abstained item is neither a discovery " *
            "nor a non-discovery and appears in neither the numerator nor the denominator, so " *
            "batch-wide control is not what was computed and must not be what is claimed."))

        # THE DENOMINATOR AND ITS TWO PARTS, CROSS-CHECKED RATHER THAN DERIVED FROM ONE ANOTHER.
        n_total > 0 || throw(ArgumentError(
            "P14BatchDecision: `n_total` must be > 0; got $n_total. A decided fraction over an " *
            "empty batch is not a number."))
        n_decided + n_abstained == n_total || throw(ArgumentError(
            "P14BatchDecision: n_decided ($n_decided) + n_abstained ($n_abstained) must equal " *
            "n_total ($n_total). Every item is decided or abstained; there is no third state " *
            "and no item may be counted twice."))
        length(results) == n_total || throw(ArgumentError(
            "P14BatchDecision: `results` holds $(length(results)) items but n_total is " *
            "$n_total. `results` is the WHOLE batch, abstentions included; a mismatch means the " *
            "quoted denominator is not the batch the results came from."))
        decided_fraction == n_decided / n_total || throw(ArgumentError(
            "P14BatchDecision: `decided_fraction` ($decided_fraction) must equal n_decided / " *
            "n_total ($(n_decided / n_total)). The decided fraction is the companion number the " *
            "alpha is meaningless without, so it may not be an independently supplied figure."))

        # THE ABSTENTION LEDGER MUST BALANCE. A reason set that does not account for exactly the
        # silenced items means some abstention was recorded without its trigger -- which is the
        # failure D-05 and D-06 both exist to prevent, arriving one level up at the batch.
        all(k -> k in p14_abstain_reasons(), keys(abstain_reason_counts)) || throw(ArgumentError(
            "P14BatchDecision: `abstain_reason_counts` may only carry keys from " *
            "$(p14_abstain_reasons()); got $(keys(abstain_reason_counts))."))
        counted = sum(values(abstain_reason_counts); init = 0)
        counted == n_abstained || throw(ArgumentError(
            "P14BatchDecision: `abstain_reason_counts` sums to $counted but n_abstained is " *
            "$n_abstained. Every ABSTAIN carries which trigger fired; a batch whose reasons do " *
            "not add up has silenced items nobody can explain."))

        # THE ACCEPTED SET IS INDEXED INTO THE DECIDED SUBSET, NOT INTO THE BATCH. This catches
        # the one confusion that would silently report the wrong items as discoveries.
        k_star == length(accepted) || throw(ArgumentError(
            "P14BatchDecision: k_star ($k_star) must equal length(accepted) " *
            "($(length(accepted)))."))
        all(i -> 1 <= i <= n_decided, accepted) || throw(ArgumentError(
            "P14BatchDecision: every entry of `accepted` must be a position INTO THE DECIDED " *
            "SUBSET, i.e. in 1:$n_decided; got $(accepted). These are not positions into " *
            "`results`, and indexing the batch with them would report the wrong items."))
        (0.0 <= posterior_expected_fdp <= alpha_fdr + 1e-12) || throw(ArgumentError(
            "P14BatchDecision: the realized posterior expected FDP ($posterior_expected_fdp) " *
            "must lie in [0, alpha_fdr] ($alpha_fdr); a value outside it did not come from the " *
            "prefix rule this batch claims to have applied."))

        # THE PRIOR, AND WHERE IT CAME FROM.
        class_prior_source in P14_CLASS_PRIOR_SOURCES || throw(ArgumentError(
            "P14BatchDecision: `class_prior_source` must be one of $(P14_CLASS_PRIOR_SOURCES); " *
            "got $class_prior_source. The DEFAULT is the MEASURED class mass of the prior the " *
            "net was trained under, so the default is measured rather than chosen."))
        Set(keys(class_prior_used)) == Set(P14_CLASS_KEYS) || throw(ArgumentError(
            "P14BatchDecision: `class_prior_used` must carry exactly the keys $(P14_CLASS_KEYS) " *
            "(in any order); got $(keys(class_prior_used))."))

        # THE SENSITIVITY CURVE IS REQUIRED OUTPUT, NEVER OPTIONAL. Prevalence mismatch is the
        # biggest way this guarantee fails silently: the class prior is the SIMULATOR's class mix,
        # not any real batch's prevalence, and a batch that is 90 % coloc makes the rule too
        # conservative while one that is 1 % coloc makes the claim false. SC1-c is REPORTED, NOT
        # GATED -- no bar is attached to it, and failing to report it IS the failure.
        length(pi_sensitivity) == length(P14_PI_COLOC_GRID) || throw(ArgumentError(
            "P14BatchDecision: `pi_sensitivity` must carry one rung per P14_PI_COLOC_GRID value " *
            "($(length(P14_PI_COLOC_GRID))); got $(length(pi_sensitivity))."))
        Tuple(e.pi_coloc for e in pi_sensitivity) == P14_PI_COLOC_GRID || throw(ArgumentError(
            "P14BatchDecision: `pi_sensitivity` must sweep exactly $(P14_PI_COLOC_GRID), in that " *
            "order; got $(Tuple(e.pi_coloc for e in pi_sensitivity)). The grid is pre-registered " *
            "so the sweep cannot be chosen after the fact."))
        for (i, e) in enumerate(pi_sensitivity)
            all(k -> haskey(e, k), P14_PI_SENSITIVITY_KEYS) || throw(ArgumentError(
                "P14BatchDecision: rung $i of `pi_sensitivity` must carry " *
                "$(P14_PI_SENSITIVITY_KEYS); got $(keys(e))."))
        end

        return new(results, Float64(alpha_fdr), Int(k_star), Float64(threshold_t_star),
                   Float64(implied_cost_ratio), Float64(posterior_expected_fdp), accepted,
                   Int(n_total), Int(n_decided), Int(n_abstained), Float64(decided_fraction),
                   abstain_reason_counts, class_prior_used, class_prior_source,
                   pi_sensitivity, conformal,
                   fdr_scope,   # proven === :decided_subset_only above
                   meta)
    end
end

# The const and the literal the constructor compares against are pinned to each other MECHANICALLY,
# so the two spellings cannot drift, and the headline token is derived from the symbol rather than
# typed a second time.
@assert P14_HEADLINE_SCOPE_TOKEN ==
        "FDR SCOPE: " * uppercase(replace(String(P14_FDR_SCOPE), "_" => " "))

# --- The named limits, enumerable FROM CODE -----------------------------------------------------
# EVERY ONE OF THESE TRAVELS WITH EVERY PHASE-14 RESULT. A report that drops one is then a report
# that failed to iterate a tuple, which is visible in a diff, rather than a report that forgot a
# sentence, which is not. See the banner for why item 2 says "real-data" rather than naming the
# withdrawn substrate: this file is scanned by the lane guard, a string literal is code to that
# scan, and the full record lives in the allowlisted spike/p14/consts.jl section H.
const P14_NAMED_LIMITS = (
    "The decision layer is built in the spike research lane and is NOT shipped in v2.0 (D-01): it has a src-shaped signature and subtypes the shipped AbstractColocResult, but src-shaped is not in-src, and no byte of src/ is edited by this phase.",
    "The conformal guarantee is SIMULATOR-DERIVED (D-02/D-03) and inherits the simulator's misspecification in full; the intended real-data bound is ABSENT, not merely loose, because the substrate it was to be computed on holds no unsealed physical ground truth and its bytes are unfetched (D-03a).",
    "The real-data check is six committed microscopy TIFFs in two conditions (D-03a) -- an ILLUSTRATION, not a coverage claim. It bounds nothing tightly and must never be quoted as though it did.",
    "FDR is controlled over the DECIDED SUBSET ONLY -- the batch minus abstentions -- and the decided fraction must be quoted beside alpha, because a rule that abstains on most of a batch hits any alpha trivially. The guarantee is a posterior expectation conditional on the model, not a frequentist long-run rate.",
    "The per-tile map is UNCONTROLLED DISPLAY, never an FDR-controlled call (D-04): Phase 12 returned NO on calibrated per-region uncertainty and the shipped per-tile map carries no uncertainty field, so there is nothing per-tile to control against.",
    "SC1 is AMENDED by D-02 (ConformalPrediction.jl is not added; conformal is met in substance, hand-rolled) and SC2 is AMENDED by D-05 (the literal OR of three triggers is replaced by an asymmetric rule in which cross-method disagreement alone DECIDES). Any report citing a Phase-14 result must cite the original ROADMAP wording alongside the amendment.",
    "The four SC3 bars and the 0.20 coverage floor are JUDGEMENT CALLS with no derivation (P14_JUDGEMENT_CALL_BARS). They were ratified by the user before any Phase-14 result existed, and a missed bar is REPORTED, not re-tuned.",
    "The underlying net is an epoch-4 checkpoint of an early-overfitting run (13-REPORT limit C), and every Phase-14 number inherits that limit unchanged.",
)

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

# --- The batch constructor ----------------------------------------------------------------------

"""
    p14_batch_decision(results; alpha_fdr, fdr, n_total, abstain_reason_counts,
                       class_prior_used, class_prior_source, pi_sensitivity, conformal,
                       meta = NamedTuple()) -> P14BatchDecision

Assemble the batch object from the whole batch of per-pair results and the output of
[`p14_fdr_over_decided`](@ref).

**EVERY KEYWORD EXCEPT `meta` IS REQUIRED, AND `n_total` MOST OF ALL.** Deriving the denominator
from `length(results)` is forbidden here, and the reason is the whole point of the field:
`(alpha_fdr, n_decided/n_total)` is ONE quantity, not two, because a rule that abstains on 95 % of
a batch hits any alpha trivially. Making the denominator a REQUIRED INPUT is the structural version
of that sentence -- the same trick as `p13_calibration_meta`, where `auc` is a required keyword
precisely so a calibration verdict can never ship without its discrimination number
(`spike/p13/result.jl:313-315`). It is required AND cross-checked against `length(results)`, so it
must be stated and cannot be misstated.

`fdr` must be the tuple `p14_fdr_over_decided` returns -- not the bare `p14_bayes_fdr` one. That is
enforced by requiring its scope and denominator fields, because the scoped rule is the only one
whose output is a CLAIM rather than a fragment: it already carries `n_total`, `n_decided`,
`decided_fraction` and the scope label, and this constructor asserts all four agree with what the
caller stated.

`n_abstained` is NOT a keyword: it is derived by summing `abstain_reason_counts`, and the identity
`n_decided + n_abstained == n_total` is then a real check rather than a tautology. A batch whose
reasons do not add up has silenced items nobody can explain, which is the D-05/D-06 failure
arriving at batch level.

The scope is set BY THIS CONSTRUCTOR and is not a parameter: `:decided_subset_only` is the only
legal value in this phase.
"""
function p14_batch_decision(results::Vector{P14Result};
                            alpha_fdr::Real,
                            fdr::NamedTuple,
                            n_total::Integer,
                            abstain_reason_counts::NamedTuple,
                            class_prior_used::NamedTuple,
                            class_prior_source::Symbol,
                            pi_sensitivity,
                            conformal::NamedTuple,
                            meta::NamedTuple = NamedTuple())
    # THE RULE'S OWN OUTPUT MUST BE THE SCOPED ONE. Accepting the bare `p14_bayes_fdr` tuple here
    # would let a batch-wide number be stored under a decided-subset label.
    for k in (:accepted, :k_star, :fdp, :t_star, :cost_ratio, :n_total, :n_decided,
              :decided_fraction, :fdr_scope)   # the last must read :decided_subset_only
        haskey(fdr, k) || throw(ArgumentError(
            "p14_batch_decision: `fdr` must be the NamedTuple `p14_fdr_over_decided` returns " *
            "(missing key `$k`). The bare `p14_bayes_fdr` output carries no scope and no " *
            "denominator, so storing it here would label a batch-wide number as a " *
            "decided-subset claim."))
    end
    fdr.fdr_scope === :decided_subset_only || throw(ArgumentError(
        "p14_batch_decision: fdr_scope must be :decided_subset_only; got $(fdr.fdr_scope)."))
    fdr.n_total == n_total || throw(ArgumentError(
        "p14_batch_decision: the rule was run over a batch of $(fdr.n_total) but `n_total` says " *
        "$n_total. The denominator is required precisely so it is stated once and checked, not " *
        "carried twice with two values."))

    n_abstained = sum(values(abstain_reason_counts); init = 0)
    return P14BatchDecision(results, Float64(alpha_fdr), Int(fdr.k_star), Float64(fdr.t_star),
                            Float64(fdr.cost_ratio), Float64(fdr.fdp), Vector{Int}(fdr.accepted),
                            Int(n_total), Int(fdr.n_decided), Int(n_abstained),
                            fdr.decided_fraction, abstain_reason_counts, class_prior_used,
                            class_prior_source, [e for e in pi_sensitivity], conformal,
                            P14_FDR_SCOPE, meta)
end

# --- The headline -------------------------------------------------------------------------------

"""
    _p14_headline_line(alpha_fdr, posterior_expected_fdp, n_decided, n_total) -> String

The single headline line, built from the four numbers that may never be separated.

Factored out of [`p14_headline`](@ref) so that [`p14_headline_template_probe`](@ref) exercises the
SAME code path with placeholder numbers, rather than a second format string that could drift from
the real one.
"""
function _p14_headline_line(alpha_fdr::Real, posterior_expected_fdp::Real,
                            n_decided::Integer, n_total::Integer)
    frac = n_total == 0 ? NaN : n_decided / n_total
    return string("alpha_FDR = ", alpha_fdr,
                  " | realized posterior expected FDP = ", round(posterior_expected_fdp; digits = 4),
                  " | decided ", n_decided, "/", n_total,
                  " (", round(100 * frac; digits = 1), "%)",
                  " | ", P14_HEADLINE_SCOPE_TOKEN)
end

"""
    p14_headline(b::P14BatchDecision) -> String

ONE line carrying the alpha, the realized posterior expected FDP, the decided fraction as both a
ratio and a percentage, and the literal token `FDR SCOPE: DECIDED SUBSET ONLY`.

**The report and every runner print THIS rather than assembling their own sentence.** That is the
whole design: `(alpha_fdr, n_decided/n_total)` is one quantity, and two numbers that live in one
string cannot be separated by a copy-paste. A rule that abstains on most of a batch hits any alpha
trivially, so an alpha quoted without its denominator is not a weaker claim -- it is a misleading
one.
"""
p14_headline(b::P14BatchDecision) =
    _p14_headline_line(b.alpha_fdr, b.posterior_expected_fdp, b.n_decided, b.n_total)

"""
    p14_headline_template_probe() -> String

The headline that [`p14_headline`](@ref) emits, computed through the SAME code path on placeholder
numbers, so a smoke check can assert the mandatory tokens are present without first assembling a
whole batch.

It is a PROBE, not a template string: if the real headline ever stopped carrying the scope token or
the decided fraction, this would stop carrying them too. A second hand-written format string would
have kept passing.
"""
p14_headline_template_probe() = _p14_headline_line(0.10, 0.0432, 37, 100)

# --- The named limits ---------------------------------------------------------------------------

"""
    p14_named_limits() -> NTuple{8,String}

The honesty items that must travel with EVERY Phase-14 result, enumerable from code:

 1. the decision layer is NOT shipped in v2.0 (D-01);
 2. the conformal guarantee is simulator-derived and the intended real-data bound is ABSENT
    (D-03/D-03a);
 3. the real-data check is six committed TIFFs -- an illustration, not a coverage claim (D-03a);
 4. FDR is controlled over the decided subset only, and the decided fraction travels with alpha;
 5. per-tile output is uncontrolled display (D-04);
 6. SC1 is amended by D-02 and SC2 by D-05, and the original wording must be cited alongside;
 7. the four SC3 bars and the coverage floor are judgement calls with no derivation;
 8. the underlying net is an epoch-4 checkpoint of an early-overfitting run (13-REPORT limit C).

Enumerated rather than remembered. A report that drops one is a report that failed to iterate a
tuple, which shows up in a diff; a report that forgot a sentence does not.
"""
p14_named_limits() = P14_NAMED_LIMITS
