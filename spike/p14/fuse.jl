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

# spike/p14/fuse.jl --- the abstention rule: D-05's ASYMMETRIC fusion over D-06's THREE-VALUED OOD.
#
# SC2 IS AMENDED BY D-05, AND THE AMENDMENT IS CITED HERE BECAUSE A READER OF THIS FILE IS EXACTLY
# THE READER WHO NEEDS IT. The ROADMAP's SC2 asks for
#
#     OOD  OR  cross-method disagreement  OR  ambiguous conformal set   =>   be silent
#
# and that literal OR is REPLACED by the rule below. Two reasons, both recorded so the original is
# never planned against again:
#   1. Phase 9 exists so that v2.0 can be positioned as "knows when the classics are wrong, not
#      merely agrees with them". Making disagreement with Costes/Manders a trigger for SILENCE would
#      mute the tool precisely in the cases that demonstrate its headline differentiator. The
#      disagreement is the PRODUCT THESIS, not a failure mode.
#   2. Three OR'd triggers COMPOUND. A tool that is usually silent is not useful, and SC2's
#      protective intent is fully preserved by keeping the two triggers that really are evidence
#      against ourselves (OOD, and a conformal set that does not separate the classes).
# ANY REPORT QUOTING A DECISION FROM THIS FILE MUST CITE THE ORIGINAL SC2 WORDING ALONGSIDE THE
# AMENDMENT.
#
# CONSIDERED AND NOT CHOSEN, recorded so a later reader knows it was weighed rather than missed: a
# FOUR-STATE output adding an explicit `DECIDED-CONTRA-CLASSICAL` action beside {decide, abstain}.
# It is the more transparent design and is the one to reach for if reviewers push back on the
# asymmetry. It was rejected for v2.0 as a non-standard output state needing extra defence; the same
# information is carried here as a REASON on a two-state action, which every existing consumer can
# read. See 14-CONTEXT.md Deferred Ideas.
#
# D-06, AND WHY IT IS A SCIENTIFIC FAILURE RATHER THAN A STYLE ONE. `ood_verdict`
# (src/amortized/ood.jl:394-396) computes
#
#     thr  = haskey(ood_nulls, :thr) ? ood_nulls.thr : nothing
#     flag = thr === nothing ? false : fused > thr
#
# so when `:thr` is absent it returns a FALSE flag FOR EVERY INPUT, UNCONDITIONALLY. The persisted
# `ood_nulls_G.jld2` carries only the frozen density fit and no threshold, so this is the DEFAULT
# shape, not an edge case. `src/amortized/local_map.jl:70-72` states it outright: "a `false` flag on
# a non-sentinel tile means NOT CHECKED, not in distribution", and `:110-119` names the consequence
# -- the flag is VACUOUS, fires only structurally, and SILENTLY UNDERSTATES RISK. Reading it as "in
# distribution" makes the tool speak confidently in exactly the regime where it has no evidence it
# should. The failure this prevents is silent and scientific; the cost it imposes is loud and
# operational, which is the right way round.
#
# `p14_ood_state` IS THE ONLY FUNCTION IN THE PHASE-14 SOURCES PERMITTED TO INSPECT A RAW VERDICT
# FLAG, because turning a raw flag into the three-valued state is its entire job.
# spike/test/test_p14_decoupling.jl's SC2-c scan is lane-wide and its exemption is POSITIONAL: it
# locates `function p14_ood_state` and its matching `end` at the same indentation and skips exactly
# those lines. THE LONG `function ... end` FORM IS THEREFORE REQUIRED -- a short-form definition
# would be un-delimitable, and the exemption would have to become a whole-file pass, which is not an
# exemption but a hole.
#
# DECOUPLING (hard constraint, CLAUDE.md and D-01): spike-local. Reaches src/ READ-ONLY for the
# `OODVerdict` type through a guarded include, touches no src/ byte, adds no dependency, writes
# nothing, opens no image and consumes no random stream at all.
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :p14_fuse) || include(joinpath(@__DIR__, "fuse.jl"))

# --- Guarded includes, in dependency order ------------------------------------------------------
# The Phase-14 Tier-1 pre-registration. Nothing in this file reads a THRESHOLD from it: the fusion
# rule is pure symbol arithmetic and has no tunable number in it at all, which is deliberate --
# there is no knob here to be adjusted after seeing a result.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
# READ-ONLY reach into src/ for `OODVerdict`. TWO levels up, not three.
isdefined(@__MODULE__, :OODVerdict) ||
    include(joinpath(@__DIR__, "..", "..", "src", "results.jl"))

# THE STATUS VOCABULARY IS DELIBERATELY DUPLICATED FROM spike/p14/conformal.jl, AND THE DUPLICATION
# IS SELF-POLICING. This file must stay loadable WITHOUT conformal.jl, because conformal.jl pulls
# the whole Phase-13 stack (Flux, NeuralEstimators, Images) through posterior.jl -- roughly sixteen
# seconds of package load -- while the rule below is pure symbol arithmetic that a caller must be
# able to test in a second. So both files declare the same three symbols under the same name behind
# the same guard, and whichever loads SECOND asserts agreement rather than silently losing. A
# textual divergence therefore fails loudly the moment both files are in one process, which
# `decide.jl` guarantees. Recorded here rather than resolved by an include, because the include
# would trade a loud, cheap check for a slow one.
if isdefined(@__MODULE__, :P14_CONFORMAL_STATUSES)
    @assert P14_CONFORMAL_STATUSES === (:singleton, :ambiguous, :empty) "spike/p14/fuse.jl: P14_CONFORMAL_STATUSES was already defined as $(P14_CONFORMAL_STATUSES), which disagrees with this file's (:singleton, :ambiguous, :empty). The conformal status vocabulary is duplicated across spike/p14/conformal.jl and this file on purpose (see the note above); a divergence must fail loudly here rather than resolve to whichever file loaded first."
else
    const P14_CONFORMAL_STATUSES = (:singleton, :ambiguous, :empty)
end

# THE GUARD COVERS ONLY THE `const`s. Julia 1.12 DROPS every docstring written inside an
# `if ... end` block, so the documented functions sit at top level; method redefinition under a
# re-include is silent and harmless, while `const` redefinition is what would warn.
if !isdefined(@__MODULE__, :P14_FUSE_LOADED)
    # THE THREE-VALUED OOD STATE, frozen as a tuple so the validation below enumerates it
    # mechanically rather than by hand. NOTE THAT IT IS NOT A Bool AND CANNOT BE COERCED TO ONE:
    # that is the whole of D-06.
    const P14_OOD_STATES = (:fired, :clear, :not_checked)

    # The CLOSED set of abstention reasons. A new reason cannot appear untracked, because the result
    # type and the report both enumerate this tuple rather than collecting whatever turns up.
    const P14_ABSTAIN_REASONS = (:ood_fired, :ood_not_checked, :conformal_ambiguous,
                                 :conformal_empty, :disagreement_and_ood)

    # The two DECIDE reasons. `:decided_contra_classical` is the one D-05 exists to make reachable.
    const P14_DECIDE_REASONS = (:decided, :decided_contra_classical)

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else can
    # make this body skip.
    const P14_FUSE_LOADED = true
end

# --- The three-valued OOD resolver ----------------------------------------------------------------

"""
    p14_ood_state(ood_nulls, verdict::OODVerdict) -> Symbol

Resolve the OOD channel into ONE OF THREE STATES — `:fired`, `:clear` or `:not_checked` — and
**never a `Bool`**.

# The mechanic this exists for (D-06)

`ood_verdict` computes `flag = thr === nothing ? false : fused > thr` (`src/amortized/ood.jl:396`).
So when `:thr` is absent from `ood_nulls`, it returns a **false verdict for EVERY input,
unconditionally** — and that is the *default* shape, because the persisted `ood_nulls_G.jld2`
carries only the frozen density fit and no operating point. A false verdict is therefore genuinely
two-valued at the source, and `src/amortized/local_map.jl:70-72` says so outright: *"a `false` flag
on a non-sentinel tile means NOT CHECKED, not in distribution"*, with `:110-119` naming the
consequence — the channel fires only structurally, never as a detector, **which silently understates
risk**.

# The acceptance predicate is byte-for-byte the shipped one

`(t isa Real && isfinite(t))` is exactly `_recorded_ood_threshold`'s test at
`src/amortized/local_map.jl:107`, and matching it is not cosmetic: a non-finite threshold is *also*
not-checked there, so a Phase-14 resolver that accepted `Inf` would silently disagree with the
shipped semantics while looking correct. `-Inf`, `NaN`, `nothing` and a non-`Real` all resolve to
`:not_checked` for the same reason.

# This is the ONLY place permitted to inspect the raw verdict

`spike/test/test_p14_decoupling.jl`'s SC2-c scan is lane-wide and exempts exactly this function
body, located POSITIONALLY from `function p14_ood_state` to its matching `end`. The long
`function ... end` form is required for that exemption to be delimitable at all.

# What is NOT wired in this lane, stated rather than papered over

The Phase-14 null is a single-channel (density) fit on the Phase-13 net's own pool. The noise and
posterior-predictive channels are **not** wired here, so `:clear` means "one channel says clear",
not "three channels say clear".
"""
function p14_ood_state(ood_nulls, verdict::OODVerdict)
    haskey(ood_nulls, :thr) || return :not_checked
    t = ood_nulls.thr
    (t isa Real && isfinite(t)) || return :not_checked
    return verdict.flag ? :fired : :clear
end

# --- The fusion -----------------------------------------------------------------------------------

"""
    p14_abstain_reasons() -> NTuple{5, Symbol}

The CLOSED set of abstention reasons: `(:ood_fired, :ood_not_checked, :conformal_ambiguous,
:conformal_empty, :disagreement_and_ood)`.

Enumerated rather than collected, so a new reason cannot appear untracked in an artifact and a
reason that stopped being reachable shows up as a dead member. Every ABSTAIN carries one of these:
**a silent abstention is nearly as bad as a wrong call**, because a reader cannot tell a tool that
declined from a tool that failed.
"""
p14_abstain_reasons() = P14_ABSTAIN_REASONS

"""
    p14_fuse(; ood_state, conformal_status, disagree, allow_unchecked_ood = false) -> NamedTuple

The D-05 ASYMMETRIC abstention rule. Returns `(action, reason)` with `action` in
`(:decide, :abstain)` and `reason` in `p14_abstain_reasons()` or `(:decided,
:decided_contra_classical)`.

# SC2's literal `∨` is REPLACED by this rule

ROADMAP SC2 asks for `OOD ∨ cross-method disagreement ∨ ambiguous set ⇒ be silent`. Under that OR,
the tool would abstain **exactly where Phase 9 says it should speak** — Phase 9 exists so that v2.0
can be positioned as *"knows when the classics are wrong, not merely agrees with them"*, so making
disagreement with Costes/Manders a trigger for silence mutes the headline differentiator. And three
OR'd triggers compound: a tool that is usually silent is not useful. D-05 keeps the two triggers
that really are evidence against ourselves and turns the third into a RECORDED FIELD:

| input | action |
|---|---|
| OOD fired | ABSTAIN |
| ambiguous conformal set | ABSTAIN |
| empty conformal set | ABSTAIN (a distinct channel — see below) |
| cross-method disagreement **alone** | **DECIDE**, recorded as `:decided_contra_classical` |
| cross-method disagreement **∧** OOD not clear | ABSTAIN (an independent reason to distrust ourselves) |
| OOD not-checked, no opt-out | ABSTAIN (D-06) |

*Considered and not chosen:* a FOUR-STATE output adding `DECIDED-CONTRA-CLASSICAL` as its own
action. It is the more transparent design and is the one to reach for if reviewers push back on the
asymmetry; it was rejected for v2.0 as a non-standard output state needing extra defence. The same
information is carried here as a reason on a two-state action.

# The order of the branches is LOAD-BEARING

Each branch cites the decision it implements, and they are evaluated in the order written:

 1. `:fired` — the detector answered and said no. Nothing downstream can overturn that.
 2. `:not_checked` without the explicit opt-out — D-06. Placed second so that a caller who opted out
    still gets every other trigger applied, rather than opting out of the whole rule.
 3. `:ambiguous` — the evidence does not separate two or more classes at the conformal level.
 4. `:empty` — **deliberately not folded into `:ambiguous`.** Both abstain, so folding them would
    change no decision and would destroy the most informative single output of the hedge: an empty
    set means the point is more nonconforming than `(1 − α)` of the entire calibration set, a
    distribution-free misspecification signal on a *completely different channel* from the
    Mahalanobis/PP/noise detector, and the only signal in this phase that does not need the
    simulator to be right about densities — only about exchangeability.
 5. `disagree && ood_state !== :clear` — reachable ONLY under the explicit opt-out, since
    `:fired` and un-opted-out `:not_checked` have already returned. Disagreement is no longer
    *alone* here: there is an independent reason to distrust ourselves.
 6. otherwise DECIDE, carrying the disagreement as `:decided_contra_classical` when there was one.

# Every ABSTAIN carries which trigger fired

That is a requirement, not a convenience. A silent abstention is nearly as bad as a wrong call: the
reader cannot distinguish a tool that declined from a tool that broke, and cannot tell an OOD
abstention (fix the input) from a conformal-empty one (question the simulator).

# Counter-examples that make the rule load-bearing

    p14_fuse(ood_state = :clear, conformal_status = :singleton, disagree = true).reason
        === :decided_contra_classical      # NEVER an abstention -- this cell IS D-05's amendment
    p14_fuse(ood_state = :not_checked, conformal_status = :singleton, disagree = false).reason
        === :ood_not_checked               # NEVER :decided -- a not-checked channel is not "clear"
"""
function p14_fuse(; ood_state::Symbol, conformal_status::Symbol, disagree::Bool,
                    allow_unchecked_ood::Bool = false)
    # VALIDATE FIRST, then decide. A bad symbol reaching the branch chain would fall through every
    # early return and be reported as a confident `:decided`, which is the worst available failure.
    ood_state in P14_OOD_STATES || throw(ArgumentError(
        "p14_fuse (D-06): ood_state must be one of $(P14_OOD_STATES) -- the OOD input is " *
        "THREE-valued, and a Bool cannot represent it, because a false verdict on an unchecked " *
        "channel means NOT CHECKED rather than in distribution. Got $ood_state."))
    conformal_status in P14_CONFORMAL_STATUSES || throw(ArgumentError(
        "p14_fuse: conformal_status must be one of $(P14_CONFORMAL_STATUSES); got " *
        "$conformal_status. The empty set is a status of its own and is never reported as " *
        "ambiguous (D-05)."))

    ood_state === :fired && return (action = :abstain, reason = :ood_fired)                 # D-05
    (ood_state === :not_checked && !allow_unchecked_ood) &&
        return (action = :abstain, reason = :ood_not_checked)                               # D-06
    conformal_status === :ambiguous &&
        return (action = :abstain, reason = :conformal_ambiguous)                           # D-05
    conformal_status === :empty &&
        return (action = :abstain, reason = :conformal_empty)   # D-05, a DISTINCT channel
    (disagree && ood_state !== :clear) &&
        return (action = :abstain, reason = :disagreement_and_ood)  # D-05, opt-out only
    return (action = :decide,
            reason = disagree ? :decided_contra_classical : :decided)  # D-05: the Phase-9 thesis
end
