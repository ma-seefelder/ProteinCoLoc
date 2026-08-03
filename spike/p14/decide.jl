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

# spike/p14/decide.jl --- the composition point: `decide_coloc`, the entry point the phase is
# named for.
#
# "SRC-SHAPED" IS NOT "IN SRC", AND THIS FILE IS WHERE THAT DISTINCTION IS EASIEST TO LOSE.
# D-01 says the decision layer is built in the SPIKE RESEARCH LANE against a signature it WOULD
# have in src/ -- a bundle plus images in, an AbstractColocResult subtype out -- so that a later
# promotion would be mechanical rather than a rewrite. It says nothing about shipping it, and the
# phase makes ZERO edits to src/. THE DECISION LAYER IS NOT SHIPPED IN v2.0. Any report quoting a
# number from this file must say so plainly; `spike/test/test_p14_decoupling.jl` asserts src/ is
# byte-unchanged and untracked-clean on every run, which is what makes that sentence executable
# rather than a promise.
#
# THE ORDERING CLAIM THIS FILE EXISTS TO IMPLEMENT: ABSTAIN FIRST, THEN FDR-SORT THE SURVIVORS.
# See `decide_coloc`'s docstring for Assumption A2 and its two-line derivation. The reverse order
# is wrong and the reason is recorded there so it is not re-derived at some later date by someone
# who assumes the order was arbitrary.
#
# WHAT THIS FILE DOES NOT DO. It does not train, does not fit a threshold, does not choose a
# number and does not write anything. Every constant it consumes was frozen before any Phase-14
# result existed (`spike/p14/consts.jl`), and the ONE number that is neither frozen here nor
# chosen -- tau -- is LOADED from the Phase-13 artifact with its provenance asserted (D-07).
#
# THE THREE SEEDS, KEPT APART. `P14_DEV_SEED` is the Phase-14 development stream. The comparator's
# own master seed IS `NPE_MASTER_SEED`, a FORBIDDEN stream (14-RESEARCH Pitfall 7), and an
# executor copying the Phase-9 Costes call site verbatim would consume it silently --
# `test_p14_consts.jl` cannot catch that, because it checks CONSTANTS and not CALL SITES. So the
# Costes call site here is guarded twice: `_p14_assert_costes_seed` refuses every seed in the
# frozen forbidden inventory at RUN TIME, and the lane grep in `test_p14_decoupling.jl` refuses
# the forbidden spelling at REVIEW TIME.
#
# THE SCALE, KEPT STRAIGHT. tau is defined on rho_true; `patch_correlation` returns the induced
# mean per-patch correlation mu. Comparing them directly would be a criterion applied in a unit it
# was not derived for, which is one of the four recorded "the bar was wrong" families in this
# project. `ghat` is the frozen monotone bridge and is applied on BOTH channels before either
# reaches the cut; `cross_method.basis` records which scale the comparison was made on, so the
# choice is auditable after the fact rather than implicit.
#
# THE ACCESSOR THE DECISION LAYER NEVER ROUTES THROUGH. The shipped Bool accessor `is_ood(` is
# lossy in a direction that matters (D-06): a `false` return collapses :clear and :not_checked,
# and reading the second as the first makes the tool speak confidently in exactly the regime where
# it has no evidence it should. This file reads the THREE-VALUED state through `p14_ood_state` and
# nothing else; the lane guard asserts the Bool accessor's call spelling appears here zero times.
# (Named in this comment block on purpose: the guard scans the COMMENT-STRIPPED source, so a
# docstring naming it would be code to that scan while a `#` block is not.)
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three.
#
# INCLUDING THIS FILE DOES NOTHING. The script entry point at the foot is guarded by
# `abspath(PROGRAM_FILE) == @__FILE__`, so a test or a runner that includes it gets definitions
# and no side effects.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :decide_coloc) || include(joinpath(@__DIR__, "decide.jl"))
#
# RUN THE ACCOMPANYING TEST PER FILE, never through the aggregate suite (which exits 1 early at
# the Phase-4 speedup gate and masks every later include block):
#
#     julia --project=spike spike/test/test_p14_decide.jl

# --- Guarded includes, in dependency order -------------------------------------------------------
# The Phase-14 pre-registration and the D-07 threshold loader come first; the rule modules next;
# the Phase-13 read surface and the classical comparator last, because they are the heaviest and
# because `spike/contract.jl` (which `spike/comparator/classical.jl` requires and deliberately
# does NOT include itself) is already in scope by then through the Phase-13 chain.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED)   || include(joinpath(@__DIR__, "provenance.jl"))
isdefined(@__MODULE__, :P14_POSTERIOR_LOADED)    || include(joinpath(@__DIR__, "posterior.jl"))
isdefined(@__MODULE__, :p14_bayes_fdr)           || include(joinpath(@__DIR__, "fdr.jl"))
isdefined(@__MODULE__, :p14_conformal_set)       || include(joinpath(@__DIR__, "conformal.jl"))
isdefined(@__MODULE__, :P14_FUSE_LOADED)         || include(joinpath(@__DIR__, "fuse.jl"))
isdefined(@__MODULE__, :P14Result)               || include(joinpath(@__DIR__, "result.jl"))
isdefined(@__MODULE__, :ThreeWayEvidenceNet) ||
    include(joinpath(@__DIR__, "..", "p13", "net.jl"))
isdefined(@__MODULE__, :three_way_label) ||
    include(joinpath(@__DIR__, "..", "p13", "labels.jl"))
isdefined(@__MODULE__, :load_p13_basis) ||
    include(joinpath(@__DIR__, "..", "p13", "preconditions.jl"))
# READ-ONLY reach into src/ for the SHIPPED OOD channel. Consumed unchanged; never reimplemented.
isdefined(@__MODULE__, :ood_verdict) ||
    include(joinpath(@__DIR__, "..", "..", "src", "amortized", "ood.jl"))
# The comparator's frozen pre-registration (COSTES_ALPHA, the scramble count and block size). The
# sentinel is deliberately the alpha rather than the comparator's master seed: naming that seed
# here would put a forbidden stream's name one keystroke from a call site (Pitfall 7).
isdefined(@__MODULE__, :COSTES_ALPHA) ||
    include(joinpath(@__DIR__, "..", "comparator", "config.jl"))
isdefined(@__MODULE__, :manders) ||
    include(joinpath(@__DIR__, "..", "comparator", "classical.jl"))
isdefined(@__MODULE__, :ghat) ||
    include(joinpath(@__DIR__, "..", "simulator", "ghat.jl"))

# THE GUARD BLOCK COVERS ONLY THE `const`s. Julia 1.12 DROPS every docstring written inside an
# `if ... end` block -- the parser emits the `Core.@doc` call but the docsystem never registers it
# -- so every documented function below sits at TOP LEVEL. Method redefinition under a re-include
# is silent and harmless; `const` redefinition is what would warn. Same split, same reason, as
# spike/p13/result.jl, spike/p14/posterior.jl and spike/p14/result.jl.
if !isdefined(@__MODULE__, :P14_DECIDE_LOADED)
    # The keys ONE item of a batch must carry, and the keys it MAY carry. Enumerated so a
    # misspelled optional key fails loudly instead of silently disabling the cross-method channel
    # -- a silent `false` disagreement reads as AGREEMENT, which is the one failure this phase
    # must not produce quietly.
    const P14_PAIR_REQUIRED_KEYS = (:Zs, :Zc)
    const P14_PAIR_OPTIONAL_KEYS = (:lambda, :mci_s, :mci_c, :idx, :local_map_uncontrolled)

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else
    # can make this body skip.
    const P14_DECIDE_LOADED = true
end

# --- The bundle: everything inherited, loaded once ------------------------------------------------

"""
    p14_load_bundle(; net_path, probe_path, gate_path) -> NamedTuple

Load the Phase-13 inheritance ONCE and return it as `(h, tau, prov, prior, calibration, grid)`.

**This is the only place in Phase 14 that obtains tau (D-07).** The threshold is LOADED from the
Phase-13 artifacts by [`p14_load_tau`](@ref), which asserts the four on-disk values agree and pins
the sha table on both sides, and the net's own recorded value is asserted equal to it again here.
There is deliberately **no keyword by which a caller can override tau**: copying the value is how
two numbers silently diverge, re-deriving it guarantees a second different cut on the same
hypothesis space, and no legitimate caller in this phase needs a third option.

# Fields

- `h`: the loaded three-way net handle (`net`, `head_log_odds`, `input_dim`, `num_summaries`,
  `width`, `tau`, `consts_sha`, `cut_variant`, `schema_version`, `meta`).
- `tau`: the inherited three-way cut, measured by Phase 13 and never tuned here.
- `prov`: the full D-07 provenance block, embedded in every artifact this phase writes.
- `prior`: `h.meta.pi_class_masses` -- the **MEASURED** realized class mix of the prior draws the
  net was trained under. Measured rather than chosen, which is why it is the default and why
  `class_prior_source` records when a caller replaced it.
- `calibration`: Phase 13's own `CalibrationMeta` for the coloc head, CARRIED FORWARD unchanged.
- `grid`: the patch grid the net was trained for, read off the artifact.

# The prior is the simulator's class mix, not any batch's prevalence

Stated here because it is the biggest way the downstream FDR guarantee fails silently. See
[`decide_coloc`](@ref) for the sensitivity curve that exists to surface it.
"""
function p14_load_bundle(; net_path = joinpath(P14_P13_DIR, "three_way_net.jld2"),
                           probe_path = joinpath(P14_P13_DIR, "tau_probe_report.jld2"),
                           gate_path = joinpath(P14_P13_DIR, "three_way_gate_report.jld2"))
    h = load_three_way(net_path)
    prov = p14_load_tau(; net_path = net_path, probe_path = probe_path, gate_path = gate_path)

    # ASSERTED AGAIN, ON THIS SIDE. `p14_load_tau` reopens the net itself, so this catches the one
    # case its four-way check cannot: two different net paths reaching one bundle.
    @assert h.tau == prov.tau "D-07 breach: the net handle this bundle carries records tau = $(h.tau) while the provenance block loaded from $(net_path) records $(prov.tau); the bundle would then be quoting one artifact and cutting with another"

    hasproperty(h.meta, :pi_class_masses) || error(
        "p14_load_bundle: the net at $net_path records no meta.pi_class_masses. The class prior " *
        "MUST be the measured class mix of the prior the net was trained under; substituting a " *
        "plausible one here would make it a chosen number, which is precisely what this " *
        "project's discipline forbids.")

    grid = hasproperty(h.meta, :grid) ? Int(h.meta.grid) : isqrt(Int(h.num_summaries))
    return (h = h, tau = prov.tau, prov = prov,
            prior = h.meta.pi_class_masses,
            calibration = _p14_carry_p13_calibration(prov.gate_path, grid),
            grid = grid)
end

"""
    _p14_carry_p13_calibration(gate_path, grid) -> CalibrationMeta

Phase 13's coloc-head calibration, **carried forward unchanged** from the persisted gate report.

# Why this is a carry and not a recomputation

The FDR guarantee is an identity GIVEN correct posteriors; nothing rescues it if the probabilities
are wrong, so every Phase-14 number must ship with the calibration evidence beside it. That
evidence already exists and was already gated: re-running the binning on Phase-14 draws would
produce a SECOND, differently-sourced ECE that a reader would inevitably conflate with the gated
one, and no reader could tell which number backed which claim. So the reliability curve, the ECE
and the MCE are read back off `three_way_gate_report.jld2` and pushed through the SHARED
`p13_calibration_meta` builder, which recomputes the empty-bin count, the vacuity label and the
band verdict from that carried curve.

Those three recomputed values are then asserted equal to the ones Phase 13 persisted. That is a
real check rather than a tautology: the curve is carried, but the verdicts are re-derived from it,
so a carry that landed on the wrong head, the wrong bins or the wrong AUC fails here instead of
travelling.

Only the COLOC head is carried, and deliberately: the composite null this phase sorts on is
`1 - P(coloc)`, so the coloc head is the one whose calibration the rule rests on. The exclusion
head's numbers stay in the Phase-13 report where they were gated.
"""
function _p14_carry_p13_calibration(gate_path, grid::Integer)
    rep = _p14_load_report(gate_path, "three-way gate")
    key(k) = _p14_report_key(rep, k, "three-way gate", gate_path)

    cal = CalibrationResult(Vector{Float64}(key("bin_midpoints_coloc")),
                            Vector{Float64}(key("predicted_rate_coloc")),
                            Vector{Float64}(key("observed_rate_coloc")),
                            Vector{Int}(key("bin_counts_coloc")),
                            Float64(key("ece_coloc")),
                            Float64(key("mce_coloc")))
    meta = p13_calibration_meta(cal; grid = grid, auc = Float64(key("head_auc_coloc")),
                                gate = (source = :p13_three_way_gate_report,
                                        head = :coloc,
                                        carried_forward = true,
                                        recomputed_in_phase14 = false,
                                        gate_path = gate_path))

    @assert meta.gate.empty_bins == Int(key("empty_bins_coloc")) "the carried Phase-13 reliability curve re-derives $(meta.gate.empty_bins) empty bins but the gate report recorded $(key("empty_bins_coloc")); the curve carried here is not the curve that was gated"
    @assert meta.gate.vacuous == key("vacuous_coloc") "the carried Phase-13 curve re-derives vacuous = $(meta.gate.vacuous) but the gate report recorded $(key("vacuous_coloc")); the ECE and the AUC carried here did not come from one another"
    @assert meta.gate.verdict === Symbol(key("verdict_coloc")) "the carried Phase-13 curve re-derives the band verdict $(meta.gate.verdict) but the gate report recorded $(key("verdict_coloc"))"
    return meta
end

"""
    p14_summary_length(bundle) -> Int

The length one standardized acquisition summary must have for this net: `2 * G^2`, DERIVED from
the net's own recorded `num_summaries`.

**Never written as a literal.** The realized input width is a function of the grid and of Phase
11's conditioning width; a hard-coded number would be silently wrong the day either moves, with
nothing thrown. `encode_conditioned_pair` re-asserts the assembled width against
`three_way_input_dim` for the same reason.
"""
p14_summary_length(bundle::NamedTuple) = 2 * Int(bundle.h.num_summaries)

# --- The OOD input shape --------------------------------------------------------------------------

"""
    p14_ood_input(ref) -> NamedTuple

Adapt what [`p14_ood_reference`](@ref) returns into the shape the SHIPPED `ood_verdict` reads.

# The mismatch this exists to bridge, recorded rather than patched upstream

`ood_verdict` reads its density fit as `ood_nulls.density` (`src/amortized/ood.jl:369`), while
`p14_ood_reference` returns the raw `fit_ood_nulls` tuple `(cont, muS, C)` merged with the
operating point. Those two shapes are not the same, and the difference is SILENT until the first
call: the fit is present under no key `ood_verdict` looks for. `spike/p14/pools.jl` is frozen and
is not edited to fix this; the adapter lives here, at the composition point, where the two
surfaces actually meet.

Accepts either a full reference NamedTuple or a NamedTuple that is already in verdict shape, and
returns the latter unchanged. `:thr` is carried ONLY when the reference actually fitted an
operating point -- its absence is what makes `p14_ood_state` resolve `:not_checked`, which is the
honest answer and the one that ABSTAINS by default (D-06). A substituted threshold would convert
a loud, correct abstention into a silent, confident answer.
"""
function p14_ood_input(ref::NamedTuple)
    haskey(ref, :density) && return ref
    haskey(ref, :nulls) || throw(ArgumentError(
        "p14_ood_input: expected either a verdict-shaped NamedTuple carrying `:density` or a " *
        "`p14_ood_reference` return carrying `:nulls`; got $(keys(ref))."))
    ref.nulls === nothing && return NamedTuple()
    base = (density = ref.nulls,)
    return (haskey(ref, :thr) && ref.thr !== nothing) ? merge(base, (thr = ref.thr,)) : base
end

"""
    _p14_ood_verdict(ood_nulls, Zs, Zc) -> Union{OODVerdict, Nothing}

The SHIPPED `ood_verdict`, consumed unchanged, or `nothing` when no density fit is present.

`ood_verdict` requires `ood_nulls.density` and would raise on a null that was never fitted. A
missing fit is not an error in this phase -- it is the NOT-CHECKED state -- so it is resolved to
`nothing` here and turned into `:not_checked` by the caller, rather than crashing a batch whose
correct answer is "abstain, because nothing checked this".

The posterior-predictive channel is switched OFF (`with_pp = false`): it needs re-simulation and a
model handle this lane does not carry. Only the summary-density channel is wired, so a `:clear`
here means ONE channel says clear, not three.
"""
function _p14_ood_verdict(ood_nulls, Zs, Zc)
    (ood_nulls isa NamedTuple && haskey(ood_nulls, :density) && ood_nulls.density !== nothing) ||
        return nothing
    return ood_verdict(ood_nulls, Zs, Zc; with_pp = false)
end

# --- The amortized call, BY NAME ------------------------------------------------------------------

"""
    _p14_amortized_call(p) -> ThreeWayClass

The class the amortized posterior calls, computed as a THREE-TERM comparison over NAMED fields --
never `argmax` over a positional collection.

# Why a positional read is the live hazard here

FIVE class orderings coexist in this codebase: `three_way_probs` emits `(coloc, exclusion)`,
`three_way_log_bf` emits `(coloc, exclusion, random)`, `ThreeWayLogBF` stores
`(coloc, random, exclusion)`, the gate's confusion order is `["exclusion", "random", "coloc"]`,
and `pi_class_masses` is keyed `(exclusion, random, coloc)`. Every aggregate statistic this phase
reports -- ECE, AUC, FDR -- is invariant to a CONSISTENT relabelling, which is exactly what makes
an inconsistent one so dangerous: it produces numbers that look right. Phase 12 lost a headline
finding to an orthogonal permutation for the same reason.

# The tie rule, stated rather than left to the comparison operator

A class wins only on a STRICT majority; an exact tie falls to `RANDOM`. That mirrors the frozen
`three_way_label`, whose dead zone absorbs everything ambiguous, and it errs toward the null,
which is the right direction for a rule whose false discoveries are the controlled quantity. Ties
are measure-zero for a continuous posterior, so this costs nothing and removes a silent choice.
"""
function _p14_amortized_call(p::NamedTuple)
    _p14_check_class_keys(p, P14_CLASS_KEYS, "_p14_amortized_call: `p`")
    (p.coloc > p.random && p.coloc > p.exclusion) && return COLOC
    (p.exclusion > p.random && p.exclusion > p.coloc) && return EXCLUSION
    return RANDOM
end

# --- The cross-method channels (D-05) -------------------------------------------------------------

"""
    _p14_assert_costes_seed(seed) -> nothing

Refuse, at RUN TIME, any Costes seed that lies in the frozen forbidden-seed inventory.

The comparator's own master seed IS `NPE_MASTER_SEED`, one of the pre-registered streams this
project's seed discipline forbids a research lane from consuming, and it is the FIRST POSITIONAL
ARGUMENT of the Costes p-value function -- so copying the Phase-9 call site verbatim consumes it
with nothing thrown. `spike/test/test_p14_consts.jl` cannot catch that: it asserts CONSTANTS
disjoint, and the constant it asserts disjoint is exactly the one that would be passed here.

This assertion is the run-time half of the guard; the lane grep in
`spike/test/test_p14_decoupling.jl` is the review-time half. Two halves because the failure is
invisible in the output: a forbidden stream produces perfectly ordinary-looking p-values.
"""
function _p14_assert_costes_seed(seed)
    UInt64(seed) in _p13_forbidden() && throw(ArgumentError(
        "p14_cross_method (Pitfall 7): the Costes seed $(repr(seed)) is in the frozen " *
        "forbidden-seed inventory. A research lane may not consume a pre-registered stream; " *
        "pass P14_DEV_SEED, which is the default and is recorded on every result as " *
        "`cross_method.costes_seed`."))
    return nothing
end

"""
    p14_cross_method(mci_s, mci_c, amortized_call, tau; idx, seed = P14_DEV_SEED) -> NamedTuple

The classical comparison, on the SAME SCALE as the amortized call, recorded as a named field
rather than suppressed (D-05).

# Two INDEPENDENT channels, and why each is constructed the way it is

**Channel 1 -- the classical three-way call.** `ghat(patch_correlation(mci))` on both the sample
and the control, then the FROZEN `three_way_label` with the INHERITED cut. Each link was checked:
`patch_correlation` is the mean per-patch correlation over the frozen 8x8 grid; `mbar_from_summary`
-- the statistic the cut was MEASURED on -- is the mean over present patches, i.e. **the same
quantity**; and `ghat` is the frozen, monotone, clamped piecewise-linear map from that mu to
rho_true, already running in the direction needed, so no inversion is written here. Applying the
cut directly to mu would be a criterion applied in a unit it was not derived for, which is one of
the four recorded "the bar was wrong" families in this project. No second threshold is invented.

**Channel 2 -- Costes significance.** A block-permutation test, genuinely independent of channel 1
rather than a second correlation cut, at the comparator's own frozen alpha. The seed is the
Phase-14 development stream and is asserted against the forbidden inventory before the call.

**Manders is recorded and NOT cut on.** `M1`/`M2` are fraction-of-signal-above-threshold
quantities on `[0, 1]` with no map to rho_true, so any cut on them would be a NEW CHOSEN CONSTANT
-- exactly what this phase must avoid. They travel as a descriptive column so a reader can see
them; no call is derived from them.

# Disagreement DECIDES, it does not silence (D-05, and this is the product thesis)

The literal reading of the original SC2 makes disagreement with the classics a trigger for
silence, which would mute the tool precisely in the cases that demonstrate its headline
differentiator. D-05 amends that: disagreement alone DECIDES and is recorded here; only
disagreement TOGETHER WITH a non-clear OOD state abstains. Any report quoting a decision must cite
the original SC2 wording alongside the amendment.

# `basis` is not decoration

It records WHICH SCALE the comparison was made on. That is the thing that is otherwise forgotten,
and a cross-method record without it is a comparison nobody can reconstruct later; `P14Result`'s
constructor refuses one that lacks it.
"""
function p14_cross_method(mci_s, mci_c, amortized_call::ThreeWayClass, tau::Real;
                          idx::Integer, seed = P14_DEV_SEED)
    _p14_assert_costes_seed(seed)

    # CHANNEL 1. `ghat` FIRST, on both acquisitions, then the inherited cut -- one scale, one cut.
    rho_hat_s = ghat(patch_correlation(mci_s))
    rho_hat_c = ghat(patch_correlation(mci_c))
    classical_call = three_way_label(rho_hat_s, rho_hat_c; tau = tau)
    disagree_label = classical_call !== amortized_call

    # CHANNEL 2. `seed` defaults to P14_DEV_SEED and was asserted non-forbidden above (Pitfall 7).
    p_s = costes_p(seed, idx, mci_s)  # `seed` defaults to P14_DEV_SEED; asserted above (Pitfall 7)
    costes_call = (isfinite(p_s) && p_s < COSTES_ALPHA) ? :coloc : :not_coloc
    disagree_costes = (costes_call === :coloc) != (amortized_call === COLOC)

    m1, m2 = manders(mci_s)   # DESCRIPTIVE ONLY -- no call is derived from a Manders coefficient

    return (classical_call   = classical_call,
            costes_p_sample  = p_s,
            costes_call      = costes_call,
            rho_hat_sample   = rho_hat_s,
            rho_hat_control  = rho_hat_c,
            manders_sample   = (M1 = m1, M2 = m2),
            disagree_label   = disagree_label,
            disagree_costes  = disagree_costes,
            disagree_any     = disagree_label || disagree_costes,
            basis            = :ghat_rho_true_scale,
            tau              = Float64(tau),
            costes_seed      = UInt64(seed),
            computed         = true)
end

"""
    _p14_no_cross_method(note) -> NamedTuple

The cross-method record for a pair whose images were not supplied, carrying `basis =
:not_computed` and the reason.

**Never a bare `false`.** A `disagree_any = false` with no note reads as AGREEMENT -- the classics
were consulted and concurred -- when in fact nothing was consulted at all. That is the `bad(note)`
early-return idiom of the Phase-13 composition lane, and it exists so a batch summary cannot
report "no disagreements" for a run in which no comparison happened.
"""
_p14_no_cross_method(note::AbstractString) =
    (classical_call = nothing, costes_p_sample = NaN, costes_call = nothing,
     rho_hat_sample = NaN, rho_hat_control = NaN, manders_sample = nothing,
     disagree_label = false, disagree_costes = false, disagree_any = false,
     basis = :not_computed, tau = NaN, costes_seed = nothing, computed = false, note = note)

"""
    _p14_recall_disagreement(cm, amortized_call) -> Bool

Re-derive `disagree_any` from a stored cross-method record against a (possibly different)
amortized call.

Needed by the prior-sensitivity sweep: the classical call and the Costes verdict do NOT depend on
the class prior, but the AMORTIZED call does, so a sweep that reused the stored `disagree_any`
would be comparing the classics against a call the rule no longer makes. Recomputing from the two
stored classical verdicts keeps the sweep a genuine re-run of the whole rule rather than a partial
one.
"""
function _p14_recall_disagreement(cm::NamedTuple, amortized_call::ThreeWayClass)
    cm.basis === :not_computed && return false
    dl = cm.classical_call !== amortized_call
    dc = (cm.costes_call === :coloc) != (amortized_call === COLOC)
    return dl || dc
end

# --- The per-pair path ----------------------------------------------------------------------------

"""
    p14_decide_one(bundle, Zs, Zc, ood_nulls, qhat; lambda, idx, mci_s = nothing,
                   mci_c = nothing, prior = bundle.prior, allow_unchecked_ood = false,
                   local_map_uncontrolled = nothing, posterior = zeros(0, 0),
                   costes_seed = P14_DEV_SEED, meta = NamedTuple()) -> P14Result

One image pair, end to end: the evidence net, the class posterior, the conformal hedge, the
three-valued OOD state, the two cross-method channels, and the D-05 fusion, composed into ONE
result constructed once at the end.

# Order of operations

 1. the pair is encoded through Phase 11's OWN conditioning encoder and pushed through the
    evidence net in a single forward pass; the returned triple is normalised through
    `_p13_as_three_way_logbf` rather than built positionally, because the read surface emits
    `(coloc, exclusion, random)` while the stored type is `(coloc, random, exclusion)`;
 2. the three-class posterior and the composite null `v = 1 - P(coloc)`;
 3. the OOD channel, resolved to the THREE-VALUED state -- fired / clear / not-checked;
 4. the split-conformal set at the pre-registered miscoverage level;
 5. the amortized call, by name;
 6. the two cross-method channels, when both images were supplied;
 7. the D-05 asymmetric fusion;
 8. the result.

# THE COLOC SIDE OF `decision` IS PROVISIONAL

`decision` is `:abstain` when the rule abstained. Otherwise it is `:coloc` when the amortized call
is the coloc class and `:not_coloc` when it is not -- **but on the coloc side that is a provisional
call, not a discovery.** Bayesian FDR is a BATCH quantity: which `:coloc` calls are ACCEPTED is
decided by [`decide_coloc`](@ref)'s prefix rule over the decided subset, which rewrites this field
and records the provisional value in `meta.provisional_decision`. A single result read outside a
batch therefore carries an evidence call, not an FDR-controlled discovery, and the two must not be
quoted as the same thing.

# Named limits (non-negotiable honesty — carried into the manuscript)

`p14_named_limits()` returns the eight items verbatim and a report enumerates that tuple rather
than remembering it. The four that bind THIS function: the decision layer is built in the spike
research lane and is NOT shipped in v2.0 (D-01); the conformal guarantee is simulator-derived and
inherits the simulator's misspecification in full (D-02/D-03/D-03a); the per-tile map carried in
`local_map_uncontrolled` is UNCONTROLLED DISPLAY and never an FDR-controlled call (D-04); and the
underlying net is an epoch-4 checkpoint of an early-overfitting run, a limit every number here
inherits unchanged.

Two more that are properties of this call: only the summary-density OOD channel is wired in this
lane, so `:clear` means one channel says clear rather than three; and the evidence path produces
NO posterior draws, so `posterior` defaults to an EMPTY matrix, which means "no draws were
produced" and never "zero draws were interesting".

# Returns

A [`P14Result`](@ref). Read it through the accessors rather than the fields: `decision`,
`abstain_reason`, `ood_state`, `class_posterior`, `null_posterior`, `null_split`,
`conformal_set` / `conformal_status`, `cross_method`, `log_bf_vs_random`, and `local_map` --
whose companion `p14_local_map_is_controlled` is always `false`, by design.
"""
function p14_decide_one(bundle::NamedTuple, Zs::AbstractVector, Zc::AbstractVector,
                        ood_nulls, qhat::Real;
                        lambda, idx::Integer,
                        mci_s = nothing, mci_c = nothing,
                        prior::NamedTuple = bundle.prior,
                        allow_unchecked_ood::Bool = false,
                        local_map_uncontrolled = nothing,
                        posterior::AbstractMatrix{<:Real} = zeros(Float64, 0, 0),
                        costes_seed = P14_DEV_SEED,
                        meta::NamedTuple = NamedTuple())
    # --- validation first (the house order: name, requirement, observed value) ---------------
    n_expected = p14_summary_length(bundle)
    for (nm, Z) in (("Zs", Zs), ("Zc", Zc))
        length(Z) == n_expected || throw(ArgumentError(
            "p14_decide_one: `$nm` must hold $n_expected entries (2*G^2 for this net, DERIVED " *
            "from its recorded num_summaries); got $(length(Z))."))
        all(isfinite, Z) || throw(ArgumentError(
            "p14_decide_one: `$nm` carries a non-finite entry. A degenerate summary scored " *
            "against the density null yields a number with no interpretation."))
    end
    isfinite(qhat) || throw(ArgumentError(
        "p14_decide_one: `qhat` must be finite; got $qhat. A non-finite conformal threshold " *
        "silently admits every class or none, and both look like a legitimate regime."))

    # 1. ONE forward pass. The conditioning rides Phase 11's OWN encoder -- never a raw lambda
    #    appended by hand, which would feed the net a value it never saw at a row it reads, with
    #    nothing thrown.
    Z = p13_encode_pair(Zs, Zc, lambda)
    lbf = _p13_as_three_way_logbf(three_way_log_bf(bundle.h.net, Z, bundle.h.head_log_odds))

    # 2. The posterior and the composite null.
    p = p14_class_posterior(lbf, prior)
    v = p14_null_posterior(p)
    split = p14_null_split(p)

    # 3. The OOD channel, THREE-VALUED. `p14_ood_state` is the only route (D-06).
    verdict = _p14_ood_verdict(ood_nulls, Zs, Zc)
    state = verdict === nothing ? :not_checked : p14_ood_state(ood_nulls, verdict)

    # 4. The distribution-free hedge.
    cset = p14_conformal_set(p, qhat)

    # 5. The call, by name.
    call = _p14_amortized_call(p)

    # 6. The classics -- both images or neither, and a NOTE rather than a silent `false`.
    cm = (mci_s !== nothing && mci_c !== nothing) ?
         p14_cross_method(mci_s, mci_c, call, bundle.tau; idx = idx, seed = costes_seed) :
         _p14_no_cross_method("no image pair was supplied to this call, so no classical " *
                              "comparison was made; this is NOT agreement with the classics")

    # 7. The D-05 asymmetric fusion, consumed unchanged.
    fused = p14_fuse(; ood_state = state, conformal_status = cset.status,
                       disagree = cm.disagree_any, allow_unchecked_ood = allow_unchecked_ood)

    dec = fused.action === :abstain ? :abstain : (call === COLOC ? :coloc : :not_coloc)
    reason = fused.action === :abstain ? fused.reason : nothing

    tw = ThreeHypothesisColocResult(bundle.grid, posterior, lbf,
                                    verdict === nothing ? OODVerdict(NaN, false, NamedTuple()) :
                                                          verdict,
                                    bundle.calibration,
                                    (; ood_state = state,
                                       ood_channels_wired = (:density,),
                                       ood_channels_not_wired = (:pp, :noise)))

    # 8. Constructed ONCE, at the end, with a provenance meta.
    return P14Result(tw, p, v, split, dec, reason, state, cset, cm, local_map_uncontrolled,
                     merge(meta, (idx = Int(idx),
                                  lambda = lambda,
                                  amortized_call = call,
                                  fuse_reason = fused.reason,
                                  allow_unchecked_ood = allow_unchecked_ood,
                                  class_prior_used = prior,
                                  decision_is_provisional_on_the_coloc_side = true,
                                  provenance = p14_provenance_record(bundle.prov))))
end

# --- The src-shaped per-pair entry point ----------------------------------------------------------

"""
    decide_coloc(bundle, img::MultiChannelImage, control::MultiChannelImage,
                 channels::AbstractVector{<:Integer}; ood_nulls, qhat,
                 lambda = P13_PHASE11_REFERENCE_LAMBDA, idx = 1, basis = nothing,
                 prior = bundle.prior, allow_unchecked_ood = false,
                 local_map_uncontrolled = nothing) -> P14Result

The **src-shaped** per-pair entry point: a bundle plus two images and the two channels to compare,
in; an `AbstractColocResult` subtype, out.

# "src-shaped" is not "in src"

The signature deliberately mirrors the shipped `colocalization_amortized(img, control, channels;
...)` so that a later promotion would be mechanical rather than a rewrite (D-01). **It is built in
the spike research lane and is NOT shipped in v2.0**, no byte of `src/` is edited by this phase,
and the lane guard asserts that on every run. A report must say so plainly.

# Arguments

- `bundle`: the [`p14_load_bundle`](@ref) handle -- the net, the inherited cut, its provenance,
  the measured class prior and Phase 13's carried calibration.
- `img`: the sample image.
- `control`: the control image.
- `channels`: the TWO channels to compare (e.g. `[1, 2]`).

# Keywords

- `ood_nulls`: the density fit plus its operating point, in the shape [`p14_ood_input`](@ref)
  produces. Omitting the operating point resolves to `:not_checked` and therefore ABSTAINS, which
  is the honest default (D-06).
- `qhat`: the split-conformal threshold, cut at the pre-registered miscoverage level. An ARGUMENT,
  never read from a constant here, and never derived from the FDR level -- the two are distinct
  quantities that happen to share a value at one grid point, which is exactly why (D-07).
- `lambda`: the registration-uncertainty level. Defaults to the reference value the cut was
  measured at. **Phase 11 closed negative on inferring registration from the summary**, so this is
  a quantity a caller measures externally rather than one this tool recovers; the default is a
  stated reference, not an estimate.
- `basis`: the frozen Phase-11 standardizer. Resolved AFTER validation so a bad `channels` throws
  without first reading a five-megabyte artifact.

# Named limits (non-negotiable honesty — carried into the manuscript)

`p14_named_limits()` returns all eight verbatim. The ones this entry point must not be quoted
without: the decision layer is NOT shipped in v2.0 (D-01); the conformal guarantee is
SIMULATOR-DERIVED and the intended real-data bound is absent rather than merely loose
(D-02/D-03/D-03a); the per-tile map is UNCONTROLLED DISPLAY, never an FDR-controlled call (D-04);
the original SC1 and SC2 are AMENDED by D-02 and D-05 and the original wording must be cited
alongside; and the underlying net is an epoch-4 checkpoint of an early-overfitting run. To those:
a single pair carries an evidence call, and a DISCOVERY is only defined by the batch rule.

# Returns

A [`P14Result`](@ref) -- a subtype of the shipped `AbstractColocResult`. Read it through
`decision`, `abstain_reason`, `ood_state`, `class_posterior`, `null_posterior`, `null_split`,
`conformal_set`, `conformal_status`, `cross_method`, `log_bf_vs_random` and `local_map`.
"""
function decide_coloc(bundle::NamedTuple, img::MultiChannelImage, control::MultiChannelImage,
                      channels::AbstractVector{<:Integer};
                      ood_nulls, qhat::Real,
                      lambda = P13_PHASE11_REFERENCE_LAMBDA,
                      idx::Integer = 1,
                      basis = nothing,
                      prior::NamedTuple = bundle.prior,
                      allow_unchecked_ood::Bool = false,
                      local_map_uncontrolled = nothing,
                      costes_seed = P14_DEV_SEED,
                      meta::NamedTuple = NamedTuple())
    # --- input validation (V5), BEFORE any artifact is opened -------------------------------
    length(channels) == 2 || throw(ArgumentError(
        "decide_coloc: `channels` must select exactly 2 channels, got $(length(channels))."))
    for (nm, m) in (("img", img), ("control", control))
        for c in channels
            (1 <= c <= length(m.data)) || throw(ArgumentError(
                "decide_coloc: channel $c out of range for $nm (has $(length(m.data)) channels)."))
        end
    end
    channels[1] == channels[2] && throw(ArgumentError(
        "decide_coloc: `channels` must select two DISTINCT channels; got $(collect(channels)). " *
        "A channel compared against itself has correlation 1 by construction."))

    b = basis === nothing ? load_p13_basis() : basis

    # The frozen ingestion, byte-for-byte the path the net was trained on and the path the
    # calibration was measured through: select the pair, encode, standardize through the
    # INHERITED Phase-11 standardizer. No crop, resize or pad step exists anywhere in it.
    pair_s = _p14_channel_pair(img, channels)
    pair_c = _p14_channel_pair(control, channels)
    Zs = standardize_summary(encode_d01(patch_summary(pair_s)), b.zt, b.variant)
    Zc = standardize_summary(encode_d01(patch_summary(pair_c)), b.zt, b.variant)

    return p14_decide_one(bundle, Zs, Zc, ood_nulls, qhat;
                          lambda = lambda, idx = idx, mci_s = pair_s, mci_c = pair_c,
                          prior = prior, allow_unchecked_ood = allow_unchecked_ood,
                          local_map_uncontrolled = local_map_uncontrolled,
                          costes_seed = costes_seed,
                          meta = merge(meta, (channels = collect(Int, channels),
                                              basis_provenance = :p13_inherited_phase11_zt)))
end

"""
    _p14_channel_pair(m::MultiChannelImage, channels) -> MultiChannelImage

The two selected channels of an image, as the 2-channel container every frozen statistic in this
lane reads.

`patch_summary`, `patch_correlation`, the Manders coefficients and the Costes test all read
`data[1]` and `data[2]` unconditionally, so the channel SELECTION has to happen once, here, rather
than being re-expressed at four call sites where three of them could drift. Built through the
lane's own constructor, so the Otsu thresholds the Manders coefficients need are computed the same
way the frozen encoder computes them.
"""
_p14_channel_pair(m::MultiChannelImage, channels::AbstractVector{<:Integer}) =
    build_mci(Vector{Matrix{Float64}}([Matrix{Float64}(m.data[channels[1]]),
                                       Matrix{Float64}(m.data[channels[2]])]);
              name = m.name)

# --- The batch path: ABSTAIN FIRST, THEN SORT -----------------------------------------------------

"""
    _p14_rule_at_prior(recs, qhat, prior, alpha_fdr, allow_unchecked_ood, n_total) -> NamedTuple

The WHOLE rule -- posterior, hedge, fusion, abstention partition, FDR prefix -- at one class
prior.

Factored out so the headline run and every rung of the prior-sensitivity curve go through **one
code path**. A sweep written as a second, simpler loop would keep agreeing with the headline right
up until the day it stopped, and nothing would say which of the two was wrong.

The prior enters in more places than the null posterior: it moves the class posterior, hence the
conformal set (hence the abstention), hence the amortized call (hence the cross-method
disagreement), hence the fusion. All of it is recomputed here; only the OOD state and the two
CLASSICAL verdicts are carried, because those are the only quantities in the rule that do not
depend on the prior at all.
"""
function _p14_rule_at_prior(recs, qhat::Real, prior::NamedTuple, alpha_fdr::Real,
                            allow_unchecked_ood::Bool, n_total::Integer)
    ps    = [p14_class_posterior(r.lbf, prior) for r in recs]
    calls = [_p14_amortized_call(p) for p in ps]
    csets = [p14_conformal_set(p, qhat) for p in ps]
    dis   = [_p14_recall_disagreement(recs[i].cm, calls[i]) for i in eachindex(recs)]
    acts  = [p14_fuse(; ood_state = recs[i].ood, conformal_status = csets[i].status,
                        disagree = dis[i], allow_unchecked_ood = allow_unchecked_ood)
             for i in eachindex(recs)]

    # STEP A -- ABSTAIN FIRST. The partition happens before anything is sorted.
    decided_idx = [i for i in eachindex(acts) if acts[i].action === :decide]
    # STEP B -- THEN SORT, over the DECIDED SUBSET ONLY.
    v = [p14_null_posterior(ps[i]) for i in decided_idx]
    fdr = p14_fdr_over_decided(v, alpha_fdr; n_total = n_total)
    return (ps = ps, calls = calls, csets = csets, disagree = dis, acts = acts,
            decided_idx = decided_idx, v = v, fdr = fdr)
end

"""
    _p14_reason_counts(results) -> NamedTuple

How many items each trigger silenced, keyed in the frozen `p14_abstain_reasons()` order.

Built by ENUMERATING the closed reason set rather than by collecting whatever turned up, so a
reason that stopped being reachable shows up as an absent key and a reason that is not in the set
cannot appear at all. The batch constructor then derives `n_abstained` from this ledger, which
turns the identity `n_decided + n_abstained == n_total` into a real check that every silenced item
was explained.
"""
function _p14_reason_counts(results)
    counts = Dict{Symbol, Int}()
    for r in results
        r.decision === :abstain || continue
        counts[r.abstain_reason] = get(counts, r.abstain_reason, 0) + 1
    end
    ks = Tuple(k for k in p14_abstain_reasons() if haskey(counts, k))
    return NamedTuple{ks}(Tuple(counts[k] for k in ks))
end

"""
    _p14_prior_at(prior, pi_c) -> NamedTuple

The measured class prior with its coloc mass replaced by `pi_c`, the other two renormalised IN
THEIR MEASURED RATIO.

The sweep asks "what if the batch's coloc prevalence were different", not "what if the whole class
structure were different", so the exclusion-to-random ratio the simulator realized is held fixed
and only the coloc mass moves. Rebuilding the two remaining masses in equal parts instead would
silently sweep two things at once and the curve would answer a question nobody asked.
"""
function _p14_prior_at(prior::NamedTuple, pi_c::Real)
    (0.0 < pi_c < 1.0) || throw(ArgumentError(
        "_p14_prior_at: the swept coloc mass must lie strictly inside (0,1); got $pi_c."))
    rest = prior.random + prior.exclusion
    rest > 0 || throw(ArgumentError(
        "_p14_prior_at: the measured non-coloc mass is $rest; the ratio it is renormalised in " *
        "is not defined."))
    s = 1.0 - pi_c
    return (coloc = Float64(pi_c),
            random = s * prior.random / rest,
            exclusion = s * prior.exclusion / rest)
end

"""
    decide_coloc(bundle, pairs, ood_nulls, qhat; alpha_fdr, prior = bundle.prior,
                 allow_unchecked_ood = false, pi_grid = P14_PI_COLOC_GRID, idx_offset = 0,
                 meta = NamedTuple()) -> P14BatchDecision

The batch entry point the phase is named for: `{coloc / not_coloc / ABSTAIN}` over a batch at a
**user-set** Bayesian FDR, with the decided fraction and the prior-sensitivity curve as required
output rather than optional extras.

`pairs` is a vector of NamedTuples carrying `Zs` and `Zc` (required) and any of `lambda`, `mci_s`,
`mci_c`, `idx` and `local_map_uncontrolled` (optional). An unknown key throws rather than being
ignored, because a misspelled `mci_s` would silently disable the cross-method channel and the
resulting `disagree_any = false` reads as agreement.

`alpha_fdr` is a REQUIRED keyword. SC1 calls it *"a user-set Bayesian FDR"*, so it is a per-call
user parameter: it never defaults to a constant of the pre-registration and is never derived from
the conformal miscoverage level, which is a different quantity that happens to share a value at
one grid point of the pre-registered sweep (D-07).

# ABSTAIN FIRST, THEN FDR-SORT THE SURVIVORS. The order is the point.

 1. every pair is decided individually;
 2. **the batch is partitioned into decided and abstained, and the reason ledger is tallied;**
 3. **only then** are the survivors' null posteriors sorted and cut by the prefix rule;
 4. the accepted positions are mapped back to ORIGINAL batch indices and asserted non-abstaining
    -- an off-by-one there is silent and would report the wrong items as discoveries;
 5. the decided calls are FINALISED: an accepted item is `:coloc`, a decided-but-not-accepted item
    is `:not_coloc`. The results are REBUILT, not mutated, and each carries its provisional call in
    `meta.provisional_decision` so a divergence between the argmax and the rule is auditable rather
    than silent.

## Assumption A2 -- why abstaining first needs no selective-inference correction

The controlled quantity is `E[FDP | data]`. The rejection set `R` is a deterministic function of
the observed data -- the abstention filter and the sort are both data-measurable -- so `R` is
sigma(data)-measurable and pulls straight out of the conditional expectation:

    E[FDP | data] = (1/|R|) * sum over i in R of P(H0_i | data) = (1/|R|) * sum over i in R of v_i

and that holds **for ANY data-dependent selection of `R`**, including one that filtered on OOD
status or conformal ambiguity. Unlike frequentist FDR, where a data-dependent pre-selection is a
genuine selective-inference problem requiring correction, the Bayesian posterior quantity is
immune because the conditioning has already happened.

This is a two-line derivation stated as such, not a quoted theorem. **It is load-bearing** -- it
is the entire justification for abstain-then-sort and for the decided-subset scope -- so the
report must STATE it explicitly and let a referee check it, rather than assert it and be believed.

## Why the REVERSE order is wrong, recorded so it is not re-derived

Sorting the whole batch and then abstaining from some of the ACCEPTED items removes items from a
set whose running mean was computed over all of them, changing both the numerator and the
denominator in an uncontrolled way; the residual accepted set can then have a running mean that
EXCEEDS the level. Abstaining first also improves the premise, because it routes away precisely
the items whose posteriors are least likely to be right.

# The honest claim, and the mandatory companion number

> The posterior expected false discovery proportion is controlled at the user-set level **over the
> DECIDED subset only** -- the batch minus abstentions. It is **not** controlled over the whole
> batch: an abstained item is neither a discovery nor a non-discovery and appears in neither the
> numerator nor the denominator. The guarantee is conditional on the model (calibrated posteriors
> and a correct class prior) and is a posterior expectation, not a frequentist long-run rate.

A rule that abstains on most of a batch hits any level trivially, so the DECIDED FRACTION travels
with the level as one quantity. `p14_headline` emits both in a single string; quote that rather
than assembling a sentence.

# The prior-sensitivity curve is REQUIRED output

`v_i` depends on the class prior, and the default prior is the **SIMULATOR's** class mix, not any
real batch's prevalence. A batch that is 90 % coloc makes every `v_i` too large and the rule too
conservative; a batch that is 1 % coloc makes them too small and **the FDR claim FALSE**. So the
whole rule is re-run at every rung of the pre-registered grid and the curve is a field of the
result. Empirical-Bayes estimation of the prior from the batch is deliberately NOT offered as a
default: it is a second estimated quantity with its own failure modes, and estimating the prior
from the same batch the guarantee is quoted over would make that guarantee CIRCULAR.

The grid is pre-registered, and `P14BatchDecision` pins the curve to it, so passing a different
`pi_grid` fails loudly at construction rather than quietly producing an off-grid curve.

# Returns

A `P14BatchDecision` whose `fdr_scope` is always `:decided_subset_only`. Read `p14_headline(b)`
for the one line that may never be split, `b.pi_sensitivity` for the curve, `b.results` for the
finalised per-pair results, and `p14_named_limits()` for the eight items that travel with any
quoted number.
"""
function decide_coloc(bundle::NamedTuple, pairs::AbstractVector, ood_nulls, qhat::Real;
                      alpha_fdr::Real,
                      prior::NamedTuple = bundle.prior,
                      allow_unchecked_ood::Bool = false,
                      pi_grid = P14_PI_COLOC_GRID,
                      idx_offset::Integer = 0,
                      costes_seed = P14_DEV_SEED,
                      meta::NamedTuple = NamedTuple())
    # --- validation first --------------------------------------------------------------------
    0.0 < alpha_fdr < 1.0 || throw(ArgumentError(
        "decide_coloc: `alpha_fdr` must lie strictly inside (0,1); got $alpha_fdr. It is a " *
        "per-call USER parameter (SC1), never a constant of the pre-registration and never " *
        "derived from the conformal miscoverage level."))
    n_total = length(pairs)
    n_total > 0 || throw(ArgumentError(
        "decide_coloc: the batch is empty. A decided fraction over an empty batch is not a " *
        "number, and the fraction is not optional."))
    for (i, q) in enumerate(pairs)
        q isa NamedTuple || throw(ArgumentError(
            "decide_coloc: item $i of `pairs` must be a NamedTuple carrying at least " *
            "$(P14_PAIR_REQUIRED_KEYS); got $(typeof(q))."))
        for k in P14_PAIR_REQUIRED_KEYS
            haskey(q, k) || throw(ArgumentError(
                "decide_coloc: item $i of `pairs` is missing the required key `$k`."))
        end
        for k in keys(q)
            (k in P14_PAIR_REQUIRED_KEYS || k in P14_PAIR_OPTIONAL_KEYS) || throw(ArgumentError(
                "decide_coloc: item $i of `pairs` carries the unknown key `$k`; the permitted " *
                "keys are $(P14_PAIR_REQUIRED_KEYS) and $(P14_PAIR_OPTIONAL_KEYS). An ignored " *
                "misspelling would silently disable a channel, and a disabled cross-method " *
                "channel reports `disagree_any = false`, which reads as AGREEMENT."))
        end
    end

    # 1. EVERY PAIR IS DECIDED INDIVIDUALLY.
    provisional = P14Result[
        p14_decide_one(bundle, q.Zs, q.Zc, ood_nulls, qhat;
                       lambda = get(q, :lambda, P13_PHASE11_REFERENCE_LAMBDA),
                       idx = get(q, :idx, idx_offset + i),
                       mci_s = get(q, :mci_s, nothing), mci_c = get(q, :mci_c, nothing),
                       prior = prior, allow_unchecked_ood = allow_unchecked_ood,
                       local_map_uncontrolled = get(q, :local_map_uncontrolled, nothing),
                       costes_seed = costes_seed)
        for (i, q) in enumerate(pairs)]

    # The rule inputs that do NOT depend on the prior, read back off the results themselves --
    # which also proves a result carries enough to reproduce the rule that produced it.
    recs = [(lbf = log_bf_vs_random(r), ood = ood_state(r), cm = cross_method(r))
            for r in provisional]

    # 2-3. ABSTAIN FIRST, THEN SORT -- through the same code path the sensitivity curve uses.
    run = _p14_rule_at_prior(recs, qhat, prior, alpha_fdr, allow_unchecked_ood, n_total)

    # THE TWO PATHS MUST AGREE. `p14_decide_one` and `_p14_rule_at_prior` compute the same fusion
    # from the same inputs; a divergence means one of them drifted, and it would be invisible in
    # the output because both produce a perfectly well-formed batch.
    for i in eachindex(provisional)
        @assert (provisional[i].decision === :abstain) == (run.acts[i].action === :abstain) "decide_coloc: item $i is $(provisional[i].decision) on the per-pair path but $(run.acts[i].action) on the batch path; the two computations of the D-05 fusion have drifted apart"
    end

    # 4. BACK TO ORIGINAL BATCH INDICES. `fdr.accepted` holds positions into the DECIDED SUBSET,
    #    never into the batch, so the mapping is explicit and then asserted.
    accepted_original = [run.decided_idx[j] for j in run.fdr.accepted]
    @assert all(i -> provisional[i].decision !== :abstain, accepted_original) "decide_coloc: an ABSTAINED item reached the accepted set; the decided-subset index mapping is off and the wrong items would be reported as discoveries"
    @assert length(accepted_original) == run.fdr.k_star "decide_coloc: the mapped accepted set holds $(length(accepted_original)) items but k_star is $(run.fdr.k_star)"

    # 5. FINALISE. Rebuilt, never mutated -- the type is immutable and that is deliberate.
    acc = Set(accepted_original)
    results = P14Result[_p14_finalise(provisional[i], i in acc) for i in eachindex(provisional)]

    # 6. The REQUIRED prior-sensitivity curve, every rung a full re-run of the whole rule.
    sens = [begin
                sr = _p14_rule_at_prior(recs, qhat, _p14_prior_at(prior, pc), alpha_fdr,
                                        allow_unchecked_ood, n_total)
                (pi_coloc = pc, k_star = sr.fdr.k_star, fdp = sr.fdr.fdp,
                 n_accepted = length(sr.fdr.accepted), n_decided = sr.fdr.n_decided,
                 decided_fraction = sr.fdr.decided_fraction)
            end
            for pc in pi_grid]

    return p14_batch_decision(results;
        alpha_fdr = alpha_fdr,
        fdr = run.fdr,
        n_total = n_total,
        abstain_reason_counts = _p14_reason_counts(results),
        class_prior_used = prior,
        class_prior_source = _p14_prior_source(prior, bundle),
        pi_sensitivity = sens,
        conformal = p14_hedge_diagnostics([r.conformal for r in results], qhat),
        meta = merge(meta, (accepted_batch_indices = accepted_original,
                            decided_batch_indices = run.decided_idx,
                            allow_unchecked_ood = allow_unchecked_ood,
                            idx_offset = Int(idx_offset),
                            ood_channels_wired = (:density,),
                            ood_channels_not_wired = (:pp, :noise),
                            named_limits = p14_named_limits(),
                            provenance = p14_provenance_record(bundle.prov))))
end

"""
    _p14_finalise(r, accepted) -> P14Result

Rebuild one per-pair result with its FINAL decision: `:coloc` when the batch rule accepted it,
`:not_coloc` when the rule decided but did not accept it, and `:abstain` untouched.

**A discovery is defined by the rule, not by the argmax.** The two agree in the overwhelming
majority of cases -- a low null posterior is what makes an item both the argmax coloc call and an
accepted one -- but at a permissive level they can part company, and when they do the RULE is
authoritative, because the controlled quantity is a property of the accepted set. The provisional
call is kept in `meta.provisional_decision` so the disagreement is visible rather than erased.

Rebuilt rather than mutated: `P14Result` is immutable, its constructor is where every invariant of
the type is enforced, and going through it again is what proves the finalised object is still a
legal one.
"""
function _p14_finalise(r::P14Result, accepted::Bool)
    r.decision === :abstain && return r
    final = accepted ? :coloc : :not_coloc
    final === r.decision && r.meta isa NamedTuple && haskey(r.meta, :provisional_decision) &&
        return r
    return P14Result(r.three_way, r.class_posterior, r.null_posterior, r.null_split,
                     final, nothing, r.ood_state, r.conformal, r.cross_method,
                     r.local_map_uncontrolled,
                     merge(r.meta, (provisional_decision = r.decision,
                                    accepted_by_fdr_rule = accepted,
                                    decision_is_provisional_on_the_coloc_side = false)))
end

"""
    _p14_prior_source(prior, bundle) -> Symbol

`:measured_pi_class_masses` when the prior in use IS the measured class mix of the prior the net
was trained under, `:caller_supplied` otherwise.

Compared BY NAME rather than by object identity, so a caller who rebuilt the same three masses in
a different key order is correctly recorded as having used the measured prior; the label is about
which NUMBERS were used, not about which object was passed.
"""
function _p14_prior_source(prior::NamedTuple, bundle::NamedTuple)
    m = bundle.prior
    same = prior.coloc == m.coloc && prior.random == m.random && prior.exclusion == m.exclusion
    return same ? :measured_pi_class_masses : :caller_supplied
end

# --- Script entry point ---------------------------------------------------------------------------
# Including this file does NOTHING. Running it prints what was inherited and nothing else: no
# training, no simulation, no figure, no artifact, and no random stream consumed.
if abspath(PROGRAM_FILE) == @__FILE__
    let b = p14_load_bundle()
        println("spike/p14/decide.jl -- the inherited Phase-13 bundle")
        println("  tau (LOADED, D-07)      = ", b.tau)
        println("  four-way agreed         = ", b.prov.four_way_agreed)
        println("  grid / summary length   = ", b.grid, " / ", p14_summary_length(b))
        println("  class prior (MEASURED)  = ", b.prior)
        println("  carried calibration     = ECE ", round(b.calibration.ece; sigdigits = 4),
                ", AUC ", round(b.calibration.gate.auc; sigdigits = 4),
                ", verdict ", b.calibration.gate.verdict)
        println()
        println("NAMED LIMITS (all eight travel with every number this file produces):")
        for (i, s) in enumerate(p14_named_limits())
            println("  ", i, ". ", s)
        end
    end
end
