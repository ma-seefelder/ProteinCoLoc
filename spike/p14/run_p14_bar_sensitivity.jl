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

# spike/p14/run_p14_bar_sensitivity.jl --- THE VERDICT-SENSITIVITY STRIP. READ-ONLY.
#
# =============================================================================================
# NO BAR WAS CHANGED BY THIS FILE, AND NO BAR CAN BE.
# =============================================================================================
# THIS RUNNER GATES NOTHING. It re-simulates nothing, re-trains nothing, and re-measures
# nothing. It loads the ALREADY-MEASURED Phase-14 statistics out of the artifacts the gated
# runners wrote, and compares each frozen statistic against a SWEEP OF HYPOTHETICAL bar values
# so that a reader can see whether a verdict sits next to a boundary or far from one.
#
# WHY IT EXISTS. `spike/p14/consts.jl` section E states at length that the four SC3 bars and the
# coverage floor are JUDGEMENT CALLS WITH NO DERIVATION, and this project has recorded four
# separate cases of an underived bar measuring something other than what it named. The user, on
# 2026-08-04, ruled the Phase-14 report's presentation of those five bars as follows:
#
#     "Add a verdict-sensitivity strip. Keep every bar frozen and every verdict as-is, but add a
#      read-only table showing what each SC3 verdict would have been across a range of nearby bar
#      values. This changes nothing measured -- it lets a reader see whether a conclusion hinges
#      on the arbitrary number or is robust to it, which is the actual defence against an
#      underived bar. It is presentation, not relaxation."
#
# THAT IS THE ENTIRE MANDATE OF THIS FILE: PRESENTATION, NOT RELAXATION.
#
# THE FROZEN BARS ARE READ FROM `consts.jl` BY NAME AND ARE NEVER RETYPED HERE. The sweep grids
# below are HYPOTHETICAL values with the `P14_BS_` prefix -- a namespace `consts.jl` does not use
# and no gated runner reads -- so a sweep rung can never be mistaken for, or promoted into, a
# pre-registered bar. Each sweep is asserted to CONTAIN its frozen bar, so the frozen value is
# always visible inside its own strip rather than sitting outside the picture.
#
# WHAT A ONE-SIDED BAR'S SENSITIVITY ACTUALLY IS, STATED PLAINLY RATHER THAN DRESSED UP. Every
# SC3 gate has the form `statistic >= bar`, so the value at which the verdict flips IS the
# measured statistic, by construction. The informative content of the strip is therefore the
# DISTANCE between the frozen bar and that flip point -- how much room the verdict had -- and the
# rungs make that distance legible instead of asking a reader to do the subtraction. This is said
# out loud in the artifact's own honesty block so nobody mistakes the strip for new evidence.
#
# THE COVERAGE FLOOR IS DIFFERENT IN KIND, AND IS TREATED DIFFERENTLY. `P14_COVERAGE_FLOOR` is
# not a pass/fail bar: it is the coverage at which the risk-coverage curve is TRUNCATED before
# SC3-b's Spearman is taken. Its sensitivity is therefore a genuine re-derivation of the
# statistic at other truncations, computed from the RAW FULL-RANGE curve that
# `p14_riskcoverage_report.jld2` persists (`curve_coverage`, `curve_selective_risk`). The
# recomputation is PINNED at the frozen floor to the persisted `spearman` and `spearman_n_points`
# before any other rung is reported, so a strip built on a different curve, a different
# truncation convention or a different correlation function fails loudly instead of quietly
# disagreeing with the number it sits beside.
#
# THE SECOND TABLE IS SC3-d's SCOPE, NOT SC3-d's VERDICT. SC3-d FAILED at a margin of -0.129 and
# nothing here changes that. The user's SC3-d ruling (2026-08-04,
# `.planning/phases/14-decision-and-abstention-layer/14-SC3D-RULING.md`) requires the report to
# state, beside the failure, that the measurement was taken on a DENSITY-ONLY detector
# (`ood_channels_wired = (:density,)`) while the shipped detector fuses density with noise. That
# comparison is loaded here from the SHIPPED gate report,
# `artifacts/amended_v2/grid_8/gate_report_8.jld2` key `report[:ood]`, rather than retyped from
# the ruling document -- the ruling is a document, the artifact is the evidence.
#
# READ-ONLY IN THE STRONGEST SENSE. This runner opens `artifacts/` for READING only and writes
# exactly one file, its own gitignored `p14_bar_sensitivity.jld2`. It touches no gated artifact,
# no figure and no constant.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY AND NOT A TEST. Including it does nothing; run it:
#
#     julia --project=spike spike/p14/run_p14_bar_sensitivity.jl
#
# DECOUPLING: spike-local, no dependency added (JLD2 and StatsBase are among the frozen sixteen).

using JLD2
using Printf
import StatsBase      # corspearman, never hand-rolled -- the same function SC3-b was scored with

# --- Guarded includes, in dependency order --------------------------------------------------------
# Only two are needed. This runner deliberately does NOT include `decide.jl`: it loads no net,
# draws no sample and needs no posterior. A runner that pulls the whole inference chain in order to
# read six numbers out of a `.jld2` invites the suspicion that it recomputed one of them.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED)   || include(joinpath(@__DIR__, "provenance.jl"))

if !isdefined(@__MODULE__, :P14_BS_RUNNER_LOADED)

    "The strip's own artifact. Gitignored and regenerable, like every other `spike/p14/*.jld2`."
    const P14_BS_REPORT_PATH = joinpath(@__DIR__, "p14_bar_sensitivity.jld2")

    "The GATED SC3-a/b/c artifact, opened READ-ONLY. Never rewritten by this runner."
    const P14_BS_RC_PATH = joinpath(@__DIR__, "p14_riskcoverage_report.jld2")

    "The GATED SC3-d artifact, opened READ-ONLY. Never rewritten by this runner."
    const P14_BS_OOD_PATH = joinpath(@__DIR__, "p14_ood_arm_report.jld2")

    """
    The SHIPPED grid-8 gate report, opened READ-ONLY.

    This is the v2.0 ship-gate artifact, not a Phase-14 output. It is read for exactly one
    purpose: `report[:ood]` carries `density_auc`, `noise_auc` and `fused_auc` per
    misspecification family, which is the evidence that SC3-d's failing family is the family the
    UNWIRED channel exists for. Nothing is written to `artifacts/`.
    """
    const P14_BS_GATE8_PATH =
        normpath(joinpath(@__DIR__, "..", "..", "artifacts", "amended_v2", "grid_8",
                          "gate_report_8.jld2"))

    # =========================================================================================
    # THE HYPOTHETICAL SWEEP GRIDS -- NOT BARS, AND NAMED SO THEY CANNOT BE MISTAKEN FOR BARS
    # =========================================================================================
    # `P14_BS_` is a namespace `consts.jl` never uses and no gated runner reads. Each grid is
    # asserted below to CONTAIN the frozen bar it sweeps around, so the frozen value always
    # appears inside its own strip. None of these values gates anything, here or anywhere.
    const P14_BS_SKILL_SWEEP    = (0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.91, 0.95)
    const P14_BS_SPEARMAN_SWEEP = (0.80, 0.85, 0.90, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00)
    const P14_BS_AUC_SWEEP      = (0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.91, 0.92, 0.95, 0.99)
    const P14_BS_OOD_SWEEP      = (0.00, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)

    """
    The coverage truncations SC3-b's Spearman is RE-DERIVED at.

    Not a bar and not a candidate bar: a set of alternative truncation points for the raw curve,
    spanning from a floor low enough to include the one-item-denominator noise the frozen floor
    exists to exclude, up to a floor that keeps only the right half of the curve. The frozen
    `P14_COVERAGE_FLOOR` is asserted to be a member so the reported statistic appears inside its
    own sweep.
    """
    const P14_BS_COVERAGE_SWEEP = (0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70)

    """
    The tolerance at which the RE-DERIVED Spearman at the frozen floor must equal the persisted one.

    A reproduction check, never a threshold on a result. If the strip's truncation-and-correlate
    path cannot reproduce the number `run_p14_riskcoverage.jl` persisted, then every other rung of
    the coverage sweep is being produced by a different computation from the one SC3-b was scored
    with, and the strip would be quietly comparing two things.
    """
    const P14_BS_PIN_TOL = 1e-12

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else
    # can make this body skip.
    const P14_BS_RUNNER_LOADED = true
end

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    _p14_bs_save_report(path, required; kwargs...) -> String

Atomically persist an artifact: write `path * ".tmp"`, REOPEN it read-only and integrity-check
every key in `required`, then `mv(...; force = true)`.

The `_p14_rc_save_report` idiom (`spike/p14/run_p14_riskcoverage.jl:198`), carried per runner
exactly as every Phase-13 and Phase-14 runner carries it, and deliberately NOT imported from a
sibling runner: importing would mean including that runner and inheriting its `main`.
"""
function _p14_bs_save_report(path, required; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        for k in required
            @assert haskey(f, k) "_p14_bs_save_report: integrity check failed -- $tmp is missing the required key `$k`"
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    _p14_bs_load(path, what) -> Dict

Load a `.jld2` READ-ONLY, failing loudly and by name if it is absent.

An absent artifact is never silently substituted, and never regenerated by this runner: the
strip's whole claim is that it reports numbers someone else measured, so producing one itself
would destroy the only property it has.
"""
function _p14_bs_load(path::AbstractString, what::AbstractString)
    isfile(path) || error("""
        run_p14_bar_sensitivity: $what is ABSENT at
            $path

        This runner does NOT produce it and will not substitute for it. Re-run the gated runner
        that writes it, then run this strip again.""")
    return JLD2.load(path)
end

"""
    _p14_bs_verdict(measured, bar) -> Symbol

The verdict a one-sided `statistic >= bar` gate would return, and nothing else.

Written as one function used by every row of the strip -- including the row at the FROZEN bar,
which is then pinned to the verdict the gated artifact itself persisted. Two spellings of "did it
pass" is how a strip ends up disagreeing with the table above it.
"""
_p14_bs_verdict(measured::Real, bar::Real) = measured >= bar ? :MET : :NOT_MET

"""
    _p14_bs_strip(measured, frozen_bar, sweep) -> NamedTuple

One bar's sensitivity strip: the verdict at every hypothetical rung, the frozen rung marked, and
the flip point.

THE FLIP POINT IS THE MEASURED STATISTIC, BY CONSTRUCTION, because the gate is one-sided. That is
not a discovery and is not presented as one; what the strip carries is the DISTANCE from the
frozen bar to it -- `measured - frozen_bar` -- which is the quantity a reader asking "does this
conclusion hinge on the arbitrary number?" actually needs.
"""
function _p14_bs_strip(measured::Real, frozen_bar::Real, sweep)
    bars     = collect(Float64, sweep)
    verdicts = [_p14_bs_verdict(measured, b) for b in bars]
    return (measured = Float64(measured),
            frozen_bar = Float64(frozen_bar),
            bars = bars,
            verdicts = verdicts,
            frozen_index = findfirst(==(Float64(frozen_bar)), bars),
            frozen_verdict = _p14_bs_verdict(measured, frozen_bar),
            flip_point = Float64(measured),
            distance_to_flip = Float64(measured) - Float64(frozen_bar))
end

"""
    _p14_bs_coverage_spearman(cov, risk, floor) -> NamedTuple

SC3-b's statistic RE-DERIVED at an alternative truncation of the raw full-range curve.

Identical in construction to `run_p14_riskcoverage.jl:747-752` -- select the curve points with
`coverage >= floor` and `<= 1.0`, then `StatsBase.corspearman` -- so that at the frozen floor it
must reproduce the persisted number exactly. Fewer than two points in range is reported as `NaN`
with the count, never as a correlation.
"""
function _p14_bs_coverage_spearman(cov::AbstractVector{<:Real},
                                   risk::AbstractVector{<:Real},
                                   floor::Real)
    in_range = [(c >= floor) & (c <= 1.0) for c in cov]
    n = count(in_range)
    n >= 2 || return (spearman = NaN, n_points = n)
    return (spearman = StatsBase.corspearman(cov[in_range], risk[in_range]), n_points = n)
end

"""
    _p14_bs_claims() -> NamedTuple

The honesty block, written INTO the artifact as string keys so it travels WITH the strip and
cannot be dropped by a report writer.

The strip is the single most misreadable object this phase produces: a table of "what the verdict
would have been" is one careless sentence away from being read as a table of verdicts. Every
sentence below exists to make that misreading impossible from the artifact alone.
"""
function _p14_bs_claims()
    return (
        strip_note =
            "NO BAR WAS CHANGED. This artifact GATES NOTHING. The only bar values that gate " *
            "anything in Phase 14 are the frozen ones in `spike/p14/consts.jl`, which is " *
            "byte-unchanged against HEAD, and every reported Phase-14 verdict is the verdict at " *
            "those frozen values. The rungs here are HYPOTHETICAL: they exist so a reader can " *
            "see whether a conclusion hinges on an underived number or is robust to it. This is " *
            "PRESENTATION, NOT RELAXATION.",

        mandate_note =
            "USER RULING, 2026-08-04, on the presentation of the five judgement-call bars: `Add " *
            "a verdict-sensitivity strip. Keep every bar frozen and every verdict as-is, but add " *
            "a read-only table showing what each SC3 verdict would have been across a range of " *
            "nearby bar values. This changes nothing measured -- it lets a reader see whether a " *
            "conclusion hinges on the arbitrary number or is robust to it, which is the actual " *
            "defence against an underived bar. It is presentation, not relaxation.'",

        one_sided_note =
            "WHAT A ONE-SIDED BAR'S SENSITIVITY IS, SAID PLAINLY. Every SC3 gate has the form " *
            "`statistic >= bar`, so the bar value at which the verdict flips IS the measured " *
            "statistic, by construction. The strip's informative content is therefore the " *
            "DISTANCE from the frozen bar to that flip point, i.e. the margin the verdict had. " *
            "It is not new evidence and must not be quoted as any.",

        coverage_floor_note =
            "THE COVERAGE FLOOR IS DIFFERENT IN KIND FROM THE OTHER FOUR. It is not a pass/fail " *
            "bar: it is the coverage at which the raw risk-coverage curve is TRUNCATED before " *
            "SC3-b's Spearman is computed. Its sensitivity is therefore a genuine re-derivation " *
            "of the statistic at other truncations, taken from the RAW FULL-RANGE curve " *
            "persisted in `p14_riskcoverage_report.jld2`. The re-derivation is PINNED at the " *
            "frozen floor to the persisted `spearman` and `spearman_n_points` before any other " *
            "rung is reported. The reason a floor exists at all is unchanged: at coverage k/n a " *
            "single hard item moves the selective risk by 1/k, so at the left edge the " *
            "denominator is one item and the risk is exactly 0 or 1.",

        sc3d_scope_note =
            "SC3-d FAILED, at a margin of -0.129 against the frozen floor, and NOTHING IN THIS " *
            "ARTIFACT SOFTENS THAT. The `density_auc` / `noise_auc` / `fused_auc` table read " *
            "here from the SHIPPED gate report is a statement of SCOPE, not of verdict: the " *
            "Phase-14 decision layer wired ONE of the shipped flag's three channels " *
            "(`ood_channels_wired = (:density,)`), and the channel it did not wire is the one " *
            "whose family SC3-d failed on. THE FUSED DETECTOR'S SC3-d MARGIN IS UNMEASURED -- " *
            "Phase 14 never ran the fused detector through `decide_coloc` -- and no reader may " *
            "infer from this table that SC3-d would pass with all channels wired.",

        not_recomputed_note =
            "NOTHING WAS RE-SIMULATED, RE-TRAINED OR RE-MEASURED. `skill`, `spearman`, " *
            "`auc_hard` and the SC3-d margin are LOADED from the gated artifacts and compared " *
            "against hypothetical thresholds. The single quantity re-derived here is SC3-b's " *
            "Spearman at alternative TRUNCATIONS of the already-persisted raw curve, and it is " *
            "pinned to the persisted value at the frozen truncation.",

        amendment = P14_AMENDMENT_NOTICE,
    )
end

# =============================================================================================
# The read-only run
# =============================================================================================

"""
    main(; report_path = P14_BS_REPORT_PATH, verbose = true) -> NamedTuple

The verdict-sensitivity strip: five rows, all read-only, none of them a gate.
"""
function main(; report_path = P14_BS_REPORT_PATH, verbose::Bool = true)

    t_start = time()

    # --- 0. THE DECOUPLING PROOF AND THE PRE-REGISTRATION PROOF, AT RUN TIME ------------------
    # The shared lane guard covers `src/` and the frozen spike environment. A THIRD check is
    # added here that no other runner makes and that this runner in particular owes the reader:
    # `spike/p14/consts.jl` -- the file holding every bar this strip sweeps around -- must be
    # byte-unchanged against HEAD WHILE the strip is being built. A strip built beside an edited
    # pre-registration would be the exact move it exists to make unnecessary.
    lane = p14_assert_lane_clean(; verbose = verbose)
    consts_clean = if Sys.which("git") === nothing
        missing
    else
        ok = success(Cmd(`git diff --quiet HEAD -- spike/p14/consts.jl`; dir = P14_REPO_ROOT))
        @assert ok "PRE-REGISTRATION BREACH: `git diff --quiet HEAD -- spike/p14/consts.jl` FAILED. The frozen bars moved under the run; a sensitivity strip built beside an edited pre-registration is worthless"
        ok
    end

    # --- 1. THE SWEEPS CONTAIN THEIR FROZEN BARS ---------------------------------------------
    # Asserted rather than eyeballed: a sweep that does not contain its frozen bar would show a
    # reader every value EXCEPT the one that actually gates, which is the one presentation
    # failure this strip cannot afford.
    @assert Float64(P14_SKILL_FLOOR)      in P14_BS_SKILL_SWEEP    "the skill sweep does not contain the frozen P14_SKILL_FLOOR"
    @assert Float64(P14_SPEARMAN_FLOOR)   in P14_BS_SPEARMAN_SWEEP "the spearman sweep does not contain the frozen P14_SPEARMAN_FLOOR"
    @assert Float64(P14_AUC_HARD_FLOOR)   in P14_BS_AUC_SWEEP      "the AUC sweep does not contain the frozen P14_AUC_HARD_FLOOR"
    @assert Float64(P14_OOD_MARGIN_FLOOR) in P14_BS_OOD_SWEEP      "the OOD-margin sweep does not contain the frozen P14_OOD_MARGIN_FLOOR"
    @assert Float64(P14_COVERAGE_FLOOR)   in P14_BS_COVERAGE_SWEEP "the coverage sweep does not contain the frozen P14_COVERAGE_FLOOR"

    # --- 2. LOAD THE ALREADY-MEASURED STATISTICS ---------------------------------------------
    rc  = _p14_bs_load(P14_BS_RC_PATH,  "the GATED SC3-a/b/c artifact")
    ood = _p14_bs_load(P14_BS_OOD_PATH, "the GATED SC3-d artifact")

    skill    = Float64(rc["skill"])
    spearman = Float64(rc["spearman"])
    auc_hard = Float64(rc["auc_hard"])
    margin   = Float64(ood["margin"])

    cov  = Vector{Float64}(rc["curve_coverage"])
    risk = Vector{Float64}(rc["curve_selective_risk"])
    @assert rc["curve_is_raw_full_range"] == true "run_p14_bar_sensitivity: the persisted curve is not flagged raw-full-range; the coverage sweep must be taken on the untrimmed curve"

    if verbose
        println("="^78)
        println("PHASE-14 VERDICT-SENSITIVITY STRIP -- READ-ONLY, GATES NOTHING")
        println("="^78)
        println("NO BAR WAS CHANGED. The frozen bars in spike/p14/consts.jl are the ONLY values")
        println("that gate anything, and consts.jl is byte-unchanged against HEAD ($consts_clean).")
        println("This strip is PRESENTATION, NOT RELAXATION -- the user's ruling of 2026-08-04.")
        println("FROZEN BARS, read from consts.jl BY NAME (never retyped here):")
        println("  P14_SKILL_FLOOR      = $P14_SKILL_FLOOR")
        println("  P14_SPEARMAN_FLOOR   = $P14_SPEARMAN_FLOOR")
        println("  P14_COVERAGE_FLOOR   = $P14_COVERAGE_FLOOR   (a truncation, NOT a pass/fail bar)")
        println("  P14_AUC_HARD_FLOOR   = $P14_AUC_HARD_FLOOR")
        println("  P14_OOD_MARGIN_FLOOR = $P14_OOD_MARGIN_FLOOR")
        println("  P14_JUDGEMENT_CALL_BARS = $P14_JUDGEMENT_CALL_BARS")
        println("MEASURED STATISTICS, LOADED (not recomputed):")
        println("  skill    = $skill")
        println("  spearman = $spearman")
        println("  auc_hard = $auc_hard")
        println("  SC3-d margin = $margin   (weakest family: $(ood["weakest_family"]))")
        println("="^78)
    end

    # --- 3. THE FOUR ONE-SIDED STRIPS ---------------------------------------------------------
    # The verdict at the FROZEN rung is pinned to the verdict the gated artifact itself persisted.
    # This is the check that makes the strip trustworthy: if this runner's verdict function
    # disagreed with the gated runner's at the one value that matters, every hypothetical rung
    # beside it would be produced by a different rule from the one Phase 14 was scored under.
    strip_skill    = _p14_bs_strip(skill,    P14_SKILL_FLOOR,      P14_BS_SKILL_SWEEP)
    strip_spearman = _p14_bs_strip(spearman, P14_SPEARMAN_FLOOR,   P14_BS_SPEARMAN_SWEEP)
    strip_auc      = _p14_bs_strip(auc_hard, P14_AUC_HARD_FLOOR,   P14_BS_AUC_SWEEP)
    strip_ood      = _p14_bs_strip(margin,   P14_OOD_MARGIN_FLOOR, P14_BS_OOD_SWEEP)

    @assert (strip_skill.frozen_verdict === :MET)    == Bool(rc["sc3a_met"])  "run_p14_bar_sensitivity: the strip's verdict at the frozen skill floor disagrees with the persisted sc3a_met"
    @assert (strip_spearman.frozen_verdict === :MET) == Bool(rc["sc3b_met"])  "run_p14_bar_sensitivity: the strip's verdict at the frozen spearman floor disagrees with the persisted sc3b_met"
    @assert (strip_auc.frozen_verdict === :MET)      == Bool(rc["sc3c_met"])  "run_p14_bar_sensitivity: the strip's verdict at the frozen AUC floor disagrees with the persisted sc3c_met"
    @assert (strip_ood.frozen_verdict === :MET)      == Bool(ood["sc3d_met"]) "run_p14_bar_sensitivity: the strip's verdict at the frozen OOD-margin floor disagrees with the persisted sc3d_met"

    # --- 4. THE COVERAGE-FLOOR SWEEP, PINNED AT THE FROZEN TRUNCATION -------------------------
    pin = _p14_bs_coverage_spearman(cov, risk, P14_COVERAGE_FLOOR)
    @assert abs(pin.spearman - spearman) < P14_BS_PIN_TOL "run_p14_bar_sensitivity: re-deriving SC3-b's Spearman at the FROZEN floor gives $(pin.spearman) against the persisted $spearman; the sweep is not the computation SC3-b was scored with"
    @assert pin.n_points == Int(rc["spearman_n_points"]) "run_p14_bar_sensitivity: the frozen truncation selects $(pin.n_points) curve points against the persisted $(rc["spearman_n_points"])"

    cov_floors     = collect(Float64, P14_BS_COVERAGE_SWEEP)
    cov_spearmans  = Float64[]
    cov_npoints    = Int[]
    cov_verdicts   = Symbol[]
    for f in cov_floors
        r = _p14_bs_coverage_spearman(cov, risk, f)
        push!(cov_spearmans, r.spearman)
        push!(cov_npoints, r.n_points)
        push!(cov_verdicts, isnan(r.spearman) ? :UNDEFINED : _p14_bs_verdict(r.spearman, P14_SPEARMAN_FLOOR))
    end

    # --- 5. SC3-d's SCOPE, LOADED FROM THE SHIPPED GATE REPORT --------------------------------
    # A statement of scope, never of verdict. See `sc3d_scope_note`.
    gate8 = _p14_bs_load(P14_BS_GATE8_PATH, "the SHIPPED grid-8 gate report")
    g_ood = gate8["report"][:ood]
    fams  = collect(Symbol, g_ood.families)
    density_auc = Dict(f => collect(Float64, g_ood.density_auc[f]) for f in fams)
    noise_auc   = Dict(f => collect(Float64, g_ood.noise_auc[f])   for f in fams)
    fused_auc   = Dict(f => collect(Float64, g_ood.fused_auc[f])   for f in fams)

    # --- 6. PRINT THE STRIP -------------------------------------------------------------------
    if verbose
        println()
        println("### Verdict sensitivity -- what each verdict would have been at a different bar")
        println("### (NO BAR WAS CHANGED. Hypothetical rungs; the frozen bar is marked *.)")
        println()
        for (name, sym, st) in (("SC3-a selective skill", :P14_SKILL_FLOOR, strip_skill),
                                ("SC3-b Spearman",        :P14_SPEARMAN_FLOOR, strip_spearman),
                                ("SC3-c AUC_hard",        :P14_AUC_HARD_FLOOR, strip_auc),
                                ("SC3-d abstain margin",  :P14_OOD_MARGIN_FLOOR, strip_ood))
            @printf("| %-22s | measured %-12s |", name, @sprintf("%.4f", st.measured))
            for (b, v) in zip(st.bars, st.verdicts)
                mark = (b == st.frozen_bar) ? "*" : " "
                @printf(" %s%.2f:%s |", mark, b, v === :MET ? "MET" : "NOT")
            end
            @printf("  frozen %s = %.2f -> %s, flips at %.4f (distance %+.4f)\n",
                    String(sym), st.frozen_bar, String(st.frozen_verdict),
                    st.flip_point, st.distance_to_flip)
        end

        println()
        println("### P14_COVERAGE_FLOOR -- different in kind: a TRUNCATION, not a pass/fail bar")
        println("### SC3-b's Spearman RE-DERIVED from the persisted raw full-range curve.")
        println()
        println("| truncation | curve points | Spearman | vs the frozen spearman floor |")
        for i in eachindex(cov_floors)
            mark = (cov_floors[i] == Float64(P14_COVERAGE_FLOOR)) ? "*" : " "
            @printf("| %s%.2f | %5d | %s | %s |\n", mark, cov_floors[i], cov_npoints[i],
                    isnan(cov_spearmans[i]) ? "undefined" : @sprintf("%.4f", cov_spearmans[i]),
                    String(cov_verdicts[i]))
        end

        println()
        println("### SC3-d SCOPE (NOT verdict) -- artifacts/amended_v2/grid_8/gate_report_8.jld2, report[:ood]")
        println("### SC3-d FAILED at margin $margin and nothing here changes that.")
        println()
        println("| family | density_auc (what P14 wired) | noise_auc | fused_auc (SHIPPED) |")
        for f in fams
            @printf("| %-11s | %s | %s | %s |\n", String(f),
                    join(density_auc[f], " / "), join(noise_auc[f], " / "),
                    join(fused_auc[f], " / "))
        end
        println("id_threshold = $(g_ood.id_threshold)   fused_threshold = $(g_ood.fused_threshold)   " *
                "levels = $(g_ood.levels)   n = $(g_ood.n)")
        println("P14 wired: ", ood["channels_wired"], "   NOT wired: ", ood["channels_not_wired"])
        println()
    end

    # --- 7. PERSIST ---------------------------------------------------------------------------
    elapsed = time() - t_start
    _p14_bs_save_report(report_path,
        ("strip_note", "mandate_note", "one_sided_note", "coverage_floor_note",
         "sc3d_scope_note", "not_recomputed_note", "amendment",
         "skill", "spearman", "auc_hard", "sc3d_margin",
         "skill_strip", "spearman_strip", "auc_hard_strip", "ood_margin_strip",
         "coverage_floors", "coverage_spearmans", "coverage_n_points", "coverage_verdicts",
         "coverage_pin_ok", "gate8_density_auc", "gate8_noise_auc", "gate8_fused_auc",
         "gates_nothing", "no_bar_changed", "provenance");
        _p14_bs_claims()...,
        # --- what this artifact is, machine-readable ---
        gates_nothing = true,
        no_bar_changed = true,
        recomputed_anything_gated = false,
        consts_byte_unchanged_at_run = consts_clean,
        lane_checked = lane.checked,
        # --- the LOADED measured statistics, and where they came from ---
        skill = skill,
        spearman = spearman,
        auc_hard = auc_hard,
        sc3d_margin = margin,
        sc3d_weakest_family = ood["weakest_family"],
        source_riskcoverage = P14_BS_RC_PATH,
        source_ood_arm = P14_BS_OOD_PATH,
        source_shipped_gate = P14_BS_GATE8_PATH,
        # --- the four one-sided strips ---
        skill_strip = strip_skill,
        spearman_strip = strip_spearman,
        auc_hard_strip = strip_auc,
        ood_margin_strip = strip_ood,
        # --- the coverage-floor re-derivation ---
        coverage_floors = cov_floors,
        coverage_spearmans = cov_spearmans,
        coverage_n_points = cov_npoints,
        coverage_verdicts = cov_verdicts,
        coverage_pin_ok = true,
        coverage_pin_spearman = pin.spearman,
        coverage_pin_n_points = pin.n_points,
        coverage_pin_tolerance = P14_BS_PIN_TOL,
        # --- SC3-d's scope, from the SHIPPED artifact ---
        gate8_families = fams,
        gate8_density_auc = density_auc,
        gate8_noise_auc = noise_auc,
        gate8_fused_auc = fused_auc,
        gate8_id_threshold = g_ood.id_threshold,
        gate8_fused_threshold = g_ood.fused_threshold,
        gate8_levels = g_ood.levels,
        gate8_n = g_ood.n,
        p14_channels_wired = ood["channels_wired"],
        p14_channels_not_wired = ood["channels_not_wired"],
        fused_sc3d_margin_measured = false,
        # --- the frozen bars, written INTO the artifact so the strip carries them ---
        P14_SKILL_FLOOR = P14_SKILL_FLOOR,
        P14_SPEARMAN_FLOOR = P14_SPEARMAN_FLOOR,
        P14_COVERAGE_FLOOR = P14_COVERAGE_FLOOR,
        P14_AUC_HARD_FLOOR = P14_AUC_HARD_FLOOR,
        P14_OOD_MARGIN_FLOOR = P14_OOD_MARGIN_FLOOR,
        P14_JUDGEMENT_CALL_BARS = P14_JUDGEMENT_CALL_BARS,
        # --- the hypothetical sweeps, labelled as such ---
        sweep_skill = collect(Float64, P14_BS_SKILL_SWEEP),
        sweep_spearman = collect(Float64, P14_BS_SPEARMAN_SWEEP),
        sweep_auc_hard = collect(Float64, P14_BS_AUC_SWEEP),
        sweep_ood_margin = collect(Float64, P14_BS_OOD_SWEEP),
        sweep_coverage = collect(Float64, P14_BS_COVERAGE_SWEEP),
        sweeps_are_hypothetical = true,
        elapsed_s = elapsed,
        provenance = p14_provenance_record(p14_load_tau();
            extra = (runner = "run_p14_bar_sensitivity.jl", gated = false, elapsed_s = elapsed)),
    )

    verbose && println("artifact -> $report_path   ($(round(elapsed; digits = 2)) s)")

    return (skill_strip = strip_skill, spearman_strip = strip_spearman,
            auc_hard_strip = strip_auc, ood_margin_strip = strip_ood,
            coverage_floors = cov_floors, coverage_spearmans = cov_spearmans,
            report_path = report_path)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so including this file can NEVER trigger the run. Run it deliberately:
#
#     julia --project=spike spike/p14/run_p14_bar_sensitivity.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
