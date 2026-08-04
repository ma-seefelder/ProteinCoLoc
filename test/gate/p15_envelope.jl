#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/p15_envelope.jl --- the Phase-15 DECISION RULES, as pure functions.
#
# This file contains NO inference, NO forward model and NO net. It is arithmetic over records
# that some other file produced. That separation is the point: the break criterion and the SC2
# verdict are the two most easily-corrupted decisions in Phase 15, so they are written down and
# fixture-proven HERE, BEFORE any rung exists to argue with them
# (`test/gate/test_p15_break.jl`).
#
# ============================ WHY ECE IS THE BREAK STATISTIC (D-04) ===========================
# ECE is an EFFECT SIZE. It measures how far empirical central-interval coverage sits from
# nominal, in coverage units, and it does NOT grow teeth as M rises. The coverage band and the
# KS / chi-square p-values are computed at EVERY rung and gate NOTHING. Both alternatives were
# rejected for reasons already paid for in this project:
#
#   * A SIGNIFICANCE-BASED break measures `SBC_M`, not the tool. At M = 2000 the chi-square test
#     rejects on drift far below what matters — the documented over-powered-chi-square artifact
#     behind part of the Phase-7 amended-gate failure. Reporting the p-values but not gating on
#     them keeps them auditable instead of hidden.
#   * A PRE-REGISTERED COVERAGE BAND is structurally the criterion the 2026-07-31 Phase-12 ruling
#     named as "a criterion whose pass/fail tracks something other than what it was meant to
#     measure". Same shape, same failure. So the band is reported and never gates.
#
# ========================= WHY THE ANCHOR IS RUNG 0, NOT THE SHIPPED ECE (D-04a) ==============
# The break threshold is `ECE(rung 0, same M, same stream) + P15_ECE_MARGIN`. It is NEVER the
# shipped gate report's ECE. This is arithmetic, not preference: the shipped Delta-rho ECE is
# 0.00429 at M = 2000, while the null mean of a PERFECTLY CALIBRATED estimator at that same M is
# 0.00751. Anchoring at 0.00429 would therefore declare a break on roughly 90 % of perfectly
# calibrated rungs — which is D-04's own failure mode reappearing inside the anchor. Anchoring at
# rung 0 at the sweep's own M makes the null's noise floor COMMON-MODE between anchor and rung,
# so it cancels instead of accumulating.
#
# ================================ MCE IS FORBIDDEN IN THIS PHASE ==============================
# `_bin_calibration` (`sbc.jl:60-108`) walks all 50 reliability bins INCLUDING the 31 empty ones,
# where the predicted rate is the bin midpoint and the observed rate is 0. The maximum gap is
# therefore always attained on an empty bin, and the statistic is pinned at exactly 0.99 on all
# eight columns of the shipped report REGARDLESS OF THE NET. It is an artifact of the binning,
# not a measurement of anything. No Phase-15 code may read that field — `p15_report_only` below
# does not even pass it through, and `test_p15_break.jl` enforces the ban with a source guard
# that has been shown to fire (`P15_MCE_FORBIDDEN`, D-04a).
#
# ============================== THE DETECTOR IS READ-ONLY HERE (D-14) =========================
# This file CONSUMES `gate_ood_roc`'s output (`fire_rate[fam][lvl]` and the MEASURED
# `id_fire_rate`) and nothing else. It never reads, writes or derives a detector parameter.
# Lowering the detector's in-distribution quantile until a LATE axis passes would be tuning a
# detector until it passes its own gate; that response is forbidden by D-14 and is enforced by a
# second source guard in `test_p15_break.jl`.
#
# LOADING: guarded includes, mirroring the rest of `test/gate/`. `p15_consts.jl` supplies the
# frozen Tier-1 bars; `sbc.jl` supplies the PRODUCTION statistics (`sbc_calibration`,
# `sbc_coverage`, `SBC_PARAM_LABELS`) this file is defined against.

isdefined(@__MODULE__, :SBC_M)    || include(joinpath(@__DIR__, "p15_consts.jl"))
isdefined(@__MODULE__, :sbc_ranks) || include(joinpath(@__DIR__, "sbc.jl"))

import JLD2
import Random123: Philox4x
import Statistics

"The repository root, resolved from THIS file rather than from the process working directory."
const P15_ENVELOPE_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# =============================================================================================
# The ECE identity — the NULL-MODEL TWIN of the production statistic, never its replacement
# =============================================================================================

"""
    p15_ece(ranks_col, L; levels = P15_COVERAGE_LEVELS) -> Float64

The closed form of the Phase-15 break statistic for ONE rank column:

    ECE = mean over a in {0.05, 0.10, ..., 0.95} of |a - c(a)|,
    c(a) = mean( |(rank + 0.5)/(L + 1) - 0.5| <= a/2 )

i.e. the mean absolute deviation of empirical central-interval coverage from nominal, over the 19
nominal levels. Derived from `sbc.jl:302-315` together with `_bin_calibration` (`sbc.jl:60-108`):
`sbc_calibration` hands 19*M (predicted, positive) pairs to 50 equal-width bins, each nominal
level lands in its OWN bin holding exactly M samples, and the 31 empty bins carry weight 0. The
statistic therefore does not depend on `SBC_BINS` at all.

THIS FUNCTION HAS EXACTLY TWO ROLES.
 1. It drives `p15_ece_null` — a null model that has no net, no forward model and no data, so it
    cannot use the production path.
 2. It is cross-checked against `sbc_calibration(...).ece` in `test_p15_break.jl`.

THE PRODUCTION PATH IS `sbc_calibration(ranks; L, n_bins).ece`, AND THIS FUNCTION MUST NEVER
REPLACE IT IN THE SWEEP. A second implementation on the reported path would make the envelope
incomparable to the shipped gate report, which is the one comparison the envelope exists to make.
"""
function p15_ece(ranks_col::AbstractVector{<:Integer}, L::Integer;
                 levels = P15_COVERAGE_LEVELS)
    n = length(ranks_col)
    n > 0 || throw(ArgumentError("p15_ece: empty rank vector"))
    L > 0 || throw(ArgumentError("p15_ece: L must be positive, got $L"))
    # The PIT transform and the central-interval distance, formed exactly as `sbc_coverage` and
    # `sbc_calibration` form them (elementwise division by `L + 1`, not by a reciprocal).
    u = (ranks_col .+ 0.5) ./ (L + 1)
    d = abs.(u .- 0.5)
    acc  = 0.0
    nlev = 0
    for alpha in levels
        half = alpha / 2
        hits = 0
        @inbounds for x in d
            x <= half && (hits += 1)
        end
        acc  += abs(alpha - hits / n)
        nlev += 1
    end
    nlev > 0 || throw(ArgumentError("p15_ece: empty level grid"))
    return acc / nlev
end

"""
    p15_ece_null(M, L; B, q, seed, rng) -> (mean = ..., q = ...)

The Monte-Carlo null distribution of `p15_ece` under PERFECT CALIBRATION: draw `B` independent
vectors of `M` uniform ranks on `0:L`, evaluate the identity on each, and return the mean and the
`q` quantile.

THIS IS THE REGENERATOR FOR THE TIER-1 MARGIN. `P15_ECE_MARGIN` is `P15_ECE_NULL_Q95 -
P15_ECE_NULL_MEAN`, and both terms come from this function alone — no net, no forward model, no
outcome. The bar can therefore be re-derived from the recorded test design (`P15_ECE_NULL_B`,
`P15_ECE_NULL_SEED`, `P15_SWEEP_M`, `SBC_L`) by anyone, at any time, with the artifacts tree
deleted. That is what makes it a bar rather than a preference.

ECE HAS A POSITIVE EXPECTATION EVEN UNDER PERFECT CALIBRATION — it is a mean of 19 absolute
deviations of binomial proportions — with the closed form `E[ECE_0] ~= 0.336/sqrt(M)`. That is
the whole reason the break threshold needs a measured anchor and a null-derived margin instead of
a number someone liked.

THE DEFAULT STREAM IS `Philox4x(UInt64, (seed, 0))`, because that is the design
`p15_consts.jl:628` RECORDS as having produced `P15_ECE_NULL_MEAN` / `P15_ECE_NULL_Q95`. A
regenerator that used a different RNG family would be a second estimate of the same quantity, not
a regeneration of the pinned one.
"""
function p15_ece_null(M::Integer, L::Integer;
                      B::Integer = P15_ECE_NULL_B,
                      q::Real    = 0.95,
                      seed::Integer = P15_ECE_NULL_SEED,
                      rng = Philox4x(UInt64, (UInt64(seed), UInt64(0))))
    B > 0 || throw(ArgumentError("p15_ece_null: B must be positive, got $B"))
    v = Vector{Float64}(undef, B)
    for b in 1:B
        v[b] = p15_ece(rand(rng, 0:L, M), L)
    end
    return (mean = Statistics.mean(v), q = Statistics.quantile(v, q))
end

# =============================================================================================
# The shipped anchor — REPORTING CONTEXT, and explicitly NOT the break anchor
# =============================================================================================

"""
    p15_read_shipped_anchor(path) -> NamedTuple or nothing

Read the shipped gate report and return its per-column ECEs as a `Dict` from parameter label to
ECE, together with the `M`, the `L` and the seed the report was produced at:

    (path = ..., ece = Dict(label => ece), M = ..., L = ..., seed = ...)

Returns `nothing` — never an error, never a fabricated value — when the report is absent. That
matches the honest `:not_trained` behaviour of `run_gate.jl:16-19`, and it matters here because
`artifacts/` is git-ignored, so the report is absent in a fresh clone and in every worktree.

**THIS IS REPORTING CONTEXT. IT IS NOT THE BREAK ANCHOR (D-04a).** The shipped values are
rho_true 0.03713 and Delta-rho 0.00429, both measured at M = 2000. The null mean of a PERFECTLY
CALIBRATED estimator at M = 2000 is 0.00751 — so the shipped Delta-rho value sits BELOW its own
null mean, and a threshold anchored there would fire on roughly 90 % of perfectly calibrated
rungs. The break anchor is `ECE(rung 0)` measured at `P15_SWEEP_M` on the Phase-15 stream and
appended to `p15_consts.jl` under the `:P15_ECE_ANCHOR_MEASURED` sentinel.
"""
function p15_read_shipped_anchor(path::AbstractString =
        joinpath(P15_ENVELOPE_REPO_ROOT, P15_ARTIFACTS_ROOT, "grid_8", "gate_report_8.jld2"))
    isfile(path) || return nothing
    rep = JLD2.jldopen(path, "r") do f
        haskey(f, "report") ? f["report"] : nothing
    end
    rep === nothing && return nothing
    return (path  = String(path),
            ece   = Dict{String,Float64}(String(p.label) => float(p.ece)
                                         for p in rep.sbc.per_param),
            M     = rep.sbc.M,
            L     = rep.sbc.L,
            seed  = rep.seed)
end

# =============================================================================================
# The break reduction — columns 1 and 8 ONLY (D-05)
# =============================================================================================

"""
    p15_break_threshold_for(ece_rung0) -> Float64

Thin delegate to the Tier-1 `p15_break_threshold` (`p15_consts.jl:674`), present so that the
sweep runner and the report writer never inline `rung0 + margin` arithmetic of their own. The
threshold FORM lives in the frozen pre-registration; this file only applies it.
"""
p15_break_threshold_for(ece_rung0::Real) = p15_break_threshold(ece_rung0)

"""
    p15_break_rung(eces, thr) -> Int or nothing

The first index `k` in `eachindex(eces)` with `eces[k] > thr`, else `nothing`.

`nothing` MEANS "DID NOT BREAK WITHIN THE SWEPT LADDER" — D-13's UNBOUNDED reading — and it is
not a pass. The honest report in that case is the BOUND ("calibration holds at least out to the
last rung on this axis"), never a wider search.
"""
function p15_break_rung(eces::AbstractVector{<:Real}, thr::Real)
    return findfirst(>(float(thr)), eces)
end

"Resolve a per-column threshold from a scalar, a callable, or anything indexable by column."
function _p15_col_thr(thr, col::Integer)
    thr isa Real     && return float(thr)
    thr isa Function && return float(thr(col))
    return float(thr[col])
end

"""
    p15_axis_break(rows_for_axis, thr_by_column) -> (L_break = ..., per_column = ...)

The break rung for ONE axis. `rows_for_axis` is that axis's per-rung records — each a NamedTuple
carrying at least `rung` and `ece` (an 8-vector in `SBC_PARAM_LABELS` order) — and
`thr_by_column` is a scalar, a callable, or anything indexable by column index.

**ONLY COLUMNS 1 AND 8 ARE READ (D-05).** Those are the two SBC-calibrated targets (randomized-
rank atom handling) and the quantities a user acts on. The six nuisance marginals are measured and
reported at every rung and gate NOTHING: their ~0.08 SD drift is an already-disclosed named v2.0
limit (`docs/amortized.md`), so gating on it would re-fail the tool for a reason already on the
record and would collapse the envelope to near-zero on every axis, producing no map at all.

The returned `L_break` is the EARLIEST crossing across the two gating columns — the envelope ends
where EITHER target stops being calibrated — expressed as the record's own `rung`, not as a
positional index. `nothing` means neither gating column crossed within the ladder.

`P15_BREAK_COLUMNS == (1, 8)` is asserted at call time, so a later edit to the frozen tuple cannot
silently widen the gate through this reduction.
"""
function p15_axis_break(rows_for_axis, thr_by_column)
    @assert P15_BREAK_COLUMNS == (1, 8) "P15_BREAK_COLUMNS was widened away from the frozen " *
        "(1, 8) gating pair; the D-05 envelope definition is a pre-registration matter, not an " *
        "in-phase edit."
    rows = sort(collect(rows_for_axis); by = r -> r.rung)
    @assert all(r -> r.rung != P15_RUNG0, rows) "p15_axis_break: rung $(P15_RUNG0) is the " *
        "shared ANCHOR, measured once per net; it is not a ladder rung and must not be swept in."
    per_column = Dict{Int,Union{Nothing,Int}}()
    for col in P15_BREAK_COLUMNS
        eces = [float(r.ece[col]) for r in rows]
        k    = p15_break_rung(eces, _p15_col_thr(thr_by_column, col))
        per_column[col] = k === nothing ? nothing : rows[k].rung
    end
    crossed = [v for v in values(per_column) if v !== nothing]
    return (L_break = isempty(crossed) ? nothing : minimum(crossed),
            per_column = per_column)
end

# =============================================================================================
# The reporting-only record — everything that is measured and gates nothing
# =============================================================================================

"""
    p15_report_only(row) -> NamedTuple

Assemble the NON-GATING record for one rung: the 8-vector of KS p-values, the 8-vector of
chi-square p-values, the `sbc_coverage` nominal/empirical curves for the two gating columns, the
shrinkage and vacuous diagnostics, the per-column ECE for the six `P15_NUISANCE_COLUMNS`, and the
realised image-size counts. Every field here is REPORTED and gates NOTHING (D-04, D-05).

Fields absent from `row` come back as `nothing` rather than raising, because the sweep runner and
the fixtures populate different subsets and a reporting record must never be the thing that stops
a run.

The maximum-calibration-error field is NOT in the returned record — not even as a pass-through.
It is pinned at 0.99 on all eight columns by empty reliability bins and is an artifact, not a
measurement (D-04a).
"""
function p15_report_only(row)
    _f(name) = hasproperty(row, name) ? getproperty(row, name) : nothing

    ece = _f(:ece)
    nuisance_ece = ece === nothing ? nothing :
        Dict{String,Float64}(SBC_PARAM_LABELS[c] => float(ece[c]) for c in P15_NUISANCE_COLUMNS)

    ranks = _f(:ranks)
    coverage = nothing
    if ranks !== nothing
        L = _f(:L)
        L = L === nothing ? SBC_L : L
        coverage = Dict{String,NamedTuple}(
            SBC_PARAM_LABELS[c] => sbc_coverage(vec(view(ranks, :, c));
                                                levels = P15_COVERAGE_LEVELS, L = L)
            for c in P15_BREAK_COLUMNS)
    elseif _f(:coverage) !== nothing
        coverage = _f(:coverage)
    end

    return (axis                   = _f(:axis),
            rung                   = _f(:rung),
            mechanism              = _f(:mechanism),
            prior_edge             = _f(:prior_edge),
            ks_p                   = _f(:ks_p),
            chi2_p                 = _f(:chi2_p),
            coverage               = coverage,
            nuisance_ece           = nuisance_ece,
            shrinkage              = _f(:shrinkage),
            vacuous                = _f(:vacuous),
            realised_imsize_counts = _f(:realised_imsize_counts))
end

# =============================================================================================
# SC2 — the OOD crossing rule (D-07) and the THREE-VALUED per-axis verdict (D-06)
# =============================================================================================

"""
    p15_ood_rung(fire_rates, id_fire_rate; margin = P15_OOD_FIRE_MARGIN) -> Int or nothing

The first rung index `k` with `fire_rates[k] > id_fire_rate + margin`, else `nothing`.

**THE BASELINE IS THE MEASURED ONE, NEVER THE DESIGN VALUE (D-07).** `gate_ood_roc` returns
`id_fire_rate = mean(id_fused .> fused_thr)` — what this net, on this stream, actually does on
in-distribution images. Both alternatives were rejected before any rung existed:

  * **"ANY IMAGE FIRES"** is not a rule. The detector's in-distribution operating quantile is 0.95
    (`src/amortized/ood.jl:61`), so about 5 % of IN-DISTRIBUTION images fire BY CONSTRUCTION.
    An any-image rule therefore trips at rung 1 on every axis and proves nothing at all.
  * **A CHOSEN FIXED RATE** (50 %, say) would be an underived bar — exactly the Stage-1-ceiling
    pattern this project has already had to adjudicate, and the reason `p15_consts.jl` requires
    every bar to carry its derivation next to it.

The margin is binomial design arithmetic frozen in the pre-registration:
`P15_OOD_FIRE_MARGIN = q95(Binomial(P15_OOD_N_POS, 0.05))/P15_OOD_N_POS - 0.05 = 0.025`, derived
at load time by `_p15_binom_quantile` rather than typed.
"""
function p15_ood_rung(fire_rates::AbstractVector{<:Real}, id_fire_rate::Real;
                      margin::Real = P15_OOD_FIRE_MARGIN)
    thr = float(id_fire_rate) + float(margin)
    return findfirst(>(thr), fire_rates)
end

"""
    p15_axis_verdict(L_break, L_ood) -> Symbol

The SC2 verdict for ONE axis, drawn from `P15_SC2_STATES`. Three-valued, and pre-registered
before any rung existed (D-06):

  * `:protective`      — the flag fires AT OR BEFORE the rung where coverage breaks.
  * `:late`            — coverage breaks while the flag is still quiet. **THIS IS THE FAILURE SC2
                         EXISTS TO CATCH.** It includes the case where the flag never fires at all
                         while coverage did break.
  * `:silent_but_safe` — coverage never breaks within the swept ladder, so no warning was owed.

**WHY THE THIRD STATE EXISTS.** "The flag never fired" has two OPPOSITE meanings — working
correctly where nothing was wrong, versus failing to warn where something was — and a two-valued
pass/fail collapses them. Collapsing them is what would force a fifth after-the-fact amendment on
this milestone, so the third state is declared before results exist rather than discovered after.

**THE ASYMMETRY THAT MOTIVATES IT.** Spillover, registration and autofluorescence are IN-PRIOR
axes (`P15_AXIS_MECHANISM`): the detector's null was fitted on in-distribution data that already
CONTAINS them, so it may correctly never fire on them. Scoring that as a failure would penalise
the detector for being right.
"""
function p15_axis_verdict(L_break, L_ood)
    L_break === nothing && return :silent_but_safe
    (L_ood !== nothing && L_ood <= L_break) && return :protective
    return :late
end

"""
    p15_sc2_pass(verdicts) -> Bool

`true` iff no axis is `:late`. `verdicts` may be a Dict, a NamedTuple or any iterable of symbols.

**WHAT A LATE AXIS LICENSES, AND WHAT IT DOES NOT (D-14).** A LATE axis is reported as a NAMED
LIMIT — "the OOD flag does not protect against X" — and the phase CLOSES on it. That is a real,
publishable finding on the Phase-11/12 negative-but-useful precedent. Two responses are
EXPLICITLY FORBIDDEN: lowering the detector's in-distribution operating quantile until the flag
fires early enough (that is tuning a detector until it passes its own gate), and building a new
detector channel for the failing mechanism (that is new capability, with its own phase). No
retuning, no retraining.
"""
p15_sc2_pass(verdicts) = !any(==(:late), values(verdicts))

"""
    p15_axis_verdicts(break_by_axis, ood_by_axis) -> NamedTuple

Assemble the whole SC2 result over `P15_AXES`. `break_by_axis` and `ood_by_axis` map each axis to
its `L_break` / `L_ood` (`nothing` where there was no crossing); every axis in `P15_AXES` must be
present in both, because a missing key would silently read as "no break" and turn an INCOMPLETE
sweep into an `:unbounded_envelope` conclusion.

Each axis record carries its MECHANISM LABEL (`P15_AXIS_MECHANISM`) and its PRIOR BOUNDARY RUNG
(`P15_PRIOR_BOUNDARY_RUNG`), so the domain map reads as two overlaid stories rather than one
number: an IN-PRIOR axis that breaks is a TRAINING failure — the net saw those theta and still
miscalibrates — while an OUT-OF-MODEL axis that breaks is a SCOPE limit, expected and honest, and
precisely what the flag exists for. A boundary rung of `nothing` means "this axis has no prior
support to mark", never "unmeasured".

`degenerate` attaches D-13's two PRE-DECLARED readings mechanically instead of leaving them to be
recognised by eye:
  * `:empty_envelope`     — every axis breaks at rung 1. The tool has no usable operating range at
    the swept resolution. A reportable negative, not a phase failure, and it does NOT license
    widening the margin, lowering M or re-anchoring the threshold.
  * `:unbounded_envelope` — no axis breaks at any rung. The envelope exceeds the swept range; the
    honest report is the BOUND, not a wider search. This is the ONLY case in which
    `P15_ITERATION_ALLOWANCE` may be spent.
  * `nothing`             — neither.
"""
function p15_axis_verdicts(break_by_axis, ood_by_axis)
    missing_break = [a for a in P15_AXES if !haskey(break_by_axis, a)]
    missing_ood   = [a for a in P15_AXES if !haskey(ood_by_axis, a)]
    isempty(missing_break) || error("p15_axis_verdicts: no L_break for axes $(missing_break); " *
        "a missing axis would read as 'no break' and fake an unbounded envelope.")
    isempty(missing_ood) || error("p15_axis_verdicts: no L_ood for axes $(missing_ood); " *
        "a missing axis would read as 'flag never fired' and fake a LATE verdict.")

    per_axis = Dict{Symbol,NamedTuple}()
    verdicts = Dict{Symbol,Symbol}()
    for axis in P15_AXES
        lb = break_by_axis[axis]
        lo = ood_by_axis[axis]
        v  = p15_axis_verdict(lb, lo)
        v in P15_SC2_STATES || error("p15_axis_verdicts: verdict $v for $axis is outside the " *
            "pre-registered P15_SC2_STATES $(P15_SC2_STATES).")
        verdicts[axis] = v
        per_axis[axis] = (axis                = axis,
                          verdict             = v,
                          L_break             = lb,
                          L_ood               = lo,
                          mechanism           = P15_AXIS_MECHANISM[axis],
                          prior_boundary_rung = P15_PRIOR_BOUNDARY_RUNG[axis])
    end

    breaks     = [break_by_axis[a] for a in P15_AXES]
    degenerate = all(b -> b == 1, breaks)         ? :empty_envelope :
                 all(b -> b === nothing, breaks)  ? :unbounded_envelope : nothing

    return (axes = per_axis, verdicts = verdicts,
            sc2_pass = p15_sc2_pass(verdicts), degenerate = degenerate)
end

