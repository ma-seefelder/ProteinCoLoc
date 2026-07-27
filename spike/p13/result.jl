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

# spike/p13/result.jl --- ThreeHypothesisColocResult, the spike-lane three-way result type (D-08).
#
# (a) THE NAME THIS FILE EXISTS TO REPLACE. The sketch at src/results.jl:161-192 gives the
#     Phase-13 result a field named `log_bf_simplex :: NTuple{3,Float64}`. That name is a
#     MISNOMER and D-08 replaces it. The three values are
#         log BF(coloc : random),   0.0,   log BF(exclusion : random)
#     i.e. TWO FREE NUMBERS PLUS A STRUCTURAL ZERO -- the zero is the reference class scored
#     against itself, not a third measurement. A point on a 2-simplex has non-negative entries
#     that SUM TO ONE; these are unbounded real log Bayes factors on the whole line, one of which
#     is pinned at zero by construction. A reader who sees "simplex" and never reads the comment
#     will quote them as model probabilities, which is exactly the reading D-08 rejected when it
#     declined to invent a prior over the three hypotheses (turning evidence into a chosen
#     hypothesis is Phase 14's job, not this type's). The replacement name states its own
#     reference class: `log_bf_vs_random`.
#
# (b) WHERE THE EXECUTABLE TYPE LIVES, AND WHY THE SKETCH IS NOT EDITED. D-01 scopes all
#     Phase-13 work to spike/, and CLAUDE.md requires src/ to stay provably untouched during the
#     spike. So the EXECUTABLE subtype is defined HERE, in spike/p13/, reached to its supertype
#     by a READ-ONLY include of src/results.jl -- the same read-only reach spike/contract.jl:46-47
#     already uses for src/LoadImages.jl and src/colocalization.jl. The comment-only sketch block
#     at src/results.jl:161-192 stays BYTE-UNCHANGED. The rename there is RECORDED as a one-line
#     productionization item (P13_RESULTS_RENAME_DEFERRED, spike/p13/consts.jl section J) and
#     carried in the phase report, NOT performed. It is a src/ edit with no scientific benefit
#     inside this phase; if the user wants it, it is a separate docs-only commit and a user
#     decision, not a planner inference (13-RESEARCH G3).
#
# (c) PHASE 7 D-02 IS STILL SATISFIED. D-02 requires every new result variant to slot in as a NEW
#     subtype of the real `AbstractColocResult`, never as bolted-on fields on an existing struct.
#     That is structurally what happens here: `ThreeHypothesisColocResult <: AbstractColocResult`
#     against the ACTUAL supertype from src/results.jl, and the four shared accessors are
#     implemented for it. Defining it in spike/ changes nothing about that relationship.
#
# (d) DECOUPLING (hard constraint, CLAUDE.md): spike-local code. The two reaches into src/ are
#     read-only `include`s; no src/ byte changes, no package is added, and nothing here is
#     imported by the shipped package.
#
# WHAT THIS FILE DELIBERATELY DOES NOT CONTAIN. No TOST, no equivalence test, no p-value
# machinery (D-14): the D-13 calibration verdict is a traffic-light BAND over ECE, not a
# hypothesis test, so there is no over-power failure mode of the kind the M = 2000 SBC
# point-null had and no equivalence machinery is needed. And no accessor carrying the simplex
# misnomer forward -- the sketch's own accessor name does not survive into the spike lane either.

# ORDER MATTERS (spike/contract.jl:34-39, the documented load order): StatsBase + Statistics must
# be in scope before src/colocalization.jl is reached transitively, and Images before
# src/LoadImages.jl (its convenience constructor calls Images.otsu_threshold). Both are reached
# through spike/validation/sbc.jl -> harness.jl -> contract.jl below, so the `using` lines come
# FIRST and the includes after.
using StatsBase      # corspearman/corkendall for the transitively-reached correlation() Dict
using Statistics     # mean, quantile
using Images         # otsu_threshold (transitive src/LoadImages.jl requirement)

# --- Guarded includes, in dependency order --------------------------------------------------
# READ-ONLY reach into src/ for the supertype, `_iface_error`, `OODVerdict` and `CalibrationMeta`.
# NOTE THE PATH DEPTH: from spike/p13/ the repo root is TWO levels up, not three.
isdefined(@__MODULE__, :AbstractColocResult) ||
    include(joinpath(@__DIR__, "..", "..", "src", "results.jl"))
# The Tier-1 pre-registration: P13_ECE_GREEN / P13_ECE_YELLOW / P13_ECE_NBINS /
# P13_GATE_STATISTIC / P13_VACUOUS_AUC_FLOOR / P13_RESULTS_RENAME_DEFERRED.
isdefined(@__MODULE__, :P13_DEV_SEED) || include(joinpath(@__DIR__, "consts.jl"))
# The shared, gate-lineage calibration surface: CalibrationResult + _bin_calibration. Guarded on
# the FUNCTION rather than the struct, because the function is what this file's meta-builder
# consumes. NEVER hand-modified here -- see the empty-bin note on p13_calibration_meta.
isdefined(@__MODULE__, :_bin_calibration) ||
    include(joinpath(@__DIR__, "..", "validation", "sbc.jl"))

# THE GUARD BLOCK COVERS ONLY THE `const`s AND THE `struct`, NOT THE FUNCTIONS, AND THAT SPLIT IS
# DELIBERATE. Julia 1.12 DROPS every docstring written inside an `if ... end` block -- the parser
# emits the `Core.@doc` call but the docsystem never registers it, so `@doc f` returns `nothing`
# for a function documented inside a guard (verified on 1.12.6). D-08 requires the `bayes_factor`
# choice and the do-not-fabricate-Delta-rho rule to be DOCUMENTED, and the accompanying test
# asserts the docstring is retrievable, so the functions live at TOP LEVEL where the docsystem
# can see them. Method redefinition under a re-include is silent and harmless in Julia; only
# `const` and `struct` redefinition would warn or throw, and those are what the guard protects.
if !isdefined(@__MODULE__, :ThreeHypothesisColocResult)

# --- The reference-class-carrying evidence triple -------------------------------------------
#
# The KEY ORDER of the stored NamedTuple is (:coloc, :random, :exclusion). The sketch's positional
# `NTuple{3,Float64}` was documented as `{null, coloc, anti-coloc}` -- REFERENCE CLASS FIRST --
# which no reader would guess from a name containing "simplex", and which a caller indexing `[1]`
# for "the colocalization evidence" would get silently wrong. A NamedTuple removes the hazard
# outright: every read site is by NAME, so a reordering cannot change a number.
const P13_LOG_BF_KEYS = (:coloc, :random, :exclusion)

"""
    ThreeWayLogBF

The stored type of [`ThreeHypothesisColocResult`](@ref)'s evidence field:
`NamedTuple{(:coloc, :random, :exclusion), Tuple{Float64,Float64,Float64}}`.

`coloc` is `log BF(coloc : random)`, `exclusion` is `log BF(exclusion : random)`, and `random`
is the reference class scored against itself, hence identically `0.0` (D-08). These are
unbounded log Bayes factors, not probabilities and not a point on a simplex.
"""
const ThreeWayLogBF = NamedTuple{P13_LOG_BF_KEYS, Tuple{Float64,Float64,Float64}}

"""
    ThreeHypothesisColocResult(grid, posterior, log_bf_vs_random, ood, calibration, meta)

Result of the Phase-13 three-hypothesis amortized evidence path: one forward pass through the
two-head evidence net yields log Bayes factors for `coloc` and `exclusion`, each against the
`random` reference, with an OOD verdict and the D-13 calibration report.

A new subtype of the real `AbstractColocResult` (Phase 7 D-02), defined in the spike lane per
D-01; `src/results.jl` is byte-unchanged and its comment-only sketch is superseded by this file
rather than edited.

# Fields
- `grid::Int`: patch grid `G` the evidence net was trained for.
- `posterior::Matrix{Float64}`: posterior draws backing the result (parameter x draw layout).
- `log_bf_vs_random::ThreeWayLogBF`: `(coloc, random, exclusion)` log Bayes factors against the
  random reference; `random` is a STRUCTURAL zero and the constructor enforces it.
- `ood::OODVerdict`: out-of-distribution / misspecification verdict.
- `calibration::CalibrationMeta`: the D-13 calibration report (build it with
  [`p13_calibration_meta`](@ref), which carries the head AUC and the empty-bin count).
- `meta::NamedTuple`: run metadata. If -- and ONLY if -- the run actually carried control
  posterior draws, put them in `meta.delta_rho_draws`; see [`delta_rho`](@ref).

The constructor accepts the evidence triple with its three keys in any order and throws an
`ArgumentError` naming D-08 if `random` is not exactly `0.0`.
"""
struct ThreeHypothesisColocResult <: AbstractColocResult
    grid             :: Int
    posterior        :: Matrix{Float64}
    log_bf_vs_random :: ThreeWayLogBF
    ood              :: OODVerdict
    calibration      :: CalibrationMeta
    meta             :: NamedTuple

    function ThreeHypothesisColocResult(grid::Integer, posterior::AbstractMatrix{<:Real},
                                        log_bf_vs_random::NamedTuple, ood::OODVerdict,
                                        calibration::CalibrationMeta, meta::NamedTuple)
        lbf = _p13_as_three_way_logbf(log_bf_vs_random)
        # STRUCTURAL, NOT MERELY CONVENTIONAL. D-08 emits evidence AGAINST the random reference,
        # so the reference entry is the reference scored against itself: identically zero, for
        # every input, by construction. A non-zero value here means the caller built the triple
        # against some other reference (or fabricated a third measurement), and every downstream
        # comparison would silently be on a different scale.
        lbf.random === 0.0 || throw(ArgumentError(
            "ThreeHypothesisColocResult (D-08): log_bf_vs_random.random must be exactly 0.0 " *
            "-- the random class is the REFERENCE, scored against itself. Got $(lbf.random)."))
        return new(Int(grid), Matrix{Float64}(posterior), lbf, ood, calibration, meta)
    end
end

end # if !isdefined(@__MODULE__, :ThreeHypothesisColocResult) -- consts + struct only

"""
    _p13_as_three_way_logbf(nt::NamedTuple) -> ThreeWayLogBF

Normalise any NamedTuple carrying exactly the keys `:coloc`, `:random` and `:exclusion` (in ANY
order) into the canonical `(:coloc, :random, :exclusion)` layout.

This exists so the read surface's own output ordering cannot become a silent coupling:
`three_way_log_bf` (`spike/p13/net.jl`) returns `(coloc, exclusion, random)`, and a positional
struct would have mis-stored it. Missing or extra keys throw an `ArgumentError` naming D-08
rather than being silently dropped.
"""
function _p13_as_three_way_logbf(nt::NamedTuple)
    Set(keys(nt)) == Set(P13_LOG_BF_KEYS) || throw(ArgumentError(
        "ThreeHypothesisColocResult (D-08): the evidence triple must carry exactly the keys " *
        "$(P13_LOG_BF_KEYS); got $(keys(nt))."))
    return ThreeWayLogBF((Float64(nt.coloc), Float64(nt.random), Float64(nt.exclusion)))
end

# --- The four shared accessors (the AbstractColocResult interface) ---------------------------

"""
    posterior_draws(r::ThreeHypothesisColocResult)

The posterior draws backing the result (parameter x draw layout), unchanged from the shipped
`AmortizedColocResult` semantics.
"""
posterior_draws(r::ThreeHypothesisColocResult) = r.posterior

"""
    is_ood(r::ThreeHypothesisColocResult)

Whether the input was flagged out-of-distribution / misspecified, unchanged from the shipped
`AmortizedColocResult` semantics.
"""
is_ood(r::ThreeHypothesisColocResult) = r.ood.flag

"""
    bayes_factor(r::ThreeHypothesisColocResult) -> Float64

Return `r.log_bf_vs_random.coloc`, i.e. `log BF(coloc : random)`.

THIS CHOICE IS DOCUMENTED RATHER THAN SILENT, because a single-number accessor on a three-way
result is inherently lossy and the reader deserves to know which number they got. The shared
interface's `bayes_factor` has exactly one shipped meaning -- "the colocalization Bayes factor
evidence" (`src/results.jl:59-65`) -- and the coloc-vs-random entry is the semantically closest
quantity this result carries. Returning it means a consumer written against the binary
`AmortizedColocResult` keeps working and gets the number it meant, rather than an error.

A caller who wants THE THREE-WAY EVIDENCE must use [`log_bf_vs_random`](@ref), which returns all
three entries by name. `bayes_factor` deliberately says nothing about the exclusion hypothesis;
reading a three-way verdict off this one number would be wrong.
"""
bayes_factor(r::ThreeHypothesisColocResult) = r.log_bf_vs_random.coloc

"""
    delta_rho(r::ThreeHypothesisColocResult)

Return the carried Delta rho draws if -- and only if -- the run actually produced them, i.e. if
`meta.delta_rho_draws` is present. Otherwise fall through to `_iface_error`.

DO NOT FABRICATE A DELTA RHO. Delta rho is a sample-minus-control contrast that requires a
CONTROL posterior; the three-way evidence path does not need one and in general does not carry
one. Synthesising a plausible-looking value (zeros, the coloc draws, a difference against a
notional zero control) would put a number in front of a reader that no measurement backs. Letting
the documented interface error fire is the honest answer, and it is exactly what
`_iface_error` (`src/results.jl:47-49`) exists for: it names the accessor, the type, and the four
accessors the interface requires.

D-05's own counter-example is why this matters: a sample at rho = 0.8 under a control at rho = 0.9
has Delta rho < 0, which a naive three-way reading would publish as "mutually exclusive" while
the sample is strongly colocalized.
"""
function delta_rho(r::ThreeHypothesisColocResult)
    haskey(r.meta, :delta_rho_draws) || _iface_error(r, :delta_rho)
    return r.meta.delta_rho_draws
end

# --- The named three-way accessor ------------------------------------------------------------

"""
    log_bf_vs_random(r::ThreeHypothesisColocResult) -> ThreeWayLogBF

The full evidence triple by name: `(coloc = log BF(C:R), random = 0.0, exclusion = log BF(E:R))`.

This is the accessor a three-way consumer should use. There is deliberately NO alias carrying the
sketch's simplex naming: the misnomer does not survive into the spike lane, and nothing consumes
the sketch today (it is a comment), so a clean name costs nothing now and everything later
(13-RESEARCH G2).
"""
log_bf_vs_random(r::ThreeHypothesisColocResult) = r.log_bf_vs_random

# --- The D-13 calibration surface ------------------------------------------------------------

"""
    p13_traffic_light(ece; green = P13_ECE_GREEN, yellow = P13_ECE_YELLOW) -> Symbol

The Phase-13 calibration verdict: `:green` if `ece <= green`, `:yellow` if `ece <= yellow`,
else `:red`.

RE-DECLARED LOCALLY rather than imported from `spike/validation/sbc.jl`'s `sbc_traffic_light`,
matching the established pattern that every pre-registration in this project is SELF-CONTAINED
(each `test/gate/gate_consts_*.jl` re-declares its own band). The defaults come from the frozen
`spike/p13/consts.jl`, so this function can be read and audited without chasing another phase's
current value -- and if Phase 5's band ever moved, Phase 13's reported verdicts would not move
with it.

It is a BAND, not a hypothesis test (D-14).
"""
function p13_traffic_light(ece::Real; green::Real = P13_ECE_GREEN,
                           yellow::Real = P13_ECE_YELLOW)
    ece <= green  && return :green
    ece <= yellow && return :yellow
    return :red
end

"""
    vacuous_pass(ece, auc; ece_green = P13_ECE_GREEN, auc_floor = P13_VACUOUS_AUC_FLOOR) -> Bool

`true` when the ECE is green but the AUC is at or below `auc_floor` -- a VACUOUS PASS.

A head that learns nothing produces well-calibrated-looking probabilities by simply reproducing
the base rate: "0.5, always" is perfectly calibrated and perfectly useless. That is the F3
vacuous-SBC-pass failure of `07-CALIBRATION-FINDINGS.md` one level up (13-RESEARCH Pitfall 9), and
it is why every calibration claim in this phase ships with its discrimination number beside it.
`p13_calibration_meta` makes that structural rather than a reporting habit.

A `true` return is NOT a failure verdict on its own -- it is a LABEL. The result must be reported
as a vacuous pass rather than quoted as calibration evidence.
"""
vacuous_pass(ece::Real, auc::Real; ece_green::Real = P13_ECE_GREEN,
             auc_floor::Real = P13_VACUOUS_AUC_FLOOR) = (ece <= ece_green) && (auc <= auc_floor)

"""
    p13_calibration_meta(cal::CalibrationResult; grid, auc, gate = NamedTuple()) -> CalibrationMeta

Build the `CalibrationMeta` a [`ThreeHypothesisColocResult`](@ref) carries from a spike
`CalibrationResult`, the grid, and the head's AUC.

FIELD-ORDER CONTRACT, STATED BECAUSE THE TWO TYPES DIFFER. The spike's `CalibrationResult`
(`spike/validation/sbc.jl:65-72`) has SIX fields: `bin_midpoints`, `predicted_rate`,
`observed_rate`, `bin_counts`, `ece`, `mce`. `CalibrationMeta` (`src/results.jl:117-126`) has
those same six IN THE SAME ORDER plus `grid::Int` and `gate::NamedTuple`. THE CARRIED VALUE IS
THE EIGHT-FIELD `CalibrationMeta`: the two extra fields are the provenance that lets a result
prove which pre-registered gate it was calibrated under, which is the whole point of attaching it.

`auc` IS A REQUIRED KEYWORD, NOT AN OPTION. The vacuous-pass guard is only real if it cannot be
omitted, so the discrimination number is part of the constructor's contract: there is no way to
build a Phase-13 calibration verdict that does not carry the AUC it must be read against.

The returned `gate` NamedTuple always carries, at minimum:

  - `statistic = P13_GATE_STATISTIC` (`:ece`) -- ECE IS THE GATE STATISTIC AND MCE IS NOT.
  - `ece_green`, `ece_yellow`, `n_bins` -- the frozen band and binning.
  - `mce_reported_not_gated = true` and `empty_bins` -- THE EMPTY-BIN TRAP, recorded beside the
    MCE it invalidates. `_bin_calibration` scores an EMPTY bin as `predicted_rate = midpoint`,
    `observed_rate = 0.0`, which contributes ZERO ECE weight (the ECE term is weighted by
    `count_i/total`) but FULL MCE weight (MCE is an unweighted max over bins). On a
    well-separated three-way problem many of the ten bins are empty, so MCE would be DOMINATED BY
    EMPTY BINS and would fail the gate for a reason that has nothing to do with calibration.
    `_bin_calibration` is shared, gate-lineage code and is NOT hand-modified to work around this;
    the workaround is the CHOICE OF GATE STATISTIC, and the empty-bin count travels with the
    number so a reader can see how much of the MCE is an artifact.
  - `auc` and `vacuous` -- the discrimination number and the Pitfall-9 label.
  - `verdict` -- the `p13_traffic_light` band verdict on the ECE.

Any caller-supplied `gate` (run provenance such as `(; seed, consts_hash, passed)`) is MERGED
UNDER these, so caller provenance is preserved but the required keys can never be overwritten or
omitted.
"""
function p13_calibration_meta(cal::CalibrationResult; grid::Integer, auc::Real,
                              gate::NamedTuple = NamedTuple())
    empty_bins = count(iszero, cal.bin_counts)
    required = (statistic              = P13_GATE_STATISTIC,
                ece_green              = Float64(P13_ECE_GREEN),
                ece_yellow             = Float64(P13_ECE_YELLOW),
                n_bins                 = Int(P13_ECE_NBINS),
                mce_reported_not_gated = true,
                empty_bins             = empty_bins,
                auc                    = Float64(auc),
                vacuous                = vacuous_pass(cal.ece, auc),
                verdict                = p13_traffic_light(cal.ece))
    return CalibrationMeta(cal.bin_midpoints, cal.predicted_rate, cal.observed_rate,
                           cal.bin_counts, cal.ece, cal.mce, Int(grid),
                           merge(gate, required))
end
