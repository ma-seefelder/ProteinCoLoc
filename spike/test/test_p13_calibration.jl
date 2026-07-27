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

# spike/test/test_p13_calibration.jl --- the D-13 / D-14 evidence-calibration gate.
#
# WHAT THIS FILE IS FOR. D-13 gates on the reliability of the heads' implied class probabilities,
# using the SHARED `_bin_calibration` / `CalibrationResult` surface. Three properties of that
# arrangement are load-bearing and none of them is obvious from reading the constant that records
# them, so each is DEMONSTRATED here rather than asserted in prose:
#
#   1. THE STATISTIC WORKS IN BOTH DIRECTIONS. A calibration check that can only return "green"
#      proves nothing. The Phase-5 both-directions principle (test_sbc.jl:77-88) is copied against
#      the Phase-13 band: perfectly calibrated input must be green, maximally miscalibrated input
#      must be red.
#   2. THE EMPTY-BIN TRAP IS REAL (13-RESEARCH Pitfall 7). `_bin_calibration` scores an EMPTY bin
#      as `predicted_rate = midpoint`, `observed_rate = 0.0`. That carries ZERO ECE weight but FULL
#      MCE weight, so on a well-separated problem MCE is dominated by bins that contain no data.
#      This testset IS the documentation of why `P13_GATE_STATISTIC === :ece`, and it exists so
#      nobody "fixes" the shared, gate-lineage `_bin_calibration` to make MCE usable.
#   3. A GREEN ECE IS NOT EVIDENCE ON ITS OWN (Pitfall 9). A head that learns nothing reproduces
#      the base rate and looks perfectly calibrated. Every Phase-13 calibration verdict therefore
#      ships with its discrimination number, and `p13_calibration_meta` makes that structural.
#
# D-14: the traffic light is a BAND over ECE, not a hypothesis test, so there is no over-power
# failure mode of the kind the M = 2000 SBC point-null had and NO equivalence machinery is
# introduced anywhere in the phase. That absence is asserted, not just intended.
#
# `spike/validation/sbc.jl` is the unit under test and is NEVER modified by this phase.

using Test
using Statistics

# Units under test: the Phase-13 calibration surface (which pulls consts.jl, src/results.jl
# read-only and validation/sbc.jl), the class enum, and the hand-rolled ROC/AUC.
isdefined(@__MODULE__, :p13_calibration_meta) ||
    include(joinpath(@__DIR__, "..", "p13", "result.jl"))
isdefined(@__MODULE__, :_bin_calibration) ||
    include(joinpath(@__DIR__, "..", "validation", "sbc.jl"))
isdefined(@__MODULE__, :ThreeWayClass) ||
    include(joinpath(@__DIR__, "..", "p13", "labels.jl"))
# roc_auc: the in-repo, hand-rolled, tie-aware Mann-Whitney AUC. REUSED, never re-implemented --
# the same rule tau_probe.jl follows, and the reason no ROC package may enter the environment.
isdefined(@__MODULE__, :roc_auc) ||
    include(joinpath(@__DIR__, "..", "validation", "ood.jl"))

# --- The three-class synthetic reliability fixture, built ONCE -------------------------------
#
# 300 labelled samples, 100 per class, with each head's probability set to a value that is
# well-calibrated ON ITS OWN CLASS PAIR and deliberately WRONG on the third class -- which is
# exactly what a head trained under D-11's per-head restriction does, because it never saw the
# third class. The numbers are chosen so the restricted curve is green and the pooled curve is
# red, making the restriction's effect impossible to miss. Nothing here is random: the fixture is
# deterministic, so no Phase-13 stream (reported or fixture) is consumed at all.
const CAL_FIX_LABELS = vcat(fill(EXCLUSION, 100), fill(RANDOM, 100), fill(COLOC, 100))
const CAL_FIX_P_COLOC     = [l === COLOC ? 0.96 : (l === RANDOM ? 0.04 : 0.96)
                             for l in CAL_FIX_LABELS]
const CAL_FIX_P_EXCLUSION = [l === EXCLUSION ? 0.96 : (l === RANDOM ? 0.04 : 0.96)
                             for l in CAL_FIX_LABELS]

"The D-11 per-head evaluation mask: a head is scored ONLY on its own class pair."
_head_mask(labels, positive) = [l === positive || l === RANDOM for l in labels]

# The executable source of every Phase-13 file, comments stripped, for the D-14 absence grep.
const P13_SRC_FILES = sort(filter(f -> endswith(f, ".jl"),
                                  readdir(joinpath(@__DIR__, "..", "p13"); join = true)))
const P13_ALL_CODE = Dict(
    basename(f) => join(filter(l -> !startswith(strip(l), "#"),
                               split(read(f, String), '\n')), '\n')
    for f in P13_SRC_FILES)

@testset "P13 evidence calibration (D-13, D-14)" verbose = true begin

    @testset "both directions" begin
        # COPIED PRINCIPLE (test_sbc.jl:77-88), re-run against the PHASE-13 band and the
        # PHASE-13 bin count, because Phase 13 re-declares its own pre-registration rather than
        # importing Phase 5's -- if Phase 5's band ever moved, Phase 13's verdicts must not move
        # with it.
        cal_good = _bin_calibration(fill(0.5, 100), vcat(trues(50), falses(50));
                                    n_bins = P13_ECE_NBINS)
        @test cal_good.ece < P13_ECE_GREEN
        @test p13_traffic_light(cal_good.ece) === :green
        @test p13_traffic_light(cal_good.ece) isa Symbol

        cal_bad = _bin_calibration(fill(1.0, 100), falses(100); n_bins = P13_ECE_NBINS)
        @test cal_bad.ece > P13_ECE_YELLOW
        @test p13_traffic_light(cal_bad.ece) === :red

        # The band is ordered and the middle rung is reachable, so ":green" is not the only
        # verdict the function can return.
        @test p13_traffic_light((P13_ECE_GREEN + P13_ECE_YELLOW) / 2) === :yellow
        @test P13_ECE_GREEN < P13_ECE_YELLOW
    end

    @testset "THE EMPTY-BIN TRAP: MCE is not a usable gate statistic" begin
        # THE CONSTRUCTION. All 100 predicted probabilities sit at 0.05 and exactly 5 of the 100
        # samples are positive, so the single occupied bin is PERFECTLY calibrated: predicted 0.05,
        # observed 0.05, gap 0. The other nine bins are EMPTY.
        #
        # WHAT `_bin_calibration` DOES WITH AN EMPTY BIN (spike/validation/sbc.jl:109-115): it
        # records `predicted_rate = midpoint` and `observed_rate = 0.0`. The per-bin gap for the
        # top bin is therefore |0.95 - 0.0| = 0.95 -- computed entirely from a bin that contains
        # NO DATA.
        #
        # WHY THAT MATTERS. ECE weights each gap by `count_i / total`, so an empty bin contributes
        # EXACTLY ZERO to the ECE (sbc.jl:118-126). MCE is an UNWEIGHTED max over bins, so the
        # empty bin contributes its full 0.95. The result below is a reliability input that is
        # flawless by construction and yet scores MCE ~ 0.95: an MCE gate would have failed it for
        # a reason that has nothing whatsoever to do with calibration.
        #
        # On a WELL-SEPARATED three-way problem -- which is exactly what a working evidence net
        # produces, probabilities piled up near 0 and near 1 -- most of the ten bins are empty, so
        # this is the normal case for Phase 13, not a contrived corner.
        #
        # THE FIX IS THE CHOICE OF GATE STATISTIC, NOT A CODE CHANGE. `_bin_calibration` is shared,
        # gate-lineage code used by the Phase-5 SBC gate; hand-modifying it to special-case empty
        # bins would silently change a number that has already been reported. Phase 13 instead
        # pre-registers `P13_GATE_STATISTIC = :ece` and REPORTS the MCE beside its empty-bin count.
        # THIS TESTSET IS THE DOCUMENTATION OF THAT DECISION.
        cal = _bin_calibration(fill(0.05, 100), vcat(trues(5), falses(95));
                               n_bins = P13_ECE_NBINS)

        n_empty = count(iszero, cal.bin_counts)
        @test n_empty > 0
        @test n_empty == P13_ECE_NBINS - 1        # exactly one occupied bin, nine empty

        # The occupied bin really is perfectly calibrated, so the MCE below cannot be blamed on it.
        occupied = findall(!iszero, cal.bin_counts)
        @test length(occupied) == 1
        @test cal.predicted_rate[occupied[1]] ≈ cal.observed_rate[occupied[1]] atol = 1e-12

        # ECE is green; MCE is enormous. Same input, same function, same call.
        @test cal.ece < P13_ECE_GREEN
        @test p13_traffic_light(cal.ece) === :green
        @test cal.mce > P13_ECE_YELLOW
        @test cal.mce > 0.9

        # THE ATTRIBUTION, ASSERTED AND NOT ASSUMED: the bin that MAXIMISES the gap -- i.e. the bin
        # that IS the MCE -- is one of the empty ones. Without this the testset would only show
        # "MCE is big", not "MCE is big BECAUSE of empty bins".
        worst = argmax(abs.(cal.predicted_rate .- cal.observed_rate))
        @test iszero(cal.bin_counts[worst])
        @test cal.predicted_rate[worst] == cal.bin_midpoints[worst]
        @test cal.observed_rate[worst] == 0.0

        # And the pre-registration says so, in the same words.
        @test P13_GATE_STATISTIC === :ece
        @test P13_GATE_STATISTIC !== :mce
    end

    @testset "the gate statistic is ECE by pre-registration" begin
        # Asserted as LITERALS, deliberately. A test that reads the constant it is checking cannot
        # notice the constant changing, and spike/p13/consts.jl is a frozen pre-registration whose
        # values are the whole basis for the phase's numbers being citable.
        @test P13_GATE_STATISTIC === :ece
        @test P13_ECE_NBINS  == 10
        @test P13_ECE_GREEN  == 0.05
        @test P13_ECE_YELLOW == 0.10
        @test P13_VACUOUS_AUC_FLOOR == 0.60
        # The meta builder records that MCE was reported rather than gated, so a downstream reader
        # of a result never has to reconstruct this decision from a document.
        meta = p13_calibration_meta(
            _bin_calibration(fill(0.5, 100), vcat(trues(50), falses(50));
                             n_bins = P13_ECE_NBINS);
            grid = 8, auc = 0.97)
        @test meta.gate.statistic === :ece
        @test meta.gate.mce_reported_not_gated === true
        @test meta.gate.n_bins == P13_ECE_NBINS
        @test isfinite(meta.mce)                 # reported, and carried
        @test meta isa CalibrationMeta
    end

    @testset "per-head restriction" begin
        # D-11: each head trains only on its own class pair, so each reliability curve must be
        # evaluated only on that pair. Feeding the coloc head an exclusion-class sample asks it
        # about a class it never saw -- the answer is extrapolation, and scoring it as calibration
        # error attributes a training-distribution question to the calibration statistic.
        coloc_mask = _head_mask(CAL_FIX_LABELS, COLOC)
        excl_mask  = _head_mask(CAL_FIX_LABELS, EXCLUSION)
        @test count(coloc_mask) == 200 && count(excl_mask) == 200
        @test !any(coloc_mask[CAL_FIX_LABELS .=== EXCLUSION])
        @test !any(excl_mask[CAL_FIX_LABELS .=== COLOC])

        cal_coloc_restricted = _bin_calibration(CAL_FIX_P_COLOC[coloc_mask],
                                                CAL_FIX_LABELS[coloc_mask] .=== COLOC;
                                                n_bins = P13_ECE_NBINS)
        cal_excl_restricted  = _bin_calibration(CAL_FIX_P_EXCLUSION[excl_mask],
                                                CAL_FIX_LABELS[excl_mask] .=== EXCLUSION;
                                                n_bins = P13_ECE_NBINS)
        @test p13_traffic_light(cal_coloc_restricted.ece) === :green
        @test p13_traffic_light(cal_excl_restricted.ece)  === :green

        # The SAME head scored on ALL THREE classes -- the mistake this restriction prevents.
        cal_coloc_pooled = _bin_calibration(CAL_FIX_P_COLOC,
                                            CAL_FIX_LABELS .=== COLOC;
                                            n_bins = P13_ECE_NBINS)
        cal_excl_pooled  = _bin_calibration(CAL_FIX_P_EXCLUSION,
                                            CAL_FIX_LABELS .=== EXCLUSION;
                                            n_bins = P13_ECE_NBINS)

        # LOAD-BEARING, NOT COSMETIC: including the third class changes the ECE, and changes the
        # VERDICT. If these were equal the restriction would be a convention with no consequence
        # and this testset would be decorative.
        @test cal_coloc_pooled.ece != cal_coloc_restricted.ece
        @test cal_excl_pooled.ece  != cal_excl_restricted.ece
        @test cal_coloc_pooled.ece > cal_coloc_restricted.ece
        @test cal_excl_pooled.ece  > cal_excl_restricted.ece
        @test p13_traffic_light(cal_coloc_pooled.ece) === :red
        @test p13_traffic_light(cal_excl_pooled.ece)  === :red
    end

    @testset "vacuous-pass guard (Pitfall 9)" begin
        # A HEAD THAT LEARNS NOTHING. It emits the base rate for every input, so its reliability
        # diagram is flawless -- ECE 0, green -- and it carries no information at all: its AUC is
        # exactly 0.5. This is the F3 vacuous-SBC-pass failure one level up, and a report that
        # quoted the green ECE without the AUC beside it would be quoting nothing.
        n = P13_MIN_EVAL_PER_HEAD
        flat_p = fill(0.5, 2n)
        flat_y = vcat(trues(n), falses(n))
        cal_flat = _bin_calibration(flat_p, flat_y; n_bins = P13_ECE_NBINS)
        auc_flat = roc_auc(flat_p[.!flat_y], flat_p[flat_y])[3]
        @test p13_traffic_light(cal_flat.ece) === :green
        @test auc_flat == 0.5
        @test vacuous_pass(cal_flat.ece, auc_flat)

        # A DISCRIMINATING HEAD, calibrated to the SAME green band, so the only thing separating
        # the two cases is the discrimination number.
        sharp_p = [y ? 0.96 : 0.04 for y in flat_y]
        cal_sharp = _bin_calibration(sharp_p, flat_y; n_bins = P13_ECE_NBINS)
        auc_sharp = roc_auc(sharp_p[.!flat_y], sharp_p[flat_y])[3]
        @test p13_traffic_light(cal_sharp.ece) === :green
        @test auc_sharp == 1.0
        @test !vacuous_pass(cal_sharp.ece, auc_sharp)

        # The floor is what separates them, and it is the pre-registered one.
        @test vacuous_pass(cal_sharp.ece, P13_VACUOUS_AUC_FLOOR)
        @test !vacuous_pass(cal_sharp.ece, nextfloat(P13_VACUOUS_AUC_FLOOR))
        # A RED ece is never a "vacuous pass" -- it is not a pass at all.
        @test !vacuous_pass(0.5, 0.5)

        # STRUCTURAL, NOT A REPORTING HABIT: the gate NamedTuple a result carries always holds the
        # head AUC and the empty-bin count, so no Phase-13 calibration verdict can be quoted
        # without the two numbers it must be read against.
        meta_flat  = p13_calibration_meta(cal_flat;  grid = 8, auc = auc_flat)
        meta_sharp = p13_calibration_meta(cal_sharp; grid = 8, auc = auc_sharp)
        @test haskey(meta_flat.gate, :auc)
        @test haskey(meta_flat.gate, :empty_bins)
        @test meta_flat.gate.auc == 0.5
        @test meta_flat.gate.vacuous === true
        @test meta_sharp.gate.vacuous === false
        @test meta_sharp.gate.empty_bins == count(iszero, cal_sharp.bin_counts)
        @test meta_sharp.gate.empty_bins > 0     # a sharp head DOES empty most of its bins
        # Caller provenance is preserved, but the required keys cannot be dropped or overwritten.
        meta_prov = p13_calibration_meta(cal_sharp; grid = 8, auc = auc_sharp,
                                         gate = (; seed = 0x0B13DE71, auc = -1.0))
        @test meta_prov.gate.seed == 0x0B13DE71
        @test meta_prov.gate.auc == auc_sharp
    end

    @testset "no equivalence machinery (D-14)" begin
        # The traffic light is a BAND, not a p-value. There is no null hypothesis, so there is no
        # power to be over-supplied and nothing an equivalence test could add -- unlike Phase 11's
        # D-08, where a point-null WAS being tested. Asserted as an ABSENCE across every Phase-13
        # source file, comments stripped, so a later plan cannot quietly introduce one.
        #
        # THE ONE DOCUMENTED EXEMPTION: `P13_TOST_REQUIRED = false` in the frozen consts file. That
        # binding IS the pre-registered record that no such machinery is required, so it is removed
        # from the text before the grep rather than allowed to satisfy it.
        @test P13_TOST_REQUIRED == false
        @test length(P13_SRC_FILES) >= 7          # the grep is not scanning an empty directory
        for (name, code) in P13_ALL_CODE
            scrubbed = replace(code, "P13_TOST_REQUIRED" => "")
            @test !occursin("tost", lowercase(scrubbed))
            @test !occursin("equivalence", lowercase(scrubbed))
        end
        # ... and result.jl in particular carries none of it.
        @test haskey(P13_ALL_CODE, "result.jl")
    end

    @testset "minimum evaluation count is pre-registered" begin
        # ECE is biased DOWNWARD at small n and when bins are sparse, so a green ECE on a small
        # evaluation set is not evidence of calibration -- it is evidence that the bins were not
        # populated enough to disagree. The minimum per-head count is therefore part of the frozen
        # pre-registration, and P13_GATE_M was sized so each head's RESTRICTED evaluation set
        # clears it with margin.
        @test P13_MIN_EVAL_PER_HEAD == 1000
        @test P13_GATE_M > 2 * P13_MIN_EVAL_PER_HEAD
    end

    @testset "calibration ran CPU-only (D-04)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
