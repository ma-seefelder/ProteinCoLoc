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

# spike/test/test_p13_consts.jl --- the Phase-13 Tier-1 pre-registration gate (D-01, D-04).
#
# WHAT THIS FILE IS FOR. `spike/p13/consts.jl` is the scientific-integrity artifact of
# Phase 13: every threshold, ladder, seed, statistic and design choice that a reported
# Phase-13 number will be scored against. Its value is entirely in the fact that it was
# committed FIRST. This suite asserts every locked value AS A LITERAL, so a post-hoc
# edit -- the one failure mode that would silently invalidate every downstream number --
# breaks the suite loudly instead of passing quietly.
#
# It also re-runs the D-01 seed-disjointness proof against the RECOMPUTED ship-gate seed
# families (a comment naming a derived seed is not evidence), and proves the Tier-2 tau
# slot is still EMPTY.
#
# Mirrors test_sbc.jl / test_bf.jl: license header -> using Test -> guarded include of
# the unit under test -> source read ONCE at module level -> one outer @testset. No
# training, no simulation, no figure call -- this gate is milliseconds.

using Test

# Unit under test: the Tier-1 pre-registration. Guarded so a re-include under the harness
# is a silent no-op.
isdefined(@__MODULE__, :P13_DEV_SEED) || include(joinpath(@__DIR__, "..", "p13", "consts.jl"))

# The pre-registration source, read ONCE (the test_bf.jl:67 idiom), for the source-grep
# assertions below. COMMENT LINES ARE STRIPPED before any count/occursin assertion, so the
# file's own explanatory prose -- which necessarily NAMES the forbidden forms in order to
# forbid them -- cannot invalidate its own gate.
const P13_CONSTS_SRC = read(joinpath(@__DIR__, "..", "p13", "consts.jl"), String)
const P13_CONSTS_CODE = join(
    filter(l -> !startswith(strip(l), "#"), split(P13_CONSTS_SRC, '\n')), '\n')

@testset "P13 Tier-1 pre-registration (D-01, D-04)" verbose = true begin

    @testset "seeds are fresh and disjoint" begin
        # The locked literals. A changed seed is a changed experiment.
        @test P13_DEV_SEED == 0x0000_0000_0B13_DE71
        @test P13_FIX_SEED == 0x0000_0000_0B13_F1F7
        @test P13_SALT     == 0x2545_F491_4F6C_DD1D

        # D-01: neither fresh stream is any reserved, burned or RECOMPUTED gate seed.
        @test !(UInt64(P13_DEV_SEED) in _p13_forbidden())
        @test !(UInt64(P13_FIX_SEED) in _p13_forbidden())

        # Pairwise, by name, so a failure says WHICH stream was collided with.
        for reserved in (VAL_MASTER_SEED, NPE_MASTER_SEED, VAL_FIX_SEED, RATIO_PAIR_SEED,
                         CORPUS_MASTER_SEED, DEFAULT_MASTER_SEED, P11_DEV_SEED)
            @test UInt64(P13_DEV_SEED) != UInt64(reserved)
            @test UInt64(P13_FIX_SEED) != UInt64(reserved)
        end
        @test UInt64(P13_DEV_SEED) != UInt64(P13_FIX_SEED)

        # A repeated (seed, salt) pair is a repeated Philox key, hence a repeated stream.
        @test !(P13_SALT in P13_REPO_SALTS)

        # PROOF THAT THE RECOMPUTE PATH RAN, not that a comment claims it did: both frozen
        # gate families are DERIVED, and at least one member of each must be in the set.
        @test any(v -> UInt64(v) in _p13_forbidden(), values(_GC.PROD_SEED))
        @test any(v -> UInt64(v) in _p13_forbidden(), values(_GC.PROD_SEED_V2))
        # A duplicate would silently shrink the forbidden set.
        @test length(_p13_forbidden()) == length(_p13_forbidden_list())
    end

    @testset "streams are counter-separated" begin
        # Distinct reserved counters must not share a sub-stream.
        @test rand(p13_rng(P13_TAU_COUNTER), UInt64) != rand(p13_rng(P13_DATAGEN_COUNTER), UInt64)
        # The fixture stream never overlaps the reported one, even at the same counter.
        @test rand(p13_fix_rng(P13_FIXTURE_COUNTER), UInt64) !=
              rand(p13_rng(P13_FIXTURE_COUNTER), UInt64)
        # Reproducibility: same construction, same draw.
        @test rand(p13_rng(P13_TAU_COUNTER), UInt64) == rand(p13_rng(P13_TAU_COUNTER), UInt64)
        # Fixtures must never ride a reported counter.
        @test P13_FIXTURE_COUNTER ∉ (P13_TAU_COUNTER, P13_DATAGEN_COUNTER, P13_GATE_COUNTER,
                                     P13_ALPHA_COUNTER, P13_CONTINUITY_COUNTER)
    end

    @testset "D-05 cut variant is locked" begin
        # Variant (a): one tau governs the level AND the control-contrast factor. Switching
        # variants would relabel every training example, so it is locked, not inferred.
        @test P13_CUT_VARIANT === :tau_contrast
    end

    @testset "D-06 probe spec is locked" begin
        @test P13_TAU_DELTA_GRID == (0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20)
        @test P13_TAU_AUC == 0.90
        @test P13_TAU_R == 400
        @test P13_TAU_STATISTIC === :mbar_mask_weighted
        # The single most important difference from the Phase-11 D-06 probe: a PAIRED design
        # would measure a counterfactual the deployed tool never has and understate tau.
        @test P13_TAU_DESIGN === :unpaired
        @test P13_TAU_BOTH_DIRECTIONS == true
        @test P13_TAU_BOOTSTRAP_B == 1000
        # The abort criterion must be immune to its own result.
        @test P13_TAU_ABORT_EXTEND_GRID == false
        @test P13_TAU_REQUIRE_POST_P11_SIMULATOR == true
        # A RULE plus an expected value; the realized number is asserted against Phase-11's
        # LAMBDA_MAX at probe time (plan 13-10).
        @test P13_TAU_REFERENCE_LAMBDA_RULE === :widest_rung
        @test P13_TAU_REFERENCE_LAMBDA_EXPECTED == 3.0
    end

    @testset "tau is measured and locked (Tier 2)" begin
        # FLIPPED BY PLAN 13-10, which is the plan the previous version of this testset named
        # as the one that would flip it. Until 13-10 ran, tau NOT existing was the
        # pre-registration guarantee and was asserted here. 13-10 made the measurement in two
        # commits, in this order: 581602a (the runner plus the persisted A(delta) curve, with
        # no Tier-2 value anywhere) and then the Tier-2 append to spike/p13/consts.jl. This
        # edit belongs to the SECOND of those two commits, so the guarantee that survives is
        # now the ORDERING, visible in git history, rather than tau's absence.
        @test isdefined(@__MODULE__, :P13_TAU)
        # A grid point, never an interpolated or rounded value: tau is read OFF the frozen
        # grid, so a tau that is not on it would mean the reading rule was bypassed.
        @test P13_TAU in P13_TAU_DELTA_GRID
        # The measured discriminability actually cleared the frozen bar.
        @test P13_TAU_MEASURED_AUC >= P13_TAU_AUC
        # Both provenance shas are full 40-character hex object names, so the simulator the
        # probe measured against and the tree it ran on are both pinned rather than described.
        @test occursin(r"^[0-9a-f]{40}$", P13_TAU_PROBE_SHA)
        @test occursin(r"^[0-9a-f]{40}$", P13_TAU_SIMULATOR_SHA)
        # tau is really tau-of-lambda; the reference must travel with the number.
        @test P13_TAU_REFERENCE_LAMBDA == P13_TAU_REFERENCE_LAMBDA_EXPECTED
        @test P13_TAU_PROBE_ARTIFACT == "spike/p13/tau_probe_report.jld2"
        # The Tier-1 accessor now resolves instead of throwing -- same function, no edit.
        @test p13_tau() === P13_TAU
    end

    @testset "D-07 design and bars are locked" begin
        @test P13_STRATIFICATION === :class_frequency
        @test P13_STRATIFICATION_FALLBACK === :importance_weighted
        @test P13_TARGET_CLASS_FREQ == (exclusion = 1//3, random = 1//3, coloc = 1//3)
        # ONE documented iteration, authorised in advance, with ONE named trigger.
        @test P13_ITERATION_ALLOWANCE == 1
        @test P13_ITERATION_TRIGGER isa AbstractString
        @test !isempty(strip(P13_ITERATION_TRIGGER))
        @test occursin("high", lowercase(P13_ITERATION_TRIGGER))
        # F5 verification pass bars, and the control without which the test cannot tell
        # "the correction works" from "the correction is unnecessary here".
        @test P13_F5_CORR_MIN == 0.99
        @test P13_F5_MAXABS_TOL == 0.25
        @test P13_F5_CENTRAL_FRAC == 0.90
        @test P13_F5_NEGCTRL_REQUIRED == true
        @test P13_F5_SKEW_FREQ == (0.60, 0.25, 0.15)
    end

    @testset "D-09/D-10 recipe is copied verbatim" begin
        # Copied, never chosen: this is what D-10's "difference attributable to the head
        # alone" argument rests on. Any drift here silently breaks the attribution.
        @test P13_SUMMARY_WIDTH == 256
        @test P13_NUM_SUMMARIES == 64
        @test P13_TRAIN_N == 48_000
        @test P13_EPOCHS == 300
        @test P13_BATCHSIZE == 128
        @test P13_LR == 2.5e-4
        @test P13_WEIGHT_DECAY == 1e-4
        @test P13_VAL_FRAC == 0.15
        @test P13_STOPPING_EPOCHS == 40
        @test P13_USE_GPU == false
        @test P13_SMOKE_N == 200
        @test P13_SMOKE_EPOCHS == 2
        @test P13_LAMBDA_PLACEMENT === :append_after_pair_encode
        @test P13_INPUT_WIDTH_RULE === :ratio_input_dim_plus_ncond
        # NEVER 321 AS A LITERAL: Phase 11 may deliver n_cond != 1.
        @test !occursin("321", P13_CONSTS_CODE)
    end

    @testset "D-12/D-13 gate numbers are locked" begin
        @test P13_GATE_M == 4000
        @test P13_AUC_FLOOR_COLOC == 0.90
        @test P13_AUC_FLOOR_EXCLUSION == 0.90
        # The confusion matrix is DESCRIPTIVE; decisions/abstention are Phase 14's scope.
        @test P13_CONFUSION_RULE === :argmax_descriptive
        @test P13_CONTINUITY_GATED == false
        @test P13_ECE_GREEN == 0.05
        @test P13_ECE_YELLOW == 0.10
        @test P13_ECE_NBINS == 10
        # ECE IS THE GATE STATISTIC AND MCE IS EXPLICITLY NOT: _bin_calibration scores an
        # EMPTY bin as |midpoint - 0.0| with ZERO ECE weight but FULL MCE weight, so on a
        # well-separated three-way problem MCE is dominated by empty bins.
        @test P13_GATE_STATISTIC === :ece
        @test P13_MIN_EVAL_PER_HEAD == 1000
        @test P13_VACUOUS_AUC_FLOOR == 0.60
        # D-14: a traffic-light band is not a hypothesis test, so there is no over-power
        # pathology and no equivalence machinery is needed.
        @test P13_TOST_REQUIRED == false
    end

    @testset "D-16 alpha ladder is locked" begin
        @test P13_ALPHA_LADDER == (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0)
        @test P13_ALPHA_BG_QUANTILE == 0.05
        @test P13_ALPHA_MASK_RULE === :calculate_mask_ch1_once
        @test P13_ALPHA_N_IMAGES == 64
        @test P13_ALPHA_INVARIANT_RTOL == 1e-10
        @test P13_ALPHA_GATED == false
        # TWO-VALUED, asserted as a literal, so a silent reduction to one arm breaks the
        # suite rather than quietly turning two experiments into one.
        @test P13_ALPHA_SUBSTRATE == (:simulated, :real)
        # PITFALL 2b: the alternative all-pixels-strictly-positive invariant is FALSE on the
        # committed real fixtures at alpha = 0 (positive_c2 carries 2 source zeros,
        # negative_c2 carries 3, of 1 414 528) and is DELIBERATELY not the locked rule.
        # The locked rule is "alpha introduces no NEW zero".
        @test P13_ALPHA_ZERO_INVARIANT === :count_iszero_preserved
        @test P13_ALPHA_MASK_FRACTION_BOUNDS == (0.01, 0.40)
        @test P13_ALPHA_MAX_VALUE_BOUND == 1.0
    end

    @testset "D-15 real-image arm knobs are locked" begin
        # The three sentinels that make "qualitative only, reported never gated, read-only
        # input" machine-checked rather than merely written.
        @test P13_REAL_IS_GATED == false
        @test P13_REAL_QUALITATIVE_ONLY == true
        @test P13_REAL_READ_ONLY == true

        @test P13_REAL_CONDITIONS == ("positive", "negative")
        @test P13_REAL_SAMPLE == "positive"
        @test P13_REAL_CONTROL == "negative"
        @test P13_REAL_CHANNEL_PAIR == (1, 2)
        @test P13_REAL_REDUNDANCY_PAIR == (1, 3)
        @test P13_REAL_IMSIZE == (1028, 1376)
        @test P13_REAL_GRID_TRUNCATION_ROWS == 4

        # The frozen ghat lineage anchors the ingestion regresses against.
        @test P13_REAL_ANCHOR_MBAR.positive == 0.3292
        @test P13_REAL_ANCHOR_MBAR.negative == 0.2481
        @test P13_REAL_ANCHOR_TOL == 1e-3

        @test P13_REAL_ALPHA_GRID == P13_ALPHA_LADDER
        @test P13_REAL_LAMBDA_READS_RULE === :full_phase11_ladder
        @test P13_REAL_LAMBDA_HEADLINE_RULE === :widest_rung
        @test P13_REAL_LAMBDA_HEADLINE_EXPECTED == 3.0

        @test P13_REAL_OOD_COMPARISON == true
        @test P13_REAL_OOD_SHIPPED_DENSITY == 433.69
        @test P13_REAL_OOD_SHIPPED_THRESHOLD == 179.14

        @test P13_REAL_NAMING_CORRECTION isa AbstractString
        @test !isempty(strip(P13_REAL_NAMING_CORRECTION))
        @test occursin("biological", P13_REAL_NAMING_CORRECTION)
        @test occursin("0.2481", P13_REAL_NAMING_CORRECTION)
        @test P13_REAL_SUBSTITUTION_RECORD isa AbstractString
        @test !isempty(strip(P13_REAL_SUBSTITUTION_RECORD))
        @test occursin("Phase 16", P13_REAL_SUBSTITUTION_RECORD)

        # THE NOT-A-GATE SOURCE ASSERTION. Over the comment-stripped source, collect every
        # P13_REAL_* binding whose name is bar-shaped and assert the set is EXACTLY the two
        # DOCUMENTED EXEMPTIONS:
        #   P13_REAL_ANCHOR_TOL          -- an ingestion REGRESSION tolerance, not a P13 bar
        #   P13_REAL_OOD_SHIPPED_THRESHOLD -- a RECORDED reference value belonging to the
        #                                     SHIPPED net, not a P13 bar
        # SET EQUALITY, not a count: a count would let the exemption list silently grow by
        # one while another entry was removed. The name pattern is \b-anchored on BOTH ends
        # so a bar-shaped SUBSTRING inside a non-bar name (e.g. the "MIN" inside
        # P13_REAL_NAMING_CORRECTION) is not counted as a threshold.
        bar_names = Set(m.match for m in eachmatch(
            r"\bP13_REAL_[A-Z_]*(?:FLOOR|THRESHOLD|MIN|MAX|TOL)\b", P13_CONSTS_CODE))
        @test bar_names == Set(["P13_REAL_ANCHOR_TOL", "P13_REAL_OOD_SHIPPED_THRESHOLD"])
    end

    @testset "F5 imsize mixture was READ not retyped" begin
        # The binding invariant is train-joint == eval-joint. Reading the mixture through the
        # isolated module is what makes drift impossible; retyping it is what makes drift
        # invisible.
        @test P13_IMSIZE_SET == _GC.SBC_IMSIZE_SET
        @test P13_IMSIZE_WEIGHTS == _GC.SBC_IMSIZE_WEIGHTS
        @test sum(P13_IMSIZE_WEIGHTS) ≈ 1.0
        @test length(P13_IMSIZE_SET) == length(P13_IMSIZE_WEIGHTS)
        # The anchor frame size must not appear as a retyped literal outside comments.
        @test !occursin("(1376, 1028)", P13_CONSTS_CODE)
        # ...and the read must actually go through the frozen gate pre-registration.
        @test occursin("gate_consts_8_v2", P13_CONSTS_SRC)
    end

    @testset "deferrals and deviations are on the record" begin
        @test P13_SRC_UNTOUCHED == true
        @test P13_PHYSICAL_ANCHOR_DEFERRED_TO == "Phase 16"
        @test P13_RESULTS_RENAME_DEFERRED == true
        @test length(P13_DECLARED_DEVIATIONS) == 4
        @test all(d -> d isa AbstractString && !isempty(strip(d)), P13_DECLARED_DEVIATIONS)
    end

    @testset "P13 ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
