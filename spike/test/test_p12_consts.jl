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

# spike/test/test_p12_consts.jl --- the Phase-12 Tier-1 pre-registration gate (R-5).
#
# THIS FILE EXISTS TO MAKE A POST-HOC EDIT LOUD. Every locked value is asserted AS A
# LITERAL here, so changing a constant in `p12_consts.jl` after a run breaks the suite
# rather than silently rewriting the pre-registration. The assertions are deliberately
# redundant with the `@assert` self-checks inside the consts file: the self-checks prove
# internal CONSISTENCY, these prove the values are the SPECIFIC values that were committed.
#
# IT ALSO GUARDS THE WIRING. Phase 12 inherits an include-ordering trap -- a testset placed
# after `test_p13_correction.jl` never runs, because that file's outer testset throws -- so
# testset 9 asserts the ORDER of the includes in `runtests.jl` at source level, and asserts
# that every Phase-12 test file has a line in the aggregator. A test that cannot run is worse
# than a test that fails.
#
# DECOUPLING (CLAUDE.md): spike-local; reaches `test/gate/` only READ-ONLY, transitively
# through the isolated `GateV2P12` module the consts file opens.

using Test

# Unit under test: the Tier-1 pre-registration. Guarded for idempotency (S2).
isdefined(@__MODULE__, :P12_DEV_SEED) ||
    include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))

# Source text with whole-line comments stripped, so a filename MENTIONED IN PROSE cannot be
# mistaken for an include. Ported from `spike/test/test_stage6_regression.jl:79-80`.
_strip_comment_lines(src::AbstractString) =
    join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')

# The ten Phase-12 test filenames, in aggregator order, kept in ONE tuple so the count cannot
# drift between the wiring assertion below and `test_p12_suite.jl` itself.
const P12_TEST_FILES = ("test_p12_consts.jl", "test_p12_lattice.jl", "test_p12_prior.jl",
                        "test_p12_architecture.jl", "test_p12_datagen.jl", "test_p12_train.jl",
                        "test_p12_result.jl", "test_p12_sbc.jl", "test_p12_coverage.jl",
                        "test_p12_decoupling.jl")

@testset "P12 Tier-1 pre-registration (R-5)" verbose = true begin

    @testset "Tier-1 seed literals (D-01 / R-5)" begin
        @test P12_DEV_SEED     === 0x0000_0000_0B12_DE71
        @test P12_FIX_SEED     === 0x0000_0000_0B12_F1F7
        @test P12_SALT         === 0xFF51_AFD7_ED55_8CCD
        @test P12_DATAGEN_SALT === 0xC2B2_AE3D_27D4_EB4F
        @test P12_DEV_SEED     isa UInt64
        @test P12_FIX_SEED     isa UInt64
        @test P12_SALT         isa UInt64
        @test P12_DATAGEN_SALT isa UInt64
    end

    @testset "seed disjointness is EXECUTABLE, not documentary" begin
        @test !(UInt64(P12_DEV_SEED) in _p12_forbidden())
        @test !(UInt64(P12_FIX_SEED) in _p12_forbidden())
        # Pairwise, BY NAME, so a failure says WHICH stream collided.
        @test P12_DEV_SEED != NPE_MASTER_SEED
        @test P12_DEV_SEED != VAL_MASTER_SEED
        @test P12_DEV_SEED != VAL_FIX_SEED
        @test P12_DEV_SEED != DEFAULT_MASTER_SEED
        @test P12_DEV_SEED != RATIO_PAIR_SEED
        @test P12_DEV_SEED != CORPUS_MASTER_SEED
        @test P12_DEV_SEED != P11_DEV_SEED
        @test P12_DEV_SEED != P13_DEV_SEED
        @test P12_DEV_SEED != P13_FIX_SEED
        @test P12_FIX_SEED != P13_FIX_SEED
        # The reserved-stream literals themselves are the locked ones (a typo here would
        # silently shrink the forbidden set to a set the DEV seed trivially misses).
        @test NPE_MASTER_SEED     === 0x0000_0000_00C0_FFEE
        @test VAL_MASTER_SEED     === 0x0000_0000_5BC0_FFEE
        @test VAL_FIX_SEED        === 0x0000_0000_00F1_F7ED
        @test DEFAULT_MASTER_SEED === 0x0000_0000_0000_0001
        @test RATIO_PAIR_SEED     === 0x0000_0000_004A_7107
        @test CORPUS_MASTER_SEED  === 0x0000_0000_00C0_5EED
        @test P11_DEV_SEED        === 0x0000_0000_0B11_DE71
        @test P13_DEV_SEED        === 0x0000_0000_0B13_DE71
        @test P13_FIX_SEED        === 0x0000_0000_0B13_F1F7
        # Phase 13's burned-key inventory is carried WHOLE, and the three Phase-11 tuples it
        # subsumes are carried too -- so a reader can check either source.
        @test length(P13_BURNED_DEV_SEEDS) == 25
        @test issubset(Set(UInt64.(F2_DEV_SEEDS)),        Set(UInt64.(P13_BURNED_DEV_SEEDS)))
        @test issubset(Set(UInt64.(SPIKE_DEV_SEEDS)),     Set(UInt64.(P13_BURNED_DEV_SEEDS)))
        @test issubset(Set(UInt64.(P11_RESEARCH_BURNED)), Set(UInt64.(P13_BURNED_DEV_SEEDS)))
        @test all(s -> P12_DEV_SEED != s, P13_BURNED_DEV_SEEDS)
        @test all(s -> P12_FIX_SEED != s, P13_BURNED_DEV_SEEDS)
        # The eight ship-gate seeds are RECOMPUTED, and the recomputation reproduces the frozen
        # file's own values -- so "forbidden" covers the real gate streams, not a comment.
        @test length(PROD_SEED) == 4 && length(PROD_SEED_V2) == 4
        @test PROD_SEED    == GateV2P12.PROD_SEED
        @test PROD_SEED_V2 == GateV2P12.PROD_SEED_V2
        @test all(s -> P12_DEV_SEED != s, values(PROD_SEED))
        @test all(s -> P12_DEV_SEED != s, values(PROD_SEED_V2))
        # No duplicate may hide inside the forbidden list.
        @test length(_p12_forbidden_list()) == length(unique(_p12_forbidden_list()))
        @test length(_p12_forbidden()) == length(_p12_forbidden_list())
    end

    @testset "salt inventory covers the P11_DATAGEN_SALT / P13_SALT collision" begin
        # ONE literal, TWO phases: spike/p13/consts.jl:192 and spike/data/p11_generate.jl:90.
        # Phase 13's own inventory omits it, so nothing in the repository would have caught a
        # third phase reaching for the same value. This assertion is what closes that.
        @test 0x2545_F491_4F6C_DD1D in P12_REPO_SALTS
        @test 0xA24B_AED4_663E_E121 in P12_REPO_SALTS    # P11_SALT, carried by Phase 13 already
        @test P12_SALT         ∉ P12_REPO_SALTS
        @test P12_DATAGEN_SALT ∉ P12_REPO_SALTS
        @test P12_SALT != P12_DATAGEN_SALT
        @test length(unique(P12_REPO_SALTS)) == length(P12_REPO_SALTS)
    end

    @testset "counter streams are behaviourally distinct" begin
        @test P12_FIXTURE_COUNTER == 99
        @test P12_FIXTURE_COUNTER ∉ (P12_SIM02_COUNTER, P12_STAGE1_COUNTER,
                                     P12_ELL_RIDGE_COUNTER, P12_EPS_RIDGE_COUNTER,
                                     P12_MINISPIKE_COUNTER, P12_DATAGEN_COUNTER,
                                     P12_SBC_COUNTER, P12_COVERAGE_COUNTER, P12_GUARDS_COUNTER)
        # Distinct counters ⇒ distinct Philox sub-streams off the same key.
        @test rand(p12_rng(P12_STAGE1_COUNTER), UInt64) != rand(p12_rng(P12_SBC_COUNTER), UInt64)
        # ... the stream is deterministic under the fixed seed (R-5 reproducibility) ...
        @test rand(p12_fix_rng(P12_FIXTURE_COUNTER), UInt64) ==
              rand(p12_fix_rng(P12_FIXTURE_COUNTER), UInt64)
        # ... and the FIXTURE key is a genuinely different key, not merely a different counter.
        @test rand(p12_rng(P12_STAGE1_COUNTER), UInt64) !=
              rand(p12_fix_rng(P12_STAGE1_COUNTER), UInt64)
    end

    @testset "recorded derivations reproduce" begin
        # The INPUTS are asserted as literals first, so the derivation checks below cannot pass
        # by both sides drifting together.
        @test P12_STAGE1_N                       == 25_000
        @test P12_STAGE1_TEST_FRACTION           == 0.25
        @test P12_STAGE1_IMSIZE                  == (512, 512)
        @test P12_STAGE1_RATIO_CEILING           == 0.95
        @test P12_STAGE1_CONTROL_CEILING         == 0.5
        @test P12_STAGE1_MIN_RUNGS               == 1
        @test P12_STAGE1_N_TEST_PER_RUNG_DERIVED == 1250
        @test P12_R1_LADDER                      == (0.05, 0.25, 0.50, 0.75, 0.95)
        @test P12_STAGE2_N_MIN                   == 271
        @test P12_STAGE2_COVERAGE_TOST_DELTA     == 0.03
        @test P12_COVERAGE_NOMINAL               == 0.90
        @test P12_STAGE2_LOGSCORE_MIN            == 0.02
        @test Z_TWO_SIDED_90                     == 1.644854
        @test P12_SBC_M                          == 2000
        @test P12_SBC_L                          == 999
        @test P12_SBC_BINS                       == 50
        @test P12_SBC_TARGET_KS_ALPHA            == 0.01
        @test P12_SBC_TOST_DELTA                 == 0.10
        @test P12_SBC_TOST_ALPHA                 == 0.05
        @test P12_SIM02_W1_TOL_PERREGION         == 0.10
        @test P12_N_PAIRS                        == 50_000
        @test P12_MINISPIKE_N                    == 10_000
        @test P12_DATAGEN_WALLCLOCK_CEILING_MIN  == 150
        @test P12_ITERATION_ALLOWANCE            == 1
        # N_MIN is SIZED with the Wald half-width; the interval REPORTED is Wilson
        # (p11_stats.jl:49-64). 0.49346/√N ≤ δ ⇒ N ≥ 270.6 ⇒ 271.
        @test ceil(Int, (Z_TWO_SIDED_90 * sqrt(0.09) / P12_STAGE2_COVERAGE_TOST_DELTA)^2) ==
              P12_STAGE2_N_MIN
        # The Stage-1 ceiling is 2.5 sampling sd below 1.0 at n_test per rung.
        @test round(1 - 2.5 / sqrt(2 * P12_STAGE1_N_TEST_PER_RUNG_DERIVED); digits = 3) ==
              P12_STAGE1_RATIO_CEILING
        # n_test counts DATASETS, not region-draws: 0.25 · 25 000 / 5 rungs = 1250.
        @test P12_STAGE1_N_TEST_PER_RUNG_DERIVED ==
              round(Int, P12_STAGE1_TEST_FRACTION * P12_STAGE1_N / length(P12_R1_LADDER))
        # The theta row budget: global c0 + 63 deviations + 7 nuisances + the appended r1 row.
        @test P12_D_MINISPIKE == 1 + P12_K_DEV + 7 + 1 == 72
        @test P12_K_DEV == 63
        @test P12_G     == 8
    end

    @testset "gated and reporting-only constants are disjoint (reported != gated)" begin
        @test isempty(intersect(P12_GATING_CONSTANTS, P12_REPORTING_ONLY_CONSTANTS))
        @test :P12_RADIAL_ENERGY_CEILING in P12_REPORTING_ONLY_CONSTANTS
        @test :P12_STAGE1_RATIO_CEILING  in P12_GATING_CONSTANTS
        @test all(s -> isdefined(@__MODULE__, s),
                  (P12_GATING_CONSTANTS..., P12_REPORTING_ONLY_CONSTANTS...))
        # The reporting-only bars are locked literals too -- "not a gate" is not "not frozen".
        @test P12_RADIAL_ENERGY_CEILING   == 0.50
        @test P12_OFFSET_GRID_TOL         == 0.10
        @test P12_VACUOUS_SHRINKAGE_FLOOR == 0.90
        @test length(unique(P12_GATING_CONSTANTS)) == length(P12_GATING_CONSTANTS)
        @test length(unique(P12_REPORTING_ONLY_CONSTANTS)) ==
              length(P12_REPORTING_ONLY_CONSTANTS)
    end

    @testset "declared deviations are the pre-registered five (D-02/D-08/R-3/R-7)" begin
        @test length(P12_RESEARCH_NET_DEVIATIONS) == 5
        @test all(!isempty, P12_RESEARCH_NET_DEVIATIONS)
        @test all(d -> d isa AbstractString, P12_RESEARCH_NET_DEVIATIONS)
    end

    @testset "Tier 2 does not exist yet" begin
        # THESE THREE ASSERTIONS ARE REMOVED IN THE SAME COMMIT THAT APPENDS THE CORRESPONDING
        # TIER-2 BLOCK, AND ONLY THEN -- 12-15 opens :P12_CHOSEN_PRIOR and :P12_N_LOW, 12-16
        # opens :P12_FISHERZ_NEFF. Until then they prove Tier 2 has not been pre-empted, which is
        # what makes the two-tier structure a pre-registration rather than a filing convention.
        # (Same discipline as `test_p11_consts.jl:135-139`.)
        @test !isdefined(@__MODULE__, :P12_CHOSEN_PRIOR)
        @test !isdefined(@__MODULE__, :P12_N_LOW)
        @test !isdefined(@__MODULE__, :P12_FISHERZ_NEFF)
    end

    @testset "Phase-12 wiring manifest — a late include cannot silently never run" begin
        rt        = _strip_comment_lines(read(joinpath(@__DIR__, "runtests.jl"), String))
        suite_src = _strip_comment_lines(read(joinpath(@__DIR__, "test_p12_suite.jl"), String))

        # The aggregator is wired into the suite at all ...
        @test findfirst("test_p12_suite.jl", rt) !== nothing

        # ... and it precedes EVERY OTHER PHASE-TEST INCLUDE, which is the guarantee this
        # testset exists to encode.
        #
        # THE WEAKER ASSERTION THIS REPLACES ASSERTED ONLY "before `test_p13_correction.jl`",
        # AND THAT WAS NOT ENOUGH -- it was true while the aggregator never ran. A thrown
        # `@testset` aborts every remaining include, so "before the file that throws" is
        # contingent on knowing which file throws TODAY, and on 2026-07-29 that was
        # `test_npe.jl` (Phase-4 SC3, `median_speedup > SPEEDUP_GATE`) rather than the Phase-13
        # correction arm the wiring had been reasoned against. FIRST POSITION IS STRUCTURAL: no
        # other phase's failure can mask Phase 12, whichever file throws.
        #
        # The sibling filenames are DERIVED FROM THE SOURCE TEXT, never hardcoded -- a literal
        # list here would go stale the moment another phase adds an include, which is precisely
        # the class of staleness this assertion replaces. The pattern matches the phase-test
        # include form `include(joinpath(@__DIR__, "<name>.jl"))`; the smoke include inside the
        # outer testset uses a different, two-segment path and is deliberately not in scope.
        _inc_pat  = r"include\(joinpath\(@__DIR__, \"([A-Za-z0-9_]+\.jl)\"\)\)"
        _includes = [(m.offset, m.captures[1]) for m in eachmatch(_inc_pat, rt)]
        @test !isempty(_includes)
        @test any(t -> t[2] == "test_p12_suite.jl", _includes)
        _p12_at = first(t[1] for t in _includes if t[2] == "test_p12_suite.jl")
        for (off, name) in _includes
            name == "test_p12_suite.jl" && continue
            # Written as a pair so a regression names WHICH include jumped ahead of us.
            @test (name => _p12_at < off) == (name => true)
        end
        @test _p12_at == minimum(t[1] for t in _includes)
        # Corollary, kept explicit because it is a stated acceptance criterion of plan 12-01:
        # the aggregator is ahead of the throwing Phase-13 correction arm in particular.
        @test _p12_at < first(findfirst("test_p13_correction.jl", rt))

        # Every Phase-12 test file exists AND has a line in the aggregator, so adding a file
        # without wiring it is a failure here rather than a silent omission at run time.
        @test length(P12_TEST_FILES) == 10
        for name in P12_TEST_FILES
            @test isfile(joinpath(@__DIR__, name))
            @test occursin(name, suite_src)
        end
    end

    @testset "model surface ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
