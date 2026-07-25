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

# spike/test/test_p11_consts.jl --- the Phase-11 Tier-1 pre-registration gate (D-04, D-01).
#
# THIS FILE EXISTS TO MAKE A POST-HOC EDIT LOUD. Every locked value is asserted AS A
# LITERAL here, so changing a constant in `p11_consts.jl` after a run breaks the suite
# rather than silently rewriting the pre-registration. The assertions are deliberately
# redundant with the `@assert` self-checks inside the consts file: the self-checks prove
# internal consistency, these prove the values are the SPECIFIC values that were committed.
#
# DECOUPLING (CLAUDE.md): spike-local; reaches `test/gate/` only READ-ONLY, transitively
# through the isolated `GateV2` module the consts file opens.

using Test

# Unit under test: the Tier-1 pre-registration. Guarded for idempotency (S2).
isdefined(@__MODULE__, :P11_DEV_SEED) ||
    include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))

@testset "P11 Tier-1 pre-registration (D-04)" verbose = true begin

    @testset "the DEV seed and salt are the locked literals (D-01)" begin
        @test P11_DEV_SEED === 0x0000_0000_0B11_DE71
        @test P11_SALT     === 0xA24B_AED4_663E_E121
        @test P11_DEV_SEED isa UInt64
        @test P11_SALT     isa UInt64
    end

    @testset "the DEV seed is disjoint from every forbidden stream (D-01)" begin
        @test !(UInt64(P11_DEV_SEED) in _p11_forbidden())
        # Pairwise, by name, so a failure says WHICH stream collided.
        @test P11_DEV_SEED != NPE_MASTER_SEED
        @test P11_DEV_SEED != VAL_MASTER_SEED
        @test P11_DEV_SEED != VAL_FIX_SEED
        @test P11_DEV_SEED != CORPUS_MASTER_SEED
        @test P11_DEV_SEED != DEFAULT_MASTER_SEED
        @test all(s -> P11_DEV_SEED != s, F2_DEV_SEEDS)
        @test all(s -> P11_DEV_SEED != s, SPIKE_DEV_SEEDS)
        @test all(s -> P11_DEV_SEED != s, P11_RESEARCH_BURNED)
        @test all(s -> P11_DEV_SEED != s, values(PROD_SEED))
        @test all(s -> P11_DEV_SEED != s, values(PROD_SEED_V2))
        # The reserved-stream literals themselves are the locked ones (a typo here would
        # silently shrink the forbidden set to a set the DEV seed trivially misses).
        @test VAL_FIX_SEED       === 0x0000_0000_00F1_F7ED   # = 0xF1F7ED (consts.jl:79)
        @test VAL_MASTER_SEED    === 0x0000_0000_5BC0_FFEE
        @test NPE_MASTER_SEED    === 0x0000_0000_00C0_FFEE
        @test CORPUS_MASTER_SEED === 0x0000_0000_00C0_5EED
        # No duplicate may hide inside the forbidden set.
        @test length(_p11_forbidden()) == length(unique(_p11_forbidden()))
        # The eight derived gate seeds are RECOMPUTED, and the recomputation reproduces the
        # frozen file's own values -- so "forbidden" covers the real gate streams.
        @test length(PROD_SEED) == 4 && length(PROD_SEED_V2) == 4
        @test PROD_SEED    == GateV2.PROD_SEED
        @test PROD_SEED_V2 == GateV2.PROD_SEED_V2
        @test isempty(intersect(Set(values(PROD_SEED)), Set(values(PROD_SEED_V2))))
    end

    @testset "the reserved counters keep fixtures off the reported streams" begin
        @test P11_FIXTURE_COUNTER == 99
        @test P11_FIXTURE_COUNTER ∉ (P11_PROBE_COUNTER, P11_DATAGEN_COUNTER,
                                     P11_LADDER_COUNTER, P11_BREAKDOWN_COUNTER,
                                     P11_ATTENUATION_COUNTER, P11_REALIMAGE_COUNTER)
        # Distinct counters ⇒ distinct Philox sub-streams off the same key.
        @test rand(p11_rng(P11_PROBE_COUNTER), UInt64) != rand(p11_rng(P11_LADDER_COUNTER), UInt64)
        # ... and the stream is deterministic under the fixed seed (D-01 reproducibility).
        @test rand(p11_rng(P11_FIXTURE_COUNTER), UInt64) == rand(p11_rng(P11_FIXTURE_COUNTER), UInt64)
    end

    @testset "the SC2 ladder and λ range are the locked values (D-03, D-07)" begin
        @test LAMBDA_MIN == 0.25
        @test LAMBDA_MAX == 3.0
        @test SC2_RUNGS == (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
        @test first(SC2_RUNGS) > 0.0        # the no-λ-zero rule (interpolation-onset trap)
        @test first(SC2_RUNGS) == LAMBDA_MIN
        @test last(SC2_RUNGS)  == LAMBDA_MAX
        @test issorted(SC2_RUNGS)
    end

    @testset "the SC3 beyond-prior ladder is locked and REPORTED-not-gated (D-14)" begin
        @test SC3_SHIFT_RUNGS == (3.5, 4.0, 5.0, 6.0, 8.0)
        @test SC3_EPS_RUNGS   == (0.025, 0.03, 0.04, 0.05)
        @test issorted(SC3_SHIFT_RUNGS)
        @test issorted(SC3_EPS_RUNGS)
        @test first(SC3_SHIFT_RUNGS) > LAMBDA_MAX     # genuinely beyond the training prior
        @test SC3_READ_LAMBDA_PRIMARY   == 3.0
        @test SC3_READ_LAMBDA_SECONDARY == 1.0
        @test SC3_READ_LAMBDA_PRIMARY   in SC2_RUNGS
        @test SC3_READ_LAMBDA_SECONDARY in SC2_RUNGS
    end

    @testset "the equivalence / coverage constants are the locked values (D-07, D-08)" begin
        @test SC2_COVERAGE_NOMINAL == 0.90
        @test SC2_TOST_DELTA       == 0.03
        @test SC2_TOST_ALPHA       == 0.05
        @test P11_HOLM_FWER        == 0.05
        @test PERM_B               == 10_000
        @test N_PER_RUNG           == 500
        @test N_MIN_DERIVED        == 271
        @test N_PER_RUNG >= 271
        @test N_PER_RUNG >= N_MIN_DERIVED
        @test Z_TWO_SIDED_90 ≈ 1.644853627
        # The recorded N_min derivation reproduces: 0.49346/√N ≤ δ ⇒ N ≥ 270.6 ⇒ 271.
        @test ceil(Int, (Z_TWO_SIDED_90 * sqrt(0.09) / SC2_TOST_DELTA)^2) == N_MIN_DERIVED
    end

    @testset "the probe spec and its abort criterion are locked BEFORE the probe (D-06)" begin
        @test P11_PROBE_THETA_BASE_N == 5
        @test P11_PROBE_R            == 32
        @test P11_PROBE_SHIFT_RUNGS  == (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0)
        @test P11_PROBE_EPS_RUNGS    == (0.0, 0.0025, 0.005, 0.01, 0.02)
        @test P11_PROBE_DRHO_RUNGS   == (0.02, 0.05, 0.10, 0.20)
        @test P11_PROBE_METRIC === :paired_l2_rows_1_64
        # The abort criterion must be immune to the probe's own result.
        @test P11_PROBE_S_FLOOR        == 0.9
        @test P11_PROBE_SPAN_FLOOR     == 0.02
        @test SC2_SPEARMAN_ATTENUATION == 0.5
        # Tier 2 is NOT here yet -- it is appended in a second guard block after the probe.
        @test !isdefined(@__MODULE__, :SC2_SPEARMAN_FLOOR)
    end

    @testset "the F5 imsize mixture is READ from the frozen gate file (F5)" begin
        # These literals are the ASSERTION THAT THE READ SUCCEEDED, not a second source of
        # truth: `p11_consts.jl` never types the mixture, it binds GateV2's constants.
        @test P11_IMSIZE_SET == ((512, 512), (1024, 1024), (1376, 1028), (2048, 2048))
        @test P11_IMSIZE_WEIGHTS == (0.40, 0.25, 0.25, 0.10)
        @test P11_IMSIZE_SET     === GateV2.SBC_IMSIZE_SET
        @test P11_IMSIZE_WEIGHTS === GateV2.SBC_IMSIZE_WEIGHTS
        @test sum(P11_IMSIZE_WEIGHTS) ≈ 1.0
        @test P11_PROBE_COMPARABILITY_IMSIZE == (256, 256)
        @test !(P11_PROBE_COMPARABILITY_IMSIZE in P11_IMSIZE_SET)
    end

    @testset "the budget and its ceiling are locked (no silent downgrade)" begin
        @test P11_N_PAIRS == 50_000
        @test P11_DATAGEN_WALLCLOCK_CEILING_MIN == 150
    end

    @testset "the stage-6 regression contract is exact equality (D-10)" begin
        @test P11_STAGE6_EXACT === true
        @test P11_STAGE6_TOLERANCE_FALLBACK == 1e-13
    end

    @testset "the D-04 novelties: one iteration, four named deviations" begin
        @test P11_ITERATION_ALLOWANCE == 1
        @test length(P11_RESEARCH_NET_DEVIATIONS) == 4
        @test all(d -> d isa AbstractString, P11_RESEARCH_NET_DEVIATIONS)
        # Each declared deviation names the thing it changes, so the report can be checked
        # against this list rather than against prose.
        @test any(d -> occursin("chromatic_eps", d), P11_RESEARCH_NET_DEVIATIONS)
        @test any(d -> occursin("SHIFT_PRIOR", d),   P11_RESEARCH_NET_DEVIATIONS)
        @test any(d -> occursin("129", d),           P11_RESEARCH_NET_DEVIATIONS)
        @test any(d -> occursin("BoundedThetaTransform", d), P11_RESEARCH_NET_DEVIATIONS)
    end

    @testset "P11 consts ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
