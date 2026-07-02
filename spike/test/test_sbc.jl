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

# spike/test/test_sbc.jl --- Phase-5 SBC calibration gate (SBC-01..04).
#
# The FAST per-task gate: it exercises the SBC machinery on a SMALL fixture
# (SBC_FIX_M/SBC_FIX_L at a 64×64 imsize) drawn from the FIXTURE stream
# val_rng(VAL_FIX_SEED) -- so the quick gate NEVER consumes the reserved reported
# VAL_MASTER_SEED stream (D-02, keeping the reported run never-tuned-against). The
# REPORTED M=SBC_M / L=SBC_L run is owned by Wave-3 (05-04 run_sbc.jl), not here.
#
# Mirrors test_npe.jl: license header → using Test → guarded includes of the units
# under test → const fixture computed ONCE before the gate → one outer @testset.
# NO figure call appears here (figures are script artifacts, not gate assertions).

using Test
using Random

# Unit under test: the SBC diagnostics (pulls harness.jl + consts.jl transitively).
isdefined(@__MODULE__, :sbc_ranks) || include(joinpath(@__DIR__, "..", "validation", "sbc.jl"))

# --- Fixture computed ONCE (small, 64×64, FIXTURE stream — never the reported one) --
const SBC_FIX_M_MODEL = load_frozen_model()
const SBC_FIX_RANKS   = sbc_ranks(SBC_FIX_M_MODEL; M = SBC_FIX_M, L = SBC_FIX_L,
                                  imsize = (64, 64), rng = val_rng(VAL_FIX_SEED))

@testset "Phase 5 — SBC Calibration (SBC-01..04)" verbose = true begin

    @testset "SC1 (SBC-01) rank table shape + range" begin
        # M×8 integer table: 7 θ marginal ranks + 1 dedicated paired Δρ rank (D-01),
        # every rank in 0:L. The paired Δρ column is computed from differenced draws,
        # not inherited from the ρ_true marginal ranks.
        @test size(SBC_FIX_RANKS) == (SBC_FIX_M, 8)
        @test eltype(SBC_FIX_RANKS) == Int
        @test all(0 .<= SBC_FIX_RANKS .<= SBC_FIX_L)
        # the Δρ column is a genuine independent column, not a copy of the ρ_true column
        @test SBC_FIX_RANKS[:, 8] != SBC_FIX_RANKS[:, 1]
    end

    @testset "SC2 (SBC-02) uniformity + coverage + ECE traffic-light" begin
        # (a) uniformity p-values are finite, in [0,1], on the fixture ρ_true ranks.
        uρ = sbc_uniformity(SBC_FIX_RANKS[:, 1]; L = SBC_FIX_L, bins = SBC_FIX_BINS)
        @test isfinite(uρ.ks_p)   && 0.0 <= uρ.ks_p   <= 1.0
        @test isfinite(uρ.chi2_p) && 0.0 <= uρ.chi2_p <= 1.0

        # (b) a DELIBERATELY biased rank vector (all zeros) drives both p ≈ 0 below the
        #     pre-registered alphas — the test exercises BOTH directions (calibrated vs
        #     miscalibrated), so a stuck "always-pass" p-value would be caught.
        ub = sbc_uniformity(zeros(Int, SBC_FIX_M); L = SBC_FIX_L, bins = SBC_FIX_BINS)
        @test ub.ks_p   < SBC_KS_ALPHA
        @test ub.chi2_p < SBC_CHI2_ALPHA

        # (c) coverage curve: empirical coverages are valid probabilities and the curve
        #     is monotone non-decreasing in the nominal level.
        cov = sbc_coverage(SBC_FIX_RANKS[:, 1]; L = SBC_FIX_L)
        @test length(cov.nominal) == length(cov.empirical)
        @test all(0.0 .<= cov.empirical .<= 1.0)
        @test issorted(cov.empirical)

        # (d) _bin_calibration on PERFECTLY-calibrated synthetic input → ECE ≈ 0 (green);
        #     on MAXIMALLY-miscalibrated input → ECE near its max (red). Traffic-light
        #     returns a Symbol on both.
        cal_good = _bin_calibration(fill(0.5, 100),
                                    vcat(trues(50), falses(50)); n_bins = 10)
        @test cal_good.ece < SBC_ECE_GREEN
        @test sbc_traffic_light(cal_good.ece) === :green
        @test sbc_traffic_light(cal_good.ece) isa Symbol

        cal_bad = _bin_calibration(fill(1.0, 100), falses(100); n_bins = 10)
        @test cal_bad.ece > SBC_ECE_YELLOW
        @test sbc_traffic_light(cal_bad.ece) === :red

        # (e) the SBC-routed calibration returns a finite-ECE CalibrationResult and a
        #     Symbol verdict on the fixture ranks.
        cal_sbc = sbc_calibration(SBC_FIX_RANKS[:, 1]; L = SBC_FIX_L, n_bins = SBC_FIX_BINS)
        @test cal_sbc isa CalibrationResult
        @test isfinite(cal_sbc.ece) && isfinite(cal_sbc.mce)
        @test sbc_traffic_light(cal_sbc.ece) isa Symbol
    end

    @testset "SC3 (SBC-03) pre-registration + disjoint seed (anti-snooping)" begin
        # The reported stream MUST be disjoint from the training/holdout stream (D-02).
        @test VAL_MASTER_SEED != NPE_MASTER_SEED
        @test VAL_MASTER_SEED == 0x5BC0FFEE     # the pre-registered locked value
        @test VAL_FIX_SEED    != VAL_MASTER_SEED # fixtures never consume the reported stream
        # the pass/fail thresholds are the locked pre-registered values.
        @test SBC_M == 2000
        @test SBC_L == 999
        @test (SBC_L + 1) % SBC_BINS == 0        # rank-bin evenness (Pitfall 4)
    end

    @testset "SC4 (SBC-04) caption pairs calibration with OOD" begin
        @test occursin("under the simulator", SBC_CAPTION)
        @test occursin("OOD", SBC_CAPTION)
    end

    @testset "SBC ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
