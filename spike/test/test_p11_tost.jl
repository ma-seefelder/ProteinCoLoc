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

# spike/test/test_p11_tost.jl --- the Phase-11 equivalence statistics gate (D-07, D-08).
#
# The load-bearing test in this file is the HOLM DIRECTION one. `sbc_holm_pass` passes when
# Holm rejects NOTHING; the equivalence family passes when Holm rejects EVERYTHING. A verbatim
# copy of the SBC verdict would invert SC2(b) silently and plausibly, so the direction is
# pinned against a HAND-COMPUTED adjusted vector, in both the passing and the failing case.
#
# Every statistic is exercised in BOTH directions (equivalent AND non-equivalent, covered AND
# broken-down), mirroring `test_sbc.jl:63-68`, so a stuck always-pass statistic is caught.

using Test
using Random

# Unit under test. Guarded for idempotency (S2); pulls p11_consts.jl transitively.
isdefined(@__MODULE__, :p11_tost_pass) ||
    include(joinpath(@__DIR__, "..", "validation", "p11_stats.jl"))

@testset "P11 equivalence / monotonicity statistics (D-07, D-08)" verbose = true begin

    @testset "Wilson CI" begin
        ci = wilson_ci(450, 500)                      # p̂ = SC2_COVERAGE_NOMINAL exactly
        @test ci.lo < SC2_COVERAGE_NOMINAL < ci.hi    # the interval brackets nominal
        # The pre-registered N makes equivalence DECIDABLE: the half-width fits in the band.
        @test (ci.hi - ci.lo) / 2 < SC2_TOST_DELTA
        # The knife-edge N_MIN_DERIVED is the boundary case it was derived as.
        ci_min = wilson_ci(round(Int, 0.90 * N_MIN_DERIVED), N_MIN_DERIVED)
        @test (ci_min.hi - ci_min.lo) / 2 <= SC2_TOST_DELTA + 1e-3
        # Degenerate p̂ = 1 stays inside [0, 1] up to floating slop (Wald would not).
        ci1 = wilson_ci(500, 500)
        @test ci1.hi <= 1.0 + 1e-12
        @test ci1.lo > 0.0
        # A low coverage sits entirely below nominal.
        @test wilson_ci(400, 500).hi < SC2_COVERAGE_NOMINAL
    end

    @testset "TOST p-value (both directions)" begin
        # Coverage exactly nominal at the pre-registered N ⇒ equivalence is established.
        @test tost_pvalue(450, 500) < SC2_TOST_ALPHA
        # Coverage 0.80 is a full 10 pp off nominal ⇒ equivalence must NOT be declared.
        @test tost_pvalue(400, 500) > SC2_TOST_ALPHA
        # A small N cannot establish equivalence even at exactly nominal coverage (the test
        # is conservative in the right direction; this is the reason N_PER_RUNG is large).
        @test tost_pvalue(45, 50) > SC2_TOST_ALPHA
        # The degenerate p̂ = 1 guard keeps the statistic finite rather than NaN.
        @test isfinite(tost_pvalue(500, 500))
        @test isfinite(tost_pvalue(0, 500))
    end

    @testset "Holm DIRECTION (the mirror-image trap)" begin
        # HAND COMPUTED for p = [0.001, 0.002, 0.04], m = 3:
        #   ascending order 0.001, 0.002, 0.04 with multipliers 3, 2, 1
        #   0.001·3 = 0.003            running max = 0.003
        #   0.002·2 = 0.004            running max = 0.004
        #   0.040·1 = 0.040            running max = 0.040
        p = [0.001, 0.002, 0.04]
        adj = p11_holm_adjusted(p)
        @test adj ≈ [0.003, 0.004, 0.04]
        # Every adjusted p is ≤ FWER ⇒ Holm rejects EVERY non-equivalence null ⇒ SC2(b) PASSES.
        @test p11_tost_pass(p; fwer = 0.05) == true
        # ... and the SBC reading of the SAME vector is the OPPOSITE verdict. This line is the
        # whole point of the testset: it proves the two directions genuinely disagree here, so
        # a verbatim copy of `sbc_holm_pass` would have silently inverted SC2(b).
        @test all(>(0.05), p11_holm_adjusted(p)) == false

        # A family where one rung is nowhere near equivalent ⇒ SC2(b) must FAIL.
        q = [0.001, 0.002, 0.30]
        @test p11_holm_adjusted(q) ≈ [0.003, 0.004, 0.30]
        @test p11_tost_pass(q; fwer = 0.05) == false

        # Holm is monotone non-decreasing in the sorted order and capped at 1.
        r = p11_holm_adjusted([0.9, 0.8, 0.7])
        @test all(<=(1.0), r)
        @test p11_tost_pass(r) == false
    end

    @testset "agreement with the frozen adjuster" begin
        # `test/gate/sbc.jl` CANNOT be loaded standalone under `--project=spike`: it includes
        # `test/gate/harness.jl:35`, which does `using ProteinCoLoc`, and that package does not
        # resolve in the spike environment (MEASURED — the load raises
        # `ArgumentError: Package ProteinCoLoc not found in current path`). No dependency is
        # added to make it load; instead the agreement is asserted against the FROZEN
        # `holm_adjusted` in `test/gate/gate_consts_8_v2.jl:361-371` — which is the
        # authoritative implementation that `sbc_holm` itself defers to, and which the isolated
        # `GateV2` module already loaded for the imsize mixture.
        @test_skip GateSBC.sbc_holm_adjusted   # unavailable under --project=spike, see above

        pv = rand(Random.Xoshiro(2026), 8) .* 0.2      # a fixed 8-vector, deterministic
        @test p11_holm_adjusted(pv) ≈ GateV2.holm_adjusted(pv)
        @test p11_holm_adjusted([0.001, 0.002, 0.04]) ≈ GateV2.holm_adjusted([0.001, 0.002, 0.04])
        # And the frozen SBC-direction verdict really is the opposite one on that vector.
        @test GateV2.holm_pass([0.001, 0.002, 0.04]; fwer = 0.05) == false
        @test p11_tost_pass([0.001, 0.002, 0.04]; fwer = 0.05)    == true
    end

    @testset "permutation Spearman over rung labels" begin
        rng = p11_rng(P11_FIXTURE_COUNTER)
        lam = repeat(collect(SC2_RUNGS), inner = 20)
        # A width that increases with λ plus modest noise ⇒ strong positive S, tiny p_perm.
        w   = lam .+ 0.05 .* randn(Random.Xoshiro(11), length(lam))
        res = p11_perm_spearman(lam, w; B = 200, rng = rng)
        @test res.S > 0.9
        @test res.p_perm <= 0.05
        # A width INDEPENDENT of λ must not be declared monotone (the other direction).
        w0   = randn(Random.Xoshiro(12), length(lam))
        res0 = p11_perm_spearman(lam, w0; B = 200, rng = p11_rng(P11_FIXTURE_COUNTER))
        @test res0.p_perm > 0.05
        @test 0.0 < res0.p_perm <= 1.0
        @test_throws ArgumentError p11_perm_spearman(lam, w[1:end-1]; B = 10, rng = rng)
    end

    @testset "shrinkage vacuity diagnostic (reporting only)" begin
        # Column 1 informative (posterior sd ≪ prior sd), column 2 vacuous (posterior = prior).
        prior = hcat(randn(Random.Xoshiro(21), 200), randn(Random.Xoshiro(22), 200))
        post  = hcat(fill(0.05, 200), fill(1.0, 200))
        s = p11_shrinkage(post, prior)
        @test length(s.shrinkage) == 2
        @test s.vacuous[1] == false
        @test s.vacuous[2] == true
        @test s.shrinkage[1] < s.shrinkage[2]
    end

    @testset "breakdown_point" begin
        n  = N_PER_RUNG
        # Coverage holds at rungs 1-2 and collapses at rung 3 ⇒ the 3rd magnitude is reported.
        ks = [450, 450, 400, 380, 350]
        @test breakdown_point(ks, n, SC3_SHIFT_RUNGS) == SC3_SHIFT_RUNGS[3]
        # Right-censored: nothing on the ladder drops evidentially below nominal.
        @test breakdown_point(fill(490, 5), n, SC3_SHIFT_RUNGS) === nothing
        # FIRST, not last: a non-monotone curve still reports the first evidential drop.
        @test breakdown_point([450, 400, 450, 380, 350], n, SC3_SHIFT_RUNGS) == SC3_SHIFT_RUNGS[2]
        @test_throws ArgumentError breakdown_point([450, 450], n, SC3_SHIFT_RUNGS)
    end

    @testset "P11 statistics ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
