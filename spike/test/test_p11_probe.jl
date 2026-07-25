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

# spike/test/test_p11_probe.jl --- FIXTURE-SCALE unit test of the D-06 probe's helpers.
#
# THIS IS NOT THE PROBE. It runs at 128x128 with 2 replicates so it stays a seconds-scale unit
# test, mirroring the `SBC_FIX_*` discipline in `spike/validation/consts.jl`: the fixture arm
# exercises the machinery, the RESERVED arm produces the reported number. Nothing asserted here is
# a measurement, a threshold, or evidence about any net.
#
# THE INCLUDE BELOW MUST NOT TRIGGER THE REPORTED RUN. `P11_PROBE_LOAD_ONLY` is defined first, so
# `run_p11_probe.jl` loads its helpers and stops short of `main()`. Without that guard this file
# would re-consume the reserved probe stream (T-11-21) and take minutes.
#
# DECOUPLING (CLAUDE.md): spike-local; reaches src/ only READ-ONLY through `spike/contract.jl`.

using Test
using Statistics

const P11_PROBE_LOAD_ONLY = true
isdefined(@__MODULE__, :_delta_norm) ||
    include(joinpath(@__DIR__, "..", "validation", "run_p11_probe.jl"))

# Fixture-scale knobs. NOT pre-registered values; they exist only to keep this file fast.
const PROBE_FIX_IMSIZE = (128, 128)
const PROBE_FIX_R      = 2
const PROBE_FIX_THETA  = (ρ_true = 0.4, spillover = 0.1, autofluorescence = 0.05,
                          label_efficiency = 0.8, shift_dx = 0.0, shift_dy = 0.0,
                          noise = 0.5, chromatic_eps = 0.0)

@testset "P11 D-06 probe helpers (fixture scale)" verbose = true begin

    @testset "dRho_eq exactly inverts the fitted line (the headline unit)" begin
        slope, intercept = 6.48, 0.041
        drho = [0.02, 0.05, 0.10, 0.20]
        norms = slope .* drho .+ intercept
        for (d, n) in zip(drho, norms)
            @test isapprox(_drho_eq(n, slope, intercept), d; atol = 1e-12)
        end
        # ... and the fit that supplies those coefficients recovers them exactly on noiseless data.
        fit = _ols(drho, norms)
        @test isapprox(fit.slope,     slope;     atol = 1e-10)
        @test isapprox(fit.intercept, intercept; atol = 1e-10)
        @test isapprox(fit.r2,        1.0;       atol = 1e-12)
        # Round trip through the pair, which is how the probe actually uses them.
        @test isapprox(_drho_eq(fit.slope * 0.137 + fit.intercept, fit.slope, fit.intercept),
                       0.137; atol = 1e-10)
        # BOTH DIRECTIONS: a fit with genuine scatter must NOT report R^2 = 1.
        noisy = _ols(drho, norms .+ [0.3, -0.3, 0.3, -0.3])
        @test noisy.r2 < 0.99
    end

    @testset "the metric reads rows 1:64 ONLY (mask rows cannot move it)" begin
        base = vcat(zeros(64), ones(64))
        # Rows 65:128 differ wildly; the metric must be exactly 0.
        wild = vcat(zeros(64), fill(1.0e6, 64))
        @test _delta_norm(wild, base) == 0.0
        # A single continuous row moving by 3 gives exactly 3.
        one_row = copy(base); one_row[1] += 3.0
        @test _delta_norm(one_row, base) == 3.0
        # ... and it is a Euclidean norm over the continuous block, not a sum or a max.
        two_rows = copy(base); two_rows[1] += 3.0; two_rows[2] += 4.0
        @test _delta_norm(two_rows, base) == 5.0
    end

    @testset "diagonal injection preserves the rung label as the true displacement" begin
        for m in (0.0, 0.25, 1.5, 3.0, 8.0)
            dx, dy = _diag_shift(m)
            @test isapprox(sqrt(dx^2 + dy^2), m; atol = 1e-12)
            @test dx == dy
        end
    end

    @testset "the two key families are disjoint and reproducible" begin
        pk = _pair_key.(1:P11_PROBE_R)
        ik = _imsize_key.(1:P11_PROBE_R)
        @test length(unique(pk)) == P11_PROBE_R
        @test isempty(intersect(Set(pk), Set(ik)))
        @test P11_PROBE_COUNTER ∉ pk
        @test _pair_key(7) == _pair_key(7)          # pure function of the replicate index
    end

    @testset "the imsize sampler respects the F5 weights (Monte-Carlo tolerance)" begin
        n = 4000
        rng = p11_rng(P11_FIXTURE_COUNTER)   # FIXTURES ONLY -- never a reported counter
        draws = [_draw_imsize(rng) for _ in 1:n]
        for (i, sz) in enumerate(P11_IMSIZE_SET)
            @test isapprox(count(==(sz), draws) / n, P11_IMSIZE_WEIGHTS[i]; atol = 0.03)
        end
        @test all(d -> d in P11_IMSIZE_SET, draws)
    end

    @testset "the paired key re-derivation is BIT-IDENTICAL (the design's load-bearing claim)" begin
        k = _pair_key(1)
        a = simulate_pair(p11_rng(k), PROBE_FIX_THETA; imsize = PROBE_FIX_IMSIZE)
        b = simulate_pair(p11_rng(k), PROBE_FIX_THETA; imsize = PROBE_FIX_IMSIZE)
        @test a == b                     # exact, not isapprox
        @test a[1] == b[1] && a[2] == b[2]
        # The whole summary chain therefore reproduces exactly, so the zero rung is exactly zero.
        s1 = _summary(k, PROBE_FIX_THETA; imsize = PROBE_FIX_IMSIZE)
        s2 = _summary(k, merge(PROBE_FIX_THETA, (shift_dx = 0.0, shift_dy = 0.0)); imsize = PROBE_FIX_IMSIZE)
        @test _delta_norm(s2, s1) == 0.0
        # A DIFFERENT key on the same theta does NOT reproduce -- so the equality above is a
        # property of the pairing, not of a degenerate simulator.
        s3 = _summary(_pair_key(2), PROBE_FIX_THETA; imsize = PROBE_FIX_IMSIZE)
        @test _delta_norm(s3, s1) > 0.0
    end

    @testset "pairing beats the between-key floor (why the design is paired at all)" begin
        # At R = 2 this is an existence check, not a measurement: the paired displacement under a
        # real misalignment must be SMALLER than the displacement two independent keys produce at
        # IDENTICAL theta. If it were not, the probe would be reporting Monte-Carlo noise.
        dx, dy = _diag_shift(3.0)
        paired   = Float64[]
        unpaired = Float64[]
        for r in 1:PROBE_FIX_R
            k  = _pair_key(r)
            s0 = _summary(k, PROBE_FIX_THETA; imsize = PROBE_FIX_IMSIZE)
            sp = _summary(k, merge(PROBE_FIX_THETA, (shift_dx = dx, shift_dy = dy));
                          imsize = PROBE_FIX_IMSIZE)
            su = _summary(_imsize_key(r), PROBE_FIX_THETA; imsize = PROBE_FIX_IMSIZE)
            push!(paired,   _delta_norm(sp, s0))
            push!(unpaired, _delta_norm(su, s0))
        end
        @test all(>(0.0), paired)                       # the effect is present
        @test mean(paired) < mean(unpaired)             # and it is below the between-key floor
    end

    @testset "the chromatic axis moves the summary on its own (SEPARATE from registration)" begin
        k  = _pair_key(3)
        s0 = _summary(k, PROBE_FIX_THETA; imsize = PROBE_FIX_IMSIZE)
        se = _summary(k, merge(PROBE_FIX_THETA, (chromatic_eps = 0.02,)); imsize = PROBE_FIX_IMSIZE)
        @test _delta_norm(se, s0) > 0.0
        # chromatic_eps = 0 is the reference itself, so it must land exactly on zero.
        sz = _summary(k, merge(PROBE_FIX_THETA, (chromatic_eps = 0.0,)); imsize = PROBE_FIX_IMSIZE)
        @test _delta_norm(sz, s0) == 0.0
    end

    @testset "the rung sets are the pre-registered ones, unioned only where required" begin
        # The shift arm must contain every pre-registered probe rung, every extension rung AND
        # every SC2 rung -- the last because S_probe is defined across SC2_RUNGS.
        @test all(x -> x in PROBE_SHIFT_ARM, P11_PROBE_SHIFT_RUNGS)
        @test all(x -> x in PROBE_SHIFT_ARM, P11_PROBE_SHIFT_EXT)
        @test all(x -> x in PROBE_SHIFT_ARM, SC2_RUNGS)
        @test all(x -> x in PROBE_EPS_ARM,   P11_PROBE_EPS_RUNGS)
        @test all(x -> x in PROBE_EPS_ARM,   P11_PROBE_EPS_EXT)
        @test PROBE_DRHO_ARM == Float64.(P11_PROBE_DRHO_RUNGS)
        @test issorted(PROBE_SHIFT_ARM) && issorted(PROBE_EPS_ARM)
        @test first(PROBE_SHIFT_ARM) == 0.0 && first(PROBE_EPS_ARM) == 0.0
        # Nothing was invented: the arms are exactly the union, no extra rung.
        @test length(PROBE_SHIFT_ARM) ==
              length(unique(Float64[P11_PROBE_SHIFT_RUNGS..., P11_PROBE_SHIFT_EXT..., SC2_RUNGS...]))
    end

    @testset "the probe cannot write the pre-registration (T-11-19)" begin
        src = read(joinpath(@__DIR__, "..", "validation", "run_p11_probe.jl"), String)
        # Booleans are computed BEFORE @test so a failure prints `false` instead of dumping the
        # whole source file into the test output (11-04 deviation: a naive occursin does that).
        writes_consts = occursin(r"(open|write|jldsave)\([^\n]*p11_consts", src)
        is_a_gate     = occursin(string("@", "testset"), src)
        prior_draws   = count(x -> occursin("sample_prior", x), split(src, '\n'))
        @test writes_consts == false     # no write path from the probe to the pre-registration
        @test is_a_gate     == false     # the probe produces thresholds; it never gates on them
        @test prior_draws   <= 1         # exactly the one frozen theta-base draw, nowhere else
        # The report path is the only artifact it names for writing.
        @test normpath(P11_PROBE_REPORT_PATH) ==
              normpath(joinpath(@__DIR__, "..", "validation", "p11_probe_report.jld2"))
    end

    @testset "probe helpers ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
