#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/test_p15_families.jl --- the assertions on the THREE Phase-15 sweep axes.
#
# `test/gate/p15_misspec.jl` adds `spillover`, `registration` and `autofluorescence` to the four
# shipped OOD families so that all seven can be swept through the EXISTING gate machinery. Every
# property the sweep depends on is asserted here, and every one of them runs in seconds — no net
# is loaded, no inference runs, nothing here burns a rung of compute.
#
# WHAT "MATCHED PAIRS" ACTUALLY BUYS, AND WHERE IT STOPS (read before trusting a cross-rung delta).
# The three new generators are PURE θ-OVERRIDES on `ProteinCoLoc.simulate_pair`: at a fixed seed,
# rung k and rung 0 share the prior draw θ*, the three latent Gaussian fields and both Bernoulli
# thinning masks, because `simulate_pair` consumes those in stages 1-2 — BEFORE it ever reads
# `spillover` (stage 4), `autofluorescence` (stage 5) or the shift (stage 6), neither of which
# consumes any rng at all. Testset 3 asserts that as an exact image equality.
#
# It does NOT follow — and it is MEASURABLY FALSE (testset 3, second half) — that the two runs
# leave the stream in the same state. Stage 7 draws `Poisson(GAIN·intensity)`, whose per-pixel rng
# consumption depends on the RATE, which is exactly what the override changed. So the streams
# diverge at the end of the FIRST draw, and `sbc_gate` threads ONE rng through all M draws
# (`harness.jl:188-204`). The matched-pairs guarantee therefore covers draw 1 and no further; from
# draw 2 on, two rungs see different θ* sequences. Cross-rung ECEs are neither independent
# replicates nor exactly paired ones, and a sweep report must not claim either.
#
# ISOLATED MODULE, as everywhere in this phase: `p15_misspec.jl` pulls in `p15_consts.jl`, which
# defines the same const names (`SBC_M`, `prod_seed`, …) `run_gate.jl` has already loaded into the
# test module, and every gate consts file is guarded on `:SBC_M` — so a direct include would be
# silently skipped and every assertion below would read the wrong file.
#
# RUNS BOTH WAYS: standalone (`julia --project=. --threads=auto test/gate/test_p15_families.jl`)
# and `include`d from `test/runtests.jl`. Every path is built from `@__DIR__`.

using Test
using ProteinCoLoc
import Random123
import Statistics: mean

module P15Fam
    include(joinpath(@__DIR__, "p15_misspec.jl"))
end

# =============================================================================================
# Fixtures
# =============================================================================================
#
# FIXTURE SEEDS, NOT GATE SEEDS. Nothing below is a Phase-15 measurement: no bar is scored, no
# rung is decided, no number here reaches a report. These streams are deliberately NOT taken from
# `prod_rng` so that running the test suite can never consume, advance or otherwise entangle the
# pre-registered production stream (D-11). They are plain literals for the same reason a unit-test
# fixture is: reproducibility of the test, nothing more.
const P15FAM_STREAM_SEED = UInt64(0x0000_0000_0F15_0003)   # drives simulate_pair
const P15FAM_THETA_SEED  = UInt64(0x0000_0000_0F15_00C0)   # drives sample_prior

"A fresh Philox stream at `seed` — the generators mutate their rng, so every comparison needs its own."
p15fam_rng(seed) = Random123.Philox4x(UInt64, (UInt64(seed), UInt64(0)))

# 128×128 at G = 8 gives 16 px per patch (≥ 4 is the guard's floor), so the guard passes while each
# image stays ~4× cheaper than the 256² gate size. Small enough that the whole file runs in seconds.
const P15FAM_IMSIZE = (128, 128)
const P15FAM_G      = 8

"The three axes this file covers, paired with the θ override each rung is supposed to apply."
const P15FAM_NEW_AXES = (
    (:spillover,        (θ, k) -> (; spillover = P15Fam.P15_SPILLOVER_LADDER[k])),
    (:registration,     (θ, k) -> (; shift_dx = P15Fam.P15_SHIFT_LADDER[k] / sqrt(2.0),
                                     shift_dy = P15Fam.P15_SHIFT_LADDER[k] / sqrt(2.0))),
    (:autofluorescence, (θ, k) -> (; autofluorescence = P15Fam.P15_AF_LADDER[k])),
)

@testset "Phase-15 sweep axes (p15_misspec.jl)" begin

    # =========================================================================================
    @testset "generator contract (D-02)" begin
        # Every one of the seven axes must satisfy the contract `gate_ood_roc` consumes:
        # `gen(rng, θ; imsize, level, G) -> Vector` of exactly 2 `Matrix{Float64}`, all finite, ≥ 0.
        θ = ProteinCoLoc.sample_prior(p15fam_rng(P15FAM_THETA_SEED))
        @test length(keys(P15Fam.P15_FAMILIES)) == 7

        for axis in keys(P15Fam.P15_FAMILIES)
            gen = P15Fam.P15_FAMILIES[axis]
            img = gen(p15fam_rng(P15FAM_STREAM_SEED), θ;
                      imsize = P15FAM_IMSIZE, level = 1, G = P15FAM_G)
            @test img isa Vector
            @test length(img) == 2
            @test all(c -> c isa Matrix{Float64}, img)
            @test all(c -> size(c) == P15FAM_IMSIZE, img)
            @test all(c -> all(isfinite, c), img)
            @test all(c -> all(≥(0.0), c), img)
        end

        # The imsize guard is not decorative: a 32×32 image gives 4 px per patch on the 8×8 grid,
        # which `patch_summary` cannot summarize. All three new generators must refuse it, and must
        # refuse it BEFORE simulating (the guard is the first statement, including at level 0).
        for (axis, _) in P15FAM_NEW_AXES
            gen = P15Fam.P15_FAMILIES[axis]
            @test_throws ArgumentError gen(p15fam_rng(P15FAM_STREAM_SEED), θ;
                                           imsize = (32, 32), level = 1, G = P15FAM_G)
            @test_throws ArgumentError gen(p15fam_rng(P15FAM_STREAM_SEED), θ;
                                           imsize = (32, 32), level = 0, G = P15FAM_G)
        end
    end

    # =========================================================================================
    @testset "rung-0 identity (D-03)" begin
        # THE RUNG-0 ANCHOR IS ONLY AN ANCHOR IF IT IS THE UNMISSPECIFIED MODEL. Exact equality,
        # no tolerance: level 0 must be byte-for-byte `simulate_pair` on an identically-seeded
        # stream. This is the cheap unit test that catches a whole class of wiring errors (a
        # dropped `level == 0` branch, an off-by-one into the ladder, a stray extra rng draw)
        # before any rung burns compute.
        θ = ProteinCoLoc.sample_prior(p15fam_rng(P15FAM_THETA_SEED))
        ref = ProteinCoLoc.simulate_pair(p15fam_rng(P15FAM_STREAM_SEED), θ; imsize = P15FAM_IMSIZE)

        for (axis, _) in P15FAM_NEW_AXES
            gen = P15Fam.P15_FAMILIES[axis]
            got = gen(p15fam_rng(P15FAM_STREAM_SEED), θ;
                      imsize = P15FAM_IMSIZE, level = 0, G = P15FAM_G)
            @test got[1] == ref[1]
            @test got[2] == ref[2]
        end
    end

    # =========================================================================================
    @testset "pure θ-override, and the measured limit of matched pairs (D-03)" begin
        # PART 1 — what IS true. Each rung must be EXACTLY `simulate_pair` on θ with one field
        # overridden by the frozen ladder, on the same stream. That is the whole content of the
        # matched-pairs design at the level of a single draw: the prior draw, both latent field
        # triples and both thinning masks are shared, because stages 1-2 run before any overridden
        # field is read and stages 4/6 consume no rng at all.
        θ = ProteinCoLoc.sample_prior(p15fam_rng(P15FAM_THETA_SEED))

        for (axis, override) in P15FAM_NEW_AXES
            gen = P15Fam.P15_FAMILIES[axis]
            for k in 1:P15Fam.P15_RUNGS
                got = gen(p15fam_rng(P15FAM_STREAM_SEED), θ;
                          imsize = P15FAM_IMSIZE, level = k, G = P15FAM_G)
                ref = ProteinCoLoc.simulate_pair(p15fam_rng(P15FAM_STREAM_SEED),
                                                 merge(θ, override(θ, k)); imsize = P15FAM_IMSIZE)
                @test got[1] == ref[1]
                @test got[2] == ref[2]
            end
        end

        # The registration ladder is defined on the shift NORM, applied at a fixed 45°. The norm
        # reproduces the ladder value to within one ulp — `r/√2` then `hypot` is not an exact
        # round trip in Float64 (measured: 0.5 → 0.49999999999999994). Stated as a tolerance rather
        # than claimed exact, because a sub-ulp difference in a sub-pixel shift is physically
        # nothing but an `==` here would be a false claim.
        for r in P15Fam.P15_SHIFT_LADDER
            d = r / sqrt(2.0)
            @test hypot(d, d) ≈ r rtol = 4eps(Float64)
        end

        # PART 2 — what is NOT true, MEASURED rather than assumed. The stream state after a
        # generator call is NOT level-invariant: stage 7's `Poisson(GAIN·intensity)` consumes a
        # rate-dependent number of rng values, and the override is precisely a change of rate.
        # These assertions exist so the false claim ("identical rng consumption across levels")
        # cannot silently re-enter, and so the sweep runner (15-05) treats the divergence as a
        # known property rather than discovering it in a report.
        #
        # THE FOUR CONSECUTIVE WORDS ARE DELIBERATE. Philox4x buffers four `UInt64` per counter
        # block, so comparing a single post-draw would miss nothing in principle but reads as if it
        # might; four makes the comparison obviously insensitive to buffer position.
        "The first four words the NEXT draw of an M-draw arm would see, per rung 0..R."
        post_words(gen, θ) = map(0:P15Fam.P15_RUNGS) do k
            r = p15fam_rng(P15FAM_STREAM_SEED)
            gen(r, θ; imsize = P15FAM_IMSIZE, level = k, G = P15FAM_G)
            [rand(r, UInt64) for _ in 1:4]
        end
        diverging_rungs(post) = findall(k -> post[k + 1] != post[1], 1:P15Fam.P15_RUNGS)

        # `registration` and `autofluorescence` diverge at EVERY rung. MEASURED over six fixture
        # θ (this one and five others): 5 of 5 rungs, every time. Both move the rate over the whole
        # frame — the shift changes every pixel's intensity, the offset adds up to 0.8 to it.
        for axis in (:registration, :autofluorescence)
            post = post_words(P15Fam.P15_FAMILIES[axis], θ)
            @info "post-generator stream state by rung (matched pairs end after draw 1)" axis =
                axis diverging_rungs = string(diverging_rungs(post))
            @test diverging_rungs(post) == collect(1:P15Fam.P15_RUNGS)
        end

        # `spillover` is θ-DEPENDENT and is therefore REPORTED, NOT ASSERTED at this θ. It perturbs
        # only channel 1 (`ch1 .+= θ.spillover .* ch2`), so whether the Poisson consumption count
        # actually changes depends on where the resulting rates sit. MEASURED across six fixture θ:
        # 0, 1, 2 and 5 of 5 rungs diverged. Asserting either invariance or divergence here would
        # be asserting a coincidence of one θ.
        post_a = post_words(P15Fam.P15_FAMILIES.spillover, θ)
        @info "post-generator stream state by rung — spillover is θ-dependent, reported not asserted" diverging_rungs =
            string(diverging_rungs(post_a))

        # But the CLAIM "identical rng consumption across levels" must be falsified for all three
        # axes, not only two — so the diverging case is SHOWN on a second fixture θ (measured 5 of
        # 5), rather than left as an assurance that it can happen.
        θ_b = ProteinCoLoc.sample_prior(p15fam_rng(P15FAM_THETA_SEED + 1))
        @test !isempty(diverging_rungs(post_words(P15Fam.P15_FAMILIES.spillover, θ_b)))
    end

    # =========================================================================================
    @testset "monotonicity in level (D-02)" begin
        # A LADDER THAT IS NOT MONOTONE IS NOT A LADDER — an envelope reported against it would be
        # meaningless. Magnitude is measured as the mean absolute deviation from the axis's OWN
        # rung-0 image at the same seed and the same θ, averaged over prior draws.
        #
        # HONESTY ABOUT WHAT THIS NUMBER CONTAINS: it is the deviation of the FINAL image, so it
        # includes stage 7's resampled photon noise. Once the override moves the Poisson rate
        # enough for the draw to decorrelate, a decorrelation component (~0.16 per channel at
        # GAIN = 50) enters and never leaves. That is why rung 1 already sits near 0.25 on the
        # offset axes rather than near 0 — the reported values are TOTAL image change, not the
        # isolated effect of the parameter.
        #
        # 24 DRAWS, NOT 8: at 8 draws the `spillover` rung-1→rung-2 step (0.1 → 0.2, both inside
        # the U(0, 0.2) prior) is smaller than the across-draw variability and the sequence was
        # MEASURED non-monotone at 2 of 3 fixture seeds. At 24 it held at 3 of 3. The test is
        # deterministic either way; the draw count is what decides whether it asserts a real
        # property or a lucky seed.
        nrep = 24
        acc = Dict(axis => zeros(Float64, P15Fam.P15_RUNGS) for (axis, _) in P15FAM_NEW_AXES)

        for rep in 1:nrep
            θ = ProteinCoLoc.sample_prior(p15fam_rng(P15FAM_THETA_SEED + rep))
            seed = P15FAM_STREAM_SEED + rep
            img0 = ProteinCoLoc.simulate_pair(p15fam_rng(seed), θ; imsize = P15FAM_IMSIZE)
            for (axis, _) in P15FAM_NEW_AXES
                gen = P15Fam.P15_FAMILIES[axis]
                for k in 1:P15Fam.P15_RUNGS
                    imgk = gen(p15fam_rng(seed), θ; imsize = P15FAM_IMSIZE, level = k,
                               G = P15FAM_G)
                    acc[axis][k] += mean(abs.(imgk[1] .- img0[1]) .+ abs.(imgk[2] .- img0[2]))
                end
            end
        end

        for (axis, _) in P15FAM_NEW_AXES
            mad = acc[axis] ./ nrep
            # Reported in the log: these ARE the ladders' effective magnitudes, and they belong in
            # the domain map's honesty section, not only in a passing assertion.
            # `string(...)`: the logger elides long vectors, and a ladder reported as "0.2471,
            # 0.26147, ⋮" is not reported.
            @info "ladder effective magnitude (mean |Δ| from rung 0, $(nrep) prior draws)" axis =
                axis mad = join(round.(mad; digits = 5), ", ")
            @test all(k -> mad[k] < mad[k + 1], 1:(P15Fam.P15_RUNGS - 1))
        end
    end
end
