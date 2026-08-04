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

"The repository root, resolved from this file rather than from the process working directory."
const P15FAM_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# =============================================================================================
# A SELF-CONTAINED SHA-256, and why this file does not `import SHA`
# =============================================================================================
#
# `SHA` is a Julia stdlib and loads fine under `julia --project=.`, but `Pkg.test()` runs the suite
# in a SANDBOX built from `[deps]` plus the `[targets] test` extras, where `import SHA` fails
# outright and aborts the whole suite (observed by plan 15-01). The obvious fix — adding `SHA` to
# `[extras]` — would EDIT `Project.toml`, which is one of the eight files this phase pins as
# frozen. Pinning a file and editing it in the same phase is not a pre-registration.
#
# `test_p15_consts.jl` carries the same ~50 lines for the same reason. The duplication is
# deliberate rather than factored into a shared helper: `runtests.jl`'s Phase-15 wiring is a FIXED
# five-file list of `test_p15_*.jl` names (15-01 built it so later plans never edit that file
# again), so a shared helper would be a sixth file outside the phase's declared surface. The two
# implementations are CROSS-CHECKED against each other below whenever both are loaded, which turns
# the duplication into an independent second opinion on every pinned digest.

const _P15FAM_SHA256_K = UInt32[
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]

_p15fam_rotr(x::UInt32, n::Int) = (x >> n) | (x << (32 - n))

"SHA-256 (FIPS 180-4) as 64 lowercase hex chars. Dependency-free; verified against NIST vectors below."
function _p15fam_sha256(data::AbstractVector{UInt8})
    msg    = Vector{UInt8}(data)
    bitlen = UInt64(length(data)) * 8
    push!(msg, 0x80)
    while length(msg) % 64 != 56
        push!(msg, 0x00)
    end
    for i in 7:-1:0
        push!(msg, UInt8((bitlen >> (8i)) & 0xff))
    end

    h = UInt32[0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
               0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
    w = Vector{UInt32}(undef, 64)
    for base in 1:64:length(msg)
        @inbounds for i in 0:15
            o = base + 4i
            w[i + 1] = (UInt32(msg[o]) << 24) | (UInt32(msg[o + 1]) << 16) |
                       (UInt32(msg[o + 2]) << 8) | UInt32(msg[o + 3])
        end
        @inbounds for i in 17:64
            s0 = _p15fam_rotr(w[i - 15], 7) ⊻ _p15fam_rotr(w[i - 15], 18) ⊻ (w[i - 15] >> 3)
            s1 = _p15fam_rotr(w[i - 2], 17) ⊻ _p15fam_rotr(w[i - 2], 19) ⊻ (w[i - 2] >> 10)
            w[i] = w[i - 16] + s0 + w[i - 7] + s1
        end
        a, b, c, d, e, f, g, hh = h[1], h[2], h[3], h[4], h[5], h[6], h[7], h[8]
        @inbounds for i in 1:64
            S1  = _p15fam_rotr(e, 6) ⊻ _p15fam_rotr(e, 11) ⊻ _p15fam_rotr(e, 25)
            ch  = (e & f) ⊻ (~e & g)
            t1  = hh + S1 + ch + _P15FAM_SHA256_K[i] + w[i]
            S0  = _p15fam_rotr(a, 2) ⊻ _p15fam_rotr(a, 13) ⊻ _p15fam_rotr(a, 22)
            maj = (a & b) ⊻ (a & c) ⊻ (b & c)
            t2  = S0 + maj
            hh, g, f, e = g, f, e, d + t1
            d, c, b, a  = c, b, a, t1 + t2
        end
        h[1] += a; h[2] += b; h[3] += c; h[4] += d
        h[5] += e; h[6] += f; h[7] += g; h[8] += hh
    end
    return join(string(x; base = 16, pad = 8) for x in h)
end

_p15fam_sha256(s::AbstractString) = _p15fam_sha256(Vector{UInt8}(codeunits(s)))

"""
    _p15fam_lf_sha256(path) -> String

SHA-256 of a text file with CRLF newlines normalized to LF — the normalization
`P15_FROZEN_HASH_NORMALIZATION = :lf` records. `.gitattributes` sets `* text=auto`, so the
working-tree bytes of every text file here are CRLF on Windows and LF on Linux; a raw-bytes pin
would fail on the ubuntu-latest CI this phase exists to build, for a reason that has nothing to do
with the file's contents. A real content change is still caught.
"""
_p15fam_lf_sha256(path::AbstractString) =
    _p15fam_sha256(Vector{UInt8}(codeunits(replace(String(read(path)), "\r\n" => "\n"))))

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

    # =========================================================================================
    @testset "OOD_FAMILIES is byte-unchanged (D-02a)" begin
        # WHY THIS IS A TEST AND NOT A CONVENTION. `OOD_FAMILIES` is an INPUT TO A RECORDED
        # VERDICT: `gate_ood_roc` pools the strongest-level scores over ALL families into
        # `combined_auc`, and `ood_gate` returns `passed = auc_ok && all(values(fam_pass))`. Adding
        # the three Phase-15 axes to that table would therefore change the pooled AUC and add three
        # terms to the pass conjunction — so the shipped, recorded `OOD — PASS (pooled AUC 1.0)`
        # would silently become a DIFFERENT number reported under the old one's name. Not
        # "would need re-running": non-comparable, with nothing to signal it.
        @test keys(P15Fam.OOD_FAMILIES) == (:texture, :noise, :optics, :background)
        @test P15Fam.OOD_FAMILIES.texture    === P15Fam.misspec_texture
        @test P15Fam.OOD_FAMILIES.noise      === P15Fam.misspec_noise
        @test P15Fam.OOD_FAMILIES.optics     === P15Fam.misspec_optics
        @test P15Fam.OOD_FAMILIES.background === P15Fam.misspec_background

        # The in-memory table can only be right if the file it comes from is unchanged, so the pin
        # is checked too — the same digest `p15_consts.jl` froze before any Phase-15 code existed.
        # NIST vectors first: a hand-rolled digest that has not been shown correct proves nothing
        # about the file it is applied to. These four exercise the empty, single-block, two-block
        # and three-block padding paths, which is where such an implementation actually goes wrong.
        @test _p15fam_sha256("") ==
              "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        @test _p15fam_sha256("abc") ==
              "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        @test _p15fam_sha256("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
              "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        @test _p15fam_sha256("abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn" *
                             "hijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu") ==
              "cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac45037afee9d1"

        misspec_path = joinpath(P15FAM_REPO_ROOT, P15Fam.P15_FROZEN_FILE_SHA256.misspec.path)
        got  = _p15fam_lf_sha256(misspec_path)
        want = P15Fam.P15_FROZEN_FILE_SHA256.misspec.sha256
        got == want || @error(
            "test/gate/misspec.jl HAS CHANGED. Mutating OOD_FAMILIES silently rewrites the " *
            "shipped, PASSED OOD arm's combined_auc and its pass conjunction, so the recorded " *
            "verdict stops being reproducible from the code that claims it. New axes must ride " *
            "P15_FAMILIES = merge(OOD_FAMILIES, ...) through the existing `families =` keyword.",
            path = misspec_path, expected = want, actual = got)
        @test got == want

        # An INDEPENDENT SECOND OPINION when one is available. `test_p15_consts.jl` carries its own
        # copy of this digest (see the header); when both files are loaded — which is the case
        # under `runtests.jl` — the two implementations are required to agree, so a shared bug in
        # one of them cannot pass unnoticed.
        if isdefined(@__MODULE__, :_p15_lf_sha256)
            @test _p15_lf_sha256(misspec_path) == got
        else
            @info "cross-check against test_p15_consts.jl's digest skipped (standalone run)"
        end
    end

    # =========================================================================================
    @testset "mechanism typing is complete (D-02, D-02a)" begin
        # THE TYPING IS THE FINDING. An in-prior axis that breaks coverage is a TRAINING failure —
        # the net saw those θ and still miscalibrates. An out-of-model axis that breaks is a SCOPE
        # limit, which is what the OOD flag exists to report. An axis with no label would land in
        # the domain map as an unreadable number, so "every axis is typed" is an assertion.
        @test Tuple(keys(P15Fam.P15_FAMILIES)) == P15Fam.P15_AXES   # shipped four FIRST

        for axis in keys(P15Fam.P15_FAMILIES)
            @test haskey(P15Fam.P15_AXIS_MECHANISM, axis)
            @test haskey(P15Fam.P15_PRIOR_BOUNDARY_RUNG, axis)
            @test P15Fam.p15_axis_mechanism(axis) === P15Fam.P15_AXIS_MECHANISM[axis]
            @test P15Fam.p15_axis_boundary_rung(axis) === P15Fam.P15_PRIOR_BOUNDARY_RUNG[axis]

            if P15Fam.p15_axis_mechanism(axis) === :in_prior
                # The boundary must be an EXACT rung, not a value the ladder straddles (D-03).
                @test P15Fam.p15_axis_boundary_rung(axis) == 2
            else
                @test P15Fam.p15_axis_mechanism(axis) === :out_of_model
                # `nothing` means "no prior support to mark", never "unmeasured".
                @test P15Fam.p15_axis_boundary_rung(axis) === nothing
            end
        end

        # Exactly the three axes this file adds are the in-prior ones — which is why they had to be
        # built: without them the envelope would have no in-prior half at all.
        in_prior = filter(a -> P15Fam.p15_axis_mechanism(a) === :in_prior,
                          collect(keys(P15Fam.P15_FAMILIES)))
        @test Set(in_prior) == Set((:spillover, :registration, :autofluorescence))

        # D-02a's reclassification, asserted rather than left in prose: the shipped `background`
        # family composes a gradient × vignette the simulator has no parameter for and its level-1
        # bleed is already 4× the prior maximum, so it cannot express "at the prior edge" at all.
        @test P15Fam.P15_AXIS_MECHANISM.background === :out_of_model
        @test P15Fam.p15_axis_boundary_rung(:background) === nothing

        # The boundary rungs ARE the simulator's prior extrema — the ladders index the frozen
        # tuples rather than re-typing magnitudes, and this is what stops the two drifting apart.
        @test P15Fam.P15_SPILLOVER_LADDER[2] == maximum(ProteinCoLoc.SPILLOVER_PRIOR)
        @test P15Fam.P15_SHIFT_LADDER[2]     == maximum(ProteinCoLoc.SHIFT_PRIOR)
        @test P15Fam.P15_AF_LADDER[2]        == maximum(ProteinCoLoc.AUTOFLUORESCENCE_PRIOR)
    end

    # =========================================================================================
    @testset "misspec_simulator wrapping (D-02)" begin
        # THE POINT OF THIS TESTSET: the new axes must drop into `sbc_gate`'s `sim` slot with NO
        # adapter, because "the sweep is an orchestration loop over existing calls" is the entire
        # reason this phase is affordable. If a new generator needed a bespoke wrapper, the sweep
        # would be new harness code and would need its own correctness argument.
        θ   = ProteinCoLoc.sample_prior(p15fam_rng(P15FAM_THETA_SEED))
        gen = P15Fam.P15_FAMILIES.autofluorescence
        sim = P15Fam.misspec_simulator(gen; level = 2, G = P15FAM_G)

        @test sim isa NamedTuple
        @test keys(sim) == (:sample_prior, :simulate_pair, :build_mci)
        @test sim.sample_prior === ProteinCoLoc.sample_prior   # the prior stays π(θ): the
        @test sim.build_mci    === ProteinCoLoc.build_mci      # misspecification is in the FORWARD
                                                               # model, which is what the term means

        wrapped = sim.simulate_pair(p15fam_rng(P15FAM_STREAM_SEED), θ; imsize = P15FAM_IMSIZE)
        direct  = gen(p15fam_rng(P15FAM_STREAM_SEED), θ;
                      imsize = P15FAM_IMSIZE, level = 2, G = P15FAM_G)
        @test wrapped[1] == direct[1]
        @test wrapped[2] == direct[2]
    end
end
