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

# spike/test/test_p13_preconditions.jl --- the Phase-11 binding gate (D-02, D-03).
#
# THE ABSENT-ARTIFACT PATH IS THE IMPORTANT TEST, and it is deliberately the first testset.
# D-02 says Phase 13 trains on Phase 11's registration-aware FROZEN summary transform, not on the
# shipped grid-8 one. The catastrophic failure is not a crash; it is a Phase-13 net that trains
# happily on the WRONG BASIS and reports numbers that mean something other than what the phase
# claims. Nothing downstream can detect that from the trained net, so the only place it can be
# caught is the load boundary -- which makes "does the block actually fire, with a message that
# tells the reader what to do" the assertion this whole file exists for.
#
# THIS FILE RUNS IN ANY STATE OF THE REPOSITORY. The block assertion, the no-fallback source
# assertions and the derived-conditioning-width assertions need no artifact at all. Only the
# six-item contract testset needs the Phase-11 net, and when that is absent it SKIPS EXPLICITLY
# with Phase 11 named -- it never fabricates a stand-in artifact, because an artifact invented to
# make a gate green is the same wrong-basis failure wearing a different hat.
#
# Mirrors test_bf.jl / test_p13_net.jl: license header -> using Test -> guarded include of the
# unit under test -> module-level fixtures computed ONCE -> one outer @testset. No image
# simulation, no training and no figure call.

using Test

# Unit under test: the Phase-11 binding surface (pulls p11_architecture.jl -> p13/consts.jl ->
# train_npe.jl -> infer.jl -> p13/net.jl and the simulator/contract chain).
isdefined(@__MODULE__, :p13_require_phase11) ||
    include(joinpath(@__DIR__, "..", "p13", "preconditions.jl"))

# The EXECUTABLE source of the unit under test, with whole-line comments stripped BEFORE any
# count, so a grep assertion cannot be satisfied (or defeated) by prose. This is the same strip
# the plan's acceptance criteria apply (`grep -v '^\s*#'`), reproduced in-suite so the check is
# executable rather than a one-off shell invocation. Docstrings SURVIVE the strip on purpose: a
# forbidden path quoted in a docstring is still a path a reader may copy, so the docstrings in
# preconditions.jl deliberately name none of them.
const P13_PRE_SRC  = read(joinpath(@__DIR__, "..", "p13", "preconditions.jl"), String)
const P13_PRE_CODE = join(filter(l -> !startswith(strip(l), "#"),
                                 split(P13_PRE_SRC, '\n')), '\n')

# A lambda strictly inside the Phase-11 trained range, so the encoding is an interpolation rather
# than the extrapolation the encoder deliberately declines to clamp.
const P13_PRE_LAMBDA = sum(P13_PHASE11_LAMBDA_RANGE) / 2

# The patch grid the 128-row frozen summary implies: 128 = 2*G^2, so G = 8. Written as the
# derivation rather than the number so the assertion below still reads correctly at another grid.
const P13_PRE_G = isqrt(128 ÷ 2)

@testset "P13 Phase-11 preconditions (D-02, D-03)" verbose = true begin

    @testset "absence is a loud, instructional block" begin
        # A path that cannot exist, inside a scratch directory, so this testset's outcome is
        # INDEPENDENT of whether the real Phase-11 artifact is on disk. This is the one Phase-13
        # assertion that is runnable in every state of the repository, and it must stay that way.
        missing_path = joinpath(mktempdir(), "no_such_phase11_net.jld2")
        @test !isfile(missing_path)

        # A HARD error, not a status return and not a skip. `run_gate.jl`'s `:not_trained` shape
        # is right for a RUNNER script that may legitimately decline to score; it is wrong here,
        # because D-02 makes a missing Phase-11 net a block on the phase.
        @test_throws ErrorException p13_require_phase11(missing_path)

        # The message must tell the reader WHAT is blocked, WHERE the decision is recorded, and
        # WHAT NOT TO DO about it. The grid-8 prohibition is the load-bearing sentence: without
        # it, the obvious "fix" for a reader hitting this error is to point the constant at the
        # shipped bundle, which is precisely the failure the block exists to prevent.
        @test_throws "BLOCKED on Phase 11" p13_require_phase11(missing_path)
        @test_throws "13-CONTEXT D-02" p13_require_phase11(missing_path)
        @test_throws "grid-8" p13_require_phase11(missing_path)
        # And it must name the artifact it wanted, so the reader can check the right directory.
        @test_throws "research NPE" p13_require_phase11(missing_path)
    end

    @testset "no fallback path exists" begin
        # A `try` ANYWHERE in this file's executable code would let a caller convert the D-02
        # block into a warning and carry on with the wrong basis -- the exact failure D-02
        # forbids, and the reason this is asserted on the source rather than left to review.
        @test !occursin("try", P13_PRE_CODE)
        # No route to the pre-Phase-11 net, to any shipped bundle directory, or to the grid-8
        # basis. Checked as substrings of the comment-stripped source, so a fallback cannot be
        # added later without turning this red.
        @test !occursin("trained_npe.jld2", P13_PRE_CODE)
        @test !occursin("artifacts/", P13_PRE_CODE)
        @test !occursin("grid_8", P13_PRE_CODE)
        # The bound artifact IS the Phase-11 research net, by name, and the constant is built
        # from `@__DIR__` rather than from an absolute or user-supplied path (T-13-04).
        @test occursin("p11_research_npe.jld2", P13_PRE_CODE)
        @test occursin("@__DIR__", P13_PRE_CODE)
        @test basename(P13_PHASE11_NET) == "p11_research_npe.jld2"
        @test basename(P13_PHASE11_NET) != "trained_npe.jld2"
        # DEVIATION (a) REGRESSION GUARD, at source level: the concrete theta-transform type must
        # be asserted, because a load without p11_architecture.jl in scope only WARNS and hands
        # back a reconstructed placeholder with the right field names.
        @test occursin("P11BoundedThetaTransform", P13_PRE_CODE)
        # And the include preamble must order Phase 11 BEFORE the loader, for the same reason.
        @test findfirst("p11_architecture.jl", P13_PRE_CODE)[1] <
              findfirst("train_npe.jl", P13_PRE_CODE)[1]
    end

    @testset "the conditioning length is derived, not literal" begin
        # DERIVED from Phase 11's own encoder, never typed. `length` of a scalar is 1 in Julia,
        # so the same expression covers a scalar and a vector encoding without a branch.
        @test occursin("length(encode_lambda(", P13_PRE_CODE)
        @test occursin("p13_conditioning_length", P13_PRE_CODE)
        @test occursin("three_way_input_dim", P13_PRE_CODE)
        # The realized width must appear NOWHERE as a literal: at G = 8 with a scalar encoding
        # it is 5*8^2 + 1, and hard-coding it would be silently wrong the moment Phase 11
        # delivered a conditioning vector of another length, with nothing thrown.
        @test !occursin(string(three_way_input_dim(P13_PRE_G, 1)), P13_PRE_CODE)

        n_cond = p13_conditioning_length()
        @test n_cond isa Int
        @test n_cond == length(encode_lambda(P13_PHASE11_REFERENCE_LAMBDA))
        @test n_cond >= 1
        # Memoized: the second call returns the same number without re-deriving it.
        @test p13_conditioning_length() == n_cond
        # The pass-through returns exactly that many components at any lambda in range.
        @test length(p13_encode_lambda(P13_PRE_LAMBDA)) == n_cond
        @test length(p13_encode_lambda(P13_PHASE11_LAMBDA_RANGE[1])) == n_cond
        @test length(p13_encode_lambda(P13_PHASE11_LAMBDA_RANGE[2])) == n_cond
        # The encoding is INHERITED, not re-derived: the endpoints map to 0 and 1 exactly as
        # Phase 11's affine map does, and out-of-range values are NOT clamped (an extrapolation
        # the caller must own -- silent clamping would hide it).
        @test p13_encode_lambda(P13_PHASE11_LAMBDA_RANGE[1])[1] == 0.0
        @test p13_encode_lambda(P13_PHASE11_LAMBDA_RANGE[2])[1] == 1.0
        @test p13_encode_lambda(P13_PHASE11_LAMBDA_RANGE[2] + 1)[1] > 1.0
    end

    @testset "the :widest_rung reference lambda resolves to Phase 11's LAMBDA_MAX" begin
        # tau is really tau-of-lambda (13-10), so the reference must travel with the number. If
        # Phase 11's ladder ever moved, `:widest_rung` would silently re-point tau at a
        # different measurement -- which is why the frozen constant is asserted equal to the
        # live Phase-11 range rather than trusted.
        @test P13_TAU_REFERENCE_LAMBDA_RULE === :widest_rung
        @test P13_PHASE11_LAMBDA_RANGE isa Tuple{<:Real,<:Real}
        @test P13_PHASE11_LAMBDA_RANGE == (LAMBDA_MIN, LAMBDA_MAX)
        @test P13_PHASE11_LAMBDA_RANGE[2] == P13_TAU_REFERENCE_LAMBDA
        @test P13_PHASE11_REFERENCE_LAMBDA == P13_TAU_REFERENCE_LAMBDA
        # LAMBDA_MIN is deliberately never 0: the 0 -> 0.25 step is interpolation onset, a
        # qualitatively different regime, not misalignment sensitivity.
        @test P13_PHASE11_LAMBDA_RANGE[1] > 0
    end

    if isfile(P13_PHASE11_NET)
        @testset "the real artifact satisfies the six-item contract" begin
            h = p13_require_phase11()

            # item 1 -- the trained net, and the RESEARCH lane rather than a shipped bundle.
            @test h.estimator !== nothing
            @test h.meta.research_lane === true
            @test h.meta.shipped === false

            # item 2 -- the FROZEN summary transform, present and a real ZScoreTransform.
            @test h.zt isa StatsBase.ZScoreTransform

            # item 3 -- the FROZEN theta transform. THIS IS THE DEVIATION (a) REGRESSION GUARD:
            # loaded without p11_architecture.jl in scope, JLD2 emits a WARNING and substitutes
            # a `JLD2.ReconstructedMutable` whose field names are (:lo, :hi, :zt) -- so every
            # field-name check still passes while the frozen theta basis is silently not the
            # real one. Only the concrete type distinguishes them.
            @test h.θzt isa P11BoundedThetaTransform

            # item 6 -- theta arity and POSITIONAL field order, against the LIVE prior box so
            # the check moves with prior.jl and a reordered row cannot pass by coincidence.
            bounds = p11_theta_prior_bounds()
            @test length(h.θzt.lo) == 8
            @test length(h.θzt.lo) == length(bounds)
            @test h.θzt.lo == Float64[b[1] for b in bounds]
            @test h.θzt.hi == Float64[b[2] for b in bounds]
            # Rows 5-6 are shift_dx/shift_dy at Phase 11's WIDENED half-width; row 8 is the
            # chromatic epsilon appended LAST so every pre-existing theta index stays valid.
            @test h.θzt.lo[5] == h.θzt.lo[6] == -LAMBDA_MAX
            @test h.θzt.hi[5] == h.θzt.hi[6] == LAMBDA_MAX
            @test h.θzt.hi[8] == -h.θzt.lo[8] > 0

            # item 4 -- the lambda range, bound from Phase 11's Tier-1 pre-registration rather
            # than from net metadata (deviation (b): the handle carries no lambda_range field,
            # and a frozen pre-registered consts file cannot drift with a retrain).
            @test P13_PHASE11_LAMBDA_RANGE isa Tuple{<:Real,<:Real}
            @test P13_PHASE11_LAMBDA_RANGE[2] == P13_TAU_REFERENCE_LAMBDA

            # item 5 -- the F5 imsize provenance: train-joint MUST equal eval-joint, or every
            # rank/coverage/calibration claim is made under the wrong joint.
            @test P11_IMSIZE_SET == P13_IMSIZE_SET
            @test P11_IMSIZE_WEIGHTS == P13_IMSIZE_WEIGHTS

            # The row arithmetic closes WITHOUT a literal: the summary is 2 * (continuous rows)
            # because standardize_summary bypasses the mask block, and the net's input is that
            # plus the conditioning rows.
            @test h.variant === :min
            @test h.d_in == 2 * length(h.zt.mean) + p13_conditioning_length()
            @test isqrt(length(h.zt.mean))^2 == length(h.zt.mean)
            @test isqrt(length(h.zt.mean)) == P13_PRE_G
        end

        @testset "the basis is loaded ONCE and the standardizer is inherited" begin
            b = load_p13_basis()
            # Memoized: every later Phase-13 file shares ONE handle and therefore ONE zt object.
            @test load_p13_basis() === b
            @test load_p13_basis().zt === b.zt

            # The HARD half of assert_frozen_zt: a look-alike handle carrying a transform that
            # did not come from the artifact load fails outright.
            n_cont = length(b.zt.mean)
            fake   = (estimator = b.estimator, θzt = b.θzt,
                      zt = StatsBase.ZScoreTransform(n_cont, 2, copy(b.zt.mean),
                                                     copy(b.zt.scale)),
                      variant = b.variant, d_in = b.d_in, meta = b.meta)
            @test_throws ErrorException assert_frozen_zt(fake, zeros(Float32, 2 * n_cont, 400))
            @test_throws "INHERITED" assert_frozen_zt(fake, zeros(Float32, 2 * n_cont, 400))

            # The SOFT half (Pitfall 5). A pool that is NOT centred and unit-scaled reads as
            # inherited, which is the expected shape of a frozen transform applied to a
            # different pool.
            shifted = Float32.(randn(2 * n_cont, 400) .* 2.5 .+ 3)
            @test assert_frozen_zt(b, shifted).verdict === :inherited
            # A pool that IS exactly centred and unit-scaled is the RE-FIT WARNING SIGN: it is
            # REPORTED (verdict plus a warning), never silently accepted and never an error,
            # because two similar joints can legitimately give similar moments.
            centred = Float32.(randn(2 * n_cont, 4000))
            r = (@test_logs (:warn,) match_mode = :any assert_frozen_zt(b, centred))
            @test r.verdict === :refit_suspected
            # Too few columns to judge is `:indeterminate`, never a false clean bill.
            @test assert_frozen_zt(b, zeros(Float32, 2 * n_cont, 3)).verdict === :indeterminate
        end
    else
        @testset "the real artifact satisfies the six-item contract" begin
            # SKIPPED BECAUSE PHASE 11 HAS NOT DELIVERED ITS RESEARCH NPE. Do NOT fabricate a
            # stand-in artifact and do NOT point the constant at the shipped basis to make this
            # green: an artifact invented to satisfy a gate is the same wrong-basis failure the
            # first testset exists to block. The block assertion above still runs and still
            # passes, so the suite reports the truth in this state rather than going silent.
            @test_skip isfile(P13_PHASE11_NET)
        end
        @testset "the basis is loaded ONCE and the standardizer is inherited" begin
            @test_skip isfile(P13_PHASE11_NET)
        end
    end

    @testset "lambda placement" begin
        # The CORRECT form: the conditioning goes LAST, after the pair encoding, and the
        # realized width equals the DERIVED width. This runs without the artifact, because the
        # encoding surface and the encoder are both pure functions.
        n_cond = p13_conditioning_length()
        Zs = zeros(Float64, 2 * P13_PRE_G^2)
        Zc = zeros(Float64, 2 * P13_PRE_G^2)
        Zp = p13_encode_pair(Zs, Zc, P13_PRE_LAMBDA)
        @test length(Zp) == three_way_input_dim(P13_PRE_G, n_cond)
        @test eltype(Zp) === Float32
        # The conditioning value sits in the LAST n_cond rows, unstandardized and unsmeared.
        @test Zp[end - n_cond + 1:end] == Float32.(p13_encode_lambda(P13_PRE_LAMBDA))
        # ONE lambda governs the pair (pre-registered modelling choice), so the width grows by
        # n_cond, not 2 * n_cond.
        @test length(Zp) - three_way_input_dim(P13_PRE_G, 0) == n_cond

        # THE WRONG FORM -- lambda folded INTO each summary before the pair encoding -- makes
        # each summary odd-length and throws on the kept `iseven` guard. THAT LOUD FAILURE IS A
        # GIFT: without it the conditioning value would be pushed through the contrast block and
        # the net would train happily on a mis-shaped input.
        @test_throws ArgumentError p13_encode_pair(zeros(Float64, 2 * P13_PRE_G^2 + n_cond),
                                                  zeros(Float64, 2 * P13_PRE_G^2 + n_cond),
                                                  P13_PRE_LAMBDA)
        @test_throws "even" p13_encode_pair(zeros(Float64, 2 * P13_PRE_G^2 + n_cond),
                                            zeros(Float64, 2 * P13_PRE_G^2 + n_cond),
                                            P13_PRE_LAMBDA)
        @test_throws "2*G^2" p13_encode_pair(zeros(Float64, 2 * P13_PRE_G^2 + n_cond),
                                             zeros(Float64, 2 * P13_PRE_G^2 + n_cond),
                                             P13_PRE_LAMBDA)
        # Mismatched sample/control lengths are a different, equally loud failure.
        @test_throws DimensionMismatch p13_encode_pair(Zs, zeros(Float64, 2 * P13_PRE_G^2 - 2),
                                                      P13_PRE_LAMBDA)
    end

    @testset "D-04 CPU-only (CUDA never enters this path)" begin
        # Binding the Phase-11 basis is a JLD2 read plus arithmetic; nothing here touches a
        # device. CUDA must not be loaded as a module by this include chain.
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
        @test P13_USE_GPU == false
    end

end
