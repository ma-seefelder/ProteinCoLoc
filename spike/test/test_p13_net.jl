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

# spike/test/test_p13_net.jl --- the Phase-13 three-way evidence-net gate (D-08..D-11).
#
# WHAT THIS FILE IS FOR. The three-way head is the place this phase forks away from the shipped
# binary net's gate lineage, so every claim the fork rests on is asserted here rather than
# assumed:
#
#   D-09 / A5 -- the whole route depends on NeuralEstimators' `train` accepting a CUSTOM
#         estimator subtype AND honouring a PASSED loss. Both are properties of unexported
#         internals of a pre-1.0 package, so they are proven by an actual two-epoch training run
#         on a 200-sample toy BEFORE any expensive run. A failure here is the trigger for the
#         pre-declared route-R2 hand-rolled loop, which would have to be recorded as a DECLARED
#         deviation -- never silently substituted, because R2 weakens D-10's attribution claim.
#   D-10 -- the trunk must be the shipped topology VERBATIM, so the widths and activations are
#         walked layer by layer instead of trusted to a builder that looks right.
#   D-11 -- each head must train ONLY on its own class pair. The distinguishing property is that
#         perturbing an INACTIVE logit leaves the loss BIT-IDENTICAL; a one-vs-rest scheme would
#         fail exactly that assertion and pass every shape check.
#   D-08 -- the read surface must return `random = 0.0` by construction and must subtract the
#         PER-HEAD measured correction, never the binary read surface's own term. The absence of
#         the binary term is asserted by grepping the executable source, because a wrong copy is
#         numerically invisible in the file it would have been copied from.
#   D-03 -- the input width must be DERIVED. The literal is asserted absent from executable code.
#
# Mirrors test_bf.jl / test_p13_labels.jl: license header -> using Test -> guarded include of the
# unit under test -> fixtures computed ONCE at module level -> one outer @testset. No image
# simulation and no figure call; the only training is the A5 smoke.
#
# THE SUITE IS RUNNABLE BEFORE THE TIER-2 TAU COMMIT. The net does not depend on tau (labelling
# does), so nothing here calls the tau accessor and the last testset asserts tau is still absent.

using Test
using Random
using Statistics

# Unit under test: the three-way net surface (pulls p13/consts.jl -> p13/labels.jl).
isdefined(@__MODULE__, :ThreeWayEvidenceNet) ||
    include(joinpath(@__DIR__, "..", "p13", "net.jl"))

# The EXECUTABLE source of the unit under test, with whole-line comments stripped BEFORE any
# count, so a grep assertion cannot be satisfied (or defeated) by prose. Docstrings survive the
# strip on purpose: a forbidden identifier quoted in a docstring is still a forbidden identifier
# a reader may copy, so the docstrings in net.jl deliberately name none of them.
const P13_NET_SRC  = read(joinpath(@__DIR__, "..", "p13", "net.jl"), String)
const P13_NET_CODE = join(filter(l -> !startswith(strip(l), "#"),
                                 split(P13_NET_SRC, '\n')), '\n')

# --- The A5 toy, built ONCE -----------------------------------------------------------------
# A CLASS-SEPARATED GAUSSIAN TOY: exclusion samples sit at -NET_FIX_DELTA, coloc at
# +NET_FIX_DELTA, random at 0, so both heads have a real signal to find and two epochs at the
# COPIED learning rate move the loss measurably. The toy is small and synthetic on purpose --
# its job is to exercise the training PLUMBING (route R1), not to teach anything.
#
# ALL fixtures ride p13_fix_rng(P13_FIXTURE_COUNTER), so this gate NEVER consumes (and so never
# pre-observes) a REPORTED Phase-13 stream (D-01).
const NET_FIX_DIM   = 32
const NET_FIX_DELTA = 2.0f0

# THE SPLIT IS 160/40, NOT the pre-registered val_frac. NeuralEstimators' `_DataLoader` is built
# with `partial = false` (train.jl:668-676), so it DROPS an incomplete final batch: a validation
# set smaller than `batchsize` would yield ZERO batches and a 0/0 = NaN validation risk. 40 >= 32
# keeps at least one full validation batch. This is a property of the toy's size, not a change to
# any pre-registered constant.
const NET_FIX_NTR = 160
const NET_FIX_NVA = P13_SMOKE_N - NET_FIX_NTR

function _net_toy(n::Integer; input_dim::Integer = NET_FIX_DIM,
                  delta::Float32 = NET_FIX_DELTA)
    rng     = p13_fix_rng(P13_FIXTURE_COUNTER)
    classes = Vector{ThreeWayClass}(undef, n)
    Z       = Float32.(randn(rng, input_dim, n))
    for j in 1:n
        k = rand(rng, 1:3)
        classes[j] = k == 1 ? EXCLUSION : (k == 2 ? RANDOM : COLOC)
        classes[j] === COLOC     && (Z[:, j] .+= delta)
        classes[j] === EXCLUSION && (Z[:, j] .-= delta)
    end
    return (Z = Z, targets = target_matrix(classes), classes = classes)
end

const NET_FIX_TOY = _net_toy(P13_SMOKE_N)

# The hand-computed per-cell logit-BCE, written out in Float64 from its definition rather than
# called from Flux, so testset 5 checks the loss against arithmetic and not against itself.
_hand_bce(logit::Real, target::Real) =
    target == 1 ? log1p(exp(-Float64(logit))) : log1p(exp(Float64(logit)))

@testset "P13 three-way evidence net (D-08, D-09, D-10, D-11)" verbose = true begin

    @testset "SC1 shape: one forward pass, 2-by-n logits" begin
        net = build_three_way_net(320)
        @test net isa ThreeWayEvidenceNet
        @test size(net(Float32.(randn(320, 5)))) == (2, 5)
        # The trunk's own output is the learned-summary bottleneck the head reads.
        @test size(net.trunk(Float32.(randn(320, 5)))) == (P13_NUM_SUMMARIES, 5)
        # SC1's bindable content is LITERAL IN THE SOURCE: one head applied to one trunk pass.
        @test occursin("m.head(m.trunk(", P13_NET_SRC)
        # A vector input is a legitimate single column for the read surface.
        @test size(net(reshape(Float32.(randn(320)), :, 1))) == (2, 1)
    end

    @testset "D-10 topology is the shipped one verbatim" begin
        # The pre-registration's copied constants are the shipped net's constants.
        @test P13_SUMMARY_WIDTH == 256
        @test P13_NUM_SUMMARIES == 64

        input_dim = 320
        net = build_three_way_net(input_dim)
        L = net.trunk.layers
        @test length(L) == 4
        # Dense weights are (out, in): input_dim -> 256 -> 256 -> 256 -> 64.
        @test size(L[1].weight) == (P13_SUMMARY_WIDTH, input_dim)
        @test size(L[2].weight) == (P13_SUMMARY_WIDTH, P13_SUMMARY_WIDTH)
        @test size(L[3].weight) == (P13_SUMMARY_WIDTH, P13_SUMMARY_WIDTH)
        @test size(L[4].weight) == (P13_NUM_SUMMARIES, P13_SUMMARY_WIDTH)
        # gelu on the first three, identity on the bottleneck -- exactly the shipped trunk.
        @test L[1].σ === gelu && L[2].σ === gelu && L[3].σ === gelu
        @test L[4].σ === identity
        # The ONLY structural change is the head: two logits, not one model-index unit.
        @test size(net.head.weight) == (2, P13_NUM_SUMMARIES)
        # The architecture is READ BACK off the layers, so an artifact cannot disagree with it.
        @test three_way_arch(net) ==
              (input_dim = input_dim, width = P13_SUMMARY_WIDTH,
               num_summaries = P13_NUM_SUMMARIES)
        # The shipped guard is kept verbatim.
        @test_throws ArgumentError build_three_way_net(0)
    end

    @testset "input width is derived" begin
        # ratio_input_dim(G) = 5*G^2, plus the conditioning length.
        @test three_way_input_dim(8, 1) == 5 * 8^2 + 1
        @test three_way_input_dim(8, 2) == 5 * 8^2 + 2
        @test three_way_input_dim(4, 1) == 5 * 4^2 + 1
        @test three_way_input_dim(8, 0) == p13_ratio_input_dim(8)
        @test_throws ArgumentError three_way_input_dim(8, -1)
        @test_throws ArgumentError p13_ratio_input_dim(0)
        # THE TRIPWIRE: the width must never be a literal in executable code, because Phase 11
        # may deliver a conditioning vector whose length is not 1.
        @test !occursin("321", P13_NET_CODE)
        @test P13_INPUT_WIDTH_RULE === :ratio_input_dim_plus_ncond
    end

    @testset "D-11 masked loss ignores inactive cells (BIT-IDENTICAL)" begin
        # THIS IS THE TESTSET THAT DISTINGUISHES D-11 FROM ONE-VS-REST. Under one-vs-rest the
        # coloc head would train on the exclusion sample and the exclusion head on the coloc
        # sample, so perturbing those cells WOULD move the loss -- and every shape test would
        # still pass. Exact equality is the point: an approximate comparison would also accept a
        # tiny leak from an inactive cell, and a tiny leak is still a broken head boundary.
        y = target_matrix([COLOC, EXCLUSION, RANDOM])
        @test size(y) == (4, 3)
        yhat = Float32[0.7 -0.4 0.1; -0.2 0.9 -0.6]

        base = masked_two_head_bce(yhat, y)
        @test isfinite(base)

        # Column 1 is COLOC: w_E == 0, so row 2 is INACTIVE. Column 2 is EXCLUSION: w_C == 0,
        # so row 1 is INACTIVE. Perturb ONLY those two cells, by a large amount.
        pert = copy(yhat)
        pert[2, 1] += 50.0f0
        pert[1, 2] -= 50.0f0
        # BIT-IDENTICAL: exact equality, deliberately NOT an approximate comparison.
        @test masked_two_head_bce(pert, y) == base
        @test masked_two_head_bce(pert, y) === base

        # And the inverse: an ACTIVE cell must move the loss.
        act = copy(yhat)
        act[1, 1] += 1.0f0                      # the coloc head on a COLOC sample: active
        @test masked_two_head_bce(act, y) != base
        act2 = copy(yhat)
        act2[2, 3] += 1.0f0                     # the exclusion head on a RANDOM sample: active
        @test masked_two_head_bce(act2, y) != base
    end

    @testset "masked loss denominator counts active cells" begin
        # Two samples: one COLOC (coloc head active, exclusion head OFF) and one RANDOM (BOTH
        # heads active, since random is the shared negative). Active cells = 3, so the
        # denominator is 3 -- not n = 2 and not 2n = 4.
        y = target_matrix([COLOC, RANDOM])
        @test y[:, 1] == Float32[1, 1, 0, 0]
        @test y[:, 2] == Float32[0, 1, 0, 1]
        yhat = Float32[1.0 -0.5; 3.0 0.25]       # yhat[2,1] = 3.0 is the INACTIVE cell
        expected = (_hand_bce(1.0, 1) +          # coloc head, COLOC sample, target 1
                    _hand_bce(-0.5, 0) +         # coloc head, RANDOM sample, target 0
                    _hand_bce(0.25, 0)) / 3      # exclusion head, RANDOM sample, target 0
        @test isapprox(Float64(masked_two_head_bce(yhat, y)), expected; atol = 1e-6)
    end

    @testset "A5 smoke: NeuralEstimators train accepts the subtype AND the custom loss" begin
        # THE A5 MITIGATION. Route R1 rests on three unexported internals of a pre-1.0 package
        # (`_loss`, `_inputoutput`, `_construct_train_state`). This testset runs a real two-epoch
        # training pass on a 200-sample toy so that dependency is proven cheaply BEFORE any
        # expensive run. IF IT FAILS, the pre-declared route-R2 hand-rolled Flux loop is the
        # fallback -- and using it must be RECORDED AS A DECLARED DEVIATION, because R2 must
        # re-implement CosAnneal, patience and best-checkpointing identically or D-10's
        # "difference attributable to the head alone" claim is no longer true.
        @test P13_SMOKE_N == 200 && P13_SMOKE_EPOCHS == 2

        Z, tg = NET_FIX_TOY.Z, NET_FIX_TOY.targets
        @test size(Z) == (NET_FIX_DIM, P13_SMOKE_N)
        @test size(tg) == (4, P13_SMOKE_N)
        # The toy is non-degenerate: every class is present, so both heads see both of theirs.
        @test all(c -> count(==(c), NET_FIX_TOY.classes) > 0, (EXCLUSION, RANDOM, COLOC))

        Z_tr, Z_va   = Z[:, 1:NET_FIX_NTR], Z[:, NET_FIX_NTR+1:end]
        tg_tr, tg_va = tg[:, 1:NET_FIX_NTR], tg[:, NET_FIX_NTR+1:end]
        @test size(Z_va, 2) == NET_FIX_NVA >= 32     # at least one full validation batch

        fresh  = build_three_way_net(NET_FIX_DIM)
        before = masked_two_head_bce(fresh(Z), tg)
        trained = train_three_way(fresh, tg_tr, tg_va, Z_tr, Z_va;
                                  epochs = P13_SMOKE_EPOCHS, batchsize = 32,
                                  savepath = nothing, verbose = false)
        after = masked_two_head_bce(trained(Z), tg)

        @test trained isa ThreeWayEvidenceNet
        @test size(trained(Z)) == (2, P13_SMOKE_N)
        @test isfinite(before) && isfinite(after)
        # STRICTLY below: `train` returns the BEST-VALIDATION state, which is initialized to the
        # UNTRAINED state, so an equal loss would mean no epoch ever improved -- i.e. the loss
        # never reached the optimizer.
        @test after < before
    end

    @testset "custom loss is actually honoured" begin
        # THE LANDMINE THIS GUARDS AGAINST: the shipped binary estimator type hard-codes its
        # loss, so a passed loss is SILENTLY discarded -- no warning, no error. If that happened
        # here, two runs from the SAME seed under DIFFERENT losses would be bit-identical.
        # A variant that switches the exclusion head OFF entirely: its gradients must differ.
        coloc_only_bce(yhat, y) = begin
            lC = logitbinarycrossentropy.(view(yhat, 1, :), view(y, 1, :); agg = identity)
            wC = view(y, 2, :)
            sum(wC .* lC) / (sum(wC) + eps(Float32))
        end

        Z, tg        = NET_FIX_TOY.Z, NET_FIX_TOY.targets
        Z_tr, Z_va   = Z[:, 1:NET_FIX_NTR], Z[:, NET_FIX_NTR+1:end]
        tg_tr, tg_va = tg[:, 1:NET_FIX_NTR], tg[:, NET_FIX_NTR+1:end]

        # ONE initialization, DEEP-COPIED, so the two runs differ in the loss and in NOTHING
        # else. Two freshly-built nets would differ by their random initialization alone and the
        # comparison below would be vacuous.
        net0 = build_three_way_net(NET_FIX_DIM)
        a = train_three_way(deepcopy(net0), tg_tr, tg_va, Z_tr, Z_va;
                            epochs = P13_SMOKE_EPOCHS, batchsize = 32, savepath = nothing,
                            verbose = false, loss = masked_two_head_bce)
        b = train_three_way(deepcopy(net0), tg_tr, tg_va, Z_tr, Z_va;
                            epochs = P13_SMOKE_EPOCHS, batchsize = 32, savepath = nothing,
                            verbose = false, loss = coloc_only_bce)
        # Same seed, same data, same initialization -- only the loss differs. If the loss were
        # silently discarded, these two nets would be bit-identical.
        @test a.head.weight != b.head.weight
    end

    @testset "lambda placement (Pitfall 4)" begin
        zs = zeros(Float32, 128)
        zc = zeros(Float32, 128)
        Zp = encode_conditioned_pair(zs, zc, [0.5f0])
        @test length(Zp) == three_way_input_dim(8, 1)
        @test length(Zp) == length(p13_pair_encode(zs, zc)) + 1
        @test Zp[end] == 0.5f0                      # the conditioning row is LAST
        @test P13_LAMBDA_PLACEMENT === :append_after_pair_encode

        # WRONG FORM (i): folding the conditioning row into the summaries makes them odd-length,
        # and the kept `iseven` guard throws. This is the GOOD failure.
        err = try
            p13_pair_encode(zeros(Float32, 129), zeros(Float32, 129))
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("even", sprint(showerror, err))
        @test_throws ArgumentError encode_conditioned_pair(zeros(Float32, 129),
                                                          zeros(Float32, 129), [0.5f0])
        # Mismatched summary lengths are a DimensionMismatch, not a silent truncation.
        @test_throws DimensionMismatch p13_pair_encode(zeros(Float32, 128), zeros(Float32, 126))

        # WRONG FORM (ii): the conditioning row counted twice. SILENT unless the length is
        # asserted, which is why it is asserted.
        twice = vcat(p13_pair_encode(zs, zc), Float32[0.5, 0.5])
        @test length(twice) == three_way_input_dim(8, 2)
        @test length(twice) != three_way_input_dim(8, 1)

        # A conditioning vector longer than 1 is supported by construction (Phase 11 may deliver
        # one), and the realized width tracks it.
        @test length(encode_conditioned_pair(zs, zc, [0.5f0, 0.25f0])) ==
              three_way_input_dim(8, 2)
    end

    @testset "D-08 read surface" begin
        net = build_three_way_net(320)
        Zv  = Float32.(randn(320))
        Zm  = reshape(Zv, :, 1)
        s   = net(Zm)
        hlo = (coloc = 0.3, exclusion = -0.2)

        r = three_way_log_bf(net, Zv, hlo)
        # THE REFERENCE IS IDENTICALLY ZERO BY CONSTRUCTION, not by measurement (D-08).
        @test r.random === 0.0
        # The PER-HEAD correction is subtracted -- and only it.
        @test isapprox(r.coloc,     Float64(s[1, 1]) - 0.3; atol = 1e-12)
        @test isapprox(r.exclusion, Float64(s[2, 1]) + 0.2; atol = 1e-12)
        # A vector and the equivalent single-column matrix are the same read.
        rm = three_way_log_bf(net, Zm, hlo)
        @test rm.coloc == r.coloc && rm.exclusion == r.exclusion && rm.random === 0.0
        # A zero correction leaves the raw logits, which is the sanity anchor for the above.
        r0 = three_way_log_bf(net, Zv, (coloc = 0.0, exclusion = 0.0))
        @test isapprox(r0.coloc, Float64(s[1, 1]); atol = 1e-12)

        # The correction consumed here is EXACTLY the shape labels.jl produces.
        hlo_measured = head_log_odds([COLOC, COLOC, RANDOM, RANDOM, RANDOM, EXCLUSION])
        @test keys(hlo_measured) == (:coloc, :exclusion)
        @test isfinite(three_way_log_bf(net, Zv, hlo_measured).coloc)

        # D-13 scores the UNCORRECTED probabilities: calibration is a property of the classifier
        # under its own training frequencies.
        p = three_way_probs(net, Zv)
        @test 0.0 < p.coloc < 1.0 && 0.0 < p.exclusion < 1.0
        @test isapprox(p.coloc, 1 / (1 + exp(-Float64(s[1, 1]))); atol = 1e-6)

        # CPU-only, on the read path too.
        @test_throws ArgumentError three_way_log_bf(net, Zv, hlo; use_gpu = true)
    end

    @testset "the binary correction is not copied (Pitfall 1)" begin
        # The binary read surface's term is ~0.01 nat because its labels are balanced, so a
        # wrong copy is INVISIBLE in the file it would be copied from. Here the per-head terms
        # are 15-50x larger. The absence is therefore asserted on the executable source.
        @test !occursin("log_prior_odds", P13_NET_CODE)
        @test occursin("head_log_odds", P13_NET_CODE)
        # And the two rejected architectures are absent from executable code as well.
        @test !occursin("RatioEstimator(", P13_NET_CODE)
        @test !occursin("softmax", P13_NET_CODE)
        # Named limit 3: the retired KDE baseline is not even reachable from the spike env.
        @test !occursin("kde(", P13_NET_CODE)
        @test !occursin("quadgk(", P13_NET_CODE)
    end

    @testset "persistence round-trips atomically" begin
        net = build_three_way_net(64; num_summaries = 8, width = 16)
        hlo = (coloc = -0.5333, exclusion = -0.1654)
        dir = mktempdir()
        path = joinpath(dir, "three_way_net.jld2")

        @test save_three_way(path, net, hlo; meta = (note = "gate fixture",)) == path
        @test isfile(path)
        # The `.tmp` is the atomicity mechanism, so it must not survive a successful commit.
        @test !isfile(path * ".tmp")

        h = load_three_way(path)
        @test h.net isa ThreeWayEvidenceNet
        @test (h.input_dim, h.num_summaries, h.width) == (64, 8, 16)
        @test h.head_log_odds == hlo
        @test h.schema_version == P13_NET_SCHEMA
        @test h.cut_variant === P13_CUT_VARIANT
        # tau is still Tier-2, so the artifact records its absence rather than inventing a value.
        @test h.tau === nothing
        # The provenance sha ties the artifact to the frozen pre-registration.
        @test h.consts_sha == p13_consts_sha() && length(h.consts_sha) == 64
        @test h.meta.note == "gate fixture"

        # REBUILD-FROM-STATE reproduces the saved net's output (the narrow-surface contract).
        Zf = Float32.(randn(64, 4))
        @test maximum(abs.(net(Zf) .- h.net(Zf))) <= 1e-6

        @test_throws ErrorException load_three_way(joinpath(dir, "absent.jld2"))
    end

    @testset "shape confusion throws at the input boundary" begin
        net = build_three_way_net(NET_FIX_DIM)
        Z   = Float32.(randn(NET_FIX_DIM, 8))
        tg  = target_matrix(fill(RANDOM, 8))
        # A 3-row target matrix is the shape error that would otherwise train a head against the
        # wrong row of the legend.
        @test_throws ArgumentError train_three_way(net, tg[1:3, :], tg, Z, Z;
                                                   epochs = 1, batchsize = 4, savepath = nothing,
                                                   verbose = false)
        @test_throws ArgumentError train_three_way(net, tg, tg[1:3, :], Z, Z;
                                                   epochs = 1, batchsize = 4, savepath = nothing,
                                                   verbose = false)
        # A column-count disagreement between summaries and targets.
        @test_throws ArgumentError train_three_way(net, tg, tg, Z[:, 1:4], Z;
                                                   epochs = 1, batchsize = 4, savepath = nothing,
                                                   verbose = false)
    end

    @testset "CPU-only gate (D-10)" begin
        net = build_three_way_net(NET_FIX_DIM)
        Z   = Float32.(randn(NET_FIX_DIM, 8))
        tg  = target_matrix(fill(RANDOM, 8))
        @test P13_USE_GPU == false
        @test_throws ArgumentError train_three_way(net, tg, tg, Z, Z; use_gpu = true,
                                                   epochs = 1, batchsize = 4,
                                                   savepath = nothing, verbose = false)
        # CUDA is never imported by this path.
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

    @testset "the net does not depend on tau (still Tier-2)" begin
        # D-06 keeps tau out of the pre-registration until it is MEASURED. The net surface must
        # therefore be fully testable before that commit -- and it is, because tau enters through
        # the LABEL boundary, not through the architecture.
        @test !isdefined(@__MODULE__, :P13_TAU)
    end

end
