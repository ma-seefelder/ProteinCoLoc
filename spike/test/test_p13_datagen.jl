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

# spike/test/test_p13_datagen.jl --- the Phase-13 labelled-pool gate (D-02, D-03, D-07, D-11).
#
# THE FOUR SILENT FAILURES THIS FILE EXISTS TO CATCH, one testset each:
#
#   D-02 / Pitfall 5  the standardizer RE-FIT on the Phase-13 pool instead of INHERITED from
#                     Phase 11. The net would still build, train and converge, and every claim
#                     about it would silently refer to a basis it was never trained on.
#   D-03 / Pitfall 4  the conditioning value appended in the wrong place, or twice. The
#                     wrong-place form fails LOUDLY on the `iseven` guard (that is the good
#                     failure); the double-append is SILENT and only a width check catches it.
#   D-07 / F2         stratification that reshaped theta WITHIN a class rather than only
#                     changing class FREQUENCIES. The scalar per-head correction is then
#                     provably unable to recover the pi-scale Bayes factor -- plan 13-07's
#                     NEGATIVE CONTROL 2 demonstrated exactly that on a closed-form toy, and
#                     testset 2 below is the PRODUCTION-SIDE counterpart of that control.
#   Pitfall 6         a global-RNG draw anywhere in generation, which makes the pool
#                     thread-count dependent and the run unreproducible.
#
# EVERYTHING RUNS AT FIXTURE SCALE ON THE FIXTURE STREAM. Every draw here rides
# `p13_fix_datagen_rng` / `p13_fix_rng(P13_FIXTURE_COUNTER)`, never `p13_datagen_rng` and never
# `p13_rng`, so the suite provably cannot pre-observe any part of the stream a reported pool is
# drawn from. NOTHING IN THIS FILE IS DRAWN FROM AN UNSEEDED `rand`/`randn`: an unseeded fixture
# whose assertion depends on the draw's measured moments produced a real intermittent suite-
# reddening flake in plan 13-09, and reseeding from the fixture stream is what repaired it.
#
# Testsets that need the Phase-11 artifact are GUARDED on `isfile(P13_PHASE11_NET)` and
# `@test_skip` explicitly otherwise, naming Phase 11. NO STAND-IN BASIS IS EVER FABRICATED: a
# hand-rolled `zt` would make every assertion here pass while testing nothing, which is the exact
# failure mode D-02 exists to prevent.

using Test
using Statistics
using StatsBase                # the re-fit NEGATIVE CONTROL in testset 3 (test-only, see there)
using HypothesisTests          # the two-sample within-class distribution test (testset 2)
using Random

isdefined(@__MODULE__, :stratify_by_class) ||
    include(joinpath(@__DIR__, "..", "p13", "datagen.jl"))

# --- Fixture-scale knobs (test-local; no reported number depends on any of them) -------------
const TPD_TAU          = p13_tau()      # the frozen Tier-2 tau; never inlined
const TPD_N_SIM        = 150            # simulated fixture pool: 50 items per class
const TPD_N_SEL        = 600            # theta-only stratified draw: 200 per class
const TPD_N_REF        = 1_200          # theta-only UNSTRATIFIED reference draw
const TPD_REF_FIRST    = 1_000_000      # reference index block, disjoint from the stratified one
const TPD_N_TWIN       = 6              # twin-run / thread-independence pool (3 passes over it,
                                        # and the gate runs single-threaded, so keep it small)
const TPD_KS_ALPHA     = 0.001          # deliberately GENEROUS: see testset 2
const TPD_HAVE_P11     = isfile(P13_PHASE11_NET)

# Executable source of the two units under test, comment lines stripped, for the source greps.
_tpd_exec(path) = join(filter(l -> !occursin(r"^\s*#", l), readlines(path)), "\n")
const TPD_DATAGEN_SRC = _tpd_exec(joinpath(@__DIR__, "..", "p13", "datagen.jl"))
const TPD_TRAIN_SRC   = isfile(joinpath(@__DIR__, "..", "p13", "train_three_way.jl")) ?
                        _tpd_exec(joinpath(@__DIR__, "..", "p13", "train_three_way.jl")) : ""

# ONE simulated fixture pool, built once and shared by testsets 1, 3, 4 and 5. Simulating a
# separate pool per testset would quadruple the gate's wall clock for no extra coverage.
# Deliberately drawn under the F5 image-size MIXTURE (`p13_sample_imsize`, the default) rather
# than at one cheap size: a single-size sub-pool is a covariate-shifted slice of the training
# joint, and testset 3's whole subject is what the frozen transform does on the REAL joint.
const TPD_BASIS = TPD_HAVE_P11 ? load_p13_basis() : nothing
const TPD_POOL  = TPD_HAVE_P11 ?
    stratify_by_class(TPD_N_SIM, TPD_BASIS; tau = TPD_TAU, rng_for = p13_fix_datagen_rng) :
    nothing

@testset "P13 datagen and training pool (D-02, D-03, D-07, D-11)" verbose = true begin

    # =========================================================================================
    @testset "stratification hits the pre-registered frequencies" begin
    # =========================================================================================
        # The theta-only arm: `p13_select_indices` is the accept/reject rule itself, and the
        # label is a function of theta alone, so this exercises the whole decision without
        # simulating an image.
        sel = p13_select_indices(300; tau = TPD_TAU, rng_for = p13_fix_datagen_rng)
        m   = class_masses(sel.classes)
        for k in (:exclusion, :random, :coloc)
            @test isapprox(Float64(getproperty(m, k)),
                           Float64(getproperty(P13_TARGET_CLASS_FREQ, k)); atol = 0.05)
        end
        # Whole-class accept/reject to an integer quota makes the realized counts EXACT, not
        # merely close -- the 0.05 tolerance is headroom the design does not need.
        @test sel.quota == (exclusion = 100, random = 100, coloc = 100)
        @test count(==(EXCLUSION), sel.classes) == 100
        @test count(==(RANDOM),    sel.classes) == 100
        @test count(==(COLOC),     sel.classes) == 100
        @test length(sel.accepted) == 300
        @test length(unique(sel.accepted)) == 300          # no index accepted twice
        # Rejection happens, and its overhead is the documented ~1.5x at this tau -- if this
        # were 1.0x, nothing was ever rejected and the stratifier is a no-op.
        @test sel.n_draws > 300
        @test 1.2 <= sel.n_draws / 300 <= 2.5

        # The simulated arm: the SAME frequencies survive the full generate-and-summarize path.
        if TPD_HAVE_P11
            ms = class_masses([it.class for it in TPD_POOL.items])
            for k in (:exclusion, :random, :coloc)
                @test isapprox(Float64(getproperty(ms, k)),
                               Float64(getproperty(P13_TARGET_CLASS_FREQ, k)); atol = 0.05)
            end
            @test length(TPD_POOL.items) == TPD_N_SIM
        else
            @test_skip "simulated stratification arm needs the Phase-11 research NPE"
        end
    end

    # =========================================================================================
    @testset "stratification is by WHOLE CLASS (D-07-i)" begin
    # =========================================================================================
        # THE PRODUCTION-SIDE COUNTERPART OF PLAN 13-07's NEGATIVE CONTROL 2. That control
        # demonstrated on a closed-form toy that a WITHIN-CLASS reshape of theta cannot be
        # repaired by a scalar per-head correction -- the change enters INSIDE the evidence
        # integral and depends on Z (13-RESEARCH F2). So the production requirement is not
        # "stratify carefully", it is `q(theta | class) == pi(theta | class)` EXACTLY, and this
        # testset is the empirical check of that equality on the real generator.
        #
        # Design: compare the within-class rho_sample distribution of the STRATIFIED pool
        # against an UNSTRATIFIED reference drawn from a DISJOINT index block of the same
        # stream. The reference block must be disjoint, or the two samples would be the same
        # items and the test would be vacuous.
        sel = p13_select_indices(TPD_N_SEL; tau = TPD_TAU, rng_for = p13_fix_datagen_rng,
                                 shuffle = false)
        rho_strat = Dict(EXCLUSION => Float64[], RANDOM => Float64[], COLOC => Float64[])
        for (i, idx) in enumerate(sel.accepted)
            d = p13_draw_labelled_item(idx; tau = TPD_TAU, rng_for = p13_fix_datagen_rng)
            @test d.class === sel.classes[i]            # the prescreen is a pure function of idx
            push!(rho_strat[d.class], d.rho_s)
        end

        rho_ref = Dict(EXCLUSION => Float64[], RANDOM => Float64[], COLOC => Float64[])
        for idx in TPD_REF_FIRST:(TPD_REF_FIRST + TPD_N_REF - 1)
            d = p13_draw_labelled_item(idx; tau = TPD_TAU, rng_for = p13_fix_datagen_rng)
            push!(rho_ref[d.class], d.rho_s)
        end

        for c in (EXCLUSION, RANDOM, COLOC)
            @test length(rho_strat[c]) >= 50
            @test length(rho_ref[c])   >= 50
            # GENEROUS THRESHOLD ON PURPOSE. The claim being tested is NOT "the two samples are
            # identical" (they are finite draws) but "no DETECTABLE within-class reshape". A
            # strict alpha would make this testset a coin flip on the fixture stream; a reshape
            # large enough to break the scalar correction would move the whole distribution and
            # be rejected at any sane alpha.
            p = pvalue(ApproximateTwoSampleKSTest(rho_strat[c], rho_ref[c]))
            @test p > TPD_KS_ALPHA
        end

        # The class BOUNDARIES are still respected -- a within-class reshape could preserve the
        # marginal while moving mass across the cut, so the D-05 rule is asserted directly too.
        @test all(r ->  r >  TPD_TAU, rho_strat[COLOC])
        @test all(r ->  r < -TPD_TAU, rho_strat[EXCLUSION])
        # And the D-05 counter-example never appears: no exclusion-class member has rho_s > 0
        # ("less colocalized than the control" is NOT segregation; 13-RESEARCH Pitfall 3).
        @test !any(r -> r > 0, rho_strat[EXCLUSION])
    end

    # =========================================================================================
    @testset "the standardizer is inherited, not re-fit (Pitfall 5)" begin
    # =========================================================================================
        # SOURCE HALF: no z-score transform is CONSTRUCTED anywhere in the Phase-13 data path.
        @test !occursin("ZScoreTransform(", TPD_DATAGEN_SRC)
        @test !occursin("fit(ZScoreTransform", TPD_DATAGEN_SRC)
        @test !occursin("Random.seed!", TPD_DATAGEN_SRC)   # Pitfall 6: layer 1 only, see below

        if TPD_HAVE_P11
            chk = check_frozen_zt(TPD_POOL.items, TPD_BASIS)
            # HARD HALF: `assert_frozen_zt` throws unless `basis.zt` IS the object the memoized
            # single load returned. Getting a verdict at all means that identity held.
            @test chk.verdict in (:inherited, :refit_suspected, :indeterminate)
            @test chk.n_cols == 2 * TPD_N_SIM
            @test chk.n_rows == length(TPD_BASIS.zt.mean)
            @test isfinite(chk.max_abs_mean) && isfinite(chk.max_abs_sd_dev)

            # A HANDLE CARRYING A LOOK-ALIKE TRANSFORM MUST FAIL THE HARD CHECK. This is what
            # actually enforces D-02; everything else here is corroboration.
            faux = (zt = deepcopy(TPD_BASIS.zt),)
            @test_throws ErrorException assert_frozen_zt(faux, TPD_POOL.items[1].Zs)

            # THE DISCRIMINATING ASSERTION, WITH ITS OWN NEGATIVE CONTROL.
            #
            # The soft moment heuristic reports `:refit_suspected` when the standardized pool is
            # centred and unit-scaled to within a COARSE 0.05 tolerance. On this pool that can
            # legitimately happen WITHOUT any re-fit, because D-02 requires the Phase-13 pool to
            # be drawn from the SAME joint the transform was fitted on (the F5 train-joint ==
            # eval-joint invariant) -- so "close to 0 and 1" is the EXPECTED behaviour of a
            # correctly inherited transform, not evidence against it. The heuristic's stated
            # premise ("a frozen transform applied to a DIFFERENT pool should not reproduce
            # that") simply does not hold in this phase.
            #
            # What DOES separate the two hypotheses is the SCALE of the residual. A genuinely
            # re-fit transform reproduces mean 0 and sd 1 to MACHINE PRECISION, because it was
            # solved for on these very columns. An inherited one lands merely nearby. So the
            # negative control below fits a transform on the pool itself -- IN THE TEST, never
            # in `spike/p13/datagen.jl`, which the source half above asserts -- and shows the
            # measured residual is orders of magnitude smaller than the inherited one's.
            cont   = length(TPD_BASIS.zt.mean)
            M      = reduce(hcat, [reduce(hcat, [it.Zs, it.Zc]) for it in TPD_POOL.items])
            C      = Float64.(M[1:cont, :])
            zt_new = StatsBase.fit(ZScoreTransform, C; dims = 2)
            Cnew   = StatsBase.transform(zt_new, C)
            refit_mean = maximum(abs, vec(mean(Cnew, dims = 2)))
            refit_sd   = maximum(abs, vec(std(Cnew, dims = 2)) .- 1)
            @test refit_mean < 1e-8                        # a re-fit IS at machine zero
            @test refit_sd   < 1e-8
            @test chk.max_abs_mean   > 1e-4                # the inherited one is NOT
            @test chk.max_abs_sd_dev > 1e-4
            @test chk.max_abs_mean   > 1e3 * refit_mean    # and by a wide margin
        else
            @test_skip "frozen-zt inheritance check needs the Phase-11 research NPE"
        end
    end

    # =========================================================================================
    @testset "pair width is derived and lambda is appended once (Pitfall 4)" begin
    # =========================================================================================
        if TPD_HAVE_P11
            G        = isqrt(length(TPD_BASIS.zt.mean))
            n_cond   = p13_conditioning_length()
            expected = three_way_input_dim(G, n_cond)

            a = assemble_conditioned_pairs(TPD_POOL.items)
            @test size(a.Z, 1) == expected
            @test size(a.Z, 2) == TPD_N_SIM
            @test all(isfinite, a.Z)

            # The width is DERIVED, and it closes as 5*G^2 + n_cond without either number being
            # written down here either.
            @test expected == p13_ratio_input_dim(G) + n_cond
            @test expected == 2 * (2 * G^2) + G^2 + n_cond

            # THE WRONG FORM FAILS LOUDLY. Folding the conditioning value INTO each summary
            # makes them odd-length; the `iseven` guard is what turns that into an exception
            # instead of a net that trains happily on a smeared conditioning row.
            rng   = p13_fix_rng(P13_FIXTURE_COUNTER)
            odd_s = randn(rng, Float32, 2 * G^2 + 1)
            odd_c = randn(rng, Float32, 2 * G^2 + 1)
            bad   = [(Zs = odd_s, Zc = odd_c, lambda = 1.0, class = RANDOM,
                      rho_s = 0.0, rho_c = 0.0, imsize = (512, 512), idx = 1)]
            err = try
                assemble_conditioned_pairs(bad); nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("even", err.msg) || occursin("2*G^2", err.msg)

            # THE DOUBLE-APPEND IS SILENT, so it is caught by width alone. A hand-built column
            # with the conditioning counted twice must NOT have the derived width.
            gs   = randn(rng, Float32, 2 * G^2)
            gc   = randn(rng, Float32, 2 * G^2)
            cond = p13_encode_lambda(1.0)
            @test length(p13_encode_pair(gs, gc, 1.0)) == expected
            @test length(vcat(p13_pair_encode(gs, gc), Float32.(cond), Float32.(cond))) != expected
            @test length(vcat(p13_pair_encode(gs, gc), Float32.(cond), Float32.(cond))) ==
                  expected + n_cond

            # ONE lambda per pair, and it is the LAST rows of the column -- the conditioning
            # goes after the pair encoding, never inside it.
            col = p13_encode_pair(gs, gc, 1.0)
            @test col[(expected - n_cond + 1):expected] == Float32.(cond)
            @test col[1:(expected - n_cond)] == p13_pair_encode(gs, gc)
        else
            @test_skip "pair width / lambda placement checks need the Phase-11 research NPE"
        end
    end

    # =========================================================================================
    @testset "targets are the four-row matrix (D-11)" begin
    # =========================================================================================
        if TPD_HAVE_P11
            a = assemble_conditioned_pairs(TPD_POOL.items)
            @test size(a.targets, 1) == 4
            @test size(a.targets, 2) == length(TPD_POOL.items)
            @test length(a.classes)  == length(TPD_POOL.items)
            @test length(a.lambdas)  == length(TPD_POOL.items)

            # A RANDOM-class column is [0, 1, 0, 1]: the shared NEGATIVE of BOTH heads, which is
            # what makes the two logits ratios against the SAME reference (D-08).
            jr = findfirst(==(RANDOM), a.classes)
            @test jr !== nothing
            @test a.targets[:, jr] == Float32[0, 1, 0, 1]
            # A COLOC column switches the exclusion head OFF, and vice versa (D-11). Under
            # one-vs-rest those participation weights would both be 1 and each head's logit
            # would silently become a MIXTURE Bayes factor.
            jc = findfirst(==(COLOC), a.classes)
            je = findfirst(==(EXCLUSION), a.classes)
            @test a.targets[:, jc] == Float32[1, 1, 0, 0]
            @test a.targets[:, je] == Float32[0, 0, 1, 1]

            # Every lambda is inside the trained Phase-11 range, and ONE lambda governs a pair.
            @test all(l -> P13_PHASE11_LAMBDA_RANGE[1] <= l <= P13_PHASE11_LAMBDA_RANGE[2],
                      a.lambdas)
        else
            @test_skip "four-row target checks need the Phase-11 research NPE"
        end
    end

    # =========================================================================================
    @testset "reproducibility twin run" begin
    # =========================================================================================
        # THETA-ONLY ARM (runs with or without the artifact): the accept/reject walk AND the
        # deterministic permutation reproduce exactly.
        s1 = p13_select_indices(300; tau = TPD_TAU, rng_for = p13_fix_datagen_rng)
        s2 = p13_select_indices(300; tau = TPD_TAU, rng_for = p13_fix_datagen_rng)
        @test s1.accepted == s2.accepted
        @test s1.classes  == s2.classes
        @test s1.n_draws  == s2.n_draws
        # Index 0 is reserved for the permutation stream, so it is never an accepted item.
        @test !(0 in s1.accepted)

        if TPD_HAVE_P11
            # SIMULATED ARM, copying the shape of the existing thread-independence gate
            # (`spike/test/test_data_pipeline.jl:253-274`): the PARALLEL path must equal the
            # SERIAL fallback BITWISE, because every column is a pure function of its global
            # index -- no shared mutable RNG, no locks, no `Random.seed!` in the loop.
            idxs = s1.accepted[1:TPD_N_TWIN]
            a = p13_simulate_indices(idxs, TPD_BASIS; tau = TPD_TAU,
                                     rng_for = p13_fix_datagen_rng, parallel = false)
            b = p13_simulate_indices(idxs, TPD_BASIS; tau = TPD_TAU,
                                     rng_for = p13_fix_datagen_rng, parallel = true)
            Za = assemble_conditioned_pairs(a).Z
            Zb = assemble_conditioned_pairs(b).Z
            @test Za == Zb                                    # BITWISE, not approximately
            @test [x.class for x in a] == [x.class for x in b]
            @test [x.lambda for x in a] == [x.lambda for x in b]

            # OUT-OF-ORDER generation re-sorts to identical columns -- the same property the
            # Phase-3 gate asserts, restated for the Phase-13 item.
            perm = shuffle(MersenneTwister(99), collect(1:TPD_N_TWIN))
            c    = p13_simulate_indices(idxs[perm], TPD_BASIS; tau = TPD_TAU,
                                        rng_for = p13_fix_datagen_rng, parallel = false)
            invp = sortperm(perm)
            @test assemble_conditioned_pairs([c[k] for k in invp]).Z == Za
        else
            @test_skip "simulated twin-run / thread-independence arm needs the Phase-11 net"
        end
    end

    # =========================================================================================
    @testset "no forbidden basis is reachable" begin
    # =========================================================================================
        # D-02 admits NO fallback. The shipped grid-8 bundle and the pre-Phase-11 net are a
        # DIFFERENT input surface (128 rows, no conditioning row, the narrow shift prior, 7
        # theta fields), and training on either would invalidate D-02, D-03 and D-12 while
        # looking entirely normal.
        for src in (TPD_DATAGEN_SRC, TPD_TRAIN_SRC)
            @test !occursin("trained_npe.jld2", src)
            @test !occursin("grid_8", src)
            @test !occursin("artifacts/", src)
        end
        # The single binding site is `preconditions.jl`, reached by include, not by path.
        @test occursin("preconditions.jl", TPD_DATAGEN_SRC)
        @test occursin("load_p13_basis", TPD_DATAGEN_SRC)
        # And the realized input width is never written down.
        @test !occursin("321", TPD_DATAGEN_SRC)
        @test !occursin("321", TPD_TRAIN_SRC)
    end

    # =========================================================================================
    @testset "the recipe is referenced, never inlined" begin
    # =========================================================================================
        # D-10's attribution argument holds only if NOTHING in the recipe was re-chosen. A
        # literal left behind in the trainer is how a "temporary" tuning survives into the
        # reported run, so the absence of the literals is asserted rather than reviewed.
        @test TPD_TRAIN_SRC != ""                    # the trainer exists and was read
        for lit in ("48_000", "2.5e-4", "1e-4", "300")
            @test !occursin(lit, TPD_TRAIN_SRC)
        end
        # Each value is REFERENCED from the frozen pre-registration by name.
        for name in ("P13_TRAIN_N", "P13_EPOCHS", "P13_BATCHSIZE", "P13_LR",
                     "P13_WEIGHT_DECAY", "P13_VAL_FRAC", "P13_STOPPING_EPOCHS")
            @test occursin(name, TPD_TRAIN_SRC)
        end
        # And the pre-registration still holds the shipped values, so "referenced" is not
        # merely referencing something that drifted.
        @test P13_TRAIN_N         == 48_000
        @test P13_EPOCHS          == 300
        @test P13_BATCHSIZE       == 128
        @test P13_LR              == 2.5e-4
        @test P13_WEIGHT_DECAY    == 1e-4
        @test P13_VAL_FRAC        == 0.15
        @test P13_STOPPING_EPOCHS == 40
        @test P13_USE_GPU         == false

        # Layer 2 of the RNG discipline lives in the TRAINER and nowhere else, and it sits
        # IMMEDIATELY BEFORE the train call (Pitfall 6: Flux's shuffling DataLoader draws from
        # the global RNG, so a counter-based stream alone is not enough).
        @test occursin("Random.seed!", TPD_TRAIN_SRC)
        lines = split(TPD_TRAIN_SRC, "\n")
        # The CALL, not a mention: `using Random` carries an inline `# Random.seed!` note that a
        # bare `occursin` would match first, and an adjacency check anchored on a comment proves
        # nothing about where the seeding actually happens.
        k = findall(l -> occursin(r"^\s*Random\.seed!\(", l), lines)
        @test length(k) == 1
        @test occursin("train_three_way(", lines[k[1] + 1])
    end

    # =========================================================================================
    @testset "CPU-only: CUDA is never loaded" begin
    # =========================================================================================
        # CPU-only is the reproducible baseline (CLAUDE.md, D-01) and `P13_USE_GPU` is
        # pre-registered false. CUDA is never imported by this path, and no package was
        # installed to run it.
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
        @test !occursin("CUDA", TPD_DATAGEN_SRC)
        @test !occursin("CUDA", TPD_TRAIN_SRC)
        @test P13_USE_GPU == false
    end

end
