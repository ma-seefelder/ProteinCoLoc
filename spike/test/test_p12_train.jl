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

# spike/test/test_p12_train.jl --- the 12-14 trainer's two silent-failure modes, covered by tests
# that fail when the failure is introduced.
#
# THE TWO FAILURE MODES, both of which destroy an amortized calibration claim without raising
# anything at run time:
#   (1) A LEAKED STANDARDIZER. A transform fitted over the validation columns makes every later
#       coverage number optimistic by an amount nothing downstream can detect.
#   (2) AN OUT-OF-DISTRIBUTION HELD-OUT ENCODING. If masking happened AFTER standardization, a
#       held-out region would arrive at the net carrying the ROW MEAN -- "this region is perfectly
#       average" -- where read time says "this region is absent". Testset 3 is the single
#       assertion separating those two encodings.
#
# Everything runs on a 64-sample pool at (128, 128) inside `mktempdir()`, so the suite stays well
# inside its per-commit budget and NOTHING is written under `spike/data/cache/` or `spike/npe/`.

using Test
using Statistics
using StatsBase
using Random

isdefined(@__MODULE__, :train_p12_npe) ||
    include(joinpath(@__DIR__, "..", "npe", "train_p12_npe.jl"))

# 12-15's `select_prior` and its scoring runner. Loaded WITHOUT running the job: the runner's last
# line is `isdefined(@__MODULE__, :P12_MINISPIKE_LOAD_ONLY) || main()`.
#
# WHY THESE TESTS LIVE HERE AND NOT IN `test_p12_consts.jl`, which is where 12-15's `files_modified`
# would have put them: the runner transitively loads `train_p12_npe.jl` and NeuralEstimators, and
# `test_p12_consts.jl` is the Tier-1 PRE-REGISTRATION GATE and the FIRST file the aggregator runs.
# Making the pre-registration gate depend on the neural stack would mean a NeuralEstimators problem
# breaks the check that the frozen constants are intact -- a coupling regression. This file already
# loads that stack, so the tests are free here.
P12_MINISPIKE_LOAD_ONLY = true
isdefined(@__MODULE__, :select_prior) ||
    include(joinpath(@__DIR__, "..", "validation", "run_p12_minispike.jl"))

"""
Strip whole-line comments before making a source-text assertion, so a claim recorded IN a comment
cannot satisfy a check about the CODE. Mirrors `test_stage6_regression.jl:88-100`.
"""
_p12train_strip_comments(src) =
    join(filter(l -> !startswith(strip(l), "#"), split(src, "\n")), "\n")

const _P12TRAIN_SRC = read(joinpath(@__DIR__, "..", "npe", "train_p12_npe.jl"), String)
const _P12TRAIN_CODE = _p12train_strip_comments(_P12TRAIN_SRC)

@testset "P12 mini-spike trainer and its spend guard (12-14)" verbose = true begin

    # A shared 64-sample pool at (128,128) in a temp cache root. Generated ONCE.
    tmp_cache = mktempdir()
    pool64 = generate_p12_pool(64; shard_size = 32, imsize_sampler = _ -> (128, 128),
                               cache_root = tmp_cache, imsize_tag = :pinned_128x128,
                               verbose = false)

    @testset "both transforms are fit on the TRAIN split only (leak-freedom)" begin
        # THIS TESTS THE PROPERTY, NOT THE LINE OF CODE. Fit on a known train block; then
        # perturb ONLY the validation columns to extremes and refit on the same train block. If
        # any validation column could reach the fit, the two transforms would differ.
        pool = load_p12_pool(pool64)
        Z = Float64.(pool.summary_min)
        TH = Float64.(pool.theta)
        ntr = 48
        tr = 1:ntr
        va = (ntr + 1):size(Z, 2)

        zt_a, tzt_a = fit_p12_transforms(view(Z, :, tr), view(TH, :, tr))

        Z2 = copy(Z); TH2 = copy(TH)
        Z2[:, va]  .= 1.0e6          # extremes that would visibly move any mean or scale
        TH2[:, va] .= -1.0e6
        zt_b, tzt_b = fit_p12_transforms(view(Z2, :, tr), view(TH2, :, tr))

        @test zt_a.mean  == zt_b.mean
        @test zt_a.scale == zt_b.scale
        @test tzt_a.mean  == tzt_b.mean
        @test tzt_a.scale == tzt_b.scale

        # ... and the perturbation really is capable of moving a fit, so the equality above is a
        # property of the SPLIT and not of an inert perturbation.
        zt_all, _ = fit_p12_transforms(Z2, TH2)
        @test zt_all.mean != zt_a.mean
    end

    @testset "mask rows are never standardized" begin
        # A standardized 0/1 mask re-couples folds through the mask mean -- the recorded reason in
        # `src/amortized/summary.jl:85-98`, mirrored by `standardize_p12`'s contract.
        n = P12_G^2
        Zraw = vcat(randn(n, 20) .+ 5.0, ones(n, 20))
        TH   = randn(P12_D_MINISPIKE, 20)
        zt, _ = fit_p12_transforms(Zraw, TH)
        out = standardize_p12(Zraw, zt)
        @test all(out[(n + 1):2n, :] .== 1.0)

        Zraw0 = vcat(randn(n, 20) .+ 5.0, zeros(n, 20))
        out0  = standardize_p12(Zraw0, zt)
        @test all(out0[(n + 1):2n, :] .== 0.0)

        # The continuous block, by contrast, IS transformed -- otherwise the two assertions above
        # would pass on a function that does nothing at all.
        @test out[1:n, :] != Zraw[1:n, :]
    end

    @testset "masking happens on RAW rows, before standardization (Pitfall 4)" begin
        # THE ORDERING TEST. A held-out region must be encoded byte-identically to a region
        # `correlation()` could not score: raw 0.0 with mask 0, which after per-row z-scoring is
        # (0 - mu_r)/sd_r -- NOT zero. This is the one assertion separating an in-distribution
        # held-out encoding from an out-of-distribution one, and it is why 12-14 departs from
        # `12-09-PLAN.md:69`'s "to the standardized input" wording.
        n = P12_G^2
        Zraw = vcat(randn(n, 40) .+ 5.0, ones(n, 40))   # train means are far from zero
        TH   = randn(P12_D_MINISPIKE, 40)
        zt, _ = fit_p12_transforms(Zraw, TH)

        r = 7
        masked_raw = mask_regions(Zraw, [r])            # the trainer's own masking semantics
        out = standardize_p12(masked_raw, zt)

        expected = (0.0 - zt.mean[r]) / zt.scale[r]
        @test out[r, 1] ≈ expected
        @test out[n + r, 1] == 0.0                      # the mask row itself stays a hard 0

        # ... and it is NOT zero, which is what a standardize-then-mask ordering would produce.
        @test !isapprox(out[r, 1], 0.0; atol = 1e-8)
        @test abs(zt.mean[r]) > 1.0                     # the fixture really does have a non-zero
                                                        # row mean, so the check above can bite
    end

    @testset "augment_mask! respects P12_MASK_K_SET and reports the realized rate" begin
        # NOTE this exercises 12-04's `augment_mask!` (p12_architecture.jl:190) -- the trainer
        # CALLS that function and deliberately does not redefine it. On an all-ones mask block
        # there are no natural zeros, so a recount of zeroed mask entries IS the drawn k, which is
        # what makes the k-set assertion below both true and falsifiable.
        n = P12_G^2
        ncol = 4000
        Zraw = vcat(randn(n, ncol) .+ 3.0, ones(n, ncol))
        rng = p12_rng(P12_MINISPIKE_COUNTER)
        Zaug = copy(Zraw)
        augment_mask!(Zaug, rng)

        ks = _p12train_masked_counts(Zaug)
        @test all(k -> k in P12_MASK_K_SET, ks)
        @test minimum(ks) == minimum(P12_MASK_K_SET)     # both endpoints occur over enough draws
        @test maximum(ks) == maximum(P12_MASK_K_SET)

        # No column may mask the same region twice -- i.e. the count of zeroed CONTINUOUS rows
        # matches the count of zeroed MASK rows, which a with-replacement draw would break.
        cont_zeroed = [count(iszero, view(Zaug, 1:n, c)) for c in 1:ncol]
        @test cont_zeroed == ks

        # The realized rate the trainer reports is a recount of this same quantity.
        @test mean(ks) / n ≈ mean(_p12train_masked_counts(Zaug)) / n

        # An unaugmented all-ones fixture has rate exactly 0, so the rate is measuring the
        # augmentation rather than the fixture.
        @test all(iszero, _p12train_masked_counts(Zraw))
    end

    @testset "the bundle round-trips with frozen transforms and matching field names" begin
        p = joinpath(mktempdir(), "smoke.jld2")
        b = train_p12_npe(arm = :car, n = 64, pool_dir = pool64, epochs = 1, batchsize = 16,
                          out_path = p, verbose = false)
        l = load_p12_npe(p)
        @test l.schema_version == 1
        @test length(l.theta_rows) == l.D
        @test l.zscore_arm === P12_ZSCORE_ARM
        # `harness.jl` consumers depend on these exact names; a rename breaks every one of them.
        @test haskey(l, :zt)
        @test haskey(l, :theta_zt)
        @test l.D == P12_D_MINISPIKE
        @test l.realized_mask_rate isa Real && 0 <= l.realized_mask_rate <= 1
        # The FROZEN transform round-trips bit-for-bit -- a refit on load would break this.
        @test l.zt.mean == b.zt.mean
        @test l.zt.scale == b.zt.scale
    end

    @testset "sampleposterior returns D x N from the trained smoke net" begin
        b = train_p12_npe(arm = :car, n = 64, pool_dir = pool64, epochs = 1, batchsize = 16,
                          verbose = false)
        pool = load_p12_pool(pool64)
        Z1 = reshape_summary(standardize_p12(Float64.(pool.summary_min)[:, 1:1], b.zt), P12_G)
        # `N` IS A KEYWORD IN NeuralEstimators v0.2.1, NOT POSITIONAL. The positional form
        # `sampleposterior(est, Z, 32)` resolves to a DIFFERENT three-argument method and throws a
        # MethodError. Verified against the installed signature at
        # `NeuralEstimators/src/Estimators/PosteriorEstimator.jl:138`, not from memory — this is
        # the young-API hazard CLAUDE.md names, and it is the second v0.2.1 default/signature
        # surprise in this phase after `use_gpu` defaulting to true.
        smp = sampleposterior(b.estimator, Z1; N = 32, use_gpu = false)
        M = smp isa AbstractVector ? first(smp) : smp
        @test size(M, 1) == b.D
        @test size(M, 2) == 32
        @test all(isfinite, M)
    end

    @testset "the Stage-1 gate guard blocks a real spend and permits the smoke path" begin
        # EVERY call here passes `verdict_dir = td`, so the real .planning/ tree is never read and
        # never written. A testset that omitted it would silently test the repository's actual
        # verdict state and pass or fail by accident.
        td = mktempdir()

        # :absent -- a real spend is refused, and the message names both missing paths.
        @test_throws Exception train_p12_npe(arm = :car, epochs = 18, n = 10_000, verdict_dir = td)
        err = try
            train_p12_npe(arm = :car, epochs = 18, n = 10_000, verdict_dir = td)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("p12_stage1_report.jld2", err)
        @test occursin("12-STAGE1-VERDICT.md", err)

        # The smoke path is exempt and must NOT throw with the same empty verdict dir.
        @test (train_p12_npe(arm = :car, n = 64, pool_dir = pool64, epochs = 1, batchsize = 16,
                             verdict_dir = td, verbose = false); true)

        # :descope -- the arm asymmetry 12-11's routing requires.
        write(joinpath(td, "12-STAGE1-VERDICT.md"), "VERDICT: DESCOPE\n")
        @test p12_stage1_verdict(; dir = td) === :descope
        @test_throws Exception train_p12_npe(arm = :car, epochs = 18, n = 10_000, verdict_dir = td)
        # `arm = :none` is permitted -- it gets PAST the guard and fails later, on the absent
        # :none pool, which is a DIFFERENT error. Asserting the message is what distinguishes
        # "permitted by the guard" from "blocked by the guard".
        derr = try
            train_p12_npe(arm = :none, epochs = 18, n = 10_000, verdict_dir = td)
            ""
        catch e
            sprint(showerror, e)
        end
        @test !occursin("may not be trained", derr)      # NOT the guard's refusal
        @test occursin("no COMPLETE pool", derr)         # the pool check, i.e. past the guard

        # :proceed -- a real spend is permitted through the guard.
        write(joinpath(td, "12-STAGE1-VERDICT.md"), "VERDICT: PROCEED\n")
        @test p12_stage1_verdict(; dir = td) === :proceed
        perr = try
            train_p12_npe(arm = :car, epochs = 18, n = 10_000, verdict_dir = td)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("no COMPLETE pool", perr)         # past the guard, stopped by the pool
    end

    @testset "the trainer's structural contracts, asserted on the source" begin
        # Exactly ONE standardizer fit site in the whole file.
        @test count("fit(StatsBase.ZScoreTransform", _P12TRAIN_CODE) +
              count("fit(ZScoreTransform", _P12TRAIN_CODE) == 2   # zt and theta_zt, same function
        # THE TRAINER READS POOLS AND NEVER MAKES THEM.
        @test !occursin("generate_p12_pool", _P12TRAIN_CODE)
        # CPU-only, and no dependency added.
        @test !occursin("use_gpu = true", _P12TRAIN_CODE)
        @test occursin("use_gpu = false", _P12TRAIN_CODE)
        @test !occursin("Pkg.add", _P12TRAIN_CODE)
        # `augment_mask!` is CALLED, never redefined -- the 12-04 function must stay single-source.
        @test !occursin("function augment_mask!", _P12TRAIN_CODE)
        @test occursin("augment_mask!(Zraw_train_masked", _P12TRAIN_CODE)
        # ... and `p12_stage1_verdict` likewise belongs to Tier 1.
        @test !occursin("function p12_stage1_verdict", _P12TRAIN_CODE)
        # The mask -> standardize -> reshape ORDER, visible in one expression.
        @test occursin("reshape_summary(standardize_p12(Zraw_train_masked, zt)", _P12TRAIN_CODE)
        # `theta_rows` is written WITH the width argument, never bare.
        @test occursin("P12_THETA_ROWS(K_dev)", _P12TRAIN_CODE)
    end

    @testset "p12_truncate_theta is the one width operator and is a no-op at full width" begin
        th = randn(P12_D_MINISPIKE, 5)
        @test p12_truncate_theta(th, P12_K_DEV) == th          # provable no-op
        t8 = p12_truncate_theta(th, 8)
        @test size(t8, 1) == 1 + 8 + 7 + 1 == 17
        @test t8[1, :] == th[1, :]                              # c0 stays first
        @test t8[(end - 7):end, :] == th[(end - 7):end, :]      # the 7 nuisances + r1 stay last
        @test length(P12_THETA_ROWS(8)) == 17
        @test first(P12_THETA_ROWS(8)) === :c0 && last(P12_THETA_ROWS(8)) === :r1
        # A head WIDER than the pool layout is refused, not silently padded.
        @test_throws ArgumentError p12_truncate_theta(th, P12_K_DEV + 1)
        # ... and a pool whose theta is not the mini-spike width is refused too.
        @test_throws DimensionMismatch p12_truncate_theta(randn(17, 5), 8)
    end

    @testset "a truncated production head is self-consistent and loadable" begin
        # THE ASSERTION THAT WOULD HAVE CAUGHT THE MISSING TRUNCATION OPERATOR. A K_dev < 63
        # bundle must round-trip through the same loader every downstream consumer uses.
        pt = joinpath(mktempdir(), "prod.jld2")
        bt = train_p12_npe(arm = :car, n = 64, K_dev = 8, pool_dir = pool64, epochs = 1,
                           batchsize = 16, out_path = pt, verbose = false)
        @test bt.D == 17
        @test bt.K_dev == 8
        @test length(bt.theta_rows) == bt.D
        lt = load_p12_npe(pt)
        @test lt.D == 17
        @test length(lt.theta_rows) == 17
    end

    @testset "pool_indices — the default is unchanged, and a subset is leak-free" begin
        # WHY THIS KEYWORD EXISTS. A scoring block cannot come from a SECOND pool:
        # `generate_p12_sample(idx)` keys on the GLOBAL INDEX alone and `generate_p12_pool(n)`
        # always generates 1:n, so two pools at the same arm/imsize share samples 1:min(n,m)
        # byte-identically. Asserted here rather than argued, because the whole leak-free
        # construction rests on it.
        s7a = generate_p12_sample(7; arm = :car, imsize = (128, 128))
        s7b = generate_p12_sample(7; arm = :car, imsize = (128, 128))
        s8  = generate_p12_sample(8; arm = :car, imsize = (128, 128))
        @test s7a.summary_min == s7b.summary_min      # index-keyed, so reproducible ...
        @test s7a.summary_min != s8.summary_min       # ... and index-DISTINCT
        # => a "separate scoring pool" is a SUPERSET, never a held-out set.

        # The DEFAULT path is provably unchanged: nothing selected means the whole pool.
        b_all = train_p12_npe(arm = :car, n = 64, pool_dir = pool64, epochs = 1, batchsize = 16,
                              verbose = false)
        @test b_all.pool_indices_given === nothing
        @test sort(vcat(b_all.train_indices, b_all.val_indices)) == collect(1:64)

        # A SUBSET leaves the complement untouched, which is what makes it a held-out set.
        held_out = 1:16                     # the caller's scoring block (the pool HEAD)
        remainder = 17:64
        b_sub = train_p12_npe(arm = :car, n = 64, pool_dir = pool64, pool_indices = remainder,
                              epochs = 1, batchsize = 16, verbose = false)
        seen = vcat(b_sub.train_indices, b_sub.val_indices)
        @test sort(seen) == collect(remainder)
        @test isempty(intersect(seen, held_out))          # THE LEAK CHECK
        @test isempty(intersect(b_sub.train_indices, b_sub.val_indices))
        # The val block is contaminated too (it drives early stopping), so a scorer must exclude
        # BOTH — which is why both sets are recorded.
        @test !isempty(b_sub.val_indices)

        # THE LEAK CHECK IS A POSITIVE CONTROL FOR ITSELF. An `isempty(intersect(...))` that has
        # never been shown non-empty is an assurance, not a test -- and this phase has already
        # found six verifies that could not fail. Hand the SAME expression a deliberately
        # overlapping selection and confirm it fires.
        b_leak = train_p12_npe(arm = :car, n = 64, pool_dir = pool64, pool_indices = 1:64,
                               epochs = 1, batchsize = 16, verbose = false)
        seen_leak = vcat(b_leak.train_indices, b_leak.val_indices)
        @test !isempty(intersect(seen_leak, held_out))   # FIRES when the block IS overlapped
        @test length(intersect(seen_leak, held_out)) == length(held_out)

        # Malformed selections are refused, not silently repaired.
        @test_throws Exception train_p12_npe(arm = :car, n = 64, pool_dir = pool64,
                                             pool_indices = [1, 1, 2], epochs = 1, batchsize = 16,
                                             verbose = false)
        @test_throws Exception train_p12_npe(arm = :car, n = 64, pool_dir = pool64,
                                             pool_indices = 60:70, epochs = 1, batchsize = 16,
                                             verbose = false)
    end

    @testset "training ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

    # =====================================================================================
    # 12-15's SELECTION RULE, exercised on SYNTHETIC score tables.
    #
    # These were an ACCEPTANCE CRITERION of 12-15 Task 1 that commit 4287862 did not actually
    # satisfy: its message claimed all five branches were exercised, but nothing in spike/test/
    # referenced `select_prior`, so the branches had been demonstrated ad hoc in a shell command
    # and then lost. A demonstration that is not persisted is not a test.
    #
    # TWO OF THESE CASES ARE THE ONLY EVIDENCE THE MACHINERY WOULD ACT AGAINST THE MORE
    # ATTRACTIVE NUMBER, which is the entire reason a rule written in advance is worth anything:
    #   - the TIE case, where GP's RMSE is NOMINALLY LOWER and the answer is still CAR;
    #   - the DISQUALIFICATION case, where GP has the LOWER RMSE and loses anyway for being
    #     miscalibrated.
    # A pre-declared tie-break that only ever fires when it agrees with the numbers is not a
    # tie-break, it is a preference wearing a rule's name. In the real 2026-07-31 run neither
    # case arose -- GP was both miscalibrated AND worse -- so these synthetic tables are where
    # that property is actually established.
    # =====================================================================================
    @testset "select_prior — the pre-registered rule, all branches (12-15)" begin
        # Admissible band is P12_COVERAGE_NOMINAL ± P12_STAGE2_COVERAGE_TOST_DELTA = [0.87, 0.93].
        _sc(; car, gp, none, n = 1000) = Dict(
            :car  => (coverage = car[1],  rmse = car[2],  n_test = n),
            :gp   => (coverage = gp[1],   rmse = gp[2],   n_test = n),
            :none => (coverage = none[1], rmse = none[2], n_test = n))

        # -- both admissible, CAR genuinely better -------------------------------------
        r = select_prior(_sc(car = (0.90, 0.50), gp = (0.90, 0.70), none = (0.90, 0.90)))
        @test r.selected == "CAR"
        @test r.beats_ablation
        @test Set(r.admissible) == Set([:car, :gp])

        # -- both admissible, GP genuinely better: the rule CAN pick GP ------------------
        # Without this case the tie-break test below would be consistent with a rule that
        # simply always answers CAR.
        r = select_prior(_sc(car = (0.90, 0.70), gp = (0.90, 0.50), none = (0.90, 0.90)))
        @test r.selected == "GP"
        @test r.beats_ablation

        # -- EVIDENTIAL CASE 1: a TIE resolves to CAR even though GP's RMSE is LOWER ------
        # 0.4975 vs 0.5000 is a 0.5 % gap, inside the pre-declared 1 % tie band.
        r = select_prior(_sc(car = (0.90, 0.5000), gp = (0.90, 0.4975), none = (0.90, 0.90)))
        @test r.selected == "CAR"
        @test occursin("TIE", r.reason)

        # -- EVIDENTIAL CASE 2: a MISCALIBRATED arm with the LOWER RMSE is DISQUALIFIED ---
        # GP's 0.30 beats CAR's 0.70 by a wide margin and is still not selected, because
        # coverage 0.99 is outside the band. "Better" is not a property a miscalibrated
        # posterior can have.
        r = select_prior(_sc(car = (0.90, 0.70), gp = (0.99, 0.30), none = (0.90, 0.90)))
        @test r.selected == "CAR"
        @test r.admissible == [:car]

        # -- both miscalibrated -> NONE-BEATS-ABLATION, regardless of RMSE ----------------
        r = select_prior(_sc(car = (0.99, 0.10), gp = (0.60, 0.10), none = (0.90, 0.90)))
        @test r.selected == "NONE-BEATS-ABLATION"
        @test !r.beats_ablation
        @test isempty(r.admissible)

        # -- the ablation wins -> NONE-BEATS-ABLATION, not "pick a winner anyway" ---------
        r = select_prior(_sc(car = (0.90, 0.95), gp = (0.90, 0.97), none = (0.90, 0.90)))
        @test r.selected == "NONE-BEATS-ABLATION"
        @test !r.beats_ablation
        @test occursin("UNSUPPORTED", r.reason)

        # -- the n-guard: the 0.03 tolerance may not be applied at an n it was not derived
        #    for. P12_STAGE2_N_MIN = 271 DATASETS.
        @test_throws Exception select_prior(
            _sc(car = (0.90, 0.50), gp = (0.90, 0.70), none = (0.90, 0.90), n = 100))

        # -- a missing arm is refused, not defaulted -------------------------------------
        @test_throws Exception select_prior(Dict(
            :car => (coverage = 0.90, rmse = 0.5, n_test = 1000),
            :gp  => (coverage = 0.90, rmse = 0.7, n_test = 1000)))

        # -- PURITY: same input, same answer ---------------------------------------------
        s = _sc(car = (0.90, 0.50), gp = (0.91, 0.70), none = (0.90, 0.90))
        @test select_prior(s).selected == select_prior(s).selected
        @test select_prior(s).reason == select_prior(s).reason
    end

    # =====================================================================================
    # PITFALL 5 REGRESSION GUARD. The 2026-07-31T09:51Z mini-spike run scored STANDARDIZED
    # posterior draws against RAW truth because it never applied `bundle.theta_zt`, and it
    # produced a report that passed every structural check -- index disjointness, per-region
    # pooling, budget accounting -- because none of those look at SCALE. It was caught only by
    # noticing the winning arm scored 0.45 % better than PREDICTING A CONSTANT ZERO. Its
    # `selected_prior`, `k_prod_recommended` and `n_low_recommended` were all invalid, and both
    # are append-only Tier-2 constants, so harvesting it would have been permanent.
    # =====================================================================================
    @testset "the mini-spike un-standardizes theta before scoring (Pitfall 5)" begin
        src = _p12train_strip_comments(
            read(joinpath(@__DIR__, "..", "validation", "run_p12_minispike.jl"), String))
        # The inverse transform is applied, to the BUNDLE'S OWN frozen theta_zt.
        @test occursin("StatsBase.reconstruct(bundle.theta_zt", src)
        # And the raw `sampleposterior` output is not what gets scored.
        @test occursin("sampleposterior", src)

        # The scoring path reports a null-predictor baseline, so "the net learned nothing"
        # is VISIBLE in the artifact rather than looking like a measurement.
        @test occursin("trivial_rmse", src)
        @test occursin("skill", src)

        # REPORTED IS NOT GATED, asserted on the rule's own body rather than in prose. `skill`
        # was added AFTER the numbers existed, so it must not be able to reach the selection --
        # that would be exactly the post-hoc bar the two-tier pre-registration forbids. This
        # FIRES if anyone ever wires it in.
        i = findfirst("function select_prior(scores)", src)
        @test i !== nothing
        body = src[first(i):last(findnext("\nend", src, last(i)))]
        @test !occursin("skill", body)
        @test !occursin("trivial", body)
    end

    @testset "the mini-spike inverts the SMOOTHNESS PERMUTATION before the iDCT (Pitfall 6)" begin
        # THIS TESTSET IS THE FALSIFIABLE HALF OF THE 2026-07-31 REPAIR. Its first block goes RED
        # the moment the scorer reverts to a bare `p12_idct_vec` on a theta coefficient column, and
        # the numbers below are what make it red: they are measured, not asserted.
        #
        # WHY A NUMERIC TEST AND NOT ONLY A SOURCE-TEXT ONE. `run_p12_minispike.jl` legitimately
        # CONTAINS `p12_idct_vec` -- the truncation curve uses it correctly, in FLAT index space --
        # so "does the file mention the right function" cannot separate the correct use from the
        # wrong one. Only arithmetic can.

        # --- 1. The two reconstructions are DIFFERENT OBJECTS, and only one is the inverse. ------
        # Built from a REAL prior draw through the REAL assembler, so this exercises the actual
        # generation path rather than a hand-made coefficient vector.
        draw  = sample_p12_prior(p12_datagen_rng(7); arm = :car)
        th    = p12_theta_column(draw)
        truth = vec(draw.z_field)
        nreg  = P12_G^2

        correct = vec(p12_region_field(th[1:nreg]; G = P12_G))   # inverts the permutation
        naive   = p12_idct_vec(th[1:nreg]; G = P12_G)            # does NOT -- the defect

        @test maximum(abs, correct .- truth) < 1e-12             # ~3.1e-15 measured
        @test maximum(abs, naive   .- truth) > 1.0               # ~3.72 measured

        # The consequence, stated as the quantity that actually misled the run: against a SCRAMBLED
        # truth even a perfect posterior scores WORSE than predicting a constant zero.
        naive_rmse   = sqrt(mean(abs2, naive .- truth))
        trivial_rmse = sqrt(mean(abs2, truth))
        @test naive_rmse / trivial_rmse > 1.2                    # 1.302 measured

        # The permutation is not the identity -- if it ever became so, the two paths would agree
        # and the tests above would pass VACUOUSLY.
        @test p12_dct_order(P12_G) != collect(1:nreg)
        @test p12_dct_order(P12_G)[1] == 1                       # why c0 survived the defect

        # --- 2. The scorer calls the ONE inverse and does not re-implement it. -------------------
        src = _p12train_strip_comments(
            read(joinpath(@__DIR__, "..", "validation", "run_p12_minispike.jl"), String))
        @test occursin("p12_region_field(", src)
        # The scoring loop must not scatter coefficients itself: a second inverse in this repo is
        # what produced the defect in the first place.
        fn_start = findfirst("function _p12ms_score_arm", src)
        @test fn_start !== nothing
        fn_end = findnext("\nend", src, last(fn_start))
        score_body = src[first(fn_start):last(fn_end)]
        @test !occursin("p12_dct_order", score_body)

        # --- 3. The head-width guard (H1) exists and names what it protects. ---------------------
        # Rows 1:64 are the field ONLY at the full-rank head; at K_dev < 63 they run into the
        # nuisances. Asserted on the source because the run that would trip it is 12-17's, not one
        # this suite can afford to execute.
        @test occursin("bundle.K_dev == P12_K_DEV", src)
        @test occursin("bundle.D == P12_D_MINISPIKE", src)
    end

    # =====================================================================================
    # THE θ-SPACE CONTRACT, MADE MECHANICAL — a SWEEP over every Phase-12 reporting runner.
    #
    # WHY A SWEEP AND NOT ONE MORE PER-FILE TEST. The contract already existed in prose at
    # `spike/npe/infer.jl:35-38`, named itself Pitfall 5, and said in as many words that reading
    # a standardized draw as raw "would be silently wrong". It was honoured by `infer.jl:110,122`
    # and `benchmark.jl:188` -- and violated, silently, by the ONE consumer that did not read it.
    # Prose did not hold. This is the same class as the repo's include-guard sweep: solve it once,
    # structurally, rather than once per plan.
    #
    # SCOPE IS `spike/validation/run_p12_*.jl` -- the phase's REPORTING RUNNERS -- and the scope is
    # chosen, not incidental. A blanket "every file calling sampleposterior must call reconstruct"
    # would be RED on files that are correct: `test_p12_architecture.jl:264` asserts only
    # `size(draws)`, and Phase 11's `test_lambda_ablation.jl` reports RATIOS of two posterior sds
    # of the same row, which an affine ZScoreTransform leaves invariant. Scoring against RAW TRUTH
    # is the thing that needs the inverse, and the reporting runners are where that happens.
    #
    # AT THE TIME OF WRITING THIS COVERS EXACTLY ONE FILE, and that is the point: `run_p12_minispike.jl`
    # is the only run_p12_* runner that samples a posterior at all -- the other four (ell_ridge,
    # eps_ridge, sim02, stage1_ridge) are ridge estimators -- which is precisely why the defect had
    # no sibling to be caught against. 12-16, 12-18, 12-19 and 12-20 add runners that DO score
    # posteriors, and they produce this phase's SBC and coverage claims. This sweep covers them the
    # moment they exist, and goes red until they un-standardize.
    # =====================================================================================
    @testset "EVERY p12 reporting runner un-standardizes theta and reports a null baseline" begin
        vdir = joinpath(@__DIR__, "..", "validation")
        runners = sort(filter(f -> startswith(f, "run_p12_") && endswith(f, ".jl"),
                              readdir(vdir)))

        # THE SWEEP MUST NOT BE VACUOUS. If a rename ever made the glob match nothing, or match
        # only files that never sample a posterior, this testset would pass by covering NOTHING --
        # the "0 tests reported as a pass" failure 12-05 already hit once in this phase.
        @test !isempty(runners)
        scoring = filter(runners) do f
            occursin("sampleposterior", _p12train_strip_comments(read(joinpath(vdir, f), String)))
        end
        @test !isempty(scoring)
        @test "run_p12_minispike.jl" in scoring

        for f in scoring
            code = _p12train_strip_comments(read(joinpath(vdir, f), String))
            # PITFALL 5: the estimator emits STANDARDIZED theta. A runner that samples a posterior
            # and compares it to raw truth MUST invert the bundle's frozen theta transform first.
            # THIS HALF IS UNIVERSAL -- it applies to every runner that scores a posterior at all,
            # whatever statistic it reports.
            @test occursin("reconstruct(", code)

            # THE NULL BASELINE, REPORTED-ONLY, AND CONDITIONAL **BY DESIGN**. An RMSE reported alone
            # is uninterpretable; reported beside what a CONSTANT predictor achieves, it is
            # self-checking. That is the diagnostic that turned the wrong-space bug from silent into
            # obvious -- GP scoring 60 % WORSE than predicting zero is impossible to overlook, but
            # only if something compares.
            #
            # IT IS REQUIRED ONLY OF RUNNERS THAT ACTUALLY REPORT AN RMSE, and the condition is the
            # whole point. An unconditional demand would go RED on runners that are CORRECT: 12-16's
            # coverage runner reports coverage and a log score and has no RMSE at all, and 12-18's SBC
            # runner reports RANKS, whose null -- uniformity -- is already built into the statistic.
            # Demanding an RMSE baseline where there is no RMSE would create a requirement nobody can
            # satisfy: the same shape as the defect this sweep exists to catch, pointing the other way.
            # A rule that fires on correct code teaches people to suppress it.
            occursin("rmse", lowercase(code)) && @test occursin("trivial_rmse", code)
        end
    end
end
