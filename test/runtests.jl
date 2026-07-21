#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
import Images
import Statistics: cor, mean
import StatsBase
import JLD2
import Pkg
import Random

using ProteinCoLoc
using Random123
using Test
Random123.seed!(1234)

# Test fixtures are referenced by package-root-relative paths (e.g.
# "test/test_images/..."). Under `Pkg.test()` the process CWD is not guaranteed to be the
# package root, so anchor it explicitly to `<pkgroot>` (the parent of this test dir).
cd(dirname(@__DIR__))

##########################################################################################
### CO-RESOLUTION HARD GATE (Phase 7 / Finding 1)  — must run FIRST
###
### Adding NeuralEstimators/Flux to the root package alongside the legacy plot/Turing stack
### historically capped NeuralEstimators below 0.2.1 and silently DOWNGRADED it to 0.1.4,
### breaking the v0.2.1 API. This gate asserts the RESOLVED environment pins the correct
### versions (installed set via `Pkg.dependencies()`, mirroring the spike's resolve-risk
### gate — NOT a Manifest regex). If this fails, no amortized inference code can be trusted.
##########################################################################################
@testset "co-resolution gate (Finding 1)" begin
    deps = Pkg.dependencies()
    function _installed_version(name)
        for (_, info) in deps
            info.name == name && return info.version
        end
        return nothing
    end

    ne = _installed_version("NeuralEstimators")
    fl = _installed_version("Flux")

    @test ne !== nothing            # NeuralEstimators must be a hard dependency
    @test fl !== nothing            # Flux must be a hard dependency
    # The Finding-1 downgrade target: NeuralEstimators must NOT be 0.1.4.
    @test string(ne) == "0.2.1"
    @test string(fl) == "0.16.10"
end

##########################################################################################
### D-02 result-type hierarchy + shared accessor interface
##########################################################################################
@testset "D-02 result hierarchy" begin
    @test isabstracttype(ProteinCoLoc.AbstractColocResult)
    @test ProteinCoLoc.AmortizedColocResult <: ProteinCoLoc.AbstractColocResult

    # The internal Turing result type must NOT be exported (D-01).
    @test !(:AdviColocResult in names(ProteinCoLoc))
    # The Turing-path entry points must NOT be exported (breaking release, D-01).
    @test !(:colocalization in names(ProteinCoLoc))
    @test !(:compute_BayesFactor in names(ProteinCoLoc))

    # Accessor interface is defined on the supertype and overridden by the shipped subtype.
    cal = ProteinCoLoc.CalibrationMeta(
        [0.5], [0.5], [0.5], [1], 0.0, 0.0, 8, (; seed = 0, passed = true))
    ood = ProteinCoLoc.OODVerdict(0.1, false, (; density = 0.1))
    post = reshape(collect(1.0:14.0), 7, 2)      # 7×N physical-θ draws
    r = ProteinCoLoc.AmortizedColocResult(
        8, post, [0.2, 0.3], -1.5, ood, cal, (; N = 2))

    @test delta_rho(r) == [0.2, 0.3]
    @test bayes_factor(r) == -1.5
    @test is_ood(r) == false
    @test posterior_draws(r) === post

    # The supertype fallback errors for an unimplemented subtype.
    struct _Unimpl <: ProteinCoLoc.AbstractColocResult end
    @test_throws ErrorException delta_rho(_Unimpl())
end

##########################################################################################
### Test for Image Loading functionality
##########################################################################################
@testset "LoadImages" verbose = true begin
    path = ["test/test_images/positive/positive_c1.tif", "test/test_images/positive/positive_c2.tif", "test/test_images/positive/positive_c3.tif"]
    # define a test function for load_tiff
    @testset "load_tiff" begin
        img = ProteinCoLoc.load_tiff(path[1])
        @test size(img) == (1028, 1376)

        # check that image is converted to grayscale
        img_ref = Images.load(path[1])
        img_ref = Images.Gray.(img_ref)
        @test img == img_ref

        # check that image is converted to matrix
        @test typeof(img) == Matrix{Float64}
    end

    # define a test set for the MultiChannelImage constructor
    @testset "MultiChannelImage constructor" begin
        name = "test_image"
        channels = ["blue", "green", "red"]

        # test that the constructor returns a MultiChannelImage object
        img = MultiChannelImage(name, path, channels)
        @test typeof(img) == MultiChannelImage{Float64, String, Float64, Int64}

        # test that the constructor sets the fields correctly
        @test img.channels == channels
        @test img.name == name
        @test img.path == path
        @test img.pixel_size == (1028, 1376)
        @test length(img.otsu_threshold) == size(channels)[1]
    end

    # test with five channel image
    @testset "MultiChannelImage constructor (5 channels)" begin
        name = "test_image"
        channels = ["blue", "green", "red","red_2","red_3"]
        img = MultiChannelImage(name, [path[1],path[2],path[3],path[3],path[3]],channels)
        @test typeof(img) == MultiChannelImage{Float64, String, Float64, Int64}

        # test that the constructor sets the fields correctly
        @test img.channels == channels
        @test img.name == name
        @test img.pixel_size == (1028, 1376)
        @test length(img.otsu_threshold) == size(channels)[1]
    end

    # test for mask calculation
    @testset "calculate and apply mask" begin
        # test with image where all pixels are below the threshold
        # value of mask should be 0 everywhere
        img = MultiChannelImage(
            fill(fill(1.0, 128, 128),3),
            ["blue", "green", "red"],
            "all_below_otsu",
            path,(128, 128),[2.5, 2.5, 2.5]
            )

        @test ProteinCoLoc._calculate_mask(img) == fill(fill(0, 128, 128),3)
        img = ProteinCoLoc._apply_mask!(img, ProteinCoLoc._calculate_mask(img))
        @test img.data == fill(fill(0, 128, 128),3)

        # test with image where all pixels are above the threshold
        # value of mask should be 1 everywhere
        img = MultiChannelImage(
            fill(fill(1.0, 128, 128),3),
            ["blue", "green", "red"],
            "all_below_otsu",
            path,(128, 128),[0.5, 0.5, 0.5]
            )

        @test ProteinCoLoc._calculate_mask(img) == fill(fill(1, 128, 128),3)
        img = ProteinCoLoc._apply_mask!(img, ProteinCoLoc._calculate_mask(img))
        @test img.data == fill(fill(1, 128, 128),3)
        # confirm that the image size is not changed by applying the mask
        @test size(img.data[1]) == (128, 128)
        @test size(img.data[2]) == (128, 128)
        @test size(img.data[3]) == (128, 128)
    end

    @testset "MultiChannelImageStack" verbose = true begin
        path = ["test/test_images/positive/positive_c1.tif", "test/test_images/positive/positive_c2.tif", "test/test_images/positive/positive_c3.tif"]
        # load image
        img = MultiChannelImage("test_image", path, ["blue", "green", "red"])
        # create stack from image and test that the constructor returns a MultiChannelImageStack object
        img_stack = MultiChannelImageStack([img, img, img], "test_stack")
        @test typeof(img_stack) == MultiChannelImageStack{MultiChannelImage{Float64, String, Float64, Int64}, String}

        # test that the constructor sets the fields correctly
        @test img_stack.name == "test_stack"
        @test length(img_stack) == 3

        # test that indexing works as expected
        @test img_stack[1] == img
        @test_throws BoundsError img_stack[4] # check that indexing out of bounds throws an error

        # test that iterations over the stack work as expected
        for (image, idx, i) ∈ zip(img_stack, 1:length(img_stack), 1:3)
            @test image == img
            @test idx == i
        end
    end
end

##########################################################################################
### Test for patching functionality
##########################################################################################
@testset "patch" verbose = true begin
    # define a test case and call patch function
    img = rand(100, 100)  # Create a 100x100 matrix of random numbers
    num_patches = 10  # We want to divide the image into 10x10 patches
    patches = ProteinCoLoc.patch(img, num_patches)

    # check that the output is a 4D array
    @test typeof(patches) == Array{Union{Float64, Missing}, 4}

    # check that the output has the expected dimensions
    @test size(patches) == (num_patches, num_patches, size(img, 1) ÷ num_patches, size(img, 2) ÷ num_patches)

    # check that the values in the patches match the corresponding values in the original image
    patches_correct::Vector{Bool} = []
    for i in 1:num_patches
        for j in 1:num_patches
            push!(
                patches_correct,
                patches[i, j, :, :] == img[(i-1)*size(img, 1) ÷ num_patches+1:i*size(img, 1) ÷ num_patches, (j-1)*size(img, 2) ÷ num_patches+1:j*size(img, 2) ÷ num_patches]
                )
        end
    end
    @test all(patches_correct)
end

###########################################################################################
### Test for correlation calculation
###########################################################################################
@testset "Colocalization" verbose = true begin
    # define a test set for the _exclude_zero! function
    @testset "_exclude_zero" begin
        # define a test case
        a = [1, 2, 0, 3, 4, 0, 5, missing]
        b = [0, 2, 3, 0, 4, 5, 6, missing]

        # _exclude_zero returns the filtered vectors (it does not mutate in place)
        a, b = ProteinCoLoc._exclude_zero(a, b)

        # check that the output is as expected
        @test a == [2,4,5]
        @test b == [2,4,6]
    end

    # define a test set for the correlation function
    @testset "correlation" begin
        # define a test case by defining two 5x5x10x10 4D arrays of random numbers
        # corresponding to two images with 10x10 patches
        x = rand(5, 5, 10, 10)  # Create a 5x5x10x10 4D array of random numbers
        y = rand(5, 5, 10, 10)  # Create another 5x5x10x10 4D array of random numbers

        # call correlation function
        ρ = ProteinCoLoc.correlation(x, y)

        # check that the output is a 2D array
        @test typeof(ρ) == Array{Union{Float64, Missing}, 2}
        # check that the output has the expected dimensions
        @test size(ρ) == (5, 5)

        # check that the correlation is calculated correctly
        for i in 1:5
            for j in 1:5
                a = x[i, j, :, :][:]
                b = y[i, j, :, :][:]
                @test ρ[i, j] == cor(a, b)
            end
        end
    end
end

##########################################################################################
### Amortized grid-parametric summary/encoder + dimension helpers (07-01 Task 1)
###
### D-04 couples the summary-vector dimension to the patch grid G (dim 2·G²). These tests
### assert the summary/encoder and the single-source dimension helpers are grid-general for
### G ∈ {4,8,16,32} and reproduce the proven 8×8 spike constants (128 / 64 / 320).
##########################################################################################
@testset "amortized summary (grid-parametric)" begin
    sz = 256
    data = [rand(sz, sz) .+ 0.5, rand(sz, sz) .+ 0.5]   # strictly positive ⇒ clears ≥15 floor
    mci = MultiChannelImage(data, ["c1", "c2"], "synth", ["p1", "p2"], (sz, sz), [0.5, 0.5])

    # patch_summary(mci, G) is a G×G correlation matrix; encode_d01 has length 2·G².
    for G in (4, 8, 16)
        M = ProteinCoLoc.patch_summary(mci, G)
        @test size(M) == (G, G)
        @test length(ProteinCoLoc.encode_d01(M)) == 2 * G^2
    end

    # Row-partition split is grid-general: cont = 1:G², mask = G²+1:2G².
    for G in (4, 8, 16, 32)
        cont, mask = ProteinCoLoc._summary_row_partition(:min, 2 * G^2)
        @test cont == collect(1:G^2)
        @test mask == collect(G^2+1:2*G^2)
    end

    # Dimension helpers are the single source of truth and reproduce the 8×8 constants.
    @test ProteinCoLoc.summary_dim(8) == 128
    @test ProteinCoLoc.cont_rows(8) == 64
    @test ProteinCoLoc.ratio_input_dim(8) == 320
    @test ProteinCoLoc.summary_dim(16) == 512
    @test ProteinCoLoc.cont_rows(16) == 256
    @test ProteinCoLoc.ratio_input_dim(16) == 1280
end

##########################################################################################
### Grid-keyed estimator registry skeleton (07-01 Task 3, PROD-02 / D-04)
###
### The registry keys estimators by patch grid, returns a bundle for a registered grid, and
### raises a clear `train_and_register`-pointing error for an unregistered grid. The shipped
### family is CAPPED to (4,8,16,32) — 64 is DROPPED (D-04).
##########################################################################################
@testset "estimator registry (PROD-02)" begin
    # Capped shipped family — 64 dropped (D-04).
    @test ProteinCoLoc._SHIPPED_GRIDS == (4, 8, 16, 32)
    @test !(64 in ProteinCoLoc._SHIPPED_GRIDS)

    # Unregistered, non-shipped grid → ArgumentError naming train_and_register.
    @test_throws ArgumentError estimator_for(7)
    err = try
        estimator_for(7)
    catch e
        e
    end
    @test occursin("train_and_register", err.msg)

    # GPU probe is weakdep-safe: no CUDA loaded ⇒ false, never an error.
    @test ProteinCoLoc.has_cuda_device() == false

    # register! makes estimator_for return the bundle.
    cal = ProteinCoLoc.CalibrationMeta(
        [0.5], [0.5], [0.5], [1], 0.0, 0.0, 8, (; seed = 0, passed = true))
    b = ProteinCoLoc.EstimatorBundle(8, nothing, nothing, nothing, nothing, nothing, cal)
    register!(b)
    @test estimator_for(8) === b
    @test b.grid == 8
end

##########################################################################################
### Amortized NPE read surface (07-02 Task 1, PROD-01)
###
### standardize_summary applies the FROZEN zt to the continuous rows and passes the binary
### mask rows through UNCHANGED; the read functions carry a CPU-default (use_gpu=false).
##########################################################################################
@testset "amortized NPE read surface (infer.jl)" begin
    # The read functions are promoted and available on the module.
    for f in (:standardize_summary, :posterior_for, :rho_hat, :rho_draws, :delta_rho)
        @test isdefined(ProteinCoLoc, f)
    end

    # standardize_summary: frozen-zt on continuous rows, mask rows bypassed (G = 4 ⇒ d = 32).
    G = 4; nc = G^2; K = 64
    Scont = randn(nc, K) .* 3.0 .+ 7.0
    zt    = StatsBase.fit(StatsBase.ZScoreTransform, Scont; dims = 2)
    mask  = Float64.(rand(Bool, nc, K))
    S     = vcat(Scont, mask)                          # 2·G² × K
    out   = ProteinCoLoc.standardize_summary(S, zt, :min)

    @test eltype(out) == Float32
    @test size(out) == (2 * nc, K)
    # mask rows (G²+1 : 2G²) pass through UNCHANGED (never re-z-scored).
    @test out[nc+1:2*nc, :] ≈ Float32.(mask)
    # continuous rows are z-scored ⇒ per-row mean ≈ 0.
    @test all(abs.(vec(mean(Float64.(out[1:nc, :]); dims = 2)))  .< 1e-4)

    # Vector overload coerces to a column and matches the matrix result.
    outv = ProteinCoLoc.standardize_summary(S[:, 1], zt, :min)
    @test outv == out[:, 1]
    @test outv[nc+1:2*nc] ≈ Float32.(mask[:, 1])
end

##########################################################################################
### Amortized Bayes factor + non-clamped KDE baseline (07-02 Task 2, PROD-01 / T-7-06)
###
### pair_encode is the grid-general A7 difference encoding (length ratio_input_dim(G)=5G²);
### kde_log_bf_unclamped mirrors the compute_BayesFactor KDE math WITHOUT the 1e-8 floor, so
### max|Δ logBF| is measurable free of the clamp artifact.
##########################################################################################
@testset "amortized Bayes factor (bf.jl)" begin
    @test isdefined(ProteinCoLoc, :amortized_log_bf)
    @test isdefined(ProteinCoLoc, :pair_encode)
    @test isdefined(ProteinCoLoc, :kde_log_bf_unclamped)

    # pair_encode: length == ratio_input_dim(G) = 5·G², layout = vcat(Zs, Zc, contrast).
    for G in (4, 8, 16)
        d  = 2 * G^2
        Zs = rand(Float32, d); Zc = rand(Float32, d)
        enc = ProteinCoLoc.pair_encode(Zs, Zc)
        @test length(enc) == ProteinCoLoc.ratio_input_dim(G)
        @test enc[1:d] == Zs
        @test enc[d+1:2d] == Zc
        @test enc[2d+1:2d+G^2] ≈ Float32.(Zs[1:G^2] .- Zc[1:G^2])
    end
    # mismatched / odd-length inputs are rejected.
    @test_throws DimensionMismatch ProteinCoLoc.pair_encode(rand(Float32, 32), rand(Float32, 8))

    # Non-clamped KDE baseline: a FULLY one-sided posterior (all draws ≫ 0) with a prior that
    # straddles 0 drives p_post → 1 with NO 1e-8 floor, so |logBF| exceeds the clamped ceiling
    # log((1−1e-8)/1e-8) ≈ 18.42 (a clamped baseline would cap there). This is the artifact-free
    # measurement the ship-gate needs (Memo §5 / T-7-06).
    post  = fill(5.0, 400) .+ 0.01 .* randn(400)     # fully positive ⇒ one-sided
    prior = randn(4000)                              # straddles 0 ⇒ p_prior ≈ 0.5
    lb = ProteinCoLoc.kde_log_bf_unclamped([post], prior)
    @test length(lb) == 1
    @test lb[1] > 18.42                              # beyond the 1e-8 clamp ceiling ⇒ no floor

    # Numerical-domain safety: a fully one-sided NEGATIVE posterior (all draws ≪ 0) drives
    # p_post → 0. QuadGK's adaptive integral of the KDE can overshoot `p_le` a few ulp past 1.0,
    # so the un-guarded `1 - p_le` would be a tiny NEGATIVE probability and crash `log(·)` with a
    # DomainError. The [0,1] domain clamp must instead yield the honest boundary value (a large
    # negative or -Inf logBF, which the gate drops as non-finite) WITHOUT throwing.
    post_neg = fill(-5.0, 400) .+ 0.01 .* randn(400)   # fully negative ⇒ one-sided, p_post → 0
    lbn = ProteinCoLoc.kde_log_bf_unclamped([post_neg], prior)
    @test length(lbn) == 1
    @test !isnan(lbn[1])                              # no DomainError / NaN — domain-safe
    @test lbn[1] < -18.42 || isinf(lbn[1])           # honest one-sided value (large-neg or -Inf)
end

##########################################################################################
### Amortized OOD flag: density + noise + re-enabled posterior-predictive (07-02 Task 3)
###
### fit_ood_nulls restricts to the continuous rows via _summary_row_partition (1e-6·I ridge);
### the iter1 finite-guard keeps the PP channel computable on non-finite θ̂ (Memo §5 / T-7-04);
### ood_verdict OR-fuses the available channels and returns an OODVerdict.
##########################################################################################
@testset "amortized OOD flag (ood.jl)" begin
    for f in (:fit_ood_nulls, :maha_score, :noise_features, :fit_noise_null, :noise_score,
              :roc_auc, :id_threshold, :youden_j, :ood_flag, :ood_verdict,
              :pp_mismatch_score)
        @test isdefined(ProteinCoLoc, f)
    end

    # fit_ood_nulls: continuous rows only (1:G²) + 1e-6 ridge (Σ well-conditioned, cholesky OK).
    G = 4; nc = G^2; d = 2 * nc
    Zcont = randn(nc, 80)
    mask  = Float64.(rand(Bool, nc, 80))
    Ztr   = vcat(Zcont, mask)                          # d × 80 standardized-summary fixture
    nulls = ProteinCoLoc.fit_ood_nulls(Ztr; variant = :min)
    @test nulls.cont == collect(1:nc)
    @test length(nulls.μS) == nc
    s = ProteinCoLoc.maha_score(nulls, Ztr[:, 1])
    @test s ≥ 0 && isfinite(s)

    # PP finite-guard: a non-finite θ̂ maps to a finite, prior-valid tuple (no crash, no NaN).
    θ = ProteinCoLoc._theta_tuple([NaN, Inf, -Inf, NaN, 5.0, -3.0, NaN])
    @test all(isfinite, values(θ))
    @test -1.0 ≤ θ.ρ_true ≤ 1.0
    @test 0.0 ≤ θ.spillover ≤ 1.0
    @test θ.autofluorescence ≥ 0.0
    @test 0.0 ≤ θ.label_efficiency ≤ 1.0
    @test θ.noise ≥ 0.0
    @test θ.shift_dx == 5.0 && θ.shift_dy == -3.0

    # noise channel: 10 invariant features (5 per channel), all finite; null + score round-trip.
    pair = [randn(64, 64), randn(64, 64)]
    nf = ProteinCoLoc.noise_features(pair)
    @test length(nf) == 10
    @test all(isfinite, nf)
    F  = reduce(hcat, [ProteinCoLoc.noise_features([randn(64, 64), randn(64, 64)]) for _ in 1:30])
    nn = ProteinCoLoc.fit_noise_null(F)
    @test isfinite(ProteinCoLoc.noise_score(nn, pair))

    # hand-rolled ROC/AUC: separable ⇒ 1.0, identical ⇒ 0.5.
    @test ProteinCoLoc.roc_auc([0.0, 0.1, 0.2], [1.0, 1.1, 1.2])[3] == 1.0
    @test ProteinCoLoc.roc_auc([0.0, 1.0], [0.0, 1.0])[3] == 0.5

    # ood_verdict: density-only fusion (with_pp=false, no image pair) returns an OODVerdict.
    v = ProteinCoLoc.ood_verdict((; density = nulls), Ztr[:, 1]; with_pp = false)
    @test v isa ProteinCoLoc.OODVerdict
    @test isfinite(v.score)
    @test v.flag == false                              # no pre-registered threshold ⇒ no fire
    @test haskey(v.per_channel, :density)

    # With a pre-registered ID-quantile threshold the flag FIRES on an off-manifold summary.
    id_scores = [ProteinCoLoc.maha_score(nulls, Ztr[:, k]) for k in 1:80]
    thr = ProteinCoLoc.id_threshold(id_scores)
    off = 50.0 .* ones(d)                              # far off the train manifold ⇒ high density
    v2  = ProteinCoLoc.ood_verdict((; density = nulls, thr = thr), off; with_pp = false)
    @test v2.flag == true
end

##########################################################################################
### Amortized NPE architecture + training, GPU-plumbed with CPU fallback (07-03 Task 1)
###
### build_estimator is input-width-agnostic (reads d_in) and constructs q as a NormalisingFlow
### INSTANCE passed positionally; train_npe defaults use_gpu=has_cuda_device() (CPU here) with no
### throw guard, keeps the AdamW LR/decay Float64, and returns a PosteriorEstimator + frozen θzt.
##########################################################################################
@testset "amortized NPE training (train_npe.jl)" begin
    for f in (:build_estimator, :train_npe, :fit_theta_transform, :fit_summary_transform)
        @test isdefined(ProteinCoLoc, f)
    end

    # build_estimator is input-width-agnostic: two different d_in both construct an estimator
    # (a working estimator yields 7-row posterior draws through the frozen read surface).
    for d_in in (32, 72)
        e = ProteinCoLoc.build_estimator(d_in; dstar = 8, depth = 1, width = 16,
                                         num_coupling_layers = 2, flow_depth = 1, flow_width = 8)
        dz = ProteinCoLoc.posterior_for(e, Float32.(randn(d_in)); N = 4, use_gpu = false)
        @test size(dz, 1) == 7
    end

    # fit_summary_transform touches the CONTINUOUS rows only (dim G²), mask rows never z-scored.
    G = 4; nc = G^2; d = 2 * nc; n = 40
    Zraw = vcat(randn(nc, n) .* 2 .+ 3, Float64.(rand(Bool, nc, n)))
    zt   = ProteinCoLoc.fit_summary_transform(Zraw; variant = :min)
    @test length(zt.mean) == nc

    # A tiny CPU NPE train (few epochs, small N, use_gpu=false) completes and returns an estimator.
    Ztr = Float32.(ProteinCoLoc.standardize_summary(Zraw, zt, :min))
    Zva = Float32.(ProteinCoLoc.standardize_summary(Zraw, zt, :min))
    θtr = randn(7, n); θva = randn(7, n)
    res = ProteinCoLoc.train_npe(Ztr, Zva, θtr, θva; use_gpu = false, epochs = 2,
                                 batchsize = 16, dstar = 8, depth = 1, width = 16,
                                 num_coupling_layers = 2, flow_depth = 1, flow_width = 8,
                                 stopping_epochs = 2, verbose = false)
    @test res.d_in == d
    @test length(res.θzt.zt.mean) == 7          # bounded θ-space: z-score lives inside (F2)
    @test res.arch.dstar == 8
    # posterior draws are reachable through the frozen read surface (7×N un-standardized).
    draws = ProteinCoLoc.posterior_for(res.estimator, Ztr[:, 1]; N = 8, use_gpu = false)
    @test size(draws, 1) == 7
end

##########################################################################################
### Bounded θ-space (07-CALIBRATION-FINDINGS F2)
###
### The prior box is the SINGLE SOURCE OF TRUTH (derived from the prior objects `sample_prior`
### draws from), the logit bijection round-trips, and — the point of the fix — EVERY value the
### inverse map can produce lies inside the prior support, so a flow draw can never leak
### out-of-support the way the plain-ZScoreTransform θ-space did.
##########################################################################################
@testset "bounded θ-space (F2)" begin
    for f in (:theta_prior_bounds, :theta_to_unbounded, :theta_from_unbounded,
              :BoundedThetaTransform)
        @test isdefined(ProteinCoLoc, f)
    end

    b = ProteinCoLoc.theta_prior_bounds()
    @test length(b) == 7
    # SINGLE SOURCE OF TRUTH: every bound is the prior object's own support, not a literal.
    @test b[2] == (minimum(ProteinCoLoc.SPILLOVER_PRIOR), maximum(ProteinCoLoc.SPILLOVER_PRIOR))
    @test b[4] == (minimum(ProteinCoLoc.LABEL_EFFICIENCY_PRIOR),
                   maximum(ProteinCoLoc.LABEL_EFFICIENCY_PRIOR))
    @test b[7] == (minimum(ProteinCoLoc.NOISE_PRIOR), maximum(ProteinCoLoc.NOISE_PRIOR))
    # ρ_true = ghat(μ*) is a CLAMPED map ⇒ its exact range is the first/last ρ knot.
    @test b[1] == (first(ProteinCoLoc.GHAT_RHO_KNOTS), last(ProteinCoLoc.GHAT_RHO_KNOTS))
    @test all(lo < hi for (lo, hi) in b)

    # Real prior draws round-trip through the fitted transform (ε-clamp tolerance at the atoms).
    rngp = Random123.Philox4x(UInt64, (UInt64(20260721), UInt64(0)))
    θp   = hcat([collect(values(ProteinCoLoc.sample_prior(rngp))) for _ in 1:200]...)
    t    = ProteinCoLoc.fit_theta_transform(θp)
    @test t isa ProteinCoLoc.BoundedThetaTransform
    @test length(t.zt.mean) == 7
    Yp = StatsBase.transform(t, θp)
    @test all(isfinite, Yp)
    @test maximum(abs.(StatsBase.reconstruct(t, Yp) .- θp)) < 1e-4

    # THE GUARANTEE: arbitrary (even extreme) flow-space values map back INSIDE the prior box.
    wild = 50.0 .* randn(7, 500)
    back = StatsBase.reconstruct(t, wild)
    @test all(b[p][1] <= back[p, j] <= b[p][2] for p in 1:7, j in 1:size(back, 2))
    # Float32 draws (what `sampleposterior` returns) keep their element type.
    @test eltype(StatsBase.reconstruct(t, Float32.(wild))) == Float32
    # The map is strictly monotone per parameter ⇒ SBC ranks are invariant to the change of space.
    mono = StatsBase.reconstruct(t, repeat(collect(-6.0:0.5:6.0)', 7, 1))
    @test all(issorted(mono[p, :]) for p in 1:7)
end

##########################################################################################
### Amortized NRE ratio training, grid-general + GPU-plumbed with CPU fallback (07-03 Task 2)
###
### The conditioner input width derives from ratio_input_dim(G)=5G² (not literal 320); use_gpu
### defaults has_cuda_device() (CPU here) with no throw guard; no custom loss is passed.
##########################################################################################
@testset "amortized NRE training (train_ratio.jl)" begin
    for f in (:build_ratio_estimator, :train_ratio, :assemble_ratio_pairs, :measure_log_prior_odds)
        @test isdefined(ProteinCoLoc, f)
    end

    # build_ratio_estimator input width == ratio_input_dim(G) = 5G² (grid-general), not 320:
    # a working RatioEstimator yields a finite amortized log BF (one forward pass) at that width.
    for G in (4, 8)
        e   = ProteinCoLoc.build_ratio_estimator(ProteinCoLoc.ratio_input_dim(G);
                                                 num_summaries = 8, width = 16)
        Zp  = Float32.(randn(ProteinCoLoc.ratio_input_dim(G)))
        @test isfinite(ProteinCoLoc.amortized_log_bf(e, Zp, 0.0))
    end

    # assemble_ratio_pairs: grid-general pair encoding + balanced-by-construction labels.
    G = 4; nc = G^2; d = 2 * nc; N = 60
    Zstd = Float32.(randn(d, N))
    ρ    = randn(N)
    Zp, y = ProteinCoLoc.assemble_ratio_pairs(Zstd, ρ, 40)
    @test size(Zp, 1) == ProteinCoLoc.ratio_input_dim(G)
    @test all(v -> v == 0f0 || v == 1f0, y)
    @test isfinite(ProteinCoLoc.measure_log_prior_odds(y))

    # A tiny CPU ratio train (few epochs, small n, use_gpu=false) returns a usable handle.
    res = ProteinCoLoc.train_ratio(Zstd, ρ; n = 60, use_gpu = false, epochs = 2, batchsize = 16,
                                   val_frac = 0.2, stopping_epochs = 2, num_summaries = 8,
                                   summary_width = 16, verbose = false)
    @test res.input_dim == ProteinCoLoc.ratio_input_dim(G)
    @test isfinite(res.log_prior_odds)
    # the amortized log BF read runs on the trained handle (one forward pass, no quadgk/KDE).
    Zpair = ProteinCoLoc.pair_encode(Zstd[:, 1], Zstd[:, 2])
    @test isfinite(ProteinCoLoc.amortized_log_bf(res.estimator, Zpair, res.log_prior_odds))
end

##########################################################################################
### CPU-resident Flux.state persistence for estimators AND OOD nulls (07-03 Task 3, Pitfall 4)
###
### The persisted estimator JLD2 holds a Flux.state-derived key (NOT a whole PosteriorEstimator)
### + architecture metadata; load rebuilds the arch and applies Flux.loadmodel!. Save→load→CPU
### inference reproduces the pre-save output. OOD nulls round-trip through the same atomic wrapper.
##########################################################################################
@testset "CPU-resident persistence (persist.jl)" begin
    for f in (:save_estimator, :load_estimator, :save_ratio, :load_ratio,
              :save_ood_nulls, :load_ood_nulls)
        @test isdefined(ProteinCoLoc, f)
    end

    tmpdir = mktempdir()
    G = 4; nc = G^2; d = 2 * nc; n = 40
    Zraw = vcat(randn(nc, n) .* 2 .+ 3, Float64.(rand(Bool, nc, n)))
    zt   = ProteinCoLoc.fit_summary_transform(Zraw; variant = :min)
    Zstd = Float32.(ProteinCoLoc.standardize_summary(Zraw, zt, :min))
    θ    = randn(7, n)
    res  = ProteinCoLoc.train_npe(Zstd, Zstd, θ, θ; use_gpu = false, epochs = 2, batchsize = 16,
                                  dstar = 8, depth = 1, width = 16, num_coupling_layers = 2,
                                  flow_depth = 1, flow_width = 8, stopping_epochs = 2)

    # The persisted JLD2 holds a Flux.state-derived key, NOT a whole PosteriorEstimator object.
    npe_path = joinpath(tmpdir, "npe_4.jld2")
    ProteinCoLoc.save_estimator(npe_path, res.estimator, res.θzt, zt, res.arch)
    keys_on_disk = JLD2.load(npe_path)
    @test haskey(keys_on_disk, "model_state")
    @test haskey(keys_on_disk, "arch")
    @test !haskey(keys_on_disk, "estimator")           # NOT the whole object

    # Round-trip: save → load → identical CPU inference output (seeded ⇒ bitwise-equal draws).
    loaded = ProteinCoLoc.load_estimator(npe_path)
    Zq = Zstd[:, 1]
    Random123.seed!(4242)
    p_before = ProteinCoLoc.posterior_for(res.estimator, Zq; N = 32, use_gpu = false)
    Random123.seed!(4242)
    p_after  = ProteinCoLoc.posterior_for(loaded.estimator, Zq; N = 32, use_gpu = false)
    @test p_before ≈ p_after
    @test loaded.arch.dstar == 8
    @test loaded.zt.mean == zt.mean

    # Ratio handle round-trips through the CPU-resident wrapper.
    rres = ProteinCoLoc.train_ratio(Zstd, θ[1, :]; n = 60, use_gpu = false, epochs = 2,
                                    batchsize = 16, val_frac = 0.2, stopping_epochs = 2,
                                    num_summaries = 8, summary_width = 16)
    ratio_path = joinpath(tmpdir, "ratio_4.jld2")
    ProteinCoLoc.save_ratio(ratio_path, rres)
    rloaded = ProteinCoLoc.load_ratio(ratio_path)
    @test rloaded.input_dim == ProteinCoLoc.ratio_input_dim(G)
    @test rloaded.log_prior_odds == rres.log_prior_odds
    Zpair = ProteinCoLoc.pair_encode(Zstd[:, 1], Zstd[:, 2])
    @test ProteinCoLoc.amortized_log_bf(rres.estimator, Zpair, rres.log_prior_odds) ≈
          ProteinCoLoc.amortized_log_bf(rloaded.estimator, Zpair, rloaded.log_prior_odds)

    # OOD nulls (density + noise) round-trip through the same atomic .tmp→integrity→mv wrapper.
    dens  = ProteinCoLoc.fit_ood_nulls(Zstd; variant = :min)
    F     = reduce(hcat, [ProteinCoLoc.noise_features([randn(64, 64), randn(64, 64)]) for _ in 1:20])
    noise = ProteinCoLoc.fit_noise_null(F)
    nulls = (; density = dens, noise = noise)
    ood_path = joinpath(tmpdir, "ood_nulls_4.jld2")
    ProteinCoLoc.save_ood_nulls(ood_path, nulls)
    nloaded = ProteinCoLoc.load_ood_nulls(ood_path)
    @test nloaded.density.μS == dens.μS
    zq = Zstd[:, 1]
    @test ProteinCoLoc.maha_score(nloaded.density, zq) ≈ ProteinCoLoc.maha_score(dens, zq)
end

##########################################################################################
### Reusable per-grid pipeline wired into train_and_register (07-03 Task 4, PROD-01/02, D-04)
###
### _train_grid_pipeline factors datagen→train_npe→train_ratio→fit_ood_nulls→persist into ONE
### function returning an EstimatorBundle; train_and_register(grid) runs it then register!s. A
### tiny fixture run (injected synthetic datagen, no simulator) produces the three grid artifacts
### and a loadable bundle; a second call SKIPS-IF-DONE. Per-grid image-size defaults bias upward.
##########################################################################################
@testset "per-grid pipeline (pipeline.jl)" begin
    # Structural: the pipeline is a defined function and train_and_register has an Int method.
    @test isdefined(ProteinCoLoc, Symbol("_train_grid_pipeline"))
    @test hasmethod(ProteinCoLoc.train_and_register, Tuple{Int})
    @test isdefined(ProteinCoLoc, :default_imsize_for)
    @test isdefined(ProteinCoLoc, :default_npairs_for)

    # Per-grid image-size bias: 16 → ≥512², 32 → ≥1024² (finer grids draw larger images).
    s16 = ProteinCoLoc.default_imsize_for(16)
    @test all(sz -> min(sz[1], sz[2]) >= 512, s16)
    s32 = ProteinCoLoc.default_imsize_for(32)
    @test all(sz -> min(sz[1], sz[2]) >= 1024, s32)
    @test ProteinCoLoc.default_npairs_for(4) < ProteinCoLoc.default_npairs_for(32)

    # A tiny fixture run: inject a synthetic datagen (grid-4 ⇒ 2·G² = 32-row summary), tiny nets,
    # into a temp artifacts root — no forward simulator needed.
    root = mktempdir()
    G = 4; nc = G^2; d = 2 * nc; N = 80
    synth = () -> (theta = randn(7, N),
                   summary_min = vcat(randn(nc, N) .* 2 .+ 3, Float64.(rand(Bool, nc, N))))
    npe_kw = (; batchsize = 16, dstar = 8, depth = 1, width = 16,
              num_coupling_layers = 2, flow_depth = 1, flow_width = 8, stopping_epochs = 2)
    ratio_kw = (; batchsize = 16, num_summaries = 8, summary_width = 16, val_frac = 0.2,
                stopping_epochs = 2)

    b = train_and_register(G; datagen = synth, artifacts_root = root, use_gpu = false,
                           n_pairs = 60, ratio_n = 60, npe_epochs = 2, ratio_epochs = 2,
                           npe_kwargs = npe_kw, ratio_kwargs = ratio_kw)

    # A fully-populated EstimatorBundle for the grid, registered so estimator_for returns it.
    @test b isa ProteinCoLoc.EstimatorBundle
    @test b.grid == G
    @test estimator_for(G) === b
    @test b.zt !== nothing && b.θzt !== nothing
    @test haskey(b.ood_nulls, :density)
    @test isfinite(b.ratio.log_prior_odds)

    # The three grid artifacts exist under artifacts_root/grid_G/.
    gdir = joinpath(root, "grid_$(G)")
    @test isfile(joinpath(gdir, "npe_$(G).jld2"))
    @test isfile(joinpath(gdir, "ratio_$(G).jld2"))
    @test isfile(joinpath(gdir, "ood_nulls_$(G).jld2"))

    # SKIP-IF-DONE: a second call loads the persisted artifacts (a datagen that would ERROR is
    # never invoked ⇒ the skip path ran), returning a loadable bundle with a CPU-reproducible net.
    boom = () -> error("datagen must NOT run when artifacts are present (skip-if-done)")
    b2 = ProteinCoLoc._train_grid_pipeline(G; datagen = boom, artifacts_root = root,
                                           skip_if_done = true)
    @test b2 isa ProteinCoLoc.EstimatorBundle
    @test b2.grid == G
    Zq = Float32.(ProteinCoLoc.standardize_summary(
        vcat(randn(nc) .* 2 .+ 3, Float64.(rand(Bool, nc))), b2.zt, :min))
    @test size(ProteinCoLoc.posterior_for(b2.npe, Zq; N = 8, use_gpu = false), 1) == 7

    # --- TRAINING IMAGE-SIZE PROVENANCE (F6): persisted in meta, never guessed ----------------
    # This fixture INJECTED `datagen`, so the imsize_set keyword does NOT describe the pool: the
    # provenance must be recorded as :unknown rather than back-filled with default_imsize_for(G).
    @test isdefined(ProteinCoLoc, :training_imsize_provenance)
    for f in ("npe_$(G).jld2", "ratio_$(G).jld2", "ood_nulls_$(G).jld2")
        meta = JLD2.load(joinpath(gdir, f))["meta"]
        pv   = ProteinCoLoc.training_imsize_provenance(meta)
        @test hasproperty(meta, :imsize_set)          # the key EXISTS in every artifact
        @test pv.imsize_set == :unknown               # injected datagen ⇒ honestly unknown
        @test pv.imsize_source == :injected_datagen
        @test pv.recorded == false
    end

    # The real (non-injected) path records the ACTUAL set and weights it drew from. Exercising it
    # here would need the forward simulator + a full pool, so assert the accessor's contract on a
    # recorded meta of exactly the shape `_train_grid_pipeline` writes on the generate_samples path.
    rec = ProteinCoLoc.training_imsize_provenance(
              (; imsize_set = ((256, 256), (512, 512)), imsize_weights = (0.5, 0.5),
                 imsize_source = :generate_samples))
    @test rec.recorded == true
    @test rec.imsize_set == ((256, 256), (512, 512))
    @test rec.imsize_weights == (0.5, 0.5)

    # A PRE-PROVENANCE artifact (the frozen grid_4/8/16 bundles) reports :unknown — the accessor
    # must NEVER substitute today's default_imsize_for(G) for an unrecorded historical bundle.
    old = ProteinCoLoc.training_imsize_provenance((; grid = 16, n_pairs = 80_000, use_gpu = false))
    @test old.recorded == false
    @test old.imsize_set == :unknown
    @test ProteinCoLoc.training_imsize_provenance(nothing).recorded == false
end

##########################################################################################
### Per-grid CPU ship-gate harness + SBC + fresh pre-registration (07-04 Tasks 1-2, D-05)
###
### The per-grid ship-gate rides ONE grid-parametrized CPU (use_gpu=false) harness against a
### CPU-resident frozen net keyed by a FRESH disjoint PROD_SEED[G] (never the spike
### VAL_MASTER_SEED / training NPE_MASTER_SEED — Pitfall 3 / T-7-08). This smoke asserts the
### seed disjointness and drives the SBC machinery (rank table + KS/χ² + ECE) at fixture scale
### through the harness with an INJECTED fake simulator (the forward model is promoted by a later
### per-grid plan), plus that run_gate is invocable before any grid is trained.
##########################################################################################
include(joinpath(@__DIR__, "gate", "run_gate.jl"))   # brings gate_consts + harness + sbc + run_gate

@testset "ship-gate harness + SBC (gate/, 07-04)" begin
    # Fresh disjoint per-grid gate seeds (T-7-08): PROD_SEED[G] ∉ {VAL_MASTER_SEED, NPE_MASTER_SEED}.
    for G in (4, 8, 16, 32)
        @test PROD_SEED[G] != VAL_MASTER_SEED
        @test PROD_SEED[G] != NPE_MASTER_SEED
        @test PROD_SEED[G] != 0
    end
    @test length(unique(values(PROD_SEED))) == 4        # distinct per grid
    @test PROD_SALT != 0xBF58476D1CE4E5B9                # ≠ spike VAL_SALT
    @test prod_rng(8) isa Random123.Philox4x            # reproducible keyed stream

    # A tiny trained net + an INJECTED fake simulator (forward model not yet promoted) so the
    # harness + SBC machinery run end-to-end without the real simulator.
    G = 4; nc = G^2; d = 2 * nc; n = 40
    Zraw = vcat(randn(nc, n) .* 2 .+ 3, Float64.(rand(Bool, nc, n)))
    zt   = ProteinCoLoc.fit_summary_transform(Zraw; variant = :min)
    Zstd = Float32.(ProteinCoLoc.standardize_summary(Zraw, zt, :min))
    θ    = randn(7, n)
    res  = ProteinCoLoc.train_npe(Zstd, Zstd, θ, θ; use_gpu = false, epochs = 2, batchsize = 16,
                                  dstar = 8, depth = 1, width = 16, num_coupling_layers = 2,
                                  flow_depth = 1, flow_width = 8, stopping_epochs = 2)
    m = (estimator = res.estimator, zt = zt, θzt = res.θzt)

    fake_prior(rng) = (ρ_true = 2 * rand(rng) - 1, spillover = rand(rng),
                       autofluorescence = rand(rng), label_efficiency = rand(rng),
                       shift_dx = 0.0, shift_dy = 0.0, noise = rand(rng))
    fake_pair(rng, θ; imsize = (64, 64)) =
        [rand(rng, imsize...) .+ 0.5, rand(rng, imsize...) .+ 0.5]
    fake_mci(data) = MultiChannelImage(data, ["c1", "c2"], "gate_synth", ["p1", "p2"],
                                       size(data[1]), [0.5, 0.5])
    sim = (; sample_prior = fake_prior, simulate_pair = fake_pair, build_mci = fake_mci)

    # draw_simulate_infer is grid-parametrized and CPU-only; returns 7×N physical draws.
    t = draw_simulate_infer(m, prod_rng(G); G = G, imsize = (64, 64), N = SBC_FIX_L, sim = sim)
    @test size(t.draws, 1) == 7
    @test size(t.draws, 2) == SBC_FIX_L

    # tiny-M SBC rank table through the harness: M×8, every rank ∈ 0:L.
    ranks = sbc_ranks(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L, imsize = (64, 64),
                      sim = sim, rng = prod_rng(G))
    @test size(ranks) == (SBC_FIX_M, 8)
    @test all(r -> 0 <= r <= SBC_FIX_L, ranks)

    # KS + χ² uniformity via HypothesisTests (not hand-rolled) + ECE/MCE CalibrationResult.
    u = sbc_uniformity(ranks[:, 1]; L = SBC_FIX_L, bins = SBC_FIX_BINS)
    @test 0.0 <= u.ks_p <= 1.0
    @test isfinite(u.chi2_p)
    cal = sbc_calibration(collect(ranks[:, 1]); L = SBC_FIX_L, n_bins = SBC_FIX_BINS)
    @test cal isa CalibrationResult
    @test isfinite(cal.ece) && isfinite(cal.mce)
    @test sbc_traffic_light(cal.ece) in (:green, :yellow, :red)

    # The aggregated per-grid verdict runs and reports 8 columns + a pass/fail boolean.
    v = sbc_gate(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L, bins = SBC_FIX_BINS,
                 imsize = (64, 64), sim = sim, rng = prod_rng(G))
    @test length(v.per_param) == 8
    @test v.passed isa Bool
    @test v.caption == SBC_CAPTION

    # --- SHRINKAGE (vacuous-pass) DIAGNOSTIC: reporting only, verdicts untouched (F3) ---------
    # A rank-uniformity pass is NECESSARY but not SUFFICIENT: a posterior that merely reproduces
    # the prior gives uniform ranks BY CONSTRUCTION. Every column must therefore report
    # post_sd/prior_sd so such a column cannot be read as calibration evidence.
    @test all(x -> hasproperty(x, :shrinkage) && hasproperty(x, :post_sd) &&
                   hasproperty(x, :prior_sd) && hasproperty(x, :vacuous), v.per_param)
    @test all(x -> isfinite(x.post_sd) && x.post_sd >= 0, v.per_param)
    @test all(x -> x.vacuous isa Bool, v.per_param)
    @test hasproperty(v, :vacuous_params) && v.vacuous_cutoff == SBC_VACUOUS_SHRINKAGE
    # The flag is exactly the cutoff applied to the reported ratio (no hidden second rule).
    @test all(x -> !isfinite(x.shrinkage) ||
                   (x.vacuous == (x.shrinkage >= SBC_VACUOUS_SHRINKAGE)), v.per_param)

    # The diagnostic must NOT perturb the gate: `sbc_ranks_and_spread` consumes the passed rng in
    # the SAME order and amount as the old `sbc_ranks` body did (the extra outputs are std/
    # bookkeeping over samples that were ALREADY drawn), so the rank table — and hence every
    # verdict computed from it — is ==-identical.
    #
    # SEEDING NOTE: `posterior_for` → `NeuralEstimators.sampleposterior` takes no rng, so the flow's
    # posterior draws come from the GLOBAL RNG; the passed `prod_rng(G)` only drives the prior draw
    # and the forward simulation. Reproducing a rank table therefore requires pinning BOTH streams.
    # Since 07-10 the gate arms pin the global stream THEMSELVES (`seed_gate_global!(G)`); the
    # explicit seed here is retained only to keep this A/B comparison independent of that fix.
    _ranks_via(f) = (Random.seed!(20260721); f())
    base = _ranks_via(() -> sbc_ranks(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L,
                                      imsize = (64, 64), sim = sim, rng = prod_rng(G)))
    sr = _ranks_via(() -> sbc_ranks_and_spread(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L,
                                               imsize = (64, 64), sim = sim, rng = prod_rng(G)))
    @test sr.ranks == base                       # spread-returning chain ⇒ identical ranks
    @test size(sr.post_sd) == (SBC_FIX_M, 8)
    @test size(sr.prior_draws) == (SBC_FIX_M, 8)

    vA = _ranks_via(() -> sbc_gate(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L, bins = SBC_FIX_BINS,
                                   imsize = (64, 64), sim = sim, rng = prod_rng(G)))
    vB = _ranks_via(() -> sbc_gate(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L, bins = SBC_FIX_BINS,
                                   imsize = (64, 64), sim = sim, rng = prod_rng(G)))
    @test vA.ranks == base                       # the gate's ranks are the unchanged rank table
    @test vB.ranks == vA.ranks
    @test vB.passed == vA.passed && vB.ks_pass == vA.ks_pass && vB.ece_pass == vA.ece_pass
    @test [x.ks_p for x in vB.per_param] == [x.ks_p for x in vA.per_param]
    @test [x.ece for x in vB.per_param] == [x.ece for x in vA.per_param]
    # The verdict fields are a pure function of the ranks — the shrinkage columns feed nothing.
    @test vA.passed == (vA.ks_pass && vA.ece_pass)

    # A DEGENERATE posterior==prior column is flagged vacuous; an informative one is not.
    prior_col = randn(200, 1)
    vac = sbc_shrinkage(fill(StatsBase.std(prior_col), 200, 1), prior_col)
    @test vac.vacuous[1] == true                     # posterior spread == prior spread
    inf = sbc_shrinkage(fill(0.05 * StatsBase.std(prior_col), 200, 1), prior_col)
    @test inf.vacuous[1] == false                    # 20x tighter than the prior

    # run_gate is invocable before any grid is trained: honest :not_trained status, no crash.
    rep = run_gate(G; sbc = true, artifacts_root = mktempdir(), sim = sim, write_report = false)
    @test rep.status == :not_trained
    @test rep.grid == G

    # --- 07-10 REGRESSION: the gate arm is self-seeding and BIT-REPRODUCIBLE ------------------
    #
    # THE DEFECT THIS GUARDS: `posterior_for` → `NeuralEstimators.sampleposterior` accepts no
    # `rng`, so the flow's posterior draws come from the GLOBAL RNG. Passing `prod_rng(G)` pins
    # ONLY the prior draw + forward simulation. The spike pinned both streams
    # (`spike/demo.jl`: `Random.seed!(seed)` before `draw_simulate_infer`); the promotion into
    # `test/gate/` dropped the global half, making the recorded rank tables irreproducible.
    #
    # WHY THIS TEST IS STRONG: the global RNG is deliberately DESYNCHRONISED before each run (a
    # different `Random.seed!` plus a random-length burn), and NO external seeding is applied.
    # If `seed_gate_global!(G)` is ever removed from `sbc_ranks_and_spread` again, run A and run B
    # start from different global states and the rank tables diverge → this test FAILS.
    #
    # HISTORICAL CAVEAT (deliberate, do not "fix"): the grid-4/8/16 rank tables already recorded
    # in artifacts/ were produced under the UNSEEDED global stream and are NOT retroactively
    # reproducible. This test proves FUTURE runs are.
    _desync!(s) = (Random.seed!(s); rand(1 + (s % 97)); nothing)
    _gate_arm() = sbc_ranks_and_spread(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L,
                                       imsize = (64, 64), sim = sim, rng = prod_rng(G))

    _desync!(11)
    runA = _gate_arm()
    _desync!(9973)                       # global stream now in a totally different state
    runB = _gate_arm()
    @test runB.ranks == runA.ranks       # BIT-IDENTICAL rank table across independent runs
    @test runB.post_sd == runA.post_sd   # …and the posterior draws behind them are identical
    @test runB.prior_draws == runA.prior_draws

    # Same at the aggregated-verdict level: every reported p-value/ECE is bit-identical.
    _desync!(4242)
    gA = sbc_gate(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L, bins = SBC_FIX_BINS,
                  imsize = (64, 64), sim = sim, rng = prod_rng(G))
    _desync!(777)
    gB = sbc_gate(m; G = G, M = SBC_FIX_M, L = SBC_FIX_L, bins = SBC_FIX_BINS,
                  imsize = (64, 64), sim = sim, rng = prod_rng(G))
    @test gB.ranks == gA.ranks
    @test [x.ks_p for x in gB.per_param] == [x.ks_p for x in gA.per_param]
    @test [x.chi2_p for x in gB.per_param] == [x.chi2_p for x in gA.per_param]
    @test [x.ece for x in gB.per_param] == [x.ece for x in gA.per_param]
    @test gB.passed == gA.passed

    # CONTROL: without the arm's internal seeding the tables WOULD diverge — i.e. the assertions
    # above are not vacuously true because the estimator happens to be deterministic. Draw the
    # posterior directly (bypassing the arm, hence the seeding) from two different global states.
    Zprobe = ProteinCoLoc.standardize_summary(
                 ProteinCoLoc.encode_d01(ProteinCoLoc.patch_summary(
                     fake_mci(fake_pair(prod_rng(G), nothing)), G)), m.zt, :min)
    _desync!(11)
    dA = ProteinCoLoc.posterior_for(m.estimator, Zprobe; N = SBC_FIX_L, use_gpu = false)
    _desync!(9973)
    dB = ProteinCoLoc.posterior_for(m.estimator, Zprobe; N = SBC_FIX_L, use_gpu = false)
    @test dA != dB                       # posterior draws DO depend on the global RNG state

    # The global seed is derived from the frozen pre-registration, not a new constant.
    @test gate_global_seed(G) == gate_global_seed(G)          # deterministic
    @test gate_global_seed(4) != gate_global_seed(8)          # per-grid
    @test gate_global_seed(G) != prod_seed(G)                 # disjoint from the arm's own stream
    @test seed_gate_global!(G) == gate_global_seed(G)
end

##########################################################################################
### OOD control simulators for the ship-gate's --ood arm (07-09b, D-03/D-04)
###
### `test/gate/misspec.jl` supplies the POSITIVE controls (the four Phase-5 misspecification
### families) and the summary-orthogonal NEGATIVE controls the `--ood` arm needs in order to
### score a separability AUC at all. This smoke asserts the family surface, the imsize guard,
### the simulator-contract wrapper, and that the negative-control transforms really do leave the
### patch-Pearson summary (near-)invariant — the D-04 named blind spot, measured.
##########################################################################################
@testset "OOD control simulators (gate/misspec.jl, 07-09b)" begin
    G = 4; ims = (64, 64)
    @test length(OOD_FAMILIES) == 4
    @test collect(keys(OOD_FAMILIES)) == [:texture, :noise, :optics, :background]

    rng = prod_rng(G)
    θ   = ProteinCoLoc.sample_prior(rng)
    for (fam, gen) in pairs(OOD_FAMILIES)
        img = gen(rng, θ; imsize = ims, level = OOD_GRID_LEVELS, G = G)
        @test length(img) == 2
        @test all(ch -> size(ch) == ims, img)
        @test all(ch -> all(isfinite, ch) && all(>=(0.0), ch), img)
        # Grid-generalized imsize guard (mirrors simulate_pair's entry validation).
        @test_throws ArgumentError gen(rng, θ; imsize = (32, 32), level = 1, G = G)
    end

    # The simulator-contract wrapper drops straight into the harness seam.
    ps = misspec_simulator(OOD_FAMILIES.texture; level = 2, G = G)
    for f in (:sample_prior, :simulate_pair, :build_mci)
        @test haskey(ps, f)
    end
    @test length(ps.simulate_pair(rng, θ; imsize = ims)) == 2

    # Negative controls: correlation-preserving ⇒ (near-)KS-invariant summary (D-04 blind spot).
    base = ProteinCoLoc.simulate_pair(rng, θ; imsize = ims)
    negs = gate_negctrls(G; rng_factory = () -> prod_rng(G))
    for (name, tf) in pairs(negs)
        out = tf(base)
        @test length(out) == 2 && all(ch -> size(ch) == ims, out)
        δ = verify_summary_invariance(base, tf, G)
        @test isfinite(δ)
        # A G×G summary has only G² values, so the KS statistic is quantized in 1/G² steps;
        # a correlation-preserving transform must stay within a couple of those steps.
        @test δ <= 3 / G^2 + 1e-9
    end
end

##########################################################################################
### GPU-train smoke — graceful CPU fallback + CPU-resident persistence (07-04 Task 3, D-06)
##########################################################################################
include(joinpath(@__DIR__, "gpu_smoke.jl"))

##########################################################################################
### Windowed sub-tile local colocalization map (07-08, PROD-02)
###
### The in-process registration bridge (`_ensure_grid_registered`, wave 6 — Artifacts is not
### wired yet) plus `local_coloc_map`: r×c tile grid, per-tile Δρ finiteness, degenerate-tile
### sentinel, and the Phase-12 extension point left clean.
##########################################################################################
include(joinpath(@__DIR__, "test_local_map.jl"))

##########################################################################################
### RETIRED (v2.0 breaking release, D-01):
###
### The former public-API tests exercising the Turing/ADVI path — `colocalization()`,
### `compute_BayesFactor()`, `CoLocResult`, and the `plot_posterior`/`bayesplot`/
### `bayes_rangeplot` ADVI-result plotters — are removed from the core suite. That path is
### now the INTERNAL, weakdep-gated reference implementation in ext/ProteinCoLocTuringExt.jl
### and is validated by the per-grid D-05 ship-gate (later Phase-7 plans), which loads Turing
### to activate the extension. The GLMakie-backed plotting tests (`plot`,
### `local_correlation_plot`, `plot_mask`) are deferred pending the Wave-0 GLMakie/Makie
### co-resolution decision (see 07-00-SUMMARY.md § Blocker).
##########################################################################################
