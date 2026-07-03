#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
import Images
import Statistics: cor
import Pkg

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
