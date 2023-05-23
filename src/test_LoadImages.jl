using Test
include("../src/LoadImages.jl")
using .LoadImages

@testset "LoadImages" verbose = true begin 
    # define a test set for the load_tiff function
    @testset "load_tiff" begin
        # define test cases
        test_cases = [
            ("test_images/c1.tif", 1028, 1376),
            ("test_images/c2.tif", 1028, 1376),
            ("test_images/c3.tif", 1028, 1376)
        ]

        # loop through test cases
        for (path, expected_rows, expected_cols) in test_cases
            # define test name
            test_name = "loading $path"

            # call load_tiff function
            img = load_tiff(path)

            # check that the output is a matrix
            @test typeof(img) == Matrix{Float64}

            # check that the output has the expected dimensions
            @test size(img, 1) == expected_rows
            @test size(img, 2) == expected_cols
        end
    end

    # define a test set for the MultiChannelImage constructor
    @testset "MultiChannelImage constructor" begin
        path = ["test_images/c1.tif", "test_images/c2.tif", "test_images/c3.tif"]
        name = "test_image"
        channels = ["blue", "green", "red"]

        # test that the constructor returns a MultiChannelImage object
        img = MultiChannelImage(name, path, channels)
        @test typeof(img) == MultiChannelImage{Float64, String, Float64}

        # test that the constructor sets the fields correctly
        @test img.channels == channels
        @test img.name == name
        @test img.path == path
        @test img.pixel_size == (1028, 1376)
        @test length(img.otsu_threshold) == 3
    end

    @testset "Mask" begin
        path = ["test_images/c1.tif", "test_images/c2.tif", "test_images/c3.tif"]
        name = "test_image"
        channels = ["blue", "green", "red"]
        # test that the constructor returns a MultiChannelImage object
        img = MultiChannelImage(name, path, channels)
        img = apply_mask!(img)
    end
end