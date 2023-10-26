#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
using Test

##########################################################################################
### Test for Image Loading functionality
##########################################################################################
#include("./src/LoadImages.jl")
#include("./src/colocalization.jl")
using .LoadImages
using .Colocalization
import Images

@testset "LoadImages" verbose = true begin 
    path = ["test/test_images/c1.tif", "test/test_images/c2.tif", "test/test_images/c3.tif"]
    # define a test function for load_tiff
    @testset "load_tiff" begin
        img = LoadImages.load_tiff(path[1])
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
        @test typeof(img) == MultiChannelImage{Float64, String, Float64}

        # test that the constructor sets the fields correctly
        @test img.channels == channels
        @test img.name == name
        @test img.path == path
        @test img.pixel_size == (1028, 1376)
        @test length(img.otsu_threshold) == 3
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

        @test LoadImages._calculate_mask(img) == fill(fill(0, 128, 128),3)
        img = LoadImages._apply_mask!(img, LoadImages._calculate_mask(img))
        @test img.data == fill(fill(0, 128, 128),3)

        # test with image where all pixels are above the threshold
        # value of mask should be 1 everywhere
        img = MultiChannelImage(
            fill(fill(1.0, 128, 128),3),
            ["blue", "green", "red"],
            "all_below_otsu",
            path,(128, 128),[0.5, 0.5, 0.5]
            )

        @test LoadImages._calculate_mask(img) == fill(fill(1, 128, 128),3)
        img = LoadImages._apply_mask!(img, LoadImages._calculate_mask(img))
        @test img.data == fill(fill(1, 128, 128),3)
        # confirm that the image size is not changed by applying the mask
        @test size(img.data[1]) == (128, 128) 
        @test size(img.data[2]) == (128, 128)
        @test size(img.data[3]) == (128, 128) 
    end

    @testset "MultiChannelImageStack" verbose = true begin
        # load image
        img = LoadImages.MultiChannelImage("test_image", path, ["blue", "green", "red"])
        # create stack from image and test that the constructor returns a MultiChannelImageStack object
        img_stack = LoadImages.MultiChannelImageStack([img, img, img], "test_stack")
        @test typeof(img_stack) == LoadImages.MultiChannelImageStack{LoadImages.MultiChannelImage{Float64, String, Float64}, String}

        # test that the constructor sets the fields correctly
        @test img_stack.name == "test_stack"
        @test img_stack.num_images == 3

        # test that indexing works as expected
        @test img_stack[1] == img
        @test_throws BoundsError img_stack[4] # check that indexing out of bounds throws an error  

        # test that iterations over the stack work as expected
        for (image, idx, i) ∈ zip(img_stack, 1:img_stack.num_images, 1:3)
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
    patches = Colocalization.patch(img, num_patches)

    # check that the output is a 4D array
    @test typeof(patches) == Array{Union{Float64, Missing}, 4}

    # check that the output has the expected dimensions
    @test size(patches) == (num_patches, num_patches, size(img, 1) ÷ num_patches, size(img, 2) ÷ num_patches)

    # check that the values in the patches match the corresponding values in the original image
    for i in 1:num_patches
        for j in 1:num_patches
            @test patches[i, j, :, :] == img[(i-1)*size(img, 1) ÷ num_patches+1:i*size(img, 1) ÷ num_patches, (j-1)*size(img, 2) ÷ num_patches+1:j*size(img, 2) ÷ num_patches]
        end
    end
end

###########################################################################################
### Test for correlation calculation
###########################################################################################
@testset "Colocalization" verbose = true begin 
    # define a test set for the _exclude_zero! function
    @testset "_exclude_zero!" begin
        # define a test case
        a = [1, 2, 0, 3, 4, 0, 5]
        b = [0, 2, 3, 0, 4, 5, 6]

        # call _exclude_zero! function
        Colocalization._exclude_zero!(a, b)

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
        ρ = Colocalization.correlation(x, y)

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