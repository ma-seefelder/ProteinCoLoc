## colocalization.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder

module Colocalization
include("LoadImages.jl")
import .LoadImages
using Turing

path = ["test_images/c1.tif", "test_images/c2.tif", "test_images/c3.tif"]
img = LoadImages.MultiChannelImage("positive_sample", path, ["blue", "green", "red"])
mask = LoadImages._calculate_mask(img)
LoadImages._apply_mask!(img, mask)

######################################################################
# function to patch the image

"""
    patch(img::Array{Float64, 2}, num_patches::Int64)
    Patches the image into num_patches x num_patches patches.
    Return a 4D array with the patches.
"""
function patch(img::Array{Float64, 2}, num_patches::Int64)
    # calculate the size of the patches
    patch_size_x = Int64(floor(size(img, 1) / num_patches))
    patch_size_y = Int64(floor(size(img, 2) / num_patches))
    # calculate the number of patches
    num_patches_x = Int64(floor(size(img, 1) / patch_size_x))
    num_patches_y = Int64(floor(size(img, 2) / patch_size_y))
    # initialize the patches
    patches = zeros(Union{Float64, Missing}, num_patches_x, num_patches_y, patch_size_x, patch_size_y)
    # loop over the patches
    for i in 1:num_patches_x
        for j in 1:num_patches_y
            patches[i, j, :, :] = img[(i-1)*patch_size_x+1:i*patch_size_x, (j-1)*patch_size_y+1:j*patch_size_y]
        end
    end
    return patches
end

######################################################################
# function to calculate the correlation

"""
    _exclude_zero!(a::Vector{T},b::Vector{T}) where T <: Number
    Exclude all zero values from the vectors a and b.
"""
function _exclude_zero!(a::Vector{Union{T, Missing}},b::Vector{Union{T, Missing}}) where T <: Number
    # get the indices of the zero values
    a_zero = append!(findall(a .== 0), findall(isnan.(a)))
    b_zero = append!(findall(b .== 0), findall(isnan.(b)))
    # make union of the indices
    zero_indices = sort(union(a_zero, b_zero))
    # exclude all values at the indices in zero_indices
    deleteat!(a, zero_indices)
    deleteat!(b, zero_indices)
end

"""
    correlation(x::Array{Float64, 4}, y::Array{Float64, 4})
    Calculate the correlation between two 4D arrays.
    Return a 2D array with the correlation for each patch.
"""
function correlation(x::Array{T, 4}, y::Array{T, 4}) where T <: Union{Float64, Missing}
    # get the number of patches
    num_patches = size(x, 1)
    # initialize the correlation
    ρ = zeros(Union{Float64, Missing},num_patches, num_patches)
    # loop over the patches
    for i in 1:num_patches
        for j in 1:num_patches
            a = x[i, j, :, :][:] 
            b = y[i, j, :, :][:]
            _exclude_zero!(a,b)
            length(a) <= 3 ? ρ[i, j] = missing : ρ[i, j] = cor(a, b)
        end
    end
    return ρ
end

struct CoLocalizeResult
    img::LoadImages.MultiChannelImage
    control::LoadImages.MultiChannelImage
    channels::Vector{Int64}
    num_patches::Int64
    num_chains::Int64
    ρ_control::Array{Float64, 2}
    ρ_sample::Array{Float64, 2}
    ρ_diff::Array{Float64, 2}
    ρ_control_post::Array{Float64, 2}
    ρ_sample_post::Array{Float64, 2}
    ρ_diff_post::Array{Float64, 2}
end


function colocalization_Model(
    img::LoadImages.MultiChannelImage, 
    control::LoadImages.MultiChannelImage, 
    channels::Vector{Int64},
    num_patches::Int64 = 1;
    num_chains::Int64 = num_chains
    )
    # calculate colocalization
    x = img.data[channels[1]]
    y = img.data[channels[2]]
    x_control = control.data[channels[1]]
    y_control = control.data[channels[2]]

    # assert that the images are the same size
    @assert size(x) == size(y) == size(x_control) == size(y_control) "Images are not the same size"

    # patch the image
    x, y, x_control, y_control = patch.([x, y, x_control, y_control], num_patches)

    # calculate the correlation for each patch
    ρ_sample = correlation(x, y)
    ρ_control = correlation(x_control, y_control)

    # remove missing values
    ρ_sample = ρ_sample[.!ismissing.(ρ_sample)]
    ρ_control = ρ_control[.!ismissing.(ρ_control)]


    # Fisher's z-transformation of the correlation coefficients
    ρ_sample = atanh.(ρ_sample)
    ρ_control = atanh.(ρ_control)

    ###############################
    # Turing model
    ###############################

    @model function model_singular_image(control = ρ_control, sample = ρ_sample)
        # get the number of patches
        num_patches_x = size(control, 1)
        num_patches_y = size(sample, 1)

        # hyperpriors
        μ_control ~ Normal(0, 1)
        ν_control ~ Gamma(2, 0.1)
        σ_control ~ InverseGamma(ν_control, 0.1)

        μ_sample ~ Normal(0, 1)
        ν_sample ~ Gamma(2, 0.1)
        σ_sample ~ InverseGamma(ν_sample, 0.1)

        # priors
        # get prior for each patch
        μ_control_patch ~ filldist(Normal(μ_control, σ_control), num_patches)
        ν_control_patch ~ filldist(InverseGamma(ν_control, 0.1), num_patches)

        μ_sample_patch ~ filldist(Normal(μ_sample, σ_sample), num_patches)
        ν_sample_patch ~ filldist(InverseGamma(ν_sample, 0.1), num_patches)

        # likelihood
        for p ∈ 1:num_patches_x
            control[p] ~ TDist(ν_control_patch[p]) + μ_control_patch[p]
        end

        for p ∈ 1:num_patches_y
            sample[p] ~ TDist(ν_sample_patch[p]) + μ_sample_patch[p]
        end           
    end

    # sample
    chains = sample(model_singular_image(ρ_control,ρ_control), NUTS(), MCMCThreads(), 1000, num_chains; discard_adapt=false)

    # get the samples
    ρ_control_post = tanh.(chains[:μ_control])
    ρ_sample_post = tanh.(chains[:μ_sample])

    # calculate the difference
    ρ_diff = ρ_sample_post - ρ_control_post

    return CoLocalizeResult(
        img, 
        control, 
        channels, 
        num_patches, 
        num_chains, 
        ρ_control, 
        ρ_sample, 
        ρ_diff,
        ρ_control_post, 
        ρ_sample_post, 
        ρ_diff_post
        )

end


function colocalization(img::LoadImages.MultiChannelImage, channels::Vector{Int64})
    mask = LoadImages._calculate_mask(img)
    LoadImages._apply_mask!(img, mask)

    # calculate colocalization
    channel_1 = img.data[channels[1]]
    channel_2 = img.data[channels[2]]

    # calculate colocalization



    return img
end

end