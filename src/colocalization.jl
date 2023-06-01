## colocalization.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder

#module Colocalization
include("LoadImages.jl")
import .LoadImages
import DataFrames: DataFrame
using Turing
using Turing: Variational

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
######################################################################
# function to convert the posterior samples
######################################################################

function convert_posterior_samples(samples::Array{Float64, 2}, n_control::Int64, n_sample::Int64)
    # check input arguments
    6 + 2*n_control + 2*n_sample == size(samples, 1) || error("The number of parameters does not match the number of samples.")

    # get parameter_names
    parameter_names = [:μ_control, :ν_control, :σ_control, :μ_sample, :ν_sample, :σ_sample]
    [push!(parameter_names, Symbol("μ_control_", i)) for i in 1:n_control]
    [push!(parameter_names, Symbol("σ_control_", i)) for i in 1:n_control]
    [push!(parameter_names, Symbol("μ_sample_", i)) for i in 1:n_sample]
    [push!(parameter_names, Symbol("σ_sample_", i)) for i in 1:n_sample]
        
    # permute the samples
    samples = DataFrame(permutedims(samples, [2, 1]), parameter_names)

    # undo Fisher z transformation
    #samples = tanh.(samples)
end

######################################################################
# plot posterior
######################################################################
struct CoLocResult
    img::LoadImages.MultiChannelImage
    control::LoadImages.MultiChannelImage
    channels::Vector{Int64}
    num_patches::Int64
    posterior::DataFrame
    advi_result
end

function plot_posterior(posterior::CoLocResult)
    Plots.histogram(posterior.posterior.μ_control, legend = true, label = "μ_control")
    Plots.histogram!(posterior.posterior.μ_sample, label = "μ_sample")
    Plots.histogram(posterior.posterior.μ_control_1, label = "μ_control_1")

    Plots.histogram(posterior.posterior.ν_control, label = "ν_control", legend = true)
    Plots.histogram!(posterior.posterior.ν_sample, label = "ν_sample")

    Plots.histogram(posterior.posterior.σ_control, label = "σ_control", legend = true)
    Plots.histogram!(posterior.posterior.σ_sample, label = "σ_sample")

    sum(ρ_diff .> 0) / length(ρ_diff) / sum(rand(Normal(0, 1), 10_000) .> 0) / 10_000
end


function colocalization(
    img::LoadImages.MultiChannelImageStack, 
    control::LoadImages.MultiChannelImageStack, 
    channels::Vector{Int64},
    num_patches::Int64 = 1;
    num_chains::Int64 = num_chains
    )

    sample = fill(0.0, length(img), num_patches)
    for (image,idx) ∈ zip(img, 1:length(img))
        x = image.data[channels[1]]
        y = image.data[channels[2]]
        x,y = patch.([x, y], num_patches)
        sample[idx,:] = correlation(x, y)
    end

    ctrl = fill(0.0, length(img), num_patches)
    for (image,idx) ∈ zip(control, 1:length(control))
        x = image.data[channels[1]]
        y = image.data[channels[2]]
        x,y = patch.([x, y], num_patches)
        ctrl[idx,:] = correlation(x, y)
    end

    ################################
    # Turing model
    ################################
    # input data for the model are the correlation coefficients of the sample and control for each image 
    # type: Array{Float64, 2} with size (num_images, num_patches)

    @model function model(control::Array{Float64,2} = ctrl, sample::Array{Float64,2} = sample)
        # get the number of patches
        num_control = size(control, 1)
        num_sample = size(sample, 1)

        # ============= gloabal priors (per biological condition) ============= #
        # mean, degrees of freedom and standard deviation of the control 
        μ_control ~ Truncated(Normal(0, 1),-1,1)
        ν_control ~ Exponential()
        σ_control ~ Truncated(Cauchy(std(control), 0.5 * std(control)),0,1)
        # mean, degrees of freedom and standard deviation of the sample
        μ_sample ~ Truncated(Normal(0, 1),-1,1)
        ν_sample ~ Exponential()
        σ_sample ~ Truncated(Cauchy(std(sample), 0.5 * std(sample)),0,1)

        # ============= local priors (per image) ============= #
        μ_control_image ~ filldist(Truncated(Normal(μ_control, σ_control),-1,1), num_control) # mean of the control for each patch
        ν_control_image ~ filldist(Exponential(ν_control), num_control) # degress of freedom of the control for each patch

        μ_sample_image ~ filldist(Truncated(Normal(μ_sample, σ_sample),-1,1), num_sample) # mean of the sample for each patch
        ν_sample_image ~ filldist(Exponential(ν_control), num_sample) # degress of freedom of the sample for each patch

        # likelihood
        for idx ∈ 1:num_control
            control[idx] ~ TDist(ν_control_image[idx]) + μ_control_image[idx]
        end

        for idx ∈ 1:num_sample
            sample[idx] ~ TDist(ν_sample_image[idx]) + μ_sample_image[idx]
        end
    end

    # sample
    m = model_singular_image(x, y)
    q = vi(m, ADVI(10, 1000))

    # get the posterior samples
    q_samples = rand(q, 10_000)
    q_samples = convert_posterior_samples(q_samples, length(ρ_control), length(ρ_sample))
    
    # get the samples

    # calculate the difference
    ρ_diff = q_samples.μ_sample - q_samples.μ_control

    return CoLocResult(
        img, 
        control, 
        channels, 
        num_patches, 
        q_samples,
        q
        )

end

