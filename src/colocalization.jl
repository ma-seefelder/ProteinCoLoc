## colocalization.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder
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
# function to calculate the Pearson's correlation
######################################################################

"""
    _exclude_zero!(a::Vector{T},b::Vector{T}) where T <: Number
    Exclude all values that are zero in at least one of the vectors a and b.
    Return the vectors a and b without the zero and missing values as a Vector{Float64}.
"""
function _exclude_zero!(a::Vector{T},b::Vector{T}) where T <: Number  
    # get the indices of the zero values
    a_zero = append!(findall(a .== 0), findall(isnan.(a)))
    b_zero = append!(findall(b .== 0), findall(isnan.(b)))
    # make union of the indices
    zero_indices = sort(union(a_zero, b_zero))
    # exclude all values at the indices in zero_indices
    deleteat!(a, zero_indices)
    deleteat!(b, zero_indices)
end

function _exclude_zero!(a::Vector{Union{T,Missing}}, b::Vector{Union{T,Missing}}) where T <: Number
    # change missing values to zero
    a[ismissing.(a)] .= 0
    b[ismissing.(b)] .= 0

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
# function to calculate the Fractional Overlap
######################################################################
# implemented after https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3074624/

"""
    struct Fractional Overlap()
    Structure to save the results of the fractional overlap analysis. 
    The function fractional_overlap() returns this  structure that
    possesses the following fields:
        x: Fractional overlap of channel_1 / channel_2
        y: Fractional overlap of channel_2 / channel_1
        x_above_threshold: proportion of pixels above the threshold for channel 1
        y_above_threshold: proportion of pixels above the threshold for channel 2
        channel_1: name of channel 1
        channel_2: name of channel 2
"""
struct FractionalOverlap
    x::T where {T<:Real}
    y::T where {T<:Real}
    x_above_threshold::T where {T<:Real}
    y_above_threshold::T where {T<:Real}
    channel_1::T where {T<:Int}
    channel_2::T where {T<:Int}
end

"""
    print(IO::Core.IO, x::FractionalOverlap)
    Print statement for the FractionalOverlap structure. 
"""
function print(IO::Core.IO, x::FractionalOverlap)
    expression = "$(round(100*x.x; digits = 3))% of the $(x.channel_1) channel overlap with the $(x.channel_2) channel and $(round(100*x.y; digits = 3))% of the $(x.channel_2) channel overlap with the $(x.channel_1) channel." 
    return println(IO, expression)
end

"""
    _clip_to_zero(x::Matrix{T}) where T<: Real
    Helper function to clip non-zero values to zero.
"""
function _clip_to_zero(x::Matrix{T}) where T <: Union{Float64, Missing}
    Threads.@threads for i in CartesianIndices(x)
        if x[i] <= 0.0
            x[i] = 0
        end
    end
    return x    
end

"""
    prop_above_threshold(x::Matrix{T}, threshold::T) where T <: Real
    Calculates the proportion of pixels above the defined threshold. 

"""
function prop_above_threshold(x::Matrix{T}, threshold::T) where T <: Union{Float64, Missing}
    return sum(_clip_to_zero(x .- threshold)) / (size(x)[1] * size(x)[2])
end

"""
mcc(
x::Matrix{Union{Float64, Missing}}, 
y::Matrix{Union{Float64, Missing}}, 
threshold_x::Float64, 
threshold_y::Float64
)

mcc(x, y, threshold_x::Float64, threshold_y::Float64)
Low level function to compute Manders' Colocalization Coefficients (MCC) based
on the matrices x and y and the thresholds for both values. 
"""
function mcc(
    x::Matrix{T},
    y::Matrix{T},
    threshold_x::T, 
    threshold_y::T
    ) where T <: Union{Float64, Missing}

    if ismissing(threshold_x)
        threshold_x = 0.0
    end

    if ismissing(threshold_y)
        threshold_y = 0.0
    end

    x = _clip_to_zero(x .- threshold_x)
    y = _clip_to_zero(y .- threshold_y)
    x_colocal = deepcopy(x)

    # return missing if sum(x) == 0.0
    if sum(x) == 0.0
        return missing
    end

    # calculate colocalization
    Threads.@threads for i in CartesianIndices(x_colocal)
        if y[i] <= 0.0
            x_colocal[i] = 0
        end
    end  
    return sum(x_colocal) / sum(x)
end

"""
    fractional_overlap(
    img::MultiChannelImage,
    channels::Vector{I}
    ) where {I <: Int}


    fractional_overlap(
    img::MultiChannelImage,
    channels::Vector{I};
    use_otsu::Bool = true,
    ) where {I <: Int}

    This function calculates the fractional overlap of two probes and
    requires the following input arguments:
    img::MultiChannelImage: a MultiChannelImage object
    channels: Vector of channel indices to be used for the calculation
"""
function fractional_overlap(
    img::MultiChannelImage,
    channels::Vector{I}
    ) where {I <: Int}

    # calculate and apply mask
    img = _apply_mask!(img, _calculate_mask(img))
    
    # calculate MCC values
    M_1 = mcc(img.data[channels[1]], img.data[channels[2]], missing, missing)
    M_2 = mcc(img.data[channels[2]], img.data[channels[1]], missing, missing)

    fr_1 = prop_above_threshold(img.data[channels[1]], 0) # fraction of pixels in channel 1 above threshold
    fr_2 = prop_above_threshold(img.data[channels[2]], 0) # fraction of pixels in channel 2 above threshold
    return FractionalOverlap(M_1, M_2, fr_1, fr_2, channels[1], channels[2])
end

function fractional_overlap(
    x::Matrix{T},
    y::Matrix{T},
    threshold_x::T,
    threshold_y::T
    ) where {T <: Union{Float64, Missing}}

    # calculate MCC values
    M_1 = mcc(x, y, threshold_x, threshold_y)
    M_2 = mcc(y, x, threshold_x, threshold_y)

    return M_1, M_2
end

"""
    fractional_overlap(
    img::MultiChannelImage,
    control::MultiChannelImage,
    channels::Vector{I},
    num_patches::I;
    method::String = "quantile",
    quantile_level::Float64 = 0.95
    ) where {I <: Int}

    This function calculates the fractional overlap of two probes and
    requires the following input arguments:
    img::MultiChannelImage: a MultiChannelImage object
    control::MultiChannelImage: a MultiChannelImage object to be used as control
    channels: Vector of channel indices to be used for the calculation
    num_patches: number of patches to be used for the calculation
    method: method to calculate the threshold (either "quantile","max", "median")
    quantile_level: quantile level to be used for the calculation of the threshold
"""
function fractional_overlap(
    img::MultiChannelImage,
    control::MultiChannelImage,
    channels::Vector{I},
    num_patches::I;
    method::String = "quantile",
    quantile_level::Float64 = 0.95
    ) where {I <: Int}

    # calculate and apply mask
    img = _apply_mask!(img, _calculate_mask(img))
    control = _apply_mask!(control, _calculate_mask(control))

    # calculate the thresholds based on the control
    if method == "quantile"
        threshold_x = quantile(control.data[channels[1]][:], quantile_level)
        threshold_y = quantile(control.data[channels[2]][:], quantile_level)
    elseif method == "max"
        threshold_x = maximum(control.data[channels[1]][:])
        threshold_y = maximum(control.data[channels[2]][:])
    elseif method == "median"
        threshold_x = median(control.data[channels[1]][:])
        threshold_y = median(control.data[channels[2]][:])
    else
        error("Method not implemented.")
    end

    # patch the image
    x = patch(img.data[channels[1]], num_patches)
    y = patch(img.data[channels[1]], num_patches)

    # preallocate the matrices
    M_1 = zeros(Union{Float64, Missing}, num_patches, num_patches)
    M_2 = zeros(Union{Float64, Missing}, num_patches, num_patches)

    # calculate MCC values
    for i in 1:num_patches
        for j in 1:num_patches
            M_1[i, j], M_2[i, j] = fractional_overlap(x[i, j, :, :], y[i, j, :, :], threshold_x, threshold_y)
        end
    end
    return M_1, M_2
end