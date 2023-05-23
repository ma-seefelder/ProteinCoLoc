## LoadImages.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder

module LoadImages
import Images

export MultiChannelImage, load_tiff, apply_mask!

################################################################################
# Data types
struct MultiChannelImage{T <: Union{Missing, Float64}, S <: AbstractString, F <: AbstractFloat}
    data::Vector{Matrix{T}}
    channels::Vector{S}
    name::S
    path::Vector{S}
    pixel_size::Tuple{F, F}
    otsu_threshold::Vector{F}

    function MultiChannelImage(data::Vector{Matrix{T}}, channels::Vector{S}, name::S, path::Vector{S}, pixel_size::Tuple{Int64, Int64}, otsu_threshold::Vector{F}) where {T <: Union{Missing, Float64}, S <: AbstractString, F <: AbstractFloat}
        @assert length(data) == length(channels) "Length of data and channels vectors must match"
        @assert length(otsu_threshold) == length(channels) "Length of otsu_threshold vector must match number of channels"
        new{T, S, F}(data, channels, name, path, pixel_size, otsu_threshold)
    end
end

################################################################################
# Load images from tiff files

"""
    load_tiff(path::S) where {S<:AbstractString}

    Load a tiff image from a file and convert it to a matrix of Float64 values where 
    the values represent the intensity on grayscale.
"""
function load_tiff(path::S) where {S<:AbstractString}
    img = Images.load(path) # load image
    img = Images.Gray.(img) # convert to grayscale
    img = Images.convert(Matrix{Float64}, img) # convert to matrix
    return img
end

################################################################################
# Calculate mask
function _calculate_mask(img)
    mask = fill(ones(size(img.data[1])), length(img.channels))
    for i ∈ 1:length(img.channels)
        mask[i] = img.data[i] .> img.otsu_threshold[i]
    end
    return mask
end


# Apply mask
function _apply_mask!(img, mask)
    for channel ∈ 1:length(img.channels)
        img.data[channel] = img.data[channel] .* mask[channel]
    end
    return img
end

function apply_mask!(img)
    mask = _calculate_mask(img)
    img = _apply_mask!(img, mask)
    return img
end
################################################################################


################################################################################
# constructor for MultiChannelImage
function MultiChannelImage(name::S, path::Vector{S}, channels::Vector{S} =[]) where {S <: AbstractString}
    data = load_tiff.(path)

    if isempty(channels)
        @warn "No channel names provided, using default names"
        channels = ["channel_$i" for i in 1:length(data)]
    end

    return MultiChannelImage(
        data, 
        channels, 
        name, 
        path, 
        size(data[1]), 
        Images.otsu_threshold.(data)
        )
end

end # module LoadImages