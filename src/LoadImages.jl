## LoadImages.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder

module LoadImages
import CSV
import DataFrames
import Images
import PyCall
import Base: iterate, getindex

export MultiChannelImage, load_tiff, apply_mask!, MultiChannelImageStack, MultiChannelImage

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

struct MultiChannelImageStack{T <:MultiChannelImage, S <: AbstractString}
    img::Vector{T}
    name::S
    num_images::Int64 # number of images in the stack

    function MultiChannelImageStack(img::Vector{T}, name::S) where {T <: MultiChannelImage, S <: AbstractString}
        new{T, S}(img, name, length(img))
    end
end

# Method definition for struct MultiChannelImageStack
function Base.getindex(stack::MultiChannelImageStack, i::Int64)
    i > stack.num_images && throw(BoundsError(stack, i))
    return stack.img[i]
end

function Base.iterate(stack::MultiChannelImageStack, state=1)
    state > stack.num_images && return nothing
    return stack.img[state], state + 1
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
end
################################################################################
# convert lif file to tiff files

"""
    convert_lif_to_tiff(path::S) where {S<:AbstractString}

    Convert a lif file to tiff files. The lif file is expected to contain multiple
    images. Each image is saved in a separate folder. Each channel of the image is
    saved as a separate tiff file.

    Nota Bene: the current function cannot be used to convert Z-stacks in a meaningful way.
"""
function convert_lif_to_tiff(path::S) where {S<:AbstractString}
    # load python module for reading lif files
    readlif = PyCall.pyimport("readlif.reader")
    # read lif file
    f = readlif.LifFile(path)
    # get image meta_data
    meta_data = f.image_list
    # base path for saving images
    base_path = (dirname(path))

    # convert to tiff
    for img ∈ 1:f.num_images
        name = meta_data[img]["name"]
        image = f.get_image(img-1)
        image_path = "$base_path/$name"
        mkdir(image_path)
        for channel ∈ 1:meta_data[img]["channels"]
            image.get_frame(z=0, t=0, c=channel-1).save("$image_path/$name-channel-$channel.tiff")
        end
    end

    # convert meta data to csv file
    convert_lif_meta_data_to_csv(meta_data, base_path)
end

function convert_lif_meta_data_to_csv(meta_data::Vector{Dict{Any,Any}}, path::S) where {S<:AbstractString}
    df_meta_data = hcat([meta_data[1]["name"]], DataFrames.DataFrame(meta_data[1]["settings"]))
    for i ∈ 2:length(meta_data)
        try
            df_meta_data = vcat(
                df_meta_data, 
                hcat(
                    [meta_data[i]["name"]], 
                    DataFrames.DataFrame(meta_data[i]["settings"])
                    )
                )
        catch
            @warn "Could not convert meta data for image $i. This image will be skipped. Currently, z-stacks are not supported."
        end
    end

    # convert to csv
    CSV.write("$path/meta_data.csv", df_meta_data)
    
    return df_meta_data
end

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

end