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

################################################################################
# Data types
################################################################################
"""
    mutable struct MultiChannelImage{T <: Union{Missing, Float64}, S <: AbstractString, F <: AbstractFloat}
        data::Vector{Matrix{T}}
        channels::Vector{S}
        name::S
        path::Vector{S}
        pixel_size::Tuple{F, F}
        otsu_threshold::Vector{F}

    This mutable struct represents a multi-channel image.

    # Fields
    - `data`: A Vector of Matrix{T} representing the image data for each channel.
    - `channels`: A Vector of strings representing the names of the channels.
    - `name`: A string representing the name of the image.
    - `path`: A Vector of strings representing the paths to the image files for each channel.
    - `pixel_size`: A Tuple of two Float64 representing the width and height of the pixels in micrometers.
    - `otsu_threshold`: A Vector of Float64 representing the Otsu threshold for each channel.

    # Constructor
    The constructor checks that the lengths of the `data`, `channels`, and `otsu_threshold` vectors match.

    # Notes
    This struct is used to store and manipulate multi-channel images in the ProteinCoLoc package.
    """
mutable struct MultiChannelImage{T <: Union{Missing, Float64}, S <: AbstractString, F <: AbstractFloat}
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

"""
    mutable struct MultiChannelImageStack{T <: MultiChannelImage, S <: AbstractString}
        img::Vector{T}
        name::S
        num_images::Int64

    This mutable struct represents a stack of multi-channel images.

    # Fields
    - `img`: A Vector of MultiChannelImage representing the images in the stack.
    - `name`: A string representing the name of the image stack.
    - `num_images`: An Int64 representing the number of images in the stack.

    # Constructor
    The constructor takes a Vector of MultiChannelImage and a string as arguments, and sets the `num_images` field to the length of the Vector.

    # Notes
    This struct is used to store and manipulate stacks of multi-channel images in the ProteinCoLoc package.
    """
mutable struct MultiChannelImageStack{T <:MultiChannelImage, S <: AbstractString}
    img::Vector{T}
    name::S
    num_images::Int64 # number of images in the stack

    function MultiChannelImageStack(img::Vector{T}, name::S) where {T <: MultiChannelImage, S <: AbstractString}
        new{T, S}(img, name, length(img))
    end
end

# Method definition for struct MultiChannelImageStack
"""
    Base.getindex(stack::MultiChannelImageStack, i::Int64)

This function overrides the base getindex function to retrieve an image from a MultiChannelImageStack.

# Arguments
- `stack`: A MultiChannelImageStack from which to retrieve an image.
- `i`: An Int64 representing the index of the image to retrieve.

# Returns
- `img`: A MultiChannelImage representing the image at index `i` in the stack.

# Errors
- Throws a BoundsError if `i` is greater than the number of images in the stack.

# Notes
This function checks that `i` is within the bounds of the stack, and returns the image at index `i`.
"""
function Base.getindex(stack::MultiChannelImageStack, i::Int64)
    i > stack.num_images && throw(BoundsError(stack, i))
    return stack.img[i]
end

"""
    Base.iterate(stack::MultiChannelImageStack, state=1)

This function overrides the base iterate function to iterate over the images in a MultiChannelImageStack.

# Arguments
- `stack`: A MultiChannelImageStack over which to iterate.
- `state`: An Int64 representing the current state of the iteration. Default is 1.

# Returns
- A tuple containing the current image and the next state, or nothing if the end of the stack has been reached.

# Notes
This function checks if the state is greater than the number of images in the stack, and if so, returns nothing. Otherwise, it returns the current image and the next state.
"""
function Base.iterate(stack::MultiChannelImageStack, state=1)
    state > stack.num_images && return nothing
    return stack.img[state], state + 1
end
################################################################################
# Load images from tiff files

"""
    load_tiff(path::S) where {S<:AbstractString}

This function loads a TIFF image from a file and converts it to a matrix of Float64 values, where the values represent the intensity in grayscale.

# Arguments
- `path`: A string representing the path to the TIFF image file.

# Returns
- `img`: A Matrix of Float64 representing the grayscale intensity of the image.

# Notes
This function uses the Images.jl package to load the image, convert it to grayscale, and then convert it to a matrix.
"""
function load_tiff(path::S) where {S<:AbstractString}
    img = Images.load(path) # load image
    img = Images.Gray.(img) # convert to grayscale
    img = Images.convert(Matrix{Float64}, img) # convert to matrix
    return img
end

################################################################################
# Calculate mask
"""
    _calculate_mask(img::MultiChannelImage)

This function calculates a binary mask for each channel of a MultiChannelImage, where pixels with intensity greater than the Otsu threshold are set to 1 and others are set to 0.

# Arguments
- `img`: A MultiChannelImage for which to calculate the mask.

# Returns
- `mask`: A Vector of Matrix{Bool} representing the binary mask for each channel.

# Notes
This function iterates over the channels of the image, applies the Otsu threshold to the data of each channel, and stores the resulting binary mask.
"""
function _calculate_mask(img)
    mask = fill(ones(size(img.data[1])), length(img.channels))
    for i ∈ 1:length(img.channels)
        mask[i] = img.data[i] .> img.otsu_threshold[i]
    end
    return mask
end


# Apply mask
"""
    _apply_mask!(img::MultiChannelImage, mask::Vector{Matrix{Bool}})

This function applies a binary mask to each channel of a MultiChannelImage, where pixels in the mask that are set to 0 are set to 0 in the image.

# Arguments
- `img`: A MultiChannelImage to which to apply the mask.
- `mask`: A Vector of Matrix{Bool} representing the binary mask for each channel.

# Returns
- `img`: The MultiChannelImage with the mask applied.

# Notes
This function iterates over the channels of the image, multiplies the data of each channel by the corresponding mask, and stores the result back in the image.

---

    apply_mask!(img::MultiChannelImage)

This function calculates a binary mask for a MultiChannelImage and applies it to the image.

# Arguments
- `img`: A MultiChannelImage to which to apply the mask.

# Returns
- `img`: The MultiChannelImage with the mask applied.

# Notes
This function calls the `_calculate_mask` function to calculate the mask, and then calls the `_apply_mask!` function to apply the mask to the image.
"""
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
# shuffle individual pixels in each channel of a MultiChannelImage
"""
    shuffle_pixels(img::MultiChannelImage)

This function shuffles the pixel values in each channel of a MultiChannelImage.

# Arguments
- `img`: A MultiChannelImage whose pixel values are to be shuffled.

# Returns
- `img`: The MultiChannelImage with shuffled pixel values.

# Notes
This function iterates over the channels of the image, shuffles the pixel values in each channel using the `shuffle!` function from the Julia standard library, and stores the result back in the image.
"""
function shuffle_pixels(img::MultiChannelImage)
    for channel ∈ 1:length(img.channels)
        # retrieve pixel values
        img.data[channel] = shuffle!(img.data[channel])
    end
    return img
end

# shuffle blocks of pixels in each channel of a MultiChannelImage 
# with block size block_size to accomodate for autocorrelation
# use the patch function to create blocks
"""
    shuffle_blocks(img::MultiChannelImage, block_size::Int64)

This function shuffles the blocks of pixel values in each channel of a MultiChannelImage.

# Arguments
- `img`: A MultiChannelImage whose pixel blocks are to be shuffled.
- `block_size`: An Int64 representing the size of the blocks to be shuffled.

# Returns
- `img`: The MultiChannelImage with shuffled pixel blocks.

# Notes
This function calculates the number of blocks in each channel, retrieves the pixel values, shuffles the blocks, and stores the result back in the image. The shuffling is done by creating a new array with the first two dimensions shuffled. The blocks are then unpatched and the image data is updated.
"""
function shuffle_blocks(img::MultiChannelImage, block_size::Int64)
    # calculate number of blocks
    img_size = [size(img.data[1])[1], size(img.data[1])[2]]
    num_blocks = Int.(div.(img_size, sqrt(block_size)))


    for channel ∈ 1:length(img.channels)
        # retrieve pixel values
        data = patch(img.data[channel], num_blocks[1], num_blocks[2])
        ############################################################
        # shuffle blocks, i.e. only shuffle the order of the blocks
        ############################################################
        # Assume data is a 4D array and shuffle_idx is a 1D array of shuffled indices
        dims = size(data)
        indices = CartesianIndices((dims[1], dims[2]))
        shuffle_idx = indices[randperm(length(indices))]

        # Create a new array with the first two dimensions shuffled
        shuffled_data = similar(data)
        for i in eachindex(shuffle_idx)
            shuffled_data[shuffle_idx[i], :, :, :] = data[indices[i], :, :, :]
        end

        ############################################################
        # unpatch
        ############################################################
        img.data[channel] = unpatch(shuffled_data, Int.(img.pixel_size))
    end
    return img
end


################################################################################
#= Excluded in the compiled version
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

=#

################################################################################
# constructor for MultiChannelImage
"""
    MultiChannelImage(name::S, path::Vector{S}, channels::Vector{S} =[]) where {S <: AbstractString}

This constructor creates a new MultiChannelImage.

# Arguments
- `name`: A string representing the name of the image.
- `path`: A Vector of strings representing the paths to the image files for each channel.
- `channels`: A Vector of strings representing the names of the channels. Default is an empty Vector.

# Returns
- `img`: A MultiChannelImage representing the loaded image.

# Notes
This constructor loads the image data from the files at the specified paths, checks if channel names are provided and generates default names if not, calculates the pixel size and Otsu threshold for each channel, and creates a new MultiChannelImage with these values.
"""
function MultiChannelImage(name::S, path::Vector{S}, channels::Vector{S} =[]) where {S <: AbstractString}
    data = load_tiff.(path)

    if isempty(channels)
        @warn "No channel names provided, using default names"
        channels = ["channel_$i" for i in 1:length(data)]
    end

    @assert length(data) == length(channels) "Number of channels must match the number of image files"
    pixel_size = size(data[1])
    otsu_threshold = Images.otsu_threshold.(data)
    
    return MultiChannelImage(
        data, 
        channels, 
        name, 
        path, 
        pixel_size, 
        otsu_threshold
        )
end