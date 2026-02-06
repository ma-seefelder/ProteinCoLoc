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
abstract type AbstractMultiChannelImage end

# === Interface-Funktionen (erforderlich für alle Subtypen) ===
"""
    image_data(img::AbstractMultiChannelImage) -> Vector{Matrix}
Gibt die Bilddaten als Vector von Matrix pro Kanal zurück.
"""
function image_data(img::AbstractMultiChannelImage)
    hasfield(typeof(img), :data) || throw(ArgumentError(
        "Typ $(typeof(img)) muss `image_data(img)` implementieren oder ein `data` Feld haben"))
    return img.data
end

"""
    channel_names(img::AbstractMultiChannelImage) -> Vector{<:AbstractString}
Gibt die Namen der Kanäle zurück.
"""
function channel_names(img::AbstractMultiChannelImage)
    hasfield(typeof(img), :channels) || throw(ArgumentError(
        "Typ $(typeof(img)) muss `channel_names(img)` implementieren oder ein `channels` Feld haben"))
    return img.channels
end

"""
    image_name(img::AbstractMultiChannelImage) -> AbstractString
Gibt den Namen des Bildes zurück.
"""
function image_name(img::AbstractMultiChannelImage)
    hasfield(typeof(img), :name) || throw(ArgumentError(
        "Typ $(typeof(img)) muss `image_name(img)` implementieren oder ein `name` Feld haben"))
    return img.name
end

"""
    image_paths(img::AbstractMultiChannelImage) -> Vector{<:AbstractString}
Gibt die Dateipfade für jeden Kanal zurück.
"""
function image_paths(img::AbstractMultiChannelImage)
    hasfield(typeof(img), :path) || throw(ArgumentError(
        "Typ $(typeof(img)) muss `image_paths(img)` implementieren oder ein `path` Feld haben"))
    return img.path
end

"""
    pixel_dimensions(img::AbstractMultiChannelImage) -> Tuple{Integer, Integer}
Gibt die Pixeldimensionen (Breite, Höhe) zurück.
"""
function pixel_dimensions(img::AbstractMultiChannelImage)
    hasfield(typeof(img), :pixel_size) || throw(ArgumentError(
        "Typ $(typeof(img)) muss `pixel_dimensions(img)` implementieren oder ein `pixel_size` Feld haben"))
    return img.pixel_size
end

"""
    otsu_thresholds(img::AbstractMultiChannelImage) -> Vector{<:AbstractFloat}
Gibt die Otsu-Schwellenwerte für jeden Kanal zurück.
"""
function otsu_thresholds(img::AbstractMultiChannelImage)
    hasfield(typeof(img), :otsu_threshold) || throw(ArgumentError(
        "Typ $(typeof(img)) muss `otsu_thresholds(img)` implementieren oder ein `otsu_threshold` Feld haben"))
    return img.otsu_threshold
end

"""
    num_channels(img::AbstractMultiChannelImage) -> Int
Gibt die Anzahl der Kanäle zurück.
"""
num_channels(img::AbstractMultiChannelImage) = length(channel_names(img))

"""
    mutable struct MultiChannelImage{T <: Union{Missing, Float64}, S <: AbstractString, F <: AbstractFloat}
        data::Vector{Matrix{T}}
        channels::Vector{S}
        name::S
        path::Vector{S}
        pixel_size::Tuple{I, I}
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
mutable struct MultiChannelImage{T <: Union{Missing, Float64}, S <: AbstractString, F <: AbstractFloat, I <: Integer} <: AbstractMultiChannelImage
    data::Vector{Matrix{T}}
    channels::Vector{S}
    name::S
    path::Vector{S}
    pixel_size::Tuple{I, I}
    otsu_threshold::Vector{F}

    function MultiChannelImage(data::Vector{Matrix{T}}, channels::Vector{S}, name::S, path::Vector{S}, pixel_size::Tuple{I, I}, otsu_threshold::Vector{F}) where {T <: Union{Missing, Float64}, S <: AbstractString, F <: AbstractFloat, I <: Int}
        @assert length(data) == length(channels) "Length of data and channels vectors must match"
        @assert length(otsu_threshold) == length(channels) "Length of otsu_threshold vector must match number of channels"
        new{T, S, F, I}(data, channels, name, path, pixel_size, otsu_threshold)
    end
end

"""
    mutable struct MultiChannelImageStack{T <: AbstractMultiChannelImage, S <: AbstractString}
        img::Vector{T}
        name::S

    This mutable struct represents a stack of multi-channel images.

    # Fields
    - `img`: A Vector of AbstractMultiChannelImage representing the images in the stack.
    - `name`: A string representing the name of the image stack.

    # Notes
    This struct is used to store and manipulate stacks of multi-channel images in the ProteinCoLoc package.
    Use `length(stack)` to get the number of images.
    """
mutable struct MultiChannelImageStack{T <: AbstractMultiChannelImage, S <: AbstractString}
    img::Vector{T}
    name::S
end

Base.length(stack::MultiChannelImageStack) = length(stack.img)

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
function Base.getindex(stack::MultiChannelImageStack, i::Integer)
    @boundscheck checkbounds(
      Bool, stack.img, i
    ) || throw(BoundsError(stack, i))
    return @inbounds stack.img[i]
end

"""
    Base.iterate(stack::MultiChannelImageStack, state=1)

This function overrides the base iterate function to iterate over the images in a MultiChannelImageStack.

# Arguments
- `stack`: A MultiChannelImageStack over which to iterate.
- `state`: An Int64 representing the current state of the iteration. Default is 1.

# Returns
- A tuple containing the current image and the next state, or nothing if the end of the stack has been reached.
"""
function Base.iterate(stack::MultiChannelImageStack, state=1)
    state > length(stack) && return nothing
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
    _calculate_mask(img::AbstractMultiChannelImage)

This function calculates a binary mask for each channel of an AbstractMultiChannelImage, where pixels with intensity greater than the Otsu threshold are set to 1 and others are set to 0.

# Arguments
- `img`: An AbstractMultiChannelImage for which to calculate the mask.

# Returns
- `mask`: A Vector of Matrix{Bool} representing the binary mask for each channel.

# Notes
This function iterates over the channels of the image, applies the Otsu threshold to the data of each channel, and stores the resulting binary mask.
"""
function _calculate_mask(img::AbstractMultiChannelImage)
    return [
        ch .> thr
        for (ch, thr) in
            zip(image_data(img), otsu_thresholds(img))
    ]
end


# Apply mask
"""
    apply_mask!(img::AbstractMultiChannelImage)

This function calculates a binary mask for an AbstractMultiChannelImage and applies it to the image.

# Arguments
- `img`: An AbstractMultiChannelImage to which to apply the mask.

# Returns
- `img`: The AbstractMultiChannelImage with the mask applied.

# Notes
This function calls the `_calculate_mask` function to calculate the mask, and then calls the `_apply_mask!` function to apply the mask to the image.
"""
function _apply_mask!(img::AbstractMultiChannelImage, mask)
    for channel ∈ 1:length(img.channels)
        img.data[channel] = img.data[channel] .* mask[channel]
    end
    return img
end

function apply_mask!(img::AbstractMultiChannelImage)
    mask = _calculate_mask(img)
    img = _apply_mask!(img, mask)
end

################################################################################
# shuffle individual pixels in each channel of a MultiChannelImage
"""
    shuffle_pixels(img::AbstractMultiChannelImage)

This function shuffles the pixel values in each channel of an AbstractMultiChannelImage.

# Arguments
- `img`: An AbstractMultiChannelImage whose pixel values are to be shuffled.

# Returns
- `img`: The AbstractMultiChannelImage with shuffled pixel values.

# Notes
This function iterates over the channels of the image, shuffles the pixel values in each channel using the `shuffle!` function from the Julia standard library, and stores the result back in the image.
"""
function shuffle_pixels(img::AbstractMultiChannelImage)
    data = image_data(img)
    for channel ∈ 1:num_channels(img)
        # retrieve pixel values
        data[channel] = shuffle!(data[channel])
    end
    return img
end

# shuffle blocks of pixels in each channel of a MultiChannelImage 
# with block size block_size to accomodate for autocorrelation
# use the patch function to create blocks

"""
    shuffle_blocks!(img::AbstractMultiChannelImage, block_size::Integer)

Shuffle blocks of pixels in-place to preserve spatial autocorrelation.
`block_size` is the size length of each square block in pixels.
Blocks at the image edges that don't fit a full block are included
in the last block of that row/column (variable-size edge blocks).
"""
function shuffle_blocks!(img::AbstractMultiChannelImage, block_size::Integer)
    img_data = image_data(img)
    rows, cols = size(img_data[1])

    # — Compute block boundary ranges once —
    function make_ranges(total, bs)
        starts = 1:bs:total
        return [s:min(s + bs - 1, total) for s in starts]
    end

    row_ranges = make_ranges(rows, block_size)
    col_ranges = make_ranges(cols, block_size)

    # — Flat list of (row_range, col_range) block coordinates —
    blocks = [(rr, cr) for rr in row_ranges for cr in col_ranges]
    n_blocks = length(blocks)

    # — Shuffle once, reuse permutation across channels —
    perm = randperm(n_blocks)

    for ch in eachindex(img_data)
        dest = similar(img_data[ch])  # single allocation per channel

        for (dst_idx, src_idx) in enumerate(perm)
            dst_rr, dst_cr = blocks[dst_idx]
            src_rr, src_cr = blocks[src_idx]
            # Copy block — @views avoids any temporary arrays
            @views dest[dst_rr, dst_cr] .= img_data[ch][src_rr, src_cr]
        end

        img_data[ch] = dest
    end

    return img
end

# Non-mutating convenience wrapper
function shuffle_blocks(img::AbstractMultiChannelImage, block_size::Integer)
    img_copy = deepcopy(img)
    return shuffle_blocks!(img_copy, block_size)
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


    if length(data) != length(channels)
         throw(ArgumentError("Number of channels must match the number of image files"))
    end
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