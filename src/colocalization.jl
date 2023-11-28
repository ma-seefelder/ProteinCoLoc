#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel
E-Mail: proteincoloc@protonmail.com
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

# function to patch the image
"""
    patch(img::Array{Float64, 2}, num_patches::Int64)

This function divides an image into a square grid of patches.

# Arguments
- `img`: A 2D array representing the image to be patched.
- `num_patches`: An integer representing the number of patches along each axis. The image will be divided into `num_patches` x `num_patches` patches.

# Returns
- `patches`: A 4D array representing the patches of the image. The first two dimensions are the number of patches along the x and y axes, and the last two dimensions are the size of each patch.

# Notes
This function calculates the size of the patches, trims the image to fit an exact number of patches, initializes a 4D array for the patches, loops over the patches, and assigns the corresponding part of the image to each patch.
"""
function patch(img::Array{Float64, 2}, num_patches::Int64)
    # calculate the size of the patches
    patch_size_x = Int64(floor(size(img, 1) / num_patches))
    patch_size_y = Int64(floor(size(img, 2) / num_patches))
    
    # cut the image to the correct size
    img = img[1:num_patches*patch_size_x, 1:num_patches*patch_size_y]

    # initialize the patches
    patches = zeros(Union{Float64, Missing}, num_patches, num_patches, patch_size_x, patch_size_y)
    # loop over the patches
    for i in 1:num_patches
        for j in 1:num_patches
            patches[i, j, :, :] = img[(i-1)*patch_size_x+1:i*patch_size_x, (j-1)*patch_size_y+1:j*patch_size_y]
        end
    end
    return patches
end

"""
    patch(img::Array{Float64,2}, num_patches_x::Int64, num_patches_y::Int64)

This function divides an image into a specified number of patches.

# Arguments
- `img`: A 2D array representing the image to be patched.
- `num_patches_x`: An integer representing the number of patches along the x-axis.
- `num_patches_y`: An integer representing the number of patches along the y-axis.

# Returns
- `patches`: A 4D array representing the patches of the image. The first two dimensions are the number of patches along the x and y axes, and the last two dimensions are the size of each patch.

# Notes
This function calculates the size of the patches, initializes a 4D array for the patches, loops over the patches, and assigns the corresponding part of the image to each patch.
"""
function patch(img::Array{Float64,2}, num_patches_x::Int64, num_patches_y::Int64)
    # calculate the size of the patches
    patch_size_x = Int64(floor(size(img, 1) / num_patches_x))
    patch_size_y = Int64(floor(size(img, 2) / num_patches_y))
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

"""
    unpatch(patches::Array{Union{Float64, Missing}, 4}, img_size::Tuple{Int64, Int64})

This function reconstructs an image from its patches.

# Arguments
- `patches`: A 4D array representing the patches of the image. The first two dimensions are the number of patches along the x and y axes, and the last two dimensions are the size of each patch.
- `img_size`: A tuple representing the size of the original image.

# Returns
- `img`: A 2D array representing the reconstructed image.

# Notes
This function iterates over the patches, places each patch in the correct position in the image, and returns the reconstructed image.
"""
function unpatch(patches::Array{Union{Float64, Missing}, 4}, img_size::Tuple{Int64, Int64})
    # Get the dimensions of the patches
    num_patches_x, num_patches_y, patch_size_x, patch_size_y = size(patches)
    # Initialize the image
    img = zeros(Float64, img_size[1], img_size[2])
    # Loop over the patches
    for i in 1:num_patches_x
        for j in 1:num_patches_y
            img[(i-1)*patch_size_x+1:i*patch_size_x, (j-1)*patch_size_y+1:j*patch_size_y] = patches[i, j, :, :]
        end
    end
    return img
end

######################################################################
# function to calculate the Pearson's correlation
######################################################################

"""
    _exclude_zero!(a::Vector{T}, b::Vector{T}) where T <: Number

This function excludes zero and NaN values from two vectors.

# Arguments
- `a`: A vector of numbers.
- `b`: A vector of numbers.

# Returns
- This function does not return a value. It modifies the input vectors `a` and `b` in-place.

# Notes
This function finds the indices of zero and NaN values in both vectors, creates a union of these indices, and deletes the values at these indices from both vectors.
"""
function _exclude_zero!(a::Vector{T}, b::Vector{T}) where T <: Number
    # get the indices of the zero and NaN values
    zero_indices = sort(union(findall(iszero, a), findall(isnan, a), findall(iszero, b), findall(isnan, b)))

    # exclude all values at the indices in zero_indices
    deleteat!(a, zero_indices)
    deleteat!(b, zero_indices)
end

"""
    _exclude_zero!(a::Vector{Union{T,Missing}}, b::Vector{Union{T,Missing}}) where T <: Number

This function excludes zero and NaN values from two vectors.

# Arguments
- `a`: A vector of numbers or missing values.
- `b`: A vector of numbers or missing values.

# Returns
- This function does not return a value.

# Notes
This function replaces missing values with zero, finds the indices of zero and NaN values in both vectors, creates a union of these indices, and deletes the values at these indices from both vectors.
"""
function _exclude_zero!(a::Vector{Union{T,Missing}}, b::Vector{Union{T,Missing}}) where T <: Number
    # change missing values to zero
    replace!(a, missing => 0)
    replace!(b, missing => 0)

    # get the indices of the zero and NaN values
    zero_indices = sort(union(findall(iszero, a), findall(isnan, a), findall(iszero, b), findall(isnan, b)))

    # exclude all values at the indices in zero_indices
    deleteat!(a, zero_indices)
    deleteat!(b, zero_indices)
end

"""
    correlation(x::Array{Float64, 4}, y::Array{Float64, 4})

This function calculates the correlation between two 4D arrays.

# Arguments
- `x`: A 4D array representing the first set of data.
- `y`: A 4D array representing the second set of data.

# Returns
- `ρ`: A 2D array representing the correlation for each patch.

# Notes
This function iterates over the patches in the 4D arrays, excludes zeros from the data, and calculates the correlation between the patches. If the length of the data in a patch is less than or equal to 15, the correlation is set to missing.
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
            length(a) <= 15 ? ρ[i, j] = missing : ρ[i, j] = cor(a, b)
        end
    end
    return ρ
end

