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

# function to patch the image
"""
    patch(img::AbstractMatrix{T}, num_patches::Integer) where T <: Union{Float64, Missing}

This function divides an image into a square grid of patches.

# Arguments
- `img`: A 2D array representing the image to be patched.
- `num_patches`: An integer representing the number of patches along each axis. The image will be divided into `num_patches` x `num_patches` patches.

# Returns
- `patches`: A 4D array representing the patches of the image. The first two dimensions are the number of patches along the x and y axes, and the last two dimensions are the size of each patch.

# Notes
This function calculates the size of the patches, trims the image to fit an exact number of patches, initializes a 4D array for the patches, loops over the patches, and assigns the corresponding part of the image to each patch.
"""
function patch(img::AbstractMatrix{T}, num_patches::Integer) where T <: Union{Float64, Missing}
    rows, cols = size(img)
    patch_size_x = rows ÷ num_patches
    patch_size_y = cols ÷ num_patches

    # View statt Kopie für Trimming
    trimmed = @view img[1:num_patches*patch_size_x, 1:num_patches*patch_size_y]

    # Union-Typ beibehalten für Kompatibilität
    patches = Array{Union{Float64, Missing}, 4}(undef, num_patches, num_patches, patch_size_x, patch_size_y)

    # Elementweise Kopie mit @inbounds - vermeidet temporäre Arrays
    @inbounds for j in 1:num_patches
        y_start = (j-1)*patch_size_y + 1
        for i in 1:num_patches
            x_start = (i-1)*patch_size_x + 1
            for py in 1:patch_size_y, px in 1:patch_size_x
                patches[i, j, px, py] = trimmed[x_start + px - 1, y_start + py - 1]
            end
        end
    end

    return patches
end

"""
    patch(img::AbstractMatrix{T}, num_patches_x::Integer, num_patches_y::Integer) where T <: Union{Float64, Missing}

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
function patch(img::AbstractMatrix{T}, num_patches_x::Integer, num_patches_y::Integer) where T <: Union{Float64, Missing}
    rows, cols = size(img)
    patch_size_x = rows ÷ num_patches_x
    patch_size_y = cols ÷ num_patches_y

    # View statt Kopie
    trimmed = @view img[1:num_patches_x*patch_size_x, 1:num_patches_y*patch_size_y]

    patches = Array{Union{Float64, Missing}, 4}(undef, num_patches_x, num_patches_y, patch_size_x, patch_size_y)

    @inbounds for j in 1:num_patches_y
        y_start = (j-1)*patch_size_y + 1
        for i in 1:num_patches_x
            x_start = (i-1)*patch_size_x + 1
            for py in 1:patch_size_y, px in 1:patch_size_x
                patches[i, j, px, py] = trimmed[x_start + px - 1, y_start + py - 1]
            end
        end
    end

    return patches
end

"""
    unpatch(patches::AbstractArray{T, 4}, img_size::Tuple{<:Integer, <:Integer}) where T

This function reconstructs an image from its patches.

# Arguments
- `patches`: A 4D array representing the patches of the image. The first two dimensions are the number of patches along the x and y axes, and the last two dimensions are the size of each patch.
- `img_size`: A tuple representing the size of the original image.

# Returns
- `img`: A 2D array representing the reconstructed image.

# Notes
This function iterates over the patches, places each patch in the correct position in the image, and returns the reconstructed image.
"""
function unpatch(patches::AbstractArray{T, 4}, img_size::Tuple{<:Integer, <:Integer}) where T
    num_patches_x, num_patches_y, patch_size_x, patch_size_y = size(patches)

    # Float64 für Ausgabe (kein Missing nötig)
    img = Matrix{Float64}(undef, img_size[1], img_size[2])

    @inbounds for j in 1:num_patches_y
        y_start = (j-1)*patch_size_y + 1
        for i in 1:num_patches_x
            x_start = (i-1)*patch_size_x + 1
            for py in 1:patch_size_y, px in 1:patch_size_x
                img[x_start + px - 1, y_start + py - 1] = patches[i, j, px, py]
            end
        end
    end

    return img
end

######################################################################
# function to calculate the Pearson's correlation
######################################################################

"""
    _exclude_zero(a::Vector{T}, b::Vector{T}) where T <: Number

This function excludes zero and NaN values from two vectors.

# Arguments
- `a`: A vector of numbers.
- `b`: A vector of numbers.

# Returns
- This function returns two vectors of numbers without zero and NaN values.

# Notes
This function finds the indices of zero and NaN values in both vectors, creates a union of these indices, and deletes the values at these indices from both vectors.
"""
function _exclude_zero(a::Vector{T}, b::Vector{T}) where T <: Number
    # get the indices of the zero and NaN values
    zero_indices = sort(union(findall(iszero, a), findall(isnan, a), findall(iszero, b), findall(isnan, b)))

    # exclude all values at the indices in zero_indices
    deleteat!(a, zero_indices)
    deleteat!(b, zero_indices)

    return a,b
end

"""
    _exclude_zero(a::Vector{Union{T,Missing}}, b::Vector{Union{T,Missing}}) where T <: Number

This function excludes zero and NaN values from two vectors.

# Arguments
- `a`: A vector of numbers or missing values.
- `b`: A vector of numbers or missing values.

# Returns
- `a`: A vector of numbers without zero and NaN values.
- `b`: A vector of numbers without zero and NaN values.

# Notes
This function replaces missing values with zero, finds the indices of zero and NaN values in both vectors, creates a union of these indices, and deletes the values at these indices from both vectors.
"""
function _exclude_zero(a::Vector{Union{T,Missing}}, b::Vector{Union{T,Missing}}) where T <: Number
    # change missing values to zero
    replace!(a, missing => 0)
    replace!(b, missing => 0)

    # get the indices of the zero and NaN values
    zero_indices = sort(union(findall(iszero, a), findall(isnan, a), findall(iszero, b), findall(isnan, b)))

    # exclude all values at the indices in zero_indices
    deleteat!(a, zero_indices)
    deleteat!(b, zero_indices)

    a = convert(Vector{T}, a)
    b = convert(Vector{T}, b)

    return a,b
end

"""
    correlation(x::Array{Float64, 4}, y::Array{Float64, 4})

This function calculates the correlation between two 4D arrays.

# Arguments
- `x`: A 4D array representing the first set of data.
- `y`: A 4D array representing the second set of data.
- `method`: A symbol representing the correlation method to be used. The default is `:pearson`. The other options are `:spearman` and `:kendall`.

# Returns
- `ρ`: A 2D array representing the correlation for each patch.

# Notes
This function iterates over the patches in the 4D arrays, excludes zeros from the data, and calculates the correlation between the patches. If the length of the data in a patch is less than or equal to 15, the correlation is set to missing.
"""
function correlation(x::Array{T, 4}, y::Array{T, 4}; method::Symbol = :pearson) where T <: Union{Float64, Missing}
    cor_dict = Dict(:pearson => cor, :spearman => corspearman, :kendall => corkendall)
    haskey(cor_dict, method) || throw(ArgumentError("method must be one of $(keys(cor_dict))"))
    cor_func = cor_dict[method]
    # get the number of patches
    num_patches = size(x, 1)
    # initialize the correlation
    ρ = zeros(Union{Float64, Missing},num_patches, num_patches)
    # loop over the patches
    for i in 1:num_patches
        for j in 1:num_patches
            a = x[i, j, :, :][:] 
            b = y[i, j, :, :][:]
            a,b = _exclude_zero(a,b)
            length(a) <= 15 ? ρ[i, j] = missing : ρ[i, j] = cor_func(a, b)
        end
    end
    return ρ
end

