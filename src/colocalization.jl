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
    Patches the image into num_patches x num_patches patches.
    Return a 4D array with the patches.
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
    Patches the image into num_patches_x x num_patches_y patches.

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
    unpatch(patches::Array{Union{Float64, Missing}, 4})
    Unpatches the patches into an image.
    Return the image.
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

"""
    _apply_mask!(img::MultiChannelImage, mask::Matrix{Bool})
    Apply the mask to the image.
    Return the masked image.
"""

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
            length(a) <= 15 ? ρ[i, j] = missing : ρ[i, j] = cor(a, b)
        end
    end
    return ρ
end

