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

module DataTypes
export MultiChannelImage

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

# constructor for MultiChannelImage
function MultiChannelImage(name::S, path::Vector{S}, channels::Vector{S} =[]) where {S <: AbstractString}
    data = LoadImages.load_tiff.(path)

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

end # module DataTypes
