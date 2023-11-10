## utils.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder

# get_files

function get_images(path::S, nchannels::I, stack_name::S) where {I<:Integer, S<:AbstractString}
    # get files in image path
    files = readdir(path, join=false)
    # remove subdirectories
    files = filter(x -> contains(x,"."), files)
    # determine number of images in folder
    nimages = length(files) / nchannels
    # check if the number of images is an integer 
    nimages != floor(nimages) && error("The number of images in the image folder is not a multiple of the number of channels.")
    # check if the number of images is at least 1
    nimages < 1 && error("The number of images in the image folder must be at least 1.")

    # generate vector wit channel names
    channel_names = string.(collect(1:nchannels))
    # group images
    images = Vector{MultiChannelImage{Float64,String, Float64}}(undef, Int(nimages))

    # retrieve image names
    image_names = Vector{String}()
    for i in split.(files, "_")
        !in(i[1],image_names) ? push!(image_names,i[1]) : nothing
    end

    # load images
    for i in eachindex(images) 
        image_name = image_names[i]
        channel_files = fill("", nchannels)
        for j in 1:nchannels
            file_index = findfirst(file -> occursin(image_name * "_c" * string(j), file), files)
            !isnothing(file_index) || error("The image path is not valid.")
            channel_files[j] = path * "/" * files[file_index]
        end
        images[i] = MultiChannelImage(image_name, channel_files, channel_names)
        # mask image 
        images[i] = ProteinCoLoc._apply_mask!(images[i], ProteinCoLoc._calculate_mask(images[i]))
    end
    # return files
    return MultiChannelImageStack(images, stack_name)
end 