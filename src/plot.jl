## plot.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder
# module PlotImage
using Images
using GLMakie

include("LoadImages.jl")
import .LoadImages

include("colocalization.jl")
import .Colocalization: patch, correlation

import Plots: plot, vline!, hline!, annotate!, savefig, text


minmax_norm!(img::Matrix{Float64}) = (img .- minimum(img)) ./ (maximum(img) - minimum(img))


# load test images
path = ["test_images/c1.tif", "test_images/c2.tif", "test_images/c3.tif"]
img = LoadImages.MultiChannelImage("positive_sample", path, ["blue", "green", "red"])
mask = LoadImages._calculate_mask(img)
LoadImages._apply_mask!(img, mask)
cor_channel::Vector{Int64} = [2, 3]
scale_channels::Bool = true
num_patches::Int64 = 16

"""
    plot(
    img::LoadImages.MultiChannelImage,
    num_patches::Int64,
    cor_channel::Vector{Int64} = [2, 3];
    scale_channels::Bool = true,
    file::String = "test.png"
    )

    Plot the image with the patches and the correlation of each patch.

    # Arguments
    - `img::LoadImages.MultiChannelImage`: The image to plot.
    - `num_patches::Int64`: The number of patches to use.
    - `cor_channel::Vector{Int64}`: The channels to use for the calculation of correlation.
    - `scale_channels::Bool`: Whether to scale the channels intensities to the range [0, 1].
    - `file::String`: The file to save the plot to. Must end with .png or .svg.

    # Example
    ```julia
    path = ["test_images/c1.tif", "test_images/c2.tif", "test_images/c3.tif"]
    img = LoadImages.MultiChannelImage("positive_sample", path, ["blue", "green", "red"])
    mask = LoadImages._calculate_mask(img)
    LoadImages._apply_mask!(img, mask)
    cor_channel = [2, 3]
    scale_channels = true
    num_patches = 16

    plot(img, num_patches, cor_channel, scale_channels, "test.png")
    ```
"""
function plot(
    img::LoadImages.MultiChannelImage,
    num_patches::Int64,
    cor_channel::Vector{Int64} = [2, 3];
    scale_channels::Bool = true,
    file::String = "test.png"
    )
    # check that file ends with .png or .svg
    endswith(file, ".png") || endswith(file, ".svg") || @error "File must end with .png or .svg"

    # check that the number of channels is 3
    channel = length(img.channels)
    channel <= 3 || @warn "This function is currently only implemented for 3 channels"

    # get patch size
    patches = patch(img.data[1], num_patches)
    patch_size = size(patches)[3:4]

    # extract the individual channels and scale them
    scale_channels ? c = [minmax_norm!(img.data[i]) for i in 1:channel] : c = img.data

    # plot the image
    img_view = Images.colorview(Images.RGB, c[3], c[2], c[1])

    # plot image with Plots.jl package
    # without x and y axis ticks and labels and with a black background
    plt = plot(
        img_view,
        ticks = false,
        border = :none,
        legend = false,
        grid = false,
        size = (size(img.data[1])[2,],size(img.data[1])[1,]),
        background_color = :black,
        dpi = 96
        )


    # add lines to the image to separate the patches
    for i in 0:num_patches
        hline!([i*patch_size[1]], color = :white, linealpha = 0.5, linestyle = :dash)
        vline!([i*patch_size[2]], color = :white, linealpha = 0.5, linestyle = :dash)
    end

    ######################################
    # add correlation values to the image
    ######################################
    x = img.data[cor_channel[1]]
    y = img.data[cor_channel[2]]
    # assert that the images are the same size
    @assert size(x) == size(y) "Images are not the same size"
    # patch the image
    x, y = patch.([x, y], num_patches)
    # calculate the correlation for each patch
    ρ = correlation(x, y)

    # add the correlation values to the image
    for patch ∈ CartesianIndices(ρ)
        ismissing(ρ[patch]) && continue
        annotate!(
            (patch[2]-1)*patch_size[2] + patch_size[2]/2,
            (patch[1]-1)*patch_size[1] + patch_size[1]/2,
            text(round(ρ[patch],digits = 2), :white, :center, 8)
            )
    end

    # save the image
    savefig(file)

    return plt
end

"""
    plot_mask(img::LoadImages.MultiChannelImage,file::String = "test.png")
    Plot the mask of the image with the channels separated by lines.
    Save the plot to the file and also return it.

    # Arguments
    - `img::LoadImages.MultiChannelImage`: The image to plot.
    - `file::String`: The file to save the plot to. Must end with .png or .svg.
"""
function plot_mask(img::LoadImages.MultiChannelImage,file::String = "test.png")
    mask = LoadImages._calculate_mask(img)

    # initialize subplots
    plt = []

    for channel_mask ∈ mask
        push!(plt, Images.colorview(Images.Gray, channel_mask))
    end

    # combine the subplots
    plt = Images.hcat(plt...)

    # plot image with Plots.jl package
    # without x and y axis ticks and labels and with a black background
    plt = plot(
        plt,
        ticks = false,
        border = :none,
        legend = false,
        grid = false,
        size = (3*size(img.data[1])[2,],size(img.data[1])[1,]),
        background_color = :black,
        dpi = 96
        )

    # add lines to the image to separate the channels
    for i in 0:length(img.channels)
        vline!([i*size(img.data[1])[2,]], color = :white, width = 2)
    end

    # save
    savefig(plt, file)

    return plt
end

# local correlation plot
function local_correlation_plot(
    img::LoadImages.MultiChannelImage,
    num_patches::Int64,
    cor_channel::Vector{Int64} = [2, 3]
    )

    # check that the number of channels is 3
    channel = length(img.channels)
    channel <= 3 || @warn "This function is currently only implemented for 3 channels"

    # get patch size
    patches = patch(img.data[1], num_patches);
    patch_size = size(patches)[3:4]

    # check that the patch size is reasonable
    patch_size[1] * patch_size[2] > 10 || @warn "Patch size is too small for local correlation. A size between 10 and 100 px is recommended."
    patch_size[1] * patch_size[2] < 100 || @warn "Patch size is too big for local correlation. A size between 10 and 100 px is recommended."

    # patch the image and calculate the correlation
    x = img.data[cor_channel[1]]
    y = img.data[cor_channel[2]]
    # assert that the images are the same size
    @assert size(x) == size(y) "Images are not the same size"
    # patch the image
    x, y = patch.([x, y], num_patches)
    # calculate the correlation for each patch
    ρ = correlation(x, y)

    # plot the local correlation
    x = [i for i in 1 : (num_patches-1)] .* patch_size[2]
    y = [i for i in 1 : (num_patches-1)] .* patch_size[1]
    z = [ρ[i, j] for i in 1 : (num_patches-1), j in 1 : (num_patches-1)]

    # replace missing values with 0
    z = replace(z, missing => 0)

    # replace all values below 0.1 with 0
    for i in 1 : (num_patches-1)
        for j in 1 : (num_patches-1)
            if z[i, j] < 0.1
                z[i, j] = 0
            end
        end
    end

    # plot the local correlation with GLMakie
    GLMakie.activate!()
    fig = GLMakie.Figure(resolution = (size(img.data[1])[2,],size(img.data[1])[1,]))
    ax = GLMakie.Axis(
        fig[1, 1],
        xlabel = "x position [px]",
        ylabel = "y position [px]",
        title = "Colocalization of $(img.channels[cor_channel[1]]) and $(img.channels[cor_channel[2]])",
        subtitle = "Only ρ ≥ 0.1 are shown.",
        aspect = DataAspect(),
        yreversed = false
        )

    co = GLMakie.contourf!(ax, x,y,rotr90(z), levels = range(0.1,1,20))
    Colorbar(fig[1,2], co, tellheight = false, label = "|ρ|")

    # show the plot
    GLMakie.display(fig)



end


