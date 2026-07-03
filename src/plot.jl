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

"""
    minmax_norm!(img::Matrix{Float64})

This function normalizes a 2D image matrix to the range [0, 1].

# Arguments
- `img`: A 2D matrix of Float64 representing the image to be normalized.

# Returns
- `img`: The input image matrix, modified in-place, where each pixel value has been normalized to the range [0, 1].

# Notes
This function normalizes the image by subtracting the minimum pixel value from each pixel, and then dividing each pixel by the range of pixel values (maximum - minimum). The normalization is performed in-place, modifying the input image matrix directly.
"""
function minmax_norm!(img::Matrix{Float64})
    lo, hi = extrema(img)  # Single pass for min and max
    img .= (img .- lo) ./ (hi - lo)  # True in-place operation
    return img
end


"""
    cm_to_px(cm::Float64, dpi::Int)::Int

This function converts a measurement from centimeters to pixels, given a specific dots per inch (dpi) value.

# Arguments
- `cm`: A Float64 representing the measurement in centimeters.
- `dpi`: An integer representing the dots per inch (dpi) value.

# Returns
- `px`: An integer representing the measurement in pixels.

# Notes
This function first converts the measurement from centimeters to inches, then converts the measurement from inches to pixels using the dpi value. The result is rounded to the nearest integer.
"""
function cm_to_px(cm, dpi)
    inch = cm / 2.54  # Convert cm to inches
    px = round(Int, inch * dpi)  # Convert inches to pixels
    return px
end

"""
    calculate_font_size(
    resolution::Tuple{Int, Int}, 
    scale_factor::Float64
    )

This function calculates the font size based on the resolution of the plot and a scale factor.

# Arguments
- `resolution`: A tuple of two integers representing the width and height of the plot in pixels.
- `scale_factor`: A Float64 representing the scale factor to be applied to the resolution to calculate the font size.

# Returns
- `font_size`: An integer representing the calculated font size.

# Notes
This function calculates the font size by taking the minimum of the width and height of the plot, multiplying it by the scale factor, and rounding the result to the nearest integer.
"""
function calculate_font_size(resolution, scale_factor)
    return round(Int, min(resolution...) * scale_factor)
end

"""
    plot(
    img::AbstractMultiChannelImage,
    num_patches::T,
    cor_channel::Vector{T} = [2, 3];
    scale_channels::Bool = true,
    file::String = "test.png";
    channel_for_plot::Vector{T} = [1, 2, 3],
    save_to_file::Bool = true,
    cor_method::Symbol = :pearson
    ) where {T <: Int}

    Plot the image with the patches and the correlation of each patch.

    # Arguments
    - `img::AbstractMultiChannelImage`: The image to plot.
    - `num_patches::Int64`: The number of patches to use.
    - `cor_channel::Vector{Int64}`: The channels to use for the calculation of correlation.
    - `scale_channels::Bool`: Whether to scale the channels intensities to the range [0, 1].
    - `file::String`: The file to save the plot to. Must end with .png or .svg.
    - `channel_for_plot::Vector{Int64}`: The channels to use for plotting. Only 1-3 channels are supported.
        The first channel is plotted in blue, the second in green and the third in red.
    - `save_to_file::Bool`: Whether to save the plot to the file.
    - `cor_method::Symbol`: The method to use for the calculation of correlation. Default is :pearson. Other options are :spearman and :kendall.

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
    img::AbstractMultiChannelImage,
    num_patches::T,
    cor_channel::Vector{T} = [2, 3];
    scale_channels::Bool = true,
    file::String = "test.png",
    channel_for_plot::Vector{T} = [1, 2, 3],
    save_to_file::Bool = true,
    cor_method::Symbol = :pearson
    ) where {T <: Int}
    # check that file ends with .png or .svg
    endswith(file, ".png") || endswith(file, ".svg") || @error "File must end with .png or .svg"
    # check that the number of channels to plot is between 1 and 3
    length(channel_for_plot) <= 3 || @error "Only 1-3 channels are supported for plotting"
    # check that the number of channels_to_plot is <= the number of channels
    length(channel_for_plot) <= num_channels(img) || @error "The number of channels to plot must be <= the number of channels in the image"

    # get image data
    img_data = image_data(img)

    # get patch size
    patches = patch(img_data[1], num_patches)
    patch_size = size(patches)[3:4]

    # extract the individual channels and scale them
    n_channels = num_channels(img)
    scale_channels ? c = [minmax_norm!(img_data[i]) for i in 1:n_channels] : c = img_data

    # plot the image
    if length(channel_for_plot) == 3
        img_view = Images.colorview(
            Images.RGB,
            c[channel_for_plot[3]],
            c[channel_for_plot[2]],
            c[channel_for_plot[1]]
            )

    elseif length(channel_for_plot)  == 2
        img_view = Images.colorview(
            Images.RGB,
            c[channel_for_plot[2]],
            c[channel_for_plot[1]],
            Images.zeroarray
            )
    elseif length(channel_for_plot)  == 1
        img_view = Images.colorview(
            Images.Gray,
            c[channel_for_plot[1]],
            Images.zeroarray
            )
    end

    fig = GLMakie.Figure(size = (size(img_data[1])[2,],size(img_data[1])[1,]), backgroundcolor = :black)
    ax1 = GLMakie.Axis(fig[1, 1], aspect = GLMakie.DataAspect(), yreversed = true)
    GLMakie.image!(ax1, img_view')

    # add lines to the image to separate the patches
    for i in 0:num_patches
        GLMakie.hlines!([i*patch_size[1]], color = :white, alpha = 0.5, linestyle = :dash, linewidth = 1)
        GLMakie.vlines!([i*patch_size[2]], color = :white, alpha = 0.5, linestyle = :dash, linewidth = 1)
    end

    ######################################
    # add correlation values to the image
    ######################################
    x = img_data[cor_channel[1]]
    y = img_data[cor_channel[2]]
    # assert that the images are the same size
    @assert size(x) == size(y) "Images are not the same size"
    # patch the image
    x, y = patch.([x, y], num_patches)
    # calculate the correlation for each patch
    ρ = correlation(x, y, method = cor_method)

    # add the correlation values to the image
    for patch ∈ CartesianIndices(ρ)
        ismissing(ρ[patch]) && continue
        GLMakie.text!(
            ax1,
            (patch[2]-1)*patch_size[2] + patch_size[2]/2,
            (patch[1]-1)*patch_size[1] + patch_size[1]/2,
            text = string(round(ρ[patch],digits = 2)),
            color = :white, fontsize = 9	
            )
    end

    # save the image
    save_to_file && GLMakie.save(file, fig)
    return fig
end

"""
    plot_mask(img::AbstractMultiChannelImage,file::String = "test.png")
    Plot the mask of the image with the channels separated by lines.
    Save the plot to the file and also return it.

    # Arguments
    - `img::AbstractMultiChannelImage`: The image to plot.
    - `file::String`: The file to save the plot to. Must end with .png or .svg.
"""
function plot_mask(img::AbstractMultiChannelImage,file::String = "mask.png")
    mask = _calculate_mask(img)

    # initialize subplots
    plt = []

    for channel_mask ∈ mask
        push!(plt, Images.colorview(Images.Gray, channel_mask))
    end

    # combine the subplots
    plt = Images.hcat(plt...)

    # without x and y axis ticks and labels and with a black background
    fig = GLMakie.Figure(
        size = (0.3*size(plt)[2,],0.3*size(plt)[1,]),
        background_color = :black
        )

    ax1 = GLMakie.Axis(fig[1, 1], aspect = GLMakie.DataAspect(), yreversed = true)
    GLMakie.image!(ax1, plt')

    # get image data and channel names
    img_data = image_data(img)
    channels = channel_names(img)

    # add lines to the image to separate the channels
    for i in 0:length(channels)
        GLMakie.vlines!([i*size(img_data[1])[2,]], color = :white, width = 2)
    end

    # add the channel names to the bottom corner of each channel
    for i in 0:(length(channels)-1)
        GLMakie.text!(
            ax1,
            (i*size(img_data[1])[2,] +100), 100,
            text = string(channels[i+1]),
            color = :yellow, fontsize = 20.0
            )
    end
   
    # save
    GLMakie.save(file, fig)

    return fig
end

"""
    _local_correlation_plot(
    img::AbstractMultiChannelImage,
    channel_for_plot::Vector{Int},
    num_patches::Int,
    cor_channel::Vector{Int},
    cor_method::Symbol = :pearson
    )::Tuple{Matrix{Float64}, Tuple{Int, Int}}

This is a low-level function that calculates and returns the local correlation between two channels of a multi-channel image.

# Arguments
- `img`: An AbstractMultiChannelImage representing the image.
- `channel_for_plot`: A Vector of integers representing the channels to be plotted.
- `num_patches`: An integer representing the number of patches to be analyzed.
- `cor_channel`: A Vector of two integers representing the channels for which the correlation is to be calculated.
- `cor_method`: A Symbol representing the method to be used for the calculation of correlation. Default is :pearson. Other options are :spearman and :kendall.

# Returns
- `ρ`: A matrix of Float64 representing the local correlation values for each patch.
- `patch_size`: A tuple of two integers representing the width and height of the patches in pixels.

# Errors
- Throws a warning if the number of channels is not 3.
- Throws a warning if the patch size is too small or too big for local correlation.
- Throws an error if the images are not the same size.

# Notes
This function checks the number of channels, calculates the patch size, checks the patch size, patches the image, and calculates the correlation for each patch. This function should not be called directly. Use `local_correlation_plot()` instead.
"""
function _local_correlation_plot(img, channel_for_plot, num_patches, cor_channel; cor_method = :pearson)
    # check that the number of channels is 3
    length(channel_for_plot) <= 3 || @warn "This function is currently only implemented for 3 channels"

    # get patch size
    patches = patch(img.data[1], num_patches)
    patch_size = size(patches)[3:4]

    # check that the patch size is reasonable
    #patch_size[1] * patch_size[2] > 10 || @warn "Patch size is too small for local correlation. A size between 10 and 100 px is recommended."
    #patch_size[1] * patch_size[2] < 100 || @warn "Patch size is too big for local correlation. A size between 10 and 100 px is recommended."
    #@info "Patch size is $(patch_size[1]) x $(patch_size[2]) px = $(patch_size[1] * patch_size[2]) px²"    

    # patch the image and calculate the correlation
    x = img.data[cor_channel[1]]
    y = img.data[cor_channel[2]]
    # assert that the images are the same size
    @assert size(x) == size(y) "Images are not the same size"
    # patch the image
    x, y = patch.([x, y], num_patches)
    # calculate the correlation for each patch
    ρ = correlation(x, y, method = cor_method)

    return ρ, patch_size
end


"""
    local_correlation_plot(
    img::AbstractMultiChannelImage,
    num_patches::Int,
    cor_channel::Vector{Int} = [2, 3];
    channel_for_plot::Vector{Int} = [1, 2, 3],
    save::Bool = true,
    file::String = "local_correlation.png",
    cor_method::Symbol = :pearson
    )

This function generates a local correlation plot for a multi-channel image.

# Arguments
- `img`: An AbstractMultiChannelImage representing the image.
- `num_patches`: An integer representing the number of patches to be analyzed.
- `cor_channel`: A Vector of two integers representing the channels for which the correlation is to be calculated. Default is [2, 3].
- `channel_for_plot`: A Vector of integers representing the channels to be plotted. Default is [1, 2, 3].
- `save`: A boolean indicating whether to save the plot to a file. Default is true.
- `file`: A string representing the filename for the output file. Default is "local_correlation.png".
- `cor_method`: A Symbol representing the method to be used for the calculation of correlation. Default is :pearson. Other options are :spearman and :kendall.

# Returns
- `fig`: A Figure object representing the generated plot.

# Errors
- Throws a warning if no patches with a successful correlation calculation exist.
- Throws a warning if no local correlation plot could be generated.

# Notes
This function calculates the local correlation for each patch, checks that patches with a successful correlation calculation exist, tries a different patch size if necessary, plots the local correlation, and saves the plot to a file if specified.
"""
function local_correlation_plot(
    img::AbstractMultiChannelImage,
    num_patches::I,
    cor_channel::Vector{I} = [2, 3];
    channel_for_plot::Vector{I} = [1, 2, 3],
    save::Bool = true,
    file::String = "local_correlation.png",
    cor_method::Symbol = :pearson
    ) where {I <: Int}

    ρ, patch_size = _local_correlation_plot(img, channel_for_plot, num_patches, cor_channel, cor_method = cor_method)

    # check that patches with a successful correlation calculation exist
    # Cache ismissing count to avoid repeated calculations
    num_missing = count(ismissing, ρ)
    total_patches = size(ρ, 1) * size(ρ, 2)
    if num_missing == total_patches
        @warn "No patches with a successful correlation calculation exist at $num_patches patches. A different patch size is tried."
        # try a different patch
        patch_number_range = reverse(collect(10:10:num_patches))
        for pn ∈ patch_number_range
            ρ, patch_size = ProteinCoLoc._local_correlation_plot(img, channel_for_plot, pn, cor_channel)
            num_missing = count(ismissing, ρ)
            total_patches = size(ρ, 1) * size(ρ, 2)
            if num_missing < total_patches
                @info "Found patches with a successful correlation calculation at $pn patches"
                num_patches = pn
                break
            end
            if num_missing == total_patches && pn == patch_number_range[end]
                @warn "No local correlation plot could be generated."
                return nothing
            end
        end    
    end

    # plot the local correlation - direct range broadcast instead of comprehension
    x = (1:num_patches) .* patch_size[2]
    y = (1:num_patches) .* patch_size[1]
    z = [ρ[i, j] for i in 1 : (num_patches), j in 1 : (num_patches)]

    # replace missing values with 0
    z = replace(z, missing => -2.0)

    # plot the local correlation with GLMakie
    GLMakie.activate!()
    fig = GLMakie.Figure(fontsize = 12)
    ax = GLMakie.Axis(
        fig[1, 1],
        xlabel = "x position [px]",
        ylabel = "y position [px]",
        title = "Colocalization between the $(channel_names(img)[cor_channel[1]]) and $(channel_names(img)[cor_channel[2]]) channel",
        yreversed = false,
        limits = (0.0, nothing, 0.0, nothing)
        )

    co = GLMakie.contourf!(
        ax, x,y,rotr90(z), levels = range(-1.001,1.0001,50),
        colormap = :BrBG_10
        )
        
    GLMakie.Colorbar(fig[1,2], co, tellheight = false, label = "ρ")
    GLMakie.resize_to_layout!(fig)

    if save == true
        # save the plot
        GLMakie.save(file,fig)
    end
    # return the plot
    return fig
end

