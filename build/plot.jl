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
minmax_norm!(img::Matrix{Float64}) = (img .- minimum(img)) ./ (maximum(img) - minimum(img))


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
    img::MultiChannelImage,
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
    - `img::MultiChannelImage`: The image to plot.
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
    img = MultiChannelImage("positive_sample", path, ["blue", "green", "red"])
    mask = _calculate_mask(img)
    _apply_mask!(img, mask)
    cor_channel = [2, 3]
    scale_channels = true
    num_patches = 16

    plot(img, num_patches, cor_channel, scale_channels, "test.png")
    ```
"""
function plot(
    img::MultiChannelImage,
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
    length(channel_for_plot) <= length(img.channels) || @error "The number of channels to plot must be <= the number of channels in the image"

    # get patch size
    patches = patch(img.data[1], num_patches)
    patch_size = size(patches)[3:4]

    # extract the individual channels and scale them
    channel = length(img.channels)
    scale_channels ? c = [minmax_norm!(img.data[i]) for i in 1:channel] : c = img.data

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

    fig = GLMakie.Figure(size = (size(img.data[1])[2,],size(img.data[1])[1,]), backgroundcolor = :black)
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
    x = img.data[cor_channel[1]]
    y = img.data[cor_channel[2]]
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
    plot_mask(img::MultiChannelImage,file::String = "test.png")
    Plot the mask of the image with the channels separated by lines.
    Save the plot to the file and also return it.

    # Arguments
    - `img::MultiChannelImage`: The image to plot.
    - `file::String`: The file to save the plot to. Must end with .png or .svg.
"""
function plot_mask(img::MultiChannelImage,file::String = "mask.png")
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

    # add lines to the image to separate the channels
    for i in 0:length(img.channels)
        GLMakie.vlines!([i*size(img.data[1])[2,]], color = :white, width = 2)
    end

    # add the channel names to the bottom corner of each channel
    for i in 0:(length(img.channels)-1)
        GLMakie.text!(
            ax1, 
            (i*size(img.data[1])[2,] +100), 100,
            text = string(img.channels[i+1]),
            color = :yellow, fontsize = 20.0
            )
    end
   
    # save
    GLMakie.save(file, fig)

    return fig
end

"""
    _local_correlation_plot(
    img::MultiChannelImageStack, 
    channel_for_plot::Vector{Int}, 
    num_patches::Int, 
    cor_channel::Vector{Int},
    cor_method::Symbol = :pearson
    )::Tuple{Matrix{Float64}, Tuple{Int, Int}}

This is a low-level function that calculates and returns the local correlation between two channels of a multi-channel image.

# Arguments
- `img`: A MultiChannelImageStack representing the image.
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
    img::MultiChannelImage,
    num_patches::Int,
    cor_channel::Vector{Int} = [2, 3];
    channel_for_plot::Vector{Int} = [1, 2, 3],
    save::Bool = true,
    file::String = "local_correlation.png",
    cor_method::Symbol = :pearson
    )

This function generates a local correlation plot for a multi-channel image.

# Arguments
- `img`: A MultiChannelImage representing the image.
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
    img::MultiChannelImage,
    num_patches::I,
    cor_channel::Vector{I} = [2, 3];
    channel_for_plot::Vector{I} = [1, 2, 3],
    save::Bool = true,
    file::String = "local_correlation.png",
    cor_method::Symbol = :pearson
    ) where {I <: Int}

    ρ, patch_size = _local_correlation_plot(img, channel_for_plot, num_patches, cor_channel, cor_method = cor_method)
    
    # check that patches with a successful correlation calculation exist
    if sum(ismissing.(ρ)) == size(ρ)[1] * size(ρ)[2]
        @warn "No patches with a successful correlation calculation exist at $num_patches patches. A different patch size is tried."
        # try a different patch
        patch_number_range = reverse(collect(10:10:num_patches))
        for pn ∈ patch_number_range
            ρ, patch_size = ProteinCoLoc._local_correlation_plot(img, channel_for_plot, pn, cor_channel)
            if sum(ismissing.(ρ)) < size(ρ)[1] * size(ρ)[2]
                @info "Found patches with a successful correlation calculation at $pn patches"
                num_patches = pn
                break
            end
            if sum(ismissing.(ρ)) == size(ρ)[1] * size(ρ)[2] && pn == patch_number_range[end]
                @warn "No local correlation plot could be generated."
                return nothing
            end
        end    
    end

    # plot the local correlation
    x = [i for i in 1 : (num_patches)] .* patch_size[2]
    y = [i for i in 1 : (num_patches)] .* patch_size[1]
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
        title = "Colocalization between the $(img.channels[cor_channel[1]]) and $(img.channels[cor_channel[2]]) channel",
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

###############################################################################
# Plotting functions for the Bayesian inference
###############################################################################

"""
    plot_posterior(
    posterior::CoLocResult; 
    file::String = "posterior.png", 
    save::Bool = true
    )

This function generates diagnostic plots of a CoLocResult struct obtained from ProteinCoLoc.colocalization().

# Arguments
- `posterior`: A CoLocResult struct representing the posterior to be plotted.
- `file`: A string representing the filename to save the plot to. Must end with .png. Default is "posterior.png".
- `save`: A boolean indicating whether to save the plot to a file. Default is true.
- `fig_size`: A Vector of two Real numbers representing the width and height of the plot in centimeters. Default is [16.0,12.0].
- `dpi`: An integer representing the dots per inch (dpi) value. Default is 300.


# Returns
- `fig`: A Figure object representing the generated plot.

# Notes
This function generates five plots: four for the posterior distributions of ρ, ν, σ, and τ, and one for the posterior distribution of Δρ. Each plot includes a density plot for the control and sample. The plots are saved to a file if specified.
"""
function plot_posterior(
    posterior::CoLocResult; file::String = "posterior.png", 
    save::Bool = true, fig_size::Vector{T} = [16.0,12.0], dpi::I= 300
    ) where {T <: AbstractFloat, I <: Integer}

    fig_size_1 = 72*fig_size[1]/2.54
    fig_size_2 = 72*fig_size[2]/2.54
    fig = GLMakie.Figure(size = (fig_size_1, fig_size_2), backgroundcolor = :white, fontsize = 8, figure_padding = (1,8,1,1))

    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "ρ", ylabel = "P(ρ|data)", title = "P(ρ|data)",
        xticklabelrotation = deg2rad(60)
        )
    #
    # alpha and normalise no longer supported, instead colour can be used
    GLMakie.density!(
        ax1, posterior.posterior.μ_control,
        colormap = (:viridis, 0.3), label = "control",
        )
        
    GLMakie.density!(
        ax1, posterior.posterior.μ_sample,
        colormap = (:viridis, 0.3), label = "sample"
        )
        
    ax2 = GLMakie.Axis(fig[1, 2], xlabel = "ν", ylabel = "P(ν|data)", title = "P(ν|data)", xticklabelrotation = deg2rad(60))
    
    GLMakie.density!(
        ax2, posterior.posterior.ν_control, 
        colormap = (:viridis, 0.3), label = "ν_control"
        )

    GLMakie.density!(
        ax2, posterior.posterior.ν_sample, 
        colormap = (:viridis, 0.3), label = "ν_sample"
    )

    ax3 = GLMakie.Axis(fig[1, 3], xlabel = "σ", ylabel = "P(σ|data)", title = "P(σ|data)", xticklabelrotation = deg2rad(60))

    GLMakie.density!(
        ax3, posterior.posterior.σ_control, 
        colormap = (:viridis, 0.3), label = "σ_control"
        )

    GLMakie.density!(
        ax3, posterior.posterior.σ_sample, 
        colormap = (:viridis, 0.3), label = "σ_sample"
    )

    ax4 = GLMakie.Axis(fig[1, 4], xlabel = "τ", ylabel = "P(τ|data)", title = "P(τ|data)", xticklabelrotation = deg2rad(60))

    GLMakie.density!(
        ax4, posterior.posterior.τ_sample, 
        colormap = (:viridis, 0.3), label = "τ_sample"
        )

    GLMakie.density!(
        ax4, posterior.posterior.τ_control, 
        colormap = (:viridis, 0.3), label = "τ_control"
    )

    GLMakie.Legend(
        fig[1,5], ax1, orientation = :vertical, framevisible = false, 
        labelsize = 8, rowgap = 0, titlegap = 0,
        tellheight = false, valign = :top,
        patchsize = (5,5)
        )
    #########
    # Δρ 
    Δρ = posterior.posterior.μ_sample .- posterior.posterior.μ_control

    ax5 = GLMakie.Axis(
        fig[2, 1:5], xlabel = "Δρ", ylabel = "P(Δρ|data)", title = "P(Δρ|data)",
        limits = (-2, 2, nothing, nothing), 
        xticks = (collect(-2.0:0.1:2.0), string.(collect(-2.0:0.1:2.0))), 
        xticklabelrotation = deg2rad(60) 
        )

    GLMakie.density!(
        ax5, Δρ, label = "Δρ",
        color = (:darkgrey, 0.3)
    )

    GLMakie.vlines!(ax5, 0, color = :black, linestyle = :dash)

    # colgaps
    GLMakie.colgap!(fig.layout,1,GLMakie.Relative(0.02))
    GLMakie.colgap!(fig.layout,2,GLMakie.Relative(0.02))
    GLMakie.colgap!(fig.layout,3,GLMakie.Relative(0.02))
    GLMakie.colgap!(fig.layout,4,GLMakie.Relative(0.02))

    # colsize of legend
    GLMakie.colsize!(fig.layout,5,GLMakie.Relative(0.11))
    # rowgap between upper and lower plot
    GLMakie.rowgap!(fig.layout,1,GLMakie.Relative(0.02))

    save && GLMakie.save(file, fig , px_per_unit = dpi/72)
    return(fig)
end

"""
    bayesplot(
    prior::CoLocResult, 
    posterior::CoLocResult, 
    bf::Float64; 
    file::String = "bayesplot.png", 
    save::Bool = true,
    ρ_threshold::Float64 = 0.0
    )

This function generates a plot of the resulting prior and posterior distribution of Δρ and annotates the plot with the Bayes factor.

# Arguments
- `prior`: A CoLocResult struct representing the prior to be plotted.
- `posterior`: A CoLocResult struct representing the posterior to be plotted.
- `bf`: A Float64 representing the Bayes factor to be annotated on the plot.
- `file`: A string representing the filename to save the plot to. Must end with .png. Default is "bayesplot.png".
- `save`: A boolean indicating whether to save the plot to a file. Default is true.
- `ρ_threshold`: A Float64 representing the threshold for Δρ to be plotted as a vertical line. Default is 0.0.
- `fig_size`: A Vector of two Real numbers representing the width and height of the plot in centimeters. Default is [8.0,10.0].
- `dpi`: An integer representing the dots per inch (dpi) value. Default is 300.

# Returns
- `fig`: A Figure object representing the generated plot.

# Notes
This function calculates Δρ for the prior and posterior, creates a figure and axis, plots the density of the prior and posterior, plots a vertical line at the ρ threshold, annotates the plot with the Bayes factor, and saves the plot to a file if specified.
"""
function bayesplot(
    prior::CoLocResult, 
    posterior::CoLocResult, 
    bf::Float64; 
    file::String = "bayesplot.png", 
    save::Bool = true,
    ρ_threshold::Float64 = 0.0,
    fig_size::Vector{T} = [8.0,10.0],
    dpi::I= 300
    ) where {T <: AbstractFloat, I <: Integer}

    Δρ_post = posterior.posterior.μ_sample .- posterior.posterior.μ_control
    Δρ_prior = prior.posterior.μ_sample .- prior.posterior.μ_control

    fig_size_1 = 72*fig_size[1]/2.54
    fig_size_2 = 72*fig_size[2]/2.54
    fig = GLMakie.Figure(
        size = (fig_size_1, fig_size_2), backgroundcolor = :white, 
        fontsize = 8, figure_padding = (1,8,1,1)  
        )

    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "Δρ", ylabel = "PDF", title = "P(Δρ|data)",
        xticks = (collect(-2.0:0.4:2.0), string.(collect(-2.0:0.4:2.0))), 
        xticklabelrotation = deg2rad(60),
        limits = (-2, 2, nothing, nothing)
    )

    hist1a = GLMakie.density!(
        ax1, Δρ_prior, colormap = (:viridis, 0.3),
        label = "prior",
    )

    hist1b = GLMakie.density!(
        ax1, Δρ_post,
        label = "posterior", colormap = (:viridis, 0.3)
    )

    GLMakie.vlines!(
        ax1, ρ_threshold, color = :black, 
        linestyle = :dash, label = "Δρ = $ρ_threshold"
        )

    t = "BF = $(round(bf; digits = 4))"
    # Add BF to the plot
    GLMakie.text!(-1.9,0.05, text = t)
    GLMakie.Legend(
        fig[2,1], ax1, orientation = :horizontal, framevisible = false, 
        labelsize = 8, rowgap = 0, titlegap = 0,
        tellheight = false, valign = :top,
        patchsize = (5,5)
        )

    GLMakie.rowsize!(fig.layout,1,GLMakie.Relative(0.90))
    GLMakie.rowgap!(fig.layout,1,GLMakie.Relative(0.01))
    save && GLMakie.save(file, fig, px_per_unit = dpi/72)
    return fig
end

"""
    bayes_rangeplot(
    prior::CoLocResult, 
    posterior::CoLocResult; 
    Δρ::Vector{T} = collect(range(-0.5,0.5;step =0.05)), 
    save::Bool = true, 
    file::String = "bayes_rangeplot.png"
    ) where {T <: Real}

This function generates a plot of the Bayes factor for a range of Δρ thresholds.

# Arguments
- `prior`: A CoLocResult struct representing the prior to be plotted.
- `posterior`: A CoLocResult struct representing the posterior to be plotted.
- `Δρ`: A Vector of Real numbers representing the range of Δρ thresholds to be used. Default is collect(range(-0.5,0.5;step =0.05)).
- `save`: A boolean indicating whether to save the plot to a file. Default is true.
- `file`: A string representing the filename to save the plot to. Must end with .png. Default is "bayes_rangeplot.png".

# Returns
- `fig`: A Figure object representing the generated plot.

# Notes
This function calculates the Bayes factor for each Δρ threshold, creates a figure and axis, plots the Bayes factor for each Δρ threshold, plots a horizontal line at BF = 1, and saves the plot to a file if specified.
"""
function bayes_rangeplot(
    prior::CoLocResult,
    posterior::CoLocResult;
    Δρ::Vector{T} = collect(range(-0.5,0.5;step =0.05)),
    save::Bool = true,
    file::String = "bayes_rangeplot.png",
    fig_size::Vector{T} = [8.0,10.0],
    dpi::I= 300
    ) where {T <: AbstractFloat, I <: Integer}

    # calculate the bayes factor for each Δρ threshold
    bf = fill(0.0, length(Δρ))
    for idx in eachindex(Δρ)
        a, _, _ = compute_BayesFactor(posterior, prior; ρ_threshold = Δρ[idx])
        a > 0 ? bf[idx] = a : bf[idx] = NaN
    end 

    # return if all values are NaN or Inf
    if all(isnan.(bf))
        @warn "All values are bayes factors are zero. No plot is generated."
        return nothing
    end

    if all(isinf.(bf))
        @warn "All values are bayes factors are infinite. No plot is generated."
        return nothing
    end

    # plot the results
    ticks = collect(range(Δρ[1], Δρ[end], length = 11))
    fig_size_1 = 72*fig_size[1]/2.54
    fig_size_2 = 72*fig_size[2]/2.54

    fig = GLMakie.Figure(
        size = (fig_size_1, fig_size_2), backgroundcolor = :white,
        fontsize = 8, figure_padding = (1,8,1,1)
        )

    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "Δρ0", 
        ylabel = "log10(BF[Δρ > Δρ0 : Δρ ≤ Δρ0])", 
        title = "Bayes factor vs Δρ",
        limits = (Δρ[1], Δρ[end], nothing, nothing), 
        xticks = (ticks, string.(ticks)),
        xticklabelrotation = deg2rad(60)
    )

    GLMakie.lines!(ax1, Δρ, log10.(bf))

    GLMakie.hlines!(ax1, [0], color = :black, linestyle = :dash, label = "BF = 1")
    save && GLMakie.save(file, fig, px_per_unit = dpi/72)
    return fig
end
