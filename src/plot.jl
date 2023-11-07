## plot.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
using GLMakie

"""
    minmax_norm!(img::Matrix{Float64})
    Normalize the image to the range [0, 1].
"""
minmax_norm!(img::Matrix{Float64}) = (img .- minimum(img)) ./ (maximum(img) - minimum(img))


"""
    cm_to_px(cm, dpi)
    Convert centimeters to pixels with a given dpi to determine the size of the plot.
"""
function cm_to_px(cm, dpi)
    inch = cm / 2.54  # Convert cm to inches
    px = round(Int, inch * dpi)  # Convert inches to pixels
    return px
end

# Function to calculate font size
function calculate_font_size(resolution, scale_factor)
    return round(Int, min(resolution...) * scale_factor)
end

"""
    plot(
    img::LoadImages.MultiChannelImage,
    num_patches::T,
    cor_channel::Vector{T} = [2, 3];
    scale_channels::Bool = true,
    file::String = "test.png";
    channel_for_plot::Vector{T} = [1, 2, 3],
    save_to_file::Bool = true
    ) where {T <: Int}

    Plot the image with the patches and the correlation of each patch.

    # Arguments
    - `img::LoadImages.MultiChannelImage`: The image to plot.
    - `num_patches::Int64`: The number of patches to use.
    - `cor_channel::Vector{Int64}`: The channels to use for the calculation of correlation.
    - `scale_channels::Bool`: Whether to scale the channels intensities to the range [0, 1].
    - `file::String`: The file to save the plot to. Must end with .png or .svg.
    - `channel_for_plot::Vector{Int64}`: The channels to use for plotting. Only 1-3 channels are supported.
        The first channel is plotted in blue, the second in green and the third in red.
    - `save_to_file::Bool`: Whether to save the plot to the file.

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
    img::MultiChannelImage,
    num_patches::T,
    cor_channel::Vector{T} = [2, 3];
    scale_channels::Bool = true,
    file::String = "test.png",
    channel_for_plot::Vector{T} = [1, 2, 3],
    save_to_file::Bool = true
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
    if save_to_file == true
        savefig(plt, file)
    end
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
function plot_mask(img::MultiChannelImage,file::String = "mask.png")
    mask = _calculate_mask(img)

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

    # add the channel names to the bottom corner of each channel
    for i in 1:length(img.channels)
        annotate!(
            (i-1)*size(img.data[1])[2,] + 10,
            10,
            text(img.channels[i], :yellow, :left, 40)
            )
    end
   

    # save
    savefig(plt, file)

    return plt
end

# local correlation plot
function local_correlation_plot(
    img::MultiChannelImage,
    num_patches::I,
    cor_channel::Vector{I} = [2, 3];
    channel_for_plot::Vector{I} = [1, 2, 3],
    save::Bool = true,
    file::String = "local_correlation.png"
    ) where {I <: Int}

    # check that the number of channels is 3
    length(channel_for_plot) <= 3 || @warn "This function is currently only implemented for 3 channels"

    # get patch size
    patches = patch(img.data[1], num_patches);
    patch_size = size(patches)[3:4]

    # check that the patch size is reasonable
    patch_size[1] * patch_size[2] > 10 || @warn "Patch size is too small for local correlation. A size between 10 and 100 px is recommended."
    patch_size[1] * patch_size[2] < 100 || @warn "Patch size is too big for local correlation. A size between 10 and 100 px is recommended."
    @info "Patch size is $(patch_size[1]) x $(patch_size[2]) px = $(patch_size[1] * patch_size[2]) px²"    

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
    # z = replace(z, missing => 0)

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

#########################################################################################
### Fractional overlap Plots                                                          ###
#########################################################################################
"""
    plot_fractional_overlap(
    img::MultiChannelImage,
    num_patches::I = 64,
    channels::Vector{I};
    output_folder::String = "./test_images/"
    ) where {I <: Int}


    Generates a graphical representation of the fractional overlap between two channels. 
    Requires an input image of type MultiChannelImage, a vector of the color channels too
    be analysed and a path to the output folder to save the plot. 
"""
function plot_fractional_overlap(
    img::MultiChannelImage,
    control::MultiChannelImage,
    num_patches::I,
    channels::Vector{I};
    file::String = "fractional_overlap.png",
    save::Bool = true,
    method::String = "quantile",
    quantile_level::Float64 = 0.95
    ) where I <: Int

    # calculate and apply mask
    img = _apply_mask!(img, _calculate_mask(img))
    # calculate fractional overlap
    FO_1, FO_2 = fractional_overlap(
        img, control, channels, 
        num_patches; method = method, 
        quantile_level = quantile_level
        )

    # get patch size
    patches = patch(img.data[1], num_patches);
    patch_size = size(patches)[3:4]

    # plot the local fractional overlap
    x = [i for i in 1 : (num_patches-1)] .* patch_size[2]
    y = [i for i in 1 : (num_patches-1)] .* patch_size[1]
    z_1 = [FO_1[i, j] for i in 1 : (num_patches-1), j in 1 : (num_patches-1)]
    z_2 = [FO_2[i, j] for i in 1 : (num_patches-1), j in 1 : (num_patches-1)]

    # replace missing values with 0
    z_1 = replace(z_1, missing => 0)
    z_2 = replace(z_2, missing => 0)

    # plot the fractional overlap
    GLMakie.activate!()
    fig = GLMakie.Figure(fontsize = 12, resolution = (1200,600))
    ax = GLMakie.Axis(
        fig[1, 1],
        xlabel = "x position [px]",
        ylabel = "y position [px]",
        title = "Fractional overlap: $(channels[1]) : $(channels[2])",
        yreversed = false,
        limits = (0.0, nothing, 0.0, nothing)
        )

    co = GLMakie.contourf!(
        ax, x,y,rotr90(z_1), levels = range(0,1.001,100),
        colormap = :Blues
        )

    ax_2 = GLMakie.Axis(
        fig[1, 2],
        xlabel = "x position [px]",
        ylabel = "y position [px]",
        title = "Fractional overlap: $(channels[2]) : $(channels[1])",
        yreversed = false,
        limits = (0.0, nothing, 0.0, nothing)
        )
    
    co_2 = GLMakie.contourf!(
        ax_2, x,y,rotr90(z_2), levels = range(0,1.001,100),
        colormap = :Blues
        )
        
    GLMakie.Colorbar(fig[1,3], co, tellheight = false, label = "Fractional Overlap")
    GLMakie.resize_to_layout!(fig)

    if save == true
        # save the plot
        GLMakie.save(file,fig)
    end
    return fig
end


###############################################################################
# Plotting functions for the Bayesian inference
###############################################################################

"""
    plot_posterior(posterior::CoLocResult; file::String = "posterior.png", save::Bool = true)
    Generates diagnostic plots of a CoLocResult struct obtained from ProteinCoLoc.colocalization().

    # Arguments
    - `posterior::CoLocResult`: The CoLocResult struct to plot.
    - `file::String`: The file to save the plot to. Must end with .png.
    - `save::Bool`: Whether to save the plot to the file.
"""
function plot_posterior(posterior::CoLocResult; file::String = "posterior.png", save::Bool = true)
    hist1 = Plots.histogram(
        posterior.posterior.μ_control, legend = true, label = "μ_control", 
        title = "P(ρ|data)", alpha = 0.5, normalize = :pdf
        )
    Plots.histogram!(
        hist1, posterior.posterior.μ_sample, label = "μ_sample", 
        alpha = 0.5, normalize = :pdf
        )

    hist2 = Plots.histogram(
        posterior.posterior.ν_control, label = "ν_control", 
        legend = true, title = "P(ν|data)", alpha = 0.5,
        normalize = :pdf
        )

    Plots.histogram!(
        hist2, posterior.posterior.ν_sample, 
        label = "ν_sample", alpha = 0.5, normalize = :pdf
        )

    hist3 = Plots.histogram(
        posterior.posterior.σ_control, label = "σ_control", 
        legend = true, title = "P(σ|data)", alpha = 0.5,
        normalize = :pdf
        )

    Plots.histogram!(
        hist3, posterior.posterior.σ_sample, label = "σ_sample", 
        alpha = 0.5, normalize = :pdf
        )

    Δρ = posterior.posterior.μ_sample .- posterior.posterior.μ_control

    hist4 = Plots.histogram(
        Δρ,legend = true, label = "Δρ", normalize = :pdf,
        title = "P(Δρ|data)", alpha = 0.5)

    hist5 = Plots.histogram(
        posterior.posterior.τ_sample,legend = true,
        label = "τ_sample", title = "P(τ|data)", alpha = 0.5,
        normalize = :pdf
        )

    Plots.histogram!(
        hist5, posterior.posterior.τ_control, label = "τ_control", 
        alpha = 0.5, normalize = :pdf
        )

    p = Plots.plot(
        hist1, hist2, hist3, hist5, hist4, 
        layout = (3, 2), size = (800, 800)
        )

    save && Plots.savefig(p, file)
    return(p)
end

"""
    bayesplot(
    prior::CoLocResult, 
    posterior::CoLocResult, 
    bf::Float64; 
    file::String = "bayesplot.png", 
    save::Bool = true
    )
    Generates a plot of the resulitng prior and posterior distribution
    of Δ̢ρ and annotates the plot with the Bayes factor.
"""
function bayesplot(
    prior::CoLocResult, 
    posterior::CoLocResult, 
    bf::Float64; 
    file::String = "bayesplot.png", 
    save::Bool = true,
    ρ_threshold::Float64 = 0.0
    )

    Δρ_post = posterior.posterior.μ_sample .- posterior.posterior.μ_control
    Δρ_prior = prior.posterior.μ_sample .- prior.posterior.μ_control

    hist1 = Plots.histogram(
        Δρ_prior, legend = true, label = "prior", 
        title = "P(Δρ|data)", alpha = 0.35,
        xlabel = "Δρ", ylabel = "PDF",
        normalize = :pdf
        )

    Plots.histogram!(hist1, Δρ_post, label = "posterior", alpha = 0.35, normalize = :pdf)
    Plots.vline!(hist1, [0], label = "Δρ = $ρ_threshold", color = :black, linestyle = :dash)

    # add BF to the plot
    Plots.annotate!(
        hist1, 
        [(-2.0, 0.05, Plots.text("BF[Δρ>$ρ_threshold : Δρ ≤ $ρ_threshold] = $(round(bf; digits = 2))", 10, :left, :bottom, :black))], 
        font = "Helvetica", 
        color = :black
        )

    if save
        Plots.savefig(hist1, file)
    end
    return hist1
end


function bayes_rangeplot(
    prior::CoLocResult,
    posterior::CoLocResult;
    Δ̢ρ::Vector{T} = collect(range(-0.5,0.5;step =0.05)),
    save::Bool = true,
    file::String = "bayes_rangeplot.png"
    ) where {T <: Real}

    # calculate the bayes factor for each Δρ threshold
    bf = fill(0.0, length(Δ̢ρ))
    for idx in eachindex(Δ̢ρ)
        a, _, _ = compute_BayesFactor(posterior, prior; ρ_threshold = Δ̢ρ[idx])
        a > 0 ? bf[idx] = a : bf[idx] = NaN
    end 

    # plot the results
    p1 = Plots.plot(
        Δ̢ρ, log10.(bf), 
        xlabel = "Δρ0", 
        ylabel = "log10(BF[Δρ > Δρ0 : Δρ ≤ Δρ0])", 
        title = "Bayes factor vs Δρ",
        legend = false
        )

    Plots.hline!(p1, [0], label = "BF = 1", color = :black, linestyle = :dash)
    save && Plots.savefig(p1, file)
    return p1
end
