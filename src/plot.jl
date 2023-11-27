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

    fig = GLMakie.Figure(resolution = (size(img.data[1])[2,],size(img.data[1])[1,]), backgroundcolor = :black)
    ax1 = GLMakie.Axis(fig[1, 1], aspect = GLMakie.DataAspect(), yreversed = true)
    GLMakie.image!(ax1, img_view')

    # add lines to the image to separate the patches
    for i in 0:num_patches
        GLMakie.hlines!([i*patch_size[1]], color = :white, linealpha = 0.5, linestyle = :dash, linewidth = 1)
        GLMakie.vlines!([i*patch_size[2]], color = :white, linealpha = 0.5, linestyle = :dash, linewidth = 1)
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

    # without x and y axis ticks and labels and with a black background
    fig = GLMakie.Figure(
        resolution = (0.3*size(plt)[2,],0.3*size(plt)[1,]), 
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
    _local_correlation_plot(img, channel_for_plot, num_patches, cor_channel)
    Low level function to plot the local correlation between two channels.
    should not be called directly. Use local_correlation_plot() instead.
    Return a matrix of the local correlation values.
"""
function _local_correlation_plot(img, channel_for_plot, num_patches, cor_channel)
    # check that the number of channels is 3
    length(channel_for_plot) <= 3 || @warn "This function is currently only implemented for 3 channels"

    # get patch size
    patches = patch(img.data[1], num_patches)
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

    return ρ, patch_size
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

    ρ, patch_size = _local_correlation_plot(img, channel_for_plot, num_patches, cor_channel)
    
    # check that patches with a successful correlation calculation exist
    if sum(ismissing.(ρ)) == size(ρ)[1] * size(ρ)[2]
        @warn "No patches with a successful correlation calculation exist at $num_patches patches. A different patch size is tried."
        # try a different patch
        patch_number_range = reverse(collect(10:10:num_patches))
        for pn ∈ patch_number_range
            @info "Trying $pn patches"
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
    fig = GLMakie.Figure(resolution = (1200, 800))
    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "ρ", ylabel = "P(ρ|data)", title = "P(ρ|data)",
        limits = (-1, 1, nothing, nothing), 
        xticks = (collect(-1.0:0.2:1.0), string.(collect(-1.0:0.2:1.0))), 
        xticklabelrotation = deg2rad(60)
        )

    hist1a = GLMakie.density!(
        ax1, posterior.posterior.μ_control,
        alpha = 0.5, label = "control",
        normalization = :pdf
        )
        
    hist1b = GLMakie.density!(
        ax1, posterior.posterior.μ_sample,
        alpha = 0.5, label = "sample",
        normalization = :pdf
        )
        
    ax2 = GLMakie.Axis(fig[1, 2], xlabel = "ν", ylabel = "P(ν|data)", title = "P(ν|data)")
    
    hist2a = GLMakie.density!(
        ax2, posterior.posterior.ν_control, normalization = :pdf,
        alpha = 0.5, label = "ν_control"
        )

    hist2b = GLMakie.density!(
        ax2, posterior.posterior.ν_sample, normalization = :pdf,
        alpha = 0.5, label = "ν_sample"
    )

    ax3 = GLMakie.Axis(fig[1, 3], xlabel = "σ", ylabel = "P(σ|data)", title = "P(σ|data)")

    hist3a = GLMakie.density!(
        ax3, posterior.posterior.σ_control, normalization = :pdf,
        alpha = 0.5, label = "σ_control"
        )

    hist3b = GLMakie.density!(
        ax3, posterior.posterior.σ_sample, normalization = :pdf,
        alpha = 0.5, label = "σ_sample"
    )

    ax4 = GLMakie.Axis(fig[1, 4], xlabel = "τ", ylabel = "P(τ|data)", title = "P(τ|data)")

    hist4a = GLMakie.density!(
        ax4, posterior.posterior.τ_sample, normalization = :pdf,
        alpha = 0.5, label = "τ_sample"
        )

    hist4b = GLMakie.density!(
        ax4, posterior.posterior.τ_control, normalization = :pdf,
        alpha = 0.5, label = "τ_control"
    )

    lgd = GLMakie.Legend(
        fig[2,1:4], ax1, orientation = :horizontal, 
        framevisible = false
        )
    #########
    # Δ̢ρ 
    Δρ = posterior.posterior.μ_sample .- posterior.posterior.μ_control

    ax5 = GLMakie.Axis(
        fig[3, 1:4], xlabel = "Δρ", ylabel = "P(Δρ|data)", title = "P(Δρ|data)",
        limits = (-2, 2, nothing, nothing), 
        xticks = (collect(-2.0:0.1:2.0), string.(collect(-2.0:0.1:2.0))), 
        xticklabelrotation = deg2rad(60)
        
        )

    hist5 = GLMakie.density!(
        ax5, Δρ, normalization = :pdf,
        alpha = 0.5, label = "Δρ"
    )

    save && GLMakie.save(file, fig)
    return(fig)
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

    fig = GLMakie.Figure(resolution = (600, 600))
    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "Δρ", ylabel = "PDF", title = "P(Δρ|data)",
        xticks = (collect(-2.0:0.2:2.0), string.(collect(-2.0:0.2:2.0))), 
        xticklabelrotation = deg2rad(60)
    )

    hist1a = GLMakie.density!(
        ax1, Δρ_prior, normalize = :pdf,
        label = "prior", alpha = 0.30
    )

    hist1b = GLMakie.density!(
        ax1, Δρ_post, normalize = :pdf,
        label = "posterior", alpha = 0.30
    )

    GLMakie.vlines!(
        ax1, ρ_threshold, color = :black, 
        linestyle = :dash, label = "Δρ = $ρ_threshold"
        )

    t = "BF[Δρ>$ρ_threshold : Δρ ≤ $ρ_threshold] = $(round(bf; digits = 4))"
    # Add BF to the plot
    GLMakie.text!(-2,0.05, text = t)

    lgd = GLMakie.Legend(fig[2,1], ax1, orientation = :horizontal, framevisible = false)

    save && GLMakie.save(file, fig)

    return fig
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

    ticks = collect(range(Δ̢ρ[1], Δ̢ρ[end], length = 11))
    # plot the results
    fig = GLMakie.Figure(resolution = (600, 600))
    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "Δρ0", 
        ylabel = "log10(BF[Δρ > Δρ0 : Δρ ≤ Δρ0])", 
        title = "Bayes factor vs Δρ",
        limits = (Δ̢ρ[1], Δ̢ρ[end], nothing, nothing), 
        xticks = (ticks, string.(ticks)),
        xticklabelrotation = deg2rad(60)
    )

    GLMakie.lines!(ax1, Δ̢ρ, log10.(bf))

    GLMakie.hlines!(ax1, [0], color = :black, linestyle = :dash, label = "BF = 1")
    save && GLMakie.save(file, fig)
    return fig
end
