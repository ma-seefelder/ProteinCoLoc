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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
=#

"""
    get_images(
    path::AbstractString, 
    nchannels::Integer, 
    stack_name::AbstractString; 
    mask::Bool = true
    )

This function retrieves and loads multi-channel images from a specified path.

# Arguments
- `path`: A string representing the path to the image files.
- `nchannels`: An integer representing the number of channels in the images.
- `stack_name`: A string representing the name of the image stack.
- `mask`: A boolean indicating whether to apply a mask to the images. Default is true.

# Returns
- `MultiChannelImageStack`: A MultiChannelImageStack object representing the loaded images.

# Errors
- Throws an error if the number of images in the image folder is not a multiple of the number of channels.
- Throws an error if the number of images in the image folder is less than 1.
- Throws an error if the image path is not valid.

# Notes
This function retrieves the image files from the specified path, checks the number of images, generates channel names, groups the images, retrieves the image names, loads the images, applies a mask to the images if specified, and returns a MultiChannelImageStack object representing the loaded images.
"""
function get_images(path::S, nchannels::I, stack_name::S; mask::Bool = true) where {I<:Integer, S<:AbstractString}
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
        length(i) > 2 ? img_name = join(i[1:end-1],"_") : img_name = i[1]
        !in(img_name,image_names) ? push!(image_names,img_name) : nothing
    end

    # load images
    for i in eachindex(images) 
        image_name = image_names[i]
        channel_files = fill("", nchannels)
        for j in 1:nchannels
            file_index = findfirst(file -> occursin(image_name * "_c" * string(j), file), files)
            !isnothing(file_index) || error("The image path is not valid.")
            channel_files[j] = joinpath(path, files[file_index])
        end
        images[i] = MultiChannelImage(image_name, channel_files, channel_names)
        # mask image
        if mask
            images[i] = ProteinCoLoc._apply_mask!(images[i], ProteinCoLoc._calculate_mask(images[i]))
        end
    end
    # return files
    return MultiChannelImageStack(images, stack_name)
end 

"""
    plot_images(
    plot_type::Symbol,
    image_stack::MultiChannelImageStack, 
    num_patches::Int, 
    channel_indices::Vector{Int}, 
    base_filename::AbstractString, 
    suffix::AbstractString="",
    cor_method::Symbol = :pearson
    )

This function generates plots for a stack of multi-channel images.

# Arguments
- `plot_type`: A Symbol representing the type of plot to be generated. Can be either :local_correlation or :patched_correlation.
- `image_stack`: A MultiChannelImageStack representing the images to be plotted.
- `num_patches`: An integer representing the number of patches to be analyzed.
- `channel_indices`: A Vector of integers representing the channels to be analyzed.
- `base_filename`: A string representing the base filename for the output files.
- `suffix`: A string representing the suffix to be added to the base filename. Default is an empty string.
- `cor_method`: A Symbol representing the correlation method to be used. Can be either :pearson, :kendall, or :spearman. Default is :pearson.

# Returns
- None. The function saves the generated plots to files with the specified filenames.

# Errors
- Throws a warning if a plot could not be generated for an image.
- Throws an error if the plot type is not valid.

# Notes
This function iterates over the images in the image stack and generates a plot for each image based on the specified plot type. The plots are saved to files with filenames constructed from the base filename, the suffix, and the name of the image.
"""
function plot_images(
    plot_type::Symbol,
    image_stack::MultiChannelImageStack, 
    num_patches::Int, 
    channel_indices::Vector{Int}, 
    base_filename::AbstractString, 
    suffix::AbstractString="",
    cor_method::Symbol = :pearson
    )

    for image in image_stack
        # plot local correlation
        if plot_type == :local_correlation
            try 
                local_correlation_plot(
                    image, num_patches, channel_indices, file = "$base_filename$suffix"*"_$(image.:name).png",
                    cor_method = cor_method
                    )
            catch
                @warn "The local correlation plot for image $(image.:name) could not be generated."
            end
        elseif plot_type == :patched_correlation
            try 
                plot(
                    image, num_patches, channel_indices, file = "$base_filename$suffix"*"_$(image.:name).png",
                    cor_method = cor_method
                    )
            catch
                @warn "The patched correlation plot for image $(image.:name) could not be generated."
            end
        else
            error("The plot type is not valid.")
        end
    end
end

# helper function to plot all the plots for a given image and control image stack
"""
    generate_plots(
    images, control_images, channel_selection_two, number_patches, number_patches_loc, 
    number_iterations, number_posterior_samples, ρ_threshold, 
    ρ_range, ρ_range_step, output_folder_path, patched_correlation_plt, local_correlation_plt, 
    bayes_factor_plt, bayes_range_plt, posterior_plt, cor_method
    )

This function generates all plots for a given image and control image stack.

# Arguments
- `images`: A MultiChannelImageStack representing the sample images.
- `control_images`: A MultiChannelImageStack representing the control images.
- `channel_selection_two`: A Vector of two integers representing the channels to be analyzed.
- `number_patches`: An integer representing the number of patches to be analyzed.
- `number_patches_loc`: An integer representing the number of local patches to be analyzed.
- `number_iterations`: An integer representing the number of iterations for the ADVI algorithm.
- `number_posterior_samples`: An integer representing the number of posterior samples to be generated.
- `ρ_threshold`: A Float64 representing the threshold for the difference in correlation.
- `ρ_range`: A Vector of two Float64 representing the range of ρ values for the Bayes factor range plot.
- `ρ_range_step`: A Float64 representing the step size for the ρ values for the Bayes factor range plot.
- `output_folder_path`: A string representing the path to the output folder.
- `patched_correlation_plt`: A boolean indicating whether to generate the patched correlation plot.
- `local_correlation_plt`: A boolean indicating whether to generate the local correlation plot.
- `bayes_factor_plt`: A boolean indicating whether to generate the Bayes factor plot.
- `bayes_range_plt`: A boolean indicating whether to generate the Bayes factor range plot.
- `posterior_plt`: A boolean indicating whether to generate the posterior plot.
- `cor_method`: A Symbol representing the correlation method to be used. Can be either :pearson, :kendall, or :spearman.

# Returns
- None. The function saves the generated plots to the specified output folder.

# Errors
- Throws a warning if any of the plots could not be generated.

# Notes
This function performs a Bayesian colocalization analysis on the images and control images, and generates various plots based on the analysis. The plots include the patched correlation plot, the local correlation plot, the Bayes factor plot, the Bayes factor range plot, and the posterior plot. The plots are saved to the specified output folder.
"""
function generate_plots(
    images, control_images, channel_selection_two, number_patches, number_patches_loc, 
    number_iterations, number_posterior_samples, ρ_threshold, 
    ρ_range, ρ_range_step, output_folder_path, patched_correlation_plt, local_correlation_plt, 
    bayes_factor_plt, bayes_range_plt, posterior_plt, cor_method
    )
    # generate patch plot
    if patched_correlation_plt
        base_file = "$output_folder_path/patched_correlation_c$(channel_selection_two[1])_c$(channel_selection_two[2])"
        plot_images(
            :patched_correlation, images, number_patches, 
            channel_selection_two, base_file, "", cor_method
            )
        plot_images(
            :patched_correlation, control_images, number_patches, 
            channel_selection_two, base_file, "_control", cor_method
            )
    end

    # generate local correlation plot
    if local_correlation_plt
        base_file = "$output_folder_path/local_correlation_c$(channel_selection_two[1])_c$(channel_selection_two[2])"
        plot_images(
            :local_correlation, images, number_patches_loc, 
            channel_selection_two, base_file, "", cor_method
            )

        plot_images(
            :local_correlation, control_images, number_patches_loc, 
            channel_selection_two, base_file, "_control", cor_method
            )
    end

    # perform colocalization analysis
    prior, posterior = colocalization(
        images, control_images, channel_selection_two, 
        number_patches; iter = number_iterations, posterior_samples = number_posterior_samples,
        cor_method = cor_method
        )

    if ismissing(prior) || ismissing(posterior)
        @warn "The colocalization analysis could not be performed."   
        return missing, missing, missing
    end

    bf, _, _ = compute_BayesFactor(posterior, prior; ρ_threshold = ρ_threshold)
    if bayes_factor_plt
        base_file = joinpath(output_folder_path, "bayes_factor_c$(channel_selection_two[1])_c$(channel_selection_two[2])")
        try 
            bayesplot(prior, posterior, bf; file = base_file * ".png", ρ_threshold = ρ_threshold)
        catch
            @warn "The bayes factor plot could not be generated."
        end
    end

    if bayes_range_plt
        base_file = joinpath(output_folder_path, "bayes_factor_range_c$(channel_selection_two[1])_c$(channel_selection_two[2])")
        try 
            bayes_rangeplot(prior, posterior; file = base_file * ".png", Δρ = collect(range(ρ_range[1],ρ_range[2];step =ρ_range_step)))
        catch
            @warn "The bayes factor range plot could not be generated."
        end
    end

    if posterior_plt
        base_file = joinpath(output_folder_path, "posterior_c$(channel_selection_two[1])_c$(channel_selection_two[2])")
        try 
            plot_posterior(posterior; file = base_file * ".png")
        catch
            @warn "The posterior plot could not be generated."
        end
    end

    return prior, posterior, bf
end

"""
    combinations2(n::Int)

This function generates all possible combinations of two elements from 1 to n without repetition.

# Arguments
- `n`: An integer representing the range of numbers from which combinations are to be generated.

# Returns
- `combination`: A Vector of Vectors where each inner vector represents a unique combination of two numbers from 1 to n.

# Notes
This function uses two nested loops to generate all possible combinations of two numbers from 1 to n. It checks to ensure that the combination has not already been added to the list and that the two numbers are not the same.
"""
function combinations2(n)
    # generate all possible combinations of two elements from 1:n without repetition 
    combination = Vector{Vector{Int}}()
    for i in 1:n
        for j in 1:n
            if [j,i] ∉ combination && i != j
                push!(combination, [i,j])
            end
        end
    end
    return combination
end

"""
    compute_stats(samples::Vector{T}) where {T<:AbstractFloat}

This function computes the mean, median, and 95% credible interval for a vector of samples.

# Arguments
- `samples`: A vector of floating point numbers representing the samples.

# Returns
- `mean_value`: The mean of the samples.
- `median_value`: The median of the samples.
- `credible_interval_value`: A 2-element array representing the 2.5% and 97.5% quantiles of the samples, which form a 95% credible interval.

# Notes
This function uses the `mean`, `median`, and `quantile` functions from the Statistics module in Julia to compute the statistics.
"""
function compute_stats(samples::Vector{T}) where {T<:AbstractFloat}
    # compute mean, median, and 95% credible interval for prior and posterior
    mean_value = mean(samples)
    median_value = median(samples)
    credible_interval_value = quantile(samples, [0.025, 0.975])
    return mean_value, median_value, credible_interval_value
end

# generate a function to generate txt file with results
"""
    generate_txt(
        prior::CoLocResult, posterior::CoLocResult, 
        bf::T, channels::Vector{I}, ρ_threshold::T; 
        file::S = "result.csv"
    ) where {T<:AbstractFloat, I<:Integer, S<:AbstractString}

This function generates a text file with statistical results for prior and posterior samples.

# Arguments
- `prior`: A `CoLocResult` object representing the prior samples.
- `posterior`: A `CoLocResult` object representing the posterior samples.
- `bf`: A floating point number representing the Bayes factor.
- `channels`: A vector of integers representing the channels.
- `ρ_threshold`: A floating point number representing the ρ threshold.
- `file`: An optional string representing the file name. If not provided, defaults to "result.csv".

# Returns
- This function does not return a value.

# Notes
This function computes the mean, median, and 95% credible interval for the prior and posterior samples, as well as for the control images. It also computes these statistics for Δρ. It then writes these statistics to a text file, with each statistic separated by a "|". If the file does not exist, it is created; otherwise, the statistics are appended to the existing file.
"""
function generate_txt(
    prior::CoLocResult, posterior::CoLocResult, 
    bf::T, channels::Vector{I}, ρ_threshold::T; 
    file::S = "result.csv"
    ) where {T<:AbstractFloat, I<:Integer, S<:AbstractString}
    prior_samples = prior.posterior
    posterior_samples = posterior.posterior
    # compute mean, median, and 95% credible interval for prior and posterior
    prior_mean, prior_median, prior_credible_interval = compute_stats(prior_samples.:μ_sample)
    posterior_mean, posterior_median, posterior_credible_interval = compute_stats(posterior_samples.:μ_sample)
    # compute mean, median, and 95% credible interval for control images
    control_prior_mean, control_prior_median, control_prior_credible_interval = compute_stats(prior_samples.:μ_control)
    control_posterior_mean, control_posterior_median, control_posterior_credible_interval = compute_stats(posterior_samples.:μ_control)

    #  compute stats for Δρ
    Δρ_mean, Δρ_median, Δρ_credible_interval = compute_stats(posterior_samples.:μ_sample .- prior_samples.:μ_control)

    # generate file if it does not exist, else append to file
    if !isfile(file)
        # generate empty file
        open(file, "w") do f
            write(
                f,
                join([             
                "channel_1", "channel_2","ρ_threshold","bf",
                "mean_prior","median_prior","lower_credible_interval_prior","upper_credible_interval_prior",
                "mean_posterior","median_posterior","lower_credible_interval_posterior","upper_credible_interval_posterior",
                "mean_control_prior","median_control_prior","lower_credible_interval_control_prior","upper_credible_interval_control_prior",
                "mean_control_posterior","median_control_posterior","lower_credible_interval_control_posterior","upper_credible_interval_control_posterior",
                "mean_Δρ","median_Δρ","lower_credible_interval_Δρ","upper_credible_interval_Δρ"], "|"),"\n"
                )
        end
    end

    # write data
    open(file, "a") do f
        write(
            f,
            join(string.(round.([             
            channels[1], channels[2], ρ_threshold, bf,
            prior_mean, prior_median, prior_credible_interval[1], prior_credible_interval[2],
            posterior_mean, posterior_median, posterior_credible_interval[1], posterior_credible_interval[2],
            control_prior_mean, control_prior_median, control_prior_credible_interval[1], control_prior_credible_interval[2],
            control_posterior_mean, control_posterior_median, control_posterior_credible_interval[1], control_posterior_credible_interval[2],
            Δρ_mean, Δρ_median, Δρ_credible_interval[1], Δρ_credible_interval[2]], digits =6), "|")),"\n"
            )
    end
    return nothing
end
