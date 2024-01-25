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

# This file contains the main API functions for the ProteinCoLoc.jl package.
# The functions are exported and can be used by the user.
# Additionally, the function is called by the GUI to perform the analysis.

"""
    start_analysis(
    image_path::S, 
    control_image_path::S, 
    output_folder_path::S, 
    number_patches::I, 
    number_patches_loc::I, 
    number_channels::I, 
    channel_selection::Bool, 
    channel_selection_two::Vector{I}, 
    patched_correlation_plt::Bool, 
    local_correlation_plt::Bool, 
    bayes_factor_plt::Bool, 
    bayes_range_plt::Bool, 
    posterior_plt::Bool, 
    mask_plt::Bool, 
    number_iterations::I = 1000, 
    number_posterior_samples::I = 100_000,
    ρ_threshold::Float64 = 0.1, 
    ρ_range::Vector{Float64} = [-0.8, 0.8], 
    ρ_range_step::Float64 = 0.01, 
    shuffle::Bool = false, 
    shuffle_method::Symbol = :block,
    cor_method::Symbol = :pearson 
    ) where {S<:AbstractString, I<:Integer}

This function starts the analysis of multi-channel images.

# Arguments
- `image_path`: A string representing the path to the images.
- `control_image_path`: A string representing the path to the control images.
- `output_folder_path`: A string representing the path to the output folder.
- `number_patches`: An integer representing the number of patches.
- `number_patches_loc`: An integer representing the number of patches for local correlation.
- `number_channels`: An integer representing the number of channels.
- `channel_selection`: A boolean indicating whether to select channels.
- `channel_selection_two`: A Vector of integers representing the two selected channels.
- `patched_correlation_plt`: A boolean indicating whether to plot the patched correlation.
- `local_correlation_plt`: A boolean indicating whether to plot the local correlation.
- `bayes_factor_plt`: A boolean indicating whether to plot the Bayes factor.
- `bayes_range_plt`: A boolean indicating whether to plot the Bayes range.
- `posterior_plt`: A boolean indicating whether to plot the posterior.
- `mask_plt`: A boolean indicating whether to plot the mask.
- `number_iterations`: An integer representing the number of iterations. Default is 1000.
- `number_posterior_samples`: An integer representing the number of posterior samples. Default is 100_000.
- `ρ_threshold`: A Float64 representing the threshold for the Bayes factor. Default is 0.1.
- `ρ_range`: A Vector of Float64 representing the range for the Bayes factor range plot. Default is [-0.8, 0.8].
- `ρ_range_step`: A Float64 representing the step size for the Bayes factor range plot. Default is 0.01.
- `shuffle`: A boolean indicating whether to shuffle patches. Default is false.
- `shuffle_method`: A Symbol representing the shuffle method. Default is :block.
- `cor_method`: A Symbol representing the correlation method. Default is :pearson. Can be :pearson, :spearman, or :kendall.

# Returns
- Nothing. The function saves the plots to the specified output folder.

# Notes
This function checks the input parameters, loads and preprocesses the images, performs the analysis, and generates the plots.
"""
function start_analysis(
    image_path::S, # path to the images
    control_image_path::S, # path to the control images
    output_folder_path::S, # path to the output folder
    number_patches::I, # number of patchesm
    number_patches_loc::I, # number of patches for local correlation
    number_channels::I, # number of channels
    channel_selection::Bool, # channel selection
    channel_selection_two::Vector{I}, # channel selection two
    patched_correlation_plt::Bool, # patched correlation plot
    local_correlation_plt::Bool, # local correlation plot
    bayes_factor_plt::Bool, # bayes factor plot
    bayes_range_plt::Bool, # bayes range plot
    posterior_plt::Bool, # posterior plot
    mask_plt::Bool, # mask plot
    number_iterations::I = 1000, # number of iterations
    number_posterior_samples::I = 100_000,# number of posterior samples
    ρ_threshold::Float64 = 0.1, # threshold for the bayes factor
    ρ_range::Vector{Float64} = [-0.8, 0.8], # range for the bayes factor range plot
    ρ_range_step::Float64 = 0.01, # step size for the bayes factor range plot ,
    shuffle::Bool = false, # shuffle patches
    shuffle_method::Symbol = :block, # shuffle method
    cor_method::Symbol = :pearson # correlation method  
    ) where {S<:AbstractString, I<:Integer}

    ###########################################################################
    # check input parameters
    ###########################################################################

    # check if the number of channels is valid
    number_channels < 2 && error("The number of channels must be at least 2.")
    # check if the image path is valid
    !isdir(image_path) && error("The image path is not valid.")
    # check if the output folder path is valid
    !isdir(output_folder_path) && error("The output folder path is not valid.")
    # check if the output folder contains a log file
    isfile(joinpath(output_folder_path, "log.txt")) && error("The output folder contains a log file. Please choose another output folder or delete the log file and results from a previous analysis.")
    # check if the number of patches is valid
    number_patches < 1 && error("The number of patches must be at least 1.")
    # check that the number of selected channels is 2
    length(channel_selection_two) != 2 && error("The number of selected channels must be 2.")
    # check that ρ_range[1] < ρ_range[2]
    ρ_range[1] >= ρ_range[2] && error("ρ_range[1] must be smaller than ρ_range[2].")

    ###########################################################################
    # load and preprocess images
    ###########################################################################
    # load images
    images = get_images(image_path, number_channels, "images")
    # load control images
    if !shuffle 
        # check if the control image path is valid
        !isdir(control_image_path) && error("The control image path is not valid.")
        # load control images
        control_images = get_images(control_image_path, number_channels, "control images")
    else
        control_unshuffled = get_images(image_path, number_channels, "shuffled_control", mask = false)
        control = Vector{MultiChannelImage}(undef, control_unshuffled.num_images)
        for (img, idx) in zip(control_unshuffled, 1:control_unshuffled.num_images)
            img.:name = "shuffled_control_" * img.:name
            if shuffle_method == :pixel
                control[idx] = shuffle_pixels(img)
            elseif shuffle_method == :block
                control[idx] = shuffle_blocks(img, 9)
            else 
                @error "This shuffling is not supported. Please choose either :pixel or :block."
            end
        end
        control_images = MultiChannelImageStack(control, "shuffled_control")
        control_unshuffled = nothing # delete the original images to save RAM
        println("Shuffled control images.")
    end


    ###########################################################################
    # perform analysis and plotting
    ###########################################################################
    if !channel_selection
        prior, posterior, bf = generate_plots(
            images, control_images, channel_selection_two, number_patches, number_patches_loc, 
            number_iterations, number_posterior_samples, ρ_threshold, 
            ρ_range, ρ_range_step, output_folder_path, patched_correlation_plt, local_correlation_plt, 
            bayes_factor_plt, bayes_range_plt, posterior_plt, cor_method
            )

        if !ismissing(prior)
            generate_txt(prior, posterior, bf, channel_selection_two, ρ_threshold, file = joinpath(output_folder_path, "result.txt"))
            # write the DataFrames prior_samples and posterior_samples to csv
            CSV.write(
                joinpath(output_folder_path,"prior_samples_channel$(string(channel_selection_two[1]))_$(channel_selection_two[2]).csv"),
                prior.posterior
                )

            CSV.write(
                joinpath(output_folder_path,"posterior_samples_channel$(string(channel_selection_two[1]))_$(channel_selection_two[2]).csv"),
                posterior.posterior
                )
        end
    else
        # extract all possible combinations of channels 
        channel_combinations = combinations2(number_channels)
        # iterate over all channel combinations
        for channels ∈ channel_combinations
            prior, posterior, bf = generate_plots(
                images, control_images, channels, number_patches, number_patches_loc, 
                number_iterations, number_posterior_samples, ρ_threshold, 
                ρ_range, ρ_range_step, output_folder_path, patched_correlation_plt, local_correlation_plt, 
                bayes_factor_plt, bayes_range_plt, posterior_plt, cor_method)

            if !ismissing(prior)
                generate_txt(prior, posterior, bf, channels, ρ_threshold, file = joinpath(output_folder_path, "result.txt"))
                # write the DataFrames prior_samples and posterior_samples to csv
                CSV.write(
                    joinpath(output_folder_path,"prior_samples_channel$(string(channels[1]))_$(channels[2]).csv"),
                    prior.posterior
                )

                CSV.write(
                    joinpath(output_folder_path,"posterior_samples_channel$(string(channels[1]))_$(channels[2]).csv"),
                    posterior.posterior
                )
            end
        end
    end

    # mask plot
    if mask_plt
        for img_set in [images, control_images]
            for img in img_set
                plot_mask(img, joinpath(output_folder_path, "mask" * string(img.:name)  * ".png"))
            end
        end
    end
end