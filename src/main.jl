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

# This file contains the main API functions for the ProteinCoLoc.jl package.
# The functions are exported and can be used by the user.
# Additionally, the function is called by the GUI to perform the analysis.

function start_analysis(
    image_path::S, # path to the images
    control_image_path::S, # path to the control images
    output_folder_path::S, # path to the output folder
    number_patches::I, # number of patchesm
    number_patches_loc::I, # number of patches for local correlation
    number_patches_bfrobustness::Vector{I}, # number of patches for bayes factor robustness
    number_channels::I, # number of channels
    channel_selection::Bool, # channel selection
    channel_selection_two::Vector{I}, # channel selection two
    patched_correlation_plt::Bool, # patched correlation plot
    local_correlation_plt::Bool, # local correlation plot
    fractional_overlap_plt::Bool, # fractional overlap plot
    bayes_factor_plt::Bool, # bayes factor plot
    bayes_range_plt::Bool, # bayes range plot
    bayes_factor_robustness_plt::Bool, # bayes factor robustness plot,
    posterior_plt::Bool, # posterior plot
    mask_plt::Bool, # mask plot
    number_iterations::I = 1000, # number of iterations
    number_posterior_samples::I = 100_000,# number of posterior samples
    ρ_threshold::Float64 = 0.1, # threshold for the bayes factor
    ρ_range::Vector{Float64} = [-0.8, 0.8], # range for the bayes factor range plot
    ρ_range_step::Float64 = 0.01 # step size for the bayes factor range plot   
    ) where {S<:AbstractString, I<:Integer}

    ###########################################################################
    # check input parameters
    ###########################################################################

    # check if the number of channels is valid
    number_channels < 2 && error("The number of channels must be at least 2.")
    # check if the image path is valid
    !isdir(image_path) && error("The image path is not valid.")
    # check if the control image path is valid
    !isdir(control_image_path) && error("The control image path is not valid.")
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
    control_images = get_images(control_image_path, number_channels, "control images")

    ###########################################################################
    # perform analysis and plotting
    ###########################################################################
    if !channel_selection
        generate_plots(
            images, control_images, channel_selection_two, number_patches, number_patches_loc, 
            number_patches_bfrobustness, number_iterations, number_posterior_samples, ρ_threshold, 
            ρ_range, ρ_range_step, output_folder_path, patched_correlation_plt, local_correlation_plt, 
            fractional_overlap_plt, bayes_factor_plt, bayes_range_plt, bayes_factor_robustness_plt, 
            posterior_plt
            )
    else
        # extract all possible combinations of channels 
        channel_combinations = combinations2(number_channels)
        # iterate over all channel combinations
        for channels ∈ channel_combinations
            generate_plots(
                images, control_images, channels, number_patches, number_patches_loc, 
                number_patches_bfrobustness, number_iterations, number_posterior_samples, ρ_threshold, 
                ρ_range, ρ_range_step, output_folder_path, patched_correlation_plt, local_correlation_plt, 
                fractional_overlap_plt, bayes_factor_plt, bayes_range_plt, bayes_factor_robustness_plt, 
                posterior_plt
            )
        end
    end

    # mask plot
    if mask_plt
        for img_set in [images, control_images]
            for img in img_set
                plot_mask(img, output_folder_path * "/mask" * string(img.:name)  * ".png")
            end
        end
    end
end