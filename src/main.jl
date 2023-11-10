## main.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder

# This file contains the main API functions for the ProteinCoLoc.jl package.
# The functions are exported and can be used by the user.
# Additionally, the function is called by the GUI to perform the analysis.

function start_analysis(
    image_path::S, # path to the images
    control_image_path::S, # path to the control images
    output_folder_path::S, # path to the output folder
    number_patches::I, # number of patches
    number_channels::I, # number of channels
    channel_selection::Bool, # channel selection
    channel_selection_two::Vector{I}, # channel selection two
    patched_correlation_plot::Bool, # patched correlation plot
    local_correlation_plot::Bool, # local correlation plot
    fractional_overlap_plot::Bool, # fractional overlap plot
    bayes_factor_plot::Bool, # bayes factor plot
    bayes_range_plot::Bool, # bayes range plot
    bayes_factor_robustness_plot::Bool, # bayes factor robustness plot
    number_iterations::I, # number of iterations
    number_posterior_samples::I # number of posterior samples
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

    ###########################################################################
    # load and preprocess images
    ###########################################################################
    # load images
    images = get_images(image_path, number_channels, "images")
    # load control images
    control_images = get_images(control_image_path, number_channels, "control images")

    #! continue from here

end