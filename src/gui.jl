## gui.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder
#import .ProteinCoLoc
# load gui formatting file (gui_css.jl)
include("gui_css.jl")

function gui() 
    main() do app::Application
        ########################################################################
        #### generate window and main box ######################################
        ########################################################################
        #set_current_theme!(app,THEME_DEFAULT_DARK)
        win = Window(app)
        set_title!(win, "ProteinCoLoc.jl")
        set_size_request!(win, Vector2f(1400, 1000))
    
        main_box = Mousetrap.Box(ORIENTATION_VERTICAL)
        set_margin!(main_box, 50)
        set_horizontal_alignment!(main_box, ALIGNMENT_START)

    ########################################################################
    #### generate box for the introduction of the program ##################
    ########################################################################
    intro_box = Mousetrap.Box(ORIENTATION_VERTICAL)
    set_horizontal_alignment!(intro_box, ALIGNMENT_START)
    set_expand!(intro_box, true)

    ########################################################################
    # header label
    header_1 = Mousetrap.Label("ProteinCoLoc.jl")
    add_css_class!(header_1, "titlebody")
    header_2 = Mousetrap.Label("A program to analyze the colocalization of proteins by Bayesian inference")
    add_css_class!(header_2, "sub_title")
    header_3 = Mousetrap.Label(
        "
        <b>Author</b>: Dr. rer. nat. Manuel Seefelder
        <b>Email</b>: proteincoloc@protonmail.de
        <b>Version</b>: 0.1.0
        ")

    add_css_class!(header_3, "text")
    set_horizontal_alignment!(header_3, ALIGNMENT_START)
    header = vbox(header_1, header_2, header_3)
    set_margin!(header, 5)

    # add text with description of the program to the intro box
    text= Mousetrap.Label(
        "This software is designed for analyzing colocalization in immunofluorescence images, a process where the presence of multiple proteins in the same location is examined. To use the program, the user loads several target images and corresponding control images (e.g., secondary-antibody only control). The program then calculates the degree of colocalization among the proteins (Pearson's correlation and Fractional overlap) and presents the results in various figures. Additionally, it computes a Bayes Factor and the posterior probability of colocalization, only considering areas in the images where signals were detected. This is done by automatically masking the images prior to all calculations.",
        )

    add_css_class!(text, "text")
    set_justify_mode!(text, JUSTIFY_MODE_FILL)
    set_wrap_mode!(text, LABEL_WRAP_MODE_WORD_OR_CHAR)
    
    # add header and text to the intro box
    push_back!(intro_box, header)
    push_back!(intro_box, text)

    # add separator after the intro box
    separator = Separator()
    set_margin!(separator, 20)
    set_expand!(separator, true)
    set_size_request!(separator, Vector2f(0,5))

    push_back!(intro_box, separator)
    push_back!(main_box, intro_box)
    
    ########################################################################
    #### generate box for the input of the images and other seetings    ####
    ########################################################################
    main_setting_box = CenterBox(ORIENTATION_HORIZONTAL)
    setting_box = Mousetrap.Box(ORIENTATION_VERTICAL)
    set_spacing!(setting_box, 10)

    ##############################
    ## number of patches
    ##############################
    path_slider_label = Mousetrap.Label("<b>Number of NxN patches to be used for the analysis</b>")
    add_css_class!(path_slider_label, "text")
    set_justify_mode!(path_slider_label, JUSTIFY_MODE_LEFT)
    add_css_class!(path_slider_label, "text")

    patch_slider = SpinButton(1.0,100,1.0)
    set_value!(patch_slider, 16)
    set_tooltip_text!(
        patch_slider, 
        "<b>Number of NxN patches to be used for the analysis </b>: The number of patches to be used for the the calculation of the patchwise correlation coefficients. The higher the number, the longer the anaylsis will take. The default value is 16."
        )
    slider_box = vbox(path_slider_label, patch_slider) # create box for the slider
    set_spacing!(slider_box, 10) # set spacing between label and slider
    push_back!(setting_box, slider_box) # add box to the setting box

    ##############################
    ## number of patches for local correlation plot
    ##############################
    slider_loc_label = Mousetrap.Label("<b>Number of NxN patches (local correlation plot)</b>")
    set_justify_mode!(slider_loc_label, JUSTIFY_MODE_LEFT)
    add_css_class!(slider_loc_label, "text")

    slider_loc = SpinButton(20,500,1.0)
    set_value!(slider_loc, 200)
    set_tooltip_text!(
        slider_loc, 
        "<b>Number of NxN patches to be used for the local correlation plot </b>: The number of patches for the local correlation plot should be higher than the number of patches for the patched correlation plot or the bayesian inference. The higher the number of patches, the more fine grained the local correlation plot will be. The default value is 200."
        )

    slider_loc_box = vbox(slider_loc_label, slider_loc) # create box for the slider
    set_spacing!(slider_loc_box, 10) # set spacing between label and slider
    push_back!(setting_box, slider_loc_box) # add box to the setting box

    ###############################################################
    ## number of patches for bayes factor robustness plot
    ###############################################################
    slider_bfrobustness_label = Mousetrap.Label("<b>Number of NxN patches (bayes factor robustness plot)</b>")
    set_justify_mode!(slider_bfrobustness_label, JUSTIFY_MODE_LEFT)
    add_css_class!(slider_bfrobustness_label, "text")

    slider_bfrobustness_min = SpinButton(1,500,1.0)
    slider_bf_robustness_min_label = Mousetrap.Label("Min")
    set_value!(slider_bfrobustness_min, 10)
    set_tooltip_text!(slider_bfrobustness_min, "Minimal number patches that are explored in the bayes factor robustness plot.")

    slider_bfrobustness_max = SpinButton(1,100,1.0)
    slider_bf_robustness_max_label = Mousetrap.Label("Max")
    set_value!(slider_bfrobustness_max, 50)
    set_tooltip_text!(slider_bfrobustness_max, "Maximal number patches that are explored in the bayes factor robustness plot.")

    slider_bfrobustness_step = SpinButton(1,5,1.0)
    slider_bf_robustness_step_label = Mousetrap.Label("Step size")
    set_value!(slider_bfrobustness_step, 10)
    set_tooltip_text!(slider_bfrobustness_step, "Step size of the number of patches that are explored in the bayes factor robustness plot. The higher the step size, the faster the analysis will be. CAVE: A too small step size can lead to memory issues or long runtimes.")

    slider_bfrobustness_box_child = hbox(
        vbox(slider_bf_robustness_min_label, slider_bfrobustness_min), 
        vbox(slider_bf_robustness_max_label,slider_bfrobustness_max), 
        vbox(slider_bf_robustness_step_label,slider_bfrobustness_step)
        )

    # center the child box
    set_horizontal_alignment!(slider_bfrobustness_box_child, ALIGNMENT_CENTER)

    # generate separator that is displayed after the slider_bfrobustness_box_child
    separator_bfrobustness = Separator()
    set_margin!(separator_bfrobustness, 20)
    set_expand!(separator_bfrobustness, true)
    set_size_request!(separator_bfrobustness, Vector2f(0,5))

    slider_bfrobustness_box = vbox(slider_bfrobustness_label, slider_bfrobustness_box_child, separator_bfrobustness)
    set_spacing!(slider_bfrobustness_box, 10)
    set_spacing!(slider_bfrobustness_box_child, 10)
    push_back!(setting_box, slider_bfrobustness_box) # add box to the setting box

    ##############################
    ## file chooser: images
    ##############################
    image_selector_label = Mousetrap.Label("<b>Enter the path to the folder containing the sample images </b>")
    add_css_class!(image_selector_label, "text")
    # create a text field in which the path to the folder can be entered
    image_selector = Entry()
    set_text!(image_selector, "Sample Image Folder path")
    set_size_request!(image_selector, Vector2f(200, 20))

    image_selector_box = vbox(
        image_selector_label, 
        image_selector
        )

    set_spacing!(image_selector_box, 10) # set spacing between label and slider

    push_back!(setting_box, image_selector_box)

    ###############################
    ## file chooser: control images
    ###############################
    control_image_selector_label = Mousetrap.Label("<b>Enter the path to the folder containing the control images </b>")
    add_css_class!(control_image_selector_label, "text")
    # create a text field in which the path to the folder can be entered
    control_image_selector = Entry()
    set_text!(control_image_selector, "Control Imge Folder path")
    set_size_request!(control_image_selector, Vector2f(200, 20))

    control_image_selector_box = vbox(
        control_image_selector_label, 
        control_image_selector
        )

    set_spacing!(control_image_selector_box, 10) # set spacing between label and slider
    push_back!(setting_box, control_image_selector_box)

    ###############################
    ## file chooser: output folder
    ###############################
    output_folder_selector_label = Mousetrap.Label("<b>Enter the path to the output folder </b>")
    add_css_class!(output_folder_selector_label, "text")
    # create a text field in which the path to the folder can be entered
    output_folder_selector = Entry()
    set_text!(output_folder_selector, "Output Folder path")
    set_size_request!(output_folder_selector, Vector2f(200, 20))

    # generate separator that is displayed after output_folder_selector_box
    separator_outputfolder = Separator()
    set_margin!(separator_outputfolder, 20)
    set_expand!(separator_outputfolder, true)
    set_size_request!(separator_outputfolder, Vector2f(0,5))

    output_folder_selector_box = vbox(
        output_folder_selector_label, 
        output_folder_selector,
        separator_outputfolder
        )

    set_spacing!(output_folder_selector_box, 10) # set spacing between label and slider
    push_back!(setting_box, output_folder_selector_box)
    
    ###############################
    ## Analysis settings
    ###############################
    analysis_settings_label = Mousetrap.Label("General settings")
    add_css_class!(analysis_settings_label, "sub_title")
    push_front!(setting_box, analysis_settings_label)

    ######################################
    # channel selection
    number_channels_label = Mousetrap.Label("<b>How many color channels have been recorded: </b>")
    add_css_class!(number_channels_label, "text")
    number_channels = DropDown()
    number_channels_item_2 = push_back!(number_channels, "2 channels")
    number_channels_item_3 = push_back!(number_channels, "3 channels")
    number_channels_item_4 = push_back!(number_channels, "4 channels")
    number_channels_item_5 = push_back!(number_channels, "5 channels")
    
    number_channels_box = hbox(number_channels_label,number_channels)
    set_spacing!(number_channels_box, 10) # set spacing between label and slider
    push_back!(setting_box, number_channels_box)

    ######################################
    # select channels for selection
    channel_selection_label = Mousetrap.Label("<b>Select the channels to be used for the analysis: </b>")
    add_css_class!(channel_selection_label, "text")
    # all combinations
    channel_selection_button_all = Switch()
    set_is_active!(channel_selection_button_all, true)
    set_tooltip_text!(
        channel_selection_button_all, 
        "<b>All combinations </b>: If toggled, all possible combinations of the selected channels will be used for the analysis. For example, if three channels are recorded, the combinations 1-2, 1-3, and 2-3 will be used for the analysis."
        )

    channel_selection_button_all_label = Mousetrap.Label("All combinations")
    add_css_class!(channel_selection_button_all_label, "text")

    channel_selection_button_all_box = vbox(channel_selection_button_all_label, channel_selection_button_all)
    set_spacing!(channel_selection_button_all_box, 10) # set spacing between label and slider

    # only two channels
    channel_selection_entry_two_entry = Entry()
    channel_selection_entry_two_label = Mousetrap.Label("Only two channels")
    add_css_class!(channel_selection_entry_two_label, "text")
    set_text!(channel_selection_entry_two_entry, "1,2")
    
    set_tooltip_text!(
        channel_selection_entry_two_entry, 
        "<b> Enter the ID numeric) of the two color chnannels used for the colocalization analysis. </b> This field is ignored if the button 'All combinations' is selected. The two channels should be separated by a comma. For example, if the first and the third channel should be used, enter '1,3' in the field."
        ) # set tooltip text

    channel_selection_entry_two_box = vbox(channel_selection_entry_two_label, channel_selection_entry_two_entry)
    set_spacing!(channel_selection_entry_two_box, 10)
    
    # combine all and only two channels box
    channel_selection_box = hbox(channel_selection_button_all_box, channel_selection_entry_two_box)
    set_spacing!(channel_selection_box, 10)
    push_back!(setting_box, channel_selection_box)

    set_start_child!(main_setting_box, setting_box)

    ######################################
    ## Generated Plots (switches)
    ######################################
    generated_plots_label = Mousetrap.Label("Generated plots")
    add_css_class!(generated_plots_label, "sub_title")
    
    # patched correlation plot
    patched_correlation_plot_button = Switch()
    set_is_active!(patched_correlation_plot_button, true)
    patched_correlation_plot_button_label = Mousetrap.Label("Patched correlation plot")
    patched_correlation_plot_button_box = vbox(patched_correlation_plot_button_label, patched_correlation_plot_button)
    set_spacing!(patched_correlation_plot_button_box, 10)

    # local correlation plot
    local_correlation_plot_button = Switch()
    set_is_active!(local_correlation_plot_button, true)
    local_correlation_plot_button_label = Mousetrap.Label("Local correlation plot")
    local_correlation_plot_button_box = vbox(local_correlation_plot_button_label, local_correlation_plot_button)
    set_spacing!(local_correlation_plot_button_box, 10)

    # fractional overlap plot
    fractional_overlap_plot_button = Switch()
    set_is_active!(fractional_overlap_plot_button, true)
    fractional_overlap_plot_button_label = Mousetrap.Label("Fractional overlap plot")
    fractional_overlap_plot_button_box = vbox(fractional_overlap_plot_button_label, fractional_overlap_plot_button)
    set_spacing!(fractional_overlap_plot_button_box, 10)

    # Bayes factor plot
    bayes_factor_plot_button = Switch()
    set_is_active!(bayes_factor_plot_button, true)
    bayes_factor_plot_button_label = Mousetrap.Label("Bayes factor plot")
    bayes_factor_plot_button_box = vbox(bayes_factor_plot_button_label, bayes_factor_plot_button)
    set_spacing!(bayes_factor_plot_button_box, 10)

    # Bayes range plot
    bayes_range_plot_button = Switch()
    set_is_active!(bayes_range_plot_button, true)
    bayes_range_plot_button_label = Mousetrap.Label("Bayes range plot")
    bayes_range_plot_button_box = vbox(bayes_range_plot_button_label, bayes_range_plot_button)
    set_spacing!(bayes_range_plot_button_box, 10)

    # Bayes factor robustness plot
    bayes_factor_robustness_plot_button = Switch()
    set_is_active!(bayes_factor_robustness_plot_button, false)
    bayes_factor_robustness_plot_button_label = Mousetrap.Label("Bayes factor robustness plot")
    bayes_factor_robustness_plot_button_box = vbox(bayes_factor_robustness_plot_button_label, bayes_factor_robustness_plot_button)
    set_spacing!(bayes_factor_robustness_plot_button_box, 10)

    # mask plot
    mask_plot_button = Switch()
    set_is_active!(mask_plot_button, true)
    mask_plot_button_label = Mousetrap.Label("Mask plot")
    mask_plot_button_box = vbox(mask_plot_button_label, mask_plot_button)
    set_spacing!(mask_plot_button_box, 10)

    # posterior plot
    posterior_plot_button = Switch()
    set_is_active!(posterior_plot_button, true)
    posterior_plot_button_label = Mousetrap.Label("Posterior plot")
    posterior_plot_button_box = vbox(posterior_plot_button_label, posterior_plot_button)
    set_spacing!(posterior_plot_button_box, 10)
    
    # combine all options into one box
    plot_box = Mousetrap.Box(ORIENTATION_VERTICAL)
    push_back!(plot_box, generated_plots_label)
    set_spacing!(plot_box, 10)

    generated_plots_box = FlowBox(ORIENTATION_HORIZONTAL)
    push_back!(generated_plots_box, patched_correlation_plot_button_box)
    push_back!(generated_plots_box, local_correlation_plot_button_box)
    push_back!(generated_plots_box, fractional_overlap_plot_button_box)
    push_back!(generated_plots_box, bayes_factor_plot_button_box)
    push_back!(generated_plots_box, bayes_range_plot_button_box)
    push_back!(posterior_plot_button_box, mask_plot_button_box)
    push_back!(generated_plots_box, posterior_plot_button_box)
    push_back!(generated_plots_box, bayes_factor_robustness_plot_button_box)
    push_back!(plot_box, generated_plots_box)

    set_horizontal_alignment!(generated_plots_box, ALIGNMENT_CENTER)
    # add generated plots box to setting box
    set_center_child!(main_setting_box, plot_box)
    
    ###################################
    ## Settings for Bayesian analysis
    ###################################
    bayes_box = Mousetrap.Box(ORIENTATION_VERTICAL)
    set_spacing!(bayes_box, 10)

    # heading
    bayes_heading_label = Mousetrap.Label("Settings: Bayesian Inference")
    add_css_class!(bayes_heading_label, "sub_title")
    push_front!(bayes_box, bayes_heading_label)

    # number of iterations
    number_iterations_label = Mousetrap.Label("<b>Number of iterations: </b>")
    add_css_class!(number_iterations_label, "text")
    number_iterations = Entry()
    set_text!(number_iterations, "10000")

    set_tooltip_text!(
        number_iterations, 
        "<b>Number of iterations </b>: The number of iterations of the ADVI sampler. The higher the number, the longer the anaylsis will take and the more accurate the analyiss will be. The default value is 1000."
        )

    number_iterations_box = vbox(number_iterations_label, number_iterations)
    set_spacing!(number_iterations_box, 10)
    push_back!(bayes_box, number_iterations_box)

    # number of posterior samples
    number_posterior_samples_label = Mousetrap.Label("<b>Number of prior and posterior samples: </b>")
    add_css_class!(number_posterior_samples_label, "text")
    number_posterior_samples = Entry()
    set_text!(number_posterior_samples, "100000")

    set_tooltip_text!(
        number_posterior_samples, 
        "<b>Number of posterior samples </b>: The number of prior and posterior samples to be drawn from the ADVI sampler. The higher the number the more accurate the distributions will be. The default value is 100000. A too high number of posterior samples can lead to memory issues."
        )

    number_posterior_samples_box = vbox(number_posterior_samples_label, number_posterior_samples)
    set_spacing!(number_posterior_samples_box, 10)
    push_back!(bayes_box, number_posterior_samples_box)
    set_horizontal_alignment!(bayes_box, ALIGNMENT_CENTER)
    set_end_child!(main_setting_box, bayes_box)

    #############################
    ## Δρ threshold
    #############################
    ρ_threshold_label = Mousetrap.Label("<b>Δρ threshold: </b>")
    add_css_class!(ρ_threshold_label, "text")
    ρ_threshold_entry = Entry()
    set_text!(ρ_threshold_entry, "0.1")
    set_tooltip_text!(
        ρ_threshold_entry,
        "The threshold for the difference between correlation coefficient of sample and control image that is considered meaningful.This threshold is used to compute the Bayes Factor and to generate the Bayes Factor plot. The Bayes Factor compares the two hypothesis H0: &#916;&#961; &lt;= &#961;_threshold and H1: &#916;&#961; &gt; &#961;_threshold. The default value is 0.1 as small differences in the correlation coefficient are biologically difficult to interpret."
        )

    ρ_threshold_box = vbox(ρ_threshold_label, ρ_threshold_entry)
    set_spacing!(ρ_threshold_box, 10)
    push_back!(bayes_box, ρ_threshold_box)

    ########################################################################
    #### bottom box                                                     ####
    ########################################################################
    bottom_box = Mousetrap.Box(ORIENTATION_VERTICAL)
    set_spacing!(bottom_box, 10)

    # add horizontal separator
    bottom_separator = Separator()
    set_margin!(bottom_separator, 20)
    set_expand!(bottom_separator, true)
    set_size_request!(bottom_separator, Vector2f(0,5))

    # add button to start the analysis
    start_button = Mousetrap.Button()
    start_button_label = Mousetrap.Label("Start analysis")
    add_css_class!(start_button_label, "text")
    set_child!(start_button, start_button_label)

    # add spinner
    spinner = Spinner()
    set_is_spinning!(spinner, false)

    # combine
    bottom_button_spinner = hbox(start_button, spinner)
    set_spacing!(bottom_button_spinner, 10)
    set_alignment!(bottom_button_spinner, ALIGNMENT_END)
    # add to bottom box
    push_back!(bottom_box, bottom_separator)
    push_back!(bottom_box, bottom_button_spinner)

    # add frame to display log message
    log_frame = Frame()
    add_css_class!(log_frame, "log_frame")
    log_frame_label = Mousetrap.Label("Log")
    add_css_class!(log_frame_label, "sub_title")
    set_label_widget!(log_frame, log_frame_label)

    log = Mousetrap.Label("")
    set_child!(log_frame, log)

    set_horizontal_alignment!(log, ALIGNMENT_START)

    push_back!(bottom_box, log_frame)

    ########################################################################
    #### merge major boxes into main_box                                ####
    ########################################################################
    push_back!(main_box, main_setting_box)
    push_back!(main_box, bottom_box)

    ########################################################################
    ### generate viewport                                                ###
    ########################################################################
    viewport = Viewport()
    set_child!(viewport, main_box)
    
    ########################################################################
    set_child!(win, viewport)
    present!(win)

    ########################################################################
    ### define actions                                                   ###
    ########################################################################
    # define action for the start button
    connect_signal_clicked!(start_button) do self::Mousetrap.Button
        # activate spinner
        set_is_spinning!(spinner, true)

        # generate log string
        log_string = "Starting analysis...\n"
        set_text!(log, log_string)
        add_css_class!(log, "log")
        set_justify_mode!(path_slider_label, JUSTIFY_MODE_LEFT)

        # get the paths to the images
        image_path = get_text(image_selector)
        control_image_path = get_text(control_image_selector)
        output_folder_path = get_text(output_folder_selector)

        # check if the paths are valid
        if !isdir(image_path)
            log_string = log_string * "Error: The path to the image folder is not valid. Please enter a valid path.\n"
            set_text!(log, log_string)
            set_is_spinning!(spinner, false)
            return nothing
        end

        if !isdir(control_image_path)
            log_string = log_string * "Error: The path to the control image folder is not valid. Please enter a valid path.\n"
            set_text!(log, log_string)
            set_is_spinning!(spinner, false)
            return nothing
        end

        if !isdir(output_folder_path)
            log_string = log_string * "Error: The path to the output folder is not valid. Please enter a valid path.\n"
            set_text!(log, log_string)
            set_is_spinning!(spinner, false)
            return nothing
        end

        # check that the output folder does not contain a log file, i.e. contains the results of a previous analysis
        if isfile(joinpath(output_folder_path,"log.txt"))
            log_string = log_string * "Error: The output folder already contains a log file. Please choose a different output folder or delete the results of the previous analysis.\n"
            set_text!(log, log_string)
            set_is_spinning!(spinner, false)
            return nothing
        end

        # get the number of patches
        number_patches = Int64(get_value(patch_slider))

        # get the number of patches for the local correlation plot
        number_patches_loc = Int64(get_value(slider_loc))

        # get the number of patches for the bayes factor robustness plot
        numb_patches_bfrobustness = collect(
            range(
                Int64(get_value(slider_bfrobustness_min)),
                Int64(get_value(slider_bfrobustness_max));
                step = get_value(slider_bfrobustness_step)
                )
            )
        numb_patches_bfrobustness = Int64.(numb_patches_bfrobustness)
        
        # get the number of channels
        numb_channels = get_selected(number_channels)

        if numb_channels == number_channels_item_2
            numb_channels = 2
        elseif numb_channels == number_channels_item_3
            numb_channels = 3
        elseif numb_channels == number_channels_item_4
            numb_channels = 4
        elseif numb_channels == number_channels_item_5
            numb_channels = 5
        else
            log_string = log_string * "Error: Please select the number of channels.\n"
            set_text!(log, log_string)
            set_is_spinning!(spinner, false)
            return nothing
        end

        # get the channel selection
        channel_selection = get_is_active(channel_selection_button_all)
        channel_selection_two = get_text(channel_selection_entry_two_entry)
        channel_selection_two = split(channel_selection_two, ",")
        channel_selection_two = parse.(Int64, channel_selection_two)


        # get the generated plots
        patched_correlation_plot = get_is_active(patched_correlation_plot_button)
        local_correlation_plot = get_is_active(local_correlation_plot_button)
        fractional_overlap_plot = get_is_active(fractional_overlap_plot_button)
        bayes_factor_plot = get_is_active(bayes_factor_plot_button)
        bayes_range_plot = get_is_active(bayes_range_plot_button)
        bayes_factor_robustness_plot = get_is_active(bayes_factor_robustness_plot_button)
        posterior_plot = get_is_active(posterior_plot_button)
        mask_plot = get_is_active(mask_plot_button)

        # get the number of iterations
        numb_iterations = parse(Int64,get_text(number_iterations))

        # get the number of posterior samples
        numb_posterior_samples = parse(Int64,get_text(number_posterior_samples))

        # get the threshold for the correlation coefficient
        ρ_threshold = parse(Float64,get_text(ρ_threshold_entry))

        # print the settings to the log
        log_string = log_string * "Settings: \n"*
            "Number of NxN patches: "*string(number_patches)*"\n"*
            "Number of NxN patches for local correlation plot: "*string(number_patches_loc)*"\n"*
            "Number of NxN patches for bayes factor robustness plot: " *string(numb_patches_bfrobustness)*"\n"*
            "Number of channels: "* string(numb_channels)*"\n"*
            "Channel selection: "*string(channel_selection)*"\n"*
            "Channel selection two: "*string(channel_selection_two)*"\n"*
            "Patched correlation plot: "*string(patched_correlation_plot)*"\n"*
            "Local correlation plot: "*string(local_correlation_plot)*"\n"*
            "Fractional overlap plot: "*string(fractional_overlap_plot)*"\n"*
            "Bayes factor plot: "*string(bayes_factor_plot)*"\n"*
            "Bayes range plot: "*string(bayes_range_plot)*"\n"*
            "Mask plot: "*string(mask_plot)*"\n"*
            "Posterior plot: "*string(posterior_plot)*"\n"*
            "Bayes factor robustness plot: "*string(bayes_factor_robustness_plot)*"\n"*
            "Number of iterations: "*string(numb_iterations)*"\n"*
            "Number of posterior samples: "*string(numb_posterior_samples)*"\n"*
            "Δρ threshold: "*string(ρ_threshold)*"\n"
            
        set_text!(log, log_string)

        # define the range of the correlation coefficient
        ρ_range = [-0.8,0.8]
        ρ_range_step = 0.05

        # start the analysis
        start_analysis(
            image_path, # path to the images
            control_image_path, # path to the control images
            output_folder_path, # path to the output folder
            number_patches, # number of patches
            number_patches_loc, # number of patches for local correlation plot
            numb_patches_bfrobustness, # number of patches for bayes factor robustness plot 
            numb_channels, # number of channels
            channel_selection, # channel selection
            channel_selection_two, # channel selection two
            patched_correlation_plot, # patched correlation plot
            local_correlation_plot, # local correlation plot
            fractional_overlap_plot, # fractional overlap plot
            bayes_factor_plot, # bayes factor plot
            bayes_range_plot, # bayes range plot
            bayes_factor_robustness_plot, # bayes factor robustness plot
            posterior_plot, # posterior plot 
            mask_plot, # mask plot 
            numb_iterations, # number of iterations
            numb_posterior_samples, # number of posterior samples,
            ρ_threshold, # threshold for the correlation coefficient
            ρ_range, # range of the correlation coefficient
            ρ_range_step # step size of the correlation coefficient
            )

        # deactivate spinner
        set_is_spinning!(spinner, false)

        # print finished to log
        log_string = log_string * "Finished analysis."
        set_text!(log, log_string)

        # save log to file in output folder
        open(joinpath(output_folder_path,"log.txt"), "w") do log_file
            write(log_file, log_string)
        end
        
        return nothing
    end
end
end


