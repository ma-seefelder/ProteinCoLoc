## gui.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder

#import .ProteinCoLoc
using Mousetrap

main() do app::Application
    ########################################################################
    #### generate window and main box ######################################
    ########################################################################
    win = Window(app)
    set_title!(win, "ProteinCoLoc.jl")
    set_size_request!(win, Vector2f(800, 600))
    
    main_box = Box(ORIENTATION_VERTICAL)
    set_margin!(main_box, 50)

    ########################################################################
    #### generate box for the introduction of the program ##################
    ########################################################################
    intro_box = Box(ORIENTATION_VERTICAL)
    set_horizontal_alignment!(intro_box, ALIGNMENT_START)
    set_expand!(intro_box, true)

    ########################################################################
    # header label
    header = Label(
        "
        <span weight='heavy' color='#4472C4' size = '300%'> ProteinCoLoc.jl </span>
        <span weight='bold'> A program to analyze the colocalization of proteins in cells </span>
        <b> Author </b>: Dr. rer. nat. Manuel Seefelder
        <b> Email </b>: proteincoloc@protonmail.de
        <b> Version </b>: 0.1.0
        ")

    # adjust display of the header label
    set_justify_mode!(header, JUSTIFY_MODE_FILL)
    set_horizontal_alignment!(header, ALIGNMENT_START)

    # add text with description of the program to the intro box
    text= Label(
        "This program is used to analyze the colocalization of proteins in cells. Therfore, the user has to load several images as well as control images (e.g., secondary-antibody only control). The program will then calculate the colocalization of the proteins and display the results in several figures. Additionally, a Bayes Factor as well as the posterior probability of colocalization ρ(images) - ρ(reference) > 0 will be calculated.",
        )
    set_justify_mode!(text, JUSTIFY_MODE_FILL)
    set_wrap_mode!(text, LABEL_WRAP_MODE_WORD_OR_CHAR)
    
    # add header and text to the intro box
    push_back!(intro_box, header)
    push_back!(intro_box, text)
    push_back!(intro_box, Separator(ORIENTATION_HORIZONTAL))
    push_back!(main_box, intro_box)

    ########################################################################
    #### generate box for the input of the images and other seetings    ####
    ########################################################################
    setting_box = Box(ORIENTATION_VERTICAL)

    ##############################
    ## slider: number of patches
    ##############################
    path_slider_label = Label("<b>Number of NxN patches to be used for the analysis</b>")
    patch_slider = Scale(1,100,1)
    set_value!(patch_slider, 16)
    set_tooltip_text!(
        patch_slider, 
        "<b> Number of NxN patches to be used for the analysis </b>: The number of patches to be used for the the calculation of the patchwise correlation coefficients. The higher the number, the longer the anaylsis will take. The default value is 16."
        )
    set_should_draw_value!(patch_slider, true)
    slider_box = vbox(path_slider_label, patch_slider) # create box for the slider
    push_back!(setting_box, slider_box) # add box to the setting box

    ##############################
    ## file chooser: images
    ##############################
    image_selector_label = Label("<b> Select the folder containing the sample images </b> \n")
    #=
    #! the FileChooser is not working yet with FILE_CHOOSER_ACTION_SELECT_FOLDER #
    #! reported an bug to the developers of Mousetrap.jl                         #
    image_selector = FileChooser(FILE_CHOOSER_ACTION_SELECT_FOLDER, "Select")
    # register callbacks for the file chooser
    #on_accept!(image_selector) do self::FileChooser, files::Vector{FileDescriptor}
    #    println("User chose files at $files")
    #end
    #on_cancel!(image_selector) do self::FileChooser
    #    println("User cancelled the dialog")
    #end

    present!(image_selector)
    =#

    # create a text field in which the path to the folder can be entered
    image_selector = Entry()
    image_selector_button = Button()
    image_selector_selected_folder = Label("")
    set_child!(image_selector_button, Label("Enter"))
    set_text!(image_selector, "Folder path")
    set_size_request!(image_selector, Vector2f(500, 50))
    
    # save the text of the text field in a variable when the button is pressed
    connect_signal_clicked!(image_selector_button) do self::Button
        #IMAGE_PATH = get_text(image_selector)
        set_text!(
            image_selector_selected_folder, 
            "<b>The currently selected folder is: </b>
            $(get_text(image_selector))"
            )
        return nothing
    end
    
    image_selector_box = vbox(
        image_selector_label, 
        hbox(image_selector, image_selector_button), 
        image_selector_selected_folder
        )

    push_back!(setting_box, image_selector_box)

    ########################################################################
    #### generate box for the image display                             ####
    ########################################################################
    image_box = Box(ORIENTATION_VERTICAL)
    image_display = ImageDisplay()
    set_size_request!(image_display, Vector2f(800, 600))
    create_from_file!(image_display, "src/dummy_image.png")
    push_back!(image_box, image_display)

    ########################################################################
    #### merge setting and image box and add it to main_box             ####
    ########################################################################
    body_setting_image_box = hbox(setting_box, image_box)
    push_back!(main_box, body_setting_image_box)

    ########################################################################
    set_child!(win, main_box)
    present!(win)
end



#=

    main = GtkBox(:v)
    ########################################################################

    
    
    

    ########################################################################
    #### GUI body ##########################################################
    ########################################################################
    body = GtkBox(:h)

    ########################################################################
    #### generate box for the input of the images and the control images ###
    ########################################################################	
    # generate box 
    input_box = GtkBox(:v)

    ####### Folder selector for the images #################################
    # add label to the input box
    image_selector_label = GtkLabel("", margin = 10)
    GAccessor.markup(image_selector_label, """<b> Select the folder containing the sample images </b>""")
    GAccessor.line_wrap(image_selector_label,true)

    # file chooser 
    image_selector_button = GtkButton("", margin = 10)
    image_selector = GtkFileChooserDialog(
        "Select the folder containing the sample images", 
        image_selector_button, GtkFileChooserAction.SELECT_FOLDER
        )

    # add label and file chooser to the input box
    push!(input_box, image_selector_label)
    push!(input_box, image_selector)

    ####### number of patches to be used for the analysis ##################

    # make slider for the number of patches to be used for the analysis
    patch_slider = GtkScale(false,1:100, margin = 10)
    # set default value
    GAccessor.value(patch_slider, 16)
    # add label to the slider
    patch_slider_label = GtkLabel("Number of NxN patches to be used for the analysis", margin = 10)
    GAccessor.markup(
        patch_slider_label,
        """
        <b> Number of NxN patches to be used for the analysis </b>: The number of patches to be used for the the calculation of the patchwise correlation coefficients. The higher the number, the longer the anaylsis will take. The default value is 16.
        """
    )
    GAccessor.line_wrap(patch_slider_label,true)
    # add slider and label to the input box
    push!(input_box, patch_slider_label)
    push!(input_box, patch_slider)

    # add box to the body
    push!(body, input_box)

    ########################################################################
    #### generate box to display images                                  ###
    ########################################################################	
    # generate box 
    image_box = GtkBox(:v)
    push!(body, image_box)
    
    ########################################################################
    #### add body to the window                                          ###
    ########################################################################
    push!(main, intro_box)
    push!(main, body)
    push!(win, main)

    showall(win)
   
end

gui()
=#