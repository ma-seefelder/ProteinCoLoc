## gui.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder

import .ProteinCoLoc
using Gtk

function gui()
    win = GtkWindow("ProteinCoLoc.jl", 400, 200)

    # generate notebook to store pages
    notebook = GtkNotebook()

    #######################################################################
    ################### Page 1: Load images ###############################
    #######################################################################
    box = GtkBox(:v)

    # create a box,button and label for tab
    tab_box = GtkBox(:h)
    button = GtkButton("Load images")
    push!(tab_box, button)
    push!(tab_box, GtkLabel("Load images"))

    # append page to notebook
    push!(notebook, box)

    #######################################################################
    ################### Page 2: Load images ###############################
    #######################################################################

    box = GtkBox(:v)
    # 

    #######################################################################
    # add notebook to window
    push!(win, notebook)

    # display window
    showall(win)
end