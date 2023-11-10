## gui_css.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder

# This file contains the CSS for the GUI.

#####################################################################################

# Title
add_css!(
    """
    .titlebody{
        font-size: 40px;
        font-weight: bold;
        font-family: sans-serif;
        margin-bottom: 10px;
        color: #4472C4;     
    }
    """
    )
   
    
# Subtitle
add_css!(
    """
    .sub_title{
        font-size: 24px;
        font-weight: bold;
        font-family: sans-serif;
        margin-bottom: 10px;
        color: #4472C4;      
    }
    """
    )

# Text
add_css!(
    """
    .text{
        font-size: 16px;
        font-family: sans-serif;
        margin-bottom: 10px;   
    }
    """
    )

# Log
add_css!(
    """
    .log{
        font-size: 12px;
        font-family: sans-serif;
        margin-bottom: 12px; 
        color: #000000; 
    }
    """
    )

# log frame
add_css!(
    """
    .log_frame{
        border: 1px solid #E7E6E6;
        padding: 10px;
        margin-bottom: 10px;
        background-color: #E7E6E6;
    }
    """
    )
