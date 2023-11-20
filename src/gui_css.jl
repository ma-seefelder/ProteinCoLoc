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

# This file contains the CSS for the GUI.

#####################################################################################

# Title
add_css!(
    """
    .titlebody{
        font-size: 32px;
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
        font-size: 20px;
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
        font-size: 14px;
        font-family: sans-serif;
        margin-bottom: 10px;   
    }
    """
    )

# Log
add_css!(
    """
    .log{
        font-size: 10px;
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
