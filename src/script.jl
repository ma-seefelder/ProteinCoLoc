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
"""
    julia_main()::Cint

This function serves as the entry point for the application when it is run as a standalone executable.

# Returns
- `0`: A Cint representing a successful execution of the application.

# Notes
This function calls the `ProteinCoLoc.gui` function, which launches the graphical user interface (GUI) of the application. The function returns 0 to indicate a successful execution of the application.
"""
function julia_main()::Cint 
    ProteinCoLoc.gui()
    return 0
end
