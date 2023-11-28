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

using PackageCompiler
"""
    create_app(source_dir::String, target_dir::String, script::String; force::Bool = false)

This function compiles a Julia project into a standalone application.

# Arguments
- `source_dir`: A string representing the path to the source directory of the Julia project.
- `target_dir`: A string representing the path to the target directory where the application will be created.
- `script`: A string representing the path to the script file that will be run when the application is launched.
- `force`: A boolean indicating whether to force the creation of the application if the target directory already exists. Default is false.

# Returns
- This function does not return a value.

# Notes
This function uses the PackageCompiler.jl package to compile the Julia project. The `force` argument can be used to overwrite an existing application in the target directory.
"""
create_app(
    "C:/Users/Manuel/Documents/GitHub/ProteinCoLoc", 
    "C:/Users/Manuel/Desktop/ProteinCoLoc", 
    script="C:/Users/Manuel/Documents/GitHub/ProteinCoLoc/src/script.jl";
    force = true
    )