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

module ProteinCoLoc

# Dependencies
import Base: getindex, iterate
import CSV
import DataFrames
import DataFrames: DataFrame
import Distributions: pdf
import GLMakie
import Images
import KernelDensity: kde
import QuadGK: quadgk
import Statistics: quantile
import Base: Cint

using Turing
using Turing: Variational
using Mousetrap


include("LoadImages.jl")
include("colocalization.jl")
include("bayes.jl")
include("plot.jl")
include("utils.jl")
include("main.jl")
include("gui_css.jl")
include("gui.jl")
include("script.jl")

export MultiChannelImage, MultiChannelImageStack, colocalization
export correlation, patch, compute_BayesFactor, plot_posterior, CoLocResult
export fractional_overlap, plot, plot_fractional_overlap, local_correlation_plot
export plot_mask, bayesplot, bayes_rangeplot, bayesfactor_robustness
end


