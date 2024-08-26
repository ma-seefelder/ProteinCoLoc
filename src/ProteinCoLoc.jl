#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
E-Mail: manuel.seefelder@uni-ulm.de
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
 any later version.

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
import StatsBase: corspearman, corkendall
import CSV
import DataFrames
import DataFrames: DataFrame
import Distributions: pdf
import GLMakie
import Images
import KernelDensity: kde
import QuadGK: quadgk
import Statistics: quantile, mean, median
import Base: Cint
import Random: shuffle!, randperm
using Turing
using Turing: Variational
#using Mousetrap

include("load_images.jl")
include("colocalization.jl")
include("bayes.jl")
include("plot.jl")
include("utils.jl")
include("main.jl")
#include("gui_css.jl")
#include("gui.jl")
include("script.jl")

export MultiChannelImage, MultiChannelImageStack, colocalization
export correlation, patch, compute_BayesFactor, plot_posterior, CoLocResult
export plot, local_correlation_plot, plot_mask, bayesplot, bayes_rangeplot
export start_analysis, load_tiff, get_samples, apply_mask!
end

 