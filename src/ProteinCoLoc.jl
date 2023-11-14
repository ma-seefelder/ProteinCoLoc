## ProteinCoLoc.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder

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
import Plots
import Plots: annotate!, hline!, plot, savefig, text, vline!
import PyCall
import QuadGK: quadgk
import Statistics: quantile
import StatsPlots

using Turing
using Turing: Variational
using Mousetrap


include("LoadImages.jl")
include("colocalization.jl")
include("bayes.jl")
include("plot.jl")
include("utils.jl")
include("main.jl")
include("gui.jl")

export MultiChannelImage, MultiChannelImageStack, colocalization
export correlation, patch, compute_BayesFactor, plot_posterior, CoLocResult
export fractional_overlap, plot, plot_fractional_overlap, local_correlation_plot
export plot_mask, bayesplot, bayes_rangeplot, bayesfactor_robustness
end


