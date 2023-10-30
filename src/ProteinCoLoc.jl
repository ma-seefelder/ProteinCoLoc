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
import Images
import KernelDensity: kde
import Plots
import PyCall
import QuadGK: quadgk

using Turing
using Turing: Variational

include("LoadImages.jl")
include("colocalization.jl")

export MultiChannelImage, MultiChannelImageStack, colocalization
export correlation, patch, compute_BayesFactor, plot_posterior, CoLocResult
end