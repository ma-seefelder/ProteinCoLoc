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
import Statistics: quantile, mean, median, cor
import Base: Cint
import Random: shuffle!, randperm
#using Mousetrap

include("LoadImages.jl")
include("colocalization.jl")
include("results.jl")
# --- Amortized inference subsystem (v2.0, PROD-02) --------------------------------------
# Grid-parametric summary/encoder + single-source dimension helpers (the ONE grid-coupling
# site), then the grid-parametric data-gen/cache, then the grid-keyed estimator registry.
include("amortized/summary.jl")
include("amortized/datagen.jl")
include("registry.jl")
# Amortized READ surfaces (PROD-01): NPE posterior inference (infer.jl), NRE amortized Bayes
# factor + non-clamped KDE baseline (bf.jl), and the density+noise+posterior-predictive OOD flag
# (ood.jl). Ordered after registry.jl; ood.jl uses standardize_summary (infer.jl) and
# posterior_for, so the order is infer → bf → ood.
include("amortized/infer.jl")
include("amortized/bf.jl")
include("amortized/ood.jl")
# Amortized TRAINING layer (PROD-01/PROD-02, D-06): the NPE architecture + training, NRE ratio
# training, CPU-resident Flux.state persistence, and the single reusable per-grid pipeline.
# Ordered after ood.jl (uses build_estimator, fit_ood_nulls, standardize_summary) and after
# registry.jl (line above): pipeline.jl fills the `_train_grid_pipeline` referenced by
# `train_and_register`. Order: architecture → train_npe → train_ratio → persist → pipeline.
include("amortized/architecture.jl")
include("amortized/train_npe.jl")
include("amortized/train_ratio.jl")
include("amortized/persist.jl")
include("plot.jl")
include("utils.jl")
include("main.jl")
#include("gui_css.jl")
#include("gui.jl")
include("script.jl")

# --- Internal Turing/ADVI reference path (D-01 / D-03) ------------------------------------
# The ADVI reference implementation lives in ext/ProteinCoLocTuringExt.jl and loads ONLY when
# Turing is present (weakdep-gated). These generic-function stubs give core code
# (utils.jl driver, exported plotters) a binding to reference; the extension adds the methods.
# `colocalization` / `compute_BayesFactor` are intentionally NOT exported (breaking release):
# the amortized path is the sole supported public API. The ADVI-result plotters keep their
# exported names but only gain methods when Turing is loaded.
function colocalization end
function compute_BayesFactor end
function plot_posterior end
function bayesplot end
function bayes_rangeplot end

export AbstractMultiChannelImage, MultiChannelImage, MultiChannelImageStack
export image_data, channel_names, image_name, image_paths, pixel_dimensions, otsu_thresholds, num_channels
export correlation, patch, plot_posterior
export plot, local_correlation_plot,plot_mask, bayesplot, bayes_rangeplot
# D-02 result-type hierarchy (v2.0): abstract supertype + shipped amortized subtype +
# shared accessor interface. `AdviColocResult` (the internal Turing/ADVI path) is
# deliberately NOT exported (D-01).
export AbstractColocResult, AmortizedColocResult, OODVerdict, CalibrationMeta
export delta_rho, bayes_factor, is_ood, posterior_draws
# Grid-keyed estimator registry (PROD-02, D-04): the public num_patches-as-key entry points.
export estimator_for, register!, train_and_register
end

 