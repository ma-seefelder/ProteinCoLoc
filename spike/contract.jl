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

# spike/contract.jl --- the include() coupling boundary (Phase-1 D-01 fallback,
# NOTES §4) + the SIM-03 summary-contract helpers.
#
# This is the FROZEN summary contract the Wave-2 simulator targets: a synthetic
# 2-channel image flows through the UNCHANGED src/ patch()/correlation() to a
# fixed 8x8 per-patch correlation matrix and an induced mean mu. No physics here
# -- it only proves the contract end-to-end so the simulator builds against a
# known-good target rather than a moving one.
#
# DECOUPLING (hard constraint, CLAUDE.md): the two include()s below are READ-ONLY.
# src/ and the root manifests stay byte-identical to baseline f581d95; this file
# and all new code live entirely under spike/.
#
# ORDER MATTERS (NOTES §4, runtests.jl evidence): StatsBase + Statistics MUST be
# in scope BEFORE including src/colocalization.jl. correlation() builds
# `Dict(:pearson=>cor, :spearman=>corspearman, :kendall=>corkendall)` at CALL
# time, so even :pearson needs corspearman/corkendall (StatsBase) resolvable, and
# cor/mean/median/quantile come from Statistics. Images must be in scope before
# src/LoadImages.jl (its convenience constructor calls Images.otsu_threshold).

using StatsBase      # corspearman, corkendall (correlation()'s method Dict)
using Statistics     # cor, mean, median, quantile
using Images         # otsu_threshold, imfilter, Kernel (also re-exports ImageFiltering)

# --- include() coupling boundary -- READ-ONLY, never edit src/ -----------------
include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))     # MultiChannelImage(Stack)
include(joinpath(@__DIR__, "..", "src", "colocalization.jl")) # patch, correlation, _exclude_zero

"""
    build_mci(data::Vector{Matrix{Float64}}; name="sim") -> MultiChannelImage

Build a `MultiChannelImage` from a 2-element vector of channel matrices, mirroring
the convenience constructor at `src/LoadImages.jl:444-454`.

Per **D-11**, `pixel_size = size(data[1])` is the (W,H) **pixel-dimension tuple**
(`Tuple{Int,Int}`, as the struct's `I<:Int` constraint requires) -- NOT micrometers.
`otsu_threshold = Images.otsu_threshold.(data)` mirrors the source exactly. The
struct constructor asserts `length(data) == length(channels) == length(otsu_threshold)`.
"""
function build_mci(data::Vector{Matrix{Float64}}; name::String = "sim")
    @assert length(data) == 2 "contract expects exactly 2 channels, got $(length(data))"
    return MultiChannelImage(
        data,
        ["ch1", "ch2"],
        name,
        ["", ""],
        size(data[1]),                 # D-11: pixel-dim tuple, not micrometers
        Images.otsu_threshold.(data),  # mirrors src/LoadImages.jl:445
    )
end

"""
    patch_summary(mci::MultiChannelImage) -> Matrix{Union{Float64,Missing}}

The SIM-03 summary statistic: the fixed **8x8** per-patch Pearson correlation
matrix (D-10), computed by the UNCHANGED `src/colocalization.jl` `patch()` and
`correlation()`. Per patch, `_exclude_zero` drops 0.0/NaN/missing and a patch with
<=15 survivors becomes `missing` (the two contract traps).

A plain module-local function (CR-01): do NOT extend `Base.summary`, whose
contract is `(io, x)::String` for display dispatch -- returning a Matrix there
silently poisons Julia's `show`/REPL machinery and Test.jl failure formatting.
`patch_summary` gives the plan's call-site syntax without piggybacking on Base.
"""
function patch_summary(mci::MultiChannelImage)
    x = mci.data[1]
    y = mci.data[2]
    xp, yp = patch.([x, y], 8)                       # src/colocalization.jl:37, UNCHANGED
    return correlation(xp, yp; method = :pearson)    # src/colocalization.jl:221, UNCHANGED
end

"""
    induced_mu(mci::MultiChannelImage) -> Float64

The induced summary mean: `mean(skipmissing(patch_summary(mci)))`. This replicates
`src/bayes.jl::_prepare_data`'s reshape/!ismissing-filter + mean WITHOUT including
`bayes.jl` (which carries Turing/GLMakie module deps absent from the lean spike
env) -- decoupling-faithful, not scope reduction (see PLAN <interfaces>).
"""
function induced_mu(mci::MultiChannelImage)
    # CR-02: a fully-missing 8×8 summary (every patch below the ≥15-px floor) makes
    # `mean(skipmissing(...))` throw on an empty collection. Return NaN instead so
    # callers (calibration/test/demo all filter on isfinite) degrade gracefully.
    vals = collect(skipmissing(patch_summary(mci)))
    isempty(vals) && return NaN
    return Statistics.mean(vals)
end
