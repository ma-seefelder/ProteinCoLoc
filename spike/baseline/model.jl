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

# spike/baseline/model.jl --- 04-02 Task 2: the read-only lift of the hierarchical
# colocalization() @model into the ISOLATED baseline env (D-01), plus the
# build_coloc_model(...) constructor the ADVI port (04-04) will drive.
#
# DECOUPLING (hard constraint, CLAUDE.md; T-04-DECOUPLE): the three include()s
# below are READ-ONLY. src/ and the root manifests stay byte-identical to baseline
# f581d95; this file and all new code live entirely under spike/. This is the same
# include()-only coupling discipline established in spike/contract.jl -- one level
# deeper, so the reach is "..","..","src" instead of "..","src".
#
# ISOLATION (D-01, the Phase-1 landmine): Turing cannot co-resolve with the pinned
# NeuralEstimators 0.2.1 in the spike env, and the parent env's Turing is broken.
# This env (spike/baseline/, its own pinned Project/Manifest) is the ONLY safe home
# for the ground-truth ADVI baseline.
#
# ORDER MATTERS (contract.jl NOTES §4): StatsBase + Statistics must be in scope
# BEFORE src/colocalization.jl (correlation() builds a method Dict over
# cor/corspearman/corkendall at call time); Images before src/LoadImages.jl (its
# convenience constructor calls Images.otsu_threshold); Turing/DataFrames/DynamicPPL
# before src/bayes.jl (the @model macro, the DataFrame-typed CoLocResult field, and
# the DynamicPPL.Model-typed convert_posterior_samples signature).
#
# SCOPE: this file only LOADS + CONSTRUCTS the model. No variational-inference
# call lives here -- that ADVI port is 04-04. Reaching the nested @model directly is
# impossible (it is local to colocalization()'s body), so the block is lifted
# VERBATIM as a top-level `coloc_model` (spike 003 MODEL-1 confirmed it captures no
# enclosing variable, so a verbatim copy is faithful).

# --- deps that must be in scope before the src/ includes --------------------------
using StatsBase      # corspearman, corkendall (correlation()'s method Dict)
using Statistics     # cor, mean, median, quantile
using Images         # otsu_threshold (LoadImages convenience constructor)
using DataFrames     # DataFrame (CoLocResult field type in bayes.jl)
using Turing         # @model, filldist, Truncated/Cauchy/Exponential/Normal/TDist, DynamicPPL
using Random         # Xoshiro (smoke_mci_pair reproducible fixture)

# bayes.jl references the bare module name `DynamicPPL` (convert_posterior_samples
# signature, line ~24). `using Turing` loads DynamicPPL but does NOT bind that name
# in Main, and `using DynamicPPL` fails because it is a transitive (not direct) dep.
# Alias it from Turing so bayes.jl includes cleanly -- no dep added, no src/ edit.
const DynamicPPL = Turing.DynamicPPL

# --- include() coupling boundary -- READ-ONLY, never edit src/ --------------------
include(joinpath(@__DIR__, "..", "..", "src", "LoadImages.jl"))     # MultiChannelImage(Stack)
include(joinpath(@__DIR__, "..", "..", "src", "colocalization.jl")) # patch, correlation, _exclude_zero
include(joinpath(@__DIR__, "..", "..", "src", "bayes.jl"))          # colocalization, CoLocResult, _prepare_data, compute_BayesFactor

# --- verbatim lift of the hierarchical @model (src/bayes.jl:267-303) --------------
# Copied byte-for-byte from the nested `@model function model(control, sample)` so
# the baseline can instantiate + (in 04-04) run ADVI on it as a top-level model.
# The likelihood is a location-scale Student-t; priors are the exact Truncated
# Cauchy/Normal + Exponential hierarchy of src/bayes.jl. DO NOT diverge from src/.
@model function coloc_model(control, sample)
    # get the number of images
    num_control = size(control, 1)
    num_sample = size(sample, 1)

    # ============= gloabal priors (per biological condition) ============= #
    # mean, degrees of freedom and standard deviation of the control
    μ_control ~ Truncated(Cauchy(0, 0.3),-1,1)
    ν_control ~ Exponential()
    σ_control ~ Truncated(Cauchy(0, 0.3),0.0001,1)
    τ_control ~ Truncated(Cauchy(0.1, 0.3),0.0001,1)
    # mean, degrees of freedom and standard deviation of the sample, and the patch heterogeneity
    μ_sample ~ Truncated(Cauchy(0, 0.3),-1,1)
    ν_sample ~ Exponential()
    σ_sample ~ Truncated(Cauchy(0.1, 0.3),0.0001,1)
    τ_sample ~ Truncated(Cauchy(0.1, 0.3),0.0001,1)

    # ============= local priors (per image) ============= #
    μ_control_image ~ filldist(Truncated(Normal(μ_control, σ_control),-1,1), num_control) # mean of the control for each patch
    ν_control_image ~ filldist(Exponential(ν_control), num_control) # degress of freedom of the control for each patch
    σ_control_image ~ filldist(Truncated(Normal(σ_control, τ_control),0,1), num_control) # standard deviation of the control for each patch

    μ_sample_image ~ filldist(Truncated(Normal(μ_sample, σ_sample),-1,1), num_sample) # mean of the sample for each patch
    ν_sample_image ~ filldist(Exponential(ν_sample), num_sample) # degress of freedom of the sample for each patch
    σ_sample_image ~ filldist(Truncated(Normal(σ_sample, τ_sample),0,1), num_sample) # standard deviation of the control for each patch

    # likelihood
    for idx ∈ 1:num_control
        control[idx] ~ TDist(ν_control_image[idx]) * σ_control_image[idx] + μ_control_image[idx]
    end

    for idx ∈ 1:num_sample
        sample[idx] ~ TDist(ν_sample_image[idx]) * σ_sample_image[idx] + μ_sample_image[idx]
    end
end

# --- guarded re-simulation hook (04-01 spike/npe/resimulate.jl) --------------------
# 04-01 (same wave) provides `resimulate_holdout` to reproduce byte-identical raw
# holdout MultiChannelImages from stored θ + index (Pitfall 4). It runs the physics
# simulator (spike/simulator/forward.jl), whose deps (ImageTransformations,
# CoordinateTransformations, Interpolations) are deliberately NOT in this env
# (Task 1 dep set / D-11 minimalism). So the include is guarded on BOTH file
# presence AND load success: when 04-01 has merged and its deps resolve here, the
# real re-simulation path lights up; otherwise the self-contained synthetic smoke
# pair below stands in. Either way model.jl always loads.
const RESIMULATE_PATH = joinpath(@__DIR__, "..", "npe", "resimulate.jl")
const RESIMULATE_AVAILABLE = Ref(false)
if isfile(RESIMULATE_PATH)
    try
        include(RESIMULATE_PATH)
        RESIMULATE_AVAILABLE[] = true
    catch err
        @warn "spike/npe/resimulate.jl present but did not load in the baseline env \
               (likely simulator deps absent by design); falling back to the synthetic \
               smoke pair." exception = (err, catch_backtrace())
    end
end

"""
    smoke_mci_pair(; imsize=(128,128), seed=0xC0FFEE)
        -> (sample::MultiChannelImage, control::MultiChannelImage)

A self-contained, dependency-light `MultiChannelImage` pair for smoke-constructing
the model WITHOUT the physics simulator (whose deps are excluded from this env).

Two shared-latent correlated non-negative channels per stack (a linear mix of a
shared base field and independent noise, small-positive background floor) so the
UNCHANGED `patch()`/`correlation()` yield non-`missing` 8×8 summaries and
`_prepare_data` returns signal (never `nothing`). This is a construction fixture,
NOT the calibrated Phase-2 forward model -- use `resimulate_holdout` (04-01) for a
faithful re-simulated holdout pair once its env deps are present.
"""
function smoke_mci_pair(; imsize::Tuple{Int,Int} = (128, 128), seed::Integer = 0xC0FFEE)
    rng = Random.Xoshiro(seed)
    function _stack(name::String)
        base = rand(rng, imsize...)
        c1 = base .+ 0.20 .* rand(rng, imsize...) .+ 0.05
        c2 = 0.60 .* base .+ 0.40 .* rand(rng, imsize...) .+ 0.05
        return MultiChannelImage(
            [c1, c2],
            ["ch1", "ch2"],
            name,
            ["", ""],
            size(c1),
            Images.otsu_threshold.([c1, c2]),
        )
    end
    return _stack("sample"), _stack("control")
end

"""
    build_coloc_model(mci_sample, mci_control, channels, num_patches) -> DynamicPPL.Model

Instantiate the lifted `coloc_model` for a sample/control `MultiChannelImage` pair.

Each image is wrapped in a single-image `MultiChannelImageStack` and pushed through
the UNCHANGED `_prepare_data` (from the read-only `src/bayes.jl` include) to produce
the array-of-correlation-vectors the model consumes -- exactly the transformation
`colocalization()` performs before building its own (nested) model, so the returned
object is API-identical to the one 04-04 will run ADVI on. No inference is run here.

`channels` is the 2-element channel selector (e.g. `[1, 2]`); `num_patches` is the
per-axis patch count (8 for the fixed 8×8 spike grid). Errors if `_prepare_data`
finds no signal above background in either stack.
"""
function build_coloc_model(
    mci_sample::MultiChannelImage,
    mci_control::MultiChannelImage,
    channels::Vector{<:Integer},
    num_patches::Integer,
)
    sample_stack  = MultiChannelImageStack([mci_sample], "sample")
    control_stack = MultiChannelImageStack([mci_control], "control")

    sample_data = _prepare_data(sample_stack, collect(Int, channels), Int(num_patches))
    ctrl_data   = _prepare_data(control_stack, collect(Int, channels), Int(num_patches))

    (sample_data === nothing || ctrl_data === nothing) && error(
        "build_coloc_model: _prepare_data returned nothing -- no signal above " *
        "background for channels=$channels, num_patches=$num_patches.",
    )

    return coloc_model(ctrl_data, sample_data)
end
