#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder

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

#############################################################################################
# ProteinCoLocTuringExt — the internal Turing/ADVI reference path (D-01 / D-03).
#
# Weakdep-gated package extension: this whole file loads ONLY when `Turing` is present in the
# environment. It provides the (unexported, internal-only) ADVI reference implementation used
# by the D-05 ship-gate as the non-amortized accuracy baseline — NOT as supported public API.
# The public v2.0 surface is the amortized NPE/NRE path in `src/`.
#
# It defines `AdviColocResult <: ProteinCoLoc.AbstractColocResult` (the demoted `CoLocResult`,
# fields verbatim) and adds methods to the core generic-function stubs `colocalization`,
# `compute_BayesFactor`, `plot_posterior`, `bayesplot`, `bayes_rangeplot` plus the D-02
# accessor interface. Code moved verbatim from the former `src/bayes.jl` + the ADVI-result
# plotters formerly in `src/plot.jl`.
#############################################################################################
module ProteinCoLocTuringExt

using ProteinCoLoc
using ProteinCoLoc: MultiChannelImageStack, correlation, patch, AbstractColocResult
import ProteinCoLoc: colocalization, compute_BayesFactor, plot_posterior, bayesplot,
    bayes_rangeplot, delta_rho, bayes_factor, is_ood, posterior_draws

using Turing
using Turing: Variational
import DataFrames
import DataFrames: DataFrame
import KernelDensity: kde
import QuadGK: quadgk
import Distributions: pdf
import GLMakie

######################################################################
# AdviColocResult — the demoted CoLocResult (internal Turing result)
######################################################################
"""
    struct AdviColocResult <: ProteinCoLoc.AbstractColocResult

Internal result of the Turing/ADVI reference colocalization path (formerly `CoLocResult`).
Retained as internal reference code for the D-05 ship-gate; NOT exported (D-01). Fields are
verbatim from the pre-v2.0 `CoLocResult`.

# Fields
- `img`: sample `MultiChannelImageStack`.
- `control`: control `MultiChannelImageStack`.
- `channels`: analyzed channel indices.
- `num_patches`: patch grid used.
- `posterior`: posterior/prior parameter `DataFrame`.
- `advi_result`: the raw ADVI/chain object.
"""
struct AdviColocResult <: ProteinCoLoc.AbstractColocResult
    img::MultiChannelImageStack
    control::MultiChannelImageStack
    channels::Vector{Int64}
    num_patches::Int64
    posterior::DataFrame
    advi_result
end

# --- D-02 accessor interface for the internal ADVI result ---------------------------------
"""
    delta_rho(r::AdviColocResult)

Δρ posterior draws for the ADVI path: `μ_sample − μ_control` (the same columns the KDE
Bayes factor consumes).
"""
delta_rho(r::AdviColocResult) = r.posterior.μ_sample .- r.posterior.μ_control

posterior_draws(r::AdviColocResult) = r.posterior

# The ADVI path carries no OOD detector; report `false` explicitly rather than erroring.
is_ood(::AdviColocResult) = false

# The ADVI Bayes factor is intrinsically a posterior-vs-prior comparison, so the single-arg
# interface accessor is not meaningful; delegate to `compute_BayesFactor` via a 2-arg method.
bayes_factor(posterior::AdviColocResult, prior::AdviColocResult; kwargs...) =
    first(compute_BayesFactor(posterior, prior; kwargs...))
bayes_factor(::AdviColocResult) = error(
    "The ADVI/Turing Bayes factor requires BOTH a posterior and a prior AdviColocResult; " *
    "call `bayes_factor(posterior, prior)` or `compute_BayesFactor(posterior, prior)`.")

######################################################################
# function to convert the posterior samples
######################################################################
function convert_posterior_samples(samples::Array{Float64, 2}, m::T) where {T <: Turing.DynamicPPL.Model}
    # get parameter_names
    parameter_names = Turing.DynamicPPL.syms(Turing.DynamicPPL.VarInfo(m))
    # select only the necessary parameters
    parameter_names = collect(parameter_names[1:10])

    # permute the samples
    samples = samples[1:10,:]
    samples = DataFrame(permutedims(samples, [2, 1]), parameter_names)

    # undo Fisher z transformation
    #samples = tanh.(samples)
end

######################################################################
# compute_BayesFactor (KDE + numerical integration baseline)
######################################################################
"""
    compute_BayesFactor(posterior::AdviColocResult, prior::AdviColocResult; ρ_threshold = 0.0)

KDE + numerical-integration Bayes factor for the ADVI reference path. Internal reference only.
"""
function compute_BayesFactor(posterior::AdviColocResult, prior::AdviColocResult; ρ_threshold::Float64 = 0.0)
    Δρ_post = posterior.posterior.μ_sample .- posterior.posterior.μ_control
    Δρ_prior = prior.posterior.μ_sample .- prior.posterior.μ_control
    # computing the probability of Δρ <=0
    posterior_dist = kde(Δρ_post)
    p_post, ϵ_post = quadgk(x -> pdf(posterior_dist, x), -Inf, ρ_threshold)
    p_post = 1 - p_post

    ϵ_post > 1e-5 &&
        @warn "CDF of the posterior distribution is approximated by numerical integration
        with an error of $ϵ_post that is unusually large. "

    # prior
    prior_dist = kde(Δρ_prior)
    p_prior, ϵ_prior = quadgk(x -> pdf(prior_dist, x), -Inf, ρ_threshold)
    p_prior = 1 - p_prior

    ϵ_prior > 1e-5 &&
        @warn "CDF of the prior distribution is approximated by numerical integration
        with an error of $ϵ_prior that is unusually large. "
    ##################### compute the Bayes factor #####################
    # compute prior odds
    prior_odds = p_prior / (1 - p_prior)
    # compute posterior odds
    posterior_odds = p_post / (1 - p_post)
    # compute Bayes factor
    bayes_factor = posterior_odds / prior_odds
    return(bayes_factor, p_post, p_prior)
end

######################################################################
# _prepare_data (patch/correlation summary input for the ADVI model)
######################################################################
function _prepare_data(img::MultiChannelImageStack, channels::Vector{T}, num_patches::T = 1; cor_method::Symbol = :pearson) where T <: Int
    # extract channels and patch
    n_images = length(img)
    sample_image::Array{Union{Float64, Missing}, 3} = fill(0.0, n_images, num_patches, num_patches)
    for (image,idx) ∈ zip(img, 1:n_images)
        x = image.data[channels[1]]
        y = image.data[channels[2]]
        x,y = patch.([x, y], num_patches)
        sample_image[idx,:,:] = correlation(x, y, method = cor_method)
    end

    # reshape the data
    sample_img = reshape(sample_image, n_images, num_patches^2)

    sample_data = fill(Vector{Float64}(), n_images)
    for (row, idx) in zip(eachrow(sample_img), 1:n_images)
        # Single-pass filter: collect non-missing, non-NaN values
        sample_data[idx] = filter(x -> !isnan(x), collect(skipmissing(row)))
    end

    # remove images with no signal above background (empty after filtering)
    delete_idx = Int[]
    for idx ∈ 1:n_images
        if isempty(sample_data[idx])
            @warn "Image $idx has no signal above background for the selected channels.
            The image is removed from the analysis."
            push!(delete_idx, idx)
        end
    end
    length(delete_idx) >= 1 && deleteat!(sample_data, delete_idx)


    if length(sample_data) == 0
        return nothing
    end

    return(sample_data)
end

######################################################################
# colocalization (Turing @model + ADVI variational inference)
######################################################################
"""
    colocalization(img, control, channels, num_patches = 1; iter, posterior_samples, cor_method)

Turing/ADVI reference colocalization. Returns `(prior, posterior)` as `AdviColocResult`.
Internal reference only (D-01); not exported.
"""
function colocalization(
    img::MultiChannelImageStack,
    control::MultiChannelImageStack,
    channels::Vector{T},
    num_patches::T = 1;
    iter::T = 1000,
    posterior_samples::T = 100_000,
    cor_method::Symbol = :pearson
    ) where T <: Int

    ctrl_data = _prepare_data(control, channels, num_patches, cor_method = cor_method)
    sample_data = _prepare_data(img, channels, num_patches, cor_method = cor_method)

    isnothing(ctrl_data) && return missing, missing
    isnothing(sample_data) && return missing, missing

    @model function model(control, sample)
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

    # define model
    m = model(ctrl_data, sample_data)
    # get prior
    prior_chain = sample(m, Prior(), posterior_samples)

    prior = AdviColocResult(
        img, control,
        channels, num_patches,
        DataFrames.DataFrame(prior_chain[[
            :μ_control,:ν_control,:σ_control,:τ_control,
            :μ_sample,:ν_sample,:σ_sample,:τ_sample
            ]]),
        prior_chain
    )

    ######################################################
    # get the posterior samples
    # calculate number of latent variables
    num_latent = size(DataFrames.DataFrame(prior_chain))[2] - 3
    # sample
    q = vi(m, ADVI(num_latent, iter))

    # get the posterior samples
    q_samples = rand(q, posterior_samples)
    q_samples = convert_posterior_samples(q_samples, m)

    posterior = AdviColocResult(
        img, control,
        channels,num_patches,
        q_samples, q
        )

    ######################################################
    return prior, posterior
end

######################################################################
# ADVI-result plotting (moved from src/plot.jl; dispatch on AdviColocResult)
######################################################################
"""
    plot_posterior(posterior::AdviColocResult; file, save, fig_size, dpi)

Diagnostic posterior plots for the ADVI reference result. Internal reference only.
"""
function plot_posterior(
    posterior::AdviColocResult; file::String = "posterior.png",
    save::Bool = true, fig_size::Vector{T} = [16.0,12.0], dpi::I= 300
    ) where {T <: AbstractFloat, I <: Integer}

    fig_size_1 = 72*fig_size[1]/2.54
    fig_size_2 = 72*fig_size[2]/2.54
    fig = GLMakie.Figure(size = (fig_size_1, fig_size_2), backgroundcolor = :white, fontsize = 8, figure_padding = (1,8,1,1))

    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "ρ", ylabel = "P(ρ|data)", title = "P(ρ|data)",
        xticklabelrotation = deg2rad(60)
        )
    #
    # alpha and normalise no longer supported, instead colour can be used
    GLMakie.density!(
        ax1, posterior.posterior.μ_control,
        colormap = (:viridis, 0.3), label = "control",
        )

    GLMakie.density!(
        ax1, posterior.posterior.μ_sample,
        colormap = (:viridis, 0.3), label = "sample"
        )

    ax2 = GLMakie.Axis(fig[1, 2], xlabel = "ν", ylabel = "P(ν|data)", title = "P(ν|data)", xticklabelrotation = deg2rad(60))

    GLMakie.density!(
        ax2, posterior.posterior.ν_control,
        colormap = (:viridis, 0.3), label = "ν_control"
        )

    GLMakie.density!(
        ax2, posterior.posterior.ν_sample,
        colormap = (:viridis, 0.3), label = "ν_sample"
    )

    ax3 = GLMakie.Axis(fig[1, 3], xlabel = "σ", ylabel = "P(σ|data)", title = "P(σ|data)", xticklabelrotation = deg2rad(60))

    GLMakie.density!(
        ax3, posterior.posterior.σ_control,
        colormap = (:viridis, 0.3), label = "σ_control"
        )

    GLMakie.density!(
        ax3, posterior.posterior.σ_sample,
        colormap = (:viridis, 0.3), label = "σ_sample"
    )

    ax4 = GLMakie.Axis(fig[1, 4], xlabel = "τ", ylabel = "P(τ|data)", title = "P(τ|data)", xticklabelrotation = deg2rad(60))

    GLMakie.density!(
        ax4, posterior.posterior.τ_sample,
        colormap = (:viridis, 0.3), label = "τ_sample"
        )

    GLMakie.density!(
        ax4, posterior.posterior.τ_control,
        colormap = (:viridis, 0.3), label = "τ_control"
    )

    GLMakie.Legend(
        fig[1,5], ax1, orientation = :vertical, framevisible = false,
        labelsize = 8, rowgap = 0, titlegap = 0,
        tellheight = false, valign = :top,
        patchsize = (5,5)
        )
    #########
    # Δρ
    Δρ = posterior.posterior.μ_sample .- posterior.posterior.μ_control

    ax5 = GLMakie.Axis(
        fig[2, 1:5], xlabel = "Δρ", ylabel = "P(Δρ|data)", title = "P(Δρ|data)",
        limits = (-2, 2, nothing, nothing),
        xticks = (collect(-2.0:0.1:2.0), string.(collect(-2.0:0.1:2.0))),
        xticklabelrotation = deg2rad(60)
        )

    GLMakie.density!(
        ax5, Δρ, label = "Δρ",
        color = (:darkgrey, 0.3)
    )

    GLMakie.vlines!(ax5, 0, color = :black, linestyle = :dash)

    # colgaps
    GLMakie.colgap!(fig.layout,1,GLMakie.Relative(0.02))
    GLMakie.colgap!(fig.layout,2,GLMakie.Relative(0.02))
    GLMakie.colgap!(fig.layout,3,GLMakie.Relative(0.02))
    GLMakie.colgap!(fig.layout,4,GLMakie.Relative(0.02))

    # colsize of legend
    GLMakie.colsize!(fig.layout,5,GLMakie.Relative(0.11))
    # rowgap between upper and lower plot
    GLMakie.rowgap!(fig.layout,1,GLMakie.Relative(0.02))

    save && GLMakie.save(file, fig , px_per_unit = dpi/72)
    return(fig)
end

"""
    bayesplot(prior::AdviColocResult, posterior::AdviColocResult, bf::Float64; kwargs...)

Plot the prior/posterior Δρ distributions annotated with the Bayes factor. Internal reference.
"""
function bayesplot(
    prior::AdviColocResult,
    posterior::AdviColocResult,
    bf::Float64;
    file::String = "bayesplot.png",
    save::Bool = true,
    ρ_threshold::Float64 = 0.0,
    fig_size::Vector{T} = [8.0,10.0],
    dpi::I= 300
    ) where {T <: AbstractFloat, I <: Integer}

    Δρ_post = posterior.posterior.μ_sample .- posterior.posterior.μ_control
    Δρ_prior = prior.posterior.μ_sample .- prior.posterior.μ_control

    fig_size_1 = 72*fig_size[1]/2.54
    fig_size_2 = 72*fig_size[2]/2.54
    fig = GLMakie.Figure(
        size = (fig_size_1, fig_size_2), backgroundcolor = :white,
        fontsize = 8, figure_padding = (1,8,1,1)
        )

    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "Δρ", ylabel = "PDF", title = "P(Δρ|data)",
        xticks = (collect(-2.0:0.4:2.0), string.(collect(-2.0:0.4:2.0))),
        xticklabelrotation = deg2rad(60),
        limits = (-2, 2, nothing, nothing)
    )

    hist1a = GLMakie.density!(
        ax1, Δρ_prior, colormap = (:viridis, 0.3),
        label = "prior",
    )

    hist1b = GLMakie.density!(
        ax1, Δρ_post,
        label = "posterior", colormap = (:viridis, 0.3)
    )

    GLMakie.vlines!(
        ax1, ρ_threshold, color = :black,
        linestyle = :dash, label = "Δρ = $ρ_threshold"
        )

    t = "BF = $(round(bf; digits = 4))"
    # Add BF to the plot
    GLMakie.text!(-1.9,0.05, text = t)
    GLMakie.Legend(
        fig[2,1], ax1, orientation = :horizontal, framevisible = false,
        labelsize = 8, rowgap = 0, titlegap = 0,
        tellheight = false, valign = :top,
        patchsize = (5,5)
        )

    GLMakie.rowsize!(fig.layout,1,GLMakie.Relative(0.90))
    GLMakie.rowgap!(fig.layout,1,GLMakie.Relative(0.01))
    save && GLMakie.save(file, fig, px_per_unit = dpi/72)
    return fig
end

"""
    bayes_rangeplot(prior::AdviColocResult, posterior::AdviColocResult; Δρ, save, file, fig_size, dpi)

Plot the Bayes factor across a range of Δρ thresholds. Internal reference.
"""
function bayes_rangeplot(
    prior::AdviColocResult,
    posterior::AdviColocResult;
    Δρ::Vector{T} = collect(range(-0.5,0.5;step =0.05)),
    save::Bool = true,
    file::String = "bayes_rangeplot.png",
    fig_size::Vector{T} = [8.0,10.0],
    dpi::I= 300
    ) where {T <: AbstractFloat, I <: Integer}

    # Pre-compute KDE distributions once (avoiding repeated computation in loop)
    Δρ_post = posterior.posterior.μ_sample .- posterior.posterior.μ_control
    Δρ_prior = prior.posterior.μ_sample .- prior.posterior.μ_control
    posterior_dist = kde(Δρ_post)
    prior_dist = kde(Δρ_prior)

    # calculate the bayes factor for each Δρ threshold using cached KDEs
    bf = fill(0.0, length(Δρ))
    for idx in eachindex(Δρ)
        ρ_threshold = Δρ[idx]
        # Compute probability using cached KDE
        p_post, _ = quadgk(x -> pdf(posterior_dist, x), -Inf, ρ_threshold)
        p_post = 1 - p_post
        p_prior, _ = quadgk(x -> pdf(prior_dist, x), -Inf, ρ_threshold)
        p_prior = 1 - p_prior
        # Compute Bayes factor
        prior_odds = p_prior / (1 - p_prior)
        posterior_odds = p_post / (1 - p_post)
        a = posterior_odds / prior_odds
        a > 0 ? bf[idx] = a : bf[idx] = NaN
    end

    # return if all values are NaN or Inf
    if all(isnan.(bf))
        @warn "All values are bayes factors are zero. No plot is generated."
        return nothing
    end

    if all(isinf.(bf))
        @warn "All values are bayes factors are infinite. No plot is generated."
        return nothing
    end

    # plot the results
    ticks = collect(range(Δρ[1], Δρ[end], length = 11))
    fig_size_1 = 72*fig_size[1]/2.54
    fig_size_2 = 72*fig_size[2]/2.54

    fig = GLMakie.Figure(
        size = (fig_size_1, fig_size_2), backgroundcolor = :white,
        fontsize = 8, figure_padding = (1,8,1,1)
        )

    ax1 = GLMakie.Axis(
        fig[1, 1], xlabel = "Δρ0",
        ylabel = "log10(BF[Δρ > Δρ0 : Δρ ≤ Δρ0])",
        title = "Bayes factor vs Δρ",
        limits = (Δρ[1], Δρ[end], nothing, nothing),
        xticks = (ticks, string.(ticks)),
        xticklabelrotation = deg2rad(60)
    )

    GLMakie.lines!(ax1, Δρ, log10.(bf))

    GLMakie.hlines!(ax1, [0], color = :black, linestyle = :dash, label = "BF = 1")
    save && GLMakie.save(file, fig, px_per_unit = dpi/72)
    return fig
end

end # module ProteinCoLocTuringExt
