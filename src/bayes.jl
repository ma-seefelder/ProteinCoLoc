#########################################################################################
## bayes.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder
##########################################################################################

######################################################################
# function to convert the posterior samples
######################################################################

function convert_posterior_samples(samples::Array{Float64, 2}, m::T) where {T <: DynamicPPL.Model}
    # get parameter_names
    parameter_names = DynamicPPL.syms(DynamicPPL.VarInfo(m))
    # select only the necessary parameters
    parameter_names = collect(parameter_names[1:10])

    # permute the samples
    samples = samples[1:10,:]
    samples = DataFrame(permutedims(samples, [2, 1]), parameter_names)

    # undo Fisher z transformation
    #samples = tanh.(samples)
end

######################################################################
# plot posterior
######################################################################
struct CoLocResult
    img::MultiChannelImageStack
    control::MultiChannelImageStack
    channels::Vector{Int64}
    num_patches::Int64
    posterior::DataFrame
    advi_result
end

"""
    compute_BayesFactor(posterior::CoLocResult, prior::CoLocResult; ρ_theshold::Float64 = 0.0)

    Compute the Bayes factor for the colocalization of two proteins.
    H1: Δ̢ > ρ_theshold; H0: Δ̢ <= ρ_threshold
    The Bayes factor is computed as the ratio of the posterior odds and the prior odds.
    The prior odds are computed as the ratio of the probability of Δρ <= 0 under the alternative hypothesis and the null hypothesis.
    The posterior odds are computed as the ratio of the probability of Δρ <= 0 under the alternative hypothesis and the null hypothesis.
    The prior and posterior distributions are approximated by a kernel density estimation (KDE) and the 
    probability of Δρ <= ρ_threshold is approximated by numerical integration.
"""
function compute_BayesFactor(posterior::CoLocResult, prior::CoLocResult; ρ_threshold::Float64 = 0.0)
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

function colocalization(
    img::MultiChannelImageStack, 
    control::MultiChannelImageStack, 
    channels::Vector{T},
    num_patches::T = 1;
    iter::T = 1000, 
    posterior_samples::T = 100_000
    ) where T <: Int

    sample_image::Array{Union{Float64, Missing}, 3} = fill(0.0, img.num_images, num_patches, num_patches)
    for (image,idx) ∈ zip(img, 1:img.num_images)
        x = image.data[channels[1]]
        y = image.data[channels[2]]
        x,y = patch.([x, y], num_patches)
        sample_image[idx,:,:] = correlation(x, y)
    end

    ctrl_image::Array{Union{Float64, Missing}, 3} = fill(0.0, control.num_images, num_patches, num_patches)
    for (image,idx) ∈ zip(control, 1:control.num_images)
        x = image.data[channels[1]]
        y = image.data[channels[2]]
        x,y = patch.([x, y], num_patches)
        ctrl_image[idx,:, :] = correlation(x, y)
    end

    ################################
    # Turing model
    ################################
    # input data for the model are the correlation coefficients of the sample and control for each image 
    # type: Array{Float64, 2} with size (num_images, num_patches^2)
    sample_img = reshape(sample_image, img.num_images, num_patches^2)
    ctrl_img = reshape(ctrl_image, control.num_images, num_patches^2)

    sample_data = fill(Vector{Float64}(), img.num_images)
    for (row, idx) in zip(eachrow(sample_img), 1:img.num_images)
        sample_data[idx] = collect(skipmissing(row))
    end

    ctrl_data = fill(Vector{Float64}(), control.num_images)
    for (row, idx) in zip(eachrow(ctrl_img), 1:control.num_images)
        ctrl_data[idx] = collect(skipmissing(row))
    end

    @model function model(control, sample)
        # get the number of images
        num_control = size(control, 1)
        num_sample = size(sample, 1)

        # ============= gloabal priors (per biological condition) ============= #
        # mean, degrees of freedom and standard deviation of the control 
        μ_control ~ Truncated(Cauchy(0, 0.3),-1,1)
        ν_control ~ Exponential()
        σ_control ~ Truncated(Cauchy(0, 0.3),0,1)
        τ_control ~ Truncated(Cauchy(0, 0.3),0,1)
        # mean, degrees of freedom and standard deviation of the sample, and the patch heterogeneity
        μ_sample ~ Truncated(Cauchy(0, 0.3),-1,1)
        ν_sample ~ Exponential()
        σ_sample ~ Truncated(Cauchy(0, 0.3),0,1)
        τ_sample ~ Truncated(Cauchy(0, 0.3),0,1)

        # ============= local priors (per image) ============= #
        μ_control_image ~ filldist(Truncated(Normal(μ_control, σ_control),-1,1), num_control) # mean of the control for each patch
        ν_control_image ~ filldist(Exponential(ν_control), num_control) # degress of freedom of the control for each patch
        σ_control_image ~ filldist(Truncated(Normal(σ_control, τ_control),0,1), num_control) # standard deviation of the control for each patch

        μ_sample_image ~ filldist(Truncated(Normal(μ_sample, σ_sample),-1,1), num_sample) # mean of the sample for each patch
        ν_sample_image ~ filldist(Exponential(ν_sample), num_sample) # degress of freedom of the sample for each patch
        σ_sample_image ~ filldist(Truncated(Normal(σ_sample, τ_sample),0,1), num_sample) # standard deviation of the control for each patch

        # likelihood
        for idx ∈ 1:num_control
            #c = collect(skipmissing(control[idx, :]))
            control[idx] ~ TDist(ν_control_image[idx]) * σ_control_image[idx] + μ_control_image[idx]
        end

        for idx ∈ 1:num_sample
            #s = collect(skipmissing(sample[idx, :]))
            sample[idx] ~ TDist(ν_sample_image[idx]) * σ_sample_image[idx] + μ_sample_image[idx]
        end      
    end

    # define model
    m = model(ctrl_data, sample_data)
    # get prior
    prior_chain = sample(m, Prior(), posterior_samples)

    prior = CoLocResult(
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
    num_latent = 8 + size(ctrl_img)[1] * 3 + size(sample_img)[1] * 3
    # sample
    q = vi(m, ADVI(num_latent, iter))

    # get the posterior samples
    q_samples = rand(q, posterior_samples)
    q_samples = convert_posterior_samples(q_samples, m)

    posterior = CoLocResult(
        img, control, 
        channels,num_patches, 
        q_samples, q
        )

    ######################################################
    return prior, posterior
end

function bayesfactor_robustness(
    img::MultiChannelImageStack, 
    control::MultiChannelImageStack, 
    channels::Vector{T},
    num_patches::Vector{T} = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    iter::T = 1000, 
    posterior_samples::T = 100_000,
    save::Bool = true,
    file::String = "bayesfactor_robustness.png"
    ) where T <: Int

    bf = fill(0.0, length(num_patches))
    p_post = fill(0.0, length(num_patches))
    p_prior = fill(0.0, length(num_patches))

    for (n_patch, idx) ∈ zip(num_patches, 1:length(num_patches))
        println("n_patch: $n_patch; idx: $idx")
        prior, posterior = colocalization(
            img, control, channels, n_patch; 
            iter = iter, posterior_samples = posterior_samples)
        bf[idx], p_post[idx], p_prior[idx] = compute_BayesFactor(posterior, prior)
    end

    # plot the results
    fig = Figure(resolution = (800, 800))
    ax1 = Axis(fig[1, 1], xlabel = "number of patches", ylabel = "log10(Bayes factor)")
    lines!(ax1, num_patches.^2, log10.(bf))

    ax2 = Axis(fig[2, 1], xlabel = "number of patches", ylabel = "P(Δρ > 0| data)")
    lines!(ax2, num_patches.^2, p_post)

    ax3 = Axis(fig[3, 1], xlabel = "number of patches", ylabel = "P(Δρ > 0)")
    lines!(ax3, num_patches.^2, p_prior)

    # Save the figure
    save && GLMakie.save(file, fig)
    return fig, bf, p_post, p_prior 
end