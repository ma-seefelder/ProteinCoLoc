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

######################################################################
# function to convert the posterior samples
######################################################################
"""
    convert_posterior_samples(
    samples::Array{Float64, 2}, 
    m::T
    ) where {T <: DynamicPPL.Model}

This function converts the posterior samples from the ADVI algorithm into a DataFrame.

# Arguments
- `samples`: A 2D array of Float64 representing the posterior samples from the ADVI algorithm.
- `m`: A DynamicPPL.Model representing the Bayesian model.

# Returns
- `samples`: A DataFrame where each column represents a parameter of the Bayesian model and each row represents a sample from the posterior distribution.

# Notes
This function gets the parameter names from the Bayesian model, selects the necessary parameters, permutes the samples, and converts them into a DataFrame. The Fisher z transformation is undone by applying the hyperbolic tangent function to the samples.
"""
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
"""
    struct CoLocResult

This struct represents the result of a Bayesian colocalization analysis.

# Fields
- `img`: A MultiChannelImageStack representing the sample images.
- `control`: A MultiChannelImageStack representing the control images.
- `channels`: A Vector of integers representing the channels that were analyzed.
- `num_patches`: An integer representing the number of patches that were analyzed.
- `posterior`: A DataFrame representing the posterior distribution.
- `advi_result`: The result of the ADVI algorithm.

# Notes
This struct is used to store the result of a Bayesian colocalization analysis. It includes the sample and control images, the channels and number of patches that were analyzed, the posterior distribution, and the result of the ADVI algorithm.
"""
struct CoLocResult
    img::MultiChannelImageStack
    control::MultiChannelImageStack
    channels::Vector{Int64}
    num_patches::Int64
    posterior::DataFrame
    advi_result
end

"""
    compute_BayesFactor(
    posterior::CoLocResult, 
    prior::CoLocResult; 
    ρ_threshold::Float64 = 0.0
    )

This function computes the Bayes factor for the colocalization of two proteins.

# Arguments
- `posterior`: A CoLocResult object representing the posterior distribution.
- `prior`: A CoLocResult object representing the prior distribution.
- `ρ_threshold`: A Float64 representing the threshold for the difference in correlation. Default is 0.0.

# Returns
- `bayes_factor`: A Float64 representing the Bayes factor.
- `p_post`: A Float64 representing the probability of Δρ <= ρ_threshold under the alternative hypothesis.
- `p_prior`: A Float64 representing the probability of Δρ <= ρ_threshold under the null hypothesis.

# Errors
- Throws a warning if the error of the numerical integration for the CDF of the prior or posterior distribution is unusually large.

# Notes
The Bayes factor is computed as the ratio of the posterior odds and the prior odds. The prior and posterior odds are computed as the ratio of the probability of Δρ <= ρ_threshold under the alternative hypothesis and the null hypothesis. The prior and posterior distributions are approximated by a kernel density estimation (KDE), and the probability of Δρ <= ρ_threshold is approximated by numerical integration.
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


"""
    _prepare_data(
    img::MultiChannelImageStack, 
    channels::Vector{T}, 
    num_patches::T = 1
    ) where T <: Int

This is a low-level function that prepares the data for Bayesian analysis.

# Arguments
- `img`: A MultiChannelImageStack representing the sample images.
- `channels`: A Vector of integers representing the channels to be analyzed.
- `num_patches`: An integer representing the number of patches to be analyzed. Default is 1.
- `cor_method`: A Symbol representing the correlation method to be used. Default is :pearson. The other options are :spearman and :kendall.

# Returns
- `sample_data`: An array of vectors, where each vector contains the correlation values of a specific image. Images with no signal above the background for the selected channels are removed from the analysis.

# Errors
- Returns `nothing` if all images have no signal above the background for the selected channels.

# Notes
This function extracts the specified channels from the images, applies patching, and calculates the correlation between the channels. The resulting data is reshaped and missing values are skipped. Images with no signal above the background for the selected channels are removed from the analysis.
"""
function _prepare_data(img::MultiChannelImageStack, channels::Vector{T}, num_patches::T = 1; cor_method::Symbol = :pearson) where T <: Int
    # extract channels and patch 
    sample_image::Array{Union{Float64, Missing}, 3} = fill(0.0, img.num_images, num_patches, num_patches)
    for (image,idx) ∈ zip(img, 1:img.num_images)
        x = image.data[channels[1]]
        y = image.data[channels[2]]
        x,y = patch.([x, y], num_patches)
        sample_image[idx,:,:] = correlation(x, y, method = cor_method)
    end

    # reshape the data
    sample_img = reshape(sample_image, img.num_images, num_patches^2)

    sample_data = fill(Vector{Float64}(), img.num_images)
    for (row, idx) in zip(eachrow(sample_img), 1:img.num_images)
        sample_data[idx] = collect(skipmissing(row))
        # remove NaN values
        sample_data[idx] = sample_data[idx][.!isnan.(sample_data[idx])]
    end

    # remove images with no signal above background
    for idx ∈ 1:img.num_images
        size_array = size(sample_data[idx])
        if ismissing.(sample_data[idx]) == fill(true, size_array)
            @warn "Image $idx has no signal above background for the selected channels.
            The image is removed from the analysis."
            deleteat!(sample_data, idx)
        end
    end

    if length(sample_data) == 0
        return nothing
    end

    return(sample_data)
end


"""
    colocalization(
    img::MultiChannelImageStack, 
    control::MultiChannelImageStack, 
    channels::Vector{T},
    num_patches::T = 1;
    iter::T = 1000, 
    posterior_samples::T = 100_000
    ) where T <: Int

This function performs a Bayesian colocalization analysis on multi-channel image stacks. 

# Arguments
- `img`: A MultiChannelImageStack representing the sample images.
- `control`: A MultiChannelImageStack representing the control images.
- `channels`: A Vector of integers representing the channels to be analyzed.
- `num_patches`: An integer representing the number of patches to be analyzed. Default is 1.
- `iter`: An integer representing the number of iterations for the ADVI algorithm. Default is 1000.
- `posterior_samples`: An integer representing the number of posterior samples to be generated. Default is 100,000.

# Returns
- `prior`: A CoLocResult object representing the prior distribution.
- `posterior`: A CoLocResult object representing the posterior distribution.

# Errors
- Throws an error if no control or sample images with signal above background for the selected channels are found.

# Notes
This function uses a Bayesian model with global and local priors for the control and sample images. The model is defined using the Turing.jl package. The ADVI algorithm is used for variational inference, and the prior and posterior distributions are returned as CoLocResult objects.
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

    """
    @model function model(control, sample)

    This function defines a Bayesian model for analyzing multi-channel image stacks.

    # Arguments
    - `control`: An array of vectors, where each vector contains the correlation values of a specific control image.
    - `sample`: An array of vectors, where each vector contains the correlation values of a specific sample image.

    # Model
    The model includes global priors for the control and sample images, and local priors for each image. The global priors include the mean, degrees of freedom, and standard deviation of the control and sample images, and the patch heterogeneity. The local priors include the mean, degrees of freedom, and standard deviation of each image.

    The likelihood of the model is defined by the t-distribution, with the degrees of freedom, standard deviation, and mean of each image.

    # Notes
    This function uses the Turing.jl package to define the model. The model is defined using the `@model` macro, which allows for a flexible specification of probabilistic models in Julia.
    """
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
    num_latent = size(DataFrames.DataFrame(prior_chain))[2] - 3
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
