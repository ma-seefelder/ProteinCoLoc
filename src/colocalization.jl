## colocalization.jl
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Manuel Seefelder
######################################################################
# function to patch the image
"""
    patch(img::Array{Float64, 2}, num_patches::Int64)
    Patches the image into num_patches x num_patches patches.
    Return a 4D array with the patches.
"""
function patch(img::Array{Float64, 2}, num_patches::Int64)
    # calculate the size of the patches
    patch_size_x = Int64(floor(size(img, 1) / num_patches))
    patch_size_y = Int64(floor(size(img, 2) / num_patches))
    # calculate the number of patches
    num_patches_x = Int64(floor(size(img, 1) / patch_size_x))
    num_patches_y = Int64(floor(size(img, 2) / patch_size_y))
    # initialize the patches
    patches = zeros(Union{Float64, Missing}, num_patches_x, num_patches_y, patch_size_x, patch_size_y)
    # loop over the patches
    for i in 1:num_patches_x
        for j in 1:num_patches_y
            patches[i, j, :, :] = img[(i-1)*patch_size_x+1:i*patch_size_x, (j-1)*patch_size_y+1:j*patch_size_y]
        end
    end
    return patches
end

######################################################################
# function to calculate the correlation

"""
    _exclude_zero!(a::Vector{T},b::Vector{T}) where T <: Number
    Exclude all values that are zero in at least one of the vectors a and b.
    Return the vectors a and b without the zero and missing values as a Vector{Float64}.
"""
function _exclude_zero!(a::Vector{T},b::Vector{T}) where T <: Number  
    # get the indices of the zero values
    a_zero = append!(findall(a .== 0), findall(isnan.(a)))
    b_zero = append!(findall(b .== 0), findall(isnan.(b)))
    # make union of the indices
    zero_indices = sort(union(a_zero, b_zero))
    # exclude all values at the indices in zero_indices
    deleteat!(a, zero_indices)
    deleteat!(b, zero_indices)
end

function _exclude_zero!(a::Vector{Union{T,Missing}}, b::Vector{Union{T,Missing}}) where T <: Number
    # change missing values to zero
    a[ismissing.(a)] .= 0
    b[ismissing.(b)] .= 0

    # get the indices of the zero values
    a_zero = append!(findall(a .== 0), findall(isnan.(a)))
    b_zero = append!(findall(b .== 0), findall(isnan.(b)))
    # make union of the indices
    zero_indices = sort(union(a_zero, b_zero))
    # exclude all values at the indices in zero_indices
    deleteat!(a, zero_indices)
    deleteat!(b, zero_indices)
end

"""
    correlation(x::Array{Float64, 4}, y::Array{Float64, 4})
    Calculate the correlation between two 4D arrays.
    Return a 2D array with the correlation for each patch.
"""
function correlation(x::Array{T, 4}, y::Array{T, 4}) where T <: Union{Float64, Missing}
    # get the number of patches
    num_patches = size(x, 1)
    # initialize the correlation
    ρ = zeros(Union{Float64, Missing},num_patches, num_patches)
    # loop over the patches
    for i in 1:num_patches
        for j in 1:num_patches
            a = x[i, j, :, :][:] 
            b = y[i, j, :, :][:]
            _exclude_zero!(a,b)
            length(a) <= 3 ? ρ[i, j] = missing : ρ[i, j] = cor(a, b)
        end
    end
    return ρ
end
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

function plot_posterior(posterior::CoLocResult, prior::CoLocResult)
    hist1 = Plots.histogram(
        posterior.posterior.μ_control, legend = true, label = "μ_control", 
        title = "P(ρ|data)", alpha = 0.5
        )
    Plots.histogram!(hist1, posterior.posterior.μ_sample, label = "μ_sample", alpha = 0.5)

    hist2 = Plots.histogram(
        posterior.posterior.ν_control, label = "ν_control", 
        legend = true, title = "P(ν|data)", alpha = 0.5
        )

    Plots.histogram!(hist2, posterior.posterior.ν_sample, label = "ν_sample", alpha = 0.5)

    hist3 = Plots.histogram(
        posterior.posterior.σ_control, label = "σ_control", 
        legend = true, title = "P(σ|data)", alpha = 0.5
        )

    Plots.histogram!(hist3, posterior.posterior.σ_sample, label = "σ_sample", alpha = 0.5)

    Δρ = posterior.posterior.μ_sample .- posterior.posterior.μ_control

    hist4 = Plots.histogram(
        Δρ,legend = true, label = "Δρ", 
        title = "P(Δρ|data)", alpha = 0.5)

    hist5 = Plots.histogram(
        posterior.posterior.τ_sample,legend = true,
        label = "τ_sample", title = "P(τ|data)", alpha = 0.5
        )

    Plots.histogram!(hist5, posterior.posterior.τ_control, label = "τ_control", alpha = 0.5)

    p = Plots.plot(
        hist1, hist2, hist3, hist5, hist4, 
        layout = (3, 2), size = (800, 800)
        )

    Plots.display(p)
    return(p)
end


"""
    compute_BayesFactor(posterior::CoLocResult, prior::CoLocResult)
"""
function compute_BayesFactor(posterior::CoLocResult, prior::CoLocResult)
    Δρ_post = posterior.posterior.μ_sample .- posterior.posterior.μ_control
    Δρ_prior = prior.posterior.μ_sample .- prior.posterior.μ_control
    # computing the probability of Δρ <=0
    posterior_dist = kde(Δρ_post)
    p_post, ϵ_post = quadgk(x -> pdf(posterior_dist, x), -Inf, 0.0)
    p_post = 1 - p_post

    ϵ_post > 1e-5 && 
        @warn "CDF of the posterior distribution is approximated by numerical integration 
        with an error of $ϵ_post that is unusually large. " 
    
    # prior
    prior_dist = kde(Δρ_prior)
    p_prior, ϵ_prior = quadgk(x -> pdf(prior_dist, x), -Inf, 0.0)
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

    sample_img = fill(0.0, img.num_images, num_patches, num_patches)
    for (image,idx) ∈ zip(img, 1:img.num_images)
        x = image.data[channels[1]]
        y = image.data[channels[2]]
        x,y = patch.([x, y], num_patches)
        sample_img[idx,:,:] = correlation(x, y)
    end

    ctrl_img = fill(0.0, control.num_images, num_patches, num_patches)
    for (image,idx) ∈ zip(control, 1:control.num_images)
        x = image.data[channels[1]]
        y = image.data[channels[2]]
        x,y = patch.([x, y], num_patches)
        ctrl_img[idx,:, :] = correlation(x, y)
    end

    ################################
    # Turing model
    ################################
    # input data for the model are the correlation coefficients of the sample and control for each image 
    # type: Array{Float64, 2} with size (num_images, num_patches^2)
    sample_img = reshape(sample_img, img.num_images, num_patches^2)
    ctrl_img = reshape(ctrl_img, control.num_images, num_patches^2)

    @model function model(control::Matrix{Float64}, sample::Matrix{Float64})
        
        # get the number of patches
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
            control[idx] ~ TDist(ν_control_image[idx]) * σ_control_image[idx] + μ_control_image[idx]
        end

        for idx ∈ 1:num_sample
            sample[idx] ~ TDist(ν_sample_image[idx]) * σ_sample_image[idx] + μ_sample_image[idx]
        end       
    end

    # define model
    m = model(ctrl_img, sample_img)
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
    posterior_samples::T = 100_000
    ) where T <: Int

    bf = fill(0.0, length(num_patches))
    p_post = fill(0.0, length(num_patches))
    p_prior = fill(0.0, length(num_patches))

    for (n_patch, idx) ∈ zip(num_patches, 1:length(num_patches))
        prior, posterior = colocalization(
            img, control, channels, n_patch; 
            iter = iter, posterior_samples = posterior_samples)
        bf[idx], p_post[idx], p_prior[idx] = compute_BayesFactor(posterior, prior)
    end

    # plot the results
    p1 = Plots.plot(
        num_patches.^2, log10.(bf), 
        xlabel = "number of patches", 
        ylabel = "log10(Bayes factor)", 
        title = "Bayes factor vs number of patches",
        legend = false
        )

    p2 = Plots.plot(
        num_patches.^2, p_post, 
        xlabel = "number of patches", 
        ylabel = "P(Δρ > 0| data)", 
        title = "P(Δρ > 0| data) vs number of patches",
        legend = false
        )

    p3 = Plots.plot(
        num_patches.^2, p_prior, 
        xlabel = "number of patches", 
        ylabel = "P(Δρ > 0)", 
        title = "P(Δρ > 0) vs number of patches",
        legend = false
        )

    p = Plots.plot(
        p1, p2, p3,
        layout = (3, 1), size = (800, 800)
        )
    
    Plots.display(p)
    return p, bf, p_post, p_prior 
end