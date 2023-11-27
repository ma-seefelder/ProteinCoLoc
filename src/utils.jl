#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel
E-Mail: proteincoloc@protonmail.com
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
=#

# get_files
function get_images(path::S, nchannels::I, stack_name::S; mask::Bool = true) where {I<:Integer, S<:AbstractString}
    # get files in image path
    files = readdir(path, join=false)
    # remove subdirectories
    files = filter(x -> contains(x,"."), files)
    # determine number of images in folder
    nimages = length(files) / nchannels
    # check if the number of images is an integer 
    nimages != floor(nimages) && error("The number of images in the image folder is not a multiple of the number of channels.")
    # check if the number of images is at least 1
    nimages < 1 && error("The number of images in the image folder must be at least 1.")

    # generate vector wit channel names
    channel_names = string.(collect(1:nchannels))
    # group images
    images = Vector{MultiChannelImage{Float64,String, Float64}}(undef, Int(nimages))

    # retrieve image names
    image_names = Vector{String}()
    for i in split.(files, "_")
        length(i) > 2 ? img_name = join(i[1:end-1],"_") : img_name = i[1]
        !in(img_name,image_names) ? push!(image_names,img_name) : nothing
    end

    # load images
    for i in eachindex(images) 
        image_name = image_names[i]
        channel_files = fill("", nchannels)
        for j in 1:nchannels
            file_index = findfirst(file -> occursin(image_name * "_c" * string(j), file), files)
            !isnothing(file_index) || error("The image path is not valid.")
            channel_files[j] = path * "/" * files[file_index]
        end
        images[i] = MultiChannelImage(image_name, channel_files, channel_names)
        # mask image
        if mask
            images[i] = ProteinCoLoc._apply_mask!(images[i], ProteinCoLoc._calculate_mask(images[i]))
        end
    end
    # return files
    return MultiChannelImageStack(images, stack_name)
end 

# iterate over images in image stack and plot local correlation for each image
function plot_images(
    plot_type::Symbol,
    image_stack::MultiChannelImageStack, 
    num_patches::Int, 
    channel_indices::Vector{Int}, 
    base_filename::AbstractString, 
    suffix::AbstractString=""
    )

    for image in image_stack
        # plot local correlation
        if plot_type == :local_correlation
            try 
                local_correlation_plot(image, num_patches, channel_indices, file = "$base_filename$suffix"*"_$(image.:name).png")
            catch
                @warn "The local correlation plot for image $(image.:name) could not be generated."
            end
        elseif plot_type == :patched_correlation
            try 
                plot(image, num_patches, channel_indices, file = "$base_filename$suffix"*"_$(image.:name).png")
            catch
                @warn "The patched correlation plot for image $(image.:name) could not be generated."
            end
        else
            error("The plot type is not valid.")
        end
    end
end

# helper function to plot all the plots for a given image and control image stack
"""
    generate_plots(
    images, control_images, channel_selection_two, number_patches, number_patches_loc, 
    number_iterations, number_posterior_samples, ρ_threshold, 
    ρ_range, ρ_range_step, output_folder_path, patched_correlation_plt, local_correlation_plt, 
    bayes_factor_plt, bayes_range_plt, posterior_plt
    )

    Generate all plots for a given image and control image stack.
"""
function generate_plots(
    images, control_images, channel_selection_two, number_patches, number_patches_loc, 
    number_iterations, number_posterior_samples, ρ_threshold, 
    ρ_range, ρ_range_step, output_folder_path, patched_correlation_plt, local_correlation_plt, 
    bayes_factor_plt, bayes_range_plt, posterior_plt
    )

    prior, posterior = colocalization(
        images, control_images, channel_selection_two, 
        number_patches; iter = number_iterations, posterior_samples = number_posterior_samples
    )

    if patched_correlation_plt
        base_file = "$output_folder_path/patched_correlation_c$(channel_selection_two[1])_c$(channel_selection_two[2])"
        plot_images(:patched_correlation, images, number_patches, channel_selection_two, base_file)
        plot_images(:patched_correlation, control_images, number_patches, channel_selection_two, base_file, "_control")
    end

    if local_correlation_plt
        base_file = "$output_folder_path/local_correlation_c$(channel_selection_two[1])_c$(channel_selection_two[2])"
        plot_images(:local_correlation, images, number_patches_loc, channel_selection_two, base_file)
        plot_images(:local_correlation, control_images, number_patches_loc, channel_selection_two, base_file, "_control")
    end

    bf, _, _ = compute_BayesFactor(posterior, prior; ρ_threshold = ρ_threshold)
    if bayes_factor_plt
        base_file = "$output_folder_path/bayes_factor_c$(channel_selection_two[1])_c$(channel_selection_two[2])"
        try 
            bayesplot(prior, posterior, bf; file = base_file * ".png", ρ_threshold = ρ_threshold)
        catch
            @warn "The bayes factor plot could not be generated."
        end
    end

    if bayes_range_plt
        base_file = "$output_folder_path/bayes_factor_range_c$(channel_selection_two[1])_c$(channel_selection_two[2])"
        try 
            bayes_rangeplot(prior, posterior; file = base_file * ".png", Δ̢ρ = collect(range(ρ_range[1],ρ_range[2];step =ρ_range_step)))
        catch
            @warn "The bayes factor range plot could not be generated."
        end
    end

    if posterior_plt
        base_file = output_folder_path * "/posterior_c" * string(channel_selection_two[1]) * "_c" * string(channel_selection_two[2])
        try 
            plot_posterior(posterior; file = base_file * ".png")
        catch
            @warn "The posterior plot could not be generated."
        end
    end
end

function combinations2(n)
    # generate all possible combinations of two elements from 1:n without repetition 
    combination = Vector{Vector{Int}}()
    for i in 1:n
        for j in 1:n
            if [j,i] ∉ combination && i != j
                push!(combination, [i,j])
            end
        end
    end
    return combination
end