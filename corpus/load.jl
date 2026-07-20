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

# corpus/load.jl --- two-channel TIFF anchor -> MultiChannelImage via READ-ONLY src reuse (SC1, D-03).
#
# The conversion layer that makes a physical anchor directly usable by Phase 16: two single-channel
# TIFFs (the anchor layout) load into a `MultiChannelImage` and reduce through the UNCHANGED
# `src/` summary path (`patch` + `correlation`, fixed 8x8, unchanged since Phase 2) to a finite
# correlation grid. This is the D-03 preference — anchors that convert to a MultiChannelImage so
# `src/` summary functions apply verbatim.
#
# DECOUPLING (hard constraint, CLAUDE.md): this file reaches `src/` ONLY via a guarded, READ-ONLY
# `include()` of the frozen `src/LoadImages.jl` + `src/colocalization.jl`. It makes NO edits to
# `src/` and touches no root Project.toml — it merely re-uses the frozen constructors/summary
# functions. A single composite/RGB TIFF (the CBS layout) is handled SEPARATELY in 08-04
# (`corpus/cbs.jl`); this file handles the two-single-channel-TIFF anchor layout ONLY.

# Names the frozen src/ files resolve at call time (Images qualified; cor/corspearman/corkendall
# unqualified inside `correlation`). Importing them here lets load.jl reach the src/ summary path
# standalone, mirroring runtests.jl's import block.
import Images
import Statistics: mean, cor
import StatsBase: corspearman, corkendall

# Guarded, READ-ONLY reach into the FROZEN src/ (no edits — CLAUDE.md decoupling). Skipped when a
# sibling (e.g. runtests.jl) already pulled these into scope.
isdefined(@__MODULE__, :MultiChannelImage) || include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))
isdefined(@__MODULE__, :correlation)       || include(joinpath(@__DIR__, "..", "src", "colocalization.jl"))

# Guarded include of the corpus contract consts (MIN_PATCH_VALID_PIXELS documents the Pitfall-6 floor).
isdefined(@__MODULE__, :MIN_PATCH_VALID_PIXELS) || include(joinpath(@__DIR__, "config.jl"))

"""
    load_anchor(name, ch1_path, ch2_path; channels=["ch1","ch2"]) -> MultiChannelImage

Load a two-single-channel-TIFF anchor into a `MultiChannelImage` by calling the FROZEN
convenience constructor `MultiChannelImage(name, [ch1_path, ch2_path], channels)`
(`src/LoadImages.jl`), which reads each TIFF via `load_tiff` and auto-derives `pixel_size`
and the per-channel Otsu threshold. The returned image is directly consumable by the
`src/` summary path and by the Phase-9 comparator harness (which accepts
`Vector{MultiChannelImage}`).
"""
function load_anchor(name::AbstractString, ch1_path::AbstractString, ch2_path::AbstractString;
                     channels::Vector{<:AbstractString} = ["ch1", "ch2"])
    return MultiChannelImage(String(name), String[String(ch1_path), String(ch2_path)],
                             String[String(c) for c in channels])
end

"""
    summary_grid(img; num_patches=8) -> Matrix{Union{Float64,Missing}}

Reduce a two-channel `MultiChannelImage` to the `num_patches × num_patches` per-patch Pearson
correlation grid through the UNCHANGED `src/colocalization.jl` summary path: `patch(·, num_patches)`
on each channel then `correlation(·, ·)`. `num_patches` defaults to the spike-fixed 8×8 grid.
Patches with `<= MIN_PATCH_VALID_PIXELS` non-zero valid pixels come back `missing` (the frozen
`_exclude_zero` floor).
"""
function summary_grid(img; num_patches::Integer = 8)
    p1 = patch(img.data[1], num_patches)
    p2 = patch(img.data[2], num_patches)
    return correlation(p1, p2)
end

"""
    has_signal(grid) -> Bool

Pitfall-6 guard: return `true` iff `grid` has at least one non-`missing` entry AND the mean over
its valid (non-`missing`) patches is finite. An anchor whose intensity scaling / background
starves every 8×8 patch below the `MIN_PATCH_VALID_PIXELS` floor reduces to an all-`missing`
grid (or a non-finite mean) and returns `false` — a signal to re-check the real anchor bytes
before trusting the conversion.
"""
function has_signal(grid)
    any(!ismissing, grid) || return false
    return isfinite(mean(skipmissing(grid)))
end
