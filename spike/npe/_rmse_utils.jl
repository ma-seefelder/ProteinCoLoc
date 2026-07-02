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

# spike/npe/_rmse_utils.jl --- shared post-assess helpers (IN-02).
#
# SINGLE definition of the two tiny helpers that both benchmark.jl and ablation.jl
# need after a NeuralEstimators `assess`/`rmse` call. They used to be duplicated
# byte-for-byte (`_rmse_vector`/`_theta_scale` in benchmark.jl vs `_abl_rmse_vector`/
# `_abl_theta_scale` in ablation.jl) -- a drift hazard, since a fix to one would not
# track the other (IN-02). Extracted here so there is one source of truth; both
# consumers pull it in via the existing guarded `include` pattern.
#
# No `using DataFrames` (a transitive dep): `_rmse_vector` uses plain
# `getproperty`/`getindex` on the assessment DataFrame so the spike env's direct-dep
# set is unchanged.

"""
    _rmse_vector(rmse_df) -> Vector{Float64}

Extract the per-parameter RMSE from a NeuralEstimators `rmse(assessment)` DataFrame
into a vector in θ field order (ρ_true first). Keyed by the default parameter names
`θ1..θ7` so a re-ordered DataFrame cannot scramble the parameter axis. Uses plain
`getproperty`/`getindex` on the DataFrame (no `using DataFrames`, a transitive dep).
"""
function _rmse_vector(rmse_df)
    names  = rmse_df.parameter
    values = rmse_df.rmse
    lut    = Dict(String(names[i]) => Float64(values[i]) for i in 1:length(values))
    return [lut["θ$i"] for i in 1:7]
end

"""
    _theta_scale(θzt) -> Vector{Float64}

The per-parameter scale of the frozen θ `ZScoreTransform`, used to map a
standardized-space RMSE back to physical units (RMSE_phys = RMSE_std · scale).
"""
_theta_scale(θzt) = collect(Float64.(θzt.scale))
