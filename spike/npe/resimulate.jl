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

# spike/npe/resimulate.jl --- Phase-4 keyed holdout re-simulation helper
# (RESEARCH Pitfall 4, Open Question 1; D-04).
#
# The Phase-3 cache stores SUMMARIES (θ, encoded vectors, imsize, global_index),
# NOT raw MultiChannelImages -- but `colocalization()` (the ADVI baseline, 04-04)
# needs a raw stack, and the NPE speedup clock (04-05, D-08) must start from the
# SAME raw stack (summary extraction included). This module is the SINGLE shared
# re-simulation path both consume, so the two can never diverge.
#
# REPRODUCIBILITY (Open Question 1): the reserved holdout was written by
# generate.jl::_write_holdout using `holdout_rng(master_seed, j)` -- the disjoint,
# XOR-salted key namespace (D-10) -- for holdout entry j, storing
# `global_index[j] = -j`. Because `Philox4x` is a pure function of its (key, idx)
# and `sample_prior` / `sample_imsize` / `simulate_pair` are deterministic under a
# fixed rng, replaying the exact same stream reproduces the raw stack (and hence
# its summary) BIT-FOR-BIT. The Wave-0 gate in test_npe.jl proves this round-trip;
# if it ever fails, the documented fallback is to persist raw holdout images.
#
# DECOUPLING (hard constraint, CLAUDE.md): reaches src/ ONLY transitively through
# contract.jl's read-only include(); all code lives under spike/. The includes are
# GUARDED (isdefined) so this file loads both standalone and inside runtests.jl
# after the Phase-2/3 chain is already in scope.

# --- ORDER MATTERS: contract.jl FIRST (build_mci/patch_summary + read-only src/
#     coupling), then the simulator prior + forward model, then seeding + encode,
#     then the loader (load_holdout, the sole holdout.jld2 reader). Mirrors
#     generate.jl:42-51.
isdefined(@__MODULE__, :build_mci)     || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :sample_prior)  || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair) || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :holdout_rng)   || include(joinpath(@__DIR__, "..", "data", "seeding.jl"))
isdefined(@__MODULE__, :encode_d01)    || include(joinpath(@__DIR__, "..", "data", "encode.jl"))
isdefined(@__MODULE__, :load_holdout)  || include(joinpath(@__DIR__, "..", "data", "loader.jl"))

"""
    resimulate_holdout(dir, j::Integer; master_seed::Integer) -> NamedTuple

Re-simulate the raw stack for the `j`-th reserved-holdout entry of the cache at
`dir`, keyed EXACTLY as generate.jl::_write_holdout wrote it. Returns
`(mci_sample, theta, global_index)` where:

  - `mci_sample::MultiChannelImage` is the byte-identical re-simulated raw stack,
    the direct input for `colocalization()` (04-04) and the NPE summary clock
    (04-05, D-08).
  - `theta::Vector{Float64}` is the 7-vector in prior field order
    (ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy,
    noise) -- `theta[1]` is ρ_true, the known ground truth for this stack.
  - `global_index::Int` is the stored (negative) holdout index, echoed for
    alignment against `load_holdout(dir).global_index`.

The holdout key is recovered from the STORED index (`global_index[j] = -key`,
D-10) rather than assuming the identity mapping, so the re-sim tracks whatever
holdout ordering the cache actually wrote. Uses `holdout_rng` (NOT `sample_rng`)
so the reserved set is drawn from its provably disjoint stream (D-10).
"""
function resimulate_holdout(dir, j::Integer; master_seed::Integer)
    ho = load_holdout(dir)
    H  = size(ho.theta, 2)
    1 <= j <= H || throw(ArgumentError("resimulate_holdout: j=$j out of 1:$H"))
    # generate.jl stored global_index[j] = -key for the disjoint holdout stream.
    key = -Int(ho.global_index[j])
    rng = holdout_rng(master_seed, key)          # disjoint XOR-salted stream (D-10)
    θ   = sample_prior(rng)                       # i.i.d. π(θ), UNCHANGED prior.jl
    isz = sample_imsize(rng)                      # per-sample image size (D-03)
    mci = build_mci(simulate_pair(rng, θ; imsize = isz))   # UNCHANGED simulator
    return (mci_sample   = mci,
            theta        = collect(values(θ)),    # 7-vector, prior field order
            global_index = Int(ho.global_index[j]))
end

"""
    resimulate_holdout_pair(dir, j_sample::Integer, j_control::Integer;
                            master_seed::Integer) -> NamedTuple

The paired re-simulation form for the Δρ benchmark (D-04): re-simulate two
holdout stacks as a (sample, control) pair and return
`(sample, control, delta_rho_true)` where `sample`/`control` are each a
`resimulate_holdout` NamedTuple and `delta_rho_true = ρ_true(sample) −
ρ_true(control)` is the known ground-truth contrast. Both stacks come from the
SAME re-simulation path so the NPE (difference of two passes) and the ADVI
(joint on the pair) run on identical raw inputs (D-04, apples-to-apples).
"""
function resimulate_holdout_pair(dir, j_sample::Integer, j_control::Integer;
                                 master_seed::Integer)
    s = resimulate_holdout(dir, j_sample;  master_seed = master_seed)
    c = resimulate_holdout(dir, j_control; master_seed = master_seed)
    return (sample = s, control = c,
            delta_rho_true = s.theta[1] - c.theta[1])
end
