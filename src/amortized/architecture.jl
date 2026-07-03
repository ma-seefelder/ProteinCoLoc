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
# src/amortized/architecture.jl --- the NPE `PosteriorEstimator` topology (PROD-01).
#
# Promoted (input-width-agnostic, unchanged topology) from the proven spike:
#   spike/npe/architecture.jl:60-66   NPE_* architecture constants (Phase-5 calibrated)
#   spike/npe/architecture.jl:87-115  build_estimator (MLP summary net → NormalisingFlow)
#
# `build_estimator(d_in; ...)` reads its input width from `d_in = size(Ztr, 1)`, so it is
# ALREADY grid-general: a new shipped grid changes `d_in = 2·G²` but NOT this file. The
# Phase-5-calibrated higher-capacity defaults (dstar=64, depth=3/width=256, 10 coupling layers)
# are carried verbatim so the productionized estimator matches the spike-validated topology.
#
# API GOTCHA (preserved from the ONLY verified-working v0.2.1 construction): `q` is a
# CONSTRUCTED `NormalisingFlow` INSTANCE passed POSITIONALLY — the `q=` keyword convenience form
# of `PosteriorEstimator` expects a TYPE, not an instance. Do NOT "simplify" to `q=Type`.
#
# DEVICE-AGNOSTIC: the architecture itself carries no device; `use_gpu` is decided at the
# training call site (train_npe.jl, D-06). Persisted models are always CPU-resident (persist.jl).
#############################################################################################

import Flux: Chain, Dense, gelu
import NeuralEstimators: PosteriorEstimator, NormalisingFlow

# --- Architecture constants (Phase-5 calibrated; every consumer reuses these) ---------------
# The frozen Phase-4 topology (dstar=32, depth=2/width=128, 6 coupling layers) produced
# OVER-DISPERSED ρ_true/Δρ posteriors that failed the Phase-5 SBC gates; the flow lacked
# capacity to sharpen the conditionals. These raised defaults give the SAME `build_estimator`
# a higher-capacity flow. ρ_true stays θ ROW 1 and D stays 7, so downstream index/axis
# assumptions are unchanged.
const NPE_D          = 7    # θ dimension (ρ_true is ROW 1; prior field order)
const NPE_DSTAR      = 64   # learned-summary dimension feeding the flow conditioner
const NPE_DEPTH      = 3    # hidden Dense layers in the MLP summary net
const NPE_WIDTH      = 256  # hidden width of the MLP summary net
const NPE_COUPLING   = 10   # NormalisingFlow coupling layers (was library default 6)
const NPE_FLOW_DEPTH = 2    # conditioner-MLP hidden layers per affine coupling block
const NPE_FLOW_WIDTH = 128  # conditioner-MLP hidden width per affine coupling block

"""
    build_estimator(d_in::Integer, D::Integer = NPE_D;
                    dstar::Integer = NPE_DSTAR, depth::Integer = NPE_DEPTH,
                    width::Integer = NPE_WIDTH,
                    num_coupling_layers::Integer = NPE_COUPLING,
                    flow_depth::Integer = NPE_FLOW_DEPTH,
                    flow_width::Integer = NPE_FLOW_WIDTH) -> PosteriorEstimator

Build the amortized-NPE `PosteriorEstimator`: an MLP summary net mapping a `d_in`-vector
standardized summary (`d_in = 2·G²` for a `G×G` grid) to `dstar` learned summaries, feeding a
NeuralEstimators `NormalisingFlow` over the `D`-dim (default 7) simulator θ. The flow therefore
has `D` marginals; marginal 1 is ρ_true.

`q` is constructed as a `NormalisingFlow` INSTANCE and passed POSITIONALLY (the v0.2.1 API
gotcha). `num_coupling_layers` sets the flow's coupling-stack depth; `flow_depth`/`flow_width`
are forwarded to the per-block conditioner MLPs. All knobs default to the module `const`s.
`build_estimator` is INPUT-WIDTH-AGNOSTIC (reads `d_in`), so it needs no per-grid change.
"""
function build_estimator(d_in::Integer, D::Integer = NPE_D;
                         dstar::Integer = NPE_DSTAR, depth::Integer = NPE_DEPTH,
                         width::Integer = NPE_WIDTH,
                         num_coupling_layers::Integer = NPE_COUPLING,
                         flow_depth::Integer = NPE_FLOW_DEPTH,
                         flow_width::Integer = NPE_FLOW_WIDTH)
    d_in  >= 1 || throw(ArgumentError("build_estimator: d_in must be ≥ 1, got $d_in"))
    D     >= 1 || throw(ArgumentError("build_estimator: D must be ≥ 1, got $D"))
    dstar >= D || throw(ArgumentError("build_estimator: dstar ($dstar) must be ≥ D ($D)"))
    depth >= 1 || throw(ArgumentError("build_estimator: depth must be ≥ 1, got $depth"))
    num_coupling_layers >= 1 ||
        throw(ArgumentError("build_estimator: num_coupling_layers must be ≥ 1"))

    # MLP summary net: Dense(d_in→width, gelu), (depth-1) width→width gelu blocks, then a
    # linear width→dstar projection to the learned summaries. Concretely typed (not Any[...]).
    layers = Dense[Dense(d_in, width, gelu)]
    for _ in 2:depth
        push!(layers, Dense(width, width, gelu))
    end
    push!(layers, Dense(width, dstar))
    network = Chain(layers...)

    # NormalisingFlow INSTANCE, passed POSITIONALLY. `num_coupling_layers` deepens the coupling
    # stack; `depth`/`width` widen the per-block conditioner nets (the capacity that sharpens
    # the flow).
    q = NormalisingFlow(D; num_summaries = dstar,
                        num_coupling_layers = num_coupling_layers,
                        depth = flow_depth, width = flow_width)
    return PosteriorEstimator(network, q)
end
