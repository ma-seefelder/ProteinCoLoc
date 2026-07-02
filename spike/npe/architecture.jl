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

# spike/npe/architecture.jl --- Phase-4 NPE architecture builder (NPE-01, D-05).
#
# The single source of the `PosteriorEstimator` topology for the whole phase: the
# trainer (train_npe.jl), the benchmark (04-05) and the ablation (04-06) all build
# their estimator through `build_estimator` so the ARCHITECTURE IS IDENTICAL across
# runs -- only the input width (128 :min / 142 :aug) and the data change. That is
# what makes the ablation a fair, arch-controlled comparison (D-06/D-07).
#
# TOPOLOGY (D-05 -- MLP summary net → NeuralEstimators NormalisingFlow):
#   Chain(Dense(d_in, width, gelu), Dense(width, width, gelu), Dense(width, dstar))
#       → NormalisingFlow(D=7; num_summaries = dstar)
# ρ_true is θ ROW 1 (prior NamedTuple field order), so the 7-marginal flow's first
# marginal is the colocalization knob the phase reports.
#
# API GOTCHAS (preserved from the ONLY verified-working construction, 00_smoke.jl
# L31-37, against installed NeuralEstimators v0.2.1):
#   * `q` is a CONSTRUCTED NormalisingFlow INSTANCE passed POSITIONALLY -- the `q=`
#     keyword convenience form of PosteriorEstimator expects a TYPE, not an instance.
#   * `NormalisingFlow(d; num_summaries = dstar)`.
# Depth/width/dstar are Claude's discretion (RESEARCH A4) and are declared as module
# `const`s so every downstream consumer reuses the SAME defaults.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; touches no src/. CPU-only
# is enforced at the call sites (train_npe.jl / infer.jl `use_gpu=false`); the
# architecture itself is device-agnostic.

using NeuralEstimators   # PosteriorEstimator, NormalisingFlow
using Flux               # Chain, Dense, gelu

# --- Architecture constants (RESEARCH A4; the defaults every wave reuses) ---------
# PHASE-5 CALIBRATION ITERATION (2026-07-02): the frozen Phase-4 topology
# (dstar=32, summary depth=2/width=128, 6 coupling layers) produced OVER-DISPERSED
# (too-wide) ρ_true/Δρ posteriors — SBC rank domes centred at L/2, KS D≈0.15 — that
# failed the Phase-5 gates. The flow lacked capacity to sharpen the conditionals.
# The defaults below are raised so the SAME `build_estimator` yields a higher-capacity
# flow; the extra `num_coupling_layers` / `flow_depth` / `flow_width` knobs feed the
# NeuralEstimators `NormalisingFlow` coupling stack (num_coupling_layers) and the
# per-block conditioner MLPs (depth/width, forwarded as kwargs). ρ_true stays θ ROW 1
# and D stays 7, so downstream index/axis assumptions are unchanged.
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

Build the Phase-4/5 `PosteriorEstimator`: an MLP summary net mapping a `d_in`-vector
standardized loader summary (128 `:min` / 142 `:aug`) to `dstar` learned summaries,
feeding a NeuralEstimators `NormalisingFlow` over the `D`-dim (default 7) simulator
θ. The flow therefore has `D` marginals; marginal 1 is ρ_true.

`q` is constructed as a `NormalisingFlow` INSTANCE and passed POSITIONALLY (the
v0.2.1 API gotcha; see 00_smoke.jl). `num_coupling_layers` sets the flow's coupling
stack depth; `flow_depth`/`flow_width` are forwarded to the per-block conditioner
MLPs (NeuralEstimators `MLP` `depth`/`width` kwargs). All knobs default to the
module `const`s so the ablation and benchmark share the exact same topology.
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
    num_coupling_layers >= 1 || throw(ArgumentError("build_estimator: num_coupling_layers must be ≥ 1"))

    # MLP summary net: Dense(d_in→width, gelu) then (depth-1) width→width gelu
    # blocks, then a linear width→dstar projection to the learned summaries.
    layers = Dense[Dense(d_in, width, gelu)]   # concretely-typed (IN-06), not Any[...]
    for _ in 2:depth
        push!(layers, Dense(width, width, gelu))
    end
    push!(layers, Dense(width, dstar))
    network = Chain(layers...)

    # NormalisingFlow INSTANCE, passed POSITIONALLY. `num_coupling_layers` deepens the
    # coupling stack; `depth`/`width` flow through CouplingLayer→AffineCouplingBlock→MLP
    # to widen the per-block conditioner nets (the capacity that sharpens the flow).
    q = NormalisingFlow(D; num_summaries = dstar,
                        num_coupling_layers = num_coupling_layers,
                        depth = flow_depth, width = flow_width)
    return PosteriorEstimator(network, q)
end
