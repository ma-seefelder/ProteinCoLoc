#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/harness.jl --- the per-grid CPU ship-gate harness (D-05).
#
# THE ONE θ*~π → simulate → CPU-infer chain every shipped grid's independent ship-gate rides.
# Promoted (grid-generalized, CPU-forced) from the proven spike harness
# (spike/validation/harness.jl:106-136) with the THREE critical D-05 deltas:
#
#   (1) GRID-PARAMETRIZED: `patch_summary(mci, G)` — the grid `G` is a parameter (the spike
#       hardcoded 8), so a 4×4 / 8×8 / 16×16 / 32×32 gate rides the SAME harness.
#   (2) CPU-ONLY, CPU-RESIDENT NET: every NeuralEstimators call is `use_gpu = false` against a
#       CPU-resident frozen net (T-7-07). The pre-registered gate numbers must be deterministic
#       REGARDLESS of the hardware the net was trained on (D-06 reproducibility split).
#   (3) PER-GRID NET: the net is loaded via `ProteinCoLoc.load_estimator` (NOT the fixed spike
#       `trained_npe.jld2`). `m` is the `(estimator, θzt, zt, …)` bundle it returns.
#
# FROZEN-STATS DISCIPLINE: `m.zt`/`m.θzt` are the deployed frozen transforms; standardization is
# NEVER re-fit here — every summary enters the estimator's input space with the frozen `m.zt` and
# every ρ read un-standardizes with `m.θzt` (Pitfall 5).
#
# SIMULATOR INJECTION: the Phase-2 forward model (`sample_prior`/`simulate_pair`/`build_mci`) is
# promoted into `ProteinCoLoc` by a later per-grid plan. Until then a gate caller INJECTS `sim`
# (mirroring the pipeline's injectable `datagen`); once promoted, `default_simulator()` resolves
# it automatically. This lets the harness + SBC machinery be smoke-tested today without the
# forward simulator.
#
# Standalone-loadable (spike test-file convention): a per-grid gate loads its `gate_consts_<G>.jl`
# BEFORE this file (defining SBC_IMSIZE / SBC_L / prod_rng); a standalone load falls back to the
# committed template. Guarded includes keep it idempotent under `run_gate.jl` and `runtests.jl`.

using ProteinCoLoc
using StatsBase          # reconstruct (inverse ZScoreTransform)
using Statistics         # mean

# Pre-registration consts (SBC_IMSIZE / SBC_L defaults + prod_rng). Guarded: a per-grid consts
# file included first is respected; a standalone load falls back to the template.
isdefined(@__MODULE__, :SBC_IMSIZE) || include(joinpath(@__DIR__, "gate_consts_template.jl"))

"""
    default_simulator() -> NamedTuple

Lazily resolve the promoted `ProteinCoLoc` forward-model chain
`(sample_prior, simulate_pair, build_mci)`. Until a later per-grid plan promotes the simulator
into `src/`, this errors with an instruction to INJECT `sim` — so gate machinery is testable
today (inject a fake simulator) and runs unchanged once the real simulator lands.
"""
function default_simulator()
    for s in (:sample_prior, :simulate_pair, :build_mci)
        isdefined(ProteinCoLoc, s) || error(
            "draw_simulate_infer: ProteinCoLoc.$s is not promoted yet — pass " *
            "`sim = (; sample_prior, simulate_pair, build_mci)`.")
    end
    return (sample_prior  = getfield(ProteinCoLoc, :sample_prior),
            simulate_pair = getfield(ProteinCoLoc, :simulate_pair),
            build_mci     = getfield(ProteinCoLoc, :build_mci))
end

"""
    draw_simulate_infer(m, rng; G, imsize = SBC_IMSIZE, N = SBC_L, sim = default_simulator())

ONE θ*~π → simulate → CPU-infer step for a `G×G` gate. Draws a prior θ*, simulates a channel
pair, runs it through the grid-parametrized `patch_summary(mci, G)` summary, standardizes with the
FROZEN `m.zt`, takes ONE `posterior_for` pass (`use_gpu = false`), and un-standardizes with the
frozen `m.θzt`. Returns `(θ = θ*, Z = standardized-summary, draws = 7×N physical-θ)`; row 1 of
`draws` is ρ_true. CPU-only.
"""
function draw_simulate_infer(m, rng; G::Integer, imsize = SBC_IMSIZE, N = SBC_L,
                             sim = default_simulator())
    θ   = sim.sample_prior(rng)                                             # θ* ~ π(θ)
    mci = sim.build_mci(sim.simulate_pair(rng, θ; imsize = imsize))         # forward sim → MCI
    Z   = ProteinCoLoc.standardize_summary(
              ProteinCoLoc.encode_d01(ProteinCoLoc.patch_summary(mci, G)), m.zt, :min)
    draws_std = ProteinCoLoc.posterior_for(m.estimator, Z; N = N, use_gpu = false)  # CPU gate
    draws     = StatsBase.reconstruct(m.θzt, draws_std)                    # 7×N physical (Pitfall 5)
    return (θ = θ, Z = Z, draws = draws)
end

"""
    draw_simulate_infer_paired(m, rng; G, imsize = SBC_IMSIZE, N = SBC_L, sim = default_simulator())

The paired-draw path for the Δρ SBC and the BF gate (D-01): draw TWO INDEPENDENT priors θ*_s,
θ*_c, simulate/encode/standardize both against the FROZEN `m.zt`, and return the N un-standardized
ρ_true draw vectors for each PLUS the standardized paired summaries `Zs`/`Zc` (the BF gate encodes
`pair_encode(Zs, Zc)`). Returns `(θs, θc, ρs, ρc, Zs, Zc)`. CPU-only.
"""
function draw_simulate_infer_paired(m, rng; G::Integer, imsize = SBC_IMSIZE, N = SBC_L,
                                    sim = default_simulator())
    θs   = sim.sample_prior(rng)
    mcis = sim.build_mci(sim.simulate_pair(rng, θs; imsize = imsize))
    Zs   = ProteinCoLoc.standardize_summary(
               ProteinCoLoc.encode_d01(ProteinCoLoc.patch_summary(mcis, G)), m.zt, :min)

    θc   = sim.sample_prior(rng)
    mcic = sim.build_mci(sim.simulate_pair(rng, θc; imsize = imsize))
    Zc   = ProteinCoLoc.standardize_summary(
               ProteinCoLoc.encode_d01(ProteinCoLoc.patch_summary(mcic, G)), m.zt, :min)

    ρs = ProteinCoLoc.rho_draws(m.estimator, Zs, m.θzt; N = N, use_gpu = false)  # sample ρ_true
    ρc = ProteinCoLoc.rho_draws(m.estimator, Zc, m.θzt; N = N, use_gpu = false)  # control ρ_true
    return (θs = θs, θc = θc, ρs = ρs, ρc = ρc, Zs = Zs, Zc = Zc)
end

"""
    load_gate_model(path) -> NamedTuple

Load a per-grid frozen NPE CPU-resident via `ProteinCoLoc.load_estimator` (NOT the fixed spike
`trained_npe.jld2`). Returns the `(estimator, θzt, zt, arch, meta)` bundle `draw_simulate_infer`
consumes as `m`; the estimator is on CPU, so every gate inference is device-independent.
"""
load_gate_model(path) = ProteinCoLoc.load_estimator(path)
