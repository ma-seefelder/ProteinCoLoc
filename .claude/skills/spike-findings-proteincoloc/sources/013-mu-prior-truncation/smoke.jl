# Spike 013 smoke test: validate symbols, atom-freeness of the truncated prior, datagen shapes,
# and rough per-pair timing. No artifacts written.
using ProteinCoLoc
using Random, Random123, Statistics, Printf
include(joinpath(@__DIR__, "trunc_datagen.jl"))

println("nthreads = ", Threads.nthreads())

# --- 1. atom-freeness: 20k draws of ρ_true from each prior, count exact ±0.99 ------------------
rng = Philox4x(UInt64, (UInt64(0x00A1CE02), UInt64(999)))
Ndraw = 20_000
ρ_trunc = [sample_prior_trunc(rng).ρ_true for _ in 1:Ndraw]
rng2 = Philox4x(UInt64, (UInt64(0x00A1CE02), UInt64(999)))
ρ_base  = [ProteinCoLoc.sample_prior(rng2).ρ_true for _ in 1:Ndraw]
atom_trunc = count(x -> x == -0.99 || x == 0.99, ρ_trunc)
atom_base  = count(x -> x == -0.99 || x == 0.99, ρ_base)
@printf("truncated prior: %d/%d atoms (%.3f%%)  range [%.4f, %.4f]\n",
        atom_trunc, Ndraw, 100*atom_trunc/Ndraw, minimum(ρ_trunc), maximum(ρ_trunc))
@printf("baseline  prior: %d/%d atoms (%.3f%%)  range [%.4f, %.4f]\n",
        atom_base, Ndraw, 100*atom_base/Ndraw, minimum(ρ_base), maximum(ρ_base))
@assert atom_trunc == 0 "truncated prior still produced ρ_true atoms!"

# --- 2. datagen shapes + timing on a tiny batch ------------------------------------------------
module GateV2sm
    include(joinpath(@__DIR__, "..", "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end
IMSIZE_SET     = GateV2sm.SBC_IMSIZE_SET
IMSIZE_WEIGHTS = GateV2sm.SBC_IMSIZE_WEIGHTS
println("imsize mixture: ", IMSIZE_SET, "  w=", IMSIZE_WEIGHTS)

nsm = 16
t0 = time()
out = generate_samples_trunc(nsm; grid = 8, master_seed = UInt64(0x00DEC0DE),
                             imsize_set = IMSIZE_SET, imsize_weights = IMSIZE_WEIGHTS)
dt = time() - t0
@printf("generate_samples_trunc(%d): theta %s  summary_min %s  in %.1fs (%.2fs/pair wall, threaded)\n",
        nsm, size(out.theta), size(out.summary_min), dt, dt/nsm)
@assert size(out.theta) == (7, nsm)
@assert size(out.summary_min) == (ProteinCoLoc.summary_dim(8), nsm)
@assert all(isfinite, out.theta) && all(isfinite, out.summary_min)
println("smoke OK")
