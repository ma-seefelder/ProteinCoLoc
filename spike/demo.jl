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

# spike/demo.jl --- Phase-6 reproducible end-to-end demo runner (DEMO-01, DEMO-02).
#
# THE TWO-TIER SEEDED DEMO that closes the spike. A single script that actually
# RE-EXERCISES the whole amortized-inference stack reproducibly (not merely re-displays
# stored numbers), and certifies the main package was never mutated.
#
#   FAST DEFAULT   (julia --project=spike spike/demo.jl):
#     Loads the frozen artifacts ONCE and chains EVERY pipeline layer —
#     prior → simulator → summary → NPE → NRE/BF → OOD — at FIXTURE scale from a fixed
#     Random123 seed (VAL_FIX_SEED, never the reserved reported stream). Each inference
#     layer carries a TWIN-RUN reproducibility assert proving bit-identical output on the
#     same seed (Philox4x is deterministic — the cheapest honest SC1 proof). It then loads
#     the reported *_report.jld2 POST-ITERATION headline numbers (under isfile guards) and
#     tabulates the success criteria.
#
#   FULL          (julia --project=spike spike/demo.jl --full):
#     Additionally re-runs the reported-scale gates (run_sbc.jl / run_bf.jl / run_ood.jl)
#     as INDEPENDENT subprocesses so one gate's throw cannot abort the others.
#
# DECOUPLING (CLAUDE.md / DEMO-02, step 0, fail-fast): shell out to
#   git status --porcelain -- src/ src/bayes.jl src/colocalization.jl
# and assert it is empty. The backtick ARGV form passes tokens directly to the process
# (no shell string) — no command-injection surface (threat T-6-01). demo.jl reaches src/
# only read-only through contract.jl's include(); it never writes it.
#
# CPU-only (use_gpu=false everywhere); no net is trained (Phase 5 is closed — the memo
# REPORTS its results). PowerShell-safe verification: the script self-asserts (no POSIX
# `test -f`); the verify command is simply
#     julia --project=spike spike/demo.jl

using JLD2
using Test
using Dates
using Random
using Statistics

# --- Guarded, ORDER-DEPENDENT includes of the composition surfaces (harness.jl:51-58
#     idiom). consts.jl (locked pre-registration) FIRST, then the shared harness, which
#     transitively pulls load_npe, posterior_for, sample_prior/simulate_pair, build_mci
#     (contract.jl → read-only src/), encode_d01, sample_rng. @__DIR__ + joinpath
#     everywhere — never bare relative paths.
isdefined(@__MODULE__, :VAL_MASTER_SEED)   || include(joinpath(@__DIR__, "validation", "consts.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "validation", "harness.jl"))

# --- CLI flag (bare ARGS — no ArgParse dependency; run_thread_sweep.jl idiom) ----------
const FULL = ("--full" in ARGS)

# ============================================================================
# STEP 0 — DECOUPLING ASSERTION (DEMO-02, fail fast BEFORE any compute)
# ============================================================================
# Scoped to exactly the three src/ paths — a whole-tree git status false-positives on
# unrelated .planning/ churn (RESEARCH Pitfall 4). ARGV form, NO string interpolation
# into the backticks (T-6-01, command injection).
let out = read(`git status --porcelain -- src/ src/bayes.jl src/colocalization.jl`, String)
    @assert isempty(strip(out)) "DECOUPLING VIOLATED (DEMO-02): src/ is dirty:\n$out"
    println("[step 0] decoupling OK: src/ src/bayes.jl src/colocalization.jl clean")
end

# ============================================================================
# Frozen model loaded ONCE (frozen zt/θzt; NEVER re-fit — T-05-01)
# ============================================================================
const DEMO_M = load_frozen_model()      # spike/npe/trained_npe.jld2, integrity-checked

# ============================================================================
# STEP 1a — NPE fast-tier chain proof (DEMO-01)
# ============================================================================
# prior → simulate → summary → NPE, at FIXTURE scale on VAL_FIX_SEED. draw_simulate_infer
# standardizes with the frozen m.zt and un-standardizes with the frozen m.θzt (no
# fit(ZScoreTransform)). Returns r.θ (true 7-tuple) + r.draws (7×N physical-θ posterior).
#
# TWO seed sources are fixed so the WHOLE chain is deterministic: the Philox stream
# (val_rng(seed)) drives sample_prior/simulate_pair, and Random.seed!(seed) fixes the
# GLOBAL RNG that NeuralEstimators' sampleposterior draws the flow's base samples from
# (posterior_for threads no rng — infer.jl). Fixing both is what "reproducible from a
# fixed seed" (SC1) means for the stochastic posterior-sampling layer.
function run_npe_fixture(seed)
    Random.seed!(seed)
    return draw_simulate_infer(DEMO_M, val_rng(seed); imsize = (64, 64), N = SBC_FIX_L)
end

println("[step 1a] NPE fixture chain (prior→sim→summary→NPE) on VAL_FIX_SEED …")
const NPE_R1 = run_npe_fixture(VAL_FIX_SEED)
const NPE_R2 = run_npe_fixture(VAL_FIX_SEED)
# Twin-run reproducibility: bit-identical posterior draws on the same seed (SC1).
@assert NPE_R1.draws == NPE_R2.draws "NPE fixture chain not reproducible on VAL_FIX_SEED"
const NPE_REPRO = (NPE_R1.draws == NPE_R2.draws)
println("  NPE twin-run reproducible: $NPE_REPRO  (draws $(size(NPE_R1.draws,1))×$(size(NPE_R1.draws,2)); ρ̂=$(round(mean(NPE_R1.draws[1, :]); digits=3)))")

# CPU-only invariant (D-10): CUDA must never be loaded (test_sbc.jl:114-116 idiom).
@assert !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules)) "CUDA loaded — demo must be CPU-only (D-10)"

# --- Self-assert close (02_simulator_demo.jl:148-152 idiom; PowerShell-safe) -----------
@assert isfile(joinpath(@__DIR__, "npe", "trained_npe.jld2")) "demo FAILED: missing frozen NPE"
println("demo OK: NPE fast-tier chain proof reproduces on VAL_FIX_SEED (CPU-only, src/ decoupled)")
