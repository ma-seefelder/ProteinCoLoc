# ProteinCoLoc: Bayesian colocalization analysis of multi-channel fluorescence microscopy images
# Copyright (C) 2024  Manuel Seefelder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# spike/validation/run_p11_coverage.jl --- Phase-11: DOES THE NET ALREADY MARGINALISE THE SHIFT?
#
# THE QUESTION THIS SETTLES. `run_p11_hybrid.jl` adds a registration term to the net's Delta-rho
# posterior in quadrature. That is only legitimate if the net does NOT already account for the
# shift. But the training pool draws shift ~ U(-lambda, lambda), so a WELL-CALIBRATED net's
# Delta-rho posterior should ALREADY marginalise over it -- in which case the quadrature step
# double-counts and the hybrid column is an UPPER bound rather than the answer.
#
# The indirect evidence brackets it but does not settle it:
#   net-only width rises        +2.1 % across the ladder
#   exact quadrature total rises +13.4 %
# If the net marginalised perfectly it would already show ~13.4 %. Showing 2.1 % suggests it
# captures roughly a fifth. THE BRACKET IS [1.021, 1.134] AND EVERY CONCLUSION THAT MATTERS IS
# ROBUST ACROSS IT -- both ends are far below the 2.502 bar. But which end is right changes what
# Phase 11 actually CLAIMS, so it is worth settling directly rather than assuming.
#
# THE DIRECT TEST. Calibration is the arbiter. At a FIXED lambda, draw many datasets whose shift
# comes from U(-lambda, lambda), and ask whether the net's Delta-rho posterior covers the TRUE
# Delta-rho at its nominal rate.
#   - If coverage holds at BOTH lambda_min and lambda_max, the net is already marginalising
#     correctly and +2.1 % is approximately the right answer. The hybrid column double-counts.
#   - If coverage DEGRADES as lambda grows (intervals too narrow at lambda_max), the net
#     under-propagates registration uncertainty, and the shortfall is quantified by how far the
#     z-score spread exceeds 1.
# This is a property of the posterior against simulated truth -- no ratio, no bar, no gate.
#
# NOTE ON SCOPE. This measures Delta-rho interval coverage at fixed lambda. It is NOT an SBC run,
# NOT a coverage claim about the shipped bundle, and it gates nothing. It is a diagnostic that
# disambiguates one modelling assumption in a research-lane report.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_coverage.jl

using JLD2
using Dates
using Random
using Statistics

isdefined(@__MODULE__, :LAMBDA_MIN)    || include(joinpath(@__DIR__, "p11_consts.jl"))
isdefined(@__MODULE__, :encode_lambda) || include(joinpath(@__DIR__, "..", "npe", "p11_architecture.jl"))
isdefined(@__MODULE__, :posterior_for) || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))
isdefined(@__MODULE__, :simulate_pair) || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)     || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)    || include(joinpath(@__DIR__, "..", "data", "encode.jl"))

const COVERAGE_REPORT_PATH = joinpath(@__DIR__, "p11_coverage_report.jld2")
const P11_NET_PATH   = joinpath(@__DIR__, "..", "npe", "p11_research_npe.jld2")
const COV_MODEL      = load_npe(P11_NET_PATH)

const COV_N_DATASETS = 300      # per lambda rung
const COV_N_DRAWS    = 2000
const COV_IMSIZE     = P11_PROBE_COMPARABILITY_IMSIZE
const COV_LAMBDAS    = [LAMBDA_MIN, 1.0, 2.0, LAMBDA_MAX]
const COV_NOMINAL    = 0.90

# A key family disjoint from probe (+0), imsize (+1000), bottleneck (+2000) and paired (+4000).
_cov_key(i::Integer) = P11_PROBE_COUNTER * 1000 + 5000 + Int(i)

"""
    theta_at_lambda(rng, lambda) -> NamedTuple

A prior draw whose shift is conditional on a GIVEN lambda -- the generative story of the training
pool with lambda pinned instead of drawn. Every other field, `chromatic_eps` included, comes from
its existing prior exactly as `sample_prior` supplies it (D-09: chromatic_eps is NOT lambda-scaled).

The trailing comma in the merge NamedTuple is load-bearing: without it Julia parses
`(shift_dx = v)` as an assignment expression, `merge` reaches its keyword method and throws at run
time -- the gotcha `run_p11_probe.jl:164-169` documents.
"""
function theta_at_lambda(rng::AbstractRNG, lambda::Real)
    θ  = sample_prior(rng)
    dx = (2 * rand(rng) - 1) * lambda
    dy = (2 * rand(rng) - 1) * lambda
    return merge(θ, (shift_dx = dx, shift_dy = dy,))
end

"90% (or nominal) equal-tailed interval of a draw vector."
function _interval(v::AbstractVector, nominal::Real)
    a = (1 - nominal) / 2
    s = sort(collect(v))
    n = length(s)
    lo = s[max(1, floor(Int, a * n))]
    hi = s[min(n, ceil(Int, (1 - a) * n))]
    return lo, hi
end

function main()
    println("="^78)
    println("PHASE-11 DELTA-RHO COVERAGE AT FIXED LAMBDA")
    println("Does the net ALREADY marginalise the shift? Calibration is the arbiter.")
    println("="^78)
    t0 = time()

    nl = length(COV_LAMBDAS)
    cover = zeros(Float64, nl)
    zsd   = zeros(Float64, nl)
    zmean = zeros(Float64, nl)
    wmean = zeros(Float64, nl)
    rmse  = zeros(Float64, nl)

    for (j, lam) in enumerate(COV_LAMBDAS)
        hits = zeros(Bool,    COV_N_DATASETS)
        zs   = zeros(Float64, COV_N_DATASETS)
        ws   = zeros(Float64, COV_N_DATASETS)
        es   = zeros(Float64, COV_N_DATASETS)

        for i in 1:COV_N_DATASETS
            rng = p11_rng(_cov_key(j * 100_000 + i))
            θs  = theta_at_lambda(rng, lam)
            θc  = theta_at_lambda(rng, lam)
            Zs  = standardize_summary(encode_d01(patch_summary(build_mci(
                      simulate_pair(rng, θs; imsize = COV_IMSIZE)))), COV_MODEL.zt, :min)
            Zc  = standardize_summary(encode_d01(patch_summary(build_mci(
                      simulate_pair(rng, θc; imsize = COV_IMSIZE)))), COV_MODEL.zt, :min)

            Random.seed!(UInt64(1234) + UInt64(i))
            ρs = rho_draws(COV_MODEL.estimator, augment_input(Zs, lam), COV_MODEL.θzt;
                           N = COV_N_DRAWS, use_gpu = false)
            ρc = rho_draws(COV_MODEL.estimator, augment_input(Zc, lam), COV_MODEL.θzt;
                           N = COV_N_DRAWS, use_gpu = false)
            Δ  = ρs .- ρc

            truth = θs.ρ_true - θc.ρ_true
            lo, hi = _interval(Δ, COV_NOMINAL)
            hits[i] = (truth >= lo) && (truth <= hi)
            m, s    = mean(Δ), std(Δ)
            zs[i]   = s > 0 ? (truth - m) / s : NaN
            ws[i]   = s
            es[i]   = truth - m
        end

        cover[j] = mean(hits)
        zsd[j]   = std(filter(isfinite, zs))
        zmean[j] = mean(filter(isfinite, zs))
        wmean[j] = mean(ws)
        rmse[j]  = sqrt(mean(abs2, es))
        println("  lambda = $(rpad(lam, 6))  coverage = $(round(cover[j]; digits = 4)) " *
                "(nominal $(COV_NOMINAL))   z-sd = $(round(zsd[j]; digits = 4))   " *
                "z-mean = $(round(zmean[j]; digits = 4))   post sd = $(round(wmean[j]; digits = 5))" *
                "   rmse = $(round(rmse[j]; digits = 5))")
    end

    println()
    println("-"^78)
    println("READING")
    println("-"^78)
    println("  z-sd near 1.0 at every lambda  => the posterior width is right; the net ALREADY")
    println("     marginalises the shift, and the hybrid quadrature DOUBLE-COUNTS.")
    println("  z-sd growing with lambda       => intervals too narrow as registration worsens;")
    println("     the net UNDER-propagates, and the shortfall is the growth factor.")
    println()
    println("  z-sd ratio lambda_max/lambda_min = ",
            round(zsd[end] / zsd[1]; digits = 4))
    println("  posterior-width ratio            = ", round(wmean[end] / wmean[1]; digits = 4))
    println("  RMSE ratio (the width the data ACTUALLY demands) = ",
            round(rmse[end] / rmse[1]; digits = 4))

    elapsed = (time() - t0) / 60
    tmp = COVERAGE_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version = 1,
        lambdas        = COV_LAMBDAS,
        coverage       = cover,
        nominal        = COV_NOMINAL,
        z_sd           = zsd,
        z_mean         = zmean,
        posterior_sd   = wmean,
        rmse           = rmse,
        zsd_ratio      = zsd[end] / zsd[1],
        width_ratio    = wmean[end] / wmean[1],
        rmse_ratio     = rmse[end] / rmse[1],
        n_datasets     = COV_N_DATASETS,
        n_draws        = COV_N_DRAWS,
        imsize         = COV_IMSIZE,
        elapsed_min    = elapsed,
        generated      = string(Dates.now(Dates.UTC)) * "Z",
        caption        = "Phase-11 Delta-rho interval coverage at fixed lambda, used to decide " *
                         "whether the net already marginalises the registration shift. " *
                         "Diagnostic; not SBC; gates nothing.",
    )
    let d = JLD2.load(tmp); @assert haskey(d, "coverage") && haskey(d, "z_sd"); end
    mv(tmp, COVERAGE_REPORT_PATH; force = true)
    println("\nwrote ", COVERAGE_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    return nothing
end

main()
