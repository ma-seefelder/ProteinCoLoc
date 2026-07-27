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

# spike/validation/run_p11_hybrid.jl --- Phase-11 ANALYTIC UNCERTAINTY PROPAGATION (comparison only).
#
# ############################################################################################
# # THIS IS A HYBRID / ANALYTIC CORRECTION. IT IS **NOT** AMORTIZED INFERENCE.                #
# # IT IS **NOT** AN SC2 RESULT AND IT IS **NOT** A SHIP PATH.                                #
# # It exists to put ONE number beside the network's flat response so the flatness can be     #
# # judged against something, and for no other purpose.                                       #
# ############################################################################################
#
# WHY IT EXISTS. The trained research net's Delta-rho posterior width barely moves with the
# registration-uncertainty level lambda (measured ratios ~1.0). The companion diagnostics
# (`run_p11_bottleneck.jl`, `run_p11_recovery.jl`) show WHY: the 8x8 patch-correlation summary
# does not carry recoverable sub-pixel shift information, so the network has nothing to widen on.
# That is a statement about the summary, not about the physics. The FORWARD MODEL does know how
# much a registration error of a given size perturbs the colocalization estimate, and the D-06
# probe already measured it. This script propagates that known sensitivity analytically and shows
# what the ladder WOULD look like if the uncertainty were accounted for rather than inferred.
#
# WHAT IT IS NOT. Nothing here is learned, calibrated, or validated as a posterior. The inflated
# widths have no coverage guarantee. They are not produced in one forward pass. They cannot be
# reported as SC2, which is a claim about what the AMORTIZED posterior does on its own.
#
# THE CONSTRUCTION, stated plainly so it can be checked.
#
# The companion diagnostics establish that the summary carries NO recoverable information about
# the shift. That result does the work here: it means the shift POSTERIOR simply IS the shift
# PRIOR, `Uniform(-lambda, lambda)` per axis. So there is no posterior to estimate and no
# approximation to defend -- the honest propagation is to MARGINALISE Delta-rho over the shift
# prior:
#
#   sigma_reg(lambda)^2 = E_{dx,dy ~ U(-lambda,lambda)} [ drho_eq( hypot(dx, dy) )^2 ]
#
# where `drho_eq(.)` is the D-06 probe's own MEASURED shift -> Delta-rho-equivalent curve (F5 arm).
# The expectation is taken by Monte Carlo on a fixed, seeded stream.
#
# TWO COEFFICIENT TRAPS, BOTH DELIBERATELY AVOIDED HERE.
# (1) `lambda_ratio = 5.0039` is a RATIO across the ladder, not a magnitude. It CANNOT be used to
#     inflate a posterior width by multiplication. This script therefore uses the probe's ABSOLUTE
#     drho_eq magnitudes and never multiplies any width by 5. `lambda_ratio` appears only as a
#     reference column, never in the arithmetic.
# (2) The registration term must be combined with the existing width in QUADRATURE, not multiplied
#     onto it -- it is an additional independent variance component, not a scale factor.
#
#   sd_hybrid(lambda) = sqrt( sd_net(lambda)^2 + sigma_reg(lambda)^2 )
#
# The quadrature step assumes independence between the network's residual uncertainty and the
# registration-induced perturbation. That is an ASSUMPTION, stated here rather than buried: the two
# are not jointly modelled, which is precisely why this is a correction and not inference.
#
# EXPECT A SMALL NUMBER, AND READ IT AS A RESULT. Given the measured magnitudes it is entirely
# possible this adds only a little to the Delta-rho posterior. That would NOT be a disappointment
# or a failure of the method -- it would be the finding: at an 8x8 summary, registration
# uncertainty is not a major contributor to colocalization uncertainty. Nothing here is tuned to
# make the number larger.
#
# DECOUPLING (S5): spike-local; reads the frozen pre-registration, the committed research net and
# the existing probe artifact; writes one new artifact. Mutates no constant, retrains nothing,
# touches no `src/` file, and does not go near `amended_v2/grid_8`.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_hybrid.jl

using JLD2
using Dates
using Random
using Statistics

# The SAME include chain `spike/test/test_lambda_ablation.jl:61-66` uses, and for the same reason
# it documents: `spike/validation/consts.jl` and `harness.jl` are deliberately NOT included,
# because they redeclare VAL_MASTER_SEED / VAL_FIX_SEED at a different integer width than
# `p11_consts.jl` binds them, and a `const` redeclaration at a different type is a hard error.
# `infer.jl` is what supplies `load_npe` and `rho_draws`.
isdefined(@__MODULE__, :LAMBDA_MIN)       || include(joinpath(@__DIR__, "p11_consts.jl"))
isdefined(@__MODULE__, :encode_lambda)    || include(joinpath(@__DIR__, "..", "npe", "p11_architecture.jl"))
isdefined(@__MODULE__, :posterior_for)    || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))
isdefined(@__MODULE__, :simulate_pair)    || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)        || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)       || include(joinpath(@__DIR__, "..", "data", "encode.jl"))

const HYBRID_REPORT_PATH = joinpath(@__DIR__, "p11_hybrid_report.jld2")
const P11_NET_PATH       = joinpath(@__DIR__, "..", "npe", "p11_research_npe.jld2")
const PROBE_PATH         = joinpath(@__DIR__, "p11_probe_report.jld2")

const HYBRID_N_DATASETS = 12     # more than the tripwire's 3: the ladder is reported, not gated,
                                 # so a tighter estimate is worth the (cheap) extra passes
const HYBRID_N_DRAWS    = 2000   # matches ABLATION_N_DRAWS
const HYBRID_IMSIZE     = P11_PROBE_COMPARABILITY_IMSIZE

const HYBRID_MODEL = load_npe(P11_NET_PATH)

# `sampleposterior` draws from the GLOBAL stream, so the seed is re-pinned before every arm to
# keep the flow's base samples identical across lambda rungs (the tripwire's discipline).
const HYBRID_GLOBAL_SEED = rand(p11_rng(P11_FIXTURE_COUNTER), UInt64)

"""
    _interp(xs, ys, x) -> Float64

Piecewise-linear interpolation on the probe's measured curve, clamped at both ends. Deliberately
NOT a fitted parametric model: the probe's own measured points are the evidence, and fitting a
functional form here would introduce a second, unvalidated modelling step inside a diagnostic.
"""
function _interp(xs::Vector{Float64}, ys::Vector{Float64}, x::Real)
    x <= xs[1]   && return ys[1]
    x >= xs[end] && return ys[end]
    k = findfirst(i -> xs[i] <= x <= xs[i+1], 1:(length(xs)-1))
    t = (x - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] + t * (ys[k+1] - ys[k])
end

"Simulate one (sample, control) pair and freeze both 128-row summaries (tripwire construction)."
function _hybrid_pair(rng)
    θs = sample_prior(rng)
    Zs = standardize_summary(encode_d01(patch_summary(build_mci(
             simulate_pair(rng, θs; imsize = HYBRID_IMSIZE)))), HYBRID_MODEL.zt, :min)
    θc = sample_prior(rng)
    Zc = standardize_summary(encode_d01(patch_summary(build_mci(
             simulate_pair(rng, θc; imsize = HYBRID_IMSIZE)))), HYBRID_MODEL.zt, :min)
    return (Zs = Zs, Zc = Zc)
end

"Δρ posterior draws at registration level `lam` from already-computed 128-row summaries."
function _delta_rho_draws(Zs, Zc, lam; N::Integer = HYBRID_N_DRAWS)
    Random.seed!(HYBRID_GLOBAL_SEED)
    ρs = rho_draws(HYBRID_MODEL.estimator, augment_input(Zs, lam), HYBRID_MODEL.θzt;
                   N = N, use_gpu = false)
    ρc = rho_draws(HYBRID_MODEL.estimator, augment_input(Zc, lam), HYBRID_MODEL.θzt;
                   N = N, use_gpu = false)
    return ρs .- ρc
end

function main()
    println("="^78)
    println("PHASE-11 ANALYTIC UNCERTAINTY PROPAGATION")
    println("HYBRID / analytic correction -- NOT amortized inference, NOT an SC2 result.")
    println("="^78)
    t0 = time()

    probe        = JLD2.load(PROBE_PATH)
    shift_rungs  = collect(Float64.(probe["shift_rungs"]))
    shift_drho   = max.(collect(Float64.(probe["shift_drho_eq_f5"])), 0.0)  # clamp the 0-rung noise
    fwd_ideal    = Float64(probe["lambda_ratio"])                            # 5.003943268109952
    lambdas      = collect(Float64.(SC2_RUNGS))

    println("forward-model ideal ratio (probe lambda_ratio) = ", fwd_ideal)
    println("lambda ladder = ", lambdas)

    # --- Net-only Delta-rho widths per rung --------------------------------------------------
    rng  = p11_rng(P11_FIXTURE_COUNTER)
    pairs = [_hybrid_pair(rng) for _ in 1:HYBRID_N_DATASETS]

    sd_net = zeros(Float64, HYBRID_N_DATASETS, length(lambdas))
    for (i, p) in enumerate(pairs), (j, lam) in enumerate(lambdas)
        sd_net[i, j] = std(_delta_rho_draws(p.Zs, p.Zc, lam))
    end
    net_mean = [mean(view(sd_net, :, j)) for j in 1:length(lambdas)]

    # --- Analytic registration term per rung -------------------------------------------------
    # MARGINALISE over the shift PRIOR, which the recoverability result licenses: the summary
    # carries no information about the shift, so the shift posterior IS U(-lam, lam) per axis.
    #
    #   sigma_reg(lam)^2 = E_{dx,dy ~ U(-lam,lam)} [ drho_eq(hypot(dx,dy))^2 ]
    #
    # Monte Carlo on a fixed seeded stream. This replaces an earlier RMS point-estimate
    # (drho_eq evaluated at |shift|_rms), which is NOT the same thing whenever drho_eq is
    # non-linear in the displacement -- and the measured curve is visibly non-linear near zero.
    # The RMS value is retained below purely as a reported diagnostic, never in the arithmetic.
    n_mc = 200_000
    mc_rng    = p11_rng(P11_PROBE_COUNTER * 1000 + 9000)
    sigma_reg = map(lambdas) do lam
        acc = 0.0
        for _ in 1:n_mc
            dx = (2 * rand(mc_rng) - 1) * lam
            dy = (2 * rand(mc_rng) - 1) * lam
            acc += _interp(shift_rungs, shift_drho, hypot(dx, dy))^2
        end
        return sqrt(acc / n_mc)
    end
    rms_shift = [lam * sqrt(2/3) for lam in lambdas]   # reported only, not used in the arithmetic

    # --- Hybrid inflation --------------------------------------------------------------------
    hybrid = sqrt.(net_mean .^ 2 .+ sigma_reg .^ 2)

    net_ratio    = net_mean ./ net_mean[1]
    hybrid_ratio = hybrid ./ hybrid[1]
    # The forward-model column as a RATIO curve: the probe's own drho_eq shape, normalised at the
    # first rung, whose endpoint ratio is by construction the reported lambda_ratio.
    fwd_curve  = [_interp(shift_rungs, shift_drho, lam) for lam in lambdas]
    fwd_ratio  = fwd_curve ./ fwd_curve[1]

    println()
    println("-"^78)
    println("SIDE BY SIDE  --  Delta-rho posterior width, as a RATIO to the lambda = $(lambdas[1]) rung")
    println("HYBRID COLUMN IS AN ANALYTIC CORRECTION, NOT INFERENCE")
    println("-"^78)
    println(rpad("lambda", 10), rpad("net-only sd", 14), rpad("net ratio", 12),
            rpad("sigma_reg", 12), rpad("hybrid sd", 12), rpad("hybrid ratio", 14),
            rpad("fwd ratio", 12))
    for j in 1:length(lambdas)
        println(rpad(round(lambdas[j]; digits = 3), 10),
                rpad(round(net_mean[j]; digits = 5), 14),
                rpad(round(net_ratio[j]; digits = 4), 12),
                rpad(round(sigma_reg[j]; digits = 5), 12),
                rpad(round(hybrid[j]; digits = 5), 12),
                rpad(round(hybrid_ratio[j]; digits = 4), 14),
                rpad(round(fwd_ratio[j]; digits = 4), 12))
    end
    println()
    println("ENDPOINT RATIOS  (lambda_max / lambda_min)")
    println("  net-only            = ", round(net_ratio[end];    digits = 4))
    println("  analytic propagation= ", round(hybrid_ratio[end]; digits = 4))
    println("  forward-model ideal = ", round(fwd_ideal;         digits = 4),
            "   (probe drho_eq curve endpoint ratio = ", round(fwd_ratio[end]; digits = 4), ")")
    println("  pre-registered SC1g bar P11_LAMBDA_ABLATION_FACTOR = ",
            P11_LAMBDA_ABLATION_FACTOR, "  (NOT applied here -- this script gates nothing)")

    elapsed = (time() - t0) / 60
    tmp = HYBRID_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version        = 1,
        is_amortized_inference = false,
        artifact_kind         = "hybrid_analytic_correction",
        disclaimer            = "HYBRID / analytic correction. NOT amortized inference, NOT an " *
                                "SC2 result, NOT a ship path. No coverage guarantee. The net-only " *
                                "column is the only amortized quantity here.",
        lambda_rungs          = lambdas,
        sd_net_per_dataset    = sd_net,
        sd_net_mean           = net_mean,
        net_ratio             = net_ratio,
        rms_shift             = rms_shift,
        sigma_reg             = sigma_reg,
        sd_hybrid             = hybrid,
        hybrid_ratio          = hybrid_ratio,
        fwd_curve             = fwd_curve,
        fwd_ratio             = fwd_ratio,
        forward_model_ideal   = fwd_ideal,
        preregistered_bar     = P11_LAMBDA_ABLATION_FACTOR,
        n_datasets            = HYBRID_N_DATASETS,
        n_draws               = HYBRID_N_DRAWS,
        imsize                = HYBRID_IMSIZE,
        marginalisation       = "sigma_reg(lambda)^2 = E_{dx,dy~U(-lambda,lambda)}[drho_eq(hypot(dx,dy))^2], Monte Carlo, licensed by the recoverability null: the shift posterior IS the prior.",
        quadrature_assumption = "sd_hybrid = sqrt(sd_net^2 + sigma_reg^2); assumes independence " *
                                "between the network's residual uncertainty and the " *
                                "registration-induced perturbation. Stated, not validated.",
        elapsed_min           = elapsed,
        generated             = string(Dates.now(Dates.UTC)) * "Z",
    )
    let d = JLD2.load(tmp); @assert haskey(d, "hybrid_ratio") && d["is_amortized_inference"] == false; end
    mv(tmp, HYBRID_REPORT_PATH; force = true)
    println()
    println("wrote ", HYBRID_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    return nothing
end

isdefined(@__MODULE__, :P11_HYBRID_LOAD_ONLY) || main()
