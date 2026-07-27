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

# spike/validation/run_p11_bottleneck.jl --- Phase-11 summary-bottleneck diagnosis (SC1g follow-up).
#
# THIS SCRIPT ANSWERS ONE QUESTION: does the fixed 8x8 patch-correlation summary carry sub-pixel
# shift information that ANY estimator could recover, or is the bottleneck upstream of the network?
# It measures the simulator and the summary ONLY. There is NO trained network anywhere in this
# file, so a null result here CANNOT be blamed on training, capacity or optimisation.
#
# THE CENTRAL DISTINCTION THIS SCRIPT EXISTS TO MAKE.
# `run_p11_probe.jl` measures the shift effect in a PAIRED design: the same Philox key is
# re-derived for the perturbed and the reference draw, so stages 1-5 of the forward model are
# bit-identical and the reference displacement at zero injected shift is EXACTLY 0.0. That pairing
# is what let the D-06 probe see a 3 px effect at all, and its header records the reason:
#
#     independent-key ||ds||_2 (same theta, different key) = 2.653  (sd 0.132)
#     |shift| = 3.0 px, PAIRED                             = 1.281  (sd 0.162)
#
# THE NETWORK NEVER GETS THE PAIRED REFERENCE. It is handed ONE image per dataset and must infer
# the shift from that single draw, so the noise it actually faces is the INDEPENDENT-KEY floor,
# not the paired 0.0. This script re-measures that floor from scratch (the value above lives only
# in a source comment, never in an artifact) and puts it side by side with the paired effect curve
# stored in `p11_probe_report.jld2`. The ratio of those two numbers IS the bottleneck answer.
#
# DECOUPLING (S5): spike-local, reads the frozen pre-registration and the existing probe artifact,
# writes one new artifact. It mutates no constant, trains nothing, and touches no `src/` file.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_bottleneck.jl

using JLD2
using Dates
using Statistics

# Load the probe's helpers WITHOUT triggering its reported run. `P11_PROBE_LOAD_ONLY` is the
# guard `run_p11_probe.jl` itself documents for exactly this purpose; defining it here is the
# sanctioned way to reach `_summary`, `_delta_norm`, `_set_shift`, `_set_eps` and `_draw_imsize`
# instead of re-implementing them and risking a metric that silently differs from the probe's.
const P11_PROBE_LOAD_ONLY = true
include(joinpath(@__DIR__, "run_p11_probe.jl"))

const BOTTLENECK_REPORT_PATH = joinpath(@__DIR__, "p11_bottleneck_report.jld2")

# A THIRD key family, disjoint from the probe's two by construction:
#   probe pair key   = P11_PROBE_COUNTER * 1000 +    r      (r <= 999)
#   probe imsize key = P11_PROBE_COUNTER * 1000 + 1000 + r
#   THIS independent  = P11_PROBE_COUNTER * 1000 + 2000 + r
# The +2000 block cannot collide with either for any r <= 999, and P11_PROBE_R = 32.
# It is the SECOND draw of an independent pair: `_pair_key(r)` supplies the first, so the two
# differ only in their noise realisation, which is precisely the quantity being measured.
_indep_key(r::Integer) = P11_PROBE_COUNTER * 1000 + 2000 + Int(r)

"""
    _rehydrate_theta(values, fields) -> NamedTuple

`p11_probe_report.jld2` persists its frozen theta bases as plain `Vector{Float64}` rows plus a
parallel `theta_fields` name vector, not as NamedTuples. The probe's `_set_shift`/`_set_eps`
setters operate on NamedTuples, so the stored rows must be rehydrated before reuse.

Reconstructing from the artifact's OWN field-name vector (rather than re-typing the field order
here) is what keeps this script row-aligned with the probe if the theta arity ever changes.
"""
function _rehydrate_theta(values::AbstractVector{<:Real}, fields::AbstractVector)
    @assert length(values) == length(fields) "theta row and field-name vector must agree"
    return NamedTuple{Tuple(Symbol.(fields))}(Tuple(Float64.(values)))
end

"""
    independent_key_floor(thetas, imsize_of) -> (norms, mean, sd)

The noise floor a SINGLE-DRAW estimator faces: the same theta, the same image size, ZERO injected
shift and zero injected chromatic error, simulated twice under two INDEPENDENT keys. Every
stage-1..7 stochastic element differs, which is exactly the situation of a network shown one image.

Returns the per-(theta, replicate) matrix of `_delta_norm` values plus its mean and sd, on the
same `paired_l2_rows_1_64` metric the probe pre-registered, so the two are directly comparable.
"""
function independent_key_floor(thetas::Vector, imsize_of::Function)
    T = length(thetas)
    R = P11_PROBE_R
    norms = zeros(Float64, T, R)

    units = vec([(t, r) for t in 1:T, r in 1:R])
    t0    = time()
    println("[floor] $(length(units)) independent-key pairs on $(Threads.nthreads()) thread(s)")

    Threads.@threads for u in units
        t, r = u
        sz   = imsize_of(r)
        # Geometry zeroed on BOTH draws: the only difference between them is the key, i.e. the
        # noise realisation. Any displacement measured here is pure Monte-Carlo, by construction.
        theta0 = _set_eps(_set_shift(thetas[t], 0.0), 0.0)
        s_a = _summary(_pair_key(r),  theta0; imsize = sz)
        s_b = _summary(_indep_key(r), theta0; imsize = sz)
        norms[t, r] = _delta_norm(s_b, s_a)
    end

    println("[floor] complete in $(round((time() - t0) / 60; digits = 2)) min")
    return norms, mean(norms), std(norms)
end

"""
    lambda_marginalised_spread(thetas, imsize_of, lambdas) -> Dict

THE CORRECTED FOURTH TEST. Asking whether the summary distribution moves with lambda AT FIXED
THETA is null BY CONSTRUCTION: lambda has no causal path to the image except through the shift,
and theta already contains `shift_dx`/`shift_dy`. That degenerate version is run separately as a
one-line leak check.

The informative version, implemented here: hold theta fixed EXCLUDING the shift, then draw the
shift from `Uniform(-lambda, lambda)` and ask whether the summary DISTRIBUTION widens as lambda
grows. The statistic is the mean displacement from the zero-shift reference at that same key --
i.e. how far a typical draw at this lambda lands from perfect registration.

Compared against the independent-key floor, this says whether the lambda-induced spread is even
visible to a single-draw observer.
"""
function lambda_marginalised_spread(thetas::Vector, imsize_of::Function, lambdas::Vector{Float64})
    T  = length(thetas)
    R  = P11_PROBE_R
    nl = length(lambdas)
    spread = zeros(Float64, T, R, nl)

    units = vec([(t, r) for t in 1:T, r in 1:R])
    t0    = time()
    println("[spread] $(length(units)) units x $nl lambda rungs")

    Threads.@threads for u in units
        t, r   = u
        sz     = imsize_of(r)
        theta0 = _set_eps(_set_shift(thetas[t], 0.0), 0.0)
        s_base = _summary(_pair_key(r), theta0; imsize = sz)
        for (j, lam) in enumerate(lambdas)
            # One shift draw from Uniform(-lam, lam) per (theta, replicate, rung), on a stream
            # keyed off the rung so the draw is reproducible and rung-specific.
            rng = p11_rng(P11_PROBE_COUNTER * 1000 + 3000 + r * 100 + j)
            dx  = (2 * rand(rng) - 1) * lam
            dy  = (2 * rand(rng) - 1) * lam
            s   = _summary(_pair_key(r), merge(theta0, (shift_dx = dx, shift_dy = dy));
                           imsize = sz)
            spread[t, r, j] = _delta_norm(s, s_base)
        end
    end

    println("[spread] complete in $(round((time() - t0) / 60; digits = 2)) min")
    return spread
end

"""
    lambda_leak_check(theta, imsize) -> Float64

The literal fourth test, kept ONLY as a sanity check. At theta fixed INCLUDING the shift, lambda
is not an argument to the forward model at all, so re-simulating on the same key with a different
notional lambda MUST return exactly 0.0. A non-zero value would mean lambda leaks into the
simulator somewhere it should not.
"""
function lambda_leak_check(theta, imsize::Tuple{Int,Int})
    s1 = _summary(_pair_key(1), theta; imsize = imsize)
    s2 = _summary(_pair_key(1), theta; imsize = imsize)
    return _delta_norm(s2, s1)
end

function main()
    println("="^78)
    println("PHASE-11 SUMMARY-BOTTLENECK DIAGNOSIS")
    println("Simulator + summary only. No trained network. Nothing is gated by this script.")
    println("="^78)

    t_start = time()

    # Reuse the EXACT theta bases the D-06 probe froze, so the floor measured here and the effect
    # curve read from the probe artifact are the same five thetas -- not a fresh draw that would
    # make the comparison apples-to-oranges.
    report = JLD2.load(P11_PROBE_REPORT_PATH)
    theta_fields = report["theta_fields"]
    thetas = [_rehydrate_theta(row, theta_fields) for row in report["theta_bases"]]
    println("theta bases reused from the D-06 probe artifact: ", length(thetas))
    println("theta fields: ", theta_fields)
    @assert length(thetas) == P11_PROBE_THETA_BASE_N "theta base count must match the probe"

    # Paired imsize draw, identical convention to the probe's F5 arm: keyed on the replicate, so a
    # rung-to-rung difference can never be an image-size change.
    imsize_of = r -> _draw_imsize(p11_rng(_imsize_key(r)))

    # --- 1. The independent-key noise floor -------------------------------------------------
    floor_norms, floor_mean, floor_sd = independent_key_floor(thetas, imsize_of)

    # --- 2. The paired effect curve, read from the existing artifact ------------------------
    shift_rungs = collect(Float64.(report["shift_rungs"]))
    shift_mean  = collect(Float64.(report["shift_mean_f5"]))
    shift_sd    = collect(Float64.(report["shift_sd_f5"]))

    # --- 3. The corrected lambda-marginalised spread ----------------------------------------
    lambdas = collect(Float64.(SC2_RUNGS))
    spread  = lambda_marginalised_spread(thetas, imsize_of, lambdas)
    spread_mean = [mean(view(spread, :, :, j)) for j in 1:length(lambdas)]
    spread_sd   = [std(view(spread, :, :, j))  for j in 1:length(lambdas)]

    # --- 4. The degenerate leak check -------------------------------------------------------
    leak = lambda_leak_check(thetas[1], imsize_of(1))

    # --- Report -----------------------------------------------------------------------------
    println()
    println("-"^78)
    println("INDEPENDENT-KEY NOISE FLOOR (what a single-draw estimator faces)")
    println("-"^78)
    println("  mean ||ds||_2 = $(round(floor_mean; digits = 4))   sd = $(round(floor_sd; digits = 4))")
    println()
    println("-"^78)
    println("PAIRED SHIFT EFFECT (F5 arm, from the D-06 probe artifact)")
    println("-"^78)
    println(rpad("shift px", 12), rpad("mean ||ds||", 16), rpad("sd", 12),
            rpad("effect/floor", 14))
    for (j, m) in enumerate(shift_rungs)
        ratio = shift_mean[j] / floor_mean
        println(rpad(round(m; digits = 4), 12), rpad(round(shift_mean[j]; digits = 4), 16),
                rpad(round(shift_sd[j]; digits = 4), 12), rpad(round(ratio; digits = 4), 14))
    end
    println()
    println("-"^78)
    println("LAMBDA-MARGINALISED SPREAD (theta fixed EXCLUDING shift; shift ~ U(-lam, lam))")
    println("-"^78)
    println(rpad("lambda", 12), rpad("mean spread", 16), rpad("sd", 12), rpad("spread/floor", 14))
    for (j, lam) in enumerate(lambdas)
        println(rpad(round(lam; digits = 4), 12), rpad(round(spread_mean[j]; digits = 4), 16),
                rpad(round(spread_sd[j]; digits = 4), 12),
                rpad(round(spread_mean[j] / floor_mean; digits = 4), 14))
    end
    println()
    println("lambda leak check (MUST be exactly 0.0): ", leak)
    @assert leak == 0.0 "lambda leaked into the forward model at fixed theta"

    elapsed = (time() - t_start) / 60

    # The headline: the effect at LAMBDA_MAX against the floor a single-draw estimator faces.
    max_effect  = shift_mean[end]
    snr_at_max  = max_effect / floor_mean

    println()
    println("="^78)
    println("HEADLINE  effect at $(shift_rungs[end]) px = $(round(max_effect; digits = 4)), " *
            "independent-key floor = $(round(floor_mean; digits = 4))")
    println("          signal-to-noise for a SINGLE-DRAW estimator = $(round(snr_at_max; digits = 4))")
    println("="^78)

    tmp = BOTTLENECK_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version          = 1,
        floor_norms             = floor_norms,
        floor_mean              = floor_mean,
        floor_sd                = floor_sd,
        shift_rungs             = shift_rungs,
        shift_mean_f5           = shift_mean,
        shift_sd_f5             = shift_sd,
        effect_over_floor       = shift_mean ./ floor_mean,
        lambda_rungs            = lambdas,
        spread_norms            = spread,
        spread_mean             = spread_mean,
        spread_sd               = spread_sd,
        spread_over_floor       = spread_mean ./ floor_mean,
        lambda_leak             = leak,
        snr_at_lambda_max       = snr_at_max,
        theta_base_n            = length(thetas),
        replicates              = P11_PROBE_R,
        metric                  = String(P11_PROBE_METRIC),
        seed                    = P11_DEV_SEED,
        elapsed_min             = elapsed,
        generated               = string(Dates.now(Dates.UTC)) * "Z",
        caption                 = "Phase-11 summary-bottleneck diagnosis: the independent-key " *
                                  "Monte-Carlo floor a single-draw estimator faces, against the " *
                                  "paired shift effect from the D-06 probe. Simulator and summary " *
                                  "only; no trained network; gates nothing.",
    )
    let d = JLD2.load(tmp)
        @assert haskey(d, "floor_mean") && haskey(d, "snr_at_lambda_max")
    end
    mv(tmp, BOTTLENECK_REPORT_PATH; force = true)
    println("wrote ", BOTTLENECK_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")

    return nothing
end

isdefined(@__MODULE__, :P11_BOTTLENECK_LOAD_ONLY) || main()
