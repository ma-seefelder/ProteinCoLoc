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

# spike/baseline/run_bf_baseline.jl --- Phase-5 KDE Bayes-factor baseline (BF-02, path b).
#
# The reproduction TARGET for the amortized log-BF (bf.jl). This runs in the ISOLATED
# spike/baseline/ env (its own Project/Manifest, quarantined from the lean spike/ env) —
# the ONLY side allowed quadgk/KDE (05-PLAN <baseline_kde_note>).
#
# WHY THE MATH IS PORTED, NOT include()d: src/bayes.jl's compute_BayesFactor relies on the
# enclosing ProteinCoLoc module's `import KernelDensity: kde` / `import QuadGK: quadgk` and
# DynamicPPL/GLMakie; include()-ing bayes.jl standalone pulls that whole heavy stack. Mirroring
# the run_advi.jl precedent (which PORTS the src/bayes.jl @model rather than include()-ing it),
# the three compute_BayesFactor statements (src/bayes.jl:110-135) are reproduced VERBATIM below.
# src/ is NEVER edited (read-only reproduction).
#
# APPLES-TO-APPLES (D-07): this script does NOT re-simulate. It LOADS the Δρ posterior + prior
# draws bf.jl already generated (bf_sweep_draws.jld2, keyed by the sweep targets) so both sides
# score IDENTICAL inputs, and writes the per-target KDE log-BF to bf_baseline_artifact.jld2,
# which bf.jl then loads for the reproduction gate.
#
# Run:  julia --project=spike/baseline spike/baseline/run_bf_baseline.jl

using JLD2            # load sweep draws / persist the baseline artifact
using Statistics      # (reductions)
using Distributions   # pdf(::UnivariateKDE, x)
using KernelDensity   # kde  — faithful to src/bayes.jl:113,122
using QuadGK          # quadgk — faithful to src/bayes.jl:114,123
using Dates

const SWEEP_DRAWS_PATH  = joinpath(@__DIR__, "..", "validation", "bf_sweep_draws.jld2")
const BASELINE_OUT_PATH = joinpath(@__DIR__, "bf_baseline_artifact.jld2")

# src/bayes.jl compute_BayesFactor uses ρ_threshold = 0.0 (D-07, one-sided event {Δρ>0}).
const RHO_THRESHOLD = 0.0

# Numerical safety only: a KDE p never reaches exactly {0,1}, but guard the odds against a
# rounding-to-1 that would produce a non-finite log. Negligible vs. the KDE estimate; the
# three-statement math is otherwise byte-faithful to compute_BayesFactor.
_clampp(p) = clamp(p, 1e-8, 1 - 1e-8)

"""
    p_gt_threshold(draws; threshold = RHO_THRESHOLD) -> Float64

`P(Δρ > threshold)` via the compute_BayesFactor recipe (src/bayes.jl:113-115 / 122-124):
`1 − ∫_{-∞}^{threshold} kde(draws)`. KDE + numerical integration — the exact baseline math.
"""
function p_gt_threshold(draws; threshold::Real = RHO_THRESHOLD)
    dist = kde(collect(float.(draws)))
    p_le, _ = quadgk(x -> pdf(dist, x), -Inf, threshold)
    return _clampp(1 - p_le)
end

"""
    baseline_logbf(post_draws, prior_draws; threshold = RHO_THRESHOLD) -> Vector{Float64}

The KDE log Bayes factor per sweep point (src/bayes.jl:130-135, VERBATIM):
`BF = (p_post/(1−p_post)) / (p_prior/(1−p_prior))`, `logBF = log(BF)`. `p_prior` is fixed
across the sweep (the prior {Δρ>0} sample), `p_post[i]` is the posterior {Δρ>0} at target `i`.
"""
function baseline_logbf(post_draws, prior_draws; threshold::Real = RHO_THRESHOLD)
    p_prior    = p_gt_threshold(prior_draws; threshold = threshold)
    prior_odds = p_prior / (1 - p_prior)
    logbf = Vector{Float64}(undef, length(post_draws))
    for i in eachindex(post_draws)
        p_post        = p_gt_threshold(post_draws[i]; threshold = threshold)
        posterior_odds = p_post / (1 - p_post)
        logbf[i]      = log(posterior_odds / prior_odds)
    end
    return logbf
end

function main()
    isfile(SWEEP_DRAWS_PATH) ||
        error("run_bf_baseline: sweep draws $SWEEP_DRAWS_PATH absent — run generate_bf_sweep " *
              "in the spike env first (bf.jl).")
    d = JLD2.load(SWEEP_DRAWS_PATH)
    targets     = d["targets"]
    post_draws  = d["post_draws"]
    prior_draws = d["prior_draws"]

    logbf = baseline_logbf(post_draws, prior_draws)

    tmp = BASELINE_OUT_PATH * ".tmp"
    jldsave(tmp; targets = targets, logbf_kde = logbf, threshold = RHO_THRESHOLD,
            generated = string(Dates.now(Dates.UTC)) * "Z")
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "logbf_kde") "run_bf_baseline: integrity check failed"
    end
    mv(tmp, BASELINE_OUT_PATH; force = true)
    println("run_bf_baseline: wrote $BASELINE_OUT_PATH ($(length(logbf)) sweep points)")
    return BASELINE_OUT_PATH
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
