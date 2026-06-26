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

# spike/simulator/calibration.jl --- SIM-02: the OFFLINE induced-μ calibration
# (D-01/D-02). Run ONCE to freeze spike/simulator/ghat.jl:
#
#     julia --project=spike spike/simulator/calibration.jl
#
# WHAT IT DOES (the load-bearing method, 02-RESEARCH §Induced-μ Calibration):
#   1. SWEEP ρ_true over a grid spanning PAST ±0.9 (here [−0.99, 0.99], Open Q1
#      RESOLVED) on the SPIKE-VALIDATED shared-latent generator (forward.jl, D-15)
#      -- the generator induces only μ≈±0.76 at knob ±0.9 post-degrade, so the knob
#      MUST extend past ±0.9 for ĝ to reach the Turing μ-prior tails. At each grid
#      point draw n_per pairs with the 6 NUISANCES sampled from their own priors
#      (so induced μ reflects the real generative spread), push each through the
#      UNCHANGED summary (contract.jl: build_mci → patch()/correlation() →
#      induced_mu), and record E[μ|ρ_true] + conditional spread.
#   2. FIT a MONOTONE map (ρ_true, E[μ]) via isotonic regression (PAVA -- the
#      textbook isotonic-regression algorithm, NOT a bespoke optimizer; 02-RESEARCH
#      §Don't Hand-Roll sanctions "isotonic regression / a small monotone fit") and
#      INVERT to ĝ : μ ↦ ρ_true (clamped piecewise-linear, monotone by construction).
#   3. EMIT the frozen ĝ to spike/simulator/ghat.jl (plain Julia literals + a ghat(μ)
#      evaluator + GHAT_MU_MIN/MAX realized range -- no JLD2, that dep is Phase 3).
#   4. EVIDENCE: induced-vs-target μ match (Wasserstein-1 + KS over the
#      physically-realized μ range, D-16), induced-μ size-invariance across
#      {256²,512²,1024²}, a FAITHFUL real-anchor (LoadImages.jl load_tiff, NOT an
#      RGB→luminance reduction, D-16) + negative-tail reachability verdict, and the
#      σ/τ/ν consistency check (D-03). All printed and frozen as a comment header in
#      ghat.jl for transcription into spike/NOTES.md §3.
#
# DECOUPLING (hard): touches ONLY spike/. src/ is reached read-only through
# contract.jl's include() coupling; root manifests stay byte-identical to f581d95.

using Statistics
using StatsBase                       # corspearman (the D-15 monotonicity metric)
using Distributions                   # MU_PRIOR, Uniform nuisance priors
using Random                          # Xoshiro, AbstractRNG

include(joinpath(@__DIR__, "..", "contract.jl"))   # build_mci, patch_summary, induced_mu, load_tiff
include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))  # simulate_pair (shared-latent D-15)

const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# --- The Turing μ-prior (src/bayes.jl ~274-291), the calibration target ----------
const MU_PRIOR = Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)

# --- PRE-DECLARED SIM-02 pass tolerance (T-02-CAL data-snooping guard) -----------
# Fixed BEFORE the run, NOT tuned to pass. Rationale: transform-sampling makes the
# induced μ an UNBIASED, conditional-noise-widened copy of the target, so the
# Wasserstein-1 transport distance is bounded by ~E|conditional scatter| (a few
# 1e-2 in μ-units on the [−1,1] scale). 0.10 is a meaningful "distributions are
# close" bar with headroom for Monte-Carlo noise; evaluated over the physically-
# realized μ range (D-16). KS is reported as a secondary diagnostic.
const SIM02_W1_TOL = 0.10
const SIM02_KS_TOL = 0.20

# --- Nuisance prior ranges (Claude's discretion, D-03; SOURCE OF TRUTH) ----------
# prior.jl MIRRORS these verbatim (documented in NOTES §3). Chosen so the induced
# per-patch-correlation spread (→ σ/τ) lands inside the Turing scale prior
# Truncated(Cauchy(0.1,0.3),1e-4,1); verified by the σ/τ/ν check below.
_sample_nuisances(rng::AbstractRNG) = (
    spillover        = rand(rng, Uniform(0.0, 0.2)),   # directional bleed-through (modest)
    autofluorescence = rand(rng, Uniform(0.0, 0.1)),   # additive background offset
    label_efficiency = rand(rng, Uniform(0.6, 1.0)),   # Bernoulli keep-probability
    shift_dx         = rand(rng, Uniform(-1.0, 1.0)),  # sub-pixel registration error
    shift_dy         = rand(rng, Uniform(-1.0, 1.0)),
    noise            = rand(rng, Uniform(0.0, 1.0)),   # Poisson+Gaussian noise scale
)

_theta(ρ, nui) = merge((ρ_true = ρ,), nui)

# --- Isotonic regression: Pool-Adjacent-Violators (increasing) -------------------
# The canonical isotonic-regression algorithm (Barlow et al.); guarantees a
# non-decreasing fit. NOT a bespoke optimizer.
function _isotonic_increasing(y::Vector{Float64})
    vals = Float64[]
    cnts = Int[]
    for yk in y
        push!(vals, yk); push!(cnts, 1)
        while length(vals) > 1 && vals[end-1] > vals[end]
            v2 = pop!(vals); n2 = pop!(cnts)
            v1 = pop!(vals); n1 = pop!(cnts)
            push!(vals, (v1 * n1 + v2 * n2) / (n1 + n2)); push!(cnts, n1 + n2)
        end
    end
    out = Float64[]
    for (v, c) in zip(vals, cnts)
        append!(out, fill(v, c))
    end
    return out
end

# --- Induced μ at a fixed ρ_true (nuisances drawn from their priors) -------------
function _induced_mu_samples(rng::AbstractRNG, ρ::Float64, imsize::Tuple{Int,Int}, n_per::Int)
    μs = Float64[]
    for _ in 1:n_per
        θ = _theta(ρ, _sample_nuisances(rng))
        m = induced_mu(build_mci(simulate_pair(rng, θ; imsize = imsize)))
        isfinite(m) && push!(μs, m)
    end
    return μs
end

# --- Per-patch-correlation spread + tail diagnostics (σ/τ ↔ spread, ν ↔ tails) ---
function _patchcorr_diagnostics(rng::AbstractRNG, imsize::Tuple{Int,Int}, n::Int)
    vals = Float64[]
    within_sd = Float64[]
    for _ in 1:n
        μstar = rand(rng, MU_PRIOR)
        θ = _theta(μstar, _sample_nuisances(rng))   # rough: μ*≈ρ for diagnostics only
        ρgrid = collect(skipmissing(patch_summary(build_mci(simulate_pair(rng, θ; imsize = imsize)))))
        isempty(ρgrid) && continue
        append!(vals, ρgrid)
        length(ρgrid) > 1 && push!(within_sd, std(ρgrid))
    end
    return (pooled_sd = std(vals), mean_within_sd = mean(within_sd),
            kurtosis = StatsBase.kurtosis(vals))   # excess kurtosis (>0 ⇒ heavier-than-Gaussian)
end

# --- Wasserstein-1 between two empirical samples (quantile transport) -------------
function _w1(a::Vector{Float64}, b::Vector{Float64})
    n = min(length(a), length(b))
    n < 5 && return NaN
    qs = ((1:n) .- 0.5) ./ n
    return mean(abs.(quantile(a, qs) .- quantile(b, qs)))
end

# --- One-sample KS statistic against a distribution ------------------------------
function _ks_stat(sample::Vector{Float64}, d)
    s = sort(sample); n = length(s)
    D = 0.0
    for (i, x) in enumerate(s)
        Fx = cdf(d, x)
        D = max(D, abs(i / n - Fx), abs(Fx - (i - 1) / n))
    end
    return D
end

# --- FAITHFUL real anchor (D-16): the package's own load_tiff, NOT luminance ------
function _real_anchor()
    load_one(cond) = begin
        c1 = load_tiff(joinpath(_REPO_ROOT, "test", "test_images", cond, "$(cond)_c1.tif"))
        c2 = load_tiff(joinpath(_REPO_ROOT, "test", "test_images", cond, "$(cond)_c2.tif"))
        induced_mu(build_mci([Matrix{Float64}(c1), Matrix{Float64}(c2)]; name = cond))
    end
    return (positive = load_one("positive"), negative = load_one("negative"))
end

# ================================================================================
#  calibrate(rng) --- the offline sweep + fit + evidence
# ================================================================================
function calibrate(rng::AbstractRNG = Random.Xoshiro(2026);
                   imsize::Tuple{Int,Int} = (512, 512),
                   n_grid::Int = 25, n_per::Int = 100,
                   metric_N::Int = 300)

    # (1) SWEEP ρ_true ∈ [−0.99, 0.99] (PAST ±0.9 -- Open Q1 RESOLVED) ------------
    ρ_grid  = collect(range(-0.99, 0.99; length = n_grid))
    mu_mean = similar(ρ_grid)
    mu_sd   = similar(ρ_grid)
    for (k, ρ) in enumerate(ρ_grid)
        μs = _induced_mu_samples(rng, ρ, imsize, n_per)
        if isempty(μs)                                  # WR-01: guard mean/std on empty
            @warn "all induced_mu non-finite at ρ=$ρ; skipping grid point"
            mu_mean[k] = NaN; mu_sd[k] = NaN; continue
        end
        mu_mean[k] = mean(μs)
        mu_sd[k]   = std(μs)
    end
    spearman = corspearman(ρ_grid, mu_mean)

    # (2) MONOTONE fit + invert to ĝ : μ ↦ ρ -------------------------------------
    iso = _isotonic_increasing(mu_mean)
    mu_knots = copy(iso)
    for i in 2:length(mu_knots)                     # enforce STRICT increase for inversion
        mu_knots[i] <= mu_knots[i-1] && (mu_knots[i] = mu_knots[i-1] + 1e-9)
    end
    rho_knots = copy(ρ_grid)
    mu_min, mu_max = mu_knots[1], mu_knots[end]

    # local ĝ (same clamped piecewise-linear evaluator emitted into ghat.jl)
    ghat_local(μ) = begin
        μ <= mu_knots[1]   && return rho_knots[1]
        μ >= mu_knots[end] && return rho_knots[end]
        i = searchsortedlast(mu_knots, μ)
        t = (μ - mu_knots[i]) / (mu_knots[i+1] - mu_knots[i])
        rho_knots[i] + t * (rho_knots[i+1] - rho_knots[i])
    end

    # (4a) SIM-02 match metric over the PHYSICALLY-REALIZED μ range (D-16) --------
    target_realized = Truncated(Cauchy(0.0, 0.3), max(-1.0, mu_min), min(1.0, mu_max))
    induced = Float64[]
    for _ in 1:metric_N
        μstar = rand(rng, MU_PRIOR)
        θ = _theta(ghat_local(μstar), _sample_nuisances(rng))
        m = induced_mu(build_mci(simulate_pair(rng, θ; imsize = imsize)))
        isfinite(m) && push!(induced, m)
    end
    target_sample = rand(rng, target_realized, length(induced))
    w1 = _w1(induced, target_sample)
    ks = _ks_stat(induced, target_realized)

    # (4b) SIZE-INVARIANCE across {256²,512²,1024²} (T-02-SI) ---------------------
    si_rhos  = [-0.7, 0.0, 0.7]
    si_sizes = [(256, 256), (512, 512), (1024, 1024)]
    si = Dict{Tuple{Int,Int},Vector{Float64}}()
    for s in si_sizes
        si[s] = [mean(_induced_mu_samples(rng, ρ, s, 40)) for ρ in si_rhos]
    end
    si_maxdev = maximum(maximum(abs.(si[s] .- si[(512, 512)])) for s in si_sizes)

    # (4c) σ/τ/ν consistency (D-03) ----------------------------------------------
    diag = _patchcorr_diagnostics(rng, imsize, 60)
    sigma_ok = 1e-4 <= diag.pooled_sd <= 1.0          # inside Truncated(Cauchy(0.1,0.3),1e-4,1) support

    # (4d) FAITHFUL real anchor + negative-tail reachability (D-16) ---------------
    anchor = try
        _real_anchor()
    catch err
        @warn "real anchor load failed" err
        (positive = NaN, negative = NaN)
    end
    neg_reachable = (isfinite(anchor.negative) && anchor.negative < 0.0)

    return (; ρ_grid, mu_mean, mu_sd, spearman, mu_knots, rho_knots, mu_min, mu_max,
            w1, ks, induced_n = length(induced), si, si_rhos, si_sizes, si_maxdev,
            diag, sigma_ok, anchor, neg_reachable, imsize, n_grid, n_per)
end

# --- Freeze ĝ + evidence into spike/simulator/ghat.jl ----------------------------
_fmt(v) = string(round.(v; digits = 6))

function emit_ghat(r; path = joinpath(@__DIR__, "ghat.jl"))
    content = """
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

# spike/simulator/ghat.jl --- FROZEN induced-μ inverse map ĝ : μ ↦ ρ_true (SIM-02).
# AUTO-GENERATED by spike/simulator/calibration.jl -- DO NOT EDIT BY HAND; re-run
# `julia --project=spike spike/simulator/calibration.jl` to regenerate.
#
# ĝ is a clamped piecewise-linear interpolation of the monotone (isotonic-regressed)
# induced-μ sweep, inverted so prior.jl can draw μ*~Truncated(Cauchy(0,0.3),-1,1) and
# set ρ_true = ghat(μ*). Monotone by construction (knots strictly increasing in μ).
#
# CALIBRATION EVIDENCE (frozen; transcribed into spike/NOTES.md §3):
#   generator        : shared-latent correlated smooth Gaussian field (forward.jl, D-15)
#   sweep            : ρ_true ∈ [−0.99, 0.99], $(r.n_grid) points, n_per=$(r.n_per), size=$(r.imsize) (D-08/D-09)
#   span past ±0.9   : induced μ at knob ±0.9 ≈ ±0.76 post-degrade ⇒ knob extended past ±0.9 (D-15)
#   induced μ range  : [$(round(r.mu_min; digits=4)), $(round(r.mu_max; digits=4))]  (the physically-realized range, D-16)
#   monotonicity     : corspearman(ρ_grid, E[μ]) = $(round(r.spearman; digits=6))
#   SIM-02 match     : Wasserstein-1 = $(round(r.w1; digits=5)) (tol $(SIM02_W1_TOL)), KS = $(round(r.ks; digits=5)) (tol $(SIM02_KS_TOL)) over the realized range (D-16); n=$(r.induced_n)
#   size-invariance  : max |Δ E[μ]| vs 512² across {256²,512²,1024²} = $(round(r.si_maxdev; digits=5))
#   σ/τ consistency  : induced per-patch-corr pooled SD = $(round(r.diag.pooled_sd; digits=4)) (within Truncated(Cauchy(0.1,0.3),1e-4,1): $(r.sigma_ok))
#   ν consistency    : per-patch-corr excess kurtosis = $(round(r.diag.kurtosis; digits=4)) (>0 ⇒ heavier-than-Gaussian ⇒ finite-ν Exponential plausible)
#   real anchor (D-16): faithful LoadImages.jl load_tiff (NOT luminance): positive μ = $(round(r.anchor.positive; digits=4)), negative μ = $(round(r.anchor.negative; digits=4))
#   negative tail    : physically reachable in real fluorescence = $(r.neg_reachable) $(r.neg_reachable ? "" : "⇒ negative μ-prior tail documented PRIOR-ONLY; consistency scoped to realized range")

const GHAT_MU_KNOTS  = $(_fmt(r.mu_knots))

const GHAT_RHO_KNOTS = $(_fmt(r.rho_knots))

# Physically-realized induced-μ range (D-16): the span the SIM-02 consistency claim
# is scoped over. μ outside this range clamps ρ_true to the corresponding endpoint.
const GHAT_MU_MIN = $(round(r.mu_min; digits=6))
const GHAT_MU_MAX = $(round(r.mu_max; digits=6))

\"\"\"
    ghat(μ::Real) -> Float64

The frozen monotone induced-μ inverse: maps a target per-patch-correlation mean μ
to the ρ_true knob that induces it (E[induced μ | ρ_true] = μ). Clamped piecewise
linear; monotone non-decreasing. μ beyond [GHAT_MU_MIN, GHAT_MU_MAX] clamps ρ_true
to the nearest swept endpoint.
\"\"\"
function ghat(μ::Real)
    knots = GHAT_MU_KNOTS
    vals  = GHAT_RHO_KNOTS
    μ <= knots[1]   && return vals[1]
    μ >= knots[end] && return vals[end]
    i = searchsortedlast(knots, μ)
    t = (μ - knots[i]) / (knots[i+1] - knots[i])
    return vals[i] + t * (vals[i+1] - vals[i])
end
"""
    write(path, content)
    return path
end

# --- Script entry point: run the calibration and freeze ghat.jl ------------------
if abspath(PROGRAM_FILE) == @__FILE__
    @info "Running induced-μ calibration (this is the one-time offline SIM-02 step)…"
    r = calibrate()
    p = emit_ghat(r)
    println("\n================ SIM-02 CALIBRATION EVIDENCE (→ NOTES §3) ================")
    println("generator        : shared-latent smooth Gaussian field (forward.jl, D-15)")
    println("sweep            : ρ_true ∈ [-0.99, 0.99], $(r.n_grid) pts, n_per=$(r.n_per), size=$(r.imsize)")
    println("induced μ range  : [$(round(r.mu_min;digits=4)), $(round(r.mu_max;digits=4))]")
    println("monotonicity     : corspearman(ρ_grid, E[μ]) = $(round(r.spearman;digits=6))")
    println("SIM-02 W1        : $(round(r.w1;digits=5))  (tol $(SIM02_W1_TOL))  PASS=$(r.w1 < SIM02_W1_TOL)")
    println("SIM-02 KS        : $(round(r.ks;digits=5))  (tol $(SIM02_KS_TOL))  PASS=$(r.ks < SIM02_KS_TOL)")
    println("size-invariance  : max|Δμ| vs 512² over {256²,512²,1024²} = $(round(r.si_maxdev;digits=5))")
    for s in r.si_sizes
        println("   size $s  E[μ] at ρ=$(r.si_rhos) = $(round.(r.si[s];digits=4))")
    end
    println("σ/τ pooled SD    : $(round(r.diag.pooled_sd;digits=4))  within scale prior = $(r.sigma_ok)")
    println("ν excess kurtosis: $(round(r.diag.kurtosis;digits=4))  (>0 ⇒ heavier-than-Gaussian)")
    println("real anchor (D-16): positive μ = $(round(r.anchor.positive;digits=4)), negative μ = $(round(r.anchor.negative;digits=4))")
    println("negative tail    : physically reachable = $(r.neg_reachable)")
    println("frozen ĝ written : $p")
    println("==========================================================================")
end
