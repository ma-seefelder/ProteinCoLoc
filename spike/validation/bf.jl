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

# spike/validation/bf.jl --- Phase-5 amortized Bayes factor (BF-01, BF-02).
#
# The amortized log Bayes factor read in ONE forward pass from the model-comparison
# `RatioEstimator` (train_ratio.jl):
#
#     log-BF(Z) = logratio(Z; m=1) − logratio(Z; m=0) − log_prior_odds
#
# The NRE ratio identity `logratio(Z,m) = log[p(m|Z)/p(m)]` makes the difference already
# equal `log(posterior-odds) − log(prior-odds) = log-BF` (the Bayes factor); the MEASURED
# `log_prior_odds` is subtracted per D-07/Pitfall 5 (it is ≈ 0 by the symmetric i.i.d.
# model-index prior, but MEASURED, never assumed). NO quadgk / KDE / shuffle appears on
# this amortized side (BF-02) -- the training-time contrastive shuffle is internal to
# NeuralEstimators and is NOT part of this evaluation path.
#
# D-07 NULL, MIRRORED EXACTLY (src/bayes.jl:109-136, read-only reproduction target):
#   Δρ = μ_sample − μ_control ;  H1 (coloc) = {Δρ > 0} ;  H0 (null) = {Δρ ≤ 0}
#   BF = posterior-odds ÷ prior-odds for {Δρ > 0} (a one-sided interval at threshold 0).
# The model index is defined on the ρ_true contrast at threshold 0, which reproduces the
# baseline {Δμ>0} event because ghat is strictly monotone (train_ratio.jl banner).
#
# WHY src/bayes.jl IS NOT include()d HERE (decoupling-faithful, NOT scope reduction): its
# body relies on the enclosing ProteinCoLoc module's `import KernelDensity: kde` /
# `import QuadGK: quadgk` and DynamicPPL (Turing) — NONE of which are in the lean spike
# env (exactly the reason contract.jl's induced_mu re-implements bayes.jl's mean WITHOUT
# including it). The compute_BayesFactor KDE odds-ratio math (BF-02 reproduction TARGET)
# is therefore reproduced VERBATIM in the ISOLATED spike/baseline/ env by
# run_bf_baseline.jl (path b, 05-PLAN <baseline_kde_note>), which persists
# bf_baseline_artifact.jld2; bf.jl (this file) LOADS that artifact and NEVER touches src/.
#
# APPLES-TO-APPLES (D-07): both sides consume IDENTICAL Δρ inputs — bf.jl generates the
# NPE Δρ posterior draws + the Δρ prior draws for each sweep point and persists them; the
# baseline env computes its KDE log-BF on those SAME persisted draws.
#
# CPU-only (use_gpu=false on every NeuralEstimators call). Guarded includes; flat functions.

using NeuralEstimators   # logratio (the amortized one-pass ratio read)
using StatsBase          # (mean/cor come through Statistics)
using Statistics         # cor, mean, maximum
using JLD2               # persist sweep draws / load the baseline artifact

# --- ORDER MATTERS: the ratio-net surface (load_ratio + assemble/model-index framing)
#     first, then the shared harness (rho_draws / sample_prior / simulate_pair / m.zt).
#     Guarded for idempotency under runtests.jl.
isdefined(@__MODULE__, :load_ratio)          || include(joinpath(@__DIR__, "train_ratio.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))

if !isdefined(@__MODULE__, :BF_POST_DRAWS)
    # NPE Δρ-posterior draws per sweep point (the D-03 MC difference cloud the baseline
    # KDE smooths). NOT an anti-snooping threshold — a Monte-Carlo resolution knob.
    const BF_POST_DRAWS  = 2000
    # Δρ PRIOR draws (constant across the sweep — the baseline's prior-odds sample).
    const BF_PRIOR_DRAWS = 20000

    # Intermediate hand-off (bf.jl → baseline env): the Δρ posterior + prior draws keyed
    # by the sweep targets, so the baseline scores its KDE log-BF on IDENTICAL inputs.
    const BF_SWEEP_DRAWS_PATH = joinpath(@__DIR__, "bf_sweep_draws.jld2")
    # The baseline env's output (baseline env → bf.jl): the KDE log-BF per sweep target.
    const BF_BASELINE_ARTIFACT_PATH =
        joinpath(@__DIR__, "..", "baseline", "bf_baseline_artifact.jld2")
end

"""
    amortized_log_bf(est_bf, Z_pair, log_prior_odds; use_gpu = false) -> Float64

The amortized log Bayes factor for one paired summary `Z_pair` (256-dim vector or 256×1
matrix), in ONE forward pass (BF-01/BF-02):

    logratio(Z; m=1) − logratio(Z; m=0) − log_prior_odds

`grid = [0 1]` evaluates the ratio at the null (m=0) and coloc (m=1) model indices; the
MEASURED `log_prior_odds` is subtracted (Pitfall 5). Uses NeuralEstimators `logratio`
ONLY — no quadgk / KDE / shuffle. CPU-only. Returns a finite scalar.
"""
function amortized_log_bf(est_bf, Z_pair, log_prior_odds; use_gpu::Bool = false)
    use_gpu && throw(ArgumentError("amortized_log_bf: use_gpu=true out of scope (CPU-only gate)"))
    Zc = Z_pair isa AbstractVector ? reshape(Z_pair, :, 1) : Z_pair
    grid = reshape(Float32[0.0 1.0], 1, 2)          # null (m=0), coloc (m=1)
    lr = logratio(est_bf, Zc; grid = grid, use_gpu = false)   # 1×2: [logr(m=0) logr(m=1)]
    return (lr[1, 2] - lr[1, 1]) - log_prior_odds
end

"""
    build_bf_pair(m, rng, Δρ_target; imsize = SBC_IMSIZE, N = BF_POST_DRAWS) -> NamedTuple

Build a (sample, control) pair whose TRUE ρ_true contrast targets `Δρ_target`: draw a
control θ_c from the prior, then a fresh sample θ_s with its ρ_true overridden to
`clamp(θ_c.ρ_true + Δρ_target, −0.99, 0.99)` (fresh nuisances). Simulate/encode/standardize
both with the FROZEN `m.zt`, and return `(Z_pair, Δρ_post, Δρ_true)` where `Δρ_post` is the
NPE Δρ MC-difference posterior draw vector `rho_draws(s) .- rho_draws(c)` (D-03) and
`Δρ_true = θ_s.ρ_true − θ_c.ρ_true`. CPU-only.
"""
function build_bf_pair(m, rng, Δρ_target; imsize = SBC_IMSIZE, N::Integer = BF_POST_DRAWS)
    θc = sample_prior(rng)
    ρs = clamp(θc.ρ_true + Δρ_target, -0.99, 0.99)
    θs = merge(sample_prior(rng), (ρ_true = ρs,))     # target contrast, fresh nuisances

    Zs = standardize_summary(encode_d01(patch_summary(build_mci(simulate_pair(rng, θs; imsize = imsize)))), m.zt, :min)
    Zc = standardize_summary(encode_d01(patch_summary(build_mci(simulate_pair(rng, θc; imsize = imsize)))), m.zt, :min)
    Z_pair = pair_encode(Zs, Zc)   # research A7 difference encoding (concat + contrast); train_ratio.jl

    ρs_draws = rho_draws(m.estimator, Zs, m.θzt; N = N, use_gpu = false)
    ρc_draws = rho_draws(m.estimator, Zc, m.θzt; N = N, use_gpu = false)
    Δρ_post  = ρs_draws .- ρc_draws                   # D-03 MC difference
    return (Z_pair = Z_pair, Δρ_post = Δρ_post, Δρ_true = θs.ρ_true - θc.ρ_true)
end

"""
    delta_rho_prior_draws(rng, K) -> Vector{Float64}

`K` Δρ draws from the PRIOR (ρ_true_s − ρ_true_c for two independent `sample_prior`
draws) — the baseline's prior-odds sample for {Δρ>0}. Constant across the sweep (matches
`compute_BayesFactor`, whose prior term is fixed). No image simulation needed.
"""
function delta_rho_prior_draws(rng, K::Integer)
    d = Vector{Float64}(undef, K)
    for j in 1:K
        d[j] = sample_prior(rng).ρ_true - sample_prior(rng).ρ_true
    end
    return d
end

"""
    generate_bf_sweep(m, est_bf, log_prior_odds; n = BF_SWEEP_N, lo = BF_SWEEP_LO,
                      hi = BF_SWEEP_HI, N = BF_POST_DRAWS, Kprior = BF_PRIOR_DRAWS,
                      rng = val_rng(), imsize = SBC_IMSIZE,
                      draws_path = BF_SWEEP_DRAWS_PATH) -> NamedTuple

Run the held-out Δρ sweep (D-08). For each of `n` target Δρ in `[lo, hi]` it builds a pair
(`build_bf_pair`), reads the amortized log-BF (`amortized_log_bf`, one forward pass), and
collects the NPE Δρ posterior draws. It also draws the Δρ prior sample once. The Δρ
posterior + prior draws are persisted (atomic) to `draws_path` so the ISOLATED baseline env
scores its KDE log-BF on IDENTICAL inputs (apples-to-apples, D-07). Returns
`(targets, logbf_amortized, post_draws, prior_draws)`.
"""
function generate_bf_sweep(m, est_bf, log_prior_odds;
                           n::Integer = BF_SWEEP_N, lo::Real = BF_SWEEP_LO,
                           hi::Real = BF_SWEEP_HI, N::Integer = BF_POST_DRAWS,
                           Kprior::Integer = BF_PRIOR_DRAWS, rng = val_rng(),
                           imsize = SBC_IMSIZE, draws_path = BF_SWEEP_DRAWS_PATH)
    targets = collect(range(lo, hi; length = n))
    logbf_amortized = Vector{Float64}(undef, n)
    post_draws = Vector{Vector{Float64}}(undef, n)
    for (i, Δρt) in enumerate(targets)
        pr = build_bf_pair(m, rng, Δρt; imsize = imsize, N = N)
        logbf_amortized[i] = amortized_log_bf(est_bf, pr.Z_pair, log_prior_odds)
        post_draws[i] = pr.Δρ_post
    end
    prior_draws = delta_rho_prior_draws(rng, Kprior)

    if draws_path !== nothing
        _save_sweep_draws(draws_path, targets, logbf_amortized, post_draws, prior_draws)
    end
    return (targets = targets, logbf_amortized = logbf_amortized,
            post_draws = post_draws, prior_draws = prior_draws)
end

# Atomic persistence of the sweep hand-off file (save_npe idiom).
function _save_sweep_draws(path, targets, logbf_amortized, post_draws, prior_draws)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; targets = targets, logbf_amortized = logbf_amortized,
            post_draws = post_draws, prior_draws = prior_draws,
            generated = string(Dates.now(Dates.UTC)) * "Z")
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "post_draws") "_save_sweep_draws: integrity check failed"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    bf_reproduction(m, est_bf, log_prior_odds; baseline_path = BF_BASELINE_ARTIFACT_PATH,
                    sweep_path = BF_SWEEP_DRAWS_PATH, regenerate = true, kwargs...) -> NamedTuple

The D-08 reproduction gate. Produces the amortized log-BF sweep (regenerating it via
`generate_bf_sweep` when `regenerate=true`, else loading the persisted sweep at
`sweep_path`), loads the KDE baseline log-BF from `baseline_path` (produced by the
ISOLATED baseline env's run_bf_baseline.jl on IDENTICAL Δρ inputs), asserts the two are
keyed by the SAME sweep targets (apples-to-apples), and returns

    (corr = cor(logbf_amortized, logbf_kde), max_abs_err = maximum(abs.(Δ)), targets,
     logbf_amortized, logbf_kde)

The PASS rule (caller-side, D-08): `corr ≥ BF_CORR_MIN` AND `max_abs_err ≤ BF_LOGBF_TOL`
(both criteria — catches wrong ordering AND right-ordering/wrong-magnitude). The amortized
side uses ONLY `logratio` (no quadgk/KDE).
"""
function bf_reproduction(m, est_bf, log_prior_odds;
                         baseline_path = BF_BASELINE_ARTIFACT_PATH,
                         sweep_path = BF_SWEEP_DRAWS_PATH, regenerate::Bool = true,
                         kwargs...)
    if regenerate
        sw = generate_bf_sweep(m, est_bf, log_prior_odds; draws_path = sweep_path, kwargs...)
        targets = sw.targets
        amort   = sw.logbf_amortized
    else
        d = JLD2.load(sweep_path)
        targets = d["targets"]
        amort   = d["logbf_amortized"]
    end

    isfile(baseline_path) ||
        error("bf_reproduction: baseline artifact $baseline_path absent — run " *
              "`julia --project=spike/baseline spike/baseline/run_bf_baseline.jl` first (path b).")
    bl = JLD2.load(baseline_path)
    kde_logbf = bl["logbf_kde"]::AbstractVector
    bl_targets = bl["targets"]

    @assert length(kde_logbf) == length(amort) "bf_reproduction: sweep/baseline length mismatch"
    @assert isapprox(collect(Float64.(bl_targets)), collect(Float64.(targets)); atol = 1e-9) ||
            bl_targets == targets "bf_reproduction: sweep/baseline target mismatch — not apples-to-apples"

    corr = cor(amort, kde_logbf)
    max_abs_err = maximum(abs.(amort .- kde_logbf))
    return (corr = corr, max_abs_err = max_abs_err, targets = targets,
            logbf_amortized = amort, logbf_kde = kde_logbf)
end

"""
    save_bf_figure(rep; filename = "bf_agreement.png") -> String

Write the amortized-vs-KDE log-BF agreement scatter by CALLING figures.jl's
`plot_bf_agreement` (owned by Plan 05-01). Guardedly loads figures.jl (CairoMakie) only
here, so the fast test gate never pays for the plotting stack. NOT called in the gate.
"""
function save_bf_figure(rep; filename::AbstractString = "bf_agreement.png")
    isdefined(@__MODULE__, :plot_bf_agreement) ||
        include(joinpath(@__DIR__, "figures.jl"))
    return plot_bf_agreement(rep.logbf_amortized, rep.logbf_kde; filename = filename)
end
