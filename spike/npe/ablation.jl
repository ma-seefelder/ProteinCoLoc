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

# spike/npe/ablation.jl --- Phase-4 summary-statistic ablation (ABL-01 / ABL-02; D-06 / D-07).
#
# Trains the SAME NPE architecture (build_estimator) on the minimal (:min, 128-dim)
# and augmented (:aug, 142-dim) cached summaries under the loader's leak-free k=5 CV,
# scores per-parameter RMSE (all 7 θ, ρ_true headline) per fold per variant, and
# applies the pre-registered decision rule (D-07) to choose :min or :aug.
#
# ZERO RE-SIMULATION (D-06): both summary variants are already cached (Phase-3 D-02);
# `load_fold(...; variant=:min|:aug)` toggles which cached matrix is standardized. The
# ONLY thing that changes between the two arms is the `variant` argument -- the
# architecture (build_estimator), the master_seed, the fold membership, the θ-transform
# recipe and every training hyper-parameter are held IDENTICAL. That is what makes the
# comparison isolate the SUMMARY choice, not the network (D-06/D-07).
#
# RMSE VIA assess/rmse (RESEARCH "Don't Hand-Roll", mirrors benchmark.jl): per-parameter
# RMSE is taken from NeuralEstimators `rmse(assess(est, θva_std, Zva))` on the held-out
# fold, in STANDARDIZED θ-space, then scaled back to physical units by `θzt.scale`
# (Pitfall 5). Both names are QUALIFIED (`NeuralEstimators.assess` / `.rmse`) because
# `Images`/`ImageQualityIndexes` (in scope via the contract.jl chain under runtests.jl)
# also export them -- a bare call is an ambiguous-binding error (04-05 decision).
#
# PRE-REGISTERED RULE (D-07 / RESEARCH A3): `choose_summary` returns :aug ONLY iff aug
# beats min on ρ_true RMSE by the relative margin ABL_REL_MARGIN (0.05) AND improves in
# ≥ABL_FOLD_CONSISTENCY (4) of the 5 folds; otherwise :min (parsimony + OOD
# detectability). A SBC-pass-but-high-RMSE outcome is an INSUFFICIENCY signal, not a
# pass -- the rule scores accuracy, and the tie defaults to the more parsimonious /
# OOD-orthogonal-detectable minimal summary (couples to Phase 5).
#
# CPU-ONLY (D-10 / Pitfall 1): every training/assess call passes `use_gpu=false`; CUDA
# is never imported.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reaches src/ only transitively
# through train_npe.jl's guarded chain. Guarded includes keep the file loadable
# standalone AND inside runtests.jl.

using NeuralEstimators   # assess, rmse (qualified at call site; name collides with Images)
using StatsBase          # transform (frozen θ-standardization)
using Statistics         # mean

# --- ORDER MATTERS: trainer first (train_fold + build_estimator + the NPE_* consts +
#     load_fold via its guarded chain). Guarded for idempotency (standalone + runtests.jl).
isdefined(@__MODULE__, :train_fold) || include(joinpath(@__DIR__, "train_npe.jl"))

# Pre-registered ablation constants (RESEARCH A3 / test_npe.jl). Guarded so a
# redefinition to the SAME value under runtests.jl (where test_npe.jl declares them
# first) is a silent no-op; standalone `include` still gets the locked defaults.
if !isdefined(@__MODULE__, :ABL_REL_MARGIN)
    const ABL_REL_MARGIN = 0.05        # D-07 (A3): aug wins iff RMSE_aug ≤ (1−δ)·RMSE_min on ρ_true …
end
if !isdefined(@__MODULE__, :ABL_FOLD_CONSISTENCY)
    const ABL_FOLD_CONSISTENCY = 4     # … AND aug improves in ≥4 of 5 folds (else keep :min)
end

# The two ablation arms (D-06). Order is fixed so `rmse_table` slice 1 is always :min.
const ABL_VARIANTS = (:min, :aug)

"""
    _abl_rmse_vector(rmse_df) -> Vector{Float64}

Extract the per-parameter RMSE from a NeuralEstimators `rmse(assessment)` DataFrame
into a 7-vector in θ field order (ρ_true first), keyed by the default parameter names
`θ1..θ7` so a re-ordered DataFrame cannot scramble the parameter axis. Mirrors
benchmark.jl `_rmse_vector` (kept local so ablation.jl does not pull the heavier
benchmark.jl chain -- BenchmarkTools / resimulate). No `using DataFrames` (transitive).
"""
function _abl_rmse_vector(rmse_df)
    names  = rmse_df.parameter
    values = rmse_df.rmse
    lut    = Dict(String(names[i]) => Float64(values[i]) for i in 1:length(values))
    return [lut["θ$i"] for i in 1:7]
end

"""
    _abl_theta_scale(θzt) -> Vector{Float64}

The per-parameter scale of the frozen θ `ZScoreTransform`, used to map a
standardized-space RMSE back to physical units (RMSE_phys = RMSE_std · scale). Mirrors
benchmark.jl `_theta_scale`.
"""
_abl_theta_scale(θzt) = collect(Float64.(θzt.scale))

"""
    fold_rmse(dir, fold, variant; master_seed, K, N, epochs, batchsize,
              use_gpu = false) -> Vector{Float64}

Per-parameter (all 7 θ) validation RMSE for ONE fold of ONE variant, in physical units.

Trains `build_estimator` on the fold's leak-free train tensors via `train_fold`
(variant-toggled, arch-identical), then scores the HELD-OUT fold with NeuralEstimators
`rmse(assess(est, θva_std, Zva))`: `load_fold` is re-invoked (deterministic — same
`master_seed`/`variant` ⇒ byte-identical split) to recover `(Zva, θva)`, θva is
standardized with the SAME frozen `θzt` the trainer fit on θtr only, and the
standardized RMSE is scaled to physical units by `θzt.scale` (Pitfall 5). Row 1 is
ρ_true. `use_gpu=false` on every call (D-10).
"""
function fold_rmse(dir, fold::Integer, variant::Symbol; master_seed,
                   K::Integer = 5, N::Integer = 2000, epochs::Integer = 200,
                   batchsize::Integer = 64, use_gpu::Bool = false)
    use_gpu && throw(ArgumentError("fold_rmse: use_gpu=true is out of scope (D-10 CPU-only gate)"))

    # Train the arm (arch-identical; only `variant` differs → isolates the summary, D-06).
    res = train_fold(dir, fold; master_seed = master_seed, variant = variant,
                     use_gpu = false, epochs = epochs, batchsize = batchsize)

    # Recover the held-out fold tensors (deterministic re-load; same split as the trainer).
    fd = load_fold(dir, fold; K = K, master_seed = master_seed, variant = variant)

    # assess on the validation fold in standardized θ-space, then scale to physical units.
    θva_std = Float32.(StatsBase.transform(res.θzt, fd.θva))
    a       = NeuralEstimators.assess(res.estimator, θva_std, fd.Zva; use_gpu = false, N = N)
    rmse_std = _abl_rmse_vector(NeuralEstimators.rmse(a))
    return rmse_std .* _abl_theta_scale(res.θzt)      # 7-vector, physical units
end

"""
    ablate(dir; master_seed = NPE_MASTER_SEED, K = 5, N = 2000, epochs = 200,
           batchsize = 64) -> NamedTuple

Run the summary-statistic ablation (ABL-01, D-06): for each fold `f ∈ 1:K` and each
variant `∈ (:min, :aug)`, train the SAME `build_estimator` architecture on the fold's
leak-free train tensors (zero re-simulation — only the cached `variant` differs) and
score per-parameter RMSE (all 7 θ) on the held-out fold via `assess`/`rmse`.

Returns a NamedTuple:

  - `rmse_table` : `K×2×7` `Array{Float64,3}` of physical-unit RMSE, indexed
                   `[fold, variant, parameter]`; slice `[:, 1, :]` is `:min`,
                   `[:, 2, :]` is `:aug` (order = `ABL_VARIANTS`); parameter 1 is ρ_true.
  - `rho_rmse_min` / `rho_rmse_aug` : the per-variant MEAN ρ_true RMSE over the K folds
                   (the D-07 headline inputs).
  - `fold_wins_aug` : number of folds in which `:aug` beats `:min` on ρ_true RMSE
                   (the D-07 fold-consistency input).
  - `rho_per_fold` : `K×2` matrix of per-fold ρ_true RMSE (`[:,1]`=min, `[:,2]`=aug).
  - `variants` / `K` : echoed for the summary/decision consumers.

Architecture, `master_seed`, fold membership and every hyper-parameter are held
IDENTICAL across variants; only the cached summary changes (D-06). CPU-only
(`use_gpu=false`, D-10). `epochs`/`batchsize`/`N` are exposed so a small fixture gate
(SC4) runs fast; the reported ablation uses the locked training defaults.
"""
function ablate(dir; master_seed = NPE_MASTER_SEED, K::Integer = 5, N::Integer = 2000,
                epochs::Integer = 200, batchsize::Integer = 64)
    rmse_table   = Array{Float64,3}(undef, K, length(ABL_VARIANTS), 7)
    rho_per_fold = Matrix{Float64}(undef, K, length(ABL_VARIANTS))

    for f in 1:K, (vi, variant) in enumerate(ABL_VARIANTS)
        r = fold_rmse(dir, f, variant; master_seed = master_seed, K = K, N = N,
                      epochs = epochs, batchsize = batchsize, use_gpu = false)
        rmse_table[f, vi, :] = r
        rho_per_fold[f, vi]  = r[1]                    # ρ_true is row 1 (the headline)
    end

    imin, iaug   = 1, 2                                # ABL_VARIANTS = (:min, :aug)
    rho_rmse_min = mean(@view rmse_table[:, imin, 1])
    rho_rmse_aug = mean(@view rmse_table[:, iaug, 1])
    fold_wins_aug = count(f -> rmse_table[f, iaug, 1] < rmse_table[f, imin, 1], 1:K)

    return (rmse_table    = rmse_table,
            rho_rmse_min  = rho_rmse_min,
            rho_rmse_aug  = rho_rmse_aug,
            fold_wins_aug = fold_wins_aug,
            rho_per_fold  = rho_per_fold,
            variants      = ABL_VARIANTS,
            K             = K)
end
