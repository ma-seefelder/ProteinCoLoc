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

# spike/npe/train_p11_npe.jl --- Phase-11 research-net trainer (d_in = 129, D = 8). RESEARCH LANE ONLY.
#
# THIS NET IS **NOT** THE SHIPPED NET. IT IS **NOT** A DISTRIBUTABLE ARTIFACT. IT IS **NOT**
# DROP-IN COMPARABLE TO THE SHIPPED 128-ROW READ SURFACE. The shipped `EstimatorBundle` reads a
# 128-row input and a 7-marginal flow; this one reads 129 rows and returns 8 marginals, so no
# number produced through it may be quoted without naming the RESEARCH net explicitly. The
# shipped `amended_v2/grid_8` bundle, its `Artifacts.toml` pin and the Phase-7 GO stay untouched
# (D-16: no public API change, no artifact change, no gate reopened).
#
# IT CARRIES EXACTLY THE FOUR DEVIATIONS ENUMERATED IN `P11_RESEARCH_NET_DEVIATIONS`
# (spike/validation/p11_consts.jl, Tier 1, frozen before anything ran):
#   (1) chromatic_eps as the 8th theta column (D-09), so D = 8 rather than 7;
#   (2) SHIFT_PRIOR widened to Uniform(-3, 3) px -- SPIKE ONLY, deliberately NOT mirrored into
#       src/ (D-02);
#   (3) lambda appended as a 129th input row, so d_in = 129 rather than 128 (D-03);
#   (4) the BoundedThetaTransform ported from src/amortized/architecture.jl (the F2 remedy).
# ANYTHING ELSE THAT CHANGES IS AN UNDECLARED DEVIATION AND MUST BE REPORTED AS ONE.
#
# IT DELIBERATELY DOES **NOT** ROUTE THROUGH `_train_grid_pipeline` (src/amortized/pipeline.jl).
# That pipeline cannot express this net: `src/amortized/train_npe.jl:136` calls
# `build_estimator(d_in; ...)` with NO `D` and `train_npe` exposes no `D` keyword, so `D` is
# pinned to `NPE_D = 7`; and its datagen seam allocates a 7-row theta buffer. A spike-local
# trainer is what D-01 wants in any case -- the research lane must not reach into the shipped
# production path.
#
# CPU-only (D-10). Run:
#     julia --project=spike -t auto spike/npe/train_p11_npe.jl
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reaches src/ only READ-ONLY, transitively
# through the datagen module's frozen summary contract. Guarded includes keep the file loadable
# standalone AND inside a test harness.

using NeuralEstimators   # train
using Flux               # Flux.Optimisers.AdamW
using StatsBase          # ZScoreTransform, fit, transform
using JLD2               # atomic model persistence
using Dates              # UTC-labeled artifact timestamp (IN-07)
using Statistics         # mean

# --- ORDER MATTERS: the Tier-1/Tier-2 pre-registration FIRST (it binds NPE_MASTER_SEED and its
#     siblings at UInt64 width), then the research model surface (the ported theta transform,
#     the lambda encoder, the order-enforcing augmenter, build_p11_estimator), then the
#     estimator builder's own consts, then the leak-free loader (for its authoritative row
#     partition and its deterministic fold split), then the lambda-hierarchical pool.
#     Guarded for idempotency.
isdefined(@__MODULE__, :LAMBDA_MIN)         || include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :build_p11_estimator)|| include(joinpath(@__DIR__, "p11_architecture.jl"))
isdefined(@__MODULE__, :build_estimator)    || include(joinpath(@__DIR__, "architecture.jl"))
isdefined(@__MODULE__, :load_fold)          || include(joinpath(@__DIR__, "..", "data", "loader.jl"))
isdefined(@__MODULE__, :generate_p11_pool)  || include(joinpath(@__DIR__, "..", "data", "p11_generate.jl"))

"The trained Phase-11 research net. This exact path is what `test_lambda_ablation.jl` reads."
const P11_RESEARCH_NET_PATH = joinpath(@__DIR__, "p11_research_npe.jld2")

"Persistence schema for the research net (bump on any breaking layout change)."
const P11_NET_SCHEMA = 1

"Validation fraction of the deterministic k-fold split reused as the train/val split."
const P11_TRAIN_K    = 5   # fold 1 is validation (20%), folds 2..5 are training (80%)
const P11_TRAIN_FOLD = 1

# =============================================================================================
# 1. Frozen train-only transforms (leak-free), then the 129-row assembly
# =============================================================================================

"""
    fit_p11_summary_transform(Z128_train) -> ZScoreTransform

Fit the FROZEN 128-row summary transform on the TRAINING columns ONLY, exactly as
`Loader.load_fold` does: the continuous rows (1:64) are z-scored and the binary present/absent
mask rows (65:128) are BYPASSED. The row split is taken from the authoritative
`Loader._row_partition(:min, 128)` -- it is not restated here, and it must not be relaxed to
accept 129 rows (that error is the loud half of the ordering contract).

Because the fit never sees the validation columns, no validation statistic can leak into
training -- by construction, not by convention.
"""
function fit_p11_summary_transform(Z128tr::AbstractMatrix)
    size(Z128tr, 1) == 128 || throw(DimensionMismatch(
        "fit_p11_summary_transform: expected 128 raw summary rows, got $(size(Z128tr, 1))"))
    cont, _ = Loader._row_partition(:min, 128)
    zt = StatsBase.fit(StatsBase.ZScoreTransform, Z128tr[cont, :]; dims = 2)
    # ZERO-VARIANCE GUARD (loader.jl:191-193): a constant continuous feature has std 0, so the
    # z-score would divide by zero and poison that feature for every consumer. Map it to its
    # centered value by setting the scale to 1. Still leak-free: the fit saw only train columns.
    @inbounds for i in eachindex(zt.scale)
        (isfinite(zt.scale[i]) && zt.scale[i] > 0) || (zt.scale[i] = one(eltype(zt.scale)))
    end
    return zt
end

"""
    apply_p11_summary_transform(Z128, zt) -> Matrix{Float64}

Apply the frozen `zt` to a RAW 128-row summary matrix: continuous rows z-scored, mask rows
passed through UNCHANGED. The same operation `standardize_summary` performs at read time, so
train-time and read-time inputs live in exactly the same space.
"""
function apply_p11_summary_transform(Z128::AbstractMatrix, zt)
    cont, mask = Loader._row_partition(:min, size(Z128, 1))
    out = Matrix{Float64}(undef, size(Z128))
    out[cont, :] = StatsBase.transform(zt, Z128[cont, :])
    out[mask, :] = Z128[mask, :]     # binary mask rows are never re-z-scored (Anti-Pattern)
    return out
end

"""
    build_p11_inputs(Z128_std, lambda) -> Matrix{Float32}

Assemble the 129-row network inputs from the ALREADY-STANDARDIZED 128 rows and each sample's
OWN lambda, one column at a time through `augment_input`.

THE ORDER IS A CORRECTNESS CONTRACT, NOT A STYLE PREFERENCE (Pitfall 2): standardize the 128
rows FIRST, append lambda SECOND. `augment_input` asserts its input is exactly 128 rows, so a
double-append or a pre-standardization append fails loudly here rather than silently z-scoring
the user-declared conditioning variable against the training pool's lambda distribution --
which would couple a declared input to the training data and destroy the read-at-a-chosen-lambda
sweep that is the entire point of D-03.
"""
function build_p11_inputs(Z128_std::AbstractMatrix, lambda::AbstractVector)
    size(Z128_std, 1) == 128 || throw(DimensionMismatch(
        "build_p11_inputs: expected 128 STANDARDIZED rows, got $(size(Z128_std, 1))"))
    size(Z128_std, 2) == length(lambda) || throw(DimensionMismatch(
        "build_p11_inputs: $(size(Z128_std, 2)) columns vs $(length(lambda)) lambdas"))
    out = Matrix{Float32}(undef, 129, size(Z128_std, 2))
    @inbounds for j in axes(Z128_std, 2)
        out[:, j] = vec(augment_input(Z128_std[:, j], lambda[j]))
    end
    return out
end

# =============================================================================================
# 2. Atomic persistence (the save_npe idiom) and the paired reader
# =============================================================================================

"UTC, explicitly Z-labeled (IN-07): an unlabeled local timestamp was an audit defect."
_p11_now() = string(Dates.now(Dates.UTC)) * "Z"

"""
    save_p11_npe(path, result; ...) -> String

Atomically persist the trained research net: `jldsave` to a `.tmp`, reopen to integrity-check,
then `mv(...; force = true)` (a filesystem-atomic rename). A crash before the `mv` leaves only a
discardable `.tmp`, never a half-written model (T-11-30).

The key set is a SUPERSET of what `save_npe` writes, so the artifact is readable by the existing
`load_npe` (which the SC1g tripwire uses) while additionally carrying everything the Phase-11
evaluation plans need: `D`, the lambda encoding constants (so the encoder is reconstructible
without re-including the pre-registration), the four declared deviations, the seed and salt, the
risk trajectories, `elapsed_min`, and the REALIZED training image-size distribution.

`training_imsize_provenance` records realized COUNTS per size, not the intended weights: F6
established that an unrecorded training image-size distribution is a provenance defect, and the
weights a run WOULD use today are not evidence of what a given pool actually contained
(the `src/amortized/persist.jl:143-158` discipline, mirrored here).
"""
function save_p11_npe(path, result::NamedTuple; pool_dir, n_pairs, n_train, n_val,
                      elapsed_min, train_risk, val_risk, imsize_counts, imsize_meta)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp;
        schema_version  = P11_NET_SCHEMA,
        estimator       = result.estimator,
        theta_transform = result.θzt,       # the PORTED bounded transform (deviation 4)
        zt              = result.zt,        # the frozen 128-row summary transform
        variant         = :min,
        d_in            = result.d_in,      # 129 (deviation 3)
        D               = result.D,         # 8   (deviation 1)
        LAMBDA_MIN      = LAMBDA_MIN,
        LAMBDA_MAX      = LAMBDA_MAX,
        P11_RESEARCH_NET_DEVIATIONS = P11_RESEARCH_NET_DEVIATIONS,
        master_seed     = UInt64(P11_DEV_SEED),
        datagen_salt    = UInt64(P11_DATAGEN_SALT),
        train_risk      = train_risk,
        val_risk        = val_risk,
        elapsed_min     = elapsed_min,
        generated       = _p11_now(),
        training_imsize_provenance = imsize_counts,
        training_imsize_meta       = imsize_meta,
        meta = (pool_dir = string(pool_dir), n_pairs = n_pairs,
                n_train = n_train, n_val = n_val,
                dstar = NPE_DSTAR, depth = NPE_DEPTH, width = NPE_WIDTH,
                num_coupling_layers = NPE_COUPLING,
                flow_depth = NPE_FLOW_DEPTH, flow_width = NPE_FLOW_WIDTH,
                research_lane = true, shipped = false,
                generated = _p11_now()))
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "estimator") "save_p11_npe: integrity check failed, $tmp missing estimator"
        @assert haskey(f, "training_imsize_provenance") "save_p11_npe: missing provenance, $tmp"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end

"""
    load_p11_npe(path = P11_RESEARCH_NET_PATH) -> NamedTuple

The read side of the research net. EVERY Phase-11 evaluation script (plans 11-08 through 11-10)
loads the net through this one function, so the estimator, both frozen transforms and the lambda
constants can never be assembled inconsistently at two call sites.

Returns `(estimator, θzt, zt, variant, d_in, D, lambda_min, lambda_max, deviations, meta,
train_risk, val_risk, elapsed_min, generated, training_imsize_provenance)`.
"""
function load_p11_npe(path = P11_RESEARCH_NET_PATH)
    isfile(path) || error("load_p11_npe: no research net at $path (produced by plan 11-07)")
    d = JLD2.load(path)
    haskey(d, "estimator") || error("load_p11_npe: $path missing estimator key")
    return (estimator   = d["estimator"],
            θzt         = d["theta_transform"],
            zt          = d["zt"],
            variant     = d["variant"],
            d_in        = d["d_in"],
            D           = d["D"],
            lambda_min  = d["LAMBDA_MIN"],
            lambda_max  = d["LAMBDA_MAX"],
            deviations  = d["P11_RESEARCH_NET_DEVIATIONS"],
            train_risk  = get(d, "train_risk", nothing),
            val_risk    = get(d, "val_risk", nothing),
            elapsed_min = get(d, "elapsed_min", nothing),
            generated   = get(d, "generated", nothing),
            training_imsize_provenance = get(d, "training_imsize_provenance", nothing),
            meta        = get(d, "meta", nothing))
end

# =============================================================================================
# 3. The training run
# =============================================================================================

"""
    _read_risk_trajectory(logdir) -> (train::Vector{Float64}, val::Vector{Float64})

Read back the per-epoch risk trajectory NeuralEstimators writes to `loss_per_epoch.csv`
(training risk in column 1, validation risk in column 2). Parsed with `Base` only -- no
dependency is added for this (D-01). Returns empty vectors if the file is absent or unreadable;
the trajectory is provenance, never a gate, so its absence must not fail a 1.5 h run.
"""
function _read_risk_trajectory(logdir)
    path = joinpath(logdir, "loss_per_epoch.csv")
    isfile(path) || return (Float64[], Float64[])
    try
        tr = Float64[]; va = Float64[]
        for line in eachline(path)
            parts = split(strip(line), ',')
            length(parts) >= 2 || continue
            a = tryparse(Float64, parts[1]); b = tryparse(Float64, parts[2])
            (a === nothing || b === nothing) && continue
            push!(tr, a); push!(va, b)
        end
        return (tr, va)
    catch e
        e isa InterruptException && rethrow()
        return (Float64[], Float64[])
    end
end

"""
    _stall_diagnostic(train_risk, val_risk; epochs)

Print the risk R11 banner if the run shows the documented stall signature: validation risk not
improved by epoch 50, or early stopping firing before epoch 80.

THIS IS A DIAGNOSTIC, NOT A CAPACITY LEVER, AND IT CHANGES NOTHING AUTOMATICALLY. Spike 012 is
PARTIAL: raising capacity REDISTRIBUTES drift rather than reducing it, and this architecture is
already near its useful capacity ceiling. The sanctioned response to a stall is to REDUCE
`LAMBDA_MAX` and REPORT THE REDUCED LADDER as a declared deviation.
"""
function _stall_diagnostic(train_risk, val_risk; epochs::Integer)
    isempty(val_risk) && return false
    n_epochs   = length(val_risk) - 1          # row 1 is the pre-training validation risk
    improved50 = length(val_risk) > 50 ? minimum(val_risk[1:51]) < val_risk[1] : true
    stalled    = (!improved50) || (n_epochs < 80 && n_epochs < epochs)
    if stalled
        println("!"^78)
        println("RISK R11 — TRAINING STALL SIGNATURE DETECTED")
        println("!"^78)
        println("  epochs run                : $n_epochs (of a budget of $epochs)")
        println("  validation risk improved  : $(improved50 ? "yes" : "NO") by epoch 50")
        println("  initial / best validation : $(val_risk[1]) / $(minimum(val_risk))")
        println()
        println("  SANCTIONED RESPONSE: REDUCE `LAMBDA_MAX` AND REPORT THE REDUCED LADDER.")
        println("  RAISING CAPACITY IS CONTRAINDICATED — spike 012 is PARTIAL: capacity")
        println("  REDISTRIBUTES drift rather than reducing it, and this architecture is already")
        println("  near its useful capacity ceiling. Any reduction must be recorded as a declared")
        println("  deviation in the plan SUMMARY and in the Phase-11 report.")
        println("!"^78)
    end
    return stalled
end

"""
    train_p11_research_net(; n_pairs, savepath, epochs, batchsize, learning_rate,
                           weight_decay, stopping_epochs, use_gpu, verbose) -> NamedTuple

Train the Phase-11 research net on the lambda-hierarchical pool and persist it atomically.

Control flow: resolve the pool -> deterministic fold split -> fit the frozen 128-row summary
transform on TRAIN columns only -> apply to both splits -> append each sample's own lambda as a
129th row -> fit the ported bounded theta transform on TRAIN theta only -> apply to both ->
build a `d_in = 129, D = 8` estimator at UNCHANGED capacity -> train with the unchanged
Phase-5-calibrated recipe -> persist.

CPU-ONLY (D-10 / Pitfall 1). `use_gpu` DEFAULTS TO TRUE in NeuralEstimators, so `use_gpu = false`
is passed on EVERY call and `use_gpu = true` is a hard error here. Datagen dominates this phase's
budget (about 1.2 h of about 1.5 h) and is CPU-bound regardless, so a GPU would buy almost
nothing -- the CLAUDE.md CPU-reproducible baseline stands.
"""
function train_p11_research_net(; n_pairs::Integer = P11_N_PAIRS,
                                savepath = P11_RESEARCH_NET_PATH,
                                cache_root = P11_CACHE_ROOT,
                                epochs = 300,
                                batchsize = 128,
                                learning_rate = 2.5e-4,
                                weight_decay = 1e-4,
                                stopping_epochs = 40,
                                use_gpu::Bool = false,
                                verbose::Bool = true)
    use_gpu && throw(ArgumentError(
        "train_p11_research_net: use_gpu=true is out of scope (D-10 CPU-only gate)"))

    t0  = time()
    dir = p11_pool_dir(n_pairs; cache_root = cache_root)
    verbose && println("="^78)
    verbose && println("Phase-11 RESEARCH-NET training (d_in = 129, D = 8) — RESEARCH LANE ONLY")
    verbose && println("  pool     : $dir")
    verbose && println("  threads  : $(Threads.nthreads())   CPU-only (D-10)")
    verbose && println("  recipe   : epochs=$epochs batch=$batchsize lr=$learning_rate " *
                       "wd=$weight_decay patience=$stopping_epochs (Phase-5 calibrated, unchanged)")
    verbose && println("  deviations declared: $(length(P11_RESEARCH_NET_DEVIATIONS))")
    verbose && println("="^78)

    pool = load_p11_pool(dir)
    N    = size(pool.summary_min, 2)
    N == n_pairs || @warn "pool holds $N columns, expected $n_pairs"

    # --- Deterministic train/val split (the loader's own fold machinery, so the split is
    #     reproducible and independent of the generation key stream).
    folds     = Loader.all_folds(N; K = P11_TRAIN_K, master_seed = UInt64(P11_DEV_SEED))
    val_idx   = folds[P11_TRAIN_FOLD]
    train_idx = setdiff(1:N, val_idx)
    verbose && println("[1/5] split: $(length(train_idx)) train / $(length(val_idx)) val " *
                       "(K = $P11_TRAIN_K, fold $P11_TRAIN_FOLD held out)")

    # --- Frozen 128-row summary transform, fit on TRAIN columns ONLY (leak-free).
    zt       = fit_p11_summary_transform(pool.summary_min[:, train_idx])
    Ztr_128  = apply_p11_summary_transform(pool.summary_min[:, train_idx], zt)
    Zva_128  = apply_p11_summary_transform(pool.summary_min[:, val_idx],   zt)
    verbose && println("[2/5] frozen summary transform fitted train-only " *
                       "($(length(zt.mean)) continuous rows, mask rows bypassed)")

    # --- Append lambda AFTER standardization (Pitfall 2), one column at a time.
    Ztr = build_p11_inputs(Ztr_128, pool.lambda[train_idx])
    Zva = build_p11_inputs(Zva_128, pool.lambda[val_idx])
    @assert size(Ztr, 1) == 129 "training inputs must be 129 rows, got $(size(Ztr, 1))"
    @assert size(Zva, 1) == 129 "validation inputs must be 129 rows, got $(size(Zva, 1))"
    @assert length(unique(Ztr[129, :])) > 1 "the lambda row is CONSTANT across training samples " *
        "— the conditioning input carries no information (Pitfall 1 / R4)"
    verbose && println("[3/5] 129-row inputs assembled; lambda row varies over " *
                       "$(length(unique(Ztr[129, :]))) distinct values")

    # --- The PORTED bounded theta transform, fit on TRAIN theta ONLY. This is declared
    #     deviation (4) (ruling Q4). It exists so the D-02 attenuation comparison against the
    #     shipped net is not confounded by a plain z-score leaking posterior mass outside the
    #     prior box -- a known, already-diagnosed defect (F2) that would otherwise contaminate
    #     exactly the one number that comparison reports.
    θzt     = fit_p11_theta_transform(pool.theta[:, train_idx])
    θtr_std = Float32.(StatsBase.transform(θzt, pool.theta[:, train_idx]))
    θva_std = Float32.(StatsBase.transform(θzt, pool.theta[:, val_idx]))
    D_flow  = size(θtr_std, 1)
    @assert D_flow == 8 "expected 8 theta rows (D-09), got $D_flow"
    verbose && println("[4/5] ported bounded theta transform fitted train-only (D = $D_flow)")

    # --- Build and train. Capacity is HELD FIXED at the Phase-5-calibrated defaults so the
    #     comparison against the shipped net is capacity-controlled; the AdamW hyperparameters
    #     are Float64 to match NeuralEstimators' Float64 CosAnneal lr schedule (a Float32
    #     optimiser state trips `Optimisers.adjust!` when the schedule feeds a Float64 eta).
    est    = build_p11_estimator(; d_in = 129, D = 8)
    logdir = mktempdir()
    verbose && println("[5/5] training (risk log -> $logdir) ...")
    est = train(est, θtr_std, θva_std, Ztr, Zva;
                epochs = epochs,
                batchsize = batchsize,
                use_gpu = false,
                optimiser = Flux.Optimisers.AdamW(learning_rate, (0.9, 0.999), weight_decay),
                stopping_epochs = stopping_epochs,
                savepath = logdir,
                verbose = verbose)

    train_risk, val_risk = _read_risk_trajectory(logdir)
    _stall_diagnostic(train_risk, val_risk; epochs = epochs)

    imsize_counts = realized_imsize_counts(pool.imsize)
    imsize_meta   = (imsize_set = P11_IMSIZE_SET, imsize_weights = P11_IMSIZE_WEIGHTS,
                     imsize_source = :p11_consts_f5_mixture, recorded = true)
    elapsed_min   = (time() - t0) / 60

    result = (estimator = est, θzt = θzt, zt = zt, d_in = 129, D = 8)
    path   = save_p11_npe(savepath, result;
                          pool_dir = dir, n_pairs = N,
                          n_train = length(train_idx), n_val = length(val_idx),
                          elapsed_min = elapsed_min,
                          train_risk = train_risk, val_risk = val_risk,
                          imsize_counts = imsize_counts, imsize_meta = imsize_meta)

    if verbose
        println("="^78)
        println("research net persisted -> $path")
        println("  elapsed            : $(round(elapsed_min; digits = 2)) min")
        println("  epochs run         : $(max(length(val_risk) - 1, 0))")
        isempty(val_risk) ||
            println("  validation risk    : initial $(round(val_risk[1]; digits = 4)) -> " *
                    "best $(round(minimum(val_risk); digits = 4))")
        println("  realized imsizes   : $imsize_counts")
        println("="^78)
    end

    return (path = path, estimator = est, θzt = θzt, zt = zt, d_in = 129, D = 8,
            elapsed_min = elapsed_min, train_risk = train_risk, val_risk = val_risk,
            training_imsize_provenance = imsize_counts, pool_dir = dir)
end

# Bare call as the last line when run as a script (the spike house form; no PROGRAM_FILE guard,
# no arg parsing, no Pkg.activate -- the project is selected by `--project=spike`). Guarded so
# an `include` from a test or an evaluation script defines the API without launching a 20-minute
# training run.
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    train_p11_research_net()
end
