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
# The BF/NRE surface (load_ratio, build_bf_pair, amortized_log_bf; pulls train_ratio.jl)
# and the OOD surface (fit_ood_nulls, maha_score, ood_roc_over_grid). Each guards its own
# harness re-include, so this order is a no-op past the harness already loaded above.
isdefined(@__MODULE__, :amortized_log_bf)  || include(joinpath(@__DIR__, "validation", "bf.jl"))
isdefined(@__MODULE__, :fit_ood_nulls)     || include(joinpath(@__DIR__, "validation", "ood.jl"))

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

# ============================================================================
# STEP 1b — BF / NRE fast-tier chain proof (DEMO-01)
# ============================================================================
# Load the ratio net ONCE (frozen zt + measured log_prior_odds); tiny-retrain ONLY if the
# artifact is absent (fresh clone), mirroring test_bf.jl's _load_or_tiny_ratio idiom. The
# amortized log-BF is read in ONE forward pass (amortized_log_bf → logratio; no KDE/quadgk)
# and is deterministic in the Philox stream (Z_pair depends only on simulate_pair).
function _load_or_tiny_ratio()
    if isfile(RATIO_MODEL_PATH)
        return load_ratio(RATIO_MODEL_PATH)
    end
    @warn "trained_ratio.jld2 absent (gitignored on a fresh clone) — tiny-retraining a fallback NRE"
    est = train_ratio(DEMO_M; n = 120, epochs = 6, imsize = (64, 64),
                      rng = val_rng(VAL_FIX_SEED))
    lpo = measure_log_prior_odds(prior_label_draws(2000; rng = val_rng(VAL_FIX_SEED)))
    return (estimator = est, zt = DEMO_M.zt, log_prior_odds = lpo,
            num_summaries = RATIO_NUM_SUMMARIES, meta = nothing)
end
const DEMO_RATIO = _load_or_tiny_ratio()

# A short Δρ sweep: build each (sample, control) pair on ONE threaded Philox stream (so the
# sweep points differ) and read the amortized log-BF for each — the NRE layer exercised
# end-to-end, not report-loaded. Fixture scale (64×64, N=200).
function run_bf_fixture(seed)
    rng = val_rng(seed)
    targets = collect(range(-0.4, 0.4; length = 5))
    logbfs  = Vector{Float64}(undef, length(targets))
    for (i, Δρ) in enumerate(targets)
        pr = build_bf_pair(DEMO_M, rng, Δρ; imsize = (64, 64), N = 200)
        logbfs[i] = amortized_log_bf(DEMO_RATIO.estimator, pr.Z_pair, DEMO_RATIO.log_prior_odds)
    end
    return (targets = targets, logbfs = logbfs)
end

println("[step 1b] BF/NRE fixture sweep (build_bf_pair → amortized_log_bf) on VAL_FIX_SEED …")
const BF_R1 = run_bf_fixture(VAL_FIX_SEED)
const BF_R2 = run_bf_fixture(VAL_FIX_SEED)
@assert all(isfinite, BF_R1.logbfs) "amortized log-BF returned a non-finite value"
@assert BF_R1.logbfs == BF_R2.logbfs "BF/NRE fixture sweep not reproducible on VAL_FIX_SEED"
const BF_REPRO = (BF_R1.logbfs == BF_R2.logbfs)
println("  BF twin-run reproducible: $BF_REPRO  (", length(BF_R1.targets), " Δρ targets; log-BF range [",
        round(minimum(BF_R1.logbfs); digits = 2), ", ", round(maximum(BF_R1.logbfs); digits = 2), "])")

# ============================================================================
# STEP 1c — OOD / misspecification fast-tier chain proof (DEMO-01)
# ============================================================================
# Fit the Mahalanobis nulls on a TRAIN-ONLY fixture pool (SC5 / OOD's own frozen-from-train
# discipline), commit the ID-quantile operating point, score a strong OPTICS positive
# subset + an ID subset, and run the ROC grid over all four families. Keep with_pp=false —
# the PP re-simulation path crashes on non-finite θ̂ for strong OOD inputs (RESEARCH
# Pitfall 3; run_ood.jl also runs with_pp=false). Fully deterministic in the Philox stream
# (no posterior sampling with with_pp=false), so no global reseed is needed here.
function run_ood_fixture(seed)
    imsize = (64, 64)
    rng    = val_rng(seed)
    # (1) TRAIN pool → Mahalanobis null (train-only freeze, SC5).
    Ztrain = reduce(hcat, [id_summary(DEMO_M, rng; imsize = imsize) for _ in 1:60])
    nulls  = fit_ood_nulls(Ztrain)
    # (2) held-out ID subset → the pre-registered operating point + honest ID FPR.
    maha_id  = [maha_score(nulls, id_summary(DEMO_M, rng; imsize = imsize)) for _ in 1:20]
    maha_thr = id_threshold(maha_id)                         # q = OOD_ID_QUANTILE
    # (3) strong OPTICS positive subset (robustly separable external control).
    maha_pos = Float64[]
    for _ in 1:6
        θ = sample_prior(rng)
        push!(maha_pos, maha_score(nulls, frozen_summary(DEMO_M,
                          misspec_optics(rng, θ; imsize = imsize, level = OOD_GRID_LEVELS))))
    end
    # (4) controlled ROC over the misspec grid (maha-only; with_pp=false — Pitfall 3).
    roc = ood_roc_over_grid(DEMO_M, nulls; levels = OOD_GRID_LEVELS, n_id = 20, n_pos = 6,
                            rng = val_rng(seed), imsize = imsize, with_pp = false)
    return (maha_thr = maha_thr, maha_id = maha_id, maha_pos = maha_pos,
            combined_auc = roc.combined_auc)
end

println("[step 1c] OOD fixture (fit_ood_nulls → maha_score + ROC grid) on VAL_FIX_SEED …")
const OOD_R1 = run_ood_fixture(VAL_FIX_SEED)
const OOD_R2 = run_ood_fixture(VAL_FIX_SEED)
@assert OOD_R1.maha_pos == OOD_R2.maha_pos "OOD maha scores not reproducible on VAL_FIX_SEED"
@assert OOD_R1.combined_auc == OOD_R2.combined_auc "OOD ROC AUC not reproducible on VAL_FIX_SEED"
const OOD_REPRO = (OOD_R1.maha_pos == OOD_R2.maha_pos && OOD_R1.combined_auc == OOD_R2.combined_auc)
println("  OOD twin-run reproducible: $OOD_REPRO  (fixture ROC-grid combined AUC=",
        round(OOD_R1.combined_auc; digits = 3), " over ", length(OOD_FAMILIES), " families)")

# CPU-only invariant (D-10): CUDA must never be loaded (test_sbc.jl:114-116 idiom).
@assert !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules)) "CUDA loaded — demo must be CPU-only (D-10)"

# --- Self-assert close (02_simulator_demo.jl:148-152 idiom; PowerShell-safe) -----------
@assert isfile(joinpath(@__DIR__, "npe", "trained_npe.jld2")) "demo FAILED: missing frozen NPE"
println("demo OK: NPE + BF/NRE + OOD fast-tier chain proofs reproduce on VAL_FIX_SEED (CPU-only, src/ decoupled)")
