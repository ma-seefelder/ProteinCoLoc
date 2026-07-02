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

# spike/test/test_bf.jl --- Phase-5 amortized-Bayes-factor gate (BF-01/02, D-07/D-08).
#
# The FAST per-task gate. It loads the PERSISTED trained ratio net (trained_ratio.jld2)
# — it does NOT retrain in the gate (retrains a tiny net only if the artifact is absent,
# at a small imsize/epochs, mirroring test_npe.jl's fast-gate-loads-fixture idiom). All
# fixtures draw from the FIXTURE stream val_rng(VAL_FIX_SEED) so the quick gate NEVER
# consumes the reserved reported VAL_MASTER_SEED stream (D-02).
#
# The FULL reported reproduction (BF_SWEEP_N=25 at the LOCKED BF_CORR_MIN/BF_LOGBF_TOL) is
# owned by Wave-3 (05-04 run_bf.jl). Here the reproduction assertion runs only when the
# persisted sweep + baseline artifacts exist and uses a FIXTURE-RELAXED correlation floor
# (BF_FIX_CORR_MIN) — the pre-registered constants themselves are asserted-as-locked, never
# tuned. See 05-02-SUMMARY.md "Deviations": at fixture scale the magnitude criterion (D-08b)
# is NOT met (the under-trained NRE is order-correct but magnitude-inflated) — reported honestly.
#
# Mirrors test_sbc.jl: license header → using Test → guarded include of the unit under test
# → fixtures computed ONCE → one outer @testset. No figure call (figures are script artifacts).

using Test
using Random

# Unit under test: the amortized BF surface (pulls bf.jl → train_ratio.jl → harness.jl).
isdefined(@__MODULE__, :amortized_log_bf) || include(joinpath(@__DIR__, "..", "validation", "bf.jl"))

# FIXTURE-RELAXED correlation floor for the small per-task sweep. The REPORTED gate (05-04)
# asserts corr ≥ BF_CORR_MIN (0.95) AND max_abs_err ≤ BF_LOGBF_TOL (0.5) at the full sweep.
const BF_FIX_CORR_MIN = 0.75

# --- Frozen model + ratio net loaded ONCE (persisted; tiny-retrain only if absent) ------
const BF_FIX_M = load_frozen_model()

function _load_or_tiny_ratio()
    if isfile(RATIO_MODEL_PATH)
        return load_ratio(RATIO_MODEL_PATH)
    end
    est = train_ratio(BF_FIX_M; n = 120, epochs = 6, imsize = (64, 64),
                      rng = val_rng(VAL_FIX_SEED))
    lpo = measure_log_prior_odds(prior_label_draws(2000; rng = val_rng(VAL_FIX_SEED)))
    return (estimator = est, zt = BF_FIX_M.zt, log_prior_odds = lpo,
            num_summaries = RATIO_NUM_SUMMARIES, meta = nothing)
end
const BF_FIX_RATIO = _load_or_tiny_ratio()

# One fixture (sample, control) pair for the SC3 one-pass-scalar check (small imsize).
const BF_FIX_PAIR = build_bf_pair(BF_FIX_M, val_rng(VAL_FIX_SEED), 0.3; imsize = (64, 64), N = 200)

# Amortized-side source, for the BF-02 "no kde/quadgk on the amortized side" grep gate.
const BF_SRC = read(joinpath(@__DIR__, "..", "validation", "bf.jl"), String)

@testset "Phase 5 — Amortized Bayes Factor (BF-01/02, D-07/D-08)" verbose = true begin

    @testset "SC3 (BF-01) amortized log-BF is a finite scalar in ONE pass" begin
        lbf = amortized_log_bf(BF_FIX_RATIO.estimator, BF_FIX_PAIR.Z_pair,
                               BF_FIX_RATIO.log_prior_odds)
        @test lbf isa Real
        @test isfinite(lbf)
        # measured log_prior_odds is finite and ≈ 0 by the symmetric i.i.d. model-index
        # prior — but MEASURED, never assumed 0 (Pitfall 5 / D-07).
        @test isfinite(BF_FIX_RATIO.log_prior_odds)
        @test abs(BF_FIX_RATIO.log_prior_odds) < 0.5
    end

    @testset "BF-02 amortized path uses NO kde/quadgk (one forward pass)" begin
        # The amortized side is logratio-only; kde/quadgk live ONLY in the isolated
        # baseline env (spike/baseline/run_bf_baseline.jl).
        @test !occursin("kde(", BF_SRC)
        @test !occursin("quadgk(", BF_SRC)
        @test occursin("logratio(", BF_SRC)
    end

    @testset "D-08 pre-registration (both criteria locked; anti-snooping)" begin
        # The pass/fail thresholds the REPORTED reproduction (05-04) scores against are the
        # committed pre-registered values — asserted-as-locked here, NEVER tuned to pass.
        @test BF_CORR_MIN  == 0.95        # D-08a correlation floor
        @test BF_LOGBF_TOL == 0.5         # D-08b bounded |Δ log-BF|
        @test BF_SWEEP_N   == 25
        @test BF_SWEEP_LO  == -0.6
        @test BF_SWEEP_HI  == 0.8
        # fixtures never consume the reserved reported stream (D-02).
        @test VAL_FIX_SEED != VAL_MASTER_SEED
    end

    @testset "BF-02 reproduction vs KDE baseline (fixture sweep)" begin
        if isfile(BF_BASELINE_ARTIFACT_PATH) && isfile(BF_SWEEP_DRAWS_PATH)
            # Load the persisted amortized sweep + KDE-baseline artifact (computed on
            # IDENTICAL Δρ inputs by the isolated baseline env) and score the D-08 gate.
            rep = bf_reproduction(BF_FIX_M, BF_FIX_RATIO.estimator,
                                  BF_FIX_RATIO.log_prior_odds; regenerate = false)
            @test length(rep.logbf_amortized) == length(rep.logbf_kde)
            @test isfinite(rep.corr) && isfinite(rep.max_abs_err)
            # D-08a ORDERING reproduces even at fixture scale (relaxed floor; the REPORTED
            # gate asserts corr ≥ BF_CORR_MIN AND max_abs_err ≤ BF_LOGBF_TOL at 05-04).
            @test rep.corr >= BF_FIX_CORR_MIN
        else
            @test_skip "BF reproduction artifacts absent (bf_sweep_draws.jld2 / " *
                       "bf_baseline_artifact.jld2) — run generate_bf_sweep (bf.jl) then " *
                       "`julia --project=spike/baseline spike/baseline/run_bf_baseline.jl` (05-04)"
        end
    end

    @testset "BF ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
