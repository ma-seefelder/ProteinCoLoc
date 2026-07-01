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

# spike/test/test_npe.jl --- Phase-4 NPE-01/02/03 + ABL-01/02 gates.
#
# Wave-0 MISSING scaffold (this plan, 04-01): the success-criterion testsets for
# NPE training, the ADVI benchmark, and the summary ablation exist as named,
# skipped placeholders so the single `julia --project=spike spike/test/runtests.jl`
# gate already enumerates every SC the later Phase-4 waves must turn green. Each SC
# child testset holds one `@test_skip true` tagged with the requirement + wave that
# replaces it with the real gate (RESEARCH §Validation Architecture; PATTERNS
# §spike/test/test_npe.jl). Mirrors test_data_pipeline.jl exactly (license header ->
# using Test/Random/Statistics -> guarded includes -> `const`-fixtures-before-gate ->
# one outer `@testset ... verbose = true`).
#
# The six PRE-REGISTERED CONSTANTS below are committed in THIS Wave-0 plan, BEFORE
# any reported NPE/ADVI/ablation run, so the pass/fail thresholds cannot be tuned to
# the result (the "tune until calibrated" data-snooping hazard, STATE Blockers;
# RESEARCH A2/A3). They are the contract the later waves score against.
#
# ONE real gate ships now: the Wave-0 holdout raw-image reproducibility round-trip
# (Open Question 1). It proves the reserved ADVI holdout is byte-reproducible from
# its stored (θ, global_index, master_seed) so the baseline (04-04) and benchmark
# (04-05) can both re-simulate identical raw stacks. If the round-trip is not
# bit-identical this testset MUST FAIL (do NOT loosen `==` to `≈`) -- that failure
# is the signal to persist raw holdout images instead of trusting the stream.

using Test
using Random
using Statistics
using BenchmarkTools      # the >100× wall-clock harness (NPE-03); used from later waves

# Wave-2/3 generation core: brings generate_cache (the fixture builder) AND the
# UNCHANGED Phase-2 chain (contract.jl -> build_mci/patch_summary, prior/forward,
# seeding, encode_d01) into scope via guarded includes. Safe both standalone and
# after test_data_pipeline.jl already pulled generate.jl in inside runtests.jl.
isdefined(@__MODULE__, :generate_cache) || include(joinpath(@__DIR__, "..", "data", "generate.jl"))
using JLD2

# Wave-4 loader: load_holdout (the sole holdout.jld2 reader). Guarded for idempotency.
isdefined(@__MODULE__, :load_fold) || include(joinpath(@__DIR__, "..", "data", "loader.jl"))

# Phase-4 unit under test for the Wave-0 gate: the keyed holdout re-simulation helper.
isdefined(@__MODULE__, :resimulate_holdout) || include(joinpath(@__DIR__, "..", "npe", "resimulate.jl"))

# --- Pre-registered constants (declared BEFORE any reported run; RESEARCH A2/A3) --
# These are the LOCKED thresholds every later Phase-4 wave scores against. Changing
# any of them after a run would be data-snooping -- they are committed here first.
const NPE_RMSE_TOLERANCE   = 1.2       # D-09 (A2): NPE ρ_true RMSE ≤ 1.2× ADVI ρ_true RMSE
const SPEEDUP_GATE         = 100.0     # NPE-03: t_advi/t_npe > 100 at the reported thread count
const ABL_REL_MARGIN       = 0.05      # D-07 (A3): aug wins iff RMSE_aug ≤ (1−δ)·RMSE_min on ρ_true …
const ABL_FOLD_CONSISTENCY = 4         # … AND aug improves in ≥4 of 5 folds (else keep :min, parsimony)
const BENCH_THREADS        = 1         # D-13: headline thread count (sweep {1,2,4,8} reported separately)
const NPE_MASTER_SEED      = 0xC0FFEE  # Random123 reproducibility for the Phase-4 runs

# --- Wave-0 fixture: one tiny content-hash cache with a ≥20-stack reserved holdout.
# Kept small so the quick gate stays fast; the holdout is what the repro gate probes.
const NPE_REPRO_N       = 24           # small main pool (holdout is drawn separately)
const NPE_REPRO_HOLDOUT = 20           # generate.jl clamps n_holdout to ≥20 (D-10)
const NPE_REPRO_DIR     = generate_cache(mktempdir(); N = NPE_REPRO_N,
                                         master_seed = NPE_MASTER_SEED,
                                         n_holdout = NPE_REPRO_HOLDOUT)

@testset "Phase 4 — NPE Training + ADVI Benchmark + Ablation" verbose = true begin

    @testset "SC1 (NPE-01)" begin
        # Later wave (NPE training): a trained PosteriorEstimator (MLP summary net →
        # NormalisingFlow) returns a 7×N posterior in one forward pass for the ≥20
        # holdout stacks, and Δρ = MC-difference of two single-stack passes (D-03).
        @test_skip true
    end

    @testset "SC2 (NPE-02)" begin
        # Later wave (ADVI benchmark): per-parameter RMSE + interval width of the NPE
        # vs the isolated-baseline ADVI artifact (advi_artifact.jld2), in ρ-space via
        # the frozen ghat, CV-reported over the reserved holdout (D-02/D-04).
        @test_skip true
    end

    @testset "SC3 (NPE-03)" begin
        # Later wave (speedup benchmark): t_advi / t_npe > SPEEDUP_GATE at comparable
        # RMSE (≤ NPE_RMSE_TOLERANCE) at BENCH_THREADS, plus the scaling curves over N
        # and imsize (D-08/D-09/D-12/D-13), timed with BenchmarkTools.@belapsed.
        @test_skip true
    end

    @testset "SC4 (ABL-01)" begin
        # Later wave (ablation): per-parameter RMSE of the minimal (:min, 128-dim) vs
        # augmented (:aug, 142-dim) summaries under leak-free k=5 CV, zero re-simulation
        # (both variants already cached, D-06).
        @test_skip true
    end

    @testset "SC5 (ABL-02)" begin
        # Later wave (ablation decision): the pre-registered rule (ABL_REL_MARGIN +
        # ABL_FOLD_CONSISTENCY) emits :min or :aug on ρ_true, all 7 θ reported, with the
        # OOD-detectability note coupling to Phase 5 (D-07).
        @test_skip true
    end

    @testset "Wave-0 holdout raw-image reproducibility" begin
        # Open Question 1 (RESEARCH): the reserved ADVI holdout stores SUMMARIES, not
        # raw images. Prove each reserved stack is byte-reproducible from its stored
        # (θ, global_index, master_seed) via the keyed holdout_rng stream (D-10) so the
        # baseline (04-04) and the benchmark (04-05) can re-simulate identical raw
        # stacks. BIT-EXACT `==` is mandatory -- a `≈` here would hide an unreproducible
        # stream (the failure signal to persist raw images instead).
        dir = NPE_REPRO_DIR
        ho  = load_holdout(dir)
        H   = size(ho.theta, 2)
        @test H >= 20                                   # reserved set present (D-10)

        # Re-simulate ≥2 holdout entries and assert the summary round-trips exactly.
        for j in 1:3
            rs   = resimulate_holdout(dir, j; master_seed = NPE_MASTER_SEED)
            s_re = encode_d01(patch_summary(rs.mci_sample))       # 128-dim D-01 summary
            @test s_re == ho.summary_min[:, j]                    # BIT-EXACT round-trip
            @test rs.theta == ho.theta[:, j]                      # θ round-trips exactly
            @test rs.global_index == ho.global_index[j]           # negative index echoed (D-10)
        end

        # A bad holdout index must RAISE, never silently return a wrong stack (ASVS V5).
        @test_throws ArgumentError resimulate_holdout(dir, H + 1; master_seed = NPE_MASTER_SEED)
    end

end
