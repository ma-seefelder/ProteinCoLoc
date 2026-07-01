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

# Wave-2 NPE units under test (SC1): architecture + training + the inference surface.
# infer.jl pulls in train_npe.jl (train_fold/save_npe/load_npe) + architecture.jl.
isdefined(@__MODULE__, :posterior_for) || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))

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
        # NPE-01: a trained PosteriorEstimator (MLP summary net → NormalisingFlow)
        # returns a 7×N posterior in ONE forward pass for each of the ≥20 reserved
        # holdout stacks; ρ̂ = row-1 posterior mean (un-standardized) ∈ [-1,1]; and
        # Δρ = the Monte-Carlo difference of two independent single-stack passes (D-03).
        #
        # FAST GATE (RESEARCH Sampling Rate / plan): the trained model is LOADED from
        # the persisted spike/npe/trained_npe.jld2 (produced by `train_fold` at
        # NPE_MASTER_SEED in this plan) rather than retrained here. The reserved
        # holdout is read from NPE_REPRO_DIR: the holdout stream depends only on
        # (master_seed, n_holdout) -- both NPE_MASTER_SEED/20 -- so it is BYTE-IDENTICAL
        # to the holdout the persisted model was evaluated against, and was excluded
        # from that model's training pool by construction (D-10).
        SC1_N      = 1000
        model_path = joinpath(@__DIR__, "..", "npe", "trained_npe.jld2")
        @test isfile(model_path)                       # the persisted NPE-01 deliverable

        m  = load_npe(model_path)
        @test m.estimator isa NeuralEstimators.PosteriorEstimator
        @test m.d_in == 128                            # :min summary width

        ho = load_holdout(NPE_REPRO_DIR)               # reserved ≥20-stack set (D-10)
        H  = size(ho.theta, 2)
        @test H >= 20
        Zh = standardize_summary(ho.summary_min, m.zt, :min)   # frozen train-fit zt

        # (a) 7×N single-pass posterior for EACH of the ≥20 holdout stacks; and
        # (b) every ρ̂ un-standardizes to a correlation in [-1,1].
        for j in 1:H
            draws = posterior_for(m.estimator, Zh[:, j]; N = SC1_N)
            @test size(draws) == (7, SC1_N)            # one 7×N pass per stack
            r = rho_hat(m.estimator, Zh[:, j], m.θzt; N = SC1_N)
            @test isfinite(r)
            @test -1.0 <= r <= 1.0                     # row-1 (ρ_true) in range
        end

        # (c) Δρ on a (sample, control) holdout pair is a FINITE scalar and EQUALS the
        #     mean of the two independent single-stack passes (D-03). Seeding both the
        #     helper and the manual recomputation identically makes the MC draws match.
        Zs, Zc = Zh[:, 1], Zh[:, 2]
        Random.seed!(0x5C1)
        dρ = delta_rho(m.estimator, Zs, Zc, m.θzt; N = SC1_N)
        Random.seed!(0x5C1)
        ρs = rho_draws(m.estimator, Zs, m.θzt; N = SC1_N)
        ρc = rho_draws(m.estimator, Zc, m.θzt; N = SC1_N)
        @test isfinite(dρ)
        @test -2.0 <= dρ <= 2.0                        # difference of two ρ∈[-1,1]
        @test dρ ≈ mean(ρs .- ρc)                      # MC-diff of two passes (D-03)

        # (d) the whole gate ran CPU-only -- CUDA never entered the process (D-10).
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

    @testset "SC2 (NPE-02)" begin
        # NPE-02: per-parameter RMSE + 90% interval width of the trained NPE vs the
        # isolated-baseline ADVI artifact (advi_artifact.jld2), scored in ρ-space via
        # the frozen ghat and joined STRICTLY on holdout global_index (D-02/D-04).
        # The comparable-RMSE gate is the PRE-REGISTERED tolerance D-09: NPE ρ_true
        # RMSE ≤ NPE_RMSE_TOLERANCE × ADVI ρ_true RMSE (locked BEFORE this run, 04-01).
        #
        # Fixtures are the COMMITTED cross-env hand-off: the baseline holdout (whose
        # global_index -1..-20 match the artifact pairs) and advi_artifact.jld2.
        isdefined(@__MODULE__, :rmse_report) ||
            include(joinpath(@__DIR__, "..", "npe", "benchmark.jl"))

        artifact = joinpath(@__DIR__, "..", "baseline", "advi_artifact.jld2")
        holdir   = joinpath(@__DIR__, "..", "baseline", "holdout")
        @test isfile(artifact)                         # the D-02 cross-env hand-off
        @test isfile(joinpath(holdir, "holdout.jld2")) # the global_index join source

        rr = rmse_report(holdir; master_seed = NPE_MASTER_SEED,
                         artifact_path = artifact, N = 2000)

        # ≥20 stacks (10 pairs) scored on both methods (NPE-02 holdout requirement).
        @test rr.n_stacks >= 20
        @test rr.n_pairs  >= 10

        # ρ-space RMSE finite and non-negative for both methods; ADVI is the reference.
        @test isfinite(rr.npe_rho_rmse)  && rr.npe_rho_rmse  >= 0
        @test isfinite(rr.advi_rho_rmse) && rr.advi_rho_rmse >  0

        # THE GATE (D-09): NPE ρ_true RMSE within the pre-registered tolerance of ADVI.
        @test rr.npe_rho_rmse <= NPE_RMSE_TOLERANCE * rr.advi_rho_rmse

        # All 7 θ per-parameter RMSE produced (ρ_true is the headline; D-07 reports all).
        @test length(rr.npe_all7_rmse) == 7
        @test all(isfinite, rr.npe_all7_rmse)
        @test rr.npe_all7_rmse[1] == rr.npe_rho_rmse   # row 1 IS the ρ_true headline

        # 90% interval widths reported for BOTH methods (finite, positive).
        @test isfinite(rr.npe_interval_width)  && rr.npe_interval_width  > 0
        @test isfinite(rr.advi_interval_width) && rr.advi_interval_width > 0

        # Δρ RMSE scored on identical pairs for both methods (D-04).
        @test isfinite(rr.delta_rho_rmse_npe)
        @test isfinite(rr.delta_rho_rmse_advi)

        # CPU-only: CUDA never entered the process (D-10).
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

    @testset "SC3 (NPE-03)" begin
        # NPE-03: the amortized NPE is >SPEEDUP_GATE (100×) faster than per-dataset
        # ADVI at COMPARABLE RMSE (D-08/D-09). Speedup is asserted ONLY jointly with
        # the SC2 comparable-RMSE condition -- never alone (D-08). The headline is
        # stated at the pre-registered BENCH_THREADS thread count (D-13), CPU-only.
        #
        # NPE clock = summary extraction + forward pass from the RE-SIMULATED raw stack
        # (training excluded); ADVI clock = the per-pair vi() wall_clock in the artifact.
        # Both are "to-posterior": ADVI's clock excludes rand(q,100k), so the NPE clock
        # times the forward pass at a lightweight bench_N (bulk draws are the excluded
        # analog); accuracy is scored at the full N in the paired rmse result.
        isdefined(@__MODULE__, :speedup_report) ||
            include(joinpath(@__DIR__, "..", "npe", "benchmark.jl"))

        artifact = joinpath(@__DIR__, "..", "baseline", "advi_artifact.jld2")
        holdir   = joinpath(@__DIR__, "..", "baseline", "holdout")

        # The headline is only valid at the pre-registered thread count (D-13).
        @test Threads.nthreads() == BENCH_THREADS

        sr = speedup_report(holdir; master_seed = NPE_MASTER_SEED,
                            artifact_path = artifact, N = 2000, bench_N = 50,
                            bench_seconds = 0.4)

        # Recorded at the fixed, reported thread count, CPU-only (D-13/D-10).
        @test sr.threads == BENCH_THREADS
        @test sr.use_gpu == false

        # THE GATE (NPE-03): median speedup > SPEEDUP_GATE ...
        @test sr.median_speedup > SPEEDUP_GATE
        # ... asserted JOINTLY with the comparable-RMSE condition (D-08 -- never alone).
        @test sr.rmse.npe_rho_rmse <= NPE_RMSE_TOLERANCE * sr.rmse.advi_rho_rmse

        # The pairing is real: per-pair speedups and a finite ADVI wall-clock exist.
        @test length(sr.speedups) == sr.n_pairs
        @test all(isfinite, sr.t_advi) && all(>(0), sr.t_advi)
        @test isfinite(sr.full_median_speedup)   # conservative full-N lower bound reported

        # CPU-only gate: CUDA not loaded (D-10).
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
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
