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

# Wave-4 ablation units under test (SC4/SC5): ablate + choose_summary (ABL-01/02).
isdefined(@__MODULE__, :ablate) || include(joinpath(@__DIR__, "..", "npe", "ablation.jl"))

# --- Pre-registered constants (declared BEFORE any reported run; RESEARCH A2/A3) --
# These are the LOCKED thresholds every later Phase-4 wave scores against. Changing
# any of them after a run would be data-snooping -- they are committed here first.
const NPE_RMSE_TOLERANCE   = 1.2       # D-09 (A2): NPE ρ_true RMSE ≤ 1.2× ADVI ρ_true RMSE
const SPEEDUP_GATE         = 100.0     # NPE-03: t_advi/t_npe > 100 at the reported thread count
                                       # VALUE UNCHANGED. Its ASSERTION is PAUSED (not retired,
                                       # not relaxed) pending a quiet-machine re-measurement --
                                       # see DEFERRED-NPE-03-WALLCLOCK in SC3 below.
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

# --- Wave-4 ablation fixture (SC4/SC5): run the :min-vs-:aug k=5 CV ablation ONCE on
# the small fixture cache (both variants already cached, zero re-simulation, D-06) and
# share the result across SC4 (RMSE) and SC5 (decision rule) so the 10 tiny trainings
# happen a single time. Epochs/posterior-draws are kept small for the quick gate; the
# REPORTED ablation (spike/npe/summary_choice.md) uses the locked training defaults.
const ABL_FIX_EPOCHS = 30
const ABL_FIX_N      = 300
const ABL_FIX_RESULT = ablate(NPE_REPRO_DIR; master_seed = NPE_MASTER_SEED, K = 5,
                              N = ABL_FIX_N, epochs = ABL_FIX_EPOCHS, batchsize = 64)

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
        #
        # ============================ DEFERRED-NPE-03-WALLCLOCK ============================
        # PAUSED, NOT RETIRED (user ruling, 2026-07-29). `SPEEDUP_GATE = 100.0` at :74 is
        # UNCHANGED and the assertion below is UNWEAKENED -- it is `@test_skip`ed, so the
        # expression and its bar both stay in the source as the pre-registered record. This
        # is a deferral of the MEASUREMENT, not a withdrawal of the CLAIM. Nothing about the
        # threshold is relaxed; re-arming is a one-word edit (`@test_skip` -> `@test`).
        #
        # WHY, IN NUMBERS (measured, `spike/npe/flow_width_ab.jl`, artifact
        # `spike/npe/flow_width_ab.jld2`, commit ca994b9 -- 9 interleaved order-rotated reps
        # at this testset's own bench_N=50 / bench_seconds=0.4, CPU-only, 1 thread):
        #   * to clear 100 the NPE needs 4.962 ms; it measures 5.178 ms. A 1.04x MISS.
        #   * run-to-run spread across all arms is 28.8 POINTS (82.42 .. 111.19). Every arm,
        #     including an 8-marginal one, CLEARS 100 in its best repetition.
        #   * recorded history on the SAME unchanged model: 92.50, 83.97, 84.44, 89.77, 68.35
        #     -- a 24-point range with the model byte-identical throughout.
        # This is a MEASUREMENT-CONDITIONS problem, NOT a performance regression: pass/fail
        # is decided by whatever else is running on the machine (the sub-88 reps were taken
        # with 4 concurrent julia processes, the >110 reps with 2).
        #
        # THE 8-MARGINAL HYPOTHESIS IS FALSIFIED AT THE ROOT. STATE.md recorded the failure
        # as caused by Phase-11 D-09 widening the flow to 8 marginals. `speedup_report` loads
        # `DEFAULT_NPE_MODEL = spike/npe/trained_npe.jld2`, whose persisted flow is `q.d = 7`,
        # `d_in = 128`, last changed in `200e971` (the Phase-5 retrain) -- Phase 11 never
        # touched it. The 8-marginal net lives only in `p11_research_npe.jld2`, which this
        # benchmark never loads. Measured anyway, the 8th marginal costs 1.67% of forward-pass
        # latency (5.208 -> 5.295 ms) = 1.78 speedup points, against that 28.8-point spread.
        #
        # WHY THE NUMBER IS THIS NOISY -- ONLY ONE SIDE OF THE RATIO IS BENCHMARKED.
        # `BenchmarkTools` (a declared dependency, CLAUDE.md "Wall-clock NPE-ms vs
        # ADVI-minutes claim (AP4)") IS used, but only on the DENOMINATOR: `_time_npe_pair`
        # takes `@belapsed` (minimum of many samples). The NUMERATOR `t_advi` is a SINGLE
        # un-replicated `time_ns()` per pair (`spike/baseline/run_advi.jl:136-137`), frozen
        # into `advi_artifact.jld2` on 2026-07-01 under that machine's then-current load and
        # never re-measured. A ratio of a benchmarked denominator to a one-shot numerator
        # cannot be tightened by repeating this testset.
        #
        # DEFERRED OBLIGATION -- THIS MUST BE RE-RUN, NOT DROPPED. At the END of v2.0
        # development, on a QUIET machine (no concurrent julia processes), with proper
        # benchmarking tooling on BOTH sides of the ratio (re-measure `t_advi` under
        # `BenchmarkTools` rather than reusing the frozen one-shot), re-arm this assertion
        # against the UNCHANGED bar of 100.0 and report the result. If it fails under those
        # conditions, THAT is a real result about the >100x claim and must be reported as one.
        # ==================================================================================
        @test_skip sr.median_speedup > SPEEDUP_GATE
        # ... asserted JOINTLY with the comparable-RMSE condition (D-08 -- never alone).
        # The RMSE half of the joint condition is NOT paused: it is an accuracy claim, is not
        # wall-clock dependent, and stays live so a genuine accuracy regression still fails.
        @test sr.rmse.npe_rho_rmse <= NPE_RMSE_TOLERANCE * sr.rmse.advi_rho_rmse

        # The pairing is real: per-pair speedups and a finite ADVI wall-clock exist.
        @test length(sr.speedups) == sr.n_pairs
        @test all(isfinite, sr.t_advi) && all(>(0), sr.t_advi)
        @test isfinite(sr.full_median_speedup)   # conservative full-N lower bound reported

        # CPU-only gate: CUDA not loaded (D-10).
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

    @testset "SC4 (ABL-01)" begin
        # ABL-01: per-parameter RMSE of the minimal (:min, 128-dim) vs augmented
        # (:aug, 142-dim) summaries under the loader's leak-free k=5 CV, ZERO
        # re-simulation (both variants already cached, D-06). The SAME build_estimator
        # architecture + master_seed train both arms — only the cached `variant`
        # differs — so the comparison isolates the SUMMARY, not the network.
        r = ABL_FIX_RESULT

        # (a) K×2×THETA_DIM RMSE table: K folds × {:min,:aug} × θ, FINITE for BOTH arms.
        # THETA_DIM (not a literal 7) since Phase 11 D-09 appended `chromatic_eps`.
        @test size(r.rmse_table) == (5, 2, THETA_DIM)
        @test all(isfinite, r.rmse_table)
        @test r.variants == (:min, :aug)                 # slice 1 = :min, slice 2 = :aug
        @test r.K == 5

        # (b) all-7 per-parameter RMSE non-negative, per fold per variant.
        @test all(>=(0), r.rmse_table)

        # (c) the ρ_true (row-1) RMSE headline is present, finite and positive for BOTH
        #     variants (the D-07 decision inputs).
        @test isfinite(r.rho_rmse_min) && r.rho_rmse_min > 0
        @test isfinite(r.rho_rmse_aug) && r.rho_rmse_aug > 0
        @test size(r.rho_per_fold) == (5, 2)
        @test all(isfinite, r.rho_per_fold)

        # (d) fold_wins_aug counts folds where :aug beats :min on ρ_true (0..K).
        @test 0 <= r.fold_wins_aug <= r.K
        # IN-04: anchor the fold-win COUNTING PREDICATE to hand-built fixtures with a
        # KNOWN number of aug-wins (an independent oracle with fixed expected outcomes),
        # rather than re-deriving on the real result the exact formula `ablate` uses
        # (which passes by construction and would not catch a `<`/`>` sign error copied
        # into both the function and the test). rho_per_fold[:,1]=min, [:,2]=aug.
        aug_wins(m) = count(f -> m[f, 2] < m[f, 1], 1:size(m, 1))
        @test aug_wins([1.0 0.5; 1.0 0.5; 1.0 0.5; 1.0 2.0; 1.0 2.0]) == 3  # aug wins folds 1-3
        @test aug_wins([1.0 2.0; 1.0 2.0; 1.0 2.0; 1.0 2.0; 1.0 2.0]) == 0  # aug never wins
        @test aug_wins([1.0 0.5; 1.0 0.5; 1.0 0.5; 1.0 0.5; 1.0 0.5]) == 5  # aug wins all folds
        # ablate's returned fold_wins_aug must agree with that fixture-verified predicate
        # applied to its own per-fold ρ_true RMSE table (catches a wrong column / flip).
        @test r.fold_wins_aug == aug_wins(r.rho_per_fold)

        # (e) the ablation ran CPU-only — CUDA never entered the process (D-10).
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

    @testset "SC5 (ABL-02)" begin
        # ABL-02: the pre-registered decision rule (ABL_REL_MARGIN + ABL_FOLD_CONSISTENCY)
        # emits :min or :aug on the ρ_true RMSE (D-07), and the written justification
        # (spike/npe/summary_choice.md) records the choice + an OOD-detectability note
        # coupling to Phase 5.
        r      = ABL_FIX_RESULT
        chosen = choose_summary(r)

        # (a) the rule returns exactly one of the two summary variants.
        @test chosen in (:min, :aug)

        # (b) IN-04: exercise choose_summary against HAND-BUILT result fixtures with
        #     KNOWN expected outcomes (an independent oracle), instead of re-deriving the
        #     rule's own formula against the fixture result (which cannot catch a sign
        #     error present in both the function and the test). The rule returns :aug IFF
        #     aug beats min by the relative margin AND wins ≥ ABL_FOLD_CONSISTENCY folds.
        #   - clearly wins: 40% better on ρ_true (≫5% margin) AND wins all 5 folds  → :aug
        @test choose_summary((rho_rmse_min = 1.0, rho_rmse_aug = 0.6,
                              fold_wins_aug = 5)) === :aug
        #   - margin too small: only 1% better though it wins every fold             → :min
        @test choose_summary((rho_rmse_min = 1.0, rho_rmse_aug = 0.99,
                              fold_wins_aug = 5)) === :min
        #   - inconsistent: 40% better yet wins only 3 (<4) folds                    → :min
        @test choose_summary((rho_rmse_min = 1.0, rho_rmse_aug = 0.6,
                              fold_wins_aug = 3)) === :min
        #   - boundary: exactly at the margin AND exactly at the consistency threshold → :aug
        @test choose_summary((rho_rmse_min = 1.0, rho_rmse_aug = 1.0 - ABL_REL_MARGIN,
                              fold_wins_aug = ABL_FOLD_CONSISTENCY)) === :aug

        # (c) the pre-registered consts are the locked values (declared before any run).
        @test ABL_REL_MARGIN == 0.05
        @test ABL_FOLD_CONSISTENCY == 4

        # (d) the ABL-02 deliverable exists and names the OOD-detectability coupling.
        choice_doc = joinpath(@__DIR__, "..", "npe", "summary_choice.md")
        @test isfile(choice_doc)                          # the written justification
        doc = read(choice_doc, String)
        @test occursin("OOD", doc)                        # the Phase-5 coupling note
    end

    @testset "SC3-scaling (NPE-03 characterization)" begin
        # D-12/D-13: the amortized thesis IS a scaling claim. The >100× is presented as a
        # CURVE over the dataset count N and the input size (imsize), with empirical log-log
        # exponents (NPE flat/amortized vs ADVI ~linear over N), plus a CPU thread-count
        # sweep run as SEPARATE julia -t N processes (Pitfall 7) and aggregated, with the
        # headline stated at the fixed BENCH_THREADS (D-13). This is CHARACTERIZATION, not a
        # pass/fail gate -- small grids / cached timings keep it fast (the full {1,2,4,8}
        # process sweep is the reported artifact, not a CI gate).
        isdefined(@__MODULE__, :scaling_over_N) ||
            include(joinpath(@__DIR__, "..", "npe", "scaling.jl"))
        isdefined(@__MODULE__, :aggregate_thread_sweep) ||
            include(joinpath(@__DIR__, "..", "npe", "run_thread_sweep.jl"))

        artifact = joinpath(@__DIR__, "..", "baseline", "advi_artifact.jld2")
        holdir   = joinpath(@__DIR__, "..", "baseline", "holdout")
        tmp      = mktempdir()

        # (a) N-axis amortization curve: NPE total (C_train + N·t_fwd) vs ADVI (N·t_advi).
        #     train_cost supplied (the baseline/holdout dir has no CV folds to time against);
        #     the reported run measures C_train from a full cache instead.
        rN = scaling_over_N(holdir; master_seed = NPE_MASTER_SEED,
                            Ns = [1, 5, 10, 20], train_cost = 5.0,
                            artifact_path = artifact, bench_N = 20, bench_seconds = 0.2,
                            save_path = joinpath(tmp, "scaling_N.jld2"))
        @test length(rN.N_grid) == 4
        @test all(isfinite, rN.npe_time_N)  && all(>(0), rN.npe_time_N)
        @test all(isfinite, rN.advi_time_N) && all(>(0), rN.advi_time_N)
        @test isfinite(rN.npe_exponent_N) && isfinite(rN.advi_exponent_N)
        # WR-04: npe_exponent_N/advi_exponent_N are log-log slopes of the CLOSED-FORM
        # analytic curves (npe = C_train + N·t_fwd; advi = N·t_advi), so their ordering
        # is true by construction and a `<` assertion cannot fail regardless of code
        # correctness. Assert instead on the MEASURED inputs that make amortization real:
        # a positive one-time training cost, a per-dataset NPE forward pass strictly
        # cheaper than a per-dataset ADVI run, and the finite positive crossover N those
        # imply (beyond which the amortized NPE undercuts ADVI).
        @test rN.train_cost > 0                          # here supplied as 5.0
        @test rN.t_fwd_per_dataset > 0
        @test rN.t_advi_per_dataset > 0
        @test rN.t_fwd_per_dataset < rN.t_advi_per_dataset   # measured amortization premise
        @test isfinite(rN.crossover_N) && rN.crossover_N > 0
        @test isfile(joinpath(tmp, "scaling_N.jld2"))    # curve persisted atomically

        # (b) input-size curve: finite exponents for BOTH methods (imsize characterization).
        rS = scaling_over_imsize(holdir; master_seed = NPE_MASTER_SEED,
                                 sizes = [(128, 128), (256, 256)],
                                 artifact_path = artifact, bench_N = 20, bench_seconds = 0.2,
                                 save_path = joinpath(tmp, "scaling_imsize.jld2"))
        @test length(rS.imsize_grid) == 2
        @test all(isfinite, rS.npe_time_sz)  && all(>(0), rS.npe_time_sz)
        @test all(isfinite, rS.advi_time_sz) && all(>(0), rS.advi_time_sz)
        @test isfinite(rS.npe_exponent_sz) && isfinite(rS.advi_exponent_sz)
        @test isfile(joinpath(tmp, "scaling_imsize.jld2"))

        # (c) thread-sweep aggregation (D-13). The full {1,2,4,8} sweep runs as SEPARATE
        #     `julia -t N` processes in the reported artifact (run_thread_sweep spawn=true --
        #     Pitfall 7: nthreads() is immutable per process). Here the aggregation MACHINERY
        #     + atomic file production is exercised on cached per-N results written by the real
        #     worker writer path (_atomic_jldsave), keeping the gate fast (plan: cached timings).
        for (n, sp) in ((1, 300.0), (2, 170.0), (4, 95.0), (8, 60.0))
            _atomic_jldsave(joinpath(tmp, "thread_$(n).jld2");
                threads = n, median_speedup = sp, full_median_speedup = sp / 20,
                t_npe_median = 0.5 / n, t_advi_median = 0.5,
                rmse_ratio_ok = true, use_gpu = false)
        end
        agg = aggregate_thread_sweep([joinpath(tmp, "thread_$(n).jld2") for n in (1, 2, 4, 8)];
                                     outpath = joinpath(tmp, "thread_sweep.jld2"),
                                     bench_threads = BENCH_THREADS)
        @test isfile(joinpath(tmp, "thread_sweep.jld2"))   # aggregation artifact produced
        @test agg.threads == [1, 2, 4, 8]                  # sorted per-N thread-scaling table
        @test agg.headline_threads == BENCH_THREADS        # headline at the fixed count (D-13)
        @test isfinite(agg.headline_speedup)
        @test length(agg.parallel_speedup) == 4 && all(isfinite, agg.parallel_speedup)
        # IN-05: the per-thread-count D-09 RMSE-validity flag is carried through
        # aggregation (and the saved artifact), so a tolerance violation stays visible.
        @test agg.rmse_ratio_ok == [true, true, true, true]
        @test JLD2.load(joinpath(tmp, "thread_sweep.jld2"), "rmse_ratio_ok") ==
              [true, true, true, true]

        # (d) CPU-only throughout (D-10).
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
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
