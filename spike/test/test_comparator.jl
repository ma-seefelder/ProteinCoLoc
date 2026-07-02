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

# spike/test/test_comparator.jl --- Phase-9 CMP-01..08 gates (FILLED, 09-06 / D-13).
#
# The D-13 acceptance gate for the cross-method comparator harness, wired into the single
# `julia --project=spike spike/test/runtests.jl` entry point. This plan (09-06) REPLACES
# the Wave-0 skipped placeholders with the real assertions that prove SC1 (per-method
# finite table on shared inputs), SC2 (Tapqir anchor-or-skip), and SC3 (BayesInteractomics-
# pattern, seeded/reproducible).
#
# ORACLE DISCIPLINE (T-09-18): the CMP-02 Manders and CMP-04 traffic-light / divergence
# gates assert against INDEPENDENT hand-built oracles (the encode_aug M1/M2 slice; known
# divergence values with known expected colors; known ground-truth-vs-classical
# mismatches) — NEVER by re-deriving the function's own formula. A `<`/`>` sign flip or a
# no-op divergence copied into both the code and the test would slip past a self-referential
# check; the fixture oracles here catch it.
#
# The pre-declared thresholds this harness scores against were committed in Wave-0
# (spike/comparator/config.jl), BEFORE any comparison table existed — the D-14 anti-data-
# snooping guarantee. The traffic-light gate reads the SAME locked bands the artifact
# header quotes.

using Test
using Random
using JLD2

# Bring the full comparator harness into scope (guarded so this file loads standalone AND
# after runtests.jl already pulled a sibling in). run_comparator.jl transitively loads
# contract.jl (build_mci/patch_summary + the read-only src/ coupling), config.jl (the
# pre-declared D-14 consts: COSTES_*, DIVERGENCE_WARN/FAIL, MASTER_SEED, TAPQIR_*),
# classical.jl, inputs.jl, table.jl, audit.jl, tapqir_bridge.jl.
isdefined(@__MODULE__, :run_comparator) ||
    include(joinpath(@__DIR__, "..", "comparator", "run_comparator.jl"))
# encode_aug (the CMP-02 Manders-equality oracle; M1/M2 at indices 129,130 of the AUG_DIM
# vector) lives in spike/data/encode.jl and is NOT pulled in by run_comparator; guarded
# include so the equality gate can compare manders(mci) to that independent slice.
isdefined(@__MODULE__, :encode_aug) ||
    include(joinpath(@__DIR__, "..", "data", "encode.jl"))

# --- Shared fixtures (built ONCE, small + seeded, so the gate stays fast) ----------------
# A non-degenerate, positively-correlated 2-channel fixture: finite everywhere, non-zero
# per-channel sums (so Manders never divides by zero), no exact-zero pixels (so no patch is
# dropped by the frozen ≤15-survivor trap). Built via the frozen build_mci (contract.jl).
const CMP_FIX_MCI = let
    rng  = MersenneTwister(0xC0FFEE)
    base = rand(rng, 64, 64)
    ch1  = base .+ 0.10 .* rand(rng, 64, 64)
    ch2  = 0.70 .* base .+ 0.30 .* rand(rng, 64, 64)
    build_mci(Matrix{Float64}[ch1, ch2]; name = "cmp-fixture")
end

# One tiny seeded shared-input set (1 per regime) + its assembled table, reused by the
# finiteness/audit gates. Full Costes N runs once here; the per-estimator gates below use a
# reduced Costes n for speed.
const CMP_FIX_INPUTS = build_shared_inputs(; master_seed = MASTER_SEED,
                                           n_per_regime = 1, imsize = (64, 64))
const CMP_FIX_TABLE  = build_table(CMP_FIX_INPUTS; master_seed = MASTER_SEED)

# Reduced Costes scramble count for the direct-estimator gates (speed only; the range/repro
# properties hold for any n). The determinism gate keeps the FULL default N (D-05) but on a
# tiny input set so it stays fast.
const CMP_TEST_N = 8

@testset "Phase 9 — Cross-Method Comparator Harness" verbose = true begin

    @testset "CMP-01 finiteness on shared inputs" begin
        # SC1 (D-04/D-06): Costes-p, Manders M1/M2, Pearson, Spearman all emit FINITE,
        # NaN-safe scalars on a non-degenerate shared MultiChannelImage fixture.
        cp     = costes_p(MASTER_SEED, 1, CMP_FIX_MCI; n = CMP_TEST_N)
        m1, m2 = manders(CMP_FIX_MCI)
        @test isfinite(cp)
        @test isfinite(m1)
        @test isfinite(m2)
        @test isfinite(pearson_whole(CMP_FIX_MCI))
        @test isfinite(spearman_whole(CMP_FIX_MCI))

        # And on the seeded regime-labelled inputs the whole classical battery is finite.
        @test all(isfinite, CMP_FIX_TABLE.Costes_p)
        @test all(isfinite, CMP_FIX_TABLE.M1)
        @test all(isfinite, CMP_FIX_TABLE.M2)
        @test all(isfinite, CMP_FIX_TABLE.Pearson)
        @test all(isfinite, CMP_FIX_TABLE.Spearman)
    end

    @testset "CMP-02 Costes/Manders formula-equality" begin
        # (Manders ≡ encode_aug) INDEPENDENT ORACLE, D-04: the comparator's M1/M2 must
        # equal the M1/M2 slice (indices 129,130) of encode_aug's AUG_DIM vector for the
        # SAME mci — both use mci.otsu_threshold, so the ablation-augmented summary and the
        # comparator column can NEVER silently diverge, with encode.jl left byte-unchanged.
        # The formulas are identical bit-for-bit, so `==` (not `≈`) is the correct assertion.
        M   = patch_summary(CMP_FIX_MCI)
        aug = encode_aug(CMP_FIX_MCI, M)
        @test length(aug) == AUG_DIM                    # 142; indices 129,130 exist
        m1, m2 = manders(CMP_FIX_MCI)
        @test m1 == aug[129]                            # M1 slice (encode.jl:128-129)
        @test m2 == aug[130]                            # M2 slice

        # (Costes bit-reproducibility + range, D-05): the seeded scramble null reproduces
        # exactly, and the Davison–Hinkley +1/+1 correction bounds p to [1/(n+1), 1].
        p1 = costes_p(MASTER_SEED, 1, CMP_FIX_MCI; n = CMP_TEST_N)
        p2 = costes_p(MASTER_SEED, 1, CMP_FIX_MCI; n = CMP_TEST_N)
        @test p1 == p2                                  # bit-reproducible (same seed/idx)
        @test 1 / (CMP_TEST_N + 1) ≤ p1 ≤ 1.0           # exact +1/+1 bounds for this n
    end

    @testset "CMP-03 shared byte-identical inputs" begin
        # SC1 (D-01/D-02/D-12): every estimator (and, later, the NPE) consumes byte-identical
        # MultiChannelImage inputs built via build_mci over the seeded θ grid, and two
        # same-seed builds are BIT-EXACT (`==`, not `≈`) in both images and labels.
        b1 = build_shared_inputs(; master_seed = MASTER_SEED, n_per_regime = 1, imsize = (64, 64))
        b2 = build_shared_inputs(; master_seed = MASTER_SEED, n_per_regime = 1, imsize = (64, 64))
        @test length(b1.images) == 3                    # 1 per (:coloc,:random,:exclusion)
        for i in 1:length(b1.images)
            @test b1.images[i].data == b2.images[i].data # bit-exact raw channel matrices
        end
        @test b1.regime   == b2.regime                  # labels reproduce exactly
        @test b1.rho_true == b2.rho_true                # ground-truth ρ reproduces exactly
        @test all(r -> r in (:coloc, :random, :exclusion), b1.regime)   # regime ⊆ the 3 bands
    end

    @testset "CMP-04 divergence traffic-light bands (D-14)" begin
        # (i) TRAFFIC-LIGHT vs hand-built divergence values with KNOWN expected colors at the
        #     pre-declared band boundaries (an independent known-answer oracle, NOT the
        #     function's own `<`/`<` formula): green below WARN, amber AT WARN, red AT FAIL.
        @test traffic_light(0.0)                    === :green
        @test traffic_light(DIVERGENCE_WARN - 0.01) === :green
        @test traffic_light(DIVERGENCE_WARN)        === :amber   # amber begins exactly at WARN
        @test traffic_light(DIVERGENCE_FAIL - 0.01) === :amber
        @test traffic_light(DIVERGENCE_FAIL)        === :red     # red begins exactly at FAIL
        @test traffic_light(1.0)                    === :red
        @test traffic_light(missing)                === :gray    # no trusted ground truth

        # (ii) DIVERGENCE SEMANTIC correctness — the phase's core positioning payload —
        #      against known ground-truth-vs-classical MISMATCHES (independent fixtures with
        #      known expected verdicts). Catches a no-op or sign-inverted `divergence`.
        #   - a :random regime with a spuriously high SIGNIFICANT Pearson MUST flag red
        #     (classic falsely calls coloc):
        @test traffic_light(divergence(:random, classical_verdict(0.8, 0.001))) === :red
        #   - a true :coloc with a high SIGNIFICANT Pearson MUST stay green (classic right):
        @test traffic_light(divergence(:coloc, classical_verdict(0.9, 0.001))) === :green
        #   - a MISSED real :coloc (non-significant Pearson) MUST flag red (classic misses it):
        @test traffic_light(divergence(:coloc, classical_verdict(0.05, 0.9))) === :red
        #   - an :unknown regime fabricates NO divergence:
        @test divergence(:unknown, 0.5) === missing
    end

    @testset "CMP-05/08 table determinism (seeded)" begin
        # SC / D-12: two seeded run_comparator runs into SEPARATE tempdirs produce a
        # BYTE-IDENTICAL table payload (`==`/`isequal`) and resolve to the IDENTICAL
        # content-hash artifact dir, independent of thread count. run_tapqir=false keeps the
        # determinism gate subprocess-free and fully seeded (the anchor is exercised in
        # CMP-06); the FULL Costes N runs on a tiny 3-input set so the gate stays fast.
        t1 = mktempdir(); t2 = mktempdir()
        r1 = run_comparator(n_per_regime = 1, imsize = (64, 64), outdir = t1, run_tapqir = false)
        r2 = run_comparator(n_per_regime = 1, imsize = (64, 64), outdir = t2, run_tapqir = false)

        @test basename(r1.artifact_dir) == basename(r2.artifact_dir)   # identical content hash

        # Byte-identical JLD2 payload: reload each table.jld2 and compare bit-exactly
        # (isequal handles the `missing` NPE columns; `==` on the pure-Float Costes column).
        tab1 = JLD2.load(joinpath(r1.artifact_dir, "table.jld2"), "table")
        tab2 = JLD2.load(joinpath(r2.artifact_dir, "table.jld2"), "table")
        @test isequal(tab1, tab2)                                       # bit-exact payload
        @test r1.table.Costes_p == r2.table.Costes_p                    # seeded null reproduces
        @test isfile(joinpath(r1.artifact_dir, "table.csv"))           # human artifact present
    end

    @testset "CMP-06 Tapqir anchor-or-skip" begin
        # SC2 (D-07/D-08): the Tapqir bridge either reproduces the PUBLISHED anchor within
        # tolerance OR is cleanly skipped-with-flag when the PythonCall/CondaPkg/Tapqir stack
        # is unavailable — NEVER a hard failure of the classical battery. On this machine the
        # documented clean-skip path is taken (TAPQIR_PUBLISHED_VALUE=NaN, no materialized
        # Python), so :skipped is the expected-and-accepted outcome.
        anchor = tapqir_anchor()
        @test anchor.status in (:passed, :skipped)                     # anchor-or-skip
        @test hasproperty(anchor, :value) && hasproperty(anchor, :reason)
    end

    @testset "CMP-07 audit/diagnostics report" begin
        # SC3 (D-11): the BayesInteractomics-style audit renders a non-empty Markdown report
        # with the traffic-light band-count rows (the "classic is wrong here" tally) and the
        # Costes-p KS-uniformity section.
        report = comparator_audit(CMP_FIX_TABLE)
        @test report isa String
        @test !isempty(report)
        @test occursin("Traffic-Light Band Counts", report)            # the band-count section
        @test occursin("| Band | Rows | Fraction |", report)          # band-count table rows
        @test occursin("**total**", report)                            # totals row
        @test occursin("KS", report)                                   # Costes-p KS-uniformity
    end

end
