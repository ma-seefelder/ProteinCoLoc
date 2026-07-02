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

# spike/test/test_comparator.jl --- Phase-9 CMP-01..08 gates.
#
# Wave-0 MISSING scaffold (this plan, 09-01): the success-criterion testsets for the
# cross-method comparator harness exist as named, skipped placeholders so the single
# `julia --project=spike spike/test/runtests.jl` gate already enumerates every CMP
# requirement the later Phase-9 waves must turn green. Each child testset holds one
# `@test_skip true` tagged with the CMP-* requirement + the wave/plan that replaces it
# with the real gate. Mirrors test_data_pipeline.jl exactly (license header ->
# using Test -> guarded include of the pre-declared consts -> one outer
# `@testset ... verbose = true`).
#
# The pre-declared thresholds this harness scores against are committed in THIS Wave-0
# plan (spike/comparator/config.jl), BEFORE any comparison table exists — the D-14
# anti-data-snooping guarantee. This scaffold includes them so the later CMP-04
# traffic-light gate reads the SAME locked bands the artifact header quotes.

using Test

# Pre-declared comparator consts (D-14): COSTES_*, DIVERGENCE_WARN/FAIL, MASTER_SEED,
# TAPQIR_*. Guarded so the file loads both standalone and after runtests.jl already
# pulled the consts into scope; config.jl is pure `const` declarations, no logic.
isdefined(@__MODULE__, :COSTES_N_SCRAMBLE) ||
    include(joinpath(@__DIR__, "..", "comparator", "config.jl"))

@testset "Phase 9 — Cross-Method Comparator Harness" verbose = true begin

    @testset "CMP-01 finiteness on shared inputs" begin
        # Later wave: Costes-p, Manders M1/M2, Pearson, Spearman all emit FINITE,
        # NaN-safe scalars on a non-degenerate shared MultiChannelImage fixture (SC1;
        # D-04/D-06). Replaces this skip with the real per-estimator finiteness gate.
        @test_skip true
    end

    @testset "CMP-02 Costes/Manders formula-equality" begin
        # Later wave: the comparator's M1/M2 use the SAME Otsu-threshold source
        # (mci.otsu_threshold) as encode_aug so the ablation-augmented summary and the
        # comparator never silently diverge; Costes-p reproduces the seeded scramble
        # null (D-04/D-05; specifics §Manders threshold-source parity).
        @test_skip true
    end

    @testset "CMP-03 shared byte-identical inputs" begin
        # Later wave: every estimator AND the NPE consume byte-identical MultiChannelImage
        # inputs built via spike/contract.jl:build_mci over the seeded θ grid (SC1; D-01/D-02).
        @test_skip true
    end

    @testset "CMP-04 divergence traffic-light bands (D-14)" begin
        # Later wave: the divergence indicator flags rows where a classical verdict
        # diverges from the simulator ground-truth regime, presented as the
        # BayesInteractomics green/amber/red traffic light against the PRE-DECLARED
        # DIVERGENCE_WARN/FAIL bands committed in config.jl (D-10/D-14). The locked
        # bands the real gate scores against already exist and are consistent here.
        @test_skip true
    end

    @testset "CMP-05/08 table determinism (seeded)" begin
        # Later wave: two seeded runs of run_comparator.jl at MASTER_SEED produce a
        # bit-identical table across runs and thread counts (SC; D-12). Costes null,
        # table assembly, and artifact all reproduce from the Philox master seed.
        @test_skip true
    end

    @testset "CMP-06 Tapqir anchor-or-skip" begin
        # Later wave (09-04): the Tapqir bridge either reproduces the PUBLISHED anchor
        # within TAPQIR_TOL of TAPQIR_PUBLISHED_VALUE OR is cleanly skipped-with-flag
        # when the PythonCall/CondaPkg env is absent — never a hard failure of the
        # classical battery (SC2; D-07/D-08 graceful degradation).
        @test_skip true
    end

    @testset "CMP-07 audit/diagnostics report" begin
        # Later wave: reuse the BayesInteractomics comparator/audit + diagnostics
        # pattern (_bin_calibration/CalibrationResult, model_diagnostics/_ks_test_uniform
        # reporting style) for the per-method comparison table + audit summary (SC3; D-11).
        @test_skip true
    end

end
