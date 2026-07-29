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

# spike/test/runtests.jl --- re-runnable hard gate for the Phase 1 smoke (D-06).
#
# Wraps spike/00_smoke.jl in a stdlib Test harness so "green smoke is the hard
# gate for all later phases" can be re-checked on demand:
#     julia --project=spike spike/test/runtests.jl
#
# Two assertions:
#   1. ENV-02 / D-03 correctness -- the NPE recovers theta_true within tolerance.
#   2. D-04 CPU-only guarantee   -- CUDA is neither a (direct or resolved)
#      dependency nor a loaded module.
#
# NOTE on the CUDA check: Flux/NNlib/Zygote/NeuralEstimators declare CUDA as a
# *weak* dependency (an optional package extension, the CUDA-as-extension design
# since Flux v0.14). Those weakdep references appear as inert text in
# Manifest.toml, so a naive `occursin("CUDA", manifest)` regex FALSE-POSITIVES.
# The correct env-level check inspects the *resolved/installed* package set via
# Pkg.dependencies() (which excludes weakdeps), plus the project's direct deps
# and the loaded-module table.

using Test
using Pkg

@testset "NeuralEstimators CPU smoke" verbose = true begin

    # Runs the full NPE train + sample pipeline; defines mu_hat, theta_true, tol.
    include(joinpath(@__DIR__, "..", "00_smoke.jl"))

    @testset "ENV-02 correctness (recovered posterior mean)" begin
        @test abs(mu_hat - theta_true) < tol
    end

    @testset "D-04 CPU-only (no CUDA dependency, none loaded)" begin
        # (a) not a direct dependency of the spike project
        @test !haskey(Pkg.project().dependencies, "CUDA")
        # (b) not a resolved/installed package node (weakdep refs excluded)
        @test !any(p -> occursin("CUDA", p.name), values(Pkg.dependencies()))
        # (c) not loaded as a module at runtime
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
        # (d) RESOLVE-RISK GATE (Phase 2): the Phase-2 imaging/plotting stack
        #     (Images/ImageFiltering/CairoMakie/HypothesisTests) must NOT have
        #     downgraded the pinned NeuralEstimators v0.2.1 -- exactly the Phase-1
        #     GLMakie-co-resolve regression (NOTES §4). Keyed by NeuralEstimators'
        #     UUID so a rename can't silently mask it.
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
        # (e) RESOLVE-RISK GATE (Phase 3): the Wave-0 dependency additions
        #     (JLD2 the sharded-cache backend, Random123 the counter-based
        #     seeding) must be present AND must NOT have co-resolved the pinned
        #     NeuralEstimators v0.2.1 downward (same GLMakie-class regression as
        #     (d), now keyed on the Phase-3 deps). Checked by name in the
        #     resolved set, then the v0.2.1 pin is re-asserted with them present.
        _deps_by_name = Dict(p.name => p for p in values(Pkg.dependencies()))
        @test haskey(_deps_by_name, "JLD2")
        @test haskey(_deps_by_name, "Random123")
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
        # (f) RESOLVE-RISK GATE (Phase 4): the Wave-0 BenchmarkTools add (the
        #     >100× wall-clock claim, NPE-03) must be present AND must NOT have
        #     co-resolved the pinned NeuralEstimators v0.2.1 downward (RESEARCH
        #     Pitfall 2, the same GLMakie-class regression as (d)/(e), now keyed
        #     on the Phase-4 dep). Turing/AdvancedVI/ForwardDiff must NEVER enter
        #     the spike env (they live only in the isolated spike/baseline/ env,
        #     D-01) -- co-resolving Turing here would cap NeuralEstimators below
        #     the pin (the Phase-1 landmine). Asserted absent to lock the boundary.
        @test haskey(_deps_by_name, "BenchmarkTools")
        @test !haskey(Pkg.project().dependencies, "Turing")
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
        # (g) RESOLVE-RISK GATE (Phase 9): the Wave-0 tabular promotions (DataFrames the
        #     per-method comparison-table backend, CSV the human-readable artifact, D-09)
        #     must be present AND must NOT have co-resolved the pinned NeuralEstimators
        #     v0.2.1 downward (RESEARCH Pitfall; the same GLMakie-class regression as
        #     (d)/(e)/(f), now keyed on the Phase-9 deps). ADDITIONALLY the optional Tapqir
        #     bridge's Python stack (PythonCall/CondaPkg, D-08) must NEVER enter the MAIN
        #     spike env — it lives only in the isolated spike/comparator/tapqir_env (T-09-03
        #     elevation boundary). Asserted absent to lock the boundary, then the v0.2.1 pin
        #     is re-asserted with DataFrames/CSV present.
        @test haskey(_deps_by_name, "DataFrames")
        @test haskey(_deps_by_name, "CSV")
        @test !haskey(Pkg.project().dependencies, "PythonCall")
        @test !haskey(Pkg.project().dependencies, "CondaPkg")
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
        # (h) RESOLVE-RISK GATE (Phase 5): the validation bundle (SBC/BF/OOD) adds NO
        #     new package — every runtime dep (NeuralEstimators/Flux/HypothesisTests/
        #     JLD2/Random123/StatsBase/CairoMakie/Distributions) is already present from
        #     Phases 1–4. The amortized BF uses NeuralEstimators' BUILT-IN NormalisingFlow
        #     (British spelling) and the OOD ROC is hand-rolled, so NO ROC/normalizing-flow
        #     package may enter the env (CLAUDE.md "What NOT to Use"; RESEARCH Pitfall 7
        #     co-resolve landmine). Turing must stay absent (it caps NeuralEstimators below
        #     the pin). Then the v0.2.1 pin is re-asserted with the Phase-5 code loadable.
        @test !haskey(Pkg.project().dependencies, "ROCAnalysis")
        @test !haskey(Pkg.project().dependencies, "MLJ")
        @test !haskey(Pkg.project().dependencies, "NormalizingFlows")
        @test !haskey(Pkg.project().dependencies, "InvertibleNetworks")
        @test !haskey(Pkg.project().dependencies, "Turing")
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
        # (i) RESOLVE-RISK GATE (Phase 11): the registration/chromatic research lane adds NO
        #     new package at all. D-01 scopes that work to a spike-local research net trained
        #     on a fresh DEV seed, with the shipped artifact and the Phase-7 GO untouched — so
        #     the pinned, reproducible dependency set must be BYTE-IDENTICAL to the pre-Phase-11
        #     one. The research net's model surface (a ported theta transform, a lambda encoder,
        #     a 129-row input, an 8-marginal flow) is built entirely from packages Phases 1–5
        #     already resolved. Asserted as an EXPLICIT NAME SET rather than a count, so an
        #     addition names itself in the failure, and then the v0.2.1 pin is re-asserted with
        #     the Phase-11 code on disk (the same GLMakie-class regression as (d)–(h)).
        @test Set(keys(Pkg.project().dependencies)) == Set([
            "BenchmarkTools", "CSV", "CairoMakie", "CoordinateTransformations", "DataFrames",
            "Distributions", "Flux", "HypothesisTests", "ImageFiltering", "ImageTransformations",
            "Images", "Interpolations", "JLD2", "NeuralEstimators", "Random123", "StatsBase"])
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
        # (j) RESOLVE-RISK GATE (Phase 13): the three-way evidence net adds NO new package
        #     either. The two-head trunk is Flux (the shipped topology copied verbatim), the AUC
        #     is the hand-rolled `roc_auc` (validation/ood.jl:319, tie-aware Mann-Whitney, ~15
        #     lines), the calibration is the shared `_bin_calibration`, and the D-16 Otsu mask
        #     comes through `Images` 0.26.2 TRANSITIVELY -- so `ImageSegmentation` and
        #     `ImageMorphology` must NOT become DIRECT dependencies just because a segmentation
        #     word appears in the design. `ROCAnalysis`/`MLJ` would duplicate the hand-rolled ROC,
        #     `NormalizingFlows`/`InvertibleNetworks` are the CLAUDE.md "What NOT to Use" flow
        #     stacks, and `Turing` caps NeuralEstimators below the pin (the Phase-1 landmine).
        #     `QuadGK` and `KernelDensity` are the RETIRED KDE Bayes-factor baseline (docs/
        #     amortized.md named limit 3): D-12 declines to validate against it, and the D-07
        #     verification is closed-form conjugate-Gaussian precisely so no quadrature is needed
        #     -- so the retired stack must stay UNREACHABLE from the spike environment rather than
        #     merely unused. Then the v0.2.1 pin is re-asserted with the Phase-13 code on disk.
        @test !haskey(Pkg.project().dependencies, "ROCAnalysis")
        @test !haskey(Pkg.project().dependencies, "MLJ")
        @test !haskey(Pkg.project().dependencies, "NormalizingFlows")
        @test !haskey(Pkg.project().dependencies, "InvertibleNetworks")
        @test !haskey(Pkg.project().dependencies, "Turing")
        @test !haskey(Pkg.project().dependencies, "ImageSegmentation")
        @test !haskey(Pkg.project().dependencies, "ImageMorphology")
        @test !haskey(Pkg.project().dependencies, "QuadGK")
        @test !haskey(Pkg.project().dependencies, "KernelDensity")
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
    end

end

# Phase-2 Wave-0 scaffold: the SIM-03 summary-contract testset runs in the same
# harness so a single `julia --project=spike spike/test/runtests.jl` is the gate.
include(joinpath(@__DIR__, "test_simulator.jl"))

# Phase-3 Wave-0 scaffold: the training-data-pipeline SC-1..SC-4 + thread-repro
# testsets run in the same harness (as skipped placeholders until later waves
# fill them) so a single `julia --project=spike spike/test/runtests.jl` stays the gate.
include(joinpath(@__DIR__, "test_data_pipeline.jl"))

# Phase-4 Wave-0 scaffold: the NPE-training / ADVI-benchmark / ablation SC1..SC5
# testsets (skipped placeholders until later waves fill them) plus the passing
# Wave-0 holdout raw-image reproducibility gate (Open Question 1) run in the same
# harness so a single `julia --project=spike spike/test/runtests.jl` stays the gate.
include(joinpath(@__DIR__, "test_npe.jl"))

# Phase-5 scaffold: the validation bundle (SBC filled by this plan 05-01; BF/OOD
# skipped placeholders until Plans 05-02/05-03) runs in the same harness so a single
# `julia --project=spike spike/test/runtests.jl` stays the gate.
include(joinpath(@__DIR__, "test_sbc.jl"))
include(joinpath(@__DIR__, "test_bf.jl"))
include(joinpath(@__DIR__, "test_ood.jl"))

# Phase-9 Wave-0 scaffold: the cross-method comparator CMP-01..08 testsets (skipped
# placeholders until later Phase-9 waves fill them) run in the same harness so a single
# `julia --project=spike spike/test/runtests.jl` stays the gate.
include(joinpath(@__DIR__, "test_comparator.jl"))

# Phase-13 three-hypothesis evidence net: the unit testsets run in the same harness so a single
# `julia --project=spike spike/test/runtests.jl` stays the gate. Dependency order -- the frozen
# pre-registration first, then the class boundary, the alpha ladder, the real-image arm, the tau
# probe, the net, the result type and the calibration surface. The real-image arm sits directly
# after the alpha ladder because its unit includes real_images.jl, which includes alpha_series.jl;
# the guarded-include idiom makes that ordering an optimisation rather than a requirement, but
# stating it documents the dependency.
#
# THE CORRECTION ARM IS LAST, AND DELIBERATELY SO. Its two pre-registered
# `max |Delta log BF| <= P13_F5_MAXABS_TOL` assertions are MEASURED MISSES (0.5498 and 0.3783
# against a ceiling of 0.25), committed as-is by plan 13-07 rather than tuned away, so its outer
# `@testset` throws at the end of the file. A thrown testset aborts the remaining includes, so any
# sibling placed after it would silently never run. Putting it last keeps every other Phase-13
# testset executing and reporting while the honest failure still surfaces and still turns the
# suite red. Resolving that shortfall is a phase-level scientific decision (accept it as a named
# limit, amend the F5 statistic, or carry it into Phase 14's abstention layer) and is NOT resolved
# by this wiring.
include(joinpath(@__DIR__, "test_p13_consts.jl"))
include(joinpath(@__DIR__, "test_p13_labels.jl"))
include(joinpath(@__DIR__, "test_p13_alpha.jl"))
include(joinpath(@__DIR__, "test_p13_real.jl"))
include(joinpath(@__DIR__, "test_p13_tau.jl"))
include(joinpath(@__DIR__, "test_p13_net.jl"))
include(joinpath(@__DIR__, "test_p13_result.jl"))
include(joinpath(@__DIR__, "test_p13_calibration.jl"))
# The D-02/D-03 Phase-11 binding gate. Its absent-artifact assertion -- that a missing Phase-11
# research NPE is a loud, instructional BLOCK naming the grid-8 prohibition -- runs in ANY state
# of the repository, while the six-item contract assertions skip explicitly until Phase 11's net
# is on disk. Placed before the correction arm because that arm throws (see above).
include(joinpath(@__DIR__, "test_p13_preconditions.jl"))
# The D-07-i stratification / D-02 frozen-zt / D-03 lambda-placement / D-11 target gate. Like the
# binding gate above, its Phase-11-dependent testsets skip EXPLICITLY (`@test_skip`, naming
# Phase 11) until the research NPE is on disk, while the source-level assertions -- no forbidden
# basis, no re-fit standardizer, no inlined recipe value, no CUDA -- run in any state of the
# repository. Placed before the correction arm because that arm throws (see above).
include(joinpath(@__DIR__, "test_p13_datagen.jl"))
# Phase-12 spatial colocalization map: EVERY Phase-12 testset is aggregated behind the single
# include below, so adding one means adding a line to `test_p12_suite.jl` rather than to this
# file. It MUST sit before the throwing Phase-13 correction arm (see above): a thrown testset
# aborts the remaining includes, so a Phase-12 aggregator placed after it would report green by
# never running. The aggregator prints one `P12-SUITE-RAN: <file>` marker per include, because
# `Test` prints testset NAMES and none of the Phase-12 names contains its filename -- so only a
# printed marker distinguishes "the include ran" from "the include was skipped".
include(joinpath(@__DIR__, "test_p12_suite.jl"))
include(joinpath(@__DIR__, "test_p13_correction.jl"))
