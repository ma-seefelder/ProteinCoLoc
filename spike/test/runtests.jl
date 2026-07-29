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

# Phase-12 spatial colocalization map. EVERY Phase-12 testset is aggregated behind the single
# include below, so adding one means adding a line to `test_p12_suite.jl`, never to this file.
#
# IT IS FIRST AMONG THE PHASE-TEST INCLUDES, AND "FIRST" IS THE POINT. A thrown `@testset`
# aborts every remaining include, so any position other than first is contingent on knowing
# which file throws TODAY. That contingency has already failed once: this aggregator was
# originally wired immediately before `test_p13_correction.jl` on the belief that the Phase-13
# correction arm was the only throwing include, and it then never ran -- because
# `include(".../test_npe.jl")` below throws FIRST, on the Phase-4 SC3 (NPE-03) assertion
# `median_speedup > SPEEDUP_GATE`, and aborts everything after it. That shortfall is
# PRE-EXISTING, it is the >100x wall-clock claim rather than a wiring defect, re-deriving the
# bar is a pre-registration decision the user owns, and it is NOT resolved here. First position
# is structural instead of contingent: no other phase's failure can mask Phase 12, whichever
# file happens to throw. `test_p12_consts.jl` testset 9 asserts exactly that ordering.
#
# The aggregator prints one `P12-SUITE-RAN: <file>` marker per include, because `Test` prints
# testset NAMES and none of the Phase-12 names contains its filename -- so only a printed marker
# distinguishes "the include ran" from "the include was silently skipped".
include(joinpath(@__DIR__, "test_p12_suite.jl"))

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
# TWO PHASE-13 FILES ARE NOW DELIBERATELY RED, AND "PUT THE THROWING FILE LAST" NO LONGER WORKS.
#
#   * `test_p13_correction.jl` -- two pre-registered `max |Delta log BF| <= P13_F5_MAXABS_TOL`
#     MEASURED MISSES (0.5498 and 0.3783 against a ceiling of 0.25), committed as-is by plan
#     13-07 rather than tuned away.
#   * `test_p13_real.jl` -- three more, left failing by the user's ruling of 2026-07-29 on plan
#     13-17. On the amended c2/c3 channel pair the D-16 alpha ladder no longer reaches negative
#     m-bar, so the claim that negative induced mu is CONSTRUCTIBLE from real microscopy pixels is
#     RETRACTED. The assertions that encode it stay RED because this phase's own verdict says the
#     three-way Bayes factor works on SIMULATED data and is NOT SHOWN to work on real microscopy:
#     an assertion that states the real-substrate expectation and fails is telling the truth, and
#     a test that passes while the science says otherwise is worse than a red one.
#
# THE ORIGINAL WIRING PUT THE SINGLE THROWING FILE LAST, AND THAT PATTERN SUPPORTS EXACTLY ONE.
# A thrown `@testset` aborts every remaining include, so with two of them whichever runs first
# masks the other -- measured on 2026-07-29: the suite aborted at this file's `test_p13_real.jl`
# line and SEVEN Phase-13 files after it never ran. REORDERING CANNOT FIX THAT. It only chooses
# which red is hidden, and it would be contingent on knowing which file throws today -- the same
# contingency that already failed once above (`test_p12_suite.jl`).
#
# THE FIX IS STRUCTURAL, in the spirit of the first-position rule that made Phase 12 immune, and
# it has four required properties, each verified by running the suite rather than reasoned about:
#
#   1. EVERY Phase-13 sibling RUNS AND REPORTS, whichever known-red file throws. Each known-red
#      include is wrapped, so its throw stops nothing after it.
#   2. THE NAMED-LIMIT FAILURES STILL SURFACE, neither swallowed nor downgraded. `Test` prints
#      every `Test Failed at ...` line and the testset summary table BEFORE the outer testset
#      throws, so catching the throw hides nothing that was already printed. They are NOT
#      `@test_skip`-ed and NOT `@test_broken`-ed; the assertion expressions are byte-unchanged.
#   3. THE SUITE STILL EXITS NON-ZERO. The first recorded exception is RE-RAISED after the last
#      include. A green suite here would be the exact failure the ruling exists to prevent.
#   4. A KNOWN-RED FILE THAT EVER STARTS PASSING FAILS LOUDLY. The observed red set is asserted
#      EQUAL to the expected set below, so an unexpected pass is a test failure that names the
#      file -- a silently-green named limit is how a retracted claim creeps back in.
#
# And one safety property that is not negotiable: anything that is NOT a `Test.TestSetException`
# -- a genuine error, an `UndefVarError`, a load failure -- is RETHROWN IMMEDIATELY and is never
# recorded as an expected red. `include` wraps the throw in a `LoadError` (verified on Julia
# 1.12.6), so the wrapper is unwrapped before the type is tested.
#
# `P13_KNOWN_RED` is declared HERE, after the `test_p12_suite.jl` include, and that position is
# load-bearing: `test_p12_consts.jl` testset 9 locates `test_p13_correction.jl` by FIRST
# occurrence in this file's comment-stripped source and asserts the Phase-12 aggregator precedes
# it. Declaring the list above that include would break an assertion in another phase's file.
const P13_KNOWN_RED      = ["test_p13_real.jl", "test_p13_correction.jl"]
const P13_RED_OBSERVED   = String[]
const P13_RED_EXCEPTIONS = Any[]

function _record_known_red(name::AbstractString, e)
    inner = e isa LoadError ? e.error : e
    inner isa Test.TestSetException || rethrow(e)
    push!(P13_RED_OBSERVED, name)
    push!(P13_RED_EXCEPTIONS, inner)
    return nothing
end

include(joinpath(@__DIR__, "test_p13_consts.jl"))
include(joinpath(@__DIR__, "test_p13_labels.jl"))
include(joinpath(@__DIR__, "test_p13_alpha.jl"))
# KNOWN RED (3 named-limit failures: :226 once, :292 twice). Wrapped, never skipped.
try
    include(joinpath(@__DIR__, "test_p13_real.jl"))
catch e
    _record_known_red("test_p13_real.jl", e)
end
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
# KNOWN RED (2 pre-registered F5 measured misses at :245 and :246). Wrapped, never skipped.
# It stays LAST only because that is its dependency-order position; nothing depends on it being
# last any more, which is the whole point of the ledger below.
try
    include(joinpath(@__DIR__, "test_p13_correction.jl"))
catch e
    _record_known_red("test_p13_correction.jl", e)
end

# THE LEDGER. Runs after every Phase-13 include, and is the property-4 guard: if a known-red file
# ever goes green, the set equality fails and NAMES it. Written as pairs (the house idiom from
# `test_p12_consts.jl`) so the failure output identifies which file changed colour rather than
# printing two anonymous sets.
@testset "Phase-13 known-red ledger -- a named limit cannot go silently green" begin
    for name in P13_KNOWN_RED
        @test (name => name in P13_RED_OBSERVED) == (name => true)
    end
    @test Set(P13_RED_OBSERVED) == Set(P13_KNOWN_RED)
end

# PROPERTY 3. The red must survive to the exit code. Re-raise rather than `exit(1)` so the log
# ends with the real exception rather than a bare status.
if !isempty(P13_RED_EXCEPTIONS)
    println()
    println("PHASE-13 NAMED LIMITS -- THIS SUITE IS RED BY DESIGN. Red files: ",
            join(P13_RED_OBSERVED, ", "))
    println("Five failing assertions: test_p13_real.jl:226 (x1), :292 (x2) -- the RETRACTED ",
            "negative-mu construction claim; test_p13_correction.jl:245, :246 -- the two ",
            "pre-registered F5 misses. See 13-REPORT.md sections 9 (Limit E) and 9a.")
    println("DO NOT clear these to green. The red IS the finding.")
    throw(first(P13_RED_EXCEPTIONS))
end
