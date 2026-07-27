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

# spike/test/test_p13_result.jl --- the D-08 three-hypothesis result-type gate.
#
# WHAT THIS FILE IS FOR. D-08 makes three claims that are cheap to state and easy to lose, so
# each one is asserted here rather than left to the docstrings:
#
#   1. THE TYPE IS A REAL SUBTYPE. Phase 7 D-02 requires every new result variant to slot in as a
#      new `AbstractColocResult` subtype, never as bolted-on fields. Defining the executable type
#      in spike/ (D-01) must not weaken that -- so the subtype relation is asserted against the
#      ACTUAL supertype from src/results.jl, reached read-only.
#   2. THE REFERENCE CLASS IS SELF-DOCUMENTING AND STRUCTURALLY ZERO. The sketch's positional
#      triple put the reference class FIRST, which a caller indexing `[1]` for "the colocalization
#      evidence" would get silently wrong. The NamedTuple removes that hazard, and the zero on the
#      reference entry is enforced by the constructor rather than assumed by convention.
#   3. THE LOSSY ACCESSOR IS DOCUMENTED AND THE MISSING ONE DOES NOT LIE. `bayes_factor` returns
#      ONE number out of a three-way result, so the choice must be retrievable from the running
#      system (asserted through `@doc`, not by reading the source); and `delta_rho` must fall
#      through to the interface error rather than fabricate a contrast no control posterior backs.
#
# Mirrors test_p13_net.jl: license header -> using Test -> guarded include of the unit under test
# -> fixtures computed ONCE at module level -> one outer @testset. No training, no simulation, no
# figure call, and no reported stream is consumed (nothing here draws random numbers at all).

using Test

# Unit under test: the spike-lane result type (pulls src/results.jl read-only, p13/consts.jl and
# validation/sbc.jl transitively).
isdefined(@__MODULE__, :ThreeHypothesisColocResult) ||
    include(joinpath(@__DIR__, "..", "p13", "result.jl"))

# The EXECUTABLE source of the unit under test, with whole-line comments stripped BEFORE any
# count, so a grep assertion cannot be satisfied (or defeated) by prose. The banner in result.jl
# names the misnomer on purpose -- that naming is the file's documentation of what it replaced,
# and it must not be mistaken for the identifier surviving in code.
const P13_RESULT_SRC  = read(joinpath(@__DIR__, "..", "p13", "result.jl"), String)
const P13_RESULT_CODE = join(filter(l -> !startswith(strip(l), "#"),
                                    split(P13_RESULT_SRC, '\n')), '\n')

# GUARDED, NOT REDEFINED. `P13_REPO_ROOT` is already a `const` in spike/p13/real_images.jl:146,
# and the spike lane is a FLAT top-level namespace: when the whole suite runs in one process the
# two definitions land in the same module. The values are identical (both normalise to the repo
# root), but re-`const`-ing a name is exactly the collision class Phase 13's wiring plan warns
# about, so this file reuses the existing binding when one is present and defines it only when
# running standalone.
if !isdefined(@__MODULE__, :P13_REPO_ROOT)
    const P13_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
end

# --- Fixtures built ONCE --------------------------------------------------------------------
# A perfectly calibrated toy reliability input, only so the result has a real CalibrationMeta to
# carry. Nothing in this file gates on its numbers -- test_p13_calibration.jl owns that.
const RES_FIX_CAL  = _bin_calibration(fill(0.5, 100), vcat(trues(50), falses(50));
                                      n_bins = P13_ECE_NBINS)
const RES_FIX_META = p13_calibration_meta(RES_FIX_CAL; grid = 8, auc = 0.97)
const RES_FIX_OOD  = OODVerdict(1.23, false, (; density = 1.23, noise = 0.4))

"Build a valid three-hypothesis result; `meta` is where control draws live, if any."
_res(; coloc = 2.5, exclusion = -3.25, random = 0.0, meta = NamedTuple()) =
    ThreeHypothesisColocResult(8, zeros(7, 5),
                               (coloc = coloc, random = random, exclusion = exclusion),
                               RES_FIX_OOD, RES_FIX_META, meta)

@testset "P13 three-hypothesis result type (D-08)" verbose = true begin

    @testset "is a new AbstractColocResult subtype (Phase-7 D-02)" begin
        # The supertype is the REAL one from src/results.jl, reached read-only -- not a spike-local
        # look-alike. If the read-only include ever broke, this assertion is what would notice.
        @test ThreeHypothesisColocResult <: AbstractColocResult
        @test isabstracttype(AbstractColocResult)
        @test isconcretetype(ThreeHypothesisColocResult)
        @test _res() isa AbstractColocResult
        # D-02's other half: the new variant does NOT bolt fields onto the shipped struct.
        @test fieldnames(AmortizedColocResult) == (:grid, :posterior, :delta_rho_draws,
                                                   :log_bayes_factor, :ood, :calibration, :meta)
    end

    @testset "the random entry is structurally zero" begin
        r = _res()
        # `===` and not `==`: -0.0 == 0.0 is true but they are different bit patterns, and the
        # claim is that the reference class was scored against itself, giving exactly +0.0.
        @test log_bf_vs_random(r).random === 0.0
        @test r.log_bf_vs_random.random === 0.0
        # A non-zero reference entry means the triple was built against some other reference (or a
        # third measurement was fabricated), so it is rejected at construction, not downstream.
        @test_throws ArgumentError _res(random = 0.4)
        @test_throws ArgumentError _res(random = -1e-12)
        # The error names D-08 so a reader hitting it lands on the decision, not on a guess.
        err = try; _res(random = 0.4); catch e; e; end
        @test occursin("D-08", err.msg)
        # A triple missing a key is rejected too, rather than silently defaulting.
        @test_throws ArgumentError ThreeHypothesisColocResult(
            8, zeros(7, 5), (coloc = 1.0, random = 0.0), RES_FIX_OOD, RES_FIX_META, NamedTuple())
    end

    @testset "the reference class is self-documenting" begin
        r = _res(coloc = 2.5, exclusion = -3.25)
        # THE HAZARD THIS REMOVES: the sketch at src/results.jl:169-178 stores an
        # `NTuple{3,Float64}` documented as `{null, coloc, anti-coloc}` -- the REFERENCE CLASS
        # FIRST. A caller who reads a name containing "simplex" and indexes `[1]` expecting "the
        # colocalization evidence" gets the structural zero instead, silently. With a NamedTuple
        # every read site is BY NAME, so no reordering can change a number.
        @test keys(log_bf_vs_random(r)) == (:coloc, :random, :exclusion)
        @test log_bf_vs_random(r).coloc     == 2.5
        @test log_bf_vs_random(r).exclusion == -3.25
        @test log_bf_vs_random(r) isa ThreeWayLogBF
        # The read surface in spike/p13/net.jl returns its triple as (coloc, exclusion, random).
        # Storing that positionally would have swapped two entries; the constructor normalises by
        # NAME instead, so the two files cannot drift into a silent coupling.
        r2 = ThreeHypothesisColocResult(8, zeros(7, 5),
                                        (coloc = 2.5, exclusion = -3.25, random = 0.0),
                                        RES_FIX_OOD, RES_FIX_META, NamedTuple())
        @test log_bf_vs_random(r2) == log_bf_vs_random(r)
        @test log_bf_vs_random(r2).exclusion == -3.25
    end

    @testset "bayes_factor returns the coloc entry, documented" begin
        r = _res(coloc = 2.5, exclusion = -3.25)
        @test bayes_factor(r) == log_bf_vs_random(r).coloc
        @test bayes_factor(r) == 2.5
        # ... and NOT the exclusion entry, and not some fused scalar.
        @test bayes_factor(r) != log_bf_vs_random(r).exclusion
        # THE CHOICE MUST NOT BE MAKEABLE SILENTLY. A single-number accessor on a three-way result
        # is inherently lossy, so the running system -- not just the source file -- has to be able
        # to tell a caller which number they got.
        doc = string(@doc(bayes_factor))
        @test !isempty(doc)
        @test occursin("log_bf_vs_random", doc)
        @test occursin("coloc", doc)
    end

    @testset "delta_rho does not fabricate" begin
        # No control posterior was carried, so there is no Delta rho to return. Synthesising one
        # (zeros, the coloc draws, a difference against a notional zero control) would put a number
        # in front of a reader that no measurement backs; the documented interface error is the
        # honest answer.
        r = _res()
        @test_throws ErrorException delta_rho(r)
        err = try; delta_rho(r); catch e; e; end
        @test occursin("accessor interface", err.msg)
        @test occursin("delta_rho", err.msg)
        @test occursin("ThreeHypothesisColocResult", err.msg)
        # When the run DID carry control draws they are returned unchanged -- the guard is about
        # absence, not a blanket refusal.
        draws = [0.11, -0.04, 0.37]
        rc = _res(meta = (; delta_rho_draws = draws))
        @test delta_rho(rc) == draws
        # The other three accessors are implemented unconditionally.
        @test posterior_draws(r) == zeros(7, 5)
        @test is_ood(r) === false
        @test is_ood(ThreeHypothesisColocResult(8, zeros(7, 5),
                                                (coloc = 0.0, random = 0.0, exclusion = 0.0),
                                                OODVerdict(9.9, true, (; density = 9.9)),
                                                RES_FIX_META, NamedTuple())) === true
    end

    @testset "no simplex naming survives" begin
        # Comment lines are stripped first: result.jl's banner NAMES the misnomer, because
        # documenting what was replaced is the point of the banner. What must not survive is the
        # identifier in EXECUTABLE code, where a reader could copy it forward.
        @test !occursin("log_bf_simplex", P13_RESULT_CODE)
        @test !occursin("bayes_factor_simplex", P13_RESULT_CODE)
        # The replacement is actually present, so the assertions above cannot pass by the file
        # simply not mentioning evidence at all.
        @test occursin("log_bf_vs_random", P13_RESULT_CODE)
        # The banner does discuss the misnomer -- assert that too, so a future edit cannot delete
        # the explanation and leave a bare rename with no recorded reason.
        @test occursin("log_bf_simplex", P13_RESULT_SRC)
        @test occursin("MISNOMER", P13_RESULT_SRC)
    end

    @testset "src/results.jl is untouched (D-01, rename deferred)" begin
        # D-01 keeps all Phase-13 work in spike/, and CLAUDE.md requires src/ to stay provably
        # untouched during the spike. The rename of the sketch's field is therefore RECORDED
        # (P13_RESULTS_RENAME_DEFERRED) and carried as a one-line productionization item, NOT
        # performed here. This assertion is what makes "provably" mean something.
        @test P13_RESULTS_RENAME_DEFERRED == true
        @test P13_SRC_UNTOUCHED == true
        # Guarded: a sandbox without git, or a source tree that is not a checkout, skips the
        # executable half rather than failing for an unrelated reason.
        _has_git = Sys.which("git") !== nothing &&
                   (isdir(joinpath(P13_REPO_ROOT, ".git")) ||
                    isfile(joinpath(P13_REPO_ROOT, ".git")))
        if _has_git
            @test success(Cmd(`git diff --quiet HEAD -- src/results.jl`; dir = P13_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P13_REPO_ROOT))
        else
            @info "git unavailable — skipping the executable src/ byte-equality assertion"
        end
        # The sketch block is still there, still comments, still carrying the old name -- so the
        # phase report's "recorded, not performed" claim is checkable from the running suite.
        _sketch = read(joinpath(P13_REPO_ROOT, "src", "results.jl"), String)
        @test occursin("log_bf_simplex", _sketch)
        @test occursin("SKETCH ONLY", _sketch)
    end

    @testset "the result type ran CPU-only (D-04)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
