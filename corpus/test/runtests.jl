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

# corpus/test/runtests.jl --- the SINGLE offline gate for the external corpus module (SC3).
#
# One command is the whole-module gate:
#     julia --project=. corpus/test/runtests.jl
# It runs with ZERO network access. Every downstream Phase-8 plan (hash, fetch, manifest,
# load, CBS, anchors) wires its testset into this runner via the commented include stubs
# below, so the corpus always has exactly one green gate. Any live-network smoke a later
# plan adds MUST skip cleanly offline (reuse the Phase-9 tapqir_anchor :skipped contract).
#
# DECOUPLING (hard constraint, CLAUDE.md): reaches src/LoadImages.jl + src/colocalization.jl
# READ-ONLY via include() to prove the fixture converts through the UNCHANGED summary path.
# It edits no src/ and touches no root Project.toml.

using Test

# Statistics/StatsBase names that src/colocalization.jl `correlation` resolves at call time,
# and Images which src/LoadImages.jl references qualified. Bringing them into scope here lets us
# include the (module-free) src files standalone without loading the full ProteinCoLoc module
# (which would pull GLMakie and other heavy, display-bound deps).
import Images
import Statistics: mean, cor
import StatsBase: corspearman, corkendall

# Pre-declared corpus contract consts (Task 1).
include(joinpath(@__DIR__, "..", "config.jl"))

# Read-only reach into the FROZEN src/ summary path (no edits — CLAUDE.md decoupling).
include(joinpath(@__DIR__, "..", "..", "src", "LoadImages.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "colocalization.jl"))

# Deterministic offline fixture generator (Task 2).
include(joinpath(@__DIR__, "make_fixture.jl"))

const FIXTURES_DIR = joinpath(@__DIR__, "fixtures")

@testset "corpus (offline)" verbose = true begin

    @testset "fixture round-trip + determinism (SC3)" begin
        # (a) Generate into a scratch dir and assert both channels round-trip via load_tiff.
        tmp1 = mktempdir()
        p1, p2 = make_fixture(tmp1)
        @test isfile(p1) && isfile(p2)
        m1 = load_tiff(p1)
        m2 = load_tiff(p2)
        @test m1 isa Matrix{Float64}
        @test m2 isa Matrix{Float64}
        @test size(m1) == (64, 64)
        @test all(x -> x > 0, m1) && all(x -> x > 0, m2)   # strictly positive ⇒ Pitfall-6 safe

        # (b) Determinism: a second generation with the same seed is byte-identical.
        tmp2 = mktempdir()
        q1, q2 = make_fixture(tmp2)
        @test read(p1) == read(q1)
        @test read(p2) == read(q2)
    end

    @testset "fixture conversion (SC1/D-03)" begin
        # Ensure the committed fixture exists (regenerate deterministically if a fresh checkout
        # lacks it), then convert it through the UNCHANGED src/ summary path to an 8x8 grid.
        fp1 = joinpath(FIXTURES_DIR, "fixture_ch1.tif")
        fp2 = joinpath(FIXTURES_DIR, "fixture_ch2.tif")
        if !(isfile(fp1) && isfile(fp2))
            make_fixture(FIXTURES_DIR)
        end

        img = MultiChannelImage("corpus_fixture", String[fp1, fp2], String["ch1", "ch2"])
        @test length(img.data) == 2
        @test img.data[1] isa Matrix{Float64}

        # patch(·, 8) + correlation(·) from src/colocalization.jl apply verbatim.
        grid1 = patch(img.data[1], 8)
        grid2 = patch(img.data[2], 8)
        rho   = correlation(grid1, grid2)

        @test size(rho) == (8, 8)
        # Pitfall-6 guard: NOT an all-`missing` grid, and the mean over valid patches is finite.
        @test count(!ismissing, rho) > 0
        @test isfinite(mean(skipmissing(rho)))
    end

    # --- Downstream include stubs (ONE command stays the gate) ----------------------------
    # Later Phase-8 plans wire their testsets in here. Each MUST run offline (any live-network
    # smoke skips cleanly via the Phase-9 tapqir_anchor :skipped contract).
    include(joinpath(@__DIR__, "test_hash.jl"))         # wired by 08-02 (SHA-256 hashing)
    include(joinpath(@__DIR__, "test_fetch.jl"))        # wired by 08-02 (fetch: skip-with-flag / hard-fail)
    include(joinpath(@__DIR__, "test_manifest.jl"))     # wired by 08-03 (schema + tier/split guards)
    include(joinpath(@__DIR__, "test_load.jl"))         # wired by 08-03 (TIFF → MultiChannelImage)
    include(joinpath(@__DIR__, "test_cbs.jl"))          # wired by 08-04 (CBS ingestion, simulated-secondary)
    # include(joinpath(@__DIR__, "test_anchors.jl"))    # wired by 08-07 (physical anchors, sealed holdout)
end
