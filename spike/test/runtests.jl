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
    end

end
