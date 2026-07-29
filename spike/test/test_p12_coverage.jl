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

# spike/test/test_p12_coverage.jl --- P12_PENDING_SCAFFOLD, implemented by plan 12-16.
#
# WHY THE FILE EXISTS BEFORE ITS TESTS DO. It is wired into `test_p12_suite.jl` at wave 1, so it
# already has a GUARANTEED-EXECUTING home ahead of the throwing Phase-13 correction arm. Adding
# a Phase-12 test later then means filling this file, never editing `runtests.jl` -- and a file
# wired late is a file that can silently never run.
#
# IT DELIBERATELY CONTAINS NO FAILING ASSERTION. A failing `@test` throws at the end of its
# enclosing testset, and a thrown testset aborts every remaining include in the aggregator --
# the exact silent-skip trap this wiring exists to prevent. A scaffold must therefore be GREEN;
# "not yet implemented" is carried by the `@info` marker below, not by a red test.

using Test

@testset "P12 per-region coverage and the scoring rule (pending — 12-16)" begin
    @info "P12_PENDING_SCAFFOLD: spike/test/test_p12_coverage.jl is a wave-1 scaffold; plan 12-16 implements the coverage-interval and log-score tests."
    # The one invariant that holds in every state of this phase.
    @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
end
