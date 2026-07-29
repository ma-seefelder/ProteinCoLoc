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

# spike/test/test_p12_suite.jl --- THE SINGLE PHASE-12 INCLUDE POINT.
#
# ADDING A PHASE-12 TEST MEANS ADDING A LINE **HERE**, NEVER TO `runtests.jl`. One aggregator
# with one wiring line is what keeps the include ORDER a property of a single file instead of a
# thing every plan has to re-reason about.
#
# WHY THE POSITION OF THE ONE WIRING LINE MATTERS, AND WHY IT IS **FIRST**. This file is
# included from `runtests.jl` ahead of EVERY other phase-test include. A thrown `@testset`
# aborts every remaining include, so any position other than first is contingent on knowing
# which file throws TODAY -- and that contingency has already failed once, in this very phase.
#
# THE ORIGINAL WIRING PUT THIS AGGREGATOR IMMEDIATELY BEFORE `test_p13_correction.jl`, on the
# documented belief that the Phase-13 correction arm (two committed MEASURED MISSES, so its
# outer testset throws) was the include that aborts the rest. It is not the only one. Measured
# on 2026-07-29: `include(".../test_npe.jl")` throws FIRST -- its Phase-4 SC3 (NPE-03)
# assertion `median_speedup > SPEEDUP_GATE` evaluated 84.44 > 100.0 -- so the suite aborted
# there and NOTHING after it executed: this aggregator, the Phase-5 validation bundle, the
# Phase-9 comparator and the whole Phase-13 block. The aggregator was wired and still never ran.
#
# That Phase-4 shortfall is PRE-EXISTING and is the >100x wall-clock claim rather than a wiring
# defect; re-deriving the bar is a pre-registration decision the user owns, and it is NOT
# resolved here. What IS resolved here is the contingency: first position is STRUCTURAL, so no
# other phase's failure can mask Phase 12 whichever file throws. Testset 9 of
# `test_p12_consts.jl` asserts that ordering against every sibling include, derived from the
# source text rather than from a hardcoded list that would go stale the same way.
#
# THE MARKER LINES BELOW ARE WHAT MAKE THAT DIFFERENCE DETECTABLE, and they exist because the
# obvious check does not work. Julia's `Test` prints TESTSET NAMES, not filenames, and none of
# this phase's ten testset names contains its own filename -- so grepping a suite log for
# `p12_lattice` finds nothing whether or not the file ran, and a verify built that way is green
# on a suite where every Phase-12 include silently never executed. A printed marker is emitted by
# the ACT OF REACHING THE LINE, so its absence is proof the include did not run.
#
# Includes are in dependency order: the frozen pre-registration first, then the numeric kernels,
# the prior, the model surface, the data pipeline, the trainer, the result type, and finally the
# reported-diagnostic and decoupling surfaces.

println("P12-SUITE-RAN: test_p12_consts.jl")
include(joinpath(@__DIR__, "test_p12_consts.jl"))

println("P12-SUITE-RAN: test_p12_lattice.jl")
include(joinpath(@__DIR__, "test_p12_lattice.jl"))

println("P12-SUITE-RAN: test_p12_prior.jl")
include(joinpath(@__DIR__, "test_p12_prior.jl"))

println("P12-SUITE-RAN: test_p12_architecture.jl")
include(joinpath(@__DIR__, "test_p12_architecture.jl"))

println("P12-SUITE-RAN: test_p12_datagen.jl")
include(joinpath(@__DIR__, "test_p12_datagen.jl"))

println("P12-SUITE-RAN: test_p12_train.jl")
include(joinpath(@__DIR__, "test_p12_train.jl"))

println("P12-SUITE-RAN: test_p12_result.jl")
include(joinpath(@__DIR__, "test_p12_result.jl"))

println("P12-SUITE-RAN: test_p12_sbc.jl")
include(joinpath(@__DIR__, "test_p12_sbc.jl"))

println("P12-SUITE-RAN: test_p12_coverage.jl")
include(joinpath(@__DIR__, "test_p12_coverage.jl"))

println("P12-SUITE-RAN: test_p12_decoupling.jl")
include(joinpath(@__DIR__, "test_p12_decoupling.jl"))
