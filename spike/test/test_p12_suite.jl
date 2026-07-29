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
# WHY THE POSITION OF THE ONE WIRING LINE MATTERS. This file is included from `runtests.jl`
# immediately BEFORE `test_p13_correction.jl`, whose outer `@testset` carries two committed
# MEASURED MISSES and therefore THROWS. A thrown testset aborts every remaining include, so a
# sibling placed after it would report green by never running at all.
#
# THE ORDERING RULE IS "BEFORE THE FIRST THROWING INCLUDE", AND `test_p13_correction.jl` IS NOT
# THE ONLY ONE. Measured on 2026-07-29: `include(".../test_npe.jl")` at `runtests.jl:169` also
# throws today -- its SC3 (NPE-03) assertion `median_speedup > SPEEDUP_GATE` evaluates
# 84.44 > 100.0 -- so on this machine the suite aborts there and NOTHING from `runtests.jl:174`
# onward executes, this aggregator and the whole Phase-13 block included. That shortfall is a
# pre-existing, phase-level matter (the >100x wall-clock claim), it is NOT introduced by Phase 12,
# and it is NOT resolved by this wiring. It is recorded here because it is the difference between
# "these testsets are wired" and "these testsets ran".
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
