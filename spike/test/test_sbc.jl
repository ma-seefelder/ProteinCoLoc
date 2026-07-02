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

# spike/test/test_sbc.jl --- Phase-5 SBC calibration gate (SBC-01..04).
#
# Wave-1 MISSING scaffold (05-01, Task 1): a named, skipped placeholder so the single
# `julia --project=spike spike/test/runtests.jl` gate already enumerates the SBC
# success criteria. This plan's Task 3 REPLACES this scaffold with the real fast
# fixture gate (SC1 ranks M×8; SC2 KS/χ² uniformity + biased-input p≈0 + ECE
# traffic-light; SC3 VAL_MASTER_SEED ≠ NPE_MASTER_SEED anti-snooping; SC4 SBC_CAPTION).

using Test

@testset "Phase 5 — SBC Calibration (Wave 1)" verbose = true begin
    @test_skip true   # filled by Plan 05-01 Task 3 (sbc.jl + figures.jl)
end
