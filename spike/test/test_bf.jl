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

# spike/test/test_bf.jl --- Phase-5 amortized-Bayes-factor gate (BF-01/02, D-07/D-08).
#
# Wave-1 MISSING scaffold (05-01): a named, skipped placeholder so the single
# `julia --project=spike spike/test/runtests.jl` gate already enumerates the BF
# success criteria the Wave-2 plan (05-02) will turn green. Mirrors test_npe.jl's
# scaffold idiom (license header → `using Test` → one outer `@testset` holding a
# `@test_skip`). Plan 05-02 replaces the skip with the real amortized-logBF-vs-KDE
# reproduction gate scored against the pre-registered BF_CORR_MIN / BF_LOGBF_TOL
# in spike/validation/consts.jl.

using Test

@testset "Phase 5 — Amortized Bayes Factor (Wave 2)" verbose = true begin
    @test_skip true   # filled by Plan 05-02 (bf.jl + train_ratio.jl)
end
