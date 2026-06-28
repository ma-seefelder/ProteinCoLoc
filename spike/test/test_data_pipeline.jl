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

# spike/test/test_data_pipeline.jl --- Phase-3 DATA-01..03 gates.
#
# Wave-0 MISSING scaffold (this plan, 03-01): the success-criterion testsets for
# the training-data pipeline exist as named, skipped placeholders so the single
# `julia --project=spike spike/test/runtests.jl` gate already enumerates every SC
# the later waves must turn green. Each child testset holds one `@test_skip true`
# tagged with the wave that replaces it with the real gate (RESEARCH §Validation
# Architecture; PATTERNS §spike/test/test_data_pipeline.jl). Mirrors the structure
# of test_simulator.jl exactly (license header -> using Test/Random -> nested
# `@testset ... verbose = true`).

using Test
using Random

@testset "Phase 3 — Training-Data Pipeline" verbose = true begin

    @testset "SC-1 generator → 128-vector" begin
        # Wave W2: generate.jl + encode.jl produce a (128, N) summary_min whose
        # rows 65:128 are the {0,1} mask, mask=0 ⟺ value row==0, all isfinite, and
        # θ is the 7-vector. Replaces this placeholder with the real SC-1 gate.
        @test_skip true
    end

    @testset "SC-2 cache round-trip / resume / invalidation" begin
        # Wave W3: cache.jl + hashguard.jl — (a) write→reload array equality,
        # (b) delete one shard, re-run, only it regenerates and the dataset is
        # identical, (c) perturb a tracked source byte/config field → cache_hash
        # differs → loader refuses the stale cache.
        @test_skip true
    end

    @testset "SC-3 leak-free standardization" begin
        # Wave W3: loader.jl — no public `standardize_all`; Ztr rows ≈0 mean/unit
        # std, Zva NOT exactly 0/1 (fit-on-train only), binary mask rows unchanged.
        @test_skip true
    end

    @testset "SC-4 k-fold + reserved holdout" begin
        # Wave W3: loader.jl — union(folds)==1:N with pairwise ∩==∅, holdout
        # disjoint from every fold, length(holdout) ≥ 20, reproducible across two
        # loader calls.
        @test_skip true
    end

    @testset "D-11/D-12 order/thread independence" begin
        # Wave W2: seeding.jl + generate.jl — generate with nthreads()==1 and >1
        # yield byte-identical θ + summaries (counter-based Philox4x keyed per
        # global index, no shared mutable state).
        @test_skip true
    end

end
