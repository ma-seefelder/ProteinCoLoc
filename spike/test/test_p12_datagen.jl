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

# spike/test/test_p12_datagen.jl --- fixture-scale unit gate for the field-aware training pool.
#
# WHAT THIS FILE DEFENDS, IN ORDER OF WHAT IT WOULD COST TO GET WRONG:
#
#   1. THE theta ROW LAYOUT. `p12_theta_column` is the ONLY place in the phase that turns a draw
#      into a theta column. If its row order drifts from `P12_THETA_ROWS` by ONE, every SBC table,
#      every coverage number and every ablation comparison in Phase 12 is silently about the wrong
#      parameter -- and NOTHING ELSE IN THE SUITE WOULD NOTICE. Testset 1 is the only thing standing
#      between this phase and that failure, so it checks EVERY named row, by name, not a spot check.
#   2. THE STORED GROUND TRUTH. D-07 scores the DRAWN lattice. Testset 2 asserts the stored field is
#      the drawn one and asserts, positively, that it is NOT either rejected alternative.
#   3. THE PHASE-11 POOL. 54 MB of gitignored, unbacked data, already destroyed once, still needed by
#      12-15 and 12-17. Testset 3 keeps the two cache roots provably apart.
#   4. LEAK-FREEDOM. Datagen stores RAW rows; the frozen summary transform is the trainer's
#      train-only fit. Testset 4's third clause is the falsifier that makes the first two mean
#      something.
#
# FIXTURE SCALE ONLY: tiny images, tiny pools, no training, no reported number, and every pool goes
# into a `mktempdir()` so THIS SUITE NEVER WRITES INTO `spike/data/cache/`. Draws ride either the
# fixture stream (`p12_fix_rng` / `P12_FIXTURE_COUNTER`) or the datagen stream keyed by global index,
# never a reported evaluation counter.
#
# HELPERS ARE PHASE-PREFIXED. A bare `_fix_imsize`-style name has already silently overwritten a
# Phase-13 namesake in this suite through an identical signature; every local here carries `p12gen`.
#
# CPU-only. Run:
#     julia --project=spike -t auto spike/test/test_p12_datagen.jl

using Test
using Statistics
using Pkg

isdefined(@__MODULE__, :generate_p12_pool) ||
    include(joinpath(@__DIR__, "..", "data", "p12_generate.jl"))

# --- Fixture knobs. NOT pre-registered thresholds -- Monte-Carlo / wall-clock resolution only. ---
const P12GEN_FIX_IMSIZE = (128, 128)   # tiny on purpose: this file must stay seconds-scale
const P12GEN_FIX_N      = 32           # simulated pool rows for the storage / raw clauses
const P12GEN_FIX_SHARD  = 16           # two shards, so resume-by-skip has something to skip

# A fixed tiny-image sampler, injected so no clause pays the F5 mixture's cost. It is a FIXTURE
# override, not a second copy of the mixture.
_p12gen_fix_imsize(_) = P12GEN_FIX_IMSIZE

"Generate one fixture pool into a FRESH temporary cache root and return `(dir, pool)`."
function _p12gen_pool(n = P12GEN_FIX_N; shard_size = P12GEN_FIX_SHARD, kwargs...)
    root = mktempdir()
    dir  = generate_p12_pool(n; cache_root = root, shard_size = shard_size,
                             imsize_sampler = _p12gen_fix_imsize, imsize_tag = :fixture,
                             verbose = false, kwargs...)
    return (dir, load_p12_pool(dir))
end

# ONE shared pool for the storage / raw / realized clauses, so the file pays for it once.
const P12GEN_DIR, P12GEN_POOL = _p12gen_pool()

@testset "P12 field-aware datagen (D-05, D-07, D-08)" verbose = true begin

    @testset "the theta column matches P12_THETA_ROWS positionally" begin
        # THE ONLY THING BETWEEN THIS PHASE AND AN OFF-BY-ONE ROW LAYOUT. An off-by-one here would
        # fail no other test in the suite: every downstream consumer indexes theta POSITIONALLY, so
        # a shifted layout produces perfectly well-formed tables about the wrong parameters.
        # The fixture stream, never a reported counter.
        draw = sample_p12_prior(p12_fix_rng(P12_FIXTURE_COUNTER))
        col  = p12_theta_column(draw)

        @test length(col) == P12_D_MINISPIKE
        @test P12_D_MINISPIKE == 72

        # Row 1 is the DC coefficient of the GAUSSIAN field in the frozen smoothness order.
        c = p12_dct_vec(vec(draw.z_field))[p12_dct_order(P12_G)]
        @test col[p12_theta_index(:c0)] == c[1]
        # ... and the 63 deviations are the REST of that same permuted vector, in that same order.
        @test col[p12_rows_dev] == c[2:(1 + P12_K_DEV)]

        # The correlation length is APPENDED LAST (D-08), and the chromatic term sits immediately
        # before it -- the exact adjacency a one-row slip would break.
        @test col[p12_theta_index(:r1)] == draw.r1
        @test col[p12_theta_index(:chromatic_eps)] == draw.chromatic_eps
        @test p12_theta_index(:r1) == p12_theta_index(:chromatic_eps) + 1
        @test p12_theta_index(:r1) == P12_D_MINISPIKE

        # EVERY nuisance, BY NAME. Not a spot check: the point is that no single row can slip.
        for sym in (:spillover, :autofluorescence, :label_efficiency,
                    :shift_dx, :shift_dy, :noise, :chromatic_eps)
            @test col[p12_theta_index(sym)] == getproperty(draw, sym)
        end

        # The layout the indices are read from is itself the pre-registered one.
        @test length(P12_THETA_ROWS()) == P12_D_MINISPIKE
        @test P12_THETA_ROWS()[p12_row_c0] === :c0
        @test P12_THETA_ROWS()[p12_row_r1] === :r1
    end

    @testset "the stored ground truth is the DRAWN lattice (D-07)" begin
        G   = P12_G
        n   = length(P12GEN_POOL.r1)
        @test size(P12GEN_POOL.z_field) == (G, G, n)

        # (a) The stored Gaussian field is the DRAWN one, EXACTLY. Both are Float64 and no
        #     transform intervenes between the draw and the shard, so `==` is the right contract.
        for j in (1, 2, n)
            d = sample_p12_prior(p12_datagen_rng(j))
            @test P12GEN_POOL.z_field[:, :, j]   == d.z_field
            @test P12GEN_POOL.rho_field[:, :, j] == d.rho_field
            @test P12GEN_POOL.r1[j]              == d.r1
        end

        # (b) The copula step is checkable after the fact: rho is the ELEMENTWISE ghat of mu.
        @test P12GEN_POOL.rho_field == ghat.(P12GEN_POOL.mu_field)

        # (c) THE REJECTED ALTERNATIVE, ASSERTED AGAINST POSITIVELY. D-07 names the cell mean of the
        #     INTERPOLATED pixel field as the thing NOT to score: the bilinear upsample smooths
        #     across cell boundaries, so a block mean of the rendered field is a different, lossier
        #     quantity than the drawn lattice value. If datagen ever stored that instead, this
        #     clause goes red. (The other rejected alternative -- the induced per-cell patch
        #     correlation -- is a MEASURED quantity and is bounded away from the drawn field by the
        #     whole forward model; clause (a) already pins the stored field to the draw.)
        m  = P12GEN_FIX_IMSIZE[1]
        bs = m ÷ G
        ρ  = P12GEN_POOL.rho_field[:, :, 1]
        up = p12_upsample(ρ, P12GEN_FIX_IMSIZE)
        cellmean = [mean(@view up[(i - 1) * bs + 1:i * bs, (j - 1) * bs + 1:j * bs])
                    for i in 1:G, j in 1:G]
        @test cellmean != ρ
        @test maximum(abs, cellmean .- ρ) > 1e-6
    end

    @testset "the Phase-12 cache root is disjoint from Phase 11's" begin
        # Phase 11's pool is 54 MB of gitignored, unbacked data that cost ~56 min to build, has
        # already been lost once, and is still READ-ONLY input to 12-15 and 12-17.
        @test !occursin(joinpath("cache", "p11"), P12_CACHE_ROOT)
        @test basename(P12_CACHE_ROOT) == "p12"

        tmp = mktempdir()
        # The ARM is in the content hash, so 12-15's controlled CAR-vs-GP comparison cannot read one
        # arm's pool for both (T-12-34).
        @test p12_pool_dir(64; cache_root = tmp, arm = :car) !=
              p12_pool_dir(64; cache_root = tmp, arm = :gp)
        # The r1 PIN is in it too. Without this, 12-11's five ladder rungs collapse into ONE
        # directory, resume-by-skip serves rung 1 to all five, every count assertion still passes,
        # and the gate that authorises the entire training spend reports a fabricated ladder.
        @test p12_pool_dir(64; cache_root = tmp, arm = :car, r1 = 0.2) !=
              p12_pool_dir(64; cache_root = tmp, arm = :car, r1 = 0.8)
        # A pin that is ACCEPTED but not APPLIED is the other half of the same failure.
        _, pinned = _p12gen_pool(16; shard_size = 16, r1 = 0.37)
        @test all(==(0.37), pinned.r1)

        # And nothing this file generated went anywhere near the real cache tree.
        @test !occursin(P12_CACHE_ROOT, P12GEN_DIR)
    end

    @testset "summaries are stored RAW" begin
        S = P12GEN_POOL.summary_min
        @test size(S, 1) == 128
        @test size(S, 2) == length(P12GEN_POOL.r1)
        # Rows 65:128 are the binary present-mask, untouched. A z-score would have moved them.
        @test all(x -> x == 0.0 || x == 1.0, S[65:128, :])
        # THE FALSIFIER THAT MAKES THE TWO CLAUSES ABOVE MEAN SOMETHING. Per-row z-scoring of the
        # pool would drive EVERY continuous row's mean to zero to machine precision. Raw rows have
        # no reason to sit at zero, and none of them does.
        @test all(r -> abs(mean(@view S[r, :])) > 1e-6, 1:64)
    end

    @testset "generation is reproducible and resumable" begin
        # Per-global-index Philox keying means two independent runs must be BYTE-IDENTICAL.
        _, a = _p12gen_pool()
        _, b = _p12gen_pool()
        @test a.theta       == b.theta
        @test a.summary_min == b.summary_min
        @test a.z_field     == b.z_field
        @test a.r1          == b.r1
        @test a.global_index == b.global_index

        # ... and resume-by-skip must not perturb the key stream. Delete the LAST shard, regenerate,
        # and require both that the survivor was NOT rewritten (so the skip really happened) and
        # that the result is identical (so the skip cost nothing).
        root = mktempdir()
        dir  = generate_p12_pool(P12GEN_FIX_N; cache_root = root, shard_size = P12GEN_FIX_SHARD,
                                 imsize_sampler = _p12gen_fix_imsize, imsize_tag = :fixture,
                                 verbose = false)
        nsh  = cld(P12GEN_FIX_N, P12GEN_FIX_SHARD)
        @test nsh >= 2                                    # there IS a shard to skip
        keep_mtime = stat(shard_path(dir, 1)).mtime
        rm(shard_path(dir, nsh))
        @test !p12_pool_complete(dir, P12GEN_FIX_N; shard_size = P12GEN_FIX_SHARD)
        dir2 = generate_p12_pool(P12GEN_FIX_N; cache_root = root, shard_size = P12GEN_FIX_SHARD,
                                 imsize_sampler = _p12gen_fix_imsize, imsize_tag = :fixture,
                                 verbose = false)
        @test dir2 == dir
        @test stat(shard_path(dir, 1)).mtime == keep_mtime   # shard 1 was SKIPPED, not rebuilt
        @test p12_pool_complete(dir, P12GEN_FIX_N; shard_size = P12GEN_FIX_SHARD)
        r = load_p12_pool(dir2)
        @test r.theta == a.theta && r.summary_min == a.summary_min && r.z_field == a.z_field
    end

    @testset "the realized distributions are recorded, not assumed (F6)" begin
        n = length(P12GEN_POOL.r1)
        counts = realized_imsize_counts(P12GEN_POOL.imsize)
        @test sum(values(counts)) == n

        q = realized_r1_quantiles(P12GEN_POOL.r1)
        @test issorted(q.values)
        @test length(q.values) == length(q.probs)
        @test q.values[1] >= P12_R1_MIN && q.values[end] <= P12_R1_MAX

        # THE INTENDED PRIOR IS NOT EVIDENCE OF WHAT A POOL CONTAINS. Both realized records must
        # survive in the manifest, or the provenance claim is only a runtime print.
        meta = JLD2.load(joinpath(P12GEN_DIR, "meta.jld2"))
        @test haskey(meta, "realized_imsize_counts")
        @test haskey(meta, "realized_r1_quantiles")
        @test sum(values(meta["realized_imsize_counts"])) == n
        @test issorted(meta["realized_r1_quantiles"].values)
        # The manifest still carries what `open_or_invalidate` re-checks on every reopen.
        @test haskey(meta, "hash") && haskey(meta, "subhashes") && haskey(meta, "config")
        @test meta["config"].arm === :car && meta["config"].r1_pin === nothing
    end

    @testset "the wall-clock ceiling is a BLOCKER that throws" begin
        # ADDED BEYOND THE PLAN'S SEVEN, because a guard verified only by grepping its own source
        # for the word BLOCKER is not a guard that has been shown to fire. A ceiling of zero minutes
        # must ABORT the run rather than shrink the image-size arm to fit the budget.
        @test_throws ErrorException generate_p12_pool(
            P12GEN_FIX_N; cache_root = mktempdir(), shard_size = P12GEN_FIX_SHARD,
            imsize_sampler = _p12gen_fix_imsize, imsize_tag = :fixture,
            wallclock_ceiling_min = 0.0, verbose = false)
    end

    @testset "datagen ran CPU-only" begin
        @test !haskey(Pkg.project().dependencies, "CUDA")
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
