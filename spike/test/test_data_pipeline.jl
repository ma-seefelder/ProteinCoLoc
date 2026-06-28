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

# Wave W2 unit under test: the in-memory generation core. include() brings the
# UNCHANGED Phase-2 chain (guarded), seeding.jl, encode.jl, and the generator into
# scope. Guarded includes inside generate.jl make this safe both standalone and
# after test_simulator.jl already loaded contract.jl/forward.jl in runtests.jl.
include(joinpath(@__DIR__, "..", "data", "generate.jl"))
using JLD2

# Wave-3 unit under test: the content-hash guard + sharded cache (Tasks 1-2).
# Guarded so this file loads both BEFORE (RED) and AFTER (GREEN) hashguard.jl
# exists, and stays idempotent once generate.jl pulls cache.jl→hashguard.jl in.
const _HASHGUARD = joinpath(@__DIR__, "..", "data", "hashguard.jl")
isfile(_HASHGUARD) && !(@isdefined cache_hash) && include(_HASHGUARD)

# Pre-declared fixture seeds (fixed BEFORE the gate, NOT tuned to pass).
const DP_SC1_SEED    = 7
const DP_REPRO_SEED  = 11
const DP_N_FIXTURE   = 64

@testset "Phase 3 — Training-Data Pipeline" verbose = true begin

    @testset "SC-1 generator → 128-vector" begin
        # Wave W2: generate.jl + encode.jl produce a (128, N) summary_min whose
        # rows 65:128 are the {0,1} mask, mask=0 ⟺ value row==0, all isfinite, and
        # θ is the 7-vector.
        out = generate_samples(DP_N_FIXTURE; master_seed = DP_SC1_SEED, parallel = false)
        sm  = out.summary_min
        @test size(sm) == (128, DP_N_FIXTURE)
        @test size(out.theta) == (7, DP_N_FIXTURE)
        mask = sm[65:128, :]
        vals = sm[1:64, :]
        @test all(x -> x == 0.0 || x == 1.0, mask)        # rows 65:128 ∈ {0,1}
        @test all((mask .== 0.0) .== (vals .== 0.0))       # mask=0 ⟺ value row==0
        @test all(isfinite, sm)                            # every D-01 entry finite
        @test size(out.summary_aug) == (AUG_DIM, DP_N_FIXTURE)
        @test length(out.global_index) == DP_N_FIXTURE
        # surface-don't-swallow (ASVS V5): a bad θ must raise, not corrupt a column.
        badθ = (ρ_true = 2.0, spillover = 0.0, autofluorescence = 0.0,
                label_efficiency = 1.0, shift_dx = 0.0, shift_dy = 0.0, noise = 0.0)
        @test_throws ArgumentError simulate_pair(sample_rng(1, 1), badθ)
    end

    @testset "SC-2 cache round-trip / resume / invalidation" begin
        # Wave W3: cache.jl + hashguard.jl — (a) write→reload array equality,
        # (b) delete one shard, re-run, only it regenerates and the dataset is
        # identical, (c) corrupt the manifest hash → open_or_invalidate refuses it.
        Ncache = 12
        seed   = 7
        ssize  = 5                                   # force 3 shards (5,5,2)
        root   = mktempdir()

        # in-memory reference for the same keyed seed (1:N, column-major)
        ref = generate_samples(Ncache; master_seed = seed, parallel = false)

        # reload every shard in dir and hcat in file order
        _load_pool = function (dir)
            files = sort(filter(f -> startswith(basename(f), "shard_") &&
                                     endswith(f, ".jld2"),
                                readdir(dir; join = true)))
            th = reduce(hcat, [JLD2.load(f, "theta")        for f in files])
            sm = reduce(hcat, [JLD2.load(f, "summary_min")  for f in files])
            sa = reduce(hcat, [JLD2.load(f, "summary_aug")  for f in files])
            gi = reduce(vcat, [JLD2.load(f, "global_index") for f in files])
            return (theta = th, summary_min = sm, summary_aug = sa, global_index = gi)
        end

        # (a) round-trip: reloaded shards == in-memory generator (byte-identical)
        d = generate_cache(root; N = Ncache, master_seed = seed,
                           n_holdout = 20, shard_size = ssize)
        pool = _load_pool(d)
        p = sortperm(pool.global_index)
        @test pool.global_index[p] == collect(1:Ncache)
        @test isequal(pool.theta[:, p],       ref.theta)         # SC-2a
        @test isequal(pool.summary_min[:, p], ref.summary_min)
        @test isequal(pool.summary_aug[:, p], ref.summary_aug)

        # holdout (D-10): separate file, ≥20 samples, global_index disjoint from pool
        @test isfile(joinpath(d, "holdout.jld2"))
        hold = JLD2.load(joinpath(d, "holdout.jld2"))
        @test size(hold["theta"], 2) >= 20
        @test isempty(intersect(Set(pool.global_index), Set(hold["global_index"])))

        # (b) resume-by-skip: delete one shard, re-run, ONLY it regenerates
        victim     = shard_path(d, 1)
        keep       = shard_path(d, 2)
        keep_mtime = mtime(keep)
        rm(victim)
        sleep(0.05)
        d2 = generate_cache(root; N = Ncache, master_seed = seed,
                            n_holdout = 20, shard_size = ssize)
        @test d2 == d                                            # same content-hash dir
        @test isfile(victim)                                     # regenerated
        @test mtime(keep) == keep_mtime                          # untouched shard NOT rewritten
        pool2 = _load_pool(d)
        p2 = sortperm(pool2.global_index)
        @test isequal(pool2.theta[:, p2],       ref.theta)       # SC-2b: dataset identical
        @test isequal(pool2.summary_min[:, p2], ref.summary_min)
        @test isequal(pool2.summary_aug[:, p2], ref.summary_aug)

        # (c) invalidation: cache_hash config-sensitive + corrupted manifest refused
        cfg_a = generating_config(Ncache; master_seed = seed,     k = 5, shard_size = ssize)
        cfg_b = generating_config(Ncache; master_seed = seed + 1, k = 5, shard_size = ssize)
        @test cache_hash(cfg_a) != cache_hash(cfg_b)             # SC-2c precondition
        meta_path = joinpath(d, "meta.jld2")
        rm(meta_path)
        jldsave(meta_path; hash = "deadbeef", subhashes = Dict{String,String}(),
                config = cfg_a, N = Ncache, shard_size = ssize, schema_version = 1)
        @test_throws ErrorException open_or_invalidate(root, cfg_a)   # SC-2c: stale refused
    end

    @testset "committed cache fixture (on-disk schema stability)" begin
        # A tiny fixture cache (un-ignored via .gitignore negation) pins the
        # column-major on-disk schema across runs; bulk cache stays ignored.
        fshard = joinpath(@__DIR__, "..", "data", "cache", "fixture", "shard_0001.jld2")
        @test isfile(fshard)
        s = JLD2.load(fshard)
        @test size(s["theta"], 1)       == 7
        @test size(s["summary_min"], 1) == 128
        @test size(s["summary_aug"], 1) == AUG_DIM
        @test s["schema_version"]       == SCHEMA_VERSION
        sm = s["summary_min"]
        @test all(x -> x == 0.0 || x == 1.0, sm[65:128, :])      # mask rows ∈ {0,1}
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
        # Wave W2: seeding.jl + generate.jl — the parallel path equals the serial
        # fallback byte-identically, and a shuffled generation order re-sorts to
        # identical columns (Philox4x keyed per global index, no shared state).
        a = generate_samples(DP_N_FIXTURE; master_seed = DP_REPRO_SEED, parallel = false)
        b = generate_samples(DP_N_FIXTURE; master_seed = DP_REPRO_SEED, parallel = true)
        @test isequal(a.theta,        b.theta)
        @test isequal(a.summary_min,  b.summary_min)
        @test isequal(a.summary_aug,  b.summary_aug)
        @test a.global_index == b.global_index
        @test a.imsize       == b.imsize

        # Shuffled-index generation re-sorted by global index == in-order generation.
        perm = shuffle(MersenneTwister(99), collect(1:DP_N_FIXTURE))
        c    = generate_samples(DP_N_FIXTURE; master_seed = DP_REPRO_SEED,
                                indices = perm, parallel = false)
        invp = sortperm(perm)                                # invp[k] = position of global idx k
        @test c.global_index[invp] == a.global_index
        @test isequal(c.theta[:, invp],       a.theta)
        @test isequal(c.summary_min[:, invp], a.summary_min)
        @test isequal(c.summary_aug[:, invp], a.summary_aug)
    end

    @testset "hashguard content-hash version guard (D-05)" begin
        # Wave W3 / Task 1: cache_hash is deterministic, config-sensitive, and
        # source-sensitive with per-file sub-hashes for diagnosability.
        if !(@isdefined cache_hash)
            @test false  # RED: spike/data/hashguard.jl not yet implemented
        else
            c1 = (N = 64, master_seed = 1, k = 5)
            c2 = (N = 64, master_seed = 2, k = 5)
            @test cache_hash(c1) == cache_hash(c1)         # deterministic
            @test cache_hash(c1) != cache_hash(c2)         # any config field flips it
            sh = subhashes()
            @test length(sh) == length(HASH_SRC_FILES)     # one sub-hash per source
            # Source-byte sensitivity + diagnosability on an ISOLATED temp file
            # (never perturbs the real frozen sources).
            mktempdir() do td
                f = joinpath(td, "src.jl"); write(f, "alpha")
                cfg = (N = 1,)
                h1 = cache_hash(cfg; src_files = [f]); s1 = subhashes([f])
                write(f, "alphabet")                        # perturb source bytes
                h2 = cache_hash(cfg; src_files = [f]); s2 = subhashes([f])
                @test h1 != h2                              # source byte flips the hash
                @test s1[f] != s2[f]                        # the changed file's sub-hash moved
            end
        end
    end

end
