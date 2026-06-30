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
using Statistics      # mean, std (SC-3 standardization assertions)

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

# Wave-4 unit under test: the leak-free k-fold loader (loader.jl, Plan 03-04). It is
# a `module Loader` with a trailing `using .Loader`, so load_fold/load_main_pool/
# load_holdout/all_folds are directly callable AND `names(Loader)` is the controlled
# D-08 public surface the SC-3a assertion introspects. Guarded for idempotent loading.
isdefined(@__MODULE__, :load_fold) || include(joinpath(@__DIR__, "..", "data", "loader.jl"))

# Pre-declared fixture seeds (fixed BEFORE the gate, NOT tuned to pass).
const DP_SC1_SEED    = 7
const DP_REPRO_SEED  = 11
const DP_N_FIXTURE   = 64

# Loader (SC-3/SC-4) fixture parameters + tolerances (pre-declared, NOT tuned).
const DP_LOAD_SEED = 3
const DP_N_LOAD    = 40      # 40 / K=5 ⇒ 8-column val folds; one shard at SHARD_SIZE
const DP_K         = 5
const DP_MEAN_TOL  = 1e-6    # Float32 standardized train mean ≈ 0 (worst case ~3.6e-7)
const DP_STD_TOL   = 1e-3    # corrected std of standardized train ≈ 1
# One tiny cache shared by the loader gates (keeps the quick gate sub-30s).
const DP_LOAD_DIR  = generate_cache(mktempdir(); N = DP_N_LOAD,
                                    master_seed = DP_LOAD_SEED, n_holdout = 20)

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
        # Wave W4: loader.jl — no public `standardize_all` (D-08); Ztr continuous
        # rows ≈0 mean/unit std (fit-on-train, D-07), Zva NOT exactly 0/1 (val never
        # fit), binary mask rows (65:128) passed through UNCHANGED (mask bypass).
        d    = DP_LOAD_DIR
        f    = load_fold(d, 1; K = DP_K, master_seed = DP_LOAD_SEED)
        pool = load_main_pool(d)

        # SC-3a / D-08: NO public global-standardize symbol — leakage impossible by
        # construction. `names(Loader)` is the controlled public surface.
        @test !(:standardize_all in names(Loader))

        # SC-3b: fit-on-train ⇒ Ztr continuous rows (1:64) ≈ 0 mean / unit std per row
        @test all(abs.(mean(f.Ztr[1:64, :], dims = 2)) .< DP_MEAN_TOL)
        @test all(abs.(std(Float64.(f.Ztr[1:64, :]), dims = 2) .- 1) .< DP_STD_TOL)

        # SC-3b: val never used to fit ⇒ Zva continuous rows NOT all exactly 0/1
        @test !all(x -> x == 0f0 || x == 1f0, f.Zva[1:64, :])

        # SC-3b: mask bypass — rows 65:128 byte-identical to the RAW mask rows. Recompute
        # the train/val split with the SAME keyed fold RNG to address the exact columns.
        N         = size(pool.summary_min, 2)
        perm      = randperm(fold_rng(DP_LOAD_SEED), N)
        val_idx   = perm[1:DP_K:end]
        train_idx = setdiff(1:N, val_idx)
        @test f.Ztr[65:128, :] == Float32.(pool.summary_min[65:128, train_idx])
        @test f.Zva[65:128, :] == Float32.(pool.summary_min[65:128, val_idx])
    end

    @testset "SC-4 k-fold + reserved holdout" begin
        # Wave W4: loader.jl — union(folds)==1:N with pairwise ∩==∅ (D-09), the
        # reserved ≥20 holdout excluded from every fold STRUCTURALLY (D-10), and the
        # split reproducible across two loader calls.
        d     = DP_LOAD_DIR
        folds = all_folds(DP_N_LOAD; K = DP_K, master_seed = DP_LOAD_SEED)

        # SC-4a: disjoint + complete — union is exactly 1:N, every pair empty ∩
        @test length(folds) == DP_K
        @test sort(vcat(folds...)) == collect(1:DP_N_LOAD)
        for i in 1:DP_K, j in (i + 1):DP_K
            @test isempty(intersect(folds[i], folds[j]))
        end

        # SC-4a: reproducible (D-09) — same master_seed ⇒ identical membership
        @test all_folds(DP_N_LOAD; K = DP_K, master_seed = DP_LOAD_SEED) == folds
        # Two INDEPENDENT load_fold calls (saved to separate vars) — a genuine D-09
        # reproducibility check, not the tautological f(x)==f(x): same master_seed
        # must yield identical fold membership, a different seed a different split.
        fa = load_fold(d, 2; K = DP_K, master_seed = DP_LOAD_SEED)
        fb = load_fold(d, 2; K = DP_K, master_seed = DP_LOAD_SEED)
        fc = load_fold(d, 2; K = DP_K, master_seed = DP_LOAD_SEED + 1)
        @test fa.θva == fb.θva     # same seed → identical fold membership (D-09)
        @test fa.θva != fc.θva     # different seed → different split

        # SC-4b: reserved holdout present, ≥20 stacks, in its OWN file
        hold = load_holdout(d)
        @test size(hold.theta, 2) >= 20

        # SC-4b STRUCTURAL exclusion (D-10) — NOT an integer-set intersection (holdout
        # and fold indices are both 1:N-style ranges that overlap by value):
        #   (1) the main pool counts EXACTLY the N non-holdout samples, never N+H
        @test size(load_main_pool(d).theta, 2) == DP_N_LOAD
        #   (2) load_main_pool never READS holdout.jld2 — hiding the file leaves the
        #       pool column count unchanged (it only globs shard_*.jld2)
        hp = joinpath(d, "holdout.jld2")
        @test isfile(hp)
        n_before = size(load_main_pool(d).theta, 2)
        mv(hp, hp * ".hidden")
        n_after = size(load_main_pool(d).theta, 2)
        mv(hp * ".hidden", hp)
        @test n_before == n_after == DP_N_LOAD
        #   (3) the holdout RNG key namespace is XOR-salted disjoint from the main pool
        @test HOLDOUT_SALT != 0
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
