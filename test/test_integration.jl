#########################################################################################
# test/test_integration.jl --- end-to-end public-path + registry integration (07-10, PROD-01/02)
#
# Exercises the v2.0 primary public API `colocalization_amortized` on the SHIPPED 8×8 grid,
# driving the content-hashed `Artifacts.toml` lazy-load (tree-sha1 verified, T-7-01) through
# `estimator_for(8)` → `_lazy_load_from_artifact!(8)`, and asserts the D-02 accessor interface
# dispatches on the returned `AmortizedColocResult`. Structural assertions (Artifacts.toml shape,
# exports, unshipped/unregistered error paths) always run; the load-dependent assertions run
# whenever the 8×8 bundle is resolvable (artifact store OR the in-repo dev bundle).
#
# Shipped family is {8} ONLY (Go/No-Go Option A). 4/16/32/64 are NOT shipped.
#########################################################################################

# Pkg is imported by runtests.jl; use its bundled TOML + Artifacts so the test needs no extra
# declared dep in the Pkg.test sandbox.
import Pkg

@testset "public API + registry integration (07-10)" begin
    root = dirname(@__DIR__)

    # --- Artifacts.toml: content-hashed, lazy, shipped grids ONLY ---------------------------
    atoml = joinpath(root, "Artifacts.toml")
    @test isfile(atoml)
    a = Pkg.TOML.parsefile(atoml)
    @test haskey(a, "grid_8")
    @test haskey(a["grid_8"], "git-tree-sha1")
    @test a["grid_8"]["lazy"] == true
    for k in ("grid_4", "grid_16", "grid_32", "grid_64")
        @test !haskey(a, k)                        # only 8×8 ships
    end
    # Every entry is content-hashed (git-tree-sha1 present).
    for (k, v) in a
        @test haskey(v, "git-tree-sha1")
    end

    # --- exports / breaking-release surface (D-01) ------------------------------------------
    @test :colocalization_amortized in names(ProteinCoLoc)
    @test ProteinCoLoc._SHIPPED_GRIDS == (8,)
    # The Turing/ADVI reference path stays INTERNAL (not exported).
    @test !(:compute_BayesFactor in names(ProteinCoLoc))
    @test !(:colocalization in names(ProteinCoLoc))

    # --- unshipped / unregistered grid error paths (throw BEFORE any load, T-7-05) ----------
    sz = 256
    mkimg(name) = MultiChannelImage(
        [rand(sz, sz) .+ 0.5, rand(sz, sz) .+ 0.5, rand(sz, sz) .+ 0.5],
        ["c1", "c2", "c3"], name, ["p1", "p2", "p3"], (sz, sz), [0.5, 0.5, 0.5])
    img = mkimg("sample"); ctrl = mkimg("control")

    # grid 64 is not shipped ⇒ ArgumentError (no silent default).
    @test_throws ArgumentError colocalization_amortized(img, ctrl, [1, 2]; num_patches = 64)
    # an unregistered, non-shipped grid points at train_and_register.
    err = try
        colocalization_amortized(img, ctrl, [1, 2]; num_patches = 7)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("train_and_register", err.msg)
    # channel-count validation.
    @test_throws ArgumentError colocalization_amortized(img, ctrl, [1, 2, 3]; num_patches = 8)

    # --- LOAD-DEPENDENT: the full 8×8 public path via the Artifacts lazy-load ----------------
    hash = Pkg.Artifacts.artifact_hash("grid_8", atoml)
    devdir = ProteinCoLoc._dev_bundle_dir(8)
    bundle_available = hash !== nothing && (Pkg.Artifacts.artifact_exists(hash) || isdir(devdir))

    if bundle_available
        # Force the lazy-load-from-artifact path (a prior testset may have registered a stand-in
        # 8×8 bundle; clearing it makes estimator_for(8) resolve the real content-hashed artifact).
        delete!(ProteinCoLoc._REGISTRY, 8)
        b = estimator_for(8)                        # → _lazy_load_from_artifact!(8), tree-sha1 verified
        @test b isa ProteinCoLoc.EstimatorBundle
        @test b.grid == 8
        @test b.zt !== nothing && b.θzt !== nothing
        @test haskey(b.ood_nulls, :density)
        # the recorded pre-registered OOD operating point is wired in (operative density-channel flag).
        @test haskey(b.ood_nulls, :thr)
        # finite single-pass inference off the frozen net.
        Zprobe = ProteinCoLoc.standardize_summary(
            ProteinCoLoc.encode_d01(ProteinCoLoc.patch_summary(
                ProteinCoLoc._pair_mci(img, [1, 2]), 8)), b.zt, :min)
        # 7, NOT the current θ arity: `b` is the SHIPPED grid-8 bundle, whose flow marginal count
        # is frozen in the artifact and read off disk (`persist.load_estimator` uses `arch.D`).
        # Extending the simulator prior does not — and must not — move this number.
        draws = ProteinCoLoc.posterior_for(b.npe, Zprobe; N = 16, use_gpu = false)
        @test size(draws, 1) == 7 && all(isfinite, draws)

        # full public path → AmortizedColocResult with finite outputs + a boolean OOD flag.
        r = colocalization_amortized(img, ctrl, [1, 2]; num_patches = 8, N = 200)
        @test r isa ProteinCoLoc.AmortizedColocResult
        @test r.grid == 8

        # D-02 shared accessor interface dispatches on the result.
        @test size(posterior_draws(r)) == (7, 200)   # frozen shipped-bundle D (see above)
        @test all(isfinite, posterior_draws(r))
        @test length(delta_rho(r)) == 200
        @test all(isfinite, delta_rho(r))
        @test isfinite(bayes_factor(r))
        @test is_ood(r) isa Bool
        @test r.calibration isa ProteinCoLoc.CalibrationMeta
        @test r.calibration.grid == 8
    else
        @info "test_integration: 8×8 bundle not resolvable (no artifact store entry, no in-repo " *
              "dev bundle) — skipping the load-dependent assertions." devdir
        @test_skip false
    end
end
