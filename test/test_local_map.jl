#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# 07-08 — windowed sub-tile local colocalization map (`src/amortized/local_map.jl`).
#
# Two things are asserted here:
#   1. `_ensure_grid_registered` populates `_REGISTRY[G]` from a local artifact directory
#      WITHOUT touching `_lazy_load_from_artifact!`/Artifacts (the wave-6 bridge), and is a
#      no-op once the grid is registered.
#   2. `local_coloc_map` returns an r×c map whose Δρ entries are all finite (measured or the
#      documented sentinel) and whose degenerate tiles are flagged rather than fatal.
#
# The CI-portable body builds a TINY grid-4 bundle into a temp dir (the real 8×8 artifacts are
# gitignored regenerable caches). The real 8×8 path is exercised as an OPTIONAL smoke that runs
# only when `artifacts/grid_8/` is present in this checkout.
#########################################################################################

# --- helpers -----------------------------------------------------------------------------

# Persist a tiny, fast-to-train bundle for `G` into `dir`, in exactly the layout
# `_train_grid_pipeline` writes (npe_G / ratio_G / ood_nulls_G).
function _tiny_bundle_artifacts(G::Int, dir::String; n::Int = 48)
    nc   = G^2
    Zraw = vcat(randn(nc, n) .* 2 .+ 3, Float64.(rand(Bool, nc, n)))
    zt   = ProteinCoLoc.fit_summary_transform(Zraw; variant = :min)
    Zstd = Float32.(ProteinCoLoc.standardize_summary(Zraw, zt, :min))
    θ    = randn(7, n)
    npe  = ProteinCoLoc.train_npe(Zstd, Zstd, θ, θ; use_gpu = false, epochs = 2, batchsize = 16,
                                  dstar = 8, depth = 1, width = 16, num_coupling_layers = 2,
                                  flow_depth = 1, flow_width = 8, stopping_epochs = 2)
    ratio = ProteinCoLoc.train_ratio(Zstd, vec(θ[1, :]); n = 32, use_gpu = false, epochs = 2,
                                     batchsize = 16, num_summaries = 8, summary_width = 16,
                                     stopping_epochs = 2)
    nulls = (; density = ProteinCoLoc.fit_ood_nulls(Zstd; variant = :min))
    ProteinCoLoc.save_estimator(joinpath(dir, "npe_$(G).jld2"), npe.estimator, npe.θzt, zt, npe.arch)
    ProteinCoLoc.save_ratio(joinpath(dir, "ratio_$(G).jld2"), ratio)
    ProteinCoLoc.save_ood_nulls(joinpath(dir, "ood_nulls_$(G).jld2"), nulls)
    return dir
end

# EXECUTABLE code of a source file: docstrings (`\"\"\"…\"\"\"`), `#=…=#` blocks and `#` lines are
# stripped, so a mention of a symbol in prose is not mistaken for a definition of it.
function _code_only(path::String)
    s     = read(path, String)
    parts = split(s, "\"\"\"")
    body  = join(parts[1:2:end], "\n")                      # docstrings live in the odd chunks
    body  = replace(body, r"(?s)#=.*?=#" => "")
    return join(filter(l -> !startswith(lstrip(l), "#"), split(body, '\n')), '\n')
end

# A 2-channel synthetic image; strictly positive so every patch clears the ≥15-survivor floor.
function _synth_mci(sz::Tuple{Int,Int}, name::String)
    a = rand(sz...) .+ 0.5
    b = 0.7 .* a .+ 0.3 .* rand(sz...) .+ 0.5
    return MultiChannelImage([a, b], ["c1", "c2"], name, ["", ""], sz, [0.5, 0.5])
end

@testset "windowed sub-tile local map (07-08)" begin
    G   = 4
    dir = mktempdir()
    _tiny_bundle_artifacts(G, dir)

    # --- Task 1: in-process registration bridge -------------------------------------------
    # Earlier testsets in this suite register their own grid-4 bundles; start from an EMPTY
    # slot so the bridge is proven to populate it (and clear it again at the end).
    delete!(ProteinCoLoc._REGISTRY, G)
    @test !haskey(ProteinCoLoc._REGISTRY, G)
    b = ProteinCoLoc._ensure_grid_registered(G; artifact_dir = dir)
    @test haskey(ProteinCoLoc._REGISTRY, G)
    @test b.grid == G
    # estimator_for now resolves WITHOUT the (still-unwired) Artifacts lazy path.
    @test estimator_for(G) === b
    # Idempotent: a second call is a no-op returning the SAME bundle (no reload).
    @test ProteinCoLoc._ensure_grid_registered(G; artifact_dir = dir) === b
    # The bridge must not CALL the Artifacts path (docstrings may name it to explain why not).
    code = _code_only(joinpath(dirname(@__DIR__), "src", "amortized", "local_map.jl"))
    @test !occursin("_lazy_load_from_artifact!", code)
    @test !occursin("Artifacts", code)
    # A missing artifact errors loudly rather than silently falling back.
    @test_throws ErrorException ProteinCoLoc._ensure_grid_registered(99; artifact_dir = dir)

    # --- Task 2: tile-grid shape + per-tile Δρ finiteness ---------------------------------
    img  = _synth_mci((64, 64), "sample")
    ctrl = _synth_mci((64, 64), "control")

    m = local_coloc_map(img, ctrl, [1, 2]; grid = G, tiles = (2, 3), N = 16, artifact_dir = dir)
    @test m isa LocalColocMap
    @test m.grid == G
    @test m.tiles == (2, 3)
    @test size(m.delta_rho) == (2, 3)
    @test size(m.ood_flag)  == (2, 3)
    @test length(m.delta_rho) == 2 * 3            # exactly r×c entries
    @test all(isfinite, m.delta_rho)              # measured or documented sentinel, never NaN/Inf
    @test m.meta.N == 16
    @test size(m.meta.ood_score) == (2, 3)
    @test m.meta.tile_size == (32, 21)            # 64÷2 × 64÷3 (patch() trims the remainder)

    # A degenerate layout (sub-tiles smaller than the patch grid) is SENTINELLED + flagged,
    # never fatal (T-7-04).
    md = local_coloc_map(img, ctrl, [1, 2]; grid = G, tiles = (32, 32), N = 8, artifact_dir = dir)
    @test size(md.delta_rho) == (32, 32)
    @test all(isfinite, md.delta_rho)
    @test all(==(ProteinCoLoc.LOCAL_MAP_SENTINEL), md.delta_rho)
    @test all(md.ood_flag)                        # every sentinel tile is flagged
    @test md.meta.n_sentinel == 32 * 32

    # `channels` must select exactly two channels.
    @test_throws ArgumentError local_coloc_map(img, ctrl, [1]; grid = G, artifact_dir = dir)

    # --- Phase-12 extension point stays clean ---------------------------------------------
    @test !isdefined(ProteinCoLoc, :SpatialColocResult)
    @test !isdefined(ProteinCoLoc, :delta_rho_map)
    @test !isdefined(ProteinCoLoc, :uncertainty_map)
    @test !occursin("SpatialColocResult", code)      # not even referenced in executable code
    @test !occursin(r"delta_rho_map\s*\(", code)     # no definition and no call
    @test !occursin(r"uncertainty_map\s*\(", code)

    delete!(ProteinCoLoc._REGISTRY, G)             # leave the registry as we found it
end

# --- OPTIONAL: the real 8×8 bundle, when this checkout has it -------------------------------
# `artifacts/` is gitignored (regenerable cache), so this smoke is conditional by design: it
# proves the SHIPPED must-have — the gated 8×8 bundle registered in-process from
# artifacts/grid_8/*.jld2 and read per sub-tile — on any machine that has run 07-05.
let root = joinpath(dirname(@__DIR__), "artifacts", "grid_8"),
    have = all(f -> isfile(joinpath(root, f)),
               ("npe_8.jld2", "ratio_8.jld2", "ood_nulls_8.jld2"))
    if have
        @testset "windowed sub-tile local map — real 8×8 bundle (07-08)" begin
            # Earlier testsets register SYNTHETIC grid-8 bundles; drop the slot so the bridge
            # really loads artifacts/grid_8/*.jld2, then restore whatever was there.
            prior = get(ProteinCoLoc._REGISTRY, 8, nothing)
            delete!(ProteinCoLoc._REGISTRY, 8)
            ProteinCoLoc._ensure_grid_registered(8)
            @test haskey(ProteinCoLoc._REGISTRY, 8)
            img  = _synth_mci((256, 256), "sample8")
            ctrl = _synth_mci((256, 256), "control8")
            m = local_coloc_map(img, ctrl, [1, 2]; tiles = (2, 2), N = 32)   # grid = 8 default
            @test m.grid == 8
            @test size(m.delta_rho) == (2, 2)
            @test all(isfinite, m.delta_rho)
            @test m.meta.n_sentinel == 0            # 128² tiles clear the ≥15-survivor floor
            prior === nothing ? delete!(ProteinCoLoc._REGISTRY, 8) :
                                (ProteinCoLoc._REGISTRY[8] = prior)
        end
    else
        @info "test_local_map: artifacts/grid_8 absent — real-8×8 smoke skipped (gitignored cache)"
    end
end
