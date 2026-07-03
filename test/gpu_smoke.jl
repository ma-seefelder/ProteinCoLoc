#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gpu_smoke.jl --- GPU-train smoke with graceful CPU fallback (D-06, T-7-02 A7).
#
# Mirrors the ENV-02 CPU smoke but for the TRAINING path. It proves the D-06 reproducibility
# split's two halves:
#
#   (a) TRAIN MAY USE GPU: `train_npe(...; use_gpu = true)` runs on a CUDA device when one is
#       visible, and — per the CLAUDE.md graceful-fallback constraint — DEGRADES CLEANLY to CPU
#       when CUDA is absent (Flux's `gpu` is a no-op with no CUDA extension loaded), with NO error.
#       `has_cuda_device()` branches the assertion so a CPU-only machine exercises ONLY the
#       fallback path and still passes.
#
#   (b) THE FROZEN NET PERSISTS CPU-RESIDENT: `save_estimator` writes `Flux.state(cpu(est))`
#       (a device-independent nested-array snapshot, NOT a whole PosteriorEstimator), so the
#       artifact reloads on a machine with no CUDA and CPU inference runs — the pre-registered
#       ship-gate + shipped inference stay reproducible regardless of the training hardware.
#
# Standalone-runnable: `julia --project -e 'include("test/gpu_smoke.jl")'`; also included by
# runtests.jl. NOTE: this machine has no CUDA device (has_cuda_device() == false), so the reported
# path here is the graceful `use_gpu=true → CPU fallback`.

using ProteinCoLoc
using Test
import StatsBase
import JLD2

@testset "GPU-train smoke (graceful CPU fallback, D-06)" begin
    # A tiny grid-4 fixture (2·G² = 32-row summary), no forward simulator needed.
    G = 4; nc = G^2; d = 2 * nc; n = 40
    Zraw = vcat(randn(nc, n) .* 2 .+ 3, Float64.(rand(Bool, nc, n)))
    zt   = ProteinCoLoc.fit_summary_transform(Zraw; variant = :min)
    Zstd = Float32.(ProteinCoLoc.standardize_summary(Zraw, zt, :min))
    θ    = randn(7, n)

    tiny = (; batchsize = 16, dstar = 8, depth = 1, width = 16,
            num_coupling_layers = 2, flow_depth = 1, flow_width = 8, stopping_epochs = 2)

    # (a) use_gpu = true: runs on GPU if present, else degrades to CPU with NO error.
    res = ProteinCoLoc.train_npe(Zstd, Zstd, θ, θ; use_gpu = true, epochs = 2, tiny...)
    @test res.d_in == d
    @test length(res.θzt.mean) == 7
    # posterior draws are reachable through the CPU read surface (7×N).
    @test size(ProteinCoLoc.posterior_for(res.estimator, Zstd[:, 1]; N = 8, use_gpu = false), 1) == 7

    # (b) the frozen net persists CPU-RESIDENT: Flux.state key (not a whole object), reloads with
    #     no device state, CPU inference runs.
    tmpdir   = mktempdir()
    npe_path = joinpath(tmpdir, "npe_gpu_$(G).jld2")
    ProteinCoLoc.save_estimator(npe_path, res.estimator, res.θzt, zt, res.arch;
                                meta = (; grid = G, trained_use_gpu = true))
    on_disk = JLD2.load(npe_path)
    @test haskey(on_disk, "model_state")     # CPU-resident Flux.state snapshot
    @test !haskey(on_disk, "estimator")      # NOT the whole PosteriorEstimator object

    loaded = ProteinCoLoc.load_estimator(npe_path)
    draws  = ProteinCoLoc.posterior_for(loaded.estimator, Zstd[:, 1]; N = 8, use_gpu = false)
    @test size(draws, 1) == 7
    @test loaded.zt.mean == zt.mean

    # The CUDA-present branch is guarded so CPU-only CI exercises only the graceful fallback.
    if ProteinCoLoc.has_cuda_device()
        rg = ProteinCoLoc.train_npe(Zstd, Zstd, θ, θ; use_gpu = true, epochs = 2, tiny...)
        @test rg.d_in == d
        # a GPU-trained net still persists CPU-resident and reloads CPU-side.
        gpath = joinpath(tmpdir, "npe_gpu_dev_$(G).jld2")
        ProteinCoLoc.save_estimator(gpath, rg.estimator, rg.θzt, zt, rg.arch)
        lg = ProteinCoLoc.load_estimator(gpath)
        @test size(ProteinCoLoc.posterior_for(lg.estimator, Zstd[:, 1]; N = 8, use_gpu = false), 1) == 7
    else
        # Document the path taken: this machine has no CUDA, so use_gpu=true degraded to CPU
        # WITHOUT error — the graceful-fallback constraint (CLAUDE.md) held.
        @test ProteinCoLoc.has_cuda_device() == false
    end
end
