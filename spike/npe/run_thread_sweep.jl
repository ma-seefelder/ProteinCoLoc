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

# spike/npe/run_thread_sweep.jl --- Phase-4 CPU core-scaling sweep (NPE-03 / D-13; Pitfall 7).
#
# `Threads.nthreads()` is FIXED at process launch (`julia -t N`) and IMMUTABLE at runtime
# (RESEARCH Pitfall 7): a thread-count sweep therefore CANNOT be done by mutating the
# thread count inside one session. This harness SPAWNS a SEPARATE `julia --project=spike
# -t N` process for each N ∈ {1,2,4,8}; each subprocess runs the 04-05 `speedup_report`
# at its own fixed thread count and writes its NPE/ADVI wall-clock to a per-N JLD2. The
# driver then AGGREGATES the per-N files into a single thread-scaling table, identifying
# the parallel speedup / Amdahl saturation of the NPE clock.
#
# THE HEADLINE IS STATED AT A FIXED, REPORTED THREAD COUNT (D-13): the >100× claim is
# read at `BENCH_THREADS` (the pre-registered headline thread count, 04-01), with the full
# {1,2,4,8} sweep reported ALONGSIDE as characterization -- NOT at the best-performing
# thread count. Core-scaling is a reported characterization, not a pass/fail gate.
#
# HONEST CAVEAT (04-05 carry-forward): the realized ADVI baseline is fast (~0.5 s) and the
# NPE clock is the amortized to-posterior forward pass (bench_N draws, D-08 NPE-conservative
# clock: summary extraction is timed for the NPE but excluded for ADVI; WR-03).
# The sweep characterizes how the NPE clock parallelizes; the ADVI wall_clock is the frozen
# artifact reference (04-04), so the ratio's thread dependence is driven by the NPE side.
#
# CPU-ONLY (D-10): every spawned process is CPU-only (`use_gpu=false` inside speedup_report;
# CUDA is never imported). DECOUPLING (CLAUDE.md): spike-local; reaches src/ only
# transitively through the guarded benchmark.jl chain.

using JLD2            # per-N result + aggregation persistence (atomic .tmp+reopen+mv)
using Statistics      # median

# benchmark.jl brings speedup_report + the whole timing chain into scope. Guarded so the
# file is loadable standalone, as a spawned worker, AND inside runtests.jl.
isdefined(@__MODULE__, :speedup_report) || include(joinpath(@__DIR__, "benchmark.jl"))

# Pre-registered headline thread count (declared in test_npe.jl at 04-01; guarded fallback
# so a standalone include / spawned worker still gets the locked value).
if !isdefined(@__MODULE__, :BENCH_THREADS)
    const BENCH_THREADS = 1        # D-13: the fixed, reported headline thread count
end

const DEFAULT_THREAD_SWEEP = [1, 2, 4, 8]                       # D-13 sweep grid
const DEFAULT_SWEEP_OUTFILE = "thread_sweep.jld2"              # aggregation artifact name

"""
    _atomic_jldsave(path; kwargs...) -> String

Persist `kwargs` to `path` via the cache.jl atomic idiom: `jldsave` to a `.tmp`, REOPEN
to integrity-check a sentinel key, then `mv(...; force = true)`. A crash before the `mv`
leaves only a discardable `.tmp`, never a torn result file.
"""
function _atomic_jldsave(path; kwargs...)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "threads") "thread_sweep: integrity check failed: $tmp"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    thread_sweep_measure(outpath, holdir; master_seed, artifact_path, model_path,
                         bench_N, bench_seconds, expected_threads = nothing) -> String

Run `speedup_report` ONCE at the CURRENT process thread count and write the per-N result
to `outpath` (atomic). Records the fixed `Threads.nthreads()`, the median/full-N speedup,
the median NPE and ADVI wall-clocks, and whether the paired RMSE stayed within the
pre-registered tolerance. If `expected_threads` is given it is ASSERTED equal to
`Threads.nthreads()` -- the guard that this measurement really ran at the intended thread
count (a spawned `julia -t N`, never a mutated in-session count; Pitfall 7).
"""
function thread_sweep_measure(outpath, holdir; master_seed, artifact_path, model_path,
                              bench_N::Integer = 50, bench_seconds = 0.4,
                              expected_threads::Union{Nothing,Integer} = nothing)
    nt = Threads.nthreads()
    expected_threads === nothing || nt == expected_threads ||
        error("thread_sweep_measure: process has $nt threads, expected $expected_threads " *
              "(Pitfall 7 -- launch a separate `julia -t $expected_threads`).")
    sr = speedup_report(holdir; master_seed = master_seed, artifact_path = artifact_path,
                        model_path = model_path, bench_N = bench_N,
                        bench_seconds = bench_seconds)
    rmse_ratio_ok = sr.rmse.npe_rho_rmse <= 1.2 * sr.rmse.advi_rho_rmse    # D-09 tolerance
    _atomic_jldsave(outpath;
        threads             = nt,
        median_speedup      = sr.median_speedup,
        full_median_speedup = sr.full_median_speedup,
        t_npe_median        = median(sr.t_npe),
        t_advi_median       = median(sr.t_advi),
        rmse_ratio_ok       = rmse_ratio_ok,
        use_gpu             = sr.use_gpu)
    return outpath
end

"""
    aggregate_thread_sweep(files; outpath, bench_threads = BENCH_THREADS) -> NamedTuple

Aggregate per-N result files (each written by `thread_sweep_measure`) into ONE thread-
scaling table, sorted by thread count. Computes the NPE-clock PARALLEL SPEEDUP relative
to the smallest thread count (Amdahl view: `t_npe_median[1] / t_npe_median[i]`) and reads
the >100× HEADLINE at `bench_threads` (D-13 -- the fixed, reported thread count, NOT the
best). Persists the aggregation to `outpath` (atomic). Returns
`(threads, median_speedup, full_median_speedup, t_npe_median, t_advi_median,
parallel_speedup, headline_threads, headline_speedup, outpath)`.
"""
function aggregate_thread_sweep(files; outpath, bench_threads::Integer = BENCH_THREADS)
    isempty(files) && error("aggregate_thread_sweep: no per-N result files given.")
    rows = [JLD2.load(f) for f in files]
    perm = sortperm([Int(r["threads"]) for r in rows])
    rows = rows[perm]

    threads             = [Int(r["threads"])            for r in rows]
    median_speedup      = [Float64(r["median_speedup"]) for r in rows]
    full_median_speedup = [Float64(r["full_median_speedup"]) for r in rows]
    t_npe_median        = [Float64(r["t_npe_median"])   for r in rows]
    t_advi_median       = [Float64(r["t_advi_median"])  for r in rows]

    # Amdahl view: NPE-clock speedup vs the smallest thread count in the sweep.
    base_t = t_npe_median[1]
    parallel_speedup = base_t ./ t_npe_median

    hi = findfirst(==(bench_threads), threads)
    headline_speedup = hi === nothing ? NaN : median_speedup[hi]

    _atomic_jldsave(outpath;
        threads = threads, median_speedup = median_speedup,
        full_median_speedup = full_median_speedup,
        t_npe_median = t_npe_median, t_advi_median = t_advi_median,
        parallel_speedup = parallel_speedup,
        headline_threads = bench_threads, headline_speedup = headline_speedup,
        schema_version = 1)

    return (threads = threads, median_speedup = median_speedup,
            full_median_speedup = full_median_speedup,
            t_npe_median = t_npe_median, t_advi_median = t_advi_median,
            parallel_speedup = parallel_speedup,
            headline_threads = bench_threads, headline_speedup = headline_speedup,
            outpath = outpath)
end

"""
    run_thread_sweep(holdir; threads = DEFAULT_THREAD_SWEEP, master_seed,
                     artifact_path = DEFAULT_ADVI_ARTIFACT, model_path = DEFAULT_NPE_MODEL,
                     bench_N = 50, bench_seconds = 0.4, outdir = mktempdir(),
                     spawn = true) -> NamedTuple

The D-13 core-scaling driver. For each `n` in `threads`:

  - `spawn = true` (the REPORTED path): launch a SEPARATE `julia --project=spike -t n`
    subprocess running THIS file in `--worker` mode (Pitfall 7 -- the ONLY correct way to
    sweep thread counts; never mutate `nthreads()` in-session). Each writes its per-N JLD2.
  - `spawn = false` (the fast in-process path, for the test at a SINGLE thread count):
    run `thread_sweep_measure` directly in the current process, asserting its
    `Threads.nthreads()` matches the requested `n`.

Then `aggregate_thread_sweep` merges the per-N files into `outdir/thread_sweep.jld2` and
returns the aggregation, with the >100× headline read at `BENCH_THREADS`.
"""
function run_thread_sweep(holdir; threads = DEFAULT_THREAD_SWEEP, master_seed,
                          artifact_path = DEFAULT_ADVI_ARTIFACT,
                          model_path = DEFAULT_NPE_MODEL,
                          bench_N::Integer = 50, bench_seconds = 0.4,
                          outdir = mktempdir(), spawn::Bool = true)
    mkpath(outdir)
    perN_files = String[]
    for n in threads
        f = joinpath(outdir, "thread_$(n).jld2")
        if spawn
            # Pitfall 7: a genuinely separate process per thread count.
            julia_exe = joinpath(Sys.BINDIR, Base.julia_exename())
            proj      = abspath(joinpath(@__DIR__, ".."))       # the spike/ env
            script    = abspath(@__FILE__)
            run(`$julia_exe --project=$proj -t $n $script --worker $f $holdir $artifact_path $model_path $bench_N $bench_seconds $(string(master_seed))`)
        else
            thread_sweep_measure(f, holdir; master_seed = master_seed,
                                 artifact_path = artifact_path, model_path = model_path,
                                 bench_N = bench_N, bench_seconds = bench_seconds,
                                 expected_threads = n)
        end
        push!(perN_files, f)
    end
    return aggregate_thread_sweep(perN_files;
        outpath = joinpath(outdir, DEFAULT_SWEEP_OUTFILE), bench_threads = BENCH_THREADS)
end

# --- Spawned-worker entry point ------------------------------------------------------
# When this file is RUN as `julia -t N run_thread_sweep.jl --worker <out> <holdir>
# <artifact> <model> <bench_N> <bench_seconds> <master_seed>` it measures one per-N
# result at the launched (fixed) thread count. Guarded so an `include(...)` (test /
# driver) never triggers the worker branch.
if abspath(PROGRAM_FILE) == (@__FILE__) && !isempty(ARGS) && ARGS[1] == "--worker"
    _out, _hol, _art, _mod = ARGS[2], ARGS[3], ARGS[4], ARGS[5]
    _bn   = parse(Int, ARGS[6])
    _bs   = parse(Float64, ARGS[7])
    _seed = parse(UInt, ARGS[8])
    thread_sweep_measure(_out, _hol; master_seed = _seed, artifact_path = _art,
                         model_path = _mod, bench_N = _bn, bench_seconds = _bs,
                         expected_threads = Threads.nthreads())
    println("thread_sweep worker: wrote $_out at $(Threads.nthreads()) threads")
end
