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

# spike/npe/flow_width_ab.jl --- MEASUREMENT ONLY: does flow width (D_flow 7 vs 8)
# explain the Phase-4 SC3 speedup shortfall?
#
# BACKGROUND. `spike/test/test_npe.jl:74` pre-registers `SPEEDUP_GATE = 100.0`
# (NPE-03: median t_advi/t_npe > 100 at BENCH_THREADS). The gate at :230 has been
# measuring 84-93 across several serial runs -- reproducible, not noise. Because a
# thrown testset aborts every later `include`, that failure also prevents the Phase-5/9/13
# testsets from running at all, which is why it matters beyond its own claim.
#
# THE RECORDED HYPOTHESIS (STATE.md, Phase-11 wave-3 blocker): Phase-11 decision D-09
# appended `chromatic_eps` as an 8th theta column, so the spike NPE now trains an
# 8-marginal flow instead of 7, and the wider flow is a slower forward pass. STATE.md
# records this as "likely cause, NOT yet confirmed by measurement" and names the missing
# evidence exactly: "an A/B of the benchmark at D_flow 7 vs 8". This script IS that A/B.
#
# WHAT THIS SCRIPT DOES NOT DO. It changes no threshold, edits no test, and retires no
# gate. It writes one JLD2 artifact and prints a summary. Retiring or re-deriving
# SPEEDUP_GATE is a pre-registration decision the user owns.
#
# NO TRAINING IS REQUIRED, AND THAT IS A METHODOLOGICAL POINT, not a shortcut. The
# forward-pass cost of a flow is fixed by its ARCHITECTURE (layer shapes, coupling-stack
# depth), not by the VALUES its weights happen to hold. A freshly-initialised D=8 flow
# therefore costs exactly what a trained D=8 flow costs, so the width question is
# answerable for free. The `:real7` arm below (the actual persisted, trained gate model)
# is measured alongside the fresh arms precisely so that claim is checked rather than
# assumed: `:real7` and `:fresh7` must agree.
#
# THREE ARMS, all timed through the PRODUCTION `_time_npe_pair` (benchmark.jl) so this
# measures the same block the gate scores, not a re-implementation of it:
#   :real7  -- the persisted `trained_npe.jld2` exactly as `speedup_report` loads it
#              (q.d = 7, d_in = 128). This is what the failing gate measures TODAY.
#   :fresh7 -- `build_estimator(128, 7)`, untrained, same NPE_* consts.
#   :fresh8 -- `build_estimator(128, 8)`, untrained, IDENTICAL to :fresh7 except D.
# :fresh7 vs :fresh8 isolates flow width as the ONLY difference (same d_in = 128, same
# dstar/depth/width/coupling), which is what "does width explain it" requires.
#
# CONTENTION CONTROL (the machine is shared with other agents' Julia jobs):
#   * arms are INTERLEAVED, and their order ROTATES per repetition, so load drift lands
#     on all three arms roughly equally instead of biasing whichever ran first;
#   * every arm is repeated REPS times and reported as a median WITH spread (min/max,
#     quartiles), never as a point estimate;
#   * `@belapsed` (inside `_time_npe_pair`) reports the MINIMUM elapsed of many samples,
#     which already filters transient contention spikes;
#   * the concurrent `julia` process count is recorded before and after the run and
#     stored in the artifact, so the machine state is part of the record.
# If the arms' spreads OVERLAP, the honest reading is INCONCLUSIVE, not "the medians
# differ by x%" -- the printed summary says so explicitly.
#
# RUN (1 thread: the headline is only defined at BENCH_THREADS = 1, D-13):
#   julia --project=spike spike/npe/flow_width_ab.jl
# Optional env overrides: AB_REPS (default 9), AB_BENCH_SECONDS (default 0.4, the value
# SC3 itself passes), AB_BENCH_N (default 50, likewise SC3's).

using Statistics
using JLD2
using Random

# Brings the production timing block (`_time_npe_pair`), the artifact reader
# (`_load_artifact`), `load_npe`/`load_holdout`, `resimulate_holdout`, and
# `build_estimator` + the NPE_* architecture consts into scope.
isdefined(@__MODULE__, :_time_npe_pair) || include(joinpath(@__DIR__, "benchmark.jl"))

const AB_HOLDOUT_DIR  = joinpath(@__DIR__, "..", "baseline", "holdout")
const AB_ARTIFACT     = joinpath(@__DIR__, "..", "baseline", "advi_artifact.jld2")
const AB_MODEL        = joinpath(@__DIR__, "trained_npe.jld2")
const AB_OUT          = joinpath(@__DIR__, "flow_width_ab.jld2")

const AB_MASTER_SEED  = 0xC0FFEE   # == NPE_MASTER_SEED (test_npe.jl:78); holdout stream key
const AB_REPS         = parse(Int,     get(ENV, "AB_REPS",          "9"))
const AB_BENCH_N      = parse(Int,     get(ENV, "AB_BENCH_N",       "50"))
const AB_BENCH_SECS   = parse(Float64, get(ENV, "AB_BENCH_SECONDS", "0.4"))

"""
    now_utc() -> String

Integer wall-clock stamp for the artifact, without pulling in a `Dates` dependency
(the spike env's package set is an asserted, byte-pinned list; see runtests.jl gate (i)).
"""
now_utc() = string(round(Int, time()))

"""
    _julia_proc_count() -> Int

Concurrent `julia` process count (this process included), or -1 if unobtainable. Recorded
before and after the run so the artifact carries the machine state the numbers were taken
under -- a wall-clock benchmark taken under unrecorded contention is not interpretable.
"""
function _julia_proc_count()
    try
        out = read(`powershell -NoProfile -Command "(@(Get-Process julia -ErrorAction SilentlyContinue).Count)"`, String)
        return parse(Int, strip(out))
    catch
        return -1
    end
end

"""
    _quantiles(v) -> NamedTuple

Median plus spread (min, q25, q75, max) of a repetition vector. Spread is reported for
every arm because a single pair of numbers taken under contention settles nothing.
"""
function _quantiles(v)
    s = sort(collect(Float64, v))
    return (median = median(s), min = first(s), max = last(s),
            q25 = quantile(s, 0.25), q75 = quantile(s, 0.75))
end

"""
    run_flow_width_ab(; reps, bench_N, bench_seconds, outpath) -> NamedTuple

Interleaved A/B of the SC3 speedup benchmark at flow width D = 7 vs D = 8.

The holdout stacks are re-simulated ONCE up front (outside every timed block, mirroring
`speedup_report`, whose clock likewise starts from an already-materialised raw stack), so
repetitions re-time the same in-memory inputs and the arms cannot differ by re-simulation
luck. Per repetition and per arm the per-pair speedup `t_advi / t_npe` is formed and
reduced by `median` -- the exact statistic `speedup_report` reports and SC3 gates on.
"""
function run_flow_width_ab(; reps::Integer = AB_REPS, bench_N::Integer = AB_BENCH_N,
                           bench_seconds::Real = AB_BENCH_SECS, outpath = AB_OUT)
    procs_before = _julia_proc_count()
    @info "flow-width A/B starting" threads=Threads.nthreads() reps bench_N bench_seconds julia_procs=procs_before

    m   = load_npe(AB_MODEL)
    art = _load_artifact(AB_ARTIFACT)
    ho  = load_holdout(AB_HOLDOUT_DIR)

    # The gate's own model must be the 7-marginal one for this A/B to be interpretable;
    # record what it actually is rather than trusting the hypothesis.
    real_D = m.estimator.q.d
    @info "persisted gate model" path=AB_MODEL flow_d=real_D d_in=m.d_in variant=m.variant

    gi_to_col = Dict(Int(ho.global_index[j]) => j for j in 1:length(ho.global_index))
    pairs   = art["pairs"]
    npairs  = length(pairs)
    t_advi  = Float64.(art["wall_clock"])

    # Re-simulate every holdout stack ONCE (setup, never timed).
    mci_s = Vector{Any}(undef, npairs)
    mci_c = Vector{Any}(undef, npairs)
    for p in 1:npairs
        js = gi_to_col[Int(pairs[p][1])]
        jc = gi_to_col[Int(pairs[p][2])]
        mci_s[p] = resimulate_holdout(AB_HOLDOUT_DIR, js; master_seed = AB_MASTER_SEED).mci_sample
        mci_c[p] = resimulate_holdout(AB_HOLDOUT_DIR, jc; master_seed = AB_MASTER_SEED).mci_sample
    end
    @info "holdout re-simulated (outside all timed blocks)" npairs

    # The three arms. :fresh7/:fresh8 differ ONLY in D -- same d_in, same NPE_* consts.
    Random.seed!(UInt32(AB_MASTER_SEED & 0xFFFFFFFF))
    est_fresh7 = build_estimator(m.d_in, 7)
    est_fresh8 = build_estimator(m.d_in, 8)
    arms = [(:real7, m.estimator), (:fresh7, est_fresh7), (:fresh8, est_fresh8)]
    names = first.(arms)

    speedup = Dict(n => Float64[] for n in names)   # per-rep median speedup
    tnpe    = Dict(n => Float64[] for n in names)   # per-rep median NPE latency (s)

    for r in 1:reps
        order = circshift(1:length(arms), r)        # rotate arm order to cancel load drift
        for i in order
            nm, est = arms[i]
            tp = Vector{Float64}(undef, npairs)
            for p in 1:npairs
                tp[p] = _time_npe_pair(est, mci_s[p], mci_c[p], m.zt, m.variant, bench_N;
                                       bench_seconds = bench_seconds)
            end
            push!(speedup[nm], median(t_advi ./ tp))
            push!(tnpe[nm],    median(tp))
        end
        @info "rep $r/$reps" real7=round(last(speedup[:real7]), digits=2) fresh7=round(last(speedup[:fresh7]), digits=2) fresh8=round(last(speedup[:fresh8]), digits=2)
    end

    procs_after = _julia_proc_count()
    stats = Dict(string(n) => _quantiles(speedup[n]) for n in names)

    # --- Report -----------------------------------------------------------------
    println("\n", "="^74)
    println("FLOW-WIDTH A/B  --  median speedup t_advi/t_npe  (gate: > 100)")
    println("="^74)
    println("threads = $(Threads.nthreads())   reps = $reps   bench_N = $bench_N   bench_seconds = $bench_seconds")
    println("julia processes: $procs_before at start, $procs_after at end (this one included)")
    println("persisted gate model $(basename(AB_MODEL)): flow d = $real_D, d_in = $(m.d_in)")
    println("-"^74)
    for nm in names
        q = stats[string(nm)]
        println(rpad(string(nm), 9),
                "median ", rpad(round(q.median, digits = 2), 9),
                "[min ", round(q.min, digits = 2), ", max ", round(q.max, digits = 2), "]  ",
                "IQR [", round(q.q25, digits = 2), ", ", round(q.q75, digits = 2), "]  ",
                "t_npe ", round(median(tnpe[nm]) * 1e3, digits = 2), " ms")
    end
    println("-"^74)

    q7, q8 = stats["fresh7"], stats["fresh8"]
    overlap = q7.min <= q8.max && q8.min <= q7.max      # do the observed ranges overlap?
    ratio   = q7.median / q8.median
    println("fresh7 vs fresh8: median ratio = ", round(ratio, digits = 3),
            "  (", round((ratio - 1) * 100, digits = 1), "% cost of the 8th marginal)")
    println("observed ranges ", overlap ? "OVERLAP -> width effect NOT resolved by this run" :
                                          "are DISJOINT -> width effect is resolved")
    println("both arms clear 100? fresh7 ", q7.median > 100 ? "YES" : "NO",
            " / fresh8 ", q8.median > 100 ? "YES" : "NO")
    println("="^74, "\n")

    jldsave(outpath;
            arms = string.(names), speedup = Dict(string(k) => v for (k, v) in speedup),
            t_npe = Dict(string(k) => v for (k, v) in tnpe), stats = stats,
            t_advi = t_advi, n_pairs = npairs, reps = reps, bench_N = bench_N,
            bench_seconds = bench_seconds, threads = Threads.nthreads(),
            gate_model_flow_d = real_D, gate_model_d_in = m.d_in,
            julia_procs_before = procs_before, julia_procs_after = procs_after,
            speedup_gate = 100.0, generated = string(now_utc()))
    @info "artifact written" outpath
    return (stats = stats, speedup = speedup, t_npe = tnpe, overlap = overlap,
            real_D = real_D, procs_before = procs_before, procs_after = procs_after)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_flow_width_ab()
end
