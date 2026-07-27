# ProteinCoLoc: Bayesian colocalization analysis of multi-channel fluorescence microscopy images
# Copyright (C) 2024  Manuel Seefelder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# spike/validation/run_p11_paired.jl --- Phase-11 redesign FEASIBILITY test (research lane only).
#
# TWO QUESTIONS, BOTH ANSWERED WITH THE SAME NOISE-FLOOR INSTRUMENT BUILT IN
# `run_p11_bottleneck.jl`, so every number here is directly comparable to the ones already
# reported for the CURRENT summary (3 px effect 0.3966 against a floor of 0.9765, i.e. 0.41).
#
# Q1 -- DOES A FINER GRID HELP?  Not assumed: MEASURED. A finer grid raises the signal (3 px is a
#      larger fraction of a smaller patch) but ALSO raises the noise floor, because fewer pixels
#      per patch means a noisier correlation estimate. The decisive quantity is the RATIO, and it
#      is entirely possible it is flat in grid size -- which would rule out the whole "finer grid"
#      family. Swept at 8x8, 16x16, 32x32.
#
# Q2 -- DOES A SELF-REFERENTIAL (PAIRED) FEATURE HELP?  The current summary evaluates the patch-
#      correlation function at ONE fixed alignment. The registration information does not live in
#      that value -- it lives in the SHAPE of the correlation function around that point, which the
#      current summary discards entirely. That alone explains the null result, and it says the fix
#      is a different QUANTITY, not more resolution.
#
#      So: displace channel 2 against channel 1 by a set of KNOWN probe offsets, recompute the
#      patch correlations at each, and read features off the resulting curve. This is classical
#      cross-correlation registration with sub-pixel peak interpolation -- standard practice, not
#      an invented statistic. The mechanism that should make it work is the SAME one that makes the
#      D-06 probe sensitive: it compares an observation against ITSELF, so the scene and draw noise
#      are largely common across offsets and cancel in the differences, whereas the measured floor
#      of 0.9765 came from comparing INDEPENDENT draws.
#
#      THIS IS A HYPOTHESIS UNDER TEST, NOT A CONCLUSION. The verdict is pre-committed below.
#
# PRE-COMMITTED VERDICT (decided BEFORE looking, so the result cannot be rounded toward the hope):
#   CLEAR SUCCESS   : peak-location SNR > 1 AND ridge ratio < 0.8 across the ladder.
#   CLEAR FAILURE   : SNR still < 1 AND ridge still ~1.0. Then 3 px is simply not observable from
#                     a single draw at ANY summary built on patch correlations -- a legitimate and
#                     publishable answer, to be reported WITHOUT softening.
#   AMBIGUOUS       : anything else. Reported as ambiguous, NOT rounded toward the hypothesis.
#
# SCOPE. Research-lane feasibility only. Changes NO shipped summary contract, touches no `src/`
# file, retrains and reships nothing, alters no constant in `p11_consts.jl`, and spends no
# iteration allowance. The probe warp reuses the simulator's OWN `warp` + `BSpline(Linear())` path
# (`spike/simulator/forward.jl:216-218`) rather than a new interpolator.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_paired.jl

using JLD2
using Dates
using Statistics
using ImageTransformations
using CoordinateTransformations
using Interpolations

const P11_PROBE_LOAD_ONLY = true
include(joinpath(@__DIR__, "run_p11_probe.jl"))

const PAIRED_REPORT_PATH = joinpath(@__DIR__, "p11_paired_report.jld2")

const PAIRED_R          = 16                       # replicates per theta base (80 pairs total)
const PAIRED_GRIDS      = (8, 16, 32)              # Q1 sweep
const PAIRED_TEST_SHIFT = 3.0                      # = LAMBDA_MAX, the effect that must be seen
const PROBE_OFFSETS     = collect(-4.0:0.25:4.0)   # spans +/- LAMBDA_MAX with sub-pixel steps

# Independent key family, disjoint from the probe's (+0), imsize (+1000) and the bottleneck
# script's (+2000) blocks.
_paired_indep_key(r::Integer) = P11_PROBE_COUNTER * 1000 + 4000 + Int(r)

"""
    _rehydrate_theta(values, fields) -> NamedTuple

The probe artifact stores its frozen theta bases as plain Float64 rows plus a `theta_fields` name
vector. Rebuild NamedTuples from the artifact's OWN field vector so this script stays row-aligned
with the probe if the theta arity ever changes.
"""
function _rehydrate_theta(values::AbstractVector{<:Real}, fields::AbstractVector)
    @assert length(values) == length(fields)
    return NamedTuple{Tuple(Symbol.(fields))}(Tuple(Float64.(values)))
end

"""
    _raw_pair(key, theta; imsize) -> Vector{Matrix{Float64}}

The OBSERVED image pair, i.e. `simulate_pair` output before any summarisation. The probe features
below operate on this, because a registration probe must act on the observation -- re-simulating
at a different shift would change the noise realisation and destroy the self-reference that is the
entire point.
"""
_raw_pair(key::Integer, theta; imsize::Tuple{Int,Int}) =
    simulate_pair(p11_rng(key), theta; imsize = imsize)

"""
    _grid_corr(ch1, ch2, n) -> Vector{Float64}

The `n x n` patch-correlation vector, with missing patches (below the >=15-px floor) coalesced to
0.0 so the vector length is fixed at `n^2` and two draws stay index-aligned. Uses the SAME
`patch` + `correlation` pair `spike/contract.jl:85-90` calls, with the grid size as a parameter
instead of the hard-coded 8.
"""
function _grid_corr(ch1::Matrix{Float64}, ch2::Matrix{Float64}, n::Integer)
    xp, yp = patch.([ch1, ch2], n)
    c = correlation(xp, yp; method = :pearson)
    return [ismissing(v) || !isfinite(v) ? 0.0 : Float64(v) for v in vec(c)]
end

"L2 norm between two equal-length summary vectors (the probe's metric, generalised to any grid)."
_vec_norm(a::Vector{Float64}, b::Vector{Float64}) = sqrt(sum(abs2, a .- b))

"""
    _shift_channel(ch, du, dv) -> Matrix{Float64}

Displace an ALREADY-OBSERVED channel by a known sub-pixel probe offset, using the simulator's own
warp path: backward-mode `Translation(dv, du)` (slot 1 = rows = dy, slot 2 = columns = dx) with
`BSpline(Linear())` and `fillvalue = BG_FLOOR`, exactly as `forward.jl:216-218` does. Reusing that
path rather than writing an interpolator keeps the probe's resampling identical in kind to the
one the forward model already applies.
"""
function _shift_channel(ch::Matrix{Float64}, du::Real, dv::Real)
    A = Translation(Float64(dv), Float64(du))
    return Matrix{Float64}(collect(warp(ch, A, axes(ch);
                                        method = BSpline(Linear()), fillvalue = BG_FLOOR)))
end

"""
    _probe_curve(ch1, ch2, offsets; n) -> (cx, cy)

The self-referential feature. For each probe offset along each axis, displace channel 2 and record
the MEAN `n x n` patch correlation against channel 1. Two 1-D sweeps (a cross), not a full 2-D
grid: the full grid costs `length(offsets)^2` warps per image for information the axis-separable
sweeps already carry, and the cost is the binding constraint on this test.
"""
function _probe_curve(ch1::Matrix{Float64}, ch2::Matrix{Float64}, offsets::Vector{Float64};
                      n::Integer = 8)
    cx = similar(offsets)
    cy = similar(offsets)
    for (i, o) in enumerate(offsets)
        cx[i] = mean(_grid_corr(ch1, _shift_channel(ch2, o, 0.0), n))
        cy[i] = mean(_grid_corr(ch1, _shift_channel(ch2, 0.0, o), n))
    end
    return cx, cy
end

"""
    _subpixel_peak(offsets, c) -> (location, peak_value, curvature)

Sub-pixel peak by parabolic interpolation on the three points around the discrete argmax -- the
standard refinement for correlation-based registration. `curvature` is the quadratic coefficient
(negative at a maximum); its magnitude is the peak sharpness. Falls back to the discrete argmax
when the peak sits on a boundary or the parabola is degenerate.
"""
function _subpixel_peak(offsets::Vector{Float64}, c::Vector{Float64})
    k = argmax(c)
    (k == 1 || k == length(c)) && return (offsets[k], c[k], 0.0)
    y0, y1, y2 = c[k-1], c[k], c[k+1]
    denom = (y0 - 2y1 + y2)
    abs(denom) < 1e-12 && return (offsets[k], c[k], 0.0)
    h = offsets[k+1] - offsets[k]
    delta = 0.5 * (y0 - y2) / denom
    abs(delta) > 1.0 && return (offsets[k], c[k], denom / (h^2))
    return (offsets[k] + delta * h, y1 - 0.25 * (y0 - y2) * delta, denom / (h^2))
end

# ============================================================================================
# Q1 -- grid-resolution sweep of the EXISTING summary
# ============================================================================================
function grid_sweep(thetas::Vector, imsize_of::Function)
    println("\n" * "="^78)
    println("Q1  GRID-RESOLUTION SWEEP OF THE CURRENT SUMMARY")
    println("Does a finer grid improve the signal-to-noise RATIO, or does it raise both equally?")
    println("="^78)

    T, R = length(thetas), PAIRED_R
    out = Dict{Int,Any}()

    for n in PAIRED_GRIDS
        sig = zeros(Float64, T, R)   # paired: same key, 0 vs PAIRED_TEST_SHIFT
        flr = zeros(Float64, T, R)   # independent keys, both at zero shift
        units = vec([(t, r) for t in 1:T, r in 1:R])
        t0 = time()
        Threads.@threads for u in units
            t, r = u
            sz = imsize_of(r)
            th0 = _set_eps(_set_shift(thetas[t], 0.0), 0.0)
            a  = _raw_pair(_pair_key(r), th0; imsize = sz)
            s0 = _grid_corr(a[1], a[2], n)
            b  = _raw_pair(_pair_key(r), _set_shift(th0, PAIRED_TEST_SHIFT); imsize = sz)
            sig[t, r] = _vec_norm(_grid_corr(b[1], b[2], n), s0)
            c  = _raw_pair(_paired_indep_key(r), th0; imsize = sz)
            flr[t, r] = _vec_norm(_grid_corr(c[1], c[2], n), s0)
        end
        ms, mf = mean(sig), mean(flr)
        out[n] = (signal = ms, floor = mf, snr = ms / mf,
                  signal_sd = std(sig), floor_sd = std(flr))
        println("  grid $(n)x$(n):  signal = $(round(ms; digits = 4))   " *
                "floor = $(round(mf; digits = 4))   SNR = $(round(ms/mf; digits = 4))   " *
                "($(round((time()-t0)/60; digits = 2)) min)")
    end
    return out
end

# ============================================================================================
# Q2 -- the self-referential paired feature
# ============================================================================================
function paired_feature_test(thetas::Vector, imsize_of::Function)
    println("\n" * "="^78)
    println("Q2  SELF-REFERENTIAL PROBE-SHIFT FEATURE")
    println("Peak location of the correlation-vs-probe-offset curve, in PIXELS.")
    println("="^78)

    T, R = length(thetas), PAIRED_R
    px0, py0 = zeros(Float64, T, R), zeros(Float64, T, R)   # true shift 0
    pxs, pys = zeros(Float64, T, R), zeros(Float64, T, R)   # true shift PAIRED_TEST_SHIFT
    pxi, pyi = zeros(Float64, T, R), zeros(Float64, T, R)   # independent draw, true shift 0
    sharp0   = zeros(Float64, T, R)

    units = vec([(t, r) for t in 1:T, r in 1:R])
    t0 = time()
    done = Threads.Atomic{Int}(0)
    Threads.@threads for u in units
        t, r = u
        sz  = imsize_of(r)
        th0 = _set_eps(_set_shift(thetas[t], 0.0), 0.0)

        a = _raw_pair(_pair_key(r), th0; imsize = sz)
        cx, cy = _probe_curve(a[1], a[2], PROBE_OFFSETS)
        lx, _, kx = _subpixel_peak(PROBE_OFFSETS, cx)
        ly, _, _  = _subpixel_peak(PROBE_OFFSETS, cy)
        px0[t, r], py0[t, r], sharp0[t, r] = lx, ly, abs(kx)

        # The SAME key, only the geometry differs -- the signal arm.
        b = _raw_pair(_pair_key(r), _set_shift(th0, PAIRED_TEST_SHIFT); imsize = sz)
        cxs, cys = _probe_curve(b[1], b[2], PROBE_OFFSETS)
        pxs[t, r], _, _ = _subpixel_peak(PROBE_OFFSETS, cxs)
        pys[t, r], _, _ = _subpixel_peak(PROBE_OFFSETS, cys)

        # An INDEPENDENT draw at the same theta and zero shift -- the noise arm. The scatter of
        # (this - px0) is the per-draw uncertainty of the feature, in pixels.
        c = _raw_pair(_paired_indep_key(r), th0; imsize = sz)
        cxi, cyi = _probe_curve(c[1], c[2], PROBE_OFFSETS)
        pxi[t, r], _, _ = _subpixel_peak(PROBE_OFFSETS, cxi)
        pyi[t, r], _, _ = _subpixel_peak(PROBE_OFFSETS, cyi)

        n = Threads.atomic_add!(done, 1) + 1
        n % 8 == 0 && println("  [Q2] $n/$(length(units)) units " *
                              "($(round((time()-t0)/60; digits = 2)) min)")
    end

    # SIGNAL: how far the estimated peak moves when the true shift changes by PAIRED_TEST_SHIFT
    # along the diagonal, i.e. PAIRED_TEST_SHIFT/sqrt(2) per axis (the `_diag_shift` convention).
    dx_move = vec(pxs .- px0)
    dy_move = vec(pys .- py0)
    # FLOOR: scatter of the SAME estimate between two independent draws of the same scene.
    dx_noise = vec(pxi .- px0)
    dy_noise = vec(pyi .- py0)

    expected_per_axis = PAIRED_TEST_SHIFT / sqrt(2.0)
    sig_x, sig_y = mean(abs.(dx_move)), mean(abs.(dy_move))
    flr_x, flr_y = std(dx_noise),       std(dy_noise)

    println()
    println("  expected per-axis movement for a $(PAIRED_TEST_SHIFT) px diagonal shift = " *
            "$(round(expected_per_axis; digits = 4)) px")
    println("  x-axis: mean |peak move| = $(round(sig_x; digits = 4)) px   " *
            "per-draw noise sd = $(round(flr_x; digits = 4)) px   " *
            "SNR = $(round(sig_x/flr_x; digits = 4))")
    println("  y-axis: mean |peak move| = $(round(sig_y; digits = 4)) px   " *
            "per-draw noise sd = $(round(flr_y; digits = 4)) px   " *
            "SNR = $(round(sig_y/flr_y; digits = 4))")
    println("  peak-location bias at zero shift: x = $(round(mean(vec(px0)); digits = 4)) px, " *
            "y = $(round(mean(vec(py0)); digits = 4)) px  (should be near 0)")

    return (px0 = px0, py0 = py0, pxs = pxs, pys = pys, pxi = pxi, pyi = pyi,
            sharp0 = sharp0,
            signal_x = sig_x, signal_y = sig_y, floor_x = flr_x, floor_y = flr_y,
            snr_x = sig_x / flr_x, snr_y = sig_y / flr_y,
            expected_per_axis = expected_per_axis)
end

function main()
    println("="^78)
    println("PHASE-11 REDESIGN FEASIBILITY  (research lane; gates nothing; changes no contract)")
    println("="^78)
    t_start = time()

    report = JLD2.load(P11_PROBE_REPORT_PATH)
    tf     = report["theta_fields"]
    thetas = [_rehydrate_theta(row, tf) for row in report["theta_bases"]]
    imsize_of = r -> _draw_imsize(p11_rng(_imsize_key(r)))
    println("theta bases: $(length(thetas))   replicates: $PAIRED_R   " *
            "probe offsets: $(length(PROBE_OFFSETS)) over " *
            "[$(first(PROBE_OFFSETS)), $(last(PROBE_OFFSETS))] px")

    q1 = grid_sweep(thetas, imsize_of)
    q2 = paired_feature_test(thetas, imsize_of)

    println("\n" * "="^78)
    println("VERDICT INPUTS")
    println("  current 8x8 summary SNR (from run_p11_bottleneck.jl) = 0.3966 / 0.9765 = 0.4061")
    for n in PAIRED_GRIDS
        println("  grid $(n)x$(n) SNR = $(round(q1[n].snr; digits = 4))")
    end
    println("  paired probe-shift feature SNR: x = $(round(q2.snr_x; digits = 4)), " *
            "y = $(round(q2.snr_y; digits = 4))")
    println("="^78)

    elapsed = (time() - t_start) / 60
    tmp = PAIRED_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version    = 1,
        grid_sizes        = collect(PAIRED_GRIDS),
        grid_signal       = [q1[n].signal for n in PAIRED_GRIDS],
        grid_floor        = [q1[n].floor  for n in PAIRED_GRIDS],
        grid_snr          = [q1[n].snr    for n in PAIRED_GRIDS],
        probe_offsets     = PROBE_OFFSETS,
        test_shift        = PAIRED_TEST_SHIFT,
        expected_per_axis = q2.expected_per_axis,
        peak_x_zero       = q2.px0,  peak_y_zero  = q2.py0,
        peak_x_shift      = q2.pxs,  peak_y_shift = q2.pys,
        peak_x_indep      = q2.pxi,  peak_y_indep = q2.pyi,
        peak_sharpness    = q2.sharp0,
        paired_signal_x   = q2.signal_x, paired_signal_y = q2.signal_y,
        paired_floor_x    = q2.floor_x,  paired_floor_y  = q2.floor_y,
        paired_snr_x      = q2.snr_x,    paired_snr_y    = q2.snr_y,
        baseline_snr_8x8  = 0.3966 / 0.9765,
        theta_base_n      = length(thetas),
        replicates        = PAIRED_R,
        elapsed_min       = elapsed,
        generated         = string(Dates.now(Dates.UTC)) * "Z",
        caption           = "Phase-11 redesign feasibility: grid-resolution SNR sweep of the " *
                            "current summary, and the self-referential probe-shift peak-location " *
                            "feature, both against the independent-draw noise floor. Research " *
                            "lane; gates nothing; changes no shipped contract.",
    )
    let d = JLD2.load(tmp); @assert haskey(d, "paired_snr_x") && haskey(d, "grid_snr"); end
    mv(tmp, PAIRED_REPORT_PATH; force = true)
    println("wrote ", PAIRED_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    return nothing
end

isdefined(@__MODULE__, :P11_PAIRED_LOAD_ONLY) || main()
