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

# spike/validation/run_p11_peakdecomp.jl --- Phase-11: WHY is the paired peak imprecise?
#
# The paired probe-shift feature is UNBIASED and correctly scaled (peak moves 2.141 px against
# 2.121 px expected, bias ~0.1 px at zero shift). The entire negative verdict therefore rests on a
# single number: the per-draw NOISE of the peak location, measured at 2.0-2.8 px.
#
# Peak-location precision factorises, to first order, as
#
#     sigma_peak  ~=  (peak width)  /  (signal-to-noise of the individual curve points)
#
# and those two are fixed by DIFFERENT means, so reporting only the combined 2.8 px hides which
# problem we actually have:
#   - BROAD PEAK, clean points  -> the curve inherits the coarse patch structure. Patch-averaged
#     correlation is an extremely smooth function of a sub-pixel offset, so the peak is wide and
#     localises poorly no matter how clean it is. FIXABLE by computing the curve at finer
#     granularity, where the width is set by the PSF (a few px) rather than by the patch size.
#   - SHARP PEAK, noisy points  -> the design is limited by photon/draw noise. No granularity
#     change rescues it, and the candidate is dead.
#
# THIS IS NOT THE 32x32 TEST ALREADY RULED OUT. That measured the STATIC summary at finer
# resolution and found SNR flat (0.403 / 0.338 / 0.351). This measures the granularity at which the
# PAIRED CURVE is computed -- an untested combination, and the one the decomposition points at.
#
# It also raises n for the x/y anisotropy question. At n = 80 the observed log-ratio
# ln(2.835/2.065) = 0.317 sat at ~2.8 standard errors (rel. SE of a sample sd ~= 1/sqrt(2(n-1))
# = 8.0 %; for a log ratio of two independent sds, sqrt(2) x that = 11.3 %). That is the zone where
# people fool themselves, so n is raised to 200 units, where the SE falls to ~7.1 % and the same
# effect would sit at ~4.5 SE.
#
# Research lane. Changes no contract, touches no `src/`, alters no constant, gates nothing.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_peakdecomp.jl

using JLD2
using Dates
using Statistics

const P11_PAIRED_LOAD_ONLY = true
include(joinpath(@__DIR__, "run_p11_paired.jl"))

const PEAK_REPORT_PATH = joinpath(@__DIR__, "p11_peakdecomp_report.jld2")

const PD_R          = 40                        # 5 thetas x 40 = 200 units (was 80)
const PD_OFFSETS    = collect(-4.0:0.5:4.0)
const PD_GRIDS      = (8, 16, 32)               # granularity of the CURVE, not of the static summary
const PD_GRID_R     = 12                        # 5 x 12 = 60 units for the granularity sweep
const PD_TEST_SHIFT = 3.0

_pd_indep_key(r::Integer) = P11_PROBE_COUNTER * 1000 + 6000 + Int(r)

"""
    _fwhm(offsets, c) -> Float64

Full width at half maximum of the correlation peak, measured on the peak PROMINENCE
(peak minus the curve's baseline) so a curve sitting on a high pedestal is not reported as
artificially narrow. Linear interpolation to the half-prominence crossings on each side; returns
the span of the offset grid when the curve never falls to half prominence inside it (i.e. the peak
is at least that wide).
"""
function _fwhm(offsets::Vector{Float64}, c::Vector{Float64})
    k    = argmax(c)
    base = minimum(c)
    half = base + 0.5 * (c[k] - base)
    (c[k] - base) <= 0 && return offsets[end] - offsets[1]

    lo = offsets[1]
    for i in k:-1:2
        if c[i-1] <= half
            t = (half - c[i-1]) / (c[i] - c[i-1])
            lo = offsets[i-1] + t * (offsets[i] - offsets[i-1]); break
        end
    end
    hi = offsets[end]
    for i in k:(length(c)-1)
        if c[i+1] <= half
            t = (half - c[i+1]) / (c[i] - c[i+1])
            hi = offsets[i+1] - t * (offsets[i+1] - offsets[i]); break
        end
    end
    return hi - lo
end

"The probe curve along the x axis only, at an arbitrary patch grid."
function _curve_x(ch1, ch2, offsets::Vector{Float64}, n::Integer)
    c = similar(offsets)
    for (i, o) in enumerate(offsets)
        c[i] = mean(_grid_corr(ch1, _shift_channel(ch2, o, 0.0), n))
    end
    return c
end

function main()
    println("="^78)
    println("PHASE-11 PEAK-NOISE DECOMPOSITION  (is the peak BROAD, or are the POINTS noisy?)")
    println("="^78)
    t0 = time()

    report = JLD2.load(P11_PROBE_REPORT_PATH)
    tf     = report["theta_fields"]
    thetas = [_rehydrate_theta(row, tf) for row in report["theta_bases"]]
    imsize_of = r -> _draw_imsize(p11_rng(_imsize_key(r)))
    T = length(thetas)

    # ---------------------------------------------------------------------------------------
    # PART A -- decomposition + anisotropy at the current 8x8 granularity, at raised n
    # ---------------------------------------------------------------------------------------
    println("\nPART A  n = $(T * PD_R) units, grid 8x8")
    px0 = zeros(Float64, T, PD_R); py0 = zeros(Float64, T, PD_R)
    pxi = zeros(Float64, T, PD_R); pyi = zeros(Float64, T, PD_R)
    fwx = zeros(Float64, T, PD_R); fwy = zeros(Float64, T, PD_R)
    ptn = zeros(Float64, T, PD_R)   # per-point curve noise
    prm = zeros(Float64, T, PD_R)   # peak prominence

    units = vec([(t, r) for t in 1:T, r in 1:PD_R])
    done  = Threads.Atomic{Int}(0)
    Threads.@threads for u in units
        t, r = u
        sz   = imsize_of(r)
        th0  = _set_eps(_set_shift(thetas[t], 0.0), 0.0)

        a  = _raw_pair(_pair_key(r), th0; imsize = sz)
        cx, cy = _probe_curve(a[1], a[2], PD_OFFSETS)
        px0[t, r], _, _ = _subpixel_peak(PD_OFFSETS, cx)
        py0[t, r], _, _ = _subpixel_peak(PD_OFFSETS, cy)
        fwx[t, r] = _fwhm(PD_OFFSETS, cx)
        fwy[t, r] = _fwhm(PD_OFFSETS, cy)
        prm[t, r] = maximum(cx) - minimum(cx)

        b  = _raw_pair(_pd_indep_key(r), th0; imsize = sz)
        cxi, cyi = _probe_curve(b[1], b[2], PD_OFFSETS)
        pxi[t, r], _, _ = _subpixel_peak(PD_OFFSETS, cxi)
        pyi[t, r], _, _ = _subpixel_peak(PD_OFFSETS, cyi)
        # PER-POINT curve noise: the scatter between two independent realisations of the SAME
        # curve, offset by offset. This is the denominator of the factorisation.
        ptn[t, r] = std(cxi .- cx)

        n = Threads.atomic_add!(done, 1) + 1
        n % 40 == 0 && println("  [A] $n/$(length(units)) ($(round((time()-t0)/60; digits=2)) min)")
    end

    noise_x = std(vec(pxi .- px0))
    noise_y = std(vec(pyi .- py0))
    fwhm_x  = mean(vec(fwx));  fwhm_y = mean(vec(fwy))
    pt_noise = mean(vec(ptn)); prom = mean(vec(prm))

    # SE of a log ratio of two independent sample sds, n units each.
    nA = length(units)
    se_log = sqrt(2.0) / sqrt(2.0 * (nA - 1))
    obs_log = log(noise_x / noise_y)

    println("\n  peak-location per-draw noise:  x = $(round(noise_x; digits = 4)) px   " *
            "y = $(round(noise_y; digits = 4)) px")
    println("  ANISOTROPY  log-ratio = $(round(obs_log; digits = 4))   " *
            "SE = $(round(se_log; digits = 4))   " *
            "=> $(round(abs(obs_log)/se_log; digits = 2)) SE  (n = $nA)")
    println("  curve FWHM:  x = $(round(fwhm_x; digits = 4)) px   y = $(round(fwhm_y; digits = 4)) px")
    println("  peak prominence = $(round(prom; digits = 6))   " *
            "per-point curve noise = $(round(pt_noise; digits = 6))")
    println("  curve SNR (prominence / point noise) = $(round(prom / pt_noise; digits = 4))")
    println("  FACTORISATION CHECK  FWHM / curve-SNR = " *
            "$(round(fwhm_x / (prom / pt_noise); digits = 4)) px  " *
            "vs measured $(round(noise_x; digits = 4)) px")

    # ---------------------------------------------------------------------------------------
    # PART B -- does computing the CURVE at finer granularity sharpen the peak?
    # ---------------------------------------------------------------------------------------
    println("\nPART B  curve granularity sweep, x axis, n = $(T * PD_GRID_R) units")
    gres = Dict{Int,Any}()
    for n in PD_GRIDS
        p0 = zeros(Float64, T, PD_GRID_R); pi_ = zeros(Float64, T, PD_GRID_R)
        fw = zeros(Float64, T, PD_GRID_R); pn = zeros(Float64, T, PD_GRID_R)
        pr = zeros(Float64, T, PD_GRID_R)
        gu = vec([(t, r) for t in 1:T, r in 1:PD_GRID_R])
        tg = time()
        Threads.@threads for u in gu
            t, r = u
            sz  = imsize_of(r)
            th0 = _set_eps(_set_shift(thetas[t], 0.0), 0.0)
            a   = _raw_pair(_pair_key(r), th0; imsize = sz)
            ca  = _curve_x(a[1], a[2], PD_OFFSETS, n)
            p0[t, r], _, _ = _subpixel_peak(PD_OFFSETS, ca)
            fw[t, r] = _fwhm(PD_OFFSETS, ca)
            pr[t, r] = maximum(ca) - minimum(ca)
            b   = _raw_pair(_pd_indep_key(r), th0; imsize = sz)
            cb  = _curve_x(b[1], b[2], PD_OFFSETS, n)
            pi_[t, r], _, _ = _subpixel_peak(PD_OFFSETS, cb)
            pn[t, r] = std(cb .- ca)
        end
        nz = std(vec(pi_ .- p0)); fm = mean(vec(fw))
        pnn = mean(vec(pn)); prr = mean(vec(pr))
        gres[n] = (peak_noise = nz, fwhm = fm, point_noise = pnn, prominence = prr,
                   snr = prr / pnn)
        println("  grid $(n)x$(n):  peak noise = $(round(nz; digits = 4)) px   " *
                "FWHM = $(round(fm; digits = 4)) px   " *
                "curve SNR = $(round(prr/pnn; digits = 4))   " *
                "($(round((time()-tg)/60; digits = 2)) min)")
    end

    println("\n" * "="^78)
    println("READING")
    println("  FWHM shrinking with finer granularity AND peak noise falling => the peak was BROAD;")
    println("     the paired curve at finer granularity is worth one bounded look.")
    println("  FWHM flat, or noise not improving => limited by photon/draw noise; no granularity")
    println("     change rescues it, and the candidate is dead.")
    println("="^78)

    elapsed = (time() - t0) / 60
    tmp = PEAK_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version   = 1,
        n_units_partA    = nA,
        offsets          = PD_OFFSETS,
        peak_noise_x     = noise_x, peak_noise_y = noise_y,
        aniso_log_ratio  = obs_log, aniso_se = se_log, aniso_sigmas = abs(obs_log)/se_log,
        fwhm_x           = fwhm_x,  fwhm_y = fwhm_y,
        point_noise      = pt_noise, prominence = prom, curve_snr = prom / pt_noise,
        grid_sizes       = collect(PD_GRIDS),
        grid_peak_noise  = [gres[n].peak_noise  for n in PD_GRIDS],
        grid_fwhm        = [gres[n].fwhm        for n in PD_GRIDS],
        grid_point_noise = [gres[n].point_noise for n in PD_GRIDS],
        grid_curve_snr   = [gres[n].snr         for n in PD_GRIDS],
        elapsed_min      = elapsed,
        generated        = string(Dates.now(Dates.UTC)) * "Z",
        caption          = "Phase-11 peak-noise decomposition: curve FWHM vs per-point curve " *
                           "noise, the x/y anisotropy at raised n, and whether computing the " *
                           "PAIRED CURVE at finer granularity sharpens the peak. Gates nothing.",
    )
    let d = JLD2.load(tmp); @assert haskey(d, "fwhm_x") && haskey(d, "grid_fwhm"); end
    mv(tmp, PEAK_REPORT_PATH; force = true)
    println("wrote ", PEAK_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    return nothing
end

main()
