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

# spike/validation/run_p11_probe.jl --- Phase-11 D-06 pre-flight probe (simulator only, no training).
#
# WHAT THIS PROVES, AND WHAT IT DELIBERATELY DOES NOT. At the fixed 8x8 patch grid a patch is a
# LARGE fraction of the image (32 px at 256^2, 256 px at 2048^2), so a few px of misalignment is a
# small relative displacement and, after the PSF blur, the 128-dim summary might barely move. If it
# barely moves, SC2 would fail FLAT -- below the resolution of the fixed summary, NOT because the
# design is wrong. This script measures, with NO trained net and in minutes, how far the summary
# actually moves under injected misalignment, expressed in the units of a KNOWN dRho. It makes NO
# calibration, coverage, width or identifiability claim; nothing here is evidence about any net.
#
# IT PRODUCES THRESHOLDS -- IT DOES NOT GATE ON THEM. `run_sbc.jl` / `run_ood.jl` end in a Test.jl
# test-set block and exit nonzero; this script does not, and must not. Its output is the measurement
# that the Tier-2 append of the pre-registration is DERIVED from (D-06: "the probe supplies the SC2
# threshold; the threshold is not guessed"; D-08: the tolerance band must be justified from this
# probe, not chosen after seeing downstream results). The abort DECISION that reads these numbers
# belongs to a later plan and to the Tier-1 criterion, which was locked before this ran.
#
# IT MUST NEVER WRITE THE PRE-REGISTRATION FILE. The derived values are PRINTED at the end in
# copy-paste-ready `const` form, for a human-authored, append-only second guard block. A script
# that wrote its own thresholds into the file it is scored against would make the pre-registration
# circular (T-11-19). There is no write path from here to that file, by construction.
#
# RESERVED STREAM (D-01). Every draw rides `p11_rng(P11_PROBE_COUNTER)` and the two derived
# per-replicate key families below. That counter is reserved for this probe alone; no reported
# number elsewhere in Phase 11 rides on it, and none of the forbidden reserved streams
# (PROD_SEED_V2, VAL_MASTER_SEED, NPE_MASTER_SEED, VAL_FIX_SEED, CORPUS_MASTER_SEED) is touched.
# This is the SINGLE reported probe run; re-running it for nicer numbers would be a snooped
# stream (T-11-21), and the D-04 iteration allowance is a separate, explicitly declared mechanism.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reaches src/ only READ-ONLY, through
# `spike/contract.jl`'s `patch()` / `correlation()` calls. No file under src/ is written.
#
# CPU-only (D-10). Run:
#     julia --project=spike -t auto spike/validation/run_p11_probe.jl

using JLD2
using Dates
using Statistics
using StatsBase                      # corspearman (never hand-rolled, PATTERNS S6)

# Units under test. ORDER MATTERS: the pre-registration first (it binds every knob below and the
# `p11_rng` stream), then the prior/simulator halves, then the contract and the encoder. Guarded
# for idempotency (PATTERNS S2). NOTE the prior guard is keyed on CHROMATIC_PRIOR, not on the
# prior-sampling function: that function name is asserted to occur exactly ONCE in this file, at
# the frozen theta-base draw, so that no sweep can silently re-draw a theta mid-run.
isdefined(@__MODULE__, :P11_DEV_SEED)     || include(joinpath(@__DIR__, "p11_consts.jl"))
isdefined(@__MODULE__, :CHROMATIC_PRIOR)  || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair)    || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)        || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)       || include(joinpath(@__DIR__, "..", "data", "encode.jl"))

const P11_PROBE_REPORT_PATH = joinpath(@__DIR__, "p11_probe_report.jld2")

# --- Resolution / bookkeeping knobs (NOT pre-registered thresholds) -----------------------
# `run_ood.jl:70-76` voice: these govern how the run is executed, never what it is scored
# against. Every scored knob lives in the pre-registration and is read, never re-typed, below.
const PROBE_PROGRESS_EVERY = 8       # units between progress lines
const PROBE_ARM_F5         = :f5_mixture
const PROBE_ARM_256        = :imsize_256

# --- The rung sets actually evaluated ------------------------------------------------------
# The shift arm is the pre-registered probe ladder UNION the pre-registered SC2 ladder. The union
# is required, not a widening: S_probe is defined as the rank correlation across SC2_RUNGS, so the
# shift arm has to be evaluated at exactly those magnitudes for that statistic to exist at all.
# (SC2_RUNGS contributes 2.5 px; the probe ladder contributes the 0.0 px reference rung.)
const PROBE_SHIFT_ARM = Tuple(sort(unique(Float64[P11_PROBE_SHIFT_RUNGS...,
                                                  P11_PROBE_SHIFT_EXT...,
                                                  SC2_RUNGS...])))
const PROBE_EPS_ARM   = Tuple(sort(unique(Float64[P11_PROBE_EPS_RUNGS...,
                                                  P11_PROBE_EPS_EXT...])))
const PROBE_DRHO_ARM  = Tuple(Float64.(P11_PROBE_DRHO_RUNGS))

# --- Paired keying (LOAD-BEARING, NOT A NICETY) --------------------------------------------
# EVERY RUNG OF A REPLICATE IS SIMULATED FROM A FRESHLY RE-DERIVED RNG OBJECT ON THE SAME KEY, so
# stages 1-5 (the shared latent field, the two private fields, the Bernoulli thinning, the PSF)
# are bit-identical across the sweep and ONLY the stage-6 geometry differs. The measured reason:
#
#     independent-key ||ds||_2 (same theta, different key) = 2.653  (sd 0.132)
#     |shift| = 3.0 px, PAIRED                             = 1.281  (sd 0.162)
#
# THE BETWEEN-KEY MONTE-CARLO FLOOR IS ~2.1x LARGER THAN THE ENTIRE 3 px EFFECT. An unpaired probe
# would report noise and nothing else. (Honest scope of the pairing: it is EXACT through stage 5
# and exact end-to-end for channel 1; channel 2's stage-7 Poisson draw consumes a
# value-dependent amount of the stream, so its read-noise realisation necessarily diverges once
# the warp has changed channel 2. That residual is inside the measured 0.16 sd above.)
#
# Two disjoint key families off the reserved counter, asserted disjoint at run time:
#   pair key   = P11_PROBE_COUNTER * 1000 + r          -> the simulation stream for replicate r
#   imsize key = P11_PROBE_COUNTER * 1000 + 1000 + r   -> the imsize draw for replicate r
# The imsize draw is keyed on the REPLICATE, so it too is paired: one replicate is one image size
# held fixed across every rung, and a rung-to-rung difference can never be an image-size change.
_pair_key(r::Integer)   = P11_PROBE_COUNTER * 1000 + Int(r)
_imsize_key(r::Integer) = P11_PROBE_COUNTER * 1000 + 1000 + Int(r)

"""
    _draw_imsize(rng) -> Tuple{Int,Int}

Weighted draw from the F5 mixture, READ from the pre-registration (never re-typed).
"""
function _draw_imsize(rng)
    u   = rand(rng)
    acc = 0.0
    for (i, w) in enumerate(P11_IMSIZE_WEIGHTS)
        acc += w
        u <= acc && return P11_IMSIZE_SET[i]
    end
    return last(P11_IMSIZE_SET)
end

"""
    _summary(key, theta; imsize) -> Vector{Float64}

The 128-dim summary producer: `simulate_pair` -> `build_mci` -> `patch_summary` -> `encode_d01`.
The RNG object is re-derived from `key` on EVERY call; that is what makes the design paired.
"""
_summary(key::Integer, theta; imsize::Tuple{Int,Int}) =
    encode_d01(patch_summary(build_mci(simulate_pair(p11_rng(key), theta; imsize = imsize))))

"""
    _delta_norm(s_pert, s_base) -> Float64

The pre-registered metric `P11_PROBE_METRIC = :paired_l2_rows_1_64`: the Euclidean norm of the
paired summary difference over the 64 CONTINUOUS rows ONLY.

Rows 1:64 are the imputed per-patch correlations; rows 65:128 are the binary present-mask
(`spike/data/encode.jl:62-64`). The mask rows are near-constant and are EXACTLY the rows
`_summary_row_partition` excludes from standardization -- including them would dilute the
statistic with a near-zero-variance block. All 64 retained rows are the same physical quantity
(a Pearson correlation in [-1,1]) on the same scale, so no z-scoring is applied: that would
introduce a fitted nuisance for no gain, and no fitted pool exists before training anyway.
"""
_delta_norm(s_pert::AbstractVector, s_base::AbstractVector) =
    sqrt(sum(abs2, view(s_pert, 1:64) .- view(s_base, 1:64)))

"""
    _diag_shift(magnitude) -> (dx, dy)

Diagonal injection: `dx = dy = magnitude / sqrt(2)`, so `hypot(dx, dy) == magnitude` and the rung
label is the true displacement in px (the same convention `SC3_SHIFT_RUNGS` documents).
"""
_diag_shift(magnitude::Real) = (magnitude / sqrt(2.0), magnitude / sqrt(2.0))

# --- The ONLY sanctioned theta mutations ---------------------------------------------------
# NAMED SETTERS, NOT INLINE `merge` CALLS, FOR A CORRECTNESS REASON. In Julia a one-field
# `(field = value)` without a trailing comma is NOT a NamedTuple -- it is an assignment
# expression, and `merge(theta, (field = value))` therefore reaches the keyword-argument method
# and throws at run time, inside a thread, after the arm has already started. These three
# one-liners are the only places a rung value enters theta, they carry the comma once, and each
# is unit-tested to change EXACTLY the field it names and nothing else.
_set_shift(theta, magnitude::Real) =
    merge(theta, (shift_dx = magnitude / sqrt(2.0), shift_dy = magnitude / sqrt(2.0)))
_set_eps(theta, e::Real)   = merge(theta, (chromatic_eps = e,))
_set_rho(theta, rho::Real) = merge(theta, (ρ_true = rho,))

"""
    _ols(x, y) -> (slope, intercept, r2)

Ordinary least squares of `y` on `x`, plus the coefficient of determination.
"""
function _ols(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(ArgumentError("_ols: length mismatch"))
    xbar = mean(x); ybar = mean(y)
    sxx  = sum(abs2, x .- xbar)
    sxy  = sum((x .- xbar) .* (y .- ybar))
    slope     = sxy / sxx
    intercept = ybar - slope * xbar
    yhat      = intercept .+ slope .* x
    sstot     = sum(abs2, y .- ybar)
    r2        = sstot == 0 ? 1.0 : 1.0 - sum(abs2, y .- yhat) / sstot
    return (slope = slope, intercept = intercept, r2 = r2)
end

"""
    _drho_eq(norm_value, slope, intercept) -> Float64

THE HEADLINE UNIT. Inverts the calibration line `||ds||_2 ~ slope * dRho + intercept` to answer
"which dRho would displace the summary this far?". This is literally what D-06 asks for -- the
movement measured "relative to the movement induced by a known dRho" -- and it is what makes the
number reviewer-legible ("3 px of misalignment looks like a dRho of about 0.19") instead of an
uncalibrated norm with no units.
"""
_drho_eq(norm_value::Real, slope::Real, intercept::Real) = (norm_value - intercept) / slope

"""
    _rung_stats(A, j) -> (mean, median, sd)

Pooled statistics over every (theta base, replicate) cell at rung index `j`. The MEDIAN is
reported alongside the mean because one degenerate patch (mapped to 0.0 by `encode_d01`) can move
a single replicate's norm by order 1.
"""
function _rung_stats(A::Array{Float64,3}, j::Integer)
    v = vec(@view A[:, :, j])
    return (mean = mean(v), median = median(v), sd = std(v))
end

# --- Atomic persistence of the reported artifact (the save_npe idiom, PATTERNS S4) ----------
function _save_probe_report(path; kwargs...)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "shift_norm_f5") "_save_probe_report: integrity check failed ($tmp)"
        @assert haskey(f, "generated")     "_save_probe_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end

"""
    _run_arm(label, thetas, imsize_of) -> NamedTuple

Sweep all three axes for one image-size arm. `imsize_of(r)` returns the (paired) image size for
replicate `r`. Returns the raw per-(theta, replicate, rung) norm arrays.

THE THREE AXES ARE SWEPT SEPARATELY, ON PURPOSE (ruling Q3). The chromatic term is not
lambda-scaled and gets no conditioning input, so the registration axis and the chromatic axis are
never mixed into one ladder: SC2's monotone-widening claim covers REGISTRATION ONLY, and the
chromatic dose-response is reported on its own axis.
"""
function _run_arm(label::Symbol, thetas::Vector, imsize_of::Function)
    T  = length(thetas)
    R  = P11_PROBE_R
    ns = length(PROBE_SHIFT_ARM)
    ne = length(PROBE_EPS_ARM)
    nd = length(PROBE_DRHO_ARM)

    shift_norm = zeros(Float64, T, R, ns)
    eps_norm   = zeros(Float64, T, R, ne)
    drho_norm  = zeros(Float64, T, R, nd)
    used_w     = zeros(Int, T, R)
    used_h     = zeros(Int, T, R)

    # dRho DIRECTION, fixed per theta base BEFORE any evaluation. rho_true is already near the
    # edge of [-1, 1] for some prior draws, so a blind +dRho would CLAMP at the largest rung and
    # silently flatten the calibration curve. Step toward the interior instead, once per theta,
    # for every rung -- the norm is a magnitude, so the sign of the step does not bias it, but a
    # clamped step would.
    dirs = [ (thetas[t].ρ_true + maximum(PROBE_DRHO_ARM) <= 1.0) ? 1.0 : -1.0 for t in 1:T ]

    units = vec([(t, r) for t in 1:T, r in 1:R])
    done  = Threads.Atomic{Int}(0)
    t_arm = time()
    println("[$label] $(length(units)) paired units x $(1 + ns + ne + nd) evaluations " *
            "on $(Threads.nthreads()) thread(s)")

    Threads.@threads for u in units
        t, r  = u
        theta = thetas[t]
        key   = _pair_key(r)
        sz    = imsize_of(r)
        used_w[t, r] = sz[1]
        used_h[t, r] = sz[2]

        # The paired REFERENCE: the same theta with the geometry zeroed, so every displacement
        # below is caused purely by the injected misalignment and not by the theta's own
        # prior-drawn shift/chromatic values.
        theta0 = _set_eps(_set_shift(theta, 0.0), 0.0)
        s_base = _summary(key, theta0; imsize = sz)

        for (j, m) in enumerate(PROBE_SHIFT_ARM)
            s = _summary(key, _set_shift(theta0, m); imsize = sz)
            shift_norm[t, r, j] = _delta_norm(s, s_base)
        end
        for (j, e) in enumerate(PROBE_EPS_ARM)
            s = _summary(key, _set_eps(theta0, e); imsize = sz)
            eps_norm[t, r, j] = _delta_norm(s, s_base)
        end
        for (j, d) in enumerate(PROBE_DRHO_ARM)
            s = _summary(key, _set_rho(theta0, theta0.ρ_true + dirs[t] * d); imsize = sz)
            drho_norm[t, r, j] = _delta_norm(s, s_base)
        end

        n = Threads.atomic_add!(done, 1) + 1
        if n % PROBE_PROGRESS_EVERY == 0
            println("  [$label] $n/$(length(units)) units  " *
                    "($(round((time() - t_arm) / 60; digits = 2)) min)")
        end
    end

    println("[$label] arm complete in $(round((time() - t_arm) / 60; digits = 2)) min")
    return (shift_norm = shift_norm, eps_norm = eps_norm, drho_norm = drho_norm,
            used_w = used_w, used_h = used_h, dirs = dirs)
end

"""
    _arm_summary(arm) -> NamedTuple

Per-rung statistics, the dRho calibration fit and the dRho_eq column for one arm.
"""
function _arm_summary(arm)
    shift_m = [_rung_stats(arm.shift_norm, j).mean   for j in eachindex(PROBE_SHIFT_ARM)]
    shift_d = [_rung_stats(arm.shift_norm, j).median for j in eachindex(PROBE_SHIFT_ARM)]
    shift_s = [_rung_stats(arm.shift_norm, j).sd     for j in eachindex(PROBE_SHIFT_ARM)]
    eps_m   = [_rung_stats(arm.eps_norm,   j).mean   for j in eachindex(PROBE_EPS_ARM)]
    eps_d   = [_rung_stats(arm.eps_norm,   j).median for j in eachindex(PROBE_EPS_ARM)]
    eps_s   = [_rung_stats(arm.eps_norm,   j).sd     for j in eachindex(PROBE_EPS_ARM)]
    drho_m  = [_rung_stats(arm.drho_norm,  j).mean   for j in eachindex(PROBE_DRHO_ARM)]
    drho_d  = [_rung_stats(arm.drho_norm,  j).median for j in eachindex(PROBE_DRHO_ARM)]
    drho_s  = [_rung_stats(arm.drho_norm,  j).sd     for j in eachindex(PROBE_DRHO_ARM)]

    fit = _ols(collect(PROBE_DRHO_ARM), drho_m)
    eq(v) = _drho_eq(v, fit.slope, fit.intercept)
    return (shift_mean = shift_m, shift_median = shift_d, shift_sd = shift_s,
            eps_mean = eps_m, eps_median = eps_d, eps_sd = eps_s,
            drho_mean = drho_m, drho_median = drho_d, drho_sd = drho_s,
            fit = fit,
            shift_drho_eq = eq.(shift_m), eps_drho_eq = eq.(eps_m),
            drho_drho_eq  = eq.(drho_m))
end

function _print_table(title, rungs, mean_v, median_v, sd_v, eq_v; unit = "")
    println("\n", "-"^78)
    println(title)
    println(rpad("rung$unit", 12), rpad("mean ||ds||", 14), rpad("median", 14),
            rpad("sd", 14), "dRho_eq")
    println("-"^78)
    for (j, x) in enumerate(rungs)
        println(rpad(x, 12),
                rpad(round(mean_v[j];   sigdigits = 4), 14),
                rpad(round(median_v[j]; sigdigits = 4), 14),
                rpad(round(sd_v[j];     sigdigits = 3), 14),
                round(eq_v[j]; sigdigits = 4))
    end
    println("-"^78)
end

function main()
    t_start = time()
    println("="^78)
    println("Phase-11 D-06 PRE-FLIGHT PROBE — simulator only, NO training, NO verdict")
    println("  theta bases=$P11_PROBE_THETA_BASE_N  R=$P11_PROBE_R  metric=$P11_PROBE_METRIC")
    println("  shift arm (px) = $PROBE_SHIFT_ARM")
    println("  eps arm        = $PROBE_EPS_ARM")
    println("  dRho arm       = $PROBE_DRHO_ARM")
    println("  arms: F5 mixture $P11_IMSIZE_SET w=$P11_IMSIZE_WEIGHTS  +  " *
            "comparability $P11_PROBE_COMPARABILITY_IMSIZE")
    println("  stream = p11_rng(P11_PROBE_COUNTER=$P11_PROBE_COUNTER) off " *
            "P11_DEV_SEED=$(repr(P11_DEV_SEED))")
    println("="^78)

    # Seed discipline, asserted rather than asserted-in-prose (D-01).
    @assert !(UInt64(P11_DEV_SEED) in _p11_forbidden()) "probe must not ride a forbidden stream"
    @assert P11_PROBE_COUNTER ∉ (P11_DATAGEN_COUNTER, P11_LADDER_COUNTER, P11_BREAKDOWN_COUNTER,
                                 P11_ATTENUATION_COUNTER, P11_REALIMAGE_COUNTER,
                                 P11_FIXTURE_COUNTER) "probe counter must be reserved"
    @assert isempty(intersect(Set(_pair_key.(1:P11_PROBE_R)),
                              Set(_imsize_key.(1:P11_PROBE_R)))) "key families must be disjoint"
    @assert P11_PROBE_COUNTER ∉ _pair_key.(1:P11_PROBE_R)
    @assert P11_PROBE_METRIC === :paired_l2_rows_1_64 "the metric is pre-registered, not chosen here"

    # --- (1) The frozen theta base points ---------------------------------------------------
    # Drawn ONCE on the reserved counter, then FROZEN as explicit NamedTuples for the whole run.
    # Every sweep below mutates only shift_dx / shift_dy / chromatic_eps / rho_true through
    # `merge`, never by drawing again: a fresh draw consumes the stream differently and its
    # nuisances would no longer be held fixed, which is the one thing this design requires.
    rng_theta = p11_rng(P11_PROBE_COUNTER)
    thetas = [sample_prior(rng_theta) for _ in 1:P11_PROBE_THETA_BASE_N]
    thetas = [(ρ_true = t.ρ_true, spillover = t.spillover, autofluorescence = t.autofluorescence,
               label_efficiency = t.label_efficiency, shift_dx = 0.0, shift_dy = 0.0,
               noise = t.noise, chromatic_eps = 0.0) for t in thetas]
    println("\nfrozen theta base points (geometry zeroed; nuisances held fixed for the whole run):")
    for (i, t) in enumerate(thetas)
        println("  [$i] rho_true=$(round(t.ρ_true; digits=4))  spill=$(round(t.spillover; digits=4))  " *
                "autofl=$(round(t.autofluorescence; digits=4))  " *
                "label_eff=$(round(t.label_efficiency; digits=4))  noise=$(round(t.noise; digits=4))")
    end

    # --- (2) Both arms ----------------------------------------------------------------------
    imsizes_f5 = [_draw_imsize(p11_rng(_imsize_key(r))) for r in 1:P11_PROBE_R]
    arm_f5  = _run_arm(PROBE_ARM_F5,  thetas, r -> imsizes_f5[r])
    arm_256 = _run_arm(PROBE_ARM_256, thetas, _ -> P11_PROBE_COMPARABILITY_IMSIZE)

    sum_f5  = _arm_summary(arm_f5)
    sum_256 = _arm_summary(arm_256)

    # --- (3) The two derived quantities the Tier-2 append needs -----------------------------
    # Evaluated on the F5 MIXTURE arm: the ladder trains and evaluates on that mixture (F5
    # binding invariant), so its dRho_eq axis is the one the ladder's numbers will live on.
    sc2_raw      = [findfirst(==(x), PROBE_SHIFT_ARM) for x in SC2_RUNGS]
    @assert !any(isnothing, sc2_raw) "every SC2 rung must be evaluated on the shift arm"
    sc2_idx      = Int[i for i in sc2_raw]
    sc2_eq_f5    = sum_f5.shift_drho_eq[sc2_idx]
    sc2_eq_256   = sum_256.shift_drho_eq[sc2_idx]
    s_probe      = corspearman(collect(Float64.(SC2_RUNGS)), sc2_eq_f5)
    ladder_span  = maximum(sc2_eq_f5) - minimum(sc2_eq_f5)
    lambda_ratio = sc2_eq_f5[end] / sc2_eq_f5[1]

    sc2_floor       = SC2_SPEARMAN_ATTENUATION * s_probe
    ablation_raw    = SC2_SPEARMAN_ATTENUATION * lambda_ratio
    ablation_factor = max(1.05, ablation_raw)

    # Pairing determinism: the 0.0 rung of each arm is the reference itself, re-simulated from a
    # freshly re-derived RNG object. It MUST come back exactly 0.0; anything else means the paired
    # re-derivation is not bit-identical and the whole design is invalid. Reported, never thrown.
    zi_s = findfirst(==(0.0), PROBE_SHIFT_ARM)
    zi_e = findfirst(==(0.0), PROBE_EPS_ARM)
    pairing_zero_max = maximum((maximum(arm_f5.shift_norm[:, :, zi_s]),
                               maximum(arm_f5.eps_norm[:, :, zi_e]),
                               maximum(arm_256.shift_norm[:, :, zi_s]),
                               maximum(arm_256.eps_norm[:, :, zi_e])))

    elapsed_min = (time() - t_start) / 60

    # --- (4) Persist BEFORE anything is interpreted ------------------------------------------
    # The artifact must survive an honest "below resolution" result exactly as it survives a good
    # one, so it is written before a single number is compared against a floor.
    _save_probe_report(P11_PROBE_REPORT_PATH;
        # raw per-(theta, replicate, rung) tables
        shift_norm_f5 = arm_f5.shift_norm,  eps_norm_f5 = arm_f5.eps_norm,
        drho_norm_f5 = arm_f5.drho_norm,
        shift_norm_256 = arm_256.shift_norm, eps_norm_256 = arm_256.eps_norm,
        drho_norm_256 = arm_256.drho_norm,
        imsize_w_f5 = arm_f5.used_w, imsize_h_f5 = arm_f5.used_h,
        drho_direction = arm_f5.dirs,
        theta_bases = [collect(values(t)) for t in thetas],
        theta_fields = collect(String.(keys(thetas[1]))),
        # rung labels
        shift_rungs = collect(PROBE_SHIFT_ARM), eps_rungs = collect(PROBE_EPS_ARM),
        drho_rungs  = collect(PROBE_DRHO_ARM),  sc2_rungs = collect(Float64.(SC2_RUNGS)),
        # per-rung statistics, as parallel Vectors
        shift_mean_f5 = sum_f5.shift_mean, shift_median_f5 = sum_f5.shift_median,
        shift_sd_f5 = sum_f5.shift_sd, shift_drho_eq_f5 = sum_f5.shift_drho_eq,
        eps_mean_f5 = sum_f5.eps_mean, eps_median_f5 = sum_f5.eps_median,
        eps_sd_f5 = sum_f5.eps_sd, eps_drho_eq_f5 = sum_f5.eps_drho_eq,
        drho_mean_f5 = sum_f5.drho_mean, drho_median_f5 = sum_f5.drho_median,
        drho_sd_f5 = sum_f5.drho_sd, drho_drho_eq_f5 = sum_f5.drho_drho_eq,
        shift_mean_256 = sum_256.shift_mean, shift_median_256 = sum_256.shift_median,
        shift_sd_256 = sum_256.shift_sd, shift_drho_eq_256 = sum_256.shift_drho_eq,
        eps_mean_256 = sum_256.eps_mean, eps_median_256 = sum_256.eps_median,
        eps_sd_256 = sum_256.eps_sd, eps_drho_eq_256 = sum_256.eps_drho_eq,
        drho_mean_256 = sum_256.drho_mean, drho_median_256 = sum_256.drho_median,
        drho_sd_256 = sum_256.drho_sd, drho_drho_eq_256 = sum_256.drho_drho_eq,
        # the calibration fits
        drho_eq_slope_f5 = sum_f5.fit.slope, drho_eq_intercept_f5 = sum_f5.fit.intercept,
        drho_eq_r2_f5 = sum_f5.fit.r2,
        drho_eq_slope_256 = sum_256.fit.slope, drho_eq_intercept_256 = sum_256.fit.intercept,
        drho_eq_r2_256 = sum_256.fit.r2,
        # the headline derived quantities
        S_probe = s_probe, ladder_span = ladder_span,
        sc2_drho_eq_f5 = sc2_eq_f5, sc2_drho_eq_256 = sc2_eq_256,
        lambda_ratio = lambda_ratio,
        derived_SC2_SPEARMAN_FLOOR = sc2_floor,
        derived_P11_LAMBDA_ABLATION_FACTOR = ablation_factor,
        pairing_zero_max = pairing_zero_max,
        # a copy of every governing constant
        P11_PROBE_THETA_BASE_N = P11_PROBE_THETA_BASE_N,
        P11_PROBE_R = P11_PROBE_R,
        P11_PROBE_SHIFT_RUNGS = collect(P11_PROBE_SHIFT_RUNGS),
        P11_PROBE_SHIFT_EXT   = collect(P11_PROBE_SHIFT_EXT),
        P11_PROBE_EPS_RUNGS   = collect(P11_PROBE_EPS_RUNGS),
        P11_PROBE_EPS_EXT     = collect(P11_PROBE_EPS_EXT),
        P11_PROBE_DRHO_RUNGS  = collect(P11_PROBE_DRHO_RUNGS),
        P11_PROBE_METRIC = String(P11_PROBE_METRIC),
        P11_PROBE_S_FLOOR = P11_PROBE_S_FLOOR, P11_PROBE_SPAN_FLOOR = P11_PROBE_SPAN_FLOOR,
        SC2_SPEARMAN_ATTENUATION = SC2_SPEARMAN_ATTENUATION,
        P11_IMSIZE_SET = [collect(s) for s in P11_IMSIZE_SET],
        P11_IMSIZE_WEIGHTS = collect(P11_IMSIZE_WEIGHTS),
        P11_PROBE_COMPARABILITY_IMSIZE = collect(P11_PROBE_COMPARABILITY_IMSIZE),
        P11_PROBE_COUNTER = P11_PROBE_COUNTER,
        seed = UInt64(P11_DEV_SEED), salt = UInt64(P11_SALT),
        nthreads = Threads.nthreads(),
        caption = "D-06 pre-flight probe: paired summary displacement under injected " *
                  "registration and chromatic misalignment, expressed as the equivalent dRho. " *
                  "Simulator only, no trained net; produces thresholds, gates nothing.",
        elapsed_min = elapsed_min,
        generated = string(Dates.now(Dates.UTC)) * "Z")
    println("\npersisted reported artifact -> $P11_PROBE_REPORT_PATH")

    # --- (5) Report ---------------------------------------------------------------------------
    for (nm, s) in ((PROBE_ARM_F5, sum_f5), (PROBE_ARM_256, sum_256))
        _print_table("[$nm] SHIFT sweep (diagonal, dx = dy = |s|/sqrt(2))",
                     PROBE_SHIFT_ARM, s.shift_mean, s.shift_median, s.shift_sd, s.shift_drho_eq;
                     unit = " (px)")
        _print_table("[$nm] CHROMATIC sweep (shift = 0; SEPARATE AXIS, not lambda-scaled)",
                     PROBE_EPS_ARM, s.eps_mean, s.eps_median, s.eps_sd, s.eps_drho_eq)
        _print_table("[$nm] dRho CALIBRATION sweep (shift = 0, chromatic_eps = 0)",
                     PROBE_DRHO_ARM, s.drho_mean, s.drho_median, s.drho_sd, s.drho_drho_eq)
        println("[$nm] calibration fit: ||ds||_2 = $(round(s.fit.slope; sigdigits=6)) * dRho + " *
                "$(round(s.fit.intercept; sigdigits=6))   R^2 = $(round(s.fit.r2; sigdigits=6))")
    end

    println("\n", "="^78)
    println("READING NOTES (pre-registered annotations, not post-hoc caveats)")
    println("="^78)
    println("* THE 0 -> 0.25 px STEP IS INTERPOLATION ONSET, NOT SIGNAL. At exactly integer")
    println("  offsets the backward map samples on-grid and BSpline(Linear()) is an identity")
    println("  lookup; any sub-pixel offset engages the interpolation kernel. That is why")
    println("  LAMBDA_MIN = 0.25 and why no lambda = 0 rung exists in SC2_RUNGS.")
    println("* SC2's MONOTONE-WIDENING CLAIM COVERS REGISTRATION ONLY. The chromatic term is not")
    println("  lambda-scaled and gets no conditioning input, so its dose-response above is a")
    println("  SEPARATE reported axis and is never part of the SC2 ladder.")
    println("* The shift/chromatic columns of theta stay VACUOUS by design (D-05). Nothing here")
    println("  claims dx/dy/chromatic_eps are identified; the probe measures DISPLACEMENT of the")
    println("  summary, not recoverability of the parameter.")
    println("* PAIRING DETERMINISM: max ||ds|| at the zero rungs = $pairing_zero_max " *
            "(must be exactly 0.0)")

    println("\n", "="^78)
    println("DERIVED QUANTITIES vs THE TIER-1 FLOORS (comparison only — the branch decision and")
    println("the abort verdict belong to a later plan, not to this script)")
    println("="^78)
    println(rpad("quantity", 34), rpad("measured", 20), rpad("Tier-1 floor", 16), "comparison")
    println("-"^78)
    println(rpad("S_probe (Spearman, F5 arm)", 34), rpad(round(s_probe; sigdigits = 6), 20),
            rpad(P11_PROBE_S_FLOOR, 16), s_probe < P11_PROBE_S_FLOOR ? "BELOW FLOOR" : "at or above")
    println(rpad("ladder span in dRho_eq (F5)", 34), rpad(round(ladder_span; sigdigits = 6), 20),
            rpad(P11_PROBE_SPAN_FLOOR, 16),
            ladder_span < P11_PROBE_SPAN_FLOOR ? "BELOW FLOOR" : "at or above")
    println("-"^78)
    println("Tier-1 abort criterion (S_probe < floor OR span < floor): ",
            (s_probe < P11_PROBE_S_FLOOR || ladder_span < P11_PROBE_SPAN_FLOOR) ?
            "WOULD FIRE" : "would NOT fire")
    println("dRho_eq across SC2_RUNGS (F5 arm) : ", round.(sc2_eq_f5; sigdigits = 4))
    println("dRho_eq across SC2_RUNGS (256 arm): ", round.(sc2_eq_256; sigdigits = 4))
    println("lambda_max / lambda_min dRho_eq ratio (F5): ", round(lambda_ratio; sigdigits = 6))

    # --- (6) Copy-paste-ready Tier-2 block ----------------------------------------------------
    # PRINTED, NEVER WRITTEN. A human appends these to the pre-registration as a SECOND guard
    # block; this script has no write path to that file (T-11-19).
    println("\n", "="^78)
    println("TIER-2 CONSTANTS — copy-paste into a SECOND guard block; APPEND, NEVER EDIT")
    println("="^78)
    println("    const P11_PROBE_S_MEASURED     = $(s_probe)")
    println("    const P11_PROBE_SPAN_MEASURED  = $(ladder_span)")
    println("    const SC2_SPEARMAN_FLOOR       = SC2_SPEARMAN_ATTENUATION * P11_PROBE_S_MEASURED")
    println("    const P11_DRHO_EQ_SLOPE        = $(sum_f5.fit.slope)")
    println("    const P11_DRHO_EQ_INTERCEPT    = $(sum_f5.fit.intercept)")
    println("    const P11_DRHO_EQ_R2           = $(sum_f5.fit.r2)")
    println("    const P11_DRHO_EQ_SLOPE_256    = $(sum_256.fit.slope)")
    println("    const P11_DRHO_EQ_INTERCEPT_256 = $(sum_256.fit.intercept)")
    println("    const P11_DRHO_EQ_R2_256       = $(sum_256.fit.r2)")
    println("    const P11_LADDER_IMSIZE_ARM    = :f5_mixture")
    println("    const P11_LAMBDA_ABLATION_FACTOR = $(ablation_factor)")
    println("    #   raw ratio dRho_eq(LAMBDA_MAX)/dRho_eq(LAMBDA_MIN) = $(lambda_ratio)")
    println("    #   attenuated: SC2_SPEARMAN_ATTENUATION * ratio      = $(ablation_raw)")
    println("    #   floored at 1.05                                   = $(ablation_factor)")
    println("="^78)
    println("elapsed_min = $(round(elapsed_min; digits = 2))   " *
            "generated (UTC) recorded in the artifact")
    println("="^78)
    return nothing
end

# --- The reported run ----------------------------------------------------------------------
# GUARDED, DELIBERATELY, AGAINST THE HOUSE bare-`main()` SHAPE (run_sbc.jl:177). The helpers above
# are unit-tested by `spike/test/test_p11_probe.jl`, which must include THIS file to reach them; a
# bare call would make that unit test re-run the reported probe on the RESERVED stream (T-11-21)
# and take ~10 minutes. The test defines `P11_PROBE_LOAD_ONLY` before including; nothing else in
# the repository defines that name, so the reported invocation below is otherwise unconditional.
isdefined(@__MODULE__, :P11_PROBE_LOAD_ONLY) || main()
