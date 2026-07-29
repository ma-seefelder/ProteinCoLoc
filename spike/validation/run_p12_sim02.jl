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

# spike/validation/run_p12_sim02.jl --- the REPORTED per-region SIM-02 claim (D-05).
#
# THE CLAIM. Under D-05's copula -- z (marginally standard per cell) -> Phi -> the MU_PRIOR
# quantile -> the frozen ghat, elementwise -- EVERY ONE of the G^2 = 64 regions carries EXACTLY
# the `MU_PRIOR` marginal that `prior.jl:40` puts on the scalar knob. The lattice prior therefore
# contributes ONLY SPATIAL DEPENDENCE, and the CLAUDE.md constraint that pi(theta) stay consistent
# with the Turing mu/nu/sigma/tau ranges -- hence comparable with the ADVI baseline -- is untouched
# by anything Phase 12 adds. This runner is where that stops being a screening test
# (`test_p12_prior.jl` at M = 2 000) and becomes the phase's committed evidence.
#
# THE BAR, AND WHERE IT CAME FROM. `P12_SIM02_W1_TOL_PERREGION` is declared in TIER 1
# (`p12_consts.jl` section 5), is IDENTICAL to `prior.jl:46`'s `SIM02_W1_TOL = 0.10`, and was fixed
# BEFORE the field simulator existed. It is NOT reverse-engineered from any measured value, and no
# numeric tolerance literal appears anywhere in this file: the bar is READ, never re-derived.
# 12-RESEARCH.md Pattern 2's measured 0.00732 is a FORECAST of what this run will find -- a thing
# the run either confirms or refutes, never the thing that sets the bar.
#
# THE ASSERTION IS ON THE MAXIMUM OVER THE 64 REGIONS, NEVER THE MEAN. A mean lets one
# badly-mismatched region hide behind 63 good ones, which is exactly the boundary failure R-8
# identified: without the per-cell `./ sqrt.(diag(Sigma))` rescale in `field_sampler` the copula
# pushes a NON-standard normal through Phi, so the per-region prior on rho becomes
# POSITION-DEPENDENT -- the radial/edge artifact class this phase exists to detect. The mean is
# computed and reported; it gates nothing.
#
# WHY THE WHOLE LADDER, NOT ONE RUNG. 12-07 measured, BY MUTATION, that at `r1 = 0.50` a
# rescale-deleted copula reads a max per-region W1 of 0.0789 -- badly wrong, yet INSIDE the 0.10
# bar. The same deletion breaches the bar at the BOTTOM `P12_R1_LADDER` rung (0.1264). Sweeping
# every pre-registered rung of both spatial arms is therefore what makes this assertion able to
# FAIL. Nothing was moved to achieve that: the rungs are pre-registered, the bar is pre-registered,
# and only the decision to screen at all of them rather than one is this runner's.
#
# TWO SCOPINGS ARE REPORTED, AND BOTH ARE HELD TO THE SAME PRE-REGISTERED BAR.
#   `w1_*`      -- SCOPED to [GHAT_MU_MIN, GHAT_MU_MAX], the D-16 realized range, exactly as
#                  `prior.jl:42-46` scopes its own induced-mu claim and as `test_p12_prior.jl:72-88`
#                  implements that scoping. This is the pre-registered quantity and the one the
#                  research forecast describes.
#   `w1_full_*` -- the SAME estimator over the FULL [-1, 1] support of `MU_PRIOR`. It is an
#                  ADDITION, not a second threshold: it is held to the SAME
#                  `P12_SIM02_W1_TOL_PERREGION`. It exists because 12-07 recorded that the scoping
#                  DAMPS the very signal this claim is about -- a compressed marginal moves most of
#                  its mass into precisely the tails the scoping removes. Measured on the
#                  rescale-deleted mutant (fixture stream, M = 2 000): the scoped max reads
#                  0.0789 at r1 = 0.50 and PASSES, while the full-range max reads 0.1152 and FAILS,
#                  and at r1 = 0.95 the two read 0.1032 against 0.3328.
#
# WHAT ELSE IS REPORTED. The per-region `ghat` atom mass -- the fraction of `rho_field` entries
# sitting exactly on the frozen clamp endpoints `GHAT_RHO_KNOTS[1]` / `[end]`. It is
# REPORTING-ONLY AND GATES NOTHING, and the artifact caption says so, because recording a merely
# comparable number as a gate is precisely how the Phase-7 over-powered-chi-squared artifact was
# manufactured. Also reported: the empirical per-cell sd of the Gaussian field the copula consumes
# (must be ~1 everywhere -- it is the rescale, made visible), and the four corner regions by name.
#
# DECOUPLING (CLAUDE.md, hard constraint). PRIOR-ONLY: no image simulation, no `simulate_pair`, no
# network, no training, no `src/` file, no package added, no constant mutated, and nothing written
# outside this file's own artifact. `Statistics` is a stdlib reached through the default `@stdlib`
# LOAD_PATH entry exactly as `run_p11_recovery.jl:50-53` reaches `LinearAlgebra`; `Distributions`
# arrives transitively through `p12_prior.jl` -> `prior.jl`, which is where `MU_PRIOR` itself comes
# from, so this file cannot fork the prior by restating it.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p12_sim02.jl

using JLD2
using Dates
using Statistics

# ORDER MATTERS: the Tier-1 pre-registration first (it carries M, the bar, the counter, the ladder
# and the RNG), then the copula under test. Guarded, the `run_p11_recovery.jl:55-56` idiom.
isdefined(@__MODULE__, :P12_SIM02_M)      || include(joinpath(@__DIR__, "p12_consts.jl"))
isdefined(@__MODULE__, :sample_p12_prior) || include(joinpath(@__DIR__, "..", "simulator", "p12_prior.jl"))

const SIM02_REPORT_PATH = joinpath(@__DIR__, "p12_sim02_report.jld2")

# --- The scoped reference, COPIED rather than re-derived --------------------------------------
# Verbatim from `test_p12_prior.jl:77`, which in turn implements `prior.jl:42-46`'s D-16 scoping.
# Conditioning a truncated Cauchy on a sub-interval is again a truncated Cauchy, so restricting
# BOTH sides to the realized range is exact rather than an approximation and the comparison stays
# like-for-like.
const _P12S_REF_SCOPED = Truncated(Cauchy(0.0, 0.3), GHAT_MU_MIN, GHAT_MU_MAX)
# ... and the restatement cannot silently FORK pi(theta), which is the one thing `p12_prior.jl`'s
# header forbids: the base distribution is asserted identical to `MU_PRIOR`'s own, so a Tier-1 or
# `prior.jl` change that moved the mu-prior would fail this file at load rather than quietly leave
# the reported W1 measured against a stale target.
@assert _P12S_REF_SCOPED.untruncated == MU_PRIOR.untruncated

# --- Wasserstein-1 by quantile transport ------------------------------------------------------
# A DELIBERATE SECOND COPY of `calibration.jl:137-143`'s estimator, with the scoping copied from
# `test_p12_prior.jl:72-88` rather than re-derived. `_w1` is not factored out into anything
# includable: reaching it means including `calibration.jl` whole, which re-binds its OWN `MU_PRIOR`
# and `SIM02_W1_TOL` and pulls `contract.jl` -> `src/` and `forward.jl` into a prior-only runner.
# Forking pi(theta) to save nine lines is the trade this file refuses.
#
# The reference is given ANALYTICALLY instead of as a second finite sample (12-07's precedent): the
# target is known in closed form, so drawing one would only add Monte-Carlo noise to the reported
# statistic.
function _p12s_w1_scoped(x::AbstractVector{<:Real})
    s = sort([v for v in x if GHAT_MU_MIN <= v <= GHAT_MU_MAX])
    n = length(s)
    n < 5 && return NaN
    return mean(abs(s[k] - quantile(_P12S_REF_SCOPED, (k - 0.5) / n)) for k in 1:n)
end

"The same estimator against the UNSCOPED `MU_PRIOR`, over its full [-1, 1] support."
function _p12s_w1_full(x::AbstractVector{<:Real})
    s = sort(collect(float.(x)))
    n = length(s)
    n < 5 && return NaN
    return mean(abs(s[k] - quantile(MU_PRIOR, (k - 0.5) / n)) for k in 1:n)
end

# The four corner regions, by name, in `p12_idx` (column-major) order.
const _P12S_CORNER_NAMES = ["(1,1)", "(1,$(P12_G))", "($(P12_G),1)", "($(P12_G),$(P12_G))"]
const _P12S_CORNERS      = [p12_idx(1, 1), p12_idx(1, P12_G), p12_idx(P12_G, 1),
                            p12_idx(P12_G, P12_G)]

"""
    _p12s_configs() -> Vector{Tuple{Symbol,Float64}}

Every (arm, r1) this runner reports: both spatial arms across the five pre-registered
`P12_R1_LADDER` rungs, plus the D-10 matched ablation ONCE at `P12_ABLATION_R1`. Running all three
arms here is cheap and it is what makes 12-15's CAR-vs-GP comparison a comparison of SPATIAL
DEPENDENCE rather than of two differently-mis-specified marginals.
"""
function _p12s_configs()
    cfg = Tuple{Symbol,Float64}[]
    for arm in (:car, :gp), r1 in P12_R1_LADDER
        push!(cfg, (arm, float(r1)))
    end
    push!(cfg, (:none, float(P12_ABLATION_R1)))
    return cfg
end

_p12s_label(arm::Symbol, r1::Real) = "$(arm)@r1=$(r1)"

"""
    _p12s_draw(rng, arm, r1, M) -> NamedTuple

`M` draws of the D-05 copula at one (arm, r1), returning the per-region induced-mu matrix, the
per-cell sd of the Gaussian field the copula consumed, and the per-region atom count.

The sampler is HOISTED: `sample_p12_prior` re-solves the r1 bisection on every call, which is right
for a datagen loop of independent r1 draws and pointless for a fixed-r1 Monte Carlo. The object
under test -- `rho_field` -- is identical either way (12-07's testsets 1-3 make the same move).

ATOMS ARE COUNTED BY EXACT EQUALITY AGAINST THE FROZEN CLAMP ENDPOINTS, not by a hand-written
endpoint literal with an epsilon: `ghat` RETURNS `GHAT_RHO_KNOTS[1]` / `[end]` unchanged when mu
falls outside the realized range, so the atoms are those Float64s exactly, and reading the knots
means this count cannot drift from the frozen map.
"""
function _p12s_draw(rng, arm::Symbol, r1::Real, M::Integer)
    n = P12_G^2
    f = field_sampler(lattice_sigma(arm, r1))
    MU = Matrix{Float64}(undef, n, M)
    ZZ = Matrix{Float64}(undef, n, M)
    atoms = zeros(Int, n)
    lo, hi = GHAT_RHO_KNOTS[1], GHAT_RHO_KNOTS[end]
    for k in 1:M
        d = rho_field(rng, f)
        rv = vec(d.rho)
        MU[:, k] = vec(d.mu)
        ZZ[:, k] = vec(d.z)
        @inbounds for r in 1:n
            (rv[r] == lo || rv[r] == hi) && (atoms[r] += 1)
        end
    end
    return (mu = MU, zsd = vec(std(ZZ; dims = 2)), atom_mass = atoms ./ M)
end

function main()
    println("="^100)
    println("PHASE-12 PER-REGION SIM-02  (D-05 copula; M = $(P12_SIM02_M) draws per arm x rung)")
    println("GATED: max over the $(P12_G^2) regions of the per-region W1 vs MU_PRIOR, against the")
    println("       TIER-1 bar P12_SIM02_W1_TOL_PERREGION = $(P12_SIM02_W1_TOL_PERREGION), fixed")
    println("       before the field simulator existed. REPORTED, gating nothing: the atom mass.")
    println("="^100)
    t_start = time()

    # ONE STREAM, ONE COUNTER. Every draw below rides `p12_rng(P12_SIM02_COUNTER)`, the reserved
    # REPORTED sub-stream for this runner (p12_consts.jl section 2), consumed sequentially -- so the
    # run is deterministic and independent of the thread count.
    rng = p12_rng(P12_SIM02_COUNTER)

    configs   = _p12s_configs()
    labels    = String[]
    cfg_arms  = String[]
    cfg_r1    = Float64[]
    w1_pr     = Dict{String,Vector{Float64}}()
    w1full_pr = Dict{String,Vector{Float64}}()
    sd_pr     = Dict{String,Vector{Float64}}()
    atom_pr   = Dict{String,Vector{Float64}}()
    corners   = Dict{String,Vector{Float64}}()
    w1_max    = Float64[]
    w1_mean   = Float64[]
    w1f_max   = Float64[]
    w1f_mean  = Float64[]
    sd_min    = Float64[]
    sd_max    = Float64[]
    atom_mean = Float64[]
    atom_max  = Float64[]

    for (arm, r1) in configs
        lab = _p12s_label(arm, r1)
        d   = _p12s_draw(rng, arm, r1, P12_SIM02_M)
        wv  = [_p12s_w1_scoped(view(d.mu, r, :)) for r in 1:size(d.mu, 1)]
        wf  = [_p12s_w1_full(view(d.mu, r, :))   for r in 1:size(d.mu, 1)]

        push!(labels, lab); push!(cfg_arms, string(arm)); push!(cfg_r1, float(r1))
        w1_pr[lab]     = wv
        w1full_pr[lab] = wf
        sd_pr[lab]     = d.zsd
        atom_pr[lab]   = d.atom_mass
        corners[lab]   = wv[_P12S_CORNERS]
        push!(w1_max,  maximum(wv));  push!(w1_mean,  mean(wv))
        push!(w1f_max, maximum(wf));  push!(w1f_mean, mean(wf))
        push!(sd_min,  minimum(d.zsd)); push!(sd_max, maximum(d.zsd))
        push!(atom_mean, mean(d.atom_mass)); push!(atom_max, maximum(d.atom_mass))
    end

    # --- The reported table, one block per arm, one row per rung -------------------------------
    for arm in ("car", "gp", "none")
        idx = findall(==(arm), cfg_arms)
        isempty(idx) && continue
        println()
        println("-"^100)
        println("ARM = :$arm")
        println("-"^100)
        println(rpad("r1", 8), rpad("w1_max", 12), rpad("w1_mean", 12), rpad("w1full_max", 13),
                rpad("sd_min", 10), rpad("sd_max", 10), rpad("atom_mean", 12), "verdict")
        for c in idx
            pass = (w1_max[c] <= P12_SIM02_W1_TOL_PERREGION) &&
                   (w1f_max[c] <= P12_SIM02_W1_TOL_PERREGION)
            println(rpad(cfg_r1[c], 8),
                    rpad(round(w1_max[c];   digits = 6), 12),
                    rpad(round(w1_mean[c];  digits = 6), 12),
                    rpad(round(w1f_max[c];  digits = 6), 13),
                    rpad(round(sd_min[c];   digits = 4), 10),
                    rpad(round(sd_max[c];   digits = 4), 10),
                    rpad(round(atom_mean[c]; digits = 6), 12),
                    pass ? "PASS" : "FAIL")
        end
        for c in idx
            println("  corners ", labels[c], "  ",
                    join(("$(_P12S_CORNER_NAMES[i]) = $(round(corners[labels[c]][i]; digits = 6))"
                          for i in eachindex(_P12S_CORNER_NAMES)), "   "))
        end
    end

    # --- FAIL LOUDLY. This runner is one of the two that may fail the phase -----------------------
    # The gate is `maximum` over the 64 regions, per arm x rung, against the Tier-1 bar. It is an
    # @assert and NOT a printed warning on purpose: a reported claim that can only be disproved by
    # someone reading the console is not a claim.
    for (c, lab) in enumerate(labels)
        wv = w1_pr[lab]
        @assert maximum(wv) <= P12_SIM02_W1_TOL_PERREGION (
            "SIM-02 per-region FAILURE (scoped): arm = $(cfg_arms[c]), r1 = $(cfg_r1[c]), " *
            "region $(argmax(wv)) reads W1 = $(maximum(wv)) against the Tier-1 bar " *
            "P12_SIM02_W1_TOL_PERREGION = $(P12_SIM02_W1_TOL_PERREGION). Do NOT relax the bar: a " *
            "failure at the CORNERS points at the sqrt.(diag(Sigma)) rescale in field_sampler " *
            "(R-8, the known mechanism), and a failure UNIFORM across regions points at the " *
            "quantile/Phi composition in rho_field, not at the lattice.")
        wf = w1full_pr[lab]
        @assert maximum(wf) <= P12_SIM02_W1_TOL_PERREGION (
            "SIM-02 per-region FAILURE (full range): arm = $(cfg_arms[c]), r1 = $(cfg_r1[c]), " *
            "region $(argmax(wf)) reads W1 = $(maximum(wf)) against the SAME Tier-1 bar " *
            "P12_SIM02_W1_TOL_PERREGION = $(P12_SIM02_W1_TOL_PERREGION). This is the same claim " *
            "over the full MU_PRIOR support; see the header on why both scopings are reported.")
    end

    elapsed = (time() - t_start) / 60
    tmp = SIM02_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version       = 1,
        generated            = string(Dates.now(Dates.UTC)) * "Z",
        elapsed_min          = elapsed,
        M                    = P12_SIM02_M,
        arms                 = ["car", "gp", "none"],
        r1_ladder            = collect(Float64, P12_R1_LADDER),
        ablation_r1          = float(P12_ABLATION_R1),
        counter              = P12_SIM02_COUNTER,
        tolerance            = P12_SIM02_W1_TOL_PERREGION,
        tolerance_source     = "Tier 1, p12_consts.jl, fixed before the field simulator existed",
        config_labels        = labels,
        config_arms          = cfg_arms,
        config_r1            = cfg_r1,
        w1_per_region        = w1_pr,
        w1_max               = w1_max,
        w1_mean              = w1_mean,
        w1_full_per_region   = w1full_pr,
        w1_full_max          = w1f_max,
        w1_full_mean         = w1f_mean,
        w1_corners           = corners,
        corner_names         = _P12S_CORNER_NAMES,
        corner_regions       = _P12S_CORNERS,
        field_sd_per_region  = sd_pr,
        field_sd_min         = sd_min,
        field_sd_max         = sd_max,
        atom_mass_per_region = atom_pr,
        atom_mass_mean       = atom_mean,
        atom_mass_max        = atom_max,
        caption              = "Phase-12 per-region SIM-02 (D-05): the max over the 64 regions of " *
                               "the per-region induced-mu Wasserstein-1 distance to MU_PRIOR is " *
                               "the GATED quantity, against P12_SIM02_W1_TOL_PERREGION, a Tier-1 " *
                               "bar fixed before the field simulator existed and identical to " *
                               "prior.jl's SIM02_W1_TOL. `w1_*` is scoped to [GHAT_MU_MIN, " *
                               "GHAT_MU_MAX] (D-16); `w1_full_*` is the same estimator over the " *
                               "full MU_PRIOR support and is held to the SAME bar. The per-region " *
                               "ghat atom mass is REPORTED and gates nothing, as are the mean W1, " *
                               "the corner regions and the per-cell field sd.",
    )
    let d = JLD2.load(tmp)
        @assert haskey(d, "w1_per_region")
        @assert haskey(d, "w1_max")
        @assert haskey(d, "atom_mass_per_region")
        @assert haskey(d, "tolerance")
    end
    mv(tmp, SIM02_REPORT_PATH; force = true)

    println()
    println("="^100)
    println("GATED   max over all arms x rungs of the per-region W1 (scoped)     = ",
            maximum(w1_max), "   bar = ", P12_SIM02_W1_TOL_PERREGION)
    println("GATED   max over all arms x rungs of the per-region W1 (full range) = ",
            maximum(w1f_max))
    println("REPORTED, GATES NOTHING: per-region atom mass, mean over arms x rungs = ",
            mean(atom_mean), "  max = ", maximum(atom_max))
    println("wrote ", SIM02_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    println("="^100)
    return nothing
end

isdefined(@__MODULE__, :P12_SIM02_LOAD_ONLY) || main()
