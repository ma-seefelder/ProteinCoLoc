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

# spike/comparator/table.jl --- CMP-01/CMP-04/CMP-05: comparison table + divergence (D-09/D-10)
#
# The phase's positioning payload. Assembles one tidy `DataFrame` row per shared
# input carrying the classical estimator columns (Costes-p, Manders M1/M2, Pearson,
# Spearman), the simulator GROUND-TRUTH regime, optional NPE columns, and — the point
# of the whole harness — a `divergence` measure plus a green/amber/red `light` flag
# marking the rows where a classical verdict DISAGREES with the simulator ground truth
# ("here the classic is wrong, and v2.0 catches it"). Task 2 appends the atomic
# content-addressed CSV+JLD2 writer (D-05/D-07).
#
# ANTI-DATA-SNOOPING (D-14): `traffic_light` reads the pre-declared `DIVERGENCE_WARN`/
# `DIVERGENCE_FAIL` bands from comparator/config.jl — NEVER an inline literal. The
# divergence semantics are pinned DETERMINISTICALLY (not a loose mapping): a `:random`
# regime with a spuriously significant Pearson yields a LARGE divergence (classic
# falsely calls coloc), and a `:coloc` regime with a non-significant Pearson also
# yields a large divergence (classic misses real coloc).
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. The frozen src/ math is
# reached only transitively through the guarded includes below; no src/ is edited.

using DataFrames   # tidy per-input table assembly (D-09)

# Guarded includes so table.jl loads both standalone and after a sibling already
# pulled the config/estimators/inputs into scope (mirror inputs.jl:51-57).
isdefined(@__MODULE__, :DIVERGENCE_WARN)      || include(joinpath(@__DIR__, "config.jl"))
isdefined(@__MODULE__, :manders)              || include(joinpath(@__DIR__, "classical.jl"))
isdefined(@__MODULE__, :build_shared_inputs)  || include(joinpath(@__DIR__, "inputs.jl"))

# --- Traffic-light bands (D-10/D-14) --------------------------------------------

"""
    traffic_light(d) -> Symbol

Render a divergence value `d` as a green/amber/red traffic light using the
PRE-DECLARED comparator bands (`DIVERGENCE_WARN`/`DIVERGENCE_FAIL` from
comparator/config.jl — no inline literals, D-14). Monotone at the exact band
boundaries: `:green` for `d < WARN`, `:amber` for `WARN ≤ d < FAIL`, `:red` for
`d ≥ FAIL`. A `missing` divergence (e.g. an `:unknown` external input with no
trusted ground truth) is `:gray`.
"""
traffic_light(d) = ismissing(d) ? :gray :
                   d < DIVERGENCE_WARN ? :green :
                   d < DIVERGENCE_FAIL ? :amber : :red

# --- Divergence semantics (D-10) ------------------------------------------------

"""
    _regime_score(regime::Symbol) -> Union{Float64,Missing}

Map a simulator ground-truth regime to a scalar coloc score on `[-1, +1]`:
`:coloc → +1.0`, `:random → 0.0`, `:exclusion → -1.0`. An `:unknown` regime
(untrusted external input, no ground truth) returns `missing` so no divergence is
fabricated (T-09-07).
"""
_regime_score(regime::Symbol) =
    regime === :coloc     ? 1.0 :
    regime === :random    ? 0.0 :
    regime === :exclusion ? -1.0 :
    missing

"""
    classical_verdict(pearson, costes_p; alpha=COSTES_ALPHA) -> Float64

The classical estimator's coloc VERDICT on `[-1, +1]`: `clamp(pearson, -1, 1)`
when the correlation is Costes-significant (`costes_p < alpha`), else `0.0` — no
significant coloc detected. A non-finite `costes_p` (NaN from a degenerate input)
is treated as NON-significant (verdict `0.0`), never as a spurious call.
"""
function classical_verdict(pearson, costes_p; alpha = COSTES_ALPHA)
    significant = isfinite(costes_p) && costes_p < alpha
    return significant ? clamp(pearson, -1.0, 1.0) : 0.0
end

"""
    divergence(regime, verdict) -> Union{Float64,Missing}

The positioning payload: `abs(_regime_score(regime) - verdict)`, the absolute gap
between what the simulator KNOWS is true and what the classic CLAIMS. A `:random`
regime (score 0) with a significant high Pearson (verdict ≈ 0.8) gives a LARGE
divergence (classic falsely calls coloc); a `:coloc` regime (score +1) with a
non-significant Pearson (verdict 0) also gives a large divergence (classic misses
real coloc). Returns `missing` for `:unknown` (no trusted ground truth).
"""
function divergence(regime, verdict)
    rs = _regime_score(regime)
    return ismissing(rs) ? missing : abs(rs - verdict)
end

# --- Tidy per-input table (D-09) ------------------------------------------------

"""
    build_table(inputs; master_seed=inputs.master_seed, npe=nothing) -> DataFrame

Assemble the tidy comparison table (D-09): one row per shared input in `inputs`
(the `build_shared_inputs` NamedTuple), with columns

    input_id, regime, rho_true,
    Costes_p, M1, M2, Pearson, Spearman,
    NPE_rho_hat, NPE_OOD_flag,
    divergence, light

Each classical column is computed from the estimator battery in classical.jl
(`costes_p(master_seed, idx, mci)`, `manders`, `pearson_whole`, `spearman_whole`).
The `divergence`/`light` columns are the pre-declared-threshold positioning payload
(D-10/D-14). NPE columns are `missing` unless an `npe` NamedTuple with `rho_hat`
and `ood_flag` vectors is supplied (D-09 optional; OQ3). Estimators are already
NaN-safe, so no row throws on degenerate input.
"""
function build_table(inputs; master_seed = inputs.master_seed, npe = nothing)
    images   = inputs.images
    regime   = inputs.regime
    rho_true = inputs.rho_true
    n        = length(images)

    input_id = collect(1:n)
    Costes_p = Vector{Float64}(undef, n)
    M1v      = Vector{Float64}(undef, n)
    M2v      = Vector{Float64}(undef, n)
    Pearson  = Vector{Float64}(undef, n)
    Spearman = Vector{Float64}(undef, n)
    divc     = Vector{Union{Missing,Float64}}(undef, n)
    lightc   = Vector{Symbol}(undef, n)

    for i in 1:n
        mci    = images[i]
        cp     = costes_p(master_seed, i, mci)
        m1, m2 = manders(mci)
        pe     = pearson_whole(mci)
        sp     = spearman_whole(mci)

        Costes_p[i] = cp
        M1v[i]      = m1
        M2v[i]      = m2
        Pearson[i]  = pe
        Spearman[i] = sp

        verdict   = classical_verdict(pe, cp)
        d         = divergence(regime[i], verdict)
        divc[i]   = d
        lightc[i] = traffic_light(d)
    end

    # NPE columns: missing unless an NPE artifact is supplied (D-09 optional).
    if npe === nothing
        NPE_rho_hat  = fill(missing, n)
        NPE_OOD_flag = fill(missing, n)
    else
        NPE_rho_hat  = collect(npe.rho_hat)
        NPE_OOD_flag = collect(npe.ood_flag)
    end

    return DataFrame(
        input_id     = input_id,
        regime       = collect(regime),
        rho_true     = collect(rho_true),
        Costes_p     = Costes_p,
        M1           = M1v,
        M2           = M2v,
        Pearson      = Pearson,
        Spearman     = Spearman,
        NPE_rho_hat  = NPE_rho_hat,
        NPE_OOD_flag = NPE_OOD_flag,
        divergence   = divc,
        light        = lightc,
    )
end
