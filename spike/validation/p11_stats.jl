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

# spike/validation/p11_stats.jl --- Phase-11 equivalence / monotonicity statistics (D-07, D-08).
#
# The pure statistical layer the SC2 ladder and the SC3 breakdown curve score against:
# a Wilson score interval, the two one-sided equivalence tests, a Holm step-down whose
# comparison direction is the MIRROR IMAGE of the SBC one, a rung-label permutation
# Spearman, the vacuity (shrinkage) diagnostic, and the breakdown-point definition.
# Nothing here runs inference, simulates, or writes a file.
#
# DECOUPLING (CLAUDE.md): spike-local; reaches `src/` only read-only, and `test/gate/`
# only as COPIED-WITH-ATTRIBUTION source (see the two copied functions below) because
# `spike/Project.toml` cannot load the gate file's dependency chain.

# Unit dependency: the locked thresholds. Guarded for idempotency (S2).
isdefined(@__MODULE__, :P11_DEV_SEED) ||
    include(joinpath(@__DIR__, "p11_consts.jl"))

using Statistics                     # mean, std
using StatsBase                      # corspearman
using Random                         # AbstractRNG, shuffle!
import Distributions: Normal, cdf    # the standard-normal CDF; erf is NOT hand-rolled

# =========================================================================================
# Coverage interval and the two one-sided tests (D-08)
# =========================================================================================

"""
    wilson_ci(k, n; z = Z_TWO_SIDED_90) -> (lo, hi)

Wilson score interval for a Binomial proportion `k/n`, at the two-sided level implied by
`z` (`Z_TWO_SIDED_90` gives the 90 % limb, which is the `1 − 2α` interval matching a TOST
at `SC2_TOST_ALPHA` per side).

WALD IS DELIBERATELY NOT USED. At `p ≈ SC2_COVERAGE_NOMINAL = 0.90` with `n` in the low
hundreds, the Wald interval's actual coverage is poor near the boundary and it can extend
past 1; Wilson is the standard remedy and is what the pre-registered `N_MIN_DERIVED`
arithmetic assumes.
"""
function wilson_ci(k::Integer, n::Integer; z::Real = Z_TWO_SIDED_90)
    p̂ = k / n
    d = 1 + z^2 / n
    c = (p̂ + z^2 / (2n)) / d
    h = (z / d) * sqrt(p̂ * (1 - p̂) / n + z^2 / (4n^2))
    return (lo = c - h, hi = c + h)
end

"""
    tost_pvalue(k, n; p0 = SC2_COVERAGE_NOMINAL, delta = SC2_TOST_DELTA) -> Float64

The two-one-sided-tests p-value for equivalence of the empirical coverage `k/n` to `p0`
within the pre-registered band `±delta`:

    H₀₁ : p ≤ p0 − delta   vs   H₁₁ : p > p0 − delta
    H₀₂ : p ≥ p0 + delta   vs   H₁₂ : p < p0 + delta

Equivalence is declared iff BOTH are rejected, so the TOST p-value is `max(p₁, p₂)` and a
SMALL value is the equivalent (good) outcome. Normal-approximation one-sided tests, using
`cdf(Normal(), ·)` rather than a hand-rolled `erf`.

The standard error is guarded against a degenerate `p̂ ∈ {0, 1}` (which would make `se = 0`
and the ratio non-finite) by flooring the variance at `1e-12`.
"""
function tost_pvalue(k::Integer, n::Integer;
                     p0::Real = SC2_COVERAGE_NOMINAL, delta::Real = SC2_TOST_DELTA)
    p̂  = k / n
    se = sqrt(max(p̂ * (1 - p̂), 1e-12) / n)
    p1 = 1 - cdf(Normal(), (p̂ - (p0 - delta)) / se)   # H₀₁: p ≤ p0 − delta
    p2 =     cdf(Normal(), (p̂ - (p0 + delta)) / se)   # H₀₂: p ≥ p0 + delta
    return max(p1, p2)
end

# =========================================================================================
# Holm step-down --- COPIED, and then INVERTED (D-08)
# =========================================================================================

"""
    p11_holm_adjusted(p) -> Vector{Float64}

Holm–Bonferroni adjusted p-values for a family of `m = length(p)` tests, in INPUT order:
sort ascending, multiply p₍ⱼ₎ by `(m − j + 1)`, take the running maximum (monotonicity) and
cap at 1. Valid under arbitrary dependence between the tests.

COPIED VERBATIM from `test/gate/sbc.jl:352-362` (identical to the frozen
`holm_adjusted` in `test/gate/gate_consts_8_v2.jl:361-371`, which that function defers to).
It is COPIED rather than depended on because `test/gate/sbc.jl` transitively includes
`test/gate/harness.jl:35`, which does `using ProteinCoLoc` — a package that cannot resolve
under `--project=spike` (MEASURED, not assumed). Re-implementing Holm from scratch is
forbidden by the "don't hand-roll" rule, so this is a copy with attribution; a unit test
asserts it agrees with the frozen adjuster.
"""
function p11_holm_adjusted(p::AbstractVector{<:Real})
    m   = length(p)
    ord = sortperm(collect(float.(p)))
    adj = Vector{Float64}(undef, m)
    running = 0.0
    for (j, i) in enumerate(ord)
        running = max(running, (m - j + 1) * float(p[i]))
        adj[i]  = min(1.0, running)
    end
    return adj
end

# -----------------------------------------------------------------------------------------
# ***  THE DIRECTION TRAP.  READ BEFORE TOUCHING THE LINE BELOW.  ***
#
# This is the MIRROR IMAGE of `sbc_holm_pass` (`test/gate/sbc.jl:374`), which is
#
#     sbc_holm_pass(p; fwer) = all(>(fwer), sbc_holm(p))     # pass iff Holm rejects NOTHING
#
# and is CORRECT THERE, because the SBC null is UNIFORMITY and a rejection is BAD news.
# The Phase-11 family's null is NON-EQUIVALENCE, so a rejection is GOOD news: the phase must
# reject every rung's non-equivalence null. Hence `all(<=(fwer), ...)`, NOT `all(>(fwer), ...)`.
# COPYING `holm_pass` VERBATIM WOULD COMPUTE THE EXACT OPPOSITE SC2(b) VERDICT, silently, with
# no error and a plausible-looking number. For that reason the name `holm_pass` is NEVER reused
# in this file, and `spike/test/test_p11_tost.jl` pins the direction with a hand-computed vector.
#
# HONESTY NOTE for the report (this does not contradict D-08, it strengthens it). "All rungs
# equivalent" is an INTERSECTION–UNION claim: requiring every component TOST to reject at level
# α already has family-wise size ≤ α, so no multiplicity correction is REQUIRED and Holm makes
# the test strictly harder. D-08's choice is therefore deliberately conservative in the right
# direction — the same direction `gate_consts_8_v2.jl:213-218` argues for on the BF arm — and
# must be described as deliberately conservative rather than as necessary.
# -----------------------------------------------------------------------------------------

"""
    p11_tost_pass(p_tost; fwer = P11_HOLM_FWER) -> Bool

`true` iff Holm–Bonferroni at family-wise error rate `fwer` REJECTS EVERY non-equivalence
null in the family — i.e. every adjusted TOST p-value is at most `fwer`. Note `<=`, not `>`:
see the direction-trap block immediately above.
"""
p11_tost_pass(p_tost::AbstractVector{<:Real}; fwer::Real = P11_HOLM_FWER) =
    all(<=(fwer), p11_holm_adjusted(p_tost))

# =========================================================================================
# SC2(a): monotone widening, tested by a rung-label permutation (D-07a)
# =========================================================================================

"""
    p11_perm_spearman(lam, w; B = PERM_B, rng) -> (S, p_perm)

One-sided permutation test of "λ has no effect on posterior width". `S` is
`corspearman(lam, w)` over the pooled `R × N_PER_RUNG` pairs; `p_perm` is
`(1 + #{S_b ≥ S}) / (B + 1)` over `B` shuffles.

THE PERMUTATION IS OVER THE RUNG LABELS, NOT THE WIDTHS. The widths are NOT independent
within a rung (they share a rung's λ and its draw pool), which is exactly the dependence an
asymptotic Spearman p-value would misstate. Shuffling the labels makes the null "λ carries
no information about width" exact under that within-rung structure.

`rng` is a REQUIRED keyword and is threaded explicitly — `Random.seed!` is never called in
a loop anywhere in this codebase.
"""
function p11_perm_spearman(lam::AbstractVector, w::AbstractVector;
                           B::Integer = PERM_B, rng::AbstractRNG)
    length(lam) == length(w) ||
        throw(ArgumentError("p11_perm_spearman: lam and w must be the same length, " *
                            "got $(length(lam)) and $(length(w))"))
    S = corspearman(collect(float.(lam)), collect(float.(w)))
    perm  = collect(float.(lam))
    ge    = 0
    for _ in 1:B
        shuffle!(rng, perm)
        corspearman(perm, collect(float.(w))) >= S && (ge += 1)
    end
    return (S = S, p_perm = (1 + ge) / (B + 1))
end

# =========================================================================================
# Vacuity diagnostic (D-05, F3) --- REPORTING ONLY, never a pass/fail input
# =========================================================================================

# The vacuity cutoff is a REPORTING LABEL, so it lives here and deliberately NOT in the
# pre-registered `p11_consts.jl` — exactly the reasoning `test/gate/sbc.jl:226-229` gives for
# keeping `SBC_VACUOUS_SHRINKAGE` out of the frozen consts.
const P11_VACUOUS_SHRINKAGE = 0.95

"""
    p11_shrinkage(post_sd, prior_draws; cutoff = P11_VACUOUS_SHRINKAGE)
        -> (shrinkage, post_sd_mean, prior_sd, vacuous)

Per-parameter posterior-to-prior spread ratio. `prior_sd` is the standard deviation of the
θ* values ACTUALLY drawn (an empirical prior sample — no analytic form is assumed, which
matters because `ρ_true = ghat(μ*)` has no closed form); `post_sd_mean` is the mean posterior
sd; `shrinkage` is their ratio. `shrinkage ≈ 1` ⇒ the posterior retains essentially all the
prior spread ⇒ that column's rank statistics are VACUOUS. Diagnostic only.

COPIED from `test/gate/sbc.jl:243-250` (already generic in the column count, so the 7 → 8 θ
growth needs no change), for the same load-path reason given on `p11_holm_adjusted`.

*** PITFALL: PASS ONE RUNG'S PRIOR DRAWS, NOT THE POOL. *** Under the λ-conditioned
hierarchical sampler the shift columns' prior SD WITHIN a rung is `λ_r/√3`, not the pooled
marginal. Pooling across rungs inflates `prior_sd`, deflates `shrinkage`, and would make the
shift columns look INFORMATIVE — contradicting D-05, which says they are marginalized, not
identified. Compute shrinkage PER RUNG. (`chromatic_eps`'s prior is λ-independent, so its
`prior_sd` is the same pooled or per-rung.)
"""
function p11_shrinkage(post_sd::AbstractMatrix, prior_draws::AbstractMatrix;
                       cutoff::Real = P11_VACUOUS_SHRINKAGE)
    n = size(post_sd, 2)
    psd_mean = [Statistics.mean(view(post_sd, :, p)) for p in 1:n]
    prior_sd = [Statistics.std(view(prior_draws, :, p)) for p in 1:n]
    shrink   = [prior_sd[p] > 0 ? psd_mean[p] / prior_sd[p] : NaN for p in 1:n]
    vacuous  = [isfinite(s) && s >= cutoff for s in shrink]
    return (shrinkage = shrink, post_sd_mean = psd_mean, prior_sd = prior_sd, vacuous = vacuous)
end

# =========================================================================================
# SC3 breakdown point (D-14) --- REPORTED, never gated
# =========================================================================================

"""
    breakdown_point(ks, n, rungs) -> eltype(rungs) or nothing

The FIRST rung magnitude whose empirical coverage has evidentially dropped below nominal:
the smallest `rungs[i]` with `wilson_ci(ks[i], n).hi < SC2_COVERAGE_NOMINAL`. Returns
`nothing` for the right-censored case (no rung on the ladder satisfies it), which is
reported as "coverage holds past the end of the ladder".

WHY THE UPPER LIMB AND NOT THE POINT ESTIMATE: it makes "coverage has genuinely dropped" an
evidentiary claim rather than a coin flip on Monte-Carlo noise, matching the burden-of-proof
logic `gate_consts_8_v2.jl:211-218` already uses for the BF arm.

WHY THE FIRST AND NOT THE LAST: monotone degradation is expected but not guaranteed, and
"first" is unambiguous under non-monotonicity. The full curve is reported regardless.
"""
function breakdown_point(ks::AbstractVector{<:Integer}, n::Integer, rungs)
    length(ks) == length(rungs) ||
        throw(ArgumentError("breakdown_point: ks and rungs must be the same length, " *
                            "got $(length(ks)) and $(length(rungs))"))
    for (k, x) in zip(ks, rungs)
        wilson_ci(k, n).hi < SC2_COVERAGE_NOMINAL && return x
    end
    return nothing
end
