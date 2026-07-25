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

# spike/test/test_lambda_ablation.jl --- SC1g: the R4 lambda-conditioning tripwire. BLOCKS Wave 4.
#
# WHAT A FAILURE HERE MEANS, IN PLAIN TERMS: THE LAMBDA CONDITIONING IS DEAD. The 129th input
# row carries no information about the registration uncertainty the user declared, the net has
# correctly learned to ignore it, and EVERY SUBSEQUENT SC2 NUMBER WOULD BE A NULL THAT LOOKS
# EXACTLY LIKE "the effect is below the fixed 8x8 summary's resolution" WHILE ACTUALLY BEING AN
# IMPLEMENTATION BUG. Those two readings have completely different consequences -- one is an
# honest scientific limit that gets reported, the other is a defect that gets fixed -- and this
# file is the ONLY thing in the phase that tells them apart.
#
# The usual cause is R4 / Pitfall 1: the training sampler drew lambda INDEPENDENTLY of the shift
# instead of drawing lambda FIRST and then shift ~ Uniform(-lambda, lambda). If this file fails,
# inspect the datagen sampler before touching anything else.
#
# WHEN TO RUN IT: immediately after the research net is trained (plan 11-07) and BEFORE any
# ladder run (plan 11-08). Run:
#     julia --project=spike spike/test/test_lambda_ablation.jl
#
# IT REFUSES TO PASS VACUOUSLY. If the trained research net is absent this file prints a clear
# message and EXITS NON-ZERO. A tripwire that reports success because it had nothing to check is
# precisely the failure mode it exists to prevent, so there is no skip path here by design.
#
# CPU-only (D-10). Reads the net and the simulator; writes nothing.

using Test
using Statistics
using StatsBase
using Random
using Random123
using Pkg

# ORDER MATTERS: the Tier-1 (and, once it exists, Tier-2) pre-registration first, then the
# research model surface (augment_input / encode_lambda), then the read surface (load_npe,
# rho_draws, interval_width, standardize_summary), then the forward model and the FROZEN summary
# contract. Guarded for idempotency.
#
# `spike/validation/consts.jl` and `harness.jl` are DELIBERATELY NOT included: they redeclare
# VAL_MASTER_SEED / VAL_FIX_SEED at a different integer width than `p11_consts.jl` binds them,
# and a `const` redeclaration at a different type is a hard error. The three simulator/contract
# includes below are the same ones `spike/validation/harness.jl:55-57` uses.
isdefined(@__MODULE__, :LAMBDA_MIN)     || include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :encode_lambda)  || include(joinpath(@__DIR__, "..", "npe", "p11_architecture.jl"))
isdefined(@__MODULE__, :posterior_for)  || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))
isdefined(@__MODULE__, :simulate_pair)  || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)      || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)     || include(joinpath(@__DIR__, "..", "data", "encode.jl"))

# The trained Phase-11 research net. Produced by plan 11-07; it does not exist before then.
const P11_NET_PATH = joinpath(@__DIR__, "..", "npe", "p11_research_npe.jld2")

if !isfile(P11_NET_PATH)
    println("="^78)
    println("SC1g lambda-ablation tripwire — CANNOT RUN")
    println("="^78)
    println("research net not trained yet — this tripwire cannot pass vacuously.")
    println()
    println("  expected trained research net at:")
    println("      $(P11_NET_PATH)")
    println()
    println("  That file is produced by plan 11-07 (the Phase-11 research-net training run).")
    println("  Until it exists there is nothing to ablate, and reporting success would be a")
    println("  FALSE PASS on the one check that distinguishes 'below resolution' from 'the")
    println("  lambda conditioning is dead'. Exiting non-zero on purpose.")
    println("="^78)
    exit(1)
end

# --- Resolution knobs. NOT pre-registered thresholds -- Monte-Carlo resolution only. The one
#     pass/fail threshold in this file is the ablation FACTOR below.
const ABLATION_N_DATASETS = 3      # independent (sample, control) pairs; ALL must clear the bar
const ABLATION_N_DRAWS    = 2000   # posterior draws per arm
const ABLATION_IMSIZE     = P11_PROBE_COMPARABILITY_IMSIZE  # a fixed arm; this file makes no
                                                            # F5 train-joint / eval-joint claim

# THE MATERIAL-INEQUALITY THRESHOLD. This is a tripwire, so the bar must be a REAL EFFECT SIZE,
# not "the two numbers differ" — floating-point inequality would pass on Monte-Carlo noise alone
# and the tripwire would be decorative.
#
# `P11_LAMBDA_ABLATION_FACTOR` belongs in the TIER-2 block of `spike/validation/p11_consts.jl`,
# derived from the D-06 probe's measured effect scale. Tier 2 does not exist at authoring time,
# so it is read defensively with a documented fallback.
#
# 1.15 IS A PLACEHOLDER THAT PLAN 11-05 MUST SUPERSEDE. It is deliberately conservative: the
# ladder is designed to move the posterior width by a factor of order 3.7x across
# [LAMBDA_MIN, LAMBDA_MAX] (p11_consts.jl §5), so a 15% floor is far below the effect being
# claimed and above Monte-Carlo noise at 2000 draws. It is a "the conditioning is alive at all"
# bar, NOT the SC2 criterion. Once Tier 2 lands, the `isdefined` branch below picks it up with
# no edit to this file.
const P11_LAMBDA_ABLATION_FACTOR_FALLBACK = 1.15
const ABLATION_FACTOR = isdefined(@__MODULE__, :P11_LAMBDA_ABLATION_FACTOR) ?
                        P11_LAMBDA_ABLATION_FACTOR : P11_LAMBDA_ABLATION_FACTOR_FALLBACK
const ABLATION_FACTOR_IS_TIER2 = isdefined(@__MODULE__, :P11_LAMBDA_ABLATION_FACTOR)

# GLOBAL-RNG PIN (the `test/gate/harness.jl:85-96` `seed_gate_global!` discipline).
# `sampleposterior` threads NO `rng` — it draws the flow's base samples from the GLOBAL stream —
# so without a pin this tripwire cannot reproduce its own numbers, and a marginal failure could
# not be told from a lucky draw. The seed is DERIVED from the Phase-11 DEV stream on the FIXTURE
# counter (never a reported counter), so no new seed constant is introduced.
const ABLATION_GLOBAL_SEED = rand(p11_rng(P11_FIXTURE_COUNTER), UInt64)

const ABLATION_MODEL = load_npe(P11_NET_PATH)

"""
    _ablation_pair(rng) -> (Zs, Zc, θs, θc)

Simulate ONE (sample, control) dataset and push both stacks through the FROZEN path
`standardize_summary(encode_d01(patch_summary(build_mci(...))), zt, :min)` — the same chain
`spike/validation/harness.jl:126-131` uses. The two 128-row summaries are computed ONCE and
reused by every lambda arm, so the arms differ ONLY in the 129th row and there is ZERO
simulation noise between them. That is what makes this the controlled contrast D-03 promises.
"""
function _ablation_pair(rng)
    θs = sample_prior(rng)
    Zs = standardize_summary(encode_d01(patch_summary(build_mci(
             simulate_pair(rng, θs; imsize = ABLATION_IMSIZE)))), ABLATION_MODEL.zt, :min)
    θc = sample_prior(rng)
    Zc = standardize_summary(encode_d01(patch_summary(build_mci(
             simulate_pair(rng, θc; imsize = ABLATION_IMSIZE)))), ABLATION_MODEL.zt, :min)
    return (Zs = Zs, Zc = Zc, θs = θs, θc = θc)
end

"""
    _delta_rho_draws(Zs, Zc, lam; N = ABLATION_N_DRAWS) -> Vector

The Δρ posterior draw vector read at registration-uncertainty level `lam` from ALREADY-computed
128-row summaries: append the lambda row with `augment_input`, take one `sampleposterior` pass
per stack, and form the Monte-Carlo difference (the `spike/npe/infer.jl:134-138` `delta_rho`
construction, kept as DRAWS rather than collapsed to a mean so a spread statistic is available).

The global RNG is re-pinned to the SAME seed at the top of every call, so the flow's base sample
stream is IDENTICAL across lambda arms. Any change in spread is therefore attributable to the
conditioning input and to nothing else.
"""
function _delta_rho_draws(Zs, Zc, lam; N::Integer = ABLATION_N_DRAWS)
    Random.seed!(ABLATION_GLOBAL_SEED)
    ρs = rho_draws(ABLATION_MODEL.estimator, augment_input(Zs, lam), ABLATION_MODEL.θzt;
                   N = N, use_gpu = false)
    ρc = rho_draws(ABLATION_MODEL.estimator, augment_input(Zc, lam), ABLATION_MODEL.θzt;
                   N = N, use_gpu = false)
    return ρs .- ρc
end

"90% credible-interval width of a Δρ draw vector, via the sanctioned `interval_width` helper."
_delta_rho_hdi(Δ) = interval_width(reshape(collect(Δ), 1, :))[1]

# --- Measure every dataset BEFORE asserting anything, so a failure is diagnosable from the log
#     alone rather than from a bisect.
const ABLATION_RNG = p11_rng(P11_FIXTURE_COUNTER)
const ABLATION_ROWS = map(1:ABLATION_N_DATASETS) do i
    p     = _ablation_pair(ABLATION_RNG)
    Δ_lo  = _delta_rho_draws(p.Zs, p.Zc, LAMBDA_MIN)
    Δ_hi  = _delta_rho_draws(p.Zs, p.Zc, LAMBDA_MAX)
    sd_lo, sd_hi = std(Δ_lo), std(Δ_hi)
    w_lo,  w_hi  = _delta_rho_hdi(Δ_lo), _delta_rho_hdi(Δ_hi)
    (i = i, sd_lo = sd_lo, sd_hi = sd_hi, sd_ratio = sd_hi / sd_lo,
     w_lo = w_lo, w_hi = w_hi, w_ratio = w_hi / w_lo, pair = p)
end

_f(x) = rpad(string(round(x; digits = 5)), 12)

println("="^78)
println("SC1g — lambda-conditioning tripwire (R4). BLOCKS the SC2 ladder.")
println("  net            : $(P11_NET_PATH)")
println("  lambda arms    : LAMBDA_MIN = $LAMBDA_MIN  vs  LAMBDA_MAX = $LAMBDA_MAX")
println("  required factor: $(ABLATION_FACTOR) " *
        (ABLATION_FACTOR_IS_TIER2 ? "(Tier-2, pre-registered)" :
                                    "(PLACEHOLDER — plan 11-05 must supersede with Tier 2)"))
println("  draws/arm      : $ABLATION_N_DRAWS   datasets: $ABLATION_N_DATASETS   " *
        "global seed: $(repr(ABLATION_GLOBAL_SEED))")
println("-"^78)
println(rpad("dataset", 9), _f("sd(Δρ|min)"), _f("sd(Δρ|max)"), _f("sd ratio"),
        _f("hdi(min)"), _f("hdi(max)"), _f("hdi ratio"))
for r in ABLATION_ROWS
    println(rpad(r.i, 9), _f(r.sd_lo), _f(r.sd_hi), _f(r.sd_ratio),
            _f(r.w_lo), _f(r.w_hi), _f(r.w_ratio))
end
println("="^78)

@testset "SC1g — lambda conditioning is alive (R4 tripwire)" verbose = true begin

    @testset "posterior SD widens materially from LAMBDA_MIN to LAMBDA_MAX" begin
        # A MATERIAL inequality against a pre-registered factor — never a bare difference of
        # floats, which would pass on Monte-Carlo noise and make the tripwire decorative.
        for r in ABLATION_ROWS
            @test r.sd_hi > ABLATION_FACTOR * r.sd_lo
        end
        # EVERY dataset must clear the bar, so one lucky draw cannot carry the verdict.
        @test all(r.sd_hi > ABLATION_FACTOR * r.sd_lo for r in ABLATION_ROWS)
    end

    @testset "90% HDI width widens materially too (second spread statistic)" begin
        # Repeating the claim on `interval_width` means the verdict does not hinge on a single
        # spread statistic — an SD-only response could be an artifact of a heavy tail.
        for r in ABLATION_ROWS
            @test r.w_hi > ABLATION_FACTOR * r.w_lo
        end
        @test all(r.w_hi > ABLATION_FACTOR * r.w_lo for r in ABLATION_ROWS)
    end

    @testset "the statistic can FAIL (deliberately dead conditioning)" begin
        # The `test_sbc.jl:63-68` both-directions discipline. Feed the SAME lambda to both arms:
        # the conditioning input is then constant by construction, i.e. DEAD, and the inequality
        # MUST NOT hold. If this passes, the statistic is stuck at "pass" and every assertion
        # above is worthless.
        p = first(ABLATION_ROWS).pair
        dead_a = _delta_rho_draws(p.Zs, p.Zc, LAMBDA_MIN)
        dead_b = _delta_rho_draws(p.Zs, p.Zc, LAMBDA_MIN)
        @test std(dead_b) <= ABLATION_FACTOR * std(dead_a)
        @test _delta_rho_hdi(dead_b) <= ABLATION_FACTOR * _delta_rho_hdi(dead_a)
        # Same lambda + same pinned global seed => byte-identical draws, which is also the
        # reproducibility proof for the global-RNG pin itself.
        @test dead_a == dead_b
    end

    @testset "the ablation factor is sourced, and its provenance is visible" begin
        @test ABLATION_FACTOR > 1.0
        @test isfinite(ABLATION_FACTOR)
        # Until Tier 2 lands the placeholder is in force and says so in the printed banner.
        @test ABLATION_FACTOR_IS_TIER2 || ABLATION_FACTOR == P11_LAMBDA_ABLATION_FACTOR_FALLBACK
    end

    @testset "tripwire ran CPU-only (D-10)" begin
        @test !haskey(Pkg.project().dependencies, "CUDA")
        @test !any(p -> occursin("CUDA", p.name), values(Pkg.dependencies()))
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
