#########################################################################################
# Spike 012 --- DEV-seed side-by-side evaluation: HIGH-CAP flow vs amended_v2.
#
# QUESTION. Does the raised-capacity flow reduce the residual F2 marginal LOCATION drift on the
# non-identified nuisances (`autofluorescence`, `label_efficiency`) below the 0.10-SD equivalence
# margin (07-NUISANCE-SBC-SPEC §3), while keeping the coloc TARGETS (ρ_true, Δρ) clean?
#
# METHOD. BOTH nets are evaluated on the SAME fresh DEV seed over the SAME F5 image-size mixture.
# The prior draws + forward simulations come from a DEV-keyed Philox stream that does NOT depend on
# the net, so θ*/Z are BYTE-IDENTICAL across the two nets — a paired comparison of the two
# posteriors on identical simulated datasets. The flow's base draws come from the GLOBAL RNG
# (sampleposterior threads no rng); it is re-seeded to the SAME DEV value before each net, so the
# base randomness is matched too (the two nets still consume the global stream at different rates
# because their flows differ in size — unavoidable and immaterial to θ*/Z).
#
# For each of the 8 SBC columns (7 θ + the paired-draw Δρ, D-01) and each net:
#   • SBC mean-rank z         : (mean(u) − 0.5) / sqrt((1/12)/M),  u = (rank+0.5)/(L+1)
#   • marginal drift (§3)     : |mean(u) − 0.5| · sqrt(12)   in posterior-SD units (the spec's
#                               equivalence quantity; margin δ = 0.10 SD)
#   • standardized bias b     : mean((E[θ̂] − θ*)/post_sd)   (the task's requested quantity;
#                               = spike-011's "posterior-SD drift"; opposite sign to s by construction)
#   • shrinkage               : mean(post_sd) / std(θ*)      (≈1 ⇒ non-identified)
#   • KS uniformity p         : ExactOneSampleKSTest(u, Uniform(0,1))
# ρ_true carries PRIOR ATOMS (ghat(μ*) clamps at ±0.99; Spike 010 addendum) — SBC ranks are
# undefined on atoms, so for ρ_true we ALSO report atom-excluded KS and randomized-rank KS (the
# standard mixed-distribution correction). Applied ONLY to ρ_true (spec §2: the only atom carrier).
#
# DIAGNOSTIC ONLY. DEV seed asserted disjoint from every PROD/VAL/NPE/DEFAULT seed and every dev
# seed burned so far. NO gate report is written; NO PROD seed is consumed.
#
# Run:  julia --project -t auto .planning/spikes/012-capacity-lever/eval_highcap_vs_amended.jl
#########################################################################################
using Random, Random123, Statistics, Printf, JLD2
using HypothesisTests, Distributions
using ProteinCoLoc

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const GATE = joinpath(REPO, "test", "gate")

# Load the frozen amended pre-registration (imsize mixture + PROD_SEED_V2/PROD_SEED for the
# forbidden set), then the gate harness (draw_simulate_infer[_paired]). Exactly the mixture_replicate
# pattern from Spike 010c.
include(joinpath(GATE, "gate_consts_8_v2.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(GATE, "harness.jl"))

# ---- fresh DEV seed, asserted disjoint from EVERYTHING pre-registered / already burned ----------
const DEV = UInt64(0x00CA11B3)          # next in the 0x00CA11Bx dev series; fresh
const NAMED_FORBIDDEN = UInt64[
    0x5BC0FFEE,      # VAL_MASTER_SEED
    0xC0FFEE,        # NPE_MASTER_SEED
    0x1,             # DEFAULT_MASTER_SEED
    0xDE7C0DE, 0xDE7C0DE2, 0x0DE7C0D3,               # burned dev seeds
    0x00B7A771, 0x00B7A772, 0x00CA11B0, 0x00CA11B1,  # burned dev seeds
    0x00A1CE01, 0x00CA11B2,                          # burned dev seeds
]
const FORBIDDEN = Set{UInt64}(vcat(NAMED_FORBIDDEN,
                                   collect(values(PROD_SEED_V2)),   # every amended gate seed
                                   collect(values(PROD_SEED))))     # every v1 gate seed
@assert !(DEV in FORBIDDEN) "DEV seed $(repr(DEV)) collides with a forbidden seed"
@assert !(DEV in values(PROD_SEED_V2)) "DEV collides with a PROD_SEED_V2 value"
println("DEV seed = ", repr(DEV), "  (asserted disjoint from ", length(FORBIDDEN), " forbidden seeds)")

const M = parse(Int, get(ENV, "SPIKE012_M", "1000"))
const L = SBC_L                         # 999
const LAB = ("ρ_true", "spillover", "autofluorescence", "label_efficiency",
             "shift_dx", "shift_dy", "noise", "Δρ")

const NETS = (
    (name = "high-cap (spike012)", path = joinpath(REPO, "artifacts", "spike012_highcap", "grid_8", "npe_8.jld2")),
    (name = "amended_v2",          path = joinpath(REPO, "artifacts", "amended_v2",        "grid_8", "npe_8.jld2")),
)

sim = default_simulator()

"""
Run the M×L draw→simulate→infer chain for one net `m` on the DEV stream, returning the 8-column
(7 θ + Δρ) rank / posterior-mean / posterior-sd / θ* tables. The DEV `rng` (prior + simulation)
is net-independent, so θ*/Z are identical across nets.
"""
function run_eval(m)
    rng = Random123.Philox4x(UInt64, (DEV, UInt64(8)))
    Random.seed!(DEV)                    # GLOBAL stream (sampleposterior base draws)
    ranks = zeros(Int, M, 8); pmean = zeros(M, 8); psd = zeros(M, 8); θ = zeros(M, 8)
    szs = Vector{Tuple{Int,Int}}(undef, M)
    for i in 1:M
        d = draw_simulate_infer(m, rng; G = 8, imsize = :mixture,
                                imsize_set = SBC_IMSIZE_SET, imsize_weights = SBC_IMSIZE_WEIGHTS,
                                N = L, sim = sim)
        θv = (d.θ.ρ_true, d.θ.spillover, d.θ.autofluorescence, d.θ.label_efficiency,
              d.θ.shift_dx, d.θ.shift_dy, d.θ.noise)
        for p in 1:7
            ranks[i, p] = count(<(θv[p]), view(d.draws, p, :))
            θ[i, p]     = θv[p]
            pmean[i, p] = mean(view(d.draws, p, :))
            psd[i, p]   = std(view(d.draws, p, :))
        end
        szs[i] = d.imsize
        # Δρ paired draw (D-01), independent of the marginal ρ_true rank.
        pr = draw_simulate_infer_paired(m, rng; G = 8, imsize = :mixture,
                                        imsize_set = SBC_IMSIZE_SET,
                                        imsize_weights = SBC_IMSIZE_WEIGHTS, N = L, sim = sim)
        Δd = pr.ρs .- pr.ρc
        Δs = pr.θs.ρ_true - pr.θc.ρ_true
        ranks[i, 8] = count(<(Δs), Δd)
        θ[i, 8]     = Δs
        pmean[i, 8] = mean(Δd)
        psd[i, 8]   = std(Δd)
        i % 100 == 0 && @printf("    %4d/%d\n", i, M)
    end
    return (; ranks, pmean, psd, θ, szs)
end

const se_u = sqrt((1 / 12) / M)
ks_p(u) = pvalue(ExactOneSampleKSTest(u, Uniform(0.0, 1.0)))

results = Dict{String,Any}()
for net in NETS
    println("\n==================================================================")
    if !isfile(net.path)
        println("MISSING: ", net.path, "  (skipping ", net.name, ")")
        continue
    end
    println("Evaluating ", net.name, "  M=", M, " L=", L)
    m = ProteinCoLoc.load_estimator(net.path)
    @printf("  flow capacity: coupling=%d flow_depth=%d flow_width=%d  (summary dstar=%d depth=%d width=%d)\n",
            m.arch.num_coupling_layers, m.arch.flow_depth, m.arch.flow_width,
            m.arch.dstar, m.arch.depth, m.arch.width)
    r = run_eval(m)
    results[net.name] = r

    @printf("\n  %-18s %8s %8s %9s %10s %10s %12s\n",
            "parameter", "mean(u)", "z", "drift_SD", "std_bias", "shrink", "KS p")
    for p in 1:8
        u  = (r.ranks[:, p] .+ 0.5) ./ (L + 1)
        s  = mean(u) - 0.5
        z  = s / se_u
        drift = abs(s) * sqrt(12)
        b  = mean((r.pmean[:, p] .- r.θ[:, p]) ./ max.(r.psd[:, p], eps()))
        shr = mean(r.psd[:, p]) / max(std(r.θ[:, p]), eps())
        @printf("  %-18s %8.4f %8.2f %9.4f %+10.4f %10.3f %12.2e\n",
                LAB[p], mean(u), z, drift, b, shr, ks_p(u))
    end

    # ρ_true atom-aware KS (Spike 010 addendum): atoms are θ* exactly at the clamp endpoints.
    x = r.θ[:, 1]; lo, hi = minimum(x), maximum(x)
    atom = (x .== lo) .| (x .== hi)
    u1 = (r.ranks[:, 1] .+ 0.5) ./ (L + 1)
    Random.seed!(0xA70115)               # fixed, non-gate seed for the rank randomization
    ur = copy(u1); ur[atom] = rand(count(atom))
    @printf("\n  ρ_true atom handling: %d/%d atom draws (%.1f%%)\n",
            count(atom), M, 100 * count(atom) / M)
    @printf("    KS all draws       = %.3e\n", ks_p(u1))
    @printf("    KS atoms excluded  = %.3e   (n=%d)  <- valid SBC test\n",
            ks_p(u1[.!atom]), count(.!atom))
    @printf("    KS randomized rank = %.3e   (n=%d)\n", ks_p(ur), M)
end

# ---- realised mixture sanity (from the last evaluated net; identical DEV stream for both) -------
if !isempty(results)
    any_net = first(values(results))
    println("\n-- realised image-size mixture (DEV stream) --")
    for isz in SBC_IMSIZE_SET
        c = count(==(isz), any_net.szs)
        @printf("  %-12s %4d  (%.3f, target %.3f)\n", string(isz), c, c / M,
                SBC_IMSIZE_WEIGHTS[findfirst(==(isz), SBC_IMSIZE_SET)])
    end
end

# ---- persist every intermediate (store-the-intermediates convention, Spike 006) ----------------
save_kw = Dict{String,Any}("M" => M, "L" => L, "dev_seed" => DEV, "labels" => collect(LAB))
for (nm, r) in results
    key = replace(nm, r"[^a-zA-Z0-9]" => "_")
    save_kw["ranks_"*key] = r.ranks; save_kw["pmean_"*key] = r.pmean
    save_kw["psd_"*key]   = r.psd;   save_kw["theta_"*key] = r.θ
    save_kw["sizes_"*key] = r.szs
end
JLD2.save(joinpath(@__DIR__, "eval_data.jld2"), save_kw)
println("\nsaved eval_data.jld2")
