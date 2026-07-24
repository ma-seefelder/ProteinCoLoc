#########################################################################################
# Spike 013 --- DEV-seed evaluation: TRUNCATED-μ-prior net vs amended_v2.
#
# HONESTY UP FRONT. Truncating the μ-prior makes the truncated net a DIFFERENT statistical model
# (a different π(θ) ⇒ a different ρ_true support). So the amended_v2 comparison on ρ_true is NOT
# apples-to-apples: each net is evaluated under ITS OWN prior — the truncated net under the
# truncated prior, amended_v2 under the baseline prior. That is the CORRECT SBC for each. Only the
# COST measurement (§C below) is paired on identical images and is directly comparable.
#
# WHAT THIS MEASURES:
#   §A  ATOM FRACTION       — θ* draws at EXACTLY ±0.99: ~0% (truncated) vs ~6.5% (amended_v2).
#   §B  SBC (per own prior) — for each of the 8 columns (7 θ + paired-draw Δρ): mean-rank z,
#                             marginal drift (SD units), std bias, shrinkage, RAW KS p (NO
#                             randomized ranks — the whole point for ρ_true). ρ_true also reported
#                             atom-excluded / randomized for amended_v2 (Spike 010 correction).
#   §C  THE COST            — a COMMON high-|ρ| test set: both nets run on the SAME simulated
#                             images, compare |E[ρ̂]−ρ*| / RMSE / coverage / posterior-sd, binned by
#                             ρ_true. Does removing the ±0.99 training mass degrade strong-corr
#                             inference? This is the scientific trade-off the spike exists to quantify.
#
# DIAGNOSTIC ONLY. Fresh DEV seed asserted disjoint from every PROD/VAL/NPE/DEFAULT seed and every
# dev seed burned so far. NO gate report is written; NO PROD seed is consumed.
#
# Run:  julia --project -t auto .planning/spikes/013-mu-prior-truncation/eval_trunc_vs_amended.jl
#########################################################################################
using Random, Random123, Statistics, Printf, JLD2
using HypothesisTests, Distributions, StatsBase
using ProteinCoLoc
include(joinpath(@__DIR__, "trunc_datagen.jl"))     # sample_prior_trunc + MU_PRIOR_TRUNC

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const GATE = joinpath(REPO, "test", "gate")
include(joinpath(GATE, "gate_consts_8_v2.jl"))      # SBC mixture + PROD_SEED_V2 / PROD_SEED
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(GATE, "harness.jl"))

# ---- fresh DEV seed, asserted disjoint from EVERYTHING pre-registered / already burned ----------
const DEV = UInt64(0x00A1CE02)          # fresh; distinct from the 0x00DEC0DE training seed too
const NAMED_FORBIDDEN = UInt64[
    0x5BC0FFEE, 0xC0FFEE, 0x1,                                  # VAL / NPE / DEFAULT master seeds
    0xDE7C0DE, 0xDE7C0DE2, 0x0DE7C0D3,                          # burned dev seeds
    0x00B7A771, 0x00B7A772, 0x00CA11B0, 0x00CA11B1,             # burned dev seeds
    0x00A1CE01, 0x00CA11B2, 0x00CA11B3,                         # burned dev seeds
    0x00DEC0DE,                                                 # this spike's TRAINING seed
]
const FORBIDDEN = Set{UInt64}(vcat(NAMED_FORBIDDEN,
                                   collect(values(PROD_SEED_V2)),
                                   collect(values(PROD_SEED))))
@assert !(DEV in FORBIDDEN) "DEV seed $(repr(DEV)) collides with a forbidden seed"
@assert !(DEV in values(PROD_SEED_V2)) "DEV collides with a PROD_SEED_V2 value"
println("DEV seed = ", repr(DEV), "  (asserted disjoint from ", length(FORBIDDEN), " forbidden seeds)")

const M   = parse(Int, get(ENV, "SPIKE013_M", "1000"))
const L   = SBC_L                        # 999
const LAB = ("ρ_true", "spillover", "autofluorescence", "label_efficiency",
             "shift_dx", "shift_dy", "noise", "Δρ")

# Each net is evaluated under ITS OWN prior (correct SBC). Same forward model + build_mci; only the
# prior draw differs. The imsize mixture (F5) is the SAME for both, read from the frozen consts.
const SIM_BASE  = default_simulator()    # (; sample_prior, simulate_pair, build_mci)  -- baseline π
const SIM_TRUNC = (; sample_prior  = sample_prior_trunc,
                     simulate_pair = ProteinCoLoc.simulate_pair,
                     build_mci     = ProteinCoLoc.build_mci)   # truncated π

const NETS = (
    (name = "trunc (spike013)", prior = "truncated",
     path = joinpath(REPO, "artifacts", "spike013_trunc", "grid_8", "npe_8.jld2"), sim = SIM_TRUNC),
    (name = "amended_v2",       prior = "baseline",
     path = joinpath(REPO, "artifacts", "amended_v2",     "grid_8", "npe_8.jld2"), sim = SIM_BASE),
)

# ================================================================================================
# §B  SBC: M×L draw→simulate→infer under the net's OWN prior
# ================================================================================================
function run_eval(m, sim)
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
        println("MISSING: ", net.path, "  (skipping ", net.name, ")"); continue
    end
    println("Evaluating ", net.name, "  [prior=", net.prior, "]  M=", M, " L=", L)
    m = ProteinCoLoc.load_estimator(net.path)
    @printf("  flow: coupling=%d flow_depth=%d flow_width=%d  (summary dstar=%d depth=%d width=%d)\n",
            m.arch.num_coupling_layers, m.arch.flow_depth, m.arch.flow_width,
            m.arch.dstar, m.arch.depth, m.arch.width)
    r = run_eval(m, net.sim)
    results[net.name] = r

    # §A atom fraction of θ* draws (the point of the truncation)
    x = r.θ[:, 1]
    atom = (x .== -0.99) .| (x .== 0.99)
    @printf("\n  §A ATOM FRACTION (θ* ρ_true at exactly ±0.99): %d/%d = %.2f%%   range=[%.4f, %.4f]\n",
            count(atom), M, 100 * count(atom) / M, minimum(x), maximum(x))

    @printf("\n  §B  %-18s %8s %8s %9s %10s %10s %12s\n",
            "parameter", "mean(u)", "z", "drift_SD", "std_bias", "shrink", "KS_raw p")
    for p in 1:8
        u  = (r.ranks[:, p] .+ 0.5) ./ (L + 1)
        s  = mean(u) - 0.5
        z  = s / se_u
        drift = abs(s) * sqrt(12)
        bmn = mean((r.pmean[:, p] .- r.θ[:, p]) ./ max.(r.psd[:, p], eps()))
        shr = mean(r.psd[:, p]) / max(std(r.θ[:, p]), eps())
        @printf("  %-18s %8.4f %8.2f %9.4f %+10.4f %10.3f %12.2e\n",
                LAB[p], mean(u), z, drift, bmn, shr, ks_p(u))
    end

    # ρ_true atom-aware KS context (for amended_v2 the raw KS is the invalid-with-atoms number).
    u1 = (r.ranks[:, 1] .+ 0.5) ./ (L + 1)
    Random.seed!(0xA70213)               # fixed, non-gate seed for the rank randomization
    ur = copy(u1); ur[atom] = rand(count(atom))
    @printf("\n  ρ_true KS raw (all draws, NO randomization) = %.3e   <- the headline number\n", ks_p(u1))
    @printf("  ρ_true KS atoms-excluded                    = %.3e   (n=%d)\n",
            ks_p(u1[.!atom]), count(.!atom))
    @printf("  ρ_true KS randomized-rank                   = %.3e   (n=%d)\n", ks_p(ur), M)
end

# ================================================================================================
# §C  THE COST — a COMMON high-|ρ| test set, both nets on the SAME simulated images (paired).
# ================================================================================================
println("\n==================================================================")
println("§C  COST measurement — high-|ρ| inference on a COMMON simulated test set (paired)")

const RHO_TARGETS = (-0.99, -0.95, -0.90, -0.85, -0.80, 0.80, 0.85, 0.90, 0.95, 0.99)
const REPS = parse(Int, get(ENV, "SPIKE013_COSTREPS", "60"))

# Build a test θ with a PRESCRIBED ρ_true and nuisances drawn from the shared src/ nuisance priors
# (identical to both sample_prior and sample_prior_trunc — only ρ_true differs from the priors here).
function make_test_theta(rng, ρ0)
    return (ρ_true = ρ0,
            spillover        = rand(rng, ProteinCoLoc.SPILLOVER_PRIOR),
            autofluorescence = rand(rng, ProteinCoLoc.AUTOFLUORESCENCE_PRIOR),
            label_efficiency = rand(rng, ProteinCoLoc.LABEL_EFFICIENCY_PRIOR),
            shift_dx         = rand(rng, ProteinCoLoc.SHIFT_PRIOR),
            shift_dy         = rand(rng, ProteinCoLoc.SHIFT_PRIOR),
            noise            = rand(rng, ProteinCoLoc.NOISE_PRIOR))
end

# Load both nets (each carries its own frozen zt / θzt).
const COST_NETS = filter(n -> isfile(n.path), NETS)
cost_models = Dict(n.name => ProteinCoLoc.load_estimator(n.path) for n in COST_NETS)

# One shared DEV stream for the test θ + forward sims (net-independent ⇒ identical images per net).
rng_cost = Random123.Philox4x(UInt64, (DEV, UInt64(13)))
Random.seed!(DEV ⊻ UInt64(0x13))         # global stream for posterior base draws

# accumulate per (net, ρ0): |E[ρ̂]-ρ0| list, sq error, coverage90 count, sd list
cost = Dict{String,Any}()
for n in COST_NETS
    cost[n.name] = Dict(ρ0 => (abserr = Float64[], sqerr = Float64[], cov = Int[], sd = Float64[])
                        for ρ0 in RHO_TARGETS)
end

for ρ0 in RHO_TARGETS
    for _ in 1:REPS
        θt  = make_test_theta(rng_cost, ρ0)
        isz = ProteinCoLoc.sample_imsize(rng_cost; imsize_set = SBC_IMSIZE_SET,
                                         imsize_weights = SBC_IMSIZE_WEIGHTS)
        data = ProteinCoLoc.simulate_pair(rng_cost, θt; imsize = isz)   # SHARED image for both nets
        mci  = ProteinCoLoc.build_mci(data)
        summ = ProteinCoLoc.encode_d01(ProteinCoLoc.patch_summary(mci, 8))
        for n in COST_NETS
            m  = cost_models[n.name]
            Z  = ProteinCoLoc.standardize_summary(summ, m.zt, :min)      # each net's frozen zt
            dr = StatsBase.reconstruct(m.θzt,
                     ProteinCoLoc.posterior_for(m.estimator, Z; N = L, use_gpu = false))
            ρh = view(dr, 1, :)
            μh = mean(ρh)
            push!(cost[n.name][ρ0].abserr, abs(μh - ρ0))
            push!(cost[n.name][ρ0].sqerr, (μh - ρ0)^2)
            lo, hi = quantile(ρh, 0.05), quantile(ρh, 0.95)
            push!(cost[n.name][ρ0].cov, (lo <= ρ0 <= hi) ? 1 : 0)
            push!(cost[n.name][ρ0].sd, std(ρh))
        end
    end
    @printf("  ρ0=%+.2f done (%d reps)\n", ρ0, REPS)
end

println("\n  COST TABLE — mean |E[ρ̂]−ρ*| (accuracy), RMSE, 90%-coverage, mean posterior-sd")
for n in COST_NETS
    println("\n  net = ", n.name, "  [trained under ", n.prior, " prior]")
    @printf("    %8s %12s %12s %12s %12s\n", "ρ*", "mean|bias|", "RMSE", "cov90", "post_sd")
    for ρ0 in RHO_TARGETS
        c = cost[n.name][ρ0]
        @printf("    %+8.2f %12.4f %12.4f %12.3f %12.4f\n",
                ρ0, mean(c.abserr), sqrt(mean(c.sqerr)), mean(c.cov), mean(c.sd))
    end
end

# Head-to-head Δ(mean|bias|) = trunc − amended per ρ0 (positive ⇒ truncation is WORSE there).
if length(COST_NETS) == 2
    a, b = COST_NETS[1].name, COST_NETS[2].name
    tname = occursin("trunc", a) ? a : b
    aname = tname == a ? b : a
    println("\n  Δ(mean|bias|) = ", tname, " − ", aname, "   (positive ⇒ truncation WORSE at that ρ*)")
    for ρ0 in RHO_TARGETS
        dt = mean(cost[tname][ρ0].abserr) - mean(cost[aname][ρ0].abserr)
        @printf("    ρ0=%+.2f   Δ|bias| = %+.4f\n", ρ0, dt)
    end
end

# ---- realised mixture sanity (SBC stream) ------------------------------------------------------
if !isempty(results)
    any_net = first(values(results))
    println("\n-- realised image-size mixture (SBC DEV stream) --")
    for isz in SBC_IMSIZE_SET
        c = count(==(isz), any_net.szs)
        @printf("  %-12s %4d  (%.3f, target %.3f)\n", string(isz), c, c / M,
                SBC_IMSIZE_WEIGHTS[findfirst(==(isz), SBC_IMSIZE_SET)])
    end
end

# ---- persist every intermediate (store-the-intermediates convention, Spike 006) ----------------
save_kw = Dict{String,Any}("M" => M, "L" => L, "dev_seed" => DEV, "labels" => collect(LAB),
                           "rho_targets" => collect(RHO_TARGETS), "cost_reps" => REPS)
for (nm, r) in results
    key = replace(nm, r"[^a-zA-Z0-9]" => "_")
    save_kw["ranks_"*key] = r.ranks; save_kw["pmean_"*key] = r.pmean
    save_kw["psd_"*key]   = r.psd;   save_kw["theta_"*key] = r.θ
    save_kw["sizes_"*key] = r.szs
end
for (nm, d) in cost
    key = replace(nm, r"[^a-zA-Z0-9]" => "_")
    for ρ0 in RHO_TARGETS
        rk = replace(@sprintf("%+.2f", ρ0), "." => "p", "+" => "P", "-" => "N")
        save_kw["cost_$(key)_$(rk)_abserr"] = d[ρ0].abserr
        save_kw["cost_$(key)_$(rk)_sqerr"]  = d[ρ0].sqerr
        save_kw["cost_$(key)_$(rk)_cov"]    = d[ρ0].cov
        save_kw["cost_$(key)_$(rk)_sd"]     = d[ρ0].sd
    end
end
JLD2.save(joinpath(@__DIR__, "eval_data.jld2"), save_kw)
println("\nsaved eval_data.jld2")
