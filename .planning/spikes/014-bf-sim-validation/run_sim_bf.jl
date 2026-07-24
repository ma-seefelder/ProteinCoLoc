#=
Spike 014 --- Simulation-based Bayes-factor validation (proof of concept).

Proves that the amortized NRE log-BF can be validated the way SBC validates a posterior --
by DISCRIMINATION + DECISION-CALIBRATION over labelled simulated pairs -- WITHOUT the
per-pair KDE-tail baseline (which attrits/saturates; spikes 006-008).

READ-ONLY: src/ and spike/ are never edited. Uses the FROZEN amended_v2 grid-8 NRE
(artifacts/amended_v2/grid_8/ratio_8.jld2) + the frozen zt/θzt from npe_8.jld2. NO retraining.

COLOC/NULL DEFINITION (extracted, src/amortized/train_ratio.jl:107,131 + :55):
  the RatioEstimator was trained on pairs pair_encode(Zstd[:,i], Zstd[:,j]) of two INDEPENDENT
  prior draws, labelled m = Float32((ρ[i] − ρ[j]) > RATIO_SPLIT_THRESHOLD) with threshold 0.0.
  ⇒ H1 (coloc) ⟺ ρ_true_i > ρ_true_j ;  H0 (null) ⟺ ρ_true_i ≤ ρ_true_j.
  Evidence strength = signed Δρ = ρ_i − ρ_j. Our test pairs are constructed IDENTICALLY
  (two i.i.d. prior draws, sign-of-contrast label), so the validation is on-definition.

IMAGE SIZES: the F5 mixture read from the FROZEN gate_consts_8_v2.jl (SBC_IMSIZE_SET/WEIGHTS),
which EQUALS the training-joint imsize provenance recorded in the amended_v2 artifacts' meta
(verified: ((512,512),(1024,1024),(1376,1028),(2048,2048)) w=(0.4,0.25,0.25,0.1)). SBC rank
uniformity / NRE calibration hold only under the training joint, so we simulate at that mixture.

SEED: DEV seed 0x00BF5014 ("BF-Sim-014"), @assert-ed disjoint from every forbidden seed.
=#

import ProteinCoLoc as PC
import Random123: Philox4x
import Random
import Statistics: mean, median, quantile
import StatsBase: tiedrank, corspearman
import KernelDensity: kde
import QuadGK: quadgk
import Distributions: pdf
import JLD2

# ------------------------------------------------------------------------------------------
# SEED DISCIPLINE (CONVENTIONS.md: DEV seed, @assert disjoint from the full forbidden set)
# ------------------------------------------------------------------------------------------
const DEV_SEED = 0x0000_0000_00BF_5014   # fresh DEV seed for spike 014

# Load the FROZEN amended pre-registration into an ISOLATED module purely to READ its
# PROD_SEED / PROD_SEED_V2 / SBC_IMSIZE_* (no gate is run; the file is not modified).
module _GC
    include(joinpath(@__DIR__, "..", "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end

const _FORBIDDEN = Set{UInt64}(vcat(
    UInt64(0x0000_0000_0000_0001),                 # DEFAULT_MASTER_SEED
    UInt64(0x0000_0000_00C0_FFEE),                 # NPE_MASTER_SEED
    UInt64(0x0000_0000_5BC0_FFEE),                 # VAL_MASTER_SEED
    UInt64(0x0000_0000_0000_0001),
    # dev seeds burned so far (enumerated in the task brief)
    UInt64(0x0000_0000_0DE7_C0DE), UInt64(0x0000_0000_DE7C_0DE2),
    UInt64(0x0000_0000_0DE7_C0D3), UInt64(0x0000_0000_00DE_C0DE),
    UInt64(0x0000_00C0_FFEE_5BC0), UInt64(0x0000_0000_C0FF_EE5B),
    UInt64(0x0000_0000_00B7_A770), UInt64(0x0000_0000_00B7_A771),
    UInt64(0x0000_0000_00B7_A772), UInt64(0x0000_0000_00B7_A773),
    UInt64(0x0000_0000_00CA_11B0), UInt64(0x0000_0000_00CA_11B1),
    UInt64(0x0000_0000_00CA_11B2), UInt64(0x0000_0000_00CA_11B3),
    UInt64(0x0000_0000_00A1_CE01), UInt64(0x0000_0000_00A1_CE02),
    # ratio-pair / datagen salts
    UInt64(0x0000_0000_004A_7107),
    # every v1 + v2 pre-registered gate seed (read from the frozen file)
    collect(UInt64, values(_GC.PROD_SEED))...,
    collect(UInt64, values(_GC.PROD_SEED_V2))...,
))
@assert !(UInt64(DEV_SEED) in _FORBIDDEN) "DEV_SEED collides with a forbidden seed"
@assert !(UInt64(DEV_SEED) in Set(values(_GC.PROD_SEED_V2))) "DEV_SEED collides with PROD_SEED_V2"

# F5 image-size mixture, READ from the frozen amendment (single source of truth).
const IMSET = _GC.SBC_IMSIZE_SET
const IMWTS = _GC.SBC_IMSIZE_WEIGHTS
@assert IMSET == ((512, 512), (1024, 1024), (1376, 1028), (2048, 2048))

# ------------------------------------------------------------------------------------------
# Load the FROZEN nets (NO retraining).
# ------------------------------------------------------------------------------------------
const ART = joinpath(@__DIR__, "..", "..", "..", "artifacts", "amended_v2", "grid_8")
println("Loading frozen NRE + NPE (zt) from ", ART)
const RATIO = PC.load_ratio(joinpath(ART, "ratio_8.jld2"))
const NPE   = PC.load_estimator(joinpath(ART, "npe_8.jld2"))
const ZT    = NPE.zt
const ΘZT   = NPE.θzt
# Verify training-joint consistency (F5): the net's recorded imsize provenance == our sim mixture.
let prov = PC.training_imsize_provenance(NPE)
    @assert prov.imsize_set == IMSET "training imsize_set != F5 mixture -- covariate shift!"
    @assert prov.imsize_weights == IMWTS "training imsize_weights != F5 weights"
    println("training-joint imsize provenance OK: ", prov.imsize_set, " w=", prov.imsize_weights)
end
println("log_prior_odds (deployed) = ", RATIO.log_prior_odds)

const GRID = 8

# ------------------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------------------
const M       = parse(Int, get(ENV, "SPIKE014_M", "400"))   # number of test pairs
const N_POST  = parse(Int, get(ENV, "SPIKE014_NPOST", "1000")) # NPE Δρ posterior draws (KDE side)
const K_PRIOR = 20_000                                        # Δρ prior draws (KDE prior-odds)
const CLAMP_F = 1e-8                                          # spike baseline probability floor
const SAT_BOUND = log((1 - CLAMP_F) / CLAMP_F)               # ≈ 18.4207 saturation rail
println("M=$M  N_POST=$N_POST  threads=$(Threads.nthreads())")

# ------------------------------------------------------------------------------------------
# Phase A --- simulate 2M INDEPENDENT members (Philox-per-index ⇒ threadable & reproducible).
# Each member mirrors generate_sample: prior → imsize → simulate → summary → frozen-zt standardize.
# ------------------------------------------------------------------------------------------
function sim_member(idx::Int)
    rng = Philox4x(UInt64, (UInt64(DEV_SEED), UInt64(idx)))
    θ   = PC.sample_prior(rng)
    isz = PC.sample_imsize(rng; imsize_set = IMSET, imsize_weights = IMWTS)
    mci = PC.build_mci(PC.simulate_pair(rng, θ; imsize = isz))
    Zraw = PC.encode_d01(PC.patch_summary(mci, GRID))
    Z    = PC.standardize_summary(Zraw, ZT, :min)     # 2·G²=128-dim Float32 (deployed input space)
    return (Z = Z, ρ = θ.ρ_true, imsize = isz)
end

nmem = 2M
Zmem  = Vector{Vector{Float32}}(undef, nmem)
ρmem  = Vector{Float64}(undef, nmem)
szmem = Vector{Tuple{Int,Int}}(undef, nmem)
println("Phase A: simulating $nmem members ...")
tA = @elapsed begin
    Threads.@threads for idx in 1:nmem
        m = sim_member(idx)
        Zmem[idx]  = m.Z
        ρmem[idx]  = m.ρ
        szmem[idx] = m.imsize
    end
end
println("  Phase A done in $(round(tA; digits=1)) s")

# Pair consecutive members: pair k = (sample 2k-1, control 2k). Both i.i.d. ⇒ training-consistent.
sidx = collect(1:2:nmem)
cidx = collect(2:2:nmem)
ρs   = ρmem[sidx]
ρc   = ρmem[cidx]
drho = ρs .- ρc                       # signed Δρ = ρ_sample − ρ_control (evidence)
label = drho .> 0.0                   # H1 (coloc) ⟺ Δρ > 0  (train_ratio.jl:131, thr 0.0)
println("H1 (coloc) fraction = ", round(mean(label); digits=3), "  (balanced by i.i.d. construction)")

# ------------------------------------------------------------------------------------------
# Phase B --- per pair: (1) amortized NRE log-BF (deterministic forward pass), and
#             (2) NPE Δρ posterior draws for the KDE-baseline contrast (global RNG, serial).
# ------------------------------------------------------------------------------------------
Random.seed!(UInt64(DEV_SEED))        # pin the GLOBAL stream sampleposterior draws from (KDE side)

logbf_nre = Vector{Float64}(undef, M)
post_draws = Vector{Vector{Float64}}(undef, M)
println("Phase B: NRE log-BF + NPE posterior draws for $M pairs ...")
tB = @elapsed begin
    for k in 1:M
        Zs = Zmem[sidx[k]]; Zc = Zmem[cidx[k]]
        Zp = PC.pair_encode(Zs, Zc)                                   # A7 diff encoding (320-dim)
        logbf_nre[k] = PC.amortized_log_bf(RATIO.estimator, Zp, RATIO.log_prior_odds)
        ρs_d = PC.rho_draws(NPE.estimator, Zs, ΘZT; N = N_POST)
        ρc_d = PC.rho_draws(NPE.estimator, Zc, ΘZT; N = N_POST)
        post_draws[k] = ρs_d .- ρc_d                                  # D-03 MC Δρ cloud
    end
end
println("  Phase B done in $(round(tB; digits=1)) s")

# Δρ PRIOR draws (fixed; baseline prior-odds sample), on a disjoint sub-stream of DEV_SEED.
prng = Philox4x(UInt64, (UInt64(DEV_SEED) ⊻ 0xA5A5_A5A5_A5A5_A5A5, UInt64(0)))
prior_draws = Vector{Float64}(undef, K_PRIOR)
for j in 1:K_PRIOR
    prior_draws[j] = PC.sample_prior(prng).ρ_true - PC.sample_prior(prng).ρ_true
end

# ------------------------------------------------------------------------------------------
# KDE baseline on the SAME pairs: unclamped (src) + a clamped variant (1e-8 floor).
# ------------------------------------------------------------------------------------------
# Unclamped: src/amortized/bf.jl kde_log_bf_unclamped (may return ±Inf = the dropped value).
logbf_kde_unclamped = PC.kde_log_bf_unclamped(post_draws, prior_draws; threshold = 0.0)

# Clamped variant (mirrors the demoted spike baseline's _clampp=1e-8 floor → ±SAT_BOUND rail).
_p_gt_clamped(draws; thr = 0.0) = begin
    dist = kde(collect(float.(draws)))
    p_le, _ = quadgk(x -> pdf(dist, x), -Inf, thr)
    clamp(1 - p_le, CLAMP_F, 1 - CLAMP_F)
end
pp_prior_c = _p_gt_clamped(prior_draws)
prior_odds_c = pp_prior_c / (1 - pp_prior_c)
logbf_kde_clamped = map(post_draws) do d
    pp = _p_gt_clamped(d)
    log((pp / (1 - pp)) / prior_odds_c)
end

# ------------------------------------------------------------------------------------------
# METRICS
# ------------------------------------------------------------------------------------------
# (1) Discrimination: AUC(H1 vs H0) via the Mann–Whitney rank identity.
function auc_score(scores::AbstractVector{<:Real}, labels::AbstractVector{Bool})
    pos = scores[labels]; neg = scores[.!labels]
    n1 = length(pos); n0 = length(neg)
    (n1 == 0 || n0 == 0) && return NaN
    r = tiedrank(vcat(pos, neg))
    (sum(@view r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
end
auc_nre = auc_score(logbf_nre, label)
# AUC of the KDE baseline on the pairs it can score (finite only) -- for context.
finite_mask = isfinite.(logbf_kde_unclamped)
auc_kde_finite = auc_score(logbf_kde_unclamped[finite_mask], label[finite_mask])

# ROC curve for the NRE score (for the figure).
function roc_points(scores, labels)
    thr = sort(unique(scores); rev = true)
    pushfirst!(thr, Inf)
    n1 = count(labels); n0 = count(.!labels)
    tpr = Float64[]; fpr = Float64[]
    for t in thr
        pred = scores .>= t
        push!(tpr, count(pred .& labels) / n1)
        push!(fpr, count(pred .& .!labels) / n0)
    end
    (fpr = fpr, tpr = tpr)
end
roc = roc_points(logbf_nre, label)

# (2) Monotonicity of log-BF in evidence Δρ.
spear = corspearman(logbf_nre, drho)
# Binned median log-BF vs signed Δρ (equal-count bins).
nbin = 9
edges = quantile(drho, range(0, 1; length = nbin + 1))
bin_centers = Float64[]; bin_med = Float64[]; bin_lo = Float64[]; bin_hi = Float64[]
for b in 1:nbin
    lo = edges[b]; hi = edges[b+1]
    m = b < nbin ? (drho .>= lo) .& (drho .< hi) : (drho .>= lo) .& (drho .<= hi)
    any(m) || continue
    push!(bin_centers, median(drho[m]))
    push!(bin_med, median(logbf_nre[m]))
    push!(bin_lo, quantile(logbf_nre[m], 0.25))
    push!(bin_hi, quantile(logbf_nre[m], 0.75))
end
mono_ok = issorted(bin_med)   # strictly nondecreasing binned medians?

# (3) Decision calibration at the logBF>0 threshold (declare coloc).
pred0 = logbf_nre .> 0.0
acc0  = mean(pred0 .== label)
fpr0  = count(pred0 .& .!label) / max(count(.!label), 1)   # declare H1 | true H0
fnr0  = count(.!pred0 .& label) / max(count(label), 1)     # declare H0 | true H1
# Reliability of p̂(coloc)=sigmoid(posterior log-odds); lpo≈0 so ≈ sigmoid(logBF).
sig(x) = 1 / (1 + exp(-x))
phat = sig.(logbf_nre .+ RATIO.log_prior_odds)   # posterior_log_odds = logBF + 2·lpo ≈ logBF
relb = 10
rel_edges = range(0, 1; length = relb + 1)
rel_p = Float64[]; rel_emp = Float64[]; rel_n = Int[]
ece = 0.0
for b in 1:relb
    lo = rel_edges[b]; hi = rel_edges[b+1]
    m = b < relb ? (phat .>= lo) .& (phat .< hi) : (phat .>= lo) .& (phat .<= hi)
    nb = count(m); nb == 0 && continue
    mp = mean(phat[m]); me = mean(label[m])
    push!(rel_p, mp); push!(rel_emp, me); push!(rel_n, nb)
    global ece += nb / M * abs(mp - me)
end

# (4) KDE-baseline pathology on the SAME pairs (what the sim-based route sidesteps).
frac_nonfinite = mean(.!finite_mask)                       # dropped by the gate's isfinite filter
frac_saturated = mean(abs.(logbf_kde_clamped) .>= SAT_BOUND - 1e-6)
# Spike-006 signature: are the DROPPED pairs the high-|Δρ| (high-signal) ones?
absd = abs.(drho)
med_absd_dropped = any(.!finite_mask) ? median(absd[.!finite_mask]) : NaN
med_absd_kept    = any(finite_mask)  ? median(absd[finite_mask])  : NaN
# NRE keeps 100% of pairs (bounded network output) -- confirm.
frac_nre_nonfinite = mean(.!isfinite.(logbf_nre))

# ------------------------------------------------------------------------------------------
# REPORT
# ------------------------------------------------------------------------------------------
println("\n==================== SPIKE 014 RESULTS ====================")
println("DEV_SEED               = 0x", string(DEV_SEED; base = 16, pad = 8))
println("M pairs                = ", M, "   (H1 frac = ", round(mean(label); digits = 3), ")")
println("--- Discrimination ---")
println("AUC(NRE logBF, H1vsH0) = ", round(auc_nre; digits = 4))
println("AUC(KDE finite only)   = ", round(auc_kde_finite; digits = 4),
        "   (on ", count(finite_mask), "/", M, " finite pairs)")
println("--- Monotonicity in evidence Δρ ---")
println("Spearman(logBF, Δρ)    = ", round(spear; digits = 4))
println("binned medians nondecr = ", mono_ok)
println("--- Decision calibration @ logBF>0 ---")
println("accuracy               = ", round(acc0; digits = 4))
println("FPR (H1|H0)            = ", round(fpr0; digits = 4))
println("FNR (H0|H1)            = ", round(fnr0; digits = 4))
println("reliability ECE (p̂)    = ", round(ece; digits = 4))
println("--- KDE-baseline pathology on the SAME pairs ---")
println("non-finite (attrited)  = ", round(frac_nonfinite; digits = 4),
        "   (", count(.!finite_mask), "/", M, ")")
println("saturated (clamped)    = ", round(frac_saturated; digits = 4))
println("median|Δρ| dropped     = ", round(med_absd_dropped; digits = 3))
println("median|Δρ| kept        = ", round(med_absd_kept; digits = 3))
println("NRE non-finite         = ", round(frac_nre_nonfinite; digits = 4), "  (keeps ALL pairs)")
println("===========================================================\n")

# ------------------------------------------------------------------------------------------
# Persist plain arrays (for figenv/CairoMakie + reproducibility). artifacts/ is gitignored.
# ------------------------------------------------------------------------------------------
outdir = joinpath(@__DIR__, "..", "..", "..", "artifacts", "spike014")
mkpath(outdir)
outpath = joinpath(outdir, "results.jld2")
JLD2.jldsave(outpath;
    dev_seed = UInt64(DEV_SEED), M = M, n_post = N_POST,
    drho = drho, label = collect(label),
    logbf_nre = logbf_nre,
    logbf_kde_unclamped = logbf_kde_unclamped,
    logbf_kde_clamped = logbf_kde_clamped,
    roc_fpr = roc.fpr, roc_tpr = roc.tpr,
    bin_centers = bin_centers, bin_med = bin_med, bin_lo = bin_lo, bin_hi = bin_hi,
    rel_p = rel_p, rel_emp = rel_emp, rel_n = rel_n,
    auc_nre = auc_nre, auc_kde_finite = auc_kde_finite, spear = spear,
    acc0 = acc0, fpr0 = fpr0, fnr0 = fnr0, ece = ece,
    frac_nonfinite = frac_nonfinite, frac_saturated = frac_saturated,
    med_absd_dropped = med_absd_dropped, med_absd_kept = med_absd_kept,
    sat_bound = SAT_BOUND,
)
println("saved: ", outpath)
