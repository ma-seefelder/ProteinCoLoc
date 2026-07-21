#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/run_gate.jl --- the per-grid ship-gate CLI (D-05).
#
#   run_gate.jl --grid G [--sbc] [--bf] [--ood] [--artifacts-root PATH] [--consts FILE]
#
# Loads the grid's FRESH pre-registration (`gate_consts_<G>.jl`, else the template), loads the
# grid's CPU-resident frozen net via `ProteinCoLoc.load_estimator`, runs the SELECTED reported-
# scale gates through the CPU harness (every inference `use_gpu = false`), and writes a per-grid
# gate report through the atomic `.tmp` → reopen-integrity → `mv(...; force=true)` wrapper.
#
# SCOPE (07-04): this plan delivers the gate MACHINERY the per-grid plans invoke — no grid is
# trained yet. So `run_gate.jl --grid 8 --sbc` is INVOCABLE: if the grid's net artifact is absent
# it returns an honest `:not_trained` status (no crash); the full reported run happens once a
# per-grid plan trains the grid. The BF gate scores the amortized log-BF against the NON-CLAMPED
# KDE baseline (self-contained, no Turing); the OOD gate computes the ID operating point and, when
# the per-grid plan injects positive/negative control simulators, the control-separability AUC.

using ProteinCoLoc
import JLD2
import Statistics
import Statistics: cor

# Ordered, guarded includes: consts → harness → sbc. A per-grid `gate_consts_<G>.jl` included
# BEFORE this file is respected. When invoked as `run_gate.jl --grid G`, PREFER the committed
# per-grid pre-registration `gate_consts_<G>.jl` (D-05: the reported gate must be driven by the
# grid's OWN locked consts + fresh PROD_SEED[G], not the shared template) and fall back to the
# committed template only if no per-grid file exists. A standalone/runtests include (empty ARGS)
# falls back to the template.
"""
    _peek_grid_arg(args) -> Union{Int,Nothing}

Scan CLI `args` for `--grid G` before the consts are loaded, so the correct per-grid
pre-registration file can be selected at include time (the grid is otherwise parsed later).
"""
function _peek_grid_arg(args)
    for i in 1:(length(args) - 1)
        args[i] == "--grid" && return tryparse(Int, args[i + 1])
    end
    return nothing
end

"""
    _peek_consts_arg(args) -> Union{String,Nothing}

Scan CLI `args` for an EXPLICIT `--consts FILE` pre-registration override (resolved inside
`test/gate/`). The amended grid-8 pre-registration (`gate_consts_8_v2.jl`) is deliberately NOT
selected automatically by `--grid 8`: the amendment binds it to exactly ONE run on the fresh
`PROD_SEED_V2[8]` (protocol §6.3), so consuming it must be an explicit, auditable invocation
(`--grid 8 --consts gate_consts_8_v2.jl`) and never a side effect of the default path.
"""
function _peek_consts_arg(args)
    for i in 1:(length(args) - 1)
        args[i] == "--consts" && return args[i + 1]
    end
    return nothing
end

if !isdefined(@__MODULE__, :SBC_M)
    let g = _peek_grid_arg(ARGS), c = _peek_consts_arg(ARGS)
        pergrid = c !== nothing ? joinpath(@__DIR__, basename(c)) :
                  g === nothing ? nothing : joinpath(@__DIR__, "gate_consts_$(g).jl")
        if pergrid !== nothing && isfile(pergrid)
            @info "run_gate: loading per-grid pre-registration" file = basename(pergrid)
            include(pergrid)
        else
            include(joinpath(@__DIR__, "gate_consts_template.jl"))
        end
    end
end
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))
isdefined(@__MODULE__, :sbc_ranks)           || include(joinpath(@__DIR__, "sbc.jl"))
isdefined(@__MODULE__, :OOD_FAMILIES)        || include(joinpath(@__DIR__, "misspec.jl"))

# --- artifact paths (mirror src/amortized/pipeline.jl's grid store layout) -------------------
_gate_grid_dir(root, G)  = joinpath(root, "grid_$(G)")
_gate_npe_path(root, G)  = joinpath(_gate_grid_dir(root, G), "npe_$(G).jld2")
_gate_ratio_path(root, G)= joinpath(_gate_grid_dir(root, G), "ratio_$(G).jld2")
_gate_ood_path(root, G)  = joinpath(_gate_grid_dir(root, G), "ood_nulls_$(G).jld2")
_gate_report_path(root, G) = joinpath(_gate_grid_dir(root, G), "gate_report_$(G).jld2")

# On-disk schema tag for the gate report (bump on any breaking layout change).
const GATE_REPORT_SCHEMA = 1

"""
    write_gate_report(path, report) -> String

Atomically persist a per-grid gate `report` (a NamedTuple) through the SAME `.tmp` → reopen-
integrity-`@assert` → `mv(...; force=true)` wrapper the estimator/cache artifacts use: write to
`path.tmp`, reopen to integrity-check the `report` key, then filesystem-atomically rename. A crash
before the `mv` leaves only a discardable `.tmp`. Returns `path`.
"""
function write_gate_report(path, report)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    JLD2.jldsave(tmp; schema_version = GATE_REPORT_SCHEMA, report = report)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "report") "write_gate_report: integrity check failed, $tmp missing report"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    read_gate_report(path) -> Union{NamedTuple,Nothing}

Reload a previously written per-grid gate report, or `nothing` if absent/unreadable. Used by the
ARM-PRESERVING MERGE in `run_gate`: an arm that was not re-run keeps its recorded verdict.
"""
function read_gate_report(path)
    isfile(path) || return nothing
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "report") ? f["report"] : nothing
        end
    catch
        return nothing
    end
end

# Previously recorded verdict for arm `k` (`nothing` when there is no readable prior report).
_prev_arm(prev, k::Symbol) = (prev === nothing || !hasproperty(prev, k)) ? nothing :
                             getproperty(prev, k)

# =============================================================================================
# A2 — the BF correlation rule (07-GATE-AMENDMENT §3)
# =============================================================================================
#
# DEFECT 1 (precision). v1 ran `n = BF_SWEEP_N = 25` paired draws (recorded n 15–20 after dropping
# non-finite pairs) and compared a POINT estimate to 0.95. Var(atanh r̂) ≈ 1/(n−3) makes the 95 %
# Fisher-z CI ±0.07–0.10 wide there: the gate cannot discriminate ρ = 0.94 from ρ = 0.95. Derivable
# from n and the threshold alone — no result required.
#
# DEFECT 2 (extreme-value instability). `max|Δ logBF|` is a MAXIMUM, and E[max] is non-decreasing
# in n for any non-degenerate distribution. A fixed tolerance on a maximum therefore silently
# TIGHTENS whenever n is raised — which fixing defect 1 must do — for reasons unrelated to
# agreement quality.
#
# THE RULES. n = BF_GATE_N = 100 (derived from a pre-registered half-width ≤ 0.03 target at
# ρ = 0.95 ⇒ n_min = 69); the verdict is a ONE-SIDED 95 % LOWER CONFIDENCE BOUND on ρ, which is
# STRICTLY HARDER than v1 (it demands r̂ ≥ 0.9639 at n = 100); the gated tail statistic becomes the
# n-stable `quantile(|Δ logBF|, 0.95)` while `max|Δ logBF|` stays REPORTED. BF_CORR_MIN = 0.95 and
# BF_LOGBF_TOL = 0.5 are NOT re-tuned.
#
# Both rule sets are computed on the SAME paired draws and reported side by side, for the same
# attribution reason as the SBC arm.

"Fisher-z CI fallback (used when a v1 consts file is loaded); the frozen `fisher_z_ci` wins."
function gate_fisher_z_ci(r::Real, n::Integer; tail::Symbol = :two)
    (isfinite(r) && n > 3) || return (lo = -1.0, hi = 1.0)
    zc = tail === :two ? 1.959963985 : 1.644853627
    z  = atanh(clamp(float(r), -0.999999, 0.999999))
    d  = zc / sqrt(n - 3)
    return (lo = tanh(z - d), hi = tanh(z + d))
end

"Three-way LB verdict fallback; the frozen `bf_corr_verdict` of the v2 consts wins when loaded."
function gate_bf_corr_verdict(r::Real, n::Integer; thr::Real = BF_CORR_MIN, nmin::Integer = 69)
    n >= nmin || return :invalid
    isfinite(r) || return :invalid
    ci = gate_fisher_z_ci(r, n; tail = :one)
    ci.lo >= thr && return :pass
    ci.hi <  thr && return :fail
    return :inconclusive
end

# Prefer the FROZEN implementations of the loaded pre-registration; fall back only under v1 consts.
_bf_ci(r, n; tail = :two) = isdefined(@__MODULE__, :fisher_z_ci) ?
    getfield(@__MODULE__, :fisher_z_ci)(r, n; tail = tail) : gate_fisher_z_ci(r, n; tail = tail)
_bf_verdict(r, n) = isdefined(@__MODULE__, :bf_corr_verdict) ?
    getfield(@__MODULE__, :bf_corr_verdict)(r, n) : gate_bf_corr_verdict(r, n)

# The amended constants, with v1 fallbacks so this file loads under either pre-registration.
const GATE_BF_N       = isdefined(@__MODULE__, :BF_GATE_N)     ? BF_GATE_N     : BF_SWEEP_N
const GATE_BF_N_MIN   = isdefined(@__MODULE__, :BF_GATE_N_MIN) ? BF_GATE_N_MIN : 0
const GATE_BF_LOGBF_Q = isdefined(@__MODULE__, :BF_LOGBF_Q)    ? BF_LOGBF_Q    : 0.95
const BF_RULES_VERSION = isdefined(@__MODULE__, :BF_DECISION_RULE) &&
                         BF_DECISION_RULE === :lower_confidence_bound ? 2 : 1

"""
    bf_gate(m, ratio; G, sim, rng = prod_rng(G), n = GATE_BF_N, L = SBC_L,
            imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
            imsize_weights = GATE_IMSIZE_WEIGHTS, prior_n = 4000) -> NamedTuple

The amortized-Bayes-factor ship-gate: over `n` fresh-seeded paired draws (ONE shared image size per
pair under the F5 mixture) it reads the amortized log-BF (`amortized_log_bf(ratio.estimator,
pair_encode(Zs, Zc), ratio.log_prior_odds)`, one CPU forward pass) and the NON-CLAMPED KDE baseline
log-BF (`kde_log_bf_unclamped`, over the Δρ posterior vs a prior Δρ sample), non-finite pairs
dropped. CPU-only, self-contained (no Turing).

Reports BOTH rule sets on the SAME draws — `v1_verdict` (point estimate `r̂ ≥ BF_CORR_MIN` and
`max|Δ| ≤ BF_LOGBF_TOL`) and `v2_verdict` (one-sided 95 % lower confidence bound on ρ, and
`q95(|Δ|) ≤ BF_LOGBF_TOL`) — with `rules_version` naming the binding one. `max_abs_delta` is
retained as a REPORTED quantity under both, alongside `q50`/`q95`, so the KDE tail artifact stays
inspectable rather than defined away.
"""
function bf_gate(m, ratio; G::Integer, sim, rng = prod_rng(G), n::Integer = GATE_BF_N,
                 L::Integer = SBC_L, imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
                 imsize_weights = GATE_IMSIZE_WEIGHTS, prior_n::Integer = 4000)
    # REPRODUCIBILITY (07-10): the ρ draws behind every log-BF come from the GLOBAL stream
    # (`sampleposterior` threads no rng — see harness.jl). Pin it from the frozen PROD_SEED[G].
    seed_gate_global!(G)
    # Prior Δρ sample (independent θ*.ρ_true differences) for the KDE baseline denominator.
    prior_draws = Float64[sim.sample_prior(rng).ρ_true - sim.sample_prior(rng).ρ_true
                          for _ in 1:prior_n]

    amort = Float64[]
    base  = Float64[]
    sizes = Tuple{Int,Int}[]
    for _ in 1:n
        pr       = draw_simulate_infer_paired(m, rng; G = G, imsize = imsize,
                                              imsize_set = imsize_set,
                                              imsize_weights = imsize_weights, N = L, sim = sim)
        push!(sizes, pr.imsize)
        a        = ProteinCoLoc.amortized_log_bf(ratio.estimator,
                       ProteinCoLoc.pair_encode(pr.Zs, pr.Zc), ratio.log_prior_odds)
        b        = ProteinCoLoc.kde_log_bf_unclamped([pr.ρs .- pr.ρc], prior_draws)[1]
        if isfinite(a) && isfinite(b)
            push!(amort, a)
            push!(base, b)
        end
    end

    nfin    = length(amort)
    corr    = (nfin >= 2) ? cor(amort, base) : NaN
    absΔ    = isempty(amort) ? Float64[] : abs.(amort .- base)
    max_abs = isempty(absΔ) ? NaN : maximum(absΔ)
    q50     = isempty(absΔ) ? NaN : Statistics.quantile(absΔ, 0.50)
    q95     = isempty(absΔ) ? NaN : Statistics.quantile(absΔ, GATE_BF_LOGBF_Q)

    # --- v1 (original) rules, recomputed on THESE draws purely for attribution ----------------
    v1_corr_ok = isfinite(corr) && corr >= BF_CORR_MIN
    v1_tol_ok  = isfinite(max_abs) && max_abs <= BF_LOGBF_TOL
    v1 = (rules = 1, corr_pass = v1_corr_ok, tol_pass = v1_tol_ok,
          statistic = :max_abs_delta, passed = v1_corr_ok && v1_tol_ok)

    # --- v2 (amended) rules: one-sided 95 % LOWER confidence bound + n-stable q95 tail ---------
    ci_two     = _bf_ci(corr, nfin; tail = :two)
    ci_one     = _bf_ci(corr, nfin; tail = :one)
    verdict    = _bf_verdict(corr, nfin)
    v2_corr_ok = verdict === :pass
    v2_tol_ok  = isfinite(q95) && q95 <= BF_LOGBF_TOL
    v2 = (rules = 2, decision_rule = :lower_confidence_bound, corr_verdict = verdict,
          corr_lb = ci_one.lo, corr_ub = ci_one.hi, ci_two_sided = ci_two,
          n_min = GATE_BF_N_MIN, statistic = :q95_abs_delta, q = GATE_BF_LOGBF_Q,
          corr_pass = v2_corr_ok, tol_pass = v2_tol_ok,
          passed = v2_corr_ok && v2_tol_ok)

    binding = BF_RULES_VERSION == 2 ? v2 : v1
    return (grid = Int(G), n = nfin, n_attempted = Int(n), corr = corr,
            max_abs_delta = max_abs, q50_abs_delta = q50, q95_abs_delta = q95,
            imsize_set = imsize_set, imsize_weights = imsize_weights,
            realised_imsize_counts = realised_imsize_counts(sizes),
            rules_version = BF_RULES_VERSION, v1_verdict = v1, v2_verdict = v2,
            corr_pass = binding.corr_pass, tol_pass = binding.tol_pass,
            passed = binding.passed)
end

"""
    ood_gate(m, ood_nulls; G, sim, rng = prod_rng(G), n = 200, L = SBC_L, imsize = SBC_IMSIZE,
             pos_sim = nothing, neg_sim = nothing, families = OOD_FAMILIES,
             levels = OOD_GRID_LEVELS, n_pos = 40, negctrl_reps = 30) -> NamedTuple

The OOD / misspecification ship-gate — the SCORED arm (07-09b).

The positive controls are the ESTABLISHED Phase-5 misspecification families (`OOD_FAMILIES` in
`misspec.jl`: texture / noise / optics / background, four magnitude levels each), swept through the
grid-generalized controlled ROC experiment `gate_ood_roc` and scored with the OR-FUSED
density∨noise detector (the Phase-5 iter2 design; the posterior-predictive channel is excluded from
the reported fusion exactly as the reported spike run did). The gated statistic is the pooled
strongest-level separability AUC against the PRE-REGISTERED `OOD_AUC_MIN`, together with each
family's best-level AUC — the same conjunction the spike's reported OOD gate used.

The summary-orthogonal NEGATIVE controls (affine / rotate-flip / block-permute, D-04) are also
measured: their KS-invariance statistic against `OOD_KS_EPS` and their flag fire-rate against the
ID false-positive rate `1 − OOD_ID_QUANTILE`. These transforms are PROVABLY invisible to a fixed
patch-Pearson correlation summary — a near-0.5 separability on them is the NAMED structural blind
spot, reported, never hidden.

`pos_sim` (a single misspecified simulator obeying the `(sample_prior, simulate_pair, build_mci)`
contract) remains supported as an override for a one-family probe; when given it replaces the family
grid with a single pooled AUC. CPU-only. Returns a NamedTuple whose `passed` is the conjunction of
the pooled-AUC gate and the per-family AUC gates; the negative-control measurements are recorded
alongside but are NOT folded into `passed` (see the comment at the verdict).
"""
function ood_gate(m, ood_nulls; G::Integer, sim, rng = prod_rng(G), n::Integer = 200,
                  L::Integer = SBC_L, imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
                  imsize_weights = GATE_IMSIZE_WEIGHTS,
                  pos_sim = nothing, neg_sim = nothing,
                  families = OOD_FAMILIES, levels::Integer = OOD_GRID_LEVELS,
                  n_pos::Integer = 40, negctrl_reps::Integer = 30)
    # REPRODUCIBILITY (07-10): the OOD arm's scores ride draw_simulate_infer / gate_ood_roc, whose
    # posterior draws come from the GLOBAL stream (`sampleposterior` threads no rng — see
    # harness.jl). Pin it from the frozen PROD_SEED[G] so a re-run reproduces its own AUC table.
    seed_gate_global!(G)
    density =haskey(ood_nulls, :density) ? ood_nulls.density : ood_nulls
    sizes = Tuple{Int,Int}[]      # realised image sizes, for the report's realised_imsize_counts

    # --- single-simulator override (one pooled AUC, density channel) -------------------------
    if pos_sim !== nothing
        id = Float64[]
        for _ in 1:n
            t = draw_simulate_infer(m, rng; G = G, imsize = imsize, imsize_set = imsize_set,
                                    imsize_weights = imsize_weights, N = L, sim = sim)
            push!(sizes, t.imsize)
            push!(id, ProteinCoLoc.maha_score(density, t.Z))
        end
        thr = ProteinCoLoc.id_threshold(id; q = OOD_ID_QUANTILE)
        pos = Float64[]
        for _ in 1:n
            t = draw_simulate_infer(m, rng; G = G, imsize = imsize, imsize_set = imsize_set,
                                    imsize_weights = imsize_weights, N = L, sim = pos_sim)
            push!(sizes, t.imsize)
            push!(pos, ProteinCoLoc.maha_score(density, t.Z))
        end
        auc    = ProteinCoLoc.roc_auc(id, pos)[3]   # (fpr, tpr, auc) — AUC is element 3
        auc_ok = isfinite(auc) && auc >= OOD_AUC_MIN
        return (grid = Int(G), n = n, id_threshold = thr, auc = auc,
                auc_pass = auc_ok, passed = auc_ok, mode = :pos_sim,
                imsize_set = imsize_set, imsize_weights = imsize_weights,
                realised_imsize_counts = realised_imsize_counts(sizes))
    end

    # --- the reported family × level ROC grid ------------------------------------------------
    roc = gate_ood_roc(m, density; G = G, rng = rng, families = families, levels = levels,
                       n_id = n, n_pos = n_pos, imsize = imsize, imsize_set = imsize_set,
                       imsize_weights = imsize_weights)
    append!(sizes, roc.imsizes)

    fam_best = Dict(fam => maximum(roc.fused_auc[fam]) for fam in roc.families)
    fam_pass = Dict(fam => (isfinite(fam_best[fam]) && fam_best[fam] >= OOD_AUC_MIN)
                    for fam in roc.families)
    auc      = roc.combined_auc
    auc_ok   = isfinite(auc) && auc >= OOD_AUC_MIN

    # --- summary-orthogonal negative controls (D-04 named blind spot) ------------------------
    negs   = gate_negctrls(G; rng_factory = () -> prod_rng(G))
    neg_ks_max    = Dict{Symbol,Float64}()
    neg_fire_rate = Dict{Symbol,Float64}()
    neg_pass      = Dict{Symbol,Bool}()
    for (name, tf) in pairs(negs)
        ksvals = Float64[]; fires = 0
        for _ in 1:negctrl_reps
            θn   = ProteinCoLoc.sample_prior(rng)
            isz  = gate_imsize(rng; imsize = imsize, imsize_set = imsize_set,
                               imsize_weights = imsize_weights)
            push!(sizes, isz)
            base = ProteinCoLoc.simulate_pair(rng, θn; imsize = isz)
            push!(ksvals, verify_summary_invariance(base, tf, G))
            fires += (ProteinCoLoc.maha_score(density, gate_summary(m, tf(base), G)) >
                      roc.id_threshold) ? 1 : 0
        end
        ksmax = maximum(ksvals)
        fr    = fires / negctrl_reps
        neg_ks_max[name]    = ksmax
        neg_fire_rate[name] = fr
        neg_pass[name]      = (ksmax < OOD_KS_EPS) && (fr <= (1 - OOD_ID_QUANTILE) + 1e-9)
    end

    # The GATED verdict is the pre-registered SEPARABILITY criterion: the pooled strongest-level
    # AUC and every family's best-level AUC against OOD_AUC_MIN. The negative controls are
    # MEASURED and RECORDED (`neg_ks_max`/`neg_fire_rate`/`neg_pass`) but deliberately kept OUT of
    # `passed`: they quantify the D-04 structural blind spot of a correlation-only summary, and
    # their KS statistic is granularity-limited on coarse grids (a G×G summary has only G² values,
    # so the smallest non-zero KS statistic is 1/G² = 0.0625 at G=4 — already above OOD_KS_EPS
    # before any real discrepancy exists). Reported, not hidden, not folded into the AUC verdict.
    passed = auc_ok && all(values(fam_pass))
    return (grid = Int(G), n = n, n_pos = n_pos, levels = levels, mode = :family_grid,
            id_threshold = roc.id_threshold, fused_threshold = roc.fused_threshold,
            id_fire_rate = roc.id_fire_rate,
            auc = auc, auc_pass = auc_ok,
            families = roc.families, fam_best_auc = fam_best, fam_pass = fam_pass,
            fused_auc = roc.fused_auc, density_auc = roc.density_auc,
            noise_auc = roc.noise_auc, fire_rate = roc.fire_rate,
            neg_ks_max = neg_ks_max, neg_fire_rate = neg_fire_rate, neg_pass = neg_pass,
            youden_j = roc.youden.j,
            imsize_set = imsize_set, imsize_weights = imsize_weights,
            realised_imsize_counts = realised_imsize_counts(sizes),
            passed = passed)
end

"""
    run_gate(grid; sbc = true, bf = false, ood = false, artifacts_root = default root,
             sim = nothing, write_report = true, M = SBC_M, L = SBC_L, ...) -> NamedTuple

The per-grid ship-gate dispatcher. Resolves the grid's CPU-resident artifacts under
`artifacts_root/grid_G/`; if the NPE artifact is ABSENT returns `(status = :not_trained, grid)`
(invocable before any grid is trained — the full run happens in the per-grid plans). Otherwise
loads the frozen net (and ratio / OOD nulls for the selected gates) and runs the SELECTED gates
(`use_gpu = false` throughout), optionally writing the report atomically. `sim` defaults to the
promoted `default_simulator()` (inject a fake simulator until the forward model is promoted).

Under a pre-registration that sets `SBC_REQUIRE_IMSIZE_PROVENANCE` (the v2 amendment) the gate
first asserts that the net's PERSISTED training image-size distribution equals the gate mixture and
returns `(status = :provenance_mismatch, …)` — running NO arm — when it is absent or different.

The report carries `imsize_set`, `imsize_weights` and `realised_imsize_counts` (pooled over the
arms that ran) alongside the net's `training_imsize_provenance`.
"""
function run_gate(grid::Integer; sbc::Bool = true, bf::Bool = false, ood::Bool = false,
                  artifacts_root = ProteinCoLoc.default_artifacts_root(),
                  sim = nothing, write_report::Bool = true,
                  M::Integer = SBC_M, L::Integer = SBC_L, bins::Integer = SBC_BINS,
                  chi2_bins::Integer = _gate_chi2_bins(bins),
                  imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
                  imsize_weights = GATE_IMSIZE_WEIGHTS, ood_n::Integer = 200,
                  require_provenance::Bool = GATE_REQUIRE_IMSIZE_PROVENANCE,
                  pos_sim = nothing, neg_sim = nothing)
    npe_path = _gate_npe_path(artifacts_root, grid)
    if !isfile(npe_path)
        @info "run_gate: grid not trained yet — no NPE artifact; the full gate runs once a " *
              "per-grid plan trains this grid" grid npe_path
        return (status = :not_trained, grid = Int(grid), npe_path = npe_path)
    end

    simr = sim === nothing ? default_simulator() : sim
    m    = load_gate_model(npe_path)

    # --- BINDING INVARIANT (07-GATE-AMENDMENT §4): gate joint == training joint ---------------
    # SBC rank uniformity is a theorem about the joint the estimator was TRAINED under. If the
    # gate's image-size distribution differs from the net's persisted training distribution — or
    # the net has no recorded provenance at all, which is the case for every v1-frozen bundle —
    # the SBC arm would prove calibration on a distribution the net will never see. The gate
    # therefore REFUSES, before any arm executes. This is the mechanical enforcement of F5;
    # without it the validity condition is a comment.
    prov = assert_imsize_provenance(m; imsize_set = imsize_set,
                                    imsize_weights = imsize_weights,
                                    require = require_provenance)
    if !prov.ok
        @warn "run_gate: REFUSING to run — training image-size provenance absent or different " *
              "from the gate mixture (07-GATE-AMENDMENT §4)" grid training = prov.training gate = prov.gate
        return (status = :provenance_mismatch, grid = Int(grid), npe_path = npe_path,
                training_imsize_provenance = prov.training,
                imsize_set = imsize_set, imsize_weights = imsize_weights)
    end

    sbc_report = sbc ? sbc_gate(m; G = grid, M = M, L = L, bins = bins, chi2_bins = chi2_bins,
                                imsize = imsize, imsize_set = imsize_set,
                                imsize_weights = imsize_weights,
                                sim = simr, rng = prod_rng(grid)) : nothing

    bf_report = nothing
    if bf
        rpath = _gate_ratio_path(artifacts_root, grid)
        isfile(rpath) || error("run_gate: --bf needs the ratio artifact $rpath")
        ratio = ProteinCoLoc.load_ratio(rpath)
        bf_report = bf_gate(m, ratio; G = grid, sim = simr, rng = prod_rng(grid),
                            L = L, imsize = imsize, imsize_set = imsize_set,
                            imsize_weights = imsize_weights)
    end

    ood_report = nothing
    if ood
        opath = _gate_ood_path(artifacts_root, grid)
        isfile(opath) || error("run_gate: --ood needs the OOD-nulls artifact $opath")
        nulls = ProteinCoLoc.load_ood_nulls(opath)
        ood_report = ood_gate(m, nulls; G = grid, sim = simr, rng = prod_rng(grid),
                              n = ood_n, L = L, imsize = imsize, imsize_set = imsize_set,
                              imsize_weights = imsize_weights,
                              pos_sim = pos_sim, neg_sim = neg_sim)
    end

    # ARM-PRESERVING MERGE (07-09b): a partial run (e.g. `--ood` alone) must NOT erase the arms it
    # did not re-run. Any UNSELECTED arm's previously recorded verdict is carried forward verbatim
    # from the existing report; a selected arm always overwrites. Without this, re-running one arm
    # would silently null out the other two recorded verdicts.
    rpath = _gate_report_path(artifacts_root, grid)
    prev  = read_gate_report(rpath)
    sbc || (sbc_report = _prev_arm(prev, :sbc))
    bf  || (bf_report  = _prev_arm(prev, :bf))
    ood || (ood_report = _prev_arm(prev, :ood))

    # The realised mixture, pooled over the arms that actually ran (F6: record, never reconstruct).
    realised = Dict{Tuple{Int,Int},Int}()
    for arm in (sbc_report, bf_report, ood_report)
        (arm isa NamedTuple && hasproperty(arm, :realised_imsize_counts)) || continue
        for (k, v) in arm.realised_imsize_counts
            realised[k] = get(realised, k, 0) + v
        end
    end

    report = (status = :ran, grid = Int(grid), seed = prod_seed(grid),
              gate_consts_version = isdefined(@__MODULE__, :GATE_CONSTS_VERSION) ?
                                    GATE_CONSTS_VERSION : 1,
              imsize_set = imsize_set, imsize_weights = imsize_weights,
              realised_imsize_counts = realised,
              training_imsize_provenance = prov.training,
              sbc = sbc_report, bf = bf_report, ood = ood_report)
    if write_report
        write_gate_report(rpath, report)
    end
    return report
end

# --- CLI --------------------------------------------------------------------------------------

"""
    _parse_gate_args(args) -> NamedTuple

Parse `--grid G [--sbc] [--bf] [--ood] [--artifacts-root PATH]`. `--grid` is required; if no gate
flag is given, `--sbc` is the default.
"""
function _parse_gate_args(args)
    grid = nothing
    sbc = false; bf = false; ood = false
    root = ProteinCoLoc.default_artifacts_root()
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--grid"
            i < length(args) || error("run_gate: --grid needs a value")
            grid = parse(Int, args[i + 1]); i += 2
        elseif a == "--sbc"; sbc = true; i += 1
        elseif a == "--bf";  bf  = true; i += 1
        elseif a == "--ood"; ood = true; i += 1
        elseif a == "--artifacts-root"
            i < length(args) || error("run_gate: --artifacts-root needs a value")
            root = args[i + 1]; i += 2
        elseif a == "--consts"
            # Consumed at include time by `_peek_consts_arg`; accepted (and ignored) here.
            i < length(args) || error("run_gate: --consts needs a value")
            i += 2
        else
            error("run_gate: unknown argument $a")
        end
    end
    grid === nothing && error("run_gate: --grid G is required")
    (sbc || bf || ood) || (sbc = true)   # default to the SBC gate
    return (grid = grid, sbc = sbc, bf = bf, ood = ood, artifacts_root = root)
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = _parse_gate_args(ARGS)
    rep  = run_gate(opts.grid; sbc = opts.sbc, bf = opts.bf, ood = opts.ood,
                    artifacts_root = opts.artifacts_root)
    @info "run_gate complete" grid = opts.grid status = rep.status
end
