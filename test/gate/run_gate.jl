#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/run_gate.jl --- the per-grid ship-gate CLI (D-05).
#
#   run_gate.jl --grid G [--sbc] [--bf] [--ood] [--artifacts-root PATH]
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

if !isdefined(@__MODULE__, :SBC_M)
    let g = _peek_grid_arg(ARGS)
        pergrid = g === nothing ? nothing : joinpath(@__DIR__, "gate_consts_$(g).jl")
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

"""
    bf_gate(m, ratio; G, sim, rng = prod_rng(G), n = BF_SWEEP_N, L = SBC_L,
            imsize = SBC_IMSIZE, prior_n = 4000) -> NamedTuple

The amortized-Bayes-factor ship-gate: over `n` fresh-seeded paired draws it reads the amortized
log-BF (`amortized_log_bf(ratio.estimator, pair_encode(Zs, Zc), ratio.log_prior_odds)`, one CPU
forward pass) and the NON-CLAMPED KDE baseline log-BF (`kde_log_bf_unclamped`, over the Δρ
posterior vs a prior Δρ sample), then scores their correlation against `BF_CORR_MIN` and the
`max|Δ logBF|` against `BF_LOGBF_TOL` (non-finite pairs dropped). CPU-only, self-contained (no
Turing). Returns `(grid, n, corr, max_abs_delta, corr_pass, tol_pass, passed)`.
"""
function bf_gate(m, ratio; G::Integer, sim, rng = prod_rng(G), n::Integer = BF_SWEEP_N,
                 L::Integer = SBC_L, imsize = SBC_IMSIZE, prior_n::Integer = 4000)
    # Prior Δρ sample (independent θ*.ρ_true differences) for the KDE baseline denominator.
    prior_draws = Float64[sim.sample_prior(rng).ρ_true - sim.sample_prior(rng).ρ_true
                          for _ in 1:prior_n]

    amort = Float64[]
    base  = Float64[]
    for _ in 1:n
        pr       = draw_simulate_infer_paired(m, rng; G = G, imsize = imsize, N = L, sim = sim)
        a        = ProteinCoLoc.amortized_log_bf(ratio.estimator,
                       ProteinCoLoc.pair_encode(pr.Zs, pr.Zc), ratio.log_prior_odds)
        b        = ProteinCoLoc.kde_log_bf_unclamped([pr.ρs .- pr.ρc], prior_draws)[1]
        if isfinite(a) && isfinite(b)
            push!(amort, a)
            push!(base, b)
        end
    end

    corr    = (length(amort) >= 2) ? cor(amort, base) : NaN
    max_abs = isempty(amort) ? NaN : maximum(abs.(amort .- base))
    corr_ok = isfinite(corr) && corr >= BF_CORR_MIN
    tol_ok  = isfinite(max_abs) && max_abs <= BF_LOGBF_TOL
    return (grid = Int(G), n = length(amort), corr = corr, max_abs_delta = max_abs,
            corr_pass = corr_ok, tol_pass = tol_ok, passed = corr_ok && tol_ok)
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
                  L::Integer = SBC_L, imsize = SBC_IMSIZE,
                  pos_sim = nothing, neg_sim = nothing,
                  families = OOD_FAMILIES, levels::Integer = OOD_GRID_LEVELS,
                  n_pos::Integer = 40, negctrl_reps::Integer = 30)
    density = haskey(ood_nulls, :density) ? ood_nulls.density : ood_nulls

    # --- single-simulator override (one pooled AUC, density channel) -------------------------
    if pos_sim !== nothing
        id = Float64[]
        for _ in 1:n
            t = draw_simulate_infer(m, rng; G = G, imsize = imsize, N = L, sim = sim)
            push!(id, ProteinCoLoc.maha_score(density, t.Z))
        end
        thr = ProteinCoLoc.id_threshold(id; q = OOD_ID_QUANTILE)
        pos = Float64[]
        for _ in 1:n
            t = draw_simulate_infer(m, rng; G = G, imsize = imsize, N = L, sim = pos_sim)
            push!(pos, ProteinCoLoc.maha_score(density, t.Z))
        end
        auc    = ProteinCoLoc.roc_auc(id, pos)[3]   # (fpr, tpr, auc) — AUC is element 3
        auc_ok = isfinite(auc) && auc >= OOD_AUC_MIN
        return (grid = Int(G), n = n, id_threshold = thr, auc = auc,
                auc_pass = auc_ok, passed = auc_ok, mode = :pos_sim)
    end

    # --- the reported family × level ROC grid ------------------------------------------------
    roc = gate_ood_roc(m, density; G = G, rng = rng, families = families, levels = levels,
                       n_id = n, n_pos = n_pos, imsize = imsize)

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
            base = ProteinCoLoc.simulate_pair(rng, ProteinCoLoc.sample_prior(rng);
                                              imsize = imsize)
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
"""
function run_gate(grid::Integer; sbc::Bool = true, bf::Bool = false, ood::Bool = false,
                  artifacts_root = ProteinCoLoc.default_artifacts_root(),
                  sim = nothing, write_report::Bool = true,
                  M::Integer = SBC_M, L::Integer = SBC_L, bins::Integer = SBC_BINS,
                  imsize = SBC_IMSIZE, ood_n::Integer = 200,
                  pos_sim = nothing, neg_sim = nothing)
    npe_path = _gate_npe_path(artifacts_root, grid)
    if !isfile(npe_path)
        @info "run_gate: grid not trained yet — no NPE artifact; the full gate runs once a " *
              "per-grid plan trains this grid" grid npe_path
        return (status = :not_trained, grid = Int(grid), npe_path = npe_path)
    end

    simr = sim === nothing ? default_simulator() : sim
    m    = load_gate_model(npe_path)

    sbc_report = sbc ? sbc_gate(m; G = grid, M = M, L = L, bins = bins, imsize = imsize,
                                sim = simr, rng = prod_rng(grid)) : nothing

    bf_report = nothing
    if bf
        rpath = _gate_ratio_path(artifacts_root, grid)
        isfile(rpath) || error("run_gate: --bf needs the ratio artifact $rpath")
        ratio = ProteinCoLoc.load_ratio(rpath)
        bf_report = bf_gate(m, ratio; G = grid, sim = simr, rng = prod_rng(grid),
                            L = L, imsize = imsize)
    end

    ood_report = nothing
    if ood
        opath = _gate_ood_path(artifacts_root, grid)
        isfile(opath) || error("run_gate: --ood needs the OOD-nulls artifact $opath")
        nulls = ProteinCoLoc.load_ood_nulls(opath)
        ood_report = ood_gate(m, nulls; G = grid, sim = simr, rng = prod_rng(grid),
                              n = ood_n, L = L, imsize = imsize,
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

    report = (status = :ran, grid = Int(grid), seed = prod_seed(grid),
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
