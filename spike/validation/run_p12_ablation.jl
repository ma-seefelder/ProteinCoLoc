# =============================================================================================
# D-13 ABLATION — THE DESCOPE DELIVERABLE, PRODUCED UNDER A STANDALONE USER AUTHORISATION
# =============================================================================================
#
# WHAT THIS IS. The D-13 fallback: a per-region Δρ map WITH per-region uncertainty, trained with
# the spatial prior NEUTRALIZED (`arm = :none`, so `lattice_sigma` returns the identity and the
# regions do not borrow). It is strictly more than `LocalColocMap` offers today, which carries no
# uncertainty at all.
#
# WHY THIS FILE EXISTS AT ALL, AND WHY IT IS NOT 12-14 TASK 3.
# 12-14 Task 3 is the only place in the plans that defines this recipe, and it is keyed on
# `p12_stage1_verdict() === :descope`. **The Stage-1 verdict is `:proceed` and this file does not
# change that.** The descope trigger fired at the *prior-selection* step (12-15,
# `NONE-BEATS-ABLATION` in three runs), which is a state the plans' executable routing does not
# model:
#
#   * `12-CONTEXT.md:310-312` defines D-12 **Stage 1** as the CAR-vs-GP mini-spike itself —
#     "If neither does, descope **before** the full training run" — which is verbatim what
#     happened; but
#   * the **implemented** Stage-1 gate is 12-11's ridge borrowing probe, which returned PROCEED.
#
# The specified gate and the implemented gate diverged, and the trigger fired on the specified one.
# The route is therefore a NEW RECORD, never an edited one: `12-STAGE1-VERDICT.md` is NOT touched,
# and this run is authorised by
# `.planning/phases/12-spatial-colocalization-map/12-D13-AUTHORISATION.md`.
#
# WRITING `VERDICT: DESCOPE` INTO THE STAGE-1 VERDICT TO MAKE TASK 3 FIRE WOULD BE FALSIFYING A
# GATE RECORD TO OBTAIN A ROUTE. It is not done here, and `stage1_verdict` in the artifact below
# is written from `p12_stage1_verdict()` at run time — i.e. it will read `:proceed`, the truth,
# not the `:descope` that 12-14 Task 3's schema anticipated. See the DEVIATIONS block.
#
# WHY 12-16 WILL STILL FIND THIS. 12-16 resolves its bundle by **file existence, not by the
# verdict** (`12-16-PLAN.md:77-81`): it scores `p12_train_full_report.jld2` if present, else the
# bundle recorded in `p12_ablation_report.jld2`, else throws. So writing this artifact under its
# contracted name and key (`bundle_path`) is sufficient, and NO gate value has to be invented.
#
# SCALE LIMIT, NAMED HERE RATHER THAN LEFT TO BE DISCOVERED: this is the **SPIKE-SCALE** D-13
# fallback at `P12_MINISPIKE_N` = 10,000 pairs, **not** the 50,000-pair version 12-17 would have
# built. `12-CONTEXT.md:323-324` already warns that the descope ablation "exists only at spike
# scale" and that the fallback is "fully available only at the Stage-2 gate".
#
# DEVIATIONS FROM 12-14 TASK 3, DECLARED RATHER THAN SILENT
# ---------------------------------------------------------------------------------------------
# (1) `stage1_verdict` is `:proceed`, not `:descope`. Task 3's schema hard-codes `:descope`
#     because it can only run on that branch. Recording `:descope` here would be false. The true
#     value is written, plus `authorisation` / `authorisation_doc` naming the actual route.
#
# (2) TASK 3'S `realized_r1_quantiles` ACCEPTANCE CRITERION CANNOT BE SATISFIED BY TASK 3'S OWN
#     INSTRUCTION, AND THIS IS A PLAN DEFECT, NOT AN EXECUTOR CHOICE. The criterion requires the
#     quantiles be "degenerate at `P12_ABLATION_R1`" (= 0.0), "proving the ablation really was
#     trained on the neutralized prior". But the instruction is
#     `generate_p12_pool(P12_MINISPIKE_N; arm = :none)`, which leaves `r1 = nothing` — i.e. r1 is
#     **DRAWN** from `Uniform(P12_R1_MIN, P12_R1_MAX)`. Neutralization in this generator does NOT
#     happen by pinning r1: `p12_lattice.jl:311` returns `Symmetric(Matrix(1.0I))` for
#     `arm === :none` **whatever r1 is**. The criterion tests a proxy this generator does not use.
#
#     RESOLVED BY FOLLOWING THE INSTRUCTION AND RECORDING THE TRUTH, because the instruction also
#     matches the precedent that actually ran: 12-15's `:none` arm drew r1 too, so a drawn-r1
#     ablation is the one whose numbers were measured and reported. `realized_r1_quantiles` is
#     therefore written honestly (it will be ~uniform, NOT degenerate) and
#     `neutralization_mechanism` records how neutralization is actually achieved, so the
#     criterion's PURPOSE is evidenced even though its literal test is not met.
#
# (3) EPOCHS = 100, not 12-15's 18. The 12-15 re-run established that at an 18-epoch budget every
#     arm ends ON ITS BUDGET rather than on convergence, while at a 100-epoch budget the `:none`
#     arm stops early at 36 — i.e. it converges. A SHIPPED DELIVERABLE SHOULD BE TRAINED TO
#     CONVERGENCE, not stopped by its budget. Recorded in the artifact as `epochs` (the budget)
#     alongside `risk_trace`, from which the epoch it actually stopped at is readable.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p12_ablation.jl

using JLD2
using Dates
using SHA
using Statistics

isdefined(@__MODULE__, :P12_DEV_SEED)  || include(joinpath(@__DIR__, "p12_consts.jl"))
isdefined(@__MODULE__, :train_p12_npe) || include(joinpath(@__DIR__, "..", "npe", "train_p12_npe.jl"))

const ABLATION_REPORT_PATH = joinpath(@__DIR__, "p12_ablation_report.jld2")

# See DEVIATION (3). The budget, not the number of epochs that will run.
const ABLATION_EPOCHS = 100

function main()
    println("="^78)
    println("D-13 ABLATION — the descope deliverable (arm = :none, spatial prior NEUTRALIZED)")
    println("SPIKE SCALE: n = $(P12_MINISPIKE_N), NOT the 50,000-pair 12-17 version")
    println("="^78)
    t0 = time()

    # THE GATE IS READ AND REPORTED, NEVER WRITTEN. This run does not require `:descope` and does
    # not pretend to have it.
    verdict = p12_stage1_verdict()
    println("Stage-1 verdict (READ, not written): :$verdict")
    verdict === :absent && error("run_p12_ablation: the Stage-1 gate is :absent. This run is " *
        "authorised relative to an ADJUDICATED gate, not in place of one.")
    println("Authorised by 12-D13-AUTHORISATION.md (user ruling, 2026-07-31) — NOT by a " *
            ":descope verdict, which would have to be forged to exist.")

    # ---- POOL: arm = :none on the F5 mixture (both are `generate_p12_pool` defaults) ----------
    td  = time()
    dir = generate_p12_pool(P12_MINISPIKE_N; arm = :none)     # imsize_tag = :f5_mixture by default
    @assert p12_pool_complete(dir, P12_MINISPIKE_N)
    datagen_min = (time() - td) / 60
    println("pool: $dir  ($(round(datagen_min; digits = 2)) min)")

    # ---- TRAIN on the WHOLE pool. 12-16 draws its own fresh datasets (12-16-PLAN:306), so no
    #      held-out block is carved here; that was 12-15's need, not this one. -----------------
    tt = time()
    bundle = train_p12_npe(; arm = :none, n = P12_MINISPIKE_N, pool_dir = dir,
                           epochs = ABLATION_EPOCHS, verbose = true)
    train_min = (time() - tt) / 60

    # ---- DURABLE SAVE under P12_PRIMARY_CHECKOUT (12-17's rule: survives a worktree cleanup).
    #      No `worktree` substring test — the git-derived root IS the check (p12_consts.jl:540-544).
    bundle_path = joinpath(P12_PRIMARY_CHECKOUT, "spike", "artifacts", "p12", "p12_ablation_none.jld2")
    @assert startswith(normpath(bundle_path), normpath(P12_PRIMARY_CHECKOUT))
    save_p12_npe(bundle_path, bundle)
    bundle_sha256 = bytes2hex(open(sha256, bundle_path))
    println("bundle: $bundle_path\n  sha256 = $bundle_sha256")

    # ---- The realized r1 distribution, recorded HONESTLY. See DEVIATION (2). -----------------
    pool = load_p12_pool(dir)
    r1vals = Float64.(pool.r1)                       # the pool carries r1 as its own field
    r1q = quantile(r1vals, [0.0, 0.25, 0.5, 0.75, 1.0])
    r1_degenerate = maximum(abs.(r1vals .- P12_ABLATION_R1)) < 1e-12
    println("realized_r1_quantiles = ", round.(r1q; digits = 5),
            "   degenerate at P12_ABLATION_R1 = $r1_degenerate  ",
            "(FALSE is EXPECTED here — see DEVIATION 2 in this file's header)")

    imcounts = Dict{String,Int}()
    for s in pool.imsize                             # Vector{Tuple{Int,Int}}
        imcounts[string(s)] = get(imcounts, string(s), 0) + 1
    end

    elapsed_min = (time() - t0) / 60

    tmp = ABLATION_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version = 1,
        arm = :none,
        r1_prior = P12_ABLATION_R1,
        bundle_path = bundle_path,          # THE KEY 12-16 READS. Name and meaning must not drift.
        bundle_sha256 = bundle_sha256,
        pool_dir = dir,
        n_pool = P12_MINISPIKE_N,
        epochs = ABLATION_EPOCHS,
        risk_trace = bundle.risk_trace,     # the epoch it ACTUALLY stopped at is readable here
        train_seed = bundle.train_seed,
        realized_mask_rate = bundle.realized_mask_rate,
        realized_imsize_counts = imcounts,
        realized_r1_quantiles = r1q,
        realized_r1_degenerate_at_ablation_r1 = r1_degenerate,
        neutralization_mechanism = "arm = :none => lattice_sigma returns Symmetric(Matrix(1.0I)) " *
            "at p12_lattice.jl:311, INDEPENDENT of r1. Neutralization is structural (the " *
            "covariance is the identity, so regions do not borrow), NOT achieved by pinning r1. " *
            "12-14 Task 3's acceptance criterion asks for realized_r1_quantiles DEGENERATE at " *
            "P12_ABLATION_R1, but its own instruction (generate_p12_pool(n; arm = :none)) leaves " *
            "r1 DRAWN. The criterion tests a proxy this generator does not use; the instruction " *
            "was followed and the truth recorded. 12-15's :none arm drew r1 the same way, so " *
            "this deliverable matches the ablation whose numbers were reported.",
        # THE GATE, AS READ -- never as needed.
        stage1_verdict = verdict,
        stage1_verdict_note = "TRUE value, read from p12_stage1_verdict() at run time. 12-14 " *
            "Task 3's schema hard-codes :descope because it can only run on that branch; this " *
            "run is NOT on that branch and recording :descope would be false. " *
            "12-STAGE1-VERDICT.md is byte-unchanged.",
        authorisation = :user_standalone_2026_07_31,
        authorisation_doc = ".planning/phases/12-spatial-colocalization-map/12-D13-AUTHORISATION.md",
        scale_limit = "SPIKE SCALE (n = $(P12_MINISPIKE_N)), not the 50,000-pair 12-17 version. " *
            "12-CONTEXT.md:323-324 warns the descope ablation exists only at spike scale.",
        D = bundle.D, G = P12_G, K_dev = P12_K_DEV,
        datagen_min = datagen_min,
        train_min = train_min,
        elapsed_min = elapsed_min,
        datagen_ceiling_min = P12_DATAGEN_WALLCLOCK_CEILING_MIN,
        generated = string(Dates.now(Dates.UTC)) * "Z",
        julia_version = string(VERSION),
        caption = "THE D-13 DESCOPE DELIVERABLE: the neutralized-prior (arm = :none) ablation, " *
            "trained at spike scale on the F5 image-size mixture, saved durably under " *
            "P12_PRIMARY_CHECKOUT. Produced under a standalone user authorisation after 12-15 " *
            "returned NONE-BEATS-ABLATION in three runs; the Stage-1 verdict remains :proceed " *
            "and was not edited. This is the fallback release, not the spatial map.",
    )
    let d = JLD2.load(tmp)
        for k in ("arm", "bundle_path", "bundle_sha256", "realized_r1_quantiles", "pool_dir")
            @assert haskey(d, k) "artifact integrity: $k missing"
        end
        @assert isfile(d["bundle_path"]) "artifact integrity: bundle_path does not resolve"
    end
    mv(tmp, ABLATION_REPORT_PATH; force = true)

    println("\n" * "="^78)
    println("D-13 ABLATION WRITTEN")
    println("="^78)
    println("  datagen $(round(datagen_min; digits = 2)) min against " *
            "P12_DATAGEN_WALLCLOCK_CEILING_MIN = $(P12_DATAGEN_WALLCLOCK_CEILING_MIN) (PER POOL)")
    println("  train   $(round(train_min; digits = 2)) min   total $(round(elapsed_min; digits = 2)) min")
    println("  wrote ", ABLATION_REPORT_PATH)
    return nothing
end

isdefined(@__MODULE__, :P12_ABLATION_LOAD_ONLY) || main()
