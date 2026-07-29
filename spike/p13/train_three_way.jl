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

# spike/p13/train_three_way.jl --- REPORTED Phase-13 training run (D-07, D-09, D-10, D-11).
#
# =============================================================================================
# ANTI-SNOOPING BANNER. READ BEFORE CHANGING ANY NUMBER IN THIS FILE.
# =============================================================================================
# THE RECIPE WAS COPIED VERBATIM, NOT CHOSEN. Every training hyperparameter is READ from the
# frozen pre-registration `spike/p13/consts.jl`, which transcribed them from the SHIPPED binary
# trainer's own defaults at `src/amortized/train_ratio.jl:155-163`. Not one of them was selected
# by looking at a Phase-13 result. THAT IS THE ENTIRE BASIS OF D-10's ATTRIBUTION ARGUMENT: if
# the three-way net's two-way overlap differs from the shipped binary net's, the difference is
# attributable TO THE HEAD, because nothing else moved. Re-tuning even one value -- the learning
# rate, the batch size, the patience -- silently converts a comparison into a coincidence.
# NO LITERAL RECIPE VALUE APPEARS IN THIS FILE'S EXECUTABLE SOURCE, and the test suite asserts
# that, so a value cannot be "temporarily" inlined and left behind.
#
# THIS IS **ONE** RUN. A second run with different hyperparameters would spend
# `P13_ITERATION_ALLOWANCE = 1`, which `P13_ITERATION_TRIGGER` RESERVES for switching the
# stratification design to `P13_STRATIFICATION_FALLBACK` (D-07-ii) on one pre-declared condition,
# and which is currently UNSPENT. IF THIS RUN FAILS TO CONVERGE, THAT IS A FINDING FOR THE
# REPORT, NOT A LICENCE TO RETRAIN. Phase 7 was amended twice after seeing results and the
# credibility cost of that is why the allowance is a constant rather than a judgement call.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY. It composes surfaces that already exist and are already
# tested -- `stratify_by_class`/`p13_generate_pool` (datagen.jl), `head_log_odds` (labels.jl),
# `build_three_way_net`/`train_three_way`/`save_three_way` (net.jl) -- and adds no reusable
# abstraction of its own. Including it does nothing; it must be run deliberately:
#
#     julia --project=spike -t auto spike/p13/train_three_way.jl
#
# THE THREE THINGS THIS RUNNER IS RESPONSIBLE FOR GETTING RIGHT, EACH OF WHICH IS SILENT IF WRONG:
#   1. THE PER-HEAD CORRECTION IS MEASURED ON THE **TRAINING** SUBSET (D-07). Measuring it on the
#      evaluation set would tune the reported evidence scale to the data it is scored against;
#      measuring it on the whole pool would count columns the net never trained on. It is also
#      never ASSUMED to be zero -- see `measure_head_log_odds`' docstring.
#   2. `Random.seed!` IS CALLED IMMEDIATELY BEFORE `train_three_way` (Pitfall 6). Flux's shuffling
#      `DataLoader`, which NeuralEstimators' `_DataLoader` builds with `shuffle = true`
#      (`train.jl:668-676`), draws from the GLOBAL RNG. Threading a counter-based stream through
#      generation is layer 1 and is NOT sufficient on its own; this is layer 2
#      (`P13_GLOBAL_RNG_DISCIPLINE`), and it has already caused a real regression in this repo.
#   3. THE ARTIFACT IS PERSISTED **BEFORE** ANY QUALITY VERDICT IS PRINTED. A crash between
#      training and reporting must not lose hours of CPU, and a verdict must not be able to
#      retroactively shape what was recorded (the persist-before-verdict ordering plan 13-10
#      used for the tau probe).
#
# CPU-ONLY, NO NEW DEPENDENCY (CLAUDE.md, D-01, D-10). `P13_USE_GPU = false` is pre-registered,
# `train_three_way` hard-throws on `use_gpu = true`, CUDA is never imported, and `savepath` is
# passed EXPLICITLY so NeuralEstimators' loss history lands in a named artifact directory rather
# than in `tempdir()` (its default), where the run's own history would be indistinguishable from
# any other run's.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. Reaches `src/` only transitively and
# READ-ONLY through the datagen chain. `src/`, `spike/Project.toml` and `spike/Manifest.toml`
# stay byte-unchanged; nothing is installed.

using Statistics         # mean (the realized-frequency reporting)
using Dates              # UTC-labelled artifact timestamp
using Random             # Random.seed! -- layer 2 of the RNG discipline
using SHA                # stdlib: the git blob sha of the frozen pre-registration
import Flux              # Flux.state, via save_three_way

# ORDER MATTERS: the datagen surface first (it carries the whole Phase-11 precondition preamble
# through preconditions.jl, plus the pre-registration, the labels and the net). Guarded for
# idempotency, the house guarded-include idiom.
isdefined(@__MODULE__, :p13_generate_pool)   || include(joinpath(@__DIR__, "datagen.jl"))
isdefined(@__MODULE__, :ThreeWayEvidenceNet) || include(joinpath(@__DIR__, "net.jl"))

if !isdefined(@__MODULE__, :P13_NET_PATH)
    "Default on-disk path of the trained three-way evidence net."
    const P13_NET_PATH = joinpath(@__DIR__, "three_way_net.jld2")

    """
        P13_TRAIN_SAVEPATH

    Explicit `savepath` for NeuralEstimators' own training side effects. `train` otherwise
    defaults it to `tempdir()` and UNCONDITIONALLY writes `best_optimizer.bson`,
    `final_optimizer.bson` and `loss_per_epoch.csv` there, so the run's loss history would be
    anonymous and unattributable. Naming it makes the history an artifact of THIS run.
    """
    const P13_TRAIN_SAVEPATH = joinpath(@__DIR__, "three_way_train")
end

"""
    p13_git_blob_sha(path) -> String

The GIT BLOB sha1 of `path`: `sha1("blob " * filesize * "\\0" * contents)`, i.e. exactly what
`git hash-object <path>` prints.

WHY A GIT BLOB SHA AND NOT ONLY A PLAIN CONTENT HASH. `save_three_way` already records the
sha256 of `spike/p13/consts.jl`, which pins the bytes. The git blob sha additionally lets a
reader run `git cat-file -p <sha>` and get the exact pre-registration TEXT out of history, or
`git log --find-object=<sha>` and get the commit that introduced it -- so the artifact points
into the audit trail rather than merely proving the bytes were some particular bytes.
"""
function p13_git_blob_sha(path)
    bytes = read(path)
    ctx   = SHA.SHA1_CTX()
    SHA.update!(ctx, Vector{UInt8}("blob $(length(bytes))\0"))
    SHA.update!(ctx, bytes)
    return bytes2hex(SHA.digest!(ctx))
end

"""
    p13_read_risk_history(savepath) -> Union{Matrix{Float64},Nothing}

Read NeuralEstimators' `loss_per_epoch.csv` out of `savepath`: column 1 is the training risk and
column 2 the validation risk, and ROW 1 IS THE PRE-TRAINING INITIAL VALIDATION RISK DUPLICATED
INTO BOTH COLUMNS (`train.jl:544`), not epoch 1. So the best-checkpoint EPOCH is
`argmin(column 2) - 1`, and an `argmin` of 1 means no epoch ever beat the initial risk.

Parsed by hand rather than through `NeuralEstimators.loadrisk` so a runner does not pull
`DataFrames` into its load path for two columns of floats. Returns `nothing` when the file is
absent, so a caller can report "not captured" instead of crashing after a successful train.
"""
function p13_read_risk_history(savepath)
    path = joinpath(savepath, "loss_per_epoch.csv")
    isfile(path) || return nothing
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        vals = [parse(Float64, strip(x)) for x in split(s, ',')]
        push!(rows, vals)
    end
    isempty(rows) && return nothing
    return reduce(vcat, (permutedims(r) for r in rows))
end

"""
    p13_pi_head_log_odds(; n, tau) -> NamedTuple

The pi-LEVEL per-head log-odds, MEASURED from `n` cheap prior class draws (no simulation --
`prior_class_draws` needs only two `sample_prior` calls per draw, because the D-05 label is a
function of theta alone).

PRINTED BESIDE THE TRAINING-SUBSET VALUES SO A LARGE DISCREPANCY IS IMMEDIATELY VISIBLE. The two
numbers answer different questions and BOTH are wanted in the run log:

  * the pi-level pair says how imbalanced the three hypotheses are UNDER THE PRIOR. 13-RESEARCH
    E1's measured masses at tau = 0.15 are E 0.3119 / R 0.4713 / C 0.2168, so the pi-level
    log-odds are about log(0.2168/0.4713) = -0.776 and log(0.3119/0.4713) = -0.413 -- and the
    shipped BINARY net's analogous term is -0.0102, i.e. these are 40-75x larger. That is the
    concrete reason Pitfall 1 forbids copying the binary correction by analogy.
  * the TRAINING-SUBSET pair is what `three_way_log_bf` actually subtracts, and under
    `P13_TARGET_CLASS_FREQ` (equal thirds) it is small BY CONSTRUCTION -- stratification is what
    removes the imbalance. It is still MEASURED rather than set to zero, because the contiguous
    train/val split does not divide a stratified pool into exactly equal thirds.

The draws ride `p13_rng(P13_DATAGEN_COUNTER)` (`prior_class_draws`' own default), whose Philox KEY
WORD is `P13_DEV_SEED xor P13_SALT` -- provably DIFFERENT from the pool's
`P13_DEV_SEED xor P13_SALT xor P13_DATAGEN_COUNTER`, so this diagnostic cannot consume any part
of the stream the reported pool was drawn from.
"""
function p13_pi_head_log_odds(; n::Integer = 200_000, tau::Real = p13_tau())
    classes = prior_class_draws(n; tau = tau)
    return (masses = class_masses(classes), log_odds = head_log_odds(classes), n = Int(n))
end

"""
    main(; n = P13_TRAIN_N, net_path = P13_NET_PATH, verbose = true) -> NamedTuple

The reported Phase-13 training run, end to end and once.
"""
function main(; n::Integer = P13_TRAIN_N, net_path = P13_NET_PATH,
              savepath = P13_TRAIN_SAVEPATH, verbose::Bool = true)

    # --- 1. THE HARD PRECONDITION, FIRST, BEFORE ANYTHING EXPENSIVE ------------------------
    # D-02: a missing Phase-11 research NPE is a BLOCK on the phase, not a skip and not a
    # fallback to the shipped grid-8 basis. `load_p13_basis` routes through
    # `p13_require_phase11`, which asserts the whole six-item contract and throws otherwise.
    p13_require_phase11()
    basis = load_p13_basis()

    # The grid is DERIVED from the frozen standardizer's own continuous-row count, never typed.
    G         = isqrt(length(basis.zt.mean))
    n_cond    = p13_conditioning_length()
    input_dim = three_way_input_dim(G, n_cond)

    # --- 2. THE LOCKED-RECIPE BANNER -------------------------------------------------------
    # Every value is ECHOED BY NAME from the pre-registration, so the run log itself is the
    # evidence that nothing was tuned. Nothing here is a literal in this file.
    if verbose
        println("="^78)
        println("PHASE-13 THREE-WAY EVIDENCE NET -- REPORTED TRAINING RUN (D-07, D-09, D-10, D-11)")
        println("="^78)
        println("RECIPE (copied verbatim from the shipped binary trainer; NOT tuned):")
        println("  P13_TRAIN_N         = $P13_TRAIN_N        (requested n = $n)")
        println("  P13_EPOCHS          = $P13_EPOCHS")
        println("  P13_BATCHSIZE       = $P13_BATCHSIZE")
        println("  P13_LR              = $P13_LR")
        println("  P13_WEIGHT_DECAY    = $P13_WEIGHT_DECAY")
        println("  P13_VAL_FRAC        = $P13_VAL_FRAC")
        println("  P13_STOPPING_EPOCHS = $P13_STOPPING_EPOCHS")
        println("  P13_SUMMARY_WIDTH   = $P13_SUMMARY_WIDTH")
        println("  P13_NUM_SUMMARIES   = $P13_NUM_SUMMARIES")
        println("  P13_USE_GPU         = $P13_USE_GPU        (CPU-only is the reproducible baseline)")
        println("HYPOTHESIS BOUNDARY (Tier 2, MEASURED not chosen):")
        println("  P13_TAU                    = $P13_TAU")
        println("  P13_TAU_REFERENCE_LAMBDA   = $P13_TAU_REFERENCE_LAMBDA  " *
                "(tau is tau(lambda) -- never quote tau without it)")
        println("  P13_TAU_MEASURED_AUC       = $P13_TAU_MEASURED_AUC  (bar $P13_TAU_AUC)")
        println("  P13_TAU_PROBE_SHA          = $P13_TAU_PROBE_SHA")
        println("  P13_TAU_SIMULATOR_SHA      = $P13_TAU_SIMULATOR_SHA")
        println("  P13_TAU_PROBE_ARTIFACT     = $P13_TAU_PROBE_ARTIFACT")
        println("  P13_CUT_VARIANT            = $P13_CUT_VARIANT")
        println("STRATIFICATION (D-07-i):")
        println("  P13_STRATIFICATION         = $P13_STRATIFICATION   (WHOLE CLASS, never within)")
        println("  P13_TARGET_CLASS_FREQ      = $P13_TARGET_CLASS_FREQ")
        println("  P13_ITERATION_ALLOWANCE    = $P13_ITERATION_ALLOWANCE  " *
                "(reserved for the stratification switch; this run does NOT spend it)")
        println("BASIS (D-02) and INPUT WIDTH (D-03, derived never typed):")
        println("  Phase-11 artifact   = $P13_PHASE11_NET")
        println("  frozen zt           = $(length(basis.zt.mean)) continuous rows, INHERITED")
        println("  G                   = $G   n_cond = $n_cond")
        println("  input width         = three_way_input_dim($G, $n_cond) = $input_dim")
        println("  lambda range        = $P13_PHASE11_LAMBDA_RANGE")
        println("STREAM (D-01):")
        println("  P13_DEV_SEED        = $(repr(UInt64(P13_DEV_SEED))) " *
                "at counter P13_DATAGEN_COUNTER = $P13_DATAGEN_COUNTER")
        println("="^78)
    end

    # Stream hygiene, asserted rather than trusted: the reported seed must not be any seed a
    # prior lane has already observed, and must not be the fixture seed (a fixture stream that
    # governed a reported number would have been pre-observed by the test suite).
    @assert !(UInt64(P13_DEV_SEED) in _p13_forbidden()) "P13_DEV_SEED is a forbidden (pre-observed) seed"
    @assert UInt64(P13_DEV_SEED) != UInt64(P13_FIX_SEED) "the reported stream must differ from the fixture stream"
    @assert P13_USE_GPU == false "P13_USE_GPU must be false: CPU-only is the reproducible baseline"

    # --- 3. THE POOL: load, or generate through the stratifier if absent -------------------
    dir  = p13_generate_pool(n; basis = basis, verbose = verbose)
    pool = p13_load_pool_dir(dir)
    length(pool.lambda) == n || error(
        "main: the pool at $dir holds $(length(pool.lambda)) items but $n were requested")
    items    = p13_pool_to_items(pool)
    realized = class_masses(ThreeWayClass.(pool.class))

    if verbose
        println("-"^78)
        println("REALIZED CLASS FREQUENCIES vs TARGET (D-07-i)")
        for k in (:exclusion, :random, :coloc)
            r = Float64(getproperty(realized, k))
            t = Float64(getproperty(P13_TARGET_CLASS_FREQ, k))
            println("  $(rpad(k, 10)) realized $(round(r; digits = 6))   " *
                    "target $(round(t; digits = 6))   delta $(round(r - t; digits = 6))")
        end
    end

    # The Pitfall-5 warning sign, MEASURED and reported with an explicit verdict line. Read
    # `check_frozen_zt`'s docstring: the HARD half (zt object identity against the memoized
    # single load) has already passed by the time this returns at all, and the SOFT moment
    # heuristic has limited power here BY DESIGN OF D-02, because the Phase-13 pool is drawn
    # from the SAME joint the frozen transform was fitted on.
    ztchk = check_frozen_zt(items, basis)
    if verbose
        tag = ztchk.verdict === :refit_suspected ? "WARN" : "PASS"
        println("-"^78)
        println("FROZEN-zt CHECK (Pitfall 5) [$tag]: verdict = $(ztchk.verdict)")
        println("  max|per-row mean| = $(ztchk.max_abs_mean)  (diagnostic tol $P13_FROZEN_ZT_MEAN_TOL)")
        println("  max|per-row sd-1| = $(ztchk.max_abs_sd_dev)  (diagnostic tol $P13_FROZEN_ZT_SD_TOL)")
        println("  over $(ztchk.n_cols) standardized columns x $(ztchk.n_rows) continuous rows")
        if ztchk.verdict === :refit_suspected
            println("  WARN is REPORTED, NOT SUPPRESSED, and it is NOT evidence of a re-fit here:")
            println("  the SOFT heuristic assumes the Phase-13 pool is a DIFFERENT distribution")
            println("  from the one zt was fitted on, but D-02 requires it to be the SAME joint")
            println("  (F5 train-joint == eval-joint), so centred/unit-scaled moments are the")
            println("  EXPECTED outcome of a correctly INHERITED transform. The unambiguous")
            println("  evidence is the HARD zt object-identity check, which passed, plus the")
            println("  source assertion that no z-score transform is constructed in the")
            println("  Phase-13 data path (spike/test/test_p13_datagen.jl).")
        end
    end

    # --- 4. THE CONTIGUOUS TRAIN/VAL SPLIT -------------------------------------------------
    # Contiguous, exactly as the analog does it (`spike/validation/train_ratio.jl:296-302`:
    # "the stream is already i.i.d. across columns"). THAT PREMISE IS TRUE HERE ONLY BECAUSE
    # `p13_select_indices` DETERMINISTICALLY PERMUTES THE ACCEPTED ORDER: whole-class
    # accept/reject fills the classes at different times, so the raw accepted list ends in a
    # nearly pure COLOC tail and the contiguous validation slice measured E/R/C = 1266/0/5934
    # -- ZERO RANDOM, which is the shared negative of BOTH heads. See that function's docstring.
    assembled = assemble_conditioned_pairs(items)
    Z         = assembled.Z
    targets   = assembled.targets
    classes   = assembled.classes

    nval = clamp(round(Int, P13_VAL_FRAC * n), 1, n - 1)
    ntr  = n - nval
    Z_tr, Z_va             = Z[:, 1:ntr], Z[:, ntr+1:end]
    targets_tr, targets_va = targets[:, 1:ntr], targets[:, ntr+1:end]
    classes_tr             = classes[1:ntr]
    classes_va             = classes[ntr+1:end]

    if verbose
        println("-"^78)
        println("SPLIT (contiguous at P13_VAL_FRAC = $P13_VAL_FRAC): train $ntr / val $nval")
        println("  train class counts  E=$(count(==(EXCLUSION), classes_tr)) " *
                "R=$(count(==(RANDOM), classes_tr)) C=$(count(==(COLOC), classes_tr))")
        println("  val   class counts  E=$(count(==(EXCLUSION), classes_va)) " *
                "R=$(count(==(RANDOM), classes_va)) C=$(count(==(COLOC), classes_va))")
        println("  Z_tr = $(size(Z_tr))  Z_va = $(size(Z_va))  " *
                "targets_tr = $(size(targets_tr))  targets_va = $(size(targets_va))")
    end
    # RANDOM is the shared negative class of BOTH heads (D-11), so a validation split without it
    # would make the validation risk -- which drives early stopping AND the best-checkpoint rule
    # -- minimisable by predicting positive everywhere. This is the guard for that.
    count(==(RANDOM), classes_va) > 0 || error(
        "main: the validation split contains NO RANDOM items. RANDOM is the shared negative " *
        "class of both heads, so early stopping and best-checkpoint selection would be driven " *
        "by a degenerate signal. The accepted order must be permuted (p13_select_indices).")

    # --- 5. THE PER-HEAD CORRECTION, MEASURED ON THE **TRAINING** SUBSET ONLY (D-07) -------
    # Never on the evaluation set (that would tune the reported evidence scale to the data it is
    # scored against) and never on the whole pool (that would count columns the net never saw).
    hlo = head_log_odds(classes_tr)
    isfinite(hlo.coloc) && isfinite(hlo.exclusion) || error(
        "main: measured head log-odds are not finite ($(hlo)) -- degenerate head balance in the " *
        "training subset")
    piref = p13_pi_head_log_odds(; tau = p13_tau())
    if verbose
        println("-"^78)
        println("PER-HEAD EVIDENCE-SCALE CORRECTION (D-07), MEASURED ON THE TRAINING SUBSET")
        println("  head_log_odds.coloc     = $(hlo.coloc)      " *
                "[log(n_COLOC / n_RANDOM) over the training subset]")
        println("  head_log_odds.exclusion = $(hlo.exclusion)  " *
                "[log(n_EXCLUSION / n_RANDOM) over the training subset]")
        println("  pi-level reference (measured from $(piref.n) prior class draws, NO simulation):")
        println("    pi masses          = $(piref.masses)")
        println("    pi log-odds coloc  = $(piref.log_odds.coloc)")
        println("    pi log-odds exclus = $(piref.log_odds.exclusion)")
        println("  The training-subset pair is SMALL BY CONSTRUCTION: P13_TARGET_CLASS_FREQ is")
        println("  equal thirds, so stratification removes the pi-level imbalance and the")
        println("  correction goes to ~0. IT IS STILL MEASURED, NOT ASSUMED -- the contiguous")
        println("  split does not divide a stratified pool into exactly equal thirds, and D-07")
        println("  makes the measurement the contract rather than the number.")
    end

    # --- 6. THE NET ------------------------------------------------------------------------
    net = build_three_way_net(input_dim)
    verbose && println("-"^78)
    verbose && println("NET: $(three_way_arch(net))  (trunk copied from the shipped binary trainer)")

    # --- 7. TRAIN. ONE RUN. -----------------------------------------------------------------
    # LAYER 2 OF THE RNG DISCIPLINE, IMMEDIATELY BEFORE THE TRAIN CALL. Flux's shuffling
    # `DataLoader` (NeuralEstimators `_DataLoader`, `train.jl:668-676`, `shuffle = true`) draws
    # from the GLOBAL RNG, so the counter-based streams threaded through generation are NOT
    # enough on their own (Pitfall 6; `P13_GLOBAL_RNG_DISCIPLINE` layer 2). `train_three_way`
    # re-seeds identically from its own `seed` keyword just before it calls `train`, so the two
    # seedings agree by construction rather than competing; this one additionally pins anything
    # the pool assembly above might have drawn from the global RNG.
    mkpath(savepath)
    verbose && println("TRAINING (savepath = $savepath; loss history captured as an artifact)")
    t_train = time()
    Random.seed!(UInt64(P13_DEV_SEED))
    trained = train_three_way(net, targets_tr, targets_va, Z_tr, Z_va;
                              savepath = savepath, verbose = verbose)
    train_min = (time() - t_train) / 60

    # --- 8. PERSIST BEFORE ANY VERDICT ------------------------------------------------------
    consts_path = joinpath(@__DIR__, "consts.jl")
    save_three_way(net_path, trained, hlo;
        meta = (
            n                    = n,
            n_train              = ntr,
            n_val                = nval,
            epochs               = P13_EPOCHS,
            batchsize            = P13_BATCHSIZE,
            learning_rate        = P13_LR,
            weight_decay         = P13_WEIGHT_DECAY,
            val_frac             = P13_VAL_FRAC,
            stopping_epochs      = P13_STOPPING_EPOCHS,
            summary_width        = P13_SUMMARY_WIDTH,
            num_summaries        = P13_NUM_SUMMARIES,
            use_gpu              = P13_USE_GPU,
            input_dim            = input_dim,
            grid                 = G,
            n_cond               = n_cond,
            tau                  = Float64(p13_tau()),
            tau_reference_lambda = P13_TAU_REFERENCE_LAMBDA,
            tau_measured_auc     = P13_TAU_MEASURED_AUC,
            tau_probe_sha        = P13_TAU_PROBE_SHA,
            tau_probe_artifact   = P13_TAU_PROBE_ARTIFACT,
            phase11_net          = P13_PHASE11_NET,
            phase11_sha          = p13_phase11_sha(),
            phase11_lambda_range = P13_PHASE11_LAMBDA_RANGE,
            zt_continuous_rows   = length(basis.zt.mean),
            zt_provenance        = "INHERITED from the Phase-11 research NPE, never re-fit (D-02)",
            frozen_zt_verdict    = ztchk.verdict,
            frozen_zt_max_mean   = Float64(ztchk.max_abs_mean),
            frozen_zt_max_sd_dev = Float64(ztchk.max_abs_sd_dev),
            realized_class_freq  = realized,
            target_class_freq    = P13_TARGET_CLASS_FREQ,
            stratification       = P13_STRATIFICATION,
            head_log_odds_scope  = "MEASURED on the TRAINING subset only (D-07)",
            pi_class_masses      = piref.masses,
            pi_head_log_odds     = piref.log_odds,
            pool_dir             = dir,
            master_seed          = UInt64(P13_DEV_SEED),
            salt                 = UInt64(P13_SALT),
            datagen_counter      = Int(P13_DATAGEN_COUNTER),
            consts_git_blob_sha  = p13_git_blob_sha(consts_path),
            iteration_allowance  = P13_ITERATION_ALLOWANCE,
            train_minutes        = train_min,
            savepath             = savepath,
            generated            = string(Dates.now(Dates.UTC)) * "Z",
        ))
    verbose && println("PERSISTED (before any verdict): $net_path")

    # --- 9. THE CLOSING SUMMARY -------------------------------------------------------------
    risk = p13_read_risk_history(savepath)
    final_train = final_val = best_val = NaN
    best_epoch  = nothing
    if risk !== nothing && size(risk, 1) >= 1
        final_train = risk[end, 1]
        final_val   = risk[end, 2]
        bi          = argmin(view(risk, :, 2))
        best_val    = risk[bi, 2]
        # Row 1 is the PRE-TRAINING initial validation risk duplicated into both columns
        # (train.jl:544), so epoch numbering starts at row 2.
        best_epoch  = bi - 1
    end
    if verbose
        println("="^78)
        println("RUN COMPLETE -- ONE RUN, recipe unchanged, D-04 allowance UNSPENT")
        println("  wall clock (training only) = $(round(train_min; digits = 2)) min")
        println("  epochs recorded            = " *
                "$(risk === nothing ? "not captured" : string(size(risk, 1) - 1)) " *
                "of a maximum $P13_EPOCHS (patience $P13_STOPPING_EPOCHS)")
        println("  final training risk        = $final_train")
        println("  final validation risk      = $final_val")
        println("  best validation risk       = $best_val at epoch $(best_epoch)")
        println("  head_log_odds persisted    = $hlo")
        println("  artifact                   = $net_path")
        println("  loss history               = $(joinpath(savepath, "loss_per_epoch.csv"))")
        println("="^78)
    end

    return (net = trained, head_log_odds = hlo, realized_class_freq = realized,
            frozen_zt = ztchk, pool_dir = dir, n_train = ntr, n_val = nval,
            final_train_risk = final_train, final_val_risk = final_val,
            best_val_risk = best_val, best_epoch = best_epoch,
            train_minutes = train_min, pi_reference = piref, net_path = net_path)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so `include`-ing this file can NEVER trigger a training run. Run it deliberately:
#
#     julia --project=spike -t auto spike/p13/train_three_way.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
