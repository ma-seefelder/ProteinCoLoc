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

# spike/p13/net.jl --- Phase-13 three-way evidence net: shared trunk, two BCE heads (D-08..D-11).
#
# (a) WHY THIS LIVES OUTSIDE THE SHIPPED BINARY ESTIMATOR TYPE.
#     NeuralEstimators' `RatioEstimator` is structurally binary for TWO independent reasons,
#     both read out of the pinned v0.2.1 source rather than from memory:
#       1. its inference network always ends in a SINGLE output unit
#          (RatioEstimator.jl:86), so even a 3-dim model index would score a 3-vector with one
#          logit rather than emit three; and
#       2. its loss is HARD-CODED -- `_loss(estimator::RatioEstimator, loss = nothing) =
#          logitbinarycrossentropy` (RatioEstimator.jl:124) -- so a loss PASSED to `train` is
#          SILENTLY IGNORED. No warning, no error: the concrete-type method simply wins over
#          the generic fallback. A two-head masked loss handed to that type would look
#          installed and would never run.
#     That is why D-09 calls this an architectural fork and not a tuning knob.
#
# (b) WHY ROUTE R1 (A `<: NeuralEstimator` SUBTYPE) IS THE ONE TAKEN.
#     A custom subtype takes the GENERIC hooks one file over: `_loss(estimator, loss) = loss`
#     (train.jl:654) and `_inputoutput(estimator, Z, theta) = (Z, theta)` (train.jl:657), and
#     the Flux extension calls the loss as `loss(trainstate.model(input), output)`
#     (NeuralEstimatorsFluxExt.jl:67, :80) -- a two-argument function of (network output,
#     target array), BOTH of which may be multi-row matrices. So the passed loss IS honoured,
#     and the optimiser, the Float64 CosAnneal LR schedule, the early-stopping patience and the
#     best-validation checkpoint rule stay byte-for-byte the ones the shipped binary net was
#     trained with. THAT is what makes D-10's attribution argument true: any difference in the
#     two-way overlap is attributable to the HEAD, not to a re-tuned recipe.
#
# (c) THE PRE-DECLARED FALLBACK. Route R2 -- a hand-rolled Flux epoch loop -- is permitted ONLY
#     if the A5 smoke test (spike/test/test_p13_net.jl, "A5 smoke") fails. It would have to
#     re-implement CosAnneal, patience and best-checkpointing IDENTICALLY or D-10's attribution
#     claim is no longer true, so if it is ever used it must be recorded as a DECLARED
#     deviation and never silently substituted. Route R3 -- patching or forking the shipped
#     estimator type -- is rejected outright by D-09.
#
# (d) TWO THINGS THAT MUST NOT BE COPIED FROM THE BINARY READ SURFACE.
#     1. Its subtracted prior-odds term. The shipped binary path's correction belongs to a
#        DIFFERENT construction: the shuffled-theta joint-vs-marginal separation of
#        RatioEstimator.jl:94-122 already CONTAINS the ratio, so that term is a numerically
#        negligible (~0.01 nat) redundancy over ONE denominator spanning ALL labels. Here there
#        are TWO scalars, each counted over its own head's restricted class pair, and each is
#        15-50x larger. Because the binary term is ~0 a wrong copy is INVISIBLE in the file it
#        was copied from; it is not invisible here. The correction consumed below is
#        `head_log_odds` from spike/p13/labels.jl, derived from scratch in 13-RESEARCH F1/F4.
#     2. A three-class softmax head. Its logit differences are POSTERIOR odds over three
#        classes and would need a model prior over all three to become a Bayes factor --
#        exactly what D-08 refuses to invent. Two independent binary heads give two independent
#        two-model comparisons, each with its own well-defined correction.
#
# (e) THIS FILE TRAINS NOTHING. The script entry point at the foot is guarded and prints a
#     usage note; the only training this surface ever does at gate time is the 200-sample,
#     two-epoch A5 smoke in the test file. The real training run is plan 13-11's.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. Touches no src/ byte, adds no
# dependency, installs nothing, and never imports CUDA. Guarded includes keep the file loadable
# standalone AND idempotent under a runner.

import NeuralEstimators: NeuralEstimator, train   # the subtype supertype + the inherited train loop
import Flux                                        # Flux.@layer, Flux.Optimisers.AdamW, Flux.state
import Flux: Chain, Dense, gelu,                   # the D-10 trunk, lifted verbatim
             logitbinarycrossentropy,              # the per-cell loss of both heads
             sigmoid                                # the UNCORRECTED head probabilities (D-13)
import Random                                      # Random.seed! -- the global-RNG discipline
using JLD2                                         # atomic model persistence (pure Julia)
using SHA                                          # stdlib: provenance sha of the frozen consts
using Dates                                        # UTC-labelled artifact timestamp

# ORDER MATTERS: the pre-registration first (P13_NUM_SUMMARIES, P13_SUMMARY_WIDTH, the copied
# recipe, P13_USE_GPU, P13_DEV_SEED), then the label surface (target_matrix, head_log_odds).
# Guarded for idempotency, the house guarded-include idiom.
isdefined(@__MODULE__, :P13_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :three_way_label) || include(joinpath(@__DIR__, "labels.jl"))

if !isdefined(@__MODULE__, :ThreeWayEvidenceNet)

# Persistence schema for the three-way net artifact (bump on any breaking layout change).
const P13_NET_SCHEMA = 1

# --- The DERIVED input width ----------------------------------------------------------------
#
# NEVER WRITE THE INPUT WIDTH AS A LITERAL. Phase 11 may deliver a conditioning vector whose
# length is not 1, so the width is DERIVED as `p13_ratio_input_dim(G) + n_cond` with
# `n_cond = length(encode_lambda(lambda))` supplied by the caller. A hard-coded width would be
# silently wrong the moment the conditioning length changed, and nothing would throw.
# `P13_INPUT_WIDTH_RULE = :ratio_input_dim_plus_ncond` records the rule in the pre-registration.
"""
    p13_ratio_input_dim(G::Integer) -> Int

Input dimension of the paired-summary encoding for a `G x G` patch grid: `5*G^2`.

ATTRIBUTED COPY of `ratio_input_dim` (`src/amortized/summary.jl:127`), reproduced here rather
than reached through an `include` because `src/amortized/summary.jl` is not loadable standalone
(its method signatures mention `MultiChannelImage`, so including it would drag the whole frozen
image-loading chain into a file that needs one arithmetic identity). `src/` stays byte-unchanged.

The identity: the A7 difference encoding concatenates the two `2*G^2` summaries plus their
`G^2`-dim continuous-row contrast, i.e. `2*(2*G^2) + G^2 = 5*G^2`.
"""
function p13_ratio_input_dim(G::Integer)
    G >= 1 || throw(ArgumentError("p13_ratio_input_dim: G must be >= 1 (got $G)"))
    return 5 * G^2
end

"""
    three_way_input_dim(G::Integer, n_cond::Integer) -> Int

The three-way evidence net's input width: `p13_ratio_input_dim(G) + n_cond`, where `n_cond` is
the length of the registration-uncertainty conditioning vector appended AFTER the pair encoding
(D-03; `P13_LAMBDA_PLACEMENT = :append_after_pair_encode`).

DERIVED, NEVER A LITERAL -- see the block comment above this function. `n_cond = 0` is allowed
so an unconditioned net can be built for a toy or an ablation; a negative `n_cond` throws.
"""
function three_way_input_dim(G::Integer, n_cond::Integer)
    n_cond >= 0 || throw(ArgumentError("three_way_input_dim: n_cond must be >= 0 (got $n_cond)"))
    return p13_ratio_input_dim(G) + n_cond
end

# --- The architecture (D-09, D-10) ----------------------------------------------------------
"""
    ThreeWayEvidenceNet(trunk, head)

The D-09 architecture: ONE shared Flux trunk feeding a two-logit head.

    net(Z) -> 2 x n logit matrix, produced in ONE forward pass

Row 1 is the coloc-vs-random logit and row 2 the exclusion-vs-random logit. The single forward
pass is the surviving, bindable content of SC1: one network, one pass, three-way evidence
against a random reference.

A single `Dense(num_summaries, 2)` is used rather than two `Dense(num_summaries, 1)` layers.
The two are mathematically identical -- the rows do not interact -- and the single layer keeps
the one-pass claim literal in the source instead of leaving it to a reader's inference.

Subtyping `NeuralEstimator` is what buys the inherited training loop AND an honoured custom
loss (see this file's banner (b)); `Flux.@layer` is REQUIRED so `Optimisers.setup` can traverse
the struct when the training state is constructed.
"""
struct ThreeWayEvidenceNet <: NeuralEstimator
    trunk
    head
end

Flux.@layer ThreeWayEvidenceNet

(m::ThreeWayEvidenceNet)(Z) = m.head(m.trunk(Z))

"""
    build_three_way_net(input_dim; num_summaries = P13_NUM_SUMMARIES,
                        width = P13_SUMMARY_WIDTH) -> ThreeWayEvidenceNet

Build the D-10 trunk plus the D-09 two-logit head:

    Chain(Dense(input_dim, W, gelu), Dense(W, W, gelu), Dense(W, W, gelu), Dense(W, num_summaries))

with `W = P13_SUMMARY_WIDTH = 256` and `num_summaries = P13_NUM_SUMMARIES = 64`.

EVERY CONSTANT HERE IS COPIED, NOT CHOSEN (D-10). The four `Dense` lines, the two widths and
the `input_dim >= 1` guard are lifted from the shipped binary trunk builder
(`src/amortized/train_ratio.jl:72-80`); only the estimator wrapper and the head differ. Copying
rather than re-tuning is the whole basis of the attribution argument: a difference in the
two-way overlap is then attributable to the head alone.

`input_dim` comes from `three_way_input_dim(G, n_cond)` and is never a literal. SHARED BY THE
TRAINER AND THE LOADER, so training and reload build the identical topology.
"""
function build_three_way_net(input_dim::Integer;
                             num_summaries::Integer = P13_NUM_SUMMARIES,
                             width::Integer = P13_SUMMARY_WIDTH)
    input_dim >= 1 || throw(ArgumentError("build_three_way_net: input_dim must be >= 1"))
    W = width
    trunk = Chain(Dense(input_dim, W, gelu), Dense(W, W, gelu),
                  Dense(W, W, gelu), Dense(W, num_summaries))
    return ThreeWayEvidenceNet(trunk, Dense(num_summaries, 2))
end

"""
    three_way_arch(net::ThreeWayEvidenceNet) -> NamedTuple

The architecture metadata `(input_dim, width, num_summaries)` READ BACK OFF THE LAYERS rather
than carried alongside them, so a persisted artifact can never disagree with the net it stores.
`Dense` weights are `(out, in)`, hence `input_dim = size(first weight, 2)`.
"""
function three_way_arch(net::ThreeWayEvidenceNet)
    W1 = net.trunk.layers[1].weight
    Wl = net.trunk.layers[end].weight
    return (input_dim = size(W1, 2), width = size(W1, 1), num_summaries = size(Wl, 1))
end

# --- The D-11 masked two-head loss ----------------------------------------------------------
"""
    masked_two_head_bce(yhat, y) -> Float32

The D-11 masked two-term logit-BCE. `yhat` is the `2 x n` logit matrix the net returns; `y` is
the `4 x n` target matrix `target_matrix` builds, with row legend

    row 1 : y_C   BCE target of the coloc head       (1 = coloc, 0 = random)
    row 2 : w_C   participation weight of that head  (1 iff the sample is coloc OR random)
    row 3 : y_E   BCE target of the exclusion head   (1 = exclusion, 0 = random)
    row 4 : w_E   participation weight of that head  (1 iff the sample is exclusion OR random)

The denominator is the number of ACTIVE (sample, head) cells -- `sum(w_C) + sum(w_E)` -- not
`n` and not `2n`, so head imbalance cannot silently reweight one comparison against the other.

WHY THE RESTRICTION IS THE POINT (D-11). Each head trains ONLY on its own class pair: a coloc
sample switches the exclusion head OFF and vice versa, while the random class is the negative
of BOTH heads. That shared negative is what makes the two logits ratios against the SAME
reference, which is D-08's premise. Under one-vs-rest (`w_C == 1` everywhere) the coloc head's
negative class would be the MIXTURE of random and exclusion and its logit a MIXTURE Bayes
factor -- a different quantity published under D-08's name, by code that looks identical. The
TRUNK still sees every sample through at least one active head, so it trains jointly on all
three classes, which is what makes the two evidences structurally consistent (D-09).

A PASSED LOSS IS ACTUALLY HONOURED HERE. `ThreeWayEvidenceNet` is not the shipped binary
estimator type, so it takes the generic `_loss(estimator, loss) = loss` fallback
(`NeuralEstimators/src/train.jl:654`) instead of the concrete-type method that discards the
argument. The inactive cells are multiplied by an exact `0`, so the loss is BIT-IDENTICAL under
any finite perturbation of an inactive logit -- asserted with exact equality in the test suite,
because that assertion is what distinguishes this loss from one-vs-rest.
"""
function masked_two_head_bce(yhat, y)
    lC = logitbinarycrossentropy.(view(yhat, 1, :), view(y, 1, :); agg = identity)
    lE = logitbinarycrossentropy.(view(yhat, 2, :), view(y, 3, :); agg = identity)
    wC = view(y, 2, :)
    wE = view(y, 4, :)
    return (sum(wC .* lC) + sum(wE .* lE)) / (sum(wC) + sum(wE) + eps(Float32))
end

# --- Training (the recipe is COPIED, D-10) --------------------------------------------------
"""
    train_three_way(net, targets_tr, targets_va, Z_tr, Z_va; kwargs...) -> ThreeWayEvidenceNet

Train the two-head evidence net through NeuralEstimators' OWN `train` loop (route R1) on fixed
data: `Z_*` are `input_dim x n` summary matrices and `targets_*` the matching `4 x n` target
matrices from `target_matrix`.

THE RECIPE IS COPIED, NOT CHOSEN (D-10): `P13_EPOCHS`, `P13_BATCHSIZE`, `P13_LR`,
`P13_WEIGHT_DECAY` and `P13_STOPPING_EPOCHS` are the pre-registration's verbatim transcription
of the shipped binary trainer's defaults (`src/amortized/train_ratio.jl:155-163`), and the
AdamW arguments are deliberately left Float64 to match NeuralEstimators' Float64 CosAnneal
`lr_schedule` -- an `Optimisers.adjust!` gotcha, not decoration.

`loss` defaults to `masked_two_head_bce` and IS passed through (legal here, see that function's
docstring). Overriding it is what the suite's two-different-losses check exercises.

`Random.seed!` derived from `seed` is called IMMEDIATELY BEFORE `train` because Flux's shuffling
`DataLoader` -- which NeuralEstimators' `_DataLoader` builds with `shuffle = true`
(`train.jl:668-676`) -- draws from the GLOBAL RNG. Threading a counter-based stream is not
enough on its own; this has already caused a real reproducibility regression in this repository
(`P13_GLOBAL_RNG_DISCIPLINE`, layer 2).

`savepath` defaults to `nothing` so a run is side-effect-free. NeuralEstimators otherwise
defaults it to `tempdir()` and UNCONDITIONALLY writes `best_optimizer.bson`,
`final_optimizer.bson` and `loss_per_epoch.csv` there; pass an explicit directory to capture the
loss history as an artifact.

`use_gpu = true` HARD-THROWS: CPU-only is the reproducible baseline (CLAUDE.md, D-01), CUDA is
never imported, and `P13_USE_GPU = false` is pre-registered.

Errors with an `ArgumentError` when a target matrix does not have exactly 4 rows or when a
summary matrix and its target matrix disagree on the column count (shape confusion at the input
boundary is the failure that would otherwise train a head against the wrong labels).
"""
function train_three_way(net::ThreeWayEvidenceNet,
                         targets_tr::AbstractMatrix, targets_va::AbstractMatrix,
                         Z_tr::AbstractMatrix, Z_va::AbstractMatrix;
                         epochs::Integer = P13_EPOCHS,
                         batchsize::Integer = P13_BATCHSIZE,
                         learning_rate::Real = P13_LR,
                         weight_decay::Real = P13_WEIGHT_DECAY,
                         stopping_epochs::Integer = P13_STOPPING_EPOCHS,
                         use_gpu::Bool = P13_USE_GPU,
                         seed = P13_DEV_SEED,
                         loss = masked_two_head_bce,
                         savepath = nothing,
                         verbose::Bool = true)
    use_gpu && throw(ArgumentError(
        "train_three_way: use_gpu=true out of scope (CPU-only gate, D-01/D-10)"))

    size(targets_tr, 1) == 4 || throw(ArgumentError(
        "train_three_way: targets_tr must have 4 rows [y_C; w_C; y_E; w_E] " *
        "(got $(size(targets_tr, 1)))"))
    size(targets_va, 1) == 4 || throw(ArgumentError(
        "train_three_way: targets_va must have 4 rows [y_C; w_C; y_E; w_E] " *
        "(got $(size(targets_va, 1)))"))
    size(Z_tr, 2) == size(targets_tr, 2) || throw(ArgumentError(
        "train_three_way: Z_tr has $(size(Z_tr, 2)) columns but targets_tr has " *
        "$(size(targets_tr, 2))"))
    size(Z_va, 2) == size(targets_va, 2) || throw(ArgumentError(
        "train_three_way: Z_va has $(size(Z_va, 2)) columns but targets_va has " *
        "$(size(targets_va, 2))"))

    # LAYER 2 OF THE RNG DISCIPLINE, IMMEDIATELY BEFORE `train` (see docstring).
    Random.seed!(UInt64(seed))

    # FIXED-DATA train form. AdamW args stay Float64 to match the Float64 CosAnneal schedule.
    return train(net, targets_tr, targets_va, Z_tr, Z_va;
                 epochs = epochs, batchsize = batchsize, use_gpu = false,
                 optimiser = Flux.Optimisers.AdamW(learning_rate, (0.9, 0.999), weight_decay),
                 stopping_epochs = stopping_epochs, loss = loss,
                 savepath = savepath, verbose = verbose)
end

# --- The D-08 read surface ------------------------------------------------------------------
"""
    three_way_log_bf(net, Z_pair, head_log_odds; use_gpu = false) -> NamedTuple

The three-way evidence read for ONE conditioned paired summary, in ONE forward pass:

    (coloc = logit_C - head_log_odds.coloc,
     exclusion = logit_E - head_log_odds.exclusion,
     random = 0.0)

`Z_pair` is a `three_way_input_dim(G, n_cond)`-vector or the equivalent single-column matrix.

THE RANDOM ENTRY IS IDENTICALLY `0.0` BY CONSTRUCTION (D-08), not by measurement: both reported
numbers are log Bayes factors AGAINST the random reference, so the reference against itself is
zero. That is also why no model prior over the three hypotheses has to be invented -- turning
evidence into a chosen hypothesis is a later phase's job.

WHAT IS SUBTRACTED, AND WHAT MUST NOT BE. The subtracted term is the PER-HEAD MEASURED
`head_log_odds` from `spike/p13/labels.jl`, counted over that head's OWN class pair. For a plain
BCE head separating A from B at training frequencies q_A, q_B the Bayes-optimal logit is
`log[p_q(Z|A)/p_q(Z|B)] + log(q_A/q_B)`, so the log Bayes factor is the logit MINUS that
log-ratio. The binary read surface's own correction term is a DIFFERENT quantity over a
DIFFERENT denominator and must not be copied by analogy -- see this file's banner (d.1).

`use_gpu = true` HARD-THROWS (CPU-only, D-01); the keyword exists only to mirror the shipped
read surface's signature so the two are diffable.
"""
function three_way_log_bf(net::ThreeWayEvidenceNet, Z_pair, head_log_odds;
                          use_gpu::Bool = false)
    use_gpu && throw(ArgumentError(
        "three_way_log_bf: use_gpu=true out of scope (CPU-only gate, D-01/D-10)"))
    Zc = Z_pair isa AbstractVector ? reshape(Z_pair, :, 1) : Z_pair
    s  = net(Zc)                                   # 2 x 1 logits, ONE forward pass
    return (coloc     = Float64(s[1, 1]) - head_log_odds.coloc,
            exclusion = Float64(s[2, 1]) - head_log_odds.exclusion,
            random    = 0.0)
end

"""
    three_way_probs(net, Z_pair) -> NamedTuple

The two head probabilities `(coloc = sigmoid(logit_C), exclusion = sigmoid(logit_E))` from the
UNCORRECTED logits, in one forward pass.

DELIBERATELY UNCORRECTED. Calibration is a property of the classifier UNDER ITS OWN TRAINING
FREQUENCIES: the reliability of `q(A|Z)` is scored against the labels the head actually saw, so
subtracting the evidence-scale correction first would score a quantity that is not a probability
of anything. This is why D-13 gates ECE on these numbers while D-08 publishes the corrected log
Bayes factors from `three_way_log_bf`.
"""
function three_way_probs(net::ThreeWayEvidenceNet, Z_pair)
    Zc = Z_pair isa AbstractVector ? reshape(Z_pair, :, 1) : Z_pair
    s  = net(Zc)
    return (coloc = Float64(sigmoid(s[1, 1])), exclusion = Float64(sigmoid(s[2, 1])))
end

# --- The input encoding (D-03 lambda placement) ---------------------------------------------
"""
    p13_pair_encode(Zs, Zc) -> Vector{Float32}

The paired-summary encoding `vcat(Zs, Zc, Zs[1:G^2] - Zc[1:G^2])`: the two `2*G^2`-dim
frozen-`zt` summaries concatenated PLUS the `G^2`-dim continuous-row contrast. Length is
`p13_ratio_input_dim(G) = 5*G^2`, grid-general because the continuous-row count is DERIVED from
the summary length.

BYTE-FAITHFUL COPY of `pair_encode` (`src/amortized/bf.jl:91-100`), reproduced under a
`p13_` name rather than reached through an `include` because the spike-lane twin
(`spike/validation/train_ratio.jl:185`) pulls the whole Phase-5 validation harness -- a frozen
model load and a `RATIO_INPUT_DIM` literal -- into a file that needs neither, and because the
`src/` original is not loadable standalone. `src/` stays byte-unchanged.

THE `iseven` GUARD IS KEPT EXACTLY AS WRITTEN. It is the tripwire that catches the wrong
lambda placement: the wrong form `p13_pair_encode(Zs_129, Zc_129)` passes odd-length summaries
and throws here instead of silently producing a mis-shaped input. That loud failure is a gift;
do not weaken the guard to accommodate a caller.
"""
function p13_pair_encode(Zs::AbstractVector, Zc::AbstractVector)
    length(Zs) == length(Zc) ||
        throw(DimensionMismatch("p13_pair_encode: Zs ($(length(Zs))) and Zc ($(length(Zc))) " *
                                "must be the same length"))
    iseven(length(Zs)) ||
        throw(ArgumentError("p13_pair_encode: summary length $(length(Zs)) must be 2*G^2 (even)"))
    nc = length(Zs) ÷ 2                                   # continuous rows = G^2
    d  = @view(Zs[1:nc]) .- @view(Zc[1:nc])
    return Float32.(vcat(Zs, Zc, d))
end

"""
    encode_conditioned_pair(Zs, Zc, cond) -> Vector{Float32}

The net's input for one acquisition: `vcat(p13_pair_encode(Zs, Zc), Float32.(cond))`.

THE CORRECT FORM AND BOTH WRONG ONES (D-03, `P13_LAMBDA_PLACEMENT`):

    CORRECT     vcat(p13_pair_encode(Zs_2G2, Zc_2G2), cond)   # the conditioning goes LAST
    WRONG       p13_pair_encode(Zs_2G2p1, Zc_2G2p1)           # throws on the `iseven` guard
    ALSO WRONG  vcat(p13_pair_encode(Zs, Zc), cond, cond)     # conditioning counted twice

The first wrong form fails LOUDLY, which is why the `iseven` guard is worth relying on: the
conditioning row is appended AFTER standardization and after the pair encoding, so folding it
into the summaries would also push it through the contrast block. The second wrong form is
silent, which is why the realized length is asserted here against
`three_way_input_dim(G, length(cond))` rather than trusted.

`cond` may be any real iterable, including a length-1 vector for a scalar uncertainty level;
its length is whatever the upstream conditioning encoder returns and is never assumed to be 1.
"""
function encode_conditioned_pair(Zs::AbstractVector, Zc::AbstractVector, cond)
    c = Float32.(collect(cond))
    Z = vcat(p13_pair_encode(Zs, Zc), c)
    G = isqrt(length(Zs) ÷ 2)
    expected = three_way_input_dim(G, length(c))
    length(Z) == expected || error(
        "encode_conditioned_pair: built length $(length(Z)) but three_way_input_dim($G, " *
        "$(length(c))) is $expected -- the conditioning vector is misplaced or double-counted")
    return Z
end

# --- Persistence (atomic, narrow deserialization surface) -----------------------------------
"""
    save_three_way(path, net, head_log_odds; meta = (;)) -> String

Atomically persist the trained net: write to `path * ".tmp"`, REOPEN to integrity-check the
keys, then `mv(...; force = true)` (the filesystem-atomic idiom of `save_ratio`,
`spike/validation/train_ratio.jl:331-347`). A crash mid-write therefore leaves the previous
artifact intact and never a torn one; the suite asserts no `.tmp` survives a successful save.

THE NARROW DESERIALIZATION SURFACE IS DELIBERATE. What is stored is `Flux.state(net)` plus the
architecture metadata (`input_dim`, `num_summaries`, `width`) -- plain arrays and integers --
NOT the estimator object. `load_three_way` REBUILDS the topology with `build_three_way_net` and
loads the state into it, so reading an artifact never reconstructs an arbitrary object graph.

Also persisted for provenance: the schema version, the per-head correction, the tau in force
(`nothing` while tau is still unmeasured -- this function must not call the tau accessor, which
errors by design), and the sha256 of the frozen pre-registration file, so an artifact can be
tied to the exact constants it was produced under.
"""
function save_three_way(path, net::ThreeWayEvidenceNet, head_log_odds; meta = (;))
    d = dirname(path)
    isempty(d) || mkpath(d)
    arch = three_way_arch(net)
    tmp  = path * ".tmp"
    jldsave(tmp;
        schema_version = P13_NET_SCHEMA,
        model_state    = Flux.state(net),
        input_dim      = arch.input_dim,
        num_summaries  = arch.num_summaries,
        width          = arch.width,
        head_log_odds  = head_log_odds,
        tau            = isdefined(@__MODULE__, :P13_TAU) ? Float64(P13_TAU) : nothing,
        consts_sha     = p13_consts_sha(),
        cut_variant    = P13_CUT_VARIANT,
        meta           = merge((generated = string(Dates.now(Dates.UTC)) * "Z",), meta))
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "model_state") "save_three_way: integrity check failed, $tmp missing model_state"
        @assert haskey(f, "head_log_odds") "save_three_way: integrity check failed, $tmp missing head_log_odds"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end

"""
    load_three_way(path) -> NamedTuple

Reload a net persisted by `save_three_way`, REBUILDING the topology from the stored metadata and
loading the stored state into it (never deserializing an estimator object -- see
`save_three_way`). Returns
`(net, head_log_odds, input_dim, num_summaries, width, tau, consts_sha, cut_variant,
schema_version, meta)`.

Errors when the file is absent, when a required key is missing, or when the stored schema
version is not the one this file writes -- a silently-accepted stale layout would load a net
whose head semantics are not the ones documented here.
"""
function load_three_way(path)
    isfile(path) || error("load_three_way: no model file at $path")
    d = JLD2.load(path)
    for k in ("model_state", "input_dim", "num_summaries", "width", "head_log_odds")
        haskey(d, k) || error("load_three_way: $path missing $k key")
    end
    schema = get(d, "schema_version", nothing)
    schema == P13_NET_SCHEMA ||
        error("load_three_way: $path has schema_version $schema, expected $P13_NET_SCHEMA")
    net = build_three_way_net(d["input_dim"]; num_summaries = d["num_summaries"],
                              width = d["width"])
    Flux.loadmodel!(net, d["model_state"])
    return (net = net, head_log_odds = d["head_log_odds"],
            input_dim = d["input_dim"], num_summaries = d["num_summaries"],
            width = d["width"], tau = get(d, "tau", nothing),
            consts_sha = get(d, "consts_sha", nothing),
            cut_variant = get(d, "cut_variant", nothing),
            schema_version = schema, meta = get(d, "meta", nothing))
end

"""
    p13_consts_sha() -> String

The sha256 of `spike/p13/consts.jl`, hex-encoded. Recorded in every persisted artifact so a
reported number can be tied to the exact frozen pre-registration it was produced under: the
pre-registration forbids EDITING a constant, so a changed sha on a file that should only ever
have been APPENDED to is itself the audit signal.
"""
p13_consts_sha() = bytes2hex(SHA.sha256(read(joinpath(@__DIR__, "consts.jl"))))

"""
    P13_CONSTS_SHA :: NamedTuple

The TWO sha256 values `spike/p13/consts.jl` has legitimately held. **A third value is drift.**

  * `pre_amendment`  -- `c42cc8e` + `bfba6ac`; the value EVERY 13-10..13-16 artifact recorded at
    train time, including `spike/p13/three_way_net.jld2`'s `consts_sha`.
  * `post_amendment` -- after `13-D15-AMENDMENT.md` (CHANGE A: the c2/c3 operative pair, the
    dropped redundancy arm, the unmasked anchors; CHANGE B: the include-guard sentinel),
    2026-07-29.

Both literals were copied from live command output
(`bytes2hex(SHA.sha256(read("spike/p13/consts.jl")))`), never retyped and never copied out of a
document. **This amendment authorises no further amendment**, so no third value is legitimate.

WHY PINNING BOTH SIDES IS LEGITIMATE, AND WHY IT IS STRICTLY STRONGER THAN WHAT IT REPLACES.
Training reads NO `P13_REAL_*` constant (grep-verified, `13-D15-AMENDMENT.md` section 3), so the
amended values are ones the trained net never consumed -- the divergence between the artifact's
recorded sha and the current file's is by DESIGN, not drift. The assertion it replaces,
`h.consts_sha == p13_consts_sha()`, only required the two sides to AGREE WITH EACH OTHER: a
retrain plus an undisclosed edit would have satisfied it silently. Pinning both sides to named,
dated literals catches that.

REJECTED, and recorded so a later reader knows it was considered: the widened form
`h.consts_sha == p13_consts_sha() || h.consts_sha == PRE_AMENDMENT_SHA`. Once that disjunction
exists it accepts ANY future drift silently -- a second, third or undisclosed edit to
`consts.jl` also passes, because the artifact side alone satisfies the first branch. That would
convert a real integrity check into decoration, which is the one outcome the amendment must not
produce.
"""
const P13_CONSTS_SHA = (
    pre_amendment  = "100e97a37bb470bb7cb4fbdbd8e33f219d83adcf61fe398a90cf3ad3391f3fa6",
    post_amendment = "d6a6e63bf19bb13ad5c61eeaaf9fc18fab27c5aa7ba879bad5c493a5c71ae247",
)

end  # guard: !isdefined(:ThreeWayEvidenceNet)

# --- Script entry point ---------------------------------------------------------------------
# GUARDED so `include`-ing this file (the test suite, the trainer, the gate) can NEVER trigger
# training. This file deliberately has no runnable training path of its own.
if abspath(PROGRAM_FILE) == @__FILE__
    println("""
    spike/p13/net.jl defines the three-way evidence net and trains nothing.

      architecture   build_three_way_net(three_way_input_dim(G, n_cond))
      loss           masked_two_head_bce   (D-11 per-head class restriction)
      training       train_three_way       (route R1: NeuralEstimators' own train loop)
      read surface   three_way_log_bf / three_way_probs
      persistence    save_three_way / load_three_way

    The real training run is plan 13-11's trainer, which supplies the labelled three-way
    pairs, the Phase-11 conditioning length and the measured per-head corrections. The only
    training this file's gate performs is the 200-sample two-epoch A5 smoke in
    spike/test/test_p13_net.jl.

      julia --project=spike spike/test/test_p13_net.jl
    """)
end
