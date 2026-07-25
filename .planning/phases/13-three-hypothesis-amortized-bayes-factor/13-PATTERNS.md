# Phase 13: Three-Hypothesis Amortized Bayes Factor - Pattern Map

**Mapped:** 2026-07-25
**Files analyzed:** 12 (10 new, 1 modified, 1 new planning doc)
**Analogs found:** 12 / 12 (5 exact, 6 role-match, 1 partial-composite)

---

## Reading rules for the planner (read before the tables)

1. **`src/` is READ-ONLY REFERENCE in this phase (D-01 + CLAUDE.md).** Several of the strongest
   analogs below live in `src/`. Every one of them is marked **[READ-ONLY REF]**. The executor
   *copies the pattern into `spike/p13/`*; it never edits the analog. `src/results.jl:172`'s
   `log_bf_simplex` misnomer is the highest-temptation edit in the whole phase — 13-RESEARCH §G3
   Open Question 1 recommends **NOT** editing it; carry the rename as a productionization item.
2. **`spike/Project.toml` / `spike/Manifest.toml` must NOT change.** No analog below implies a new
   package. If any task drafts a `Pkg.add`, that is a planning error (13-RESEARCH §Standard Stack).
   `roc_auc`, `_bin_calibration` and the Otsu mask are all already in-repo, hand-rolled or frozen.
3. **Spike-lane analogs beat `src/` analogs** where both exist. `spike/validation/train_ratio.jl`
   is the spike-lane twin of `src/amortized/train_ratio.jl`; prefer the spike copy for style,
   the `src/` copy for the D-10 verbatim numbers (they agree — 13-RESEARCH §B2 verified all eight).
4. **RNG discipline is a hard, repeated pattern.** Two layers, always both:
   - Philox counter-based per-index: `Philox4x(UInt64, (master ⊻ SALT, UInt64(idx)))`
     (`spike/data/seeding.jl:75-96`, `spike/validation/harness.jl:93-94`).
   - `Random.seed!(derived)` immediately before any `train` / `sampleposterior` call, because
     Flux's `DataLoader(shuffle=true)` draws from the **global** RNG (13-RESEARCH Pitfall 6).
   Never `Random.seed!(i)` inside a loop; never a bare global-RNG draw in generation code.

---

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `spike/p13/consts.jl` | config (pre-registration) | static / no I/O | `spike/validation/consts.jl` (whole file) + `test/gate/gate_consts_8_v2.jl:260-345` | **exact** (two-part) |
| `spike/p13/preconditions.jl` | config guard / validator | file-I/O + assert | `spike/validation/harness.jl:44-94` + `.planning/spikes/014-bf-sim-validation/run_sim_bf.jl:40-89` | role-match |
| `spike/p13/labels.jl` | utility (pure transform) | transform | `src/amortized/train_ratio.jl:82-134` **[READ-ONLY REF]** / `spike/validation/train_ratio.jl:231-262` | **exact** |
| `spike/p13/net.jl` | model (architecture + loss + train + read surface) | batch training + request-response | `spike/validation/train_ratio.jl:264-364` + `src/amortized/train_ratio.jl:62-188` **[READ-ONLY REF]** + `src/amortized/bf.jl:61-100` **[READ-ONLY REF]** | **exact** (topology), role-match (head/loss) |
| `spike/p13/tau_probe.jl` | service (measurement probe) | batch simulate → reduce | `spike/validation/harness.jl:96-136` + `spike/validation/ood.jl:305-345` + `spike/data/seeding.jl:98-114` | role-match (composite) |
| `spike/p13/alpha_series.jl` | utility (image transform + invariants) | file-I/O / image transform | `spike/validation/ood.jl:412-502` (graded `misspec_*` families) + `:513-582` (invariance controls) | role-match |
| `spike/p13/result.jl` | model (result type) | request-response | `src/results.jl:128-159` (`AmortizedColocResult` + 4 accessors) **[READ-ONLY REF]** | **exact** |
| `spike/p13/run_three_way_gate.jl` | script (reported gate runner) | batch, exits nonzero | `spike/validation/run_bf.jl` (whole file) | **exact** |
| `spike/p13/run_alpha_series.jl` | script (reported runner, no gate) | batch, artifact-only | `spike/validation/run_ood.jl` / `run_bf.jl:59-120` (artifact + figure, no `@test`) | role-match |
| `spike/test/test_p13.jl` | test | unit fixtures | `spike/test/test_bf.jl` (whole file) + `spike/test/test_sbc.jl:36-95` | **exact** |
| `spike/test/runtests.jl` (**modified**) | test harness / resolve gate | integration | its own clauses (d)–(h) at `:59-113` + the include block at `:118-143` | **exact** (self-analog) |
| `.planning/phases/13-.../13-SC2-AMENDMENT.md` | doc (pre-registration amendment) | — | `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md:1-80` | **exact** |

---

## Pattern Assignments

### `spike/p13/consts.jl` (config, static)

**Primary analog:** `spike/validation/consts.jl` — structure, banner, guard idiom.
**Secondary analog:** `test/gate/gate_consts_8_v2.jl:260-345` — the *executable* forbidden-seed
recompute, which the primary analog lacks. 13-RESEARCH §K3 requires both.

**Copy: the anti-snooping banner + single-key re-inclusion guard** (`consts.jl:21-42, 80`):

```julia
# spike/validation/consts.jl --- Phase-5 pre-registered constants (SBC/BF/OOD).
#
# THE ANTI-SNOOPING CONTRACT (SBC-03 / D-02). Every M/L/bin-count and every
# pass/fail threshold the Phase-5 validation bundle scores against is LOCKED in
# this ONE file and committed BEFORE any reported run. ...
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local constants; touches no src/.
#
# Guarded as ONE block keyed on :SBC_M so re-inclusion under runtests.jl ... is a
# silent no-op -- a redefinition to the same value would otherwise warn on a `const`.

if !isdefined(@__MODULE__, :SBC_M)
    ...
end
```

**Change:** key the guard on a Phase-13 sentinel (e.g. `:P13_DEV_SEED`), and group sections by
decision ID (τ probe / cut variant / stratification / D-12 gate / D-13 calibration / α ladder)
rather than by SBC/BF/OOD.

**Copy: grouped constants with inline rationale + explicit disjointness comment**
(`consts.jl:44-51, 72-79`):

```julia
    const SBC_M          = 2000     # pre-registered SBC draws (θ*~π→simulate→infer)
    const SBC_L          = 999      # posterior draws per SBC draw; L+1 = 1000 …
    const SBC_ECE_GREEN  = 0.05     # ECE ≤ GREEN  → :green traffic-light verdict
    const SBC_ECE_YELLOW = 0.10     # ECE ≤ YELLOW → :yellow, else :red
    ...
    # RESERVED reported stream, consumed ONLY by the Wave-3 run (05-04). DISJOINT
    # from NPE_MASTER_SEED = 0xC0FFEE (train_npe.jl:65) so the reported SBC/BF/OOD
    # numbers never reuse the stream the net was trained/held-out on (D-02).
    const VAL_MASTER_SEED = 0x5BC0FFEE   # ≠ 0xC0FFEE (NPE_MASTER_SEED): fresh disjoint stream
    const VAL_FIX_SEED    = 0xF1F7ED
```

**Change:** re-declare `P13_ECE_GREEN = 0.05` / `P13_ECE_YELLOW = 0.10` **locally** rather than
importing `SBC_ECE_*` — 13-RESEARCH §I2 confirms every `test/gate/gate_consts_*.jl` re-declares
them; each pre-registration is self-contained.

**Copy: the recompute-don't-trust-comments forbidden-seed machinery**
(`test/gate/gate_consts_8_v2.jl:268-285, 300-330` — this is the part `spike/validation/consts.jl`
does NOT have):

```julia
    # v1 mixing constants, replicated verbatim so the v1 PROD_SEED[G] values can be RECOMPUTED
    # here and added to the forbidden set ...
    const PROD_SALT   = 0x94D0_49BB_1331_11EB
    const PROD_MASTER = 0x0000_0000_09E3_779B

    function _derive_prod_seed(G::Integer)
        rng = Philox4x(UInt64, (UInt64(PROD_MASTER) ⊻ PROD_SALT, UInt64(G)))
        s = rand(rng, UInt64)
        while s in (VAL_MASTER_SEED, NPE_MASTER_SEED, UInt64(0)); s = rand(rng, UInt64); end
        return s
    end
    const PROD_SEED = Dict{Int,UInt64}(G => _derive_prod_seed(G) for G in (4, 8, 16, 32))

    _forbidden_seeds() = (UInt64(0), NPE_MASTER_SEED, VAL_MASTER_SEED, DEFAULT_MASTER_SEED,
                          DEV_SEEDS..., values(PROD_SEED)...)
```

**Change:** name it `_p13_forbidden()`, populate it from the 13-RESEARCH §K1 inventory, and pull
`PROD_SEED` / `PROD_SEED_V2` in via the **isolated read-only module** idiom (next excerpt) instead
of re-deriving them here.

**Copy: the isolated read-only include of a frozen pre-registration**
(`.planning/spikes/014-bf-sim-validation/run_sim_bf.jl:40-68`):

```julia
# Load the FROZEN amended pre-registration into an ISOLATED module purely to READ its
# PROD_SEED / PROD_SEED_V2 / SBC_IMSIZE_* (no gate is run; the file is not modified).
module _GC
    include(joinpath(@__DIR__, "..", "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end

const _FORBIDDEN = Set{UInt64}(vcat(
    UInt64(0x0000_0000_0000_0001),                 # DEFAULT_MASTER_SEED
    UInt64(0x0000_0000_00C0_FFEE),                 # NPE_MASTER_SEED
    ...
    collect(UInt64, values(_GC.PROD_SEED))...,
    collect(UInt64, values(_GC.PROD_SEED_V2))...,
))
@assert !(UInt64(DEV_SEED) in _FORBIDDEN) "DEV_SEED collides with a forbidden seed"
```

**Change:** path is `joinpath(@__DIR__, "..", "..", "test", "gate", "gate_consts_8_v2.jl")` from
`spike/p13/` (two levels, not three). Add `P13_SALT` to a parallel salt-disjointness assert.

**Phase-13-specific addition with no analog:** the **two-commit τ discipline** (13-RESEARCH §K3) —
commit with `P13_TAU = nothing` + the probe spec frozen, run the probe, then commit
`P13_TAU = <measured>` in a separately-marked block. No existing consts file has a Tier-2
"measured after the probe" slot; the executor must invent the marker comment and the plan must
require the two commits.

---

### `spike/p13/preconditions.jl` (config guard, file-I/O + assert)

**Analog A — the guarded, order-dependent include block** (`spike/validation/harness.jl:48-58`):

```julia
# --- ORDER MATTERS: consts first, then the trainer/inference surface (load_npe +
#     posterior_for + rho_draws), then the simulator halves, contract, encode, and
#     the seeding primitives whose salt idiom val_rng mirrors. Guarded for idempotency.
isdefined(@__MODULE__, :VAL_MASTER_SEED) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :load_npe)        || include(joinpath(@__DIR__, "..", "npe", "train_npe.jl"))
isdefined(@__MODULE__, :posterior_for)   || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))
isdefined(@__MODULE__, :sample_prior)    || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair)   || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)       || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)      || include(joinpath(@__DIR__, "..", "data", "encode.jl"))
isdefined(@__MODULE__, :sample_rng)      || include(joinpath(@__DIR__, "..", "data", "seeding.jl"))
```

**Copy verbatim as the p13 include preamble** (paths become `"..", "simulator", ...` from
`spike/p13/`). This is the file that "binds the Phase-11 names once" per 13-RESEARCH §D3 —
every later p13 file includes *this*, not the eight paths.

**Analog B — the single-load frozen-model accessor** (`harness.jl:70-80`):

```julia
function load_frozen_model(path = joinpath(@__DIR__, "..", "npe", "trained_npe.jld2"))
    return load_npe(path)
end
```

**Change:** point at the **Phase-11** artifact, and make absence a loud, instructional `error`
(13-RESEARCH §D2 supplies the message text verbatim, including "Phase 13 must NOT fall back to the
shipped grid-8 basis").

**Analog C — provenance assertion against the training joint**
(`.planning/spikes/014-bf-sim-validation/run_sim_bf.jl:83-89`):

```julia
let prov = PC.training_imsize_provenance(NPE)
    @assert prov.imsize_set == IMSET "training imsize_set != F5 mixture -- covariate shift!"
    @assert prov.imsize_weights == IMWTS "training imsize_weights != F5 weights"
    println("training-joint imsize provenance OK: ", prov.imsize_set, " w=", prov.imsize_weights)
end
```

**Copy:** the assert-then-print shape. **Change:** add the 13-RESEARCH §D2 asserts —
`hasproperty(h, :zt)`, `hasproperty(h, :θzt)`, `h.lambda_range isa Tuple`,
`size(h.θzt.mean, 1) == 8`.

**Analog D — honest `:not_trained` status instead of a crash**
(`test/gate/run_gate.jl:391, 417`): `if the NPE artifact is ABSENT returns (status = :not_trained, grid)`.
Use this shape for the *runner scripts*; use the hard `error` for the precondition file itself
(D-02 makes the missing Phase-11 net a block, not a skip).

---

### `spike/p13/labels.jl` (utility, transform)

**Analog:** `src/amortized/train_ratio.jl:82-134` **[READ-ONLY REF]** (twin at
`spike/validation/train_ratio.jl:231-262`).

**Copy: the measured-never-assumed log-odds function, docstring style and error included**
(`train_ratio.jl:82-96`):

```julia
"""
    measure_log_prior_odds(model_index) -> Float64

The MEASURED log prior-odds of the coloc event (Pitfall 5 — measured, NEVER assumed 0):
`log(p / (1 − p))` with `p = mean(model_index)` the prior fraction of coloc labels. ...
Returns a finite scalar; errors on degenerate balance.
"""
function measure_log_prior_odds(model_index::AbstractVector{<:Real})
    p = mean(model_index)
    (0.0 < p < 1.0) ||
        error("measure_log_prior_odds: degenerate label balance p=$p (need 0<p<1)")
    return log(p / (1 - p))
end
```

**Change (load-bearing, 13-RESEARCH §F4):** becomes
`measure_head_log_odds(class_labels; positive, negative)` returning `log(n_pos / n_neg)` counted
**only over the head's own two classes**. Keep the degenerate-balance `error`. The docstring must
state *why it is not the binary formula* (the shuffled-θ NRE construction already contains the
ratio; a plain BCE head does not).

**Copy: the label rule lives beside its threshold constant, and the threshold is a named const**
(`train_ratio.jl:55`, `spike/validation/train_ratio.jl:97-100`):

```julia
const RATIO_SPLIT_THRESHOLD = 0.0     # D-07 coloc/null split on the ρ_true CONTRAST (≡ {Δμ>0})
...
        model_index[k]    = Float32((ρ[i1] - ρ[i2]) > threshold)
```

**Change (this is the entire point of D-05):** the one-factor `Δρ > 0` rule is **replaced** by the
two-factor cut on `ρ_sample` level × control contrast (13-RESEARCH §E1 variant (a)). `P13_TAU` is
**read from `consts.jl`, never inlined**. Add the executable counter-example assertion
`three_way_label(0.8, 0.9) != EXCLUSION` (Pitfall 3).

**Copy: the cheap prior-only label draw** (`spike/validation/train_ratio.jl:246-262`) — labels that
depend only on θ can be drawn without simulating images, which is how `log_prior_odds` is measured
at large n:

```julia
function prior_label_draws(n::Integer; rng = val_rng(), threshold = RATIO_SPLIT_THRESHOLD)
    labels = Vector{Float32}(undef, n)
    for j in 1:n
        ρs = sample_prior(rng).ρ_true
        ρc = sample_prior(rng).ρ_true
        labels[j] = Float32((ρs - ρc) > threshold)
    end
    return labels
end
```

**Change:** returns a `Vector{ThreeWayClass}` (or the 4-row target matrix from `head_targets`), and
the two per-head log-odds are measured from it. This function is also what makes the §E1 class-mass
table reproducible in-repo.

---

### `spike/p13/net.jl` (model: architecture + loss + train + read surface)

**Analog A — the trunk, to be lifted VERBATIM (D-10)**
(`src/amortized/train_ratio.jl:52-80` **[READ-ONLY REF]**; identical in
`spike/validation/train_ratio.jl:83-87, 304-311`):

```julia
const RATIO_NUM_SUMMARIES   = 64      # learned-summary bottleneck feeding the model-index head
const RATIO_SUMMARY_WIDTH   = 256     # summary-conditioner hidden width (mirrors the NPE lever)
...
function build_ratio_estimator(input_dim::Integer;
                               num_summaries::Integer = RATIO_NUM_SUMMARIES,
                               width::Integer = RATIO_SUMMARY_WIDTH)
    input_dim >= 1 || throw(ArgumentError("build_ratio_estimator: input_dim must be ≥ 1"))
    W = width
    net = Chain(Dense(input_dim, W, gelu), Dense(W, W, gelu),
                Dense(W, W, gelu), Dense(W, num_summaries))
    return RatioEstimator(net, 1; num_summaries = num_summaries)
end
```

**Copy:** the four `Dense` lines, the two constants (renamed `P13_*`), the `input_dim >= 1`
`ArgumentError` guard, and the "shared by trainer and loader so both build the identical topology"
contract in the docstring.
**Change:** drop the `RatioEstimator(...)` wrap (D-09); return
`ThreeWayEvidenceNet(trunk, Dense(num_summaries, 2))` with `Flux.@layer`. **Never write the input
width as a literal** — derive `ratio_input_dim(G) + n_cond` (13-RESEARCH §D2/Pitfall 4).

**Analog B — the training recipe, copied verbatim so D-10's attribution argument holds**
(`src/amortized/train_ratio.jl:155-183`; spike twin at `spike/validation/train_ratio.jl:278-320`):

```julia
    # Train/val split (contiguous — the stream is already i.i.d. across columns).
    nval = clamp(round(Int, val_frac * n), 1, n - 1)
    ntr  = n - nval
    Z_tr = Z_pair[:, 1:ntr]
    Z_va = Z_pair[:, ntr+1:end]
    midx_tr = reshape(model_index[1:ntr], 1, :)         # 1×ntr (num_parameters = 1)
    midx_va = reshape(model_index[ntr+1:end], 1, :)
    ...
    # FIXED-DATA train form. D-06: `use_gpu` plumbed. NO custom loss (hard-coded logit-BCE).
    # AdamW args Float64 to match NeuralEstimators' Float64 CosAnneal lr_schedule.
    est = train(est, midx_tr, midx_va, Z_tr, Z_va;
                epochs = epochs, batchsize = batchsize, use_gpu = use_gpu,
                optimiser = Flux.Optimisers.AdamW(learning_rate, (0.9, 0.999), weight_decay),
                stopping_epochs = stopping_epochs, verbose = verbose)
```

Defaults to copy exactly (`train_ratio.jl:155-163`): `n = 48_000, epochs = 300, batchsize = 128,
learning_rate = 2.5e-4, weight_decay = 1e-4, val_frac = 0.15, stopping_epochs = 40`.
**Copy also the AdamW-Float64 comment** — it is a real gotcha, not decoration.
**Change:** `midx_*` becomes the **4-row** target matrix `[y_C; w_C; y_E; w_E]` (13-RESEARCH §B3);
`loss = masked_two_head_bce` **is** passed (legal here — the generic `_loss(estimator, loss) = loss`
fallback at `NeuralEstimators/src/train.jl:654` applies to a non-`RatioEstimator` subtype).

**Analog C — the CPU-only gate** (`spike/validation/train_ratio.jl:286`):

```julia
    use_gpu && throw(ArgumentError("train_ratio: use_gpu=true out of scope (CPU-only gate, D-10)"))
```

**Copy verbatim in spirit** — the spike lane keeps the hard throw (the `src/` version relaxed it to
`has_cuda_device()` only on promotion). `spike/test/runtests.jl` asserts CUDA absent.

**Analog D — the read surface** (`src/amortized/bf.jl:61-79` **[READ-ONLY REF]**):

```julia
function amortized_log_bf(est_bf, Z_pair, log_prior_odds; use_gpu::Bool = false)
    Zc   = Z_pair isa AbstractVector ? reshape(Z_pair, :, 1) : Z_pair
    grid = reshape(Float32[0.0 1.0], 1, 2)                       # null (m=0), coloc (m=1)
    lr   = logratio(est_bf, Zc; grid = grid, use_gpu = use_gpu)  # 1×2: [logr(m=0) logr(m=1)]
    return (lr[1, 2] - lr[1, 1]) - log_prior_odds
end
```

**Copy:** the vector-or-matrix reshape guard, the `use_gpu = false` CPU default, the one-line
"ONE forward pass" contract.
**Change:** `logratio(...)` is replaced by a direct `net(Zc)` call returning 2×1 logits; the return
is the NamedTuple `(coloc = s[1,1] - head_log_odds.coloc, exclusion = s[2,1] - head_log_odds.exclusion,
random = 0.0)` (13-RESEARCH §F4). **Do NOT copy the `- log_prior_odds` term by analogy** —
Pitfall 1; the per-head correction is a different quantity with a different denominator.

**Analog E — `pair_encode`, reused UNCHANGED** (`src/amortized/bf.jl:91-100` **[READ-ONLY REF]**;
spike twin `spike/validation/train_ratio.jl:185-188`):

```julia
function pair_encode(Zs::AbstractVector, Zc::AbstractVector)
    length(Zs) == length(Zc) ||
        throw(DimensionMismatch("pair_encode: Zs ($(length(Zs))) and Zc ($(length(Zc))) must be the same length"))
    iseven(length(Zs)) ||
        throw(ArgumentError("pair_encode: summary length $(length(Zs)) must be 2·G² (even)"))
    nc = length(Zs) ÷ 2                                   # continuous rows = G² (cont_rows(G))
    d  = @view(Zs[1:nc]) .- @view(Zc[1:nc])
    return Float32.(vcat(Zs, Zc, d))
end
```

**Copy the call, not the code** — reach it through the read-only include chain. The λ row is
appended **outside**: `Z_pair = vcat(pair_encode(Zs_128, Zc_128), λ)`. The `iseven` guard is the
tripwire that catches the wrong form (Pitfall 4); do not weaken it.

**Analog F — atomic persistence** (`spike/validation/train_ratio.jl:331-347`):

```julia
function save_ratio(path, est, zt, log_prior_odds; meta = (;))
    mkpath(dirname(path)); tmp = path * ".tmp"
    jldsave(tmp; schema_version = RATIO_MODEL_SCHEMA, estimator = est, zt = zt,
        log_prior_odds = log_prior_odds, num_summaries = RATIO_NUM_SUMMARIES,
        split_threshold = RATIO_SPLIT_THRESHOLD,
        meta = merge((generated = string(Dates.now(Dates.UTC)) * "Z",), meta))
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "estimator") "save_ratio: integrity check failed, $tmp missing estimator"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end
```

**Copy the whole `.tmp` → reopen-integrity-check → `mv(force=true)` idiom** and the paired
`load_ratio` reader (`:356-364`), including the `schema_version` bump discipline.
**Change:** persist `head_log_odds` (a 2-field NamedTuple) instead of the scalar
`log_prior_odds`, plus `P13_TAU` and the consts-file sha for provenance.

**Analog G — the guarded script entry point** (`spike/validation/train_ratio.jl:369-378`):

```julia
# GUARDED so `include`-ing this file (bf.jl / test_bf.jl) never triggers training.
if abspath(PROGRAM_FILE) == @__FILE__
    m   = load_frozen_model()
    est = train_ratio(m; n = RATIO_TRAIN_N, epochs = RATIO_EPOCHS, verbose = true)
    ...
end
```

**Copy verbatim.** This is what lets `test_p13.jl` include `net.jl` without training.

---

### `spike/p13/tau_probe.jl` (service, batch simulate → reduce)

No single analog; compose three.

**Analog A — the draw→simulate→summary chain** (`spike/validation/harness.jl:106-113`):

```julia
function draw_simulate_infer(m, rng; imsize = SBC_IMSIZE, N = SBC_L)
    θ   = sample_prior(rng)                                          # θ* ~ π(θ)
    mci = build_mci(simulate_pair(rng, θ; imsize = imsize))          # forward sim → MCI
    Z   = standardize_summary(encode_d01(patch_summary(mci)), m.zt, :min)  # FROZEN m.zt
    ...
end
```

**Copy:** the exact four-call chain `sample_prior → simulate_pair → build_mci → patch_summary`.
**Change (13-RESEARCH §E4, the sequencing win):** the τ probe **stops before** `standardize_summary`
— it reads `encode_d01(patch_summary(mci))` raw, so it needs **no `zt` and no trained net** and can
run as soon as Phase 11's *simulator* merges. `ρ_true` is **pinned**, not drawn:
`merge(sample_prior(rng), (; ρ_true = ρ))`. `imsize` comes from `sample_imsize(rng)`
(the F5 mixture), **not** the literal `SBC_IMSIZE`.

**Analog B — the summary row partition the `m̄` statistic depends on**
(`spike/data/encode.jl:56-64`):

```julia
function encode_d01(M::AbstractMatrix)
    vals = vec(coalesce.(M, 0.0))                 # 64-dim, missing->0 (column-major)
    mask = vec(Float64.(.!ismissing.(M)))         # 64-dim, 1=present / 0=missing
    return vcat(vals, mask)                        # 128-dim; rows 1:64 vals, 65:128 mask
end
```

**Consequence to copy into the probe:** `m̄` must average rows `1:G²` **weighted by the mask rows
`G²+1:2G²`**, or the missing→0 imputation silently biases the statistic toward 0 (13-RESEARCH §E3
"one degenerate all-missing grid can move `m̄` by O(1)").

**Analog C — the hand-rolled tie-aware AUC** (`spike/validation/ood.jl:309-335`):

```julia
"""
    roc_auc(neg_scores, pos_scores) -> (fpr, tpr, auc)

Hand-rolled ROC curve + AUC ... the AUC is the tie-aware Mann–Whitney U statistic (exact,
robust to threshold degeneracies) `= P(score_pos > score_neg) + ½·P(=)`. ... No dependency (~15 lines).
"""
function roc_auc(neg_scores::AbstractVector, pos_scores::AbstractVector)
    P = length(pos_scores); N = length(neg_scores)
    (P == 0 || N == 0) && return (Float64[0.0, 1.0], Float64[0.0, 1.0], NaN)
    ...
    u = 0.0
    for p in pos_scores, n in neg_scores
        u += p > n ? 1.0 : (p == n ? 0.5 : 0.0)
    end
    return (fpr, tpr, u / (P * N))
end
```

**Reuse UNCHANGED** (include `spike/validation/ood.jl` or lift the function; do **not** add
`ROCAnalysis` — `runtests.jl:108-109` asserts it absent). The same function serves the D-12
per-class AUC gate.

**Analog D — disjoint per-index streams** (`spike/data/seeding.jl:75-96`):

```julia
sample_rng(master_seed::Integer, idx::Integer) =
    Philox4x(UInt64, (UInt64(master_seed), UInt64(idx)))

holdout_rng(master_seed::Integer, idx::Integer) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ HOLDOUT_SALT, UInt64(idx)))
```

**Copy the salt-XOR-into-the-first-key-word idiom** to carve the probe's two **unpaired** arms
(13-RESEARCH §E3 — unpaired is the single most important difference from the Phase-11 probe, and
the consts file must say so).

---

### `spike/p13/alpha_series.jl` (utility, image transform)

**Analog A — graded perturbation families with a `level`/α knob**
(`spike/validation/ood.jl:477-502`):

```julia
"""
    misspec_background(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS)
        -> Vector{Matrix{Float64}}

Family 4 (D-03) — BACKGROUND/illumination mismatch. Applies a multiplicative illumination
gradient + radial vignette and an autofluorescence bleed BEYOND the simulator's
`Uniform(0,0.1)` offset; `level` scales the gradient, vignette and bleed.
"""
function misspec_background(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS)
    _guard_misspec_imsize(imsize)
    base = simulate_pair(rng, θ; imsize = imsize)
    ...
    return [Matrix{Float64}(max.(ch .* max.(grad, 0.0) .* max.(vign, 0.0) .+ bleed, 0.0))
            for ch in base]
end

# The four positive-control families, in fixed order (the reported ROC grid, D-03).
const OOD_FAMILIES = (texture = misspec_texture, noise = misspec_noise, ...)
```

**Copy:** the signature shape `(rng, input; knob)`, the `Vector{Matrix{Float64}}` return, the
`max.(·, 0.0)` clamp on every output channel, and the named-tuple registry of families in fixed
order (the α ladder should be an equivalent frozen tuple in `consts.jl`).
**Change:** α is continuous in `[0,1]`, not an integer level, and the transform is applied to an
**existing** pair rather than regenerating from `simulate_pair`, so `α = 0` can be bitwise-exact.

**Analog B — the zero-set-preserving transform, which is the exact trap D-16 must dodge**
(`spike/validation/ood.jl:513-528`):

```julia
"""
    negctrl_affine(pair; a = 2.0, b = 0.5) -> Vector{Matrix{Float64}}

... Preserving the zero set is load-bearing: a blanket `+b` would turn the `max(·,0)` read-noise
zeros into signal and change which pixels the src `_exclude_zero` drops per patch, perturbing a
few patch correlations (empirically ~3/64). ...
"""
function negctrl_affine(pair; a::Real = 2.0, b::Real = 0.5)
    a > 0 || throw(ArgumentError("negctrl_affine: slope a must be > 0, got $a"))
    return [Matrix{Float64}(map(v -> v > 0 ? a * v + b : v, ch)) for ch in pair]
end
```

**This docstring is the single most valuable thing to copy in the whole α-series task.** It is the
in-repo precedent proving that changing which pixels are zero changes the summary via
`_exclude_zero`. The 13-RESEARCH §J3 rule (strictly-positive background floor `b`, assert
`all(y' .> 0)`) is the same insight applied in the opposite direction.

**Supporting reference — `_exclude_zero` itself** (`src/colocalization.jl:187-203`
**[READ-ONLY REF]**, do not modify, do not re-implement):

```julia
function _exclude_zero(a::Vector{Union{T,Missing}}, b::Vector{Union{T,Missing}}) where T <: Number
    ...
        (ismissing(ai) || ismissing(bi)) && continue
        (iszero(ai) || isnan(ai) || iszero(bi) || isnan(bi)) && continue
```

**Analog C — the invariance verifier** (`spike/validation/ood.jl:567-581`):

```julia
"""
    verify_summary_invariance(pair, transform) -> Float64

Empirically verify (D-04) that `transform` leaves the 8×8 patch-Pearson summary statistically
UNCHANGED ... Returns `NaN` if either summary is fully missing.
"""
function verify_summary_invariance(pair, transform)
    v0 = Float64.(collect(skipmissing(patch_summary(build_mci(pair)))))
    v1 = Float64.(collect(skipmissing(patch_summary(build_mci(transform(pair))))))
    (isempty(v0) || isempty(v1)) && return NaN
    return ApproximateTwoSampleKSTest(v0, v1).δ
end
```

**Copy the shape** for the α-series invariants. **Change:** the α=0 check is **stronger than
statistical** — assert bitwise `alpha_segregate(x, y, M, 0.0; b) == y` (exact `==`, not `≈`), plus
`sum(y') ≈ sum(y)` at `rtol = 1e-10` and mask-α-invariance.

**Analog D — the mask rule, frozen and already shipped** (`src/LoadImages.jl:235-241`
**[READ-ONLY REF]**):

```julia
function _calculate_mask(img::AbstractMultiChannelImage)
    return [
        ch .> thr
        for (ch, thr) in
            zip(image_data(img), otsu_thresholds(img))
    ]
end
```

**Reuse via the read-only include** (`spike/contract.jl:46-47` already includes `LoadImages.jl`);
take `M = _calculate_mask(mci)[1]`, computed **once** from the unmodified ch1 and held fixed across
the whole α ladder. This makes D-16's segmentation rule a pre-registered *existing* choice at zero
cost. The MCI construction idiom to mirror is `spike/contract.jl:60-70` (`build_mci`), which already
calls `Images.otsu_threshold.(data)`.

**Analog E — the background floor constant** (`spike/simulator/forward.jl:37, 66`):

```julia
#   (5) autofluorescence offset + BG_FLOOR  (background stays small-positive, never hard 0.0)
const BG_FLOOR   = 0.02
```

Cite this in the α-series docstring as the in-repo precedent for "background stays small-positive,
never hard 0.0"; the D-16 floor `b = quantile(vec(y), 0.05)` is its data-driven analogue and must be
asserted `> 0`.

---

### `spike/p13/result.jl` (model, request-response)

**Analog:** `src/results.jl:128-159` **[READ-ONLY REF]**.

**Copy: the struct + four-accessor block, exactly this shape** (`results.jl:146-159`):

```julia
struct AmortizedColocResult <: AbstractColocResult
    grid            :: Int
    posterior       :: Matrix{Float64}
    delta_rho_draws :: Vector{Float64}
    log_bayes_factor:: Float64
    ood             :: OODVerdict
    calibration     :: CalibrationMeta
    meta            :: NamedTuple
end

delta_rho(r::AmortizedColocResult)       = r.delta_rho_draws
bayes_factor(r::AmortizedColocResult)    = r.log_bayes_factor
is_ood(r::AmortizedColocResult)          = r.ood.flag
posterior_draws(r::AmortizedColocResult) = r.posterior
```

**Copy also the interface-error contract it must satisfy** (`results.jl:47-49`):

```julia
_iface_error(r::AbstractColocResult, f::Symbol) = error(
    "`$(f)` is not implemented for $(typeof(r)). Every `AbstractColocResult` subtype must " *
    "implement the accessor interface: delta_rho, bayes_factor, is_ood, posterior_draws.")
```

**Changes (all decision-driven):**
- The struct is defined in **`spike/p13/result.jl`**, subtyping the real `AbstractColocResult`
  reached by the read-only `include(joinpath(@__DIR__, "..", "..", "src", "results.jl"))` —
  the same read-only reach `spike/contract.jl:46-47` uses. **Zero `src/` bytes change.**
- `log_bf_simplex :: NTuple{3,Float64}` (the sketch at `src/results.jl:169-178`) becomes
  `log_bf_vs_random :: NamedTuple{(:coloc, :random, :exclusion), ...}` (D-08 + 13-RESEARCH §G1).
  The sketch's positional order is `{null, coloc, anti-coloc}` — reference class **first** — which
  is exactly the hazard the NamedTuple removes.
- `bayes_factor(r)` returns `r.log_bf_vs_random.coloc` and **must document that choice**;
  `delta_rho(r)` should fall through to `_iface_error` unless control draws are actually carried
  (13-RESEARCH §G1: do not fabricate a Δρ).
- `calibration :: CalibrationMeta` is populated from the D-13 run — note the field-order contract at
  `src/results.jl:117-126` (six `CalibrationResult` fields **plus** `grid::Int` and
  `gate::NamedTuple`), whereas the spike's `CalibrationResult` (`spike/validation/sbc.jl:65-72`) has
  only the six. The spike-lane result may carry either; state which.

**Do NOT touch `src/results.jl:161-192`.** It is comment-only, and 13-RESEARCH Open Question 1
recommends leaving it — surface as a `checkpoint:human-verify` if the user wants the rename now.

---

### `spike/p13/run_three_way_gate.jl` (script, reported gate)

**Analog:** `spike/validation/run_bf.jl` (whole file) — **exact**. The ordering below is the pattern
and it is load-bearing: *artifact persists BEFORE the gate can throw.*

**Copy: the banner that prints the locked thresholds and the stream** (`run_bf.jl:71-79`):

```julia
function main()
    println("="^78)
    println("Phase-5 REPORTED amortized Bayes factor (BF-01/02) — LOCKED consts, RESERVED stream")
    println("  sweep: N=$BF_SWEEP_N  Δρ∈[$BF_SWEEP_LO, $BF_SWEEP_HI]  imsize=$SBC_IMSIZE")
    println("  gate:  corr ≥ BF_CORR_MIN=$BF_CORR_MIN  AND  max|ΔlogBF| ≤ BF_LOGBF_TOL=$BF_LOGBF_TOL")
    println("  stream=VAL_MASTER_SEED=$(repr(VAL_MASTER_SEED))")
    println("="^78)
    @assert VAL_MASTER_SEED != NPE_MASTER_SEED "reported stream must be disjoint from training stream (D-02)"
    @assert VAL_MASTER_SEED != VAL_FIX_SEED    "reported stream must be disjoint from fixture stream (D-02)"
```

**Copy: atomic artifact write + integrity check** (`run_bf.jl:59-69`):

```julia
function _save_bf_report(path; kwargs...)
    mkpath(dirname(path)); tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "logbf_amortized") "_save_bf_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)
    return path
end
```

**Copy: persist-before-gate, with the pre-registered thresholds written INTO the artifact**
(`run_bf.jl:105-116`):

```julia
    # --- Persist the reported artifact BEFORE the gate (survives an honest fail) --------
    _save_bf_report(BF_REPORT_PATH;
        ...
        BF_CORR_MIN = BF_CORR_MIN, BF_LOGBF_TOL = BF_LOGBF_TOL,
        BF_SWEEP_LO = BF_SWEEP_LO, BF_SWEEP_HI = BF_SWEEP_HI, BF_SWEEP_N = BF_SWEEP_N,
        VAL_MASTER_SEED = UInt64(VAL_MASTER_SEED),
        generated = string(Dates.now(Dates.UTC)) * "Z")
```

**Copy: headline printed before the gate throws, then the real `@testset` gate**
(`run_bf.jl:122-140`):

```julia
    corr_pass = rep.corr >= BF_CORR_MIN
    err_pass  = rep.max_abs_err <= BF_LOGBF_TOL
    println("Overall reported BF gate: ", (corr_pass && err_pass) ? "PASS" : "FAIL")
    ...
    # --- THE REAL GATE (D-08): exits nonzero if either pre-registered criterion fails.
    # Thresholds referenced from consts.jl (never tuned). An honest failure → Phase-6.
    @testset "Reported BF pre-registered gate (BF-01/02, D-08)" begin
        @test rep.corr        >= BF_CORR_MIN      # D-08a ordering
        @test rep.max_abs_err <= BF_LOGBF_TOL     # D-08b bounded magnitude
    end
    return nothing
end

main()
```

**Changes:** gate on per-class AUC floors + confusion-matrix requirements + per-head ECE
(D-12/D-13), not on `corr`/`max|Δ|`. The binary-NRE continuity block is a **separate printed
section with no `@test`** (reported, not gated). Per Open Question 5, report the argmax confusion
matrix explicitly labelled *descriptive, not a decision rule*. Insert
`@assert success(`git diff --quiet HEAD -- src/`)`-style step-0 check (the `spike/demo.jl` precedent
named in 13-RESEARCH's test map) to prove `src/` is untouched at run time.

---

### `spike/p13/run_alpha_series.jl` (script, reported, non-gating)

**Analog:** the same `run_bf.jl` skeleton **minus** the `@testset`; closest live example of a
"artifact + figure, no assertion" runner is `run_bf.jl:105-120`:

```julia
    println("persisted reported artifact -> $BF_REPORT_PATH")
    # --- Figure (script artifact — never a gate assertion) ------------------------------
    figpath = plot_bf_agreement(rep.logbf_amortized, rep.logbf_kde; filename = "bf_agreement.png")
    println("figure -> $figpath")
```

**Copy:** the "figures are script artifacts, never gate assertions" rule (stated identically in
`spike/test/test_bf.jl:37` and `spike/validation/sbc.jl` banners), and the `plot_*` helpers living
in `spike/validation/figures.jl` (CairoMakie already resolved).
**Change:** emit the `log BF(E:R)(α)` / `log BF(C:R)(α)` / `m̄(α)` curves and the crossing point
`α*`; **no gate**, D-15/D-16 are supporting evidence. State the physical-anchor deferral to Phase 16
in the printed header (13-RESEARCH §J2 — do **not** call `open_sealed_holdout`).

---

### `spike/test/test_p13.jl` (test, unit fixtures)

**Analog:** `spike/test/test_bf.jl` (whole file) — **exact**, plus `test_sbc.jl:36-95` for the
calibration testset.

**Copy: the file skeleton, stated explicitly in the analog's own banner** (`test_bf.jl:36-37`):

```
# Mirrors test_sbc.jl: license header → using Test → guarded include of the unit under test
# → fixtures computed ONCE → one outer @testset. No figure call (figures are script artifacts).
```

**Copy: guarded include + fixtures computed once at module level** (`test_bf.jl:39-65`):

```julia
using Test
using Random

# Unit under test: the amortized BF surface (pulls bf.jl → train_ratio.jl → harness.jl).
isdefined(@__MODULE__, :amortized_log_bf) || include(joinpath(@__DIR__, "..", "validation", "bf.jl"))

# --- Frozen model + ratio net loaded ONCE (persisted; tiny-retrain only if absent) ------
const BF_FIX_M = load_frozen_model()

function _load_or_tiny_ratio()
    if isfile(RATIO_MODEL_PATH)
        return load_ratio(RATIO_MODEL_PATH)
    end
    est = train_ratio(BF_FIX_M; n = 120, epochs = 6, imsize = (64, 64),
                      rng = val_rng(VAL_FIX_SEED))
    ...
end
const BF_FIX_RATIO = _load_or_tiny_ratio()
```

**Copy this exactly** — the "load the persisted artifact, tiny-retrain only if absent, at
`n = 120, epochs = 6, imsize = (64,64)`" idiom is what keeps the gate fast. The Phase-13 A5 smoke
test (2 epochs on a 200-sample toy) slots straight into this shape.

**Copy: the fixture-stream rule** (`test_bf.jl:26-27, 99-100`):

```julia
# All fixtures draw from the FIXTURE stream val_rng(VAL_FIX_SEED) so the quick gate NEVER
# consumes the reserved reported VAL_MASTER_SEED stream (D-02).
...
        # fixtures never consume the reserved reported stream (D-02).
        @test VAL_FIX_SEED != VAL_MASTER_SEED
```

**Change:** `P13_FIX_SEED` (a second fresh seed) vs `P13_DEV_SEED`; assert both against
`_p13_forbidden()`.

**Copy: the assert-the-pre-registration-is-locked testset** (`test_bf.jl:91-101`):

```julia
    @testset "D-08 pre-registration (both criteria locked; anti-snooping)" begin
        # The pass/fail thresholds the REPORTED reproduction (05-04) scores against are the
        # committed pre-registered values — asserted-as-locked here, NEVER tuned to pass.
        @test BF_CORR_MIN  == 0.95        # D-08a correlation floor
        @test BF_LOGBF_TOL == 0.5         # D-08b bounded |Δ log-BF|
```

**Change:** assert the Phase-13 Tier-1 constants (τ grid, `P13_TAU_AUC`, cut variant,
stratification design, ECE band, AUC floors, iteration allowance = 1).

**Copy: the source-grep gate** (`test_bf.jl:67-68, 83-89`):

```julia
const BF_SRC = read(joinpath(@__DIR__, "..", "validation", "bf.jl"), String)
...
        @test !occursin("kde(", BF_SRC)
        @test !occursin("quadgk(", BF_SRC)
        @test occursin("logratio(", BF_SRC)
```

**Change/keep:** grep `spike/p13/net.jl` for the absence of `kde(`/`quadgk(` (named limit #3 — those
packages are not even in the spike env) and for the presence of the single-pass `m.head(m.trunk(`.

**Copy: the both-directions calibration test** (`test_sbc.jl:77-88`):

```julia
        # (d) _bin_calibration on PERFECTLY-calibrated synthetic input → ECE ≈ 0 (green);
        #     on MAXIMALLY-miscalibrated input → ECE near its max (red). ...
        cal_good = _bin_calibration(fill(0.5, 100), vcat(trues(50), falses(50)); n_bins = 10)
        @test cal_good.ece < SBC_ECE_GREEN
        @test sbc_traffic_light(cal_good.ece) === :green
        cal_bad = _bin_calibration(fill(1.0, 100), falses(100); n_bins = 10)
        @test cal_bad.ece > SBC_ECE_YELLOW
        @test sbc_traffic_light(cal_bad.ece) === :red
```

**Copy the exercise-both-directions principle verbatim** (it is what stops a stuck always-pass
statistic). **Add** the Phase-13-specific empty-bin MCE assertion documenting why ECE, not MCE, is
the gate (13-RESEARCH §I1 Pitfall 7).

**Copy: the CPU-only postscript** (`test_bf.jl:121-123`):

```julia
    @testset "BF ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end
```

---

### `spike/test/runtests.jl` (**modified**) — resolve-risk clause (i) + include

**Analog: the file's own clauses (d)–(h)** — clause (h) is the closest template because it is the
"this phase adds NO new package" clause (`runtests.jl:100-113`):

```julia
        # (h) RESOLVE-RISK GATE (Phase 5): the validation bundle (SBC/BF/OOD) adds NO
        #     new package — every runtime dep (NeuralEstimators/Flux/HypothesisTests/
        #     JLD2/Random123/StatsBase/CairoMakie/Distributions) is already present from
        #     Phases 1–4. ... NO ROC/normalizing-flow package may enter the env (CLAUDE.md
        #     "What NOT to Use"; RESEARCH Pitfall 7 co-resolve landmine). Turing must stay
        #     absent (it caps NeuralEstimators below the pin). Then the v0.2.1 pin is
        #     re-asserted with the Phase-5 code loadable.
        @test !haskey(Pkg.project().dependencies, "ROCAnalysis")
        @test !haskey(Pkg.project().dependencies, "MLJ")
        @test !haskey(Pkg.project().dependencies, "NormalizingFlows")
        @test !haskey(Pkg.project().dependencies, "InvertibleNetworks")
        @test !haskey(Pkg.project().dependencies, "Turing")
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
```

**Copy this block wholesale as clause (i)**, changing only the comment (Phase 13: three-way
evidence net adds NO new package; `ImageSegmentation`/`ImageMorphology` must NOT become direct deps
— `Images` 0.26.2 already provides `otsu_threshold` transitively) and re-asserting the same
NeuralEstimators UUID pin. **Also add** `@test !haskey(Pkg.project().dependencies, "CUDA")` is
already covered by clause (a).

**Copy: the bottom-of-file include convention with a one-comment rationale**
(`runtests.jl:133-138`):

```julia
# Phase-5 scaffold: the validation bundle (SBC filled by this plan 05-01; BF/OOD
# skipped placeholders until Plans 05-02/05-03) runs in the same harness so a single
# `julia --project=spike spike/test/runtests.jl` stays the gate.
include(joinpath(@__DIR__, "test_sbc.jl"))
include(joinpath(@__DIR__, "test_bf.jl"))
include(joinpath(@__DIR__, "test_ood.jl"))
```

**Change:** append `include(joinpath(@__DIR__, "test_p13.jl"))` with the equivalent Phase-13
comment. **Do not** restructure the file; append only.

---

### `.planning/phases/13-.../13-SC2-AMENDMENT.md` (doc)

**Analog:** `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md:1-80` —
**exact**, and D-12 requires it be written **before any Phase-13 result exists**.

**Copy: the header block** (`07-GATE-AMENDMENT.md:1-7`):

```markdown
# 07 — PRE-REGISTRATION AMENDMENT to the grid-8 ship-gate

**Date:** 2026-07-21
**Status:** **FROZEN SPECIFICATION. Nothing has been run against it.**
**Implements:** `test/gate/gate_consts_8_v2.jl` (new, frozen)
**Supplements — does NOT replace:** `test/gate/gate_consts_8.jl` (and `_4` / `_16`)
**Depends on:** `07-CALIBRATION-FINDINGS.md` (F1–F6), Phase-2 decisions **D-08**, **D-09**
```

**Copy: §0 provenance disclosure + the legitimacy test** (`:11-27`):

```markdown
## 0. Provenance disclosure — read this first

**This amendment was written AFTER the original grid-8 ship-gate results were seen.** ...
That is exactly why every change below is argued **only** from properties of the test design ...

The test of legitimacy applied throughout is:

> *Would this change have been made, with this justification, by someone who had seen the gate's
> code but none of its outputs?*
```

**Change — and this is the phase's strongest honesty asset:** Phase 13 inverts the disclosure.
State plainly that **no Phase-13 result exists yet**, with date + commit sha, and that the
justification is a *prior published* finding (named limit #3), not a Phase-13 outcome
(13-RESEARCH §H2). Do not blur the two positions.

**Copy: the "original is not erased" clause** (`:28-34`):

```markdown
**The original pre-registration is not erased.** ... remain **byte-unchanged** and their recorded
FAIL verdicts ... remain citable as-is. The amended gate is a **second, separately named**
pre-registration ... scored on a **fresh, disjoint seed**. Any report that cites the amended result
must cite the original FAIL alongside it. Reporting only the amended run would be exactly the thing
this document exists to prevent.
```

**Copy: the numbered-defect structure with per-item outcome-independence argument**
(`:38-68` — defect → "### Outcome-independent justification" → "*Would a blind reviewer have made
this change?*" → "### The new rule").
**Change:** two defects (13-RESEARCH §H3): (1) the KDE reference's validity ceiling
`|logBF| ≈ log(L) ≈ 6.9`; (2) the D-02 input-surface mismatch (Phase-11 129-row basis is not shared
with `compute_BayesFactor`'s ADVI path or with the shipped binary NRE) — the second is architectural
and needs no results at all.

**Also required:** annotate `.planning/ROADMAP.md:261` with a pointer to this amendment, or
`/gsd:verify-work` will score Phase 13 against superseded SC2 text (13-RESEARCH §H1).

---

## Shared Patterns

### Pattern S1 — License header + file-purpose banner (ALL new `.jl` files)
**Source:** `spike/validation/consts.jl:1-40`, `spike/data/seeding.jl:1-36`, every spike file.
**Apply to:** all 7 `spike/p13/*.jl` + `spike/test/test_p13.jl`.

```julia
#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
E-Mail: manuel.seefelder@uni-ulm.de
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany
... AGPL v3 boilerplate ...
=#

# spike/<path>.jl --- <one-line purpose> (<REQ/decision IDs>).
#
# <2-6 paragraphs: what it does, WHY the design is what it is, which decision drives it,
#  the named traps it dodges>
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; touches no src/.
```

The `DECOUPLING` line appears in essentially every spike file (`consts.jl:36`, `seeding.jl:36`,
`train_ratio.jl:57`, `encode.jl:~42`). It is not decoration — it is the phase's audit trail.

### Pattern S2 — Guarded, order-dependent includes (idempotent under `runtests.jl`)
**Source:** `spike/validation/harness.jl:48-58`, `spike/validation/train_ratio.jl:72-73`,
`spike/test/test_bf.jl:43`.
**Apply to:** every `spike/p13/*.jl` and `test_p13.jl`.

```julia
isdefined(@__MODULE__, :VAL_MASTER_SEED)   || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :load_frozen_model) || include(joinpath(@__DIR__, "harness.jl"))
```

Guard on a **symbol the included file defines**, never on a path. Paths are always
`joinpath(@__DIR__, ...)` — never interpolated, never relative (also the V5/path-traversal control
in 13-RESEARCH §Security).

### Pattern S3 — Read-only reach into `src/`
**Source:** `spike/contract.jl:30-47`.

```julia
# DECOUPLING (hard constraint, CLAUDE.md): the two include()s below are READ-ONLY.
# src/ and the root manifests stay byte-identical to baseline f581d95; this file
# and all new code live entirely under spike/.
#
# ORDER MATTERS (NOTES §4, runtests.jl evidence): StatsBase + Statistics MUST be
# in scope BEFORE including src/colocalization.jl. ... Images must be in scope before
# src/LoadImages.jl (its convenience constructor calls Images.otsu_threshold).

using StatsBase; using Statistics; using Images
include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))
include(joinpath(@__DIR__, "..", "src", "colocalization.jl"))
```

**Apply to:** `spike/p13/result.jl` (to reach `AbstractColocResult` in `src/results.jl`) and
transitively to `alpha_series.jl` (via `contract.jl`). **The `using` lines must precede the
`include`** — this is a real, documented load-order landmine.

### Pattern S4 — Two-layer RNG discipline
**Source:** `spike/data/seeding.jl:41-47, 75-96`; `spike/validation/harness.jl:60-94`.

```julia
const HOLDOUT_SALT = 0x9E3779B97F4A7C15   # reserved ADVI holdout stream (D-10)
const FOLD_SALT    = 0xD1B54A32D192ED03   # loader k-fold permutation stream (Wave-4)
...
sample_rng(master_seed, idx) = Philox4x(UInt64, (UInt64(master_seed), UInt64(idx)))
val_rng(master_seed = VAL_MASTER_SEED) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ VAL_SALT, UInt64(0)))
```

**Apply to:** `consts.jl` (define `P13_SALT`, `p13_rng(counter)`), `tau_probe.jl`, the label
generator, both runners.
**Plus the layer the analogs do NOT have:** `Random.seed!(derived_from(P13_DEV_SEED))` immediately
before every `train` and every posterior-draw loop — Flux's `DataLoader(shuffle=true)` and
`sampleposterior` use the global RNG and have already caused a real regression in this repo
(13-RESEARCH Pitfall 6). This is a Phase-13 *addition* to the established pattern; name it in the plan.

### Pattern S5 — Atomic artifact write (`.tmp` → integrity check → `mv(force=true)`)
**Source:** `spike/validation/train_ratio.jl:331-347`, `spike/validation/run_bf.jl:59-69`.
**Apply to:** persisted net, label tables, both reported reports, and the frozen `P13_TAU` block.
Never `jldsave` directly to the final path.

### Pattern S6 — Frozen-`zt` discipline (never re-fit a standardizer)
**Source:** `spike/validation/train_ratio.jl:45-47, 212`; `src/amortized/train_ratio.jl:144-146`.

```julia
# FROZEN STATS (SC5 / T-05-01): BOTH summaries in every pair are standardized with the
# FROZEN train `m.zt` (via standardize_summary, the harness surface). NO ZScoreTransform is
# re-fit in this file -- the ratio net must see summaries in the deployed NPE's input space.
...
    Zstd = standardize_summary(pool.summary_min, m.zt, variant)   # 128×N in the frozen space
```

**Apply to:** all Phase-13 labelled-data generation. The `zt` is Phase 11's, inherited via
`preconditions.jl`. Warning sign that it was re-fit: standardized Phase-13 summaries with mean ≈ 0 /
sd ≈ 1 (Pitfall 5).

### Pattern S7 — Pre-registered threshold, referenced never inlined
**Source:** `spike/validation/run_bf.jl:135-138` ("Thresholds referenced from consts.jl (never
tuned)"), `spike/test/test_bf.jl:91-98` (asserted-as-locked in the test suite).
**Apply to:** every number in the D-12/D-13 gate, τ, the α ladder, the iteration allowance.
Two obligations, both present in the analog: the runner *reads* the constant, and the test suite
*asserts its value*.

### Pattern S8 — Honest-failure reporting (never loosen a locked threshold)
**Source:** `spike/validation/run_bf.jl:36-40`, `spike/test/test_bf.jl:29-34`.

```
# ANTI-SNOOPING / HONESTY CONTRACT (D-08, CLAUDE.md HARD_GATE): BF_CORR_MIN/BF_LOGBF_TOL/
# BF_SWEEP_* were committed to consts.jl BEFORE this run. 05-02 reported an HONEST fixture
# finding — order-correct (corr=0.961) but magnitude-inflated (max|Δ|=3.18 > 0.5). Whether
# the reported-scale run clears D-08b is an honest empirical question; a shortfall is a
# Phase-6 finding, NEVER a reason to loosen a locked threshold or reseed.
```

**Apply to:** the Phase-13 report and both runners. Phase 13 additionally has a **pre-declared
single-iteration trigger** (D-04 + 13-RESEARCH §F3: exclusion-class AUC below floor *specifically at
high |ρ|* ⇒ switch to D-07-ii) — that trigger must be in `consts.jl` before the run, which is a
strengthening of this pattern with no existing analog.

### Pattern S9 — Vacuous-pass guard (report AUC beside ECE)
**Source:** the F3 precedent named in `07-CALIBRATION-FINDINGS.md`; mechanically,
`spike/test/test_sbc.jl:63-68` (a deliberately-biased input must FAIL).
**Apply to:** D-13. A green ECE with AUC ≈ 0.5 is a vacuous pass and must be labelled as such
(13-RESEARCH Pitfall 9). Every calibration claim in Phase 13 ships with its discrimination number.

---

## No Analog Found

Nothing in the file list is analog-free, but four **sub-components** have no in-repo precedent and
must be built from 13-RESEARCH's derivations rather than copied:

| Component | File | Why no analog | Planner action |
|-----------|------|---------------|----------------|
| `ThreeWayEvidenceNet <: NeuralEstimator` + `Flux.@layer` + custom loss through `train` | `spike/p13/net.jl` | Every existing net in the repo is a stock `RatioEstimator` / `PosteriorEstimator`; no custom `NeuralEstimator` subtype exists | Use 13-RESEARCH §A4/§B4 verbatim; gate it behind the A5 two-epoch smoke test **before** any real training (route R2, a hand-rolled loop, is the pre-declared fallback) |
| `masked_two_head_bce` (4-row target matrix, participation weights) | `spike/p13/net.jl` | No multi-head or masked loss anywhere in the repo | 13-RESEARCH §B3 gives the exact construction; unit-test that perturbing an *inactive* logit leaves the loss bit-identical |
| Closed-form Gaussian toy verification of the D-07 correction (with negative control) | `spike/test/test_p13.jl` | No analytic-ground-truth verification exists in the repo; every existing check is empirical | 13-RESEARCH §F5 gives runnable code + the pre-registered pass bars; the **negative control** (uncorrected logit must FAIL) is non-optional |
| `alpha_segregate` intensity-preserving redistribution | `spike/p13/alpha_series.jl` | `misspec_*` perturbs; nothing in the repo *conserves* total intensity while moving it spatially | 13-RESEARCH §J3 gives the function and its four testable invariants; the `negctrl_affine` docstring (above) is the reasoning template |

**Two blocked/deferred surfaces (not files to write):** the physical segregated anchor
(`corpus/data/` empty, `sha256 = "PENDING-FETCH"`, `split = sealed_holdout`) — do **not** call
`open_sealed_holdout`; defer to Phase 16 with a named note. And the Phase-11 research NPE bundle —
a hard execution block with no fallback (13-RESEARCH §Environment Availability).

---

## Conventions

**Convention derivation skipped** (`gsd-tools verify conventions --derive --scope spike` returned
`{"skipped": true, "reason": "no-readable-files"}` — the deriver does not cover `.jl` sources). The
table below is therefore **read off the analogs by hand**, not tool-derived; treat shares as
qualitative ("all files read" vs "mixed"), not measured.

| Axis | Dominant | Share | Entropy | Status |
|------|----------|-------|---------|--------|
| File-name casing | `snake_case.jl` (`train_ratio.jl`, `run_bf.jl`, `test_bf.jl`, `alpha_series.jl`) | all files read | low | **named contract** |
| Identifier casing | `snake_case` functions, `SCREAMING_SNAKE` consts, `CamelCase` types (`CalibrationResult`, `AmortizedColocResult`, `MultiChannelImage`) | all files read | low | **named contract** |
| Export style | **no modules, no `export`** in `spike/` — flat top-level functions reached by guarded `include` (`harness.jl:41` "Flat top-level functions (sibling style of infer.jl — no module wrapper)"); `src/` is the opposite (one `ProteinCoLoc` module, `PC.`-qualified) | split by directory | moderate repo-wide, low per-directory | **contested hotspot (author's choice)** |
| Import style | `using X  # what for` for broad surfaces in `spike/`; `import X: sym1, sym2` for narrow ones in `src/amortized/`; every `using`/`import` line carries a trailing comment naming what it is for | split by directory | moderate repo-wide, low per-directory | **contested hotspot (author's choice)** |

**Contested hotspots (author's choice).** The two contested axes are the same **intentional
directory split** the project uses everywhere: `spike/**` is a *script lane* (no module wrapper,
flat functions, guarded `include`, `using` with a purpose comment — see
`spike/validation/harness.jl:41-58`), while `src/**` is a *package lane* (one module, narrow
`import X: sym` lists — see `src/amortized/train_ratio.jl:44-49`, `src/amortized/bf.jl:54-57`). Each
half is internally consistent per-directory and contested only when measured repo-wide. This is the
same shape as the CJS↔SDK dual-resolver prototype: reviewers and planners must **match the
directory's local style**, not impose a repo-wide winner. Phase 13 lives entirely in `spike/`, so
every new file takes the **spike lane**: no module, no `export`, `using X  # comment`, guarded
`include`s, `joinpath(@__DIR__, ...)` paths.

Two further conventions worth stating because Phase 13 leans on them:

- **Docstrings are contracts, not summaries.** Every analog docstring names the decision ID
  (`D-07`, `Pitfall 5`, `SC5 / T-05-01`) and the source line it was promoted from
  (`src/amortized/bf.jl:22-24` lists three spike line ranges). Phase-13 docstrings must cite
  D-05/D-07/D-09/D-11 and the `13-RESEARCH` section that derived them.
- **Comments carry the *why*, especially the traps.** `bf.jl:117-123` spends seven lines explaining
  a `clamp(·,0,1)`. That density is the house style; a bare Phase-13 implementation with thin
  comments would be off-convention.

---

## Metadata

**Analog search scope:** `spike/` (validation, test, data, simulator, npe, contract),
`src/` (amortized/, results.jl, colocalization.jl, LoadImages.jl — read-only),
`test/gate/`, `.planning/spikes/014-bf-sim-validation/`, `.planning/phases/07-*/`.
**Files read this session:** 22 (14 full, 8 targeted ranges).
**Pattern extraction date:** 2026-07-25
**Invalidated by:** Phase 11 landing (re-read its artifact contract and `encode_lambda` form before
`preconditions.jl` is written), or any `spike/Manifest.toml` change.

---

## Amendment: real-image arm (D-15 amended, commit `5b4da6d`)

**Appended:** 2026-07-25 — **narrow amendment. Nothing above this line was edited.**
**Trigger:** `5b4da6d docs(12,13): amend real-image decisions -- corpus anchors are sealed_holdout for
Phase 16`. The amended D-15 adds a **real-microscopy arm** on the six committed TIFFs at
`test/test_images/{positive,negative}/{positive,negative}_c{1,2,3}.tif` (1028×1376). The map above
predates that amendment and contains no mapping for it (`grep -c test_images` over the pre-amendment
file = 0).

**Scope of this amendment:** 2 new files + 1 delta on an already-mapped file + 2 new shared patterns
+ 1 **correction** to a pattern asserted above. Everything else above stands unchanged.

**Hard constraints re-asserted for every analog below.**

- **`corpus/` is OFF LIMITS as a data source.** No Phase-13 executable may reference `corpus/` or call
  `open_sealed_holdout`. `corpus/load.jl:44-45` may be cited **only** as a code-shape example of the
  read-only `src/` include idiom — and even there `spike/contract.jl:30-47` is the better analog, so
  prefer it and never name `corpus/` in code.
- **`src/` stays byte-unchanged** (D-01 research lane). Every `src/` excerpt below is **[READ-ONLY
  REF]** — copy the *call*, never the file.
- **No new packages.** `Images` 0.26.2 already re-exports `otsu_threshold`; `CairoMakie`,
  `Statistics`, `StatsBase`, `JLD2` are already resolved. `spike/Project.toml` /
  `spike/Manifest.toml` stay byte-identical.

**Sequencing note the planner must carry (13-RESEARCH §J4.6).** The whole real-image arm —
ingestion, path resolution, the α-ladder on real pairs, every invariant test — is
**Phase-11-independent** and belongs in **Wave 0**. It was executed end-to-end against the shipped
grid-8 bundle this session. Only the final *three-way* read blocks on the Phase-11 net. Do not
schedule this arm behind the D-02 block.

---

### Amended File Classification (append these rows to the table above)

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `spike/p13/real_images.jl` | adapter / loader (real fixtures → summary) | file-I/O → transform | `spike/simulator/calibration.jl:59, 156-164` (`_real_anchor`) + `spike/contract.jl:30-70` | **exact** |
| `spike/p13/alpha_series.jl` (**real-image delta only**) | utility (image transform) | transform | `spike/validation/ood.jl:82-83` (`frozen_summary`, the substrate-agnostic boundary) + `src/LoadImages.jl:235-241` **[READ-ONLY REF]** | role-match |
| `spike/p13/run_real_check.jl` | script (**reported, NOT gated**) | batch → artifact + figure | `spike/validation/run_ood.jl` (whole file) + `11-10-PLAN.md:161-195` (posture, plan-level) | role-match (composite) |

---

### `spike/p13/real_images.jl` (adapter, file-I/O → transform)

**Primary analog: `spike/simulator/calibration.jl:156-164` — this is an *exact* analog.** It already
loads **these exact six files**, from a `spike/<subdir>/` depth, through the frozen loader, and feeds
them to the summary. The Phase-13 loader is this function generalized, not a new capability.

```julia
# spike/simulator/calibration.jl:156-164  --- FAITHFUL real anchor (Phase-02 D-16)
function _real_anchor()
    load_one(cond) = begin
        c1 = load_tiff(joinpath(_REPO_ROOT, "test", "test_images", cond, "$(cond)_c1.tif"))
        c2 = load_tiff(joinpath(_REPO_ROOT, "test", "test_images", cond, "$(cond)_c2.tif"))
        induced_mu(build_mci([Matrix{Float64}(c1), Matrix{Float64}(c2)]; name = cond))
    end
    return (positive = load_one("positive"), negative = load_one("negative"))
end
```

**Copy:** the `load_one(cond)` closure keyed on the condition string, the
`joinpath(_REPO_ROOT, "test", "test_images", cond, "$(cond)_c$(i).tif")` path builder, the explicit
`Matrix{Float64}(...)` conversion before `build_mci`, and the `NamedTuple` return keyed
`(positive = …, negative = …)`.
**Lineage the plan should state (13-RESEARCH §J4.0):** `spike/simulator/ghat.jl:39` carries the frozen
comment `real anchor (D-16): faithful LoadImages.jl load_tiff (NOT luminance): positive μ = 0.3292,
negative μ = 0.2481`, and both numbers reproduced **exactly** this session. Phase 13 therefore runs on
the same files, through the same loader, that the frozen `ghat` calibration used. That is real
lineage — say so; it is a far stronger position than "we picked some TIFFs".
*(Naming collision to avoid: Phase **02**'s D-16 is the real-anchor decision; Phase **13**'s D-16 is
the α-series. Different decisions, same label.)*

**Path resolution — copy this, do NOT copy the other one.**

```julia
# ✅ CORRECT — spike/simulator/calibration.jl:59. Depth from spike/p13/ is IDENTICAL.
const P13_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))     # spike/p13/x.jl -> spike/ -> root
real_tif(cond, i) = joinpath(P13_REPO_ROOT, "test", "test_images", cond, "$(cond)_c$(i).tif")
```

```julia
# ❌ DO NOT COPY — test/runtests.jl:18-21. Correct for Pkg.test(), WRONG for a spike runner:
#    mutating process CWD inside a reported run is a side effect the harness does not expect.
# Test fixtures are referenced by package-root-relative paths (e.g.
# "test/test_images/..."). Under `Pkg.test()` the process CWD is not guaranteed to be the
# package root, so anchor it explicitly to `<pkgroot>` (the parent of this test dir).
cd(dirname(@__DIR__))                      # test/runtests.jl:21
```

Add `@assert isfile(real_tif(c, i))` with a message naming the expected path (the explicit-error
house style, `src/amortized/api.jl:111-120`). Paths are always `joinpath(@__DIR__, ...)`, never
interpolated, never relative — Pattern S2 above, and the V5 path-traversal control.

**The constructor and the read chain** (`src/LoadImages.jl:432-455`, `:214-217` **[READ-ONLY REF]**):

```julia
function MultiChannelImage(name::S, path::Vector{S}, channels::Vector{S} =[]) where {S <: AbstractString}
    data = load_tiff.(path)
    if isempty(channels)
        @warn "No channel names provided, using default names"       # :436 — ALWAYS pass names
        channels = ["channel_$i" for i in 1:length(data)]
    end
    if length(data) != length(channels)
         throw(ArgumentError("Number of channels must match the number of image files"))
    end
    pixel_size = size(data[1])                       # :444 — PIXEL DIMS (1028,1376), not micrometres
    otsu_threshold = Images.otsu_threshold.(data)    # :445 — per-channel scalar, stored on the struct
    return MultiChannelImage(data, channels, name, path, pixel_size, otsu_threshold)
end

load_tiff(path::S) where {S<:AbstractString} = Float64.(Images.Gray.(Images.load(path)))   # :214-217
```

**Argument order is `(name, paths, channels)` — not `(paths, channels, name)`.** Always pass channel
names, or the `@warn` at `:436` lands in the report log.

**What `load_tiff` actually yields — measured, not assumed (13-RESEARCH §J4.1):** on-disk element
type is `RGB{N0f16}` (each `_cN.tif` is itself a 16-bit RGB file); `Images.Gray` applies **Rec-601
luma** (`0.299R + 0.587G + 0.114B`); `Float64.` gives `Matrix{Float64}` in **[0,1]** — **not**
z-scored, **not** background-subtracted. `positive_c1` measures min `6.1e-5`, max `0.198`, mean
`0.0150` — the fixtures are *dim*. Exact zeros: `positive_c1` 0, `positive_c2` **2**, `negative_c2`
**3** (this is the fact that falsifies the invariant corrected below).
`spike/NOTES.md:239`'s "**NOT** a hand-rolled RGB→luminance reduction" means *use the package's own
`load_tiff`*, not *avoid luminance*. Do not "improve" it by summing or maxing channels — that
silently moves off the frozen anchor.

**Read-only reach into `src/` — Pattern S3 above, restated because this file is the one that needs
it most** (`spike/contract.jl:30-47`, already includes both files this loader needs):

```julia
# ORDER MATTERS: StatsBase + Statistics before src/colocalization.jl (correlation() builds its
# method Dict at CALL time); Images before src/LoadImages.jl (the constructor calls otsu_threshold).
using StatsBase; using Statistics; using Images
include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))     # MultiChannelImage(Stack), load_tiff
include(joinpath(@__DIR__, "..", "src", "colocalization.jl")) # patch, correlation, _exclude_zero
```

From `spike/p13/` the include path is `joinpath(@__DIR__, "..", "..", "src", ...)`, **or** — strongly
preferred — reach it transitively with the guarded include of the existing contract file:

```julia
isdefined(@__MODULE__, :build_mci) || include(joinpath(@__DIR__, "..", "contract.jl"))
```

`corpus/load.jl:44-45` uses the same read-only include shape; it is cited here **only** as
corroboration that this idiom is the house standard. **Do not include, read, or name `corpus/` from
any Phase-13 file.**

**MCI construction from raw matrices** (`spike/contract.jl:60-70`) — the shape to mirror when
rebuilding an MCI after the α transform:

```julia
function build_mci(data::Vector{Matrix{Float64}}; name::String = "sim")
    @assert length(data) == 2 "contract expects exactly 2 channels, got $(length(data))"
    return MultiChannelImage(data, ["ch1", "ch2"], name, ["", ""],
        size(data[1]),                 # D-11: pixel-dim tuple, not micrometers
        Images.otsu_threshold.(data),  # mirrors src/LoadImages.jl:445
    )
end
```

**`patch()` SILENTLY TRUNCATES — do not plan a crop or resize step** (`src/colocalization.jl:37-43`
**[READ-ONLY REF]**):

```julia
function patch(img::AbstractMatrix{T}, num_patches::Integer) where T <: Union{Float64, Missing}
    rows, cols = size(img)
    patch_size_x = rows ÷ num_patches                                 # :39  integer division
    patch_size_y = cols ÷ num_patches                                 # :40
    trimmed = @view img[1:num_patches*patch_size_x, 1:num_patches*patch_size_y]      # :43
```

Measured at 1028×1376 (13-RESEARCH §J4.4, [VERIFIED]): `1376 ÷ 8 = 172` exactly; `1028 ÷ 8 = 128.5 →
128`, so **4 of 1028 rows (0.39 %) are dropped**, deterministically, always the same 4. Each patch
holds **22 016 px** — three orders of magnitude above the ≤15-survivor `missing` floor; measured
survivors 22 015–22 016, and **0/64 patches went `missing` on either fixture**.
**Planner consequences:** (1) **no crop, no resize, no pre-registration of a resize rule** — adding
one would be a new preprocessing choice under D-04 and would diverge from the frozen `_real_anchor()`
path; (2) the report must state the truncation in one sentence, because a reader who computes
`1028/8` will expect 128.5-px patches; (3) the `missing`/mask-row machinery is effectively **inert**
on this substrate — assert that, do not design around it.
**Footnote worth one line:** `IMSIZE_SET` (`spike/data/seeding.jl:52-58`) already contains
`(1376, 1028)` at weight 0.03, annotated `# ~21.6× -- real-data anchor (D-08)`. The simulator emits
1376 rows × 1028 cols while `load_tiff` yields 1028 × 1376; per-patch pixel count and truncation are
identical, only the column-major `vec()` ordering corresponds to a transposed layout. Statistically
immaterial (isotropic fields, symmetric `SHIFT_PRIOR`) — a footnote, not a silent difference.

**NO MASKING, NO BACKGROUND SUBTRACTION on the amortized path** (13-RESEARCH §J4.3, verified call
chain `src/amortized/api.jl:107-148`). `_apply_mask!` / `apply_mask!` appear **only** in
`src/plot.jl:114-115,224` and `src/utils.jl:86` (the legacy Turing/plotting lane).
`_apply_mask!` multiplies sub-threshold pixels by `0` (`LoadImages.jl:261`) and `_exclude_zero` then
**deletes** them — masking before the summary would change the statistic the net was trained on.
**Phase 13 must not mask.** Otsu is used only to *define* the α-series object mask, never to modify
pixels that feed `patch_summary`.

**The summary boundary — arity fork the planner must pick deliberately:**

| Call | Location | Grid |
|---|---|---|
| `patch_summary(mci)` | `spike/contract.jl:85-90` (spike lane) | **hard-wired G = 8** |
| `patch_summary(mci, G)` | `src/amortized/summary.jl:50-55` **[READ-ONLY REF]** | grid-general |

Phase 13 runs at G = 8, so the spike-lane 1-arg call is sufficient and keeps the file in the spike
lane (Conventions above). If a G-parameterized call is wanted, take the 2-arg `src/` one read-only —
do not add a third definition.

---

### `spike/p13/alpha_series.jl` — **REAL-IMAGE EXTENSION ONLY**

> This file is already fully mapped above (Analogs A–E, simulator path). **Everything there stands.**
> Only the delta is below, plus one **correction**.

#### ⚠ CORRECTION to the map above (and to `13-04-PLAN.md`)

The invariant **`all(y′ .> 0)`** — asserted at **`13-PATTERNS.md:548`** ("assert `all(y' .> 0)`") and
required at **`13-04-PLAN.md:229`** and **`13-04-PLAN.md:270`** ("`all(z .> 0)`") — **is FALSE on the
real fixtures and must not be written.**

The source TIFFs already contain exact-zero pixels **before any α**: `positive_c2` has **2** and
`negative_c2` has **3** (of 1 414 528), so `minimum(y′) = 0.0` at every α *including α = 0, where
`y′` is bitwise `y`* [VERIFIED, 13-RESEARCH §J5.3 finding 4 / Pitfall 2b].

```julia
# ❌ WRONG (fails Wave 0 on real input for a reason unrelated to the algorithm):
@assert all(y_prime .> 0)

# ✅ CORRECT — the zero SET is what matters, not zero-freeness:
z0 = count(iszero, y)                              # measured: 2 (positive_c2), 3 (negative_c2)
@assert count(iszero, y_prime) == z0               # no NEW zero created
#   …and report z0 in the artifact + report, so the reader knows the source zeros exist.
```

The *reasoning* behind the original invariant is still exactly right, and the analog that carries it
(`negctrl_affine`'s docstring, `spike/validation/ood.jl:513-528`, quoted in the map above) remains the
best template: **changing which pixels are zero changes the summary via `_exclude_zero`.** The fix is
to assert *zero-set preservation* rather than *zero absence*. `13-04-PLAN.md`'s acceptance criteria
at `:229` and its invariant list at `:270` must be amended before that plan executes.
`13-04-PLAN.md:88-93`'s `alpha_segregate` body itself is correct and needs no change; the `b > 0` and
`S_out > 0` guards stay.

#### The interface boundary that makes ONE code path serve both substrates

**Answer: `Vector{Matrix{Float64}}` (the pair) / `AbstractMatrix{Float64}` (one channel).** Both
producers already emit exactly this — `simulate_pair` returns `Vector{Matrix{Float64}}`
(`spike/simulator/forward.jl:103`) and `load_tiff.(path)` returns the same
(`src/LoadImages.jl:433`) — and the existing scorer is already typed on it:

```julia
# spike/validation/ood.jl:82-83 — THE substrate-agnostic boundary, already shipped.
# "Put a 2-channel external image `pair::Vector{Matrix{Float64}}` into the estimator's input space
#  through the FROZEN read chain … so a misspecified image is scored byte-identically to training
#  data with NO re-fitting (T-05-01 / SC5)."
frozen_summary(m, pair) =
    standardize_summary(encode_d01(patch_summary(build_mci(pair))), m.zt, :min)
```

```
     ┌─ simulate_pair(rng, θ; imsize) ──┐
     │      (spike/simulator)           │
[x, y] :: Vector{Matrix{Float64}} ──────┼──► alpha_segregate(y, M, α; b) ──► build_mci([x, y′])
     │                                  │         └─► patch_summary ─► encode_d01 ─► frozen zt
     └─ load_tiff.(paths)[1:2] ─────────┘
              (src/LoadImages)
```

**Therefore:** type `alpha_segregate` on `AbstractMatrix{Float64}` and give it **no knowledge of
provenance**. Put the three substrate-dependent knobs in a small provenance `NamedTuple`/struct
carried alongside, never inside the algorithm (13-RESEARCH §J5.4):

| Knob | Simulated | Real |
|---|---|---|
| background floor `b` | use the **same** `quantile(vec(y), 0.05)` rule for comparability (not `BG_FLOOR = 0.02`) | `quantile(vec(y), 0.05)` ≈ 0.0062 / 0.0059 [measured] |
| frame size | drawn from `IMSIZE_SET` | fixed 1028 × 1376 |
| **α = 0 semantics** | a *genuinely random* pair by construction (drawn at ρ ≈ 0) | **NOT random** — measured `m̄ = +0.329 / +0.248` |

That last row is the one real asymmetry, and it is a **reporting** obligation, not a code branch: on
simulated substrate the ladder runs **random → exclusion**; on real substrate it runs **moderately
colocalized → near-exclusion**. Both are informative; **the report must not average them into one
curve.**

#### Mask source — already frozen, already shipped, no new dependency

`Images.otsu_threshold` is re-exported by **`Images` 0.26.2, already in `spike/Manifest.toml`** — it
is the house standard, called at `src/LoadImages.jl:445`, `src/amortized/api.jl:58`,
`src/amortized/local_map.jl:208`, `src/amortized/simulator.jl:268`, `spike/contract.jl:68`.
**No `Pkg.add`.** `ImageSegmentation` / `ImageMorphology` are transitively present but **must not be
added** to `spike/Project.toml` — the resolve-risk gate (`spike/test/runtests.jl` clauses (d)–(h))
exists to catch exactly that.

The rule itself is Analog D above (`src/LoadImages.jl:235-241` `_calculate_mask`, **[READ-ONLY REF]**)
— unchanged, still the right choice, and it makes D-16's segmentation a *frozen, already-shipped*
decision at zero D-04 cost. On a real MCI the accessor path is the cheapest form:

```julia
M = _calculate_mask(mci)[1]           # == mci.data[1] .> mci.otsu_threshold[1]
#   the constructor ALREADY populated otsu_threshold (LoadImages.jl:445) — no recomputation.
#   Compute ONCE from the UNMODIFIED ch1 and hold FIXED across the whole α ladder
#   (the local_map.jl:208 "tile-local Otsu makes tiles incomparable" precedent).
```

**Measured mask behaviour** (13-RESEARCH §J5.1): positive Otsu `0.028661` → mask **13.49 %**;
negative Otsu `0.035710` → mask **22.92 %**. Both leave a large complement, so the
"mask-covers-the-frame" failure mode does not occur here. Pre-register the guard anyway:
`@assert 0.01 <= mean(M) <= 0.40`, plus `@assert b > 0` (b = 0 re-opens the `_exclude_zero` trap) and
`@assert maximum(y′) <= 1.0` recorded, **never silently clamped**. Measured α = 1 boosts: **1.509**
(positive) / **1.598** (negative) — comfortable.

#### The measured real ladders — cite these, do not re-derive them

13-RESEARCH §J5.3 executed the full ladder on both fixtures at G = 8 [VERIFIED 2026-07-25]:

| α | positive `m̄` → `ghat(m̄)` | negative `m̄` → `ghat(m̄)` | `missing` | Σ preserved |
|---|---|---|---|---|
| 0.00 | +0.3292 → +0.3305 | +0.2481 → +0.2149 | 0/64 | bitwise identity |
| 0.50 | +0.0869 → +0.0005 | −0.0187 → −0.1450 | 0/64 | ✓ (rtol 1e-10) |
| 1.00 | −0.1679 → −0.3470 | −0.2520 → −0.4419 | 0/64 | ✓ |

Sign crossing at **α ≈ 0.49** (positive) / **α ≈ 0.46** (negative) — well inside the ladder, which is
what makes `α*` informative. **Scope limit to state:** the ladder spans `ρ_true ∈ [−0.44, +0.33]`, a
densely-covered interior region — **it does NOT probe the deep-exclusion tail near the ±0.99 clamp
atoms** (§E2). A reader who assumes otherwise will over-read the result.
**Assert `missing`-count α-invariance** (a rise with α means the zero-free property broke), and
assert the mask is byte-identical across arms.

---

### `spike/p13/run_real_check.jl` (script, **REPORTED — NOT GATED**)

**Analog A (code shape, primary): `spike/validation/run_ood.jl` — the whole file.** It is the repo's
existing "score *external* inputs through the frozen chain → persist a `.jld2` → write a figure →
print tables" runner. Reuse this rather than inventing a shape; the map above already assigns it to
`run_alpha_series.jl`, and `run_real_check.jl` is its sibling.

**Copy: resolution knobs explicitly labelled NOT anti-snooping thresholds** (`run_ood.jl:70-76`) —
this is the mechanism that keeps a reported arm from mutating into a gate:

```julia
# Monte-Carlo resolution knobs for the reported grid (NOT anti-snooping thresholds —
# the pass/fail constants OOD_AUC_MIN/OOD_KS_EPS/OOD_ID_QUANTILE/OOD_GRID_LEVELS are
# locked in consts.jl). Sized for a stable reported AUC.
const OOD_N_ID  = 100
const OOD_N_POS = 40
```

**Copy: guarded includes + report path const** (`run_ood.jl:63-68`):

```julia
# Units under test: the OOD surface (pulls ood.jl → harness.jl → consts.jl transitively)
# and the figure surface (CairoMakie — figures are script artifacts, not gate assertions).
isdefined(@__MODULE__, :ood_roc_over_grid) || include(joinpath(@__DIR__, "ood.jl"))
isdefined(@__MODULE__, :plot_roc)          || include(joinpath(@__DIR__, "figures.jl"))
const OOD_REPORT_PATH = joinpath(@__DIR__, "ood_report.jld2")
```

**Copy: atomic write with a keyed integrity check** (`run_ood.jl:78-87` — the same Pattern S5 shape as
`run_bf.jl:59-69`; here the integrity key becomes a real-arm key, e.g. `"m_bar"` or `"alpha_grid"`):

```julia
function _save_ood_report(path; kwargs...)
    mkpath(dirname(path)); tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "combined_auc") "_save_ood_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)
    return path
end
```

**Copy: persist-before-anything-can-throw, with the pre-registered constants written INTO the
artifact** (`run_ood.jl:176-198`) — including `generated = string(Dates.now(Dates.UTC)) * "Z"`.

**Change vs. the analog — and this is the load-bearing difference: `run_real_check.jl` has NO
`@testset`, NO `@test`, and NO threshold of any kind.** `run_ood.jl:246-265` ends in a real gate;
delete that block entirely. 13-RESEARCH §J6.3 gives the strongest form of the constraint: *"No
pass/fail threshold is defined for any real-image quantity"* — if no threshold exists, it cannot
accidentally become a gate.

---

**Analog B (posture / file structure): `11-10-PLAN.md:161-195`** — the Phase-11 D-17 real-image
check, one phase earlier, same fixtures, same honesty problem. **Mirror this file structure rather
than inventing one**; reusing it makes the two phases' honesty claims verifiably consistent. It is a
*plan-level* analog (the file `spike/validation/run_p11_realimage.jl` does not exist yet), so copy the
four structural requirements, not code:

1. **ALL-CAPS banner** stating (a) this is CHARACTERIZATION, not a pass/fail gate, and explicitly why
   it *cannot* be a correctness test (no ground-truth ρ on real data);
2. a written **target-substitution record** — D-15 originally named the corpus anchor; `corpus/data/`
   is empty, both physical anchors are `sha256 = PENDING-FETCH`, `bytes = 0`, `split =
   sealed_holdout`, reserved for Phase 16, so opening them would irreversibly burn Phase 16;
   `test/test_images/` holds six real committed 1028 × 1376 TIFFs already exercised by
   `test/runtests.jl:86,166`;
3. **read-only discipline enforced structurally** — record a `git ls-files -s test/test_images/`
   digest at the start of the run and re-check it at the end, asserting it unchanged;
4. a **source-text assertion**: `read(@__FILE__, String)` must not contain `"corpus/"` outside the
   explanatory header, and `open_sealed_holdout` must never appear — so the abstention claim is
   *machine-checked*, not merely intended.

Also carry the `src/`-untouched step-0 check already assigned to `run_three_way_gate.jl` above
(`@assert success(git diff --quiet HEAD -- src/)`, the `spike/demo.jl` precedent).

---

**Analog C (OOD call surface) — two surfaces, pick per lane; the runner must print the verdict beside
every real result.**

*Spike lane (D-01 research lane, preferred for the Phase-13 net)* — `spike/validation/ood.jl:112-133`:

```julia
function fit_ood_nulls(Ztrain::AbstractMatrix; variant::Symbol = :min)
    cont, _ = Loader._row_partition(variant, size(Ztrain, 1))
    Zc = Float64.(Ztrain[cont, :])
    μS = vec(mean(Zc; dims = 2)); ΣS = cov(Zc; dims = 2)
    C  = cholesky(Symmetric(ΣS + 1e-6 * I))     # ridge 1e-6 (Pitfall 2)
    return (cont = cont, μS = μS, C = C)
end

maha_score(nulls, z::AbstractVector) = sum(abs2, nulls.C.L \ (Float64.(z[nulls.cont]) .- nulls.μS))
```

Threshold construction is the pre-registered ID quantile (`ood.jl:636` `maha_thr =
id_threshold(id_maha)` at `OOD_ID_QUANTILE`, fit on a **held-out ID pool** — never on the real
images). Real-image scoring goes through `frozen_summary(m, pair)` (`ood.jl:82-83`), i.e. **the exact
same call the misspecification families use** — a real pair is just another external input.

*Shipped lane (to reproduce the §J4.6 numbers)* — `src/amortized/ood.jl:354-388` **[READ-ONLY REF]**,
invoked as at `src/amortized/api.jl:143`:

```julia
ood = ood_verdict(b.ood_nulls, Zs, Zc)      # api.jl:143 — density-channel OOD flag
# returns OODVerdict(score, flag, per_channel); flag = score > ood_nulls.thr  (ood.jl:385-387)
```

**MEASURED, and this is the single most consequential finding of the amendment (13-RESEARCH §J4.6 /
§J6.2): the shipped OOD detector FLAGS BOTH FIXTURES.** Mahalanobis density **433.688** against the
frozen in-distribution threshold **`b.ood_nulls.thr = 179.1368`** — **2.4× over**, `is_ood = true` in
both read directions. This is not a bug and not a reason to skip the check; it is the OOD channel
doing its job, and it is arguably the most honest single number the arm will produce.
**Requirement: every real-image result in the artifact, in the printed table, and in the report is
reported next to its OOD score, the threshold, and its λ.** The Phase-13 net is a *different* net and
may score differently — measuring and reporting that comparison is itself a deliverable.
Also reproduce the §J4.6 sanity reads as a smoke check: `logBF = +1.988` (positive as sample) /
`−1.085` (negative as sample), mean `ρ_sample` 0.302 / 0.182 vs. the independent `ghat` anchors
0.331 / 0.215.

---

**Analog D (figures): `spike/validation/figures.jl:36-45` + the call sites at `run_ood.jl:200-202`.**
The repo is **CairoMakie**, headless, already resolved — not Plots.

```julia
using CairoMakie
CairoMakie.activate!()   # headless raster/vector backend (no OpenGL/display)

# --- shared output directory (created on demand) --------------------------------
const VAL_FIG_DIR = joinpath(@__DIR__, "figures")
_val_fig_path(name) = (isdir(VAL_FIG_DIR) || mkpath(VAL_FIG_DIR); joinpath(VAL_FIG_DIR, name))
```

Every plotter follows the same three lines: build `Figure(size = (…, …))` → `CairoMakie.Axis(fig[1,1];
title, subtitle = <CAPTION const>, xlabel, ylabel)` → `path = _val_fig_path(filename); save(path,
fig); return path` (`figures.jl:56-70, 78-88, 97-112`). **Convention: figures live in
`<runner-dir>/figures/`** — `spike/validation/figures/` for the Phase-5 runners
(`run_ood.jl:202` prints `figure -> $(joinpath(@__DIR__, "figures"))/…`), so Phase 13 writes to
**`spike/p13/figures/`** via an identical `_p13_fig_path` helper. (`spike/figures/` holds the Phase-2
`02_simulator_demo.jl` output — a different runner's directory, not a global figure dump.)
PNGs are gitignored regenerable artifacts. **Figures are script artifacts, never gate assertions**
(stated identically at `figures.jl:31` and `spike/test/test_bf.jl:37`); `test_p13.jl` must call none
of them.
Subtitle every real-image figure with the honesty caption (a `P13_REAL_CAPTION` const in
`consts.jl`, mirroring `SBC_CAPTION`): *qualitative, no ground-truth label, OOD-flagged*.

---

**λ on a real image — a decision the plan must make explicit, not a code detail** (13-RESEARCH §J4.5).
A real image has **no known λ**; there is no correct value to plug in. Recommended and to be
pre-registered: **(a) report the whole λ-sweep across the Phase-11 ladder and headline
`λ = LAMBDA_MAX = 3.0`** — a conservative λ makes the evidence *weaker*, so an exclusion verdict that
survives it is the credible one. `11-10-PLAN.md` already reads at two λ levels rather than one, so
this matches the sibling phase. **Never read a real image at an invented single λ** (Pitfall 12).

---

### Amended Shared Patterns (append to `## Shared Patterns` above)

#### Pattern S10 — Reported-not-gated real-data arm (the honesty envelope)
**Source:** `11-10-PLAN.md:161-195` (posture, plan-level); mechanically `run_ood.jl:70-76` (knobs are
not thresholds) and `run_ood.jl:47-53` (NAMED BLIND SPOT banner).
**Apply to:** `spike/p13/run_real_check.jl`, the real-image half of `run_alpha_series.jl`, the phase
report header, and every plan's acceptance criteria.

Five obligations, all machine-checkable:
1. ALL-CAPS banner: QUALITATIVE / REPORTED-NOT-GATED, and *why* it cannot be a correctness check.
2. **No pass/fail threshold exists for any real-image quantity** — the strongest form of (1).
3. Target-substitution record naming `corpus/` and why it was not used (sealed for Phase 16).
4. `read(@__FILE__, String)` must not contain `"corpus/"` outside the header; `open_sealed_holdout`
   never appears anywhere in `spike/p13/`.
5. Read-only digest of `git ls-files -s test/test_images/` asserted unchanged across the run.

Report wording is pre-drafted in 13-RESEARCH §J6.3 (named limit, report header, acceptance criteria)
— copy it rather than re-authoring, so Phases 11/13 stay verbally consistent.

#### Pattern S11 — Every real-image number ships with its OOD verdict and its λ
**Source:** `src/amortized/api.jl:143` + `src/amortized/ood.jl:385-387`; the measured
433.69 vs. 179.14 flag (13-RESEARCH §J4.6).
**Apply to:** every printed table, every JLD2 field, every figure caption in the real arm.
This is the same discipline as Pattern S9 (never a calibration number without its discrimination
number): never a real-image evidence number without its misspecification flag.
State the two framing facts alongside: (i) `positive`/`negative` are **biological conditions, not
colocalization labels** — the "negative" pair measures `m̄ = +0.248`; (ii) `spike/simulator/ghat.jl:40`
already records *"real anti-correlation never observed … documented PRIOR-ONLY"*, so **the exclusion
hypothesis has no observed real-data instance anywhere in this project**. The α-series produces the
first negative induced μ this project has ever obtained from real pixels (−0.35 / −0.44) — it does
not refute the Phase-02 verdict, it sharpens it.

---

### Amended "No Analog Found" (append to the table above)

| Component | File | Why no analog | Planner action |
|-----------|------|---------------|----------------|
| A runner with **zero** gate assertions | `spike/p13/run_real_check.jl` | every existing `run_*.jl` (`run_bf`, `run_ood`, `run_sbc`) ends in a real `@testset` that exits nonzero; none is purely descriptive | Take `run_ood.jl` and **delete** its `@testset` block; enforce absence with Pattern S10 obligation (2) — a grep for `@test` in the file must return 0 |
| λ-sweep read on an **unlabelled** real image | `spike/p13/run_real_check.jl` | the only precedent (`11-10-PLAN.md`) is a **plan**, not shipped code — `spike/validation/run_p11_realimage.jl` does not exist yet | Copy the *structure* from `11-10-PLAN.md:161-195`; if Phase 11 lands first, re-read its realized file and prefer it as the code analog |
| Substrate-provenance record travelling with an α ladder | `spike/p13/alpha_series.jl` | nothing in the repo carries "where did this pair come from" metadata alongside a transform | Small `NamedTuple` per 13-RESEARCH §J5.4 (`substrate = :real|:simulated`, `b`, `imsize`, `mask_fraction`, `source_zero_count`); persist it in the artifact so the two ladders are never averaged |

**Unchanged from above:** `alpha_segregate` itself still has no in-repo analog (nothing *conserves*
total intensity while moving it spatially); `negctrl_affine`'s docstring remains the reasoning
template. **Still blocked, still deferred:** the physical segregated anchor — `corpus/data/` empty,
`sha256 = PENDING-FETCH`, `split = sealed_holdout`. **Do not call `open_sealed_holdout`.**

---

### Amendment metadata

**Files read for this amendment:** `13-CONTEXT.md` (full), `13-RESEARCH.md` §J3/§J4/§J5/§J6 +
Code Example 5, `13-PATTERNS.md` (targeted ranges: header, classification, `alpha_series`, both
runners, shared patterns, conventions), `src/LoadImages.jl:205-269,425-455`,
`src/colocalization.jl:25-60`, `src/amortized/api.jl:100-148`, `src/amortized/ood.jl:330-388`,
`src/amortized/summary.jl:44-56`, `spike/contract.jl:25-104`, `spike/simulator/calibration.jl:40-90,
140-184`, `spike/validation/run_ood.jl` (full), `spike/validation/ood.jl:70-140`,
`spike/validation/figures.jl:36-115`, `test/runtests.jl:1-189`,
`11-10-PLAN.md:140-205`, `13-04-PLAN.md` (invariant lines).
**New analogs mapped:** 3 files, 11 concrete excerpts, 2 new shared patterns, 1 correction.
**Amendment date:** 2026-07-25
**Invalidated by:** Phase 11 landing (prefer its realized `run_p11_realimage.jl` over the plan-level
analog, and re-read `encode_lambda`/`LAMBDA_MAX` before locking the λ sweep), any change to
`test/test_images/` (the digest assertion will catch it), or any `spike/Manifest.toml` change.
