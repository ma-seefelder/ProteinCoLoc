# Phase 11: Registration and Chromatic Uncertainty as Latent — Pattern Map

**Mapped:** 2026-07-25
**Files analyzed:** 22 (14 new, 8 modified)
**Analogs found:** 22 / 22 (every new file has an in-repo shape to copy)

> **Scope of this document.** `11-RESEARCH.md` already owns the *mechanism* analysis (what the
> composed affine must be, why the probe must be paired, where the θ-arity ripples). This file
> owns only the **file shape / house style**: for each file, the closest existing analog, the
> verbatim structural excerpt to imitate, the conventions it encodes, and the naive-pattern-match
> trap. **Do not re-derive mechanism here — read RESEARCH §A–§G for that.**

---

## File Classification

### NEW files

| New file | Role | Data flow | Closest analog | Match |
|---|---|---|---|---|
| `spike/validation/p11_consts.jl` | config / pre-registration | (declarative) | `spike/validation/consts.jl` (+ `test/gate/gate_consts_8_v2.jl` for the seed-derivation + executable-rule half) | **exact** |
| `spike/test/test_p11_consts.jl` | test (unit) | assertion | `spike/test/test_sbc.jl` §"SC3 (SBC-03) pre-registration + disjoint seed" | **exact** |
| `spike/test/test_stage6_regression.jl` | test (regression) | file-I/O (golden fixture) | `spike/test/test_simulator.jl` §"SIM-01 forward pipeline" | **exact** |
| `spike/test/test_p11_smoke.jl` | test (smoke) | request-response | `spike/test/runtests.jl` §"NeuralEstimators CPU smoke" | **exact** |
| `spike/test/test_lambda_ablation.jl` | test (behavioural tripwire) | transform | `spike/test/test_sbc.jl` §SC2 clause (b) "deliberately biased … exercises BOTH directions" | role-match |
| `spike/test/test_p11_tost.jl` | test (unit, pure fn) | transform | `test/gate/gate_consts_8_v2.jl:354-390` (`holm_adjusted`/`holm_pass`) + its `runtests.jl` twin-agreement test | **exact** |
| `spike/test/test_p11_forced_theta.jl` | test (unit, invariant) | transform | `test/runtests.jl:534-563` (θ-bounds / `collect(values(...))` row-order testset) | **exact** |
| `spike/test/fixtures/p11_stage6_golden.jld2` | fixture (artifact) | file-I/O | `spike/npe/train_npe.jl:160-180` `save_npe` `.tmp`→reopen-assert→`mv` idiom | role-match |
| `spike/validation/run_p11_probe.jl` | script (experiment) | batch | `spike/validation/run_sbc.jl` | **exact** (mode differs — see §Experiment-script shape) |
| `spike/validation/run_p11_ladder.jl` | script (experiment) | batch | `spike/validation/run_sbc.jl` | **exact** |
| `spike/validation/run_p11_breakdown.jl` | script (experiment, reported-not-gated) | batch | `spike/npe/scaling.jl` header voice + `run_ood.jl` main() body | role-match |
| `spike/validation/run_p11_attenuation.jl` | script (paired A/B experiment) | batch | `spike/npe/ablation.jl` (the "single variable held identical" arm design) | role-match |
| `spike/validation/run_p11_realimage.jl` | script (experiment, read-only external data) | file-I/O | `spike/validation/run_ood.jl` (external-input scoring via `frozen_summary`) | role-match |
| `.planning/phases/11-…/11-REPORT.md` (name TBD) | report artifact | (document) | `.planning/phases/07-…/gate-8x8-amended.md` | **exact** |
| `BoundedThetaTransform` port into `spike/` (if needed) | model / transform | transform | `src/amortized/architecture.jl:43-104` | **exact** |

### MODIFIED files

| Modified file | Role | Data flow | In-file pattern to preserve |
|---|---|---|---|
| `spike/simulator/prior.jl` | model (prior) | transform | banner-comment block + one-line-per-prior `const` column, `sample_prior` NamedTuple field order |
| `spike/simulator/forward.jl` | simulator (stage pipeline) | transform | numbered `# --- (n) stage ---` comment spine; entry-guard `\|\|  throw(ArgumentError(...))` wall |
| `src/amortized/simulator.jl` | model (promoted mirror) | transform | `#####` banner; "copied UNCHANGED from the spike" provenance sentences; `theta_prior_bounds()` derive-never-literal rule |
| `src/amortized/datagen.jl` | data pipeline | batch | `DATAGEN_HASH_SRC_FILES` fixed-order list (do **not** reorder); literal `7`→`8` at `:170,:274,:393` |
| `src/amortized/ood.jl` | service | request-response | `_theta_tuple(v)` positional θ reconstruction (`:124-132`) must grow an 8th entry |
| `test/runtests.jl` | test suite | assertion | `####` banner-per-testset; replace ~17 literal `7`s with `length(ProteinCoLoc.theta_prior_bounds())` |
| `spike/npe/train_npe.jl` | trainer | batch | docstring-then-function; `use_gpu=false` on every call; `save_npe` atomic write |
| `docs/amortized.md` | docs | (document) | numbered **Named limits** list voice (see §docs shape) |

---

## Pattern Assignments

### 1. `spike/validation/p11_consts.jl` (config, pre-registration)

**Analog A (primary shape):** `spike/validation/consts.jl` — whole file, 80 lines.

```julia
# spike/validation/consts.jl --- Phase-5 pre-registered constants (SBC/BF/OOD).
#
# THE ANTI-SNOOPING CONTRACT (SBC-03 / D-02). Every M/L/bin-count and every
# pass/fail threshold the Phase-5 validation bundle scores against is LOCKED in
# this ONE file and committed BEFORE any reported run. ...
#
# SEED DISCIPLINE (D-02): VAL_MASTER_SEED is a FRESH Random123 stream, PROVABLY
# DISJOINT from the training/holdout NPE_MASTER_SEED = 0xC0FFEE (train_npe.jl:65).
# ...
# Guarded as ONE block keyed on :SBC_M so re-inclusion under runtests.jl (each of
# harness.jl / sbc.jl / test_sbc.jl guards this include) is a silent no-op --
# a redefinition to the same value would otherwise warn on a `const`.

if !isdefined(@__MODULE__, :SBC_M)
    # --- SBC (SBC-01/02/03) : reported run --------------------------------------
    const SBC_M          = 2000     # pre-registered SBC draws (θ*~π→simulate→infer)
    const SBC_L          = 999      # posterior draws per SBC draw; L+1 = 1000 …
    const SBC_BINS       = 50       # … divisible by SBC_BINS (rank-bin evenness, Pitfall 4)
    ...
    # --- Pre-registered seeds (D-02) --------------------------------------------
    # RESERVED reported stream, consumed ONLY by the Wave-3 run (05-04). DISJOINT
    # from NPE_MASTER_SEED = 0xC0FFEE (train_npe.jl:65) ...
    const VAL_MASTER_SEED = 0x5BC0FFEE   # ≠ 0xC0FFEE (NPE_MASTER_SEED): fresh disjoint stream
    const VAL_FIX_SEED    = 0xF1F7ED
end
```

**Analog B (the executable-rule + derived-seed half):** `test/gate/gate_consts_8_v2.jl:254-343`.
Copy this exact three-part shape — forbidden set → distinct (MASTER, SALT) pair → redraw loop:

```julia
    const AMEND_SALT   = 0xC4CE_B9FE_1A85_EC53
    const AMEND_MASTER = 0x0000_0000_C0DE_2026

    _forbidden_seeds() = (UInt64(0), NPE_MASTER_SEED, VAL_MASTER_SEED, DEFAULT_MASTER_SEED,
                          DEV_SEEDS..., values(PROD_SEED)...)

    function _derive_prod_seed_v2(G::Integer)
        rng = Philox4x(UInt64, (UInt64(AMEND_MASTER) ⊻ AMEND_SALT, UInt64(G)))
        s = rand(rng, UInt64)
        while s in _forbidden_seeds()
            s = rand(rng, UInt64)
        end
        return s
    end
```

**Conventions encoded**
- License block `#= … =#` (spike/`src`) or `#####` MIT banner (`test/gate/`, `scripts/`) — **spike files use the AGPL `#=…=#` block; `test/gate/`, `scripts/`, `test/runtests.jl`, `src/amortized/*` use the `#####` MIT banner.** Match the directory, not the phase.
- First code comment line is always `# <relative/path.jl> --- <one-line purpose>.`
- Whole file wrapped in `if !isdefined(@__MODULE__, :<FIRST_CONST>) … end` so re-`include` under `runtests.jl` is a silent no-op.
- Constants are **column-aligned** with `=` and carry a trailing `#` rationale on the same line.
- Section separators are `# --- <Title> ---------------` padded to ~80 cols (spike) or
  `# ===========` (gate/src).
- Seeds are written in the `0x0000_0000_XXXX_XXXX` underscore-grouped form in `gate_consts_*`,
  bare `0x5BC0FFEE` in spike consts. Either is house; be internally consistent.
- Every seed constant carries an inline `≠ <other seed>` disjointness claim in its comment, and the
  claim is **separately asserted** in a test (never only in prose).

**Naive-pattern-match trap**
`consts.jl` puts *everything* in one `if !isdefined` block committed at once. Phase 11 is **two
tiers**: Tier-1 locked at W0-1 *before the probe*, Tier-2 appended at W2-2 *after the probe*. If you
copy the single-block shape literally you will either (a) have to edit the guarded block to append —
which reads in `git log` as editing a pre-registration — or (b) hit a `const` redefinition warning.
**Use two separate guard blocks keyed on distinct sentinel consts** (e.g. `:P11_DEV_SEED` and
`:SC2_SPEARMAN_FLOOR`), and write a comment at the Tier-2 block saying it is **appended, never
edited**, with the probe artifact path it was derived from. `gate_consts_8_v2.jl:15-26`
("PROVENANCE DISCLOSURE … NOTHING HAS BEEN RUN AGAINST THIS FILE") is the voice for that note.

---

### 2. `spike/test/test_p11_consts.jl` (test, unit)

**Analog:** `spike/test/test_sbc.jl:98-107`.

```julia
    @testset "SC3 (SBC-03) pre-registration + disjoint seed (anti-snooping)" begin
        # The reported stream MUST be disjoint from the training/holdout stream (D-02).
        @test VAL_MASTER_SEED != NPE_MASTER_SEED
        @test VAL_MASTER_SEED == 0x5BC0FFEE     # the pre-registered locked value
        @test VAL_FIX_SEED    != VAL_MASTER_SEED # fixtures never consume the reported stream
        # the pass/fail thresholds are the locked pre-registered values.
        @test SBC_M == 2000
        @test SBC_L == 999
        @test (SBC_L + 1) % SBC_BINS == 0        # rank-bin evenness (Pitfall 4)
    end
```

**Conventions encoded** — the locked value is asserted **as a literal** in the test (`== 0x5BC0FFEE`),
so editing the consts file after the fact breaks the suite loudly. Disjointness is asserted
pairwise, not by a helper. Divisibility/structural invariants (`(L+1) % bins == 0`) live here too.

**Naive-pattern-match trap** — `test_sbc.jl` asserts disjointness against *two* named seeds. Phase 11
must assert `P11_DEV_SEED ∉ _p11_forbidden()` against **six** (`PROD_SEED_V2`, `VAL_MASTER_SEED`,
`NPE_MASTER_SEED`, `VAL_FIX_SEED`, `CORPUS_MASTER_SEED`, `DEV_SEEDS`, `DEFAULT_MASTER_SEED`). Follow
`gate_consts_8_v2.jl`'s `_forbidden_seeds()` set form, and note that `gate_consts_8_v2.jl` **must be
loaded into an isolated `module`** if you need its constants — it defines the same names as
`gate_consts_8.jl` and is guarded on `:SBC_M` (see `scripts/train_grid8_amended_v2.jl:25-29`, the
`module GateV2` idiom, verbatim below).

```julia
# The frozen amended pre-registration, loaded into an ISOLATED module: it defines the same const
# names as `gate_consts_8.jl`, and we want ONLY its imsize mixture here. Nothing is run against it.
module GateV2
    include(joinpath(@__DIR__, "..", "test", "gate", "gate_consts_8_v2.jl"))
end
```

---

### 3. `spike/test/test_stage6_regression.jl` + `spike/test/fixtures/p11_stage6_golden.jld2`

**Analog (test shape):** `spike/test/test_simulator.jl:90-110`.

```julia
# --- Wave-2 unit under test: the forward physics simulator (SIM-01) -----------
# include() AFTER contract.jl so build_mci/summary/induced_mu are already in scope.
include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))

# A literal 7-field θ (D-04). Helper rebuilds it with a swapped ρ_true for sweeps.
const θ_BASE = (ρ_true = 0.7, spillover = 0.1, autofluorescence = 0.05,
                label_efficiency = 0.85, shift_dx = 0.3, shift_dy = -0.4,
                noise = 0.5)
_θ(ρ) = merge(θ_BASE, (ρ_true = ρ,))

@testset "SIM-01 forward pipeline" verbose = true begin

    out = simulate_pair(Random.Xoshiro(2026), θ_BASE; imsize = (256, 256))

    @testset "shape / type / finiteness (returns 2× 256×256 Matrix{Float64})" begin
        @test out isa Vector{Matrix{Float64}}
        @test length(out) == 2
        @test all(c -> size(c) == (256, 256), out)
        @test all(c -> all(isfinite, c), out)     # no NaN/Inf (warp fillvalue trap held)
```

**Analog (fixture write idiom):** `spike/npe/train_npe.jl:160-180` — the only sanctioned artifact write.

```julia
function save_npe(path, result::NamedTuple; master_seed, fold)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp;
        schema_version = NPE_MODEL_SCHEMA,
        estimator = result.estimator,
        ...
        meta = (master_seed = master_seed, fold = fold, ...,
                generated = string(Dates_now())))
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "estimator") "save_npe: integrity check failed, $tmp missing estimator"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end
```

**Conventions encoded**
- Named `const θ_BASE` literal + `merge(θ_BASE, (field = v,))` for variants — never re-typing the tuple.
- Outer `@testset "<ID> <name>" verbose = true begin` with nested `@testset` per property; the
  nested name states the *property*, not the function.
- Fixture RNG is an explicit constructor (`Random.Xoshiro(2026)`), never `seed!` global state.
- Artifact writes always carry a `generated = string(Dates.now(Dates.UTC)) * "Z"` field
  (`train_npe.jl:182-186` names the reason: unlabeled local time was an audit defect, IN-07).

**Naive-pattern-match trap**
The golden fixture is the *only* artifact in this repo whose write must happen **before** the code
change it protects. `save_npe`'s shape is right but its *timing discipline* is not encoded anywhere —
add an explicit header line in the capture script stating the commit sha at capture time
(`PHASE11_BASE_SHA`) and store it **inside** the JLD2, so a fixture regenerated after the ε edit is
self-evidently invalid. Also: `test_simulator.jl` uses `isapprox`-style tolerance checks elsewhere;
SC1d requires **exact `==`** (RESEARCH §A3 proved bit-for-bit equality is achievable) — do not
soften it to `isapprox` when it first fails.

---

### 4. `spike/test/test_p11_smoke.jl` and `spike/test/test_lambda_ablation.jl`

**Analog (smoke):** `spike/test/runtests.jl:43-58` — resolve-risk + CPU-only gate shape.

```julia
@testset "NeuralEstimators CPU smoke" verbose = true begin

    # Runs the full NPE train + sample pipeline; defines mu_hat, theta_true, tol.
    include(joinpath(@__DIR__, "..", "00_smoke.jl"))

    @testset "ENV-02 correctness (recovered posterior mean)" begin
        @test abs(mu_hat - theta_true) < tol
    end

    @testset "D-04 CPU-only (no CUDA dependency, none loaded)" begin
        @test !haskey(Pkg.project().dependencies, "CUDA")
        @test !any(p -> occursin("CUDA", p.name), values(Pkg.dependencies()))
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
```

**Analog (ablation tripwire):** `spike/test/test_sbc.jl:63-68` — the **both-directions** discipline.

```julia
        # (b) a DELIBERATELY biased rank vector (all zeros) drives both p ≈ 0 below the
        #     pre-registered alphas — the test exercises BOTH directions (calibrated vs
        #     miscalibrated), so a stuck "always-pass" p-value would be caught.
        ub = sbc_uniformity(zeros(Int, SBC_FIX_M); L = SBC_FIX_L, bins = SBC_FIX_BINS)
        @test ub.ks_p   < SBC_KS_ALPHA
        @test ub.chi2_p < SBC_CHI2_ALPHA
```

**Conventions encoded** — every spike test file ends with a `@testset "… ran CPU-only (D-10)"`
clause asserting CUDA is absent from `Base.loaded_modules` (`test_sbc.jl:114-116`). Copy it.
Every new dependency-free phase re-asserts the NeuralEstimators v0.2.1 UUID pin in `runtests.jl`
(clauses (d)–(h), `spike/test/runtests.jl:59-113`) — **Phase 11 adds no dependency, so add clause
(i) that asserts exactly that** (no new names, pin re-asserted).

**Naive-pattern-match trap** — `test_lambda_ablation.jl` is *not* a shape test; it is the R4
behavioural tripwire that blocks Wave 4. Write it so a failure reads as "the λ conditioning is
dead", i.e. assert a **material inequality** (`sd(Δρ | λ=3.0) > sd(Δρ | λ=0.25) * <pre-registered
factor>`) against a Tier-2 const, not `!=`. `test_sbc.jl`'s biased-vector clause is the
right *spirit* (prove the statistic can fail), the wrong *statistic*.

---

### 5. `spike/test/test_p11_tost.jl` and `spike/test/test_p11_forced_theta.jl`

**Analog (pure-function test + hand-computed vector):** `test/gate/gate_consts_8_v2.jl:354-371`
is the function to copy into `spike/` (with attribution); its unit test lives in `test/runtests.jl`
and asserts the two in-repo Holm implementations agree.

```julia
    function holm_adjusted(p::AbstractVector{<:Real})
        m   = length(p)
        ord = sortperm(collect(float.(p)))
        adj = Vector{Float64}(undef, m)
        running = 0.0
        for (j, i) in enumerate(ord)
            running = max(running, (m - j + 1) * float(p[i]))
            adj[i]  = min(1.0, running)
        end
        return adj
    end
```

**Analog (field-order / row-alignment invariant):** `test/runtests.jl:534-563`.

```julia
    b = ProteinCoLoc.theta_prior_bounds()
    @test length(b) == 7
    # SINGLE SOURCE OF TRUTH: every bound is the prior object's own support, not a literal.
    @test b[2] == (minimum(ProteinCoLoc.SPILLOVER_PRIOR), maximum(ProteinCoLoc.SPILLOVER_PRIOR))
    ...
    rngp = Random123.Philox4x(UInt64, (UInt64(20260721), UInt64(0)))
    θp   = hcat([collect(values(ProteinCoLoc.sample_prior(rngp))) for _ in 1:200]...)
```

**Conventions encoded** — a copied helper carries an attribution comment naming the source
file:line and the reason it was copied rather than depended on (`spike/Project.toml` lacks the dep).
`collect(values(θ))` is the house θ→vector map, so **NamedTuple field order is a load-bearing
contract** and gets its own test.

**Naive-pattern-match trap (the single most dangerous one in this phase)** — `sbc_holm_pass`
returns `true` when Holm rejects **nothing** (`all(p̃ > fwer)`, gate_consts_8_v2.jl:373-379). The
TOST family needs the **inverted** direction: `all(<=(fwer), holm_adjusted(p_tost))` — pass means
Holm *does* reject the non-equivalence nulls. Copying `holm_pass` verbatim silently inverts the SC2b
verdict. Give the new function a **different name** (`p11_tost_pass`) and hand-compute a fixture
vector in the test (RESEARCH Pitfall 3, §E16).

---

### 6. `spike/simulator/prior.jl` (modified — add `ε`, widen `SHIFT_PRIOR`)

**Analog:** itself, `:48-82`.

```julia
# --- Nuisance prior ranges (Claude's discretion, D-03) ---------------------------
# MIRROR of calibration.jl's `_sample_nuisances` (the calibration measured E[μ|ρ]
# averaging over EXACTLY these ranges; documented in spike/NOTES.md §3). ...
const SPILLOVER_PRIOR        = Uniform(0.0, 0.2)   # directional bleed-through (modest)
const AUTOFLUORESCENCE_PRIOR = Uniform(0.0, 0.1)   # additive background offset
const LABEL_EFFICIENCY_PRIOR = Uniform(0.6, 1.0)   # Bernoulli keep-probability
const SHIFT_PRIOR            = Uniform(-1.0, 1.0)  # sub-pixel registration error (dx, dy)
const NOISE_PRIOR            = Uniform(0.0, 1.0)   # Poisson+Gaussian noise scale

function sample_prior(rng::AbstractRNG)
    μ_star = rand(rng, MU_PRIOR)
    return (
        ρ_true           = ghat(μ_star),
        spillover        = rand(rng, SPILLOVER_PRIOR),
        ...
        shift_dy         = rand(rng, SHIFT_PRIOR),
        noise            = rand(rng, NOISE_PRIOR),
    )
end
```

**Conventions encoded** — one `const <NAME>_PRIOR = <Distribution>(...)` per line, aligned `=`,
trailing physical-meaning comment; `sample_prior` mirrors that order 1:1 and every field is
`rand(rng, <NAME>_PRIOR)` (no inline distribution literals).

**Naive-pattern-match trap** — the *spike* prior and the *src* prior are byte-mirrors of each other
by deliberate policy (`src/amortized/simulator.jl:89-95` says "copied UNCHANGED"). D-11 mirrors
**`ε` only**; D-02's `SHIFT_PRIOR` widening to `Uniform(-3, 3)` is **spike-only**. If you "keep the
mirror" naively you desynchronise `theta_prior_bounds()[5:6]` from the shipped bundle's frozen
`θzt.lo/hi` (RESEARCH Anti-Patterns). Add an explicit comment on the widened spike line saying
*why* it deliberately diverges from src.

---

### 7. `spike/simulator/forward.jl` stage 6 (modified — composed affine)

**Analog:** itself, `:24-45` (header stage list) and `:168-174` (the stage being replaced).

```julia
#   (6) sub-pixel shift ch2 only   warp(ch2, Translation(dy,dx); fillvalue=BG_FLOOR)
...
    # --- (6) sub-pixel registration shift on channel 2 only (fillvalue > 0) ------
    # WR-06: Translation(a, b) shifts the FIRST array axis (rows = vertical = dy) by a
    # and the SECOND (columns = horizontal = dx) by b. Pass (dy, dx) so the named
    # fields map to their conventional physical axes (dx horizontal, dy vertical).
    shifted = warp(ch2, Translation(θ.shift_dy, θ.shift_dx), axes(ch2);
                   method = BSpline(Linear()), fillvalue = BG_FLOOR)
    ch2 = Matrix{Float64}(collect(shifted))
```

And the entry-guard wall it must extend (`:119-134`):

```julia
    all(isfinite, (θ.ρ_true, θ.spillover, θ.autofluorescence, θ.label_efficiency,
                   θ.shift_dx, θ.shift_dy, θ.noise)) ||
        throw(ArgumentError("all θ fields must be finite, got $θ"))
    (-1.0 ≤ θ.ρ_true ≤ 1.0) ||
        throw(ArgumentError("ρ_true must be in [-1, 1], got $(θ.ρ_true)"))
    # WR-07: reject physically out-of-range nuisances rather than silently clamping
```

**Conventions encoded**
- Every stage carries a `# --- (n) <name> ----` marker **and** a matching one-liner in the file
  header's numbered stage list. Both must be updated together.
- Comments carry a **ticket/finding prefix** (`WR-06`, `WR-07`, `CR-02`, `IN-03`, `D-15`,
  `Pitfall 4`) naming the review item that forced the line. A new mechanism comment without such a
  tag reads as untraceable in this codebase — tag the re-derived axis-order comment (D-10 says
  "the `WR-06` axis-order comment must be **re-derived**", so keep the `WR-06` tag and extend it).
- Guards are `<predicate> ||\n    throw(ArgumentError("<field> must be …, got $(θ.<field>)"))`,
  one per field, in θ field order.
- Fixed non-θ constants live at the top of the file as `const` with a rationale (`σ_psf`, `BG_FLOOR`).

**Naive-pattern-match trap** — the existing guard block uses a **hard-coded 7-tuple** in
`all(isfinite, (...))`. Adding `ε` requires touching that tuple *and* the header stage list *and*
the `θ is the 7-field NamedTuple (D-04)` line at `:44-45`. Three places, none of which the compiler
checks. Also: `simulate_pair` must keep accepting a 7-field θ (SC1e backward compat) — read `ε`
defensively (`get(θ, :chromatic_eps, 0.0)`-style), which is *not* the existing style but is required.

---

### 8. `src/amortized/simulator.jl` (modified — mirror `ε`)

**Analog:** itself, `:19-40` (banner) and `:122-148` (`theta_prior_bounds`).

```julia
#############################################################################################
# src/amortized/simulator.jl --- the Phase-2 forward physics simulator, PROMOTED into src/.
#
# Promoted VERBATIM (prior ranges + physics unchanged) from the proven spike:
#   spike/simulator/prior.jl   -> sample_prior (π(θ), the 7-field D-04 NamedTuple)
...
# HARD CONSTRAINT (CLAUDE.md prior-consistency): the simulator prior π(θ) is kept BYTE-CONSISTENT
# with the existing Turing @model prior ranges ...
#############################################################################################
...
**SINGLE SOURCE OF TRUTH.** Every bound is DERIVED from the prior objects defined above ...
Nothing here is a duplicated literal: change a prior above and this box follows.
"""
theta_prior_bounds() = (
    (first(GHAT_RHO_KNOTS), last(GHAT_RHO_KNOTS)),                       # ρ_true = ghat(μ*)
    (minimum(SPILLOVER_PRIOR),        maximum(SPILLOVER_PRIOR)),
    ...
    (minimum(NOISE_PRIOR),            maximum(NOISE_PRIOR)),
)
```

**Conventions encoded** — `src/amortized/*` uses the `#####`-bordered MIT banner + a
`# <path> --- <purpose>.` line + ALL-CAPS labelled paragraphs (`HARD CONSTRAINT:`,
`REPRODUCIBILITY (D-11/D-14):`, `API GOTCHA`, `WHY PROMOTED HERE`). Docstrings use `"""` with a
signature line, then prose, then a bolded invariant. Bounds are **derived, never literal** — the
new `ε` row must be `(minimum(CHROMATIC_PRIOR), maximum(CHROMATIC_PRIOR))`.

**Naive-pattern-match trap** — the file's own header says "7-field D-04 NamedTuple" in three places
and `theta_prior_bounds`'s docstring enumerates the field order in prose. All of those are now
stale. More importantly: **this file's bytes are hashed** (`datagen.jl:208-213`), so a
comment-only touch is still a provenance event. Batch every edit into the single D-11/D-12 commit
(W1-1) — a follow-up "fix the comment" commit re-digests the cache a second time.

---

### 9. `src/amortized/datagen.jl` and `src/amortized/ood.jl` (modified)

**Analog:** `datagen.jl:199-213`.

```julia
# ============================ content-hash version guard ====================================
# A hash over the data-DEFINING src source bytes PLUS a canonical serialization of the config.
# The digest names the cache directory, so ANY change to the summary source OR to a hashed config
# field ... yields a different digest and auto-separates stale/cross-grid data.

# FIXED-order data-defining src files (@__DIR__ = src/amortized/). ...
const DATAGEN_HASH_SRC_FILES = String[
    joinpath(@__DIR__, "..", "colocalization.jl"),   # patch() / correlation() (defines the summary)
    joinpath(@__DIR__, "..", "LoadImages.jl"),        # MultiChannelImage container
    joinpath(@__DIR__, "summary.jl"),                 # patch_summary(mci,G) / encode_d01
    joinpath(@__DIR__, "simulator.jl"),               # sample_prior / simulate_pair / build_mci (07-05)
]
```

**Analog:** `ood.jl:120-132`.

```julia
# Reconstruct a simulate_pair-valid θ NamedTuple from a physical-θ posterior-mean vector, first
# mapping any non-finite component to an in-range fallback, then clamping into the prior-valid
# ranges (mirrors simulate_pair's entry guard: ρ∈[-1,1], spillover∈[0,1], ...).
_theta_tuple(v) = (
    ρ_true           = clamp(_finite_or(v[1], 0.0), -1.0, 1.0),
    spillover        = clamp(_finite_or(v[2], 0.0),  0.0, 1.0),
    autofluorescence = max(_finite_or(v[3], 0.0), 0.0),
    label_efficiency = clamp(_finite_or(v[4], 1.0),  0.0, 1.0),
    shift_dx         = _finite_or(v[5], 0.0),
    shift_dy         = _finite_or(v[6], 0.0),
    noise            = max(_finite_or(v[7], 0.0), 0.0),
)
```

**Conventions encoded** — positional θ indexing (`v[1]`…`v[7]`) mirroring `sample_prior` field
order; each component gets its own clamp matching `simulate_pair`'s entry guard, with the mirroring
stated in the comment.

**Naive-pattern-match trap** — `_theta_tuple` indexes positionally into a **7-row** draw matrix.
The 8th row must be appended at the **end** (matching `sample_prior`'s new field order) or every
existing index shifts. Do **not** insert `ε` next to `shift_dx`/`shift_dy` "because it's also
geometric" — RESEARCH §A4 enumerates the ripple and the append-at-end choice is what keeps
`v[1]..v[7]` valid. Likewise **do not reorder `DATAGEN_HASH_SRC_FILES`** (fixed order is the hash
input) and do not delete orphaned caches (they are the byte record D-12 preserves by reference).

---

### 10. `test/runtests.jl` (modified)

**Analog:** `test/runtests.jl:23-32` (banner-per-testset) and `:534-543` (the literal-7 site).

```julia
##########################################################################################
### CO-RESOLUTION HARD GATE (Phase 7 / Finding 1)  — must run FIRST
###
### Adding NeuralEstimators/Flux to the root package alongside the legacy plot/Turing stack
### historically capped NeuralEstimators below 0.2.1 and silently DOWNGRADED it to 0.1.4 ...
##########################################################################################
@testset "co-resolution gate (Finding 1)" begin
```

**Conventions encoded** — each top-level `@testset` is preceded by a `####`-bordered `###`-prefixed
rationale banner naming the finding/decision ID it defends. Test names are lowercase prose with the
governing ID in parentheses.

**Naive-pattern-match trap** — RESEARCH counts **~17 literal `7`s** to change. Replace with
`length(ProteinCoLoc.theta_prior_bounds())` (the declared single source of truth) rather than
retyping `8`; otherwise Phase 12/13 repeats this chore. But note `:556` (`wild = 50.0 .* randn(7, 500)`)
and `:562` (`repeat(collect(-6.0:0.5:6.0)', 7, 1)`) are **matrix constructions**, not assertions —
they must become the derived length too, or `reconstruct` throws a `DimensionMismatch`.

---

### 11. `BoundedThetaTransform` port into `spike/` (if D-03 needs it)

**Analog:** `src/amortized/architecture.jl:43-104` — copy the comment block wholesale, it is the
rationale.

```julia
# =============================================================================================
# BOUNDED θ-SPACE (07-CALIBRATION-FINDINGS F2)
# =============================================================================================
#
# THE DEFECT: π(θ) is a hard-truncated BOX ... but `NeuralEstimators.NormalisingFlow` models θ
# in an UNBOUNDED space after a plain z-score. ...
#
# THE FIX: interpose a per-parameter LOGIT bijection from the prior box onto ℝ BEFORE the
# z-score, and invert it on read-out. ...
#
#   forward :  u = (θ − lo)/(hi − lo) ∈ [0,1]  →  y = logit(clamp(u, ε, 1−ε))  →  z-score
#   inverse :  un-z-score  →  θ = lo + (hi − lo)·logistic(y)                    ∈ [lo, hi]
#
# BOUNDS PROVENANCE: `lo`/`hi` come from `theta_prior_bounds()` (simulator.jl) ...
# No bound is ever hardcoded here.

struct BoundedThetaTransform{T<:Real,Z}
    lo::Vector{T}
    hi::Vector{T}
    zt::Z
end
```

**Naive-pattern-match trap** — the ε-clamp comment says the atoms come from `ghat`'s clamped
endpoints. Under Phase 11 the **`ε` (chromatic) column has no atoms**, so do not copy the
atom-handling prose onto it. Also beware the name collision: `THETA_LOGIT_EPS`'s `ε` and D-09's
chromatic `ε` are different things — RESEARCH §D uses `chromatic_eps` as the field name; keep that
and never write a bare `ε` in the new code.

---

## Experiment-script shape (the house form for a long-running reproducible run)

**Closest analog: `spike/validation/run_sbc.jl` (177 lines).** Its skeleton, in order:

1. AGPL `#=…=#` license block.
2. `# spike/validation/run_sbc.jl --- Phase-5 REPORTED SBC run (SBC-01..04).` then ALL-CAPS
   labelled paragraphs: what this proves, the anti-snooping contract, the caption, and finally
   the **literal run line**:
   ```julia
   # CPU-only (D-10), read-only src/ (via the harness). Run:
   #     julia --project=spike spike/validation/run_sbc.jl
   ```
   **There is no `Pkg.activate` line and no arg parsing anywhere in `spike/`** — the project is
   selected by `--project=spike` on the command line, and knobs are `const`s at the top of the
   script. `-t auto` is added on the command line for threaded runs
   (`scripts/train_grid8_amended_v2.jl:21`: `julia --project -t auto scripts/...`). Do not add
   `ArgParse`; do not call `Pkg.activate`.
3. `using Test / JLD2 / Statistics / Dates`, then **guarded includes** of the units under test:
   ```julia
   isdefined(@__MODULE__, :sbc_ranks)           || include(joinpath(@__DIR__, "sbc.jl"))
   isdefined(@__MODULE__, :plot_rank_histogram) || include(joinpath(@__DIR__, "figures.jl"))
   ```
4. `const <X>_REPORT_PATH = joinpath(@__DIR__, "<x>_report.jld2")` — **the artifact sits next to
   the script**, named `<experiment>_report.jld2`.
5. Resolution knobs as file-local `const`s, explicitly labelled as *not* pre-registered thresholds
   (`run_ood.jl:70-76`: "Monte-Carlo resolution knobs … NOT anti-snooping thresholds").
6. A private `_save_<x>_report(path; kwargs...)` using the atomic idiom:
   ```julia
   function _save_sbc_report(path; kwargs...)
       mkpath(dirname(path))
       tmp = path * ".tmp"
       jldsave(tmp; kwargs...)
       JLD2.jldopen(tmp, "r") do f
           @assert haskey(f, "rank_table") "_save_sbc_report: integrity check failed ($tmp)"
       end
       mv(tmp, path; force = true)
       return path
   end
   ```
7. `function main()` opening with a **banner block and seed assertions**:
   ```julia
       println("="^78)
       println("Phase-5 REPORTED SBC run (SBC-01..04) — LOCKED consts, RESERVED VAL stream")
       println("  M=$SBC_M  L=$SBC_L  bins=$SBC_BINS  imsize=$SBC_IMSIZE")
       println("  stream=VAL_MASTER_SEED=$(repr(VAL_MASTER_SEED))  (≠ NPE_MASTER_SEED=$(repr(NPE_MASTER_SEED)))")
       println("="^78)
       @assert VAL_MASTER_SEED != NPE_MASTER_SEED "reported stream must be disjoint from training stream (D-02)"
   ```
8. **Progress logging** = `println("[1/2] running the density-channel ROC over the full grid …")`
   plus `t0 = time()` / `elapsed = time() - t0` and a
   `println("… in $(round(elapsed/60; digits=2)) min")`. (`@info` with keyword payloads is the
   `scripts/`-under-`src/` style — `scripts/train_grid8_amended_v2.jl:43,54,59`; `println` is the
   `spike/` style. Match the directory.)
9. **Persist BEFORE the verdict** — with the comment stating why:
   ```julia
       # --- Persist the reported artifact BEFORE the gate (so it survives an honest fail) ---
   ```
   The saved payload includes the **raw data, the labels, every locked const, the seed, the caption,
   `elapsed_min`, and `generated`**. Copy that field list shape.
10. Figures written after persistence, explicitly labelled "script artifacts — never gate assertions".
11. A printed `rpad`-aligned PASS/FAIL table **before** any throw.
12. Bare `main()` as the last line (no `if abspath(PROGRAM_FILE) == @__FILE__` guard in this repo).

**Result-artifact shape:** a single `<name>_report.jld2` written atomically, keyed by plain strings,
containing (a) the raw per-draw table, (b) every derived statistic as a parallel `Vector`, (c) a
copy of every pre-registered constant that governed the run, (d) the seed as `UInt64`, (e) prose
caption, (f) `elapsed_min`, (g) UTC `generated` stamp.

**The Phase-11 divergence you must make deliberately.** `run_sbc.jl`/`run_ood.jl` end with a
`@testset` that **exits nonzero** — they are gates. Phase-11 `run_p11_breakdown.jl`,
`run_p11_attenuation.jl` and `run_p11_realimage.jl` are **reported, not gated** (D-14: "the
breakdown point beyond it is REPORTED, not gated"), and RESEARCH §Sampling Rate says experiments
must not be wired into CI. For those, keep steps 1–11 and **replace step 12's `@testset` with a
printed summary table only** — and say so in the header banner, in the voice of
`spike/npe/scaling.jl:26-30`:

```julia
# This is CHARACTERIZATION (reported), not a pass/fail
# gate -- the honest, defensible replacement for a single cherry-picked speedup point.
```

`run_p11_ladder.jl` **is** gated (SC2a/SC2b have pre-registered criteria) — it keeps the
`@testset` ending. `run_p11_probe.jl` is neither: it *produces* thresholds, so it must end by
printing the derived values in copy-paste-ready `const` form for the Tier-2 append, and must
**never** write to `p11_consts.jl` itself.

---

## Report-artifact shape (`.planning/phases/11-…/`)

**Closest analog: `.planning/phases/07-productionization-conditional-on-go/gate-8x8-amended.md`
(137 lines).** Full outline, verbatim heading text:

```
# Grid 8×8 — Amended Ship-Gate, Confirmatory Run

**Run:** 2026-07-21 17:11–17:42 (30.9 min), CPU-only, single binding invocation
**Net:** `artifacts/amended_v2/grid_8/` (bounded logit-θ space, trained on …, 50k pairs)
**Pre-registration:** `test/gate/gate_consts_8_v2.jl`, frozen at commit `c69d381` —
                     **unchanged between freezing and this run** (verified by `git diff`)
**Seed:** `PROD_SEED_V2` → 1135605683775656488 (fresh, asserted disjoint from all prior seeds)
**Protocol:** amendment §6 — ONE run, no retry, no re-seed, no threshold change. Honoured.

---

## Overall verdict: **FAIL**                     ← verdict FIRST, before any interpretation
<3-row arm summary table>
The gate FAILED. This is reported as the primary result; interpretation is separated below.

---

## 1. SBC — FAIL under **both** rule sets       ← one ## per arm, verdict in the heading
<per-parameter table: value | threshold-derived stat | shrinkage | vacuous>
### 1a. Five of eight parameters are VACUOUS    ### sub-sections carry the caveats
### 1b. ρ_true calibration BROKE under the realistic image-size regime
### 1c. The bounded-θ fix did help, but not enough
### 1d. Realised image-size mixture

## 2. BF — FAIL, and the correlation verdict is `:invalid`, not merely failing
## 3. OOD — PASS
## 4. Limitation the amendment conceded, restated
## 5. What this run establishes                 ← numbered, quotable one-claim-per-item list
## 6. Open items (none actioned in this run)
```

**How pre-registered thresholds vs measured values are presented** — always in the *same table row*,
measured first, threshold in a parenthetical or its own column:

> `| one-sided 95% lower bound | 0.9255 (needs ≥ 0.96394) |`
> `| n_min required | 69 |` … `| n (finite pairs) | **58** |`

**How verdicts are worded** — bolded uppercase `**FAIL**` / `**PASS**` in the heading itself; the
*reason class* is distinguished from the *statistic*:

> "This is a **process failure needing investigation**, not a clean statistical fail, and it must
> not be reported as one."

**How named limits are worded** — a limit is stated as a *positive property of the evidence*, then
its scope, then what would close it:

> "**5/8 parameters are vacuous.** The correlation-only summary constrains ρ_true, Δρ and
> (partially) label_efficiency; it learns essentially nothing about spillover, autofluorescence,
> shift_dx, shift_dy and noise. **Any future SBC "pass" on those columns must be read as vacuous.**"

**Companion analog for the "findings" register** (if Phase 11's report separates measurement from
inference): `07-CALIBRATION-FINDINGS.md` uses explicit epistemic labels — `**MEASURED:**`,
`**INFERRED:**`, `**NOT CONCLUSIVE —**`, `**What would close it:**` — plus a `## Severity ordering`
section up front and a `## What was deliberately NOT done` section at the end. Phase 11's report
should carry **both** the gate-memo skeleton (verdict-first, per-SC sections) and these epistemic
labels, because D-05/D-14 mix gated claims (SC2, SC3-in-prior) with reported-only ones (breakdown
point, attenuation, real-image transfer).

**Phase-11-specific sections the analog does not have (add them):**
- a `## The SC2/SC3 tension` section in **plain language** (D-13 + `<specifics>` explicitly ask that
  it not be buried) — that widening the prior to demonstrate honest uncertainty is precisely what
  makes mis-registration undetectable by the density channel;
- a vacuous-column caption on every `shift_*`/`ε` table, per-rung `prior_sd` (Pitfall 5);
- the five §F20c decoupling/provenance command outputs, stated as run and empty (D-17 requires the
  report to *say* the manuscript pipelines are unmodified);
- the pinned pre-`ε` sha, quoted (D-12).

---

## `docs/amortized.md` named-limit shape

**Limit #4, VERBATIM** (`docs/amortized.md:128-135`) — write limit #8 in exactly this voice:

```markdown
4. **Nuisance-parameter marginals carry a ~0.08-SD residual drift** not certifiable as negligible at
   M=2000. Under a pre-registered δ=0.10-SD equivalence (TOST) test, 4/6 nuisances pass; 2/6
   (autofluorescence, label_efficiency) fail — not because the drift is large (point drift ~0.08 SD
   < 0.10) but because the 90% CI edges just past 0.10. **No model lever removes it** (capacity
   redistributes it, spike 012; μ-truncation fixes only the atoms at a real cost to |ρ|≥0.95
   inference and ADVI comparability, spike 013). This is a **marginal-reproduction limit on
   parameters the fixed patch-correlation summary cannot constrain** (5/8 columns are vacuous —
   posterior ≈ prior), **not** a failure of the coloc inference. It is a named calibration limit.
```

**The voice, decomposed** (all five moves appear in #4 and in #2, #3, #6, #7):
1. **Bold lead clause naming the defect quantitatively** — "**Nuisance-parameter marginals carry a
   ~0.08-SD residual drift**".
2. **The measurement conditions** — the test, its pre-registered threshold, the count that passed.
3. **The mechanism / why it is not fixable**, with parenthetical evidence citations
   (`spike 012`, `spike 013`, `07-GO-NO-GO-UPDATE.md §4`, `gate-8x8-amended.md`).
4. **A bolded scope sentence that says what the limit is NOT** — "**not** a failure of the coloc
   inference".
5. **A closing classification sentence** — "It is a named calibration limit."

Second person is never used; the subject is always the artifact or the evidence. Em-dashes carry
the concessions. Emphasis is `**bold**`, never italics-for-emphasis.

**Where limit #8 goes** — the list is a plain numbered `1.`–`7.` markdown list under
`## Named limits (non-negotiable honesty — carry these into the manuscript)` (line 106), which is
introduced by:

> The v2.0 GO is **GO with named limits** (`07-GO-NO-GO-UPDATE.md §5`). These are properties of the
> shipped 8×8 estimator and its calibration evidence; they must be stated with any published result.

Append `8.` **after item 7 and before** the closing `**Findings and scripts:**` line (`:152-154`),
and extend that closing line's path list with the Phase-11 report. Because the intro scopes the list
to "properties of the shipped 8×8 estimator", limit #8 (D-12 provenance-by-reference) fits
naturally; it must quote the pre-`ε` sha explicitly.

**Where the D-16 interpretation section goes** — current `##` heading order in `docs/amortized.md`:

| Line | Heading |
|---|---|
| 1 | `# Amortized colocalization inference (v2.0 "AmortizedColoc")` |
| 15 | `## Public entry point` |
| 47 | `## The grid registry (num_patches as key, PROD-02)` |
| 71 | `### Content-hashed, lazy artifact loading` |
| 87 | `## CPU-reproducible by default (D-06)` |
| 95 | `## Windowed sub-tile local map (coarse)` |
| **←** | **`## Reading posterior width under registration uncertainty (Phase 11, research net)`** |
| 106 | `## Named limits (non-negotiable honesty …)` |

**Insert the D-16 section as a new `##` between line 105 and line 106** — i.e. after
"Windowed sub-tile local map (coarse)" and immediately before "Named limits". Rationale from the
existing structure: everything before line 106 is *how to read the tool*; line 106 onward is *what
not to claim*. D-16's section is an interpretation guide, so it belongs in the first group, and
sitting last in that group puts it adjacent to the limit (#8) it cross-references.

The nearest tonal analog for that section is `## Windowed sub-tile local map (coarse)`
(`:95-104`), which does exactly the job D-16 needs — describe a capability, then immediately
and bluntly bound what it is not:

> "It is a **coarse windowed readout, not a calibrated per-region posterior**: each tile is an
> independent global-ρ read on a smaller window, with no spatial prior, no borrowing of strength
> between neighbours, and no per-region uncertainty. … The calibrated per-region map (GP/CAR
> lattice, `SpatialColocResult`) is Phase 12 and is **not** implemented here."

D-16's section must do the same, with the additional burden of naming **which net every number came
from** — the shipped `amended_v2/grid_8` bundle vs the Phase-11 research net (`d_in = 129, D = 8`,
widened `SHIFT_PRIOR`, fresh DEV seed), which is *not* the shipped one and is *not* downloadable.

---

## Shared Patterns

### S1. License / banner header (applies to every new `.jl`)
**Source:** `spike/validation/consts.jl:1-19` (AGPL `#=…=#`, all `spike/` files) vs
`test/gate/gate_consts_8_v2.jl:1-5` / `test/runtests.jl:1-5` / `scripts/*.jl:1-5` /
`src/amortized/*.jl` (`#####` MIT banner).
**Apply to:** all new files. **Match the directory, not the phase.** New `spike/` files → AGPL
block. (Note the inconsistency is pre-existing and intentional-by-directory; do not "fix" it.)

### S2. Guarded include / guarded const (applies to every new `spike/` file)
**Source:** `spike/validation/harness.jl:48-58`
```julia
# --- ORDER MATTERS: consts first, then the trainer/inference surface (load_npe +
#     posterior_for + rho_draws), then the simulator halves, contract, encode, and
#     the seeding primitives whose salt idiom val_rng mirrors. Guarded for idempotency.
isdefined(@__MODULE__, :VAL_MASTER_SEED) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :load_npe)        || include(joinpath(@__DIR__, "..", "npe", "train_npe.jl"))
isdefined(@__MODULE__, :posterior_for)   || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))
isdefined(@__MODULE__, :sample_prior)    || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
```
**Apply to:** every new `spike/validation/*` and `spike/test/*` file. Includes are
`joinpath(@__DIR__, ...)`-relative, never bare strings; order is documented in a preceding comment
because it is load-bearing.

### S3. Keyed-RNG derivation (never `Random.seed!` in a loop)
**Source:** `spike/data/seeding.jl:75-96` + `spike/validation/harness.jl:93-94`
```julia
val_rng(master_seed = VAL_MASTER_SEED) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ VAL_SALT, UInt64(0)))
```
**Apply to:** every Phase-11 script and harness extension. A new stream = a new **salt**, with a
comment naming the other salts it is distinct from (`HOLDOUT_SALT`, `FOLD_SALT`, `VAL_SALT`,
`PROD_SALT`, `AMEND_SALT`).

### S4. Atomic artifact write
**Source:** `spike/npe/train_npe.jl:160-180`, `spike/validation/run_sbc.jl:56-65`.
**Apply to:** the golden fixture and all five `run_p11_*.jl` report artifacts. `mkpath` → `.tmp` →
`jldsave` → reopen + `@assert haskey(...)` → `mv(...; force = true)` → `return path`.

### S5. Read-only-`src/` decoupling declaration
**Source:** `spike/npe/train_npe.jl:48-49`
```julia
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; touches no src/. Guarded
# includes keep the file loadable standalone AND inside runtests.jl.
```
**Apply to:** every new `spike/` file — **but Phase 11 must reword it.** D-11 deliberately edits
`src/`, so the honest line is "spike-local; reaches `src/` only read-only through
`contract.jl`" (the `harness.jl:38-39` wording), and the `src/`-touching commit says so explicitly.

### S6. Never hand-roll (from RESEARCH §Don't Hand-Roll — reproduced here only as call sites)
`interval_width` (`spike/npe/infer.jl:148-151`), `sbc_holm_adjusted` (`test/gate/sbc.jl:352-362`),
`sbc_shrinkage` (`test/gate/sbc.jl:243-250`), `roc_auc`/`id_threshold`/`youden_j`
(`src/amortized/ood.jl:268-305`), `SBC_IMSIZE_SET`/`WEIGHTS` read (never retyped) from
`test/gate/gate_consts_8_v2.jl:152-153`.

---

## No Analog Found

| File / concern | Role | Data flow | Reason |
|---|---|---|---|
| `spike/test/fixtures/` **directory** | fixture store | file-I/O | No `fixtures/` directory exists anywhere in the repo. Existing golden data lives alongside its producer (`spike/npe/trained_npe.jld2`, `spike/validation/*_report.jld2`, `spike/baseline/*.jld2`). The directory itself is new; the **write idiom** (S4) and the naming convention (`<subject>_<kind>.jld2`) are established. Either follow RESEARCH's `spike/test/fixtures/` path or place it as `spike/test/p11_stage6_golden.jld2` — the repo precedent favours "next to its consumer". |
| Two-tier (append-only) pre-registration | config | — | Every existing pre-registration is a **single frozen block**. The Tier-1/Tier-2 append pattern is genuinely new; use two guard blocks + the `gate_consts_8_v2.jl:15-26` provenance-disclosure voice (see §1 trap). |
| λ-conditioned (129-row) input assembly | data pipeline | transform | No existing code appends a non-summary row after standardization; `Loader._row_partition` hard-requires exactly 128 rows (`spike/data/loader.jl:126`). This is net-new spike code — follow RESEARCH §D12/§D13, not an analog. |

---

## Conventions

Automated derivation **skipped** — `gsd-tools verify conventions --derive --scope spike` returned
`{"skipped": true, "reason": "no-readable-files"}` (the deriver does not read `.jl`). The table
below is derived by hand from the files read for this map (n = 14 Julia files across `spike/`,
`src/amortized/`, `test/`, `test/gate/`, `scripts/`).

| Axis | Dominant | Share | Entropy | Status |
|---|---|---|---|---|
| File-name casing | `snake_case.jl` (`train_npe.jl`, `gate_consts_8_v2.jl`, `run_sbc.jl`) | ~100% (exceptions: `LoadImages.jl`, the v1 legacy file) | very low | **named contract** |
| Identifier casing | `snake_case` functions/consts-lowercase; `SCREAMING_SNAKE` consts; `CamelCase` types; leading `_` for file-private helpers (`_save_sbc_report`, `_theta_tuple`, `_derive_prod_seed_v2`) | ~95% | low | **named contract** |
| Export style | Flat top-level functions, **no `module` wrapper** in `spike/` (`harness.jl:41` states this explicitly); `src/` files are `include`d into the `ProteinCoLoc` module with no per-file `module` | ~100% within each half | low | **named contract** |
| Import style | `spike/` uses `using <Pkg>` + guarded `include(joinpath(@__DIR__, …))`; `src/amortized/` uses **explicit `import Pkg: name1, name2`** (`simulator.jl:42-48`, `architecture.jl:39-41`) | ~50 / 50 **by directory** | high repo-wide, ~0 per-directory | **contested hotspot** |
| License header | `#=…=#` AGPL block in `spike/`; `#####` MIT banner in `src/amortized/`, `test/`, `test/gate/`, `scripts/` | ~50 / 50 by directory | high repo-wide, ~0 per-directory | **contested hotspot** |
| Progress logging | `println("…")` + `"="^78` banners in `spike/`; `@info "msg" key = val` in `scripts/` | ~50 / 50 by directory | high repo-wide, ~0 per-directory | **contested hotspot** |

**Contested hotspots (author's choice).** Three axes are contested repo-wide and internally
consistent **per directory** — exactly the prototype intentional-contested split (the CJS↔ESM dual
resolver pattern: `bin/lib/**` CJS vs `sdk/src/**` ESM, each half self-consistent, contested only
when pooled). Here the split is `spike/**` (exploratory, AGPL header, `using` + guarded `include`,
`println` banners) vs `src/**` + `test/**` + `scripts/**` (productionized, MIT banner, explicit
`import X: y`, `@info` with keyword payloads). **Match the directory's local style; do not
normalise across the boundary.** Phase 11 straddles it (it edits both `spike/simulator/*` and
`src/amortized/*` in one commit, W1-1) — so the *same* commit will legitimately contain both
styles. That is correct, not an inconsistency.

Two further named contracts worth stating because Phase 11 can violate them silently:
- **Every threshold constant carries an inline rationale comment and a governing decision ID**
  (`D-02`, `SBC-03`, `F5`, `WR-06`, `Pitfall 4`). A bare `const X = 0.9` is off-style.
- **Every derived quantity is derived, never retyped** — `theta_prior_bounds()` from the prior
  objects, `SBC_IMSIZE_SET` read from the frozen consts file, `length(theta_prior_bounds())` instead
  of a literal arity. This is the convention Phase 11 most needs to honour (the ~17 literal `7`s).

---

## Metadata

**Analog search scope:** `spike/` (validation, test, simulator, npe, data), `src/amortized/`,
`test/`, `test/gate/`, `scripts/`, `docs/`, `.planning/phases/07-*/`
**Files read in full or in targeted ranges:** 16
**Pattern extraction date:** 2026-07-25
**Upstream:** `11-CONTEXT.md` (17 locked decisions), `11-RESEARCH.md` (2373 lines — mechanism,
not duplicated here)
