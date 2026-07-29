# Phase 12: Spatial Colocalization Map (GP/CAR) — Pattern Map

**Mapped:** 2026-07-27
**Files analyzed:** 24 (22 new, 2 modified)
**Analogs found:** 22 / 24 (exact 15, role-match 7, none 2)

> Consumed by `gsd-planner`. Every excerpt below carries a file path and line range that was
> **read this session**. Where an analog does not exist, the file is listed under
> "No Analog Found" and the planner should fall back to `12-RESEARCH.md` § Architecture Patterns.

---

## 0. Three constraints that override any analog

1. **`spike/Project.toml` / `spike/Manifest.toml` must stay byte-unchanged.** Verified this
   session: `spike/Project.toml` has exactly 16 deps and contains **no** `LinearAlgebra`,
   `SparseArrays`, `Statistics` or `Dates`. Yet `spike/validation/run_p11_recovery.jl:50-53`
   already does `using JLD2 / Dates / Statistics / LinearAlgebra` and runs green — Julia's default
   `LOAD_PATH` carries `@stdlib`, so **stdlib `using` needs no `Project.toml` entry**. That is the
   whole mechanism for the research's zero-dependency recommendation; `SparseArrays` rides it
   identically. **Any analog that would imply a `Pkg.add` is disqualified** — and
   `spike/test/runtests.jl:123-126` asserts the *exact dependency name set*, so an addition fails
   the suite by name.
2. **Do not edit `src/`.** `src/results.jl:180-191`, `src/amortized/summary.jl`,
   `src/amortized/architecture.jl`, `src/amortized/local_map.jl` are **read-only references**. The
   executable Phase-12 result type goes in `spike/`, following `spike/p13/result.jl` (§ Pattern
   Assignments → result type).
3. **Include-ordering trap (silent-pass).** `spike/test/runtests.jl:208` is
   `include(joinpath(@__DIR__, "test_p13_correction.jl"))`, whose outer `@testset` **throws** on two
   committed measured misses. Lines 191-199 of that file record the consequence verbatim: *"A thrown
   testset aborts the remaining includes, so any sibling placed after it would silently never run."*
   **Every Phase-12 include must be inserted between line 207 (`test_p13_calibration.jl`) and line
   208.**

---

## 1. File Classification

| New / Modified File | Role | Data Flow | Closest Analog | Match |
|---|---|---|---|---|
| `spike/validation/p12_consts.jl` | config (pre-registration) | declarative / seed-derivation | `spike/validation/p11_consts.jl` | exact |
| `spike/test/test_p12_consts.jl` | test (literal gate) | assertion | `spike/test/test_p11_consts.jl` | exact |
| `spike/simulator/p12_lattice.jl` | utility (numeric kernel) | transform | `spike/simulator/ghat.jl` + `spike/simulator/prior.jl:39-58` | role-match |
| `spike/test/test_p12_lattice.jl` | test (unit) | assertion | `spike/test/test_simulator.jl` | role-match |
| `spike/simulator/p12_prior.jl` | model (prior π(θ)) | transform / sampling | `spike/simulator/prior.jl` | exact |
| `spike/test/test_p12_prior.jl` | test (unit + golden regression) | assertion | `spike/test/test_stage6_regression.jl` | exact |
| `spike/test/capture_p12_golden.jl` | script (one-shot capture) | file-I/O | `spike/test/capture_p11_golden.jl` | exact |
| `spike/test/fixtures/p12_stage1_golden.jld2` | fixture | file-I/O | `spike/test/fixtures/p11_stage6_golden.jld2` | exact |
| `spike/simulator/forward.jl` **(MODIFIED, stage 1)** | model (forward sim) | transform | its own stage 6 (`forward.jl:192-218`, the D-09/D-10 edit) | exact |
| `spike/npe/p12_architecture.jl` | model (net builder) | transform | `spike/npe/p11_architecture.jl` + `spike/npe/architecture.jl:68-115` | exact |
| `spike/test/test_p12_architecture.jl` | test (unit + smoke) | assertion | `spike/test/test_p11_architecture.jl` | exact |
| `spike/test/test_p12_decoupling.jl` | test (integration) | process / git | `spike/test/test_p13_result.jl:190-213` + `runtests.jl:114-151` | exact |
| `spike/data/p12_generate.jl` | service (datagen) | batch | `spike/data/p11_generate.jl` | exact |
| `spike/npe/train_p12_npe.jl` | service (trainer) | batch | `spike/npe/train_p11_npe.jl` | exact |
| `spike/validation/run_p12_sim02.jl` | runner (reported) | batch → `.jld2` | `spike/validation/run_p11_recovery.jl` (shape) + `spike/simulator/prior.jl:42-46` (bar) | exact |
| `spike/validation/run_p12_stage1_ridge.jl` | runner (reported, D-12 Stage 1) | batch → `.jld2` | `spike/validation/run_p11_recovery.jl` | exact |
| `spike/validation/run_p12_ell_ridge.jl` | runner (reported) | batch → `.jld2` | `spike/validation/run_p11_recovery.jl` | exact |
| `spike/validation/run_p12_eps_ridge.jl` | runner (reported, read-only on P11 pool) | batch → `.jld2` | `spike/validation/run_p11_recovery.jl` | exact |
| `spike/validation/p12_sbc.jl` + `run_p12_sbc.jl` | runner (reported) | batch → `.jld2` | `spike/validation/sbc.jl` + `run_sbc.jl` | exact |
| `spike/validation/p12_coverage.jl` + `run_p12_coverage.jl` | runner (reported) | batch → `.jld2` | `spike/validation/run_p11_coverage.jl` + `spike/validation/harness.jl` | exact |
| `spike/validation/run_p12_guards.jl` | runner (reported, S-4 / D-06 guards) | batch → `.jld2` | `spike/validation/run_p11_probe.jl` | role-match |
| `spike/p12/result.jl` (spike-local `SpatialColocResult`-shaped) | model (result type) | container | `spike/p13/result.jl` | exact |
| `.../12-SC3-AMENDMENT.md` | doc (frozen pre-declaration) | — | `.../13-SC2-AMENDMENT.md` | exact |
| `.../12-STAGE1-VERDICT.md` | doc (manual verdict) | — | `.../11-PROBE-VERDICT.md` | exact |
| `spike/test/runtests.jl` **(MODIFIED)** | test wiring | — | `runtests.jl:183-208` (the Phase-13 wiring block) | exact |

---

## 2. Pattern Assignments

### 2.1 `spike/validation/p12_consts.jl` (config, declarative)

**Analog:** `spike/validation/p11_consts.jl` — this is the file R-5 says to mirror **exactly**.
Read in full this session (487 lines, two tiers).

**Two-tier header contract** (`p11_consts.jl:30-48`) — copy the structure and the voice:

```julia
# TWO-TIER STRUCTURE -- READ THIS BEFORE EDITING ANYTHING.
#   TIER 1 is the single guard block in THIS file, committed before the probe.
#   TIER 2 (SC2_SPEARMAN_FLOOR and the Δρ_eq calibration coefficients) is APPENDED
#     LATER as a SECOND guard block keyed on a DIFFERENT sentinel const, carrying a
#     one-line provenance comment per constant naming the probe artifact it came from.
#   NO CONSTANT IN EITHER BLOCK IS EVER *EDITED*. Constants are only ever APPENDED,
#   so the git history of this file is itself the pre-registration audit trail. A
#   diff that MODIFIES a line below is, by construction, a pre-registration breach.
```

**Guard-block idiom** (`p11_consts.jl:79-80`, closing `end` at `:379`; Tier 2 opens at `:416`):

```julia
if !isdefined(@__MODULE__, :P11_DEV_SEED)
    import Random123: Philox4x     # counter-based RNG; the fresh disjoint DEV stream
    ...
end
# ... Tier 2 keyed on a DIFFERENT sentinel:
if !isdefined(@__MODULE__, :SC2_SPEARMAN_FLOOR)
    ...
end
```
Phase 12 uses `:P12_DEV_SEED` for Tier 1 and a distinct Tier-2 sentinel (e.g. `:P12_STAGE1_FLOOR`).

**Forbidden-seed set + executable disjointness** (`p11_consts.jl:86-135`, assertions at `:347-353`):

```julia
const NPE_MASTER_SEED     = 0x0000_0000_00C0_FFEE  # spike TRAINING stream (npe/train_npe.jl:65)
const VAL_MASTER_SEED     = 0x0000_0000_5BC0_FFEE  # spike VALIDATION stream (consts.jl:76)
const VAL_FIX_SEED        = 0x0000_0000_00F1_F7ED  # spike FIXTURE stream (consts.jl:79 = 0xF1F7ED)
const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001  # productionization DATAGEN stream
const CORPUS_MASTER_SEED  = 0x0000_0000_00C0_5EED  # corpus stream (corpus/manifest.csv, D-17)
...
# The v1 and v2 ship-gate seeds are DERIVED, not literal. RECOMPUTE them ... a comment
# naming a derived seed is not evidence.
_derive(master, salt, G) =
    rand(Philox4x(UInt64, (UInt64(master) ⊻ UInt64(salt), UInt64(G))), UInt64)
const PROD_SEED_V2 = Dict{Int,UInt64}(G => _derive(AMEND_MASTER, AMEND_SALT, G) for G in (4,8,16,32))

_p11_forbidden() = (UInt64(0), NPE_MASTER_SEED, VAL_MASTER_SEED, VAL_FIX_SEED,
                    DEFAULT_MASTER_SEED, CORPUS_MASTER_SEED,
                    F2_DEV_SEEDS..., SPIKE_DEV_SEEDS..., P11_RESEARCH_BURNED...,
                    values(PROD_SEED)..., values(PROD_SEED_V2)...)

@assert !(UInt64(P11_DEV_SEED) in _p11_forbidden())
@assert length(unique(_p11_forbidden())) == length(_p11_forbidden())
@assert PROD_SEED_V2 == GateV2.PROD_SEED_V2      # the recomputation reproduces the frozen file
```

**PHASE-12 DELTA THE PLANNER MUST WRITE IN — Phase 13's seeds.** `p11_consts.jl` predates Phase 13
and therefore does **not** forbid it. Read this session from `spike/p13/consts.jl:188-192` and
`:120-159`; all four must enter `_p12_forbidden()`:

```julia
const P13_DEV_SEED = 0x0000_0000_0B13_DE71  # spike/p13/consts.jl:188
const P13_FIX_SEED = 0x0000_0000_0B13_F1F7  # spike/p13/consts.jl:189
const P13_SALT     = 0x2545_F491_4F6C_DD1D  # spike/p13/consts.jl:192 — salt, not seed
const P11_DEV_SEED = 0x0000_0000_0B11_DE71  # spike/validation/p11_consts.jl:144
const RATIO_PAIR_SEED = 0x0000_0000_004A_7107  # src/amortized/train_ratio.jl:60 (p13/consts.jl:111)
# plus the 25-entry P13_BURNED_DEV_SEEDS tuple, spike/p13/consts.jl:120-146
```
Phase 13 additionally maintains a **salt inventory** (`spike/p13/consts.jl:150-159`,
`P13_REPO_SALTS`) and asserts `@assert !(P13_SALT in P13_REPO_SALTS)`. Phase 12 should copy that
list and append `P13_SALT` to it. Naming convention for the fresh seed follows both phases:
`0x0000_0000_0B12_DE71` ("0B12" = phase 12, "DE71" = DEV-1) and a fixture seed `0x0000_0000_0B12_F1F7`.

**Reserved counters, fixtures off the reported stream** (`p11_consts.jl:146-159`):

```julia
p11_rng(counter::Integer = 0) =
    Philox4x(UInt64, (UInt64(P11_DEV_SEED) ⊻ P11_SALT, UInt64(counter)))

const P11_PROBE_COUNTER       = 1   # D-06 simulator-only pre-flight probe
const P11_DATAGEN_COUNTER     = 2   # research-net training-pair generation
...
const P11_FIXTURE_COUNTER     = 99  # FIXTURES ONLY -- never a reported counter
```
Phase-12 counters map 1:1 onto the VALIDATION runners: sim02, stage1_ridge, ell_ridge, eps_ridge,
datagen, sbc, coverage, guards, + `P12_FIXTURE_COUNTER = 99`.

**Threshold-derivation voice** (`p11_consts.jl:213-241`) — a Tier-1 threshold is argued from the
*design*, never from a result, and the arithmetic is recorded verbatim:

```julia
# DERIVATION OF N_MIN (recorded verbatim ...). Require the 90 % Wilson interval at
# p̂ = SC2_COVERAGE_NOMINAL to fit inside a ±SC2_TOST_DELTA band:
#     half-width ≈ z·√(p(1−p)/N) = 1.644854·√(0.09/N) = 0.49346/√N ≤ 0.03
#     ⇒ √N ≥ 16.449 ⇒ N ≥ 270.6 ⇒ N_min = 271
```

**Tier-2 append discipline** (`p11_consts.jl:381-414`, and the product form at `:433`):

```julia
# WRITTEN AS A PRODUCT ON PURPOSE. SC2_SPEARMAN_ATTENUATION (Tier 1, 0.5) is the judgment
# allowance ... P11_PROBE_S_MEASURED is the measurement. Collapsing the two into a single
# decimal would hide which half is which.
const SC2_SPEARMAN_FLOOR = SC2_SPEARMAN_ATTENUATION * P11_PROBE_S_MEASURED
```

**Iteration allowance, declared in advance** (`p11_consts.jl:319-328`) — copy verbatim in spirit
for `P12_ITERATION_ALLOWANCE = 1` (D-12 Stage 2 allows exactly one documented iteration).

**Anti-pattern flagged from the analog:** `p11_consts.jl:75-77` opens a bare `module GateV2 ... end`
**outside** the guard block, with the reason stated at `:70-74` (*"Julia rejects a `module`
expression that is not at top level"*). If `p12_consts.jl` needs the frozen F5 imsize mixture the
same way (`P11_IMSIZE_SET = GateV2.SBC_IMSIZE_SET`, `:287-288`), it must repeat that placement —
and must **read** the mixture rather than retype it (`:283-286`: *"a retyped mixture is a mixture
that can silently drift"*). Note `spike/p13/consts.jl` uses the same idea under the name `_GC`.

---

### 2.2 `spike/test/test_p12_consts.jl` (test, assertion)

**Analog:** `spike/test/test_p11_consts.jl` (252 lines, read in full).

**Purpose header** (`test_p11_consts.jl:21-30`):

```julia
# THIS FILE EXISTS TO MAKE A POST-HOC EDIT LOUD. Every locked value is asserted AS A
# LITERAL here, so changing a constant in `p11_consts.jl` after a run breaks the suite
# rather than silently rewriting the pre-registration. The assertions are deliberately
# redundant with the `@assert` self-checks inside the consts file: the self-checks prove
# internal consistency, these prove the values are the SPECIFIC values that were committed.
```

**Guarded include + literal + executable disjointness** (`:35-36`, `:41-73`):

```julia
isdefined(@__MODULE__, :P11_DEV_SEED) ||
    include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))

@test P11_DEV_SEED === 0x0000_0000_0B11_DE71
@test !(UInt64(P11_DEV_SEED) in _p11_forbidden())
# Pairwise, by name, so a failure says WHICH stream collided.
@test P11_DEV_SEED != NPE_MASTER_SEED
@test all(s -> P11_DEV_SEED != s, values(PROD_SEED_V2))
# No duplicate may hide inside the forbidden set.
@test length(_p11_forbidden()) == length(unique(_p11_forbidden()))
```

**Counter-stream separation is *behavioural*, not documentary** (`:76-85`):

```julia
@test P11_FIXTURE_COUNTER ∉ (P11_PROBE_COUNTER, P11_DATAGEN_COUNTER, ...)
# Distinct counters ⇒ distinct Philox sub-streams off the same key.
@test rand(p11_rng(P11_PROBE_COUNTER), UInt64) != rand(p11_rng(P11_LADDER_COUNTER), UInt64)
@test rand(p11_rng(P11_FIXTURE_COUNTER), UInt64) == rand(p11_rng(P11_FIXTURE_COUNTER), UInt64)
```

**Re-derive the recorded derivation inside the test** (`:120-121`):

```julia
# The recorded N_min derivation reproduces: 0.49346/√N ≤ δ ⇒ N ≥ 270.6 ⇒ 271.
@test ceil(Int, (Z_TWO_SIDED_90 * sqrt(0.09) / SC2_TOST_DELTA)^2) == N_MIN_DERIVED
```

**"Tier 1 has not moved" testset** (`:224-239`) — re-assert the Tier-1 literals inside the Tier-2
testset, so an accidental Tier-1 edit during an append is caught by this file. Copy this whole
testset shape.

**Every testset in the analog ends with a CPU-only check** (`:176-178`, `:248-250`):

```julia
@test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
```

**Wave-0 note:** at Wave 0 Tier 2 does not yet exist. The analog handled this by asserting the
*opposite* until the probe ran — recorded at `test_p11_consts.jl:135-139`: the
`!isdefined(..., :SC2_SPEARMAN_FLOOR)` assertion *"was removed in the SAME commit that appended the
Tier-2 block, and only then."* Phase 12 should open with the negative assertion.

---

### 2.3 `spike/simulator/p12_lattice.jl` (utility, transform) — CAR/GP kernels + DCT-II

**Analog:** partial. No lattice/GP code exists in this repo. The closest **role** analogs are:

- `spike/simulator/ghat.jl` — a frozen numeric map with named knot constants and a block comment
  recording the calibration evidence that produced them (`ghat.jl:44-67`). Copy the
  "frozen-numeric-artifact" documentation discipline, **not** the arithmetic.
- `spike/simulator/prior.jl:39-58` — the const-declaration + prose-justification style for the
  distributional knobs:

```julia
# --- The Turing μ-prior (src/bayes.jl ~274-291), the SIM-02 consistency target ---
const MU_PRIOR = Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)

# --- PRE-DECLARED SIM-02 pass bar (T-02-CAL) -------------------------------------
# Identical to calibration.jl's tolerance; fixed BEFORE the calibration run, NOT
# tuned to pass.
const SIM02_W1_TOL = 0.10
```

**Import pattern (the decoupling-critical part).** From `run_p11_recovery.jl:50-53` — stdlibs are
imported plainly and are *absent from `spike/Project.toml`*:

```julia
using JLD2
using Dates
using Statistics
using LinearAlgebra
```
`p12_lattice.jl` adds `using SparseArrays` on the same basis. **Do not add either to
`spike/Project.toml`** — `runtests.jl:123-126` asserts the exact 16-name dependency set.

**Column-major indexing contract.** Derived from `src/amortized/summary.jl:72-76` (read this
session), which is the reason `idx(i,j) = (j-1)*G + i` is mandatory:

```julia
function encode_d01(M::AbstractMatrix)
    vals = vec(coalesce.(M, 0.0))                 # G²-dim, missing->0 (column-major)
    mask = vec(Float64.(.!ismissing.(M)))         # G²-dim, 1=present / 0=missing
    return vcat(vals, mask)                        # 2·G²-dim; rows 1:G² vals, G²+1:2G² mask
end
```

**Concrete code for the kernels and the mandatory `sqrt.(diag(Σ))` rescale:** `12-RESEARCH.md`
Pattern 1 (measured in the live spike env this session by the researcher). Use it verbatim; there
is no in-repo analog to prefer over it.

---

### 2.4 `spike/simulator/p12_prior.jl` (model, sampling/transform) — the D-05 copula

**Analog:** `spike/simulator/prior.jl` (106 lines, read in full) — **exact**.

**Header voice** (`prior.jl:21-32`) — states the consistency claim and the trap in the same breath:

```julia
# spike/simulator/prior.jl --- SIM-02: the simulator prior π(θ), made provably
# consistent with the Turing μ-prior via the frozen induced-μ inverse ĝ.
#
# Do NOT assert ρ_true == μ (02-PATTERNS Anti-Pattern / Pitfall 2): the nonneg
# softplus mixing + PSF + noise attenuate ρ_true into μ nonlinearly -- ĝ inverts
# exactly that measured map.
```

**Guarded include of the frozen `ghat`** (`prior.jl:37`) — D-05 applies `ghat` **elementwise**, so
the include is unchanged:

```julia
include(joinpath(@__DIR__, "ghat.jl"))   # ghat, GHAT_MU_KNOTS/RHO_KNOTS, GHAT_MU_MIN/MAX
```

**Sampler shape and the append-last θ rule** (`prior.jl:87-105`) — the rule D-08's correlation-length
column and D-01's DCT coefficients must respect:

```julia
# `chromatic_eps` is APPENDED AT THE END of the tuple, never slotted next to
# `shift_dx`/`shift_dy`, because downstream code indexes θ POSITIONALLY
# (`collect(values(θ))` row order; `src/amortized/ood.jl` `_theta_tuple` reads
# `v[1]..v[7]`). Appending keeps every existing index valid.
function sample_prior(rng::AbstractRNG)
    μ_star = rand(rng, MU_PRIOR)
    return (
        ρ_true           = ghat(μ_star),
        ...
        chromatic_eps    = rand(rng, CHROMATIC_PRIOR),   # D-09: appended LAST (positional θ)
    )
end
```

**Pre-declared-tolerance pattern for the per-region SIM-02 claim** — `prior.jl:42-46`, quoted above
in §2.3. Phase 12 declares `P12_SIM02_W1_TOL_PERREGION` in `p12_consts.jl` **before** the field
simulator runs, and asserts `maximum(w1_per_region) ≤ tol` (the **max**, not the mean).

**R-8 correction the analog cannot supply:** the per-cell rescale by `sqrt.(diag(Σ))` and the
lag-1-correlation reparametrization are new — take them from `12-RESEARCH.md` Patterns 1-2.

---

### 2.5 `spike/simulator/forward.jl` **(MODIFIED — stage 1, D-06)**

**Analog: the file's own stage 6.** The D-09/D-10 chromatic edit (commit `ca02b0e`) is the exact
precedent for "add a physics term to a stage without breaking anything downstream".

**Backward-compatible θ read** (`forward.jl:131-135`) — the pattern a ρ-**field** θ must copy so a
scalar-ρ θ still works:

```julia
# D-09 backward compatibility (SC1e): chromatic_eps is read through ONE defensive local
# binding, computed here so both the guard below and stage 6 share it. A legacy 7-field θ
# ... is therefore a valid input meaning "no chromatic aberration".
chromatic_eps_val = hasproperty(θ, :chromatic_eps) ? θ.chromatic_eps : 0.0
```

**Stage 1 as it stands today** (`forward.jl:162-171`) — the exact lines D-06 replaces with the
per-pixel `A * ρfield * B'` form:

```julia
    # --- (1) shared-latent correlated smooth Gaussian fields (D-15) -------------
    L  = _smooth_field(rng, imsize)
    ε1 = _smooth_field(rng, imsize)
    ε2 = _smooth_field(rng, imsize)
    a  = sqrt(abs(ρ))           # shared-component weight
    b  = sqrt(1.0 - abs(ρ))     # private-component weight  (a² + b² = 1)
    # sign(ρ) flips channel-2's shared component so negative ρ ⇒ anti-correlation.
    ch1 = _softplus.(a .* L .+ b .* ε1)
    ch2 = _softplus.(sign(ρ) .* a .* L .+ b .* ε2)
```

**Stage 6 is untouched by Phase 12** (`forward.jl:214-218`) — D-06 stacks *before* it, and the S-4
radial confound originates here:

```julia
    s = 1.0 / (1.0 + chromatic_eps_val)
    c = map(ax -> (first(ax) + last(ax)) / 2, axes(ch2))
    A = Translation(θ.shift_dy, θ.shift_dx) ∘ recenter(LinearMap([s 0.0; 0.0 s]), c)
    shifted = warp(ch2, A, axes(ch2); method = BSpline(Linear()), fillvalue = BG_FLOOR)
    ch2 = Matrix{Float64}(collect(shifted))
```

**Two hard "do not touch" constants** (`forward.jl:66-76`), restated by `12-RESEARCH.md` Pattern 1:

```julia
const σ_psf      = 1.3
# STRUCT_σ: field-smoothing scale. MUST be ≫ σ_psf so the induced spatial
# correlation survives PSF+noise (D-15 ...).
const STRUCT_σ   = 6.0
const BG_FLOOR   = 0.02
```
`_smooth_field` (`forward.jl:94-99`) is the **pixel-resolution texture generator**, not a lattice
prior. Do not repurpose it; do not change `STRUCT_σ`.

**Source-level tripwire that the modification must not break** (`test_stage6_regression.jl:94-96`):

```julia
        # ONE resampling pass, asserted on the simulator source itself.
        code = _strip_comment_lines(read(FORWARD_SRC_PATH, String))
        @test count("warp(", code) == 1
```

---

### 2.6 `spike/test/capture_p12_golden.jl` + `test_p12_prior.jl` (golden regression)

**Analog:** `spike/test/capture_p11_golden.jl` (136 lines) and `spike/test/test_stage6_regression.jl`
(197 lines) — both read in full. This pair is the **capture-BEFORE-you-edit** pattern D-06 needs.

**Timing discipline header** (`capture_p11_golden.jl:22-34`):

```julia
# TIMING DISCIPLINE. THIS SCRIPT MUST RUN, AND ITS ARTIFACT MUST BE COMMITTED, BEFORE A SINGLE
# BYTE OF `spike/simulator/forward.jl` STAGE 6 ... CHANGES. The whole value of the fixture is
# that it records what the simulator produced BEFORE ... a fixture regenerated AFTERWARDS would
# compare the new code against itself and pass vacuously. That failure mode is silent, so it is
# made DETECTABLE two ways: the capture-time commit sha is embedded INSIDE the JLD2 as
# `phase11_base_sha`, and the script REFUSES TO RUN if the working tree is dirty ...
```

**The dirty-tree refusal + sha embed** (`capture_p11_golden.jl:77-93`):

```julia
    dirty = readchomp(`git status --porcelain -- spike/simulator src/amortized/simulator.jl`)
    @assert isempty(dirty) """
        capture_p11_golden: the working tree is DIRTY for the simulator files this fixture
        claims to precede: ..."""
    base_sha = readchomp(`git rev-parse HEAD`)
    @assert length(base_sha) == 40 "capture_p11_golden: expected a 40-char sha, got $base_sha"

    channels = [simulate_pair(p11_rng(P11_FIXTURE_COUNTER + k), THETA_GOLDEN;
                              imsize = P11_GOLDEN_IMSIZE) for k in P11_GOLDEN_KEYS]
```

**Explicit θ literal, never drawn from the prior** (`capture_p11_golden.jl:54-65`) — the reason
transfers exactly to Phase 12, whose prior sampler's RNG consumption changes the moment a field is
drawn:

```julia
# An EXPLICIT 7-field θ literal. It is deliberately NOT drawn from the prior sampler: that
# sampler's RNG consumption changes the instant the 8th (chromatic) field is added, so a golden
# captured through it would differ after the edit for a reason that has nothing to do with
# stage 6, and the regression check would be meaningless.
const THETA_GOLDEN = (ρ_true = 0.7, spillover = 0.1, autofluorescence = 0.05,
                      label_efficiency = 0.85, shift_dx = 0.37, shift_dy = -0.82, noise = 0.5)
```

**Atomic write + integrity reopen** (`capture_p11_golden.jl:96-117`):

```julia
    mkpath(dirname(path)); tmp = path * ".tmp"
    jldsave(tmp; schema_version = 1, theta = THETA_GOLDEN, ..., phase11_base_sha = base_sha,
            generated = string(Dates.now(Dates.UTC)) * "Z", caption = "...")
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "channels") "... integrity check failed"
    end
    mv(tmp, path; force = true)   # atomic commit
```

**Consumer side — exact `==`, and refuse a regenerated fixture** (`test_stage6_regression.jl:23-34`,
`:112`, `:129-134`):

```julia
# EXACT EQUALITY IS THE PRE-REGISTERED CRITERION AND IT IS NOT NEGOTIABLE. ... DO NOT SOFTEN A
# FAILING ASSERTION TO A TOLERANCE COMPARISON IN THIS FILE.
            @test collect(composed) == collect(legacy)     # exact `==` (P11_STAGE6_EXACT)
...
        # A fixture REGENERATED after the ε edit would carry the 8th field and would compare
        # the new code against itself, passing vacuously. Refuse that fixture outright.
        @test !hasproperty(golden_theta, :chromatic_eps)
        @test length(golden_theta) == 7
```
Phase-12 analogue: the golden θ must **not** carry the field row, and the constant-ρ path must
reproduce the golden bytes with `==`. The escape hatch is a pre-registered
`P12_STAGE1_TOLERANCE_FALLBACK` whose **use is a recorded deviation** (`p11_consts.jl:314-317`).

---

### 2.7 `spike/npe/p12_architecture.jl` (model, transform)

**Analogs:** `spike/npe/architecture.jl:68-115` (the flow builder) and `spike/npe/p11_architecture.jl`
(the spike-lane wrapper).

**The v0.2.1 API gotcha and the flow construction** (`spike/npe/architecture.jl:104-114`):

```julia
    # MLP summary net: Dense(d_in→width, gelu) then (depth-1) width→width gelu
    # blocks, then a linear width→dstar projection to the learned summaries.
    layers = Dense[Dense(d_in, width, gelu)]   # concretely-typed (IN-06), not Any[...]
    for _ in 2:depth; push!(layers, Dense(width, width, gelu)); end
    push!(layers, Dense(width, dstar))
    network = Chain(layers...)

    # NormalisingFlow INSTANCE, passed POSITIONALLY. `num_coupling_layers` deepens the
    # coupling stack; `depth`/`width` flow through CouplingLayer→AffineCouplingBlock→MLP
    q = NormalisingFlow(D; num_summaries = dstar, num_coupling_layers = num_coupling_layers,
                        depth = flow_depth, width = flow_width)
    return PosteriorEstimator(network, q)
```

**Guard style** (`architecture.jl:92-99`) — copy these `ArgumentError` guards into the CNN builder:

```julia
    d_in  >= 1 || throw(ArgumentError("build_estimator: d_in must be ≥ 1, got $d_in"))
    dstar >= D || throw(ArgumentError("build_estimator: dstar ($dstar) must be ≥ D ($D)"))
```

**Capacity-held-fixed wrapper** (`p11_architecture.jl:341-373`) — Phase 12 **deviates** here (a CNN
replaces the MLP), so R-7's lean topology must be written literally into the plan, and the deviation
must be named. The analog's own words (`p11_architecture.jl:350-359`):

```julia
**NO ARCHITECTURE CODE CHANGES.** `build_estimator` is input-width-agnostic and its `D` is
already POSITIONAL ... **CAPACITY IS HELD FIXED ON PURPOSE.** Every knob defaults to the
unchanged Phase-5-calibrated value ... so the comparison against the shipped net is
CAPACITY-CONTROLLED and the only things that moved are the four declared deviations.
```
Phase-12 equivalent: `P12_RESEARCH_NET_DEVIATIONS` in `p12_consts.jl`, enumerated with a fixed
length and asserted (`p11_consts.jl:336-341`, `:373`):

```julia
    const P11_RESEARCH_NET_DEVIATIONS = (
        "chromatic_eps appended as the 8th theta column (D-09), so D = 8 rather than 7",
        ...
    )
    @assert length(P11_RESEARCH_NET_DEVIATIONS) == 4
```
Phase-12 declared deviations, from CONTEXT/RESEARCH/R-#: (1) CNN over `(G,G,2,K)` replaces the MLP
over 128 rows (D-02); (2) θ widens to `1 + K + 7 + 1` with DCT coefficients (D-01/D-08); (3) the
lattice-prior field replaces scalar ρ in stage 1 (D-05/D-06); (4) random region masking as a
training augmentation (Pitfall 4); (5) the recorded z-scoring arm (R-3).

**The reshape helper and its ordering contract** — take the code from `12-RESEARCH.md` § Code
Examples (`reshape_summary`), and enforce the mask-row rule from
`src/amortized/summary.jl:85-87` (read-only reference):

```julia
z-scoring a 0/1 mask would re-couple folds through the mask mean, so the mask rows are NEVER
standardized.
```

**Ordering-invariant assertion style for a new input contract** (`p11_architecture.jl:330-333`) —
mirror it for `reshape_summary` (standardize FIRST, then reshape):

```julia
function augment_input(Z128::AbstractMatrix, lam)
    @assert size(Z128, 1) == 128 "augment_input: expected the frozen 128-row standardized " *
        "summary, got $(size(Z128, 1)) rows — standardize FIRST, then append λ (Pitfall 2)"
```

---

### 2.8 `spike/test/test_p12_architecture.jl` (test, unit + smoke)

**Analog:** `spike/test/test_p11_architecture.jl` — testset layout read this session:
`"ported transform agrees with src"` (`:79`), `"bounded transform keeps mass inside the prior box"`
(`:122`), `"encode_lambda"` (`:164`), `"augment_input ORDER invariant"` (`:182`),
`"build_p11_estimator constructs at d_in = 129, D = 8"` (`:208`), `"the four declared deviations
are the pre-registered four (D-04)"` (`:229`), `"model surface ran CPU-only (D-10)"` (`:248`).

Constructor-guard assertions to copy (`:224-226`):

```julia
        @test_throws ArgumentError build_p11_estimator(d_in = 0)
        @test_throws ArgumentError build_p11_estimator(D = 0)
        @test_throws ArgumentError build_p11_estimator(D = 65)   # dstar (64) >= D violated
```

Phase-12 testsets map to VALIDATION rows: `reshape_summary` round-trip (no transpose), mask rows
never z-scored / mask is channel 2, CNN + wide `NormalisingFlow` builds → trains 1 epoch →
`sampleposterior` returns `D×N`, declared-deviation enumeration, CPU-only.

---

### 2.9 `spike/test/test_p12_decoupling.jl` (test, integration)

**Analog:** `spike/test/test_p13_result.jl:190-213` — **exact**, including the git-absent guard:

```julia
    @testset "src/results.jl is untouched (D-01, rename deferred)" begin
        @test P13_RESULTS_RENAME_DEFERRED == true
        @test P13_SRC_UNTOUCHED == true
        # Guarded: a sandbox without git, or a source tree that is not a checkout, skips the
        # executable half rather than failing for an unrelated reason.
        _has_git = Sys.which("git") !== nothing &&
                   (isdir(joinpath(P13_REPO_ROOT, ".git")) ||
                    isfile(joinpath(P13_REPO_ROOT, ".git")))
        if _has_git
            @test success(Cmd(`git diff --quiet HEAD -- src/results.jl`; dir = P13_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P13_REPO_ROOT))
        else
            @info "git unavailable — skipping the executable src/ byte-equality assertion"
        end
```

**Second half — the dependency-set gate.** `spike/test/runtests.jl:114-127` is the Phase-11
precedent, asserted as an **explicit name set** so an addition names itself in the failure:

```julia
        @test Set(keys(Pkg.project().dependencies)) == Set([
            "BenchmarkTools", "CSV", "CairoMakie", "CoordinateTransformations", "DataFrames",
            "Distributions", "Flux", "HypothesisTests", "ImageFiltering", "ImageTransformations",
            "Images", "Interpolations", "JLD2", "NeuralEstimators", "Random123", "StatsBase"])
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
```
Phase 12 adds a `(k)` block in the same voice ("the spatial map adds NO new package; the lattice is
64×64 dense `LinearAlgebra`; `AbstractGPs`/`KernelFunctions`/`GaussianRandomFields` must NOT become
direct deps"), plus `git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml corpus`.

---

### 2.10 `spike/validation/run_p12_stage1_ridge.jl`, `run_p12_ell_ridge.jl`, `run_p12_eps_ridge.jl`

**Analog:** `spike/validation/run_p11_recovery.jl` (263 lines, read in full) — **exact**. This is
the file `12-RESEARCH.md` § Code Examples explicitly adapts.

**Header states the exclusion, the baseline and the "gates nothing" scope** (`:17-48`):

```julia
# THE QUESTION: can ANY estimator recover the sub-pixel shift from the 128-row summary?
#
# ROW 129 IS EXCLUDED FROM THE PREDICTORS, AND THIS IS THE WHOLE POINT OF THE TEST.
# ... the prior-SD ratio across the ladder is exactly LAMBDA_MAX / LAMBDA_MIN = 12.0 by
# construction ... Any diagnostic that admits row 129 therefore measures PRIOR ECHO, not learning.
#
# THE BASELINE: predicting the conditional prior mean ... It is computed EMPIRICALLY from the
# realised shifts in each bin rather than from the analytic formula, so the comparison cannot be
# flattered by a mismatch ...
#
# Ridge is the right tool here precisely because it is weak and transparent: closed-form, no
# tuning theatre, no capacity to memorise.
#
# DECOUPLING (S5): spike-local; reads the frozen pre-registration and the existing training pool;
# writes one new artifact. Trains no network, mutates no constant, touches no `src/` file.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_recovery.jl
```

**Guarded includes + stdlib imports** (`:50-56`):

```julia
using JLD2
using Dates
using Statistics
using LinearAlgebra

isdefined(@__MODULE__, :P11_DEV_SEED)        || include(joinpath(@__DIR__, "p11_consts.jl"))
isdefined(@__MODULE__, :generate_p11_pool)   || include(joinpath(@__DIR__, "..", "data", "p11_generate.jl"))
```

**Leak-free split, train-only standardizer, validation-selected penalty** (`:63-66`, `:75-100`,
`:129-158`):

```julia
const TEST_FRACTION = 0.25
const RIDGE_GRID = [1e-6, 1e-4, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]

function _fit_standardizer(X::AbstractMatrix)
    mu = vec(mean(X; dims = 2)); sd = vec(std(X; dims = 2))
    @inbounds for i in eachindex(sd)
        (isfinite(sd[i]) && sd[i] > 1e-12) || (sd[i] = 1.0)   # zero-variance mask rows pass through
    end
    return mu, sd
end

function ridge_fit(Xs::AbstractMatrix, y::AbstractVector, penalty::Real)
    b  = mean(y); yc = y .- b
    G  = Symmetric(Xs * Xs' + penalty * I)
    w  = G \ (Xs * yc)
    return w, b
end
_rmse(pred, truth) = sqrt(mean(abs2, pred .- truth))

    ntest  = round(Int, TEST_FRACTION * N)
    test_i = (N - ntest + 1):N
    tr_all = 1:(N - ntest)
    nval   = round(Int, 0.2 * length(tr_all))
    ...
    mu, sd = _fit_standardizer(view(S, :, fit_i))          # TRAIN-ONLY fit
        # Select the penalty on VALIDATION, never on test.
        for p in RIDGE_GRID
            w, b = ridge_fit(Xfit, yfit, p)
            v = _rmse(ridge_predict(Xval, w, b), yval)
            v < best_v && ((best_p, best_v) = (p, v))
        end
```

**The positive control that proves the harness is live** (`:209-231`):

```julia
    # --- A positive control -------------------------------------------------------------------
    # rho_true IS strongly encoded in the patch-correlation summary ... If this control failed, the
    # null on the shift would be uninformative -- it would just mean the harness is broken. This is
    # what separates "the summary lacks shift information" from "this script is wrong".
    yfit = vec(TH[1, fit_i]); ...
    ctrl_prior = sqrt(mean(abs2, ytst .- mean(yfit)))
```

**Atomic reported-artifact write** (`:233-259`) — the shape every `run_p12_*.jl` reuses:

```julia
    elapsed = (time() - t_start) / 60
    tmp = RECOVERY_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version   = 1,
        n_total = N, n_fit = length(fit_i), n_val = length(val_i), n_test = length(test_i),
        predictor_rows   = 128,
        row_129_excluded = true,
        ridge_grid       = RIDGE_GRID,
        ...
        elapsed_min      = elapsed,
        generated        = string(Dates.now(Dates.UTC)) * "Z",
        caption          = "... Diagnostic; gates nothing.",
    )
    let d = JLD2.load(tmp); @assert haskey(d, "shift_dx") && haskey(d, "control_rho_ratio"); end
    mv(tmp, RECOVERY_REPORT_PATH; force = true)
```

**Runner-not-imported guard** (`:263`) — makes the file safe to `include` from a test:

```julia
isdefined(@__MODULE__, :P11_RECOVERY_LOAD_ONLY) || main()
```

**Per-file deltas the planner must specify:**
- `run_p12_stage1_ridge.jl` — predictors are the 128 rows **with row r zeroed AND masked**; target
  is the drawn lattice value at region r (D-07 truth); the positive control is the **same ridge with
  row r included**; the pooled table is replaced by a **correlation-length ladder** breakdown
  (the analog's λ-bin loop at `:176-193` is the template). **This one gates (D-12 Stage 1)** — it is
  the only runner whose header must *not* say "gates nothing".
- `run_p12_ell_ridge.jl` — target is ℓ (never conditioned on, so there is no row to exclude);
  **add Moran's I of the 64 continuous rows as an engineered predictor** and report raw-rows and
  raw+Moran arms separately (Pitfall 3). Baseline = predict the prior mean of ℓ. Reported only.
- `run_p12_eps_ridge.jl` — read-only against the **existing** Phase-11 pool
  (`p11_pool_dir(P11_N_PAIRS)` / `load_p11_pool`, `:111-116`); target is θ row 8
  (`chromatic_eps`). Reported only, never gated (R-6). **Must not regenerate or overwrite that
  pool.**

---

### 2.11 `spike/validation/run_p12_coverage.jl` + `p12_coverage.jl`

**Analogs:** `spike/validation/run_p11_coverage.jl` (the coverage runner) and
`spike/validation/harness.jl` (the shared draw→simulate→infer chain that CONTEXT says to *extend*,
not duplicate).

**Harness invariants every consumer inherits** (`harness.jl:28-42`):

```julia
# CENTRAL INVARIANTS (every consumer inherits them):
#   • FROZEN STATS (T-05-01): the net + BOTH transforms (zt, θzt) are loaded ONCE
#     via load_npe; standardization is NEVER re-fit here (no fit(ZScoreTransform)).
#   • CPU-ONLY (D-10 / Pitfall 1): use_gpu = false on EVERY NeuralEstimators call.
#   • DISJOINT SEED (D-02): val_rng threads a Random123 stream keyed off
#     VAL_MASTER_SEED ⊻ VAL_SALT ...
#   • READ-ONLY src/ (hard constraint, CLAUDE.md) ...
```

**The one chain to extend** (`harness.jl:106-113`) — Phase 12 adds a `reshape_summary` step between
`standardize_summary` and `posterior_for`, and a region-masking step before it:

```julia
function draw_simulate_infer(m, rng; imsize = SBC_IMSIZE, N = SBC_L)
    θ   = sample_prior(rng)                                          # θ* ~ π(θ)
    mci = build_mci(simulate_pair(rng, θ; imsize = imsize))          # forward sim → MCI
    Z   = standardize_summary(encode_d01(patch_summary(mci)), m.zt, :min)  # FROZEN m.zt
    draws_std = posterior_for(m.estimator, Z; N = N, use_gpu = false)      # 7×N standardized
    draws     = StatsBase.reconstruct(m.θzt, draws_std)             # 7×N physical (Pitfall 5)
    return (θ = θ, Z = Z, draws = draws)
end
```

**Paired-draw path** (`harness.jl:124-136`) — the template for the sample/control pass that yields
per-region Δρ.

**Interval + coverage helper** (`run_p11_coverage.jl:93-95`, and the print/store block at `:153-197`):

```julia
"90% (or nominal) equal-tailed interval of a draw vector."
function _interval(v::AbstractVector, nominal::Real)
    a = (1 - nominal) / 2
```

**Scope-limiting header** (`run_p11_coverage.jl:43-46`) — copy the voice, then extend it per
Pitfall 5: Phase 12 gates on **calibration AND a proper scoring rule**, and must pre-register that a
coverage-only win with no scoring-rule improvement is **not** a pass.

---

### 2.12 `spike/validation/p12_sbc.jl` + `run_p12_sbc.jl`

**Analog:** `spike/validation/sbc.jl` (read via targeted grep).

**Rank table and PIT/uniformity** (`sbc.jl:156-190`):

```julia
function sbc_ranks(m; M::Integer = SBC_M, L::Integer = SBC_L, ...)
    rank_table = Matrix{Int}(undef, M, 8)
        ...
            rank_table[i, p] = count(<(t.θ[p]), @view t.draws[p, :])
        # Δρ rank from a paired draw (D-01), independent of the marginal ρ_true rank.
        rank_table[i, 8] = count(<(Δρ_star), Δρ_draws)

function sbc_uniformity(ranks::AbstractVector{<:Integer}; L = SBC_L, bins = SBC_BINS)
    u = (ranks .+ 0.5) ./ (L + 1)
```
`HypothesisTests` is imported at `sbc.jl:44` — **never hand-roll KS/χ²**.

**Phase-12 delta (R-2 / S-3 / Pitfall 7-8):** score the *Gaussian-space* DCT coefficients (atom-free)
with the point-null, produce a **second** ρ-space table with randomized ranks, and use a
**nuisance-appropriate equivalence test (TOST)** for the nuisances and high-index deviation
coefficients. The randomized-rank recipe is in `12-RESEARCH.md` Pitfall 7; the TOST constants
pattern is `p11_consts.jl:203-241`.

---

### 2.13 `spike/data/p12_generate.jl` (service, batch)

**Analog:** `spike/data/p11_generate.jl` — **exact**.

**Generative-order header** (`p11_generate.jl:20-33`) — the ordering discipline transfers directly
(Phase 12: draw ℓ / r₁ **first**, then the field conditional on it):

```julia
# THE GENERATIVE STORY, IN EXACTLY THIS ORDER. FOR EVERY TRAINING SAMPLE i:
#     lambda_i        ~ Uniform(LAMBDA_MIN, LAMBDA_MAX)   <- THE UNCERTAINTY LEVEL, DRAWN **FIRST**
#     shift_dx_i      ~ Uniform(-lambda_i, +lambda_i)     <- DRAWN **CONDITIONAL ON** lambda_i
#     ...
# THIS ORDER IS THE ENTIRE POINT OF THE FILE. IT EXISTS TO PREVENT **PITFALL 1 / RISK R4** ...
# A NULL THAT LOOKS EXACTLY LIKE "BELOW THE FIXED 8x8 SUMMARY'S RESOLUTION" WHILE ACTUALLY BEING
# AN IMPLEMENTATION BUG.
```

**Raw storage / no z-scoring in datagen** (`p11_generate.jl:43-51`):

```julia
# WHAT THIS FILE DELIBERATELY DOES **NOT** DO. It stores the RAW 128-row summary and the raw
# lambda separately. It does not z-score anything and it does not build the 129-row input.
# Fitting the frozen summary transform is the loader's train-only job (leak-freedom) ...
```
Phase-12 equivalent: store the **drawn G×G lattice values** (D-07 truth) and the raw 128 rows;
reshape and z-score happen downstream.

**Include-order comment + guarded includes** (`p11_generate.jl:53-63`) and the salt-derivation block
(`P11_DATAGEN_SALT`, `:64-80`) transfer verbatim in shape.

**New pool dir, do not overwrite P11's** — `12-RESEARCH.md` § Runtime State Inventory.

---

### 2.14 `spike/p12/result.jl` (model, container)

**Analog:** `spike/p13/result.jl` — **exact**, including the reason the `src/` sketch is left alone
(`result.jl:35-43`):

```julia
# (b) WHERE THE EXECUTABLE TYPE LIVES, AND WHY THE SKETCH IS NOT EDITED. D-01 scopes all
#     Phase-13 work to spike/, and CLAUDE.md requires src/ to stay provably untouched during the
#     spike. So the EXECUTABLE subtype is defined HERE, in spike/p13/, reached to its supertype
#     by a READ-ONLY include of src/results.jl -- the same read-only reach spike/contract.jl:46-47
#     already uses ... The comment-only sketch block at src/results.jl:161-192 stays BYTE-UNCHANGED.
#
# (c) PHASE 7 D-02 IS STILL SATISFIED. D-02 requires every new result variant to slot in as a NEW
#     subtype of the real `AbstractColocResult`, never as bolted-on fields on an existing struct.
```
Struct declaration form: `struct ThreeHypothesisColocResult <: AbstractColocResult` (`:142`).
Phase-12: `struct SpatialColocResultSpike <: AbstractColocResult` with `region_delta_rho`,
`region_sd`, `ood`, `calibration`, `meta` per the `src/results.jl:180-191` sketch, plus
`delta_rho_map` / `uncertainty_map` accessors.

**Degenerate-region discipline** — from the read-only reference `src/amortized/local_map.jl:42-53`:

```julia
The documented finite sentinel written into `LocalColocMap.delta_rho` for a sub-tile that could
NOT be scored ... The value is `0.0` ("no evidence of a sample-vs-control difference"), and EVERY
sentinel tile also carries `ood_flag = true`, so a sentinel is never silently indistinguishable
from a measured Δρ of zero (T-7-04: a degenerate tile must never crash the map, and must never
masquerade as a finding).
const LOCAL_MAP_SENTINEL = 0.0
```
And the "missing operating point stays honestly inert" rule (`local_map.jl:88-91`): a `false` flag
with no recorded threshold means **NOT CHECKED**, not "in distribution".
**`LocalColocMap` is NOT a coverage baseline** (`local_map.jl:58-61`: *"carries no posterior draws,
no Bayes factor and no per-region uncertainty"*).

---

### 2.15 `12-SC3-AMENDMENT.md` and `12-STAGE1-VERDICT.md` (docs)

**Analog:** `.planning/phases/13-.../13-SC2-AMENDMENT.md` — **exact**. Copy its §0 provenance
disclosure table, which makes the ordering checkable afterwards:

```markdown
**Status:** **FROZEN SPECIFICATION. No Phase-13 result exists.**

| Field | Value |
|---|---|
| Date written | **2026-07-27** |
| Commit sha at time of writing | **`df16e76941cefd28147efb4d4c91278597836cbd`** |
| Phase-13 results in existence at that sha | **none** — no `spike/p13/*report*.jld2`, no trained net |

The test of legitimacy applied throughout is the precedent's, quoted verbatim:
> *Would this change have been made, with this justification, by someone who had seen the gate's
> code but none of its outputs?*
```
Phase 12's SC3 amendment (D-09: predictive, not parameter, coverage; matched ablation, not
`LocalColocMap`) must be committed **before `run_p12_coverage.jl` is ever executed**.

**Verdict doc analog:** `.planning/phases/11-.../11-PROBE-VERDICT.md` — the D-12 Stage-1
PROCEED/DESCOPE call, with the numbers quoted against the Tier-1 threshold.

---

### 2.16 `spike/test/runtests.jl` (MODIFIED)

**Analog:** its own Phase-13 wiring block (`runtests.jl:183-208`), which documents both the
dependency ordering and the abort hazard:

```julia
# Phase-13 three-hypothesis evidence net: the unit testsets run in the same harness ...
# THE CORRECTION ARM IS LAST, AND DELIBERATELY SO. Its two pre-registered
# `max |Delta log BF| <= P13_F5_MAXABS_TOL` assertions are MEASURED MISSES ... so its outer
# `@testset` throws at the end of the file. A thrown testset aborts the remaining includes, so any
# sibling placed after it would silently never run.
include(joinpath(@__DIR__, "test_p13_consts.jl"))
...
include(joinpath(@__DIR__, "test_p13_calibration.jl"))   # line 207
include(joinpath(@__DIR__, "test_p13_correction.jl"))    # line 208 — THROWS
```

**Insertion point: between line 207 and line 208.** Suggested Phase-12 order (dependency-first,
mirroring the Phase-13 comment style): `test_p12_consts.jl` → `test_p12_lattice.jl` →
`test_p12_prior.jl` → `test_p12_architecture.jl` → `test_p12_decoupling.jl`.

Add a `(k)` resolve-risk block inside the existing `"D-04 CPU-only"` testset
(`runtests.jl:114-151` is the template) rather than a new top-level testset.

---

## 3. Shared Patterns

### 3.1 Guarded, order-documented includes
**Source:** `spike/validation/harness.jl:48-58`; identical idiom in `p11_generate.jl:57-63`,
`capture_p11_golden.jl:45-48`, `test_stage6_regression.jl:62-67`.
**Apply to:** every new `.jl` file in this phase.

```julia
# --- ORDER MATTERS: consts first, then the trainer/inference surface ..., then the simulator
#     halves, contract, encode, and the seeding primitives whose salt idiom val_rng mirrors.
#     Guarded for idempotency.
isdefined(@__MODULE__, :VAL_MASTER_SEED) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :load_npe)        || include(joinpath(@__DIR__, "..", "npe", "train_npe.jl"))
```
Note the *reason* recorded at `p11_generate.jl:53-56`: the consts file binds the reserved seeds at
**UInt64 width**, and a later file rebinding one at a narrower width is a hard `const` error — so
consts must come first.

### 3.2 Single-guard-block const files
**Source:** `spike/validation/consts.jl:42` (`if !isdefined(@__MODULE__, :SBC_M)`), with the reason at
`:38-40`: *"Guarded as ONE block keyed on :SBC_M so re-inclusion under runtests.jl ... is a silent
no-op -- a redefinition to the same value would otherwise warn on a `const`."*
**Apply to:** `p12_consts.jl` (**three** blocks, three distinct sentinels: `:P12_CHOSEN_PRIOR` and
`:P12_N_LOW` opened by 12-15, `:P12_FISHERZ_NEFF` by 12-16 — `:P12_N_LOW` added when the user confirmed
`n_low` as a Tier-2 measurement on 2026-07-29).

### 3.3 CPU-only assertion in every testset
**Source:** `test_p11_consts.jl:176-178`, `:248-250`; `test_p13_result.jl:215-217`.
**Apply to:** every Phase-12 testset.

```julia
    @testset "... ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end
```

### 3.4 Atomic artifact write with integrity reopen
**Source:** `run_p11_recovery.jl:234-257`; `capture_p11_golden.jl:96-117`; same in
`run_p11_coverage.jl:197`.
**Apply to:** every `run_p12_*.jl` and every fixture capture.
Required keys: `schema_version`, `generated` (UTC ISO + `"Z"`), `elapsed_min`, `caption`, plus the
run's own fields. Then `let d = JLD2.load(tmp); @assert haskey(d, "<key>"); end` and
`mv(tmp, path; force = true)`.

### 3.5 Reported vs gated
**Source:** `run_p11_recovery.jl:42` (*"This is a DIAGNOSTIC, not a deliverable, and it gates
nothing"*) and `run_p11_coverage.jl:43-46` (*"no ratio, no bar, no gate"*).
**Apply to:** every Phase-12 runner **except** `run_p12_stage1_ridge.jl`. VALIDATION § "Reported vs
gated" is explicit that atom mass, ε identifiability and vacuous-column shrinkage are
reporting-only — recording them as gates recreates the Phase-7 over-powered-χ² failure.

### 3.6 Vacuous-column reporting (F3)
**Source:** CONTEXT `<code_context>` § Established Patterns — `shrinkage = post_sd / prior_sd` with a
boolean `vacuous` flag, reporting-only.
**Apply to:** the correlation length (D-08), every deviation coefficient, and `chromatic_eps`
(R-6 / Pitfall 1).

### 3.7 Licence header
**Source:** the AGPL `#= ... =#` block at the head of `spike/simulator/prior.jl:1-19`,
`spike/validation/consts.jl:1-19`, `spike/test/test_p11_consts.jl:1-19`. Newer runners use the
`#`-prefixed variant (`run_p11_recovery.jl:1-15`). Either is in-repo; match the neighbouring file.

### 3.8 Pathspec-scoped commits
**Source:** `12-RESEARCH.md` § Runtime State Inventory — Phase 13 is executing concurrently, and
Phase 11's closure records that concurrent agents race on the **git index**.
**Apply to:** every Phase-12 commit: `git commit -- <paths>`, never a bare `git commit`.

---

## 4. No Analog Found

| File / concern | Role | Data Flow | Reason |
|---|---|---|---|
| `spike/simulator/p12_lattice.jl` — CAR/GP kernel + DCT-II basis arithmetic | utility | transform | No lattice, GP, CAR or DCT code exists anywhere in the repo. The *file shape* copies `prior.jl`/`ghat.jl` (§2.3); the *arithmetic* must come from `12-RESEARCH.md` Pattern 1 (measured in the live spike env) and Pattern 5. |
| `spike/validation/run_p12_guards.jl` — radial-energy / offset-grid / ε=0 guards | runner | batch | Only a role-match: `run_p11_probe.jl` is a simulator-only pre-flight sweep with pre-registered rungs and a frozen metric symbol (`P11_PROBE_METRIC = :paired_l2_rows_1_64`, `p11_consts.jl:260`) — copy the *frozen-metric-symbol* discipline. The radial-basis regression and the half-cell offset (`range(1-0.5, G+0.5; length=m)`) are new; see `12-RESEARCH.md` Pattern 3 and Pitfall 1. |
| Leave-region-out predictive scoring rule (log predictive density / CRPS) | — | — | No proper scoring rule exists in the repo (the SBC/coverage stack computes ranks, ECE and interval coverage only). Pitfall 5 requires one. New code; pre-register the Fisher-z observation-noise model per `12-RESEARCH.md` § D-09 construction. |

---

## 5. Conventions

**Convention derivation skipped** (`gsd-tools verify conventions --derive` returned
`{"skipped": true, "reason": "no-readable-files"}` both repo-wide and at `--scope spike`; the
deriver does not read `.jl`). The table below is **hand-derived** from the ~15 files read this
session and is therefore evidence-based but not entropy-scored.

| Axis | Dominant | Share (hand-counted) | Entropy | Status |
|---|---|---|---|---|
| File-name casing | `snake_case.jl`, phase-prefixed `p11_*` / `p13_*` / `run_p11_*` / `test_p11_*` | ~100% of spike files | low | **named contract** |
| Identifier casing | `snake_case` functions, `SCREAMING_SNAKE` consts, `CamelCase` types; leading `_` = file-private helper | ~100% | low | **named contract** |
| Export style | **No `module`, no `export`** in spike code — flat top-level functions loaded by guarded `include` (`harness.jl:41-42`: *"Flat top-level functions (sibling style of infer.jl -- no module wrapper)"*). `module` appears only as a deliberate isolation wrapper (`p11_consts.jl:75-77` `GateV2`, `p13/consts.jl` `_GC`). | ~100% of spike | low | **named contract** |
| Import style | `using <Pkg>` at file head (never `import` except `import Random123: Philox4x`), then `isdefined(@__MODULE__, :sym) \|\| include(joinpath(@__DIR__, ...))` for repo-local code | ~100% | low | **named contract** |

**Contested hotspots (author's choice).** One genuine split exists and it is **intentional**: the
**licence-header dialect**. Older files use the AGPL `#= ... =#` block (`prior.jl:1-19`,
`consts.jl:1-19`, `test_p11_consts.jl:1-19`); newer Phase-11 runners use the `#`-prefixed short form
(`run_p11_recovery.jl:1-15`). Each half is internally consistent *per vintage*, contested only when
the whole `spike/` tree is pooled. This is the same class as the CJS↔SDK dual-resolver split named
in the GSD prototype (`bin/lib/**` CJS vs `sdk/src/**` ESM): consistent per-directory, contested
repo-wide. **Reviewers and planners should match the local neighbourhood** — a new
`spike/validation/run_p12_*.jl` takes the `#` short form of its `run_p11_*` siblings; a new
`spike/simulator/p12_*.jl` takes the `#= =#` block of `prior.jl`/`forward.jl`.

A second, weaker split: `spike/p13/*.jl` lives in a **phase-named subdirectory** while Phase 11 used
**phase-name-prefixed files in the shared directories** (`spike/validation/p11_consts.jl`). Phase 12
should follow the **Phase-11 prefix form**, because VALIDATION § Wave 0 names
`spike/validation/p12_consts.jl`, `spike/simulator/p12_lattice.jl`, `spike/npe/p12_architecture.jl`
explicitly. The one exception is the result type, where the Phase-13 form (`spike/p12/result.jl`)
is the better analog — record the mixed placement as a deliberate choice.

---

## 6. Metadata

**Analog search scope:** `spike/validation/`, `spike/simulator/`, `spike/npe/`, `spike/test/`,
`spike/data/`, `spike/p13/`, `src/amortized/`, `src/results.jl`,
`.planning/phases/11-*/`, `.planning/phases/13-*/`.

**Files read in full:** `spike/validation/consts.jl`, `spike/validation/p11_consts.jl`,
`spike/validation/run_p11_recovery.jl`, `spike/validation/harness.jl`, `spike/simulator/prior.jl`,
`spike/test/runtests.jl`, `spike/test/test_p11_consts.jl`, `spike/test/capture_p11_golden.jl`,
`spike/Project.toml`, `12-CONTEXT.md`, `12-VALIDATION.md`.
**Files read in part:** `spike/simulator/forward.jl:55-230`, `spike/simulator/ghat.jl:38-75`,
`spike/npe/architecture.jl:55-115`, `spike/npe/p11_architecture.jl:283-373`,
`spike/test/test_stage6_regression.jl:20-140`, `spike/test/test_p13_result.jl:185-219`,
`spike/p13/consts.jl:100-229`, `spike/p13/result.jl:20-60`, `spike/data/p11_generate.jl:20-80`,
`src/amortized/summary.jl:40-104`, `src/amortized/local_map.jl:19-108`,
`spike/validation/run_p11_coverage.jl:17-50`, `13-SC2-AMENDMENT.md:1-60`, `12-RESEARCH.md` (full).

**Pattern extraction date:** 2026-07-27
