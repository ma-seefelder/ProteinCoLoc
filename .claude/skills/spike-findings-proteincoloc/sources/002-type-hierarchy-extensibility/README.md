---
spike: 002
name: type-hierarchy-extensibility
type: analysis
validates: "Given AbstractMultiChannelImage + MultiChannelImageStack + CoLocResult, when audited, then a real interface contract and structure hierarchy are proposed that make the types substitutable and reusable"
verdict: VALIDATED
related: [001, 003]
tags: [extensibility, types, reuse]
---

# Spike 002: Type hierarchy & extensibility

## What This Validates

Given the image/result type system, when audited against Julia interface idioms, then the
half-finished abstractions are identified and a concrete, substitutable type/interface hierarchy
is proposed (priors for the v2.0 productionization, which will need synthetic/alternative stacks).

## Research / Method

Read `LoadImages.jl` (abstract type + accessors + the two structs), `bayes.jl` (`CoLocResult`),
`ProteinCoLoc.jl` (exports), `utils.jl`/`main.jl` (construction sites). Each finding adversarially
verified; several finder citations to `colocalization.jl` were **corrected to `bayes.jl`**.

## Findings (all verified)

| ID | Finding | Verdict | Sev | Effort |
|----|---------|---------|-----|--------|
| TYPE-2 | `hasfield`/`throw` accessors are a fake interface that leaks field names into the supertype | ✓ confirmed | Med | S |
| TYPE-3 | The accessor interface is bypassed by direct `img.data`/`img.:name` access everywhere | ✓ confirmed | Med | M |
| TYPE-8 | `MultiChannelImageStack` is `mutable` but never reassigned | ✓ confirmed | Low | S |
| TYPE-1 | `CoLocResult.advi_result::Any` + abstract container fields; one struct for prior **and** posterior | ⚠ needs-nuance | High→clarity | M |
| TYPE-4 | `MultiChannelImageStack` roots in no abstract type; partial iteration protocol | ⚠ needs-nuance | Med | M |
| TYPE-5 | `get_images` builds a Vector with a stale 3-of-4 parametric spelling (abstract eltype) | ⚠ needs-nuance | Med | S |
| TYPE-6 | `correlation` rebuilds a `Dict` of functions per call instead of dispatching | ⚠ needs-nuance | Med | S |
| TYPE-7 | Pairwise 2-channel analysis is hardcoded across the data/result layer | ⚠ needs-nuance | Low | L |

### TYPE-2 — fake `hasfield`/`throw` interface — ✓ confirmed (med)
`LoadImages.jl:31-85`. Six accessors are defined **on the abstract type**, each doing
`hasfield(typeof(img),:data) || throw(...); return img.data`. This inverts the Julia contract: the
supertype hard-codes the concrete field names, `MultiChannelImage` implements **none** of them (it
free-rides), and the six are near-identical boilerplate. Fix: erroring generic fallbacks +
concrete impls on `MultiChannelImage` (optionally `@eval`-generated). **Keep `throw(ArgumentError)`**
(not `error`) to preserve the exception type. `MultiChannelImage` is currently the only subtype.

### TYPE-3 — interface bypassed by direct field access — ✓ confirmed (med)
`bayes.jl:169-170`, `LoadImages.jl:260-261`, `plot.jl:305,314-315`, `utils.jl:140,144,149,153`,
`main.jl:142,212`. Most hot sites use `image.data[…]` / `image.:name` directly; even `_apply_mask!`,
which dispatches on `AbstractMultiChannelImage`, reaches into `img.channels`/`img.data` — so a
conforming new subtype would still break. Fix: route every field touch through accessors; add an
`image_name!` setter for the one mutation at `main.jl:142`. (Necessary but not alone sufficient for
full substitutability — e.g. `_local_correlation_plot` takes an untyped `img`.)

### TYPE-1 — `CoLocResult` typing — ⚠ needs-nuance (reframe to clarity)
`bayes.jl:75-82` (struct), `bayes.jl:310-318`/`331-335` (construction — **not** colocalization.jl).
`advi_result` is unannotated (`::Any`) and holds a `Chains` for the prior object but an ADVI `q`
for the posterior; the field literally named `posterior` holds **prior** draws in the prior object.
**Correction:** the struct is already *concrete*; `advi_result`/`img`/`control` are **never read on
any hot path** (the hot field `.posterior` is already a concrete `DataFrame`), so this is a
**clarity/extensibility** fix, not a speedup. Behavior-preserving form (no field renames):
```julia
abstract type AbstractCoLocResult end
struct CoLocResult{St<:MultiChannelImageStack, A} <: AbstractCoLocResult
    img::St; control::St; channels::Vector{Int64}; num_patches::Int64
    posterior::DataFrame; advi_result::A      # A = Chains (prior) or ADVI q (posterior)
end
```
Renaming `posterior`→`samples` would touch ~20 sites and break `typeof(x)==CoLocResult` tests
(`test/runtests.jl:203-298`) — defer or stage with `isa`.

### TYPE-4 — `MultiChannelImageStack` not abstract-rooted — ⚠ needs-nuance (med)
`LoadImages.jl:148-151,153,174-196`. The stack subtypes nothing and hand-rolls
`length`/`getindex`/`iterate` but omits `eltype`/`firstindex`/`lastindex`/`eachindex`/`size`.
**Correction:** `map`/`collect` actually still work (collect yields `Vector{Any}`); only
`stack[end]` and the index helpers break. Recommended fix: `MultiChannelImageStack{T} <:
AbstractVector{T}` (define `size`+`getindex`, delete the three hand-written methods) — restores the
full protocol for free. An `AbstractImageStack` supertype + relaxing `::MultiChannelImageStack`
signatures enables synthetic stacks (directly useful for v2.0 productionization, deferred per the
decoupling constraint). Second hardcoded signature is at `bayes.jl:234-236` (not colocalization.jl).

### TYPE-5 — stale 3-of-4 parametric spelling — ⚠ needs-nuance (med)
`utils.jl:63` writes `Vector{MultiChannelImage{Float64,String,Float64}}` (3 of 4 params → UnionAll
eltype); `main.jl:140` is the fully-abstract `Vector{MultiChannelImage}`. Fix: `Vector{
MultiChannelImage{Float64,String,Float64,Int}}` (or an inferred comprehension). **Correction:** the
proposed `eltype(control_unshuffled)` is wrong — there is no `eltype(::MultiChannelImageStack)`, so
it returns `Any`; use `eltype(control_unshuffled.img)`. Benefit is robustness, not measurable speed
(I/O-bound, few images). Same target as PERF-7.

### TYPE-6 — `correlation` method via `Dict` — ⚠ needs-nuance (med)
`colocalization.jl:221-239`. The per-call `Dict` allocation is real; **correction:** it infers as
`Dict{Symbol,Function}` (not `Any`), and the proposed `_corfun(Val(method))` sketch is **worse**
(runtime `Symbol`→`Val` is type-unstable → still dynamic). To truly specialize, use a **function
barrier**: `correlation(x,y;method)= _correlation(x,y,_corfun(Val(method)))` with
`_correlation(x,y,f::F) where F`. Cheapest behavior-preserving win: hoist a module-level
`const COR_FUNCS = Dict{Symbol,Function}(...)`. This is the *same* fix as PERF-1 — do them together.

### TYPE-7 — 2-channel hardcoding — ⚠ needs-nuance (low, L effort)
`bayes.jl:169-172`, `main.jl:123`, `bayes.jl:267-303` (`@model` — **not** colocalization.jl),
`utils.jl:290-303`. **Correction:** pairwise N-channel analysis **already works** via the
`channel_selection=false` fan-out (`main.jl:181-205`); the true gap is the absence of a **joint**
model. Retitle to "blocks joint N-channel modeling." The `.μ_sample`/`.μ_control` column convention
is consumed by `compute_BayesFactor`, `generate_txt`, and three plot functions, so a per-pair-keyed
`CoLocResult` is a larger refactor than the finding implied. Defer the joint-model rewrite.

### TYPE-8 — drop `mutable` on the stack — ✓ confirmed (low)
`LoadImages.jl:148`. No `stack.img`/`stack.name` reassignment exists anywhere; the contained
`Vector` stays mutable, so in-place pixel/mask ops are unaffected. Drop `mutable`. Pairs with TYPE-4.

## Proposed target hierarchy (sketch)
```julia
abstract type AbstractMultiChannelImage end          # + erroring-fallback interface (TYPE-2)
mutable struct MultiChannelImage{T,S,F,I} <: AbstractMultiChannelImage  # mutable: name is reassigned
abstract type AbstractImageStack{T<:AbstractMultiChannelImage} end       # TYPE-4
struct MultiChannelImageStack{T} <: AbstractVector{T}                     # immutable (TYPE-8)
abstract type AbstractCoLocResult end                                     # TYPE-1
struct CoLocResult{St,A} <: AbstractCoLocResult
```

## Signal for the Build
TYPE-2/-3 (real interface + route through accessors), TYPE-8 (drop mutable), and TYPE-1
(parametrize **without** renaming) are the high-leverage, low-risk structural moves. TYPE-6 is the
PERF-1 dispatch fix. TYPE-4's `AbstractVector` route is the cleanest extensibility lever for v2.0.
TYPE-7 (joint N-channel) is a real but large, deferrable redesign.
