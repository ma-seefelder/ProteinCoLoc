---
phase: 14-decision-and-abstention-layer
reviewed: 2026-08-04T00:00:00Z
depth: standard
files_reviewed: 25
files_reviewed_list:
  - spike/p14/conformal.jl
  - spike/p14/consts.jl
  - spike/p14/decide.jl
  - spike/p14/fdr.jl
  - spike/p14/fuse.jl
  - spike/p14/pools.jl
  - spike/p14/posterior.jl
  - spike/p14/provenance.jl
  - spike/p14/result.jl
  - spike/p14/run_p14_bar_sensitivity.jl
  - spike/p14/run_p14_conformal.jl
  - spike/p14/run_p14_fdr_check.jl
  - spike/p14/run_p14_ood_arm.jl
  - spike/p14/run_p14_real_images.jl
  - spike/p14/run_p14_riskcoverage.jl
  - spike/test/test_p14_conformal.jl
  - spike/test/test_p14_consts.jl
  - spike/test/test_p14_decide.jl
  - spike/test/test_p14_decoupling.jl
  - spike/test/test_p14_fdr.jl
  - spike/test/test_p14_fuse.jl
  - spike/test/test_p14_pools.jl
  - spike/test/test_p14_posterior.jl
  - spike/test/test_p14_provenance.jl
  - spike/test/test_p14_result.jl
findings:
  critical: 2
  warning: 11
  info: 6
  total: 19
status: issues_found
---

# Phase 14: Code Review Report

**Reviewed:** 2026-08-04
**Depth:** standard
**Files Reviewed:** 25
**Status:** issues_found

## Summary

Phase 14 implements the decision/abstention layer: a three-class posterior (`posterior.jl`),
a hand-rolled split-conformal hedge (`conformal.jl`), a Bayesian-FDR prefix rule (`fdr.jl`),
the D-05 asymmetric fusion (`fuse.jl`), the composition point (`decide.jl`), a result type
(`result.jl`), a shared draw/null layer (`pools.jl`), and six runners.

The core numerics are, on inspection, **correct**. The conformal quantile is the
`ceil((n+1)(1-α))`-th order statistic (not `Statistics.quantile`), the LAC set uses `<=` at the
threshold, `sort` is used rather than `sort!`, the FDR rule is a running **mean** (not a running
max) with the prefix property asserted, `sortperm`/`sort` are stable, the degenerate `n == 0`
branch precedes the arithmetic, log-sum-exp subtracts the max before any `exp`, every class read
is by name, and the abstain-then-sort ordering matches the stated Assumption A2. The seed
discipline is genuinely executable (recomputed forbidden inventory, full Philox key-word
cross-product, membership/non-membership asserted at draw time). The tests are largely
non-vacuous and several are deliberately adversarial.

The defects found are concentrated in three places: **(1)** a public return field whose
docstring is provably wrong and whose misuse fails *silently into the honest-looking
`:not_checked` state*; **(2)** guards and closed-vocabulary constants that are decorative rather
than enforcing; and **(3)** a set of small robustness/consistency gaps (memoization that ignores
its own keyword, a Costes key word outside the disjointness proof, two different argmax tie
rules, a Dict-ordered artifact column, `@assert` used as the gate mechanism).

Notably, **nothing was found that alters a reported Phase-14 number as the code currently
runs** — all five runners use the correct `p14_ood_input(ref)` path, so CR-01/CR-02 are latent
traps rather than active corruptions. They are classified BLOCKER because the failure they
enable is silent and produces a plausible-looking result, which is the exact failure class this
phase's own documentation says must be loud.

No recommendation in this report touches `spike/p14/consts.jl`'s frozen values.

---

## Critical Issues

### CR-01: `p14_ood_reference().ood_nulls` is documented as the shape `ood_verdict` reads — it is not, and using it silences the whole batch

**Classification:** BLOCKER
**File:** `spike/p14/pools.jl:553-557` (and the return-shape docstring at `spike/p14/pools.jl:457`)
**Also affects:** `spike/test/test_p14_pools.jl:258-259`

**Issue:**
`p14_ood_reference` returns

```julia
            # The shape `p14_ood_state` and the shipped verdict both read: the fit PLUS the
            # operating point. ...
            ood_nulls = merge(nulls, (thr = thr,)),
```

`nulls` is `fit_ood_nulls(...)`, which returns `(cont = ..., μS = ..., C = ...)`
(`spike/validation/ood.jl:112-119`). The merged tuple therefore carries
`(cont, μS, C, thr)` and **has no `:density` key**.

The shipped verdict reads its fit as `ood_nulls.density` (`src/amortized/ood.jl:368,370`), and
`_p14_ood_verdict` gates on exactly that (`spike/p14/decide.jl:291`). So the comment's claim
that this is "the shape ... the shipped verdict [reads]" is false, and `decide.jl:246-274`'s own
docstring for `p14_ood_input` says as much ("Those two shapes are not the same, and the
difference is SILENT until the first call").

The consequence of believing the comment is not an exception. `_p14_ood_verdict` returns
`nothing`, `p14_decide_one` resolves `:not_checked`, `p14_fuse` abstains under D-06, and **every
item in the batch abstains** — which is exactly the reportable, honest-looking outcome the
runners are written to print without alarm. SC3-d's margin would then be exactly zero and SC1-b's
decided fraction zero, for a reason that has nothing to do with the detector.

`spike/test/test_p14_pools.jl:258-259` asserts `haskey(OOD_REF.ood_nulls, :thr)` and
`OOD_REF.ood_nulls.thr === OOD_REF.thr`. Both pass, because `p14_ood_state` only needs `:thr` —
so the test **endorses the false shape claim without ever exercising it against `ood_verdict`**.

**Fix:** make the field actually be the shape it claims, and pin it in the test.

```julia
# spike/p14/pools.jl -- p14_ood_reference return
            # The shape BOTH `p14_ood_state` (which needs `:thr`) and the shipped
            # `ood_verdict` (which needs `:density`) read. Built through the same adapter
            # `decide.jl` uses, so there is ONE spelling of this shape in the phase.
            ood_nulls = (density = nulls, thr = thr),
```

and in `spike/test/test_p14_pools.jl`:

```julia
        @test haskey(OOD_REF.ood_nulls, :density)      # what `ood_verdict` reads
        @test haskey(OOD_REF.ood_nulls, :thr)          # what `p14_ood_state` reads
        # the shape is EXERCISED, not merely inspected
        @test ood_verdict(OOD_REF.ood_nulls, FIX_A[1].Zs, FIX_A[1].Zc; with_pp = false) isa OODVerdict
```

If the flat shape must be kept for another consumer, then delete the misleading comment and
rename the field (e.g. `ood_nulls_flat`), and make `p14_ood_input` the only sanctioned producer
of the verdict-shaped tuple.

---

### CR-02: `_p14_ood_verdict` cannot distinguish "no null was fitted" from "a null was handed over in the wrong shape"

**Classification:** BLOCKER
**File:** `spike/p14/decide.jl:290-294`

**Issue:**

```julia
function _p14_ood_verdict(ood_nulls, Zs, Zc)
    (ood_nulls isa NamedTuple && haskey(ood_nulls, :density) && ood_nulls.density !== nothing) ||
        return nothing
    return ood_verdict(ood_nulls, Zs, Zc; with_pp = false)
end
```

The `|| return nothing` arm collapses three distinct situations into one:

1. `NamedTuple()` — no null was fitted (the legitimate D-06 `:not_checked` case);
2. a **fitted** null passed in the wrong shape (CR-01, or any future caller error);
3. a non-`NamedTuple` argument entirely (e.g. a raw `fit_ood_nulls` return, a `Dict`).

Only (1) is the honest not-checked state. (2) and (3) are programming errors that this phase's
own doctrine says must fail loudly — `fuse.jl:59-60`: *"The failure this prevents is silent and
scientific; the cost it imposes is loud and operational, which is the right way round"*, and
`result.jl:213-217`: *"a silent abstention is nearly as bad as a wrong call"*. As written, a
shape error is reported as the most defensible-looking answer the layer can produce.

**Fix:** keep the *absent* case silent, make the *malformed* case throw.

```julia
function _p14_ood_verdict(ood_nulls, Zs, Zc)
    # THE ONLY SILENT CASE: nothing was fitted at all. That is D-06's :not_checked.
    (ood_nulls isa NamedTuple && isempty(ood_nulls)) && return nothing
    ood_nulls isa NamedTuple || throw(ArgumentError(
        "_p14_ood_verdict: `ood_nulls` must be a NamedTuple in the shape `p14_ood_input` " *
        "produces; got $(typeof(ood_nulls)). A malformed null is a caller error and must NOT " *
        "be reported as the honest :not_checked state (D-06)."))
    haskey(ood_nulls, :density) || throw(ArgumentError(
        "_p14_ood_verdict: `ood_nulls` carries $(keys(ood_nulls)) but the shipped `ood_verdict` " *
        "reads its fit as `ood_nulls.density` (src/amortized/ood.jl:368). A fitted null in the " *
        "wrong shape would otherwise resolve to :not_checked and silence the entire batch."))
    ood_nulls.density === nothing && return nothing
    return ood_verdict(ood_nulls, Zs, Zc; with_pp = false)
end
```

---

## Warnings

### WR-01: the Pitfall-7 "positive half" guard is satisfied by a trailing comment, not by code

**Classification:** WARNING
**File:** `spike/p14/decide.jl:405`, `spike/test/test_p14_decoupling.jl:589-595`

**Issue:** the guard requires every `costes_p(` call site in `decide.jl` to name `P14_DEV_SEED`
on the same line:

```julia
                occursin("costes_p(", s) || continue
                occursin("P14_DEV_SEED", s) || push!(bad, ...)
```

The only such call site is

```julia
    p_s = costes_p(seed, idx, mci_s)  # `seed` defaults to P14_DEV_SEED; asserted above (Pitfall 7)
```

`_p14_body` / `_P14_TABLE` (`test_p14_decoupling.jl:78-79`, `167-171`) strip only **whole-line**
comments (`startswith(strip(l), "#")`), so the trailing comment survives and is what satisfies
`occursin("P14_DEV_SEED", s)`. The *code* on that line passes a `seed` variable. Delete the
comment and the guard goes red; change the code to pass a forbidden seed and the guard stays
green. The guard is therefore not testing the property it names. (The run-time
`_p14_assert_costes_seed` is the real protection here, so this is a guard defect, not a live
seed leak.)

**Fix:** strip trailing comments too, and/or make the predicate look at the argument.

```julia
# test_p14_decoupling.jl -- strip trailing comments as well as whole-line ones
_strip_inline_comment(l) = first(split(l, '#'))   # no `#` appears inside a Phase-14 string literal
...
    lns = String[_strip_inline_comment(l) for l in split(raw, '\n') if !startswith(strip(l), "#")]
```

and, in `decide.jl`, make the call site self-evidencing:

```julia
    p_s = costes_p(seed, idx, mci_s)   # seed :: default P14_DEV_SEED, checked by _p14_assert_costes_seed
```
becomes
```julia
    @assert seed === P14_DEV_SEED || _p14_assert_costes_seed(seed) === nothing
    p_s = costes_p(seed, idx, mci_s)
```
(or simply name `P14_DEV_SEED` in the default-resolution expression on the same line).

---

### WR-02: `P14_DECIDE_REASONS` is dead — the "closed set" of DECIDE reasons enforces nothing

**Classification:** WARNING
**File:** `spike/p14/fuse.jl:117-118`

**Issue:** `const P14_DECIDE_REASONS = (:decided, :decided_contra_classical)` is declared and
**never read** anywhere in `spike/p14/` or `spike/test/` (the fuse test declares its own local
`FUSE_FIX_DECIDE_REASONS` literal at `test_p14_fuse.jl:129`). `p14_fuse` returns the two symbols
as bare literals at `fuse.jl:270-271`, and `P14Result` validates only the *abstain* reason
(`result.jl:218-230`). A typo in the decide branch (`:decided_contra_classicial`) would travel
into `meta.fuse_reason` and into every artifact with nothing to catch it — the exact asymmetry
`P14_ABSTAIN_REASONS` exists to prevent on the other branch.

**Fix:** either use it or delete it. Using it is cheap:

```julia
# spike/p14/fuse.jl, end of p14_fuse
    r = disagree ? :decided_contra_classical : :decided
    @assert r in P14_DECIDE_REASONS "p14_fuse: $r is not a member of the closed decide-reason set $(P14_DECIDE_REASONS)"
    return (action = :decide, reason = r)
```

and add a test asserting `p14_fuse(...).reason in P14_DECIDE_REASONS` for the decide rows of the
existing truth table.

---

### WR-03: `P14ClassPosterior` is dead, and the docstring claiming the return "IS" one is unchecked

**Classification:** WARNING
**File:** `spike/p14/posterior.jl:87` (declaration), `spike/p14/posterior.jl:202-203` (claim)

**Issue:** `const P14ClassPosterior = NamedTuple{P14_CLASS_KEYS, Tuple{Float64,Float64,Float64}}`
is declared inside the guard block and never referenced again — no annotation, no `convert`, no
test. The comment at `:202` asserts in prose *"this IS a `P14ClassPosterior` without a positional
construction"* while nothing verifies it; if any input were `Float32`, the returned tuple would
have a different type and the claim would silently be false.

**Fix:** make the type the return annotation, which is one token and makes the claim executable.

```julia
function p14_class_posterior(logbf::NamedTuple, prior::NamedTuple)::P14ClassPosterior
```

(or delete the constant if the type is genuinely not wanted).

---

### WR-04: `p14_tau()` and `p14_net_meta()` memoize on nothing, so their own keyword is silently ignored after the first call

**Classification:** WARNING
**File:** `spike/p14/pools.jl:134-137`, `spike/p14/pools.jl:148-154`

**Issue:**

```julia
function p14_net_meta(; net_path = joinpath(P14_P13_DIR, "three_way_net.jld2"))
    if _P14_META_CACHE[] === nothing
        isfile(net_path) || error(...)
        _P14_META_CACHE[] = load_three_way(net_path).meta
    end
    return _P14_META_CACHE[]
end
```

The cache is keyed on nothing. A second call with a *different* `net_path` returns the first
net's metadata, and neither the existence check nor the load runs. `p14_pi_class_masses()` and
`p14_assert_unstratified` sit downstream of that, so a test or runner that pointed at an
alternative net would silently be validated against the wrong measured class prior. The same
applies to `_P14_TAU_CACHE`.

**Fix:** key the cache on the path (or refuse a non-default path when the cache is warm).

```julia
const _P14_META_CACHE = Ref{Any}(nothing)
const _P14_META_PATH  = Ref{String}("")

function p14_net_meta(; net_path = joinpath(P14_P13_DIR, "three_way_net.jld2"))
    if _P14_META_CACHE[] === nothing || _P14_META_PATH[] != String(net_path)
        isfile(net_path) || error("p14_net_meta: no three-way net artifact at $net_path")
        _P14_META_CACHE[] = load_three_way(net_path).meta
        _P14_META_PATH[]  = String(net_path)
    end
    return _P14_META_CACHE[]
end
```

---

### WR-05: every reported gate and every D-01 breach check is an `@assert`, which Julia documents as removable

**Classification:** WARNING
**Files:** `spike/p14/run_p14_conformal.jl:332-334,615`; `spike/p14/run_p14_fdr_check.jl:599-601,923`;
`spike/p14/run_p14_riskcoverage.jl:600-602,947,965,979`; `spike/p14/run_p14_ood_arm.jl:468-470,883`;
`spike/p14/provenance.jl:344-345`

**Issue:** the Julia manual states of `@assert`: *"An assert might be disabled at various
optimization levels. Assert should therefore only be used as a debugging tool"*. Every SC1-b /
SC1-d / SC3-a/b/c / SC3-d verdict, and every "D-01 decoupling breach" / "D-02 decoupling breach"
check inside the runners, is implemented with `@assert`. The `error(...)` construct is used
correctly elsewhere in the same files for hard preconditions (`_p14_load_report`,
`p14_ood_reference`'s hard stop), so the inconsistency is internal.

The artifacts do remain honest — `sc3a_met`, `sc1b_met`, `iteration_trigger_fired` etc. are
persisted before the assertion — so the blast radius is limited to "the process exits 0 on a
failed gate". That is still the wrong exit code for a pre-registered gate.

**Fix:** use `error` for the verdicts and the breach checks; keep `@assert` for genuine internal
invariants.

```julia
    if reported && !(coverage >= band_lower)
        error("""
        SC1-d NOT MET. ...
        """)
    end
```

---

### WR-06: the Costes Philox key word is outside the disjointness cross-product `consts.jl` proves

**Classification:** WARNING
**File:** `spike/p14/decide.jl:345-352,405`; `spike/p14/consts.jl:168-184`

**Issue:** `consts.jl` states the standard explicitly — *"THE KEY WORDS THEMSELVES, not merely
the seeds ... `seed xor salt xor counter` is what actually keys Philox in this project's
runners, so THAT is the quantity that must be proven distinct"* — and `_p14_key_words()`
enumerates only `P14_{DEV,FIX}_SEED ⊻ P14_SALT [⊻ counter]`.

The one random stream `decide_coloc` actually consumes on real images is the Costes null, keyed
`Philox4x((master_seed ⊻ COSTES_SALT, idx))` (`spike/comparator/classical.jl:135-136`, with
`COSTES_SALT = 0xA5A5A5A5DEADBEEF`, `config.jl:39`). That key word,
`P14_DEV_SEED ⊻ COSTES_SALT`, is **not** a member of `_p14_key_words()` and is therefore covered
by no assertion. `_p14_assert_costes_seed` checks only the *seed value* against
`_p13_forbidden()`, not the resulting key word.

There is no live collision today (`P14_DEV_SEED ≠ MASTER_SEED`, so the two Costes key words
differ), but the stated proof does not cover the stream it most needs to.

**Fix:** extend the run-time guard to the key word rather than the seed alone.

```julia
function _p14_assert_costes_seed(seed)
    UInt64(seed) in _p13_forbidden() && throw(ArgumentError(...))
    kw = UInt64(seed) ⊻ COSTES_SALT
    kw in Set(_p14_p13_key_words()) && throw(ArgumentError(
        "p14_cross_method: the Costes key word 0x$(string(kw, base = 16, pad = 16)) collides " *
        "with a Phase-13 Philox key word; the null would ride draws a Phase-13 number already rode."))
    return nothing
end
```

(Adding the Costes key words to `_p14_key_words()` in `consts.jl` is *not* proposed — that file
is frozen and may only be appended to; the guard belongs at the call site.)

---

### WR-07: `_p14_rc_curve` can return `thresholds` longer than `coverage`/`selective_risk`, and the three are persisted as parallel columns

**Classification:** WARNING
**File:** `spike/p14/run_p14_riskcoverage.jl:282-303`, persisted at `:820-822`

**Issue:**

```julia
        sel == 0 && continue
        push!(cov, sel / n)
        push!(risk, nloss / sel)
    end
    return (coverage = cov, selective_risk = risk, thresholds = collect(Float64, thr))
```

`thresholds` is built from *all* unique kappa values while `coverage`/`selective_risk` skip any
threshold that selects nothing. With a finite kappa the skip is unreachable (`t` is always an
attained value), but a single `NaN` in `kappa` sorts to the end of `unique`, makes
`kappa[i] >= NaN` false for every `i`, and desynchronises the three vectors — which are then
written to the artifact as `curve_coverage` / `curve_selective_risk` / `curve_thresholds` with
no length relation asserted. Nothing upstream validates that the pool's posterior columns are
finite (`_p14_rc_load_pool` checks presence and counts only).

**Fix:** assert the invariant rather than relying on it, and validate the input.

```julia
    kept = Float64[]
    for t in thr
        ...
        sel == 0 && continue
        push!(kept, Float64(t)); push!(cov, sel / n); push!(risk, nloss / sel)
    end
    @assert length(kept) == length(cov) == length(risk) "_p14_rc_curve: the three curve columns desynchronised"
    return (coverage = cov, selective_risk = risk, thresholds = kept)
```

plus, in `_p14_rc_load_pool`, `@assert all(isfinite, ...)` on the three posterior columns.

---

### WR-08: `P14Result`'s constructor validates every field except `null_split`

**Classification:** WARNING
**File:** `spike/p14/result.jl:195-274` (constructor), field declared at `:184`

**Issue:** the constructor is described as the place *"where every invariant of the type is
enforced"* (`decide.jl:996-999`) and it does validate the decision, the OOD state, the
abstain-reason biconditional, the class-posterior key set and sum, the composite null against
`1 - p.coloc`, the conformal tuple, and `cross_method.basis`. It performs **no check at all** on
`null_split`, even though `p_random + p_exclusion == v` is precisely the identity
`_p14_fdr_partition_sanity` asserts downstream on the persisted columns
(`run_p14_fdr_check.jl:379-382`). A caller can construct a legal-looking `P14Result` whose
`null_split` describes a different item.

**Fix:**

```julia
        Set(keys(null_split)) == Set((:p_random, :p_exclusion, :v)) || throw(ArgumentError(
            "P14Result: `null_split` must carry exactly (:p_random, :p_exclusion, :v); got $(keys(null_split))."))
        abs((null_split.p_random + null_split.p_exclusion) - v) < 1e-12 || throw(ArgumentError(
            "P14Result: null_split.p_random + null_split.p_exclusion must equal the composite " *
            "null within 1e-12; got $(null_split.p_random + null_split.p_exclusion) against $v. " *
            "The composite null is additive only over a genuine partition."))
```

---

### WR-09: `abstain_reason_counts` is built by iterating a `Dict`, so the persisted key order is not deterministic

**Classification:** WARNING
**File:** `spike/p14/run_p14_fdr_check.jl:673-677, 808`

**Issue:**

```julia
    reasons = Dict{Symbol, Int}()
    ...
        abstain_reason_counts = NamedTuple{Tuple(keys(reasons))}(Tuple(values(reasons))),
```

`keys`/`values` on a `Dict` iterate in slot order, which depends on hash placement and insertion
history. The *pairing* is correct, but the resulting NamedTuple's **field order** — and hence the
shape of the persisted artifact column and the `$reasons` printout at `:679` and `:912` — is not
a stable function of the data. `decide.jl:764-772` already solves this the right way, by
enumerating the frozen `p14_abstain_reasons()` order and filtering; two spellings of one ledger
now exist.

**Fix:** reuse the frozen enumeration.

```julia
        abstain_reason_counts = let ks = Tuple(k for k in p14_abstain_reasons() if haskey(reasons, k))
            NamedTuple{ks}(Tuple(reasons[k] for k in ks))
        end,
```

---

### WR-10: two different argmax tie rules coexist, so SC3-c's `hard` and the layer's own call can disagree

**Classification:** WARNING
**Files:** `spike/p14/decide.jl:321-326`, `spike/p14/run_p14_riskcoverage.jl:254-265`

**Issue:** `_p14_amortized_call` requires a **strict** majority and falls to `RANDOM` on any tie
(*"it errs toward the null, which is the right direction for a rule whose false discoveries are
the controlled quantity"*). `_p14_rc_argmax` breaks ties toward the **first** `P14_CLASS_KEYS`
entry, i.e. `:coloc` (*"Ties resolve in the frozen `P14_CLASS_KEYS` order"*). SC3-c's `loss` /
`hard` vector is built from the second rule while the decision the layer actually emitted came
from the first. Both docstrings argue ties are measure-zero, and they are — but the phase then
has two functions named "the argmax, by name" that are not the same function, and the risk-
coverage curve is claimed to describe *"ONE ordering"* with the decision layer.

**Fix:** have `_p14_rc_argmax` delegate, so there is one rule:

```julia
_p14_rc_argmax(p::NamedTuple) = _p14_class_key(_p14_amortized_call(p))
```

(or, if the runner must stay independent of `decide.jl`, copy the strict-majority/`:random`
tie-break verbatim and say so).

---

### WR-11: `ceil(Int, (n + 1) * (1 - alpha))` is evaluated in floating point and can throw a spurious infeasibility error

**Classification:** WARNING
**File:** `spike/p14/conformal.jl:222-228`

**Issue:** the index is computed as `k = ceil(Int, (n + 1) * (1 - alpha))` on `Float64`. When
`(n+1)(1-α)` is mathematically an exact integer but the product rounds *up* by one ulp, `ceil`
returns `k+1`. The direction is always conservative for coverage (a larger `qhat` admits more
classes, so the `>= 1-α` bound is never broken), **but** it can make `k == n + 1` and trip the
infeasibility branch, throwing `"alpha = ... is infeasible at n = ..."` for an alpha that is in
fact exactly feasible. That is a hard failure on a legitimate input.

The frozen operating point (`n = 2000`, `α = 0.10` → `1800.9…` → `1801`) and the fixtures
(`n = 9`, `n = 19`) are unaffected; the hazard is at the boundary `α = 1/(n+1)`, which
`test_p14_conformal.jl:173-175` approaches deliberately.

**Fix:** compute the index without a floating-point product straddling an integer.

```julia
    # k = ceil((n+1)(1-alpha)) computed from the complement, so an exact-integer level cannot
    # round UP past its own integer and trip the feasibility branch.
    k = (n + 1) - floor(Int, (n + 1) * float(alpha))
    k = min(k, n + 1)   # keep the infeasibility branch reachable for a genuinely infeasible alpha
```

(and keep the existing `k <= n ||` throw unchanged).

---

## Info

### IN-01: `coverage_pin_ok = true` is persisted as a literal, not as the computed comparison

**Classification:** CONVENTION
**File:** `spike/p14/run_p14_bar_sensitivity.jl:522`
**Deviation:** the artifact records a verification outcome as a hard-coded `true` rather than
the value of the check.
**Convention it violates:** every other verdict field in this lane is the computed predicate
(`sc3a_met = ... >= P14_SKILL_FLOOR`, `coverage_meets_band = !iteration_trigger_fired`,
`readonly_digest_equal = digest_before == digest_after`).
**Suggested fix (recommend, non-blocking):** compute it, and combine with WR-05 so it does not
depend on `@assert` having run.

```julia
    pin_ok = abs(pin.spearman - spearman) < P14_BS_PIN_TOL &&
             pin.n_points == Int(rc["spearman_n_points"])
    pin_ok || error("run_p14_bar_sensitivity: ...")
    ...
        coverage_pin_ok = pin_ok,
```

### IN-02: the tau literal `0.15` is copied into two Phase-14 test files, and the literal guard only catches assignments

**Classification:** CONVENTION
**File:** `spike/test/test_p14_provenance.jl:78,89`; `spike/test/test_p14_decide.jl:129`
**Deviation:** D-07's stated rule is that Phase 14 *"never copies the literal"*; the guard is
`P14_TAU_LITERAL_RE = r"\w*[Tt][Aa][Uu]\w*\s*=\s*0\.15"`, so `@test P14_PROV.tau === 0.15` and
`@test DEC_B.tau === 0.15` both pass (three `=` characters do not match `\s*=\s*`), and so would
a binding with no `tau` in its name, e.g. `const P14_CUT = 0.15`.
**Convention it violates:** "the value is read from the one place it was measured, never
mirrored into this namespace" (`consts.jl:398-404`).
**Suggested fix (recommend, non-blocking):** assert against the inherited constant instead of
the literal — `@test P14_PROV.tau === P13_TAU` (already in scope via the p13 include, and
already asserted at `test_p14_provenance.jl:129`) — and widen the regex to catch any binding of
the value, e.g. `r"^\s*(const\s+)?\w+\s*=\s*0\.15\b"`.

### IN-03: the SC2-d guard is per-line and does not see a two-line `patch_correlation` → `tau` comparison

**Classification:** CONVENTION
**File:** `spike/test/test_p14_decoupling.jl:551-554`
**Deviation:** the predicate requires `patch_correlation(`, a `tau` mention and no `ghat(` on the
**same** line. `mu = patch_correlation(x)` followed by `three_way_label(mu, ...; tau = tau)` on
the next line is not caught, which is the more idiomatic Julia spelling of the very bug.
**Convention it violates:** the file's own standard that a guard must be circumventable only by
a deliberate act, not by ordinary formatting.
**Suggested fix (recommend, non-blocking):** track the local bound to `patch_correlation(` and
flag a later `tau` comparison against that name within the same function body; or, more cheaply,
also flag any line matching `patch_correlation(` whose result is assigned to a name that is later
compared against `tau` without an intervening `ghat(`.

### IN-04: `:disagreement_and_ood` names a state that is unreachable for the reason its name implies

**Classification:** CONVENTION
**File:** `spike/p14/fuse.jl:114-115, 268-269`
**Deviation:** branch 5 is reached only when `ood_state === :not_checked` under
`allow_unchecked_ood = true` (`:fired` and un-opted-out `:not_checked` have already returned), so
the reason recorded is `:disagreement_and_ood` for a case in which the detector **never ran**.
The name reads as "the detector fired *and* the classics disagreed".
**Convention it violates:** D-06's own rule, restated three times in this phase, that
`ood_fired` and `ood_not_checked` must never be summed into an "OOD" total
(`run_p14_ood_arm.jl:383-393`).
**Suggested fix (recommend, non-blocking):** the reason set is a closed vocabulary carried in
artifacts, so renaming it is a schema change — instead, document the reachability explicitly in
`p14_abstain_reasons()`'s docstring and in the `_p14_ood_breakdown` table so a reader of the
artifact cannot read this column as detector performance.

### IN-05: five byte-identical `_p14_*_save_report` helpers and two identical cross-tab helpers

**Classification:** CONVENTION
**Files:** `run_p14_conformal.jl:249-261`, `run_p14_fdr_check.jl:206-218`,
`run_p14_riskcoverage.jl:198-210`, `run_p14_ood_arm.jl:164-176`,
`run_p14_real_images.jl:250-262`; and `run_p14_conformal.jl:275-291` vs
`run_p14_ood_arm.jl:261-277`.
**Deviation:** seven copies of two helpers, each with a comment explaining that importing from a
sibling runner would inherit its `main`.
**Convention it violates:** the phase's own repeatedly stated rule that *"a local
re-implementation is how two versions of one statistic silently diverge"*
(`pools.jl`, `posterior.jl`, `test_p14_decoupling.jl:_NO_REIMPLEMENT`).
**Suggested fix (recommend, non-blocking):** the stated obstacle is real but avoidable — put the
two helpers in a small `spike/p14/artifact_io.jl` behind the standard guarded include
(`isdefined(@__MODULE__, :p14_save_report) || include(...)`), which defines no `main` and so
cannot be inherited.

### IN-06: `_p14_headline_line`'s `n_total == 0` branch is unreachable from `p14_headline`

**Classification:** CONVENTION
**File:** `spike/p14/result.jl:801`
**Deviation:** `frac = n_total == 0 ? NaN : n_decided / n_total`, but `P14BatchDecision`'s
constructor rejects `n_total <= 0` (`result.jl:379-381`), so the guard is dead on the only real
call path; it exists for `p14_headline_template_probe()`, which passes `100`.
**Convention it violates:** dead-branch avoidance; the file's own note that a check which cannot
fire is not a check.
**Suggested fix (recommend, non-blocking):** either drop the branch, or `@assert n_total > 0` in
`_p14_headline_line` so the precondition is stated once where the division happens.

---

## Verified-correct (no finding)

Recorded so a later reader knows these were checked rather than skipped:

- `p14_conformal_quantile` uses the order statistic, not `Statistics.quantile`; `sort` (not
  `sort!`) so the caller's calibration vector is not reordered; the `n = 0` and non-finite-score
  paths throw rather than return.
- `p14_conformal_set` uses `<=` at the threshold (matching the quantile's definition) and builds
  the set by iterating the frozen key tuple, so the emitted order is independent of the caller's
  declaration order.
- `p14_bayes_fdr`: running **mean** (not running max); `n == 0` handled before `cumsum`/
  `findlast`; `accepted` holds original indices; `sortperm` is stable; the prefix property is
  asserted rather than assumed; `t_star`/`cost_ratio` are `NaN` (not `0`) when nothing is
  accepted.
- `p14_class_posterior`: log-sum-exp with the max subtracted before any `exp`; three-term `max`
  by name; the D-08 structural zero checked with `===`; non-positive prior masses rejected before
  `log`.
- `_p14_rule_at_prior` is the single code path for both the headline run and every sensitivity
  rung, and `decide_coloc` asserts the per-pair and batch fusions agree item by item.
- The decided-subset → batch index mapping is explicit and asserted on both sides
  (`decide.jl:948-950`, `run_p14_fdr_check.jl:696-698`).
- `OOD_FAMILIES` is a NamedTuple, so `collect(keys(...))` in `run_p14_ood_arm.jl` is
  order-deterministic; the arm tables and `weakest_family` are therefore reproducible.
- `p14_misspec_pool`'s item-1 summariser agreement check is genuinely non-vacuous: both sides
  construct their RNG from the same index and `p13_simulate_labelled` composes exactly the same
  four calls `_p14_summarize` does.
- `p14_assert_unstratified`'s two checks are not redundant: the balance check requires **all
  three** masses within `tol_sd` SE of 1/3 simultaneously, which the measured prior does not
  satisfy.
- `test_p14_fdr.jl` and `test_p14_conformal.jl` both carry the directed discriminating fixtures
  they claim (running-max vs running-mean; order statistic vs interpolated quantile), and neither
  is satisfiable by the wrong implementation.

---

_Reviewed: 2026-08-04_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
