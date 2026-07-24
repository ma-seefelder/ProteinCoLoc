---
spike: 003
name: bayes-model-modularity
type: analysis
validates: "Given the inline @model + ADVI inference flow + post-processing, when audited, then modularity defects and (active/latent) correctness bugs are identified and a runnable, testable redesign is proposed"
verdict: VALIDATED
related: [001, 002, 005]
tags: [bayes, modularity, correctness]
---

# Spike 003: Bayesian model & inference modularity

## What This Validates

Given the Turing model and the ADVI inference + post-processing flow, when audited, then the
separation-of-concerns defects and correctness issues are identified and a modular, runnable,
testable redesign is proposed — including reuse of the prior ranges for the v2.0 SBI pipeline.

## Research / Method

Read `bayes.jl` (`@model`, `colocalization()`, `convert_posterior_samples`, `compute_BayesFactor`,
`CoLocResult`) and `colocalization.jl` (correlation gate). Verified the resolved/installed Turing
stack against `Manifest.toml` and the package cache. **Important correction:** the `@model`,
`colocalization()`, the `sample/vi/rand` calls, and `CoLocResult` construction all live in
**`bayes.jl`**, not `colocalization.jl` (which ends at line 240); the analysis brief mislabeled
this and the per-finding citations are corrected below.

## Findings (all verified)

| ID | Finding | Verdict | Sev | Bug? | Effort |
|----|---------|---------|-----|------|--------|
| MODEL-1 | Turing `@model` defined **inline** inside `colocalization()` — not testable/reusable | ⚠ needs-nuance | High | no | M |
| INFER-1 | Posterior path uses **removed** `ADVI(n,iter)`/`vi(m,alg)`/`DynamicPPL.syms` API | ⚠ needs-nuance | High | **yes** | M |
| POST-3 | `generate_txt` Δρ mixes posterior μ_sample with **prior** μ_control | ✓ confirmed | High | **yes** | S |
| POST-1 | `convert_posterior_samples` magic `1:10` mislabels 2 local params as globals | ⚠ needs-nuance | Med | **yes** | S |
| POST-2 | Fisher-z docstring claims a transform that is never applied (dead `#tanh`) | ⚠ needs-nuance | Med | no* | S |
| INFER-2 | `num_latent = ncols-3` fragile magic, misused as ADVI's sample arg | ⚠ needs-nuance | Med | no | S |
| DRY-1 | BF odds / KDE-integration logic duplicated in `compute_BayesFactor` & `bayes_rangeplot` | ✓ confirmed | Med | no | S |
| MAGIC-1 | Hardcoded `length(a) <= 15` patch-inclusion threshold | ✓ confirmed | Low | no | S |
| INFER-3 | Prior drawn as full Chains over all locals, then 99% discarded | ⚠ needs-nuance | Low | no | M |

\* POST-2 reclassified by the verifier from bug → stale-doc / dead-code trap (no wrong runtime output today).

### MODEL-1 — lift the `@model` to module scope — ⚠ needs-nuance (high)
`bayes.jl:267-303` (`@model`), `bayes.jl:234-339` (`colocalization()`), instantiation `bayes.jl:306`.
The model is nested inside `colocalization()`, so it can't be unit-tested, used for
prior-predictive checks, or have its prior ranges reused by the v2.0 SBI simulator (CLAUDE.md
requires the simulator prior stay consistent with these ranges, `bayes.jl:274-291`). **Verified the
`@model` body captures no enclosing variable**, so lifting it to a top-level `@model coloc_model`
is behavior-preserving. Also relocate the **orphan docstring** at `bayes.jl:250-266` (a no-op
string literal sitting inside the function) onto the lifted model. Split `fit_prior`/`fit_posterior`
helpers so definition / inference / packaging are separable.

### INFER-1 — posterior path uses a removed API (may not run) — ⚠ needs-nuance (high, BUG)
`bayes.jl:323-329`, `bayes.jl:44`. Written against pre-0.35 Turing/AdvancedVI:
- `ADVI(num_latent, iter)` → `ADVI` is now `const ADVI = KLMinRepGradDescent` whose sole positional
  arg is an AD backend, not two ints → **MethodError**.
- `vi(m, alg)` → `vi` now needs `vi(model, q, max_iter; …)`; the 2-arg form fails.
- `DynamicPPL.syms(VarInfo(m))` → no longer public in DynamicPPL 0.40/0.41 → **UndefVarError**.

**Critical caveat (verifier):** the proposed fix `q,_info,_state = vi(...)` only works on Turing
**0.43.x**; on **0.44.5/0.45.0** (installed) `vi` returns a non-iterable `VIResult` struct, so use
`res = vi(m, q0, iter; adtype=AutoForwardDiff()); q = res.q`. Also resolve the
**Manifest/cache mismatch** (pinned Turing 0.42.8 / DynamicPPL 0.39.13 are not even installed;
0.43–0.45 / 0.40–0.41 are) before porting. This is the single biggest *runnability* risk in `src/`.

### POST-3 — Δρ prior/posterior mix-up — ✓ confirmed (high, BUG)
`utils.jl:368` computes `posterior.μ_sample .- prior.μ_control` — everywhere else Δρ is
within-distribution (`bayes.jl:110`, `plot.jl:527,597,680`). The `result.txt` Δρ mean/median/CI are
wrong and inconsistent with the Bayes factor in the same row. Fix: use `posterior_samples.:μ_control`
(or the `delta_rho` helper from VIZ-1). **The same bug was independently found by 4 of 5 dimensions.**
Helper-name caveat: `Δρ` collides with locals in `plot.jl` — name the helper `delta_rho`.

### POST-1 — `convert_posterior_samples` magic `1:10` — ⚠ needs-nuance (med, BUG)
`bayes.jl:42-50`. The model has 8 scalar globals then 6 `filldist` groups; `rand(q,N)` returns rows
in declaration order, so rows 9–10 are `μ_control_image[1..2]`, but `syms[1:10]` labels them
`:μ_control_image`/`:ν_control_image`. Columns 9–10 carry **per-image local values under wrong
global-sounding names**. Fix: select the 8 globals explicitly (`1:8` + named vector). **Correction:**
these mislabeled columns **do escape** to the user-facing `posterior_samples` CSV (`main.jl:176-179`)
— so the posterior CSV has 10 cols while the prior CSV has 8. The `1:8` fix changes the posterior
CSV column count 10→8 (call this out in the PR); no in-code consumer is affected.

### POST-2 — Fisher-z dead code + false docstring — ⚠ needs-nuance (med, not a live bug)
`bayes.jl:40,52-53`. The docstring says "Fisher z … undone by tanh," but the only such line is
commented out and **no `atanh`/`tanh` is applied anywhere** — the model lives in raw `[-1,1]`
correlation space (priors truncated to `[-1,1]`). Current output is already correct, so reclassify
`is_correctness_bug=false`: it's a stale docstring + a **trap** (uncommenting `tanh.(samples)` would
corrupt the unbounded ν and positive σ/τ columns). Fix: delete the dead line, correct the docstring.
(The "implicit return" sub-claim is wrong — the trailing assignment already returns the DataFrame.)

### INFER-2 / DRY-1 / MAGIC-1 / INFER-3
- **INFER-2** (`bayes.jl:323,325`): `num_latent = DataFrame(prior_chain) ncols - 3` is an
  undocumented, version-fragile count, and the legacy `ADVI`'s first arg is the MC sample count,
  **not** the latent dim — a misnomer/misconfiguration. Folds into INFER-1; prefer
  `length(DynamicPPL.VarInfo(m)[:])` if a count is still needed.
- **DRY-1** (`bayes.jl:109-137`, `plot.jl:680-699`): the tail-prob-via-KDE/`quadgk` + odds→BF math
  is copy-pasted and has **already diverged** (rangeplot drops the `ε>1e-5` warning, adds a NaN
  guard). Extract `_tail_prob`/`_bayes_factor`; `_tail_prob` must return `(1-p, ϵ)` so the warning
  stays at the `compute_BayesFactor` site. (Overlaps VIZ-4.)
- **MAGIC-1** (`colocalization.jl:235`): promote `length(a) <= 15` to `min_pairs::Int = 15` (default
  preserves behavior). It *is* documented in the docstring but not in the signature / not tunable.
- **INFER-3** (`bayes.jl:308-318,323`): `sample(m, Prior(), 100_000)` materializes a Chains over all
  locals, kept in `advi_result` (dead — never read), to extract only the 8 global marginals. Draw
  the 8 globals directly with `rand(dist,N)` (sharing the prior specs with the lifted model).

## Signal for the Build
**INFER-1 first** — confirm the code even runs on the resolved stack before any other bayes work
(resolve the Manifest mismatch, port to the modern `vi`). Then MODEL-1 (lift the model; unlocks
testing **and** v2.0 prior reuse) and the bug fixes POST-3/POST-1. DRY-1 + MAGIC-1 + INFER-2/3 are
clean follow-ups. All of this touches `src/` → schedule outside the decoupled spike.
