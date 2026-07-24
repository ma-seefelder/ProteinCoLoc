# Reference: `src/` Refactoring Analysis (spikes 001–005)

**Scope:** A *separate, earlier* spike session — an **analysis-only** refactoring roadmap for the
ProteinCoLoc v1.0 `src/` core. Nothing here was implemented. Every proposal is behavior-preserving
unless explicitly flagged as a correctness-bug fix. 5 dimensions, every finding adversarially
verified by an independent skeptic (46 agents total); hot-path claims micro-benchmarked.

> **This section is intentionally brief — it is an index, not the detail.** The consolidated,
> prioritized, effort-tagged roadmap lives in `.planning/spikes/PROPOSAL.md`. Per-finding evidence
> (stable IDs `PERF-n`/`TYPE-n`/`MODEL-n`/`INFER-n`/`POST-n`/`VIZ-n`/`PROV-n`/`OUT-n`) lives in the
> copied `sources/00N-*/README.md`. Read those before touching `src/`.

## Constraints (read before building)

- **Hard gate — do NOT touch `src/` during the active v2.0 spike.** CLAUDE.md forbids it. These
  refactors are a *separate effort*: sequence them **before** resuming the spike or **during the
  v2.0 productionization** phase. Do not interleave.
- **Stack reality:** the resolved/installed Turing stack (0.43–0.45 / DynamicPPL 0.40–0.41) differs
  from the legacy `ADVI(n,iter)` / `vi(m,alg)` / `DynamicPPL.syms` API the posterior path is written
  against. The pinned Manifest (Turing 0.42.8) is *not even installed*. **Confirm the posterior path
  runs at all (INFER-1 / B2) before any other bayes work** — it is the biggest runnability risk.
- **File-map correction (trust the verifiers, not the brief):** the analysis brief mislabeled
  `bayes.jl` code as `colocalization.jl`. `colocalization.jl` (≤ line 240) contains **only**
  `patch`/`unpatch`/`_exclude_zero`/`correlation`. The `@model`, `colocalization()`,
  `convert_posterior_samples`, `compute_BayesFactor`, `_prepare_data`, and `CoLocResult`
  construction all live in **`bayes.jl`**.

## Correctness bugs (Tier 0 — fix first, independent of any refactor)

The single most important finding of the whole session, found independently by **4 of 5**
dimensions:

- **B1 / Δρ bug (`utils.jl:368`)** — the reported headline effect size subtracts a **prior** draw
  of `μ_control` from a **posterior** draw of `μ_sample` (two independent distributions). The
  `result.txt` Δρ mean/median/CI are wrong and inconsistent with the Bayes factor in the same row.
  1-line fix: `posterior_samples.:μ_sample .- posterior_samples.:μ_control`. Regenerate any cached
  result files. Triangulation across 4 dimensions = highest confidence.

Other confirmed bugs: **B2** removed-Turing-API (may not run), **B3** `convert_posterior_samples`
magic `1:10` mislabels 2 local params (posterior CSV has 10 cols vs prior's 8), **B4**
`minmax_norm!` mutates shared image data so inference silently depends on plot flags, **B5** false
Fisher-z docstring / dead `#tanh` trap (not a live bug). Detail + fixes in `PROPOSAL.md` Tier 0.

## The five verified refactoring blueprints

1. **Hot-path type-stability via function barriers (spike 001 / PERF-1 ≡ TYPE-6).** The measured
   win. `correlation()` (`colocalization.jl:221-239`) rebuilds a `Dict` of functions per call
   (→ dynamic dispatch on an inferred `::Function`) and does `vec(collect(view(...)))` Union
   materializations per patch pair. **Fix: resolve the function once, pass it through a barrier**
   `_corr(f::F, x, y) where F` that filters zero/NaN/missing in a single pass into two reused
   `Float64` buffers. **Measured 3.4–6.0× faster, 135–1394× less memory** (`bench_results.txt`).
   A bare `Val(runtime_symbol)` is **not** a fix (itself type-unstable). Also **PERF-5**: hoist
   `InterpKDE` out of the `quadgk` integrand (`bayes.jl:113-124`) — it currently rebuilds a
   2048-point spline on every integration node. Both numerically identical.

2. **Real type hierarchy (spike 002).** The image accessors are a *fake* interface: six methods on
   the abstract type each do `hasfield(...) || throw(...); return img.data`, hard-coding concrete
   field names into the supertype (`MultiChannelImage` implements none, free-rides). And the
   interface is bypassed by direct `img.data`/`img.:name` access everywhere. **Fix: erroring generic
   fallbacks + concrete impls on the subtype (keep `throw(ArgumentError)`); route every field touch
   through accessors.** `MultiChannelImageStack{T} <: AbstractVector{T}` (define `size`+`getindex`,
   delete the hand-rolled `length`/`getindex`/`iterate`) restores the full container protocol for
   free — the cleanest extensibility lever for v2.0 synthetic stacks. Parametrize
   `CoLocResult{St,A} <: AbstractCoLocResult` **without renaming fields** (a `posterior`→`samples`
   rename touches ~20 sites and breaks `typeof(x)==CoLocResult` tests). Drop unused `mutable` on the
   stack.

3. **`bayes` `@model` modularity (spike 003).** The Turing `@model` is defined **inline inside
   `colocalization()`** (`bayes.jl:267-303`), so it cannot be unit-tested, used for
   prior-predictive checks, or have its prior ranges reused by the v2.0 SBI simulator. **Verified the
   body captures no enclosing variable → lifting to a top-level `@model coloc_model` is
   behavior-preserving.** This directly *enables v2.0* (CLAUDE.md requires the simulator prior stay
   consistent with these μ/ν/σ/τ ranges). Split `fit_prior`/`fit_posterior`; relocate the orphan
   docstring onto the lifted model.

4. **Headless viz reuse (spike 004).** Every figure is built/saved with **GLMakie** (needs a live
   OpenGL context — dead on headless servers/CI/GUI-batch); the advertised `.svg` export throws.
   **Fix: render file output with CairoMakie** (PNG+SVG+PDF, no GL context; already the `spike/`
   pattern), change `endswith||@error` guards to `throw(ArgumentError)`. DRY: 7–8 hand-coded Δρ
   copies → one `delta_rho(::CoLocResult)` helper (kills the B1 bug class structurally); shared
   `bf_at` for the duplicated-and-already-diverged BF math; **scoped `with_theme(COLOC_THEME) do…end`
   — never a global `set_theme!`** (it would bleed into the black-background image plots). Replace 5
   bare `catch` blocks with `@warn … exception=(e, catch_backtrace())`.

5. **Provenance by construction (spike 005).** No RNG/seed is threaded through any stochastic step
   → **runs are non-reproducible** (a published Bayes factor cannot be reproduced). No run-manifest:
   zero versions/seed/timestamp/git-commit/input-hashes in any output. **Fix, in dependency order:**
   R1 thread `seed=nothing` → one `rng = Xoshiro(seed)` into `sample`/`rand`/shuffles (`Random.seed!`
   before `vi`, which takes no rng); R2 `write_manifest()` → `run_manifest.toml` (VERSION, pkg
   versions, seed, timestamp, git commit, args, per-input SHA-256, stdlib TOML/SHA/Dates/Pkg); R3
   embed `provenance::NamedTuple` in `CoLocResult` via a defaulted constructor. Also harden the
   hand-rolled `|`-delimited writer (drift-prone parallel header/value lists) and the append-mode
   `result.txt` contamination.

## What to avoid

- Don't implement any of this inside the decoupled v2.0 spike — it edits `src/`.
- Don't use `Val(runtime_symbol)` as the "type-stability fix" — it is itself unstable; use a
  function barrier that resolves the concrete function once.
- Don't rename struct fields as part of parametrizing types — stage behind property aliases; a
  rename breaks `typeof(x)==CoLocResult` tests (`test/runtests.jl:203-298`).
- Don't use a global `set_theme!` for the Bayesian plots — scope it, or you corrupt the image plots.
- Don't trust the analysis brief's file map — several `colocalization.jl` citations are actually
  `bayes.jl` (the verifiers corrected them).

## Origin

`.planning/spikes/MANIFEST.md` (idea + scope + requirements), `.planning/spikes/PROPOSAL.md`
(consolidated roadmap, tiers 0–4, execution order), `sources/001-hotpath-perf-memory/` …
`sources/005-provenance-reproducibility/`. Also in user auto-memory as
`src-refactoring-proposal-spike.md`.
