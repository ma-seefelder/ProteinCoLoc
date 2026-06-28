# Spike Conventions

Patterns established during the refactoring-analysis spike session. These are recommendations for
the eventual `src/` refactor, plus conventions for any future spikes in this repo.

## Stack (analysis spikes)

- **Micro-benchmarks:** standalone Julia scripts using only **stdlib** (`Statistics`, `Random`) +
  built-in `@allocated`/`@elapsed` (best-of-N with a warmup) — no `BenchmarkTools`, no heavy deps,
  no real images. Keeps benchmarks fast, reproducible, and decoupled from the `Turing`/`GLMakie`
  load. See `001-hotpath-perf-memory/bench_hotpath.jl`.
- **Synthetic data:** seed `Random.seed!`; simulate post-Otsu masking with a fraction of zeros.

## Structure

- One spike directory per analysis dimension: `NNN-<kebab-name>/` with `README.md` (frontmatter +
  findings table + per-finding detail) and any scripts/results.
- Findings carry stable IDs (`PERF-1`, `TYPE-2`, …) reused across the MANIFEST and PROPOSAL so they
  can be cross-referenced.

## Patterns (analysis method)

- **Find → adversarially verify.** Every proposed change is checked against the actual code by an
  independent skeptic that records a verdict (`confirmed` / `needs-nuance` / `overstated` /
  `incorrect`) + the evidence it found at the cited lines. Only verified, corrected claims go into
  the deliverables.
- **Cite `file.jl:line-range` for every claim**, and trust the verifier's file map over the brief —
  the analysis brief mislabeled `bayes.jl` code as `colocalization.jl`; the verifiers caught it.
- **Triangulation = signal.** The Δρ bug was found independently by 4 of 5 dimensions → highest
  confidence; treat cross-dimension agreement as a priority booster.
- **Separate correctness bugs from refactors** (the user asked for this): refactors are
  behavior-preserving; bugs are flagged with `is_correctness_bug` and listed separately.

## Recommended `src/` refactor conventions (for the build)

- **Type stability via function barriers:** resolve a `Symbol`/`Dict` selection to a *concrete*
  function once, then pass it into an inner `f(x,y, g::G) where G`. A bare `Val(runtime_symbol)` is
  **not** a fix (it is itself type-unstable).
- **Single source of truth for derived quantities:** e.g. one `delta_rho(::CoLocResult)` instead of
  7 inline copies — duplication is what let the Δρ bug drift in.
- **Real Julia interfaces:** abstract type declares erroring-fallback methods; concrete subtypes
  implement them; callers go through accessors, never `obj.field`. Preserve exception types
  (`ArgumentError`) when converting `hasfield`/`throw` checks.
- **Behavior-preserving first:** parametrize/annotate structs **without** renaming fields; stage
  renames behind property aliases. Update the `typeof(x)==CoLocResult` tests if you parametrize.
- **Provenance by construction:** thread an explicit `rng`/`seed`; embed run parameters in the
  result type; emit a structured `run_manifest.toml` (versions, seed, params, input hashes).
- **Headless plotting:** use `CairoMakie` for file output (PNG/SVG/PDF, no GL context); use scoped
  `with_theme(...) do … end`, never a global `set_theme!` that bleeds into the image plots.
- **Validate against the resolved stack:** the installed Turing (0.43–0.45) ≠ the legacy ADVI/`vi`
  API in `bayes.jl`; resolve `Manifest.toml` and port before relying on the posterior path.
