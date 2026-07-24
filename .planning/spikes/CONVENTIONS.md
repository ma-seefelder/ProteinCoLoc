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

## Numerics / gate-diagnostic spikes (added in session 006-008)

- **Isolated figure environment.** `CairoMakie` is NOT a root dependency and the root
  co-resolution is fragile (the GLMakie/Makie 0.21-vs-0.24 conflict cost a Phase-7 plan). Render
  figures from `.planning/spikes/figenv` (`julia --project=.planning/spikes/figenv fig.jl`) and
  never add plotting deps to the root `Project.toml`/`Manifest.toml`.
- **Avoid new root deps in spikes generally** — a 3-line local `logsumexp` beats adding
  `LogExpFunctions` to the root manifest.
- **DEV seeds only, asserted disjoint.** Any spike touching the ship-gate machinery must run on a
  DEV seed and `@assert` it is disjoint from `PROD_SEED*`, `VAL_MASTER_SEED` (0x5BC0FFEE),
  `NPE_MASTER_SEED` (0xC0FFEE) and previously used dev seeds. Never consume a pre-registered seed
  in a spike, and never write a gate report from one.
- **Store the intermediate quantities, not just the verdict.** Spike 006 stored `p_post` per pair,
  which let Spike 007 answer a different question with *zero* additional inference cost. Cheap
  insurance against re-running expensive draws.
- **Validate a replacement estimator against the incumbent where both are well-defined** before
  interpreting where they differ. Spike 008's median agreement of 0.002 in the non-degenerate
  regime is what made its tail disagreement interpretable rather than suspect.
- **Report a refuted hypothesis as the finding.** Spike 008's proposed remedy failed; the failure
  identified the real problem (the baseline is unresolvable in the tail). Verdict PARTIAL, not a
  quiet rewrite of the question.

## Calibration / retrain-comparison spikes (added in session 009-013)

- **DEV-seed disjointness, made explicit and enumerated.** Every SBC/calibration spike ran on a DEV
  seed `@assert`ed disjoint from the *full enumerated* forbidden set — `PROD_SEED_V2`, every v1
  `PROD_SEED[G]`, `VAL_MASTER_SEED` (0x5BC0FFEE), `NPE_MASTER_SEED` (0xC0FFEE), and all prior dev
  seeds (012 asserted against 20 forbidden seeds by name). No production seed consumed, no gate
  report written, from any spike.
- **Separate artifact root + byte-identical baseline proof for any retrain.** A spike that retrains
  a net writes to its **own** artifact root (`artifacts/spike012_highcap/`, `artifacts/spike013_trunc/`)
  and captures `sha256_before.txt`, then re-verifies the frozen bundles are **byte-identical** after
  (amended_v2 `npe_8` sha256 `198bb078…` unchanged). This is how a retrain-comparison spike proves it
  did not perturb the pre-registered model.
- **Inject model changes through a seam, never by editing `src/`.** Spike 013 changed the prior via a
  custom `sample_prior_trunc`/`trunc_datagen` passed into `_train_grid_pipeline`'s `datagen` seam —
  `src/` stayed untouched while still training a genuinely different statistical model. Prefer an
  existing injection seam over any edit to frozen code.
- **Compare nets on the SAME DEV seed, and trust the band not the cell.** 012/013 evaluated both nets
  on one shared DEV seed (M=1000–2000, L=999) with atom-aware ρ_true handling. Per-parameter pass/fail
  carries ~±1.5 z-units of seed noise (Spike 010), so the reportable result is the **drift-magnitude
  band** (0.02–0.09 SD) and the **direction of change**, not an individual parameter's verdict on one
  seed. State the single-seed limit every time.
- **Randomized ranks are the standard fix for a mixed (discrete+continuous) prior.** When θ* can land
  on a prior atom (ρ_true = ghat(μ) clamped to ±0.99), plain SBC ranks are invalid. Replace atom-case
  ranks with fresh uniforms (or exclude atoms and report the excluded fraction) before the KS test —
  prefer this test-side correction over changing the model.
- **Read shrinkage as width, not location; state both readings of `vacuous`.** `vacuous = shrinkage
  ≥ 0.95` flags *uninformative about spread*, not *posterior = prior*. A vacuous parameter can still
  reject SBC on a location shift. Always report location (mean-u z) and shape (tail/middle) separately.
- **Justify a test-design change with an outcome-independent argument.** A practical-equivalence band
  for non-identified nuisances is defensible only via a power analysis that does **not** depend on the
  results (Spike 011) — the amendment §5 refused a results-dependent band to avoid a p-hacking
  appearance. Keep the justification outcome-independent, and route any change through §6.
