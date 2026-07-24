# 07-10 Summary — productionize the amortized-only public API (shipped family = {8})

**Plan:** 07-10 (wave 8) — public entry point + content-hashed lazy Artifacts registry + integration
test + docs. Adapted to Go/No-Go Option A (`07-GO-NO-GO-UPDATE.md`): `_SHIPPED_GRIDS = (8,)`.

## What shipped

### Task 1 — public entry point `colocalization_amortized` (`src/amortized/api.jl`, new)
`colocalization_amortized(img, control, channels; num_patches = 8, N = 2000, use_gpu = false)` →
`AmortizedColocResult`. It composes the already-promoted read surfaces (defines no new inference
math): `estimator_for(G)` → frozen summary (`_pair_mci` → `patch_summary` → `encode_d01` →
`standardize_summary`, mask rows bypassed) → `posterior_for` (un-standardized with `θzt`, Pitfall 5)
→ Δρ draw cloud (`rho_draws` difference, D-03) → `amortized_log_bf(pair_encode(Zs,Zc))` →
`ood_verdict`. `use_gpu = false` on every call (D-06). Input validation: exactly 2 channels,
in-range channels, finite standardized summaries, `N ≥ 1`. Exported from `src/ProteinCoLoc.jl`
(included after `local_map.jl`). The 4 accessors (`posterior_draws`/`delta_rho`/`bayes_factor`/
`is_ood`) dispatch on the result via `AbstractColocResult`.

### Task 2 — content-hashed Artifacts + registry population (`Artifacts.toml`, `src/registry.jl`)
- `_SHIPPED_GRIDS = (4,8,16,32)` → **`(8,)`**, with a comment citing the Go/No-Go decision and the
  per-grid exclusion reasons (4 never post-hoc re-analysed; 16 gate FAILED; 32 CAPPED; 64 dropped).
- `Artifacts.toml` written via the canonical `bind_artifact!` flow: a single `[grid_8]` entry,
  `lazy = true`, `git-tree-sha1 = 90e6b63a8a234d067b407fefd7914f2ae4845448`, with a `[[grid_8.download]]`
  stanza (GitHub-Release URL `…/releases/download/v2.0.0/grid_8.tar.gz` + real tarball
  `sha256 = 17162904…b6d5`). The tarball asset to upload is `artifacts/grid_8.tar.gz` (gitignored).
- `_lazy_load_from_artifact!(8)` wired: resolves the bundle from (1) the artifact store, (2) an
  in-repo dev-bundle fallback verified against the Artifacts.toml `git-tree-sha1`, or (3) a lazy
  Release download — **every path tree-sha1 verified before load (T-7-01)**. Then
  `_bundle_from_artifacts` + `_activate_recorded_ood_threshold` (wires the recorded density-channel
  ID operating point 179.14 from `gate_report_8.jld2`) + `register!`.
- **Shipping bundle = `artifacts/amended_v2/grid_8/`** — the retrained bounded-θ / realistic-imsize
  8×8 net the GO rests on (npe/ratio/ood_nulls + gate_report; `imsize_source = :generate_samples`).
  Verified: loads CPU-only, produces finite 7×N draws; `estimator_for(8)` returns it with `:thr`
  wired; full `colocalization_amortized` returns finite posterior/Δρ/logBF + a boolean OOD flag.
- Stdlibs `Artifacts` + `Pkg` added to `Project.toml [deps]` (`Pkg` moved out of `[extras]`); resolve
  clean, NeuralEstimators **0.2.1** / Flux **0.16.10** unchanged (Manifest delta is only the 2 stdlib
  entries — no downgrade).

### Task 3 — integration test + docs (`test/test_integration.jl`, `docs/amortized.md`)
- `test/test_integration.jl` (wired into `runtests.jl`): Artifacts.toml shape (git-tree-sha1 + lazy,
  no grid_4/16/32/64); exports (`colocalization_amortized` present; Turing path unexported); unshipped
  (64) + unregistered (7, mentions `train_and_register`) + channel-count error paths; and — forcing
  the lazy-load path by clearing `_REGISTRY[8]` — the full 8×8 public path with finite outputs and
  D-02 accessor dispatch.
- `docs/amortized.md`: the public API, the registry + `train_and_register`, the shipped-`{8}` family
  with exclusion reasons, the content-hashed lazy Artifacts load, the CPU-default note, the windowed
  sub-tile local map, and — verbatim from `07-GO-NO-GO-UPDATE.md §5` — the **NAMED LIMITS**
  (twice-amended gate / post-hoc verdict; ρ_true randomized-rank atom handling; simulation-based BF
  AUC 0.994, not a per-pair KDE reference; ~0.08-SD nuisance drift; ±1.5-z seed variance; F5
  image-size regime; density-only shipped OOD flag).

## Deviations from the plan (adapted to the GO)

- **`_SHIPPED_GRIDS = (8,)` not `(4,8,16,32)`** — the plan predates the calibration investigation and
  assumed grids "pass their gates". No grid passed literally; the GO rests on 8×8 (post-hoc randomized
  ranks + simulation-based BF). Only 8×8 ships.
- **32 absent from the registry** — 07-09 CAP decision (no `grid_32` entry).
- The plan's "only gate-PASSED grids populate `_SHIPPED_GRIDS`" is honored in spirit: **no grid passed
  its literal gate**, so the shipped grid is the one whose corrected post-hoc re-analysis carries the
  GO — 8×8 — and its named limits are documented rather than hidden behind a "pass" claim.

## Honesty notes

- The shipped OOD flag is the **density Mahalanobis channel** at the recorded ID operating point; the
  fused noise/PP channels that reached gate AUC 1.0 are gate-time constructs not carried in the frozen
  bundle. Documented in `docs/amortized.md` (limit 7) and the `api.jl` docstring.
- The `AmortizedColocResult.calibration` is the `_bundle_from_artifacts` placeholder
  (`gated = false`); the gate provenance and named limits live in `docs/amortized.md` and the memo.

## Constraints honored

`spike/` byte-untouched. No `test/gate/gate_consts_*.jl` touched. No retraining, no gate run, no
`PROD_SEED_V2` consumed. Existing frozen artifacts (`amended_v2`, `grid_4/8/16`, `spike01x`)
byte-identical (`create_artifact` copies; it never mutates the source). Full `Pkg.test()` run at the
end.
