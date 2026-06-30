---
phase: 03-training-data-pipeline
fixed_at: 2026-06-30T00:00:00Z
review_path: .planning/phases/03-training-data-pipeline/03-REVIEW.md
iteration: 1
findings_in_scope: 5
fixed: 5
skipped: 0
status: all_fixed
---

# Phase 03: Code Review Fix Report

**Fixed at:** 2026-06-30
**Source review:** .planning/phases/03-training-data-pipeline/03-REVIEW.md
**Iteration:** 1

**Summary:**
- Findings in scope: 5 (CR-01, WR-01, WR-02, WR-03, WR-04 — `critical_warning` scope; IN-01 deferred)
- Fixed: 5
- Skipped: 0

**Verification:** Full spike suite re-run in an isolated git worktree —
`julia --project=spike spike/test/runtests.jl` exits 0 (NeuralEstimators CPU
smoke 8/8, all simulator sets, Phase 3 Training-Data Pipeline **69/69** — one
more than the prior 68 because WR-04 adds a genuine different-seed assertion).
Decoupling gate re-checked: `git diff --quiet f581d95 -- src Project.toml
Manifest.toml` is clean; only files under `spike/` were modified.

## Fixed Issues

### CR-01: Bare `catch` in `_loads_ok` swallows `InterruptException`

**Files modified:** `spike/data/cache.jl`
**Commit:** c9f870b
**Applied fix:** Rebound the bare `catch` to `catch e` and added
`e isa InterruptException && rethrow()` before `return false`, so Ctrl-C (and
other interrupt-class exceptions) during the resume integrity check aborts the
`generate_cache` run instead of being silently consumed and treated as a
"not-done" shard. Torn-`.tmp` open/read errors still resolve to `false` as
intended.

### WR-01: `load_fold` returned `theta` as `Float64` (dtype inconsistency)

**Files modified:** `spike/data/loader.jl`
**Commit:** 78d764d
**Applied fix:** Wrapped both theta returns in `Float32.(...)` —
`θtr = Float32.(θ[:, train_idx])`, `θva = Float32.(θ[:, val_idx])` — so the
returned `(Ztr, θtr, Zva, θva)` tensors are uniformly Float32, matching the
docstring ("returns Float32 d×K tensors") and the Phase-4 NeuralEstimators
`train` input contract. SC-3/SC-4 still pass; SC-4's `θva` equality/inequality
assertions operate on the new Float32 matrices.

### WR-02: `seeding.jl` / `encode.jl` included unconditionally in `generate.jl`

**Files modified:** `spike/data/generate.jl`
**Commit:** 935800f
**Applied fix:** Added `isdefined` guards matching the neighbouring includes —
`isdefined(@__MODULE__, :HOLDOUT_SALT)  || include(".../seeding.jl")` and
`isdefined(@__MODULE__, :N_AUG_MOMENTS) || include(".../encode.jl")`. Guard
symbols verified present (`HOLDOUT_SALT` at seeding.jl:45, `N_AUG_MOMENTS` at
encode.jl:53), making repeated `include`/Revise re-includes idempotent and
const-redefinition-warning-free (and `--depwarn=error`-safe).

### WR-03: `write_meta` committed without a pre-commit integrity check

**Files modified:** `spike/data/cache.jl`
**Commit:** cb8431d
**Applied fix:** Reopened the `.tmp` and asserted the `hash` key is present
before the atomic `mv`, mirroring `write_shard`'s defensive
write-`.tmp` → reopen-verify → `mv` pattern. A partial/corrupt `meta.jld2` is
now caught before it can replace the committed manifest and brick
`open_or_invalidate` on the next run.

### WR-04: SC-4 `load_fold` reproducibility assertion was tautological

**Files modified:** `spike/test/test_data_pipeline.jl`
**Commit:** 67c21c0
**Applied fix:** Replaced the `f(x) == f(x)` self-comparison with two
independent `load_fold` calls saved to separate variables (`fa`, `fb`) plus a
third with a perturbed seed (`fc`, `master_seed + 1`), then asserted
`fa.θva == fb.θva` (same seed → identical D-09 fold membership) and
`fa.θva != fc.θva` (different seed → different split). This is now a genuine
reproducibility gate that would catch an off-by-one/stride bug in fold
membership rather than passing vacuously. Both assertions pass (Phase 3:
69/69).

---

_Fixed: 2026-06-30_
_Fixer: Claude (gsd-code-fixer)_
_Iteration: 1_
