---
phase: 12-spatial-colocalization-map
plan: 03
status: complete
subsystem: lattice prior arithmetic
tags: [car, gp, matern32, dct-ii, rescale, r1-reparametrization, column-major, zero-deps]
requires: ["12-01"]
provides:
  - "spike/simulator/p12_lattice.jl — CAR and GP covariance arms, sqrt(diag) rescaled field sampler, induced-lag-1 reparametrization of both arms, exact orthonormal DCT-II basis with a frozen smoothness ordering, radial design matrix for the S-4 guard"
  - "spike/test/test_p12_lattice.jl — nine testsets, 137 assertions, replacing the 12-01 pending scaffold"
affects: []
tech-stack:
  added: []
  patterns:
    - "stdlib-only numeric kernel (LinearAlgebra + SparseArrays via the default @stdlib LOAD_PATH), zero Project.toml change"
    - "falsifier-paired assertions: every property that could pass vacuously carries a companion assertion that fails when the property is removed"
    - "the marginal rescale lives in the SAMPLER, not in Sigma, so only a DRAW can verify it"
key-files:
  created:
    - spike/simulator/p12_lattice.jl
  modified:
    - spike/test/test_p12_lattice.jl
decisions:
  - "The unit-marginal claim is asserted on 20 000 draws, never on `Sigma ./ (sd*sd')` — that expression has a unit diagonal for ANY SPD Sigma and is green whether or not the sampler ever divides. Confirmed load-bearing by mutation: max |sd − 1| = 0.0259 with the rescale, 0.3635 without."
  - "The header comment does NOT name the forbidden third-party GP/random-field packages, so their absence stays greppable. The first draft declared their absence by naming them, which would have tripped the plan's own acceptance grep."
  - "`gp_sigma` attempts the Cholesky itself and throws with the failing (ell, kernel, jitter) rather than returning a matrix that fails later somewhere unrelated — the T-12-16 mitigation is only real if the throw happens where the jitter is chosen."
metrics:
  duration: ~55 min
  completed: 2026-07-29
  tasks: 2
  commits: 3
  files: 2
---

# Phase 12 Plan 03: Lattice Kernels, Rescale, r₁ Reparametrization and the DCT-II Basis — Summary

Both spatial prior arms now exist as 64×64 dense covariances addressed by ONE physically
meaningful knob (induced lag-1 correlation, solved to 1e-6 by bisection over the pre-registered
brackets), draws from either are marginally standard **per cell** rather than only on average,
and the deviation basis is an exact orthonormal bijection with a frozen smoothness ordering —
all on stdlib `LinearAlgebra`/`SparseArrays` with `spike/Project.toml` and `spike/Manifest.toml`
byte-unchanged.

## What Was Built

**Task 1 — `spike/simulator/p12_lattice.jl`** (commit `4c6e13a`, 441 lines, new file)

Flat top-level functions, no `module`, `#= =#` AGPL header, guarded include of the frozen
Tier-1 pre-registration:

| Function | What it delivers |
|---|---|
| `p12_idx(i, j; G)` | the single statement of the column-major contract, `(j - 1) * G + i` |
| `lattice_adjacency(G)` | 4-neighbour sparse adjacency in `p12_idx` order |
| `car_sigma(G, α)` | `Σ = (D − αW)⁻¹`, `Symmetric`, guarded to `0 ≤ α < 1` |
| `gp_sigma(G, ℓ; jitter, kernel)` | Matérn-3/2 **default**, `:se` reachable; pre-registered jitter; throws on a failed Cholesky with the failing parameters |
| `field_sampler(Σ)` | the closure `rng -> (L*randn(rng,n)) ./ sqrt.(diag(Σ))` |
| `induced_lag1(Σ; G)` | mean Pearson correlation over 4-neighbour pairs of `Σ ./ (sd*sd')` |
| `car_alpha_for_r1`, `gp_ell_for_r1` | bisection to `P12_R1_SOLVE_TOL`, monotonicity asserted on the bracket endpoints, out-of-range requests refused |
| `lattice_sigma(arm, r1; G, kernel)` | `:car` / `:gp` / `:none`; `:none` (or `r1 == P12_ABLATION_R1`) is `Symmetric(Matrix(1.0I, 64, 64))` |
| `p12_dct_matrix`, `p12_dct`, `p12_idct`, `p12_dct_vec`, `p12_idct_vec`, `p12_dct_order` | the exact orthonormal basis and its frozen smoothness ordering |
| `radial_basis(G)` | `[1 radius]` design matrix for the S-4 radial-energy guard |

Nothing hard-codes 8 in a body; every function takes `G` as an argument or a keyword defaulting
to `P12_G`.

**Why the GP default is Matérn-3/2 and not SE, measured rather than asserted.** At the top of
the useful range the SE kernel is badly conditioned and Matérn-3/2 is not — measured this
session at ℓ = 8: `cond(K) = 5.49e9` for `:se` against `9.89e4` for `:matern32`, a factor of
5.6e4. `:se` stays reachable by keyword so the mini-spike can report both.

**Task 2 — `spike/test/test_p12_lattice.jl`** (commit `02ddaf3`, 198 lines, scaffold replaced)

Nine testsets, 137 assertions, all randomness on `p12_fix_rng(P12_FIXTURE_COUNTER)` — the
fixture seed on the fixture counter. `p12_rng` (the reported stream) is not touched anywhere in
the file. The `P12_PENDING_SCAFFOLD` marker is gone (`grep -c` returns 0).

## Verification Results — what was actually run and observed

| Check | Command | Observed |
|---|---|---|
| Task 1 verify (verbatim from the plan) | the plan's `julia -e` one-liner | **PASS** — printed `0.5 lattice-ok` |
| induced r₁ solved at the printed rung | same run | `0.499999140` — \|err\| = 8.6e-7, inside `P12_R1_SOLVE_TOL` = 1e-6 |
| DCT orthonormality residual | `maximum(abs, C*C' - I)` | **6.66e-16**, below the 1e-12 acceptance bar |
| per-cell marginal sd, 4000 draws, **with** rescale | `maximum(abs, std(X;dims=2) .- 1)` | **0.0259** against the 0.06 bar |
| per-cell marginal sd, 4000 draws, **rescale deleted** (scratch copy) | same assertion | **0.3635 — verify FAILS.** The assertion is load-bearing |
| Task 2 verify | `julia --project=spike -e 'using Test; include("spike/test/test_p12_lattice.jl")'` | **PASS** — 137 passed, 0 failed, 0 errored, exit 0 |
| testset 2 falsifier bites | real testsets run against the rescale-deleted scratch copy | **PASS** — 2 passed / **1 failed** at `test_p12_lattice.jl:95`, and the suite threw |
| **nine testsets green under `runtests.jl`** | `julia --project=spike spike/test/runtests.jl` | **PASS** — 29 / 3 / 4 / 13 / 4 / 10 / 2 / 71 / 1 = **137, zero failures, zero errors**, all after the `P12-SUITE-RAN: test_p12_lattice.jl` marker |
| `spike/Project.toml` / `Manifest.toml` unchanged | `git diff --quiet HEAD --` | **PASS** |
| `src/` and `corpus/` unchanged | `git diff --quiet HEAD -- src corpus` | **PASS** |
| `.planning/STATE.md`, `.planning/ROADMAP.md`, `p12_consts.jl` unchanged | `git diff --quiet HEAD --` | **PASS** |
| `spike/data/cache/p11` intact | `du -sh` | **PASS** — 54 MB |
| no `Pkg.add`, no `push!(LOAD_PATH`, no third-party GP/random-field package name | `grep -E` over the file | **PASS** — zero matches (see Deviation 1) |
| `(j - 1) * G + i` present; `(i - 1) * G + j` absent | `grep` | **PASS** — 1 and 0 |
| `lattice_sigma(:none, ·)` is exactly `Symmetric(Matrix(1.0I,64,64))` | `==` | **PASS** |
| `P12_PENDING_SCAFFOLD` gone | `grep -c` | **PASS** — 0 |
| r₁ ladder, both arms | all five `P12_R1_LADDER` rungs | **PASS** — max \|err\| 1.3e-6, bar 1e-5 |

**Suite exit code is still non-zero and it is not ours.** `runtests.jl` aborts at
`test_npe.jl:230` on the pre-existing Phase-4 `SPEEDUP_GATE` (`median_speedup` measured
**88.60** this run against the pre-registered bar 100.0 — a third independent measurement, after
84.44 and 89.77 during 12-01). That gate was not touched: re-deriving the bar is a
pre-registration decision the user owns (standing project ruling). Every Phase-12 testset runs
and reports **before** that abort, exactly as 12-01's hardened wiring was built to guarantee.

**Not run:** `julia --project=. -e 'using Pkg; Pkg.test()'`. No pre-plan baseline for it exists
(12-01 recorded the same gap), so "unchanged" is not something this run can assert. `src/` is
byte-unchanged by `git diff --quiet` and this plan adds nothing to the main package and loads
nothing from it, so the main-package suite cannot have been affected. Recorded rather than
implied.

## Deviations from Plan

### 1. [Rule 2 — Correctness] The header could not declare the forbidden package names by naming them

- **Found during:** Task 1, while checking the acceptance criteria before committing.
- **Issue:** The first draft's dependency-discipline header read "there is deliberately no
  AbstractGPs / KernelFunctions / GaussianRandomFields anywhere". That sentence makes the
  acceptance criterion *"the file contains … no reference to `AbstractGPs`, `KernelFunctions`
  or `GaussianRandomFields`"* **false as written**, and would trip any grep-based decoupling
  assertion 12-05 later adds — a comment claiming a property while breaking the check for it.
- **Fix:** The paragraph now states the same prohibition without the literal names ("nothing
  below … reaches for a third-party Gaussian-process / random-field library … every such name is
  deliberately absent from this file so its absence stays greppable"), and the same for
  `Pkg.add` / `push!(LOAD_PATH`. Verified: `grep -E "Pkg\.add|push!\(LOAD_PATH|AbstractGPs|KernelFunctions|GaussianRandomFields"` returns zero matches.
- **Commit:** `4c6e13a`.

### 2. Assertions beyond the plan's enumerated list, in four testsets

- **Found during:** Task 2.
- **What:** testset 3 adds `0 < r < 1` on both ladders; testset 4 adds
  `@test_throws ArgumentError car_alpha_for_r1(0.99999999)` and
  `lattice_sigma(:nonsense, 0.5)`; testset 5 adds `lattice_sigma(:car, P12_ABLATION_R1)` and the
  ablation's own marginal band; testset 6 adds the `p12_dct_vec`/`p12_idct_vec` flat round-trip
  and `sort(p12_dct_order(8)) == 1:64`; testset 7 adds `p12_dct(ones)[1,1] ≈ 8`; testset 8 adds
  the adjacency and radial-basis indexing checks and `p12_idx(3,6) != p12_idx(6,3)`.
- **Why:** each covers a surface the plan's item list creates but does not assert. The flat
  `_vec` forms are the shape θ rows are actually stored in, so a matrix-only round-trip would
  leave the *used* path untested; `p12_idx(3,6) != p12_idx(6,3)` is what gives the column-major
  assertions content (on a symmetric fixture they would pass under a transpose). No testset was
  added or removed — the plan's nine are exactly what shipped.

### 3. `M::Matrix{Union{Float64,Missing}}` is compared with `isequal`, not `==`

- **Found during:** Task 2, testset 8.
- **Issue:** the plan asks for `vec(M)[p12_idx(i,j)] == M[i,j]` and `reshape(vec(M),8,8) == M` on
  a `Union{Float64,Missing}` matrix. With a real `missing` present both expressions evaluate to
  `missing`, and `@test missing` **errors** rather than passing or failing.
- **Fix:** `isequal` throughout, and the fixture deliberately carries two `missing` entries so
  the Union eltype is exercised rather than merely declared. Writing `M` with no missings would
  have satisfied `==` while testing nothing about the eltype `encode_d01` actually receives.
- **Commit:** `02ddaf3`.

## Assumption Drift (advisory)

**The research's CAR lag-1 table is not reproduced by the estimator the plan defines, and the GP
one is — exactly.**

- **Planned:** `12-RESEARCH.md` Pattern 1's measured table (CAR α = 0.5/0.9/0.95/0.99/0.999 →
  r₁ = 0.136/0.351/0.438/0.692/0.947) describes the same quantity `induced_lag1` computes, and
  `p12_consts.jl:301-304` quotes its endpoints (0.136 and 0.947) as the R-8 motivation.
- **Actual, measured this session** with the plan's definition ("mean Pearson correlation over
  all 4-neighbour pairs of `Σ ./ (sd*sd')`"):

  | α | plan's estimator | `12-RESEARCH.md` |
  |---|---|---|
  | 0.5 | **0.1561** | 0.136 |
  | 0.9 | **0.4117** | 0.351 |
  | 0.95 | **0.5064** | 0.438 |
  | 0.99 | **0.7222** | 0.692 |
  | 0.999 | **0.9447** | 0.947 |

  The GP row, by contrast, reproduces to four decimals on the same code path
  (`:se` at ℓ = 0.5/1/2/4/8 → 0.1353 / 0.6065 / 0.8825 / 0.9692 / 0.9922 against the document's
  0.135 / 0.607 / 0.882 / 0.969 / 0.992).
- **Why:** a stationary GP kernel gives every 4-neighbour pair the *same* correlation, so every
  reasonable pooling of "lag-1" coincides — which is why that row matches exactly. The CAR arm
  is **not** stationary (boundary cells have degree 2 or 3), so the answer depends on how pairs
  are pooled, and `12-RESEARCH.md` records the numbers without recording the pooling. Two
  alternative poolings were tried and neither reproduces the table either (variance-weighted
  pooled Pearson: 0.1596/0.419/0.5132/0.723/0.9445; interior-pairs-only:
  0.1376/0.3743/0.4716/0.7119/0.9475).
- **Why this is advisory and not a blocker.** Nothing gated depends on the α↔r₁ map: the phase
  parametrizes by r₁, so α is always *derived* by bisection through whichever `induced_lag1` the
  code defines, and no pre-registered constant asserts an α value. The R-8 conclusion is
  unchanged and if anything is restated more conservatively by the new numbers (α = 0.5 still
  buys r₁ = 0.156, and the plan's own falsifier bar `< 0.2` holds). The `p12_consts.jl:301-304`
  endpoints are **prose motivation in a comment, not an assertion**, so Tier 1 is neither
  violated nor in need of an append. Recorded here so that no later document quotes the
  research table and the executed code as if they were the same measurement.

## Known Stubs

None. Both files are complete implementations; nothing in this plan is deferred to a later one.

## Threat Flags

None. No new network, auth, file-access or schema surface: the plan adds pure in-memory linear
algebra and a test file, writes no artifact, and reads only the frozen Tier-1 constants.

## For the Next Plan

- **12-04** (`p12_architecture.jl`) inherits the column-major contract through `p12_idx`; the
  `reshape(rows[1:64], 8, 8)` it needs is asserted transpose-free in testset 8 here.
- **12-07** (`p12_prior.jl`) is the file that owns the copula. `field_sampler` already delivers
  the marginally-standard `z` the copula assumes, so 12-07 should NOT re-standardize; SIM-02
  per region is testable directly against `MU_PRIOR ∘ ĝ`.
- **`P12_K_DEV = 63` deviation rows are `p12_dct_order(8)[2:end]`**, and the global term is
  coefficient index 1. That ordering is frozen by `p12_dct_order`, so a later truncation study
  (12-15) can take "the first K modes" literally.
- **`radial_basis(8)` is the S-4 guard's geometry.** Do not re-derive it in the guard runner — a
  second derivation is a second chance to transpose it.
- The Phase-4 `SPEEDUP_GATE` still aborts `runtests.jl` (88.60 vs 100.0 this run). Read the
  `P12-SUITE-RAN` markers and the Phase-12 testset summaries; do not read a red suite as a
  Phase-12 failure.

## Self-Check: PASSED

Both files verified present on disk (`spike/simulator/p12_lattice.jl`,
`spike/test/test_p12_lattice.jl`); both task commits verified present in `git log`: `4c6e13a`
(lattice arithmetic), `02ddaf3` (test file). `spike/Project.toml`, `spike/Manifest.toml`, `src/`
and `corpus/` byte-unchanged; `spike/validation/p12_consts.jl` byte-unchanged (Tier 1 not
appended to, as this plan opens no sentinel); `.planning/STATE.md` and `.planning/ROADMAP.md`
untouched; `spike/data/cache/p11` intact at 54 MB; `test_npe.jl` and `SPEEDUP_GATE` untouched.
All commits used explicit pathspecs; no `git stash`, `rebase`, `amend`, `reset` or `clean` was
run, and no index.lock race occurred despite the concurrent Phase-13 executor on the same branch.
