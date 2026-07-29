---
phase: 12-spatial-colocalization-map
plan: 04
status: complete
subsystem: spatial input contract + the pinned convolutional head
tags: [cnn, reshape, column-major, mcar-masking, normalising-flow, theta-layout, r-7, r-3, cpu-only]
requires: ["12-01"]
provides:
  - "spike/npe/p12_architecture.jl — reshape_summary (transpose-free (G,G,2,K) view of the UNCHANGED 128-row summary), mask_regions / sample_mask_k / augment_mask! (the Pitfall-4 MCAR augmentation), the LEAN pre-registered CNN summary net (70 872 params), build_p12_estimator (CNN + a 72-marginal NormalisingFlow), P12_THETA_ROWS + p12_theta_index + the four named row positions, and P12_ZSCORE_ARM = :per_row"
  - "spike/test/test_p12_architecture.jl — nine testsets, 1141 assertions, replacing the 12-01 pending scaffold"
affects: []
tech-stack:
  added: []
  patterns:
    - "the summary STATISTIC is unchanged; only its VIEW changes — reshape_summary is a pure re-view of the same numbers"
    - "region-unobserved is expressed in encode_d01's EXISTING value-0/mask-0 vocabulary, never in a new sentinel"
    - "anti-pattern assertions run on the COMMENT-STRIPPED source, so the prohibition can be documented in the file it constrains"
    - "use_gpu = false is passed EXPLICITLY on every NeuralEstimators call because the v0.2.1 default is true"
key-files:
  created:
    - spike/npe/p12_architecture.jl
    - .planning/phases/12-spatial-colocalization-map/deferred-items.md
  modified:
    - spike/test/test_p12_architecture.jl
decisions:
  - "The rejected fat topology is described in PROSE (a 2 -> 32 -> 64 -> 64 stack flattening to 4096 into a 4096 -> 256 Dense) rather than transcribed as Julia, because the plan's own acceptance criterion forbids the literals `Dense(4096` and `=> 64,` ANYWHERE in the file. Writing the comment literally would have made the file fail the criterion the comment exists to justify — the same trap 12-03 hit with the forbidden GP package names."
  - "`sampleposterior` is called with an EXPLICIT `use_gpu = false`, which the plan's task body omits. The pinned v0.2.1 signature (`PosteriorEstimator.jl:130`) defaults `use_gpu` to TRUE, so the omission would opt the call into the GPU path on any machine that has one — and the plan's own acceptance criterion demands `use_gpu = false` on every NeuralEstimators call."
  - "The T-12-18 pooling-absence assertion strips `#=` blocks, `\"\"\"` docstrings AND `#` lines before searching. The plan requires the anti-pattern to be NAMED in the file's prose and the docstring to record the `Unet` fact; an assertion on the raw text would therefore forbid documenting the very rule it enforces. Stripping first keeps the check pointed at the executable topology, which is the only place pooling could do harm."
metrics:
  duration: ~70 min
  completed: 2026-07-29
  tasks: 2
  commits: 2
  files: 3
---

# Phase 12 Plan 04: Spatial View, the Lean CNN and the 72-Row θ Layout — Summary

The unchanged 128-row summary is now readable as an 8×8×2 lattice with the column-major
contract asserted rather than assumed, the pre-registered LEAN convolutional head is pinned at
70 872 parameters (the rejected alternative is 1 137 760), a 72-marginal `NormalisingFlow`
trains one epoch and samples 72×32 on CPU in ~44 s with no custom `ApproximateDistribution`, and
the θ row order is fixed, named and asserted at load — `spike/Project.toml` and
`spike/Manifest.toml` byte-unchanged.

## What Was Built

**Task 1 — `spike/npe/p12_architecture.jl`** (commit `0b36c16`, 399 lines, new file)

Flat top-level functions, no `module`, `#= =#` AGPL header, guarded includes of the frozen
Tier-1 pre-registration and then `p12_lattice.jl` (for the shared `p12_idx` contract).

| Function / constant | What it delivers |
|---|---|
| `reshape_summary(Zraw, G = P12_G)` | `(2G²)×K → (G,G,2,K)`; channel 1 correlations, channel 2 present-mask; `DimensionMismatch` on a wrong row count; a vector method for the single-data-set case |
| `mask_regions(Zraw, idxs; G)` | a COPY with the continuous AND mask row of each named region zeroed — the SAME encoding `encode_d01` already emits below the ≥15-pixel floor |
| `sample_mask_k(rng; k_set)` | one draw from the Tier-1 `P12_MASK_K_SET = 0:8` |
| `augment_mask!(Zraw, rng; G, k_set)` | in-place MCAR augmentation, an independent `k` per column |
| `build_p12_summary_net(; G, dstar)` | the LEAN topology, layer stack transcribed literally, head width `Dense(8 * G^2, dstar)` DERIVED from `G` |
| `build_p12_estimator(; G, D, dstar, num_coupling_layers, flow_depth, flow_width)` | `PosteriorEstimator(network, q)` with `q` a `NormalisingFlow` INSTANCE passed POSITIONALLY; `ArgumentError` guards on `D ≥ 2`, `dstar ≥ 1`, `G ≥ 2`, `num_coupling_layers ≥ 1` |
| `P12_THETA_ROWS(K)`, `p12_theta_index(sym; K)` | the row order and its name→position map, `@assert`ed at load |
| `p12_row_c0`, `p12_rows_dev`, `p12_rows_nuisance`, `p12_row_r1` | `1`, `2:64`, `65:71`, `72` |
| `P12_ZSCORE_ARM` | `:per_row` — the recorded R-3 arm (T-12-19) |

Nothing hard-codes 8 or 512 in a body: `build_p12_summary_net(; G = 6, dstar = 32)` builds and
runs, asserted in testset 3.

**The parameter count is not a claim, it is arithmetic that reproduces the research table
exactly.** `304 + 4640 + 264 + 65 664 = 70 872`, which is the `12-RESEARCH.md:1176` figure to
the digit. This is the first Phase-12 measurement that has reproduced a research-table row on
the nose (12-03's CAR lag-1 row did not).

**Task 2 — `spike/test/test_p12_architecture.jl`** (commit `17e798e`, 334 lines, scaffold
replaced)

Nine testsets, 1141 assertions. All DATA randomness rides `p12_fix_rng(P12_FIXTURE_COUNTER)`;
`p12_rng` is not touched anywhere in the file. Flux weight init rides the task-local RNG, seeded
from a locally constructed `Xoshiro` — a different RNG family entirely, so it cannot collide
with the Philox key space. `P12_PENDING_SCAFFOLD` is gone (`grep -c` returns 0). The outer
testset carries `verbose = true` so all nine names print with their own pass counts even when
green: a collapsed summary cannot distinguish "nine testsets passed" from "one passed and eight
were never written".

## Verification Results — what was actually run and observed

| Check | Command | Observed |
|---|---|---|
| Task 1 verify (verbatim from the plan) | the plan's `julia -e` one-liner | **PASS** — printed `(128, 4) 70872 arch-ok`, exit 0 |
| parameter count in the lean band | same run | **70 872**, inside the plan's 60 000 – 90 000 acceptance band |
| literal layer specs present | `grep -c` for `2  => 16` / `16 => 32` / `1,1` / `32 =>  8` | **2 / 2 / 2 / 2** (code + the docstring's rendering) |
| `Dense(8 * G^2` present, `Dense(4096` and `=> 64,` absent | `grep -c` | **1 / 0 / 0** |
| pooling absent from the comment-stripped source | `#=`, `"""` and `#` stripped, then searched | `GlobalMeanPool` 0, `DeepSet` 0, `Unet` 0, `mean(` 0 |
| θ layout | `P12_THETA_ROWS()` | length **72**, first `:c0`, last `:r1`, `[65:71]` the seven nuisances in the stated order; `p12_theta_index(:r1) = 72`, `(:c0) = 1`, `length(p12_rows_dev) = 63` |
| Task 2 verify | `julia --project=spike -e 'using Test; include("spike/test/test_p12_architecture.jl")'` | **PASS — 1141 passed, 0 failed, 0 errored**, exit 0 |
| all nine testsets visible | same run | 32 / 1042 / 21 / 7 / 3 / 23 / 10 / 2 / 1 = **1141** |
| testset 5 (train + sample) wall clock | same run | **44.2 s**, under the 60 s bar; `train` one epoch then `sampleposterior` → `(72, 32)`, all finite |
| **testset 1 falsifier bites** | `reshape_summary` mutated to a transposed variant in place, real testsets re-run | **PASS — 16 passed / 16 FAILED** in testset 1 alone (lines 101,102,105,106,110-113,122,123 and the per-index loop), plus 2 knock-on failures in testset 2. The round-trip assertion is load-bearing |
| falsifier reverted | restored from a scratchpad backup, then `git diff --quiet HEAD --` | **PASS — byte-identical to `0b36c16`** |
| nine testsets green **under `runtests.jl`** | `julia --project=spike spike/test/runtests.jl` | **PASS** — 1141 / 0 / 0, printed after the `P12-SUITE-RAN: test_p12_architecture.jl` marker |
| `spike/Project.toml` / `Manifest.toml` unchanged | `git diff --quiet HEAD --` | **PASS** |
| `src/` unchanged | `git diff --quiet HEAD -- src` | **PASS** |
| `spike/validation/p12_consts.jl` unchanged | `git diff --quiet HEAD --` | **PASS** (this plan opens no Tier-2 sentinel) |
| `.planning/STATE.md`, `.planning/ROADMAP.md` unchanged | `git diff --quiet HEAD --` | **PASS** (after `git update-index --refresh`; the first call returned 1 on stale index stat data, not on content) |
| `spike/data/cache/p11` intact | `du -sh` | **PASS** — 54 MB |
| no stray files written by `train` | `git status --short` | **PASS** — no new untracked files; NeuralEstimators' default `savepath` writes nothing |

### The full-suite exit code is non-zero, and this time it is a DIFFERENT reason than the briefing predicted

The carry-forward I was briefed on — `test_npe.jl` aborting on the Phase-4 `SPEEDUP_GATE` —
**did not happen this run. Phase 4 SC3 (NPE-03) passed 9/9.** Execution therefore reached
`test_sbc.jl` (five includes later) for the first time in this phase, and hit a **latent
sentinel collision that has been present since 12-01** and was being masked by that very abort:

`spike/validation/p12_consts.jl:103` re-declares `const VAL_MASTER_SEED` (same VALUE as
`consts.jl:76` — `0x5BC0FFEE` — so no seed changed and no pre-registration was breached), and
`spike/validation/harness.jl:51` guards its include of the Phase-5 pre-registration on exactly
that NAME. With the Phase-12 aggregator wired FIRST at `runtests.jl:174`, `consts.jl` is
skipped, `SBC_FIX_M` never exists, and `test_sbc.jl:41` throws.

Reproduced in one process:

```
before:            VAL_MASTER_SEED = false
after p12_consts:  VAL_MASTER_SEED = true    SBC_M = false
after harness:     SBC_M = false   SBC_FIX_M = false   load_frozen_model = true
```

`load_frozen_model` is defined while `SBC_FIX_M` is not — i.e. `harness.jl` loaded and the
Phase-5 constants file specifically did not.

**This is not 12-04's.** The two files this plan touches are not in that chain: the Phase-12
aggregator provably leaves `SBC_M`, `SBC_FIX_M` and `sbc_ranks` undefined (checked directly),
and both offending lines predate this plan. Fixing it would require editing FROZEN Tier-1
`p12_consts.jl` (append-only; this plan opens no sentinel) or `harness.jl` / `runtests.jl`,
neither of which is in this plan's `files_modified` — and the FIRST-position ordering is 12-01's
asserted structural decision. Logged in full, with the one-line preferred fix, as **DEF-12-01**
in `.planning/phases/12-spatial-colocalization-map/deferred-items.md`. It should be resolved
before any plan depends on a green full-suite run.

**Not run:** `julia --project=. -e 'using Pkg; Pkg.test()'`. No pre-plan baseline for it exists
(12-01 and 12-03 recorded the same gap), so "unchanged" is not something this run can assert.
`src/` is byte-unchanged by `git diff --quiet` and this plan adds nothing to the main package and
loads nothing from it, so the main-package suite cannot have been affected. Recorded rather than
implied.

## Deviations from Plan

### 1. [Rule 2 — Correctness] `sampleposterior` needed an explicit `use_gpu = false`

- **Found during:** Task 2, checking the acceptance criteria before committing.
- **Issue:** The plan's Task-2 body writes `sampleposterior(est, Z_train[:, :, :, 1:1]; N = 32)`
  with no `use_gpu`, while its own acceptance criteria say *"Every NeuralEstimators call in the
  file passes `use_gpu = false`"* and the action text says *"`use_gpu = false` on EVERY
  NeuralEstimators call, no exceptions."* Reading the pinned source settles which is right:
  `PosteriorEstimator.jl:130` is
  `sampleposterior(estimator::PosteriorEstimator, Z; N = 1000, device = nothing, use_gpu::Bool = true, kwargs...)`
  — **the library default is TRUE.** Omitting it opts the call into the GPU path on any machine
  that has one, which is the D-10 / Pitfall-1 trap, and would make testset 9's CPU-only claim
  true only by accident of hardware.
- **Fix:** `sampleposterior(...; N = 32, use_gpu = false)`, with the signature quoted in a
  comment so the next reader does not "simplify" it away.
- **Commit:** `17e798e`.

### 2. [Rule 3 — Blocking] The rejected fat topology is described in prose, not transcribed

- **Found during:** Task 1, checking the acceptance criteria before committing.
- **Issue:** The plan asks the comment to record *"the fat alternative (`2→32→64→64`, flatten
  4096, `Dense(4096,256)`)"*, and its own acceptance criterion says the file *"does NOT contain
  `Dense(4096` or `=> 64,`"*. Transcribing the comment literally makes the file fail the
  criterion the comment exists to justify. This is the same trap 12-03 hit with the forbidden
  GP package names (12-03 Deviation 1).
- **Fix:** the cost comment names the same architecture in prose — "a 2 -> 32 -> 64 -> 64
  convolution stack with no channel bottleneck, flattening 8*8*64 = 4096 into a 4096 -> 256
  Dense" — carrying every number (1 137 760 params, 12.5 h vs 1.8 h, the 92 % figure) without
  the two forbidden literals. Verified: `grep -c` returns 0 for both over the whole raw file.
- **Commit:** `0b36c16`.

### 3. A vector method on `reshape_summary`, required by the plan's own test text

- **Found during:** Task 2.
- **Issue:** the plan's signature is `reshape_summary(Zraw::AbstractMatrix, G)`, but its Task-2
  text calls `reshape_summary(Z, 8)` where `Z = vcat(vec(...), vec(...))` is a 128-element
  VECTOR. As written the two do not compose.
- **Fix:** `reshape_summary(v::AbstractVector, G) = reshape_summary(reshape(v, :, 1), G)`, which
  is exactly the precedent `p11_architecture.jl:335` sets for `augment_input`. `mask_regions`
  gets the matching vector method for the same reason.
- **Commit:** `0b36c16`.

### 4. The T-12-18 assertion strips docstrings as well as `#` lines

- **Found during:** Task 2.
- **Issue:** the plan requires the anti-pattern to be NAMED in the file (*"record the
  ANTI-PATTERN explicitly: no `GlobalMeanPool`, no `DeepSet` …"*) and requires the DOCSTRING to
  record that `Unet` is not exported by v0.2.1 — while T-12-18's mitigation is a source
  assertion that those names appear nowhere. A raw-text assertion would forbid documenting the
  rule it enforces.
- **Fix:** `_p12_strip_noncode` removes `#= =#` blocks, `"""…"""` docstrings and `#`-prefixed
  lines before searching, so the check is pointed at the executable topology — the only place
  the pooling could do harm. The reason is written into the helper.
- **Commit:** `17e798e`.

### 5. Assertions beyond the plan's enumerated list, in four testsets

- **Found during:** Task 2.
- **What:** testset 1 adds the per-index `p12_idx` cross-check, the explicit transposed-variant
  inequality, and a multi-column no-bleed check; testset 2 adds the not-mutated /
  bit-identical-elsewhere pair and the out-of-range `ArgumentError`s; testset 3 adds a
  `G = 6, dstar = 32` build proving the head width follows `G`, and the raw-file `Dense(4096` /
  `=> 64,` absence; testset 6 adds `!(:ρ_true in rows)` and a `K = 3` genericity check.
- **Why:** each covers a surface the plan's item list creates but does not assert. The
  transposed-variant inequality in particular is what stops the two equality assertions being
  jointly satisfiable on an accidentally symmetric fixture — the fixture's missingness pattern
  is deliberately asymmetric for the same reason. **No testset was added or removed — the plan's
  nine are exactly what shipped.**

## Assumption Drift (advisory)

**The full-suite abort point moved, and the reason the briefing gave for it is no longer true.**

- **Planned / briefed:** *"`runtests.jl` exits NON-ZERO for a pre-existing Phase-4 reason Phase
  12 does not own (`test_npe.jl` SPEEDUP_GATE, measured 84.44 / 88.60 / 89.77 against a bar of
  100)."* Three prior sessions measured the same shortfall, so the working assumption was that
  the abort point is stable and everything after `test_npe.jl` is simply unreached.
- **Actual, this session:** Phase 4 SC3 (NPE-03) **passed 9/9** — the speedup cleared 100 — and
  the suite instead aborted five includes later at `test_sbc.jl:41`.
- **Why it matters:** the two facts are causally linked. The `SPEEDUP_GATE` abort was
  *masking* DEF-12-01, so "the suite is red for a known Phase-4 reason" was quietly doing double
  duty as "nothing after `test_npe.jl` is broken". It is not: `SPEEDUP_GATE` is a *measurement*
  against a threshold and therefore a coin-flip on machine load, which makes any inference drawn
  from where the suite stops unreliable. Per-file runs (as Phase 13 does, and as this plan did)
  remain the usable signal.
- Advisory only. Nothing in 12-04 is gated on the full suite, and both of this plan's files are
  green under `runtests.jl` and standalone.

## Known Stubs

None. Both files are complete implementations; nothing in this plan is deferred to a later one.
`P12_ZSCORE_ARM` records `:per_row` as the arm in force — that is R-3's DEFAULT, not a
placeholder; `:shared_scalar` is a declared alternative whose use would be a recorded deviation.

## Threat Flags

None. No new network, auth, file-access or schema surface: the plan adds a pure in-memory array
re-view, an in-memory augmentation, two network builders and a test file. Nothing is written to
disk (confirmed — `train` produced no artifacts), nothing is read beyond the frozen Tier-1
constants and the sibling lattice module.

Threat-register coverage: T-12-14 (testsets 1-2 + the executed transposed falsifier), T-12-17
(testset 3's `< 100_000` and `< 200_000` ceilings plus the exact `== 70_872`), T-12-18 (testset
3's comment-stripped source assertion), T-12-19 (testset 8 + the `P12_ZSCORE_ARM` constant),
T-12-09 (zero `Pkg.add`; `git diff --quiet` green on `spike/Project.toml`/`Manifest.toml`).

## For the Next Plan

- **The θ row layout is now NAMED, so nobody should re-derive it.** Take `p12_row_c0`,
  `p12_rows_dev` (`2:64`), `p12_rows_nuisance` (`65:71`) and `p12_row_r1` (`72`), or
  `p12_theta_index(:sym)`. A second derivation is a second chance to be off by one, exactly as
  a second `radial_basis` derivation would be a second chance to transpose.
- **ρ_true is NOT a θ row and testset 6 asserts it.** Any plan that reports a scalar ρ must
  label it DERIVED (`ghat(quantile(MU_PRIOR, Φ(c0 / G)))`, R-1).
- **`P12_ZSCORE_ARM` must be persisted into every Phase-12 artifact** (T-12-19). It is a plain
  `Symbol`; add it to each `jldsave` alongside `schema_version`.
- **Standardize FIRST, then reshape.** `reshape_summary`'s docstring carries the contract but,
  unlike `augment_input`, it cannot assert it at runtime — a 128-row matrix looks the same
  standardized or not. 12-11 (`train_p12_npe.jl`) and 12-16 (coverage) own that ordering.
- **`augment_mask!` is TRAINING-TIME ONLY.** If leave-region-out coverage later comes out wildly
  off nominal for BOTH the spatial arm and the D-10 ablation, suspect this augmentation before
  suspecting the prior — a shared failure across two arms that differ only in spatial
  correlation cannot be caused by the thing that differs between them.
- **DEF-12-01 blocks any plan that wants a green full-suite run.** One line in `harness.jl:51`
  (re-guard on `:SBC_M`) is the preferred fix; do NOT fix it by editing `p12_consts.jl` or by
  moving the Phase-12 include.
- **Do not read the full-suite exit code as a Phase-12 signal.** Read the `P12-SUITE-RAN`
  markers and the per-testset summaries, as 12-03 already advised — and note that the abort
  point is not stable between runs.

## Self-Check: PASSED

All three files verified present on disk: `spike/npe/p12_architecture.jl` (399 lines),
`spike/test/test_p12_architecture.jl` (334 lines),
`.planning/phases/12-spatial-colocalization-map/deferred-items.md` (69 lines). Both task commits
verified present AND verified ancestors of `HEAD` (`git merge-base --is-ancestor`): `0b36c16`
(architecture), `17e798e` (tests). `spike/Project.toml`, `spike/Manifest.toml`, `src/` and
`spike/validation/p12_consts.jl` byte-unchanged; `.planning/STATE.md` and `.planning/ROADMAP.md`
untouched; `spike/data/cache/p11` intact at 54 MB; `artifacts/amended_v2/grid_8` untouched. Every
commit used an explicit pathspec; no `git stash`, `rebase`, `amend`, `reset`, `clean` or
`checkout --` was run, and no `index.lock` race occurred despite the concurrent Phase-13 executor
committing to the same branch during this session.
