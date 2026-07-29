---
phase: 12-spatial-colocalization-map
plan: 12
status: complete
subsystem: spike/p12 (result type)
tags: [result-type, delta-rho, structural-identity, sentinel-discipline, decoupling]
requires:
  - "spike/npe/p12_architecture.jl (12-04) — P12_ZSCORE_ARM, the theta row map"
  - "spike/simulator/p12_lattice.jl (12-03) — p12_idct_vec, p12_dct_order"
  - "spike/validation/p12_consts.jl (12-01) — P12_REPO_ROOT, P12_G"
  - "spike/validation/sbc.jl — CalibrationResult (included, never redefined)"
  - "src/results.jl — AbstractColocResult, read-only include"
provides:
  - "SpatialColocResultSpike <: AbstractColocResult — three maps + three sd maps + three OOD maps"
  - "delta_rho_map / uncertainty_map — the two names the src/ sketch specifies"
  - "p12_region_maps / p12_mark_unscorable! — the working bundle and its single sentinel writer"
  - "p12_result_meta / p12_ood_checked / p12_region_field"
  - "P12_SPATIAL_SENTINEL, P12_IDENTITY_ATOL, P12_RESULT_META_KEYS"
affects:
  - "12-16 (p12_coloc_map writes into this type)"
  - "12-19 (the real arm builds the delta_rho_available = false mode)"
tech-stack:
  added: []
  patterns:
    - "spike/p13/result.jl mirrored structurally: executable subtype in spike/, sketch untouched"
    - "guard block covers consts + struct only; functions at top level so Julia 1.12 keeps docstrings"
key-files:
  created:
    - spike/p12/result.jl
  modified:
    - spike/test/test_p12_result.jl
decisions:
  - "The Delta-rho guard is a STRUCTURAL IDENTITY, elementwise; no range check was reinstated"
  - "p12_region_maps initialises FAIL-CLOSED (every region unscorable) — plan left this open"
  - "bayes_factor routes to _iface_error rather than fabricating an evidence quantity"
metrics:
  duration_min: 42
  tasks: 2
  files: 2
  assertions: 120
  completed: 2026-07-29
---

# Phase 12 Plan 12: Per-Region Result Type Summary

A per-region Delta-rho map **with** a per-region uncertainty map now exists as a genuine
`AbstractColocResult` subtype in `spike/p12/result.jl`, carrying all three maps
(`region_delta_rho` primary, plus `region_rho_sample` and `region_rho_control`), whose constructor
enforces the Delta-rho semantics as a **structural identity** rather than a range and makes an
unscorable region impossible to confuse with a measured zero — with `src/results.jl` byte-unchanged
and its `SpatialColocResult` identifier still unclaimed.

## What was built

| Task | Name | Commit | Files |
|---|---|---|---|
| 1 | `spike/p12/result.jl` — subtype, accessors, sentinel discipline | `82ddc8b` | `spike/p12/result.jl` (698 lines, new) |
| 2 | `test_p12_result.jl` — eleven testsets replacing the 12-01 scaffold | `299237f` | `spike/test/test_p12_result.jl` (469 lines) |
| — | Deviation: phase-prefix the test helpers (Rule 2) | `160e219` | `spike/test/test_p12_result.jl` |

**Observed verification results** (each command was actually run; results are reported as observed,
not paraphrased from the plan):

- Task 1 `<verify>`: the plan's one-liner exits 0 and prints `result-type-ok`. It constructs the
  thirteen-argument three-map struct, asserts the identity and the union rule, asserts a violated
  identity throws `ArgumentError`, and constructs the single-stack mode.
- Task 2 `<verify>`: `julia --project=spike -e 'using Test; include("spike/test/test_p12_result.jl")'`
  → **120 passed, 0 failed, 0 errored**, all eleven testsets present.
- In-process integration: `include("spike/test/test_p12_suite.jl")` runs the whole Phase-12
  aggregator in one process and reaches the end. The result testsets are green there
  (120/120), and the **live 12-05 decoupling test is 24/24 with the new `spike/p12/result.jl`
  present in its scan surface** — so the new file passes the sealed-holdout token scan and the
  `src/` byte-equality assertions rather than merely not being looked at.

**Not run, and stated rather than implied:** the plan's third `<verification>` line,
`julia --project=. -e 'using Pkg; Pkg.test()'` ("unchanged"). It was not executed. What *was*
verified is the input it depends on: `src`, root `Project.toml` and root `Manifest.toml` are all
byte-unchanged against HEAD and carry no untracked files, and nothing in this plan is imported by
the shipped package — so the root suite is unaffected by construction. That is an argument, not a
measurement, and it is labelled as one. (`MEMORY.md` separately records that the root ADVI path may
not run at all on the current Manifest, which is a pre-existing condition this plan does not touch.)

## The Delta-rho ruling, as implemented

`region_delta_rho` is the primary named deliverable and the two single-stack maps are additionally
exposed, at zero extra forward passes. The guard is the elementwise identity

    region_delta_rho ≈ region_rho_sample − region_rho_control      (atol = P12_IDENTITY_ATOL = 1e-9)

scoped to regions where none of the three is a sentinel. **No `[-2, 2]` range guard was
reinstated**; `[-1, 1]` on the single-stack maps and `[-2, 2]` on the difference survive only as
sanity checks, labelled as such in the code, in the error messages and in the tests. Delta-rho is
NOT a theta row — it is nowhere near `P12_THETA_ROWS`, which stays at 72.

### The identity guard was proved capable of failing

Required by the brief, and done three ways rather than asserted:

1. **The mutant matrix.** A scratch copy of `result.jl` with the identity block deleted was run
   against two mislabels. Both **constructed successfully**: the sample map handed over as the
   difference, and the difference negated. With the identity block present, both throw
   `ArgumentError` whose message names `STRUCTURAL IDENTITY FAILED`. The scratch copy was deleted;
   `git diff --quiet HEAD -- spike/p12/result.jl` confirms the shipped file is byte-unchanged.
2. **The range's blindness is on the record.** In the same breath the test asserts that the
   mislabelled map lies inside `[-2, 2]` *and* inside `[-1, 1]` — so the retired bound would have
   admitted it silently. That is the ruling's argument made executable rather than restated.
3. **Testset 3's three falsifiers were each proved load-bearing**, one constructor check at a time,
   in three scratch mutants. The result is a clean diagonal — each falsifier goes silent for its own
   check and for no other:

   | scratch mutation | (a) finite sd under `ood` | (b) `1e-9` Δρ under `ood` | (c) `NaN` sd under `!ood` |
   |---|---|---|---|
   | sd-pairing `NaN` branch deleted | **CONSTRUCTS** | throws | throws |
   | sentinel-pairing deleted | throws | **CONSTRUCTS** | throws |
   | sd-pairing positive branch deleted | throws | throws | **CONSTRUCTS** |

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 2 — Missing critical safeguard] The test helper `_res` would silently overwrite Phase 13's**

- **Found during:** Task 2, after the testsets were green.
- **Issue:** `spike/test/test_p13_result.jl:77` binds a helper named `_res`, declared
  `_res(; coloc = 2.5, ...)`. The plan's implied fixture shape, `_res(maps = default; ...)`,
  generates a zero-positional method with the **identical** signature `Tuple{typeof(_res)}` — so in
  the spike lane's flat top-level namespace the later include silently replaces the earlier one.
  Verified empirically in a scratch process, not inferred. Today's include order
  (`test_p12_suite.jl` runs first) hides it, which is precisely what makes it worth removing: the
  defect would surface only on a suite reorder and would then read as a Phase-13 failure. This is
  the same collision class `test_p13_result.jl:58-66` already guards `P13_REPO_ROOT` against.
- **Fix:** `_res → _p12res`, `_meta → _p12meta`, `_measured_maps → _p12_measured_maps`,
  `RES_FIX → P12_RES_FIX`, with the reasoning recorded in the fixtures banner. 120/120 unchanged.
- **Commit:** `160e219`

**2. [Rule 1 — Bug in my own first draft] A float-literal comparison that was simply wrong**

- **Found during:** Task 2, first run — one failing assertion out of 118.
- **Issue:** I asserted `region_delta_rho[1] == 0.4` on a fixture built as `0.6 - 0.2`, which is
  `0.39999999999999997` in binary floating point.
- **Fix:** compare with `≈ 0.6 - 0.2`, and assert the two single-stack values exactly (they are
  filled, not computed), with the reason in a comment.
- **Commit:** folded into `299237f`.

### Decisions the plan left open

**`p12_region_maps` initial state — FAIL-CLOSED.** The plan specifies the bundle's contents and
order but not its initial values. It is initialised fully unscorable (every rho map at
`P12_SPATIAL_SENTINEL`, every sd `NaN`, every flag `true`), so a write path that skips a region
yields an honest sentinel rather than a confident zero. The initial bundle is itself a legal
construction — it is exactly the no-control-stack shape with the sample stack unmeasured too.
Recorded in the docstring.

**`bayes_factor` refuses rather than fabricates.** The plan requires "the four shared accessors,
matching the analog's signatures exactly" without saying what `bayes_factor` should return. This
phase computes no evidence quantity, so it routes to `_iface_error`, mirroring `spike/p13/result.jl`'s
treatment of the Delta-rho *it* does not carry. The test asserts the refusal and its message.

## Plan/reality discrepancies found (reported, none blocking)

Every wave so far has found a defect; these are this wave's, and all are citation drift rather than
a false premise. None changed what was built, and none warranted `status: blocked`.

| Cited | Plan says | Measured | Material? |
|---|---|---|---|
| `src/results.jl:181` for `region_delta_rho` | "the primary named deliverable, per `src/results.jl:181`" | `:181` is the `struct SpatialColocResult <: AbstractColocResult` line; `region_delta_rho` is at **`:183`**. The block is `:180-191`, which the plan states correctly elsewhere. | No — the sketch was read in full and mirrored as read. |
| `test_p13_result.jl:185-219` for the "src untouched" testset | Task 2 item 10 | That testset is **`:190-213`**; `:185-219` over-spans into the preceding testset and to EOF. The companion citation `:215-217` (CPU-only) is exact. | No — the testset was copied in shape from the real range. |
| `spike/p13/result.jl:150-175` for the validating inner constructor | Task 1 read_first | The inner constructor is **`:150-163`**; `:166` is the guard block's `end`. | No. |

The plan's load-bearing citations were checked and are **correct**: `src/amortized/infer.jl:117-121`
is the MC difference (the frozen amendment §5's citation, which I mirrored, rather than the brief's
`:108`, which is the docstring signature line); `local_map.jl:42-53`, `:56-79` and `:88-91` are the
sentinel, the "no per-region uncertainty" docstring and the NOT-CHECKED rule respectively;
`p13/result.jl:179` and `:337` are the keyed-NamedTuple guard and the meta builder;
`sbc.jl:57-155` does contain `CalibrationResult`, `_bin_calibration` and `sbc_traffic_light`.

## Assumption Drift (advisory)

**The `delta_rho(r)` accessor returns a rho, not a Delta-rho.** Planned and implemented exactly as
the plan specifies (`return r.meta.derived_rho`, labelled DERIVED in the docstring's first
sentence). What turned out true on contact with the code is that the *name* still says "delta_rho"
while the *value* is a single-stack derived rho — a mismatch a docstring can warn about but cannot
remove, and the shared interface's method name is not this plan's to change. Recorded because a
reader who skips the docstring gets a plausible number under a wrong name, and because it is the
kind of thing a future productionization phase should fix at the interface rather than in a
docstring. Advisory only; nothing was gated on it and nothing was changed.

## Proactive audit: include-guard poisoning

A concurrent commit (`6d3bc8a`, Phase-13 lane) recorded a spike-wide latent defect — 29 of 220
include guards keyed on a sentinel declared by more than one file, so the guard short-circuits and
the included body silently never runs. Since this plan writes **five** new guards, each was audited
against that defect class before this summary was written. All five sentinels are singly-owned:

| Sentinel | Sole declaring file |
|---|---|
| `AbstractColocResult` | `src/results.jl` |
| `_bin_calibration` | `spike/validation/sbc.jl` |
| `P12_ZSCORE_ARM` | `spike/npe/p12_architecture.jl` |
| `p12_idx` | `spike/simulator/p12_lattice.jl` |
| `SpatialColocResultSpike` | `spike/p12/result.jl` |

Relatedly, and recorded in the file: `sbc.jl` is included **before** `p12_architecture.jl` (hence
before `p12_consts.jl`), because DEF-12-02/DEF-12-03 say a process that loads the Phase-12
pre-registration first can silently skip a foreign one — `p12_consts.jl` binds `:P11_DEV_SEED` and
`:P13_DEV_SEED`. Loading the foreign surfaces first is the standing workaround, and the plan's
prescribed include order already had it right.

## Known Stubs

None. Every function in `spike/p12/result.jl` is fully implemented; nothing returns a placeholder,
an empty container or a hardcoded constant standing in for a computation. The one deliberate
refusal, `bayes_factor`, raises the documented interface error rather than returning a stub value —
which is the opposite of a stub.

## Threat Flags

None. The two files introduce no network endpoint, no auth path, no file write and no schema at a
trust boundary. The only file access is read-only: `src/results.jl` (`include`) and
`src/amortized/local_map.jl` (`read`, as text). The plan's threat register entries T-12-20, T-12-43,
T-12-44, T-12-45, T-12-46 and T-12-47 are each mitigated and each has a corresponding testset
(10, 3, 8, 7, 4 and 1 respectively).

## Read-only surfaces, verified byte-unchanged after the run

`git diff --quiet HEAD` succeeds and `git status --porcelain` is empty for every one of:
`src` (including the untracked case), `spike/Project.toml`, `spike/Manifest.toml`, root
`Project.toml`, root `Manifest.toml`, `corpus`, `spike/data/cache/p11`, `.planning/STATE.md` and
`.planning/ROADMAP.md`. `spike/data/cache/p12` does not exist and was not created. Every commit was
made with an explicit pathspec (`git commit -- <path>`), as the concurrency rule requires; no
`git add -A`, no `git stash`, no bare commit. No index-lock race occurred, though a concurrent
Phase-13 commit (`6d3bc8a`, `.planning/STATE.md` only) did land between Task 1 and Task 2 — it
touches nothing in this plan's surface and `82ddc8b` remains an ancestor of `HEAD`.

## Self-Check: PASSED

- `spike/p12/result.jl` — FOUND
- `spike/test/test_p12_result.jl` — FOUND
- Commit `82ddc8b` — FOUND
- Commit `299237f` — FOUND
- Commit `160e219` — FOUND
