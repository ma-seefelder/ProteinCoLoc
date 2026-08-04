---
phase: 15-calibration-operating-envelope-and-ci-gate
plan: 03
subsystem: gate-sweep-axes
tags: [misspecification, ladders, prior-boundary, matched-pairs, ood, frozen-tables]
requires:
  - test/gate/misspec.jl (the generator contract, OOD_FAMILIES, misspec_simulator — read, never written)
  - test/gate/p15_consts.jl (P15_SPILLOVER_LADDER / P15_SHIFT_LADDER / P15_AF_LADDER, the mechanism typing, the frozen-file pins)
  - src/amortized/simulator.jl (sample_prior's θ fields, simulate_pair's range guards and rng-consuming stages)
provides:
  - test/gate/p15_misspec.jl (misspec_spillover / misspec_registration / misspec_autofluorescence)
  - P15_FAMILIES — the seven-axis table, a merge of OOD_FAMILIES, never a mutation
  - p15_axis_mechanism / p15_axis_boundary_rung accessors
  - test/gate/test_p15_families.jl (158 standalone / 159 in-suite assertions)
  - the MEASURED per-axis ladder magnitudes the domain map's honesty section needs
affects:
  - 15-05 (the sweep runner: consumes P15_FAMILIES through the existing `families =` keyword; must read the matched-pairs limit below before framing any cross-rung delta)
  - 15-07 (the rung-0 ECE anchor; rung 0 here is a proven byte-identity with the unmisspecified simulator)
  - 15-09 (the domain map: the mechanism labels and the effective-magnitude table are its inputs)
tech-stack:
  added: []
  patterns:
    - pure θ-override generators that delegate to simulate_pair rather than forking it
    - local merge of a frozen table, passed through an existing keyword
    - guarded-include preamble for standalone-loadable gate files
    - dependency-free NIST-verified SHA-256 in a test (Project.toml is pinned frozen)
key-files:
  created:
    - test/gate/p15_misspec.jl
    - test/gate/test_p15_families.jl
  modified: []
decisions:
  - "The plan's 'identical rng consumption across level 0..5' is MEASURABLY FALSE — stage 7's Poisson consumption is rate-dependent, so matched pairs cover draw 1 and no further; the test now asserts the property that is true (exact pure-θ-override equality) and asserts the divergence so the false claim cannot re-enter"
  - "The 45° shift reproduces the ladder norm to within ONE ULP, not exactly; asserted with a tolerance and stated in the docstring rather than claimed exact"
  - "Monotonicity uses 24 prior draws, not the planned 8: at 8 the spillover rung-1→rung-2 step was measured non-monotone at 2 of 3 fixture seeds"
  - "The misspec.jl pin is checked with a dependency-free LF-normalized SHA-256, inheriting 15-01's two findings (SHA is unavailable under Pkg.test(); .gitattributes makes raw bytes platform-dependent)"
metrics:
  duration_minutes: 38
  tasks: 2
  files_created: 2
  files_modified: 0
  tests_added: 158
  completed: 2026-08-04
---

# Phase 15 Plan 03: The Three Missing Sweep Axes Summary

Built `spillover`, `registration` and `autofluorescence` as pure θ-overrides on the unmodified
forward model, so all seven axes are addressable through one `P15_FAMILIES` NamedTuple while
`test/gate/misspec.jl` stays byte-unchanged — and, in the course of testing them, **measured two of
the plan's stated properties to be false and corrected them rather than asserting them**.

## What Was Built

### `test/gate/p15_misspec.jl` (210 lines)

Three generators, each satisfying the shipped `(rng, θ; imsize, level, G)` contract, each opening
with `_guard_misspec_imsize` unconditionally (including at level 0), and each **delegating to
`ProteinCoLoc.simulate_pair`** rather than reimplementing any part of it — a forked forward model
would make the rung-0 identity unverifiable by construction, which is the one property these axes
exist to have.

| Axis | Override | Boundary rung | Ladder endpoint is |
|------|----------|---------------|--------------------|
| `spillover` | `θ.spillover = P15_SPILLOVER_LADDER[k]` | 2 (`= 0.2 = maximum(SPILLOVER_PRIOR)`) | the **simulator's own guard** (`simulator.jl:232`, `0 ≤ spillover ≤ 1`) |
| `registration` | `θ.shift_dx = θ.shift_dy = P15_SHIFT_LADDER[k]/√2` | 2 (`= 1.0 px = maximum(SHIFT_PRIOR)`) | a stated choice (no range guard exists) |
| `autofluorescence` | `θ.autofluorescence = P15_AF_LADDER[k]` | 2 (`= 0.1 = maximum(AUTOFLUORESCENCE_PRIOR)`) | a stated choice (only `≥ 0` is guarded) |

`const P15_FAMILIES = merge(OOD_FAMILIES, (...))` — seven entries with the shipped four **first**,
so `Tuple(keys(P15_FAMILIES)) == P15_AXES` and the existing report ordering is preserved. There is
no assignment to `OOD_FAMILIES` anywhere in the file (`grep -nE '^\s*(const\s+)?OOD_FAMILIES\s*='`
returns nothing) and `git diff --quiet HEAD -- test/gate/misspec.jl` exits **0**.

Plus `p15_axis_mechanism` / `p15_axis_boundary_rung`, each documenting that `nothing` means "this
axis has no prior support to mark", never "unmeasured".

**No level-0 entry point for the four shipped families**, deliberately: rung 0 is the
unmisspecified model, so it is the same measurement for every axis and 15-05 measures it once
through `default_simulator()`.

### `test/gate/test_p15_families.jl` (454 lines, 158 assertions standalone / 159 in-suite)

Seven testsets. No net is loaded, no inference runs; the whole file is **20.9 s wall-clock warm**
(9.8 s of testset time), against a 60 s acceptance budget.

## The measured ladder magnitudes (the deliverable `<output>` asks for)

Mean absolute deviation of the final image from the axis's **own rung-0 image** at the same seed
and the same θ, averaged over 24 prior draws at 128×128:

| Axis | rung 1 | rung 2 (**prior boundary**) | rung 3 | rung 4 | rung 5 |
|------|--------|-----------------------------|--------|--------|--------|
| `spillover` | 0.24710 | **0.26147** | 0.33706 | 0.47759 | 0.72997 |
| `registration` | 0.14592 | **0.15310** | 0.18202 | 0.23879 | 0.33880 |
| `autofluorescence` | 0.26650 | **0.28401** | 0.37845 | 0.70348 | 1.47986 |

**What this number contains, stated so the domain map does not over-read it.** It is the deviation
of the *final* image, so it includes stage 7's resampled photon noise. Once an override moves the
Poisson rate far enough for the draw to decorrelate, a decorrelation component (~0.16 per channel
at `GAIN = 50`) enters and never leaves. That is why rung 1 already sits near 0.25 on the two
offset axes rather than near 0: **these are total image-change magnitudes, not the isolated effect
of the parameter.** The rung-1→rung-2 step is correspondingly small on `spillover` (0.0144) and
`registration` (0.0072) — both steps stay inside or at the prior boundary, where the parameter
change is small relative to the resampling.

Measured, separately, that this is not an artificial floor: a physically negligible override
(`autofluorescence + 1e-12`) produces a **byte-identical** image, mean |Δ| exactly **0.0**.

## Verification (observed, not assumed)

| Check | Command | Observed |
|-------|---------|----------|
| Task 1 verify (7 vs 4 keys, exact key set) | the plan's `module M; include(...)` one-liner | printed `ok` |
| Task 2 verify, standalone | `julia --project=. --threads=auto test/gate/test_p15_families.jl` | **158 pass, 0 fail, 0 error**, 9.8 s |
| Task 2 wall-clock | `time` around the same command | **20.9 s** (budget 60 s) |
| Full suite (the `Pkg.test()` sandbox path) | `julia --project=. --threads=auto -e 'using Pkg; Pkg.test()'` | `Testing ProteinCoLoc tests passed`; the Phase-15 axes testset reports **159 pass** (the extra one is the cross-check against `test_p15_consts.jl`'s digest, which only exists when both files are loaded) |
| 15-01's testset still green | `julia --project=. --threads=auto test/gate/test_p15_consts.jl` | 8 testsets, **303 pass, 0 fail** (50 + 6 + 9 + 103 + 34 + 41 + 41 + 19) |
| Frozen files untouched | `git status --porcelain test/gate/misspec.jl test/gate/sbc.jl test/gate/harness.jl test/runtests.jl` | **empty** |
| `misspec.jl` unmodified | `git diff --quiet HEAD -- test/gate/misspec.jl` | exit **0** |
| No `OOD_FAMILIES` assignment | `grep -nE '^\s*(const\s+)?OOD_FAMILIES\s*=' test/gate/p15_misspec.jl` | no match |
| `merge(OOD_FAMILIES` present | `grep -c` | **2** (the code and its header explanation) |
| No accidental deletions | `git diff --diff-filter=D --name-only HEAD~3 HEAD` | **empty** |

RED was observed before GREEN: the first run of `test_p15_families.jl` exited **1** with
`SystemError: opening file ... p15_misspec.jl: No such file or directory`, and that failing state is
commit `9ad5d7d`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] "Identical rng consumption across `level = 0..5`" is measurably FALSE**
- **Found during:** Task 1, probing the plan's stated behaviour before asserting it
- **Issue:** The plan's Task-1 `<behavior>` and Task-2 testset 3 both require asserting that the
  number of rng values each generator consumes is identical across levels, on the stated ground
  that "stages 4 and 6 consume no rng". Stages 4/5/6 indeed consume none — but **stage 7 does**, and
  `rand(rng, Poisson(GAIN·intensity))` consumes a **rate-dependent** number of values, while the
  override is precisely a change of rate. Measured post-generator stream state versus rung 0, over
  six fixture θ:

  | Axis | rungs that diverge from rung 0 |
  |------|--------------------------------|
  | `registration` | **5 of 5, at every one of the six θ** |
  | `autofluorescence` | **5 of 5, at every one of the six θ** |
  | `spillover` | **θ-dependent: 0, 5, 5, 1, 5, 2 of 5** across the six θ, in seed order |

  Asserting the plan's version would have produced a test that fails (it did — 1 failure on the
  first GREEN run), and worse, the *claim* is what a sweep report would have leaned on.
- **Consequence, which matters more than the test:** `sbc_gate` threads **one** rng through all M
  draws (`harness.jl:188-204`), so two rungs agree on draw 1 and then see **different θ\*
  sequences from draw 2 on**. The matched-pairs variance reduction the plan and research both cite
  covers exactly one draw. Cross-rung ECEs are neither independent replicates nor exactly paired
  ones, and 15-05/15-09 must claim neither.
- **Fix:** Replaced the false assertion with the two that are true and load-bearing: (a) each rung
  is **exactly** `simulate_pair` on θ with one field overridden by the frozen ladder, on the same
  stream — which is the entire content of the matched-pairs design *per draw*, and pins the wiring
  just as tightly; (b) the divergence itself is asserted where it is robust (`registration`,
  `autofluorescence`), **reported and not asserted** where it is θ-dependent (`spillover`), and
  **shown** on a second fixture θ where spillover does diverge — so "identical consumption" is
  falsified for all three axes rather than only two.
- **Files modified:** `test/gate/test_p15_families.jl`, `test/gate/p15_misspec.jl` (header)
- **Commit:** `8e9e32a`

**2. [Rule 1 - Bug] The 45° shift norm is exact to one ULP, not exact**
- **Found during:** Task 1
- **Issue:** The plan (and the research it quotes) states the override makes `|(dx, dy)|` equal the
  ladder value **EXACTLY**. Measured: `r/√2` followed by `hypot` is not an exact Float64 round trip
  — `0.5 → 0.49999999999999994`, `1.0 → 0.9999999999999999`, and so on for all five rungs.
- **Fix:** Asserted with `rtol = 4eps(Float64)` and stated in the `misspec_registration` docstring
  that the norm reproduces the ladder value "to within one ulp", with the measured example. A
  sub-ulp difference in a sub-pixel shift is physically nothing; an `==` in the test would have
  been a false claim about the code, which is the class of defect this phase's discipline exists to
  prevent.
- **Commit:** `8e9e32a`

**3. [Rule 1 - Bug] 8 prior draws is not enough to assert monotonicity; 24 is**
- **Found during:** Task 2
- **Issue:** The plan specifies averaging the magnitude over **8** prior draws and asserting the
  five-element sequence strictly increasing. Measured at 8 draws across three independent fixture
  seeds, `spillover` was **non-monotone at 2 of 3** — the rung-1→rung-2 step (0.1 → 0.2, both
  inside the `U(0, 0.2)` prior) is smaller than the across-draw variability. At 24 draws it held at
  **3 of 3**.
- **Fix:** `nrep = 24`, with the measurement written next to it. The test is deterministic at any
  draw count, so this is not a flakiness fix — it is the difference between asserting a real
  property and asserting a lucky seed.
- **Commit:** `9ef05ae`

**4. [Rule 3 - Blocking] `import SHA` would abort the suite; the pin is checked LF-normalized**
- **Found during:** Task 2
- **Issue:** The plan's Task-2 action specifies `import SHA` and
  `bytes2hex(SHA.sha256(read("test/gate/misspec.jl")))`. Both halves are wrong here, for reasons
  15-01 already recorded and pinned into the repository: `SHA` is unavailable inside `Pkg.test()`'s
  sandbox (it is not in `[deps]`, and adding it to `[extras]` would edit the `Project.toml` this
  phase pins as frozen), and the pin is over **LF-normalized** bytes because `.gitattributes`
  `* text=auto` makes raw working-tree bytes CRLF on Windows and LF on the Linux CI this phase
  exists to build. A raw-bytes comparison would have failed on CI for a reason unrelated to the
  file's contents.
- **Fix:** A dependency-free SHA-256 in the test file, verified against **four** published FIPS
  180-4 / NIST vectors covering the empty, single-block, two-block and three-block padding paths,
  applied after `replace(s, "\r\n" => "\n")`.
- **Duplication, and why it was chosen over a shared helper:** `test_p15_consts.jl` carries the same
  ~50 lines. A shared `test/gate/p15_sha256.jl` would be a **sixth** file outside this plan's
  declared surface, and `runtests.jl`'s Phase-15 wiring is a fixed five-file list of `test_p15_*.jl`
  names that 15-01 built specifically so later plans never edit that file again. Instead the
  duplication was turned into a feature: when both files are loaded (every `runtests.jl` run) the
  two **independent implementations are asserted to agree** on `misspec.jl`, so a bug in either
  cannot pass unnoticed. That cross-check is the 159th assertion.
- **Commit:** `9ef05ae`

### Mechanics decisions taken autonomously

**5. The TDD RED for Task 1 was written into Task 2's file.** Task 1 carries `tdd="true"` with a
`<behavior>` block, but the plan assigns the only test file to Task 2. Rather than create a
throwaway test, the RED commit (`9ad5d7d`) put Task 1's behaviours — contract, guard, rung-0
identity, θ-override, monotonicity — into `test/gate/test_p15_families.jl` and was observed to fail
with the file-not-found error; Task 2 then appended testsets 5-7. Gate sequence in `git log`:
`test(15-03)` → `feat(15-03)` → `test(15-03)`.

**6. Fixture seeds are literals, deliberately not `prod_rng`.** `P15FAM_STREAM_SEED` and
`P15FAM_THETA_SEED` are plain constants so that running the test suite can never consume, advance
or otherwise entangle the pre-registered production stream (D-11). Nothing in this file is a
Phase-15 measurement — no bar is scored and no rung is decided — so these are unit-test fixtures in
the ordinary sense, and are labelled as such in the file.

## Assumption Drift (advisory)

**A. "Stages 4 and 6 consume no rng" was treated as equivalent to "the generator's rng consumption
is level-independent".** *Planned* (plan `<interfaces>`, research Pattern 2 point 2): the two are
stated as the same fact, and the matched-pairs claim rests on it. *Actual:* the first is true and
the second is false, because stage **7** runs after the overridden stages and its consumption
depends on the intensities they produced. *Why it matters:* the same inference — "this parameter is
consumed by a non-sampling stage, therefore the stream is unaffected" — is available for any future
θ-override axis, and it will be wrong for the same reason every time. The general rule is that any
override that changes a **rate** changes the stream, regardless of which stage reads it.

**B. The magnitude of a ladder rung was assumed to be readable from the image difference.**
*Planned:* the monotonicity testset is framed as reporting "the ladders' effective magnitudes".
*Actual:* the quantity is dominated at the low rungs by resampled photon noise, not by the
parameter — rung 1 sits near 0.25 on the offset axes where the parameter's own contribution is a
fraction of that. *Why it matters:* 15-09 should present these as *total image change*, and should
not compare them across axes as if they were parameter-effect sizes.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no schema at a trust boundary, and no file
access beyond reading two already-committed repository files inside a test.

`T-15-09` (tampering with the frozen `OOD_FAMILIES`) is mitigated as planned and then some: the new
axes enter only through `merge`, the key set and the four function identities are asserted, the
file's LF-normalized SHA-256 is checked against `P15_FROZEN_FILE_SHA256.misspec` on every test run
with a failure message that states the consequence, and a second independent digest implementation
cross-checks the first whenever both are loaded. `T-15-15` (ladder literals drifting from the
prior) is mitigated by indexing the frozen ladders rather than re-typing magnitudes, plus three
assertions that rung 2 **is** `maximum(<PRIOR>)` read from `ProteinCoLoc` itself.

## Known Stubs

None.

## Notes for Later Plans

- **15-05 (sweep runner): the matched-pairs guarantee covers draw 1 only.** If the sweep's design
  depends on cross-rung noise reduction beyond that, it needs an explicit mechanism — e.g. deriving
  a per-draw substream keyed on (seed, draw index) so every rung sees the same θ\* sequence. That is
  a design decision with a real cost and it was deliberately **not** taken here: re-syncing the
  stream inside a Phase-15 simulator wrapper would make the misspecified arms' draw sequence differ
  from `default_simulator()`'s, which is exactly the shared rung-0 anchor they are compared against.
  Rule on it explicitly rather than inheriting it.
- **The `families =` keyword is the only supported route.** Pass `P15_FAMILIES` to `gate_ood_roc`;
  never append to `OOD_FAMILIES`.
- **`SBC_IMSIZE` is `:mixture` under `p15_consts.jl`**, not a tuple, so the generators' declared
  `imsize = SBC_IMSIZE` default is decorative exactly as it is in `misspec.jl`. Every call site must
  pass `imsize` explicitly — the guard will throw on the symbol otherwise.
- **A dependency-free digest is now the phase's pattern.** Any later Phase-15 test needing a hash
  should reuse `_p15_lf_sha256` (loaded from `test_p15_consts.jl` under `runtests.jl`) or copy the
  same ~50 lines with a distinct name; `import SHA` will abort the suite in the sandbox.

## Self-Check: PASSED

- `test/gate/p15_misspec.jl` — FOUND (210 lines, ≥ 120 required)
- `test/gate/test_p15_families.jl` — FOUND (454 lines, ≥ 110 required)
- `.planning/phases/15-calibration-operating-envelope-and-ci-gate/15-03-SUMMARY.md` — FOUND
- commits `9ad5d7d`, `8e9e32a`, `9ef05ae` — all present in `git log`
- `git status --porcelain` on `misspec.jl` / `sbc.jl` / `harness.jl` / `runtests.jl` — empty
- Every number in this summary was printed by a command run in this session; none is transcribed
  from the plan or the research document.
