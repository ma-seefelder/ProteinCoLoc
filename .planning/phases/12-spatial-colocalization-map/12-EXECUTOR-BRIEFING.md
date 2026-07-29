# Phase 12 — Standing Executor Briefing

Written and maintained by the phase-12 orchestrator. **Every phase-12 executor reads this file
before touching anything.** It carries what is true across all plans, so each dispatch does not have
to restate it. Facts here were measured in this execution run, not assumed.

---

## 1. Execution mode — NO WORKTREES

You run **sequentially on the MAIN working tree** at `C:\Users\Manuel\Documents\GitHub\ProteinCoLoc`.
Deliberate ruling (`3e31d8e`), for two reasons:

- 12-09's training pool lands in gitignored `spike/data/cache/p12/` under an append-only 150-minute
  ceiling, and is consumed by 12-11/13/14/17 in *later waves*. A between-wave worktree cleanup would
  destroy it — which is exactly how Phase 11 lost 54 MB and ~56 min.
- `spike/data/cache/p11` (54 MB, read-only, needed by 12-15 and 12-17) is gitignored and therefore
  **does not exist inside a fresh worktree at all**.

## 2. Concurrency — you are not alone on this branch

A **phase-13 executor** and **the user** both commit to `gsd/v2.0-milestone` while you work.

- ALWAYS commit with explicit paths: `git commit -m "..." -- <path> <path>`.
  **NEVER** a bare `git commit`, **NEVER** `-a`, **NEVER** `git add -A` / `git add .`.
- **NEVER** `git stash`, `git rebase`, `git commit --amend`, `git reset --hard`, `git checkout --`.
  Revert a mutation by restoring a byte-exact copy you saved yourself.
- Scope every `git status` / `git diff` assertion to the paths you mean
  (`git diff --quiet HEAD -- src`), never a repo-wide clean-tree check — it will see another
  agent's in-flight files and be flaky for reasons unrelated to what you are testing.
- On an `index.lock` or a race: WAIT and retry. Document it. Never force.
- **Do NOT modify `.planning/STATE.md` or `.planning/ROADMAP.md`.** The orchestrator owns those and
  writes them additively. Skip the usual `state.advance-plan` / `roadmap.update-plan-progress` steps.
- Hooks run. **Do NOT pass `--no-verify`.** If a hook fails, fix the cause.
- **Order:** Write SUMMARY.md → commit → *only then* any narration. Nothing between Write and commit.

## 3. Frozen surfaces — read, obey, never edit

| File | Status |
|---|---|
| `spike/validation/p12_consts.jl` | **Tier-1 pre-registration, APPEND-ONLY.** Every threshold you need is here. Do not invent or re-derive one. Only 12-15 and 12-16 may append, under reserved Tier-2 sentinels. |
| `12-SC3-AMENDMENT.md` | **Frozen** success criteria. Cite it; never edit it. |
| `src/` | Untouchable. The spike decoupling constraint is absolute, and 12-05's test is LIVE and will catch an edit. Read `src/results.jl` / `src/amortized/infer.jl` as sketches to mirror. |
| `artifacts/amended_v2/grid_8` | Shipped. Do not touch; do not reopen the Phase-7 ship gate. |
| `spike/Project.toml`, `spike/Manifest.toml` | Byte-unchanged is an asserted invariant. |
| `spike/data/cache/p11/` | 54 MB, READ-ONLY. 12-15 and 12-17 depend on it. Never delete, overwrite, move or regenerate. |
| `corpus/` sealed holdout | Reserved for Phase 16. Do not fetch, do not open. |

## 4. Built so far — reuse, do not reimplement

| From | File |
|---|---|
| 12-01 | `spike/validation/p12_consts.jl`, `spike/test/test_p12_suite.jl`, nine scaffolds |
| 12-03 | `spike/simulator/p12_lattice.jl` — CAR + GP arms, r1→α bisection, per-cell rescale |
| 12-04 | `spike/npe/p12_architecture.jl` — 70,872 params |
| 12-05 | `spike/test/test_p12_decoupling.jl` — LIVE, self-scanning |
| 12-06 | `spike/validation/run_p12_eps_ridge.jl` + report |
| 12-07 | `spike/simulator/p12_prior.jl` |
| 12-12 | `spike/p12/result.jl` — the three maps + the structural identity guard |

Your own `spike/test/test_p12_*.jl` file usually **already exists** as a green pending scaffold from
12-01. **REPLACE its body**; do not create the file. Several are appended to by more than one plan.

## 5. Gate signal, and traps that have already bitten

- `runtests.jl` includes `test_p12_suite.jl` **FIRST** among the phase includes. Do not move it.
- **Per-file runs are the gate signal**, for Phase 12 and Phase 13 alike. A full-suite green is not
  reliably available.
- **The suite's abort point is NOT a stable fact.** It has moved three times during this phase:
  Phase-4 `SPEEDUP_GATE` (measured 84.44 / 88.60 / 89.77 / passing) → `test_sbc.jl` → Phase-13
  consts. Never build a check on "the file that throws is X".
- **Do NOT touch `SPEEDUP_GATE`.** The user paused it in `a494247` — their ruling, not yours.
- **Include-guard poisoning is a REPO-WIDE latent defect, now the user's.** `6d3bc8a` audited 220
  guards over 118 files and found **29 poisoned across eight sentinel families**; 23 guard on a seed
  or prior bound — exactly the names other phases MUST mirror to assert stream disjointness. Fixed
  so far: `:SBC_M` (`fb76b84`). Still live and relevant to you: `p12_consts.jl` binds
  `:P11_DEV_SEED` and `:P13_DEV_SEED`, the sentinels `p11_consts.jl:79` and `p13/consts.jl:100`
  guard their Tier-1 blocks on. **If you load more than one pre-registration in a process, load the
  FOREIGN ones FIRST.** Do not fix this class yourself — it is coupled to an open user ruling.
- **NeuralEstimators v0.2.1 defaults `use_gpu = true`** on `train` / `sampleposterior`. CPU-only is
  this project's baseline: pass `use_gpu = false` **explicitly**.
- **Phase-prefix your test helpers.** A bare `_res`-style helper silently overwrote its Phase-13
  namesake via an identical zero-positional method signature (found by 12-12).

## 6. Report; do not guess

If anything in your plan does not match reality — a cited file/line/symbol that does not exist, a
constant whose stated value and derivation disagree, a **verify that cannot fail**, an instruction
contradicting `p12_consts.jl` or the frozen amendment — **STOP**, write it into SUMMARY.md with
`status: blocked`, and name exactly what disagreed with what. Do not commit a guess.

The plans passed six self-review gates and four independent review passes ending CLEAN, so the prior
that an odd-looking instruction is **deliberate** is high. Report what looks wrong rather than
silently "improving" it. But real defects have surfaced in **every wave so far**:

| Plan | Defect found in execution |
|---|---|
| 12-01 | Plan premise factually false — the include it was positioned against was not the one that throws |
| 12-02 | Three stale line citations, an understated section, a line-wrap that broke an acceptance literal |
| 12-03 | A research table the plan's own estimator does not reproduce |
| 12-04 | A library default (`use_gpu = true`) silently contradicting the plan's own acceptance criterion |
| 12-05 | A testset inert because its target directory is gitignored — reported as 0 tests, not as a pass |
| 12-06 | An `@assert isdir(dir)` that cannot fire, because the resolver creates the directory it checks |
| 12-07 | **A testset specified on the four cells its own mutation does not break — as written it could never fail** |
| 12-12 | A fixture helper whose signature silently overwrote its Phase-13 namesake |

**Prove your assertions can fail.** Mutate the mechanism, watch the assertion go red, revert
byte-exactly. Where a plan's chosen operating point does not falsify a mutation, keep the plan's
assertions, **add** falsifiable ones, and screen at a point where the **same pre-registered bar**
actually bites — 12-07's precedent. **Never move a threshold to make a test bite.**

## 7. Settled decisions — do NOT reopen, do NOT re-derive

- **Δρ: ship ALL THREE maps** (`region_rho_sample`, `region_rho_control`, `region_delta_rho`).
  `region_delta_rho` is the PRIMARY named deliverable. The guard is a **STRUCTURAL IDENTITY**
  `region_delta_rho ≈ region_rho_sample − region_rho_control`, asserted **elementwise** — NOT a
  range check. A `[-2,2]` guard was REPLACED, not tightened; reinstating it is a regression, because
  a true Δρ genuinely spans [-2,2] and no range can distinguish it from a ρ map.
- **Δρ is NOT a θ row.** θ is 72 rows. Δρ is derived, carried as a separate appended rank column.
- **The training pool is SINGLE-STACK.** No paired control draw — it would push datagen past the
  append-only `P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150`. The Δρ* truth is built at SCORING time from
  two independent prior draws.
- **Real-image arm: REPORTED, not gated.** `effective_independent_n = 2`; the 128 region-draws are
  pseudo-replicates. Δρ is NOT scored on real data, and the reason is **EXCHANGEABILITY, not shared
  nuisances**: the simulated pair is two exchangeable draws from one prior (a WITHIN-population
  difference), while the two real anchors are deliberately non-exchangeable (BETWEEN-population).
  `corpus/data/` supplies **ZERO** images. The Stage-2 gate rests on 12-16's simulated arm (N ≥ 271).
- **`n_low` is a TIER-2 APPEND** from 12-15's truncation curve (`n_low_source = :tier2_from_minispike`).
  **THREE** Tier-2 sentinels, not two.
- **Guard 4 is PURELY DESCRIPTIVE** — no pass/fail fraction. Bars exist for Guards 1 and 3 only;
  Guards 2 and 4 carry a written-attribution duty.
- **`ridge_residual_shrinkage`** (12-06, 12-13) and **`shrinkage`** (12-15, 12-18) are **TWO
  ESTIMATORS, deliberately named differently.** 12-13:131-152 fixes a four-branch interpretation
  rule in advance. Never average them; never rename one to the other; a disagreement is a **FINDING**.
- `P12_GUARD_METRICS` is runner-local by choice.

## 8. Measurements already on the record

- **12-06 (ε-ridge):** `eps_ratio = 0.96307` — 5.8 sampling-sd below 1.0 at n_test = 12500, so ε is
  genuinely identified, the first nuisance in this project distinguishable from its prior. But
  `ridge_residual_shrinkage = 0.96870` against the frozen floor 0.90 → **`vacuous = true`**; 96.9% of
  the ε prior spread survives. Positive control 0.15742 (matches Phase 11's 0.157), so the harness is
  live and the result is informative. Consequence, under rules frozen before the number was seen:
  **S-4's radial chromatic confound is LIVE**, the 12-20 guards carry the SC3 claim, and named limit
  #4 is *sharpened*, not revised.
- **12-03 (advisory):** the research's CAR lag-1 table is NOT reproduced by the estimator the plan
  defines (0.1561/0.4117/0.5064/0.7222/0.9447 vs 0.136/0.351/0.438/0.692/0.947); the GP row matches
  exactly. Cause: CAR is non-stationary and the research did not record its pooling. Nothing gated
  depends on it. **Do not quote the two as the same measurement.**
- **12-07 (advisory):** the CAR marginal-sd spread runs OPPOSITE to what `12-RESEARCH.md` Pattern 1
  and `p12_lattice.jl:41-48` record — corner (1,1) = 0.9922 vs interior (4,4) = 0.6575 at α = 0.95,
  as it must be, since `Σ = (D − αW)⁻¹` gives low-degree cells the largest variance.

## 9. Compute discipline

- Any pool or cache under `spike/data/cache/` is **gitignored and unbacked**. If you generate one,
  ensure it regenerates byte-identically from the fixed counter-based (Philox) seed, and say so.
  **NEVER** record a pool in a SUMMARY as "available for reuse at no further compute cost."
- **The measurement is the deliverable.** Report what the numbers say, including when they say a
  thing is unidentifiable — Phase 11 closed negative-but-useful on exactly that and it was correct.
  Do NOT tune a threshold, widen a grid, or re-run for a nicer number. Do NOT spend the phase's
  iteration allowance. If a pre-registered bar is missed, **report the miss**.
- Do not shorten a pre-registered sample size to save time.
