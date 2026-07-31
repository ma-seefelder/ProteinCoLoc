---
phase: 12-spatial-colocalization-map
plan: 14
status: complete
subsystem: spike-npe
tags: [training-surface, leak-freedom, mask-ordering, stage1-gate-guard, theta-truncation, D-10]
requires:
  - "spike/npe/p12_architecture.jl (build_p12_estimator, reshape_summary, mask_regions, augment_mask!, P12_THETA_ROWS, p12_theta_index)"
  - "spike/data/p12_generate.jl (load_p12_pool, p12_pool_dir, p12_pool_complete)"
  - "spike/validation/p12_consts.jl (p12_stage1_verdict, P12_MASK_K_SET, P12_ZSCORE_ARM, P12_MINISPIKE_N/COUNTER)"
provides:
  - "spike/npe/train_p12_npe.jl — THE single Phase-12 training surface (train/save/load + p12_truncate_theta)"
  - "spike/test/test_p12_train.jl — 64 tests over leak-freedom, mask ordering, the gate guard and truncation"
affects:
  - "12-15: trains three arms through train_p12_npe; K_dev defaults keep the mini-spike full-width"
  - "12-17: trains the production head via K_dev = P12_K_PROD — UNBUILDABLE without p12_truncate_theta"
  - "12-16/12-18/12-19/12-20: all read bundles through load_p12_npe, which REFUSES a D/theta_rows mismatch"
tech-stack:
  added: []
  patterns:
    - "arm-as-keyword so CAR / GP / ablation are one code path (what makes SC3's claim attributable)"
    - "gate guard with a smoke-path exemption and a verdict_dir seam for testing all three branches"
---

# 12-14 — the single Phase-12 training surface, and the two failure modes it makes impossible

`train_p12_npe.jl` exists, its transforms cannot leak, a held-out region is encoded byte-identically
to a naturally unusable patch, and it will not spend compute before the gate that exists to stop it.
**64 tests pass, 0 fail, 0 error, in 91.3 s** (budget 120 s).

**Task 3 was a NO-OP, as designed.** The Stage-1 verdict is `:proceed`, so the D-13 ablation branch
did not fire, `p12_ablation_report.jld2` does not exist, and no compute was spent there. That branch
trains only on a DESCOPE, where 12-11 routes the fallback through this plan.

## Declared deviations — the plan says one thing and this implementation does another

### DEVIATION 1 — `augment_mask!` is CALLED, never redefined (a namespace collision the plan would have created)

`12-14-PLAN.md:154` specifies `augment_mask!(Zraw, rng; G, kset)` as one of this file's functions.
**That function already exists**: `spike/npe/p12_architecture.jl:190`, written by 12-04, with the
keyword spelled `k_set` — and this trainer `include`s that file into the same scope. Defining it
again would have **silently overwritten the pre-registered 12-04 augmentation**, with the winner
decided by include order, and `test_p12_architecture.jl` would then have been exercising the
trainer's copy rather than the one it was written for.

This is the **second namespace collision in this phase**, after 12-12's `_res` helper overwrote its
Phase-13 namesake. It was not escalated as a blocker because the plan's *intent* is unambiguous — it
wants the realized masking rate recorded (F6: "the intended distribution is not evidence of what a
run contained") — no threshold or pre-registered value moved, and the repair is forced rather than
discretionary.

**What replaces it.** 12-04's `augment_mask!` returns `Zraw`, not the drawn `k`, so the trainer
MEASURES what the net saw from the data instead: `_p12train_masked_counts` is called before and
after the augmentation. That is strictly better evidence than the drawn `k` — it counts what
actually reached the net — and it keeps the augmentation single-sourced.

**A trap this fix created, and how it was closed.** The plan's testset 4 asks that "the returned
realized rate matches a recount of the zeroed mask entries". Under this change that would be an
**identity — a test that cannot fail**, which is precisely the defect class this phase has already
found six times. It is kept falsifiable by exercising 12-04's `augment_mask!` on an **all-ones-mask
fixture**, where no natural zeros exist and the recount therefore IS the drawn `k`, so
"every realized k ∈ P12_MASK_K_SET" is a real constraint. A source-level assertion additionally
requires that `function augment_mask!` never appear in the trainer, closing the reintroduction path.

**THE THREE MASK RATES ARE NAMED FOR WHAT THEY COUNT, AND NONE IS THE DRAWN `k`.**

| Key | What it is |
|---|---|
| `natural_mask_rate` | zeroed mask rows per region **BEFORE** augmentation — the forward model's own ≥15-surviving-pixel floor |
| `realized_mask_rate` | zeroed mask rows per region **AFTER** augmentation. A **TOTAL, not a delta**. (The key name 12-14 mandates.) |
| `augmented_only_mask_rate` | the difference, i.e. **NEWLY** zeroed rows — which **UNDER-COUNTS the drawn `k`** by any overlap with an already-masked region |

`realized_mask_rate * G²` equals the mean drawn `k` **only when the pre-augmentation mask has no
zeros**. On the F5 mixture natural masking is rare, so the two are usually close — which is exactly
why the distinction is written into the bundle as `mask_rate_naming_note` rather than left to be
inferred from a number that looks about right. Comparing any of them against `P12_MASK_K_SET`
without that condition holding would be the same unit error as conflating a ridge residual ratio
with a posterior width.

### DEVIATION 2 — the pool-generation command is named by FILE, not by symbol

The plan asks the refusal message to "name the generation command", *and* requires that after
stripping comments the file contain **zero** occurrences of `generate_p12_pool`. **Those cannot both
be met literally** — the first put the symbol in a string literal and the second's text scan found
it there. Resolved in favour of the checkable one: the message names
`spike/data/p12_generate.jl` and the plans that own the spend (12-15, 12-17), so a reader can still
act on it, while "the trainer never generates" stays verifiable by a text scan rather than by
review. The reasoning is recorded in the message itself.

### DEVIATION 3 (recorded, from the plan's own `<context>`) — masking on RAW rows

`12-09-PLAN.md:69` says the MCAR masking is applied "to the standardized input". **This file does
not follow it**, as 12-14 instructs. Masking happens on the RAW rows so a held-out region is encoded
exactly as a patch `correlation()` could not score: raw `0.0`, mask `0`, which after per-row
z-scoring is `(0 − μ_r)/σ_r` and **not** zero. Masking after standardization would place the row
*mean* at the held-out region — "this region is perfectly average" where read time says "this region
is absent". This is a deviation from a sibling plan's PROSE, not from any decision.

## A library-API surprise, the second in this phase

**`sampleposterior`'s `N` is a KEYWORD in NeuralEstimators v0.2.1, not positional.** The natural
`sampleposterior(est, Z, 32)` resolves to a *different* three-argument method and throws a
`MethodError`. Verified against the installed signature at
`NeuralEstimators/src/Estimators/PosteriorEstimator.jl:138`, not from memory. This is the second
v0.2.1 signature/default surprise here, after `use_gpu` defaulting to `true` (which 12-04 hit), and
it is exactly the hazard CLAUDE.md flags: *"API surface is young (v0.x); verify signatures against
`dev` docs, not memory."* The correct call is `sampleposterior(est, Z; N = 32, use_gpu = false)`.

## `p12_truncate_theta` — the operator whose absence survived three plan-checker passes

The pool's θ is always 72 rows; a production head at `K < 63` is that layout with the deviation block
shortened. Without the operator **12-17 is unbuildable**: the trainer would fit `D = 72`, write a
72-entry `theta_rows`, and `load_p12_npe` — which *refuses*, never repairs, a bundle whose
`theta_rows` disagrees with its `D` — would reject the production bundle for every downstream
consumer. `K == P12_K_DEV` returns the input unchanged, so the default path is a **provable no-op**,
and a `K` wider than the pool layout is refused rather than silently padded.

## Verification

| Check | Result |
|---|---|
| `test_p12_train.jl` per-file | **64 pass / 0 fail / 0 error**, 91.3 s (budget 120 s) |
| Plan's Task-1 automated verify | PASS — prints exactly `72 17 per_row train-ok` |
| Plan's Task-3 verify (PROCEED route) | PASS — `no-op by design`; no `p12_ablation_report.jld2` exists |
| **Testset 3 fails when the trainer standardizes before masking** | **CONFIRMED RED** (property true→false under mutation) |
| **Testset 1 fails when `fit_p12_transforms` sees the full pool** | **CONFIRMED RED** (property true→false under mutation) |
| `grep -c P12_PENDING_SCAFFOLD spike/test/test_p12_train.jl` | 0 |
| Exactly one standardizer fit site (2 `fit` calls, both inside `fit_p12_transforms`) | PASS |
| Zero `generate_p12_pool`, zero `use_gpu = true`, zero `Pkg.add` in code | PASS |
| `function augment_mask!` / `function p12_stage1_verdict` absent from the trainer | PASS |
| mask→standardize→reshape order visible in one expression | PASS (`reshape_summary(standardize_p12(Zraw_train_masked, zt), …)`) |
| `git diff` on `src`, `Project.toml`, `Manifest.toml`, `corpus`, `p12_consts.jl` | PASS — 0 changed |
| `git status --porcelain -- spike/data/cache spike/npe` | clean apart from this plan's own new trainer |

**The falsification harness never wrote to the tracked trainer.** A phase-13 executor and the user
also commit to this branch, so mutating a tracked file in place — even briefly — risks another agent
committing the mutation. Each mutant was written to a *sibling scratch file* in the same directory
(so `@__DIR__`-relative includes still resolve), loaded in a fresh process where the real trainer was
never included, probed, and deleted. Byte-identity of the real file was asserted afterwards, and no
scratch file was left behind.

## The gate guard, and why all three branches are testable

`_assert_spend_allowed` exempts the smoke path (`epochs <= 1 && n <= 512`) and says so in its own
message, so a test run never depends on the repository's live verdict state. The `verdict_dir`
keyword threads to Tier-1's `p12_stage1_verdict(; dir = …)`, which is what lets testset 7 exercise
**absent / DESCOPE / PROCEED** against a temp dir — necessary because 12-11 (wave 6) always writes a
verdict before this plan (wave 7) runs, so a bare call can never reach the `:absent` branch.

The DESCOPE asymmetry is asserted by *message*, not merely by "it threw": `arm = :none` is confirmed
to get **past** the guard and fail later on the absent `:none` pool (`no COMPLETE pool`), which is
what distinguishes "permitted by the guard" from "blocked by the guard". `arm = :car` on DESCOPE
fails with the guard's own refusal.

## Commits

| Commit | Contents |
|---|---|
| (this) | `spike/npe/train_p12_npe.jl`, `spike/test/test_p12_train.jl`, this summary |
