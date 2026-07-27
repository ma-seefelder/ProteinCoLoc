---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 06
subsystem: amortized-inference
tags: [neuralestimators, flux, three-way-bayes-factor, masked-bce, custom-estimator-subtype, jld2]

# Dependency graph
requires:
  - phase: 13-three-hypothesis-amortized-bayes-factor (plan 13-01)
    provides: "spike/p13/consts.jl -- the frozen Tier-1 pre-registration (P13_SUMMARY_WIDTH, P13_NUM_SUMMARIES, the copied recipe, P13_SMOKE_N/EPOCHS, P13_USE_GPU, P13_DEV_SEED/P13_FIX_SEED, P13_INPUT_WIDTH_RULE, P13_LAMBDA_PLACEMENT)"
  - phase: 13-three-hypothesis-amortized-bayes-factor (plan 13-03)
    provides: "spike/p13/labels.jl -- ThreeWayClass, target_matrix (the 4-row [y_C; w_C; y_E; w_E] legend), head_log_odds"
provides:
  - "ThreeWayEvidenceNet <: NeuralEstimator -- one shared trunk, one Dense(num_summaries, 2) head, 2-by-n logits in ONE forward pass"
  - "build_three_way_net / three_way_arch -- the D-10 trunk verbatim, with architecture metadata read back off the layers"
  - "masked_two_head_bce -- the D-11 per-head class restriction as participation weights"
  - "train_three_way -- route R1: NeuralEstimators' own train loop with an HONOURED custom loss, copied recipe, CPU-only hard gate, global-RNG reseed"
  - "three_way_input_dim(G, n_cond) / p13_ratio_input_dim(G) -- the DERIVED input width (never a literal)"
  - "p13_pair_encode / encode_conditioned_pair -- the pair encoding with the iseven tripwire and the conditioning row appended LAST"
  - "three_way_log_bf / three_way_probs -- the D-08 read surface (random == 0.0, per-head correction) and the UNCORRECTED probabilities D-13 scores"
  - "save_three_way / load_three_way -- atomic, narrow-surface persistence (Flux.state + arch metadata, not the object)"
  - "spike/test/test_p13_net.jl -- 91 assertions incl. the A5 two-epoch train smoke and the D-11 inactive-logit bit-identity"
affects: [13-08 test wiring, 13-09 preconditions binding, 13-11 the real training run, 13-12 the D-12/D-13 gate, 13-13 alpha series read, 13-16 real-image read, 14-decision-and-abstention]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "custom NeuralEstimator subtype (route R1) to obtain a multi-logit head with an honoured custom loss"
    - "masked multi-head BCE with per-sample participation weights and an active-cell denominator"
    - "input width derived from a grid helper plus a conditioning length, asserted never-literal by a source grep"
    - "narrow-surface model persistence: rebuild topology from metadata, load state into it"

key-files:
  created:
    - spike/p13/net.jl
    - spike/test/test_p13_net.jl
  modified: []

key-decisions:
  - "Route R1 SUCCEEDED: NeuralEstimators 0.2.1's train accepts ThreeWayEvidenceNet and honours the passed masked loss; the pre-declared route-R2 hand-rolled loop was NOT used, so D-10's attribution argument holds unweakened"
  - "p13_ratio_input_dim and p13_pair_encode are attributed local copies (src/amortized/summary.jl:127, src/amortized/bf.jl:91-100) because neither src/ file is loadable standalone and the spike twin drags the whole Phase-5 harness; src/ stays byte-unchanged"
  - "train_three_way exposes a `loss` keyword (defaulting to masked_two_head_bce) so the two-different-losses check can prove the loss is not silently discarded"
  - "The 'never write the width as a literal' and 'do not copy the binary prior-odds term' warnings live in # comments, not docstrings, because the machine-checked source greps strip only comment lines"
  - "The A5 smoke splits 160/40 rather than at the pre-registered val_frac: NeuralEstimators' DataLoader uses partial=false, so a validation set smaller than batchsize yields zero batches and a 0/0 validation risk"

patterns-established:
  - "Pattern: prove a pre-1.0 package's internal extension hooks with a real two-epoch training run before any expensive run, and name the fallback as a declared deviation in advance"
  - "Pattern: assert a head-boundary invariant with EXACT equality (perturbing an inactive logit must not move the loss), because an approximate comparison would accept a leak"
  - "Pattern: forbid an identifier from executable code with a comment-stripped source grep, and keep docstrings free of the forbidden token"

requirements-completed: [D-08, D-09, D-10, D-11, D-01, D-03]

# Metrics
duration: 41min
completed: 2026-07-27
---

# Phase 13 Plan 06: Three-Way Evidence Network Summary

**A `ThreeWayEvidenceNet <: NeuralEstimator` with the shipped 256/64 trunk verbatim and a single `Dense(64,2)` head, trained through NeuralEstimators' own loop with a masked two-term logit-BCE that provably restricts each head to its own class pair, reading out `log BF(coloc:random)` and `log BF(exclusion:random)` against an identically-zero random reference in one forward pass.**

## Performance

- **Duration:** ~41 min
- **Started:** 2026-07-27T00:00:00Z (approximate; worktree session)
- **Completed:** 2026-07-27
- **Tasks:** 2
- **Files modified:** 2 created, 0 modified

## Accomplishments

- **Route R1 is proven, not assumed.** The A5 smoke runs a real two-epoch training pass on a 200-sample toy: NeuralEstimators 0.2.1's `train` accepts the custom subtype (generic `_construct_train_state`, `_inputoutput`) and **honours the passed loss** (generic `_loss(estimator, loss) = loss`). The pre-declared route-R2 hand-rolled loop was **not** needed, so the optimizer, the Float64 CosAnneal schedule, the patience rule and the best-validation checkpoint are byte-for-byte the shipped binary net's — which is exactly what D-10's "any difference is attributable to the head alone" claim requires.
- **D-11 is machine-checked, not documented.** Perturbing an INACTIVE logit by 50 leaves `masked_two_head_bce` **bit-identical** (asserted with `==` and `===`, never an approximate comparison); perturbing an ACTIVE logit moves it. That single pair of assertions is what distinguishes this loss from one-vs-rest, which would pass every shape test while silently turning each logit into a mixture Bayes factor.
- **The evidence-scale correction cannot be confused with the binary one.** `three_way_log_bf` subtracts the per-head `head_log_odds` NamedTuple and nothing else; the binary read surface's own term is asserted **absent from executable source** by a comment-stripped grep, because a wrong copy is numerically invisible in the file it would have been copied from.
- **The input width is derived and the literal is forbidden.** `three_way_input_dim(G, n_cond)` = `5*G^2 + n_cond`, with `n_cond` supplied by the caller so plan 13-09 can bind the Phase-11 conditioning length later. The literal never appears in executable code (grep count 0), and the wrong lambda placement still fails loudly on the kept `iseven` guard.
- **91 assertions, ~25 s, CPU-only, tau-independent.** The suite passes while `P13_TAU` is still Tier-2 (the net depends on tau only through the label boundary), and asserts CUDA absent.

### The A5 smoke's measured numbers (the plan's required record)

| Quantity | Value |
|---|---|
| Route taken | **R1** (`<: NeuralEstimator` subtype through NeuralEstimators' `train`) — SUCCEEDED |
| Route R2 (hand-rolled loop) used? | **No** |
| Toy | `P13_SMOKE_N = 200` samples, `input_dim = 32`, class means at −2/0/+2, split 160/40 |
| Realized class counts | exclusion 64, random 70, coloc 66 |
| Epochs | `P13_SMOKE_EPOCHS = 2`, batchsize 32, `P13_LR = 2.5e-4`, AdamW(Float64) |
| Masked loss, freshly-initialized net | **0.58271** |
| Masked loss, returned (best-validation) net | **0.35380** |
| Validation risk trajectory | 0.56248 (initial) → 0.44 (epoch 1) → 0.353 (epoch 2) |

The validation risk improved at **both** epochs, so the returned net is genuinely the epoch-2 state and the `after < before` assertion is not marginal (`train` initializes its best-state to the untrained net, so an equal loss would have meant the loss never reached the optimizer).

## Task Commits

1. **Task 1: the two-head evidence net, the masked loss and the read surface** — `41d647f` (feat)
2. **Task 2: the net test incl. the A5 smoke and the D-11 inactive-logit invariance** — `1184db7` (test)

**Plan metadata:** committed with this SUMMARY (docs).

## Files Created/Modified

- `spike/p13/net.jl` (548 lines) — `ThreeWayEvidenceNet` (+ `Flux.@layer`), `build_three_way_net`, `three_way_arch`, `p13_ratio_input_dim`, `three_way_input_dim`, `masked_two_head_bce`, `train_three_way`, `three_way_log_bf`, `three_way_probs`, `p13_pair_encode`, `encode_conditioned_pair`, `save_three_way`, `load_three_way`, `p13_consts_sha`, and a guarded script entry point that prints a usage note and trains nothing.
- `spike/test/test_p13_net.jl` (424 lines) — 14 testsets / 91 assertions: SC1 one-pass shape, the walked D-10 topology, the derived width plus the never-literal tripwire, the D-11 bit-identity, the hand-computed active-cell denominator, the A5 train smoke, the two-different-losses honoured-loss check, lambda placement (both wrong forms), the D-08 read surface, the Pitfall-1 source greps, the atomic round-trip, input-boundary shape throws, the CPU-only gate, and the tau-still-Tier-2 clause.

## Decisions Made

1. **`p13_ratio_input_dim` and `p13_pair_encode` are attributed local copies, not reached through an include.** The plan permitted either the spike-lane twin or an attributed copy. The twin (`spike/validation/train_ratio.jl`) includes the whole Phase-5 harness (a frozen-model load) and carries a `RATIO_INPUT_DIM = 320` **literal** rather than the grid-general helper; the `src/` originals are not loadable standalone (`src/amortized/summary.jl` mentions `MultiChannelImage` in a method signature, so including it would drag the frozen image chain in for one arithmetic identity). Both copies carry an attribution comment naming the source file and line range and the reason, and the `iseven` guard is kept exactly as written. `src/` is byte-unchanged (asserted).
2. **`train_three_way` gained a `loss` keyword** (defaulting to `masked_two_head_bce`). The plan's signature list omitted it, but task 2's "custom loss is actually honoured" check requires training the same toy under two different losses — which is the only available proof, short of touching package internals, that the shipped type's silent-discard landmine has not been inherited.
3. **`three_way_log_bf` hard-throws on `use_gpu = true`** rather than merely defaulting it false. The keyword exists only so the signature is diffable against the shipped read surface; CPU-only is a hard constraint, so silently accepting `true` would be worse than throwing.
4. **The two-different-losses check trains from ONE deep-copied initialization.** Two freshly built nets differ by their random initialization alone, which would have made the comparison vacuous.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Two plan clauses conflicted; the machine-checked one won**

- **Found during:** Task 1 (writing `net.jl`)
- **Issue:** The plan's `<action>` asked for the docstring of `three_way_input_dim` to say, in capitals, "NEVER write 321 as a literal", and for `three_way_log_bf`'s docstring to name the binary `- log_prior_odds` term as the thing not to copy. Both instructions place the forbidden tokens on **non-comment lines**, which the plan's own acceptance criteria and `<verification>` block assert must contain **zero** occurrences of `321` and `log_prior_odds` (`grep -v '^\s*#' | grep -c ...` → `0`). Docstrings are not comment lines, so following the prose literally would have failed the gate — and would also have left a copy-pasteable forbidden literal in the rendered help text.
- **Fix:** Both warnings are preserved verbatim in intent and in capitals, but placed in `#` comment blocks immediately adjacent to the definitions (the derived-width block comment above `p13_ratio_input_dim`, and banner clause (d.1) above the read surface). The docstrings state the same rule without naming the token ("DERIVED, NEVER A LITERAL"; "the binary read surface's own correction term"). The test asserts both absences on the comment-stripped source.
- **Files modified:** `spike/p13/net.jl`
- **Verification:** `grep -v '^\s*#' spike/p13/net.jl | grep -c '321'` → `0`; same for `log_prior_odds`, `RatioEstimator(`, `softmax`, `kde(`/`quadgk(` → `0`. Testsets "input width is derived" and "the binary correction is not copied (Pitfall 1)" assert the same on `P13_NET_CODE`.
- **Committed in:** `41d647f` (Task 1 commit)

**2. [Rule 3 - Blocking] The A5 smoke split had to be 160/40, not the pre-registered `val_frac`**

- **Found during:** Task 2 (writing the A5 smoke)
- **Issue:** NeuralEstimators builds its loaders with `partial = false` (`train.jl:668-676`), so an incomplete final batch is **dropped**. At `P13_SMOKE_N = 200` with `P13_VAL_FRAC = 0.15` the validation set would hold 30 samples < `batchsize = 32`, giving **zero** validation batches and a `0/0 = NaN` validation risk — the early-stopping and best-checkpoint logic would then compare against `NaN` and the smoke would be meaningless rather than failing loudly.
- **Fix:** The toy splits 160/40 so the validation set holds at least one full batch, with an in-test assertion (`size(Z_va, 2) == NET_FIX_NVA >= 32`) and a comment naming the `partial = false` cause. This is a property of the toy's size only; **no pre-registered constant was changed** (`P13_VAL_FRAC` is untouched and is what plan 13-11's real run will use at n = 48 000, where 15% is 7 200 ≫ 128).
- **Files modified:** `spike/test/test_p13_net.jl`
- **Verification:** The A5 testset passes with a finite, monotonically improving validation risk (0.562 → 0.44 → 0.353) and `after < before`.
- **Committed in:** `1184db7` (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (both Rule 3 - blocking). No architectural change, no Rule 4 escalation, no pre-registered constant touched.
**Impact on plan:** Both were forced by facts on the ground (the plan's own machine-checked criteria in the first case, an installed-package behaviour in the second). Every substantive instruction — the topology, the recipe, the masked loss, the per-head correction, the derived width, the atomic persistence, the CPU-only gate — was implemented as written. No scope creep.

## Assumption Drift (advisory)

**1. `train`'s returned net is the BEST-VALIDATION state, not the final-epoch state**

- **Found during:** Task 2 (designing the A5 pass criterion)
- **Planned:** The plan asks that "the training loss on the toy at the end is strictly below the loss of the freshly-initialized net", implicitly reading the returned object as the end-of-training net.
- **Actual:** `train` initializes `min_val_risk` to the **initial** validation risk and `trainstate_best` to the **untrained** state (`train.jl:227-237`), returning `trainstate_best`. If no epoch improved validation risk, the returned net *is* the fresh net and the two losses would be exactly equal.
- **Why it matters to a reader:** the assertion is therefore stronger than it looks — it also proves at least one epoch improved validation risk, i.e. that the loss actually reached the optimizer. It also means a marginal toy could produce a false failure, which is why the toy's class means are separated by 2σ and why the realized validation trajectory (improving at both epochs) is recorded above rather than left implicit.

Advisory only; nothing was gated on it and no plan instruction was overridden.

## Issues Encountered

- **`RatioEstimator`'s hard-coded loss** — the premise of the whole plan — was re-verified against the installed source (`RatioEstimator.jl:124` vs the generic `train.jl:654`) before writing a line, and the generic-hook path was probed end-to-end in a scratch script before committing. No surprises.
- **Broadcasting a keyword-argument loss** (`logitbinarycrossentropy.(...; agg = identity)`) was verified to work under Julia 1.12.6 / Flux 0.16.10 rather than assumed.
- **`Flux.state` + JLD2 + `loadmodel!`** round-trip for a custom `@layer` struct was probed before being written into `save_three_way`/`load_three_way`; the rebuilt net reproduces the saved net's output exactly (`max |Δ| = 0.0` in the probe, asserted `<= 1e-6` in the suite).

## Verification (observed, not inferred)

Every command below was RUN from the worktree root and its real output observed:

| Check | Command | Observed |
|---|---|---|
| Suite (script) | `julia --project=spike spike/test/test_p13_net.jl` | **91 passed / 91 total**, 0 fail, 0 error, 25.3 s test time, 31 s wall |
| Suite (includable) | `julia --project=spike -e 'include("spike/test/test_p13_net.jl")'` | **91 passed / 91 total** |
| Shape | `... build_three_way_net(320)(randn(Float32,320,7))` | `shape ok` (size `(2, 7)`) |
| Derived width | `three_way_input_dim(8,1) == 321`, `(8,2) == 322` | `width ok` |
| Never-literal | `grep -v '^\s*#' spike/p13/net.jl \| grep -c '321'` | `0` |
| Rejected architectures | same grep for `RatioEstimator(` / `softmax` | `0` / `0` |
| Named limit 3 | same grep for `kde(\|quadgk(` | `0` |
| Pitfall 1 | same grep for `log_prior_odds` | `0` |
| One-pass literal | `grep -c 'm.head(m.trunk(' spike/p13/net.jl` | `1` |
| RNG discipline | `grep -c 'Random.seed!' spike/p13/net.jl` | `3` (one call inside `train_three_way`, two in prose) |
| Exact-equality gate | `grep -A6 'BIT-IDENTICAL' ... \| grep -c isapprox` | `0`; the window shows `==` and `===` |
| A5 named | `grep -c 'A5' spike/test/test_p13_net.jl` | `5` |
| Decoupling | `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | exit `0` |
| Frozen pre-registration | `git status --short` after both commits | clean; `spike/p13/consts.jl` and `spike/test/runtests.jl` untouched |

**Not verified here (out of this plan's scope):** the net has been trained on nothing but the 200-sample toy. Its discrimination (D-12), its calibration (D-13), the per-head correction's numerical exactness on the §F5 conjugate-Gaussian toy, the lambda-response of the log-BF pair, and the realized conditioning length all belong to later plans. Nothing in this plan measures or claims any of them.

## User Setup Required

None — no packages installed, no external service, no environment variable. `spike/Project.toml` and `spike/Manifest.toml` are byte-unchanged and CUDA remains absent.

## Next Phase Readiness

**Ready:**
- The architecture, loss, trainer, read surface and persistence are in place and gated. Plan 13-11's trainer can call `build_three_way_net(three_way_input_dim(G, n_cond))` and `train_three_way(...)` directly.
- The A5 route-R1 risk is **retired** before any expensive run; D-10's attribution argument is intact.
- `three_way_probs` gives D-13 the uncorrected probabilities it must score; `three_way_log_bf` gives D-08/D-12 the corrected pair.

**Blockers / concerns for downstream plans:**
- `spike/test/test_p13_net.jl` is **not yet wired into `spike/test/runtests.jl`** — that is plan 13-08's job (this plan did not touch `runtests.jl`).
- `n_cond` is a **caller-supplied parameter**, deliberately unbound here. Plan 13-09's `spike/p13/preconditions.jl` must bind it from the Phase-11 handle (`length(encode_lambda(lambda))`); the width helper is ready for a value other than 1.
- Training still blocks on the Phase-11 research NPE (D-02) for the frozen `zt`, and on the Tier-2 tau measurement (plan 13-10) for any labelled datum.
- `save_three_way` records `tau = nothing` while tau is Tier-2. After plan 13-10 appends the measured value, artifacts produced by 13-11 will carry it automatically — but any artifact saved before that point is, correctly, marked tau-less.

## Self-Check: PASSED

- `spike/p13/net.jl` — FOUND (30 167 bytes)
- `spike/test/test_p13_net.jl` — FOUND (23 289 bytes)
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-06-SUMMARY.md` — FOUND
- commit `41d647f` — FOUND in `git log`
- commit `1184db7` — FOUND in `git log`
- `git status --short` — clean (no untracked or uncommitted file left behind)
- `spike/p13/consts.jl`, `spike/test/runtests.jl`, `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `.planning/STATE.md`, `.planning/ROADMAP.md` — untouched

---
*Phase: 13-three-hypothesis-amortized-bayes-factor*
*Completed: 2026-07-27*
