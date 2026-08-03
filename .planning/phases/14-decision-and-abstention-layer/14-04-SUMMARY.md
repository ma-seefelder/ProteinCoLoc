---
phase: 14
plan: 04
subsystem: decision-primitives
status: COMPLETE
tags: [posterior, bayesian-fdr, class-order, d-04, d-07, d-08, wave-2]
requires:
  - "14-01 — spike/p14/consts.jl frozen at a6c8258 (P14_CLASS_KEYS, P14_FIX_SEED, P14_FIXTURE_COUNTER)"
  - "14-03 — spike/test/test_p14_decoupling.jl, the lane guard the two new modules pass"
  - "spike/p14/provenance.jl (473aac3) — p14_load_tau(), the only D-07 route to the three-way cut"
  - "spike/p13/{net,result,labels}.jl — read-only: load_three_way, P13_LOG_BF_KEYS, prior_class_draws, class_masses"
provides:
  - "p14_class_posterior — the three-class posterior, composed BY NAME from the MEASURED prior"
  - "p14_null_posterior / p14_null_split — the composite null {random ∪ exclusion}, with its decomposition"
  - "p14_confidence — the one confidence ordering the risk-coverage curve and the conformal hedge share"
  - "p14_rederive_class_masses — pi_class_masses re-derived, not trusted (Assumptions Log A7)"
  - "p14_bayes_fdr — the running-MEAN prefix rule, prefix property asserted"
  - "p14_fdr_over_decided — the batch wrapper that cannot report an FDR without its decided fraction"
affects:
  - "14-05+ (conformal, fuse, decide): every downstream layer sorts on p14_null_posterior and orders on p14_confidence"
  - "SC1-b/SC1-c integration runners consume p14_fdr_over_decided directly"
tech-stack:
  added: []          # zero packages; sortperm + cumsum + findlast are Base
  patterns:
    - "no positional index anywhere near a class — including a three-term `max` instead of `maximum`"
    - "the fixture is built by INVERTING the arithmetic, never by recording what the unit said"
    - "the equal-mass fixture is kept but LABELLED as a sanity check, not the class-order proof"
    - "the rule it must NOT be is written once, in the TEST, as the reference the fixture discriminates against"
    - "required keyword as a structural honesty device (n_total, following p13_calibration_meta's auc)"
    - "degenerate branch before the arithmetic (the roc_auc house shape)"
key-files:
  created:
    - spike/p14/posterior.jl              # 309 lines
    - spike/p14/fdr.jl                    # 240 lines
    - spike/test/test_p14_posterior.jl    # 257 lines, 48 assertions
    - spike/test/test_p14_fdr.jl          # 237 lines, 55 assertions
  modified: []
decisions: [D-04, D-07, D-08]
validates: ["SC1-a", "Class order", "Prior re-derivation"]
metrics:
  duration: ~55 min
  completed: 2026-08-03
  tasks_completed: 2
  tasks_total: 2
  per_file_pass_count: 103   # 48 posterior + 55 fdr
---

# Phase 14 Plan 04: The Three-Class Posterior and the Bayesian-FDR Prefix Rule — Summary

**The two numeric primitives the whole decision layer rests on now exist, and both of them fail
when subverted — proved by running the subversion, not by arguing about it.** A deliberately
permuted class assignment fails the posterior gate in 9 places; a running maximum substituted for
the running mean fails the FDR gate in 11. Neither substitution is hypothetical: the first is the
exact shape of the Phase-12 permutation artifact, and the second is the more conservative rule
that would never produce an obviously wrong answer.

## Commits, in order

| # | Task | Gate | Commit | Files |
|---|---|---|---|---|
| 1 | Three-class posterior | RED | `3ddbf8c` | `spike/test/test_p14_posterior.jl` (+257) |
| 1 | Three-class posterior | GREEN | `4f2ccc1` | `spike/p14/posterior.jl` (+309) |
| 2 | Bayesian-FDR prefix rule | RED | `56e2c3d` | `spike/test/test_p14_fdr.jl` (+237) |
| 2 | Bayesian-FDR prefix rule | GREEN | `6a45b7a` | `spike/p14/fdr.jl` (+240) |

No REFACTOR commit: neither implementation needed one after its GREEN gate.

## Per-file pass counts (the aggregate suite cannot be used, 14-PATTERNS §4.8)

| File | Assertions | Testset time | Wall time | Exit |
|---|---|---|---|---|
| `test_p14_posterior.jl` | **48** | 2.2 s | ~18 s | 0 |
| `test_p14_fdr.jl` | **55** | 1.6 s | ~12.6 s | 0 |
| `test_p14_consts.jl` (regression) | 273 | 0.2 s | — | 0 |
| `test_p14_provenance.jl` (regression) | **63** (was 59) | 2.2 s | — | 0 |
| `test_p14_decoupling.jl` (regression) | 42 | 1.3 s | — | 0 |

`test_p14_provenance.jl` gained four assertions without being edited: it greps the Phase-14 lane by
glob, so the two new modules armed checks that were previously waiting for them. That is the
grow-with-the-phase guard working as designed.

## The falsification runs, and what they showed

Both are recorded because the plan required them and because the *pattern of failure* is itself
evidence.

**Posterior — permuted class assignment** (scratch edit to line 203, verified present on disk by
`grep` before the run, and verified reverted by `grep` after):

```
p = (coloc = w_random / s, random = w_exclusion / s, exclusion = w_coloc / s)
→ 39 passed, 9 failed, exit 1
```

The instructive part is *which* testset stayed green: **the flat 1/3-each fixture passed all five
of its assertions under the permutation.** That is the file's own thesis made visible — an
equal-mass fixture cannot catch a permutation, which is why the unequal 0.7 / 0.2 / 0.1 target
under an unequal 0.2 / 0.5 / 0.3 prior is the actual proof and the flat fixture is labelled a
sanity check in both the code and the test output.

**FDR — running maximum substituted for the running mean** (scratch edit to line 156, mutation
confirmed on disk by the `cummax|accumulate(max` grep flipping 0 → 1, and confirmed reverted by it
flipping back to 0):

```
run = accumulate(max, sv)
→ 44 passed, 11 failed, exit 1
```

## What was built

### `spike/p14/posterior.jl`

`p14_class_posterior(logbf, prior)` composes `P(class | Z)` from the two Phase-13 log Bayes factors
and a class prior, by log-sum-exp over `log(prior_k) + logbf_k`. **There is no positional index
anywhere in the file** — not one `[1]`, `[2]` or `[3]`, and the log-sum-exp maximum is a
three-term `max` over named locals rather than `maximum` over a collection, precisely so no
reordering can change a number. Five class orderings coexist upstream and every headline metric is
invariant to a *consistent* relabelling, so a positional read would produce numbers that look
right.

The D-08 structural zero is checked with `===` and not `==`, because `-0.0 == 0.0` is true while
the two are different bit patterns and the claim is that the reference class was scored against
itself.

`p14_null_posterior` is the composite null `H₀ = {random ∪ exclusion}`, documented with all four
reasons that pooling is right rather than convenient — additivity over a partition needs no
approximation; a null of `{random}` alone would leave exclusion items *uncontrolled*, which is a
hole and not a tighter guarantee; it keeps the FDR unit identical to the decision unit (D-04); and
calling a segregated item a colocalization discovery is the worst available error.
`p14_null_split` carries the decomposition beside the total, because a batch whose false
discoveries are mostly exclusion mass is scientifically different from one whose are mostly random
mass.

`p14_confidence` is the three-term max by name — the κ the risk-coverage curve sweeps and the same
quantity the conformal hedge will threshold, recorded once so the two layers describe **one**
ordering rather than two that happen to agree on the fixtures.

### `spike/p14/fdr.jl`

`p14_bayes_fdr(v, alpha)` is `sortperm` + `cumsum` + `findlast`. **Zero packages added.** The
degenerate branch answers an empty batch *before* the arithmetic, `@assert issorted(run)` states
the prefix property rather than assuming it, `accepted` holds original indices, and `t_star` plus
its implied `cost_ratio` are reported so the decision-theoretic reading is available to a referee
without the implementation committing to a 3×3 loss matrix this project has no basis to elicit.
The docstring states the controlled quantity **with** its qualifier: `E[FDP | data] ≤ α`,
conditional on the model, **not** a frequentist long-run rate.

`p14_fdr_over_decided(v_decided, alpha; n_total)` makes `n_total` a **required** keyword, copying
the structural trick from `p13_calibration_meta` where `auc` is required so a calibration verdict
cannot ship without its discrimination number. `fdr_scope = :decided_subset_only` carries the scope
in machine-readable form. Assumption A2 (abstain-first, then sort) is written out as the two-line
σ(data)-measurability derivation it is, **and the reverse order is recorded as wrong** so it is not
re-derived later.

## The prior, re-derived rather than trusted (Assumptions Log A7)

`p14_rederive_class_masses` draws prior pairs through the existing `prior_class_draws` path,
labels them with the frozen `three_way_label` cut at the **inherited** τ (obtained only through
`p14_load_tau()`, never a literal), and reduces with the existing `class_masses`. Nothing is
re-implemented. It rides the fixture stream `p14_fix_rng(P14_FIXTURE_COUNTER)`, so a re-derivation
can never pre-observe a reported stream.

**Realized masses at n = 20,000, beside `h.meta.pi_class_masses`:**

| class | `h.meta.pi_class_masses` | re-derived | \|diff\| | 3 binomial SE | within band |
|---|---|---|---|---|---|
| coloc | 0.21655 | **0.2163** | 0.00025 | 0.008738 | ✅ |
| random | 0.47073 | **0.4713** | 0.00057 | 0.010588 | ✅ |
| exclusion | 0.31272 | **0.3124** | 0.00032 | 0.009834 | ✅ |

Every class agrees to well inside a *tenth* of its 3-SE band. A7 is corroborated: the measured
prior the FDR guarantee rests on is reproducible from the simulator prior and the frozen label
rule, and the numbers are printed by the test on every run rather than only asserted.

## Verification (observed, not inferred)

| Check | Command | Result |
|---|---|---|
| SC1-a + the FDR unit | `julia --project=spike spike/test/test_p14_fdr.jl` | **exit 0**, 55/55 |
| Class order + prior re-derivation | `julia --project=spike spike/test/test_p14_posterior.jl` | **exit 0**, 48/48 |
| Lane guard | `julia --project=spike spike/test/test_p14_decoupling.jl` | **exit 0**, 42/42 |
| Pre-registration | `julia --project=spike spike/test/test_p14_consts.jl` | **exit 0**, 273/273 |
| D-07 provenance | `julia --project=spike spike/test/test_p14_provenance.jl` | **exit 0**, 63/63 |
| Frozen surfaces | `git diff --exit-code HEAD -- src spike/Project.toml spike/Manifest.toml corpus spike/p14/consts.jl spike/p14/provenance.jl` | **exit 0** |
| No stray files under the frozen roots | `git status --porcelain -- src spike/Project.toml spike/Manifest.toml corpus artifacts` | empty |
| No file deleted by any of the four commits | `git diff --diff-filter=D --name-only HEAD~4 HEAD` | empty |
| No positional class read | `grep -v '^\s*#' spike/p14/posterior.jl \| grep -E '\[1\]\|\[2\]\|\[3\]'` | **0 lines** (0 even before the class filter) |
| The running max is absent from the unit | `grep -c 'cummax\|accumulate(max' spike/p14/fdr.jl` | **0** |
| `fdr_scope` is machine-readable and single-valued | `grep -n 'fdr_scope\s*=' spike/p14/fdr.jl` | 3 hits, the only assigned value is `:decided_subset_only` |

## Deviations from Plan

### Auto-fixed / judgement additions

**1. [Rule 2 — missing validation] `p14_bayes_fdr` validates that `v` holds finite probabilities**

- **Found during:** Task 2, writing the degenerate-input branch.
- **Issue:** the plan specified validation of `alpha` only. A `NaN` in `v` sorts to the end of the
  ascending order and turns every running mean after it into `NaN`, so the sequence the prefix
  property is asserted on is silently no longer ordered — the `@assert issorted(run)` would then
  fire with a message about *sorting* rather than about the caller's input, sending a reader to the
  wrong file. A value outside [0,1] is not a posterior probability at all.
- **Fix:** `_p14_check_null_posteriors` runs before the arithmetic and throws an `ArgumentError`
  naming the function, the offending index and the offending value. Three assertions cover it.
- **Files:** `spike/p14/fdr.jl`, `spike/test/test_p14_fdr.jl`. **Commit:** `6a45b7a` / `56e2c3d`.

**2. [Rule 2 — missing validation] `p14_class_posterior` rejects a non-positive prior mass**

- **Found during:** Task 1.
- **Issue:** `log(0.0)` is `-Inf`, which silently zeroes a whole class instead of erroring. Zeroing
  a class is a *decision*, not an arithmetic accident, and it must not be reachable by accident.
- **Fix:** explicit positivity check with a message naming the function and all three observed
  masses. **Commit:** `4f2ccc1`.

**3. [Rule 2 — missing validation] `p14_fdr_over_decided` rejects `n_total ≤ 0` and
`n_decided > n_total`**

- **Found during:** Task 2. A decided subset larger than its batch is a bookkeeping error, not a
  decided fraction above 1, and the whole point of the required keyword is that the fraction means
  something. **Commit:** `6a45b7a`.

**4. [Rule 3 — blocking] `posterior.jl` includes `provenance.jl` in addition to the four includes
the plan named**

- **Found during:** Task 1. `p14_rederive_class_masses` obtains τ through `p14_load_tau()`, which
  lives in `provenance.jl`; the plan named `consts.jl`, `p13/net.jl`, `p13/result.jl` and
  `p13/labels.jl` only. The added include is guarded on `:P14_PROVENANCE_LOADED` and is a superset
  of two of the four (provenance already pulls `consts.jl` and `p13/net.jl`), so nothing is
  loosened — and routing τ through the D-07 loader rather than a keyword is *stricter* than the
  alternative.

**5. [Design choice within Claude's discretion] `p14_rederive_class_masses` takes `n` as a
REQUIRED keyword and exposes no `tau` keyword**

- The plan's signature sketch was `(; n, rng)`. No sample size for this activity is frozen in the
  pre-registration, and inventing a default would put a chosen-looking number into a phase whose
  discipline is that numbers are measured or frozen first. The agreement band is
  `3·√(p(1−p)/n)`, so a caller who does not state `n` has not stated what the comparison means.
  Likewise τ is fetched inside the function body rather than offered as an overridable keyword, so
  there is no route to the cut that bypasses D-07.

### Assumption Drift (advisory)

**1. "Both test files complete in under 10 seconds each" — true of the test work, not of wall
clock, and the gap is package load.**

- **Planned:** `test_p14_fdr.jl` runs in under 10 s.
- **Actual:** testset time **1.6 s**; wall time **≈ 12.6 s**. `test_p14_posterior.jl` is 2.2 s of
  testset inside ≈ 18 s of wall.
- **Why:** both files reach the Phase-13 surfaces through guarded includes, which load Flux,
  NeuralEstimators and Images. That fixed cost is ≈ 16 s in this environment and is paid before a
  single assertion runs. It is not new to this plan: the already-committed
  `test_p14_provenance.jl` shows the same shape (2.2 s of testset inside 18.5 s of wall), because
  it includes the same stack.
- **Not repairable without breaking the plan:** the plan's own `key_links` require
  `fdr.jl → posterior.jl` and `posterior.jl → p13/net.jl`, and D-07 requires τ to come from the
  loader. Dropping the includes to save load time would sever exactly the links the plan asks for.
- **Non-blocking, recorded so the number is not later quoted as a regression.** The substantive
  claim — *no simulation ran, the primitives are pure arithmetic, feedback is fast* — holds:
  103 assertions in 3.8 s of combined testset time.

**2. The plan's `p14_class_posterior` sketch used `w[1]`, `w[2]`, `w[3]`; the implementation uses
named locals.**

- **Planned:** 14-RESEARCH's code example builds the result as
  `(coloc = w[1]/s, random = w[2]/s, exclusion = w[3]/s)`.
- **Actual:** three named locals and no indexing at all.
- **Why:** the plan's own acceptance criterion forbids a positional index on any line that also
  names a class — and the RESEARCH sketch's return line names all three *and* indexes all three.
  The named-local form satisfies the criterion, and it is the form the file's whole argument
  demands. The arithmetic is identical.

## Known Stubs

None. Both modules are complete implementations of their specified behaviour; nothing returns a
placeholder and nothing is wired to an empty data source.

## Threat Flags

None. No file created here opens a network endpoint, reads an image, touches an auth path, or
changes a schema at a trust boundary. Both modules are pure arithmetic over already-loaded numbers;
the only file-system reach is the read-only `p14_load_tau()` route that `14-02` already shipped and
`test_p14_provenance.jl` already guards.

## Threat register, discharged

| ID | Threat | How it is now mitigated |
|---|---|---|
| T-14-01 | Class order permuted in `p14_class_posterior` | Every read and write by name; zero positional indices in the file; unequal 0.7/0.2/0.1 fixture asserted field by field; shuffled-declaration-order fixtures; **falsification run shows 9 failures and, tellingly, the flat fixture still passing** |
| T-14-06 | FDR reported without the decided fraction | `n_total` is a required keyword (`UndefKeywordError` asserted); `fdr_scope = :decided_subset_only` is machine-readable |
| T-14-18 | Running max substituted for running mean | Discriminating fixture asserts `k* = 3` **and** asserts the max rule would say 2; `grep` for the max idiom returns 0; falsification run fails 11 |
| T-14-19 | Prior silently chosen rather than measured | `prior` is an explicit argument, documented as `h.meta.pi_class_masses`; `p14_rederive_class_masses` reproduces it inside 3 binomial SEs and prints both tables |
| T-14-20 | Log-sum-exp omitted, ±22-nat log BFs overflow | Max subtracted before `exp`; `logbf.coloc = 800.0` and the mirrored `−800.0` underflow case both asserted finite and summing to 1 |
| T-14-SC | Package installs | Zero packages; `spike/Project.toml` and `spike/Manifest.toml` byte-unchanged, asserted by the running lane guard |

## What this does NOT establish

- **No FDR number has been produced.** This plan built the *rule*; SC1-b/SC1-c measure whether the
  realized FDP tracks the level, and that runs on fresh simulator draws in a later plan.
- **The prior sensitivity is unmeasured.** `pi_class_masses` is the *simulator prior's* class mix,
  not any real batch's prevalence. A batch that is 90 % coloc makes every null posterior too large
  and the rule too conservative; a batch that is 1 % coloc makes them too small and the claim
  false. The sweep over `P14_PI_COLOC_GRID` is what bounds this, and it is reported, not gated.
- **Assumption A2 remains an assumption.** It is now written out as a derivation in the code, which
  is what the validation contract asks for, but the manual-verification row still stands: the
  report must state it explicitly so a referee can check it rather than take it on trust.
- **D-04 is unchanged.** FDR is controlled per image pair. Per-region/per-tile FDR is out of scope
  and is not merely unimplemented — Phase 12 returned NO, and `LocalColocMap` carries no
  uncertainty field to control against.

## Self-Check: PASSED

| Claim | Check | Result |
|---|---|---|
| `spike/p14/posterior.jl` exists | `[ -f ]` | FOUND |
| `spike/p14/fdr.jl` exists | `[ -f ]` | FOUND |
| `spike/test/test_p14_posterior.jl` exists | `[ -f ]` | FOUND |
| `spike/test/test_p14_fdr.jl` exists | `[ -f ]` | FOUND |
| `3ddbf8c` in history | `git log` | FOUND |
| `4f2ccc1` in history | `git log` | FOUND |
| `56e2c3d` in history | `git log` | FOUND |
| `6a45b7a` in history | `git log` | FOUND |
