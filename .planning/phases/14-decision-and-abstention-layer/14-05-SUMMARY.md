---
phase: 14
plan: 05
subsystem: hedge-and-abstention
status: COMPLETE
tags: [conformal, abstention, d-02, d-05, d-06, sc1-amended, sc2-amended, wave-2]
requires:
  - "14-01 — spike/p14/consts.jl frozen at a6c8258 (P14_CLASS_KEYS, P14_AMBIGUOUS_RATE_FLOOR, P14_FIXTURE_COUNTER)"
  - "14-03 — spike/test/test_p14_decoupling.jl, the lane guard both new modules pass (its SC2-c exemption is now armed)"
  - "14-04 — spike/p14/posterior.jl (4f2ccc1): p14_class_posterior / p14_null_posterior / p14_confidence / _p14_check_class_keys"
  - "src/results.jl — read-only: the OODVerdict type"
provides:
  - "p14_lac_score — the LAC/THR nonconformity score, on the SAME posterior the FDR rule sorts on"
  - "p14_conformal_quantile — the ceil(Int, (n + 1) * (1 - alpha))-th ORDER STATISTIC, never an interpolated quantile"
  - "p14_conformal_set — the three set-size regimes, with :empty kept as its own signal"
  - "p14_hedge_diagnostics — the realized set-size distribution with vacuous_hedge as a LABEL"
  - "p14_conformal_calibrate — the labelled-pair fit, zero dependencies"
  - "p14_ood_state — the THREE-valued OOD resolver, the only Phase-14 function that inspects a raw verdict flag"
  - "p14_fuse — the D-05 asymmetric abstention rule, replacing SC2's literal OR"
  - "p14_abstain_reasons — the closed five-member abstention vocabulary"
affects:
  - "14-06+ (decide.jl): composes p14_fuse over p14_ood_state, p14_conformal_set and the cross-method arm"
  - "SC1-d / SC3 integration runners consume p14_conformal_calibrate + p14_hedge_diagnostics"
  - "Every Phase-14 report: SC1 and SC2 are AMENDED and the amendment must be cited alongside any result"
tech-stack:
  added: []          # zero packages; sort + ceil + integer indexing are Base
  patterns:
    - "the WRONG answer is constructed in the test (Statistics.quantile) as the reference the right one is asserted against"
    - "a docstring is CODE to a comment-stripping guard -- the corpus token ban caught exactly that"
    - "the feasibility bound is an expression in the argument n, never a literal, so a withdrawn figure cannot become a default"
    - "a duplicated vocabulary made SELF-POLICING by an isdefined-guard whose else-branch asserts agreement"
    - "the long function...end form kept deliberately, because the lane guard's exemption is POSITIONAL"
    - "every ABSTAIN carries its trigger; the reason set is closed AND asserted fully reachable"
key-files:
  created:
    - spike/p14/conformal.jl              # 401 lines
    - spike/p14/fuse.jl                   # 272 lines
    - spike/test/test_p14_conformal.jl    # 312 lines, 83 assertions
    - spike/test/test_p14_fuse.jl         # 314 lines, 222 assertions
  modified: []
decisions: [D-02, D-05, D-06]
validates: ["SC1-e", "SC2-a", "SC2-b", "SC2-c"]
metrics:
  duration: ~70 min
  completed: 2026-08-03
  tasks_completed: 2
  tasks_total: 2
  per_file_pass_count: 305   # 83 conformal + 222 fuse
---

# Phase 14 Plan 05: The Conformal Hedge and the Abstention Rule — Summary

**Both of the ROADMAP's amended success criteria are now code, and both amendments are bound by a
test that fails when the amendment is reverted — proved by reverting them and watching the failures,
not by arguing about it.** Substituting `Statistics.quantile` for the order statistic fails the
conformal gate in 6 places. Restoring SC2's literal `∨` in branch 5 fails the fusion gate in 6
places, two of them the single cell that *is* D-05.

## READ THIS BEFORE QUOTING ANY NUMBER FROM THIS PLAN

**SC1 and SC2 are AMENDED.** Any report citing a Phase-14 result produced by these two modules must
cite the **original ROADMAP wording alongside the amendment**, in both directions:

- **SC1 originally names `ConformalPrediction.jl`.** That library is **NOT added** (D-02). Conformal
  is met **in substance**, hand-rolled: `sort` + `ceil` + integer indexing. Two verified reasons —
  adding it forces a `Pkg.resolve()` against an environment a *running test* asserts byte-frozen
  (and this milestone already paid for that once, when the Phase-7 Wave-0 co-resolution gate failed
  on a transitive conflict that pushed NeuralEstimators below the pin the whole spike rests on); and
  it wraps MLJ models rather than a NeuralEstimators posterior, so it would need adapting anyway to
  supply a construction that is three lines.
- **SC2 originally reads `OOD ∨ cross-method disagreement ∨ ambiguous set ⇒ be silent`.** That
  literal OR is **replaced** (D-05) by an asymmetric rule in which **cross-method disagreement alone
  DECIDES** and is recorded as `:decided_contra_classical`. The OR would silence the tool exactly
  where Phase 9 says it should speak — disagreeing with Costes/Manders is the product thesis, not a
  failure mode — and three OR'd triggers compound into a tool that is usually silent.

"Met in substance" and "amended" are claims a reader is entitled to check. Both are cited in the
source of the file that implements them, per `14-PATTERNS` §4.6 and the `P14_DECLARED_DEVIATIONS`
enumeration frozen in `spike/p14/consts.jl`.

## Commits, in order

| # | Task | Gate | Commit | Files |
|---|---|---|---|---|
| 1 | Split conformal (LAC) | RED | `9dd9d10` | `spike/test/test_p14_conformal.jl` (+312) |
| 1 | Split conformal (LAC) | GREEN | `bf7888a` | `spike/p14/conformal.jl` (+401) |
| 2 | D-05 fusion / D-06 OOD | RED | `2c8453b` | `spike/test/test_p14_fuse.jl` (+314) |
| 2 | D-05 fusion / D-06 OOD | GREEN | `b5e1f22` | `spike/p14/fuse.jl` (+272) |

No REFACTOR commit: neither implementation needed one after its GREEN gate.

## Per-file pass counts (the aggregate suite cannot be used, 14-PATTERNS §4.8)

| File | Assertions | Testset time | Wall time | Exit |
|---|---|---|---|---|
| `test_p14_conformal.jl` | **83** | 2.7 s | 13.5 s | 0 |
| `test_p14_fuse.jl` | **222** | 1.2 s | **3.1 s** | 0 |
| `test_p14_consts.jl` (regression) | 273 | 0.2 s | — | 0 |
| `test_p14_provenance.jl` (regression) | **67** (was 63) | 2.3 s | — | 0 |
| `test_p14_decoupling.jl` (regression) | **45** (was 42) | 1.4 s | — | 0 |
| `test_p14_posterior.jl` (regression) | 48 | 2.1 s | — | 0 |
| `test_p14_fdr.jl` (regression) | 55 | 1.7 s | — | 0 |

Both lane-wide guards grew without being edited, which is the grow-with-the-phase design working:
`test_p14_provenance.jl` gained 4 assertions and `test_p14_decoupling.jl` gained 3, of which the
important ones are **SC2-c going from 1 assertion to 4** — its positive half armed itself the moment
`fuse.jl` landed, located `function p14_ood_state`'s body positionally, and confirmed the exemption
is delimitable.

## The three falsification runs, and what they showed

All three were run, and in every case the mutation was confirmed present on disk before the run and
confirmed absent after — the 14-02 lesson that a silent no-op substitution produces a check that
never ran.

**1. Conformal — the interpolated quantile substituted for the order statistic**
(`return Statistics.quantile(collect(Float64, scores), 1 - alpha)`):

```
83 assertions -> 77 passed, 6 failed, exit 1
```

The instructive part is *where* it failed: 5 of the 6 in the n = 9 discrimination testset and 1 in
the feasibility testset — and **nothing else in the file noticed at all**. The three set-size
regimes, the by-name reads, the vacuity label and the calibration shape were all still green. That
is Pitfall 6 made visible: the substitution is invisible to every test that does not target it
directly, which is exactly why it survives to production at n = 2000.

**2. Fusion — SC2's literal OR restored** (branch 5's guard reduced from
`disagree && ood_state !== :clear` to `disagree`):

```
222 assertions -> 216 passed, 6 failed, exit 1
```

Two truth-table cells, three assertions in the named-cell testset, and one reachability assertion in
the closed-reason-set testset. The amendment is bound, not described.

**3. The duplicated status vocabulary — a deliberate divergence** (`:empty` → `:vacant` in
`fuse.jl`, then both files loaded in one process):

```
ERROR: LoadError: AssertionError: spike/p14/conformal.jl: P14_CONFORMAL_STATUSES was already
defined as (:singleton, :ambiguous, :vacant), which disagrees with this file's
(:singleton, :ambiguous, :empty). ...
```

Both load orders were also verified clean before and after. The duplication is self-policing rather
than merely commented (see Deviations 5).

## What was built

### `spike/p14/conformal.jl`

`p14_conformal_quantile(scores, alpha)` is `sort(collect(Float64, scores))[ceil(Int, (n + 1) * (1 - alpha))]`.
**`Statistics.quantile` is not called anywhere in the file** — verified by grep (see Verification),
and the name is nonetheless *in scope* here, because `spike/p13/result.jl` does `using Statistics`
and the spike lane is one flat namespace. That is recorded in the banner rather than papered over:
the guard is the directed n = 9 test and the grep, **not** the import. `Statistics` is reached
qualified, for `median` in the score summary, which is the honest minimum.

It **throws** rather than clamps when `ceil((n + 1)(1 − α)) > n`, which is the deliberate difference
from `ghat` — the closest in-repo analog and one that *does* clamp. Clamping is right for a
calibration curve and wrong here: a clamped value would look like a quantile and carry no coverage
guarantee at all. The feasibility bound in the message is computed **from the argument `n`**, never
from a literal, so D-03a's withdrawn figure cannot become an executable default.

`p14_conformal_set` keeps `:empty` as its own status and the docstring says why in the terms the
report needs: an empty set means the point is more nonconforming than `(1 − α)` of the entire
calibration set — a distribution-free misspecification signal **on a completely different channel**
from the Mahalanobis/PP/noise detector, and the one signal in this phase that does not need the
simulator to be right about *densities*, only about *exchangeability*. Both regimes abstain, so
folding them would change no decision and would destroy the more informative one silently. The
arithmetic is recorded too: with three classes an empty set needs `q̂ < 2/3`, which is reachable for
this net rather than a textbook curiosity.

`p14_hedge_diagnostics` emits `vacuous_hedge` as a **LABEL beside the rates**, copying
`vacuous_pass`'s design (`spike/p13/result.jl:284-298`). It is explicitly not a failure condition,
because abstaining rarely can also be the *correct* behaviour on well-specified data; what must not
happen is a coverage number quoted without it. An **empty** outcome collection is refused rather
than answered — `0 + 0 < floor` would label a hedge that was never exercised as vacuous, the right
word for entirely the wrong reason.

The score is **LAC**, and the rejection of the BF-margin alternative is written into the docstring
because it is load-bearing: a margin score's `(1 − α)` quantile sits *by construction* in the
confident tail, and `13-REPORT.md` hands this phase a named finding that the corrected log Bayes
factors lose sub-nat precision exactly there. The posterior probability compresses that tail into
`[0,1]`, where both Phase-13 heads are ECE-green with zero empty bins.

### `spike/p14/fuse.jl`

`p14_ood_state(ood_nulls, verdict)` returns `:fired` / `:clear` / `:not_checked` and **never a
`Bool`**. Its acceptance predicate is `(t isa Real && isfinite(t))` — byte-for-byte
`_recorded_ood_threshold`'s test at `src/amortized/local_map.jl:107` — so `Inf`, `-Inf`, `NaN`,
`nothing` and a non-`Real` all resolve to `:not_checked` in agreement with the shipped semantics. A
resolver that accepted `Inf` would silently disagree while looking correct.

It is the **only** function in the Phase-14 sources that inspects a raw verdict flag, and the long
`function ... end` form is kept deliberately: `test_p14_decoupling.jl`'s SC2-c exemption is located
*positionally*, so a short-form definition would be un-delimitable and the exemption would have to
become a whole-file pass — which is not an exemption but a hole.

`p14_fuse` implements D-05's asymmetry with one early return per branch, each citing its decision in
a trailing comment, and the branch order stated as load-bearing. Every ABSTAIN carries which trigger
fired, because a reader cannot otherwise distinguish a tool that declined from a tool that broke —
nor an OOD abstention (fix the input) from a conformal-empty one (question the simulator). The
four-state alternative that would have added `DECIDED-CONTRA-CLASSICAL` as its own *action* is
recorded in the docstring as considered-and-not-chosen, with the note that it is the design to reach
for if reviewers push back on the asymmetry.

`p14_abstain_reasons()` is closed at five members, and the test asserts not only that every produced
reason is inside it but that **every member is reachable** — a closed set with a dead member is a set
that has stopped describing the rule.

## The 36-cell truth table, and why it is written out by hand

`test_p14_fuse.jl` carries all 36 rows as literal data. A test that computes its expectation by
re-implementing the rule proves only that two copies of the same mistake agree, so no expected value
in that table came from calling `p14_fuse`. The table is then asserted to **be** the complete
3 × 3 × 2 × 2 product — no cell missing, no cell twice — so it cannot pass by quietly omitting the
row that hurts. Each cell runs in its own named `@testset`, so a failure names the cell.

| `ood_state` | `conformal_status` | cells | verdict |
|---|---|---|---|
| `:fired` | any | 12 | `(:abstain, :ood_fired)` — the detector answered and said no |
| `:clear` | `:singleton` | 4 | `(:decide, :decided)` / `(:decide, :decided_contra_classical)` |
| `:clear` | `:ambiguous` / `:empty` | 8 | `(:abstain, :conformal_ambiguous)` / `(:abstain, :conformal_empty)` |
| `:not_checked`, no opt-out | any | 6 | `(:abstain, :ood_not_checked)` — D-06 |
| `:not_checked`, opt-out | `:singleton` | 2 | `(:decide, :decided)` / `(:abstain, :disagreement_and_ood)` |
| `:not_checked`, opt-out | `:ambiguous` / `:empty` | 4 | the conformal trigger, which precedes branch 5 |

**The one cell that is the whole amendment** is asserted separately by name:
`(:clear, :singleton, disagree = true, allow_unchecked_ood = false) ⇒ (:decide,
:decided_contra_classical)`, together with its no-disagreement twin, so the reason field is shown to
be *carrying* the disagreement rather than being the constant answer for a clear singleton.

## Verification (observed, not inferred)

| Check | Command | Result |
|---|---|---|
| SC1-e | `julia --project=spike spike/test/test_p14_conformal.jl` | **exit 0**, 83/83 |
| SC2-a + SC2-b | `julia --project=spike spike/test/test_p14_fuse.jl` | **exit 0**, 222/222 |
| SC2-c armed | `julia --project=spike spike/test/test_p14_decoupling.jl` | **exit 0**, 45/45 (SC2-c 1 → 4) |
| Pre-registration | `julia --project=spike spike/test/test_p14_consts.jl` | **exit 0**, 273/273 |
| D-07 provenance | `julia --project=spike spike/test/test_p14_provenance.jl` | **exit 0**, 67/67 |
| Class order | `julia --project=spike spike/test/test_p14_posterior.jl` | **exit 0**, 48/48 |
| SC1-a | `julia --project=spike spike/test/test_p14_fdr.jl` | **exit 0**, 55/55 |
| Frozen surfaces | `git diff --exit-code HEAD -- src spike/Project.toml spike/Manifest.toml corpus` | **exit 0** |
| No frozen p14 module touched | `git diff --exit-code HEAD~4 HEAD -- spike/p14/{consts,provenance,posterior,fdr}.jl spike/test/test_p14_decoupling.jl` | **exit 0** |
| Only the four allowed files changed | `git diff --name-only HEAD~4 HEAD` | exactly the four |
| No file deleted by any of the four commits | `git diff --diff-filter=D --name-only HEAD~4 HEAD` | empty |
| No stray files under the lane | `git status --porcelain -- spike src corpus artifacts` | empty |
| The interpolating quantile is absent | `grep -v '^\s*#' spike/p14/conformal.jl \| grep -nE '(^\|[^_[:alnum:]])quantile\(\|Statistics\.quantile'` | **0 hits** (see Deviation 1) |
| The withdrawn D-03a figure is absent | `grep -v '^\s*#' spike/p14/conformal.jl \| grep -cE '0\.032\|1/31'` | **0** |
| Only `p14_ood_state` reads a raw flag | `grep -v '^\s*#' spike/p14/fuse.jl \| grep -c '\.flag'` | **1**, and it is `return verdict.flag ? :fired : :clear` |
| The three-valued state is named | `grep -c 'not_checked' spike/p14/fuse.jl` | **13** (≥ 3 required) |
| The named amendment cell is asserted | `grep -n 'decided_contra_classical' spike/test/test_p14_fuse.jl` | 6 hits incl. the table row and the standalone assertion |
| Both modules co-load in EITHER order | `include` both, both orders | OK; `P14_CONFORMAL_STATUSES` identical |

## Deviations from Plan

### 1. [Plan defect — the acceptance grep is self-defeating] The Task-1 quantile grep cannot pass as written

- **Found during:** Task 1, running the acceptance criteria.
- **Issue:** the criterion is
  `grep -v '^\s*#' spike/p14/conformal.jl | grep -c 'Statistics.quantile\|quantile('` **returns 0**.
  But the same plan's `must_haves` require the file to *provide* `p14_conformal_quantile`, and the
  string `p14_conformal_quantile(` **contains** `quantile(`. The literal criterion is therefore
  unsatisfiable by any file that satisfies the plan's own artifact contract. Observed: the literal
  grep returns **3**, and all three hits are the function's own docstring signature, its `function`
  line, and its one internal call site.
- **This is the third instance of the recorded pattern** — a guard that must NAME what it forbids
  trips its own scan, because strings in code are not comment-stripped. `test_p14_decoupling.jl`
  answers it by assembling needles from fragments; a *name* cannot be assembled, so the answer here
  has to be on the grep side.
- **Fix:** the substantive check — no `quantile(` call that is not part of `p14_conformal_quantile`
  — run as
  `grep -v '^\s*#' spike/p14/conformal.jl | grep -nE '(^|[^_[:alnum:]])quantile\(|Statistics\.quantile'`,
  which returns **0 hits**. The binding property is unchanged and is additionally proved by
  falsification run 1. **No code was weakened to accommodate the grep.**
- **Commit:** `bf7888a`.

### 2. [Rule 1 — caught by a running guard] A docstring named a banned token, and the lane guard failed

- **Found during:** Task 1, first regression run of `test_p14_decoupling.jl`.
- **Issue:** `test_p14_decoupling.jl` failed with
  `token_hits == ["spike/p14/conformal.jl: corpus"]`. The first draft of
  `p14_conformal_quantile`'s docstring explained the withdrawn D-03a floor and named the token doing
  it. **Comments are stripped by that scan; docstrings are not** — a docstring is code text.
- **Why this is worth recording rather than quietly fixing:** it is Pitfall 8 ("the token ban is a
  code ban, not a prose ban") landing on this plan, and it is the *inverse* of the failure mode the
  pitfall warns about. An executor who reacted by deleting every mention would have broken the
  guard's own preconditions; the correct move is to put the explanation where the scan permits it.
- **Fix:** the full explanation moved into a `#` comment block immediately above the function — a
  more prominent place, not a lesser one — and the docstring now cross-references it. The comment
  itself records why it lives there.
- **Files:** `spike/p14/conformal.jl`. **Commit:** `bf7888a`.

### 3. [Test correction during the GREEN gate] Two test-side claims about the resolver were wrong

- **Found during:** Task 2, first GREEN run — 1 fail and 1 error against a 36/36 green truth table.
- **(a)** `p14_ood_state(Dict(:thr => 3.5), …)` errored with `FieldError: type Dict has no field
  thr`. The **implementation is correct**: it mirrors the shipped call site's
  `haskey(ood_nulls, :thr) ? ood_nulls.thr : nothing` byte-for-byte, and `ood_nulls` is NamedTuple-
  shaped in the shipped path. Broadening the resolver to `get(ood_nulls, :thr, nothing)` so a Dict
  would work would be a *silent divergence from the shipped semantics* — precisely the class of
  divergence D-06 exists to prevent. The test was corrected to assert the Dict fails **loudly**, and
  that a container with no threshold at all still resolves conservatively to `:not_checked`.
- **(b)** A `Bool` passed as `ood_state` raises `TypeError`, not `MethodError` — it is refused by the
  *signature*, before a single branch runs, which is a stronger guarantee than any body check. The
  test was corrected and a second signature-level case (`disagree = :yes`) added.
- **Files:** `spike/test/test_p14_fuse.jl`. **Commit:** `b5e1f22`.

### 4. [Rule 2 — missing validation] Six input checks the plan did not specify

- `_p14_check_scores` rejects a non-finite calibration score. The damage is specific and silent:
  `sort` puts `NaN` at the END of the ascending order, so one `NaN` makes the top order statistic
  `NaN`, and `NaN <= qhat` is false for every class — **the conformal set comes back EMPTY for every
  input**, and an empty set is a result this phase reports as a misspecification *finding*. An
  arithmetic fault would be indistinguishable from the most interesting output the hedge can produce.
- `p14_conformal_quantile` rejects `n = 0` explicitly, so an empty calibration set is not reported
  through the infeasibility message (which would send a reader to the wrong question).
- `p14_conformal_set` and `p14_hedge_diagnostics` reject a non-finite `qhat`: it silently admits
  every class or none, and both look like a legitimate regime.
- `p14_conformal_calibrate` rejects a length mismatch between posteriors and labels — `zip` would
  otherwise truncate silently and shift the quantile.
- `_p14_set_status` validates against `P14_CONFORMAL_STATUSES`: a mistyped symbol would be counted
  into none of the three rates, and the rates would sum to less than 1 with nothing saying why.
- `p14_hedge_diagnostics` refuses an empty collection (see above).
- **Commit:** `bf7888a`.

### 5. [Design decision within Claude's discretion] `P14_CONFORMAL_STATUSES` is duplicated, and made self-policing

- **The problem:** `fuse.jl` must validate `conformal_status`, and `conformal.jl` produces it. One
  vocabulary, two files. The plan's include list for `fuse.jl` is `consts.jl` + `src/results.jl`
  only, and `consts.jl` is FROZEN so the tuple cannot be added there.
- **Why not simply include `conformal.jl` from `fuse.jl`:** `conformal.jl` pulls the entire
  Phase-13 stack (Flux, NeuralEstimators, Images) through `posterior.jl` — about 10 s of package
  load in this environment. `test_p14_fuse.jl` currently runs in **3.1 s wall**; the include would
  have made it ~14 s for a three-symbol tuple, trading a fast, cheap check for a slow one on a rule
  that is pure symbol arithmetic.
- **What was done instead:** both files declare the tuple under the same name behind the same
  `isdefined` guard, and whichever loads **second asserts agreement** rather than silently losing.
  A textual divergence therefore fails loudly the moment both are in one process — which
  `decide.jl` guarantees. Falsification run 3 confirms the assertion fires and names both tuples;
  both load orders were verified clean.
- **Commits:** `bf7888a`, `b5e1f22`.

### 6. [Shape choice within Claude's discretion] `p14_fuse` returns a NamedTuple, not a positional Tuple

- The plan's action text specifies `-> NamedTuple` while `14-RESEARCH`'s sketch returns a bare
  `(Symbol, Symbol)`. The NamedTuple `(action = …, reason = …)` was chosen: it is the phase's
  discipline (every read by name, no positional hazard) and it matches every other Phase-14
  primitive. Consequence for readers: `r == (:decide, :decided)` is **false**; the test compares
  `(r.action, r.reason) === (…, …)`.

### 7. [Additive] `fuse.jl` declares three `const` tuples the plan listed only as function bodies

`P14_OOD_STATES`, `P14_ABSTAIN_REASONS` and `P14_DECIDE_REASONS` are named constants inside the
guard block, so the validation message enumerates the legal states mechanically and
`p14_abstain_reasons()` returns a frozen binding rather than a fresh literal per call. No behaviour
differs; the enumerations become greppable.

### Assumption Drift (advisory)

**1. "Both test files complete in under 10 seconds" — true for one, and the reason the other is
slower is structural, not incidental.**

- **Planned:** each unit file under 10 s.
- **Actual:** `test_p14_fuse.jl` **3.1 s wall** (1.2 s of testset). `test_p14_conformal.jl`
  **13.5 s wall** (2.7 s of testset).
- **Why:** `conformal.jl` reaches the Phase-13 surfaces through `posterior.jl`, which loads Flux,
  NeuralEstimators and Images — a fixed ~10 s paid before a single assertion runs, and the same
  shape `14-04` already recorded. It is not repairable without severing the plan's own `key_links`
  requirement that `conformal.jl → posterior.jl`.
- **The advisory part:** `fuse.jl` was *deliberately kept out of that dependency* (Deviation 5), and
  the measured 4× difference is the evidence that the choice mattered. The substantive claim — no
  simulation ran, the primitives are pure arithmetic, feedback is fast — holds: **305 assertions in
  3.9 s of combined testset time.**

**2. The plan's behaviour examples read as bare tuples; the implementation returns a NamedTuple.**

- **Planned:** `p14_fuse(...) returns (:abstain, :ood_not_checked)`.
- **Actual:** `(action = :abstain, reason = :ood_not_checked)`, per the same plan's `-> NamedTuple`
  signature line.
- **Why it is recorded:** the two statements in the plan are not literally consistent, and a later
  reader comparing this file against the plan's prose would otherwise see a mismatch. The signature
  line was followed because it agrees with the phase's by-name discipline.

## Known Stubs

None. Both modules are complete implementations of their specified behaviour; nothing returns a
placeholder and nothing is wired to an empty data source. The two functions that *look* like they
await data — `p14_hedge_diagnostics` and `p14_conformal_calibrate` — are complete transforms whose
inputs the SC1-d integration runner will supply.

## Threat Flags

None. No file created here opens a network endpoint, reads an image, touches an auth path, or
changes a schema at a trust boundary. `conformal.jl` is pure arithmetic over already-loaded numbers;
`fuse.jl` is pure symbol arithmetic with no file-system reach at all.

## Threat register, discharged

| ID | Threat | How it is now mitigated |
|---|---|---|
| T-14-08 | `verdict.flag == false` read as an in-distribution proxy | `p14_ood_state` is three-valued and never returns a `Bool`; it is the ONLY function that touches a raw flag (grep = 1, inside its own body); SC2-c's positional exemption located it and armed 3 further assertions; `:not_checked` abstains by default across 6 of the 36 cells, with an opt-out a caller must set deliberately |
| T-14-09 | `Statistics.quantile` substituted for the order statistic | Implementation is `sort(...)[ceil(Int, (n + 1) * (1 - alpha))]`; boundary-anchored grep finds zero interpolating calls; the n = 9 test asserts `== maximum(s)` AND `!= Statistics.quantile(s, 0.9)` AND the direction; **falsification run fails 6** |
| T-14-05 | adding `ConformalPrediction.jl` | Zero packages added; `spike/Project.toml` and `Manifest.toml` byte-unchanged, asserted by the running lane guard, which also asserts `!haskey(..., "ConformalPrediction")` and `!haskey(..., "MLJ")` |
| T-14-21 | SC2's literal OR reinstated, silencing the Phase-9 thesis | 36-cell truth table from an INDEPENDENT hand-written expectation table, asserted to be the complete product once each; the `(:clear, :singleton, true, false)` cell asserted by name with its twin; **falsification run fails 6** |
| T-14-22 | the empty conformal set folded into `:ambiguous` | Separate `:empty` status and separate `:conformal_empty` reason; a dedicated fixture asserts `:empty` is not reported as `:ambiguous`; `p14_abstain_reasons()` is closed at five AND every member asserted reachable |
| T-14-23 | a vacuous hedge reported as a pass | `p14_hedge_diagnostics` emits `vacuous_hedge` as a LABEL beside the three rates; the test asserts it is returned and not thrown, and that a hedge which DOES fire is not labelled |
| T-14-SC | package installs | Zero installs. `sort` + `ceil` + integer indexing for the hedge; symbol comparison for the fusion |

## What this does NOT establish

- **No coverage number exists yet.** This plan built the *construction*; SC1-d measures whether the
  realized marginal coverage clears `1 − α_conf − 3·√(α(1−α)/n_eval) = 0.880`, and that runs on
  fresh unstratified simulator draws in a later plan.
- **`p14_conformal_calibrate` cannot check its own precondition.** The guarantee needs the
  calibration points to be exchangeable with the test point. A class-BALANCED calibration set drawn
  against a prior-distributed deployment stream is **not** exchangeable, and the realized coverage
  will then miss `1 − α` silently. Nothing in that signature can see how the scores were drawn; the
  unstratified-draw assertion belongs to the runner, and this is written into the docstring so it is
  not assumed to have been handled here.
- **`:clear` means ONE channel says clear.** The Phase-14 null is a single-channel density fit; the
  noise and posterior-predictive channels are not wired in this lane. Stated in the `p14_ood_state`
  docstring rather than left to be discovered.
- **No decision has been made about any real or simulated dataset.** `p14_fuse` is the rule;
  `decide.jl` is the composition point, and the cross-method `disagree` input it consumes does not
  exist yet.
- **The headline conformal guarantee will be simulator-derived** and therefore inherits the
  simulator's misspecification (D-03, unchanged by this plan). That is a named manuscript limit.

## Self-Check: PASSED

| Claim | Check | Result |
|---|---|---|
| `spike/p14/conformal.jl` exists | `[ -f ]` | FOUND |
| `spike/p14/fuse.jl` exists | `[ -f ]` | FOUND |
| `spike/test/test_p14_conformal.jl` exists | `[ -f ]` | FOUND |
| `spike/test/test_p14_fuse.jl` exists | `[ -f ]` | FOUND |
| `9dd9d10` in history | `git log` | FOUND |
| `bf7888a` in history | `git log` | FOUND |
| `2c8453b` in history | `git log` | FOUND |
| `b5e1f22` in history | `git log` | FOUND |
