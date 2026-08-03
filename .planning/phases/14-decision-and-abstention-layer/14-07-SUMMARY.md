---
phase: 14
plan: 07
subsystem: decision-composition
status: COMPLETE
tags: [decide-coloc, abstain-first, assumption-a2, d-01, d-04, d-05, d-06, d-07, sc2-d, pitfall-7, wave-4]
requires:
  - "14-01 — spike/p14/consts.jl (P14_DEV_SEED, P14_PI_COLOC_GRID, P14_CLASS_KEYS, P14_FIXTURE_COUNTER)"
  - "14-02 — spike/p14/provenance.jl (p14_load_tau, p14_provenance_record, _p14_load_report)"
  - "14-03 — spike/test/test_p14_decoupling.jl (the lane guard whose SC2-d and costes-seed positive halves this plan ARMS)"
  - "14-04 — spike/p14/{posterior,fdr,conformal}.jl (p14_class_posterior, p14_fdr_over_decided, p14_conformal_set, p14_hedge_diagnostics)"
  - "14-05 — spike/p14/fuse.jl (p14_ood_state, p14_fuse, p14_abstain_reasons)"
  - "14-06 — spike/p14/result.jl (P14Result, P14BatchDecision, p14_batch_decision, p14_headline, p14_named_limits)"
  - "14-08 — spike/p14/pools.jl (p14_ood_reference; consumed through a NEW adapter, see deviation 1)"
  - "spike/p13/{net,labels,preconditions,result}.jl + three_way_net.jld2 + three_way_gate_report.jld2 — READ-ONLY"
  - "spike/comparator/{config,classical}.jl, spike/simulator/ghat.jl — frozen, READ-ONLY"
  - "src/{results.jl,registry.jl,amortized/{ood,local_map}.jl} — READ-ONLY, zero bytes edited"
provides:
  - "p14_load_bundle — the inherited Phase-13 handle: tau LOADED once (D-07, no override keyword), the MEASURED class prior, and Phase 13's own CalibrationMeta carried forward"
  - "p14_cross_method — the two D-05 channels on the ghat rho_true scale, Costes on P14_DEV_SEED, Manders descriptive only"
  - "p14_decide_one — one image pair to one P14Result"
  - "decide_coloc(bundle, img, control, channels; ...) -> P14Result — the SRC-SHAPED per-pair entry (D-01)"
  - "decide_coloc(bundle, pairs, ood_nulls, qhat; alpha_fdr, ...) -> P14BatchDecision — ABSTAIN FIRST, THEN SORT"
  - "p14_ood_input — the adapter between p14_ood_reference's return and the shape the SHIPPED ood_verdict reads"
  - "Assumption A2 named, derived and bound by a test that fails if the order is reversed"
affects:
  - "every later Phase-14 runner: they call decide_coloc rather than assembling the rule"
  - "the phase report: it prints p14_headline(b) and enumerates b.meta.named_limits"
tech-stack:
  added: []          # zero packages; every call is into already-frozen in-repo code
  patterns:
    - "one code path for the headline and the sensitivity sweep, with the two paths' agreement ASSERTED — a second simpler loop would agree until the day it stopped"
    - "a fixture whose every threshold is DERIVED from its own draws, so nothing goes stale silently"
    - "falsify the guard, then close what the falsification found: the SC2-d value assertion exists because the label assertion was proven insufficient"
    - "a run-time forbidden-seed refusal beside the review-time grep — the grep checks spellings, the assertion checks values"
    - "the batch rule is authoritative for the coloc call; the provisional argmax call is PRESERVED so a divergence is auditable rather than erased"
key-files:
  created:
    - spike/p14/decide.jl             # 1049 lines
    - spike/test/test_p14_decide.jl   # 445 lines, 146 assertions
  modified: []
decisions: [D-01, D-04, D-05, D-06, D-07]
validates: ["SC1-a", "SC2-a", "SC2-b", "SC2-d"]
metrics:
  duration: ~80 min
  completed: 2026-08-03
  tasks_completed: 3
  tasks_total: 3
  per_file_pass_count: 146
  test_wall_clock_s: 25.4
  test_testing_time_s: 5.6
---

# Phase 14 Plan 07: `decide_coloc` — the Decision and Abstention Entry Point — Summary

**The ordering claim is now code, and it is bound by a test that a reordering turns red.** The
batch is partitioned into decided and abstained **before anything is sorted**; only the survivors'
null posteriors reach the prefix rule. Assumption A2 — the rejection set is a deterministic
function of the observed data, hence σ(data)-measurable, hence it pulls out of the conditional
expectation for *any* data-dependent selection — is named as such in the source, derived in two
lines, and stated alongside why the reverse order is wrong, so neither has to be re-derived by
someone who assumes the order was arbitrary.

## SAY THIS PLAINLY IN ANY REPORT: "SRC-SHAPED" IS NOT "IN SRC"

`decide_coloc` takes the signature it **would** have in `src/` — a bundle plus two images and the
two channels to compare, returning a subtype of the shipped `AbstractColocResult` — so that a
later promotion would be mechanical rather than a rewrite (D-01). **It is built entirely in
`spike/` and is NOT shipped in v2.0.** This plan edits zero bytes of `src/`;
`git diff --exit-code HEAD -- src spike/Project.toml spike/Manifest.toml corpus spike/validation`
exits 0, `git status --porcelain -- src` is empty, and both facts are asserted inside the test
file as well as in the lane guard. `spike/test/test_p14_decide.jl` also asserts
`!isfile("src/amortized/decide.jl")` — the file this layer would become if anyone confused the two
readings.

## What was built

### `spike/p14/decide.jl` (1049 lines)

**`p14_load_bundle`** — the inheritance, loaded once. τ arrives through `p14_load_tau` with its
four-way on-disk agreement and its two-sided sha table asserted, and the net's own recorded value
is asserted equal to it again on this side. **There is no keyword by which a caller can override
τ**, and the test asserts the keyword set is exactly `(:net_path, :probe_path, :gate_path)`:
copying the value is how two numbers silently diverge, re-deriving it guarantees a second
different cut on the same hypothesis space, and no legitimate caller in this phase needs a third
option (D-07). The bundle also carries the **measured** `pi_class_masses` and Phase 13's own
`CalibrationMeta`.

**Phase 13's calibration is CARRIED, not recomputed.** `_p14_carry_p13_calibration` reads the
coloc head's reliability curve, ECE, MCE and AUC back off `three_way_gate_report.jld2` and pushes
them through the *shared* `p13_calibration_meta` builder. Re-running `_bin_calibration` on
Phase-14 draws would have produced a second, differently-sourced ECE that a reader would
inevitably conflate with the gated one. The carry is not taken on trust: the builder re-derives
the empty-bin count, the vacuity label and the band verdict from the carried curve, and all three
are asserted equal to the values Phase 13 persisted — so a carry that landed on the wrong head or
the wrong bins fails here instead of travelling.

**`p14_cross_method`** — the two D-05 channels. Channel 1 applies `ghat` to the patch-correlation
of **both** acquisitions before the inherited cut, because the cut is defined on ρ_true while the
classical statistic returns the induced mean per-patch correlation μ, and comparing them directly
would be a criterion applied in a unit it was not derived for — one of the four recorded "the bar
was wrong" families in this project. Channel 2 is the Costes block-permutation test, genuinely
independent rather than a second correlation cut. Manders travels as a **descriptive column with
no call derived from it**, because a cut on a fraction-of-signal quantity would be a new chosen
constant. `basis = :ghat_rho_true_scale` is recorded on every computed record.

**`p14_decide_one`** — net forward pass → class posterior → three-valued OOD state → conformal set
→ amortized call **by name** → D-05 fusion → one `P14Result` built once at the end with a
provenance meta. The evidence triple is normalised through `_p13_as_three_way_logbf` rather than
built positionally, and the argmax is a three-term comparison over named fields with an explicit
tie rule (a tie falls to `RANDOM`, mirroring the frozen label rule's dead zone and erring toward
the null).

**`decide_coloc`, batch** — abstain first, then sort, then map the accepted positions back to
original batch indices and **assert** none of them abstained. The prior-sensitivity curve is
required output and re-runs the *whole* rule at every pre-registered rung through the **same code
path** as the headline — necessary rather than fastidious, because the prior moves the conformal
set (hence the abstention) and the amortized call (hence the cross-method disagreement) as well as
the null posterior. The two computations of the fusion, one on the per-pair path and one on the
batch path, are asserted to agree item by item.

### `spike/test/test_p14_decide.jl` (445 lines, 146 assertions, 5.6 s testing / 25.4 s wall)

Eleven testsets. Two of them exist in their present form only because a weaker version was
**falsified against a deliberately broken implementation** (see *Falsification*, below).

## Verification (observed, not assumed)

Every command below was run and its exit status or output observed.

| Command | Result |
|---|---|
| `julia --project=spike spike/test/test_p14_decide.jl` | **EXIT 0** — 146 passed, 0 failed, 25.4 s wall (5.6 s testing) |
| `julia --project=spike spike/test/test_p14_decoupling.jl` | EXIT 0 — **48/48** (was 45/45; the SC2-d and costes-seed positive halves are now ARMED) |
| `julia --project=spike spike/test/test_p14_consts.jl` | EXIT 0 — 273 |
| `julia --project=spike spike/test/test_p14_provenance.jl` | EXIT 0 — 72 |
| `julia --project=spike spike/test/test_p14_posterior.jl` | EXIT 0 |
| `julia --project=spike spike/test/test_p14_fdr.jl` | EXIT 0 — 55 |
| `julia --project=spike spike/test/test_p14_conformal.jl` | EXIT 0 — 83 |
| `julia --project=spike spike/test/test_p14_fuse.jl` | EXIT 0 — 222 |
| `julia --project=spike spike/test/test_p14_result.jl` | EXIT 0 — 121 |
| `julia --project=spike spike/test/test_p14_pools.jl` | EXIT 0 |
| plan Task-1 `<automated>` | `bundle ok, prior=(exclusion = 0.31272, random = 0.47073, coloc = 0.21655)` |
| plan Task-2 `<automated>` (substituted, deviation 5) | both `decide_coloc` methods present; `P14Result <: AbstractColocResult` |
| `git diff --exit-code HEAD -- src spike/Project.toml spike/Manifest.toml corpus spike/validation` | **EXIT 0** |
| `grep -c 'ghat(patch_correlation(' spike/p14/decide.jl` | **3** (≥ 2 required) |
| `grep -v '^\s*#' … \| grep -c 'costes_p(MASTER_SEED\|costes_p(COSTES'` | **0** |
| `grep -v '^\s*#' … \| grep -c '321'` | **0** — the input width is derived |
| `grep -v '^\s*#' … \| grep -c 'is_ood('` | **0** — the lossy Bool accessor is never routed through |
| `grep -v '^\s*#' … \| grep -c 'alpha_fdr = P14_ALPHA_CONFORMAL'` | **0** (D-07) |
| `grep -c 'A2' spike/p14/decide.jl` | **2** |
| `grep -c 'pi_sensitivity\|pi_grid' spike/p14/decide.jl` | **6** |
| FDR call sites | exactly ONE, `p14_fdr_over_decided(v, …)` with `v` the decided-subset vector |

`STATE.md` and `ROADMAP.md` were **not** touched: the orchestrator owns those writes.

## Falsification (every guard was broken on purpose, and the mutation was verified to land)

| Mutation | Verified landed | Result |
|---|---|---|
| `_p14_rule_at_prior` reordered to **sort-then-abstain** (prefix over the whole batch, then drop abstained) | `git diff` showed +5/−1 | `test_p14_decide.jl` **EXIT 1** — the ordering testset failed at the k\* equality, the accepted-index mapping and every discriminating level |
| `p14_decide_one`'s `allow_unchecked_ood` default flipped to `true` | `git diff --stat` 1/1, line 525 confirmed | D-06 testset **3 failures** |
| `ghat` dropped from the **sample** channel of `p14_cross_method` | `git diff` showed the exact line | **First attempt: BOTH guards stayed green.** See below. |

**The third mutation is the finding of this plan.** Dropping the rescaling from *one* of the two
classical channels was caught by **neither** the lane guard **nor** the first draft of this test:

- the lane's SC2-d rule is per-line and fires only on a line mentioning both the statistic and the
  cut — and here they live on different lines, so a stripped `ghat` is invisible to it;
- its positive half (`occursin("ghat(patch_correlation(", code)`) still matched, because the
  *control* channel kept its call;
- a range check passes, because a raw induced mean also lies in `[-1, 1]`;
- and `basis = :ghat_rho_true_scale` is a **literal**, which a wrong computation carries just as
  happily as a right one.

The test now recomputes the rescaled value from the same frozen ingestion on **both** channels,
asserts equality, asserts the rescaling actually **moved** the number (`|ĝ(μ) − μ| > 1e-6`; the
observed shift on the fixture is 0.114 / 0.124), and asserts the classical call is the frozen rule
applied at the inherited cut on that scale. Re-running the mutation against the closed version:
**EXIT 1**, 2 failures. The lane guard's own weakness is left as it is and recorded here rather
than widened — editing a frozen Phase-14 module is out of this plan's scope, and the value-level
assertion is the stronger check anyway.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 — Blocking] `p14_ood_reference(...).ood_nulls` is NOT the shape the shipped
`ood_verdict` reads, and the mismatch is silent**
- **Found during:** Task 1
- **Issue:** `ood_verdict` reads its density fit as `ood_nulls.density` (`src/amortized/ood.jl:369`),
  while `p14_ood_reference` (14-08, frozen) returns the raw `fit_ood_nulls` tuple `(cont, μS, C)`
  merged with `(thr = …)`. There is no `:density` key on it. `test_p14_pools.jl` cannot catch this:
  it asserts the `:thr` key and calls `p14_ood_state`, which reads only `:thr` — it never calls
  `ood_verdict` with the object.
- **Fix:** `p14_ood_input(ref)` in `decide.jl` adapts the reference into verdict shape, carrying
  `:thr` **only** when an operating point was actually fitted, so an unfitted null keeps resolving
  to `:not_checked` and therefore ABSTAINS. `spike/p14/pools.jl` is **not edited** — it is a frozen
  Phase-14 module, and the adapter belongs at the composition point where the two surfaces meet.
- **Verified:** the decide test drives both shapes (`DEC_OOD` with the operating point, and
  `DEC_OOD_UNCHECKED` without) and asserts `:fired`/`:clear` in the first case and `:not_checked`
  in the second.
- **Commit:** 2034573

**2. [Rule 1 — Bug in the plan's own suggested code] `encode_conditioned_pair(Zs, Zc, [lambda])`
bypasses Phase 11's conditioning encoder**
- **Found during:** Task 1
- **Issue:** the plan's action text says to build the net input as
  `encode_conditioned_pair(Zs, Zc, [lambda])`. That appends the **raw** λ where the net expects
  Phase 11's affine `[0,1]` map of it. `p13_encode_lambda`'s own docstring names this exact hazard:
  *"re-deriving it here — rescaling, clamping, standardizing, or 'improving' it — would feed the net
  a conditioning value it never saw at a row it reads, and nothing would throw."* At the reference
  λ = 3.0 against a trained range of (0.25, 3.0) the raw value is 3.0 where the encoder yields 1.0.
- **Fix:** `p13_encode_pair(Zs, Zc, lambda)` is used instead — the sanctioned route, which routes
  through the inherited encoder and re-asserts the assembled width against `three_way_input_dim`.
- **Commit:** 2034573

**3. [Rule 2 — Missing critical] `ood_verdict`'s unchecked precondition is guarded**
- **Found during:** Task 1
- **Issue:** `ood_verdict` raises on a null that carries no density fit. In this phase a missing fit
  is not an error — it is the NOT-CHECKED state — so a batch whose correct answer is "abstain,
  because nothing checked this" would instead crash.
- **Fix:** `_p14_ood_verdict` returns `nothing` when no fit is present and the caller resolves that
  to `:not_checked`. The shipped function is still consumed unchanged; only its precondition is
  checked before the call.
- **Commit:** 2034573

**4. [Rule 2 — Missing critical] a run-time refusal beside the Costes seed grep**
- **Found during:** Task 1
- **Issue:** the lane grep bans the forbidden *spellings* `costes_p(MASTER_SEED` and
  `costes_p(COSTES`. It cannot catch a forbidden *value* arriving through a variable, which is
  exactly how a caller would pass one.
- **Fix:** `_p14_assert_costes_seed` throws on any seed in `_p13_forbidden()` before the call.
  Asserted live in the test with `first(_p13_forbidden())`.
- **Commit:** 2034573

**5. [Rule 2 — Missing critical] the batch rule is authoritative for the coloc call, and the
provisional call is PRESERVED**
- **Found during:** Task 2
- **Issue:** the plan says an accepted item "keeps `:coloc`". At a permissive level the FDR prefix
  can accept an item whose argmax is not the coloc class, so "keeps" is not always what happens —
  and silently overwriting the argmax call would erase a real disagreement.
- **Fix:** the rule is authoritative (a discovery is defined by the accepted set, which is what the
  controlled quantity is a property of), and `_p14_finalise` records the provisional call in
  `meta.provisional_decision` beside `meta.accepted_by_fdr_rule`. Results are **rebuilt through the
  constructor**, never mutated, so the finalised object is proven still legal.
- **Commit:** bc3197b

**6. [Rule 2 — Missing critical] one code path for the headline and the sensitivity sweep**
- **Found during:** Task 2
- **Issue:** the plan describes the sweep as "recompute every `class_posterior` and re-run steps
  2-3". Written as a second loop, that would drift from the headline path and agree with it right
  up until the day it stopped, with nothing to say which was wrong. It also under-states the work:
  the prior moves the conformal set and the amortized call (hence the cross-method disagreement),
  not only the null posterior.
- **Fix:** `_p14_rule_at_prior` computes the whole rule at one prior and is called by both the
  headline and every rung; `_p14_recall_disagreement` re-derives disagreement from the two stored
  classical verdicts against the sweep's own amortized call. The per-pair and batch computations of
  the fusion are asserted to agree item by item.
- **Commit:** bc3197b

**7. [Rule 2] the `local_map` keyword is named `local_map_uncontrolled`**
- **Found during:** Task 1
- **Issue:** the plan's signature names it `local_map`, which both shadows the accessor of that name
  in the lane's flat namespace and drops the D-04 label at the one place a caller types it.
- **Fix:** the keyword matches the field name. D-04's whole design is that the label lives in the
  type rather than only in a docstring; extending that to the call site costs nothing.
- **Commit:** 2034573

### Plan text that could not be satisfied literally

**8. Task 1's behaviour list and its action text contradict each other on `cross_method.basis`.**
The behaviour says `basis === :ghat_rho_true_scale` on **every** returned result; the action says a
result with no images gets `(disagree_any = false, basis = :not_computed, …)` with a note. The
action was implemented: recording `:ghat_rho_true_scale` for a comparison that never happened would
be a false provenance record, which is precisely what `basis` exists to prevent. The test asserts
both arms — `:ghat_rho_true_scale` on the live image path, `:not_computed` with a `note` when no
images were supplied.

**9. Task 2's `<automated>` verify cannot pass with typed signatures, and the types are worth
more.** `hasmethod(decide_coloc, Tuple{Any,Any,Any,Any})` is `false` unless all four positional
arguments are untyped. Typing them (`::NamedTuple`, `::AbstractVector` vs `::MultiChannelImage`) is
what makes the two entry points dispatch unambiguously and what gives a caller a real error rather
than a validation message from inside the body. The check was replaced with the strictly stronger
pair — `hasmethod(decide_coloc, Tuple{NamedTuple, MultiChannelImage, MultiChannelImage, Vector{Int}})`
and `P14Result <: AbstractColocResult` — which asserts the src-shape claim itself rather than the
existence of a four-argument method.

**10. The Costes call line carries `P14_DEV_SEED` in a trailing comment.** The lane guard requires
every `costes_p(` line to name the Phase-14 stream, and the call passes the `seed` parameter (whose
default *is* `P14_DEV_SEED`). Whole-line comments are stripped by the scan but trailing ones are
not, so the trailing note satisfies it truthfully. The real guard on that call is the run-time
refusal of deviation 4; the comment is the review-time signpost.

**11. `pi_grid` is a keyword the batch constructor pins.** The plan exposes `pi_grid`, but
`P14BatchDecision` asserts the curve sweeps exactly `P14_PI_COLOC_GRID` in order. Passing anything
else therefore fails loudly at construction rather than quietly producing an off-grid curve. Kept
as written and documented, because failing loudly is the correct behaviour for a pre-registered
grid.

## Assumption Drift (advisory)

**1. "src-shaped" lands on a per-pair method, not on the batch one.**
- **Planned:** D-01 and the plan objective describe `decide_coloc` as taking a bundle plus images
  and returning an `AbstractColocResult` subtype; Task 2's signature is
  `decide_coloc(bundle, pairs, ood_nulls, qhat; alpha_fdr, …) -> P14BatchDecision`.
- **Actual:** `P14BatchDecision` is deliberately **not** an `AbstractColocResult` (14-06 built it as
  a batch object), and Task 1's own behaviour list plus Task 3's testset 7 both require a `channels`
  argument that only an image-taking method can have. So `decide_coloc` has **two** methods: the
  src-shaped per-pair one returning `P14Result <: AbstractColocResult`, and the batch one returning
  `P14BatchDecision`. FDR is a batch quantity and a single pair cannot express it, so this is the
  only shape that satisfies both D-01 and SC1.
- **Why it matters:** a reader checking the D-01 claim must look at the per-pair method. The test
  asserts exactly that, so the claim is checkable rather than inferable.

**2. λ defaults rather than being required.** The plan makes `lambda` a required keyword of
`p14_decide_one`; it stays required there, but the src-shaped entry defaults it to
`P13_PHASE11_REFERENCE_LAMBDA` — the value the cut was measured at. Recorded because Phase 11 closed
**negative** on inferring registration from the summary: this is a quantity a caller measures
externally (fiducial beads), not one this tool recovers, so the default is a stated reference and
not an estimate. It is documented as such in the `# Named limits` section.

**3. TDD gates could not run in the order the task types imply.** Tasks 1 and 2 carry
`tdd="true"`, but the plan's own `<files>` list puts `test_p14_decide.jl` in Task 3, so no failing
test file could precede them. The RED step used was each task's own `<automated>` verify, run
before the file existed; the git log therefore reads `feat` → `feat` → `test`. Identical to what
14-06 recorded, and it is the plan's task decomposition rather than a skipped gate. The
falsification table above is the stronger substitute: three deliberate mutations, each verified to
land, two of them red and the third leading to a new assertion.

## Surprises worth recording

**The fixture needed no pinned magic numbers, and that turned out to matter.** The discriminating
property of testset 4 — a level at which abstain-first and sort-then-abstain give different answers
— is *searched for inside the test* over `0.01:0.01:0.99` and asserted to exist (11 such levels on
this fixture), rather than transcribed. The same for the operating point that splits the batch (the
sixth-largest Mahalanobis score of the drawn items). Since the Philox stream is counter-based the
draws are identical on every machine, so the test is fully deterministic without a single number
that could go stale when something upstream moves.

## Known Stubs

None. Every function in `spike/p14/decide.jl` is fully implemented and exercised by the test file.
The `posterior` field of the composed Phase-13 result defaults to an **empty** matrix, which is not
a stub but a documented fact: the three-way evidence path produces no posterior draws, and an empty
matrix means "no draws were produced", never "zero draws were interesting".

## Self-Check: PASSED

- `spike/p14/decide.jl` — FOUND (1049 lines)
- `spike/test/test_p14_decide.jl` — FOUND (445 lines)
- `2034573` — FOUND
- `bc3197b` — FOUND
- `c7a5de7` — FOUND
