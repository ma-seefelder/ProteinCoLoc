---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 11
status: complete
subsystem: amortized-inference
tags: [datagen, training, stratification, frozen-basis, three-way-evidence]
requires:
  - "spike/p13/preconditions.jl (13-09) — the Phase-11 binding, frozen zt, lambda encoder"
  - "spike/p13/labels.jl (13-03) — three_way_label, target_matrix, head_log_odds, class_masses"
  - "spike/p13/net.jl (13-06) — build_three_way_net, masked_two_head_bce, train_three_way"
  - "spike/p13/consts.jl (13-01/13-10) — the frozen pre-registration incl. MEASURED P13_TAU = 0.15"
  - "spike/npe/p11_research_npe.jld2 (Phase 11) — the bound research basis"
provides:
  - "spike/p13/datagen.jl — p13_generate_pool / stratify_by_class / assemble_conditioned_pairs / check_frozen_zt"
  - "spike/p13/three_way_net.jld2 — the trained three-way evidence net + measured per-head log-odds + full provenance"
  - "spike/test/test_p13_datagen.jl — the D-02/D-03/D-07/D-11 gate"
affects:
  - "13-12 (lambda response), 13-13, 13-14, 13-16 — all read the trained net and its head_log_odds"
tech-stack:
  added: []
  patterns:
    - "Content-hash-named sharded pool with atomic .tmp -> reopen-integrity -> mv(force) commit"
    - "Resume-by-skip via a Phase-13-local shard-done predicate (cache.jl's _loads_ok would reject every finished shard)"
    - "Two-layer RNG discipline: Philox per global index in generation, Random.seed! immediately before train"
key-files:
  created:
    - spike/p13/train_three_way.jl
    - spike/p13/three_way_net.jld2
    - spike/p13/three_way_train/loss_per_epoch.csv
    - spike/test/test_p13_datagen.jl
    - .planning/phases/13-three-hypothesis-amortized-bayes-factor/deferred-items.md
  modified:
    - spike/test/runtests.jl
    - .gitignore
decisions:
  - "Task order was inverted (3 before 2) to de-risk ~53 KB of inherited uncommitted work before a multi-hour compute — a verification-ordering choice, not a scope change"
  - "Optimiser/train-state BSONs gitignored; the loss history CSVs stay tracked as the run's evidence"
metrics:
  duration: ~115 min wall clock (of which ~20 min pool generation, 2.37 min training)
  completed: 2026-07-29
---

# Phase 13 Plan 11: Labelled Pool + Three-Way Evidence Net Training Summary

The three-way evidence net is trained, once, on the Phase-11 frozen basis from a 48,000-item
class-stratified conditioned pool — realized class frequencies are **exactly** equal thirds, the
standardizer is provably inherited rather than re-fit, and both per-head corrections are measured
on the training subset and persisted with the artifact. The run **overfits early** (best validation
risk at epoch 4 of 45), which is recorded as an honest finding; the recipe was not tuned.

## Execution context: this was a RESUMED plan

This plan was executed by a **third** agent. State inherited from predecessors, verified on disk
before any work began:

| Task | Inherited state | What this run did |
|------|-----------------|-------------------|
| Task 1 `spike/p13/datagen.jl` | **Already committed** as `4839c66` (1117 lines) | Not redone. Its behaviour is exercised by Task 3's gate, which now passes 706/706. |
| Task 2 `spike/p13/train_three_way.jl` | Present but **untracked** (27,417 bytes), written by a predecessor that died on a session limit | Reviewed against the plan's acceptance criteria, found compliant, **not rewritten**. Ran it; committed it with the artifact it produced. |
| Task 3 `spike/test/test_p13_datagen.jl` | Present but **untracked** (25,639 bytes); `runtests.jl` modified +6/-0 | Reviewed, found compliant, **not rewritten**. Ran it, then committed. |
| Training pool | 5 of 12 shards on disk from the dead predecessor | Resumed by skip; 7 further shards generated. |
| `spike/p13/three_way_net.jld2` | **Absent** — the training run had never happened | Produced by this run. |

**Task order was deliberately inverted: Task 3 was executed and committed before Task 2.** ~53 KB
of unreviewed, uncommitted predecessor work was at risk while a ~70-minute compute ran; Task 3's
gate asserts source-level properties of *both* files and needs no trained net, so running and
committing it first converted that risk into a commit before the long spend started. This is a
verification-ordering choice. No task content, scope or acceptance criterion changed.

## What was built

**The pool** (`spike/p13/datagen.jl`, Task 1, pre-committed) generated 48,000 labelled items in 12
shards of 4,000, each standardized with the **inherited** Phase-11 `zt`, each carrying **one**
lambda appended once after `pair_encode`, and each labelled from simulator ground truth alone via
the D-05 three-way cut at the measured `P13_TAU = 0.15`.

**The runner** (`spike/p13/train_three_way.jl`, Task 2) calls `p13_require_phase11()` first, echoes
every recipe constant by name, measures the per-head corrections on the training subset only,
seeds the global RNG immediately before `train_three_way`, and persists **before** printing any
verdict.

**The gate** (`spike/test/test_p13_datagen.jl`, Task 3) is nine testsets wired into `runtests.jl`
additively (+6/-0).

## Measured results

### Realized class frequencies vs the pre-registered target (D-07-i)

| Class | Realized | Target | Delta |
|-------|----------|--------|-------|
| exclusion | 0.333333 | 1//3 | **0.0** |
| random | 0.333333 | 1//3 | **0.0** |
| coloc | 0.333333 | 1//3 | **0.0** |

Exact, not approximate: whole-class accept/reject to an integer quota (16,000 each) makes the
counts exact by construction. Acceptance overhead **1.545×** (48,000 accepted from 74,149 draws),
consistent with the class masses under pi.

### The frozen-zt measurement (Pitfall 5) — **PASS**

```
verdict           = inherited
max|per-row mean| = 0.057353463   (diagnostic tol 0.05)
max|per-row sd-1| = 0.035997152   (diagnostic tol 0.05)
over 96,000 standardized columns x 64 continuous rows
```

The Pitfall-5 warning sign is a **reported number**, not a hope. The discriminating evidence is not
the tolerance comparison but the *scale* of the residual: the gate's negative control fits a
transform on the pool itself and confirms a genuine re-fit lands at machine zero (< 1e-8), while
the inherited transform's residual is > 1e-4 — more than 1e3× larger. The hard half (zt object
identity against the memoized single load) passed as a precondition of getting a verdict at all.

### The two per-head corrections (D-07), measured on the TRAINING subset only

| Head | Measured (training subset) | pi-level reference | Ratio |
|------|---------------------------|--------------------|-------|
| coloc | **+0.0029420438990930215** | −0.7764632135784509 | ~264× smaller |
| exclusion | **+0.0023543271616353607** | −0.4089764597543872 | ~174× smaller |

pi-level class masses measured from 200,000 prior class draws (no simulation): exclusion 0.31272 /
random 0.47073 / coloc 0.21655 — matching 13-RESEARCH E1's 0.3119 / 0.4713 / 0.2168 to ~3 decimals.

The training-subset pair is near zero **by construction**: `P13_TARGET_CLASS_FREQ` is equal thirds,
so stratification is what removes the prior imbalance. It is still **measured, not assumed** — the
contiguous split does not divide a stratified pool into exactly equal thirds (train E/R/C =
13608/13576/13616, val 2392/2424/2384). This is also the concrete reason Pitfall 1 forbids copying
the shipped binary net's correction (−0.0102) by analogy: the pi-level terms here are 40–75× larger.

### The training run — ONE run, recipe unchanged

| Quantity | Value |
|----------|-------|
| Input width | `three_way_input_dim(8, 1)` = **321** (derived, never typed) |
| Split | 40,800 train / 7,200 validation (contiguous at `P13_VAL_FRAC = 0.15`) |
| Initial validation risk | 0.7684933 |
| **Best validation risk** | **0.13625102 at epoch 4** |
| Final training risk | 0.036551084 |
| Final validation risk | 0.38632473 |
| Epochs recorded | 45 of a maximum 300 (patience 40, early stopping fired) |
| Wall clock (training only) | 2.37 min (132.52 s), CPU-only, 32 threads |

**HONEST FINDING — the run overfits early and the recipe was NOT tuned.** Validation risk reached
its minimum at **epoch 4** and then rose monotonically to 0.386 while training risk fell to 0.037 —
a ~10× train/validation gap. Early stopping fired at epoch 45.

The verbatim recipe (`P13_TRAIN_N`, `P13_EPOCHS`, `P13_BATCHSIZE`, `P13_LR`, `P13_WEIGHT_DECAY`,
`P13_VAL_FRAC`, `P13_STOPPING_EPOCHS`) was copied from `src/amortized/train_ratio.jl:155-163` for
the shipped *binary* net; D-10's attribution argument requires it to stay untouched, so any
difference from the binary net is attributable to the head alone. Re-tuning even one value would
convert that comparison into a coincidence and would spend `P13_ITERATION_ALLOWANCE = 1`, which
`P13_ITERATION_TRIGGER` reserves for the stratification switch. **The allowance remains UNSPENT and
`consts.jl` is byte-unchanged.** Per the plan's own instruction ("If it fails to converge, record
the finding — do not tune the recipe"), this is recorded, not repaired.

**What is persisted is the epoch-4 best checkpoint, not the overfit final net.** Verified in the
dependency source: NeuralEstimators' `train` returns `_trainstate_to_device(trainstate_best,
cpu_device())` (`train.jl:643`), i.e. the best-validation-risk state moved to CPU. The overfitting
therefore affects how much of the 300-epoch budget was usable, not which weights shipped.

**Carry-forward for 13-12/13/14/16:** whether a best-at-epoch-4 checkpoint is *sufficient* for the
downstream log-BF and calibration reads is an open question this plan does not answer, because
answering it by retraining is exactly what the iteration allowance forbids without the pre-declared
trigger. It should be judged on the downstream gates.

### Artifact provenance (`spike/p13/three_way_net.jld2`, 953,027 bytes)

Verified present by direct read: `head_log_odds`, `tau = 0.15`, `cut_variant = tau_contrast`,
`consts_sha` (sha256 `100e97a37bb470bb…`), plus in `meta`: `consts_git_blob_sha`
(`5a4ea2223ce11bcae2f974c643ec1e523458e8ac`), `phase11_net` + `phase11_sha`
(`659cda09112c4f20…`), `phase11_lambda_range = (0.25, 3.0)`, `realized_class_freq`,
`target_class_freq`, `stratification`, `head_log_odds_scope`, `zt_provenance`,
`frozen_zt_verdict/max_mean/max_sd_dev`, `pi_class_masses`, `pi_head_log_odds`, `pool_dir`,
`master_seed`, `salt`, `datagen_counter`, every recipe constant, `train_minutes`, `savepath` and a
UTC `generated` stamp. No `.tmp` sibling remains.

## The training pool: where it lives and what losing it costs

`spike/data/cache/p13/7c65a1a90caeb23b4511e104fb553f4c30245f66378845686f791d03520a29a1/`
— 12 shards + `meta.jld2`, **50 MB**.

**It is GITIGNORED** (`.gitignore:407`, `spike/data/cache/*`). It is therefore **not protected by
version control and has no backup**. This project has already destroyed a comparable artifact once:
`git worktree remove --force` deleted Phase 11's 54 MB pool. This plan was executed on the **main
working tree with no worktree** specifically so that could not recur.

**Regenerating it costs ~2 hours of CPU** (measured ~10 min per 4,000-item shard, 12 shards).
Because generation is counter-based Philox keyed per global index, regeneration is
**byte-identical** — so a loss is a **compute cost, not a re-seed** and would not invalidate any
reported number. It must not be described as "available for reuse at no further cost".

The resume-by-skip design was exercised for real twice in this run and worked both times: the
5 predecessor shards were skipped (their Jul-28 mtimes verified unchanged, and the run reported
`[resume] 10 of 12 shards already complete`), and when the process was killed mid-run at 10 shards
it resumed with no torn `.tmp` and no lost shard.

## Verification performed

| Check | Result |
|-------|--------|
| `julia --project=spike spike/test/test_p13_datagen.jl` | **exit 0, 706 pass / 706 total**, zero skips (Phase-11 net present, so every `@test_skip` arm ran for real) |
| `julia --project=. -e 'using Pkg; Pkg.test()'` (root suite) | **exit 0, green** |
| `julia --project=spike spike/test/runtests.jl` (full spike suite) | **exit 1 — pre-existing, NOT a regression.** See below. |
| Task-2 artifact verify (corrected form) | prints `net artifact ok` |
| Rebuilt net output shape | `(2, 3)` |
| No `.tmp` sibling of the artifact | confirmed |
| `grep -c 'include(joinpath(@__DIR__, "test_p13' spike/test/runtests.jl` | **11** (as required) |
| `grep -c 'test_p13_datagen.jl' spike/test/runtests.jl` | **1** (as required) |
| `git diff --numstat -- spike/test/runtests.jl` | **6 insertions, 0 deletions** |
| `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | **exit 0 — byte-unchanged** |
| Recipe literals absent from trainer source | `48_000` / `2.5e-4` → 0 occurrences |
| `Random.seed!` in `datagen.jl` | 0 (Philox per index only) |
| `ZScoreTransform` constructed in `datagen.jl` | 0 |
| Forbidden basis (`trained_npe.jld2`, `grid_8`, `artifacts/`) reachable | 0 |
| `use_gpu = true` in trainer | 0; CUDA never loaded |

**On the full-suite exit 1.** The suite aborts at `spike/test/runtests.jl:169` →
`test_npe.jl:230`, the Phase-4 `SPEEDUP_GATE` (median speedup **68.35** against the pre-registered
bar 100.0). This is the blocker STATE.md already records from the Phase-11 wave-3 post-merge gate
(then 92.50× and 83.97×). It is **not** caused by this plan and `test_npe.jl` is not a file this
plan owns. Two details are new and are logged in
`13-three-hypothesis-amortized-bayes-factor/deferred-items.md`: the number has drifted further
down, and the abort is ~44 lines **earlier** than STATE.md states — so the masked region is the
entire Phase-13 include block, not just what follows `test_p13_correction.jl` at line 213.
Per-file runs are the available gate signal, which is the practice STATE.md already prescribes for
the rest of Phase 13.

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 1 — Bug] Plan 13-11's own Task-2 `<verify>` command is defective as written**
- **Found during:** Task 2 verification.
- **Issue:** the plan's command is
  `@assert haskey(h, "head_log_odds") || hasproperty(h, :head_log_odds)`. `load_three_way` returns
  a `NamedTuple`, and `haskey(::NamedTuple, ::String)` **throws a `MethodError`** rather than
  returning `false`, so `||` never short-circuits to the `hasproperty` branch. The command fails
  on a perfectly valid artifact.
- **Fix:** operands swapped (`hasproperty(...) || haskey(...)`), which then passes and prints
  `net artifact ok`. The artifact was never in question — every claimed provenance key was
  independently confirmed present by direct read.
- **Files modified:** none (the defect is in the PLAN text, which this executor does not own).
- **Also logged in:** `deferred-items.md` (D-13-B).

**2. [Rule 3 — Blocking] `.gitignore` rule for the training side effects**
- **Found during:** Task 2 commit.
- **Issue:** NeuralEstimators writes `best_optimizer.bson`, `best_trainstate.bson`,
  `final_optimizer.bson` and `final_trainstate.bson` (~9.4 MB total) into the run's `savepath`
  alongside the loss history. Leaving them untracked violates the "never leave generated files
  untracked" rule; committing them adds 9.4 MB of machine state redundant with the persisted net.
- **Fix:** a scoped, commented `.gitignore` rule following the house convention — the BSONs are
  ignored, and `loss_per_epoch.csv` + `train_time.csv` are explicitly negated back in, because the
  loss history is precisely why the plan requires `savepath` to be named rather than left in
  `tempdir()`.
- **Files modified:** `.gitignore` (additive, 8 lines).
- **Commit:** `41a24f4`.

### Process deviations (not code)

**3. Task order inverted — Task 3 executed and committed before Task 2.** Rationale in
"Execution context" above. A verification-ordering choice, not a scope change; no task content or
acceptance criterion was altered.

**4. The first training process was killed by the harness mid-run** (after 10 of 12 shards, with
the block-buffered stdout for that attempt lost). Recovered by relaunching detached; the atomic
shard design meant **zero** shards were lost and none were regenerated. Recorded because it is the
second time in this project that a long compute has been interrupted, and the resume design is
what made it a non-event both times.

## Assumption Drift (advisory)

**A. The training run was budgeted as the long pole; it was not.**
- **Planned:** the plan and the dispatch brief treat the training run as the phase's "one large
  compute spend", following Phase 11 where training took 16.08 min after a 56.37 min datagen.
- **Actual:** training took **2.37 minutes**. Pool generation took ~120 minutes across two
  sessions. The spend is ~98% data generation.
- **Why:** the three-way net is a small MLP trunk plus two scalar heads, not a normalising flow —
  epochs cost ~1.8 s each after the first. The expensive part is simulating and summarizing 48,000
  image pairs under the F5 image-size mixture.
- **Materiality:** it changes how later Phase-13 plans should be scheduled and, more importantly,
  it means the *marginal* cost of a second training run is minutes — which makes the
  `P13_ITERATION_ALLOWANCE` discipline a matter of pre-registration integrity rather than of
  compute budget. That distinction should be stated explicitly wherever the allowance is discussed,
  so nobody argues for a retrain on the grounds that it is cheap. **It is cheap; that is not the
  reason it is forbidden.**

**B. The full spike suite was expected to reach the Phase-13 block.**
- **Planned:** the plan's acceptance criteria assume `runtests.jl` runs to completion and exits 0,
  and STATE.md attributes the known abort to `test_p13_correction.jl` at line 213.
- **Actual:** the suite aborts at line 169 in Phase-4 `test_npe.jl`, so no `test_p13_*` file runs
  in a full-suite invocation at all.
- **Why:** a second, earlier pre-existing red (`SPEEDUP_GATE`) that STATE.md records as a blocker
  but does not connect to the masking effect.
- **Materiality:** any future plan that writes "the spike suite is green" as an acceptance
  criterion will be unsatisfiable until one of the two reds is resolved by a user decision.

## Requirements satisfied

D-02 (frozen inherited zt), D-03 (one lambda appended once, width derived), D-07 (whole-class
stratification + measured per-head correction), D-11 (four-row target matrix), D-10 (verbatim
recipe), D-01 (reproducible from the reserved dev stream, global RNG seeded before train), D-04
(iteration allowance untouched).

## Commits

| Commit | Task | Description |
|--------|------|-------------|
| `4839c66` | Task 1 | *(pre-existing — landed by a predecessor agent, not this run)* labelled class-stratified conditioned datagen on the Phase-11 basis |
| `133558d` | Task 3 | assert the stratification design, the inherited zt and the derived pair width |
| `41a24f4` | Task 2 | train the three-way evidence net once on the Phase-11 basis |

Commits from a parallel agent working on Phase 12 are interleaved in the branch history between
`133558d` and `41a24f4`. All commits here used explicit paths only; no race occurred on the index
and nothing was rebased or amended.

## Self-Check: PASSED

- `spike/p13/datagen.jl` — FOUND
- `spike/p13/train_three_way.jl` — FOUND
- `spike/p13/three_way_net.jld2` — FOUND (953,027 bytes)
- `spike/p13/three_way_train/loss_per_epoch.csv` — FOUND
- `spike/test/test_p13_datagen.jl` — FOUND
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/deferred-items.md` — FOUND
- commits `4839c66`, `133558d`, `41a24f4` — all FOUND in history
