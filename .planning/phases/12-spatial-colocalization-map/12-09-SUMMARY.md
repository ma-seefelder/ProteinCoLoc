---
phase: 12-spatial-colocalization-map
plan: 09
status: complete
subsystem: training-pool generation (field-aware datagen)
tags: [datagen, content-hash, resume-by-skip, theta-layout, single-stack, wall-clock-blocker, mutation-tested]
requires: ["12-01", "12-03", "12-04", "12-05", "12-07", "12-08"]
provides:
  - "spike/data/p12_generate.jl — p12_datagen_rng, p12_theta_column (the phase's ONLY theta assembler), sample_p12_imsize, generate_p12_sample/_block, p12_shard_loads_ok/_done, write_p12_shard, p12_generating_config, write_p12_meta, p12_pool_dir, p12_pool_complete, generate_p12_pool, load_p12_pool, realized_imsize_counts, realized_r1_quantiles, P12_CACHE_ROOT, P12_POOL_SCHEMA, P12_SUMMARY_MIN_DIM"
  - "spike/test/test_p12_datagen.jl — 8 testsets / 65 assertions replacing the 12-01 scaffold"
affects: ["12-11", "12-13", "12-14", "12-15", "12-17"]
tech-stack:
  added: []
  patterns:
    - "a keyword a SIBLING plan depends on must be visible in the signature that declares it — the r1 pin reaches both the draw and the content hash, and both halves are asserted"
    - "a phase-local meta writer is what lets the REALIZED distributions ride in the manifest without entering the hash the manifest names"
    - "a defensive runtime assertion can violate a stricter STRUCTURAL ban — 12-05 forbids the p11 literal outright, which the guard needed to spell"
key-files:
  created:
    - spike/data/p12_generate.jl
    - .planning/phases/12-spatial-colocalization-map/12-09-SUMMARY.md
  modified:
    - spike/test/test_p12_datagen.jl
decisions:
  - "NO LARGE POOL WAS GENERATED, AND THE REASON IS THAT ITS KEYWORDS ARE NOT KNOWABLE TODAY. Every keyword that enters the content hash for the 50 000-pair pool is owned by a later plan: 12-17 generates it with `arm = P12_CHOSEN_PRIOR`, a TIER-2 constant 12-15 has not yet appended. 12-11 generates its own five 5 000-sample rung pools, 12-15 its own 10 000-sample per-arm pools. 12-14:375 states the scope verbatim: ' 12-09 only builds 64 samples into a mktempdir()'. Guessing an arm would have burned ~73 min into a directory nothing reads."
  - "MEASURED: 87-90 ms/sample on the real F5 mixture at 32 threads. P12_N_PAIRS = 50 000 projects to ~73 min against the append-only 150-min ceiling; P12_MINISPIKE_N to ~15 min; one 12-11 rung (5 000) to ~7.3 min. The research's 65-80 min/50k estimate is confirmed. The ceiling is not at risk for any pre-registered size."
  - "The draw is consumed BEFORE the image size, so the header's 'r1 DRAWN **FIRST**' is literally true of the key stream. The plan's prose listed imsize first; the analog (p11_generate.jl:226-238) draws theta first, and the header's ordering claim is the stronger constraint."
  - "ADDED p12_pool_complete: p12_pool_dir CREATES the directory it resolves (cache.jl:162), so `isdir(p12_pool_dir(...))` is TRUE for a pool that has never been generated. Three later plans currently plan an isdir presence check that cannot fail."
  - "ADDED an 8th testset exercising the wall-clock BLOCKER, because a guard verified only by grepping its own source for the word BLOCKER has not been shown to fire."
metrics:
  duration: ~75 min
  completed: 2026-07-29
  tasks: 2
  commits: 3
  files: 2
---

# Phase 12 Plan 09: The Field-Aware Training Pool — Summary

`spike/data/p12_generate.jl` turns 12-07's prior and 12-08's field-accepting simulator into the
data every later Phase-12 plan consumes: r₁ drawn first, the field conditional on it, the RAW
128-row summary and the **drawn 8×8 lattice** stored side by side, sharded and content-hash
addressed under a cache root that cannot collide with Phase 11's.

---

## 1. The thing the orchestrator asked for that I did not do, and why

The dispatch's success criteria include *"Pool generated at the pre-registered size, in the PRIMARY
checkout, wall-clock MEASURED and reported against the 150-min ceiling."* **I generated no large
pool.** This is reported rather than silently done or silently skipped.

**The plan of record contains no such task.** Its two tasks build the file and its test; its only
generation is 64 + 32 samples. And three sibling plans each own a pool of their own, with keywords
that **enter the content hash**:

| Plan | Pool it generates itself | Keywords in the hash |
|---|---|---|
| 12-11:131 | five pools of `P12_STAGE1_N ÷ 5` = 5 000 | `arm = :car`, **`r1 = rung`** (one per rung), imsize pinned to `P12_STAGE1_IMSIZE` |
| 12-15:147 | `P12_MINISPIKE_N` = 10 000, per arm | `arm ∈ (:car, :gp, :none)`, imsize pinned to `P12_MINISPIKE_IMSIZE` |
| 12-14:377 | `P12_MINISPIKE_N`, `arm = :none`, **descope branch only** | F5 mixture |
| 12-17:137 | `P12_N_PAIRS` = 50 000 | **`arm = P12_CHOSEN_PRIOR`** |

`P12_CHOSEN_PRIOR` is a **Tier-2 constant that does not exist yet** — 12-15 appends it under the
reserved `:P12_CHOSEN_PRIOR` sentinel after measuring CAR against GP. So the arm of the 50 000-pair
pool is *unknowable today*. Generating `arm = :car` would be a guess; if 12-15 selects `:gp` the
~73 minutes lands in a directory nothing will ever read. The same is true, one keyword down, for
12-11's rungs: its `imsize_tag` is not specified in its plan text, and a mismatched tag hashes to a
different directory.

12-14:375 states the settled scope in its own words, in a plan that passed the same review gates:

> *"12-09 only builds 64 samples into a `mktempdir()`, 12-11 builds five 5,000-sample `arm = :car`
> pools whose content hash differs, and 12-15 is on the skip list."*

**What I did instead** — because the measurement is the deliverable — is measure the per-sample cost
on the real F5 mixture and report the projection for every pre-registered size (§4). That is the
fact the 150-minute ceiling exists to establish, and it now exists without spending an hour on a
guessed hash.

**If the orchestrator wants a pool pre-built anyway,** the only one that is fully determined today
is 12-11's ladder, and it still needs `imsize_tag` pinned by a ruling before it is safe to spend.

---

## 2. What was built

### `spike/data/p12_generate.jl`

**The generative order is enforced, not just documented.** The header states r₁ → field → nuisances
→ imsize, and the code cannot contradict it: this file never re-implements a draw, it calls
`sample_p12_prior`, which is the single place the order lives (`p12_prior.jl:184-207`).

**`p12_theta_column` is the phase's only θ assembler.** It builds
`[c₀ ; c₁..c₆₃ ; 7 nuisances ; r₁]` from the **Gaussian** field's DCT coefficients in the frozen
`p12_dct_order`, and guards itself with in-function `@assert`s read off `p12_theta_index` — the
layout, not arithmetic — so a layout change cannot pass silently here.

**Both fields are stored, deliberately.** `z_field` is the atom-free SBC target (D-07 ground truth);
`rho_field` is the deliverable's units for 12-19; `mu_field` makes the copula step checkable after
the fact. Fields are `G×G×N` arrays, not vectors of matrices.

**The pool is SINGLE-STACK.** No paired control draw. The header records why (Δρ is a read-time
construction from two independent passes; a paired pool would take 65-80 min/50k to ~130-160 min
against an append-only 150-min ceiling), so the decision is legible at the point where someone
would be tempted to re-add it.

**`arm`, the r₁ prior bounds AND the r₁ pin all enter the content hash.** The `r1_pin` field is the
one that keeps 12-11's five rungs in five directories. Its absence was mutation-proved to be
undetectable by anything else (§5).

**Storage is raw.** Comment-stripped, the file contains zero occurrences of `fit(ZScoreTransform`,
`standardize`, `reshape_summary`, `augment_mask!`, and exactly one `simulate_pair(` and one
`patch_summary(`.

### `spike/test/test_p12_datagen.jl`

The 12-01 scaffold is replaced. **8 testsets, 65 assertions, 0 failures, 0 errors, 5.1 s** (25 s
including package load — the plan's 90 s budget). Every pool goes into a `mktempdir()`; the real
cache tree's directory listing is identical before and after a run.

---

## 3. Verification actually run (not inferred)

| Check | Command | Observed |
|---|---|---|
| Plan Task-1 verify, verbatim | `julia --project=spike -t auto -e '…generate_p12_pool(64;…)…'` | exits 0, prints `(72, 64) datagen-ok`; dir under `spike/data/cache/p12/` |
| Plan Task-2 verify | `julia --project=spike -e 'using Test; include("spike/test/test_p12_datagen.jl")'` | **65 pass / 65 total**, 8 testsets |
| Phase-12 aggregator | `include("spike/test/test_p12_suite.jl")` | all ten files ran; datagen green **inside** the suite (65/65); no failures anywhere |
| 12-05 decoupling (LIVE) | `include("spike/test/test_p12_decoupling.jl")` | **26 pass / 26** — see §5, it caught a real defect first |
| Resume-by-skip | delete last shard, regenerate | `[resume] 1 of 2 shards already complete — skipping them`; shard 1 **mtime unchanged**; reloaded pool `==` the original |
| Thread independence | `parallel = false` vs `true`, separate roots | `theta`, `summary_min`, `z_field` all `==` |
| Wall-clock BLOCKER | `wallclock_ceiling_min = 0.0` | throws `ErrorException`, message contains the literal `BLOCKER` |
| Scaffold marker gone | `grep -c P12_PENDING_SCAFFOLD` | `0` |
| Frozen surfaces | `git diff --quiet HEAD -- src corpus spike/Project.toml spike/Manifest.toml spike/simulator spike/npe` | **CLEAN** |
| `.planning` files I must not touch | `git status --porcelain -- .planning/STATE.md .planning/ROADMAP.md` | empty |

### The Phase-11 pool, before and after

Byte-identical. Sizes **and** mtimes, all six files:

```
BEFORE and AFTER (identical):
  meta.jld2         11 764 bytes   mtime 1785165536
  shard_0001.jld2   11 206 103     mtime 1785163029
  shard_0002.jld2   11 206 103     mtime 1785163648
  shard_0003.jld2   11 206 103     mtime 1785164275
  shard_0004.jld2   11 206 103     mtime 1785164906
  shard_0005.jld2   11 206 103     mtime 1785165536
  du -sh spike/data/cache/p11 -> 54M
git status --porcelain -- spike/data/cache/p11 -> empty
```

---

## 4. Measurements on the record

Measured this session, 32 threads, CPU-only, on the **real F5 mixture**
(`((512,512), (1024,1024), (1376,1028), (2048,2048))` at `(0.4, 0.25, 0.25, 0.1)`):

| Quantity | Measured |
|---|---|
| Prior draw alone (CAR bisection + 64×64 Cholesky), serial | **2.18 ms/sample** |
| Full sample on the F5 mixture, n = 64 | 5.8 s → **89.9 ms/sample** |
| Full sample on the F5 mixture, n = 128 | 11.1 s → **87.0 ms/sample** |

Projected against the append-only `P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150`:

| Pre-registered size | Projection | Margin |
|---|---|---|
| one 12-11 rung, 5 000 | **~7.3 min** | 20× |
| `P12_MINISPIKE_N` = 10 000 | **~14.5 min** | 10× |
| `P12_N_PAIRS` = 50 000 | **~72.5 min** | ~2× |

This **confirms the research's 65-80 min/50k estimate** and says the ceiling is not at risk for any
pre-registered size. Realized image-size counts at n = 128 (49 / 35 / 28 / 16 = 0.383 / 0.273 /
0.219 / 0.125) track the F5 weights.

### The regeneration cost, stated honestly

`spike/data/cache/p12/` is **gitignored and unbacked**, exactly like the Phase-11 pool that this
milestone already lost once to a worktree cleanup and paid ~56 minutes to rebuild. Every pool this
file writes is a **pure function of its global index** through the counter-based Philox stream
`P12_DEV_SEED ⊻ P12_DATAGEN_SALT`, so a destroyed pool **regenerates byte-identically** — verified
here, both across independent cache roots and across a parallel/serial split.

**Regenerating is byte-identical; it is NOT free.** At the measured 87 ms/sample, rebuilding the
50 000-pair pool costs **~73 minutes of CPU**, the mini-spike pool ~15 minutes, and 12-11's five
rungs ~37 minutes together. Nothing in this phase may be described as "available for reuse at no
further compute cost" — that sentence is what made Phase 11's loss expensive.

The 372 KB currently under `spike/data/cache/p12/` is verification scratch (64- and 32-sample
pools plus four empty hash-probe directories), not a deliverable.

---

## 5. Deviations from plan

### `[Rule 1 — Bug] The datagen source violated 12-05's live structural ban on the `p11` literal`

- **Found during:** post-Task-2 run of `test_p12_decoupling.jl`
- **Issue:** `test_p12_decoupling.jl:337` requires that the Phase-11 cache-root literal appear
  **nowhere** in `p12_generate.jl`'s comment-stripped source, so that no code path can construct
  that path at all. My load-time guard `@assert !occursin(joinpath("cache", "p11"), P12_CACHE_ROOT)`
  had to *spell* the banned literal in order to check for it, and failed the live test on its own
  defensive assertion (25 pass / 1 fail).
- **Why 12-05 is right and the guard was wrong:** a structural absence is a stronger guarantee than
  a runtime check that a constant is not something. The ban makes the p11 path unconstructible; the
  assertion only made it un-equal.
- **Fix:** removed the assertion, kept `@assert basename(P12_CACHE_ROOT) == "p12"`, and recorded the
  reasoning at the site. The negative half is asserted from the *test* files, where naming another
  phase's root is what they are for (my testset 3 and 12-05 itself).
- **Verified:** `test_p12_decoupling.jl` → **26 pass / 26**.
- **Commit:** `8880c29`

### `[Rule 2 — Missing critical functionality] `p12_pool_complete`, because `isdir` cannot fail`

- **Found during:** Task 1, while probing arm/r₁ hash separation — four *empty* directories appeared.
- **Issue:** `p12_pool_dir` resolves through `open_or_invalidate`, which **creates the directory it
  resolves** (`cache.jl:162`). So `isdir(p12_pool_dir(n; arm = arm))` is `true` the first time it is
  ever asked, for a pool that has never been generated. Demonstrated live:

  ```
  isdir right after resolve (the vacuous check) = true
  p12_pool_complete before generating           = false
  p12_pool_complete after generating            = true
  ```

  **This is 12-06's defect class exactly** (*"an `@assert isdir(dir)` that cannot fire, because the
  resolver creates the directory it checks"*), and it is aimed at three later plans: **12-14:203**
  plans to "assert `isdir`"; **12-15:147** and **12-17:137** plan to generate "only if absent".
- **Fix:** added `p12_pool_complete(dir, n; shard_size)` — true iff every shard passes the reopen
  integrity check *and* the manifest exists. Documented at the site as *the* presence test.
- **Commit:** `13dc121`

### `[Deliberate addition] An 8th testset exercising the wall-clock BLOCKER`

The plan specifies seven testsets and verifies the wall-clock guard by grepping the source for the
word `BLOCKER`. A guard verified by grepping for its own message has not been shown to **fire**, so
testset 8 drives `wallclock_ceiling_min = 0.0` and requires an `ErrorException`. Mirrors
`test_p11_generate.jl:179-184`. All seven plan-specified testsets are present and unmodified.

### Falsification proofs (mutate → red → revert byte-exactly)

Per §6 of the briefing, every load-bearing assertion was proved able to fail. All three reverts were
confirmed by SHA-256 against a pre-mutation copy.

| Mutation | Assertion under test | Result |
|---|---|---|
| Drop `r1_pin` from `p12_generating_config` | the r₁-pin hash separation | `AssertionError("the r1 PIN must enter the content hash")` — **red** |
| Append `r1` before `chromatic_eps` in `p12_theta_column` (with its in-function guard silenced, so the *testset* is what is exercised) | testset 1 | **3 failures** at lines 98, 99, 106 — the `:r1`, adjacency and by-name-loop clauses |
| Introduce per-row z-scoring into `generate_p12_block` | testset 4 | **1 failure** at line 182 — the non-zero-row-mean falsifier. The shape and mask clauses still passed, confirming the falsifier is the load-bearing one |

---

## 6. Advisories for later plans (nothing gated depends on these)

1. **`imsize_tag` is the one keyword a pinned-size caller must not forget.** The image-size *sampler*
   is a closure and cannot be hashed, so `imsize_tag` is the only thing separating an F5-mixture pool
   from a pool pinned to one size. 12-11, 12-14 and 12-15 all pin `imsize_sampler` in their plan text
   and **none of them names an `imsize_tag`**. Two consequences:
   - a pinned pool generated without a tag hashes to the *F5-mixture* directory and is
     indistinguishable from one;
   - **12-14's descope branch** (10 000, `arm = :none`, F5) and **12-15's `:none` arm** (10 000,
     `arm = :none`, pinned 512²) would then collide on *the same directory*. They are mutually
     exclusive by routing today, so this is latent rather than live — but it is T-12-34 in a
     different guise. The realized counts are written into `meta.jld2` (F6), so the mistake is at
     least diagnosable after the fact.
2. **`p12_pool_dir` is not a pure query — it creates.** Use `p12_pool_complete`, not `isdir`. (§5.)
3. **The prior draw costs 2.18 ms/sample** in CAR bisection + Cholesky, recomputed per sample even
   when r₁ is pinned. That is ~2.5 % of the 87 ms sample cost on the F5 mixture, but it is ~30 % of
   the cost at 512² alone — relevant to 12-11's five pinned-r₁ rungs, where the *same* Σ is solved
   5 000 times per rung. Not optimised here: the plan mandates going through `sample_p12_prior`, and
   the projection clears the ceiling comfortably either way.

## 7. Assumption Drift (advisory)

- **Draw order vs. image-size order.** *Planned:* the plan's `<action>` lists `imsize` before
  `draw = sample_p12_prior(...)`. *Actual:* the draw is consumed first, the image size second.
  *Why:* the plan's own header specification says r₁ is *"DRAWN **FIRST**"*, and the analog
  (`p11_generate.jl:226-238`) draws θ before imsize. Consuming an image size ahead of r₁ would make
  the header's ordering claim literally false of the key stream. Nothing downstream depends on which
  ordering is used — both are deterministic — but the two produce different (equally valid) pools,
  so the choice is recorded rather than left implicit.

---

## 8. Threat Flags

None. This plan adds no network endpoint, no auth path, no file access outside the gitignored
Phase-12 cache root, and no schema change at a trust boundary. The `<threat_model>` dispositions
were applied: T-12-32 (single assembler + layout-keyed `@assert` + by-name testset + mutation proof),
T-12-11 (distinct cache root, testset 3, before/after p11 listing), T-12-33 (raw storage, zero
`fit(ZScoreTransform` in comment-stripped source, the non-zero-mean falsifier), T-12-34 (`arm` and
`r1_pin` in the config — see the §6 advisory for the residual `imsize_tag` gap), T-12-17
(resume-by-skip + the BLOCKER, both exercised), T-12-35 (realized records in `meta.jld2`), T-12-07
(all three commits pathspec-scoped to my own files, on a branch a Phase-13 executor is
concurrently committing to).

---

## 9. Known Stubs

None.

---

## Self-Check: PASSED

- `spike/data/p12_generate.jl` — FOUND
- `spike/test/test_p12_datagen.jl` — FOUND (scaffold marker count `0`)
- `13dc121` `feat(12-09): field-aware sharded datagen with its own cache root` — FOUND
- `4513a48` `test(12-09): replace the datagen scaffold with the row-layout and ground-truth gate` — FOUND
- `8880c29` `fix(12-09): drop the p11 literal 12-05 bans from the datagen source` — FOUND
- `spike/data/cache/p11` — 54 MB, six files, sizes and mtimes unchanged
- `src`, `spike/Project.toml`, `spike/Manifest.toml`, `spike/simulator`, `spike/npe`, `corpus` — byte-unchanged
- `.planning/STATE.md`, `.planning/ROADMAP.md` — untouched
