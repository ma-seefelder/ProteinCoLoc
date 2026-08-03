---
phase: 14
plan: 08
subsystem: draw-and-null-machinery
status: COMPLETE
tags: [pools, unstratified, exchangeability, ood-null, basis-provenance, d-03, d-06, d-07, wave-4]
requires:
  - "14-01 — spike/p14/consts.jl frozen: P14_DEV_SEED / P14_FIX_SEED / P14_SALT, the reserved counters, and the key-word cross-product `_p14_key_words()` / `_p14_p13_key_words()`"
  - "14-02 — spike/p14/provenance.jl: p14_load_tau (D-07), P14_P13_DIR, the provenance record"
  - "14-04 — spike/p14/posterior.jl: pulled for the class surface and the Phase-13 net handle"
  - "14-05 — spike/p14/fuse.jl: p14_ood_state, p14_fuse (the D-06 chain this file proves end to end)"
  - "spike/p13/datagen.jl, labels.jl, preconditions.jl — read-only: p13_simulate_indices, p13_draw_labelled_item, p13_assert_classes_agree, class_masses, load_pool, load_p13_basis"
  - "spike/validation/ood.jl — read-only: fit_ood_nulls, maha_score, id_threshold, OOD_ID_QUANTILE, OOD_FAMILIES, OOD_GRID_LEVELS"
  - "spike/p13/three_way_net.jld2 + the 51 MB gitignored pool at meta.pool_dir — read-only"
provides:
  - "p14_draw_pool — the UNSTRATIFIED, pi-distributed draw every reported Phase-14 number rides"
  - "p14_draw_classes — the theta-only twin at the same indices, pinned to the pool by p13_assert_classes_agree"
  - "p14_assert_unstratified — the equal-thirds detector, throwing with the CAUSE before the symptom"
  - "p14_datagen_key_word / p14_datagen_rng — the per-item stream, with 14-01's cross-product membership asserted at every draw"
  - "p14_ood_reference — the density null on the Phase-13 net's OWN pool, degrading to available=false rather than substituting a threshold"
  - "p14_misspec_pool — matched in-distribution vs misspecified arms from the pre-registered families"
  - "p14_pool_provenance — the pool block every runner embeds"
  - "p14_gate_dev_substrate — the Phase-13 gate set, labelled development-only as FIELDS"
  - "p14_cache_root — the new Phase-14 cache root; cache/p13 stays read-only"
affects:
  - "every Phase-14 runner (SC1-b, SC1-d, SC3-a..d): they draw through p14_draw_pool and persist p14_assert_unstratified's realized masses"
  - "the SC3-d misspecification arm consumes p14_misspec_pool and p14_ood_reference"
  - "the phase report: `channels_not_wired = (:pp, :noise)` is a NAMED LIMIT, not a footnote"
tech-stack:
  added: []          # zero packages; every draw and null primitive already existed in-repo
  patterns:
    - "the plan's own suggested lambda was the bug: a stream constructor that ignores the item index returns n identical items, and every reproducibility assertion still passes on them"
    - "a cheap theta-only twin of an expensive draw, PINNED to it by the existing prescreen/simulation agreement check rather than trusted"
    - "the balance check runs BEFORE the pi check so a balanced set reports its cause, not its symptom -- and the mutation that removes it is caught by the MESSAGE, not the throw"
    - "a docstring is CODE to a comment-stripping guard: every forbidden name in this module lives in a `#` block"
    - "a failure path returns basis_provenance = :unavailable, because a null that was never fitted has no basis to claim"
    - "a listing digest named `pool_listing_sha256` so it cannot be mistaken for a content hash"
key-files:
  created:
    - spike/p14/pools.jl               # 699 lines
    - spike/test/test_p14_pools.jl     # 362 lines, 102 assertions
  modified: []
decisions: [D-03, D-06, D-07]
validates: ["SC1-b", "SC1-d", "SC3-a", "SC3-d", "Prior re-derivation"]
metrics:
  duration: ~95 min
  completed: 2026-08-03
  tasks_completed: 3
  tasks_total: 3
  per_file_pass_count: 102
  test_wall_clock_s: 74
  mutations_run: 6
---

# Phase 14 Plan 08: The Draw and Null-Fitting Machinery — Summary

**Both silent failure modes are now closed by assertions that were watched to FAIL, and the OOD
reference reproduces Phase 13's own reported in-distribution threshold to the last digit —
167.54463293780395 on 96,000 columns across 12 shards — which is the strongest available evidence
that the construction really was reused rather than re-derived.**

## READ THIS BEFORE QUOTING ANY NUMBER FROM THIS PLAN

Nothing here is a result. This plan builds the machinery five later runners draw through; the only
numbers below are fixture-scale draws on the FIXTURE stream and a null fit on an existing pool. No
reported stream was consumed, and no bar was scored against anything.

## What was built

`spike/p14/pools.jl` (699 lines) and `spike/test/test_p14_pools.jl` (362 lines, 102 assertions,
74 s, CPU-only).

### Pitfall 1 — the calibration set is pi-distributed, and a balanced one fails LOUDLY

`p14_draw_pool` draws through `p13_simulate_indices`, the unstratified path. The Phase-13
accept/reject balancer is named only inside `#` comment blocks, because the guard that forbids it
greps the comment-stripped source — a docstring is code to that grep, and this project has now hit
that trap repeatedly.

`p14_assert_unstratified` runs TWO checks in a deliberate order:

1. the **balance** check (all three masses simultaneously within 3 binomial SE of 1/3) — first,
   because it is the more specific diagnosis;
2. the **pi-agreement** check against `meta.pi_class_masses`.

They are not redundant. At n = 600 the EXCLUSION mass alone sits 1.2 SE from 1/3, so a "some mass is
near a third" test would fire on a perfectly good draw; only the simultaneous condition identifies a
balanced set. And a balanced set would also fail check 2 — but as *"it rode a stream or a cut it was
not supposed to ride"*, which sends the reader hunting the wrong thing. The mutation that disables
check 1 is caught by the **message**, not by the throw.

Realized masses of the reference fixture draw (n = 600, `P14_FIXTURE_COUNTER`, fixture stream):

| class | realized | expected (`meta.pi_class_masses`) | deviation |
|---|---|---|---|
| exclusion | 0.30500 | 0.31272 | 0.41 SE |
| random    | 0.49000 | 0.47073 | 0.95 SE |
| coloc     | 0.20500 | 0.21655 | 0.69 SE |

Distance from equal thirds: coloc 2.2 SE and random 2.7 SE **outside** the 3-SE balance band, which
is what makes the detector non-vacuous at this sample size.

### Pitfall 4 — the null is fitted on the net's OWN basis

`p14_ood_reference()` resolves `pool_dir` from the trained net's `meta.pool_dir` (never guessed),
reads all 12 `shard_*.jld2`, stacks BOTH the `Zs` and `Zc` halves, fits the density null and takes
the pre-registered in-distribution quantile:

| quantity | value |
|---|---|
| `thr` | **167.54463293780395** |
| `n_pool` | 96,000 columns (48,000 pairs, both halves) |
| `n_shards` | 12 |
| `pool_dir` | `spike/p13/../data/cache/p13/7c65a1a90caeb23b4511e104fb553f4c30245f66378845686f791d03520a29a1` |
| pool bytes | 51,972,300 |
| `pool_listing_sha256` | `ec2e3f8b9ecbde25f81818a9e3de19737d40216af8edd81f118bf4d33235e93e` |
| `q` | 0.95 (`OOD_ID_QUANTILE`) |
| `basis_provenance` | `:p13_pool_phase11_zt` |
| score quantiles | q50 40.17 / q90 149.57 / q95 167.54 / q99 199.47 / max 290.69 |

**The pool was found and used on the main working tree.** It is gitignored bulk data and does not
exist inside a worktree; this plan would have produced a silently wrong (or absent) null there.

**Independent cross-check:** `spike/p13/realimage_report.jld2` records
`ood_phase13.threshold = 167.54463293780395`, `n_id = 96000`, `n_shards = 12` — identical in every
digit. Two implementations of one construction, run months apart, agreeing exactly.

### D-06, proved end to end

`p14_ood_reference(pool_dir = <nonexistent>)` returns `available = false`, `thr = nothing`, an
`ood_nulls` **without** a `:thr` key, and a note — it never throws and never substitutes a value.
The test then walks the full chain:

```
p14_ood_reference (missing pool) -> p14_ood_state == :not_checked
                                 -> p14_fuse == (:abstain, :ood_not_checked)
```

with the **positive half** beside it (present pool -> `:clear` -> `:decide`), so the chain cannot
pass by abstaining unconditionally.

### The named limit this lane must carry

`channels_wired = (:density,)`, `channels_not_wired = (:pp, :noise)`. The shipped flag is an
OR-fusion over three channels; only the density channel is wired here. A `:clear` resolved
downstream means **one channel says clear**, not three. This is recorded as a FIELD on every return
so it travels into the artifact.

## Verification

| check | result |
|---|---|
| `julia --project=spike spike/test/test_p14_pools.jl` | **exit 0**, 102/102, 74 s |
| every other `test_p14_*.jl`, run per file | **all exit 0** (9 files total) |
| Task 1 `<verify>` (600-item draw + `p14_assert_unstratified`) | exit 0, masses printed above |
| Task 2 `<verify>` (`p14_ood_reference`, basis assertion) | exit 0, `available=true`, finite `thr` |
| `grep -v '^\s*#' pools.jl \| grep -c 'stratify_by_class'` | 0 |
| `grep -v '^\s*#' pools.jl \| grep -c 'youden_j'` | 0 |
| `grep -v '^\s*#' pools.jl \| grep -c 'ood_nulls_8.jld2\|amended_v2'` | 0 |
| `git status --porcelain -- spike/data/cache/p13 artifacts` | empty |
| `git diff --exit-code HEAD -- src spike/Project.toml spike/Manifest.toml corpus spike/validation` | **exit 0** |

### Six mutations, each verified to land in the diff before the run

| mutation | caught by | failures |
|---|---|---|
| per-item RNG ignores its index | testsets 1 + 2 | 5 fail + 1 error |
| equal-thirds check disabled | testset 2 | 1 fail (on the MESSAGE) |
| `basis_provenance` relabelled | testset 4 | 2 fail |
| a threshold substituted on the failure path | testset 5 | 6 fail |
| the balancing helper called from `pools.jl` | testset 3 | 1 fail |
| the gate set marked `reported = true` | testset 6 | 1 fail |

Every mutation was reverted with `git checkout --` and the clean run re-confirmed at 102/102.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 — Bug] The plan's suggested RNG lambda would have produced `n` identical items**

- **Found during:** Task 1
- **Issue:** the plan's `<action>` suggests
  `p13_simulate_indices(1:n, basis; rng_for = i -> rng_for(counter))`. `p14_rng(counter)` returns a
  *fresh* Philox keyed by `(seed ⊻ salt, counter)` on every call, so that lambda hands **every item
  the same stream** and the draw is `n` copies of one item. It is silent: identical items are valid
  items, and reproducibility, counter-separation and even `p13_assert_classes_agree` all still pass.
- **Fix:** mirrored `p13_datagen_rng` (`spike/p13/datagen.jl:180`) instead — the activity counter is
  folded into the KEY word and the per-item GLOBAL INDEX goes in the second word, which is what
  makes a column a pure function of its index. Implemented as named functions
  (`p14_datagen_key_word` / `p14_datagen_rng`) rather than an inline lambda, so 14-01's key-word
  cross-product membership can be asserted on every draw.
- **Test added:** testset 1 now asserts within-draw distinctness; the mutation reproducing the
  plan's lambda was run and fails 5 assertions + 1 error.
- **Files modified:** `spike/p14/pools.jl`, `spike/test/test_p14_pools.jl`
- **Commits:** efc4334, c8c9130

**2. [Rule 3 — Blocking] Task 1's `<verify>` command needs ~11 minutes single-threaded**

- **Found during:** Task 1
- **Issue:** an image item costs ~1.13 s; the 600-item draw the verify command performs is ~680 s on
  one thread.
- **Fix:** ran it with `julia --project=spike -t auto` (32 threads, 125.5 s). This is
  byte-identical by construction — each item is a pure function of its index, which is the stated
  guarantee of `p13_simulate_indices` — and it was confirmed rather than assumed: the first three
  classes match the earlier single-threaded run exactly, and the realized masses match the
  theta-only twin exactly.
- **Commit:** efc4334

**3. [Rule 3 — Blocking] `p14_draw_classes` added, so the test can meet BOTH of the plan's bars**

- **Found during:** Task 3
- **Issue:** the plan wants a ~600-item pi-distribution check AND a test that finishes in under
  120 s. Those are incompatible with image draws.
- **Fix:** added `p14_draw_classes(n; counter, rng_for)`, the theta-only half of the same draw at
  the same indices — the label is a function of theta alone, so this is the same computation, not an
  approximation. It is pinned to the pool twice: by `p13_assert_classes_agree` inside
  `p14_draw_pool`, and by a test assertion that the first items agree.
- **Commit:** efc4334

**4. [Rule 2 — Missing] `basis_provenance = :unavailable` on failure paths**

- **Found during:** Task 2
- **Issue:** the plan pins `basis_provenance === :p13_pool_phase11_zt` only on successful returns
  and does not say what a failure returns.
- **Fix:** failure paths return `:unavailable`. A null that was never fitted has no basis, and
  labelling one anyway would be Pitfall 4 in miniature.
- **Commit:** efc4334

**5. [Rule 2 — Missing] An eighth testset, so `p14_misspec_pool` is actually executed**

- **Found during:** Task 3
- **Issue:** the plan's seven testsets never call `p14_misspec_pool`; a function nobody runs is a
  function that does not work.
- **Fix:** added a matched-arms testset (n = 1, fixture stream, `:noise` at the strongest rung),
  asserting equal arm lengths, identical lambda/imsize/class, differing summaries, finiteness, and
  that an unknown family throws.
- **Commit:** c8c9130

**6. [process] Tasks 1 and 2 share one file and are covered by ONE commit**

Both tasks write `spike/p14/pools.jl` and were written in a single pass. Each was verified
separately by its own `<verify>` command before the commit, and both are described in it. Task 3 is
its own commit.

## Assumption Drift (advisory)

**1. "Reuse `p13_real_id_reference` verbatim" — reused as a CONSTRUCTION, not as a call**

- **Planned:** the plan says to reuse `p13_real_id_reference`
  (`spike/p13/run_p13_realimage.jl:412`) *verbatim* rather than reimplementing it.
- **Actual:** that function lives in a RUNNER which drags `real_images.jl`, `alpha_series.jl` and
  CairoMakie into a file that needs one covariance — and its own docstring records that it is an
  *attributed copy* of the shipped `fit_ood_nulls` / `maha_score`, written that way for exactly this
  reason. So the same construction is expressed here through the canonical in-repo
  `fit_ood_nulls` / `maha_score` / `id_threshold`, which are the originals those copies attribute
  to. Nothing is re-implemented: three existing functions are called.
- **Why it is material, and why it is settled:** "verbatim" could have meant the call. The
  equivalence is not argued, it is measured — the threshold, the column count and the shard count
  all reproduce Phase 13's reported values exactly (167.54463293780395 / 96000 / 12).

**2. The misspecified fixture item scored LOWER than its in-distribution twin — an observation, not
a result**

- **Planned:** nothing in the plan says which way a misspecification arm should move.
- **Actual:** the single fixture item (family `:noise`, level 4) scored 0.203 on the density channel
  against 11.13 for its matched in-distribution twin, both far below `thr = 167.54`.
- **Why it is recorded here and NOWHERE else:** n = 1 measures nothing and this plan pre-registers
  no bar for it — SC3-d is a later runner's business at 1,000 items per arm. It is noted because it
  is consistent with the recorded history that the density channel alone is weak and the NOISE
  channel carried OOD detection, which is precisely why `channels_not_wired = (:pp, :noise)` is
  carried as a named limit rather than a footnote. **Do not quote this as a finding.**

## Known Stubs

None. Every function in `spike/p14/pools.jl` is executed by the test file, including
`p14_misspec_pool` and `p14_pool_provenance`.

## Self-Check: PASSED

- `spike/p14/pools.jl` — FOUND (699 lines)
- `spike/test/test_p14_pools.jl` — FOUND (362 lines)
- commit efc4334 — FOUND
- commit c8c9130 — FOUND
