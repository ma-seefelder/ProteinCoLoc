---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 17
status: complete
subsystem: amortized-inference
tags: [d-15-amended, d-16, pre-registration-amendment, real-images, include-guard, def-12-03, consts-sha, reported-not-gated, named-limits, retraction, red-by-design]
requires:
  - "spike/p13/three_way_net.jld2 (13-11) -- the trained net, READ ONLY; its train-time consts_sha is now the pinned pre_amendment literal"
  - "spike/p13/real_images.jl (13-15) -- pair-agnostic ingestion; needed no executable change"
  - "spike/p13/consts.jl (13-01/13-10) -- the Tier-1 pre-registration, AMENDED by this plan"
  - ".planning/phases/13-three-hypothesis-amortized-bayes-factor/13-D15-AMENDMENT.md -- the frozen authorisation"
provides:
  - "P13_REAL_CHANNEL_PAIR = (2, 3) -- the operative green/red target pair"
  - "P13_REAL_ANCHOR_MBAR = (positive = 0.4603, negative = 0.3815) -- unmasked c2/c3 anchors"
  - "P13_CONSTS_SHA (spike/p13/net.jl) -- the named, dated pre/post-amendment sha pair asserted at all three runner sites"
  - "P13_DECLARED_DEVIATIONS as the Tier-1 include-guard sentinel -- DEF-12-03 closed"
  - "spike/p13/realimage_report.jld2 -- the re-run on the amended pair, carrying superseded_* comparison fields"
  - "13-REPORT.md Limit E -- the negative-mu regime is NOT DEMONSTRABLE on this substrate; the constructibility claim is RETRACTED"
  - "spike/test/runtests.jl P13_KNOWN_RED ledger -- the Phase-13 block survives MORE THAN ONE deliberately-red file"
affects:
  - ".planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md (sections 8b, 8c, 9, 13, 14)"
  - ".planning/ROADMAP.md (Phase-13 row, 13-16 bullet, 13-17 bullet, Verdict real-arm clause)"
tech-stack:
  added: []
  patterns:
    - "guard an include on a name the target file EXCLUSIVELY owns -- never on a seed or a prior bound (CONVENTIONS.md C-01)"
    - "a train-time provenance guard pinned to NAMED, DATED literals on BOTH sides, never a disjunction"
key-files:
  created: []
  modified:
    - spike/p13/consts.jl
    - spike/p13/net.jl
    - spike/p13/labels.jl
    - spike/p13/preconditions.jl
    - spike/p13/result.jl
    - spike/p13/tau_probe.jl
    - spike/p13/toy_gaussian.jl
    - spike/p13/run_three_way_gate.jl
    - spike/p13/run_alpha_series.jl
    - spike/p13/run_p13_realimage.jl
    - spike/p13/real_images.jl
    - spike/p13/realimage_report.jld2
    - spike/test/test_p13_consts.jl
    - spike/test/test_p13_real.jl
    - spike/test/runtests.jl
    - .planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md
    - .planning/phases/13-three-hypothesis-amortized-bayes-factor/13-15-SUMMARY.md
    - .planning/phases/13-three-hypothesis-amortized-bayes-factor/13-16-SUMMARY.md
    - .planning/ROADMAP.md
    - .planning/STATE.md
decisions:
  - "M1 (amend consts.jl in place) rather than M2 (isolated-module read) -- ruled by the user; the byte-lock breaks either way so CHANGE B rides along at zero additional pre-registration cost"
  - "The consts_sha guard is a NAMED, DATED pre/post literal pair asserted on both sides, NOT a widened disjunction that would accept future drift"
  - "Three test_p13_real.jl assertions broken by the corrected pair are LEFT FAILING rather than rewritten -- the disposition was escalated, and the user RULED on 2026-07-29 that they STAY RED as named limits"
  - "The D-16 real-substrate expectation is NOT amended: a green test standing beside a 'not shown' conclusion is worse than a red one"
  - "The masking is fixed STRUCTURALLY (wrap + ledger), not by re-ordering -- the last-position pattern supports exactly ONE throwing file and there are now two"
metrics:
  duration: "~2 h (13-17 as planned) + ~1.5 h (the ruling: retraction, masking fix, close)"
  completed: 2026-07-29
---

# Phase 13 Plan 17: The D-15 Pre-Registration Amendment and the 13-16 Re-Run — Summary

Amended the byte-locked Phase-13 pre-registration in **one edit carrying two separately justified
changes**, re-derived the train-time `consts_sha` guard so it means *more* than before, and re-ran
the real-image arm once on the corrected substrate — which **left the conclusion unchanged and
retracted two claims the phase had made**. Under the user's ruling of 2026-07-29 those retractions
are now **named limits carried as RED TESTS**: the assertions stay failing because the phase's own
verdict says the real arm is *not shown* to work, and the harness was fixed structurally so that
**more than one deliberately-red file can coexist without masking its siblings**.

---

## ✅ STATUS: `complete` — and what had to be true first

This plan stood at `partial` for one reason: it surfaced a pre-registration question it deliberately
did not resolve, and three test assertions were failing as a result. **The user ruled on
2026-07-29** and the ruling has been executed. `complete` is claimed only because all three of the
following are now true, each verified rather than asserted:

1. **The named limits are RECORDED** — the retraction written plainly in `13-REPORT.md` §8c, carried
   as **Limit E** in §9, inventoried assertion-by-assertion in **§9a**, with the ruling verbatim in
   §14.9 and an additive `STATE.md` block. See **§8** below.
2. **The masking is FIXED** — structurally, verified by a full-suite run each side. See **§7a**.
3. **The verification can be stated honestly** — and it is stated as *"the suite runs to completion,
   every Phase-13 file executes and reports, and it still exits non-zero on exactly the five
   named-limit assertions."* **No full-suite green is claimed, because there must not be one.**

The plan index counts a plan complete on a file's mere existence, and this project has been damaged
by that before. Nothing here rests on existence: every claim below carries the command output it
came from.

---

## 1. The ruling this plan implements, recorded verbatim

The user ruled on 2026-07-29:

> **M1 — fix at the source, inside the same authorised amendment. M2 (the isolated-module read) is
> the tool for files that CANNOT be edited because no amendment authorises opening them.**

The reasoning, reproduced because *a mechanism chosen for a reason is auditable and one chosen by
default is not*:

> *The byte-lock breaks either way, and `h.consts_sha == p13_consts_sha()` has to be re-derived in
> the gate, alpha and realimage runners regardless — so the include-guard fix rides along at **ZERO
> additional pre-registration cost**. That coupling is precisely why the answer is M1 and not M2.*

---

## 2. TWO CHANGES, TWO JUSTIFICATIONS — disclosed separately, on the page

Both rode in one edit to `spike/p13/consts.jl`. Neither is a threshold moved after proving
inconvenient. **No seed, bar, floor, band, ladder, tolerance or allowance changed in either**, and a
`-U0` diff over changed `const` **declaration** lines shows exactly three: the pair, the dropped
redundancy pair, and the anchors. No seed literal appears in any changed line at all.

| | **CHANGE A — the channel pair** | **CHANGE B — the include-guard sentinel** |
|---|---|---|
| What moved | `P13_REAL_CHANNEL_PAIR` `(1,2)` → **`(2,3)`**; `P13_REAL_REDUNDANCY_PAIR` **DROPPED**; `P13_REAL_ANCHOR_MBAR` → **`0.4603 / 0.3815`** (UNMASKED) | the Tier-1 body wrapper and **all eight callers** guard on **`:P13_DECLARED_DEVIATIONS`** instead of `:P13_DEV_SEED` |
| What was wrong | it **named the WRONG PHYSICAL OBJECT** — `c1` is the DAPI/Hoechst nuclear counterstain (`test/runtests.jl:105`), so both pairs measured counterstain-versus-protein overlap, which is not colocalization | it guarded on a name **another phase must LEGITIMATELY MIRROR** — `p12_consts.jl:109` declares `P13_DEV_SEED` to assert seed disjointness, so a full-suite run skipped the whole body and ~104 Tier-1 constants went undefined (DEF-12-03) |
| Any bar move? | **No** | **No — it changes no value at all**, only which symbol an `if` tests |
| Independent justification | the upstream `ghat.jl:59` had **already** been corrected to `(2,3)` on 2026-07-25, **39 h before `consts.jl` was written**; this closes Phase 13's drift from its own frozen calibration source | it is the `fb76b84` precedent (*guard on the name it actually owns*) applied to the one case where the owner file is the file needing the edit; diagnosed from a **suite abort**, not from any Phase-13 number |

**`spike/validation/p12_consts.jl` was NOT touched and is byte-unchanged**; its `P13_DEV_SEED`
mirror is **correct**, because asserting seed disjointness requires naming the seeds.
`P13_DEV_SEED = 0x0000_0000_0B13_DE71` is byte-unchanged at its declaration.

**The amendment took the WORSE number, and it survives.** Unmasked, the corrected pair separates the
fixtures by **0.0788** — against **0.0811** for the superseded pair, and against **0.8576** for the
masked read of the *same* corrected pair, which `ghat.jl:69-71` itself calls the sharper separator
and which no prior source forbids. It is refused on a code property: `patch_summary` applies no Otsu
mask, so the net was trained on unmasked summaries and masked ones would be off-distribution input.
**The masked `+0.8238 / −0.0338` appear nowhere in the arm except in prose stating they are
forbidden and why** — verified: 0 occurrences in the comment-stripped source of `consts.jl`.

---

## 3. Provenance — every value copied from live command output

| | Pre-amendment | Post-amendment |
|---|---|---|
| `sha256(read("spike/p13/consts.jl"))` | `100e97a37bb470bb7cb4fbdbd8e33f219d83adcf61fe398a90cf3ad3391f3fa6` | `d6a6e63bf19bb13ad5c61eeaaf9fc18fab27c5aa7ba879bad5c493a5c71ae247` |
| `git hash-object --no-filters` | `5a4ea2223ce11bcae2f974c643ec1e523458e8ac` | `ef117e6ebc42095b777320ae0ceb680c13d7c6b9` |
| `git rev-parse HEAD:spike/p13/consts.jl` | `70fe66df33832af2e648ed3f81cc3ad356ce4cd8` | `d82caea64e19e8509f3abc1c09140cc7c2e2d497` |
| `git log --oneline -- spike/p13/consts.jl` | two commits: `c42cc8e`, `bfba6ac` | **three** — the amendment commit is **`9fe9b02`** |

**Cross-checked before any byte moved:** the 13-11 net's own recorded `h.consts_sha` equals
`100e97a3…`, so the plan's premise held.

**The guard was re-derived to mean MORE, not less.** `spike/p13/net.jl` defines
`P13_CONSTS_SHA = (pre_amendment = …, post_amendment = …)` and each of the three runners now asserts
**both sides** against those named literals. The old form (`h.consts_sha == p13_consts_sha()`) only
required the two sides to *agree with each other*, so a retrain plus an undisclosed edit would have
passed silently; a **third** sha now fails. The widened `||` form is **rejected in a code comment**.
No guard was deleted or blanket-disabled. `gate_report.jld2` and `alpha_report.jld2` are
byte-unchanged, and the two gate-arm runners changed in **two hunks each — the banner and the
assertion, both the sha re-derivation and nothing else** (verified by hunk inspection).

---

## 4. CHANGE B verified — and the honest qualification

| | BEFORE | AFTER |
|---|---|---|
| **Full suite** abort | **`runtests.jl:174`** → `test_p12_suite.jl:86` → `test_p12_decoupling.jl:150` (P12 decoupling, 25 pass / 1 fail). **Phase 13 never reached at all.** | **`runtests.jl:223`** → `test_p13_real.jl`. Phase-13 block **executes**: `consts` **137/137, 0 errors**, `labels` 164/164, `alpha` 129/129 |
| suite exit | 1 | 1 |
| **Targeted DEF-12-03 repro** (mirror loaded, then `test_p13_consts.jl` — the exact suite condition) | **22 pass / 1 fail / 114 error**, every error an `UndefVarError` | **137 pass / 0 error** |

**The raw before/after suite comparison is confounded and is NOT claimed.** The BEFORE run did not
abort where the amendment predicted — it aborted in the **Phase-12** block — and a concurrently
executing Phase-12 agent fixed that between the two runs (`8880c29`). **That movement is Phase 12's,
not CHANGE B's.** The targeted reproduction is the clean evidence, and it is exact.
`P13_DECLARED_DEVIATIONS` was confirmed **not** defined by the Phase-12 mirror.

**The §5B.5 knock-on was VERIFIED, not assumed.** Loading the narrow (`UInt32`) declarations first
and then letting the now-un-skipped Tier-1 body re-declare them wide (`UInt64`) executed **without
error** on Julia 1.12.6; values agree, only widths differ. **The STOP RULE did not trigger and no
seed literal was edited.**

---

## 5. The re-run — one run, and what it says

`julia --project=spike spike/p13/run_p13_realimage.jl`, exit 0, 1.33 min, on the main working tree.
**No RNG stream consumed**, so no reserved counter (1 τ, 2 datagen, 3 gate, 4 α, 5 continuity,
99 fixtures) could collide; `rng_stream = "NONE"` is persisted and re-asserted from the runner's own
source text.

### DOES THE CONCLUSION CHANGE? — answered in the same words the report uses

> **The channel pair was wrong, it has been corrected under an authorised amendment, the arm was
> re-run once — and the conclusion is unchanged.** The numbers moved, one verdict moved, and the
> overall reading did not: this arm is still **qualitative, n = 2 specimens, no ground-truth label,
> OOD-bound, weak evidence — coherence and never correctness.**

**That the limitation survived its own correction is the point.** It shows the weakness is
**structural** — two unlabelled specimens outside the training distribution — rather than an
artefact of having read the counterstain channel.

**This is not a rescue, and three measured facts say so:** the binding OOD limit got **worse**
(4.198× against 2.491×); the amendment adopted the **narrower** separation; and the phase now
reports **one arm where it reported two**.

### The delta table — every row labelled with the pair it was measured on

| Quantity | SUPERSEDED — pair `(1, 2)` | **AMENDED — pair `(2, 3)`** |
|---|---|---|
| `m-bar` positive / negative | +0.32916 / +0.24805 | **+0.46027 / +0.38147** |
| headline `log BF(C:R)` / `(E:R)`, positive as sample | −0.4658 / −10.2727 | **+5.70311 / −12.9744** |
| headline `log BF(C:R)` / `(E:R)`, negative as sample | −3.9432 / −9.0534 | **−8.65723 / −12.8627** |
| descriptive argmax, positive / negative as sample | RANDOM / RANDOM | **COLOC / RANDOM** |
| Phase-13 OOD density / threshold / ratio | 417.2974 / 167.5446 / 2.491× | **703.2995 / 167.5446 / 4.198×** |
| per-acquisition densities (positive / negative) | 417.2974 / 321.1505 | **703.2995 / 398.7601** |
| shipped OOD density / threshold / ratio | 433.6884 / 179.1368 / 2.421× | **948.9745 / 179.1368 / 5.297×** |
| A12 outcome | `agreement` | **`agreement`** |
| λ response span, exclusion / coloc heads | 0.0117 / 0.0585 ; 0.2595 / 0.1750 | **0.0577 / 0.0731 ; 0.1747 / 0.1155** |
| `alpha_star_real` positive / negative | 0.875 / 0.625 | **`nothing` / `nothing`** |
| ladder `m-bar` span, positive | +0.32916 → **−0.16787** | **+0.46027 → +0.13806** (stays positive) |
| ladder `m-bar` span, negative | +0.24805 → **−0.25198** | **+0.38147 → +0.20414** (stays positive) |
| realized mask fraction positive / negative | 0.134933 / 0.229189 (**on c1**) | **0.051150 / 0.032858** (**on c2**) |
| eight α invariants | all hold, worst residual 2.86e-16 | **all hold, worst residual 1.80e-16** |
| ρ through the frozen `ghat`, positive / negative | +0.33042 / +0.21480 | **+0.48253 / +0.38434** |
| D-05 contrast vs `tau = 0.15` | 0.11562 — inside | **0.09819 — inside** |
| redundancy arm | reported, 607.5728 (3.63×) | **DROPPED — no second pair exists** |

**The headline OOD figure is reported as a CHANGE, never restated.** It was already demonstrably
pair-dependent: the dropped `(1, 3)` pair had landed at 607.5728, a 46% move.

**The mask fraction stayed inside `P13_ALPHA_MASK_FRACTION_BOUNDS = (0.01, 0.40)`** on the c2 mask,
so the runner did **not** take its `INVARIANT VIOLATION` branch. The band was not widened and was
never at risk of being.

---

## 6. Artifacts, and what regenerating them costs

- **`spike/p13/realimage_report.jld2`** — TRACKED and committed (`652692c`). 65 keys. Carries
  `is_gated = false`, `qualitative_only = true`, the amended pair `[2, 3]`, the amended anchors,
  `rng_stream` beginning `NONE`, the `superseded_*` comparison fields and the dropped-arm record.
  Its `consts_sha` is the **post**-amendment value and its `net_consts_sha` the **pre**-amendment
  one — the deliberate, asserted divergence.
- **`spike/figures/p13_realimage.png`** — **GITIGNORED** (`.gitignore:401`, the house `*.png` rule),
  275,259 bytes, regenerated by the re-run but **NOT committed and NOT recoverable from git**.
  Regenerating it costs a **full ~1.33 min re-run** of `run_p13_realimage.jl`. The substrate is
  deterministic and consumes no RNG stream, so the regeneration is byte-identical — but it is a
  re-run, not a free copy.
- **Read-only proof holds:** the `git ls-files -s test/test_images` digest is
  `eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276` before **and** after, and
  `git status --porcelain test/test_images/` is empty. `spike/data/cache/p11/` was neither read nor
  written.

---

## 7. Verification — which is suite-level and which is per-file

**Be precise about this, because a signed-off checklist asserting something that had not happened
has cost this project before.**

**Suite-level (observed, `julia --project=spike spike/test/runtests.jl`):** the Phase-13 block now
**executes** rather than being masked — `test_p13_consts.jl` 137/137 with **0 errors**,
`test_p13_labels.jl` 164/164, `test_p13_alpha.jl` 129/129. **The suite still exits 1** and that is
correct (the Phase-4 `SPEEDUP_GATE` is `@test_skip`-paused, `test_p13_correction.jl` carries two
committed pre-registered measured misses). **A full-suite green was neither achieved nor claimed.**

**As of the ruling this is HISTORY: the suite no longer aborts there.** What follows is the state
the ruling inherited — the suite aborted at `runtests.jl:223`, on `test_p13_real.jl`, which this
plan caused (three deliberate measured misses, §8). It masked the **seven** Phase-13 files after it,
so **all seven were run individually** and are reported here rather than assumed. **§7a records what
replaced it, and every one of these files now also runs inside the suite.**

| File | Observed | Exit |
|---|---|---|
| `test_p13_consts.jl` | 137 / 137 | 0 |
| `test_p13_real.jl` | 273 pass / **3 fail** | **1** — the deliberate misses |
| `test_p13_tau.jl` | 99 / 99 | 0 |
| `test_p13_net.jl` | 93 / 93 | 0 |
| `test_p13_result.jl` | 43 / 43 | 0 |
| `test_p13_calibration.jl` | 95 / 95 | 0 |
| `test_p13_preconditions.jl` | 74 / 74 | 0 |
| `test_p13_datagen.jl` | 706 / 706 | 0 |
| `test_p13_correction.jl` | 88 pass / **2 fail** | **1** — pre-existing, committed, pre-registered |

---

## 7a. The masking fix — measured, and stated without a win it did not earn

**The house pattern was documented, and it is no longer sufficient.** `spike/test/runtests.jl` said,
in so many words, that the correction arm is last *"deliberately so … a thrown testset aborts the
remaining includes, so any sibling placed after it would silently never run."* **That pattern
supports exactly ONE throwing file.** There are now two — `test_p13_real.jl` (3 named-limit
failures) and `test_p13_correction.jl` (2 pre-registered F5 misses) — so whichever runs first masks
the other. **Re-ordering cannot fix that**; it only chooses which red is hidden, and it re-creates
the *"before the file that throws today"* contingency that already failed once in this same harness
when `test_npe.jl` threw ahead of the file the wiring had been reasoned against.

**The fix is structural.** Each known-red include is wrapped; a **ledger testset** then asserts the
observed red set equals the expected one, and the recorded exception is re-raised after the last
include. Four required properties, each **verified empirically rather than trusted**:

| Property | Verification |
|---|---|
| 1. Every Phase-13 sibling runs and reports | full-suite run: **11** Phase-13 file-level testsets report, against **4** before |
| 2. The named limits still surface, not swallowed, not downgraded | the log carries **exactly 5 `Test Failed at` lines** with unchanged `Evaluated:` values; `Test` prints each failure and its summary table *before* the testset throws |
| 3. The suite still exits non-zero | **exit 1**, ending on `ERROR: LoadError: Some tests did not pass: 273 passed, 3 failed` |
| 4. An unexpected pass fails loudly | exercised on a mock with one known-red file swapped for a green one: the ledger failed with `"second.jl" => false == "second.jl" => true`, exit 1 |

Plus the non-negotiable safety property: **anything that is not a `Test.TestSetException` is
rethrown immediately and never recorded as an expected red.** `include` wraps the throw in a
`LoadError` — **verified on Julia 1.12.6, not assumed** — so the wrapper is unwrapped before the
type test. Exercised on a mock raising a genuine `error(...)`: it propagated, was not recorded, and
aborted the run, which is correct for a real defect.

### One full-suite run each side

| | BEFORE | AFTER |
|---|---|---|
| Phase-13 file-level testsets reporting | **4** (`consts`, `labels`, `alpha`, `real`) | **11** — all of them, plus the ledger |
| Phase-13 assertions reported | 703 pass / 3 fail | **1,901 pass / 5 fail** + the ledger's 3 pass |
| Phase-13 files masked | **7** | **0** |
| aborts at | `runtests.jl:223` | runs to completion; re-raises after the ledger |
| `Test Failed at` lines in the whole log | 3 | **5** |
| non-Phase-13 failures | 0 | 0 |
| `test_p12_consts.jl` testset 9 (the ordering assertion) | 44 / 44 | **44 / 44** |
| suite exit | **1** | **1** |

**Which verification is suite-level and which is per-file, stated precisely.** Everything in the
table above is **suite-level** (`julia --project=spike spike/test/runtests.jl`, one run each side).
The two known-red files were **additionally** run per-file and agree exactly with their in-suite
numbers: `test_p13_real.jl` **273 / 3** exit 1, `test_p13_correction.jl` **88 / 2** exit 1. The
four ledger properties were verified **suite-level for 1–3** and **on an isolated mock for 4 and for
the rethrow-a-genuine-error path**, because neither can be provoked on the real suite without
editing an assertion the ruling forbids touching.

**NO FULL-SUITE GREEN IS CLAIMED, AND THERE MUST NOT BE ONE.** The honest statement is: *the suite
runs to completion, every Phase-13 file executes and reports, and the suite still exits non-zero on
exactly the five named-limit assertions.*

**The concurrency is disclosed rather than assumed away.** A Phase-12 executor was live throughout
and landed five commits (`986d518` … `406015e`) **between** the two runs. **None touched any file
under `spike/test/` or `spike/p13/`** — they touched `.planning/` documents,
`spike/validation/run_p12_stage1_ridge.jl` and `spike/validation/p12_stage1_report.jld2` — and both
runs report zero Phase-12 failures with identical Phase-13 per-file numbers. The measured delta is
attributable to `spike/test/runtests.jl`, the only file this change edited. **That is a narrower
claim than "the fix caused the whole difference", and it is the one the evidence supports.**

**What did NOT change, checked line by line:** no assertion expression, no `@test_skip`, no
`@test_broken`, no threshold, no seed, no bar, no floor, no band, no allowance. **No include was
re-ordered** — `test_p12_suite.jl` keeps FIRST position and the two known-red files keep their
existing positions. The non-comment diff is 3 `const`s, one 6-line helper, two `try`/`catch`
wrappers, one ledger testset and one re-raise block.

---

**Decoupling, re-verified after every task:** `git diff --quiet HEAD --` on `src/`, `test/`, `docs/`,
`corpus/`, `spike/Project.toml`, `spike/Manifest.toml`, `spike/validation/`,
`spike/p13/gate_report.jld2`, `spike/p13/alpha_report.jld2` all exit 0. `spike/test/runtests.jl` is
byte-unchanged. `P13_ITERATION_ALLOWANCE` reads **1 of 1, UNSPENT**. The Phase-16 seal is shut. The
not-a-gate exemption set is still exactly `{P13_REAL_ANCHOR_TOL, P13_REAL_OOD_SHIPPED_THRESHOLD}`.

---

## 8. ✅ THE BLOCKER, RULED AND CLOSED — the red stays, and it is the finding

**The user ruled on 2026-07-29: accept the three misses as NAMED LIMITS, leave them RED, and do NOT
amend the D-16 real-substrate expectation.** Recorded verbatim, because the reasoning *is* the
deliverable:

> *Phase 13's own two-sentence verdict already says the three-way Bayes factor works on simulated
> data and is NOT SHOWN to work on real microscopy. A test that asserts the real-substrate
> expectation and FAILS is therefore TELLING THE TRUTH. Amending it to expect the new measurement
> would produce a GREEN TEST STANDING NEXT TO A CONCLUSION THAT SAYS "NOT SHOWN" — a test that
> passes while the science says otherwise is worse than a red one, and this project's credibility
> rests on exactly that not happening. **The red is not a defect to be cleared; it is the finding,
> encoded where a future reader will trip over it.***

**Verified at source, not taken on trust:** that two-sentence verdict is in `STATE.md` (the 13-14
entry) and in `13-REPORT.md`'s header — *"works on simulated data … and it is not shown to work on
real microscopy."* The tests and the conclusion now agree.

### 8.1 THE RETRACTION, written plainly

**The manuscript-bound sentence *"negative induced μ is CONSTRUCTIBLE from real microscopy pixels
via the D-16 mask-based reassignment"* is RETRACTED.**

**The reason: it was an ARTEFACT of segregating a target channel against a nuclear counterstain.**
The α mask is the Otsu mask of the pair's FIRST channel — nuclear `c1` before, green target `c2`
now. Moving green out of nuclei reaches negative correlation; moving red out of green, two channels
that co-occur on a shared bright background, does not.

**The evidence, read out of `spike/p13/realimage_report.jld2` rather than quoted from a document:**
`alpha_star_real` = `(positive = nothing, negative = nothing)` against `superseded_alpha_star_real`
= `(positive = 0.875, negative = 0.625)`; and `m-bar` runs **+0.4602658 → +0.1380565** and
**+0.3814746 → +0.2041410**, **staying POSITIVE on both fixtures at all nine rungs**.

**The correction COST the claim rather than revealing a new one** — the same register as
`13-D15-AMENDMENT.md` §2.1, where the amendment took the worse number at the one point a better and
permitted one was on the table. It removed a claim and put nothing in its place. **This is not a
discovery and it is not presented as one.** Phase 13 does *not* assert the converse either: two
fixtures cannot establish that the negative regime is physically unreachable, only that it is **not
demonstrable on this substrate**. The phase claims strictly less than it did.

### 8.2 THE CORROBORATION — two independent routes, verified at source

This is what makes the retraction robust rather than one surprising measurement. Both were checked
against `.planning/STATE.md`'s *Quick Tasks Completed* table **and** against
`spike/simulator/ghat.jl` at HEAD, not cited from memory:

| Quick task | Date | Commit | What it established |
|---|---|---|---|
| **260725-vl8** | 2026-07-25 | `151ad79` | The unqualified *"negative tail not physically reachable / PRIOR-ONLY"* claim in `ghat.jl` is **WITHDRAWN — not inverted.** |
| **260725-wb7** | 2026-07-25 | `78dc37f` | `neg_reachable` flipping to `true` on a masked reading of **−0.0338 ≈ 0** is a **PREDICATE ARTEFACT, not evidence**; the vl8 withdrawal **STANDS**. |

`spike/simulator/ghat.jl:84-91` still carries both today: *"STILL NOT ESTABLISHED either way … the
former unqualified 'physically reachable = false ⇒ PRIOR-ONLY' reading stays WITHDRAWN and is NOT
reinstated; equally, the clean +0.82 / −0.03 control separation is NOT evidence for the opposite
claim."*

**Not overclaimed:** neither route says the negative regime is unreachable. They say **this substrate
carries no evidence that it is reachable**, and that the one number that looked like such evidence
was a predicate artefact on a value indistinguishable from zero. The α-ladder is the third route and
points the same way. **Three routes, one direction** — and the agreement is worth more than any one
of them alone.

### 8.3 The D-05 disagreement stays RED, and it is an INFORMATIVE red

Contrast **0.09819** is inside `tau = 0.15`, so the label rule assigns **RANDOM** while the net's
argmax says **COLOC**. Left standing: the net is not retrained, `tau` is not moved, the label rule is
not re-derived.

**A disagreement on a fixture flagged OOD at 4.198× its own threshold is close to what one should
EXPECT.** Both fixtures score above the maximum density of all 96,000 simulated acquisitions the null
was fit on. A net reading an input that far outside its training distribution and *still* agreeing
with a rule derived from the same summary would have been the surprising outcome. **It is the OOD
flag surfacing a second time in an independent place** — an informative red, not an anomalous one.

### 8.4 The three assertions, as they stood when the ruling arrived

The three `test_p13_real.jl` assertions were **left FAILING rather than rewritten**, because
rewriting a substrate expectation to match what was measured *after* the measurement is exactly the
act the amendment exists not to be. **The ruling confirms that disposition and makes it permanent.**

1. **The D-16 α-ladder no longer crosses zero.** `alpha_star_real` is `nothing` on **both**
   conditions (was 0.875 / 0.625); `m-bar` runs +0.46027 → +0.13806 and +0.38147 → +0.20414,
   entirely positive. Strict monotone decrease still holds — only the reach into negative
   correlation is gone. **So the manuscript-bound sentence *"Negative induced μ is CONSTRUCTIBLE
   from real microscopy pixels via the D-16 mask-based reassignment"* is RETRACTED.** It was
   measured on the counterstain pair. The mechanism is not mysterious: the α mask is the Otsu mask
   of the pair's **first** channel — nuclear `c1` before, green target `c2` now — and moving red out
   of green objects, when the two co-occur on a shared bright background, cannot drive the
   correlation negative the way moving green out of nuclei could. The c2 mask fractions are also far
   smaller (0.0511 / 0.0329 against 0.1349 / 0.2292).
2. **The D-05 coherence check that previously PASSED now DISAGREES.** Contrast **0.09819** is inside
   `tau = 0.15`, so the phase's own label rule assigns RANDOM — while the net's descriptive argmax
   says **COLOC** on the positive direction. `13-REPORT.md` previously called this *"a consistency
   check that passed"*.
3. **A third assertion, `REAL_ZEROS[c] > 0`, is a non-degeneracy guard** that now reads `0 > 0`
   because the transformed channel is `c3`, which carries no source zeros on one condition. The
   invariant it guards (zero count preserved) still **passes** on both conditions.

**Where the ruling landed each of them.** Nothing was moved to resolve any of it:
`P13_ITERATION_ALLOWANCE` reads **1 of 1 UNSPENT**, `spike/p13/consts.jl` is byte-unchanged at git
blob `d82caea6`, and no seed, bar, floor, band or tolerance was touched.

| # | Assertion | Ruling |
|---|---|---|
| 1 | `test_p13_real.jl:292` (fails twice, once per condition) — the α-ladder crossing | **STAYS RED.** Named as **Limit E**; the claim it encodes is retracted (§8.1). |
| 2 | The D-05 coherence disagreement — **not a test assertion**, a report-level finding | **STAYS as a disagreement**, recorded as an informative red (§8.3). |
| 3 | `test_p13_real.jl:226` — `REAL_ZEROS[c] > 0`, the non-degeneracy guard | **STAYS RED.** Reads `0 > 0` on the negative condition because `c3` carries no source zeros there. The invariant it guards still **PASSES** on both conditions. |

---

## 9. DEF-12-03, and the 28 guards NOT fixed

**DEF-12-03 is CLOSED BY REFERENCE.** `.planning/phases/12-spatial-colocalization-map/deferred-items.md`
was **NOT edited**: a live Phase-12 executor owns that file and a cross-phase write for a one-line
status change is not worth the merge. The Phase-12 owner can mark it when it next opens the file.

**The other 28 poisoned include guards found by the same mechanical sweep (220 guards, 118 files, 29
poisoned, 8 families) were DOCUMENTED and NOT FIXED**, with a remedy each in
`.planning/CONVENTIONS.md` **C-01 appendix**. Fixing another phase's frozen file on this executor's
own initiative is precisely the act this project's amendment discipline exists to prevent. No
thirtieth was found during execution.

---

## 10. Deviations from plan

### Auto-fixed / judgement calls

**1. [Rule 2 — completeness] Two further stale guard-describing comments in `consts.jl` were
corrected.** *Found during:* Task 2. The plan assigned only the commented illustration at `:93`.
`:79` (*"Guarded as ONE Tier-1 block keyed on `:P13_DEV_SEED`"*) and `:840` (*"a DIFFERENT sentinel
from the Tier-1 block's `:P13_DEV_SEED`"*) described the guard and became **false** the moment the
sentinel moved. Leaving them would have gone on teaching the bug the edit just fixed — the same
drift class the amendment repairs. Comment-only; no value changed. *Commit:* `9fe9b02`.

**2. [Rule 2 — honesty] `P13_REAL_OOD_SHIPPED_DENSITY` was labelled in `consts.jl` as a
superseded-pair measurement.** *Found during:* Task 2. §6.5 says it stays as a frozen historical
reference; a future reader had no way to know from the file that 433.69 was measured on `(1,2)`.
Comment-only; the constant is byte-unchanged. *Commit:* `9fe9b02`.

**3. [Rule 3 — blocking] Line endings.** The three runners were rewritten with LF by a scripted
patch and were restored to the repo's CRLF working-tree convention. `core.autocrlf=true` and
`* text=auto` mean git normalizes the blob either way, so no spurious change landed (verified: 21
changed lines per runner, not 875). `consts.jl` was never converted, so its post-amendment sha256 is
valid.

**4. [judgement] Three failing assertions left failing rather than rewritten.** Documented in place
with the reason, and escalated (§8). This is a decision *not* to act; the assertion expressions are
byte-unchanged. **The user ruled on the escalation and the disposition is now permanent: they stay
RED as named limits and the D-16 real-substrate expectation is NOT amended.**

**9. [Rule 2 — honesty] `13-REPORT.md` §9's *"A fifth entry was considered and NOT added"* paragraph
was re-worded to *"A further named-limit entry …"*.** *Found during:* the ruling. Adding **Limit E**
made a literal ordinal go stale the moment it landed — the same class of stale describing-text as
deviations 1 and 2 above. **Wording only; the evaluated non-finding it records (exclusion AUC
0.997257 in the deep tail, so the anticipated prior-atom limit's condition did not hold) is
unchanged**, and no measurement was touched.

**10. [judgement] Property 4 of the masking fix, and the rethrow-a-genuine-error path, were verified
on an ISOLATED MOCK rather than on the real suite.** *Found during:* the ruling. Neither can be
provoked in the real suite without editing an assertion the ruling forbids touching — making a
known-red file pass, or making one raise a non-test error. The mock reproduces the exact ledger code
against a green file and against a file that calls `error(...)`. **Reported as a mock, never as a
suite-level observation** (§7a).

### Plan-criterion defects, reported rather than overridden (binding constraint 10)

**5. `M` must be `0` for the two gate-arm runners.** The whitelist regex
(`consts_sha|p13_consts_sha|P13_CONSTS_SHA|13-D15-AMENDMENT|^[+-]\s*(#|$)`) does not match three
banner-continuation `println` lines that say `consts.jl` rather than `consts_sha`. Observed `M = 3`.
**The intended property was verified directly and holds:** each gate-arm runner's diff is exactly
**two hunks** — the sha banner and the sha assertion — and nothing else.

**6. `grep -c 'channels = P13_REAL_CHANNEL_PAIR' spike/p13/real_images.jl` must output `5`.**
It outputs **10**, because five are function-signature keyword defaults and five are docstring
echoes of those signatures. **It also output 10 at HEAD**, so nothing moved; the five signature
defaults are intact and byte-identical.

**7. `git diff -U0 … 13-REPORT.md | grep -c <gate figures>` must be `0` for `.planning/ROADMAP.md`'s
Verdict block.** Unsatisfiable as written: the Verdict is a **single line** carrying both the gate
clause and the real-arm clause, and the same plan instructs that the real-arm clause be amended.
**Verified directly instead:** the gate clause is **byte-identical** across the two revisions (418
bytes, equal).

**8. The BEFORE suite did not abort where the plan predicted** (`runtests.jl:220` with
`test_p13_consts.jl` at 22/1/114). It aborted earlier, at `runtests.jl:174`, in the **Phase-12**
block, so Phase 13 was never reached. Recorded as observed. A targeted reproduction was added — and
captured **before** any edit — to give CHANGE B clean evidence independent of the Phase-12 state.

### Concurrency

A **Phase-12 executor was live on this branch throughout** and interleaved six commits, including
`61c293a` between this plan's Task-4 and Task-5 commits and `8880c29`, which changed the suite
baseline between the BEFORE and AFTER runs. **No race was rebased or amended**; both are documented.

---

## 11. Commits

| Commit | What |
|---|---|
| `9fe9b02` | Tasks 2+3 — the ONE `consts.jl` edit (CHANGE A + CHANGE B), all nine guard sites, the six siblings, the tests, and the `P13_CONSTS_SHA` re-derivation at all three runner sites. **Landed BEFORE the re-run.** |
| `5441afb` | Task 4 — the corrected runner, the deleted redundancy arm, the de-staled ingestion header |
| `652692c` | Task 5 — the re-run artifact |
| `cb8c587` | Task 6 — `13-REPORT.md` (§8b, §8c, §9, **§13 corrected**, §14.6–§14.8) and the two superseding notes |
| `3257783` | Task 6 — `.planning/ROADMAP.md`, on its own pathspec |
| `156ca4e` | **The ruling** — `spike/test/runtests.jl`: the known-red ledger, the structural masking fix |
| `dfb5991` | **The ruling** — `13-REPORT.md`: the retraction (§8c), **Limit E** (§9), **§9a**, §14.9 |
| *(this commit)* | **The ruling** — `13-17-SUMMARY.md` flipped to `complete`, plus the additive `STATE.md` block |

All used `git commit … -- <explicit paths>`. Never `-a`, never `--no-verify`.

---

## Self-Check: PASSED

Re-run after the ruling, against disk and git rather than asserted:

- **Files exist:** `13-17-SUMMARY.md`, `13-REPORT.md`, `spike/test/runtests.jl`, `.planning/STATE.md`.
- **Commits exist:** `156ca4e`, `dfb5991` (this ruling); `151ad79`, `78dc37f` (the corroborating
  quick tasks, verified present rather than cited); `986d518` … `406015e` (the concurrent Phase-12
  window).
- **`spike/p13/consts.jl` git blob `d82caea64e19e8509f3abc1c09140cc7c2e2d497`** — the amendment
  blob, unchanged. `P13_ITERATION_ALLOWANCE = 1`, **UNSPENT**.
- **`git diff --quiet HEAD --` exits 0** for `src/`, `spike/Project.toml`, `spike/Manifest.toml`,
  `spike/validation/p12_consts.jl`, **all of `spike/p13/`**, and **both known-red test files** —
  so the five failing assertions are byte-identical to HEAD. The one `@test_skip` in
  `test_p13_real.jl:319` is a **pre-existing, unrelated** placeholder, present at HEAD, and is not
  one of the five.
- **`test_p12_suite.jl` still at `runtests.jl:174`, FIRST among the phase-test includes**, and
  `test_p12_consts.jl` testset 9 passes **44 / 44**.
- **The suite exits 1**, on exactly the five named-limit assertions. **That is the intended state.**
