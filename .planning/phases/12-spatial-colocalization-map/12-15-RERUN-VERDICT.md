# 12-15 RE-RUN — VERDICT

**Every number below is re-derived by loading the JLD2 artifacts, never copied from stdout.** The
training traces are the exception and are read from the run logs, which is stated at their table.

Reading fixed in advance by `12-15-RERUN-PREDECLARATION.md`, committed `ea1d2b7` at
**2026-07-31T15:37:40+02:00**. The artifacts were generated at **16:59:48Z** (control) and
**17:29:28Z** (treatment). **The pre-declaration's commit precedes both, and that ordering is the
point.**

| artifact | epochs | generated | in repo as |
|---|---|---|---|
| unseeded baseline | 18 | 2026-07-31T12:57:29Z | `spike/validation/p12_minispike_report.jld2` (unchanged, `2fd75352…`) |
| **control** | 18 | 2026-07-31T16:59:48Z | `spike/validation/p12_minispike_report.rerun_control_18ep.jld2` |
| **treatment** | 100 | 2026-07-31T17:29:28Z | `spike/validation/p12_minispike_report.rerun_treatment_100ep.jld2` |

---

## 1. THE BRANCH THAT FIRED

**`select_prior` returned `NONE-BEATS-ABLATION` in the control AND in the treatment.** It returned
`NONE-BEATS-ABLATION` in the unseeded baseline too. **Three runs, three times the same branch.**

Stated as the branch, per the pre-declaration §2.2, first bullet:

> *`select_prior` returns `NONE-BEATS-ABLATION` … **This is the answer. It is reported as the answer,
> and nothing re-runs.***

**§2.1 did not fire** (no arm was admissible, so no arm could be selected).
**§2.3(a) did not fire** on the treatment — it requires skill to turn positive, and at 100 epochs
every arm has skill ≤ 0.
**§2.3(b) did not fire** — rule 1 returned before rule 3, so no RMSE margin was ever compared and no
TIE was declared.
**§2.3(c) DID fire.** See §5.

### 1.1 The rule returned on a DIFFERENT rule than last time

This is the one substantive change in the branch trace and it must not be glossed:

| run | returning rule | why |
|---|---|---|
| unseeded 18 ep | **rule 3** | `:gp` was admissible but did not beat the ablation |
| **control 18 ep** | **rule 1** | **BOTH** spatial arms miscalibrated → `adm` empty |
| **treatment 100 ep** | **rule 1** | **BOTH** spatial arms miscalibrated → `adm` empty |

Recorded `selection_reason`, treatment, verbatim from the artifact: *"BOTH spatial arms are
miscalibrated: coverage 0.7982 (car) and 0.8003 (gp) against 0.9 +/- 0.03. A miscalibrated arm is
DISQUALIFIED, not 'better' (rule 1)."*

### 1.2 BOTH pre-declared edge cases fired for real

§2.2 closed two edge cases **by execution, in advance**, so that no branch was left to be inferred
after the numbers existed. **Both of them then happened.** This is the strongest available evidence
that the closure was worth doing, and it is why nothing here required a judgement call:

**Edge case 1 — "neither spatial arm admissible → `NONE-BEATS-ABLATION`, even when both arms crush
the ablation on RMSE."** In the treatment, **both** spatial arms beat the ablation on RMSE — `:car`
0.99811 and `:gp` 1.10068 against `:none` 1.14663, i.e. by **12.95 %** and **4.01 %** — and both are
disqualified on coverage regardless.

**Edge case 2 — "an INADMISSIBLE arm with positive skill that beats the ablation can never be
selected, whatever its skill or RMSE."** In the control, `:car` has **positive skill (+0.08629)**
**and** beats the ablation on RMSE (0.90850 against 1.02352, by 10.28 %) — **and is still
disqualified**, coverage 0.99734 against the band [0.87, 0.93].

**A miscalibrated arm is disqualified, not "better."** That is the rule as written before any number
existed, and it is doing exactly the work it was written to do.

---

## 2. THE NUMBERS

Band: `P12_COVERAGE_NOMINAL = 0.90` ± `P12_STAGE2_COVERAGE_TOST_DELTA = 0.03` → **[0.87, 0.93]**,
asserted at `n_test = 1000 ≥ P12_STAGE2_N_MIN = 271` in every run.

### CONTROL — 18 epochs, seeded

| arm | RMSE | trivial | **skill** | coverage | admissible | `c0_rmse` |
|---|---|---|---|---|---|---|
| `:car`  | 0.90850 | 0.99430 | **+0.08629** | 0.99734 | **no** | 1.94103 |
| `:gp`   | 1.60322 | 0.99331 | **−0.61402** | 0.99998 | **no** | 6.03903 |
| `:none` | 1.02352 | 0.99634 | **−0.02728** | 0.97544 | — | 0.82373 |

`admissible = Symbol[]` — **empty**. `beats_ablation = false`.

### TREATMENT — 100 epochs, same seed

| arm | RMSE | trivial | **skill** | coverage | admissible | `c0_rmse` |
|---|---|---|---|---|---|---|
| `:car`  | 0.99811 | 0.99430 | **−0.00383** | 0.79816 | **no** | 1.21123 |
| `:gp`   | 1.10068 | 0.99331 | **−0.10810** | 0.80028 | **no** | 1.14184 |
| `:none` | 1.14663 | 0.99634 | **−0.15084** | 0.82403 | — | 0.81130 |

`admissible = Symbol[]` — **empty**. `beats_ablation = false`.

---

## 3. THE ONE AXIS: what 82 more epochs bought

**The seeding worked, and it is verifiable rather than asserted.** Control and treatment have
**byte-identical epoch-1 traces for every arm** (`:car` initial validation risk 421798.8, epoch 1
train 1328.974 / val 145.889 — identical in both runs). **The two runs start from literally the same
weights.** The epoch budget is therefore the only thing that differs, which is what the control
exists to establish.

### 3.1 EARLY STOPPING FIRED — the confound the run was built to remove is GONE

This is the decisive operational result. Read from the run logs, not the artifacts, because the
bundles are not persisted (`out_path = nothing`); the per-epoch trace is now persisted **inside the
bundle** by `faa1871`, but the mini-spike does not save its bundles to disk.

| arm | control: epochs run | treatment: epochs run | final val risk 18 → 100 | change |
|---|---|---|---|---|
| `:car`  | 18 of 18 | **56 of 100** | 109.915 → 91.731 | **−16.54 %** |
| `:gp`   | 18 of 18 | **99 of 100** | 119.920 → 88.145 | **−26.50 %** |
| `:none` | 18 of 18 | **36 of 100** | 108.361 → 96.282 | **−11.15 %** |

**At 18 epochs all three arms ran the full budget — the run ended on its budget, not on convergence.
At 100 epochs two of three arms stopped EARLY, at 56 and 36.** `12-MINISPIKE-VERDICT.md:112` recorded
that the experiment as designed *"CONFOUNDS 'converged' with 'annealed to a stop'"*. **That confound
is now removed by measurement.**

**So the under-training diagnosis was RIGHT, and it does not rescue the premise.** The 18-epoch run
genuinely was under-trained; training to convergence genuinely does reduce validation risk, by 11–27 %;
and **no arm becomes admissible and no arm reaches positive skill.** The leading alternative
explanation for the original negative has been tested and eliminated. **That makes this negative
stronger than the one it replaces, not weaker.**

### 3.2 Skill did not improve with training — it moved in both directions

| arm | skill @ 18 ep | skill @ 100 ep | Δ |
|---|---|---|---|
| `:car`  | **+0.08629** | −0.00383 | **−0.09012 (worse)** |
| `:gp`   | −0.61402 | −0.10810 | +0.50592 (better) |
| `:none` | −0.02728 | −0.15084 | −0.12356 (worse) |

**Two of three arms got worse.** The one arm that had positive skill lost it. More training is not a
monotone improvement in skill on this task, and at the converged budget **no arm beats a constant-zero
predictor.**

### 3.3 Coverage crossed the band without landing in it — REPORTED, and deliberately NOT acted on

| run | `:car` | `:gp` | `:none` | relative to band [0.87, 0.93] |
|---|---|---|---|---|
| control, 18 ep | 0.99734 | 0.99998 | 0.97544 | **all three OVER-cover** |
| treatment, 100 ep | 0.79816 | 0.80028 | 0.82403 | **all three UNDER-cover** |

Longer training makes the posteriors monotonically narrower, and between 18 and 100 epochs **every
arm crosses the admissible band from above to below without stopping inside it.** On this evidence
admissibility tracks *training duration*, not the prior.

**This is reported as an observation and nothing is done with it, for two independent reasons, and
either alone would be sufficient.** First, searching for the epoch count that lands coverage inside
the band is **a search**, and the user authorised **one run**. Second, and worse, it would be
**tuning a model until it passes its own calibration gate** — which does not produce a calibrated
model, it produces a gate that has stopped measuring anything. That is precisely the failure this
phase has already recorded three families of. Acting on it requires a fresh user decision.

---

## 4. THE INIT-VARIANCE SMELL TEST (§3) — THEY DIVERGE, AND SUBSTANTIALLY

Seeded 18-epoch control against the preserved **unseeded** 18-epoch baseline. Same file, same pool,
same 18 epochs, same held-out block; **only the weight init differs.**

| arm | RMSE unseeded | RMSE seeded | Δ | **relative** | skill unseeded | skill seeded |
|---|---|---|---|---|---|---|
| `:car`  | 1.43483 | 0.90850 | −0.52632 | **−36.68 %** | −0.44305 | **+0.08629** |
| `:gp`   | 1.06141 | 1.60322 | +0.54182 | **+51.05 %** | −0.06856 | −0.61402 |
| `:none` | 1.05793 | 1.02352 | −0.03441 | −3.25 % | −0.06181 | −0.02728 |

Initial validation risks differ too (`:car` 852,488.1 → 421,798.8; `:gp` 1,306,758.8 → 620,645.4;
`:none` 97,665.2 → 59,427.9), confirming a genuinely different draw rather than a scoring difference.

**The two spatial arms swap which one is worse.** `:car` goes from worst arm to best-and-positive-skill;
`:gp` goes from the only admissible arm to the worst arm in any run.

**Read against the margin it has to be read against: the gap that decided rule 3 in the unseeded run
was 0.33 %. The init spread measured here is 37–51 % — two orders of magnitude larger.**

Per the pre-declaration §3, the "diverge substantially" branch:

> *If they diverge substantially — direct evidence that no single run settles a 0.33 % margin, which
> retro-justifies refusing to read a ranking out of the existing run.*

**It fires.** `12-MINISPIKE-VERDICT.md:62-76` argued from the *mechanism* (unseeded init) that the
0.33 % margin was not a stable quantity and declined to read a ranking out of it. **That refusal is
now supported by direct observation and not only by an argument.**

**THE LIMIT, CARRIED AS DECLARED IN ADVANCE: this is n = 1 versus n = 1. IT BOUNDS NOTHING.** It is a
smell test, not a variance estimate. It licenses no conclusion and refuses none, and it is not used
below to support anything.

**What is nonetheless robust across all three runs: the branch.** Every per-arm number moved, some by
half their own size; `NONE-BEATS-ABLATION` did not.

---

## 5. §2.3(c) FIRES — REPORTED, NOT ACTED ON

> *Train/val separates (overfitting appears) with skill still ≤ 0. That is evidence for a **new**
> hypothesis — that 9,000 training samples is the binding constraint. **REPORTED, NOT ACTED ON.** It
> requires a fresh user decision and must not be folded into this run's conclusion.*

| arm | val/train @ 18 ep | val/train @ 100 ep |
|---|---|---|
| `:car`  | 1.0168 | **1.1448** |
| `:gp`   | 1.0054 | **1.0742** |
| `:none` | 1.0291 | **1.1581** |

Train and validation **have begun to separate** where at 18 epochs they were essentially on top of
each other (1.005–1.029), and **skill is ≤ 0 for every arm at 100 epochs.** The condition is met
exactly as written.

**This is not folded into the conclusion below, and no pool was changed.** The pre-declaration fixed
the pool at 10,000 in advance and recorded that 50k × 100 epochs projects to ~220 min, which would
**breach** `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150`. **The larger-pool experiment is unavailable
without a fresh user decision, not merely declined.**

---

## 6. BUDGETS — EACH AGAINST ITS OWN CONSTANT, BY NAME

**Two separate budgets that happen to carry the same number.** Never summed, never crossed. The
mini-spike ceiling is enforced **per run** by `_p12ms_check_budget`; the datagen ceiling is enforced
**per `generate_p12_pool` call**.

| run | `elapsed_min` | against | % of ceiling | BLOCKER? |
|---|---|---|---|---|
| control | **11.29** | `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150` | 7.5 % | no |
| treatment | **29.29** | `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150` | 19.5 % | no |

**`datagen_min` is PER ARM**, each pool judged independently inside `generate_p12_pool` against
**`P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150`**: `:car` 0.09, `:gp` 0.00, `:none` 0.00 — in **both**
runs. All near zero because the pools were already complete and resume-by-skip cost nothing, exactly
as the pre-declaration predicted (epochs are not in the pool's content hash). **Their sum is not a
quantity anything compares to 150.**

Projection check: the pre-declaration projected 10.56 + 40.28 = 50.84 min. **Actual 11.29 + 29.29 =
40.58 min.** The treatment came in **under** projection because early stopping ended two arms before
the budget — the projection assumed all 100 epochs would run. No ceiling was approached in either run.

---

## 7. DOES THE CONCLUSION CHANGE?

**No. The headline conclusion is unchanged, and it is now better supported.**

**`SELECTED PRIOR: NONE-BEATS-ABLATION`. The spatial premise is UNSUPPORTED at mini-spike scale.**

What is **strengthened**:

- The **under-training explanation is tested and eliminated.** Early stopping fires at 56/99/36 of
  100; validation risk improves 11–27 %; no arm becomes admissible; no arm reaches positive skill.
  The original run could not separate "under-trained" from "cannot learn this". **This one can, and
  the answer is not under-training.**

What is **corrected**, and this is a real correction to the previous verdict's framing:

- **`12-MINISPIKE-VERDICT.md:28` called it "the stronger finding" that no arm learned the field at
  all — every arm worse than a constant-zero predictor. THAT CLAIM DOES NOT SURVIVE A RE-SEED.** At
  the same 18 epochs with a different init, `:car` reaches skill **+0.08629**. "Every arm has negative
  skill" is an artifact of one initialisation, not a property of the method.
- **What IS robust across all three runs is narrower and should replace it: no spatial arm is ever
  both calibrated and better than the ablation.** In three runs, the set of admissible arms was
  `[:gp]`, then `[]`, then `[]`, and `beats_ablation` was `false` every time.
- The **0.33 % margin** is now demonstrably noise against a 37–51 % init spread. No ranking among the
  three arms is established by any of these runs.

What does **not** drift, carried verbatim as required:

- **The qualifier "at mini-spike scale" stays on the result.**
- > **A comparison between two failed arms is not evidence about priors — the CAR-vs-GP question is
  > UNRESOLVED, not answered.**

  This re-run **does not** vindicate or refute a comparison the earlier run never made. In the control
  and the treatment **neither** arm was even admissible, so the two runs added here compare two arms
  that both failed rule 1 — which is, if anything, further from a prior comparison than the original.
- **The two negatives stay reported as two, and they are the ORIGINAL run's two** — GP
  calibrated-but-uninformative (coverage 0.89666, skill −0.069) versus CAR mis-centred AND over-wide
  (coverage 0.99981, `c0_rmse` 5.206 against the ablation's 0.798). Collapsing them loses the more
  diagnostic half. **In the new runs the pairing does not recur in that form** — at 18 epochs seeded,
  both spatial arms are over-wide (0.997, 1.000); at 100 epochs both are under-wide (0.798, 0.800).
  That is a third and a fourth pattern, not a reproduction of the first two.
- **`select_prior` does not read `skill` at all.** It did not need the §2.0 disposition this time: the
  rule returned `NONE-BEATS-ABLATION` rather than naming a winner, so there was no winner whose skill
  had to be checked. **The rule was not touched.**

---

## 8. WHAT WAS NOT DONE

- **No Tier-2 constant appended.** `spike/validation/p12_consts.jl` verified **byte-unchanged**
  (`git diff --quiet HEAD` clean). `P12_CHOSEN_PRIOR`, `P12_K_PROD`, `P12_D_PROD` and `P12_N_LOW` all
  **absent** — all three sentinels intact. The branch that fired forbids the append.
- **No further seeds.** §2.3(b) never fired, but the standing instruction holds regardless: **one run,
  not a search.**
- **`PRODUCTION RANK: K = 35`** is reported by both new runs (against 24 in the unseeded run) and is
  **not** appended, for the same reason as before — it is scaled against `winner_rmse`, which is a
  posterior RMSE from a net that is not admissible. The **truncation curve itself is BYTE-IDENTICAL
  across all three runs** — it is arithmetic on the drawn field ensemble and touches no net. That is
  the control confirming nothing structural changed between runs; only the nets did.
- **Wave 9 (12-17) not started.** Still impossible as written: it generates its pool with
  `arm = P12_CHOSEN_PRIOR`, which does not exist on this branch.
- **`.planning/STATE.md` not modified** — the orchestrator owns it (`12-EXECUTOR-BRIEFING.md:32`).
- `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `artifacts/` all verified clean;
  `spike/data/cache/p11/` untouched.
- **`spike/validation/p12_minispike_report.jld2` restored byte-exactly** to its committed blob
  `2fd75352…` after the runner overwrote it (the report path is a fixed constant, so the treatment
  clobbered it). Restored by copying the preserved artifact, **not** by `git checkout --`, and
  verified by `git hash-object`. The two new runs are committed under their own names so that no
  existing citation silently changes meaning.
