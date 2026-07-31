# 12-15 RE-RUN — PRE-DECLARED INTERPRETATION

**COMMITTED BEFORE THE RUN LAUNCHES. The git timestamp of this file's commit precedes the run's
artifact, and that ordering is the point: a reading fixed after the numbers exist is not a
pre-declaration.**

Authorised by the user (via the orchestrator) as **ONE pre-declared larger run, with its own
baseline**. Not a search. If it comes back negative, that is the answer and it is reported as the
answer.

---

## 1. What is being changed, and what is NOT

| | value | changed? |
|---|---|---|
| pool size `P12_MINISPIKE_N` | 10,000 | **NO** |
| held-out block | 1,000 (pool head) | **NO** |
| arms | `(:car, :gp, :none)` | **NO** |
| head rank `K = P12_K_DEV` | 63 (full rank) | **NO** |
| `select_prior` | byte-identical | **NO** |
| every Tier-1 constant | — | **NO** |
| **epoch budget** | **18 → 100** | **YES — THE ONE AXIS** |
| **weight-init seeding** | unseeded → `p12_rng(P12_MINISPIKE_COUNTER + 2)` | **YES — auditability fix** |
| **per-epoch risk trace** | not persisted → persisted in bundle | **YES — auditability fix** |

**EPOCHS AND LR-SCHEDULE-LENGTH ARE ONE AXIS, NOT TWO.** NeuralEstimators anneals the schedule *over*
`epochs`, so they cannot move independently. That coupling is precisely the confound being removed —
in the 18-epoch run the LR annealed 33× to 1.51E-05 across exactly the budget, so the flat tail was
the schedule stopping, not the model saturating.

**THE POOL IS DELIBERATELY NOT CHANGED.** Train/val ratios of 1.007 / 1.036 / 1.038 at epoch 18 show
there is no overfitting, so data is **not** the binding constraint. Changing an axis the evidence does
not implicate would confound the axis it does. (Also: 50k × 100 epochs projects to **~220 min**, which
would **breach** `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150` — so the expensive option is not merely
declined, it is unavailable without a fresh user decision.)

**THE TWO SEEDING/TRACE CHANGES ARE AUDITABILITY FIXES TO RECORDED DEFECTS, NOT TUNING.** Neither
changes what the model learns; both change what can be checked afterwards. Seeding uses counter `+2`
off the **existing** reserved `P12_MINISPIKE_COUNTER` (alongside the current `+0` and `+1` mask
streams), so **no new Tier-1 constant and no pre-registration edit**.

### The design: a CONTROL and a TREATMENT, at the same seed

| run | epochs | seed | purpose | cost |
|---|---|---|---|---|
| **control** | 18 | `+2` | the same configuration as the existing result, but seeded | 10.56 min |
| **treatment** | 100 | `+2` | the one axis moved | 40.28 min |

**Why the control is not optional:** seeding fixes the init, so the existing unseeded artifact
**cannot** serve as the baseline. Comparing a seeded 100-epoch run against an unseeded 18-epoch one
would confound seeding with epochs — the exact error this run exists to avoid.

**Projected total 50.84 min against `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150` — NO BREACH** (33.9 %
of ceiling). Datagen is 0.00 min: pool directories are content-hash keyed on
`(n, arm, imsize_tag, shard_size)` and **epochs are not in the hash**. A breach would be a `BLOCKER`
and a fresh user decision, never absorbed by trimming the run.

---

## 2. THE BRANCHES, FIXED BEFORE ANY NUMBER EXISTS

### 2.0 A necessary precision: `select_prior` does NOT read `skill`

**Verified by execution, not by reading** (`select_prior` returns `CAR` on a score table that has no
`skill` field at all). The rule's inputs are `coverage`, `rmse` and `n_test` only, and
`test_p12_train.jl` asserts that `skill` and `trivial` never appear in its body.

**Consequence, stated so it cannot be read as a moved bar:** the rule **can** return `CAR` or `GP`
with negative skill — it nearly did on the invalid run. The skill condition in §2.1 is an
**interpretive layer declared in advance, sitting on top of an unchanged rule**. It adds no bar to
`select_prior`, changes no threshold, and is not a gate. If the rule names a winner and that winner
has skill ≤ 0, **the rule's output is reported verbatim AND the premise is recorded as unsupported** —
both, without touching the rule.

### 2.1 SUPPORTS the spatial premise

**All three must hold:**
1. `select_prior` returns `CAR` or `GP` (so the arm is admissible, `beats_ablation = true`); **and**
2. that arm's **`skill > 0`** — it beats the constant-zero predictor; **and**
3. its RMSE beats `:none` by **≥ 1 %** (see §2.3 for < 1 %).

Only then does the CAR-vs-GP comparison become meaningful **for the first time**.

### 2.2 LEAVES IT UNSUPPORTED at the larger scale

Any of:
- `select_prior` returns `NONE-BEATS-ABLATION`; **or**
- it names a winner whose `skill ≤ 0`.

**This is the answer. It is reported as the answer, and nothing re-runs.**

**Edge cases closed against the code, by execution rather than reasoning:**
- **Neither spatial arm admissible** → `adm` is empty → `NONE-BEATS-ABLATION`, `beats_ablation =
  false`, **even when both arms crush the ablation on RMSE** (verified: coverage 0.99/0.60 with RMSE
  0.10/0.05 against an ablation at 0.90 still returns `NONE-BEATS-ABLATION`). Rule 1 is a hard gate.
- **An INADMISSIBLE arm with positive skill that beats the ablation** → it is never placed in `adm`,
  so it can **never** be selected, whatever its skill or RMSE (verified: an inadmissible GP at RMSE
  0.05 loses to an admissible CAR at 0.95, and the branch is `NONE-BEATS-ABLATION` because CAR then
  fails rule 3). **A miscalibrated arm is disqualified, not "better" — no branch is left to be
  inferred after the numbers exist.**

### 2.3 AMBIGUITY, PRE-COMMITTED

**(a) Skill turns positive but no spatial arm beats the ablation.**
**NOT ambiguous, and the most informative outcome available: the method works and the spatial prior
adds nothing.** Reported as a **clean negative on the premise**. This is the outcome most open to
being spun as "promising", and pre-declaring it as a clean negative is what stops that.

**(b) A spatial arm beats the ablation by < 1 %.**
**Declared a TIE; the comparison is reported UNRESOLVED.** One seeded run cannot measure init
variance, so a sub-1 % margin is not distinguishable from luck. **1 % is `select_prior`'s existing
tie-break threshold, reused rather than invented — it is not a new bar.**
**I will NOT run further seeds to resolve it. That is a search, and the user authorised one run.**

**(c) Train/val separates (overfitting appears) with skill still ≤ 0.**
That is evidence for a **new** hypothesis — that 9,000 training samples is the binding constraint.
**REPORTED, NOT ACTED ON.** It requires a fresh user decision and must not be folded into this run's
conclusion.

---

## 3. THE INIT-VARIANCE SMELL TEST — free information, and its limit stated up front

The seeded 18-epoch control is directly comparable to the existing **unseeded** 18-epoch numbers.
This will be reported:

- **If they land close** — single-run results are more trustworthy than assumed.
- **If they diverge substantially** — direct evidence that no single run settles a 0.33 % margin,
  which retro-justifies refusing to read a ranking out of the existing run.

**LIMIT, DECLARED IN ADVANCE: this is n = 1 versus n = 1. It bounds nothing.** It is a smell test, not
a variance estimate, and it will be reported with that caveat attached. It will **not** be used to
license or refuse any conclusion.

---

## 4. WHAT DOES NOT DRIFT, WHATEVER THE RUN SAYS

- **The qualifier "at mini-spike scale" stays on the existing result.**
- **Verbatim, into the verdict: _a comparison between two failed arms is not evidence about priors —
  the CAR-vs-GP question is UNRESOLVED, not answered._** The new run will **not** be presented as
  vindicating or refuting a comparison the existing run never made.
- **The two negative results stay reported as two**: GP calibrated-but-uninformative (coverage
  0.89666, skill −0.069) versus CAR mis-centred AND over-wide (coverage 0.99981, `c0_rmse` **5.206**
  against the ablation's 0.798). Collapsing them loses the more diagnostic half.
- **No Tier-2 constant is appended unless the branch that fired licenses it.** On
  `NONE-BEATS-ABLATION`, `12-15-PLAN:265` forbids both blocks and all three sentinels stay.
- **One run.** If §2.2 fires, it is reported and the phase stops for the user.
