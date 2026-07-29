# 13 — PRE-REGISTRATION AMENDMENT to the D-15 real-image arm: the channel pair

**Date:** 2026-07-29
**Status:** **PROPOSED.** Frozen on approval of §7 (the one open question is *where the amended
constant lives*, not *what it is*). Every value in §5 is settled by the user's ruling and is not
reopened by §7.
**Implements:** the user ruling of 2026-07-29 (`f411ee7`), restating and completing the factual
correction first recorded 2026-07-25 (`ea4a7d3`)
**Amends:** `spike/p13/consts.jl` section I2 — `P13_REAL_CHANNEL_PAIR`,
`P13_REAL_REDUNDANCY_PAIR`, `P13_REAL_ANCHOR_MBAR`, and the two prose constants that quote them
**Supplements — does NOT replace:** `13-SC2-AMENDMENT.md` §5 (the SC3 three-arm scope) and
`.planning/ROADMAP.md` Phase-13 SC3
**Executed by:** plan `13-17-PLAN.md`

---

## 0. Provenance disclosure — read this first, because the ordering is against us

**This amendment lands AFTER the results it affects exist. `13-SC2-AMENDMENT.md` could open by
saying "no Phase-13 result exists"; this one cannot.** Plan 13-16 ran on 2026-07-29 at 17:07
(`2112bed`), the phase report was written the same day (`0c1a0bd`), and this document is dated
after both. That is, on its face, the exact shape of the pattern `07-GATE-AMENDMENT.md` §6.4 named
and that `P13_ITERATION_ALLOWANCE = 1` exists to prevent. Saying so first is the only honest way to
open.

The claim of this document is **not** that the ordering is innocent. It is that the *act* is a
different act, and that the difference is checkable rather than asserted. Four facts carry it, all
verified against the working tree at `d7195cc4c5ba5a4aaceef8f2d59ddf9bd4422be1`:

1. **The constant names the wrong physical object** (§1). `test/runtests.jl:105` records the
   fixture channels as `["blue", "green", "red"]`, so channel 1 is the DAPI/Hoechst nuclear
   counterstain. A colocalization arm pre-registered on the counterstain measures
   counterstain-versus-protein overlap, which is not colocalization. Correcting a factual
   mis-designation is not the same act as relaxing a bar that proved inconvenient.
2. **The upstream source was already corrected, and Phase 13 drifted from it** (§2). This is the
   strongest fact and it is the one to lead with. `spike/simulator/ghat.jl:58-68` and
   `spike/simulator/calibration.jl:351,433` already mark `0.3292 / 0.2481` **SUPERSEDED** and
   already record the corrected unmasked `0.4603 / 0.3815` — done on 2026-07-25, **thirty-nine
   hours before `spike/p13/consts.jl` was written**. Phase 13's pre-registration never followed its
   own upstream calibration source. This is drift being repaired, not a result being chased.
3. **No success criterion moves** (§3). No AUC floor, no ECE band, no `P13_ITERATION_ALLOWANCE`, no
   gating threshold of any kind is touched. The phase's only gating arm (13-12, simulated) does not
   read `P13_REAL_CHANNEL_PAIR` at all — grep-verified, §3.
4. **The amendment picks the WORSE number, and a better one was available and unforbidden**
   (§2.1, §4). Unmasked, the two fixtures barely separate: `0.4603` against `0.3815` — a separation
   of **0.0788**, narrower even than the **0.0811** of the superseded pair it replaces. The masked
   read of the *same corrected pair* separates them by **0.8576** (`+0.8238 / −0.0338`), and
   `ghat.jl:69-71` — the frozen upstream this amendment adopts — says in its own words that the
   masked pair "separates the controls far more sharply", while `:79-83` **retracts** the earlier
   "MASKED is biased NEGATIVE" framing as having OVERSTATED the case. **Upstream therefore neither
   forbids the masked read nor prefers the unmasked one.** This amendment forbids it anyway (§5.3),
   on a property of the code rather than a preference, and takes the ten-times-narrower number — and
   then deletes an arm on top of it (§5.2). **An amendment that selects against its own arm's
   interest, when the better number was there for the taking, is structurally not a result being
   chased.** §2.1 states that check in mechanical form.

### 0.1 The authorisation timeline, stated without softening and without overstatement

| When (verified) | Commit | What |
|---|---|---|
| 2026-07-25 22:54 | `151ad79` | `ghat.jl` real-anchor header corrected to record MASKED as well as unmasked reads |
| 2026-07-25 23:28 | `78dc37f` | `ghat.jl` header corrected to the **c2/c3** pair; `0.3292 / 0.2481` marked **SUPERSEDED**; unmasked `0.4603 / 0.3815` recorded |
| 2026-07-25 23:56 | `ea4a7d3` | STATE.md blocker recorded: the Phase-13 plans "pre-register the real-image arm on the DAPI counterstain and **MUST be corrected before Phase 13 executes**"; `11-10-PLAN.md` fixed in the same commit |
| **2026-07-27 14:43** | **`c42cc8e`** | **Phase-13 Tier-1 `consts.jl` LOCKED — still `(1,2)` / `(1,3)` / `0.3292` / `0.2481`. **38 h 46 min 18 s** after the blocker (23:56:54 → 14:43:12, both copied from `git log`).** |
| 2026-07-27 21:46 | `bfba6ac` | Tier-2 `P13_TAU` appended (the only other commit to `consts.jl`) |
| 2026-07-29 17:07 | `2112bed` | **13-16 executed** on `(1,2)` primary and `(1,3)` redundancy |
| 2026-07-29 17:48 | `f411ee7` | The miss recorded; the user's **full disposition** received: pair `(2,3)`, redundancy dropped, unmasked values, document as an amendment |
| 2026-07-29 18:09 | `d7195cc` | HEAD at the time of writing this document |

Read that table exactly as it reads, in two parts:

- **The factual correction predates everything.** Which pair is the colocalization pair was settled
  and committed on 2026-07-25, **before `consts.jl` existed** and **3 days 17 hours before 13-16
  ran**. Nothing about *which channels are the target proteins* was learned from a Phase-13 result.
- **The specific disposition postdates the run by 41 minutes.** Dropping the redundancy arm,
  choosing unmasked over masked, and handling it as an amendment rather than a named limit were
  ruled at 17:48 on 2026-07-29, after 13-16 finished at 17:07. That is not obscured here. What was
  ruled after the run is *how to dispose of an already-known defect*, not *whether the defect
  exists*.

Neither half of that is improved by hiding the other.

### 0.2 The executors did not err, and that matters for reading this document

13-15 and 13-16 both detected the problem independently and both **declined to fix it**.
`13-16-SUMMARY.md:276-279` names `c1` as the counterstain and records that choosing a new channel
pair after the fixtures had been measured would be the wrong act, so it named the limit instead of
editing a frozen constant. `spike/p13/run_p13_realimage.jl:49-57` says the same thing in the
runner's own banner. **That was the correct call under the pre-registration as it stood**, and it is
the reason the correction requires a document rather than a patch: the executors' refusal is what
makes an authorised, dated, argued amendment the only legitimate route.

### 0.3 The precedent's test, applied

`13-SC2-AMENDMENT.md` §0 quotes the legitimacy test verbatim:

> *Would this change have been made, with this justification, by someone who had seen the gate's
> code but none of its outputs?*

**Yes — demonstrably, because someone did.** The 2026-07-25 blocker (`ea4a7d3`) was written before
any Phase-13 code existed, from `test/runtests.jl:105` and the corrected `ghat.jl` alone, and it
specifies the same pair and the same two unmasked numbers this document adopts. The identical
correction was applied to `11-10-PLAN.md` in that same commit and has stood since. Phase 13's
failure was not of judgement but of application: the ruling was on the record and was not carried
into `consts.jl`.

**Do not file this alongside `07-GATE-AMENDMENT.md` as "the third time they moved a goalpost."** No
goalpost is in this document. §3 is the greppable form of that claim.

---

## 1. Fact 1 — the constant names the wrong physical object

### The defect

`spike/p13/consts.jl:552-553` freeze

```
const P13_REAL_CHANNEL_PAIR    = (1, 2)
const P13_REAL_REDUNDANCY_PAIR = (1, 3)
```

`test/runtests.jl:105` records the fixture channel map:

```julia
channels = ["blue", "green", "red"]
```

So `c1` is the **blue channel — the DAPI/Hoechst nuclear counterstain**. It is not a target
protein and it was never intended to be one: a counterstain marks nuclei so the operator can find
cells. Both pre-registered pairs therefore contain it, and the primary arm as executed measured
**counterstain-versus-green-protein overlap**.

That quantity is not colocalization under any definition this project uses. It is dominated by the
trivial fact that both channels are bright where there is a cell and dark where there is not — the
shared-background covariance that `ghat.jl:72-75` names explicitly as the unmasked estimator's
positive bias. Calling it "colocalization of two targets" is a category error about the physical
object, not an optimistic reading of a number.

### Why this is a correction and not a relaxation

A relaxation changes *how well* the measured object must score. This changes *which object is
measured*. The test that separates them is mechanical: after this amendment, is there any quantity
that previously had to clear a bar and now has to clear a lower one? **No** — because no real-image
quantity has ever had a bar at all (`P13_REAL_IS_GATED = false`; `13-SC2-AMENDMENT.md` §5(g), stated
in greppable words: *No pass/fail threshold is defined for any real-image quantity*). See §3.

---

## 2. Fact 2 — the upstream source was corrected first, and Phase 13 drifted from it

**This is the strongest fact in the document.** The correction Phase 13 is applying is not new
information. It is Phase 13 catching up with its own frozen calibration source, which was corrected
two days before Phase 13's pre-registration was written.

`spike/simulator/ghat.jl:58-68` [VERIFIED at HEAD]:

```
#   real anchor (D-16): faithful LoadImages.jl load_tiff (NOT luminance), measured BOTH ways
#     channel pair     : c2/c3 = green/red  (ANCHOR_CHANNEL_PAIR = (2, 3))
#     channel map (test/runtests.jl:105): c1 = blue (DAPI/Hoechst nuclear counterstain), c2 = green, c3 = red
#     why the pair matters: c1 is a nuclear COUNTERSTAIN, not a target protein, so any pair
#       involving it measures counterstain-vs-protein overlap and NOT colocalization. The earlier
#       figures (positive 0.3292 / negative 0.2481, unmasked c1/c2) were measured against the
#       counterstain and are SUPERSEDED, not merely updated.
#     UNMASKED (what patch_summary / the v2.0 training + inference path sees -- no Otsu mask):
#       positive μ = 0.4603 (64/64 patches), negative μ = 0.3815 (64/64 patches)
#     MASKED (the v1.0 analysis pipeline: per-channel Otsu _apply_mask! then _exclude_zero):
#       positive μ = 0.8238 (22/64 patches), negative μ = -0.0338 (20/64 patches)
```

`spike/simulator/calibration.jl:433` [VERIFIED at HEAD] emits the same ruling from the generator,
so a regeneration stays a no-op on it:

> `why: c1 is a nuclear COUNTERSTAIN, not a target protein -- any pair involving it measures
> counterstain-vs-protein overlap, NOT colocalization; the earlier unmasked c1/c2 figures
> (positive 0.3292 / negative 0.2481) are SUPERSEDED`

Both landed in `78dc37f` on **2026-07-25 23:28**. `spike/p13/consts.jl` was first committed in
`c42cc8e` on **2026-07-27 14:43** — and it quotes `positive = 0.3292, negative = 0.2481` at line
576 as "the frozen lineage anchors recorded at `spike/simulator/ghat.jl:39`", citing a line whose
content had already been superseded at the moment of citation.

### The consequence for how this amendment must be read

Phase 13 is not choosing new numbers after seeing a disappointing result. It is adopting numbers
that were already frozen upstream, already published in the simulator's own calibration header,
already applied to a sibling phase's plan (`11-10-PLAN.md`, `ea4a7d3`), and already recorded as a
blocker against Phase 13 specifically. The `consts.jl` values are the outlier, not the amendment.

*Would a blind reviewer have made this change?* **Yes, and one did, twice, before Phase 13 existed
— once in `ghat.jl`, once in `11-10-PLAN.md`.**

### 2.1 The check that separates a correction from a relaxation: this amendment picks the WORSE number

The strongest test of whether an amendment is a correction or a rescue is not what it argues. It is
whether, at the one point where a **better number was available and permitted**, it took it.
**It did not — and the record shows the better number was there for the taking.**

`spike/simulator/ghat.jl:69-71` [VERIFIED at HEAD] — the same frozen upstream header this amendment
adopts in §2 — states, unprompted:

```
#   anchor caveat    : on this pair the two estimators AGREE in order -- positive > negative under
#     BOTH -- and the MASKED pair separates the controls far more sharply (0.8238 vs
#     -0.0338) than the unmasked pair does (0.4603 vs 0.3815).
```

and `:79-83` goes further, **retracting the one argument that had previously been used against
masking**:

```
#     The earlier "MASKED is biased NEGATIVE" framing
#       OVERSTATED it: here the masked POSITIVE control reads 0.8238, so the selection does
#       not prevent detecting strong colocalization -- the large negative readings that motivated
#       that wording were predominantly the wrong-channel (counterstain) artifact, not the mask.
```

Read together, the upstream position is unambiguous: **it neither forbids the masked read nor
prefers the unmasked one.** It records both on the corrected pair, rates the masked pair the sharper
separator, and withdraws its own prior objection to masking. Meanwhile `:72-75` says the *unmasked*
estimator is biased POSITIVE and that on this pair the bias "COMPRESSES the contrast between the
controls, which is why the unmasked pair barely separates them."

So at the moment of choosing, the options were:

| Candidate anchor | separation | upstream's own verdict on it | adopted? |
|---|---|---|---|
| **masked** c2/c3 `+0.8238 / −0.0338` | **0.8576** | *"separates the controls far more sharply"*; the anti-masked framing **retracted as overstated** | **NO — forbidden by §5.3** |
| **unmasked** c2/c3 `+0.4603 / +0.3815` | **0.0788** | biased POSITIVE, *"COMPRESSES the contrast … barely separates them"* | **YES** |
| unmasked c1/c2 `+0.3292 / +0.2481` | 0.0811 | measured against the counterstain — **SUPERSEDED** | no (§1) |

**This amendment took the 0.0788.** It adopted a separation roughly **eleven times narrower** than
the one it was permitted to adopt, and narrower even than the superseded number it replaces. The
reason is a property of the code and not a preference: `patch_summary` (`spike/contract.jl:85`)
applies **no Otsu mask**, so the net was trained on unmasked summaries and masked summaries would be
off-distribution input to it (§5.3 in full). Then, on top of that, the amendment **deletes one of
its own two arms** (§5.2).

An amendment that is chasing a result does not do this. Handed a permitted, upstream-endorsed,
order-of-magnitude-sharper number, it takes it. The mechanical form of the check, for a reader who
wants a test rather than a paragraph:

> **Was a better number available, unforbidden by any prior source, at the moment of amendment?
> Yes — 0.8576 against 0.0788, on the very pair being adopted.
> Was it taken? No.**

That is the difference between a correction and a relaxation, and this document would rather have it
on the page than leave it implicit. It is also falsifiable: if a later reader finds the masked
figures adopted anywhere in this arm, this paragraph is a lie and §5.3's prohibition was theatre.
`13-17-PLAN.md` makes that a grepped acceptance criterion rather than a promise.

---

## 3. Fact 3 — no gating threshold moves, and the gating arm cannot see this change

### The greppable claim

**This amendment changes no floor, no band, no tolerance, no allowance and no bar of any kind.**
Specifically untouched, by name:

`P13_AUC_FLOOR_COLOC`, `P13_AUC_FLOOR_EXCLUSION`, `P13_GATE_M`, `P13_MIN_EVAL_PER_HEAD`,
`P13_VACUOUS_AUC_FLOOR`, `P13_ECE_GREEN`, `P13_ECE_YELLOW`, `P13_ECE_NBINS`, `P13_GATE_STATISTIC`,
`P13_TOST_REQUIRED`, `P13_CONFUSION_RULE`, `P13_CONTINUITY_GATED`, `P13_TAU`,
`P13_ITERATION_ALLOWANCE`, `P13_ALPHA_GATED`, `P13_REAL_IS_GATED`, `P13_REAL_QUALITATIVE_ONLY`,
`P13_REAL_READ_ONLY`, every `P13_*_SEED`, every `P13_*_COUNTER`, `P13_ALPHA_LADDER`,
`P13_ALPHA_MASK_FRACTION_BOUNDS`, `P13_ALPHA_MAX_VALUE_BOUND`, `P13_ALPHA_INVARIANT_RTOL`,
`P13_REAL_ANCHOR_TOL`, `P13_REAL_OOD_SHIPPED_THRESHOLD`, `P13_REAL_LAMBDA_*`.

`P13_ITERATION_ALLOWANCE` stays **1 of 1, UNSPENT**, and this amendment does not spend it and
cannot: `13-SC2-AMENDMENT.md` §7 scopes the allowance to **the gate** and states that it "can
**never** be spent in response to a real-image observation." Nothing here is a re-run of the gate,
nothing retrains a net, and no weight changes.

### The gating arm is grep-provably blind to this

`13-SC2-AMENDMENT.md` §5 records exactly one gating arm: simulator ground truth, plan 13-12.
Counted occurrences of `P13_REAL_CHANNEL_PAIR` or `P13_REAL_REDUNDANCY_PAIR` [VERIFIED at HEAD]:

| File | Arm | Occurrences |
|---|---|---|
| `spike/p13/run_three_way_gate.jl` | **the gate (13-12)** | **0** |
| `spike/p13/train_three_way.jl` | training (13-11) | **0** |
| `spike/p13/run_alpha_series.jl` | simulated α-ladder (13-13) | **0** |
| `spike/p13/labels.jl` | D-05 label rule | **0** |
| `spike/p13/datagen.jl` | labelled pool | **0** |
| `spike/p13/net.jl` | architecture / persistence | **0** |
| `spike/p13/alpha_series.jl` | the shared α transform | **0** |
| `spike/p13/run_p13_realimage.jl` | **the real arm (13-16)** | **16 + 8** (on 14 + 6 lines; `:588` carries both) |
| `spike/p13/real_images.jl` | ingestion (13-15) | **11** — 5 keyword defaults (`:206, 237, 290, 371, 423`), 5 docstring echoes (`:183, 225, 273, 348, 393`), 1 header prose line (`:40`) |

**Count provenance, because these two rows were wrong in the first draft of this document.** An
independent recount on 2026-07-29 found the two CONSUMING rows undercounted — they read "10 + 4" and
"8, all as keyword defaults" — and they are corrected above. The correction makes the amendment's
blast radius **larger**, not smaller, and the `real_images.jl` row is no longer "all as keyword
defaults": 5 of its 11 are docstring echoes and 1 is header prose, which is why §6.4a now assigns
that file remediation instead of listing it. **The ZERO rows — the load-bearing ones, the ones that
carry the claim that the gating arm cannot see this change — were exact and are re-verified exact.**
Commands, so a reader can repeat them:
`grep -o 'P13_REAL_CHANNEL_PAIR' <file> | wc -l` and the same for `P13_REAL_REDUNDANCY_PAIR`.

**The phase's gate verdict is untouched by this amendment**: per-class AUC 0.990169 (coloc) /
0.988262 (exclusion) against the 0.90 floors, ECE 0.0119808 / 0.012886 against the 0.05 green band,
6/6 on one run, all on simulator ground truth. Those numbers are not re-read, not re-run, and not
re-interpreted here.

### The 13-12 / 13-13 artifacts are not touched

Their reported artifacts (`spike/p13/gate_report.jld2`, `spike/p13/alpha_report.jld2`), plans and
summaries are out of scope for plan 13-17 and are named as such in its constraints. The one way
this amendment could have reached them — the shared `h.consts_sha == p13_consts_sha()` train-time
guard in all three runners — is the open question escalated in §7, and it is escalated **precisely
because** resolving it the wrong way would widen the blast radius to the gating arm's
reproducibility.

---

## 4. Fact 4 — the honest consequence is worse-or-equal, and is stated in those terms

Unmasked, on the corrected pair, the two fixtures **barely separate**:

| Fixture | Superseded, unmasked c1/c2 | **Amended, unmasked c2/c3** |
|---|---|---|
| positive | +0.3292 | **+0.4603** |
| negative | +0.2481 | **+0.3815** |
| separation | 0.0811 | **0.0788** |

**And a third option existed that this amendment declined.** The masked read of the *same amended
pair* separates the fixtures by **0.8576** (`+0.8238 / −0.0338`), which `ghat.jl:69-71` calls the
sharper separator and which no prior source forbids. §5.3 forbids it here, on a code property; §2.1
states why that refusal — taking 0.0788 when 0.8576 was available and permitted — is the
load-bearing check that this is a correction rather than a relaxation.

The separation does not improve. Both fixtures move *further into positive correlation*, which is
if anything a stronger illustration of `13-SC2-AMENDMENT.md` §5(d)(ii): **every unmodified real
image this project has ever measured sits at positive induced μ, including the "negative"
biological control.** The corrected pair sharpens that statement rather than rescuing anything.

**The correction is therefore expected to leave the arm's conclusion where it was:** qualitative,
n = 2 specimens, no ground-truth label, OOD-bound, weak evidence. Plan 13-17 is required to say so
in words if that is what it finds (§10), and is forbidden from presenting the correction as a
rescue.

**What is expected to change is the headline OOD number**, and that is a change, not a restatement
— see §6.4.

---

## 5. The amendment, exactly

### 5.1 The constants

| Constant | Frozen (SUPERSEDED) | **Amended** | Authority |
|---|---|---|---|
| `P13_REAL_CHANNEL_PAIR` | `(1, 2)` | **`(2, 3)`** | user ruling `f411ee7`; `ghat.jl:59` `ANCHOR_CHANNEL_PAIR = (2, 3)` |
| `P13_REAL_REDUNDANCY_PAIR` | `(1, 3)` | **DROPPED — no replacement** | §5.2 |
| `P13_REAL_ANCHOR_MBAR.positive` | `0.3292` | **`0.4603`** | `ghat.jl:66` (unmasked, 64/64 patches) |
| `P13_REAL_ANCHOR_MBAR.negative` | `0.2481` | **`0.3815`** | `ghat.jl:66` (unmasked, 64/64 patches) |
| `P13_REAL_ANCHOR_TOL` | `1e-3` | **`1e-3` — unchanged** | an ingestion regression tolerance, not a bar |

`P13_REAL_ANCHOR_MBAR` keeps its meaning exactly: a **read-chain identity check** proving the
Phase-13 ingestion is the same load path the frozen simulator calibration used. It is not, and
after this amendment still is not, a Phase-13 pass/fail bar. `P13_REAL_ANCHOR_TOL` remains one of
the two documented exemptions in the not-a-gate source assertion; the second,
`P13_REAL_OOD_SHIPPED_THRESHOLD`, is also unchanged. **The exemption set stays exactly those two**,
and plan 13-17 must not introduce a third bar-shaped `P13_REAL_*` name — `test_p13_consts.jl:272-274`
enforces that by set equality and must keep passing.

### 5.2 Why the redundancy arm is dropped rather than repointed

With `c1` excluded and the fixtures carrying exactly three channels, **there is no second pair**.
The frozen `(1,3)` was DAPI-versus-red: it measured `−0.0925` masked on the *positive* fixture,
which is noise, not a corroborating read. A redundancy arm that corroborates nothing is worse than
no redundancy arm, because it is quoted as a second observation.

Dropping it **reduces** what the phase reports. It removes a table, a printed section, a persisted
artifact field and a paragraph of the report. An amendment that deletes one of its own arms is not
an amendment that is buying itself evidence.

The remaining arm is a single channel pair on two specimens. §10 requires that to be stated as the
weakening it is.

### 5.3 UNMASKED, not masked — and the masked figures are the trap

**The amended anchors are the UNMASKED values `0.4603 / 0.3815`. The masked values
`+0.8238 / −0.0338` are FORBIDDEN in this arm, however much better they separate.**

The reason is a property of the code, not a preference:

- The v2.0 amortized path applies **no Otsu mask**. `patch_summary` (`spike/contract.jl:85`, matching
  `src/amortized/summary.jl:54`) reads raw per-patch correlations. The net was **trained** on
  unmasked summaries.
- Feeding it masked summaries would present the net with inputs from a distribution it never saw.
  Given that the fixtures are **already** flagged out-of-distribution at ~2.5× on unmasked input,
  masking would push them further out and the resulting log-BFs would be measuring the mask, not
  the specimens.
- `ghat.jl:76-83` additionally records the masked estimator's own range-restriction caveat: each
  channel is thresholded independently and `_exclude_zero` then keeps only pixels bright in **both**
  — a selection on both variables — leaving only ~20-22 of 64 patches above the ≥15-survivor floor.

So the sharper-looking number is the *wrong* number twice over: wrong distribution for the net, and
a selection-biased statistic. **Choosing the worse-separating figure because it is the correct one
is the whole content of §0.** It is stated here so a reader can check that the choice went against
the arm's apparent interest.

**Be precise about whose prohibition this is.** `ghat.jl` does **not** forbid the masked read — it
records both reads on the corrected pair, calls the masked pair the sharper separator (`:69-71`) and
retracts its own earlier anti-masked framing as overstated (`:79-83`). The prohibition in this
section is therefore **this amendment's own choice, made against its own arm's interest**, not an
inherited constraint it had no say in. §2.1 states that as the load-bearing check. A reader who
wants to attack this document should attack §2.1, because if the masked figures were adopted
anywhere the whole legitimacy argument collapses — and that is greppable.

### 5.4 What follows mechanically from the pair change, and must not be silently absorbed

The α-ladder's object mask is `P13_ALPHA_MASK_RULE = :calculate_mask_ch1_once` — the Otsu mask of
**the first channel of the pair**, held fixed across the ladder. On the amended pair that is **c2
(green)**, not c1 (blue). Two consequences, both of which plan 13-17 must report rather than
absorb:

1. **The measured mask fractions change.** The frozen record (13-RESEARCH §J5.1) measured 13.49%
   (positive) and 22.92% (negative) **on c1**. The c2 mask fraction is unmeasured.
2. **`P13_ALPHA_MASK_FRACTION_BOUNDS = (0.01, 0.40)` does not move.** If the realized c2 mask
   fraction falls outside that band, `verify_alpha_invariants` returns `mask_fraction_ok = false`,
   the runner takes its `INVARIANT VIOLATION` branch, the ladder curves are **not read**, and the
   artifact is still persisted with the violation on the record. **That outcome is a reportable
   finding and is NOT a licence to widen the bounds.** Widening a pre-registered bound after a
   measurement is exactly the act this document is not.

---

## 6. Affected-site inventory (VERIFIED at `d7195cc`)

### 6.1 Sites that make reproducing the SUPERSEDED numbers a pass condition

These currently assert the wrong values and are the reason a hand-patch was started and
deliberately reverted:

| Site | Current assertion |
|---|---|
| `spike/test/test_p13_consts.jl:235` | `@test P13_REAL_CHANNEL_PAIR == (1, 2)` |
| `spike/test/test_p13_consts.jl:236` | `@test P13_REAL_REDUNDANCY_PAIR == (1, 3)` |
| `spike/test/test_p13_consts.jl:241-242` | `@test P13_REAL_ANCHOR_MBAR.positive == 0.3292` / `.negative == 0.2481` |
| `spike/test/test_p13_consts.jl:257` | `@test occursin("0.2481", P13_REAL_NAMING_CORRECTION)` |
| `spike/test/test_p13_real.jl:144-146` | `real_mbar("positive"/"negative") ≈ P13_REAL_ANCHOR_MBAR` at the DEFAULT pair |
| `spike/test/test_p13_real.jl:314` | `@test occursin("+0.2481", P13_REAL_NAMING_CORRECTION)` |
| `spike/p13/consts.jl:708-709` | `@assert occursin("biological", …) && occursin("0.2481", P13_REAL_NAMING_CORRECTION)` |
| `13-01-PLAN.md:485` | acceptance criterion pinning the bar-name exemption set (must keep passing unchanged) |
| `13-01-PLAN.md:558-561, 568` | verify criteria asserting `(1,2)`, `(1,3)`, `0.3292`, `0.2481` |
| `13-15-PLAN.md:17` | a `must_have` asserting the SUPERSEDED `ghat` anchors |
| `13-15-PLAN.md:354, 515` | acceptance criteria asserting `0.3292` / `0.2481` |

### 6.2 Constants and prose in `spike/p13/consts.jl`

`:482` (α-substrate comment quoting `+0.3292 / +0.2481`), `:552` pair, `:553` redundancy pair,
`:576` `P13_REAL_ANCHOR_MBAR`, `:630` (`P13_REAL_NAMING_CORRECTION` body quoting `+0.2481`),
`:709` (the assertion requiring the string to contain `"0.2481"`).

Note `:546-548` — the section comment claiming the primary pair "is the EXACT pair the frozen
`ghat` calibration anchor was measured on." After 2026-07-25 that sentence was already false; it is
the drift of §2 written down in one line.

### 6.3 Code that breaks when the redundancy pair is dropped

`spike/p13/run_p13_realimage.jl`:

- `:50` header prose freezing both pairs
- `:588` prints `P13_REAL_REDUNDANCY_PAIR` in the locked-value banner
- `:724-727` the `sweep_rows("primary", P13_REAL_CHANNEL_PAIR)` call and its section label
- `:734-735` the shipped-OOD comparison read
- `:797`, `:830-836`, `:852` the α-ladder reads and the persisted `channels` field
- `:905` prints `P13_REAL_ANCHOR_MBAR` in the never-averaged note
- `:926-940` **SECTION 5 in full** — the `[5/6]` progress line, the `sweep_rows("redundancy", …)`
  call, `print_sweep`, and the paragraph naming the counterstain limit
- `:1005` the persisted `redundancy_channels` artifact field
- `:1069` the returned `redundancy_sweep` field
- `:1007-1008` `real_provenance(…; channels = P13_REAL_CHANNEL_PAIR)`

`sweep_rows(arm, channels)` at `:665` is already parameterized on the pair, so dropping the arm is a
deletion of the second call and its section, not a refactor.

### 6.4 Prose quoting the superseded anchors

`spike/p13/real_images.jl:35, 44, 68-69, 277-283`; `spike/p13/run_p13_realimage.jl:45`;
`spike/p13/run_alpha_series.jl:97, 399-400, 674`; `spike/test/test_p13_real.jl:35, 137`;
`13-REPORT.md:30, 568, 595, 684, 857, 992`; `13-16-SUMMARY.md:68, 132, 215, 267, 276, 406`;
`13-15-SUMMARY.md:69-72, 150, 222-224, 230, 236, 242`.

#### 6.4a REMEDIATED — `spike/p13/real_images.jl` is assigned to a task, not merely listed

**This file was in no task's file list until 2026-07-29, and that was a defect in this document.**
Listing a stale site in an inventory without assigning it remediation is the same drift class this
amendment exists to repair, one level down: §2's whole argument is that a corrected upstream was
recorded and then **not carried into the file that reads it**.

Three blocks here are load-bearing prose and each goes stale, differently under each §7 mechanism:

| Site | Block | Under **M1** | Under **M2** |
|---|---|---|---|
| `:39-49` | (b2) the channel-pair caveat; `:40` opens `P13_REAL_CHANNEL_PAIR IS FROZEN AT c1/c2.` | **literally FALSE in shipped source** | literally true, but **misleading about what is operative** — the reader learns what is frozen, not what is read |
| `:66-71` | (e) THE LINEAGE CLAIM — records the c1/c2 anchors `0.3292 / 0.2481` as *the* claim, unqualified | anchors no longer bound to those values | anchors still bound, but no longer the colocalization reference |
| `:277-283` | the `real_mbar` docstring, restating the frozen lineage expectation | same as (e) | same as (e) |

**Assigned to `13-17-PLAN.md` Task 3, with explicit M1/M2-conditional wording, as a COMMENT-ONLY
edit.** No executable line changes: the five `channels = P13_REAL_CHANNEL_PAIR` keyword defaults are
**not** re-pointed — under M2 they cannot be, and under M1 they follow automatically — and the
runner instead passes `channels` explicitly at every call site, which is what makes the operative
pair visible at the point of use. The lineage claim in block (e) is **kept, not withdrawn**: what it
asserts is the identity of the read chain, and that identity holds on whichever pair it is exercised
on. Only its status as *the colocalization reference* is superseded.

#### 6.4b LEFT DELIBERATELY — the frozen sibling documents

**`13-SC2-AMENDMENT.md:343, 381, 413, 423, 478` also quote `+0.3292 / +0.2481`. Those are NOT
edited.** That document is FROZEN and its numbers were correct as prior measurements at the time it
was written; this document supersedes them by reference, which is how the two amendments stack.
The same applies to `13-RESEARCH.md` and `13-VALIDATION.md`: prior-measurement records, cited, not
rewritten.

### 6.5 The headline OOD number is pair-dependent and WILL move

The reported binding limit of the whole arm is the OOD flag. Measured on the superseded pairs
[`13-16-SUMMARY.md`, VERIFIED]:

| Arm / pair | Phase-13 density | ID threshold | ratio | shipped density | shipped threshold | ratio |
|---|---|---|---|---|---|---|
| primary `(1,2)` | **417.2974** | 167.5446 | **2.491×** | 433.6884 | 179.1368 | 2.421× |
| redundancy `(1,3)` | **607.5728** | 167.5446 | **3.63×** | — | — | — |

Per-acquisition, primary: positive 417.2974 / negative 321.1505. Redundancy: positive 410.24 /
negative 607.57. **A different channel pair already moved this number by 46%.** The amended `(2,3)`
value is unmeasured and must be **reported as a change**, never restated from the old run.

Two knock-on notes for the re-run:

- `P13_REAL_OOD_SHIPPED_DENSITY = 433.69` is a **measurement on the superseded pair**. It stays in
  `consts.jl` as a frozen historical reference, but the corrected run's shipped-detector reading on
  `(2,3)` will differ from it and **that difference is a result, not a discrepancy to reconcile.**
  `P13_REAL_OOD_SHIPPED_THRESHOLD = 179.14` is a property of the *bundle*, not of the pair, and is
  unaffected — it remains the runner's bundle selector (`run_p13_realimage.jl:352`), so the
  pre-registration is still what picks the comparison reference.
- `run_p13_realimage.jl:741` falls back to `P13_REAL_OOD_SHIPPED_DENSITY` when no shipped bundle is
  readable. After this amendment that fallback would quote a superseded-pair number; if it fires it
  must be labelled as such.

### 6.6 `.planning/ROADMAP.md` — the phase's most-read summary surface, and it was missing here

**`.planning/ROADMAP.md` quotes the superseded-pair result in four places and was absent from this
inventory.** It is the surface a reader hits first — before the report, before any summary — and §9's
rule ("any report that quotes an amended real-image number must quote the superseded one beside it")
applies to it with more force than to anything else listed above, not less.

| Site | What it quotes | Status |
|---|---|---|
| `:176` | Phase-13 table row: *"13-16: real-image arm REPORTED (RANDOM on both pairs, both OOD-flagged 2.49x/2.42x, agreement)"* | measured on the **SUPERSEDED** pairs; "both pairs" refers to the dropped redundancy arm |
| `:275` | `**Plans**: 16 plans (9 waves …)` | stale — becomes **17 plans (10 waves)** |
| `:291` | the `13-16` plan bullet: `417.30 vs 167.54 (2.49×)`, `433.69 vs 179.14 (2.42×)`, `alpha_star_real = 0.875 / 0.625`, *"RANDOM on both unmodified pairs"*, contrast ±0.1156 | every figure measured on the **SUPERSEDED** pair |
| `:293` | the `**Verdict:**` block, which correctly reports the **gate** but folds the real arm into the same sentence | the gate half is untouched by this amendment and must stay so; the real-arm half is superseded |

**Split across two owners, deliberately:**

1. **Now, by the planner** — the part that is knowable without the re-run: the plan count and wave
   count, a `13-17-PLAN.md` bullet, and a SUPERSEDED marker on every quoted real-image figure. Done
   as the smallest additive edit, pathspec-scoped, with the file re-read immediately before writing
   **because a Phase-12 executor is live on this same file**.
2. **After the re-run, by `13-17-PLAN.md` Task 5** — the amended figures, each printed beside its
   superseded counterpart with the pair named, per §9. `.planning/ROADMAP.md` is in that task's file
   list, and the same re-read-immediately-before-writing discipline applies.

**The gate row and the `**Verdict:**` block's gate sentences are NOT touched by either owner.** The
gate's AUC/ECE figures are simulator-ground-truth results that this amendment cannot reach (§3).

---

## 7. THE ONE OPEN QUESTION — escalated, not resolved here

**Where the amended constants live is NOT settled by the user's ruling, and this planner does not
settle it.** It is escalated as a blocking decision at Task 1 of `13-17-PLAN.md`.

### The finding that forces the question

`spike/p13/run_p13_realimage.jl:634` asserts:

```julia
@assert h.consts_sha == p13_consts_sha() "the frozen pre-registration changed since the net was trained: consts.jl sha256 does not match the artifact's"
```

`h` is the **trained net artifact from 13-11**, whose `consts_sha` was recorded at train time
(`spike/p13/net.jl:472, 524`). The identical assertion exists at `run_three_way_gate.jl:353` (the
**gate** runner) and `run_alpha_series.jl:451`.

**Therefore: editing `spike/p13/consts.jl` in place makes the 13-16 re-run abort at load**, and
simultaneously makes the gate and α-ladder runners un-reproducible against their own net. Retraining
is not available (no further training is authorised; retraining would invalidate 13-12's gate).

`spike/test/test_p13_net.jl:380` is **not** affected — it round-trips a net saved inside the test,
so its sha is recomputed at save time. Verified.

### The two options, both defensible

**M1 — amend in place.** Edit `consts.jl` section I2; widen the three train-time guards to accept
the recorded pre-amendment sha alongside the current one, justified by the fact that **training
never read any `P13_REAL_*` constant** (grep-verified, §3), so the amended values are ones the net
did not consume.
*Cost:* touches an integrity guard in order to produce a result — the shape this project is most
careful about — and retires the standing claim that `consts.jl` has "two commits in its whole
history, neither postdating a result" (`13-REPORT.md` §13).

**M2 — amend beside.** Leave `consts.jl` **byte-unchanged**; put the amended values in a new,
separately named, separately committed `spike/p13/consts_d15_amendment.jl` with its own provenance
header and its own not-a-gate assertion. The runner reads the amended names and prints the
superseded frozen values beside them.
*Cost:* two files must be read to know the operative pair, and the module must carry an explicit
"this supersedes `consts.jl:552-576`" statement or it is a trap of its own.
*Benefit:* every sha guard keeps holding untouched; the gating arm's reproducibility is untouched;
the superseded pre-registration stays byte-present and citable, which is exactly the discipline
`13-SC2-AMENDMENT.md` §6 requires ("The original is not erased"); and
`test_p13_consts.jl:235-242` keeps asserting the frozen values **truthfully**, with the new values
asserted in a new testset.

**This planner's reading is that M2 dominates on every stated constraint, and this planner does not
adopt it.** The choice is a provenance-narrative decision about a frozen pre-registration, and this
project's credibility is the thing being spent. It is the user's.

---

## 8. What this amendment does NOT do

- It does **not** reopen the Phase-7 ship gate, and does not touch `amended_v2` or `grid_8`.
- It does **not** open the sealed holdout. `open_sealed_holdout` is not called; no path under the
  sealed provenance tree is constructed; the Phase-16 seal stays SHUT. The one physically segregated
  real anchor remains deferred to Phase 16.
- It does **not** re-run, re-read or re-interpret 13-12 (the gate) or 13-13 (the simulated α-ladder).
- It does **not** spend `P13_ITERATION_ALLOWANCE`, retrain any net, or change any weight.
- It does **not** touch `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `test/`, `docs/` or
  `corpus/`.
- It does **not** edit `13-SC2-AMENDMENT.md`, `13-RESEARCH.md` or `13-VALIDATION.md`. Those record
  prior measurements that were correct as prior measurements; this document supersedes them by
  reference.
- It authorises **no further amendment**. A third amendment is not authorised by this document and
  cannot be authorised by amending this document — the same clause `13-SC2-AMENDMENT.md` §7 carries,
  for the same reason.

---

## 9. The original is not erased

`spike/p13/consts.jl`'s frozen values remain citable in git at `c42cc8e` and `bfba6ac` whichever
mechanism §7 selects, and under M2 they remain byte-present in the working tree. `13-SC2-AMENDMENT.md`
keeps its numbers. `13-16-SUMMARY.md` keeps its measured tables in full.

**Any report that quotes an amended real-image number must quote the superseded one beside it**, and
must say which pair each was measured on. Plan 13-17 makes that a checked property of `13-REPORT.md`,
not an intention.

**"Any report" includes `.planning/ROADMAP.md`, and that surface matters most, not least.** It is
the first thing a reader reaches — before the report, before any summary — and it was missing from
this inventory until 2026-07-29. §6.6 lists its four sites and splits the work: the planner marks
them SUPERSEDED and corrects the plan/wave count now; plan 13-17 Task 5 prints the amended figures
beside them after the re-run, on its own pathspec, with the file re-read immediately before writing
because a Phase-12 executor is live on it.

---

## 10. What the re-run is required to answer, in words

The prior reading of the real arm is: **qualitative, n = 2 specimens, no ground-truth label,
OOD-bound, weak evidence, coherence and never correctness.**

Plan 13-17 must require `13-REPORT.md` to answer, in plain prose and near the top of the real-image
section: **does the conclusion change?**

- **If the corrected pair lands in the same place, the report must SAY SO PLAINLY** — *"we corrected
  the channel pair and the conclusion is unchanged"*. That is a **good** result. It tells a reader
  the limitation is **structural** — a property of n = 2 unlabelled specimens outside the training
  distribution — rather than an artefact of having read the wrong channels. A limitation that
  survives its own correction is a stronger statement than one that was never tested.
- **If it changes, the change must be stated as a change**, with both numbers and both pairs, and
  the reader told which is superseded.
- **Overselling the correction as a rescue is forbidden.** The separation got marginally *worse*
  (§4). No sentence may imply the corrected arm is strong evidence, and the report may not drop or
  soften any of: n = 2, no label, both fixtures OOD-flagged, the α-ladder's construction-parameter
  caveat, or the Phase-16 deferral.

---

## Appendix — provenance of every number quoted above

| Number | Source | Kind |
|---|---|---|
| `["blue", "green", "red"]` | `test/runtests.jl:105` | property of the committed fixtures |
| `(1,2)`, `(1,3)`, `0.3292`, `0.2481` | `spike/p13/consts.jl:552, 553, 576` | the SUPERSEDED pre-registration |
| `(2,3)`, unmasked `0.4603 / 0.3815` | `spike/simulator/ghat.jl:59, 66` (frozen, `78dc37f`) | prior measurement, frozen simulator calibration |
| masked `+0.8238 / −0.0338`, `−0.0925` | `spike/simulator/ghat.jl:68`; STATE.md quick task `260725-wb7` | prior measurement — **FORBIDDEN in this arm**, §5.3 |
| 417.2974 / 167.5446 / 2.491×; 433.6884 / 179.1368 / 2.421×; 321.1505 | `13-16-SUMMARY.md` §1, `realimage_report.jld2` | Phase-13 result on the SUPERSEDED pair |
| 607.5728 / 3.63×; 410.24 / 607.57 | `13-16-SUMMARY.md:139-144` | Phase-13 result on the SUPERSEDED redundancy pair |
| AUC 0.990169 / 0.988262; ECE 0.0119808 / 0.012886 | `13-12-SUMMARY.md`, `gate_report.jld2` | Phase-13 GATE result — **not re-read, not re-run** |
| mask fractions 13.49% / 22.92% | 13-RESEARCH §J5.1 | prior measurement **on c1** — superseded as a c2 expectation, §5.4 |
| every commit sha and timestamp in §0.1 | `git log`, copied from output | repository history |
| `d7195cc4c5ba5a4aaceef8f2d59ddf9bd4422be1` | `git rev-parse HEAD` at the time of writing | the sha this inventory was verified at |

| separations `0.0811`, `0.0788`, `0.8576` | subtraction of the paired figures in the rows above | **arithmetic on prior measurements**, not a measurement |
| `39 h 15 min 06 s` (`78dc37f` → `c42cc8e`), `38 h 46 min 18 s` (`ea4a7d3` → `c42cc8e`), `40 min 49 s` (`2112bed` → `f411ee7`) | differences of the `git log` timestamps in §0.1 | arithmetic on repository history |
| occurrence counts in §3 | `grep -o '<name>' <file> \| wc -l`, re-verified 2026-07-29 | property of the source tree |

**No number in this document was measured by this document.** Every Phase-13 figure quoted is a
prior result, labelled with the pair it was measured on, and the amended anchors are prior
measurements of the frozen simulator calibration — not new readings taken to justify the change.
