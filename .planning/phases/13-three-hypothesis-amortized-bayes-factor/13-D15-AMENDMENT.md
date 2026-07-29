# 13 — PRE-REGISTRATION AMENDMENT: the real-image channel pair (CHANGE A) and the include-guard sentinel (CHANGE B)

**Date:** 2026-07-29
**Status:** **FROZEN.** §7 is RESOLVED — the user ruled **M1, amend in place**, on 2026-07-29.
Every value in §5 was already settled by the user's ruling and is not reopened.
**Implements:** the user ruling of 2026-07-29 (`f411ee7`) on the channel pair, restating and
completing the factual correction first recorded 2026-07-25 (`ea4a7d3`); and the user ruling of
2026-07-29 on the include-guard mechanism (M1, §7.1)
**Amends `spike/p13/consts.jl` in ONE edit, with TWO separately justified changes (§0.4):**
**CHANGE A** — section I2: `P13_REAL_CHANNEL_PAIR`, `P13_REAL_REDUNDANCY_PAIR`,
`P13_REAL_ANCHOR_MBAR` and the prose constants that quote them.
**CHANGE B** — the body wrapper at `:100`: the include-guard sentinel, together with
`spike/test/test_p13_consts.jl:42` and six sibling callers.
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

### 0.4 TWO CHANGES, TWO JUSTIFICATIONS, DISCLOSED SEPARATELY — read this before the diff

**This amendment opens `spike/p13/consts.jl` once and makes TWO changes to it. They are different
acts with different reasons, and neither is a threshold moved after proving inconvenient.** They
are separated here, on the page, so that a reader who was not present does not have to infer the
second one from a diff.

| | **CHANGE A — the channel pair** (§5) | **CHANGE B — the include-guard sentinel** (§5B) |
|---|---|---|
| What moves | `P13_REAL_CHANNEL_PAIR` `(1,2)` → **`(2,3)`**; `P13_REAL_REDUNDANCY_PAIR` **DROPPED**; `P13_REAL_ANCHOR_MBAR` `0.3292 / 0.2481` → **`0.4603 / 0.3815`** (UNMASKED) | the body wrapper at `:100` guards on **`:P13_DECLARED_DEVIATIONS`** instead of `:P13_DEV_SEED`, and `spike/test/test_p13_consts.jl:42` plus six sibling callers move with it |
| What is WRONG with the text today | it **names the WRONG PHYSICAL OBJECT.** `test/runtests.jl:105` records the fixture channels as `["blue", "green", "red"]`, so `c1` is the DAPI/Hoechst nuclear counterstain — **not a target protein at all** — and the arm as executed measured counterstain-versus-protein overlap, which is not colocalization | it guards on a name **another phase must LEGITIMATELY MIRROR.** `spike/validation/p12_consts.jl:109` declares `P13_DEV_SEED` in order to assert seed disjointness. **That mirror is CORRECT and `p12_consts.jl` must NOT be touched** (append-only Tier 1). The defect is on the guard side, every time |
| Kind of error | a factual mis-designation of the object being measured | a latent loading defect: in a full-suite run the wrapper sees the mirror, skips its own body, and ~114 Tier-1 constants are never defined (`runtests.jl:220`, 22 pass / 1 fail / **114 `UndefVarError`**) |
| Does any bar move? | **NO** — §3 enumerates every untouched threshold by name | **NO** — CHANGE B changes no value at all. `P13_DEV_SEED = 0x0000_0000_0B13_DE71` stays byte-present and byte-unchanged at `:188`; only *which name the `if` tests* moves |
| Would a blind reviewer have made it? | **Yes, and one did, twice, before Phase 13 existed** (§2, §0.3) | **Yes** — it is the `fb76b84` precedent (*"guard on the name it actually owns"*) applied to the one case where the owner file is the file that needs editing |
| Why it rides in THIS amendment | — | **ZERO additional pre-registration cost.** `consts.jl`'s sha256 changes for CHANGE A regardless, and `h.consts_sha == p13_consts_sha()` has to be re-derived in the gate, alpha and real-image runners regardless (§7.3). CHANGE B adds no new consequence to that. Holding it back would mean opening the byte-locked file a SECOND time later, under a SECOND amendment — strictly worse |

**Neither change is a threshold moved after proving inconvenient, and both are written so they read
that way to a hostile reviewer.** No seed, bar, floor, band, ladder, tolerance or allowance changes
in either. §3 is the greppable form of that claim and it covers both.

**The coupling is the reason §7 resolves to M1 and not M2.** A "leave the frozen file alone"
mechanism cannot fix CHANGE B at all: the broken end **is** `consts.jl:100`, and a caller-only fix
is provably useless — the caller would correctly call `include`, and the wrapper would still see
`P13_DEV_SEED` defined by Phase 12 and skip the entire body.

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

## 5. CHANGE A, exactly - the real-image channel pair

*(CHANGE B - the include-guard sentinel - is §5B. The two are disclosed separately on purpose;
see §0.4.)*

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

## 5B. CHANGE B, exactly — the include-guard sentinel

### 5B.1 The defect

`spike/p13/consts.jl:100` wraps the file's **entire Tier-1 body** (lines 100-731, ~104 `const`s) in

```julia
if !isdefined(@__MODULE__, :P13_DEV_SEED)
```

and `:188` is where `P13_DEV_SEED` is itself declared. The sentinel is the file's own seed — and
`spike/validation/p12_consts.jl:109` **also** declares it:

```julia
const P13_DEV_SEED        = 0x0000_0000_0B13_DE71  # spike/p13/consts.jl:188
```

**That mirror is CORRECT and `spike/validation/p12_consts.jl` MUST NOT be touched.** It is
append-only Tier 1, and asserting stream disjointness *requires* naming the streams:
`spike/test/test_p12_consts.jl:80` asserts `P12_DEV_SEED != P13_DEV_SEED`, which is only possible
if Phase 12 names Phase 13's seed. **The defect is on the guard side, every time.**

The consequence, traced end to end and reproduced empirically (`.planning/STATE.md`, 2026-07-29):
`runtests.jl:174` includes `test_p12_suite.jl` FIRST — deliberately, so Phase 12 cannot be masked
by anyone else's abort — and that chain defines `P13_DEV_SEED` in `Main`. At `runtests.jl:220`,
`test_p13_consts.jl:42`'s guard sees the mirrored symbol, **skips its `include` entirely**, the
wrapper at `:100` never runs, and Phase 13 reports **22 pass / 1 fail / 114 error**, every error an
`UndefVarError` on a Tier-1 constant. The thrown testset then aborts the **eight** remaining
Phase-13 includes. **Per-file runs pass**, which is exactly why this stayed invisible: in a
single-file process nothing else has defined the sentinel yet.

### 5B.2 Why this is not a threshold move — and the check that separates them

The check §1 uses applies unchanged: *a relaxation changes how well a measured quantity must score;
this changes nothing that is measured at all.* CHANGE B alters **which symbol an `if` tests**. It
changes no value, reads no result, and cannot be reached by any measurement.

Mechanically:

- `P13_DEV_SEED = 0x0000_0000_0B13_DE71` stays **byte-present and byte-unchanged** at `:188`;
- no constant is added, removed or re-valued by CHANGE B;
- §3's enumerated untouched-threshold list is unaffected;
- the change is **outcome-independent in the strongest sense available**: it was diagnosed from a
  suite abort rather than from any Phase-13 number, and it makes the suite report **more**
  failures, not fewer — it un-masks eight Phase-13 test files that are currently skipped.

**A change that increases the number of live assertions is not a change that is buying itself a
result.** After it lands, `test_p13_correction.jl`'s two pre-registered MEASURED MISSES still fail
and the full suite still exits 1. That is the intended state and it is not repaired here.

### 5B.3 The sentinel, chosen and ARGUED

`spike/p13/consts.jl` exclusively declares **100** constants (each verified by a fresh repo-wide
`grep -rn "const NAME\b" --include=*.jl .` returning exactly one declaring file). The wrapper is
re-pointed to:

> **`P13_DECLARED_DEVIATIONS`** — `spike/p13/consts.jl:662`

against the four criteria in `.planning/CONVENTIONS.md` C-01:

1. **Exclusively owned.** Fresh repo-wide grep: `const P13_DECLARED_DEVIATIONS` appears at
   `spike/p13/consts.jl:662` and **nowhere else**; the only other file that so much as mentions the
   name is `spike/test/test_p13_consts.jl`. It does **not** appear in
   `spike/validation/p12_consts.jl`.
2. **Tier-1 and declared UNCONDITIONALLY inside the guarded body.** `:662` sits inside `100-731` at
   guard-block indent, with no nested `if` / `begin` between `:100` and it.
3. **NOT a seed and NOT a prior bound.** That is the whole lesson of the sweep: **23 of the 29
   poisoned guards guard on a seed or a prior bound**, because those are exactly the names other
   phases are *required* to mirror in order to prove disjointness.
4. **About the FILE'S IDENTITY, not a value another phase might quote.** It is the
   pre-registration's own register of its four declared model deviations — a meta-property of the
   document. A Phase-14 deviation register would be `P14_DECLARED_DEVIATIONS`; there is no
   requirement anywhere, and no motive, for another phase to declare Phase 13's.

It is also the name `.planning/STATE.md` had already named as option one when it recorded the
blocker, so this is not an invention produced at the moment of fixing.

**The candidates weighed and REJECTED, on the page:**

| Candidate | Exclusively owned? | Rejected because |
|---|---|---|
| `P13_TAU` (`:794`) | yes | **Tier-2, and legitimately SOMETIMES ABSENT.** `consts.jl:745` explicitly tolerates its absence (`!isdefined(:P13_TAU) \|\| …`) and it lives in a *separately guarded* block outside the Tier-1 body. **A sometimes-absent sentinel is worse than the bug** — it would let the Tier-1 body run twice |
| `P13_GATE_STATISTIC` (`:442`) | yes | It is a **gate knob** (`:ece`). Making the loading of the whole pre-registration depend on a gate constant is exactly the coupling this project should not create, and a later phase asserting *"we do not use Phase 13's statistic"* is a plausible mirror |
| `P13_CUT_VARIANT` (`:255`) | yes | **Persisted into the trained-net artifact** (`net.jl`, `h.cut_variant`; asserted at `test_p13_net.jl:373`). Any later phase reading a Phase-13 net has an obvious motive to declare its own expected `P13_CUT_VARIANT` to assert artifact compatibility — the mirror hazard again |
| `P13_LAMBDA_PLACEMENT` (`:398`) | yes | Referenced from `preconditions.jl`, which is the **cross-phase compatibility surface**. A constant whose job is to be compared against another phase's is the wrong shape for a private sentinel |
| `P13_INPUT_WIDTH_RULE` (`:402`) | yes | Same class: an architecture rule a Phase-14 net would have to match, so a plausible thing for Phase 14 to re-declare |
| `P13_STRATIFICATION` (`:334`) | yes | A **design knob** (`:class_frequency`) read by five files; semantically about the data-generation design rather than the file's identity. Not unsafe today, but strictly weaker than a deviations register on criterion (4) |

None of the rejected candidates is *unsafe today*. They are rejected because criterion (4) is what
makes a sentinel durable, and only `P13_DECLARED_DEVIATIONS` satisfies it outright.

### 5B.4 BOTH ENDS, IN ONE EDIT — and the six sibling callers

**A caller-only fix is provably useless here**, and saying why is the point: the caller would
correctly call `include`, and the wrapper at `:100` would still see `P13_DEV_SEED` defined by
Phase 12 and skip the entire body. `fb76b84` (`:SBC_M`) only had to fix callers because
`validation/consts.jl` already guarded its own body on the name it owns. Here **both ends use the
poisoned name**, and this project's signature failure is a fix applied at one end of a pair. So
both ends move together, in one edit:

| # | Site | Today | After |
|---|---|---|---|
| 1 | `spike/p13/consts.jl:100` — **the owner-side wrapper** | `if !isdefined(@__MODULE__, :P13_DEV_SEED)` | `if !isdefined(@__MODULE__, :P13_DECLARED_DEVIATIONS)` |
| 2 | `spike/test/test_p13_consts.jl:42` — **the aborting caller** | `isdefined(@__MODULE__, :P13_DEV_SEED) \|\| include(…)` | `isdefined(@__MODULE__, :P13_DECLARED_DEVIATIONS) \|\| include(…)` |

**And the six sibling callers, IDENTIFIED rather than left as a later surprise.** Each guards its
own `include` of `consts.jl` on the same poisoned name, and each is a Phase-13 file carrying no sha
guard of its own:

`spike/p13/labels.jl:73`, `spike/p13/net.jl:89`, `spike/p13/preconditions.jl:191`,
`spike/p13/result.jl:79`, `spike/p13/tau_probe.jl:92`, `spike/p13/toy_gaussian.jl:79`.

They are **harmless today given that `test_p13_consts.jl` loads first and, after this amendment,
succeeds** — a skipped include is correct once the constants are already present. They are fixed
anyway, because leaving them makes correctness depend on load ORDER, and order-independence is the
entire property the guarded-include idiom exists to buy. `spike/p13/consts.jl:93` is a **commented
illustration** of the idiom, not a live guard; its text is updated to match so the file does not go
on teaching the bug it just fixed.

**NOT touched:** `spike/validation/p12_consts.jl` (§5B.1), and the 28 other poisoned guards the
same sweep found (§11).

### 5B.5 The knock-on CHANGE B makes LIVE, and that must be VERIFIED rather than assumed

While the guard was skipping the body, every `const` in that body was never executed in a
full-suite process. Repairing the guard makes ~104 of them execute into `Main` **for the first
time**, alongside declarations of the same names made earlier by other files. Ten names collide.
Seven agree in value **and** type. Three disagree in **type**, because in Julia a hex literal's
width is its type (`0xC0FFEE` is `UInt32`; `0x0000_0000_00C0_FFEE` is `UInt64`):

| Name | narrow (`UInt32`) | wide (`UInt64`) |
|---|---|---|
| `NPE_MASTER_SEED` | `npe/train_npe.jl:65`, `test/test_npe.jl:81`, `baseline/run_advi.jl:74` | `p13/consts.jl:108` and the three pre-registrations |
| `VAL_MASTER_SEED` | `validation/consts.jl:76` | `p13/consts.jl:109` and the three pre-registrations |
| `VAL_FIX_SEED` | `validation/consts.jl:79` | `p13/consts.jl:110` and the three pre-registrations |

The **values** agree; only the widths differ, and all three are FORBIDDEN foreign seeds that
Phase 13 *reserves* rather than uses. Julia 1.12.6 — this project's pinned version — permits
constant redefinition, so the expectation is that this is a non-event. **But it is an expectation,
and it has never been executed.** `13-17-PLAN.md` therefore makes it a verified step with a stop
rule rather than a hope.

**No seed literal may be edited to resolve it.** If a redefinition does fail, the remedy is on the
loading side (isolated-module read, or load order) and the executor **STOPS and escalates** rather
than touching a seed. Recorded as a standing rule in `.planning/CONVENTIONS.md` C-02.

---

## 6. Affected-site inventory (VERIFIED at `d7195cc`)

*(§6.1-§6.6 are CHANGE A's sites. CHANGE B's sites are §5B.4 and §6.7.)*

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

### 6.7 CHANGE B's sites, for inventory completeness

The nine guard sites are enumerated in §5B.4 (one owner-side wrapper, one aborting caller, six
sibling callers, one commented illustration at `:93`). The three `consts_sha` assertion sites that
CHANGE A and CHANGE B jointly trip are enumerated in §7.3. Nothing else in the repository reads the
sentinel as a sentinel — verified by
`grep -rn "isdefined(.*:P13_DEV_SEED" --include=*.jl .`, which returns exactly those nine lines.

---

## 7. RESOLVED — M1 (amend in place). The ruling, its reasoning, and the sha re-derivation

**Status: SETTLED by the user on 2026-07-29. The mechanism is M1 — edit `spike/p13/consts.jl` in
place — and BOTH changes of §0.4 ride in that ONE edit.** This question is not reopened. The text
below replaces the escalation that previously stood here; the two options as they were put to the
user are preserved verbatim in git and in `13-17-PLAN.md`'s history.

### 7.1 The ruling, and the reasoning that produced it

> **M1 — fix at the source, inside the same authorised amendment. M2 (the isolated-module read) is
> the tool for files that CANNOT be edited because no amendment authorises opening them.**

The reasoning is recorded because a mechanism chosen for a reason is auditable and a mechanism
chosen by default is not:

**The byte-lock breaks either way.** CHANGE A on its own changes `consts.jl`'s sha256, and
`h.consts_sha == p13_consts_sha()` has to be re-derived in the gate, alpha and real-image runners
regardless of how CHANGE B is handled. **The include-guard fix therefore rides along at ZERO
additional pre-registration cost.** That coupling is precisely why the answer is M1 and not M2:

- Under M2, `consts.jl` would stay byte-unchanged — but **CHANGE B could not be fixed at all**,
  because the broken end *is* `consts.jl:100`. A caller-only fix is provably useless (§5B.4). The
  defect would instead have to be carried by an isolated-module read at every caller, and
  `spike/test/test_p13_consts.jl` asserts ~114 constants UNQUALIFIED, so that read would need all
  of them re-bound into the test module — mechanical, but not free, and strictly more machinery
  than the one-line fix it replaces.
- Under M2 the phase would also end holding a *second* source of truth for the operative channel
  pair while still carrying the guard defect: two files to read, one bug outstanding.
- Under M1, one file is the pre-registration, a reader finds the operative pair where they expect
  it, and the guard defect closes inside the same authorised, dated, argued edit.

**M2 remains the right tool where it applies** — a frozen file that no amendment opens — and §11
assigns exactly that remedy to those of the other 28 hits whose owner file cannot be edited. The
ruling is not "M1 is better than M2"; it is "M1 here, M2 there, and the distinction is whether an
amendment authorises opening the file."

### 7.2 What M1 costs, stated rather than minimised

M1's price was recorded before the ruling and is not softened now.

1. **It retires a standing claim in `13-REPORT.md` §13.** That section states that `consts.jl` has
   two commits in its whole history and that *"No commit that produced a Phase-13 result touched
   it"*, and quotes `git hash-object` → `5a4ea222…` and `sha256` → `100e97a3…` as *"matches
   `consts_sha` and `net_consts_sha` in every artifact"*. **All three become FALSE the moment the
   file is edited.** `13-17-PLAN.md` therefore rewrites §13 with pre- and post-amendment values
   side by side. **The section whose entire function is proving nothing was tampered with must
   itself be corrected — and that is said plainly here rather than done quietly.**
2. **It touches an integrity guard in order to produce a result** — the shape this project is most
   careful about. §7.3 specifies the mechanism so the guard means MORE afterwards, not less.
3. **It makes `spike/test/test_p13_consts.jl:235-242` assert values that no longer exist**, so
   those assertions are rewritten to the amended values rather than relabelled. The superseded
   values stay citable in git at `c42cc8e` / `bfba6ac` (§9).
4. **It makes a lineage assertion silently change meaning unless it is pinned.**
   `spike/test/test_p13_real.jl:144-146` calls `real_mbar("positive")` with **no `channels`
   argument**, riding the default, and compares against `P13_REAL_ANCHOR_MBAR`. Under M1 **both
   sides move together**: the call becomes a c2/c3 read compared against the c2/c3 anchors, the
   assertion **still passes**, and the c1/c2 read-chain identity proof is **silently destroyed**
   while its own comment goes on claiming to test it. `13-17-PLAN.md` pins it to
   `channels = (1, 2)` against the literals `0.3292` / `0.2481`, so it keeps testing what it says.

### 7.3 The `consts_sha` re-derivation — first-class work, not a side effect

`spike/p13/net.jl:524` defines
`p13_consts_sha() = bytes2hex(SHA.sha256(read(joinpath(@__DIR__, "consts.jl"))))`, and three
runners assert it against the **13-11 trained net's** recorded value:

| Site | Runner | Arm |
|---|---|---|
| `spike/p13/run_three_way_gate.jl:353` | the gate | **the phase's ONLY gating arm** |
| `spike/p13/run_alpha_series.jl:451` | the simulated α-ladder | reported, not gated |
| `spike/p13/run_p13_realimage.jl:634` | the real arm | reported, not gated |

`spike/test/test_p13_net.jl:380` is **not** affected — it round-trips a net saved inside the test,
so its sha is recomputed at save time. Verified. `spike/p13/datagen.jl:870` and
`spike/p13/run_tau_probe.jl:96` *record* a sha but assert nothing against it. **Retraining is NOT
authorised** and would invalidate 13-12's gate result.

**The mechanism chosen, and why the guard still means something afterwards.**

A widened `h.consts_sha == p13_consts_sha() || h.consts_sha == PRE_AMENDMENT_SHA` is **rejected**.
Once that disjunction exists it accepts **any future drift silently**: a second, third or
undisclosed edit to `consts.jl` also passes, because the artifact side alone satisfies the first
branch. That converts a real integrity check into decoration, which is the one outcome this
amendment must not produce.

Instead the single equality is replaced by **two named, dated, amendment-aware equalities**, both
pinned to literals copied from command output:

```julia
# spike/p13/net.jl, beside p13_consts_sha()
# The TWO sha256 values spike/p13/consts.jl has legitimately held. A third value is drift.
#   pre_amendment  -- c42cc8e + bfba6ac; the value every 13-10..13-16 artifact recorded
#   post_amendment -- after 13-D15-AMENDMENT.md (CHANGE A + CHANGE B), 2026-07-29
const P13_CONSTS_SHA = (pre_amendment = "<copied>", post_amendment = "<copied>")
```

and, at each of the three sites, in place of the single assertion:

```julia
@assert h.consts_sha     == P13_CONSTS_SHA.pre_amendment  "the net was not trained under the frozen pre-registration"
@assert p13_consts_sha() == P13_CONSTS_SHA.post_amendment "consts.jl is not at either sha this amendment authorises"
```

**This is STRICTLY STRONGER than what it replaces, and that is the point.** The original assertion
only required the two sides to *agree with each other* — a retrain plus an edit would have passed
it silently. The replacement pins **both sides to named literals**, so:

- a net trained under any other pre-registration fails (the guard's original purpose, kept intact);
- `consts.jl` at any sha other than the two this amendment authorises fails — **including a future
  undisclosed edit**, which the original form would have accepted;
- the two values *differing* becomes an **asserted, documented fact** instead of a silent
  tolerance, and `13-REPORT.md` §13 prints both.

The pre-amendment literal is captured **before any byte changes** (`13-17-PLAN.md` Task 1) from
live command output and cross-checked against `h.consts_sha` as the net itself reports it; the
post-amendment literal is captured after the edit, from the same command. Neither is retyped and
neither is copied out of a document — including out of this one.

**Nothing else is authorised.** No net is retrained, no weight changes, no gate is re-run,
`spike/p13/gate_report.jld2` and `spike/p13/alpha_report.jld2` stay byte-unchanged, and
`P13_ITERATION_ALLOWANCE` stays 1 of 1, UNSPENT.

### 7.4 Ordering: the guard fix lands BEFORE the 13-16 re-run

**Until the guard is fixed, nobody can use the suite to verify anything — including the corrected
13-16.** `runtests.jl` aborts at `:220` and eight Phase-13 test files never execute, so a re-run
verified only by per-file runs would be another per-file-only result. `13-17-PLAN.md` therefore
orders the guard fix and the sha re-derivation **before** the re-run, and requires the suite to be
observed reaching and executing the whole Phase-13 block first. That ordering is what makes the
re-run verifiable rather than merely repeated.

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
- It does **not** touch `spike/validation/p12_consts.jl`. That file's `P13_DEV_SEED` mirror is
  **CORRECT** — asserting seed disjointness requires naming the seeds — and the file is append-only
  Tier 1. CHANGE B's defect is on the guard side (§5B.1).
- It does **not** change any seed, bar, floor, band, ladder, tolerance or allowance. Neither
  CHANGE A nor CHANGE B moves a value that anything is scored against (§3, §5B.2).
- It does **not** fix the 28 other poisoned include guards the same sweep found. They are
  documented with a remedy each and left alone (§11).
- It does **not** retrain, re-run or re-read the gate, and it does **not** widen the
  `consts_sha` guard into a disjunction that would accept future drift (§7.3).
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

## 11. The 28 OTHER poisoned include guards — DOCUMENTED, NOT FIXED

The sweep that found CHANGE B's defect found **29** of them: **220 guards over 118 files** under
`spike/`, **29** whose sentinel is declared as a `const` by more than one file, in **eight
families**. This amendment fixes **one** family — the one whose file it already had open, for
CHANGE A, under this authorisation.

**The other 28 are recorded with a recommended remedy each and are NOT touched.** Editing another
phase's frozen file on a planner's own initiative is precisely the act this project's amendment
discipline exists to prevent, and a fix applied without an authorising document has the shape §0
spends four facts distinguishing this one from. The remedies split cleanly:

- **re-point at the source (M1-style)** where the owner file is editable — not a frozen
  pre-registration, or already opened by an authorised amendment;
- **isolated-module read (M2-style)** where the owner file cannot be edited, which is not a
  departure from the precedent but the only form of it available when the owner is frozen. This
  repo already uses it three times: `module _P11C` (`spike/p13/preconditions.jl:160`), `module _GC`
  (`spike/p13/consts.jl:96`) and `module GateV2` (`spike/validation/p11_consts.jl`).

**The full eight-family table, with a remedy per family and the owner-editability finding that
selects it, lives in `.planning/CONVENTIONS.md` (C-01 appendix).** It is there rather than here
because it is not a Phase-13 fact.

**The running record of individual observations belongs to Phase 12 and is accumulating in the
right place:** `.planning/phases/12-spatial-colocalization-map/deferred-items.md` logs this exact
collision as **DEF-12-03** (`p12_consts.jl:109` → `p13/consts.jl:100`, status OPEN), alongside
`13-09`, **DEF-12-01** (FIXED by `fb76b84`) and **DEF-12-02** (worked around per-runner).
**This amendment CLOSES DEF-12-03 rather than starting a parallel ledger**, and `13-17-PLAN.md`
records the closure by reference.

### 11.1 The standing rule, stated ONCE, for every future phase

> **Never guard an include on a seed or a prior bound, because those are exactly the names other
> phases are required to mirror.**

**23 of the 29 hits guard on a seed or a prior bound.** Four separate phases independently reached
for one, because R-5 *obliges* every Phase-N pre-registration to re-declare prior phases' seeds in
order to assert stream disjointness. The collision is structural, not careless — and that one
sentence would have prevented all 29.

It is recorded as **`.planning/CONVENTIONS.md` C-01**, a new durable file registered from
`.planning/PROJECT.md`'s `## Constraints` block. **That home was chosen deliberately:** every GSD
plan's `<context>` references `@.planning/PROJECT.md`, so a Phase-14 planner reaches it without
being told to. `STATE.md` was rejected as the home because it is a 900-line running log whose
blocker entries get closed — a rule buried in it is a record, not a convention; Phase 12's
`deferred-items.md` was rejected because Phase 14 will not open Phase 12's deferred items; and a
Phase-13 plan was rejected for the same reason. `CLAUDE.md`'s `## Conventions` section — which is
auto-loaded into every session and currently reads *"Conventions not yet established"* — is the
natural second home and is **recommended to the user**, not written by this planner.

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
| `220` guards / `118` files / `29` poisoned / `8` families / `23` on a seed-or-bound | the mechanical sweep of 2026-07-29, recorded in `.planning/STATE.md` | property of the source tree |
| `100` constants exclusively declared by `spike/p13/consts.jl`; `~104` declared inside the guard body | `grep -rn "const NAME\b" --include=*.jl .` per name, re-verified 2026-07-29 | property of the source tree |
| `22 pass / 1 fail / 114 error` at `runtests.jl:220` | observed full-suite run, recorded in `.planning/STATE.md` 2026-07-29 | prior observation |
| `0xC0FFEE` (`UInt32`) vs `0x0000_0000_00C0_FFEE` (`UInt64`) and the two siblings | read from the declaring lines named in §5B.5 | property of the source tree |
| Julia `1.12.6` | `spike/Manifest.toml:3` `julia_version` | property of the pinned environment |
| `5a4ea222…`, `100e97a3…`, `70fe66df…` (§7.2) | `13-REPORT.md` §13 as it reads today | prior observation — **superseded by this amendment; re-derived from live output by `13-17-PLAN.md`, never copied from here** |

**No number in this document was measured by this document.** Every Phase-13 figure quoted is a
prior result, labelled with the pair it was measured on, and the amended anchors are prior
measurements of the frozen simulator calibration — not new readings taken to justify the change.
