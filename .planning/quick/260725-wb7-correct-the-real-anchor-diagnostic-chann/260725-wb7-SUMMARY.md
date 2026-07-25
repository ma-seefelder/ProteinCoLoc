---
phase: quick-260725-wb7
plan: 01
subsystem: spike/simulator (SIM-02 induced-μ calibration diagnostic)
tags: [calibration, diagnostic, honesty-fix, d-16, real-anchor, channel-pair, otsu-mask]
requires:
  - test/runtests.jl:105 (channel-name ground truth ["blue", "green", "red"])
  - src/LoadImages.jl::_calculate_mask / _apply_mask! (called read-only via contract.jl)
  - spike/contract.jl::build_mci / patch_summary / induced_mu / load_tiff (unchanged)
provides:
  - "ANCHOR_CHANNEL_PAIR = (2, 3) / ANCHOR_CHANNEL_NAMES: the measured pair is a named, defaulted parameter"
  - "_real_anchor(; pair): masked AND unmasked induced μ + surviving-patch count per fixture, for any pair"
  - "calibrate(; anchor_pair) returns anchor_pair so the provenance can state which channels produced a number"
  - "corrected real-anchor provenance (pair + channel map + SUPERSEDED note) in the frozen ghat.jl header"
affects:
  - spike/NOTES.md §3 (still quotes the c1/c2 counterstain figures — out of scope, see Follow-ups)
  - spike/test/test_simulator.jl:225 (same — out of scope)
  - .planning/phases/02-.../02-03-SUMMARY.md, 13-RESEARCH.md, 13-PATTERNS.md, 13-15-PLAN.md, 13-01-PLAN.md (same — orchestrator's job)
tech-stack:
  added: []
  patterns:
    - "the channel pair behind a published number is a named constant, threaded into the return value and rendered into the provenance"
    - "diagnostic reports BOTH estimators rather than silently choosing one (carried over from vl8)"
key-files:
  created: []
  modified:
    - spike/simulator/calibration.jl
    - spike/simulator/ghat.jl
decisions:
  - "the real anchor measures channels 2/3 (green/red) — the colocalization pair; c1 is the DAPI/Hoechst nuclear counterstain, so every c1 pair measured counterstain-vs-protein overlap, not colocalization"
  - "the c1/c2 figures (positive 0.3292 / negative 0.2481) are SUPERSEDED, not merely updated"
  - "neg_reachable flipping false → true is a predicate artifact of a ≈0 reading (−0.0338), NOT evidence; the vl8 WITHDRAWAL of the PRIOR-ONLY claim STANDS and is not inverted"
  - "the 'MASKED is biased NEGATIVE' framing is downgraded to a range-restriction caveat; the large negative readings were predominantly the wrong-channel artifact"
  - "ghat.jl's anchor block is hand-patched again, not regenerated — re-running calibration.jl would re-fit the frozen knots"
metrics:
  duration: ~20 min
  completed: 2026-07-25
  tasks: 3
  files: 2
---

# Quick 260725-wb7: Correct the Real-Anchor Diagnostic Channel Pair Summary

The SIM-02 real anchor now measures the colocalization pair (c2/c3 = green/red) instead of the
DAPI counterstain pair (c1/c2), the pair is a named parameter recorded in both the emitted
provenance and the printed report together with the `test/runtests.jl:105` channel map, and the two
vl8 sentences that the corrected measurement falsified are fixed — without reinstating any withdrawn
claim.

## What Was Measured

All three pairs, via the package's own `load_tiff` (D-16 — faithful load, not an RGB→luminance
reduction), on the two fixtures in `test/test_images/`. The c1/c2 and c2/c3 rows were reproduced by
this task's Task-1 verify run; the c1/c3 row is the planner's pre-dispatch probe (not re-run here).

| pair             | condition | unmasked μ | patches | masked μ | patches |
|------------------|-----------|-----------:|--------:|---------:|--------:|
| c1-c2 blue/green | positive  | 0.3292     | 64/64   | −0.0729  | 21/64   |
| c1-c2 blue/green | negative  | 0.2481     | 64/64   |  0.1149  | 20/64   |
| c1-c3 blue/red   | positive  | 0.2350     | 64/64   | −0.0925  | 20/64   |
| c1-c3 blue/red   | negative  | 0.3756     | 64/64   | −0.1568  | 44/64   |
| **c2-c3 green/red** | **positive** | **0.4603** | **64/64** | **+0.8238** | **22/64** |
| **c2-c3 green/red** | **negative** | **0.3815** | **64/64** | **−0.0338** | **20/64** |

Exact values observed in the Task-1 verify run for the shipped default pair `(2, 3)`:
`positive = (unmasked_mu = 0.4602657619986303, unmasked_n = 64, masked_mu = 0.8238196296329328, masked_n = 22)`,
`negative = (unmasked_mu = 0.3814745646648433, unmasked_n = 64, masked_mu = -0.0337877825290027, masked_n = 20)`.

The same run with `pair = (1, 2)` reproduced the vl8 numbers exactly
(`0.329163298494035 / -0.07293570883861245 / n=21`), proving the parameter really drives the load
rather than a coincidental re-derivation.

**Why c1 was the wrong channel.** `test/runtests.jl:105` records the fixture channels as
`["blue", "green", "red"]`. Channel 1 (blue) is the DAPI/Hoechst nuclear counterstain, not a target
protein — the user, a domain expert, confirmed this directly. Correlating a counterstain against a
target protein yields noise, which is exactly what every c1 pair shows (−0.16 to +0.11 masked). The
colocalization pair is green/red = channels 2 and 3.

**On the correct pair, masked separates the controls as their labels promise:** positive +0.8238,
negative −0.0338. Unmasked, even the correct pair barely separates (0.4603 vs 0.3815) — the
positive background bias compresses the contrast. So the correct pair *and* masking are both
required to see the separation.

## `neg_reachable` — flipped to `true`, and why that changes no verdict

`neg_reachable` is defined as `isfinite(masked_mu) && masked_mu < 0.0` on the negative fixture. It
was `false` (c1/c2 masked negative = +0.1149); on the correct pair it is **`true`**
(−0.0338 < 0).

**This is a predicate artifact, not evidence.** −0.0338 is approximately zero, not convincingly
negative, and n = 2 fixtures cannot settle whether a negative induced-μ tail is physically
reachable in either direction. Both the code comment and the emitted provenance state this
explicitly. Concretely, the following remain true in the emitted text:

- negative-tail reachability is **STILL NOT ESTABLISHED** either way;
- the original unqualified `"physically reachable = false ⇒ PRIOR-ONLY"` reading stays
  **WITHDRAWN** (vl8) — it is neither reinstated nor inverted;
- the clean +0.82 / −0.03 control separation is explicitly recorded as **not** evidence for the
  opposite claim;
- the SIM-02 consistency claim stays scoped to `[GHAT_MU_MIN, GHAT_MU_MAX]`, unchanged and not
  dependent on this anchor.

## What Was Built

**Task 1 — `spike/simulator/calibration.jl` (commit `d51e7b4`)**
`ANCHOR_CHANNEL_PAIR = (2, 3)` and `ANCHOR_CHANNEL_NAMES = ("blue (DAPI/Hoechst nuclear
counterstain)", "green", "red")` are module-level constants with a comment recording that c1 is a
counterstain and why the pair is now explicit; `_channel_label(i)` renders `"c2 = green"`-style
labels. `_real_anchor(; pair::Tuple{Int,Int} = ANCHOR_CHANNEL_PAIR)` interpolates the pair into the
fixture filenames (`$(cond)_c$(pair[1]).tif`) and into the MCI name (`positive_c2c3`), so even the
`MultiChannelImage` carries the pair. Everything vl8 established is untouched: the unmasked
measurement is still taken **before** `_apply_mask!` (which mutates `mci.data` in place — the
load-bearing-ordering comment is verbatim), both measurements still go through `_anchor_measure`,
and `_anchor_measure` / `_anchor_na` are unchanged. `calibrate(; anchor_pair = ANCHOR_CHANNEL_PAIR)`
passes it to `_real_anchor` inside the existing `try` and returns `anchor_pair` in the NamedTuple;
the `catch` fallback is unchanged. `neg_reachable` keeps its definition and reads the masked
negative μ — only its comment changed.

**Task 2 — emit template + printed report (commit `15b85ac`)**
Both now open the anchor block with the pair (`c2/c3 = green/red (ANCHOR_CHANNEL_PAIR = (2, 3))`,
rendered from `r.anchor_pair` and `ANCHOR_CHANNEL_NAMES`, not hardcoded), the full
`test/runtests.jl:105` channel map, and the counterstain explanation ending in the SUPERSEDED note.
Three prose corrections:
- **deleted** "the two estimators disagree in SIGN and in ORDER" — on c2/c3 they **agree** in order
  (positive > negative under both), and masked separates the controls far more sharply than
  unmasked. The text now says that.
- **softened** the mask bias: the range-restriction caveat stays in full (independent per-channel
  Otsu + `_exclude_zero` keeps only pixels bright in BOTH, a selection on both variables that can
  attenuate or distort the correlation — the known thresholded-Pearson weakness behind
  Manders/Costes — leaving only ~20–22 of 64 patches above the ≥15-survivor floor), but the text
  now states plainly that "MASKED is biased NEGATIVE" **overstated** it: the masked positive control
  reads +0.8238, so the selection does not prevent detecting strong colocalization, and the large
  negative readings were predominantly the wrong-channel (counterstain) artifact.
- **kept and extended** the unmasked positive-bias sentence: on this pair the shared dark background
  also *compresses* the control contrast, which is why unmasked barely separates them.

The four μ values and four counts still interpolate `r.anchor.*`, so they became the corrected
numbers automatically. All other evidence lines (generator, sweep, induced-μ range, monotonicity,
SIM-02 match, size-invariance, σ/τ, ν) and the `GHAT_*` / `ghat()` emission are byte-identical.

**Task 3 — `spike/simulator/ghat.jl` header (commit `78dc37f`)**
Comment-only hand patch. The corrected block is copied **verbatim from the rendered Task-2
template output** (captured from the temp-path render), so a future regeneration is a no-op on this
block. A second dated `HAND-PATCH PROVENANCE (2026-07-25, second entry)` note records the pair
change, the `test/runtests.jl:105` justification, and the unchanged reason for hand-patching
(re-running `calibration.jl` re-runs a stochastic sweep and re-fits the frozen knots).

## Verification (all run, output observed)

| Check | Command | Observed |
|-------|---------|----------|
| Task 1 | `julia --project=spike -e 'include("spike/simulator/calibration.jl"); …'` | printed both pairs' NamedTuples (values quoted above); all asserts held — `ANCHOR_CHANNEL_PAIR == (2,3)`, 64/64 unmasked, 22 and 20 masked, all four μ within 5e-4 of the reference, and `pair=(1,2)` reproduced 0.3292 / −0.0729 / n=21; **`REAL-ANCHOR-PAIR OK`** |
| Task 2 (template) | probe rendering `emit_ghat(r; path=tempdir()/ghat_probe_wb7.jl)` with a stub `r` | rendered block contains `UNMASKED`, `MASKED`, `c2/c3`, `green`, `red`, `counterstain`, `test/runtests.jl:105`, `SUPERSEDED`, `NOT ESTABLISHED`, `WITHDRAWN`, `Manders`, all four rounded μ, `64/64`, `22/64`, `20/64`; `disagree in SIGN` absent; **`EMIT-TEMPLATE OK`** |
| Task 2 (report) | probe evaluating the 8 report `println` source lines (430–437) against a stub `r` | all eight rendered (full text observed); same substrings present incl. `ANCHOR_CHANNEL_PAIR = (2, 3)`; **`REPORT-LINES OK`** |
| Task 3 (diff) | `git diff -U0 -- spike/simulator/ghat.jl \| grep -E '^[+-][^+-]' \| grep -vE '^[+-]\s*#'` | **empty** — every added/removed line is a comment |
| Task 3 (load) | `julia --project=spike -e 'include("spike/simulator/ghat.jl"); …'` | `GHAT_MU_MIN == -0.67976`, `GHAT_MU_MAX == 0.847149`, 25 + 25 knots, `ghat(0.0) = -0.12004360606355068` to 1e-12; header substrings all present; **`GHAT-HEADER OK`** |
| Decoupling | `git status --porcelain src/ spike/contract.jl` | **empty** |
| Scope | `git diff --name-only HEAD~3 HEAD` | exactly `spike/simulator/calibration.jl`, `spike/simulator/ghat.jl` |
| Deletions | `git diff --diff-filter=D --name-only HEAD~3 HEAD` | **empty** |
| Planning artifacts | `git status --porcelain .planning/phases/` | only the two pre-existing untracked `*-PATTERNS.md` files present before dispatch — not touched by this plan |

`calibrate()` itself was deliberately **not** run (it would re-run the stochastic sweep and re-fit
the frozen knots). Every verify `include`s the file; the `PROGRAM_FILE` guard keeps the script entry
inert. `calibrate()`'s new `anchor_pair` keyword and return field are therefore verified at
parse/compile level (the file loads) and by the report probe consuming `r.anchor_pair`, not by a
full sweep — recorded here rather than implied.

## Deviations from Plan

**None affecting scope or behavior.** Two implementation details worth recording:

1. **[Rule 3 — blocking] The caveat numbers are interpolated, not hardcoded.** The plan's Task-2
   text quotes literal values ("+0.8238 vs −0.0338", "0.4603 vs 0.3815") in the caveat prose.
   Hardcoding them would drift silently if the anchor were ever re-measured, so they are
   interpolated from `r.anchor.*` — the rendered output is byte-for-byte what the plan specified
   (verified above), with the drift hazard removed. The one genuinely historical pair
   (0.3292 / 0.2481) *is* a literal, correctly so: it names a superseded past measurement.
2. **Added verification beyond the plan.** The plan's Task-2 verify covered only the emit template.
   Following vl8's precedent, a second probe evaluates the eight report `println` source lines
   directly against a stub `r`, so the "mirror the same content in the printed report" requirement
   is backed by observed output rather than by inspection.

Also updated (in-scope, same file, same claim): the `_real_anchor` header comment in
`calibration.jl` carried its own copy of the "they disagree in SIGN and in ORDER" sentence. It would
have been left factually wrong by an emit-template-only edit, so it was corrected too.

## Assumption Drift (advisory)

None material. Every number the plan predicted for the (2, 3) pair — and every number it predicted
for the (1, 2) control run — was reproduced exactly, so no re-numbering of the ghat.jl header was
needed.

## Out-of-Scope Follow-ups (identified, NOT edited)

Seven artifacts still quote the c1/c2 counterstain figures (0.3292 / 0.2481) and/or the withdrawn
PRIOR-ONLY reading:

1. **`spike/NOTES.md` §3 "Real anchor (D-16)"** — the μ figures plus "the negative μ-prior tail is
   NOT physically reachable" / "real fluorescence data never showed anti-correlation".
2. **`spike/test/test_simulator.jl:225`** — comment: "the negative μ-prior tail past `GHAT_MU_MIN`
   is prior-only (real anchor > 0)".
3. **`.planning/phases/02-forward-simulator-summary-contract/02-03-SUMMARY.md`** — records
   "positive μ=0.329, negative μ=0.248 (both >0)" as a frozen decision and headline claim.
4. **`.planning/phases/13-.../13-RESEARCH.md`**
5. **`.planning/phases/13-.../13-PATTERNS.md`**
6. **`.planning/phases/13-.../13-15-PLAN.md`**
7. **`.planning/phases/13-.../13-01-PLAN.md`**

Items 3–7 are `.planning/` artifacts: the Phase 11 and Phase 13 plan corrections are the
**orchestrator's separate job** and were explicitly out of scope here (hard constraint 1). Items 1–2
are spike source/test files, out of `files_modified` for this plan; item 3 is a historical phase
summary, so correcting it is a records question (amend vs. append a correction note), not a code
change.

## Known Stubs

None. The probe scripts live in the session scratchpad, not in the repo; the temp-rendered
`ghat_probe_wb7.jl` was written to `tempdir()`, never to `spike/simulator/ghat.jl`.

## Threat Flags

None. No new network endpoint, auth path, file-access pattern, or schema change. Against the plan's
threat register: T-wb7-01 mitigated (ghat.jl diff proven comment-only, all four `GHAT_*` constants
and `ghat(0.0)` re-checked to 1e-12); T-wb7-02 mitigated (`src/`, `spike/contract.jl` and
`src/amortized/summary.jl` untouched — `git status --porcelain` empty; the Otsu mask stays confined
to `_real_anchor()`, no masking added to any training or inference path); T-wb7-03 mitigated (the
pair is a named constant, threaded into `calibrate()`'s return, and rendered into both the emitted
provenance and the printed report alongside the channel map); T-wb7-04 mitigated (NOT ESTABLISHED +
WITHDRAWN asserted present after `neg_reachable` flipped to `true`); T-wb7-05 accepted (no
concurrent-executor file was touched or staged).

## Commits

| Task | Commit | Message |
|------|--------|---------|
| 1 | `d51e7b4` | `fix(wb7-01): measure the real anchor on the colocalization pair c2/c3, not the DAPI counterstain` |
| 2 | `15b85ac` | `docs(wb7-02): record the anchor channel pair and correct the two wrong-pair claims in the emit template and report` |
| 3 | `78dc37f` | `docs(wb7-03): correct the frozen ghat.jl real-anchor header to the c2/c3 pair (comment-only)` |

## Self-Check: PASSED

- Files exist: `spike/simulator/calibration.jl`, `spike/simulator/ghat.jl`, this SUMMARY.
- Commits exist: `d51e7b4`, `15b85ac`, `78dc37f`.
- `git diff --name-only HEAD~3 HEAD` = exactly the two files in `files_modified`; no deletions.
- `git status --porcelain src/ spike/contract.jl spike/simulator/` empty; the only other dirty paths
  are the pre-existing `.planning/HANDOFF.json` deletion, two untracked `*-PATTERNS.md` files and
  `decisions/` — none staged by this plan.
