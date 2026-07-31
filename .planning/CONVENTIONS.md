# Project conventions — standing rules that apply to EVERY phase

Rules recorded here are **not phase-scoped**. They exist because the same defect was found more
than once, in more than one phase, by more than one author. A rule earns a place here only when it
has been observed to recur; anything observed once belongs in that phase's `deferred-items.md` or
in `STATE.md`.

This file is registered from `.planning/PROJECT.md` (`## Constraints`) so that every GSD plan,
whose `<context>` block references `@.planning/PROJECT.md`, reaches it.

---

## C-01 — NEVER guard an `include` on a seed or a prior bound

> **Never guard an include on a seed or a prior bound, because those are exactly the names other
> phases are required to mirror.**

That one sentence, applied from the beginning, would have prevented all 29 hits in the sweep
below.

### The idiom, and how it breaks

The house idiom for making an `include` idempotent is:

```julia
isdefined(@__MODULE__, :SENTINEL) || include(joinpath(@__DIR__, "target.jl"))   # caller side
if !isdefined(@__MODULE__, :SENTINEL)                                           # owner side
    # ... the whole body ...
end
```

The idiom is only sound if `:SENTINEL` is a name the target file **exclusively owns**. If any
other file can define that name first, the guard reads "already loaded" when nothing was loaded,
the `include` becomes a silent no-op, and the target's constants never exist. The failure surfaces
far away, as `UndefVarError` on an unrelated name, and — worst of all — **only in the suite**:
per-file runs pass, because in a single-file process nothing else has defined the sentinel yet.

### Why seeds are the worst possible sentinel *in this project specifically*

Requirement R-5 obliges every Phase-N pre-registration to **re-declare prior phases' seeds** in
order to assert stream disjointness. A seed name is therefore, by design, a name that several
files are *required* to declare. Guarding on one is guaranteeing the collision. Prior bounds
(`LAMBDA_MIN`, `MU_PRIOR`, …) are the same class for the same reason: they are re-declared so a
later phase can assert it stayed inside an earlier phase's range.

**The empirical shape of the defect: 23 of the 29 poisoned guards guard on a seed or a prior
bound.** Four separate phases independently reached for one. This is structural, not careless.

### What TO guard on

A name that is:

1. **exclusively owned** — verify with a fresh repo-wide `grep -rn "const NAME\b" --include=*.jl .`
   and confirm exactly one declaring file. Do not trust a list; re-run the grep;
2. **declared unconditionally inside the guarded body** — never a Tier-2 / appended / optional
   constant whose absence the file itself tolerates. A sometimes-absent sentinel is worse than the
   bug it was meant to fix;
3. **not a seed and not a prior bound** (rule C-01);
4. **about the file's IDENTITY, not a value another phase might quote** — a declared-deviations
   register, a pre-registration id, a schema marker. Not a design knob, not a gate statistic, not
   an architecture rule: each of those is a value a later phase may legitimately mirror in order to
   assert compatibility with it.

### The two remedies, and when each applies

| Remedy | Use when | Form |
|---|---|---|
| **Re-point at the source** | the owner file **can** be edited — it is not frozen, or an authorised amendment already opens it | change the owner-side wrapper AND every caller-side guard, in ONE edit, to a name satisfying (1)-(4) |
| **Isolated-module read** | the owner file **cannot** be edited — a frozen pre-registration with no amendment open | `module _XC; include(".../frozen.jl"); end`, then re-bind what is needed from `_XC.NAME`. The file is read, not modified, and no gate is run |

**A caller-only fix is provably useless when both ends use the poisoned name.** The caller would
correctly call `include`, and the owner-side wrapper would still see the foreign definition and
skip its entire body. Fix both ends or fix neither. `fb76b84` (`:SBC_M`) was a caller-only fix and
it was *correct there* only because `validation/consts.jl` already guarded its own body on the
name it owns.

**Never "fix" this by editing the mirroring file.** The mirror is correct — asserting seed
disjointness requires naming the seeds. `spike/validation/p12_consts.jl` is append-only Tier 1 and
must not be touched for this. **The defect is on the guard side, every time.**

### In-repo precedents

- Isolated-module read: `module _P11C` (`spike/p13/preconditions.jl:160`, whose comment states
  *"NO GATE IS RUN AND THE FILE IS NOT MODIFIED -- this is a read"*), `module _GC`
  (`spike/p13/consts.jl:96`), `module GateV2` (`spike/validation/p11_consts.jl`).
- Caller re-point: `fb76b84` — *"guard validation/consts.jl on :SBC_M, the name it actually owns"*.
- Both-ends re-point: `13-D15-AMENDMENT.md` §5B (`spike/p13/consts.jl:100` +
  `spike/test/test_p13_consts.jl:42` + six sibling callers), executed by `13-17-PLAN.md`.

### THE GENERALISATION: guard poisoning and name collisions are ONE root cause, not two

*Added 2026-07-31, authorised by the phase-12 orchestrator, after the second function-name collision
of the milestone. Scope of that authorisation: extend C-01 with this subsection only.*

Everything above is about **guard sentinels** colliding. The same milestone has now twice had
**function names** collide, in a way that produces a different symptom and needs the same fix. They
are two faces of one root cause, and treating them as separate problems means fixing each twice.

**The root cause: this project `include`s files into SHARED SCOPE.** Almost every spike file is
`include`d into `Main` (or into a single test module) rather than being read behind a namespace. In
that arrangement two different files may bind the same name — and **whichever loads last wins,
silently, with no error and no warning.** What the collision costs depends only on what kind of name
it is:

| Colliding name | Symptom | Loudness |
|---|---|---|
| a **guard sentinel** | the owner-side `if !isdefined(...)` body is SKIPPED entirely, so ~100 constants are never defined and the failure surfaces far away as `UndefVarError` | loud, but misattributed |
| a **function or helper** | the later definition **overwrites** the earlier one for every caller, including the tests written for the earlier one | **completely silent** — a green suite may be testing the wrong function |

The function case is the more dangerous of the two precisely because nothing fails. A test file that
believes it is exercising the owner's implementation may be exercising a caller's replacement.

**Instances so far (all three found in execution, none by review):**

| # | Kind | What happened |
|---|---|---|
| 1 | guards | the 2026-07-29 sweep: **29 poisoned guards over eight sentinel families**, 220 guards across 118 files (appendix below) |
| 2 | function | **12-12**: a bare `_res`-style test helper with an identical zero-positional signature silently overwrote its **Phase-13 namesake** |
| 3 | function | **12-14**: the plan specified `augment_mask!` as one of `spike/npe/train_p12_npe.jl`'s functions, but 12-04 already defines it at `spike/npe/p12_architecture.jl:190` — which that trainer `include`s. Defining it again would have overwritten the pre-registered augmentation, with include order picking the winner, and `test_p12_architecture.jl` would then have been testing the trainer's copy. Avoided by calling 12-04's function and never redefining it |

**THE ONE REMEDY, IN BOTH DIRECTIONS: read a foreign file through a PRIVATE MODULE, never into
shared scope.** The repo already does this in three places and they are the pattern to copy:
`module _P11C` (`spike/p13/preconditions.jl:160`), `module _GC` (`spike/p13/consts.jl:96`),
`module GateV2` (`spike/validation/p11_consts.jl:90`). A private module makes the import explicit —
you re-bind exactly the names you meant to take — so neither a sentinel nor a helper can be captured
by accident.

**Two rules that follow, and that every future phase inherits:**

1. **Before defining ANY function in a file that `include`s another, check whether that name already
   exists in the included tree.** `grep -rn "function <name>" spike/` costs seconds. A plan
   instructing you to define a function is not evidence that the function does not already exist —
   in instance 3 the plan explicitly specified one that did.
2. **Phase-prefix helpers you own** (`_p12train_…`, `_p12s1_…`), so that even in shared scope a new
   helper cannot capture a name another phase depends on. This is the cheap half of the remedy and
   costs nothing to apply universally.

**A plan that asks you to define a name that already exists is reporting a defect, not granting
permission.** Call the existing one; record the deviation.

---

## C-01 appendix — the mechanical sweep (2026-07-29): 220 guards, 118 files, 29 poisoned

Question asked of every `isdefined(..., :NAME) || include(...)` and `if !isdefined(..., :NAME)`
under `spike/`: **is `:NAME` declared as a `const` by more than one file?** Recorded in full in
`.planning/STATE.md` (Blockers/Concerns, 2026-07-29). Reproduced here with a remedy per family,
because STATE.md is a running log and this is a standing rule.

| # | Sentinel | Guards | Declared by | Owner-side block | Owner editable? | Remedy |
|---|---|---|---|---|---|---|
| 1 | `:P13_DEV_SEED` | **9** | `p13/consts.jl:188`, `validation/p12_consts.jl:109` | `p13/consts.jl:100` | **yes** — `13-D15-AMENDMENT.md` opens it | **RE-POINT AT SOURCE → `:P13_DECLARED_DEVIATIONS`. FIXED by `13-17-PLAN.md`.** Both ends plus six sibling callers (`labels.jl:73`, `net.jl:89`, `preconditions.jl:191`, `result.jl:79`, `tau_probe.jl:92`, `toy_gaussian.jl:79`) |
| 2 | `:LAMBDA_MIN` | **9** | `validation/p11_consts.jl`, `p13/preconditions.jl:171` (a re-bind of Phase 11's value) | `p11_consts.jl:79` (itself guarded on the poisoned `:P11_DEV_SEED`) | **no** — Phase-11 Tier-1, frozen, no amendment open | **ISOLATED-MODULE READ** at each of the 8 callers, or re-point them to a name `p11_consts.jl` alone owns (e.g. `SC2_SPEARMAN_ATTENUATION`). Not fixable at the source until Phase 11 opens an amendment |
| 3 | `:P11_DEV_SEED` | **4** | `p13/consts.jl:116`, `p11_consts.jl:144`, `p12_consts.jl:108` | `p11_consts.jl:79` | **no** — same as #2 | **ISOLATED-MODULE READ** — already the `_P11C` precedent (`p13/preconditions.jl:160`) and already worked around per-runner as DEF-12-02 |
| 4 | `:NPE_MASTER_SEED` | 1 | **six** files (`baseline/run_advi.jl`, `npe/train_npe.jl`, `p13/consts.jl`, `test/test_npe.jl`, `p11_consts.jl`, `p12_consts.jl`) | `npe/train_npe.jl:64` | **yes** — not a frozen pre-registration | **RE-POINT AT SOURCE** to a name `train_npe.jl` alone owns. Note the literals differ in TYPE (`0xC0FFEE` is `UInt32`; `0x0000_0000_00C0_FFEE` is `UInt64`) — see C-02 |
| 5 | `:ABL_REL_MARGIN` | 1 | `npe/ablation.jl`, `test/test_npe.jl` | `npe/ablation.jl:67` | **yes** | **RE-POINT AT SOURCE** |
| 6 | `:ABL_FOLD_CONSISTENCY` | 1 | `npe/ablation.jl`, `test/test_npe.jl` | `npe/ablation.jl:70` | **yes** | **RE-POINT AT SOURCE** |
| 7 | `:BENCH_THREADS` | 1 | `npe/run_thread_sweep.jl`, `test/test_npe.jl` | `npe/run_thread_sweep.jl:55` | **yes** | **RE-POINT AT SOURCE** |
| 8 | `:MU_PRIOR` | 1 | `simulator/calibration.jl`, `simulator/prior.jl` | caller `simulator/p12_prior.jl:68` | **yes**, but see caveat | **RE-POINT AT SOURCE.** CAVEAT: `spike/p13/consts.jl:797` freezes `P13_TAU_SIMULATOR_SHA = 78dc37f…` as `git log -1 -- spike/simulator`; any commit touching `spike/simulator/` moves that pointer, so this one needs its own disclosure |
| 9 | `:P11_IMSIZE_SET` | 1 | `p13/preconditions.jl:176`, `p11_consts.jl` | `p13/preconditions.jl:176` | **yes** | **RE-POINT AT SOURCE** (Phase-13-owned) |
| 10 | `:P13_REPO_ROOT` | 1 | `p13/real_images.jl`, `test/test_p13_result.jl` | `test/test_p13_result.jl:64` | **yes** | **RE-POINT AT SOURCE** (Phase-13-owned; deliberately NOT done by `13-17`, whose `real_images.jl` edit is comment-only) |

**Totals:** 9 + 9 + 4 + 1 + 3 + 1 + 1 + 1 = **29**. Family 1 is fixed by `13-17-PLAN.md`. **The
other 28 are DOCUMENTED and NOT fixed** — fixing another phase's frozen file on a planner's own
initiative is exactly the act this project's amendment discipline exists to prevent. Each is
recorded above with the remedy that phase's owner should apply when it next legitimately opens
the file.

**`p13/consts.jl:93` is a FALSE POSITIVE of the scan** — a commented illustration of the idiom,
not a live guard. The live owner-side wrapper is `:100`.

**BEFORE TRUSTING ANY SWEEP, RUN IT AGAINST A SUBSET WHOSE ANSWER YOU ALREADY KNOW.** An over-broad
filter produces exactly the false alarm that trains you to distrust sweeps. (2026-07-31: a
ceiling-sweep exclusion pattern skipping any path *containing* `p12_consts.jl` also silently
swallowed `test_p12_consts.jl`, falsely reporting five gating constants as unread — caught only
because constants the author had just written tests for came back "unread". The same lesson cost the
plan-checker's `n_low` regex three tries, fixed only once it was RUN against the forms it claimed to
forbid.)

**Accumulating record:** the same class is logged as `13-09`, `DEF-12-01` (FIXED, `fb76b84`),
`DEF-12-02` (worked around per-runner) and `DEF-12-03` (family 1, fixed by `13-17`) in
`.planning/phases/12-spatial-colocalization-map/deferred-items.md`. That file is the right place
for new observations; this file is the rule.

---

## C-02 — a hex literal's WIDTH is its type, and a guard fix can make that suddenly matter

Julia sizes an unsigned hex literal by its digit count: `0xC0FFEE` (6 digits) is a `UInt32`;
`0x0000_0000_00C0_FFEE` (16 digits) is a `UInt64`. `==` compares them equal; `===` does not.

While a poisoned guard is skipping a file's body, that file's `const` re-declarations never
execute, so any type disagreement with an already-loaded mirror is **invisible**. Repairing the
guard makes the body run, and every mirrored constant is then re-declared into the shared module
for the first time.

**So: before repairing an include guard, diff the literal TYPES of every name the newly-running
body declares against every other declaration of that name that may already be loaded.** Known
disagreements at the time of writing:

| Name | narrow form | wide form |
|---|---|---|
| `NPE_MASTER_SEED` | `0xC0FFEE` (`UInt32`) — `npe/train_npe.jl:65`, `test/test_npe.jl:81`, `baseline/run_advi.jl:74` | `0x0000_0000_00C0_FFEE` (`UInt64`) — `p13/consts.jl:108`, `p11_consts.jl:86`, `p12_consts.jl:102` |
| `VAL_MASTER_SEED` | `0x5BC0FFEE` (`UInt32`) — `validation/consts.jl:76` | `0x0000_0000_5BC0_FFEE` (`UInt64`) — the three pre-registrations |
| `VAL_FIX_SEED` | `0xF1F7ED` (`UInt32`) — `validation/consts.jl:79` | `0x0000_0000_00F1_F7ED` (`UInt64`) — the three pre-registrations |

The **values** agree; only the widths differ. This is a hazard to VERIFY, not a licence to change
a seed: **no seed literal may be edited to resolve it.** If a redefinition does fail, the remedy
is on the guard/loading side (isolated-module read, or load order), never on the seed.

---

## C-03 — `git diff` context lines are not changes; gate on `-U0` and on `^[+-]`

A negative acceptance criterion of the form

```
git diff HEAD -- <file> | grep -cE '<forbidden identifiers>'   # must be 0
```

**fails on correct work.** Default diff context is three lines, so editing line *N* prints lines
*N−3 … N+3* as context, and a forbidden identifier sitting anywhere in that window is counted even
though nothing about it changed. A criterion that fails a correct executor trains that executor to
override the one check that proves nothing forbidden moved — the worst available incentive.

**Correct form:**

```
git diff -U0 HEAD -- <file> | grep -E '^[+-]' | grep -v '^[+-][+-]' | grep -cE '<forbidden>'
```

- `-U0` removes context lines;
- `grep -E '^[+-]'` keeps only added/removed lines (hunk headers start with `@`);
- `grep -v '^[+-][+-]'` drops the `+++`/`---` file headers.

**And narrow the target further when the identifier can legitimately appear in a changed line for
an unrelated reason.** Gating on changed **declaration** lines only —
`… | grep -E '^[+-]\s*const ' | grep -cE '<forbidden>'` — is what distinguishes *"a bar moved"*
from *"a line that happens to mention a bar was touched"*.

Note also: `grep -c` exits **1** when the count is `0`, so a `… && …` chain around a
must-be-zero grep will report failure on success. Capture the count and compare it.

---

## C-04 — index-keyed generation has TWO faces: reproducible pools, NON-DISJOINT pools

*Added 2026-07-31, authorised by the phase-12 orchestrator, after a wave-8 finding.*

Every pool generator in this repository keys its per-sample RNG on the **global sample index**, in
the counter word of a counter-based Philox stream. **Verified at all three defining sites, not
inferred from family resemblance:**

| Generator | Constructor |
|---|---|
| Phase 11 | `Philox4x(UInt64, (P11_DEV_SEED ⊻ P11_DATAGEN_SALT, UInt64(idx)))` — `p11_generate.jl:113-114` |
| Phase 12 | `Philox4x(UInt64, (P12_DEV_SEED ⊻ P12_DATAGEN_SALT, UInt64(idx)))` — `p12_generate.jl:147-148` |
| Phase 13 | `Philox4x(UInt64, (P13_DEV_SEED ⊻ P13_SALT ⊻ P13_DATAGEN_COUNTER, UInt64(idx)))` — `p13/datagen.jl:181` |

This is not one phase's design choice. It is the counter-based Random123 pattern the project adopted
at the top level and that CLAUDE.md mandates (*"counter-based, reproducible seeding of the
simulator"*), so **every future phase that builds a pool inherits both faces of it.**

### The two faces

| Face | Consequence |
|---|---|
| **THE GUARANTEE** | A lost pool regenerates **byte-identically**. This is why rebuilding Phase 11's destroyed 54 MB pool was a *compute cost, not a re-seed*, and cost no iteration allowance. |
| **THE TRAP** | Pool builders generate indices `1:n`. So two pools at the **same configuration** share samples `1:min(n,m)` **byte-identically**. **A "separate scoring pool" is a SUPERSET of the training pool, never a held-out set.** |

**The Phase-11 story teaches only the guarantee.** A reader who knows only that half will conclude
that generating a second pool yields fresh samples. It yields an overlap. Both faces are the same
property.

### The rules that follow

1. **A held-out block must be CARVED OUT OF THE ONE POOL, never generated as a second pool.** There
   is no index-offset keyword on any of the three builders, and changing the configuration to force
   a different content hash makes it *a different experiment*, not a held-out set.
2. **EXCLUDE THE VALIDATION BLOCK TOO, not just the training block.** It drives early stopping inside
   `NeuralEstimators.train`, so it is contaminated for model selection — **a leak that passes every
   index check**, because those indices genuinely *are* disjoint from the training set. This is the
   half that survives the obvious test.
3. **Assert disjointness against the trainer's RECORDED index sets, not against re-derived
   arithmetic** — and then show the assertion FIRES on a deliberately overlapping fixture. An
   `isempty(intersect(...))` never shown non-empty is an assurance, not a test.
4. **Drawing fresh through the VALIDATION stream is structurally safe and is the preferred pattern
   for evaluation.** Pools ride `DEV_SEED ⊻ DATAGEN_SALT` keyed by index; the `harness.jl`
   draw-simulate path rides `DEV_SEED ⊻ SALT` keyed by counter. The two salts are asserted distinct
   at load time. An SBC or coverage number computed on training samples cannot arise on that path.

### In-repo instances

- **Phase 12 / 12-15** — the finding. `train_p12_npe`'s `pool_indices` keyword exists so the scoring
  block is carved from the one pool; the bundle records `train_indices` / `val_indices` for the
  disjointness assertion. `12-EXECUTOR-BRIEFING.md` §10 carries the concrete case.
- **Phase 12 / 12-11** — safe: pinning `r1` changes the *sample* at a given index, and each rung is
  scored within its own split, never across rungs.
- **Phase 12 / 12-16, 12-18, 12-19** — safe by rule 4: they never read cached pools at all.

---

## C-05 — a WORDING rule is checked on the RAW source; a STRUCTURAL ban on the COMMENT-STRIPPED source

**Authorised 2026-07-31 by the phase-12 orchestrator, after the rule was arrived at twice: once as a
recorded plan defect (12-02) and once as a live test failure in 12-16 that repeated it.**

Several plans in this project assert things about source text. Two DIFFERENT kinds of assertion get
made that way, and **they need opposite treatments of comments.** Applying one treatment to both is
what breaks them, and it breaks them in the direction that looks fine.

| kind of assertion | example | check against | why |
|---|---|---|---|
| **STRUCTURAL BAN** — "this construct does not appear" | `fit(ZScoreTransform` absent from a read-time file; `use_gpu = true` absent; no `GlobalMeanPool` in the Phase-12 net; `generate_p12_pool` absent from the trainer | **COMMENT-STRIPPED** source | A comment that merely *mentions* the banned construct would satisfy an `occursin` on the raw text. **A ban a comment can satisfy is not a ban.** Every one of these files legitimately discusses the thing it forbids — that discussion is the point — so the check must not see it. |
| **WORDING RULE** — "this sentence is present" | Pitfall 6's `no spatial borrowing, full nuisance and global borrowing`; `coverage-only win`; the R-1 `derived` label | **RAW** source | **A wording rule lives precisely in the prose, so stripping comments is exactly what deletes the thing being checked.** These rules exist because a phrase is a falsifiable claim about the method; the assertion pins the phrase where a reader will meet it. |

### The failure this prevents, in both directions

- **Ban checked on raw source** → passes on a file that only *talks about* the construct, or fails on
  a file that is compliant but honest about what it forbids. Silent either way.
- **Wording rule checked on stripped source** → fails on a compliant file (the rule is in a header
  comment and got stripped), or, worse, **passes for the wrong reason** if the phrase happens to
  survive in a docstring somewhere unrelated to where the rule belongs.

### And a corollary that bit twice: A REQUIRED LITERAL MUST NOT BE LINE-WRAPPED

An `occursin("...")` is a substring test, so a phrase broken across a newline **does not match** —
even though a human reads the file as containing it. Where a file is required to carry an exact
phrase, write that phrase **on one line**, unwrapped, even if it makes the line long or forces an
indented display block. `12-02-SUMMARY.md` records "a line-wrap that broke an acceptance literal" as
a plan defect; 12-16 reproduced it independently in `p12_coverage.jl`, where the Pitfall-6 wording
survived unwrapped in the file header and was wrapped in a docstring.

**The fix in both cases was to correct the mechanism, never to loosen the check.** A wording
assertion that gets relaxed until it passes has stopped asserting the wording.

### In-repo instances

- **Phase 12 / 12-02** — the first instance, recorded among that plan's line-wrap and stale-citation
  defects.
- **Phase 12 / 12-16** — `spike/test/test_p12_coverage.jl` testset 8 applies the split explicitly:
  `raw` for the two wording rules, `_p12cov_tstrip(raw)` for `fit(ZScoreTransform` and
  `use_gpu = true`. `spike/validation/p12_coverage.jl` carries the Pitfall-6 phrase unwrapped in
  BOTH its header comment and an indented docstring block, so the rule survives either treatment.
- **Existing consumers to read this way, not to change:** `test_p12_consts.jl` and
  `test_p12_decoupling.jl` both define `_strip_comment_lines` and use it for STRUCTURAL bans, which
  is the correct side of this rule and needs no edit.
