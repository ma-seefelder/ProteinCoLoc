# Phase 12 — Deferred Items

Out-of-scope discoveries logged during execution. Nothing here is fixed by the plan that found
it; each names the plan that found it, why it is out of scope, and what the fix would touch.

---

## DEF-12-01 — `runtests.jl` aborts at `test_sbc.jl:41` on `UndefVarError: SBC_FIX_M`

**Found during:** 12-04, running the plan-level verification `julia --project=spike spike/test/runtests.jl`.
**Caused by:** NOT 12-04. The two files 12-04 touches (`spike/npe/p12_architecture.jl`,
`spike/test/test_p12_architecture.jl`) are not in the failing chain — see the proof below.

### The mechanism, confirmed executably

1. `spike/validation/p12_consts.jl:103` re-declares
   `const VAL_MASTER_SEED = 0x0000_0000_5BC0_FFEE` as one entry of its forbidden-seed
   inventory. That is a faithful copy of the `p11_consts.jl:107` pattern and the VALUE IS
   IDENTICAL to `spike/validation/consts.jl:76` (`0x5BC0FFEE`), so **no seed changed and no
   pre-registration was breached**. The collision is in the NAME, used as a guard sentinel.
2. `spike/validation/harness.jl:51` guards its include of the Phase-5 pre-registration on
   exactly that name: `isdefined(@__MODULE__, :VAL_MASTER_SEED) || include(".../consts.jl")`.
3. `spike/test/runtests.jl:174` runs the Phase-12 aggregator **FIRST** (12-01's deliberate
   structural choice), which loads `p12_consts.jl` into `Main`.
4. So by the time `test_sbc.jl` includes `sbc.jl → harness.jl`, `VAL_MASTER_SEED` is already
   defined, `consts.jl` is SKIPPED, and `SBC_M` / `SBC_FIX_M` / `SBC_FIX_L` / `SBC_FIX_BINS` and
   the BF/OOD constants never come into existence.
5. `test_sbc.jl:41` then throws `UndefVarError: SBC_FIX_M`, aborting every include after it.

Reproduced in one process (12-04 session, 2026-07-29):

```
before:            VAL_MASTER_SEED = false
after p12_consts:  VAL_MASTER_SEED = true    SBC_M = false
after harness:     SBC_M = false   SBC_FIX_M = false   load_frozen_model = true
```

`load_frozen_model` IS defined while `SBC_FIX_M` is not — i.e. `harness.jl` loaded and the
Phase-5 constants file specifically did not.

### Why it only appeared now

It was MASKED by the Phase-4 `SPEEDUP_GATE`. `test_npe.jl` is included at `runtests.jl:189`,
five lines before `test_sbc.jl`, and in the 12-01 / 12-02 / 12-03 sessions it threw first
(measured `median_speedup` 84.44 / 88.60 / 89.77 against a bar of 100). This run **Phase 4 SC3
(NPE-03) passed 9/9**, so execution reached line 194 for the first time and the latent collision
surfaced. The defect has been present since 12-01 wired the aggregator first.

### Why 12-04 does not fix it

- `spike/validation/p12_consts.jl` is **FROZEN Tier-1, append-only**. 12-04 opens no Tier-2
  sentinel, and the offending line is an existing declaration — changing it would be a
  pre-registration breach by construction (the file's own header says a diff that MODIFIES a
  line is the breach).
- `spike/test/runtests.jl` and `spike/validation/harness.jl` are outside 12-04's
  `files_modified`. The FIRST-position ordering is 12-01's asserted structural decision
  (`test_p12_consts.jl` testset 9 asserts it), so moving the include is a decision that plan
  owns, not this one.

### Candidate fixes, for whoever owns it

- Re-guard `harness.jl:51` on a sentinel that only Phase-5's `consts.jl` defines (`:SBC_M` is
  the obvious one, and is already the sentinel `consts.jl:41` guards its own body on). One-line,
  no pre-registration surface touched. **Preferred.**
- Or have `test_sbc.jl` guard on `:SBC_FIX_M` rather than `:sbc_ranks` (`test_sbc.jl:37`).
- Do NOT "fix" it by editing `p12_consts.jl` or by moving the Phase-12 include.

**Suggested owner:** 12-05 (`test_p12_decoupling.jl`) already owns the suite-integrity surface,
or a standalone `/gsd-quick`.

### RESOLVED 2026-07-29 by the orchestrator — commit `fb76b84`

The preferred fix was taken, and widened: **three** callers guarded on `:VAL_MASTER_SEED`, not one
— `spike/demo.jl:62`, `spike/validation/harness.jl:51`, `spike/validation/train_ratio.jl:72`. All
three now guard on `:SBC_M`, the name `validation/consts.jl` uniquely owns and already guards its
own block on (`consts.jl:38-42`, which documents that sentinel in prose).

Not a Phase-12 defect. FOUR files declare `VAL_MASTER_SEED`, all to the same value: the owner
`validation/consts.jl:76` plus `p11_consts.jl:87`, `p12_consts.jl:103` and `p13/consts.jl:109` —
each of the latter three *must* name it, because each lists it as a forbidden seed. Phase 11 and
Phase 13 planted the same landmine; Phase 12 only stepped on it first.

No seed, threshold, bar or pre-registered constant changed; no frozen consts file edited. Verified
empirically in both load orders — p12-first yields `SBC_M = 2000`, `SBC_FIX_M = 50`,
`load_frozen_model` defined; harness-alone yields the same and double-include is idempotent.
`VAL_MASTER_SEED` keeps its value (1539375086) and resolves to the owner's `UInt32`.
Regression confirmed cleared: `test_sbc.jl` runs **28/28** with the Phase-12 aggregator loaded first.

---

## DEF-12-02 — `p12_consts.jl` binds `:P11_DEV_SEED`, the sentinel `p11_consts.jl` guards its Tier-1 block on

**Found during:** 12-06, building the ε-ridge runner against the Phase-11 pool.
**Status:** WORKED AROUND per-runner, NOT structurally fixed. Live for every later plan.

Same class as DEF-12-01, opposite direction. `p11_consts.jl:79` guards its whole Tier-1 block on
`:P11_DEV_SEED`; `p12_consts.jl:108` legitimately binds that same name because R-5 requires Phase 12
to forbid the Phase-11 stream **by name**. So loading p12 first silently no-ops p11's Tier-1 block,
after which p11's Tier-2 block — guarded on a *different* sentinel, so it still runs — dies on
`UndefVarError: SC2_SPEARMAN_ATTENUATION`.

**Only `p11`-then-`p12` loads.** 12-06 ran both orders and fixed the order in its own runner header.

**This will bite 12-13, 12-15 and 12-17**, each of which reaches for a Phase-11 file. Any new
`run_p12_*` runner that includes both pre-registrations MUST include `p11_consts.jl` BEFORE
`p12_consts.jl`.

Not repaired structurally because both files are FROZEN Tier-1 and append-only, so neither the
sentinel nor the binding can be changed. The durable repair is the 13-09 precedent — read the
foreign consts file through a private module (`module _P11C`, cf. `module _GC` in
`spike/p13/consts.jl:96-98` and `module GateV2P12` in `p12_consts.jl`) — which a later plan may
adopt if the load-order convention proves too fragile to carry.

**General lesson, now observed three times in this repo** (13-09, DEF-12-01, DEF-12-02):
*"guarded include ⇒ order-free" is FALSE whenever two frozen consts files reserve the same name* —
and every Phase-N pre-registration is REQUIRED to reserve the prior phases' seed names, so this
collision is structural, not accidental. Guard each block on a sentinel the block alone owns.

---

## DEF-12-03 — the SAME collision against Phase 13: `runtests.jl` now dies with 114 `UndefVarError`s in `test_p13_consts.jl`

**Found during:** 12-07, running the plan-level verification `julia --project=spike spike/test/runtests.jl`.
**Caused by:** NOT 12-07. Neither file this plan touches (`spike/simulator/p12_prior.jl`,
`spike/test/test_p12_prior.jl`) appears anywhere in the failing chain — proof below.

### The mechanism

Exactly DEF-12-02's class, third direction:

1. `spike/p13/consts.jl:100` guards its **entire** Tier-1 block on `if !isdefined(@__MODULE__, :P13_DEV_SEED)`.
2. `spike/validation/p12_consts.jl:109` legitimately binds `const P13_DEV_SEED = 0x0000_0000_0B13_DE71`,
   because R-5 requires Phase 12 to forbid the Phase-13 stream **by name**. (The value is identical
   to `p13/consts.jl:188`; no seed changed and no pre-registration was breached — the collision is
   in the NAME, used as a sentinel.)
3. `runtests.jl:174` loads the Phase-12 aggregator FIRST, so `test_p12_consts.jl` puts
   `P13_DEV_SEED` into `Main` before anything Phase-13 runs.
4. `test_p13_consts.jl:42`'s `isdefined(@__MODULE__, :P13_DEV_SEED) || include(".../p13/consts.jl")`
   therefore SKIPS the include, and `P13_SALT`, `_p13_forbidden`, `P13_DECLARED_DEVIATIONS` and every
   other Phase-13 constant never comes into existence.
5. `test_p13_consts.jl` reports **22 passed, 1 failed, 114 errored** and throws, aborting
   `runtests.jl` (exit 1) at line 220.

Reproduced in one process with **neither 12-07 file loaded** (2026-07-29):

```
include("spike/test/test_p12_consts.jl")   # 12-01's file, alone
>>> after p12_consts: P13_DEV_SEED defined? true   P13_SALT defined? false
```

`P13_DEV_SEED` is defined while `P13_SALT` is not — i.e. the sentinel is set and the block that
would define the rest was skipped. Same signature as DEF-12-01's proof.

### Why it appeared only now

It was MASKED by the Phase-4 `SPEEDUP_GATE`, which threw at `test_npe.jl` and aborted everything
after it in the 12-01/12-02/12-03 sessions (measured `median_speedup` 84.44 / 88.60 / 89.77 against
a bar of 100). Commit `a494247` — *"test(04): PAUSE the NPE-03 wall-clock gate — not retire it, not
relax it"* — converted that assertion to `@test_broken`, so Phase 4 now reports `149 passed,
1 broken` instead of throwing, and `runtests.jl` reached the Phase-13 block for the **first time**.
The defect has been present since 12-01 committed `p12_consts.jl`.

### Why 12-07 does not fix it

- `spike/validation/p12_consts.jl` is **FROZEN Tier-1, append-only**; the offending line is an
  existing declaration, so editing it is a pre-registration breach by construction.
- `spike/test/test_p13_consts.jl` and `spike/p13/consts.jl` are Phase-13 files and **a Phase-13
  executor is live on this same branch**. Editing them from a Phase-12 plan is both out of scope
  and a concurrent-write hazard.
- 12-07's `files_modified` is two files, neither of which is in the chain.

### Candidate fixes, for whoever owns it

- Re-guard `test_p13_consts.jl:42` on a sentinel `spike/p13/consts.jl` **alone** owns —
  `:P13_SALT` is the obvious choice, exactly as DEF-12-01 was resolved by re-guarding three
  callers on `:SBC_M`. One line, no pre-registration surface touched. **Preferred.**
- Or have `spike/p13/consts.jl:100` guard on `:P13_SALT` instead of `:P13_DEV_SEED` (same
  one-line shape, but touches a frozen file's guard line — prefer the caller-side fix).
- Do NOT "fix" it by editing `p12_consts.jl` or by moving the Phase-12 include; `test_p12_consts.jl`
  testset 9 asserts the first-position ordering and 12-01 owns that decision.

**Suggested owner:** the Phase-13 executor (it owns those files and is live), or a standalone
`/gsd-quick` once Phase 13's wave settles.

**Fourth observation of the same class** (13-09, DEF-12-01, DEF-12-02, DEF-12-03). The prediction in
DEF-12-02's closing paragraph — that every Phase-N pre-registration is *required* to reserve prior
phases' seed names, so the collision is structural — is now confirmed against every phase
`p12_consts.jl` reserves a seed name for: Phase 5 (DEF-12-01), Phase 11 (DEF-12-02), Phase 13
(DEF-12-03). The durable repair is one convention, applied once: **guard each block on a sentinel
that block alone owns** — a salt or a bar, never a seed name, because seed names are exactly the
names other phases are obliged to re-declare.
