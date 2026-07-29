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
