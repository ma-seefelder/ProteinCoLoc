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
