# Phase 13 — Deferred / out-of-scope items

Discoveries made while executing Phase-13 plans that are **not** caused by Phase-13 code and are
therefore out of scope for the plan that found them (executor SCOPE BOUNDARY rule). Recorded so
they are visible, not fixed here.

---

## D-13-A — `spike/test/runtests.jl` aborts at `test_npe.jl`, BEFORE the Phase-13 block

**Found during:** plan 13-11 (wave 7), verification step.
**Status:** pre-existing, NOT a regression, NOT fixed.

The full spike suite exits 1 at `spike/test/runtests.jl:169` → `spike/test/test_npe.jl:230`:

```
SC3 (NPE-03): Test Failed at spike/test/test_npe.jl:230
  Expression: sr.median_speedup > SPEEDUP_GATE
   Evaluated: 68.34638733442779 > 100.0
Phase 4 — NPE Training + ADVI Benchmark + Ablation |  149     1    150
```

**Two things this adds to the existing record.**

1. **The measured speedup has drifted further down.** The STATE.md blocker
   (`[Phase 11 wave-3 post-merge gate, 2026-07-25]`) records **92.50×** and **83.97×** on two
   consecutive serial runs. This run measured **68.35×**. The gate bar `SPEEDUP_GATE = 100.0` is
   pre-registered and was NOT touched.

2. **The masking point is EARLIER than STATE.md states.** STATE.md attributes the loss of a
   full-suite green signal to `test_p13_correction.jl` at `runtests.jl:213`. In fact the suite
   never reaches line 213 — it aborts at line 169, ~44 lines earlier. Consequently the ENTIRE
   Phase-13 include block (lines ~200-215, eleven `test_p13_*` files) is masked in a full-suite
   run, including `test_p13_datagen.jl` added by plan 13-11.

**Consequence for Phase 13 verification.** Plan 13-11's acceptance criterion
"`julia --project=spike spike/test/runtests.jl` exits 0 across all phases" is **unmeetable for
reasons that predate the plan**. Per-file runs are the available gate signal, which is the
practice STATE.md already prescribes for the rest of Phase 13. Plan 13-11 verified
`julia --project=spike spike/test/test_p13_datagen.jl` → **706 pass / 706 total, exit 0**.

**Not fixed here because** `SPEEDUP_GATE` is a pre-registered threshold and re-deriving or
relaxing it is a user decision, not an executor call — the same discipline that kept
`P13_TAU`, `P11_LAMBDA_ABLATION_FACTOR` and the Phase-7 gate constants untouched. `test_npe.jl`
is also Phase-4 code that plan 13-11 does not own.

**Owed decision (user):** the two options already named in the STATE.md blocker stand —
(a) confirm the 7-vs-8 flow-marginal hypothesis and record it as a named limit scoped to the
research net, or (b) explicitly retire the legacy Phase-4 speedup gate as not meaningful against
an 8-column prior. A third question is now also open: whether the drift 92.50 → 83.97 → 68.35
across runs is machine-load noise or a real trend, which needs repeated measurement to answer.

---

## D-13-B — plan 13-11's Task-2 `<verify>` command is defective as written

**Found during:** plan 13-11, Task 2 verification.
**Status:** worked around in-run; the PLAN text still carries the defective command.

The plan's automated verify is:

```julia
h = load_three_way("spike/p13/three_way_net.jld2")
@assert haskey(h, "head_log_odds") || hasproperty(h, :head_log_odds)
```

`load_three_way` returns a `NamedTuple`, and `haskey(::NamedTuple, ::String)` **throws a
`MethodError`** rather than returning `false`, so `||` never reaches the `hasproperty` branch.
The command fails on a perfectly valid artifact. Swapping the operands
(`hasproperty(h, :head_log_odds) || haskey(h, "head_log_odds")`) makes it behave as intended and
it then passes. The artifact itself was never in question.
