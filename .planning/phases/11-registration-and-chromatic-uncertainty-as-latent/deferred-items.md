# Phase 11 — deferred items

Out-of-scope discoveries logged during execution. Not fixed here.

## From plan 11-02 (chromatic ε as an 8th θ column)

- **`scripts/diag_bounded_theta_dev.jl:62-65` still hardcodes a θ arity of 7**
  (`Matrix{Int}(undef, M, 7)`, `Matrix{Float64}(undef, M, 7)`, `zeros(Float64, 7)`,
  `Matrix{Float64}(undef, M, 7)`).
  It is a developer diagnostic script, not part of the package or of `Pkg.test()`, so nothing
  currently breaks — but it would under-count columns if re-run against the extended prior.
  Not fixed here: `scripts/` is outside plan 11-02's file scope, and touching it would widen the
  D-11/D-12 provenance commit for no correctness gain.

- **`test/gate/sbc.jl` `SBC_PARAM_LABELS` is still 8 labels (7 θ + Δρ)** and
  `sbc_ranks_and_spread` still builds an `M×8` rank table over `for p in 1:7`.
  This is deliberate and matches `11-RESEARCH.md` §A4 row 25: the src-side ship gate is not
  re-run in Phase 11 (D-01), the frozen pre-registration must not be edited, and the loop is
  arity-safe against a wider posterior (it reads the first 7 rows of an 8-row draw matrix).
  Revisit only in a phase that actually re-runs the src gate.

## From plan 11-04 (research-net scaffold)

- **`spike/test/test_simulator.jl:199-200` still pins the PRE-11-02 seven-field θ**, so
  `julia --project=spike spike/test/runtests.jl` currently exits 1:

  ```
  sample_prior: 7-field θ, in-range ρ_true, deterministic (D-14): Test Failed at
    spike/test/test_simulator.jl:200
    Expression: keys(θ) == (:ρ_true, :spillover, :autofluorescence, :label_efficiency,
                            :shift_dx, :shift_dy, :noise)
    Evaluated:  (:ρ_true, …, :noise, :chromatic_eps) == (:ρ_true, …, :noise)
  ```

  **Pre-existing, and not caused by plan 11-04.** `spike/simulator/prior.jl` gained
  `chromatic_eps` in `ca02b0e` (plan 11-02); `spike/test/test_simulator.jl` was last touched in
  `02c3bee` (Phase 2) and was never updated to the post-edit arity the way its sibling
  `spike/test/test_p11_forced_theta.jl` was (11-02 deviation 3). The failure is present at plan
  11-04's worktree base `2aae21f` with no 11-04 file loaded.

  **Blast radius is larger than one assertion:** `test_simulator.jl` THROWS, so
  `runtests.jl` aborts there and the five later includes (`test_data_pipeline.jl`,
  `test_npe.jl`, `test_sbc.jl`, `test_bf.jl`, `test_ood.jl`, `test_comparator.jl`) do not run
  at all.

  Not fixed here: `spike/test/test_simulator.jl` is outside plan 11-04's declared file set, and
  plan 11-03 was executing in parallel against the same phase, so an undeclared edit risked a
  merge collision on a file neither plan owns. The fix itself is mechanical and mirrors what
  11-02 already did for `test_p11_forced_theta.jl` — extend the expected tuple with
  `:chromatic_eps`, keep the assertion a tripwire rather than a restatement, and update the
  testset name, which still says "7-field θ".

  The Phase-11 clause (i) dependency gate that plan 11-04 added to `runtests.jl` is
  **unaffected and green** (`D-04 CPU-only (no CUDA dependency, none loaded)`: 23/23 pass,
  up from 21 assertions).
