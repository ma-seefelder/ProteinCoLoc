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
