// claims.typ — the claim spine as DATA + a render function (D-06, RESEARCH Pattern 1).
//
// `#let`-only (no bare top-level content) so `#import "claims.typ": claims, claim_table`
// pulls definitions and emits NOTHING at the import site (Pitfall 4). The rendered table is
// produced by main.typ via `#figure(claim_table(claims), ...) <tab:claims>`.
//
// The `claims` array is itself the diffable plain-source artifact; `claim_table` renders it.
// Schema per row (all string-valued for a clean table + easy diff):
//   id     — claim identifier (e.g. "C1")
//   claim  — one-line claim statement (carry honesty caveats in-line, e.g. ">100x (caveated)")
//   axis   — differentiator axis (Amortization | Calibration | Speed | ...)
//   phase  — supporting / owning phase(s) (e.g. "4,7")
//   fig    — figure/experiment ID that substantiates it (e.g. "F1")
//   status — "supported (measured)" | "supported (caveat)" | "reserved slot" (D-07 honesty)
//
// This is a compile-safe stub seeded with a SINGLE placeholder row (Plan 10-02 seeds the full
// headline claim set: amortized ms inference, >100x speedup caveated, SBC, OOD, BF, registration
// UQ, spatial map, cross-method, external validation).

#let claims = (
  (
    id: "C0",
    claim: "Placeholder — full claim set seeded in Plan 10-02",
    axis: "—",
    phase: "10",
    fig: "—",
    status: "reserved slot",
  ),
)

#let claim_table(claims) = table(
  columns: 6,
  table.header([ID], [Claim], [Axis], [Phase], [Fig], [Status]),
  ..claims.map(c => (c.id, c.claim, c.axis, c.phase, c.fig, c.status)).flatten(),
)
