// figures/specs.typ — figure specification entries (D-10). Each spec reserves a figure's slot
// with its ID, working title, panel list, the phase that produces its data, the claim row(s) it
// supports, and a status placeholder — so later phases swap in the rendered figure without
// renumbering. Every spec (from Plan 04 onward) names its future figures/<name>.yaml + <name>.typ
// target so later phases clone the F1 pattern (D-16).
//
// This is a compile-safe stub — the top-level `= Figure Specifications` heading is emitted by
// main.typ. The full spec set + the worked F1 (speedup-vs-RMSE) figure land in Plan 04.
// No citekey/label references here yet (compile stays green, Pitfall 2).

// TODO(phase 10 plan 04): figure spec entries (F1–F7) + per-figure YAML targets; build the worked
// F1 speedup-vs-RMSE figure end-to-end (figures/f1_speedup.{yaml,typ} + f1_speedup_data/*.csv)
// using figstyle + colours.yaml roles, seeded with the honest Phase-4 numbers.
_Draft slot._ Figure specifications (F1–F7) are enumerated here in Plan 04, each naming its
per-figure YAML + `.typ` target.
