#!/usr/bin/env bash
# build.sh — fallback / CI REGRESSION gate for the ProteinCoLoc v2.0 manuscript (D-12, D-17).
# Mirrors build.ps1 (primary, Windows). Run from the repo root:
#   bash manuscript/build.sh
#
# --root manuscript makes root-absolute lookups of /colours.yaml and /figures/*.yaml (the
# per-figure YAML convention, D-14/D-15) resolve. --font-path bundles the Libertinus OTFs.
#
# Exit-code correctness (D-17, Pitfall 1): the `typst compile` line has NO pipe, so the shell
# reports typst's own exit status — a broken build cannot silently pass. `set -euo pipefail`
# aborts on the first failure. A non-fatal "did not converge / did not stabilize" warning from
# lib/tiered-bib.typ is EXPECTED and does NOT fail the gate: we read the compiler's exit code,
# we do NOT grep stderr for the word "warning".
#
# Typst version: the source template pinned 0.14; this machine runs 0.15.0. Default to the
# `typst` on PATH; override with the TYPST env var to match build.ps1's hook.
#
# Beyond compiling, this gate ASSERTS the Wave-2/3 content is present so later phases
# regression-check the skeleton (D-12): claim table (SC1), related-work matrix + four Tapqir
# axes (SC2), >=7 figure specs (SC3), the D-07 Supports->claims cross-check, and the D-14/D-16
# per-figure YAML machinery (figstyle loader, F1 YAML, F1 typ free of hard-coded geometry/hex).
# grep counts filter comment lines (^//) so a token in a header comment cannot self-satisfy an
# assertion nor falsely trip a negative check.

set -euo pipefail

TYPST="${TYPST:-typst}"

# Repo-relative, forward-slash paths from the repo root (Pitfall 3).
MAN="manuscript"
PDF="$MAN/main.pdf"

fail() { echo "GATE FAIL: $1" >&2; exit 1; }

# --- Step 1: compile (NO pipe — preserve typst's exit code; warnings are non-fatal) ----------
"$TYPST" compile --root "$MAN" --font-path "$MAN/fonts" "$MAN/main.typ" "$PDF"

# --- Step 2: content assertions (each exits non-zero on failure) ------------------------------
# Strip full-line Typst comments (// only — in Typst `#` is CODE, e.g. `#let`, `#table(`, not a
# comment) so a token appearing only in a header comment cannot satisfy (or, for the negative F1
# checks, falsely trip) an assertion.
noc() { grep -vE '^[[:space:]]*//' "$1"; }

# SC1 — compiles + the claim table renders.
[ -f "$PDF" ] || fail "SC1: $PDF was not produced by the compile step"
noc "$MAN/claims.typ" | grep -q "claim_table" \
  || fail "SC1: claim_table render function missing from claims.typ"

# SC2 — related-work matrix + all four Tapqir differentiation axes.
RW="$MAN/sections/related_work.typ"
noc "$RW" | grep -qi "amortiz"        || fail "SC2: 'amortization' axis missing from related_work.typ"
noc "$RW" | grep -qiE "sbc|calibrat"  || fail "SC2: 'SBC/calibration' axis missing from related_work.typ"
noc "$RW" | grep -qi "registration"   || fail "SC2: 'registration' axis missing from related_work.typ"
noc "$RW" | grep -qi "spatial"        || fail "SC2: 'spatial' axis missing from related_work.typ"
noc "$RW" | grep -q  "table("         || fail "SC2: related-work comparison matrix (table(...)) missing"

# SC3 — >=7 figure specs enumerating their producing phases.
SPEC_COUNT=$(noc "$MAN/figures/specs.typ" | grep -cE "SPEC F[0-9]")
[ "$SPEC_COUNT" -ge 7 ] || fail "SC3: only $SPEC_COUNT figure specs found (<7)"

# D-07 — every 'Supports: Cx' token in specs.typ resolves to a claim id in claims.typ.
for cid in $(noc "$MAN/figures/specs.typ" | grep -oE "Supports:[ ]*[C0-9, ]+" \
                | grep -oE "C[0-9]+" | sort -u); do
  grep -q "id: \"$cid\"" "$MAN/claims.typ" \
    || fail "D-07: figure spec references $cid but no matching 'id: \"$cid\"' in claims.typ"
done

# D-14/D-16 — the per-figure YAML machinery holds.
[ -f "$MAN/lib/figstyle.typ" ] || fail "D-14: lib/figstyle.typ missing"
grep -q "figstyle" "$MAN/lib/figstyle.typ" || fail "D-14: figstyle loader not defined in lib/figstyle.typ"
[ -f "$MAN/figures/f1_speedup.yaml" ] || fail "D-14: figures/f1_speedup.yaml missing"
grep -q 'figstyle("f1_speedup")' "$MAN/figures/f1_speedup.typ" \
  || fail "D-16: f1_speedup.typ does not load its YAML via figstyle(\"f1_speedup\")"
# F1 must be fully YAML-driven: no hard-coded geometry (cm) and no inline hex in code lines.
if noc "$MAN/figures/f1_speedup.typ" | grep -qE '[0-9]+(\.[0-9]+)?cm'; then
  fail "D-14: f1_speedup.typ contains hard-coded cm geometry (must come from the YAML)"
fi
if noc "$MAN/figures/f1_speedup.typ" | grep -qE 'rgb\("#'; then
  fail "D-14: f1_speedup.typ contains an inline hex colour (must come from the palette role)"
fi

echo "OK: DoD met — compiles + claim table + matrix + >=7 specs ($SPEC_COUNT) + YAML-driven F1"
