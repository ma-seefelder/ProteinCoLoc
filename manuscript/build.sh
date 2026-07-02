#!/usr/bin/env bash
# build.sh — fallback / CI compile gate for the ProteinCoLoc v2.0 manuscript (D-17).
# Mirrors build.ps1 (primary, Windows). Run from the repo root:
#   bash manuscript/build.sh
#
# --root manuscript makes root-absolute lookups of /colours.yaml and /figures/*.yaml (the
# per-figure YAML convention, D-14/D-15) resolve. --font-path bundles the Libertinus OTFs.
#
# Exit-code correctness (D-17, Pitfall 1): the `typst compile` line has NO pipe, so the shell
# reports typst's own exit status — a broken build cannot silently pass. `set -euo pipefail`
# aborts on the first failure.
#
# Typst version: the source template pinned 0.14; this machine runs 0.15.0. Default to the
# `typst` on PATH; override with the TYPST env var to match build.ps1's hook.
#
# Content assertions (grep for rendered claim table / related-work matrix, etc.) are added in
# Plan 05 — they would fail against Wave-1 stubs, so this gate is compile-only for now.

set -euo pipefail

TYPST="${TYPST:-typst}"

# Forward-slash, repo-relative paths from the repo root (Pitfall 3).
"$TYPST" compile --root manuscript --font-path manuscript/fonts "manuscript/main.typ" "manuscript/main.pdf"

echo "OK: manuscript/main.pdf built (exit 0)"
