# Phase 1 — Root Baseline Record (ENV-01 / DEMO-02)

**Recorded:** 2026-06-26
**Plan:** 01-01

## Reconciliation decision

**Chosen option:** `commit-as-baseline`

The root `Project.toml` and `Manifest.toml` were modified vs the pre-planning HEAD
(`9055a13`): `Project.toml` gained a `[compat]` block pinning `GLMakie = "0.10.5"`
(3 insertions) and `Manifest.toml` was fully re-resolved (~1165 insertions / ~1073
deletions). The package author (user) confirmed these edits were **intentional** and
elected to freeze them as the new spike baseline rather than restore to HEAD or adopt a
non-committed snapshot.

## Baseline reference

| Property | Value |
|----------|-------|
| **Baseline ref (commit)** | `f581d95` (`f581d95dbc3596c5d9a064d6325ed70ca43a79c0`) |
| **Commit subject** | `chore(baseline): freeze root manifests as v2.0 spike baseline` |
| **Protected files** | `Project.toml`, `Manifest.toml`, `src/` |
| **Prior publication tags (unchanged)** | `v1.0.2-scirep`, `v2.0-baseline` |

## Untouched-root verify command (use in all later phases)

Run from the repo root; exits non-zero if any protected file drifted from the baseline:

```sh
git diff --quiet f581d95 -- Project.toml Manifest.toml src/ && echo "root untouched" || { echo "ROOT DRIFT vs baseline f581d95"; exit 1; }
```

Plan 01-04 (decoupling proof) and any future ENV-01 / DEMO-02 check MUST assert
byte-identity against `f581d95` using the command above.

## Notes

- `.planning/` artifacts are **not** protected by this baseline — only `Project.toml`,
  `Manifest.toml`, and `src/`. Planning-doc commits do not constitute "root drift".
- The spike lives entirely under `spike/`; nothing in this phase edits `src/` or the
  root manifests after this commit.
