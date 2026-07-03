# Phase 6: Reproducible Demo + Go/No-Go Memo - Discussion Log

> **Audit trail only.** Not consumed by downstream agents. Decisions live in 06-CONTEXT.md.

**Date:** 2026-07-03
**Phase:** 06-reproducible-demo-go-no-go-memo
**Mode:** discuss (interactive)
**Areas discussed:** Memo verdict · Pre-registration honesty · demo.jl repro scope · Full-build-out call

## Note on process

The first two AskUserQuestion prompts timed out (60s, no response) and Claude auto-drafted
defaults. The user then reported never having seen the questions. The prompts were re-issued,
reached the user, and **all four areas were answered interactively.** The final decisions below
are the user's actual answers — two of which override the earlier auto-drafts.

## Questions & answers

| Area | Options offered | User's answer | vs draft |
|------|-----------------|---------------|----------|
| Memo verdict | Conditional Go / Clean Go / iterate-first / No-Go | **Clean Go** | **override** (draft was Conditional Go) |
| Pre-registration honesty | report both+name / both+independent confirm gate / post-iteration only | **Both + independent confirm gate** | **override** (draft was report-both only) |
| demo.jl repro scope | two-tier / full-only / fast-only | **Two-tier (fast + --full)** | matches draft |
| Full-build-out call | DAG P11∥P13 lead / spatial-first / productionize-only | **DAG, P11∥P13 lead** | matches draft |

## What the overrides changed

- **Clean Go (D-01):** the memo recommends full productionization now, not a Go gated behind a
  blocking calibration close-out. It must still be explicit that the literal pre-registered
  gates did not pass as written, justifying the Go by reading each residual failure to a
  characterized non-method cause (χ²-over-power, data-scale, clamped baseline).
- **Independent confirm gate (D-05):** the Clean Go's credibility is backed by a hard Phase-7
  ship-gate — a fresh-seed, re-pre-registered SBC/BF/OOD confirmation run that must reproduce
  the story before amortized inference ships into `src/`. This replaces the earlier
  "close-out blocks feature work" framing; feature work (the DAG) is not blocked, but shipping
  is.

## Deferred ideas captured

- Independent confirmation ship-gate (fresh-seed re-pre-registered run; optional larger-cache
  retrain + non-clamped BF baseline) → executed in Phase 7 before src/ integration.
- Re-enabling OOD PP channel in the reported OR-fusion → Phase 7 hardening.
- RxInfer independent ADVI cross-check → post-Go, paper nice-to-have (BACK-01).
