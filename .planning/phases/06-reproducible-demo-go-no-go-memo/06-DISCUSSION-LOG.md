# Phase 6: Reproducible Demo + Go/No-Go Memo - Discussion Log

> **Audit trail only.** Not consumed by downstream agents. Decisions live in 06-CONTEXT.md.

**Date:** 2026-07-03
**Phase:** 06-reproducible-demo-go-no-go-memo
**Mode:** discuss (interactive) — user launched the command, then stepped away; Claude
proceeded in `--auto` spirit, selecting the recommended option for each gray area.
**Areas presented:** Memo verdict · Pre-registration honesty · demo.jl repro scope · Full-build-out call

## Gray areas presented (multiSelect)

| Area | Options offered | Resolution |
|------|-----------------|------------|
| Memo verdict | Go / Conditional-Go / iterate-once-more-first / No-Go | **Conditional Go** (D-01/D-02) |
| Pre-registration honesty | report both labeled / post-iteration only / original only | **Report both, labeled; name the snooping exposure + 3 mitigations** (D-03/D-04/D-05) |
| demo.jl repro scope | full reported-scale re-run / fast frozen-artifact chain / two-tier | **Two-tier: fast default + `--full`** (D-06) |
| Full-build-out call | roadmap DAG as-is / prioritized axis / data-scale gate first | **DAG, gated on a re-pre-registered calibration close-out; P11∥P13 lead, P12 follows** (D-07/D-08) |

## Auto-Resolved (user away)

All four areas auto-resolved with the recommended option, grounded in:
- `05-04-SUMMARY.md` (original pre-registered all-FAIL run)
- `STATE.md` iter1/iter2 accumulated-context notes (OOD→PASS, SBC ECE-green, BF corr 0.936 + clamped-baseline diagnosis)
- `consts.jl` (byte-unchanged pre-registration)
- `ROADMAP.md` Phase 6 SC1–3 + Phase 8–16 DAG

Three decisions flagged **[USER-OWNED]** in CONTEXT.md (verdict, honesty framing, build-out
call) as the ones most worth human confirmation before planning.

## Rationale highlights

- **Conditional Go, not clean Go:** the literal pre-registered gates did not pass as written;
  honesty forbids claiming they did. But OOD passes, SBC is calibrated (ECE green), BF fails
  only via a clamped-baseline tail artifact, and the >100× thesis was proven in Phase 4 — so
  the evidence supports proceeding, gated on a bounded close-out.
- **Both number sets reported:** the retrains happened after the pre-registered result was
  seen (a real data-snooping exposure); the memo must name that and lean on consts-unchanged +
  disjoint-DEV-seed + single-confirmatory-VAL-run rather than hide it.
- **Two-tier demo:** SC1 needs CPU-only reproducibility; a fast default keeps demo.jl in the
  method's own milliseconds-to-minutes spirit, `--full` gives the faithful path.

## Deferred ideas captured

- Calibration close-out (200k retrain, non-clamped BF baseline, fresh pre-registration) → Phase 7.
- Re-enabling OOD PP channel in the reported OR-fusion → Phase 7 hardening.
- RxInfer independent ADVI cross-check → post-Go, paper nice-to-have (BACK-01).
