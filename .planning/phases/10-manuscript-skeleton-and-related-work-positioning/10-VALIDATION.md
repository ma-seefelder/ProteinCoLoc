---
phase: 10
slug: manuscript-skeleton-and-related-work-positioning
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-02
---

# Phase 10 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> This is a documentation phase: the **Typst compile gate is the test framework** — the compiler *is* the test.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | `typst compile` exit-code gate (no unit-test framework; Typst 0.15.0 compiler is the test) |
| **Config file** | `manuscript/main.typ` (entry) + `manuscript/build.sh` (compile + content assertions) |
| **Quick run command** | `typst compile manuscript/main.typ manuscript/main.pdf` |
| **Full suite command** | `bash manuscript/build.sh` |
| **Estimated runtime** | ~2 seconds (incremental compile sub-second) |

**Compile-gate gotchas (measured against installed 0.15.0):**
- An undefined `@citekey` or `<label>` makes `typst compile` **fail exit 1** — seeded stubs must NOT reference not-yet-existing keys/labels. (Free regression check.)
- Piping `typst compile | head` masks the exit code (reports the pipe tail's 0). The build script MUST read Typst's own exit status directly.

---

## Sampling Rate

- **After every task commit:** Run `typst compile manuscript/main.typ manuscript/main.pdf` (must stay exit 0 — cheap, sub-second incremental).
- **After every plan wave:** Run `bash manuscript/build.sh` (compile + content assertions).
- **Before `/gsd:verify-work`:** `typst compile` green with the claim table AND related-work matrix rendering in the PDF (D-12).
- **Max feedback latency:** ~2 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 10-01-xx | 01 | 1 | SC1 (compile) | — | N/A (docs phase) | build/smoke | `typst compile manuscript/main.typ manuscript/main.pdf; echo $?` (==0) | ❌ W0 | ⬜ pending |
| 10-02-xx | 02 | 2 | D-06 claim table | — | N/A | build + source | `test -f manuscript/claims.typ` + compile | ❌ W0 | ⬜ pending |
| 10-03-xx | 03 | 2 | SC2 Tapqir 4-axis | — | N/A | content check | `grep -qi "amortiz\|SBC\|registration\|spatial" manuscript/sections/related_work.typ` | ❌ W0 | ⬜ pending |
| 10-04-xx | 04 | 2 | SC3 figure specs | — | N/A | content check | `grep -c "F[0-9]" manuscript/figures/specs.typ` ≥ 7 | ❌ W0 | ⬜ pending |
| 10-05-xx | 05 | 3 | D-12 gate script | — | N/A | build | `bash manuscript/build.sh; echo $?` (==0) | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky. Exact task IDs finalized by the planner.*

---

## Wave 0 Requirements

- [ ] `manuscript/main.typ` — entry with set-rules + `#include`s + `#bibliography("refs.bib")` (SC1)
- [ ] `manuscript/claims.typ` — claim array-of-dicts + `claim_table()` render function (D-06)
- [ ] `manuscript/refs.bib` — seeded citations, only keys actually referenced (D-11)
- [ ] `manuscript/sections/*.typ` — stub files with `// TODO(phase N)` markers (D-05)
- [ ] `manuscript/sections/related_work.typ` — comparison matrix + Tapqir 4-axis prose (D-08/D-09)
- [ ] `manuscript/figures/specs.typ` — figure spec entries (D-10)
- [ ] `manuscript/build.sh` — exit-code-correct compile gate (D-12)
- [ ] Framework install: none — Typst 0.15.0 already on PATH

*No test framework install needed — the Typst compiler is present and is the validation surface.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Claim table + related-work matrix visually render in the PDF | D-12 | Compile exit 0 proves it builds; visual confirmation that tables render (not empty) needs a human/PDF glance | Open `manuscript/main.pdf`, confirm the claim table and related-work matrix appear with rows populated |
| Tapqir per-axis cells are factually defensible | D-09 | Three matrix cells (formal SBC, explicit BF, registration-UQ) flagged ASSUMED in research — accuracy is a human judgment against Ordabayev et al. 2022 | Verify each ASSUMED cell against the eLife paper before locking prose |

---

## Validation Sign-Off

- [ ] All tasks have an automated compile/content verify or a Wave 0 dependency
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags (`typst watch` not used in gate)
- [ ] Feedback latency < 5s
- [ ] `nyquist_compliant: true` set in frontmatter once plans satisfy the above

**Approval:** pending
