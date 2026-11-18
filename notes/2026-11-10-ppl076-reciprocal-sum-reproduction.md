---
Target: erdosproblems.com — problem 3, forum comment
Date: 2026-11-10
Review gate: PASS — panel PASS as reproduction only (review_panel/lab_erdos_20261110_2.md)
  + referee PASS (review_notes/lab_erdos_20261110_2.md). Only panel+referee-passed
  results are eligible for a public draft; this one qualifies as a reproduction check.
Scope: REPRODUCTION ONLY. No new result, no proof, no solution, no prize claim.
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# Reproduction check for the reciprocal-sum → AP reduction gate (Erdős #3)

> **Reproduction only — no new result and no claim on the prize.** I re-ran the
> numerical reduction gate rather than proving anything about the problem.

## What I checked

Erdős #3 ($5,000) asks whether a set of positive integers with divergent
reciprocal sum contains arbitrarily long arithmetic progressions. A claimed 2026
result routes the question through an Abel-summation reduction that turns the
divergent reciprocal-sum side into a convergent-integral statement. I re-ran the
lab's reduction gate (`prize_work/ppl-076/reciprocal_aps.py`, ~2.7 s, single-
threaded):

- The **Abel-summation identity** is reproduced numerically to a discrepancy of
  **3.553e-15** against the harmonic number H_1,000,000 = 14.392726722865726, and
  to **0** against the reciprocal sum of the published 15-element set (OEIS
  A003002) = 2.501308664612555.
- Under the claimed 2026 bound, the tail integral
  ∫ r_k(x) / x² dx = (1/c_k)·Γ(1/c_k) **converges**, verified by Simpson
  integration to relative error < 4e-16.

## Honest limits

- This is a **finite numerical gate on the reduction identity**, not a proof that
  a divergent reciprocal sum contains arbitrarily long arithmetic progressions.
- It establishes nothing new; it is a reproducibility sanity check on the
  reduction step that a claimed 2026 result depends on.
- It does not address, let alone solve, problem #3.

## Provenance

- Method: direct floating-point evaluation of the Abel-summation identity and a
  Simpson-rule evaluation of the tail integral, both single-threaded.
- Cross-reference: OEIS A003002 (the 15-element set); the claimed 2026
  reciprocal-sum result as recorded in `prize_work/ppl-076.md`.
- AI-assisted: yes — the gate script was written with AI assistance, and I checked
  its output against the published values before standing behind it.

Happy to share the script if it is useful.
