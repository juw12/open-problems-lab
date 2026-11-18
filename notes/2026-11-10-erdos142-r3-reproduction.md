---
Target: erdosproblems.com — problem 142, forum comment
Date: 2026-11-10
Review gate: PASS — panel PASS (review_panel/lab_erdos_20261110_1.md) + referee
  PASS (review_notes/lab_erdos_20261110_1.md). Only panel+referee-passed results are
  eligible for a public draft, and this is the only one.
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

# Reproduction check for the r_3(N) data (Erdős #142)

> **Reproduction only — no new result and no claim on the prize.** I independently
> re-enumerated the known small values of r_3(N) rather than re-reading the paper.

## What I checked

r_3(N) is the size of the largest subset of {1, …, N} with no 3-term arithmetic
progression. I recomputed it for N = 0 … 40 with an exact branch-and-bound search
(single-threaded, ~115 s over the whole range):

- **All 16 published record values of r_3(N) (OEIS A003002) are reproduced exactly**,
  including r_3(40) = 15.
- The published optimal 15-set at N = 40 —
  {1, 2, 4, 5, 10, 11, 13, 14, 28, 29, 31, 32, 37, 38, 40} —
  is verified 3-AP-free by direct check.
- At N = 40 the value 15 is a step-function plateau (15/40 = 0.375); it is the largest
  3-AP-free size at that N, not an asymptotic statement.

## Honest limits

- This is a **finite reproduction** for N ≤ 40. It checks the *table*, and it says
  nothing about the asymptotic formula for r_k(N) that problem #142 asks for.
- It establishes nothing new; it is a reproducibility sanity check on known data.
- Because the search is exact, each value is the true maximum for that N, not a
  heuristic lower bound.

## Provenance

- Method: exact branch-and-bound over 3-AP-free subsets, one N at a time.
- Cross-reference: OEIS A003002 (r_3(N)).
- AI-assisted: yes — the enumeration script was written with AI assistance, and I
  checked its output against the published record values and the explicit N = 40 set
  before standing behind it.

Happy to share the script if it is useful.
