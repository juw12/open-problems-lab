---
Target: erdosproblems.com — problem 710, forum comment
Date: 2026-11-11
Review gate: PASS — panel PASS (review_panel/lab_erdos_fcc0a5910ca140653d8472fd044f0d6a5c109e2d57ebd6cee295932121852265_8.md)
  + referee PASS (review_notes/lab_erdos_fcc0a5910ca140653d8472fd044f0d6a5c109e2d57ebd6cee295932121852265_8.md).
  Only panel+referee-passed results are eligible for a public draft; this one
  qualifies as a computed extension.
Scope: COMPUTED EXTENSION ONLY — the exact values f(10000), f(15000),
  f(20000) plus the flat-ratio observation on n ≤ 20000. No proof, no
  solution, no prize claim — problem 710 remains open.
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# Exact values of f(n) for Erdős–Pomerance problem 710 up to n = 20,000

> **Computed extension only — no asymptotic proof.** I computed exact values
> of f(n) well past the range of the only published exact value and checked
> what the ratio does; I did not prove the asymptotics.

## What I checked

Erdős–Pomerance #710 (£2,000, ≈$78 at 1992 rates; open since 1980): f(n) is
the least L for which there are distinct a₁,…,aₙ in (n, n+L) with k | aₖ for
each k. The only published exact value is f(10) ≤ 15 (Erdős–Pomerance 1980),
with bounds (2/√e + o(1))·n·√(log n / log log n) ≤ f(n) ≤ 1.7398·n·√(log n).

I computed f(n) exactly for n = 1…5000 by bipartite matching (Hopcroft–Karp +
binary search), cross-checked by an independent backtracking for n ≤ 30, and
pushed to n = 20,000 with 2–3 independent methods (including the referee's
own from-scratch Hopcroft–Karp):

| n | f(n) |
|---|---|
| 10 | 15 (the published value — reproduced exactly) |
| 30 | 46 |
| 100 | 161 |
| 500 | 878 |
| 10,000 | **20,626** |
| 15,000 | **31,576** |
| 20,000 | **42,561** |

- Gates: f(1)=2, f(2)=3, f(3)=4, f(4)=6, f(10)=15; backtracking and
  Hopcroft–Karp agree for all n ≤ 30.
- The ratio r_t = f(n) / (n·√log n) is **flat at 0.676–0.680** across the
  push (not drifting toward 0.5203), and a 30-point OLS fit gives
  A = 0.41089 < √2/e = 0.52026 — against the narrative that f(n)/(n√log n)
  converges to √2/e (an unverified, unreviewed AI-assisted 2026-09-23 claim,
  not cited here as a result).
- The Erdős–Pomerance lower bound at n = 20,000 is 50,421.5; the exact value
  sits at 84.4% of it (f/EP = 0.844105).

## Honest limits

- These are **exact finite values and a finite-range ratio observation**,
  not an asymptotic proof; problem 710 remains open and the
  Erdős–Pomerance 1980 bounds are the standing theoretical result.
- The flat-ratio observation is a finite-sample statement on n ≤ 20,000: it
  is evidence against one specific convergence narrative, not a replacement
  for the lower bound.
- (A later lab push to n = 50,000 exists but has not passed the review
  panel and is not cited here.)

## Provenance

- Method: bipartite matching (Hopcroft–Karp) between positions k and
  candidate values a ≡ 0 mod k in (n, n+L), binary search on L; independent
  backtracking cross-check for n ≤ 30; independent re-implementations for
  the push.
- Scripts: `prize_work/ppl-129/` (exact files in the dossier
  `prize_work/ppl-129.md`).
- Review cross-references: panel `review_panel/lab_erdos_fcc0a5910ca140653d8472fd044f0d6a5c109e2d57ebd6cee295932121852265_8.md`
  (4-method gate re-run incl. from-scratch Hopcroft–Karp) + referee
  `review_notes/lab_erdos_fcc0a5910ca140653d8472fd044f0d6a5c109e2d57ebd6cee295932121852265_8.md`
  (3 independent methods; every committed number reproduced exactly).
- AI-assisted: yes — the solvers were written with AI assistance, and I
  checked their output against independent re-runs (the panel's 4-method
  gate and the referee's own implementation) before standing behind it.

Happy to share the script if it is useful.
