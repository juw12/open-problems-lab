---
Target: erdosproblems.com — problem 50, forum comment
Date: 2026-11-11
Review gate: PASS — panel PASS (review_panel/lab_erdos_20261110_29.md) + referee
  PASS (review_notes/lab_erdos_20261110_29.md — the referee re-ran the script
  at N=1e6 and N=1e7 and reproduced every committed number exactly). Only
  panel+referee-passed results are eligible for a public draft.
Scope: COMPUTED EXTENSION ONLY. First high-resolution interior table of
  f(c) = density of { n ≤ N : φ(n) < c·n } for 0 < c < 1, at N = 10⁷.
  No proof, no solution, no counterexample, no prize claim. Problem 50
  remains OPEN.
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# Problem 50 — first high-resolution interior table of the φ(n) < c n density (N = 10⁷)

> Computed extension only — finite-N numerical values of a known (singularity-proven) distribution function. No theorem, no claim on the prize. Problem 50 remains open.

## What I checked

- **Gates (reproduced exactly at N = 10⁷):** mean(φ(n)/n) = 0.607927101854 = 6/π²; (1/N²)·Σφ(n) = 0.303963550927 = 3/π²; f(0) = 0, f(1) = 1 − 1/N.
- **Interior CDF table (N = 10⁷), reproduced exactly by an independent re-run:**

  | c | 0.1 | 0.26 | 0.34 | 0.42 | 0.5 | 0.58 | 0.66 | 0.74 | 0.82 | 0.9 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | f(c) | 0 | 0.017957 | 0.17235 | 0.251825 | 0.475023 | 0.555692 | 0.610535 | 0.687583 | 0.745559 | 0.785871 |

- **Singularity diagnostic** (max window-averaged slope under refinement): 13.677 → 48.293 → 82.494 → 183.671 → 366.152 → 728.74 → 1445.811 as the window shrinks from 0.00999 to 5×10⁻⁵ (a 200× shrink), i.e. growth ×105.711 at N = 10⁷ (×115.70 at N = 10⁶) — steepest at a singular point, consistent with f being purely singular.
- **Stable-window-average points: 5250** — a resolution artifact: the count grows with the number of sample windows, so it is *not* a witness of differentiability.

## Honest limits

- Finite-N (N ≤ 10⁷), single-threaded and deterministic, but **no asymptotic statement**: the tail c→0 (Weingartner, Proc. AMS 135 (2007) + arXiv:1011.4262) and the endpoint c→1 (Tenenbaum–Toulmonde, Funct. Approx. 35 (2006)) are covered by the literature; this table covers the interior 0.1 ≤ c ≤ 0.9.
- Purely singular (Erdős, Resenhas 1995) is a theorem; the slope divergence above is a **diagnostic consistent with** that theorem, not a re-proof of it.
- **No verified point of finite derivative** — problem 50 (pointwise differentiability?) remains open.
- Not a proof, not a solution, not a counterexample.

## Provenance

- Method: segmented totient sieve, `prize_work/ppl-003/phi_dist.py` (deterministic, single-threaded; N=10⁶ run 3.2 s, N=10⁷ run 25.3 s); dossier `prize_work/ppl-003.md`.
- Cross-reference: review panel `review_panel/lab_erdos_20261110_29.md`, referee `review_notes/lab_erdos_20261110_29.md` (independent re-run, full-precision match).
- AI-assisted: yes — the computation and this write-up were produced by an AI research agent, and I checked its output against an independent re-run (referee, both N=10⁶ and N=10⁷, every committed number reproduced exactly) before standing behind it.

Happy to share the script if it is useful.
