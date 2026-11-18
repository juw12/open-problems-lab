---
Target: Erdős #50 — kill test for the certified-rational-enclosure route (atom sum + published tail bound) at interior c
Date: 2026-11-13
Review gate: panel PASS `review_panel/lab_erdos_261113_43.md` + referee PASS `review_notes/lab_erdos_261113_43.md` (2026-11-13)
Scope: negative — the specific certification route (exact atom sum + published tail bound) cannot certify an interior f(c); NOT a proof, NOT a solution; Erdős #50 remains open
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# Erdős #50: the atom-sum + tail-bound certification route fails at interior c

> Problem 50 (erdosproblems.com/50, $250): is there no x for which f′(x) exists and is positive? f(c) is the natural density of {n : φ(n) < c n} (Schoenberg proved existence and continuity). Live-checked 2026-11-13: OPEN, 0 comments, 0 proof claims.

## What I checked

I attempted to certify an interior value f(c) with rational error bounds using the two pieces that are each exactly computable or already published:

- the **exact atom sum**: φ(n)/n = ∏_{q|n}(1 − 1/q) depends only on rad(n), and the one-sided-below atom sum is an exact finite-N integer identity. Verified 490/490 (point, window) pairs at N = 10⁷, reproduced at N = 2·10⁷ by an exact-integer independent verifier (3/3; atoms 1,192,434 / 307,540 / 432,414; 0 monotone-in-N violations).
- a **published tail bound** for the m-truncation remainder T(M) = (1/N)#{n ≤ N : φ(n)/n < 1/2, rad(n) > M}.

The negative: at N = 2·10⁷ the **measured** T(M) stays ≥ 0.104 for **every** M ≤ N/2 (curve 0.505, 0.498, 0.482, 0.451, 0.386, 0.281, 0.104 at M = 10⁴, 3·10⁴, 10⁵, 3·10⁵, 10⁶, 3·10⁶, 10⁷), while the observed f_N(1/2) movement 10⁷ → 2·10⁷ is only 1.55·10⁻⁵ (0.5111904 → 0.51120590; 0.5112406 at 10⁸). The remainder exceeds the movement by a factor of **6705.7 (≈ 6706×)** at the tightest rung.

**Kill, and its honest framing** (per the review panel): 6706× is measured against the single smallest 10⁷ → 2·10⁷ rung movement. Against the dataset's own 10⁷ → 10⁸ spread (max 2.12·10⁻³, median 8.5·10⁻⁵, 40/99 points rising) the ratio is ≈ **49×** — still four orders of magnitude above what a certified enclosure would need (≫ 1). The KILL is unaffected by which comparison is used.

Why the remainder stays Θ(1) at interior c:

- each atom class {n : rad(n) = m} has natural density 0;
- the families n = m·P (P large prime, φ(m)/m < c) sit at **every** radical scale — e.g. n = 2P at c = 1/2 — keeping the m-truncation remainder Θ(1) up to M = N/2;
- the published tail estimate (Weingartner 2007, log A(t) as t → ∞) is a c → 0 asymptotic: it bounds neither the m-truncation remainder at fixed interior c nor |f_N − f|.

Also negative: f_N is **not monotone in N** (40/99 grid points rise between 10⁷ and 10⁸), so the "enclosures nest when N doubles" premise fails as well.

Known results in the neighbourhood (all cited in the dossier): Erdős 1995 (f is purely singular — f′ is not finite-positive a.e.); Weingartner 2007 + part II arXiv:1011.4262 (c → 0 tail only); Tenenbaum–Toulmonde 2006 (c → 1 endpoint). Closest published interior number: Delégise 1989, 0.2474 < A(2) < 0.2480 — but for the σ(n)/n density, a different function. No published computational table of f(c) at interior c was found in our search (2026-11-13).

## Verification

- **Referee PASS:** every number reproduces exactly — KILL min ratio 6705.7, 490/490 gate re-run at N = 2·10⁷, N = 10⁸ grid 20 rungs byte-equal, cross-lineage ALL_OK.
- **Panel PASS (independent):** the panel's own SPF-sieve recomputation at N = 8·10⁶ gave T(N/2) = 0.103935 and min ratio **2631** — KILL FIRED independently of the lab's rad-sieve code.
- **Dataset gates (cross-lineage):** Σφ/N² = 0.3039635516 vs 3/π²; mean φ/n = 0.6079271026 vs 6/π²; #{n ≤ 10⁸ : φ(n) = n − 1} = π(10⁸) = 5,761,455; f_N(1/2) at 10⁷ / 2·10⁷ = 5,111,904 / 10,224,118.
- Two non-blocking precision notes carried into this draft: (i) the 6706× vs ≈49× framing above; (ii) the problem's date on erdosproblems.com (1934) differs from the citation in Weingartner 2007 (1928) — hence "Schoenberg" without a year here.

## Honest limits

- This is a **negative result about one specific certification route**, not a statement about the problem: the route (exact atom sum + published tail bound) fails at interior c by ≥ 4 orders of magnitude. It does not say the problem is hopeless.
- The 6706× / 49× ratios measure this route against the empirical N-movement of one dataset; they are not theorems.
- What a certified interior value would need is a tail bound at interior c, or a modulus-of-continuity argument (cf. Toulmonde 2006) giving a certified |f_N − f| at interior c — that is exactly what is missing, and such an approach is welcome.
- No proof is claimed. **Erdős #50 remains open.**

## Provenance

- **Scripts:** `prize_work/erdos-50-certified-fc-enclosures/` — `killtest.py` (T(M) curve, N = 2·10⁷), `gate490_N2e7.py` (card's gate + exact-integer verifier at N = 2·10⁷), `grid_N1e8.py` + `check_grid.py` (N = 10⁸ dataset, 20 rungs × 99 points, cross-lineage checks), `extrap.py` (convergence-law fit).
- **Cross-check:** panel's independent SPF-sieve recomputation at N = 8·10⁶ (KILL FIRED, min ratio 2631).
- **Review:** panel PASS + referee PASS (2026-11-13).
- **AI-assisted:** performed by an AI agent under human supervision. All numbers are reproducible from the committed scripts.

Not posted to the problem site. The scripts are in this repository.
