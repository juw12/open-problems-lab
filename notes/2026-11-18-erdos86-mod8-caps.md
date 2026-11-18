---
Target: erdosproblems.com problem 86 (C4-free subgraphs of the n-cube) forum comment; reusable as a Kaggle discussion write-up
Status: DRAFT only — not posted, not submitted; review gate PASSED (Referee PASS + Review Panel PASS, both sides)
Scope: REPRODUCTION ONLY — a corollary of a published lemma (Minamoto, arXiv:2603.29127v5, Lemma 12), applied at dimensions that paper does not evaluate; no proof, no frontier move
AI-assisted: yes (OpenCode)
Gate: `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d/`
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# The mod-8 quantised two-point cap on Q_n: it bites at exactly the even n with S_pl(n) ≡ 2 (mod 4)

*(Shorter forum version at the end. Nothing below is a proof or a solution to Problem 86.)*

For the odd-square / fully-frustrated subfamily of Problem 86, the two-point inequality gives an upper
bound on `g(n)` of the form `cap(n) = (E + S_pl(n))/2` with `E = n·2^(n−1)` and `S_pl(n)` the largest
admissible even value of `Σ h`. Minamoto's **Lemma 12** (arXiv:2603.29127v5) adds the congruence
`(h−a)(h−a−2) ≡ 0 (mod 8)`, which quantises that bound. The paper applies the lemma only at n = 6, 7, 8.

Two things came out of applying it at every n, with exact integer arithmetic:

1. **At the five open dimensions where the quantisation bites, the cap drops by exactly 1:**

   | n | plain cap (DPTV79) | quantised cap (Lemma 12) | S_Q |
   |---|---|---|---|
   | 12 | 15701 | **15700** | 6824 |
   | 18 | 1454899 | **1454898** | 550500 |
   | 20 | 6396313 | **6396312** | 2306864 |
   | 22 | 27892121 | **27892120** | 9646896 |
   | 24 | 120795955 | **120795954** | 40265316 |

   It gives nothing at n = 14 and n = 25, and nothing at any odd n.

2. **Which n sharpen has a one-line description.** With `ρ(n) = min_a (E + a(a+2)·2^(n−1))/(2a+2)` and
   `S_pl(n) = 2⌊ρ(n)/2⌋`:

   > `cap*(n) = cap(n) − 1` **iff** `n` is even **and** `S_pl(n) ≡ 2 (mod 4)`.

   For n ≥ 4, `N(a) ≡ 0 (mod 8)` for every admissible `a`, so the congruence reduces to `W·S ≡ 0 (mod 8)`
   with `W = 2a+2`: `a` even forces `S ≡ 0 (mod 4)` (the only case that bites, costing exactly 2 when
   `S_pl ≡ 2 mod 4`), `a ≡ 1 (mod 4)` forces only the DPTV parity, `a ≡ 3 (mod 4)` is vacuous. Checked
   against brute-force minimisation at **every n = 3…600: 598/598 agree**, giving **131 sharpened
   dimensions, 125 of them open dimensions > 26**.

Why this is worth writing down even though it is a corollary: it lets a reader recompute any row of the
table — or any n beyond 600 — from one residue, instead of running a search. It also explains the two
open dimensions that do *not* sharpen (n = 25 is odd; at n = 14 the binding `a = 2` gives
`S_pl = 30036 ≡ 0 mod 4`).

## What was checked (and how)

- **Exact integers only.** Python `int` / `Fraction` throughout, no floating point in either script.
- **Three routes to the table.** A from-scratch CRT/closed-form implementation (92 gates, exit 0); an
  independent verifier with no shared code that gets `S_max` by naive downward scan and the plain cap by
  cross-multiplied fractions, enumerating the programme optimum exhaustively at n = 4…8 (32 checks, exit 0);
  and the earlier lab job's gate, re-run unmodified.
- **The gate is anchored on published values.** It reproduces all 17 published exact rows (equality at the
  16 attained dimensions; n = 9 keeps its realisability slack of 4), and reproduces Minamoto's three
  published uses of Lemma 12: `S_max(6,2)=72`, `S_max(7,1)=160`, `S_max(8,2)=340`. At n = 6, where the lemma
  *is* published, the quantised cap is 132 and `g(6) = 132` is attained.
- **The multiset witnesses are exact.** At each open dimension the extremal multiset has `2^(n−1)` entries,
  `Σh = S_Q` and `Σh² = E`; each is re-checked by both scripts, and the Δ identity is evaluated both ways
  (`Σ(h−a)(h−a−2)` over the multiset vs `N(a) − W·S`).

## Honest limits (please keep these if you post)

- **This is a reproduction, not a frontier move.** The cap refinement is the residue consequence of
  Minamoto's Lemma 12; Wang (arXiv:2611.14906) already contains this cap family in its `parity-bound` form
  with `r = a+1`. The best **published** rigorous upper bound at n = 12 therefore remains the plain
  **15701**, and the sharpened 15700 is a corollary, not a published improvement.
- **The −1 changes no bracket at any open dimension.** At n = 12 the lab's own bracket is
  **[15616, 15700] — 84 wide**; moving the top by one does not narrow it materially.
- **The multisets are not fields.** Each multiset witnesses the *relaxed programme* optimum `S_Q`. Whether a
  fully-frustrated field attains `cap*` (i.e. whether `g(n)` equals it) is a different question and is open
  at every one of these dimensions. **n = 9 is the standing counterexample**: `g(9) = 1532` (Wang,
  arXiv:2611.14906, Thm 1.4) against a cap of 1536, a slack of 4 that no two-point or higher-moment
  argument visible at the profile level sees — the lab's own n = 12 field has `S = 6652`, not 6824.
- **No `g(n)` is determined anywhere here, no proof is claimed, and nothing was submitted to anyone.**
- Ties to the conjecture itself (`f(n) ≤ (½+o(1))·n·2^(n−1)`) are nil: this is an upper-bound bookkeeping
  refinement inside the odd-square subfamily.

## Provenance

- Dossier: `prize_work/erdos-86-confirm-the-mod-8-sharpened-caps-at-the-five-open-d.md`; jobs
  `lab_erdos_261117_2` (table + 92 gates + independent verifier) and `lab_erdos_261117_5` (the criterion,
  598/598 at n = 3…600), both rated `reproduction`.
- Scripts (exact, single-threaded): `fresh_mod8_table.py`, `verify_caps_independent.py`,
  `j261117_5_criterion.py`, `criterion_check_600.py`; artifacts `cap_table.json`, `cap_table.md`,
  `j261117_5_criterion.json`, plus the three gate logs.
- Reviewed from both sides: a referee re-ran all five scripts (leaf-identical output) and scanned the
  relations independently; a simulated panel re-derived the values with its own engines, including an
  independent DP at n ≤ 10 (at n = 6: DP 72 against the plain cap's 74). Their conditions are the four
  bullets in "Honest limits" above.

---

## Shorter forum version (erdosproblems.com/86)

> A small exact addendum on the odd-square side of this problem, using Lemma 12 of Minamoto
> (arXiv:2603.29127v5) — the mod-8 quantisation of the two-point inequality — at the dimensions the paper
> does not evaluate. Where it bites, the two-point cap drops by exactly 1: n = 12: 15701 → 15700,
> n = 18: 1454899 → 1454898, n = 20: 6396313 → 6396312, n = 22: 27892121 → 27892120,
> n = 24: 120795955 → 120795954. It gives nothing at n = 14, at n = 25, or at any odd n, and the pattern is
> one line: `cap* = cap − 1` iff n is even and `S_pl(n) ≡ 2 (mod 4)`, verified against brute force at every
> n = 3…600 (598/598; 131 sharpened dimensions, 125 of them open > 26).
>
> This is a corollary, not an improvement: the best published bound at n = 12 stays the plain 15701, and
> Wang's parity-bound (arXiv:2611.14906) already contains this cap family. The −1 moves no bracket of
> practical width (the n = 12 bracket [15616, 15700] is 84 wide). The extremal multisets certify the relaxed
> programme optimum only — they are not fields, no `g(n)` is determined, and n = 9 (g(9) = 1532 against cap
> 1536) shows the realisability gap is real. Scripts and gate logs are available on request; the whole
> computation is exact-integer and runs in under a minute per leg.
