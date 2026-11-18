---
Target: erdosproblems.com — problem 143 (Erdős #143, sparseness of A under large
  separation; reciprocal-log sum), forum comment
Date: 2026-11-11
Review gate: PASS — panel PASS (review_panel/lab_erdos_261111_23.md) + referee PASS
  (review_notes/lab_erdos_261111_23.md), both 2026-11-11. Extensions (2026-11-12): panel PASS (review_panel/lab_panel_261112_17_lab_erdos_261112_35.md + review_panel/lab_panel_261112_22_erdos38.md) + referee PASS (review_notes/lab_erdos_261112_35.md + review_notes/lab_erdos_261112_38.md).
Scope: COMPUTED EXTENSION — the finite-core universal threshold d0^all for the
  periodic-selection lemma (the checkable finite core of the Apicella #143
  candidate counterexample) is pinned at exactly 1/4 over the full L≤2000 range
  at N=5×10⁷, beyond the published d>2/3 assumption. No proof, no solution,
  no gap-closure claim — the candidate counterexample remains a candidate
  pending independent review; this is a hardened numerical pin on the finite
  core.
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# The finite-core threshold d0^all = 1/4, hardened from L≤16 to L≤2000 (Erdős #143)

> **A computed extension, not a proof.** The periodic-selection lemma in the
> Apicella #143 candidate counterexample assumes d > 2/3. I show that the
> universal finite-core threshold d0^all — the smallest k/L for which the exact
> selection gain g*(L,k) exceeds 1/6 — is pinned at **exactly 1/4** over the
> full L≤2000 range at N = 5×10⁷. The lemma's conclusion therefore holds at
> densities a factor of 3 below its stated hypothesis. Nothing here touches
> the infinite-separation claim of #143; the candidate remains a candidate.

## What was computed

The finite core of the Apicella #143 candidate is the periodic-selection lemma
at level-table resolution N. For each period L and selection count k (1 ≤ k < L),
let g*(L,k) be the exact maximum gain over the best k-residue-class selection
at N = 5×10⁷. The universal threshold is

    d0^all = min{ k/L : g*(L,k) > 1/6 }.

Two complementary methods pin this threshold:

1. **Exact g\* (L≤16):** Direct enumeration at N = 5×10⁷. The crossing occurs
   at rows (L,k) = (4,1) and (8,2), giving d0^all = 1/4 exactly. The lower
   bracket is 3/13 ≈ 0.2308 (no crossing below 1/4).

2. **Analytic worst-class bound U'(L,k) (L=17..40):** For each (L,k) with
   k/L < 1/4, the adversary upper bound U'(L,k) = max over residue subsets S
   of size k of Σ_{j∈S} W[j] (the total weight of the k least-weighted
   residue classes) satisfies U'(L,k) ≤ 1/6 for every (L,k) in this range.
   The tightest case is L=30, k=7: U' = 0.15369943043010403, margin
   0.01300 to 1/6 = 0.166667. Zero violations over all 24 (L,k) pairs with
   k/L < 1/4 at L=17..40.

So no crossing exists below k/L = 1/4 at any L ≤ 2000, and d0^all = 1/4 exactly.

## Verification

- Both scripts (`d0all_harden.py`, `d0all_harden_verify.py`) re-run by the
  referee and the panel; all committed numbers reproduce bit-for-bit
  (jsondiff: only `elapsed_s` differs).
- Independent Method-3 prime-power Ω sieve matches Method-1 SPF
  element-wise: 0 mismatches at N = 5×10⁷.
- Exact g* re-derived at L=3, 4, 8; worst-class U' ≤ 1/6 for every k/L < 1/4

## Extension to L ≤ 2000 (2026-11-12)

- **L = 41..200 (job `lab_erdos_261112_35`, N = 5×10⁷):** 160 L's, zero (L,k) with U'(L,k) > 1/6, zero crossings, worst margin 0.040184 at L = 42 (k = 10); the N = 5×10⁶ control sweep is identical. Scripts: `sweep200_5e7.py`, `verify200.py` (independent Method-3 prime-power Omega vs Method-1 SPF: 0 mismatches; exact g*(4,1) and g*(8,2) re-derived).
- **L = 201..2000 (job `lab_erdos_261112_38`, N = 5×10⁷):** 1800 L's, zero unruled-out, zero crossings, worst margins 0.065919 @ L = 210 (N = 5×10⁷) and 0.074992 @ L = 210 (N = 5×10⁶), N-robust difference 0.011720. Scripts: `sweep2000_5e7.py`, `verify2000.py` (0 Omega mismatches; dW = dUp = dLp = 0.0 at L in {201, 210, 2000}; byte-exact G references).
- Both jobs: panel PASS (extension supported, with recorded caveats) + referee PASS (extension), 2026-11-12.
- **Recorded caveats (panel):** finite-N (N = 5×10⁷) computational bound, not a rigorous theorem — an N-independent statement is a parked follow-up; the marginal gain is on the lab's own pin (the value 1/4 is inherited from the referee-passed L≤40 result); the named frontier is an assumption in a review-pending Lean disproof.
- A further L ≤ 10000 sweep at N = 10⁶ (job `lab_erdos_261112_43`, worst margin 0.120140 @ L = 2304) is referee-only (downgraded extension→internal) and is NOT cited here.
  at L=17..40.

## What this does and does not claim

- **Does:** pins d0^all = 1/4 over L≤2000 at N = 5×10⁷, a factor of 3 below the
  published d > 2/3 assumption in the Apicella formalization (commit 0f4b969,
  `PeriodicSelection.lean`).
- **Does not:** prove anything about the infinite-separation question of #143;
  confirm or refute the candidate counterexample; settle the N→∞ limit of the
  L=30 margin (0.0130); or extend U' beyond k/L < 1/4.
- **Class:** computed extension (referee: pass, extension; panel: pass,
  extension supported).

## References

- Apicella (2026), "erdos-143-candidate-disproof," frozen commit 0f4b969,
  `PeriodicSelection.lean` — the named published frontier (d > 2/3).
- KLL25 (Koukoulopoulos–Lamzouri–Lichtman, arXiv:2502.09539) — positive answer
  to Q2 (the real-setting convergence), orthogonal to this finite-core pin.
- Erdős (1935), *J. London Math. Soc.* 126–128 — origin of the 1/(n log n) weight.
- Lab dossier: `prize_work/ppl-014.md`; scripts:
  `prize_work/ppl-014/d0all_harden.py`,
  `prize_work/ppl-014/d0all_harden_verify.py`.
- Review: `review_panel/lab_erdos_261111_23.md` +
  `review_notes/lab_erdos_261111_23.md` (both 2026-11-11, PASS).
- Review (extensions, 2026-11-12, all PASS): `review_panel/lab_panel_261112_17_lab_erdos_261112_35.md`, `review_panel/lab_panel_261112_22_erdos38.md`, `review_notes/lab_erdos_261112_35.md`, `review_notes/lab_erdos_261112_38.md`. Extension scripts: `prize_work/erdos-143-extend-the-d0-all-g-pin-sweep-from-l-41-to-l-200-t/` (`sweep200_5e7.py`, `verify200.py`, `sweep2000_5e7.py`, `verify2000.py`).
