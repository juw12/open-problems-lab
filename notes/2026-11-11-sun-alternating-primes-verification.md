---
Target: prizeproblems.org — problem 043 (Sun's personal prize), forum comment
Date: 2026-11-11
Review gate: PASS —
  - m ≤ 10⁸: panel PASS (review_panel/lab_prize_20261110_5.md) + referee PASS (review_notes/lab_prize_20261110_5.md);
  - m ≤ 3×10⁸: panel PASS (review_panel/lab_prize_261111_5.md) + referee PASS (review_notes/lab_prize_261111_5.md).
  - m ≤ 5×10⁸: panel PASS (review_panel/lab_prize_261111_9.md) + referee PASS (review_notes/lab_prize_261111_9.md).
  - m ≤ 8×10⁸: panel PASS (review_panel/lab_panel_261112_15_prize4.md) + referee PASS (review_notes/lab_prize_261112_4.md).
  - m ≤ 9×10⁸: panel PASS (review_panel/lab_panel_261112_18_prize9.md) + referee PASS (review_notes/lab_prize_261112_9.md).
  - m ≤ 10⁹: panel PASS (review_panel/lab_panel_261112_23_lab_prize_261112_11.md) + referee PASS (review_notes/lab_prize_261112_11.md).
  Only panel+referee-passed results are eligible for a public draft; this one
  qualifies as a computed extension.
Scope: COMPUTED EXTENSION ONLY (empirical verification m ≤ 10⁹, including
  m ≤ 10⁸, m ≤ 3×10⁸, m ≤ 5×10⁸, m ≤ 8×10⁸ and m ≤ 9×10⁸). No proof, no solution, no prize claim — the $1,000 is for the
  first proof, and the conjecture remains open.
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# Extending the verification of Sun's alternating-primes conjecture to m ≤ 10⁹

> **Verification only — no proof, no claim on the prize.** I verified the
> conjecture on a hundred times the published range; I did not prove it.

## What I checked

Zhi-Wei Sun's conjecture (Prime-AP.pdf, slide 11, April 2012; $1,000 for the
first proof): for every positive integer m there are consecutive primes
pₖ,…,pₙ whose alternating sum

    m = pₙ − pₙ₋₁ + ⋯ + (−1)ⁿ⁻ᵏ pₖ

equals m, with the primes not exceeding 2m + 2.2√m (m even) or m + 4.6√m
(m odd). Sun verified it for m ≤ 10⁷. I verified it for **m = 1 … 10⁹ —
0 violations out of 1,000,000,000** (a hundredfold extension of the
published frontier), using exact integer arithmetic throughout.

- The four worked examples from the slides reproduce exactly: 10 = 17−13+11−7+5−3
  (top 17), 20 = 41−37+31−29+23−19+17−13+11−7+5−3 (top 41), 303 → top 383,
  2382 → top 4871.
- The bound is tight (exact min top prime equals the bound) at exactly two
  values across all of 1…10⁹: **m = 303** (383 = ⌊303 + 4.6√303⌋) and
  **m = 2382** (4871 = ⌊2·2382 + 2.2√2382⌋) — the two published examples.
  No new bound-tight or near-tight case appears at scale.
- For m > 10⁶ the ratio t(m)/B(m) (min top over the bound) never exceeds
  0.9999327 (worst at m = 257,787,004); the smallest bound gap B−t = 21,676
  at m = 101,641,678 (t = 203,283,859, B = 203,305,535). The smallest gap in
  the new 7×10⁸+1…8×10⁸ range is 119,105 at m = 700,542,151, in the new
  8×10⁸+1…9×10⁸ range is 127,073 at m = 800,094,929, and in the new
  9×10⁸+1…10⁹ range is 135,201 at m = 902,523,441, so both global records
  are unchanged. Largest top prime used: 1,999,781,011 at
  m = 999,999,762 (even m, bound 2m + 2.2√m = 2,000,069,570).
- Structural pattern: t(m) ≈ 2m − O(√m) for even m (gap to 2m grows like
  ≈ 4–9·√m near 10⁸); t(m) ≈ m for odd m.
- The alternating prefix sums s₁..s₁₆ = 2, 1, 4, 3, 8, 5, 12, 7, 16, 13,
  18, 19, 22, 21, 26, 27 agree with OEIS A008347.

## How it was checked (independent re-runs)

- Gate: `verify_all.py 1 200000 gate` (~1 s): 0 violations, all four slide
  examples reproduced, and 1,013 (m, min-top, bound) triples from an earlier
  independent exact run match to the value.
- Independent verifier (`verify.py`, a separate implementation: its own
  odd-only sieve to 199,993,643, direct alternating-sum recomputation, no
  prefix sums): 1,019 representations (19 extreme m values + 1,000 random)
  re-checked, 0 failures.
- Review: the panel re-ran the full 10⁸ verification (byte-identical
  result) and the 1,019-rep verifier (0 failures); the referee re-ran the
  200k gate (0 violations) and 1,000 reps, and confirmed the committed
  99,000,001–100,000,000 min-top slice is byte-identical to the stored
  array.
- m ≤ 3×10⁸ gate 3 (SHA-256 byte-identical re-verify of the three committed
  10⁸ min-top arrays) and gate 3B (the new 3×10⁸ pipeline `verify_all3.py`
  reproduces the committed 10⁸ result byte-identically: 0 violations, ties = 2
  at m = 303/2382, all four slide examples reproduced) — both PASS.
- m ≤ 3×10⁸ independent verifier (`verify3b.py`: its own odd-only sieve to
  600,048,105, direct alternating sums, no prefix sums, no data arrays): 600
  reps (4×150 per slice, seed 7) plus whole-range checks over all 2×10⁸ stored
  t-values — 600/600 reps, 0 failures.
- Review of the 3×10⁸ result: the panel re-ran both gates (byte-identical)
  and confirmed every headline number reproduces from the stored t-arrays and
  `agg3b.json`; the extension is grounded against Sun's published m ≤ 10⁷
  (Prime-AP.pdf slide 11, April 2012) with no later published verification
  found.
 - m ≤ 5×10⁸ gate 6 (SHA-256 byte-identical re-verify of the two committed 5×10⁸
   min-top slices, 4×10⁸+1…5×10⁸) and independent verifier (`verify6.py`: its own
   odd-only sieve, direct alternating-sum recomputation, no prefix sums, no data
   arrays): 360 reps (per-slice, seed 7) plus whole-range checks over all 1×10⁸ new
   t-values — 360/360 reps, 0 failures; aggregation in `agg6.json` reproduces every
   headline number — PASS.
 - Review of the 5×10⁸ result: the panel and the independent referee both re-ran
   the 5×10⁸ pipeline and verifier (PASS, 2026-11-12) and confirmed every headline
   number reproduces from the stored t-slices and `agg6.json`.
 - m ≤ 8×10⁸ gate 9 (`gate9.py`): SHA-256 byte-identical re-verify of the committed
   1×10⁸ min-top t-arrays (3 ranges), cross-gate byte-match, 0 violations, slide
   examples t(303)=383 and t(2382)=4871, and new-array integrity — PASS.
 - m ≤ 8×10⁸ independent verifier (`verify9.py`: its own odd-only sieve to
   1,599,890,023, direct alternating-sum recomputation, no prefix sums, no data
   arrays): 360 reps (per-slice, seed 7) plus whole-range checks over both new
   5×10⁷ slices (0 violations, 0 over-bound, 0 below m, 0 non-prime) — 360/360
   reps, 0 failures; aggregation in `agg9.json` reproduces every headline number —
   PASS.
 - Review of the 8×10⁸ result: the panel re-ran gate 9 (PASS) and the independent
   referee re-ran the 8×10⁸ verifier (360/360 reps, 0 failures), 2026-11-12; every
   headline number was confirmed against `agg9.json` / `verify9.json`.
  - m ≤ 9×10⁸ gate 10 (`gate10.py`): SHA-256 byte-identical re-verify of the committed 1×10⁸ min-top t-arrays (3 ranges), cross-gate byte-match, 0 violations, slide examples t(303)=383 and t(2382)=4871, new-array integrity — PASS.
  - m ≤ 9×10⁸ independent verifier (`verify10.py`: its own odd-only sieve to 1,799,879,071, direct alternating-sum recomputation, no prefix sums, no data arrays): 360 reps (per-slice, seed 7) over both new 5×10⁷ slices (0 violations, 0 over-bound, 0 below m, 0 non-prime) — 360/360 reps, 0 failures; aggregation in `agg10.json` reproduces every headline number — PASS.
  - Review of the 9×10⁸ result (2026-11-12): the panel re-ran gate 10 (PASS); the independent referee re-ran verify10 (360/360 reps, 0 failures) and additionally re-derived every value with a from-scratch bound checker over all 100,000,000 new m (0 violations, 0 over-bound, 0 tie-bound, max top 1,799,879,071, smallest gap 127,073) — every headline number confirmed against `agg10.json` / `verify10.json`.
  - m ≤ 10⁹ gate 11 (`gate11.py`): SHA-256 byte-identical re-verify of the committed 1×10⁸ min-top t-arrays (3 ranges), cross-gate byte-match, 0 violations, slide examples t(303)=383 and t(2382)=4871, new-array integrity — PASS.
  - m ≤ 10⁹ independent verifier (`verify11.py`: its own odd-only sieve to 1,999,781,011, direct alternating-sum recomputation, no prefix sums, no data arrays): 360 reps (per-slice, seed 7) over both new 5×10⁷ slices (0 violations, 0 over-bound, 0 below m, 0 non-prime) — 360/360 reps, 0 failures; aggregation in `agg11.json` reproduces every headline number — PASS.
  - Review of the 10⁹ result (2026-11-12): the panel re-ran gate 11 (PASS) and a fresh independent base-range re-derivation (byte-identical: n_viol=0, t(303)=383, t(2382)=4871); the independent referee re-ran gate 11, the aggregator over both new 5×10⁷ slices (100,000,000 new m: 0 violations, 0 over-bound, 0 tie-bound) and verify11 (360/360 reps, 0 failures), and confirmed `agg11.json` leaf-identical — every headline number confirmed.

## Honest limits

- This is **empirical verification, not a proof**: no counterexample was
  found in 10⁷ < m ≤ 10⁹, and the conjecture for all m remains open. The
  $1,000 prize is explicitly for the first proof.
- This post cites the full panel- and referee-cleared frontier m ≤ 10⁹.
  10⁹ is the last int32 milestone: extending further needs 64-bit indices.
  The 100× is cumulative over the sequence of verification jobs (not one run),
  and the baseline is Sun's slide-deck verification of m ≤ 10⁷.
- The bound is tight at the two small published examples (m = 303, m = 2382),
  so any proof needs fine-grained control, not just density.
- The 5×10⁸ range was run as one continuous pipeline; at the boundary value
  m = 5×10⁸ the bound is B = 1,000,049,193 and the min top is t = 999,834,067
  (gap 215,126).

## Provenance

- Method: a per-m exact min-top "window method" built on the parity
  structure of the alternating prefix sums (sᵢ = pᵢ − sᵢ₋₁): every
  alternating sum over a block is sᵢ − sⱼ (same parity) or sᵢ + sⱼ (opposite
  parity), giving O(N) passes per block length; sieve to 2·10⁸, prefix sums,
  membership tables. Exact integer arithmetic everywhere (the bounds via
  integer square roots: ⌊2m + 2.2√m⌋ = 2m + isqrt(121m)//5 and
  ⌊m + 4.6√m⌋ = m + isqrt(529m)//5).
- Scripts: `prize_work/ppl-043/verify_all.py` (the 10⁸ result),
  `prize_work/ppl-043/verify_all3.py` (the 3×10⁸ result),
  `prize_work/ppl-043/verify_fb5.py` + `gate6.py` (the 5×10⁸ pipeline + gate),
  `prize_work/ppl-043/verify_fb8.py` + `gate9.py` (the 8×10⁸ pipeline + gate), `prize_work/ppl-043/verify_fb10.py` + `gate10.py` (the 9×10⁸ pipeline + gate), `prize_work/ppl-043/verify_fb11.py` + `gate11.py` (the 10⁹ pipeline + gate);
  `verify.py`, `verify3b.py`, `verify6.py`, `verify9.py`, `verify10.py` and `verify11.py` (independent
  verifiers); statistics in `prize_work/ppl-043/agg3b.json`, `agg6.json`, `agg9.json`,
  `agg10.json` and `agg11.json`; full record in `prize_work/ppl-043.md`.
- Review cross-references — m ≤ 10⁸: panel `review_panel/lab_prize_20261110_5.md`,
  referee `review_notes/lab_prize_20261110_5.md`; m ≤ 3×10⁸: panel
  `review_panel/lab_prize_261111_5.md`, referee `review_notes/lab_prize_261111_5.md`;
  m ≤ 5×10⁸: panel `review_panel/lab_prize_261111_9.md`, referee
  `review_notes/lab_prize_261111_9.md`; m ≤ 8×10⁸: panel
  `review_panel/lab_panel_261112_15_prize4.md`, referee
  `review_notes/lab_prize_261112_4.md` (2026-11-12).
- Review cross-reference, m ≤ 9×10⁸: panel `review_panel/lab_panel_261112_18_prize9.md`, referee `review_notes/lab_prize_261112_9.md` (2026-11-12).
- Review cross-reference, m ≤ 10⁹: panel `review_panel/lab_panel_261112_23_lab_prize_261112_11.md`, referee `review_notes/lab_prize_261112_11.md` (2026-11-12).
- AI-assisted: yes — the scripts were written with AI assistance, and I
  checked their output against independent re-runs (the panel's full
  re-run and the referee's own verifier) before standing behind it.

Happy to share the script if it is useful.
