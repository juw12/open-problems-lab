# Referee note — lab_prize_261112_9 (referee job lab_review_261112_26)

**Verdict: pass** — extension reconfirmed; every headline number reproduces independently; no overclaim.

**One line:** Sun's alternating-sums-of-consecutive-primes conjecture holds for m = 1..9e8 (90× the published 2012 frontier of 1e7); the new 8e8+1..9e8 range (1e8 values) has 0 violations, 0 over-bound, 0 tie-bound, max top prime 1799879071 — all independently re-derived.

**Re-run / evidence (this session, single-threaded):**
- `gate10.py` fresh → **GATE PASS**: (a) 3 committed 1e8 base t-arrays SHA-256 match; (b) cross-gate `_gate10/t_1_200000.npy` byte-matches committed `t_1_200000.npy`, 0 violations, and reproduces the worked examples t(303)=383, t(2382)=4871; (c) new-array SHAs recorded.
- **My own independent checker** (own from-scratch bound `bnd(m)=(2m+isqrt(121m))//5` even / `(m+isqrt(529m))//5` odd, no shared code) over all 1e8 new m: `n_zero=0, n_over_bound=0, n_tie_bound=0, n_below_m=0`; `max_t=1799879071`; `smallest (bnd−t) gap=127073`; ratios p50/p99/p99.9/max = 0.9998552/0.999896/0.9998973/0.9998988 → rounds to the claimed 0.99986/0.99990/0.99990/0.99990.
- Job's own independent verifier `verify10.py` (own sieve to 1799879071, ~88.9M primes) run fresh on both slices: slice 800000001–850000000 zero/over/below/not_prime = 0, 180 reps; slice 850000001–900000000 zero/over/below/not_prime = 0, 180 reps; **360 direct alternating-sum recomputations, 0 failed → VERIFY10 PASS.**
- The 303/2382 examples and the all-zero violation counts are the decisive independent checks; the argmax t-value (1799879071) is confirmed by my own argmax of the committed t-arrays.

**Class check (extension):** the published frontier is fixed at m ≤ 1e7 (Sun, *Primes in Arithmetic Progressions*, Apr 2012, slide 11, fetched via web.archive.org 2026-11-10). The Sun conjecture is unproven, so verification to 9e8 does **not** follow from any published theorem; it is a genuine range extension (90× the published frontier) and consistent with the lab's standing convention (job lab_prize_261112_4 at 8e8 was already referee-passed as "extension"). supported_class = **extension**.

**Wording:** the dossier (`ppl-043.md`) is honest — "Empirical only — no proof claim; a counterexample anywhere in 1..9e8 would have been recorded automatically (none found)." No "proof"/"solved"/"proven". No FINDINGS.md fix needed.

**Note (hidden-number hygiene):** the argmax m and the smallest-gap m show as 64-hex in my raw re-run output (8-digit runs); the job states them as 899,999,880 and 800,094,929. I independently confirmed the t-VALUES (max 1799879071, gap 127073) and the zero-violation counts across all 1e8 m in-script, so the headline does not depend on reading those m-pointers directly.
