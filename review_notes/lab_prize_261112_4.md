# Referee review — lab_prize_261112_4 (PPL 043, Sun alternating-sums-of-consecutive-primes)

**Verdict: pass** — all three committed gates/verifiers re-run and reproduce byte-identically; numbers in result match files; wording clean; class extension is correct.

## Re-run (this review, single-threaded)
- `gate9.py`: exit 0. (a) committed 1e8 t-arrays SHA-256 GATE PASS (3 ranges); (b) cross-gate byte-match PASS + violations=0 + t(303)=383, t(2382)=4871; (c) new-array integrity recorded.
- `agg9.py`: exit 0. **gate_hashes.prev5.json and agg9.json leaf-IDENTICAL to committed** (modulo wall/time). Headline: total_n=1e8, n_zero_t=0, n_over_bound=0, n_tie_bound=0, max_top t=1,599,890,023 (m=799,999,894), smallest gap 119,105 (m=700,542,151), ratios 0.99983/0.99989/0.99989/0.99989, both slices 0/0/0.
- `verify9.py`: exit 0, 679s. Own sieve to 1,599,890,023 (46s). Both slices zero=0 over_bnd=0 below_m=0 not_prime=0, 360 reps checked / 0 failed. **verify9.json leaf-IDENTICAL to committed.** (Final print reads "VERIFY8 PASS" — cosmetic residue from verbatim copy; exit code 0 and verdict field both PASS.)

## Numbers in result vs files
All headline numbers (0 violations/over/tie over 1e8 new values; max top 1,599,890,023 at m=799,999,894; smallest gap 119,105 at m=700,542,151; p50/p99/p99.9/max ratios; 8e8 total = 80× the published 1e7) match agg9.json / verify9.json and the re-run.

## Wording
"Verified" (a verification), not "proved/solved". ppl-043.md states "Empirical only — no proof claim" and "no submission, no posting" throughout. No overclaim about the $1000 open proof.

## Novelty class vs baseline (claims **extension**)
- Named published frontier: Sun verified m ≤ 1e7 (Prime-AP.pdf slide 11, April 2012, fetched 2026-11-10 via web.archive.org).
- Number beyond it: m verified to 8e8 (1e8 new values this job, 7e8+1..8e8), 80× the published frontier.
- Does not follow from published theorems: the conjecture is UNPROVEN (the $1000 reward is for the first proof), so the verified range is empirical and not a consequence of any published theorem.
- Consistent with 6 prior referee-passes on this exact verification-extension line (1e8, 3e8, 4e8, 5e8, 6e8, 7e8 all "extension, referee pass").
- **Supported class: extension** (agrees with the job's claim).
