# Referee — lab_prize_261111_5 (ppl-043, Sun alternating sums of consecutive primes)

Verdict: **pass** — 30x verification extension properly grounded; all gates/verifiers re-run single-threaded and every claimed number re-derives byte-identically; no overclaim.

Re-run (referee, single-threaded, `prize_work/ppl-043/`; worker JSONs backed up, re-runs diffed, originals restored):
- Gate `gate3.py` → **PASS**: SHA-256 of the three 1..1e8 t-arrays (`t_1_200000.npy`, `t_1_50000000.npy`, `t_50000001_100000000.npy`) all match `gate_hashes.prev.json` → the 100,000,000 prior values are byte-identical/intact (re-gated, not re-computed).
- Independent verifier `verify3b.py` (own odd-only sieve, own `isqrt` Sun bound, direct alternating sums over the whole 2e8 new range, 4 slices) → **PASS in 629 s**: all 4 slices `zero=0 over_bnd=0 below_m=0 not_prime=0`, 150 direct-sum reps each, `reps checked: 600, failed: 0`; own sieve `p_N = 599967871` (matches stored `max_top`).
- Aggregator re-derivation `agg3b.py` (own `isqrt`) → **exact match** to committed `agg3b.json` on every headline number: min gap `B−t = 21676 @ m=101641678` (t=203283859, B=203305535); worst ratio `0.9999327 @ m=257787004`; `max_top = 599967871 @ m=299999932`; boundary `m=3e8`: B=600038105, t=599967829 (gap 70276); violations = 0. Only `time_s` differed from the worker's run (ref 98.6 s vs worker 101.2 s) — all substantive numbers identical.

All numbers in `lab_results/lab_prize_261111_5.json` match the stored files exactly.

Wording: dossier/result say "Empirical only — no proof claim", "Verification alone is not the prize", "No submission, no posting", "empirical" throughout. No "proof"/"solved". A proof is required for the prize, so the empirical framing is honest.

Novelty class: claims **extension**; supported. Named published frontier = Zhi-Wei Sun, *Prime-AP.pdf*, slide 11 (April 2012), verification m ≤ 10^7 (re-verified by panel `lab_panel_20261110_28` against the archived slide; no later published verification found). The lab reaches m ≤ 3e8 = a 30x range extension beyond that frontier (and 3x the lab's prior 1e8 win); the conjecture is open, so 3e8 does not follow from any published theorem. **extension** (range).
