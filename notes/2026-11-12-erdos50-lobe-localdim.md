---
Target: Erdős #50 — the |S|=7 c~0.65 interior lobe: exact prime fingerprint + local-complexity gap
Date: 2026-11-12
Review gate: panel PASS `review_panel/lab_panel_261112_13_erdos25.md` + referee PASS `review_notes/lab_erdos_261112_25.md` (2026-11-12)
Scope: extension — first measurement of the densest interior lobe of Schoenberg's f(c) (never measured before: no local exponent, no prime content); NOT a proof, NOT a solution
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# Erdős #50: the densest interior lobe of f(c) — exact prime fingerprint and local-complexity gap

> This is a numerical extension. The published interior results (Erdős 1995: f is purely
> singular, f′ = 0 a.e.; Weingartner 2007: the tail as c→0; Tenenbaum–Toulmonde 2006: c→1)
> do not measure the interior. This run measures the densest interior region — the
> |S| = 7, c ≈ 0.65 lobe of Schoenberg's f(c) = dens{n: φ(n) < c n} — for the first
> time: its exact prime-content fingerprint and its local-complexity contrast with the
> |S| ≤ 4 core.

## What I checked

The |S| = 7, c ≈ 0.65 lobe holds **162,608 points at N = 10⁷** — the densest interior
region of f(c), never measured before for any local exponent or prime content.

**Exact prime fingerprint** (fraction of the 162,608 lobe points whose index set S
contains the prime p; exact integer counts, reproduced exactly by the review re-runs):

- most common: **3 → 56.9%** (92,465 pts); **5 → 39.3%** (63,862); **11 → 19.1%**
  (31,079); **23 → 16.6%** (27,004)
- least common: **7 → 11.9%** (19,336); **17 → 13.0%** (21,120); then ~13.2% each for
  191 / 193 / 197 / 199

That the *small* primes 7 and 17 are rarer in the lobe than 191–199 is the open probe
this fingerprint raises: the lobe is not simply "dense in small primes".

**Local-complexity contrast.** A finite-N local-complexity diagnostic (log–log fit of
the one-sided-below profile around each of the top 100 lobe points) gives a snapshot
value h ≈ 1.485 (sd 0.007) for the lobe, versus h ≈ 0.0804 (clean, N-stable) for the
|S| ≤ 4 core — a **lobe-minus-core gap of 1.4042**. The densest lobe is locally far
more complex than the core, and is still resolving as N grows (below).

The densest single point is **S = [5, 11, 19, 79, 109, 163, 193]**, c ≈ 0.666
(windowed-slope growth 70.902) — a useful anchor for any follow-up.

## Honest framing (read before the numbers)

- The h ≈ 1.485 figure is a **finite-N local-complexity snapshot, not a certified
  local Hölder exponent or limit**: 0/100 of the top points fit a clean power law
  (r² > 0.9; mean r² = 0.68), and the value is **not N-stable** (5×10⁶ → 10⁷ shifts
  it by 0.26; 0/100 points stable within 0.01; window-shift 0.09; robust_to_N_and_W =
  False). Treat it as "the lobe is locally more complex and still-resolving", not as
  a measured exponent.
- The **prime fingerprint is exact at N = 10⁷** (integer counts), but its N-stability
  is not certified.
- The load-bearing new facts are the **exact fingerprint** and the **lobe-vs-core gap
  1.4042**, not the absolute 1.485.

## Verification

- The panel re-ran `prize_work/ppl-003/lobe_localdim.py` on a scratch copy (committed
  artifacts untouched): **G1** n_divergent = **1,289,274** (exact), per-|S| divergent
  counts match; **K1** lobe |S| = 7 points = **162,608**, kill test not fired;
  **G2** 60-core cross-check reproduced exactly (mean h = 0.0804, sd 0.0147).
- Every headline number above reproduced exactly from the re-run; data-leaf diff of the
  re-run JSON vs the committed `lobe_localdim_N10000000.json`: **0** (bit-identical).
- The independent referee re-verified the same gates and the fingerprint counts
  (2026-11-12).

## Honest limits

- This is a **measurement**, not a proof. It characterizes one interior region of f(c);
  it does not prove non-differentiability anywhere, nor compute f(c) at any c.
- The lobe value h ≈ 1.485 is a still-resolving finite-N snapshot (see above), not a
  certified exponent.
- Erdős #50 remains open.

## Provenance

- **Scripts:** `prize_work/ppl-003/lobe_localdim.py`, output
  `prize_work/ppl-003/lobe_localdim_N10000000.json`; full record in
  `prize_work/ppl-003.md`.
- **Review:** panel PASS + referee PASS (2026-11-12) —
  `review_panel/lab_panel_261112_13_erdos25.md`, `review_notes/lab_erdos_261112_25.md`.
- **AI-assisted:** performed by an AI agent under human supervision. All numbers are
  reproducible from the committed scripts.

Not posted to the problem site. The scripts are in this repository.
