# Open Problems Lab

**Small, checkable computations on open maths problems. No proofs, no prize claims.**

Reproducible computations on open problems in number theory and combinatorics
(mostly from [erdosproblems.com](https://www.erdosproblems.com)), produced by an
AI-agent research lab that I run on my own PC.

> Nothing in this repository is a proof, a solution or a prize claim. Every note is one
> of: a reproduction of published values, a finite computation that goes past the
> published range, or a negative result. The work was computed and re-run by AI agents.
> **No human mathematician has reviewed it.** I am a mechanical engineer, not a
> mathematician; I built the lab and the checking process, and I am publishing the
> results so that people who know these problems can check them, use them or tell me
> what is wrong.

## Highlights

- Sun's alternating-primes conjecture checked for every m up to 10⁹ with 0 violations; the published check stops at 10⁷.
- Exact values of f(n) for Erdős problem 710 for every n up to 5,000 and at n = 10,000, 15,000 and 20,000; the lab found only n = 10 published.
- A negative result on Erdős problem 50: one natural way to certify the density f(c) cannot work, with the measurement that shows it.
- A one-line rule for where a published congruence lowers the cap in Erdős problem 86, checked for every n from 3 to 600.
- Every result must first reproduce published values, then pass a second independent method and a separate AI referee re-run.
- Drafts that turned out wrong or already known are listed openly instead of being deleted.

## What is here

| Note | Problem | Kind | What it says |
|---|---|---|---|
| [erdos142-r3-reproduction](notes/2026-11-10-erdos142-r3-reproduction.md) | Erdős #142 | Reproduction | Re-enumerates r_3(N) for N ≤ 40; all 16 record values in OEIS A003002 match. |
| [ppl076-reciprocal-sum-reproduction](notes/2026-11-10-ppl076-reciprocal-sum-reproduction.md) | Erdős #3 | Reproduction | Numerical check of the Abel-summation reduction step only (not of any proof). |
| [s5-off-diagonal-ramsey-frontier](notes/2026-11-11-s5-off-diagonal-ramsey-frontier.md) | Erdős #986 | Dated reference | Re-states and dates the published lower and upper log-exponents for r(s,t), s ≥ 5. No new bound. |
| [sun-alternating-primes-verification](notes/2026-11-11-sun-alternating-primes-verification.md) | Sun's alternating-primes conjecture | Verification | 0 violations for m ≤ 10⁹ (published range: m ≤ 10⁷). The bound is tight only at m = 303 and m = 2382. |
| [erdos710-f-n-exact-push](notes/2026-11-11-erdos710-f-n-exact-push.md) | Erdős #710 | Computed extension | Exact f(n) by bipartite matching: f(10000) = 20626, f(15000) = 31576, f(20000) = 42561. |
| [erdos143-d0all-1-4-hardened](notes/2026-11-11-erdos143-d0all-1-4-hardened.md) | Erdős #143 | Computed extension | A finite-N threshold in a candidate counterexample's selection lemma stays at 1/4 for all periods L ≤ 2000 (N = 5×10⁷). |
| [erdos50-phi-distribution-interior-cdf](notes/2026-11-11-erdos50-phi-distribution-interior-cdf.md) | Erdős #50 | Computed extension | Table of the density of {n ≤ N : φ(n) < c·n} at interior c, N = 10⁷. |
| [erdos50-band-edges](notes/2026-11-12-erdos50-band-edges.md) | Erdős #50 | Computed extension | Rational end points 8/35 and 42367104/44709361 of the band where the finite-N windowed slope diverges. |
| [erdos50-sharpest-interior-region](notes/2026-11-12-erdos50-sharpest-interior-region.md) | Erdős #50 | Computed extension | The sharpest finite-N slope growth sits in a cluster near c ≈ 0.6664. |
| [erdos50-lobe-localdim](notes/2026-11-12-erdos50-lobe-localdim.md) | Erdős #50 | Computed extension | Prime content and local-complexity measurements of the densest interior region near c ≈ 0.65. |
| [erdos50-fc-enclosure-kill-test](notes/2026-11-13-erdos50-fc-enclosure-kill-test.md) | Erdős #50 | **Negative result** | One natural way to certify f(c) at interior c (exact atom sum + published tail bound) cannot work: the remainder stays ≥ 0.104. |
| [erdos86-mod8-caps](notes/2026-11-18-erdos86-mod8-caps.md) | Erdős #86 | Corollary + table | Applies a published congruence (Minamoto, Lemma 12) at every n: the two-point cap drops by 1 at n = 12, 18, 20, 22, 24, with a one-line rule checked for n = 3…600. |

Each note lists its own limits. Where a note says "first" or "not published", that is the
result of the lab's literature search on that date, not a guarantee.

## Layout

- `notes/` - the 12 notes above.
- `prize_work/` - the scripts and output files behind each note, plus the lab's working
  log for each problem (long, unedited).
- `review_notes/` - the AI referee's re-run report for each result.
- [METHOD.md](METHOD.md) - how a result is checked before it gets a note.

## Reproducing a result

Python 3.11+. Most scripts need `numpy`; a few use `scipy`, `sympy` or `pypdf`. Quick example,
under one second:

```
python prize_work/ppl-076/reciprocal_aps.py
```

Some scripts were written for one machine and contain the absolute path
`E:\Projects\MillenniumResearch`. Replace it with the path of your clone. A few large
output files (over 400 KB) and third-party PDFs are not included; the scripts regenerate
the outputs.

## Corrections and retractions

The lab was wrong several times. These drafts were written and then pulled before
publication, and they are not in this repository:

- **Erdős #86, C4-free certificates for Q15 to Q25 (five drafts).** Superseded: Wang
  (arXiv:2611.14906) and Minamoto (arXiv:2603.29127v5) publish exact values for the same
  family, so the lab's lower bounds were no longer new.
- **Erdős #241, a B₃ set with 10 elements.** The lab bracketed N₁₀ in [310, 359], then
  found that N₁₀ = 310 was already settled (T. Oliveira e Silva, erdosproblems.com forum
  thread 241, 21 March 2026). The literature search had missed a forum post.
- **Erdős #50, Rényi dimension spectrum.** Held: the estimate changes with the bin width
  and has not been shown to converge.
- **Separating words, nS(n) through n = 25.** Held: an update is still waiting for review.

Also withheld: 17 notes that claim status errors in a third-party problem catalogue
(not yet checked by a person), and work on a competition that is still running.

## Found a mistake?

Please open an issue. A wrong number, a missed reference or a claim that is already
known are all useful to me, and I will record the correction here.

## Licence

Code: MIT. Text and data: CC BY 4.0. See [LICENSE](LICENSE).

Juwen Lie - [github.com/juw12](https://github.com/juw12)
