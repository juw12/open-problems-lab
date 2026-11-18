---
Target: erdosproblems.com — problem 986 (Erdős #986, the s ≥ 5 off-diagonal Ramsey
  sharp-log-exponent frontier), forum comment
Date: 2026-11-11
Review gate: PASS — panel PASS (review_panel/lab_erdos_261111_10.md) + referee PASS
  (review_notes/lab_erdos_261111_10.md). Only panel+referee-passed results are
  eligible for a public draft; this one qualifies as a dated reference /
  reproduction of the published frontier (not a new bound).
Scope: DATED REFERENCE / REPRODUCTION ONLY. No new bound, no proof, no solution,
  no gap-closure claim — the s ≥ 5 frontier (lower exponent 2s−4 vs upper
  exponent s−2) is the published standing result; I independently re-state,
  date, and re-check it.
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# A dated reference card for the s ≥ 5 off-diagonal Ramsey sharp-log-exponent frontier (Erdős #986)

> **A dated reference / reproduction of the published frontier — not a new bound.**
> I independently re-stated the lower and upper sharp-log exponents, dated them,
> and checked that nothing published after the lower bound (2026-05-27) narrows
> the gap. I did not close the gap, and I am not crediting an unverified preprint
> that claims to.

## What I checked

The standing frontier for the off-diagonal Ramsey number r(s,t) as t → ∞, in the
sharp-log-exponent form r(s,t) = t^(s−1) / (log t)^(c(s)+o(1)), for s ≥ 5:

- **Lower exponent 2s−4:** r(s,t) ≥ t^(s−1) / (log t)^(2s−4). Bradač (2026),
  arXiv:2605.28793 (v1, 2026-05-27). The bound actually holds for s ≥ 3 (so it
  already covers Erdős #165 / #166); it is a corollary of a published
  deletion-threshold / dependent-random-choice barrier, not a new method.
- **Upper exponent s−2:** r(s,t) ≤ (1+o(1))·t^(s−1) / (log t)^(s−2).
  Ajtai–Komlós–Szemerédi (1980, *Combinatorica* 1:353–370), with the (1+o(1))
  constant from Li–Rousseau–Zang (2001). This is the standing upper since 1980.
- **At s = 5 the gap is exponent 6 (lower) vs 3 (upper).** The exponent gap
  (2s−4) − (s−2) = s−2 grows with s; the s ≥ 5 regime is the live frontier
  because the *upper* side is open, not because the lower bound is new there.
- **Context (s = 4):** Chen (arXiv:2610.15293) and Lin–Zhang (arXiv:2610.17891)
  prove r(4,t) = t³ / (log t)^(2+o(1)) using the s ≥ 5 lower as an input, so
  s = 4 has reached the s−2 upper up to o(1) — that line is essentially closed
  and is **not** a narrowing of the s ≥ 5 gap.
- **Standing-frontier page:** erdosproblems.com/986 (live) confirms this frontier
  as of 2026-11-11.

## How it was checked (independent re-runs)

- **Gate:** `prize_work/erdos-dated-reference-card-for-the-s-5-off-diagonal-ramsey-l/gate.py`
  — re-run 2026-11-11, **13/13 PASS**, with every live anchor re-fetched (Bradač
  v1 2026-05-27; the absence record; the OpenAI preprint finding; the `sorry`
  Lean stub; erdosproblems.com/986). The stored result (`gate_result.json`)
  re-runs byte-identical.
- **Absence scan (post-cutoff 2026-05-27 → 2026-11-11):** arXiv API, 9 AND-queries,
  **44 post-cutoff hits, 0 that narrow the frontier** (the referee independently
  verified the triage: n_beat = 0, all 44 verdicts non_narrowing).
- **Review:** the panel re-ran the gate (13/13) and re-fetched every live anchor;
  the referee re-ran it fresh and independently verified the 0-narrowing count.

## Honest limits

- This is a **dated reference / reproduction, not a new bound**: it re-states and
  dates the published 2s−4 (lower) and s−2 (upper) exponents and checks that
  nothing published post-cutoff narrows the gap. It does **not** close the gap.
- The gap (exponent 6 vs 3 at s = 5) is open.
- The web route of the absence scan surfaced **two unverified non-arXiv preprints
  (2026-09-24) that claim to close the gap** (exponent s−2 for all s ≥ 5). They
  are not credited: their Lean comparators are `theorem main … := by sorry`
  stubs — no build, no axiom audit, no peer review, no arXiv listing. Because of
  that unverified closure claim, I deliberately **withhold the word "confirmed"**
  (the card stays needs_review) rather than assert a closed frontier.
- Precision: the routes the absence scan actually executed were arXiv, web, and
  primary source (a third "Semantic Scholar" route is listed in the card but was
  not run).

## Provenance

- **Sources:** Bradač (2026) arXiv:2605.28793 (v1 2026-05-27);
  Ajtai–Komlós–Szemerédi (1980) *Combinatorica* 1:353–370; Li–Rousseau–Zang
  (2001); prior s ≥ 5 lowers Spencer (1977) *Theoret. Comput. Sci.* 11:55–64 and
  Bohman–Keevash (2010); in-window s = 4 Chen arXiv:2610.15293, Lin–Zhang
  arXiv:2610.17891; standing-frontier page erdosproblems.com/986.
- **Card + record:** the dated reference card under
  `prize_work/erdos-dated-reference-card-for-the-s-5-off-diagonal-ramsey-l/`
  (gate `gate.py`, triage `triage_notes.json`).
- **Review cross-references:** panel `review_panel/lab_erdos_261111_10.md`,
  referee `review_notes/lab_erdos_261111_10.md` (both PASS as a dated reference).
- **AI-assisted:** yes — produced with AI assistance and independently panel- and
  referee-reviewed (both PASS as a dated reference) before standing behind it.

Happy to share the gate script and the triage record if it is useful.
