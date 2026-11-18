---
Target: Erdős #50 — first interior localization of the sharpest divergent-slope region
Date: 2026-11-12
Review gate: panel PASS `review_panel/lab_panel_261112_10_erdos50.md` + referee PASS `review_notes/lab_erdos_261112_13.md` (2026-11-12)
Scope: extension — first interior localization of the divergent-slope set past the Erdős-1995 a.e.-only frontier; NOT a proof, NOT a solution
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# Erdős #50: first interior localization of the sharpest divergent-slope region

> This is a numerical extension. Erdős (1995) established the interior singular measure is a.e. singular but cited no location of the sharpest divergence. This run localizes it to a specific interior cluster.

## What I checked

Inside the interior divergent band (pinned separately: LEFT = 8/35, RIGHT ≈ 0.94761), the sharpest windowed-slope divergence is concentrated in a **cluster near c ≈ 0.6664**:

- cluster c ∈ [0.666417, 0.666499]
- max windowed-slope growth = **70.90**, at set-sizes |S| = 4..7
- the per-|S| divergent fraction is **flat ≈ 0.020 for |S| ≥ 4** — the localization persists as the set-size grows, so it is not a finite-|S| artifact

## Verification

- **j13 gate:** all checks pass. A normalization check gives mean(r) = 0.607927115, matching 6/π² = 0.607927115 to 6 digits (unbiased sieve). 179,446 towers with |S| ≤ 4; total divergent towers 1,289,274.
- **V1 / V2 / V3:** three independent verifiers agree **byte-identically** (max relative difference 0.0).

This is a **new result**: the first localization of *where* the sharpest divergence sits inside the interior band. The prior frontier (Erdős 1995) was a.e.-only, with no location.

## Honest limits

- This is a **localization**, not a proof. It says the sharpest divergence is near c ≈ 0.6664; it does not prove f(c) is non-differentiable there, nor that no other points diverge.
- The cluster interval [0.666417, 0.666499] is a numerical resolution bound, not an exact rational.
- Erdős #50 remains open.

## Provenance

- **Scripts:** `prize_work/ppl-003/` — the j13 gate + three independent verifiers.
- **Cross-check:** V1/V2/V3 byte-identical, max relative difference 0.0.
- **Review:** panel PASS + referee PASS (2026-11-12).
- **AI-assisted:** performed by an AI agent under human supervision. All numbers are reproducible from the committed scripts.

Not posted to the problem site. The scripts are in this repository.
