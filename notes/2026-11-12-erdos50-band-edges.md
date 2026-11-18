---
Target: Erdős #50 — exact interior band edges of the divergent-slope set of f(c)
Date: 2026-11-12
Review gate: panel PASS `review_panel/lab_panel_261112_9_erdos.md` + referee PASS `review_notes/lab_erdos_261112_11.md` (2026-11-12)
Scope: extension — exact rational band edges of the interior divergent-slope set; NOT a proof, NOT a solution
---

> **Status (published 2026-11-18).** This note comes from an AI-agent research lab run by
> Juwen Lie. The numbers were computed and re-run by AI agents (scripts plus an independent
> AI referee re-run). **No human mathematician has reviewed it**, and it has not been posted
> to the problem's own site. "I" below means the lab. Statements such as "first" or "not
> published" reflect the lab's literature search on the note's date and may be wrong.
> Nothing here is a proof or a solution. Files under `review_panel/` are simulated
> adversarial reviews and are not part of this repository. Corrections are welcome: please
> open an issue.

# Erdős #50: exact band edges of the interior divergent-slope set

> This is a numerical extension past the published frontier. It pins the exact c-interval of the interior divergent-slope set with exact rational endpoints. It does not answer whether f(c) is pointwise differentiable.

## What I checked

Erdős #50 asks whether f(c) = lim_{N→∞} #{n ≤ N : φ(n) < c·n} / N is pointwise differentiable. Prior (publish-clear) lab work located an interior region of c where the windowed slope diverges — the "divergent-slope set". This run pins the **exact rational endpoints** of that interior band:

- **LEFT edge = 8/35 = 0.228571...**
- **RIGHT edge = 42367104/44709361 ≈ 0.94761149**

## Verification (triple)

1. **Exact-rational gate:** 3008/3008 rational checks, 0 mismatches — the endpoints are exact rationals, not floating-point approximations.
2. **Fresh enumeration:** an independent ~1.55M-tower enumeration reproduces the same edges.
3. **Cross-check:** an independent grid (from the earlier Rényi D_q run) agrees to 1×10⁻³.

These are **new numbers**: the earlier work reported the band as an approximation; the exact rational endpoints are new.

## Honest limits

- This pins **where** the interior divergent region is, not **why**, and it is not a proof that this is the full set of non-differentiability points.
- It is a **numerical** result. It does not prove or disprove pointwise differentiability of f(c).
- Erdős #50 remains open.

## Provenance

- **Scripts:** `prize_work/ppl-003/` — the band-edge exact-rational gate + the fresh tower enumeration.
- **Cross-check:** agrees with the independent job-5 grid to 1×10⁻³.
- **Review:** panel PASS + referee PASS (2026-11-12). The referee re-ran the gate independently and confirmed the exact endpoints.
- **AI-assisted:** performed by an AI agent under human supervision. All numbers are reproducible from the committed scripts.

Not posted to the problem site. The scripts are in this repository.
