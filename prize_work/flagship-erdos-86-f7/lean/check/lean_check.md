# Lean check of this folder

Every file below was accepted by Lean: no error, no `sorry`, and no axiom beyond the three standard ones (propext, Classical.choice, Quot.sound).

| File | Lean verdict | Theorems | Axioms used | Seconds |
|---|---|---|---|---|
| `Averaging.lean` | checked | 18 of 18 | Classical.choice, Quot.sound, propext | 3.9 |
| `Q8Bound689.lean` | checked | 82 of 82 | Classical.choice, Quot.sound, propext | 21.2 |
| `Q8Parity.lean` | checked | 56 of 56 | Classical.choice, Quot.sound, propext | 8.2 |
| `Stars.lean` | checked | 42 of 42 | Classical.choice, Quot.sound, propext | 5.2 |

Run: https://github.com/juw12/open-problems-lab/actions/runs/45804583610 (commit 7375a52370a430932aa289c65572650ffc99b1a2), Lean and Mathlib as pinned in `project/`.
