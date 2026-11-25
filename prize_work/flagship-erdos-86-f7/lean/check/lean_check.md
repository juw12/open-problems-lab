# Lean check of this folder

Not every file is checked - see the table.

| File | Lean verdict | Theorems | Axioms used | Seconds |
|---|---|---|---|---|
| `Averaging.lean` | checked | 18 of 18 | Classical.choice, Quot.sound, propext | 5.4 |
| `Q8Bound689.lean` | failed | 81 of 82 | Classical.choice, Quot.sound, propext | 17.4 |
| `Q8Parity.lean` | failed | 55 of 56 | Classical.choice, Quot.sound, propext | 6.6 |
| `Stars.lean` | failed | 41 of 42 | Classical.choice, Quot.sound, propext | 6.9 |

Run: https://github.com/juw12/open-problems-lab/actions/runs/45602990202 (commit 13abe9d359afef47623efd3d157043c9b7cbf71d), Lean and Mathlib as pinned in `project/`.
