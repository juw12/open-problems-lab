# Lean check of this folder

Every file below was accepted by Lean: no error, no `sorry`, and no axiom beyond the three standard ones (propext, Classical.choice, Quot.sound).

| File | Lean verdict | Theorems | Axioms used | Seconds |
|---|---|---|---|---|
| `Averaging.lean` | checked | 18 of 18 | Classical.choice, Quot.sound, propext | 4.1 |
| `Q8Bound689.lean` | checked | 82 of 82 | Classical.choice, Quot.sound, propext | 17.4 |
| `Q8Parity.lean` | checked | 56 of 56 | Classical.choice, Quot.sound, propext | 9.0 |
| `Stars.lean` | checked | 42 of 42 | Classical.choice, Quot.sound, propext | 7.8 |

Run: https://github.com/juw12/open-problems-lab/actions/runs/45603989943 (commit e6b51b3ae5ceae1052e625dde01aaae1ea7cc16b), Lean and Mathlib as pinned in `project/`.
