# Q5 inequality: HOLDS

Solvers that cover every case at target 801: cpsat, highs.  Control (800 is reachable): d4 by cpsat, highs, d5 by cpsat, highs.

| solver | case | top | target | verdict | seconds | limit |
|---|---|---|---|---|---|---|
| cpsat | all plain |  | 801 | UNKNOWN | 19200.4 | 19200 |
| cpsat | all |  | 801 | UNKNOWN | 19200.5 | 19200 |
| cpsat | d4 |  | 800 | FEASIBLE | 1.1 | 1500 |
| cpsat | d4 | 1 | 801 | INFEASIBLE | 0.3 | 19200 |
| cpsat | d4 | 2 | 801 | INFEASIBLE | 1.9 | 19200 |
| cpsat | d4 | 3 | 801 | INFEASIBLE | 1.8 | 19200 |
| cpsat | d4 | 4 | 801 | INFEASIBLE | 25.2 | 19200 |
| cpsat | d4 |  | 801 | INFEASIBLE | 19.7 | 19200 |
| cpsat | d5 |  | 800 | FEASIBLE | 8.1 | 1500 |
| cpsat | d5 | 1 | 801 | INFEASIBLE | 0.4 | 19200 |
| cpsat | d5 | 2 | 801 | INFEASIBLE | 15.2 | 19200 |
| cpsat | d5 | 3 | 801 | INFEASIBLE | 30.0 | 19200 |
| cpsat | d5 | 4 | 801 | INFEASIBLE | 1224.8 | 19200 |
| cpsat | d5 | 5 | 801 | INFEASIBLE | 92.7 | 19200 |
| cpsat | d5 |  | 801 | INFEASIBLE | 1199.0 | 19200 |
| highs | d4 |  | 800 | FEASIBLE | 5.4 | 1500 |
| highs | d4 | 1 | 801 | INFEASIBLE | 0.6 | 19200 |
| highs | d4 | 2 | 801 | INFEASIBLE | 2.5 | 19200 |
| highs | d4 | 3 | 801 | INFEASIBLE | 9.0 | 19200 |
| highs | d4 | 4 | 801 | INFEASIBLE | 159.6 | 19200 |
| highs | d4 |  | 801 | INFEASIBLE | 616.5 | 19200 |
| highs | d5 |  | 800 | FEASIBLE | 87.7 | 1500 |
| highs | d5 | 1 | 801 | INFEASIBLE | 0.4 | 19200 |
| highs | d5 | 2 | 801 | INFEASIBLE | 313.3 | 19200 |
| highs | d5 | 3 | 801 | INFEASIBLE | 369.2 | 19200 |
| highs | d5 | 4 | 801 | INFEASIBLE | 13763.4 | 19200 |
| highs | d5 | 5 | 801 | INFEASIBLE | 299.2 | 19200 |
| highs | d5 |  | 801 | UNKNOWN | 19200.7 | 19200 |

Run: https://github.com/juw12/open-problems-lab/actions/runs/45592147713 (commit 07c42e2eff18a9f4528cb1cd3690cd53eca01e75). A solver verdict, not a formal proof.
