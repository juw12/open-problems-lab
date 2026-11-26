# Q5 inequality: HOLDS

Solvers that cover every case at target 801: cadical+drat-trim.  Control (800 is reachable): d4 by cadical+drat-trim, d5 by cadical+drat-trim.

| solver | case | neighbour degrees | target | verdict | seconds | limit |
|---|---|---|---|---|---|---|
| cadical+drat-trim | all plain |  | 801 | UNKNOWN | 9000.0 | 9000 |
| cadical+drat-trim | d4 plain |  | 800 | FEASIBLE | 8.0 | 7200 |
| cadical+drat-trim | d4 plain | 1 | 801 | INFEASIBLE | 3.0 | 7200 |
| cadical+drat-trim | d4 plain | 2 | 801 | INFEASIBLE | 46.6 | 7200 |
| cadical+drat-trim | d4 plain | 3 | 801 | INFEASIBLE | 120.9 | 7200 |
| cadical+drat-trim | d4 plain | 4 | 801 | INFEASIBLE | 378.9 | 7200 |
| cadical+drat-trim | d5 plain |  | 800 | FEASIBLE | 38.7 | 7200 |
| cadical+drat-trim | d5 plain | 4,1,1 | 801 | INFEASIBLE | 31.1 | 7200 |
| cadical+drat-trim | d5 plain | 4,2,1 | 801 | INFEASIBLE | 92.2 | 7200 |
| cadical+drat-trim | d5 plain | 4,2,2 | 801 | INFEASIBLE | 1543.1 | 7200 |
| cadical+drat-trim | d5 plain | 4,3,1 | 801 | INFEASIBLE | 86.0 | 7200 |
| cadical+drat-trim | d5 plain | 4,3,2 | 801 | INFEASIBLE | 1775.8 | 7200 |
| cadical+drat-trim | d5 plain | 4,3,3 | 801 | INFEASIBLE | 2263.9 | 7200 |
| cadical+drat-trim | d5 plain | 4,4,1 | 801 | INFEASIBLE | 111.1 | 7200 |
| cadical+drat-trim | d5 plain | 4,4,2 | 801 | INFEASIBLE | 1278.4 | 7200 |
| cadical+drat-trim | d5 plain | 4,4,3 | 801 | INFEASIBLE | 1300.9 | 7200 |
| cadical+drat-trim | d5 plain | 4,4,4 | 801 | INFEASIBLE | 76.2 | 7200 |
| cadical+drat-trim | d5 plain | 1 | 801 | INFEASIBLE | 19.1 | 7200 |
| cadical+drat-trim | d5 plain | 2 | 801 | INFEASIBLE | 907.4 | 7200 |
| cadical+drat-trim | d5 plain | 3 | 801 | INFEASIBLE | 2118.4 | 7200 |
| cadical+drat-trim | d5 plain | 5 | 801 | INFEASIBLE | 394.4 | 7200 |

Runs: https://github.com/juw12/open-problems-lab/actions/runs/45604846118 (the case split; 18 questions and 2 controls) and https://github.com/juw12/open-problems-lab/actions/runs/45606831102 (the whole problem with no case split: time limit hit, UNKNOWN, decides nothing).
INFEASIBLE here means: CaDiCaL 2.1.3 found no solution and drat-trim verified its DRAT proof against the formula (each JSON file has the checker's last lines, the sizes and the SHA-256 of formula and proof; the proofs, 9.4 GB together, were not kept).
What this table rests on: the encoder `q5_sat_proof.py`, the case split (vertex 0 has the largest degree, its edges go in the first directions, its neighbours are in order of degree) and drat-trim. Not on the solver.
The first run's own last step could not commit these files (a merge conflict with the record of a cancelled run), so they were added by hand from the run's artifacts on 2026-11-26; the table is `q5_verdict.py` over this folder.
