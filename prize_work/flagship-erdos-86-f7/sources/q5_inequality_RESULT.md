# Solver check of the Q5 inequality (independent of Thomason-Wagner's Table 1)

Claim (*): 12 s1 - 5 s2 + 11 s4 - 20 s5 <= 800 for every C4-free subgraph of Q5; the same thing per vertex:
sum_v g(d(v)) <= 800 with g = (0, 12, 19, 21, 29, 45). It gives f(8) <= 689 by averaging over the 448 sub-5-cubes.

## What rests on (*) since 2026-11-25 18:20

`../lean/Q8Bound689.lean`, theorem `bound_8_689_of_Q5Ineq : Q5Ineq -> Bound 8 689`, checked by Lean
(https://github.com/juw12/open-problems-lab/actions/runs/45603989943; record in `../lean/check/lean_check.md`). `Q5Ineq` is (*) word for word: for every C4-free
edge set K of Q5, the sum of g(degree) over the 32 vertices is at most 800. So f(8) <= 689 now needs exactly
this one statement and nothing else. Lean does not check the statement; the runs below do.

## SAT run with checked proofs (finished 2026-11-25 19:14): HOLDS, every case verified

`q5_sat_proof.py` writes each case as a CNF formula, CaDiCaL 2.1.3 answers, and for each "no" drat-trim
checks the solver's DRAT proof (only VERIFIED counts). Workflow `.github/workflows/q5-sat-proof.yml`,
run https://github.com/juw12/open-problems-lab/actions/runs/45604846118.

- **All 18 cases at target 801: INFEASIBLE, and drat-trim says VERIFIED for each proof.** Both controls at
  target 800 are FEASIBLE. The encoder's self-test passed (brute force on Q3, CP-SAT on Q4).
- Sizes: about 5800 s of solving and 9.4 GB of proofs in total. Hardest case, a vertex of degree 5 with
  neighbour degrees 4, 3, 3: 990 s to solve, a 1.47 GB proof, 1274 s to check.
- The proofs were not kept (too large). Each JSON file in `q5_sat_runs/` has the SHA-256 of formula and
  proof, the sizes and the checker's last lines. Table: `q5_sat_runs/q5_sat_verdict.md`.
- The whole problem with no case split and no symmetry (https://github.com/juw12/open-problems-lab/actions/runs/45606831102,
  workflow `q5-sat-proof-all.yml`): no answer in 9000 s. It decides nothing.
- **What this rests on:** the encoder `q5_sat_proof.py`, the case split (vertex 0 has the largest degree, its
  edges go in the first directions, its neighbours are in order of degree; largest degree <= 3 gives at most
  32 x 21 = 672) and drat-trim. Not on the solver. It is not a Lean proof of (*), drat-trim is not a
  formally verified checker, and no person has read the encoder or the case split.
- History: the first run of the workflow (45604007728) failed before solving anything, because drat-trim did
  not compile with its own makefile on the runner (fixed in commit c9f722a). The last job of run 45604846118
  failed on a merge conflict while committing its files; they were added from the run's artifacts in
  commit 7cdf09a.

## Exact solvers (final, 2026-11-25 22:30): HOLDS, by CP-SAT and by HiGHS. Solver verdicts, no proof file.

- **How.** `q5_inequality_cases.py` turns (*) into feasibility questions: is there a C4-free subgraph of Q5 with
  sum >= 801 in this case? Cases: largest degree <= 3 (no run: sum <= 672), `d4` (largest degree 4), `d5` (a vertex
  of degree 5), each cut further by the degrees of the neighbours of that vertex. `q5_verdict.py` reads all the
  JSON files and states one verdict. Every FEASIBLE graph is re-checked from its edge list alone.
- **CP-SAT (OR-Tools 9.15, 4 threads): INFEASIBLE in every case.** d4 whole (20 s) and split in four; d5 whole
  (1199 s), split by d(e_0) = 1..5 (the slow one is d(e_0) = 4: 1225 s), and that one split again in ten.
- **HiGHS (SciPy 1.18): INFEASIBLE in every case of the split** - all of d4 (whole: 617 s) and d5 with
  d(e_0) = 1, 2, 3, 4, 5 (the slow one is d(e_0) = 4: 13763 s). HiGHS on d5 as one question hit the 5.3 h limit
  (UNKNOWN); the five parts cover it.
- **Control.** With target 800 instead of 801 both solvers find a graph in d4 and in d5 (re-checked), so the model
  is not infeasible by a modelling error.
- **No case split: no answer.** CP-SAT on the whole problem with no case split (`all`) and with no case split
  and none of the added rows (`all plain`, the run that trusts the least) both hit the 5.3 h limit: UNKNOWN.
  So the case split is still trusted.
- **Where it ran.** GitHub Actions on the public repository, not this PC (the PC was over its Cairn load limit):
  https://github.com/juw12/open-problems-lab/actions/runs/45592147713 (28 questions, commit 07c42e2) and
  https://github.com/juw12/open-problems-lab/actions/runs/45596181027 (the finer split, commit 17dfd6a).
  Workflows `.github/workflows/q5-inequality.yml` and `q5-inequality-split.yml` (started by hand, no schedule).
  The first run commits its files and verdict to `prize_work/flagship-erdos-86-f7/sources/q5_runs/` of the public
  repository when its last job ends.
- **What it does not give.** No proof file was produced: the verdict is as good as CP-SAT. The case split uses
  two symmetry facts (translations and permutations of directions are automorphisms of Q5; the neighbours of the
  chosen vertex can be put in order of degree) and four families of added rows; each is argued in the docstring
  of `q5_inequality_cases.py` and none has been read by a person. `all plain` would remove both.
- **To repeat on any PC:** `python q5_inequality_cases.py cpsat d4 600 out.json 4` (20 s on 4 threads);
  d5 needs about 20 minutes on 4 threads.

## One parity class: the maximum is 404, not 400

`q5_class_bound.py` asks for the largest sum of g over ONE parity class (16 vertices). CP-SAT: optimum 404, proved
(11 s; class degrees 2 x6, 4 x10; the other class then has 382, total 786; graph re-checked, 52 edges;
`q5_class_bound_cpsat.json`). So (*) cannot be proved one class at a time (a class can be above 400), and the
trivial bound from this is 808.

## History

- 2026-11-25 15:10 gauge run on this PC (`q5_inequality_ilp.py`, HiGHS, 45 s, no cutoff): best found 800, proven
  upper bound 812 (`q5_inequality_ilp_gauge.json`). Not decided.
- 2026-11-25 16:25 cutoff run started on this PC (Windows task "Millennium Q5 inequality check",
  `q5_inequality_ilp_guarded.py`). **Stopped at 17:03 with no verdict**: the PC went over its load limit and the
  run was starved. The task was deleted; `q5_inequality_ilp_guarded.py` is kept only as a record.
- An earlier attempt with CP-SAT and no case split found 800 in 240 s and could not decide 801.
- 2026-11-25 17:08 moved to GitHub Actions with the case split; decided by CP-SAT at 17:30.

## Table at 17:50 (from `q5_verdict.py` over `q5_runs/`; the final table is `q5_runs/q5_verdict.md`)

Solvers that cover every case at target 801: cpsat.  Control (800 is reachable): d4 by cpsat, highs, d5 by cpsat, highs.

| solver | case | neighbour degrees | target | verdict | seconds | limit |
|---|---|---|---|---|---|---|
| cpsat | d4 |  | 800 | FEASIBLE | 1.1 | 1500 |
| cpsat | d4 | 1 | 801 | INFEASIBLE | 0.3 | 19200 |
| cpsat | d4 | 2 | 801 | INFEASIBLE | 1.9 | 19200 |
| cpsat | d4 | 3 | 801 | INFEASIBLE | 1.8 | 19200 |
| cpsat | d4 | 4 | 801 | INFEASIBLE | 25.2 | 19200 |
| cpsat | d4 |  | 801 | INFEASIBLE | 19.7 | 19200 |
| cpsat | d5 |  | 800 | FEASIBLE | 8.1 | 1500 |
| cpsat | d5 | 4,1,1 | 801 | INFEASIBLE | 0.6 | 19200 |
| cpsat | d5 | 4,2,1 | 801 | INFEASIBLE | 0.6 | 19200 |
| cpsat | d5 | 4,2,2 | 801 | INFEASIBLE | 84.0 | 19200 |
| cpsat | d5 | 4,3,1 | 801 | INFEASIBLE | 0.6 | 19200 |
| cpsat | d5 | 4,3,2 | 801 | INFEASIBLE | 280.3 | 19200 |
| cpsat | d5 | 4,3,3 | 801 | INFEASIBLE | 190.9 | 19200 |
| cpsat | d5 | 4,4,1 | 801 | INFEASIBLE | 0.6 | 19200 |
| cpsat | d5 | 4,4,2 | 801 | INFEASIBLE | 334.8 | 19200 |
| cpsat | d5 | 4,4,3 | 801 | INFEASIBLE | 68.6 | 19200 |
| cpsat | d5 | 4,4,4 | 801 | INFEASIBLE | 1.4 | 19200 |
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
| highs | d5 | 4,1,1 | 801 | INFEASIBLE | 0.6 | 19200 |
| highs | d5 | 4,2,1 | 801 | INFEASIBLE | 0.9 | 19200 |
| highs | d5 | 4,3,1 | 801 | INFEASIBLE | 0.8 | 19200 |
| highs | d5 | 4,4,1 | 801 | INFEASIBLE | 0.7 | 19200 |
| highs | d5 | 4,4,4 | 801 | INFEASIBLE | 4.0 | 19200 |
| highs | d5 | 1 | 801 | INFEASIBLE | 0.4 | 19200 |
| highs | d5 | 2 | 801 | INFEASIBLE | 313.3 | 19200 |
| highs | d5 | 3 | 801 | INFEASIBLE | 369.2 | 19200 |
| highs | d5 | 5 | 801 | INFEASIBLE | 299.2 | 19200 |
