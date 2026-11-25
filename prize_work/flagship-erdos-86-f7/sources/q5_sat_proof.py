"""The Q5 inequality as SAT problems whose "no" answers come with a proof a small checker verifies.

Claim (*): every C4-free subgraph H of Q5 has  sum_v g(d(v)) <= 800,  g = (0, 12, 19, 21, 29, 45).

For one case of q5_inequality_cases.py (same case split, same symmetry argument) this script writes
a formula in conjunctive normal form that has a solution exactly when the case holds a C4-free
graph with sum >= target.  A SAT solver (CaDiCaL) then answers; when the answer is "no solution" it
also writes a DRAT proof, and drat-trim - a separate, much smaller program - checks that proof
against the formula.  So for a "no" one trusts this encoder and drat-trim, not the solver.

Variables.
  x_e       edge e is in H                                   (80)
  z_{v,d}   vertex v has degree d                            (32 * 6)
  r_{i,s}   the degrees of vertices 0..i-1 add up to s       (only sums that can still reach target)
Clauses.
  C4:       for each of the 80 squares, not all four edges.
  degree:   for each vertex and each of the 32 patterns of its five edges: pattern -> z_{v,|pattern|};
            and at most one z_{v,d} per vertex.
  sum:      r_{0,0};  r_{i,s} and z_{i,d} -> r_{i+1,s+g(d)}  while s + g(d) < target;
            r_{i,s} and z_{i,d} -> false  when s + g(d) + 45 * (31 - i) < target (target out of reach).
  case:     as in q5_inequality_cases.py (edges at vertex 0, no degree 5 in case d4, neighbour
            degrees in order, the fixed first neighbour degrees).
Why "no solution" proves the case: take a graph of the case with sum >= target, set x and z from
it, and set r_{i,s} true exactly when s is the sum of the first i degrees and s < target.  Every
clause holds (a partial sum of such a graph can always still reach target).  So a formula with no
solution means no such graph.  Only this direction is needed.  The other direction is tested, not
used: `selftest` compares with brute force on Q3 and on random small targets, and the control runs
(target 800) must return a graph, which is re-checked from its edge list alone.

Usage:
  python q5_sat_proof.py selftest
  python q5_sat_proof.py case <d4|d5|all> <out.json> <solve seconds> [control] [top=k | nbr=a,b,...]
The second form needs `cadical` and `drat-trim` on the PATH (the workflow builds both).
"""
import hashlib
import json
import os
import subprocess
import sys
import time
from itertools import combinations, product
from math import comb


# solving and proof checking together stay below this many seconds (one job has 6 hours)
TOTAL_SECONDS = 19000


def weights(n):
    return [12 * d - 5 * comb(d, 2) + 11 * comb(d, 4) - 20 * comb(d, 5) for d in range(n + 1)]


class Cube:
    def __init__(self, n):
        self.n = n
        self.vertices = list(range(2 ** n))
        self.edges = [(v, v | (1 << i)) for v in self.vertices for i in range(n) if not v & (1 << i)]
        self.edge_number = {edge: k + 1 for k, edge in enumerate(self.edges)}
        self.squares = []
        for i, j in combinations(range(n), 2):
            for v in self.vertices:
                if not v & (1 << i) and not v & (1 << j):
                    a, b, c, d = v, v | (1 << i), v | (1 << i) | (1 << j), v | (1 << j)
                    self.squares.append([self.x(a, b), self.x(b, c), self.x(c, d), self.x(d, a)])

    def x(self, a, b):
        return self.edge_number[(min(a, b), max(a, b))]

    def z(self, v, d):
        return len(self.edges) + v * (self.n + 1) + d + 1


def formula(n, g, target, case="all", prefix=()):
    """Returns (number of variables, clauses, cube).  Clauses are lists of non-zero integers."""
    cube = Cube(n)
    clauses = [[-e for e in square] for square in cube.squares]
    for v in cube.vertices:
        incident = [cube.x(v, v ^ (1 << i)) for i in range(n)]
        for pattern in product((0, 1), repeat=n):
            clause = [-e if bit else e for e, bit in zip(incident, pattern)]
            clauses.append(clause + [cube.z(v, sum(pattern))])
        for d1, d2 in combinations(range(n + 1), 2):
            clauses.append([-cube.z(v, d1), -cube.z(v, d2)])

    next_variable = len(cube.edges) + len(cube.vertices) * (n + 1) + 1
    best = max(g)
    count = len(cube.vertices)
    layer = {0: next_variable}
    next_variable += 1
    clauses.append([layer[0]])
    for i in range(count):
        following = {}
        for s, r in sorted(layer.items()):
            for d in range(n + 1):
                s2 = s + g[d]
                if s2 >= target:
                    continue
                if s2 + best * (count - 1 - i) < target:
                    clauses.append([-r, -cube.z(i, d)])
                    continue
                if s2 not in following:
                    following[s2] = next_variable
                    next_variable += 1
                clauses.append([-r, -cube.z(i, d), following[s2]])
        layer = following
    assert not layer, "after the last vertex every sum below target is out of reach"

    def degree_is(v, d):
        return cube.z(v, d)

    if case in ("d4", "d5"):
        assert n == 5
        top = 4 if case == "d4" else 5
        for i in range(n):
            e = cube.x(0, 1 << i)
            clauses.append([e] if i < top else [-e])
        if case == "d4":
            for v in cube.vertices:
                clauses.append([-degree_is(v, 5)])
        for i in range(top - 1):
            for low, high in combinations(range(n + 1), 2):
                clauses.append([-degree_is(1 << i, low), -degree_is(1 << (i + 1), high)])
        assert len(prefix) <= top and all(1 <= value <= top for value in prefix)
        assert all(prefix[i] >= prefix[i + 1] for i in range(len(prefix) - 1))
        for i, value in enumerate(prefix):
            clauses.append([degree_is(1 << i, value)])
    else:
        assert case == "all" and not prefix
    return next_variable - 1, clauses, cube


def graph_from_model(cube, g, model, target):
    """Re-checks a solver's graph from its edge list alone: C4-free and the sum."""
    true = {literal for literal in model if literal > 0}
    chosen = {edge for edge, k in cube.edge_number.items() if k in true}
    for square in cube.squares:
        assert sum(e in true for e in square) <= 3, "the solver's graph has a 4-cycle"
    degrees = [sum((min(v, v ^ (1 << i)), max(v, v ^ (1 << i))) in chosen for i in range(cube.n))
               for v in cube.vertices]
    total = sum(g[d] for d in degrees)
    assert total >= target, "the solver's graph misses the target"
    return {"edges": len(chosen), "sum": total, "degrees": sorted(degrees), "edge_list": sorted(chosen)}


def write_dimacs(path, variables, clauses):
    with open(path, "w", encoding="ascii", newline="\n") as handle:
        handle.write(f"p cnf {variables} {len(clauses)}\n")
        for clause in clauses:
            handle.write(" ".join(map(str, clause)) + " 0\n")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def selftest():
    from pysat.solvers import Cadical153

    def has_solution(n, g, target, case="all", prefix=()):
        variables, clauses, cube = formula(n, g, target, case, prefix)
        with Cadical153(bootstrap_with=clauses) as solver:
            if not solver.solve():
                return None
            return graph_from_model(cube, g, solver.get_model(), target)

    report = {}
    # Q3, every graph by brute force, three weight vectors, every target around the maximum
    cube = Cube(3)
    for g in ([0, 12, 19, 21], [0, 5, 1, 9], [3, 0, 7, 2]):
        best = -1
        for mask in range(2 ** len(cube.edges)):
            chosen = {k + 1 for k in range(len(cube.edges)) if mask >> k & 1}
            if any(all(e in chosen for e in square) for square in cube.squares):
                continue
            degrees = [sum(cube.x(v, v ^ (1 << i)) in chosen for i in range(3)) for v in cube.vertices]
            best = max(best, sum(g[d] for d in degrees))
        for target in range(max(1, best - 3), best + 4):
            found = has_solution(3, g, target)
            assert (found is not None) == (target <= best), (g, target, best)
        report[f"Q3 g={g}"] = f"maximum {best} by brute force; formula agrees at targets {best - 3}..{best + 3}"
    # Q4 with the real weights against a second method (CP-SAT), if it is installed
    try:
        from ortools.sat.python import cp_model
        g = weights(4)
        cube = Cube(4)
        model = cp_model.CpModel()
        edge = {k: model.NewBoolVar(f"x{k}") for k in cube.edge_number.values()}
        for square in cube.squares:
            model.Add(sum(edge[e] for e in square) <= 3)
        total = []
        for v in cube.vertices:
            one_hot = [model.NewBoolVar(f"z{v}_{d}") for d in range(5)]
            model.Add(sum(one_hot) == 1)
            model.Add(sum(d * one_hot[d] for d in range(5)) == sum(edge[cube.x(v, v ^ (1 << i))] for i in range(4)))
            total.append(sum(g[d] * one_hot[d] for d in range(5)))
        model.Maximize(sum(total))
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = 60
        assert solver.StatusName(solver.Solve(model)) == "OPTIMAL"
        best = int(round(solver.ObjectiveValue()))
        assert has_solution(4, g, best) is not None and has_solution(4, g, best + 1) is None
        report["Q4 real weights"] = f"maximum {best} by CP-SAT; formula has a solution at {best}, none at {best + 1}"
    except ImportError:
        report["Q4 real weights"] = "skipped (ortools not installed)"
    # Q5, the real problem: the graphs with sum 800 are found in both cases, and easy cases are refuted
    g = weights(5)
    assert g == [0, 12, 19, 21, 29, 45]
    for case in ("d4", "d5"):
        found = has_solution(5, g, 800, case)
        assert found is not None and found["sum"] == 800, case
        report[f"Q5 {case} target 800"] = f"graph found, sum {found['sum']}, {found['edges']} edges, re-checked"
    for case, prefix in (("d4", (1,)), ("d5", (1,)), ("d5", (4, 1, 1))):
        assert has_solution(5, g, 801, case, prefix) is None
        report[f"Q5 {case} nbr={list(prefix)} target 801"] = "no solution (in-process solver, no proof file)"
    print(json.dumps(report, indent=1))


def run_case(case, out, seconds, flags):
    target = 800 if "control" in flags else 801
    prefix = next((tuple(int(part) for part in flag[4:].split(",")) for flag in flags
                   if flag.startswith(("top=", "nbr="))), ())
    g = weights(5)
    variables, clauses, cube = formula(5, g, target, case, prefix)
    stem = os.path.splitext(out)[0]
    cnf, proof = stem + ".cnf", stem + ".drat"
    write_dimacs(cnf, variables, clauses)
    report = {
        "question": f"is there a C4-free subgraph of Q5 in case {case}"
                    + (f", neighbour degrees of 0 start {list(prefix)}," if prefix else "") + f" with sum g(d) >= {target}?",
        "method": "CNF + CaDiCaL + DRAT proof checked by drat-trim", "case": case, "nbr": list(prefix),
        "top_degree": prefix[0] if len(prefix) == 1 else None, "target": target,
        "variables": variables, "clauses": len(clauses), "cnf_sha256": sha256(cnf),
    }
    started = time.time()
    solved = subprocess.run(["cadical", "-t", str(int(seconds)), cnf, proof], capture_output=True, text=True)
    report["solve_seconds"] = round(time.time() - started, 1)
    report["versions"] = {
        "python": sys.version.split()[0],
        "cadical": subprocess.run(["cadical", "--version"], capture_output=True, text=True).stdout.strip(),
        "drat-trim": os.environ.get("DRAT_TRIM_COMMIT", "unknown commit"),
    }
    answer = next((line[2:].strip() for line in solved.stdout.splitlines() if line.startswith("s ")), "UNKNOWN")
    report["solver_answer"] = answer
    if answer == "SATISFIABLE":
        model = [int(token) for line in solved.stdout.splitlines() if line.startswith("v ")
                 for token in line[2:].split() if token != "0"]
        report["graph"] = graph_from_model(cube, g, model, target)
        report["verdict"] = "FEASIBLE"
    elif answer == "UNSATISFIABLE":
        report["proof_bytes"] = os.path.getsize(proof)
        started = time.time()
        check_limit = int(max(600, TOTAL_SECONDS - report["solve_seconds"]))
        checked = subprocess.run(["drat-trim", cnf, proof, "-t", str(check_limit)], capture_output=True, text=True)
        report["check_limit"] = check_limit
        report["check_seconds"] = round(time.time() - started, 1)
        report["checker_says"] = next((line[2:].strip() for line in checked.stdout.splitlines()
                                       if line.startswith("s ")), "no verdict line")
        report["checker_tail"] = checked.stdout.splitlines()[-6:]
        report["proof_sha256"] = sha256(proof)
        report["verdict"] = "INFEASIBLE" if report["checker_says"] == "VERIFIED" else "UNKNOWN"
        report["proof_checked"] = report["checker_says"] == "VERIFIED"
    else:
        report["verdict"] = "UNKNOWN"
    report["seconds"] = round(report["solve_seconds"] + report.get("check_seconds", 0), 1)
    report["solver"] = "cadical+drat-trim"
    report["status"] = answer + (", proof " + report.get("checker_says", "") if answer == "UNSATISFIABLE" else "")
    report["plain"] = True
    report["time_limit"] = seconds
    printable = {k: v for k, v in report.items() if k != "graph"}
    if "graph" in report:
        printable["graph_sum"] = report["graph"]["sum"]
    print(json.dumps(printable, indent=1))
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=1)
    if os.path.exists(proof):
        os.remove(proof)


def main():
    if sys.argv[1] == "selftest":
        selftest()
    else:
        assert sys.argv[1] == "case"
        run_case(sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5:])


if __name__ == "__main__":
    main()
