"""One parity class at a time: can the Q5 inequality be proved class by class?  No.

Q5 is bipartite: A = vertices of even weight, B = odd weight (16 each).  A square of Q5 has two
A-vertices and two B-vertices, and it lies in H exactly when both A-vertices are joined to both
B-vertices.  So "H is C4-free" says only: two A-vertices have at most one common neighbour in H.

Question (one class): is  sum over v in A of g(d(v)) <= 400,  g = (0, 12, 19, 21, 29, 45)?
If yes, the same would hold for B (x -> x xor 1 swaps the classes) and adding the two would give
12 s1 - 5 s2 + 11 s4 - 20 s5 <= 800  for every C4-free subgraph of Q5.

Answer: NO.  CP-SAT finds a C4-free subgraph whose A-sum is 402 (A-degrees 2 x7, 3 x3, 4 x4, 5 x2;
run of 2026-11-25, `q5_class_bound_cpsat.json`).  So the inequality over both classes cannot be
proved one class at a time: a class that is above 400 must be paid for by the other class.  The
inequality itself is not touched by this - see `q5_inequality_cases.py`.

This script asks a solver for the maximum of the class sum (CP-SAT or HiGHS, chosen on the command
line), writes the graph it found, and re-checks that graph from its edge list alone.
It does not use Table 1 of Thomason-Wagner.

Usage: python q5_class_bound.py <cpsat|highs> <seconds> <out.json> [workers]
"""
import json
import sys
import time
from itertools import combinations

G = [0, 12, 19, 21, 29, 45]
N = 5
A = [v for v in range(32) if bin(v).count("1") % 2 == 0]
NEIGHBOURS = {v: [v ^ (1 << i) for i in range(N)] for v in A}
PAIRS = [(u, v) for u, v in combinations(A, 2) if bin(u ^ v).count("1") == 2]
assert len(A) == 16 and len(PAIRS) == 80


def common(u, v):
    both = sorted(set(NEIGHBOURS[u]) & set(NEIGHBOURS[v]))
    assert len(both) == 2
    return both


def solve_cpsat(seconds, workers):
    from ortools.sat.python import cp_model

    model = cp_model.CpModel()
    x = {(v, b): model.NewBoolVar(f"x_{v}_{b}") for v in A for b in NEIGHBOURS[v]}
    gain = []
    for v in A:
        degree = model.NewIntVar(0, N, f"d_{v}")
        model.Add(degree == sum(x[v, b] for b in NEIGHBOURS[v]))
        value = model.NewIntVar(0, G[N], f"g_{v}")
        model.AddElement(degree, G, value)
        gain.append(value)
    for u, v in PAIRS:
        a, b = common(u, v)
        model.Add(x[u, a] + x[u, b] + x[v, a] + x[v, b] <= 3)
    model.Maximize(sum(gain))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    status = solver.Solve(model)
    return {
        "status": solver.StatusName(status),
        "best_found": solver.ObjectiveValue() if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None,
        "proven_upper_bound": solver.BestObjectiveBound(),
        "edge_list": sorted([v, b] for (v, b), variable in x.items() if solver.Value(variable))
        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None,
    }


def check_graph(edge_list):
    """Re-checks a graph from its edge list alone: edges of Q5, C4-free, the A-sum and B-sum."""
    chosen = {(min(a, b), max(a, b)) for a, b in edge_list}
    assert all(bin(a ^ b).count("1") == 1 and 0 <= a < b < 32 for a, b in chosen)
    for i, j in combinations(range(N), 2):
        for v in range(32):
            if not v & (1 << i) and not v & (1 << j):
                square = [(v, v | (1 << i)), (v | (1 << i), v | (1 << i) | (1 << j)),
                          (v | (1 << j), v | (1 << i) | (1 << j)), (v, v | (1 << j))]
                assert sum(edge in chosen for edge in square) <= 3, "the graph has a 4-cycle"
    degrees = [sum((min(v, v ^ (1 << i)), max(v, v ^ (1 << i))) in chosen for i in range(N)) for v in range(32)]
    even = [v for v in range(32) if bin(v).count("1") % 2 == 0]
    odd = [v for v in range(32) if bin(v).count("1") % 2 == 1]
    return {"edges": len(chosen), "degrees_A": sorted(degrees[v] for v in even),
            "degrees_B": sorted(degrees[v] for v in odd),
            "sum_A": sum(G[degrees[v]] for v in even), "sum_B": sum(G[degrees[v]] for v in odd)}


def solve_highs(seconds):
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix

    edge = {}
    for v in A:
        for b in NEIGHBOURS[v]:
            edge[v, b] = len(edge)
    z = {(v, d): len(edge) + k * (N + 1) + d for k, v in enumerate(A) for d in range(N + 1)}
    width = len(edge) + len(z)
    rows, lower, upper = [], [], []
    for u, v in PAIRS:
        a, b = common(u, v)
        rows.append({edge[u, a]: 1, edge[u, b]: 1, edge[v, a]: 1, edge[v, b]: 1})
        lower.append(-np.inf)
        upper.append(3)
    for v in A:
        rows.append({z[v, d]: 1 for d in range(N + 1)})
        lower.append(1)
        upper.append(1)
        link = {z[v, d]: d for d in range(1, N + 1)}
        link.update({edge[v, b]: -1 for b in NEIGHBOURS[v]})
        rows.append(link)
        lower.append(0)
        upper.append(0)
    matrix = lil_matrix((len(rows), width))
    for r, row in enumerate(rows):
        for column, value in row.items():
            matrix[r, column] = value
    objective = np.zeros(width)
    for (v, d), column in z.items():
        objective[column] = -G[d]
    result = milp(objective, constraints=[LinearConstraint(matrix.tocsr(), np.array(lower), np.array(upper))],
                  integrality=np.ones(width), bounds=Bounds(0, 1), options={"time_limit": seconds, "disp": False})
    return {
        "status": f"{int(result.status)}: {result.message}",
        "best_found": None if result.fun is None else -float(result.fun),
        "proven_upper_bound": None if getattr(result, "mip_dual_bound", None) is None else -float(result.mip_dual_bound),
    }


def main():
    solver_name = sys.argv[1]
    seconds = float(sys.argv[2])
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    started = time.time()
    report = solve_cpsat(seconds, workers) if solver_name == "cpsat" else solve_highs(seconds)
    checked = check_graph(report["edge_list"]) if report.get("edge_list") else None
    report.update({
        "question": "is the sum over one parity class of Q5 of g(d(v)) <= 400 for every C4-free subgraph?",
        "solver": solver_name,
        "seconds": round(time.time() - started, 1),
        "checked_graph": checked,
        "answer": "NO" if checked and checked["sum_A"] >= 401 else
                  "YES" if report["proven_upper_bound"] is not None and report["proven_upper_bound"] < 401 - 1e-6
                  else "NOT DECIDED",
    })
    print(json.dumps(report, indent=1))
    with open(sys.argv[3], "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=1)


if __name__ == "__main__":
    main()
