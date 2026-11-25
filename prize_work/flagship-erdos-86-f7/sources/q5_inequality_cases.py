"""The Q5 inequality as a proof by cases, each case a feasibility question for an exact solver.

Claim (*): every C4-free subgraph H of Q5 has  sum_v g(d(v)) <= 800,  g = (0, 12, 19, 21, 29, 45).
(This is 12 s1 - 5 s2 + 11 s4 - 20 s5 <= 800 with s_r = sum_v C(d(v), r); f(8) <= 689 follows from it
by averaging over the 448 sub-5-cubes of Q8.)  Nothing here uses Table 1 of Thomason-Wagner (2009).

Cases.  Let D be the largest degree of H.  Translations and permutations of the five directions are
automorphisms of Q5, so we may move a vertex of degree D to 0 and its edges to directions 0..D-1.
  D <= 3: the sum is at most 32 * 21 = 672.  Nothing to compute.
  d4:     no vertex has degree 5; vertex 0 has exactly the edges in directions 0, 1, 2, 3.
  d5:     vertex 0 has all five edges.
The directions used at 0 can still be permuted, so the neighbours of 0 may be put in order of
degree: d(e_0) >= d(e_1) >= ...  Each case asks: is there such an H with sum >= 801?
"INFEASIBLE" in both cases proves (*), given the solver.  Case `all` has no case split and no
ordering (slower; fewer things to trust).

Finer split: `top=k` adds d(e_0) = k, the largest degree among the neighbours of 0.  k runs over
1..4 in case d4 and 1..5 in case d5; the nine sub-cases together cover d4 and d5.
Finer still: `nbr=a,b,...` fixes d(e_0) = a, d(e_1) = b, ... (a >= b >= ...; `top=k` is `nbr=k`).
The sub-cases nbr=a,b,k for k = 1..b together cover nbr=a,b, and so on.

Control: with `control` the target is 800 instead of 801; the answer must be FEASIBLE in every case
(the graphs exist), which shows the model is not infeasible by a modelling error.

Rows added to speed the solvers up; each is true for every C4-free subgraph of Q5:
  - same class, distance 2: the two vertices are not both of degree 5 (they would share 2 neighbours);
  - around a vertex u: its H-neighbours have pairwise disjoint sets of other neighbours among the 10
    vertices at distance 2 from u, so  sum_{v in S} (d(v) + 4 x_uv - 5) <= 10  for every set S of
    neighbours of u;
  - per class: sum C(d,2) <= 80 (a pair at distance 2 has at most one common neighbour) and
    sum C(d,3) <= 40 (one 3-star per class in each of the 40 sub-3-cubes).
`plain` leaves all of these out.

Usage: python q5_inequality_cases.py <cpsat|highs> <d4|d5|all> <seconds> <out.json> [workers] [control] [plain] [top=k | nbr=a,b,...]
"""
import json
import sys
import time
from itertools import combinations
from math import comb

N = 5
G = [12 * d - 5 * comb(d, 2) + 11 * comb(d, 4) - 20 * comb(d, 5) for d in range(N + 1)]
assert G == [0, 12, 19, 21, 29, 45]
VERTICES = list(range(2 ** N))
EDGES = [(v, v | (1 << i)) for v in VERTICES for i in range(N) if not v & (1 << i)]
EDGE = {edge: k for k, edge in enumerate(EDGES)}
WIDTH = len(EDGES) + len(VERTICES) * (N + 1)
INF = float("inf")


def x(a, b):
    return EDGE[(min(a, b), max(a, b))]


def z(v, d):
    return len(EDGES) + v * (N + 1) + d


def degree(v, sign=1):
    return {z(v, d): sign * d for d in range(1, N + 1)}


def rows_for(case, target, plain, prefix=()):
    """Rows (coefficients, low, high) over 0/1 variables: 80 edges, then one-hot degrees."""
    rows = []
    for i, j in combinations(range(N), 2):
        for v in VERTICES:
            if v & (1 << i) or v & (1 << j):
                continue
            a, b, c, d = v, v | (1 << i), v | (1 << i) | (1 << j), v | (1 << j)
            rows.append(({x(a, b): 1, x(b, c): 1, x(c, d): 1, x(d, a): 1}, -INF, 3))
    assert len(rows) == 80
    for v in VERTICES:
        rows.append(({z(v, d): 1 for d in range(N + 1)}, 1, 1))
        link = degree(v)
        for i in range(N):
            link[x(v, v ^ (1 << i))] = -1
        rows.append((link, 0, 0))
    rows.append(({z(v, d): G[d] for v in VERTICES for d in range(1, N + 1)}, target, INF))

    if case in ("d4", "d5"):
        top = 4 if case == "d4" else 5
        for i in range(N):
            present = 1 if i < top else 0
            rows.append(({x(0, 1 << i): 1}, present, present))
        if case == "d4":
            for v in VERTICES:
                rows.append(({z(v, 5): 1}, 0, 0))
        for i in range(top - 1):
            order = degree(1 << i)
            order.update(degree(1 << (i + 1), -1))
            rows.append((order, 0, INF))
        assert len(prefix) <= top and all(1 <= value <= top for value in prefix)
        assert all(prefix[i] >= prefix[i + 1] for i in range(len(prefix) - 1))
        for i, value in enumerate(prefix):
            rows.append((degree(1 << i), value, value))
    else:
        assert case == "all" and not prefix

    if not plain:
        for u, v in combinations(VERTICES, 2):
            if bin(u ^ v).count("1") == 2:
                rows.append(({z(u, 5): 1, z(v, 5): 1}, -INF, 1))
        for u in VERTICES:
            for size in range(1, N + 1):
                for dims in combinations(range(N), size):
                    cut = {}
                    for i in dims:
                        v = u ^ (1 << i)
                        cut.update(degree(v))
                        cut[x(u, v)] = 4
                    rows.append((cut, -INF, 10 + 5 * size))
        for parity in (0, 1):
            side = [v for v in VERTICES if bin(v).count("1") % 2 == parity]
            rows.append(({z(v, d): comb(d, 2) for v in side for d in range(2, N + 1)}, -INF, 80))
            rows.append(({z(v, d): comb(d, 3) for v in side for d in range(3, N + 1)}, -INF, 40))
    return rows


def check_graph(values, case, target):
    """Re-checks a solver's graph from its edge list alone: C4-free, the case, the sum."""
    chosen = {edge for edge, k in EDGE.items() if values[k] > 0.5}
    degrees = [sum((min(v, v ^ (1 << i)), max(v, v ^ (1 << i))) in chosen for i in range(N)) for v in VERTICES]
    for i, j in combinations(range(N), 2):
        for v in VERTICES:
            if not v & (1 << i) and not v & (1 << j):
                square = [(v, v | (1 << i)), (v | (1 << i), v | (1 << i) | (1 << j)),
                          (v | (1 << j), v | (1 << i) | (1 << j)), (v, v | (1 << j))]
                assert sum(edge in chosen for edge in square) <= 3, "the solver's graph has a 4-cycle"
    total = sum(G[d] for d in degrees)
    assert total >= target, "the solver's graph misses the target"
    return {"edges": len(chosen), "sum": total, "degrees": sorted(degrees),
            "edge_list": sorted(chosen)}


def solve_cpsat(rows, seconds, workers):
    from ortools.sat.python import cp_model

    model = cp_model.CpModel()
    variables = [model.NewBoolVar(f"b{k}") for k in range(WIDTH)]
    for coefficients, low, high in rows:
        expression = sum(int(value) * variables[column] for column, value in coefficients.items())
        if low == high:
            model.Add(expression == int(low))
        else:
            if low != -INF:
                model.Add(expression >= int(low))
            if high != INF:
                model.Add(expression <= int(high))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    status = solver.Solve(model)
    name = solver.StatusName(status)
    values = [solver.Value(variable) for variable in variables] if name in ("OPTIMAL", "FEASIBLE") else None
    verdict = {"INFEASIBLE": "INFEASIBLE", "OPTIMAL": "FEASIBLE", "FEASIBLE": "FEASIBLE"}.get(name, "UNKNOWN")
    return verdict, name, values, {"branches": solver.NumBranches(), "conflicts": solver.NumConflicts()}


def solve_highs(rows, seconds):
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix

    matrix = lil_matrix((len(rows), WIDTH))
    for r, (coefficients, _, _) in enumerate(rows):
        for column, value in coefficients.items():
            matrix[r, column] = value
    result = milp(np.zeros(WIDTH),
                  constraints=[LinearConstraint(matrix.tocsr(), np.array([row[1] for row in rows]),
                                                np.array([row[2] for row in rows]))],
                  integrality=np.ones(WIDTH), bounds=Bounds(0, 1), options={"time_limit": seconds, "disp": False})
    status = int(result.status)
    verdict = {0: "FEASIBLE", 2: "INFEASIBLE"}.get(status, "UNKNOWN")
    values = list(result.x) if status == 0 and result.x is not None else None
    return verdict, f"{status}: {result.message}", values, {"nodes": int(getattr(result, "mip_node_count", -1) or -1)}


def versions(solver_name):
    found = {"python": sys.version.split()[0]}
    if solver_name == "cpsat":
        import ortools
        found["ortools"] = ortools.__version__
    else:
        import scipy
        found["scipy"] = scipy.__version__
    return found


def main():
    solver_name, case, seconds, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
    flags = sys.argv[5:]
    workers = next((int(flag) for flag in flags if flag.isdigit()), 1)
    target = 800 if "control" in flags else 801
    prefix = next((tuple(int(part) for part in flag[4:].split(",")) for flag in flags
                   if flag.startswith(("top=", "nbr="))), ())
    top_degree = prefix[0] if len(prefix) == 1 else None
    rows = rows_for(case, target, "plain" in flags, prefix)
    started = time.time()
    if solver_name == "cpsat":
        verdict, status, values, effort = solve_cpsat(rows, seconds, workers)
    else:
        verdict, status, values, effort = solve_highs(rows, seconds)
    report = {
        "question": f"is there a C4-free subgraph of Q5 in case {case}"
                    + (f", neighbour degrees of 0 start {list(prefix)}," if prefix else "") + f" with sum g(d) >= {target}?",
        "solver": solver_name, "case": case, "top_degree": top_degree, "nbr": list(prefix), "target": target,
        "plain": "plain" in flags, "workers": workers,
        "verdict": verdict, "status": status, "rows": len(rows), "seconds": round(time.time() - started, 1),
        "time_limit": seconds, **effort, "versions": versions(solver_name),
    }
    if values is not None:
        report["graph"] = check_graph(values, case, target)
    print(json.dumps({k: v for k, v in report.items() if k != "graph"} |
                     ({"graph_sum": report["graph"]["sum"], "graph_degrees": report["graph"]["degrees"]}
                      if "graph" in report else {}), indent=1))
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=1)


if __name__ == "__main__":
    main()
