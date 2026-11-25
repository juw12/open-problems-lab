"""Independent check of one inequality on the 5-cube, by integer programming.

Claim (read off Table 1 of Thomason-Wagner 2009, see thomason_wagner_certificate.py):
for every C4-free subgraph H of Q5, with s_r = sum_v C(d(v), r),

    12 s1 - 5 s2 + 11 s4 - 20 s5 <= 800.

Per vertex this is g(d) = 12 d - 5 C(d,2) + 11 C(d,4) - 20 C(d,5) = 0, 12, 19, 21, 29, 45
for d = 0..5.  This script maximises sum_v g(d(v)) over C4-free subgraphs of Q5 with
HiGHS, without using Table 1.  If the proven upper bound is below 801, the claim holds
(subject to trusting the solver), and f(8) <= 689 follows by averaging over sub-5-cubes.

Usage: python q5_inequality_ilp.py <time_limit_seconds> <out.json> [cutoff] [nosym]
  cutoff: ask only for a graph with gain >= 801; status 2 (infeasible) proves the claim.
  nosym:  switch off the symmetry-breaking rows (slower, fewer things to trust).
"""
import json
import sys
import time
from itertools import combinations
from math import comb

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

N = 5
SYMMETRY = "nosym" not in sys.argv
CUTOFF = "cutoff" in sys.argv
GAIN = [12 * d - 5 * comb(d, 2) + 11 * comb(d, 4) - 20 * comb(d, 5) for d in range(N + 1)]
assert GAIN == [0, 12, 19, 21, 29, 45]

vertices = list(range(2 ** N))
edges = [(v, v | (1 << i)) for v in vertices for i in range(N) if not v & (1 << i)]
edge_index = {edge: k for k, edge in enumerate(edges)}
assert len(edges) == 80


def edge_id(a, b):
    return edge_index[(min(a, b), max(a, b))]


cycles = []
for i, j in combinations(range(N), 2):
    for v in vertices:
        if v & (1 << i) or v & (1 << j):
            continue
        a, b, c, d = v, v | (1 << i), v | (1 << i) | (1 << j), v | (1 << j)
        cycles.append([edge_id(a, b), edge_id(b, c), edge_id(c, d), edge_id(d, a)])
assert len(cycles) == 80

# Variables: x_e (80), then z_{v,d} (32 * 6).
n_x = len(edges)
n_z = len(vertices) * (N + 1)
width = n_x + n_z


def z(v, d):
    return n_x + v * (N + 1) + d


rows = []
lower = []
upper = []


def add(coefficients, lo, hi):
    rows.append(coefficients)
    lower.append(lo)
    upper.append(hi)


for cycle in cycles:
    add({e: 1 for e in cycle}, -np.inf, 3)
for v in vertices:
    add({z(v, d): 1 for d in range(N + 1)}, 1, 1)
    link = {z(v, d): d for d in range(1, N + 1)}
    for i in range(N):
        link[edge_id(v, v ^ (1 << i))] = -1
    add(link, 0, 0)
# Strengthening cuts that are true for every C4-free graph:
# a sub-3-cube holds at most 9 edges, a sub-4-cube at most 24, the whole cube at most 56.
for size, cap in ((3, 9), (4, 24), (5, 56)):
    for dims in combinations(range(N), size):
        mask = sum(1 << i for i in dims)
        for base in vertices:
            if base & mask:
                continue
            inside = {edge_id(a, b): 1 for (a, b) in edges
                      if (a & ~mask) == base and (b & ~mask) == base}
            assert len(inside) == size * 2 ** (size - 1)
            add(inside, -np.inf, cap)

# Pair counting (true for every C4-free graph): two vertices at distance 2 have at most
# one common neighbour.  (a) By parity class: sum over one class of C(d,2) <= 80, and
# sum of C(d,3) <= 40 (one 3-star per class in each of the 40 sub-3-cubes, Chung 1992).
for parity in (0, 1):
    side = [v for v in vertices if bin(v).count("1") % 2 == parity]
    add({z(v, d): comb(d, 2) for v in side for d in range(2, N + 1)}, -np.inf, 80)
    add({z(v, d): comb(d, 3) for v in side for d in range(3, N + 1)}, -np.inf, 40)
# (b) Around one vertex u: the neighbours v of u in H have disjoint sets of other
# neighbours, all among the 10 vertices at distance 2 from u.  So for every set S of
# neighbours, sum_{v in S} (d(v) - 1) [uv in H] <= 10, and (d - 1)[x] >= d + 4x - 5.
for u in vertices:
    for size in range(1, N + 1):
        for dims in combinations(range(N), size):
            cut = {}
            for i in dims:
                v = u ^ (1 << i)
                for d in range(1, N + 1):
                    cut[z(v, d)] = d
                cut[edge_id(u, v)] = 4
            add(cut, -np.inf, 10 + 5 * size)
# Symmetry breaking (the cube's symmetries move any vertex to 0 and permute directions):
# vertex 0 has the largest degree, and its edges use the first directions.
if SYMMETRY:
    for v in vertices[1:]:
        cut = {z(0, d): d for d in range(1, N + 1)}
        for d in range(1, N + 1):
            cut[z(v, d)] = -d
        add(cut, 0, np.inf)
    for i in range(N - 1):
        add({edge_id(0, 1 << i): 1, edge_id(0, 1 << (i + 1)): -1}, 0, np.inf)
# Only look for a counterexample: total gain at least 801.
if CUTOFF:
    add({z(v, d): GAIN[d] for v in vertices for d in range(1, N + 1)}, 801, np.inf)

matrix = lil_matrix((len(rows), width))
for r, coefficients in enumerate(rows):
    for column, value in coefficients.items():
        matrix[r, column] = value

objective = np.zeros(width)
for v in vertices:
    for d in range(N + 1):
        objective[z(v, d)] = -GAIN[d]

time_limit = float(sys.argv[1]) if len(sys.argv) > 1 else 45.0
started = time.time()
result = milp(
    objective,
    constraints=[LinearConstraint(matrix.tocsr(), np.array(lower), np.array(upper))],
    integrality=np.ones(width),
    bounds=Bounds(0, 1),
    options={"time_limit": time_limit, "disp": False},
)
report = {
    "claim": "12 s1 - 5 s2 + 11 s4 - 20 s5 <= 800 for every C4-free subgraph of Q5",
    "status": int(result.status),
    "message": str(result.message),
    "best_found": None if result.fun is None else -float(result.fun),
    "proven_upper_bound": None if getattr(result, "mip_dual_bound", None) is None else -float(result.mip_dual_bound),
    "nodes": int(getattr(result, "mip_node_count", -1) or -1),
    "seconds": round(time.time() - started, 1),
    "time_limit": time_limit,
    "symmetry_breaking": SYMMETRY,
    "cutoff_801": CUTOFF,
}
report["claim_proved_by_solver"] = bool(
    (report["proven_upper_bound"] is not None and report["proven_upper_bound"] < 801 - 1e-6)
    or (CUTOFF and int(result.status) == 2)
)
print(json.dumps(report, indent=1))
if len(sys.argv) > 2:
    with open(sys.argv[2], "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=1)
