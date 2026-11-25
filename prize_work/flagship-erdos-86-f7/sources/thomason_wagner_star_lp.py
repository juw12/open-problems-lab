"""Finite-n bound from Thomason-Wagner (2009) Table 1.

Table 1 of "Bounding the size of square-free subgraphs of the hypercube"
(Discrete Math. 309 (2009) 1730-1735) lists the 33 extreme points of the convex
hull of the star 5-tuples (s1..s5) of C4-free subgraphs of Q5, where
s_r = sum over vertices of C(degree, r).

For a C4-free G in Q_n, every r-star lies in C(n-r, 5-r) sub-5-cubes, and there
are C(n,5) * 2^(n-5) sub-5-cubes. So the average star tuple over the sub-5-cubes
lies in that hull:

    sum_v C(d_v, r) * C(n-r, 5-r) = C(n,5) 2^(n-5) * sum_j lambda_j p_j[r],  r = 1..5.

Maximising the number of edges over degree counts n_d (relaxed to reals, then as
integers) gives an upper bound on f(n) that follows from their published table.
The question for our note: is that bound for n = 8 below 691 / 690?
"""
from math import comb

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

TABLE_1 = [
    (0, 0, 0, 0, 0), (32, 0, 0, 0, 0), (36, 40, 40, 20, 4), (48, 40, 40, 20, 4),
    (48, 52, 40, 20, 4), (58, 64, 56, 24, 4), (60, 60, 40, 20, 4), (60, 64, 56, 24, 4),
    (64, 32, 0, 0, 0), (64, 68, 44, 20, 4), (76, 68, 40, 20, 4), (80, 96, 64, 16, 0),
    (84, 100, 68, 24, 4), (84, 100, 72, 28, 4), (90, 110, 76, 28, 4), (92, 108, 72, 28, 4),
    (96, 112, 64, 16, 0), (96, 116, 69, 17, 0), (96, 116, 71, 20, 1), (96, 116, 76, 28, 4),
    (96, 118, 74, 25, 4), (96, 96, 32, 0, 0), (98, 118, 69, 17, 0), (98, 120, 74, 25, 4),
    (102, 124, 72, 22, 4), (102, 126, 76, 24, 4), (104, 124, 68, 20, 4), (104, 128, 76, 24, 4),
    (104, 130, 76, 18, 0), (106, 133, 79, 23, 3), (108, 136, 79, 22, 3), (110, 140, 80, 20, 2),
    (112, 144, 80, 16, 0),
]
assert len(TABLE_1) == 33


def check_table():
    """Each point must be a possible star tuple shape: s1 even, s1 <= 112 (56 edges)."""
    for point in TABLE_1:
        assert point[0] % 2 == 0 and point[0] <= 112, point
    # Two points we can recompute by hand.
    assert (112, 144, 80, 16, 0) == tuple(16 * comb(3, r) + 16 * comb(4, r) for r in range(1, 6))
    assert (96, 96, 32, 0, 0) == tuple(32 * comb(3, r) for r in range(1, 6))


def build(n, with_chung_q3=False):
    """Variables: n_0..n_n (vertex counts by degree), then 33 lambdas."""
    degrees = list(range(n + 1))
    cubes = comb(n, 5) * 2 ** (n - 5)
    width = len(degrees) + len(TABLE_1)
    rows, rhs = [], []
    for r in range(1, 6):
        row = np.zeros(width)
        for d in degrees:
            row[d] = comb(d, r) * comb(n - r, 5 - r)
        for j, point in enumerate(TABLE_1):
            row[len(degrees) + j] = -cubes * point[r - 1]
        rows.append(row)
        rhs.append(0.0)
    row = np.zeros(width)
    row[: len(degrees)] = 1
    rows.append(row)
    rhs.append(2 ** n)
    row = np.zeros(width)
    row[len(degrees):] = 1
    rows.append(row)
    rhs.append(1.0)
    upper_rows, upper_rhs = [], []
    if with_chung_q3:
        # Chung: a C4-free subgraph of Q3 has at most two 3-stars.
        row = np.zeros(width)
        for d in degrees:
            row[d] = comb(d, 3)
        upper_rows.append(row)
        upper_rhs.append(2 * comb(n, 3) * 2 ** (n - 3))
    objective = np.zeros(width)
    for d in degrees:
        objective[d] = -d / 2.0
    return objective, np.array(rows), np.array(rhs), upper_rows, upper_rhs, len(degrees)


def real_bound(n, with_chung_q3=False):
    objective, a_eq, b_eq, a_ub, b_ub, _ = build(n, with_chung_q3)
    result = linprog(
        objective, A_eq=a_eq, b_eq=b_eq,
        A_ub=np.array(a_ub) if a_ub else None, b_ub=np.array(b_ub) if b_ub else None,
        bounds=(0, None), method="highs",
    )
    assert result.status == 0, result.message
    return -result.fun


def integer_bound(n, with_chung_q3=False):
    objective, a_eq, b_eq, a_ub, b_ub, count = build(n, with_chung_q3)
    constraints = [LinearConstraint(a_eq, b_eq - 1e-7, b_eq + 1e-7)]
    if a_ub:
        constraints.append(LinearConstraint(np.array(a_ub), -np.inf, np.array(b_ub)))
    integrality = np.zeros(len(objective))
    integrality[:count] = 1
    result = milp(objective, constraints=constraints, integrality=integrality, bounds=Bounds(0, np.inf))
    assert result.status == 0, result.message
    degrees = {d: int(round(result.x[d])) for d in range(count) if round(result.x[d])}
    return -result.fun, degrees


if __name__ == "__main__":
    check_table()
    known = {5: "f(5) = 56", 6: "f(6) = 132", 7: "f(7) >= 304", 8: "682 <= f(8) <= 690 (ours)", 9: "", 10: ""}
    for n in range(5, 11):
        real = real_bound(n)
        both = real_bound(n, with_chung_q3=True)
        whole, degrees = integer_bound(n, with_chung_q3=True)
        print(f"n={n}: Q5-star LP {real:.4f}; with Chung Q3 {both:.4f}; integer degrees {whole:.1f} "
              f"-> f({n}) <= {int(whole + 1e-6)}   [{known[n]}]   degrees {degrees}")
