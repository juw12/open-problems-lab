"""Exact certificate that Table 1 of Thomason-Wagner (2009) gives f(n) <= B for small n.

Setting. G is a C4-free subgraph of Q_n, s_r = sum over all vertices of C(d(v), r).
Each r-star lies in m_r = C(n-r, 5-r) sub-5-cubes and there are N = C(n,5) 2^(n-5)
of them, so the point (m_r s_r / N), r = 1..5, lies in the convex hull of the 33
points of Table 1.

Certificate. Rational y_1..y_5 and alpha, beta with
    (A)  d/2 <= alpha + sum_r y_r m_r C(d, r)      for every degree d = 0..n
    (B)  N * sum_r y_r p[r] <= beta                for every one of the 33 points p
Then  edges = sum_v d(v)/2 <= 2^n alpha + sum_r y_r m_r s_r <= 2^n alpha + beta.

The script finds y with a float LP, rounds it to fractions, takes the smallest
alpha and beta that make (A) and (B) true, and re-checks every row in exact fractions.
"""
import json
import sys
from fractions import Fraction
from math import comb, floor

import numpy as np
from scipy.optimize import linprog

from thomason_wagner_star_lp import TABLE_1


def rows_for(n):
    """Rows over the variables (y1..y5, alpha, beta); each row means row . x >= rhs."""
    cubes = comb(n, 5) * 2 ** (n - 5)
    mult = [comb(n - r, 5 - r) for r in range(1, 6)]
    rows, rhs, names = [], [], []
    for d in range(n + 1):
        rows.append([mult[r - 1] * comb(d, r) for r in range(1, 6)] + [1, 0])
        rhs.append(Fraction(d, 2))
        names.append(f"degree {d}")
    for index, point in enumerate(TABLE_1):
        rows.append([-cubes * point[r] for r in range(5)] + [0, 1])
        rhs.append(Fraction(0))
        names.append(f"s{index} {point}")
    return rows, rhs, names, cubes, mult


def exact_bound(n, y):
    """For rational y, the best alpha and beta are maxima over the rows; always a valid certificate."""
    rows, rhs, names, cubes, mult = rows_for(n)
    alpha = max(Fraction(d, 2) - sum(y[r] * mult[r] * comb(d, r + 1) for r in range(5)) for d in range(n + 1))
    beta = max(cubes * sum(y[r] * point[r] for r in range(5)) for point in TABLE_1)
    return 2 ** n * alpha + beta, alpha, beta


def certificate(n):
    rows, rhs, names, cubes, mult = rows_for(n)
    objective = [0, 0, 0, 0, 0, 2 ** n, 1]
    result = linprog(
        np.array(objective, dtype=float),
        A_ub=-np.array(rows, dtype=float), b_ub=-np.array([float(v) for v in rhs]),
        bounds=[(None, None)] * 7, method="highs",
    )
    assert result.status == 0, result.message
    best = None
    for limit in (10, 20, 50, 100, 200, 500, 1000, 5000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 8):
        y = [Fraction(float(v)).limit_denominator(limit) for v in result.x[:5]]
        value, alpha, beta = exact_bound(n, y)
        if best is None or value < best[0]:
            best = (value, y, alpha, beta)
        if abs(float(value) - result.fun) < 1e-9:
            break
    value, y, alpha, beta = best
    # Exact re-check of every row, independent of how y was found.
    for d in range(n + 1):
        assert Fraction(d, 2) <= alpha + sum(y[r] * mult[r] * comb(d, r + 1) for r in range(5))
    for point in TABLE_1:
        assert cubes * sum(y[r] * point[r] for r in range(5)) <= beta
    tight_degrees = [d for d in range(n + 1)
                     if Fraction(d, 2) == alpha + sum(y[r] * mult[r] * comb(d, r + 1) for r in range(5))]
    tight_points = [f"s{j}" for j, point in enumerate(TABLE_1)
                    if cubes * sum(y[r] * point[r] for r in range(5)) == beta]
    return {
        "n": n,
        "sub_5_cubes": cubes,
        "multiplicity_m_r": mult,
        "y": [str(v) for v in y],
        "alpha": str(alpha),
        "beta": str(beta),
        "bound_exact": str(value),
        "bound_float": float(value),
        "lp_optimum_float": float(result.fun),
        "f_n_at_most": floor(value),
        "tight_rows": [f"degree {d}" for d in tight_degrees] + tight_points,
        "checked": "all degree rows (A) and all 33 table rows (B) hold in exact fractions",
    }


if __name__ == "__main__":
    out = [certificate(n) for n in range(5, 11)]
    for item in out:
        print(f"n={item['n']}: bound {item['bound_exact']} = {item['bound_float']:.4f} -> f({item['n']}) <= {item['f_n_at_most']}")
        print(f"   y = {item['y']}  alpha = {item['alpha']}  beta = {item['beta']}")
        print(f"   tight: {item['tight_rows']}")
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w", encoding="utf-8") as handle:
            json.dump(out, handle, indent=1)
