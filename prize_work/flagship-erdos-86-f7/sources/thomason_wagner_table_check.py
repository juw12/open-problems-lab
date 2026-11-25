"""Check that every row of Thomason-Wagner Table 1 (as we transcribed it) is the star
tuple of a real degree multiset on the 32 vertices of Q5.

s_r = sum_d n_d C(d, r) is triangular in (n_1..n_5), so the degree counts are
determined by the tuple. A transcription slip would normally give a negative or
impossible count. Also reports the simple necessary conditions per row.
"""
from math import comb

from thomason_wagner_star_lp import TABLE_1

bad = 0
for point in TABLE_1:
    counts = {}
    remaining = list(point)
    for degree in range(5, 0, -1):
        counts[degree] = remaining[degree - 1]
        for r in range(1, degree):
            remaining[r - 1] -= counts[degree] * comb(degree, r)
    counts[0] = 32 - sum(counts.values())
    edges = point[0] // 2
    ok = all(value >= 0 for value in counts.values()) and edges <= 56 and point[2] <= 160
    bad += not ok
    shown = {degree: counts[degree] for degree in range(6) if counts[degree]}
    print(f"{str(point):24} edges {edges:2}  degrees {shown}  {'ok' if ok else 'IMPOSSIBLE'}")
print(f"rows: {len(TABLE_1)}, impossible rows: {bad}")
