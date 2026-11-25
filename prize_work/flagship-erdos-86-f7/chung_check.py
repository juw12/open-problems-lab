"""Checks on Chung (1992), "Subgraphs of a hypercube containing no small even cycles", Section 4.

Written by the Vector session of 2026-11-25 after reading the paper (scan from the author's site,
sources/chung_1992_small_even_cycles.pdf, journal pages 280-285).

Part A  Edge bound that the *printed* Lemma 11 would give if it were exact for every n:
            sum_v C(d_v,3) + (1/16) sum_v C(n-d_v,3) <= 2 C(n,3) 2^(n-3).
Part B  Goodman's bound in the printed form of Lemma 7, mono >= (1-3p+3p^2) C(t,3), tested on all
        graphs with t <= 6 nodes. (It is an asymptotic form; the exact Goodman bound has (t-1)^2.)
Part C  Lemmas 8, 9, 10 (proofs omitted in the paper) checked on all 4-cycle-free subgraphs of Q3,
        and the per-cube inequality  (#degree-3 nodes) + (1/4)(#admissible nodes) <= 2.
Part D  Both inequalities evaluated on the 682-edge witness in Q8 (q8_odd_square_682.json) if present.
Part E  A 14-edge forest in Q4 for which the printed Lemma 11 fails (so it is not exact for small n).

Standard library only. Run: python chung_check.py
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------- Part A
def max_edges(n: int, weights: list[Fraction]) -> int:
    """Largest e with min over degree sequences of sum weights[d_v] <= 2 C(n,3) 2^(n-3).

    weights is convex in d, so the minimum is at the most equal degree sequence.
    """
    capacity = Fraction(2 * comb(n, 3) * 2 ** (n - 3))
    vertices = 2 ** n
    best = -1
    for edges in range(n * 2 ** (n - 1) + 1):
        low, extra = divmod(2 * edges, vertices)
        value = (vertices - extra) * weights[low]
        if extra:
            value += extra * weights[low + 1]
        if value <= capacity:
            best = edges
    return best


def part_a() -> None:
    print("Part A  edge bounds (n, star capacity only, printed Lemma 11 read as exact, known value)")
    known = {4: "24", 5: "56", 6: "132", 7: "304", 8: ">= 682"}
    for n in range(4, 11):
        star = [Fraction(comb(d, 3)) for d in range(n + 1)]
        lemma11 = [Fraction(comb(d, 3)) + Fraction(comb(n - d, 3), 16) for d in range(n + 1)]
        convex = all(
            lemma11[d + 1] - lemma11[d] >= lemma11[d] - lemma11[d - 1] for d in range(1, n)
        )
        assert convex, n
        print(f"  n={n:2d}  star={max_edges(n, star):5d}  lemma11={max_edges(n, lemma11):5d}"
              f"  known f(n) = {known.get(n, '?')}")
    # the two n = 8 lines used in the notes
    h = lambda d: Fraction(comb(d, 3)) + Fraction(comb(8 - d, 3), 16)
    print("  n=8, 691 edges: 154 nodes of degree 5, 102 of degree 6:",
          float(154 * h(5) + 102 * h(6)), "against capacity 3584")
    print("  n=8, 690 edges: 156 nodes of degree 5, 100 of degree 6:",
          float(156 * h(5) + 100 * h(6)), "against capacity 3584")
    print("  n=8, 691 edges, star part only:", 154 * 10 + 102 * 20, "(slack 4)")


# ---------------------------------------------------------------- Part B
def part_b() -> None:
    print("Part B  printed Lemma 7 (Goodman) on all graphs with t nodes")
    for t in range(3, 7):
        pairs = list(itertools.combinations(range(t), 2))
        triples = list(itertools.combinations(range(t), 3))
        worst = None
        failures = 0
        for mask in range(1 << len(pairs)):
            edge = {pair: (mask >> i) & 1 for i, pair in enumerate(pairs)}
            mono = sum(
                1 for a, b, c in triples
                if edge[(a, b)] == edge[(a, c)] == edge[(b, c)]
            )
            p = Fraction(bin(mask).count("1"), len(pairs))
            printed = (1 - 3 * p + 3 * p * p) * comb(t, 3)
            if mono < printed:
                failures += 1
                gap = printed - mono
                if worst is None or gap > worst[0]:
                    worst = (gap, mono, p)
            # weaker use made in Lemma 11: mono >= C(t,3)/4
        quarter_fail = 0
        for mask in range(1 << len(pairs)):
            edge = {pair: (mask >> i) & 1 for i, pair in enumerate(pairs)}
            mono = sum(
                1 for a, b, c in triples
                if edge[(a, b)] == edge[(a, c)] == edge[(b, c)]
            )
            if 4 * mono < comb(t, 3):
                quarter_fail += 1
        print(f"  t={t}: printed bound fails on {failures} of {1 << len(pairs)} labelled graphs"
              f" (worst: mono={worst[1] if worst else '-'}, p={worst[2] if worst else '-'});"
              f" 'mono >= C(t,3)/4' fails on {quarter_fail}")


# ---------------------------------------------------------------- Part C
Q3_NODES = list(range(8))
Q3_EDGES = [(u, u | (1 << b)) for u in range(8) for b in range(3) if not u & (1 << b)]
Q3_FACES = []
for axis in range(3):
    for value in (0, 1):
        face = [u for u in range(8) if (u >> axis) & 1 == value]
        Q3_FACES.append(frozenset(e for e in Q3_EDGES if e[0] in face and e[1] in face))


def ordered(u: int, v: int) -> tuple[int, int]:
    return (u, v) if u < v else (v, u)


def admissible(graph: frozenset, v: int) -> bool:
    """v has no edge of G in the Q3, and its three pairs of neighbours are all blue or all not blue.

    A pair {u, w} of non-G-neighbours of v is blue when the 4-cycle v, u, x, w has both edges
    u-x and x-w in G (x = the node opposite v in that face).
    """
    bits = (1, 2, 4)
    if any(ordered(v, v ^ b) in graph for b in bits):
        return False
    blue = 0
    for b1, b2 in itertools.combinations(bits, 2):
        u, w, x = v ^ b1, v ^ b2, v ^ b1 ^ b2
        if ordered(u, x) in graph and ordered(w, x) in graph:
            blue += 1
    return blue in (0, 3)


def part_c() -> None:
    print("Part C  Lemmas 8-10 on all 4-cycle-free subgraphs of Q3")
    table: dict[int, int] = {}
    count = 0
    worst_weight = Fraction(0)
    for mask in range(1 << 12):
        graph = frozenset(e for i, e in enumerate(Q3_EDGES) if (mask >> i) & 1)
        if any(face <= graph for face in Q3_FACES):
            continue
        count += 1
        degree = {u: 0 for u in Q3_NODES}
        for a, b in graph:
            degree[a] += 1
            degree[b] += 1
        full = sum(1 for u in Q3_NODES if degree[u] == 3)
        adm = sum(1 for u in Q3_NODES if admissible(graph, u))
        table[full] = max(table.get(full, 0), adm)
        worst_weight = max(worst_weight, full + Fraction(adm, 4))
    print(f"  {count} 4-cycle-free subgraphs of Q3")
    for full in sorted(table):
        print(f"  degree-3 nodes = {full}: most admissible nodes = {table[full]}")
    print("  max of (#degree-3) + (#admissible)/4 over all of them:", worst_weight)
    assert max(table) <= 2 and table.get(2, 0) == 0 and table.get(1, 0) <= 1 and table[0] <= 8
    assert worst_weight <= 2


# ---------------------------------------------------------------- Part D
def part_d() -> None:
    path = HERE / "q8_odd_square_682.json"
    if not path.exists():
        print("Part D  skipped (no witness file)")
        return
    data = json.loads(path.read_text(encoding="utf-8"))
    edges = data.get("edges") if isinstance(data, dict) else data
    if edges is None:
        print("Part D  skipped (witness file has no 'edges' key; keys:", list(data)[:8], ")")
        return
    n = 8
    graph = set()
    for a, b in edges:
        graph.add(ordered(int(a), int(b)))
    degree = [0] * (1 << n)
    for a, b in graph:
        degree[a] += 1
        degree[b] += 1
    star = sum(comb(d, 3) for d in degree)
    printed = star + Fraction(sum(comb(n - d, 3) for d in degree), 16)
    mono = 0
    for v in range(1 << n):
        missing = [1 << b for b in range(n) if ordered(v, v ^ (1 << b)) not in graph]
        for b1, b2, b3 in itertools.combinations(missing, 3):
            blue = 0
            for c1, c2 in ((b1, b2), (b1, b3), (b2, b3)):
                x = v ^ c1 ^ c2
                if ordered(v ^ c1, x) in graph and ordered(v ^ c2, x) in graph:
                    blue += 1
            if blue in (0, 3):
                mono += 1
    capacity = 2 * comb(n, 3) * 2 ** (n - 3)
    print(f"Part D  682-edge witness: edges={len(graph)}, degrees={sorted(set(degree))},"
          f" star sum={star}, printed Lemma 11 left side={float(printed)},"
          f" exact left side (star + mono/4)={float(star + Fraction(mono, 4))}, capacity={capacity}")
    assert star + Fraction(mono, 4) <= capacity

# ---------------------------------------------------------------- Part E
Q4_COUNTEREXAMPLE = [
    (0, 1), (0, 2), (0, 4), (0, 8), (1, 3), (1, 5), (1, 9),
    (6, 14), (7, 15), (10, 14), (11, 15), (12, 14), (13, 15), (14, 15),
]


def part_e() -> None:
    """The printed Lemma 11 is not an exact statement: a 14-edge forest in Q4 violates it.

    Found by chung_lemma11_search.py (n = 4, HiGHS optimum 16.75 against capacity 16); re-checked
    here without the solver.
    """
    n = 4
    graph = {ordered(a, b) for a, b in Q4_COUNTEREXAMPLE}
    assert len(graph) == 14 and all(bin(a ^ b).count("1") == 1 for a, b in graph)
    # 4-cycle-free: no 2-face has all four edges
    for b1, b2 in itertools.combinations(range(n), 2):
        for u in range(1 << n):
            if u & (1 << b1) or u & (1 << b2):
                continue
            a, b, c = u | (1 << b1), u | (1 << b2), u | (1 << b1) | (1 << b2)
            face = [ordered(u, a), ordered(u, b), ordered(a, c), ordered(b, c)]
            assert not all(e in graph for e in face)
    degree = [0] * (1 << n)
    for a, b in graph:
        degree[a] += 1
        degree[b] += 1
    # a_i: fraction of 3-subcubes with i nodes of degree 3 inside the subcube
    counts = {0: 0, 1: 0, 2: 0}
    for fixed in range(n):
        for value in (0, 1):
            nodes = [u for u in range(1 << n) if (u >> fixed) & 1 == value]
            full = 0
            for u in nodes:
                inside = sum(1 for b in range(n) if b != fixed and ordered(u, u ^ (1 << b)) in graph)
                full += inside == 3
            counts[full] += 1
    subcubes = comb(n, 3) * 2 ** (n - 3)
    left = Fraction(counts[1] + 2 * counts[0])            # (a1 + 2 a0) C(n,3) 2^(n-3)
    right = Fraction(sum(comb(n - d, 3) for d in degree), 16)
    print("Part E  Q4 forest with 14 edges, degree histogram",
          {d: degree.count(d) for d in sorted(set(degree))})
    print(f"  3-subcubes by number of degree-3 nodes: {counts} (of {subcubes})")
    print(f"  printed Lemma 11: left side (a1 + 2 a0) C(4,3) 2 = {left}, right side = {right}"
          f"  ->  {'holds' if left >= right else 'FAILS'}")
    assert left < right


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print("DONE")
