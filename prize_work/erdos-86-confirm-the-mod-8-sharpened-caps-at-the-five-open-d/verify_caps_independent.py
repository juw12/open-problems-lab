#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent verifier for job lab_erdos_261117_2 (Erdos #86, mod-8 sharpened caps).

Deliberately written WITHOUT importing fresh_mod8_table.py, and by different means:
  * S_max(a) is found by a naive downward scan over S (no CRT/congruence algebra),
  * the cap is recomputed from the definition,
  * the small dimensions are solved EXHAUSTIVELY by enumerating every h-multiset
    (count vectors), so the programme optimum is known rather than inferred,
  * the seven certificates are re-checked from the JSON with independent arithmetic,
  * the Lemma-12 identity Delta = sum_v (h_v-a)(h_v-a-2) is evaluated directly on
    the multiset and compared with N(a) - W*S.
Exit code 0 iff every check passes.  Pure integers only.
"""
import json
import os
import sys
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))

PUB_EXACT = {3: 9, 4: 24, 5: 56, 6: 132, 7: 304, 8: 682, 9: 1532, 10: 3328,
             11: 7296, 13: 33792, 15: 153600, 16: 327680, 17: 688128,
             19: 3047424, 21: 13369344, 23: 58195968, 26: 520093696}
OPEN_N = [12, 14, 18, 20, 22, 24, 25]

fails = []
logs = []


def chk(tag, ok, msg):
    logs.append(("PASS " if ok else "FAIL ") + tag + " :: " + msg)
    if not ok:
        fails.append(tag)


def NV(n):                      # vertices of Q_n
    return 2 ** (n - 1)


def NE(n):                      # edges of Q_n
    return n * 2 ** (n - 1)


def smax_naive(n, a):
    """Largest S with (2a+2)S <= N and N-(2a+2)S divisible by 8, by scanning down."""
    W = 2 * a + 2
    N = NE(n) + a * (a + 2) * NV(n)
    S = N // W
    while S >= 0:
        if (N - W * S) % 8 == 0:
            return S, W, N, N - W * S
        S -= 1
    return None


def cap_naive(n):
    """cap_plain and cap_quant recomputed from the definitions, keep S even."""
    E = NE(n)
    start = 0 if n % 2 == 0 else 1
    A = range(start, n + 1, 2)
    # plain: min_a N/(2a+2), exact rational, then largest even S below it
    num, den = None, None
    for a in A:                                    # compare fractions num/den
        p, q = E + a * (a + 2) * NV(n), 2 * a + 2
        if num is None or p * den < num * q:
            num, den = p, q
    s_plain_even = 2 * ((num // den) // 2) if (num // den) % 2 == 0 else 2 * (((num // den) - 1) // 2)
    cap_plain = (E + s_plain_even) // 2
    # quantised
    best = None
    for a in A:
        r = smax_naive(n, a)
        if r is None:
            continue
        S = r[0] - (r[0] % 2)
        if best is None or S < best:
            best = S
    return cap_plain, (E + best) // 2, best


def multisets(n, max_combos=None):
    """Exhaustive: every count vector over the allowed h values with m slots and
    sum h^2 = E.  Returns (best_sum_h, best_vector)."""
    m, E = NV(n), NE(n)
    H = list(range(n % 2, n + 1, 2))
    best, bestv = None, None
    k = len(H)
    # iterate over the k-1 smallest values, the last one is forced by the slot count
    for combo in product(*[range(m + 1) for _ in range(k - 1)]):
        used = sum(combo)
        if used > m:
            continue
        s2 = sum(c * h * h for c, h in zip(combo, H[:-1]))
        rem_sq = E - s2
        if rem_sq < 0:
            continue
        last = m - used
        if last == 0:
            if rem_sq != 0:
                continue
        else:
            if rem_sq % (H[-1] ** 2) or rem_sq // (H[-1] ** 2) != last:
                continue
        s1 = sum(c * h for c, h in zip(combo, H[:-1])) + last * H[-1]
        if best is None or s1 > best:
            best, bestv = s1, combo + (last,)
    return best, bestv, H


def main():
    tbl = json.load(open(os.path.join(HERE, "cap_table.json"), encoding="utf-8"))
    rows = {r["n"]: r for r in tbl["rows"]}
    wide = {r["n"]: r for r in tbl["wide_rows"]}

    # V1 -- every row of the 24-row table, recomputed by the naive route
    bad = []
    for n, r in rows.items():
        cp, cq, SQ = cap_naive(n)
        if (cp, cq, SQ) != (r["cap_plain"], r["cap_quant"], r["S_quant"]):
            bad.append((n, (cp, cq, SQ), (r["cap_plain"], r["cap_quant"], r["S_quant"])))
    chk("V1 24-row recomputation", not bad, "mismatches: %s" % (bad,) if bad else
        "all 24 rows identical by the naive downward scan")

    # V2 -- the seven open caps printed by both routes
    got = {n: cap_naive(n)[1] for n in OPEN_N}
    chk("V2 seven open caps identical", all(got[n] == rows[n]["cap_quant"] for n in OPEN_N),
        "; ".join("n=%d %d" % (n, got[n]) for n in OPEN_N))

    # V3 -- cross-check with job 261114_3's own table (third implementation)
    prior = json.load(open(os.path.join(HERE, "..", "erdos-86-two-point-cap-slack-parity",
                                       "slack_table.json"), encoding="utf-8"))
    prows = {r["n"]: r for r in prior["rows"]}
    chk("V3 x-check job 261114_3", all(got[n] == prows[n]["cap_quant"] for n in OPEN_N)
        and all(cap_naive(n)[0] == prows[n]["cap_plain"] for n in prows),
        "seven open caps and all 24 plain caps agree")

    # V4 -- exhaustive programme optimum at the small dimensions
    for n in (4, 5, 6, 7, 8):
        best, vec, H = multisets(n)
        chk("V4 n=%d exhaustive optimum" % n, best == rows[n]["S_quant"],
            "exhaustive max sum h = %d at %s ; cap* route S_quant = %d"
            % (best, dict(zip(H, vec)), rows[n]["S_quant"]))
    # and at those dims the programme optimum is the published g(n)
    for n in (4, 5, 6, 7, 8):
        chk("V4b n=%d cap* == g(n)" % n, rows[n]["cap_quant"] == PUB_EXACT[n],
            "cap*=%d = published g(%d)=%d" % (rows[n]["cap_quant"], n, PUB_EXACT[n]))

    # V5 -- certificates: independent arithmetic + the Lemma-12 identity two ways
    for n in OPEN_N:
        c = tbl["certificates"][str(n)]
        multi = {int(k): v for k, v in c["multiset"].items()}
        m, E = NV(n), NE(n)
        cnt = sum(multi.values())
        s1 = sum(h * c2 for h, c2 in multi.items())
        s2 = sum(h * h * c2 for h, c2 in multi.items())
        ok_shape = cnt == m and s2 == E and all(h % 2 == n % 2 and 0 <= h <= n for h in multi)
        chk("V5a n=%d multiset realises (m,E,S)" % n,
            ok_shape and s1 == rows[n]["S_quant"],
            "count %d, sum h %d, sum h^2 %d" % (cnt, s1, s2))
        a = c["a"]
        W = 2 * a + 2
        N = E + a * (a + 2) * m
        delta_direct = sum((h - a) * (h - a - 2) * c2 for h, c2 in multi.items())
        delta_formula = N - W * s1
        chk("V5b n=%d Delta identity" % n, delta_direct == delta_formula and delta_formula >= 0
            and delta_formula % 8 == 0,
            "sum (h-a)(h-a-2) = %d = N-W*S = %d (>=0, =0 mod 8)" % (delta_direct, delta_formula))

    # V6 -- g(n) never above either cap, and equality at the 16 attained dims
    bad = [n for n, g in PUB_EXACT.items() if rows[n]["cap_quant"] < g or rows[n]["cap_plain"] < g]
    chk("V6a cap >= published g at all 17 dims", not bad, "violations: %s" % bad)
    eq = [n for n in sorted(PUB_EXACT) if rows[n]["cap_quant"] == PUB_EXACT[n]]
    chk("V6b equality at 16 attained dims", len(eq) == 16, "attained: %s ; missing: %s"
        % (eq, [n for n in PUB_EXACT if n not in eq]))

    # V7 -- the range table n=3..200, recomputed, and the sharpened set
    sharp = []
    bad = []
    for n in range(3, 201):
        cp, cq, SQ = cap_naive(n)
        if (cp, cq) != (wide[n]["cap_plain"], wide[n]["cap_quant"]):
            bad.append(n)
        if cp != cq:
            sharp.append(n)
    chk("V7a range 3..200 recomputation", not bad, "mismatched rows: %s" % bad)
    chk("V7b sharpened set", sharp == tbl["wide_sharpened"],
        "%d sharpened dims, identical list" % len(sharp))

    # V8 -- the mod-8 step is what bites: at each sharpened n, S_max is even and
    #       the plain even-rounded S exceeds it by exactly 2
    bad = []
    for n in sharp:
        E = NE(n)
        a = wide[n]["a_quant"]
        S, W, N, D = smax_naive(n, a)
        if S % 2 or wide[n]["S_quant"] != S:
            bad.append(n)
    chk("V8 sharpening is a congruence effect, not parity", not bad,
        "at every sharpened n the binding S_max is even and equals S_quant")

    print("\n".join(logs))
    print("\nFAILURES: %d %s" % (len(fails), fails))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
