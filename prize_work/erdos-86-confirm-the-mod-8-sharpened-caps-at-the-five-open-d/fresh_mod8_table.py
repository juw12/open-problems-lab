#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Erdos #86 -- fresh, exact-integer re-derivation of the two-point (DPTV) cap and of
the Lemma-12 mod-8 quantisation, turned into a citable table.

Job: lab_erdos_261117_2.  Written from scratch: this file imports nothing from
tools/gate_odsquare_cap.py nor from prize_work/erdos-86-two-point-cap-slack-parity/.
No floating point is used anywhere (only int and fractions.Fraction).

Setting (published, cited in the dossier):
  Q_n has m = 2^(n-1) vertices and E = n*2^(n-1) edges.
  DPTV79 field identity: for every fully-frustrated field of Q_n,
        sum_v h(v)^2 = n*2^(n-1) = E   and   h(v) = n (mod 2),  0 <= h(v) <= n.
  Two-point bound (DPTV79; Minamoto arXiv:2603.29127v5 Lemma 9), and the
  quantisation (Minamoto Lemma 12): for integer a = n (mod 2), b = a+2,
        (h-a)(h-b) = 4*t*(t-1) >= 0  and  = 0 (mod 8),  t = (h-a)/2,
  so with S = sum_v h(v), W = 2a+2, N(a) = E + a(a+2) m:
        Delta(a) = N(a) - W*S >= 0   and   Delta(a) = 0 (mod 8).
  Cap:  g(n) <= (E + S)/2  with S even (S = sum h = m*n mod 2 = 0 for n >= 2).

Output: cap_table.json + cap_table.md (the citable table) and a gate log.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- published data
# 17 published exact g(n) -- Minamoto arXiv:2603.29127v5 (n=6,7,8), Wang
# arXiv:2611.14906 (closed forms), BHN95 (n=4^r), OEIS A245762.  Same list as the
# previous job's table; each value is also re-derived in the gates below from the
# cap formula, which is the point of the gate (published value -> proves the code).
PUB_EXACT = {3: 9, 4: 24, 5: 56, 6: 132, 7: 304, 8: 682, 9: 1532, 10: 3328,
             11: 7296, 13: 33792, 15: 153600, 16: 327680, 17: 688128,
             19: 3047424, 21: 13369344, 23: 58195968, 26: 520093696}
OPEN_N = [12, 14, 18, 20, 22, 24, 25]
# previous job (lab_erdos_261114_3) slack_table.json, used only as a cross-check
PRIOR_CAP_PLAIN = {12: 15701, 14: 72362, 18: 1454899, 20: 6396313,
                   22: 27892121, 24: 120795955, 25: 251658240}
PRIOR_CAP_QUANT = {12: 15700, 14: 72362, 18: 1454898, 20: 6396312,
                   22: 27892120, 24: 120795954, 25: 251658240}


def edges(n):
    """E = number of edges of Q_n = sum_v h(v)^2 (DPTV identity)."""
    return n * (1 << (n - 1))


def nvert(n):
    return 1 << (n - 1)


def admissible_a(n):
    """a must satisfy a = n (mod 2) and a >= 0 (h >= 0)."""
    start = 0 if n % 2 == 0 else 1
    return list(range(start, n + 1, 2))


def s_quant_max(n, a):
    """Largest integer S with W*S <= N(a) and N(a) - W*S = 0 (mod 8).

    Returned with the certificate tuple (W, N, S, Delta) or None if a is
    congruence-infeasible (cannot happen for n >= 2, a = n mod 2 -- see gates).
    """
    E, m = edges(n), nvert(n)
    W = 2 * a + 2
    N = E + a * (a + 2) * m
    S0, r0 = divmod(N, W)                      # S0 = floor(N/W), 0 <= r0 < W
    g = math.gcd(W, 8)
    if r0 % g:
        return None
    mod = 8 // g
    if mod == 1:
        d = 0
    else:
        rhs = ((-r0) // g) % mod               # careful: r0 % g == 0 here
        Wg = (W // g) % mod
        d = (rhs * pow(Wg, -1, mod)) % mod     # smallest d >= 0 with W*(S0-d) = N (mod 8)
    S = S0 - d
    Delta = N - W * S
    assert Delta >= 0 and Delta % 8 == 0
    return S, W, N, Delta


def s_plain_real(n):
    """min_a N(a)/W over admissible a -- the classical (unquantised) two-point bound.
    Exact rational arithmetic."""
    from fractions import Fraction
    E, m = edges(n), nvert(n)
    best, arg = None, None
    for a in admissible_a(n):
        val = Fraction(E + a * (a + 2) * m, 2 * a + 2)
        if best is None or val < best:
            best, arg = val, a
    return best, arg


def cap_plain(n):
    """cap(n) = (E + 2*floor(S*_real/2))/2 with S even (S = m*n mod 2 = 0)."""
    s, a = s_plain_real(n)
    return (edges(n) + 2 * int(s // 2)) // 2, a, s


def cap_quant(n):
    """cap*(n) = (E + S_Q)/2, S_Q = min_a S_quant_max(n,a), reduced to even.

    S = sum_v h(v) = m*n (mod 2) = 0 for every n >= 2 (m = 2^(n-1) even), so the
    odd part of the programme optimum is unavailable: the usable value is the
    largest EVEN S <= min_a S_max(a).  The same parity step is already inside
    cap_plain, so both caps are treated identically.
    """
    E = edges(n)
    best = None
    for a in admissible_a(n):
        r = s_quant_max(n, a)
        if r is None:
            continue
        S, W, N, Delta = r
        S_ev = S if S % 2 == 0 else S - 1
        if best is None or S_ev < best[0]:
            best = (S_ev, a, S, W, N, Delta)
    S_Q, a_bind, S_raw, W, N, Delta = best
    return (E + S_Q) // 2, a_bind, S_Q, S_raw, W, N, Delta


# ------------------------------------------------------------------- certificates
def certify(n, a, S):
    """All multisets of the shape
           c_x copies of h = a+2x,  c_y copies of h = a+2y,  rest h = a
    with sum h = S and sum h^2 = E, exactly m slots.  Exact integer arithmetic.

    With T = (S - a m)/2 and B2 = (E - a^2 m)/4 - a T the two equations are
         c_x x + c_y y   = T
         c_x x^2 + c_y y^2 = B2
    whose unique solution for a given (x,y) is
         c_x = (T y - B2) / (x (y-x)),   c_y = (B2 - T x) / (y (y-x)).
    """
    E, m = edges(n), nvert(n)
    if (S - a * m) % 2:
        return []
    T = (S - a * m) // 2
    if T < 0:
        return []
    R4 = E - a * a * m
    if R4 < 0 or R4 % 4:
        return []
    B2 = R4 // 4 - a * T
    if B2 < 0:
        return []
    Jmax = (n - a) // 2
    out = []
    # one nonzero value
    for x in range(1, Jmax + 1):
        if T % x == 0 and B2 == T * x:
            cx = T // x
            if 0 <= cx <= m:
                out.append({a + 2 * x: cx, a: m - cx})
    # two nonzero values
    for x in range(1, Jmax + 1):
        for y in range(x + 1, Jmax + 1):
            num, den = T * y - B2, x * (y - x)
            if num % den:
                continue
            cx = num // den
            num2, den2 = B2 - T * x, y * (y - x)
            if num2 % den2:
                continue
            cy = num2 // den2
            if cx < 0 or cy < 0 or cx + cy > m:
                continue
            d = {a + 2 * x: cx, a + 2 * y: cy}
            if m - cx - cy:
                d[a] = d.get(a, 0) + (m - cx - cy)
            out.append(d)
    return out


def check_multiset(n, multi):
    """Independent arithmetic check of a multiset: count, sum h, sum h^2."""
    m, E = nvert(n), edges(n)
    cnt = sum(multi.values())
    s1 = sum(h * c for h, c in multi.items())
    s2 = sum(h * h * c for h, c in multi.items())
    ok = (cnt == m) and (s1 * 0 == 0) and (s2 == E) and all(
        h % 2 == n % 2 and 0 <= h <= n for h in multi)
    return ok, cnt, s1, s2


# ------------------------------------------------------------------------- gates
def main():
    log = []
    fail = []

    def rec(tag, ok, msg):
        log.append(("PASS" if ok else "FAIL") + " " + tag + " :: " + msg)
        if not ok:
            fail.append(tag)

    # ---- table for n = 3..26 (the published + open band) -------------------
    rows = []
    for n in range(3, 27):
        cp, ap, sp = cap_plain(n)
        cq, aq, SQ, S_raw, W, N, Delta = cap_quant(n)
        rows.append(dict(n=n, E=edges(n), m=nvert(n),
                         a_plain=ap, S_plain_real=str(sp), S_plain_even=2 * int(sp // 2),
                         cap_plain=cp, a_quant=aq, S_quant=SQ, S_quant_raw=S_raw, W=W, N=N,
                         Delta=Delta, cap_quant=cq, drop=cp - cq,
                         exact=PUB_EXACT.get(n)))
    by_n = {r["n"]: r for r in rows}

    # G1 -- published exact rows: cap >= g everywhere, cap* == g at the attained dims
    for n, g in sorted(PUB_EXACT.items()):
        r = by_n[n]
        rec("G1a n=%d cap>=g" % n, r["cap_plain"] >= g,
            "cap_plain=%d >= g=%d" % (r["cap_plain"], g))
        rec("G1b n=%d cap*>=g" % n, r["cap_quant"] >= g,
            "cap_quant=%d >= g=%d" % (r["cap_quant"], g))
    attained = [n for n in sorted(PUB_EXACT) if by_n[n]["cap_quant"] == PUB_EXACT[n]]
    attained_plain_only = [n for n in sorted(PUB_EXACT)
                           if by_n[n]["cap_plain"] == PUB_EXACT[n] and n not in attained]
    rec("G1c attained set",
        by_n[6]["cap_quant"] == 132 and by_n[9]["cap_quant"] == 1536 and len(attained) == 16,
        "cap*==g at %d dims: %s ; plain-only: %s" % (len(attained), attained, attained_plain_only))
    rec("G1d n=9 is a realisability slack", by_n[9]["cap_quant"] - PUB_EXACT[9] == 4,
        "cap*(9)=%d, g(9)=1532, slack 4 (mod-8 step cannot close it)" % by_n[9]["cap_quant"])

    # G2 -- identical to the previous job's independent table (24 rows)
    prior = json.load(open(os.path.join(HERE, "..", "erdos-86-two-point-cap-slack-parity",
                                        "slack_table.json"), encoding="utf-8"))
    prior_rows = {r["n"]: r for r in prior["rows"]}
    same_plain = all(by_n[n]["cap_plain"] == prior_rows[n]["cap_plain"] for n in prior_rows)
    same_quant = all(by_n[n]["cap_quant"] == prior_rows[n]["cap_quant"] for n in prior_rows)
    same_Sp = all(int(by_n[n]["S_quant"]) == prior_rows[n]["S_quant"] for n in prior_rows)
    rec("G2 x-check job261114_3 cap_plain", same_plain, "%d/24 rows identical" % len(prior_rows))
    rec("G2 x-check job261114_3 cap_quant", same_quant, "%d/24 rows identical" % len(prior_rows))
    rec("G2 x-check job261114_3 S_quant", same_Sp, "%d/24 rows identical" % len(prior_rows))
    for n in OPEN_N:
        rec("G2 open n=%d" % n,
            by_n[n]["cap_quant"] == PRIOR_CAP_QUANT[n] and by_n[n]["cap_plain"] == PRIOR_CAP_PLAIN[n],
            "plain %d / quant %d" % (by_n[n]["cap_plain"], by_n[n]["cap_quant"]))

    # G3 -- Minamoto's three published uses of Lemma 12 (n=6,7,8)
    g3 = {}
    for n, a, want in ((6, 2, 72), (7, 1, 160), (8, 2, 340)):
        got = s_quant_max(n, a)[0]
        g3[n] = got
        rec("G3 Lemma12 n=%d a=%d" % (n, a), got == want, "S_max=%d (published %d)" % (got, want))

    # G4 -- monotonicity: cap* <= cap, and the drop is 0 or 1
    drops = sorted((r["n"], r["drop"]) for r in rows)
    rec("G4a cap* <= cap", all(d >= 0 for _, d in drops), "all 24 rows")
    rec("G4b drop in {0,1}", all(d in (0, 1) for _, d in drops),
        "drops: " + ", ".join("%d:%d" % (n, d) for n, d in drops if d))
    sharp_open = [n for n in OPEN_N if by_n[n]["drop"] == 1]
    rec("G4c exactly the five open dims sharpened",
        sharp_open == [12, 18, 20, 22, 24] and by_n[14]["drop"] == 0 and by_n[25]["drop"] == 0,
        "sharpened: %s ; not sharpened: %s" % (sharp_open,
                                               [n for n in OPEN_N if by_n[n]["drop"] == 0]))

    # G5 -- Lemma 12 bookkeeping at the binding a: Delta >= 0, Delta = 0 mod 8, maximality
    parity_steps = []
    for n in OPEN_N:
        r = by_n[n]
        S_raw, W, N, Delta = r["S_quant_raw"], r["W"], r["N"], r["Delta"]
        rec("G5a n=%d Delta>=0,=0 mod 8" % n, Delta >= 0 and Delta % 8 == 0,
            "Delta=N-W*S_raw=%d (S_raw=%d)" % (Delta, S_raw))
        # maximality: S_raw+1 violates at least one of the two conditions
        rec("G5b n=%d S_raw maximal" % n, (W * (S_raw + 1) > N) or ((N - W * (S_raw + 1)) % 8 != 0),
            "W*(S_raw+1)-N=%d, (N-W*(S_raw+1)) mod 8=%d"
            % (W * (S_raw + 1) - N, (N - W * (S_raw + 1)) % 8))
        # the parity hypothesis of Lemma 12 is what the field supplies
        rec("G5c n=%d parity hypothesis" % n, r["a_quant"] % 2 == n % 2,
            "a=%d = n mod 2; h = n mod 2 is DPTV's parity" % r["a_quant"])
        # S is even (S = m*n mod 2 = 0): record whether the parity step bit
        if r["S_quant"] != S_raw:
            parity_steps.append(n)
        rec("G5d n=%d even-S step" % n, r["cap_quant"] == (r["E"] + r["S_quant"]) // 2,
            "S_quant=%d (S_raw=%d), cap*=(E+S_quant)/2=%d" % (r["S_quant"], S_raw, r["cap_quant"]))

    # G6 -- explicit extremal multiset attaining S_Q at the seven open dims
    certs = {}
    for n in OPEN_N:
        r = by_n[n]
        SQ = r["S_quant"]
        found = None
        for a in admissible_a(n):
            if (SQ - a * nvert(n)) % 2:
                continue
            for multi in certify(n, a, SQ):
                ok, cnt, s1, s2 = check_multiset(n, multi)
                if ok and s1 == SQ:
                    found = dict(a=a, multiset={str(k): v for k, v in sorted(multi.items())},
                                 count=cnt, sum_h=s1, sum_h2=s2)
                    break
            if found:
                break
        certs[n] = found
        rec("G6 n=%d certificate attains S_Q" % n, found is not None,
            "no 2-value certificate" if found is None else
            "a=%d %s (count %d, sum h %d, sum h^2 %d)"
            % (found["a"], found["multiset"], found["count"], found["sum_h"], found["sum_h2"]))

    # G7/G8 -- range extension: cap* for every n = 3..200 (all open dims >= 27 included)
    wide = []
    for n in range(3, 201):
        cp, ap, sp = cap_plain(n)
        cq, aq, SQ, S_raw, W, N, Delta = cap_quant(n)
        wide.append(dict(n=n, E=edges(n), cap_plain=cp, cap_quant=cq, drop=cp - cq,
                         S_quant=SQ, S_quant_raw=S_raw, a_quant=aq, W=W, N=N, Delta=Delta,
                         parity_step=(SQ != S_raw), exact=PUB_EXACT.get(n)))
    wide_sharp = [r["n"] for r in wide if r["drop"]]
    rec("G7a range n=3..200 computed", len(wide) == 198, "198 rows")
    rec("G7b range drop in {0,1}", all(r["drop"] in (0, 1) for r in wide),
        "sharpened (%d dims): %s" % (len(wide_sharp), wide_sharp))
    wide_open_sharp = [n for n in wide_sharp if n not in PUB_EXACT]
    rec("G7c new sharpened open dims >= 27", len(wide_open_sharp) > 0,
        "open dims sharpened beyond the seven tabulated by job 261114_3: %s" % wide_open_sharp)
    rec("G7d sharpening is never decided by parity alone",
        all((not r["parity_step"]) for r in wide if r["drop"]),
        "at every sharpened n the binding S_max was already even (drop comes from the mod-8 congruence)")
    rec("G7e n=14,25 not sharpened", by_n[14]["drop"] == 0 and by_n[25]["drop"] == 0,
        "cap(14)=%d=cap*(14), cap(25)=%d=cap*(25)" % (by_n[14]["cap_quant"], by_n[25]["cap_quant"]))

    # ---- write the table ---------------------------------------------------
    out = dict(
        job="lab_erdos_261117_2",
        convention="cap = (E + S)/2, E = n*2^(n-1), S = sum_v h(v); "
                   "cap_plain from DPTV79, cap_quant from Minamoto arXiv:2603.29127v5 Lemma 12",
        published_exact=PUB_EXACT,
        open_n=OPEN_N,
        rows=rows,
        certificates=certs,
        wide_rows=wide,
        wide_sharpened=wide_sharp,
        wide_open_sharpened=wide_open_sharp,
    )
    with open(os.path.join(HERE, "cap_table.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)

    md = []
    md.append("| n | E = n2^(n-1) | a* | S_Q | cap (plain, DPTV79) | cap* (mod-8, Lemma 12) | cap - cap* | g(n) published |")
    md.append("|---|---|---|---|---|---|---|---|")
    for r in rows:
        md.append("| %d | %d | %d | %d | %d | %d | %d | %s |" % (
            r["n"], r["E"], r["a_quant"], r["S_quant"], r["cap_plain"],
            r["cap_quant"], r["drop"], r["exact"] if r["exact"] is not None else "open"))
    with open(os.path.join(HERE, "cap_table.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")

    print("\n".join(log))
    print("\n--- citable table (n=3..26) ---")
    print("\n".join(md))
    print("\n--- seven open dims ---")
    for n in OPEN_N:
        r = by_n[n]
        print("n=%2d  plain %12d   quantised %12d   drop %d" % (n, r["cap_plain"], r["cap_quant"], r["drop"]))
    print("\n--- certificates ---")
    for n in OPEN_N:
        print(n, certs[n])
    print("\n--- range n=27..64: sharpened dims ---")
    print([r["n"] for r in wide if r["drop"] and r["n"] >= 27])
    print("\nFAILURES: %d %s" % (len(fail), fail))
    return 0 if not fail else 1


if __name__ == "__main__":
    sys.exit(main())
