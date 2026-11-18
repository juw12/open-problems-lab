#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Erdos #86 -- leg C: independent audit of the mod-8 sharpened caps at the open dims.

Job lab_erdos_261117_2.  Third, from-scratch implementation of the two-point
programme + Minamoto Lemma 12 quantisation, written to disagree with the other
two legs if either is wrong.  Nothing is imported from
  gate_slack_parity.py   (leg A, job lab_erdos_261114_3), or
  fresh_mod8_table.py    (leg B, this job).
Only the standard library; exact integer / Fraction arithmetic, no floats in any
comparison.

Differences in method from the other legs (deliberate):
  * the plain two-point minimum is found by CROSS-MULTIPLYING integers, not by
    Fraction comparisons;
  * S_max(a) is found by scanning DOWN from N//W until Delta == 0 (mod 8) --
    leg B instead solves the congruence with a modular inverse (CRT);
  * attainment of S_Q is tested by a COMPLETE search over supports of size <= 3
    (1, 2 and 3 distinct h-levels) using an integer Vandermonde solve, where
    leg B fixes a two-level closed form and leg A fixes a 3-point ansatz;
  * the whole apparatus is validated on a REAL field: the committed n=12
    certificate's spins are decoded here and its h-profile is recomputed, so
    cap = (E + sum_U h)/2 is checked to return the certificate's own 15614 and
    Delta(a) is checked to be 0 (mod 8) on a field that exists.

Writes audit_caps.json and prints the gate log + the seven-open-dims block.
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

# Published exact g(n) for the odd-square family (Erdos #86).  Provenance:
# Minamoto arXiv:2603.29127v5 (n=6,7,8), Wang arXiv:2611.14906 (closed forms),
# BHN95 / LNSW26 / OEIS A245762.  Used ONLY as the gate: a value the code must
# reproduce, never as an input to the bound.
PUB_EXACT = {3: 9, 4: 24, 5: 56, 6: 132, 7: 304, 8: 682, 9: 1532, 10: 3328,
             11: 7296, 13: 33792, 15: 153600, 16: 327680, 17: 688128,
             19: 3047424, 21: 13369344, 23: 58195968, 26: 520093696}
OPEN_N = [12, 14, 18, 20, 22, 24, 25]
# job lab_erdos_261114_5 (n=12 bracket) and lab_erdos_261114_3 (slack table):
# cross-checks only, loaded from their committed JSONs below.
BRACKET = os.path.join(ROOT, "prize_work", "erdos-86-n12-publish-grade-bracket-15614-15700",
                       "bracket.json")
LEGB = os.path.join(HERE, "cap_table.json")


def nvert(n):
    return 1 << (n - 1)


def edges(n):
    return n * (1 << (n - 1))


def admissible(n):
    """a >= 0, a = n (mod 2): the DPTV parity class of every local field."""
    return list(range(n % 2, n + 1, 2))


def programme_plain(n):
    """argmin_a N(a)/W(a) by integer cross-multiplication; then floor to even S."""
    E, m = edges(n), nvert(n)
    bestN, bestW = None, None
    for a in admissible(n):
        N, W = E + a * (a + 2) * m, 2 * a + 2
        if bestN is None or N * bestW < bestN * W:      # N/W < bestN/bestW
            bestN, bestW = N, W
    S_floor = bestN // bestW
    S_even = 2 * (S_floor // 2)
    return dict(a=min(a for a in admissible(n) if (edges(n) + a * (a + 2) * nvert(n)) * bestW
                      == bestN * (2 * a + 2)),
                W=bestW, N=bestN, S_floor=S_floor, S_even=S_even,
                cap=(E + S_even) // 2)


def s_max(n, a):
    """Largest integer S with W*S <= N and (N - W*S) = 0 (mod 8): downward scan."""
    E, m = edges(n), nvert(n)
    N, W = E + a * (a + 2) * m, 2 * a + 2
    S = N // W
    while (N - W * S) % 8 != 0:
        S -= 1
    return S, N, W, N - W * S


def programme_quant(n):
    """min_a S_max(a); S = sum_U h is even (m even) so floor to even."""
    per_a = {}
    best = None
    for a in admissible(n):
        S, N, W, D = s_max(n, a)
        per_a[a] = dict(S_max=S, N=N, W=W, Delta=D)
        if best is None or S < best[0]:
            best = (S, a)
    S_Q = 2 * (best[0] // 2)
    return dict(a=best[1], S_Q=S_Q, S_min_raw=best[0], cap=(edges(n) + S_Q) // 2,
                per_a=per_a)


# ---------------------------------------------------------------- attainment
def levels(n):
    """All h values a field may take: h = n mod 2, |h| <= n (DPTV)."""
    return list(range(-n, n + 1, 2)) if n % 2 == 0 else list(range(-n, n + 1, 2))


def support_solve(n, lv, S):
    """Exact counts c on the given level tuple lv with
         sum c = m,  sum c*h = S,  sum c*h^2 = E.
    Size 1 and 2 use the first two equations and CHECK the third; size 3 solves
    the integer Vandermonde system.  Returns a dict {h: count} or None."""
    m, E = nvert(n), edges(n)
    if len(lv) == 1:
        h = lv[0]
        if m * h != S or m * h * h != E:
            return None
        return {h: m}
    if len(lv) == 2:
        x, y = lv
        # c_x + c_y = m ; x c_x + y c_y = S
        num, den = S - m * y, x - y                    # c_x = (S - m y)/(x - y)
        if num % den:
            return None
        cx = num // den
        cy = m - cx
        if cx < 0 or cy < 0:
            return None
        if cx * x * x + cy * y * y != E:
            return None
        return {h: c for h, c in ((x, cx), (y, cy)) if c}
    if len(lv) == 3:
        x, y, z = lv
        M = [[1, 1, 1], [x, y, z], [x * x, y * y, z * z]]
        rhs = [m, S, E]
        det = (y - x) * (z - x) * (z - y)              # Vandermonde
        if det == 0:
            return None
        # Cramer
        def det3(col):
            A = [row[:] for row in M]
            for r in range(3):
                A[r][col] = rhs[r]
            return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                    - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                    + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
        out = {}
        for col, h in enumerate((x, y, z)):
            num = det3(col)
            if num % det:
                return None
            c = num // det
            if c < 0:
                return None
            if c:
                out[h] = c
        return out
    raise ValueError("support size")


def complete_search(n, S, maxsize=3):
    """First exact hit with a support of size <= maxsize (levels in increasing order)."""
    L = levels(n)
    k = len(L)
    for i in range(k):
        r = support_solve(n, (L[i],), S)
        if r:
            return r
    for i in range(k):
        for j in range(i + 1, k):
            r = support_solve(n, (L[i], L[j]), S)
            if r:
                return r
    if maxsize >= 3:
        for i in range(k):
            for j in range(i + 1, k):
                for t in range(j + 1, k):
                    r = support_solve(n, (L[i], L[j], L[t]), S)
                    if r:
                        return r
    return None


def check_multi(n, multi, S):
    m, E = nvert(n), edges(n)
    c = sum(multi.values())
    s1 = sum(h * v for h, v in multi.items())
    s2 = sum(h * h * v for h, v in multi.items())
    return c == m and s1 == S and s2 == E and all(h % 2 == n % 2 and abs(h) <= n
                                                 for h in multi)


# ------------------------------------------------------- the real n=12 field
def decode_field(n, spin_hex, edge_count, order):
    """Rebuild the signing J(x, x^2^d) = (-1)^popcount(x & (2^d - 1)), gauge it by
    the spins, and return (positive edges, h on the even-parity half)."""
    raw = bytes.fromhex(spin_hex)
    V = 1 << n
    bit = (lambda i: (raw[i >> 3] >> (i & 7)) & 1) if order == 0 else \
          (lambda i: (raw[i >> 3] >> (7 - (i & 7))) & 1)
    s = [1 - 2 * bit(i) for i in range(V)]
    pc = [bin(x).count("1") for x in range(V)]
    pos = 0
    h = [0] * V
    for x in range(V):
        for d in range(n):
            y = x ^ (1 << d)
            J = -1 if (pc[x & ((1 << d) - 1)] & 1) else 1
            sg = J * s[x] * s[y]
            if sg > 0 and not ((x >> d) & 1):    # count each undirected edge once
                pos += 1
            h[x] += sg                     # h(x) = sum over the n incident edges
    return pos, h, s


def field_gates(rec, fail):
    path = os.path.join(ROOT, "prize_work", "erdos-86-n12-publish-grade-bracket-15614-15700",
                        "cert_n12_15614.json")
    if not os.path.exists(path):
        rec("F0 n=12 certificate present", False, path)
        return None
    cert = json.load(open(path, encoding="utf-8"))
    n = cert["n"]
    got = None
    probes = {}
    for order in (0, 1):
        pos, h, s = decode_field(n, cert["spins"], cert["num_edges"], order)
        probes[order] = pos
        if pos == cert["num_edges"] and got is None:
            got = (order, h, pos)
    rec("F1 the spin encoding decodes to exactly %d positive edges under one bit order"
        % cert["num_edges"], got is not None,
        "probe LSB-first gave %d, MSB-first gave %d -> using order %s"
        % (probes[0], probes[1], got[0] if got else None))
    if got is None:
        return None
    order, h, pos = got
    V, m, E = 1 << n, nvert(n), edges(n)
    even = [v for v in range(V) if bin(v).count("1") % 2 == 0]
    S_U = sum(h[v] for v in even)
    Q_U = sum(h[v] * h[v] for v in even)
    rec("F2 h = n mod 2 and |h| <= n on all %d vertices" % V,
        all(h[v] % 2 == n % 2 and abs(h[v]) <= n for v in range(V)), "")
    rec("F3 field identity sum_U h^2 = n*2^(n-1) = %d" % E, Q_U == E, "got %d" % Q_U)
    rec("F4 cap formula returns the certificate's own value: (E+sum_U h)/2 = %d" % cert["num_edges"],
        (E + S_U) // 2 == cert["num_edges"] and (E + S_U) % 2 == 0,
        "sum_U h = %d -> %d" % (S_U, (E + S_U) // 2))
    q = programme_quant(n)
    rec("F5 the field respects the quantised bound: sum_U h = %d <= S_Q = %d" % (S_U, q["S_Q"]),
        S_U <= q["S_Q"], "")
    D = q["per_a"][q["a"]]["N"] - q["per_a"][q["a"]]["W"] * S_U
    rec("F6 Delta(a*=%d) on a real field is >= 0 and = 0 mod 8" % q["a"],
        D >= 0 and D % 8 == 0, "Delta = %d" % D)
    all_a_ok = all((pp["N"] - pp["W"] * S_U) % 8 == 0 and pp["N"] - pp["W"] * S_U >= 0
                   for pp in q["per_a"].values())
    rec("F7 Delta(a) >= 0 and = 0 mod 8 for EVERY admissible a on that field", all_a_ok,
        "a in %s" % sorted(q["per_a"]))
    return dict(n=n, bit_order=order, positive_edges=pos, sum_U_h=S_U, sum_U_h2=Q_U,
                delta_at_a_star=D, cap_from_field=(E + S_U) // 2)


# --------------------------------------------------------------------- main
def main():
    log, fail = [], []

    def rec(tag, ok, msg):
        log.append(("PASS " if ok else "FAIL ") + tag + (" :: " + msg if msg else ""))
        if not ok:
            fail.append(tag)

    rows = {}
    for n in range(3, 27):
        pl = programme_plain(n)
        qu = programme_quant(n)
        rows[n] = dict(n=n, E=edges(n), m=nvert(n),
                       a_plain=pl["a"], W_plain=pl["W"], N_plain=pl["N"],
                       S_plain_floor=pl["S_floor"], S_plain=pl["S_even"],
                       cap_plain=pl["cap"],
                       a_quant=qu["a"], S_quant=qu["S_Q"], cap_quant=qu["cap"],
                       drop=pl["cap"] - qu["cap"], exact=PUB_EXACT.get(n),
                       per_a=qu["per_a"])

    # A. gate: published exact rows
    bad_hi = [(n, rows[n]["cap_quant"], g) for n, g in sorted(PUB_EXACT.items())
              if rows[n]["cap_quant"] < g]
    rec("A1 cap_quant >= published g at all %d exact dims" % len(PUB_EXACT), not bad_hi,
        str(bad_hi or "none below"))
    att = [n for n in sorted(PUB_EXACT) if rows[n]["cap_quant"] == PUB_EXACT[n]]
    rec("A2 cap_quant == g at the attained dims", len(att) == 16, "attained: %s" % att)
    rec("A3 the two known slacks: n=6 -> 133-132, n=9 -> 1536-1532",
        rows[6]["cap_plain"] == 133 and rows[9]["cap_plain"] == 1536
        and rows[6]["cap_quant"] - PUB_EXACT[6] == 0
        and rows[9]["cap_quant"] - PUB_EXACT[9] == 4,
        "n=6 plain %d quant %d ; n=9 plain %d quant %d"
        % (rows[6]["cap_plain"], rows[6]["cap_quant"], rows[9]["cap_plain"],
           rows[9]["cap_quant"]))

    # B. Lemma 12's own published instances (n = 6,7,8) and its "next-best" pairs
    inst = {(6, 2): 72, (6, 0): 96, (7, 1): 160, (7, 3): 176, (8, 2): 340, (8, 4): 408}
    ok = all(s_max(n, a)[0] == v for (n, a), v in inst.items())
    rec("B1 Minamoto Lemma 12 instances n=6,7,8 (and the next-best a)", ok,
        ", ".join("n=%d a=%d -> %d" % (n, a, s_max(n, a)[0]) for (n, a) in inst))

    # C. maximality of every S_max(a): S+1 must break one of the two conditions
    viol = []
    for n in range(3, 27):
        for a, pp in rows[n]["per_a"].items():
            S, N, W = pp["S_max"], pp["N"], pp["W"]
            if not ((W * (S + 1) > N) or ((N - W * (S + 1)) % 8 != 0)):
                viol.append((n, a))
            if pp["Delta"] < 0 or pp["Delta"] % 8:
                viol.append((n, a, "delta"))
    rec("C1 S_max(a) maximal for every n=3..26 and every admissible a", not viol, str(viol))

    # D. cap* <= cap and the drop is 0 or 1
    rec("D1 cap_quant <= cap_plain and drop in {0,1}",
        all(0 <= r["drop"] <= 1 for r in rows.values()),
        "drops: " + ", ".join("%d:%d" % (n, rows[n]["drop"]) for n in rows if rows[n]["drop"]))
    sharp = [n for n in OPEN_N if rows[n]["drop"] == 1]
    rec("D2 exactly the five open dims 12,18,20,22,24 are sharpened; 14 and 25 are not",
        sharp == [12, 18, 20, 22, 24] and rows[14]["drop"] == 0 and rows[25]["drop"] == 0,
        "sharpened %s" % sharp)

    # E. attainment of S_Q by an explicit multiset (complete search, supports <= 3)
    certs, infeas = {}, {}
    for n in OPEN_N:
        S = rows[n]["S_quant"]
        r = complete_search(n, S)
        ok = bool(r) and check_multi(n, r, S)
        certs[n] = dict(multiset={str(h): c for h, c in sorted(r.items())} if r else None,
                        count=sum(r.values()) if r else 0,
                        sum_h=sum(h * c for h, c in r.items()) if r else 0,
                        sum_h2=sum(h * h * c for h, c in r.items()) if r else 0)
        rec("E1 n=%d: S_Q=%d attained by an explicit multiset" % (n, S), ok,
            str(certs[n]["multiset"]))
        # S_Q + 2 (the next admissible even sum) must have NO support of size <= 3
        r2 = complete_search(n, S + 2)
        infeas[n] = r2 is None
        rec("E2 n=%d: no support of size <= 3 attains S_Q+2 = %d" % (n, S + 2), r2 is None,
            "" if r2 is None else str(r2))

    # F. the real n=12 field
    field = field_gates(rec, fail)

    # G. cross-checks against the other two legs / other jobs
    if os.path.exists(LEGB):
        legb = json.load(open(LEGB, encoding="utf-8"))
        bb = {r["n"]: r for r in legb["rows"]}
        same = all(bb[n]["cap_plain"] == rows[n]["cap_plain"]
                   and bb[n]["cap_quant"] == rows[n]["cap_quant"]
                   and bb[n]["S_quant"] == rows[n]["S_quant"] for n in bb)
        rec("G1 identical to leg B (fresh_mod8_table.py) on all %d rows" % len(bb), same, "")
    else:
        rec("G1 leg B table present", False, LEGB)
    if os.path.exists(BRACKET):
        br = json.load(open(BRACKET, encoding="utf-8"))
        rec("G2 n=12 upper end matches job 261114_5's bracket (15700)",
            br["upper_end"]["value"] == rows[12]["cap_quant"] == 15700,
            "bracket %s, this leg %d" % (br["upper_end"]["value"], rows[12]["cap_quant"]))
        rec("G3 n=12 witness of job 261114_5 agrees with this leg's S_Q",
            br["upper_end"]["ip12_optimum"] == rows[12]["S_quant"] == 6824, "")

    # H. range by-product: sharpened dims n = 3..200 (same band as leg B)
    wide = []
    for n in range(3, 201):
        pl = programme_plain(n)
        qu = programme_quant(n)
        wide.append(dict(n=n, cap_plain=pl["cap"], cap_quant=qu["cap"],
                         drop=pl["cap"] - qu["cap"], S_quant=qu["S_Q"]))
    wide_sharp = [r["n"] for r in wide if r["drop"]]
    rec("H1 range n=3..200: 198 rows, drop in {0,1}",
        len(wide) == 198 and all(r["drop"] in (0, 1) for r in wide),
        "%d sharpened dims" % len(wide_sharp))
    beyond = [n for n in wide_sharp if n > 26]
    rec("H2 sharpened dims beyond the tabulated band (>26)", True, str(beyond))
    if os.path.exists(LEGB):
        legb_wide = legb.get("wide_sharpened")
        if legb_wide is not None:
            rec("H3 identical sharpened-dim list to leg B over n=3..200",
                sorted(legb_wide) == sorted(wide_sharp),
                "leg B: %d dims, this leg: %d dims" % (len(legb_wide), len(wide_sharp)))

    out = dict(job="lab_erdos_261117_2", leg="C (audit_mod8_caps.py)",
               convention="m=|U|=2^(n-1), E=n*2^(n-1), S=sum_U h, cap=(E+S)/2, "
                          "S even for n>=2; cap_plain = DPTV two-point, "
                          "cap_quant = Minamoto arXiv:2603.29127v5 Lemma 12",
               published_exact=PUB_EXACT, open_n=OPEN_N,
               rows=rows, certificates=certs, S_Q_plus_2_infeasible_support_le_3=infeas,
               field_n12=field, wide=wide, wide_sharpened=wide_sharp,
               wide_sharpened_beyond_26=beyond, gates_pass=not fail, failures=fail)
    with open(os.path.join(HERE, "audit_caps.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)

    print("\n".join(log))
    print("\n--- seven open dims (plain DPTV cap vs mod-8 sharpened cap) ---")
    for n in OPEN_N:
        r = rows[n]
        print("n=%2d  a*=%d  S_plain=%6d  S_Q=%6d  plain cap %12d  mod-8 cap %12d  drop %d"
              % (n, r["a_quant"], r["S_plain"], r["S_quant"], r["cap_plain"],
                 r["cap_quant"], r["drop"]))
    print("\n--- extremal multiset per open dim ---")
    for n in OPEN_N:
        print("n=%2d  %s" % (n, certs[n]["multiset"]))
    print("\n--- per-a S_max at n=12 and n=25 (the minimisation) ---")
    for n in (12, 25):
        print("n=%d  %s" % (n, {a: pp["S_max"] for a, pp in rows[n]["per_a"].items()}))
    print("\n--- range by-product ---")
    print("sharpened dims 3..64: %s" % wide_sharp)
    print("sharpened dims > 26 : %s" % beyond)
    print("\nFAILURES: %d %s" % (len(fail), fail))
    return 0 if not fail else 1


if __name__ == "__main__":
    sys.exit(main())
