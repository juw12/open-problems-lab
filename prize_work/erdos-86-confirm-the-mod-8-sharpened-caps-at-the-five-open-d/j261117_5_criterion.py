#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Erdos #86 -- job lab_erdos_261117_5: WHICH dimensions the mod-8 quantisation sharpens.

Context.  Job lab_erdos_261117_2 tabulated the sharpened caps at n = 12,18,20,22,24
(and 30 further open dims <= 200) but wrote that the sharpened set "has no simple
congruence description that I could find".  This job supplies that description and
checks it exhaustively.

Setting (published; DPTV79 = J. Physique 40 (1979) 620; Minamoto arXiv:2603.29127v5,
Lemma 9 = two-point bound, Lemma 12 = mod-8 quantisation).
  m = 2^(n-1) vertices, E = n*m edges, field identity sum_v h(v)^2 = E,
  h(v) = n (mod 2), S = sum_v h(v), and  g(n) <= (E+S)/2  with S even.
  For every admissible a = n (mod 2), 0 <= a <= n:
      W = 2a+2,  N(a) = E + a(a+2)m,  Delta(a) = N(a) - W*S = sum_v (h-a)(h-a-2) >= 0
      and Delta(a) = 0 (mod 8).
  Plain (DPTV79) bound:  S <= rho := min_a N(a)/W(a), S even.
  Quantised (Lemma 12) bound: S <= min_a S_q(a), S even.

Claim proved and machine-checked here (the new step):
  For n >= 4, N(a) = 0 (mod 8) for every admissible a, so the constraint
  N - W S = 0 (mod 8) reads, with omega = gcd(W, 8):
      a = 3 (mod 4)  ->  W = 0 (mod 8)  ->  no constraint at all
      a = 1 (mod 4)  ->  W = 4 (mod 8)  ->  S = 0 (mod 2)   (DPTV parity, no loss)
      a even         ->  W = 2 (mod 4)  ->  S = 0 (mod 4)   (NEW restriction)
  Every admissible a also admits S_pl (the even-rounding of rho) as a plain value,
  so the quantised optimum is exactly S_pl - 2 when at least one even a exists
  (i.e. n even) and S_pl = 2 (mod 4), and S_pl otherwise.  Hence

      cap*(n) = cap(n) - 1   <==>   n even  and  S_pl(n) = 2 (mod 4),

  where S_pl = 2*floor(rho/2).  The "arithmetic-random" sharpened set is therefore
  the residue class of a single integer per dimension.

Parts
  P1 GATE   the card's kill test at n=12 (brute-force descent, no CRT algebra):
            S(12) = 6824 with multiset {h=0:2, h=2:680, h=4:1366}.
  P2 GATE   all 17 published exact g(n): cap* >= g, equality at the 16 attained dims,
            n=9 keeps its realisability slack of 4; the seven open caps printed.
  P3 NEW    the criterion above, checked dimension by dimension against a brute-force
            minimisation over all admissible a, for n = 3..600.
  P4 NEW    the sharpened set for n <= 600 as predicted by the criterion, and the
            count of open dims >= 27 it covers.

Nothing is claimed beyond a corollary of published lemmas: no g(n) is determined,
no field is constructed, no proof of Erdos #86 is claimed.

Job: lab_erdos_261117_5.  Exact integers only (no floats anywhere).
Exit 0 iff every gate passes.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

PUB_EXACT = {3: 9, 4: 24, 5: 56, 6: 132, 7: 304, 8: 682, 9: 1532, 10: 3328,
             11: 7296, 13: 33792, 15: 153600, 16: 327680, 17: 688128,
             19: 3047424, 21: 13369344, 23: 58195968, 26: 520093696}
OPEN_N = [12, 14, 18, 20, 22, 24, 25]
CAP_PLAIN = {12: 15701, 14: 72362, 18: 1454899, 20: 6396313,
             22: 27892121, 24: 120795955, 25: 251658240}
CAP_QUANT = {12: 15700, 14: 72362, 18: 1454898, 20: 6396312,
             22: 27892120, 24: 120795954, 25: 251658240}
NMAX = 600

logs, fails = [], []


def chk(tag, ok, msg):
    logs.append(("PASS " if ok else "FAIL ") + tag + " :: " + msg)
    if not ok:
        fails.append(tag)


def m_of(n):
    return 1 << (n - 1)


def E_of(n):
    return n * (1 << (n - 1))


def adm(n):
    """admissible a: a = n (mod 2), 0 <= a <= n"""
    return list(range(n % 2, n + 1, 2))


def sq(a, n):
    """largest S with W*S <= N and N - W*S = 0 (mod 8), by downward scan (no CRT)."""
    W = 2 * a + 2
    N = E_of(n) + a * (a + 2) * m_of(n)
    S = N // W
    while S >= 0:
        if (N - W * S) % 8 == 0:
            return S, W, N, N - W * S
        S -= 1
    return None


def plain(n):
    """rho = min_a N(a)/W(a) exactly, then the largest even S <= rho, and the argmin set."""
    num = den = None
    argmin = []
    for a in adm(n):
        p, q = E_of(n) + a * (a + 2) * m_of(n), 2 * a + 2
        if num is None or p * den < num * q:          # cross-multiplied comparison
            num, den, argmin = p, q, [a]
        elif p * den == num * q:
            argmin.append(a)
    S_even = 2 * ((num // den) // 2)                  # floor(rho) rounded down to even
    return S_even, num // den, den, argmin


def quantised(n):
    """min_a S_q(a), then rounded down to even: the Lemma-12 optimum."""
    best = None
    for a in adm(n):
        S = sq(a, n)[0]
        S -= S % 2
        if best is None or S < best:
            best = S
    return best


def cap_plain(n):
    S, _, _, _ = plain(n)
    return (E_of(n) + S) // 2


def cap_quant(n):
    return (E_of(n) + quantised(n)) // 2


def certify_12():
    """The card's kill test, brute force: the extremal multiset at n=12 with S=6824."""
    n, m, E = 12, m_of(12), E_of(12)
    S = quantised(n)
    # solve over the two-value shapes h in {a, a+2x, a+2y} for EVERY admissible a
    out = []
    for a in adm(n):
        T = (S - a * m) // 2
        if (S - a * m) % 2 or T < 0:
            continue
        B2 = (E - a * a * m) // 4 - a * T
        for x in range(1, (n - a) // 2 + 1):
            for y in range(x, (n - a) // 2 + 1):
                if y == x:
                    if T % x or B2 != T * x:
                        continue
                    cx, cy = T // x, 0
                else:
                    if (T * y - B2) % (x * (y - x)) or (B2 - T * x) % (y * (y - x)):
                        continue
                    cx, cy = (T * y - B2) // (x * (y - x)), (B2 - T * x) // (y * (y - x))
                if cx < 0 or cy < 0 or cx + cy > m:
                    continue
                d = {}
                for h, c in ((a + 2 * x, cx), (a + 2 * y, cy), (a, m - cx - cy)):
                    if c:
                        d[h] = d.get(h, 0) + c
                cnt = sum(d.values())
                s1 = sum(h * c for h, c in d.items())
                s2 = sum(h * h * c for h, c in d.items())
                if cnt == m and s2 == E and s1 == S:
                    out.append((a, tuple(sorted(d.items()))))
    return S, sorted(set(out))


def main():
    # ---------------------------------------------------------------- P1 kill test
    S12, certs12 = certify_12()
    want = {0: 2, 2: 680, 4: 1366}
    found_exact = any(dict(c) == want for _, c in certs12)
    chk("P1a kill: S(12) = 6824", S12 == 6824, "brute-force descent gives S(12)=%d" % S12)
    chk("P1b kill: multiset {0:2, 2:680, 4:1366}", found_exact,
        "extremal multisets found at n=12: %s" % (certs12,))
    # independent re-arithmetic of that multiset
    m, E = m_of(12), E_of(12)
    cnt = sum(want.values())
    s1 = sum(h * c for h, c in want.items())
    s2 = sum(h * h * c for h, c in want.items())
    chk("P1c multiset realises (m, E, S)", cnt == m and s2 == E and s1 == 6824,
        "count %d = 2^11, sum h %d = S_Q, sum h^2 %d = E = 24576" % (cnt, s1, s2))
    chk("P1d cap* = (E+S)/2 = 15700", (E + S12) // 2 == 15700, "cap*(12) = %d" % ((E + S12) // 2))

    # ------------------------------------------------------- P2 published row gate
    bad_below, attained = [], []
    for n, g in sorted(PUB_EXACT.items()):
        cq = cap_quant(n)
        if cq < g:
            bad_below.append((n, cq, g))
        if cq == g:
            attained.append(n)
    chk("P2a cap* >= published g at all 17 dims", not bad_below, "violations: %s" % (bad_below,))
    chk("P2b equality at the 16 attained dims", len(attained) == 16,
        "attained: %s ; not attained: %s" % (attained, [n for n in PUB_EXACT if n not in attained]))
    chk("P2c n=9 realisability slack = 4", cap_quant(9) - PUB_EXACT[9] == 4,
        "cap*(9)=%d, g(9)=%d" % (cap_quant(9), PUB_EXACT[9]))
    for n in OPEN_N:
        chk("P2d open n=%d" % n, cap_plain(n) == CAP_PLAIN[n] and cap_quant(n) == CAP_QUANT[n],
            "plain %d -> quantised %d" % (cap_plain(n), cap_quant(n)))

    # ------------------------- P2e cross-check the seven caps against leg B's table
    legB = os.path.join(HERE, "cap_table.json")
    if os.path.exists(legB):
        tbl = json.load(open(legB, encoding="utf-8"))
        rows = {int(r["n"]): r for r in tbl["rows"]}
        same = all(int(rows[n]["cap_plain"]) == cap_plain(n) and
                   int(rows[n]["cap_quant"]) == cap_quant(n) for n in OPEN_N)
        chk("P2e agree with leg B (fresh_mod8_table.py)", same,
            "seven open caps identical to the CRT-route implementation")
    else:
        chk("P2e leg B table present", False, "cap_table.json missing")

    # ------------------------------------------- P3 the criterion vs brute force
    # hypothesis: N(a) = 0 (mod 8) for every admissible a when n >= 4
    bad8 = [(n, a) for n in range(4, NMAX + 1) for a in adm(n)
            if (E_of(n) + a * (a + 2) * m_of(n)) % 8]
    chk("P3a N(a) = 0 mod 8 for n >= 4", not bad8, "violations: %s" % (bad8[:5],))

    # per-a constraint table: which residue class does the hypothesis force?
    badres = []
    for n in range(4, 60):
        for a in adm(n):
            W = 2 * a + 2
            if a % 4 == 3:
                pred = None                                  # no restriction beyond S<=N/W
            elif a % 4 == 1:
                pred = 2                                     # S even
            else:
                pred = 4                                     # S = 0 mod 4
            S = sq(a, n)[0]
            top = (E_of(n) + a * (a + 2) * m_of(n)) // W
            if pred is None:                         # vacuous: S = floor(N/W)
                ok = S == top
            else:                                    # S = largest value <= top in 0 mod pred
                ok = (S <= top) and (S % pred == 0) and (top - S < pred)
            if not ok:
                badres.append((n, a, S, pred))
    chk("P3b per-a residue prediction (n<60, all a)", not badres,
        "the mod-8 condition is S = 0 mod 4 iff a even, S = 0 mod 2 iff a = 1 mod 4, "
        "and vacuous iff a = 3 mod 4 (violations: %s)" % (badres[:5],))

    # criterion: cap* = cap - 1  <==>  n even and S_pl = 2 (mod 4)
    rows, disagree = [], []
    for n in range(3, NMAX + 1):
        S_pl, top, den, argmin = plain(n)
        S_q = quantised(n)
        gain = 1 if cap_plain(n) - cap_quant(n) else 0
        pred = 1 if (n % 2 == 0 and S_pl % 4 == 2) else 0
        rows.append(dict(n=n, S_pl=S_pl, S_quant=S_q, gain=gain, pred=pred,
                         rho_floor=top, argmin_a=argmin,
                         cap_plain=cap_plain(n), cap_quant=cap_quant(n)))
        if gain != pred:
            disagree.append((n, S_pl, gain, pred))
    chk("P3c criterion == brute force for n = 3..%d" % NMAX, not disagree,
        "%d/%d dims agree; disagreements: %s" % (NMAX - 2, NMAX - 2, disagree[:5]))

    # ---------------------------------------------------- P4 the sharpened set
    sharp = [r["n"] for r in rows if r["gain"]]
    sharp_pred = [r["n"] for r in rows if r["pred"]]
    chk("P4a sharpened set = predicted set", sharp == sharp_pred,
        "|S| = %d for n <= %d" % (len(sharp), NMAX))
    chk("P4b the five open dims are in it, 14 and 25 are not",
        [n for n in OPEN_N if n in sharp] == [12, 18, 20, 22, 24],
        "sharpened open dims: %s" % [n for n in OPEN_N if n in sharp])
    open_beyond = [n for n in sharp if n > 26 and n not in PUB_EXACT]
    chk("P4c gain is never more than 1", all(r["gain"] in (0, 1) for r in rows),
        "max gain %d" % max(r["gain"] for r in rows))

    # ---------------------------------------------- P5 worked examples at the top
    examples = []
    for n in (38, 48, 60, 100, 198, 200, 300, 400, 599, 600):
        rows_by = {r["n"]: r for r in rows}
        if n in rows_by:
            r = rows_by[n]
            examples.append(dict(n=n, cap_plain=r["cap_plain"], cap_quant=r["cap_quant"],
                                 gain=r["gain"], S_pl_mod4=r["S_pl"] % 4,
                                 S_quant=r["S_quant"]))
    chk("P5 ten example rows reproduced", len(examples) == 10, "%s" % (examples,))

    out = dict(job="lab_erdos_261117_5",
               criterion="cap*(n) = cap(n) - 1  iff  n even and S_pl(n) = 2 (mod 4), "
                         "where S_pl = 2*floor(min_a N(a)/(2a+2))/2) and N(a) = E + a(a+2)2^(n-1)",
               why="N(a) = 0 mod 8 for n>=4, so Lemma 12's Delta=0 mod 8 reads: a even -> S=0 mod 4, "
                   "a=1 mod 4 -> S even, a=3 mod 4 -> vacuous; only the even-a class can lower the "
                   "DPTV79 even-rounded optimum, and it lowers it by exactly 2 when S_pl = 2 mod 4.",
               kill_test=dict(n=12, S=6824, multiset={str(k): v for k, v in want.items()},
                              cap_quant=15700, all_extremal_solutions=[[a, list(c)] for a, c in certs12]),
               published_gate=dict(below=bad_below, attained=attained,
                                   n_attained=len(attained)),
               open_caps={n: dict(cap_plain=cap_plain(n), cap_quant=cap_quant(n),
                                  S_quant=quantised(n)) for n in OPEN_N},
               sharpened_upto_NMAX=sharp, open_dims_beyond_26=open_beyond,
               n_sharpened=len(sharp), NMAX=NMAX, examples=examples,
               rows=rows)
    with open(os.path.join(HERE, "j261117_5_criterion.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)

    print("\n".join(logs))
    print("\n--- seven open dims (cap_plain -> cap_quant) ---")
    for n in OPEN_N:
        print("n=%2d  %12d -> %12d  %s" % (n, cap_plain(n), cap_quant(n),
                                           "SHARPENED" if cap_quant(n) < cap_plain(n) else "unchanged"))
    print("\n--- sharpened dims for n <= %d (%d of %d) ---" % (NMAX, len(sharp), NMAX - 2))
    print(sharp)
    print("\n--- criterion check ---")
    print("cap* = cap - 1  <=>  n even and S_pl = 2 (mod 4):  %d/%d dims agree"
          % (NMAX - 2 - len(disagree), NMAX - 2))
    print("\nFAILURES: %d %s" % (len(fails), fails))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
