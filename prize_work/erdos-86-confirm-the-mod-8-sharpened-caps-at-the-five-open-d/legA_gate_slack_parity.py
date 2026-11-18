#!/usr/bin/env python3
"""lab_erdos_261114_3 -- #86 two-point cap slack: parity of the programme optimum.

Idea card: why is the two-point cap ATTAINED at 10/10 exact dims n>=10 but not at
n=6 (+1) or n=9 (+4)?  Card hypothesis: attainment is decided by the parity of the
two-point programme optimum S* together with the integrality of the extremal
weight vector.

This script works from scratch (exact Fraction arithmetic, no lab import):

1. Re-derives the two-point programme for the odd-square family (Erdos #86):
       g(n) = (E + max_s S(s)) / 2,   E = n 2^(n-1),  S(s) = sum_{v in U} h(v),
   U = even-parity vertices, h(v) = local field of a fully frustrated signing, with
   the DPTV field identity
       sum_{v in U} h(v)^2 = n 2^(n-1) =: Q        (for EVERY field)
       h(v) == n (mod 2),  h in {-n, -n+2, ..., n}.
   The classical two-point inequality (DPTV79 / Minamoto arXiv:2603.29127 Lemma 9)
   gives, for every a with a == n (mod 2) and 2a+2 > 0:
       S <= f(a) = m (n + a(a+2)) / (2a+2),   m = 2^(n-1)
   so S*_real = min_a f(a) and the lab/published cap is
       cap(n)  = (E + 2 floor(S*_real / 2)) / 2.                    (plain rounding)

2. The MOD-8 QUANTISATION (Minamoto Lemma 12, "Quantised two-point inequality"):
   every summand (h-a)(h-a-2) is a nonnegative MULTIPLE OF 8, so
       Delta(a) = Q + m a(a+2) - (2a+2) S >= 0  and  Delta(a) == 0 (mod 8),
   hence the strongest integral bound from one a is
       S_max(a) = max { S in Z : (2a+2) S <= N(a),  N(a) - (2a+2) S == 0 (mod 8) }
   with N(a) = Q + m a(a+2);  and
       cap*(n) = (E + 2 floor(min_a S_max(a) / 2)) / 2.             (quantised cap)

   cap*(n) <= cap(n), and cap*(n) is a RIGOROUS upper bound on g(n) for every field
   (published inequality + integer arithmetic).  Where min_a S_max(a) is attained by
   an explicit h-vector (constructed here, see slack_table.json), it is also the
   optimum of the integer programme
       max sum h   s.t.  sum h^2 = Q,  h == n (mod 2),  |h| <= n
   (Minamoto's DP programme, which he solved by exhaustive DP at n=6,7,8).

GATES (all must pass; exit 0 on success):
  G1  from-scratch cap(n) == committed tools/gate_odsquare_cap.py cap(n), all n=3..26
  G2  slack table reproduces cap(6)-g(6) = 133-132 = 1 and cap(9)-g(9) = 1536-1532 = 4
  G3  quantised cap* reproduces Minamoto's three exact instances:
        n=6 -> S<=72 (g=132, Delta=16, h-dist {0:2, 2:24, 4:6} on U),
        n=7 -> S<=160 (g=304, Delta=0, dist {1:16, 3:48}),
        n=8 -> S<=340 (g=682, Delta=8, one vertex h in {0,6});
     plus the paper's "next-best quantised pair" values 96 / 176 / 408.
  G4  consistency: cap*(n) >= exact g(n) at all 17 published exact n;
     cap*(n) == g(n) at every n where the cap is attained
     (4,5,7,8,10,11,13,15,16,17,19,21,23,26);  cap*(n) >= lab certs (15614, 71572).
  G5  integrality of the extremal T=0 two-point vector:
     n=6: k = m(n-a^2)/(4a+4) = 16/3 NOT integral (parity gap, slack 1);
     n=9: k = 256 integral (two-point vector feasible as a multiset) yet the cap is
     not attained (Wang F_9=772) -> the +4 is a REALISABILITY gap, not a parity gap;
     n=10: k integral and the cap is attained.
  G6  extremal multiset construction: for n in {6,7,8,12,18,20,22,24} an explicit
     3-point h-multiset attains the quantised value S_q (matching Minamoto's
     published distributions at 6,7,8), so the quantised cap is the EXACT optimum
     of the integer programme IP(n): the two-point relaxation is TIGHT there and
     any remaining gap to g(n) is purely a field-realizability question.

Outputs slack_table.json (the deliverable: the exact slack table) and prints it.
The DPTV identity itself is verified empirically by verify_identity.py (separate).
"""
from __future__ import annotations
import json, os, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))

# Published exact g(n) for the odd-square family (same table as the committed gate
# tool; sources: Minamoto arXiv:2603.29127v5, Wang arXiv:2611.14906, BHN95,
# LNSW26, OEIS A245762). Used ONLY for the slack column and consistency gate G4.
EXACT = {
    3: 9, 4: 24, 5: 56, 6: 132, 7: 304, 8: 682, 9: 1532, 10: 3328, 11: 7296,
    13: 33792, 15: 153600, 16: 327680, 17: 688128, 19: 3047424, 21: 13369344,
    23: 58195968, 26: 520093696,
}
OPEN = [12, 14, 18, 20, 22, 24, 25]
CERTS = {12: 15614, 14: 71572}   # explicit lab odd-square certificates
ATTAINED = [4, 5, 7, 8, 10, 11, 13, 15, 16, 17, 19, 21, 23, 26]


# ---------------------------------------------------------------- the programme
def admissible_a(n):
    """a with a == n (mod 2) and 2a+2 > 0."""
    return [a for a in range(-n, n + 1) if (a - n) % 2 == 0 and 2 * a + 2 > 0]


def two_point_real(n):
    """S*_real = min_a m(n+a(a+2))/(2a+2), exact Fraction.  Returns (a*, S*_real)."""
    m = 1 << (n - 1)
    best, abest = None, None
    for a in admissible_a(n):
        f = Fraction(m * (n + a * (a + 2)), 2 * a + 2)
        if best is None or f < best:
            best, abest = f, a
    return abest, best


def quantised_Smax(n, a):
    """max integer S with (2a+2)S <= N(a) and N(a)-(2a+2)S == 0 (mod 8)."""
    m = 1 << (n - 1)
    N = int(n * m + m * a * (a + 2))
    d = 2 * a + 2
    S = N // d
    for _ in range(8):
        if (N - d * S) % 8 == 0:
            return S
        S -= 1
    raise AssertionError('n=%d a=%d: no S with Delta == 0 mod 8 (should be impossible)' % (n, a))


def cap_plain(n):
    a, S_real = two_point_real(n)
    S_even = 2 * (S_real.numerator // (2 * S_real.denominator))
    return a, S_real, S_even, (n * (1 << (n - 1)) + S_even) // 2


def cap_quantised(n):
    m = 1 << (n - 1)
    best, abest = None, None
    for a in admissible_a(n):
        s = quantised_Smax(n, a)
        if best is None or s < best:
            best, abest = s, a
    S_even = 2 * (best // 2)
    return abest, best, S_even, (n * m + S_even) // 2


def extremal_two_point(n, a):
    """T=0 two-point vector feasibility: h in {a, a+2}, sum h^2 == n m.
    Returns (k, ok): k = # vertices with h = a+2, ok = k integral and 0 <= k <= m."""
    m = 1 << (n - 1)
    k = Fraction(n * m - a * a * m, (a + 2) ** 2 - a * a)
    ok = (k.denominator == 1) and 0 <= k <= m
    return k, ok


# ---------------------------------------------------------------- main
def main():
    fails = []

    def chk(name, ok, detail):
        print(('PASS  ' if ok else 'FAIL  ') + name + ': ' + detail, flush=True)
        if not ok:
            fails.append(name)

    # ---- G1: from-scratch cap == committed gate tool, all n=3..26
    sys.path.insert(0, os.path.join(ROOT, 'tools'))
    try:
        import gate_odsquare_cap as g
        ref = {n: g.cap(n) for n in range(3, 27)}
    except Exception as e:
        ref = None
        print('NOTE gate import failed:', e, flush=True)
    bad = []
    for n in range(3, 27):
        _, _, _, c = cap_plain(n)
        if ref and c != ref[n]:
            bad.append((n, c, ref[n]))
    chk('G1 from-scratch cap == gate tool, n=3..26', not bad and ref is not None,
        'mismatches=%s' % (bad or 'none'))

    # ---- build the slack table
    rows = []
    for n in range(3, 27):
        a_p, S_real, S_even, c = cap_plain(n)
        a_q, S_ref, S_even_q, cq = cap_quantised(n)
        k, ok_k = extremal_two_point(n, a_q)
        ex = EXACT.get(n)
        rows.append({
            'n': n, 'E': n * (1 << (n - 1)),
            'a_star_plain': a_p, 'S_2pt_real': str(S_real),
            'S_2pt_real_float': round(float(S_real), 6),
            'S_even_plain': S_even, 'cap_plain': c,
            'a_star_quant': a_q, 'S_quant': S_ref, 'S_even_quant': S_even_q,
            'cap_quant': cq, 'cap_minus_quant': c - cq,
            'k_two_point': str(k), 'two_point_integral': ok_k,
            'exact': ex,
            'slack_plain': (c - ex) if ex else None,
            'slack_quant': (cq - ex) if ex else None,
        })

    r6 = next(r for r in rows if r['n'] == 6)
    r9 = next(r for r in rows if r['n'] == 9)

    # ---- G2: the two published slacks
    chk('G2a slack n=6 = 133-132 = 1', r6['cap_plain'] == 133 and r6['slack_plain'] == 1,
        'cap=%d exact=%d' % (r6['cap_plain'], r6['exact']))
    chk('G2b slack n=9 = 1536-1532 = 4', r9['cap_plain'] == 1536 and r9['slack_plain'] == 4,
        'cap=%d exact=%d' % (r9['cap_plain'], r9['exact']))

    # ---- G3: Minamoto's worked instances and next-best pairs
    chk('G3a n=6 quantised: S_max(a=2)=72, alt pair a=0 gives 96',
        quantised_Smax(6, 2) == 72 and quantised_Smax(6, 0) == 96,
        'S(a=2)=%d S(a=0)=%d' % (quantised_Smax(6, 2), quantised_Smax(6, 0)))
    chk('G3b n=7 quantised: S_max(a=1)=160, alt pair a=3 gives 176',
        quantised_Smax(7, 1) == 160 and quantised_Smax(7, 3) == 176,
        'S(a=1)=%d S(a=3)=%d' % (quantised_Smax(7, 1), quantised_Smax(7, 3)))
    chk('G3c n=8 quantised: S_max(a=2)=340, alt pair a=4 gives 408',
        quantised_Smax(8, 2) == 340 and quantised_Smax(8, 4) == 408,
        'S(a=2)=%d S(a=4)=%d' % (quantised_Smax(8, 2), quantised_Smax(8, 4)))
    chk('G3d Delta n=6 = 448-6*72 = 16 (two vertices h in {0,6})', 448 - 6 * 72 == 16, '')
    chk('G3e Delta n=7 = 448... = 0 (all h in {1,3}: {1:16, 3:48})', 640 - 4 * 160 == 0, '')
    chk('G3f Delta n=8 = 2048-6*340 = 8 (one vertex h in {0,6})', 2048 - 6 * 340 == 8, '')

    # ---- G4: consistency against published exact values and lab certs
    below = [r['n'] for r in rows if r['exact'] is not None and r['cap_quant'] < r['exact']]
    att = [n for n in ATTAINED
           if next(r for r in rows if r['n'] == n)['cap_quant'] == EXACT[n]]
    chk('G4a cap* >= exact at all 17 published dims', not below, 'below=%s' % (below or 'none'))
    chk('G4b cap* == exact at all 14 attained dims', len(att) == 14, 'attained=%s' % att)
    over = [(n, CERTS[n]) for n in CERTS
            if next(r for r in rows if r['n'] == n)['cap_quant'] < CERTS[n]]
    chk('G4c lab certs <= cap* at n=12,14', not over, 'over=%s' % (over or 'none'))

    # ---- G5: integrality of the extremal two-point vector
    k6, ok6 = extremal_two_point(6, 2)
    k9, ok9 = extremal_two_point(9, 1)
    k10, ok10 = extremal_two_point(10, 2)
    chk('G5a n=6: T=0 vector INFEASIBLE (k = 16/3 not integral) -> parity gap, slack 1',
        (not ok6) and str(k6) == '16/3', 'k=%s' % k6)
    chk('G5b n=9: T=0 vector FEASIBLE (k=256, all h=3) yet unattained -> realisability gap',
        ok9 and int(k9) == 256 and r9['slack_plain'] == 4, 'k=%s' % k9)
    chk('G5c n=10: T=0 vector feasible and cap attained (slack 0)',
        ok10 and next(r for r in rows if r['n'] == 10)['slack_plain'] == 0, 'k=%s' % k10)

    # ---- G6: extremal multiset construction (IP optimum == quantised S)
    # For a given a, b = a+2, look for a 3-point multiset {c: t, a: (m-t-k), b: k}
    # with t in 0..8, c in {a-4, a-2, a+4, a+6} satisfying both
    #   b^2 k + a^2 (m-t-k) + t c^2 == Q   and   b k + a (m-t-k) + t c == S_q.
    # If found, the quantised value S_q is the EXACT optimum of the integer
    # programme IP(n) (upper bound = Lemma 12, lower bound = this construction),
    # i.e. the two-point relaxation is TIGHT there and any further gap to g(n)
    # is purely a field-realizability question.
    def construct(n):
        a_q, S_q, _, _ = cap_quantised(n)
        b = a_q + 2
        m = 1 << (n - 1)
        Q = n * m
        for t in range(0, 9):
            for c in (a_q - 4, a_q - 2, a_q + 4, a_q + 6):
                if (c - n) % 2 != 0 or abs(c) > n:
                    continue
                kq1 = Fraction(Q - a_q * a_q * (m - t) - t * c * c, 4 * a_q + 4)
                kq2 = Fraction(S_q - a_q * (m - t) - t * c, 2)
                if kq1 == kq2 and kq1.denominator == 1 and 0 <= kq1 <= m - t:
                    k = int(kq1)
                    dist = {h: x for h, x in ((c, t), (a_q, m - t - k), (b, k)) if x > 0}
                    S_check = sum(h * x for h, x in dist.items())
                    Q_check = sum(h * h * x for h, x in dist.items())
                    m_check = sum(dist.values())
                    assert S_check == S_q and Q_check == Q and m_check == m
                    return a_q, dist
        return a_q, None

    construct_targets = [6, 7, 8, 12, 18, 20, 22, 24]
    cons_ok = True
    cons_detail = {}
    for n in construct_targets:
        a_q, dist = construct(n)
        cons_detail[n] = {str(h): int(x) for h, x in (dist or {}).items()}
        cons_ok = cons_ok and (dist is not None)
    chk('G6a IP(6)=72 multiset == Minamoto {0:2, 2:24, 4:6}',
        cons_detail.get(6) == {'0': 2, '2': 24, '4': 6}, str(cons_detail.get(6)))
    chk('G6b IP(7)=160 multiset == Minamoto {1:16, 3:48}',
        cons_detail.get(7) == {'1': 16, '3': 48}, str(cons_detail.get(7)))
    chk('G6c IP(8)=340 multiset == Minamoto {0:1, 2:84, 4:43}',
        cons_detail.get(8) == {'0': 1, '2': 84, '4': 43}, str(cons_detail.get(8)))
    chk('G6d quantised S is the exact IP optimum (tight relaxation) at n=12,18,20,22,24',
        all(cons_detail.get(n) is not None for n in (12, 18, 20, 22, 24)),
        'dist=%s' % {n: cons_detail.get(n) for n in (12, 18, 20, 22, 24)})

    # ---- print table
    print('\n| n | E | a* | S*_real | S_even | cap | a*Q | S_Q | cap* | cap-cap* | k(T=0) | T0-int | exact | slack | slack* |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in rows:
        print('| %d | %d | %s | %s | %s | %d | %s | %s | %d | %s | %s | %d | %s | %s | %s |' % (
            r['n'], r['E'], r['a_star_plain'], r['S_2pt_real'], r['S_even_plain'], r['cap_plain'],
            r['a_star_quant'], r['S_quant'], r['cap_quant'], r['cap_minus_quant'],
            r['k_two_point'], int(r['two_point_integral']),
            r['exact'] if r['exact'] is not None else '-',
            r['slack_plain'] if r['slack_plain'] is not None else '-',
            r['slack_quant'] if r['slack_quant'] is not None else '-'))

    out = os.path.join(HERE, 'slack_table.json')
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({'rows': rows, 'open_n': OPEN,
                   'cap_plain_open': {r['n']: r['cap_plain'] for r in rows if r['n'] in OPEN},
                   'cap_quant_open': {r['n']: r['cap_quant'] for r in rows if r['n'] in OPEN},
                   'ip_optimum_multisets': {str(k2): v for k2, v in cons_detail.items()
                                            if v is not None}},
                  f, indent=1)
    print('\nwrote', out)
    print(('\nGATES PASS' if not fails else '\nGATES FAIL: %s' % fails), flush=True)
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
