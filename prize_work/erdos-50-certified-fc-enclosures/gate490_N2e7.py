"""lab_erdos_261113_43 -- gate490 at a SECOND N = 2e7.

CARD GATE:  the integer atom count must reproduce job 34's 490/490 (point,window)
            exact agreement at a second N.
CARD VERIFIER (two independent implementations of the atom sum agree on the integer
            count at three (point,window) pairs):
  (1) pruned DFS on squarefree prime sets   = job 34's atoms_band_arrays (int64 M_U(N))
  (2) totient-sieve direct count            = EXACT integer comparison cd*phi(n) < cn*n
                                              and cd*b*phi(n) > L*n (no floating point)
Extras:
  * monotone in N: model(2e7) >= model(1e7) and measured(2e7) >= measured(1e7) at all 490
  * G4: the 99-grid at N=1e5 by THREE independent lineages:
        naive trial division / int64 phi cache (ppl-003) / int64 Euler sieve (this job)
"""
import sys, json, time
import numpy as np
from fractions import Fraction
from math import gcd

sys.path.insert(0, r"E:\Projects\MillenniumResearch\prize_work\erdos-confirm-the-27-239-s-5-points-are-true-singularities-v")
from atomgroups import atoms_band_arrays
from atomlib import primes_upto

ROOT = r"E:\Projects\MillenniumResearch"
J34 = ROOT + r"\prize_work\erdos-confirm-the-27-239-s-5-points-are-true-singularities-v"
OUT = ROOT + r"\prize_work\erdos-50-certified-fc-enclosures"
SW = ROOT + r"\prize_work\ppl-003"
N = 20_000_000
W0 = 5e-3
LADDER = [5e-3, 2e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
t0 = time.time()

job34 = json.load(open(J34 + r"\stage_A_N10000000.json"))
assert job34["N"] == 10_000_000 and len(job34["primary"]) == 50

# ---------------------------------------------------------------- measured side
phi = np.load(f"{SW}/phi_cache_{N}.npy")
rs = np.sort(phi[1:].astype(np.float64) / np.arange(1, N + 1, dtype=np.float64))
print(f"[measured setup] {time.time()-t0:.1f}s", flush=True)

per_group = {}
total_pairs = total_agree = 0
mono_model_bad = mono_meas_bad = 0
primary_mod = {}   # (idx, w) -> model int, for the exact verifier
primary_pts = None

for key, label in [("primary", "primary"), ("secondary_growth_2_3", "sec23"),
                   ("secondary_off_half", "secoff")]:
    rows34 = job34[key]
    pts = []
    for r in rows34:
        cn, cd = (int(x) for x in r["c"].split("/"))
        pts.append((cn, cd, tuple(r["S"])))
    lots = [cn / cd for cn, cd, _ in pts]
    c_lo = min(lots) - W0
    c_hi = max(lots)
    t = time.time()
    cfl, cnt, prd, kk, nodes = atoms_band_arrays(c_lo, c_hi, N, verbose=True)
    o = np.argsort(cfl)
    cfl_s, cnt_s = cfl[o], cnt[o]
    csum = np.concatenate(([0], np.cumsum(cnt_s)))
    g_pairs = g_agree = 0
    for i, ((cn, cd, S), lot) in enumerate(zip(pts, lots)):
        r34 = rows34[i]
        i_hi = int(np.searchsorted(rs, lot, 'left'))          # measured side (rs, size N)
        hi_c = int(np.searchsorted(cfl_s, lot, 'left'))       # model side (atoms, size ~1e6)
        meas, mod = {}, {}
        for w in LADDER:
            meas[w] = int(i_hi - np.searchsorted(rs, lot - w, 'left'))
            mod[w] = int(csum[hi_c] - csum[int(np.searchsorted(cfl_s, lot - w, 'left'))])
            if key == "primary":
                primary_mod[(i, w)] = mod[w]
        g_pairs += 7
        g_agree += sum(1 for w in LADDER if meas[w] == mod[w])
        mono_model_bad += sum(1 for w in LADDER if mod[w] < int(r34["model"][str(w)]))
        mono_meas_bad += sum(1 for w in LADDER if meas[w] < int(r34["measured"][str(w)]))
    total_pairs += g_pairs
    total_agree += g_agree
    if label == "primary":
        primary_pts = pts
    per_group[label] = {"pairs": g_pairs, "agree": g_agree,
                        "band": [c_lo, c_hi], "atoms": int(cnt.size), "nodes": int(nodes),
                        "s": round(time.time() - t, 1)}
    print(f"[{label}] {g_agree}/{g_pairs} agree  ({time.time()-t:.1f}s)", flush=True)

gate490 = total_agree == total_pairs
print(f"[GATE490] {total_agree}/{total_pairs} at N=2e7 -> {gate490};  "
      f"mono model violations={mono_model_bad}, meas violations={mono_meas_bad}", flush=True)

# ---------------------------------------------------------------- exact verifier
# three (point, window) pairs, windows as exact decimal fractions; int64-safe choices
VERIF = [(0, Fraction(1, 200)), (24, Fraction(1, 10000)), (49, Fraction(2, 1000))]
x = phi[1:].astype(np.int64)
narr = np.arange(1, N + 1, dtype=np.int64)
verif_rows = []
for idx, wf in VERIF:
    cn, cd, S = primary_pts[idx]
    a, b = wf.numerator, wf.denominator
    L = cn * b - cd * a                      # lower bound = L/(cd*b), exact rational
    assert L > 0
    # int64 overflow safety (max product <= 4e18 << 9.22e18)
    while cd * b * N > 4e18 or L * N > 4e18:
        wf = wf / 10
        a, b = wf.numerator, wf.denominator
        L = cn * b - cd * a
    below_hi = cd * x < cn * narr            # phi(n)/n < c       (strict)
    above_lo = cd * b * x > L * narr         # phi(n)/n > c - w   (strict)
    exact = int(np.count_nonzero(below_hi & above_lo))
    w = float(a / b)
    model = primary_mod[(idx, w)]
    verif_rows.append({"point_idx": idx, "S": list(S), "c": f"{cn}/{cd}",
                       "w": str(a) + "/" + str(b), "model_dfs": int(model),
                       "exact_int": exact, "agree": exact == model})
    print(f"[VERIF] S={list(S)} c={cn}/{cd} w={a}/{b}: dfs={model} exact={exact} "
          f"agree={exact == model}", flush=True)
verif_pass = all(v["agree"] for v in verif_rows)

# ---------------------------------------------------------------- G4: N=1e5 x 3 lineages
N5 = 100_000
t = time.time()
small_p = [int(p) for p in primes_upto(400)]
def phi_naive(n):
    r, m = n, n
    for d in small_p:
        if d * d > m:
            break
        if m % d == 0:
            while m % d == 0:
                m //= d
            r -= r // d
    if m > 1:
        r -= r // m
    return r
ph_n = np.array([phi_naive(i) for i in range(1, N5 + 1)], dtype=np.int64)
ph_c = phi[1:N5 + 1].astype(np.int64)                       # int64 cache lineage
ph_s = np.arange(N5 + 1, dtype=np.int64)                     # Euler sieve lineage
for i in range(2, N5 + 1):
    if ph_s[i] == i:
        ph_s[i::i] -= ph_s[i::i] // i
n5 = np.arange(1, N5 + 1, dtype=np.int64)
def grid99(ph):
    q = (100 * ph) // n5
    bc = np.bincount(q, minlength=101)
    return bc.cumsum()[:99].astype(np.int64)
g_n, g_c, g_s = grid99(ph_n), grid99(ph_c), grid99(ph_s[1:N5 + 1])
g4 = bool(np.array_equal(g_n, g_c) and np.array_equal(g_c, g_s))
print(f"[G4] 99-grid at N=1e5: naive vs cache vs sieve -> {g4} ({time.time()-t:.1f}s)", flush=True)

res = {"N": N, "ladder": [str(w) for w in LADDER],
       "gate490": {"pairs": total_pairs, "agree": total_agree, "pass": bool(gate490)},
       "per_group": per_group,
       "monotone_in_N": {"model_violations": int(mono_model_bad),
                         "meas_violations": int(mono_meas_bad)},
       "verifier_exact_vs_dfs": verif_rows,
       "verifier_pass": bool(verif_pass),
       "G4_N1e5_three_lineages": {"pass": g4},
       "elapsed_s": round(time.time() - t0, 1)}
json.dump(res, open(OUT + r"\gate490_N20000000.json", "w"), indent=1)
print(f"[done] {time.time()-t0:.0f}s", flush=True)
