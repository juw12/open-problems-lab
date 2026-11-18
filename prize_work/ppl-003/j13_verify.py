#!/usr/bin/env python3
# job lab_erdos_261112_13 -- INDEPENDENT VERIFIER (Method B, NOT the same code as j13_main):
#  (V1) N-stability: recompute |S|<=6 divergent fractions at N=5e6 (independent cache) vs N=1e7.
#  (V2) windowing/threshold robustness: |S|=5 fraction under 5 different (w_large,w_small,thresh)
#       -> is the ~0.020 plateau an artifact of the specific 0.05/0.00025/thresh-5 choice?
#  (V3) method independence: top-15 towers per |S| in {5,6,7}, recompute growth by DIRECT boolean
#       in-window count on the UNSORTED r (n-order) array -- a different code path than searchsorted
#       on r_sorted used by j13_main. Must agree to <1%.
import os, time, json, math
import numpy as np
from itertools import combinations
_T0 = time.time()

D = r'E:\Projects\MillenniumResearch\prize_work\ppl-003'
PRIMES_MAX = 199

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[0] = s[1] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if s[i]:
            s[i * i::i] = False
    return s

def load(N):
    phi = np.load(os.path.join(D, f'phi_cache_{N}.npy'))
    n = np.arange(N + 1, dtype=np.float64)
    r = phi[1:].astype(np.float64) / n[1:]
    return np.sort(r), r

def growth_ss(r_sorted, N, cvals, wL, wS):
    # Method A path: searchsorted on sorted
    g = np.empty(len(cvals))
    for st in range(0, len(cvals), 200000):
        en = min(st + 200000, len(cvals)); C = cvals[st:en]
        s0 = (np.searchsorted(r_sorted, C + wL, 'left') - np.searchsorted(r_sorted, C - wL, 'left')) / (2.0 * wL * N)
        s1 = (np.searchsorted(r_sorted, C + wS, 'left') - np.searchsorted(r_sorted, C - wS, 'left')) / (2.0 * wS * N)
        g[st:en] = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), np.nan)
    return g

def growth_brute(r_unsorted, N, c, wL, wS):
    # Method B path: direct boolean count on UNSORTED r
    s0 = float(np.count_nonzero((r_unsorted >= c - wL) & (r_unsorted <= c + wL))) / (2.0 * wL * N)
    s1 = float(np.count_nonzero((r_unsorted >= c - wS) & (r_unsorted <= c + wS))) / (2.0 * wS * N)
    return s1 / s0 if s0 > 0 else float('nan')

def enum_k(P, k):
    cvals = np.empty(math.comb(len(P), k), dtype=np.float64)
    f = {q: (q - 1) / q for q in P}; i = 0
    for S in combinations(P, k):
        c = 1.0
        for q in S: c *= f[q]
        cvals[i] = c; i += 1
    return cvals

P = [int(p) for p in primes_upto(PRIMES_MAX).nonzero()[0]]
out = {}

# ---------- V1: N-stability |S|<=6 at N=5e6 ----------
r5, r5u = load(5000000)
frac5 = {}
for k in range(1, 7):
    cv = enum_k(P, k)
    g = growth_ss(r5, 5000000, cv, 0.05, 0.00025)
    frac5[k] = round(float((g >= 5.0).mean()), 6)
# N=1e7 reference (recomputed here, same as main)
r10, r10u = load(10000000)
frac10 = {}
for k in range(1, 7):
    cv = enum_k(P, k)
    g = growth_ss(r10, 10000000, cv, 0.05, 0.00025)
    frac10[k] = round(float((g >= 5.0).mean()), 6)
out["V1_N_stability"] = {"frac_N5e6": frac5, "frac_N1e7": frac10,
                          "max_abs_diff_S4_6": round(max(abs(frac5[k]-frac10[k]) for k in (4,5,6)), 5)}
print("V1 N-stability |S|<=6:")
for k in (4,5,6):
    print(f"  |S|={k}: N5e6={frac5[k]} N1e7={frac10[k]} diff={abs(frac5[k]-frac10[k]):.5f}")

# ---------- V2: windowing/threshold robustness, |S|=5 ----------
cv5 = enum_k(P, 5)
combos = [ (0.05,0.00025,5.0,"lab-standard"), (0.03,0.0003,5.0,"w*0.6,thr5"),
           (0.02,0.0002,5.0,"w*0.4,thr5"), (0.05,0.00025,3.0,"thr3"),
           (0.05,0.00025,10.0,"thr10") ]
v2 = []
for wL,wS,th,lbl in combos:
    g = growth_ss(r10, 10000000, cv5, wL, wS)
    v2.append({"cfg": lbl, "wL": wL, "wS": wS, "thr": th, "frac_S5": round(float((g>=th).mean()),5)})
out["V2_windowing_S5"] = v2
print("\nV2 windowing/threshold robustness |S|=5 (frac should stay ~0.02):")
for v in v2: print(f"  {v['cfg']:14s} frac={v['frac_S5']}")

# ---------- V3: brute-force (unsorted) vs searchsorted, top-15 per |S| in {5,6,7} ----------
v3 = []
maxrel = 0.0
for k in (5,6,7):
    cv = enum_k(P, k)
    g_ss = growth_ss(r10, 10000000, cv, 0.05, 0.00025)
    order = np.argsort(-np.nan_to_num(g_ss, nan=-1))[:15]
    for i in order:
        c = float(cv[i]); g_b = growth_brute(r10u, 10000000, c, 0.05, 0.00025)
        rel = abs(g_b - g_ss[i]) / max(1e-9, abs(g_ss[i]))
        maxrel = max(maxrel, rel)
        v3.append({"absS": k, "c": round(c,6), "growth_searchsorted": round(float(g_ss[i]),4),
                   "growth_brute_unsorted": round(g_b,4), "rel_diff": round(rel,5)})
out["V3_bruteforce_top15"] = {"max_rel_diff": round(maxrel,5), "pass_lt_0.01": bool(maxrel < 0.01),
                               "detail": v3[:9]}
print(f"\nV3 brute-force (unsorted) vs searchsorted, top-15/|S| in {{5,6,7}}: max_rel_diff={maxrel:.5f} -> {'PASS' if maxrel<0.01 else 'FAIL'}")
for d in v3[:6]: print(f"  |S|={d['absS']} c={d['c']} ss={d['growth_searchsorted']} brute={d['growth_brute_unsorted']} rel={d['rel_diff']}")

json.dump(out, open(os.path.join(D, 'j13_verify.json'), 'w'), indent=2)
print(f"\nwrote j13_verify.json  ({time.time()-_T0:.0f}s)")
