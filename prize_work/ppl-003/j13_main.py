#!/usr/bin/env python3
# job lab_erdos_261112_13 -- MAIN: faithful |S|<=7 interior singular-set scan at N=1e7
# (my own code, reusing the cached int64 totient sieve). Produces the two universal invariants:
#   (1) the divergent-fraction plateau p*(|S|) for |S|=1..7
#   (2) the sharpness / mass structure (max growth per |S|, c=1/2 vs |S|>=5 towers, delta-f)
# GATE: |S|<=4 divergent count == 3008 (job 35) AND c=1/2 growth ~= 58.46.
import os, time, json, math
import numpy as np
from itertools import combinations

D = r'E:\Projects\MillenniumResearch\prize_work\ppl-003'
W_LARGE, W_SMALL, N_WIN = 0.05, 0.00025, 11
WINDOWS = np.geomspace(W_LARGE, W_SMALL, N_WIN)
DIVERGE_THRESH = 5.0
PRIMES_MAX, MAX_S = 199, 7

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[0] = s[1] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if s[i]:
            s[i * i::i] = False
    return s

def load_r_sorted(N):
    phi = np.load(os.path.join(D, f'phi_cache_{N}.npy'))
    n = np.arange(N + 1, dtype=np.float64)
    r = phi[1:].astype(np.float64) / n[1:]
    return np.sort(r), r

def scan_growth(r_sorted, N, cvals):
    w0, w1 = WINDOWS[0], WINDOWS[-1]
    growth = np.empty(len(cvals))
    for st in range(0, len(cvals), 200000):
        en = min(st + 200000, len(cvals)); C = cvals[st:en]
        s0 = (np.searchsorted(r_sorted, C + w0, 'left') - np.searchsorted(r_sorted, C - w0, 'left')) / (2.0 * w0 * N)
        s1 = (np.searchsorted(r_sorted, C + w1, 'left') - np.searchsorted(r_sorted, C - w1, 'left')) / (2.0 * w1 * N)
        growth[st:en] = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), np.nan)
    return growth

def df_change(r_sorted, N, c, w=0.02):
    # density change f(c+w)-f(c-w) = windowed mass
    return (np.searchsorted(r_sorted, c + w, 'left') - np.searchsorted(r_sorted, c - w, 'left')) / N

def enumerate(P, max_S):
    total = int(sum(math.comb(len(P), k) for k in range(1, max_S + 1)))
    cvals = np.empty(total, dtype=np.float64); absS = np.zeros(total, dtype=np.int8)
    has2 = np.zeros(total, dtype=bool); i = 0
    f = {q: (q - 1) / q for q in P}
    for k in range(1, max_S + 1):
        for S in combinations(P, k):
            c = 1.0
            for q in S: c *= f[q]
            cvals[i] = c; absS[i] = k; has2[i] = (S[0] == 2); i += 1
    return cvals, absS, has2

N = 10000000
t0 = time.time()
r_sorted, r = load_r_sorted(N)
P = [int(p) for p in primes_upto(PRIMES_MAX).nonzero()[0]]
t1 = time.time(); cvals, absS, has2 = enumerate(P, MAX_S); t_enum = time.time() - t1
t2 = time.time(); growth = scan_growth(r_sorted, N, cvals); t_scan = time.time() - t2
div = growth >= DIVERGE_THRESH
print(f"N={N} load={time.time()-t0:.1f}s enum={t_enum:.1f}s scan={t_scan:.1f}s total_div={int(div.sum())}")

out = {"N": N, "n_primes": len(P), "max_S": MAX_S, "t_enum_s": round(t_enum,1), "t_scan_s": round(t_scan,1),
       "windows": [W_LARGE, W_SMALL], "diverge_thresh": DIVERGE_THRESH,
       "total_divergent": int(div.sum())}

# GATE
m4 = absS <= 4
gate_count = int(div[m4].sum())
c12g = float(growth[np.argmin(np.abs(cvals - 0.5))])
out["gate"] = {"S4_divergent": gate_count, "S4_expect": 3008, "c12_growth": round(c12g,4), "c12_expect": 58.46,
               "pass": bool(gate_count == 3008 and abs(c12g - 58.46) < 1.0)}
print(f"GATE: |S|<=4 div={gate_count} (exp 3008) c12={c12g:.4f} (exp 58.46) -> {out['gate']['pass']}")

# INvariant 1: divergent fraction per |S|
frac = {}
mg = {}
cnt = {}
for k in range(1, MAX_S + 1):
    m = absS == k
    nd = int(div[m].sum())
    frac[k] = round(float(div[m].mean()), 6)
    cnt[k] = int(m.sum())
    mg[k] = round(float(np.nanmax(growth[m][div[m]])), 4) if nd else None
out["divergent_frac_by_absS"] = frac
out["max_growth_by_absS"] = mg
out["n_towers_by_absS"] = cnt
print("divergent_frac_by_absS:", frac)
print("max_growth_by_absS   :", mg)

# plateau: mean of |S|=4..7 and spread
plateau = [frac[k] for k in (4,5,6,7)]
out["plateau_S4_7"] = {"mean": round(float(np.mean(plateau)),5), "min": round(min(plateau),5),
                        "max": round(max(plateau),5), "rel_spread_pct": round(100*(max(plateau)-min(plateau))/np.mean(plateau),2)}
plateau57 = [frac[k] for k in (5,6,7)]
out["plateau_S5_7"] = {"mean": round(float(np.mean(plateau57)),5), "rel_spread_pct": round(100*(max(plateau57)-min(plateau57))/np.mean(plateau57),2)}

# INvariant 2: c=1/2 vs top |S|>=5 (mass vs sharpness)
top = {}
for k in (1,2,3,4,5,6,7):
    m = (absS == k) & div
    idx = np.where(m)[0]
    if len(idx)==0: top[k]=[]; continue
    idx = idx[np.argsort(-growth[idx])][:3]
    lst=[]
    for i in idx:
        lst.append({"c": round(float(cvals[i]),6), "growth": round(float(growth[i]),4),
                    "df_0p02": round(float(df_change(r_sorted,N,float(cvals[i]),0.02)),5),
                    "has2": bool(has2[i])})
    top[k]=lst
out["top3_by_growth_per_absS"] = top
print("top3 per |S| (c, growth, df_0.02, has2):")
for k in sorted(top):
    print(f"  |S|={k}:", top[k])

# c=1/2 detail
i12 = int(np.argmin(np.abs(cvals - 0.5)))
out["c12_detail"] = {"c": 0.5, "growth": round(float(growth[i12]),4),
                     "df_0p02": round(float(df_change(r_sorted,N,0.5,0.02)),5),
                     "max_growth_S1": mg[1], "note": "c=1/2 vs its own |S|=1 class max and vs |S|>=4 ceiling"}
print(f"\nc=1/2: growth={out['c12_detail']['growth']} df_0.02={out['c12_detail']['df_0p02']} | max |S|=1 growth={mg[1]} | |S|>=4 ceiling~{mg[4]}")

json.dump(out, open(os.path.join(D, 'j13_invariants.json'), 'w'), indent=2)
print(f"\nwrote j13_invariants.json  ({time.time()-t0:.1f}s total)")
