#!/usr/bin/env python3
# job lab_erdos_261112_13 -- KILL-TEST / GATE: confirm a faithful re-implementation of the
# windowed-slope growth machinery reproduces the lab's known |S|<=4 values at N=1e7
# (3008 divergent, c=1/2 growth 58.46). If this gate fails, the method is untrusted -> switch.
import os, time, json, math
import numpy as np
from itertools import combinations

D = r'E:\Projects\MillenniumResearch\prize_work\ppl-003'
W_LARGE, W_SMALL, N_WIN = 0.05, 0.00025, 11
WINDOWS = np.geomspace(W_LARGE, W_SMALL, N_WIN)
DIVERGE_THRESH = 5.0
PRIMES_MAX, MAX_S = 199, 4

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[0] = s[1] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if s[i]:
            s[i * i::i] = False
    return s

def load_r_sorted(N):
    phi = np.load(os.path.join(D, f'phi_cache_{N}.npy'))   # int64 totients phi[0..N]
    n = np.arange(N + 1, dtype=np.float64)
    r = phi[1:].astype(np.float64) / n[1:]                 # phi(n)/n for n=1..N
    r_sorted = np.sort(r)
    return r_sorted, r

def scan_growth(r_sorted, N, cvals):
    w0, w1 = WINDOWS[0], WINDOWS[-1]
    growth = np.empty(len(cvals))
    for st in range(0, len(cvals), 200000):
        en = min(st + 200000, len(cvals)); C = cvals[st:en]
        s0 = (np.searchsorted(r_sorted, C + w0, 'left') - np.searchsorted(r_sorted, C - w0, 'left')) / (2.0 * w0 * N)
        s1 = (np.searchsorted(r_sorted, C + w1, 'left') - np.searchsorted(r_sorted, C - w1, 'left')) / (2.0 * w1 * N)
        growth[st:en] = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), np.nan)
    return growth

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
mean_r = float(r.mean())
print(f"loaded r_sorted N={N} in {time.time()-t0:.1f}s; mean(r)={mean_r:.9f} vs 6/pi^2={6/math.pi**2:.9f}")
assert abs(mean_r - 6/math.pi**2) < 2e-5, "phi-cache sanity gate FAILED"

P = [int(p) for p in primes_upto(PRIMES_MAX).nonzero()[0]]  # [2,3,5,...,199]
print(f"n primes <=199: {len(P)} (first={P[0]}, last={P[-1]})")
t1 = time.time()
cvals, absS, has2 = enumerate(P, MAX_S)
print(f"enumerated |S|<={MAX_S}: {len(cvals)} towers in {time.time()-t1:.1f}s")
t2 = time.time()
growth = scan_growth(r_sorted, N, cvals)
div = growth >= DIVERGE_THRESH
print(f"scan done in {time.time()-t2:.1f}s; total divergent = {int(div.sum())}")

# per-|S|
for k in range(1, MAX_S + 1):
    m = absS == k
    frac = div[m].mean() if m.sum() else 0.0
    print(f"  |S|={k}: n_div={int(div[m].sum()):6d} / {int(m.sum()):7d}  frac={frac:.6f}  max_growth={growth[m][div[m]].max():.4f}" if div[m].any() else f"  |S|={k}: n_div=0")

# c=1/2 gate
c12 = (1/2)
idx = np.argmin(np.abs(cvals - c12))
print(f"\nc=1/2 tower S=[2]: growth={growth[idx]:.4f} (expect ~58.46)")
print(f"total divergent |S|<=4 = {int(div.sum())} (expect 3008)")
print(f"ALL GATES: count_ok={int(div.sum())==3008}  c12_ok={abs(growth[idx]-58.46)<1.0}")
