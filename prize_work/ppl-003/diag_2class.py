#!/usr/bin/env python3
# diag_2class.py — is the |S|=8 2-class collapse (2pqrstuvw=0) real or a threshold artifact?
# Scans ONLY the 2-class |S|=8 (combos [2]+7 odd primes, C(45,7)=45,379,620) at N=5e6 and
# reports the growth distribution (max + fraction above 1/2/3/5/10). If max growth < 5.0,
# the collapse is REAL (no 2-class point reaches the divergent threshold); if some exceed
# 5.0, the main run's 2-class count has a bug.
import sys, time, itertools
sys.path.insert(0, "E:/Projects/MillenniumResearch/prize_work/ppl-003")
import numpy as np
from scan_s8 import SLUG, PRIMES_MAX, BLOCK, W_LARGE, W_SMALL, DIVERGE_THRESH, primes_upto, chunk_growth

N = 5000000
phi = np.load(f"{SLUG}/phi_cache_{N}.npy")
nn = np.arange(1, N + 1, dtype=np.float64)
r_sorted = np.sort(phi[1:].astype(np.float64) / nn)
ps = primes_upto(PRIMES_MAX)
P = [int(q) for q in range(2, PRIMES_MAX + 1) if ps[q]]
odd = [q for q in P if q != 2]
F = np.empty(PRIMES_MAX + 1, dtype=np.float64)
for q in P:
    F[q] = (q - 1) / q

def c2_chunks(block):
    it = itertools.combinations(odd, 7)
    while True:
        buf = list(itertools.islice(it, block))
        if not buf:
            break
        a = np.array(buf, dtype=np.int32)
        yield F[a].prod(axis=1) * F[2]     # c for S = [2] + 7 odd primes

t0 = time.time()
gmax = 0.0
above = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 8: 0, 10: 0}
total = 0
for c in c2_chunks(BLOCK):
    g = chunk_growth(r_sorted, N, c)
    total += len(c)
    gm = float(np.nanmax(g))
    if gm > gmax:
        gmax = gm
    for t in above:
        above[t] += int(np.nansum(g >= t))
t = time.time() - t0
print("2class |S|=8 combos scanned: %d  (%.1fs)" % (total, t))
print("MAX growth (2class |S|=8) = %.4f   (DIVERGE_THRESH=%.1f)" % (gmax, DIVERGE_THRESH))
print("fraction above threshold: " + ", ".join(">= %d: %.5f" % (k, above[k] / total) for k in sorted(above)))
print("VERDICT: collapse real (max<5) =" , gmax < DIVERGE_THRESH)
