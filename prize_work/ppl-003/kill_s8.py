#!/usr/bin/env python3
# kill_s8.py — cheap |S|=8 kill test + rate probe (job lab_erdos_261112_3).
# Runs ONLY k=8 for a few chunks at N=5e6, stops early. Confirms the |S|=8 class is
# nonzero (kill test) and measures the live scan rate to size the real run.
import sys, time
sys.path.insert(0, "E:/Projects/MillenniumResearch/prize_work/ppl-003")
from math import comb
import numpy as np
from scan_s8 import (SLUG, PRIMES_MAX, BLOCK, DIVERGE_THRESH, W_LARGE, W_SMALL,
                     primes_upto, combo_chunks, chunk_growth)

N = 5000000
MAX_CHUNKS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
t0 = time.time()
phi = np.load(f"{SLUG}/phi_cache_{N}.npy")
nn = np.arange(1, N + 1, dtype=np.float64)
r_sorted = np.sort(phi[1:].astype(np.float64) / nn)
ps = primes_upto(PRIMES_MAX)
P = [int(q) for q in range(2, PRIMES_MAX + 1) if ps[q]]
F = np.empty(PRIMES_MAX + 1, dtype=np.float64)
for q in P:
    F[q] = (q - 1) / q
t_sieve = time.time() - t0

found = 0
total_cand = 0
chunks = 0
t1 = time.time()
for combo in combo_chunks(P, 8, BLOCK):
    c = F[combo].prod(axis=1)
    g = chunk_growth(r_sorted, N, c)
    d = int((g >= DIVERGE_THRESH).sum())
    found += d
    total_cand += len(combo)
    chunks += 1
    if chunks >= MAX_CHUNKS:
        break
t_scan = time.time() - t1
rate = total_cand / t_scan if t_scan > 0 else 0
full8 = comb(46, 8)
full_all = sum(comb(46, k) for k in range(1, 9))
print("sieve+sort+setup: %.1fs" % t_sieve)
print("scanned %d |S|=8 cands in %d chunks: %.1fs -> %.0f cand/s" % (total_cand, chunks, t_scan, rate))
print("FOUND |S|=8 divergents (first %d cands): %d" % (total_cand, found))
print("KILL_TEST_nonzero =", found > 0)
print("frac divergent (first chunks) = %.5f" % (found / total_cand if total_cand else 0))
print("est |S|=8 full (%d) = %.0fs" % (full8, full8 / rate if rate else 0))
print("est full |S|<=8 (%d) = %.0fs" % (full_all, full_all / rate if rate else 0))
