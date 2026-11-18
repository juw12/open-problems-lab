#!/usr/bin/env python3
"""Smoke test: all job-35 code paths at N=20k (seconds). Not a result."""
import sys, time, math
HERE = r"E:\Projects\MillenniumResearch\prize_work\erdos-143-extend-the-d0-all-g-pin-sweep-from-l-41-to-l-200-t"
sys.path.insert(0, r"E:\Projects\MillenniumResearch\prize_work\ppl-014")
sys.path.insert(0, HERE)
import numpy as np
from density_threshold import omega_spf, build_W, gstar, THRESH, RMAX
from per_round_weight_sharpening import level_mass, max_layer
from bounds200 import (U_maxlevel, U_sum, L_prime, kmax_below_quarter,
                       build_W_np)
from verify200 import U_prime_indep, L_prime_indep, gstar_indep, w2_direct

N = 20_000
t0 = time.time()
om = omega_spf(N)
assert om[8] == 3 and om[12] == 3 and om[360] == 6

import array
wt = array.array('d', [0.0]) * (N + 1)
for n in range(2, N + 1):
    wt[n] = 1.0 / (n * math.log(n))
for e in (2, 3, 4):
    wt[e] = 0.0

# numpy gate: build_W_np vs pure-Python build_W
wtn = np.zeros(N + 1)
wtn[2:] = np.array(wt[2:], dtype=np.float64)
for e in (2, 3, 4):
    wtn[e] = 0.0
nn = np.arange(N + 1, dtype=np.int32)   # 0,1 never enter any r>=1 mask
omn = np.frombuffer(om, dtype=np.int8)
for L in (4, 41, 200):
    Wp = build_W(om, wt, N, L, RMAX)
    Wn = build_W_np(omn, wtn, nn, L, RMAX)
    d = max(abs(Wp[j][r] - float(Wn[j, r])) for j in range(L)
            for r in range(RMAX + 1))
    print(f"L={L:3d} np-vs-py maxdiff={d:.2e}  {'OK' if d < 1e-12 else 'FAIL'}")
    assert d < 1e-12

# G1 machinery (level_mass at 2e6 needs its own sieve; test at N=20k only
# for code-path, no ref)
f1 = max_layer(level_mass(N, om, L=4, S=(1, 2, 3)))[1]
print(f"level_mass path ok: max-layer@20k = {f1}")

# bound functions: consistency identities
W = build_W(om, wt, N, 41, RMAX)
for k in (1, 2, 5, 10):
    Up = U_maxlevel(W, 41, k, RMAX)
    U = U_sum(W, 41, k, RMAX)
    Lp = L_prime(W, 41, k, RMAX)
    Upi = U_prime_indep(W, 41, k, RMAX)
    Lpi = L_prime_indep(W, 41, k, RMAX)
    ok = (Lp <= Up and abs(Upi - Up) < 1e-15 and abs(Lpi - Lp) < 1e-15
          and Lp <= U)
    print(f"k={k:2d} Lp={Lp:.6f} Up={Up:.6f} U={U:.6f} indep_ok={ok}")
    assert ok

# kmax identity
for L in (41, 44, 45, 200):
    assert kmax_below_quarter(L) * 4 < L
    assert (kmax_below_quarter(L) + 1) * 4 >= L
print("kmax_below_quarter ok")

# gstar_indep vs shared gstar at L=4,8
for L, k in ((4, 1), (8, 2)):
    a = gstar(build_W(om, wt, N, L, RMAX), L, k, RMAX)
    b = gstar_indep(build_W(om, wt, N, L, RMAX), L, k, RMAX)
    assert abs(a - b) < 1e-15, (L, k, a, b)
    print(f"gstar_indep L={L} k={k}: {a} == shared {b}")

# w2_direct vs build_W (same omega input)
W2 = w2_direct(om, wtn, N, 41, RMAX)
Wp = build_W(om, wt, N, 41, RMAX)
d = max(abs(W2[j][r] - Wp[j][r]) for j in range(41) for r in range(RMAX + 1))
print(f"w2_direct vs build_W maxdiff={d:.2e}  {'OK' if d < 1e-12 else 'FAIL'}")
assert d < 1e-12

print(f"\nSMOKE OK  {time.time()-t0:.1f}s")
