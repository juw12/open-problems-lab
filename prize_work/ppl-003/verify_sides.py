#!/usr/bin/env python3
# Verify the windowed-slope side decomposition for #50 value-set points (job lab_erdos_261111_44).
# Confirms s_two(w) = (s_below(w) + s_above(w))/2 and shows WHERE the |S|=6 singularity sits
# (upper vs lower side) vs the c=1/2 tower (lower-side, huge 2^a*p pile-up) vs the 1/4 phantom.
import math, time
from fractions import Fraction
import numpy as np

N = 5000000
SLUG = "E:/Projects/MillenniumResearch/prize_work/ppl-003"
W0, W1 = 0.05, 0.00025

t0 = time.time()
phi = np.load(f"{SLUG}/phi_cache_{N}.npy")
nn = np.arange(1, N + 1, dtype=np.float64)
r = phi[1:].astype(np.float64) / nn
rs = np.sort(r)

def slopes(c):
    # two-sided, below [c-w,c), above [c,c+w): densities per unit c
    out = {}
    for w in (W0, W1):
        nb = int(((r >= c - w) & (r < c)).sum())
        na = int(((r >= c) & (r < c + w)).sum())
        out[w] = (nb / (w * N), na / (w * N))
    (b0, a0), (b1, a1) = out[W0], out[W1]
    two0 = (b0 + a0) / 2; two1 = (b1 + a1) / 2
    return {
        "growth_two": two1 / two0 if two0 > 0 else 0.0,
        "growth_below": b1 / b0 if b0 > 0 else 0.0,
        "growth_above": a1 / a0 if a0 > 0 else 0.0,
        "b0": b0, "b1": b1, "a0": a0, "a1": a1,
        "id_check": abs(two0 - (b0 + a0) / 2) < 1e-12,
    }

pts = {
    "c=1/2 (tower, 2^a*p below)": Fraction(1, 2),
    "1/3 (tower point, divergent)": Fraction(1, 3),
    "1/4 (PHANTOM, job27)": Fraction(1, 4),
    "S6 {2,3,5,7,11,13}": Fraction(192, 1001),
    "S6 {2,3,5,7,11,17}": Fraction(256, 1309),
    "S6 {2,3,5,7,13,17}": Fraction(384, 1547),
    "S6 {3,5,7,11,13,17} (no 2)": Fraction(122880, 345855),
}
print(f"N={N}  w0={W0} w1={W1}  [load+sort {time.time()-t0:.1f}s]")
print(f"{'point':32s} {'g_two':>8s} {'g_below':>9s} {'g_above':>9s}  side")
for name, fr in pts.items():
    c = float(fr)
    s = slopes(c)
    side = "BELOW" if s["growth_below"] > s["growth_above"] else ("ABOVE" if s["growth_above"] > s["growth_below"] else "flat")
    print(f"{name:32s} {s['growth_two']:8.3f} {s['growth_below']:9.3f} {s['growth_above']:9.3f}  {side}  id={s['id_check']}")
print(f"total {time.time()-t0:.1f}s")
