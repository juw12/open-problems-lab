# Job lab_erdos_261112_11 — Erdos #50 (PPL 003) — MAIN RESULT
# Pin the two edges of the interior singular band to EXACT prime-tower values and
# characterize the band as the "percolated" (well-accumulated) core of the
# prime-tower value set. Independent cross-check: the grid-based interior_map
# (job 261112_5, different machinery) c_range, N-stable to 1e-3.
#
# All divergent points are prime-tower values c_S = prod_{p in S}(1-1/p) for the
# prime set S. The singular band = tower values that are accumulation points.
# The |S|<=4 scan is COMPLETE (all 3008 divergent points stored), so the edges
# are pinned exactly from it; the |S|>=5 scans store only top-by-growth, so we
# use the N-stable grid c_range for the |S|>=5 tail.

import json, os, math
from fractions import Fraction
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))

def tower(S):
    n = 1; d = 1
    for p in S:
        n *= (p - 1); d *= p
    return Fraction(n, d)

# ---- Load the COMPLETE |S|<=4 divergent set (the authoritative edge source) ----
d4 = json.load(open(os.path.join(HERE, "scan_N10000000.json")))
S4 = {}   # S -> (c, class, growth)
for e in d4["divergent_list"]:
    cstr = e["c"]; a, b = cstr.split("/"); cs = Fraction(int(a), int(b))
    S = tuple(sorted(e["S"]))
    S4[S] = (cs, e.get("class"), e.get("growth"))

# also load the top-by-growth subsets for |S|=5..8 (for local-density, biased but useful)
extras = {}
for fname, key in [("scan_s5_N10000000.json", "divergent_detail"),
                   ("scan_s6_N10000000.json", "top_detail"),
                   ("scan_s7_N10000000.json", "top_detail"),
                   ("scan_s8_N5000000k8_8.json", "top_detail_s8")]:
    d = json.load(open(os.path.join(HERE, fname)))
    for e in d.get(key, []):
        cstr = e.get("c") or e.get("c_reduced"); S = tuple(sorted(e["S"]))
        a, b = cstr.split("/"); cs = Fraction(int(a), int(b))
        extras[S] = (cs, e.get("class"), e.get("growth"), fname, d.get("N"))

# ---- GATE: exact c_S == stored reduced rational, |S|<=4 (complete) ----
gate_mismatch = sum(1 for S, (cs, cl, g) in S4.items() if tower(S) != cs)
print(f"GATE |S|<=4 exact c_S==stored: {len(S4)} pts, mismatch={gate_mismatch}")

# ---- EDGES from the complete |S|<=4 set ----
order4 = sorted(S4.items(), key=lambda kv: kv[1][0])
left_S, (left_c, left_cl, left_g) = order4[0]
right_S, (right_c, right_cl, right_g) = order4[-1]
print("\n=== LEFT EDGE (min, |S|<=4 complete) ===")
print(f"S={list(left_S)}  c_exact={left_c} = {float(left_c):.6f}  class={left_cl}  growth={left_g}")
print("=== RIGHT EDGE (max, |S|<=4 complete) ===")
print(f"S={list(right_S)}  c_exact={right_c} = {float(right_c):.8f}  class={right_cl}  growth={right_g}")

# ---- divergent single primes (which sit INSIDE the band, not at the edge) ----
singles = sorted([(int(S[0]), S) for S in S4 if len(S) == 1])
print("\ndivergent single primes (|S|=1):", [p for p, _ in singles])
# the first single prime ABOVE the right edge (isolated / smooth)
above = [p for p, _ in singles if Fraction(p - 1, p) > right_c]
print("single primes above right edge (isolated):", above)
# show {23} specifically
print(f"{{23}} = {float(Fraction(22,23)):.6f}  (above right edge {float(right_c):.6f} -> outside band)")

# ---- right-edge region: confirm all 12 largest are |S|=4 large-prime towers ----
print("\n12 largest |S|<=4 divergent towers (right-edge region):")
for S, (cs, cl, g) in order4[-12:]:
    print(f"  c={float(cs):.8f}  S={list(S)}  |S|={len(S)}  growth={g}")

# ---- LOCAL-DENSITY / percolation profile over the |S|<=4 complete set ----
# For each c in [0.2, 0.99] on a 1e-3 grid, count |S|<=4 divergent towers within w.
w = 5e-4
cs_vals = sorted(float(cs) for cs, cl, g in S4.values())
n_grid = 800
lo, hi = 0.20, 0.99
dens = []
for i in range(n_grid):
    c0 = lo + (hi - lo) * i / (n_grid - 1)
    k = sum(1 for v in cs_vals if abs(v - c0) <= w)
    dens.append((c0, k))
# band = contiguous c region where density >= threshold (percolated core)
import statistics
maxk = max(k for c, k in dens)
thr = max(3, int(0.02 * maxk))
inband = [c for c, k in dens if k >= thr]
if inband:
    band_lo, band_hi = min(inband), max(inband)
    print(f"\npercolation profile: w={w}, max_local_density={maxk}, thr={thr}")
    print(f"  percolated-core span (|S|<=4 complete): [{band_lo:.4f}, {band_hi:.4f}]")
# density just OUTSIDE each edge (show isolation)
def dens_at(c0):
    return sum(1 for v in cs_vals if abs(v - c0) <= w)
print(f"  local density at left edge {float(left_c):.4f}: {dens_at(float(left_c))}")
print(f"  local density 0.05 below left edge:          {dens_at(float(left_c)-0.05)}")
print(f"  local density at right edge {float(right_c):.4f}: {dens_at(float(right_c))}")
print(f"  local density 0.05 above right edge:         {dens_at(float(right_c)+0.05)}")

# ---- cross-check: N-stable grid c_range from the independent interior_map ----
gm = json.load(open(os.path.join(HERE, "interior_map_full_N10000000.json")))
h = gm.get("headline", {})
print("\n=== INDEPENDENT CROSS-CHECK (grid-based interior_map, job 261112_5) ===")
print("N=1e7 grid c_range:", h.get("c_range"), " n_div_cells:", h.get("n_div_fine_cells"), " measure:", h.get("measure_of_[0,1]"))
# N=5e6 grid
gm5 = None
for cand in ["interior_map_full_N5000000.json"]:
    p = os.path.join(HERE, cand)
    if os.path.exists(p):
        gm5 = json.load(open(p)); break
if gm5:
    print("N=5e6 grid c_range:", gm5.get("headline", {}).get("c_range"))

# diff of exact edges vs grid
print(f"\nexact left  {float(left_c):.6f}  vs grid min {h['c_range'][0]:.5f}  diff={float(left_c)-h['c_range'][0]:+.2e}")
print(f"exact right {float(right_c):.6f}  vs grid max {h['c_range'][1]:.5f}  diff={float(right_c)-h['c_range'][1]:+.2e}")

# ---- write result JSON ----
res = {
    "job": "lab_erdos_261112_11",
    "N": 10000000,
    "gate": {"n_S4": len(S4), "mismatch": gate_mismatch, "pass": gate_mismatch == 0},
    "left_edge": {"S": list(left_S), "c_exact": str(left_c), "c_float": float(left_c),
                   "class": left_cl, "growth": left_g},
    "right_edge": {"S": list(right_S), "c_exact": str(right_c), "c_float": float(right_c),
                   "class": right_cl, "growth": right_g,
                   "caveat": "|S|<=4 complete; |S|>=5 tail via N-stable grid (rightmost stored |S|<=4 tower)"},
    "divergent_single_primes_inside_band": [p for p, _ in singles if Fraction(p - 1, p) < right_c],
    "first_single_prime_above_edge": min(above) if above else None,
    "band": {"percolation_span_S4": [band_lo, band_hi] if inband else None,
             "grid_c_range_N1e7": h.get("c_range"),
             "grid_measure": h.get("measure_of_[0,1]"),
             "grid_n_div_cells": h.get("n_div_fine_cells")},
    "percolation": {"w": w, "max_local_density_S4": maxk, "thr": thr,
                    "dens_left_edge": dens_at(float(left_c)),
                    "dens_left_minus_0p05": dens_at(float(left_c) - 0.05),
                    "dens_right_edge": dens_at(float(right_c)),
                    "dens_right_plus_0p05": dens_at(float(right_c) + 0.05)},
}
out = os.path.join(HERE, "band_edges_result.json")
json.dump(res, open(out, "w"), indent=2)
print(f"\nwrote {out}")
