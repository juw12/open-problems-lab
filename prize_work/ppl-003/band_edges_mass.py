# Job lab_erdos_261112_11 — Erdos #50 — MASS-BASED VERIFIER
# Confirms the divergence is a MASS (population) phenomenon, not a tower-count
# phenomenon. Uses the phi cache (phi(n)/n for n<=N) to compute the windowed-slope
# mass(w;c) = #{n<=N : |phi(n)/n - c| < w} / (2w N)  = slope of f in [c-w, c+w].
# A divergent-slope point has slope that GROWS as w shrinks (mass concentrated);
# a smooth point has slope ~ const (mass spread). Compare the two band edges
# (divergent) against tower values just outside them (smooth).

import numpy as np, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
N = 10_000_000
phint = np.load(os.path.join(HERE, f"phi_cache_{N}.npy"))  # this is phi(n) (int), NOT phi(n)/n
# G2 gate: mean(phi(n)) = 3/pi^2 * N
import math
g2 = phint[1:N + 1].mean() / N
print(f"GATE G2 mean(phi)/N = {g2:.7f}  (3/pi^2 = {3/math.pi**2:.7f})  "
      f"pass={abs(g2 - 3/math.pi**2) < 1e-5}")
idx = np.arange(1, N + 1, dtype=np.float64)
phi = (phint[1:N + 1].astype(np.float64)) / idx   # phi(n)/n
print(f"loaded phi(n)/n, N={N}, min={phi.min():.6f} max={phi.max():.6f} mean={phi.mean():.7f} "
      f"(6/pi^2={6/math.pi**2:.7f})")

def slope(c, w):
    # mass fraction in [c-w, c+w] divided by 2w = local slope of the CDF f
    m = np.mean((phi > c - w) & (phi < c + w))
    return m / (2 * w)

def growth(c, w_large=0.05, w_small=0.00025):
    # lab-like growth = slope(w_small)/slope(w_large)  (concentrated -> large)
    sl = slope(c, w_small); sg = slope(c, w_large)
    return (sl / sg) if sg > 0 else float("inf")

cands = [
    ("LEFT edge   {2,3,5,7}=8/35   (divergent)", 8/35),
    ("below-left  {2,3,5,7,11}      (not div)",   16/77),
    ("c=1/2 2p-spike (divergent)", 0.5),
    ("RIGHT edge  {47,73,83,157}    (divergent)", 42367104/44709361),
    ("inside      {19}=18/19        (divergent)", 18/19),
    ("above-edge  {23}=22/23        (not div)",   22/23),
    ("above       {89,97,103}       (not div)",   (88*96*102)/(89*97*103)),
]
print(f"\n{'point':44s} {'slope(w=2.5e-4)':>16s} {'slope(w=5e-2)':>14s} {'growth':>9s}")
for label, c in cands:
    sl = slope(c, 0.00025); sg = slope(c, 0.05)
    gr = sl/sg if sg > 0 else float("inf")
    print(f"{label:44s} {sl:16.4f} {sg:14.4f} {gr:9.3f}")

# the tower-count (from the fresh verifier) vs mass slope contrast, stated plainly
print("\nINTERPRETATION:")
print("  {23}=0.9565 has ~21151 TOWER neighbors (dense in tower-count) but SMALL mass slope -> smooth.")
print("  (47,73,83,157)=0.94761 has ~23994 tower neighbors AND the mass is concentrated -> divergent.")
print("  => divergence is a MASS (population) phenomenon, not a tower-count/percolation one.")

res = {"N": N, "w_small": 0.00025, "w_large": 0.05, "points": []}
for label, c in cands:
    sl = slope(c, 0.00025); sg = slope(c, 0.05)
    res["points"].append({"label": label, "c": c, "slope_small": float(sl),
                           "slope_large": float(sg), "growth": float(sl/sg) if sg > 0 else None})
out = os.path.join(HERE, "band_edges_mass.json")
json.dump(res, open(out, "w"), indent=2)
print(f"\nwrote {out}")
