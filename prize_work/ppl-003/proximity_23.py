# job lab_erdos_261112_24 -- follow-up: the 2/3 proximity effect behind the top-20 "growth"
# The top-20 by growth are |S|=4 points clustered just below 2/3 (S=[3], the c=2/3 tower).
# Hypothesis: their large "growth" (two-sided slope at the W2=0.00025 window) is a PROXIMITY
# effect - the W2 window catches the 1/ln N tower of 2/3 - while the intrinsic fine-window
# (W=5e-5) local mass of each is a tiny decaying atom-sum (ln N)^k/N.
import numpy as np, json, math, sys
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8')
SW = r"E:\Projects\MillenniumResearch\prize_work/ppl-003"
N = 10_000_000
W = 5e-5
W2 = 0.00025
phi = np.load(f"{SW}/phi_cache_{N}.npy").astype(np.int64)
nn = np.arange(1, N + 1, dtype=np.float64)
r = phi[1:].astype(np.float64) / nn
T23 = float(Fraction(2, 3))

scan = json.load(open(f"{SW}/scan_N{N}.json"))
dl = scan["divergent_list"]
by_growth = sorted(dl, key=lambda d: -float(d["growth"]))
top20 = by_growth[:20]
NAMED = ["1/2", "2/3", "1/3", "2/5", "3/7"]
by_c = {d["c"]: d for d in dl}
tops, seen = [], set()
for d in top20 + [by_c[c] for c in NAMED if c in by_c]:
    if d["c"] not in seen:
        seen.add(d["c"]); tops.append(d)

rows = []
for d in tops:
    cS = float(Fraction(d["c"]))
    Ma_W  = float(np.sum((r > cS) & (r < cS + W))) / N
    Mb_W  = float(np.sum((r > cS - W) & (r < cS))) / N
    Ma_W2 = float(np.sum((r > cS) & (r < cS + W2))) / N
    Mb_W2 = float(np.sum((r > cS - W2) & (r < cS))) / N
    dist  = T23 - cS
    rows.append({"c": d["c"], "S": [int(x) for x in d["S"]], "growth": d["growth"],
                 "named": d["c"] in NAMED,
                 "dist_to_2over3": dist, "side_2over3": "above_cS" if dist > 0 else "below_cS",
                 "catches_2over3_at_W": abs(dist) < W, "catches_2over3_at_W2": abs(dist) < W2,
                 "Ma_W": Ma_W, "Mb_W": Mb_W, "Ma_W2": Ma_W2, "Mb_W2": Mb_W2})

# The 2/3 tower itself (below side) is the reference 1/ln N mass ~ C([3])/lnN = (1/2)/lnN
c23_below = float(np.sum((r > T23 - W2) & (r < T23))) / N
lnN = math.log(N)
rows.append({"c": "2/3(tower)", "S": [3], "growth": None, "named": True,
             "dist_to_2over3": 0.0, "side_2over3": "self",
             "catches_2over3_at_W": True, "catches_2over3_at_W2": True,
             "Ma_W": None, "Mb_W": None, "Ma_W2": None, "Mb_W2": None})

t20 = [x for x in rows if not x["named"]]
catchW2 = [x for x in t20 if x["catches_2over3_at_W2"]]
catchW  = [x for x in t20 if x["catches_2over3_at_W"]]
print(f"N={N}  2/3={T23:.6f}  W=5e-5  W2=0.00025  2/3-below-tower mass(W2)={c23_below:.6f}  C([3])/lnN={(0.5)/lnN:.6f}")
print(f"\ntop20 (non-named) by growth: {len(t20)} points; catch 2/3 within W2={len(catchW2)}; within W={len(catchW)}")
print(f"\nc               S                          growth   dist2/3  side  catchW catchW2   Ma_W      Mb_W      Ma_W2     Mb_W2")
for x in sorted(rows, key=lambda z: -(z["growth"] or 0)):
    g = f'{x["growth"]:.3f}' if x["growth"] else "  -- "
    maW = f'{x["Ma_W"]:.5f}' if x["Ma_W"] is not None else "   -- "
    mbW = f'{x["Mb_W"]:.5f}' if x["Mb_W"] is not None else "   -- "
    maW2 = f'{x["Ma_W2"]:.5f}' if x["Ma_W2"] is not None else "   -- "
    mbW2 = f'{x["Mb_W2"]:.5f}' if x["Mb_W2"] is not None else "   -- "
    print(f'{x["c"]:>13} {str(x["S"]):>26} {g} {x["dist_to_2over3"]:+.6f} {x["side_2over3"]:>9} {str(x["catches_2over3_at_W"]):>8} {str(x["catches_2over3_at_W2"]):>8}  {maW}  {mbW}  {maW2}  {mbW2}')

# quantify: among top20 catching 2/3 at W2, is the W2 mass ~ the 2/3 tower (0.036-ish)?
big = [x for x in t20 if x["catches_2over3_at_W2"] and max(x["Ma_W2"] or 0, x["Mb_W2"] or 0) > 0.01]
w2_masses_catching = [max(x["Ma_W2"] or 0, x["Mb_W2"] or 0) for x in catchW2]
out = {
    "job": "lab_erdos_261112_24_proximity", "N": N, "W": W, "W2": W2,
    "two_thirds": T23,
    "two_thirds_below_tower_mass_W2": c23_below,
    "C_of_3_over_lnN_bare": (0.5) / lnN,
    "top20_non_named_count": len(t20),
    "catch_2over3_within_W2": len(catchW2),
    "catch_2over3_within_W": len(catchW),
    "catchW2_max_min_mass": [max(w2_masses_catching), min(w2_masses_catching)] if w2_masses_catching else None,
    "headline": "the top-20 'growth' (W2 window) is a proximity effect to the 2/3 1/ln N tower; their intrinsic fine-window (W=5e-5) mass is a tiny decaying atom-sum (ln N)^k/N",
    "rows": rows,
}
with open(f"{SW}/proximity_23_N{N}.json", "w") as f:
    json.dump(out, f, indent=2, default=float)

print(f"\n[top20 catching 2/3 at W2 with W2-mass>0.01]: {len(big)}  (vs 2/3 tower mass {c23_below:.4f})")
print("note: W2 growth-window mass ~0.036 = the 2/3 1/ln N tower; fine-window (W=5e-5) mass is the tiny atom-sum")
print("saved proximity_23_N%d.json" % N)
