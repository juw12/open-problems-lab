# Job lab_erdos_261112_11 — Erdos #50 (PPL 003)
# Explore the EXACT prime-tower values of the interior singular set and pin its
# two EDGES (min and max c over all divergent-slope points) to exact rationals.
#
# Each divergent-slope point in the lab scans is a prime-tower value
#   c_S = prod_{p in S} (1 - 1/p) = prod(p-1)/prod(p)
# for the prime set S (the set of prime divisors of n; phi(n)/n depends only on
# this set). The singular band [0.22827, 0.94761] (job 261112_5) is the range of
# tower values that are accumulation points of the tower-value set.
#
# Here we: load every divergent point from the committed scans (|S|=1..8), recompute
# the EXACT rational c_S from S, gate it against the lab's stored reduced 'c',
# and report the global min and max tower value (the two band edges).

import json, os, math
from fractions import Fraction
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))

def load_divergent(fname):
    """Return list of (c_float, S, class, absS, growth, source) from a scan json."""
    p = os.path.join(HERE, fname)
    d = json.load(open(p))
    N = d.get("N")
    out = []
    def grab(entry):
        cstr = entry.get("c") or entry.get("c_reduced")
        S = entry.get("S")
        cl = entry.get("class")
        g = entry.get("growth")
        if cstr is None or S is None:
            return None
        a, b = cstr.split("/")
        return (Fraction(int(a), int(b)), tuple(sorted(S)), cl, g, N, fname)
    # try the known divergent-list keys in order
    for key in ["divergent_list", "divergent_detail", "top_detail",
                "top_detail_s8", "top_s7_by_growth", "top_s8_by_growth"]:
        if key in d and isinstance(d[key], list):
            for e in d[key]:
                r = grab(e)
                if r:
                    out.append(r)
    return out

def tower(S):
    num = 1; den = 1
    for p in S:
        num *= (p - 1); den *= p
    return Fraction(num, den)

allpts = []
for fname in ["scan_N10000000.json", "scan_s5_N10000000.json",
              "scan_s6_N10000000.json", "scan_s7_N10000000.json",
              "scan_s8_N5000000k8_8.json"]:
    pts = load_divergent(fname)
    allpts += pts
    print(f"{fname}: {len(pts)} divergent points loaded")

# Gate: recompute exact c_S from S and compare to the stored reduced rational.
mismatch = 0; matched = 0
for (cstored, S, cl, g, N, fname) in allpts:
    ct = tower(S)
    if ct == cstored:
        matched += 1
    else:
        mismatch += 1
        if mismatch <= 5:
            print("MISMATCH", fname, "stored", cstored, "tower(S)", ct, "S", S)
print(f"\nGATE exact c_S==stored: matched={matched} mismatch={mismatch} "
      f"({100*matched/(matched+mismatch):.4f}%)")

# Deduplicate by (S) keeping one record; compute global min/max c.
byS = {}
for (cstored, S, cl, g, N, fname) in allpts:
    ct = tower(S)
    if S not in byS:
        byS[S] = (ct, cl, g, N, fname)
print(f"unique prime sets S: {len(byS)}")

# global min and max tower value
smin = min(byS.items(), key=lambda kv: kv[1][0])
smax = max(byS.items(), key=lambda kv: kv[1][0])
print("\n=== LEFT EDGE (min c over divergent towers) ===")
print("S =", smin[0], " c_exact =", smin[1][0], "=", float(smin[1][0]),
      " class", smin[1][1], "growth", smin[1][2], "N", smin[1][3])
print("=== RIGHT EDGE (max c over divergent towers) ===")
print("S =", smax[0], " c_exact =", smax[1][0], "=", float(smax[1][0]),
      " class", smax[1][1], "growth", smax[1][2], "N", smax[1][3])

# Hypothesis checks against named towers
def named(Sname):
    return tower(tuple(Sname))
print("\n=== hypothesis values ===")
for label, S in [("{2,3,5,7}", [2,3,5,7]), ("{2,3,5,7,11}", [2,3,5,7,11]),
                 ("{19}", [19]), ("{17}", [17]), ("{23}", [23]),
                 ("{2,19}", [2,19])]:
    f = named(S)
    print(f"{label}: {f} = {float(f):.6f}")

# Sort divergent towers by c and dump the extreme 12 on each side
order = sorted(byS.items(), key=lambda kv: kv[1][0])
print("\n=== 12 smallest divergent towers (left edge region) ===")
for S, (ct, cl, g, N, fname) in order[:12]:
    print(f"  c={float(ct):.6f}  S={S}  |S|={len(S)}  growth={g}  N={N}")
print("\n=== 12 largest divergent towers (right edge region) ===")
for S, (ct, cl, g, N, fname) in order[-12:]:
    print(f"  c={float(ct):.6f}  S={S}  |S|={len(S)}  growth={g}  N={N}")

# Distribution of |S| and of max-prime among divergent towers
from collections import Counter
absS = Counter(len(S) for S in byS)
print("\n|S| distribution of divergent towers:", dict(sorted(absS.items())))
maxp = Counter(max(S) for S in byS)
print("max-prime distribution (top 12):", maxp.most_common(12))
