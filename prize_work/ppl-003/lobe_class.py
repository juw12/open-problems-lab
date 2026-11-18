#!/usr/bin/env python3
"""Job lab_erdos_261112_21 - PART 1: class composition of the c~0.65 lobe.

Decompose the divergent-slope points in [0.60,0.70] (the c~0.65 lobe, the
highest-density interior region of the Schoenberg f(c) band) by |S| class and
2p/has-2 status. Reuses the scan_s7 machinery (enumerate_arrays, scan_growth,
class_of) on the cached exact int64 totient.

GATE (proves the re-run is correct): exact match of the GLOBAL
divergent_by_absS and divergent_by_class against the committed
scan_s7_N10000000.json (job 44 numbers). If the global counts match exactly,
the [0.60,0.70] filter (a pure mask on cvals over the same div_idx) is exact.

INDEPENDENT VERIFIER (different code path, not the same code twice):
two-sided windowed-slope growth recomputed with a DIRECT boolean mask on the
UNSORTED r array (scan_growth uses searchsorted on the SORTED r) for the
top-K divergent [0.60,0.70] points by growth. Same window convention
[c-w, c+w).

KILL TEST (pre-test): is the c~0.65 lobe a BAND or a single isolated c-value
spike in the committed job-261112_5 map? Count the 1e-5 fine cells in
[0.62632,0.67684]; a single cell = sub-grid isolated -> kill.

Empirical/heuristic only. No proof or solution is claimed.
"""
import numpy as np, json, sys, time
from pathlib import Path

sys.path.insert(0, "E:/Projects/MillenniumResearch/prize_work/ppl-003")
import scan_s7 as S

N = 10_000_000
SLO = S.SLUG
OUT = Path(SLO) / f"lobe_class_N{N}.json"
W0 = S.WINDOWS[0]   # 0.05
W1 = S.WINDOWS[-1]  # 0.00025
TH = S.DIVERGE_THRESH  # 5.0

t0 = time.time()
lo_core, hi_core = 0.62632, 0.67684   # bin12 (the c~0.65 lobe core, highest profile bin)

# ------------------------------------------------------------------ LOAD r
phi = np.load(Path(SLO) / "phi_cache_10000000.npy")  # shape (N+1,), phi[0]=0
n = np.arange(1, N + 1, dtype=np.int64)
r = (phi[1:] / n).astype(np.float64)   # phi(n)/n for n=1..N
rs = np.sort(r)
print(f"[{time.time()-t0:5.1f}s] loaded r (N={N})")

# ------------------------------------------------- enumerate + growth (|S|<=7)
ps = S.primes_upto(S.PRIMES_MAX)
P = [int(q) for q in range(2, S.PRIMES_MAX + 1) if ps[q]]
cvals, Sarr, has2, absS = S.enumerate_arrays(P, 7)
M = len(cvals)
print(f"[{time.time()-t0:5.1f}s] enumerated {M} candidates |S|<=7")

growth = S.scan_growth(rs, N, cvals)
div = growth >= TH
div_idx = np.flatnonzero(div)
print(f"[{time.time()-t0:5.1f}s] divergent (growth>={TH}): {len(div_idx)}")

# --------------------------------------------------------------- KILL TEST
# Band vs isolated: count distinct 1e-5 fine cells occupied by the DIVERGENT
# points inside the lobe core [0.62632,0.67684]. A single cell = sub-grid
# isolated spike -> kill; many cells = a band -> proceed. (From our own
# div_idx, i.e. the actual divergent c-values, not a stored profile.)
cc = cvals[div_idx][(cvals[div_idx] >= lo_core) & (cvals[div_idx] <= hi_core)]
cells = np.unique(np.floor(cc / 1e-5).astype(np.int64))
n_core_cells = int(len(cells))
core_cmin, core_cmax = (float(cc.min()), float(cc.max())) if len(cc) else (0.0, 0.0)
kill_isolated = (n_core_cells <= 1)
print(f"[kill test] divergent fine cells (1e-5) in [{lo_core},{hi_core}]: "
      f"n={n_core_cells}, c in [{core_cmin:.5f},{core_cmax:.5f}]  -> "
      f"{'ISOLATED -> KILL' if kill_isolated else 'BAND -> proceed'}")

# ---------------------------------------------------------------- GATE (exact)
com = json.load(open(Path(SLO) / "scan_s7_N10000000.json"))
com_absS = {int(k): v for k, v in com["divergent_by_absS"].items()}

by_absS = {}
by_class = {}
for i in div_idx:
    k = int(absS[i]); by_absS[k] = by_absS.get(k, 0) + 1
    cl = S.class_of(absS[i], bool(has2[i])); by_class[cl] = by_class.get(cl, 0) + 1

gate_total = (len(div_idx) == com["n_divergent"])
gate_absS = (by_absS == com_absS)
gate_class = (by_class == com["divergent_by_class"])
print(f"[{time.time()-t0:5.1f}s] GATE total={gate_total} by_absS={gate_absS} by_class={gate_class}")
if not (gate_total and gate_absS and gate_class):
    print("  GATE MISMATCH detail:")
    for k in sorted(set(list(by_absS) + list(com_absS))):
        if by_absS.get(k) != com_absS.get(k):
            print(f"    absS {k}: {by_absS.get(k)} vs {com_absS.get(k)}")
    for k in sorted(set(list(by_class) + list(com["divergent_by_class"]))):
        if by_class.get(k) != com["divergent_by_class"].get(k):
            print(f"    class {k}: {by_class.get(k)} vs {com['divergent_by_class'].get(k)}")

# ------------------------------------------------------- decomposition helper
def decomp(lo, hi):
    m = (cvals[div_idx] >= lo) & (cvals[div_idx] <= hi)
    sub = div_idx[m]; nsub = int(len(sub))
    babs, bcls = {}, {}
    h2 = has2[sub]; aS = absS[sub]
    for j in range(nsub):
        k = int(aS[j]); babs[k] = babs.get(k, 0) + 1
        cl = S.class_of(aS[j], bool(h2[j])); bcls[cl] = bcls.get(cl, 0) + 1
    n_has2 = int(h2.sum()); n_not2 = nsub - n_has2
    dom = max(bcls, key=bcls.get) if bcls else None
    return {
        "range": [lo, hi], "n_divergent": nsub,
        "n_has2": n_has2, "n_not_has2": n_not2,
        "frac_has2": (n_has2 / nsub) if nsub else 0.0,
        "by_absS": {str(k): babs[k] for k in sorted(babs)},
        "by_class": dict(sorted(bcls.items(), key=lambda x: -x[1])),
        "dominant_class": dom,
        "dominant_class_count": (bcls[dom] if dom else 0),
        "dominant_class_share": (bcls[dom] / nsub if (dom and nsub) else 0.0),
        "top3_classes": [[k, v] for k, v in sorted(bcls.items(), key=lambda x: -x[1])[:3]],
    }

d070 = decomp(0.60, 0.70)
dcore = decomp(lo_core, hi_core)
print(f"\n[0.60,0.70]: n={d070['n_divergent']} dominant={d070['dominant_class']} "
      f"({d070['dominant_class_count']}, share={d070['dominant_class_share']:.4f}) "
      f"has2={d070['n_has2']} ({d070['frac_has2']:.3f})")
print("   by_class:", d070["by_class"])
print("   by_absS:", d070["by_absS"])
print(f"[{lo_core},{hi_core}]: n={dcore['n_divergent']} dominant={dcore['dominant_class']} "
      f"({dcore['dominant_class_count']}, share={dcore['dominant_class_share']:.4f})")
print("   by_class:", dcore["by_class"])

# ------------------------------------- INDEPENDENT VERIFIER (boolean on r)
sub07 = div_idx[(cvals[div_idx] >= 0.60) & (cvals[div_idx] <= 0.70)]
sub07 = sub07[np.argsort(-growth[sub07])]   # top by growth
K = 50
top = sub07[:K]
vcheck = []; mismatch = 0; n_div_indep = 0
for i in top:
    c = float(cvals[i])
    s0 = float(np.count_nonzero((r >= c - W0) & (r < c + W0))) / (2 * W0 * N)
    s1 = float(np.count_nonzero((r >= c - W1) & (r < c + W1))) / (2 * W1 * N)
    gi = (s1 / s0) if s0 > 0 else np.nan
    ok = (not np.isnan(gi)) and (not np.isnan(growth[i])) and \
         abs(gi - growth[i]) <= 0.10 * max(abs(float(growth[i])), 1.0)
    if not ok: mismatch += 1
    if (not np.isnan(gi)) and gi >= TH: n_div_indep += 1
    vcheck.append({
        "c": round(c, 6), "absS": int(absS[i]),
        "class": S.class_of(absS[i], bool(has2[i])),
        "growth_scan": round(float(growth[i]), 4),
        "growth_indep": None if np.isnan(gi) else round(float(gi), 4),
        "match10pct": bool(ok), "div_indep": bool(gi >= TH) if not np.isnan(gi) else None,
    })
print(f"[{time.time()-t0:5.1f}s] verifier: K={len(top)}, match(10%)={len(top)-mismatch}/{len(top)}, "
      f"still-divergent(indep)={n_div_indep}/{len(top)}")

out = {
    "job_id": "lab_erdos_261112_21", "part": 1,
    "N": N, "w_large": W0, "w_small": W1, "diverge_threshold": TH,
    "kill_test": {"core_range": [lo_core, hi_core],
                  "n_divergent_fine_cells_1e-5": n_core_cells,
                  "core_cmin": round(core_cmin, 6), "core_cmax": round(core_cmax, 6),
                  "isolated": kill_isolated, "proceed": (not kill_isolated)},
    "gate": {"n_candidates": int(M), "n_divergent": int(len(div_idx)),
             "n_divergent_committed": com["n_divergent"],
             "match_total": gate_total, "match_by_absS": gate_absS, "match_by_class": gate_class,
             "PASS": bool(gate_total and gate_absS and gate_class),
             "by_absS_mine": {str(k): v for k, v in sorted(by_absS.items())},
             "by_class_mine": dict(sorted(by_class.items(), key=lambda x: -x[1]))},
    "decomp_060_070": d070,
    "decomp_core_bin12": dcore,
    "independent_verifier": {"method": "direct boolean two-sided window on UNSORTED r "
        "(scan_growth uses searchsorted on sorted r)", "K": int(len(top)),
        "match_10pct": int(len(top) - mismatch), "still_divergent_indep": int(n_div_indep),
        "points": vcheck},
    "wall_s": round(time.time() - t0, 1),
}
OUT.write_text(json.dumps(out, indent=2))
print(f"\nwrote {OUT}   (wall {time.time()-t0:.1f}s)")
