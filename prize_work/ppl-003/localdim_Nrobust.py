# -*- coding: utf-8 -*-
"""
Erdos #50 / PPL 003 - N-robustness (convergence-vs-drift) test for the job-38
core local dimension h~0.08.

Question (idea card): is the core local dimension h~0.08 a LIMIT quantity
(h(c;N) -> a stable value as N->inf), or FINITE-N DRIFT (h(c;N) keeps moving
with N, no limit on the reachable range)?

Method (mirrors job-38 exactly, so the gate reproduces its 60-core):
  For each STRONG sharp singularity c, the one-sided-below local dimension
  h(c;N) = p_c = 1 + gamma_c,  gamma = slope of log s_-(w) vs log w over the
  FIXED window range W_FINE=[5e-5,5e-3] (job-38's asymptotic regime, above the
  N*eps=10 finite-N floor). M_-(c;w)=F(c)-F(c-w)~w^p_c.
  Compute h(c;N) at N = 1e6, 2e6, 5e6, 1e7 (cached + freshly sieved, SAME window
  range at every N for an apples-to-apples N-comparison) and test:
    * GATE: reproduce job-38's 60-core at N=1e7 (n=60, p mean 0.0804, sd 0.0147,
      range [0.0659,0.1177]) + the 4 scalar gates.
    * CONVERGENCE (card): a limit local dimension exists iff h(c;N) is stable
      within +/-0.01 across N=2e6..1e7 for >=1 point.
    * KILL (card): top-3 points at N=1e6 vs N=2e6 only - if ALL three h differ
      by >0.05 between the two N, no convergence on this range -> finite-N drift.
    * TREND: do the N-increments |h(N)-h(N/2)| shrink (=> convergence, limit
      needs larger N) or stay large / oscillate (=> drift, no limit)?
  Independent VERIFIER (not the same code twice):
    (i)  DIFFERENT WINDOWING: re-estimate h with W_WIDE=[5e-5,1e-2] (a different
         slope window) and check it agrees with the W_FINE estimate.
    (ii) CO-SIEVED SECOND N: freshly sieve N=8e6 (between 5e6 and 1e7, an N NOT
         in the main set) and check h(c;8e6) vs h(c;1e7) - adjacent large N must
         agree for a stable limit; disagreement = drift even at large N.

Empirical/heuristic only - NOT a proof, NOT a solution, no claim on #50, no Lean.
"""
import json, math, os
import numpy as np
from fractions import Fraction

# Canonical data home (job-35 scan list + cached exact totient sieves). The script
# is location-independent: it reads data from DATA_DIR and writes its JSON to CWD.
DATA_DIR = r"E:\Projects\MillenniumResearch\prize_work\ppl-003"

W_FINE = np.geomspace(5e-5, 5e-3, 12)   # job-38 primary regime (same at every N)
W_WIDE = np.geomspace(5e-5, 1e-2, 12)   # different-windowing verifier
STRONG = 20.0
TOPK   = 25           # "top ~20 strong singularities" (sharpest local singularities)
N_MAIN = [1_000_000, 2_000_000, 5_000_000, 10_000_000]
N_VER  = 8_000_000    # co-sieved second N (verifier)

def totient_sieve(N):
    phi = np.arange(N+1, dtype=np.int64)
    for i in range(2, N+1):
        if phi[i] == i:            # i prime
            phi[i::i] -= phi[i::i] // i
    return phi

def load_or_sieve(N, path=None):
    import os
    if path is not None and os.path.exists(path):
        phi = np.load(path, mmap_mode='r').astype(np.int64)
    else:
        phi = totient_sieve(N)
    n = np.arange(1, len(phi), dtype=np.int64)
    r = (phi[1:].astype(np.float64) / n.astype(np.float64)).astype(np.float64)
    rs = np.sort(r)
    return phi, r, rs

# ---- job-38 machinery, copied verbatim for a faithful gate ----
def F(rs, c):
    return float(np.searchsorted(rs, c, 'left')) / len(rs)

def Fvec(rs, arr):
    return np.searchsorted(rs, arr, 'left') / len(rs)

def slope_growth(rs, c):
    W = np.geomspace(0.05, 0.00025, 11)
    def s(w):
        return (F(rs, c + w) - F(rs, c - w)) / (2.0 * w)
    return s(W[-1]) / s(W[0])

def onerow(rs, c, W):
    cw = c - W
    Mcw = F(rs, c) - Fvec(rs, cw)
    Mcw = np.clip(Mcw, 0.0, None)
    return Mcw, Mcw / W

def fit_gamma(W, s):
    m = s > 0
    if m.sum() < 4:
        return None, None
    x = np.log(W[m]); y = np.log(s[m])
    A = np.polyfit(x, y, 1)
    yhat = np.polyval(A, x)
    ss_res = float(np.sum((y - yhat) ** 2)); ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None
    return float(A[0]), r2

def h_of(rs, c, W, N, mincount=20):
    Mcw, s = onerow(rs, c, W)
    m_fine = float(Mcw[-1])
    gamma, r2 = fit_gamma(W, s)
    p = (1.0 + gamma) if gamma is not None else None
    resolves = bool(m_fine * N >= mincount)
    return p, gamma, r2, resolves

def main():
    out = {"job": "lab_erdos_261112_22", "W_FINE": [float(x) for x in W_FINE],
           "W_WIDE": [float(x) for x in W_WIDE], "strong_thresh": STRONG, "topk": TOPK}

    data = {}
    for N in N_MAIN + [N_VER]:
        path = (os.path.join(DATA_DIR, "phi_cache_%d.npy" % N)) if N in (5_000_000, 10_000_000) else None
        phi, r, rs = load_or_sieve(N, path)
        data[N] = (phi, r, rs)
        out["N_%d" % N] = {"source": ("cache" if path else "sieved"),
                           "mean_phiN": round(float(np.mean(r)), 6),
                           "expect": round(6.0/math.pi**2, 6)}
    mean1e7 = data[10_000_000][1]
    out["gate_mean_phiN_1e7"] = round(float(np.mean(mean1e7)), 6)

    # ================= GATE: reproduce job-38 at N=1e7 =================
    rs1 = data[10_000_000][2]
    j35 = json.load(open(os.path.join(DATA_DIR, "scan_N10000000.json")))
    dlist = j35["divergent_list"]
    strong = [r for r in dlist if r["growth"] >= STRONG]
    out["n_strong"] = len(strong)

    M12 = F(rs1, 0.5) - F(rs1, 0.5 - 5e-5)
    g12 = slope_growth(rs1, 0.5)
    gam12, r2_12 = fit_gamma(W_FINE, onerow(rs1, 0.5, W_FINE)[1])
    out["gate"] = {
        "mean_phiN": round(float(np.mean(mean1e7)), 6), "mean_pass": abs(float(np.mean(mean1e7)) - 6.0/math.pi**2) < 5e-4,
        "M12": round(float(M12), 6), "M12_pass": abs(float(M12) - 0.072339) < 3e-4,
        "growth12": round(g12, 3), "growth12_pass": abs(g12 - 58.46) < 1.0,
        "gamma12": round(gam12, 4), "r2_12": round(r2_12, 4),
        "gamma12_pass": (-1.2 <= gam12 <= -0.6) and (r2_12 or 0) > 0.95,
    }

    # recompute the 60-core exactly as job-38 (clean = resolves & r2>0.9)
    clean1e7 = []
    for r in strong:
        cf = float(Fraction(r["c"]))
        p, g, r2, res = h_of(rs1, cf, W_FINE, 10_000_000, 20)
        if res and r2 is not None and r2 > 0.9 and p is not None:
            clean1e7.append((r["c"], r["class"], p, g, r2))
    out["gate"]["n_clean_1e7"] = len(clean1e7)
    core = [x for x in clean1e7 if x[2] <= 0.15]
    cp = np.array([x[2] for x in core])
    out["gate"]["core_n"] = int(len(cp))
    out["gate"]["core_mean"] = round(float(cp.mean()), 4)
    out["gate"]["core_sd"] = round(float(cp.std()), 4)
    out["gate"]["core_min"] = round(float(cp.min()), 4)
    out["gate"]["core_max"] = round(float(cp.max()), 4)
    out["gate"]["core_expect"] = {"n": 60, "mean": 0.0804, "sd": 0.0147, "min": 0.0659, "max": 0.1177}
    out["gate"]["core_pass"] = bool(abs(cp.mean() - 0.0804) < 0.003 and abs(cp.std() - 0.0147) < 0.005
                                     and abs(cp.min() - 0.0659) < 0.004 and abs(cp.max() - 0.1177) < 0.006
                                     and abs(len(cp) - 60) <= 2)
    out["gate"]["all_pass"] = bool(out["gate"]["mean_pass"] and out["gate"]["M12_pass"]
                                    and out["gate"]["growth12_pass"] and out["gate"]["gamma12_pass"]
                                    and out["gate"]["core_pass"])

    # ================= SELECT top sharp points at N=1e7 (smallest p) =================
    ranked = sorted(clean1e7, key=lambda x: x[2])[:TOPK]
    sel = [x[0] for x in ranked]
    out["selected_top"] = [(c, cls, round(p, 4)) for c, cls, p, g, r2 in ranked]

    # ================= h(c;N) at every N (primary W_FINE + verifier W_WIDE) =================
    hmat = {}   # hmat[c] = {N: {"fine":p, "wide":p, "resolves":.., "r2":..}}
    for c in sel:
        cf = float(Fraction(c))
        hmat[c] = {}
        for N in N_MAIN + [N_VER]:
            rs = data[N][2]
            pf, gf, r2f, resf = h_of(rs, cf, W_FINE, N, 20)
            pw, gw, r2w, resw = h_of(rs, cf, W_WIDE, N, 20)
            hmat[c][N] = {"fine": (round(pf, 5) if pf is not None else None),
                          "wide": (round(pw, 5) if pw is not None else None),
                          "r2": (round(r2f, 4) if r2f is not None else None),
                          "resolves": bool(resf)}
    out["h_matrix"] = hmat

    # ================= CONVERGENCE TEST (card): stable within 0.01 across 2e6..1e7 =================
    conv = {}
    n_stable = 0
    for c in sel:
        vals = []
        for N in (2_000_000, 5_000_000, 10_000_000):
            v = hmat[c][N]["fine"]
            if v is not None:
                vals.append(v)
        if len(vals) < 3:
            conv[c] = {"n_vals": len(vals), "spread": None, "stable0.01": None}
            continue
        spread = max(vals) - min(vals)
        stable = spread <= 0.01
        n_stable += 1 if stable else 0
        conv[c] = {"h_2e6": hmat[c][2_000_000]["fine"], "h_5e6": hmat[c][5_000_000]["fine"],
                   "h_1e7": hmat[c][10_000_000]["fine"], "spread": round(spread, 5),
                   "stable0.01": bool(stable)}
    out["convergence"] = conv
    out["convergence"]["n_stable_of_top%d" % TOPK] = n_stable
    out["convergence"]["limit_exists_per_card"] = bool(n_stable >= 1)

    # ================= KILL TEST (card): top-3, N=1e6 vs N=2e6, all >0.05? =================
    kill = {}
    alldiff = True
    for c in sel[:3]:
        a = hmat[c][1_000_000]["fine"]; b = hmat[c][2_000_000]["fine"]
        d = abs(a - b) if (a is not None and b is not None) else None
        kill[c] = {"h_1e6": a, "h_2e6": b, "diff": (round(d, 5) if d is not None else None)}
        if d is None or d <= 0.05:
            alldiff = False
    out["kill_test_top3"] = kill
    out["kill_test_top3"]["all_diff_gt_0.05"] = bool(alldiff)
    out["kill_test_top3"]["drift_negative_triggered"] = bool(alldiff)

    # ================= TREND: increment shrinkage (convergence vs drift) =================
    def seq(c):
        return [hmat[c][N]["fine"] for N in N_MAIN]
    trend = {}
    for c in sel:
        s = seq(c)
        if any(v is None for v in s):
            trend[c] = None; continue
        d1 = s[1]-s[0]; d2 = s[2]-s[1]; d3 = s[3]-s[2]
        mags = [abs(d1), abs(d2), abs(d3)]
        # shrinking => mags[2]<mags[1]<mags[0]; oscillation => sign changes in d
        signs = [np.sign(d) for d in (d1, d2, d3)]
        osc = int((signs[0]*signs[1] < 0) + (signs[1]*signs[2] < 0))
        shrinking = bool(mags[2] < mags[1] and mags[1] < mags[0])
        trend[c] = {"d12": round(d1, 5), "d25": round(d2, 5), "d5_10": round(d3, 5),
                    "mags": [round(m, 5) for m in mags], "sign_changes": int(osc),
                    "monotone_shrinking": shrinking}
    out["trend"] = trend
    shr = sum(1 for t in trend.values() if isinstance(t, dict) and t["monotone_shrinking"])
    osc = sum(1 for t in trend.values() if isinstance(t, dict) and t["sign_changes"] >= 2)
    out["trend_summary"] = {"n_shrinking": int(shr), "n_oscillating": int(osc),
                            "n_points": len([c for c in sel if isinstance(trend.get(c), dict)])}

    # ================= VERIFIER (i): different windowing W_FINE vs W_WIDE at each N =================
    vw = {}
    allw = []
    for c in sel:
        for N in N_MAIN + [N_VER]:
            a = hmat[c][N]["fine"]; b = hmat[c][N]["wide"]
            if a is not None and b is not None:
                allw.append(abs(a-b))
    out["verifier_windowing"] = {
        "n_pairs": len(allw),
        "mean_abs_diff": round(float(np.mean(allw)), 5) if allw else None,
        "max_abs_diff": round(float(np.max(allw)), 5) if allw else None,
        "frac_within_0.01": round(float(np.mean([d <= 0.01 for d in allw])), 4) if allw else None,
        "frac_within_0.05": round(float(np.mean([d <= 0.05 for d in allw])), 4) if allw else None,
    }
    # per-point windowing diff at N=1e7 (cleanest)
    out["verifier_windowing"]["per_point_1e7"] = {
        c: (round(abs(hmat[c][10_000_000]["fine"] - hmat[c][10_000_000]["wide"]), 5)
            if hmat[c][10_000_000]["fine"] is not None and hmat[c][10_000_000]["wide"] is not None else None)
        for c in sel}

    # ================= VERIFIER (ii): co-sieved N=8e6 vs N=1e7 (adjacent large N) =================
    va = []
    for c in sel:
        a = hmat[c][10_000_000]["fine"]; b = hmat[c][N_VER]["fine"]
        if a is not None and b is not None:
            va.append(abs(a-b))
    out["verifier_cosieved_8e6"] = {
        "n": len(va),
        "mean_abs_diff": round(float(np.mean(va)), 5) if va else None,
        "max_abs_diff": round(float(np.max(va)), 5) if va else None,
        "frac_within_0.01": round(float(np.mean([d <= 0.01 for d in va])), 4) if va else None,
        "frac_within_0.05": round(float(np.mean([d <= 0.05 for d in va])), 4) if va else None,
    }

    # ================= HEADLINE: overall N-to-N variation =================
    allv = []
    for c in sel:
        for N in (2_000_000, 5_000_000, 10_000_000):
            v = hmat[c][N]["fine"]
            if v is not None:
                allv.append(v)
    # per-point mean |delta| across consecutive N in 2e6..1e7
    dsum = 0.0; dn = 0
    for c in sel:
        s = [hmat[c][N]["fine"] for N in (2_000_000, 5_000_000, 10_000_000)]
        if all(v is not None for v in s):
            dsum += abs(s[1]-s[0]); dsum += abs(s[2]-s[1]); dn += 2
    out["headline"] = {
        "n_points": len(sel),
        "points_stable_within_0.01": n_stable,
        "mean_abs_consecutive_delta_2e6_1e7": round(dsum/dn, 5) if dn else None,
        "core_sd_job38": 0.0147,
        "job38_N5e6_1e7_mean_dgamma": 0.052,
    }

    def _pyo(o):
        return o.item() if hasattr(o, "item") else str(o)
    with open("localdim_Nrobust.json", "w") as f:
        json.dump(out, f, indent=2, default=_pyo)

    # ================= console =================
    g = out["gate"]
    print("=== GATE (reproduce job-38 @N=1e7) ===")
    print(" mean_phiN=%s M12=%s growth12=%s gamma12=%s r2=%s" % (g["mean_phiN"], g["M12"], g["growth12"], g["gamma12"], g["r2_12"]))
    print(" 60-core: n=%d mean=%s sd=%s range=[%s,%s]  EXPECT n=60 mean=0.0804 sd=0.0147 [0.0659,0.1177]" %
          (g["core_n"], g["core_mean"], g["core_sd"], g["core_min"], g["core_max"]))
    print(" ALL GATES PASS:", g["all_pass"], " (n_clean=%d)" % g["n_clean_1e7"])
    print("\n=== selected top-%d sharpest (by p @1e7) ===" % TOPK)
    for c, cls, p in out["selected_top"]:
        print("   %-22s %-6s p(1e7)=%s" % (c, cls, p))
    print("\n=== h(c;N) matrix (primary W_FINE) ===")
    hdr = "c".ljust(22) + "  ".join(("N=%s" % ("%.0e" % N)).ljust(10) for N in N_MAIN) + "  8e6"
    print(hdr)
    for c in sel:
        row = c.ljust(22) + "  ".join(("%s" % hmat[c][N]["fine"]).ljust(10) for N in N_MAIN)
        row += "%s" % hmat[c][N_VER]["fine"]
        print(row)
    print("\n=== CONVERGENCE (card): stable within 0.01 across 2e6..1e7? ===")
    for c, v in conv.items():
        if isinstance(v, dict) and v.get("spread") is not None:
            print("   %-22s spread=%s stable0.01=%s (2e6/5e6/1e7 = %s/%s/%s)" %
                  (c, v["spread"], v["stable0.01"], v["h_2e6"], v["h_5e6"], v["h_1e7"]))
    print(" points stable within 0.01: %d / %d  => limit exists per card: %s" %
          (n_stable, TOPK, out["convergence"]["limit_exists_per_card"]))
    print("\n=== KILL TEST (card): top-3, N=1e6 vs N=2e6, all >0.05? ===")
    for c, v in kill.items():
        if isinstance(v, dict):
            print("   %-22s h(1e6)=%s h(2e6)=%s diff=%s" % (c, v["h_1e6"], v["h_2e6"], v["diff"]))
    print(" all top-3 diff>0.05:", out["kill_test_top3"]["all_diff_gt_0.05"],
          " => drift negative triggered:", out["kill_test_top3"]["drift_negative_triggered"])
    print("\n=== TREND (convergence=shrinking increments, drift=oscillating/persistent) ===")
    print(" n_monotone_shrinking=%d  n_oscillating(>=2 sign changes)=%d" % (shr, osc))
    for c, t in trend.items():
        if not isinstance(t, dict):
            continue
        print("   %-22s d=[%s,%s,%s] mags=%s signch=%s shrink=%s" %
              (c, t["d12"], t["d25"], t["d5_10"], t["mags"], t["sign_changes"], t["monotone_shrinking"]))
    print("\n=== VERIFIER (i) different windowing W_FINE vs W_WIDE ===")
    print("  ", out["verifier_windowing"])
    print("=== VERIFIER (ii) co-sieved N=8e6 vs N=1e7 ===")
    print("  ", out["verifier_cosieved_8e6"])
    print("\n=== HEADLINE ===")
    print("  ", out["headline"])
    print("\nWROTE localdim_Nrobust.json")

if __name__ == "__main__":
    main()
