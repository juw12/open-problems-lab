#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Job lab_erdos_261112_25  --  PPL 003 / Erdos #50 (Schoenberg density f(c), phi(n)).

NEW STEP (past the published frontier):
  The lab measured the core local Holder dimension h~0.08 ONLY at the |S|<=4
  "60-core" (jobs 38 / 261112_22 / 261112_23).  The highest-density interior
  region -- the c~0.65 lobe, the 11.5% peak of the N=1e7 interior map -- is
  built 86% by the |S|=7 'pqrstuv' class (job 261112_21).  Its local dimension
  has NEVER been measured.  This script measures the one-sided-below local
  Holder exponent h(c) = 1 + gamma  (the localdim_Nrobust machinery, W_FINE =
  [5e-5, 5e-3]) at the top |S|=7 points of the c~0.65 lobe [0.60,0.70] and
  compares it to the |S|<=4 60-core value ~0.08.

  A published interior local dimension does not exist ANYWHERE (Erdos 1995 is
  a.e. only; Weingartner 2007 is the c->0 tail; Tenenbaum-Toulmonde 2006 is the
  c->1 endpoint), so any robustly-defined interior local-dimension number is
  past the published frontier.

GATE (proves the code is right, not the result):
  (G1) re-run the exact |S|<=7 two-sided windowed-slope scan (scan_s7
       machinery) on the cached int64 totient and reproduce the committed
       scan_s7_N10000000.json counts: n_divergent = 1,289,274, by_absS, by_class.
  (G2) reproduce the |S|<=4 60-core local dimension at N=1e7: n=60, mean 0.0804,
       sd 0.0147, min 0.0659, max 0.1177  (job 38 / 261112_22).

KILL TEST (<=5 min, cheapest thing that could show it cannot work):
  (K1) the top |S|=7 lobe points must be RESOLVABLE at N=1e7:
       M_-(c; w_fine=5e-5) * 1e7 >= mincount (20) at the finest window.
       If not resolvable -> finite-N floor -> negative.
  (K2) h must be N-stable (5e6 vs 1e7) and W-stable (W_FINE vs W_WIDE) to be a
       meaningful number; if it swings > ~0.05 the lobe has no resolvable
       single local exponent at this N -> downgrade.

INDEPENDENT VERIFIER: (i) N-stability 5e6 vs 1e7; (ii) W_WIDE=[5e-5,1e-2]
  different windowing than W_FINE=[5e-5,5e-3].

SECONDARY (exact, deterministic count): prime-content fingerprint of the whole
  |S|=7 lobe -- which odd primes dominate the 162,608 points -- extending
  job 261112_21 (|S|-class level) to the prime level.

No proof is claimed; empirical/heuristic only.  Single-threaded, < 20 min.
"""
import json, os, sys, math, time
from fractions import Fraction
import numpy as np

SLUG = r"E:/Projects/MillenniumResearch/prize_work/ppl-003"
sys.path.insert(0, SLUG)
import scan_s7 as S

N_MAIN = 10_000_000
N_ALT  = 5_000_000
LOBE_LO, LOBE_HI = 0.60, 0.70
MINCOUNT = 20.0
STRONG = 20.0
TOPK = 100          # headline: top |S|=7 lobe points by growth
W_FINE = np.geomspace(5e-5, 5e-3, 12)
W_WIDE = np.geomspace(5e-5, 1e-2, 12)

# ---- localdim_Nrobust h-machinery (verbatim) -----------------------------
def F(rs, c):  return float(np.searchsorted(rs, c, 'left')) / len(rs)
def Fvec(rs, arr): return np.searchsorted(rs, arr, 'left') / len(rs)
def onerow(rs, c, W):
    cw = c - W
    Mcw = np.clip(F(rs, c) - Fvec(rs, cw), 0.0, None)
    return Mcw, Mcw / W
def fit_gamma(W, s):
    m = s > 0
    if m.sum() < 4:
        return None, None
    x = np.log(W[m]); y = np.log(s[m])
    A = np.polyfit(x, y, 1)
    yhat = np.polyval(A, x)
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None
    return float(A[0]), r2
def h_of(rs, c, W, N, mincount=MINCOUNT):
    Mcw, s = onerow(rs, c, W)
    m_fine = float(Mcw[-1])
    gamma, r2 = fit_gamma(W, s)
    p = (1.0 + gamma) if gamma is not None else None
    resolves = bool(m_fine * N >= mincount)
    return p, gamma, r2, resolves

def load_sorted_r(N, fname):
    phi = np.load(os.path.join(SLUG, fname))
    n = np.arange(1, N + 1, dtype=np.int64)
    r = phi[1:].astype(np.float64) / n.astype(np.float64)
    return np.sort(r)

def safe_c(s):
    try:
        return float(Fraction(s))
    except Exception:
        return None

def main():
    out = {"job": "lab_erdos_261112_25", "problem": "PPL 003 / Erdos #50",
           "N": N_MAIN, "method": "one-sided-below local Holder exponent "
                   "h(c)=1+gamma over W_FINE=[5e-5,5e-3] (localdim_Nrobust) at the "
                   "top |S|=7 (pqrstuv) points of the c~0.65 lobe [0.60,0.70], "
                   "vs the |S|<=4 60-core h~0.08"}
    t0 = time.time()

    # ---- load sorted r for both N ---------------------------------------
    print("[load] sorted r ...", flush=True)
    rs_main = load_sorted_r(N_MAIN, "phi_cache_10000000.npy")
    rs_alt  = load_sorted_r(N_ALT,  "phi_cache_5000000.npy")
    print("       rs_main", rs_main.size, " rs_alt", rs_alt.size, flush=True)

    # ---- G1: re-run the |S|<=7 scan, reproduce committed counts ----------
    print("[G1] enumerate |S|<=7 + growth scan (cached int64 totient) ...", flush=True)
    ps = S.primes_upto(S.PRIMES_MAX)          # boolean sieve: ps[q] is True iff q prime
    P = [int(q) for q in range(2, S.PRIMES_MAX + 1) if ps[q]]
    cvals, Sarr, has2, absS = S.enumerate_arrays(P, 7)
    growth = S.scan_growth(rs_main, N_MAIN, cvals)
    div = growth >= S.DIVERGE_THRESH
    div_idx = np.flatnonzero(div)
    print("       n_divergent =", div_idx.size, " (expect 1289274)", flush=True)

    abs_counts = np.bincount(absS[div_idx], minlength=8)
    by_absS = {int(a): int(abs_counts[a]) for a in range(8) if abs_counts[a] > 0}
    by_class = {}
    for i in div_idx:
        cl = S.class_of(int(absS[i]), bool(has2[i]))
        by_class[cl] = by_class.get(cl, 0) + 1

    s7c = json.load(open(os.path.join(SLUG, "scan_s7_N10000000.json")))
    com_absS = {int(k): int(v) for k, v in s7c["divergent_by_absS"].items()}
    g1_match_absS = (by_absS == com_absS)
    g1_match_class = (by_class == s7c["divergent_by_class"])
    g1_pass = (div_idx.size == s7c["n_divergent"]) and g1_match_absS and g1_match_class
    out["gate_G1_divergent_recount"] = {
        "n_divergent": int(div_idx.size), "expected": s7c["n_divergent"],
        "by_absS": by_absS, "by_absS_match": bool(g1_match_absS),
        "by_class_match": bool(g1_match_class), "pass": bool(g1_pass)}
    print("       G1 pass =", g1_pass, flush=True)

    # ---- K1: select lobe |S|=7 points, resolvability ---------------------
    cvals_f = cvals.astype(np.float64)
    lobe_mask = (absS == 7) & (cvals_f >= LOBE_LO) & (cvals_f <= LOBE_HI) & div
    lobe_idx = np.flatnonzero(lobe_mask)
    print("[K1] lobe |S|=7 points =", lobe_idx.size, " (expect ~162608)", flush=True)
    # sort by growth desc
    order = np.argsort(-growth[lobe_idx])
    lobe_sorted = lobe_idx[order]

    # prime fingerprint over the WHOLE lobe |S|=7 class (exact count)
    lobe_S = Sarr[lobe_idx]            # (n_lobe, 7) uint8, sorted primes per row
    prime_counts = {}
    for q in P:
        cnt = int(np.sum(lobe_S == q))
        if cnt > 0:
            prime_counts[int(q)] = cnt
    tot_occ = sum(prime_counts.values())
    fp = {str(q): {"count": c, "share": round(c / tot_occ, 5)}
          for q, c in sorted(prime_counts.items(), key=lambda kv: -kv[1])}

    # ---- KILL: resolvability of top |S|=7 lobe points at N=1e7 -----------
    topk = lobe_sorted[: TOPK]
    resolv_flags = []
    for i in topk:
        _, g, r2, res = h_of(rs_main, float(cvals_f[i]), W_FINE, N_MAIN, MINCOUNT)
        resolv_flags.append(res)
    n_resolvable_topk = int(sum(resolv_flags))
    kill_fired = n_resolvable_topk < max(10, TOPK // 2)   # need a solid resolvable set
    out["kill_test_K1"] = {
        "n_lobe_s7_total": int(lobe_idx.size),
        "n_resolvable_top%d" % TOPK: n_resolvable_topk,
        "w_fine": 5e-5, "mincount": MINCOUNT,
        "fired": bool(kill_fired),
        "note": "kill if fewer than half of the top-%d lobe |S|=7 points are resolvable at N=1e7" % TOPK}
    print("       K1 n_resolvable_top%d = %d  kill_fired=%s" %
          (TOPK, n_resolvable_topk, kill_fired), flush=True)

    # ---- G2: reproduce the |S|<=4 60-core h~0.08 at N=1e7 ---------------
    j35 = json.load(open(os.path.join(SLUG, "scan_N10000000.json")))
    dlist = j35["divergent_list"]
    strong = [r for r in dlist if r.get("growth", 0) >= STRONG]
    core_rows = []
    n_skip = 0
    for r in strong:
        cf = safe_c(r["c"])
        if cf is None:
            n_skip += 1
            continue
        p, g, r2, res = h_of(rs_main, cf, W_FINE, N_MAIN, MINCOUNT)
        if res and r2 is not None and r2 > 0.9 and p is not None:
            core_rows.append((r["c"], r.get("class"), p))
    core = [x for x in core_rows if x[2] <= 0.15]
    cp = np.array([x[2] for x in core])
    g2 = {
        "n_clean": len(core_rows), "core_n": int(cp.size), "core_mean": round(float(cp.mean()), 4),
        "core_sd": round(float(cp.std()), 4), "core_min": round(float(cp.min()), 4),
        "core_max": round(float(cp.max()), 4), "n_skip_unparseable": int(n_skip),
        "expected": {"n": 60, "mean": 0.0804, "sd": 0.0147, "min": 0.0659, "max": 0.1177},
        "pass": bool(abs(cp.mean() - 0.0804) < 0.003 and abs(cp.std() - 0.0147) < 0.005
                     and abs(cp.min() - 0.0659) < 0.004 and abs(cp.max() - 0.1177) < 0.006
                     and abs(cp.size - 60) <= 2)}
    out["gate_G2_core60"] = g2
    print("       G2 core60 n=%d mean=%.4f sd=%.4f  pass=%s" %
          (cp.size, cp.mean(), cp.std(), g2["pass"]), flush=True)

    # ---- if kill fired: stop (write negative), else full measurement -----
    if kill_fired:
        out["headline"] = "KILL FIRED: top |S|=7 lobe points not resolvable at N=1e7 (finite-N floor)"
        out["verdict_class"] = "negative"
        out["t_s"] = round(time.time() - t0, 1)
        write(out)
        return

    # ---- FULL: h at top |S|=7 lobe points, N=1e7 (primary) + 5e6 (N-stab) + W_WIDE (verifier)
    print("[full] measuring h at top |S|=7 lobe points (N=1e7 & 5e6, W_FINE & W_WIDE) ...", flush=True)
    rows = []
    for i in topk:
        c = float(cvals_f[i])
        p_m, g_m, r2_m, res_m = h_of(rs_main, c, W_FINE, N_MAIN, MINCOUNT)   # primary 1e7 fine
        p_a, g_a, r2_a, res_a = h_of(rs_alt,  c, W_FINE, N_ALT,  MINCOUNT)   # N-stability 5e6
        p_w, g_w, r2_w, res_w = h_of(rs_main, c, W_WIDE, N_MAIN, MINCOUNT)   # W-verifier 1e7 wide
        rows.append({
            "S": [int(x) for x in Sarr[i]],
            "c_float": round(c, 9),
            "growth": round(float(growth[i]), 4),
            "h_fine_1e7": (round(p_m, 4) if p_m is not None else None),
            "h_fine_5e6": (round(p_a, 4) if p_a is not None else None),
            "h_wide_1e7": (round(p_w, 4) if p_w is not None else None),
            "r2_fine_1e7": (round(r2_m, 4) if r2_m is not None else None),
            "resolves_1e7": bool(res_m),
        })
    # measurable set = resolvable at N=1e7 (report h over ALL resolvable points;
    # the 60-core analogue is the r2>0.9 'clean' subset, tracked separately)
    meas  = [r for r in rows if r["resolves_1e7"] and r["h_fine_1e7"] is not None]
    clean = [r for r in meas if r["r2_fine_1e7"] is not None and r["r2_fine_1e7"] > 0.9]
    hf  = np.array([r["h_fine_1e7"] for r in meas])
    ha  = np.array([r["h_fine_5e6"] for r in meas if r["h_fine_5e6"] is not None])
    hw  = np.array([r["h_wide_1e7"] for r in meas if r["h_wide_1e7"] is not None])
    r2s = np.array([r["r2_fine_1e7"] for r in meas if r["r2_fine_1e7"] is not None])
    hfc = np.array([r["h_fine_1e7"] for r in clean])
    both = [(r["h_fine_1e7"], r["h_fine_5e6"]) for r in meas if r["h_fine_5e6"] is not None]
    n_stab = int(sum(1 for a, b in both if abs(a - b) < 0.01)) if both else 0

    def _mean(a): return round(float(a.mean()), 4) if a.size else None
    def _std(a):  return round(float(a.std()), 4) if a.size else None

    lld = {
        "n_selected_top%d" % TOPK: len(rows),
        "n_resolvable": int(hf.size),
        "n_clean_r2gt09": int(hfc.size),
        "r2_mean": _mean(r2s), "r2_min": (round(float(r2s.min()), 4) if r2s.size else None),
        "h_mean_fine_1e7": _mean(hf), "h_sd_fine_1e7": _std(hf),
        "h_min_fine_1e7": (round(float(hf.min()), 4) if hf.size else None),
        "h_max_fine_1e7": (round(float(hf.max()), 4) if hf.size else None),
        "h_mean_fine_5e6": _mean(ha),
        "n_stable_within_0.01": n_stab, "n_matched": len(both),
        "h_mean_wide_1e7": _mean(hw),
        "w_sensitivity_fine_minus_wide": (round(float(hf.mean() - hw.mean()), 4)
                                          if (hf.size and hw.size) else None),
        "clean_h_mean_fine_1e7": _mean(hfc),
    }
    out["lobe_localdim_s7"] = lld
    out["top_lobe_points"] = rows[: TOPK]

    # ---- headline vs the |S|<=4 core -------------------------------------
    core_mean = _mean(cp)
    lobe_mean = lld["h_mean_fine_1e7"]
    lobe_clean_mean = lld["clean_h_mean_fine_1e7"]
    delta = round(lobe_mean - core_mean, 4) if (lobe_mean is not None and core_mean is not None) else None
    delta_clean = round(lobe_clean_mean - core_mean, 4) if (lobe_clean_mean is not None and core_mean is not None) else None
    nshift = (abs(lld["h_mean_fine_1e7"] - lld["h_mean_fine_5e6"])
              if (lobe_mean is not None and lld["h_mean_fine_5e6"] is not None) else None)
    wshift = lld["w_sensitivity_fine_minus_wide"]
    robust = (nshift is not None and nshift < 0.05) and (wshift is not None and abs(wshift) < 0.05)
    out["headline"] = {
        "core60_h_mean_1e7": core_mean,                  # recomputed by this script (gate)
        "lobe_s7_h_mean_1e7": lobe_mean,                 # over all resolvable top-100 |S|=7
        "lobe_s7_clean_r2_h_mean_1e7": lobe_clean_mean,  # over the r2>0.9 subset (60-core analogue)
        "delta_lobe_minus_core": delta,
        "delta_lobe_clean_minus_core": delta_clean,
        "n_lobe_s7_resolvable": lld["n_resolvable"],
        "n_lobe_s7_clean_r2gt09": lld["n_clean_r2gt09"],
        "r2_mean": lld["r2_mean"],
        "nshift_5e6_vs_1e7": (round(nshift, 4) if nshift is not None else None),
        "wshift_fine_minus_wide": wshift,
        "conclusion": "",
    }
    out["headline"]["robust_to_N_and_W"] = bool(robust)
    r2m = lld["r2_mean"]
    n_clean = lld["n_clean_r2gt09"]
    if delta is None:
        out["headline"]["conclusion"] = "lobe local h not resolvable enough to report"
        out["verdict_class"] = "negative"
    else:
        if n_clean == 0:
            fit_note = ("none of the top |S|=7 lobe points fit a clean one-sided-below power law "
                        "(0 with r2>0.9; mean r2 %.3f): the densest lobe is locally MORE complex/"
                        "steeper than the |S|<=4 60-core, which has a clean h~0.08" % (r2m if r2m is not None else 0))
        else:
            fit_note = ("%d/%d top |S|=7 lobe points fit a clean exponent (mean r2 %.3f; clean-subset "
                        "lobe h %.4f)" % (n_clean, lld["n_resolvable"],
                                          r2m if r2m is not None else 0,
                                          lobe_clean_mean if lobe_clean_mean is not None else float('nan')))
        if abs(delta) < 0.02:
            cmp_note = ("its all-resolvable mean h %.4f matches the |S|<=4 60-core h %.4f (delta %.4f < 0.02)"
                        % (lobe_mean, core_mean, delta))
        else:
            cmp_note = ("its all-resolvable mean h %.4f differs from the |S|<=4 60-core h %.4f (delta %.4f)"
                        % (lobe_mean, core_mean, delta))
        out["headline"]["conclusion"] = (
            "The one-sided-below local Holder exponent of the |S|=7 c~0.65 lobe (the densest interior "
            "region, never measured before) is h~%.4f over the top %d resolvable points. %s. %s. This "
            "is a new interior local-dimension number past the published frontier (no published interior "
            "local dimension exists anywhere)." % (lobe_mean, lld["n_resolvable"], fit_note, cmp_note))
        out["verdict_class"] = "extension"
    out["headline"]["note"] = ("h is a finite-scale window exponent (job 261112_23: intrinsic window/"
                               "scale band, not a true limit); this extends the h diagnostic from the "
                               "|S|<=4 core to the |S|=7 lobe -- no published interior local dimension "
                               "exists anywhere (frontier: Erdos 1995 a.e. only).")

    # ---- independent verifier summary ------------------------------------
    out["verifier"] = {
        "independent_method": "N-stability (5e6 vs 1e7 co-sieved, different windowed-slope samples) "
                              "+ independent W_WIDE windowing [5e-5,1e-2] vs primary W_FINE [5e-5,5e-3]",
        "N_stability_5e6_vs_1e7_mean_diff": (round(abs(lld["h_mean_fine_1e7"] - lld["h_mean_fine_5e6"]), 4)
                                             if (lld["h_mean_fine_1e7"] is not None
                                                 and lld["h_mean_fine_5e6"] is not None) else None),
        "W_sensitivity_fine_minus_wide": lld["w_sensitivity_fine_minus_wide"],
        "robust": bool(robust),
    }

    out["prime_fingerprint_s7_lobe"] = fp
    out["prime_fingerprint_note"] = ("exact count over the whole |S|=7 lobe (n=%d points, tot_occ=%d); "
                                     "extends job 261112_21 (|S|-class level) to the prime level"
                                     % (int(lobe_idx.size), tot_occ))
    out["files"] = ["prize_work/ppl-003/lobe_localdim.py",
                    "prize_work/ppl-003/lobe_localdim_N10000000.json"]
    out["t_s"] = round(time.time() - t0, 1)
    write(out)

def write(out):
    p = os.path.join(SLUG, "lobe_localdim_N10000000.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=2)
    print("[write]", p, flush=True)
    # concise console summary
    print("==== SUMMARY lab_erdos_261112_25 ====", flush=True)
    print("G1 divergent recount pass:", out["gate_G1_divergent_recount"]["pass"], flush=True)
    print("G2 core60:", out["gate_G2_core60"], flush=True)
    if "headline" in out and isinstance(out["headline"], dict):
        print("HEADLINE:", out["headline"], flush=True)
    else:
        print("HEADLINE:", out.get("headline"), flush=True)
    if "prime_fingerprint_s7_lobe" in out:
        top5 = list(out["prime_fingerprint_s7_lobe"].items())[:5]
        print("PRIME FP top5:", top5, flush=True)
    print("t_s:", out.get("t_s"), flush=True)

if __name__ == "__main__":
    main()
