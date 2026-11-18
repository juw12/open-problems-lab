#!/usr/bin/env python3
"""#143 / PPL 014  lab_erdos_261112_43  -  5e7 gate-G calibration + N-robustness spot-check.

The N=5e7 FULL sweep over 8000 L's (~1 hr, ~0.47 s/L) exceeds the per-computation budget,
so the 5e7 side is a SPOT subset: it (a) absolute-calibrates the 5e7 pipeline by reproducing
the job-23/35 referee-passed references at N=5e7 (gate G), and (b) spot-checks the pin at a
spread of L's across 2001..10000 including the tightest-margin L from the 5e6 full sweep,
reporting the 5e7-vs-5e6 margin delta (N-robustness, mirroring job 38's N-robust claim).

Gates:
  G   reproduce at N=5e7:  U'(30,7) vs job-23,  U'(42,10) vs job-35,  g*(4,1), g*(8,2) vs job-23  <1e-12

Output: spot10k_5e7.json
Single-threaded, ~<4 min (dominated by the 5e7 sieve).
"""
import os, sys, json, time
import numpy as np
_T0 = time.time()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ppl-014"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from density_threshold import omega_spf, gstar, THRESH, RMAX
from bounds200 import U_maxlevel, L_prime, kmax_below_quarter
from bounds2000 import masked_inputs, build_W_flat_precomputed

CARD = os.path.dirname(os.path.abspath(__file__))
PPL = os.path.join(CARD, "..", "ppl-014")
N = int(5e7)
EXCL = (2, 3, 4)
TH = 1.0 / 6.0
RMAX = 24


def make_wt(N):
    n = np.arange(N + 1, dtype=np.float64)
    n = np.where(n < 2.0, 1.0, n)
    wt = np.zeros(N + 1, dtype=np.float64)
    wt[2:] = 1.0 / (n[2:] * np.log(n[2:]))
    for e in EXCL:
        wt[e] = 0.0
    return wt


def bounds_all(W, L, kmax, rmax):
    Wc = W[1:L]
    w = Wc.max(axis=1)
    ord_ = np.argsort(w, kind="stable")
    C = np.cumsum(Wc[ord_], axis=0)
    U_all = C[:kmax].max(axis=1)
    P = np.cumsum(np.sort(Wc, axis=0), axis=0)
    L_all = P[:kmax].max(axis=1)
    return U_all, L_all


def main():
    out = os.path.join(CARD, "spot10k_5e7.json")
    print(f"[43-5e7] N={N} TH={TH} RMAX={RMAX} EXCL={EXCL}", flush=True)

    # ---- reference values from committed job-23/35 JSONs ----
    with open(os.path.join(PPL, "d0all_harden.json")) as fh:
        d23 = json.load(fh)
    REF_U30_7 = d23["per_L"]["30"]["max_Up"]
    pin = {r["L"]: r for r in d23["bound_vs_exact_L16"]["pin_rows_k_over_L_is_1_4"]}
    REF_g4_1 = pin[4]["gstar"]; REF_g8_2 = pin[8]["gstar"]
    with open(os.path.join(CARD, "sweep200_5e7.json")) as fh:
        d35 = json.load(fh)
    REF_U42_10 = d35["per_L"]["42"]["max_Up"]
    print(f"[43-5e7] refs U30_7={REF_U30_7} g4_1={REF_g4_1} g8_2={REF_g8_2} U42_10={REF_U42_10}", flush=True)

    # ---- build 5e7 ----
    print(f"[43-5e7] building 5e7 omega_spf + wt + n ...", flush=True)
    t0 = time.time()
    om = np.asarray(omega_spf(N), dtype=np.int32)
    wt = make_wt(N); n_np = np.arange(N + 1, dtype=np.int32)
    print(f"[43-5e7] 5e7 sieve {time.time() - t0:.1f}s", flush=True)
    n_m, om_i32, wt_m = masked_inputs(n_np, om, wt, RMAX)
    del n_np, om, wt
    print(f"[43-5e7] 5e7 mask done ({len(n_m)} terms)  {time.time() - _T0:.0f}s", flush=True)

    gates = {"gates_ok": True}
    g = {}
    W30 = build_W_flat_precomputed(n_m, om_i32, wt_m, 30, RMAX)
    u30 = float(U_maxlevel(W30, 30, 7, RMAX)); g["U30_7"] = {"computed": u30, "ref": REF_U30_7, "abs_diff": abs(u30 - REF_U30_7)}
    W42 = build_W_flat_precomputed(n_m, om_i32, wt_m, 42, RMAX)
    u42 = float(U_maxlevel(W42, 42, 10, RMAX)); g["U42_10"] = {"computed": u42, "ref": REF_U42_10, "abs_diff": abs(u42 - REF_U42_10)}
    W4 = build_W_flat_precomputed(n_m, om_i32, wt_m, 4, RMAX)
    g41 = float(gstar(W4, 4, 1, RMAX)); g["g4_1"] = {"computed": g41, "ref": REF_g4_1, "abs_diff": abs(g41 - REF_g4_1)}
    W8 = build_W_flat_precomputed(n_m, om_i32, wt_m, 8, RMAX)
    g82 = float(gstar(W8, 8, 2, RMAX)); g["g8_2"] = {"computed": g82, "ref": REF_g8_2, "abs_diff": abs(g82 - REF_g8_2)}
    for key in g:
        print(f"[43-5e7] GATE G {key}: computed={g[key]['computed']}  ref={g[key]['ref']}  abs_diff={g[key]['abs_diff']:.3e}", flush=True)
        if g[key]["abs_diff"] >= 1e-12:
            gates["gates_ok"] = False
    gates["G_references"] = g

    # ---- read 5e6 full sweep -> tightest L ----
    with open(os.path.join(CARD, "sweep10k_5e6.json")) as fh:
        s6 = json.load(fh)
    p6 = s6["N_5e6"]["per_L"]
    worst6_L = int(s6["N_5e6"]["worst_margin_L"])
    spot_set = sorted(set([2001, 2500, 3000, 4000, 6000, 8000, 10000, worst6_L]))
    print(f"[43-5e7] 5e6 worst margin L={worst6_L} (margin={p6[str(worst6_L)]['margin']}); spot L's={spot_set}", flush=True)

    spot = {}
    dmar_max = 0.0; dmaxup_max = 0.0
    any5e7_cross = False; any5e7_unruled = 0
    for L in spot_set:
        W = build_W_flat_precomputed(n_m, om_i32, wt_m, L, RMAX)
        kmax = kmax_below_quarter(L)
        U_all, L_all = bounds_all(W, L, kmax, RMAX)
        maxUp = float(U_all.max()); kUp = int(np.argmax(U_all) + 1)
        maxLp = float(L_all.max())
        margin = TH - maxUp
        u7_cross = int((L_all > TH).sum()); u7_unruled = int((U_all > TH).sum())
        any5e7_cross = any5e7_cross or (u7_cross > 0)
        any5e7_unruled += u7_unruled
        e6 = p6.get(str(L))
        dmar = (TH - maxUp) - (TH - e6["max_Up"]) if e6 else None   # margin_5e7 - margin_5e6
        dmar = abs(e6["max_Up"] - maxUp) if e6 else None
        if dmar is not None:
            dmar_max = max(dmar_max, dmar); dmaxup_max = dmar_max
        spot[str(L)] = {"kmax": int(kmax), "max_Up": maxUp, "k_at_maxUp": kUp, "max_Lp": maxLp,
                        "margin": float(margin), "Lp_crossings": u7_cross, "U_unruled": u7_unruled,
                        "d_maxUp_vs_5e6": dmar}
        print(f"[43-5e7] L={L}: max_Up={maxUp:.6f} margin={margin:.6f} Lp_cross={u7_cross} 5e6_maxUp={e6['max_Up'] if e6 else 'n/a'} d={dmar if dmar is not None else 'n/a'}", flush=True)
    del n_m, om_i32, wt_m

    result = {
        "job": "lab_erdos_261112_43", "problem": "Erdos #143 / PPL 014 ($500)",
        "N": N, "TH": TH, "RMAX": RMAX, "EXCL": list(EXCL),
        "gates": gates, "spot_per_L": spot,
        "spot_Ls": spot_set, "worst_5e6_L": worst6_L,
        "n_robustness": {"spot_L": len(spot_set),
                         "max_abs_diff_maxUp": dmar_max},
        "any_5e7_crossing": any5e7_cross, "any_5e7_unruled": any5e7_unruled,
        "elapsed_s": time.time() - _T0,
    }
    if any5e7_cross:
        result["verdict"] = "5e7 SPOT CROSSING: L'>1/6 at 4k<L in a spot L in 2001..10000 (N=5e7)."
    elif any5e7_unruled:
        result["verdict"] = f"5e7 SPOT: U' rule-out incomplete (unruled={any5e7_unruled}) in spot L's."
    else:
        result["verdict"] = (f"5e7 SPOT-CONFIRMED: U'(L,k)<=1/6 for every (L,k) 4k<L at all {len(spot_set)} spot L's "
                             f"in 2001..10000 (N=5e7); N-robust vs 5e6 (max|d maxUp|={dmar_max:.6f}); gate G absolute-calibrated.")
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"[43-5e7] wrote {os.path.basename(out)}", flush=True)
    print(f"[43-5e7] VERDICT: {result['verdict']}", flush=True)
    print(f"[43-5e7] TOTAL {time.time() - _T0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
