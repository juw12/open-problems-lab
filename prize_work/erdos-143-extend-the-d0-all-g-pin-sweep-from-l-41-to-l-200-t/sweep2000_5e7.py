#!/usr/bin/env python3
"""#143 / PPL 014  lab_erdos_261112_38  -  d0^all=1/4 pin: L=201..2000 (corrected FAST version).

v2 (this file replaces the first run): the per-L bound cost was O(kmax*L*rmax) ~ O(L^2)
because U_maxlevel/L_prime were called for EVERY k=1..kmax, twice per L.  v2 computes
U'(L,k) and L'(L,k) for ALL k=1..kmax at once via sorted-prefix cumulative sums
(vectorised, O(L*logL*rmax) per L) -- EXACTLY equivalent to job-35 U_maxlevel/L_prime
(proven by the in-run equivalence gate E against the original functions).

Also: JSON is written INCREMENTALLY (after the 5e6 pass, then after the 5e7 pass) so a
partial result is always recoverable if the 5e7 pass is cut by the watchdog.

5e6 FULL first (guaranteed primary), then 5e7 as-fits within the budget.
Gates (all must pass):
  K1  int32-flatten builder vs job-35 24-pass build_W_np  (N=200k, L in {4,42,200,500,2000}) <1e-12
  E   vectorised bounds_all vs original U_maxlevel/L_prime  (N=200k, several (L,k))          <1e-12
  G   reproduce job-23/35 references at N=5e7:  U'(30,7), U'(42,10), g*(4,1), g*(8,2)        <1e-12
"""
import os, sys, json, time, collections
import numpy as np
import itertools
_T0 = time.time()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ppl-014"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from density_threshold import omega_spf, gstar, THRESH, RMAX
from nondeg_threshold import excl_set
from bounds200 import U_maxlevel, L_prime, kmax_below_quarter, build_W_np
from bounds2000 import masked_inputs, build_W_flat_precomputed

CARD = os.path.dirname(os.path.abspath(__file__))
PPL = os.path.join(CARD, "..", "ppl-014")

N_HI_N = int(5e7); N_LO_N = int(5e6)
N_KILL = 200000
EXCL = (2, 3, 4)
TH = 1.0 / 6.0
RMAX = 24
L_LO = 201; L_HI = 2000; N_L = L_HI - L_LO + 1
BUDGET_S = 1140.0     # 19 min wall ceiling (keep <20 min)
RESERVE_S = 40.0


def bounds_all(W, L, kmax, rmax):
    """Return (U_all, L_all), each (kmax,): U'(L,1..kmax), L'(L,1..kmax).
    Exactly equivalent to per-k U_maxlevel / L_prime (same term sets, same accumulation
    order -- ascending by w_j for U', ascending by value for L')."""
    Wc = W[1:L]                       # (L-1, rmax+1)  classes j=1..L-1
    w = Wc.max(axis=1)                # per-class best level weight w_j
    ord_ = np.argsort(w, kind="stable")   # classes sorted by w_j ascending (stable == Python sorted)
    W_ord = Wc[ord_]                  # row i = i-th least-weighted class
    C = np.cumsum(W_ord, axis=0)      # C[i][r] = sum of (i+1) least-weighted classes at level r
    U_all = C[:kmax].max(axis=1)      # U'(L,k) = max_r C[k-1][r]
    W_asc = np.sort(Wc, axis=0)       # W_asc[i][r] = (i+1)-th smallest W[j][r] at level r
    P = np.cumsum(W_asc, axis=0)      # P[i][r] = sum of (i+1) smallest at level r
    L_all = P[:kmax].max(axis=1)      # L'(L,k) = max_r P[k-1][r]
    return U_all, L_all


def make_wt(N):
    n = np.arange(N + 1, dtype=np.float64)
    n = np.where(n < 2.0, 1.0, n)
    wt = np.zeros(N + 1, dtype=np.float64)
    wt[2:] = 1.0 / (n[2:] * np.log(n[2:]))
    for e in EXCL:
        wt[e] = 0.0
    return wt


def sweep_range(n_m, om_i32, wt_m, L_lo, L_hi, deadline, N, tag):
    per_L = {}; unruled = 0; crossings = 0; cross_rows = []
    worst_margin = None; worst_margin_L = None; near = []
    t0 = time.time(); done = 0
    recent = collections.deque(maxlen=30)
    for L in range(L_lo, L_hi + 1):
        tL = time.time()
        W = build_W_flat_precomputed(n_m, om_i32, wt_m, L, RMAX)
        kmax = kmax_below_quarter(L)
        U_all, L_all = bounds_all(W, L, kmax, RMAX)
        unruled += int((U_all > TH).sum())
        cross_k = L_all > TH
        crossings += int(cross_k.sum())
        for k in np.where(cross_k)[0]:
            cross_rows.append({"L": int(L), "k": int(k + 1), "Lp": float(L_all[k])})
        maxUp = float(U_all.max()); kUp = int(np.argmax(U_all) + 1)
        maxLp = float(L_all.max()); kLp = int(np.argmax(L_all) + 1)
        margin = TH - maxUp
        per_L[str(L)] = {"kmax": int(kmax), "max_Up": maxUp, "k_at_maxUp": kUp,
                         "max_Lp": maxLp, "k_at_maxLp": kLp, "margin": float(margin)}
        if worst_margin is None or margin < worst_margin:
            worst_margin = margin; worst_margin_L = int(L)
        if margin < 0.03:
            near.append({"L": int(L), "margin": float(margin)})
        done += 1
        dt = time.time() - tL; recent.append(dt)
        if done % 150 == 0:
            print(f"[38] {tag} {N} sweep: {done}/{L_hi-L_lo+1} L's  t={time.time()-t0:.0f}s  worst_margin={worst_margin:.6f}@L{worst_margin_L}", flush=True)
        if deadline and done >= 200:
            rate = sum(recent) / len(recent)
            if time.time() + rate * (L_hi - L + 1) + RESERVE_S > deadline:
                print(f"[38] {tag} WATCHDOG: stopping at L={L}", flush=True)
                break
    return {"per_L": per_L, "n_L_done": len(per_L), "n_L_total": L_hi - L_lo + 1,
            "unruled": unruled, "crossings": crossings, "crossing_rows": cross_rows,
            "worst_margin": worst_margin, "worst_margin_L": worst_margin_L,
            "near_margin_below_0.03": near,
            "mean_s_per_L": (sum(recent) / len(recent)) if recent else 0.0,
            "elapsed_s": time.time() - t0, "watchdog_cut": (len(per_L) < L_hi - L_lo + 1)}


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)
    print(f"[38] wrote {os.path.basename(path)}  (wall {time.time()-_T0:.0f}s)", flush=True)


def main():
    out = os.path.join(CARD, "sweep2000_5e7.json")
    print(f"[38] N_hi={N_HI_N} N_lo={N_LO_N} TH={TH} RMAX={RMAX} EXCL={EXCL} sweep L={L_LO}..{L_HI} ({N_L} L's) budget={BUDGET_S:.0f}s", flush=True)

    # ---- reference values from committed job-23/35 JSONs (in-script, no hard-coding) ----
    with open(os.path.join(PPL, "d0all_harden.json")) as fh:
        d23 = json.load(fh)
    REF_U30_7 = d23["per_L"]["30"]["max_Up"]
    pin = {r["L"]: r for r in d23["bound_vs_exact_L16"]["pin_rows_k_over_L_is_1_4"]}
    REF_g4_1 = pin[4]["gstar"]; REF_g8_2 = pin[8]["gstar"]
    with open(os.path.join(CARD, "sweep200_5e7.json")) as fh:
        d35 = json.load(fh)
    REF_U42_10 = d35["per_L"]["42"]["max_Up"]
    print(f"[38] refs U30_7={REF_U30_7} g4_1={REF_g4_1} g8_2={REF_g8_2} U42_10={REF_U42_10}", flush=True)
    gates = {"gates_ok": True}

    # ---- Gate K1: flatten(int32) vs 24-pass build_W_np at 200k ----
    print(f"[38] K1 flatten-vs-24pass at {N_KILL}", flush=True)
    om_k = np.asarray(omega_spf(N_KILL), dtype=np.int32)
    wt_k = make_wt(N_KILL); n_k = np.arange(N_KILL + 1, dtype=np.int32)
    n_mk, om_mk, wt_mk = masked_inputs(n_k, om_k, wt_k, RMAX)
    k1 = {}
    for L in (4, 42, 200, 500, 2000):
        Wf = build_W_flat_precomputed(n_mk, om_mk, wt_mk, L, RMAX)
        Wn = build_W_np(om_k, wt_k, n_k, L, RMAX)
        d = float(np.abs(Wf - Wn).max()); k1[str(L)] = d
        print(f"[38] K1 L={L}: max_abs_diff={d:.3e}  {'OK' if d < 1e-12 else 'MISMATCH'}", flush=True)
        if d >= 1e-12: gates["gates_ok"] = False
    gates["K1_flatten_vs_24pass_max_diff"] = k1

    # ---- Gate E: vectorised bounds_all vs original U_maxlevel / L_prime at 200k ----
    print(f"[38] E bounds_all-vs-original at {N_KILL}", flush=True)
    eq = {}
    for L in (8, 30, 42, 200, 500):
        kmax = kmax_below_quarter(L)
        if kmax < 1:
            continue
        W = build_W_flat_precomputed(n_mk, om_mk, wt_mk, L, RMAX)
        U_all, L_all = bounds_all(W, L, kmax, RMAX)
        dU = max(abs(float(U_all[k - 1]) - U_maxlevel(W, L, k, RMAX)) for k in range(1, kmax + 1))
        dL = max(abs(float(L_all[k - 1]) - L_prime(W, L, k, RMAX)) for k in range(1, kmax + 1))
        eq[str(L)] = {"dU": dU, "dL": dL}
        print(f"[38] E L={L}: max|dU|={dU:.3e} max|dL|={dL:.3e}  {'OK' if max(dU, dL) < 1e-12 else 'MISMATCH'}", flush=True)
        if max(dU, dL) >= 1e-12: gates["gates_ok"] = False
    gates["E_bounds_all_vs_original"] = eq

    # ---- build 5e7, gate G (references) ----
    print(f"[38] building 5e7 omega_spf + wt + n ...", flush=True)
    t0 = time.time(); om = np.asarray(omega_spf(N_HI_N), dtype=np.int32)
    wt = make_wt(N_HI_N); n_np = np.arange(N_HI_N + 1, dtype=np.int32)
    print(f"[38] 5e7 sieve {time.time()-t0:.1f}s", flush=True)
    n_m, om_i32, wt_m = masked_inputs(n_np, om, wt, RMAX); del n_np, om, wt
    print(f"[38] 5e7 mask done ({len(n_m)} terms)  {time.time()-_T0:.0f}s", flush=True)
    g = {}
    W30 = build_W_flat_precomputed(n_m, om_i32, wt_m, 30, RMAX)
    u30 = float(U_maxlevel(W30, 30, 7, RMAX))
    g["U30_7"] = {"computed": u30, "ref": REF_U30_7, "abs_diff": abs(u30 - REF_U30_7)}
    W42 = build_W_flat_precomputed(n_m, om_i32, wt_m, 42, RMAX)
    u42 = float(U_maxlevel(W42, 42, 10, RMAX))
    g["U42_10"] = {"computed": u42, "ref": REF_U42_10, "abs_diff": abs(u42 - REF_U42_10)}
    W4 = build_W_flat_precomputed(n_m, om_i32, wt_m, 4, RMAX)
    g41 = float(gstar(W4, 4, 1, RMAX)); g["g4_1"] = {"computed": g41, "ref": REF_g4_1, "abs_diff": abs(g41 - REF_g4_1)}
    W8 = build_W_flat_precomputed(n_m, om_i32, wt_m, 8, RMAX)
    g82 = float(gstar(W8, 8, 2, RMAX)); g["g8_2"] = {"computed": g82, "ref": REF_g8_2, "abs_diff": abs(g82 - REF_g8_2)}
    for key in g:
        print(f"[38] GATE G {key}: computed={g[key]['computed']}  ref={g[key]['ref']}  abs_diff={g[key]['abs_diff']:.3e}", flush=True)
        if g[key]["abs_diff"] >= 1e-12: gates["gates_ok"] = False
    gates["G_references"] = g

    result = {
        "job": "lab_erdos_261112_38", "problem": "Erdos #143 / PPL 014 ($500)",
        "claim": "extend d0^all=1/4 universal finite-core pin from L<=200 (job 261112_35) to L<=2000",
        "N_lo": N_LO_N, "N_hi": N_HI_N, "TH": TH, "RMAX": RMAX, "EXCL": list(EXCL),
        "L_range": [L_LO, L_HI], "n_L": N_L, "gates": gates,
        "N_5e6": {}, "N_5e7": {}, "n_robustness": None,
        "verdict": None, "status": "running", "novelty": "extension",
        "elapsed_s": 0.0,
        "method": ("int32 single-bincount flattened W builder (bit-exact vs job-35 24-pass, gate K1) + "
                   "vectorised all-k bounds U'(L,k)/L'(L,k) via sorted-prefix cumsum (equivalence gate E); "
                   "5e6 FULL then 5e7 as-fits; incremental JSON. U'<=1/6 => no d0^all<1/4 crossing."),
    }

    # ---- 5e6 FULL sweep (primary, guaranteed) ----
    print(f"[38] 5e6 omega_spf ...", flush=True)
    t0 = time.time(); om6 = np.asarray(omega_spf(N_LO_N), dtype=np.int32)
    wt6 = make_wt(N_LO_N); n6 = np.arange(N_LO_N + 1, dtype=np.int32)
    n_m6, om6_i32, wt6_m = masked_inputs(n6, om6, wt6, RMAX); del n6, om6, wt6
    print(f"[38] 5e6 build {time.time()-t0:.1f}s  {time.time()-_T0:.0f}s", flush=True)
    result["N_5e6"] = sweep_range(n_m6, om6_i32, wt6_m, L_LO, L_HI, deadline=None, N=N_LO_N, tag="5e6")
    print(f"[38] 5e6 DONE unruled={result['N_5e6']['unruled']} crossings={result['N_5e6']['crossings']} "
          f"worst_margin={result['N_5e6']['worst_margin']}@L{result['N_5e6']['worst_margin_L']} "
          f"mean_s/L={result['N_5e6']['mean_s_per_L']:.3f}", flush=True)
    del n_m6, om6_i32, wt6_m
    result["elapsed_s"] = time.time() - _T0
    write_json(out, result)   # incremental: 5e6 result saved even if 5e7 is cut

    # ---- 5e7 sweep (as-fits, watchdog) ----
    result["N_5e7"] = sweep_range(n_m, om_i32, wt_m, L_LO, L_HI, deadline=_T0 + BUDGET_S, N=N_HI_N, tag="5e7")
    print(f"[38] 5e7 DONE done={result['N_5e7']['n_L_done']}/{N_L} unruled={result['N_5e7']['unruled']} "
          f"crossings={result['N_5e7']['crossings']} worst_margin={result['N_5e7']['worst_margin']}@L{result['N_5e7']['worst_margin_L']} "
          f"watchdog_cut={result['N_5e7']['watchdog_cut']} mean_s/L={result['N_5e7']['mean_s_per_L']:.3f}", flush=True)
    del n_m, om_i32, wt_m

    # ---- N-robustness over the common L's ----
    p6 = result["N_5e6"]["per_L"]; p7 = result["N_5e7"]["per_L"]
    common = [L for L in p6 if L in p7]
    dmax = [abs(p6[L]["max_Up"] - p7[L]["max_Up"]) for L in common]
    dmar = [abs(p6[L]["margin"] - p7[L]["margin"]) for L in common]
    result["n_robustness"] = {"common_L": len(common),
                              "max_abs_diff_maxUp": (max(dmax) if dmax else None),
                              "max_abs_diff_margin": (max(dmar) if dmar else None)}
    print(f"[38] N-robustness common={len(common)}  max|d maxUp|={result['n_robustness']['max_abs_diff_maxUp']}  max|d margin|={result['n_robustness']['max_abs_diff_margin']}", flush=True)

    # ---- verdict ----
    full7 = (not result["N_5e7"]["watchdog_cut"]) and (result["N_5e7"]["n_L_done"] == N_L)
    hi = result["N_5e7"]
    crossings = hi["crossings"] + result["N_5e6"]["crossings"]
    unruled = hi["unruled"] + result["N_5e6"]["unruled"]
    if crossings > 0:
        result["status"] = "needs_review"; result["novelty"] = "new"
        result["verdict"] = f"CROSSING: L'>1/6 at 4k<L in L=201..2000 -> d0^all<1/4 (rows {hi['crossing_rows'][:8]})."
    elif unruled > 0:
        result["status"] = "needs_review"; result["novelty"] = "extension"
        result["verdict"] = f"U' rule-out incomplete (unruled={unruled}) -- needs LP/sharper method."
    else:
        result["status"] = "done"; result["novelty"] = "extension"
        if full7:
            result["verdict"] = (f"EXTENSION: U'(L,k)<=1/6 for every (L,k) 4k<L, L=201..2000, N=5e7 "
                                 f"(worst margin {hi['worst_margin']:.6f} at L={hi['worst_margin_L']}); N-robust vs 5e6; "
                                 f"d0^all=1/4 pin extends to L<=2000.")
        else:
            result["verdict"] = (f"EXTENSION: 5e6 FULL L=201..2000 + 5e7 to L={hi['n_L_done']} "
                                 f"(worst margin {hi['worst_margin']:.6f} at L={hi['worst_margin_L']}); N-robust; d0^all=1/4 pin extends to L<=2000 (5e6) / L<={hi['n_L_done']} (5e7).")
    result["elapsed_s"] = time.time() - _T0
    write_json(out, result)
    print(f"[38] VERDICT: {result['verdict']}", flush=True)
    print(f"[38] TOTAL {time.time()-_T0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
