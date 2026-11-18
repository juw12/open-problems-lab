#!/usr/bin/env python3
"""
Job lab_erdos_261112_38 - INDEPENDENT VERIFIER (N=5e7) for the L=201..2000 extension.

Independent re-derivation (a second implementation, not the sweep code re-run):
  O1  OMEGA element-wise: Method 1 (omega_spf, SPF sieve) vs Method 3
      (omega_prime_power, prime-power sieve) at the full N=5e7 -> 0 mismatches.
      (vectorised via np.frombuffer)
  O2  W-table: for representative L (tightest-margin L from the sweep, the low
      edge L=201, the high edge L=2000) an INDEPENDENT pure-Python double-loop
      W2 (Method-3 Omega, no shared build code) compared entry-wise to the
      numpy pipeline W (build_W_np, Method-1 Omega): < 1e-12.
  O3  U'(L,k) and L'(L,k) re-derived with FRESH implementations (U_prime_indep,
      L_prime_indep) on W2, at the sweep's argmax k, compared to the sweep's
      stored per-L max_Up / max_Lp: < 1e-12.
  O4  pin rows g*(4,1), g*(8,2) re-derived by fresh enumeration (Method-3
      Omega) vs the job-23 referee-passed references (d0all_harden.json).
  O5  wt table spot check against the literal formula 1/(n log n).

All ok -> the L=201..2000 headline (U'<=1/6 for every (L,k) with 4k<L, N=5e7)
machinery is independently confirmed.  Single-threaded, < 20 min.
Output: verify2000.json
"""
import sys, time, json, math
HERE = r"E:\Projects\MillenniumResearch\prize_work\erdos-143-extend-the-d0-all-g-pin-sweep-from-l-41-to-l-200-t"
sys.path.insert(0, r"E:\Projects\MillenniumResearch\prize_work\ppl-014")
sys.path.insert(0, HERE)
import numpy as np
from density_threshold import omega_spf, THRESH, RMAX
from nondeg_threshold_verify import omega_prime_power
from bounds200 import kmax_below_quarter, build_W_np

N = 50_000_000
JOB23 = r"E:\Projects\MillenniumResearch\prize_work\ppl-014\d0all_harden.json"
SWEEP = HERE + r"\sweep2000_5e7.json"


def U_prime_indep(W, L, k, rmax):
    best = [(max(W[j][r] for r in range(1, rmax + 1)), j) for j in range(1, L)]
    best.sort(key=lambda t: (t[0], t[1]))
    Sstar = [j for _, j in best[:k]]
    out = 0.0
    for r in range(1, rmax + 1):
        s = sum(W[j][r] for j in Sstar)
        if s > out:
            out = s
    return out


def L_prime_indep(W, L, k, rmax):
    import heapq
    out = 0.0
    for r in range(1, rmax + 1):
        s = sum(heapq.nsmallest(k, (W[j][r] for j in range(1, L))))
        if s > out:
            out = s
    return out


def gstar_indep(W, L, k, rmax):
    from itertools import combinations
    worst = float("inf")
    for S in combinations(range(1, L), k):
        bl = 0.0
        for r in range(1, rmax + 1):
            s = 0.0
            for j in S:
                s += W[j][r]
            if s > bl:
                bl = s
        if bl < worst:
            worst = bl
    return worst


def w2_direct(om3, wt, N, L, rmax):
    W = [[0.0] * (rmax + 1) for _ in range(L)]
    for n in range(2, N + 1):
        r = om3[n]
        if 1 <= r <= rmax:
            W[n % L][r] += wt[n]
    return W


def main():
    t0 = time.time()
    res = {"job": "lab_erdos_261112_38_verify", "N": N, "threshold": THRESH, "rmax": RMAX}
    print(f"[start] verify N={N} t=0.0s", flush=True)

    # ---- O1: element-wise Omega, two independent methods (vectorised) ----
    om1 = omega_spf(N); print(f"[m1] SPF omega {time.time()-t0:.1f}s", flush=True)
    om3 = omega_prime_power(N); print(f"[m3] prime-power omega {time.time()-t0:.1f}s", flush=True)
    a = np.frombuffer(om1, dtype=np.int8); b = np.frombuffer(om3, dtype=np.int8)
    mm = int((a != b).sum()); mx = int(np.abs(a.astype(np.int16) - b.astype(np.int16)).max())
    res["omega_mismatches"] = mm; res["omega_max_absdiff"] = mx; res["omega_agree"] = (mm == 0)
    print(f"[O1] omega mismatches={mm} (max {mx})  {'AGREE' if mm == 0 else 'DISAGREE'}  {time.time()-t0:.1f}s", flush=True)

    # ---- weights + O5 spot check ----
    t_ = np.arange(2, N + 1, dtype=np.float64)
    wt = np.zeros(N + 1, dtype=np.float64); wt[2:] = 1.0 / (t_ * np.log(t_))
    for e in (2, 3, 4): wt[e] = 0.0
    n_arr = np.arange(N + 1, dtype=np.int32)
    om1_np = np.frombuffer(om1, dtype=np.int8)
    spot = {}
    for n in (1000003, 49999999, 123457, 2, 4):
        ref = 0.0 if n in (2, 3, 4) else 1.0 / (n * math.log(n))
        spot[n] = bool(abs(float(wt[n]) - ref) < 1e-15)
    res["wt_spot_ok"] = all(spot.values()) and bool(spot)
    print(f"[O5] wt spot {spot}  {time.time()-t0:.1f}s", flush=True)

    # ---- sweep + job-23 data ----
    with open(SWEEP, "r", encoding="utf-8") as fh:
        runB = json.load(fh)
    with open(JOB23, "r", encoding="utf-8") as fh:
        job23 = json.load(fh)
    pin = {(r["L"], r["k"]): r["gstar"] for r in job23["bound_vs_exact_L16"]["pin_rows_k_over_L_is_1_4"]}
    N7 = runB.get("N_5e7") or runB.get("N_5e6") or {}
    perL = N7.get("per_L", {})
    tightest = str(N7.get("worst_margin_L", 2000))
    L_check = sorted({201, 2000, int(tightest)})
    res["N_used_for_perL"] = "N_5e7" if runB.get("N_5e7") else "N_5e6"
    res["L_check"] = L_check
    print(f"[L_check] {L_check}  (tightest={tightest})  {time.time()-t0:.1f}s", flush=True)

    # ---- O2 + O3 per checked L ----
    o2, o3 = {}, {}; all_o2 = all_o3 = True
    for L in L_check:
        tb = time.time()
        W2 = w2_direct(om3, wt, N, L, RMAX)
        Wn = build_W_np(om1_np, wt, n_arr, L, RMAX)
        dW = max(abs(W2[j][r] - float(Wn[j, r])) for j in range(L) for r in range(RMAX + 1))
        o2[str(L)] = dW; all_o2 = all_o2 and dW < 1e-12
        st = perL.get(str(L))
        if st is not None:
            dUp = abs(U_prime_indep(W2, L, st["k_at_maxUp"], RMAX) - st["max_Up"])
            dLp = abs(L_prime_indep(W2, L, st["k_at_maxLp"], RMAX) - st["max_Lp"])
            o3[str(L)] = {"dUp": dUp, "dLp": dLp,
                          "maxUp_indep": U_prime_indep(W2, L, st["k_at_maxUp"], RMAX),
                          "maxLp_indep": L_prime_indep(W2, L, st["k_at_maxLp"], RMAX),
                          "stored_maxUp": st["max_Up"], "stored_maxLp": st["max_Lp"]}
        else:
            dUp = dLp = float("nan")
            o3[str(L)] = {"note": f"L={L} not in sweep per_L (partial?)"}
        all_o3 = all_o3 and (dW < 1e-12) and (st is not None) and (dUp < 1e-12) and (dLp < 1e-12)
        print(f"[O2/O3] L={L} dW={dW:.2e} dUp={dUp:.2e} dLp={dLp:.2e}  {time.time()-tb:.1f}s  {time.time()-t0:.1f}s", flush=True)
    res["o2_W_entry_maxdiff"] = o2
    res["o3_bounds_maxdiff"] = o3

    # ---- O4: pin rows, fresh enumeration, Method-3 Omega ----
    W4 = w2_direct(om3, wt, N, 4, RMAX); W8 = w2_direct(om3, wt, N, 8, RMAX)
    g41 = gstar_indep(W4, 4, 1, RMAX); g82 = gstar_indep(W8, 8, 2, RMAX)
    o4 = {
        "gstar_4_1": {"value": g41, "ref_job23": pin[(4, 1)], "ok": abs(g41 - pin[(4, 1)]) < 1e-12},
        "gstar_8_2": {"value": g82, "ref_job23": pin[(8, 2)], "ok": abs(g82 - pin[(8, 2)]) < 1e-12},
        "pin_crosses_1_6": bool(g41 > THRESH and g82 > THRESH),
    }
    res["o4_pin_rows"] = o4
    print(f"[O4] g*(4,1)={g41} g*(8,2)={g82} ok={o4['gstar_4_1']['ok'] and o4['gstar_8_2']['ok']}  {time.time()-t0:.1f}s", flush=True)

    # ---- headline cross-check ----
    gates_ok = (runB.get("gates") or {}).get("gates_ok") if isinstance(runB.get("gates"), dict) else runB.get("gates_ok")
    res["runB_status"] = runB.get("status")
    res["runB_verdict"] = runB.get("verdict")
    res["runB_gates_ok"] = gates_ok
    res["all_ok"] = bool(res["omega_agree"] and res["wt_spot_ok"] and all_o2 and all_o3
                         and o4["gstar_4_1"]["ok"] and o4["gstar_8_2"]["ok"]
                         and o4["pin_crosses_1_6"] and gates_ok)
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(HERE + r"\verify2000.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(f"\n[done] all_ok={res['all_ok']} omega_agree={res['omega_agree']} O2_ok={all_o2} O3_ok={all_o3} "
          f"O4_ok={o4['gstar_4_1']['ok'] and o4['gstar_8_2']['ok']} gates_ok={gates_ok} elapsed={res['elapsed_s']}s", flush=True)


if __name__ == "__main__":
    main()
