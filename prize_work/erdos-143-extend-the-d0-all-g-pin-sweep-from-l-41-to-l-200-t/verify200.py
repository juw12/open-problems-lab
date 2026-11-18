#!/usr/bin/env python3
"""
Job lab_erdos_261112_35 - INDEPENDENT VERIFIER (N=5e7).

Independent part (per the card's verifier spec: "element-wise Omega
cross-check (0 mismatches), re-derived by an independent script"):
  O1  OMEGA: Method 1 (omega_spf, SPF sieve) vs Method 3
      (omega_prime_power, prime-power sieve from job-37's verifier) checked
      ELEMENT-WISE at the full N=5e7 -> 0 mismatches expected (job-23's
      verifier achieved 0 at this N, run 183 s).
  O2  W-table: for L_check = {41, 200, 30, tightest-L from run B} an
      INDEPENDENT pure-Python double-loop W2 (no shared build_W code,
      Method-3 Omega) compared entry-wise to the numpy pipeline's W
      (Method-1 Omega): < 1e-12.  This catches level/level-assignment bugs
      beyond the pointwise Omega check.
  O3  U' and L' re-derived with FRESH implementations (U_prime_indep,
      L_prime_indep) on W2 for every k <= (L-1)//4 at each checked L;
      compared to run B's stored per-L max_Up / max_Lp (< 1e-12).
  O4  pin rows g*(4,1), g*(8,2) re-derived by fresh enumeration (Method-3
      Omega) vs run B's gate values AND vs the job-23 referee-passed
      references (d0all_harden.json bound_vs_exact_L16.pin_rows).
  O5  wt table spot check against the literal formula 1/(n log n).

All ok -> run B's headline (U' <= 1/6 for all (L,k) with 4k<L, L=41..200,
N=5e7) is independently confirmed.

Single-threaded, < 20 min.  Output: verify200.json
"""
import sys, time, json, math, array, heapq

HERE = r"E:\Projects\MillenniumResearch\prize_work\erdos-143-extend-the-d0-all-g-pin-sweep-from-l-41-to-l-200-t"
sys.path.insert(0, r"E:\Projects\MillenniumResearch\prize_work\ppl-014")
sys.path.insert(0, HERE)

import numpy as np
from density_threshold import omega_spf, THRESH, RMAX
from nondeg_threshold_verify import omega_prime_power
from bounds200 import kmax_below_quarter, build_W_np

N = 50_000_000
JOB23 = r"E:\Projects\MillenniumResearch\prize_work\ppl-014\d0all_harden.json"


def U_prime_indep(W, L, k, rmax):
    """Fresh re-implementation of the worst-class upper bound U'(L,k)."""
    best = [(max(W[j][r] for r in range(1, rmax + 1)), j)
            for j in range(1, L)]
    best.sort(key=lambda t: (t[0], t[1]))
    Sstar = [j for _, j in best[:k]]
    out = 0.0
    for r in range(1, rmax + 1):
        s = sum(W[j][r] for j in Sstar)
        if s > out:
            out = s
    return out


def L_prime_indep(W, L, k, rmax):
    """Fresh re-implementation of the lower bound L'(L,k)."""
    out = 0.0
    for r in range(1, rmax + 1):
        s = sum(heapq.nsmallest(k, (W[j][r] for j in range(1, L))))
        if s > out:
            out = s
    return out


def gstar_indep(W, L, k, rmax):
    """Fresh exact min-S max-r enumeration (combinations, no shared code)."""
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
    """Independent W table: plain Python double loop, Method-3 Omega."""
    W = [[0.0] * (rmax + 1) for _ in range(L)]
    for n in range(2, N + 1):
        r = om3[n]
        if 1 <= r <= rmax:
            W[n % L][r] += wt[n]
    return W


def main():
    t0 = time.time()
    res = {"job": "lab_erdos_261112_35_verify", "N": N,
           "threshold": THRESH, "rmax": RMAX}
    print(f"[start] verify N={N} t=0.0s", flush=True)

    # ---------- O1: element-wise Omega, two independent methods ----------
    om1 = omega_spf(N)
    print(f"[m1] SPF omega done  {time.time()-t0:.1f}s", flush=True)
    om3 = omega_prime_power(N)
    print(f"[m3] prime-power omega done  {time.time()-t0:.1f}s", flush=True)
    mm, mx = 0, 0.0
    for x in range(N + 1):
        if om1[x] != om3[x]:
            mm += 1
            d = abs(int(om1[x]) - int(om3[x]))
            if d > mx:
                mx = d
    res["omega_mismatches"] = mm
    res["omega_max_absdiff"] = mx
    res["omega_agree"] = (mm == 0)
    print(f"[O1] omega mismatches={mm} (max {mx})  "
          f"{'AGREE' if mm == 0 else 'DISAGREE'}  {time.time()-t0:.1f}s",
          flush=True)

    # ---------- weights (numpy, same formula) + O5 spot check ----------
    t_ = np.arange(2, N + 1, dtype=np.float64)
    wt = np.zeros(N + 1, dtype=np.float64)
    wt[2:] = 1.0 / (t_ * np.log(t_))
    for e in (2, 3, 4):
        wt[e] = 0.0
    n_arr = np.arange(N + 1, dtype=np.int32)
    om1_np = np.frombuffer(om1, dtype=np.int8)
    spot = {}
    for n in (1000003, 49999999, 123457, 2, 4):
        ref = 0.0 if n in (2, 3, 4) else 1.0 / (n * math.log(n))
        spot[n] = bool(abs(float(wt[n]) - ref) < 1e-15)
    res["wt_spot_ok"] = all(spot.values()) and bool(spot)
    print(f"[O5] wt spot check {spot}  {time.time()-t0:.1f}s", flush=True)

    # ---------- run B + job-23 data ----------
    with open(HERE + r"\sweep200_5e7.json", "r", encoding="utf-8") as fh:
        runB = json.load(fh)
    with open(JOB23, "r", encoding="utf-8") as fh:
        job23 = json.load(fh)
    pin = {(r["L"], r["k"]): r["gstar"]
           for r in job23["bound_vs_exact_L16"]["pin_rows_k_over_L_is_1_4"]}
    perL = runB["per_L"]
    tightest = str(runB["worst_margin_L"])
    L_check = sorted({41, 200, 30, int(tightest)})
    res["L_check"] = L_check

    # ---------- O2 + O3 per checked L ----------
    o2, o3 = {}, {}
    all_o2 = all_o3 = True
    for L in L_check:
        tb = time.time()
        W2 = w2_direct(om3, wt, N, L, RMAX)
        Wn = build_W_np(om1_np, wt, n_arr, L, RMAX)
        dW = max(abs(W2[j][r] - float(Wn[j, r]))
                 for j in range(L) for r in range(RMAX + 1))
        o2[str(L)] = dW
        all_o2 = all_o2 and dW < 1e-12
        kmax = kmax_below_quarter(L)
        up_rows, lp_rows = {}, {}
        for k in range(1, kmax + 1):
            up_rows[k] = U_prime_indep(W2, L, k, RMAX)
            lp_rows[k] = L_prime_indep(W2, L, k, RMAX)
        st = perL.get(str(L))
        if st is not None:
            dUp = max(abs(up_rows[k] - row["Up"])
                      for k, row in enumerate(st["rows"], start=1))
            dLp = max(abs(lp_rows[k] - row["Lp"])
                      for k, row in enumerate(st["rows"], start=1))
        else:  # L=30: run B's ref is the G3 gate value U'(30,7) (= max over k)
            ref_up = runB["gate_g3_5e7"]["Up307"]["value"]
            dUp = abs(max(up_rows.values()) - ref_up)
            dLp = 0.0
        o3[str(L)] = {"dUp_max": dUp, "dLp_max": dLp,
                      "maxUp_indep": max(up_rows.values()),
                      "maxLp_indep": max(lp_rows.values())}
        all_o3 = all_o3 and dUp < 1e-12 and dLp < 1e-12
        print(f"[O2/O3] L={L} dW={dW:.2e} dUp={dUp:.2e} dLp={dLp:.2e} "
              f"  {time.time()-t0:.1f}s", flush=True)
    res["o2_W_entry_maxdiff"] = o2
    res["o3_bounds_maxdiff"] = o3

    # ---------- O4: pin rows, fresh enumeration, Method-3 Omega ----------
    W4 = w2_direct(om3, wt, N, 4, RMAX)
    W8 = w2_direct(om3, wt, N, 8, RMAX)
    g41 = gstar_indep(W4, 4, 1, RMAX)
    g82 = gstar_indep(W8, 8, 2, RMAX)
    o4 = {
        "gstar_4_1": {"value": g41, "ref_job23": pin[(4, 1)],
                      "ref_runB": runB["gate_g3_5e7"]["gstar_4_1"]["value"],
                      "ok": abs(g41 - pin[(4, 1)]) < 1e-12
                      and abs(g41 - runB["gate_g3_5e7"]["gstar_4_1"]["value"]) < 1e-12},
        "gstar_8_2": {"value": g82, "ref_job23": pin[(8, 2)],
                      "ref_runB": runB["gate_g3_5e7"]["gstar_8_2"]["value"],
                      "ok": abs(g82 - pin[(8, 2)]) < 1e-12
                      and abs(g82 - runB["gate_g3_5e7"]["gstar_8_2"]["value"]) < 1e-12},
        "pin_crosses_1_6": bool(g41 > THRESH and g82 > THRESH),
    }
    res["o4_pin_rows"] = o4
    print(f"[O4] g*(4,1)={g41} g*(8,2)={g82} ok={o4['gstar_4_1']['ok']}"
          f" and {o4['gstar_8_2']['ok']}  {time.time()-t0:.1f}s", flush=True)

    # ---------- headline cross-check ----------
    res["runB_gates_ok"] = runB.get("gates_ok")
    res["runB_verdict"] = runB.get("verdict")
    res["all_ok"] = bool(res["omega_agree"] and res["wt_spot_ok"]
                         and all_o2 and all_o3 and o4["gstar_4_1"]["ok"]
                         and o4["gstar_8_2"]["ok"] and o4["pin_crosses_1_6"]
                         and runB.get("gates_ok")
                         and not runB.get("kill_triggered", True)
                         and runB.get("n_L_done") == 160
                         and not runB.get("unruled_Up_over_1_6")
                         and not runB.get("proven_crossings_Lp_over_1_6"))
    res["elapsed_s"] = round(time.time() - t0, 1)
    out = HERE + r"\verify200.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1,
                  default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(f"\n[done] all_ok={res['all_ok']} omega_agree={res['omega_agree']} "
          f"O2_ok={all_o2} O3_ok={all_o3} O4_ok={o4['gstar_4_1']['ok'] and o4['gstar_8_2']['ok']} "
          f"elapsed={res['elapsed_s']}s", flush=True)


if __name__ == "__main__":
    main()
