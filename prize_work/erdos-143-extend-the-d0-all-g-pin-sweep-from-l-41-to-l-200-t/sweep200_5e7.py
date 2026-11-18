#!/usr/bin/env python3
"""
Job lab_erdos_261112_35 - RUN B: full sweep L=41..200 at N=5e7 (numpy build_W).

Headline run: extends the job-261111_23 pin (d0^all = 1/4 for L<=40 at N=5e7)
to L<=200 at the SAME N, using a vectorized build_W that must first pass:

  KILL-TEST / GATE KP  build_W_np vs pure-Python build_W at N=200k,
      L in {4,41,200}, all entries < 1e-12.  (numpy path dead -> abort,
       report partial; the card's kill criterion.)
  GATE G3  (my 5e7 pipeline vs job-261111_23 referee-passed values):
      U'(30,7) == 0.15369943043010403  (job-23 per_L["30"]["max_Up"])
      g*(4,1), g*(8,2) == job-23 exact_rows values  (full precision from
      prize_work/ppl-014/d0all_harden.json)

SWEEP: L=41..200, k=1..(L-1)//4 (4k<L):  U'(L,k), U(L,k), L'(L,k).
  - all U' <= 1/6  -> pin extends to L<=200 (combined with job-23 L<=40)
  - any U' > 1/6   -> that L is NOT ruled out (report which, pin to L<=min-1)
  - any L' > 1/6   -> PROVEN crossing below 1/4 (g* >= L' > 1/6): a win.

Single-threaded.  Output: sweep200_5e7.json
"""
import sys, time, json, math, array

HERE = r"E:\Projects\MillenniumResearch\prize_work\erdos-143-extend-the-d0-all-g-pin-sweep-from-l-41-to-l-200-t"
sys.path.insert(0, r"E:\Projects\MillenniumResearch\prize_work\ppl-014")
sys.path.insert(0, HERE)

import numpy as np
from density_threshold import (omega_spf, build_W, gstar, THRESH, RMAX)
from nondeg_threshold import excl_set
from bounds200 import (U_maxlevel, U_sum, L_prime, kmax_below_quarter,
                       build_W_np)

N = 50_000_000
NGATE_KP = 200_000
UP307_REF = 0.15369943043010403   # job-23 per_L["30"]["max_Up"] @5e7
LS_SWEEP = list(range(41, 201))
KILL_PER_L_S = 30.0               # card kill: cheap-eval budget
JOB23 = r"E:\Projects\MillenniumResearch\prize_work\ppl-014\d0all_harden.json"


def wt_np(nmax, exclude):
    t = np.arange(2, nmax + 1, dtype=np.float64)
    w = np.zeros(nmax + 1, dtype=np.float64)
    w[2:] = 1.0 / (t * np.log(t))
    for e in exclude:
        w[e] = 0.0
    return w


def main():
    t0 = time.time()
    res = {"job": "lab_erdos_261112_35_runB", "N": N, "exclude": list(excl_set(10)),
           "threshold": THRESH, "rmax": RMAX}
    print(f"[start] runB N={N} t=0.0s", flush=True)

    # ---------- KILL-TEST: numpy build_W gate at 200k ----------
    om_g = omega_spf(NGATE_KP)
    wt_g = np.zeros(NGATE_KP + 1)
    t_ = np.arange(2, NGATE_KP + 1, dtype=np.float64)
    wt_g[2:] = 1.0 / (t_ * np.log(t_))
    for e in (2, 3, 4):
        wt_g[e] = 0.0
    n_g = np.arange(NGATE_KP + 1, dtype=np.int32)
    om_g_np = np.frombuffer(om_g, dtype=np.int8)
    kp = {}
    kp_ok = True
    for L in (4, 41, 200):
        W_py = build_W(om_g, wt_g.tolist(), NGATE_KP, L, RMAX)  # shared ref impl
        W_np = build_W_np(om_g_np, wt_g, n_g, L, RMAX)
        d = float(np.max(np.abs(np.array(W_py) - W_np)))
        kp[str(L)] = d
        kp_ok = kp_ok and d < 1e-12
    res["killtest_np_vs_py_200k"] = kp
    res["killtest_ok"] = bool(kp_ok)
    print(f"[KP] np vs py max diff at 200k: {kp} ok={kp_ok}  "
          f"{time.time()-t0:.1f}s", flush=True)
    if not kp_ok:
        res["verdict"] = ("KILL: numpy build_W fails 200k gate vs pure-Python "
                          "reference; 5e7 sweep aborted. See run A (5e6).")
        res["elapsed_s"] = round(time.time() - t0, 1)
        with open(HERE + r"\sweep200_5e7.json", "w", encoding="utf-8") as fh:
            json.dump(res, fh, indent=1)
        print("[ABORT] numpy gate failed", flush=True)
        return

    # ---------- sieve + weights at 5e7 ----------
    omega = omega_spf(N)
    wt = wt_np(N, excl_set(10))
    n_arr = np.arange(N + 1, dtype=np.int32)
    om_np = np.frombuffer(omega, dtype=np.int8)
    del omega
    print(f"[sieve] omega+wt done  {time.time()-t0:.1f}s", flush=True)

    # ---------- GATE G3: 5e7 pipeline vs job-23 referee-passed values ----------
    with open(JOB23, "r", encoding="utf-8") as fh:
        job23 = json.load(fh)
    # full-precision job-23 refs (referee-passed), single source of truth:
    pin = {(r["L"], r["k"]): r["gstar"]
           for r in job23["bound_vs_exact_L16"]["pin_rows_k_over_L_is_1_4"]}
    g41_ref, g82_ref = pin[(4, 1)], pin[(8, 2)]
    g3_parts = {}
    W30 = build_W_np(om_np, wt, n_arr, 30, RMAX)
    up307 = float(U_maxlevel(W30, 30, 7, RMAX))
    g3_parts["Up307"] = {"value": up307, "ref": UP307_REF,
                         "ok": abs(up307 - UP307_REF) < 1e-12}
    W4 = build_W_np(om_np, wt, n_arr, 4, RMAX)
    W8 = build_W_np(om_np, wt, n_arr, 8, RMAX)
    g41 = float(gstar(W4, 4, 1, RMAX))
    g82 = float(gstar(W8, 8, 2, RMAX))
    g3_parts["gstar_4_1"] = {"value": g41, "ref": g41_ref,
                             "ok": abs(g41 - g41_ref) < 1e-12}
    g3_parts["gstar_8_2"] = {"value": g82, "ref": g82_ref,
                             "ok": abs(g82 - g82_ref) < 1e-12}
    g3_parts["pin_crosses"] = bool(g41 > THRESH and g82 > THRESH)
    res["gate_g3_5e7"] = g3_parts
    res["gates_ok"] = bool(g3_parts["Up307"]["ok"] and g3_parts["gstar_4_1"]["ok"]
                           and g3_parts["gstar_8_2"]["ok"]
                           and g3_parts["pin_crosses"])
    print(f"[G3] U'(30,7)={up307} ok={g3_parts['Up307']['ok']}; "
          f"g*(4,1)={g41} ok={g3_parts['gstar_4_1']['ok']}; "
          f"g*(8,2)={g82} ok={g3_parts['gstar_8_2']['ok']}  "
          f"{time.time()-t0:.1f}s", flush=True)

    # ---------- SWEEP L=41..200 ----------
    per_L = {}
    unruled, proven_crossings = [], []
    build_times = []
    kill = False
    for i, L in enumerate(LS_SWEEP):
        tb = time.time()
        W = build_W_np(om_np, wt, n_arr, L, RMAX)
        dt = time.time() - tb
        build_times.append(dt)
        if i == 4:   # probe after 5 L's
            mt5 = sum(build_times) / len(build_times)
            res["build_W_np_probe_s"] = mt5
            print(f"[probe] build_W_np mean {mt5:.2f}s/L -> projected "
                  f"{mt5 * len(LS_SWEEP) / 60.0:.1f} min for 160 L", flush=True)
            if mt5 > KILL_PER_L_S:
                kill = True
                break
        kmax = kmax_below_quarter(L)
        rowk = []
        for k in range(1, kmax + 1):
            rowk.append({"k": k, "Up": float(U_maxlevel(W, L, k, RMAX)),
                         "U": float(U_sum(W, L, k, RMAX)),
                         "Lp": float(L_prime(W, L, k, RMAX))})
        maxUp = max(rowk, key=lambda r: r["Up"])
        maxLp = max(rowk, key=lambda r: r["Lp"])
        per_L[str(L)] = {"kmax": kmax, "rows": rowk,
                         "max_Up": maxUp["Up"], "k_at_maxUp": maxUp["k"],
                         "max_Lp": maxLp["Lp"], "k_at_maxLp": maxLp["k"],
                         "margin": THRESH - maxUp["Up"]}
        if maxUp["Up"] > THRESH:
            unruled.append({"L": L, "k": maxUp["k"], "Up": maxUp["Up"]})
        if maxLp["Lp"] > THRESH:
            proven_crossings.append({"L": L, "k": maxLp["k"], "Lp": maxLp["Lp"],
                                     "d": maxLp["k"] / L})
        if L % 40 == 1:
            print(f"  [L={L:3d}] maxUp={maxUp['Up']:.6f} "
                  f"margin={THRESH - maxUp['Up']:.6f}  {time.time()-t0:7.1f}s",
                  flush=True)
    res["kill_triggered"] = bool(kill)
    res["n_L_done"] = len(per_L)
    res["build_W_np_mean_s"] = sum(build_times) / len(build_times) if build_times else None

    if not kill:
        worst_L = min(per_L, key=lambda L: per_L[L]["margin"])
        res.update({
            "per_L": per_L,
            "unruled_Up_over_1_6": unruled,
            "proven_crossings_Lp_over_1_6": proven_crossings,
            "worst_margin_L": worst_L,
            "worst_margin": per_L[worst_L]["margin"],
            "near_zero_margin_lt_0_03": sorted(
                [int(L) for L in per_L if per_L[L]["margin"] < 0.03]),
        })
        if proven_crossings:
            res["verdict"] = (
                f"PROVEN CROSSING BELOW 1/4: L'(L,k)>1/6 (lower bound on g*) at "
                f"{proven_crossings} -> g*(L,k) > 1/6 for k/L < 1/4, so the "
                f"pin d0^all drops below 1/4.  needs_review.")
        elif unruled:
            ruled = [int(L) for L in per_L
                     if not any(u["L"] == int(L) for u in unruled)]
            res["verdict"] = (
                f"PARTIAL: U' > 1/6 (bound too loose, NOT a proven crossing) at "
                f"{unruled}; pin 1/4 verified by U' to L<={max(ruled) if ruled else 40} "
                f"plus job-23 L<=40; the listed L's need exact g* or an LP "
                f"(game-value) sharpening.")
        else:
            res["verdict"] = (
                f"EXTENSION: U'(L,k) <= 1/6 for every (L,k) with 4k<L, L=41..200 "
                f"(worst margin {res['worst_margin']:.6f} at L={worst_L}); "
                f"combined with job-261111_23 (L<=40, exact L<=16), the "
                f"d0^all = 1/4 pin extends to L<=200 at N=5e7.")
    else:
        res["verdict"] = (
            "KILL: numpy build_W > 30 s/L at 5e7 (probe); 5e7 sweep incomplete. "
            f"See run A (N=5e6, pure Python) + probe mean {res['build_W_np_probe_s']:.1f}s/L.")
    res["elapsed_s"] = round(time.time() - t0, 1)
    out = HERE + r"\sweep200_5e7.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n[done] gates_ok={res['gates_ok']} kill={kill} n_L={res['n_L_done']} "
          f"unruled={len(unruled)} crossings={len(proven_crossings)} "
          f"elapsed={res['elapsed_s']}s", flush=True)
    print("VERDICT:", res["verdict"], flush=True)


if __name__ == "__main__":
    main()
