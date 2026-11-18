#!/usr/bin/env python3
"""
Independent verifier for d0all_harden.py (job lab_erdos_261111_23).

Independent part = a SECOND Omega implementation (Method 3, prime-power sieve),
cross-checked ELEMENT-WISE against the SPF sieve (Method 1) at the FULL N=5e7.
That is the load-bearing check: the main headline is a deterministic function of
Omega (via build_W -> g*/U'), so 0 element-wise Omega mismatches validates the
whole headline by transitivity, exactly the pattern job 37's verifier used.

On top of that it RE-DERIVES (with the independent Method-3 Omega):
  V1  the pin: g*(4,1) > 1/6 and g*(8,2) > 1/6  ->  d0^all = 1/4 at L<=16;
  V2  the crossing rule-out at the tightest L (L=30) and at the range endpoints
      (L=17, L=40):  U'(L,k) <= 1/6 for every k with k/L < 1/4;
  V3  the bound is EXACT at k=1 (U' == g*, gap < 1e-12) and a VALID upper bound
      (g* <= U') at k>=2, at L<=16.

It re-implements U' independently (U_maxlevel_indep) and compares it to the
main script's U_maxlevel to check the bound code, not just the Omega.

Single-threaded, N=5e7, <20 min.
"""
import sys, time, json, math, array

sys.path.insert(0, r"E:\Projects\MillenniumResearch\prize_work\ppl-014")

from density_threshold import omega_spf, build_W, gstar, THRESH, RMAX
from nondeg_threshold import excl_set
from nondeg_threshold_verify import omega_prime_power
from d0all_harden import (weight_table_arr, per_class_best, U_sum,
                          U_maxlevel)

N = 50_000_000
L_EXACT = [3, 4, 8]              # pin rows (L=4,k=1 ; L=8,k=2) + L=3 (no lower xing)
L_ANALYTIC = [17, 24, 30, 35, 40]  # endpoints + tightest (L=30)


def U_maxlevel_indep(W, L, k, rmax):
    """Independent re-implementation of the worst-class upper bound U'(L,k)."""
    wm = []
    for j in range(1, L):
        m = 0.0
        for r in range(1, rmax + 1):
            if W[j][r] > m:
                m = W[j][r]
        wm.append((m, j))
    wm.sort()
    Sstar = [j for m, j in wm[:k]]
    best = 0.0
    for r in range(1, rmax + 1):
        s = 0.0
        for j in Sstar:
            s += W[j][r]
        if s > best:
            best = s
    return best


def compare_omega(a, b, nmax):
    mm = 0
    mx = 0.0
    for x in range(nmax + 1):
        if a[x] != b[x]:
            mm += 1
            d = abs(a[x] - b[x])
            if d > mx:
                mx = d
    return mm, mx


def main():
    t0 = time.time()
    result = {"job": "lab_erdos_261111_23_verify", "N": N,
              "exclude": list(excl_set(10)), "threshold": THRESH, "rmax": RMAX}
    print(f"[start] verifier N={N}  t=0.0s", flush=True)

    # ---------- KEY: two independent Omega methods, element-wise at 5e7 ----------
    om1 = omega_spf(N)
    print(f"[m1] omega_spf (SPF) done   {time.time()-t0:.1f}s", flush=True)
    om3 = omega_prime_power(N)
    print(f"[m3] omega_prime_power done {time.time()-t0:.1f}s", flush=True)
    mm, mx = compare_omega(om1, om3, N)
    result["omega_mismatches"] = mm
    result["omega_max_absdiff"] = mx
    result["omega_agree"] = (mm == 0)
    print(f"[omega] mismatches={mm} (max {mx})  ->  "
          f"{'AGREE' if mm == 0 else 'DISAGREE'}   {time.time()-t0:.1f}s",
          flush=True)

    wt = weight_table_arr(N, excl_set(10))

    # ---------- V1 + V3: exact g* at L in L_EXACT ----------
    exact = {}
    for L in L_EXACT:
        W = build_W(om3, wt, N, L, RMAX)          # independent Omega
        for k in range(1, L):
            ge = gstar(W, L, k, RMAX)             # shared combinatorics
            up_main = U_maxlevel(W, L, k, RMAX)   # main bound
            up_indep = U_maxlevel_indep(W, L, k, RMAX)  # independent bound
            exact.setdefault(L, []).append({
                "k": k, "gstar": ge, "Up_main": up_main,
                "Up_indep": up_indep,
                "bound_code_ok": abs(up_main - up_indep) < 1e-12,
                "upper_bound_ok": ge <= up_main + 1e-15,   # g* <= U'
                "gap_k1_exact": (abs(up_main - ge) < 1e-12) if k == 1 else None,
            })
        print(f"  [L={L}] exact g* re-derived (Method-3 omega)  "
              f"{time.time()-t0:.1f}s", flush=True)

    g41 = next(r["gstar"] for r in exact[4] if r["k"] == 1)
    g82 = next(r["gstar"] for r in exact[8] if r["k"] == 2)
    g31 = next(r["gstar"] for r in exact[3] if r["k"] == 1)
    v1_pin = (g41 > THRESH) and (g82 > THRESH)
    v1_no_lower = all(r["gstar"] <= THRESH
                      for L in L_EXACT for r in exact[L]
                      if r["k"] / L < 0.25 - 1e-12)
    v3_upper = all(r["upper_bound_ok"] for L in L_EXACT for r in exact[L])
    v3_k1exact = all(r["gap_k1_exact"] for L in L_EXACT
                     for r in exact[L] if r["k"] == 1)
    v3_code = all(r["bound_code_ok"] for L in L_EXACT for r in exact[L])
    result["v1_pin"] = {
        "gstar_4_1": g41, "gstar_8_2": g82, "gstar_3_1": g31,
        "pin_crosses_1_6": v1_pin, "no_lower_crossing_Lsmall": v1_no_lower,
    }
    result["v3_bound"] = {
        "g_leq_Uall_rows": v3_upper, "U_eq_g_at_k1": v3_k1exact,
        "bound_code_matches_indep": v3_code,
    }

    # ---------- V2: crossing rule-out at L_ANALYTIC (incl. tightest L=30) ----------
    cross = {}
    all_ruled = True
    for L in L_ANALYTIC:
        W = build_W(om3, wt, N, L, RMAX)
        kmax_below = (L - 1) // 4
        worst = None
        for k in range(1, kmax_below + 1):
            up = U_maxlevel_indep(W, L, k, RMAX)   # independent bound
            if worst is None or up > worst["Up"]:
                worst = {"k": k, "Up": up}
        cross[str(L)] = {"kmax_below_1_4": kmax_below, "k_at_maxUp": worst["k"],
                         "max_Up": worst["Up"],
                         "ruled_out": worst["Up"] <= THRESH}
        if worst["Up"] > THRESH:
            all_ruled = False
        print(f"  [L={L}] max_Up(k< L/4) = {worst['Up']:.6f}  "
              f"{'<= 1/6 OK' if worst['Up'] <= THRESH else 'VIOLATION'}   "
              f"{time.time()-t0:.1f}s", flush=True)
    result["v2_crossing_ruleout"] = cross
    result["v2_all_ruled"] = all_ruled

    # ---------- cross-check the MAIN headline JSON ----------
    with open(r"E:\Projects\MillenniumResearch\prize_work\ppl-014\d0all_harden.json",
              "r", encoding="utf-8") as fh:
        main = json.load(fh)
    main30 = main["per_L"]["30"]["max_Up"]
    this30 = cross["30"]["max_Up"]
    xc = {
        "main_crossing_ruled_out": main["crossing_ruled_out_all_L"],
        "main_pin_L40": main["pin_d0all_L40"],
        "main_d0all_L16": main["gate_exact_d0all_L16"]["d_frac"],
        "this_d0all_pin_crosses": v1_pin,
        "this_all_ruled": all_ruled,
        "L30_main_maxUp": main30,
        "L30_this_maxUp": this30,
        "L30_match": abs(main30 - this30) < 1e-12,
    }
    result["crosscheck_main"] = xc

    result["all_ok"] = bool(result["omega_agree"] and v1_pin and v1_no_lower
                            and all_ruled and v3_upper and v3_k1exact
                            and v3_code and xc["L30_match"]
                            and main["crossing_ruled_out_all_L"]
                            and main["pin_d0all_L40"] == "1/4")
    result["elapsed_s"] = round(time.time() - t0, 2)

    out = (r"E:\Projects\MillenniumResearch\prize_work\ppl-014"
           r"\d0all_harden_verify.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
    print(f"\n[done] all_ok={result['all_ok']}  omega_agree="
          f"{result['omega_agree']}  v2_all_ruled={all_ruled}  "
          f"L30_match={xc['L30_match']}  elapsed={result['elapsed_s']}s",
          flush=True)
    assert result["all_ok"], "INDEPENDENT VERIFIER FAILED"


if __name__ == "__main__":
    main()
