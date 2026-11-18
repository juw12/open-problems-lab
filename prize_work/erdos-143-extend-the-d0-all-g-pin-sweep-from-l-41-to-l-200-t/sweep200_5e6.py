#!/usr/bin/env python3
"""
Job lab_erdos_261112_35 - RUN A: triage sweep L=41..200 at N=5e6 (pure Python).

Reuses job-261111_23 machinery (prize_work/ppl-014/): omega_spf, build_W,
gstar, d0_of, level_mass gate.  Pure-Python build_W (~0.9 s/L at 5e6).

GATES (before any new claim):
  G1  level_mass(N=2e6, L=4, S=(1,2,3)) max-layer == 1.56769507127912
  G2  exact g* at L=3..16 (N=5e6) -> d0^all == 1/4
      (cross rows (4,1),(8,2); bracket_lo 3/13)  [pin robustness 1e6..5e7
       was established by earlier jobs; this re-derives it at the sweep N]

SWEEP: L=41..200, k=1..(L-1)//4 (4k<L):  U'(L,k), U(L,k), L'(L,k)
  - U' > 1/6  -> 'un-ruled' (upper bound too loose at that L)
  - L' > 1/6  -> PROVEN crossing below 1/4 (a win, needs_review)
  - near-zero: margin 1/6 - max_U' < 0.03  (candidates for 5e7 re-check;
    run B re-checks ALL L at 5e7 anyway, so this is informational)

Kill-test output: per-L build_W timing -> projected 5e7 sweep cost.
Single-threaded.  Output: sweep200_5e6.json
"""
import sys, time, json, math, array

HERE = r"E:\Projects\MillenniumResearch\prize_work\erdos-143-extend-the-d0-all-g-pin-sweep-from-l-41-to-l-200-t"
sys.path.insert(0, r"E:\Projects\MillenniumResearch\prize_work\ppl-014")
sys.path.insert(0, HERE)

from density_threshold import (omega_spf, build_W, gstar, d0_of, THRESH, RMAX)
from per_round_weight_sharpening import level_mass, max_layer
from nondeg_threshold import excl_set
from bounds200 import (U_maxlevel, U_sum, L_prime, kmax_below_quarter)

N = 5_000_000
NGATE = 2_000_000
F1_2E6_REF = 1.56769507127912
LS_SWEEP = list(range(41, 201))
LS_EXACT = list(range(3, 17))


def weight_table_arr(nmax, exclude):
    wt = array.array('d', [0.0]) * (nmax + 1)
    for n in range(2, nmax + 1):
        wt[n] = 1.0 / (n * math.log(n))
    for e in exclude:
        wt[e] = 0.0
    return wt


def main():
    t0 = time.time()
    res = {"job": "lab_erdos_261112_35_runA", "N": N, "exclude": list(excl_set(10)),
           "threshold": THRESH, "rmax": RMAX}
    print(f"[start] runA N={N} t=0.0s", flush=True)

    # ---------- GATE 1 ----------
    om_gate = omega_spf(NGATE)
    f1 = max_layer(level_mass(NGATE, om_gate, L=4, S=(1, 2, 3)))[1]
    g1 = abs(f1 - F1_2E6_REF) < 1e-12
    res["gate_level_mass"] = {"value": f1, "ref": F1_2E6_REF, "ok": g1}
    print(f"[G1] level_mass(2e6) max-layer={f1} ok={g1}  {time.time()-t0:.1f}s",
          flush=True)

    # ---------- sieve + weights ----------
    omega = omega_spf(N)
    wt = weight_table_arr(N, excl_set(10))
    print(f"[sieve] omega+wt done  {time.time()-t0:.1f}s", flush=True)

    # ---------- GATE 2: exact pin at L<=16, N=5e6 ----------
    rows = []
    for L in LS_EXACT:
        W = build_W(omega, wt, N, L, RMAX)
        for k in range(1, L):
            rows.append({"L": L, "k": k, "gstar": gstar(W, L, k, RMAX),
                         "d": k / L, "d_frac": f"{k}/{L}"})
    d0all = d0_of(rows, "gstar", THRESH)
    g2 = (d0all is not None and d0all["d_frac"] == "1/4"
          and sorted(d0all["cross_rows"]) == [[4, 1], [8, 2]]
          and abs(d0all["bracket_lo"] - 3 / 13) < 1e-12)
    res["gate_exact_d0all_L16_5e6"] = d0all
    res["gate2_ok"] = bool(g2)
    print(f"[G2] exact d0^all(L<=16)@5e6 = {d0all and d0all['d_frac']} "
          f"cross={d0all and d0all['cross_rows']} ok={g2}  {time.time()-t0:.1f}s",
          flush=True)

    # ---------- SWEEP L=41..200 ----------
    per_L = {}
    unruled, proven_crossings = [], []
    build_times = []
    for L in LS_SWEEP:
        tb = time.time()
        W = build_W(omega, wt, N, L, RMAX)
        build_times.append(time.time() - tb)
        kmax = kmax_below_quarter(L)
        rowk = []
        for k in range(1, kmax + 1):
            Up = U_maxlevel(W, L, k, RMAX)
            Lp = L_prime(W, L, k, RMAX)
            rowk.append({"k": k, "Up": Up, "Lp": Lp})
        maxUp = max(rowk, key=lambda r: r["Up"])
        maxLp = max(rowk, key=lambda r: r["Lp"])
        margin = THRESH - maxUp["Up"]
        per_L[str(L)] = {"kmax": kmax, "rows": rowk,
                         "max_Up": maxUp["Up"], "k_at_maxUp": maxUp["k"],
                         "max_Lp": maxLp["Lp"], "k_at_maxLp": maxLp["k"],
                         "margin": margin}
        if maxUp["Up"] > THRESH:
            unruled.append({"L": L, "k": maxUp["k"], "Up": maxUp["Up"]})
        if maxLp["Lp"] > THRESH:
            proven_crossings.append({"L": L, "k": maxLp["k"], "Lp": maxLp["Lp"],
                                     "d": maxLp["k"] / L})
        if L % 20 == 1:
            print(f"  [L={L:3d}] maxUp={maxUp['Up']:.6f} margin={margin:.6f} "
                  f"{time.time()-t0:7.1f}s", flush=True)

    mt = sum(build_times) / len(build_times)
    res.update({
        "gates_ok": bool(g1 and g2),
        "n_L": len(LS_SWEEP),
        "per_L": per_L,
        "unruled_Up_over_1_6": unruled,
        "proven_crossings_Lp_over_1_6": proven_crossings,
        "near_zero_margin_lt_0_03": sorted(
            [L for L in LS_SWEEP if per_L[str(L)]["margin"] < 0.03]),
        "worst_margin_L": min(LS_SWEEP, key=lambda L: per_L[str(L)]["margin"]),
        "build_W_mean_s_5e6": mt,
        "projected_5e7_sweep_min": round(mt * (N * 10) / 60.0, 1),
        "elapsed_s": round(time.time() - t0, 1),
    })
    out = HERE + r"\sweep200_5e6.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n[done] gates_ok={res['gates_ok']} unruled={len(unruled)} "
          f"proven_crossings={len(proven_crossings)} "
          f"worst_margin_L={res['worst_margin_L']} "
          f"(margin={per_L[str(res['worst_margin_L'])]['margin']:.6f}) "
          f"buildW_mean={mt:.2f}s proj5e7={res['projected_5e7_sweep_min']}min "
          f"elapsed={res['elapsed_s']}s", flush=True)


if __name__ == "__main__":
    main()
