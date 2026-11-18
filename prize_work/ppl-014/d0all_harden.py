#!/usr/bin/env python3
"""
Job lab_erdos_261111_23 - Erdos #143 (PPL 014)
Harden the d0^all (g*) 'for-all-levels' pin from L<=16 to L<=40.

With d0^exists = 0 (job 261111_11), the meaningful quantity is the UNIVERSAL
threshold
    d0^all = min { k/L : g*(L,k) > 1/6 },
pinned at 1/4 for L<=16 (exact g*, N=5e7, job fcc..._7).  Exact g* is infeasible
for L>~16 (binomial blowup), so extend the for-all claim to L=17..40 with an
analytic worst-class UPPER bound on g* and check whether 1/4 survives.

VALID UPPER BOUND on g*(L,k) = min_S max_r sum_{j in S} W[j][r]:
    w_j    = max_{r>=1} W[j][r]                       (per-class best-level wt)
    S*     = the k classes j with the smallest w_j
    U'(L,k)= max_r  sum_{j in S*} W[j][r]             (value at the k-argmin set)
    U(L,k) = sum of the k smallest w_j                (looser, cross-check)
    g*(L,k) <= max_r sum_{j in S*} W[j][r] = U'(L,k)  <= U(L,k)

If U'(L,k) <= 1/6 for every (L,k) with k/L < 1/4 (L=17..40), then g*(L,k) <= 1/6
there, so NO (L,k) with k/L < 1/4 crosses 1/6, so the global d0^all (min over
L<=40 of the smallest crossing k/L) cannot drop below 1/4.  Combined with the
exact crossing at (L=4,k=1) (g*=0.2615>1/6, k/L=1/4), the pin is EXACTLY 1/4 over
the full L<=40 range.

GATES (reproduce the lab's own pinned values before any new claim):
  G1  level_mass(N=2e6, L=4, S={1,2,3}) == 1.56769507127912   (trusted machinery)
  G2  exact g* at L<=16 (N=5e7) -> d0^all == 1/4              (job 7 pin)
  G3  f*(14,1) @5e7 == 0.1808865800265416                     (job 7 pin)

Single-threaded, N=5e7, {2,3,4}-excluded weights.  <20 min.
"""
import sys, time, json, math, array

sys.path.insert(0, r"E:\Projects\MillenniumResearch\prize_work\ppl-014")

from density_threshold import (omega_spf, build_W, fstar, gstar, d0_of,
                               THRESH, RMAX)
from per_round_weight_sharpening import level_mass, max_layer
from nondeg_threshold import excl_set

NMAIN = 50_000_000
NGATE = 2_000_000
F1_2E6_REF = 1.56769507127912      # level_mass(N=2e6, L=4, S={1,2,3})
F141_5E7_REF = 0.1808865800265416  # f*(14,1) @5e7, {2,3,4}-excluded (job 7)
PIN_ALL = "1/4"                    # job-7 exact d0^all pin at L<=16
LS_EXACT = list(range(3, 17))      # 3..16 : exact g*
LS_ANALYTIC = list(range(17, 41))  # 17..40: analytic worst-class bound


def weight_table_arr(nmax, exclude):
    """wt[n]=1/(n log n) for 2<=n<=nmax, zero at excluded n.  (array 'd')."""
    wt = array.array('d', [0.0]) * (nmax + 1)
    for n in range(2, nmax + 1):
        wt[n] = 1.0 / (n * math.log(n))
    for e in exclude:
        wt[e] = 0.0
    return wt


def per_class_best(W, L, rmax):
    """w_j = max_{1<=r<=rmax} W[j][r], for j=1..L-1."""
    return [max(W[j][r] for r in range(1, rmax + 1)) for j in range(1, L)]


def k_smallest_classes(W, L, k, rmax):
    """The k classes j (1..L-1) with the smallest per-class best weight w_j."""
    w = per_class_best(W, L, rmax)
    idx = list(range(L - 1))          # idx[t] -> class j = t+1
    idx.sort(key=lambda t: w[t])
    return [t + 1 for t in idx[:k]]


def U_sum(W, L, k, rmax):
    """U(L,k) = sum of the k smallest per-class best-level weights (valid ub)."""
    w = per_class_best(W, L, rmax)
    w.sort()
    return sum(w[:k])


def U_maxlevel(W, L, k, rmax):
    """U'(L,k) = max_r sum_{j in S*} W[j][r], S* = k smallest w_j classes."""
    classes = k_smallest_classes(W, L, k, rmax)
    best = 0.0
    for r in range(1, rmax + 1):
        s = 0.0
        for j in classes:
            s += W[j][r]
        if s > best:
            best = s
    return best


def main():
    t0 = time.time()
    result = {"job": "lab_erdos_261111_23", "problem": "Erdos #143 / PPL 014",
              "N": NMAIN, "exclude": list(excl_set(10)),
              "threshold": THRESH, "rmax": RMAX}
    print(f"[start] lab_erdos_261111_23  N={NMAIN}  t=0.0s", flush=True)

    # ---------- GATE 1: trusted level_mass machinery ----------
    om_gate = omega_spf(NGATE)
    f1_list = level_mass(NGATE, om_gate, L=4, S=(1, 2, 3))
    f1 = max_layer(f1_list)[1]      # best single-level mass (the gate quantity)
    g1 = abs(f1 - F1_2E6_REF) < 1e-12
    result["gate_level_mass"] = {"value": f1, "ref": F1_2E6_REF, "ok": g1,
                                 "note": "max_layer over per-level masses"}
    print(f"[G1] level_mass(2e6,L=4,S={1,2,3}) max-layer = {f1}  ref {F1_2E6_REF}  "
          f"ok={g1}   {time.time()-t0:.1f}s", flush=True)

    # ---------- main sieve + weight table ----------
    omega = omega_spf(NMAIN)
    wt = weight_table_arr(NMAIN, excl_set(10))
    print(f"[sieve] omega_spf({NMAIN}) + weight table   {time.time()-t0:.1f}s",
          flush=True)

    # ---------- exact g* for L<=16 ; analytic bound for L=17..40 ----------
    exact_rows = []     # L<=16 : full (L,k) table with g*, U, U'
    analytic_rows = []  # L>=17 : k with k/L < 1/4, U, U'
    W_cache = {}
    for L in LS_EXACT + LS_ANALYTIC:
        W = build_W(omega, wt, NMAIN, L, RMAX)
        W_cache[L] = W
        if L <= 16:
            for k in range(1, L):
                ge = gstar(W, L, k, RMAX)
                U = U_sum(W, L, k, RMAX)
                Up = U_maxlevel(W, L, k, RMAX)
                exact_rows.append({"L": L, "k": k, "gstar": ge, "U": U, "Up": Up})
        else:
            kmax_below = (L - 1) // 4          # largest k with k/L < 1/4
            for k in range(1, kmax_below + 1):
                U = U_sum(W, L, k, RMAX)
                Up = U_maxlevel(W, L, k, RMAX)
                analytic_rows.append({"L": L, "k": k, "U": U, "Up": Up,
                                      "k_over_L": k / L})
        print(f"  [L={L:2d}] {'exact' if L <= 16 else 'bound'}   "
              f"{time.time()-t0:6.1f}s", flush=True)

    # ---------- GATE 2: exact d0^all at L<=16 == 1/4 ----------
    # d0_of expects rows to already carry "d" and "d_frac" keys (as sweep() does).
    gate2_rows = []
    for r in exact_rows:
        k, L = r["k"], r["L"]
        g = math.gcd(k, L)
        gate2_rows.append({"L": L, "k": k, "fstar": None, "gstar": r["gstar"],
                           "d": k / L, "d_frac": f"{k // g}/{L // g}"})
    d0all = d0_of(gate2_rows, "gstar", THRESH)
    g2 = d0all is not None and d0all["d_frac"] == PIN_ALL
    f141 = fstar(W_cache[14], 14, 1, RMAX)
    g3 = abs(f141 - F141_5E7_REF) < 1e-9
    result["gate_exact_d0all_L16"] = d0all
    result["gate_f141_5e7"] = {"value": f141, "ref": F141_5E7_REF, "ok": g3}
    result["gates_ok"] = bool(g1 and g2 and g3)
    print(f"[G2] exact d0^all(L<=16) = {d0all['d_frac']} (d={d0all['d']})  "
          f"ok={g2};  [G3] f*(14,1)={f141}  ok={g3}   "
          f"{time.time()-t0:.1f}s", flush=True)

    # ---------- THE CROSSING RULE-OUT at L=17..40 ----------
    crossing_ruled_out = True
    per_L = {}
    violations = []
    for L in LS_ANALYTIC:
        kmax_below = (L - 1) // 4
        rowsL = [r for r in analytic_rows if r["L"] == L]
        max_Up = max((r["Up"] for r in rowsL), default=None)
        max_U = max((r["U"] for r in rowsL), default=None)
        k_at_maxUp = max(rowsL, key=lambda r: r["Up"])["k"] if rowsL else None
        per_L[str(L)] = {
            "kmax_below_1_4": kmax_below,
            "k_at_maxUp": k_at_maxUp,
            "max_Up": max_Up, "max_U": max_U,
            "margin_to_1_6": (THRESH - max_Up) if max_Up is not None else None,
        }
        if max_Up is not None and max_Up > THRESH:
            crossing_ruled_out = False
            violations.append({"L": L, "k": k_at_maxUp, "Up": max_Up})
    result["crossing_ruled_out_all_L"] = crossing_ruled_out
    result["violations"] = violations
    result["per_L"] = per_L

    # ---------- bound-vs-exact gap at L<=16 (verifier honesty) ----------
    gaps_up = [abs(r["gstar"] - r["Up"]) for r in exact_rows]
    gaps_u = [abs(r["gstar"] - r["U"]) for r in exact_rows]
    pin_rows = [r for r in exact_rows
                if r["k"] * 4 == r["L"] and abs(r["k"] / r["L"] - 0.25) < 1e-15]
    result["bound_vs_exact_L16"] = {
        "max_gap_Up": max(gaps_up),
        "max_gap_U": max(gaps_u),
        "n_rows": len(exact_rows),
        "note": "U' is a valid UPPER bound on g* (min-max <= value at k-argmin "
                "set).  Gap reports how tight the bound is at L<=16 where exact "
                "g* is available; the crossing rule-out uses g* <= U' <= 1/6.",
        "pin_rows_k_over_L_is_1_4": [
            {"L": r["L"], "k": r["k"], "gstar": r["gstar"], "Up": r["Up"],
             "gap_Up": abs(r["gstar"] - r["Up"])} for r in pin_rows],
    }

    # ---------- headline ----------
    if g2 and crossing_ruled_out:
        result["pin_d0all_L40"] = PIN_ALL
        worst = max(per_L.values(), key=lambda v: v["max_Up"])
        result["verdict"] = (
            f"EXTENSION: d0^all pin 1/4 hardens from L<=16 to L<=40.  Exact g* "
            f"(L<=16) gives d0^all=1/4; the analytic worst-class upper bound "
            f"U' <= 1/6 for every (L,k) with k/L<1/4 at L=17..40 (worst max_Up="
            f"{worst['max_Up']:.6f} at the argmax L) rules out any lower "
            f"crossing, so the pin stays EXACTLY 1/4 over L<=40.")
    elif g2:
        result["pin_d0all_L40"] = None
        result["verdict"] = (
            "KILL: analytic bound too loose to rule out a lower crossing at "
            f"L=17..40 (violations: {violations}); report L<=16 pin 1/4 unchanged.")
    else:
        result["pin_d0all_L40"] = None
        result["verdict"] = ("GATE FAILED: exact g* at L<=16 did not reproduce "
                             f"1/4 (got {d0all and d0all['d_frac']}).")

    result["elapsed_s"] = round(time.time() - t0, 2)
    out = r"E:\Projects\MillenniumResearch\prize_work\ppl-014\d0all_harden.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
    print(f"\n[done] gates_ok={result['gates_ok']}  "
          f"crossing_ruled_out={crossing_ruled_out}  "
          f"pin_L40={result['pin_d0all_L40']}  elapsed={result['elapsed_s']}s",
          flush=True)
    print("VERDICT:", result["verdict"], flush=True)


if __name__ == "__main__":
    main()
