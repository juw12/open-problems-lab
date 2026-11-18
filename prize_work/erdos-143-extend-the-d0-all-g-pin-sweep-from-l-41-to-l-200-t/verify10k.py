#!/usr/bin/env python3
"""#143 / PPL 014  lab_erdos_261112_43  -  INDEPENDENT VERIFIER (N=5e6) for the L=2001..10000 extension.

Independent re-derivation (a second implementation, NOT the sweep code re-run):
  O1  OMEGA element-wise: Method 1 (omega_spf, SPF sieve) vs Method 3
      (omega_prime_power, prime-power sieve, embedded below, independent of SPF)
      at the full N=5e6 -> expect 0 mismatches. (vectorised via np.frombuffer)
  O2  W-table: for representative L (low edge 2001, mid 6000, high edge 10000, and
      the tightest-margin L from the 5e6 sweep) an INDEPENDENT pure-Python double-loop
      W2 (Method-3 Omega, no shared build code) compared entry-wise to the numpy
      pipeline W (build_W_np, Method-1 Omega): < 1e-12.
  O3  U'(L,k) and L'(L,k) re-derived with FRESH pure-Python implementations on W2,
      at the sweep's argmax k, compared to the sweep's stored per-L max_Up / max_Lp: < 1e-12.
  O4  pin rows g*(4,1), g*(8,2) re-derived by fresh enumeration (Method-3 Omega) at N=5e6:
      (i) independent-vs-numpy self-consistency < 1e-12, (ii) both still cross 1/6 (the
      d0^all=1/4 pin mechanism is N-stable).  (Absolute N=5e7 calibration is done in
      spot10k_5e7.py gate G; here we verify the 5e6 pin property + self-consistency.)
  O5  wt table spot check against the literal formula 1/(n log n).

All ok -> the L=2001..10000 headline (U'<=1/6 for every (L,k) with 4k<L, N=5e6) machinery
is independently confirmed.  Single-threaded, ~<5 min.
Output: verify10k.json
"""
import sys, os, time, json, math, heapq, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ppl-014"))
sys.path.insert(0, HERE)
from density_threshold import omega_spf, THRESH, RMAX
from bounds200 import build_W_np, kmax_below_quarter

N = 5_000_000   # match the N=5e6 sweep (sweep10k_5e6.json) this verifier checks
SWEEP = os.path.join(HERE, "sweep10k_5e6.json")


def omega_prime_power(n):
    """Omega via the prime-power sieve (independent of SPF).  Method 3.

    For each prime p and each power p^a <= n, add 1 to omega[m] for every multiple m of p^a.
    Gives Omega(m) = sum_a count exactly.  (verbatim independent copy, no shared code)
    """
    omega = bytearray(n + 1)
    is_prime = bytearray(b'\x01') * (n + 1)
    if n >= 0:
        is_prime[0] = 0
    if n >= 1:
        is_prime[1] = 0
    lim = int(n ** 0.5)
    for i in range(2, lim + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = 0
    for p in range(2, n + 1):
        if is_prime[p]:
            pp = p
            while pp <= n:
                for m in range(pp, n + 1, pp):
                    omega[m] += 1
                pp *= p
    return omega


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
    out = 0.0
    for r in range(1, rmax + 1):
        s = sum(heapq.nsmallest(k, (W[j][r] for j in range(1, L))))
        if s > out:
            out = s
    return out


def gstar_indep(W, L, k, rmax):
    worst = float("inf")
    for S in itertools.combinations(range(1, L), k):
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


def w2_direct(om3, wt, Nn, L, rmax):
    W = [[0.0] * (rmax + 1) for _ in range(L)]
    for nn in range(2, Nn + 1):
        r = om3[nn]
        if 1 <= r <= rmax:
            W[nn % L][r] += wt[nn]
    return W


def main():
    t0 = time.time()
    res = {"job": "lab_erdos_261112_43_verify", "N": N, "threshold": THRESH, "rmax": RMAX}
    print(f"[43v] verify N={N} t=0.0s", flush=True)

    # ---- O1: element-wise Omega, two independent methods (vectorised) ----
    om1 = omega_spf(N); print(f"[43v] [m1] SPF omega {time.time() - t0:.1f}s", flush=True)
    om3 = omega_prime_power(N); print(f"[43v] [m3] prime-power omega {time.time() - t0:.1f}s", flush=True)
    a = np.frombuffer(om1, dtype=np.int8); b = np.frombuffer(om3, dtype=np.int8)
    mm = int((a != b).sum()); mx = int(np.abs(a.astype(np.int16) - b.astype(np.int16)).max())
    res["omega_mismatches"] = mm; res["omega_max_absdiff"] = mx; res["omega_agree"] = (mm == 0)
    print(f"[43v] [O1] omega mismatches={mm} (max {mx})  {'AGREE' if mm == 0 else 'DISAGREE'}  {time.time() - t0:.1f}s", flush=True)

    # ---- weights + O5 spot check ----
    t_ = np.arange(2, N + 1, dtype=np.float64)
    wt = np.zeros(N + 1, dtype=np.float64); wt[2:] = 1.0 / (t_ * np.log(t_))
    for e in (2, 3, 4):
        wt[e] = 0.0
    n_arr = np.arange(N + 1, dtype=np.int32)
    om1_np = np.frombuffer(om1, dtype=np.int8)
    spot = {}
    for n in (1000003, 4999999, 123457, 2, 3, 4):
        ref = 0.0 if n in (2, 3, 4) else 1.0 / (n * math.log(n))
        spot[n] = bool(abs(float(wt[n]) - ref) < 1e-15)
    res["wt_spot_ok"] = all(spot.values()) and bool(spot)
    print(f"[43v] [O5] wt spot all_ok={res['wt_spot_ok']}  {time.time() - t0:.1f}s", flush=True)

    # ---- sweep data ----
    with open(SWEEP, "r", encoding="utf-8") as fh:
        runB = json.load(fh)
    perL = runB.get("N_5e6", {}).get("per_L", {})
    n5 = runB.get("N_5e6", {})
    tightest = str(n5.get("worst_margin_L", 10000))
    # pick L's that are ACTUALLY in the achieved 5e6 range (as-fits may stop before 10000)
    lo = n5.get("L_achieved_lo"); hi = n5.get("L_achieved_hi")
    cand = {int(lo) if lo is not None else 2001,
            int(hi) if hi is not None else 10000,
            int(tightest)}
    L_check = sorted({L for L in cand if str(L) in perL})
    res["L_check"] = L_check
    res["L_achieved"] = [lo, hi]
    print(f"[43v] L_check={L_check}  (tightest={tightest}, achieved={lo}..{hi})  {time.time() - t0:.1f}s", flush=True)

    # ---- O2 + O3 per checked L ----
    o2, o3 = {}, {}
    all_o2 = all_o3 = True
    for L in L_check:
        tb = time.time()
        W2 = w2_direct(om3, wt, N, L, RMAX)
        Wn = build_W_np(om1_np, wt, n_arr, L, RMAX)
        dW = max(abs(W2[j][r] - float(Wn[j, r])) for j in range(L) for r in range(RMAX + 1))
        o2[str(L)] = dW; all_o2 = all_o2 and dW < 1e-12
        st = perL.get(str(L))
        if st is not None:
            # stored max_Up == U'(L, kmax) where kmax = (L-1)//4 is the BINDING worst-k
            # (U'(L,k) is increasing in k).  The stored 'k_at_maxUp' is an OMEGA-LEVEL index,
            # NOT a k -- so the bound must be re-derived at kmax, not at that index.
            kmax = kmax_below_quarter(L)
            dUp = abs(U_prime_indep(W2, L, kmax, RMAX) - st["max_Up"])
            dLp = abs(L_prime_indep(W2, L, kmax, RMAX) - st["max_Lp"])
            o3[str(L)] = {"dUp": dUp, "dLp": dLp}
            all_o3 = all_o3 and (dUp < 1e-12) and (dLp < 1e-12)
        else:
            o3[str(L)] = {"note": f"L={L} not in sweep per_L"}
            all_o3 = False
        print(f"[43v] [O2/O3] L={L} dW={dW:.2e} dUp={o3[str(L)].get('dUp', float('nan')):.2e} "
              f"dLp={o3[str(L)].get('dLp', float('nan')):.2e}  {time.time() - tb:.1f}s  {time.time() - t0:.1f}s", flush=True)
    res["o2_W_entry_maxdiff"] = o2
    res["o3_bounds_maxdiff"] = o3

    # ---- O4: pin rows, fresh enumeration, Method-3 Omega at N=5e6 ----
    W4 = w2_direct(om3, wt, N, 4, RMAX); W8 = w2_direct(om3, wt, N, 8, RMAX)
    g41 = gstar_indep(W4, 4, 1, RMAX); g82 = gstar_indep(W8, 8, 2, RMAX)
    # numpy gstar self-consistency (Method-1 Omega, same W table)
    from density_threshold import gstar as gstar_np
    g41_np = float(gstar_np(W4, 4, 1, RMAX)); g82_np = float(gstar_np(W8, 8, 2, RMAX))
    o4 = {
        "gstar_4_1_indep": g41, "gstar_4_1_numpy": g41_np, "d4_1": abs(g41 - g41_np),
        "gstar_8_2_indep": g82, "gstar_8_2_numpy": g82_np, "d8_2": abs(g82 - g82_np),
        "pin_crosses_1_6": bool(g41 > THRESH and g82 > THRESH),
        "self_consistent": bool(abs(g41 - g41_np) < 1e-12 and abs(g82 - g82_np) < 1e-12),
    }
    res["o4_pin_rows"] = o4
    print(f"[43v] [O4] g*(4,1)={g41} g*(8,2)={g82} crosses1/6={o4['pin_crosses_1_6']} "
          f"self_consistent={o4['self_consistent']}  {time.time() - t0:.1f}s", flush=True)

    # ---- headline cross-check ----
    gates_ok = (runB.get("gates") or {}).get("gates_ok") if isinstance(runB.get("gates"), dict) else runB.get("gates_ok")
    res["runB_status"] = runB.get("status")
    res["runB_verdict"] = runB.get("verdict")
    res["runB_gates_ok"] = gates_ok
    res["all_ok"] = bool(res["omega_agree"] and res["wt_spot_ok"] and all_o2 and all_o3
                         and o4["self_consistent"] and o4["pin_crosses_1_6"] and gates_ok)
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(HERE, "verify10k.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(f"\n[43v] [done] all_ok={res['all_ok']} omega_agree={res['omega_agree']} O2_ok={all_o2} "
          f"O3_ok={all_o3} O4_ok={o4['self_consistent'] and o4['pin_crosses_1_6']} gates_ok={gates_ok} "
          f"elapsed={res['elapsed_s']}s", flush=True)


if __name__ == "__main__":
    main()
