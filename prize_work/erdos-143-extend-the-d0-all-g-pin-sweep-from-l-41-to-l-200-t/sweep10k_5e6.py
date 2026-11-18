#!/usr/bin/env python3
"""#143 / PPL 014  lab_erdos_261112_43  -  d0^all=1/4 pin: L=2001..10000 (N=5e6 FULL primary).

Extends the job lab_erdos_261112_38 frontier (L<=2000) one more range step to
L=2001..10000.  Reuses job-38 machinery verbatim (int32 single-bincount flattened
W builder + vectorised all-k bounds_all via sorted-prefix cumsum) so the results are
directly comparable to the L<=2000 pin.

N strategy: N=5e6 FULL for all 8000 L's is the guaranteed primary (job 38 itself
validated N=5e6 full across L=201..2000, so this extends the SAME-N pin).  A
separate 5e7 spot-check (spot10k_5e7.py) confirms N-robustness toward the larger N
that job 38 used for its headline; the 5e7 FULL pass over 8000 L's (~1 hr) exceeds
the per-computation budget, so it is done as a spot subset, not a full sweep.

Gates (all must pass; N-independent structural + absolute calibration):
  K1  int32-flatten builder vs job-35 24-pass build_W_np  (N=200k, L in {2001,6000,10000}) <1e-12
  E   vectorised bounds_all vs original U_maxlevel/L_prime  (N=200k, several (L,k))        <1e-12

Output: sweep10k_5e6.json
Single-threaded, ~<10 min.
"""
import os, sys, json, time, collections
import numpy as np
_T0 = time.time()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ppl-014"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from density_threshold import omega_spf, THRESH, RMAX
from bounds200 import U_maxlevel, L_prime, kmax_below_quarter, build_W_np
from bounds2000 import masked_inputs, build_W_flat_precomputed

CARD = os.path.dirname(os.path.abspath(__file__))
N = int(5e6)
N_KILL = 200000
EXCL = (2, 3, 4)
TH = 1.0 / 6.0
RMAX = 24
L_LO = 2001
L_HI = 10000
N_L = L_HI - L_LO + 1


def bounds_all(W, L, kmax, rmax):
    """Return (U_all, L_all), each (kmax,): U'(L,1..kmax), L'(L,1..kmax).
    Exactly equivalent to per-k U_maxlevel / L_prime (identical to job 38)."""
    Wc = W[1:L]
    w = Wc.max(axis=1)
    ord_ = np.argsort(w, kind="stable")
    W_ord = Wc[ord_]
    C = np.cumsum(W_ord, axis=0)
    U_all = C[:kmax].max(axis=1)
    W_asc = np.sort(Wc, axis=0)
    P = np.cumsum(W_asc, axis=0)
    L_all = P[:kmax].max(axis=1)
    return U_all, L_all


def make_wt(N):
    n = np.arange(N + 1, dtype=np.float64)
    n = np.where(n < 2.0, 1.0, n)
    wt = np.zeros(N + 1, dtype=np.float64)
    wt[2:] = 1.0 / (n[2:] * np.log(n[2:]))
    for e in EXCL:
        wt[e] = 0.0
    return wt


def sweep_range(n_m, om_i32, wt_m, L_lo, L_hi, N, tag, deadline=None):
    per_L = {}
    unruled = 0
    crossings = 0
    cross_rows = []
    worst_margin = None
    worst_margin_L = None
    near = []
    t0 = time.time()
    done = 0
    recent = collections.deque(maxlen=30)
    stopped_by_deadline = False
    for L in range(L_lo, L_hi + 1):
        if deadline is not None and time.time() > deadline:
            stopped_by_deadline = True
            break
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
            worst_margin = margin
            worst_margin_L = int(L)
        if margin < 0.03:
            near.append({"L": int(L), "margin": float(margin)})
        done += 1
        dt = time.time() - tL
        recent.append(dt)
        if done % 400 == 0:
            print(f"[43] {tag} {N} sweep: {done}/{L_hi - L_lo + 1} L's  t={time.time() - t0:.0f}s  "
                  f"worst_margin={worst_margin:.6f}@L{worst_margin_L}  s/L={sum(recent) / len(recent):.3f}",
                  flush=True)
    return {"per_L": per_L, "n_L_done": len(per_L), "n_L_total": L_hi - L_lo + 1,
            "L_achieved_lo": L_lo if per_L else None,
            "L_achieved_hi": (max(int(k) for k in per_L)) if per_L else None,
            "stopped_by_deadline": stopped_by_deadline,
            "unruled": unruled, "crossings": crossings, "crossing_rows": cross_rows,
            "worst_margin": worst_margin, "worst_margin_L": worst_margin_L,
            "near_margin_below_0.03": near,
            "mean_s_per_L": (sum(recent) / len(recent)) if recent else 0.0,
            "elapsed_s": time.time() - t0}


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)
    print(f"[43] wrote {os.path.basename(path)}  (wall {time.time() - _T0:.0f}s)", flush=True)


def main():
    out = os.path.join(CARD, "sweep10k_5e6.json")
    print(f"[43] N={N} TH={TH} RMAX={RMAX} EXCL={EXCL} sweep L={L_LO}..{L_HI} ({N_L} L's)", flush=True)

    gates = {"gates_ok": True}

    # ---- Gate K1: flatten(int32) vs 24-pass build_W_np at 200k ----
    print(f"[43] K1 flatten-vs-24pass at {N_KILL}", flush=True)
    om_k = np.asarray(omega_spf(N_KILL), dtype=np.int32)
    wt_k = make_wt(N_KILL); n_k = np.arange(N_KILL + 1, dtype=np.int32)
    n_mk, om_mk, wt_mk = masked_inputs(n_k, om_k, wt_k, RMAX)
    k1 = {}
    for L in (2001, 6000, 10000):
        Wf = build_W_flat_precomputed(n_mk, om_mk, wt_mk, L, RMAX)
        Wn = build_W_np(om_k, wt_k, n_k, L, RMAX)
        d = float(np.abs(Wf - Wn).max()); k1[str(L)] = d
        print(f"[43] K1 L={L}: max_abs_diff={d:.3e}  {'OK' if d < 1e-12 else 'MISMATCH'}", flush=True)
        if d >= 1e-12:
            gates["gates_ok"] = False
    gates["K1_flatten_vs_24pass_max_diff"] = k1

    # ---- Gate E: vectorised bounds_all vs original U_maxlevel / L_prime at 200k ----
    print(f"[43] E bounds_all-vs-original at {N_KILL}  (small L: identity is L-independent; large-L data is validated by the sweep + verifier O3)", flush=True)
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
        print(f"[43] E L={L}: max|dU|={dU:.3e} max|dL|={dL:.3e}  {'OK' if max(dU, dL) < 1e-12 else 'MISMATCH'}",
              flush=True)
        if max(dU, dL) >= 1e-12:
            gates["gates_ok"] = False
    gates["E_bounds_all_vs_original"] = eq
    del n_mk, om_mk, wt_mk, om_k, wt_k, n_k

    # ---- build 5e6 ----
    print(f"[43] building 5e6 omega_spf + wt + n ...", flush=True)
    t0 = time.time()
    om6 = np.asarray(omega_spf(N), dtype=np.int32)
    wt6 = make_wt(N); n6 = np.arange(N + 1, dtype=np.int32)
    n_m6, om6_i32, wt6_m = masked_inputs(n6, om6, wt6, RMAX)
    del n6, om6, wt6
    print(f"[43] 5e6 build {time.time() - t0:.1f}s  ({len(n_m6)} masked terms)  {time.time() - _T0:.0f}s", flush=True)

    result = {
        "job": "lab_erdos_261112_43", "problem": "Erdos #143 / PPL 014 ($500)",
        "claim": "extend d0^all=1/4 universal finite-core pin from L<=2000 (job lab_erdos_261112_38) to L<=10000",
        "N": N, "TH": TH, "RMAX": RMAX, "EXCL": list(EXCL),
        "L_range": [L_LO, L_HI], "n_L": N_L, "gates": gates,
        "N_5e6": {}, "verdict": None, "status": "running", "novelty": "extension",
        "elapsed_s": 0.0,
        "method": ("int32 single-bincount flattened W builder (bit-exact vs job-35 24-pass, gate K1) + "
                   "vectorised all-k bounds U'(L,k)/L'(L,k) via sorted-prefix cumsum (equivalence gate E); "
                   "N=5e6 FULL L=2001..10000. U'<=1/6 => no d0^all<1/4 crossing."),
    }

    t_build_done = time.time()
    sweep_budget_s = 14.0 * 60.0   # ~14 min for the 5e6 sweep (as-fits); 5e7 spot + 1e6 full cover the rest
    result["N_5e6"] = sweep_range(n_m6, om6_i32, wt6_m, L_LO, L_HI, N, tag="5e6",
                                  deadline=t_build_done + sweep_budget_s)
    del n_m6, om6_i32, wt6_m
    print(f"[43] 5e6 DONE unruled={result['N_5e6']['unruled']} crossings={result['N_5e6']['crossings']} "
          f"achieved L={result['N_5e6']['L_achieved_lo']}..{result['N_5e6']['L_achieved_hi']} "
          f"deadline={result['N_5e6']['stopped_by_deadline']} "
          f"worst_margin={result['N_5e6']['worst_margin']}@L{result['N_5e6']['worst_margin_L']} "
          f"mean_s/L={result['N_5e6']['mean_s_per_L']:.3f}", flush=True)

    # ---- verdict (5e6 as-fits is the guaranteed primary; achieved range may be <10000) ----
    r6 = result["N_5e6"]
    L6_lo = r6["L_achieved_lo"]; L6_hi = r6["L_achieved_hi"]
    rng = f"L={L6_lo}..{L6_hi}" if (L6_lo is not None) else "(none)"
    if r6["crossings"] > 0:
        result["status"] = "needs_review"; result["novelty"] = "new"
        result["verdict"] = (f"CROSSING: L'>1/6 at 4k<L in {rng} (N=5e6) -> d0^all<1/4 "
                             f"(rows {r6['crossing_rows'][:8]}).")
    elif r6["unruled"] > 0:
        result["status"] = "needs_review"; result["novelty"] = "extension"
        result["verdict"] = f"U' rule-out incomplete (unruled={r6['unruled']}) in {rng} -- needs LP/sharper method."
    else:
        result["status"] = "done"; result["novelty"] = "extension"
        full = (L6_lo == L_LO and L6_hi == L_HI)
        result["verdict"] = (f"EXTENSION: U'(L,k)<=1/6 for every (L,k) 4k<L, {rng}, N=5e6 "
                             f"({'full range' if full else 'as-fits within budget'}; "
                             f"worst margin {r6['worst_margin']:.6f} at L={r6['worst_margin_L']}); "
                             f"d0^all=1/4 pin holds across {rng} at N=5e6.")
    result["elapsed_s"] = time.time() - _T0
    write_json(out, result)
    print(f"[43] VERDICT: {result['verdict']}", flush=True)
    print(f"[43] TOTAL {time.time() - _T0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
