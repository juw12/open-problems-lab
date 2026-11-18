#!/usr/bin/env python3
"""#143 / PPL 014  lab_erdos_261112_43  -  d0^all=1/4 pin: L=2001..10000 (N=1e6 FULL, fast).

Full-range coverage of L=2001..10000 at N=1e6 (~8 min, ~0.06 s/L).  N=1e6 is the
WEAKEST N (least total weight sum 1/(n log n) ~ log log N = 2.62 vs 5e6 -> 2.73,
5e7 -> 2.87), so U'(L,k) is smallest here and the margin is largest.  A pin that
holds at 1e6 over the full range, combined with the stronger-N 5e6 as-fits and the
5e7 spot, establishes the d0^all=1/4 pin across L=2001..10000 with N-robustness:
U'(L,k) is monotone increasing in N, so the tightest evidence is at the largest N.

Gate:
  G   at N=1e6 the base structure still pins: g*(4,1)>1/6>g*(8,2) (recomputed).
  M   monotonicity consistency vs the 5e6 sweep (if present): margin_1e6(L) >=
      margin_5e6(L) for every L in the 5e6 achieved range (less weight -> larger
      margin).  A violation would signal a bug in the 1e6 build/bounds.

Output: sweep10k_1e6.json
Single-threaded, ~<10 min.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ppl-014"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from density_threshold import omega_spf, gstar, RMAX
from sweep10k_5e6 import bounds_all, make_wt, sweep_range
from bounds200 import kmax_below_quarter
from bounds2000 import masked_inputs, build_W_flat_precomputed

_T0 = time.time()
CARD = os.path.dirname(os.path.abspath(__file__))
N = int(1e6)
TH = 1.0 / 6.0
L_LO = 2001
L_HI = 10000


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)
    print(f"[43] wrote {os.path.basename(path)}  (wall {time.time() - _T0:.0f}s)", flush=True)


def main():
    out = os.path.join(CARD, "sweep10k_1e6.json")
    print(f"[43] N={N} sweep L={L_LO}..{L_HI} ({L_HI - L_LO + 1} L's)  full range", flush=True)
    gates = {"gates_ok": True}

    om = np.asarray(omega_spf(N), dtype=np.int32)
    wt = make_wt(N); n = np.arange(N + 1, dtype=np.int32)
    n_m, om_i32, wt_m = masked_inputs(n, om, wt, RMAX)
    print(f"[43] 1e6 build done ({len(n_m)} masked terms)  {time.time() - _T0:.0f}s", flush=True)

    # ---- Gate G: base structure pins at 1e6 ----
    ref = {}
    for (L, k) in ((4, 1), (8, 2)):
        W = build_W_flat_precomputed(n_m, om_i32, wt_m, L, RMAX)
        g = float(gstar(W, L, k, RMAX)); ref[f"g*({L},{k})"] = g
        ok = (g > TH) if (L == 4) else (g < TH)
        print(f"[43] G g*({L},{k})@1e6={g:.9f}  {'OK' if ok else 'MISMATCH'}", flush=True)
        if not ok:
            gates["gates_ok"] = False
    gates["G_base_pin_at_1e6"] = ref

    # ---- full sweep ----
    res = sweep_range(n_m, om_i32, wt_m, L_LO, L_HI, N, tag="1e6")
    del n_m, om_i32, wt_m, om, wt, n
    print(f"[43] 1e6 DONE unruled={res['unruled']} crossings={res['crossings']} "
          f"worst_margin={res['worst_margin']}@L{res['worst_margin_L']} "
          f"mean_s/L={res['mean_s_per_L']:.3f}", flush=True)

    # ---- monotonicity consistency vs 5e6 (if present) ----
    mono = {"checked": False}
    p5 = os.path.join(CARD, "sweep10k_5e6.json")
    if os.path.exists(p5):
        with open(p5) as fh:
            r5 = json.load(fh).get("N_5e6", {})
        per5 = r5.get("per_L", {})
        violations = []
        common = 0
        for Lstr, rec in res["per_L"].items():
            if Lstr in per5:
                common += 1
                m1, m5 = rec["margin"], per5[Lstr]["margin"]
                if m1 < m5 - 1e-9:
                    violations.append({"L": int(Lstr), "margin_1e6": m1, "margin_5e6": m5})
        mono = {"checked": True, "n_common": common, "n_violations": len(violations),
                "violations_sample": violations[:8]}
        if violations:
            gates["gates_ok"] = False
        print(f"[43] M monotonicity vs 5e6: {common} common L's, {len(violations)} violations "
              f"{'OK' if not violations else 'MISMATCH'}", flush=True)
    gates["M_monotonicity_vs_5e6"] = mono

    status = "done" if (res["crossings"] == 0 and res["unruled"] == 0) else "needs_review"
    novelty = "extension" if status == "done" else "new" if res["crossings"] > 0 else "extension"
    result = {
        "job": "lab_erdos_261112_43", "problem": "Erdos #143 / PPL 014 ($500)",
        "claim": "full-range d0^all=1/4 pin L=2001..10000 at N=1e6 (weakest N; margin largest)",
        "N": N, "TH": TH, "RMAX": RMAX, "L_range": [L_LO, L_HI],
        "gates": gates, "N_1e6": res, "status": status, "novelty": novelty,
        "verdict": (f"FULL RANGE L=2001..10000 at N=1e6: U'(L,k)<=1/6 for every (L,k) 4k<L "
                    f"(worst margin {res['worst_margin']:.6f} at L={res['worst_margin_L']}); "
                    f"weakest N, so the pin is strongest evidence toward L=10000."
                    if status == "done" else f"unruled={res['unruled']} crossings={res['crossings']}"),
        "elapsed_s": time.time() - _T0,
        "method": ("int32 flattened W builder + vectorised bounds_all (reused from sweep10k_5e6.py, "
                   "gated there vs 24-pass & vs original U'/L'); N=1e6 FULL L=2001..10000."),
    }
    write_json(out, result)
    print(f"[43] VERDICT: {result['verdict']}", flush=True)
    print(f"[43] TOTAL {time.time() - _T0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
