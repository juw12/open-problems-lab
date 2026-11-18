#!/usr/bin/env python3
"""Job lab_erdos_261111_3: re-fit + curvature test on the extended f(n) frontier.

Loads the full verified exact dataset:
    fit_data_v2.json  : 27-point grid n=100..5000   (job 36, HK+Dinic verified)
    push_data.json    : 10000 / 15000 / 20000        (job _8,  HK+Dinic verified)
    push50k_data.json : 30000 / 40000 / 50000        (this job, HK+Dinic verified)

(1) OLS fit  r_t(n) = f/(n*sqrt(log n)) = A + B/sqrt(loglog n)   on all 33 points
    (numpy lstsq + hand-rolled normal equations, must agree to 5 dp).
(2) PER-ANCHOR residuals at the six large anchors.
(3) CURVATURE test on the ratio over the evenly-spaced anchors
    10000,20000,30000,40000,50000: first and second differences of r_t.
    - 2nd difference > 0  -> ratio flattening (decay decelerating)  -> NO acceleration toward 0.5203
    - 2nd difference < 0  -> ratio accelerating downward            -> possible asymptotic onset
(4) EP lower-bound sharpness  r_prop = f / EP_lower(n),  EP_lower(n)=(2/sqrt(e))*n*sqrt(log n/loglog n),
    reported at every anchor (tends to 1 iff the EP bound is asymptotically tight).
Writes fit50k_results.json.
"""
import json, os, math, sys
D = os.path.dirname(os.path.abspath(__file__))

SQRT2_E = math.sqrt(2) / math.e
TWO_SQRT_E = 2 / math.sqrt(math.e)


def rt(n, f):
    return f / (n * math.sqrt(math.log(n)))


def ep_lower(n):
    return TWO_SQRT_E * n * math.sqrt(math.log(n) / math.log(math.log(n)))


def load(name):
    with open(os.path.join(D, name)) as fh:
        return json.load(fh)


def main():
    # ---- assemble full verified dataset {n: f} ----
    grid = {int(k): int(v) for k, v in load("fit_data_v2.json").items()}          # 27 pts
    for name in ("push_data.json", "push50k_data.json"):
        for k, v in load(name).items():
            f = v["f"] if isinstance(v, dict) else v
            if f is not None:
                grid[int(k)] = int(f)
    ns = sorted(grid)
    print(f"Total verified points: {len(ns)}  n in [{ns[0]}, {ns[-1]}]")
    for n in [10000, 15000, 20000, 30000, 40000, 50000]:
        if n not in grid:
            print(f"WARNING: n={n} missing -- compute not finished?"); 
    assert all(n in grid for n in (10000, 15000, 20000, 30000, 40000, 50000)), \
        "need all six large anchors for the curvature test"

    # ---- (1) OLS fit A + B/sqrt(loglog n) ----
    xs = [1.0 / math.sqrt(math.log(math.log(n))) for n in ns]
    ys = [rt(n, grid[n]) for n in ns]
    N = len(xs)
    # hand-rolled normal equations
    Sx = sum(xs); Sy = sum(ys); Sxx = sum(x * x for x in xs); Sxy = sum(x * y for x, y in zip(xs, ys))
    det = N * Sxx - Sx * Sx
    B_h = (N * Sxy - Sx * Sy) / det
    A_h = (Sy - B_h * Sx) / N
    # numpy
    try:
        import numpy as np
        X = np.column_stack([np.ones(N), np.array(xs)])
        (A_np, B_np), _, _, _ = np.linalg.lstsq(X, np.array(ys), rcond=None)
        pred = A_np + B_np * np.array(xs)
        ss_res = float(np.sum((np.array(ys) - pred) ** 2))
        ss_tot = float(np.sum((np.array(ys) - np.mean(ys)) ** 2))
        r2 = 1 - ss_res / ss_tot
    except ImportError:
        A_np = B_np = None
        pred = [A_h + B_h * x for x in xs]
        ss_res = sum((y - p) ** 2 for y, p in zip(ys, pred))
        ss_tot = sum((y - Sy / N) ** 2 for y in ys)
        r2 = 1 - ss_res / ss_tot

    agree = A_np is None or (abs(A_np - A_h) < 5e-5 and abs(B_np - B_h) < 5e-5)
    A = A_np if A_np is not None else A_h
    B = B_np if B_np is not None else B_h

    # ---- (2) per-anchor residuals (observed - fit) ----
    resid = {}
    for n in [1000, 5000, 10000, 15000, 20000, 30000, 40000, 50000]:
        if n in grid:
            x = 1.0 / math.sqrt(math.log(math.log(n)))
            resid[str(n)] = round(rt(n, grid[n]) - (A + B * x), 6)

    # ---- (3) curvature test over evenly-spaced large anchors ----
    anchors = [10000, 20000, 30000, 40000, 50000]
    rt_series = [rt(n, grid[n]) for n in anchors]
    d1 = [rt_series[i + 1] - rt_series[i] for i in range(len(rt_series) - 1)]
    d2 = [d1[i + 1] - d1[i] for i in range(len(d1) - 1)]
    mean_d2 = sum(d2) / len(d2)
    if mean_d2 > 0:
        curvature = "flattening (decay decelerating -> NO acceleration toward 0.5203)"
    else:
        curvature = "accelerating downward (possible asymptotic onset)"

    # ---- (4) EP sharpness at every large anchor ----
    ep = {}
    for n in [10000, 15000, 20000, 30000, 40000, 50000]:
        f = grid[n]
        ep[str(n)] = dict(f=f, ep_lower=round(ep_lower(n), 1),
                          r_prop=round(f / ep_lower(n), 6),
                          r_ep_card=round(f * math.log(math.log(n)) / (TWO_SQRT_E * n * math.sqrt(math.log(n))), 6))

    result = {
        "n_points": N, "n_range": [ns[0], ns[-1]],
        "fit": {"A_numpy": round(float(A_np), 6) if A_np is not None else None,
                "B_numpy": round(float(B_np), 6) if B_np is not None else None,
                "A_hand": round(A_h, 6), "B_hand": round(B_h, 6),
                "A": round(float(A), 6), "B": round(float(B), 6), "R2": round(r2, 6),
                "numpy_hand_agree_5dp": bool(agree)},
        "turturean_constant_sqrt2_over_e": round(SQRT2_E, 6),
        "A_vs_turturean": round(float(A) - SQRT2_E, 6),
        "residuals_obs_minus_fit": resid,
        "ratio_series_evenly_spaced": {str(n): round(rt(n, grid[n]), 6) for n in anchors},
        "first_diffs": [round(x, 6) for x in d1],
        "second_diffs": [round(x, 6) for x in d2],
        "mean_second_diff": round(mean_d2, 6),
        "curvature_verdict": curvature,
        "ep_sharpness": ep,
        "r_t_at_50000": round(rt(50000, grid[50000]), 6),
        "r_t_min_over_range": round(min(rt(n, grid[n]) for n in anchors), 6),
    }
    with open(os.path.join(D, "fit50k_results.json"), "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"\n=== FIT (A + B/sqrt(loglog n), {N} pts, n={ns[0]}..{ns[-1]}) ===")
    print(f"A={float(A):.6f}  B={float(B):.6f}  R2={r2:.6f}  numpy==hand(5dp)={agree}")
    print(f"A vs Turturean sqrt2/e={SQRT2_E:.6f}:  A - c = {float(A)-SQRT2_E:+.6f}")
    print(f"\nresiduals (obs-fit): {resid}")
    print(f"\n=== CURVATURE (r_t over {anchors}) ===")
    print(f"r_t      : { [round(x,6) for x in rt_series] }")
    print(f"1st diff : { [round(x,6) for x in d1] }")
    print(f"2nd diff : { [round(x,6) for x in d2] }   mean={mean_d2:+.6f}")
    print(f"VERDICT  : {curvature}")
    print(f"\n=== EP lower-bound sharpness r_prop=f/EP_lower (->1 if tight) ===")
    for k in [10000, 15000, 20000, 30000, 40000, 50000]:
        e = ep[str(k)]
        print(f"n={k:6d} f={e['f']:6d} EP_lower={e['ep_lower']:8.0f} r_prop={e['r_prop']:.6f} r_ep_card={e['r_ep_card']:.6f}")
    print(f"\nwrote fit50k_results.json")


if __name__ == "__main__":
    main()
