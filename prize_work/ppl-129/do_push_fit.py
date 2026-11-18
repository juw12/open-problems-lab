#!/usr/bin/env python3
"""Combined fit after pushing f(n) to n=10000/15000/20000 (job ..._8).
Dataset = 27-point grid (fit_data_v2.json, job 36) + 3 new anchors (push_data.json) = 30 pts.
Fit model M1:  ratio(n) = f(n)/(n sqrt(log n)) = A + B/sqrt(loglog n).
Reports A, B, R^2, per-point residuals, and both convergence ratios at the new anchors.
Writes push_results.json.
"""
import json, os, math
D = os.path.dirname(os.path.abspath(__file__))
grid = {int(k): v for k, v in json.load(open(os.path.join(D, "fit_data_v2.json"))).items()}
push = json.load(open(os.path.join(D, "push_data.json")))
for k, v in push.items():
    grid[int(k)] = v["f"]

pts = sorted(grid.items())
def ratio(n, f): return f / (n * math.sqrt(math.log(n)))
xs = [1.0 / math.sqrt(math.log(math.log(n))) for n, _ in pts]
ys = [ratio(n, f) for n, f in pts]
N = len(pts)
Sx, Sy = sum(xs), sum(ys)
Sxx, Sxy = sum(x*x for x in xs), sum(x*y for x, y in zip(xs, ys))
det = N*Sxx - Sx*Sx
B = (N*Sxy - Sx*Sy) / det
A = (Sy - B*Sx) / N
pred = [A + B*x for x in xs]
ss_res = sum((y-p)**2 for y, p in zip(ys, pred))
ss_tot = sum((y - Sy/N)**2 for y in ys)
r2 = 1 - ss_res/ss_tot
import numpy as np
X = np.column_stack([np.ones(N), np.array(xs)])
(A_np, B_np), *_ = np.linalg.lstsq(X, np.array(ys), rcond=None)

resid = {str(n): round(ratio(n, f) - (A + B*(1.0/math.sqrt(math.log(math.log(n))))), 6) for n, f in pts}

SQRT2_E = math.sqrt(2)/math.e
TWO_SQRT_E = 2/math.sqrt(math.e)
anchors = {}
for k, v in push.items():
    n = int(k); f = v["f"]
    anchors[str(n)] = dict(f=f, dinic_ok=v.get("dinic_ok"),
        ratio_t=round(ratio(n, f), 6),
        ratio_ep=round(f*math.log(math.log(n))/(TWO_SQRT_E*n*math.sqrt(math.log(n))), 6))

out = dict(
    n_points=N, x_range=[round(min(xs), 5), round(max(xs), 5)],
    A=round(A, 5), B=round(B, 5), r2=round(r2, 5),
    A_numpy=round(float(A_np), 5), B_numpy=round(float(B_np), 5),
    sqrt2_over_e=round(SQRT2_E, 5), two_over_sqrt_e=round(TWO_SQRT_E, 5),
    anchors=anchors, residuals_M1=resid,
    ratios_all={str(n): round(ratio(n, f), 5) for n, f in pts},
)
json.dump(out, open(os.path.join(D, "push_results.json"), "w"), indent=2)
print(json.dumps(out, indent=2))
