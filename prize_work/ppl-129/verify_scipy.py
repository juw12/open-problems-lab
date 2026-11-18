"""Job lab_erdos_261111_4 - INDEPENDENT VERIFIER.

Second optimizer: scipy.optimize.least_squares doing a JOINT nonlinear least-squares
fit of  y_n = f(n)/(n*sqrt(ln n))  =  C + D*(ln ln n)^(-p)  over (C, D, p) directly
(NO grid, NO per-p linear OLS) - a genuinely different algorithm from the grid+OLS
primary fit, same objective (ratio-space SSE), same 27-point grid.

Two configurations are reported:
  A (card-compliant): p constrained to the card's stated range [0.2, 1.5]. This is the
     apples-to-apples check the card asks for: does a second optimizer recover the same
     (C, D, p) as the grid+OLS fit? PASS = |C_A - C_ols| <= 0.02.
  B (unconstrained p): p in [0.02, 3.0]. Diagnostic only - it exposes the flat valley /
     degenerate direction (p->0, C->-inf, D->+inf) that makes C unidentifiable.
Multi-start on both to avoid local minima; the lowest-cost solution is kept.
"""
import math, json
import numpy as np
from scipy.optimize import least_squares

data = {int(k): v for k, v in json.load(open("fit_data_v2.json")).items()}
pts = sorted(data.items())
N = len(pts)
ys = np.array([f / (n * math.sqrt(math.log(n))) for n, f in pts])
lnln = np.array([math.log(math.log(n)) for n, _ in pts])
TSS = float(np.sum((ys - ys.mean()) ** 2))


def resid(theta):
    C, D, p = theta
    return (C + D * lnln ** (-p)) - ys


def multistart(plo, phi):
    best = None
    for p0 in (0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.5):
        for C0 in (-0.5, 0.35, 0.5, 0.6, 0.7, 0.75):
            for D0 in (0.1, 0.3, 0.5, 0.7):
                p0c = min(max(p0, plo), phi)
                try:
                    r = least_squares(resid, x0=[C0, D0, p0c],
                                      bounds=([-1.0, -3.0, plo], [1.6, 3.0, phi]),
                                      method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                                      max_nfev=40000)
                except Exception:
                    continue
                if best is None or r.cost < best.cost:
                    best = r
    C, D, p = [float(v) for v in best.x]
    sse = 2.0 * best.cost
    r2 = 1.0 - sse / TSS
    aic = N * math.log(sse / N) + 4.0
    bic = N * math.log(sse / N) + 2.0 * math.log(N)
    return {"C": C, "D": D, "p": p, "r2": r2, "sse": sse, "aic": aic, "bic": bic}


A = multistart(0.2, 1.5)      # card-compliant
B = multistart(0.02, 3.0)     # diagnostic

prim = json.load(open("free_p_scan.json"))
C_ols = prim["best_by_r2"]["C"]
p_ols = prim["best_by_r2"]["p"]
r2_ols = prim["best_by_r2"]["r2"]
aic_ols = prim["best_by_r2"]["aic"]
sse_ols = (1.0 - r2_ols) * TSS

agree_A = abs(A["C"] - C_ols) <= 0.02
rel_sse_A = abs(A["sse"] - sse_ols) / sse_ols
sse_ok_A = A["sse"] <= sse_ols * (1.0 + 1e-3)   # OLS best is a near-exhaustive lower bound

out = {
    "n_points": N,
    "A_card_compliant_p_0.2_1.5": {k: round(v, 6) for k, v in A.items()},
    "B_unconstrained_p_0.02_3.0": {k: round(v, 6) for k, v in B.items()},
    "primary_ols": {"C": C_ols, "p": p_ols, "r2": r2_ols, "aic": aic_ols,
                    "sse": round(sse_ols, 10)},
    "A_delta_C_vs_ols": round(abs(A["C"] - C_ols), 6),
    "A_delta_p_vs_ols": round(abs(A["p"] - p_ols), 6),
    "A_C_agrees_within_0.02": bool(agree_A),
    "A_sse_rel_diff": round(rel_sse_A, 8),
    "A_sse_matches_ols_best": bool(sse_ok_A),
    "diagnosis": ("unconstrained fit (B) runs to p=%.3f, C=%.3f, D=%.3f - the p->0/"
                  "C->-inf degenerate valley; C is not identifiable, two correct "
                  "optimizers land at different C with ~equal SSE" % (B["p"], B["C"], B["D"])),
}
with open("verify_scipy.json", "w") as fh:
    json.dump(out, fh, indent=2)

print(f"A (p in [0.2,1.5]): C={A['C']:.5f} D={A['D']:.5f} p={A['p']:.5f} R2={A['r2']:.6f} SSE={A['sse']:.10f}")
print(f"B (p in [0.02,3.0]):C={B['C']:.5f} D={B['D']:.5f} p={B['p']:.5f} R2={B['r2']:.6f} SSE={B['sse']:.10f}")
print(f"primary grid+OLS    : C={C_ols:.5f} p={p_ols:.3f} R2={r2_ols:.6f} SSE={sse_ols:.10f}")
print(f"A |C-C_ols|={abs(A['C']-C_ols):.6f} (<=0.02? {agree_A}); A SSE rel-diff={rel_sse_A:.3e} (matches/beat? {sse_ok_A})")
print("wrote verify_scipy.json")
