"""Fast fit + analysis. Reads fit_data_v2.json, writes fit_results.json. ~seconds."""
import math, json
import numpy as np

out = {int(k): v for k, v in json.load(open("fit_data_v2.json")).items()}
pts = sorted(out.items())

def ratio(n, f):
    return f / (n * math.sqrt(math.log(n)))

xs = [1.0 / math.sqrt(math.log(math.log(n))) for n, _ in pts]
ys = [ratio(n, f) for n, f in pts]
N = len(pts)

X = np.column_stack([np.ones(N), np.array(xs)])
(A, B), res, rank, sv = np.linalg.lstsq(X, np.array(ys), rcond=None)
pred = A + B * np.array(xs)
ss_res = float(np.sum((np.array(ys) - pred) ** 2))
ss_tot = float(np.sum((np.array(ys) - np.mean(ys)) ** 2))
r2 = 1 - ss_res / ss_tot

# independent route: hand-rolled normal equations
Sx = sum(xs); Sy = sum(ys); Sxx = sum(x * x for x in xs); Sxy = sum(x * y for x, y in zip(xs, ys))
det = N * Sxx - Sx * Sx
Bh = (N * Sxy - Sx * Sy) / det
Ah = (Sy - Bh * Sx) / N

# model selection
C = Sy / N
ss_a = sum((y - C) ** 2 for y in ys)
D = Sxy / Sxx
ss_b = sum((y - D * x) ** 2 for x, y in zip(xs, ys))
aic_a = N * math.log(ss_a / N) + 2
aic_b = N * math.log(ss_b / N) + 2
aic_c = N * math.log(ss_res / N) + 4

nonmono = sum(1 for i in range(N - 1) if ys[i + 1] > ys[i] + 1e-9)
SQRT2_OVER_E = math.sqrt(2) / math.e
TWO_OVER_SQRTE = 2 / math.sqrt(math.e)

# residual detail (a few anchor points)
anchors = {}
for n, f in pts:
    if n in (100, 500, 1000, 2000, 3000, 5000):
        r = ratio(n, f)
        x = 1.0 / math.sqrt(math.log(math.log(n)))
        anchors[str(n)] = {"f": f, "ratio": round(r, 5), "model": round(float(A + B * x), 5),
                           "resid": round(r - float(A + B * x), 5)}

result = {
    "n_points": N,
    "x_range": [round(min(xs), 5), round(max(xs), 5)],
    "A_numpy": round(float(A), 5), "B_numpy": round(float(B), 5), "r2": round(r2, 5),
    "A_hand": round(float(Ah), 5), "B_hand": round(float(Bh), 5),
    "turturean_C": round(float(C), 5), "eplb_D": round(float(D), 5),
    "aic": {"turturean_const": round(aic_a, 2), "eplb_origin": round(aic_b, 2), "combined": round(aic_c, 2)},
    "r2_origin_model": round(1 - ss_b / ss_tot, 5),
    "num_nonmonotone_increases": int(nonmono),
    "A_in_[0.50,0.75]": bool(0.50 <= float(A) <= 0.75),
    "ratio_at_5000": round(ratio(5000, out[5000]), 5),
    "ratio_at_1000": round(ratio(1000, out[1000]), 5),
    "sqrt2_over_e": round(SQRT2_OVER_E, 5),
    "two_over_sqrt_e": round(TWO_OVER_SQRTE, 5),
    "anchors": anchors,
    "data": {str(k): v for k, v in out.items()},
}
with open("fit_results.json", "w") as fh:
    json.dump(result, fh, indent=2)

print(f"A={A:.5f} B={B:.5f} R2={r2:.5f}   (hand: A={Ah:.5f} B={Bh:.5f})")
print(f"Turturean C={C:.5f} | EP-lb D={D:.5f} R2_origin={1-ss_b/ss_tot:.4f}")
print(f"AIC const={aic_a:.2f} origin={aic_b:.2f} combined={aic_c:.2f}")
print(f"nonmonotone_increases={nonmono}  A in [0.5,0.75]={0.50<=float(A)<=0.75}")
print(f"ratio(1000)={ratio(1000,out[1000]):.5f} ratio(5000)={ratio(5000,out[5000]):.5f}")
print(f"target A sqrt(2)/e={SQRT2_OVER_E:.5f}; target B 2/sqrt(e)={TWO_OVER_SQRTE:.5f}")
