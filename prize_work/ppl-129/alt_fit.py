"""
Alternative asymptotic fits for #710 (lab_erdos_20261110_36), run after the
card's 1/sqrt(loglog n) model is killed (A out of [0.50, 0.75]).

Models (ratio r(n) = f(n)/(n*sqrt(log n)), linear in the free params):
  M1 card     : r = A + B/sqrt(loglog n)          [the idea card]
  M2 const    : r = A                              [best constant multiple]
  M3 EP-orig  : r = B/sqrt(loglog n)              [A=0, pure EP-lb decay]
  M4 alt-log  : r = C + D/sqrt(log n)             [alternative correction]
  M5 alt-ll   : r = C + D/loglog n                [slower loglog correction]

For each: OLS slope/intercept, R^2, RMSE, and AIC (k = number of free params,
including the noise term). Best AIC = smallest.
"""
import json, math

p = r"E:\Projects\MillenniumResearch\prize_work\ppl-129\fit_data_v2.json"
data = json.load(open(p, encoding="utf-8"))

ns  = [int(k) for k in data]
fs  = [int(data[str(n)]) for n in ns]
N   = len(ns)

import numpy as np

def ratio(n, f): return f / (n * math.sqrt(math.log(n)))
ys = [ratio(n, f) for n, f in zip(ns, fs)]
y  = np.array(ys, dtype=float)

def aic(n, r2, k):
    # OLS with normal errors: AIC = n*ln(RSS/n) + 2k, k = #free params
    rss = (1.0 - r2) * float(np.sum((y - y.mean())**2))
    return n * math.log(rss / n) + 2.0 * k

def fit(X, label, k):
    X = np.asarray(X, dtype=float)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot
    rmse = math.sqrt(ss_res / len(y))
    a = aic(len(y), r2, k)
    print(f"{label:12s} coef={list(round(float(c),5) for c in coef)} R2={r2:.4f} RMSE={rmse:.5f} AIC={a:.2f}")
    return coef, r2, rmse, a

logn    = np.array([math.log(n) for n in ns])
loglogn = np.array([math.log(math.log(n)) for n in ns])
ones    = np.ones(N)

print("n points =", N, " ratio range", round(float(y.min()),5), "-", round(float(y.max()),5))
print("-" * 78)
fit(np.column_stack([ones, 1.0/np.sqrt(loglogn)]), "M1 card A+B/ll", 2)
fit(ones.reshape(-1,1), "M2 const", 1)
fit((1.0/np.sqrt(loglogn)).reshape(-1,1), "M3 EP B/ll", 1)
fit(np.column_stack([ones, 1.0/np.sqrt(logn)]), "M4 C+D/sqrt(log)", 2)
fit(np.column_stack([ones, 1.0/loglogn]), "M5 C+D/loglog", 2)
print("-" * 78)
print("targets: Turturean A=sqrt(2)/e=0.52026 ; EP-lb B=2/sqrt(e)=1.21306")
