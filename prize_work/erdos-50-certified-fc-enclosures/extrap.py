# -*- coding: utf-8 -*-
"""Job lab_erdos_261113_43, option 2 analysis: convergence law for f_N(k/100).

Data: grid_N100000000.json — exact counts C_k(N_end) for 20 rungs (5e6..1e8,
step 5e6), 99 points k/100, C_k = #{n<=N : 100*phi(n)//n <= k-1}.

Per k: fit f_N = f0 + A/ln(N) on rungs [1e7 .. 5e7] (9 rungs), predict
rungs [5.5e7 .. 1e8] (10 rungs), RMSE = honest error bar. Compare with the
f0 + B/sqrt(N) law. Nesting check: monotone decrease in N.
No published model is assumed; this is an empirical law, verified by the
held-out half of the rungs (the data was computed, not generated).
"""
import json
import time
import numpy as np

t0 = time.time()
HERE = __file__
g = json.load(open(HERE.replace("extrap.py", "grid_N100000000.json"), encoding="utf-8"))
rungs = g["counts_per_rung"]                 # [[N_end, [C_1..C_99]]] x20
Nends = np.array([r[0] for r in rungs], float)
C = np.array([r[1] for r in rungs], "int64") # 20 x 99
assert C.shape == (20, 99)
F = C / Nends[:, None].astype("float64")     # f_N(k/100)

# fit rungs: N in [1e7, 5e7] -> indices 1..9 ; predict: indices 10..19
fit_idx = np.arange(1, 10)
pred_idx = np.arange(10, 20)
lnN = np.log(Nends)

per_k = []
for k in range(1, 100):
    y = F[:, k - 1]
    # 1/ln N law
    Xf = 1.0 / lnN[fit_idx]
    Mf = np.column_stack([np.ones_like(Xf), Xf])
    b_ln, *_ = np.linalg.lstsq(Mf, y[fit_idx], rcond=None)
    res_ln = Mf @ b_ln - y[fit_idx]
    rmse_ln = float(np.sqrt(np.mean(res_ln ** 2)))
    pred_ln = Mf @ b_ln
    # sqrt N law
    Xs = 1.0 / np.sqrt(Nends[fit_idx])
    Ms = np.column_stack([np.ones_like(Xs), Xs])
    b_sq, *_ = np.linalg.lstsq(Ms, y[fit_idx], rcond=None)
    rmse_sq = float(np.sqrt(np.mean((Ms @ b_sq - y[fit_idx]) ** 2)))
    # hold-out error at the last predict rung (N=1e8), both laws
    p_ln_1e8 = float(b_ln[0] + b_ln[1] / lnN[19])
    p_sq_1e8 = float(b_sq[0] + b_sq[1] / np.sqrt(Nends[19]))
    a1e8 = float(y[19])
    # monotone decrease in N over ALL rungs (nesting)
    mono_full = bool(np.all(np.diff(y) < 0))
    mono_fitA = bool(b_ln[1] < 0)  # 1/ln law decreasing
    # hold-out band: |pred - actual| at every predicted rung (1/ln law)
    hold_errs = [abs(float(b_ln[0] + b_ln[1] / lnN[i]) - y[i]) for i in pred_idx]
    within = sum(1 for e in hold_errs if e <= 2.0 * rmse_ln)
    per_k.append({
        "k": k, "f_1e7": float(y[1]), "f_1e8": a1e8,
        "delta_1e7_1e8": a1e8 - float(y[1]),
        "f0_ln": float(b_ln[0]), "A_ln": b_ln[1], "rmse_ln": rmse_ln,
        "rmse_sqrt": rmse_sq,
        "pred_1e8_ln": p_ln_1e8, "pred_1e8_sqrt": p_sq_1e8,
        "abs_err_1e8_ln": abs(p_ln_1e8 - a1e8), "abs_err_1e8_sqrt": abs(p_sq_1e8 - a1e8),
        "holdout_within_2rmse": f"{within}/10",
        "mono_decrease_full": mono_full, "fit_A_negative": mono_fitA,
    })

n_ln_better = sum(1 for r in per_k if r["rmse_ln"] < r["rmse_sqrt"])
n_mono = sum(1 for r in per_k if r["mono_decrease_full"])
n_aneg = sum(1 for r in per_k if r["fit_A_negative"])
# headline point k=50 (c=1/2)
h50 = per_k[49]
out = {
    "N": 100000000, "fit_window": [int(Nends[1]), int(Nends[9])],
    "pred_window": [int(Nends[10]), int(Nends[19])],
    "law": "f_N = f0 + A/ln(N) vs f0 + B/sqrt(N); RMSE on fit window = error bar",
    "rmse_ln_better_count": f"{n_ln_better}/99",
    "mono_decrease_full_count": f"{n_mono}/99",
    "fit_A_negative_count": f"{n_aneg}/99",
    "c_half": {key: h50[key] for key in
               ("f_1e7", "f_1e8", "delta_1e7_1e8", "f0_ln", "A_ln", "rmse_ln",
                "rmse_sqrt", "pred_1e8_ln", "pred_1e8_sqrt",
                "abs_err_1e8_ln", "abs_err_1e8_sqrt", "holdout_within_2rmse",
                "mono_decrease_full", "fit_A_negative")},
    "per_k": per_k,
    "elapsed_s": round(time.time() - t0, 2),
}
with open(HERE.replace("extrap.py", "extrap_N100000000.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1)
print(f"[fit] 1/lnN beats sqrt(N) RMSE: {n_ln_better}/99 ; mono decrease in N: {n_mono}/99 ; A<0: {n_aneg}/99", flush=True)
print(f"[c=1/2] f(1e7)={h50['f_1e7']:.7f} f(1e8)={h50['f_1e8']:.7f} delta={h50['delta_1e7_1e8']:.2e} "
      f"f0={h50['f0_ln']:.7f} A={h50['A_ln']:.3e} rmse_ln={h50['rmse_ln']:.2e} rmse_sqrt={h50['rmse_sqrt']:.2e} "
      f"pred_1e8={h50['pred_1e8_ln']:.7f} err={h50['abs_err_1e8_ln']:.2e} ({h50['holdout_within_2rmse']} holdout)", flush=True)
worst = max(per_k, key=lambda r: r["abs_err_1e8_ln"] / max(r["rmse_ln"], 1e-300))
print(f"[worst] k={worst['k']} pred_err={worst['abs_err_1e8_ln']:.2e} rmse={worst['rmse_ln']:.2e} "
      f"ratio={worst['abs_err_1e8_ln']/max(worst['rmse_ln'],1e-300):.1f}", flush=True)
print(f"[done] {time.time()-t0:.1f}s", flush=True)
