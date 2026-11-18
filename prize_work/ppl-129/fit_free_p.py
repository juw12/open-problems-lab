"""Job lab_erdos_261111_4.

Free decay-exponent fit of  f(n) = n*sqrt(ln n) * (C + D/(ln ln n)^p)
on the existing 27-point exact-f grid (fit_data_v2.json, job-36 data).

Idea card: no new f(n) computation. Fit the FREE-p form by OLS over a p-grid in
[0.2, 1.5]; report the best (C, D, p), the full p-scan, and whether C can reach the
claimed Turturean constant sqrt(2)/e = 0.5203 (unverified 2026-09-23 claim).

For fixed p the model is exactly LINEAR in (C, D):
    y_n = f(n)/(n*sqrt(ln n))  =  C + D * x_n(p),   x_n(p) = 1/(ln ln n)^p
so each p is a plain 2-parameter OLS. The p-scan is then a 1-D search over the
nonlinear parameter p. AIC/BIC use the job-36 convention  AIC = N*ln(SSE/N) + 2k.

GATE (code correctness, not the result): OLS at p=0.5 must reproduce job-36 M1
    A = 0.41581,  B = 0.39749,  R2 = 0.7062,  AIC = -246.86
to 4-5 digits. If the gate passes, the OLS implementation is the same one job-36 used.

KILL TEST (pre-registered, from the card):  best R2 < 0.72,  OR  C shifts > 0.05
across adjacent p values  ->  'correction form unidentifiable on 27 points'.
'Adjacent' is operationalized three ways and all are reported: a 0.1 step (natural
exponent step), a 0.05 step, and the credible-AIC region {AIC(p)-min AIC <= 2}.
"""
import math, json
import numpy as np

data = {int(k): v for k, v in json.load(open("fit_data_v2.json")).items()}
pts = sorted(data.items())
N = len(pts)
SQRT2_OVER_E = math.sqrt(2.0) / math.e          # 0.52029822...
assert N == 27, f"expected 27 points, got {N}"

ys = np.array([f / (n * math.sqrt(math.log(n))) for n, f in pts])       # ratio space
lnln = np.array([math.log(math.log(n)) for n, _ in pts])                # ln ln n
TSS = float(np.sum((ys - ys.mean()) ** 2))


def ols(p):
    """Exact 2-param OLS of y on [1, (ln ln n)^-p]. Returns (C, D, r2, aic, bic, sse)."""
    x = lnln ** (-p)
    X = np.column_stack([np.ones(N), x])
    (C, D), _res, _rank, _sv = np.linalg.lstsq(X, ys, rcond=None)
    sse = float(np.sum((ys - (C + D * x)) ** 2))
    r2 = 1.0 - sse / TSS
    aic = N * math.log(sse / N) + 4.0           # k = 2
    bic = N * math.log(sse / N) + 2.0 * math.log(N)
    return float(C), float(D), r2, aic, bic, sse


# ---------------------------------------------------------------- GATE
gC, gD, gR2, gAIC, _gBIC, _gsse = ols(0.5)
gate = (abs(gC - 0.41581) < 2e-4 and abs(gD - 0.39749) < 2e-4
        and abs(gR2 - 0.7062) < 5e-4 and abs(gAIC - (-246.86)) < 0.20)
print(f"GATE p=0.5: C={gC:.5f} D={gD:.5f} R2={gR2:.5f} AIC={gAIC:.2f}"
      f"   [job36: A=0.41581 B=0.39749 R2=0.7062 AIC=-246.86]  -> {'PASS' if gate else 'FAIL'}")
if not gate:
    raise SystemExit("GATE FAILED: OLS does not reproduce job-36 M1; code is wrong, abort.")

# ---------------------------------------------------------------- p-scan
pmin, pmax, step = 0.2, 1.5, 0.005
ps = np.round(np.arange(pmin, pmax + step / 2, step), 6)
scan = []
for p in ps:
    C, D, r2, aic, bic, sse = ols(float(p))
    scan.append({"p": float(p), "C": C, "D": D, "r2": r2, "aic": aic, "bic": bic, "sse": sse})

best_r2  = max(scan, key=lambda s: s["r2"])
best_aic = min(scan, key=lambda s: s["aic"])
best_bic = min(scan, key=lambda s: s["bic"])

Ccurve = {s["p"]: s["C"] for s in scan}
AIC_min = best_aic["aic"]


def c_at(p):
    p = min(max(p, pmin), pmax)
    p = round(p, 6)
    return Ccurve.get(p, ols(p)[0])


# ---------------- sensitivity of C near the R2-best p (operationalize the kill)
pbest = best_r2["p"]
Cbest = best_r2["C"]

# (a) 0.1-step "adjacent" shift, and (b) 0.05-step, over the whole grid
shift_010 = max(abs(scan[i + 1]["C"] - scan[i]["C"]) for i in range(0, len(scan) - 1, 20))  # 0.1 apart
shift_005 = max(abs(scan[i + 1]["C"] - scan[i]["C"]) for i in range(0, len(scan) - 1, 10))  # 0.05 apart
shift_005_grid = max(abs(scan[i + 1]["C"] - scan[i]["C"]) for i in range(len(scan) - 1))     # true 0.005 grid step

# (c) credible-AIC region: { p : AIC(p) - minAIC <= 2 }
cred = [s for s in scan if s["aic"] - AIC_min <= 2.0]
cred_p = [s["p"] for s in cred]
cred_C = [s["C"] for s in cred]
cred_range = (min(cred_p), max(cred_p), min(cred_C), max(cred_C))
cred_Cspan = max(cred_C) - min(cred_C)

# full-grid C range
full_Crange = max(s["C"] for s in scan) - min(s["C"] for s in scan)

# ---------------- can C reach the Turturean constant?
target = SQRT2_OVER_E
hits = [s["p"] for s in scan if abs(s["C"] - target) < 1e-9]
# nearest p to target C
nearest = min(scan, key=lambda s: abs(s["C"] - target))
# C at the best-p and best-AIC p
c_reach_at_bestR2 = Cbest
c_reach_at_bestAIC = best_aic["C"]
# min |C(p) - target| over the credible region
min_gap_cred = min(abs(c - target) for c in cred_C) if cred_C else float("nan")

# ---------------- kill test (pre-registered)
kill_r2 = best_r2["r2"] < 0.72
kill_shift = shift_010 > 0.05
kill = kill_r2 or kill_shift

result = {
    "n_points": N,
    "gate": {"p0.5_C": round(gC, 5), "p0.5_D": round(gD, 5), "p0.5_R2": round(gR2, 5),
             "p0.5_AIC": round(gAIC, 2), "pass": bool(gate)},
    "turturean_target_C": round(target, 6),
    "best_by_r2": {k: round(v, 6) for k, v in best_r2.items()},
    "best_by_aic": {k: round(v, 6) for k, v in best_aic.items()},
    "best_by_bic": {k: round(v, 6) for k, v in best_bic.items()},
    "sensitivity": {
        "pbest_R2": pbest, "Cbest_R2": round(Cbest, 6),
        "shift_C_0.10_window": round(shift_010, 5),
        "shift_C_0.05_window": round(shift_005, 5),
        "shift_C_grid_0.005": round(shift_005_grid, 5),
        "credible_AIC_le2_p_range": [round(cred_p[0], 4), round(cred_p[-1], 4)],
        "credible_AIC_le2_C_range": [round(min(cred_C), 5), round(max(cred_C), 5)],
        "credible_AIC_le2_C_span": round(cred_Cspan, 5),
        "full_grid_C_range": round(full_Crange, 5),
        "central_dC_dp_at_pbest": round((c_at(pbest + 0.05) - c_at(pbest - 0.05)) / 0.10, 4),
    },
    "c_reaches_0.5203": {
        "C_at_bestR2_p": round(c_reach_at_bestR2, 6),
        "C_at_bestAIC_p": round(c_reach_at_bestAIC, 6),
        "p_where_C_nearest_target": nearest["p"],
        "min_|C_target|_over_grid": round(abs(nearest["C"] - target), 6),
        "min_|C_target|_over_credible_region": round(min_gap_cred, 6),
        "abs_gap_at_bestR2_p": round(abs(c_reach_at_bestR2 - target), 6),
    },
    "kill_test": {
        "best_R2": round(best_r2["r2"], 5),
        "kill_R2_lt_0.72": bool(kill_r2),
        "shift_C_0.10_gt_0.05": bool(kill_shift),
        "kill_fires": bool(kill),
        "verdict": ("correction form unidentifiable on 27 points" if kill
                    else "correction form identifiable on 27 points"),
    },
    "scan": scan,
}
with open("free_p_scan.json", "w") as fh:
    json.dump(result, fh, indent=2)

# ---------------- human-readable log
with open("fit_free_p.log", "w") as fh:
    fh.write(f"Free decay-exponent fit  f(n)=n*sqrt(ln n)*(C + D/(ln ln n)^p)  "
             f"on {N}-point grid (n=100..5000)\n")
    fh.write(f"Turturean target C = sqrt(2)/e = {target:.6f}\n\n")
    fh.write(f"GATE p=0.5: C={gC:.5f} D={gD:.5f} R2={gR2:.5f} AIC={gAIC:.2f} "
             f"{'PASS' if gate else 'FAIL'}\n\n")
    fh.write(f"Best by R2 : p={best_r2['p']:.3f}  C={best_r2['C']:.5f}  D={best_r2['D']:.5f}  "
             f"R2={best_r2['r2']:.5f}  AIC={best_r2['aic']:.2f}\n")
    fh.write(f"Best by AIC: p={best_aic['p']:.3f}  C={best_aic['C']:.5f}  D={best_aic['D']:.5f}  "
             f"R2={best_aic['r2']:.5f}  AIC={best_aic['aic']:.2f}\n")
    fh.write(f"Best by BIC: p={best_bic['p']:.3f}  C={best_bic['C']:.5f}  D={best_bic['D']:.5f}  "
             f"R2={best_bic['r2']:.5f}  BIC={best_bic['bic']:.2f}\n\n")
    fh.write(f"Sensitivity of C to p:\n")
    fh.write(f"  C at R2-best p={pbest:.3f}: {Cbest:.5f}\n")
    fh.write(f"  max C-shift over 0.10 p-window : {shift_010:.5f}\n")
    fh.write(f"  max C-shift over 0.05 p-window : {shift_005:.5f}\n")
    fh.write(f"  max C-shift over 0.005 grid step: {shift_005_grid:.5f}\n")
    fh.write(f"  central dC/dp at pbest          : "
             f"{(c_at(pbest+0.05)-c_at(pbest-0.05))/0.10:.4f}\n")
    fh.write(f"  credible region (AIC-min<=2): p in [{cred_p[0]:.3f},{cred_p[-1]:.3f}], "
             f"C in [{min(cred_C):.5f},{max(cred_C):.5f}]  (span {cred_Cspan:.5f})\n")
    fh.write(f"  full-grid C range              : {full_Crange:.5f}\n\n")
    fh.write(f"Can C reach the Turturean constant {target:.5f}?\n")
    fh.write(f"  C at best-R2 p                 : {c_reach_at_bestR2:.5f}  "
             f"(gap {abs(c_reach_at_bestR2-target):.5f})\n")
    fh.write(f"  C at best-AIC p                : {c_reach_at_bestAIC:.5f}\n")
    fh.write(f"  nearest C on grid to target    : p={nearest['p']:.3f}, |C-target|="
             f"{abs(nearest['C']-target):.6f}\n")
    fh.write(f"  min |C-target| in credible region: {min_gap_cred:.6f}\n\n")
    fh.write(f"KILL TEST: best R2={best_r2['r2']:.5f} (<0.72? {kill_r2}); "
             f"max 0.10 C-shift={shift_010:.5f} (>0.05? {kill_shift})\n")
    fh.write(f"  -> {'KILL FIRES: ' if kill else 'kill does not fire: '}{result['kill_test']['verdict']}\n")
    fh.write(f"\nFull p-scan (p, C, D, R2, AIC):\n")
    fh.write(f"{'p':>7} {'C':>10} {'D':>10} {'R2':>9} {'AIC':>9}\n")
    for s in scan:
        fh.write(f"{s['p']:7.3f} {s['C']:10.5f} {s['D']:10.5f} {s['r2']:9.5f} {s['aic']:9.2f}\n")

print(f"Best R2 : p={best_r2['p']:.3f} C={best_r2['C']:.5f} D={best_r2['D']:.5f} R2={best_r2['r2']:.5f}")
print(f"Best AIC: p={best_aic['p']:.3f} C={best_aic['C']:.5f} D={best_aic['D']:.5f} AIC={best_aic['aic']:.2f}")
print(f"Best BIC: p={best_bic['p']:.3f} C={best_bic['C']:.5f} D={best_bic['D']:.5f} BIC={best_bic['bic']:.2f}")
print(f"C-sens: shift0.1={shift_010:.4f} shift0.05={shift_005:.4f} credCspan={cred_Cspan:.4f} "
      f"fullCspan={full_Crange:.4f}")
print(f"C reach target {target:.5f}: at bestR2 p gap={abs(c_reach_at_bestR2-target):.5f}; "
      f"min gap grid={abs(nearest['C']-target):.6f}; min gap cred={min_gap_cred:.6f}")
print(f"KILL: R2<0.72={kill_r2}  shift0.1>0.05={kill_shift}  -> fires={kill}: {result['kill_test']['verdict']}")
print("wrote free_p_scan.json and fit_free_p.log")
