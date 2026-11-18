#!/usr/bin/env python3
"""Job lab_erdos_261112_21 - PART 2: local D_2 (correlation dimension) contrast
between the c~0.65 lobe (highest-density interior region) and the c~0.5 lobe.

METHOD. For a window W=[a,b) and scale h, take the CONDITIONAL sample
r_w = { r : a <= r < b } and partition W into h-width sub-bins. The local
Renyi D_2 (correlation dimension) of the conditional measure mu|_W is
  D_2 = log( sum_i m_i^2 ) / log(h),   m_i = (count in bin i)/(total in W),
computed with the job-261112_8 renyi_D (mass normalized WITHIN the window).
Two independent windowings are used (the job-8 V-A/V-B pair):
  V-A: searchsorted on the restricted (sorted) sample
  V-B: bincount on the restricted (unsorted) sample
Their agreement is the independent check. N-stability is checked at N=1e7
vs N=5e6. The two lobes are compared at MATCHED h and MATCHED window width.

The result is the lobe difference D_2(0.65) - D_2(0.5) at matched h.
Kill (per card): |diff| < 0.02 -> lobes dimensionally indistinguishable.

Empirical/heuristic only. No proof or solution is claimed.
"""
import numpy as np, json, time, math
from pathlib import Path

SLO = "E:/Projects/MillenniumResearch/prize_work/ppl-003"
OUT = Path(SLO) / "lobe_localD2.json"

# Inlined from renyi_limit_N1e8.py (path-independent; avoids importing the
# module, whose top-level calls compute_N and fails outside its CWD).
def renyi_D(counts, h, q):
    m = counts[counts > 0].astype(np.float64)
    m = m / m.sum()
    if q == 0.0:
        return math.log(len(m)) / math.log(1.0 / h)
    if q == 1.0:
        return float(-np.sum(m * np.log(m)) / math.log(1.0 / h))
    return math.log(np.sum(m ** q)) / ((q - 1.0) * math.log(h))
CACHE = {5_000_000: "phi_cache_5000000.npy", 10_000_000: "phi_cache_10000000.npy"}

# matched windows (same width), from the job-261112_5 profile bins
W_LOBE065 = (0.62632, 0.67684)   # bin12, c~0.65 lobe (11.5%, HIGHEST)
W_LOBE050 = (0.47474, 0.52526)   # bin9,  c~0.5  lobe (10.6%)
W_LOBE065_w = (0.60, 0.70)       # wider symmetric (robustness)
W_LOBE050_w = (0.45, 0.55)       # wider symmetric (robustness)
H_GRID = [1e-5, 2e-5, 5e-5, 1e-4, 2e-4]
NS = [10_000_000, 5_000_000]

def load_r(N):
    phi = np.load(Path(SLO) / CACHE[N])
    n = np.arange(1, N + 1, dtype=np.float64)
    return (phi[1:].astype(np.float64) / n)

def local_D2(r, a, b, h):
    rw = r[(r >= a) & (r < b)]
    nw = rw.shape[0]
    nb = int(np.floor((b - a) / h + 1e-9))
    edges = a + np.arange(nb + 1) * h
    # V-A: searchsorted on the restricted SORTED sample (searchsorted needs sorted input)
    rws = np.sort(rw)
    LA = np.searchsorted(rws, edges[:-1], side="left")
    RA = np.searchsorted(rws, edges[1:], side="right")
    cA = np.maximum(RA - LA, 0)
    # V-B: bincount (direct floor) on the UNRESTRICTED-ORDER sample, independent windowing
    idx = np.clip(np.floor((rw - a) / h).astype(np.int64), 0, nb - 1)
    cB = np.bincount(idx, minlength=nb)
    D2_A = renyi_D(cA, h, 2.0)
    D2_B = renyi_D(cB, h, 2.0)
    # point-mass dominance: fraction of window mass in the single largest h-bin
    ctop = np.sort(cB)[::-1]
    top1 = float(ctop[0]) / nw
    top3 = float(ctop[:3].sum()) / nw
    return D2_A, D2_B, int(nw), int(nb), top1, top3

t0 = time.time()
out = {"job_id": "lab_erdos_261112_21", "part": 2, "method":
       "local Renyi D_2 (correlation dimension) of conditional measure mu|_window, "
       "mass normalized within window; V-A searchsorted vs V-B bincount independent "
       "windowing; N-stability 1e7 vs 5e6; matched h and matched window width",
       "windows": {"lobe_065_bin12": list(W_LOBE065), "lobe_050_bin9": list(W_LOBE050),
                   "lobe_065_wide": list(W_LOBE065_w), "lobe_050_wide": list(W_LOBE050_w)},
       "h_grid": H_GRID, "N_list": NS, "by_N": {}}

for N in NS:
    r = load_r(N)
    per_N = {}
    # local mass in each window (fraction of N)
    per_N["local_mass"] = {
        "lobe_065_bin12": float(((r >= W_LOBE065[0]) & (r < W_LOBE065[1])).mean()),
        "lobe_050_bin9": float(((r >= W_LOBE050[0]) & (r < W_LOBE050[1])).mean()),
    }
    per_N["D2_by_h"] = {}
    for h in H_GRID:
        row = {}
        for tag, (a, b) in {"lobe_065": W_LOBE065, "lobe_050": W_LOBE050,
                            "lobe_065_wide": W_LOBE065_w, "lobe_050_wide": W_LOBE050_w}.items():
            dA, dB, nw, nb, top1, top3 = local_D2(r, a, b, h)
            row[tag] = {"D2_A": round(float(dA), 5), "D2_B": round(float(dB), 5),
                        "absVA_VB": round(float(abs(dA - dB)), 6), "n_win": nw, "n_bins": nb,
                        "top1_bin_frac": round(top1, 5), "top3_bin_frac": round(top3, 5)}
        per_N["D2_by_h"][f"{h:.0e}"] = row
    out["by_N"][str(N)] = per_N
    print(f"[{time.time()-t0:5.1f}s] N={N} local_mass 0.65={per_N['local_mass']['lobe_065_bin12']:.4f} "
          f"0.50={per_N['local_mass']['lobe_050_bin9']:.4f}")

# ---- summary: lobe D_2 difference at matched h (primary bin12 vs bin9) ----
N = NS[0]  # primary N=1e7
D2b = out["by_N"][str(N)]["D2_by_h"]
summary = {"primary_N": N, "diff_by_h": {}, "wide_diff_by_h": {}}
max_va = 0.0
for h in H_GRID:
    hk = f"{h:.0e}"
    d65 = D2b[hk]["lobe_065"]["D2_A"]
    d50 = D2b[hk]["lobe_050"]["D2_A"]
    max_va = max(max_va, D2b[hk]["lobe_065"]["absVA_VB"], D2b[hk]["lobe_050"]["absVA_VB"])
    summary["diff_by_h"][hk] = round(d65 - d50, 5)
    summary["wide_diff_by_h"][hk] = round(D2b[hk]["lobe_065_wide"]["D2_A"] -
                                          D2b[hk]["lobe_050_wide"]["D2_A"], 5)
# primary h = finest = 1e-5
summary["diff_primary_h1e5"] = summary["diff_by_h"]["1e-05"]
summary["abs_diff_primary"] = abs(summary["diff_primary_h1e5"])
summary["kill_threshold"] = 0.02
summary["kill_fired"] = bool(abs(summary["diff_primary_h1e5"]) < 0.02)
summary["vA_vB_max_abs"] = round(float(max_va), 6)
# point-mass dominance at primary h (why local D_2 is a finite-N artifact)
summary["point_mass_top1_bin_frac_h1e5"] = {
    "lobe_065": D2b["1e-05"]["lobe_065"]["top1_bin_frac"],
    "lobe_050": D2b["1e-05"]["lobe_050"]["top1_bin_frac"],
}
# N-stability of the primary-h difference
d65_5m = out["by_N"][str(NS[1])]["D2_by_h"]["1e-05"]["lobe_065"]["D2_A"]
d50_5m = out["by_N"][str(NS[1])]["D2_by_h"]["1e-05"]["lobe_050"]["D2_A"]
summary["diff_5e6_h1e5"] = round(d65_5m - d50_5m, 5)
summary["N_stability_delta"] = round(summary["diff_primary_h1e5"] - summary["diff_5e6_h1e5"], 5)
out["summary"] = summary
out["wall_s"] = round(time.time() - t0, 1)

OUT.write_text(json.dumps(out, indent=2))
print(f"\n== SUMMARY (primary N={N}, matched window bin12 vs bin9) ==")
for hk, dv in summary["diff_by_h"].items():
    print(f"  h={hk}: D2(0.65)-D2(0.5) = {dv:+.5f}")
print(f"  wide windows: " + ", ".join(f"{k}:{v:+.5f}" for k, v in summary["wide_diff_by_h"].items()))
print(f"  primary h=1e-5 diff = {summary['diff_primary_h1e5']:+.5f}  |diff| = {summary['abs_diff_primary']:.5f}")
print(f"  V-A vs V-B max abs = {summary['vA_vB_max_abs']:.6f}  (independent windowing agreement)")
print(f"  N=5e6 h=1e-5 diff = {summary['diff_5e6_h1e5']:+.5f}  (N-stability delta {summary['N_stability_delta']:+.5f})")
print(f"  KILL (|diff|<0.02) fired = {summary['kill_fired']}")
print(f"\nwrote {OUT}   (wall {time.time()-t0:.1f}s)")
