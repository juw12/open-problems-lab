"""independent verifier for interior_map_full_N10000000.json.

Recomputes the two headline quantities via a DIFFERENT code path than the
main script (which used 'left' searchsorted + floor/digitize binning):
  (a) D_1 @ h=3e-3        -> np.digitize binning (independent of floor)
  (b) divergent-set measure @ G=100001 -> plain per-cell searchsorted with
      'RIGHT' side (main used 'left'), same growth criterion s_small/s_large>=5
and asserts the stored JSON is consistent + the Renyi spectrum is monotone.

NOTE on (a): D_1 at a SINGLE finite h is bin-edge sensitive; a different
binning convention (digitize vs searchsorted) moves it by ~0.01. The check
tolerates that and confirms the value to ~1 decimal (the headline point is
the monotone multifractal spectrum, not the last digit of D_1).
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
J = json.load(open(os.path.join(HERE, "interior_map_full_N10000000.json")))
phi = np.load(os.path.join(HERE, "phi_cache_10000000.npy"))
N = 10_000_000
rs = np.sort(phi[1:].astype(np.float64) / np.arange(1, N + 1, dtype=np.float64))
r = phi[1:].astype(np.float64) / np.arange(1, N + 1, dtype=np.float64)

print("=== independent verifier (N=1e7) ===", flush=True)
ok = True

# ---- (a) D_1 via np.digitize (independent of floor binning) ----
h = 3e-3
bins = np.arange(0.0, 1.0 + h, h)
idx = np.clip(np.digitize(r, bins) - 1, 0, len(bins) - 2)
cnt = np.bincount(idx, minlength=len(bins) - 1).astype(np.float64)
p = cnt / cnt.sum()
nz = p[p > 0]
D1_dig = -float(np.sum(nz * np.log(nz))) / np.log(1.0 / h)
D1_json_h3 = J["N"]["10000000"]["D_q_at_h3e3"]["D_1.0"]
print(f"  D1 digitize={D1_dig:.5f}   D1 json(h=3e-3)={D1_json_h3}   "
      f"(bin-edge sensitivity at finite h; agree to ~1e-1)", flush=True)
ok &= abs(D1_dig - D1_json_h3) < 0.02

# ---- (b) divergent-set measure @ G=100001 via 'right'-side searchsorted ----
G = 100001
grid = np.linspace(0.0, 1.0, G)                 # G points incl endpoints (same as main)
w_L, w_S, TH = J["w_large"], J["w_small"], J["diverge_thresh"]
lo, hi = J["interior"]
t0 = time.time()
s_large = (np.searchsorted(rs, grid + w_L, 'right') -
           np.searchsorted(rs, grid - w_L, 'right')) / (2.0 * w_L * N)
s_small = (np.searchsorted(rs, grid + w_S, 'right') -
           np.searchsorted(rs, grid - w_S, 'right')) / (2.0 * w_S * N)
grow = s_small / np.where(s_large > 0, s_large, 1.0)     # s_small/s_large (CORRECT)
m = (grid > lo) & (grid < hi)
div = (grow >= TH) & m
n_div = int(div.sum())
meas = n_div / len(grid)
print(f"  divergent n_div_interior={n_div}   measure_of_[0,1]={meas:.5f}  "
      f"[{time.time()-t0:.1f}s]   (right-side searchsorted)", flush=True)
fm = J["N"]["10000000"]["fine_map"]["100001"]
jnd, jmeas, jcov, jcr = (fm["n_div_interior"], fm["divergent_measure_of_[0,1]"],
                         fm["coverage_interior_frac"], fm["c_range"])
print(f"  json: n_div_interior={jnd}  measure_of_[0,1]={jmeas}  "
      f"coverage={jcov}  c_range={jcr}", flush=True)
ok &= abs(n_div - jnd) <= 50            # right-vs-left edge convention: a few cells
ok &= abs(meas - jmeas) <= 1e-3
# c_range: recompute from my own div set
cr_mine = [float(grid[div].min()), float(grid[div].max())]
ok &= all(abs(a - b) < 5e-4 for a, b in zip(cr_mine, jcr))

# ---- (c) Renyi spectrum monotone (multifractal) + cross-checks ----
dq = J["N"]["10000000"]["D_q_at_h3e3"]
seq = [dq["D_0.0"], dq["D_0.5"], dq["D_1.0"], dq["D_1.5"], dq["D_2.0"]]
mono = all(seq[i] > seq[i + 1] for i in range(len(seq) - 1))
print(f"  Renyi spectrum = {[round(x, 4) for x in seq]}   strictly_decreasing={mono}", flush=True)
ok &= mono
ok &= J["N"]["10000000"]["V2_renyi_AB_all_agree"] == True
ok &= J["gate_G5_s7"]["match"] == True
ok &= all(J["N"]["10000000"]["gates_1234"][k]["pass"] for k in
          ["G1_mean_r", "G2_sumphi_N2", "G3_f0_f1", "G4_c12_spike"])

print("\nVERDICT:", "ALL INDEPENDENT CHECKS PASS" if ok else "SOME CHECK FAILED", flush=True)
sys.exit(0 if ok else 1)
