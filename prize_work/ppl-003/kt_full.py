#!/usr/bin/env python3
# KILL TEST (job lab_erdos_261112_5 RETRY) for the fine-resolution interior map.
# Question: does a 1e-5-resolution grid at N=1e7 reveal substantially more divergent
# structure than the existing coarse 4001-grid backstop (72 cells)?  If the fine grid
# gives only ~72 cells, the "fine map" is a trivial re-sample -> card killed -> do Renyi only.
# Reuses the EXACT referee-passed machinery from scan_s7 (two-sided windowed slope).
import os, sys, time, math
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
from scan_s7 import (W_LARGE, W_SMALL, WINDOWS, DIVERGE_THRESH, PRIMES_MAX,
                     primes_upto)

def grid_growth(r_sorted, N, G):
    grid = np.linspace(0.0, 1.0, G)
    w0, w1 = WINDOWS[0], WINDOWS[-1]
    s0 = (np.searchsorted(r_sorted, grid + w0, 'left') -
          np.searchsorted(r_sorted, grid - w0, 'left')) / (2.0 * w0 * N)
    s1 = (np.searchsorted(r_sorted, grid + w1, 'left') -
          np.searchsorted(r_sorted, grid - w1, 'left')) / (2.0 * w1 * N)
    g = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), 0.0)
    return grid, g

def main():
    N = 10_000_000
    phi = np.load(os.path.join(BASE, f"phi_cache_{N}.npy"))
    nn = np.arange(1, N + 1, dtype=np.float64)
    r_sorted = np.sort(phi[1:].astype(np.float64) / nn)
    # sanity gate: mean(r) = 6/pi^2
    mean_r = float(r_sorted.mean())
    g1 = abs(mean_r - 6 / math.pi ** 2) < 5e-4
    print(f"gate mean_r=6/pi^2: {g1} (measured {mean_r:.6f} vs {6/math.pi**2:.6f})", flush=True)

    for G in [4001, 100001]:
        t0 = time.time()
        grid, g = grid_growth(r_sorted, N, G)
        # interior only: exclude edge bands [0,0.02] and [0.98,1]
        m = (grid > 0.02) & (grid < 0.98)
        div = np.flatnonzero((g >= DIVERGE_THRESH) & m)
        gd = grid[div]
        # coverage = fraction of interior covered by divergent cells
        cov = float(np.mean((g >= DIVERGE_THRESH) & m))
        # 20-bin spatial profile (interior [0.02,0.98])
        print(f"G={G}: n_div(interior)={len(div)}  coverage_interior={cov:.4f}  "
              f"[{time.time()-t0:.1f}s]  c-range of div: "
              f"[{gd.min():.4f},{gd.max():.4f}]" if len(div) else
              f"G={G}: n_div(interior)=0", flush=True)
        if G == 100001 and len(div):
            edges = np.linspace(0.02, 0.98, 21)
            hist, _ = np.histogram(gd, bins=edges)
            tot = hist.sum()
            prof = [f"{h/tot*100:.0f}" for h in hist]
            print("  19-bin divergence density (% of divergent cells, interior):",
                  " ".join(prof), flush=True)

if __name__ == "__main__":
    main()
