#!/usr/bin/env python3
# job lab_erdos_261112_13 -- quick robustness: is the c=2/3 > c=1/2 sharpness ranking window-stable?
import os, math
import numpy as np
D = r'E:\Projects\MillenniumResearch\prize_work\ppl-003'
phi = np.load(os.path.join(D, 'phi_cache_10000000.npy'))
n = np.arange(10000001, dtype=np.float64)
r = phi[1:].astype(np.float64) / n[1:]; rs = np.sort(r); N = 10000000
def growth(c, wL, wS):
    s0 = (np.searchsorted(rs, c + wL, 'left') - np.searchsorted(rs, c - wL, 'left')) / (2*wL*N)
    s1 = (np.searchsorted(rs, c + wS, 'left') - np.searchsorted(rs, c - wS, 'left')) / (2*wS*N)
    return s1/s0 if s0>0 else float('nan')
def df(c, w=0.02):
    return (np.searchsorted(rs, c + w, 'left') - np.searchsorted(rs, c - w, 'left'))/N
c23, c12 = 2/3, 0.5
print(f"{'windows':16s} {'g(2/3)':>9s} {'g(1/2)':>9s} {'2/3>1/2':>8s}  {'df(2/3)':>8s} {'df(1/2)':>8s}")
for wL,wS in [(0.05,0.00025),(0.03,0.0003),(0.02,0.0002),(0.1,0.001)]:
    g23, g12 = growth(c23,wL,wS), growth(c12,wL,wS)
    print(f"{f'{wL}/{wS}':16s} {g23:9.3f} {g12:9.3f} {str(g23>g12):>8s}  {df(c23):8.4f} {df(c12):8.4f}")
