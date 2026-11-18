#!/usr/bin/env python3
"""#143 lab_erdos_261112_43 diagnostic: isolate the O3 dUp~1e-2 discrepancy.

Compares, on the SAME W at N=5e6 for L in {2001,2304}:
  W builders : build_W_flat_precomputed (sweep)  vs  build_W_np (24-pass)
  U' bounds  : bounds_all (sweep, vectorised)  vs  U_maxlevel (ORIGINAL reference, bounds200)
               vs  U_prime_indep (verifier pure-python)
Whichever matches U_maxlevel (the job-38-validated reference) is correct.
"""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ppl-014"))
sys.path.insert(0, HERE)
from density_threshold import omega_spf, RMAX
from bounds200 import build_W_np, U_maxlevel, L_prime, kmax_below_quarter
from bounds2000 import masked_inputs, build_W_flat_precomputed
from sweep10k_5e6 import bounds_all, make_wt
from verify10k import U_prime_indep, L_prime_indep, w2_direct, omega_prime_power

_T0 = time.time()
N = 5_000_000
print(f"[diag] N={N}", flush=True)
om = np.asarray(omega_spf(N), dtype=np.int32)
om3 = omega_prime_power(N)
wt = make_wt(N); n = np.arange(N + 1, dtype=np.int32)
n_m, om_i32, wt_m = masked_inputs(n, om, wt, RMAX)
om3_np = np.frombuffer(om3, dtype=np.int8)

for L in (2001, 2304):
    Wf = build_W_flat_precomputed(n_m, om_i32, wt_m, L, RMAX)   # sweep W
    Wn = build_W_np(om, wt, n, L, RMAX)                          # 24-pass W
    W2 = w2_direct(om3, wt, N, L, RMAX)                           # pure-python W (Method-3)
    W2n = np.array(W2, dtype=np.float64)
    d_fn = float(np.abs(Wf - Wn).max())
    d_n2 = float(np.abs(Wn - W2n).max())
    kmax = kmax_below_quarter(L)
    Uf, _ = bounds_all(Wf, L, kmax, RMAX)
    Un, _ = bounds_all(Wn, L, kmax, RMAX)
    U_ref = U_maxlevel(Wn, L, kmax, RMAX)        # ORIGINAL reference (per-k)
    U_flast = float(Uf[-1]); Un_last = float(Un[-1])
    U_indep = U_prime_indep(Wn, L, kmax, RMAX)   # verifier pure-python
    print(f"[diag] L={L} kmax={kmax}", flush=True)
    print(f"   W  |Wf-Wn|={d_fn:.3e}  |Wn-W2|={d_n2:.3e}   (W builders)", flush=True)
    print(f"   U'(L,kmax) bounds_all(Wf)={U_flast:.10f}", flush=True)
    print(f"                bounds_all(Wn)={Un_last:.10f}", flush=True)
    print(f"                U_maxlevel(Wn) [REFERENCE]={U_ref:.10f}", flush=True)
    print(f"                U_prime_indep(Wn)={U_indep:.10f}", flush=True)
    print(f"   d(bounds_all(Wn), ref)={abs(Un_last-U_ref):.3e}   d(U_prime_indep, ref)={abs(U_indep-U_ref):.3e}", flush=True)
    print(f"   t={time.time()-_T0:.1f}s", flush=True)
