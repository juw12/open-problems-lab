"""lab_erdos_261112_38 - flattened single-bincount W-builder for the L>200 sweep.

build_W_np (bounds200.py, job-35 referee-passed) does 24 separate bincounts
(one per Omega-level r=1..RMAX). build_W_flat here does ONE bincount over the
flattened index (j*(RMAX+1)+r), which is ~24x cheaper per L.

The index max is (L-1)*(RMAX+1)+RMAX < 50000 for L<=2000, so it fits int32; we
keep the index int32 (no int64 allocation). Bit-exactness: both methods
accumulate each (j,r) bucket as  sum_{n: n%L==j, omega(n)==r} wt[n]  in
increasing-n order (numpy bincount adds weights in input order); the flattened
form feeds the same per-bucket weight sequence, so the float64 sums are
identical. Confirmed: the sweep gates this (K1) and the verifier re-derives it
from an independent Method-3 pure-Python build (O2a/O2b). Timetest showed
max|W_i32 - W_i64| = 0.0 at N=5e7.
"""
import numpy as np


def masked_inputs(n_np, om_np, wt_np, rmax):
    """Precompute the (1<=omega<=rmax) masked arrays once for reuse across L's.
    n_m is int32, om_i32 is int32, wt_m is float64."""
    m = (om_np >= 1) & (om_np <= rmax)
    n_m = n_np[m]
    om_i32 = om_np[m].astype(np.int32)
    wt_m = wt_np[m]
    return n_m, om_i32, wt_m


def build_W_flat_precomputed(n_m, om_i32, wt_m, L, rmax):
    """One flattened int32 bincount. n_m/om_i32/wt_m precomputed once."""
    idx = n_m % L
    idx *= (rmax + 1)
    idx += om_i32
    return np.bincount(idx, weights=wt_m, minlength=L * (rmax + 1)).reshape(L, rmax + 1)


def build_W_flat(n_np, om_np, wt_np, L, rmax):
    """Convenience: recompute the mask each call (fine for a handful of L's)."""
    n_m, om_i32, wt_m = masked_inputs(n_np, om_np, wt_np, rmax)
    return build_W_flat_precomputed(n_m, om_i32, wt_m, L, rmax)
