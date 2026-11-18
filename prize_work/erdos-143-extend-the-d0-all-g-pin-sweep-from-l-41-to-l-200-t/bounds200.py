#!/usr/bin/env python3
"""
Shared quantities for job lab_erdos_261112_35 (Erdos #143 / PPL-014).

Extends the job-261111_23 machinery (prize_work/ppl-014/d0all_harden.py):
  W[j][r]  = sum_{2<=n<=N, n mod L == j, Omega(n) == r} wt[n],  j=0..L-1, r=0..rmax
  wt[n]    = 1/(n log n), zeroed at n in {2,3,4} (excl_set(10))
  (exhibit residues S are subsets of {1..L-1}, same convention as jobs 21/31/37)

Quantities per (L,k), k with 4k < L (density k/L < 1/4):
  U'(L,k) = max_r sum_{j in S*} W[j][r],  S* = k classes with smallest
           per-class best weight w_j = max_r W[j][r]
           -> VALID UPPER bound on g*(L,k)  (min-max <= value at k-argmin set)
  U(L,k)  = sum of the k smallest w_j      (looser upper bound, cross-check)
  L'(L,k) = max_r sum of the k smallest W[j][r]  (j=1..L-1, at level r)
           -> VALID LOWER bound on g*(L,k)  (min_S max_r >= max_r min_S)
           L' > 1/6 would PROVE a genuine crossing below density 1/4.

g*(L,k) = min_{|S|=k} max_r sum_{j in S} W[j][r]  (worst-exhibit best level).

The numpy build_W (build_W_np) is a vectorized re-implementation of the
pure-Python build_W (density_threshold.py) and must pass the gate in
sweep200_5e7.py before any 5e7 sweep.
"""
import numpy as np

THRESH = 1.0 / 6.0
RMAX = 24
EXCL = (2, 3, 4)   # excl_set(10)


def kmax_below_quarter(L):
    """Largest k with k/L < 1/4, i.e. 4k < L."""
    return (L - 1) // 4


def per_class_best(W, L, rmax=RMAX):
    """w_j = max_{1<=r<=rmax} W[j][r], for j=1..L-1."""
    return [max(W[j][r] for r in range(1, rmax + 1)) for j in range(1, L)]


def U_maxlevel(W, L, k, rmax=RMAX):
    """U'(L,k) = max_r sum_{j in S*} W[j][r], S* = k smallest w_j classes.
    (identical to job-261111_23 U_maxlevel)"""
    w = per_class_best(W, L, rmax)
    idx = list(range(L - 1))
    idx.sort(key=lambda t: w[t])
    Sstar = [t + 1 for t in idx[:k]]
    best = 0.0
    for r in range(1, rmax + 1):
        s = 0.0
        for j in Sstar:
            s += W[j][r]
        if s > best:
            best = s
    return best


def U_sum(W, L, k, rmax=RMAX):
    """U(L,k) = sum of the k smallest per-class best-level weights."""
    w = per_class_best(W, L, rmax)
    w.sort()
    return sum(w[:k])


def L_prime(W, L, k, rmax=RMAX):
    """L'(L,k) = max_r sum of the k smallest W[j][r] at level r (j=1..L-1).
    Lower bound on g*(L,k): min_S max_r F(S,r) >= max_r min_S F(S,r)."""
    best = 0.0
    for r in range(1, rmax + 1):
        vals = sorted(W[j][r] for j in range(1, L))
        s = sum(vals[:k])
        if s > best:
            best = s
    return best


def build_W_np(om_np, wt_np, n_np, L, rmax=RMAX):
    """Vectorized W table (numpy).  om_np int8 Omega, wt_np float64 weights,
    n_np int32 arange(2..N).  Returns float64 (L, rmax+1) array, W[j][r]."""
    W = np.zeros((L, rmax + 1), dtype=np.float64)
    for r in range(1, rmax + 1):
        m = om_np == r
        if m.any():
            W[:, r] = np.bincount(n_np[m] % L, weights=wt_np[m], minlength=L)
    return W
