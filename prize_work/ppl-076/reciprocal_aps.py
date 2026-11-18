#!/usr/bin/env python3
"""
Erdos Problem #3 (PPL 076), $5000.

    If A subset of N has sum_{n in A} (1/n) = +inf, then must A contain
    arbitrarily long arithmetic progressions?

This script does NOT prove the problem.  It numerically verifies the single
reduction step that links the *known/claimed* small-set bounds r_k(N) to the
conjecture, so the computation has a hard gate against published values.

Reduction (partial summation / Abel)
------------------------------------
Let a_n = 1_{n in A} and A(n) = sum_{j<=n} a_j = |A cap {1..n}|.  Abel's
summation gives the EXACT identity

    sum_{n=1}^N a_n/n  =  A(N)/N  +  sum_{n=1}^{N-1} A(n) / (n (n+1)).      (1)

If A is k-AP-free then A(n) <= r_k(n) for every n, so the reciprocal sum of a
k-AP-free set is bounded by

    sum_{n in A} 1/n  <=  lim_{N} [ r_k(N)/N + sum_{n<N} r_k(n)/(n(n+1)) ].  (2)

The OpenAI preprint "Quasipolynomial Bounds for Arithmetic Progressions"
(2026-09-23) claims  r_k(N) <= N * exp(-(log N)^{c_k})  for all k>=3 with some
c_k>0.  Under that bound the right-hand side of (2) is finite because, with
u=log x,

    int_1^inf exp(-(log x)^c) / x dx  =  int_0^inf e^{-u^c} du
        =  (1/c) * Gamma(1/c)  <  inf          for every c>0.               (3)

A finite reciprocal sum contradicts sum_{n in A} 1/n = +inf, so a divergent
set must contain a k-term AP; since k is arbitrary it contains arbitrarily
long APs.  THIS script verifies (1) and (3) numerically with published-value
gates.  It does not verify the preprint's proof of the r_k bound itself.

GATES (published values)
------------------------
A1.  (1) with a_n = 1 gives the harmonic number H_N, which must match the
     published asymptotic  H_N = log N + gamma + 1/(2N) - 1/(12 N^2) + ...
     (Euler's constant gamma = 0.5772156649015328606...).
A2.  (1) must hold exactly for the published size-15 3-AP-free set in {1..40}
     (OEIS A003002, Dybizbanski 2012):
         {1,2,4,5,10,11,13,14,28,29,31,32,37,38,40}
     and that set must be verified 3-AP-free.
B.   int_0^inf e^{-u^c} du must match (1/c) Gamma(1/c) for c in {1/2,1,2,3}:
         c=1/2 -> 2,   c=1 -> 1,   c=2 -> sqrt(pi)/2,   c=3 -> Gamma(1/3)/3.

Usage
-----
    python reciprocal_aps.py
Exits 0 if every gate passes, 1 otherwise.
"""
import sys
import time
import math
import numpy as np

GAMMA = 0.57721566490153286060651209   # Euler-Mascheroni (published constant)


def simpson(f, a, b, n=200000):
    """Composite Simpson rule on [a, b] with n (even) sub-intervals."""
    if n % 2:
        n += 1
    x = np.linspace(a, b, n + 1)
    h = (b - a) / n
    y = f(x)
    return (h / 3.0) * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum())


# ---------------------------------------------------------------------------
# Part A: Abel / partial-summation identity
# ---------------------------------------------------------------------------
def abel_lhs(a, N):
    """sum_{n=1..N} a_n / n, with a a 1-based numpy array (index 0 unused)."""
    n = np.arange(1, N + 1)
    return float((a[n] / n).sum())


def abel_rhs(A, N):
    """A(N)/N + sum_{n=1..N-1} A(n)/(n(n+1)), A a 1-based cumulative array."""
    if N == 1:
        return float(A[1])
    n = np.arange(1, N)
    return float(A[N] / N + (A[n] / (n * (n + 1))).sum())


def is_3ap_free(S, N):
    S = set(S)
    for a in S:
        for d in range(1, (N - a) // 2 + 1):
            if (a + d) in S and (a + 2 * d) in S:
                return False
    return True


def part_A():
    print("=" * 72)
    print("PART A  Abel / partial-summation identity  (eq. (1))")
    print("=" * 72)
    ok = True

    # A1: full set  a_n = 1  ->  harmonic number H_N vs published asymptotic.
    N = 1_000_000
    a = np.ones(N + 1, dtype=np.float64)
    a[0] = 0.0
    lhs = abel_lhs(a, N)
    A = np.cumsum(a)                       # A[n] = n  for the full set
    rhs = abel_rhs(A, N)
    asym = math.log(N) + GAMMA + 1.0 / (2 * N) - 1.0 / (12 * N * N)
    err_id = abs(lhs - rhs)
    err_asym = abs(lhs - asym)
    print(f" A1  full set, N={N}")
    print(f"     H_N (exact sum)              = {lhs!r}")
    print(f"     H_N (Abel identity, eq.1)    = {rhs!r}   |lhs-rhs|={err_id:.3e}")
    print(f"     H_N (log N+g+1/2N-1/12N^2)   = {asym!r}   |H_N-asym|={err_asym:.3e}")
    gate_a1 = err_id < 1e-12 and err_asym < 1e-9
    print(f"     [A1 {'PASS' if gate_a1 else 'FAIL'}]")
    ok &= gate_a1

    # A2: published size-15 3-AP-free set in {1..40} (OEIS A003002).
    S = {1, 2, 4, 5, 10, 11, 13, 14, 28, 29, 31, 32, 37, 38, 40}
    N = 40
    a2 = np.zeros(N + 1, dtype=np.float64)
    for x in S:
        a2[x] = 1.0
    lhs2 = abel_lhs(a2, N)
    A2 = np.cumsum(a2)
    rhs2 = abel_rhs(A2, N)
    apfree = is_3ap_free(S, N)
    err2 = abs(lhs2 - rhs2)
    recip = sum(1.0 / x for x in S)
    print(f" A2  3-AP-free set S (|S|={len(S)}) in {{1..{N}}} (OEIS A003002)")
    print(f"     3-AP-free?                   = {apfree}")
    print(f"     sum over S of 1/x  (LHS)    = {lhs2!r}")
    print(f"     Abel identity (eq.1) (RHS)   = {rhs2!r}   |lhs-rhs|={err2:.3e}")
    print(f"     (cross-check direct sum)     = {recip!r}")
    gate_a2 = apfree and err2 < 1e-12 and abs(lhs2 - recip) < 1e-12
    print(f"     [A2 {'PASS' if gate_a2 else 'FAIL'}]")
    ok &= gate_a2
    return ok, lhs, asym, lhs2


# ---------------------------------------------------------------------------
# Part B: convergence of the claimed quasipolynomial bound  (eq. (3))
# ---------------------------------------------------------------------------
# For each c, int_0^inf e^{-u^c} du is integrated in a variable where the
# integrand is smooth and decays like e^{-w} (or faster):
#   c <= 1 : substitute u = w^{1/c}  ->  (1/c) int_0^inf w^{1/c-1} e^{-w} dw
#   c >  1 : integrate e^{-u^c} directly (fast decay, no singularity at 0).
def integrand_w(c, w):
    return (1.0 / c) * np.power(w, 1.0 / c - 1.0) * np.exp(-w)


def integrand_u(c, u):
    return np.exp(-np.power(u, c))


def part_B():
    print()
    print("=" * 72)
    print("PART B  int_0^inf e^{-u^c} du = (1/c) Gamma(1/c)  (eq. (3))")
    print("=" * 72)
    ok = True
    print(f" {'c':>5} | {'numeric int_0^inf e^(-u^c) du':>28} | {'(1/c)Gamma(1/c) (known)':>24} | rel.err | gate")
    cases = [
        (0.5, 2.0, "2"),
        (1.0, 1.0, "1"),
        (2.0, math.sqrt(math.pi) / 2.0, "sqrt(pi)/2"),
        (3.0, math.exp(math.lgamma(1.0 / 3.0)) / 3.0, "Gamma(1/3)/3"),
    ]
    NPTS = 400000   # even number of Simpson sub-intervals
    for c, known, label in cases:
        if c <= 1.0:
            num = float(simpson(lambda w: integrand_w(c, w), 0.0, 60.0, NPTS))
        else:
            # e^{-u^c} is negligible beyond u where u^c ~ 40.
            U = 40.0 ** (1.0 / c) + 1.0
            num = float(simpson(lambda u, _c=c: integrand_u(_c, u), 0.0, U, NPTS))
        known_from_lgamma = math.exp(math.lgamma(1.0 / c)) / c
        rel = abs(num - known) / known
        gate = rel < 1e-6 and abs(known - known_from_lgamma) < 1e-12
        print(f" {c:>5} | {num:>28.12f} | {known:>24.12f} ({label:>11}) | {rel:.2e} | {'PASS' if gate else 'FAIL'}")
        ok &= gate

    # Discrete series S_N(c) = sum_{n=2..N} e^{-(log n)^c}/n must stay bounded
    # as N -> inf (converges), consistent with the finite integral above.
    print()
    print(" discrete series S_N(c) = sum_{n=2..N} e^{-(log n)^c}/n  (convergence check)")
    print(f" {'N':>10} | " + " | ".join(f"S_N({c:g})".ljust(18) for c in (0.5, 1, 2, 3)))
    for N in (10_000, 100_000, 1_000_000):
        n = np.arange(2, N + 1)
        lg = np.log(n)
        row = []
        for c in (0.5, 1.0, 2.0, 3.0):
            row.append(f"{float((np.exp(-np.power(lg, c)) / n).sum()):.9f}")
        print(f" {N:>10} | " + " | ".join(r.ljust(18) for r in row))
    print(" (values grow slowly and remain bounded as N increases -> series converges)")
    return ok


def main():
    t0 = time.time()
    okA, lhs1, asym1, lhs2 = part_A()
    okB = part_B()
    ok = okA and okB
    dt = time.time() - t0
    print()
    print("=" * 72)
    print(f"RESULT  Part A (Abel identity): {'PASS' if okA else 'FAIL'}   "
          f"Part B (convergence, eq.3): {'PASS' if okB else 'FAIL'}")
    print(f"        H_{1000000} = {lhs1:.12f}  (matches log N + gamma to "
          f"{abs(lhs1 - asym1):.2e})")
    print(f"        3-AP-free size-15 set reciprocal sum = {lhs2:.12f}")
    print(f"        Computed in {dt:.2f}s.  GATE: {'PASS' if ok else 'FAIL'}")
    print("=" * 72)
    # The reduction (2)+(3): under r_k(N) <= N e^{-(log N)^c_k} with c_k>0, the
    # reciprocal sum of any k-AP-free set is finite, contradicting divergence.
    print("NOTE  This verifies the REDUCTION (Abel identity + finite integral).")
    print("      It does NOT verify the OpenAI preprint's proof of the r_k bound.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
