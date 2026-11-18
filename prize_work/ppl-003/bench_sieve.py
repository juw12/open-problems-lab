"""Benchmark the vectorized totient sieve at N=1e8 to confirm it fits the 20-min budget.
Times: (1) bool prime sieve, (2) vectorized int32 phi sieve, (3) float multiplicative r sieve.
Empirical only. Single-threaded.
"""
import time
import numpy as np

def prime_bool(N):
    sp = np.ones(N + 1, dtype=bool); sp[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if sp[i]: sp[i*i::i] = False
    return sp

def phi_sieve_int32(N, primes):
    phi = np.arange(N + 1, dtype=np.int32)
    for p in primes:
        phi[p::p] = (phi[p::p] // p) * (p - 1)
    return phi

def r_sieve_float(N, primes):
    r = np.ones(N + 1, dtype=np.float64)
    for p in primes:
        r[p::p] *= (1.0 - 1.0 / p)
    return r

N = 100_000_000
t0 = time.time()
sp = prime_bool(N)
primes = np.nonzero(sp)[0]
t1 = time.time()
print("prime bool sieve: %.2fs  (nprimes=%d)" % (t1 - t0, len(primes)))

t0 = time.time()
phi = phi_sieve_int32(N, primes)
t1 = time.time()
mean_r = phi[1:].astype(np.float64).mean() / 1e8  # rough mean phi / N
print("int32 phi sieve: %.2fs  (phi[N]=%d, spot phi[100]=%d)" % (t1 - t0, int(phi[N]), int(phi[100])))

# free
del phi

t0 = time.time()
r = r_sieve_float(N, primes)
t1 = time.time()
print("float r sieve: %.2fs  (r[N]=%.9f)" % (t1 - t0, float(r[N])))
del r, sp, primes
print("TOTAL elapsed ~ (see above)")
