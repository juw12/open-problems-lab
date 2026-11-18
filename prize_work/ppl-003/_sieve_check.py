# -*- coding: utf-8 -*-
"""Quick check: cache dtype + my totient sieve matches the cached sieve."""
import numpy as np, math, time

def totient_sieve(N):
    phi = np.arange(N+1, dtype=np.int64)
    for i in range(2, N+1):
        if phi[i] == i:            # i is prime
            phi[i::i] -= phi[i::i] // i
    return phi

# cache dtype
c = np.load("phi_cache_10000000.npy")
print("cache 1e7 dtype:", c.dtype, "shape:", c.shape, "len:", len(c))
print("cache[0:6] =", c[0:6])
print("cache phi[1000] =", c[1000], " expected phi(1000)=400")
print("cache phi[999999] =", c[999999])
print("cache mean phi/n:", float(np.mean(c[1:].astype(float)/np.arange(1,len(c),dtype=float))))

# my sieve at small N to verify correctness + speed
t0=time.time()
phi = totient_sieve(2_000_000)
t1=time.time()
print("\nmy sieve N=2e6 time: %.2fs" % (t1-t0))
print("my phi[1000] =", phi[1000], " expected 400")
print("my phi[2000000] =", phi[2000000])
print("my mean phi/n:", float(np.mean(phi[1:].astype(float)/np.arange(1,2000001,dtype=float))), " expect 6/pi^2=", round(6.0/math.pi**2,6))
# spot-check a few known totients
for n in [1,2,6,12,60,100,1000,999999,2000000]:
    print("  phi(%d) = %d" % (n, phi[n]))
