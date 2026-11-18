"""Diagnose the 48-n residual at N=1e7, eps=5e-5: what are the n in [0.5-eps,0.5) that are NOT of the
form 2^a * p (a>=1, p odd prime)?  Classify by radical / odd-part factorization and by float-boundary
proximity, to decide whether the residual is a ulp boundary artifact or a genuine family."""
import math
import numpy as np
from fractions import Fraction

N = 10_000_000
eps = 5e-5
lo = 0.5 - eps

def is_prime_array(N):
    isp = np.ones(N + 1, dtype=bool); isp[:2] = False
    for i in range(2, int(math.isqrt(N)) + 1):
        if isp[i]: isp[i*i::i] = False
    return isp

isp = is_prime_array(N)
phi = np.arange(N + 1, dtype=np.int64)
for p in np.nonzero(isp)[0]:
    for j in range(p, N + 1, p):
        phi[j] = phi[j] // p * (p - 1)
n = np.arange(1, N + 1)
r = phi[1:].astype(np.float64) / n.astype(np.float64)
in_win = (r >= lo) & (r < 0.5)
lowbit = n & -n
v2 = np.log2(lowbit.astype(np.float64)).astype(np.int64)
oddpart = n // lowbit
is_tower = (v2 >= 1) & (oddpart >= 3) & isp[oddpart]

resid_n = n[in_win & ~is_tower]
print("residual count:", len(resid_n))
if len(resid_n):
    # classify odd parts: prime power? semiprime? radical size?
    from collections import Counter
    op = oddpart[in_win & ~is_tower]
    # exact rational value for each residual n and its distance to lo
    print("first 40 residual n (n, v2, oddpart, is_oddpart_pp, exact r, dist_to_lo):")
    def factor_small(m):
        fs = {}
        d = 2
        while d*d <= m and d <= 100000:
            while m % d == 0:
                fs[d] = fs.get(d, 0) + 1
                m //= d
            d += 1 if d == 2 else 2
        if m > 1: fs[m] = fs.get(m, 0) + 1
        return fs
    cnt = Counter()
    for nn in resid_n[:200]:
        f = factor_small(int(oddpart[nn-1]))
        key = tuple(sorted(f.items()))
        cnt[key] += 1
    print("odd-part factor patterns (top):")
    for k, v in cnt.most_common(10):
        print("   ", k, "x", v)
    # exact distance to lo for the smallest |r-lo|
    dd = []
    for nn in resid_n:
        exact = Fraction(int(phi[nn]), int(nn))
        d = exact - Fraction(1,2) + Fraction(1,20000)  # r - lo  (lo = 1/2 - 1/20000)
        dd.append((abs(float(d)), int(nn), int(oddpart[nn-1]), int(v2[nn-1])))
    dd.sort()
    print("closest-to-boundary residual (|r-lo|, n, oddpart, v2):")
    for x in dd[:10]:
        print("   ", x)
