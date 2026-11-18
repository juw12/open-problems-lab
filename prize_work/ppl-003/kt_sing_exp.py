"""Kill test (lab_erdos_261111_38): local windowed-slope divergence EXPOENT at c=1/2.

Goal: for c=1/2 (cached phi, N=1e7) measure the one-sided-below windowed slope
  s_-(w) = (F(c) - F(c-w)) / w
over a fine window range above the finite-N floor (N*eps=10 -> floor ~1e-6), fit the
log-log exponent gamma = d log s_-(w) / d log w (job 24 measured ~ -0.9 for the
two-sided slope at c=1/2), and cross-check the one-sided-below MASS M_-(w) against the
sieve-independent exact 2^a*p tower (job 22). If it reproduces job 24's exponent and the
tower matches, the local-exponent machinery is sound -> proceed to all strong points.

Gates:
  (a) M_-(0.5; 5e-5) reproduces job 8's 2p-spike mass 0.0723393 (|d|<1e-5).
  (b) sieve M_-(w) == exact-tower M(w) for every fine window (|d|<1e-4) -> independent.
  (c) log-log fit over w in [5e-5, 0.005] has R^2 > 0.99 and gamma in [-1.2, -0.6].
Empirical only; no proof/Lean/posting. Single-threaded, cached phi (instant).
"""
import numpy as np
from scipy.optimize import curve_fit  # may not be present; fall back to polyfit

N = 10_000_000
phi = np.load("phi_cache_%d.npy" % N)
n = np.arange(1, len(phi), dtype=np.int64)
r = phi[1:].astype(np.float64) / n.astype(np.float64)
rs = np.sort(r)

def F(c):
    return float(np.searchsorted(rs, c, 'left')) / N

def M_below(c, w):
    return F(c) - F(c - w)

def slope_below(c, w):
    return M_below(c, w) / w

# exact 2^a*p tower (sieve-independent): density of n=2^a*p (a>=1, p odd prime)
# with (p-1)/(2p) in (0.5-eps, 0.5)  <=>  1/(2p) < eps  <=>  p > 1/(2eps), p <= N/2^a.
import math
psieve = np.ones(N // 2 + 1, dtype=bool); psieve[:2] = False
for i in range(2, int(math.isqrt(N // 2)) + 1):
    if psieve[i]:
        psieve[i * i::i] = False
P = np.flatnonzero(psieve)
oddP = P[P % 2 == 1]

def tower_mass(eps):
    lo = 1.0 / (2.0 * eps)
    tot = 0
    a = 1
    while (1 << a) <= N:
        cap = N >> a
        if cap < 3:
            break
        # odd primes p with lo < p <= cap
        tot += int(np.searchsorted(oddP, cap, 'right') - np.searchsorted(oddP, lo, 'left'))
        a += 1
    return tot / N

WINDOWS = np.array([5e-3, 2.5e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2.5e-5, 1e-5, 5e-6, 2e-6, 1e-6])

print("=== c=1/2 one-sided-below slope/mass vs window (N=1e7) ===")
print("w        M_sieve   M_tower   dM        slope     slope_tower")
rows = []
for w in WINDOWS:
    ms = M_below(0.5, float(w))
    mt = tower_mass(float(w))
    ss = slope_below(0.5, float(w))
    st = mt / float(w)
    rows.append((float(w), ms, mt, ss))
    print("%.1e  %.6f  %.6f  %+.2e   %.4f   %.4f" % (w, ms, mt, ms - mt, ss, st))

# gate (a): reproduce job-8 spike mass at eps=5e-5
m5e5 = M_below(0.5, 5e-5)
print("\nGATE(a) M_-(0.5;5e-5)=%.6f  job8=0.0723393  match=%s" %
      (m5e5, abs(m5e5 - 0.0723393) < 1e-5))

# gate (b): sieve vs tower over the fit range
dm = [abs(M_below(0.5, float(w)) - tower_mass(float(w))) for w in WINDOWS]
print("GATE(b) max |sieveM - towerM| over all windows = %.2e  (all<1e-4: %s)" %
      (max(dm), all(d < 1e-4 for d in dm)))

# gate (c): log-log fit of the SLOPE over w in [5e-5, 0.005]
fitW = WINDOWS[(WINDOWS >= 5e-5) & (WINDOWS <= 5e-3)]
xs = np.log(fitW.astype(float))
ys = np.log([M_below(0.5, float(w)) / float(w) for w in fitW])
A = np.polyfit(xs, ys, 1)
gamma = A[0]
pred = np.polyval(A, xs)
ss_res = float(np.sum((ys - pred) ** 2)); ss_tot = float(np.sum((ys - ys.mean()) ** 2))
r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
print("\nGATE(c) slope log-log fit over w in [5e-5,5e-3]:")
print("  gamma = dlog(s_+)/dlog(w) = %.4f   R^2 = %.4f" % (gamma, r2))
print("  in-range? gamma in [-1.2,-0.6]: %s   R^2>0.99: %s" %
      (-1.2 <= gamma <= -0.6, r2 > 0.99))
print("\nKILL-TEST VERDICT: pass = (a) and (b) and (c) ->",
      (abs(m5e5 - 0.0723393) < 1e-5) and (all(d < 1e-4 for d in dm)) and (-1.2 <= gamma <= -0.6) and (r2 > 0.99))
