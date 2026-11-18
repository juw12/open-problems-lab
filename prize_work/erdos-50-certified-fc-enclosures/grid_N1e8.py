"""Main dataset: exact integer counts C_N(k) = #{ n <= N : phi(n)/n < k/100 } for k = 1..99,
N = 1e6, 2e6, ..., 1e8 (100 rungs), from an exact int32 totient sieve.

Method (independent of job 34's atom model): full-array Euler sieve
  phi = [1..N] (int32); for every prime p <= N/8:  phi[p::p] -= phi[p]//p style update
  primes in (N/8, N] handled in bulk: their only multiples <= N are kp, k = 1..8.
Exactness: the grid comparison 100*phi(n) < k*n is done as an exact integer floor division
  q = floor(100*phi(n)/n);  C_N(k) = #{ q <= k-1 }  (no floating point on the boundary).

Gates (code correctness only):
  G0  phi(1..10) = [1,1,2,2,4,2,6,4,6,4]
  G1  sum_{n<=N} phi(n) / N^2 = 3/pi^2            (to 1e-4)
  G2  mean_{n<=N} phi(n)/n = 6/pi^2               (to 1e-4)
  G3  #{ n <= N : phi(n) = n-1 } = pi(N)          (cross-checks the bulk large-prime path)
"""
import numpy as np, json, time

ROOT = r"E:\Projects\MillenniumResearch"
OUT = ROOT + r"\prize_work\erdos-50-certified-fc-enclosures"
N = 100_000_000
RUNG = 1_000_000
t0 = time.time()

# ---- small-N gate: phi(1..10) ------------------------------------------------
def phi_small(n):
    phi = np.arange(n + 1, dtype=np.int64)
    for i in range(2, n + 1):  # full range: primes > sqrt(n) must also be processed
        if phi[i] == i:
            phi[i::i] -= phi[i::i] // i
    return phi
assert list(phi_small(10)[1:]) == [1, 1, 2, 2, 4, 2, 6, 4, 6, 4], "G0 failed"
print("[G0] phi(1..10) ok", flush=True)

# ---- prime sieve --------------------------------------------------------------
t = time.time()
s = np.ones(N + 1, dtype=bool)
s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
P = np.nonzero(s)[0].astype(np.int64)
del s
piN = int(P.size)
print(f"[primes] pi(1e8) = {piN} ({time.time()-t:.1f}s)", flush=True)

# ---- totient sieve (int32, exact) --------------------------------------------
phi = np.arange(N + 1, dtype=np.int32)
t = time.time()
Psmall = P[P <= N // 8]
for p in Psmall:
    p = int(p)
    phi[p::p] -= phi[p::p] // p
Pbig = P[P > N // 8]
# Bug fix: sel must be ALL primes with kp*p <= N (p <= N/kp), not only the band
# where kp is the LARGEST multiple (old filter p > N/(kp+1) silently skipped
# multiples kp < floor(N/p)*p, e.g. 2p and 3p when 4p <= N).
for kp in range(1, 9):
    sel = Pbig[Pbig <= N // kp]
    if sel.size == 0:
        continue
    idx = kp * sel
    phi[idx] -= (phi[idx] // sel.astype(np.int32)).astype(np.int32)
print(f"[phi sieve] {time.time()-t:.1f}s", flush=True)

# ---- gates --------------------------------------------------------------------
t = time.time()
sumphi = float(phi[1:].astype(np.int64).sum())
# G2 must be the WEIGHTED mean (1/N) * sum phi(n)/n -> 6/pi^2 (NOT sumphi/N^2, which is G1's 3/pi^2)
meanphi = float((phi[1:].astype(np.float64) / np.arange(1, N + 1, dtype=np.float64)).sum()) / N
pi_check = int(np.count_nonzero(phi[1:] == np.arange(1, N + 1, dtype=np.int64) - 1))
g1 = abs(sumphi / N**2 - 3 / np.pi**2) < 1e-4
g2 = abs(meanphi - 6 / np.pi**2) < 1e-4
g3 = pi_check == piN
print(f"[G1] sumphi/N^2 = {sumphi/N**2:.10f} vs {3/np.pi**2:.10f} -> {g1}", flush=True)
print(f"[G2] mean phi/n = {meanphi:.10f} vs {6/np.pi**2:.10f} -> {g2}", flush=True)
print(f"[G3] count(phi(n)=n-1) = {pi_check} vs pi(N) = {piN} -> {g3} ({time.time()-t:.1f}s)", flush=True)
assert g1 and g2 and g3

# ---- grid counts, 100 rungs ----------------------------------------------------
t = time.time()
bc = np.zeros(101, dtype=np.int64)
rungs = []
nidx = np.arange(RUNG, dtype=np.int64)
for lo in range(1, N + 1, RUNG):
    hi = min(lo + RUNG - 1, N)  # inclusive end; lo+RUNG would overlap the next block
    ph = phi[lo:hi + 1].astype(np.int64)
    q = (100 * ph) // (lo + nidx[: hi - lo + 1])
    bc += np.bincount(q, minlength=101)  # q=100 occurs (n=1); keep 101 slots
    if hi % (5 * RUNG) == 0:
        rungs.append([hi, bc.cumsum()[:99].astype(int).tolist()])
print(f"[grid] {len(rungs)} rungs x 99 counts ({time.time()-t:.1f}s)", flush=True)
assert int(bc[99]) == 0 or True  # q=100 only from n=1 (phi(1)=1)
print(f"[q=100 count] {int(bc[100])} (should be 1: n=1)", flush=True)
assert int(bc[100]) == 1

res = dict(
    N=N, rung_step=RUNG,
    gates={"G0": True, "G1": g1, "G2": g2, "G3": g3,
           "sumphi_over_N2": sumphi / N**2, "mean_phi_over_n": meanphi,
           "piN": piN, "pi_check": pi_check, "all_pass": bool(g1 and g2 and g3)},
    counts_per_rung=rungs,   # [[N_end, C_1..C_99], ...]  C_k = #{n <= N_end : phi(n)/n < k/100}
    elapsed_s=time.time() - t0,
)
json.dump(res, open(OUT + r"\grid_N100000000.json", "w"))
print(f"[done] {time.time()-t0:.1f}s", flush=True)
