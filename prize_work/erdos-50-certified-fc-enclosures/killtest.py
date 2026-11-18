"""Kill test for card erdos-50-certified-fc-enclosures (job lab_erdos_261113_43).

CARD CLAIM: f(c) = (exact finite atom sum over squarefree prime sets U with c_U < c)
                     + (truncation tail bounded by the Weingartner tail estimate),
so that f(c) can be enclosed between two exact rationals at c = k/100.

KILL CRITERION (pre-registered in the card): at c = 1/2 the tail bound must be SMALLER
than the gap between f(1/2) at N=1e7 and N=1e8.  If not, the enclosure is vacuous.

This script measures, with cached exact totients (no new sieve):
  (1) the gap  D(N1,N2) = f_N2(1/2) - f_N1(1/2)  at the two N we can afford cheaply;
  (2) the truncation remainder that the card's "tail" has to bound:
        T(M;N,c) = #{n <= N : phi(n)/n < c, rad(n) > M} / N,
      i.e. exactly the mass of the atoms the truncation at radical size M drops,
      as a function of M (M = 1e4 ... N/2).  T(M) is the SMALLEST interval width the
      card's scheme can have at truncation M, since it is a measured, not a bounded,
      quantity -- and any *bound* on it is >= it.
Gate: sum_{n<=N} phi(n)/N^2 = 3/pi^2 (code correctness only).
"""
import numpy as np, json, time

ROOT = r"E:\Projects\MillenniumResearch"
OUT = ROOT + r"\prize_work\erdos-50-certified-fc-enclosures"
t0 = time.time()

phi1 = np.load(ROOT + r"\prize_work\ppl-003\phi_cache_10000000.npy")   # phi(n), n = index
phi2 = np.load(ROOT + r"\prize_work\ppl-003\phi_cache_20000000.npy")
N1, N2 = 10_000_000, 20_000_000
gate = {}
gate["N1"] = dict(sumphi_over_N2=float(phi1.sum()) / N1**2, target=3 / np.pi**2)
gate["N2"] = dict(sumphi_over_N2=float(phi2.sum()) / N2**2, target=3 / np.pi**2)
for k, v in gate.items():
    v["pass"] = abs(v["sumphi_over_N2"] - v["target"]) < 1e-4
print("[gate]", gate, flush=True)

GRID = [(1, 3), (1, 2), (2, 3), (3, 4)]
fN = {}
for name, ph, N in (("N1e7", phi1, N1), ("N2e7", phi2, N2)):
    n = np.arange(1, N + 1, dtype=np.int64)
    row = {}
    for (p, q) in GRID:
        row[f"{p}/{q}"] = int(np.count_nonzero(q * ph[1:] < p * n))
    fN[name] = row
    print(f"[f_N {name}]", row, f"({time.time()-t0:.1f}s)", flush=True)
    del n

gap_half = (fN["N2e7"]["1/2"] / N2) - (fN["N1e7"]["1/2"] / N1)
gap_third = (fN["N2e7"]["1/3"] / N2) - (fN["N1e7"]["1/3"] / N1)
print(f"[gap 1/2] {gap_half:.8g}   [gap 1/3] {gap_third:.8g}", flush=True)

# --- truncation remainder T(M;N,1/2) at N=2e7 -------------------------------
N = N2
primes = np.nonzero(np.ones(N + 1, dtype=bool))[0][:0]  # placeholder
s = np.ones(N + 1, dtype=bool)
s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
P = np.nonzero(s)[0]
print(f"[primes] {len(P)} ({time.time()-t0:.1f}s)", flush=True)
del s
rad = np.ones(N + 1, dtype=np.int64)
rad[0] = 0
t = time.time()
for p in P:
    p = int(p)
    rad[p::p] *= p
print(f"[rad sieve] {time.time()-t:.1f}s", flush=True)

n = np.arange(1, N + 1, dtype=np.int64)
below = (2 * phi2[1:] < n)          # phi(n)/n < 1/2, exact integer test
tot = int(np.count_nonzero(below))
print(f"[f_N(1/2) total] {tot} / {N} = {tot/N:.8f}", flush=True)
tail = {}
for M in (10**4, 3 * 10**4, 10**5, 3 * 10**5, 10**6, 3 * 10**6, 10**7):
    c = int(np.count_nonzero(below & (rad[1:] > M)))
    tail[M] = dict(count=c, density=c / N, ratio_to_gap=(c / N) / abs(gap_half))
    print(f"[tail] M={M:>9} count={c:>10} density={c/N:.6f} "
          f"density/|gap|={(c/N)/abs(gap_half):.1f}  ({time.time()-t0:.0f}s)", flush=True)

res = dict(gate=gate, fN=fN, gap_half=gap_half, gap_third=gap_third,
           tail_half=tail, N=N, elapsed=time.time() - t0,
           verdict="KILL FIRED" if min(v["ratio_to_gap"] for v in tail.values()) > 1 else "not fired",
           min_ratio=min(v["ratio_to_gap"] for v in tail.values()))
json.dump(res, open(OUT + r"\killtest.json", "w"), indent=1)
print("[verdict]", res["verdict"], "min ratio to gap:", res["min_ratio"], flush=True)
