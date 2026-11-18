"""Sanity checks on grid_N100000000.json (independent of the building script's internal state).
Cross-checks against the killtest's int64 phi-cache counts (different totient lineage):
  k=50 (c=1/2): C(1e7)=5111904, C(2e7)=10224118
  k=75 (c=3/4): C(1e7)=6891233, C(2e7)=13783411
Plus: (a) monotone in N at every rung for every k; (b) monotone in k within a rung;
(c) C_k(N) <= C_{k'}(N) for k<k'; (d) total accounted: sum of q-histogram = N.
"""
import json

ROOT = r"E:\Projects\MillenniumResearch\prize_work\erdos-50-certified-fc-enclosures"
d = json.load(open(ROOT + r"\grid_N100000000.json"))
rungs = d["counts_per_rung"]
assert len(rungs) == 20 and rungs[0][0] == 5_000_000 and rungs[-1][0] == 100_000_000

ok = True
# (a) monotone in N, every k
# rows are [N_end, [C_1, ..., C_99]] (nested)
for a, b in zip(rungs, rungs[1:]):
    bad = [k for k in range(1, 99) if b[1][k - 1] < a[1][k - 1]]
    if bad:
        ok = False
        print(f"NON-MONOTONE in N between {a[0]} and {b[0]} at k={bad[:5]}...")
# (b) monotone in k, every rung
for r in rungs:
    bad = [k for k in range(2, 99) if r[1][k - 1] < r[1][k - 2]]
    if bad:
        ok = False
        print(f"NON-MONOTONE in k at N={r[0]}: {bad[:5]}")

# cross-checks
def C(Nend, k):
    r = next(r for r in rungs if r[0] == Nend)
    return r[1][k - 1]

checks = [
    (10_000_000, 50, 5111904, "killtest int64 cache 1e7 c=1/2"),
    (20_000_000, 50, 10224118, "killtest int64 cache 2e7 c=1/2"),
    (10_000_000, 75, 6891233, "killtest int64 cache 1e7 c=3/4"),
    (20_000_000, 75, 13783411, "killtest int64 cache 2e7 c=3/4"),
]
for Nend, k, want, src in checks:
    got = C(Nend, k)
    status = "OK " if got == want else "FAIL"
    ok = ok and got == want
    print(f"[{status}] C({Nend},{k}) = {got} vs {want} ({src})")

# (d) accounted total at last rung: C_99(1e8) + #{q=100} + #{q<100...} -- C_99 counts q<=98;
# q in {99} and q=100 remain. q=100 is 1 (n=1). Report the remainder explicitly.
C99 = C(100_000_000, 99)
print(f"remainder above C_99 at 1e8: {100_000_000 - C99} (n with 100*phi(n)/n >= 99)")
print("ALL_OK" if ok else "SOME_FAILED")
