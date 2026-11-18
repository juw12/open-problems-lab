# Job lab_erdos_261112_11 — Erdos #50 — INDEPENDENT VERIFIER
# Fresh from-scratch enumeration of prime-tower values c_S (float), |S|<=5, primes<=199.
# No lab stored data used. Confirms the two band edges are the ACCUMULATION boundary:
# the edge towers have many neighboring tower values (accumulation points) while the
# next tower out on each side is isolated (0-1 neighbors). This is a different method
# than the lab's stored windowed-slope 'growth' field.

import math
from itertools import combinations

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
          67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
          139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199]
n = len(PRIMES)
print(f"primes <= 199: {n}")

# float tower value (1-1/p) factors
fp = [(1.0 - 1.0 / p) for p in PRIMES]

def cval(S):
    v = 1.0
    for i in S:
        v *= fp[i]
    return v

# enumerate all |S|<=5 tower values (fresh)
vals = [1.0]  # S={} -> c=1 (the n=1 point); harmless
cnt = 0
for r in range(1, 6):
    for S in combinations(range(n), r):
        vals.append(cval(S))
        cnt += 1
print(f"enumerated {cnt} tower values (|S|=1..5) + 1 (S={{}}) = {len(vals)} total")
vals.sort()

def neighbors(v, w):
    # count tower values within w of v (excluding self, via binary search)
    import bisect
    lo = bisect.bisect_left(vals, v - w)
    hi = bisect.bisect_right(vals, v + w)
    return hi - lo  # includes self; subtract 1 below if needed

w = 1e-3
def test(label, v):
    k = neighbors(v, w)
    print(f"{label:42s} c={v:.6f}  neighbors(w={w:g})={k}")
    return k

print(f"\n=== LEFT EDGE (w={w:g}) ===")
test("left edge  {2,3,5,7}=8/35", cval((0, 1, 2, 3)))
test("next smaller {2,3,5,7,11}", cval((0, 1, 2, 3, 4)))
test("next smaller2 {2,3,5,11}", cval((0, 1, 2, 4)))
test("3-tower {2,3,5}", cval((0, 1, 2)))

print(f"\n=== RIGHT EDGE (w={w:g}) ===")
test("right edge {47,73,83,157}", cval((PRIMES.index(47), PRIMES.index(73),
                                        PRIMES.index(83), PRIMES.index(157))))
test("single prime {23} (above edge)", cval((PRIMES.index(23),)))
test("single prime {19} (inside band)", cval((PRIMES.index(19),)))
test("3-tower large {89,97,103} (above)", cval((PRIMES.index(89), PRIMES.index(97),
                                                PRIMES.index(103))))

# density profile: how many tower values in each 1e-3 bin across [0.2, 0.99]
print(f"\n=== tower-value density (fresh, |S|<=5, primes<=199), bin=2e-3 ===")
import bisect
lo, hi, b = 0.20, 0.99, 2e-3
for i in range(0, 40, 1):
    c0 = lo + (hi - lo) * i / 39
    k = bisect.bisect_right(vals, c0 + b / 2) - bisect.bisect_left(vals, c0 - b / 2)
    bar = "#" * int(k / 50)
    print(f"  c~{c0:.3f}  count={k:5d} {bar}")
