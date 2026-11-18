# Test the j261117_5 criterion against this job's independently computed 198 rows:
#   cap* = cap - 1  iff  n even and S_pl = 2 (mod 4),
#   S_pl = 2*floor(S*_real/2) (the even-rounded DPTV79 optimum).
import json, io
from fractions import Fraction
t = json.load(io.open("cap_table.json", encoding="utf-8"))
w = {r["n"]: r for r in t["wide_rows"]}
def Spl(n):
    E = n * 2 ** (n - 1); m = 2 ** (n - 1)
    start = 0 if n % 2 == 0 else 1
    best = None
    for a in range(start, n + 1, 2):
        v = Fraction(E + a * (a + 2) * m, 2 * a + 2)
        if best is None or v < best: best = v
    return 2 * int(best // 2)
bad = []; rows = 0
for n in sorted(w):
    pred = (n % 2 == 0) and (Spl(n) % 4 == 2)
    act = w[n]["cap_plain"] - w[n]["cap_quant"] == 1
    rows += 1
    if pred != act: bad.append((n, pred, act, Spl(n) % 4))
print("rows", rows, "mismatches", bad)
print("sharpened predicted:", [n for n in sorted(w) if (n % 2 == 0) and Spl(n) % 4 == 2])
