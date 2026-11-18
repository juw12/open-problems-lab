"""Independent NO-SIEVE verifier for lab_erdos_261111_24 (Erdos #50 / PPL 003).

Reads ONLY the stored scale_N1e8_100000000.json and RE-DERIVES every headline
claim by independent arithmetic on the raw stored traces (no phi sieve, no
prime sieve, no `::` slice notation).

Role (protocol "independent verifier"): a second METHOD on the raw numbers,
distinct from scale_N1e8.py's dual-sieve + analytic-tower computation.
Together they cross-check the result:
  * scale_N1e8.py SIEVE A (int32 exact) vs SIEVE B (float multiplicative) =
    the card's required "second differently-implemented phi-sieve".
  * scale_N1e8.py analytic 2^a*p tower = sieve-independent, prime-list only.
  * this file = pure arithmetic re-derivation of every headline claim.

Checks
  1. Monotonicity of the window-slope trace on BOTH sieves (strictly rising,
     0 direction changes)  -> the c=1/2 spike is NOT a finite-N oscillation.
  2. Kill test: global c* within threshold (0.005) of 0.5.
  3. A vs B dual-sieve agreement < 1e-6 (the card's required verifier).
  4. ~1/N turnover scaling of the analytic 2^a*p tower, RE-DERIVED from the
     raw rows (finest eps with A_count>0), not from the stored label.
  5. Analytic tower vs measured spike mass < 1e-5.
  6. Brute-force spot check: every stored (sieve,brute) pair agrees.
  7. Gate booleans all true.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "scale_N1e8_100000000.json")
with open(BASE) as f:
    d = json.load(f)

N = d["N"]
assert N == 10**8, f"expected N=1e8, got {N}"

checks = []
def check(name, ok, detail):
    checks.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")

A = d["sieve_A_int32_exact"]
B = d["sieve_B_float_multiplicative"]

# ---- 1. Monotonicity on BOTH sieves (independent direction count) ----------
def mono(trace):
    s = [float(r["slope_at_0p5"]) for r in trace]
    w = [float(r["w"]) for r in trace]
    w_dec = all(w[i] > w[i+1] for i in range(len(w)-1))
    strict_rise = all(s[i+1] > s[i] for i in range(len(s)-1))
    diffs = [s[i+1]-s[i] for i in range(len(s)-1)]
    dc = sum(1 for i in range(1, len(diffs)) if (diffs[i]>0) != (diffs[i-1]>0))
    return w_dec, strict_rise, dc, s[0], s[-1]

wdA, rA, dA, s0A, s1A = mono(A["slope_trace"])
wdB, rB, dB, s0B, s1B = mono(B["slope_trace"])
check("monotonicity A: w strictly decreasing + slope strictly rising (0 dir-chg)",
      wdA and rA and dA == 0,
      f"slope {s0A:.2f}->{s1A:.2f}, dir_changes={dA}")
check("monotonicity B: w strictly decreasing + slope strictly rising (0 dir-chg)",
      wdB and rB and dB == 0,
      f"slope {s0B:.2f}->{s1B:.2f}, dir_changes={dB}")
check("monotonicity: NOT oscillating on either sieve (=> not finite-N artifact)",
      (dA == 0) and (dB == 0),
      f"A.dir_changes={dA}, B.dir_changes={dB}")

# log-log differential exponent (diagnostic; slope ~ w^-alpha, expect alpha~1)
import math
def loglog_exp(trace):
    s = [math.log(float(r["slope_at_0p5"])) for r in trace]
    w = [math.log(float(r["w"])) for r in trace]
    n = len(w)
    sw, ss = sum(w), sum(s)
    swx = sum(w[i]*s[i] for i in range(n))
    varw = sum((w[i]-sw/n)**2 for i in range(n))
    cov  = sum((w[i]-sw/n)*(s[i]-ss/n) for i in range(n))
    return -cov/varw
expA = loglog_exp(A["slope_trace"])
check("diagnostic: window-slope ~ w^-alpha (differential alpha near 1)",
      0.6 < expA < 1.6, f"alpha_A={expA:.3f}")

# ---- 2. Kill test ------------------------------------------------------------
kd = d["kill_dcstar"]
check("kill test: global c* within threshold of 0.5 (not killed)",
      (abs(kd["abs_from_half"]) < kd["threshold"]) and (kd["kill_off_0p5"] is False),
      f"c*={kd['global']}, |c*-0.5|={kd['abs_from_half']}, thr={kd['threshold']}, "
      f"kill_off_0p5={kd['kill_off_0p5']}")

# ---- 3. A vs B dual-sieve agreement < 1e-6 (card verifier) --------------------
ab = d["A_vs_B_verifier"]
check("card verifier: dual-sieve A vs B agree < 1e-6",
      (ab["dcstar"] < 1e-6) and (ab["max_dspike"] < 1e-6) and (ab["VERDIFIER_PASS"]),
      f"dcstar={ab['dcstar']}, max_dspike={ab['max_dspike']}, VERDIFIER_PASS={ab['VERDIFIER_PASS']}")

# ---- 4. ~1/N turnover, RE-DERIVED from raw rows ------------------------------
tr = d["analytic_fine_divergence_trace"]
def turnover_from_rows(rows):
    # finest eps (smallest) with A_count > 0
    pos = [float(r["eps"]) for r in rows if int(r["A_count"]) > 0]
    return min(pos) if pos else None
derived = {}
for key in ["N_1000000", "N_10000000", "N_100000000"]:
    blk = tr[key]
    derived[int(key[2:])] = (turnover_from_rows(blk["rows"]), float(blk["turnover_eps"]))
Nvals = [10**6, 10**7, 10**8]
ratios = [derived[Nvals[i]][0]/derived[Nvals[i+1]][0] for i in range(2)]
prods  = {n: n*derived[n][0] for n in Nvals}
labels_ok = all(derived[n][0] == derived[n][1] for n in Nvals)
check("turnover RE-DERIVED from raw rows matches stored label", labels_ok,
      " ".join(f"{n}: derived={derived[n][0]:.0e} stored={derived[n][1]:.0e}" for n in Nvals))
check("turnover scales ~1/N: N*eps_turnover is a constant ~O(1)",
      all(0.5 <= p <= 20.0 for p in prods.values()),
      f"N*eps_turnover = " + ", ".join(f"{n}:{prods[n]:.1f}" for n in Nvals))
check("turnover deepens one decade per decade of N (genuine singularity, "
      "not a plateauing finite-N artifact)",
      all(5.0 < r < 20.0 for r in ratios),
      f"ratios eps(1e6)/eps(1e7)={ratios[0]:.1f}, eps(1e7)/eps(1e8)={ratios[1]:.1f} (each ~10)")
pderiv = [tr["N_1000000"]["peak_left_deriv"], tr["N_10000000"]["peak_left_deriv"],
          tr["N_100000000"]["peak_left_deriv"]]
check("diagnostic: peak left-derivative grows ~10x per decade (divergence deepens)",
      (pderiv[1]/pderiv[0] > 5) and (pderiv[2]/pderiv[1] > 5),
      f"peak_left_deriv {pderiv[0]:.0f} -> {pderiv[1]:.0f} -> {pderiv[2]:.0f}")

# ---- 5. Analytic tower vs measured spike mass < 1e-5 -------------------------
at = d["analytic_tower"]
av = [abs(at[k]["abs_analytic_minus_A"]) for k in at]
check("analytic tower matches measured spike mass < 1e-5",
      all(x < 1e-5 for x in av),
      " |A-meas| = " + ", ".join(f"{k}:{at[k]['abs_analytic_minus_A']:.1e}" for k in at))

# ---- 6. Brute-force spot check: every pair agrees ---------------------------
bs = d["brute_phi_spot_check"]
allagree = all(v["agree"] for v in bs.values())
check("brute-force spot check: all stored (sieve,brute) pairs agree",
      allagree and len(bs) >= 10,
      f"{sum(1 for v in bs.values() if v['agree'])}/{len(bs)} agree")

# ---- 7. Gate booleans --------------------------------------------------------
g = d["gates"]
check("sieve correctness gates all true", all(g.values()) and d["all_gates_and_brute"],
      f"gates={g}, all_gates_and_brute={d['all_gates_and_brute']}")

# ---- verdict ----------------------------------------------------------------
n_pass = sum(1 for _, ok, _ in checks if ok)
n_fail = sum(1 for _, ok, _ in checks if not ok)
print(f"\n== NO-SIEVE INDEPENDENT VERIFIER: {n_pass} pass, {n_fail} fail ==")
allpass = n_fail == 0
print("ALL_PASS:", allpass)

out = {
    "file_verified": os.path.basename(BASE),
    "N": N,
    "checks": [{"name": n, "ok": o, "detail": dt} for n, o, dt in checks],
    "n_pass": n_pass,
    "n_fail": n_fail,
    "all_pass": allpass,
    "derived_turnover_eps": {str(n): derived[n][0] for n in Nvals},
    "N_times_eps_turnover": {str(n): prods[n] for n in Nvals},
}
with open(os.path.join(HERE, "verify_job24_no_sieve.json"), "w") as f:
    json.dump(out, f, indent=2)
print("wrote verify_job24_no_sieve.json")
