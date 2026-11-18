# job lab_erdos_261112_20 — Erdos #50: generalize job-22 c=1/2 tower to arbitrary tower top c_S
# Derive the 1/ln N mass + C(S); verify top-20 mapped tower tops at N=1e7.
import numpy as np, json, math, time, sys
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8')

ROOT = r"E:\Projects\MillenniumResearch"
SW = f"{ROOT}/prize_work/ppl-003"
N = 10_000_000
W = 5e-5          # one-sided window (job-8 c=1/2 window)
W2 = 0.00025      # job-35 two-sided window (for in_window repro gate)
t_start = time.time()

# ---------- 1. exact totient sieve (load cached int64) + independent gate ----------
phi = np.load(f"{SW}/phi_cache_{N}.npy").astype(np.int64)
assert phi.shape == (N + 1,)

def phi_naive(n):
    res, x, p = n, n, 2
    while p * p <= x:
        if x % p == 0:
            while x % p == 0:
                x //= p
            res -= res // p
        p += 1
    if x > 1:
        res -= res // x
    return res

chk = np.array([phi_naive(i) for i in range(1, 100001)], dtype=np.int64)
phi_gate_small = bool(np.array_equal(phi[1:100001], chk))

def sieve_primes(Nn):
    isp = np.ones(Nn + 1, dtype=bool)
    isp[:2] = False
    for i in range(2, int(Nn ** 0.5) + 1):
        if isp[i]:
            isp[i * i::i] = False
    return isp

isp = sieve_primes(N)
pla = np.cumsum(isp).astype(np.int64)

def pi(x):
    x = int(x)
    if x < 2:
        return 0
    if x > N:
        x = N
    return int(pla[x])

n_arr = np.arange(1, N + 1, dtype=np.float64)
r_ratio = phi[1:].astype(np.float64) / n_arr

gates = {
    "phi_small_naive_match": phi_gate_small,
    "mean_ratio": float(np.mean(r_ratio)),
    "mean_ratio_expected_6pi2": 6 / math.pi ** 2,
    "sumphi_over_N2": float(np.sum(phi[1:].astype(np.int64))) / float(N * N),
    "sumphi_over_N2_expected_3pi2": 3 / math.pi ** 2,
}
gates["all_pass"] = (phi_gate_small
                     and abs(gates["mean_ratio"] - gates["mean_ratio_expected_6pi2"]) < 1e-5
                     and abs(gates["sumphi_over_N2"] - gates["sumphi_over_N2_expected_3pi2"]) < 1e-5)

# ---------- 2. tower tops to verify ----------
# "top-20 mapped tower tops" = top-20 by growth, UNION the card's named top-5
# (the named top-5 have nonempty below towers; the top-20-by-growth are dominated by
#  large-prime |S|=4 points whose below tower is EMPTY / above-dominant).
NAMED_TOP5 = ["1/3", "2/5", "3/7", "1/2", "2/3"]
with open(f"{SW}/scan_N{N}.json") as f:
    scan = json.load(f)
dl = scan["divergent_list"]
by_growth = sorted(dl, key=lambda d: -float(d["growth"]))
top20 = by_growth[:20]
by_c = {d["c"]: d for d in dl}
tops, seen = [], set()
for d in top20 + [by_c[c] for c in NAMED_TOP5 if c in by_c]:
    if d["c"] not in seen:
        seen.add(d["c"])
        tops.append(d)

def s_smooth_Ms(S, Mmax):
    Ms = []
    S = [int(q) for q in S]
    def rec(i, cur):
        if i == len(S):
            Ms.append(cur)
            return
        q = S[i]
        p = q
        while cur * p <= Mmax:
            rec(i + 1, cur * p)
            p *= q
    rec(0, 1)
    return Ms

def tower_A_primesum(S, cS):
    S = [int(q) for q in S]
    mS = 1
    for q in S:
        mS *= q
    rmin = int(cS / W)          # floor(cS/W)
    Mmax = int(N * W / cS)
    Ms = s_smooth_Ms(S, Mmax)
    cnt = 0
    for M in Ms:
        if M < mS:
            continue
        hi = N // M
        if hi <= rmin:
            continue
        cnt += pi(hi) - pi(rmin)
    return cnt

def tower_B_strip(S, cS):
    S = [int(q) for q in S]
    mS = 1
    for q in S:
        mS *= q
    rmin = cS / W
    kmax = N // mS
    k = np.arange(1, kmax + 1, dtype=np.int64)
    rem = k.copy()
    for q in S:
        while True:
            mask = (rem % q == 0)
            if not mask.any():
                break
            rem[mask] = rem[mask] // q
    in_range = (rem >= 2) & (rem <= N)
    mask = in_range & isp[rem] & (rem > rmin)
    for q in S:
        mask &= (rem != q)
    return int(mask.sum())

def D_S(S, cutoff=1_000_000):
    S = [int(q) for q in S]
    mS = 1
    for q in S:
        mS *= q
    Ms = s_smooth_Ms(S, cutoff)
    s = 0.0
    for M in Ms:
        if M >= mS:
            s += (math.log(M) + 1) / M
    return s

results = []
for rank, d in enumerate(tops, 1):
    S = [int(x) for x in d["S"]]
    cS = float(Fraction(d["c"]))
    mS = 1
    for q in S:
        mS *= q
    C_S = 1.0
    for q in S:
        C_S /= (q - 1)
    M_below = float(np.sum((r_ratio > cS - W) & (r_ratio < cS))) / N  # full one-sided-below mass
    A = tower_A_primesum(S, cS)   # derived tower mass (prime-sum, from pi)
    B = tower_B_strip(S, cS)      # measured tower mass (sieve S-strip)
    Ds = D_S(S)
    lnN = math.log(N)
    bare = C_S / lnN
    pnt2 = (C_S + Ds / lnN) / lnN
    Bm = B / N
    rec = {
        "rank": rank, "c": d["c"], "S": S, "m_S": mS, "cS_float": cS, "growth": d["growth"],
        "in_named_top5": d["c"] in NAMED_TOP5,
        "measured_full_below": M_below,
        "tower_A_count": A, "tower_A": A / N,
        "tower_B_count": B, "tower_B": Bm,
        "AB_integer_equal": (A == B),
        "tower_nonempty": (B > 0),
        "contamination_below_minus_tower": M_below - Bm,  # other prime sets with c_T in (c_S-W,c_S)
        "C_S": C_S, "D_S": Ds,
        "bare_C_over_lnN": bare, "pnt_2term": pnt2,
        "relerr_A_vs_B_pct": (100 * abs(A / N - Bm) / Bm) if Bm > 0 else None,
        "relerr_pnt2_vs_B_pct": (100 * abs(pnt2 - Bm) / Bm) if Bm > 0 else None,
    }
    inw_repro = int(np.sum((r_ratio >= cS - W2) & (r_ratio <= cS + W2)))
    rec["job35_in_window_stored"] = d["in_window"]
    rec["job35_in_window_repro"] = inw_repro
    rec["job35_gate_pass"] = (inw_repro == d["in_window"])
    results.append(rec)

# KILL TEST: card's named top-5 (all nonempty below towers); compare derived (pnt2, A) to measured tower mass B
t5 = [r for r in results if r["in_named_top5"]]
kill_vals = [r["relerr_pnt2_vs_B_pct"] for r in t5 if r["relerr_pnt2_vs_B_pct"] is not None]
kill_max_err = max(kill_vals) if kill_vals else None
kill_pass = all((r["relerr_pnt2_vs_B_pct"] is None or r["relerr_pnt2_vs_B_pct"] <= 10.0) for r in t5)
kill_AB_pass = all(r["AB_integer_equal"] for r in t5)

out = {
    "job": "lab_erdos_261112_20", "N": N, "W": W, "W_two_sided_gate": W2,
    "gates": gates,
    "all_gates_pass": bool(gates["all_pass"]
                           and all(r["AB_integer_equal"] for r in results)
                           and all(r["job35_gate_pass"] for r in results)),
    "headline": {
        "tower_mass_identity": "T_S(N,W)=(1/N)*sum_{M S-smooth, m_S<=M<=N*W/c_S}[pi(N/M)-pi(c_S/W)] ~ C(S)/lnN, C(S)=prod_{q in S}1/(q-1) (generalizes job-22 C({2})=1)",
        "C_S_named_top5": {r["c"]: r["C_S"] for r in t5},
        "recovers_job22_S2": (abs(next(r["C_S"] for r in results if r["S"] == [2]) - 1.0) < 1e-12),
    },
    "top20": results,
    "kill_test": {
        "note": "derived tower mass (pnt2 / A) vs measured tower mass B at the card's named top-5 (all nonempty below towers)",
        "points": [r["c"] for r in t5],
        "per_point_relerr_A_vs_B_pct": {r["c"]: r["relerr_A_vs_B_pct"] for r in t5},
        "per_point_relerr_pnt2_vs_B_pct": {r["c"]: r["relerr_pnt2_vs_B_pct"] for r in t5},
        "max_relerr_pnt2_vs_B_pct": kill_max_err,
        "AB_integer_equal_all_top5": bool(kill_AB_pass),
        "threshold_pct": 10.0,
        "pass": bool(kill_pass and kill_AB_pass),
    },
    "elapsed_s": time.time() - t_start,
}

with open(f"{SW}/tower_general_N{N}.json", "w") as f:
    json.dump(out, f, indent=2, default=float)

print("gates:", json.dumps(gates, indent=1))
print("all_gates_pass:", out["all_gates_pass"])
print("C_S named top5:", out["headline"]["C_S_named_top5"])
print("kill:", json.dumps(out["kill_test"], indent=1))
for r in results:
    print(f"  {r['c']:>10} S={r['S']} Mb={r['measured_full_below']:.6f} A={r['tower_A']:.6f} B={r['tower_B']:.6f} ABeq={r['AB_integer_equal']} nonempty={r['tower_nonempty']} contam={r['contamination_below_minus_tower']:.6f} C={r['C_S']:.4f} pnt2={r['pnt_2term']:.6f} errA%={r['relerr_A_vs_B_pct']} errPNT%={r['relerr_pnt2_vs_B_pct']} j35={r['job35_gate_pass']} name5={r['in_named_top5']}")
print("elapsed_s:", round(out["elapsed_s"], 1))
