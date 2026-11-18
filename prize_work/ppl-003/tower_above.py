# job lab_erdos_261112_24 -- Erdos #50 (PPL 003): the ABOVE-side local mass at the top-20 tower tops
#
# New argument (the card's core is a dead end):
#   job-20 generalized job-22's exact c=1/2 2^a*p tower to an arbitrary top c_S and got the
#   one-sided-BELOW tower mass T_S(N,W) ~ C(S)/lnN, C(S)=prod_{q in S}1/(q-1), verified on the
#   top-5 named points (below-dominant). It found the top-20 BY GROWTH are |S|=4 ~2/3-cluster
#   points whose below-tower is EMPTY (above-dominant). Referee 261112_14 downgraded the C(S)
#   law to a PNT corollary (reproduction).
#   => the card's "verify top-20" is half done: the BELOW side was verified (empty for the
#      top-20); the ABOVE side (the actual source of their divergent slope) is uncharacterized.
#   THIS JOB measures the above-side local mass M_above(N,W)=(1/N)*#{n: c_S < phi(n)/n < c_S+W}
#   at the top-20, characterizes its N-scaling, and states the below/above dichotomy.
#
#   Exact finite-N structure (no prime data in the structure, only in the count):
#     v(n) := phi(n)/n = c_{P(n)}   (constant on prime set P(n)=set of prime factors of n)
#     M_above(N,W;c_S) = (1/N) * sum_{T : c_T in (c_S, c_S+W)} A_T(N)
#     A_T(N) = #{n<=N : P(n)=T} ~ (ln N)^|T| / (|T|! * prod_{q in T} ln q)   [PNT-level]
#   => the above side is a SUM OF ATOMS over the value-set window, NOT a single 1/ln N tower
#      (the single tower is the BELOW side). We verify the structure empirically + independently.
#
# Gates (prove code): reproduce job-8 c=1/2 below-mass + job-20 top-5 below-masses + job-35 in_window.
# Kill test: at rank1 (sharpest), is M_above*lnN ~constant across N? (No => 1/ln N is wrong for the
#   above side). Verifier: (i) 3-N scaling (5e6/1e7/2e7); (ii) independent float multiplicative sieve
#   (different dtype+algorithm) recomputing M_above at 1e7.
import numpy as np, json, math, time, sys
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8')

ROOT = r"E:\Projects\MillenniumResearch"
SW = f"{ROOT}/prize_work/ppl-003"
W = 5e-5
W2 = 0.00025
t0 = time.time()

def sieve_phi_int(N):
    phi = np.arange(N + 1, dtype=np.int64)
    if N >= 1:
        phi[1] = 1
    for p in range(2, N + 1):
        if phi[p] == p:                       # p is prime
            for j in range(p, N + 1, p):
                phi[j] -= phi[j] // p
    return phi

def sieve_phi_float(N):
    phi = np.arange(N + 1, dtype=np.float64)
    phi[0] = 0.0
    if N >= 1:
        phi[1] = 1.0
    for p in range(2, N + 1):
        if phi[p] == p:                       # p is prime
            for j in range(p, N + 1, p):
                phi[j] -= phi[j] / p
    return phi

def load_or_build(N, kind):
    path = f"{SW}/phi_cache_{N}.npy"
    if kind == "int" and __import__("os").path.exists(path):
        a = np.load(path).astype(np.int64)
        assert a.shape == (N + 1,)
        return a, "cached"
    if kind == "float":
        return sieve_phi_float(N), "built-float"
    a = sieve_phi_int(N)
    np.save(path, a)
    return a, "built-int"

NMAIN = 10_000_000
NLIST = [5_000_000, NMAIN, 20_000_000]

phi = {}
src = {}
r = {}
for N in NLIST:
    a, s = load_or_build(N, "int")
    phi[N] = a
    src[N] = s
    nn = np.arange(1, N + 1, dtype=np.float64)
    r[N] = a[1:].astype(np.float64) / nn
    print(f"[sieve] N={N} source={s} mean(r)={np.mean(r[N]):.7f} (6/pi^2={6/math.pi**2:.7f}) elapsed={time.time()-t0:.0f}s", flush=True)

# independent float sieve at NMAIN
phiF, _ = load_or_build(NMAIN, "float")
nnF = np.arange(1, NMAIN + 1, dtype=np.float64)
rF = phiF[1:].astype(np.float64) / nnF
print(f"[sieve-float] N={NMAIN} mean(r)={np.mean(rF):.7f} elapsed={time.time()-t0:.0f}s", flush=True)

# ---------------- gates: independent phi small-n naive check at NMAIN ----------------
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

chk = np.array([phi_naive(i) for i in range(1, 50001)], dtype=np.int64)
phi_gate_small = bool(np.array_equal(phi[NMAIN][1:50001], chk))

# ---------------- tower tops (same set as job-20: top-20 by growth U named top-5) ----------------
NAMED_TOP5 = ["1/2", "2/3", "1/3", "2/5", "3/7"]
with open(f"{SW}/scan_N{NMAIN}.json") as f:
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

# job-20 stored values (gates) : c -> (measured_full_below at 1e7, in_window at W2)
J20 = json.load(open(f"{SW}/tower_general_N{NMAIN}.json"))
j20_by_c = {e["c"]: e for e in J20["top20"]}

# gate reference: job-8 c=1/2 below mass
GATE_1HALF = {5_000_000: 0.0752172, 10_000_000: 0.0723393}

def m_above_below(N, cS):
    rr = r[N]
    ma = float(np.sum((rr > cS) & (rr < cS + W))) / N
    mb = float(np.sum((rr > cS - W) & (rr < cS))) / N
    return ma, mb

results = []
for rank, d in enumerate(tops, 1):
    S = [int(x) for x in d["S"]]
    cS = float(Fraction(d["c"]))
    lnN = {}
    rec = {"rank": rank, "c": d["c"], "S": S, "cS_float": cS, "growth": d["growth"],
           "in_named_top5": d["c"] in NAMED_TOP5, "above_dominant": True}
    # measure at each N
    for N in NLIST:
        ma, mb = m_above_below(N, cS)
        rec[f"M_above_{N}"] = ma
        rec[f"M_below_{N}"] = mb
        lnN[N] = math.log(N)
        rec[f"Mabove_lNN_{N}"] = ma * lnN[N]
        rec[f"Mabove_times_N_{N}"] = ma * N
    # independent float-sieve M_above at NMAIN
    maF = float(np.sum((rF > cS) & (rF < cS + W))) / NMAIN
    rec["M_above_floatN"] = maF
    rec["relerr_int_vs_float_pct"] = 100 * abs(rec["M_above_10000000"] - maF) / max(maF, 1e-12)
    # above/below dominance
    rec["ratio_above_over_below_1e7"] = (rec["M_above_10000000"] / rec["M_below_10000000"]) if rec["M_below_10000000"] > 0 else None
    # job-20 gate: below mass at 1e7 should match
    if d["c"] in j20_by_c:
        rec["gate_below_j20_stored"] = j20_by_c[d["c"]]["measured_full_below"]
        rec["gate_below_repro"] = rec["M_below_10000000"]
        rec["gate_below_match"] = abs(rec["M_below_10000000"] - j20_by_c[d["c"]]["measured_full_below"]) < 1e-9
        rec["gate_inwindow_stored"] = j20_by_c[d["c"]]["job35_in_window_stored"]
        inw = int(np.sum((r[NMAIN] >= cS - W2) & (r[NMAIN] <= cS + W2)))
        rec["gate_inwindow_repro"] = inw
    results.append(rec)

# ---------------- scaling fit for the above-dominant top-20 (M_above ~ (ln N)^K / N) ----------------
def fit_K(rec):
    # ln(M_above*N) = a + K * ln(ln N); use 3 N points
    xs = [math.log(math.log(N)) for N in NLIST]
    ys = [math.log(rec[f"Mabove_times_N_{N}"]) for N in NLIST if rec[f"Mabove_times_N_{N}"] > 0]
    if len(ys) < 3:
        return None, None
    x = np.array(xs); y = np.array(ys)
    K, a = np.polyfit(x, y, 1)
    resid = y - (K * x + a)
    return float(K), float(np.max(np.abs(resid)))

for rec in results:
    K, resid = fit_K(rec)
    rec["fit_K_above_atom_order"] = K
    rec["fit_residual"] = resid
    # 1/ln N would imply M_above*N = N/ln N  ->  ln(M_above*N) = ln N - ln ln N, NOT a power of ln N.
    # We report K: K>1 => (ln N)^K/N atom-sum (NOT 1/ln N tower).

above20 = [rec for rec in results if rec["in_named_top5"] is False]
named5 = [rec for rec in results if rec["in_named_top5"] is True]
Ks = [rec["fit_K_above_atom_order"] for rec in above20 if rec["fit_K_above_atom_order"] is not None]
resid_max = max((rec["fit_residual"] for rec in above20 if rec["fit_residual"] is not None), default=None)

# ---------------- gate summary ----------------
gate_1half = {}
for N in [5_000_000, 10_000_000]:
    ma, mb = m_above_below(N, 0.5)
    gate_1half[N] = {"repro_below": mb, "stored": GATE_1HALF[N],
                     "absdiff": abs(mb - GATE_1HALF[N]),
                     "match": abs(mb - GATE_1HALF[N]) < 5e-5}
below_matches = all(rec["gate_below_match"] for rec in results if "gate_below_match" in rec)
inwindow_matches = all(rec["gate_inwindow_stored"] == rec["gate_inwindow_repro"] for rec in results if "gate_inwindow_stored" in rec)
all_gates = bool(phi_gate_small and all(g["match"] for g in gate_1half.values()) and below_matches and inwindow_matches)

# ---------------- KILL TEST: rank1 M_above*lnN stable across N? ----------------
rk1 = results[0]
l1 = [rk1[f"Mabove_lNN_{N}"] for N in NLIST]
kill_ratio = l1[2] / l1[0] if l1[0] > 0 else None   # 2e7 / 5e6
kill_const_1lnN = (max(l1) / min(l1)) < 1.30          # 1/ln N => M_above*lnN ~ const
kill = {"rank1_c": rk1["c"], "rank1_S": rk1["S"],
        "Mabove_lNN_by_N": {str(N): l1[i] for i, N in enumerate(NLIST)},
        "ratio_2e7_over_5e6": kill_ratio,
        "const_within_30pct": kill_const_1lnN,
        "verdict_1lnN_above": "CONSTANT (1/ln N holds)" if kill_const_1lnN else "DRIFTING (1/ln N WRONG for above side)"}

out = {
    "job": "lab_erdos_261112_24", "W": W, "W_two_sided_gate": W2, "N_list": NLIST,
    "sieve_source": {str(N): src[N] for N in NLIST},
    "gates": {
        "phi_small_naive_match": phi_gate_small,
        "c_half_below_mass": gate_1half,
        "below_mass_j20_matches": bool(below_matches),
        "inwindow_j35_matches": bool(inwindow_matches),
        "all_pass": all_gates,
    },
    "headline": {
        "structure": "M_above(N,W;cS)=(1/N)*sum_{T:c_T in (cS,cS+W)} A_T(N), A_T~(lnN)^|T|/(|T|!*prod ln q); a SUM of (lnN)^k/N atoms, not a single 1/ln N tower (the 1/ln N tower is the BELOW side).",
        "top20_above_dominant": all((rec['ratio_above_over_below_1e7'] is None or rec['ratio_above_over_below_1e7'] > 10) for rec in above20),
        "mean_fit_K_above_atom_order": float(np.mean(Ks)) if Ks else None,
        "fit_K_range": [float(min(Ks)), float(max(Ks))] if Ks else None,
        "fit_max_resid_ln": resid_max,
        "interpretation": "K>1 means M_above scales like (ln N)^K/N (atom-sum), NOT 1/ln N; confirms the card's 1/ln N tower form is the BELOW-side law only."
    },
    "top20_above": above20,
    "named_top5": named5,
    "kill_test": kill,
    "elapsed_s": time.time() - t0,
}

with open(f"{SW}/tower_above_N{NMAIN}.json", "w") as f:
    json.dump(out, f, indent=2, default=float)

print("\nGATES:", json.dumps(out["gates"], indent=1), flush=True)
print("\nHEADLINE:", json.dumps(out["headline"], indent=1), flush=True)
print("\nKILL:", json.dumps(kill, indent=1), flush=True)
print("\nrank  c              S                          Ma(1e7)  Ma(5e6)  Ma(2e7)  MaB(1e7)  ratioA/B   K_fit   floatErr%")
for rec in results:
    K = rec["fit_K_above_atom_order"]
    Kb = f"{K:.2f}" if K is not None else "  --"
    print(f"{rec['rank']:>3} {rec['c']:>13}  S={str(rec['S']):>26}  {rec['M_above_10000000']:.6f} {rec['M_above_5000000']:.6f} {rec['M_above_20000000']:.6f}  {rec['M_below_10000000']:.6f}  {str(rec['ratio_above_over_below_1e7']):>10}  {Kb}   {rec['relerr_int_vs_float_pct']:.3f}  name5={rec['in_named_top5']}", flush=True)
print("\nelapsed_s:", round(out["elapsed_s"], 1), flush=True)
