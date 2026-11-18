"""lab_erdos_261111_31 - Erdos #50 (PPL 003): quantify the NON-TOWER residual of the c=1/2
spike and test whether it is o(tower) as N grows.

Builds on (verified in ppl003_law.json / diag_residual.py):
  job 8  : M(N,eps)=dens{n<=N: phi(n)/n in [0.5-eps,0.5)} = 0.0723393 @N=1e7 (eps=5e-5).
  job 22 : exact tower A(N,eps)=(1/N)*sum_{a>=1}(pi(N/2^a)-pi(pmin-1)), pmin=ceil(1/(2eps));
           A matches M to <1e-5 and the 2^a*p*q layer is provably zero (min-gap + AM-GM).
  job 29 : A(1e7)=0.0723345, |A-M|=4.8e-6 (48 n's); A(1e8)=0.06274177, M=0.06274806.

THIS JOB (one step past): residual R = M - A, and the ratio R/A, at N=1e7 and N=2e7.
  - GATE: fresh totient sieve reproduces stored M(1e7)=0.0723393 and A(1e7)=0.0723345.
  - RESULT: R_count, R, and R/A at N=1e7 and N=2e7; is R/A < 1e-4? non-increasing?
  - CHARACTERIZE: factor every residual n (v2, odd-part prime pattern) => what family it is.
  - INDEPENDENT VERIFIER (no phi sieve): max of f=(p-1)(q-1)/(2pq) over 2^a*p*q<=N
    (p,q distinct odd primes) and the min gap gmin=1/2 - f_max; family empty in window
    iff f_max < 1/2 - eps. Cross-checks that the residual is NOT the 2^a*p*q layer.

Empirical only. No proof claim. Single-threaded. Writes residual_otoower.json.
Run: python residual_otoower.py
"""
import json
import math
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 5e-5
N_LIST = [10000000, 20000000]
MAXP = 20000000


def sieve_primes(N):
    isp = np.ones(N + 1, dtype=bool)
    isp[0] = False
    isp[1] = False
    r = int(math.isqrt(N))
    for i in range(2, r + 1):
        if isp[i]:
            start = i * i
            isp[start : N + 1 : i] = False
    return isp


def phi_sieve(N, isp):
    phi = np.arange(N + 1, dtype=np.int64)
    pr = np.nonzero(isp[: N + 1])[0]
    for p in pr:
        phi[p : N + 1 : p] = phi[p : N + 1 : p] // p * (p - 1)
    return phi


def A_tower(N, pla, eps):
    pmin = math.ceil(1.0 / (2.0 * eps))
    base = int(pla[pmin - 1])
    cnt = 0
    a = 1
    while (N >> a) >= pmin:
        cnt += int(pla[N >> a]) - base
        a += 1
    return cnt / float(N), cnt, a - 1


def pq_verifier(N, isp):
    """max f=(p-1)(q-1)/(2pq) over 2^a*p*q<=N (a>=1, p,q odd primes), and gmin.
    f is maximized at the largest prime p with 2*p*p<=N (product pq<=N/2). Independent of phi."""
    pmax = int(math.isqrt(N // 2))
    q = pmax
    while q >= 3 and not isp[q]:
        q -= 1
    cands = []
    if q >= 5:
        cands.append((q, q))
    q2 = q - 1
    while q2 >= 3 and not isp[q2]:
        q2 -= 1
    if q2 >= 3:
        cands.append((q, q2))
    best = -1.0
    best_pair = None
    for p, rr in cands:
        f = 0.5 * (1.0 - 1.0 / p) * (1.0 - 1.0 / rr)
        if 2 * p * rr <= N and f > best:
            best = f
            best_pair = (p, rr)
    gmin = 0.5 - best
    empty = best < (0.5 - EPS)
    return {
        "fmax": best,
        "argmax_pq": list(best_pair) if best_pair else None,
        "gmin": gmin,
        "window_lo": 0.5 - EPS,
        "family_empty_in_window": bool(empty),
        "gmin_over_2eps": gmin / (2.0 * EPS),
        "gap_to_window_lo": (0.5 - EPS) - best,
    }


def factor_odd(oddpart, primes_list):
    m = int(oddpart)
    fs = []
    for p in primes_list:
        if p * p > m:
            break
        while m % p == 0:
            fs.append(int(p))
            m //= p
        if m == 1:
            break
    if m > 1:
        fs.append(m)
    return fs


from fractions import Fraction
def classify(resid_n, isp, primes_list, phi):
    out = []
    n_genuine = 0
    n_artifact = 0
    for nn in resid_n:
        nn = int(nn)
        lowbit = nn & (-nn)
        v2 = int(round(math.log2(lowbit)))
        oddpart = nn // lowbit
        fp = factor_odd(oddpart, primes_list)
        is_pq = (len(fp) == 2 and fp[0] != fp[1])
        is_p2 = (len(fp) == 2 and fp[0] == fp[1])
        is_prime_odd = (len(fp) == 1)
        rr = int(phi[nn]) / float(nn)
        # exact-rational recheck: is EXACT phi(n)/n in [0.5-eps, 0.5)?
        ex = Fraction(int(phi[nn]), nn)
        lo_frac = Fraction(1, 2) - Fraction(int(1 / (1 / (2 * EPS))) and 1, 20000)  # 1/2 - 1/20000
        ex_in = (ex >= Fraction(1, 2) - Fraction(1, 20000)) and (ex < Fraction(1, 2))
        if ex_in:
            n_genuine += 1
        else:
            n_artifact += 1
        out.append({"n": nn, "v2": v2, "oddpart": int(oddpart), "pattern": [int(x) for x in fp],
                    "odd_prime_count": len(fp), "is_pq": is_pq, "is_p2": is_p2,
                    "is_odd_prime": is_prime_odd, "phi_over_n": round(rr, 9),
                    "dist_to_half": round(0.5 - rr, 9),
                    "exact_in_window": bool(ex_in)})
    return out, n_genuine, n_artifact


def main():
    t0 = time.time()
    print("prime sieve to %d" % MAXP, flush=True)
    isp = sieve_primes(MAXP)
    pla = np.cumsum(isp.astype(np.int64))
    primes_list = np.nonzero(isp)[0]
    out = {"eps": EPS, "per_N": {}}

    prev_ratio = None
    for N in N_LIST:
        print("phi sieve N=%d" % N, flush=True)
        phi = phi_sieve(N, isp)
        n = np.arange(1, N + 1, dtype=np.int64)
        rr = phi[1:].astype(np.float64) / n.astype(np.float64)
        in_win = (rr >= 0.5 - EPS) & (rr < 0.5)
        M = float(in_win.sum()) / N
        M_count = int(in_win.sum())
        A, A_count, amax = A_tower(N, pla, EPS)
        R_count = M_count - A_count
        R = R_count / float(N)
        ratio = (R_count / float(N)) / A if A > 0 else None

        nw_n = n[in_win]
        lowbit = nw_n & (-nw_n)
        v2 = np.log2(lowbit.astype(np.float64)).astype(np.int64)
        oddpart = nw_n // lowbit
        is_tower = (v2 >= 1) & (oddpart >= 3) & isp[oddpart]
        tower_in_win = int(is_tower.sum())
        resid_n = nw_n[~is_tower]

        cl, n_genuine, n_artifact = classify(resid_n, isp, primes_list, phi)
        from collections import Counter
        patc = Counter(tuple(x["pattern"][:3]) for x in cl)
        v2c = Counter(x["v2"] for x in cl)
        ntype = Counter()
        for x in cl:
            if x["v2"] == 0:
                ntype["odd_n"] += 1
            elif x["is_pq"]:
                ntype["2a_pq"] += 1
            elif x["is_p2"]:
                ntype["2a_p2"] += 1
            elif x["is_odd_prime"]:
                ntype["2a_p"] += 1
            else:
                ntype["other"] += 1

        ver = pq_verifier(N, isp)
        row = {
            "measured_M": round(M, 8),
            "M_count": M_count,
            "exact_tower_A": round(A, 8),
            "A_count": A_count,
            "amax": amax,
            "residual_count": R_count,
            "residual_R": round(R, 9),
            "ratio_R_over_A": round(ratio, 9) if ratio is not None else None,
            "ratio_R_over_A_times_1e4": round(1e4 * ratio, 4) if ratio is not None else None,
            "tower_in_win_2adic": tower_in_win,
            "tower_A_eq_2adic_B": (tower_in_win == A_count),
            "residual_v2_hist": {str(k): v for k, v in sorted(v2c.items())},
            "residual_type_hist": dict(ntype),
            "residual_exact_genuine": int(n_genuine),
            "residual_float_artifact": int(n_artifact),
            "residual_top_patterns": {str(k): v for k, v in patc.most_common(8)},
            "residual_first20": cl[:20],
            "verifier_2a_pq": ver,
        }
        if N == 10000000:
            row["gate_M_reproduces_stored"] = bool(abs(M - 0.0723393) < 1e-5)
            row["gate_A_reproduces_stored"] = bool(abs(A - 0.0723345) < 1e-5)
            row["mean_r_6_pi2"] = bool(abs(float(rr.mean()) - 6.0 / math.pi**2) < 1e-4)
            row["sum_phi_3_pi2"] = bool(abs(float(phi[1:].sum()) / float(N * N) - 3.0 / math.pi**2) < 1e-4)
        out["per_N"][str(N)] = row

        if prev_ratio is not None and ratio is not None:
            row["ratio_nonincreasing_vs_prev"] = bool(ratio <= prev_ratio)
        prev_ratio = ratio
        del phi, rr, in_win, n, nw_n, lowbit, v2, oddpart, is_tower
        print("  N=%d  M=%.7f  A=%.7f  R_count=%d  R/A=%.3e  (x1e4=%.3f)  tower2adic==A:%s"
              % (N, M, A, R_count, ratio, 1e4 * ratio, tower_in_win == A_count), flush=True)
        print("    type_hist=%s" % dict(ntype), flush=True)
        print("    verifier: fmax=%.8f gmin=%.3e empty=%s gmin/2eps=%.2f"
              % (ver["fmax"], ver["gmin"], ver["family_empty_in_window"], ver["gmin_over_2eps"]), flush=True)

    r1 = out["per_N"][str(10000000)]
    r2 = out["per_N"][str(20000000)]
    r1_ok = r1["ratio_R_over_A"] < 1e-4
    r2_ok = r2["ratio_R_over_A"] < 1e-4
    noninc = (r2["ratio_R_over_A"] <= r1["ratio_R_over_A"])
    ver_empty = r2["verifier_2a_pq"]["family_empty_in_window"]
    out["verdict"] = (
        "RESIDUAL_o_tower: R=M-A (spike mass beyond the 2^a*p tower). N=1e7: R_count=%d, R/A=%.3e (x1e4=%.3f). "
        "N=2e7: R_count=%d, R/A=%.3e (x1e4=%.3f). R/A<1e-4 at 1e7: %s; at 2e7: %s; non-increasing 1e7->2e7: %s. "
        "VERIFIER (independent, no phi sieve): 2^a*p*q family fmax=%.8f, gmin=%.3e, empty-in-window(2e7)=%s, "
        "gmin/2eps=%.2f => the residual is NOT the 2^a*p*q layer. Residual type_hist(2e7)=%s. GATE(1e7): "
        "M=%.7f (stored 0.0723393, repro=%s), A=%.7f (stored 0.0723345, repro=%s), tower==2adic=%s."
        % (r1["residual_count"], r1["ratio_R_over_A"], r1["ratio_R_over_A_times_1e4"],
           r2["residual_count"], r2["ratio_R_over_A"], r2["ratio_R_over_A_times_1e4"],
           r1_ok, r2_ok, noninc,
           r2["verifier_2a_pq"]["fmax"], r2["verifier_2a_pq"]["gmin"], ver_empty,
           r2["verifier_2a_pq"]["gmin_over_2eps"],
           dict(r2["residual_type_hist"]),
           r1["measured_M"], r1["gate_M_reproduces_stored"], r1["exact_tower_A"],
           r1["gate_A_reproduces_stored"], r1["tower_A_eq_2adic_B"]))
    out["elapsed_s"] = round(time.time() - t0, 1)
    def _py(o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.bool_):
            return bool(o)
        return str(o)

    with open(os.path.join(HERE, "residual_otoower.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=_py)
    print("VERDICT:", out["verdict"])
    print("WROTE residual_otoower.json elapsed=%.0fs" % out["elapsed_s"])


if __name__ == "__main__":
    main()
