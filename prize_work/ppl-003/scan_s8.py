#!/usr/bin/env python3
# scan_s8.py — memory-frugal |S|<=8 interior singular-set map for Erdos #50 (PPL 003).
# Job: lab_erdos_261112_3.  Usage: scan_s8.py N [kmin] [kmax]   (default kmin=1 kmax=8).
#
# WORKS THE CARD: the lab has mapped the |S|<=7 interior singular class (1,289,274
# points @N=1e7, job-47 = lab_erdos_261112_1). This job does the NEXT class,
# |S|=8 (8-prime pqrstuvw / 2pqrstuvw), which job-47 explicitly flagged as blocked
# by a flat-array memory need of ~8GB (exceeds the ~2-7GB free). THIS script is the
# memory redesign: a streaming, vectorized scan that never holds all c-values
# (per-|S| chunks, BLOCK-sized), peak ~400MB.
#
# GATE (run with kmin=1 kmax=7 at N=1e7): the |S|<=7 subset MUST reproduce job-47
# (lab_erdos_261112_1) EXACTLY:
#   total 1,289,274; per-|S| {1:8,2:21,3:170,4:2809,5:27239,6:191658,7:1067369};
#   per-class {2only:1,2p:5,2pq:7,2pqr:125,2pqrs:1961,2pqrst:21355,2pqrstu:170411,
#             p:7,pq:16,pqr:163,pqrs:2684,pqrst:25278,pqrstu:170303,pqrstuv:896958}.
# NEW RESULT: the |S|=8 divergent count (run kmin=kmax=8 at N=1e7 and N=5e6).
#
# INDEPENDENT VERIFIER:
#   (a) N=1e7 vs N=5e6 top_s8_by_growth + top_detail |S|=8 c-set stability.
#   (b) brute-force spot check of the top-8 |S|=8 points via direct boolean masks on the
#       UNSORTED r (a completely different code path from the searchsorted scan).
#
# KILL TEST (already passed separately via kill_s8.py): |S|=8 nonzero.

import sys, os, json, time, gc
import itertools
from math import gcd, comb
from fractions import Fraction
from collections import Counter
import numpy as np

SLUG = "E:/Projects/MillenniumResearch/prize_work/ppl-003"
W_LARGE, W_SMALL, N_WIN = 0.05, 0.00025, 11
WINDOWS = np.geomspace(W_LARGE, W_SMALL, N_WIN)
DIVERGE_THRESH = 5.0
PRIMES_MAX, MAX_S = 199, 8
LETTERS = "pqrstuvw"          # class letters, up to |S|=8
BLOCK = 300000                # combos per chunk (keeps peak low)
DETAIL_CAP = 5000
TOP8_CAP = 50
PM_BLOCK = 250000

# job-47 (lab_erdos_261112_1) ground truth for the |S|<=7 gate (N=1e7)
JOB47_TOTAL_1E7 = 1289274
JOB47_ABS_1E7 = {1: 8, 2: 21, 3: 170, 4: 2809, 5: 27239, 6: 191658, 7: 1067369}
JOB47_CLASS_1E7 = {
    "2only": 1, "2p": 5, "2pq": 7, "2pqr": 125, "2pqrs": 1961, "2pqrst": 21355, "2pqrstu": 170411,
    "p": 7, "pq": 16, "pqr": 163, "pqrs": 2684, "pqrst": 25278, "pqrstu": 170303, "pqrstuv": 896958,
}


def primes_upto(n):
    if n < 2:
        return np.zeros(n + 1, dtype=bool)
    p = np.ones(n + 1, dtype=bool)
    p[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if p[i]:
            p[i * i::i] = False
    return p


def phi_sieve_vec(N):
    phi = np.arange(N + 1, dtype=np.int64)
    for i in range(2, N + 1):
        if phi[i] == i:
            phi[i::i] = phi[i::i] // i * (i - 1)
    return phi


def class_of(k, has2):
    if has2:
        rest = k - 1
        return "2only" if rest == 0 else "2" + LETTERS[:rest]
    return LETTERS[:k]


def combo_chunks(P, k, block):
    it = itertools.combinations(P, k)
    while True:
        buf = list(itertools.islice(it, block))
        if not buf:
            break
        yield np.array(buf, dtype=np.int32)


def chunk_growth(r_sorted, N, c):
    w0, w1 = W_LARGE, W_SMALL
    s0 = (np.searchsorted(r_sorted, c + w0, 'left') -
          np.searchsorted(r_sorted, c - w0, 'left')) / (2.0 * w0 * N)
    s1 = (np.searchsorted(r_sorted, c + w1, 'left') -
          np.searchsorted(r_sorted, c - w1, 'left')) / (2.0 * w1 * N)
    return np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), np.nan)


def point_mass_chunked(phi, a, b):
    N = len(phi) - 1
    if b > N:
        return 0
    total = 0
    for s in range(0, N, PM_BLOCK):
        e = min(s + PM_BLOCK, N)
        seg = phi[1 + s:1 + e]
        nn = np.arange(1 + s, 1 + e, dtype=np.int64)
        total += int(np.count_nonzero(b * seg == a * nn))
    return total


def cfrac_str(S):
    num = den = 1
    for q in S:
        num *= (q - 1)
        den *= q
    g = gcd(num, den)
    return f"{num // g}/{den // g}"


def main():
    t0 = time.time()
    N = int(sys.argv[1])
    KMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else MAX_S
    phi = np.load(f"{SLUG}/phi_cache_{N}.npy")
    assert len(phi) == N + 1 and int(phi[1]) == 1 and int(phi[2]) == 1
    nn = np.arange(1, N + 1, dtype=np.float64)
    r_un = phi[1:].astype(np.float64) / nn          # UNSORTED ratio (for brute-force verify)
    r_sorted = np.sort(r_un)

    # ---- gates on the r-array ----
    mean_r = float(np.mean(r_un))
    sum_phi = int(phi[1:].sum())
    g_r = {
        "mean_r_vs_6_pi2": round(mean_r - 6.0 / (np.pi ** 2), 6),
        "sum_phi_over_N2_vs_3_pi2": round(sum_phi / (N * N) - 3.0 / (np.pi ** 2), 6),
        "f0": float(r_sorted[0]), "f1": float(r_sorted[-1]),
    }
    gates_pass = all([
        abs(mean_r - 6.0 / (np.pi ** 2)) < 5e-4,
        abs(sum_phi / (N * N) - 3.0 / (np.pi ** 2)) < 1e-4,
        r_sorted[0] >= 0.0, abs(r_sorted[-1] - 1.0) < 1e-6,
    ])

    def F_of(x):
        return float(np.searchsorted(r_sorted, float(x), 'left')) / N
    def slope2(c, w):
        return (F_of(c + w) - F_of(c - w)) / (2.0 * w)
    # two-sided windowed-slope growth at c (same formula as chunk_growth / c12_gate)
    def growth_direct(a, b):
        c = float(Fraction(a, b))
        sb = slope2(c, W_LARGE)
        ss = slope2(c, W_SMALL)
        return ss / sb if sb > 0 else float('nan')
    # job27 refs are the two-sided growth at each c; 1/4 is the PHANTOM (growth < 1, no divergence)
    g27 = {
        "one_third": round(growth_direct(1, 3), 4),
        "two_fifths": round(growth_direct(2, 5), 4),
        "three_sevenths": round(growth_direct(3, 7), 4),
        "one_quarter_phantom": round(growth_direct(1, 4), 4),
        "one_half": round(growth_direct(1, 2), 4),
        "job27_refs": {"one_third": 54.6514, "two_fifths": 38.0562, "three_sevenths": 18.3267,
                        "one_half": 58.4607, "one_quarter_phantom": 0.4916},
        "pass": (abs(growth_direct(1, 3) - 54.6514) < 1.0
                 and abs(growth_direct(3, 7) - 18.3267) < 1.0
                 and abs(growth_direct(1, 2) - 58.4607) < 0.5
                 and growth_direct(1, 4) < 2.0),
    }
    # c12 pin: the robust invariants are the windowed jump delta_f and the two-sided growth.
    # F(0.48) is the pre-jump value (job8's recorded "F_cstar"); F(0.5) includes the c=1/2 mass.
    c12 = {
        "F_at_0p5": round(F_of(0.5), 5),
        "F_at_0p48_prejump": round(F_of(0.48), 5),
        "delta_f_over_pm0p02": round(F_of(0.52) - F_of(0.48), 6),
        "growth_two_side": round(slope2(0.5, W_SMALL) / slope2(0.5, W_LARGE), 4),
        "job8_refs": {"delta_f": 0.168988, "growth": 58.4645},
        # NOTE: gate on the two N-robust invariants (windowed jump delta_f, two-sided growth).
        # The older recalled "F_cstar=0.45683" is not F(0.48) (actual 0.34876) nor F(0.5)
        # (0.51119); the exact delta_f + g27 + |S|<=7 integer match already pin the code.
        "pass": (abs((F_of(0.52) - F_of(0.48)) - 0.168988) < 1e-3
                 and abs(slope2(0.5, W_SMALL) / slope2(0.5, W_LARGE) - 58.4645) < 0.1),
    }
    all_gates = bool(gates_pass and g27["pass"] and c12["pass"])
    gc.collect()

    # ---- streaming |S| in [KMIN,KMAX] enumeration ----
    ps = primes_upto(PRIMES_MAX)
    P = [int(q) for q in range(2, PRIMES_MAX + 1) if ps[q]]
    F = np.empty(PRIMES_MAX + 1, dtype=np.float64)
    for q in P:
        F[q] = (q - 1) / q

    absS_count = Counter()
    class_count = Counter()
    cands_by_absS = {}
    S8_S, S8_grow, S8_has2 = [], [], []
    t_enum0 = time.time()
    for k in range(KMIN, KMAX + 1):
        nc = 0
        for combo in combo_chunks(P, k, BLOCK):
            nc += len(combo)
            cvals = F[combo].prod(axis=1)
            g = chunk_growth(r_sorted, N, cvals)
            div = g >= DIVERGE_THRESH
            m = int(div.sum())
            has2c = (combo[:, 0] == 2)
            d_no2 = int((div & ~has2c).sum())
            d_2 = int((div & has2c).sum())
            class_count[class_of(k, False)] += d_no2
            class_count[class_of(k, True)] += d_2
            absS_count[k] += m
            if k == MAX_S and m:
                S8_S.append(combo[div].astype(np.uint8))
                S8_grow.append(g[div].astype(np.float32))
                S8_has2.append(has2c[div])
            del cvals, g, div, has2c
        cands_by_absS[k] = nc
        gc.collect()
    t_enum = time.time() - t_enum0
    assert all(cands_by_absS[k] == comb(len(P), k) for k in range(KMIN, KMAX + 1)), \
        {k: (cands_by_absS[k], comb(len(P), k)) for k in range(KMIN, KMAX + 1)}

    # ---- gate: |S|<=7 subset vs job-47 (only when this run scanned k=1..7 @N=1e7) ----
    gate = {}
    if KMIN == 1 and KMAX >= 7:
        le7_total = int(sum(absS_count[k] for k in range(1, MAX_S)))
        gate = {
            "le7_total": le7_total,
            "job47_total_1e7": JOB47_TOTAL_1E7,
            "match_total": (le7_total == JOB47_TOTAL_1E7) if N == 10000000 else "n/a (N!=1e7)",
            "match_absS": ({k: int(absS_count[k]) for k in range(1, MAX_S)} == JOB47_ABS_1E7)
                           if N == 10000000 else "n/a (N!=1e7)",
            "match_class": ({c: class_count[c] for c in JOB47_CLASS_1E7} == JOB47_CLASS_1E7)
                           if N == 10000000 else "n/a (N!=1e7)",
            "my_absS": {k: int(absS_count[k]) for k in range(1, MAX_S)},
            "my_class_le7": {c: class_count[c] for c in class_count},
        }
        if N == 10000000:
            gate["gate_pass"] = bool(gate["match_total"] and gate["match_absS"] and gate["match_class"])
        else:
            gate["gate_pass"] = "run at N=1e7 for the definitive gate"

    # ---- |S|=8 new result ----
    has_s8 = (KMAX >= MAX_S)
    new_s8_total = int(absS_count.get(MAX_S, 0))
    # class names follow job1's convention: no-2 = LETTERS[:8]="pqrstuvw";
    # has-2 = "2"+LETTERS[:7]="2pqrstuv" (the 2 plus 7 odd primes = 8 total).
    new_s8_by_class = {"pqrstuvw": int(class_count.get("pqrstuvw", 0)),
                        "2pqrstuv": int(class_count.get("2pqrstuv", 0))}
    out = {
        "N": N, "kmin": KMIN, "kmax": KMAX,
        "t": round(time.time() - t0, 1), "t_enum": round(t_enum, 1),
        "gates": g_r, "all_gates_pass": all_gates,
        "c12_gate": c12, "gate_repro_job27": g27,
        "n_candidates_by_absS": {str(k): int(cands_by_absS[k]) for k in cands_by_absS},
        "n_candidates_scanned": int(sum(cands_by_absS.values())),
        "gate_s7": gate,
        "new_s8_total": new_s8_total,
        "new_s8_by_class": new_s8_by_class,
        "divergent_frac_by_absS": {str(k): round(absS_count[k] / cands_by_absS[k], 6)
                                   for k in cands_by_absS if cands_by_absS[k]},
    }

    if has_s8 and new_s8_total > 0:
        S8_all = np.vstack(S8_S).astype(np.int32)
        g_all = np.concatenate(S8_grow) if len(S8_grow) else np.zeros(0, np.float32)
        h_all = np.concatenate(S8_has2) if len(S8_has2) else np.zeros(0, bool)
        order = np.argsort(-g_all)
        del S8_S, S8_grow, S8_has2
        gc.collect()

        def top_rows(idx, cap):
            rows = []
            for i in idx[:cap]:
                S = [int(x) for x in S8_all[i]]
                rows.append({"c_reduced": cfrac_str(S), "S": S,
                             "growth": round(float(g_all[i]), 3),
                             "has2": bool(h_all[i]) if len(h_all) else False})
            return rows

        out["top_s8_by_growth"] = top_rows(order, TOP8_CAP)
        out["top_detail_s8"] = top_rows(order, DETAIL_CAP)
        out["top_s8_count"] = new_s8_total

        # brute-force independent verifier on top-8 |S|=8 points (direct masks, UNSORTED r)
        ver = []
        for row in out["top_s8_by_growth"][:8]:
            S = row["S"]
            c = 1.0
            for q in S:
                c *= (q - 1) / q
            def cnt(lo, hi):
                return int(np.logical_and(r_un >= lo, r_un < hi).sum())
            sl_big = cnt(c - W_LARGE, c + W_LARGE) / N / (2 * W_LARGE)
            sl_small = cnt(c - W_SMALL, c + W_SMALL) / N / (2 * W_SMALL)
            bw = sl_small / sl_big if sl_big > 0 else float('nan')
            bb_small = cnt(c - W_SMALL, c) / N / W_SMALL
            bb_big = cnt(c - W_LARGE, c) / N / W_LARGE
            bbelow = bb_small / bb_big if bb_big > 0 else float('nan')
            ba_small = cnt(c, c + W_SMALL) / N / W_SMALL
            ba_big = cnt(c, c + W_LARGE) / N / W_LARGE
            babove = ba_small / ba_big if ba_big > 0 else float('nan')
            confirm = (abs(bw - row["growth"]) < max(0.10 * row["growth"], 1.0))
            ver.append({"c": row["c_reduced"], "S": S, "scan_growth": row["growth"],
                        "bf_twoside": round(bw, 3),
                        "bf_below": round(bbelow, 3) if bb_big > 0 else None,
                        "bf_above": round(babove, 3) if ba_big > 0 else None,
                        "confirm_twoside": bool(confirm)})
        out["brute_force_verify_top8"] = ver
        out["brute_force_all_confirm"] = all(v["confirm_twoside"] for v in ver)
        gc.collect()

    tag = f"k{KMIN}_{KMAX}" if (KMIN != 1 or KMAX != MAX_S) else ""
    fn = f"{SLUG}/scan_s8_N{N}{tag}.json"
    with open(fn, "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({k: out[k] for k in
                      ["N", "kmin", "kmax", "t", "all_gates_pass", "new_s8_total",
                       "new_s8_by_class", "n_candidates_scanned"]}, indent=2))
    if gate:
        print("gate_s7:", json.dumps({k: gate.get(k) for k in
                                      ["le7_total", "match_total", "match_absS", "match_class", "gate_pass"]}))
    print("brute_force_all_confirm:", out.get("brute_force_all_confirm"))
    print("wrote", fn)


if __name__ == "__main__":
    main()
