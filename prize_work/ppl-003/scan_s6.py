#!/usr/bin/env python3
# #50 (PPL 003, $250) -- memory-frugal |S|<=6 interior singular-set map (job lab_erdos_261111_44).
# Continues programme card (job 35 = |S|<=4 map, job 38 = local-dim + mass, job 40 = |S|<=5 map).
# Job 40's own next idea (dossier) was exactly this: push |S| to 6 (C(46,6)=9366819 candidates)
# to complete the per-class trend. This job does the |S|=6 step (6-prime class pqrstu / 2pqrst).
#
# WHY |S|=6 extends the frontier (not a re-run):
#   Job 40 enumerated |S|<=5 (1,550,200 candidates -> 30,247 divergent @N=1e7). The |S|=6 class
#   is a NEW prime-set family: 9,366,819 additional candidate c-values (c_S = prod_{q in S}
#   (1-1/q) for 6 distinct primes <=199). If they show divergent windowed-slope (growth>=5) they
#   are NEW interior singular points the published frontier (Erdos 1995: f purely singular,
#   f'=0 a.e. -- a.e. statement only, no interior map) and the lab |S|<=5 map never covered.
#
# MECHANISM (reproduced exactly as the gate): c_S is a divergent-slope interior singular point IFF
#   it is an accumulation point of the value set V: neighbors c_S(1-1/r) -> c_S (r prime, r not in
#   S) pile into [c-w, c] from below, so the two-sided windowed slope s(w)=(F(c+w)-F(c-w))/(2w)
#   DIVERGES as w shrinks (growth=s(w_small)/s(w_large)>=5). A phantom has no pile-up.
#   f'=0 a.e. => a.e. c is a phantom; only the countable value set V (accumulation closure) is
#   singular, a null set (consistent with #50: no x with f' in (0,inf)).
#
# METHOD (reuses job 8/27/35/40 machinery verbatim; identical float64 c-product + WINDOWS +
#   searchsorted, so the |S|<=5 subset is BIT-IDENTICAL to job 40):
#   exact int64 totient sieve (cached) -> r=phi[1:]/n -> SORTED r so F(c)=#n:phi(n)/n<c =
#   searchsorted(r_sorted,c)/N -> for every |S|<=6 value-set c-value a two-sided windowed slope
#   over a 200x shrink (w:0.05 -> 0.00025); growth=s(w_small)/s(w_large); DIVERGE iff growth>=5.
#
# GATE (proves the code is right): the |S|<=5 subset of this scan reproduces job 40's EXACT
#   counts at N=1e7: n_divergent(|S|<=5)=30247, per-|S| {1:8,2:21,3:170,4:2809,5:27239}, |S|=5
#   class pqrst=25278 + 2pqrs=1961; and the |S|<=4 subset = job 35's 3008.
# NEW RESULT (past the frontier): the |S|=6 divergent count @N=1e7 (pqrstu + 2pqrst) + full
#   per-|S| trend (job 40 ended at |S|=5; no published interior map at any |S| or N).
# INDEPENDENT VERIFIER (not the same code twice):
#   (a) second finite-N N=5e6: the |S|=6 count + top-5000 |S|=6 c-set must be stable.
#   (b) brute-force spot-check: top-8 strongest |S|=6 points, recompute the two-sided growth AND
#       the one-sided-below pile-up growth by DIRECT boolean in-window counts on the UNSORTED r
#       array (different code path than searchsorted on r_sorted).
# KILL TEST (analogous to job 40): the N=5e6 run first; if it finds ZERO divergent 6-prime
#   points (growth>=5) the whole |S|=6 class is phantom -> negative -> stop.
# MEMORY-FRUGAL (the |S|<=5 first version OOM'd at 1.26GB): flat numpy arrays + chunked growth
#   scan (only growth kept) + brute-force spot-check after freeing Sarr/cvals. Peak ~490MB.
#
# Usage:  python scan_s6.py 10000000   (N=1e7 primary; run 5000000 first for the kill test)

import sys, os, json, time, math, gc
from math import gcd
from fractions import Fraction
from itertools import combinations
from collections import Counter
import numpy as np

SLUG = "E:/Projects/MillenniumResearch/prize_work/ppl-003"
W_LARGE, W_SMALL, N_WIN = 0.05, 0.00025, 11
WINDOWS = np.geomspace(W_LARGE, W_SMALL, N_WIN)
DIVERGE_THRESH = 5.0
GRID_N, M_FINE, PRIMES_MAX, MAX_S = 4001, 400001, 199, 6
S6_CLASSES = ("pqrstu", "2pqrst")
BLOCK_SCAN = 200000
BLOCK_PM = 1000000
DETAIL_CAP = 5000

# Job-40 ground truth (|S|<=5 subset gate), from scan_s5_N10000000.json:
JOB40_TOTAL_1E7 = 30247
JOB40_ABS_1E7 = {1: 8, 2: 21, 3: 170, 4: 2809, 5: 27239}
JOB40_PQRST_1E7 = 25278
JOB40_2PQRS_1E7 = 1961
# Job-35 ground truth (|S|<=4 subset gate):
JOB35_CLASS_1E7 = {"2only": 1, "2p": 5, "2pq": 7, "2pqr": 125,
                   "p": 7, "pq": 16, "pqr": 163, "pqrs": 2684}
JOB35_TOTAL_1E7 = 3008


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[0] = s[1] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


def phi_sieve_vec(N):
    phi = np.arange(N + 1, dtype=np.int64)
    for p in range(2, N + 1):
        if phi[p] == p:
            for m in range(p, N + 1, p):
                phi[m] = (phi[m] // p) * (p - 1)
    return phi


def class_of(absS, has2):
    L = "pqrstuv"
    if has2:
        rest = absS - 1
        return "2only" if rest == 0 else "2" + L[:rest]
    return L[:absS]


def enumerate_arrays(P, max_S):
    total = int(sum(math.comb(len(P), k) for k in range(1, max_S + 1)))
    cvals = np.empty(total, dtype=np.float64)
    Sarr = np.zeros((total, max_S), dtype=np.int16)
    has2 = np.zeros(total, dtype=bool)
    absS = np.zeros(total, dtype=np.int8)
    f = {q: (q - 1) / q for q in P}
    i = 0
    for k in range(1, max_S + 1):
        for S in combinations(P, k):
            c = 1.0
            for q in S:
                c *= f[q]
            cvals[i] = c
            Sarr[i, :k] = S
            has2[i] = (S[0] == 2)
            absS[i] = k
            i += 1
    return cvals, Sarr, has2, absS


def scan_growth(r_sorted, N, cvals):
    w0, w1 = WINDOWS[0], WINDOWS[-1]
    growth = np.empty(len(cvals))
    for st in range(0, len(cvals), BLOCK_SCAN):
        en = min(st + BLOCK_SCAN, len(cvals))
        C = cvals[st:en]
        s0 = (np.searchsorted(r_sorted, C + w0, 'left') -
              np.searchsorted(r_sorted, C - w0, 'left')) / (2.0 * w0 * N)
        s1 = (np.searchsorted(r_sorted, C + w1, 'left') -
              np.searchsorted(r_sorted, C - w1, 'left')) / (2.0 * w1 * N)
        growth[st:en] = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), np.nan)
    return growth


def growth_direct(r_sorted, N, c):
    s0 = (np.searchsorted(r_sorted, c + W_LARGE, 'left') -
          np.searchsorted(r_sorted, c - W_LARGE, 'left')) / (2.0 * W_LARGE * N)
    s1 = (np.searchsorted(r_sorted, c + W_SMALL, 'left') -
          np.searchsorted(r_sorted, c - W_SMALL, 'left')) / (2.0 * W_SMALL * N)
    return (s1 / s0) if s0 > 0 else 0.0


def point_mass_chunked(phi, a, b):
    N = len(phi) - 1
    cnt = 0
    for st in range(2, N + 1, BLOCK_PM):
        en = min(st + BLOCK_PM, N + 1)
        idx = np.arange(st, en, dtype=np.int64)
        cnt += int(np.sum(b * phi[idx] == a * idx))
    return cnt


def cfrac_str(S, k):
    num = den = 1
    for q in S[:k]:
        num *= (q - 1); den *= q
    g = gcd(num, den)
    return f"{num // g}/{den // g}"


def cstar_pin(r_sorted, N, win=0.002, digits=3):
    cs = np.linspace(0.30, 0.70, 401)
    Fb = np.searchsorted(r_sorted, cs, 'left') / N
    best, best_j = 0.0, 0
    for j in range(1, len(cs) - 1):
        slope = (Fb[j + 1] - Fb[j - 1]) / (2 * (cs[j + 1] - cs[j]))
        if slope > best:
            best, best_j = slope, j
    c0 = float(cs[best_j])
    for _ in range(3):
        lo, hi = max(c0 - win, 0.01), min(c0 + win, 0.99)
        cf = np.linspace(lo, hi, 801)
        Ff = np.searchsorted(r_sorted, cf, 'left') / N
        bestf, bestf_j = 0.0, 0
        for j in range(1, len(cf) - 1):
            slope = (Ff[j + 1] - Ff[j - 1]) / (2 * (cf[j + 1] - cf[j]))
            if slope > bestf:
                bestf, bestf_j = slope, j
        c0 = float(cf[bestf_j]); win *= 0.5
    Fc = float(np.searchsorted(r_sorted, c0, 'left') / N)
    Fl = float(np.searchsorted(r_sorted, c0 - 0.02, 'left') / N)
    Fr = float(np.searchsorted(r_sorted, c0 + 0.02, 'left') / N)
    return {"cstar_pin_%ddp" % digits: round(c0, digits),
            "delta_f_over_pm0p02": round(Fr - Fl, 6),
            "F_cstar": round(Fc, 6),
            "growth_two_side_N": round(growth_direct(r_sorted, N, c0), 4)}


def grid_backstop(r_sorted, N, cvals_sorted):
    grid = np.linspace(0.0, 1.0, GRID_N)
    w0, w1 = WINDOWS[0], WINDOWS[-1]
    s0 = (np.searchsorted(r_sorted, grid + w0, 'left') -
          np.searchsorted(r_sorted, grid - w0, 'left')) / (2.0 * w0 * N)
    s1 = (np.searchsorted(r_sorted, grid + w1, 'left') -
          np.searchsorted(r_sorted, grid - w1, 'left')) / (2.0 * w1 * N)
    gg = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), 0.0)
    div_g = np.flatnonzero(gg >= DIVERGE_THRESH)
    explained = 0
    for gi in div_g:
        g = float(grid[gi])
        if g > 0.985:
            explained += 1
            continue
        j = int(np.searchsorted(cvals_sorted, g))
        nears = []
        if j < len(cvals_sorted):
            nears.append(abs(float(cvals_sorted[j]) - g))
        if j > 0:
            nears.append(abs(float(cvals_sorted[j - 1]) - g))
        if nears and min(nears) < 1e-3:
            explained += 1
    return {"n_div_grid_cells": int(len(div_g)), "n_explained": int(explained),
            "n_unexplained": int(len(div_g) - explained),
            "all_explained": bool(len(div_g) == explained)}


def subset_counts(idx_arr, absS, has2):
    cls = Counter(); abss = Counter()
    for i in idx_arr:
        k = int(absS[i]); h2 = bool(has2[i])
        cls[class_of(k, h2)] += 1; abss[k] += 1
    return cls, abss


def main():
    t0 = time.time()
    if len(sys.argv) < 2:
        print("usage: scan_s6.py N"); sys.exit(1)
    N = int(sys.argv[1])
    cache = f"{SLUG}/phi_cache_{N}.npy"
    if os.path.exists(cache):
        phi = np.load(cache)
    else:
        phi = phi_sieve_vec(N)
        np.save(cache, phi)
    nn = np.arange(1, N + 1, dtype=np.float64)
    r_sorted = np.sort(phi[1:].astype(np.float64) / nn)

    out = {"N": N, "max_S": MAX_S, "primes_max": PRIMES_MAX,
           "w_large": W_LARGE, "w_small": W_SMALL,
           "windows": [round(float(w), 6) for w in WINDOWS],
           "diverge_thresh": DIVERGE_THRESH}

    # --- gates on the r-array (job 40 identical) ---
    gates = {}
    gates["mean_r_6_pi2"] = {"measured": round(float(r_sorted.mean()), 6),
                             "theory": round(6 / math.pi ** 2, 6),
                             "pass": bool(abs(float(r_sorted.mean()) - 6 / math.pi ** 2) < 5e-4)}
    gates["sum_phi_N2_3_pi2"] = {"measured": round(float(phi[1:].sum()) / (N * N), 6),
                                 "theory": round(3 / math.pi ** 2, 6),
                                 "pass": bool(abs(float(phi[1:].sum()) / (N * N) - 3 / math.pi ** 2) < 1e-4)}
    gates["f0_is_0"] = {"pass": bool(r_sorted[0] >= 0.0)}
    gates["f1_is_1"] = {"pass": bool((r_sorted < 1.0).sum() > N - 2)}
    out["gates"] = gates
    out["all_gates_pass"] = all(g["pass"] for g in gates.values())
    out["c12_gate"] = cstar_pin(r_sorted, N)

    job27_gate = {}
    for name, fr, gexp in [("1/3", Fraction(1, 3), 54.6514), ("2/5", Fraction(2, 5), 38.0562),
                           ("3/7", Fraction(3, 7), 18.3267), ("1/2", Fraction(1, 2), 58.4607),
                           ("1/4", Fraction(1, 4), 0.4916)]:
        g = growth_direct(r_sorted, N, float(fr))
        a, b = fr.numerator, fr.denominator
        pm = 0 if b > N else point_mass_chunked(phi, a, b)
        ok = (abs(g - gexp) <= max(1e-4, abs(gexp) * 0.02)) if gexp > 1 else (abs(g - gexp) <= 0.02)
        job27_gate[name] = {"growth": round(float(g), 4), "pm": int(pm),
                            "exp_growth": gexp, "match": bool(ok)}
    out["gate_repro_job27"] = job27_gate
    out["gate_job27_all_match"] = all(v["match"] for v in job27_gate.values())

    # --- enumerate |S|<=6 (identical to job 40 for |S|<=5) ---
    ps = primes_upto(PRIMES_MAX)
    P = [int(q) for q in range(2, PRIMES_MAX + 1) if ps[q]]
    t1 = time.time()
    cvals, Sarr, has2, absS = enumerate_arrays(P, MAX_S)
    out["n_candidates"] = int(len(cvals))
    out["n_candidates_by_absS"] = {str(k): int((absS == k).sum()) for k in range(1, MAX_S + 1)}
    out["t_enum_s"] = round(time.time() - t1, 1)
    print("candidates:", out["n_candidates"],
          "| |S|=6:", int((absS == 6).sum()), "| |S|<=5:", int((absS <= 5).sum()),
          f"[{out['t_enum_s']}s]")

    # --- windowed-slope growth scan ---
    t2 = time.time()
    growth = scan_growth(r_sorted, N, cvals)
    out["t_scan_s"] = round(time.time() - t2, 1)
    div = growth >= DIVERGE_THRESH
    div_idx = np.flatnonzero(div)
    out["n_divergent"] = int(len(div_idx))
    print("scan", out["t_scan_s"], "s | divergent(|S|<=6):", out["n_divergent"])

    # --- GATE: |S|<=5 subset must reproduce job 40 exactly (and |S|<=4 = job 35) @N=1e7 ---
    full_cls, full_abs = subset_counts(div_idx, absS, has2)
    out["divergent_by_absS"] = {str(k): int(full_abs.get(k, 0)) for k in range(1, MAX_S + 1)}
    out["divergent_by_class"] = {c: int(full_cls.get(c, 0))
                                 for c in sorted(set(full_cls.keys()) | set(S6_CLASSES))}

    is1e7 = (N == 10000000)
    s5_idx = div_idx[absS[div_idx] <= 5]
    s5_cls, s5_abs = subset_counts(s5_idx, absS, has2)
    out["gate_s5"] = {
        "s5_total_this": int(len(s5_idx)),
        "s5_by_absS_this": {str(k): int(s5_abs.get(k, 0)) for k in range(1, 6)},
        "s5_by_class_this": {c: int(s5_cls.get(c, 0)) for c in ("pqrst", "2pqrs")},
        "job40_total": JOB40_TOTAL_1E7,
        "job40_by_absS": {str(k): JOB40_ABS_1E7[k] for k in range(1, 6)},
        "job40_pqrst": JOB40_PQRST_1E7, "job40_2pqrs": JOB40_2PQRS_1E7,
        "match_total_job40": (int(len(s5_idx)) == JOB40_TOTAL_1E7) if is1e7 else None,
        "match_absS_job40": (all(s5_abs.get(k, 0) == JOB40_ABS_1E7[k] for k in range(1, 6))) if is1e7 else None,
        "match_s5_class_job40": (s5_cls.get("pqrst", 0) == JOB40_PQRST_1E7 and
                                 s5_cls.get("2pqrs", 0) == JOB40_2PQRS_1E7) if is1e7 else None,
    }
    s4_idx = div_idx[absS[div_idx] <= 4]
    s4_cls, s4_abs = subset_counts(s4_idx, absS, has2)
    out["gate_s4"] = {
        "s4_total_this": int(len(s4_idx)),
        "s4_by_class_this": {c: int(s4_cls.get(c, 0)) for c in JOB35_CLASS_1E7},
        "job35_total": JOB35_TOTAL_1E7,
        "job35_by_class": JOB35_CLASS_1E7,
        "match_total_job35": (int(len(s4_idx)) == JOB35_TOTAL_1E7) if is1e7 else None,
        "match_class_job35": (all(s4_cls.get(c, 0) == v for c, v in JOB35_CLASS_1E7.items())) if is1e7 else None,
    }
    out["gate_pass"] = (bool(out["gate_s5"]["match_total_job40"]) and
                        bool(out["gate_s5"]["match_absS_job40"]) and
                        bool(out["gate_s5"]["match_s5_class_job40"]) and
                        bool(out["gate_s4"]["match_total_job35"])) if is1e7 else None
    if is1e7:
        print("GATE |S|<=5 vs job40:", out["gate_s5"]["match_total_job40"],
              out["gate_s5"]["match_absS_job40"], out["gate_s5"]["match_s5_class_job40"],
              "| |S|<=4 vs job35:", out["gate_s4"]["match_total_job35"],
              "-> gate_pass:", out["gate_pass"])

    # --- NEW RESULT: the |S|=6 divergent class ---
    s6_idx = div_idx[absS[div_idx] == 6]
    s6_cls, _ = subset_counts(s6_idx, absS, has2)
    out["new_s6_total"] = int(len(s6_idx))
    out["new_s6_by_class"] = {c: int(s6_cls.get(c, 0)) for c in S6_CLASSES}
    out["s6_trend_per_absS"] = out["divergent_by_absS"]
    print("NEW |S|=6 divergent:", out["new_s6_total"],
          "by_class:", out["new_s6_by_class"])

    # --- kill test flag (for the N=5e6 run) ---
    out["kill_test_s6_nonzero"] = bool(len(s6_idx) > 0)

    # --- detail: top-5000 divergent by growth ---
    order = div_idx[np.argsort(-growth[div_idx])]
    cap = min(DETAIL_CAP, len(order))
    detail = []
    for rank, ii in enumerate(order[:cap], start=1):
        i = int(ii)
        k = int(absS[i]); h2 = bool(has2[i])
        c = float(cvals[i]); S = [int(x) for x in Sarr[i, :k]]
        cs = cfrac_str(S, k)
        cf = Fraction(cs)
        a, b = cf.numerator, cf.denominator
        pm = 0 if b > N else point_mass_chunked(phi, a, b)
        inw = int(np.searchsorted(r_sorted, c + W_SMALL, 'left') -
                  np.searchsorted(r_sorted, c - W_SMALL, 'left'))
        detail.append({"rank": rank, "c": cs, "class": class_of(k, h2), "absS": k,
                       "S": S, "growth": round(float(growth[i]), 4),
                       "point_mass": int(pm), "in_window_w_small": int(inw)})
    out["top_detail"] = detail

    # --- grid backstop (lightweight: are all divergent grid cells near a |S|<=6 c_S?) ---
    cvals_sorted = np.sort(cvals)
    out["grid_backstop"] = grid_backstop(r_sorted, N, cvals_sorted)
    del cvals_sorted
    print("grid_backstop:", out["grid_backstop"])

    # --- free big arrays before the brute-force spot check ---
    del Sarr, cvals, absS, has2, growth
    gc.collect()

    # --- INDEPENDENT VERIFIER (b): brute-force spot check on top-8 |S|=6 (direct boolean,
    #     UNSORTED r -- a different code path than searchsorted on r_sorted) ---
    s6_detail = [d for d in detail if d["absS"] == 6][:8]
    if s6_detail:
        r_un = phi[1:].astype(np.float64) / nn   # UNSORTED r, transient
        w0, w1 = WINDOWS[0], WINDOWS[-1]
        spot = []
        for d in s6_detail:
            c = float(Fraction(d["c"]))
            cnt0 = int(((r_un >= c - w0) & (r_un < c + w0)).sum()); s0 = cnt0 / (2 * w0 * N)
            cnt1 = int(((r_un >= c - w1) & (r_un < c + w1)).sum()); s1 = cnt1 / (2 * w1 * N)
            g_twoside = (s1 / s0) if s0 > 0 else 0.0
            b0 = int(((r_un >= c - w0) & (r_un < c)).sum()); sb0 = b0 / (w0 * N)
            b1 = int(((r_un >= c - w1) & (r_un < c)).sum()); sb1 = b1 / (w1 * N)
            g_below = (sb1 / sb0) if sb0 > 0 else 0.0
            a0 = int(((r_un > c) & (r_un <= c + w0)).sum()); sa0 = a0 / (w0 * N)
            a1 = int(((r_un > c) & (r_un <= c + w1)).sum()); sa1 = a1 / (w1 * N)
            g_above = (sa1 / sa0) if sa0 > 0 else 0.0
            g_ss = float(d["growth"])
            spot.append({"c": d["c"], "growth_searchsorted": g_ss,
                         "growth_brute_twoside": round(float(g_twoside), 4),
                         "growth_brute_below": round(float(g_below), 4),
                         "growth_brute_above": round(float(g_above), 4),
                         "side": "below" if g_below >= g_above else "above",
                         "confirm_twoside": bool(g_twoside >= DIVERGE_THRESH and abs(g_twoside - g_ss) <= 0.1 * g_ss + 1.0),
                         "confirm_below": bool(g_below >= DIVERGE_THRESH),
                         "confirm_above": bool(g_above >= DIVERGE_THRESH)})
        out["brute_force_spot_check_top8_s6"] = spot
        out["spot_confirm_all_twoside"] = all(x["confirm_twoside"] for x in spot)
        out["spot_sides"] = Counter(x["side"] for x in spot)
        out["spot_all_above"] = bool(out["spot_sides"].get("above", 0) == len(spot))
        out["spot_all_below"] = bool(out["spot_sides"].get("below", 0) == len(spot))
        del r_un
        gc.collect()
        print("spot_check top8 |S|=6 two-sided confirm:", out["spot_confirm_all_twoside"],
              "sides:", dict(out["spot_sides"]), "all_above:", out["spot_all_above"])

    # --- INDEPENDENT VERIFIER (a): stability vs the second N (read the other run's JSON) ---
    other = 10000000 if N == 5000000 else 5000000
    oth = f"{SLUG}/scan_s6_N{other}.json"
    if os.path.exists(oth):
        with open(oth) as fh:
            o = json.load(fh)
        o6 = set(d["c"] for d in o.get("top_detail", []) if d.get("absS") == 6)
        m6 = set(d["c"] for d in detail if d["absS"] == 6)
        inter6 = len(o6 & m6)
        out["stability_vs_N"] = {
            "other_N": other,
            "s6_total_this": int(out["new_s6_total"]),
            "s6_total_other": int(o.get("new_s6_total", -1)),
            "s6_topK_this": len(m6),
            "s6_topK_inter": inter6,
            "s6_topK_stable_frac": round(inter6 / len(m6), 4) if m6 else 1.0,
            "s6_count_ratio_this_over_other":
                (round(out["new_s6_total"] / o["new_s6_total"], 4)
                 if o.get("new_s6_total") else None),
        }
        print("stability vs N", other, ":", out["stability_vs_N"])
    else:
        out["stability_vs_N"] = {"other_N": other, "note": "other run JSON not present yet"}

    out["t_total_s"] = round(time.time() - t0, 1)
    outp = f"{SLUG}/scan_s6_N{N}.json"
    with open(outp, "w") as fh:
        json.dump(out, fh, indent=2, default=_js)
    print("wrote", outp, "| total", out["t_total_s"], "s")


def _js(o):
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


if __name__ == "__main__":
    main()
