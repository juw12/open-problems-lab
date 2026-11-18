"""#50 (PPL 003, $250) -- MEMORY-FRUGAL |S|<=5 interior singular-set map (job lab_erdos_261111_40).

The first version (scan_singular_set_s5.py) OOM-killed (~1.26GB peak on a ~3GB-free shared PC):
it built a Python list of 1.55M heavy Fraction tuples + a 1.55M-key dict. This version:
  * enumerates |S|<=5 as FLAT numpy arrays (cvals float64, Sarr int16, has2 bool, absS int8)
    preallocated -- no Fraction objects, no big lists.
  * CHUNKS the two-sided windowed-slope growth scan (block=200k) over r_sorted via searchsorted.
  * CHUNKS the exact point-mass count (block=1M) so its transient is ~26MB regardless of N.
  * drops the unsorted r array after sorting (only r_sorted kept); uses searchsorted for F(c).
  * loads the exact phi cache as int64 only where needed; frees intermediates.
Everything else (windows, DIVERGE_THRESH, analytic + job-27 + c=1/2 + |S|<=4-repro job-35 gates,
grid backstop, stability) reuses the job 8/27/35 machinery verbatim so the code is bit-identical.

RESULT: the NEW |S|=5 class (pqrst, 2pqrs) divergent points + the full |S|<=5 map.
GATE  (code right only if all pass): analytic + job-27 growth/mass + c=1/2 pin + the |S|<=4
      subset reproducing job-35's 3008 per-class EXACTLY at N=1e7.
VERIFIER: |S|=5 divergent c-set stable across N (5e6 vs 1e7), c-locations exact rationals.
Run: python scan_s5.py <N>   (N=5000000 kill-test, N=10000000 primary)
Writes: scan_s5_N<N>.json
"""
import sys, os, json, time, gc
import numpy as np
from math import comb, gcd
from itertools import combinations
from collections import Counter

SLUG = "E:/Projects/MillenniumResearch/prize_work/ppl-003"
W_LARGE, W_SMALL, N_WIN = 0.05, 0.00025, 11
WINDOWS = np.geomspace(W_LARGE, W_SMALL, N_WIN)
DIVERGE_THRESH = 5.0
GRID_N = 4001
M_FINE = 400001
PRIMES_MAX = 199
MAX_S = 5
S5_CLASSES = ("pqrst", "2pqrs")
BLOCK_SCAN = 200000
BLOCK_PM = 1000000

JOB35_CLASS_1E7 = {"2only": 1, "2p": 5, "2pq": 7, "2pqr": 125,
                   "p": 7, "pq": 16, "pqr": 163, "pqrs": 2684}
JOB35_TOTAL_1E7 = 3008
JOB35_TOTAL_5E6 = 3005


def primes_upto(n):
    if n < 2: return []
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s).tolist()


def phi_sieve_vec(N):
    phi = np.arange(N + 1, dtype=np.int64)
    for p in primes_upto(N):
        idx = np.arange(p, N + 1, p)
        phi[idx] = (phi[idx] // p) * (p - 1)
    return phi


def class_of(absS, has2):
    L = "pqrst"
    if has2:
        rest = absS - 1
        return "2only" if rest == 0 else "2" + L[:rest]
    return L[:absS]


def enumerate_arrays(P, max_S):
    """flat, memory-frugal: cvals(float64), Sarr(int16, n x max_S, 0-padded), has2(bool), absS(int8)."""
    total = sum(comb(len(P), k) for k in range(1, max_S + 1))
    cvals = np.empty(total, dtype=np.float64)
    Sarr = np.zeros((total, max_S), dtype=np.int16)
    has2 = np.zeros(total, dtype=bool)
    absS = np.zeros(total, dtype=np.int8)
    f = {q: (q - 1) / q for q in P}
    i = 0
    for k in range(1, max_S + 1):
        for S in combinations(P, k):
            c = 1.0
            for q in S: c *= f[q]
            cvals[i] = c
            Sarr[i, :k] = S
            has2[i] = (2 in S)
            absS[i] = k
            i += 1
    print("enumerated %d candidates (|S|<=%d)" % (i, max_S))
    return cvals, Sarr, has2, absS


def scan_growth(r_sorted, N, cvals, windows):
    """two-sided slope growth s(w_small)/s(w_large), CHUNKED to keep transient memory flat."""
    w0, w1 = windows[0], windows[-1]
    s0 = np.empty(len(cvals)); s1 = np.empty(len(cvals))
    for st in range(0, len(cvals), BLOCK_SCAN):
        en = min(st + BLOCK_SCAN, len(cvals))
        C = cvals[st:en]
        s0[st:en] = (np.searchsorted(r_sorted, C + w0, 'left') -
                     np.searchsorted(r_sorted, C - w0, 'left')) / (2.0 * w0 * N)
        s1[st:en] = (np.searchsorted(r_sorted, C + w1, 'left') -
                     np.searchsorted(r_sorted, C - w1, 'left')) / (2.0 * w1 * N)
    with np.errstate(divide='ignore', invalid='ignore'):
        growth = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), np.nan)
    return growth, s0, s1


def growth_direct(r_sorted, N, c, windows):
    def s_at(w):
        return (np.searchsorted(r_sorted, c + w, 'left') -
                np.searchsorted(r_sorted, c - w, 'left')) / (2.0 * w * N)
    s = [s_at(w) for w in windows]
    return s[-1] / s[0] if s[0] > 0 else float('nan')


def point_mass_chunked(phi, a, b):
    """# n in [1,N] with b*phi[n] == a*n, CHUNKED (transient ~ a*block)."""
    total = 0
    for st in range(1, len(phi), BLOCK_PM):
        en = min(st + BLOCK_PM, len(phi))
        nch = np.arange(st, en, dtype=np.int64)
        total += int((b * phi[st:en] == a * nch).sum())
    return total


def fine_cdf(r_sorted):
    edges = np.linspace(0.0, 1.0, M_FINE + 1)
    hist, _ = np.histogram(r_sorted, bins=edges)
    F = np.zeros(len(edges)); F[1:] = np.cumsum(hist) / float(len(r_sorted))
    return edges, F


def cstar_pin(edges, F, window):
    dc = float(edges[1] - edges[0]); h = max(1, int(round(window / (2.0 * dc))))
    s = (F[2 * h:] - F[:-2 * h]) / (2.0 * h * dc)
    return float(edges[h + 1 + int(np.argmax(s))])


def cfrac_str(S, k):
    a, b = 1, 1
    for q in S[:k]:
        a *= (q - 1); b *= q
    g = gcd(a, b); a //= g; b //= g
    return "%d/%d" % (a, b)


def main():
    t0 = time.time()
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000
    out = {"N": N, "max_S": MAX_S, "w_large": W_LARGE, "w_small": W_SMALL,
           "diverge_thresh": DIVERGE_THRESH, "elapsed_s": 0.0}
    pfile = SLUG + "/phi_cache_%s.npy" % N
    phi = np.load(pfile) if os.path.exists(pfile) else None
    if phi is None:
        phi = phi_sieve_vec(N); np.save(pfile, phi); print("built+cached phi")
    else:
        print("loaded cached phi (N=%d)" % N)

    pf = phi[1:].astype(np.float64)
    nf = np.arange(1, N + 1, dtype=np.float64)
    r_sorted = np.sort(pf / nf)
    del pf, nf; gc.collect()
    mean_r = float(r_sorted.mean()); sum_phi = int(phi[1:].sum())
    gates = {
        "mean_r_6_pi2":  abs(mean_r - 6.0 / np.pi ** 2) < 1e-4,
        "sum_phi_3_pi2": abs(sum_phi / N ** 2 - 3.0 / np.pi ** 2) < 1e-4,
        "f0_is_0":       bool((r_sorted < 0.0).sum() == 0),
        "f1_is_1":       bool((r_sorted < 1.0).sum() > N - 2),
    }
    out["gates"] = gates; out["all_gates_pass"] = bool(all(gates.values()))

    # c=1/2 gate (job 8)
    edges, F = fine_cdf(r_sorted); del edges, F; gc.collect()
    cstar02 = cstar_pin(*fine_cdf(r_sorted), 0.002)
    Fb = lambda c: np.searchsorted(r_sorted, c, 'left') / N
    delta_f = Fb(cstar02 + 0.02) - Fb(cstar02 - 0.02)
    growth05 = growth_direct(r_sorted, N, 0.5, WINDOWS)
    out["c12_gate"] = {
        "cstar_pin_3dp": round(cstar02, 3),
        "delta_f_over_pm0p02": round(delta_f, 6),
        "growth_two_side_N": round(growth05, 4),
        "gate_cstar_0p500": abs(cstar02 - 0.5) < 0.005,
        "gate_delta_f_0p169": abs(delta_f - 0.169) < 0.005,
        "gate_growth05": abs(growth05 - 58.4607) < 0.5,
    }

    # job-27 two-sided growth + exact point mass (reduced a/b for each gate value)
    job27 = [("1/3", 1, 3, 54.6514, 152), ("2/5", 2, 5, 38.0562, 100),
             ("3/7", 3, 7, 18.3267, 81), ("1/2", 1, 2, 58.4607, 23),
             ("1/4", 1, 4, 0.4916, 0)]
    gr = {}
    for lab, a, b, g27, m27 in job27:
        c = a / b
        g = growth_direct(r_sorted, N, c, WINDOWS)
        m = point_mass_chunked(phi, a, b)
        gr[lab] = {"c": lab, "growth": round(g, 4), "point_mass": m,
                   "match_growth_job27": abs(g - g27) < 0.1,
                   "match_mass_job27": (m == m27) if N == 10_000_000 else None}
    out["gate_repro_job27"] = gr
    out["gate_repro_all_match"] = bool(all(v["match_growth_job27"] for v in gr.values()))
    print("job27 gate all-match:", out["gate_repro_all_match"],
          {k: (v["growth"], v["point_mass"]) for k, v in gr.items()})

    # enumerate |S|<=5 candidates (frugal) + controls
    P = primes_upto(PRIMES_MAX)
    cvals, Sarr, has2, absS = enumerate_arrays(P, MAX_S)
    control = {"0.25_phantom": 0.25, "0.30_gen": 0.30, "0.45_gen": 0.45,
               "0.60_gen": 0.60, "0.70_gen": 0.70, "0.85_gen": 0.85}
    all_c = np.concatenate([cvals, np.array(list(control.values()))])
    t_scan = time.time()
    growth_all, s0_all, s1_all = scan_growth(r_sorted, N, all_c, WINDOWS)
    growth_cand = growth_all[:len(cvals)]
    growth_ctrl = {k: (None if not np.isfinite(g) else round(float(g), 4))
                   for k, g in zip(control, growth_all[len(cvals):])}
    print("scan growth %.1fs over %d cands" % (time.time() - t_scan, len(cvals)))
    del s0_all, s1_all, all_c; gc.collect()
    diverge_cand = growth_cand >= DIVERGE_THRESH

    # THE KEY GATE: |S|<=4 subset must reproduce job-35
    s4 = diverge_cand & (absS <= 4)
    s4_by_class = Counter(class_of(int(absS[i]), bool(has2[i])) for i in np.flatnonzero(s4))
    gate_s4 = {
        "s4_total": int(s4.sum()),
        "s4_by_class": dict(s4_by_class),
        "job35_total": JOB35_TOTAL_1E7, "job35_by_class": JOB35_CLASS_1E7,
        "job35_total_5e6": JOB35_TOTAL_5E6,
        "match_total_job35": (int(s4.sum()) == JOB35_TOTAL_1E7) if N == 10_000_000
                             else (int(s4.sum()) == JOB35_TOTAL_5E6) if N == 5_000_000 else None,
        "match_byclass_job35": (dict(s4_by_class) == JOB35_CLASS_1E7) if N == 10_000_000 else None,
    }
    out["gate_s4_repro_job35"] = gate_s4
    print("GATE |S|<=4: total=%d  byclass=%s  match_total=%s"
          % (gate_s4["s4_total"], dict(s4_by_class), gate_s4["match_total_job35"]))

    # NEW RESULT: |S|=5 + full |S|<=5 divergent set
    div_idx = np.flatnonzero(diverge_cand)
    by_class = Counter(class_of(int(absS[i]), bool(has2[i])) for i in div_idx)
    by_absS = Counter(int(absS[i]) for i in div_idx)
    out["n_candidates"] = int(len(cvals))
    out["n_divergent"] = int(diverge_cand.sum())
    out["divergent_by_class"] = {k: int(v) for k, v in sorted(by_class.items())}
    out["divergent_by_absS"] = {str(k): int(v) for k, v in sorted(by_absS.items())}
    out["new_s5_by_class"] = {k: int(v) for k, v in by_class.items() if k in S5_CLASSES}
    out["new_s5_total"] = int(by_absS.get(5, 0))
    out["controls_growth"] = growth_ctrl
    out["controls_all_phantom"] = bool(all((v is None) or (v < DIVERGE_THRESH) for v in growth_ctrl.values()))

    # per-point detail for top-2500 divergent by growth (exact mass/in-window)
    order = div_idx[np.argsort(-growth_cand[div_idx])][:2500]
    detail = []
    for i in order:
        S_i = [int(x) for x in Sarr[i][:absS[i]]]
        a_full, b_full = 1, 1
        for q in S_i:
            a_full *= (q - 1); b_full *= q
        g = gcd(a_full, b_full)
        a, b = a_full // g, b_full // g
        cf = float(cvals[i])
        # smallest n whose prime-support is exactly S is prod(S); if > N, point-mass is 0
        pm = 0 if b_full > N else int(point_mass_chunked(phi, a, b))
        rec = {"c": "%d/%d" % (a, b), "class": class_of(int(absS[i]), bool(has2[i])),
               "absS": int(absS[i]), "S": S_i, "growth": round(float(growth_cand[i]), 3),
               "point_mass": pm,
               "in_window": int(np.searchsorted(r_sorted, cf + W_SMALL, 'left') -
                                 np.searchsorted(r_sorted, cf - W_SMALL, 'left'))}
        detail.append(rec)
    out["divergent_detail"] = detail
    s5_detail = [d for d in detail if d["absS"] == 5]
    out["s5_top20_by_growth"] = s5_detail[:20]
    out["s5_count_in_detail"] = len(s5_detail)
    print("NEW |S|<=5 @N=%d: total divergent=%d by_absS=%s new|S|=5=%d byclass(s5)=%s"
          % (N, out["n_divergent"], out["divergent_by_absS"], out["new_s5_total"], out["new_s5_by_class"]))

    # STABILITY vs the other N (card verifier)
    other_N = 5_000_000 if N == 10_000_000 else 10_000_000
    ofile = SLUG + "/scan_s5_N%d.json" % other_N
    if os.path.exists(ofile):
        od = json.load(open(ofile))
        o_s5 = {d["c"] for d in od.get("divergent_detail", []) if d.get("absS") == 5}
        o_all = {d["c"] for d in od.get("divergent_detail", [])}
        my_s5 = {d["c"] for d in detail if d["absS"] == 5}
        my_all = {d["c"] for d in detail}
        stab = {
            "other_N": other_N,
            "s5_intersect": len(my_s5 & o_s5), "s5_this": len(my_s5), "s5_other": len(o_s5),
            "s5_stable_frac": round(len(my_s5 & o_s5) / max(1, len(my_s5)), 4),
            "all_intersect": len(my_all & o_all), "all_this": len(my_all), "all_other": len(o_all),
            "all_stable_frac": round(len(my_all & o_all) / max(1, len(my_all)), 4),
        }
        out["stability_vs_other_N"] = stab
        print("STABILITY vs N=%d: |S|=5 %d/%d (%.3f)  all %d/%d (%.3f)"
              % (other_N, stab["s5_intersect"], stab["s5_this"], stab["s5_stable_frac"],
                 stab["all_intersect"], stab["all_this"], stab["all_stable_frac"]))

    # COARSE GRID completeness backstop
    g = np.linspace(0.0, 1.0, GRID_N)
    gg, _, _ = scan_growth(r_sorted, N, g, WINDOWS)
    grid_div = np.flatnonzero(np.isfinite(gg) & (gg >= DIVERGE_THRESH))
    strong_cS = sorted(float(cvals[i]) for i in div_idx if growth_cand[i] >= 30.0)
    def explained(gv):
        if gv > 0.985: return "endpoint"
        if float(np.min(np.abs(cvals - gv))) < 1e-3: return "cS"
        if any(0.0 <= (sc - gv) <= 0.02 for sc in strong_cS): return "accum-edge"
        return "UNEXPLAINED"
    lab_count = Counter(explained(float(g[gi])) for gi in grid_div)
    out["grid_completeness"] = {
        "grid_n": GRID_N, "grid_div_count": int(len(grid_div)),
        "label_counts": dict(lab_count),
        "unexplained": [round(float(g[gi]), 5) for gi in grid_div if explained(float(g[gi])) == "UNEXPLAINED"],
        "all_explained": lab_count.get("UNEXPLAINED", 0) == 0,
    }
    print("grid: %d div labels=%s all_explained=%s"
          % (out["grid_completeness"]["grid_div_count"], out["grid_completeness"]["label_counts"],
             out["grid_completeness"]["all_explained"]))

    out["elapsed_s"] = round(time.time() - t0, 1)
    fp = SLUG + "/scan_s5_N%d.json" % N
    def _js(o):
        if isinstance(o, np.bool_): return bool(o)
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        return str(o)
    with open(fp, "w") as fh:
        json.dump(out, fh, indent=2, default=_js)
    print("WROTE", fp, " elapsed=%.1fs" % out["elapsed_s"])


if __name__ == "__main__":
    main()
