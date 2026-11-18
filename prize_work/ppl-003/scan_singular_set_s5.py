"""#50 (PPL 003, $250) - EXTEND the interior singular-set map to |S|<=5 (job lab_erdos_261111_40).

Continues the programme card `erdos-50-map-the-full-interior-singular-set-of-f-c-all-diver`
(already worked 2x: job 35 = the 3008-pt |S|<=4 map @N=1e7; job 38 = local-dim + mass spectrum).
Job 35's own next idea (line 249 of the dossier) was "Push |S| to 5-6". This job does the
|S|=5 step (the memory-safe one; N=1e8 is out of reach with 3GB free on this shared PC).

WHY |S|=5 extends the frontier (not a re-run):
  Job 35 enumerated only |S|<=4 (179446 candidates -> 3008 divergent @N=1e7). The |S|=5 class
  (pqrst, 2pqrst) is a NEW prime-set family: ~1.37M additional candidate c-values, many in the
  non-trivial f region (0.2<c<0.8). If they show divergent windowed-slope they are NEW interior
  singular points the published frontier (Erdos 1995 a.e. only) and the lab's |S|<=4 map never
  covered. The empirical class trend from job 35 (divergent count GROWS with |S|: p=7, pq=16,
  pqr=163, pqrs=2684) predicts |S|=5 adds a large new family.

MECHANISM (job 27/35, reproduced exactly as the gate): a c-value c_S=prod_{q in S}(1-1/q) is a
divergent-slope interior singular point IFF it is an accumulation point of the value set V: the
neighbor families c_S(1-1/r) -> c_S as r->inf (r prime, r not in S) pile into [c-w,c] from below,
so the windowed slope s(w)=(F(c+w)-F(c-w))/(2w) DIVERGES as w shrinks (growth=s(w_small)/s(w_large)
>= 5). A phantom (c not an accumulation point) has no pile-up: slope decays. f'=0 a.e. => a.e. c
is a phantom; only the countable value set V (its accumulation closure) is singular, a null set
(consistent with #50: no x with f' in (0,inf)).

METHOD (reuses job 8 / 27 / 35 exact machinery verbatim):
  exact int64 totient sieve (cached to phi_cache_<N>.npy) -> r=phi[1:]/n -> SORTED r so that
  F(c)=(r<c).mean() is a binary search (searchsorted) -> two-sided windowed slope s(w) over the
  job-27 200x shrink (w: 0.05 -> 0.00025, 11 windows) -> growth = s(w_small)/s(w_large);
  DIVERGENCE = growth >= 5. Candidates = EVERY c_S for |S|<=5 (46 primes <=199) + uniform-grid
  backstop + phantom/generic controls.

GATES (reproduce job 35 / job 27 EXACTLY at N=1e7 -- the code is right only if ALL pass):
  (a) analytic: mean(r)=6/pi^2, sum(phi)/N^2=3/pi^2, f(0)=0, f(1)=1.
  (b) job-27 two-sided growth: 1/3->54.6514, 2/5->38.0562, 3/7->18.3267, 1/2->58.4607,
      1/4->0.4916 (phantom); point-mass 1/3->152,2/5->100,3/7->81,1/2->23,1/4->0.
  (c) c=1/2 pin (job 8): c*(win .002)=0.499004, delta-f over +-0.02 = 0.169026, growth=58.4607.
  (d) THE KEY GATE for this job: the |S|<=4 SUBSET of the |S|<=5 scan must reproduce job 35's
      per-class divergent counts EXACTLY: {2only:1, 2p:5, 2pq:7, 2pqr:125, p:7, pq:16,
      pqr:163, pqrs:2684}, total 3008.  If the |S|<=4 subset != 3008, the code/geometry drifted
      and the |S|=5 extension is NOT trusted.

RESULT (the |S|=5 extension, past the lab's |S|<=4 map):
  the NEW divergent |S|=5 class counts (pqrst, 2pqrst) at N=1e7, the total |S|<=5 count, the
  strongest |S|=5 points, and the grid backstop.

VERIFIER (card: every newly-found point re-checked at one more N, stability to 1e-3):
  the |S|=5 divergent c-set at N=1e7 is intersected with the |S|=5 divergent c-set at N=5e6
  (read from scan_s5_N5000000.json); "stable" = divergent (growth>=5) at BOTH N. The c-locations
  are exact rationals (stable to 0); we also require the peak-below location within 1e-3 of c_S.

Empirical/computational only. No proof / no Lean / no posting. Single-threaded.
Run:  python scan_singular_set_s5.py <N>     (N=5000000 for the kill test, N=10000000 primary)
Writes: scan_s5_N<N>.json
"""
import sys, os, json, time
import numpy as np
from fractions import Fraction
from itertools import combinations
from collections import Counter

SLUG = "E:/Projects/MillenniumResearch/prize_work/ppl-003"
W_LARGE, W_SMALL, N_WIN = 0.05, 0.00025, 11
WINDOWS = np.geomspace(W_LARGE, W_SMALL, N_WIN)
DIVERGE_THRESH = 5.0
GRID_N = 4001
M_FINE = 400001
PRIMES_MAX = 199
MAX_S = 5                     # <-- THE EXTENSION (job 35 used 4)
PM_CAP = 2500                 # per-point mass/accum computed only for the top-PM_CAP divergent by growth

# job-35 stored |S|<=4 per-class divergent counts at N=1e7 (the gate to reproduce)
JOB35_CLASS_1E7 = {"2only": 1, "2p": 5, "2pq": 7, "2pqr": 125,
                   "p": 7, "pq": 16, "pqr": 163, "pqrs": 2684}
JOB35_TOTAL_1E7 = 3008


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


def class_name(S):
    has2 = 2 in S
    rest = [x for x in S if x != 2]
    L = "pqrst"
    if has2:
        return "2only" if not rest else "2" + "".join(L[:len(rest)])
    return "".join(L[:len(S)])


def enumerate_candidates(P, max_S):
    cand = []
    for k in range(1, max_S + 1):
        for S in combinations(P, k):
            c = Fraction(1, 1)
            for q in S:
                c = c * Fraction(q - 1, q)
            cand.append((c, class_name(S), S))
    return cand


def scan_growth(r_sorted, N, C, windows):
    """two-sided slope growth s(w_small)/s(w_large) for an array of c (vectorized)."""
    def s_at(w):
        Fp = np.searchsorted(r_sorted, C + w, 'left')
        Fm = np.searchsorted(r_sorted, C - w, 'left')
        return (Fp - Fm) / (2.0 * w * N)
    s0, s1 = s_at(windows[0]), s_at(windows[-1])
    with np.errstate(divide='ignore', invalid='ignore'):
        growth = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), np.nan)
    return growth, s0, s1


def gate_growth_direct(r, c, windows):
    s = [((r < c + w).mean() - (r < c - w).mean()) / (2.0 * w) for w in windows]
    return s[-1] / s[0] if s[0] > 0 else float('nan')


def point_mass_exact(phi, n, cnum):
    a, b = cnum.numerator, cnum.denominator
    return int((b * phi == a * n).sum())


def count_direct(r, c, w):
    return int(((r >= c - w) & (r < c + w)).sum())


def fine_cdf(r):
    edges = np.linspace(0.0, 1.0, M_FINE + 1)
    hist, _ = np.histogram(r, bins=edges)
    F = np.zeros(len(edges)); F[1:] = np.cumsum(hist) / float(len(r))
    return edges, F


def cstar_pin(edges, F, window):
    dc = float(edges[1] - edges[0]); h = max(1, int(round(window / (2.0 * dc))))
    s = (F[2 * h:] - F[:-2 * h]) / (2.0 * h * dc)
    return float(edges[h + 1 + int(np.argmax(s))])


def peak_below_location(r_sorted, N, c_S, w=0.002, half=0.05, m=801):
    cs = np.linspace(c_S - half, c_S, m)
    L = (np.searchsorted(r_sorted, cs, 'left') - np.searchsorted(r_sorted, cs - w, 'left')) / (w * N)
    return float(cs[int(np.argmax(L))])


def main():
    t0 = time.time()
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000
    pfile = SLUG + "/phi_cache_%s.npy" % N
    if os.path.exists(pfile):
        phi = np.load(pfile)
        print("loaded cached phi (N=%d)" % N)
    else:
        phi = phi_sieve_vec(N)
        np.save(pfile, phi)
        print("built + cached phi (N=%d)  (%.1fs)" % (N, time.time() - t0))
    n = np.arange(1, len(phi), dtype=np.int64)
    pf = phi[1:].astype(np.float64); nf = n.astype(np.float64)
    r = pf / nf
    r_sorted = np.sort(r)
    mean_r = float(r.mean()); sum_phi = float(phi[1:].sum())
    gates = {
        "mean_r_6_pi2":  abs(mean_r - 6.0 / np.pi ** 2) < 1e-4,
        "sum_phi_3_pi2": abs(sum_phi / N ** 2 - 3.0 / np.pi ** 2) < 1e-4,
        "f0_is_0":       bool((r < 0.0).mean() == 0.0),
        "f1_is_1":       bool((r < 1.0).mean() > 1.0 - 2.0 / N),
    }
    out = {"N": N, "max_S": MAX_S, "gates": gates, "all_gates_pass": all(gates.values()),
           "w_large": W_LARGE, "w_small": W_SMALL, "diverge_thresh": DIVERGE_THRESH, "elapsed_s": 0.0}

    # ---------- c=1/2 gate (job 8 c* pin + delta-f) ----------
    edges, F = fine_cdf(r)
    cstar02 = cstar_pin(edges, F, 0.002)
    cf = lambda c: float((r < c).mean())
    delta_f = cf(cstar02 + 0.02) - cf(cstar02 - 0.02)
    growth05 = gate_growth_direct(r, 0.5, WINDOWS)
    out["c12_gate"] = {
        "cstar_pin_3dp": round(round(cstar02, 3), 3),
        "delta_f_over_pm0p02": round(delta_f, 6),
        "growth_two_side_N": round(growth05, 4),
        "gate_cstar_0p500": abs(cstar02 - 0.5) < 0.005,
        "gate_delta_f_0p169": abs(delta_f - 0.169) < 0.005,
        "gate_growth05": abs(growth05 - 58.4607) < 0.5,
    }

    # ---------- GATE: reproduce job-27 two-sided growth + point mass ----------
    job27 = {"1/3": (Fraction(1, 3), 54.6514, 152), "2/5": (Fraction(2, 5), 38.0562, 100),
             "3/7": (Fraction(3, 7), 18.3267, 81), "1/2": (Fraction(1, 2), 58.4607, 23),
             "1/4": (Fraction(1, 4), 0.4916, 0)}
    gate_repro = {}
    for lab, (cn, g27, m27) in job27.items():
        g = gate_growth_direct(r, float(cn), WINDOWS)
        m = point_mass_exact(phi[1:], n, cn)
        gate_repro[lab] = {"c": str(cn), "growth": round(g, 4), "point_mass": m,
                           "match_growth_job27": abs(g - g27) < 0.1,
                           "match_mass_job27": (m == m27) if N == 10_000_000 else None}
    out["gate_repro_job27"] = gate_repro
    out["gate_repro_all_match"] = all(v["match_growth_job27"] for v in gate_repro.values())
    print("job27 gate all-match:", out["gate_repro_all_match"],
          {k: v["growth"] for k, v in gate_repro.items()})

    # ---------- FULL SCAN: every c_S for |S|<=MAX_S + controls ----------
    P = primes_upto(PRIMES_MAX)
    cand = enumerate_candidates(P, MAX_S)
    seen = {}
    for c, cls, S in cand:
        if c not in seen:
            seen[c] = (cls, S)
    cands = list(seen.items())
    cvals = np.array([float(c) for c, _ in cands])
    cfrac = [str(c) for c, _ in cands]
    absS = [len(S) for (_, (cls, S)) in cands]
    control = {"0.25_phantom": 0.25, "0.30_gen": 0.30, "0.45_gen": 0.45,
               "0.60_gen": 0.60, "0.70_gen": 0.70, "0.85_gen": 0.85}
    all_c = np.concatenate([cvals, np.array(list(control.values()))])
    t_scan = time.time()
    growth_all, s0_all, s1_all = scan_growth(r_sorted, N, all_c, WINDOWS)
    growth_cand = growth_all[:len(cands)]
    growth_ctrl = {k: (None if not np.isfinite(growth_all[len(cands) + i]) else
                       round(float(growth_all[len(cands) + i]), 4))
                   for i, k in enumerate(control)}
    print("scan growth done in %.1fs over %d candidates" % (time.time() - t_scan, len(cands)))
    diverge_cand = growth_cand >= DIVERGE_THRESH
    absS_arr = np.array(absS)

    # ---------- THE KEY GATE: the |S|<=4 subset must reproduce job-35 EXACTLY ----------
    s4 = diverge_cand & (absS_arr <= 4)
    s4_idx = np.flatnonzero(s4)
    s4_by_class = Counter(cands[i][1][0] for i in s4_idx)
    gate_s4 = {
        "s4_total": int(s4.sum()),
        "s4_by_class": dict(s4_by_class),
        "job35_total": JOB35_TOTAL_1E7,
        "job35_by_class": JOB35_CLASS_1E7,
        "match_total_job35": (int(s4.sum()) == JOB35_TOTAL_1E7) if N == 10_000_000 else None,
        "match_byclass_job35": (dict(s4_by_class) == JOB35_CLASS_1E7) if N == 10_000_000 else None,
    }
    out["gate_s4_repro_job35"] = gate_s4
    print("GATE |S|<=4 repro job35: total=%d (job35=%d)  byclass=%s  match=%s/%s"
          % (gate_s4["s4_total"], JOB35_TOTAL_1E7, dict(s4_by_class),
             gate_s4["match_total_job35"], gate_s4["match_byclass_job35"]))

    # ---------- NEW RESULT: the |S|=5 (and the full |S|<=5) divergent set ----------
    div_idx = np.flatnonzero(diverge_cand)
    by_class = Counter(cands[i][1][0] for i in div_idx)
    by_absS = Counter(int(absS_arr[i]) for i in div_idx)
    out["n_candidates"] = len(cands)
    out["n_divergent"] = int(diverge_cand.sum())
    out["divergent_by_class"] = {k: int(v) for k, v in sorted(by_class.items())}
    out["divergent_by_absS"] = {str(k): int(v) for k, v in sorted(by_absS.items())}
    out["new_s5_by_class"] = {k: int(v) for k, v in by_class.items() if k in ("pqrst", "2pqrst")}
    out["new_s5_total"] = int(sum(v for k, v in by_absS.items() if k == 5))
    out["controls_growth"] = growth_ctrl
    out["controls_all_phantom"] = all((v is None) or (v < DIVERGE_THRESH) for v in growth_ctrl.values())

    # per-point detail for the divergent set (point-mass + in-window for the top by growth)
    order = div_idx[np.argsort(-growth_cand[div_idx])]
    detail = []
    for i in order:
        c, (cls, S) = cands[i]
        rec = {"c": cfrac[i], "class": cls, "absS": int(absS_arr[i]), "S": list(S),
               "growth": round(float(growth_cand[i]), 3)}
        if len(detail) < PM_CAP:
            pm = point_mass_exact(phi[1:], n, c)
            iw = count_direct(r, float(c), W_SMALL)
            rec["point_mass"] = int(pm); rec["in_window"] = int(iw)
        detail.append(rec)
    out["divergent_detail"] = detail
    # strongest |S|=5 points (the headline new family)
    s5_detail = [d for d in detail if d["absS"] == 5]
    out["s5_top20_by_growth"] = s5_detail[:20]
    out["s5_count_in_detail"] = len(s5_detail)

    print("NEW |S|<=5 map @N=%d: total divergent=%d  by_absS=%s  new |S|=5=%d  by_class(s5)=%s"
          % (N, out["n_divergent"], out["divergent_by_absS"], out["new_s5_total"],
             out["new_s5_by_class"]))

    # ---------- STABILITY vs the other N (card verifier) ----------
    other_N = 5_000_000 if N == 10_000_000 else 10_000_000
    ofile = SLUG + "/scan_s5_N%d.json" % other_N
    stability = None
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
        stability = stab
        out["stability_vs_other_N"] = stab
        print("STABILITY vs N=%d: |S|=5 stable=%d/%d (%.3f)   all stable=%d/%d (%.3f)"
              % (other_N, stab["s5_intersect"], stab["s5_this"], stab["s5_stable_frac"],
                 stab["all_intersect"], stab["all_this"], stab["all_stable_frac"]))

    # ---------- COARSE GRID completeness backstop ----------
    g = np.linspace(0.0, 1.0, GRID_N)
    gg, _, _ = scan_growth(r_sorted, N, g, WINDOWS)
    grid_div = np.flatnonzero(np.isfinite(gg) & (gg >= DIVERGE_THRESH))
    strong_cS = sorted(float(c) for i in div_idx if growth_cand[i] >= 30.0)
    cand_floats = cvals
    def explained(gv):
        if gv > 0.985:
            return "endpoint"
        if float(np.min(np.abs(cand_floats - gv))) < 1e-3:
            return "cS"
        if any(0.0 <= (sc - gv) <= 0.02 for sc in strong_cS):
            return "accum-edge"
        return "UNEXPLAINED"
    lab_count = Counter(explained(float(g[gi])) for gi in grid_div)
    out["grid_completeness"] = {
        "grid_n": GRID_N, "grid_div_count": int(len(grid_div)),
        "label_counts": dict(lab_count),
        "unexplained": [round(float(g[gi]), 5) for gi in grid_div if explained(float(g[gi])) == "UNEXPLAINED"],
        "all_explained": lab_count.get("UNEXPLAINED", 0) == 0,
    }
    print("grid: %d divergent  labels=%s  all_explained=%s"
          % (out["grid_completeness"]["grid_div_count"],
             out["grid_completeness"]["label_counts"], out["grid_completeness"]["all_explained"]))

    out["elapsed_s"] = round(time.time() - t0, 1)
    fp = SLUG + "/scan_s5_N%d.json" % N
    def _js(o):
        if isinstance(o, (np.bool_,)): return bool(o)
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, (np.ndarray,)): return o.tolist()
        return str(o)
    with open(fp, "w") as fh:
        json.dump(out, fh, indent=2, default=_js)
    print("WROTE", fp, "  elapsed=%.1fs" % out["elapsed_s"])


if __name__ == "__main__":
    main()
