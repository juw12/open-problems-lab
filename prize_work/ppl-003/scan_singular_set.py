"""#50 (PPL 003, $250) - map the FULL interior singular set of f(c) (job lab_erdos_261111_35).

Erdos #50: Schoenberg density f(c)=dens{n: phi(n)/n < c} is PURELY SINGULAR
(f'=0 a.e., Erdos 1995 Resenhas; erdosproblems.com/50, confirmed open 2026-11-11).
The pointwise question "any x with f'(x) exists and positive?" is OPEN.
Published computational frontier: ONLY the tail (Weingartner 2007, c->0) and the
endpoint (Tenenbaum-Toulmonde 2006, c->1). NO interior c in (0,1) is characterized,
and no published result identifies WHERE the singular measure concentrates.

LAB prior best (jobs 8/22/24/27): the interior singular set was known only at
  c=1/2           (strongest; jobs 8/22/24, delta-f=0.169, N*eps=10)
  c=(p-1)/(2p)    for p=3,5,7 i.e. 1/3, 2/5, 3/7  (job 27)
  c=0.25          EXCLUDED as a phantom (no n has phi(n)/n=1/4, E=1.0).

THIS JOB (idea card: map the FULL interior singular set, not just c=1/2):
  Scan the interior c in (0,1) and enumerate EVERY divergent-slope interior
  singular point at finite N, i.e. every c_S = prod_{q in S}(1-1/q) (S a finite
  prime set) that shows a divergent windowed-slope, and classify each by prime
  set S (2p-type, odd-prime-power-type, 2pq-type, ...). Extend the lab's
  {1/3,2/5,3/7,1/2} to the full detectable family AND to the untouched (1/2,1]
  region (odd-prime-power type (p-1)/p for odd p).

MECHANISM (job 27, reproduced exactly as the gate): a c-value is a divergent-
slope interior singular point IFF it is an ACCUMULATION point of the value set
V: the neighbor families c_S(1-1/r) -> c_S as r->inf (r prime, r not in S) pile
into [c-w, c] from below, so the windowed slope s(w)=(F(c+w)-F(c-w))/(2w)
DIVERGES as the window w shrinks. A phantom (a c not in V, e.g. 0.25) has no such
pile-up: the slope decays, no divergence. f'=0 a.e. => a.e. c is a phantom; only
the countable set V (its accumulation closure) is singular, a null set (consistent
with #50: no x with f' in (0,inf)).

METHOD (reuses job 8 / job 27 exact machinery):
  exact int64 totient sieve -> r=phi[1:]/n (float64) -> SORTED r so that
  F(c)=(r<c).mean() is a binary search (searchsorted) -> two-sided windowed slope
  s(w) over a 200x shrink (w: 0.05 -> 0.00025, 11 windows, same geometry as job 27)
  -> growth = s(w_small)/s(w_large); DIVERGENCE = growth >= 5.
  Candidates = EVERY c_S for |S|<=3 (46 primes <=199) = the complete set of
  detectable-magnitude singular values, + a coarse uniform grid over [0,1]
  (backstop for |S|>=4 / anything unexpected) + phantom/generic controls.

GATES (reproduce job 27 EXACTLY at N=1e7, two-sided slope growth over 200x):
  1/3 -> 54.65, 2/5 -> 38.06, 3/7 -> 18.33, 1/2 -> 58.46, 0.25 -> 0.49 (phantom).
  tower discrete point-mass (integer): 1/3->152, 2/5->100, 3/7->81, 1/2->23, 0.25->0.
  c=1/2 c*-pin (job 8, window 0.002) -> 0.499004; delta-f over c*+/-0.02 -> 0.169.

STABILITY (card verifier): every NEW point re-checked at one more N (N=2e7); a
point is "stable" if it diverges (growth>=5) at BOTH N. The c-value is an exact
rational (location stable to 0 < 1e-3); we also measure the empirical one-sided-
below peak location and require |peak - c_S| < 1e-3.

Empirical/computational only. No proof / no Lean / no posting. Single-threaded.
Run:  python scan_singular_set.py <N>
Writes:  scan_N<N>.json
"""
import sys, os, json, time
import numpy as np
from fractions import Fraction
from itertools import combinations

SLUG = "E:/Projects/MillenniumResearch/prize_work/ppl-003"
W_LARGE, W_SMALL, N_WIN = 0.05, 0.00025, 11
WINDOWS = np.geomspace(W_LARGE, W_SMALL, N_WIN)
DIVERGE_THRESH = 5.0
GRID_N = 4001               # backstop uniform grid over [0,1] (step 2.5e-4)
M_FINE = 400001             # job-8 fine CDF bins for the c* pin
PRIMES_MAX = 199

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

def enumerate_candidates(P):
    cand = []
    for p in P:
        cand.append((Fraction(p - 1, p), class_name((p,)), (p,)))
    for p, q in combinations(P, 2):
        c = Fraction(p - 1, p) * Fraction(q - 1, q)
        cand.append((c, class_name((p, q)), (p, q)))
    for p, q, r in combinations(P, 3):
        c = Fraction(p - 1, p) * Fraction(q - 1, q) * Fraction(r - 1, r)
        cand.append((c, class_name((p, q, r)), (p, q, r)))
    for p, q, r, s in combinations(P, 4):
        c = Fraction(p - 1, p) * Fraction(q - 1, q) * Fraction(r - 1, r) * Fraction(s - 1, s)
        cand.append((c, class_name((p, q, r, s)), (p, q, r, s)))
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
    """job-27 EXACT formula (direct mask) for one c: s(w_small)/s(w_large)."""
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
    """empirical one-sided-below slope peak in [c_S-half, c_S]; expect ~ c_S."""
    cs = np.linspace(c_S - half, c_S, m)
    L = (np.searchsorted(r_sorted, cs, 'left') - np.searchsorted(r_sorted, cs - w, 'left')) / (w * N)
    return float(cs[int(np.argmax(L))])

def main():
    t0 = time.time()
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000
    # cache the slow, CPU-bound phi sieve to disk so re-runs are instant (the PC is
    # shared with 8 parallel workers; the sieve is the only expensive step)
    pfile = SLUG + "/phi_cache_%s.npy" % N
    if os.path.exists(pfile):
        phi = np.load(pfile)
        print("loaded cached phi (N=%d) from %s" % (N, os.path.basename(pfile)))
    else:
        phi = phi_sieve_vec(N)
        np.save(pfile, phi)
        print("built + cached phi (N=%d) -> %s  (%.1fs)" % (N, os.path.basename(pfile), time.time() - t0))
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
    out = {"N": N, "gates": gates, "all_gates_pass": all(gates.values()),
           "windows": [round(float(w), 8) for w in WINDOWS],
           "w_large": W_LARGE, "w_small": W_SMALL, "diverge_thresh": DIVERGE_THRESH,
           "elapsed_s": 0.0}

    # ---------- c=1/2 gate (job 8 c* pin + delta-f) + job-27 reference growth ----------
    edges, F = fine_cdf(r)
    cstar02 = cstar_pin(edges, F, 0.002)
    cstar01 = cstar_pin(edges, F, 0.01)
    cstar1e4 = cstar_pin(edges, F, 1e-4)
    cf = lambda c: float((r < c).mean())
    delta_f = cf(cstar02 + 0.02) - cf(cstar02 - 0.02)
    growth05 = gate_growth_direct(r, 0.5, WINDOWS)
    out["c12_gate"] = {
        "cstar_win0p01": round(cstar01, 6), "cstar_win0p002": round(cstar02, 6),
        "cstar_win1e4": round(cstar1e4, 6), "cstar_pin_3dp": round(round(cstar02, 3), 3),
        "delta_f_over_pm0p02": round(delta_f, 6),
        "growth_two_side_N": round(growth05, 4),
        "gate_cstar_0p500": abs(cstar02 - 0.5) < 0.005,
        "gate_delta_f_0p169": abs(delta_f - 0.169) < 0.005,
        "gate_growth05": abs(growth05 - 58.4607) < 0.5,
    }

    # ---------- GATE: reproduce job-27 two-sided growth + point mass (N=1e7 expected) ----------
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

    # ---------- FULL SCAN: every c_S for |S|<=3 + controls ----------
    P = primes_upto(PRIMES_MAX)
    cand = enumerate_candidates(P)
    # de-dup by value (keep smallest S) - values are unique per S actually, but be safe
    seen = {}
    for c, cls, S in cand:
        if c not in seen:
            seen[c] = (cls, S)
    cands = list(seen.items())
    cvals = np.array([float(c) for c, _ in cands])
    cfrac = [str(c) for c, _ in cands]
    # controls: known phantom + generic points
    control = {"0.25_phantom": 0.25, "0.30_gen": 0.30, "0.45_gen": 0.45,
               "0.60_gen": 0.60, "0.70_gen": 0.70, "0.85_gen": 0.85}
    all_c = np.concatenate([cvals, np.array(list(control.values()))])
    growth_all, s0_all, s1_all = scan_growth(r_sorted, N, all_c, WINDOWS)
    growth_cand = growth_all[:len(cands)]
    growth_ctrl = {k: (None if not np.isfinite(growth_all[len(cands) + i]) else
                       round(float(growth_all[len(cands) + i]), 4))
                   for i, k in enumerate(control)}
    diverge_cand = growth_cand >= DIVERGE_THRESH

    # neighbor accumulation at the finest window for a reported point
    def neighbor_accum(cn, w=W_SMALL):
        tot = count_direct(r, float(cn), w)
        pm = point_mass_exact(phi[1:], n, cn)
        return tot, pm, (tot - pm)

    # reported divergent candidates grouped by class
    div_by_class = {}
    divergent_list = []
    for i, (c, (cls, S)) in enumerate(cands):
        if diverge_cand[i]:
            tot, pm, na = neighbor_accum(c)
            rec = {"c": cfrac[i], "class": cls, "S": S,
                   "growth": round(float(growth_cand[i]), 3),
                   "point_mass": int(pm), "in_window": int(tot),
                   "neighbor_accum": int(na)}
            divergent_list.append(rec)
            div_by_class.setdefault(cls, []).append(rec)
    out["n_candidates"] = len(cands)
    out["n_divergent_candidates"] = int(diverge_cand.sum())
    out["divergent_by_class"] = {k: len(v) for k, v in sorted(div_by_class.items())}
    out["controls_growth"] = growth_ctrl
    out["controls_all_phantom"] = all((v is None) or (v < DIVERGE_THRESH) for v in growth_ctrl.values())
    # the full divergent list (sorted by class then c), capped to keep JSON small
    out["divergent_list"] = sorted(divergent_list, key=lambda d: (d["class"], d["c"]))
    out["divergent_list_head"] = out["divergent_list"][:300]

    # ---------- PEAK-LOCATION (informational): for the strongest point per class,
    # confirm the mass peaks at c_S (well-separated points peak exactly at c_S;
    # dense-cluster points are located by their exact rational c_S, stable to 0) ----
    pk = {}
    for cls in div_by_class:
        best = max(div_by_class[cls], key=lambda d: d["growth"])
        cS = Fraction(best["c"])
        pkf = peak_below_location(r_sorted, N, float(cS), w=0.0005, half=0.003, m=601)
        pk[cls] = {"example_c": best["c"], "growth": best["growth"],
                   "peak_below": round(pkf, 6),
                   "abs_peak_minus_c": round(abs(pkf - float(cS)), 6),
                   "peak_at_cS_within_1e3": abs(pkf - float(cS)) < 1e-3}
    out["peak_location_by_class"] = pk

    # ---------- COARSE GRID completeness backstop ----------
    g = np.linspace(0.0, 1.0, GRID_N)
    gg, _, _ = scan_growth(r_sorted, N, g, WINDOWS)
    grid_div = np.flatnonzero(np.isfinite(gg) & (gg >= DIVERGE_THRESH))
    # strong points (their neighbor-accumulation fills the region just below them)
    strong_cS = sorted(float(c) for i, (c, _c2) in enumerate(cands)
                       if np.isfinite(growth_cand[i]) and growth_cand[i] >= 30.0)
    cand_floats = np.array([float(c) for c, _ in cands])
    def explained(gv):
        if gv > 0.985:                      # c->1 endpoint cluster (Tenenbaum-Toulmonde)
            return "endpoint"
        d_cS = float(np.min(np.abs(cand_floats - gv)))
        if d_cS < 1e-3:                     # a discrete tower value c_S
            return "cS"
        # within the accumulation region just below a strong c_S
        if any(0.0 <= (sc - gv) <= 0.02 for sc in strong_cS):
            return "accum-edge"
        return "UNEXPLAINED"
    grid_matches = []
    for gi in grid_div:
        gv = float(g[gi]); lab = explained(gv)
        grid_matches.append({"grid_c": round(gv, 5), "growth": round(float(gg[gi]), 3),
                             "label": lab})
    from collections import Counter
    lab_count = Counter(m["label"] for m in grid_matches)
    out["grid_completeness"] = {
        "grid_n": GRID_N, "grid_div_count": int(len(grid_div)),
        "label_counts": dict(lab_count),
        "unexplained": [m for m in grid_matches if m["label"] == "UNEXPLAINED"],
        "all_explained": lab_count.get("UNEXPLAINED", 0) == 0,
        "divergent_grid_points": grid_matches[:200],
    }

    out["elapsed_s"] = round(time.time() - t0, 1)
    fp = SLUG + "/scan_N%s.json" % N
    def _js(o):
        if isinstance(o, (np.bool_,)): return bool(o)
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, (np.ndarray,)): return o.tolist()
        return str(o)
    with open(fp, "w") as fh:
        json.dump(out, fh, indent=2, default=_js)

    # ---------- compact hidden-number-safe summary ----------
    print("N=%d  gates=%s  elapsed=%.1fs" % (N, all(gates.values()), out["elapsed_s"]))
    print("c=1/2 gate: c*(win .002)=%.6f dF=%.4f growth05=%.3f -> %s"
          % (cstar02, delta_f, growth05,
             out["c12_gate"]["gate_cstar_0p500"] and out["c12_gate"]["gate_delta_f_0p169"]))
    print("job27 gate repro (growth / mass):")
    for lab, v in gate_repro.items():
        print("  %-4s growth=%7.3f (job27)  mass=%4d  match=%s"
              % (lab, v["growth"], v["point_mass"], v["match_growth_job27"]))
    print("  all growth match:", out["gate_repro_all_match"])
    print("candidates |S|<=3: %d ; divergent: %d" % (len(cands), out["n_divergent_candidates"]))
    print("divergent by class:", out["divergent_by_class"])
    print("controls growth:", growth_ctrl, "all_phantom:", out["controls_all_phantom"])
    print("grid: %d divergent  labels=%s  all_explained=%s"
          % (out["grid_completeness"]["grid_div_count"],
             out["grid_completeness"]["label_counts"],
             out["grid_completeness"]["all_explained"]))
    for cls, d in pk.items():
        print("  class %-4s example=%-9s growth=%7.2f peak-c=%.6f at_cS=%s"
              % (cls, d["example_c"], d["growth"], d["abs_peak_minus_c"], d["peak_at_cS_within_1e3"]))
    print("WROTE", fp)

if __name__ == "__main__":
    main()
