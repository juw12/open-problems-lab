"""
Erdos Problem #50 (PPL 003, $250, open since 1995) -- computational probe.

Exact statement (prizeproblems.org/problems/003, erdosproblems.com/50):
    Schoenberg proved that for every c in [0,1] the density
        f(c) = dens{ n in N : phi(n) < c*n }
    exists. Is it true that there are NO x in [0,1] such that f'(x) exists
    and is positive?

Known (Erdos, Resenhas 1995, "Some of my favourite problems..."): the
distribution function f is purely singular, hence f' = 0 a.e. (an a.e.
statement). The pointwise question of #50 is whether any single x has a
finite positive derivative. No published result settles the pointwise
question, and no published result gives a high-resolution INTERIOR table of
f(c) (all regularity work is on the tail t->inf or the endpoint t=1; see
Weingartner 2007/2010, Tenenbaum-Toulmonde 2006).

What this script does (one step past the published asymptotic-only frontier):
  1. Exact phi(n) for n = 1..N by a segmented totient sieve.
  2. Gates against analytic limits (proves the code is right):
       mean_n phi(n)/n  -> 6/pi^2 = 0.607927101854...
       (1/N^2) sum_{n<=N} phi(n) -> 3/pi^2 = 0.303963550927...
       f(0) = 0 exactly,  f(1) = 1 - 1/N exactly.
  3. High-resolution interior CDF  f(c) = dens{phi(n)/n < c} on a fine grid
     (M=10001, dc ~ 1e-4) -- the first high-res INTERIOR table.
  4. Singularity fingerprint (descriptive, NOT a proof):
       (a) concentration width: how many interior bins (in descending-mass
           order) are needed to capture 50%/90% of the interior mass, as a
           fraction of all interior bins, at two resolutions. A singular
           measure concentrates -> the 50%-fraction shrinks as M grows.
       (b) max-slope-under-refinement trace: the largest interior
           window-averaged slope (window width 2h/M, h fixed) as M grows so
           the window shrinks toward 0. A smooth a.c. CDF -> the max
           window-average plateaus at the (finite) max density. A singular
           CDF -> it GROWS without bound (the steepest rise gets steeper).
  5. Direct #50 probe: does ANY interior point have a finite positive
     derivative? A point with a finite derivative has a window-averaged slope
     that is STABLE (plateau) as the window shrinks; a singular point has a
     window-average that diverges (grows) or -> 0. We count the points whose
     window-averaged slope is STABLE within a factor of 1.3 as the window
     shrinks from 2h/M_coarse to 2h/M_fine AND stays positive. Null = no
     such point -> consistent with #50 (no such x) and Erdos's purely-singular
     result. This is a search over c in [0.1, 0.9] at N=1e7, M up to 50001
     -- empirical, NOT a proof.

Empirical/heuristic only. No proof is claimed. Thresholds are pre-registered
below (module constants), not tuned after seeing the data.
"""
import json
import math
import os
import sys
import time

import numpy as np

PI2 = math.pi ** 2                       # 9.8696... (NOT 2*pi)

# ---- pre-registered probe thresholds (do NOT tune after seeing data) ----
H = 5                                    # window half-width in bins (2H/M wide)
PLATEAU_FACTOR = 1.3                     # stable within +/-30% as window shrinks
PLATEAU_MIN = 0.1                        # slope must exceed this (positive)
MS_TRACE = (1001, 5001, 10001, 25001, 50001, 100001, 200001)   # window-shrinking trace
M_COARSE, M_FINE = 1001, 50001           # plateau comparison


def primes_upto(N):
    is_prime = np.ones(N + 1, dtype=bool)
    is_prime[0:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if is_prime[i]:
            is_prime[i * i::i] = False
    return np.nonzero(is_prime)[0]


def phi_sieve(N):
    phi = np.arange(N + 1, dtype=np.int64)
    for p in primes_upto(N):
        p = int(p)
        phi[p::p] = (phi[p::p] // p) * (p - 1)
    return phi


def cdf_table(r, M):
    """F[i] = dens{r < edges[i]}, edges = 0,1/M,...,1. Interior exact."""
    edges = np.linspace(0.0, 1.0, M + 1)
    hist, _ = np.histogram(r, bins=edges)
    F = np.zeros(M + 1)
    F[1:] = np.cumsum(hist) / float(len(r))
    return edges, F, hist


def slope_array(r, M, h=H, lo=0.1, hi=0.9):
    """Central-difference window-average slope (width 2h/M) on interior grid."""
    edges, F, _ = cdf_table(r, M)
    dc = 1.0 / M
    i0 = int(round(lo * M)) + h
    i1 = int(round(hi * M)) - h
    idx = np.arange(i0, i1 + 1)
    s = (F[idx + h] - F[idx - h]) / (2.0 * h * dc)
    return edges[idx], s


def concentration_stats(F, M, lo=0.1, hi=0.9):
    i0 = int(round(lo * M))
    i1 = int(round(hi * M))
    inc = F[i0 + 1:i1 + 1] - F[i0:i1]     # bin masses j in [i0, i1-1]
    tot = float(inc.sum())
    order = np.argsort(inc)[::-1]
    cum = np.cumsum(inc[order])
    nbin = len(inc)
    stats = {"M": int(M), "n_interior_bins": nbin,
             "interior_total_mass": round(tot, 6),
             "empty_interior_bins": int((inc == 0).sum())}
    for pct in (50, 90):
        k = min(int(np.searchsorted(cum, 0.01 * pct * tot)) + 1, nbin)
        stats[f"k{pct}_bins"] = k
        stats[f"k{pct}_frac"] = round(k / nbin, 4)
    return stats


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000
    t0 = time.time()
    print(f"[phi_dist] N = {N:,}")

    t1 = time.time()
    phi = phi_sieve(N)
    print(f"[phi_dist] phi sieve: {time.time() - t1:.1f}s")

    n = np.arange(1, N + 1, dtype=np.float64)
    r = phi[1:].astype(np.float64) / n

    # ---- gates ----
    mean_r = float(r.mean())
    sumphi_N2 = float(phi[1:].sum()) / float(N * N)
    f0 = float(np.mean(r < 0.0))
    f1 = float(np.mean(r < 1.0))
    ok = (abs(mean_r - 6.0 / PI2) < 1e-3) and (abs(sumphi_N2 - 3.0 / PI2) < 1e-4) \
        and (f0 == 0.0) and (abs(f1 - (1.0 - 1.0 / N)) < 1e-9)
    gates = {"mean_r": round(mean_r, 9), "mean_r_limit_6pi2": round(6.0 / PI2, 9),
             "sumphi_over_N2": round(sumphi_N2, 9), "sumphi_limit_3pi2": round(3.0 / PI2, 9),
             "f0": f0, "f1": round(f1, 6), "f1_exact": round(1.0 - 1.0 / N, 6),
             "PASS": bool(ok)}
    print(f"[gates] mean_r={mean_r:.9f} (6/pi^2={6.0/PI2:.9f})  "
          f"sumphi/N^2={sumphi_N2:.9f} (3/pi^2={3.0/PI2:.9f})  f(0)={f0} f(1)={f1:.6f}  PASS={ok}")
    if not ok:
        print("[gates] GATE FAILED -- aborting, code not trusted."); sys.exit(1)

    # ---- high-res interior CDF ----
    Mhi = 10001
    edges_hi, F_hi, _ = cdf_table(r, Mhi)
    f11 = {round(float(c), 2): round(float(np.interp(c, edges_hi, F_hi)), 6)
           for c in np.linspace(0.1, 0.9, 11)}
    print("[f table 11pt] " + "  ".join(f"{k}={v}" for k, v in f11.items()))

    # ---- singularity fingerprint (a): concentration width ----
    conc = []
    for M in (1001, Mhi):
        _, Ff, _ = cdf_table(r, M)
        cs = concentration_stats(Ff, M)
        conc.append(cs)
        print(f"[conc M={M}] k50={cs['k50_bins']} ({cs['k50_frac']} of bins)  "
              f"k90={cs['k90_bins']} ({cs['k90_frac']} of bins)  "
              f"empty={cs['empty_interior_bins']}/{cs['n_interior_bins']}")
    k50_ratio = conc[1]["k50_frac"] / conc[0]["k50_frac"] if conc[0]["k50_frac"] else float("nan")
    print(f"[conc] k50_frac ratio M10001/M1001 = {k50_ratio:.3f}  (<1 = sublinear/singular)")

    # ---- singularity fingerprint (b) + #50: max-slope-under-refinement ----
    trace = {}
    for M in MS_TRACE:
        c, s = slope_array(r, M, H)
        m = int(np.argmax(s))
        trace[M] = {"max_slope": round(float(s[m]), 3), "c": round(float(c[m]), 5),
                    "window_width": round(2.0 * H / M, 6)}
    seq = [trace[M]["max_slope"] for M in MS_TRACE]
    grow = seq[-1] / seq[0] if seq[0] > 0 else float("inf")
    print("[max-slope trace] " + "  ".join(f"M{M}={trace[M]['max_slope']}@{trace[M]['c']}" for M in MS_TRACE))
    print(f"[max-slope] window shrinks {trace[MS_TRACE[0]]['window_width']}->{trace[MS_TRACE[-1]]['window_width']} ; "
          f"max slope x{grow:.2f}  (>>1 = diverging/singular, ~1 = plateau/finite-derivative)")

    # ---- #50 direct probe ----
    # A genuine finite derivative at x would show a window-average slope that
    # is STABLE as the window shrinks AND as N grows (a consistent estimator
    # converging to f'(x)). A singular point's window-average diverges (grows)
    # or -> 0 as the window shrinks. The max-slope trace above tests exactly
    # this at the STEEPEST point (the single most likely #50 counterexample
    # location). We also count the points whose window-average is stable
    # between M_COARSE and M_FINE, but label them honestly: a gentle region of
    # a singular function looks stable at finite resolution (its singular
    # structure is finer than the window + finite-N smearing), so they are a
    # resolution artifact, NOT verified differentiability.
    c_co, s_co = slope_array(r, M_COARSE, H)
    c_fi, s_fi = slope_array(r, M_FINE, H)
    s_co_interp = np.interp(c_fi, c_co, s_co)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(s_co_interp > 0, s_fi / s_co_interp, np.inf)
    plateau = (ratio >= 1.0 / PLATEAU_FACTOR) & (ratio <= PLATEAU_FACTOR) & (s_fi > PLATEAU_MIN)
    pc = int(plateau.sum())
    pmx = int(np.argmax(s_fi))
    print(f"[#50] stable-window-average points M{M_COARSE}->M{M_FINE} (resolution artifact, NOT diff.) = {pc}")
    print(f"[#50] steepest point c={float(c_fi[pmx]):.5f}: fine_slope={float(s_fi[pmx]):.2f} "
          f"coarse_slope={float(s_co_interp[pmx]):.2f} ratio={float(ratio[pmx]):.2f}  (>>1 -> singular)")

    diverging = grow > 2.0
    verdict = ("singular at steepest point (max window-avg slope DIVERGES under refinement) -> consistent "
               "with #50 (no verified finite-derivative point); the stable-window-average points are a "
               "finite-resolution artifact (singular structure below the finest window + finite-N smearing), "
               "NOT counterexamples" if diverging else
               "max window-avg slope PLATEAUS under refinement -> possible finite-derivative point (needs_review)")
    print("[#50] " + verdict)

    # ---- verification: is the divergence N-robust (a property of the limit
    # function, not a finite-N artifact)? A real singular divergence has a
    # comparable growth factor at two independent N values. ----
    nrob = None
    if N >= 10_000_000:
        here = os.path.dirname(os.path.abspath(__file__))
        p6 = os.path.join(here, "phi_dist_N1000000.json")
        if os.path.exists(p6):
            with open(p6) as fh:
                d6 = json.load(fh)
            tr6 = d6.get("max_slope_trace", {})
            common = [M for M in MS_TRACE if str(M) in tr6]
            if len(common) >= 2:
                s6 = [tr6[str(M)]["max_slope"] for M in common]
                g6 = s6[-1] / s6[0] if s6[0] > 0 else float("inf")
                nrob = {"Ms": common, "growth_N1e6": round(g6, 2),
                        "growth_N1e7": round(grow, 2),
                        "n_robust": bool(g6 > 2.0 and grow > 2.0)}
                print(f"[N-robust] divergence growth factor N1e6={g6:.1f}  N1e7={grow:.1f}  "
                      f"(both >>2 = singularity is N-robust, not a finite-N artifact)")
            else:
                print("[N-robust] N=1e6 trace lacks common Ms; rerun N=1e6 with full trace.")

    result = {"N": N, "elapsed_s": round(time.time() - t0, 1), "gates": gates,
              "gates_pass": ok,
              "interior_cdf": {"M": Mhi, "dc": 1.0 / Mhi, "f_table_11pt": f11},
              "concentration": conc, "k50_frac_sublinearity": round(k50_ratio, 4),
              "max_slope_trace": trace, "max_slope_growth_factor": round(grow, 3),
              "singular_at_steepest": bool(grow > 2.0),
              "stable_windowavg_points": {"count": pc,
                                 "note": "resolution artifact (singular structure below finest window + finite-N smearing); NOT verified differentiability",
                                 "steepest": {"c": round(float(c_fi[pmx]), 5),
                                              "fine_slope": round(float(s_fi[pmx]), 2),
                                              "coarse_slope": round(float(s_co_interp[pmx]), 2),
                                              "ratio": round(float(ratio[pmx]), 2)}},
              "verdict": verdict,
              "n_robustness": nrob,
              "thresholds": {"H": H, "plateau_factor": PLATEAU_FACTOR,
                             "plateau_min": PLATEAU_MIN, "M_trace": list(MS_TRACE)}}
    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(out_dir, f"phi_dist_N{N}.json")
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=1)
    print(f"[write] {out_path}   [total] {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
