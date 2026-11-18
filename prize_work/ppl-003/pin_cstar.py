"""Pin the interior transition point c* of f(c)=dens{n: phi(n)/n < c} on the N=1e7 grid,
using EXACT empirical counts f(c)=(r<c).mean() (not np.interp, which undershoots at the
c=0.5 density spike). Verifies an independent second N.  Empirical/computational only:
no proof claims.  Single-threaded, ~15-25 s.

Run:  python pin_cstar.py <N>
Writes:  pin_cstar_N<N>.json
"""
import sys, json, time
import numpy as np

M_FINE = 400001          # fine profile bins (bin width ~2.5e-6)
KILL_DELTA_F = 0.005     # kill if |f(c*+0.02)-f(c*-0.02)| < this  (=> plateau)
KILL_HALF    = 0.02
WIN_COARSE   = 0.01      # coarse window for slope divergence
WIN_FINE     = 1e-4      # fine window for slope divergence
PIN_WINDOW   = 0.002     # window used to pin c* (matches the +/-0.002 target)

def primes_upto(n):
    if n < 2: return np.array([], dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool); sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]: sieve[i*i::i] = False
    return np.nonzero(sieve)[0]

def phi_sieve(N):
    phi = np.arange(N + 1, dtype=np.int64)
    for p in primes_upto(N):
        for j in range(p, N + 1, p):
            phi[j] = phi[j] // p * (p - 1)
    return phi

def f_exact(r, c):
    """EXACT empirical f(c) = fraction of r strictly < c (the definition)."""
    return float((r < c).mean())

def fine_cdf(r):
    edges = np.linspace(0.0, 1.0, M_FINE + 1)
    hist, _ = np.histogram(r, bins=edges)
    F = np.zeros(len(edges)); F[1:] = np.cumsum(hist) / float(len(r))
    return edges, F, float(edges[1] - edges[0])

def window_slope_max(edges, F, window):
    dc = float(edges[1] - edges[0]); M = len(edges) - 1
    h = max(1, int(round(window / (2.0 * dc))))
    s = (F[2*h:] - F[:-2*h]) / (2.0 * h * dc)
    cstar = float(edges[h + 1 + int(np.argmax(s))])
    return cstar, float(s.max()), h

def hist_selfcheck(edges, F, r):
    """Prove the histogram CDF is correct: F at the edge nearest 0.5 == direct count."""
    i = int(np.searchsorted(edges, 0.5, side='right'))   # first edge > 0.5
    F_at = float(F[i])                                    # fraction r < edges[i] (~0.5)
    direct = float((r < 0.5).mean())
    return F_at, direct, abs(F_at - direct)

def main():
    t0 = time.time()
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000
    phi = phi_sieve(N)
    n = np.arange(1, N + 1)
    r = phi[1:].astype(np.float64) / n.astype(np.float64)
    mean_r = float(r.mean()); sum_phi = float(phi[1:].sum())

    edges, F, dc = fine_cdf(r)
    F_at, direct, sc_err = hist_selfcheck(edges, F, r)

    # analytic gates
    gates = {
        "mean_r_6_pi2":   abs(mean_r - 6.0/np.pi**2) < 1e-4,
        "sum_phi_3_pi2":  abs(sum_phi/N**2 - 3.0/np.pi**2) < 1e-4,
        "f0_is_0":        f_exact(r, 0.0) == 0.0,
        "f1_is_1":        f_exact(r, 1.0) > 1.0 - 2.0/N,   # n=1 has r=1.0 exactly
    }
    allg = all(gates.values())
    # independent CDF check: every even non-power-of-2 has phi(n)/n < 0.5  => f(0.5) >= 0.5
    f05_check = f_exact(r, 0.5) >= 0.5 - 1e-9

    # job29 comparison (job29's f values are np.interp @ M=10001)
    j29 = {0.10:0.0, 0.18:0.0, 0.26:0.0180, 0.34:0.17235, 0.42:0.251825,
           0.50:0.475023, 0.58:0.555692, 0.66:0.610535, 0.74:0.687583,
           0.82:0.745559, 0.90:0.785871}
    table = []
    for c in sorted(j29):
        fv = f_exact(r, c)
        table.append({"c": c, "my_f_exact": round(fv, 6),
                      "job29_f_interp": j29[c],
                      "match_abs0005": abs(fv - j29[c]) < 0.0005})
    # the three gentle points must agree with job29 (they did); f(0.5) is the flagged one
    gentle_ok = all(t["match_abs0005"] for t in table if t["c"] in (0.34, 0.58, 0.90))

    median = float(np.quantile(r, 0.5))
    cstar, slope_pin, hpin = window_slope_max(edges, F, PIN_WINDOW)
    cs_co, sl_co, _  = window_slope_max(edges, F, WIN_COARSE)
    cs_fi, sl_fi, _  = window_slope_max(edges, F, WIN_FINE)
    growth = sl_fi / sl_co if sl_co > 0 else float('nan')

    cf = f_exact(r, cstar)
    dlo = f_exact(r, cstar - KILL_HALF); dhi = f_exact(r, cstar + KILL_HALF)
    delta_f = dhi - dlo
    kill = abs(delta_f) < KILL_DELTA_F

    profile = [{"c": round(cstar + d, 4), "f_exact": round(f_exact(r, cstar + d), 6)}
               for d in (0.02, 0.01, 0.005, 0.002, 0.001, 0.0)]

    # 2p spike: mass just below 0.5 (the local density divergence)
    spike = {}
    for e in (0.001, 0.0005, 0.00025, 0.0001, 0.00005):
        spike["eps_" + format(e, 'g')] = round(f_exact(r, 0.5) - f_exact(r, 0.5 - e), 7)

    out = {
        "N": N, "M_fine": M_FINE, "bin_width": dc,
        "gates": gates, "all_gates_pass": bool(allg),
        "hist_selfcheck": {"F_at_0p5_edge": round(F_at, 7), "direct_r_lt_0p5": round(direct, 7),
                            "abs_err": sc_err, "histogram_correct": bool(sc_err < 1e-5)},
        "f05_ge_0p5_check": bool(f05_check),
        "job29_table": table, "gentle_points_match_job29": bool(gentle_ok),
        "f05_job29_interp": 0.475023, "f05_true": round(direct, 6),
        "f05_is_interp_artifact": bool(direct - 0.475023 > 0.02),
        "median_f0p5_crossing": round(median, 6),
        "cstar": round(cstar, 6), "cstar_pin_window": PIN_WINDOW, "slope_at_cstar": slope_pin,
        "slope_divergence": {"coarse_win0p01": sl_co, "fine_win1e4": sl_fi,
                              "growth": round(growth, 3), "cs_coarse": round(cs_co, 5),
                              "cs_fine": round(cs_fi, 5)},
        "kill": {"cstar": round(cstar, 5), "f_cstar": round(cf, 6),
                 "f_cstar_minus_0p02": round(dlo, 6), "f_cstar_plus_0p02": round(dhi, 6),
                 "delta_f": round(delta_f, 6), "threshold": KILL_DELTA_F,
                 "is_plateau": bool(kill),
                 "verdict": "PLATEAU-transition-not-sharp" if kill else "SHARP-transition-confirmed"},
        "local_profile_around_cstar": profile,
        "twop_spike_mass_below_0p5": spike,
        "elapsed_s": round(time.time() - t0, 1),
    }
    fp = "pin_cstar_N%s.json" % N
    with open(fp, "w") as fh: json.dump(out, fh, indent=2)

    print("N=%d  %ds  M=%d  bin=%.2e" % (N, out["elapsed_s"], M_FINE, dc))
    print("gates:", gates, "ALL=%s" % allg)
    print("hist self-check: F@0.5-edge=%.6f direct=%.6f err=%.2e correct=%s"
          % (F_at, direct, sc_err, out["hist_selfcheck"]["histogram_correct"]))
    print("f(0.5) true=%.6f  (job29 interp=0.475023)  f05>=0.5=%s  artifact=%s"
          % (direct, f05_check, out["f05_is_interp_artifact"]))
    print("median (f=0.5 crossing)=%.6f" % median)
    print("c* (steepest, win %.4f)=%.6f  slope=%.2f" % (PIN_WINDOW, cstar, slope_pin))
    print("slope divergence: coarse(win .01)=%.3f  fine(win 1e-4)=%.3f  growth=x%.2f"
          % (sl_co, sl_fi, growth))
    print("KILL: delta_f=%.6f (thr %.4f) -> %s" % (delta_f, KILL_DELTA_F, out["kill"]["verdict"]))
    print("  f(c*)=%.6f f(c*-0.02)=%.6f f(c*+0.02)=%.6f" % (cf, dlo, dhi))
    print("2p spike mass in [0.5-eps,0.5):", spike)
    print("local profile:", profile)
    print("WROTE", fp)

if __name__ == "__main__":
    main()
