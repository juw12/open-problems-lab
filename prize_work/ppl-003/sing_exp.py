# -*- coding: utf-8 -*-
"""
Erdos #50 (PPL 003) - local-singularity-strength (Hölder/local-dimension) + mass
spectrum of the FULL interior singular set of f(c)=dens{n: phi(n)/n < c}.

New sub-target (beyond job-35's 3008-point MAP): for each STRONG divergent-slope
interior point, measure the one-sided-below windowed-slope divergence EXPOENT
gamma = d log s_-(w) / d log w  with  s_-(c;w) = (F(c)-F(c-w))/w  over fine windows
w in [5e-5, 5e-3] (above the N*eps=10 finite-N floor ~1e-6, so asymptotic regime).
Define the LOCAL MASS (local-dimension) exponent p_c = 1 + gamma_c, since the
one-sided-below pile-up mass M_-(c;w)=F(c)-F(c-w) ~ w^p_c as w->0+. A smooth point
has p=1; a singular accumulation point has p<1 (sublinear, mass concentrates).
Test MONOFRACtal (p universal ~ const across all strong points) vs MULTIFRACTAL
(p stratified by prime-set size |S|). Also compute the mass spectrum (rank strong
points by pile-up mass M_-(c;w0), Pareto exponent, cumulative share, top-K rationals)
to answer the card's core question: is c=1/2 isolated or part of a structured set?

Borrowed field: physics (fractal / local Hölder dimension / multifractal spectrum +
Pareto power-law). This is the w->0 axis, distinct from job-35's two-sided GROWTH
metric (a different quantity) and from job-35-idea-2 (the N->inf axis, already in bank).

Gates (must pass before the result):
  (a) one-sided-below M_-(1/2;5e-5) == 0.072339  (job 8, N=1e7)
  (b) c=1/2 two-sided windowed-slope growth == 58.46  (job 27/35, w 0.05->0.00025)
  (c) c=1/2 one-sided-below slope exponent gamma ~= -0.9  (job 24 two-sided alpha~0.93)
Independent verifier (NOT the same code run twice):
  (i)  repeat the per-point exponent at the INDEPENDENT second N=5e6 (N-robust, above
       its ~2e-6 floor) -> gamma_c(1e7) vs gamma_c(5e6) must agree within tolerance;
  (ii) analytic 2^a*p tower anchors the c=1/2 finest-regime mass (job 22/24).

Hidden-number-safe output: print only booleans, class names, short exact rationals,
and floats rounded to <=4 significant figures. Never print raw integer counts >=1e8.
"""
import json, math, sys
import numpy as np
from fractions import Fraction

N1 = 10_000_000
N2 = 5_000_000
STRONG = 20.0          # growth threshold for "strong" points (job-35 metric)
W_FINE = np.geomspace(5e-5, 5e-3, 12)   # kill-tested primary exponent regime
W_WIDE = np.geomspace(5e-5, 1e-2, 12)   # robustness cross-check (200x shrink)
W0   = 5e-3            # pile-up mass reference window for the mass spectrum
J35  = "scan_N10000000.json"

def build(N, path):
    phi = np.load(path, mmap_mode='r').astype(np.int64)
    n = np.arange(1, len(phi), dtype=np.int64)
    r = (phi[1:].astype(np.float64) / n.astype(np.float64)).astype(np.float64)
    rs = np.sort(r)
    return r, rs

def F(rs, c):
    return float(np.searchsorted(rs, c, 'left')) / len(rs)

def Fvec(rs, arr):
    return np.searchsorted(rs, arr, 'left') / len(rs)

def slope_growth(rs, c):
    # two-sided windowed slope, job-35 windows 0.05 -> 0.00025 (11 geomspace)
    W = np.geomspace(0.05, 0.00025, 11)
    def s(w):
        return (F(rs, c + w) - F(rs, c - w)) / (2.0 * w)
    return s(W[-1]) / s(W[0])

def onerow(rs, c, W):
    # one-sided-below slope s_-(c;w)=(F(c)-F(c-w))/w and mass M_-(c;w)=F(c)-F(c-w)
    cw = c - W
    Mcw = F(rs, c) - Fvec(rs, cw)        # mass in [c-w, c)  (below c)
    Mcw = np.clip(Mcw, 0.0, None)
    s = Mcw / W
    return Mcw, s

def fit_gamma(W, s):
    m = s > 0
    if m.sum() < 4:
        return None, None
    x = np.log(W[m]); y = np.log(s[m])
    A = np.polyfit(x, y, 1)
    yhat = np.polyval(A, x)
    ss_res = float(np.sum((y - yhat) ** 2)); ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None
    return float(A[0]), r2

def main():
    out = {"job": "lab_erdos_261111_38", "N": [N1, N2], "strong_thresh": STRONG,
           "W_FINE": [float(x) for x in W_FINE], "W_WIDE": [float(x) for x in W_WIDE],
           "w0_mass": W0}

    r1, rs1 = build(N1, "phi_cache_%d.npy" % N1)
    r2, rs2 = build(N2, "phi_cache_%d.npy" % N2)
    out["n_cache"] = [int(len(r1)), int(len(r2))]

    # ---- phi-cache integrity gates ----
    mean_r = float(np.mean(r1)); mean_expect = 6.0 / (math.pi ** 2)
    out["gate_mean_phiN"] = {"value": round(mean_r, 6), "expect": round(mean_expect, 6),
                              "absdiff": round(abs(mean_r - mean_expect), 6),
                              "pass": abs(mean_r - mean_expect) < 5e-4}

    j35 = json.load(open(J35))
    dlist = j35["divergent_list"]
    out["n_divergent_job35"] = len(dlist)

    # ---- gate (a): one-sided-below mass at c=1/2, N=1e7 ----
    M12_fine = F(rs1, 0.5) - F(rs1, 0.5 - 5e-5)
    out["gate_a_mass12"] = {"value": round(float(M12_fine), 6), "expect": 0.072339,
                             "absdiff": round(abs(float(M12_fine) - 0.072339), 6),
                             "pass": abs(float(M12_fine) - 0.072339) < 3e-4}
    # ---- gate (b): c=1/2 two-sided growth (job 27/35) ----
    g12 = slope_growth(rs1, 0.5)
    out["gate_b_growth12"] = {"value": round(g12, 3), "expect": 58.46,
                               "absdiff": round(abs(g12 - 58.46), 3),
                               "pass": abs(g12 - 58.46) < 1.0}
    # ---- gate (c): c=1/2 one-sided-below slope exponent ~ -0.9 (job 24) ----
    Mcw, s = onerow(rs1, 0.5, W_FINE)
    gam12, r2_12 = fit_gamma(W_FINE, s)
    out["gate_c_gamma12"] = {"value": round(gam12, 4) if gam12 is not None else None,
                              "r2": round(r2_12, 4) if r2_12 is not None else None,
                              "expect_range": [-1.2, -0.6],
                              "pass": gam12 is not None and -1.2 <= gam12 <= -0.6 and (r2_12 or 0) > 0.95}
    # robustness: wide-window exponent at c=1/2
    Mcw_w, s_w = onerow(rs1, 0.5, W_WIDE)
    gam12_w, r2_12_w = fit_gamma(W_WIDE, s_w)
    out["gate_c_gamma12_wide"] = {"value": round(gam12_w, 4) if gam12_w is not None else None,
                                   "r2": round(r2_12_w, 4) if r2_12_w is not None else None}

    out["all_gates_pass"] = (out["gate_mean_phiN"]["pass"] and out["gate_a_mass12"]["pass"]
                             and out["gate_b_growth12"]["pass"] and out["gate_c_gamma12"]["pass"])

    # ---- analytic 2^a*p tower anchor at c=1/2 (finest regime), job 22/24 ----
    def tower_A(N, w):
        # A(N,eps)=(1/N) sum_{p odd, 1/(2eps)<=p<=N/2} floor(log2(N/p)); eps=w
        import math as _m
        lo = 1.0 / (2.0 * w); hi = N / 2.0
        total = 0
        # iterate odd integers p in [lo,hi]; only primes contribute but the 2^a*p
        # FAMILIES with phi/n=1/2-1/(2p) need p prime. Use job-22 closed form:
        # mass from 2^a*p is sum over odd primes p of (number of a) / N.
        # Approx via counting primes by simple sieve up to hi for w0 regime.
        hi_i = int(hi)
        if hi_i < 3:
            return 0.0
        s = bytearray(b'\x01') * (hi_i + 1); s[0:2] = b'\x00\x00'
        for q in range(3, int(hi_i ** 0.5) + 1, 2):
            if s[q]:
                s[q*q:hi_i+1:2*q] = b'\x00' * len(s[q*q:hi_i+1:2*q])
        A = 0.0
        p = lo
        # sum over odd primes p>=lo, <=hi of floor(log2(N/p))
        # loop over set bits is slow; use np on the boolean sieve
        bs = np.frombuffer(s, dtype=np.uint8).astype(bool)
        primes = np.where(bs)[0]
        primes = primes[(primes >= lo) & (primes <= hi) & (primes % 2 == 1)]
        if len(primes) == 0:
            return 0.0
        terms = np.floor(np.log2(N / primes.astype(np.float64)))
        A = float(np.sum(terms)) / N
        return A
    # NOTE: float64 helper
    try:
        tower12 = tower_A(N1, W_FINE[-1])
        out["tower_anchor12_fine"] = {"tower": round(tower12, 6), "sieve_mass": round(float(M12_fine), 6),
                                       "ratio": round(tower12 / float(M12_fine), 5) if M12_fine > 0 else None}
    except Exception as e:
        out["tower_anchor12_fine"] = {"error": str(e)[:80]}

    # ---- select strong points (growth >= STRONG) from job-35 divergent list ----
    strong = [r for r in dlist if r["growth"] >= STRONG]
    out["n_strong"] = len(strong)

    def process(rs, N, W, mincount_frac):
        rows = []
        for r in strong:
            c = Fraction(r["c"])
            cf = float(c)
            Mcw, s = onerow(rs, cf, W)
            m_fine = float(Mcw[-1])           # mass in finest window [c-Wmin, c)
            gamma, r2 = fit_gamma(W, s)
            p = (1.0 + gamma) if gamma is not None else None
            rows.append({
                "c": r["c"], "class": r["class"], "S": r["S"], "growth": r["growth"],
                "absS": len(r["S"]),
                "M_w0": round(float(F(rs, cf) - F(rs, cf - W0)), 6),   # pile-up mass at w0
                "M_fine": round(m_fine, 6),
                "gamma": round(gamma, 4) if gamma is not None else None,
                "p_local": round(p, 4) if p is not None else None,
                "r2": round(r2, 4) if r2 is not None else None,
                "resolves": bool(m_fine * N >= mincount_frac),
            })
        return rows

    rows1 = process(rs1, N1, W_FINE, 20)
    rows2 = process(rs2, N2, W_FINE, 10)
    out["n_resolves_1e7"] = sum(1 for x in rows1 if x["resolves"])
    out["n_resolves_5e6"] = sum(1 for x in rows2 if x["resolves"])

    # ---- robustness: wide-window exponents (N=1e7) ----
    def wide_gamma(rs, cf):
        Mcw, s = onerow(rs, cf, W_WIDE)
        g, r2 = fit_gamma(W_WIDE, s)
        return (round(g, 4) if g is not None else None), (round(r2, 4) if r2 is not None else None)
    # store for a few headline points only (keep file lean)
    head_names = ["1/2", "1/3", "2/5", "2/3"]
    out["wide_gamma_head"] = {}
    for r in strong:
        if r["c"] in head_names:
            g, r2 = wide_gamma(rs1, float(Fraction(r["c"])))
            out["wide_gamma_head"][r["c"]] = {"gamma": g, "r2": r2}

    # ---- universal vs stratified-by-|S| test on CLEAN points (r2>0.9, resolves) ----
    from collections import defaultdict
    clean = [x for x in rows1 if x["resolves"] and x["r2"] is not None and x["r2"] > 0.9
             and x["p_local"] is not None]
    out["n_clean"] = len(clean)
    by_abs = defaultdict(list)
    for x in clean:
        by_abs[x["absS"]].append(x["p_local"])
    stats = {}
    allp = []
    for a in sorted(by_abs):
        v = np.array(by_abs[a])
        allp.extend(v.tolist())
        stats[str(a)] = {"n": int(len(v)), "p_mean": round(float(v.mean()), 4),
                          "p_sd": round(float(v.std()), 4),
                          "p_min": round(float(v.min()), 4), "p_max": round(float(v.max()), 4)}
    allp = np.array(allp)
    overall = {"n": int(len(allp)), "p_mean": round(float(allp.mean()), 4),
               "p_sd": round(float(allp.std()), 4),
               "p_min": round(float(allp.min()), 4), "p_max": round(float(allp.max()), 4)}
    # stratification signal: range of per-|S| means vs overall sd
    means = np.array([stats[str(a)]["p_mean"] for a in sorted(by_abs)])
    out["p_local_by_absS"] = stats
    out["p_local_overall"] = overall
    out["absS_mean_range"] = round(float(means.max() - means.min()), 4) if len(means) >= 2 else None
    # classification rule: stratified if per-|S| mean spread > 1.5x overall per-point sd
    sd_overall = float(allp.std()) if len(allp) else 0.0
    out["stratified_by_absS"] = bool(len(means) >= 2 and (float(means.max() - means.min()) > 1.5 * sd_overall))

    # ---- independent verifier: N=5e6 gamma agreement (clean points) ----
    map2 = {x["c"]: x for x in rows2}
    agree = 0; tot = 0; diffs = []
    for x in clean:
        y = map2.get(x["c"])
        if y is None or y["gamma"] is None or x["gamma"] is None:
            continue
        if not (y["resolves"] and y["r2"] is not None and y["r2"] > 0.9):
            continue
        tot += 1
        d = abs(x["gamma"] - y["gamma"])
        diffs.append(d)
        if d <= 0.15:
            agree += 1
    out["verifier_N2"] = {
        "n_compared": int(tot), "n_agree_tol0.15": int(agree),
        "agree_frac": round(agree / tot, 4) if tot else None,
        "gamma_diff_mean": round(float(np.mean(diffs)), 4) if diffs else None,
        "gamma_diff_max": round(float(np.max(diffs)), 4) if diffs else None,
        "pass": bool(tot >= 50 and (agree / tot if tot else 0) >= 0.80),
    }

    # ---- mass spectrum (rank strong points by M_w0 pile-up mass, N=1e7) ----
    ms = sorted([x for x in rows1 if x["M_w0"] is not None], key=lambda x: -x["M_w0"])
    masses = np.array([x["M_w0"] for x in ms], dtype=float)
    ranks = np.arange(1, len(masses) + 1)
    mpos = masses > 0
    if mpos.sum() >= 10:
        beta, b = np.polyfit(np.log(ranks[mpos]), np.log(masses[mpos]), 1)
        beta = -beta  # mass ~ rank^-beta
    else:
        beta, b = None, None
    # cumulative share of top fractions
    tot_mass = float(masses.sum())
    def cumshare(k):
        return float(masses[:k].sum()) / tot_mass if tot_mass > 0 else None
    out["mass_spectrum"] = {
        "n_ranked": int(len(ms)),
        "pareto_beta": round(float(beta), 4) if beta is not None else None,
        "top1_share": round(cumshare(1), 4), "top10_share": round(cumshare(min(10, len(ms))), 4),
        "top50_share": round(cumshare(min(50, len(ms))), 4),
        "top10pct_share": round(cumshare(max(1, int(0.1 * len(ms)))), 4),
        "c12_rank": None, "c12_M_w0": None,
        "top12": [(x["c"], x["class"], x["M_w0"], x["p_local"]) for x in ms[:12]],
    }
    for i, x in enumerate(ms):
        if x["c"] == "1/2":
            out["mass_spectrum"]["c12_rank"] = int(i + 1)
            out["mass_spectrum"]["c12_M_w0"] = x["M_w0"]

    # ---- gamma (local-dimension) spectrum headline: strongest local singularities ----
    cg = [x for x in clean if x["p_local"] is not None]
    cg = sorted(cg, key=lambda x: x["p_local"])  # smallest p = strongest singularity
    out["localdim_spectrum"] = {
        "n": len(cg),
        "min_p": round(min(x["p_local"] for x in cg), 4),
        "median_p": round(float(np.median([x["p_local"] for x in cg])), 4),
        "top10_lowest_p": [(x["c"], x["class"], x["p_local"], x["gamma"]) for x in cg[:10]],
        "c12_p_local": next((x["p_local"] for x in cg if x["c"] == "1/2"), None),
    }

    # store full per-point rows (hidden-safe: fractions + rounded floats + short rationals)
    out["rows_N1e7"] = rows1
    out["rows_N5e6"] = rows2

    def _pyo(o):
        if hasattr(o, "item"):
            return o.item()
        return str(o)
    with open("sing_exp_N10000000.json", "w") as f:
        json.dump(out, f, indent=2, default=_pyo)

    # ---- hidden-number-safe console summary ----
    print("=== gates ===")
    print("mean(phi/N) ok:", out["gate_mean_phiN"]["pass"],
          "(", round(out["gate_mean_phiN"]["value"], 5), ")")
    print("a: M_-(1/2;5e-5) ok:", out["gate_a_mass12"]["pass"],
          "(", round(out["gate_a_mass12"]["value"], 5), "expect 0.072339 )")
    print("b: growth(1/2) ok:", out["gate_b_growth12"]["pass"],
          "(", round(out["gate_b_growth12"]["value"], 3), "expect 58.46 )")
    print("c: gamma(1/2) ok:", out["gate_c_gamma12"]["pass"],
          "gamma=", out["gate_c_gamma12"]["value"], "r2=", out["gate_c_gamma12"]["r2"])
    print("c-wide: gamma(1/2) wide =", out["gate_c_gamma12_wide"]["value"])
    print("ALL GATES PASS:", out["all_gates_pass"])
    print("tower anchor(1/2,finest):", out.get("tower_anchor12_fine"))
    print("=== selection ===")
    print("n_strong:", out["n_strong"], "resolves1e7:", out["n_resolves_1e7"],
          "resolves5e6:", out["n_resolves_5e6"])
    print("=== local-dimension / Hölder characterization ===")
    print("n_clean:", out["n_clean"])
    print("p_local by |S|:", json.dumps(out["p_local_by_absS"]))
    print("p_local overall:", json.dumps(out["p_local_overall"]))
    print("stratified_by_absS:", out["stratified_by_absS"])
    print("localdim spectrum min/median p:", out["localdim_spectrum"]["min_p"],
          out["localdim_spectrum"]["median_p"], "c1/2 p=", out["localdim_spectrum"]["c12_p_local"])
    print("=== independent verifier N=5e6 ===")
    print("compared:", out["verifier_N2"]["n_compared"], "agree<=0.15:", out["verifier_N2"]["n_agree_tol0.15"],
          "frac:", out["verifier_N2"]["agree_frac"], "PASS:", out["verifier_N2"]["pass"])
    print("=== mass spectrum ===")
    print("pareto_beta:", out["mass_spectrum"]["pareto_beta"],
          "top1/10/50 share:", out["mass_spectrum"]["top1_share"],
          out["mass_spectrum"]["top10_share"], out["mass_spectrum"]["top50_share"])
    print("c1/2 mass rank:", out["mass_spectrum"]["c12_rank"], "/", out["mass_spectrum"]["n_ranked"],
          "M_w0=", out["mass_spectrum"]["c12_M_w0"])
    print("top5 by mass:", [(t[0], t[1]) for t in out["mass_spectrum"]["top12"][:5]])
    print("top5 by local dim (smallest p):", [(t[0], t[1]) for t in out["localdim_spectrum"]["top10_lowest_p"][:5]])
    print("WIDE gamma head:", json.dumps(out.get("wide_gamma_head", {})))
    print("WROTE sing_exp_N10000000.json")

if __name__ == "__main__":
    main()
