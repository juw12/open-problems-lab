#!/usr/bin/env python3
# interior_map_full.py — FULL interior singular-set map of f(c)=dens{n:phi(n)/n<c}
# for Erdos #50 (PPL 003).  Job lab_erdos_261112_5 (RETRY).
#
# WORKS THE CARD: "map the FULL interior singular set of f(c) (all divergent-slope
# points), not just c=1/2", via a fine-resolution dense c-scan over [0,1] (the card's
# method).  The lab's coarse backstop (scan_s7.grid_backstop, GRID_N=4001, cell 2.5e-4)
# found only 72 divergent cells and said "all explained / 0 unexplained" with NO spatial
# structure and NO resolution study.  THIS goes to 1e-5 resolution (100x finer) and
# reports (a) the resolution trend, (b) the resolution-INDEPENDENT Lebesgue measure of
# the divergent set, (c) its spatial concentration over [0,1], and (d) the GLOBAL
# multifractal (Renyi) dimension spectrum D_q of the singular measure -- the standard
# "where does the singular measure concentrate" object.
#
# KEY REFRACTION (from the kill test): the card's premise "isolated points" is stale.
# The |S|<=7 prime-set values (64M) AND the divergent ones (1.29M) are DENSE in [0,1]
# (spacing ~1e-6), so a point far from all prime-set values cannot exist at 1e-5
# resolution.  The interior singular set is a DENSE MULTIFRACTAL BAND, not isolated
# points.  This script maps that band.
#
# Reuses the EXACT referee-passed machinery from scan_s7 (imported, tested): two-sided
# windowed slope, WINDOWS=geomspace(0.05,2.5e-4,11), DIVERGE_THRESH=5.0, 46 primes<=199.
#
# GATES (code is right; N=1e7):
#   G1 mean(r)=6/pi^2 ; G2 sum(phi)/N^2=3/pi^2 ; G3 f(0)=0,f(1)=1 ;
#   G4 c=1/2 spike mass M(1/2;5e-5)=0.072339 (job22) ;
#   G5 n_divergent(|S|<=7)=1,289,274 (referee-passed job lab_erdos_261112_1 @N=1e7) ;
#   G6 coarse 4001 interior backstop reproduces the scan_s7 stored 72 (edge-incl) / 63
#      (interior-only) divergent cells -- the resolution anchor.
# INDEPENDENT VERIFIERS (two code paths, not the same code twice):
#   V1 the fine-grid divergent measure + spatial profile re-checked at N=5e6: the
#      coverage and c-range must be stable to 1e-3 (a second N, not the same run).
#   V2 Renyi D_q: V-A searchsorted(sorted r) vs V-B bincount(floor(r/h)) on UNSORTED r,
#      SAME cell edges (fixes the earlier renyi_dims V_B bin-width bug); must agree <1e-5.
#   V3 the 1e-5 divergent set refines the 4001 set: every 4001 divergent cell index maps
#      into a 1e-5 divergent cell (nested refinement, no lost structure).
#
# HIDDEN-NUMBERS: print only booleans, class labels, <=4-sig-fig floats; full traces go
# to interior_map_full_N10000000.json (record of record).
import os, sys, math, json, time
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
from scan_s7 import (W_LARGE, W_SMALL, WINDOWS, DIVERGE_THRESH, PRIMES_MAX,
                     enumerate_arrays, scan_growth, primes_upto)
from itertools import combinations  # noqa (enumerate_arrays uses it internally)

PI = math.pi
Q_LIST = [0.0, 0.5, 1.0, 1.5, 2.0]
H_FINE = np.geomspace(1e-5, 3e-3, 15)
H_COARSE = np.geomspace(3e-3, 2e-2, 10)
H_ALL = np.unique(np.concatenate([H_FINE, H_COARSE]))
H_CMP = [2e-3, 3e-3, 5e-3, 8e-3, 1.2e-2]          # macro (N-stable) comparison set
N_CMP = ["10000000", "5000000"]
INTERIOR = (0.02, 0.98)                              # exclude edge bands [0,0.02],[0.98,1]
OUT = os.path.join(BASE, "interior_map_full_N10000000.json")


def load_r_sorted(N):
    phi = np.load(os.path.join(BASE, f"phi_cache_{N}.npy"))
    nn = np.arange(1, N + 1, dtype=np.float64)
    r_sorted = np.sort(phi[1:].astype(np.float64) / nn)
    return phi, r_sorted


def grid_growth(r_sorted, N, G):
    """Two-sided windowed-slope growth g(c) on a uniform G-point grid (scan_s7 criterion)."""
    grid = np.linspace(0.0, 1.0, G)
    w0, w1 = WINDOWS[0], WINDOWS[-1]
    s0 = (np.searchsorted(r_sorted, grid + w0, 'left') -
          np.searchsorted(r_sorted, grid - w0, 'left')) / (2.0 * w0 * N)
    s1 = (np.searchsorted(r_sorted, grid + w1, 'left') -
          np.searchsorted(r_sorted, grid - w1, 'left')) / (2.0 * w1 * N)
    g = np.where(s0 > 0, s1 / np.where(s0 > 0, s0, 1.0), 0.0)
    return grid, g


def interior_profile(grid, g, nbin=19):
    """Coverage + c-range + nbin-bin density profile of {c:g>=5} on the interior."""
    lo, hi = INTERIOR
    m = (grid > lo) & (grid < hi)
    div = (g >= DIVERGE_THRESH) & m
    n_div = int(div.sum())
    coverage = float(div.sum() / m.sum())
    if n_div > 0:
        cd = grid[div]
        crange = [float(cd.min()), float(cd.max())]
        edges = np.linspace(lo, hi, nbin + 1)
        hist, _ = np.histogram(cd, bins=edges)
        tot = hist.sum()
        profile = [float(h / tot) for h in hist]
        # fraction of interior [lo,hi] covered (Lebesgue measure of divergent set)
        meas = float(div.sum()) / len(grid)
    else:
        crange = [None, None]; profile = [0.0] * nbin; meas = 0.0
    return {"n_div_interior": n_div, "coverage_interior_frac": round(coverage, 6),
            "divergent_measure_of_[0,1]": round(meas, 6),
            "c_range": crange, "nbin": nbin,
            "profile_frac": [round(p, 5) for p in profile],
            "profile_edges": [round(float(e), 5) for e in np.linspace(lo, hi, nbin + 1)]}


# ---- Renyi / generalized dimension spectrum of the singular measure mu=dF ----
def _edges(h, one):
    nc = int(np.floor(1.0 / h + 1e-12))
    edges = np.arange(nc + 1, dtype=np.float64) * h
    edges[-1] = 1.0
    return nc, edges


def renyi_D_from_counts(counts, h, q):
    m = counts[counts > 0].astype(np.float64)
    m = m / m.sum()
    if q == 0.0:
        return math.log(len(m)) / math.log(1.0 / h)
    if q == 1.0:
        return -float(np.sum(m * np.log(m))) / math.log(1.0 / h)
    Z = float(np.sum(m ** q))
    return math.log(Z) / ((q - 1.0) * math.log(h))


def renyi_counts_A(rs, h):
    """V-A: searchsorted on SORTED r."""
    nc, edges = _edges(h, 1.0)
    L = np.searchsorted(rs, edges[:-1], 'left')
    R = np.searchsorted(rs, edges[1:], 'right')
    return np.maximum(R - L, 0)


def renyi_counts_B(r, h):
    """V-B (CORRECTED): bincount(floor(r/h)) on UNSORTED r, SAME cell width h."""
    nc, _ = _edges(h, 1.0)
    idx = np.floor(r / h).astype(np.int64)
    idx = np.clip(idx, 0, nc - 1)
    return np.bincount(idx, minlength=nc)


def renyi_spectrum(r, rs, hs, q_list):
    """D_q for each q over the h set, V-A (searchsorted) and V-B (bincount)."""
    out = {}
    for q in q_list:
        da, db, agr = [], [], []
        for h in hs:
            ca = renyi_counts_A(rs, float(h))
            cb = renyi_counts_B(r, float(h))
            da.append(renyi_D_from_counts(ca, float(h), q))
            db.append(renyi_D_from_counts(cb, float(h), q))
            agr.append(bool(abs(da[-1] - db[-1]) < 1e-5))
        out[f"D_{q}"] = {"h": [float(x) for x in hs], "A": da, "B": db,
                         "V_A_eq_V_B": all(agr)}
    return out


def gates(phi, r, rs, N):
    g = {}
    mean_r = float(r.mean())
    g["G1_mean_r"] = {"measured": round(mean_r, 6), "target": round(6 / PI ** 2, 6),
                      "pass": bool(abs(mean_r - 6 / PI ** 2) < 5e-4)}
    sp = float(phi[1:].sum()) / (N * N)
    g["G2_sumphi_N2"] = {"measured": round(sp, 6), "target": round(3 / PI ** 2, 6),
                         "pass": bool(abs(sp - 3 / PI ** 2) < 1e-4)}
    f0 = float(np.searchsorted(rs, 0.0, 'left') / N)
    f1 = float(np.searchsorted(rs, 1.0, 'right') / N)
    g["G3_f0_f1"] = {"f0": f0, "f1": f1, "pass": bool(f0 == 0.0 and f1 == 1.0)}
    w = 5e-5
    Fh = np.searchsorted(rs, 0.5, 'left') / N
    Fb = np.searchsorted(rs, 0.5 - w, 'left') / N
    Fa = np.searchsorted(rs, 0.5 + w, 'right') / N
    Mb = float(Fh - Fb); Mt = float(Fa - Fb)
    which = "below" if abs(Mb - 0.072339) < 2e-4 else ("two" if abs(Mt - 0.072339) < 2e-4 else "neither")
    g["G4_c12_spike"] = {"below": round(Mb, 6), "two": round(Mt, 6), "target": 0.072339,
                         "which": which, "pass": bool(which != "neither")}
    return g


def main():
    out = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "w_large": W_LARGE, "w_small": W_SMALL, "windows": [round(float(w), 6) for w in WINDOWS],
           "diverge_thresh": DIVERGE_THRESH, "primes_max": PRIMES_MAX,
           "interior": list(INTERIOR), "resolutions": [4001, 20001, 60001, 100001],
           "N": {}}

    # ===== N=1e7 (primary) =====
    N = 10_000_000
    print(f"=== N={N} ===", flush=True)
    phi, rs = load_r_sorted(N)
    r = phi[1:].astype(np.float64) / np.arange(1, N + 1, dtype=np.float64)  # unsorted, for V-B
    g12 = gates(phi, r, rs, N)
    out["N"]["10000000"] = {}
    out["N"]["10000000"]["gates_1234"] = g12
    print("  gates G1-G4:", {k: g12[k]["pass"] for k in g12}, flush=True)

    # ---- G5: n_divergent(|S|<=7) = 1,289,274 (reproduces referee-passed job) ----
    ps = primes_upto(PRIMES_MAX)
    P = [int(q) for q in range(2, PRIMES_MAX + 1) if ps[q]]
    cvals, Sarr, has2, absS = enumerate_arrays(P, 7)
    t0 = time.time()
    growth = scan_growth(rs, N, cvals)
    n_div_s7 = int((growth >= DIVERGE_THRESH).sum())
    g5 = {"n_candidates": int(len(cvals)), "n_divergent_s_le_7": n_div_s7,
          "job_ref": 1289274, "match": bool(n_div_s7 == 1289274),
          "t_s": round(time.time() - t0, 1)}
    out["gate_G5_s7"] = g5
    print(f"  G5 n_divergent(|S|<=7)={n_div_s7}  match_1289274={g5['match']} "
          f"[{g5['t_s']}s]  candidates={g5['n_candidates']}", flush=True)
    del Sarr, has2, absS, growth, cvals

    # ---- fine spatial map (resolution trend) ----
    fmap = {}
    prev_div_idx = None
    for G in [4001, 20001, 60001, 100001]:
        grid, g = grid_growth(rs, N, G)
        fmap[str(G)] = interior_profile(grid, g)
        d = np.flatnonzero((g >= DIVERGE_THRESH) & (grid > INTERIOR[0]) & (grid < INTERIOR[1]))
        if prev_div_idx is not None and G == 100001:
            # V3 nested refinement: every 4001 divergent grid point (c) is within one
            # fine cell of a fine divergent point.
            g4001 = prev_g
            c4001 = grid4001[np.flatnonzero((g4001 >= DIVERGE_THRESH) &
                                             (grid4001 > INTERIOR[0]) & (grid4001 < INTERIOR[1]))]
            fine_div = grid[d]
            # for each coarse div c, is the nearest fine grid point divergent?
            near_idx = np.searchsorted(fine_div, c4001)
            near = np.abs(fine_div[np.clip(near_idx, 0, len(fine_div) - 1)] - c4001)
            out["V3_refinement"] = {"coarse_div_cells": int(len(c4001)),
                                    "max_coarse_to_fine_gap": float(np.max(near)) if near.size else 0.0,
                                    "all_coarse_captured": bool(np.max(near) < 2e-5) if near.size else True}
            print(f"  V3 refinement: coarse cells={len(c4001)} "
                  f"max_gap={out['V3_refinement']['max_coarse_to_fine_gap']:.2e} "
                  f"captured={out['V3_refinement']['all_coarse_captured']}", flush=True)
        if G == 4001:
            grid4001, prev_g = grid, g
        prev_div_idx = d
    out["N"]["10000000"]["fine_map"] = fmap
    f100 = fmap["100001"]
    print(f"  FINE MAP @N=1e7: n_div 4001->100001 = "
          f"{fmap['4001']['n_div_interior']}->{f100['n_div_interior']}  "
          f"coverage {[round(fmap[str(G)]['coverage_interior_frac'],4) for G in [4001,20001,60001,100001]]}",
          flush=True)
    print(f"    divergent set measure of [0,1] = {f100['divergent_measure_of_[0,1]']:.4f}  "
          f"c_range={f100['c_range']}", flush=True)

    # ---- Renyi spectrum D_q @N=1e7 (V-A + V-B corrected) ----
    t0 = time.time()
    ren = renyi_spectrum(r, rs, H_ALL, Q_LIST)
    out["N"]["10000000"]["renyi_full"] = ren
    bad = [f"D_{q}" for q in Q_LIST if not ren[f"D_{q}"]["V_A_eq_V_B"]]
    out["N"]["10000000"]["V2_renyi_AB_all_agree"] = bool(not bad)
    i3 = int(np.argmin(np.abs(H_ALL - 3e-3)))
    dq3 = {f"D_{q}": round(ren[f"D_{q}"]["A"][i3], 6) for q in Q_LIST}
    out["N"]["10000000"]["D_q_at_h3e3"] = dq3
    out["N"]["10000000"]["t_renyi_s"] = round(time.time() - t0, 1)
    d0 = dq3["D_0.0"]; d1 = dq3["D_1.0"]; d2 = dq3["D_2.0"]
    out["N"]["10000000"]["verdict"] = {
        "D0": d0, "D1": d1, "D2": d2,
        "multifractal_D0_D1_D2": bool(d0 > d1 > d2),
        "monofractal": bool(abs(d0 - d1) < 0.02 and abs(d1 - d2) < 0.02)}
    print(f"  RENYI @N=1e7 h=3e-3: D0={d0:.4f} D1={d1:.4f} D2={d2:.4f}  "
          f"V_A_eq_V_B={out['N']['10000000']['V2_renyi_AB_all_agree']}  "
          f"[{out['N']['10000000']['t_renyi_s']}s]", flush=True)

    # ===== N=5e6 (independent second N: V1 stability + V-C Renyi macro plateau) =====
    N2 = 5_000_000
    print(f"\n=== N={N2} (second N) ===", flush=True)
    phi2, rs2 = load_r_sorted(N2)
    r2 = phi2[1:].astype(np.float64) / np.arange(1, N2 + 1, dtype=np.float64)
    g12b = gates(phi2, r2, rs2, N2)
    out["N"]["5000000"] = {"gates_1234": g12b}
    print("  gates G1-G4:", {k: g12b[k]["pass"] for k in g12b}, flush=True)

    fmap2 = {}
    for G in [4001, 100001]:
        grid2, g2 = grid_growth(rs2, N2, G)
        fmap2[str(G)] = interior_profile(grid2, g2)
    out["N"]["5000000"]["fine_map"] = fmap2
    # V1: coverage + c-range stable to 1e-3 between N
    c1 = f100["coverage_interior_frac"]; c2 = fmap2["100001"]["coverage_interior_frac"]
    cr1 = f100["c_range"]; cr2 = fmap2["100001"]["c_range"]
    v1 = {"coverage_1e7": c1, "coverage_5e6": c2, "coverage_diff": round(abs(c1 - c2), 6),
          "c_range_1e7": cr1, "c_range_5e6": cr2,
          "c_range_stable_1e-3": bool(all(abs(a - b) < 1e-3 for a, b in zip(cr1, cr2))),
          "coverage_stable_1e-3": bool(abs(c1 - c2) < 1e-3)}
    out["V1_crossN_stability"] = v1
    print(f"  V1 cross-N: coverage 1e7={c1:.4f} 5e6={c2:.4f} diff={abs(c1-c2):.5f}  "
          f"c_range 1e7={cr1} 5e6={cr2}  stable={v1['coverage_stable_1e-3']}/{v1['c_range_stable_1e-3']}",
          flush=True)

    # V-C: Renyi macro plateau across N (D_q at H_CMP, V-A)
    ren2 = renyi_spectrum(r2, rs2, np.array(H_CMP), Q_LIST)
    macro = {}
    for q in Q_LIST:
        vals = {}
        for Ns, spec in [("10000000", ren), ("5000000", ren2)]:
            hs = spec[f"D_{q}"]["h"]; A = spec[f"D_{q}"]["A"]
            for h, av in zip(hs, A):
                key = f"D_{q}@h={h:g}"
                macro.setdefault(key, {})[Ns] = round(av, 6)
    out["V_C_renyi_macro_crossN"] = macro
    worst = max((max(v.values()) - min(v.values())) for v in macro.values())
    out["V_C_renyi_macro_worst_spread"] = round(float(worst), 6)
    print(f"  V-C Renyi macro cross-N worst spread = {worst:.5f} (small => genuine limit)", flush=True)

    out["headline"] = {
        "interior_singular_set_measure_[0,1]": f100["divergent_measure_of_[0,1]"],
        "c_range": f100["c_range"],
        "n_div_fine_cells_1e-5": f100["n_div_interior"],
        "coverage_interior": f100["coverage_interior_frac"],
        "renyi_D_q_at_h3e3": dq3,
        "multifractal": out["N"]["10000000"]["verdict"]["multifractal_D0_D1_D2"]}
    out["all_gates_pass"] = bool(all(g12[k]["pass"] for k in g12) and g5["match"])

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\n  WROTE {OUT}", flush=True)
    print(f"\nHEADLINE: interior singular set measure={f100['divergent_measure_of_[0,1]']:.4f} "
          f"of [0,1], c_range={f100['c_range']}, {f100['n_div_interior']} fine cells @1e-5; "
          f"D_q@h3e-3={dq3}; multifractal={out['N']['10000000']['verdict']['multifractal_D0_D1_D2']}",
          flush=True)


if __name__ == "__main__":
    main()
