#!/usr/bin/env python3
"""Job lab_erdos_261111_3: push exact f(n) for Erdos #710 (PPL 129) past n=20000.

The card's literal goal (push to 10000/15000/20000) was completed and Dinic-verified by
lab_erdos_..._8 (f(10000)=20626, f(15000)=31576, f(20000)=42561). This job builds ON that
frontier and extends the verified exact range to n=50000, adding the three anchors
30000 / 40000 / 50000, so the Turturean-constant and EP-lower-bound convergence tests can
be re-run on a wider pre-asymptotic range, and so a curvature (second-difference) test of
the ratio r_t(n)=f/(n*sqrt(log n)) -- which a flat level test cannot detect -- can be run.

PRIMARY           : proven O(E)-divisor-adjacency Hopcroft-Karp  (fit_loglog.f_exact).
INDEPENDENT VERIFY: Dinic max-flow  (fit_loglog.dinic_maxflow), a different algorithm
                    family. For each anchor N the exact value F is certified by
                    dinic(N, F)==True  AND  dinic(N, F-1)==False  (the minimality condition).

GATE (fast, independent, run before any new anchor):
    naive HK reproduces f(1)=2,f(2)=3,f(3)=4,f(4)=6,f(10)=15,f(20)=29,f(30)=46
    optimized HK reproduces f(100)=161,f(500)=878,f(1000)=1817,f(2000)=3815,f(5000)=10081
    (last two re-confirm the prior job-36 / job-25 anchors).

Resumable: push50k_data.json flushed after each anchor. Usage:
    python push50k.py             # all anchors not yet done
    python push50k.py 30000       # just this anchor (kill test)
"""
import json, os, math, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_loglog as fl

NMAX = 50000
fl.SPF = fl.build_spf(4 * NMAX)        # divisors_le factors m up to ~4n (n=50000 -> m<200000)

D = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(D, "push50k_data.json")
ANCHORS = [30000, 40000, 50000]
SQRT2_E = math.sqrt(2) / math.e
TWO_SQRT_E = 2 / math.sqrt(math.e)

PUBLISHED = {1: 2, 2: 3, 3: 4, 4: 6, 10: 15, 20: 29, 30: 46,
             100: 161, 500: 878, 1000: 1817, 2000: 3815, 5000: 10081}


def rt(n, f):    # (a) Turturean-style ratio  -> target sqrt(2)/e = 0.52026
    return f / (n * math.sqrt(math.log(n)))


def ep_lower(n):  # EP 1980 lower bound  (2/sqrt(e)) * n * sqrt(log n / loglog n)
    return TWO_SQRT_E * n * math.sqrt(math.log(n) / math.log(math.log(n)))


def re_card(n, f):  # (b) card's literal formula  f*loglog n / ((2/sqrt(e)) n sqrt(log n))
    return f * math.log(math.log(n)) / (TWO_SQRT_E * n * math.sqrt(math.log(n)))


def r_prop(n, f):  # (c) the definition that tends to 1 iff the EP bound is asymptotically tight
    return f / ep_lower(n)


def verify(n, f):
    okF = fl.dinic_maxflow(n, f)         # must be True
    okFm1 = fl.dinic_maxflow(n, f - 1)   # must be False
    return bool(okF), bool(okFm1)


def gate():
    small_ok = all(fl.naive_f_exact(n) == v for n, v in
                   [(1, 2), (2, 3), (3, 4), (4, 6), (10, 15), (20, 29), (30, 46)])
    big_ok = all(fl.f_exact(n) == PUBLISHED[n] for n in (100, 500, 1000, 2000, 5000))
    detail = {n: fl.f_exact(n) for n in (100, 500, 1000, 2000, 5000)}
    ok = small_ok and big_ok
    print(f"GATE: {'PASS' if ok else 'FAIL'}  small={small_ok}  "
          f"big={ {k: (detail[k] == PUBLISHED[k]) for k in detail} }", flush=True)
    return ok


def main():
    data = {}
    if os.path.exists(F):
        data = {int(k): v for k, v in json.load(open(F)).items()}
    targets = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else ANCHORS

    print("GATE (independent) ...", flush=True)
    if not gate():
        print("ABORT: gate failed -- do not trust new anchors", flush=True)
        sys.exit(1)

    for n in targets:
        if n in data and data[n].get("f") is not None and data[n].get("dinic_ok"):
            print(f"n={n} already done f={data[n]['f']}", flush=True)
            continue
        t1 = time.time()
        f = fl.f_exact(n)
        thk = time.time() - t1
        if f is None:
            print(f"n={n} NOT FOUND (hi>200000)", flush=True)
            continue
        tv = time.time()
        okF, okFm1 = verify(n, f)
        tver = time.time() - tv
        dinic_ok = okF and (not okFm1)
        data[n] = dict(f=f, hk_time_s=round(thk, 1), dinic_time_s=round(tver, 1),
                       dinic_ok_at_f=okF, dinic_false_at_fm1=(not okFm1), dinic_ok=dinic_ok,
                       ratio_t=round(rt(n, f), 6),
                       ratio_ep_card=round(re_card(n, f), 6),
                       ep_lower=round(ep_lower(n), 1),
                       r_prop=round(r_prop(n, f), 6))
        with open(F, "w") as fh:
            json.dump({str(k): v for k, v in sorted(data.items())}, fh, indent=2)
        print(f"n={n} f={f} hk={thk:.1f}s dinic={tver:.1f}s dinic_ok={dinic_ok} "
              f"r_t={rt(n, f):.6f} r_ep={re_card(n, f):.6f} r_prop={r_prop(n, f):.6f} "
              f"(ep_lower={ep_lower(n):.0f})", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
