#!/usr/bin/env python3
"""Job lab_erdos_..._8: push exact f(n) for Erdos #710 (PPL 129) to n=10000/15000/20000.

PRIMARY   : proven O(E)-divisor-adjacency Hopcroft-Karp  (fit_loglog.f_exact).
INDEPENDENT VERIFIER : Dinic max-flow  (fit_loglog.dinic_maxflow), a different
algorithm family (max-flow, not matching).  For each anchor N the exact value F is
certified by the boundary check  dinic(N, F) == True  AND  dinic(N, F-1) == False
(matching exists at the found F, not one shorter) -- the exact minimality condition.

Convention (same as job-25/36, erdosproblems (n, n+f(n))):  F = minimal L such that
the window {N+1, ..., N+L-1} holds distinct a_k with k | a_k for 1<=k<=N.

Ratios reported at each anchor:
  (a) Turturean-style   r_t(N) = f / (N sqrt(log N))                       -> sqrt(2)/e = 0.52026
  (b) EP-lb sharpness   r_e(N) = f loglog(N) / ((2/sqrt(e)) N sqrt(log N)) -> 1 (asymptotically tight)

Resumable: push_data.json flushed after each anchor.  Usage:
  python push_frontier.py            # all anchors not yet done
  python push_frontier.py 10000      # just this anchor
"""
import json, os, math, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_loglog as fl

NMAX = 20000
fl.SPF = fl.build_spf(4 * NMAX)          # divisors_le factors m up to ~4n

F = os.path.join(os.path.dirname(os.path.abspath(__file__)), "push_data.json")
ANCHORS = [10000, 15000, 20000]
SQRT2_E = math.sqrt(2) / math.e
TWO_SQRT_E = 2 / math.sqrt(math.e)


def rt(n, f):    # (a) Turturean-style ratio
    return f / (n * math.sqrt(math.log(n)))


def re_(n, f):   # (b) EP lower-bound sharpness ratio
    return f * math.log(math.log(n)) / (TWO_SQRT_E * n * math.sqrt(math.log(n)))


def verify(n, f):
    okF = fl.dinic_maxflow(n, f)         # must be True
    okFm1 = fl.dinic_maxflow(n, f - 1)   # must be False
    return bool(okF), bool(okFm1)


def main():
    data = {}
    if os.path.exists(F):
        data = {int(k): v for k, v in json.load(open(F)).items()}
    targets = [int(sys.argv[1])] if len(sys.argv) > 1 else ANCHORS

    # GATE (fast, independent) -- always run before touching new anchors
    g_ok = all(fl.naive_f_exact(n) == v for n, v in
               [(1, 2), (2, 3), (3, 4), (4, 6), (10, 15), (20, 29), (30, 46)])
    g_ok = g_ok and fl.f_exact(100) == 161 and fl.f_exact(500) == 878 and fl.f_exact(1000) == 1817
    print(f"GATE: {'PASS' if g_ok else 'FAIL'}", flush=True)
    if not g_ok:
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
                       ratio_t=round(rt(n, f), 6), ratio_ep=round(re_(n, f), 6))
        with open(F, "w") as fh:
            json.dump({str(k): v for k, v in sorted(data.items())}, fh, indent=2)
        print(f"n={n} f={f} hk={thk:.1f}s dinic={tver:.1f}s dinic_ok={dinic_ok} "
              f"r_t={rt(n, f):.6f} r_ep={re_(n, f):.6f}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
