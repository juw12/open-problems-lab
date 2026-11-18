"""verify_all.py M0 M1 [gate]

Exact min-top search for EVERY m in [M0, M1] using the prefix-sum window method.

Method (see dossier prize_work/ppl-043.md):
  s_0 = 0, s_i = p_i - s_{i-1}.  Any alternating sum of consecutive primes
  p_k..p_n equals s_n - s_{k-1} (n-(k-1) even) or s_n + s_{k-1} (odd).
  Facts: s_i even <=> (i odd or i=0); s_i odd <=> i even >= 2.
  E = {s at even indices} = {0} U {odd values}; O = {s at odd indices} = {even values >= 2}.
  Every representation of m has top p_n >= m (p_i = s_i + s_{i-1}).
  ODD m: top in [m, bnd(m)]; per candidate i, at most one j (sj = m-si or si-m).
  EVEN m: exhaustive case split -> K1o (m in O, j=0), D1 (si = m+ev[c], c>=1),
          D2 (si = m+od[d], d>=0).  Sum cases with j>=2 are parity-impossible.
First hit in ascending top order = exact minimum top.
"""
import json, sys, time, heapq, random
import numpy as np

def main():
    t0 = time.time()
    m0 = int(sys.argv[1]); m1 = int(sys.argv[2])
    gate = 'gate' in sys.argv
    assert 1 <= m0 <= m1
    D = 'data'
    meta = json.load(open(D + '/meta.json'))
    Smax = meta['Smax']
    print(f'loading data (Smax={Smax}) ...', flush=True)
    pr = np.load(D + '/pr.npy')
    s  = np.load(D + '/s.npy')
    ev = np.load(D + '/ev.npy')
    od = np.load(D + '/od.npy')
    memE = open(D + '/memE.ba', 'rb').read()
    memO = open(D + '/memO.ba', 'rb').read()
    bnd_all = np.load(D + '/bnd_all.npy')
    assert len(memE) == Smax + 1 and len(memO) == Smax + 1
    assert len(bnd_all) >= m1
    pr_l = pr.tolist(); s_l = s.tolist()
    N = len(pr_l)
    invE = np.empty(Smax + 1, np.int32)
    for pos, v in enumerate(ev.tolist()): invE[v] = pos
    invO = np.empty(Smax + 1, np.int32)
    for pos, v in enumerate(od.tolist()): invO[v] = pos
    t_load = time.time() - t0
    print(f'loaded in {t_load:.1f}s; processing m in [{m0}, {m1}]', flush=True)

    t_arr = np.zeros(m1 - m0 + 1, np.int32)   # min top; 0 = violation
    violations = []
    gap_heap = []                              # (gap, m, t, b) smallest gaps
    n_tie = 0
    max_top = 0; max_top_m = 0
    rnd = random.Random(42)
    rnd_set = set(rnd.randint(m0, m1) for _ in range(1000))
    force = {10, 20, 303, 2382}
    reps = {}
    p0 = 1

    def note(mm, i, j, how, tt, b):
        nonlocal n_tie, max_top, max_top_m
        gap = b - tt
        if gap < 0:
            violations.append(mm); return
        t_arr[mm - m0] = tt
        if gap == 0: n_tie += 1
        if len(gap_heap) < 20:
            heapq.heappush(gap_heap, (gap, mm, tt, b))
        else:
            heapq.heapreplace(gap_heap, (gap, mm, tt, b))
        if tt > max_top: max_top, max_top_m = tt, mm
        if mm % 200000 == 0 or mm in rnd_set or mm in force:
            reps[mm] = [i, j, how]

    t_start = time.time()
    for mm in range(m0, m1 + 1):
        while pr_l[p0] < mm: p0 += 1
        b = int(bnd_all[mm - 1])
        if mm & 1:
            i = p0; hit = None
            while i < N and pr_l[i] <= b:
                si = s_l[i]
                if i & 1:
                    sj = mm - si
                    if 1 <= sj <= Smax and memE[sj]:
                        j = int(invE[sj]) * 2
                        if j < i: hit = (i, j, 'sum'); break
                else:
                    if si == mm: hit = (i, 0, 'k1'); break
                    sj = mm - si
                    if 2 <= sj <= Smax and memO[sj]:
                        j = int(invO[sj]) * 2 + 1
                        if j < i: hit = (i, j, 'sum'); break
                i += 1
            if hit is not None:
                i, j, how = hit
                note(mm, i, j, how, pr_l[i], b)
        else:
            bt = None
            if memO[mm]:
                i = int(invO[mm]) * 2 + 1
                bt = (pr_l[i], i, 0, 'k1')
            c = 1
            while True:
                v = mm + ev[c]
                if v > Smax: break
                if memE[v]:
                    a = int(invE[v]); i = a * 2; tt = pr_l[i]
                    if bt is None or tt < bt[0]: bt = (tt, i, 2 * c, 'diff')
                    break
                if bt is not None and v >= bt[0]: break
                c += 1
            d = 0
            while True:
                v = mm + od[d]
                if v > Smax: break
                if memO[v]:
                    bpos = int(invO[v]); i = bpos * 2 + 1; tt = pr_l[i]
                    if bt is None or tt < bt[0]: bt = (tt, i, 2 * d + 1, 'diff')
                    break
                if bt is not None and v >= bt[0]: break
                d += 1
            if bt is None:
                violations.append(mm)
            else:
                tt, i, j, how = bt
                note(mm, i, j, how, tt, b)
        if (mm - m0 + 1) % 200000 == 0:
            el = time.time() - t_start
            rate = (mm - m0 + 1) / el
            eta = (m1 - mm) / rate
            print(f'  m={mm}  elapsed={el:.0f}s  rate={rate/1000:.1f}k/s  eta={eta:.0f}s', flush=True)

    el = time.time() - t_start
    gap_heap.sort()
    out = {
        'm0': m0, 'm1': m1, 'n': m1 - m0 + 1,
        'time_s': round(el, 1), 'load_s': round(t_load, 1),
        'violations': violations[:200], 'n_viol': len(violations),
        'n_tie_bound': n_tie,
        'worst_gap20': [{'m': m, 't': t, 'b': b, 'gap': g} for (g, m, t, b) in gap_heap],
        'max_top': {'m': max_top_m, 't': max_top},
        'reps': {str(k): v for k, v in reps.items()},
        't_file': f't_{m0}_{m1}.npy',
    }
    np.save(f't_{m0}_{m1}.npy', t_arr)
    fn = f'verify_all_{m0}_{m1}.json'
    json.dump(out, open(fn, 'w'), indent=1)
    print(f'done: {out["n"]} m in {el:.0f}s; violations={len(violations)}; ties={n_tie}; wrote {fn}', flush=True)

    if gate:
        base = json.load(open('out_200000.json'))
        ok = True
        # (a) all published samples: exact min top must match
        mism = []
        for (m, t, bnd) in base['samples']:
            if t_arr[m - m0] != t:
                mism.append((m, int(t_arr[m - m0]), t))
            if int(bnd_all[m - 1]) != bnd:
                mism.append((m, 'bnd', int(bnd_all[m - 1]), bnd))
        # (b) slide examples
        ex = {10: 17, 20: 41, 303: 383, 2382: 4871}
        for m, t in ex.items():
            if t_arr[m - m0] != t:
                ok = False; print(f'EXAMPLE FAIL m={m}: got {t_arr[m-m0]} want {t}')
        if mism:
            ok = False; print(f'SAMPLE MISMATCHES ({len(mism)}):', mism[:10])
        if violations:
            ok = False; print(f'VIOLATIONS: {violations[:10]}')
        print('GATE PASS' if ok else 'GATE FAIL', flush=True)
        sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
