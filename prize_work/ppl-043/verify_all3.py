"""verify_all3.py M0 M1 — exact min-top window method, memory-bounded layout.

IDENTICAL algorithm to verify_all.py (same first-hit-in-ascending-top search,
same K1/D1/D2 even-m case split, same 'how' labels and rep bookkeeping); only
the data access differs: data3/ int32 arrays, bit-packed membership, bounds
computed per m via math.isqrt, inverse index via searchsorted on the
strictly increasing ev/od, memoryview scalar access (no tolist()).

Also records reps for the independent direct-sum verifier: every 200000th m,
1000 random m, and the chunk's boundary extremes. Writes
t_<M0>_<M1>.npy + verify_all3_<M0>_<M1>.json.
"""
import json, sys, os, time, math, heapq, random
import numpy as np

def main():
    t0 = time.time()
    m0 = int(sys.argv[1]); m1 = int(sys.argv[2])
    assert 1 <= m0 <= m1
    D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data3')
    meta = json.load(open(D + '/meta3.json'))
    Smax = meta['Smax']; SmaxE = meta['SmaxE']; SmaxO = meta['SmaxO']
    N = meta['N']
    print(f'loading data3 (Smax={Smax}, SmaxE={SmaxE}, SmaxO={SmaxO}, N={N}) ...',
          flush=True)
    pr = np.load(D + '/pr.npy')
    s  = np.load(D + '/s.npy')
    ev = np.load(D + '/ev.npy')
    od = np.load(D + '/od.npy')
    memE = bytes(np.load(D + '/memE_bits.npy'))
    memO = bytes(np.load(D + '/memO_bits.npy'))
    pr_mv = memoryview(pr); s_mv = memoryview(s)
    t_load = time.time() - t0
    print(f'loaded in {t_load:.1f}s; processing m in [{m0}, {m1}]', flush=True)

    t_arr = np.zeros(m1 - m0 + 1, np.int32)   # min top; 0 = none found
    violations = []
    gap_heap = []
    n_tie = 0
    max_top = 0; max_top_m = 0
    rnd = random.Random(42)
    rnd_set = set(rnd.randint(m0, m1) for _ in range(1000))
    force = {m0, m0 + 1, m1, m1 - 1, (m0 + m1) // 2}
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
        while pr_mv[p0] < mm: p0 += 1
        b = (2 * mm + math.isqrt(121 * mm) // 5) if not (mm & 1) \
            else (mm + math.isqrt(529 * mm) // 5)
        if mm & 1:
            i = p0; hit = None
            while i < N and pr_mv[i] <= b:
                si = s_mv[i]
                if i & 1:
                    sj = mm - si
                    if 1 <= sj <= SmaxE and memE[sj >> 3] & (1 << (sj & 7)):
                        j = int(np.searchsorted(ev, np.int32(sj))) * 2
                        if j < i: hit = (i, j, 'sum'); break
                else:
                    if si == mm: hit = (i, 0, 'k1'); break
                    sj = mm - si
                    if 2 <= sj <= SmaxO and memO[sj >> 3] & (1 << (sj & 7)):
                        j = int(np.searchsorted(od, np.int32(sj))) * 2 + 1
                        if j < i: hit = (i, j, 'sum'); break
                i += 1
            if hit is not None:
                i, j, how = hit
                note(mm, i, j, how, pr_mv[i], b)
        else:
            bt = None
            if mm <= SmaxO and memO[mm >> 3] & (1 << (mm & 7)):
                i = int(np.searchsorted(od, np.int32(mm))) * 2 + 1
                bt = (pr_mv[i], i, 0, 'k1')
            c = 1
            while True:
                v = mm + ev[c]
                if v > SmaxE: break
                if memE[v >> 3] & (1 << (v & 7)):
                    a = int(np.searchsorted(ev, np.int32(v)))
                    i = a * 2; tt = pr_mv[i]
                    if bt is None or tt < bt[0]: bt = (tt, i, 2 * c, 'diff')
                    break
                if bt is not None and v >= bt[0]: break
                c += 1
            d = 0
            while True:
                v = mm + od[d]
                if v > SmaxO: break
                if memO[v >> 3] & (1 << (v & 7)):
                    bpos = int(np.searchsorted(od, np.int32(v)))
                    i = bpos * 2 + 1; tt = pr_mv[i]
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
            print(f'  m={mm}  elapsed={el:.0f}s  rate={rate/1000:.1f}k/s  eta={eta:.0f}s',
                  flush=True)

    el = time.time() - t_start
    # any m with no found rep (t == 0) is a missed search => violation
    zero_idx = np.nonzero(t_arr == 0)[0]
    for z in zero_idx:
        violations.append(int(z) + m0)
    gap_heap.sort()
    out = {
        'm0': m0, 'm1': m1, 'n': m1 - m0 + 1,
        'time_s': round(el, 1), 'load_s': round(t_load, 1),
        'violations': sorted(set(violations))[:200],
        'n_viol': len(set(violations)),
        'n_zero_t': int(zero_idx.size),
        'n_tie_bound': n_tie,
        'worst_gap20': [{'m': m, 't': t, 'b': b, 'gap': g}
                        for (g, m, t, b) in gap_heap],
        'max_top': {'m': max_top_m, 't': max_top},
        'reps': {str(k): v for k, v in reps.items()},
        't_file': f't_{m0}_{m1}.npy',
    }
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.environ.get('PPL043_OUT', here)   # redirect outputs (gate runs)
    os.makedirs(outdir, exist_ok=True)
    np.save(os.path.join(outdir, out['t_file']), t_arr)
    fn = f'verify_all3_{m0}_{m1}.json'
    json.dump(out, open(os.path.join(outdir, fn), 'w'), indent=1)
    print(f'done: {out["n"]} m in {el:.0f}s; violations={out["n_viol"]}; '
          f'ties={n_tie}; wrote {fn}', flush=True)

if __name__ == '__main__':
    main()
