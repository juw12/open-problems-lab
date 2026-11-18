"""verify_fb5.py M0 M1 — exact min-top window method, FAST layout, M = 5e8.

SAME algorithm and SAME output schema as verify_fb.py (window method over the
alternating prefix sums; first hit in ascending top order = exact min top):
  t_<M0>_<M1>.npy + verify_all3_<M0>_<M1>.json (schema-compatible; the 'code'
  field records provenance). PPL043_OUT env var redirects outputs (gate runs).

Copied from verify_fb.py; two memory changes for M = 5e8 (only ~10.8 GB free
on the shared 32 GB box, 8 workers):
  - MERGED inverse table: ev values are odd/0, od values are even, so the two
    sets are disjoint by value-parity -> ONE int32 table of size Smax+1 holds
    both position maps (halves the ~8 GB of two separate invE/invO tables).
  - numpy indexing of pr/s directly (no pr_l/s_l Python-list copies), saving
    ~2.4 GB. Membership still via the data5 bit-packed tables.
Peak memory ~6 GB (vs ~11.4 GB for the unmodified 3e8 layout at this scale).
The gate (byte-identical cross-check on 1..200000) proves these access changes
do not alter the computed t-values.
"""
import json, sys, os, time, math, heapq, random
import numpy as np

def build(D):
    """Load data3/ and build the fast access structures (call once)."""
    meta = json.load(open(os.path.join(D, 'meta5.json')))
    SmaxE = meta['SmaxE']; SmaxO = meta['SmaxO']
    print(f'loading data5 (SmaxE={SmaxE}, SmaxO={SmaxO}) ...', flush=True)
    pr = np.load(os.path.join(D, 'pr.npy'))
    s  = np.load(os.path.join(D, 's.npy'))
    ev = np.load(os.path.join(D, 'ev.npy'))
    od = np.load(os.path.join(D, 'od.npy'))
    memE = bytes(np.load(os.path.join(D, 'memE_bits.npy')))
    memO = bytes(np.load(os.path.join(D, 'memO_bits.npy')))
    assert len(memE) * 8 > SmaxE and len(memO) * 8 > SmaxO

    # Merged inverse table: ev values are odd/0, od values are even -> the two
    # sets are DISJOINT by value-parity, so one int32 table of size Smax+1 holds
    # both position maps. Halves ~8 GB of two separate tables to ~4 GB.
    Smax = meta['Smax']
    inv = np.full(Smax + 1, -1, np.int32)
    inv[ev] = np.arange(ev.size, dtype=np.int32)
    inv[od] = np.arange(od.size, dtype=np.int32)
    # numpy indexing (no pr_l/s_l Python-list copies) saves ~2.4 GB.
    st = dict(SmaxE=SmaxE, SmaxO=SmaxO, Smax=Smax, pr=pr, s=s, Np=pr.size,
              inv=inv, ev=ev, od=od, memE=memE, memO=memO)
    return st

def _isqrt_vec(x):
    k = np.sqrt(x.astype(np.float64)).astype(np.int64)
    k = np.where((k + 1) * (k + 1) <= x, k + 1, k)
    k = np.where(k * k > x, k - 1, k)
    assert np.all(k * k <= x) and np.all((k + 1) * (k + 1) > x)
    return k

def run_range(st, m0, m1, outdir):
    SmaxE = st['SmaxE']; SmaxO = st['SmaxO']
    pr = st['pr']; s = st['s']; Np = st['Np']
    inv = st['inv']; ev = st['ev']; od = st['od']
    memE = st['memE']; memO = st['memO']

    t0 = time.time()
    mm_ = np.arange(m0, m1 + 1, dtype=np.int64)
    bnd = np.empty(mm_.size, dtype=np.int64)
    even = (mm_ & 1) == 0
    bnd[even] = 2 * mm_[even] + _isqrt_vec(121 * mm_[even]) // 5
    bnd[~even] = mm_[~even] + _isqrt_vec(529 * mm_[~even]) // 5
    del mm_

    t_arr = np.zeros(m1 - m0 + 1, np.int32)
    violations = []
    gap_heap = []
    n_tie = 0
    max_top = 0; max_top_m = 0
    rnd = random.Random(42)
    rnd_set = set(rnd.randint(m0, m1) for _ in range(1000))
    force = {m0, m0 + 1, m1, m1 - 1, (m0 + m1) // 2}
    if m0 <= 2382:
        force |= {10, 20, 303, 2382}
    reps = {}
    p0 = 1

    def note(mm, i, j, how, tt, b):
        nonlocal n_tie, max_top, max_top_m
        tt = int(tt)
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
        while pr[p0] < mm: p0 += 1
        b = int(bnd[mm - m0])
        if mm & 1:
            i = p0; hit = None
            while i < Np and pr[i] <= b:
                si = s[i]
                if i & 1:
                    sj = mm - si
                    if 1 <= sj <= SmaxE and (memE[sj >> 3] >> (sj & 7)) & 1:
                        j = int(inv[sj]) * 2
                        if j < i: hit = (i, j, 'sum'); break
                else:
                    if si == mm: hit = (i, 0, 'k1'); break
                    sj = mm - si
                    if 2 <= sj <= SmaxO and (memO[sj >> 3] >> (sj & 7)) & 1:
                        j = int(inv[sj]) * 2 + 1
                        if j < i: hit = (i, j, 'sum'); break
                i += 1
            if hit is not None:
                i, j, how = hit
                note(mm, i, j, how, pr[i], b)
        else:
            bt = None
            if mm <= SmaxO and (memO[mm >> 3] >> (mm & 7)) & 1:
                i = int(inv[mm]) * 2 + 1
                bt = (pr[i], i, 0, 'k1')
            c = 1
            while True:
                v = mm + ev[c]
                if v > SmaxE: break
                if (memE[v >> 3] >> (v & 7)) & 1:
                    a = int(inv[v]); i = a * 2; tt = pr[i]
                    if bt is None or tt < bt[0]: bt = (tt, i, 2 * c, 'diff')
                    break
                if bt is not None and v >= bt[0]: break
                c += 1
            d = 0
            while True:
                v = mm + od[d]
                if v > SmaxO: break
                if (memO[v >> 3] >> (v & 7)) & 1:
                    bpos = int(inv[v]); i = bpos * 2 + 1; tt = pr[i]
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
    zero_idx = np.nonzero(t_arr == 0)[0]
    for z in zero_idx:
        violations.append(int(z) + m0)
    gap_heap.sort()
    out = {
        'm0': m0, 'm1': m1, 'n': m1 - m0 + 1,
        'time_s': round(el, 1), 'load_s': round(time.time() - t0, 1),
        'code': 'verify_fb5.py',
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
    os.makedirs(outdir, exist_ok=True)
    np.save(os.path.join(outdir, out['t_file']), t_arr)
    fn = f'verify_all3_{m0}_{m1}.json'
    json.dump(out, open(os.path.join(outdir, fn), 'w'), indent=1)
    print(f'done: {out["n"]} m in {el:.0f}s ({out["n"]/el/1000:.1f}k/s); '
          f'violations={out["n_viol"]}; ties={n_tie}; wrote {fn}', flush=True)
    return out

def main():
    D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data5')
    m0 = int(sys.argv[1]); m1 = int(sys.argv[2])
    assert 1 <= m0 <= m1
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(os.path.abspath(__file__))
    st = build(D)
    run_range(st, m0, m1, outdir)

if __name__ == '__main__':
    main()
