"""
Job lab_erdos_20261110_36: fit the 1/sqrt(loglog n) correction to Erdos #710 / PPL 129.

Model under test (idea card):  ratio(n) = f(n)/(n*sqrt(log n)) = A + B/sqrt(loglog n)
  A compared to Turturean claimed constant sqrt(2)/e = 0.5203.
  B compared to EP lower-bound constant 2/sqrt(e) = 1.2131 (the job-25 dossier wrote
  0.797, an arithmetic error: 2/sqrt(e) = 1.2131 exactly).

Method
  1. Regenerate exact f(n) on a dense grid n=100..5000 with an OPTIMIZED Hopcroft-Karp
     (adjacency built divisor-based in O(E), not O(n*f) per call).
  2. GATE: reproduce job-25 published values AND ErPo80 f(10)=15.
  3. INDEPENDENT VERIFIER: naive O(n*f) HK (job-25 proven code) on a sample + Dinic
     max-flow (different algorithm) on n=500,1000,2000,5000.
  4. Fit via 3 least-squares routes (numpy.lstsq, hand-rolled normal eqs, two-anchored)
     + model selection (Turturean ratio=C, EP-lb ratio=D/sqrtloglog, combined).
  5. Kill test: ratio monotone in n? A in [0.50, 0.75]?
"""
import math, json, sys, time
from collections import deque

sys.setrecursionlimit(1000000)
NMAX = 5000
HI_FACTOR = 3


def build_spf(limit):
    spf = list(range(limit + 1))
    if limit >= 0: spf[0] = 0
    if limit >= 1: spf[1] = 1
    for i in range(2, int(limit**0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


SPF = build_spf(4 * NMAX)


def divisors_le(m, n):
    """All divisors d of m with 1 <= d <= n, ascending."""
    if m <= 0:
        return []
    divs = [1]
    tmp = m
    while tmp > 1:
        p = SPF[tmp]
        base = list(divs)          # divisors from previous primes
        pk = 1
        while tmp % p == 0:
            tmp //= p
            pk *= p
            divs += [d * pk for d in base]
    return [d for d in sorted(set(divs)) if d <= n]


def build_adj_full(n, hi):
    """adj_full[k] = sorted list of right indices j (right value = n+1+j, j<hi-1) with k | (n+1+j)."""
    adj = [[] for _ in range(n + 1)]
    for j in range(hi - 1):
        m = n + 1 + j
        for d in divisors_le(m, n):
            adj[d].append(j)
    return adj


def hk_match_exists(n, adj_full, f):
    """True iff there is a matching covering all left [1..n] inside (n, n+f)
    = right values {n+1, ..., n+f-1} = indices [0..f-2]."""
    cap = f - 1
    if cap < n:
        return False
    matchL = [-1] * (n + 1)
    matchR = [-1] * cap
    INF = float('inf')
    dist = [0] * (n + 1)

    def bfs():
        q = deque()
        for k in range(1, n + 1):
            if matchL[k] == -1:
                dist[k] = 0
                q.append(k)
            else:
                dist[k] = INF
        found = False
        while q:
            k = q.popleft()
            for j in adj_full[k]:
                if j >= cap:
                    break
                ml = matchR[j]
                if ml == -1:
                    found = True
                elif dist[ml] > dist[k] + 1:
                    dist[ml] = dist[k] + 1
                    q.append(ml)
        return found

    def dfs(k):
        for j in adj_full[k]:
            if j >= cap:
                break
            ml = matchR[j]
            if ml == -1 or (dist[ml] == dist[k] + 1 and dfs(ml)):
                matchL[k] = j
                matchR[j] = k
                return True
        dist[k] = INF
        return False

    size = 0
    while bfs():
        for k in range(1, n + 1):
            if matchL[k] == -1 and dfs(k):
                size += 1
    return size == n


def f_exact(n):
    hi = HI_FACTOR * n
    adj_full = build_adj_full(n, hi)
    while not hk_match_exists(n, adj_full, hi):
        hi *= 2
        adj_full = build_adj_full(n, hi)
        if hi > 200000:
            return None
    lo = n + 1
    while lo < hi:
        mid = (lo + hi) // 2
        if hk_match_exists(n, adj_full, mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def naive_hk(n, f):
    """Job-25 proven code: rebuild adjacency O(n*f) per call."""
    right_list = list(range(n + 1, n + f))
    if len(right_list) < n:
        return False
    m_to_idx = {m: i for i, m in enumerate(right_list)}
    adj = {}
    for k in range(1, n + 1):
        adj[k] = [m_to_idx[m] for m in right_list if m % k == 0]
        if not adj[k]:
            return False
    matchL, matchR = {}, {}
    INF = float('inf')
    dist = {}

    def bfs():
        queue = deque()
        for k in range(1, n + 1):
            dist[k] = 0 if k not in matchL else INF
            if k not in matchL:
                queue.append(k)
        found = False
        while queue:
            k = queue.popleft()
            for j in adj[k]:
                ml = matchR.get(j, -1)
                if ml == -1:
                    found = True
                elif dist.get(ml, INF) > dist[k] + 1:
                    dist[ml] = dist[k] + 1
                    queue.append(ml)
        return found

    def dfs(k):
        for j in adj[k]:
            ml = matchR.get(j, -1)
            if ml == -1 or (dist.get(ml, INF) == dist[k] + 1 and dfs(ml)):
                matchL[k] = j
                matchR[j] = k
                return True
        dist[k] = INF
        return False

    matching_size = 0
    while bfs():
        for k in range(1, n + 1):
            if k not in matchL and dfs(k):
                matching_size += 1
    return matching_size == n


def naive_f_exact(n, max_f=20000):
    lo, hi = n + 1, max_f
    if not naive_hk(n, hi):
        hi = max_f * 2
        while not naive_hk(n, hi):
            hi *= 2
            if hi > 100000:
                return None
    while lo < hi:
        mid = (lo + hi) // 2
        if naive_hk(n, mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def dinic_maxflow(n, cap):
    """Independent verifier. Source S, sink T. Left k->cap edge to right j iff k|(n+1+j), j<cap.
       Returns True iff max flow == n (perfect matching)."""
    if cap < n:
        return False
    S = 0
    L0 = 1            # left k occupies node 1 + (k-1)
    R0 = 1 + n        # right j occupies node R0 + j
    T = 1 + n + cap
    N = T + 1
    to, nxt, capf, head = [], [], [], [-1] * N

    def add(u, v, c):
        to.append(v); capf.append(c); nxt.append(head[u]); head[u] = len(to) - 1
        to.append(u); capf.append(0); nxt.append(head[v]); head[v] = len(to) - 1

    for k in range(1, n + 1):
        add(S, L0 + k - 1, 1)
    for j in range(cap):
        add(R0 + j, T, 1)
    for k in range(1, n + 1):
        # right values n+1 .. n+cap-1 divisible by k
        # first multiple >= n+1
        start = ((n + 1 + k - 1) // k) * k
        m = start
        while m < n + cap:
            j = m - (n + 1)
            if 0 <= j < cap:
                add(L0 + k - 1, R0 + j, 1)
            m += k
    from collections import deque as dq
    level = [0] * N
    it = [0] * N
    flow = 0

    def bfs2():
        for i in range(N):
            level[i] = -1
        q = dq([S])
        level[S] = 0
        while q:
            u = q.popleft()
            e = head[u]
            while e != -1:
                if capf[e] > 0 and level[to[e]] < 0:
                    level[to[e]] = level[u] + 1
                    q.append(to[e])
                e = nxt[e]
        return level[T] >= 0

    def dfs2(u, f):
        if u == T:
            return f
        e = it[u]
        while e != -1:
            v = to[e]
            if capf[e] > 0 and level[v] == level[u] + 1:
                ret = dfs2(v, f if f < capf[e] else capf[e])
                if ret > 0:
                    capf[e] -= ret
                    capf[e ^ 1] += ret
                    return ret
            e = nxt[e]
            it[u] = e
        return 0

    while bfs2():
        for i in range(N):
            it[i] = head[i]
        while True:
            pushed = dfs2(S, 1 << 30)
            if pushed == 0:
                break
            flow += pushed
            if flow >= n:
                return True
    return flow >= n


def ratio(n, f):
    return f / (n * math.sqrt(math.log(n)))


def main():
    t0 = time.time()
    grid = [100, 125, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 1250, 1500, 1750,
            2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000]
    published = {100: 161, 500: 878, 1000: 1817, 2000: 3815, 3000: 5841, 5000: 10081,
                 10: 15, 1: 2, 2: 3, 3: 4, 4: 6, 20: 29, 30: 46}

    data = {}
    print(f"{'n':>6} {'f(n)':>7} {'ratio':>9} {'time':>7}")
    for n in grid:
        t1 = time.time()
        f = f_exact(n)
        dt = time.time() - t1
        data[n] = f
        if f is None:
            print(f"{n:>6} NOT FOUND ({dt:.1f}s)")
        else:
            print(f"{n:>6} {f:>7} {ratio(n, f):>9.6f} {dt:>6.2f}s")
    print(f"Data gen: {time.time()-t0:.1f}s")

    # Gate
    gate_ok = True
    for n, v in published.items():
        if n in data:
            got = data[n]
            if got != v:
                print(f"GATE FAIL n={n}: got {got}, expected {v}")
                gate_ok = False
    # small-n gate via naive (f(1..4), f(10), f(20), f(30))
    for n in [1, 2, 3, 4, 10, 20, 30]:
        fn = naive_f_exact(n) if n <= 100 else None
        if fn is not None and fn != published[n]:
            print(f"SMALL-GATE FAIL n={n}: naive {fn}, expected {published[n]}")
            gate_ok = False
    print("GATE:", "PASS" if gate_ok else "FAIL")

    # Independent verifier: naive HK vs optimized HK on a sample
    print("\n--- cross-check optimized vs naive HK ---")
    xcheck_ok = True
    for n in [100, 200, 500, 1000]:
        a, b = data.get(n), naive_f_exact(n)
        ok = a == b
        xcheck_ok = xcheck_ok and ok
        print(f"n={n:>5}: optimized={a} naive={b} {'OK' if ok else 'MISMATCH'}")

    # Independent verifier: Dinic vs HK on key n
    print("\n--- cross-check HK vs Dinic max-flow ---")
    dinic_ok = True
    for n in [500, 1000, 2000, 5000]:
        f = data.get(n)
        if f is None:
            continue
        # find minimal f with Dinic by binary search (use f as upper, f-1 as false)
        lo, hi = n + 1, f
        while not dinic_maxflow(n, hi):
            hi += 1
        while lo < hi:
            mid = (lo + hi) // 2
            if dinic_maxflow(n, mid):
                hi = mid
            else:
                lo = mid + 1
        d = lo
        ok = d == f
        dinic_ok = dinic_ok and ok
        print(f"n={n:>5}: HK={f} Dinic={d} {'OK' if ok else 'MISMATCH'}")

    # Save dataset
    out = {n: data[n] for n in sorted(data) if data[n] is not None}
    with open("fit_data.json", "w") as fh:
        json.dump(out, fh)

    # ---------- FIT ----------
    pts = [(n, data[n]) for n in sorted(data) if data[n] is not None]
    xs = [1.0 / math.sqrt(math.log(math.log(n))) for n, _ in pts]
    ys = [ratio(n, f) for n, f in pts]
    print(f"\nFit on {len(pts)} points. x=1/sqrt(loglog) in [{min(xs):.4f},{max(xs):.4f}]")

    # Route 1: numpy lstsq
    import numpy as np
    X = np.column_stack([np.ones(len(xs)), np.array(xs)])
    (A_np, B_np), res, rank, sv = np.linalg.lstsq(X, np.array(ys), rcond=None)
    pred = A_np + B_np * np.array(xs)
    ss_res = float(np.sum((np.array(ys) - pred) ** 2))
    ss_tot = float(np.sum((np.array(ys) - np.mean(ys)) ** 2))
    r2_np = 1 - ss_res / ss_tot

    # Route 2: hand-rolled normal equations
    Sx = sum(xs); Sy = sum(ys); Sxx = sum(x * x for x in xs); Sxy = sum(x * y for x, y in zip(xs, ys))
    N = len(xs)
    det = N * Sxx - Sx * Sx
    B_h = (N * Sxy - Sx * Sy) / det
    A_h = (Sy - B_h * Sx) / N

    # Route 3: two-anchored (first and last point)
    B_t = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
    A_t = ys[0] - B_t * xs[0]

    # Model selection
    # (a) Turturean: ratio = C
    C = sum(ys) / N
    ss_a = sum((y - C) ** 2 for y in ys)
    # (b) EP-lb: ratio = D/sqrtloglog (no intercept, through origin in x)
    D = sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)
    ss_b = sum((y - D * x) ** 2 for x, y in zip(xs, ys))
    # (c) combined: A + B x  -> ss_res
    k_a, k_b, k_c = 1, 1, 2
    aic_a = N * math.log(ss_a / N) + 2 * k_a
    aic_b = N * math.log(ss_b / N) + 2 * k_b
    aic_c = N * math.log(ss_res / N) + 2 * k_c

    # monotonicity of ratio
    mono_dec = all(ys[i + 1] <= ys[i] + 1e-9 for i in range(len(ys) - 1))
    nonmono = sum(1 for i in range(len(ys) - 1) if ys[i + 1] > ys[i] + 1e-9)

    SQRT2_OVER_E = math.sqrt(2) / math.e
    TWO_OVER_SQRTE = 2 / math.sqrt(math.e)

    result = {
        "gate": "PASS" if gate_ok else "FAIL",
        "xcheck_opt_vs_naive": "PASS" if xcheck_ok else "FAIL",
        "xcheck_hk_vs_dinic": "PASS" if dinic_ok else "FAIL",
        "n_points": N,
        "x_range": [round(min(xs), 5), round(max(xs), 5)],
        "A_numpy": round(float(A_np), 5), "B_numpy": round(float(B_np), 5), "r2": round(r2_np, 5),
        "A_hand": round(float(A_h), 5), "B_hand": round(float(B_h), 5),
        "A_twoanchor": round(float(A_t), 5), "B_twoanchor": round(float(B_t), 5),
        "turturean_C": round(float(C), 5),
        "eplb_D": round(float(D), 5),
        "aic": {"turturean_const": round(aic_a, 2), "eplb_origin": round(aic_b, 2), "combined": round(aic_c, 2)},
        "r2_const_model": round(1 - ss_a / ss_tot, 5),
        "r2_origin_model": round(1 - ss_b / ss_tot, 5),
        "mono_decreasing": bool(mono_dec),
        "num_nonmonotone_increases": int(nonmono),
        "A_in_[0.50,0.75]": bool(0.50 <= float(A_np) <= 0.75),
        "sqrt2_over_e": round(SQRT2_OVER_E, 5),
        "two_over_sqrt_e": round(TWO_OVER_SQRTE, 5),
        "data": {str(k): v for k, v in out.items()},
        "ratios": {str(n): round(ratio(n, data[n]), 5) for n in sorted(data) if data[n] is not None},
    }
    with open("fit_results.json", "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"\n=== FIT RESULTS ===")
    print(f"numpy  A={A_np:.5f}  B={B_np:.5f}  R2={r2_np:.5f}")
    print(f"hand   A={A_h:.5f}  B={B_h:.5f}")
    print(f"2anchor A={A_t:.5f}  B={B_t:.5f}")
    print(f"Turturean const model: C={C:.5f}  R2={1-ss_a/ss_tot:.5f}")
    print(f"EP-lb origin model:    D={D:.5f}  R2={1-ss_b/ss_tot:.5f}")
    print(f"AIC: const={aic_a:.2f}  origin={aic_b:.2f}  combined={aic_c:.2f}")
    print(f"mono_decreasing={mono_dec}  nonmono_increases={nonmono}")
    print(f"A in [0.50,0.75] = {0.50 <= float(A_np) <= 0.75}")
    print(f"Target A: sqrt(2)/e={SQRT2_OVER_E:.5f};  Target B: 2/sqrt(e)={TWO_OVER_SQRTE:.5f}")
    print(f"Total: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
