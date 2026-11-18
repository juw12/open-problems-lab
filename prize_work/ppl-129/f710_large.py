"""
Extended computation of f(n) for Erdos Problem #710 for n up to 2000.
Uses Hopcroft-Karp with binary search on f.
"""
import math
import time
from collections import deque

def hopcroft_karp(n, f):
    """Returns True iff a perfect matching exists: left={1..n}, right={n+1..n+f-1}, edge k-m iff k|m."""
    right_list = list(range(n + 1, n + f))
    if len(right_list) < n:
        return False

    m_to_idx = {m: i for i, m in enumerate(right_list)}
    adj = {}
    for k in range(1, n + 1):
        adj[k] = [m_to_idx[m] for m in right_list if m % k == 0]
        if not adj[k]:
            return False

    matchL = {}
    matchR = {}
    INF = float('inf')
    dist = {}

    def bfs():
        queue = deque()
        for k in range(1, n + 1):
            if k not in matchL:
                dist[k] = 0
                queue.append(k)
            else:
                dist[k] = INF
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
            if k not in matchL:
                if dfs(k):
                    matching_size += 1
    return matching_size == n


def f_exact_hk(n, max_f=5000):
    lo, hi = n + 1, max_f
    if not hopcroft_karp(n, hi):
        hi = max_f
        while not hopcroft_karp(n, hi):
            hi *= 2
            if hi > 50000:
                return None
    while lo < hi:
        mid = (lo + hi) // 2
        if hopcroft_karp(n, mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def main():
    t0 = time.time()
    claimed = math.sqrt(2.0) / math.e
    print(f"Claimed asymptotic constant: sqrt(2)/e = {claimed:.6f}")
    print(f"{'n':>6} {'f(n)':>6} {'ratio':>10} {'delta':>10}")
    
    ns = [100, 200, 300, 400, 500, 700, 1000, 1500, 2000]
    prev_ratio = None
    results = {}
    for n in ns:
        t1 = time.time()
        f = f_exact_hk(n)
        dt = time.time() - t1
        results[n] = f
        if f is None:
            print(f"{n:>6} NOT FOUND ({dt:.1f}s)")
            continue
        ratio = f / (n * math.sqrt(math.log(n)))
        delta = ratio - prev_ratio if prev_ratio is not None else 0
        prev_ratio = ratio
        print(f"{n:>6} {f:>6} {ratio:>10.6f} {delta:>+10.6f}  ({dt:.2f}s)")

    print(f"\nTotal time: {time.time()-t0:.2f}s")

    # Compare with lower/upper bounds
    print(f"\n{'n':>6} {'f(n)':>6} {'ratio':>10} {'EP_lb':>10} {'EP_ub':>10} {'claimed':>10}")
    for n in ns:
        f = results.get(n)
        if f is None: continue
        ratio = f / (n * math.sqrt(math.log(n)))
        ln = math.log(n)
        lln = math.log(ln)
        ep_lb = (2.0 / math.sqrt(math.e)) / math.sqrt(lln)  # (2/sqrt(e))/sqrt(loglog n)
        ep_ub = 1.7398
        print(f"{n:>6} {f:>6} {ratio:>10.6f} {ep_lb:>10.6f} {ep_ub:>10.6f} {claimed:>10.6f}")


if __name__ == "__main__":
    main()
