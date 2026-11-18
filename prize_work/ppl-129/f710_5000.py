"""
Push f(n) to n=5000 to see if the ratio is converging to sqrt(2)/e or leveling off.
"""
import math
import time
from collections import deque

def hopcroft_karp(n, f):
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

def f_exact_hk(n, max_f=10000):
    lo, hi = n + 1, max_f
    if not hopcroft_karp(n, hi):
        while not hopcroft_karp(n, hi):
            hi *= 2
            if hi > 100000:
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
    print(f"Claimed: sqrt(2)/e = {claimed:.6f}")
    print(f"{'n':>6} {'f(n)':>7} {'ratio':>10} {'1/sqrt(lln)':>12} {'time':>6}")
    
    ns = [500, 1000, 2000, 3000, 5000]
    for n in ns:
        t1 = time.time()
        f = f_exact_hk(n)
        dt = time.time() - t1
        if f is None:
            print(f"{n:>6} NOT FOUND ({dt:.1f}s)")
            continue
        ratio = f / (n * math.sqrt(math.log(n)))
        lln = math.log(math.log(n))
        inv_sqrt_lln = 1.0 / math.sqrt(lln)
        print(f"{n:>6} {f:>7} {ratio:>10.6f} {inv_sqrt_lln:>12.6f} {dt:>6.2f}s")
    
    print(f"\nTotal: {time.time()-t0:.1f}s")
    
    # Fit: ratio ≈ C + A/sqrt(loglog n) to estimate C
    data = [(n, f_exact_hk(n)) for n in [500, 1000, 2000, 3000, 5000]]
    data = [(n, f) for n, f in data if f is not None]
    if len(data) >= 2:
        # Simple: check if ratio - claimed is proportional to 1/sqrt(loglog n)
        print(f"\n{'n':>6} {'ratio':>10} {'ratio-claimed':>14} {'1/sqrt(lln)':>12} {'(ratio-c)/1/sqrt(lln)':>20}")
        for n, f in data:
            ratio = f / (n * math.sqrt(math.log(n)))
            lln = math.log(math.log(n))
            inv = 1.0 / math.sqrt(lln)
            diff = ratio - claimed
            print(f"{n:>6} {ratio:>10.6f} {diff:>14.6f} {inv:>12.6f} {diff/inv:>20.6f}")

if __name__ == "__main__":
    main()
