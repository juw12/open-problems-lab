"""
Extended computation of f(n) for Erdos Problem #710 using Hopcroft-Karp matching.
Computes f(n) for n = 1..200 to see the trend of f(n)/(n*sqrt(ln n)).
"""
import math
import time
from collections import deque

def hopcroft_karp(n, f):
    """Check if a perfect matching exists using Hopcroft-Karp.
    Left: {1..n}, Right: {n+1..n+f-1}, edge k-m iff k|m.
    Returns True iff a perfect matching covering all left vertices exists."""
    right_list = list(range(n + 1, n + f))
    if len(right_list) < n:
        return False

    # Build adjacency: left k -> list of right indices (0-based)
    m_to_idx = {m: i for i, m in enumerate(right_list)}
    adj = {}
    for k in range(1, n + 1):
        adj[k] = [m_to_idx[m] for m in right_list if m % k == 0]
        if not adj[k]:
            return False

    # Hopcroft-Karp
    # matchL[k] = right index matched to left k, or -1
    # matchR[j] = left index matched to right j, or -1
    matchL = {}
    matchR = {}

    INF = float('inf')
    dist = {}

    def bfs():
        """Find all shortest augmenting paths. Returns True if one exists."""
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
        """Try to find an augmenting path from left vertex k."""
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


def f_exact_hk(n, max_f=2000):
    """Compute f(n) using Hopcroft-Karp. Binary search on f."""
    lo, hi = n + 1, max_f
    # First check if hi is big enough
    if not hopcroft_karp(n, hi):
        # Extend
        hi = max_f * 2
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

    # Gates from the backtracking run
    gates = {1: 2, 2: 3, 3: 4, 4: 6, 10: 15, 20: 29, 30: 46}

    print("n\tf(n)\tf/(n*sqrt(ln n))\tclaimed(sqrt2/e)")
    results = {}
    for n in list(range(1, 31)) + [40, 50, 60, 80, 100, 120, 150, 200]:
        f = f_exact_hk(n)
        results[n] = f
        if f is None:
            print(f"{n}\tNOT FOUND")
            continue
        claimed = math.sqrt(2.0) / math.e
        if n >= 3:
            ratio = f / (n * math.sqrt(math.log(n)))
            print(f"{n}\t{f}\t{ratio:.6f}\t{claimed:.6f}")
        else:
            print(f"{n}\t{f}\t-\t{claimed:.6f}")

    elapsed = time.time() - t0
    print(f"\nTime: {elapsed:.2f}s")

    # Gate check
    gate_pass = all(results.get(n) == v for n, v in gates.items())
    print(f"GATE: {'PASS' if gate_pass else 'FAIL'}")
    if not gate_pass:
        for n, v in gates.items():
            if results.get(n) != v:
                print(f"  MISMATCH n={n}: HK={results.get(n)}, expected={v}")

    # Trend analysis
    print("\nTrend (f(n)/(n*sqrt(ln n))):")
    for n in [10, 20, 30, 50, 100, 200]:
        if n in results and results[n] is not None:
            ratio = results[n] / (n * math.sqrt(math.log(n)))
            print(f"  n={n:4d}: ratio = {ratio:.6f}")

    claimed = math.sqrt(2.0) / math.e
    print(f"\nClaimed asymptotic constant: sqrt(2)/e = {claimed:.6f}")


if __name__ == "__main__":
    main()
