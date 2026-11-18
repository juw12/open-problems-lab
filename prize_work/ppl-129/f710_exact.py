"""
Exact computation of f(n) for Erdos Problem #710 (PPL 129).

f(n) = min f such that (n, n+f) contains distinct integers a_1..a_n
with k | a_k for all 1 <= k <= n.

Method: bipartite matching (left = {1..n}, right = {n+1..n+f-1}, edge k-m iff k|m).
f(n) is the smallest f for which a perfect matching exists.

Gates (hand-verified):
  f(1)=2, f(2)=3, f(3)=4, f(4)=6
"""
import math
import sys
import time

def has_perfect_matching(n, f):
    """Check if a perfect matching exists: left={1..n}, right={n+1..n+f-1}, edge k->m iff k|m."""
    right = list(range(n + 1, n + f))
    if len(right) < n:
        return False

    # Build adjacency
    adj = {}
    for k in range(1, n + 1):
        adj[k] = [m for m in right if m % k == 0]
        if len(adj[k]) == 0:
            return False

    # Most-constrained-first ordering
    order = sorted(range(1, n + 1), key=lambda k: len(adj[k]))

    # Backtracking with bitset for speed
    # right vertices indexed 0..(f-2) corresponding to m = n+1+i
    n_right = len(right)
    # map m -> index
    m_to_idx = {m: i for i, m in enumerate(right)}
    adj_idx = {k: [m_to_idx[m] for m in adj[k]] for k in range(1, n + 1)}

    used = [False] * n_right

    def backtrack(idx):
        if idx == n:
            return True
        k = order[idx]
        for i in adj_idx[k]:
            if not used[i]:
                used[i] = True
                if backtrack(idx + 1):
                    return True
                used[i] = False
        return False

    return backtrack(0)


def f_exact(n, max_f=500):
    for f in range(n + 1, max_f + 1):
        if has_perfect_matching(n, f):
            return f
    return None


def main():
    t0 = time.time()
    results = {}
    for n in range(1, 31):
        f = f_exact(n)
        results[n] = f
        if n >= 3 and f is not None:
            log_n = math.log(n)
            ratio = f / (n * math.sqrt(log_n))
            claimed = math.sqrt(2.0) / math.e  # ~0.4839
            print(f"n={n:3d}  f={f:4d}  f/(n*sqrt(ln n))={ratio:.6f}  claimed(sqrt2/e)={claimed:.6f}")
        else:
            print(f"n={n:3d}  f={f}")

    elapsed = time.time() - t0
    print(f"\nTime: {elapsed:.2f}s")
    print("\nFull table:")
    for n in range(1, 31):
        print(f"f({n}) = {results[n]}")

    # Gate check
    gates = {1: 2, 2: 3, 3: 4, 4: 6}
    gate_pass = all(results.get(n) == v for n, v in gates.items())
    print(f"\nGATE: {'PASS' if gate_pass else 'FAIL'} (hand values {gates})")
    if not gate_pass:
        for n, v in gates.items():
            if results.get(n) != v:
                print(f"  MISMATCH n={n}: computed={results.get(n)}, expected={v}")


if __name__ == "__main__":
    main()
