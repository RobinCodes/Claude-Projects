"""Trace spectrum of cyclic binary words with p L's and q R's.

min trace is p*q+2 (the clumped word L^p R^q); max is the balanced/Christoffel one.
Question: which integers in between actually occur?
"""
import sys
from collections import defaultdict

def spectra(N):
    """spec[(p,q)] = set of traces of words with p L's, q R's,  p+q = n <= N"""
    spec = defaultdict(set)
    # DFS over all words, carrying the 2x2 product
    def rec(n, p, a, b, c, d):
        if n:
            spec[(p, n - p)].add(a + d)
        if n == N:
            return
        # append L: M * [[1,1],[0,1]]
        rec(n + 1, p + 1, a, a + b, c, c + d)
        # append R: M * [[1,0],[1,1]]
        rec(n + 1, p, a + b, b, c + d, d)
    rec(0, 0, 1, 0, 0, 1)
    return spec

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    spec = spectra(N)
    print(f"{'p':>3} {'q':>3} {'min':>6} {'max':>8} {'#':>6} {'firstgap':>9}  missing (first few)")
    for n in range(2, N + 1):
        for p in range(1, n):
            q = n - p
            if p > q:
                continue
            S = spec[(p, q)]
            lo, hi = p * q + 2, max(S)
            miss = [t for t in range(lo, hi + 1) if t not in S]
            fg = miss[0] if miss else None
            print(f"{p:>3} {q:>3} {lo:>6} {hi:>8} {len(S):>6} {str(fg):>9}  {miss[:8]}")
        print()
