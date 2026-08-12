"""The CUT-SUM of a binary necklace.

A necklace of length n can be cut in n places; each cut gives a word, hence a matrix
M_0,...,M_{n-1}, all conjugate (same trace t).  Rotating the necklace only permutes
them, so   S(w) = sum_i M_i   is a well-defined MATRIX invariant of the necklace.

tr S = n*t always.  What is det S?  What is  Delta = (tr S)^2 - 4 det S ?
"""
import sys
from math import gcd
from collections import defaultdict

L = (1, 1, 0, 1)
R = (1, 0, 1, 1)

def mul(A, B):
    return (A[0]*B[0] + A[1]*B[2], A[0]*B[1] + A[1]*B[3],
            A[2]*B[0] + A[3]*B[2], A[2]*B[1] + A[3]*B[3])

def mat(word):
    M = (1, 0, 0, 1)
    for ch in word:
        M = mul(M, L if ch == 'L' else R)
    return M

def necklaces(n):
    seen = set()
    for x in range(1 << n):
        w = "".join('L' if (x >> i) & 1 else 'R' for i in range(n))
        c = min(w[i:] + w[:i] for i in range(n))
        if c not in seen:
            seen.add(c)
            yield c

def cutsum(w):
    n = len(w)
    S = [0, 0, 0, 0]
    for i in range(n):
        M = mat(w[i:] + w[:i])
        for k in range(4):
            S[k] += M[k]
    return tuple(S)

def sqfree_part(m):
    """m = f^2 * core ; return (core, f)"""
    if m == 0:
        return (0, 0)
    s = 1 if m > 0 else -1
    m = abs(m)
    f = 1
    d = 2
    while d * d <= m:
        while m % (d * d) == 0:
            m //= d * d
            f *= d
        d += 1
    return (s * m, f)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    print(f"{'word':<14} {'n':>2} {'t':>5} {'detS':>10} {'Delta':>12} {'Delta/D':>10} {'core':>7}")
    for n in range(1, N + 1):
        for w in necklaces(n):
            M = mat(w)
            t = M[0] + M[3]
            if t < 3:
                continue
            S = cutsum(w)
            trS = S[0] + S[3]
            detS = S[0]*S[3] - S[1]*S[2]
            D = t*t - 4
            Delta = trS*trS - 4*detS
            core, f = sqfree_part(Delta)
            q = Delta / D if D else 0
            print(f"{w:<14} {n:>2} {t:>5} {detS:>10} {Delta:>12} {q:>10.4f} {core:>7}")
        print()
