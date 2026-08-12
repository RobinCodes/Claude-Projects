"""Fast: enumerate every M in SL_2(Z>=0) with trace t, get its L/R word statistics.

M = [[a,b],[c,d]],  ad - bc = 1,  a+d = t,  all entries >= 0.
Then bc = ad-1 =: N, so for each a and each divisor b|N we get exactly one matrix.
The monoid SL_2(Z>=0) is free on L=[[1,1],[0,1]], R=[[1,0],[1,1]], so each such M
is a unique word; the run-length decomposition gives #L, #R in O(log) time.
"""
from math import isqrt

def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i != n // i:
                ds.append(n // i)
        i += 1
    return ds

INF = 1 << 62

def word_counts(a, b, c, d):
    """number of L's and R's in the unique L/R word of [[a,b],[c,d]];
    peel leading letters:  L^{-1}M = [[a-c,b-d],[c,d]],  R^{-1}M = [[a,b],[c-a,d-b]]"""
    nL = nR = 0
    while not (a == 1 and b == 0 and c == 0 and d == 1):
        if a >= c and b >= d:                       # word starts with L
            k = min(a // c if c else INF, b // d if d else INF)
            if k == 0 or k == INF:
                return None
            a -= k * c; b -= k * d; nL += k
        elif c >= a and d >= b:                     # word starts with R
            k = min(c // a if a else INF, d // b if b else INF)
            if k == 0 or k == INF:
                return None
            c -= k * a; d -= k * b; nR += k
        else:
            return None
    return nL, nR

def matrices(t):
    """yield (a,b,c,d) for all M in SL_2(Z>=0), trace t"""
    for a in range(1, t):
        d = t - a
        N = a * d - 1
        if N < 1:
            continue
        for b in divisors(N):
            yield (a, b, N // b, d)

def profile(t):
    n0 = 0
    s_ell = s_psi2 = s_psi4 = s_LR = 0
    for (a, b, c, d) in matrices(t):
        wc = word_counts(a, b, c, d)
        assert wc is not None, (a, b, c, d)
        nL, nR = wc
        n0 += 1
        s_ell += nL + nR
        p = nL - nR
        s_psi2 += p * p
        s_psi4 += p ** 4
        s_LR += nL * nR
    return n0, s_ell, s_psi2, s_psi4, s_LR

if __name__ == "__main__":
    import sys
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    print(f"{'t':>5} {'#M':>7} {'sum_ell':>10} {'sum_psi2':>12} {'sum_LR':>12}")
    for t in range(3, T + 1):
        n0, se, s2, s4, sLR = profile(t)
        print(f"{t:>5} {n0:>7} {se:>10} {s2:>12} {sLR:>12}")
