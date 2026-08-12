"""Word statistics on the necklaces of a given trace.

For a cyclic word w over {L,R} with matrix M(w) in SL_2(Z), trace t:
  ell(w) = |w|,  psi(w) = #L - #R   (the Rademacher symbol of the conjugacy class)
Tabulate the multisets over all necklaces of trace t and look for exact laws.
"""
import sys
from collections import defaultdict
from necktrace import necklaces_by_trace

def numdiv(n):
    c = 0
    i = 1
    while i * i <= n:
        if n % i == 0:
            c += 2 if i != n // i else 1
        i += 1
    return c

def stats(T):
    nt = necklaces_by_trace(T)
    out = {}
    for t in range(3, T + 1):
        data = []
        for runs in nt.get(t, []):
            nL = sum(runs[0::2])
            nR = sum(runs[1::2])
            data.append((nL + nR, nL - nR))
        out[t] = data
    return out

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    S = stats(T)
    print(f"{'t':>4} {'h':>4} {'sum_ell':>8} {'divsum':>8} {'sum_psi2':>9} {'sum|psi|':>9}  psi multiset")
    for t in range(3, T + 1):
        d = S[t]
        h = len(d)
        se = sum(x[0] for x in d)
        sp2 = sum(x[1] ** 2 for x in d)
        spa = sum(abs(x[1]) for x in d)
        ds = sum(numdiv(a * (t - a) - 1) for a in range(1, t))
        psis = sorted(x[1] for x in d)
        print(f"{t:>4} {h:>4} {se:>8} {ds:>8} {sp2:>9} {spa:>9}  {psis}")
