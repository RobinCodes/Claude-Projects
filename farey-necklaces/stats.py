"""Length statistics of the necklaces with a given SL_2 trace."""
import sys
from math import isqrt, gcd
from collections import defaultdict
from necktrace import necklaces_by_trace, form_cycles

def divisor_sum(n):
    s = 0
    i = 1
    while i * i <= n:
        if n % i == 0:
            s += i
            if i != n // i:
                s += n // i
        i += 1
    return s

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    nt = necklaces_by_trace(T)
    print(f"{'t':>4} {'h':>4} {'sumlen':>8} {'sumruns':>8} {'D=t^2-4':>8} {'sigma(D)':>9} {'lengths'}")
    rows = []
    for t in range(3, T + 1):
        nk = nt.get(t, [])
        lens = sorted(sum(w) for w in nk)
        runs = sum(len(w) for w in nk)
        D = t * t - 4
        rows.append((t, len(nk), sum(lens), runs, D, divisor_sum(D), lens))
        print(f"{t:>4} {len(nk):>4} {sum(lens):>8} {runs:>8} {D:>8} {divisor_sum(D):>9} {lens}")
