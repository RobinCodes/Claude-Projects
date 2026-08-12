"""Attack on the OEIS conjecture  A264598(t) = 2 * A257007(t).

A264598(t) = #{ X in the free monoid on L,R : tr X = t }
A257007(t) = #{ Zagier-reduced forms of discriminant t^2-4 }

Claimed proof:
  (i)   A264598(t) = sum over conjugacy classes C of trace t of p(C)  [p = primitive period]
                   = sum_C ( #L(u_C) + #R(u_C) ),  u_C = primitive root word
  (ii)  A257007(t) = sum_C #R(u_C)   [Zagier cycle length = minus-CF period length = #R]
  (iii) L <-> R (conjugation by [[0,1],[1,0]]) permutes the classes of trace t and
        swaps #L with #R, so sum_C #L = sum_C #R.
  => A264598 = 2 * A257007.
"""
import sys
from math import isqrt
from necktrace import necklaces_by_trace
from mat import matrices

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

def zagier_count(D):
    """OEIS A257007 formula: pairs (h,k), k^2<D, k^2=D mod 4, h | (D-k^2)/4, 2h > sqrt(D)-k"""
    s = isqrt(D)
    n = 0
    for k in range(-s, s + 1):
        if k * k >= D or (D - k * k) % 4:
            continue
        m = (D - k * k) // 4
        if m == 0:
            continue
        for h in divisors(m):
            # h > (sqrt(D)-k)/2  <=>  2h + k > sqrt(D)
            if 2 * h + k > 0 and (2 * h + k) ** 2 > D:
                n += 1
    return n

def prim_root(runs):
    w = "".join(("L" if i % 2 == 0 else "R") * a for i, a in enumerate(runs))
    n = len(w)
    for p in range(1, n + 1):
        if n % p == 0 and w == w[:p] * (n // p):
            return w[:p]
    return w

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    nt = necklaces_by_trace(T)
    print(f"{'t':>4} {'#matrices':>10} {'sum p':>7} {'sum#L':>7} {'sum#R':>7} {'Zagier':>7} {'2*Z':>8} {'ok':>5}")
    allok = True
    for t in range(3, T + 1):
        nk = nt.get(t, [])
        if not nk:
            continue
        roots = [prim_root(r) for r in nk]
        sL = sum(u.count('L') for u in roots)
        sR = sum(u.count('R') for u in roots)
        sp = sum(len(u) for u in roots)
        nm = sum(1 for _ in matrices(t))
        z = zagier_count(t * t - 4)
        ok = (nm == sp == sL + sR) and (sL == sR) and (z == sR) and (nm == 2 * z)
        allok &= ok
        print(f"{t:>4} {nm:>10} {sp:>7} {sL:>7} {sR:>7} {z:>7} {2*z:>8} {str(ok):>5}")
    print("\nALL CHECKS PASS" if allok else "\n*** MISMATCH ***")
