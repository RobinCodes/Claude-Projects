#!/usr/bin/env python3
"""
Re-check every claim of PROOF.md.  Exact integer arithmetic throughout.

    python3 verify.py            # default ranges, a couple of minutes
    python3 verify.py --full     # the published ranges, ~1 hour

Non-zero exit status means a claim failed.
"""
import sys, time
from math import gcd, comb
from fractions import Fraction
from core import (pmul, pmod, pdivexact, ptrim, cyclotomic, Y_at_zeta,
                  as_rational_int, orbit_sizes, cycle_sizes_affine,
                  coeff_subset_sum, order_mod, K_set, fixed_count, rad, criterion)

FAIL = 0
def report(name, bad, inst, extra=''):
    global FAIL
    ok = 'OK  ' if bad == 0 else 'FAIL'
    if bad: FAIL += 1
    print(f'  [{ok}] {name:<46} {inst:>9,} checks, {bad} failures {extra}')

# ---------------------------------------------------------------- helpers
def qint(j):  return tuple([1] * j) if j > 0 else ()

def qbinom(a, b):
    if b < 0 or b > a: return ()
    num = (1,)
    for i in range(1, a + 1): num = pmul(num, qint(i))
    den = (1,)
    for i in range(1, b + 1): den = pmul(den, qint(i))
    for i in range(1, a - b + 1): den = pmul(den, qint(i))
    return pdivexact(num, den)

def Qpoly(small):
    acc = (1,)
    for s in small: acc = pmul(acc, tuple([1] + [0] * (s - 1) + [1]))
    return acc

def multipliers(nmax, crit_only=True):
    for n in range(3, nmax + 1):
        for v in range(2, n):
            if gcd(v, n) != 1: continue
            m = order_mod(v, n)
            if m < 2: continue
            if crit_only and not criterion(v, n, m): continue
            yield n, v, m

# ---------------------------------------------------------------- Theorem A
def check_theoremA(nmax, nmax3):
    """#Fix_necklaces(v) = #Fix_words(v)/f, unconditionally, any content."""
    bad = inst = 0
    for n, v, m in multipliers(nmax, crit_only=False):     # criterion NOT assumed
        f = gcd(v - 1, n); sizes = orbit_sizes(v, n)
        for k in range(1, n):
            if gcd(k, n) != 1: continue
            inst += 1
            w = coeff_subset_sum(sizes, k)
            if w % f or w // f != fixed_count(v, n, k): bad += 1
    report('Theorem A, 2 colours (criterion not assumed)', bad, inst, f'n<={nmax}')

    def words_content(sizes, alpha):
        dp = {(0,) * len(alpha): 1}
        for s in sizes:
            nd = {}
            for st, c in dp.items():
                for i in range(len(alpha)):
                    t = list(st); t[i] += s
                    if t[i] <= alpha[i]: nd[tuple(t)] = nd.get(tuple(t), 0) + c
            dp = nd
        return dp.get(tuple(alpha), 0)

    bad = inst = 0
    for n, v, m in multipliers(nmax3, crit_only=False):
        f = gcd(v - 1, n); sizes = orbit_sizes(v, n)
        for a1 in range(1, n - 1):
            for a2 in range(1, n - a1):
                a3 = n - a1 - a2
                if a3 < 1 or gcd(gcd(a1, a2), a3) != 1: continue
                inst += 1
                w = words_content(sizes, (a1, a2, a3))
                tot = sum(words_content(cycle_sizes_affine(v, c, n), (a1, a2, a3))
                          for c in K_set(v, n, m))
                if w % f or w // f != tot // n: bad += 1
    report('Theorem A, 3 colours, all contents', bad, inst, f'n<={nmax3}')

# ---------------------------------------------------------------- Theorem B
def check_theoremB(nmax):
    """Under the criterion P is the subgroup C_L, L in {1,2,3,4,6}, v acts as +-1."""
    bad = inst = 0; shapes = {}
    for n, v, m in multipliers(nmax):
        inst += 1
        P = set()
        for y in range(n):
            z = (v * y) % n; L = 1
            while z != y: z = (v * z) % n; L += 1
            if L < m: P.add(y)
        L = len(P)
        ok = (L in (1, 2, 3, 4, 6) and m % rad(L) == 0
              and P == {(n // L) * i % n for i in range(L)}
              and (L <= 2 or v % L in (1, L - 1)))
        if not ok: bad += 1
        eps = 1 if (L <= 2 or v % L == 1) else -1
        shapes[(L, eps)] = shapes.get((L, eps), 0) + 1
    report('Theorem B, structure of P', bad, inst, f'n<={nmax}')
    print(f'         shapes (L,eps) observed: {sorted(shapes)}')

# ---------------------------------------------------------------- Theorem C
def check_theoremC(nmax_coeff, nmax_e2e):
    """3.3/3.4 coefficient identity, then the end-to-end sieving identity."""
    bad = inst = 0
    for n, v, m in multipliers(nmax_coeff):
        f = gcd(v - 1, n); sizes = orbit_sizes(v, n)
        small = [s for s in sizes if s < m]
        Q = Qpoly(small); r = n % m; u = (sum(small) - r) // m
        Pm = cyclotomic(m)
        svals = sorted({k % m for k in range(1, n) if gcd(k, n) == 1})
        for s in svals:
            for j in range(0, len(Q) // m + 3):
                idx = j * m + s
                q = Q[idx] if idx < len(Q) else 0
                inst += 1
                if r >= 1:
                    lhs = pmod(pmul((q,), qint(r)), Pm)
                    C = comb(u, j) if 0 <= j <= u else 0
                    rhs = pmod(pmul((f * C,), pmod(qbinom(r, s), Pm)), Pm)
                else:
                    if s == 0: continue
                    lhs = pmod(pmul((q,), qint(s)), Pm)
                    C = comb(u - 1, j) if 0 <= j <= u - 1 else 0
                    e = (-(s * (s - 1) // 2)) % m
                    rhs = pmod(pmul((f * C * (-1) ** (s - 1),),
                                    tuple([0] * e + [1])), Pm)
                if ptrim(lhs) != ptrim(rhs): bad += 1
    report('Theorem C, coefficient identity (3.3/3.4)', bad, inst, f'n<={nmax_coeff}')

    sb = nb = ns = nn = 0
    for n in range(3, nmax_e2e + 1):
        for v in range(2, n):
            if gcd(v, n) != 1: continue
            m = order_mod(v, n)
            if m < 2: continue
            crit = criterion(v, n, m); f = gcd(v - 1, n); sizes = orbit_sizes(v, n)
            allok = True
            for k in range(1, n):
                if gcd(k, n) != 1: continue
                Yz = as_rational_int(Y_at_zeta(n, k, m), m)
                if Yz is None or Yz != coeff_subset_sum(sizes, k) // f:
                    allok = False; break
            if crit:
                ns += 1; sb += (not allok)
            else:
                nn += 1; nb += allok
    report('Theorem C, end-to-end sufficiency', sb, ns, f'n<={nmax_e2e}')
    report('necessity (verified, not proved in general)', nb, nn, f'n<={nmax_e2e}')

# ---------------------------------------------------------------- Proposition D
def ratio_if_rational(m, r, s):
    Pm = cyclotomic(m)
    A = list(pmod(qbinom(r, s), Pm)); B = list(pmod(qint(r), Pm))
    A += [0] * (len(Pm) - 1 - len(A)); B += [0] * (len(Pm) - 1 - len(B))
    piv = next((i for i, x in enumerate(B) if x != 0), None)
    if piv is None: return None
    c = Fraction(A[piv], B[piv])
    return c if all(Fraction(A[i]) == c * B[i] for i in range(len(B))) else None

def check_propD(mmax):
    """Hypothesis 3.6 is false; the congruence is necessary; counterexamples have r=m-1."""
    ce = []; viol = 0; inst = 0
    for m in range(2, mmax + 1):
        for r in range(2, m):
            for s in sorted(set([0, r] + list(range(2, r - 1)))):
                if s > r: continue
                inst += 1
                c = ratio_if_rational(m, r, s)
                if c is not None:
                    ce.append((m, r, s, c))
                    if (r - 1 - s * (r - s)) % m: viol += 1     # congruence must hold
    report('Prop D(1): congruence necessary for rationality', viol, inst, f'm<={mmax}')
    nonlast = sum(1 for (m, r, s, c) in ce if r != m - 1)
    notpm1  = sum(1 for (m, r, s, c) in ce if abs(c) != 1)
    report('Prop D(3): every counterexample has r = m-1', nonlast, len(ce), f'm<={mmax}')
    report('Prop D(4): every rational value has modulus 1', notpm1, len(ce), f'm<={mmax}')
    print(f'         Hypothesis 3.6 is FALSE: {len(ce)} counterexamples with m<={mmax};'
          f' smallest {ce[0][:3] if ce else None} with value {ce[0][3] if ce else None}')
    # [r]_zeta irrational for 2<=r<m  (Prop D(2))
    bad = inst = 0
    for m in range(2, mmax + 1):
        for r in range(2, m):
            inst += 1
            A = list(pmod(qint(r), cyclotomic(m)))
            if len(ptrim(A)) <= 1: bad += 1
    report('Prop D(2): [r]_zeta irrational for 2<=r<m', bad, inst, f'm<={mmax}')

# ---------------------------------------------------------------- the n=99 witness
def check_witness():
    n, v, k = 99, 19, 13
    m = order_mod(v, n); f = gcd(v - 1, n)
    Yz = as_rational_int(Y_at_zeta(n, k, m), m)
    true = coeff_subset_sum(orbit_sizes(v, n), k) // f
    ok = (m == 10 and f == 9 and (n % m) == m - 1 and Yz == -9 and true == 84)
    report('witness n=99,v=19,k=13: Y=-9 but count=84', 0 if ok else 1, 1,
           f'Y={Yz}, count={true}')

# ---------------------------------------------------------------- main
if __name__ == '__main__':
    full = '--full' in sys.argv
    t0 = time.time()
    print('Verifying PROOF.md' + (' (full ranges)' if full else ' (default ranges)'))
    check_witness()
    check_theoremA(40 if not full else 60, 18 if not full else 28)
    check_theoremB(120 if not full else 400)
    check_theoremC(60 if not full else 200, 40 if not full else 130)
    check_propD(24 if not full else 60)
    print(f'\n{"ALL CLAIMS VERIFIED" if FAIL == 0 else str(FAIL) + " CLAIM(S) FAILED"}'
          f'   ({time.time() - t0:.0f}s)')
    sys.exit(1 if FAIL else 0)
