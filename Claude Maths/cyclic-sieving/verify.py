#!/usr/bin/env python3
"""
verify.py -- the headline verification, self-contained entry point.

Checks the Criterion for the multiplier cyclic sieving phenomenon on binary necklaces:

    For n >= 3 and a unit v mod n of multiplicative order m, put
        f_e = gcd(v^e - 1, n)      (# of points of Z/n held fixed by v^e)
    Then, for EVERY k coprime to n,
        ( [n choose k]_q / [n]_q )  evaluated at  q = zeta_m
            ==  #{ k-subset necklaces of Z/n fixed by the multiplier x -> v x }
    holds  <=>  for every proper divisor e of m:
        (i)  phi(f_e) <= 2      i.e.  f_e in {1,2,3,4,6}   [crystallographic restriction]
        (ii) rad(f_e) | m

All arithmetic is exact (integers + cyclotomic polynomial arithmetic). No floats.

Usage:
    python3 verify.py            # n = 4..40, ~1 minute
    python3 verify.py 4 96       # the full published range (slow; run in background)
"""
import sys
from math import gcd
from lib2 import order, fixY_alpha, phi, rad
from crit import Yat_fast


def criterion(n, v, m):
    """The predicted verdict: does sieving hold at the multiplier v?"""
    for e in range(1, m):
        if m % e:
            continue
        fe = gcd(pow(v, e, n) - 1, n)
        if phi(fe) > 2 or m % rad(fe) != 0:
            return False
    return True


def truth(n, v, m):
    """The computed verdict: test every admissible k."""
    for k in range(1, n):
        if gcd(n, k) != 1:
            continue
        if Yat_fast(n, k, m) != fixY_alpha(n, [k, n - k], v, m):
            return False
    return True


def main(lo, hi):
    tested = held = mismatches = 0
    for n in range(lo, hi + 1):
        for v in range(2, n):
            if gcd(v, n) != 1:
                continue
            m = order(v, n)
            t, p = truth(n, v, m), criterion(n, v, m)
            tested += 1
            held += t
            if t != p:
                mismatches += 1
                print(f"  COUNTEREXAMPLE  n={n} v={v} m={m} "
                      f"f={gcd(v-1,n)} computed={t} predicted={p}", flush=True)
        if n % 10 == 0:
            print(f"  n<={n}: multipliers={tested} sieving holds={held} "
                  f"mismatches={mismatches}", flush=True)
    print(f"\nRANGE n={lo}..{hi}")
    print(f"  multipliers tested : {tested}")
    print(f"  sieving holds      : {held}")
    print(f"  COUNTEREXAMPLES    : {mismatches}")
    return mismatches


if __name__ == '__main__':
    lo = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    sys.exit(1 if main(lo, hi) else 0)
