"""Per-class check: multiset of Zagier cycle lengths  ==  multiset of #R over necklaces.

Zagier-reduced form of discriminant D: (A,B,C) with A,C >= 1 and B > A+C.
They fall into cycles; consecutive forms satisfy A' = C and B' = -B (mod 2C).
Zagier's theory: the cycles <-> the classes, and a cycle's length is the period of
the MINUS (Hirzebruch-Jung) continued fraction.  The plus<->minus conversion
   [a1,a2,a3,a4,...]_+ = [a1+1, 2^(a2-1), a3+2, 2^(a4-1), ...]_-
makes that period length equal to  a2+a4+... = #R.
"""
import sys
from math import isqrt
from collections import Counter
from necktrace import necklaces_by_trace
from zagier import divisors, prim_root

def zagier_forms(D):
    """all (A,B,C), A,C>=1, B>A+C, B^2-4AC=D.
    Bound: S=A+C satisfies B>=S+1 and S^2 > B^2-D, so D > 2S+1, i.e. S < (D-1)/2,
    hence B^2 = D+4AC <= D+S^2 < D + D^2/4."""
    out = []
    Bmax = isqrt(D + D * D // 4) + 2
    s = isqrt(D)
    for B in range(s + 1, Bmax + 1):
        m = B * B - D
        if m % 4:
            continue
        AC = m // 4
        if AC < 1:
            continue
        for A in divisors(AC):
            C = AC // A
            if A + C < B:
                out.append((A, B, C))
    return out

def zagier_cycles(D):
    forms = zagier_forms(D)
    fs = set(forms)
    index = {}
    for f in forms:
        index.setdefault(f[0], []).append(f)
    cycles = []
    unseen = set(forms)
    while unseen:
        f = min(unseen)
        cyc = []
        g = f
        while g in unseen:
            unseen.discard(g)
            cyc.append(g)
            A, B, C = g
            nxt = [h for h in index.get(C, []) if (h[1] + B) % (2 * C) == 0]
            if len(nxt) != 1:
                return None          # successor not unique -> our rule is wrong
            g = nxt[0]
        cycles.append(len(cyc))
    return sorted(cycles)

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    nt = necklaces_by_trace(T)
    ok_all = True
    for t in range(3, T + 1):
        nk = nt.get(t, [])
        if not nk:
            continue
        rs = sorted(prim_root(r).count('R') for r in nk)
        zc = zagier_cycles(t * t - 4)
        ok = (zc == rs)
        ok_all &= ok
        flag = "" if ok else "   <<< MISMATCH"
        print(f"t={t:>3}  zagier cycle lengths {str(zc):<40} #R per class {str(rs):<40}{flag}")
    print("\nPER-CLASS MULTISETS MATCH EVERYWHERE" if ok_all else "\n*** MISMATCH ***")
