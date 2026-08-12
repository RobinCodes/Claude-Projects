"""The configuration of moment vectors V(C) attached to the classes of one discriminant.

For a necklace w of length n with trace t:
    V(w) = 2 S(w) - n t I          (integral, traceless)
    Delta(w) = q(V) = a^2 + bc' ... concretely for [[a,b],[c,-a]]: a^2+bc
V is equivariant, so the symmetry group <rot, rev, comp> moves V by conjugation.
Questions: does Delta separate the orbits?  Does sum_C V(C) vanish or have structure?
"""
import sys
from collections import defaultdict
from necktrace import necklaces_by_trace
from moment import mat, cut_data, runs_to_word

def orbit_rep(w):
    n = len(w)
    comp = str.maketrans("LR", "RL")
    outs = set()
    for base in (w, w[::-1], w.translate(comp), w[::-1].translate(comp)):
        for i in range(n):
            outs.add(base[i:] + base[:i])
    return min(outs)

def V(w):
    n, Ms, S, trS, detS, Delta = cut_data(w)
    t = Ms[0][0] + Ms[0][3]
    return (2*S[0] - n*t, 2*S[1], 2*S[2], 2*S[3] - n*t), Delta, t, n

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    nt = necklaces_by_trace(T)
    fails = []
    print(f"{'t':>4} {'h':>4} {'orbits':>7} {'#Delta':>7} {'sep':>5}   sum_C V(C)")
    for t in range(3, T + 1):
        ws = [runs_to_word(r) for r in nt.get(t, [])]
        if not ws:
            continue
        orb = {}
        tot = [0, 0, 0, 0]
        for w in ws:
            v, D, _, _ = V(w)
            orb.setdefault(orbit_rep(w), D)
            for k in range(4):
                tot[k] += v[k]
        ds = list(orb.values())
        sep = len(set(ds)) == len(ds)
        if not sep:
            fails.append(t)
        print(f"{t:>4} {len(ws):>4} {len(orb):>7} {len(set(ds)):>7} {str(sep):>5}   {tuple(tot)}")
    print("\ntraces where Delta fails to separate ORBITS:", fails)
