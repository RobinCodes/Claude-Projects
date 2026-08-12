"""Does det(cut-sum) separate necklaces that the trace cannot?

Group necklaces of length n by trace. Within a trace class, quotient by the
symmetry group <rotation, reversal, complement> (which det S cannot see anyway).
Then ask whether det S distinguishes the remaining distinct classes.
"""
import sys
from collections import defaultdict
from cutsum import mat, cutsum, necklaces

def orbit(w):
    n = len(w)
    comp = str.maketrans("LR", "RL")
    outs = set()
    for base in (w, w[::-1], w.translate(comp), w[::-1].translate(comp)):
        for i in range(n):
            outs.add(base[i:] + base[:i])
    return min(outs)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    for n in range(2, N + 1):
        by_t = defaultdict(dict)          # t -> {orbit_rep: detS}
        for w in necklaces(n):
            M = mat(w)
            t = M[0] + M[3]
            if t < 3:
                continue
            S = cutsum(w)
            by_t[t][orbit(w)] = S[0] * S[3] - S[1] * S[2]
        nclash = nsep = 0
        detail = []
        for t, d in sorted(by_t.items()):
            if len(d) > 1:                # trace fails to separate these classes
                nclash += len(d) - 1
                vals = set(d.values())
                nsep += len(vals) - 1
                if len(vals) > 1:
                    detail.append((t, sorted(d.items())[:3], sorted(vals)))
        print(f"n={n:>3}  classes with a trace-clash: {nclash:>4}   detS resolves: {nsep:>4}")
        for t, items, vals in detail[:3]:
            print(f"        t={t}: {[k for k, _ in items]} -> detS {vals}")
