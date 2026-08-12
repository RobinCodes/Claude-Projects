"""Find the exceptional pairs: distinct classes (mod rot/rev/comp) with the same
length, the same trace, AND the same cut-sum determinant.  What are they?"""
import sys
from collections import defaultdict
from moment import mat, cut_data
from cutsum import necklaces
from config import orbit_rep

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    for n in range(2, N + 1):
        seen = {}
        groups = defaultdict(list)
        for w in necklaces(n):
            M = mat(w)
            t = M[0] + M[3]
            if t < 3:
                continue
            r = orbit_rep(w)
            if r in seen:
                continue
            seen[r] = 1
            _, _, S, trS, detS, Delta = cut_data(w)
            groups[(t, detS)].append((r, S))
        hits = [(k, v) for k, v in groups.items() if len(v) > 1]
        print(f"n={n:>3}  Delta-collisions: {len(hits)}")
        for (t, detS), v in sorted(hits):
            print(f"    t={t} detS={detS}")
            for r, S in v:
                print(f"        {r}   S={S}")
