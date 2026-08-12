"""The cut-sum invariant  Delta(w)  of a binary necklace, and what it knows.

S(w)   = sum over the n cuts of the loop of the cut word's matrix   (well defined)
Delta  = (tr S)^2 - 4 det S                                          (an integer)

Checks:
  (1)  det S = n + sum over cut-pairs of tr[u,v]      (Fricke / commutator identity)
  (2)  Delta = -4 det(2S - nt I)/4 ... i.e. Delta = q(2S-ntI), q = -det on traceless
  (3)  does Delta separate the classes of a fixed trace t (= a fixed discriminant)?
"""
import sys
from collections import defaultdict
from necktrace import necklaces_by_trace

def mul(A, B):
    return (A[0]*B[0] + A[1]*B[2], A[0]*B[1] + A[1]*B[3],
            A[2]*B[0] + A[3]*B[2], A[2]*B[1] + A[3]*B[3])

L = (1, 1, 0, 1)
R = (1, 0, 1, 1)

def mat(w):
    M = (1, 0, 0, 1)
    for ch in w:
        M = mul(M, L if ch == 'L' else R)
    return M

def runs_to_word(runs):
    return "".join(("L" if i % 2 == 0 else "R") * a for i, a in enumerate(runs))

def cut_data(w):
    n = len(w)
    Ms = [mat(w[i:] + w[:i]) for i in range(n)]
    S = [0, 0, 0, 0]
    for M in Ms:
        for k in range(4):
            S[k] += M[k]
    detS = S[0]*S[3] - S[1]*S[2]
    trS = S[0] + S[3]
    return n, Ms, tuple(S), trS, detS, trS*trS - 4*detS

def check_identities(w):
    n, Ms, S, trS, detS, Delta = cut_data(w)
    t = Ms[0][0] + Ms[0][3]
    # (1) commutator identity:  det S = n + sum_{i<j} tr[u,v],  u,v the two arcs
    tot = 0
    for i in range(n):
        for j in range(i + 1, n):
            u = mat(w[i:j])                      # arc from cut i to cut j
            v = mat(w[j:] + w[:i])               # complementary arc
            x, y, z = u[0]+u[3], v[0]+v[3], t
            tot += x*x + y*y + z*z - x*y*z - 2   # Fricke: tr of the commutator
    ok1 = (detS == n + tot)
    # (2) Delta = q(2S - nt I) with q([[a,b],[c,-a]]) = a^2 + bc
    a = S[0] - S[3]
    ok2 = (Delta == a*a + 4*S[1]*S[2])
    return ok1, ok2, t, n, Delta

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "check"
    if mode == "check":
        from cutsum import necklaces
        bad = 0
        for n in range(2, 10):
            for w in necklaces(n):
                if mat(w)[0] + mat(w)[3] < 3:
                    continue
                o1, o2, t, nn, D = check_identities(w)
                if not (o1 and o2):
                    bad += 1
                    print("FAIL", w, o1, o2)
        print("identity check done, failures:", bad)
    else:
        T = int(sys.argv[2]) if len(sys.argv) > 2 else 60
        nt = necklaces_by_trace(T)
        print(f"{'t':>4} {'h':>4} {'#distinct Delta':>16} {'separates?':>11}   Deltas")
        bad = []
        for t in range(3, T + 1):
            ws = [runs_to_word(r) for r in nt.get(t, [])]
            if not ws:
                continue
            ds = [cut_data(w)[5] for w in ws]
            sep = len(set(ds)) == len(ds)
            if not sep:
                bad.append(t)
            print(f"{t:>4} {len(ws):>4} {len(set(ds)):>16} {str(sep):>11}   {sorted(ds)[:6]}")
        print("\ntraces where Delta FAILS to separate the classes:", bad)
