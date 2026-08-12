"""Binary necklaces over {L,R} -> trace of the SL_2(Z) matrix product.

The trace is a cyclic invariant, so it is a function on necklaces.
Question: how many necklaces have a given trace t, and what is that number?
"""
import sys
from math import isqrt, gcd
from collections import defaultdict

# ---------- enumerate necklaces by trace, via run-length (continued fraction) form ----------
# cyclic word  L^{a1} R^{a2} L^{a3} ... R^{a_{2k}}   <->   matrix prod_{i} [[a_i,1],[1,0]]
# trace = K(a_1..a_2k) + K(a_2..a_{2k-1})   (continuants); K is monotone under appending.

def matmul(A, B):
    return (A[0]*B[0] + A[1]*B[2], A[0]*B[1] + A[1]*B[3],
            A[2]*B[0] + A[3]*B[2], A[2]*B[1] + A[3]*B[3])

def word_matrix(runs):
    M = (1, 0, 0, 1)
    for a in runs:
        M = matmul(M, (a, 1, 1, 0))
    return M

def canon_runs(runs):
    """canonical form of the run sequence under cyclic rotation by 2 (keeps L/R phase)"""
    n = len(runs)
    best = None
    for s in range(0, n, 2):
        cand = tuple(runs[s:] + runs[:s])
        if best is None or cand < best:
            best = cand
    return best

def necklaces_by_trace(T):
    """all primitive-or-not necklaces (both letters present) with trace <= T"""
    out = defaultdict(list)
    seen = set()
    # DFS over run sequences of even length
    stack = [([], (1, 0, 0, 1))]
    while stack:
        runs, M = stack.pop()
        if runs and len(runs) % 2 == 0:
            t = M[0] + M[3]
            if t <= T:
                c = canon_runs(runs)
                if c not in seen:
                    seen.add(c)
                    out[t].append(c)
        # extend
        for a in range(1, T + 1):
            M2 = matmul(M, (a, 1, 1, 0))
            if M2[0] > T:          # continuant already too big
                break
            stack.append((runs + [a], M2))
    return out

# ---------- class numbers: cycles of reduced indefinite binary quadratic forms ----------

def reduced_forms(D):
    """all reduced indefinite forms (a,b,c) of discriminant D>0 nonsquare"""
    s = isqrt(D)
    res = []
    for b in range(1 if D % 2 else 2, s + 1, 2):
        if (D - b * b) % 4:
            continue
        ac = (b * b - D) // 4          # a*c = (b^2-D)/4 < 0
        m = -ac                         # a*(-c) = m > 0
        for a in divisors_signed(m):
            c = ac // a
            # reduced: 0 < b < sqrt(D) and sqrt(D)-b < 2|a| < sqrt(D)+b
            if 0 < b and b * b < D and (s_lt(D, b, 2 * abs(a))) :
                res.append((a, b, c))
    return res

def s_lt(D, b, twoa):
    # sqrt(D) - b < twoa < sqrt(D) + b   with exact integer arithmetic
    # twoa > sqrt(D) - b  <=>  twoa + b > sqrt(D)  <=> (twoa+b)^2 > D
    # twoa < sqrt(D) + b  <=>  twoa - b < sqrt(D)  <=> twoa-b<0 or (twoa-b)^2 < D
    if (twoa + b) ** 2 <= D:
        return False
    d = twoa - b
    if d < 0:
        return True
    return d * d < D

def divisors_signed(m):
    ds = []
    i = 1
    while i * i <= m:
        if m % i == 0:
            ds.append(i); ds.append(-i)
            if i != m // i:
                ds.append(m // i); ds.append(-(m // i))
        i += 1
    return ds

def rho(f, D):
    a, b, c = f
    s = isqrt(D)
    # choose b' = -b mod 2c with  sqrt(D)-2|c| < b' < sqrt(D)
    twoc = 2 * abs(c)
    # b' in (s-2|c|, s]  roughly; solve exactly
    b2 = (-b) % twoc
    # want largest b2 with b2 <= sqrt(D) i.e. b2^2 < D  (D nonsquare)
    while b2 * b2 < D:
        b2 += twoc
    while b2 * b2 > D:
        b2 -= twoc
    # now b2 is the largest value <= sqrt(D) congruent to -b mod 2c
    a2 = c
    c2 = (b2 * b2 - D) // (4 * c)
    return (a2, b2, c2)

def form_cycles(D, primitive_only=True):
    forms = [f for f in reduced_forms(D)]
    if primitive_only:
        forms = [f for f in forms if gcd(gcd(f[0], f[1]), f[2]) == 1]
    fs = set(forms)
    cycles = []
    while fs:
        f = next(iter(fs))
        cyc = []
        g = f
        while g in fs:
            fs.discard(g)
            cyc.append(g)
            g = rho(g, D)
        cycles.append(cyc)
    return cycles

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    nt = necklaces_by_trace(T)
    print(f"{'t':>4} {'#neck':>6} {'D':>7} {'h+(prim)':>9} {'h(all)':>7}   sample")
    for t in range(3, T + 1):
        D = t * t - 4
        nk = nt.get(t, [])
        cp = form_cycles(D, True)
        ca = form_cycles(D, False)
        sample = ""
        if nk:
            sample = " ".join("".join(("L" if i % 2 == 0 else "R") * 1 * a for i, a in enumerate(w))
                              for w in sorted(nk)[:3])
        print(f"{t:>4} {len(nk):>6} {D:>7} {len(cp):>9} {len(ca):>7}   {sample[:60]}")
