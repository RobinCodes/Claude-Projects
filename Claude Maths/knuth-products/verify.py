#!/usr/bin/env python3
"""
verify.py -- exact verification for "Knuth products in linear numeration systems".

All arithmetic is exact (Python ints / Z[x]/(chi)).  No floating point anywhere.

Objects
-------
A *numeration system* is U = (u_0=1 < u_1 < u_2 < ...) of positive integers
obeying  u_k = a_1 u_{k-1} + ... + a_d u_{k-d}  (k >= d),
chi(x) = x^d - a_1 x^{d-1} - ... - a_d  irreducible over Q, beta = dominant root.

Greedy expansion K(n): repeatedly subtract the largest u_k <= n (multiset of indices).
Knuth product:   n o_s m = sum_{k in K(n), l in K(m)} u_{k+l+s}.
   (U = Zeckendorf (1,2,3,5,8,...), s = 2  is Knuth's Fibonacci multiplication.)
Phi(n) = sum_{k in K(n)} beta^k  in Z[beta];   S = Phi(N).
"""
from bisect import bisect_right

# ---------- exact arithmetic in Z[beta] = Z[x]/(chi) ----------
class Zbeta:
    def __init__(self, a):                     # chi = x^d - a1 x^{d-1} - ... - ad
        self.a = list(a); self.d = len(a)
    def red(self, c):
        c = list(c)
        while len(c) > self.d:
            top = c.pop(); k = len(c)
            for i, ai in enumerate(self.a, 1): c[k-i] += ai*top
        while len(c) < self.d: c.append(0)
        return tuple(c)
    def mul(self, x, y):
        c = [0]*(len(x)+len(y)-1)
        for i, xi in enumerate(x):
            if xi:
                for j, yj in enumerate(y): c[i+j] += xi*yj
        return self.red(c)
    def pw(self, k): return self.red([0]*k + [1])
    def add(self, x, y): return tuple(p+q for p, q in zip(x, y))

def is_irreducible(a):
    """rational-root-free + no factorisation into monic integer polys of lower degree
       (brute force over divisors of the constant term for d<=4; enough for our systems)"""
    d = len(a)
    if d == 1: return True                             # linear polys are irreducible
    chi = [1] + [-ai for ai in a]                      # descending coefficients
    def ev(t):
        s = 0
        for c in chi: s = s*t + c
        return s
    const = chi[-1]
    if const == 0: return False
    divs = [t for t in range(-abs(const), abs(const)+1) if t and const % t == 0]
    if any(ev(t) == 0 for t in divs): return False     # kills linear factors
    if d <= 3: return True                             # no linear factor => irreducible
    # d == 4: rule out a product of two monic integer quadratics
    for p1 in range(-20, 21):
        for q1 in range(-20, 21):
            for p2 in range(-20, 21):
                for q2 in range(-20, 21):
                    # (x^2+p1x+q1)(x^2+p2x+q2)
                    if (p1+p2, q1+q2+p1*p2, p1*q2+p2*q1, q1*q2) == tuple(chi[1:]):
                        return False
    return True

class System:
    def __init__(self, a, init, name, terms=260):
        self.a, self.name = list(a), name
        self.R = Zbeta(a)
        U = list(init)
        while len(U) < terms: U.append(sum(c*U[-1-i] for i, c in enumerate(a)))
        assert U[0] == 1 and all(U[i] < U[i+1] for i in range(len(U)-1)), name
        self.U = U
        self.pw = [self.R.pw(k) for k in range(terms)]
    def greedy(self, n):
        out = []
        while n > 0:
            k = bisect_right(self.U, n) - 1
            out.append(k); n -= self.U[k]
        return out
    def Phi(self, n):
        v = (0,)*self.R.d
        for k in self.greedy(n): v = self.R.add(v, self.pw[k])
        return v
    def L(self, z):                                   # the functional of Theorem 1
        return sum(c*self.U[i] for i, c in enumerate(z))
    def circ(self, n, m, s):
        return sum(self.U[i+j+s] for i in self.greedy(n) for j in self.greedy(m))
    def in_S(self, z):                                # membership test for S = Phi(N)
        n = self.L(z)
        return n >= 0 and self.Phi(n) == z

# ---------- the systems ----------
def canon(a):
    d = len(a); U = [1]
    for k in range(1, d): U.append(sum(a[i]*U[k-1-i] for i in range(k)) + 1)
    return U

SYSTEMS = [
    System((2,),      [1],       "base 2                 1,2,4,8,..."),
    System((3,),      [1],       "base 3                 1,3,9,27,..."),
    System((1,1),     [1,2],     "Zeckendorf (Fibonacci) 1,2,3,5,8,..."),
    System((1,1),     [1,3],     "Lucas                  1,3,4,7,11,..."),
    System((1,1),     [1,4],     "Fib-shift              1,4,5,9,14,..."),
    System((2,1),     canon((2,1)), "Pell canonical      1,3,7,17,41,..."),
    System((2,1),     [1,2],     "Pell natural           1,2,5,12,29,..."),
    System((3,1),     canon((3,1)), "x^2=3x+1 canonical  1,4,13,43,..."),
    System((2,2),     canon((2,2)), "x^2=2x+2 canonical  1,3,8,22,..."),
    System((1,1,1),   [1,2,4],   "Tribonacci             1,2,4,7,13,..."),
    System((1,0,1),   canon((1,0,1)), "x^3=x^2+1 canon    1,2,3,4,6,9,..."),
    System((1,3),     canon((1,3)),   "x^2=x+3  (not Pisot) 1,2,5,11,..."),
    System((1,0,0,1), canon((1,0,0,1)), "x^4=x^3+1 canon  1,2,3,4,5,7,..."),
]

def hdr(t): print("\n" + "="*78 + "\n" + t + "\n" + "="*78)

fails = []
def check(cond, msg):
    if not cond: fails.append(msg); print("   *** FAIL:", msg)

# ---------- 0. irreducibility of the characteristic polynomials used ----------
hdr("0.  characteristic polynomials are irreducible (needed for Theorem 1)")
for S in SYSTEMS:
    ok = is_irreducible(S.a)
    print(f"   {S.name:34s} chi irreducible: {ok}")
    check(ok, f"chi reducible for {S.name}")

# ---------- 1. Theorem 1: transport formula ----------
hdr("1.  THEOREM 1   n o_s m = L(beta^s Phi(n) Phi(m))    and   L(Phi(n)) = n")
for S in SYSTEMS:
    ok = True
    for n in range(1, 41):
        if S.L(S.Phi(n)) != n: ok = False
    for s in range(0, 7):
        for n in range(1, 31):
            for m in range(1, 31):
                lhs = S.circ(n, m, s)
                rhs = S.L(S.R.mul(S.pw[s], S.R.mul(S.Phi(n), S.Phi(m))))
                if lhs != rhs: ok = False; break
    print(f"   {S.name:34s} {'OK' if ok else 'FAIL'}")
    check(ok, f"transport formula fails for {S.name}")

# ---------- 2. Theorem 2: criterion  beta^s S S subset S  <=> morphism ----------
hdr("2.  THEOREM 2   [beta^s Phi(n)Phi(m) in S  for all n,m]  <=>  Phi(n o_s m)=beta^s Phi(n)Phi(m)")
print("     and this implies associativity.   N = 24, s = 0..8")
print(f"   {'system':34s} {'closure holds for s in':24s} {'assoc holds for s in'}")
N = 24
for S in SYSTEMS:
    clo, ass = [], []
    for s in range(9):
        c = all(S.in_S(S.R.mul(S.pw[s], S.R.mul(S.Phi(n), S.Phi(m))))
                for n in range(1, N+1) for m in range(1, N+1))
        b = all(S.Phi(S.circ(n, m, s)) == S.R.mul(S.pw[s], S.R.mul(S.Phi(n), S.Phi(m)))
                for n in range(1, N+1) for m in range(1, N+1))
        check(c == b, f"closure<=>morphism mismatch, {S.name}, s={s}")
        a = all(S.circ(S.circ(n, m, s), p, s) == S.circ(n, S.circ(m, p, s), s)
                for n in range(1, N+1) for m in range(1, N+1) for p in range(1, N+1))
        if c: clo.append(s)
        if a: ass.append(s)
        check(not c or a, f"closure but not associative, {S.name}, s={s}")
    print(f"   {S.name:34s} {str(clo):24s} {ass}")

print(fails and f"\n{len(fails)} FAILURES" or "\nno failures so far")

# ---------- 3. the phenomena ----------
hdr("3.  PHENOMENA.  good shifts  G(U) = {s : o_s associative},  s = 0..14, n,m,p <= 40")
FOCUS = [S for S in SYSTEMS if S.name.split()[0] in
         ("Zeckendorf", "Lucas", "Fib-shift", "Pell")]
N3 = 40
for S in FOCUS:
    G = []
    gre = {n: S.greedy(n) for n in range(1, N3+1)}
    for s in range(15):
        def c(n, m, gn=None, gm=None):
            gn = gre.get(n) or S.greedy(n); gm = gre.get(m) or S.greedy(m)
            return sum(S.U[i+j+s] for i in gn for j in gm)
        ok = True
        for n in range(1, N3+1):
            for m in range(1, N3+1):
                nm = c(n, m)
                for p in range(1, N3+1):
                    if c(nm, p) != c(n, c(m, p)): ok = False; break
                if not ok: break
            if not ok: break
        if ok: G.append(s)
    print(f"   {S.name:34s} G = {G}")

hdr("4.  EXPLICIT MINIMAL COUNTEREXAMPLES (rigorous finite certificates)")
def minimal_counterexample(S, s, N=60):
    for n in range(1, N):
        for m in range(1, N):
            for p in range(1, N):
                a = S.circ(S.circ(n, m, s), p, s); b = S.circ(n, S.circ(m, p, s), s)
                if a != b: return (n, m, p, a, b)
    return None
for name, s in [("Lucas", 3), ("Fib-shift", 5), ("Zeckendorf (Fibonacci)", 1)]:
    S = next(x for x in SYSTEMS if x.name.startswith(name))
    r = minimal_counterexample(S, s)
    print(f"   {S.name:34s} s={s}: (n,m,p)=({r[0]},{r[1]},{r[2]}) -> "
          f"(n o m) o p = {r[3]} != {r[4]} = n o (m o p)")
    check(r is not None, f"expected a counterexample for {name} s={s}")

hdr("5.  Knuth's theorem (control): Zeckendorf, s=2, is associative;  n o_2 m for small n,m")
Z = next(x for x in SYSTEMS if x.name.startswith("Zeckendorf"))
ok = all(Z.circ(Z.circ(n, m, 2), p, 2) == Z.circ(n, Z.circ(m, p, 2), 2)
         for n in range(1, 46) for m in range(1, 46) for p in range(1, 46))
print(f"   associativity of Knuth's o_2 for n,m,p <= 45: {ok}")
check(ok, "Knuth's theorem failed(!)")
print("   table of n o_2 m:")
print("        " + "".join(f"{m:6d}" for m in range(1, 9)))
for n in range(1, 9):
    print(f"   {n:3d} |" + "".join(f"{Z.circ(n,m,2):6d}" for m in range(1, 9)))

hdr("SUMMARY")
print(f"   failures: {len(fails)}")
for f in fails: print("     -", f)
print("   " + ("ALL CHECKS PASSED" if not fails else "SOME CHECKS FAILED"))

hdr("6.  LONG-RANGE SCAN for the two 'never associative' systems (s = 0..20, n,m,p <= 12)")
for name in ("x^2=x+3", "x^4=x^3+1"):
    S = next(x for x in SYSTEMS if x.name.startswith(name))
    G = []
    for s in range(21):
        ok = all(S.circ(S.circ(n, m, s), p, s) == S.circ(n, S.circ(m, p, s), s)
                 for n in range(1, 13) for m in range(1, 13) for p in range(1, 13))
        if ok: G.append(s)
    print(f"   {S.name:34s} G ∩ [0,20] = {G}")
print("\n   (x^2=x+3: beta=2.302.. is NOT Pisot, conjugate -1.302..)")
print("   (x^4=x^3+1: beta=1.380.. IS Pisot, but lacks the finiteness property (F))")

hdr("7.  LEMMA: closure at s  =>  beta^s S subset S  =>  closure at every multiple ks")
for S in SYSTEMS:
    row = []
    for s in range(1, 8):
        clo = all(S.in_S(S.R.mul(S.pw[s], S.R.mul(S.Phi(n), S.Phi(m))))
                  for n in range(1, 21) for m in range(1, 21))
        if not clo: continue
        shift_ok = all(S.in_S(S.R.mul(S.pw[s], S.Phi(n))) for n in range(1, 60))
        check(shift_ok, f"closure at s={s} but beta^s S not in S, {S.name}")
        mult_ok = all(all(S.in_S(S.R.mul(S.pw[k*s], S.R.mul(S.Phi(n), S.Phi(m))))
                          for n in range(1, 16) for m in range(1, 16)) for k in (2, 3))
        check(mult_ok, f"closure at s={s} but not at multiples, {S.name}")
        row.append(s)
    print(f"   {S.name:34s} closure shifts <=7: {row}  (lemma holds)")

hdr("8.  Lucas system: 2 in S but 2*phi not in S   (so phi*S is NOT inside S)")
Lu = next(x for x in SYSTEMS if x.name.startswith("Lucas"))
two   = Lu.Phi(2)                       # 2 = u_0+u_0 -> 2*beta^0
twophi = Lu.R.mul(Lu.pw[1], two)
print(f"   Phi(2)      = {two}    in S: {Lu.in_S(two)}")
print(f"   beta*Phi(2) = {twophi}    in S: {Lu.in_S(twophi)}   (L of it = {Lu.L(twophi)}, "
      f"Phi of that = {Lu.Phi(Lu.L(twophi))})")
check(Lu.in_S(two) and not Lu.in_S(twophi), "expected 2 in S, 2*beta not in S")

hdr("9.  L(z) = Tr(lambda z):  Zeckendorf <-> lambda = phi^2/sqrt5,  Lucas <-> lambda = phi")
from fractions import Fraction as Fr
def mulQ(x, y):                          # (a,b) = a + b*sqrt5
    return (x[0]*y[0] + 5*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
phi = (Fr(1,2), Fr(1,2))
lam_Z = mulQ(mulQ(phi, phi), (0, Fr(1,5)))       # phi^2 / sqrt5 = phi^2 * sqrt5/5
lam_L = phi
for lam, S, tag in ((lam_Z, next(x for x in SYSTEMS if x.name.startswith("Zeckendorf")), "Zeckendorf"),
                    (lam_L, Lu, "Lucas")):
    p = (Fr(1), Fr(0)); ok = True
    for k in range(15):
        if 2*mulQ(lam, p)[0] != S.U[k]: ok = False
        p = mulQ(p, phi)
    print(f"   {tag:12s} Tr(lambda*phi^k) = u_k for k<15: {ok}")
    check(ok, f"trace representation wrong for {tag}")

hdr("FINAL")
print(f"   failures: {len(fails)}")
for f in fails: print("     -", f)
print("   " + ("ALL CHECKS PASSED" if not fails else "SOME CHECKS FAILED"))
