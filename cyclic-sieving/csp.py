"""General cyclic-sieving detector.
Given a finite set S, a cyclic action c of order n, and a statistic stat,
test whether (S, C_n, X(q)=sum q^stat) exhibits CSP:
   X(zeta_n^d) = #{s : c^d(s)=s}   for all d.
Exact arithmetic: work in Z[q]/(q^n-1), evaluate at zeta^d by reducing mod cyclotomic.
"""
from fractions import Fraction

def poly_from_stats(stats, n=None):
    """stats: list of ints -> coefficient dict"""
    c = {}
    for s in stats: c[s] = c.get(s,0)+1
    return c

def cyclotomic(n, cache={}):
    """Phi_n(x) as coeff list (low->high), integer coefficients."""
    if n in cache: return cache[n]
    # x^n - 1 = prod_{d|n} Phi_d
    num = [-1]+[0]*(n-1)+[1]
    for d in range(1,n):
        if n%d==0:
            num = polydiv(num, cyclotomic(d))
    cache[n]=num
    return num

def polydiv(a,b):
    a=a[:]; q=[0]*(len(a)-len(b)+1)
    for i in range(len(a)-len(b), -1, -1):
        coef = a[i+len(b)-1]//b[-1]
        q[i]=coef
        for j in range(len(b)):
            a[i+j]-=coef*b[j]
    assert all(x==0 for x in a), (a,b)
    return q

def polymod(a, m):
    a=a[:]
    while len(a)>=len(m):
        if a[-1]==0: a.pop(); continue
        coef=a[-1]//m[-1]
        assert a[-1]==coef*m[-1]
        off=len(a)-len(m)
        for j in range(len(m)):
            a[off+j]-=coef*m[j]
        while a and a[-1]==0: a.pop()
    return a

def eval_at_root(coeffs, n, d):
    """coeffs: dict exp->coef. Evaluate at zeta_n^d.  zeta_n^d is a primitive m-th root, m=n/gcd(n,d).
    Returns an integer if the value is rational, else None (value in Z[zeta_m])."""
    from math import gcd
    m = n//gcd(n,d)
    # reduce exponents: q^e -> q^(e*d mod n) then as primitive m-th root: exponent*d/gcd ... 
    # simpler: substitute q = w where w^m=1 primitive. exponent e -> e*d mod n; but zeta_n^(e*d) = zeta_m^(e*d/gcd(n,d))
    g = gcd(n,d)
    red = [0]*m
    for e,cf in coeffs.items():
        ee = (e*d) % n
        assert ee % g == 0
        red[(ee//g) % m] += cf
    # now value = sum red[i] zeta_m^i ; reduce mod Phi_m
    r = polymod(red, cyclotomic(m))
    while r and r[-1]==0: r.pop()
    if len(r)<=1:
        return r[0] if r else 0
    return ('nonrat', r)

def test_csp(objs, act, stat, n, verbose=False):
    """objs: list; act: function giving c(s); stat: dict obj->int or function"""
    idx = {o:i for i,o in enumerate(objs)}
    # build permutation
    perm = [idx[act(o)] for o in objs]
    coeffs = {}
    for o in objs:
        s = stat(o) if callable(stat) else stat[o]
        coeffs[s] = coeffs.get(s,0)+1
    ok = True; details=[]
    for d in range(n):
        # fixed points of c^d
        cnt = 0
        for i in range(len(objs)):
            j=i
            for _ in range(d): j=perm[j]
            if j==i: cnt+=1
        v = eval_at_root(coeffs, n, d)
        good = (v == cnt)
        details.append((d, v, cnt, good))
        if not good: ok=False
    return ok, details, coeffs
