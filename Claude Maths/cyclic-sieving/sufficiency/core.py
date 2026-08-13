"""Exact machinery: cyclotomics, Y(zeta_m), multiplier fixed-point counts."""
from math import gcd
from functools import lru_cache

# ---------- polynomials over Z, as tuples of coefficients (low -> high) ----------
def ptrim(a):
    a=list(a)
    while a and a[-1]==0: a.pop()
    return tuple(a)

def pmul(a,b):
    if not a or not b: return ()
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): r[i+j]+=x*y
    return ptrim(r)

def pdivexact(a,b):
    """exact division a/b (b monic-ish, leading coeff divides)"""
    a=list(a); b=list(b)
    q=[0]*(len(a)-len(b)+1)
    for i in range(len(a)-len(b), -1, -1):
        c,rem = divmod(a[i+len(b)-1], b[-1])
        assert rem==0, "not exact"
        q[i]=c
        for j,y in enumerate(b): a[i+j]-=c*y
    assert all(x==0 for x in a), "remainder nonzero"
    return ptrim(q)

def pmod(a,b):
    """a mod b, b monic"""
    a=list(a); b=list(b)
    assert b[-1]==1
    for i in range(len(a)-len(b), -1, -1):
        c=a[i+len(b)-1]
        if c:
            for j,y in enumerate(b): a[i+j]-=c*y
    return ptrim(a)

@lru_cache(maxsize=None)
def cyclotomic(d):
    num=tuple([-1]+[0]*(d-1)+[1])          # x^d - 1
    for e in range(1,d):
        if d%e==0: num=pdivexact(num, cyclotomic(e))
    return num

@lru_cache(maxsize=None)
def divisors(n): return tuple(d for d in range(1,n+1) if n%d==0)

# ---------- Y_{n,k}(q) = qbin(n,k)/[n]_q  =  prod_{d in E} Phi_d ----------
def Y_at_zeta(n,k,m):
    """Y_{n,k}(zeta_m) as element of Z[x]/Phi_m ; returns tuple of coeffs."""
    Pm=cyclotomic(m)
    acc=(1,)
    for d in range(2, n+1):
        if n%d and (k%d) > (n%d):
            acc=pmod(pmul(acc, pmod(cyclotomic(d),Pm)), Pm)
    return acc

def as_rational_int(poly, m):
    """If element of Z[x]/Phi_m is a rational integer, return it, else None."""
    p=ptrim(poly)
    if len(p)<=1:
        return p[0] if p else 0
    return None

# ---------- multiplier combinatorics ----------
def orbit_sizes(v,n):
    seen=[False]*n; sizes=[]
    for y in range(n):
        if not seen[y]:
            L=0; z=y
            while not seen[z]:
                seen[z]=True; z=(v*z)%n; L+=1
            sizes.append(L)
    return sizes

def cycle_sizes_affine(v,c,n):
    seen=[False]*n; sizes=[]
    for y in range(n):
        if not seen[y]:
            L=0; z=y
            while not seen[z]:
                seen[z]=True; z=(v*z+c)%n; L+=1
            sizes.append(L)
    return sizes

def coeff_subset_sum(sizes, k):
    """[z^k] prod (1+z^L)"""
    dp=[0]*(k+1); dp[0]=1
    for L in sizes:
        if L<=k:
            for i in range(k,L-1,-1):
                if dp[i-L]: dp[i]+=dp[i-L]
    return dp[k]

def order_mod(v,n):
    m=1; x=v%n
    while x!=1:
        x=(x*v)%n; m+=1
    return m

def K_set(v,n,m):
    S=sum(pow(v,i,n) for i in range(m))%n
    return [c for c in range(n) if (c*S)%n==0]

def fixed_count(v,n,k):
    m=order_mod(v,n)
    tot=0
    for c in K_set(v,n,m):
        tot+=coeff_subset_sum(cycle_sizes_affine(v,c,n), k)
    assert tot%n==0
    return tot//n

def rad(x):
    r=1; d=2; y=x
    while d*d<=y:
        if y%d==0:
            r*=d
            while y%d==0: y//=d
        d+=1
    if y>1: r*=y
    return r

def criterion(v,n,m):
    for e in divisors(m):
        if e==m: continue
        fe=gcd(pow(v,e,n)-1 if pow(v,e,n)>=1 else n, n)
        fe=gcd((pow(v,e,n)-1)%n, n)
        if fe==0: fe=n
        if fe not in (1,2,3,4,6): return False
        if m%rad(fe)!=0: return False
    return True
