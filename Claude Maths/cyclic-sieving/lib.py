from math import comb, factorial, gcd
from functools import lru_cache

N = 30

def mobius_table(n):
    mu = [1]*(n+1); primes=[]
    is_comp=[False]*(n+1)
    for i in range(2,n+1):
        if not is_comp[i]: primes.append(i); mu[i]=-1
        for p in primes:
            if i*p>n: break
            is_comp[i*p]=True
            if i%p==0: mu[i*p]=0; break
            else: mu[i*p]=-mu[i]
    mu[1]=1
    return mu
MU = mobius_table(200)

def divisors(n):
    d=[]
    i=1
    while i*i<=n:
        if n%i==0:
            d.append(i)
            if i!=n//i: d.append(n//i)
        i+=1
    return sorted(d)

def orbit_transform(a):
    """a is 1-indexed list a[1..N]. returns b[n] = (1/n) sum_{d|n} mu(n/d) a[d], as Fraction-free check."""
    b=[None]*(len(a))
    for n in range(1,len(a)):
        s=0
        for d in divisors(n):
            s += MU[n//d]*a[d]
        b[n] = (s, n)   # keep numerator, denominator
    return b

def classify(a, name):
    """returns (integral_upto, nonneg_upto, b list)"""
    bs=[]
    integral=True; nonneg=True; firstbad=None
    for n in range(1,len(a)):
        s=0
        for d in divisors(n):
            s += MU[n//d]*a[d]
        if s % n != 0:
            integral=False
            if firstbad is None: firstbad=('nonint',n)
            bs.append(None); continue
        v=s//n
        bs.append(v)
        if v<0:
            nonneg=False
            if firstbad is None: firstbad=('neg',n)
    return integral, nonneg, bs, firstbad
