from math import gcd
from itertools import combinations

def orbits(objs, act):
    seen=set(); orbs=[]
    for o in objs:
        if o in seen: continue
        orb=[]; x=o
        while x not in seen:
            seen.add(x); orb.append(x); x=act(x)
        orbs.append(orb)
    return orbs

def canonical_dist(objs, act, n):
    """Return length-n vector: canonical CSP polynomial coefficients mod q^n-1."""
    v=[0]*n
    for orb in orbits(objs,act):
        m=len(orb)
        for j in range(m):
            v[(j*(n//m))%n]+=1
    return v

def stat_dist(objs, stat, n):
    v=[0]*n
    for o in objs: v[stat(o)%n]+=1
    return v

def qint(n):
    return [1]*n  # [n]_q = 1+q+...+q^{n-1}

def polymul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): r[i+j]+=x*y
    return r

def polydivexact(a,b):
    a=a[:]; 
    while a and a[-1]==0: a.pop()
    q=[0]*(len(a)-len(b)+1)
    for i in range(len(a)-len(b),-1,-1):
        c=a[i+len(b)-1]//b[-1]
        q[i]=c
        for j in range(len(b)): a[i+j]-=c*b[j]
    assert all(x==0 for x in a), "not exact"
    return q

def qbinom(n,k):
    if k<0 or k>n: return [0]
    num=[1]
    for i in range(n-k+1,n+1): num=polymul(num,qint(i))
    den=[1]
    for i in range(1,k+1): den=polymul(den,qint(i))
    return polydivexact(num,den)

def reduce_mod_qn1(p,n):
    v=[0]*n
    for i,c in enumerate(p): v[i%n]+=c
    return v

