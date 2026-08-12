"""Multiplier cyclic sieving on necklaces: exact test for composite n."""
from math import gcd, comb
from canon import qbinom, qint, polydivexact
from csp import eval_at_root
import sys

def order(u,n):
    o=1; x=u%n
    while x!=1: x=x*u%n; o+=1
    return o

def rot(mask,n,j):
    return ((mask<<j)|(mask>>(n-j)))&((1<<n)-1)

def canon(mask,n):
    return min(rot(mask,n,j) for j in range(n))

def subsets_masks(n,k):
    from itertools import combinations
    for S in combinations(range(n),k):
        m=0
        for x in S: m|=1<<x
        yield m

def multmap(mask,n,u):
    out=0
    for x in range(n):
        if mask>>x & 1: out |= 1<<((u*x)%n)
    return out

def test(n,k,u,verbose=False):
    N=order(u,n)
    reps={}
    for m in subsets_masks(n,k):
        c=canon(m,n)
        reps[c]=1
    R=list(reps)
    fix=[]
    for d in range(N):
        v=pow(u,d,n)
        c=0
        for m in R:
            if canon(multmap(m,n,v),n)==m: c+=1
        fix.append(c)
    X=polydivexact(qbinom(n,k), qint(n))
    co={i:cc for i,cc in enumerate(X) if cc}
    vals=[eval_at_root(co,N,d) for d in range(N)]
    ok=all(v==f for v,f in zip(vals,fix))
    return ok,N,fix,vals,len(R)

if __name__=='__main__':
    lo,hi=int(sys.argv[1]),int(sys.argv[2])
    for n in range(lo,hi+1):
        units=[u for u in range(1,n) if gcd(u,n)==1]
        for k in range(1,n):
            if gcd(n,k)!=1: continue
            if k> n//2: continue
            for u in units:
                if u==1: continue
                ok,N,fix,vals,sz=test(n,k,u)
                if not ok:
                    print(f"FAIL n={n} k={k} u={u} ord={N} |Y|={sz}\n     Fix={fix}\n     Y  ={vals}")
        print(f"  ...n={n} done", flush=True)
