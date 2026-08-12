from math import gcd
from fast import order, fixY
from canon import qbinom, qint, polydivexact
from csp import eval_at_root
from collections import defaultdict

Yc={}
def Yat(n,k,m):
    if (n,k) not in Yc: Yc[(n,k)]=polydivexact(qbinom(n,k), qint(n))
    X=Yc[(n,k)]; co={i:c for i,c in enumerate(X) if c}
    return eval_at_root(co,m,1)

def orbsizes(v,n):
    seen=[False]*n; s=[]
    for y in range(n):
        if seen[y]: continue
        l=0;z=y
        while not seen[z]: seen[z]=True; z=v*z%n; l+=1
        s.append(l)
    return s

res=defaultdict(set)
for n in range(4,49):
    for v in range(2,n):
        if gcd(v,n)!=1: continue
        m=order(v,n); sz=orbsizes(v,n)
        if not all(x in (1,m) for x in sz): continue      # uniform only
        f=sz.count(1); t=sz.count(m)
        allok=True
        for k in range(1,n):
            if gcd(n,k)!=1: continue
            if Yat(n,k,m)!=fixY(n,k,v): allok=False; break
        res[(m,f)].add((n,t,allok))
print("uniform multipliers: cycle type (m^t, 1^f)")
print(f"{'m':>3} {'f':>3}  verdict           examples (n,t)")
for (m,f) in sorted(res):
    S=res[(m,f)]
    ok=[x for x in S if x[2]]; bad=[x for x in S if not x[2]]
    verdict = "ALWAYS HOLDS" if not bad else ("ALWAYS FAILS" if not ok else "MIXED")
    ex = sorted([(x[0],x[1]) for x in (ok if ok else bad)])[:4]
    print(f"{m:>3} {f:>3}  {verdict:14s}  ok={len(ok):2d} bad={len(bad):2d}  {ex}")
