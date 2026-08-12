from math import gcd, comb
from fast import order, fixY, Yval
from canon import qbinom, qint, polydivexact
from csp import eval_at_root
import sys

Ycache={}
def Ypoly(n,k):
    if (n,k) not in Ycache:
        Ycache[(n,k)]=polydivexact(qbinom(n,k), qint(n))
    return Ycache[(n,k)]

def Yat(n,k,m):
    X=Ypoly(n,k); co={i:c for i,c in enumerate(X) if c}
    return eval_at_root(co,m,1)

def status(n,u):
    N=order(u,n); ok=True
    for k in range(1,n):
        if gcd(n,k)!=1: continue
        for d in range(1,N):
            v=pow(u,d,n); m=order(v,n)
            if Yat(n,k,m)!=fixY(n,k,v): return False,N
    return True,N

for n in range(4,41):
    units=[u for u in range(2,n) if gcd(u,n)==1]
    for u in units:
        ok,N=status(n,u)
        prof=sorted({gcd(pow(u,d,n)-1,n) for d in range(1,N)})
        print(f"n={n:3d} u={u:3d} ord={N:2d} {'PASS' if ok else 'FAIL'}  fixedset-sizes={prof}", flush=True)
