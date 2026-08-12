from math import gcd
from fast import order
from canon import qbinom, qint, polydivexact
from csp import eval_at_root
from collections import defaultdict
import sys

Yc={}
def Yat(n,k,m):
    if (n,k) not in Yc: Yc[(n,k)]=polydivexact(qbinom(n,k), qint(n))
    X=Yc[(n,k)]; co={i:c for i,c in enumerate(X) if c}
    return eval_at_root(co,m,1)

def cyc(v,c,n):
    seen=bytearray(n); s=[]
    for y in range(n):
        if seen[y]: continue
        l=0;z=y
        while not seen[z]: seen[z]=1; z=(v*z+c)%n; l+=1
        s.append(l)
    return s
def Nk(sizes,k):
    dp=[0]*(k+1); dp[0]=1
    for s in sizes:
        if s<=k:
            for j in range(k,s-1,-1): dp[j]+=dp[j-s]
    return dp[k]
def fixY(n,k,v):
    tot=sum(Nk(cyc(v,c,n),k) for c in range(n)); assert tot%n==0
    return tot//n

res=defaultdict(lambda:[0,0,[]])
NM=int(sys.argv[1])
for n in range(4,NM+1):
    for v in range(2,n):
        if gcd(v,n)!=1: continue
        m=order(v,n)
        sz=cyc(v,0,n)
        if not all(x in (1,m) for x in sz): continue
        f=sz.count(1)
        allok=True
        for k in range(1,n):
            if gcd(n,k)!=1: continue
            if Yat(n,k,m)!=fixY(n,k,v): allok=False; break
        r=res[(m,f)]
        r[0 if allok else 1]+=1
        if len(r[2])<3: r[2].append((n,'ok' if allok else 'BAD'))
print(f"{'m':>3}{'f':>4}  hold fail   examples")
for key in sorted(res):
    h,fl,ex=res[key]
    mark = "HOLD" if fl==0 else ("FAIL" if h==0 else "***MIXED***")
    print(f"{key[0]:>3}{key[1]:>4}  {h:>4} {fl:>4}   {mark:11s} {ex}")
