from math import gcd
from lib2 import *
from crit import Yat_fast
holds=[]; tested=0
for n in range(4,61):
    for v in range(2,n):
        if gcd(v,n)!=1: continue
        m=order(v,n); sz=cyc(v,0,n)
        if all(x in (1,m) for x in sz): continue
        tested+=1; ok=True
        for k in range(1,n):
            if gcd(n,k)!=1: continue
            if Yat_fast(n,k,m)!=fixY_alpha(n,[k,n-k],v,m): ok=False; break
        if ok: holds.append((n,v,m,sorted(set(sz)),gcd(v-1,n)))
print(f"non-uniform multipliers tested: {tested}")
print(f"non-uniform that SATISFY the sieving: {len(holds)}")
for x in holds[:20]: print("   n,v,m,orbit-sizes,f =",x)
