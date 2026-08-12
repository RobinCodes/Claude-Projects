from math import gcd
from lib2 import *
from crit import Yat_fast
import sys
lo,hi=int(sys.argv[1]),int(sys.argv[2])
mis=0;tested=0;hold=0
examples=[]
for n in range(lo,hi+1):
    for v in range(2,n):
        if gcd(v,n)!=1: continue
        m=order(v,n); f=gcd(v-1,n)
        ok=True
        for k in range(1,n):
            if gcd(n,k)!=1: continue
            if Yat_fast(n,k,m)!=fixY_alpha(n,[k,n-k],v,m): ok=False; break
        tested+=1; hold+=ok
        if ok!=predict(m,f):
            mis+=1
            if len(examples)<12: examples.append((n,v,m,f,ok,predict(m,f)))
    if n%10==0: print(f"  n<={n}: tested={tested} hold={hold} mismatch={mis}",flush=True)
print(f"RANGE {lo}..{hi}: tested={tested} hold={hold} MISMATCHES={mis}")
for e in examples: print("   n,v,m,f,truth,pred =",e)
