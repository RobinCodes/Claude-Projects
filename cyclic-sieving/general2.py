from math import gcd
from lib2 import *
from crit import Yat_fast
import sys

def predict2(n,v,m):
    for e in range(1,m):
        if m%e: continue
        fe=gcd(pow(v,e,n)-1,n)
        if phi(fe)>2 or m%rad(fe)!=0: return False
    return True

lo,hi=int(sys.argv[1]),int(sys.argv[2])
mis=0;tested=0;hold=0;ex=[]
for n in range(lo,hi+1):
    for v in range(2,n):
        if gcd(v,n)!=1: continue
        m=order(v,n)
        ok=True
        for k in range(1,n):
            if gcd(n,k)!=1: continue
            if Yat_fast(n,k,m)!=fixY_alpha(n,[k,n-k],v,m): ok=False; break
        tested+=1; hold+=ok
        if ok!=predict2(n,v,m):
            mis+=1
            if len(ex)<10: ex.append((n,v,m,gcd(v-1,n),ok,predict2(n,v,m)))
    if n%12==0: print(f"  n<={n}: tested={tested} hold={hold} mismatch={mis}",flush=True)
print(f"RANGE {lo}..{hi}: tested={tested} hold={hold} MISMATCHES={mis}")
for e in ex: print("   n,v,m,f,truth,pred =",e)
