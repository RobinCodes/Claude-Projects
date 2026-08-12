from math import gcd
from crit import Yat_fast
import sys

def order(u,n):
    o=1;x=u%n
    while x!=1: x=x*u%n; o+=1
    return o
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
def fixY(n,k,v,m):
    S=0; p=1
    for _ in range(m): S=(S+p)%n; p=p*v%n
    g=gcd(S,n) if S else n          # |K| = gcd(S,n); if S==0, K = Z_n
    step=n//g
    tot=0
    for c in range(0,n,step): tot+=Nk(cyc(v,c,n),k)
    assert tot%n==0,(n,k,v)
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
def phi(x):
    r=x; d=2; y=x
    while d*d<=y:
        if y%d==0:
            while y%d==0: y//=d
            r-=r//d
        d+=1
    if y>1: r-=r//y
    return r

def predict(m,f): return phi(f)<=2 and m%rad(f)==0

if __name__ == '__main__':
    lo,hi=int(sys.argv[1]),int(sys.argv[2])
    mism=0; tested=0; holds=0
    for n in range(lo,hi+1):
        for v in range(2,n):
            if gcd(v,n)!=1: continue
            m=order(v,n); sz=cyc(v,0,n)
            if not all(x in (1,m) for x in sz): continue
            f=sz.count(1)
            truth=True
            for k in range(1,n):
                if gcd(n,k)!=1: continue
                if Yat_fast(n,k,m)!=fixY(n,k,v,m): truth=False; break
            tested+=1; holds+=truth
            if truth!=predict(m,f):
                mism+=1
                print(f"  MISMATCH n={n} v={v} m={m} f={f} truth={truth} pred={predict(m,f)}",flush=True)
        if n%10==0: print(f"n<={n}: tested={tested} hold={holds} mismatches={mism}",flush=True)
    print(f"DONE range {lo}..{hi}: tested={tested} hold={holds} MISMATCHES={mism}")

