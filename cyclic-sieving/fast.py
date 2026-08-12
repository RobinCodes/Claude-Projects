"""Exact, non-enumerative test of the multiplier CSP on binary necklaces."""
from math import gcd, comb
from canon import qbinom, qint, polydivexact
from csp import eval_at_root, cyclotomic, polymod

def order(u,n):
    o=1; x=u%n
    while x!=1: x=x*u%n; o+=1
    return o

def cycletype(v,c,n):
    seen=[False]*n; sizes=[]
    for y in range(n):
        if seen[y]: continue
        l=0; z=y
        while not seen[z]:
            seen[z]=True; z=(v*z+c)%n; l+=1
        sizes.append(l)
    return sizes

def Nk(sizes,k):
    dp=[0]*(k+1); dp[0]=1
    for s in sizes:
        for j in range(k,s-1,-1):
            dp[j]+=dp[j-s]
    return dp[k]

def fixY(n,k,v):
    tot=0
    for c in range(n):
        tot+=Nk(cycletype(v,c,n),k)
    assert tot%n==0, (n,k,v,tot)
    return tot//n

def Yval(n,k,m):
    """Y(q)=qbinom(n,k)/[n]_q evaluated at a primitive m-th root of unity, exactly.
       Returns int or ('nonrat', coeffs)."""
    X=polydivexact(qbinom(n,k), qint(n))
    co={i:c for i,c in enumerate(X) if c}
    return eval_at_root(co, m, 1)

def brute_check(n, verbose=False):
    """for each unit u and each k coprime to n, test CSP over all powers"""
    units=[u for u in range(1,n) if gcd(u,n)==1]
    res={}
    for u in units:
        if u==1: continue
        N=order(u,n)
        ok=True; detail=[]
        for k in range(1,n):
            if gcd(n,k)!=1: continue
            for d in range(N):
                v=pow(u,d,n); m=order(v,n) if v!=1 else 1
                f=fixY(n,k,v)
                y=Yval(n,k,m) if m>1 else comb(n,k)//n
                if y!=f:
                    ok=False; detail.append((k,d,m,y,f))
        res[u]=(N,ok,detail)
    return res

if __name__=='__main__':
    import sys
    for n in range(4, int(sys.argv[1])+1):
        r=brute_check(n)
        good=[u for u,(N,ok,_) in r.items() if ok]
        bad=[u for u,(N,ok,_) in r.items() if not ok]
        print(f"n={n:3d}  PASS u={good}   FAIL u={bad}", flush=True)
