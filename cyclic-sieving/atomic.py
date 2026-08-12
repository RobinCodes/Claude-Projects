from math import gcd, comb
from fast import order, fixY
from canon import qbinom, qint, polydivexact
from csp import eval_at_root
from collections import Counter, defaultdict

Yc={}
def Yat(n,k,m):
    if (n,k) not in Yc:
        Yc[(n,k)]=polydivexact(qbinom(n,k), qint(n))
    X=Yc[(n,k)]; co={i:c for i,c in enumerate(X) if c}
    return eval_at_root(co,m,1)

def orbits_of(v,n):
    seen=[False]*n; sizes=[]
    for y in range(n):
        if seen[y]: continue
        l=0;z=y
        while not seen[z]: seen[z]=True; z=v*z%n; l+=1
        sizes.append(l)
    return tuple(sorted(sizes))

rows=[]
NMAX=34
for n in range(4,NMAX+1):
    for v in range(2,n):
        if gcd(v,n)!=1: continue
        m=order(v,n)
        orb=orbits_of(v,n)
        uniform = all(x in (1,m) for x in orb)
        f=orb.count(1)
        assert f==gcd(v-1,n)
        for k in range(1,n//2+1):
            if gcd(n,k)!=1: continue
            y=Yat(n,k,m); fx=fixY(n,k,v)
            rows.append(dict(n=n,k=k,v=v,m=m,r=n%m,s=k%m,f=f,uniform=uniform,
                             orb=orb, ok=(y==fx), y=y, fx=fx))
print("total atomic instances:", len(rows))
ok=[r for r in rows if r['ok']]; bad=[r for r in rows if not r['ok']]
print("hold:",len(ok)," fail:",len(bad))

# cross-tab by (r = n mod m, s = k mod m)
tab=defaultdict(lambda:[0,0])
for r in rows:
    tab[(r['r'],r['s'])][0 if r['ok'] else 1]+=1
print("\n (n mod m, k mod m) -> [hold, fail]")
for key in sorted(tab): 
    h,fl=tab[key]
    if fl or h: print("   ",key,tab[key])
