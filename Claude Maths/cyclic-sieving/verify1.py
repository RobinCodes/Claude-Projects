from math import comb, gcd
from canon import qbinom, qint, polydivexact
from csp import eval_at_root

def check(n, kmax=None, group=None, label=""):
    """n: modulus. group: list of units generating cyclic Gamma (we use a generator).
       Verify X(q)=qbinom(n,k)/[n]_q  sieves multiplier action on k-subsets mod translation."""
    units=[u for u in range(1,n) if gcd(u,n)==1]
    # find generator of (Z/n)^*
    gen=None
    for g in units:
        s=set(); x=1
        for _ in range(len(units)): x=x*g%n; s.add(x)
        if len(s)==len(units): gen=g;break
    if gen is None: return None,"(Z/%d)* not cyclic"%n
    N=len(units)
    bad=[]
    for k in range(1,n):
        if gcd(n,k)!=1: continue
        X=polydivexact(qbinom(n,k), qint(n))
        co={i:c for i,c in enumerate(X) if c}
        # predicted fixed points via orbit-union formula
        for d in range(N):
            s_=gcd(d,N) if d>0 else N
            m=N//s_ if d>0 else 1
            h=pow(gen,d,n)
            # brute force fixed count on translation classes
            # A subset fixed iff h*A = A + t  <=> exists t
            cnt=0
            from itertools import combinations
            seen=set()
            for A in combinations(range(n),k):
                c=min(tuple(sorted((a+t)%n for a in A)) for t in range(n))
                if c in seen: continue
                seen.add(c)
                B=tuple(sorted((h*a)%n for a in A))
                cB=min(tuple(sorted((b+t)%n for b in B)) for t in range(n))
                if cB==c: cnt+=1
            v=eval_at_root(co,N,d)
            if v!=cnt: bad.append((k,d,v,cnt))
    return bad, f"n={n} N={N} gen={gen}"

for n in [5,7,9,11,13,14,17,18,19,25,27,22,23]:
    b,info=check(n)
    if b is None: print(info); continue
    print(f"{info}: {'CSP HOLDS for all k coprime to n' if not b else 'FAILS: '+str(b[:6])}")
