from itertools import combinations
from canon import qbinom, qint, polydivexact, reduce_mod_qn1
from csp import eval_at_root

def primroot(p):
    for g in range(2,p):
        s=set(); x=1
        for _ in range(p-1): x=x*g%p; s.add(x)
        if len(s)==p-1: return g

for p in [5,7,11,13,17,19]:
    g=primroot(p); N=p-1
    print(f"\n=== p={p}, generator g={g}, group order {N} ===")
    for k in range(1,p):
        # translation classes of k-subsets
        def canon(A): return min(tuple(sorted((a+t)%p for a in A)) for t in range(p))
        classes=set()
        for A in combinations(range(p),k): classes.add(canon(A))
        classes=sorted(classes)
        act=lambda A: canon(tuple(sorted((g*a)%p for a in A)))
        # fixed points of act^d
        fix=[]
        for d in range(N):
            c=0
            for A in classes:
                x=A
                for _ in range(d): x=act(x)
                if x==A: c+=1
            fix.append(c)
        # candidate X(q) = qbinom(p,k)/[p]_q
        try:
            X=polydivexact(qbinom(p,k), qint(p))
        except AssertionError:
            print(f"  k={k}: [p]_q does not divide"); continue
        co={}
        for i,c in enumerate(X):
            if c: co[i]=co.get(i,0)+c
        vals=[eval_at_root(co,N,d) for d in range(N)]
        ok = all(v==f for v,f in zip(vals,fix))
        print(f"  k={k}: |S|={len(classes)}  Fix={fix}  X(zeta^d)={vals}  {'*** CSP ***' if ok else 'no'}")
