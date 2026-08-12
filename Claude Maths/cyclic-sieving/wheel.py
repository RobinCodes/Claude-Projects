from itertools import combinations
from canon import canonical_dist

def spanning_trees(nv, edges):
    """return list of frozenset of edge-indices forming spanning trees"""
    res=[]
    m=len(edges)
    for comb in combinations(range(m), nv-1):
        par=list(range(nv))
        def find(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        ok=True
        for e in comb:
            a,b=edges[e]; ra,rb=find(a),find(b)
            if ra==rb: ok=False;break
            par[ra]=rb
        if ok: res.append(frozenset(comb))
    return res

for n in range(3,10):
    # wheel: hub = n, rim = 0..n-1
    edges=[]
    for i in range(n): edges.append((i,(i+1)%n))
    for i in range(n): edges.append((i,n))
    eidx={}
    for i,e in enumerate(edges): eidx[frozenset(e)]=i
    def rot_edge(i):
        a,b=edges[i]
        a2 = (a+1)%n if a<n else n
        b2 = (b+1)%n if b<n else n
        return eidx[frozenset((a2,b2))]
    T=spanning_trees(n+1, edges)
    Tset=set(T)
    def act(t): return frozenset(rot_edge(i) for i in t)
    assert all(act(t) in Tset for t in T)
    can=canonical_dist(T, act, n)
    # fixed point counts
    fix=[]
    for d in range(n):
        c=0
        for t in T:
            x=t
            for _ in range(d): x=act(x)
            if x==t: c+=1
        fix.append(c)
    print(f"n={n} #trees={len(T)}  Fix(sigma^d)={fix}")
    print(f"      canonical X coeffs mod q^{n}-1: {can}")
