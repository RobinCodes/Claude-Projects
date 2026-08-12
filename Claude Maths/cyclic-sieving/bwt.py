from itertools import product

def bwt(w):
    n=len(w)
    rots=sorted(tuple(w[(i+j)%n] for j in range(n)) for i in range(n))
    return tuple(r[-1] for r in rots)

def canon(w):
    n=len(w)
    return min(tuple(w[(i+j)%n] for j in range(n)) for i in range(n))

def necklaces(n):
    seen=set(); out=[]
    for w in product((0,1),repeat=n):
        c=canon(w)
        if c not in seen: seen.add(c); out.append(c)
    return out

for n in range(2,17):
    N=necklaces(n)
    f={w:canon(bwt(w)) for w in N}
    # functional graph analysis
    color={}; cycles=[]
    for start in N:
        path=[]; x=start
        while x not in color:
            color[x]=('onpath',len(path)); path.append(x); x=f[x]
        if color[x][0]=='onpath' and x in path:
            i=path.index(x); cyc=path[i:]
            cycles.append(len(cyc))
        for y in path: color[y]=('done',0)
    fixed=[w for w in N if f[w]==w]
    from collections import Counter
    print(f"n={n:2d} necklaces={len(N):5d} fixedpts={len(fixed):3d} cycles={Counter(cycles)}")
