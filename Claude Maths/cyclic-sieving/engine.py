from itertools import product, combinations
from math import gcd
from canon import canonical_dist

def dist_of(objs, stat, N):
    v=[0]*N
    for o in objs: v[stat(o)%N]+=1
    return v

def test(objs, act, N, stats, tag, out):
    can = canonical_dist(objs, act, N)
    for sname, s in stats.items():
        d = dist_of(objs, s, N)
        if d==can:
            out.append((tag, sname, 'exact'))
        else:
            for sh in range(N):
                if [d[(i-sh)%N] for i in range(N)]==can:
                    out.append((tag, sname, f'shift{sh}')); break

def STATS(n, r=2):
    def sumpos(w): return sum(i*w[i] for i in range(len(w)))
    def sumpos1(w): return sum((i+1)*w[i] for i in range(len(w)))
    def wsum(w): return sum(w)
    def maj(w): return sum(i+1 for i in range(len(w)-1) if w[i]>w[i+1])
    def comaj(w): return sum(len(w)-1-i for i in range(len(w)-1) if w[i]>w[i+1])
    def inv(w): return sum(1 for i in range(len(w)) for j in range(i+1,len(w)) if w[i]>w[j])
    def cycdes(w): return sum(1 for i in range(len(w)) if w[i]>w[(i+1)%len(w)])
    def sumval(w): return sum(i*w[i] for i in range(len(w)))
    def sumall(w): return sum((i+1)*w[i] for i in range(len(w)))
    return {'sumpos':sumpos,'sumpos1':sumpos1,'maj':maj,'comaj':comaj,'inv':inv,
            'wsum':wsum,'cycdes':cycdes}

out=[]
# F1/F3: all words over Z_r, rotation
for r in (2,3):
  for n in range(2,9):
    objs=[tuple(w) for w in product(range(r),repeat=n)]
    act=lambda w:tuple(w[(i+1)%len(w)] for i in range(len(w)))
    test(objs,act,n,STATS(n,r),f'words r={r} n={n} rot',out)

# F7: binary words, twisted rotation (complement twist), order 2n
for n in range(2,9):
    objs=[tuple(w) for w in product((0,1),repeat=n)]
    act=lambda w:tuple(1-w[(i+1)%len(w)] for i in range(len(w)))
    test(objs,act,2*n,STATS(n),f'binary n={n} TWIST-complement',out)

# F9: words over Z_r with affine twist tau(w)_i = w_{i+1}+1
from math import lcm
for r in (2,3,4):
  for n in range(2,8):
    objs=[tuple(w) for w in product(range(r),repeat=n)]
    act=lambda w:tuple((w[(i+1)%len(w)]+1)%r for i in range(len(w)))
    test(objs,act,lcm(n,r),STATS(n,r),f'words r={r} n={n} AFFINE-TWIST',out)

for t in out: print(t)
print("hits:",len(out))
