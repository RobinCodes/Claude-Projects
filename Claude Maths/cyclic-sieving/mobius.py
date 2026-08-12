from math import gcd
from canon import canonical_dist, reduce_mod_qn1, polymul, qint

def words(n,q=2):
    for i in range(q**n):
        w=[]
        x=i
        for _ in range(n): w.append(x%q); x//=q
        yield tuple(w)

def twisted_rot(w):  # rotate then complement
    n=len(w)
    r=tuple(w[(i+1)%n] for i in range(n))
    return tuple(1-x for x in r)

print("=== Mobius (twisted) necklaces of binary words: canonical 2n-th-root polynomial ===")
for n in range(2,11):
    objs=list(words(n))
    N=2*n
    can=canonical_dist(objs, twisted_rot, N)
    # candidate products
    cands={}
    cands['prod(1+q^(2i+1)), i=0..n-1'] = None
    p=[1]
    for i in range(n):
        p=polymul(p,[1]+[0]*(2*i)+[1])
    cands['prod(1+q^(2i+1)), i=0..n-1']=reduce_mod_qn1(p,N)
    p2=[1]
    for i in range(n):
        p2=polymul(p2,[1]+[0]*(2*i+1)+[1])   # 1+q^{2i+2}
    cands['prod(1+q^(2i+2))']=reduce_mod_qn1(p2,N)
    p3=[1]
    for i in range(n):
        p3=polymul(p3,[1]+[0]*(i)+[1])   # 1+q^{i+1}
    cands['prod(1+q^(i+1))']=reduce_mod_qn1(p3,N)
    out=[]
    for name,v in cands.items():
        if v==can: out.append(name+" EXACT")
        else:
            for s in range(N):
                if [v[(i-s)%N] for i in range(N)]==can: out.append(name+f" shift{s}"); break
    print(f"n={n} canon={can}  ->  {out}")
