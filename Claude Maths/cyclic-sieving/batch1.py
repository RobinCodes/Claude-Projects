from itertools import product
from collections import Counter
from math import gcd

def necklace_reps(n):
    seen=set(); reps=[]
    for bits in product((0,1),repeat=n):
        if bits in seen: continue
        orb={tuple(bits[(i+j)%n] for i in range(n)) for j in range(n)}
        seen|=orb
        reps.append(bits)
    return reps

def autocorr(w):
    n=len(w)
    return tuple(sum(w[i]*w[(i+d)%n] for i in range(n)) for d in range(n))

S={}
NMAX=18
S['necklaces']=[]
S['distinct_autocorr']=[]
S['homometric_classes_gt1']=[]
S['max_homometry_class']=[]
for n in range(1,NMAX+1):
    reps=necklace_reps(n)
    S['necklaces'].append(len(reps))
    c=Counter(autocorr(w) for w in reps)
    S['distinct_autocorr'].append(len(c))
    S['homometric_classes_gt1'].append(sum(1 for v in c.values() if v>1))
    S['max_homometry_class'].append(max(c.values()))
for k,v in S.items(): print(k, v)
