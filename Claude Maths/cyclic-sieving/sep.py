from canon import *
from itertools import combinations
from math import gcd

def cyclic_gaps_ok(s, n, m):
    k=len(s)
    for i in range(k):
        g=(s[(i+1)%k]-s[i])%n
        if i==k-1: g=(s[0]+n-s[-1])%n
        if g<m and not (k==1): return False
    return True

def gaps(s,n):
    k=len(s); return [ (s[(i+1)%k]-s[i])%n if i<k-1 else (s[0]+n-s[-1])%n for i in range(k)]

print("=== m-separated cyclic k-subsets: does [n]/[n-mk] * qbinom(n-mk,k) sieve? ===")
bad=[]
for m in range(1,5):
  for n in range(4,17):
    for k in range(1,n//m+1):
        if n-m*k<k: continue
        objs=[]
        for s in combinations(range(n),k):
            g=gaps(s,n)
            if k==1 or all(x>=m+1 for x in g): objs.append(s)
        if not objs: continue
        act=lambda s: tuple(sorted((x+1)%n for x in s))
        can=canonical_dist(objs,act,n)
        try:
            cand=polydivexact(polymul(qint(n),qbinom(n-m*k,k)),qint(n-m*k))
        except AssertionError:
            bad.append((m,n,k,'noprod')); continue
        if sum(cand)!=len(objs):
            bad.append((m,n,k,'wrongcount',sum(cand),len(objs))); continue
        red=reduce_mod_qn1(cand,n)
        if red!=can:
            sh=None
            for s in range(n):
                if [red[(i-s)%n] for i in range(n)]==can: sh=s;break
            bad.append((m,n,k,'mismatch',sh))
print("failures:",bad[:40], "total tested ok" if not bad else "")
print("count of failures:",len(bad))
