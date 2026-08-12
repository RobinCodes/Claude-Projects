from math import gcd
from lib2 import *
from csp import cyclotomic, polymod
from itertools import combinations

def qmultinom_over_n(n,alpha):
    """Y(q) = qmultinomial(n;alpha)/[n]_q as product of cyclotomics Phi_d^{e_d}"""
    E=[]
    for d in range(2,n+1):
        e = n//d - sum(a//d for a in alpha)
        e -= 1 if n%d==0 else 0
        if e<0: return None
        E += [d]*e
    return E
def Yat(E,m):
    if E is None: return None
    if E.count(m)>0: return 0
    Pm=cyclotomic(m); acc=[1]
    for d in E:
        f=polymod(cyclotomic(d)[:],Pm); r=[0]*(len(acc)+len(f)-1)
        for i,a in enumerate(acc):
            if a:
                for j,b in enumerate(f): r[i+j]+=a*b
        acc=polymod(r,Pm)
    while acc and acc[-1]==0: acc.pop()
    return acc[0] if len(acc)<=1 else ('nonrat',tuple(acc))
def predict2(n,v,m):
    for e in range(1,m):
        if m%e: continue
        fe=gcd(pow(v,e,n)-1,n)
        if phi(fe)>2 or m%rad(fe)!=0: return False
    return True

def comps(n,parts):
    if parts==1: yield [n]; return
    for i in range(1,n-parts+2):
        for rest in comps(n-i,parts-1): yield [i]+rest

mis=0;tested=0
for n in range(4,29):
    for v in range(2,n):
        if gcd(v,n)!=1: continue
        m=order(v,n); ok=True
        for alpha in comps(n,3):
            g=0
            for a in alpha: g=gcd(g,a)
            if g!=1: continue
            E=qmultinom_over_n(n,alpha)
            if Yat(E,m)!=fixY_alpha(n,alpha,v,m): ok=False; break
        tested+=1
        if ok!=predict2(n,v,m):
            mis+=1; print("  MISMATCH(3-colour)",n,v,m,ok,predict2(n,v,m))
print(f"3-COLOUR test: multipliers tested={tested} mismatches={mis}")
