from math import gcd
from csp import cyclotomic, polymod
from canon import qbinom, qint, polydivexact

def Yfactors(n,k):
    return [d for d in range(2,n+1) if n%d and (k%d)>(n%d)]

def Yat_fast(n,k,m):
    """Y(zeta_m) exactly, as element of Z[q]/Phi_m (list). Rational -> int."""
    E=Yfactors(n,k)
    if m in E: return 0
    Pm=cyclotomic(m); acc=[1]
    for d in E:
        f=polymod(cyclotomic(d)[:],Pm)
        r=[0]*(len(acc)+len(f)-1)
        for i,a in enumerate(acc):
            if a:
                for j,b in enumerate(f): r[i+j]+=a*b
        acc=polymod(r,Pm)
    while acc and acc[-1]==0: acc.pop()
    if len(acc)<=1: return acc[0] if acc else 0
    return ('nonrat',tuple(acc))

