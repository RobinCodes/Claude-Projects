from lib import *
from math import comb, factorial, isqrt

M = 26  # test n=1..M

def seq(f): return [0]+[f(n) for n in range(1,M+1)]

Z = {}
Z['C(2n,n)']        = seq(lambda n: comb(2*n,n))
Z['Catalan']        = seq(lambda n: comb(2*n,n)//(n+1))
Z['C(3n,n)']        = seq(lambda n: comb(3*n,n))
Z['C(4n,2n)']       = seq(lambda n: comb(4*n,2*n))
Z['C(2n,n)^2']      = seq(lambda n: comb(2*n,n)**2)
Z['Franel']         = seq(lambda n: sum(comb(n,k)**3 for k in range(n+1)))
Z['sumC^2=C(2n,n)'] = seq(lambda n: sum(comb(n,k)**2 for k in range(n+1)))
Z['sumC^4']         = seq(lambda n: sum(comb(n,k)**4 for k in range(n+1)))
Z['Apery A']        = seq(lambda n: sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1)))
Z['Apery B']        = seq(lambda n: sum(comb(n,k)**2*comb(n+k,k) for k in range(n+1)))
Z['Motzkin']        = seq(lambda n: sum(comb(n,2*k)*comb(2*k,k)//(k+1) for k in range(n//2+1)))
Z['Delannoy']       = seq(lambda n: sum(comb(n,k)*comb(n+k,k) for k in range(n+1)))
Z['lgSchroeder']    = seq(lambda n: sum(comb(n,k)*comb(n+k,k)//(k+1)*1 for k in range(n+1)) )
Z['Bell']           = None
Z['partitions']     = None
Z['n!']             = seq(lambda n: factorial(n))
Z['derangements']   = seq(lambda n: round(factorial(n)/2.718281828459045) if n>0 else 1)
Z['2^n']            = seq(lambda n: 2**n)
Z['2^n-1']          = seq(lambda n: 2**n-1)
Z['2^n+1']          = seq(lambda n: 2**n+1)
Z['3^n-2^n']        = seq(lambda n: 3**n-2**n)
Z['n^n']            = seq(lambda n: n**n)
Z['Fib']            = None
Z['Lucas']          = None
Z['Perrin']         = None
Z['sigma(n)']       = seq(lambda n: sum(divisors(n)))
Z['involutions']    = None
Z['Domb']           = seq(lambda n: sum(comb(n,k)**2*comb(2*k,k)*comb(2*(n-k),n-k) for k in range(n+1)))
Z['C(2n,n)*3^?']    = seq(lambda n: sum(comb(n,k)*comb(2*k,k)*comb(2*(n-k),n-k) for k in range(n+1)))
Z['centralTrinom']  = seq(lambda n: sum(comb(n,k)*comb(n-k,k) for k in range(n//2+1)))
Z['Fine']           = None
Z['Bessel/invol']   = None

# recursive ones
def bell(M):
    B=[1]
    row=[1]
    for i in range(M):
        new=[row[-1]]
        for x in row: new.append(new[-1]+x)
        row=new; B.append(row[0])
    return B
b=bell(M+1); Z['Bell']=[0]+[b[n] for n in range(1,M+1)]

def parts(M):
    p=[0]*(M+1); p[0]=1
    for k in range(1,M+1):
        for i in range(k,M+1): p[i]+=p[i-k]
    return p
p=parts(M); Z['partitions']=[0]+[p[n] for n in range(1,M+1)]

F=[0,1]
for i in range(2,M+2): F.append(F[-1]+F[-2])
Z['Fib']=[0]+[F[n] for n in range(1,M+1)]
L=[2,1]
for i in range(2,M+2): L.append(L[-1]+L[-2])
Z['Lucas']=[0]+[L[n] for n in range(1,M+1)]
P=[3,0,2]
for i in range(3,M+2): P.append(P[-2]+P[-3])
Z['Perrin']=[0]+[P[n] for n in range(1,M+1)]

I=[1,1]
for n in range(2,M+2): I.append(I[-1]+(n-1)*I[-2])
Z['involutions']=[0]+[I[n] for n in range(1,M+1)]

Fine=[1,0,1,2,6,18,57,186,622,2120,7338,25724,91144,325878,1174281,4260282,15548694,57048048,210295326,778483932,2891461209,10778209834,40312874461,151245191722,569074240533,2146336125648,8118259581354]
Z['Fine']=[0]+Fine[1:M+1]

rows=[]
for name,a in Z.items():
    if a is None: continue
    intg,nn,bs,bad = classify(a,name)
    rows.append((name,intg,nn,bad,bs[:12]))
for r in sorted(rows, key=lambda r:(not r[1], not r[2])):
    flag = "REALIZABLE " if (r[1] and r[2]) else ("gauss-only " if r[1] else "           ")
    print(f"{flag}{r[0]:16s} int={r[1]!s:5s} nonneg={r[2]!s:5s} bad={str(r[3]):12s} b={r[4]}")
