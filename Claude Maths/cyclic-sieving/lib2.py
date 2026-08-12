from math import gcd
def order(u,n):
    o=1;x=u%n
    while x!=1: x=x*u%n; o+=1
    return o
def cyc(v,c,n):
    seen=bytearray(n); s=[]
    for y in range(n):
        if seen[y]: continue
        l=0;z=y
        while not seen[z]: seen[z]=1; z=(v*z+c)%n; l+=1
        s.append(l)
    return s
def Ncontent(sizes, alpha):
    """# of ways to color cycles so colour i is used alpha_i times"""
    from functools import lru_cache
    dp={tuple([0]*len(alpha)):1}
    for s in sizes:
        nd={}
        for st,c in dp.items():
            for i in range(len(alpha)):
                if st[i]+s<=alpha[i]:
                    t=list(st); t[i]+=s; t=tuple(t)
                    nd[t]=nd.get(t,0)+c
        dp=nd
    return dp.get(tuple(alpha),0)
def fixY_alpha(n,alpha,v,m):
    S=0;p=1
    for _ in range(m): S=(S+p)%n; p=p*v%n
    g=gcd(S,n) if S else n
    step=n//g
    tot=sum(Ncontent(cyc(v,c,n),alpha) for c in range(0,n,step))
    assert tot%n==0
    return tot//n
def rad(x):
    r=1;d=2;y=x
    while d*d<=y:
        if y%d==0:
            r*=d
            while y%d==0: y//=d
        d+=1
    if y>1: r*=y
    return r
def phi(x):
    r=x;d=2;y=x
    while d*d<=y:
        if y%d==0:
            while y%d==0: y//=d
            r-=r//d
        d+=1
    if y>1: r-=r//y
    return r
def predict(m,f): return phi(f)<=2 and m%rad(f)==0
