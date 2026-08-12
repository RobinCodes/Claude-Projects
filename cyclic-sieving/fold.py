from canon import polymul, polydivexact

def cyclo(n, cache={}):
    if n in cache: return cache[n]
    num=[-1]+[0]*(n-1)+[1]
    for d in range(1,n):
        if n%d==0: num=polydivexact(num, cyclo(d))
    cache[n]=num; return num

def fold(p,m):
    v=[0]*m
    for i,c in enumerate(p): v[i%m]+=c
    return v

def s(v): return "["+",".join(f"{x:+d}" for x in v)+"]"

print("Fold of Phi_n modulo x^m - 1   (m | n)")
for n in [4,8,9,12,16,18,20,24,27,15,21,35,45,105,30,60]:
    ds=[d for d in range(2,n+1) if n%d==0]
    print(f"\n n={n}  Phi_n deg={len(cyclo(n))-1}")
    for m in ds:
        f=fold(cyclo(n),m)
        print(f"   m={m:3d}: {s(f)}")
