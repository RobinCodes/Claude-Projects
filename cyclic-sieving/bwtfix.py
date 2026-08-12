from bwt import bwt, canon, necklaces
for n in range(2,19):
    N=necklaces(n)
    fx=[w for w in N if canon(bwt(w))==w]
    print(f"n={n:2d} ({len(fx):2d}): "+"  ".join("".join(map(str,w)) for w in sorted(fx,key=lambda w:(sum(w),w))))
