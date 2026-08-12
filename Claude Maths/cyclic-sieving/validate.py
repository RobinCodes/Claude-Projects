from csp import *
from itertools import combinations
# Validate on the classic RSW theorem: k-subsets of Z_n under rotation, X = q-binomial with stat = sum(S) - (k choose 2)
n, k = 8, 3
objs = list(combinations(range(n), k))
def act(s): return tuple(sorted((x+1)%n for x in s))
def stat(s): return sum(s) - k*(k-1)//2
ok, det, co = test_csp(objs, act, stat, n)
print("RSW q-binomial CSP (n=8,k=3):", ok)
for d,v,c,g in det: print("  d=%d  X(z^d)=%s  Fix=%d  %s"%(d,v,c,"OK" if g else "FAIL"))
