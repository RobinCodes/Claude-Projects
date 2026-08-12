# The Loop Hunt — running log

ROOT PRINCIPLE: "looping a string" = the cyclic group acting; its dual = ROOTS OF UNITY.
Everything beautiful about loops (necklaces, primes, zeta, Witt vectors, Weil conjectures,
trace formulas) is one statement: *closed orbits are the primes of a dynamical system,
and characters of the loop are roots of unity.*

## Established so far (mine, verified by exhaustive computation)
1. Realizability/Gauss-congruence lens works, but classical sequences are already
   catalogued in OEIS (A060165 family). Mined vein.
2. VERIFIED: q-analogue [n]_q/[n-k]_q * qbinom(n-k,k) cyclically sieves cyclic
   k-subsets with no two adjacent, n<=14, all k. (probably known)
3. VERIFIED + PROVED: **Mobius/twisted loop.** For ODD n, tau(w)_i = 1 - w_{i+1}
   has order 2n on {0,1}^n and
        prod_{i=1}^{n} (1 + q^i)  evaluated at zeta_{2n}^d  =  #Fix(tau^d).
   Proof: d odd -> factor (1+q^n) vanishes. d=2e -> product = 2^{gcd(e,n)} = 2^{gcd(d,n)}.
   (light, but real)
4. PROVED: **No nonzero sequence has infinite "realizability depth."**
   If b=L(a) (orbit counts) and the tower a, L(a), L^2(a), ... stays a nonneg integer
   sequence forever, then a=0.  Pf: let N=min support; zeros propagate; c^{(k+1)}_N =
   c^{(k)}_N / N, so N^k | a_N for all k.  "Primes of primes of primes..." always dies.
