"""End-to-end verification of every claim.

(1) Phi(t) := #{X in free monoid <L,R> : tr X = t} = sum_{a=1}^{t-1} d(a(t-a)-1)
              [ = Technau 2023, Thm 2:  Phi(N) = Upsilon(N^2-4) ]
(2) Upsilon(D) = 2 * Z(D), Z = #Zagier-reduced forms,  via the fixed-point-free
    involution (h,k) -> (m_k/h, -k) on divisor pairs.        [two lines, elementary]
(3) Z(D) = sum over classes of #R  (per class: Zagier cycle length = #R)
(4) L<->R swap is a bijection on classes of trace t swapping #L and #R
        => sum #L = sum #R  => Phi(t) = 2 Z(t^2-4)         [A264598 = 2*A257007]
(5) TOTAL MATRIX:  sum_{tr X = t} X  =  (t/2) Phi(t) I  +  (sum_a sigma(a(t-a)-1)) J
"""
import sys
from math import isqrt
from necktrace import necklaces_by_trace
from mat import matrices
from zagier import divisors, zagier_count, prim_root

def d(n):
    return len(divisors(n))

def sigma(n):
    return sum(divisors(n))

def upsilon(D):
    """#{(lam,mu,m) in N^2 x Z : 4 lam mu + m^2 = D}"""
    s = isqrt(D)
    tot = 0
    for m in range(-s, s + 1):
        if m * m >= D or (D - m * m) % 4:
            continue
        q = (D - m * m) // 4
        if q >= 1:
            tot += d(q)
    return tot

def involution_check(D):
    """the map (h,k) -> (m_k/h, -k) is a fixed-point-free involution on
       P = {(h,k): |k|<sqrt(D), k^2=D mod 4, h | (D-k^2)/4}
       and it swaps the Zagier condition 2h+k > sqrt(D) with its negation."""
    s = isqrt(D)
    P = []
    for k in range(-s, s + 1):
        if k * k >= D or (D - k * k) % 4:
            continue
        m = (D - k * k) // 4
        if m < 1:
            continue
        for h in divisors(m):
            P.append((h, k))
    def cond(h, k):
        return (2 * h + k) > 0 and (2 * h + k) ** 2 > D
    for (h, k) in P:
        m = (D - k * k) // 4
        h2, k2 = m // h, -k
        assert (h2, k2) in set(P)
        assert (h, k) != (h2, k2), "fixed point!"
        assert cond(h, k) != cond(h2, k2), "involution does not swap the condition"
    return len(P), sum(1 for (h, k) in P if cond(h, k))

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    nt = necklaces_by_trace(T)
    allok = True
    print(f"{'t':>4} {'Phi':>7} {'sum d':>7} {'Upsilon':>8} {'Z':>6} {'sumR':>6} {'sumL':>6} "
          f"{'S11':>9} {'(t/2)Phi':>9} {'S12':>9} {'sum sig':>9}  ok")
    for t in range(3, T + 1):
        D = t * t - 4
        Ms = list(matrices(t))
        Phi = len(Ms)
        sd = sum(d(a * (t - a) - 1) for a in range(1, t))
        U = upsilon(D)
        Z = zagier_count(D)
        roots = [prim_root(r) for r in nt.get(t, [])]
        sR = sum(u.count('R') for u in roots)
        sL = sum(u.count('L') for u in roots)
        S11 = sum(m[0] for m in Ms)
        S12 = sum(m[1] for m in Ms)
        S21 = sum(m[2] for m in Ms)
        S22 = sum(m[3] for m in Ms)
        ssig = sum(sigma(a * (t - a) - 1) for a in range(1, t))
        nP, nC = involution_check(D)
        ok = (Phi == sd == U == 2 * Z and Z == sR == sL and nP == U and nC == Z
              and S11 == S22 == t * Phi // 2 and S12 == S21 == ssig)
        allok &= ok
        if t <= 20 or not ok:
            print(f"{t:>4} {Phi:>7} {sd:>7} {U:>8} {Z:>6} {sR:>6} {sL:>6} "
                  f"{S11:>9} {t*Phi//2:>9} {S12:>9} {ssig:>9}  {ok}")
    print(f"\nchecked t = 3..{T}:", "ALL CLAIMS VERIFIED" if allok else "*** FAILURE ***")
