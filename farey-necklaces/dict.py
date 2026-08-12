"""The necklace <-> quadratic-form dictionary, checked numerically.

achiral   = necklace fixed by reversal            (conjecture: ambiguous class)
selfcomp  = necklace fixed by L<->R complement    (conjecture: trace = s^2+2 tower)
"""
import sys
from necktrace import necklaces_by_trace

def runs_to_word(runs):
    return "".join(("L" if i % 2 == 0 else "R") * a for i, a in enumerate(runs))

def rotations(w):
    return {w[i:] + w[:i] for i in range(len(w))}

def canon(w):
    return min(rotations(w))

def analyse(T):
    nt = necklaces_by_trace(T)
    rows = []
    for t in range(3, T + 1):
        ws = [runs_to_word(r) for r in nt.get(t, [])]
        ach = sum(1 for w in ws if w[::-1] in rotations(w))
        sc = sum(1 for w in ws if w.translate(str.maketrans("LR", "RL")) in rotations(w))
        both = sum(1 for w in ws
                   if w[::-1] in rotations(w)
                   and w.translate(str.maketrans("LR", "RL")) in rotations(w))
        rows.append((t, len(ws), ach, sc, both))
    return rows

def omega(n):
    """number of distinct prime factors"""
    c, d = 0, 2
    while d * d <= n:
        if n % d == 0:
            c += 1
            while n % d == 0:
                n //= d
        d += 1
    return c + (1 if n > 1 else 0)

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    print(f"{'t':>4} {'h':>4} {'achiral':>8} {'selfcomp':>9} {'both':>5} {'w(t^2-4)':>9} {'2^(w-1)':>8}  {'s^2+2?':>7}")
    for (t, h, ach, sc, both) in analyse(T):
        D = t * t - 4
        w = omega(D)
        sq = ""
        for s in range(1, t):
            if s * s + 2 == t:
                sq = f"s={s}"
        print(f"{t:>4} {h:>4} {ach:>8} {sc:>9} {both:>5} {w:>9} {2**(w-1):>8}  {sq:>7}")
