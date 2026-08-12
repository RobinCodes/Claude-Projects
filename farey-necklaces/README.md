# Loop and Swap — the Farey spin chain and Zagier-reduced forms

Session of 2026-08-11. Write-up:
https://claude.ai/code/artifact/20acb34f-dbdf-4f8d-be2e-ea76681e610e

## The result

`L = [[1,1],[0,1]]`, `R = [[1,0],[1,1]]` generate a **free** monoid = all non-negative
integer matrices of determinant 1. Let

    Phi(t) = #{ words in L,R whose matrix has trace t }        (OEIS A264598)
    Z(D)   = #{ Zagier-reduced forms (A,B,C): A,C>=1, B>A+C, B^2-4AC=D }   (A257007 at D=t^2-4)

**Sloane conjectured in OEIS (Nov 2015) that `Phi(t) = 2*Z(t^2-4)`, and the entries still
carry "It would be nice to have a proof".** It is true. Two proofs:

**1. The swap.** Trace ignores rotation, so the objects are necklaces = conjugacy classes
= ideal classes of discriminant `t^2-4`. Each class contributes its primitive period
`p = #L + #R` matrices, so `Phi(t) = sum_C (#L + #R)`. Zagier's cycle of a class has length
equal to the minus-CF period, which by the plus/minus dictionary
`[a1,a2,a3,...]_+ = [a1+1, 2^(a2-1), a3+2, ...]_-` equals `#R`. So `Z = sum_C #R`.
Conjugation by `J = [[0,1],[1,0]]` swaps L and R; it has det -1 but preserves trace and
determinant, so it permutes the classes of trace t and exchanges `#L` with `#R`, giving
`sum #L = sum #R`. Hence `Phi = 2Z`.  (Per-class refinement, strictly stronger:
**Zagier cycle length of a class = #R of its word**.)

**2. Two lines from published work.** Technau (arXiv:2304.08143, Thm 2) proved
`Phi(t) = Upsilon(t^2-4)` with `Upsilon(D) = sum_{k^2<D, k^2=D mod 4} d((D-k^2)/4)`.
Zagier-reduced forms = pairs `(h,k)` with `h | (D-k^2)/4` and `2h+k > sqrt(D)`. The map
`(h,k) -> ((D-k^2)/(4h), -k)` is an involution on all such pairs, fixed-point-free
(`t^2-4` is never a square), and it flips the inequality. So exactly half are reduced:
`Upsilon = 2Z`. Technau never mentions Zagier reduction and OEIS never mentions Technau,
which is why this sat unnoticed.

## Also new here

- **Cut-sum.** `S(w) = sum of the matrices of all n cuts of a necklace` — a well-defined
  *matrix* invariant, finer than the trace. Exact identity
  `det S = n + sum over cut-pairs of tr[u,v]` (polarisation of det + Fricke).
  `Delta = (tr S)^2 - 4 det S` = squared Minkowski norm of the sum of the n pole vectors:
  a centre of mass for a closed geodesic over its Farey crossings.
  Delta separates all classes mod <rot,reversal,swap> for t<=34 and n<=14; first
  collisions at n=15 (2 pairs) and n=16 (3 pairs), where all of S agrees. Unexplained.
- **Total matrix.** `sum_{tr X = t} X = (t/2) Phi(t) I + (sum_a sigma(a(t-a)-1)) J`.
  The `d` half is Technau's theorem; the `sigma` half looks untouched — the natural next
  asymptotics question for the spin chain (Hooley / Bykovskii-Ustinov did `d`).
- Observation: a necklace is self-complementary iff its trace is `s^2+2` (it is the square
  of a determinant -1 matrix).

## Files

Flat Python 3, stdlib only, exact integer arithmetic, no floats in any decision.

    final.py     END-TO-END CHECK.  `python3 final.py 90` -> "ALL CLAIMS VERIFIED"
    zagier.py    Phi vs Zagier count vs sum #L / sum #R      (checked t=3..150)
    zcycle.py    per-class Zagier cycles vs #R multisets     (checked t=3..24)
    necktrace.py necklace enumeration by trace; Gauss-reduced form cycles
    mat.py       fast enumeration of SL_2(Z>=0) by trace, word statistics
    moment.py    cut-sum identities (`python3 moment.py check`)
    collide.py   the exceptional cut-sum coincidences at n=15,16
    cutsum.py sep.py config.py rade.py stats.py spec.py dict.py   exploration

## Open

Explain the n=15,16 cut-sum coincidences; asymptotics of `sum_a sigma(a(t-a)-1)`;
whether Delta-collisions form an infinite family. OEIS A264598/A257007/A264597 should be
updated to remove the conjecture flags.

## Paper

`paper/paper.tex` — formal write-up, 10 pp., compiles with `pdflatex paper.tex` (run twice);
amsart + amsmath/amssymb/booktabs/geometry/hyperref only, no TikZ. `paper/paper.pdf` is the
built copy, `paper/paper.html` the web version
(https://claude.ai/code/artifact/a43d56f6-7978-44a1-9a0b-e440a7e9905f).

Important correction found while writing it up: Technau's *intermediate* equation (3.2)
already reads `Phi(N) = 2 Upsilon^<(N^2-4)` — the factor 2 was proved in 2023. The only
missing step was recognising `Upsilon^<` as the count of Zagier-reduced forms, which is the
change of variables `(A,B,C) -> (A, B-2A)` of Lemma 3.1. The involution proof above is a
re-derivation of Technau's own Lemma 4 in different coordinates. The necklace/swap proof is
independent of all of it and gives the finer per-class statement.
