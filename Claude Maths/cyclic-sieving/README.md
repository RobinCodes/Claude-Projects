# Multiplier cyclic sieving on necklaces

Work from a session on 2026-08-09/10. The headline result is a criterion deciding when a
q-binomial evaluated at a root of unity counts the necklaces a *multiplier* leaves fixed.
The answer is the crystallographic restriction: **{1, 2, 3, 4, 6}**.

**2026-08-13 follow-up:** sufficiency is now **proved**, in full generality, and one of the
paper's unproved hypotheses is **false**. The paper below has been rewritten to incorporate
this. See also `sufficiency/PROOF.md` (web version:
https://claude.ai/code/artifact/4cf026e9-fef5-4f6d-8818-21224149d990) and reproduce with
`python3 sufficiency/verify.py`.

**Paper:** `writeup/paper.pdf` (16 pp.) --- build with `pdflatex paper.tex` (run twice).
**Source:** `writeup/paper.tex`, single file, standard packages only.
**Web version:** `writeup/crystal-sieving.html`
(published at https://claude.ai/code/artifact/85b442c4-382e-4fa9-8176-b4135affd510)

Both figures --- the `(m,f)` verdict grid and the symmetry-order strip --- are drawn in the
base LaTeX `picture` environment. No TikZ/PGF, no external image files: `texlive-basic`,
`texlive-latex` and the standard graphics bundle are enough. `hyperref` and T1 `fontenc` are
loaded only if installed, so the paper also compiles without them. Compiles with no warnings
and no overfull boxes. `writeup/mkfigs.py` regenerates that figure code from the verdict
data, so the picture and the table it came from cannot drift apart.

---

## The result

Colour the positions `Z/n` with two colours, `k` of them black, and quotient by rotation:
these are **necklaces**, `C(n,k)/n` of them when `gcd(n,k)=1`. The natural q-analogue is

    Y(q) = [n choose k]_q / [n]_q          (a polynomial exactly when gcd(n,k)=1)

Rotation now acts trivially, but **multipliers** still act: for a unit `v` mod `n`, relabel
positions by `x -> v x`. These are the symmetries in the normaliser of the rotation group,
i.e. the affine group `Z/n ⋊ (Z/n)^×`.

### Criterion (exhaustively verified; sufficiency not yet proved)

Let `v` have multiplicative order `m`, and for `e | m` set

    f_e = gcd(v^e - 1, n) = #{ x in Z/n : v^e x = x }

Then

    Y(zeta_m) = #{ k-subset necklaces of Z/n fixed by v }     for every k coprime to n

holds **iff** for every proper divisor `e` of `m`:

1. `phi(f_e) <= 2`, i.e. `f_e in {1,2,3,4,6}`  — the crystallographic restriction
2. `rad(f_e) | m`

For multipliers whose orbits all have size 1 or `m`, this collapses to a condition on the
pair `(m, f)` alone — independent of `n` and of the number of orbits.

### Lemma (proved)

The engine that makes everything computable. With `A_c(y) = v y + c` and
`S = 1 + v + ... + v^(m-1)`, for any colour content `alpha` with `gcd(alpha) = 1`:

    #Fix(v) = (1/n) * sum_{c : cS = 0 mod n} N_alpha(A_c)

where `N_alpha(A_c)` counts words of content `alpha` constant on the cycles of `A_c`.

*Proof.* A necklace is `v`-fixed iff some representative is constant on the cycles of some
`A_c`. Distinct `c` never share such a word (it would be fixed by a nontrivial rotation,
impossible when `gcd(alpha)=1`), so the sum needs no inclusion–exclusion. If `cS != 0` then
`A_c^m` is a nontrivial translation of order `e > 1`; every cycle length of `A_c` is then
divisible by `e`, so every colour class is too, forcing `e | gcd(alpha) = 1`. Those terms
vanish. ∎

### Why {1,2,3,4,6}

For uniform-orbit multipliers, `n = f + mt` so `n ≡ f (mod m)`. By q-Lucas,

    Y(zeta) = C(floor(n/m), floor(k/m)) * [r choose s]_zeta / [r]_zeta,
    r = n mod m, s = k mod m

which is rational only for `s = 1`, `s = r-1` (ratio 1) or `s > r` (zero). So every *bad*
residue `s in {0, 2, 3, ..., r-2, r}` must be unreachable by `k` coprime to `n`, which forces
every prime `q <= f` other than `f-1` to divide `gcd(m,n)`. If such a `q` divided `n` but not
`f`, uniformity gives `m | q-1`, so `m <= q-1`; but `q | m` gives `q <= m`. Contradiction.
Hence every prime `q <= f`, `q != f-1`, divides `f` — which permits only 1, 2, 3, 4, 6.

---

## Verification

Exact integer + cyclotomic arithmetic throughout; no floating point.

| Form of criterion | Range | Multipliers | Counterexamples |
|---|---|---|---|
| General (all divisors `e \| m`, no uniformity hypothesis) | n = 4..96 | 2,709 | 0 |
| Uniform-orbit special case, `(m,f)` form | n = 4..115 | 2,682 | 0 |

Each multiplier is tested against *every* admissible `k`, so the instance count is in the
tens of thousands. Both sides were cross-validated against independent routes: the polynomial
side against direct polynomial division for all `n < 30`, the fixed-point side against
brute-force necklace enumeration.

### Reproduce

    python3 verify.py            # n = 4..40, about a minute
    python3 verify.py 4 96       # the published range; run in background, ~40 min

Non-zero exit status means a counterexample was found.

---

## Status: what is proved and what is not

**Updated 2026-08-13 — sufficiency is now PROVED. See `sufficiency/PROOF.md`.**

- **Proved:** the fixed-point lemma above; and now, strictly stronger and unconditional,
  `#Fix_necklaces(v) = #Fix_words(v)/f` for every content with `gcd(alpha)=1` (Theorem A).
- **Proved:** under the criterion the points of period `< m` form the cyclic subgroup `C_L`
  with `L in {1,2,3,4,6}` and `v` acting as `+-1` on it (Theorem B) — a restatement of the
  criterion as "the non-generic points form a subgroup of crystallographic order".
- **Proved:** **sufficiency of the criterion, in full generality** — no uniformity hypothesis
  (Theorem C). The proof reduces to: `Q(z) = prod_{|O|<m}(1+z^|O|)` is palindromic of degree
  `L` with `[z^1]Q = f`, and the achievable residues are exactly `{1, L-1}`, which is
  precisely the condition `phi(L) <= 2`.
- **Disproved:** Hypothesis 3.6 of the paper. The smallest counterexample is `m=9, r=8, s=4`,
  where the ratio equals `-1`; there are 67 with `m <= 60`, all at `r = m-1`, all `+-1`.
  The first multiplier at which the paper's necessity proof actually depends on it is
  `n = 99, v = 19` — just past the verified range `n <= 96`. Necessity still holds there,
  but because `Y_{99,13}(zeta_10) = -9` while the true count is `84`: the evaluation is a
  rational, indeed negative, integer that is simply wrong.
- **Verified only:** necessity in the general multi-level case.

### Open

- **Three or more colours is strictly stronger.** At `n = 28` the two-colour sieving holds for
  multipliers where the three-colour version fails (`multicolor.py`). Unexplained — but
  Theorem A now applies verbatim to every content, so the discrepancy is entirely on the
  `q`-analogue side, not the counting side.
- **Necessity** in the non-uniform multi-level case (sufficiency is done).
- **Residual irrationality claim:** for `2 <= r <= m-2`, is `qbinom(r,s)|_zeta / [r]_zeta`
  always irrational at the bad residues? Verified `m <= 60`. A weaker statement suffices for
  necessity: *if the ratio is rational then it has modulus 1* (proved for `r = m-1`).
- **BWT-fixed necklaces** (`bwt.py`, `bwtfix.py`): among 4,116 necklaces of length 16 exactly
  11 are their own Burrows–Wheeler transform. Counts for n = 2..18 are
  `3,4,4,4,7,6,5,9,5,6,9,8,8,11,11,7,15` — not in the OEIS. Irregular and uncharacterised.

---

## Prior work / novelty ledger

- **New:** the criterion; fixed-point counts 3, 4, 6 are cases the literature does not cover.
  They sit exactly in the regime E. N. Stucky flags in arXiv:1812.04578 §3 as one where
  *"we do not have a general conjecture."*
- **New:** the exact fixed-point lemma.
- **Rediscovered:** sieving for prime `n` — derived and proved here by q-Lucas before being
  found as Stucky's Corollary 3.3. His.
- **Rediscovered:** `([n]_q/[n-mk]_q) * [n-mk choose k]_q` sieves k-subsets of a cycle with all
  cyclic gaps `>= m+1` (`sep.py`, verified to n=16). Known.
- **Rediscovered:** central binomials, Apéry, Franel, Domb numbers all count periodic points of
  some dynamical system; their orbit counts are the OEIS A060165 family (`zoo.py`).
- **Side result, proved** (`data/exploration-log.md`): iterated orbit-counting terminates. If
  `L` maps a periodic-point sequence to its orbit counts, then for any nonzero `a >= 0` some
  `L^j(a)` leaves the nonnegative integers — bounded by `log_N(a_N)` at the first nonzero term.

References:
- E. N. Stucky, *Parity-Unimodality and a Cyclic Sieving Phenomenon for Necklaces*, arXiv:1812.04578
- V. Reiner, D. Stanton, D. White, *The cyclic sieving phenomenon*, JCTA 108 (2004)

---

## File guide

Everything is flat Python 3, standard library only (`urllib` for the OEIS oracle). Keep the
files in one directory — they import each other by module name.

### The 2026-08-13 proof

| File | Role |
|---|---|
| `sufficiency/PROOF.md` | Theorems A, B, C and Proposition D. The write-up. |
| `sufficiency/verify.py` | **Entry point.** Re-checks every claim; `--full` for published ranges. |
| `sufficiency/core.py` | Exact cyclotomics, `Y(zeta_m)`, orbit and affine-cycle machinery. |
| `sufficiency/criterion.html` | Web version of the write-up. |

### Core (the result)

| File | Role |
|---|---|
| `verify.py` | **Entry point.** Verifies the criterion over a range of `n`. |
| `crit.py` | Fast exact evaluation of `Y(q)` at a root of unity, via `Y(q) = prod Phi_d(q)` over `d` with `d ∤ n` and `(k mod d) > (n mod d)`. |
| `lib2.py` | Fixed-point counting via the lemma; `phi`, `rad`, the predicate. Supports arbitrary colour contents. |
| `csp.py` | Cyclotomic polynomials, polynomial mod, exact root-of-unity evaluation, generic CSP tester. |
| `canon.py` | q-integers, q-binomials, exact polynomial division, canonical CSP distributions. |

### Verification runs

| File | Role |
|---|---|
| `general2.py` | The general criterion, no uniformity hypothesis. Produced `data/verify-general-n47-96.txt`. |
| `final.py` | The uniform `(m,f)` criterion. Produced `data/verify-uniform-n61-115.txt`. |
| `uniform.py`, `uniform2.py` | Build the `(m,f)` verdict table. Produced `data/uniform-mf-table-n4-84.txt`. |
| `general.py` | The *simpler* criterion using only `f_1`. Fails — kept because it documents why the multi-level form is needed. |
| `nonuniform.py` | Shows uniformity is not necessary: 38 of 280 non-uniform multipliers still sieve. |
| `multicolor.py` | The three-colour test, where the criterion becomes strictly stronger. |
| `profile.py`, `atomic.py`, `fast.py`, `mcsp.py` | Earlier pass/fail sweeps and cross-tabulations. `mcsp.py` is the slow brute-force used to validate the fast path. |

### Exploration (the search that got there)

| File | Role |
|---|---|
| `oeis.py` | OEIS novelty oracle with an on-disk cache (`oeis_cache/`). |
| `zoo.py`, `lib.py` | Realizability / Gauss-congruence sweep over classical sequences. |
| `mult.py` | First hit: the multiplier CSP for prime `n`. |
| `mobius.py` | Twisted / Möbius-loop sieving at 2n-th roots of unity. |
| `sep.py` | m-separated cyclic subsets. |
| `bwt.py`, `bwtfix.py` | Burrows–Wheeler transform as a dynamical system on necklaces. |
| `batch1.py` | Circular autocorrelation / homometry statistics ("can one hear a necklace?"). |
| `wheel.py`, `fold.py`, `engine.py`, `validate.py` | Spanning trees of wheels, folded cyclotomics, the broad CSP engine, RSW validation. |
| `verify1.py` | Superseded slow verifier; kept for provenance. |

### Data

`data/exploration-log.md` is the running log with the proved side results.
The `.txt` files are raw outputs of the verification runs.
