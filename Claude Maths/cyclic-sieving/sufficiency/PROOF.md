# Sufficiency of the multiplier sieving criterion — and a false hypothesis

Session of 2026-08-12/13. Web version:
https://claude.ai/code/artifact/4cf026e9-fef5-4f6d-8818-21224149d990

**The formal write-up is now `../writeup/paper.tex` (16 pp.), which has been rewritten to
incorporate everything below.** This file is the informal companion: same mathematics, less
ceremony. Numbering in the paper: Theorem A = Thm 3.2 (collapse), Theorem B = Thm 7.2
(structure), Theorem C = Thm 6.3 (sufficiency), Proposition D = Prop 5.1 / Cor 5.2 /
Prop 5.3 and Remarks 5.4, 9.5.

This closes the main open problem of `writeup/paper.pdf`
("Verified only: sufficiency, and the general multi-level form") and corrects one of its
unproved hypotheses, which turns out to be **false**.

Three results:

| | Statement | Status before | Status now |
|---|---|---|---|
| **A** | `#Fix_necklaces(v) = #Fix_words(v) / f`, for every content with `gcd(alpha)=1` | Lemma 3.1 gave a sum over `c in K` | **Proved**, unconditionally, in four lines; strictly stronger |
| **C** | The criterion is **sufficient** — for *all* multipliers, not just uniform ones | Conjectured, verified to `n=96` | **Proved** |
| **D** | Hypothesis 3.6 (an irrationality claim) | Unproved, "verified for all `m <= 96`" | **False.** 67 counterexamples with `m <= 60`; corrected statement proved where it matters |

Notation is that of the paper: `n >= 3`, `v` a unit mod `n` of order `m > 1`,
`f_e = gcd(v^e - 1, n) = |Fix(v^e)|`, `f = f_1`, `rad` the radical,
`Y_{n,k}(q) = qbinom(n,k)/[n]_q`, and

> **Criterion.** For every proper divisor `e | m`: `f_e in {1,2,3,4,6}` and `rad(f_e) | m`.

---

## 1. Theorem A — the collapse

Let `A_c(y) = v y + c` be the affine maps with linear part `v`, and
`S_j = 1 + v + ... + v^(j-1)`.

> **Theorem A.** Let `alpha` be any content with `gcd(alpha) = 1` and let `v` be any
> multiplier. Then every `v`-fixed `alpha`-necklace has **exactly `f` representatives fixed
> by `v` as words**, and consequently
>
>     #Fix_necklaces(v)  =  (1/f) * #Fix_words(v),
>
> where `#Fix_words(v)` counts words of content `alpha` constant on the `<v>`-orbits of `Z/n`.
> In the two-colour case, writing `n_d` for the number of `<v>`-orbits of size exactly `d`,
>
>     #Fix(v)  =  (1/f) * [z^k] prod_{d | m} (1 + z^d)^{n_d}.

*Proof.* By Steps 1–2 of the paper's Lemma 3.1, `[w]` is `v`-fixed iff `w` is constant on the
cycles of exactly one `A_c`, so `#Fix = (1/n) sum_c N_alpha(A_c)`.

Let `l` be any cycle length of `A_c`. Then `A_c^l` has a fixed point, i.e.
`(v^l - 1) y = -c S_l` is solvable, so `f_l | c S_l`. Now `f | f_l` (because `v-1 | v^l-1`),
and `v = 1 (mod f)` gives `S_l = l (mod f)`. Hence

    f | c*l   for every cycle length l of A_c.

So `e_c := f / gcd(f, c)` divides **every** cycle length of `A_c`. If `f` does not divide `c`
then `e_c > 1`; take a prime `p | e_c`. Then `p | f | n`, and `p` divides every cycle length,
so every colour class of a word constant on the cycles of `A_c` — being a union of cycles —
has size divisible by `p`. Thus `p | gcd(alpha) = 1`, a contradiction. Therefore
`N_alpha(A_c) = 0` unless `f | c`.

The set `{c : f | c}` is exactly `(v-1)(Z/n)`, of size `n/f`, and for `c = (v-1)d` the
translation `y -> y + d` conjugates `A_c` to `A_0 = ` multiplication by `v`. So all `n/f`
surviving terms are equal to `N_alpha(v)`, and
`#Fix = (1/n)(n/f) N_alpha(v) = N_alpha(v)/f`. ∎

This is strictly stronger than Lemma 3.1: it identifies *which* `c` contribute (only
`c in (v-1)Z/n`), shows they all contribute equally, and needs no hypothesis whatever — in
particular no uniformity and no criterion. It reduces the whole subject to cycle-index
bookkeeping for a single permutation.

*Verified:* 29,780 instances (`n <= 60`, two colours, all `v`, all admissible `k`, criterion
true **and** false) and 34,186 instances with three colours (`n <= 28`, all contents with
`gcd(alpha)=1`). Zero failures.

---

## 2. Theorem B — what the criterion actually says

Write `P = { y in Z/n : the <v>-orbit of y has size < m }`. Since a period `< m` divides some
maximal proper divisor `m/p`,

    P = union over primes p | m of Fix(v^(m/p)),        g_p := f_(m/p).

**Lemma 2.1.** The criterion holds iff `g_p in {1,2,3,4,6}` and `rad(g_p) | m` for every prime
`p | m`.
*Proof.* Every proper `e | m` divides some `m/p`, so `f_e | g_p`; and every divisor of an
element of `{1,2,3,4,6}` is again in `{1,2,3,4,6}`, with `rad(f_e) | rad(g_p) | m`. ∎

**Lemma 2.2 (no clash).** Assume the criterion. Then the `g_p` are totally ordered by
divisibility.
*Proof.* The unordered pairs from `{1,2,3,4,6}` that are incomparable are `{2,3}`, `{3,4}`,
`{4,6}`. Two standing facts: (i) `g_p | n`; (ii) if `2 | n` then `v` is odd, so `2 | v-1 |
v^e-1`, hence `2 | f_e` for **every** `e`.
- `{2,3}`: `g_p = 2` forces `2 | n`, so by (ii) every `f_e` is even — contradicting `g_q = 3`.
- `{3,4}`: `g_q = 4` forces `2 | n`, same contradiction with `g_p = 3`.
- `{4,6}`: `4 | n` and `6 | n`, so `12 | n`. `g_p = 4` exactly means `3` does not divide
  `v^(m/p) - 1`, so `ord_3(v) = 2` and `2` does not divide `m/p`. `g_q = 6` exactly means `4`
  does not divide `v^(m/q) - 1`, so `ord_4(v) = 2` and `2` does not divide `m/q`. But
  `3 | n` gives `ord_3(v) = 2 | m`, so `m` is even; then `2 ∤ m/p` and `2 ∤ m/q` force
  `p = q = 2`, contradicting `p != q`. ∎

> **Theorem B.** Assume the criterion, and put `L = max_p g_p`. Then `P` is precisely the
> cyclic subgroup of `Z/n` of order `L`, with
>
>     L in {1,2,3,4,6},   rad(L) | m,   and v acts on P = C_L as +1 or -1.

*Proof.* By Lemma 2.2 the subgroups `Fix(v^(m/p)) = C_{g_p}` are nested, so their union is
`C_L`. `v` restricts to an automorphism of `C_L`, and `Aut(C_L) = (Z/L)^x` has order
`phi(L) <= 2`, so `v` acts as `±1`. ∎

Equivalently: **the criterion says exactly that the non-generic points form a subgroup of
crystallographic order.** The eight possible shapes `(L, eps)` are
`(1,+), (2,+), (3,±), (4,±), (6,±)`, and all eight occur.

*Verified:* 19,820 criterion-satisfying multipliers with `n <= 400`; `P` is the subgroup of
order `L in {1,2,3,4,6}` and `v` acts as `±1` in every single case, all eight shapes present.

---

## 3. Theorem C — sufficiency

> **Theorem C.** If the criterion holds then
> `Y_{n,k}(zeta_m) = #{k-subset necklaces of Z/n fixed by v}` for every `k` coprime to `n` —
> with no uniformity hypothesis.

By Theorem A it suffices to prove the purely combinatorial identity

    [z^k] N(z)  =  f * Y_{n,k}(zeta_m),      N(z) = prod_{d|m} (1 + z^d)^{n_d}.       (*)

### 3.1 Splitting off the small part

Factor `N(z) = Q(z) * (1 + z^m)^{n_m}` where `Q(z) = prod_{orbits O with |O| < m} (1 + z^|O|)`.
By Theorem B the orbits inside `P` partition `C_L`, so

    deg Q = |P| = L,   and Q is palindromic (a product of palindromic factors),
    [z^1] Q = #{orbits of size 1} = f,   hence  [z^(L-1)] Q = [z^1] Q = f.            (**)

Write `n = a m + r`, `k = b m + s` with `0 <= r,s < m`. Since `n = L + m n_m` we get
`r = L mod m`, and `u := (L - r)/m >= 0`. Because `rad(L) | m`, the pairs with `m <= L` are
exactly `(L,m) in {(2,2), (3,3), (4,2), (4,4), (6,6)}`, and all of them have `r = 0`. So

    m > L  <=>  m does not divide n,  and then r = L, u = 0;
    m <= L <=>  m divides n,          and then r = 0, u = L/m >= 1.

### 3.2 Which residues `s` occur

`s = k mod m` is realisable by some `k` coprime to `n` only if `p` does not divide `s` for
every prime `p | gcd(n,m)`. As `L | n` and `rad(L) | m`, every prime of `L` lies in
`gcd(n,m)`, so

    every achievable s satisfies gcd(s, rad(L)) = 1.                                   (***)

For `L in {2,3,4,6}` the integers `s` with `0 <= s <= L` and `gcd(s, rad(L)) = 1` are exactly
`s in {1, L-1}`: `s = 0` and `s = L` share every prime of `L`, and each `s` with
`2 <= s <= L-2` shares one too. **That last sentence is the crystallographic restriction**: it
holds precisely when `phi(L) <= 2`, i.e. `L in {1,2,3,4,6}`. (For `L = 1`, `s in {0,1}`.)

### 3.3 The case `m` does not divide `n`  (`r = L`, `u = 0`)

By `q`-Lucas, `Y_{n,k}(zeta_m) = C(a,b) * rho_s` with `rho_s = qbinom(L,s)|_zeta / [L]_zeta`,
and `n_m = a`. Extracting coefficients,
`[z^k] N = sum_j C(a, b-j) * [z^(jm+s)] Q`. Since `deg Q = L < m`, only `j = 0` survives, so
`[z^k] N = C(a,b) * [z^s] Q`. Hence (*) reduces to `[z^s] Q = f * rho_s` for achievable `s`:

- `s > L`: `[z^s]Q = 0` and `qbinom(L,s) = 0`, so both sides vanish.
- `s in {1, L-1}` (all that remains, by 3.2): `qbinom(L,1) = qbinom(L,L-1) = [L]_q`, so
  `rho_s = 1`; and `[z^s] Q = f` by (**). Both sides equal `f`.
- `L = 1`: `s in {0,1}`, `rho_s = 1`, `Q = 1+z`, `f = 1`. Both sides equal `1`.  ∎

### 3.4 The case `m` divides `n`  (`r = 0`, `u = L/m >= 1`)

Here `[n]_q` vanishes at `zeta_m`, so use `Y_{n,k}(q) = qbinom(n-1,k-1)/[k]_q`. With `m | n`
and `gcd(k,m) = 1` we have `s >= 1`, `n-1 = m(a-1) + (m-1)`, `k-1 = mb + (s-1)`, and the
standard evaluation `qbinom(m-1,j)|_zeta = (-1)^j zeta^(-j(j+1)/2)` together with
`[k]_zeta = [s]_zeta` gives

    Y_{n,k}(zeta_m) = C(a-1, b) * sigma_s,     sigma_s = (-1)^(s-1) zeta^(-s(s-1)/2) / [s]_zeta.

Since `n_m = a - u`, Vandermonde `C(a-1,b) = sum_j C(u-1,j) C(a-u, b-j)` reduces (*) to

    [z^(jm+s)] Q  =  f * sigma_s * C(u-1, j)   for all j >= 0.

By 3.2 the achievable `s` in `[0, m-1]` are `s in {1, m-1}`. One computes `sigma_1 = 1`, and
using `[m-1]_zeta = -zeta^(-1)`, `sigma_(m-1) = (-1)^(m+1) zeta^(1 - (m-1)(m-2)/2) = 1` for
each of `m = 2,3,4,6`. So the requirement is `[z^(jm+s)] Q = f * C(u-1,j)`.

- `(L,m) in {(2,2),(3,3),(4,4),(6,6)}`: `u = 1`, so `C(0,j) = [j=0]`. For `j = 0`,
  `s in {1, m-1} = {1, L-1}` and `[z^s]Q = f` by (**). For `j >= 1`, `jm + s > L = deg Q`, so
  `[z^(jm+s)]Q = 0`. ✔
- `(L,m) = (4,2)`: `u = 2`, achievable `s = 1`, need `[z^(2j+1)]Q = f * C(1,j)`. `j = 0` gives
  `[z^1]Q = f` ✔; `j = 1` gives `[z^3]Q = [z^(L-1)]Q = f` ✔ by palindromicity; `j >= 2` gives
  `2j+1 > 4 = deg Q`, both sides `0` ✔.  ∎

So the entire proof is: *`Q` is palindromic of degree `L` with linear coefficient `f`, and the
achievable residues are exactly `{1, L-1}` — which is the statement `phi(L) <= 2`.*

*Verified:* the coefficient identity of 3.3/3.4 was checked in 1,048,609 instances
(`n <= 200`, every criterion-satisfying multiplier, every achievable `s`, every `j`), zero
failures. The end-to-end sieving identity `Y_{n,k}(zeta_m) = #Fix(v)` was checked
independently for every multiplier and every admissible `k` with `n <= 48` (11,578
instances), zero failures; `verify.py --full` extends this to `n <= 130`.

---

## 4. Hypothesis 3.6 is false

> **Hypothesis 3.6 (paper).** For `zeta` a primitive `m`-th root of unity, `2 <= r < m` and
> `s in {0,2,3,...,r-2,r}`, the number `qbinom(r,s)|_zeta / [r]_zeta` is irrational.

**It is not.** The smallest counterexample is `m = 9, r = 8, s = 4`, where the ratio is `-1`.
There are 25 counterexamples with `m <= 34` and 67 with `m <= 60`; every one has `r = m-1`
and value `±1`.

**Proposition D.**
1. *(Necessary condition, proved.)* If `x = qbinom(r,s)|_zeta / [r]_zeta` is rational then
   `m | r - 1 - s(r-s)`.
   *Proof.* `qbinom(r,s)_{q^-1} = q^{-s(r-s)} qbinom(r,s)_q` and `[r]_{q^-1} = q^{1-r}[r]_q`,
   so complex conjugation sends `x` to `zeta^(r-1-s(r-s)) x`. Also `x != 0` for `s <= r < m`,
   since `Phi_m` divides `qbinom(r,s)` only when `s > r`. ∎
2. *(The case `s in {0,r}`, proved.)* `[r]_zeta` is **irrational** for `2 <= r < m`.
   *Proof.* If `[r]_zeta = c in Q` then `c != 0` and conjugation gives `c = zeta^(1-r) c`, so
   `m | r-1`, impossible for `2 <= r < m`. ∎ Hence the paper's use of the hypothesis at
   `s in {0,r}` is sound.
3. *(The case `r = m-1`, proved.)* `qbinom(m-1,s)|_zeta / [m-1]_zeta = (-1)^(s+1)
   zeta^(1 - s(s+1)/2)`, of modulus `1`. It is rational exactly when
   `zeta^(1-s(s+1)/2) = ±1`, which happens for many `s` — these are all the counterexamples
   found.
4. *(What survives, conjecturally.)* For `2 <= r <= m-2` the ratio appears to be always
   irrational; verified for all `m <= 60` and all admissible `r,s`, using (1) to cut the
   search. This restricted statement is what the paper needs.

**Does the gap bite?** Yes. The paper's necessity proof invokes Hypothesis 3.6 whenever a bad
residue is achievable; that argument is invalid exactly when `r = n mod m = m-1` and the
residue is one of the rational ones. The first multiplier where this happens is

    n = 99,  v = 19   (m = 10, f = 9, uniform, r = m-1 = 9)

— just **outside** the paper's verified range `n <= 96`, which is why neither the proof nor
the computation caught it. Necessity nonetheless still holds there, but for a different
reason: at `k = 13` (so `s = 3`) the evaluation is

    Y_{99,13}(zeta_10) = -9        while the true fixed-necklace count is 84.

The polynomial evaluates to a perfectly rational — indeed **negative** — integer, and simply
gives the wrong answer.

**Repair.** In the uniform case with `r = f`, the sieving identity forces
`rho_s = C(f,s)/f`. For `s in {0,f}` use (2): `rho_s = 1/[f]_zeta` is irrational, done. For
`2 <= s <= f-2` (which needs `f >= 4`): if `rho_s` is irrational, done; and if `r = m-1` then
by (3) `|rho_s| = 1`, whereas `C(f,s)/f >= C(f,2)/f = (f-1)/2 > 1`, so equality fails. Hence
necessity is unconditional except for the residual case `2 <= r <= m-2` with `rho_s`
rational, which item (4) says is empty for `m <= 60`. A clean sufficient replacement for
Hypothesis 3.6 is the weaker claim **"if the ratio is rational then it has modulus 1"** —
true in every case computed, and proved above whenever `r = m-1`.

---

## 5. Status after this session

- **Proved:** Theorem A (unconditional, all contents); the structure theorem B; **sufficiency
  of the criterion in full generality** (Theorem C) — the paper's main open problem;
  Proposition D items 1–3.
- **Disproved:** Hypothesis 3.6 as stated.
- **Still open:** necessity in the general (non-uniform) multi-level case — verified only;
  the residual irrationality claim D(4); and the three-colour phenomenon at `n = 28`, on
  which Theorem A now gives direct leverage, since it holds for every content.

## 6. Files

    core.py       exact cyclotomics, Y(zeta_m), orbit/affine cycle machinery
    verify.py     ENTRY POINT.  `python3 verify.py` re-checks every claim above.

Everything is exact integer arithmetic; no floating point enters any decision.
