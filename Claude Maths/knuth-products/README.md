# Knuth products in linear numeration systems

**What is proved here.** Knuth's *Fibonacci multiplication* — the surprising fact that

$$n \circ m \;=\; \sum_{i,j} F_{k_i+l_j}, \qquad n=\sum_i F_{k_i},\; m=\sum_j F_{l_j} \ \text{(Zeckendorf)}$$

is associative — is an instance of one exact identity that holds in **every** linear numeration
system, with no hypotheses at all:

$$\boxed{\,n \circ_s m \;=\; L\!\left(\beta^{s}\,\Phi(n)\,\Phi(m)\right)\,}$$

where $\Phi$ reads a greedy expansion as an element of $\mathbb{Z}[\beta]$ and $L$ is the unique
$\mathbb{Z}$-linear functional $\mathbb{Z}[\beta]\to\mathbb{Z}$ with $L(\beta^k)=u_k$. Associativity
is then *exactly* the statement that the image $S=\Phi(\mathbb{N})$ is closed under multiplication
by $\beta^{s}$ — and the failure of associativity is measured by a closed-form obstruction
(Theorem 2(iv)).

Applied to systems with **non-canonical initial conditions** — a case the literature does not
cover — this produces phenomena nobody would guess:

* the set of good shifts **need not be an upward-closed set**. For the Lucas system
  $U=(1,3,4,7,11,18,\dots)$ the product $\circ_s$ is associative for $s=2$, **fails for $s=3$**,
  and is associative again for every $4\le s\le 14$ (verified). Minimal certificate of failure:
  $(1\circ_3 2)\circ_3 6 = 268 \ne 278 = 1\circ_3(2\circ_3 6)$.
* a non-canonical system can be **strictly better** than the canonical one: the natural Pell
  system $(1,2,5,12,29,\dots)$ is associative from $s=1$, the canonical Pell system
  $(1,3,7,17,41,\dots)$ only from $s=2$.

**Paper:** `paper/paper.pdf` (8 pp.) — build with `pdflatex paper.tex` (run twice).
**Verification:** `python3 verify.py` (exact integer arithmetic, no floating point, ~4 min).

---

## 1. Setup

Let $\chi(x)=x^d-a_1x^{d-1}-\cdots-a_d\in\mathbb{Z}[x]$ be the **minimal polynomial** of a real
number $\beta>1$, and let

$$U=(u_k)_{k\ge0},\qquad u_0=1,\qquad u_0<u_1<u_2<\cdots,\qquad
u_k=a_1u_{k-1}+\cdots+a_du_{k-d}\ (k\ge d).$$

The initial values $u_0,\dots,u_{d-1}$ are arbitrary subject to $u_0 = 1$ and monotonicity; taking
the *canonical* ones $u_k=a_1u_{k-1}+\cdots+a_ku_0+1$ gives the classical numeration system of a
simple Parry number, but nothing below needs that choice.

Because $u_0=1$, the **greedy expansion** of $n\in\mathbb{N}$ — repeatedly subtract the largest
$u_k\le n$ — terminates and yields a finite multiset $K(n)$ of indices with
$n=\sum_{k\in K(n)}u_k$. Define

$$\Phi(n)=\sum_{k\in K(n)}\beta^{k}\in\mathbb{Z}[\beta],\qquad S:=\Phi(\mathbb{N}),$$

$$n\circ_s m \;=\; \sum_{k\in K(n)}\sum_{l\in K(m)} u_{k+l+s}\qquad (s\ge 0).$$

For $U=(1,2,3,5,8,\dots)$ (Zeckendorf, $u_k=F_{k+2}$) and $s=2$ this is exactly Knuth's product,
since $k\mapsto k+2$ converts his Fibonacci indices to ours. For $u_k=b^k$ one gets
$n\circ_s m=b^s\,nm$.

## 2. Theorem 1 (transport)

> **Theorem 1.** There is a unique $\mathbb{Z}$-linear map $L:\mathbb{Z}[\beta]\to\mathbb{Z}$ with
> $L(\beta^k)=u_k$ for all $k\ge0$, and
> 1. $L(\Phi(n))=n$ for every $n\in\mathbb{N}$;
> 2. $n\circ_s m=L\!\big(\beta^{s}\Phi(n)\Phi(m)\big)$ for all $n,m\in\mathbb{N}$, $s\ge0$.

*Proof.* $\chi$ is monic and irreducible, so $\mathbb{Z}[\beta]\cong\mathbb{Z}[x]/(\chi)$ is a free
$\mathbb{Z}$-module with basis $1,\beta,\dots,\beta^{d-1}$. Define $L(\beta^i)=u_i$ for $i<d$ and
extend $\mathbb{Z}$-linearly; this is forced, hence unique. For $k\ge d$ multiply $\chi(\beta)=0$ by
$\beta^{k-d}$ to get $\beta^{k}=a_1\beta^{k-1}+\cdots+a_d\beta^{k-d}$; by induction
$L(\beta^k)=a_1u_{k-1}+\cdots+a_du_{k-d}=u_k$.

(1) is immediate: $L(\Phi(n))=\sum_{k\in K(n)}u_k=n$. For (2),
$\beta^s\Phi(n)\Phi(m)=\sum_{k\in K(n),\,l\in K(m)}\beta^{k+l+s}$; apply $L$ and use
$L(\beta^{k+l+s})=u_{k+l+s}$. $\blacksquare$

Two consequences worth stating separately.

* $\Phi:\mathbb{N}\to S$ is a **bijection** with inverse $L|_S$ (by (1)). So $S$ — call its elements
  the **$U$-integers** — is a set of representatives, and $n\circ_s m$ is nothing but the ring
  product $\beta^s\Phi(n)\Phi(m)$ *pushed back to $\mathbb{N}$ by $L$*.
* Since the trace form of a separable extension is non-degenerate, for $d\ge1$ there is a unique
  $\lambda\in\mathbb{Q}(\beta)$ with $L(z)=\operatorname{Tr}_{\mathbb{Q}(\beta)/\mathbb{Q}}(\lambda z)$.
  Choosing the initial conditions of $U$ *is* choosing $\lambda$: the Zeckendorf system is
  $\lambda=\varphi^2/\sqrt5$, the Lucas system $(1,3,4,7,\dots)$ is $\lambda=\varphi$.

## 3. Theorem 2 (associativity is closure, and the exact obstruction)

Fix $s\ge0$ and define the **defect**

$$\delta(n,m):=\Phi(n\circ_s m)-\beta^{s}\Phi(n)\Phi(m)\in\mathbb{Z}[\beta].$$

> **Theorem 2.**
> 1. $\delta(n,m)\in\ker L$ for all $n,m$.
> 2. $\beta^{s}S\cdot S\subseteq S \iff \delta\equiv 0$.
> 3. If $\beta^{s}S\cdot S\subseteq S$, then $\Psi:=\beta^{s}\Phi$ is injective and satisfies
>    $\Psi(n\circ_s m)=\Psi(n)\Psi(m)$. Hence $(\mathbb{N},\circ_s)$ is isomorphic to the
>    multiplicative semigroup $\Psi(\mathbb{N})\subset(\mathbb{Z}[\beta],\times)$; in particular
>    $\circ_s$ is associative.
> 4. $\circ_s$ is associative **iff**
>    $$L\big(\beta^{s}\,\delta(n,m)\,\Phi(p)\big)=L\big(\beta^{s}\,\Phi(n)\,\delta(m,p)\big)
>      \qquad\text{for all } n,m,p\in\mathbb{N}.$$

*Proof.* (1) $L(\Phi(n\circ_sm))=n\circ_s m=L(\beta^s\Phi(n)\Phi(m))$ by Theorem 1.

(2) ($\Rightarrow$) Put $z=\beta^s\Phi(n)\Phi(m)\in S$, say $z=\Phi(N)$. Then $L(z)=N$, while
$L(z)=n\circ_sm$ by Theorem 1(2); so $N=n\circ_sm$ and $\Phi(n\circ_sm)=z$, i.e. $\delta(n,m)=0$.
($\Leftarrow$) $\beta^s\Phi(n)\Phi(m)=\Phi(n\circ_sm)\in S$.

(3) $\Psi(n\circ_sm)=\beta^{s}\Phi(n\circ_sm)=\beta^{s}\cdot\beta^{s}\Phi(n)\Phi(m)=\Psi(n)\Psi(m)$,
and $\Psi$ is injective because $\Phi$ is and $\beta^s$ is a unit in $\mathbb{Q}(\beta)$. Therefore
$\Psi\big((n\circ_sm)\circ_sp\big)=\Psi(n)\Psi(m)\Psi(p)=\Psi\big(n\circ_s(m\circ_sp)\big)$, and
injectivity gives associativity.

(4) By Theorem 1(2) and the definition of $\delta$,
$$(n\circ_sm)\circ_sp=L\big(\beta^{s}\Phi(n\circ_sm)\Phi(p)\big)
 =L\big(\beta^{2s}\Phi(n)\Phi(m)\Phi(p)\big)+L\big(\beta^{s}\delta(n,m)\Phi(p)\big),$$
and symmetrically
$n\circ_s(m\circ_sp)=L\big(\beta^{2s}\Phi(n)\Phi(m)\Phi(p)\big)+L\big(\beta^{s}\Phi(n)\delta(m,p)\big)$.
Subtract. $\blacksquare$

Part (4) is the useful form: the two "extra" terms are *equal as integers* far more often than the
defects vanish — a priori. Experimentally they never are (§5, Observation 4), which says Knuth-type
associativity is never an accident.

> **Lemma 3.** If $\beta^{s}S\cdot S\subseteq S$ then $\beta^{s}S\subseteq S$, and consequently
> $\beta^{ks}S\cdot S\subseteq S$ for every $k\ge1$: the set of closure shifts is closed under
> multiplication by positive integers.

*Proof.* $\Phi(1)=\beta^{0}=1\in S$, so $\beta^{s}S=\beta^{s}S\cdot\{1\}\subseteq\beta^sS\cdot S
\subseteq S$. Iterating, $\beta^{js}S\subseteq S$ for all $j\ge0$, whence
$\beta^{(j+1)s}S\cdot S=\beta^{js}\big(\beta^{s}S\cdot S\big)\subseteq\beta^{js}S\subseteq S$.
$\blacksquare$

Lemma 3 is exactly why a *gap* in the good-shift set can only ever occur at non-multiples of a good
shift — see Observation 1.

**Knuth's theorem in this language.** For the Zeckendorf system, $\beta=\varphi$ and
$S=\{\sum_k d_k\varphi^k : d\ \text{has no two consecutive }1\text{'s}\}$ is the set of non-negative
golden-mean integers $\mathbb{Z}_\varphi^{\ge0}$. Knuth's theorem is the closure statement
$\varphi^{2}\,\mathbb{Z}_\varphi^{\ge0}\!\cdot\!\mathbb{Z}_\varphi^{\ge0}\subseteq
\mathbb{Z}_\varphi^{\ge0}$, i.e. "a product of two golden-mean integers has at most two golden-mean
fractional digits."

## 4. What is classical and what is new

*Classical.* Knuth introduced $\circ$ and proved associativity for Zeckendorf
(*Fibonacci multiplication*, Appl. Math. Lett. 1 (1988) 57–60). Arnoux gave the $\mathbb{Z}[\varphi]$
explanation (*Some remarks about Fibonacci multiplication*, Appl. Math. Lett. 2 (1989) 319–320).
The shifted product for general recurrences **with canonical initial values**, associative once the
shift is large enough, is Grabner et al., *Associativity of recurrence multiplication*, Appl. Math.
Lett. 7 (1994) 85–90 — that is exactly the canonical case of Theorem 2(3) plus the fact (Remark in
§3 of the paper) that there the admissible shifts form an interval. Messaoudi generalised the
construction to Pisot numeration systems (*Généralisation de la
multiplication de Fibonacci*, Math. Slovaca 50 (2000) 135–148; *Tribonacci multiplication*, Appl.
Math. Lett. 15 (2002) 981–985). For the **canonical** system of a simple Parry number, greedy
expansions are the $\beta$-admissible strings (Bertrand), so $S=\mathbb{Z}_\beta^{\ge0}$ and the
closure condition of Theorem 2 is exactly the finiteness of the quantity $L^{\otimes}$ studied by
Burdík–Frougny–Gazeau–Krejcar and Ambrož–Frougny–Masáková–Pelantová (*Arithmetics on number systems
with irrational bases*, Bull. Belg. Math. Soc. 10 (2003)): $s$ is good iff $s\ge L^{\otimes}(\beta)$.

*New here.* (a) Theorem 1 with **arbitrary** initial conditions — the numeration system is decoupled
from $\beta$, and the functional $L=\operatorname{Tr}(\lambda\,\cdot)$ becomes a free parameter, so
$S$ ranges over a whole family of "$U$-integer" sets rather than only $\mathbb{Z}_\beta$.
(b) The exact obstruction identity, Theorem 2(4). (c) The phenomena in §5, which are invisible in the
canonical setting: there $\{s:\text{good}\}=[L^{\otimes},\infty)$ is an interval by construction.

## 5. Phenomena (all machine-verified by `verify.py`)

Write $G(U)=\{s:\circ_s\ \text{is associative}\}$.

| system | $u_0,u_1,\dots$ | $\beta$ | $G(U)\cap[0,14]$ |
|---|---|---|---|
| base 2 | 1,2,4,8,… | 2 | 0,1,2,3,4,… |
| Zeckendorf | 1,2,3,5,8,… | $\varphi$ | 2,3,4,5,… (Knuth) |
| **Lucas** | **1,3,4,7,11,…** | $\varphi$ | **2, 4,5,6,…**  ← 3 missing |
| **Fib-shift** | **1,4,5,9,14,…** | $\varphi$ | **4, 6,7,8,…**  ← 5 missing |
| Pell canonical | 1,3,7,17,41,… | $1+\sqrt2$ | 2,3,4,… |
| **Pell natural** | **1,2,5,12,29,…** | $1+\sqrt2$ | **1**,2,3,… ← beats canonical |
| $x^2=2x+2$ canon. | 1,3,8,22,… | $1+\sqrt3$ | 4,5,6,… |
| Tribonacci | 1,2,4,7,13,… | 1.8392… | 3,4,5,… |
| $x^3=x^2+1$ canon. | 1,2,3,4,6,9,… | 1.4655… | 7,8,9,… |
| $x^2=x+3$ (not Pisot) | 1,2,5,11,… | 2.3027… | $\varnothing$ (checked $s\le20$) |
| $x^4=x^3+1$ (Pisot, no (F)) | 1,2,3,4,5,7,… | 1.3802… | $\varnothing$ (checked $s\le20$) |

**Observation 1 — $G(U)$ need not be upward closed.** Lucas: $2,4\in G$ but $3\notin G$, certified by
$(1\circ_3 2)\circ_3 6=268\neq278=1\circ_3(2\circ_3 6)$. Fib-shift: $4,6\in G$ but $5\notin G$,
certified by $(1\circ_5 8)\circ_5 8=8017\neq8072$. Both certificates are finite exact computations,
hence rigorous.

*Why a gap is possible.* Lemma 3 forces every **multiple** of a good shift to be good, which is why
all even shifts are good for the Lucas system once $2$ is — but $3$ is not a multiple of $2$, and
nothing propagates closure to it. What would propagate it is $\varphi S\subseteq S$, and that fails:
in the Lucas system $2=\Phi(2)\in S$ (greedy: $2=u_0+u_0$) while $2\varphi\notin S$, because
$L(2\varphi)=2u_1=6$ and $\Phi(6)=\varphi^{2}+2\varphi^{0}\ne2\varphi$ (greedy: $6=u_2+u_0+u_0$).
The odd shifts $5,7,9,\dots$ are nevertheless good, so Lemma 3 does not explain all of $G$; the
combinatorial reason those survive while $3$ does not is open.

**Observation 2 — $G(U)$ is not determined by $\beta$.** Three systems above share $\beta=\varphi$
and have $G=[2,\infty)$, $\{2\}\cup[4,\infty)$, $\{4\}\cup[6,\infty)$. Likewise the two Pell systems
differ. So the good-shift set is a genuine invariant of the *pair* $(\beta,\lambda)$, not of $\beta$.

**Observation 3 — a non-canonical system can beat the canonical one.** $G$ for natural Pell starts at
$1$, one better than canonical Pell; equivalently, $\lambda$ can be chosen so that the $U$-integers
are closed under a smaller power of $\beta$ than the $\beta$-integers are.

**Observation 4 — no accidental associativity.** In all 13 systems and all shifts $s\le8$ tested,
$$\circ_s\ \text{associative}\iff \beta^{s}S\cdot S\subseteq S .$$
Theorem 2(3) proves $\Leftarrow$. The converse is *not* proved here; by Theorem 2(4) it amounts to
showing that the bilinear obstruction cannot vanish identically without $\delta$ itself vanishing.
**Conjecture:** it always does — i.e. every associative Knuth product is a transported ring product.

**Observation 5 — two ways to have $G=\varnothing$.** $\beta$ not Pisot ($x^2=x+3$) makes even the
$\beta$-expansions of products badly behaved; but $\beta$ Pisot is not sufficient: $x^4=x^3+1$ is
Pisot (the second-smallest Pisot number, 1.3802…) yet $G=\varnothing$ up to $s=20$, in accordance
with its digit string $d_\beta(1)=1001$ violating the Frougny–Solomyak descending-digit condition for
the finiteness property (F).

## 6. Files

* `paper/paper.tex`, `paper/paper.pdf` — the formal write-up (amsart, 8 pp.).
* `verify.py` — exact verification of everything above: irreducibility of the $\chi$'s, Theorem 1,
  the equivalence of Theorem 2(2) with the morphism property, the good-shift tables, the minimal
  counterexamples, Knuth's theorem as a control, and the long-range scans. Exits with
  `ALL CHECKS PASSED`.
