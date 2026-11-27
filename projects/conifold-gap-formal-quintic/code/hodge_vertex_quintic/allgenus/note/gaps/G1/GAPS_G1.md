# Gap cluster G1 (section s02 of the note): results, part 1 of 5

Author-facing; not part of the note. Task key `gaps/G1`. Object: the five items of cluster G1 (OUTLINE_3 s02; OUTLINE_REVIEW R1,
R2, R7; DECISIONS D2). Notation of `NOTE/NOTATION.md` (+`_2`); `NOTE` = `allgenus/note/`. Pieces: this file (s.0 summary, s.1 item 1);
`GAPS_G1_2.md` and `GAPS_G1_3.md` (s.2 item 2, `prop:set-cyclicLR`: inputs, lemmas, proof, corollaries, literature, evidence);
`GAPS_G1_4.md` (s.3 item 3, `thm:set-qrr`; s.4 item 4, `lem:set-dilaton`); `GAPS_G1_5.md` (s.5 item 5, the model; s.6 consequences
for the design; s.7 reading record). Evidence script and log:
`g1_cyclic_check.py`, `g1_cyclic_check.log` (this folder). Labels: PROVED, VERIFIED, CONJECTURAL, SKETCH, OPEN. AI-written; no
human has read it.

## 0. Summary (one resolution per item)

| item | resolution | where |
|---|---|---|
| 1. (M), genus-0 mirror theorem for the formal quintic | **CITED**: Givental, alg-geom/9603021v2 (IMRN 1996, root key MR1408320), Corollary 11.7 (with Theorem 11.6, Proposition 11.4 and the definitions of s.9), case `n = 4`, `r = 1`, `l_1 = 5`, `lambda' = 0`, then `lambda_i = zeta^i`. Read in the PDF. Hypotheses checked; translation of Givental's `s_i` to the twisted J-function written out (s.1.3). Coates-Givental Thm 2 is NOT a literal pin (non-equivariant `X`). | s.1 |
| 2. `prop:set-cyclicLR` (R1) | **PROVED** (all orders in `q`) from (M) and the genus-0 axioms of the twisted theory: Z/5 grading, unit, Frobenius property, Birkhoff recursion, the Wronskian of PF4 and the inhomogeneous period `PF4 c_4 = 1`. Also FOUND in the literature as a cited statement (GJR18 s.5.1 (13); Guo-Ross 1703.06955 Props 6.1-6.2 with Zagier-Zinger; LP18 (17), (19)-(20)). VERIFIED to `q^30` with failing controls. | GAPS_G1_2, GAPS_G1_3 |
| 3(a). passage `R(q=0) = diag(h_(w_i))` -> `A^QRR_n(q=0) = h_n` | **PROVED** (evaluation at `q = 0` of the defining relation of the flat amplitude; three lines once the definition of `Phi_i` is fixed). | GAPS_G1_4 s.3.1 |
| 3(b). uniqueness: the R-matrix at `q = 0` is the vertex R-matrix | **CITED** with the argument written: Teleman 0712.0160v3 Theorem (2) in (1.15), 8.7 Proposition (i)-(ii), 8.10 Proposition (read verbatim in the v3 PDF; 8.4 Proposition read as context), LP18 eq. (26) (read), GJR18 Thm 4.1 (read), PPZ15 Prop. 2.14 (read). New step PROVED here: `R(q=0)` is diagonal of the form `exp(odd)` by the QDE alone, so only rank-one uniqueness is needed. | GAPS_G1_4 s.3.2 |
| 4. `lem:set-dilaton` | **PROVED**, with pinned inputs: PPZ15 (1303.1043v3) Def. 0.3 (iii), Prop. 2.12, Def. 2.13; Arbarello-Cornalba alg-geom/9406008 formula (1.7) with `kappa_0 = 2g - 2 + n` (read). Hypothesis on `kappa` made explicit (`kappa` in `1 + m`, m-adically complete). | GAPS_G1_4 s.4 |
| 5a. symplecticity of `M_0` | **PROVED** (constancy of the pairing of a solution and an opposite solution of the self-adjoint model system; homogeneity). Not proved in the sources (LG_MODEL s.3: VERIFIED to `zhat^40`; UR s.7 table: asserted). | GAPS_G1_5 s.5.1 |
| 5b. homogeneity of the model CohFT | **PROVED** (degree count in the R-matrix action). Stated without proof in the sources (MSG1 s.1; RD1 s.1.2 "PROVED, LG_MODEL s.3, UR s.1", which prove only the homogeneity of the rows). | GAPS_G1_5 s.5.2 |

Nothing in the cluster is left OPEN. Two residual points are flagged honestly (not gaps of the chain, but places where the note
should say what it relies on): (i) Givental 1996's proof of Cor. 11.7 has a known history of being regarded as terse; the note may add
[LLY97] (not read here) as an independent proof (s.1.5); (ii) Teleman's injectivity "the CohFT determines its classification datum" is
part of his classification and not a separately numbered statement; in our case only the rank-one case is needed, and an elementary
route to it is recorded (s.3.2, step 5).

## 1. Item 1: the genus-0 mirror theorem (M) for the formal quintic

### 1.1 The statement the note uses (proposed wording for `thm:set-mirror`, name (M), R2)

Let `T = (C^*)^5` act on `P^4` with weights `lambda_0..lambda_4` on `H^0(P^4, O(1))`, `H = c_1^T(O(1))`, `H|_(p_i) = lambda_i`, and let
`E_(0,n,d) := R^0 pi_* f^* O(5)` on `Mbar_(0,n)(P^4, d)` (a vector bundle of rank `5d + 1`, since `R^1 pi_* f^* O(5) = 0` in genus 0).
The genus-0 twisted invariants are `<psi^(k_1) g_1, .., psi^(k_n) g_n>_(0,n,d) := Int_[Mbar_(0,n)(P^4,d)] e_T(E_(0,n,d)) prod ev_j^* g_j psi_j^(k_j)`,
the twisted pairing `eta(a, b) := Int_(P^4) a b e_T(O(5))`, and the J-function
`J(t, z) := z + t + Sum_(n,d) (Q^d/n!) Sum_a <t, .., t, H^a/(z - psi)>_(0,n+1,d) H_a^vee`, `H_a^vee` the `eta`-dual basis ([CG07, s.6,
(10)] normalisation, with the twisted pairing). Then, at the specialisation `lambda_i = zeta^i` (where `H^5 = 1`, `eta(H^a,H^b) =
5 delta_(a+b = 3 mod 5)`), with `Q = 1` and `t = T H`:

**(M)** `J(T(q) H, z) = I(q, z)/C_0(q)`, where `I|_(p_i) = I_i(q, z) = z q^(zeta^i/z) Sum_d q^d prod_(k=1..5d)(5 zeta^i + kz) / prod_j
prod_(k=1..d)(zeta^i - zeta^j + kz)` (`def:set-genus0`), `C_0 = Sum (5d)!/(d!)^5 q^d`, `T(q) = I_1/C_0 = log q + Itilde_1/C_0`.

(At generic weights the right side is `exp(-tau_0/z) I/C_0` with `tau_0` proportional to `e_1(lambda)`; `e_1(zeta^i) = 0`.)

### 1.2 The pinned statement (read in the PDF of alg-geom/9603021v2)

Givental, "Equivariant Gromov-Witten invariants", s.9 sets up (p.23-26): `Y_(2,d)` = genus-0 degree-`d` stable maps `(C, x_0, x_1) ->
CP^n`; `E_d` = Euler class of the bundle with fibre `H^0(C, phi^*(H^(l_1) + .. + H^(l_r)))`; `S_d(hbar) := E_d/(hbar + c_1^(0))`,
`c_1^(0)` the first Chern class of the universal TANGENT line at `x_0`; `S(t, hbar) := e^(Pt/hbar) Sum_d e^(dt) (e_0)_*(S_d(hbar))` (for
`d = 0`, `(e_0)_* S_0 := Euler(+_j H^(l_j))`); the equivariant counterpart `S'` in `T x T'`-equivariant cohomology (`T` the torus of
`C^(n+1)`, `T'` acting on the fibres of `+_j H^(l_j)` with weights `lambda'_j`; Euler class `(l_1 p - lambda'_1)..(l_r p - lambda'_r)`),
`phi_i := prod_(j != i)(p - lambda_j)`. p.27: `Z_i := Sum_(d>=0) q^d Int_(Y_(2,d)) e_0^*(phi_i) (hbar + c)^(-1) E'_d`, `E'_d` the Euler
class of the bundle of sections vanishing at `x_0`, and `<phi_i, S'> = e^(lambda_i t/hbar)(l lambda_i - lambda')(Z_i|_(q = e^t))` (`r = 1`).
p.40, **Theorem 11.6** (CY case `l_1 + .. + l_r = n + 1`): the hypergeometric `Z^*_i = Sum_d q^d prod_r prod_(m=1..l_r d)(l_r lambda_i -
lambda'_r + m hbar)/prod_(alpha=0..n) prod_(m=1..d)(lambda_i - lambda_alpha + m hbar)` coincides with `{Z_i}` after (1) `Q = q exp{Sum_a l_a
[g_(l_a)(q) - g_1(q)]/f(q)}`, (2) multiplication of `Z_i(Q(q), hbar)` by `exp{(f(q) hbar)^(-1)[Sum_a (l_a lambda_i - lambda'_a) g_(l_a)(q) -
(Sum_alpha (lambda_i - lambda_alpha)) g_1(q)]}`, (3) multiplication by `f(q)`; here `f(q) = Sum_d prod_a (l_a d)!/(d!)^(n+1) q^d`, `g_l(q) =
Sum_(d>=1) q^d prod_a (l_a d)!/(d!)^(n+1) Sum_(m=1..ld) 1/m`. p.41, **Corollary 11.7**: for the solutions `s_i = e^(lambda_i T/hbar)
Z_i(e^T, hbar)` "to the equivariant quantum cohomology differential equations", the operations (1) `T = t + Sum_a l_a [g_(l_a)(e^t) -
g_1(e^t)]/f(e^t)`, (2) multiplication by `f(e^t) exp{[g_1(e^t)(Sum_alpha lambda_alpha) - Sum_a lambda'_a g_(l_a)(e^t)]/(hbar f(e^t))}`
transform `{s_i}` to `s^*_i = e^(pt/hbar) Sum_d e^(dt) prod_a prod_(m=1..l_a d)(l_a p - lambda'_a + m hbar)/prod_alpha prod_(m=1..d)(p -
lambda_alpha + m hbar)|_(p = lambda_i)`, solutions of `prod_alpha(hbar d/dt - lambda_alpha) s^* = e^t prod_a prod_(m=1..l_a)(hbar l_a d/dt -
lambda'_a + m hbar) s^*`; "For `lambda' = 0`, `lambda_0 + .. + lambda_n = 0` the solutions `s^*_i` have the following integral
representation". (The proof rests on Prop. 11.4 (3)-(5), p.37-39, and the polynomiality of Prop. 11.2-11.3, p.35-36.)

### 1.3 Translation of conventions (written out)

Case of the formal quintic: `n = 4`, `r = 1`, `l_1 = 5` (`l_1 = n + 1`: Givental's CY case), `lambda'_1 = 0` (LP18 (10): the fibre weight of
`O(5)` at `p_i` is `5 lambda_i`; GJR18 s.3.1 (9): twist parameter `t = 0`), `p = H`, `hbar = z`, `e^t = q`.
1. *`f` and the mirror map.* `f = Sum (5d)!/(d!)^5 q^d = C_0`; `5(g_5 - g_1) = 5 Sum q^d ((5d)!/(d!)^5) Sum_(m=d+1..5d) 1/m = Itilde_1`, so
   operation (1) is `T = t + Itilde_1/C_0 = I_1/C_0`, the mirror map of `def:set-genus0`.
2. *The factor.* At `lambda' = 0` the factor (2) is `C_0 exp{e_1(lambda) g_1/(z C_0)}`; at `lambda_i = zeta^i`, `e_1 = 0`, so it is `C_0`.
3. *`s^*_i` is the I-function.* `s^*_i(t) = q^(lambda_i/z) Sum_d q^d prod_(m=1..5d)(5 lambda_i + mz)/prod_alpha prod_(m=1..d)(lambda_i -
   lambda_alpha + mz) = I_i(q,z)/z`.
4. *`s_i` is the J-function.* By the definition of the twisted pairing, the `eta`-dual push-forward of a class `beta` on
   `Mbar_(0,1)(P^4,d)` is `ev_*(beta)/e_T(O(5))`; since `E_(0,1,d) = ev^* O(5) + E'_d` (sections vanishing at the marking; exact
   sequence of vector bundles), `ev_*(e(E_(0,1,d)) beta)/e(O(5)) = ev_*(e(E'_d) beta)`. With the divisor equation, `J(TH, z) =
   e^(TH/z)(z + Sum_(d>=1) e^(dT) ev_*(e(E'_d)/(z - psi)))`. In Givental's `S'` the class `E_d/(hbar + c)` lives on the two-pointed
   `Y_(2,d)`, `c = -psi_0`, nothing inserted at `x_1`; pushing forward along the map forgetting `x_1` (`E_d` is pulled back from
   `Mbar_(0,1)`) and using `pi_*(psi_0^k) = psi_0^(k-1)` (string equation for stable maps, `d >= 1`: `psi_0 = pi^* psi_0 + D_(0,1)`,
   `D_(0,1)` the locus where `x_0, x_1` lie on a contracted 3-pointed component; the same proof as Arbarello-Cornalba (1.9) for
   `Mbar_(g,n)`, read) gives `pi_*(1/(hbar - psi_0)) = hbar^(-1)(hbar - psi_0)^(-1)`. Hence `S'/e(E) = J(tH, hbar)/hbar` (the `d = 0` term is `e^(pt/hbar)` by Givental's convention), and
   `s_i = (S'/e(E))|_(p_i) = J(TH, z)|_(p_i)/z`, using `<phi_i, S'> = S'|_(p_i)` (`phi_i|_(p_i) = e(T_(p_i) P^4)`).
5. *Conclusion.* Corollary 11.7 reads `C_0(q) J(T(q)H, z)|_(p_i)/z = I_i(q,z)/z` at `lambda' = 0` and `e_1(lambda) = 0`, i.e. (M).

### 1.4 Hypotheses (checked)

* *CY case:* `l_1 = 5 = n + 1`. Theorem 11.6 and Cor. 11.7 are stated exactly for `l_1 + .. + l_r = n + 1`.
* *Equivariance:* Givental works `T x T'`-equivariantly with `T` the full torus of `C^(n+1)` (p.26); the formal quintic is the
  specialisation `lambda' = 0`, explicitly allowed (the integral-representation sentence of Cor. 11.7 and the proof of Cor. 11.8 use
  `lambda' = 0`); every class involved is polynomial in `(lambda, lambda', hbar)` (Prop. 11.2-11.4), so `lambda' = 0` is a specialisation.
* *Convexity (what is needed):* `E_d` must be a vector bundle on `Y_(2,d)` (Givental's `E_d` is the Euler class of the bundle with fibre
  `H^0`). This holds: `H^1(C, f^* O(5)) = 0` for every genus-0 stable map (degree `>= 0` on each component). Consequently the formal
  quintic's genus-0 Euler class "defined after localisation" (LP18 s.0.2, `c_top(R^0)/c_top(R^1)` on fixed loci) is the honest
  equivariant Euler class of `E_(0,n,d)` (`R^1 = 0`), and the twisted genus-0 invariants are polynomials in `lambda`. In genus 0 nothing
  beyond convexity is used; in higher genus (not needed for (M)) the definition after localisation needs `R^1` to carry only moving
  weights on every fixed locus: `R^1` comes from the contracted components, `E^vee (x) O(5)|_(p_i)`, whose weights are `5 lambda_i` shifted
  by nilpotent Hodge classes, so `c_top(R^1)` is invertible (`5 zeta^i != 0`). (The edge weights `(5d lambda_i + k(lambda_j - lambda_i))/d` of
  `R^0` sit in the numerator; they do not vanish at the zeta point either, since `zeta^(j-i)` is not real for `j != i`.)
* *Invertibility:* `e_T(O(5)) = 5H` is invertible after localisation (`5 lambda_i != 0`); this is what makes the twisted pairing
  non-degenerate (at the zeta point `eta(H^a,H^b) = 5 delta_(a+b = 3)`). It is not used by the identity of Cor. 11.7 itself.
* *Specialisation `lambda_i = zeta^i`:* both sides of (M), coefficientwise in `q` and `1/z` (after removing `e^(lambda_i T/z)`), are
  polynomials in `lambda` (left: equivariant push-forwards of genuine classes on proper moduli spaces; right: expand `1/(lambda_i -
  lambda_j + mz)` in `1/z`); `T(q)`, `C_0` do not depend on `lambda`. So evaluation at `lambda_i = zeta^i` is a ring homomorphism.

### 1.5 The other candidates (read) and why they are not the pin

* **Coates-Givental math/0110142v2, Theorem 2** (s.7, p.9, read): for `E` a sum of line bundles with integer Chern roots, the
  hypergeometric modification `I_E` of `J_X` lies on the cone `L_(e,E)` of the theory twisted by the `S^1`-equivariant Euler class
  `e(E) = prod (lambda + rho_i)` (`S^1` acting on the fibres only, s.1 p.4; `X` a compact projective manifold, no torus on `X`); no
  convexity hypothesis. Corollaries 1-3 (s.4-s.5, read) are quantum Riemann-Roch for the Euler class and quantum Serre; they do not
  contain the T-equivariant statement either. The formal quintic needs `X = P^4` with its `T`-action and `O(5)` with the induced
  linearisation; the `T`-equivariant version of Thm 2 is standard (the proof of CG Thm 1 is Grothendieck-Riemann-Roch, valid
  equivariantly) but it is not stated in CG. Hence CG Thm 2 is the modern formulation, not a literal pin. The project's earlier
  pointer "math/0110142 Thm 2, Cor 5.1" (PLAN_2026-11-03 "Key facts") is not accurate: CG has no Corollary 5.1.
* **LP18 = arXiv:1803.01409v2** (read, s.0.2-s.2.4, p.3-19): the formal quintic is defined (5) with the Euler class after localisation;
  s.2.4.3 (17) computes the stable-quotient S-operator `S(H^j)` from the I-function "using Birkhoff factorization, ... see [18]"
  (Kim-Lho), and (19)-(20) define `C_0..C_4`; the GW-side genus-0 statement enters only through the wall-crossing [7], [8]
  (Ciocan-Fontanine-Kim; not read). No GW-side mirror theorem is stated for the formal quintic.
* **GJR18 = arXiv:1812.11908v1** (read, s.3.1, s.4.1, s.5.1-5.2): s.5.2 (22) "by the genus zero mirror theorem `J_t(H_delta) =
  H_delta S^*(H_delta) 1/I_0`", without a pinned statement (their bibliography lists [14] = Givental alg-geom/9701016, [15] =
  Givental alg-geom/9603021, [28] = Lian-Liu-Yau).
* **Guo-Ross arXiv:1703.06955v1** (Compositio 155 (2019); read s.2.1-2.3, s.6.1-6.3): Theorem 2.1 is the non-equivariant quintic
  ("Givental [Giv98], Lian-Liu-Yau [LLY97]"); for the twisted (formal-quintic) theory, Prop. 6.1 ("standard properties of Givental's
  Lagrangian cone") uses (M) without a pin.
* **Lian-Liu-Yau, Mirror principle I** (Asian J. Math. 1 (1997)): NOT read. It is the standard second reference; a reader should
  confirm that its main theorem is stated `T`-equivariantly before it is cited for (M).

### 1.6 Recommendation for the note (input (M); decision D2)

1. Print (M) as in s.1.1 (name (M), never (S1); R2), cite "[MR1408320, Corollary 11.7 (with Theorem 11.6), case `n = 4`, `r = 1`,
   `l_1 = 5`, `lambda' = 0`]" and add one paragraph with the translation of s.1.3 (items 1-5) and the hypothesis check of s.1.4
   (convexity of `O(5)` in genus 0; polynomiality; `e_1(zeta^i) = 0`). Numbers are those of the arXiv v2 PDF (the only version;
   v1 withdrawn); `\draftnote{IMRN numbering of Givental 1996 not compared}` at the first citation, as for Teleman.
2. Mention [MR2276766, Theorem 2] as the modern formulation (non-equivariant in `X`) and, if the authors wish, LLY97 as an
   independent proof (`\draftnote{LLY97 not read by the gap phase}`).
3. The genus-0 axioms used together with (M) in `prop:set-cyclicLR` (string, divisor, QDE, unit, Frobenius; GAPS_G1_2 s.2.1) are
   properties of the twisted genus-0 CohFT; state them once next to (M) (sources there).
4. Replace "(S1)" by "(M)" in OUTLINE s.0.2, `sec:int-inputs`, `thm:set-mirror`, TITLE_ABSTRACT, LABELS l.12 and l.32 (R2).
