# DIGEST_L_4: s.2 the novelty statement; s.4 conventions of the literature

Continuation of `DIGEST_L_3.md`. Keys as in s.3 (`DIGEST_L_5.md`). Plain prose drafts below are suggestions for the writers, built only
from the claims L1-L48; every mathematical sentence in them is a transcription of a project source named in brackets.

## 2. The novelty statement

### 2.1 What the introduction may state
* **Main theorem** (transcribe from MECHANISM_STAGE5 s.2 and NATURAL_FORM s.1 item 7): for every `g >= 2`, in the natural form
  `M = M_0(tau(u) zhat)`, `R' = Psi(u) G`, of the Hodge-vertex Givental-Teleman graph sum of the formal quintic near `u = 0`, every
  stable graph of genus `g` with at least one edge is regular at `u = 0`, and the vertex graph is `c_g t_c^(2-2g)` + regular,
  `c_g = (-1)^(g-1) B_2g / (2g(2g-2)) 5^(1-g)`. Conditional on (Q) and (QRR) only (rule 8); Theorem U is proved in the note.
* **Corollary (conifold gap of the formal quintic)**: summing over the leg-free graphs, the genus-`g` potential of the formal quintic,
  expanded at the conifold point in a conifold frame, is `c_g t_c^(2-2g) + O(1)`, for every `g >= 2`. Scope: the conifold frames of the
  sources (`varpi_0^cf = f_0 + alpha g + beta f_2`, FORMAL_QUINTIC_genus2 s.1; the same family as CGYZ Remark 1.4). This corollary
  is "the gap of the total" (MECHANISM s.2); its one-step derivation from N(g) and the identification of the total graph sum with the
  formal quintic's `F_g` must be written in the note from the sources (s.7 item 4).
* **Relation to the literature** (L5, L8, L9): the gap for the formal quintic is open in the literature (CGYZ26 Remark 1.9: expected,
  not proved; GJR18 Remark 1.6 for the quintic); it is one of the two key statements of the approach of Guo-Janda-Ruan to the quintic's
  gap (CGYZ26 Remark 1.9). The quintic's gap itself is proved, by a different route, in CGYZ26 (Theorem 2).
* **The constant**: `c_g` equals the constant of CGYZ26 Conjecture 1.5 for the true quintic, in the same normalisation of the vanishing
  period (`tau_c = Delta + (7/10) Delta^2 + O(Delta^3)`, CGYZ26 Example 1.3; `t_c = u + (7/10) u^2 + ...`, FORMAL_QUINTIC_genus2 s.1).
  This is a comparison of two printed formulas; state it as a remark, with a `\draftnote` asking the authors to compare the definitions
  of the conifold expansion (CGYZ26 Defs 1.1-1.2 against the sources' map "Cf", FORMAL_QUINTIC_genus2 s.1) before submission.
* **Method, and what is new in it** (L43-L48): the factorisation through the Landau-Ginzburg model `x^4 + v/x` with primitive form
  `dx`; Theorem Gamma (no pole at `L = 1`, value `h_(w_0)`), new as a written proof as far as the literature read shows; `Ginf` is
  GJR18 Lemma 7.7 (cite as such); (QRR) quoted.

### 2.2 Draft positioning paragraphs (plain text, to be written in LaTeX by the writers)
1. "The conifold gap is a prediction of the topological B-model: near a conifold point the genus-g potential, expanded in the flat
   coordinate given by the vanishing period, has the universal leading pole of the c = 1 string at the self-dual radius and no other
   pole [Bershadsky:1993cx; Ghoshal:1995wm; HK07; HKQ09, (1.2)]. For the quintic threefold it is the one of the four structural
   predictions listed in [GJR18, abstract] that remained open after [GJR18; CGL21; CGL26] (see [GJR18, Remark 1.6]); it is an input
   of the computation of invariants of high genus [HKQ09; LR22]. It is elementary for the resolved conifold [Br25, footnote 1], it was proved for local
   P^2 in [Br25] (after a check to genus 7 in [CI21]), for the topological recursion of toric mirror curves in [FCCZ26], and for the
   quintic and the other one-parameter Fermat hypersurfaces in [CGYZ26, Theorem 2], by master-space localisation."
2. "The formal quintic [LP18, (5)] is the O(5)-twisted, torus-equivariant Gromov-Witten theory of P^4 specialised at lambda_i =
   zeta^i [ZZ08]. It satisfies the holomorphic anomaly equations of the quintic [LPfq, Theorems 1-2], and its constants of integration
   are computable by localisation [LPfq, s.0.4]. It is the main tool of the approach of Guo-Janda-Ruan to the quintic's gap, which
   reduces it to the gap for the formal quintic and to the conifold regularity of certain of its correlators [CGYZ26, Remark 1.9]. The
   gap for the formal quintic has not been proved; the authors of [CGYZ26, Remark 1.9] expect their method to give it."
3. "We prove, conditionally on (Q) and (QRR), a graph-local refinement for the formal quintic: [main theorem]. Summing over graphs
   gives the gap for the formal quintic in every genus, with the constant c_g, which agrees with the constant of [CGYZ26, Conjecture
   1.5] for the quintic in the same normalisation."
4. "The graph sum is the Givental-Teleman sum of the semisimple formal quintic CohFT [MR2917177; PPZ15; GJR18, Theorem 4.1],
   normalised at q = 0 by quantum Riemann-Roch [Giv01, Theorem 1; MR2276766, Corollary 1; GJR18, (17)]. Near the conifold point its
   R-matrix factorises through the R-matrix of the Landau-Ginzburg model x^4 + v/x with primitive form dx, a Hurwitz Frobenius manifold
   whose CohFT is computed by topological recursion on that curve [Dunin-Barkowski:2015caa, Theorem 8.1]; topological recursion is used
   only for this model curve [Eynard:2007kz; ABDKS25; Hock25], never for a mirror curve of the quintic. The normalisation at q = 0 is
   carried to the conifold point by a formal stationary-phase identity for a Mellin-Barnes integral (Theorem Gamma). Its mechanism, a
   change of chart that makes the integrand factorise into Euler integrals, is Givental's for projective spaces [MR1901075, s.10]; an
   equivalent identity for the Landau-Ginzburg integral of the formal quintic is stated without proof in [GJR18, Remark 7.2], and the
   value at the orbifold point is [GJR18, Lemma 7.7]. For Gamma factors in quantum cohomology see also [MR2553377; GGI16]."
(Paragraph 4's mathematical content: MECHANISM_STAGE5 s.1, s.3.1, s.3.6; NATURAL_FORM s.0; STAGE5_SYNTHESIS_2 s.8 item 2.)

### 2.3 What the introduction must not claim
1. The gap for the true quintic, or any statement about the quintic's own graph sums (CGYZ26 proves the quintic's gap; the note
   treats the formal quintic only).
2. That a gap has never been obtained from R-matrix methods: CGYZ26's proof factorises a master-space R-matrix, `R = R^reg R^prin`,
   with vanishing non-leading graph contributions (Prop. 3.1, s.3 opening). Say instead that the note factorises the R-matrix of the
   formal quintic itself, through a universal model. Brini's sentence (Br25 s.1.1) may be quoted only as his remark.
3. GJR's second key statement (conifold regularity of correlators with special insertions, CGYZ26 Remark 1.9): its precise form is
   unpublished; do not claim to prove it, and do not identify it with the note's Model Lemma.
4. That Lho-Pandharipande conjectured the gap for the formal quintic (they do not discuss the conifold point: L15).
5. Anything unconditional: (QRR) is quoted (L27), not proved; Theorem U is proved in the note; the inputs (Q) are listed once.
6. Novelty of the formal stationary-phase lemmas (standard in spirit; proved in the appendix for completeness), of the `Ginf` formula
   (GJR18 Lemma 7.7), or of the chart mechanism (Givental s.10). For Theorem Gamma: "we are not aware of a proof in the literature".
7. Any use of topological recursion on the quintic's mirror curve (none is made); FCCZ26 is not used.
8. The gap or N(g) for frames outside the sources' conifold frames, for `g = 1`, or for `q` away from the conifold point.
9. Priority statements with dates, or a judgement of the cited preprints (CGYZ26 is a preprint, v2; Br25 is "to appear").
10. Do not cite arXiv:2610.38473 (withdrawn) or unpublished lecture notes without the authors' decision (`\draftnote`).

## 4. Conventions of the literature (sign dictionaries; flag every difference)
* **R-matrix action.** PPZ15 Def. 2.2: legs `R^-1(psi)`, edges `(eta^-1 - R^-1(psi') eta^-1 R^-1(psi'')^t)/(psi' + psi'')`, "the
  same convention is used by Givental and Teleman"; the project uses it (LIT s.1.2). arXiv:1712.02528 uses `R^-1` relative to it (LIT
  s.1.3, footnote 11 there; not re-checked by me). GJR18 Lemma 5.2: `S C_alpha^-1 = e^(u_alpha/z) R`.
* **Characters.** Giv01 Thm 1: tangent characters, sign minus; MR1901075 9.1: cotangent characters, sign plus; equal (odd powers).
* **Twist weight.** GJR18 (17): `(-5 h_alpha)` (authors' LaTeX); GJR17 (20): `(t - 5H)`; MR2276766 Cor. 1: an Euler-class twist enters
  as `exp(+Sum b_m (z/(lambda + rho))^(2m-1))`; LPfq (26): `(5 lambda_i)/e(E^* (x) (5 lambda_i))`. Project: `w_i` contains `-5 zeta^i`
  and `h_w = exp(-Sum_k b_k z^(2k-1) Sum_(a in w) a^(1-2k))` (MECHANISM_STAGE5 s.2). All agree (LIT s.1.4).
* **Bernoulli numbers.** MR1728879 and Br25: `t/(e^t - 1) = Sum B_n t^n/n!`; MR1901075 s.4: `x/(1 - e^(-x))`; Giv01 s.2.2: `z/(e^z - 1)`;
  MR2276766 s.3 prints its own expansion. Only `B_2k`, `k >= 1`, and the polynomials `B_k(h)` (DLMF 24) enter: no conflict.
* **Normalisation of asymptotic series.** ZZ08 Thm 4(i): `Phi_s(0) = delta_(0,s)` (the I-function's); LPfq (18) has the same; (QRR)
  is a different normalisation (by `h`). GJR18 Thm 7.1 normalises by `e^(u_alpha/z) (-2 pi z)^(m/2)`.
* **Landau-Ginzburg variables and measures.** GJR18 s.7: `(x_i)_alpha = L_alpha - lambda_i`, `y_alpha = m L_alpha`, measure `dx/x`
  with `x^(k/m)`; the project's charts: `x_j = zeta^j - L`, `y = 5L` (a sign flip), with `dy` for `Gamma(1 + 5/z)` and `dx_j/x_j` for the
  others (STAGE5_SYNTHESIS_2 s.7 item 1); MR1901075 s.10: `dx`; MR2553377 Lemma 4.20: `dw/w`, at `z = -1`; MR2276766 (11): `dx` with
  `(2 pi z (lambda + rho))^(-1/2)`.
* **Gap constants.** HKQ09 (1.2): `(-1)^(g-1) B_2g/(2g(2g-2)) t_D^(2-2g)`; CGYZ26 Conj. 1.5 (quintic): `(-1)^(g-1) B_2g/(5^(g-1)
  2g(2g-2)) tau_c^(2-2g)`; CGYZ26 Def. 1.8 (gap property of a CohFT) has `kappa^(1-g) (-1)^g B_2g/(2g(2g-2)) omega^(2-2g)` and Thm 1
  concludes with `(-1)^(g-1)`: do not import Def. 1.8's sign; Br25 (1.5) (local `P^2`): `3^(g-1) B_2g/(2g(2g-2)) t_CF^(2-2g)`. Project:
  `c_g = (-1)^(g-1) B_2g/(2g(2g-2)) 5^(1-g) = |B_2g|/(2g(2g-2) 5^(g-1))` (NATURAL_FORM s.1; STAGE2_SYNTHESIS s.1); model:
  `F_h^TR = B_2h/(2h(2h-2)) v^(2-2h)`, frozen `F_h = (-5)^(1-h) F_h^TR` (STAGE2_SYNTHESIS M8). Consistent.
* **Topological recursion.** Eynard:2007kz (4-3): `dE_q(p) = (1/2) Int_q^(q-bar) B`, opposite to earlier papers (its footnote); the
  project's kernel is minus the printed (4-3), and ABDKS25 (1.1) with `omega_(1,0) = -y dx` equals the printed one (STAGE2_SYNTHESIS_2
  K1); dictionary `omega_(h,n) = (-1)^n W^print_(h,n)`, `F_h = -F^print_h` (STAGE2_SYNTHESIS M9); Eynard:2007kz (4-10) defines
  `F^(g)`; the convention of Hock25 Cor. 1.2 must be stated (STAGE2_SYNTHESIS_2 O13). Project statements, not re-derived by me.
* **Hurwitz-space names.** The model `p^4 + a_2 p^2 + a_3 p + a_4 + a_5/(p - a_1)` (NATURAL_FORM s.2) is Dunin-Barkowski:2015caa (8-1)
  with `n = 4`, `d = (n-2)/n = 1/2`; Dubrovin-type name `M_(0;3,0)` (pole orders 4 and 1; NATURAL_FORM s.0); CvdLPS23 calls it
  `H_(0;1,N-1)` with `N = 5`. One space, three names: the note should fix one (NOTATION) and mention the others once.
* **Statement numbers.** Teleman: v3 PDF numbering, not ar5iv's (L24). All other numbers above agree between the HTML renderings and the
  PDFs (checked: GJR18, Dunin-Barkowski:2015caa, ABDKS25, Eynard:2007kz, PPZ15, Hock25, Giv01, MR1901075, MR2276766, MR1728879, ZZ08,
  CvdLPS23).
