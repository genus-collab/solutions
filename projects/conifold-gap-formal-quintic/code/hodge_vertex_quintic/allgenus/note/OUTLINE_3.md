# Outline, part 3: s01 introduction; s02 set-up, items 1-9 (items 10-17 in `OUTLINE_10.md`)

Symbols as in `NOTATION.md` (renamings: charges `nu_a` (source `c_a`), cyclic coefficients `kappa_a`, `b_a` (source `c_a`, `b_a`),
vanishing period `varpi_1` (source `g`), free coefficients `frak f_m` (source `f_m`), normalisation series sans `D` (source `D`),
Theorem-U freedom sans `c` (source `c(z)`), LD constants `varrho_b` (source `c_b`), insertions `psi^k H^b` (source `tau_k(H^b)`),
vertex genera `g_v`, the PPZ translation `\ttr_R` (source `T`), the degree variable `\Qn` (source `Q`)).

## s01_intro (written last, by the assembler) — target 5 pp

Purpose: state Theorems A, C, Gamma and Corollary B with their exact hypotheses, list every input once with the conditionality
sentence, position the results against the literature without overclaiming, and give the strategy and the plan.

Covers: DL s.2.1-s.2.3 (novelty, positioning drafts, must-not-claim list) with REV R4, R5; DL s.1 L1-L19, L43-L48; MS5 s.1-s.2;
NF s.0-s.2; DA A14-A20, A24; DC1 UNI-7; OUTLINE s.0.2 (inputs, conditionality); G5_2 s.2.6 (the constant).

Contents, in order.
1. `sec:int-setting`. The formal quintic [LP18, s.0.2 (5); LP20, s.0.2 (5)] as the O(5)-twisted equivariant theory of `P^4` at
   `lambda_i = zeta^i` [ZZ08]; its semisimple CohFT, R-matrix normalised at `q = 0` by quantum Riemann-Roch, Givental-Teleman graph sum
   [MR2917177; PPZ15; GJR18, Thm 4.1]; the conifold point `u = 1 - 5^5 q = 0`, the vanishing period and `t_c`. One paragraph: the
   computations at genus 2 of appendix C show that the linear Hodge frames of the conifold case fail here (DA A18, FQ2 s.3; REV_2 M6;
   numbers only in appendix C); the analogy of the natural form with the Hodge factor (DA A20; NF s.2; UR s.7; `rmk:mai-analogy`).
2. `sec:int-gap`. Positioning paragraphs 1-2 of DL s.2.2, transcribed with REV R5: physics origin [Bershadsky:1993cx via HKQ09;
   Ghoshal:1995wm; HK07; HKQ09 (1.2)]; the quintic [GJR18 abstract, Rem. 1.6; CGL21; CGL26; LR22]; proved cases [Br25; CI21 as an
   earlier check; FCCZ26 for TR of toric mirror curves; CGYZ26, Theorem 2, for the quintic, and nothing else under that number]; the
   gap of the formal quintic open [CGYZ26, Rem. 1.9]. The other targets of CGYZ26 s.5 are not cited (not read). DL s.2.2 paragraph 4
   with REV R4: "the analogous identity for their Landau-Ginzburg integral is stated without proof in [GJR18, Remark 7.2]".
   `\draftnote{re-check the latest versions of arXiv:2609.13921 and arXiv:2509.19298 before submission}`.
3. `sec:int-results`. Statements, each with the hypotheses of its body version:
   * `thm:int-main` (Theorem A) = `thm:mai-main` (text of OUTLINE s.0.2, with `T_Gamma` as in REV R6), "for every `g >= 2` and every
     conifold frame of `def:set-frames`".
   * `cor:int-gap` (Corollary B) = `cor:mai-gap`: in every such frame the conifold expansion `F_g^cf` of the genus-`g` potential
     (`def:mai-Fcf`) is `c_g t_c^(2-2g) + Q[[u]]` (no fractional power of `u`, no `t_c^k`, `3 - 2g <= k <= -1`). No `\draftnote` (G5_2 s.2).
   * `thm:int-D` (Theorem C) = `thm:gam-D` with `cor:mai-actual`: `D = 1`; the actual R-matrix is in the canonical normalisation.
   * `thm:int-gamma` (Theorem Gamma) = `thm:gam-main` (ii) and (iv): `c(z; 1) = h_(w_0)(z)`; `c(z; 0) = Ginf(z)`, the series of [GJR18,
     Lemma 7.7, (37)] (REV_2 M4: the note's content at `L = 0` is the value of its own series, not a new formula).
   * One sentence: Theorem A is graph-local; the computations of appendix C show that a statement about the total cannot imply it
     (DC1 UNI-7; REV_2 M6).
   * `rmk:int-constant` (D10): `c_g` equals the constant of [CGYZ26, Conj. 1.5] in the same normalisation of the vanishing period, and
     the two conifold expansions agree (G5_2 s.2.6; pointer of G5R F2.4: [CGYZ26] Ex. 1.3 and s.4.2, Prop.-Def. 4.5 (1));
     `\draftnote{numbers from arXiv:2609.13921v2, read in the HTML rendering; to be re-read in the PDF}`.
4. `sec:int-inputs`. The list of OUTLINE s.0.2 ("Exact final list of inputs"), verbatim, with which steps use which (column "inputs"
   of the chain table), and the conditionality sentence of OUTLINE s.0.2. State: Theorem U and every lemma of the chain are proved in
   the note; topological recursion is used only on `x = W_v(p)`, `y = p`; no finite computation enters a proof. There is no list of
   open chain steps (none is open); the derivation of (QRR) from the CohFT is a remark (`rmk:set-qrrderiv`).
5. `sec:int-notclaimed`: DL s.2.3 items 1-10, as sentences: nothing about the true quintic (the formal quintic is a different theory:
   `F_2(0) = -19/9600` against `-5/144`, DA A24, numbers in appendix C); no claim about GJR's second key statement; nothing for frames
   outside `def:set-frames`, for `g = 1`, away from `u = 0`; (QRR) is quoted, not proved; no novelty claimed for the formal
   stationary-phase lemmas, the `Ginf` formula [GJR18, Lemma 7.7, (37)] or the change of chart [MR1901075, s.10]; for Theorem Gamma "we
   are not aware of a proof in the literature", and the analogous identity for GJR's Landau-Ginzburg integral is stated without proof
   in [GJR18, Remark 7.2] (REV R4).
6. `sec:int-strategy`: the chain C1-C13 of OUTLINE s.0.2 in prose (half a page), and where (QRR), (M), (G0), (T) and each quote enter.
7. `sec:int-plan`: plan of the note.

`\draftnote`s: title (working; in `main.tex`); D11 (`\draftnote{the authors decide whether to cite their companion notes and how to
relate the note to the sentence of [Br25, s.1.1] on R-matrix methods}`); Teleman v3 and Givental v2 numbering at the first citations
(s.0 rule 4). Dropped: every status, every date, the human-check lists. Length 5 pp.

## s02_setup — target 10 pp

Purpose: define the formal quintic with its genus-0 data and the quoted inputs (M), (G0); prove the large-radius cyclic structure of
the rows; define the rows with (QRR) and the graph sum with (T); the conifold frames and the conifold expansion; the model, `M_0` and
its symplecticity, the model CohFT, its homogeneity and the dilaton lemma; the hatted theory, the outer graph sum and Property N(g).
Proofs written here: `prop:set-cyclicLR` (about 1.5 pp), `lem:set-frobbasis`, `lem:set-M0sympl`, `lem:set-modhom`, `lem:set-dilaton`;
the others at most a paragraph.

Covers: DA s.1 D1-D20, A1-A6, A16; DB D1-D4, D6-D8, D16, S1, S2, S18 (a)-(b); DD1 D1.1-D1.3, D1.9-D1.11; DD2 D3, QRR-1..3; DC1 D1-D3;
DC2 D6, D14; FQ2 s.1-s.2; OU2 s.1; S15 s.1.3; UR s.1-s.2; LG s.0-s.3; G1 (all pieces), G1R, G1R_2; G4_2 s.5; G5_2 s.2.2; GS s.2.1-s.2.4,
s.2.10, s.2.11.

Subsections (REV_2 M14): `sec:set-fq` items 1-6; `sec:set-rows` items 7-8; `sec:set-graph` item 9; `sec:set-frames` items 10-12;
`sec:set-model` items 13-15; `sec:set-hat` item 16; `sec:set-N` item 17.

1. `def:set-fq` (defn). Torus weights `lambda_0..lambda_4` on `H^0(P^4, O(1))`; `N_(g,d) = Int_[Mbar_g(P^4,d)]^vir e(R pi_* f^* O(5))`,
   the Euler class defined after virtual localisation [MR1666787], then `lambda_i = zeta^i lambda`, `lambda = 1`; fixed points `p_i`,
   tangent weights `lambda_i - lambda_j`, fibre weight `5 lambda_i` ([LP20] (10)); `F_g := Sum_(d >= 0) N_(g,d) \Qn^d` (`d = 0` included,
   FQ2 s.2), with the degree variable `\Qn` (REV R8 (g)); along the mirror map `\Qn = e^(T(q))` (point `t = T H` at Novikov variable 1,
   G1 s.1.1), so `F_g` is a function of `q` (G5_2 s.2.2). Sources DA D1; LP18 s.0.2 (5); LP20 (5), (10); ZZ08 (GJR18 s.3.1 (9)).
2. `def:set-genus0` (defn). `D = q d/dq`; `C_0 = Sum (5d)!/(d!)^5 q^d`, `I_1 = C_0 log q + I~_1`, `T = I_1/C_0`, `C_1 = D T`; `L = (1 -
   5^5 q)^(-1/5)`, `D L = (L^6 - L)/5`; the I-function at `p_i` (DA D2; G1 s.1.1); PF5 `= (zD)^5 - 1 - q prod_(k=1..5)(5 zD + kz)` and PF4
   `= D^4 - 5q prod_(k=1..4)(5D + k)` (DB D2), their `u`-forms. The eigenvalues of `H*`, `Delta_i` and `D u_i` are consequences, stated in
   `cor:set-cyclicPF5` (a): the eigenvalues of `H*` are `zeta^i L/C_1` (the earlier outline wrote `zeta^i L`, the eigenvalue of `C_1 H*`).
   * `lem:set-pf5` (lem) `PF5 = z^5 D PF4 - 1`. Proof DB S1 steps 1-3 [W] (CYC Cor.; CF2 s.8.2; `D q = q (D + 1)`, G1_2 s.2.1).
3. `thm:set-mirror` (quoted, (M); DEC D2, Q7). Statement G1 s.1.1: twisted genus-0 invariants with `e_T(E_(0,n,d))`, twisted pairing,
   the J-function in the normalisation of [MR2276766, s.6 (10)] with the twisted pairing; at `lambda_i = zeta^i`, Novikov variable 1,
   `t = T H`: `J(T(q) H, z) = I(q, z)/C_0(q)`. Pin: [MR1408320, Theorem 11.6], Calabi-Yau case `n = 4`, `r = 1`, `l_1 = 5`, `lambda' = 0`,
   rewritten for `s_i` in two lines (GS s.2.1), with the set-up of s.9 (p.23-27) and Prop. 11.4; one paragraph with the translation of
   Givental's `s_i` to the twisted J-function (G1 s.1.3 items 1-5) and the hypotheses (G1 s.1.4: genus-0 convexity of `O(5)`,
   polynomiality, `e_1(zeta^i) = 0`); one clause: [MR1408320, Cor. 11.7] prints `l_alpha` for `lambda_alpha` (G1R s.1.1). Givental's
   letters (`s_i`, `f`, `E_i`, `g_1`, `hbar`, `l_a`) are declared local. `\draftnote{numbers from arXiv v2; IMRN numbering not compared}`
   at the first citation. Optional sentence: [MR2276766, Thm 2] is the modern formulation (non-equivariant in `X`; not the pin, G1 s.1.5).
   A second source (LLY97) only under DEC Q7. The Yukawa coupling reproducing the genus-0 numbers is evidence (appendix C).
   * `thm:set-G0` (quoted, (G0)): G1_2 s.2.1 (G0-a)-(G0-d): the genus-0 twisted classes form a CohFT with unit for `eta` ([GJR18] s.4.1;
     [PPZ15] Def. 0.3); polynomial and `S_5`-symmetric in `lambda`, homogeneous of degree `Sum a_j - n`; string and divisor equations;
     splitting (hence WDVV) and the genus-0 TRR; the quantum differential equation ([GJR18] s.5.2 eq. (14)). Stated once, next to (M).
4. `prop:set-cyclicLR` (prop; REV R1). Statement G1_2 s.2.0: assume (M), (G0); at `t = T(q) H`, in the flat basis, `H * H^a = kappa^LR_a
   H^(a+1)` (`a` in Z/5, `H^5 = 1`), `kappa^LR = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1)`; equivalently for every solution `v` of `z \partial_T v =
   H * v` the lowered components satisfy `z D v_a = C_1 kappa^LR_a v_(a+1)`; for the colour-`i` column of `def:set-rows`, `X_(a+1) = (C_1
   kappa^LR_a)^-1 (z D + zeta^i L) X_a` (`a = 0..3`), `X_a := zeta^-i Q^(i)_a`, `X_0 = L Phi_i/C_0`. The superscript LR separates
   `kappa^LR` from the frame coefficients `kappa_a` of `prop:gau-quintic`. Proof [G]: local lemmas `lem:set-grading` (G1_2 s.2.2, Lemma
   A), `lem:set-unitfrob` (s.2.3, B), `lem:set-birkhoff` (s.2.4, C), `lem:set-periods` (s.2.5, D), `lem:set-wronskian` (s.2.6, E: the
   Wronskian is `L^10`), `lem:set-order3` (s.2.7, F); the proof G1_3 s.2.8, with the sentence "only derivatives along the line `t = T H`
   are used" (G1R s.2.1); G1_3 s.2.9 says what is used. Local letters (G1R s.2.3): `\partial_T = C_1^-1 D` (not s03's `\partial_t`), the
   formal `log q` is `\logq`, the Wronskian `\Wr`; declared local letters for G1's `s_a`, `g`, `w`, `c_k`, `phi_k`, `Mop` (with one
   sentence: the large-radius analogue of `M_cyc`), `X_a`, `Fhat`, `ctilde_m`. Optional sentence (G1_3 s.2.11): the product is used in
   [LP20] (17), (19)-(20) and [GR19] Props. 6.1-6.2; [GJR18] (13) is not cited (misprinted, G1R s.2.3).
   * `cor:set-cyclicPF5` (cor): G1_3 s.2.10 (a)-(b) [G]: (a) the eigenvalues `mu_i = zeta^i L/C_1` of `H*`, the idempotents `e_i = (1/5)
     Sum_a (prod_(b<a) kappa^LR_b) mu_i^-a H^a`, `Delta_i = C_0^2 zeta^(3i)/L^2`, `\partial_T u_i = mu_i`, i.e. `D u_i = zeta^i L` (constants as
     in `lem:set-canon`); (b) for every solution `v` of the quantum differential equation `C_0 v_0` solves PF5; for `v = R e~_i e^(u_i/z)`
     the function `e^(u_i/z) z L Phi_i(zhat, u)` solves PF5 (the WKB form of colour `i`). Part (c) of G1_3 s.2.10 is `cor:gau-rowsNN` (s03).
5. `def:set-points` (defn): `u = 1 - 5^5 q`, `s = u^(1/5)`, `L = 1/s`, `D = -(1-u) d/du`, `zhat = z/s^4 = L^4 z`; large radius `u = 1`
   (`L = s = 1`, `zhat = z`), conifold `u = 0`, orbifold `q = oo` (`L = 0`). DA D3; DD1 D1.1.
6. `def:set-frob` (defn): the state space `\bbH` with flat basis `1, H, .., H^4` (REV_2 M13 (e)); `eta(H^a, H^b) = 5 delta_(a+b = 3 mod 5)`,
   `eta^-1 = eta/25`; the idempotents of `cor:set-cyclicPF5` (a); normalised canonical basis `e~_i = Delta_i^(1/2) e_i`; square roots
   `zeta^(1/2) := zeta^3`, `Delta_i^(1/2) = zeta^(4i) C_0/L`; `W_i^0 = zeta^(3i)`; the TFT `\TFT` of the formal quintic,
   `\TFT_(g,n)(e~_i, .., e~_i) = Delta_i^((2g-2+n)/2)`, zero on mixed colours (REV R6). DA D4, C2, C3; G1_4 s.3.2 Step 2.
7. `def:set-rows` (defn): Givental's convention `S ~ R e^(u/z)`, `R = 1 + O(z)`, symplectic; the columns `R(z) e~_i e^(u_i/z)` solve
   `z \partial_T v = H * v`, coefficients power series in `q` along `t = T(q) H` (G1_4 s.3; [GJR18] s.4.1); rows in lowered form `Q_a^(i)(z) =
   eta(H^a, R(z) e~_i)` (S15 s.1.3; OU2 s.1). DEFINITIONS (G1_4 s.3.1): `Phi_i(zhat, u) := C_0 Q^(i)_0(z)/(zeta^i L)`, `zhat = L^4 z`;
   `ph_n` by `Phi_0 = Sum_n zhat^n ph_n(u)`, `ph_0 = 1`; flat amplitude `A_n := [z^n](C_0 Q^(0)_0) = L^(4n+1) ph_n(L^-5)`. The WKB form of
   colour `i` is `cor:set-cyclicPF5` (b); the covariance clauses `Q^(i)_(a,n) = zeta^(i(a+1-n)) Q^(0)_(a,n)` and `Phi_i(zhat, u) =
   Phi_0(zeta^-i zhat, u)` are `lem:gau-colrows` (s03; GS s.1 item (a)); HPF and its free coefficient are `lem:gau-HPF` (s03).
8. `thm:set-qrr` (quoted input (QRR); DD2 QRR-1): in Givental's convention `R(q = 0) = diag(h_(w_0), .., h_(w_4))` in the basis
   `e~_i(q = 0)`, `w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`, `h_w(z) = exp(-Sum_k bk_k z^(2k-1) Sum_(a in w) a^(1-2k))`. Form used,
   for every colour: `Phi_i(z, 1) = h_(w_i)(z)`; in particular `A^QRR_n(q = 0) = h_n`, `h = h_(w_0)`. Proof of the form used: G1_4 s.3.1
   [G] (`q = 0` is `u = 1`, `L = 1`, `zhat = z`, `C_0 = 1`; `eta(1, e~_i) = zeta^i L/C_0` by `cor:set-cyclicPF5` (a) and the square roots
   of item 6). Sources as in OUTLINE s.0.2. No `\draftnote` (the earlier (a), (b) are deleted).
   * `rmk:set-qrrsign` (rmk; DD2 QRR-2): the twist weight is `-5 zeta^i` (Euler-class twist: minus its fibre weight); agreement of
     [GJR18 (17)], [MR2276766 Cor. 1], [LP20 (26)]; `h_1 = -3/20` (`Sum_(a in w_0) a^-1 = 9/5`), whereas the weight `+5` would give `-11/60`;
     that `ph_1 = -3/20` is forced by HPF is `lem:gau-HPF` (forward reference).
   * `rmk:set-qrrderiv` (rmk; DEC Q1; replaces `rmk:set-qrrvertex`): how (QRR) follows from the CohFT; not used by the chain. Steps 1-5
     of G1_4 s.3.2 [G] with C3b-1 (the R-matrix of [GJR18] Thm 4.1 solves the QDE: [MR2917177] 8.7 Proposition (ii), `R = E exp h`, `E`
     solving 8.2 Proposition (ii); G1R_2 s.4.2), C3b-3 (Teleman's `T_z^-1 exp h(z) T_z (I_A)` is `exp(h).TFT` in PPZ's convention: [PPZ15]
     p.20, with the genus-1 sign check `r_1 = a_1/12`; G1R_2 s.4.4), C3b-4 (Step 5 rests on the non-vanishing of `ch_(2k-1)(E)`; neither
     Teleman 8.10 nor Madsen-Weiss is cited for it; G1R_2 s.4.5). Citations: [MR2917177] Theorem (2) in (1.15), (8.6), 8.7 Proposition
     (i)-(ii), 8.10 Proposition (context); [LP20] (26); [GJR18] Thm 4.1; [PPZ15] Prop. 2.14; Mumford's relation as in [MR1728879] s.2.1.
     One `\draftnote` (O1): `\draftnote{Mumford's relation ch_{2k}(E) = 0 and, for every k >= 1, the non-vanishing of ch_{2k-1}(E) in
     some H^*(Mbar_{g,n}) are used here without a pinned statement}`. Never "Prop. 8.5(ii) (s.8.2)". GS s.2.2's route: `DRAFT_NOTES.md` only.
   * `lem:set-hrot` (lem): `h_(w_i)(z) = h_(w_0)(zeta^-i z)` (from `w_i = zeta^i w_0`; DD1 cor:glu-colours; CP s.0.2). One line [W].
9. `def:set-graphsum` (defn): leg-free stable graphs of genus `g`, vertex genera `g_v`; the R-matrix action of [PPZ15] Def. 2.2 (vertices,
   legs `R^-1(psi)`, edges `(eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi')`, a power series by symplecticity, [PPZ15] Rem. 2.1 and
   (6)), with the translation of Def. 2.5/2.13, `\ttr_R(z) = z 1 - z R^-1(z) 1`; the contraction `Cont_Gamma[E, q; Omega]` INCLUDES the
   factor `|Aut Gamma|^-1` (S15 s.1.1; G5R F2.2), so sums over graphs carry no further factor; coloured graphs; the colour-basis form
   follows from the naturality `A_*(R.Omega) = (A R A^-1).(A_* Omega)` for a constant isometry `A`, immediate from [PPZ15] Defs. 2.2, 2.13
   (G1R_2 s.7; GS s.1 item (b)); the explicit colour-basis edge `E_ij` and the leaves `Dil_i(d)` only in appendix B (REV_2 M8).
   * `thm:set-GT` (quoted, (T)): the CohFT of the formal quintic is the graph sum with the R-matrix of `thm:set-qrr` ([MR2917177] Theorem
     (2) in (1.15); [GJR18] s.4.1, Thm 4.1). Calibration numbers (DA V1) in appendix C.
