# DIGEST_C1 (piece 6 of 6): s.6 Literature; s.7 Open issues for the writers

## s.6 Literature used (as the sources give it; root `refs.bib` key where one exists, searched by eprint number and title)

| arXiv id | statement numbers | what is used (digest item) | read in the full text by the sources? | root refs.bib key |
|---|---|---|---|---|
| 1303.1043 | s.2; Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9 (S2A s.2 list (Q)) | the R-matrix action: edge formula, legs, translation (D12, B1, C1) | statement numbers located in the PDF by the Stage-2 flat referee (S2A s.2); RD1 cites "s.2" only; TRF: "standard" | none found |
| 1509.06954v2 | Thm 8.1 (Hurwitz superpotential `p^n + a_2 p^(n-2) + .. + a_n + a_(n+1)/(p - a_1)`: CohFT = TR on `(CP^1, f, p, B)`); Thm 5.1 with conditions (1)-(4); Prop. 5.2 (homogeneity of R, TRRF corr. 4); (2-7) metric; (2-26)-(2-35) TR conventions and dictionary | TR-1, TR-2 (quote (Q)) | yes: TRM s.1 (s.1-8 of the PDF); TRRF s.1 (pp. 1-3, 6-13, 16-19, 24-27; "transcription accurate"; proof of Thm 8.1 uses `p_Dub = p - a_1/(n+1)`) | `Dunin-Barkowski:2015caa` |
| 2502.18115v2 | Cor. 1.2 with formula (1.3); Thm 1.1 = Thm 3.2; s.2.1 hypotheses; kernel (2.1) (p. 6); `F_g` (1.1) = (3.1); b-points (p. 13); s.3.2, Cor. 3.3, Prop. 3.4 (printed sign slip) | TR-4 route 1 (optional quote of (Q)) | yes: TRM (s.1-3.2, Prop. 3.4-3.5); TRRF (pp. 1-21); XYRF (pp. 1-5) | none found |
| 2212.00320 | Thm 7.1 (formula (7.9)); Thm 1.8, 1.14; (1.1) kernel; (1.13); example (7.2); Remark 7.3 | TR-7, TR-8, TR-9 (quote (Q)) | yes: XY (s.1, s.6-7); XYRF (s.1.1-1.5, s.7, pp. 2-12, 46-49) | none found |
| 2211.12259 | Thm 2 (H^2 tower = rooted hypermaps; calibration (1.24)); (1.6)-(1.8), (1.18); main theorem (extended rationally (N-1)-constrained KP) | MOD-3 (superseded); context only | yes: MS2 s.8 (pp. 1-11); MRF item 5 (Thm 2 and its proof, pp. 7-10); TRM (s.1-2.1) | none found |
| math/0108100 | the Fock-space formula `(R^)G(q) = [exp(hbar E_R(d,d)/2) G](R^-1 q)` (no number given) | B6 (not used) | not stated; "conventions not re-checked" (RD3 s.7) | `MR1901075` |
| math-ph/0702045 | kernel (4-3)/(4-9); Thm 4.1; Thm 5.1 (variational formula, s.5.3 (5-27)-(5-30)); s.4 (linear loop equation) | CV9 (orientation); TR-10, TR-11 (variational identities); the chain's uses (Theorem V, M6) belong to another reader | by the Stage-2 logleg referee (S2B K1), not by my sources' authors | `Eynard:2007kz` |
| 1211.4021 | (2-30)-(2-35) of 1509.06954 derive from it | none (TRM s.1 context) | no | `DuninBarkowski:2012bw` |
| 1110.2949 | "R-action = deformation of B and y dx", quoted "from memory, not checked" (XY s.4) | none | no | `Eynard:2011ga` |
| 2312.16950; 2408.02608; 2211.08917; 2304.03032 | log TR; generalised TR; x-y relation (genus 0; general genus) | none ("not needed", XY s.1; RD3 s.5.2; S2B s.6.4 lists 2312.16950 for an optional (Q1)-free proof of (V)) | abstracts or search summaries only | none found |
| 1411.7039 | two polarisations / Fock sheaf | RD3 s.5.3 reformulation only | no statement number | none found |
| hep-th/9506122 | c = 1 at the self-dual radius = the conifold | RD1 s.2.2 interpretation only | no | `Ghoshal:1995wm` |

Standard facts used without a pinned citation in my sources (each needs a citation with statement number or a written proof): the
uniqueness of the homogeneous R-matrix of a semisimple conformal Frobenius manifold (TR-2; TRM s.2.3, s.8); "the CohFTs over the
Frobenius manifold form the shifted family" (TR-5; TRM s.8); the dilaton equation and homogeneity of the model (D3, D4; LG s.3, UR s.1);
the difference equation of the Barnes `G` function (MOD-4); Lagrange inversion (TR-4); the GUE free energies (XYRF s.2). The project
cites arXiv:0712.0160v3 (root key `MR2917177`) Theorem (2) and Prop. 8.5(ii) only for the (QRR) uniqueness (MECH5 s.2), not for TR-2.

## s.7 Open issues for the writers (conflicts, gaps, \draftnote candidates, placement)

**O-1. Notation clashes (for NOTATION.md).** `c`: charges `c_a`; LD constants `c_b = (-1,1,-1,-3,-2)`; the scalar series `c(z)`;
vertex constants `c_g`; branch points `c_j`; `c0`. `e`: Euler-dilaton vector; basis vectors `e_a` (`Y^C(0) = -e_1`); the exponent `e`;
TR's normalised canonical basis `e_j`. `B`: `B = prod b_a`; the Bergman kernel; Bernoulli `B_2g`; black polygons (MOD-3). `D`: UR's
diagonal series `D`; D(ii), D(iii); the D-directions; `D = d/dx` (LD); TRM's class degree `D = 3g - 3 + n - Sum d`; univ's dressing `D`.
`K`: the slice; UR's frame `K = (1, uH, H^2, H^3, H^4)`; `K_q` (Lemma N); the Hankel `K_0` of Stage 5; the TR kernel `K(z0, z)`.
`L`: the leg class; the Mellin-Barnes variable `L` (Stages 4-5). `Lambda`: theory's slice `Lambda_v`; TR's `Lambda = -5 kappa`; the leg
map `Lambda(s)` (dict). `P`: the dressing; the polynomial part `P(z)` in TR-4; `P_(i,d)` in TR-3; Theorem P(b). `phi`: the basis `phi_a`;
B5's covector series; TR-4's series `phi(t)`; UNI-2's coordinate change; also `Phi` (primitive of `y dx`), `Phi_h` (D5), `Phi^SP`.
`S`: the calibration; `S(t) = sinh(t/2)/(t/2)`; weight spaces `S_w`; bvar systems S1-S7; `s` (B3) versus `s^5 = u` (TRM, univ). `tau`:
descendants versus `tau(u)`. `x`: the TR function versus the LG variable of "`x^4 + v/x`". `z`: descendant variable versus XY's curve
coordinate (S2A uses `p`). `u`: conifold coordinate versus TRM's residue in `W_u`. `H`: the class, `H(k,b)` (univ), `Hc`, `H_(0;1,4)`,
the purity hypothesis (H). `omega_(g,n)` versus `omega^`. "frame": CV16. Statement names: "Lemma 0" (frames; CYC constants), "Lemma T"
(leaf; Stage-5 translation invariance), "Lemma C" (Fock form B6; BVRF coordinate UNI-2), "Lemma F" / "Corollary F" / (F), (F_g),
"Lemma E" / "(E)" / E1-E13 (SYN01) / E2-E5 (univ) / E4 (the finite identity), "Lemma K", "Theorem G" (gauge identity; gap), "Theorem R"
(conifold), "Corollary U" / Theorem U, "Corollary L", ML1/ML2 (two meanings), M1-M10 (Stage 2) / (M1), (M2) (XY). Rename consistently.

**O-2. Exact hypotheses of quintic Theorem R and of the transport lemmas.** State them as in D11 and C1: flat-regular, Z/5-graded,
`P_0` invertible, (Sp), (A); the model at `v(u)` with `v(0) = 0`, `v'(0) != 0` (used in B3(d), and implicitly in C1(b): `T_P^(2-2g)` and the convergence of `Phi_g(v + s)` need `v` with a simple zero, RD2 s.3.2 step (3)); ML_K(h) for `h <= g - 1` (a) and `h = g`
(b). For `P_can`: (U), (CYC) (S5A S9), (G), Prop. S, Lemmas V, Q. Use the precise leaf statement of B3 (`q_T in K`, `Pi_0(q_T) = T_P`,
`d q_T in -z + K`, `Pi_0(d q_T) = v + s`) instead of the shorthand "`Pi_0 = t_c`".

**O-3. Steps needing expansion ([W-c] in s.3)**: none is a gap; all were checked by TRF/MRF, but rule 2 needs sentences:
(a) Lemma 0's existence and direct-sum step for infinite sums, with the ring stated (A1 steps 2, 4); (b) the duality step "annihilated
by all `V^(k)` in a slot => no `phi_1` component" (B2 step 5); (c) the exclusion of the other `u^1` terms in `s = O(u^2)` (B3 step 6); (d)
the `h = 1` case of the homogeneity identity in Lemma F (B4 step 1); (e) the genus bound for vertices of graphs with an edge (only in
TRF s.2; C1 step 2); (f) the u-adic convergence of frame expansions inside correlators and its agreement with flat expansions (C1 step 6;
MS2 s.10 "argued, not formalised"); (g) the resummation in the vertex graph (only in TRF s.2; C1 step 9). Suggestion: one lemma "u-adic
convergence" (from RD2 s.3.1 remark, s.3.2 (1)-(3), TRF s.2, MRF item 8(c)-(d), S2A M5) placed before Theorem R.

**O-4. Theorem R(c), Prop. I, the deformation route, the Fock-space form**: (c) is a one-line outline, Prop. I is SKETCH beyond the
model, neither is used. Omit, or a remark with a \draftnote ("outline only").

**O-5. Coefficient rings**: fix once (homogeneous elements over `Q[v]`; `Q[[u]]` after `v = v(u)`; completed tensor products for edges).

**O-6. Lemma F wording** (B4): "Laurent monomials in `u`" -> "monomials in `v(u)`".

**O-7. Lowered versus vector form** (CV6): RD1 does not distinguish them; the note must fix one and give `P^vec = J P J`, `d = (P_0)_33`.

**O-8. The ring of Lemma E's expansion** (B1 step 2): TRF says the expansion of `1/(x + psi')` is legitimate "in `Q((x^-1))[[psi1]]`";
the term `R^-1(x) x^(-j-1)` has infinitely many positive powers of `x`, and the residue pairing with `alpha(x)` involves sums that converge
only u-adically (RD2 s.3.1 remark). One precise sentence fixing the ring is needed. (A presentational point; the lemma is CONFIRMED.)

**O-9. Corollary L** (C3): the input `t_c - vhat = u^2/30 + O(u^3)` is VERIFIED only; keep C3 as a remark with the number cited from the
evidence appendix, or \draftnote "derive `[u^2] s` by hand". It is not used by the chain.

**O-10. Uniqueness of the homogeneous R-matrix** (TR-2): two sentences in TRM s.2.3, labelled "standard, not re-proved" (TRM s.8),
inputs to be named (TRRF correction 4). Needed for "frozen = TR" to all orders (M1). Write the recursion argument with the named inputs,
or cite with a statement number; otherwise \draftnote.

**O-11. TR-5's "shifted family" input** (standard, not re-proved): not needed by the chain (ML_K contains the primary-sector lemma).
Keep TR-5 as a remark or omit.

**O-12. TR-8, the Borel step**: one paragraph in XY s.2.3. Route 2 is optional (it removes the need for arXiv:2502.18115); if kept,
write the Borel identity as a formal identity of coefficient sequences.

**O-13. The genus-1 base** (TR-11): XY's `F_1 = -(1/12) log v` is WEAKENED (orientation fixed by consistency). The chain's `h = 1` base
(S2A M9) uses homogeneity and the value `-1/12`, "one exact `Q(zeta5)` number"; rule 3 forbids a proof citing a computation. A written
derivation is needed (sources containing explicit data: XY s.2.2's closed `omega_(1,1)`; XY s.2.3), else a \draftnote. Owner: the
reader of the Theorem V part; flagged here because the XY derivation is one of my sources.

**O-14. Conjecture A** (UNI-3): only its flat form follows (Theorem R for non-gauge `P`, plus ML_K). The hatted form for non-gauge
dressings has no written S15-type identification: state universality in flat form only; the hatted engine form of UNI-1 is SKETCH
(S2A s.2).

**O-15. Necessity** (UNI-8, MOD-1, MOD-2): examples only, in the evidence appendix; no theorem "each hypothesis is necessary"; (Sp) has
no example; "regularity <=> annihilation at first order (genus 2)" is a finite computation. Phrase as "cannot be dropped in general,
see Appendix".

**O-16. Graph-locality** (UNI-7): a methodological remark (introduction or s.9), numbers in the appendix; it explains why the proof is
graph by graph and why no total-potential argument can work.

**O-17. Material to exclude.** TRM s.6 and XY s.4 (Delta B deformation of the model curve by the gauge, contour arguments) are SKETCH and
unchecked; XY s.5 (matrix model) SKETCH; RD3 s.5-s.6 (candidate routes, ranked Stage-2 routes) and RD2 s.3.2's comparison with the
conifold project are project planning (rule 5). TR is used only on the model curve (rule 8); none of these sketches may be presented
as an argument.

**O-18. Superseded wordings to avoid**: ML "CONJECTURAL" (RD1, MS1, MS2, SYN01 E10); Cor. N, N_9 (RD2 s.3.3); "1013 singular lists"
(RD1 s.0, s.2.4: 326 genuine cancellations); "flatness is neither sufficient nor necessary" (UN1 s.0 item 4); "genus 4 is not a sharp
test" (MS1 s.0 item 8, MS2 s.9); "`Y^C` and the TR monodromy are different functionals" (XY s.3; REFUTED, SYN01 C5); Lemma K's
"`phi_j -> 0`" (MS2 s.6); "D(ii) PROVED (critic)" (MECH1 s.1; sketch-level, proof of record dict referee s.3, S2B O26); "ML on `K_0`
IS L-regularity" (SYN01 s.0; also needs LD, D(iii), M5, S2B O27); "each hypothesis on `P` is needed in examples" (RD1 s.0, RD2 s.3.4).

**O-19. Hard-rule items in the sources**: author names in TRM, XY and MS2 s.8, and author-named citations (TRRF s.6 item 1; XYRF s.5;
MRF s.4): cite by key only. TRM's "EO/DNOPS conventions" label for the Gap Theorem is inaccurate as to the printed EO kernel (CV9).

**O-20. Placement (suggested; numbering of S2B s.6.3).**
* s.3 (tower, slice, frames; transport): D1-D10, D12, D18; A1 (+ A3 inside it or right after), A2, B1, B2, B3, B4, the convergence lemma
  of O-3, B7 for the quintic; remarks A5, A6; B5 as a remark (second proof) optional; B6 omitted or a remark.
* s.4 (quintic Theorem R): C1 (a), (b); C2 remark; C3 remark (O-9); C4, C5 omitted or one remark.
* s.5 (the model as TR): TR-1, TR-2 (O-10), TR-3 (or LD), TR-7 (hypothesis check of the quote); MOD-3 optional remark.
* s.6 (L-regularity, Theorem V, Gap; ML_K): TR-9, TR-4 (route 1, with TR-8 as an alternative proof or remark); MOD-1 and MOD-5 as
  corollaries of ML_K; MOD-2, MOD-4, MOD-6 remarks; TR-5, TR-6, TR-10-TR-12 remarks or omitted.
* s.9 (universality): UNI-1, UNI-2 (renamed), UNI-3 (flat form), UNI-4, UNI-7, UNI-8 remarks; UNI-5, UNI-6, UNI-9 in the appendix.
* Appendix A (sign conventions): CV1-CV17. Appendix B (computations): s.5. Appendix C (errata), if kept: the corrections recorded in s.2
  (as mathematical corrections only, without project bookkeeping).

**O-21. Inputs not digested here** (other readers): (G), Prop. S, Lemmas V, Q, 1 (GI s.2-s.4); Theorem U and (R1)-(R3) (UR); the
definition of Conjecture N and of the hatted outer sum (NF); Theorem S15 (flat worker); LD, D(ii), D(iii), Theorem V (Stage 2);
(CYC), Theorem P(b), Theorem Gamma, the glue (Stages 4-5). Every statement of this digest that uses them says so.

**O-22. Conflicts recorded in this digest, with their adjudications** (index): V^(k) versus Pi_i (CV1, no source conflict, a derived
sign); weight sign conventions (CV3, equal); `c_g` forms (CV4, SYN01 C2); `<H>_1`, `<e>_1` signs (CV5, SYN01 C3); lowered versus vector
(CV6, TRF s.2); `m5` versus `c_5` (CV7, SYN01 C11); curve normalisations (CV8); kernel orientation (CV9, S2B K1); leg dictionaries (CV10,
S2B K10) and the string sign (SYN01 C8); the Euler-dilaton leg (CV11, SYN01 C7); metric and psi signs (CV12, TRRF corrections 2, 6);
MRF's `Pi_0` sign (CV13); three leaf differences (CV14); `lambda`, `c0` (CV15, SYN01 C10, S2B K6); "frame" (CV16); genus-4 sharpness
(MOD-6, SYN01 C1); Lemma K's hypothesis (B5, MRF item 6); "no single V" scope (MOD-2, MRF item 3); `Y^C` versus the TR monodromy
(TR-10, SYN01 C5); the s.1.5 input (C6, S2B K8); D(ii)'s status (O-18, S2B O26).
