# DIGEST_C1 (piece 2 of 6): s.2 Statements, parts A-C (tower and slice, transport, quintic Theorem R)

Continues `DIGEST_C1.md` (abbreviations and definitions D1-D19 there). Each item: proposed label (prefix by topic: `tow` tower and
slice, `tra` transport, `red` reduction = quintic Theorem R; suggested section in brackets, numbering of the write-up plan S2B s.6.3),
type, statement in final corrected form with hypotheses and scope, source, status in the sources, corrections applied, dependencies.
"Used by the chain" refers to the final chain S5A S8 = S2A s.2 steps 1-8.

## A. The vanishing tower, the slice, frames [s.3]

**A1 `lem:tow-frames` (lem; "Lemma 0").** `Hc_+ = K_0 (+) Q[[z]] phi_1`. For every `(k, a)` with `a != 1` there is a unique
`frame(k, a) = z^k phi_a + Sum_(j<k) f_j z^j phi_1` in `K_0`, the `f_j` homogeneous monomials in `v`; the frame coordinates of
`n in K_0` are its flat `phi_a`-components, `a != 1`. There is no frame with `a = 1`: the descendants `tau_k(H)`, `k >= 1`, are
transverse to `K_0`. Scope: coefficients `Q[v]` (homogeneous elements) or `Q[[u]]` with `v = v(u)` (u-adic sums, A3).
Source RD1 s.1.3. Status PROVED (RD1); frames rebuilt independently at general `v` with homogeneity asserted (TRF s.3, s.4). Used by
the chain: S2A step 4 and M5. Correction: rename (clash with "Lemma 0" of `stage5/critic/CYC_PROOF.md`); state the coefficient ring.
Depends on: D6 (`Y^C_1 = -1 + O(z^-5)`, `Y^C_a` has no `z^0` term for `a != 1`).

**A2 `cor:tow-slice` (cor).** `-z + K = -z phi_0 + (K_0 (+) Q phi_1)`, and `Pi_0(-z phi_0 + n + s phi_1) = v + s` for `n in K_0`: on the
slice the only coordinate seen by `Pi_0` is the residue direction `H`. Source RD1 s.1.3 ("The slice"). PROVED from A1 and the values
`Pi_0(-z phi_0) = v`, `Pi_i(-z phi_0) = 0` (`i >= 1`), `Pi_0(phi_1) = 1`, `Pi_i(phi_1) = 0` (`i >= 1`) (RD1 s.1.3, "From the shape of
Y^C"). Used by the chain (steps 4, 5).

**A3 `lem:tow-coeff` (lem; S2A "M5").** Over `Q[[u]]` (`v = v(u)`), `K_0` is spanned u-adically by the frames; in `frame(k, a)` the
coefficient `f_j` is `const . v^((4(k-j) + c_a + 1)/5)` with exponent `>= 1`. Hence every insertion of quintic Theorem R is a
`Q[[u]]`-combination of homogeneous elements of `K_0` and of `tau_0(H)`, and frame expansions converge u-adically and agree with flat
expansions inside any correlator (only finitely many flat slots survive the dimension bound). Source S2A s.1 M5; RD2 s.3.2 step (1)
parenthetical; TRF s.2 (Theorem R paragraph). Status PROVED (Lemma 0 plus the exponent count); it closes the dict referee's
correction 5(iii) (S2A M5). Used by the chain (steps 4, 7.4).

**A4 `lem:tow-Kbasis` (lem; model description).** (i) `K = Q tau_0(H) + span{kappa(X, k)}`; the `kappa(X, k)` are log-free with
non-negative powers of `v` and span `K_0`. (ii) `Hc_+ = K (+) span{tau_k(H) : k >= 1}` (direct: `Pi_i(tau_k(H)) = PiC(tau_(k-i)(H))` is
triangular with diagonal 1). (iii) The membership list of D10. Source MS1 s.1 ("K and its basis"). Status PROVED ("elementary");
MRF item 2 (`rf_flat`: `Pi_i(kappa) = 0`, list confirmed), item 5 ("no resonance" -> "log-free"). Not needed by the chain (which uses
frames, A1, A3); optional equivalent basis, used in the ML_K evidence (MS1 s.5).

**A5 `rmk:tow-consistency` (rmk).** (a) String: `l_(-1) q := [q/z]_+` satisfies `l_(-1) V^(k) = -V^(k+1)` (equivalently
`Pi_0([q/z]_+) = Pi_1(q)`), so it preserves `-z + K` and `Pi_0` on it; off `K` no function of `Pi_0` alone can be the singular part:
`<tau_1(H) tau_0(1)>_h = <H>_h != 0` although both insertions are killed by `Pi_0`, and `tau_1(H)` is not in `K` (`V^(2)(tau_1 H) = -1`).
Inside `K_0` the frame `frame(2,0) = tau_2(1) + tau_1(H)` compensates: `<frame(2,0) tau_0(1)>_h = <frame(1,0)>_h = 0` (RD3 s.4.2).
(b) Dilaton: `-z + K` is a cone and `Pi_0` linear, so `Phi_h o Pi_0` has the dilaton weight of `Fb_h`. (c) Euler-dilaton: given
ML_K, `<frame(1,0)>_h = (2h-2) F_h + v F_h' = 0` (`h >= 2`), and at `h = 1` the exponent is 0 and the anomaly
`<frame(1,0)>_1 = 5/24 - 1/12 = 1/8` is allowed (regular); `<frame(1,0) H^m>_1 = 0` for `m >= 1` (TRF s.2). (d) `phi_0, phi_2, phi_3,
phi_4 in K_0`, so ML_K contains the primary-sector lemma (MOD-5, piece 3). (e) On `-z + K` the singular part `F_h(Pi_0)` is also invariant
under every `q -> [z^-m q]_+` (`Pi_m = 0` on K), which are not symmetries of the full potential for `m >= 2`; higher Virasoro operators
not examined. Source RD1 s.2.1; RD3 s.4.2; MS2 s.9(d); TRF s.2 (CONFIRMED). Status PROVED (consistency checks). Not used by the chain.

**A6 `lem:tow-isotropic` (lem or rmk).** `Sigma_C = Q((z^-1)) J^C` is isotropic: `Y^C(-z)^T eta^-1 Y^C(z)` is constant in `u` (flatness;
`eta^-1` is proportional to `eta`) and vanishes at `u = 0` because `Y^C(0) = -e_1` and `(eta^-1)_11 = 0`. Source RD1 s.2.2; TRF s.2
(CONFIRMED). Status PROVED. Not used; interpretation (c = 1 sector, cf. hep-th/9506122) only.

## B. Transport [s.3]

**B1 `lem:tra-edge` (lem; "Lemma E").** Let `R(z) = 1 + O(z)` satisfy `R(-z)^T eta R(z) = eta`, `E_R` its edge (D12), and `f in Hc_-`
(vector form). Contracting the covector `Omega(f, .)` into the first slot of `E_R` gives the vector
`(Omega(f, .) (x) 1) E_R (psi') = -R(psi')^-1 (R f)_(>=0)(psi')`. In particular `E_R` annihilates `Omega(f, .)` (in either slot, `E_R` being
symmetric) iff `R f in Hc_-`. Scope: coefficients in `Q[[u]]`, infinite sums u-adically convergent (true in every application: D6
decay). Source RD2 s.3.1 (= "Lemma (L)" of the gauge worker, VERIFIED there with random covectors). Status PROVED; TRF s.2 CONFIRMED by
hand and `rf_lemE` (exact Q, random symplectic `R`, controls fail). Used by the chain (step 4). Depends on: the edge formula of
arXiv:1303.1043 s.2 (Def. 2.2 per S2A s.2) and (Sp) only.

**B2 `cor:tra-edges` (cor; "Corollary E'").** If `P = P_0 R` satisfies (Sp) and (A), then for every `k >= 1` the edges of `P` annihilate
`V^(k)` in each slot, because `(P z^-k J^C)_(>=0) = z^-k (P J^C)_(>=k) = 0`; equivalently every `Pi_i`, `i >= 0`, annihilates both slots.
Hence `E_P in K_0 (x) K_0` (completed; coefficients in `Q[[u]]`, `v = v(u)`), as `Q[[u]]`-combinations of frames (A3). Moreover
(A) <=> the edges annihilate `V^(1)` <=> they annihilate the whole flag (TRF s.2; RD3 s.4.2 (i)). Source RD2 s.3.1; second write-up
MS2 s.6 Lemma K(a) (B5). Status PROVED; TRF s.2 CONFIRMED; MRF item 6 CONFIRMED (Lemma K). Used by the chain (step 4: "edges in
`K_0 (x) K_0`"). Depends on: B1; `E_P = E_R` because `P_0^T eta P_0 = eta` ((Sp) at `z^0`) (TRF s.2); A1 for the frame expansion.

**B3 `lem:tra-leaf` (lem; "Lemma T").** Let `P` be as in D11 with (Sp) and (A). (a) `Omega(f, q_T) = Omega(P f, -z phi_0)` for every `f`.
(b) `V^(k)(q_T) = (-1)^(k+1) [z^(k-2)] Y'_0(-z)`, which is `T_P` for `k = 1` and `0` for `k >= 2` (by (A) for `k >= 3`, by the grading
for `k = 2`); so `Pi_0(q_T) = T_P` and `q_T in K`. (c) `q_T = -(z/d) phi_0 + T_(>=2)`, `T_(>=2) in z^2 Hc_+`; hence
`qt := d q_T = -z phi_0 + d T_(>=2)` lies in `-z + K` with `Pi_0(qt) = d T_P =: v + s`. (d) `s = d T_P - v in u^2 Q[[u]]`, because
`T_P = v'(0) (P_0)_00(0) u + O(u^2)` and `(P_0)_00 (P_0)_33 = 1` by (Sp) at `z^0`. Source RD2 s.3.1; second write-up MS2 s.6
Lemma K(b). Status PROVED; TRF s.2 CONFIRMED (equivalent form `V^(k)(q_T) = -[z^(k-2)] Y'_0`). Used by the chain (step 4).
Corrections: the shorthand "the leaf point lies in `-z + K` with `Pi_0 = t_c`" (SYN01 s.0, MECH1 s.1 item 3) means (b)-(c) with
`T_P = t_c`: `q_T in K`, `Pi_0(q_T) = t_c`, and the dilaton-rescaled `d q_T in -z + K` with `Pi_0 = d t_c = v + s` (S2A step 4 has the
precise form). TRF adds (VERIFIED): `s` is frame independent, `s = tau^(1/4) varpi_1 - v` with `varpi_0 t_c = varpi_1`,
`d = tau^(1/4) varpi_0`, and `s = u^2/30 + (127/3150) u^3 + ..` in both conifold frames. Depends on: (Sp), (A), grading, D6.

**B4 `lem:tra-finite` (lem; "Lemma F").** For each stable graph `Gamma`, `Cont_Gamma` is a FINITE sum of products of model correlators
at `v(u)` (monomials `C v(u)^e`, i.e. Laurent series in `u`) and of entries of `P_k`, `(P^-1)_k`, `d^(+-1)` (power series in `u`). In the
flat frame nothing has to converge; `P(0, z) != 1` is harmless. Source RD2 s.3.1. Status PROVED; TRF s.2 CONFIRMED. Used by the chain
(step 4: "every graph term is a finite sum"). Correction: RD2 writes "Laurent monomials in `u`"; read "monomials in `v(u)`".

**B5 `rmk:tra-K` (lem, alternative proof of B2-B3; "Lemma K" of the model note).** Let `R(z)` be symplectic (`R^*(-z) R(z) = 1`) with
edge `V(w, z) = (eta^-1 - (R^-1(w) (x) R^-1(z)) eta^-1)/(w + z)`, and `phi(w) = Sum_(j>=0) phi_j w^-j` a covector series with no positive
power of `w` such that `(R^-1(w))^T phi(w)` has none either. Put `Pi_m(f) := Sum_k phi_(k-m)(f_k)`. Then (a) `Pi_m` annihilates both slots
of `V` for every `m >= 0`; (b) for the leaf point `-z R^-1(z) 1 = -z + T(z)`: `Pi_m(-z + T) = -[w^(m-1)] <(R^-1(w))^T phi(w), 1>`, which
is 0 for `m >= 2`, `-<y_0, 1>` for `m = 1`, `-[w^-1] <(R^-1)^T phi, 1>` for `m = 0`. With `phi` = the model's vanishing covector and
`R` the gauge, the hypothesis is (G), and one gets "edges in `K_0 (x) K_0`, leaf point in `-z + K` with `Pi_0 = t_c`", in the canonical
normalisation only (for `P -> P c(z)`, `c` scalar, the hypothesis fails at the first non-zero `c_k`). Source MS2 s.6. Status PROVED;
MRF item 6 CONFIRMED ("elementary"; shorter: `Pi_m` is the `m = 0` functional of `w^-m phi(w)`). Correction (MRF item 6; SYN01 s.2
"Model"): the stated hypothesis "`phi_j -> 0` u-adically" FAILS in the hatted frame (`phi_(i,m) = O(u)` exactly); the conclusion
survives because every hatted edge is a polynomial in `(psi, psi')` at each `u`-order ("Theorem D", GAUGE_IDENTITY_HAT, not read); in the
flat frame `phi_j = O(u^(j - floor(j/5)))`. Rename (SYN01 C9: "Lemma K"/"K" clash with UR's frame `K`). Placement: remark (second proof).

**B6 `rmk:tra-fock` (rmk; theory's "Lemma C", Fock-space form).** If `f` and `P f` lie in `Hc_-`, then `P^ o (mult. by Omega(f, q)) =
(mult. by Omega(P f, q)) o P^`; hence `P^(Phi(V^(1), V^(2), ..) D) = Phi(V'^(1), V'^(2), ..) P^ D` for every function `Phi` of the
`V^(k)` (`V'` built from `Y' = P Y^C`): the regular gauge transports the whole c = 1 sector. Source RD2 s.3.1. Status PROVED given the
formula `(R^)G(q) = [exp(hbar E_R(d, d)/2) G](R^-1 q)` of arXiv:math/0108100 (conventions not re-checked, RD3 s.7); TRF s.2 CONFIRMED
as a reformulation. Not used by Theorem R or the chain. Rename (clash with BVRF's Lemma C, UNI-2 in piece 3).

**B7 `prop:tra-canonical` (prop; the quintic case).** `P_can` satisfies the hypotheses of D11 with `Y' = Y^C_Q` and `T_P = t_c`. Hence
its edges lie in `K_0 (x) K_0` as `Q[[u]]`-combinations of frames, its leaf point satisfies `q_T in K`, `Pi_0(q_T) = t_c`,
`d q_T in -z + K` with `Pi_0 = d t_c = v + s`, `s in u^2 Q[[u]]`, and every graph term is a finite sum. Source S2A s.2 step 4; RD1
s.1.4 example (i). Status PROVED given (U) (UR s.3), (G) (GI s.3), Prop. S (GI s.4), Lemmas V, Q, 1 (GI s.2); Theorem U applies to the
actual quintic only through (CYC) (S5A S9; PROVED in Stage 5). VERIFIED on the actual rows (s.5 items EV-T1..EV-T4). Used by the chain
(step 4). "Three independent write-ups" (MECH1 s.1 item 3): RD2 s.3.1, MS2 s.6, GI's Lemma (L).

## C. Quintic Theorem R [s.4]

**C1 `thm:red-main` (thm; "quintic Theorem R").** Fix `g >= 2`. Let `P` be a dressing of the model at `v(u)` (D11: flat-regular, graded,
`P_0` invertible) satisfying (Sp) and (A), and let `Omega^P = P.Omega^mod_(v(u))` (D12).
(a) If ML_K(h) holds for `0 <= h <= g - 1`, then `Cont_Gamma in Q[[u]]` for every stable graph `Gamma` of genus `g` with at least one edge.
(b) If moreover ML_K(g) holds, then `Cont_vertex = c_g T_P^(2-2g) + Q[[u]]`.
The gauge equation (flatness in `u`, the u-ODE) is not used. For `P = P_can` (B7) the vertex graph is `c_g t_c^(2-2g) + Q[[u]]`: this is
Conjecture N(g) for the flat sum `P_can.Omega^mod_(v(u))`; Theorem S15 (S2A step 6) transfers it to the hatted outer sum.
Source RD2 s.3.2 (general `P`, flat frame); second write-up MS2 s.7 (canonical gauge, hatted outer sum, assuming the factorisation
"Corollary F", hence (H)); S2A s.2 step 5. Status PROVED, conditional on ML_K; refereed twice: TRF s.0 item 1 and s.2 (CONFIRMED, no
gap; conditional also on the s.1.5 identification, now S15), MRF item 8 (accepted as a conditional reduction; no gap beyond B5's
hypothesis). Final: with ML_K PROVED mod (Q) (S2A s.1), Theorem R is unconditional modulo (Q) for any `P` as in D11. Used by the chain:
(a), (b) (step 5).
Corrections: name "quintic Theorem R" (SYN01 C9: clash with the conifold theorem); RD2 s.3.3 Cor. N / N_9 omit the s.1.5 hypothesis
(TRF s.0 item 3(iv); SYN01 s.7) and are superseded (C6 below); "each of the three hypotheses on `P` is needed in examples" is WEAKENED:
there is no example for (Sp) (TRF s.0 item 3(v); SYN01 s.2); MS2 s.7's hypothesis "(iii) ML_h for every `h >= 0`" is more than needed
(`h <= g` suffices). Depends on: A1-A3, B1-B4, ML_K (D14), the dilaton equation (D4), arXiv:1303.1043 s.2.
(c) (in RD2, not needed): the ancestor potential of `Omega^P` restricted to its own slice `Lambda'` (built from `Y'`) satisfies ML(g)
with `V'^(1)`; legs `w` enter as `P^-1 w` and `V^(k)(P^-1 w) = V'^(k)(w)`. Status: one-line outline (TRF s.2 "correct one-line
outline"; RD3 s.7). Omit, or remark with a \draftnote.

**C2 `rmk:red-hyp` (rmk; what is used of ML).** (i) Only on `-z + K`: insertions from `K_0` and `tau_0(H)`; nothing about the
D-directions `tau_k(H)`, `k >= 1`, alone. (ii) Only the vanishing direction: singular correlators with at least one `K_0` insertion
vanish; no formula for a non-zero singular part except the gap. (iii) Graphs with edges: genera `h <= g - 1` only, including `h = 0, 1`.
(iv) Vertex graph: the gap `<tau_0(H)^m>_g = Phi_g^(m)(v)` + regular (equivalently `F_g = c_g v^(2-2g)`) and the `n >= 1` statement at
`h = g`. (v) For each `g` only the finitely many correlators that reach hatted order `<= 2g - 3`. Source RD2 s.3.2 ("Hypotheses
needed"); TRF s.2 (step (2) needs ML only for lists with at least one `K_0` insertion). Status PROVED (it is a reading of the proof).

**C3 `cor:red-leaf` (cor or rmk; "Corollary L", leaf bookkeeping).** Given ML_K(g): in the normalisation with leafless vertex
`c_g vhat^(2-2g)` (`vhat = v/d`), the leaves turn it into `c_g t_c^(2-2g)` + regular, whose first correction is
`-(g-1) c_g (2/vhat^(2g-1)) (t_c - vhat)`; with `t_c - vhat = u^2/30 + O(u^3)` this is `-(g-1) c_g u^(3-2g)/15`, the same in every frame
(SYN01 C6: `2 gamma_t + 2 w_1 - 4/3 = 1/15`, `gamma_t + w_1 = 7/10`). Source RD2 s.3.3; TRF s.2 (CONFIRMED given ML; genus-2 control
X1: leaves removed give the leafless `-1/900 u^-1 = -(3g-3)(4/9) c_g`, the leaves supply `-1/18000 = -(g-1) c_g/15`). Status PROVED given
ML(g), with the input `t_c - vhat = u^2/30 + O(u^3)` VERIFIED only (`gi_hat_table.log`; TRF `rf_flat`; `cr_conv` (4)). Not used by the
chain. Placement: remark; the number `1/30` goes to the evidence appendix unless a written derivation is found (s.7).

**C4 `prop:red-gaugeclass` (prop; "Proposition I").** ML(h <= g) is an invariant of the gauge class: if it holds for one system of the
shape (canonically normalised), it holds for all (apply C1(c) to `P_(A') P_A^-1`, which is flat-regular with (Sp) and (A)). Source RD3
s.4.3. Status PROVED from the model to any system; SKETCH from a system to another or back to the model (homogeneity of the source
replaced by the hatted u-filtration, NF s.6.1); TRF s.0 item 4 ("only model-outwards"). Not used. Omit or remark.

**C5 `rmk:red-deformation` (rmk; deformation route).** Join the pulled-back model `A_0` to `A_1` by a path `A_s` with the same product
and residue; `X_s := (d_s P_s) P_s^-1` satisfies (D0) `z X_s' = [A_s, X_s] + d_s A_s`; (D1) flat-regular, graded, infinitesimally
symplectic; (D2) `X_s Y^C_(A_s) = d_s Y^C_(A_s)` (no positive power: annihilates the flag); (D3) `d_s F_g = (1/2) Sum B_s(dF_(g1), dF_(g2))
+ (1/2) B_s(d^2 F_(g-1)) + L_s F_g` with `B_s in Sym^2 K_0`. The Hessian of a naive non-linear `Vn` enters only on `K_0 (x) K_0`, where
ML has no Hessian. It is the infinitesimal form of C1(c) and gives nothing beyond it. Source RD3 s.4.1-s.4.3. Status PROVED
(reformulation; TRF s.2 CONFIRMED). Not used. Optional remark.

**C6 Superseded statements (do not transcribe as such).** RD2 s.3.3 Corollary N (needs Theorem U, (H), (N), ML) and Corollary N_9
(genus `<= 9` reduced to ML(h <= g)): superseded by the assembled theorem (MECH5 s.2): N(g) for all `g >= 2` modulo (U), (Q), (QRR).
The s.1.5 identification (RD1 s.1.5) is no longer an input: Theorem S15 (S2B K8). RD1 s.0 item 3 ("ML VERIFIED on 1013 lists") and
s.7 (ML "CONJECTURAL") are superseded by ML_K PROVED mod (Q) (S2A s.1). MECH1 s.1's "four inputs" and SYN01 s.4's statement are
superseded by MECH2 s.1 and then MECH5 s.2.
